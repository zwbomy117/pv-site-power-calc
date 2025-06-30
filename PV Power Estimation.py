# -*- coding: utf-8 -*-
"""
Title: Site-level PV Power Output Estimation using pvlib
Author: Zhao Wenbo
Contact: zwbomy117@163.com
Last Modified: 2025-06-30

Description:
    This script provides core functions for estimating Direct Normal Irradiance (DNI) and photovoltaic (PV)
    power output at a site level using the pvlib library. The workflow includes solar position calculation,
    irradiance decomposition, module temperature modeling, and DC/AC conversion based on the Sandia model.

    This script offers the **foundational code** for PV power estimation.
    To implement a complete workflow (e.g., batch station simulation, file IO, capacity factor calculation),
    users should connect these functions with their own datasets and call them accordingly.

    The `multiprocessing` library is recommended for large-scale or multi-station processing,
    as it significantly speeds up the computation. For technical support or collaboration,
    please contact the author via: zwbomy117@163.com

    > This work builds upon the excellent functionality provided by the **pvlib** library.
    > Special thanks to the pvlib development team for their outstanding open-source contributions to the solar energy research community.
"""

import pandas as pd
import numpy as np
import os
import pvlib
from pvlib import location

# Set PROJ_LIB path if needed (for solar position/coordinate projections in pvlib)
os.environ['PROJ_LIB'] = r'C:\ProgramData\Anaconda3\envs\pvlib\Lib\site-packages\osgeo\data\proj'

# Load PV module and inverter parameters from the SAM database
sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
sapm_inverters = pvlib.pvsystem.retrieve_sam('cecinverter')

# Specify the PV module and inverter models
module = sandia_modules['Canadian_Solar_CS5P_220M___2009_']
inverter = sapm_inverters['ABB__MICRO_0_25_I_OUTD_US_208__208V_']

# Temperature model parameters for module mounting configuration
temperature_model_parameters = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']


# ----------------------------------------------------------------------------------------
# Function 1: Compute DNI from GHI/DHI and solar position
# ----------------------------------------------------------------------------------------
def dot_dni_series_calculate(Location, ghi_series, dhi_series, temperature_series, tz):
    """
    Estimate the Direct Normal Irradiance (DNI) time series for a single site.

    Parameters:
        Location (tuple): (latitude, longitude, altitude) in decimal degrees and meters.
        ghi_series (pd.Series): Time-indexed Global Horizontal Irradiance (W/m²).
        dhi_series (pd.Series): Time-indexed Diffuse Horizontal Irradiance (W/m²).
        temperature_series (pd.Series): Time-indexed ambient temperature (°C).
        tz (str): Timezone of the location (e.g., 'UTC').

    Returns:
        pd.Series: Time series of estimated DNI (W/m²).
    """
    latitude, longitude, altitude = Location
    loc = location.Location(latitude=latitude, longitude=longitude, altitude=altitude, tz=tz)

    # Calculate solar position for each timestamp
    solpos = pvlib.solarposition.get_solarposition(
        time=ghi_series.index,
        latitude=latitude,
        longitude=longitude,
        altitude=altitude,
        temperature=temperature_series,
        pressure=pvlib.atmosphere.alt2pres(altitude)
    )

    # Prepare data frame for clear-sky modeling
    clear_input = pd.DataFrame({
        'utc_time': ghi_series.index,
        'apparent_zenith': solpos['apparent_zenith'],
        'zenith': solpos['zenith'],
        'apparent_elevation': solpos['apparent_elevation']
    }).set_index('utc_time', drop=True)

    # Estimate clear-sky DNI for quality control
    clear_sky_dni = loc.get_clearsky(times=clear_input.index, solar_position=clear_input)

    # Compute actual DNI using pvlib’s decomposition model
    dni = pvlib.irradiance.dni(
        ghi=ghi_series,
        dhi=dhi_series,
        zenith=solpos['zenith'].astype(np.float32),
        clearsky_dni=clear_sky_dni['dni'],
        clearsky_tolerance=1.0,  # Stricter threshold for exceeding clear-sky
        zenith_threshold_for_zero_dni=85,  # High zenith → zero DNI
        zenith_threshold_for_clearsky_limit=30  # Optional: mark too-high DNI as invalid
    )

    return dni


# ----------------------------------------------------------------------------------------
# Function 2: Estimate PV power output (AC) for a single location
# ----------------------------------------------------------------------------------------
def dot_pv_series_calculate(Location, ghi_series, dni_series, dhi_series, temperature_series, windspeed_series):
    """
    Estimate the photovoltaic AC power output time series for a single site.

    Parameters:
        Location (tuple): (latitude, longitude, altitude) of the station.
        ghi_series (pd.Series): Global Horizontal Irradiance (W/m²).
        dni_series (pd.Series): Direct Normal Irradiance (W/m²).
        dhi_series (pd.Series): Diffuse Horizontal Irradiance (W/m²).
        temperature_series (pd.Series): Ambient temperature (°C).
        windspeed_series (pd.Series): Wind speed at hub height or surface (m/s).

    Returns:
        pd.Series: AC power output time series (W).
    """
    latitude, longitude, altitude = Location
    ghi = ghi_series
    dni = dni_series
    dhi = dhi_series
    temperature = temperature_series
    windspeed = windspeed_series

    # PV system configuration: south-facing, tilt = latitude
    system = {
        'module': module,
        'inverter': inverter,
        'surface_azimuth': 180,
        'surface_tilt': latitude
        # You may replace latitude with a predefined optimal tilt angle if needed
    }

    # Compute solar position
    solpos = pvlib.solarposition.get_solarposition(
        time=ghi.index,
        latitude=latitude,
        longitude=longitude,
        altitude=altitude,
        temperature=temperature,
        pressure=pvlib.atmosphere.alt2pres(altitude)
    )

    # Get extra-terrestrial DNI and airmass for irradiance modeling
    dni_extra = pvlib.irradiance.get_extra_radiation(ghi.index)
    airmass = pvlib.atmosphere.get_relative_airmass(solpos['apparent_zenith'])
    pressure = pvlib.atmosphere.alt2pres(altitude)
    am_abs = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)

    # Angle of incidence (AOI) on the module surface
    aoi = pvlib.irradiance.aoi(
        surface_tilt=system['surface_tilt'],
        surface_azimuth=system['surface_azimuth'],
        solar_zenith=solpos['apparent_zenith'],
        solar_azimuth=solpos['azimuth']
    )

    # Compute total POA (plane of array) irradiance
    total_irradiance = pvlib.irradiance.get_total_irradiance(
        surface_tilt=system['surface_tilt'],
        surface_azimuth=system['surface_azimuth'],
        solar_zenith=solpos['apparent_zenith'],
        solar_azimuth=solpos['azimuth'],
        dni=dni,
        ghi=ghi,
        dhi=dhi,
        dni_extra=dni_extra,
        model='haydavies'  # Anisotropic sky model
    )

    # Estimate cell/module temperature
    cell_temperature = pvlib.temperature.sapm_cell(
        poa_global=total_irradiance['poa_global'],
        temp_air=temperature,
        wind_speed=windspeed,
        **temperature_model_parameters
    )

    # Estimate effective irradiance received by the PV cell
    effective_irradiance = pvlib.pvsystem.sapm_effective_irradiance(
        poa_direct=total_irradiance['poa_direct'],
        poa_diffuse=total_irradiance['poa_diffuse'],
        airmass_absolute=am_abs,
        aoi=aoi,
        module=module
    )

    # Estimate DC output using the SAPM model
    dc = pvlib.pvsystem.sapm(effective_irradiance, cell_temperature, module)

    # Convert to AC output using inverter model
    ac = pvlib.inverter.sandia(v_dc=dc['v_mp'], p_dc=dc['p_mp'], inverter=inverter)

    return ac

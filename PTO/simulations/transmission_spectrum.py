#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:26:10 2023

@author: chamaeleontis
"""

#%%
import pandas as pd
from enum import Enum
from dataclasses import dataclass
import numpy as np
import astropy.units as u
import specutils as sp
import os
import scipy as sci
from expecto import get_spectrum
import PTO.spectra_manipulation as sm
import astropy
from astropy.modeling import models
import astropy.constants as con
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy as sci
import seaborn as sns
import matplotlib.gridspec as gs
import logging
from PTO.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 
#%%
@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class ESPRESSO_wavelength:
    """Class holding the wavelength range of ESPRESSO."""
    
    def __post_init__(self):
        """
        Estimated wavelength range as taken from real observation.

        Returns
        -------
        None.

        """
        self.wavelength_range = np.array([
        np.linspace(3771.2787376605725, 3829.7840287759504,9111),
        np.linspace(3771.4643227636125, 3829.7715446057596,9111),
        np.linspace(3794.8947452174707, 3853.7087472919184,9111),
        np.linspace(3794.8425844119834, 3853.7583061178107,9111),
        np.linspace(3818.734013707605, 3877.9601104946737,9111),
        np.linspace(3818.7323760406034, 3877.985853468359,9111),
        np.linspace(3842.8916578752874, 3902.4891093367537,9111),
        np.linspace(3842.878990395664, 3902.5312931681874,9111),
        np.linspace(3867.392129559423, 3927.3494595582843,9111),
        np.linspace(3867.375688541572, 3927.3841788348454,9111),
        np.linspace(3892.2116576126086, 3952.516259845867,9111),
        np.linspace(3892.2011293018563, 3952.5373701159447,9111),
        np.linspace(3917.348089087959, 3978.0193594390835,9111),
        np.linspace(3917.3175168008233, 3978.0408962098973,9111),
        np.linspace(3942.773556340861, 4003.84729094209,9111),
        np.linspace(3942.758392522059, 4003.8661996595506,9111),
        np.linspace(3968.5667515910864, 4030.0123100342607,9111),
        np.linspace(3968.554500551588, 4030.0402519983454,9111),
        np.linspace(3994.71166352614, 4056.520478527957,9111),
        np.linspace(3994.7084094329925, 4056.5395076678797,9111),
        np.linspace(4021.152737825885, 4083.400802019439,9111),
        np.linspace(4021.161605847676, 4083.4065133247295,9111),
        np.linspace(4047.964187098823, 4110.620833215921,9111),
        np.linspace(4047.9791575688255, 4110.623376891031,9111),
        np.linspace(4075.152139524644, 4138.2062390193705,9111),
        np.linspace(4075.158760668843, 4138.212937489863,9111),
        np.linspace(4102.7377338066835, 4166.148710147385,9111),
        np.linspace(4102.742617012644, 4166.1633918340995,9111),
        np.linspace(4130.685955536411, 4194.475122121532,9111),
        np.linspace(4130.708603652426, 4194.47012183829,9111),
        np.linspace(4158.974557576077, 4223.213557824446,9111),
        np.linspace(4158.985540650628, 4223.209882585936,9111),
        np.linspace(4187.718417690256, 4252.317160184001,9111),
        np.linspace(4187.723469369687, 4252.317266254806,9111),
        np.linspace(4216.815284111027, 4281.843072416275,9111),
        np.linspace(4216.8126207053765, 4281.845932830434,9111),
        np.linspace(4246.334816103123, 4311.781048994436,9111),
        np.linspace(4246.336166010289, 4311.7829533859785,9111),
        np.linspace(4276.267489255504, 4342.12368830769,9111),
        np.linspace(4276.260066529987, 4342.139017964975,9111),
        np.linspace(4306.617164207534, 4372.913822026877,9111),
        np.linspace(4306.627237001473, 4372.912142438144,9111),
        np.linspace(4337.415108152269, 4404.143758486397,9111),
        np.linspace(4337.405867239342, 4404.154172603916,9111),
        np.linspace(4368.646015983025, 4435.816763539261,9111),
        np.linspace(4368.6424339489095, 4435.831662164844,9111),
        np.linspace(4400.32979584596, 4467.950330415987,9111),
        np.linspace(4400.333093321937, 4467.956109252706,9111),
        np.linspace(4432.496780549474, 4500.547470691217,9111),
        np.linspace(4432.48615094564, 4500.55851877716,9111),
        np.linspace(4465.095317628195, 4533.638748247901,9111),
        np.linspace(4465.109462411758, 4533.647418063654,9111),
        np.linspace(4498.215180887234, 4567.208881036579,9111),
        np.linspace(4498.235624922552, 4567.2131507180275,9111),
        np.linspace(4531.831642314217, 4601.27131696573,9111),
        np.linspace(4531.830471130549, 4601.28398962037,9111),
        np.linspace(4565.946648462809, 4635.8462375713325,9111),
        np.linspace(4565.955552378357, 4635.8602927238935,9111),
        np.linspace(4600.578010724452, 4670.948257607756,9111),
        np.linspace(4600.581788340105, 4670.963301415562,9111),
        np.linspace(4635.731478076988, 4706.597847317468,9111),
        np.linspace(4635.756768476915, 4706.602565432542,9111),
        np.linspace(4671.437658850833, 4742.781160357433,9111),
        np.linspace(4671.455528811719, 4742.7893308012935,9111),
        np.linspace(4707.707058862909, 4779.525992998565,9111),
        np.linspace(4707.712056736231, 4779.543013013408,9111),
        np.linspace(4744.530994505388, 4816.853433769093,9111),
        np.linspace(4744.545583582409, 4816.856300909258,9111),
        np.linspace(4781.943564819925, 4854.761467764942,9111),
        np.linspace(4781.957715035795, 4854.763825088372,9111),
        np.linspace(4819.950406606762, 4893.256216867388,9111),
        np.linspace(4819.964376405277, 4893.26810437854,9111),
        np.linspace(4858.571811727038, 4932.384191553762,9111),
        np.linspace(4858.593223196309, 4932.384453545037,9111),
        np.linspace(4897.824173676835, 4972.12770683428,9111),
        np.linspace(4897.8425261785505, 4972.1319312744145,9111),
        np.linspace(4937.715781689066, 5012.528234585245,9111),
        np.linspace(4937.726010671691, 5012.529273781539,9111),
        np.linspace(4978.264194020407, 5053.576337829653,9111),
        np.linspace(4978.266144421772, 5053.577786403188,9111),
        np.linspace(5019.480251279039, 5095.303058344128,9111),
        np.linspace(5019.502423188249, 5095.300595221686,9111),
        np.linspace(5061.398044204499, 5137.718955527876,9111),
        np.linspace(5061.428149288104, 5137.713415802885,9111),
        np.linspace(5104.021036208722, 5180.8521344122155,9111),
        np.linspace(5104.040009148319, 5180.848331757356,9111),
        np.linspace(5147.383791067555, 5224.699961813259,9111),
        np.linspace(5147.398516981898, 5224.69572769761,9111),
        np.linspace(5191.419636147633, 5269.771277282278,9111),
        np.linspace(5191.473229341502, 5269.926193480165,9111),
        np.linspace(5193.165101663176, 5271.210731684364,9111),
        np.linspace(5193.1726471525635, 5271.222218453554,9111),
        np.linspace(5237.945746567275, 5316.644531405807,9111),
        np.linspace(5237.952960964355, 5316.655663517843,9111),
        np.linspace(5283.502716463571, 5362.8675364615865,9111),
        np.linspace(5283.509528476912, 5362.878526816698,9111),
        np.linspace(5329.861894507855, 5409.901322700275,9111),
        np.linspace(5329.867549387057, 5409.912229807305,9111),
        np.linspace(5377.041935747941, 5457.7675282508035,9111),
        np.linspace(5377.0468985671805, 5457.7784089527595,9111),
        np.linspace(5425.064514018341, 5506.488524610052,9111),
        np.linspace(5425.068972228696, 5506.499361670704,9111),
        np.linspace(5473.952645168638, 5556.087087948789,9111),
        np.linspace(5473.956423954884, 5556.0980360513795,9111),
        np.linspace(5523.729388783924, 5606.5873867712035,9111),
        np.linspace(5523.732881032934, 5606.5981670616975,9111),
        np.linspace(5574.419586840785, 5658.0140978728605,9111),
        np.linspace(5574.422692450903, 5658.0249486882885,9111),
        np.linspace(5626.048572259214, 5710.392948405479,9111),
        np.linspace(5626.051351328716, 5710.403888723457,9111),
        np.linspace(5678.642794816133, 5763.750582183552,9111),
        np.linspace(5678.645331079689, 5763.7616417701665,9111),
        np.linspace(5732.229491532979, 5818.114816027069,9111),
        np.linspace(5732.232006565761, 5818.1258171567215,9111),
        np.linspace(5786.837132880204, 5873.514336214492,9111),
        np.linspace(5786.8394946423305, 5873.525449553423,9111),
        np.linspace(5842.494904595001, 5929.978617450081,9111),
        np.linspace(5842.497433968182, 5929.989907592699,9111),
        np.linspace(5899.233913912579, 5987.538922086424,9111),
        np.linspace(5899.236556004323, 5987.550061561587,9111),
        np.linspace(5957.086176968194, 6046.227490000123,9111),
        np.linspace(5957.088934489347, 6046.238709426095,9111),
        np.linspace(6016.084177088486, 6106.078267686153,9111),
        np.linspace(6016.08725313049, 6106.0893861676395,9111),
        np.linspace(6076.262857711647, 6167.125325878146,9111),
        np.linspace(6076.266076580646, 6167.1366289928455,9111),
        np.linspace(6137.657907914385, 6229.404373678272,9111),
        np.linspace(6137.661406707827, 6229.415691771861,9111),
        np.linspace(6200.30675391242, 6292.9541900770455,9111),
        np.linspace(6200.310666998568, 6292.965482422278,9111),
        np.linspace(6264.247997746767, 6357.813632038661,9111),
        np.linspace(6264.252101411256, 6357.824986062797,9111),
        np.linspace(6329.522063862854, 6424.023730359819,9111),
        np.linspace(6329.526637929429, 6424.035000211732,9111),
        np.linspace(6396.171398644979, 6491.627130208693,9111),
        np.linspace(6396.176290581616, 6491.638354238909,9111),
        np.linspace(6464.239907748623, 6560.667756554076,9111),
        np.linspace(6464.245126639533, 6560.679116522222,9111),
        np.linspace(6533.7739144007555, 6631.191793458423,9111),
        np.linspace(6533.7795188669, 6631.202824333651,9111),
        np.linspace(6604.8209638476, 6703.247449304857,9111),
        np.linspace(6604.826851277418, 6703.258548242413,9111),
        np.linspace(6677.4307840135225, 6776.886068150898,9111),
        np.linspace(6677.436923658317, 6776.896946193677,9111),
        np.linspace(6751.655903724899, 6852.159980454181,9111),
        np.linspace(6751.662316544761, 6852.170723078734,9111),
        np.linspace(6827.551570336545, 6929.123564916879,9111),
        np.linspace(6827.558163145319, 6929.1342335409945,9111),
        np.linspace(6905.17415275321, 7007.834726179728,9111),
        np.linspace(6905.180888332385, 7007.845407905051,9111),
        np.linspace(6984.583933873233, 7088.353126205972,9111),
        np.linspace(6984.590793907555, 7088.363399711146,9111),
        np.linspace(7065.8438036309, 7170.741481398268,9111),
        np.linspace(7065.850758167719, 7170.751480893239,9111),
        np.linspace(7149.018578697368, 7255.066405834495,9111),
        np.linspace(7149.025490481778, 7255.076102379623,9111),
        np.linspace(7234.177685185124, 7341.396069814298,9111),
        np.linspace(7234.184405943058, 7341.405319379039,9111),
        np.linspace(7321.393319139038, 7429.802823133282,9111),
        np.linspace(7321.399516755308, 7429.811870288571,9111),
        np.linspace(7410.740551835409, 7520.362761154747,9111),
        np.linspace(7410.746201321001, 7520.371734220259,9111),
        np.linspace(7502.300503402569, 7613.153017058346,9111),
        np.linspace(7502.305877581938, 7613.161306250178,9111),
        np.linspace(7596.156676007851, 7708.2574877288125,9111),
        np.linspace(7596.161035747009, 7708.265607093863,9111),
        np.linspace(7692.3946690851635, 7805.7661139164575,9111),
        np.linspace(7692.398574597109, 7805.774664521214,9111),
        np.linspace(7791.110253537907, 7905.5145802111865,9111),
        np.linspace(7791.1138966871695, 7905.599289073489,9111),
        ])
        self.wavelength_range = np.sort(self.wavelength_range.flatten())
        self.wavelength_range *= u.AA
        return

#%%
class Instruments(Enum):
    """Enum class holding wavelength ranges for instruments."""
    
    ESPRESSO = ESPRESSO_wavelength().wavelength_range
    ESPRESSO_4UT = ''
    HARPS = ''
    NIRPS = ''
    # # FIXME - This is not correct, but no official ETC for ANDES is provided yet
    # ANDES = ANDES_wavelength_range.wavelength_range
    
    
class Instruments_size(Enum):
    """Enum class holding wavelength ranges for instruments."""
    
    ESPRESSO = 8.2 * u.m
    ESPRESSO_4UT = 4*8.2*u.m
    HARPS = 3.6*u.m
    NIRPS = 3.6*u.m
    ANDES = 39 * u.m

#%% Scaling by instrument
def _scale_instrument(spectrum: sp.Spectrum1D,
                      from_instrument:str,
                      to_instrument: str):
    """
    Scale the SNR between instrument, using the telescope mirror size. Instrument must be within the Instruments_size Enum class.

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Spectrum to scale
    from_instrument : str
        Original instrument, usually ESPRESSO, for which the code was implemented.
    to_instrument : str
        Target instrument, for which we want to simulate the dataset.
    """
    
    scale_factor = (Instruments_size[to_instrument].value / Instruments_size[from_instrument].value)**2 
    
    new_spectrum = sp.Spectrum1D(
        spectral_axis = spectrum.spectral_axis,
        flux = spectrum.flux * scale_factor,
        uncertainty = astropy.nddata.StdDevUncertainty(
            np.sqrt(spectrum.flux * scale_factor)
            )
        mask = np.zeros_like(spectrum.flux),
        meta= {}
        )
    
    return new_spectrum

#%% vactoair
def vactoair(wlnm: np.ndarray):
    """
    Change from vacuum to air wavelengths.

    Parameters
    ----------
    wlnm : np.ndarray
        Wavelenght grid in nanometers (actually invariant towards units)
    """
    wlA = wlnm*10.0
    s = 1e4/wlA
    f = 1.0 + 5.792105e-2/(238.0185e0 - s**2) + 1.67917e-3/( 57.362e0 - s**2)
    return(wlA/f/10.0)

#%%
def generate_mock_spectrum(system_parameters: pd.DataFrame,
                           instrument: str = 'ESPRESSO') -> sp.Spectrum1D:
    """
    Generate mock spectrum using the PHOENIX high resolution spectra.

    Parameters
    ----------
    system_parameters : pd.DataFrame
        Stellar parameters to consider as provided by NASA archive.
    instrument : str, optional
        Which instrument to consider. The default is 'ESPRESSO'.

    Returns
    -------
    best_scenarion_spectrum : sp.Spectrum1D
        Mock spectrum interpolated to instrument wavelength grid assuming the best weather conditions.
    average_scenarion_spectrum : sp.Spectrum1D
        Mock spectrum interpolated to instrument wavelength grid assuming the average weather conditions.
    worst_scenarion_spectrum : sp.Spectrum1D
        Mock spectrum interpolated to instrument wavelength grid assuming the worst weather conditions.

    """
    stellar_spectrum = get_spectrum(
        T_eff=system_parameters['st_teff'],
        log_g=system_parameters['st_logg'],
        cache=True
    )
    
    stellar_spectrum = sp.Spectrum1D(spectral_axis = vactoair(stellar_spectrum.spectral_axis.value) * u.AA,
                             flux = stellar_spectrum.flux,
                             uncertainty = astropy.nddata.StdDevUncertainty(
                                 np.zeros(stellar_spectrum.flux.shape)
                             ),
                             )
    
    stellar_spectrum = sm.interpolate2commonframe(stellar_spectrum,Instruments[instrument].value)
    
    best_scenario, average_scenario, worst_scenario = load_SNR_table_ESPRESSO(system_parameters,
                                                                     instrument = instrument)
    
    best_scenario = scale_estimated_SNR_with_Vmag(best_scenario, system_parameters)
    average_scenario = scale_estimated_SNR_with_Vmag(average_scenario, system_parameters)
    worst_scenario = scale_estimated_SNR_with_Vmag(worst_scenario, system_parameters)
    
    # FIXME: This is ESPRESSO specific, do a branch for other instrument. 
    # Right now, this is fixed by _scale_instrument function, but should be unified
    best_scenario_spectrum = scale_spectrum_with_SNR(stellar_spectrum, best_scenario)
    average_scenario_spectrum = scale_spectrum_with_SNR(stellar_spectrum, average_scenario)
    worst_scenario_spectrum = scale_spectrum_with_SNR(stellar_spectrum, worst_scenario)
    
    if instrument != 'ESPRESSO':
        best_scenario_spectrum = _scale_instrument(best_scenario_spectrum, 'ESPRESSO', instrument)
        average_scenario_spectrum = _scale_instrument(average_scenario_spectrum, 'ESPRESSO', instrument)
        worst_scenario_spectrum = _scale_instrument(worst_scenario_spectrum, 'ESPRESSO', instrument)
    
    return best_scenario_spectrum, average_scenario_spectrum, worst_scenario_spectrum

#%%
def load_ETC_file(filename: str) -> np.ndarray:
    """
    Load the ETC files based on filename.

    Parameters
    ----------
    filename : str
        Filename location.

    Returns
    -------
    SNR : np.ndarray
        Array of [wavelength, expected SNR] as loaded from the file.

    """
    f = open(filename,'r')
    
    SNR =[]
    for line in f.readlines():
        if line.startswith('#'):
            continue
        
        wavelength, snr = line.replace('\n','').split('\t')
        
        SNR.append([float(wavelength), float(snr)])
    
    return np.asarray(SNR)

#%%
def find_stellar_type_ESPRESSO(T_eff: float) -> str:
    """
    Find stellar type based on effective temperature.

    Parameters
    ----------
    T_eff : float
        Effective temperature of the star in K.

    Returns
    -------
    spectral_type : str
        Spectral type closest to the true one that is available in ESPRESSO ETC.

    """
    if T_eff > 54000:
        return 'O5'
    elif T_eff > 43300:
        return 'O9'
    elif T_eff > 29200:
        return 'B1'
    elif T_eff > 23000:
        return 'B3'
    elif T_eff > 15000:
        return 'B8'
    elif T_eff > 11800:
        return 'B9'
    elif T_eff > 10000:
        return 'A0'
    elif T_eff > 9450:
        return 'A1'
    elif T_eff > 8480:
        return 'F0'
    elif T_eff > 5900:
        return 'G0'
    elif T_eff > 5560:
        return 'G2'
    elif T_eff > 4730:
        return 'K2'
    elif T_eff > 3865:
        return 'K7'
    else:
        return 'M2'
    
#%%
def calculate_semimajor_axis(system_parameters: pd.DataFrame)->pd.DataFrame:
    """
    Calculate semimajor axis provided all parameters.

    Parameters
    ----------
    system_parameters : pd.DataFrame
        System parameters as extracted from NASA archive.

    Raises
    ------
    ValueError
        If all the necessary components for calculation are not found, raise ValueError.

    Returns
    -------
    system_parameters : pd.DataFrame
        Semimajor axis of the system in units of au.

    """
    if np.isfinite(system_parameters['pl_orbsmax']):
        return system_parameters
    
    if (np.isfinite(system_parameters['st_mass']) and
        np.isfinite(system_parameters['pl_bmasse']) and
        np.isfinite(system_parameters['pl_orbper'])
        ):
        a = ((con.G * (
                system_parameters['st_mass']*u.M_sun +
                system_parameters['pl_bmasse']*u.M_earth
                ).decompose() *
            ((system_parameters['pl_orbper'] * u.day) **2)).decompose() /
            np.pi**2 / 4
            )**(1/3)
        system_parameters['pl_orbsmax'] = a.to(u.au).value
        return system_parameters
    else:
        raise ValueError ("System parameters are not complete for calculation of semimajor axis")

#%%
def calculate_inclination(system_parameters: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate inclination of the system.

    Parameters
    ----------
    system_parameters : pd.DataFrame
        System parameters as provided by the NASA archive.

    Raises
    ------
    ValueError
        On incomplete system parameters table.

    Returns
    -------
    system_parameters : pd.DataFrame
        Value of the eccentricity.

    """
    if np.isfinite(system_parameters['pl_orbincl']):
        return system_parameters
    
    if (np.isfinite(system_parameters['pl_orbeccen']) and
        np.isfinite(system_parameters['pl_orblper'])
        ):
        eccentricity_correction = ((1 + system_parameters['pl_orbeccen'] *
                                    np.sin(system_parameters['pl_orblper'] * np.pi/180)
                                    ) /
                                   (1 - system_parameters['pl_orbeccen']**2))
    else:
        eccentricity_correction = 1
    
    if (np.isfinite(system_parameters['st_rad']) and
        np.isfinite(system_parameters['pl_imppar']) and
        np.isfinite(system_parameters['pl_orbsmax'])
        ):
        i = np.arccos(
                (system_parameters['st_rad'] * u.R_sun *
                system_parameters['pl_imppar'] *
                eccentricity_correction /
                (system_parameters['pl_orbsmax'] * u.au)).decompose()
                )
        system_parameters['pl_orbincl'] = i.to(u.deg).value
        return system_parameters
    else:
        raise ValueError ("System parameters are not complete for calculation of inclination")

#%%
def load_SNR_table_ESPRESSO(system_parameters: pd.DataFrame,
                   instrument: str = 'ESPRESSO'
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load SNR tables based on the system parameters and instrument used.

    Parameters
    ----------
    system_parameters : pd.DataFrame
        System parameters as provided by NASA archive.
    instrument : str, optional
        Which instrument to use. The default is 'ESPRESSO'.

    Returns
    -------
    best_scenario : np.ndarray
        Array of [wavelength, expected SNR] for best case scenario.
    average_scenario : np.ndarray
        Array of [wavelength, expected SNR] for average (50%) case scenario.
    worst_scenario : np.ndarray
        Array of [wavelength, expected SNR] for worst case scenario.

    """
    spectral_type = find_stellar_type_ESPRESSO(float(system_parameters['st_teff']))
    location_of_ETC_files = './ETC_{instrument}/{spectral_type}/'.format(instrument = instrument, spectral_type = spectral_type)
    
    best_scenario = load_ETC_file(location_of_ETC_files + '0.5arcsec_900s.etc')
    average_scenario = load_ETC_file(location_of_ETC_files + '0.8arcsec_900s.etc')
    worst_scenario = load_ETC_file(location_of_ETC_files + '1.3arcsec_900s.etc')
    
    return best_scenario, average_scenario, worst_scenario
#%%
def scale_estimated_SNR_with_Vmag(scenario: np.ndarray,
                                  system_parameters: pd.DataFrame) -> np.ndarray:
    """
    Scale estimated SNR with V magnitude of the target.

    Parameters
    ----------
    scenario : np.ndarray
        Scenario as loaded by load_SNR_table function.
    system_parameters : pd.DataFrame
        System parameters as loaded by NASA archive.

    Returns
    -------
    scenario : np.ndarray
        Adjusted SNR of given scenario.

    """
    scenario[:,1] *= np.sqrt(10**(0.4*(10 - system_parameters['sy_vmag'])))
    
    return scenario

#%%
def scale_spectrum_with_SNR(spectrum: sp.Spectrum1D,
                            SNR: np.ndarray) -> sp.Spectrum1D:
    """
    Scale spectrum according to estimated SNR.

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Model spectrum to scale.
    SNR : np.ndarray, optional
        Estimated SNR of single exposure as provided by load_SNR_table.

    Returns
    -------
    spectrum : sp.Spectrum1D
        Simulated exposure with given SNR.

    """
    spectrum = spectrum / spectrum.flux.mean()
    
    p_coeff = np.polyfit(SNR[:,0]*10, SNR[:,1], 4, rcond=None, full=False, w=None, cov=False)
    
    pol = np.poly1d(p_coeff)
    
    spectrum = sp.Spectrum1D(
        spectral_axis = spectrum.spectral_axis,
        flux = spectrum.flux * pol(spectrum.spectral_axis.value)**2,
        uncertainty = astropy.nddata.StdDevUncertainty(
            np.sqrt(spectrum.flux * pol(spectrum.spectral_axis.value)**2)
            ),
        wcs = spectrum.wcs
        )
    
    return spectrum
#%%
def simulate_transit_parameters(system_parameters:pd.DataFrame,
                                exptime: int = 900,
                                overheads: int = 68
                                ) -> tuple[int, int, np.ndarray, np.ndarray]:
    """
    Simulate a transit parameters.

    Parameters
    ----------
    system_parameters : pd.DataFrame
        System parameters.
    exptime : int, optional
        Exposure time to use for the observation. The default is 900.
    overheads : int, optional
        Overheads to account for each observation. The default is 68.

    Returns
    -------
    num_out_transit : int
        Number of out of transit exposures.
    num_in_transit : int
        Number of in transit exposures.
    phase : np.ndarray
        Phase array for the simulated transit.
    velocity_planet : np.ndarray
        Velocity array of the planet for the simulated transit.

    """
    
    system_parameters = calculate_semimajor_axis(system_parameters)
    system_parameters = calculate_inclination(system_parameters)
    
    P = system_parameters['pl_orbper'] * u.day
    a = system_parameters['pl_orbsmax'] *u.au
    i = system_parameters['pl_orbincl'] *u.deg
    
    if np.isnan(a) and np.isfine(P):
        a = (con.G * (
                system_parameters['st_mass']*u.M_sun +
                system_parameters['pl_bmasse']*u.M_earth
                ) /
            np.pi**2 / 4
            )**(1/3)
    
    if np.isnan(i):
        pass
    
    from transits.calculate_transit_windows import define_baseline
    baseline_length = define_baseline(system_parameters['pl_trandur'])
    
    num_out_transit = round(baseline_length * 3600 / (exptime+overheads))
    if num_out_transit%2 == 1:
        num_out_transit += 1
    
    num_in_transit = round(system_parameters['pl_trandur'] * 3600 / (exptime+overheads))
    num_exposures = num_out_transit + num_in_transit
    
    exposure_phase = (exptime+overheads) / system_parameters['pl_orbper'] / 3600 / 24
    
    phase = np.linspace(-exposure_phase * num_exposures/2,
                        exposure_phase * num_exposures/2,
                        num_exposures)
    
    velocity_planet = ((2. * np.pi * a.to(u.m) /
                       (P.to(u.s))) * \
                       np.sin(2.*np.pi * u.rad * phase) * \
                       np.sin(np.radians(i))
                       )
    
    T1 = -(system_parameters['pl_trandur']/24 / system_parameters['pl_orbper'])/2
    T4 = -T1
    
    return num_out_transit, num_in_transit, phase, velocity_planet, T1, T4

#%% simulate_dataset
def simulate_dataset(
    mock_spectrum: sp.Spectrum1D,
    system_parameters: pd.DataFrame,
    instrument: str= 'ESPRESSO',
    exptime: int = 900,
    overheads: int = 68
    ) -> sp.SpectrumList:
    """
    Simulate a dataset for given transit

    Parameters
    ----------
    mock_spectrum : sp.Spectrum1D
        Spectrum from which to create the dataset.
    system_parameters : pd.DataFrame
        System parameters as extracted from NASA archive.

    Returns
    -------
    mock_dataset : sp.SpectrumList
        Mock dataset.

    """
    # Simulate a transit
    (num_out_transit,
    num_in_transit,
    phases,
    velocity_planet,
    T1,
    T4) = simulate_transit_parameters(system_parameters,
                                      exptime= exptime,
                                      overheads= overheads)
    
    mock_dataset = sp.SpectrumList()
    # Loop through velocity list
    for ind, (velocity, phase) in enumerate(zip(velocity_planet, phases)):
        if phase > T1 and phase < T4:
            transit = True
        else:
            transit = False
        
        simulated_flux = mock_spectrum.flux + np.random.normal(0, mock_spectrum.uncertainty.array, len(mock_spectrum.flux))
        simulated_uncertainty = mock_spectrum.uncertainty.array
        
        new_spectrum =  sp.Spectrum1D(
                        spectral_axis = mock_spectrum.spectral_axis,
                        flux = simulated_flux,
                        uncertainty = astropy.nddata.StdDevUncertainty(simulated_uncertainty),
                        mask = np.zeros_like(simulated_flux),
                        meta = {
                            'vel_pl': velocity,
                            'Transit_full': transit,
                            'normalization': True,
                            'Night_num': 1,
                            'RF': '',
                            'Night': 'Simulated night',
                            'S_N': 1
                            }
                        )
        
        new_spectrum = sm.normalize_spectrum(new_spectrum, quantile=.5, linfit=False)
        
        if transit:
            new_spectrum = inject_mock_planetary_signal(new_spectrum, system_parameters)
            new_spectrum.meta.update(
                    {
                    'vel_pl': velocity,
                    'Transit_full': transit,
                    'normalization': True,
                    'Night_num': 1,
                    'RF': '',
                    'Night': 'Simulated night',
                    'S_N': 1
                    }
                )

            
        mock_dataset.append(
            new_spectrum
            )
    return mock_dataset

#%% simulate_transmission_spectrum
def simulate_transmission_spectrum(mock_dataset: sp.SpectrumList) -> sp.Spectrum1D:
    """
    Simulate transmission spectrum.

    Parameters
    ----------
    mock_dataset : sp.SpectrumList
        Mock dataset to simulate transmission spectrum from.

    Returns
    -------
    transmission_spectrum : sp.Spectrum1D
        Simulated transmission spectrum.

    """
    shifted_data = sp.SpectrumList()
    # Calculate master out
    master_out = sm.calculate_master_list(mock_dataset,
                                          key = 'Transit_full',
                                          value = False,
                                          sn_type='quadratic'
                                          )
    
    for item in mock_dataset:
        # Divide by master out
        M_out_corrected = item.divide(master_out[0], handle_meta= 'first_found')
        
        # Shift to planetary rest frame
        shifted_data.append(
            sm.shift_spectrum_multiprocessing(M_out_corrected,
                                              shift_BERV = 0,
                                              shift_v_sys = 0,
                                              shift_v_star = 0,
                                              shift_v_planet = -1,
                                              shift_constant = 0)
            )
    # Calculate transmission spectrum
    transmission_spectrum = sm.calculate_master_list(shifted_data,
                                                     key = 'Transit_full',
                                                     value = True,
                                                     sn_type='quadratic')[0]
    
    return shifted_data, transmission_spectrum
#%% get_amplitudes_signal
def get_amplitudes_signal(system_parameters: pd.DataFrame):
    """
    Calculate strength of the transmission signal of sodium given planetary equilibrium temperature.

    Parameters
    ----------
    system_parameters : pd.DataFrame
        System parameters as extracted from NASA archive.

    Returns
    -------
    H : float
        Strength of the signal in number of scale heights.

    """
    logger.warning('This is obsolete way of calculating the signal strength')
    
    if system_parameters['pl_eqt'] > 2650:
        H = 20
        return H, H
    else:
        x = system_parameters['pl_eqt']
        H = (-13510/(x-2650)+36)
        H2 = (-13387/(x-2650)+20.60254066)
        return H, H2

#%%
def scale_to_excess_absorption(H_num: float,
                               H_size: float,
                               RsRp_ratio: float,
                               Rp : float
                               ):
    """
    Scale number of scale heights to excess absorption.

    Parameters
    ----------
    H_num : float
        Number of scale heights to scale.
    H_size : float
        Size of one atmospheric scale height in km.
    RsRp_ratio : float
        Star-to-Planet radius ratio.
    Rp : float
        Radius of the planet in km.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    nominator = ((H_num * H_size +1)**2 * RsRp_ratio**(-2)) - 1
    denominator = RsRp_ratio**(-1) - 1
    H_size /= Rp
    
    return 1 - ((H_num * H_size+1)**2*RsRp_ratio**(-2)-(1))/(-1+RsRp_ratio**(-2))


#%% inject_mock_planetary_signal
def inject_mock_planetary_signal(spectrum: sp.Spectrum1D,
                                 system_parameters: pd.DataFrame,
                                 trend:u.Quantity = 0*u.km/u.s) -> sp.Spectrum1D:
    """
    Inject mock planetary signal.

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Spectrum to which to inject the signal.
    system_parameters : pd.DataFrame
        System parameters as extracted from NASA archive.

    Returns
    -------
    spectrum : sp.Spectrum1D
        Spectrum with injected signal.

    """
    v = spectrum.meta['vel_pl']
    doppler_factor = (1+v/con.c).decompose()
    D2_mean = 5889.950*u.AA 
    D1_mean = 5895.924*u.AA
    D2_mean /= doppler_factor
    D1_mean /= doppler_factor
    
    H_num_D2, H_num_D1 = get_amplitudes_signal(system_parameters)
    amp_D2 = - scale_to_excess_absorption(H_num_D2,
                                        system_parameters['H'],
                                        (system_parameters['st_rad'] *u.R_sun/ 
                                         system_parameters['pl_radj'] /u.R_jup
                                         ).decompose(),
                                        (system_parameters['pl_radj']*u.R_jup).to(u.km).value
                                        )
    amp_D1 = - scale_to_excess_absorption(H_num_D1,
                                        system_parameters['H'],
                                        (system_parameters['st_rad'] *u.R_sun/ 
                                         system_parameters['pl_radj'] /u.R_jup
                                         ).decompose(),
                                        (system_parameters['pl_radj']*u.R_jup).to(u.km).value
                                        )
    
    sodium_doublet_model = (models.Gaussian1D(amplitude = amp_D2,
                                              mean = D2_mean,
                                              stddev =.1*u.AA) +
                            models.Gaussian1D(amplitude = amp_D1,
                                              mean = D1_mean,
                                              stddev =.1*u.AA) +
                            models.Const1D(amplitude = 1)
                            )
    new_spectrum = spectrum * sodium_doublet_model(spectrum.spectral_axis)
    
    return new_spectrum


def _inject_planetary_signal_petitRADtrans(
    planet_name: str,
    system_parameters,
    spectrum: sp.Spectrum1D,
    species: list) -> sp.Spectrum1D:
    import rats.parameters as para
    
    system_parameters = para.SystemParametersComposite()
    system_parameters.load_NASA_CompositeTable_values(planet_name=planet_name,
                                                      )
    import rats.modeling_CCF as ccf
    template = ccf._create_single_template(
        SystemParameters= system_parameters,
        spectral_axis= spectrum.spectral_axis,
        species= species
        )
    
    
    
    return

#%% full_simulation
def full_simulation(mock_spectrum: sp.Spectrum1D,
                    system_parameters:pd.DataFrame,
                    instrument: str = 'ESPRESSO',
                    exptime: float = 900*u.s,
                    overheads: float= 68*u.s) -> tuple[sp.SpectrumList, sp.Spectrum1D]:
    """
    Simulate the transit observation.

    Parameters
    ----------
    mock_spectrum : sp.Spectrum1D
        Mock stellar spectrum.
    system_parameters : pd.DataFrame
        System parameters as extracted from NASA archive.

    Returns
    -------
    mock_dataset : sp.SpectrumList
        Complete mock dataset.
    transmission_spectrum : sp.Spectrum1D
        Simulated resulting transmission spectrum.

    """
    mock_dataset = simulate_dataset(mock_spectrum, system_parameters, instrument, exptime, overheads)
    master_corrected_dataset, transmission_spectrum = simulate_transmission_spectrum(mock_dataset)
    return mock_dataset, master_corrected_dataset, transmission_spectrum
#%%
def simulate_all_planets(
    filtered_table: pd.DataFrame,
    save_directory: str,
    instrument: str = 'ESPRESSO',
    exptime: float = 900*u.s,
    overheads: float= 68*u.s,
    conditions: list = [False, True, False]
    ):

    for ind, row in filtered_table.iterrows():
        logger.info(str('Calculation of: '+row['pl_name']))
        while True:
            try:
                logger.info('Loading stellar spectrum')
                best_spectrum, average_spectrum, worst_spectrum = generate_mock_spectrum(row)
            except: # Sometimes the loading of PHOENIX spectra fails due to TimeoutError, retry on exception
                logger.warning('Received Error, retrying.')
                continue
            break
        

        if (np.isfinite(row['pl_bmasse']) and
            (np.isfinite(row['pl_imppar']) or np.isfinite(row['pl_orbincl']))
            ): # Skip planets without mass or impact parameter/inclination
            
            sr = sp.SpectralRegion(5876*u.AA, 5930*u.AA)
            best_spectrum = sm.extract_region(best_spectrum, sr)
            average_spectrum = sm.extract_region(average_spectrum, sr)
            worst_spectrum = sm.extract_region(worst_spectrum, sr)
            
            (best_dataset,
            best_Mout_corrected,
            best_transmission_spectrum) = full_simulation(best_spectrum,
                                                          row,
                                                          instrument,
                                                          exptime,
                                                          overheads)
            (average_dataset,
            average_Mout_corrected,
            average_transmission_spectrum ) = full_simulation(average_spectrum,
                                                              row,
                                                          instrument,
                                                          exptime,
                                                          overheads)
            (worst_dataset,
            worst_Mout_corrected,
            worst_transmission_spectrum) = full_simulation(worst_spectrum, row,
                                                          instrument,
                                                          exptime,
                                                          overheads)
        
            fig, ax, ax2 = plot_transmission_spectrum(average_transmission_spectrum, color = 'cyan')
        
            ax.set_title(row['pl_name'])
            os.makedirs(save_directory,
                        mode = 0o777,
                        exist_ok = True) 
            fig.savefig(save_directory + row['pl_name']+'.png')
            logger.info('Saving transmission spectrum simulation in: \n' + save_directory + row['pl_name']+'.png')
        else:
            logger.warning('Skipping planet: '+row['pl_name'])
        
    return

#%%
def double_gaussian(x:float,
                    mean_0:float,
                    amplitude_0:float,
                    stddev_0:float,
                    mean_1:float,
                    amplitude_1:float,
                    stddev_1:float,
                    ) -> float:
    """
    Double gaussian model to use for the fitting.

    Parameters
    ----------
    x : float
        The x parameter to use as the input variable for the model.
    mean_0 : float
        Mean of the first Gaussian line.
    amplitude_0 : float
        Amplitude of the first Gaussian line.
    stddev_0 : float
        Standard deviation of the first Gaussian line.
    mean_1 : float
        Mean of the second Gaussian line.
    amplitude_1 : float
        Amplitude of second Gaussian line.
    stddev_1 : float
        Standard deviation of the second Gaussian line.

    Returns
    -------
    model : float
        Value of the model.

    """
    G_0 = - amplitude_0 * np.exp(-((x-mean_0)**2)/(2*stddev_0**2))
    G_1 = - amplitude_1 * np.exp(-((x-mean_1)**2)/(stddev_1**2))
    return G_1 + G_0


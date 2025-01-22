# import specutils as sp
# import astropy
# import astropy.units as u
# import numpy as np
# import pandas as pd
# import os 


# # cwd = os.getcwd()
# # os.chdir('/media/chamaeleontis/Observatory_main/Code/')
# # from rats.utilities import default_logger_format
# # import rats.parameters as para
# # import rats.spectra_manipulation as sm
# # os.chdir(cwd)

# import calculations as calc
# import instruments as inst
# import plots as plot
# import logging


# logger = logging.getLogger(__name__)
# logger = default_logger_format(logger) 

# def _load_ETC_file(filename: str) -> np.ndarray:
#     """
#     Load the ETC files based on filename.

#     Parameters
#     ----------
#     filename : str
#         Filename location.

#     Returns
#     -------
#     SNR : np.ndarray
#         Array of [wavelength, expected SNR] as loaded from the file.

#     """
#     f = open(filename,'r')
    
#     SNR =[]
#     for line in f.readlines():
#         if line.startswith('#'):
#             continue
        
#         wavelength, snr = line.replace('\n','').split('\t')
        
#         SNR.append([float(wavelength), float(snr)])
    
#     return np.asarray(SNR)

# #%%
# def _find_stellar_type_ESPRESSO(T_eff: float) -> str:
#     """
#     Find stellar type based on effective temperature.

#     Parameters
#     ----------
#     T_eff : float
#         Effective temperature of the star in K.

#     Returns
#     -------
#     spectral_type : str
#         Spectral type closest to the true one that is available in ESPRESSO ETC.

#     """
#     if T_eff > 54000:
#         return 'O5'
#     elif T_eff > 43300:
#         return 'O9'
#     elif T_eff > 29200:
#         return 'B1'
#     elif T_eff > 23000:
#         return 'B3'
#     elif T_eff > 15000:
#         return 'B8'
#     elif T_eff > 11800:
#         return 'B9'
#     elif T_eff > 10000:
#         return 'A0'
#     elif T_eff > 9450:
#         return 'A1'
#     elif T_eff > 8480:
#         return 'F0'
#     elif T_eff > 5900:
#         return 'G0'
#     elif T_eff > 5560:
#         return 'G2'
#     elif T_eff > 4730:
#         return 'K2'
#     elif T_eff > 3865:
#         return 'K7'
#     else:
#         return 'M2'
    

# def _load_SNR_table_ESPRESSO(
#     TargetParameters: para.SystemParametersComposite,
#     instrument: str = 'ESPRESSO'
#     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
#     """
#     Load SNR tables based on the system parameters and instrument used.

#     Parameters
#     ----------
#     system_parameters : pd.DataFrame
#         System parameters as provided by NASA archive.
#     instrument : str, optional
#         Which instrument to use. The default is 'ESPRESSO'.

#     Returns
#     -------
#     best_scenario : np.ndarray
#         Array of [wavelength, expected SNR] for best case scenario.
#     average_scenario : np.ndarray
#         Array of [wavelength, expected SNR] for average (50%) case scenario.
#     worst_scenario : np.ndarray
#         Array of [wavelength, expected SNR] for worst case scenario.

#     """
#     if instrument != 'ESPRESSO' or 'ESPRESSO_4UT':
#         logger.warning('Loading SNR values from ESPRESSO, even through usage of different instrument. This will increase error in the resulting dataset, which however will not show in the result.')
    
#     spectral_type = _find_stellar_type_ESPRESSO(
#         float(TargetParameters.Star.temperature.data)
#         )
    
#     location_of_ETC_files = 'PTO/./ETC_{instrument}/{spectral_type}/'.format(instrument = instrument, spectral_type = spectral_type)
    
#     best_scenario = _load_ETC_file(location_of_ETC_files + '0.5arcsec_900s.etc')
#     average_scenario = _load_ETC_file(location_of_ETC_files + '0.8arcsec_900s.etc')
#     worst_scenario = _load_ETC_file(location_of_ETC_files + '1.3arcsec_900s.etc')
    
#     return best_scenario, average_scenario, worst_scenario
# #%% Scaling by instrument
# def _scale_instrument(spectrum: sp.Spectrum1D,
#                       from_instrument:str,
#                       to_instrument: str):
#     """
#     Scale the SNR between instrument, using the telescope mirror size. Instrument must be within the Instruments_size Enum class.

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Spectrum to scale
#     from_instrument : str
#         Original instrument, usually ESPRESSO, for which the code was implemented.
#     to_instrument : str
#         Target instrument, for which we want to simulate the dataset.
#     """
    
#     scale_factor = (inst.Instruments_size[to_instrument].value / inst.Instruments_size[from_instrument].value)**2 
    
#     new_spectrum = sp.Spectrum1D(
#         spectral_axis = spectrum.spectral_axis,
#         flux = spectrum.flux * scale_factor,
#         uncertainty = astropy.nddata.StdDevUncertainty(
#             np.sqrt(spectrum.flux * scale_factor)
#             ),
#         mask = np.zeros_like(spectrum.flux),
#         meta= {}
#         )
    
#     return new_spectrum

# #%% vactoair
# def vactoair(wlnm: np.ndarray):
#     """
#     Change from vacuum to air wavelengths.

#     Parameters
#     ----------
#     wlnm : np.ndarray
#         Wavelenght grid in nanometers (actually invariant towards units)
#     """
#     wlA = wlnm*10.0
#     s = 1e4/wlA
#     f = 1.0 + 5.792105e-2/(238.0185e0 - s**2) + 1.67917e-3/( 57.362e0 - s**2)
#     return(wlA/f/10.0)

# #%%
# def scale_estimated_SNR_with_Vmag(scenario: np.ndarray,
#                                   TargetParameters: para.SystemParametersComposite) -> np.ndarray:
#     """
#     Scale estimated SNR with V magnitude of the target.

#     Parameters
#     ----------
#     scenario : np.ndarray
#         Scenario as loaded by load_SNR_table function.
#     system_parameters : para.SystemParametersComposite
#         System parameters as loaded by NASA archive using the rats.para.SystemParametersComposite class.

#     Returns
#     -------
#     scenario : np.ndarray
#         Adjusted SNR of given scenario.

#     """
#     scenario[:,1] *= np.sqrt(10**(0.4*(10 - TargetParameters.Star.magnitudes.V.data)))
#     return scenario

# #%%
# def _scale_spectrum_with_SNR(
#     spectrum: sp.Spectrum1D,
#     SNR: np.ndarray) -> sp.Spectrum1D:
#     """
#     Scale spectrum according to estimated SNR.

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Model spectrum to scale.
#     SNR : np.ndarray, optional
#         Estimated SNR of single exposure as provided by _load_SNR_table.

#     Returns
#     -------
#     spectrum : sp.Spectrum1D
#         Simulated exposure with given SNR.

#     """
#     spectrum = spectrum / spectrum.flux.mean()
    
#     p_coeff = np.polyfit(SNR[:,0]*10, SNR[:,1], 4, rcond=None, full=False, w=None, cov=False)
    
#     pol = np.poly1d(p_coeff)
    
#     spectrum = sp.Spectrum1D(
#         spectral_axis = spectrum.spectral_axis,
#         flux = spectrum.flux * pol(spectrum.spectral_axis.value)**2,
#         uncertainty = astropy.nddata.StdDevUncertainty(
#             np.sqrt(spectrum.flux * pol(spectrum.spectral_axis.value)**2)
#             )
#         )
    
#     return spectrum

# #%%
# def generate_mock_spectrum(
#     TargetParameters: para.SystemParametersComposite,
#     instrument: str = 'ESPRESSO',
#     exptime: u.Quantity = 900*u.s,
#     ) -> sp.Spectrum1D:
#     """
#     Generate mock spectrum using the PHOENIX high resolution spectra.

#     Parameters
#     ----------
#     TargetParameters : para.SystemParametersComposite
#         Target parameters as loaded by rats.para.SystemParametersComposite class
#     instrument : str, optional
#         Which instrument to consider. The default is 'ESPRESSO'.

#     Returns
#     -------
#     best_scenarion_spectrum : sp.Spectrum1D
#         Mock spectrum interpolated to instrument wavelength grid assuming the best weather conditions.
#     average_scenarion_spectrum : sp.Spectrum1D
#         Mock spectrum interpolated to instrument wavelength grid assuming the average weather conditions.
#     worst_scenarion_spectrum : sp.Spectrum1D
#         Mock spectrum interpolated to instrument wavelength grid assuming the worst weather conditions.

#     """
#     stellar_spectrum = TargetParameters.Star.stellar_model()
    
#     new_stellar_spectrum = sm.interpolate2commonframe(
#         stellar_spectrum,
#         inst.Instruments[instrument]
#         )
    
#     best_scenario, average_scenario, worst_scenario = _load_SNR_table_ESPRESSO(
#         TargetParameters,
#         instrument = instrument
#         )
    
#     best_scenario = scale_estimated_SNR_with_Vmag(best_scenario, TargetParameters)
#     average_scenario = scale_estimated_SNR_with_Vmag(average_scenario, TargetParameters)
#     worst_scenario = scale_estimated_SNR_with_Vmag(worst_scenario, TargetParameters)
    
#     # FIXME: This is ESPRESSO specific, do a branch for other instrument. 
#     # Right now, this is fixed by _scale_instrument function, but should be unified
#     best_scenario_spectrum = _scale_spectrum_with_SNR(new_stellar_spectrum, best_scenario)
#     average_scenario_spectrum = _scale_spectrum_with_SNR(new_stellar_spectrum, average_scenario)
#     worst_scenario_spectrum = _scale_spectrum_with_SNR(new_stellar_spectrum, worst_scenario)
    
#     if instrument != 'ESPRESSO' or instrument != 'ESPRESSO_4UT':
#         best_scenario_spectrum = _scale_instrument(best_scenario_spectrum, 'ESPRESSO', instrument)
#         average_scenario_spectrum = _scale_instrument(average_scenario_spectrum, 'ESPRESSO', instrument)
#         worst_scenario_spectrum = _scale_instrument(worst_scenario_spectrum, 'ESPRESSO', instrument)
    
#     if exptime != 900*u.s:
#         best_scenario_spectrum = best_scenario_spectrum.multiply(exptime/900*u.s)
#         average_scenario_spectrum = average_scenario_spectrum.multiply(exptime/900*u.s)
#         worst_scenario_spectrum = worst_scenario_spectrum.multiply(exptime/900*u.s)
    
#     best_scenario_spectrum, average_scenario_spectrum, worst_scenario_spectrum = _add_mask([best_scenario_spectrum, average_scenario_spectrum, worst_scenario_spectrum])
    
#     return best_scenario_spectrum, average_scenario_spectrum, worst_scenario_spectrum


# def _add_mask(dataset: sp.SpectrumList) -> sp.SpectrumList:
#     """
#     Adds a mask to a dataset based on its flux element.

#     Parameters
#     ----------
#     dataset : sp.SpectrumList
#         Dataset to add mask in.

#     Returns
#     -------
#     dataset : sp.SpectrumList
#         Dataset with mask argument.
#     """
#     for spectrum in dataset:
#         spectrum.mask = np.isnan(spectrum.flux)
    
#     return dataset
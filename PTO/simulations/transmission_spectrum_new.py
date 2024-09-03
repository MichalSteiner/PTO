"""
_summary_
"""
#%% Importing libraries
import os
from re import A

import astropy
import astropy.units as u
import specutils as sp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
from copy import deepcopy
from dataclasses import dataclass 
import astropy
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors



# FIXME This is ugly!
import sys
sys.path.append(os.path.dirname(__file__))
sys.path.append('/media/chamaeleontis/Observatory_main/Code/observations_transits/PTO/PTO/transits/')
sys.path.append('/media/chamaeleontis/Observatory_main/Code/observations_transits/PTO/PTO/')

# FIXME Change this to package download.
cwd = os.getcwd()
os.chdir('/media/chamaeleontis/Observatory_main/Code/')
from rats.utilities import default_logger_format
import rats.spectra_manipulation as sm
import rats.parameters as para
# import rats.plots.spectra as ps
os.chdir(cwd)

import calculations as calc
import instruments as inst
import plots as plot
import calculate_transit_windows as ctw

logger = logging.getLogger(__name__)
logger = default_logger_format(logger) 

#%%
@dataclass
class SingleConditionResult:
    CCF_results: list | None = None
    data_PRF: list | None = None
    transmission_spectrum_combined: sp.Spectrum1D | None = None
    model: sp.Spectrum1D | None = None
    pass

@dataclass
class ResultsSinglePlanet():
    best: SingleConditionResult | None = None
    average: SingleConditionResult | None = None
    worst: SingleConditionResult | None = None
    TargetParameters: para.SystemParametersComposite | None = None
    pass


def _spectrum_flattening(spectrum: sp.Spectrum1D) -> sp.Spectrum1D:
    """
    Flattens the spectrum by rolling window quantile. Useful for flattening models.

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Spectrum to flatten

    Returns
    -------
    new_spectrum : sp.Spectrum1D
        Flattened spectrum
    """
    flux_array = pd.Series(spectrum.flux.value)
    flattened_flux = flux_array.rolling(250, min_periods=1, center=True).quantile(.95)
    
    flattening_spectrum = sp.Spectrum1D(
        spectral_axis= spectrum.spectral_axis,
        flux = flattened_flux.values* spectrum.flux.unit,
        uncertainty = astropy.nddata.StdDevUncertainty(np.zeros_like(spectrum.spectral_axis.value)),
        meta = {},
        mask = np.isnan(flattened_flux)
    )
    
    # FIXME This is dumb!
    spectrum = sm.replace_flux_units_transmission([spectrum], u.dimensionless_unscaled)[0]
    
    new_spectrum = spectrum.divide(flattening_spectrum, handle_meta = 'first_found')
    
    return new_spectrum

def flatten_spectrum_list(spectrum_list: sp.SpectrumList) -> sp.SpectrumList:
    """
    Wrapper to flatten a list of spectra

    Parameters
    ----------
    spectrum_list : sp.SpectrumList
        Spectrum list to flatten

    Returns
    -------
    new_spectrum_list : sp.SpectrumList
        Flattened list
    """
    new_spectrum_list = sp.SpectrumList()
    for item in spectrum_list:
        new_spectrum_list.append(
            _spectrum_flattening(item)
        )

    return new_spectrum_list 

#%%
def simulate_transmission_spectra(target_names: list,
                                  species_list: list,
                                  instrument: str = 'ESPRESSO',
                                  exptime: u.Quantity | None = None,
                                  overheads: u.Quantity = 68*u.s,
                                  best_calculate: bool = True,
                                  average_calculate: bool = True,
                                  worst_calculate: bool = True,
                                  number_of_transits: int = 1,
                                  ) -> dict:
    
    
    # DOCUMENTME
    # TODO
    logger.info('Starting simulation of transmission dataset')
    
    all_planet_results = {}
    
    for target in target_names:
        logger.info(f'    Simulation of {target}')
        
        best_results, average_results, worst_results, TargetParameters = _simulate_transmission_single_planet(
            target= target,
            species_list= species_list,
            instrument= instrument,
            exptime= exptime,
            overheads= overheads,
            best_calculate= best_calculate,
            average_calculate= average_calculate,
            worst_calculate= worst_calculate,
            number_of_transits= number_of_transits,
        )
        
        PlanetResults = ResultsSinglePlanet(
            best= best_results,
            average= average_results,
            worst= worst_results,
            TargetParameters= TargetParameters
        )
        logger.info(f'    Simulation of {target} complete!')
        all_planet_results[target] = PlanetResults
        
        # plot_CCF(PlanetResults.average,
        #          species= 'Fe',
        #          Target_parameters= PlanetResults.TargetParameters)
        
        # plot_TS(transmission_spectrum= PlanetResults.average.transmission_spectrum_combined,
        #         model= PlanetResults.average.model)
        
        # plot_sodium_resolved(PlanetResults.average.data_PRF)
    return all_planet_results

def _load_system_parameters(target_name: str) -> para.SystemParametersComposite:
    """
    Load system parameters using the SystemParametersComposite class from rats package.

    Parameters
    ----------
    target_name : str
        Target name as defined by NASA Exoplanet archive.

    Returns
    -------
    Target_parameters : para.SystemParametersComposite
        Target parameters loaded through the TAP service of NASA composite table.
    """
    
    Target_parameters = para.SystemParametersComposite(
        filename= os.getcwd() + f'/saved_data/system_parameters_{target_name}.pkl'
        )
    Target_parameters.load_NASA_CompositeTable_values(
        planet_name = target_name,
        force_load= True
        )
    
    return Target_parameters

def _simulate_transmission_single_planet(
    target: str,
    species_list: list,
    instrument: str = 'ESPRESSO',
    exptime: u.Quantity | None = None,
    overheads: u.Quantity = 68*u.s,
    best_calculate: bool = True,
    average_calculate: bool = True,
    worst_calculate: bool = True,
    number_of_transits: int = 1,
    ) -> [SingleConditionResult, SingleConditionResult, SingleConditionResult, para.SystemParametersComposite]:
    """
    Simulates a transmission dataset for single planet

    Parameters
    ----------
    target : str
        Planet name as available on NASA Exoplanet archive.
    species_list : list
        List of species to consider for petitRADTRANS model.
    instrument : str, optional
        Instrument to create the dataset for, by default 'ESPRESSO'
    exptime : u.Quantity | None, optional
        Exposure time to consider, by default None. If None, will try to optimize the result [CURRENTLY NOT IMPLEMENTED].
    overheads : u.Quantity, optional
        Overheads length, by default 68*u.s [ESPRESSO 1UT].
    best_calculate : bool, optional
        If True, will calculate the dataset for best seeing conditions, else will skip. By default True.
    average_calculate : bool, optional
        If True, will calculate the dataset for average seeing conditions, else will skip. By default True.
    worst_calculate : bool, optional
        If True, will calculate the dataset for worst seeing conditions, else will skip. By default True.
    number_of_transits : int, optional
        Number of transits to stack, by default 1.

    Returns
    -------
    best_results : SingleConditionResult
        Class holding the best weather condition results for given planet
    average_results : SingleConditionResult
        Class holding the average weather condition results for given planet
    worst_results : SingleConditionResult
        Class holding the worst weather condition results for given planet
    Target_parameters : para.SystemParametersComposite
        Target parameters as extracted by rats pipeline.
    """
    Target_parameters = _load_system_parameters(target_name= target)
    
    best_dataset, average_dataset, worst_dataset, model = _create_dataset(
        Target_parameters= Target_parameters,
        instrument= instrument,
        exptime= exptime,
        overheads= overheads,
        species= species_list,
        number_of_transits= number_of_transits)
    
    best = SingleConditionResult()
    average = SingleConditionResult()
    worst = SingleConditionResult()
    
    if best_calculate:
        best_CCF_result, best_data_PRF, best_transmission_spectrum = _run_trasmission_spectroscopy_pipeline(
            dataset= best_dataset,
            Target_parameters= Target_parameters,
            species= species_list
            )
        best = SingleConditionResult(
            CCF_results = best_CCF_result,
            data_PRF= best_data_PRF,
            transmission_spectrum_combined= best_transmission_spectrum,
            model= model
            )
        
    if average_calculate:
        average_CCF_result, average_data_PRF, average_transmission_spectrum = _run_trasmission_spectroscopy_pipeline(
            dataset= average_dataset,
            Target_parameters= Target_parameters,
            species= species_list
            )
        average = SingleConditionResult(
            CCF_results = average_CCF_result,
            data_PRF= average_data_PRF,
            transmission_spectrum_combined= average_transmission_spectrum,
            model= model
            )
        
    if worst_calculate:
        worst_CCF_result, worst_data_PRF, worst_transmission_spectrum = _run_trasmission_spectroscopy_pipeline(
            dataset= worst_dataset,
            Target_parameters= Target_parameters,
            species= species_list
            )
        worst = SingleConditionResult(
            CCF_results = worst_CCF_result,
            data_PRF= worst_data_PRF,
            transmission_spectrum_combined= worst_transmission_spectrum,
            model= model
            )
    
    return best, average, worst, Target_parameters

def _create_dataset(Target_parameters: para.SystemParametersComposite,
                    instrument = 'ESPRESSO',
                    exptime: u.Quantity | None = None,
                    overheads: u.Quantity = 68*u.s,
                    species: list = [],
                    number_of_transits: int = 1,
                    ) -> [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList, sp.Spectrum1D]:
    """
    Create a simulated dataset for single transit observation for given planet.

    Parameters
    ----------
    Target_parameters : para.SystemParametersComposite
        Target parameters, as loaded by rats package.
    instrument : str, optional
        Which instrument to create dataset for, by default 'ESPRESSO'. Other datasets lack several functions, like proper resolution and wavelength grid. It can be easily expanded on though.
    exptime : u.Quantity | None, optional
        Exposure time for single exposure, by default None. If None, default value is taken.
            To be implemented: Estimate optimal exposure through several simulations.
    overheads : u.Quantity, optional
        Overheads to consider, by default 68*u.s. These are the corresponding overheads for 1-UT ESPRESSO dataset with 2x1 binning mode.
    number_of_transits : int, optional
        Number of transits to stack, by default 1.

    Returns
    -------
    best_dataset, average_dataset, worst_dataset : [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList]
        Set of three datasets for different seeing.
    """
    
    if exptime is None:
        exptime = 900*u.s # TODO do more optimal selection
    
    baseline_length = ctw.define_baseline(Target_parameters.Ephemeris.transit_length_partial.data)
    transit_length = Target_parameters.Ephemeris.transit_length_partial.data
    
    observation_length = (baseline_length + transit_length) * u.hour
    period = Target_parameters.Ephemeris.period.data * u.day
    
    phases = np.linspace((-(observation_length/period).decompose()/ 2).value,
                         (+(observation_length/period).decompose()/ 2).value,
                         round(+(observation_length/(exptime+overheads)).decompose().value)
                         )
    
    
    best_scenario, average_scenario, worst_scenario = calc.generate_mock_spectrum(
        TargetParameters= Target_parameters,
        instrument= instrument,
        exptime= exptime
        )
    best_dataset, average_dataset, worst_dataset = _create_spectrum_list(
        phases= phases,
        Target_parameters= Target_parameters,
        best_scenario= best_scenario,
        average_scenario= average_scenario,
        worst_scenario= worst_scenario,
        number_of_transits = number_of_transits
        )
    
    best_dataset = sm.normalize_list(best_dataset,
                                     force_multiprocessing=False)
    average_dataset = sm.normalize_list(average_dataset,
                                     force_multiprocessing=False)
    worst_dataset = sm.normalize_list(worst_dataset,
                                     force_multiprocessing=False)
    
    
    best_dataset, average_dataset, worst_dataset, model = _inject_signal(
        best_dataset= best_dataset,
        average_dataset= average_dataset,
        worst_dataset= worst_dataset,
        Target_parameters= Target_parameters,
        species= species)

    return best_dataset, average_dataset, worst_dataset, model

def _create_spectrum_list(
    phases: np.ndarray,
    Target_parameters: para.SystemParametersComposite,
    best_scenario: sp.Spectrum1D,
    average_scenario: sp.Spectrum1D,
    worst_scenario: sp.Spectrum1D,
    number_of_transits: int = 1) -> [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList]:
    """
    Expands a template simulation spectrum into a dataset, including relevant phase and velocities data.

    Parameters
    ----------
    phases : np.ndarray
        Array of phases for which to generate the data.
    Target_parameters : para.SystemParametersComposite
        Target parameters as loaded by rats package.
    best_scenario : sp.Spectrum1D
        Best scenario simulation spectrum.
    average_scenario : sp.Spectrum1D
        Average scenario simulation spectrum.
    worst_scenario : sp.Spectrum1D
        Worst scenario simulation spectrum.
    number_of_transits : int
        Number of transits to stack together.

    Returns
    -------
    best_dataset, average_dataset, worst_dataset : [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList]
        Set of three datasets for different seeing conditions.
    """
    
    best_dataset, average_dataset, worst_dataset = sp.SpectrumList(), sp.SpectrumList(), sp.SpectrumList() 
    
    for ind, phase in enumerate(phases):
        general_meta = {
            'Night_num': 1,
            'Night': 'Simulated night',
            'Spec_num': ind,
            'normalization': True,
            'RF': 'Star',
            }
        
        best_dataset.append(_generate_similar_spectrum(best_scenario, number_of_transits))
        average_dataset.append(_generate_similar_spectrum(average_scenario, number_of_transits))
        worst_dataset.append(_generate_similar_spectrum(worst_scenario, number_of_transits))
        
        best_dataset[ind].meta['Phase'] = phase
        average_dataset[ind].meta['Phase'] = phase
        worst_dataset[ind].meta['Phase'] = phase
        
        best_dataset[ind].meta.update(general_meta)
        average_dataset[ind].meta.update(general_meta)
        worst_dataset[ind].meta.update(general_meta)
        
        if (phase < Target_parameters.Ephemeris.contact_T1.data or
            phase > Target_parameters.Ephemeris.contact_T4.data):
            best_dataset[ind].meta['Transit_full'] = False
            average_dataset[ind].meta['Transit_full'] = False
            worst_dataset[ind].meta['Transit_full'] = False
        else:
            best_dataset[ind].meta['Transit_full'] = True
            average_dataset[ind].meta['Transit_full'] = True
            worst_dataset[ind].meta['Transit_full'] = True
    
    Target_parameters.calculate_velocities_list(best_dataset)
    Target_parameters.calculate_velocities_list(average_dataset)
    Target_parameters.calculate_velocities_list(worst_dataset)
    logger.info(f'        Created simulated dataset for {Target_parameters.Planet.name}')
    return best_dataset, average_dataset, worst_dataset

def _generate_similar_spectrum(spectrum: sp.Spectrum1D,
                               number_of_transits: int = 1) -> sp.Spectrum1D:
    """
    Generates a similar spectrum given a model spectrum. This is done by adding randomly sampling normal distribution with standard deviation = uncertainty of the spectrum. This sample is then added to the flux array. Uncertainty of the new spectrum is then the square-root of the new flux.

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Model simulated spectrum.
    number_of_transits : int
        Number of transits to stack together, by default 1.

    Returns
    -------
    new_spectrum : sp.Spectrum1D
        New realization of the model spectrum assuming a given noise.
    """
    new_flux = (spectrum.flux.value
                + (np.random.normal(0, spectrum.uncertainty.array, len(spectrum.flux))/ np.sqrt(number_of_transits))
                )*u.dimensionless_unscaled
    
    new_spectrum = sp.Spectrum1D(
        spectral_axis = spectrum.spectral_axis,
        flux= new_flux,
        uncertainty= astropy.nddata.StdDevUncertainty(np.sqrt(new_flux)),
        mask= np.isnan(new_flux),
    )
    return new_spectrum
    


def _model_signal(Target_parameters: para.SystemParametersComposite,
                  spectral_axis: sp.SpectralAxis,
                  species: list,) -> sp.Spectrum1D:
    """
    Model a signal to inject into the dataset.

    Parameters
    ----------
    Target_parameters : para.SystemParametersComposite
        Target parameters as loaded by rats package
    spectral_axis : sp.SpectralAxis
        Spectral axis on which to consider the model
    species : list
        Species list to consider.

    Returns
    -------
    model : sp.Spectrum1D
        Model spectrum that can be injected in the data.
    """
    import rats.modeling_CCF
    model = rats.modeling_CCF._create_single_template(
        SystemParameters= Target_parameters,
        spectral_axis= spectral_axis,
        species= species,
        MMW_value= 2.33,
        )
    model = model.add(Target_parameters.Planet.radius.data * u.R_jup,
                      handle_meta= 'first_found').divide(
                          Target_parameters.Planet.radius.data * u.R_jup,
                          handle_meta= 'first_found')
    equivalency_transmission, F_lam, R, R_plam, delta_lam, H_num = sm.custom_transmission_units(Target_parameters)
    model = sm.replace_flux_units_transmission([model], R_plam)[0]
    
    model = model.new_flux_unit(F_lam, equivalencies=equivalency_transmission)
    model = _spectrum_flattening(model)
    
    model = sp.Spectrum1D(
        spectral_axis = model.spectral_axis,
        flux = np.where(model.flux.value < (1- (0.1*(1-np.nanmin(model.flux.value)))), 1 - model.flux.value, 0)  * u.dimensionless_unscaled,
        uncertainty = astropy.nddata.StdDevUncertainty(model.uncertainty.array),
        meta = model.meta,
        mask = model.mask
    )
    
    
    return model

def _inject_signal(best_dataset: sp.SpectrumList,
                   average_dataset: sp.SpectrumList,
                   worst_dataset: sp.SpectrumList,
                   Target_parameters: para.SystemParametersComposite,
                   species: list,
                   ) -> [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList]:
    """
    Injects a planetary atmospheric model inside the datasets for given species.

    Parameters
    ----------
    best_dataset : sp.SpectrumList
        Best simulated dataset.
    average_dataset : sp.SpectrumList
        Average simulated dataset.
    worst_dataset : sp.SpectrumList
        Worst simulated dataset.
    Target_parameters : para.SystemParametersComposite
        Target parameters as loaded by rats package
    species : list
        Species list. Must correspond to petitRADtrans list name.

    Returns
    -------
    best_dataset_injected, average_dataset_injected, worst_dataset_injected : [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList]
        Set of three datasets for different seeing conditions with injected planetary signal.
    """
    
    
    logger.info(f'        Creating template model for species list: {species}')
    for ind, specimen in enumerate(species):
        specimen_model = _model_signal(
            Target_parameters= Target_parameters,
            spectral_axis= best_dataset[0].spectral_axis,
            species= [specimen])

        if specimen.startswith('Na'):
            logger.warning('Found sodium as requested species. petitRADTRANS shows significant redshift of ~2.5km/s in its model, so it is manually shifted back to nominal value. Some small effect might be still visible.')
            specimen_model = sm._shift_spectrum(specimen_model,
                                                velocities= [2.5*u.km/u.s])
        if ind == 0:
            model = specimen_model
        if ind != 0:
            model = sp.Spectrum1D(
                spectral_axis = model.spectral_axis,
                flux = (model.flux.value + (specimen_model.flux.value * 5) *u.dimensionless_unscaled),
                uncertainty = astropy.nddata.StdDevUncertainty(np.zeros_like(model.spectral_axis.value)),
                mask = np.zeros_like(model.spectral_axis, dtype=bool),
                meta= {}
            )
        
    best_dataset_injected, average_dataset_injected, worst_dataset_injected = sp.SpectrumList(), sp.SpectrumList(), sp.SpectrumList()
    
    
    for best, average, worst in zip(best_dataset, average_dataset, worst_dataset):
        shifted_model = sm._shift_spectrum(
            model,
            velocities= [-best.meta['velocity_planet'].data * best.meta['velocity_planet'].unit,
                         +best.meta['velocity_star'].data * best.meta['velocity_star'].unit
                         ]
            )
        if best.meta['Transit_full']:
            shifted_model = shifted_model.multiply(-1*u.dimensionless_unscaled, handle_meta = 'first_found'
                                                   ).add(1*u.dimensionless_unscaled, handle_meta = 'first_found')
            best_dataset_injected.append(
                best.multiply(shifted_model, handle_meta='first_found')
                )
            average_dataset_injected.append(
                average.multiply(shifted_model, handle_meta='first_found')
                )
            worst_dataset_injected.append(
                worst.multiply(shifted_model, handle_meta='first_found')
                )
        else:
            best_dataset_injected.append(best)
            average_dataset_injected.append(average)
            worst_dataset_injected.append(worst)
            
    logger.info(f'        Finished injecting signal in simulated dataset for {Target_parameters.Planet.name}')
    return best_dataset_injected, average_dataset_injected, worst_dataset_injected, model

def _run_trasmission_spectroscopy_pipeline(dataset: sp.SpectrumList,
                                           Target_parameters: para.SystemParametersComposite,
                                           species: list) -> [dict, sp.Spectrum1D]:
    """
    Run transmission spectroscopy pipeline.

    Parameters
    ----------
    dataset : sp.SpectrumList
        Original dataset to calculate transmission spectrum from.
    Target_parameters : para.SystemParametersComposite
        Target parameters as loaded by rats package

    Returns
    -------
    result_CCF : dict
        Resulting CCF between simulated dataset and each of the requested species.
    transmission_spectrum : sp.Spectrum1D
        Average transmission spectrum calculated over the dataset.
    """
    dataset = calc._add_mask(dataset=dataset)
    master_out = sm.calculate_master_list(spectrum_list=dataset,
                                          key= 'Transit_full',
                                          value= False,
                                          )
    
    data_out_corrected = sm.spec_list_master_correct(spec_list=dataset,
                                                 master=master_out,
                                                 pkl_name = 'star_corrected_noRM.pkl'
                                                 )
    data_PRF = sm.shift_list(spectrum_list=data_out_corrected,
                                 shift_BERV=0,
                                 shift_v_sys = 0,
                                 shift_v_star = -1,
                                 shift_v_planet = 1,
                                 force_multiprocessing= False,
                                 pkl_name = 'data_PRF.pkl'
                                 )
    
    transmission_spectrum = sm.calculate_master_list(spectrum_list=data_PRF,
                                                 key = 'Transit_full',
                                                 value =True,
                                                 sn_type= None,
                                                 pkl_name = 'transmission_spectrum.pkl'
                                                 )
    import rats.spectra_cross_correlate as ccf
    
    result_CCF = {}
    data_PRF = sm.sigma_clipping_list(spectrum_list=data_PRF, num_of_sigma=5)
    
    new_data_PRF = sp.SpectrumList([item.subtract(1*u.dimensionless_unscaled, handle_meta= 'first_found').multiply(-1*u.dimensionless_unscaled, handle_meta = 'first_found') for item in data_PRF])
    
    for specimen in species:
        if specimen.startswith('Na'):
            logger.warning('Sodium detected as requested species, however petitRADTRANS sodium model is not precise enough, skipping. The sodium doublet is still visible.')
            continue
        
        test_model = _model_signal(
            Target_parameters=Target_parameters,
            spectral_axis=new_data_PRF[0].spectral_axis,
            species= [specimen]
        )
        
        CCF_model = ccf.cross_correlate_list(
            spectrum_list=new_data_PRF,
            model=test_model,
            sn_type= 'quadratic_error',
            velocities = np.arange(-50, 50.5, 0.5) * u.km/u.s
            )
        result_CCF[specimen] = CCF_model
    return result_CCF, data_PRF, transmission_spectrum[0]

# import rats.plots.plot_spectra_new as ps

def plot_sodium_resolved(average_data_PRF,
                         Target_parameters,
                         fig: plt.Figure | None = None,
                         ax_first: plt.Axes | None = None,
                         ax_second: plt.Axes | None = None):

    if fig is None:
        fig, ax = plt.subplots(1, figsize= (6*2,4*2))
    
    
    test = np.linspace(0,180,1)

    t=30
    cut_data_PRF = sm.extract_region_in_list(average_data_PRF, sp.SpectralRegion(5888*u.AA,5898*u.AA))
    spectrum_list = cut_data_PRF
    phase = [item.meta['Phase'] for item in spectrum_list]
    flux = np.asarray([(item.flux.value) for item in spectrum_list])
    
    cmap = sns.diverging_palette(90, 270, s=75, l=50, center='light', as_cmap=True)
    norm = colors.CenteredNorm(vcenter=1, halfrange=0.01)
    
    for ax in [ax_first, ax_second]:
        cax = ax.pcolormesh(spectrum_list[0].spectral_axis.value, phase, flux, cmap= cmap, norm=norm)
        y_lowest, y_highest = ax.get_ylim()
        ax.vlines(5889.950,y_lowest, Target_parameters.Ephemeris.contact_T1.data,color='darkred',ls='--')
        ax.vlines(5889.950,Target_parameters.Ephemeris.contact_T4.data,y_highest,color='darkred',ls='--')
        ax.vlines(5895.924,y_lowest, Target_parameters.Ephemeris.contact_T1.data,color='darkred',ls='--')
        ax.vlines(5895.924,Target_parameters.Ephemeris.contact_T4.data,y_highest,color='darkred',ls='--')
        Target_parameters.plot_contact_points(ax)
    cbar = fig.colorbar(cax, ax=[ax_first, ax_second], orientation='horizontal', location='top')
    cbar.set_label('<-Absorption (yellow) - $F_{in}/F_{out}$ - "Emission" (blue)->')
    ax_first.set_xlim(5889,5891)
    ax_second.set_xlim(5895,5897)
    
    ax_first.spines['right'].set_visible(False)
    ax_second.spines['left'].set_visible(False)
    ax_first.yaxis.tick_left()
    ax_second.yaxis.tick_right()
    ax_first.tick_params(axis='x', rotation=45)
    ax_second.tick_params(axis='x', rotation=45)
    ax_first.set_ylabel('Phase [1]')
    return



def plot_TS(transmission_spectrum: sp.Spectrum1D,
            model: sp.Spectrum1D | None = None,
            fig: plt.Figure | None = None,
            ax_first: plt.Axes | None = None,
            ax_second: plt.Axes | None = None,
            ):
    
    
    if fig is None:
        fig, ax = plt.subplots(1, figsize= (12,4))
        
    
    transmission_spectrum = sm.extract_region_in_spectrum(transmission_spectrum,
                                                          sp.SpectralRegion(5886*u.AA,5900*u.AA))
    model = sm.extract_region_in_spectrum(model,
                                          sp.SpectralRegion(5886*u.AA,5900*u.AA))
    
    for ax in [ax_first, ax_second]:
        ax.axvline(5889.950,color='darkred',ls='--', label='_nolegend_')
        ax.axvline(5895.924,color='darkred',ls='--', label='_nolegend_')
        ps.PlotSingleSpectrum.plot_spectrum(
            ax = ax,
            spectrum = transmission_spectrum,
            binning_factor= 10,
            color_spectrum = sns.color_palette("pastel")[0],
            color_bin = sns.color_palette("dark")[0],
            )
    
    ax_first.set_xlim(5889,5891)
    ax_second.set_xlim(5895,5897)
    
    ax_first.spines['right'].set_visible(False)
    ax_second.spines['left'].set_visible(False)
    ax_first.yaxis.tick_left()
    ax_second.yaxis.tick_right()
    ax_first.tick_params(axis='x', rotation=45)
    ax_second.tick_params(axis='x', rotation=45)

    # ax.set_xlabel('Wavelength [$\AA$]')
    ax_first.set_ylabel('Excess absorption [1]')
    return

def plot_CCF(Results: SingleConditionResult,
             species: str,
             Target_parameters: para.SystemParametersComposite,
             fig: plt.Figure | None = None,
             ax: plt.Axes | None = None):
    # TODO

    
    if fig is None:
        fig, ax = plt.subplots(1, figsize=(6*2,4*2))
    
    x= Results.CCF_results[species][0].spectral_axis.value  
    
    ind = np.argwhere(abs(x) > 25)
    noise_level = np.std(np.asarray([item.flux[ind].value for item in Results.CCF_results[species]]))* 15 # The factor of thirty is emprical comparison between 1-transit data and simulations. It comes from the fact that we assume ONLY white noise, but data are commonly contaminated by systematics and red noise as well. Factor of 15 works for 4-UT.
    
    y = np.asarray([item.meta['Phase'] for item in Results.CCF_results[species]])
    z = np.asarray([(item.flux.value)/noise_level for item in Results.CCF_results[species]])
    z = z + np.random.normal(0, 1, z.shape)
    
    z = z[:]
    
    
    
    norm = colors.CenteredNorm(vcenter= 0, halfrange= z.max())
    cmap = sns.diverging_palette(270, 90, s=75, l=50, center='light', as_cmap=True)
    cax = ax.pcolormesh(x,y,z,
                        cmap=cmap,
                        norm=norm)
    cbar = fig.colorbar(cax, ax=[ax], orientation='horizontal', location='top')
    cbar.set_label('<-"Emission" (blue) - SNR [1] -Absorption (yellow)->')
    y_lowest, y_highest = ax.get_ylim()
    
    ax.vlines(0,y_lowest, Target_parameters.Ephemeris.contact_T1.data,color='darkred',ls='--')
    ax.vlines(0,Target_parameters.Ephemeris.contact_T4.data,y_highest,color='darkred',ls='--')
    Target_parameters.plot_contact_points(ax)
    
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel('Phase [1]')
    
    return

def _extract_signal():
    return

def _fit_signal():
    test= (astropy.modeling.models.Gaussian1D(amplitude= -0.00449,
                                              mean=5889.950*u.AA,
                                              stddev=0.1*u.AA
                                              ) + 
           astropy.modeling.models.Gaussian1D(amplitude= -0.00385,
                                              mean=5895.924*u.AA,
                                              stddev=0.1*u.AA) + 
           astropy.modeling.models.Const1D(amplitude=1)
           )
    
    return

def _multiplot_full(PlanetResults):
    import matplotlib.gridspec as gridspec
    # fig, axs = plt.subplots(1, 5, squeeze= False, figsize=(6.6*3, 2*3))
    
    fig = plt.figure(layout='compressed', figsize=(6.6*3, 2*3))
    subfig_l, subfig_m, subfig_r = fig.subfigures(1, 3)
    
    gs = gridspec.GridSpec(9, 2, wspace= 0.10)
    ax = subfig_l.add_subplot(gs[0:8, 0])
    ax_second = subfig_l.add_subplot(gs[0:8, 1], sharey= ax)
    ax2 = subfig_m.add_subplot(gs[0:8, 0])
    ax2_second = subfig_m.add_subplot(gs[0:8, 1], sharey= ax2)
    ax3 = subfig_r.add_subplot(gs[0:8, :])
    
    subfig_l.suptitle('Simulated transmission spectrum, sodium doublet \n combined transmission spectrum')
    subfig_m.suptitle('Simulated transmission spectrum, sodium doublet \n temporally resolved')
    subfig_r.suptitle('Injection-Recovery of Iron from the simulated dataset')
    subfig_l.supxlabel('Wavelength $[\AA]$')
    subfig_m.supxlabel('Wavelength $[\AA]$')
    
    plot_TS(transmission_spectrum= PlanetResults.average.transmission_spectrum_combined,
            model= PlanetResults.average.model,
            fig= fig,
            ax_first= ax,
            ax_second= ax_second)
        
    plot_sodium_resolved(PlanetResults.average.data_PRF,
                         PlanetResults.TargetParameters,
                         fig= fig,
                         ax_first= ax2,
                         ax_second= ax2_second)
    
    plot_CCF(PlanetResults.average,
             species= 'Fe',
             Target_parameters= PlanetResults.TargetParameters,
             fig= fig,
             ax = ax3)
    return

def _multiplot_divided():
    import matplotlib.gridspec as gridspec
    # fig, axs = plt.subplots(1, 5, squeeze= False, figsize=(6.6*3, 2*3))
    
    fig = plt.figure(layout='compressed', figsize=(6.6*3, 2*3))
    subfig_l, subfig_m = fig.subfigures(1, 2)
    
    gs = gridspec.GridSpec(9, 2, wspace= 0.05)
    ax = subfig_l.add_subplot(gs[0:8, 0])
    ax_second = subfig_l.add_subplot(gs[0:8, 1], sharey= ax)
    ax2 = subfig_m.add_subplot(gs[0:8, 0])
    ax2_second = subfig_m.add_subplot(gs[0:8, 1], sharey= ax2)
    
    subfig_l.suptitle('Simulated transmission spectrum, sodium doublet \n combined transmission spectrum')
    subfig_m.suptitle('Simulated transmission spectrum, sodium doublet \n temporally resolved')
    subfig_l.supxlabel('Wavelength $[\AA]$')
    subfig_m.supxlabel('Wavelength $[\AA]$')
    
    plot_TS(transmission_spectrum= PlanetResults.best.transmission_spectrum_combined,
            model= PlanetResults.best.model,
            fig= fig,
            ax_first= ax,
            ax_second= ax_second)
    
    plot_sodium_resolved(PlanetResults.best.data_PRF,
                         PlanetResults.TargetParameters,
                         fig= fig,
                         ax_first= ax2,
                         ax_second= ax2_second)
    
    fig, ax = plt.subplots(1, figsize=(2.2*3, 2*3))
    plot_CCF(PlanetResults.average,
             species= 'Fe',
             Target_parameters= PlanetResults.TargetParameters,
             fig= fig,
             ax = ax)
    
    
    return

if __name__ == '__main__':
    target_names = [
        # 'WASP-193 b',
        # 'TOI-615 b',
        # 'TOI-622 b',
        # 'TOI-2641 b',
        # 'TOI-132 b',
        # 'WASP-76 b',
        # 'WASP-121 b',
        'WASP-78 b'
        ]
    species_list= [
        'Na_allard_new',
        'Fe',
        # 'Cr'
        ]
    instrument = 'ESPRESSO'
    exptime= 900*u.s
    overheads = 68*u.s
    
    # instrument = 'ESPRESSO_4UT'
    # overheads = 41*u.s
    # exptime= 900*u.s
    simulate_transmission_spectra(target_names= target_names,
                                  species_list= species_list,
                                  instrument = instrument,
                                  exptime = exptime,
                                  overheads = overheads,
                                  best_calculate = True,
                                  average_calculate = False,
                                  worst_calculate = False,
                                  number_of_transits = 2
                                  )
    

"""
_summary_
"""
#%% Importing libraries
import os

import astropy
import astropy.units as u
import specutils as sp
import numpy as np
import matplotlib.pyplot as plt
import expecto
import logging
from copy import deepcopy

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
os.chdir(cwd)

import calculations as calc
import instruments as inst
import plots as plot
import calculate_transit_windows as ctw

logger = logging.getLogger(__name__)
logger = default_logger_format(logger) 

#%%
def simulate_transmission_spectra(target_names: list,
                                  species_list: list,
                                  instrument: str = 'ESPRESSO',
                                  exptime: u.Quantity | None = None,
                                  overheads: u.Quantity = 68*u.s,
                                  best_calculate: bool = True,
                                  average_calculate: bool = True,
                                  worst_calculate: bool = True
                                  ) -> dict:
    
    
    # DOCUMENTME
    # TODO
    logger.info('Starting simulation of transmission dataset')
    
    planet_results = {}
    
    for target in target_names:
        logger.info(f'    Simulation of {target}')
        
        best_results, average_results, worst_results, model = _simulate_transmission_single_planet(
            target= target,
            species_list= species_list,
            instrument= instrument,
            exptime= exptime,
            overheads= overheads,
            best_calculate= best_calculate,
            average_calculate= average_calculate,
            worst_calculate= worst_calculate,
        )
        
        planet_results[target] = {
            'best': best_results,
            'average': average_results,
            'worst': worst_results,
            'model': model
            }
        
    return planet_results

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
    worst_calculate: bool = True
    ) -> [dict, dict, dict, sp.Spectrum1D]:
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

    Returns
    -------
    best_results : dict
        Dictionary holding the best weather condition results for given planet
    average_results : dict
        Dictionary holding the average weather condition results for given planet
    worst_results, model[dict, dict, dict, sp.Spectrum1D]
        Dictionary holding the worst weather condition results for given planet
    model : sp.Spectrum1D
        A model that has been injected into the data. 
    """
    Target_parameters = _load_system_parameters(target_name= target)
    
    best_dataset, average_dataset, worst_dataset, model = _create_dataset(
        Target_parameters= Target_parameters,
        instrument= instrument,
        exptime= exptime,
        overheads= overheads,
        species= species_list)
    
    best_results = {}
    average_results = {}
    worst_results = {}
    
    if best_calculate:
        best_CCF_result, best_data_PRF, best_transmission_spectrum = _run_trasmission_spectroscopy_pipeline(
            dataset= best_dataset,
            Target_parameters= Target_parameters,
            species= species_list
            )
        best_results[target] = [best_CCF_result, best_data_PRF, best_transmission_spectrum]
    if average_calculate:
        average_CCF_result, average_data_PRF, average_transmission_spectrum = _run_trasmission_spectroscopy_pipeline(
            dataset= average_dataset,
            Target_parameters= Target_parameters,
            species= species_list
            )
        average_results[target] = [average_CCF_result, average_data_PRF, average_transmission_spectrum]
    
    if worst_calculate:
        worst_CCF_result, worst_data_PRF, worst_transmission_spectrum = _run_trasmission_spectroscopy_pipeline(
            dataset= worst_dataset,
            Target_parameters= Target_parameters,
            species= species_list
            )
        worst_results[target] = [worst_CCF_result, worst_data_PRF, worst_transmission_spectrum]
    
    return best_results, average_results, worst_results, model

def _create_dataset(Target_parameters: para.SystemParametersComposite,
                    instrument = 'ESPRESSO',
                    exptime: u.Quantity | None = None,
                    overheads: u.Quantity = 68*u.s,
                    species: list = [],
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
        worst_scenario= worst_scenario)
    
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
    worst_scenario: sp.Spectrum1D) -> [sp.SpectrumList, sp.SpectrumList, sp.SpectrumList]:
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
        
        best_dataset.append(_generate_similar_spectrum(best_scenario))
        average_dataset.append(_generate_similar_spectrum(average_scenario))
        worst_dataset.append(_generate_similar_spectrum(worst_scenario))
        
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

def _generate_similar_spectrum(spectrum: sp.Spectrum1D) -> sp.Spectrum1D:
    """
    Generates a similar spectrum given a model spectrum. This is done by adding randomly sampling normal distribution with standard deviation = uncertainty of the spectrum. This sample is then added to the flux array. Uncertainty of the new spectrum is then the square-root of the new flux.

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Model simulated spectrum.

    Returns
    -------
    new_spectrum : sp.Spectrum1D
        New realization of the model spectrum assuming a given noise.
    """
    new_flux = (spectrum.flux.value + np.random.normal(0, spectrum.uncertainty.array, len(spectrum.flux)))*u.dimensionless_unscaled
    
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
                flux = (model.flux.value - (1-specimen_model.flux.value)) *u.dimensionless_unscaled,
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
    dataset = calc._add_mask(dataset)
    master_out = sm.calculate_master_list(dataset,
                                          key= 'Transit_full',
                                          value= False,
                                          )
    
    data_out_corrected = sm.spec_list_master_correct(dataset,
                                                 master_out,
                                                 pkl_name = 'star_corrected_noRM.pkl'
                                                 )
    data_PRF = sm.shift_list(data_out_corrected,
                                 shift_BERV=0,
                                 shift_v_sys = 0,
                                 shift_v_star = -1,
                                 shift_v_planet = 1,
                                 force_multiprocessing= False,
                                 pkl_name = 'data_PRF.pkl'
                                 )
    transmission_spectrum = sm.calculate_master_list(data_PRF,
                                                 key = 'Transit_full',
                                                 value =True,
                                                 sn_type= None,
                                                 pkl_name = 'transmission_spectrum.pkl'
                                                 )
    import rats.spectra_cross_correlate as ccf
    
    result_CCF = {}
    for specimen in species:
        test_model = _model_signal(
            Target_parameters,
            data_PRF[0].spectral_axis,
            species= [specimen]
        )
        
        CCF_model = ccf.cross_correlate_list(
            data_PRF,
            test_model,
            sn_type= None
            )
        result_CCF[specimen] = CCF_model
    return result_CCF, data_PRF, transmission_spectrum[0]

# import rats.plots.plot_spectra_new as ps

def _test():
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.colors as colors

    test = np.linspace(0,180,1)

    t=30
    cut_data_PRF = sm.extract_region_in_list(average_data_PRF, sp.SpectralRegion(5886*u.AA,5900*u.AA))
    fig, ax = plt.subplots(1, figsize= (6*2,4*2))
    spectrum_list = cut_data_PRF
    # spectrum_list = sm.get_sublist(spectrum_list,'Night_num',night+1)

    cmap= sns.diverging_palette(t, t+180, s=75, l=50, sep=1, center='light', as_cmap=True)
    cmap = 'PiYG'



    phase = [item.meta['Phase'] for item in spectrum_list]
    # flux = np.asarray([(item.flux.value - np.nanmedian(item.flux).value) for item in spectrum_list])
    flux = np.asarray([(item.flux.value) for item in spectrum_list])
    vmax = 1.005
    vmin = 0.995

    # bounds = [0.95,0.97,1.03,1.05]
    # norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')
    norm = colors.CenteredNorm(vcenter=1, halfrange=0.01)
    # cax = ax.pcolormesh(spectrum_list[0].spectral_axis.value, phase, flux, cmap= cmap, norm=colors.LogNorm(vmin=0.95, vmax=1.05))
    cax = ax.pcolormesh(spectrum_list[0].spectral_axis.value, phase, flux, cmap= cmap, norm=norm)
    fig.colorbar(cax)

    y_lowest, y_highest = ax.get_ylim()
    ax.vlines(5889.950,y_lowest, Target_parameters.Ephemeris.contact_T1.data,color='black',ls='--')
    ax.vlines(5889.950,Target_parameters.Ephemeris.contact_T4.data,y_highest,color='black',ls='--')
    ax.vlines(5895.924,y_lowest, Target_parameters.Ephemeris.contact_T1.data,color='black',ls='--')
    ax.vlines(5895.924,Target_parameters.Ephemeris.contact_T4.data,y_highest,color='black',ls='--')
    Target_parameters.plot_contact_points(ax)
    ax.set_title(f'Na doublet simulation for {Target_parameters.Planet.name}, Planetary RF, night {spectrum_list[0].meta["Night_num"]}', fontsize= 20)
    ax.set_xlabel('Velocity [km/s]', fontsize= 20)
    ax.set_ylabel('$F_{in}/ F_{out}$ [1]', fontsize= 20)
    
    fig.show()



def plot_TS(transmission_spectrum: sp.Spectrum1D,
            model: sp.Spectrum1D | None = None):
    test= (astropy.modeling.models.Gaussian1D(amplitude= -0.00449,
                                            mean=5889.950*u.AA,
                                            stddev=0.1*u.AA
                                            ) + 
        astropy.modeling.models.Gaussian1D(amplitude= -0.00385,
                                            mean=5895.924*u.AA,
                                            stddev=0.1*u.AA) + 
        astropy.modeling.models.Const1D(amplitude=1)
        )
    import astropy
    import rats.plots.plot_spectra_new as ps
    fig, ax = plt.subplots(1, squeeze = False, figsize= (6*2,4*2))
    transmission_spectrum = sm.extract_region_in_spectrum(average_results['WASP-78 b'][2], sp.SpectralRegion(5886*u.AA,5900*u.AA))
    model = sp.Spectrum1D(
        spectral_axis= model.spectral_axis,
        flux= model.flux.value * u.dimensionless_unscaled,
        uncertainty= astropy.nddata.StdDevUncertainty(model.uncertainty.array),
        mask = np.isnan(model.flux.value),
        meta= {}
    )
    
    model = sm.extract_region_in_spectrum(model, sp.SpectralRegion(5886*u.AA,5900*u.AA))
    ps.PlotSingleSpectrum(
        spectrum = transmission_spectrum,
        binning_factor= 30,
        fig = fig,
        axs = ax,
        color = 'black',
        mode= 'whitemode_presentation'
    )
    
    ax[0,0].plot(model.spectral_axis, model.flux, color='darkred', zorder=3, label='petitRADtrans model')
    # ax[0,0].plot(model.spectral_axis, test(model.spectral_axis)*(np.mean(transmission_spectrum.flux)), color='darkgreen', zorder=3, label='Tabernero+2020')
    ax[0,0].axvline(5889.950,color='darkred',ls='--', label='_nolegend_')
    ax[0,0].axvline(5895.924,color='darkred',ls='--', label='_nolegend_')
    ax[0,0].set_xlim(5886,5900)
    # ax[0,0].set_ylim(0.93, 1.025)
    ax[0,0].set_xlabel('Wavelength [$\AA$]')
    ax[0,0].set_ylabel('Excess absorption [1]')
    # ax[0,0].legend(['petitRADtrans model', 'Tabernero+2020','Simulation'])
    return

def plot_CCF():
    # TODO
    species = 'Na_allard_new'
    species = 'Fe'
    
    test= average_CCF_result[species]
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, figsize=(6*2,4*2))
    x= test[0].spectral_axis.value
    y = np.asarray([item.meta['Phase'] for item in test])
    z = np.asarray([item.flux.value for item in test])
    
    z = z[:] - np.median(z, axis=1, keepdims=True)
    vmax = 1E-6
    vmin = -vmax
    cax = ax.pcolormesh(x,y,z,
                        vmin=vmin,
                        vmax=vmax,
                        cmap='PiYG')
    fig.colorbar(cax)
    
    y_lowest, y_highest = ax.get_ylim()
    ax.vlines(0,y_lowest, Target_parameters.Ephemeris.contact_T1.data,color='black',ls='--')
    ax.vlines(0,Target_parameters.Ephemeris.contact_T4.data,y_highest,color='black',ls='--')
    Target_parameters.plot_contact_points(ax)
    
    ax.set_xlabel('Velocity [km/s]', fontsize=20)
    ax.set_ylabel('Phase [1]', fontsize=20)
    ax.set_title(f'Simulation of species: {species} for planet: {Target_parameters.Planet.name}', fontsize=20)
    
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

if __name__ == '__main__':
    target_names = [
        # 'WASP-193 b',
        # 'TOI-615 b',
        # 'TOI-622 b',
        # 'TOI-2641 b',
        # 'TOI-132 b',
        'WASP-76 b',
        # 'WASP-78 b'
        ]
    species_list= [
        'Na_allard_new',
        'Fe'
        ]
    instrument = 'ESPRESSO'
    exptime= 900*u.s
    overheads = 68*u.s
    simulate_transmission_spectra(target_names= target_names,
                                  species_list= species_list,
                                  instrument = instrument,
                                  exptime = exptime,
                                  overheads = overheads,
                                  best_calculate = False,
                                  average_calculate = True,
                                  worst_calculate = False
                                  )
    
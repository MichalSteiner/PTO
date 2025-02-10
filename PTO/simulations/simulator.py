from dataclasses import dataclass
import astropy.units as u
import specutils as sp
import specutils.manipulation as spm
import pandas as pd
from expecto import get_spectrum
import numpy as np
import matplotlib.pyplot as plt
from . import spectra_manipulation as sm
import astropy
import astropy.constants as con

import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
if not logger.handlers:
    logger = logger_default(logger)

try:
    import petitRADTRANS as prt
except:
    logger.warning('petitRADTRANS instalation not found. Please install petitRADTRANS to use the PTO.simulations module.')
    pass

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



@dataclass
class Simulator():
    
    def define_target(self,
                      information: pd.Series,
                      ):
        
        self.information = information
        
        self.transit_length = information['Planet.TransitDuration'] *u.hour
        self.baseline_length = information['Planet.Baseline']
        
        self._simulate_stellar_spectrum()

        
    def define_instrument_mode(self,
                               instrument_mode):
        ind = np.where(np.logical_and(
            instrument_mode.exposure_time_calculator.start_wavelength < self.default_stellar_spectrum.wavelength,
            self.default_stellar_spectrum.wavelength < instrument_mode.exposure_time_calculator.end_wavelength
            ))
        
        self.default_stellar_spectrum = sp.Spectrum1D(
            spectral_axis= self.default_stellar_spectrum.spectral_axis[ind],
            flux= self.default_stellar_spectrum.flux[ind],
            uncertainty= self.default_stellar_spectrum.uncertainty[ind] if self.default_stellar_spectrum.uncertainty is not None else None,
        )
        
        instrument_mode.exposure_time_calculator.open_all_scenarios(self.information['Star.EffectiveTemperature'])
        
        self.best_SNR, self.mean_SNR, self.worst_SNR = instrument_mode.exposure_time_calculator.best_scenario, instrument_mode.exposure_time_calculator.mean_scenario, instrument_mode.exposure_time_calculator.worst_scenario
        
        self.normalized_stellar_spectrum = sm.normalize_spectrum(self.default_stellar_spectrum)
        
        # Calculate noise level for each point
        noise_level_array = 1/self.best_SNR(self.normalized_stellar_spectrum.spectral_axis.to(u.nm))

        # Generate Gaussian noise for each point
        noise = np.array([np.random.normal(0, sigma) 
                        for sigma in noise_level_array])
        
        self.observed_stellar_spectrum = sp.Spectrum1D(
            spectral_axis= self.normalized_stellar_spectrum.spectral_axis,
            flux= (self.normalized_stellar_spectrum.flux.value + noise) * self.normalized_stellar_spectrum.flux.unit,
            uncertainty= astropy.nddata.StdDevUncertainty(noise_level_array)
        )
        
        self._simulate_CCF_output()

    
    def _simulate_stellar_spectrum(self):
        
        self.default_stellar_spectrum = get_spectrum(T_eff=self.information['Star.EffectiveTemperature'],
                                                     log_g=self.information['Star.Logg'],
                                                     cache=False,
                                                     vacuum=False
                                                     )
    
    def _simulate_CCF_output(self):
        
        from specutils.analysis import correlation
        
        corr, lag = correlation.template_correlate(self.observed_stellar_spectrum, 1-self.normalized_stellar_spectrum)
        
        ind = np.where(np.logical_and(lag < 200*u.km/u.s, lag > -200*u.km/u.s))
        self.default_CCF = sp.Spectrum1D(
            spectral_axis= lag[ind],
            flux= corr[ind],
            uncertainty = astropy.nddata.StdDevUncertainty(np.sqrt(corr[ind]))
        )
    
    def define_planet_petitRADTRANS(self):
        from petitRADTRANS.radtrans import Radtrans
        
        cloud_species = [
            'Al2O3(s)_crystalline__DHS',
            'CaTiO3(s)_crystalline__DHS',
            'Fe(s)_crystalline__DHS',
            'MgSiO3(s)_crystalline__DHS',
            'Cr(s)__Mie',
            'MnS(s)__DHS',
            'Na2S(s)__DHS',
        ]
        
        atmosphere = Radtrans(
            pressures=np.logspace(-10,2,130),
            line_species=[
                'Na',
                'Fe',
                'K'
            ],
            rayleigh_species=['H2', 'He'],
            gas_continuum_contributors=['H2-H2', 'H2-He'],
            wavelength_boundaries=[self.observed_stellar_spectrum.spectral_axis[0].to(u.um).value,
                                   self.observed_stellar_spectrum.spectral_axis[-1].to(u.um).value],
            line_opacity_mode='lbl',
            cloud_species = cloud_species,
        )
        
        # The model includes a single species that is included in both models, otherwise pRT will throw an error
        atmosphere_default_species = Radtrans(
            pressures=np.logspace(-10,2,130),
            line_species=[
                # 'Na',
                # 'Fe',
                'K'
            ],
            rayleigh_species=['H2', 'He'],
            gas_continuum_contributors=['H2-H2', 'H2-He'],
            wavelength_boundaries=[self.observed_stellar_spectrum.spectral_axis[0].to(u.um).value,
                                   self.observed_stellar_spectrum.spectral_axis[-1].to(u.um).value],
            line_opacity_mode='lbl',
            cloud_species = cloud_species,
        )
        
        
        from petitRADTRANS import physical_constants as cst
        from petitRADTRANS.physics import temperature_profile_function_guillot_global
        
        planet_radius = self.information['Planet.RadiusJupiter'] * cst.r_jup_mean
        reference_gravity = ((self.information['Planet.MassJupiter'] *u.M_jup * con.G / (self.information['Planet.RadiusJupiter']**2 *u.R_jup**2)).decompose()).value*100
        reference_pressure = 0.01

        pressures = atmosphere.pressures*1e-6 # cgs to bar
        infrared_mean_opacity = 0.01
        gamma = 0.4
        intrinsic_temperature = 200
        equilibrium_temperature = self.information['Planet.EquilibriumTemperature']

        temperatures = temperature_profile_function_guillot_global(
            pressures=pressures,
            infrared_mean_opacity=infrared_mean_opacity,
            gamma=gamma,
            gravities=reference_gravity,
            intrinsic_temperature=intrinsic_temperature,
            equilibrium_temperature=equilibrium_temperature
        )

        mass_fractions = {
            'H2': 0.74 * np.ones_like(temperatures),
            'He': 0.24 * np.ones_like(temperatures),
            'Fe': 1e-7 * np.ones_like(temperatures),
            'Na': 1e-9 * np.ones_like(temperatures),
            'K' : 1e-9 * np.ones_like(temperatures),
        }
        
        logger.info('Mass fractions for species:')
        for species in mass_fractions:
            logger.info(f"    {species}: {mass_fractions[species][0]}")

        
        cloud_particles_mean_radii = {}
        for species in cloud_species:
            mass_fractions[species] = 0.0000005 * np.ones_like(temperatures)
            cloud_particles_mean_radii[species] = 5e-5 * np.ones_like(temperatures)

        
        cloud_particle_radius_distribution_std = 1.05  # a value of 1.0 would be a delta function, so we assume a very narrow distribtion here
        
        mean_molar_masses = 2.33 * np.ones_like(temperatures)  #  2.33 is a typical value for H2-He dominated atmospheres

        wavelengths, transit_radii, _ = atmosphere.calculate_transit_radii(
            temperatures=temperatures,
            mass_fractions=mass_fractions,
            mean_molar_masses=mean_molar_masses,
            reference_gravity=reference_gravity,
            planet_radius=planet_radius,
            reference_pressure=reference_pressure,
            cloud_particles_mean_radii=cloud_particles_mean_radii,
            cloud_particle_radius_distribution_std = cloud_particle_radius_distribution_std
        )
        wavelengths_no_species, transit_radii_no_species, _ = atmosphere_default_species.calculate_transit_radii(
            temperatures=temperatures,
            mass_fractions=mass_fractions,
            mean_molar_masses=mean_molar_masses,
            reference_gravity=reference_gravity,
            planet_radius=planet_radius,
            reference_pressure=reference_pressure,
            cloud_particles_mean_radii=cloud_particles_mean_radii,
            cloud_particle_radius_distribution_std = cloud_particle_radius_distribution_std
        )
        
        self.planet_spectrum = sp.Spectrum1D(
            spectral_axis = wavelengths*u.um,
            flux = (transit_radii_no_species/transit_radii)*u.dimensionless_unscaled
        )
        raise NotImplementedError('This method is not yet implemented. Fix the planet parameters.')

    
    def _simulate_RM_CLV_effect(self):
        ...
    
    def _get_dataset_spectra(self):
        ...
    
    def _get_dataset_ccf(self):
        ...
    
    def analyze_transmission_spectrum(self):
        ...
    
    def analyze_RM(self):
        ...
    
    

if __name__ == '__main__':
    
    from petitRADTRANS.config import petitradtrans_config_parser
    petitradtrans_config_parser.set_input_data_path(r'/media/chamaeleontis/Observatory_main/prt_test/input_data')  # input_data directory

    
    logger.warning('Debugging mode: Simulator module')
    
    from ..database.NASA_exoplanet_archive import NASA_Exoplanet_Archive_CompositeDefault
    
    Target = NASA_Exoplanet_Archive_CompositeDefault()
    logger.print('Hello there!')
    Target.load_API_table(force_load=False)
    Target.table = Target.table[Target.table['Planet.Name'] == 'WASP-76 b']
    from ..transits.windows import define_baseline
    Target.table = define_baseline(Target.table)
    
    from ..telescopes.instruments import ESPRESSO
    
    Simulation = Simulator()
    Simulation.define_target(Target.table.iloc[0])
    Simulation.define_instrument_mode(ESPRESSO.modes[0])
    Simulation.define_planet_petitRADTRANS()
    
    
    logger.warning('Debugging finished')

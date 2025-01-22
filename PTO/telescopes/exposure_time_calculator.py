from dataclasses import dataclass
import os
import numpy as np
import specutils as sp
import astropy.units as u
from scipy.interpolate import interp1d

def _find_stellar_type_ESPRESSO(T_eff: float) -> str:
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

def read_ESO_ETC_output(filename:str):
    # Initialize lists to store the data
    wavelength = []
    SNR = []

    # Read the file
    with open(filename, 'r') as file:
        for line in file:
            # Ignore commented lines
            if line.startswith('#'):
                continue
            # Split the line into wavelength and SNR
            parts = line.split()
            if len(parts) == 2:
                wavelength.append(float(parts[0]))
                SNR.append(float(parts[1]))

    # Convert lists to numpy arrays
    wavelength = np.array(wavelength)
    SNR = np.array(SNR)

    # Create an interpolation function
    interpolation_function = interp1d(wavelength, SNR, kind='linear', fill_value="extrapolate")
    return interpolation_function

class ETC():
    def __post_init__(self):
        self.spectral_axis = np.arange(self.start_wavelength.value,
                                       self.end_wavelength.value,
                                       (self.start_wavelength.value + self.end_wavelength.value)/ 2  / (self.resolution * self.spectral_resolution_sampling)) *u.nm
    
    pass


@dataclass
class ETC_ESPRESSO_1UT_HR(ETC):
    start_wavelength = 378.2 * u.nm
    end_wavelength = 788.7 * u.nm
    resolution = 70000
    spectral_resolution_sampling = 3.5
    
    def _get_stellar_type(self, stellar_temperature: float) -> str:
        return _find_stellar_type_ESPRESSO(stellar_temperature)
    
    
    
    def _load_file(self,
                   stellar_temperature:float,
                   seeing:float):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        file_path = f'{current_dir}/ETC/ETC_ESPRESSO/{self._get_stellar_type(stellar_temperature)}/{seeing}arcsec_900s.etc'

        interpolation_function = read_ESO_ETC_output(filename= file_path)

        return interpolation_function
    
    def open_all_scenarios(self, stellar_temperature: float):
        self.best_scenario = self._load_file(stellar_temperature, 0.5)
        self.mean_scenario = self._load_file(stellar_temperature, 0.8)
        self.worst_scenario = self._load_file(stellar_temperature, 1.3)


@dataclass
class ETC_ESPRESSO_4UT(ETC):
    start_wavelength = 378.2 * u.nm
    end_wavelength = 788.7 * u.nm
    resolution = 70000
    spectral_resolution_sampling = 3.5
    
    def _get_stellar_type(self,
                        stellar_temperature:float):
        return _find_stellar_type_ESPRESSO(stellar_temperature)

    def _load_file(self,
                   stellar_temperature:float,
                   seeing:float):
        
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = f'{current_dir}/ETC/ETC_ESPRESSO_4UT/{self._get_stellar_type(stellar_temperature)}/{seeing}arcsec_900s.etc'
        interpolation_function = read_ESO_ETC_output(filename= file_path)

        return interpolation_function
    
    def open_all_scenarios(self, stellar_temperature: float):
        self.best_scenario = self._load_file(stellar_temperature, 0.5)
        self.mean_scenario = self._load_file(stellar_temperature, 0.8)
        self.worst_scenario = self._load_file(stellar_temperature, 1.3)
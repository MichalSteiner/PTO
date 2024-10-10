#%% Importing libraries
from dataclasses import dataclass
import numpy as np
import astropy.units as u
from enum import Enum


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
        self.wavelength_range = np.linspace(3771.2787376605725, 7905.599289073489, 443262) * u.AA # Number of S1D pixels 
        return

@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class ESPRESSO_4UT_wavelength:
    """Class holding the wavelength range of ESPRESSO 4UT mode."""
    
    def __post_init__(self):
        """
        Estimated wavelength range as taken from real observation.

        Returns
        -------
        None.

        """
        self.wavelength_range = np.linspace(3771.2787376605725, 7905.599289073489,221570) *u.AA # Number of S1D pixels
        return
#%%
Instruments = {
    'ESPRESSO': ESPRESSO_wavelength().wavelength_range,
    'ESPRESSO_4UT': ESPRESSO_4UT_wavelength().wavelength_range,
    'HARPS': '',
    'NIRPS': '',
    }

    
class Instruments_size(Enum):
    """Enum class holding wavelength ranges for instruments."""
    
    ESPRESSO = 8.2 * u.m
    ESPRESSO_4UT = 4*8.2*u.m
    HARPS = 3.6*u.m
    NIRPS = 3.6*u.m
    ANDES = 39 * u.m

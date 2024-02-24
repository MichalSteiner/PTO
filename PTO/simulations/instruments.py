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
        self.wavelength_range = np.array([
        np.linspace(3771.2787376605725, 7905.599289073489,170*9111),
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

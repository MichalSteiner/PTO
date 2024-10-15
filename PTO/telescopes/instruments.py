import astropy.units as u
from dataclasses import dataclass, field
import numpy as np
import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 

@dataclass
class Spectrograph():
    
    name: str
    resolution_mode: str
    modes: list = field(default_factory=list) 
    
    def _add_observing_mode(self,
                            mode_name:str,
                            resolution:int,
                            simultaneous = None,
                            simulator:bool = False):
        self.modes.append(
            Mode(
                mode_name=mode_name,
                resolution=resolution,
                simultaneous=simultaneous,
                simulator= simulator
            )
        )

@dataclass
class Mode():
    mode_name: str
    resolution: int
    simultaneous: bool = None
    simulator: bool = False
    pass


#%% List of instruments available
#%%% CARMENES
CARMENES = Spectrograph(
    'CARMENES',
    resolution_mode= 'high'
    )   

CARMENES._add_observing_mode(
    'VIS',
    94600,
    True
)

CARMENES._add_observing_mode(
    'NIR',
    80400,
    True
)
#%%% CORALIE
CORALIE = Spectrograph(
    'CORALIE',
    resolution_mode= 'high',
)

CORALIE._add_observing_mode(
    'HR',
    50000,
)

#%%% ESPRESSO, 1 UT
ESPRESSO = Spectrograph(
    'ESPRESSO',
    resolution_mode='high',
)

ESPRESSO._add_observing_mode(
    'HR',
    134000,
)

ESPRESSO._add_observing_mode(
    'UHR',
    190000,
)

ESPRESSO_4UT = Spectrograph(
    'ESPRESSO_4UT',
    resolution_mode='high'
)
ESPRESSO_4UT._add_observing_mode(
    'MR',
    70000,
)
#%%% EXPRES
EXPRES = Spectrograph(
    'EXPRES',
    resolution_mode='high'
)

EXPRES._add_observing_mode(
    'HR',
    137000,
)
#%%% GIANO
GIANO = Spectrograph(
    'GIANO',
    resolution_mode='high'
)

GIANO._add_observing_mode(
    'HR',
    50000,
)

GIANO._add_observing_mode(
    'LR',
    25000,
)
#%%% HARPS
HARPS = Spectrograph(
    'HARPS',
    resolution_mode='high'
)

HARPS._add_observing_mode(
    'HR',
    115000,
)

#%%% HARPS-N
HARPS_N = Spectrograph(
    'HARPS-N',
    resolution_mode='high'
)

HARPS_N._add_observing_mode(
    'HR',
    115000,
    'with_GIANO',
)
#%%% MAROON-X
MAROON_X = Spectrograph(
    'MAROON-X',
    resolution_mode='high'
)

MAROON_X._add_observing_mode(
    'HR',
    80000,
)

MAROON_X._add_observing_mode(
    'UHR',
    100000,
)
#%%% NIRPS
NIRPS = Spectrograph(
    'NIRPS',
    resolution_mode= 'high'
)
NIRPS._add_observing_mode(
    'HR',
    100000,
    'with HARPS',
)

NIRPS._add_observing_mode(
    'MR',
    75000,
    'with HARPS',
)
#%%% SOPHIE
SOPHIE = Spectrograph(
    'SOPHIE',
    resolution_mode= 'high'
)

SOPHIE._add_observing_mode(
    'HR',
    75000,
)
SOPHIE._add_observing_mode(
    'HE',
    40000,
)

#%%% SPIRou
SPIROU = Spectrograph(
    'SPIRou',
    resolution_mode= 'high'
)
SPIROU._add_observing_mode(
    'HR',
    70000,
)
#%%% UVES
UVES = Spectrograph(
    'UVES',
    resolution_mode= 'high'
)

UVES._add_observing_mode(
    'HR',
    np.mean([40000, 110000]),
)

def print_all_spectrographs(resolution: str = 'all') -> None:
    """
    Prints all spectrographs based on 

    Parameters
    ----------
    resolution : str, optional
        Resolution mode filter, by default 'all'. If 'high' or 'low' filter is used, only high-resolution or low-resolution spectrographs are provided.
    """
    if resolution != 'all':
        raise NotImplementedError
    
    for var_name, var_value in globals().items():
        if isinstance(var_value, Spectrograph):
            if resolution == 'all' or var_value.resolution_mode == resolution:
                logger.print(f"{var_name}: {var_value.name}")

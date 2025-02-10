import astropy.units as u
from dataclasses import dataclass, field
import numpy as np
import logging
from ..utils.utilities import logger_default
from ..simulations.simulator import Simulator
from . import exposure_time_calculator as etc

logger = logging.getLogger(__name__)
if not logger.handlers:
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
                            exposure_time_calculator: None | etc.ETC = None,
                            simulator: bool = False):
        self.modes.append(
            Mode(
                mode_name=mode_name,
                resolution=resolution,
                simultaneous=simultaneous,
                simulator= simulator,
                exposure_time_calculator= exposure_time_calculator
            )
        )

@dataclass
class Mode(Simulator):
    mode_name: str
    resolution: int
    simultaneous: bool = None
    simulator: bool = False,
    exposure_time_calculator: None | etc.ETC = None
    
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
    'with NIR CARMENES mode'
)

CARMENES._add_observing_mode(
    'NIR',
    80400,
    'with VIS CARMENES mode'
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
    exposure_time_calculator= etc.ETC_ESPRESSO_1UT_HR()
)

ESPRESSO._add_observing_mode(
    'UHR',
    190000
)

ESPRESSO_4UT = Spectrograph(
    'ESPRESSO_4UT',
    resolution_mode='high'
)
ESPRESSO_4UT._add_observing_mode(
    'MR',
    70000,
    exposure_time_calculator= etc.ETC_ESPRESSO_4UT()
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
    'with NIRPS',
)

#%%% HARPS-N
HARPS_N = Spectrograph(
    'HARPS-N',
    resolution_mode='high'
)

HARPS_N._add_observing_mode(
    'HR',
    115000,
    'with GIANO',
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
                logger.print(f'='*25)
                logger.print(f"{var_name}: {var_value.name}")
                logger.print('Modes:')
                for mode in var_value.modes:
                    logger.print(f"    {mode.mode_name}: R= {mode.resolution} | Simulator: {mode.simulator} | Simultaneous to: {mode.simultaneous} | ETC: {True if mode.exposure_time_calculator is not None else False}")

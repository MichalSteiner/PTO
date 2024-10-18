import pandas as pd
from dataclasses import dataclass
import dill as pickle
import os
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from ..utils.utilities import logger_default
from .plot import PlotUtilitiesComposite
from .calculation import CalculationUtilities

logger = logging.getLogger(__name__)
logger = logger_default(logger) 

class Utilities():
    """
    This class holds basic utilities for both Catalogs Classes.
    """
    
    def save(self) -> None:
        """
        Saves the class into a predefined filename.
        """
        if not(os.path.exists('./saved_files')):
            logger.info('Directory not found, creating new directory in:')
            logger.info(f'    {os.getcwd()}/saved_files')
            os.makedirs('./saved_files', mode = 0o777, exist_ok = False) 
        
        with open(f'./saved_files/{self.filename}', 'wb') as output_file:
            pickle.dump(self.__dict__, output_file)
            logger.info(f'File saved succesfully in:')
            logger.info(f'    {os.getcwd()}/saved_files/{self.filename}')
        
        return None 
        
    def load(self) -> None:
        """
        Loads the class from predefined filename.
        
        If this functions fails to load, the standard initialization routine will run instead.
        """
        try:
            with open(f'./saved_files/{self.filename}', 'rb') as input_file:
                self.__dict__  =  pickle.load(input_file)
        except:
            logger.info(f'Failed to load file:')
            logger.info(f'    {os.getcwd()}/saved_files/{self.filename}')
        return None
    
    def _print_keys(self,
                   keytype:str) -> None:
        """
        Prints keys available in the table attribute.

        Parameters
        ----------
        keytype : str
            Keytype to look for. Keywords should start with this keytype to be printed. E.g., to select planet related keys, use keytype='Planet'
        """
        for key in self.table.keys():
            if key.startswith(keytype):
                logger.print(key)
        return None
        
    
    def print_position_keys(self) -> None:
        """
        Prints all position keys.
        """
        logger.print('='*25)
        logger.print('Position keys:')
        logger.print('='*25)
        self._print_keys(keytype='Position.')
    
    def print_system_keys(self) -> None:
        """
        Prints all system keys.
        """
        logger.print('='*25)
        logger.print('System keys:')
        logger.print('='*25)
        self._print_keys(keytype='System.')
    
    def print_star_keys(self) -> None:
        """
        Prints all star keys.
        """
        logger.print('='*25)
        logger.print('Stellar keys:')
        logger.print('='*25)
        self._print_keys(keytype='Star.')
        
    def print_planet_keys(self) -> None:
        """
        Prints all planet keys.
        """
        logger.print('='*25)
        logger.print('Planet keys:')
        logger.print('='*25)
        self._print_keys(keytype='Planet.')
    
    def print_discovery_keys(self) -> None:
        """
        Prints all discovery keys.
        """
        logger.print('='*25)
        logger.print('Discovery keys:')
        logger.print('='*25)
        self._print_keys(keytype='Discovery.')

    def print_magnitude_keys(self) -> None:
        """
        Prints all magnitude keys.
        """
        logger.print('='*25)
        logger.print('Magnitude keys:')
        logger.print('='*25)
        self._print_keys(keytype='Magnitude.')
        
    def print_flag_keys(self) -> None:
        """
        Prints all flag keys.
        """
        logger.print('='*25)
        logger.print('Flag keys:')
        logger.print('='*25)
        self._print_keys(keytype='Flag.')
    
    def print_all_keys(self) -> None:
        """
        Prints all keys.
        """
        self.print_position_keys()
        self.print_system_keys()
        self.print_star_keys()
        self.print_planet_keys()
        self.print_discovery_keys()
        self.print_magnitude_keys()
        self.print_flag_keys()
    



@dataclass
class CatalogComposite(Utilities,
                       PlotUtilitiesComposite,
                       CalculationUtilities):
    table: pd.DataFrame | None = None
    filename: str = 'CatalogComposite.pkl'
    drop_mode: str = 'drop'
    
    def _get_all(self) -> None:
        """
        This functions calculates the missing values that can be calculated.
        """
        if self.drop_mode == 'drop':
            logger.warning('    Droping all values without errorbars. To instead replace the errorbars with 0 change Catalogs "drop_mode" key to "replace"')
        if self.drop_mode == 'replace':
            logger.info('Replacing NaN errorbars with 0, if the value is defined')
        self._handle_keys_without_errors(mode=self.drop_mode)
        
        logger.info('Converting Earth and Jupiter Radius units in the catalog')
        self._add_Earth_and_Jupiter_units()
        
        logger.info('Checking for inclination and impact parameters values')
        self._calculate_impact_parameter()
        
        logger.info('Calculation of T_14 and related values')
        self._calculate_transit_length()

class CatalogFull(Utilities):
    table: pd.DataFrame | None = None
    filename: str = 'CatalogFull.pkl'

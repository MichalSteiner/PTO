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
import csv



logger = logging.getLogger(__name__)
logger = logger_default(logger) 

class _Utilities():
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
        logger.info(f'Filename: {os.getcwd()}/saved_files/{self.filename}')
        
        try:
            with open(f'{os.getcwd()}/saved_files/{self.filename}', 'rb') as input_file:
                self.__dict__  =  pickle.load(input_file)
        except:
            logger.info(f'Failed to load file:')
            logger.info(f'    {os.getcwd()}/saved_files/{self.filename}')
        return
    
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

    def create_exoplanet_csv(self,
                             filename='exoplanet_data.csv'):
        """
        Creates an empty CSV file with exoplanet-related headers if it does not already exist.

        Args:
            filename (str): Name of the CSV file to create (default: 'exoplanet_data.csv')
        """
        if not os.path.exists(filename):
            headers = [
                'Planet.Name',
                'Position.RightAscension',
                'Position.Declination',
                'Magnitude.V',
                'Planet.TransitDuration',
                'Planet.TransitMidpoint',
                'Planet.TransitMidpoint.Error',
                'Planet.Period',
                'Planet.Period.Error',
            ]

            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=';')
                writer.writerow(headers)
                logger.info(f'File {filename} created successfully.')
        else:
            logger.info(f'File {filename} already exists.')
            
    def load_exoplanet_csv(self,
                           filename:str):
        """
        Loads exoplanet data from a CSV file into the table attribute.

        Args:
            filename (str): Name of the CSV file to load
        """
        try:
            self.table = pd.read_csv(filename, delimiter=';')
            logger.info(f'File {filename} loaded successfully.')
            self._expand_error()
        except Exception as e:
            logger.error(f'Failed to load file {filename}: {e}')
            
    
    def _expand_error(self):
        """
        Expands error columns in the table by creating corresponding lower and upper error columns.
        
        This method iterates through the keys of the table and for each key that ends with '.Error',
        it creates two new keys: one for the lower error and one for the upper error. Both new keys
        are initially set to the same value as the original error key.
        
        Example:
            If the table has a key 'value.Error', this method will create 'value.Error.Lower' and
            'value.Error.Upper' with the same values as 'value.Error'.
        """
        for key in self.table.keys():
            if key.endswith('.Error'):
                self.table[f"{key}.Lower"] = self.table[key]
                self.table[f"{key}.Upper"] = self.table[key]
        
    
    
    def print_ephemeris(self):
        logger.print('='*25)
        logger.print('Used ephemeris values:')
        logger.print('='*25)
        for ind, row in self.table.iterrows():
            logger.print(f"Planet with indice {ind}: RA:{row['Position.RightAsce']}")
            



@dataclass
class CatalogComposite(_Utilities,
                       PlotUtilitiesComposite,
                       CalculationUtilities):
    """
    CatalogComposite class is a class that combines various utility classes and provides methods to handle and process exoplanet composite catalog data.
    
    Attributes:
        table (pd.DataFrame | None): The main data table containing the catalog data.
        filename (str): The filename to save or load the catalog data. Default is 'CatalogComposite.pkl'.
        drop_mode (str): The mode to handle missing error bars. Can be 'drop' to drop rows without error bars or 'replace' to replace NaN error bars with 0. Default is 'drop'.
        legacy_table (pd.DataFrame | None): An optional legacy data table for backward compatibility.
        
    Methods:
        _get_all() -> None:
            Calculates the missing values that can be calculated in the catalog.
            - Logs a warning if drop_mode is 'drop' and drops all values without error bars.
            - Logs an info message if drop_mode is 'replace' and replaces NaN error bars with 0.
            - Fills in Earth and Jupiter Radius units in the catalog if empty.
            - Fills in semimajor axis over stellar radius ratio in the catalog if empty.
            - Checks for inclination and impact parameters values.
            - Calculates T_14 and related values.
    """
    
    table: pd.DataFrame | None = None
    filename: str = 'CatalogComposite.pkl'
    drop_mode: str = 'drop'
    legacy_table: pd.DataFrame | None = None
    
    def _get_all(self) -> None:
        """
        This functions calculates the missing values that can be calculated.
        """
        if self.drop_mode == 'drop':
            logger.warning('    Droping all values without errorbars. To instead replace the errorbars with 0 change Catalogs "drop_mode" key to "replace"')
        if self.drop_mode == 'replace':
            logger.info('Replacing NaN errorbars with 0, if the value is defined')
        self._handle_keys_without_errors(mode=self.drop_mode)
        
        logger.info('Filling in Earth and Jupiter Radius units in the catalog if empty')
        self._add_Earth_and_Jupiter_units()
        

        logger.info('Filling in semimajor axis over stellar radius ratio in the catalog if empty')
        self._calculate_R_s_a()
        
        logger.info('Checking for inclination and impact parameters values')
        self._calculate_impact_parameter()
        
        logger.info('Calculation of T_14 and related values')
        self._calculate_transit_length()

class CatalogFull(_Utilities, CalculationUtilities):
    """
    CatalogComposite class is a class that combines various utility classes and provides methods to handle and process exoplanet full catalog data.
    
    Attributes:
        table (pd.DataFrame | None): The main data table containing the catalog data.
        filename (str): The filename to save or load the catalog data. Default is 'CatalogComposite.pkl'.
        drop_mode (str): The mode to handle missing error bars. Can be 'drop' to drop rows without error bars or 'replace' to replace NaN error bars with 0. Default is 'drop'.
        legacy_table (pd.DataFrame | None): An optional legacy data table for backward compatibility.
        
    Methods:
        _get_all() -> None:
            Calculates the missing values that can be calculated in the catalog.
            - Logs a warning if drop_mode is 'drop' and drops all values without error bars.
            - Logs an info message if drop_mode is 'replace' and replaces NaN error bars with 0.
            - Fills in Earth and Jupiter Radius units in the catalog if empty.
            - Fills in semimajor axis over stellar radius ratio in the catalog if empty.
            - Checks for inclination and impact parameters values.
            - Calculates T_14 and related values.
    """
    
    
    table: pd.DataFrame | None = None
    filename: str = 'CatalogFull.pkl'
    drop_mode: str = 'drop'
    legacy_table: pd.DataFrame | None = None

    def _get_all(self) -> None:
        """
        This function calculates the missing values that can be calculated.
        """
        if self.drop_mode == 'drop':
            logger.warning('    Dropping all values without errorbars. To instead replace the errorbars with 0 change Catalogs "drop_mode" key to "replace"')
        if self.drop_mode == 'replace':
            logger.info('Replacing NaN errorbars with 0, if the value is defined')
        
        self._handle_keys_without_errors(mode=self.drop_mode)
        
        logger.info('Filling in Earth and Jupiter Radius units in the catalog if empty')
        self._add_Earth_and_Jupiter_units()
        
        logger.info('Filling in semimajor axis over stellar radius ratio in the catalog if empty')
        self._calculate_R_s_a()
        
        logger.info('Checking for inclination and impact parameters values')
        self._calculate_impact_parameter()
        
        logger.info('Calculation of T_14 and related values')
        self._calculate_transit_length()
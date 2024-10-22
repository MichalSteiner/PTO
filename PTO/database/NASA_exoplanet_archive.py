from . import catalog as cat
import datetime
import pyvo as vo
import pandas as pd
import numpy as np
import logging
from ..utils.utilities import logger_default
from .mappers import _NASA_EXOPLANET_ARCHIVE_COMPOSITE_MAPPER

logger = logging.getLogger(__name__)
logger = logger_default(logger) 

class NASA_Exoplanet_Archive_CompositeDefault(cat.CatalogComposite):
    def load_API_table(self, force_load=False) -> None:
        """
        Loads the Table using the API system, in particular the TAP protocol. 
        
        This is rerun every week, but the output is saved and by default loaded instead of rerunning the TAP protocol

        Parameters
        ----------
        force_load : bool, optional
            Flag to trigger reloading of the TAP protocol, by default False. If False, the self.filename is going to loaded, and only if a week or more passed the TAP protocol is rerun. If True or if last run happened week or more ago, the TAP protocal is relaunched. 
        """
        try:
            if force_load:
                logger.info('Forced reload:')
                raise
            logger.info('Trying to load NASA Exoplanet Archive Composite table')
            self.load()
            if (datetime.datetime.now() - self.time) < datetime.datetime(days=7):
                logger.info('Too old data, reloading:')
                raise
        except:
            logger.info('Accessing NASA Exoplanet Archive')
            service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP/") 
            logger.info('Fetching table')
            self.table = pd.DataFrame(service.search("SELECT * FROM pscomppars"))
            logger.info('Table fetched successfully')
            self.time = datetime.datetime.now()
            self._rename_columns()
            self._drop_columns()
            self._absolute_errors()
            self._get_all()
            self.legacy_table = self.table
            self.save()
            
    def _rename_columns(self) -> None:
        """
        Renames the columns in the pandas dataframe
        """
        self.table = self.table.rename(columns= _NASA_EXOPLANET_ARCHIVE_COMPOSITE_MAPPER)
    
    def _drop_columns(self) -> None:
        """
        Drops all columns that are irrelevant.
        """
        _TODROP = [key for key in self.table.keys() if
                    not(key.startswith('Planet.')) and 
                    not(key.startswith('Star.')) and 
                    not(key.startswith('Magnitude.')) and 
                    not(key.startswith('Position.')) and 
                    not(key.startswith('System.')) and 
                    not(key.startswith('Flag.')) and 
                    not(key.startswith('Discovery.'))
                   ]
        self.table = self.table.drop(_TODROP, axis=1)
    
    def _absolute_errors(self) -> None:
        """
        Reverses the sign of the lower error to absolute value. This ensures further functionality.
        """
        keys = [key for key in self.table.keys() if key.endswith('.Lower')]
        
        for key in keys:
            self.table[key] = np.abs(self.table[key])

class NASA_Exoplanet_Archive_CompositeMostPrecise():
    ...
    
class NASA_Exoplanet_Archive_FullTable():
    ...
    
if __name__ == '__main__':
    import os
    os.chdir('/media/chamaeleontis/Observatory_main/Code/observations_transits/PTO/')
    test = NASA_Exoplanet_Archive_CompositeDefault()
    logger.print('Hello there!')
    test.load_API_table(force_load=True)
    test.print_all_keys()
    # fig, ax = test.plot_diagram(
    #     x_key = 'Planet.Period',
    #     y_key = 'Planet.RadiusJupiter',
    # )

    logger.print(f"Length before further filtering of the table: {test.table.shape[0]}")
    test.table = test.table[test.table['Magnitude.V'] < 10]
    test.table = test.table[test.table['Planet.RadiusEarth'] > 3]
    test.table = test.table[test.table['Planet.RadiusEarth'] < 8]
    test.table = test.table[test.table['Planet.Period'] < 30]
    logger.print(f"Length after further filtering of the table: {test.table.shape[0]}")

    fig, ax = test.highlight_sample(
        x_key= 'Planet.Period',
        y_key= 'Planet.RadiusEarth',
        marker= 'd'
    )
    ax.set_xlim(0.1,50)
    logger.print('General Kenobi!!!!')
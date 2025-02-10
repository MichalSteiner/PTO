from . import catalog as cat
import datetime
import pyvo as vo
import pandas as pd
import numpy as np
import logging
from ..utils.utilities import logger_default
from .mappers import _NASA_EXOPLANET_ARCHIVE_COMPOSITE_MAPPER, _NASA_EXOPLANET_ARCHIVE_FULL_MAPPER

logger = logging.getLogger(__name__)
if not logger.handlers:
    logger = logger_default(logger)

class NASA_Exoplanet_Archive_CompositeDefault(cat.CatalogComposite):
    """
    NASA_Exoplanet_Archive_CompositeDefault is a class that extends the CatalogComposite class to handle the loading and processing of the NASA Exoplanet Archive Composite table using the TAP protocol.

    Methods
    -------
    :meth:`.NASA_Exoplanet_Archive_CompositeDefault.load_API_table`
        Loads the Table using the API system, in particular the TAP protocol. This method is rerun every week, but the output is saved and by default loaded instead of rerunning the TAP protocol.
    :meth:`.NASA_Exoplanet_Archive_CompositeDefault._rename_columns`
        Renames the columns in the pandas dataframe according to a predefined mapping.
    :meth:`.NASA_Exoplanet_Archive_CompositeDefault._drop_columns`
        Drops all columns that are irrelevant from the pandas dataframe.
    :meth:`.NASA_Exoplanet_Archive_CompositeDefault._absolute_errors`
        Reverses the sign of the lower error to an absolute value to ensure further functionality.
    """
    
    def load_API_table(self, force_load=False) -> None:
        """
        Loads the Table using the API system, in particular the TAP protocol.

        This is rerun every week, but the output is saved and by default loaded instead of rerunning the TAP protocol.

        Parameters
        ----------
        force_load : bool
            Flag to trigger reloading of the TAP protocol.
                If False, the self.filename is going to be loaded, and only if a week 
                or more passed the TAP protocol is rerun. If True or if last run 
                happened week or more ago, the TAP protocal is relaunched.
                (default: False)
        """
        
        try:
            if force_load:
                logger.info('Forced reload:')
                raise
            logger.info('Trying to load NASA Exoplanet Archive Composite table')
            self.load()
            if (datetime.datetime.now() - self.time) > datetime.timedelta(days=7):
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


class NASA_Exoplanet_Archive_FullTable(cat.CatalogFull):
    """
    NASA_Exoplanet_Archive_FullTable is a class that handles the loading and processing of the NASA Exoplanet Archive Composite table using the TAP protocol. It extends the CatalogFull class from the Catalog module.
    
    Methods
    -------
    load_API_table(force_load=False) -> None
        Loads the Table using the API system, in particular the TAP protocol. This method is rerun every week, but the output is saved and by default loaded instead of rerunning the TAP protocol.
    _rename_columns() -> None
        Renames the columns in the pandas dataframe according to a predefined mapping.
    _drop_columns() -> None
        Drops all columns that are irrelevant from the pandas dataframe.
    get_most_precise_value(group, column)
        Returns the most precise value for a given column within a group, based on the smallest error range.
    aggregate_most_precise_values() -> None
        Aggregates the most precise values for each planet by considering the error ranges for relevant columns.
    """
    
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
            logger.info('Success!')
            if (datetime.datetime.now() - self.time) > datetime.timedelta(days=7):
                logger.info('Too old data, reloading:')
                raise
        except:
            logger.info('Accessing NASA Exoplanet Archive')
            service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP/") 
            logger.info('Fetching table')
            self.table = pd.DataFrame(service.search("SELECT * FROM ps"))
            logger.info('Table fetched successfully')
            self.time = datetime.datetime.now()
            self._rename_columns()
            self._drop_columns()
            # self._absolute_errors()
            self._get_all()
            
            self.aggregate_most_precise_values()
            
            self.legacy_table = self.table
            self.save()
            
    def _rename_columns(self) -> None:
        """
        Renames the columns in the pandas dataframe
        """
        self.table = self.table.rename(columns= _NASA_EXOPLANET_ARCHIVE_FULL_MAPPER)

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


class NASA_Exoplanet_Archive_CompositeMostPrecise(NASA_Exoplanet_Archive_FullTable):
    
    
    def get_most_precise_value(self, group, column):
        error_upper = column + '.Error.Upper'
        error_lower = column + '.Error.Lower'
        if error_upper in group.columns and error_lower in group.columns and (group[error_upper].notna() != False).all() and (group[error_lower].notna() != False).all():
            group['error_range'] = group[error_upper] - group[error_lower]
            most_precise_row = group.loc[group['error_range'].idxmin()]
            return most_precise_row[column]
        else:
            return group[column].iloc[0]

    def aggregate_most_precise_values(self):
        # List of columns to consider for precision
        columns_to_consider = [col for col in self.table.columns if not col.endswith('.Error.Upper') and not col.endswith('.Error.Lower') and f"{col}.Error.Upper" in self.table.columns and f"{col}.Error.Lower" in self.table.columns]
        
        # Group by 'Planet.Name'
        grouped = self.table.groupby('Planet.Name')
        
        # Create a new DataFrame to store the aggregated results
        aggregated_df = []
        
        for name, group in grouped:
            aggregated_row = {'Planet.Name': name}
            for column in columns_to_consider:
                if column != 'Planet.Name':
                    aggregated_row[column] = self.get_most_precise_value(group, column)
            aggregated_df.append(aggregated_row)
        
        self.table = pd.DataFrame(aggregated_df)
    
    
    ...
    


if __name__ == '__main__':
    logger.warning('Debugging Database module')
    logger.warning('='*25)
    import os
    os.chdir('/media/chamaeleontis/Observatory_main/Code/observations_transits/PTO/')
    
    test_full = NASA_Exoplanet_Archive_FullTable()
    test_full.load_API_table(force_load=True)
    logger.print('Hello there!')
    
    test = NASA_Exoplanet_Archive_CompositeDefault()
    logger.print('Hello there!')
    test.load_API_table(force_load=False)
    # test.print_all_keys()
    # fig, ax = test.plot_population_diagram(
    #     x_key = 'Planet.Period',
    #     y_key = 'Planet.RadiusJupiter',
    # )

    logger.print(f"Length before further filtering of the table: {test.table.shape[0]}")
    test.table = test.table[test.table['Planet.Name'] == 'WASP-31 b']
    logger.print(f"Length after further filtering of the table: {test.table.shape[0]}")
    
    # fig, ax = test.highlight_sample(
    #     x_key= 'Planet.Period',
    #     y_key= 'Planet.RadiusJupiter',
    # )
    # ax.set_xlim(0.1,50)
    # logger.print('General Kenobi!!!!')
    
    
    
    from . import catalog as cat
    from ..telescopes import telescopes as tel
    from ..transits.windows import Windows

    # ATREIDES = cat.CatalogComposite()
    # ATREIDES.create_exoplanet_csv('ATREIDES_custom.csv')
    # ATREIDES.load_exoplanet_csv('ATREIDES_custom.csv')



    for P in [108, 109, 110, 111]:
        Transits = Windows(
            table = test.table,
            directory= f'/media/chamaeleontis/Observatory_main/ESO_scheduling/PTO_developement/',
            observing_period= f'ESO.{P}'
        )
        
        Transits.directory = f'/media/chamaeleontis/Observatory_main/ESO_scheduling/PTO_developement/'
        Transits.generate_observability(
            location= tel.VLT,
            partial= False,
            velocity_offset= 0,
            velocity_range = 5
        )


#%%
# def transform_time_format(input_string):
#     # Split the input string
#     target, date_str, quality, time_start, time_end, _ = input_string.split(';')

#     # Parse the date
#     year = int(date_str[0:4])
#     month = int(date_str[4:6])
#     day = int(date_str[6:8])

#     # Parse the time
#     hour_start = int(time_start.split(':')[0])
#     minute_start = int(time_start.split(':')[1])

#     # If hour < 12, it's morning of the next day
#     if hour_start < 12:
#         observation_date = datetime(year, month, day) + timedelta(days=1)
#     else:
#         observation_date = datetime(year, month, day)

#     hour_end = int(time_end.split(':')[0])
#     minute_end = int(time_end.split(':')[1])

#     if hour_end < 12:
#         observation_date_end = datetime(year, month, day) + timedelta(days=1)
#     else:
#         observation_date_end = datetime(year, month, day)


#     # Create the start time
#     start_time = datetime(observation_date.year, 
#                          observation_date.month, 
#                          observation_date.day, 
#                          hour_start, 
#                          minute_start)
    
#     end_time = datetime(observation_date_end.year, 
#                         observation_date_end.month, 
#                         observation_date_end.day, 
#                         hour_end, 
#                         minute_end)
    
#     # Format to ISO time string
#     start_iso = start_time.strftime("%Y-%m-%dT%H:%M")
#     end_iso = end_time.strftime("%Y-%m-%dT%H:%M")

#     # Create the output string
#     output = f"between({start_iso},{end_iso},{1},\"{target} P117\")"

#     return output

# Example usage
# input_string = "CoRoT-22 b;20260522;10;03:37;09:54;1"
# result = transform_time_format(input_string)
# print(result)

# # %%
# list_of_strings = [
#     'TOI-3071 b;20260112;1;03:45;07:35;1',
#     'TOI-3071 b;20260131;1;03:50;07:40;1',
#     'TOI-3071 b;20260214;1;02:19;06:09;1',
#     'TOI-3071 b;20260219;1;03:56;07:46;1',
#     'TOI-3071 b;20260228;1;00:47;04:37;1',
#     'TOI-3071 b;20260305;1;02:25;06:15;1',
#     'TOI-3071 b;20260310;1;04:02;07:52;1',
#     'TOI-3071 b;20260319;1;00:53;04:43;1',
#     'TOI-3071 b;20260324;1;02:31;06:21;1',
#     'TOI-3071 b;20260329;1;04:08;07:58;1',
# ]


# list_of_strings = [
#     'TOI-3071 b;20250407;1;03:54;07:44;1',
#     'TOI-3071 b;20250416;1;00:45;04:35;1',
#     'TOI-3071 b;20250421;1;02:23;06:13;1',
#     'TOI-3071 b;20250505;1;00:51;04:41;1',
#     'TOI-3071 b;20250402;1;02:17;06:07;1',
#     'TOI-3071 b;20250524;1;00:57;04:47;1',
#     'TOI-3071 b;20250607;1;23:25;03:15;1',
# ]

# list_of_strings = [
#     'TOI-3071 b;20260407;1;00:59;04:49;1',
#     'TOI-3071 b;20260412;1;02:36;06:26;1',
#     'TOI-3071 b;20260426;1;01:05;04:55;1',
#     'TOI-3071 b;20260501;1;02:42;06:32;1',
#     'TOI-3071 b;20260510;1;23:33;03:23;1',
#     'TOI-3071 b;20260515;1;01:11;05:01;1',
#     'TOI-3071 b;20260529;1;23:39;03:29;1',
#     'TOI-3071 b;20260617;1;23:45;03:35;1',
# ]
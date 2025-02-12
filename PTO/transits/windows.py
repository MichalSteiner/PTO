from dataclasses import dataclass, field
import pandas as pd
import astropy.time as astime
import astropy.units as u
import logging
import csv
import numpy as np
import os
from ..utils.utilities import logger_default
from datetime import datetime, timedelta
from .observability import Event
from ..telescopes.telescopes import Telescope, VLT

logger = logging.getLogger(__name__)
if not logger.handlers:
    logger = logger_default(logger)


def define_baseline(table):
    array_length = len(table)
    table['Planet.Baseline'] = np.nanmax([np.nanmin([[3]*array_length, .75*table['Planet.TransitDuration']], axis = 0),[2]*array_length],                                            axis=0)
    return table


@dataclass
class Windows:
    """
    A class to handle the generation and management of transit windows for exoplanet observations.
    
    Attributes
    ----------
    table : pd.DataFrame
        DataFrame containing the exoplanet data.
    center_phase : float, optional
        The center phase of the transit (default is 0.0).
    observing_period : astime.Time | None | str, optional
        The observing period for the transits (default is None).
    baseline : u.Quantity | None, optional
        The baseline for the transit windows (default is None).
    large_program : bool, optional
        Flag indicating if the program is a large program (default is False).
    directory : str, optional
        Directory to save output files (default is '').
    Airmass_limit : float, optional
        The airmass limit for observations (default is None).
    
    Methods
    -------
    generate_windows():
        Generate transit windows for the exoplanets in the table.
    print_windows():
        Print the transit windows for each exoplanet.
    __post_init__():
        Initialize the observing period and generate transit windows.
    get_proposal_period():
        Determine the observing period based on the provided input.
    _get_period_for_next_year():
        Set the observing period to the next year.
    _get_proposal_period_from_string():
        Parse the observing period from a string.
    get_dates_for_ESO_semester(P_number: int) -> astime.Time:
        Get the start and end dates for an ESO semester.
    to_eso_format(time_start, time_end, target_name, quality=1) -> str:
    generate_observability(location: Telescope, partial: float = 1, velocity_offset: None | float = None, velocity_range: float = 5, save_figures: bool = True):
    define_baseline():
    add_header_to_csv(csv_writer: csv.writer):
        Add a header row to a CSV file using the provided CSV writer.
    add_event_to_csv(csv_writer: csv.writer, event: Event):
        Add an event's details to a CSV file using the provided CSV writer.
    """
    
    table: pd.DataFrame
    center_phase: float = 0.0
    observing_period: astime.Time | None | str = None
    baseline: u.Quantity | None = None
    large_program: bool = False
    directory: str = ''
    Airmass_limit: float = None
    
    
    def generate_windows(self):
        """
        Generate transit windows for planets in the table.
        
        This method calculates the midpoints of transit windows for each planet in the table
        within the specified observing period. It also calculates the associated errors for
        these midpoints.
        Raises:
            ValueError: If the table is empty and no calculations can be performed.
        Returns:
            None
        """
        
        if len(self.table) == 0:
            raise ValueError('Table is empty. Cannot calculate windows')
        else:
            logger.info(f'About to calculate event midpoints for {len(self.table)} planets')
        
        Ti_values = []
        Ti_sigma_values = []
        
        T1, T2 = self.observing_period[0], self.observing_period[1]
        
        for _, row in self.table.iterrows():
            T_c = row['Planet.TransitMidpoint'] * u.day
            T_c_error = np.max([row['Planet.TransitMidpoint.Error.Lower'], row['Planet.TransitMidpoint.Error.Upper']]) * u.day
            P = row['Planet.Period'] * u.day
            P_error = np.max([row['Planet.Period.Error.Lower'], row['Planet.Period.Error.Upper']]) * u.day
            
            n_min = np.ceil(((T1 - T_c).jd)*u.day / P)  # Round up to next integer
            n_max = np.floor(((T2 - T_c).jd)*u.day / P)  # Round down to previous integer

            n_values = np.arange(n_min, n_max + 1, dtype=int)

            Ti = T_c + n_values * P
            sigma_values = np.sqrt(T_c_error**2 + n_values**2 * P_error**2)
            # Shouldn't change, but to be sure we are within the observing period
            valid_Ti = Ti[(Ti > T1.jd *u.day) & (Ti < T2.jd*u.day)]
            valid_sigma = sigma_values[(Ti > T1.jd*u.day) & (Ti < T2.jd*u.day)].to(u.min)

            valid_Ti = [astime.Time(t, format='jd') for t in valid_Ti]
            
            Ti_values.append(valid_Ti)
            Ti_sigma_values.append(valid_sigma)

        self.table['Planet.TransitWindowCenter'] = Ti_values
        self.table['Planet.TransitWindowCenter.Error'] = Ti_sigma_values

        return
    
    def print_windows(self):
        """
        Prints the transit windows for each planet in the table.
        
        Iterates over each row in the table and prints the transit window centers
        along with their uncertainties for each planet. The output is formatted
        with a separator line and includes the planet name, transit window center,
        uncertainty, and the formatted date and time in UT.
        
        Returns:
            None
        """
        
        for _, row in self.table.iterrows():
            logger.print('='*25)
            logger.print(f"Transit windows centers for {row['Planet.Name']}")
            for window, uncertainty in zip(row['Planet.TransitWindowCenter'], row['Planet.TransitWindowCenter.Error']):
                logger.print(f"    {window} ± {uncertainty:.2f} | {window.strftime('%Y-%m-%d %H:%M')} UT")
    
    def __post_init__(self):
        """
        Post-initialization method to set up the observing period, generate windows, 
        and define the baseline. This method is automatically called after the 
        object's initialization.
        
        It performs the following steps:
        1. Retrieves the proposal period and logs the observing period.
        2. Generates observation windows.
        3. Defines the baseline for observations.
        """
        
        self.get_proposal_period()
        logger.info('='*25)
        logger.info('Set observing period:')
        logger.info(f'    {self.observing_period[0]}')
        logger.info(f'    {self.observing_period[1]}')
        logger.info('='*25)
        
        self.generate_windows()
        self.define_baseline()
        
        
    
    
    def get_proposal_period(self):
        """
        Determines the proposal period based on the observing period attribute.
        
        This method checks the type of the `observing_period` attribute and 
        determines the proposal period accordingly. If `observing_period` is 
        None, it calls `_get_period_for_next_year()`. If `observing_period` is 
        a string, it calls `_get_proposal_period_from_string()`. If 
        `observing_period` is of type `astime.Time`, it does nothing. For any 
        other type, it raises a `ValueError`.
        Raises:
            ValueError: If `observing_period` is not None, str, or astime.Time.
        """
        
        if self.observing_period is None:
            self._get_period_for_next_year()
        elif type(self.observing_period) == str:
            self._get_proposal_period_from_string()
        elif type(self.observing_period) == astime.Time:
            return
        else:
            raise ValueError('Invalid type for observing period.') 

    def _get_period_for_next_year(self):
        """
        Calculate the observing period for the next year.
        This method sets the `observing_period` attribute to a time range
        starting from today at 12:00 PM to the same date and time next year.
        The dates are formatted as 'YYYY-MM-DD 12:00:00.000' and converted
        to an `astime.Time` object with UTC scale.
        
        Returns:
            None
        """
        
        today = datetime.now().replace(hour=12, minute=0, second=0, microsecond=0)
        next_year = today.replace(year=today.year + 1)

        date1 = today.strftime('%Y-%m-%d 12:00:00.000')
        date2 = next_year.strftime('%Y-%m-%d 12:00:00.000')
        
        self.observing_period = astime.Time([date1, date2], scale='utc')

    def _get_proposal_period_from_string(self):
        """
        Extracts the proposal period from the observing period string and sets the observing period
        to the corresponding dates for the ESO semester.
        If the observing period string starts with 'ESO.', it extracts the numeric part of the string,
        converts it to an integer, and uses it to get the dates for the ESO semester.
        
        Raises:
            ValueError: If the observing period string does not start with 'ESO.'.
        """
        
        if self.observing_period.startswith('ESO.'):
            self.observing_period = self.get_dates_for_ESO_semester(
                P_number= int(''.join(c for c in self.observing_period if c.isdigit())),
            )
        else: 
            raise ValueError('Invalid string format.')

    
    def get_dates_for_ESO_semester(self,
                                   P_number: int,
                                   ) -> astime.Time:
        """
        Provide starting and ending date of the ESO semester.

        Parameters
        ----------
        P_number : int
            Which semester are we looking for.
            
        Returns
        -------
        time : astime.Time
            Starting and ending date of the semester.

        """
        if not(self.large_program):
            if P_number%2 ==0:
                time = astime.Time([f'{1967+ int(P_number/2)}-10-01 12:00:00.000',  
                                    f'{1968+int(P_number/2)}-04-01 12:00:00.000'], scale = 'utc')
            else:
                time = astime.Time([f'{1968 + int((P_number-1)/2)}-04-01 12:00:00.000',
                                    f'{1968+ int((P_number-1)/2)}-10-01 12:00:00.000'], scale = 'utc')
        else:
            if P_number%2 ==0:
                time = astime.Time([f'{1967+ int(P_number/2)}-10-01 12:00:00.000',  
                                    f'{1969+ int(P_number/2)}-10-01 12:00:00.000'], scale = 'utc')
            else:
                time = astime.Time([f'{1968 + int((P_number-1)/2)}-04-01 12:00:00.000',
                                    f'{1970 + int((P_number-1)/2)}-04-01 12:00:00.000'], scale = 'utc')
        return time
    
    def to_eso_format(time_start, time_end, target_name, quality=1):
        """
        Convert two astropy Time objects to ESO's between() format.

        Parameters
        ----------
        time_start : astropy.time.Time
            Start time
        time_end : astropy.time.Time
            End time
        target_name : str
            Name of the target (e.g., "K2-406 b partial")
        quality : int, optional
            Quality parameter (default=1)

        Returns
        -------
        str
            ESO formatted time string

        Examples
        --------
        >>> from astropy.time import Time
        >>> t1 = Time('2025-06-24T23:24:00')
        >>> t2 = Time('2025-06-25T06:16:00')
        >>> print(to_eso_format(t1, t2, "K2-406 b partial"))
        between(2025-06-24T23:24,2025-06-25T06:16,1,"K2-406 b partial")
        """
        # Convert to ISO format and remove seconds if present
        t1_str = time_start.iso.split('.')[0].replace(':00', '')
        t2_str = time_end.iso.split('.')[0].replace(':00', '')

        return f'between({t1_str},{t2_str},{quality},"{target_name}")'
    
    
    def generate_observability(self,
                               location: Telescope,
                               partial: float = 1,
                               velocity_offset: None | float = None,
                               velocity_range: float = 5, #km/s,
                               save_figures: bool = True,
                               ):
        """
        Generate observability windows for transiting exoplanets.
        
        Parameters:
        -----------
        location : Telescope
            The telescope location for which to generate observability windows.
        partial : float, optional
            The fraction of the transit duration to consider for partial observability (default is 1).
        velocity_offset : None or float, optional
            The velocity offset to apply (default is None).
        velocity_range : float, optional
            The range of velocities to consider in km/s (default is 5).
        save_figures : bool, optional
            Whether to save figures generated during the process (default is True).
        
        Returns:
        --------
        None
        """
        
        if not(self.Airmass_limit):
            if location.name == 'Very Large Telescope (VLT)':
                logger.info('Setting airmass limit to 2.2.')
                self.Airmass_limit = 2.2
            else:
                logger.warning('No airmass limit set. Setting to 2')
                self.Airmass_limit = 2
        
        os.makedirs(os.path.join(self.directory, location.name), exist_ok=True)
        complete_summary = os.path.join(self.directory, location.name, 'complete_summary.csv')
        csvfile_complete = open(complete_summary, 'w', newline='')
        csv_writer_complete = csv.writer(csvfile_complete)
        self.add_header_to_csv(csv_writer_complete)

        
        for _, row in self.table.iterrows():
            logger.print('='*25)
            logger.print(f"Working on {row['Planet.Name']}")
            logger.print(f"    with Tc: {row['Planet.TransitMidpoint']}, P: {row['Planet.Period']} days and T14: {row['Planet.TransitDuration']} hours")
            logger.print('='*25)
            valid_events = []
            
            for window, window_uncertainty in zip(row['Planet.TransitWindowCenter'], row['Planet.TransitWindowCenter.Error']):
                new_event = Event(
                        Location=location,
                        TransitDuration= row['Planet.TransitDuration'],
                        baseline=row['Planet.Baseline'],
                        Night= window,
                        row= row,
                        directory= self.directory,
                        Uncertainty= window_uncertainty,
                        partial= partial,
                        Airmass_limit= self.Airmass_limit,
                        velocity_offset= velocity_offset,
                        velocity_range= velocity_range,
                        save_figures= save_figures,
                    )
                
                
                
                if new_event.quality > 0:
                    valid_events.append(new_event)
            
            if valid_events:
                per_planet_summary = os.path.join(self.directory, location.name, row['Planet.Name'].replace(' ',''), 'planet_summary.csv')
                with open(per_planet_summary, 'w', newline='') as csvfile_per_planet:
                    csv_writer_per_planet = csv.writer(csvfile_per_planet)
                    self.add_header_to_csv(csv_writer_per_planet)
                
                    # Write data for all valid events
                    for event in valid_events:
                        self.add_event_to_csv(csv_writer_per_planet, event)
                        self.add_event_to_csv(csv_writer_complete, event)
        csvfile_complete.close()
        self.generate_calendar(complete_summary)

    def define_baseline(self):
        """
        Define the baseline for the transit windows.
        
        This method checks if the baseline attribute is None. If it is, it calls the 
        `define_baseline` function to set the baseline for the table. Otherwise, it 
        sets the 'Planet.Baseline' column in the table to the baseline value converted 
        to hours.
        
        Returns:
            None
        """
        
        
        if self.baseline is None:
            self.table = define_baseline(self.table)
        else:
            self.table['Planet.Baseline'] = self.baseline.to(u.hour)
        
        return

    def add_header_to_csv(self,
                          csv_writer: csv.writer):
        """
        Adds a header row to a CSV file using the provided CSV writer.
        
        Parameters:
        -----------
        csv_writer : csv.writer
            The CSV writer object used to write to the CSV file.
        
        The header row contains the following columns:
        - 'Planet Name'
        - 'Night'
        - 'Observation start'
        - 'Observation end'
        - 'Quality'
        - 'Period'
        - 'Transit midpoint'
        - 'Transit center'
        - 'Transit center error [min]'
        - 'SM mode observable'
        """
        csv_writer.writerow([
            'Planet Name',
            'Night',
            'Observation start',
            'Observation end',
            'Quality',
            'Period',
            'Transit midpoint',
            'Transit center',
            'Transit center error [min]',
            'SM mode observable'
        ])
        return

    def add_event_to_csv(self,
                         csv_writer:csv.writer,
                         event:Event):
        """
        Adds an event's details to a CSV file using the provided CSV writer.

        Parameters:
        -----------
        csv_writer : csv.writer
            The CSV writer object used to write to the CSV file.
        event : Event
            The event object containing details about the transit event.

        The following details are written to the CSV file:
        - Planet Name
        - Night (formatted as YYYYMMDD)
        - Start Time (formatted as HH:MM)
        - End Time (formatted as HH:MM)
        - Quality of the event
        - Planet's orbital period
        - Planet's transit midpoint
        - Night of the event
        - Uncertainty in minutes
        - Whether the transit is observable by Service Mode using ESO rules
        """
        
        
        csv_writer.writerow([
            event.row['Planet.Name'],  # Planet Name
            event.TimeArray.midnight.datetime.strftime('%Y%m%d'),              # Night
            event.TimeArray.time_array[event.ObservationsFull.complete[0]].strftime('%H:%M'),  # Start Time
            event.TimeArray.time_array[event.ObservationsFull.complete[-1]].strftime('%H:%M'),  # End Time
            event.quality,       # Quality
            event.row['Planet.Period'],
            event.row['Planet.TransitMidpoint'],
            event.Night,
            event.Uncertainty.to(u.min).value,
            len(event.ObservationsFull.Ingress) > 30,
        ])
        return

    
    def generate_calendar(self, complete_summary: str):
        
        import pandas as pd
        import calplot
        import numpy as np
        
        from matplotlib.colors import ListedColormap,BoundaryNorm
        
        start_calendar = self.observing_period[0].datetime
        end_calendar = self.observing_period[-1].datetime
        
        df = pd.read_csv(complete_summary)
        df['Night'] = pd.to_datetime(df['Night'], format='%Y%m%d')
        
        def create_event_calendar(events_df, date_column, start_date, end_date):

            events_series = events_df[date_column].value_counts()
            events_series.index = pd.to_datetime(events_series.index)
            date_range = pd.date_range(start=start_date, end=end_date, freq='D').normalize()
            events_series = events_series.reindex(date_range, fill_value=0)
            
            
            # Create custom colormap with grey for out-of-range dates
            colors = ['green', 'green', 'green', 'orange']  # lightgrey for out-of-range
            bounds = [-0.5, 0.9, 1.1, 5.5]  # or similar small range around 1

            custom_cmap = ListedColormap(colors)
            norm = BoundaryNorm(bounds, custom_cmap.N)
            
            # Plot the calendar
            fig, ax = calplot.calplot(
                events_series,
                cmap=custom_cmap,  # Red-Yellow-Green colormap
                fillcolor='black',  # Color for days without data
                daylabels='MTWTFSS',
                dayticks=[0, 1, 2, 3, 4, 5, 6],  # Show some day labels
                figsize=(16, 8),
                suptitle='Transits events',
                colorbar=False
            )
            return fig

        # Usage example:
        # start_date and end_date should be datetime objects
        # df is your DataFrame with events
        fig = create_event_calendar(
            df, 
            'Night',
            start_calendar,
            end_calendar
        )
        fig.savefig(complete_summary.replace('csv', 'png'))
    
    
if __name__ == '__main__':
    import os
    from ..database.NASA_exoplanet_archive import NASA_Exoplanet_Archive_CompositeDefault
    
    os.chdir('/media/chamaeleontis/Observatory_main/Code/observations_transits/PTO/')
    test = NASA_Exoplanet_Archive_CompositeDefault()
    test.load_API_table(force_load=True)
    
    logger.print(f"Length before further filtering of the table: {test.table.shape[0]}")
    test.table = test.table[test.table['Magnitude.V'] < 10]
    test.table = test.table[test.table['Planet.RadiusEarth'] > 3]
    test.table = test.table[test.table['Planet.RadiusEarth'] < 8]
    test.table = test.table[test.table['Planet.Period'] < 30]
    # test.table = test.table[test.table['Planet.Name'].isin(['HD 209458 b', 'HD 189733 b', 'WASP-76 b', ])]
    
    logger.print(f"Length after further filtering of the table: {test.table.shape[0]}")

    
    Transits = Windows(
        table = test.table,
        directory= '/media/chamaeleontis/Observatory_main/ESO_scheduling/PTO_developement/',
        large_program= False
    )
    
    Transits.print_windows()
    
    Transits.generate_observability(
        location= VLT,
        save_figures= False,
    )
    
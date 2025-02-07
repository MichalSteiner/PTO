from dataclasses import dataclass, field
import pandas as pd
import astropy.time as astime
import astropy.units as u
import logging
import numpy as np
from ..utils.utilities import logger_default
from datetime import datetime, timedelta
from .observability import Event
from ..telescopes.telescopes import Telescope, VLT

logger = logging.getLogger(__name__)
logger = logger_default(logger) 

def define_baseline(table):
    array_length = len(table)
    table['Planet.Baseline'] = np.nanmax([np.nanmin([[3]*array_length, .75*table['Planet.TransitDuration']], axis = 0),[2]*array_length],                                            axis=0)
    return table


@dataclass
class Windows:
    table: pd.DataFrame
    center_phase: float = 0.0
    observing_period: astime.Time | None | str = None
    baseline: u.Quantity | None = None
    large_program: bool = False
    directory: str = ''
    Airmass_limit: float = None
    
    windows_in_period:list =  field(default_factory=list)
    
    
    def generate_windows(self):
        
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
        for _, row in self.table.iterrows():
            logger.print('='*25)
            logger.print(f"Transit windows centers for {row['Planet.Name']}")
            for window, uncertainty in zip(row['Planet.TransitWindowCenter'], row['Planet.TransitWindowCenter.Error']):
                logger.print(f"    {window} ± {uncertainty:.2f} | {window.strftime('%Y-%m-%d %H:%M')} UT")
    
    def __post_init__(self):
        self.get_proposal_period()
        logger.info('='*25)
        logger.info('Set observing period:')
        logger.info(f'    {self.observing_period[0]}')
        logger.info(f'    {self.observing_period[1]}')
        logger.info('='*25)
        
        self.generate_windows()
        self.define_baseline()
        
        
    
    
    def get_proposal_period(self):
        if self.observing_period is None:
            self._get_period_for_next_year()
        elif type(self.observing_period) == str:
            self._get_proposal_period_from_string()
        elif type(self.observing_period) == astime.Time:
            return
        else:
            raise ValueError('Invalid type for observing period.') 

    def _get_period_for_next_year(self):
        today = datetime.now().replace(hour=12, minute=0, second=0, microsecond=0)
        next_year = today.replace(year=today.year + 1)

        date1 = today.strftime('%Y-%m-%d 12:00:00.000')
        date2 = next_year.strftime('%Y-%m-%d 12:00:00.000')
        
        self.observing_period = astime.Time([date1, date2], scale='utc')

    def _get_proposal_period_from_string(self):
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
                               velocity_range: float = 5 #km/s
                               ):
        
        if not(self.Airmass_limit):
            if location.name == 'Very Large Telescope (VLT)':
                logger.info('Setting airmass limit to 2.2.')
                self.Airmass_limit = 2.2
            else:
                logger.warning('No airmass limit set. Setting to 2')
                self.Airmass_limit = 2
        
        for _, row in self.table.iterrows():
            logger.print('='*25)
            logger.print(f"Working on {row['Planet.Name']}")
            logger.print(f"    with Tc: {row['Planet.TransitMidpoint']}, P: {row['Planet.Period']} days and T14: {row['Planet.TransitDuration']} hours")
            logger.print('='*25)
            
            
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
                        velocity_range= velocity_range
                    )
                
                self.windows_in_period.append(new_event)
                

    def define_baseline(self):
        
        if self.baseline is None:
            self.table = define_baseline(self.table)
        else:
            self.table['Planet.Baseline'] = self.baseline.to(u.hour)
        
        return



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
    logger.print(f"Length after further filtering of the table: {test.table.shape[0]}")

    
    Transits = Windows(
        table = test.table,
        observing_period = 'ESO.115',
        directory= '/media/chamaeleontis/Observatory_main/ESO_scheduling/PTO_developement/',
        large_program= False
    )
    
    Transits.print_windows()
    
    Transits.generate_observability(
        location= VLT
    )
    
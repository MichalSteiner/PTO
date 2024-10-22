from dataclasses import dataclass, field
import pandas as pd
import astropy.time as astime
import astropy.units as u
import numpy as np
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 


class Windows():
    table: pd.DataFrame
    center_phase: float = 0.0
    observing_period: astime.Time
    baseline: u.Quantity
    
    windows_in_period =  field(default_factory=list)
    
    
    def generate_windows(self):
        # Extract columns from the DataFrame
        P = self.table['Planet.Period'].values * u.day  # Orbital period
        P_error = np.max([self.table['Planet.Period.Error.Lower'], self.table['Planet.Period.Error.Upper'].values])
        
        T_C = astime.Time(self.table['Planet.TransitMidpoint'].values, format='jd') + (self.center_phase * P)
        T_C_error = np.maximum(self.table['Planet.TransitMidpoint.Error.Lower'].values, self.table['Planet.TransitMidpoint.Error.Upper'].values)
        
        date_start = self.observing_period[0]
        date_end = self.observing_period[1]

        # Initialize result lists
        all_transit_windows = [[] for _ in range(len(T_C))]

        # Step 1: Ensure that all T_C are less than the starting date
        n_pre = np.zeros(len(T_C))  # Counter for pre-adjustments
        while np.any(T_C > date_start):
            mask = T_C > date_start
            T_C[mask] -= P[mask]
            n_pre[mask] += 1

        # Step 2: Now calculate the transit windows from date_start to date_end
        n_post = np.zeros(len(T_C))  # Counter for post-adjustments
        while np.any(T_C < date_end):
            mask = T_C < date_end
            T_C[mask] += P[mask]
            n_post[mask] += 1

        # Check for valid transit windows in the desired date range
        for i in range(len(T_C)):
            if date_start < T_C[i] < date_end:
                n_total = n_pre[i] + n_post[i]
                sigma_T0 = np.sqrt(T_C_error[i]**2 + n_total**2 * P_error[i]**2)
                all_transit_windows[i].append([T_C[i], sigma_T0 * u.day])
    
    
    
    pass



if __name__ == '__main__':
    import os
    from ..database.NASA_exoplanet_archive import NASA_Exoplanet_Archive_CompositeDefault
    
    os.chdir('/media/chamaeleontis/Observatory_main/Code/observations_transits/PTO/')
    test = NASA_Exoplanet_Archive_CompositeDefault()
    test.load_API_table(force_load=True)
    test.print_all_keys()
    
    logger.print(f"Length before further filtering of the table: {test.table.shape[0]}")
    test.table = test.table[test.table['Magnitude.V'] < 10]
    test.table = test.table[test.table['Planet.RadiusEarth'] > 3]
    test.table = test.table[test.table['Planet.RadiusEarth'] < 8]
    test.table = test.table[test.table['Planet.Period'] < 30]
    logger.print(f"Length after further filtering of the table: {test.table.shape[0]}")


    Transits = Windows(
        table = test.table
        observing_period = astime.Time()
        baseline: u.Quantity
    )
    
    
    logger.print('General Kenobi!!!!')
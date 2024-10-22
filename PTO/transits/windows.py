from dataclasses import dataclass
import pandas as pd
import astropy.time as astime
import astropy.units as u

class Windows():
    table: pd.DataFrame
    center_phase: float = 0.0
    observing_period: astime.Time
    baseline: u.Quantity
    
    
    
    
    
    
    pass
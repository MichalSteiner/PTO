import specutils as sp
import pandas as pd
import numpy as np
from astropy.modeling import models, fitting

def normalize_spectrum(spectrum: sp.Spectrum1D):
    normalization_function = pd.Series(spectrum.flux).rolling(10000,min_periods=1,center=True).quantile(0.85)
    new_spectrum = spectrum.divide(normalization_function.values * spectrum.flux.unit)
    return new_spectrum
import sys
import pandas as pd
import datetime
import astropy.units as u
import numpy as np
import astropy.time as astime
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass, field
import re
import matplotlib.dates as md 
import os
from ..telescopes import telescopes as tel
from .plot import WindowPlot
import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
if not logger.handlers:
    logger = logger_default(logger)



@dataclass
class _TimeArray:
    """Class holding timings of the transit."""
    
    def _define_times(self,
                      Night: astime.Time,
                      TransitDuration: float) -> None:
        """
        Define basic times and time arrays for the plot.

        Parameters
        ----------
        Night : astime.Time
            The central time of the night.
        TransitDuration : float
            Duration of the transit in hours.

        Returns
        -------
        None
        """
        self.midnight = astime.Time(round(Night.value) - 0.5,
                                    format=Night.format)
        
        self.time_array = np.linspace(-12, 12, 24*60+1) * u.hour + self.midnight
        self.T1, self.T4 = ((Night - TransitDuration * u.hour / 2),
                            (Night + TransitDuration * u.hour / 2))
        return None
    
    def _define_sunset_and_sunrise(self,
                                   sunset: datetime.datetime,
                                   sunrise: datetime.datetime) -> None:
        """
        Define the times for sunset and sunrise.

        Parameters
        ----------
        sunset : datetime.datetime
            Time of sunset.
        sunrise : datetime.datetime
            Time of sunrise.

        Returns
        -------
        None
        """
        self.Sunset = sunset
        self.Sunrise = sunrise
        return None
        

@dataclass
class _AltitudeArray:
    """Class holding altitudes of Sun, Moon, and target."""
    
    def _define_targets(self,
                        TimeArray: _TimeArray,
                        Location: tel.Telescope,
                        Planet_row: pd.Series) -> None:
        """
        Define the altitude-azimuth frame and calculate positions of the target, Moon, and Sun.

        Parameters
        ----------
        TimeArray : _TimeArray
            Object holding the time array for the observations.
        Location : tel.Telescope
            Telescope location.
        Planet_row : pd.Series
            Series containing planetary parameters.

        Returns
        -------
        None
        """
        self.altitude_azimuth_frame = coord.AltAz(
            obstime=TimeArray.time_array,
            location=Location.location
            )
        self.Sun = coord.get_sun(TimeArray.time_array).transform_to(self.altitude_azimuth_frame)
        self.Moon = coord.get_body("moon", TimeArray.time_array, location=Location.location).transform_to(self.altitude_azimuth_frame)
        self.target = coord.SkyCoord(
            ra=coord.Angle(Planet_row['Position.RightAscension'] * u.deg),
            dec=coord.Angle(Planet_row['Position.Declination'] * u.deg),
            frame='icrs'
        ).transform_to(self.altitude_azimuth_frame)
        self.MoonSeparation = np.min(self.Moon.separation(self.target))
        
        ind = np.argwhere((self.Sun.alt.value < 0) * self.Sun.alt.value != 0)
        TimeArray._define_sunset_and_sunrise(TimeArray.time_array[ind[0]].to_value(format='datetime')[0],
                                             TimeArray.time_array[ind[-1]].to_value(format='datetime')[0])
        
        return None

@dataclass
class _Indices:
    """Class holding indices of various arrays."""
    
    def _define_indices(self,
                        TimeArray: _TimeArray,
                        AltitudeArray: _AltitudeArray,
                        Airmass_Limit: float) -> None:
        """
        Define the indices of transit, ingress, egress, and out-of-transit.

        Parameters
        ----------
        TimeArray : _TimeArray
            Object holding the time array for the observations.
        AltitudeArray : _AltitudeArray
            Object holding the altitude-azimuth frame and positions of the target, Moon, and Sun.
        Airmass_Limit : float
            The airmass limit for the observations.

        Returns
        -------
        None
        """
        self.Transit = np.where(
            np.logical_and(
                TimeArray.time_array > TimeArray.T1,
                TimeArray.time_array < TimeArray.T4
                )
            )
        self.Out_of_transit = np.where(
            np.logical_or(
                TimeArray.time_array < TimeArray.T1,
                TimeArray.time_array > TimeArray.T4
                )
            )
        self.Ingress = np.where(
                TimeArray.time_array < TimeArray.T1,
            )
        self.Egress = np.where(
                TimeArray.time_array > TimeArray.T4
            )
        self.Observable = np.where(
            np.logical_and(
                AltitudeArray.target.secz < Airmass_Limit,
                AltitudeArray.target.secz >= 1
                )
            )
        return None


class _Visibility:
    """Class holding the indices when baseline, target, transit, ingress, and egress is visible."""
    
    def _define_visibility(self,
                           NightDefiningIndices: np.ndarray,
                           Indices: _Indices) -> None:
        """
        Define the visibility depending on full/twilight mode and indices of transit.

        Parameters
        ----------
        NightDefiningIndices : np.ndarray
            Indices defining the night period.
        Indices : _Indices
            Object holding the indices of various arrays.

        Returns
        -------
        None
        """
        self.Visibility = np.intersect1d(Indices.Observable, NightDefiningIndices)
        self.Transit = np.intersect1d(self.Visibility,
                                      Indices.Transit)
        self.Baseline = np.intersect1d(self.Visibility,
                                       Indices.Out_of_transit)
        self.Ingress = np.intersect1d(self.Visibility,
                                      Indices.Ingress)
        self.Egress = np.intersect1d(self.Visibility,
                                      Indices.Egress)
        return

class _Observations:
    """Class holding the indices for the observations.
    
    This class defines the observation timings for ingress, egress, and transit phases
    based on the visibility of the target. It ensures that the proposed observation time
    cover the required baseline length for both ingress and egress phases.
    """
    
    def _define_observations(self,
                             Visibility: _Visibility,
                             baseline_length: float) -> None:
        """
        Define the observation periods for ingress, egress, and transit phases.

        Parameters
        ----------
        Visibility : _Visibility
            Visibility of the target's transit.
        baseline_length : float
            Required baseline length in minutes.

        Returns
        -------
        None
        """
        baseline_length = int(round(baseline_length * 60 / 2))
        
        def adjust_visibility(visibility, baseline_length, ingress= True):
            if len(visibility) > baseline_length and ingress:
                return visibility[-baseline_length:], 0
            elif len(visibility) > baseline_length and not(ingress):
                return visibility[:baseline_length], 0
            else:
                return visibility, baseline_length - len(visibility)
            
        self.Ingress, missing_length_ingress = adjust_visibility(Visibility.Ingress, baseline_length)
        self.Egress, missing_length_egress = adjust_visibility(Visibility.Egress, baseline_length, ingress=False)
        
        if missing_length_ingress != 0:
            if len(Visibility.Egress) > (baseline_length + missing_length_ingress):
                self.Egress = Visibility.Egress[:(baseline_length + missing_length_ingress)]
            else:
                self.Egress = Visibility.Egress
        
        if missing_length_egress != 0:
            if len(Visibility.Ingress) > (baseline_length + missing_length_egress):
                self.Ingress = Visibility.Ingress[-(baseline_length + missing_length_egress):]
            else:
                self.Ingress = Visibility.Ingress
        
        self.Transit = Visibility.Transit
        self.complete = np.concatenate([self.Ingress, self.Transit, self.Egress])

        return

#%% Flags transit
@dataclass
class _Flags_transit:
    """Flags for a given transit opportunity in terms of visibility."""
    def check(self,
              Visibility: _Visibility,
              row: pd.Series,
              baseline_length: float,
              partial: bool = False) -> None:
        """
        Check the values of flags.

        Parameters
        ----------
        Visibility : _Visibility
            Visibility of the target's transit. Can be full night or twilight-included version.
        row : pd.Series
            Series containing planetary parameters.
        baseline_length : float
            Required baseline length in minutes.
        partial : bool, optional
            Minimum required coverage for partial visibility. Default is False.

        Returns
        -------
        None
        """
        if partial == False:
            partial = 1.0
        
        # Visibility of transit
        transit_duration_minutes = row['Planet.TransitDuration'] * 60
        self.transit_coverage = len(Visibility.Transit) / transit_duration_minutes
        self.transit = self.transit_coverage > partial
        
        # Baseline visibility
        baseline_length_minutes = baseline_length * 60
        half_baseline_length_minutes = baseline_length_minutes / 2
        
        self.baseline = len(Visibility.Baseline) > baseline_length_minutes
        self.baseline_ingress = len(Visibility.Ingress) > half_baseline_length_minutes
        self.baseline_egress = len(Visibility.Egress) > half_baseline_length_minutes
        
        self.baseline_coverage = min(len(Visibility.Baseline) / baseline_length_minutes, 1.0)
        self.baseline_ingress_coverage = min(len(Visibility.Ingress) / half_baseline_length_minutes, 1.0)
        self.baseline_egress_coverage = min(len(Visibility.Egress) / half_baseline_length_minutes, 1.0)
        
        if self.baseline_coverage > partial:
            self.baseline = True
        return None
#%% Flags for a given window
@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class _Flags_window:
    """Flags for given transit window in terms of feasibility."""
    moon_angle_warning: bool = False # 30 deg < separation < 45 deg
    moon_angle_critical: bool = False # separation < 30 deg
    
    def check(self,
              FlagsTransitFull: _Flags_transit,
              FlagsTransitTwilight: _Flags_transit,
              MoonSeparation: float
              ):
        """
        Check the flags values.

        Parameters
        ----------
        FlagsTransitFull : _Flags_transit
            Flags for full night observation.
        FlagsTransitTwilight : _Flags_transit
            Flags for twilight-included observation.
        MoonSeparation : float
            Minimal separation between target and Moon.

        Returns
        -------
        None.

        """
        
        if FlagsTransitFull.transit and FlagsTransitFull.baseline:
            self.visible = True
        else:
            self.visible = False
            
        if FlagsTransitTwilight.transit and FlagsTransitTwilight.baseline:
            self.visible_twilight = True
        else:
            self.visible_twilight = False
        
        if MoonSeparation <30*u.deg:
            self.moon_angle_critical = True
            self.moon_angle_warning = False
        elif MoonSeparation < 45 * u.deg:
            self.moon_angle_critical = False
            self.moon_angle_warning = True
        else:
            self.moon_angle_critical = False
            self.moon_angle_warning = False
        return

#%% Transit plots


@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class Event(WindowPlot):
    """Class holding the settings for a single transit window plot."""
    
    Location: tel.Telescope
    TransitDuration: float
    baseline: u.Quantity
    Night: astime.Time
    Uncertainty: u.Quantity
    partial: bool | float
    directory: str
    row: pd.Series
    Airmass_limit: float
    TimeArray: _TimeArray = field(default_factory=_TimeArray)
    AltitudeArray: _AltitudeArray = field(default_factory=_AltitudeArray)
    VisibilityFull: _Visibility = field(default_factory=_Visibility)
    VisibilityTwilight: _Visibility = field(default_factory=_Visibility)
    Indices: _Indices = field(default_factory=_Indices)
    FlagsFull: _Flags_transit = field(default_factory=_Flags_transit)
    FlagsTwilight: _Flags_transit = field(default_factory=_Flags_transit)
    FlagsWindow: _Flags_window = field(default_factory=_Flags_window)
    ObservationsFull: _Observations = field(default_factory=_Observations)
    ObservationsTwilight: _Observations = field(default_factory=_Observations)
    velocity_offset: None | float = None
    velocity_range: float = 5
    save_figures: bool = True
    
    
    def __post_init__(self):
        """
        Automatic setup of the plot and quality estimate.

        Returns
        -------
        None.

        """
        self._define_arrays()
        self._calculate_BERV_offset()
        self._calculate_baseline_observations()
        self._define_observability()
        self._define_flags()
        self._estimate_quality_of_transit_window()
        
        if self.quality > 0:
            fig, _ = self.generate_plots()
            directory = os.path.join(self.directory, self.Location.name, self.row['Planet.Name'].replace(' ',''))
            os.makedirs(directory, exist_ok=True)
            self.observing_night = (self.TimeArray.midnight-1*u.day)
            
            
            filepath = os.path.join(directory, f"Q{self.quality}_{self.observing_night.strftime('%Y%m%d')}_{self.row['Planet.Name'].replace(' ','')}.png")
            
            if self.save_figures:
                fig.savefig(filepath)
            logger.info(f"        Transit on: {self.observing_night.strftime('%Y%m%d')} - Quality: {self.quality}; Uncertainty: {self.Uncertainty:.2f}")
        
        return
    
    def _calculate_BERV_offset(self,
                               ):
        if self.velocity_offset is None:
            return
        
        BERV = []
        from PyAstronomy import pyasl
        for time in self.TimeArray.time_array:
            # BERV.append(pyasl.baryvel(time.jd, deq=2000.0))
            BERV.append(pyasl.baryCorr(time.jd,
                                       self.row['Position.RightAscension'],
                                       self.row['Position.Declination'],
                                       deq=2000.0)[1]
                        )
        self.BERV = np.asarray(BERV)
        return
    
    def _calculate_baseline_observations(self):
        """Calculate how long baseline we need."""
        self.baseline_length = self.baseline
        return
    
    def _define_arrays(self):
        """Define the arrays for the plot."""
        self.TimeArray._define_times(Night = self.Night,
                                    TransitDuration = self.TransitDuration)
        
        self.AltitudeArray._define_targets(self.TimeArray,
                                           self.Location,
                                           self.row
                                           )
        
        self.Indices._define_indices(self.TimeArray,
                                     self.AltitudeArray,
                                     self.Airmass_limit,
                                     )
        full_defining_indices = np.argwhere(self.AltitudeArray.Sun.alt.value < -18)
        
        self.VisibilityFull._define_visibility(NightDefiningIndices = full_defining_indices,
                                               Indices = self.Indices)
        
        twilight_defining_indices = np.argwhere(
            (self.AltitudeArray.Sun.alt.value < 0) * self.AltitudeArray.Sun.alt.value != 0
            )
        self.VisibilityTwilight._define_visibility(NightDefiningIndices = twilight_defining_indices,
                                                   Indices = self.Indices)
        
        return
    
    def _define_flags(self):
        """
        Define the flags to estimate transit window quality.

        Returns
        -------
        None.

        """
        self.FlagsFull.check(Visibility = self.VisibilityFull,
                             row = self.row,
                             baseline_length = self.baseline_length,
                             partial= self.partial
                             )
        self.FlagsTwilight.check(Visibility = self.VisibilityTwilight,
                                 row = self.row,
                                 baseline_length = self.baseline_length,
                                 partial= self.partial
                                 )
        self.FlagsWindow.check(
            FlagsTransitFull= self.FlagsFull,
            FlagsTransitTwilight= self.FlagsTwilight,
            MoonSeparation= self.AltitudeArray.MoonSeparation)
        
        return
    
    def _define_observability(self):
        
        self.ObservationsFull._define_observations(
            Visibility= self.VisibilityFull,
            baseline_length= self.baseline_length
            )
        self.ObservationsTwilight._define_observations(
            Visibility= self.VisibilityTwilight,
            baseline_length= self.baseline_length
            )
        return
    
    def _estimate_quality_of_transit_window(self):
        """
        Estimate quality of transit window.
        
        The flagging works like this:
        If not visible (even in twilight) quality = -9999 (ignored completely)
        
        If not visible during full night quality = 5
        
        Else quality = 1
        
        Afterwards, multiple increments are added in case of bad additional condition.
        

        Returns
        -------
        None.

        """
        if not(self.FlagsWindow.visible):
            if not(self.FlagsWindow.visible_twilight):
                self.quality = -9999
                return
            else:
                self.quality = 2
        else:
            self.quality = 1
        
        if self.FlagsWindow.moon_angle_critical:
            self.quality += 10
        else:
            if self.FlagsWindow.moon_angle_warning:
                self.quality += 3
        
        if self.FlagsFull.transit_coverage > 0.8:
            pass    
        elif self.FlagsFull.transit_coverage > 0.5:
            self.quality += 2
        else:
            self.quality += 4
            
        if self.FlagsFull.baseline_coverage > 0.8:
            pass
        elif self.FlagsFull.baseline_coverage > 0.5:
            self.quality += 2
        else:
            self.quality += 4
        
        if self.FlagsFull.baseline_ingress_coverage > 0.5 and self.FlagsFull.baseline_egress_coverage > 0.5:
            pass
        else:
            self.quality += 1
        return


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 10:35:50 2023

@author: chamaeleontis
"""
#%% Importing libraries
import sys
import pandas as pd
from PTO.PTO.transits.NASA_exo import NASA_Exoplanets_archive, logger
import datetime
import astropy.units as u
import numpy as np
import astropy.time as astime
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
from PTO.PTO.transits.observatories import Observatories
import re
import matplotlib.dates as md 
import os

#%%
def define_baseline(time_t14:float) -> float:
    """
    Define baseline based on input time.

    Parameters
    ----------
    time_t14 : float
        Transit duration length in hours.

    Returns
    -------
    baseline_length : float
        Length of baseline.

    """
    baseline_length = np.min([3, .75*time_t14])
    baseline_length = np.max([baseline_length, 2])
    return baseline_length

#%% Class definitions
#%% Planet transits
@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class _Planet_transits:
    """Class holding all the transit windows and parameters for convenience."""
    name: str
    transit_windows: list
    parameters: pd.core.series.Series


@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class _TimeArray:
    """Class holding timings of the transit."""
    def _define_times(self,
                      Night,
                      Planet_transit)->None:
        """Define basic times and time arrays for the plot."""
        self.midnight = astime.Time(round(Night.value) - 0.5,
                                    format=Night.format)
        
        self.time_array = np.linspace(-12,12,24*60+1)*u.hour + self.midnight
        self.T1, self.T4 = ((Night - Planet_transit.parameters['pl_trandur']*u.hour/2) ,
                            (Night + Planet_transit.parameters['pl_trandur']*u.hour/2))
        return None
    
    def _define_sunset_and_sunrise(self,
                                   sunset,
                                   sunrise)->None:
        self.Sunset = sunset
        self.Sunrise = sunrise
        

@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class _Indices:
    """Class holding indices of various arrays."""
    def _define_indices(self,
                        TimeArray,
                        AltitudeArray)->None:
        """Define the indices of transit, ingress, egress and out-of-transit."""
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
                AltitudeArray.target.secz < 2.2,
                AltitudeArray.target.secz >= 1
                )
            )

class _Visibility:
    """Class holding the indices when baseline, target, transit, ingress and egress is visible."""
    def _define_visibility(self,
                           NightDefiningIndices,
                           Indices
                           )->None:
        """Define the visibility depending on full/twilight mode and Indices of transit."""
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
    """Class holding the indices for the observations."""
    def _define_observations(self,
                             Visibility: _Visibility,
                             baseline_length: float
                             ):
        baseline_length = int(round(baseline_length*60/2))
        
        if len(Visibility.Ingress) > baseline_length:
            self.Ingress = Visibility.Ingress[-baseline_length:]
            missing_length_ingress = 0
        else:
            self.Ingress = Visibility.Ingress
            missing_length_ingress = baseline_length - len(Visibility.Ingress)
            
        if len(Visibility.Egress) > baseline_length:
            self.Egress = Visibility.Egress[:baseline_length]
            missing_length_egress = 0
        else:
            self.Egress = Visibility.Egress
            missing_length_egress = baseline_length - len(Visibility.Egress)
        
        if missing_length_ingress != 0 and missing_length_egress != 0:
            self.Ingress, self.Egress = Visibility.Ingress, Visibility.Egress
        elif missing_length_ingress != 0:
            if len(Visibility.Egress) > (baseline_length+missing_length_egress):
                self.Egress = Visibility.Egress[:(baseline_length+missing_length_ingress)]
            else:
                self.Egress = Visibility.Egress
        elif missing_length_egress != 0:
            if len(Visibility.Ingress) > (baseline_length+missing_length_ingress):
                self.Ingress = Visibility.Ingress[-(baseline_length+missing_length_egress):]
            else:
                self.Ingress = Visibility.Ingress
        
        self.Transit = Visibility.Transit
        self.complete = np.concatenate([self.Ingress,self.Transit, self.Egress])

        return


@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class _AltitudeArray:
    """Class holding alitudes of Sun, Moon and target."""
    
    def _define_targets(self,
                        TimeArray,
                        Location:str,
                        Planet_transit,
                        )->None:
        """Define the altitude azimuth frame and calculates target, Moon and Sun position in it."""
        self.altitude_azimuth_frame = coord.AltAz(
            obstime = TimeArray.time_array,
            location = Observatories[Location].value
            )
        self.Sun = coord.get_sun(TimeArray.time_array
                                 ).transform_to(
                                     self.altitude_azimuth_frame)
        self.Moon = coord.get_body("moon",
                                TimeArray.time_array,
                                location = Observatories[Location].value
                                ).transform_to(
                                    self.altitude_azimuth_frame)
        self.target = coord.SkyCoord(
            ra = coord.Angle(Planet_transit.parameters['ra'] * u.deg),
            dec = coord.Angle(Planet_transit.parameters['dec'] * u.deg),
            frame='icrs'
            ).transform_to(
                self.altitude_azimuth_frame)
        self.MoonSeparation = np.min(self.Moon.separation(self.target))
        
        ind = np.argwhere((self.Sun.alt.value < 0) * self.Sun.alt.value != 0)
        TimeArray._define_sunset_and_sunrise(TimeArray.time_array[ind[0]].to_value(format= 'datetime')[0],
                                             TimeArray.time_array[ind[-1]].to_value(format= 'datetime')[0])
        
        return None

#%% Flags transit
@dataclass(init=True, repr=False, eq=False, order=False, unsafe_hash=False, frozen=False)
class _Flags_transit:
    """Flags for given transit opportunity in terms of visibility."""
    def check(self,
              Visibility: _Visibility,
              parameters: pd.DataFrame,
              baseline_length: float,
              partial : bool | float = False
              ):
        """
        Check the values of flags.

        Parameters
        ----------
        Visibility : _Visibility
            Visibility of target's transit. Can be full night or twilight-included version.
        parameters : pd.DataFrame
            Planetary parameters.
        baseline_length : float
            Baseline length.

        Returns
        -------
        None.

        """
        if partial == False:
            partial = 1
        
        # Visibility of transit
        if (len(Visibility.Transit) > (parameters['pl_trandur'] * 60)-1):
            self.transit = True
            self.transit_coverage = 1
        else:
            self.transit = False
            self.transit_coverage = ((len(Visibility.Transit))/ (parameters['pl_trandur'] * 60))
            
            if self.transit_coverage > partial:
                self.transit = True
            
        
        # Baseline visibility
        self.baseline = (len(Visibility.Baseline) > baseline_length * 60)
        self.baseline_ingress = (len(Visibility.Ingress) > baseline_length/2*60)
        self.baseline_egress = (len(Visibility.Egress) > baseline_length/2*60)
        
        self.baseline_coverage = np.min([( len(Visibility.Baseline)) / (baseline_length * 60),1])
        self.baseline_ingress_coverage = np.min([( len(Visibility.Ingress) / (baseline_length/2 * 60)),1])
        self.baseline_egress_coverage = np.min([(len(Visibility.Egress) / (baseline_length/2 * 60) ),1])
        
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
class Plot_transit:
    """Class holding the settings for a single transit window plot."""
    
    Location: str
    Planet_transit: _Planet_transits
    Night: astime.Time
    Uncertainty: u.quantity.Quantity
    partial: bool | float
    TimeArray: _TimeArray = _TimeArray()
    AltitudeArray: _AltitudeArray = _AltitudeArray()
    VisibilityFull: _Visibility = _Visibility()
    VisibilityTwilight: _Visibility = _Visibility()
    Indices: _Indices = _Indices()
    FlagsFull: _Flags_transit = _Flags_transit()
    FlagsTwilight: _Flags_transit = _Flags_transit()
    FlagsWindow: _Flags_window = _Flags_window()
    ObservationsFull: _Observations = _Observations()
    ObservationsTwilight: _Observations = _Observations()
    def __post_init__(self):
        """
        Automatic setup of the plot and quality estimate.

        Returns
        -------
        None.

        """
        self._define_arrays()        
        self._calculate_baseline_observations()
        self._define_observability()
        self._define_flags()
        self._estimate_quality_of_transit_window()
        
        if self.quality > 0:
            self._log_transit_info()
            self._create_plot()
        
        return
    
    def _create_plot(self):
        """
        Create transit window plot.

        Returns
        -------
        None.

        """
        with plt.ioff():
            fig, ax = plt.subplots(1,figsize=(18,12))
            self._plot_basic_setting(fig, ax)
            self._plot_twilight(ax)
            self._plot_Moon(ax)
            self._plot_target(ax)
            self._plot_observation(ax)
            self._write_system_parameters(ax)
            
            self.fig, self.ax = fig, ax
        return
    
    def _log_transit_info(self):
        """
        Log transit information and quality.

        Returns
        -------
        None.

        """
        logger.info('Transit during night: ' + str(self.TimeArray.midnight.to_value('datetime').strftime('%Y%m%d')) + ' is of quality: '+ str(self.quality) +' with uncertainty of {sigma:.2f} min'.format(sigma = self.Uncertainty.to(u.min).value))
        if self.quality != 1:
            logger.warning('This is a degraded transit because:')
            if self.FlagsWindow.moon_angle_critical:
                logger.warning('    Moon separation angle too low')
            if self.FlagsWindow.moon_angle_warning:
                logger.warning('    Moon separation angle warning (30 deg < angle < 45 deg)')
            if (not(self.FlagsWindow.visible) and
                self.FlagsWindow.visible_twilight):
                logger.warning('    Transit window is heavily visible inside the twilight.')
            if self.FlagsFull.baseline_ingress_coverage < 0.9:
                logger.warning(f'    Ingress baseline is not full with coverage {self.FlagsFull.baseline_ingress_coverage:.0%} during full night')
                logger.info(f'    Ingress baseline including twilight coverage:{self.FlagsTwilight.baseline_ingress_coverage:.0%}')
            if self.FlagsFull.baseline_egress_coverage < 0.9:
                logger.warning(f'    Egress baseline is not full with coverage {self.FlagsFull.baseline_egress_coverage:.0%} during full night')
                logger.info(f'    Egress baseline including twilight coverage: {self.FlagsTwilight.baseline_egress_coverage:.0%}')
            if self.FlagsFull.transit_coverage < 0.9:
                logger.warning(f'    Transit event is not full with coverage {self.FlagsFull.transit_coverage:.0%} during full night')
                logger.info(f'    Transit event including twilight coverage: {self.FlagsTwilight.transit_coverage:.0%}')
        
    def _plot_target(self,
                     ax:plt.Axes):
        """
        Plot the target, highlighting the transit area.

        Parameters
        ----------
        ax : plt.Axes
            Artist on which to plot.

        Returns
        -------
        None.

        """
        ax.plot(self.TimeArray.time_array.to_value(format= 'datetime'),
                self.AltitudeArray.target.alt.value,
                ls = '-',
                label='Target',
                linewidth=4,
                color='red',
                )
        ax.plot(self.TimeArray.time_array[self.Indices.Transit].to_value(format= 'datetime'),
                self.AltitudeArray.target[self.Indices.Transit].alt.value,
                ls = '-',
                label='Target',
                linewidth=10,
                color='red',
                )
        
        ax.scatter(
            self.TimeArray.time_array[self.Indices.Transit[0]].to_value(format= 'datetime'),
            self.AltitudeArray.target[self.Indices.Transit[0]].alt.value,
            color = 'white',
            )
        ax.scatter(
            self.TimeArray.time_array[self.Indices.Transit[-1]].to_value(format= 'datetime'),
            self.AltitudeArray.target[self.Indices.Transit[-1]].alt.value,
            color = 'white',
            )
        
        return
    
    def _plot_Moon(self,
                   ax:plt.Axes):
        """
        Overplot the Moon position.

        Parameters
        ----------
        ax : plt.Axes
            Artist on which to plot on.

        Returns
        -------
        None.

        """
        ax.plot(self.TimeArray.time_array.to_value(format= 'datetime'),
                self.AltitudeArray.Moon.alt.value,
                ls = '-',
                label='Moon',
                linewidth=2,
                color='yellow',
                alpha=0.7
                )
        return
    
    def _plot_twilight(self,
                       ax:plt.Axes):
        """
        Overplot the twilight area.

        Parameters
        ----------
        ax : plt.Axes
            Artist on which to plot on.

        Returns
        -------
        None.

        """
        ax.fill_between(self.TimeArray.time_array.to_value(format= 'datetime'),
                        0, 90,
                        self.AltitudeArray.Sun.alt.value < 0,
                        color='blue',
                        alpha = 0.3
                        )
        ax.fill_between(self.TimeArray.time_array.to_value(format= 'datetime'),
                        0, 90,
                        self.AltitudeArray.Sun.alt.value < -18,
                        color ='blue',
                        alpha= 0.5
                        )
        return

    def _plot_observation(self,
                          ax):
        
        
        #CLEANME
        if len(self.ObservationsFull.complete) != 0:
            ax.fill_between(self.TimeArray.time_array[self.ObservationsFull.complete[0]:
                                            self.ObservationsFull.complete[-1]].to_value(
                                                format= 'datetime'),
                            0, 10,
                            color='lime',
                            alpha = 1
                            )
            ax.text((self.TimeArray.time_array[self.ObservationsFull.complete[0]]+30*u.min).to_value(
                format= 'datetime'),
                    5,
                    'Observations (no twilight): ' + 
                    self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(
                        format= 'datetime').strftime('%H:%M') + ' - ' +
                    self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(
                        format= 'datetime').strftime('%H:%M'),
                    fontsize = 10,
                    )
            ax.text(self.TimeArray.midnight.to_value(
                format= 'datetime'),
                    -10,
                    'Observations (no twilight): ' + 
                    self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(
                        format= 'datetime').strftime('%H:%M') + ' - ' +
                    self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(
                        format= 'datetime').strftime('%H:%M'),
                    fontsize = 15,
                    bbox=dict(
                        facecolor='green',
                        edgecolor='black',
                        boxstyle='round,pad=1'
                        )
                    )
            
            self.full_observation_time = [
                self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(format= 'datetime'),
                self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(format= 'datetime')
                ]
        
        self.twilight_observation_time = None
        
        if ((self.FlagsFull.baseline_ingress) or # Ingress baseline missing
            (self.FlagsFull.baseline_egress) or # Egress baseline missing
            (not(self.FlagsWindow.visible) and self.FlagsWindow.visible_twilight)): # Transit into twilight
            
            
            
            ax.fill_between(self.TimeArray.time_array[self.ObservationsTwilight.complete[0]:
                                            self.ObservationsTwilight.complete[-1]].to_value(
                                                format= 'datetime'),
                            10, 20,
                            color='goldenrod',
                            alpha = 1
                            )

            ax.text((self.TimeArray.time_array[self.ObservationsTwilight.complete[0]]+30*u.min).to_value(
                format= 'datetime'),
                    15,
                    'Observation (with twilight): ' +
                    self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(
                             format= 'datetime').strftime('%H:%M') + ' - ' +
                    self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(
                             format= 'datetime').strftime('%H:%M'),
                    fontsize = 10,
                    )
            
            ax.text(self.TimeArray.midnight.to_value(
                format= 'datetime'),
                    -20,
                    'Observation (with twilight): ' +
                    self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(
                             format= 'datetime').strftime('%H:%M') + ' - ' +
                    self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(
                             format= 'datetime').strftime('%H:%M'),
                    fontsize = 15,
                    bbox=dict(
                        facecolor='orange',
                        edgecolor='black',
                        boxstyle='round,pad=1'
                        )
                    )
            
            ax.axvline(self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(format= 'datetime'), ls='--', lw=2, color='orange')
            ax.axvline(self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(format= 'datetime'), ls='--', lw=2, color='orange')
            self.twilight_observation_time = [
                self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(format= 'datetime'),
                self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(format= 'datetime')
                ]
        
        ax.axvline(self.TimeArray.T1.to_value(format= 'datetime'), ls='--', lw=2, color='darkred')
        ax.axvline(self.TimeArray.T4.to_value(format= 'datetime'), ls='--', lw=2, color='darkred')
        if len(self.ObservationsFull.complete) != 0:
            ax.axvline(self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(format= 'datetime'), ls='--', lw=2, color='lime')
            ax.axvline(self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(format= 'datetime'), ls='--', lw=2, color='lime')
        
        
        
        return
    
    def _plot_basic_setting(self,
                            fig: plt.Figure,
                            ax:plt.Axes):
        """
        Plot basic settings regardless of what planet is plotted.

        Parameters
        ----------
        fig : plt.Figure
            Figure on which to plot on.
        ax : plt.Axes
            Artist on which to plot on.

        Returns
        -------
        None.

        """
        Night = (self.TimeArray.midnight-1*u.day).to_value('datetime')
        ax.set_title('Transit of '+ str(self.Planet_transit.parameters['pl_name'])+
                     ' on night: '+str(Night.strftime('%Y%m%d')) + 
                     ' at Paranal with precision on Tc of: ' + 
                     '{sigma:.2f} min'.format(sigma = self.Uncertainty.to(u.min).value) +
                     '; Quality of transit:' + str(self.quality),
                     fontsize = 18
                     )

        ax.axhline(
            np.arcsin(1/2.2)* 180/np.pi,
            ls='--', lw=2,
            color= 'darkred'
            )
        ax.grid(color='black', linestyle='dashed', linewidth=0.3, alpha=1.0)
        xfmt = md.DateFormatter('%d/%m  %H')
        ax.xaxis.set_major_formatter(xfmt)
        fig.autofmt_xdate()
        hours = md.HourLocator()
        ax.xaxis.set_major_locator(hours)
        
        ax.set_ylim(0,90)
        ax2 = ax.twinx()
        ax.set_ylabel('Altitude [deg]', color='k')
        ax2.set_ylabel('Airmass', color='k')
        ax2.tick_params('y', colors='k')
        ax2.set_ylim([0,90])
        ax2.set_yticklabels(['',5.76,2.92,2.00,1.56,1.31,1.15,1.06,1.02,1.00]) 
        
        ax.set_xlim(self.TimeArray.Sunset,
                    self.TimeArray.Sunrise)
        
        
        return
    
    def _write_system_parameters(self,
                                 ax: plt.Axes):
        """
        Write the system parameters for given system, including Right Ascension, Declination,V magnitude, Moon angle separation, Transit length, Period, ingress and egress timings.

        Parameters
        ----------
        ax : plt.Axes
            DESCRIPTION.

        Returns
        -------
        None.

        """
        ax.text(self.TimeArray.Sunset,
                97,
                ('RA: {ra:.2f}'.format(ra=self.Planet_transit.parameters['ra']) + 
                '; DEC: {dec:.2f}'.format(dec=self.Planet_transit.parameters['dec']) + \
                '; V = {V_mag:.3f}'.format(V_mag = self.Planet_transit.parameters['sy_vmag']) + \
                '; Moon = {moon:.0f} deg'.format(moon = self.AltitudeArray.MoonSeparation.value) + \
                '; $T_{14}$' + ' = {t14:.3f} h'.format(t14 = self.Planet_transit.parameters['pl_trandur']) + \
                '; P = {p:.5f} d'.format(p=self.Planet_transit.parameters['pl_orbper'])
                ),
                style='italic',size=12
                )
        
        ax.text(self.TimeArray.Sunset,
                -23,
                "ingress=" + self.TimeArray.T1.to_value('datetime').strftime("%H:%M"),
                style='italic',
                bbox=dict(facecolor='none', edgecolor='black', pad=7)
                )
        
        ax.text(self.TimeArray.Sunrise,
                -23,
                "egress=" + self.TimeArray.T4.to_value('datetime').strftime("%H:%M"),
                style='italic',
                bbox=dict(facecolor='none', edgecolor='black', pad=7)
                )
        return
    
    def _calculate_baseline_observations(self):
        """Calculate how long baseline we need."""
        self.baseline_length = define_baseline(self.Planet_transit.parameters['pl_trandur'])
        return
    
    def _define_arrays(self):
        """Define the arrays for the plot."""
        self.TimeArray._define_times(Night = self.Night,
                                    Planet_transit = self.Planet_transit,)
        
        self.AltitudeArray._define_targets(self.TimeArray,
                                           self.Location,
                                           self.Planet_transit
                                           )
        
        self.Indices._define_indices(self.TimeArray,
                                     self.AltitudeArray)
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
                             parameters = self.Planet_transit.parameters,
                             baseline_length = self.baseline_length,
                             partial= self.partial
                             )
        self.FlagsTwilight.check(Visibility = self.VisibilityTwilight,
                                 parameters = self.Planet_transit.parameters,
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
                self.quality = 5
        else:
            self.quality = 1
        
        if self.FlagsWindow.moon_angle_critical:
            self.quality += 10
        else:
            if self.FlagsWindow.moon_angle_warning:
                self.quality += 1
        
        if not(self.FlagsFull.baseline):
            if self.FlagsFull.baseline_coverage > 0.8:
                self.quality += 1
            else:
                if self.FlagsTwilight.baseline_coverage != 1:
                    self.quality += 2
                else:
                    self.quality += 1
        
        if self.FlagsFull.baseline == True:
            if ((self.FlagsFull.baseline_ingress_coverage < 1) or
                (self.FlagsFull.baseline_egress_coverage < 1)):
                if ((self.FlagsTwilight.baseline_ingress_coverage < 1) or
                    (self.FlagsTwilight.baseline_egress_coverage < 1)):
                    self.quality += 2
                else:
                    self.quality += 1
        elif (self.FlagsFull.baseline_coverage > 0.8):
            self.quality += 2
        else:
            self.quality += 4
        
        
        return

    def save_transit_window(self,
                            location:str
                            ):
        """
        Save the transit window plot at given location.

        Parameters
        ----------
        location : str
            Filepath to the saved figure.

        Returns
        -------
        None.

        """
        self.filename = location
        self.fig.savefig(location)
        plt.close('all')
        return
    
    def output_observability(self):
        
        # Full observability
        # Planned observability
        
        return

#%%
class Transit_windows:
    """Class holding all transit windows calculated."""
    def __init__(
            self,
            database: NASA_Exoplanets_archive,
            semester_date: astime.Time,
            save_directory: str,
            ):
        """
        Initialize the class Transit windows. 

        Parameters
        ----------
        database : NASA_Exoplanets_archive
            Which database we take data from.
        semester_date : astime.Time
            Starting and ending date of the observations.
        save_directory : str
            Save directory for the figures.

        Returns
        -------
        None.

        """
        self.planets = []
        
        filtered_table = database.filtered_table_all
        filtered_table = self._check_parameters_in_table(filtered_table)
        
        filtered_table_composite = database.filtered_table_composite
        filtered_table_composite = self._check_parameters_in_table(filtered_table_composite,
                                                                   True)
        filtered_table_all = self._add_transit_length_if_missing(filtered_table_all= filtered_table,
                                                                 filtered_table_composite= filtered_table_composite)
        
        
        self.calculate_transit_windows(
            exoplanet_table = filtered_table,
            date = semester_date
            )
        
        self.calculate_transit_windows(
            exoplanet_table = filtered_table_composite,
            date = semester_date
            )
        
        self.save_directory = save_directory
        
        return
    
    def _check_parameters_in_table(self,
                                   filtered_table:pd.DataFrame,
                                   composite:bool = False) -> pd.DataFrame:
        """
        Filter out references without full ephemeris parameters.

        Parameters
        ----------
        filtered_table : pd.DataFrame
            Filtered table from the NASA table.
        composite : bool, optional
            Whether we use composite/complete tables. The default is False.

        Returns
        -------
        filtered_table : pd.DataFrame
            Filtered table.

        """
        original_length = len(filtered_table)
        condition_composite = (filtered_table['pl_orbper'].notna() &
                               filtered_table['pl_orbpererr1'].notna() &
                               filtered_table['pl_orbpererr2'].notna() &
                               filtered_table['pl_tranmid'].notna() &
                               filtered_table['pl_tranmiderr1'].notna() &
                               filtered_table['pl_tranmiderr2'].notna()
                               )
        dropped_targets = filtered_table[~condition_composite]
        
        # WHY IS NASA COMPOSITE TABLE INCLUDING VALUES WITHOUT ERRORBARS AS MOST PRECISE!!!!!
        filtered_table = filtered_table[condition_composite]
        
        
        if composite:
            filtered_table = filtered_table.loc[:,
                                                [
                                                'pl_name',
                                                'pl_orbper', 
                                                'pl_orbpererr1', 
                                                'pl_orbpererr2', 
                                                'pl_trandur',
                                                'pl_trandurerr1',
                                                'pl_trandurerr2',
                                                'pl_tranmid',
                                                'pl_tranmiderr1',
                                                'pl_tranmiderr2',
                                                'ra',
                                                'dec',
                                                'sy_vmag'
                                                ]
                                                ]
            filtered_table['pl_refname'] = ['Composite']*len(filtered_table)
            if len(filtered_table) != original_length:
                logger.critical('Some of the planet rows were dropped because not all necessary values were found in the table. Dropped targets:')
                for ind, row in dropped_targets.iterrows():
                    name = row['pl_name']
                    logger.critical(f'    {name}')
        else:
            filtered_table.loc[:, 'pl_refname'] = [re.sub('<[^<]+?>', '', ref_name) for ref_name in filtered_table['pl_refname'].values] 
            
        return filtered_table
    
    def _add_transit_length_if_missing(self,
                                       filtered_table_all: pd.DataFrame,
                                       filtered_table_composite:pd.DataFrame) -> pd.DataFrame:
        """
        Add transit length to rows if it is missing using the composite table.

        Parameters
        ----------
        filtered_table_all : pd.DataFrame
            Filtered NASA Exoplanet Full table.
        filtered_table_composite : pd.DataFrame
            Filtered NASA Exoplanet Composite table.

        Returns
        -------
        pd.DataFrame
            Updated table with the missing transit lengths.
        """
        
        
        condition_no_transit_length = filtered_table_all['pl_trandur'].isna()
        
        
        for ind, row in filtered_table_all.loc[condition_no_transit_length].iterrows():
            pl_name = row['pl_name']
            composite_row = filtered_table_composite[filtered_table_composite['pl_name'] == pl_name]
            T14 =  composite_row['pl_trandur'].values
            if not(np.isfinite(T14)):
                logger.warning(f'Transit length not found in composite table for planet {pl_name}')
                continue
            
            filtered_table_all[filtered_table_all['pl_name'] == pl_name].fillna(value = T14[0])
            
        
        return filtered_table_all
    
    def _calculate_transit_windows_for_single_planet(self,
                                                     row: pd.core.series.Series,
                                                     date: astime.Time
                                                     ) -> list:
        """
        Calculate transit windows for single planet.

        Parameters
        ----------
        row : pd.core.series.Series
            Row of pd.DataFrame holding all the information on the planet.
        date : astime.Time
            Semester in which we want to observe.

        Returns
        -------
        transit_windows : list
            List of all posible transit window.

        """
        P = row['pl_orbper'] * u.day
        P_error = np.max([row['pl_orbpererr1'], row['pl_orbpererr2']])
        T_C = astime.Time(row['pl_tranmid'], format = 'jd')
        T_C_error = np.max([row['pl_tranmiderr1'], row['pl_tranmiderr2']])
        transit_windows = []
        n = 0
        
        while date[0] < T_C: # To ensure we search from beginning. This will add some uncertainties if searching for past windows, but that is generally irrelevant in the purpose of this code case
            T_C -= P
            n +=1
            logger.warning('Please use T_C that is less than the starting date of observations. Otherwise, the error propagation will be overestimated.')

        while date[1] > T_C:
            T_C += P
            n +=1
            if date[0] < T_C and date[1] > T_C:
                sigma_T0 = np.sqrt(T_C_error**2 + n**2 * P_error**2)
                transit_windows.append([T_C, sigma_T0*u.day])
                
        return transit_windows
        
    
    def calculate_transit_windows(self,
                                  exoplanet_table: pd.DataFrame,
                                  date: astime.Time
                                  ):
        """
        Calculate all transit windows for all planet/references within given period.

        Parameters
        ----------
        exoplanet_table : pd.DataFrame
            Table of all references to use.
        date : astime.Time
            When to consider transit windows.

        Returns
        -------
        None.

        """
        for ind, row in exoplanet_table.iterrows():
            transit_windows = self._calculate_transit_windows_for_single_planet(
                    row = row,
                    date = date
                    )
            
            
            self.planets.append(
                _Planet_transits(
                    name = row['pl_name'],
                    transit_windows = transit_windows,
                    parameters = row
                    )
                )
        return
    
    
    def _print_info_about_window(self,
                                 Planet_transit: _Planet_transits
                                 ):
        
        logger.info('Currently working on %s'%(Planet_transit.name) + ' with reference %s'%(Planet_transit.parameters['pl_refname']))
        logger.info('P = {period:.6f} d'.format(period = Planet_transit.parameters['pl_orbper']) +
                    u" \u00B1 " + '{period_error:.6f}'.format(
                        period_error = max([Planet_transit.parameters['pl_orbpererr1'],
                                            Planet_transit.parameters['pl_orbpererr2']])
                        )
                    )
        logger.info('T14 = {t14:.3f} h'.format(t14 = Planet_transit.parameters['pl_trandur'])  +
                    u" \u00B1 " + '{t14_error:.3f}'.format(
                        t14_error = max([Planet_transit.parameters['pl_trandurerr1'],
                                            Planet_transit.parameters['pl_trandurerr2']])
                        )
                    )
        logger.info('Tc ' + '= {tc:.6f} [JD]'.format(tc = Planet_transit.parameters['pl_tranmid']) +
                    u" \u00B1 " + '{tc_error:.6f}'.format(
                        tc_error = max([Planet_transit.parameters['pl_tranmiderr1'],
                                        Planet_transit.parameters['pl_tranmiderr2']])
                        )
                    )
        logger.info('='*60)
            
        
        return

    def plot_transit_windows(self,
                             Location: str,
                             partial: bool | float = False,
                             ):
        """
        Creates plots of all valid transit windows provided Location, assuming the filtered table and all possible ephemeris sets.

        Parameters
        ----------
        Location : str
            Location where we wish to observe. Must be implemented in the observatories.Observatories class
        partial : bool
            Whether include partial transits in the planning.

        Returns
        -------
        None.

        """
        logger.info('='*60)
        
        # if partial:
        #     logger.critical('Partial windows are not implemented yet!')
        self._order_list()

        
        for Planet_transit in self.planets:
            if len(Planet_transit.transit_windows) == 0:
                logger.warning('No transit window found for planet %s'%Planet_transit.parameters['pl_name'])
                continue
            self._print_info_about_window(Planet_transit)
            best_reference = self._find_best_ephemeris(Planet_transit= Planet_transit)
            
            if best_reference == Planet_transit.parameters['pl_refname']:
                path_file_reference = Location + '/%s'%(Planet_transit.name) + '/best_%s'%(Planet_transit.parameters['pl_refname'])
                txt_summary =  Location + '/%s'%(Planet_transit.name) + '/best_%s'%(Planet_transit.parameters['pl_refname']) + '.csv'
            else:
                path_file_reference = Location + '/%s'%(Planet_transit.name) + '/%s'%(Planet_transit.parameters['pl_refname'])
                txt_summary =  Location + '/%s'%(Planet_transit.name) + '/best_%s'%(Planet_transit.parameters['pl_refname']) + '.csv'
                
            path_file_reference = path_file_reference.replace(' ', '')
            txt_summary = txt_summary.replace(' ', '')
            
            for night, uncertainty in Planet_transit.transit_windows:
                new_Plot = Plot_transit(Location = Location,
                                          Planet_transit = Planet_transit,
                                          Night = night,
                                          Uncertainty = uncertainty,
                                          partial= partial
                                          )
                
                if new_Plot.quality != -9999:
                    semester = self._check_P_ESO(new_Plot.Night)
                    file_name = self.save_directory + '/' + path_file_reference + '/Q' + str(new_Plot.quality) + '_'  + f'{semester}_'+ (new_Plot.Night-1*u.day).to_value('datetime').strftime('%Y%m%d')
                    # Make directory if it doesn't exist
                    os.makedirs(self.save_directory + '/' + path_file_reference, mode = 0o777, exist_ok = True) 
                    new_Plot.save_transit_window(file_name)
                    
                    with open(self.save_directory + '/' +txt_summary, 'a') as f:
                        line = Planet_transit.name + ';' + (new_Plot.Night-1*u.day).to_value('datetime').strftime('%Y%m%d') + ';' + str(new_Plot.quality) + ';' + new_Plot.TimeArray.time_array[new_Plot.ObservationsFull.complete[0]].to_value(
                             format= 'datetime').strftime('%H:%M') + ';' + new_Plot.TimeArray.time_array[new_Plot.ObservationsFull.complete[-1]].to_value(
                             format= 'datetime').strftime('%H:%M') + ';'+ str(int(len(new_Plot.VisibilityFull.Ingress)>30))  + '\n' 
                        f.write(line)
                    
                else:
                    night = (new_Plot.Night-1*u.day).to_value('datetime').strftime('%Y%m%d')
                    logger.warning(f'Skipping transit window {night}')
            logger.info('='*60)
        return
    
    def _find_best_ephemeris(self,
                             Planet_transit: _Planet_transits,
                             ) -> str:
        """
        Finds which of the available ephemeris is the best.

        Parameters
        ----------
        Planet_transit : _Planet_transits
            Planet transits class, including the 

        Returns
        -------
        best_reference : str
            Best reference for given transit window period
        """
        best_reference = None
        for ind, planet in enumerate(self.planets):
            if (best_reference is None and 
                Planet_transit.parameters['pl_name'] == planet.parameters['pl_name']):
                best_reference = planet.parameters['pl_refname']
                uncertainty = planet.transit_windows[0][1]
            elif Planet_transit.parameters['pl_name'] == planet.parameters['pl_name']:
                if planet.transit_windows[0][1] < uncertainty:
                    best_reference = planet.parameters['pl_refname']
                    uncertainty = planet.transit_windows[0][1]
            else:
                continue

        return best_reference
    
    def _order_list(self):
        """
        Order list of planet ephemeris by planet name. This ensures the calculation is done planet by planet.
        """
        
        order = []
        for Planet_transit in self.planets:
            order.append(Planet_transit.name)
            indices = np.argsort(np.array(order))
        new_list = []
        # CLEANME: This is dumb.
        for ind in indices:
            new_list.append(self.planets[ind])
            
        self.planets = new_list
        return
    
    
    def _check_P_ESO(self,
                     Night: astime.Time) -> str:
        """
        Calculates in which period the transit night is to order the transit based on it.,

        Parameters
        ----------
        Night : astime.Time
            Night of the transit

        Returns
        -------
        semester : str
            Which period of ESO proposal call does this window fall in.
        """
        P112 = astime.Time(['2024-04-01 12:00:00.000'],scale = 'utc')
        P113 = astime.Time(['2024-10-01 12:00:00.000'],scale = 'utc')
        semester = None
        ind = 0
        while semester is None:
            if Night > P113:
                P113 += 1*u.year
                P112 += 1*u.year
                ind += 2
            elif Night > P112:
                P112 += 1*u.year
                ind += 1
            else:
                semester = 112+ind

        return str(semester)
    
    
    def create_transit_summary(self,
                               Plots,
                               filename:str
                               ):
        information = {}
        
        for plot in Plots:
            baseline_full = ('between(' +
                             self.full_observation_time[0].strftime('%Y-%m-%dT%H:%M') + ', '+ 
                             self.full_observation_time[1].strftime('%Y-%m-%dT%H:%M') + ', '+
                             str(plot.quality) + ', "Transit of %s; '%(plot.Planet_transit.name) +
                             '{time:.2f}'.format(time = (
                                 plot.full_observation_time[1] - 
                                 plot.full_observation_time[0]).seconds/3600
                                 ) + ' hours of observation"'
                             )
            baseline_twilight = None
            if not(plot.twilight_observation_time is None):
                baseline_twilight = ('between(' +
                                 self.twilight_observation_time[0].strftime('%Y-%m-%dT%H:%M') + ', '+ 
                                 self.twilight_observation_time[1].strftime('%Y-%m-%dT%H:%M') + ', '+
                                 str(plot.quality) + ', "Transit of %s; '%(plot.Planet_transit.name) +
                                 '{time:.2f}'.format(time = (
                                     plot.twilight_observation_time[1] - 
                                     plot.twilight_observation_time[0]).seconds/3600
                                     ) + ' hours of observation"'
                                 )
            
            try:
                information[str(plot.quality)].append({
                    'Night': plot.Night.to_value('datetime').strftime('%Y%m%d'),
                    'Baseline_full': baseline_full,
                    'Baseline_twilight': baseline_twilight,
                    'filename': plot.filename,
                    'quality': plot.quality
                    })
            except:
                information[str(plot.quality)] = [{
                    'Night': plot.Night.to_value('datetime').strftime('%Y%m%d'),
                    'Baseline_full': baseline_full,
                    'Baseline_twilight': baseline_twilight,
                    'filename': plot.filename,
                    'quality': plot.quality
                    }]
            pass
        
        
        
        return
        
        
#%% Get dates for ESO semester counting
def get_dates_for_ESO_semester(
        P_number: int,
        large_program: bool = False
        ) -> astime.Time:
    """
    Provide starting and ending date of the ESO semester.

    Parameters
    ----------
    P_number : int
        Which semester are we looking for.
    large_program : bool, optional
        If we consider a large (4 semesters) program. Will give the 6 semesters by default, as ESO asks for transit windows in two more semesters. The default is False.

    Returns
    -------
    time : astime.Time
        Starting and ending date of the semester.

    """
    P112 = astime.Time(['2023-10-01 12:00:00.000', '2024-04-01 12:00:00.000'],
                       scale = 'utc',)
    P113 = astime.Time(['2024-04-01 12:00:00.000', '2024-10-01 12:00:00.000'],
                       scale = 'utc',)
    P112_large = astime.Time(['2023-10-01 12:00:00.000', '2026-10-01 12:00:00.000'],
                             scale = 'utc',)
    P113_large = astime.Time(['2024-04-01 12:00:00.000', '2027-04-01 12:00:00.000'],
                             scale = 'utc',)
    
    if P_number%2 == 0:
        if large_program:
            time = P112_large + ((P_number - 112)/2)*u.year
        else:
            time = P112 + ((P_number - 112)/2)*u.year
    else:
        if large_program:
            time = P113_large + ((P_number - 113)/2)*u.year
        else:
            time = P113 + ((P_number - 113)/2)*u.year
    return time

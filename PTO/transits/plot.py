import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
import astropy.time as astime
import datetime
import pandas as pd
from ..telescopes.telescopes import Telescope
import astropy.units as u
import numpy as np
import matplotlib.dates as md
import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
if not logger.handlers:
    logger = logger_default(logger)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.dates import date2num
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import astropy.constants as con
from PyAstronomy import pyasl
import numpy as np


def calculate_planetary_rv(time, row):
    """
    Calculates the radial velocity of a planet based on its orbital parameters.

    This function computes the radial velocity of a planet using its orbital parameters and the time of observation. It extracts the necessary parameters from the provided dictionary and performs the calculations using Kepler's laws.

    Parameters
    ----------
    time : astropy.time.Time
        The time at which to calculate the radial velocity.
    row : dict
        A dictionary containing the orbital parameters of the planet and star. Expected keys are:
        - 'Planet.Period': Orbital period of the planet.
        - 'Planet.SemiMajorAxis': Semi-major axis of the planet's orbit.
        - 'Planet.Eccentricity': Eccentricity of the planet's orbit.
        - 'Planet.ArgumentOfPeriastron': Argument of periastron (in degrees).
        - 'Planet.TransitMidpoint': Time of the planet's transit midpoint.
        - 'Star.Mass': Mass of the star.
        - 'Planet.MassJupiter': Mass of the planet in Jupiter masses.
        - 'Planet.Inclination': Inclination of the planet's orbit (in degrees).
        - 'Omega_deg' (optional): Longitude of the ascending node (in degrees). Default is 0.0.

    Returns
    -------
    numpy.ndarray
        The radial velocity of the planet at the given time.
    """
    # Extract parameters with defaults for optional ones
    period = row['Planet.Period']
    semi_major_axis = row['Planet.SemiMajorAxis']
    eccentricity = row['Planet.Eccentricity']
    omega_deg = row['Planet.ArgumentOfPeriastron']
    t_transit = row['Planet.TransitMidpoint']
    m_star = row['Star.Mass']
    m_planet = row['Planet.MassJupiter'] * u.M_jup.to(u.M_sun)
    inclination_deg = row['Planet.Inclination']
    Omega_deg = row.get('Omega_deg', 0.0)
    t = time.jd
    omega = np.radians(omega_deg)
    # Calculate time of periastron from transit time
    f_transit = np.pi/2 - omega
    E_transit = 2 * np.arctan(np.sqrt((1-eccentricity)/(1+eccentricity)) * 
                             np.tan(f_transit/2))
    M_transit = E_transit - eccentricity * np.sin(E_transit)
    dt = M_transit * period / (2*np.pi)
    t_periastron = t_transit - dt

    # Create Keplerian orbit
    ke = pyasl.KeplerEllipse(semi_major_axis, period, e=eccentricity,
                            Omega=Omega_deg, i=inclination_deg, w=omega_deg)

    # Calculate orbital position
    vel = ke.xyzVel(t-t_periastron)[:,2] *u.au/u.day
    vel = vel * (u.au/u.day).to(u.m/u.s)
    # Extract radial velocity (z-component of velocity)
    mass_ratio = m_planet / (m_star + m_planet)
    star_vel = vel * mass_ratio
    planet_vel = star_vel * (-(m_star/m_planet))
    return planet_vel

@dataclass
class WindowPlot():
    def generate_plots(self):
        with plt.ioff():
            fig, ax = plt.subplots(1, figsize=(18, 12))
            self._plot_basic_setting(fig, ax)
            self._plot_twilight(ax)
            self._plot_Moon(ax)
            self._plot_target(ax)
            self._plot_observation(ax)
            self._plot_BERV_spine(ax)
            self._write_system_parameters(ax)
            return fig, ax

    def _plot_BERV_spine(self, ax):
        if self.velocity_offset is None:
            return
        
        try:
            self.row['System.Velocity']
            pass
        except:
            logger.warning('Systemic velocity not found. Setting to 0')
            self.row['System.Velocity'] = 0
        
        if np.isnan(self.row['System.Velocity']):
            self.row['System.Velocity'] = 0
        
        planet_rv = calculate_planetary_rv(self.TimeArray.time_array, self.row)
        
        # Get the top spine's path
        x_data = date2num(self.TimeArray.time_array.to_value(format='datetime'))
        values = abs((-self.BERV + self.row['System.Velocity']) - self.velocity_offset + planet_rv)

        _, ymax = ax.get_ylim()

        # Create points for the line segments
        points = np.array([x_data, np.ones_like(x_data) * ymax+0.05]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Define colors for different ranges
        colors = ['red', 'yellow', 'green']
        cmap = ListedColormap(colors)

        # Define boundaries for color changes
        bounds = [0, self.velocity_range, self.velocity_range*2, 999]
        norm = BoundaryNorm(bounds, cmap.N)

        # Create LineCollection with discrete colors
        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=20)
        lc.set_array(values[:-1])
        
        # Remove original top spine
        ax.spines['top'].set_visible(False)

        # Add colored line collection
        ax.add_collection(lc)
        cbar = plt.colorbar(lc)
        cbar.set_ticks([2.5, 7.5, (bounds[-1]-bounds[-2])/2])  # Center of each color range
        cbar.set_ticklabels(['Contamination in range', 'Contamination in 2x range', 'No contamination'])
        
        # raise NotImplementedError('Implement planet velocity properly')
        return 

    def _plot_contamination_in_transit(self, ax):
        
        
        
        
        ...

    def _plot_basic_setting(self,
                            fig: plt.Figure,
                            ax: plt.Axes):
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

        ax.set_title(
            f"Transit of {self.row['Planet.Name']} on night: {Night.strftime('%Y%m%d')} at {self.Location.name} with precision on $T_c$ of: {self.Uncertainty.to(u.min):.2f}; Quality: {self.quality}", fontsize=14)

        ax.axhline(
            np.arcsin(1/self.Airmass_limit) * 180/np.pi,
            ls='--', lw=2,
            color='darkred'
        )

        ax.grid(color='black', linestyle='dashed', linewidth=0.3, alpha=1.0)
        xfmt = md.DateFormatter('%d/%m  %H')
        ax.xaxis.set_major_formatter(xfmt)
        fig.autofmt_xdate()
        hours = md.HourLocator()
        ax.xaxis.set_major_locator(hours)

        ax.set_ylim(0, 90)
        ax2 = ax.twinx()
        ax.set_ylabel('Altitude [deg]', color='k')
        ax2.set_ylabel('Airmass', color='k')
        ax2.tick_params('y', colors='k')
        ax2.set_ylim([0, 90])
        ax2.set_yticklabels(['', 5.76, 2.92, 2.00, 1.56,
                            1.31, 1.15, 1.06, 1.02, 1.00])

        ax.set_xlim(self.TimeArray.Sunset,
                    self.TimeArray.Sunrise)

    def _plot_twilight(self,
                       ax: plt.Axes):
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
        ax.fill_between(self.TimeArray.time_array.to_value(format='datetime'),
                        0, 90,
                        self.AltitudeArray.Sun.alt.value < 0,
                        color='blue',
                        alpha=0.3
                        )
        ax.fill_between(self.TimeArray.time_array.to_value(format='datetime'),
                        0, 90,
                        self.AltitudeArray.Sun.alt.value < -18,
                        color='blue',
                        alpha=0.5
                        )
        return

    def _plot_Moon(self,
                   ax: plt.Axes):
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
        ax.plot(self.TimeArray.time_array.to_value(format='datetime'),
                self.AltitudeArray.Moon.alt.value,
                ls='-',
                label='Moon',
                linewidth=2,
                color='yellow',
                alpha=0.7
                )
        return

    def _plot_target(self,
                     ax: plt.Axes):
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
        ax.plot(self.TimeArray.time_array.to_value(format='datetime'),
                self.AltitudeArray.target.alt.value,
                ls='-',
                label='Target',
                linewidth=4,
                color='red',
                )
        ax.plot(self.TimeArray.time_array[self.Indices.Transit].to_value(format='datetime'),
                self.AltitudeArray.target[self.Indices.Transit].alt.value,
                ls='-',
                label='Target',
                linewidth=10,
                color='red',
                )

        ax.scatter(
            self.TimeArray.time_array[self.Indices.Transit[0][0]].to_value(
                format='datetime'),
            self.AltitudeArray.target[self.Indices.Transit[0][0]].alt.value,
            color='white',
            zorder=999,
            s=200,
        )
        ax.scatter(
            self.TimeArray.time_array[self.Indices.Transit[-1]
                                      [-1]].to_value(format='datetime'),
            self.AltitudeArray.target[self.Indices.Transit[-1][-1]].alt.value,
            color='white',
            zorder=999,
            s=200,
        )

        return

    def _plot_observation(self,
                          ax: plt.Axes):

        # CLEANME
        self.twilight_observation_time = [
                    self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(
                        format='datetime'),
                    self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(
                        format='datetime')
                ]
                
        if len(self.ObservationsFull.complete) != 0:
            ax.fill_between(self.TimeArray.time_array[self.ObservationsFull.complete].to_value(format='datetime'),
                0, 10,
                color='lime',
                alpha=1,
            )

            ax.text((self.TimeArray.time_array[self.ObservationsFull.complete[0]]+30*u.min).to_value(format='datetime'),
                    5,
                    'Observations (no twilight): ' +
                    self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(
                        format='datetime').strftime('%H:%M') + ' - ' +
                    self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(
                        format='datetime').strftime('%H:%M'),
                    fontsize=10,
                    )
            
            self.full_observation_time = [
                self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(
                    format='datetime'),
                self.TimeArray.time_array[self.ObservationsFull.complete[-1]
                                          ].to_value(format='datetime')
            ]
            self.twilight_full_is_same = (self.full_observation_time == self.twilight_observation_time)
        else:
            self.twilight_full_is_same = False
            
        
        if (((self.FlagsFull.baseline_ingress) or  # Ingress baseline missing
            (self.FlagsFull.baseline_egress) or  # Egress baseline missing
                (not (self.FlagsWindow.visible) and
                 self.FlagsWindow.visible_twilight)) and # Transit into twilight
            not(self.twilight_full_is_same)): # Twilight == Full observations
            
            ax.fill_between(self.TimeArray.time_array[self.ObservationsTwilight.complete].to_value(format='datetime'),
                10, 20,
                color='orange',
                alpha=1,
            )

            ax.text((self.TimeArray.time_array[self.ObservationsTwilight.complete[0]]+30*u.min).to_value(
                format='datetime'),
                15,
                'Observation (with twilight): ' +
                self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(
                format='datetime').strftime('%H:%M') + ' - ' +
                self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(
                format='datetime').strftime('%H:%M'),
                fontsize=10,
            )

            ax.axvline(self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(format='datetime'),
                       ls='--',
                       lw=2,
                       color='orange'
                       )
            ax.axvline(self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(format='datetime'),
                       ls='--',
                       lw=2,
                       color='orange')


        ax.axvline(self.TimeArray.T1.to_value(format='datetime'),
                   ls='--', lw=2, color='black')
        ax.axvline(self.TimeArray.T4.to_value(format='datetime'),
                   ls='--', lw=2, color='black')
        if len(self.ObservationsFull.complete) != 0:
            ax.axvline(self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(
                format='datetime'), ls='--', lw=2, color='lime')
            ax.axvline(self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(
                format='datetime'), ls='--', lw=2, color='lime')

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
        from matplotlib import rc

        if self.AltitudeArray.MoonSeparation.value > 45:
            color_moon = 'white'
        elif self.AltitudeArray.MoonSeparation.value > 30:
            color_moon = 'orange'
        else:
            color_moon = 'red'

        format_text = {'style': 'italic',
                       'size': 12
                       }

        text = ax.text(self.TimeArray.Sunset, 97,
                       f"RA: {self.row['Position.RightAscension']:.2f} ",
                       **format_text
                       )

        text = ax.annotate(
            f"DEC: {self.row['Position.Declination']:.2f} ",
            xycoords=text,
            xy=(1, 0), verticalalignment="bottom",
            **format_text
        )

        text = ax.annotate(
            f"; V = {self.row['Magnitude.V']:.3f} ",
            xycoords=text,
            xy=(1, 0), verticalalignment="bottom",
            **format_text
        )
        text = ax.annotate(
            f"; Moon = {self.AltitudeArray.MoonSeparation.value:.0f} deg ",
            xycoords=text,
            xy=(1, 0), verticalalignment="bottom",
            bbox = dict(facecolor=color_moon, edgecolor=color_moon,
                    boxstyle='round,pad=1'),
            **format_text,
        )
        text = ax.annotate(
            '; $T_{14}$' + f" = {self.row['Planet.TransitDuration']:.3f} h ",
            xycoords=text,
            xy=(1, 0), verticalalignment="bottom",
            **format_text
        )

        text = ax.annotate(
            f"; P = {self.row['Planet.Period']:.5f} d ",
            xycoords=text,
            xy=(1, 0), verticalalignment="bottom",
            **format_text
        )
        
        if 'Flag.TransitTimingVariations' in self.row.keys():
            if self.row['Flag.TransitTimingVariations']:
                text = ax.annotate(
                    f"; TTV detected",
                    xycoords=text,
                    xy=(1, 0), verticalalignment="bottom",
                    bbox = dict(facecolor='red', edgecolor='red',
                    boxstyle='round,pad=1'),
                    **format_text
                )
        
        sunset_text = ax.text(self.TimeArray.Sunset, -10,
                              f"ingress= {self.TimeArray.T1.to_value('datetime').strftime('%H:%M')}",
                              style='italic',
                              bbox=dict(
                                  facecolor='none',
                                  edgecolor='black',
                                  pad=7,
                    boxstyle='round,pad=1'
                              )
                              )

        sunset_text = ax.annotate(
            (f"Observations (no twilight): {self.TimeArray.time_array[self.ObservationsFull.complete[0]].to_value(format='datetime').strftime('%H:%M')} " +
             f"- {self.TimeArray.time_array[self.ObservationsFull.complete[-1]].to_value(format='datetime').strftime('%H:%M')}"),
            xycoords=sunset_text,
            xy=(1.4, 0), verticalalignment="bottom",
            bbox=dict(
                facecolor='green',
                edgecolor='black',
                boxstyle='round, pad=1'
            )
        )

        if not(self.full_observation_time == self.twilight_observation_time):
            sunset_text = ax.annotate(
                (f"Observation (with twilight): {self.TimeArray.time_array[self.ObservationsTwilight.complete[0]].to_value(format='datetime').strftime('%H:%M')} " +
                f"- {self.TimeArray.time_array[self.ObservationsTwilight.complete[-1]].to_value(format='datetime').strftime('%H:%M')}"),
                xycoords=sunset_text,
                xy=(1.4, 0), verticalalignment="bottom",
                bbox=dict(
                    facecolor='orange',
                    edgecolor='black',
                    boxstyle='round,pad=1'
                ),
            )

        ax.text(self.TimeArray.Sunrise, -10,
                f"egress= {self.TimeArray.T4.to_value('datetime').strftime('%H:%M')}",
                style='italic',
                bbox=dict(
                    facecolor='none',
                    edgecolor='black',
                    pad=7)
                )
        return

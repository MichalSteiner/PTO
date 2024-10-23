#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 13:45:10 2023

@author: chamaeleontis
"""

#%% Importing libraries
import datetime
import pandas as pd
import pyvo as vo
import astropy
import astropy.constants as con
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from dataclasses import dataclass, asdict, field
import pickle
import logging
import numpy as np
import os
import seaborn as sns
from colorlog import ColoredFormatter


from ..utilities import logger_default
logger = logging.getLogger(__name__)
logger = logger_default(logger) 

#%% Convert_Earth_to_Jupiter
def rad_earth2jupiter(x): # Convert radius of Earth to Jupiter
    return x * con.R_earth.to(u.R_jup).value

def rad_jupiter2earth(x): # Convert radius Jupiter to Earth
    return x * con.R_jup.to(u.R_earth).value

def mass_earth2jupiter(x): # Convert mass Earth to Jupiter
    return x * con.M_earth.to(u.M_jup).value

def mass_jupiter2earth(x): # Convert mass Jupiter to Earth
    return x * con.M_jup.to(u.M_earth).value

#%% NASA_exoplanets_archive
class NASA_Exoplanets_archive:
    """Class holding information extracted from the NASA Composite table."""
    
    def _save(self):
        """
        Save class as a pickle file to predetermined location.

        Returns
        -------
        None.

        """
        os.makedirs('./saved_files/', mode = 0o777, exist_ok = True) 
        with open('./saved_files/NASA_exoplanets_archive.pkl', 'wb') as output_file:
            pickle.dump(self.__dict__, output_file)


    def _load(self):
        """
        Load class from pickle file in predetermined location.

        Returns
        -------
        None.

        """
        with open('./saved_files/NASA_exoplanets_archive.pkl', 'rb') as input_file:
            self.__dict__  =  pickle.load(input_file)
    
    def _load_API_table(self, full = False):
        """
        Load the Planetary Systems (PS) and Planetary System Composite (PSComp) tables.

        Returns
        -------
        None.

        """
        service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP/") # Initialization of TAP service
        # Searching relevant parameters
        self.nasa_table_composite = pd.DataFrame(service.search("SELECT * FROM pscomppars"))
        if full:
            self.nasa_table_all = pd.DataFrame(service.search("SELECT * FROM ps"))
        else:
            self.nasa_table_all = pd.DataFrame(service.search("SELECT pl_name, pl_tranmid,pl_tranmiderr1, pl_tranmiderr2, pl_orbper, pl_orbpererr1, pl_orbpererr2, pl_trandur, pl_trandurerr1, pl_trandurerr2, pl_refname, dec, sy_vmag, ra FROM ps"))
        self.time = datetime.datetime.now()
        self._get_most_precise()
        self._add_Earth_and_Jupiter_units()
        self._calculate_insolation_flux()
        self._calculate_scale_height()
        self._calculate_transit_length()
       
        
        
        self.filtered_table_composite = self.nasa_table_composite
        self.filtered_table_all = self.nasa_table_all
        self._save()
        
    
    
    def __init__(self,
                 force_reload: bool = False,
                 use_PS: bool = True,
                 full_PS: bool = False):
        """
        Initialize the class. This will try to load the NASA table already saved in predetermined location. If it fails, it will extract the table from NASA archive through TAP service. In case the table is loaded, but the loading was done more than week ago, the table will be reloaded automatically.
    
        Parameters
        ----------
        force_reload : bool, optional
            If turned to True it forces reload of the table. The default is False.
    
        Returns
        -------
        None.
    
        """
        try: # Test loading the function
            logger.info('Trying to load NASA table')
            self._load()
            if (datetime.datetime.now() - self.time).days > 7 or force_reload: # If not < week data or requested, reload anyway
                logger.info('Reloading NASA table [new week] or forced reload. This can take a long time if using the full (PS) table (around 10 min)')
                self._load_API_table(full_PS)
        except: # Loading failed, reload the archive - should happen only the first time running the code
            logger.warning('Weird, no table found. This should happen only for the first time running the code. This also takes a long time if loading both Composite (PSComp) and Full (PS) tables (around 10 min)')
            self._load_API_table(full_PS)
        
        
        
        if not(use_PS): # Not the most effective way of solving
            self.nasa_table_all = self.nasa_table_all[0:0]
            self.filtered_table_all = self.filtered_table_all[0:0]
        return

    
    def _add_Earth_and_Jupiter_units(self):
        """
        Add planetary radii in units of Earth/ Jupiter if only one version exists.

        Returns
        -------
        None.

        """
# =============================================================================
#         Earth to Jupiter radius
# =============================================================================
        condition_composite = (self.nasa_table_composite['pl_rade'].notna() & # Planetary radius [Earth units]
                               self.nasa_table_composite['pl_radj'].isna() # Planetary radius [Jupiter units]
                                )
        indices_composite = self.nasa_table_composite[condition_composite].index
        
        self.nasa_table_composite.loc[indices_composite, 'pl_radj'] = (
            self.nasa_table_composite.loc[indices_composite, 'pl_rade'] * con.R_jup / con.R_earth
            )
        
# =============================================================================
#         Jupiter to Earth radius
# =============================================================================
        condition_composite = (self.nasa_table_composite['pl_radj'].notna() & # Planetary radius [Jupiter units]
                               self.nasa_table_composite['pl_rade'].isna() # Planetary radius [Earth units]
                                )
        indices_composite = self.nasa_table_composite[condition_composite].index
        
        self.nasa_table_composite.loc[indices_composite, 'pl_rade'] = (
            self.nasa_table_composite.loc[indices_composite, 'pl_radj'] * con.R_earth / con.R_jup
            )
        return

# =============================================================================
#     Calculation of additional parameters
# =============================================================================
    def _calculate_transit_length(self):
        r"""
        Calculates the transit length.
        
        Equations used:
            T_{14} = \frac{P}{\pi} * \arcsin(\frac{\frac{R_s}{a} * \sqrt{(1 + \frac{R_p}{R_s})^2 - b^2)}}{\sin(i)})
            T_{14} = \frac{P}{\pi} * \arcsin{(\frac{\frac{R_s}{a} * \sqrt{1+ (\frac{R_p}{R_s})^2} - b^2}{\sin{i}})}
        
        In case of eccentricity, correction factor is used:
            \frac{\sqrt{1-e^2}}{(1+e)\sin(\omega)} 
            \frac{\sqrt{1-e^2}}{(1+e)} 

        Returns
        -------
        None.

        """
        logger.critical('Implement uncertainties when calculating transit length')
        condition_composite = (
            self.nasa_table_composite['pl_trandur'].isna() & # Transit duration
            self.nasa_table_composite['pl_orbper'].notna() & # Orbital period
            self.nasa_table_composite['pl_imppar'].notna() & # Impact parameter
            self.nasa_table_composite['pl_orbincl'].notna() & # Orbital inclination
            self.nasa_table_composite['pl_ratdor'].notna() & # Ratio of semimajor axis to stellar radius (a/Rs)
            self.nasa_table_composite['pl_ratror'].notna() # Ratio of planet to stellar ratio (Rp/Rs)
            )
        
        indices_composite = self.nasa_table_composite[condition_composite].index
        indices_composite_combined = indices_composite
        
# =============================================================================
#         $$
#         T_{14} = \frac{P}{\pi} * \arcsin(\frac{\frac{R_s}{a} * \sqrt{(1 + \frac{R_p}{R_s})^2 - b^2)}}{\sin(i)})
#         $$
# =============================================================================
        self.nasa_table_composite.loc[indices_composite, 'pl_trandur'] = (
            (
            (self.nasa_table_composite.loc[indices_composite, 'pl_orbper'] * 24 / np.pi)*
            np.arcsin(
                    (1./
                     self.nasa_table_composite.loc[indices_composite, 'pl_ratdor']
                     ) *
                    np.sqrt((1. +
                             self.nasa_table_composite.loc[indices_composite, 'pl_ratror'])**2. -
                            self.nasa_table_composite.loc[indices_composite, 'pl_imppar']**2.
                            ) /
                    np.sin((self.nasa_table_composite.loc[indices_composite, 'pl_orbincl']*np.pi/180.))
                    )
            )
            )
        
        condition_composite = (
            self.nasa_table_composite['pl_trandur'].isna() & # Transit duration
            self.nasa_table_composite['pl_orbper'].notna() & # Orbital period
            self.nasa_table_composite['pl_imppar'].notna() & # Impact parameter
            (self.nasa_table_composite['pl_imppar'] >= 0) & # Impact parameter
            (self.nasa_table_composite['pl_imppar'] <= 1) & # Impact parameter
            self.nasa_table_composite['pl_orbincl'].notna() & # Orbital inclination
            self.nasa_table_composite['pl_orbsmax'].notna() & # Semi-major axis
            
            self.nasa_table_composite['pl_rade'].notna() & # Planetary radius [Earth units]
            self.nasa_table_composite['st_rad'].notna() # Stellar radius [Solar units]
            )
        
        
        indices_composite = self.nasa_table_composite[condition_composite].index

        
# =============================================================================
#         $$
#         T_{14} = \frac{P}{\pi} * \arcsin{(\frac{\frac{R_s}{a} * \sqrt{1+ (\frac{R_p}{R_s})^2} - b^2}{\sin{i}})}
#         $$
# =============================================================================
        self.nasa_table_composite.loc[indices_composite, 'pl_trandur'] = (
            (self.nasa_table_composite.loc[indices_composite, 'pl_orbper'] * 24 / np.pi) *
            np.arcsin(
                        (
                    self.nasa_table_composite.loc[indices_composite, 'st_rad']  /
                    self.nasa_table_composite.loc[indices_composite, 'pl_orbsmax'] *
                    con.R_sun / con.au 
                            ) *
                    np.sqrt((1. +
                            (
                            (self.nasa_table_composite.loc[indices_composite, 'pl_rade'])/
                            self.nasa_table_composite.loc[indices_composite, 'st_rad'] *
                            con.R_earth / con.R_sun)**2 -
                            self.nasa_table_composite.loc[indices_composite, 'pl_imppar']**2.
                            ) /
                    np.sin((self.nasa_table_composite.loc[indices_composite, 'pl_orbincl']*np.pi/180.))
                            )
                    )
            )
        
        indices_composite_combined = indices_composite_combined.union(indices_composite)
# =============================================================================
#         Eccentricity correction
# =============================================================================
        condition_composite = (
            (self.nasa_table_composite['pl_orbeccen'] > 0.) & # Eccentricity
            self.nasa_table_composite['pl_orblper'].notna() # Argument of periastron
            )
        
        indices_composite = self.nasa_table_composite[condition_composite].index
        indices_composite = indices_composite_combined.intersection(indices_composite)
        
# =============================================================================
#         $$
#         T_{14} = T_{14} * \frac{\sqrt{1-e^2}}{(1+e)\sin(\omega)} 
#         $$
# =============================================================================
        self.nasa_table_composite.loc[indices_composite, 'pl_trandur'] *= np.sqrt(
            1. - self.nasa_table_composite.loc[indices_composite,'pl_orbeccen']**2.
            ) / (
                1. + (self.nasa_table_composite.loc[indices_composite, 'pl_orbeccen'] * 
                np.sin(self.nasa_table_composite.loc[indices_composite, 'pl_orblper'] * 
                       np.pi / 180.))
                )
        
        
        condition_composite = (
            (self.nasa_table_composite['pl_orbeccen'] > 0.) & # Eccentricity
            self.nasa_table_composite['pl_orblper'].isna() # Argument of periastron
            )
        
        indices_composite = self.nasa_table_composite[condition_composite].index
        indices_composite = indices_composite_combined.intersection(indices_composite)
        
# =============================================================================
#         $$
#         T_{14} = T_{14} * \frac{\sqrt{1-e^2}}{(1+e)} 
#         $$
# =============================================================================
        self.nasa_table_composite.loc[indices_composite, 'pl_trandur'] *= np.sqrt(
            1. - self.nasa_table_composite['pl_orbeccen']**2.
            ) / (1. + self.nasa_table_composite.loc[indices_composite, 'pl_orbeccen'])
        
        return
    
    
    def _calculate_insolation_flux(self):
        """
        Calculate insolation flux for planets with sufficient parameters.

        Returns
        -------
        None.

        """
        condition_composite = (self.nasa_table_composite['pl_insol'].isna() & 
                               self.nasa_table_composite['st_lum'].notna() &
                               self.nasa_table_composite['pl_orbsmax'].notna()
                               )
        
        indices_composite = self.nasa_table_composite[condition_composite].index
        
        self.nasa_table_composite.loc[indices_composite, 'pl_insol'] = (
            10**self.nasa_table_composite.loc[indices_composite, 'st_lum'] / 
            self.nasa_table_composite.loc[indices_composite,'pl_orbsmax']**2
            )
        
        return
    
    def _calculate_scale_height(self,
                               molecular_mass:float = 2.3):
        """
        Calculate atmospheric scale heights for all planets with sufficient parameters, assuming molecular mass of 2.3kg/mol as default.

        Parameters
        ----------
        mol_mass : float, optional
            Molecular mass in kilograms per mol [kg/mol]. The default is 2.3.

        Returns
        -------
        None.

        """
        try:
            self.nasa_table_composite['H']
        except:
            self.nasa_table_composite['H'] = [np.nan]*len(self.nasa_table_composite)
        
        condition_composite = (
                            self.nasa_table_composite['H'].isna() &
                            self.nasa_table_composite['pl_bmassj'].notna() &
                            self.nasa_table_composite['pl_radj'].notna() &
                            self.nasa_table_composite['pl_eqt'].notna()
                            )
        
        indices_composite = self.nasa_table_composite[condition_composite].index
        
        logger.critical('Atmospheric scale height calculation is not debugged yet!')
        self.nasa_table_composite.loc[indices_composite, 'H']  = (
            con.N_A *
            self.nasa_table_composite.loc[indices_composite, 'pl_eqt'] * u.K  *
            ((self.nasa_table_composite.loc[indices_composite, 'pl_radj'] * (1*u.R_jup).to(u.m))**2) *
            con.k_B /
            (
            (molecular_mass*u.kg/u.mol) *
            con.G *
            self.nasa_table_composite.loc[indices_composite, 'pl_bmassj'] * (1*u.M_jup).to(u.kg))
            )
        
        return
    
    
    def reset_filter(self):
        """
        Reset the filtered_table atributes to full table.

        Returns
        -------
        None.

        """
        self.filtered_table_composite = self.nasa_table_composite
        self.filtered_table_all = self.nasa_table_all
        return
    
    def _bind_filter_all(self):
        """
        Bind the filter used by composite table to the all table.

        Returns
        -------
        None.

        """
        condition = self.filtered_table_all['pl_name'].isin(
            self.filtered_table_composite['pl_name'].values
            )
        self.filtered_table_all = self.filtered_table_all[condition]
        return
    
    def _bind_filter_composite(self):
        """
        Bind the filter used by composite table to the all table.

        Returns
        -------
        None.

        """
        condition = self.filtered_table_composite['pl_name'].isin(
            self.filtered_table_all['pl_name'].values
            )
        self.filtered_table_composite = self.filtered_table_composite[condition]
        return
    
    def filter_table_composite(self,
                               condition: pd.core.series.Series) -> pd.DataFrame:
        """
        Filter the composite table by pre-defined condition indices.

        Parameters
        ----------
        condition : pd.core.series.Series
            Series with boolean values that filter the composite table.

        Returns
        -------
        filtered_table_composite : pd.DataFrame
            Filtered table.

        """
        self.filtered_table_composite = self.nasa_table_composite[condition]
        self._bind_filter_all()
        return self.filtered_table_composite
    
    def filter_table_all(self, condition):
        """
        Filter the complete table by pre-defined condition indices.

        Parameters
        ----------
        condition : pd.core.series.Series
            Series with boolean values that filter the composite table.

        Returns
        -------
        filtered_table_all : pd.DataFrame
            Filtered table.

        """
        self.filtered_table_all = self.nasa_table_all[condition]
        self._bind_filter_composite()
        return self.filtered_table_all
    
    def _check_if_fig_and_ax_exist(self,
                                   fig: matplotlib.figure.Figure | None,
                                   ax: plt.Axes | None) -> tuple[matplotlib.figure.Figure,
                                                                 plt.Axes]:
        """
        Checks whether figure and artist is created, and will create new one if not.

        Parameters
        ----------
        fig : matplotlib.figure.Figure | None
            Variable of the figure. If None, the figure is created with single artist
        ax : plt.Axes | None
            Variable of the artist. If None, a single artist is created

        Returns
        -------
        fig : matplotlib.figure.Figure
            Resulting figure, which is unchanged if passed, and created if None is passed instead.
        ax : plt.Axes
            Resulting artist, which is unchanged if passed, and created if None is passed instead.

        """
        
        if fig is None:
            fig,ax = plt.subplots(1)
            return fig, ax 
        else:
            return fig, ax
    
    def _default_plot_insolation_flux_diagram_settings(self,
                                                       ax):
        """Set the ax scales to log scale."""
        ax.set_xscale('log')
        ax.set_yscale('log')
        return
        
    
    def plot_radius_insolation_flux_diagram(self,
                                            fig: matplotlib.figure.Figure | None = None,
                                            ax: plt.Axes | None =None,
                                            ) -> tuple[matplotlib.figure.Figure, plt.Axes]:
        logger.critical("TODO: Rewrite this function to be more general")
        fig, ax = self._check_if_fig_and_ax_exist(fig,ax)
        
        self._default_plot_insolation_flux_diagram_settings(ax)
        
        condition_composite = (self.nasa_table_composite['pl_insol'].notna() &
                               self.nasa_table_composite['pl_rade'].notna() &
                               (self.nasa_table_composite['pl_insol'] != 0) &
                               (self.nasa_table_composite['pl_rade'] != 0)
                               )
        
        indices_composite = self.nasa_table_composite[condition_composite].index
        
        ax.invert_xaxis()
        

        sns.kdeplot(x = self.nasa_table_composite.loc[indices_composite, 'pl_insol'],
                    y = self.nasa_table_composite.loc[indices_composite, 'pl_rade'],
                    # cmap=cmap,
                    fill=True,
                    label='_nolegend_',
                    bw_adjust=0.6,
                    log_scale =[True, True],
                    levels = [0,.0001,.001,.01,.1,.3,.5,.8,1],
                    ax=ax
                    )

        ax.scatter(
                self.nasa_table_composite.loc[indices_composite, 'pl_insol'],
                self.nasa_table_composite.loc[indices_composite, 'pl_rade'],
                # color = 'black',
                # edgecolors = 'black',
                label='_nolegend_',
                alpha = 0.5,
                s = 20,
            )
        
        
        ax2 = ax.secondary_yaxis('right', functions=(rad_earth2jupiter, rad_jupiter2earth))
        ax.set_title('Radius vs Insolation flux for known exoplanets')

        ax.set_xlabel('Insolation flux [$S_{\oplus}$]')
        ax2.set_ylabel('Planetary radius $[R_{Jup}]$')
        ax.set_ylabel('Planetary radius $[R_{\oplus}]$')
        
        ax.set_xlim(10E4,10E0)
        ax.set_ylim(10E-1,25E0)
        ax.get_figure().gca().set_title("")

        # if ykey == 'pl_bmasse':
        #     ax2 = ax.secondary_yaxis('right', functions=(mass_earth2jupiter, mass_jupiter2earth))
        #     ax2.set_ylabel('Planetary mass $[M_{Jup}]$', fontsize = 36)
        return fig,ax
        
    
    def plot_mass_insolation_flux_diagram(self,
                                          fig=None,
                                          ax=None):
        fig, ax = self._check_if_fig_and_ax_exist(fig,ax)
        
        return
    
    def create_table_observation(self):
        return
    
    def add_custom_ephemeris(self,
                             Transit_midpoint: float,
                             Transit_midpoint_lower_error: float,
                             Transit_midpoint_upper_error: float,
                             Period: float,
                             Period_lower_error: float,
                             Period_upper_error: float,
                             Transit_length: float,
                             Transit_length_lower_error: float,
                             Transit_length_upper_error: float,
                             Right_Ascension: float,
                             Declination: float,
                             V_magnitude: float,
                             Reference: str,
                             Planet_name: str
                             ):
        """
        Adds a custom ephemeris to the table.

        Parameters
        ----------
        Transit_midpoint : float
            Transit midpoint of the system.
        Transit_midpoint_lower_error : float
            Transit midpoint minus uncertainty.
        Transit_midpoint_upper_error : float
            Transit midpoint plus uncertainty.
        Period : float
            Period of the system.
        Period_lower_error : float
            Period minus uncertainty.
        Period_upper_error : float
            Period plus uncertainty.
        Transit_length : float
            Transit length of the system.
        Transit_length_lower_error : float
            Transit length minus uncertainty.
        Transit_length_upper_error : float
            Transit length plus uncertainty.
        Right_Ascension : float
            Right ascension of the system.
        Declination : float
            Declination of the system.
        V_magnitude : float
            V magnitude of the system.
        Reference : str
            Reference to use for the ephemeris.
        Planet_name : str
            Planet name to use for the ephemeris.

        Returns
        -------
        None.

        """
        
        parameters = {key: '' for key in self.nasa_table_all.columns}
        
        parameters.update({
            'pl_name': Planet_name,
            'pl_orbper': Period, 
            'pl_orbpererr1': Period_upper_error, 
            'pl_orbpererr2': Period_lower_error, 
            'pl_trandur': Transit_length,
            'pl_trandurerr1': Transit_length_upper_error,
            'pl_trandurerr2': Transit_length_lower_error,
            'pl_tranmid': Transit_midpoint,
            'pl_tranmiderr1': Transit_midpoint_upper_error,
            'pl_tranmiderr2': Transit_midpoint_lower_error,
            'ra': Right_Ascension,
            'dec': Declination,
            'sy_vmag': V_magnitude,
            'pl_refname': Reference,
            })
        
        pieces = [self.nasa_table_all, pd.DataFrame(parameters,  index=[0])]
        
        self.nasa_table_all = pd.concat(pieces)
        try:
            self.nasa_table_all = self.nasa_table_all.reset_index()
        except:
            pass
        
        pieces = [self.filtered_table_all, pd.DataFrame(parameters,  index=[0])]
        self.filtered_table_all = pd.concat(pieces)
        try:
            self.filtered_table_all = self.filtered_table_all.reset_index()
        except:
            pass
        
        return

    def _get_most_precise_row(self,
                              planet_data: pd.DataFrame,
                              composite_data: pd.DataFrame):
        for key in [
            'pl_orbper',
            'pl_tranmid',
            'pl_trandur'
            ]:
            condition = (composite_data[key].isna() |
                         composite_data[key+'err1'].isna() |
                         composite_data[key+'err2'].isna())
            
            indices_composite = composite_data[condition].index
            
            average_error = (planet_data[key+'err1'] + abs(planet_data[key+'err2']))/2
            if average_error.idxmin() is np.nan:
                continue
            composite_data.loc[indices_composite, key] = planet_data.loc[average_error.idxmin()][key]
            composite_data.loc[indices_composite, key+'err2'] = planet_data.loc[average_error.idxmin()][key+'err1']
            composite_data.loc[indices_composite, key+'err2'] = planet_data.loc[average_error.idxmin()][key+'err2']
        
        
        return composite_data
    
    
    def _get_most_precise(self):
        
        condition = (
            (self.nasa_table_composite['pl_tranmid'].isna() |
            self.nasa_table_composite['pl_tranmiderr1'].isna() |
            self.nasa_table_composite['pl_tranmiderr2'].isna() |
            self.nasa_table_composite['pl_orbper'].isna() |
            self.nasa_table_composite['pl_orbpererr1'].isna() |
            self.nasa_table_composite['pl_orbpererr2'].isna() |
            self.nasa_table_composite['pl_trandur'].isna() |
            self.nasa_table_composite['pl_trandurerr1'].isna() |
            self.nasa_table_composite['pl_trandurerr2'].isna()) &
            (self.nasa_table_composite['tran_flag'] == 1) # Transiting system
            )
        indices = self.nasa_table_composite[condition].index
        
        for value in self.nasa_table_composite.loc[indices, 'pl_name']:
            updated_row = self._get_most_precise_row(planet_data= self.nasa_table_all[self.nasa_table_all['pl_name'] == value],
                                                     composite_data= self.nasa_table_composite[self.nasa_table_composite['pl_name'] == value])
            
            condition = (self.nasa_table_composite['pl_name'] == value)

            indices_composite = self.nasa_table_composite[condition].index
            self.nasa_table_composite.loc[indices_composite] = updated_row

            # self.nasa_table_composite[self.nasa_table_composite['pl_name'] == value] = updated_row
        
        return

def load_csv_ephemeris(Exoplanets: NASA_Exoplanets_archive,
                       table: str,
                       reference_name: str):
    logger.info(f'Loading ephemeris from table located: {table}')
    
    read_table = pd.read_csv(table, sep= ';')
    
    for ind, row in read_table.iterrows():
        
        Exoplanets.add_custom_ephemeris(
            Transit_midpoint = float(row['T0']),
            Transit_midpoint_lower_error = float(row['T0_err']),
            Transit_midpoint_upper_error = float(row['T0_err']),
            Period = float(row['P']),
            Period_lower_error = float(row['P_err']),
            Period_upper_error = float(row['P_err']),
            Transit_length = float(row['T14']),
            Transit_length_lower_error = 0,
            Transit_length_upper_error = 0,
            Right_Ascension = float(row['ra']),
            Declination = float(row['dec']),
            V_magnitude = float(row['mag_V']),
            Reference = reference_name,
            Planet_name = row['pl_name'],
            )
    return



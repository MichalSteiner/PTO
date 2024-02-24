#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 10:37:18 2023

@author: chamaeleontis
"""

#%% Import libraries
from enum import Enum
import astropy.coordinates as coord
import astropy.units as u
#%%
class Observatories(Enum):
    """
    Enum class object holding all the observatories we wish to use.
    
    In case observatory is not found, please add a new parameter to the class, with latitude, longitude and height.
    
    Example of how to setup custom Earth location
        Paranal = coord.EarthLocation(
            lat = coord.Angle('-24:37:38.0', unit=u.degree),
            lon = coord.Angle('-70:24:17.0', unit=u.degree),
            height = 2635 * u.m,
            )
        
        It is possible your observatory is in the astropy.coordinate list of known sites
        You can extract this list by  coord.EarthLocation.get_site_names()
    """
    
    Paranal = coord.EarthLocation.of_site('paranal')
    LaSilla = coord.EarthLocation.of_site('lasilla')
    LaPalma = coord.EarthLocation.of_site('lapalma')



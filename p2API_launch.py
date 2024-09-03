#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 10:49:08 2023

@author: chamaeleontis
"""

#%% Import API
import pprint
import p2api
import pandas as pd
import PTO.transits.NASA_exo
import astropy.units as u
from astroquery.simbad import Simbad
from dataclasses import dataclass
import os
import astropy.time as astime
import logging

from PTO.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 

p = pprint.PrettyPrinter(indent=4)

#%% NightList class
@dataclass
class NightList:
    """Wrapper class holding all the nights information"""
    
    Exoplanets_archive: PTO.transits.NASA_exo.NASA_Exoplanets_archive
    
    
    def _generate_obs_ESPRESSO(self,
                               target_information: pd.DataFrame
                               ):
        """
        Generate default observation with ESPRESSO for the purpose of Rossiter-McLaughlin effect and transmission spectroscopy observations.

        Parameters
        ----------
        observation : eso.P2.ob
            Observation template created by the p2 API.
        target : str
            Target name.

        Returns
        -------
        None.

        """
        target = target_information['Name'].values[0][:-2]
        
        try:
            sdata = simb.query_object(target)
            self.observing_block['target']['name'] = target
            rasp = sdata['RA_hms_A_ICRS_J2000_2000'][0].replace(' ', ':').split('.')
            self.observing_block['target']['ra'] = f'{rasp[0]}.{rasp[1][:3]}'
            decsp = sdata['DEC_hms_D_ICRS_J2000_2000'][0].replace(' ', ':').split('.')
            self.observing_block['target']['dec'] = f'{decsp[0]}.{decsp[1][:3]}'
            self.observing_block['target']['properMotionRa'] = sdata['PMRA'][0] / 1000
            self.observing_block['target']['properMotionDec'] = sdata['PMDEC'][0] / 1000
        except:
            logger.warning("Target", target_information['Name'], 'not found on SIMBAD')
        
        self.observing_block['instrument'] = 'ESPRESSO'
        
        self.observing_block['constraints']['name'] = target_information['Name'].values[0][:-2]
        self.observing_block['constraints']['airmass'] = 2.2
        self.observing_block['constraints']['seeing'] = 2
        self.observing_block['constraints']['fli'] = 1
        self.observing_block['constraints']['skyTransparency'] = 'Variable, thin cirrus'
        self.observing_block['constraints']['moonDistance'] = 30
        self.observing_block['constraints']['twilight'] = 0
        
        self.observing_block['obsDescription']['userComments'] = 'Transit of: ' + self.observing_block['name']
        self.observing_block['obsDescription']['name'] = self.observing_block['name']
        
        return None
        
    def generate_observation_list(self,
                                  user: str,
                                  password: str,
                                  runContainerId: int,
                                  list_of_planets: list,
                                  directory_windows: str,
                                  SM_mode: bool,
                                  quality: int
                                  ):
        
        self._connect_P2API_Paranal(user, password)
        
        for planet, SM in list_of_planets:
            
            if SM != SM_mode:
                continue
            
            csv_pathname = f'{directory_windows}/{planet.replace(" ","")}/'
            for file in os.listdir(csv_pathname):
                if file.startswith('best_') and file.endswith('.csv'):
                    csv_pathname += file
            
            
            target_info = pd.read_csv(csv_pathname, sep=';', names= ['Name', 'Night', 'Quality', 'Start', 'End','SM'])
            
            
            
            self.observing_block, self.observing_block_ID = None, None
            self.observing_block_version = None
            
            logger.info('Currently working on observations of:', target_info['Name'].values[0])
            folder, folderVersion = self.api.createFolder(runContainerId, target_info['Name'].values[0])
            folderId = folder['containerId']

# =============================================================================
#             Create observing block for given planet
# =============================================================================
            observing_block, observing_block_version = self.api.createOB(folderId, target_info['Name'].values[0])
            self.observing_block = observing_block
            self.observing_block_ID = observing_block['obId']
            self.observing_block_version = observing_block_version
# =============================================================================
#             Generate observation for given planet
# =============================================================================
            self._generate_obs_ESPRESSO(target_info)
            self._create_template(self.observing_block_ID, target_info, planet)
            self.observing_block, self.observing_block_version = self.api.saveOB(
                self.observing_block,
                self.observing_block_version)
            
            if SM == True:
                self._define_SM_time_constraint(self.observing_block_ID, target_info, quality)
            else:
                self._define_dVM_time_constraint(self.observing_block_ID, target_info, quality)
                
            
            self._verify_OB(self.observing_block, self.observing_block_ID)
            
            p.pprint(self.observing_block)
        
        return
    
    
    def _connect_P2API_Paranal(self,
                               user: str,
                               password: str):
        
        logger.info('About to connect to ESO Paranal API')
        # self.api = p2api.ApiConnection('demo', '52052', 'tutorial')
        self.api = p2api.ApiConnection('production',
                                  user,
                                  password)
        logger.info('Connected to ESO Paranal API')
        self.api.p2fc_url = 'https://www.eso.org/p2fc/'
        
        return
    
    def _verify_OB(self,
                   ob,
                   observing_block_ID
                   ):
        
        response, _ = self.api.verifyOB(observing_block_ID, True)
        if response['observable']:
          logger.info('*** Congratulations. Your OB', observing_block_ID, ob['name'], 'is observable!')
        # else:
        #   logger.info('OB', observing_block_ID, 'is >>not observable<<. See messages below.')
        print(' ', '\n  '.join(response['messages']))

    
    def _create_template(self,
                         observing_block_ID,
                         observation,
                         planet
                         ):
        
# =============================================================================
#         Create acquisition template
# =============================================================================
        self.acquisition_template,\
        self.acquisition_template_Version = self.api.createTemplate(
                                                    observing_block_ID,
                                                    'ESPRESSO_singleHR_acq_obj'
                                                    )
        try:
            star_V_mag = self.Exoplanets_archive.nasa_table_composite[
                self.Exoplanets_archive.nasa_table_composite['pl_name']==planet
                ]['sy_vmag'].values[0]
            stellar_type = self.Exoplanets_archive.nasa_table_composite[
                self.Exoplanets_archive.nasa_table_composite['pl_name']==planet
                ]['st_spectype'].values[0]
        except:
            logger.warning(f'{planet} has not been found in Exoplanet archive, replacing V mag and stellar type with 999 and "NA" respectivelly')
            star_V_mag = 999
            stellar_type = 'NA'
        
        if len(stellar_type) > 2:
            stellar_type= stellar_type[:2]
        
        self.acquisition_template, self.acquisition_template_Version = self.api.setTemplateParams(
                observing_block_ID,
                self.acquisition_template,
                {
                'OCS.OBJ.MV': star_V_mag, # V magnitude of host star
                'OCS.OBJ.RV': 0, # Guess radial velocity
                'OCS.OBJ.SP.TYPE': stellar_type
                },
                self.acquisition_template_Version
                )
        
# =============================================================================
#         Create science image template
# =============================================================================
        self.science_template, \
        self.science_template_Version = self.api.createTemplate(
                                                    observing_block_ID,
                                                    'ESPRESSO_singleHR_obs'
                                                    )
        
        self.science_template, self.science_template_Version = self.api.setTemplateParams(
                observing_block_ID,
                self.science_template,
                {
                'DET1.MODE': '2x1_SLOW',# Binning mode
                'DET1.UIT1': 900, # Exposure time
                'SEQ.CALSOURCEB': 'SKY', # Fiber B
                'SEQ.NO': 100 # Number of exposures
                },
                self.science_template_Version
                )
        
        p.pprint(self.science_template)
        
        self.api.generateFindingChart(observing_block_ID)
        
        return
    
    def _define_dVM_time_constraint(self,
                                    obId,
                                    target_info,
                                    quality):
        
        absTCs, atcVersion = self.api.getAbsoluteTimeConstraints(obId)
        

        
        night_list = []
        
        for ind, row in target_info.iterrows():
            if row['Quality'] > quality:
                continue
            
            start_date, start_time = row['Night'], row['Start']
            end_date, end_time = row['Night'], row['End']
            night_list.append(
                {
                    'from':self._format_time_dVM(start_date, start_time),
                    'to':self._format_time_dVM(end_date, end_time)
                }
            )
            
        
        absTCs, atcVersion = self.api.saveAbsoluteTimeConstraints(obId,
                                                             night_list,
                                                             atcVersion)
    
        return
    
    def _define_SM_time_constraint(self,
                                   obId,
                                   target_info,
                                   quality):
        absTCs, atcVersion = self.api.getAbsoluteTimeConstraints(obId)
        
        night_list = []
        
        for ind, row in target_info.iterrows():
            if row['Quality'] > quality:
                continue
            if row['SM'] == 0:
                continue
            
            start_date, start_time = row['Night'], row['Start']
            
            start_SM, end_SM = self._format_time_SM(start_date, start_time)
            
            night_list.append(
                {
                    'from':start_SM,
                    'to': end_SM
                }
            )
        
        absTCs, atcVersion = self.api.saveAbsoluteTimeConstraints(obId,
                                                             night_list,
                                                             atcVersion)
        
        return
    
    def _format_time_SM(self, date, time):
        date = str(date)
        date = astime.Time(f'{date[:4]}-{date[4:6]}-{date[6:]}')
        date += int(time[:2]) * u.hour
        date += int(time[3:]) *u.minute
        
        if int(time[:2])< 12:
            date += 1*u.day
        date.format = 'isot'
        start_date = date - 15*u.minute
        end_date = date + 15*u.minute
        
        return start_date.value[:-4], end_date.value[:-4]

    def _format_time_dVM(self, date, time):
        date = str(date)
        date = astime.Time(f'{date[:4]}-{date[4:6]}-{date[6:]}')
        date += int(time[:2]) * u.hour
        date += int(time[3:]) *u.minute
        
        if int(time[:2])< 12:
            date += 1*u.day
        date.format = 'isot'
        
        return date.value[:-4]
    
#%% Target information
# =============================================================================
# Load targets informations from the ATREIDES spreadsheet
# =============================================================================
simb = Simbad()
old_fields = list(simb.get_votable_fields())
simb.add_votable_fields('ra(hms;A;ICRS;J2000;2000)',
  'dec(hms;D;ICRS;J2000;2000)', 'pmra', 'pmdec', 'sptype', 'flux(V)',
  'flux(I)', 'flux(J)')
simb.remove_votable_fields(*old_fields)


# TODO: This is not optimal solution. A better handling of user/password is recommended
user = 'Cham'
password = 'MsVbPgozj3Sx2XWc'
runContainerId_SMP114 = 3953548
runContainerId_VMP114 = 3961226
runContainerId_dVMP114 = 3953991
Exoplanets = PTO.transits.NASA_exo.NASA_Exoplanets_archive(
        force_reload = False,
        # use_PS = False,
        full_PS = False
        )

list_of_planets = [
    ['K2-105 b', False],#VM
    # ['K2-27 b', False], #VM
    # ['K2-334 b', False], #dVM
    # ['K2-353 b', False], #dVM
    # ['K2-370 b', True], #SM
    # ['TOI-1231 b', False], #dVM
    # ['TOI-257 b', False], #dVM
    # ['TOI-431 d', True], #SM
    # ['TOI-2000 c', False], #VM
    # ['HD 93963 A c', False], #VM
    # ['K2-100 b', True], # SM
    # ['TOI-1853 b', True] #SM
    ]

directory_windows = '/media/chamaeleontis/Observatory_main/ESO_scheduling/ATREIDES_P114_planning_only/Paranal'
directory_windows = '/media/chamaeleontis/Observatory_main/ESO_scheduling/ATREIDES_P114__only_20240527/Paranal'



ATREIDES_P114 = NightList(
    Exoplanets_archive= Exoplanets,
    )

ATREIDES_P114.generate_observation_list(
    user= user,
    password= password,
    runContainerId= runContainerId_VMP114,
    list_of_planets= list_of_planets,
    directory_windows= directory_windows,
    SM_mode= False,
    quality= 9
    )
# ATREIDES_P114.generate_observation_list(
#     user= user,
#     password= password,
#     runContainerId= runContainerId_SMP114,
#     list_of_planets= list_of_planets,
#     directory_windows= directory_windows,
#     SM_mode= True,
#     quality= 5
#     )


print('Yipeeeeeeeeeeeee')


#%% SIMBAD
# =============================================================================
# Initialization of SIMBAD
# =============================================================================



# #%% List_of_nights
# # =============================================================================
# # Create a structure class "NightList" holding all nights
# # =============================================================================
# Nights_P112 = NightList(
#     information = targets_information
#     )


# for ind, planet in targets_observed.iterrows():
#     if planet['mode'] == 'dVM':
#         continue

#     planet_name = planet['pl_name']
#     nights_observed = planet['nights'].split(' ; ')
#     mode = planet['mode']
    
#     Nights_P112.add_planet(
#         planet_name = planet_name,
#         nights = nights_observed,
#         mode = mode
#         )



#%% Generate OBs
# =============================================================================
# Iterativelly generates obs in the P2 system according to the ATREIDES spreadsheet
# =============================================================================



# Nights_P112.generate_observation_list()


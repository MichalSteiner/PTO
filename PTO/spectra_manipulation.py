# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Fri Sep  1 15:56:11 2023

# @author: chamaeleontis
# """
# import numpy as np
# import specutils as sp
# import scipy as sci
# import astropy
# import pandas as pd
# import astropy.units as u
# from astropy.wcs import WCS
# import astropy.constants as con
# #%% interpolate2commonframe
# def interpolate2commonframe(spectrum,new_spectral_axis):
#     """
#     # TODO
#         Clean up the masking region handling
#     Interpolate spectrum to new spectral axis

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Spectrum to interpolate to new spectral axis
#     new_spectral_axis : sp.Spectrum1D.spectral_axis 
#         Spectral axis to intepolate the spectrum to

#     Returns
#     -------
#     new_spectrum : sp.Spectrum1D
#         New spectrum interpolated to given spectral_axis
        
#     TODO:
#         Check how does this function behaves with masks
#         Don't interpolate for same wavelength grid

#     """
    
#     # Mask handling
#     mask_flux = ~np.isfinite(spectrum.flux) # Ensure nans are not included
#     mask_err = ~np.isfinite(spectrum.uncertainty.array) # Sometimes y_err is NaN while flux isnt? Possible through some divide or np.sqrt(negative)
#     mask = mask_flux + mask_err # Gives zero values in each good pixel (values are False and False)
#     mask = ~mask # Indices of good pixels (both flux and error)
        
#     change_value = np.where(mask[:-1] != mask[1:])[0]
#     mask_region_list = []
#     for ind,value in enumerate(change_value):
#         if ind == len(change_value)-1:
#             break
#         next_value = change_value[ind+1]
#         if mask[value] and ~mask[value+1] and ~mask[next_value] and mask[next_value+1]:
#             mask_region_list.append(sp.SpectralRegion(
#                 np.nanmean([new_spectral_axis[value].value,new_spectral_axis[value+1].value])*new_spectral_axis.unit
#                 ,
#                 np.nanmean([new_spectral_axis[next_value].value,new_spectral_axis[next_value+1].value])*new_spectral_axis.unit
#                 ))
#         pass

#     # Interpolation function for flux - cubic spline with no extrapolation
#     flux_int = sci.interpolate.CubicSpline(spectrum.spectral_axis[mask],
#                                            spectrum.flux[mask],
#                                            extrapolate= False)
#     # Interpolation function for uncertainty - cubic spline with no extrapolation
    
#     # Calculated with square of uncertainty, than final uncertainty is np.sqrt()
#     err_int = sci.interpolate.CubicSpline(spectrum.spectral_axis[mask],
#                                            spectrum.uncertainty.array[mask]**2,
#                                            extrapolate= False)
#     new_flux = flux_int(new_spectral_axis) # Interpolate on the old wave_grid
#     import warnings
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         new_uncertainty = np.sqrt(err_int(new_spectral_axis)) # Interpolate on the old wave_grid
#     new_uncertainty = astropy.nddata.StdDevUncertainty(new_uncertainty) # Interpolate on the old wave_grid
#     new_spectrum = sp.Spectrum1D(
#         spectral_axis = new_spectral_axis,
#         flux = new_flux * spectrum.flux.unit,
#         uncertainty = new_uncertainty,
#         meta = spectrum.meta.copy(),
#         mask = np.isnan(new_flux),
#         wcs = spectrum.wcs,
#         )
    
#     if len(mask_region_list) !=0:
#         for region in mask_region_list:
#             new_spectrum = exciser_fill_with_nan(new_spectrum,region)
    
    
#     return new_spectrum

# #%% exciser_fill_with_nan
# def exciser_fill_with_nan(spectrum,region):
#     """
#     Takes a spectrum and fills given region with NaNs
#     Input:
#         spectrum ; sp.Spectrum1D - spectrum to mask
#         region ; sp.SpectralRegion - region to mask
#     Output:
#         new_spectrum ; sp.Spectrum1D - masked spectrum
    
#     """
#     spectral_axis = spectrum.spectral_axis
#     excise_indices = None

#     for subregion in region:
#         # Find the indices of the spectral_axis array corresponding to the subregion
#         region_mask = (spectral_axis >= region.lower) & (spectral_axis < region.upper)
#         region_mask = (spectral_axis >= subregion.lower) & (spectral_axis < subregion.upper)
#         temp_indices = np.nonzero(region_mask)[0]
#         if excise_indices is None:
#             excise_indices = temp_indices
#         else:
#             excise_indices = np.hstack((excise_indices, temp_indices))

#     new_spectral_axis = spectrum.spectral_axis.copy()
#     new_flux = spectrum.flux.copy()
#     modified_flux = new_flux
#     modified_flux[excise_indices] = np.nan
#     if spectrum.mask is not None:

#         new_mask = spectrum.mask
#         new_mask[excise_indices] = True
#     else:
#         new_mask = None
#     if spectrum.uncertainty is not None:

#         new_uncertainty = spectrum.uncertainty
#         # new_uncertainty[excise_indices] = np.nan
#     else:
#         new_uncertainty = None

#     # Return a new object with the regions excised.
#     return sp.Spectrum1D(flux=modified_flux,
#                       spectral_axis=new_spectral_axis,
#                       uncertainty=new_uncertainty,
#                       mask=new_mask,
#                       meta = spectrum.meta.copy(),
#                       wcs=spectrum.wcs,
#                       velocity_convention=spectrum.velocity_convention)


# #%% normalize_spectrum
# def normalize_spectrum(spectrum,quantile=.85,linfit=False):
#     """
#     Normalize spectrum depending on the size and linfit values.
#     Normalization function is either rolling quantile window (with size of 7500 pixels), or linear fit
    
#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Spectrum to normalize.
#     quantile : Float [0;1], optional
#         Quantile by which to normalize the spectrum. The default is .85.
#     linfit : bool, optional
#         Whether to fit by linear fit. Works only for spectra of length less than 10000. The default is False.

#     Returns
#     -------
#     normalized_spectrum : sp.Spectrum1D
#         Normalized spectrum based on parameters.

#     """
#     if (len(spectrum.flux) <10000) and linfit==True: 
#         p = np.polyfit(spectrum.spectral_axis.value[~np.isnan(spectrum.flux)],
#                          spectrum.flux.value[~np.isnan(spectrum.flux)],
#                          1,
#                          rcond=None,
#                          full=False,
#                          w=None,
#                          cov=False
#                          )
#         tmp_spec = sp.Spectrum1D(
#             spectral_axis = spectrum.spectral_axis,
#             flux = (np.polyval(p,spectrum.spectral_axis.value))*pd.Series(spectrum.flux.value / np.polyval(p,spectrum.spectral_axis.value)).quantile(quantile)*spectrum.flux.unit,
#             )
#     elif (len(spectrum.flux) <10000) and linfit==False:
#         tmp_spec = sp.Spectrum1D(
#             spectral_axis = spectrum.spectral_axis,
#             flux =np.full_like(spectrum.spectral_axis.value, pd.Series(spectrum.flux.value).fillna(999999).quantile(quantile))*spectrum.flux.unit,
#             )
#     else:
#         tmp_spec = sp.Spectrum1D(spectral_axis = spectrum.spectral_axis,
#                                  flux = np.array( pd.Series(spectrum.flux.value).rolling(7500
#                                    ,min_periods=1,center=True).quantile(quantile))*spectrum.flux.unit,
#                                   )
    
#     normalization = spectrum.divide(tmp_spec,
#                                     handle_mask = 'first_found',
#                                     handle_meta = 'first_found',
#                                     )
#     normalization.spectral_axis.value.put(np.arange(len(spectrum.spectral_axis)),spectrum.spectral_axis)
    
#     # Rewrite the spectral axis, as operations are done only on the flux array and WCS is dumb.
#     normalized_spectrum = sp.Spectrum1D(spectral_axis = spectrum.spectral_axis,
#                              flux = normalization.flux*u.dimensionless_unscaled,
#                              uncertainty = normalization.uncertainty,
#                              mask = spectrum.mask.copy(),
#                              meta = spectrum.meta.copy(),
#                              wcs = spectrum.wcs,
#                               )
#     normalized_spectrum.meta['normalization'] = True
#     return normalized_spectrum
# #%% calculate_master_list
# def calculate_master_list(spec_list,
#                           key = None,
#                           value = None,
#                           sn_type='quadratic',
#                           force_load = False,
#                           force_skip = False,
#                           pkl_name = ''
#                           ):
#     """
#     Calculates master list
#     Input:
#         spec_list ; sp.SpectrumList
#         key = 'Transit' ; key of meta dictionary
#         value = False ; value for which master should calculated (eg. 'Transit'==False for out-of-Transit master)
#         sn_type = 'S_N' ; type of weighting, options are: ('None','S_N','quadratic','quadratic_combined')
#     Ouput:
#         master_list ; sp.SpectrumList with length num_nights+1
#     Error:
#         When spec_list is not normalized
#     """
#     # Warning for non-normalized spectra
#     if (spec_list[0].meta['normalization'] == False) & (key != None):
#         message = 'Warning: Non-normalized spectra'
#     # Getting the right sublist based on type of master
#     if key == None:
#         # All spectra together
#         sublist = spec_list.copy()
#     else:
#         # Specified sublist 
#         sublist = get_sublist(spec_list,key,value)
#         # Necessary for master_all
#         spec_list = sublist
        
        
#     # What type of master we have
#     spec_type = get_spec_type(key,value)
#     # Creating master_list
#     master_list = sp.SpectrumList()
#     # Number of nights
#     num_nights = sublist[-1].meta['Night_num']
#     # Master of nights combined
#     master_all = get_master(sublist,
#                             spec_type,
#                             night='nights-combined',
#                             num_night = '"all"',
#                             rf = sublist[0].meta['RF'],
#                             sn_type=sn_type)
#     # Appending master_all to list
#     master_list.append(master_all)
#     # For cycle through nights
#     for ni in range(num_nights):
#         # Getting night data
#         sublist = get_sublist(spec_list,'Night_num',ni+1)
#         # Calculating master_night
#         master_night = get_master(sublist,
#                                   spec_type,
#                                   night=sublist[0].meta['Night'],
#                                   num_night = str(ni+1),
#                                   rf = sublist[0].meta['RF'],
#                                   sn_type=sn_type)
#         # Appending master_night
#         master_list.append(master_night)

#     return master_list

# #%% extract_velocity_field
# def extract_velocity_field(spectrum:sp.Spectrum1D,
#                            shift_BERV:float,
#                            shift_v_sys:float,
#                            shift_v_star:float,
#                            shift_v_planet:float,
#                            shift_constant=0,
#                            ):
#     """
#     Extracts velocity field for the shift
#     Input:
#         spectrum ; sp.Spectrum1D - spectrum from which to extract velocities
#         shift_BERV ; float - 1/0/-1, otherwise the velocity is scaled
#         shift_v_sys ; float - 1/0/-1, otherwise the velocity is scaled 
#         shift_v_star ; float - 1/0/-1, otherwise the velocity is scaled 
#         shift_v_planet ; float - 1/0/-1, otherwise the velocity is scaled
#         shift_constant ; float * u.m/u.s or equivalent
#     Output:
#         velocity_field ; list - list of velocities to shift by
#     """
#     velocity_field = []
#     if shift_BERV != 0:
#         velocity_field.append(spectrum.meta['BERV'] * shift_BERV)
#     if shift_v_sys != 0:
#         velocity_field.append(spectrum.meta['vel_sys'] * shift_v_sys)
#     if shift_v_star != 0:
#         velocity_field.append(spectrum.meta['vel_st'] * shift_v_star)
#     if shift_v_planet != 0:
#         velocity_field.append(spectrum.meta['vel_pl'] * shift_v_planet)
#     if shift_constant != 0:
#         velocity_field.append(shift_constant)
#     return velocity_field
# #%% shift_spectrum_multiprocessing
# def shift_spectrum_multiprocessing(spectrum,
#                            shift_BERV,
#                            shift_v_sys,
#                            shift_v_star,
#                            shift_v_planet,
#                            shift_constant):
#     """
#     Convenience function to pass to multiprocessing

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Spectrum to shift.
#     shift_BERV : float
#         Shift by BERV?
#     shift_v_sys : float
#         Shift by systemic velocity?
#     shift_v_star : float
#         Shift by stellar velocity?
#     shift_v_planet : float
#         Shift by planetary velocity?
#     shift_constant : float
#         Shift by arbitrary constant velocity?

#     Returns
#     -------
#     new_spectrum : sp.Spectrum1D
#         Shifted spectrum.

#     """
#     velocity_field =  extract_velocity_field(spectrum,
#                                shift_BERV = shift_BERV,
#                                shift_v_sys = shift_v_sys,
#                                shift_v_star = shift_v_star,
#                                shift_v_planet = shift_v_planet,
#                                shift_constant = shift_constant,
#                                )
#     new_spectrum = shift_spectrum(spectrum,velocity_field)
#     return new_spectrum

# #%% shift_spectrum
# def shift_spectrum(spectrum: sp.Spectrum1D | sp.SpectrumCollection, 
#                     velocities: list
#                     ) -> sp.Spectrum1D | sp.SpectrumCollection:
#     """
#     Shifts spectrum by a list of velocities.

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D | sp.SpectrumCollection
#         Spectrum to shift. The output will be the same as input
#     velocities : list
#         Velocity list to shift by. Must be a list of astropy Quantities in the units of velocity.

#     Returns
#     -------
#     new_spectrum : sp.Spectrum1D | sp.SpectrumCollection
#         Shifted spectrum, interpolated to the old wavelength grid.
#     """
#     term = 1 # Initiation of the term (1 since we are dividing/multiplying the term)
#     for velocity in velocities: # To factor all velocities (replace with relativistic one?)
#         term = term / (1+velocity.to(u.m/u.s)/con.c) # Divide as the velocities are in opposite sign
#     new_x_axis = spectrum.spectral_axis * term # Apply Doppler shift
    
#     # Mask handling
#     # Combine the mask arrays of flux and errors
#     mask_flux = ~np.isfinite(spectrum.flux)
#     mask_err = ~np.isfinite(spectrum.uncertainty.array)
#     mask = np.logical_or(mask_flux, mask_err)

#     # Interpolation function for flux - cubic spline with no extrapolation
#     interpolate_flux = sci.interpolate.CubicSpline(new_x_axis[~mask],
#                                            spectrum.flux[~mask],
#                                            extrapolate= False)
#     # Interpolation function for uncertainty - cubic spline with no extrapolation
    
#     # Calculated with square of uncertainty, than final uncertainty is np.sqrt()
#     interpolate_error = sci.interpolate.CubicSpline(new_x_axis[~mask],
#                                            spectrum.uncertainty.array[~mask]**2,
#                                            extrapolate= False)
#     # Applying interpolation functions to the old wave_grid
#     # mask = ~mask # Indices of good pixels (both flux and error)
#     new_flux = interpolate_flux(spectrum.spectral_axis) # Interpolate on the old wave_grid
    
#     import warnings
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         new_uncertainty = np.sqrt(interpolate_error(spectrum.spectral_axis)) # Interpolate on the old wave_grid
        
#     # Masking values that were NaN
#     new_flux[mask] = np.nan
#     new_uncertainty[mask] = np.nan
    
#     new_uncertainty = astropy.nddata.StdDevUncertainty(new_uncertainty) # Interpolate on the old wave_grid
    
#     # Create new spectrum
#     if type(spectrum) == sp.Spectrum1D:
#         new_spectrum = sp.Spectrum1D(
#             spectral_axis = spectrum.spectral_axis,
#             flux = new_flux * spectrum.flux.unit,
#             uncertainty = new_uncertainty,
#             meta = spectrum.meta.copy(),
#             mask = mask,
#             wcs = spectrum.wcs,
#             )
#     elif type(spectrum) == sp.SpectrumCollection:
#         logger.warning('Spectral Collection format has been untested, please verify')
#         new_spectrum = sp.SpectrumCollection(
#             spectral_axis = spectrum.spectral_axis,
#             flux = new_flux * spectrum.flux.unit,
#             uncertainty = new_uncertainty,
#             meta = spectrum.meta.copy(),
#             mask = mask,
#             wcs = spectrum.wcs,
#             )
#         logger.info('Finished correctly')
    
#     return new_spectrum
# #%% get_sublist
# def get_sublist(spec_list,key,value,mode='normal'):
#     """
#     Returns a sublist with given item.meta[key] == value
#     Input:
#         spec_list ; sp.SpectrumList
#         key ; string of the meta keyword
#         value ; what value should be in sublist (eg. Transit == True)
#         mode = 'normal' or 'equal'; Whether condition is ==, 
#             or < ('less'),
#             or > ('more'),
#             or != ('non-equal')
#     Output:
#         new_spec_list ; sp.SpectrumList sublist of spec_list
#     """
#     new_spec_list = spec_list.copy()
#     # For single value extraction
#     if (mode == 'normal') or (mode == 'equal'):
#         for item in spec_list:
#             if item.meta[key] != value:
#                 new_spec_list.remove(item)
#     # For values smaller than value
#     elif mode == 'less':
#         for item in spec_list:
#             if item.meta[key] > value:
#                 new_spec_list.remove(item)
#     # For values higher than value
#     elif mode == 'more':
#         for item in spec_list:
#             if item.meta[key] < value:
#                 new_spec_list.remove(item)
#     # For values that are non-equal to value
#     elif mode == 'non-equal':
#         for item in spec_list:
#             if item.meta[key] == value:
#                 new_spec_list.remove(item)
#     return new_spec_list
# #%% get_spec_type
# def get_spec_type(key,value):
#     """
#     Supplementary function giving spec_type value for each key,value used
#     Used for labeling plots
#     Input:
#         key ; key of meta dictionary
#         value ; value of meta dictionary
#     Output:
#         spec_type ; Type of master spectrum (eg. 'out-of-Transit master')
#     """
#     # Type of master based on in/out of transit
#     if (key == 'Transit') or (key == 'Transit_full') or (key == 'Preingress') or (key == 'Postegress'):
#         if value == False:
#             spec_type = 'Out-of-transit'
#         else:
#             spec_type = 'In-transit (transmission)'
            
#     # Type of master based on before/after telluric correction
#     if key == 'telluric_corrected':
#         if value == True:
#             spec_type = 'After-telluric-correction'
#         else:
#             spec_type = 'Before-telluric-correction'
#     # Set None master type (for debugging)
#     if key == None:
#         spec_type = 'None'
#     return spec_type


# #%% get_master
# def get_master(spec_list,spec_type = '',night = '',num_night = '',rf = '',sn_type='quadratic'):
#     """
#     Calculates master spectrum of spectrum list
#     Input:
#         spec_list ; sp.SpectrumList
#         spec_type ; Type of resulting master spectrum (for automatic labels)
#         night ; number of night/ 'all' nights (for automatic labels)
#         rf ; rest frame of spectrum (for automatic labels)
#         sn_type ; type of weighting
#             possible options are:
#                 'S_N' ; linear S_N weight
#                 'quadratic' ; quadratic S_N weight
#                 'quadratic_combined' ; quadratic S_N and light curve flux
#                 'None' ; fill with 1
#     Output:
#         master ; sp.Spectrum1D - master spectrum
#     """
#     # Unit of flux
#     unit_flux = spec_list[0].flux.unit
#     # Allocate wavelength, flux and flux_err arrays
#     spectral_axis = spec_list[0].spectral_axis
#     flux = np.zeros(spec_list[0].spectral_axis.shape)
#     flux_err = np.zeros(spec_list[0].spectral_axis.shape)
#     # Allocate weighting
#     # Since some pixels might be masked, its necessary to weight by pixel
#     weights_total = np.zeros(spec_list[0].spectral_axis.shape)
    
#     # For cycle through spec_list
#     for item in spec_list:
#         mask_flux = np.isnan(item.flux) # Masking NaNs
#         mask_err = np.isnan(item.uncertainty.array) # Masking NaNs
#         mask = mask_flux + mask_err # Getting combined mask (zero is for correct pixels)
#         mask != 0
#         # Taking flux and error from spectrum
#         tmp_flux = item.flux.value
#         tmp_err = item.uncertainty.array
#         # Assigning weights according to type of weighting
#         if sn_type == 'S_N':
#             weights = [item.meta['S_N']] * ~mask 
#         elif sn_type == 'quadratic':
#             weights = (np.asarray([item.meta['S_N']])**2) * ~mask
#         elif sn_type =='quadratic_error':
#             weights = ((np.asarray([item.meta['S_N']])**2*np.asarray([item.uncertainty.array])**2) * ~mask).flatten()
#         elif sn_type == 'quadratic_combined':
#             weights = (np.asarray([item.meta['S_N']])**2 +\
#                        np.asarray(item.meta['delta']**2)) * ~mask
#         elif sn_type == 'None':
#             weights = [1] * len(flux)
#         elif sn_type == 'quadratic_error':
#             weights = item.uncertainty.array**(-2) * ~mask
#         # Erasing NaN values with 0
#         tmp_flux = np.where(mask == False,tmp_flux,0)
#         tmp_err = np.where(mask == False,tmp_err,0)
#         weights = np.where(mask == False,weights,0)
#         # Suming flux, flux_err and weights for averaging
#         flux += tmp_flux * weights
#         flux_err += (tmp_err*weights)**2
#         weights_total += weights
#     # Averaging flux and flux_err
#     flux = flux / weights_total
#     flux_err = np.sqrt(flux_err) / weights_total
#     # Creation of master spectrum
#     master = sp.Spectrum1D(
#         flux = flux * unit_flux,
#         spectral_axis = spectral_axis,
#         uncertainty =  astropy.nddata.StdDevUncertainty(flux_err),
#         mask = np.isnan(flux),
#         )
#     # Updating type of master
#     master.meta = {'type':spec_type,
#                    'night':night,
#                    'Night_num':num_night,
#                    'RF': rf
#         }
#     return master

# #%% extract_subregion
# def extract_subregion(spectrum,subregion):
#     """
#     Convenience function that extracts indices of given subregion

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Spectrum from which to extract the subspectrum.
#     subregion : sp.SpectralRegion
#         Subregion of sp.SpectralRegion with size of 1.

#     Returns
#     -------
#     ind : array
#         Indices of which pixels to include.

#     """
#     ind = np.where(
#         np.logical_and(
#         spectrum.spectral_axis > subregion.lower,
#         spectrum.spectral_axis < subregion.upper,
#         )
#         )
#     return ind

# #%% extract_region
# def extract_region(spectrum,spectral_region):
#     """
#     Extract region from spectrum

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Spectrum from which to extract the subspectrum..
#     spectral_region : sp.SpectralRegion
#         Spectral region from which to take the subspectrum.

#     Returns
#     -------
#     cut_spectrum : sp.Spectrum1D
#         Cut spectrum with spectral_region.

#     """
#     ind = np.array([],dtype =int)
    
#     for subregion in spectral_region: # Extract all indices 
#         ind = np.append(ind,(extract_subregion(spectrum,subregion)))
#     cut_spectrum = sp.Spectrum1D( # Create new spectrum with old parameters
#         spectral_axis =spectrum.spectral_axis[ind],
#         flux = spectrum.flux[ind],
#         uncertainty = spectrum.uncertainty[ind],
#         )
#     return cut_spectrum
# #%% error_bin
# def error_bin(array):
#     """
#     Calculates error in a bin based on np.sqrt(sum(error_in_bin))
#     Input:
#         array ; np.array of values to bin together
#     Output:
#         value ; resulting value
#     """
#     # Change list to array
#     if isinstance(array,list):
#         array = np.asarray(array)
#     # For arrays longer than 1 value
#     if len(array) != 1:
#         value = np.sqrt(np.sum(array**2)/(len(array))**2)
#     # Else nothing changes
#     else:
#         value = array[0]
#     return value
# #%% binning_spectrum
# def binning_spectrum(spectrum,bin_factor = 10):
#     """
#     Bins a input spectrum by bin_factor*pixels

#     Parameters
#     ----------
#     spectrum : sp.Spectrum1D
#         Input spectrum to bin.
#     bin_factor : int, optional
#         How many pixels we want to bin by. The default is 10.

#     Returns
#     -------
#     x
#         x values of the binned spectrum.
#     y
#         y values of the binned spectrum.
#     y_err
#         y_err values of the binned spectrum.

#     """
#     num_bins = round(spectrum.spectral_axis.shape[0] / bin_factor)
#     x = sci.stats.binned_statistic(spectrum.spectral_axis, spectrum.spectral_axis, statistic='mean', bins=num_bins)
#     y = sci.stats.binned_statistic(spectrum.spectral_axis, spectrum.flux, statistic='mean', bins=num_bins)
#     y_err = sci.stats.binned_statistic(spectrum.spectral_axis, spectrum.uncertainty.array, statistic=error_bin, bins=num_bins)
#     return x.statistic,y.statistic,y_err.statistic
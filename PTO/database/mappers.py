

_NASA_EXOPLANET_ARCHIVE_COMPOSITE_MAPPER = {
    # Basic names
    'pl_name': 'Planet.Name',
    'hostname': 'Star.Name',
    'pl_letter': 'Planet.Letter',
    'hd_name': 'Star.Catalogue.HD',
    'hip_name': 'Star.Catalogue.HIP',
    'tic_id': 'Star.Catalogue.TIC',
    'gaia_id': 'Star.Catalogue.GAIA',
    # System composition
    'sy_snum': 'System.NumberOfStars',
    'sy_pnum': 'System.NumberOfPlanets',
    'sy_mnum': 'System.NumberOfMoons',
    'cb_flag': 'Flag.BinarySystem',
    # Planet discovery
    'discoverymethod': 'Discovery.Method',
    'disc_year': 'Discovery.Year',
    'disc_refname': 'Discovery.Reference',
    'disc_pubdate': 'Discovery.PublicationDate',
    'disc_locale': 'Discovery.Locale',
    'disc_facility': 'Discovery.Facility',
    'disc_telescope': 'Discovery.Telescope',
    'disc_instrument': 'Discovery.Instrument',
    # Detections
    'rv_flag': 'Flag.RadialVelocity',
    'pul_flag': 'Flag.PulsarTimingVariations',
    'ptv_flag': 'Flag.PulsationTimingVariations',
    'tran_flag': 'Flag.Transit',
    'ast_flag': 'Flag.AstrometricVariations',
    'obm_flag': 'Flag.OrbitalBrightnessModulations',
    'micro_flag': 'Flag.Microlensing',
    'etv_flag': 'Flag.EclipseTimingVariations',
    'ima_flag': 'Flag.Imaging',
    'dkin_flag': 'Flag.DiskKinematics',
    # Planet parameters
    'pl_refname': 'Planet.Reference',
    
    'pl_orbper': 'Planet.Period',
    'pl_orbpererr1': 'Planet.Period.Error.Upper',
    'pl_orbpererr2': 'Planet.Period.Error.Lower',
    'pl_orbper_reflink': 'Planet.Period.Reference',

    'pl_orbsmax': 'Planet.SemiMajorAxis',
    'pl_orbsmaxerr1': 'Planet.SemiMajorAxis.Error.Upper',
    'pl_orbsmaxerr2': 'Planet.SemiMajorAxis.Error.Lower',
    'pl_orbsmax_reflink': 'Planet.SemiMajorAxis.Reference',

    'pl_angsep': 'Planet.AngularSeparation',
    'pl_angseperr1': 'Planet.AngularSeparation.Error.Upper',
    'pl_angseperr2': 'Planet.AngularSeparation.Error.Lower',
    'pl_angsep_reflink': 'Planet.AngularSeparation.Reference',

    'pl_rade': 'Planet.RadiusEarth',
    'pl_radeerr1': 'Planet.RadiusEarth.Error.Upper',
    'pl_radeerr2': 'Planet.RadiusEarth.Error.Lower',
    'pl_rade_reflink': 'Planet.RadiusEarth.Reference',

    'pl_radj': 'Planet.RadiusJupiter',
    'pl_radjerr1': 'Planet.RadiusJupiter.Error.Upper',
    'pl_radjerr2': 'Planet.RadiusJupiter.Error.Lower',
    'pl_radj_reflink': 'Planet.RadiusJupiter.Reference',

    'pl_masse': 'Planet.MassEarth',
    'pl_masseerr1': 'Planet.MassEarth.Error.Upper',
    'pl_masseerr2': 'Planet.MassEarth.Error.Lower',
    'pl_masse_reflink': 'Planet.MassEarth.Reference',

    'pl_massj': 'Planet.MassJupiter',
    'pl_massjerr1': 'Planet.MassJupiter.Error.Upper',
    'pl_massjerr2': 'Planet.MassJupiter.Error.Lower',
    'pl_massj_reflink': 'Planet.MassJupiter.Reference',
    
    'pl_msinie': 'Planet.MinimumMassEarth',
    'pl_msinieerr1': 'Planet.MinimumMassEarth.Error.Upper',
    'pl_msinieerr2': 'Planet.MinimumMassEarth.Error.Lower',
    'pl_msinie_reflink': 'Planet.MinimumMassEarth.Reference',

    'pl_msinij': 'Planet.MinimumMassJupiter',
    'pl_msinijerr1': 'Planet.MinimumMassJupiter.Error.Upper',
    'pl_msinijerr2': 'Planet.MinimumMassJupiter.Error.Lower',
    'pl_msinij_reflink': 'Planet.MinimumMassJupiter.Reference',
    
    'pl_bmasse': 'Planet.BestMassEstimateEarth',
    'pl_bmasseerr1': 'Planet.BestMassEstimateEarth.Error.Upper',
    'pl_bmasseerr2': 'Planet.BestMassEstimateEarth.Error.Lower',
    'pl_bmasse_reflink': 'Planet.BestMassEstimateEarth.Reference',

    'pl_bmassj': 'Planet.BestMassEstimateJupiter',
    'pl_bmassjerr1': 'Planet.BestMassEstimateJupiter.Error.Upper',
    'pl_bmassjerr2': 'Planet.BestMassEstimateJupiter.Error.Lower',
    'pl_bmassj_reflink': 'Planet.BestMassEstimateJupiter.Reference',
    
    'pl_dens': 'Planet.Density',
    'pl_denserr1': 'Planet.Density.Error.Upper',
    'pl_denserr2': 'Planet.Density.Error.Lower',
    'pl_dens_reflink': 'Planet.Density.Reference',
    
    'pl_orbeccen': 'Planet.Eccentricity',
    'pl_orbeccenerr1': 'Planet.Eccentricity.Error.Upper',
    'pl_orbeccenerr2': 'Planet.Eccentricity.Error.Lower',
    'pl_orbeccen_reflink': 'Planet.Eccentricity.Reference',
    
    'pl_insol': 'Planet.InsolationFlux',
    'pl_insolerr1': 'Planet.InsolationFlux.Error.Upper',
    'pl_insolerr2': 'Planet.InsolationFlux.Error.Lower',
    'pl_insol_reflink': 'Planet.InsolationFlux.Reference',
    
    'pl_eqt': 'Planet.EquilibriumTemperature',
    'pl_eqterr1': 'Planet.EquilibriumTemperature.Error.Upper',
    'pl_eqterr2': 'Planet.EquilibriumTemperature.Error.Lower',
    'pl_eqt_reflink': 'Planet.EquilibriumTemperature.Reference',

    'pl_orbincl': 'Planet.Inclination',
    'pl_orbinclerr1': 'Planet.Inclination.Error.Upper',
    'pl_orbinclerr2': 'Planet.Inclination.Error.Lower',
    'pl_orbincl_reflink': 'Planet.Inclination.Reference',

    'pl_tranmid': 'Planet.TransitMidpoint',
    'pl_tranmiderr1': 'Planet.TransitMidpoint.Error.Upper',
    'pl_tranmiderr2': 'Planet.TransitMidpoint.Error.Lower',
    'pl_tranmid_reflink': 'Planet.TransitMidpoint.Reference',

    'pl_imppar': 'Planet.ImpactParameter',
    'pl_impparerr1': 'Planet.ImpactParameter.Error.Upper',
    'pl_impparerr2': 'Planet.ImpactParameter.Error.Lower',
    'pl_imppar_reflink': 'Planet.ImpactParameter.Reference',

    'pl_trandep': 'Planet.TransitDepth',
    'pl_trandeperr1': 'Planet.TransitDepth.Error.Upper',
    'pl_trandeperr2': 'Planet.TransitDepth.Error.Lower',
    'pl_trandep_reflink': 'Planet.TransitDepth.Reference',

    'pl_trandur': 'Planet.TransitDuration',
    'pl_trandurerr1': 'Planet.TransitDuration.Error.Upper',
    'pl_trandurerr2': 'Planet.TransitDuration.Error.Lower',
    'pl_trandur_reflink': 'Planet.TransitDuration.Reference',

    'pl_ratdor': 'Planet.RatioSemiMajorAxisToStellarRadius',
    'pl_ratdorerr1': 'Planet.RatioSemiMajorAxisToStellarRadius.Error.Upper',
    'pl_ratdorerr2': 'Planet.RatioSemiMajorAxisToStellarRadius.Error.Lower',
    'pl_ratdor_reflink': 'Planet.RatioSemiMajorAxisToStellarRadius.Reference',

    'pl_ratror': 'Planet.RatioPlanetRadiusToStellarRadius',
    'pl_ratrorerr1': 'Planet.RatioPlanetRadiusToStellarRadius.Error.Upper',
    'pl_ratrorerr2': 'Planet.RatioPlanetRadiusToStellarRadius.Error.Lower',
    'pl_ratror_reflink': 'Planet.RatioPlanetRadiusToStellarRadius.Reference',

    'pl_occdep': 'Planet.OccultationDepth',
    'pl_occdeperr1': 'Planet.OccultationDepth.Error.Upper',
    'pl_occdeperr2': 'Planet.OccultationDepth.Error.Lower',
    'pl_occdep_reflink': 'Planet.OccultationDepth.Reference',

    'pl_orbtper': 'Planet.EpochOfPeriastron',
    'pl_orbtpererr1': 'Planet.EpochOfPeriastron.Error.Upper',
    'pl_orbtpererr2': 'Planet.EpochOfPeriastron.Error.Lower',
    'pl_orbtper_reflink': 'Planet.EpochOfPeriastron.Reference',

    'pl_orblper': 'Planet.ArgumentOfPeriastron',
    'pl_orblpererr1': 'Planet.ArgumentOfPeriastron.Error.Upper',
    'pl_orblpererr2': 'Planet.ArgumentOfPeriastron.Error.Lower',
    'pl_orblper_reflink': 'Planet.ArgumentOfPeriastron.Reference',

    'pl_rvamp': 'Planet.RadialVelocityAmplitude',
    'pl_rvamperr1': 'Planet.RadialVelocityAmplitude.Error.Upper',
    'pl_rvamperr2': 'Planet.RadialVelocityAmplitude.Error.Lower',
    'pl_rvamp_reflink': 'Planet.RadialVelocityAmplitude.Reference',

    'pl_projobliq': 'Planet.ProjectedObliquity',
    'pl_projobliqerr1': 'Planet.ProjectedObliquity.Error.Upper',
    'pl_projobliqerr2': 'Planet.ProjectedObliquity.Error.Lower',
    'pl_projobliq_reflink': 'Planet.ProjectedObliquity.Reference',

    'pl_trueobliq': 'Planet.TrueObliquity',
    'pl_trueobliqerr1': 'Planet.TrueObliquity.Error.Upper',
    'pl_trueobliqerr2': 'Planet.TrueObliquity.Error.Lower',
    'pl_trueobliq_reflink': 'Planet.TrueObliquity.Reference',

    # Star parameters
    'st_refname': 'Star.Reference',
    'st_spectype': 'Star.Type',
    'st_spectype': 'Star.Type.Reference',
    
    'st_teff': 'Star.EffectiveTemperature',
    'st_tefferr1': 'Star.EffectiveTemperature.Error.Upper',
    'st_tefferr2': 'Star.EffectiveTemperature.Error.Lower',
    'st_teff_reflink': 'Star.EffectiveTemperature.Reference',
    
    'st_rad': 'Star.Radius',
    'st_raderr1': 'Star.Radius.Error.Upper',
    'st_raderr2': 'Star.Radius.Error.Lower',
    'st_rad_reflink': 'Star.Radius.Reference',
    
    'st_mass': 'Star.Mass',
    'st_masserr1': 'Star.Mass.Error.Upper',
    'st_masserr2': 'Star.Mass.Error.Lower',
    'st_mass_reflink': 'Star.Mass.Reference',

    'st_met': 'Star.Metallicity',
    'st_meterr1': 'Star.Metallicity.Error.Upper',
    'st_meterr2': 'Star.Metallicity.Error.Lower',
    'st_met_reflink': 'Star.Metallicity.Reference',

    'st_metratio': 'Star.MetallicityRatio',
    'st_metratioerr1': 'Star.MetallicityRatio.Error.Upper',
    'st_metratioerr2': 'Star.MetallicityRatio.Error.Lower',
    'st_metratio_reflink': 'Star.MetallicityRatio.Reference',

    'st_lum': 'Star.Luminosity',
    'st_lumerr1': 'Star.Luminosity.Error.Upper',
    'st_lumerr2': 'Star.Luminosity.Error.Lower',
    'st_lum_reflink': 'Star.Luminosity.Reference',
    
    'st_logg': 'Star.Logg',
    'st_loggerr1': 'Star.Logg.Error.Upper',
    'st_loggerr2': 'Star.Logg.Error.Lower',
    'st_logg_reflink': 'Star.Logg.Reference',

    'st_age': 'Star.Age',
    'st_ageerr1': 'Star.Age.Error.Upper',
    'st_ageerr2': 'Star.Age.Error.Lower',
    'st_age_reflink': 'Star.Age.Reference',

    'st_dens': 'Star.Density',
    'st_denserr1': 'Star.Density.Error.Upper',
    'st_denserr2': 'Star.Density.Error.Lower',
    'st_dens_reflink': 'Star.Density.Reference',

    'st_vsin': 'Star.RotationalVelocity',
    'st_vsinerr1': 'Star.RotationalVelocity.Error.Upper',
    'st_vsinerr2': 'Star.RotationalVelocity.Error.Lower',
    'st_vsin_reflink': 'Star.RotationalVelocity.Reference',

    'st_rotp': 'Star.RotationalPeriod',
    'st_rotperr1': 'Star.RotationalPeriod.Error.Upper',
    'st_rotperr2': 'Star.RotationalPeriod.Error.Lower',
    'st_rotp_reflink': 'Star.RotationalPeriod.Reference',

    # Technically under star, moved to System
    'st_radv': 'System.Velocity',
    'st_radverr1': 'System.Velocity.Error.Upper',
    'st_radverr2': 'System.Velocity.Error.Lower',
    'st_radv_reflink': 'System.Velocity.Reference',

    # System parameters
    'sy_refname': 'System.Reference',
    
    'sy_pm': 'System.TotalProperMotion',
    'sy_pmerr1': 'System.TotalProperMotion.Error.Upper',
    'sy_pmerr2': 'System.TotalProperMotion.Error.Lower',
    'sy_pm_reflink': 'System.TotalProperMotion.Reference',

    'sy_pmra': 'System.ProperMotionRightAscension',
    'sy_pmraerr1': 'System.ProperMotionRightAscension.Error.Upper',
    'sy_pmraerr2': 'System.ProperMotionRightAscension.Error.Lower',
    'sy_pmra_reflink': 'System.ProperMotionRightAscension.Reference',

    'sy_pmdec': 'System.ProperMotionDeclination',
    'sy_pmdecerr1': 'System.ProperMotionDeclination.Error.Upper',
    'sy_pmdecerr2': 'System.ProperMotionDeclination.Error.Lower',
    'sy_pmdec_reflink': 'System.ProperMotionDeclination.Reference',

    'sy_dist': 'System.Distance',
    'sy_disterr1': 'System.Distance.Error.Upper',
    'sy_disterr2': 'System.Distance.Error.Lower',
    'sy_dist_reflink': 'System.Distance.Reference',

    'sy_plx': 'System.Parallax',
    'sy_plxerr1': 'System.Parallax.Error.Upper',
    'sy_plxerr2': 'System.Parallax.Error.Lower',
    'sy_plx_reflink': 'System.Parallax.Reference',
    
    # Position keywords
    'rastr': 'Position.RightAscensionString',
    'decst': 'Position.DeclinationString',

    'ra': 'Position.RightAscension',
    'raerr1': 'Position.RightAscension.Error.Upper',
    'raerr2': 'Position.RightAscension.Error.Lower',
    'ra_reflink': 'Position.RightAscension.Reference',

    'dec': 'Position.Declination',
    'decerr1': 'Position.Declination.Error.Upper',
    'decerr2': 'Position.Declination.Error.Lower',
    'dec_reflink': 'Position.Declination.Reference',

    'glat': 'Position.GalacticLatitude',
    'glaterr1': 'Position.GalacticLatitude.Error.Upper',
    'glaterr2': 'Position.GalacticLatitude.Error.Lower',
    'glat_reflink': 'Position.GalacticLatitude.Reference',

    'glon': 'Position.GalacticLongitude',
    'glonerr1': 'Position.GalacticLongitude.Error.Upper',
    'glonerr2': 'Position.GalacticLongitude.Error.Lower',
    'glon_reflink': 'Position.GalacticLongitude.Reference',

    'elat': 'Position.EclipticLatitude',
    'elaterr1': 'Position.EclipticLatitude.Error.Upper',
    'elaterr2': 'Position.EclipticLatitude.Error.Lower',
    'elat_reflink': 'Position.EclipticLatitude.Reference',

    'elon': 'Position.EclipticLongitude',
    'elonerr1': 'Position.EclipticLongitude.Error.Upper',
    'elonerr2': 'Position.EclipticLongitude.Error.Lower',
    'elon_reflink': 'Position.EclipticLongitude.Reference',
    
    # Magnitude information
    'sy_bmag': 'Magnitude.B',
    'sy_bmagerr1': 'Magnitude.B.Error.Upper',
    'sy_bmagerr2': 'Magnitude.B.Error.Lower',
    'sy_bmag_reflink': 'Magnitude.B.Reference',

    'sy_vmag': 'Magnitude.V',
    'sy_vmagerr1': 'Magnitude.V.Error.Upper',
    'sy_vmagerr2': 'Magnitude.V.Error.Lower',
    'sy_vmag_reflink': 'Magnitude.V.Reference',

    'sy_jmag': 'Magnitude.J',
    'sy_jmagerr1': 'Magnitude.J.Error.Upper',
    'sy_jmagerr2': 'Magnitude.J.Error.Lower',
    'sy_jmag_reflink': 'Magnitude.J.Reference',

    'sy_hmag': 'Magnitude.H',
    'sy_hmagerr1': 'Magnitude.H.Error.Upper',
    'sy_hmagerr2': 'Magnitude.H.Error.Lower',
    'sy_hmag_reflink': 'Magnitude.H.Reference',

    'sy_kmag': 'Magnitude.K',
    'sy_kmagerr1': 'Magnitude.K.Error.Upper',
    'sy_kmagerr2': 'Magnitude.K.Error.Lower',
    'sy_kmag_reflink': 'Magnitude.K.Reference',

    'sy_umag': 'Magnitude.u',
    'sy_umagerr1': 'Magnitude.u.Error.Upper',
    'sy_umagerr2': 'Magnitude.u.Error.Lower',
    'sy_umag_reflink': 'Magnitude.u.Reference',

    'sy_gmag': 'Magnitude.g',
    'sy_gmagerr1': 'Magnitude.g.Error.Upper',
    'sy_gmagerr2': 'Magnitude.g.Error.Lower',
    'sy_gmag_reflink': 'Magnitude.g.Reference',

    'sy_rmag': 'Magnitude.r',
    'sy_rmagerr1': 'Magnitude.r.Error.Upper',
    'sy_rmagerr2': 'Magnitude.r.Error.Lower',
    'sy_rmag_reflink': 'Magnitude.r.Reference',

    'sy_imag': 'Magnitude.i',
    'sy_imagerr1': 'Magnitude.i.Error.Upper',
    'sy_imagerr2': 'Magnitude.i.Error.Lower',
    'sy_imag_reflink': 'Magnitude.i.Reference',

    'sy_zmag': 'Magnitude.z',
    'sy_zmagerr1': 'Magnitude.z.Error.Upper',
    'sy_zmagerr2': 'Magnitude.z.Error.Lower',
    'sy_zmag_reflink': 'Magnitude.z.Reference',

    'sy_w1mag': 'Magnitude.W1',
    'sy_w1magerr1': 'Magnitude.W1.Error.Upper',
    'sy_w1magerr2': 'Magnitude.W1.Error.Lower',
    'sy_w1mag_reflink': 'Magnitude.W1.Reference',

    'sy_w2mag': 'Magnitude.W2',
    'sy_w2magerr1': 'Magnitude.W2.Error.Upper',
    'sy_w2magerr2': 'Magnitude.W2.Error.Lower',
    'sy_w2mag_reflink': 'Magnitude.W2.Reference',

    'sy_w3mag': 'Magnitude.W3',
    'sy_w3magerr1': 'Magnitude.W3.Error.Upper',
    'sy_w3magerr2': 'Magnitude.W3.Error.Lower',
    'sy_w3mag_reflink': 'Magnitude.W3.Reference',

    'sy_w4mag': 'Magnitude.W4',
    'sy_w4magerr1': 'Magnitude.W4.Error.Upper',
    'sy_w4magerr2': 'Magnitude.W4.Error.Lower',
    'sy_w4mag_reflink': 'Magnitude.W4.Reference',

    'sy_gaiamag': 'Magnitude.Gaia',
    'sy_gaiamagerr1': 'Magnitude.Gaia.Error.Upper',
    'sy_gaiamagerr2': 'Magnitude.Gaia.Error.Lower',
    'sy_gaiamag_reflink': 'Magnitude.Gaia.Reference',

    'sy_icmag': 'Magnitude.IC',
    'sy_icmagerr1': 'Magnitude.IC.Error.Upper',
    'sy_icmagerr2': 'Magnitude.IC.Error.Lower',
    'sy_icmag_reflink': 'Magnitude.IC.Reference',

    'sy_tmag': 'Magnitude.TESS',
    'sy_tmagerr1': 'Magnitude.TESS.Error.Upper',
    'sy_tmagerr2': 'Magnitude.TESS.Error.Lower',
    'sy_tmag_reflink': 'Magnitude.TESS.Reference',

    'sy_kepmag': 'Magnitude.Kepler',
    'sy_kepmagerr1': 'Magnitude.Kepler.Error.Upper',
    'sy_kepmagerr2': 'Magnitude.Kepler.Error.Lower',
    'sy_kepmag_reflink': 'Magnitude.Kepler.Reference',

}
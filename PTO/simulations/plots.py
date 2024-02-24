#%% Importing libraries
import matplotlib.pyplot as plt
import specutils as sp
import seaborn as sns



#%%
def get_colors(color: str):
    """
    Get colors given a string

    Parameters
    ----------
    color : str
        Name of the color. 
        The choices are:
            'blue'
            'orange'
            'green'
            'red'
            'purple'
            'brown'
            'pink'
            'black'
            'gold'
            'cyan'

    Returns
    -------
    color_spectrum : color
        Color to use for the spectrum.
    color_bin : color
        Color of the bins to go with the spectrum.

    """
    
    if color == 'blue':
        color_spectrum = sns.color_palette("pastel")[0]
        color_bin = sns.color_palette("dark")[0]
    elif color == 'orange':
        color_spectrum = sns.color_palette("pastel")[1]
        color_bin = sns.color_palette("dark")[1]
    elif color == 'green':
        color_spectrum = sns.color_palette("pastel")[2]
        color_bin = sns.color_palette("dark")[2]
    elif color == 'red':
        color_spectrum = sns.color_palette("pastel")[3]
        color_bin = sns.color_palette("dark")[3]
    elif color == 'purple':
        color_spectrum = sns.color_palette("pastel")[4]
        color_bin = sns.color_palette("dark")[4]
    elif color == 'brown':
        color_spectrum = sns.color_palette("pastel")[5]
        color_bin = sns.color_palette("dark")[5]
    elif color == 'pink':
        color_spectrum = sns.color_palette("pastel")[6]
        color_bin = sns.color_palette("dark")[6]
    elif color == 'black':
        color_spectrum = sns.color_palette("pastel")[7]
        color_bin = sns.color_palette("dark")[7]
    elif color == 'gold':
        color_spectrum = sns.color_palette("pastel")[8]
        color_bin = sns.color_palette("dark")[8]
    elif color == 'cyan':
        color_spectrum = sns.color_palette("pastel")[9]
        color_bin = sns.color_palette("dark")[9]
    else:
        raise ValueError('Color not found.')
    return color_spectrum, color_bin


#%%
def plot_transmission_spectrum(spectrum:sp.Spectrum1D, color:str = 'green')->[]:
    """
    Plots the final transmission spectrum

    Parameters
    ----------
    spectrum : sp.Spectrum1D
        Transmission spectrum to plot.
    color : str, optional
        Color to use within the transmission spectrum. The default is 'green'.
        This will throw an error if the color is not correct.
        Refer to sts.get_color documentation for more information.

    Returns
    -------
    None.

    """
    
    
    color_spectrum, color_bin = get_colors(color)
    
    
    fig = plt.figure(constrained_layout= True, figsize= (20,6))
    spec_grid = gs.GridSpec(ncols=1,nrows=6)
    ax = fig.add_subplot(spec_grid[:4,0])
    ax2 = fig.add_subplot(spec_grid[4,0],sharex=ax)
    
    ax.errorbar(spectrum.spectral_axis.value,
                spectrum.flux.value - 1 ,
                spectrum.uncertainty.array,
                color = color_spectrum,
                label = '_nolegend_',
                fmt='.',
                markersize=10,
                elinewidth=1,
                alpha=0.3
                      )

    binning_factor = 15
    from PTO.spectra_manipulation import binning_spectrum
    x,y,yerr = binning_spectrum(spectrum,binning_factor) 
    ax.errorbar(x,
                y - 1,
                yerr,
                color = color_bin,
                label = '_nolegend_',
                fmt='.',
                markersize='20',
                markeredgecolor='black',
                markeredgewidth = 1,
                alpha = 0.5
                )
    
    from PTO.spectra_manipulation import extract_region
    sr = sp.SpectralRegion(5887*u.AA, 5897*u.AA)
    fit_spectrum = extract_region(spectrum, sr)
    try:
        popt, pcov = sci.optimize.curve_fit(double_gaussian,
                                        fit_spectrum.spectral_axis.value,
                                        fit_spectrum.flux.value - 1,
                                        p0 = [
                                            5889.950,
                                            0.01,
                                            0.1,
                                            5895.924,
                                            0.01,
                                            0.1
                                            ],
                                        sigma = fit_spectrum.uncertainty.array,
                                        )
        ax.plot(spectrum.spectral_axis.value,
                double_gaussian(spectrum.spectral_axis.value, *popt),
                color = color_bin,
                zorder = 5
                )

        perr = np.sqrt(np.diag(pcov))
        # logger.info(popt[1]/ perr[1], popt[4]/perr[4])
        # logger.info('{amplitude_D2:.4f} ± {error_D2:.4f}'.format(amplitude_D2 = popt[1], error_D2 = perr[1]))
        # logger.info('{amplitude_D1:.4f} ± {error_D1:.4f}'.format(amplitude_D1 = popt[4], error_D1 = perr[4]))

        ax2.errorbar(x,
                    y - 1 - double_gaussian(x, *popt),
                    yerr,
                    color = color_bin,
                    label = '_nolegend_',
                    fmt='.',
                    markersize='20',
                    markeredgecolor='black',
                    markeredgewidth = 1,
                    alpha = 0.5,
                    zorder = 5
                    )
        
    except RuntimeError:
        popt = [
                5889.950,
                0.,
                0.,
                5895.924,
                0.,
                0.
                ]
        logger.critical('Failed to fit the spectrum.')
    

    ax2.errorbar(spectrum.spectral_axis.value,
          (spectrum.flux.value - 1  - double_gaussian(spectrum.spectral_axis.value, *popt)),
          spectrum.uncertainty.array,
          color = color_spectrum,
          label = '_nolegend_',
          fmt='.',
          markersize=10,
          elinewidth=1,
          alpha=0.3
          )
    ax.axvline(5889.950)
    ax.axvline(5895.924)
    ax.set_xlim(5888,5898)
    ax2.set_xlabel('Wavelength [$\AA$]', fontsize = 24)
    ax.set_ylabel('Excess absorption [%]', fontsize = 24)
    ax2.set_ylabel('O-C [%]', fontsize = 24)
    
    
    ax.set_xticks([])
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0))
    ax2.xaxis.set_tick_params(labelbottom=True)
    ax2.set_xticks([5888, 5890, 5892, 5894, 5896, 5898])
    
    return fig, ax, ax2
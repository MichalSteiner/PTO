import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from dataclasses import dataclass
import logging
import pandas as pd
from ..utils.utilities import logger_default
from .mappers import LABEL_MAPPER, SCALE_MAPPER
import matplotlib.colors as mcolors

logger = logging.getLogger(__name__)
if not logger.handlers:
    logger = logger_default(logger)


@dataclass
class ColorPopulationDiagram:
    """
    A class to define a theme for a color population diagram.
    
    Attributes
    ----------
    theme : str
        The name of the theme.
    cmap : str
        The colormap used for the diagram.
    scatter : dict
        The scatter plot data kwargs.
    highlight_scatter : dict
        The highlighted scatter plot kwargs.
    """
    
    theme: str
    cmap: str

    scatter: dict
    highlight_scatter: dict


RedPopulationDiagram = ColorPopulationDiagram(
    theme='red',
    cmap=mcolors.ListedColormap(sns.color_palette(
        "Reds", as_cmap=True)(np.linspace(0.2, 1, 256))),
    scatter={
        'color': sns.color_palette('muted')[3],
        'edgecolor': 'black',
        'alpha': 0.25,
        's': 30
    },

    highlight_scatter={
        'color': sns.color_palette('bright')[0],
        'edgecolor': 'black',
        's': 200,
        'zorder': 10
    }
)

BluePopulationDiagram = ColorPopulationDiagram(
    theme='blue',
    cmap=mcolors.ListedColormap(sns.color_palette(
        "Blues", as_cmap=True)(np.linspace(0.2, 1, 256))),
    scatter={
        'color': sns.color_palette('dark')[0],
        'edgecolor': 'black',
        'alpha': 0.2,
        's': 30,
    },
    highlight_scatter={
        'color': sns.color_palette('bright')[2],
        'edgecolor': 'black',
        's': 200,
        'zorder': 10
    }
)

GreenPopulationDiagram = ColorPopulationDiagram(
    theme='green',
    cmap=mcolors.ListedColormap(sns.color_palette(
        "Greens", as_cmap=True)(np.linspace(0.2, 1, 256))),
    scatter={
        'color': sns.color_palette('dark')[2],
        'edgecolor': 'black',
        'alpha': 0.2,
        's': 30,
    },
    highlight_scatter={
        'color': sns.color_palette('bright')[0],
        'edgecolor': 'black',
        's': 200,
        'zorder': 10
    }
)

GreyScalePopulationDiagram = ColorPopulationDiagram(
    theme='grayscale',
    cmap='Greys',
    scatter={
        'color': 'black',
        'edgecolor': 'black',
        'alpha': 0.1,
        's': 30,
    },
    highlight_scatter={
        'color': sns.color_palette('bright')[0],
        'edgecolor': 'black',
        's': 200,
        'zorder': 10
    }
)

PurplePopulationDiagram = ColorPopulationDiagram(
    theme='purple',
    cmap=mcolors.ListedColormap(sns.color_palette(
        "Purples", as_cmap=True)(np.linspace(0.2, 1, 256))),
    scatter={
        'color': sns.color_palette('dark')[4],
        'edgecolor': 'black',
        'alpha': 0.2,
        's': 30,
    },

    highlight_scatter={
        'color': sns.color_palette('bright')[1],
        'edgecolor': 'black',
        's': 200,
        'zorder': 10
    }
)


def _print_PopulationDiagramTheme():
    """
    Collects and logs the themes of all instances of ColorPopulationDiagram found in the global scope.
    
    This function iterates through all global variables, identifies instances of the 
    ColorPopulationDiagram class, logs their themes using the logger, and collects these themes 
    into a list which is then returned.
    
    Returns:
        list: A list of themes from all ColorPopulationDiagram instances found in the global scope.
    """
    themes = []
    for var_value in globals().values():
        if isinstance(var_value, ColorPopulationDiagram):
            logger.info(f"{var_value.theme}")
            themes.append(var_value.theme)
    return themes


def _get_PopulationDiagramTheme(theme: str):
    """
    Returns the appropriate PopulationDiagram class based on the provided theme.
    
    Parameters
    ----------
    theme : str
        The theme for the population diagram. Valid options are:
                 'red', 'green', 'blue', 'purple', 'greyscale', 'grayscale', 'grey', 'gray'.
    
    Returns
    -------
    ColorPopulationDiagram
        The corresponding PopulationDiagram class for the given theme.
    
    Raises
    ------
    ValueError
        If the provided theme is not valid.
    """
    
    
    match theme:
        case 'red':
            return RedPopulationDiagram
        case 'green':
            return GreenPopulationDiagram
        case 'blue':
            return BluePopulationDiagram
        case 'purple':
            return PurplePopulationDiagram
        case 'greyscale' | 'grayscale' | 'grey' | 'gray':
            return GreyScalePopulationDiagram
        case _:
            logger.warning('Invalid theme. Valid options are:')
            _print_PopulationDiagramTheme()
            raise ValueError('Not a valid theme')


class PlotUtilitiesComposite():
    """
    A utility class for creating and customizing population diagrams and highlighting samples on plots.
    Methods
    -------
    _set_labels(ax: plt.Axes, x_key: str, y_key: str)
        Sets the x and y labels for the given axes based on the provided keys.
    _set_scales(ax: plt.Axes, x_key: str, y_key: str)
        Sets the x and y scales for the given axes based on the provided keys.
    plot_population_diagram(x_key: str, y_key: str, ax: plt.Axes | None = None, fig: plt.Figure | None = None, theme: str | ColorPopulationDiagram = 'blue', **kwargs) -> [plt.Figure, plt.Axes]
        Plots a population diagram on the given axes or creates new ones if not provided. Allows customization through themes and additional keyword arguments.
    highlight_sample(x_key: str, y_key: str, ax: plt.Axes | None = None, fig: plt.Figure | None = None, theme: str | ColorPopulationDiagram = 'blue', **kwargs) -> [plt.Figure, plt.Axes]
        Highlights a sample on the population diagram. If axes are not provided, it creates a new population diagram first.
    available_themes()
        Prints and returns the available themes for the population diagram.
    """
    

    def _set_labels(self,
                    ax: plt.Axes,
                    x_key: str,
                    y_key: str):
        """
        Sets the x and y labels for the given axes based on the provided keys.
        
        Parameters
        ----------
        ax : plt.Axes
            The matplotlib Axes object to set the labels on.
        x_key : str
            The key to retrieve the x-axis label from LABEL_MAPPER.
        y_key : str
            The key to retrieve the y-axis label from LABEL_MAPPER.
        """
        ax.set_xlabel(LABEL_MAPPER[x_key])
        ax.set_ylabel(LABEL_MAPPER[y_key])

    def _set_scales(self,
                    ax: plt.Axes,
                    x_key: str,
                    y_key: str):
        """
        Sets the x and y scales for the given axes based on the provided keys.

        Parameters
        ----------
        ax : plt.Axes
            The matplotlib Axes object to set the scales on.
        x_key : str
            The key to retrieve the x-axis scale from SCALE_MAPPER.
        y_key : str
            The key to retrieve the y-axis scale from SCALE_MAPPER.
        """
        ax.set_xscale(SCALE_MAPPER[x_key])
        ax.set_yscale(SCALE_MAPPER[y_key])

    def plot_population_diagram(self,
                     x_key: str,
                     y_key: str,
                     ax: plt.Axes | None = None,
                     fig: plt.Figure | None = None,
                     theme: str | ColorPopulationDiagram = 'blue',
                     **kwargs
                     ) -> [plt.Figure, plt.Axes]:
        """
        Plots a population diagram on the given axes or creates new ones if not provided.
        
        This method generates a population diagram using the specified x and y keys for the data.
        It allows customization through themes and additional keyword arguments for the scatter plot.

        Parameters
        ----------
        x_key : str
            The key to retrieve the x-axis data from the legacy_table.
        y_key : str
            The key to retrieve the y-axis data from the legacy_table.
        ax : plt.Axes | None, optional
            The matplotlib Axes object to plot on. If None, a new figure and axes are created.
        fig : plt.Figure | None, optional
            The matplotlib Figure object to plot on. If None, a new figure is created.
        theme : str | ColorPopulationDiagram, optional
            The theme for the population diagram. Can be a string representing the theme name or an instance of ColorPopulationDiagram. Default is 'blue'.
        **kwargs
            Additional keyword arguments for the scatter plot.

        Returns
        -------
        [plt.Figure, plt.Axes]
            The matplotlib Figure and Axes objects with the population diagram.
        """
        
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(12, 8))
            self._set_scales(ax, x_key, y_key)
            self._set_labels(ax, x_key, y_key)

        nan_indice = self.legacy_table[x_key].index.notna(
        ) & self.legacy_table[y_key].index.notna()

        if type(theme) == str:
            Theme = _get_PopulationDiagramTheme(theme)
        else:
            Theme = theme

        for key in Theme.scatter.keys():
            if key not in kwargs:
                kwargs[key] = Theme.scatter[key]

        sns.kdeplot(
            x=self.legacy_table[x_key][nan_indice],
            y=self.legacy_table[y_key][nan_indice],
            fill=True,
            thresh=0,
            levels=50,
            cmap=Theme.cmap,
            ax=ax,
            log_scale=(True, True)
        )

        ax.scatter(
            x=self.legacy_table[x_key][nan_indice],
            y=self.legacy_table[y_key][nan_indice],
            label= '__nolegend__',
            **kwargs
        )

        return fig, ax

    def highlight_sample(self,
                         x_key: str,
                         y_key: str,
                         ax: plt.Axes | None = None,
                         fig: plt.Figure | None = None,
                         theme: str | ColorPopulationDiagram = 'blue',
                         **kwargs
                         ) -> [plt.Figure, plt.Axes]:
        """
        Highlights a sample on the population diagram.

        This method highlights specific data points on an existing population diagram or creates a new one if axes are not provided. It uses the specified x and y keys for the data and allows customization through themes and additional keyword arguments for the scatter plot.

        Parameters
        ----------
        x_key : str
            The key to retrieve the x-axis data from the table.
        y_key : str
            The key to retrieve the y-axis data from the table.
        ax : plt.Axes | None, optional
            The matplotlib Axes object to plot on. If None, a new population diagram is created.
        fig : plt.Figure | None, optional
            The matplotlib Figure object to plot on. If None, a new figure is created.
        theme : str | ColorPopulationDiagram, optional
            The theme for highlighting the sample. Can be a string representing the theme name or an instance of ColorPopulationDiagram. Default is 'blue'.
        **kwargs
            Additional keyword arguments for the scatter plot.

        Returns
        -------
        [plt.Figure, plt.Axes]
            The matplotlib Figure and Axes objects with the highlighted sample.
        """

        if type(theme) == str:
            Theme = _get_PopulationDiagramTheme(theme)
        else:
            Theme = theme

        if ax is None:
            fig, ax = self.plot_population_diagram(
                x_key=x_key,
                y_key=y_key,
                theme=Theme
            )

        nan_indice = self.table[x_key].index.notna(
        ) & self.table[y_key].index.notna()

        for key in Theme.highlight_scatter.keys():
            if key not in kwargs:
                kwargs[key] = Theme.highlight_scatter[key]

        ax.scatter(
            x=self.table[x_key][nan_indice],
            y=self.table[y_key][nan_indice],
            **kwargs,
        )
        return fig, ax

    def available_themes(self):

        logger.print('='*25)
        logger.print('Printing themes for the plot_diagram() method')
        themes = _print_PopulationDiagramTheme()
        logger.print('='*25)
        return themes

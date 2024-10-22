import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from dataclasses import dataclass
import logging
import pandas as pd
from ..utils.utilities import logger_default

sns.set_context('talk')

logger = logging.getLogger(__name__)
logger = logger_default(logger) 



@dataclass
class ColorPopulationDiagram:
    theme: str
    cmap: str
    
    colorscatter: str
    edgecolorscatter: str
    alphascatter: float
    
    highlight_scatter: str

RedPopulationDiagram = ColorPopulationDiagram(
    theme= 'red',
    cmap='Oranges',
    colorscatter='darkred',
    edgecolorscatter='black',
    alphascatter=0.1,
    highlight_scatter= sns.color_palette('dark')[0]
)
BluePopulationDiagram = ColorPopulationDiagram(
    theme='blue',
    cmap='Blues',
    colorscatter='darkblue',
    edgecolorscatter='black',
    alphascatter=0.1,
    highlight_scatter= sns.color_palette('dark')[1]
)
GreenPopulationDiagram = ColorPopulationDiagram(
    theme='green',
    cmap='Greens',
    colorscatter='darkgreen',
    edgecolorscatter='black',
    alphascatter=0.1,
    highlight_scatter= sns.color_palette('dark')[0]
)
GreyScalePopulationDiagram = ColorPopulationDiagram(
    theme='grayscale',
    cmap='Greys',
    colorscatter='black',
    edgecolorscatter='white',
    alphascatter=0.1,
    highlight_scatter= sns.color_palette('dark')[1]
)
PurplePopulationDiagram = ColorPopulationDiagram(
    theme='purple',
    cmap='Purples',
    colorscatter='indigo',
    edgecolorscatter='black',
    alphascatter=0.1,
    highlight_scatter= sns.color_palette('dark')[0]
)

def _print_PopulationDiagramTheme():
    themes = []
    for var_value in globals().values():
        if isinstance(var_value, ColorPopulationDiagram):
            logger.info(f"{var_value.theme}")
            themes.append(var_value.theme)
    return themes

def _get_PopulationDiagramTheme(theme:str):
    match theme:
        case 'red':
            return RedPopulationDiagram
        case 'green':
            return GreenPopulationDiagram
        case 'blue':
            return BluePopulationDiagram
        case 'purple':
            return PurplePopulationDiagram
        case 'greyscale' | 'grayscale':
            return GreyScalePopulationDiagram
        case _:
            logger.warning('Invalid theme. Valid options are:')
            _print_PopulationDiagramTheme()
            raise ValueError('Not a valid theme')

def _set_log_scales(x_key: str, y_key:str, ax):
    ...


class PlotUtilitiesComposite():
    
    def plot_diagram(self,
                     x_key: str,
                     y_key: str,
                     ax: plt.Axes | None = None,
                     fig: plt.Figure | None = None,
                     theme:str | ColorPopulationDiagram = 'red'
                     ) -> [plt.Figure, plt.Axes]:
        
        ax = None
        if ax is None:
            fig, ax = plt.subplots(1, figsize=(12,8))
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        nan_indice = self.legacy_table[x_key].index.notna() & self.legacy_table[y_key].index.notna()
        
        if type(theme) == str:
            Theme = _get_PopulationDiagramTheme(theme)
        else:
            Theme = theme
        
        sns.kdeplot(
            x=self.legacy_table[x_key][nan_indice],
            y=self.legacy_table[y_key][nan_indice],
            fill=True,
            thresh=0,
            levels=40,
            cmap=Theme.cmap,
            ax= ax,
        )
        
        ax.scatter(
            x=self.legacy_table[x_key][nan_indice],
            y=self.legacy_table[y_key][nan_indice],
            color= Theme.colorscatter,
            alpha=Theme.alphascatter,
            edgecolors= Theme.edgecolorscatter
        )
        
        return fig, ax

    def highlight_sample(self,
                         x_key: str,
                         y_key: str,
                         ax: plt.Axes | None = None,
                         fig: plt.Figure | None = None,
                         theme:str | ColorPopulationDiagram = 'red',
                         **kwargs
                         ) -> [plt.Figure, plt.Axes]:
        
        if type(theme) == str:
            Theme = _get_PopulationDiagramTheme(theme)
        else:
            Theme = theme
        
        if ax is None:
            fig, ax = self.plot_diagram(
                x_key=x_key,
                y_key=y_key,
                theme=Theme
            )
        
        nan_indice = self.table[x_key].index.notna() & self.table[y_key].index.notna()
        
        ax.scatter(
            x=self.table[x_key][nan_indice],
            y=self.table[y_key][nan_indice],
            color= Theme.highlight_scatter,
            s= 100,
            **kwargs,
        )
        
        return fig, ax 

    def available_themes(self):
        
        logger.print('='*25)
        logger.print('Printing themese for the plot_diagram() method')
        themes = _print_PopulationDiagramTheme()
        logger.print('='*25)
        return themes

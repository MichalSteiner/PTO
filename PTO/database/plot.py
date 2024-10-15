import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 


@dataclass
class ColorPopulationDiagram:
    theme: str
    cmap: str
    colorscatter: str
    edgecolorscatter: str
    alphascatter: float

RedPopulationDiagram = ColorPopulationDiagram(
    theme= 'red',
    cmap='Oranges',
    colorscatter='darkred',
    edgecolorscatter='black',
    alphascatter=0.1,
)
BluePopulationDiagram = ColorPopulationDiagram(
    theme='blue',
    cmap='Blues',
    colorscatter='darkblue',
    edgecolorscatter='black',
    alphascatter=0.1,
)
GreenPopulationDiagram = ColorPopulationDiagram(
    theme='green',
    cmap='Greens',
    colorscatter='darkgreen',
    edgecolorscatter='black',
    alphascatter=0.1,
)
GreyScalePopulationDiagram = ColorPopulationDiagram(
    theme='grayscale',
    cmap='Greys',
    colorscatter='black',
    edgecolorscatter='white',
    alphascatter=0.1,
)
PurplePopulationDiagram = ColorPopulationDiagram(
    theme='purple',
    cmap='Purples',
    colorscatter='indigo',
    edgecolorscatter='black',
    alphascatter=0.1,
)

def _print_PopulationDiagramTheme():
    logger.info('red')
    logger.info('green')
    logger.info('blue')
    logger.info('purple')
    logger.info('greyscale')
    logger.info('grayscale')
    

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
        
        nan_indice = self.table[x_key].index.notna() & self.table[y_key].index.notna()
        
        if type(theme) == str:
            Theme = _get_PopulationDiagramTheme(theme)
        
        sns.kdeplot(
            x=self.table[x_key][nan_indice],
            y=self.table[y_key][nan_indice],
            fill=True,
            thresh=0,
            levels=100,
            cmap=Theme.cmap,
            ax= ax,
        )
        
        ax.scatter(
            x=self.table[x_key][nan_indice],
            y=self.table[y_key][nan_indice],
            color= Theme.colorscatter,
            alpha=Theme.alphascatter,
            edgecolors= Theme.edgecolorscatter
        )
        
        return fig, ax

    def available_themes(self):
        
        logger.print('='*25)
        logger.print('Printing themese for the plot_diagram() method')
        _print_PopulationDiagramTheme()
        logger.print('='*25)

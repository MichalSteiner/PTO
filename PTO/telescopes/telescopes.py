from . import instruments
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from dataclasses import dataclass, field
import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 


@dataclass
class Telescope():
    name: str
    location: coord.earth.EarthLocation
    instruments: list = field(default_factory=list)
    diameter: u.Quantity = 0 *u.m,
    operational: bool = True


# Calar Alto Observatory (CARMENES)
CalarAlto = Telescope(
    name='Calar Alto Observatory',
    location=coord.EarthLocation.of_site('Observatorio de Calar Alto'),
    instruments=[instruments.CARMENES],
    diameter=3.5 * u.m
)

# La Silla Observatory, Swiss Euler Telescope (CORALIE)
SwissEuler = Telescope(
    name='Swiss Euler Telescope',
    location=coord.EarthLocation.of_site('La Silla Observatory'),
    instruments=[instruments.CORALIE],
    diameter=1.2 * u.m
)

# La Silla Observatory, ESO 3.6-meter Telescope (HARPS, NIRPS)
LaSilla36m = Telescope(
    name='ESO 3.6m Telescope (La Silla)',
    location=coord.EarthLocation.of_site('La Silla Observatory (ESO)'),
    instruments=[instruments.HARPS, instruments.NIRPS],
    diameter=3.6 * u.m
)

# Paranal Observatory, VLT Unit Telescope (ESPRESSO, UVES)
VLT = Telescope(
    name='Very Large Telescope (VLT) Unit Telescope',
    location=coord.EarthLocation.of_site('Paranal'),
    instruments=[instruments.ESPRESSO, instruments.UVES],
    diameter=8.2 * u.m
)

# Paranal Observatory, VLT 4-UT (ESPRESSO 4UT)
VLT_4UT = Telescope(
    name='Very Large Telescope (VLT) 4-Unit Telescope',
    location=coord.EarthLocation.of_site('Paranal'),
    instruments=[instruments.ESPRESSO_4UT],
    diameter=2 * 8.2 * u.m  # Combined light from all 4 UTs
)

# Lowell Observatory, Lowell Discovery Telescope (EXPRES)
LowellDiscovery = Telescope(
    name='Lowell Discovery Telescope',
    location=coord.EarthLocation.of_site('Lowell Observatory'),
    instruments=[instruments.EXPRES],
    diameter=4.3 * u.m
)

# Roque de los Muchachos Observatory, Telescopio Nazionale Galileo (GIANO, HARPS-N)
TNG = Telescope(
    name='Telescopio Nazionale Galileo (TNG)',
    location=coord.EarthLocation.of_site('Roque de los Muchachos'),
    instruments=[instruments.GIANO, instruments.HARPS_N],
    diameter=3.58 * u.m
)

# Mauna Kea Observatory, Gemini North Telescope (MAROON-X)
GeminiNorth = Telescope(
    name='Gemini North Telescope',
    location=coord.EarthLocation.of_site('gemini_north'),
    instruments=[instruments.MAROON_X],
    diameter=8.1 * u.m
)

# Mauna Kea Observatory, Canada-France-Hawaii Telescope (SPIRou)
CFHT = Telescope(
    name='Canada-France-Hawaii Telescope (CFHT)',
    location=coord.EarthLocation.of_site('Canada-France-Hawaii Telescope'),
    instruments=[instruments.SPIROU],
    diameter=3.6 * u.m
)

# Haute-Provence Observatory (SOPHIE)
HauteProvence = Telescope(
    name='Haute-Provence Observatory',
    location=coord.EarthLocation.of_site('Observatoire de Haute Provence'),
    instruments=[instruments.SOPHIE],
    diameter=1.93 * u.m
)


def print_all_telescopes(instruments: str = 'all') -> None:
    if instruments != 'all':
        raise NotImplementedError
    
    telescopes_list = [telescope for telescope in globals().values() if isinstance(telescope, Telescope)]
    telescopes_list.sort(reverse= True, key=lambda telescopes_list: telescopes_list.diameter.to(u.m).value)
    
    for telescope in telescopes_list:
        logger.print(f"{telescope.name} | {telescope.diameter} | Operational: {str(telescope.operational)}")
        logger.print(f"    {[instrument.name for instrument in telescope.instruments]}")
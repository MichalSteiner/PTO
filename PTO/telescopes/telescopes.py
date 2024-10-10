import observatories
import instruments
import astropy.coordinates as coord
import astropy.units as u
from dataclasses import dataclass, field


@dataclass
class Telescope():
    name: str
    location: coord.earth.EarthLocation
    instruments: list = field(default_factory=list)
    diameter: u.Quantity = 0 *u.m

CalarAlto = Telescope(
    name= 'Calar Alto Observatory',
    location= coord.EarthLocation.of_site('paranal'),
    instruments= [instruments.CARMENES],
    diameter= 3.5 * u.m
)

# VLT = Telescope(
#     name= 'VLT',
#     location= coord.EarthLocation.of_site('paranal'),
#     instruments= [],
#     diameter= 8.2 * u.m
# )

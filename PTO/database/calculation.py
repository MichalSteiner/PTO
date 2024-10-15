import numpy as np
from dataclasses import dataclass
import astropy.units as u
import astropy.nddata as nddata
import astropy.constants as con
from sympy import symbols, Eq, solve, sin, pi, sqrt, asin, cos, acos
import logging
from ..utils.utilities import logger_default

logger = logging.getLogger(__name__)
logger = logger_default(logger) 


class CalculationUtilities():

    def __drop_errors(self,
                      key: str) -> None:
        condition = (self.table[key].notna() & (
            self.table[f'{key}.Error.Lower'].isna() |
            self.table[f'{key}.Error.Upper'].isna()
        ))
        indices = self.table[condition].index

        self.table.loc[indices, key] = np.nan
        self.table.loc[indices, f'{key}.Error.Lower'] = np.nan
        self.table.loc[indices, f'{key}.Error.Upper'] = np.nan

        return None

    def __replace_errors(self,
                         key: str) -> None:

        condition = (self.table[key].notna() & (
            self.table[f'{key}.Error.Lower'].isna() |
            self.table[f'{key}.Error.Upper'].isna()
        ))
        indices = self.table[condition].index

        self.table.loc[indices, f'{key}.Error.Lower'] = 0
        self.table.loc[indices, f'{key}.Error.Upper'] = 0

        return None

    def _drop_without_errors(self,
                             mode: str = 'drop'):

        keys_with_erros = [
            key for key in self.table.keys() if (
                f'{key}.Error.Lower' in self.table.keys() and
                f'{key}.Error.Upper' in self.table.keys() and
                not (f'{key}'.startswith('Position.'))
            )
        ]

        match mode:
            case 'drop':
                for key in keys_with_erros:
                    self.__drop_errors(key)
            case 'replace':
                for key in keys_with_erros:
                    self.__replace_errors(key)
            case _:
                raise ValueError('Invalid mode')

    def _unit_conversion(self,
                         key_original: str,
                         unit_original: u.Unit,
                         key_result: str,
                         unit_result: u.Unit
                         ) -> None:

        condition = (self.table[key_original].notna() &
                     self.table[key_result].isna()
                     )
        indices = self.table[condition].index

        if len(indices) == 0:
            return None

        ArrayValues = nddata.NDDataArray(data=self.table.loc[indices, key_original].to_numpy(),
                                         uncertainty=nddata.StdDevUncertainty(np.max(
                                             [self.table.loc[indices, f'{key_original}.Error.Lower'].to_numpy(),
                                              self.table.loc[indices, f'{key_original}.Error.Upper'].to_numpy(
                                             ),
                                             ])),
                                         unit=unit_original
                                         )
        ConvertedArray = ArrayValues.convert_unit_to(unit_result)

        self.table.loc[indices, f'{key_result}'] = ConvertedArray.data
        self.table.loc[indices,
                       f'{key_result}.Error.Lower'] = ConvertedArray.uncertainty.array
        self.table.loc[indices,
                       f'{key_result}.Error.Upper'] = ConvertedArray.uncertainty.array

        return None

    def _add_Earth_and_Jupiter_units(self) -> None:
        self._unit_conversion(key_original='Planet.RadiusEarth',
                              unit_original=u.R_earth,
                              key_result='Planet.RadiusJupiter',
                              unit_result=u.R_jupiter)

        self._unit_conversion(key_original='Planet.RadiusJupiter',
                              unit_original=u.R_jupiter,
                              key_result='Planet.RadiusEarth',
                              unit_result=u.R_earth)

        self._unit_conversion(key_original='Planet.MassEarth',
                              unit_original=u.M_earth,
                              key_result='Planet.MassJupiter',
                              unit_result=u.M_jupiter)

        self._unit_conversion(key_original='Planet.MassJupiter',
                              unit_original=u.R_jupiter,
                              key_result='Planet.MassEarth',
                              unit_result=u.R_earth)

        self._unit_conversion(key_original='Planet.MinimumMassEarth',
                              unit_original=u.M_earth,
                              key_result='Planet.MinimumMassJupiter',
                              unit_result=u.M_jupiter)

        self._unit_conversion(key_original='Planet.MinimumMassJupiter',
                              unit_original=u.R_jupiter,
                              key_result='Planet.MinimumMassEarth',
                              unit_result=u.R_earth)

        self._unit_conversion(key_original='Planet.BestMassEstimateEarth',
                              unit_original=u.M_earth,
                              key_result='Planet.BestMassEstimateJupiter',
                              unit_result=u.M_jupiter)

        self._unit_conversion(key_original='Planet.BestMassEstimateJupiter',
                              unit_original=u.R_jupiter,
                              key_result='Planet.BestMassEstimateEarth',
                              unit_result=u.R_earth)
        return None

    def _solve_equation(self,
                        Equation: Eq,
                        symbols: symbols,
                        MAPPER: dict,
                        UNIT_MAPPER: dict,
                        ):
        
        for ind, row in self.table.iterrows():
            values = {key_name: key_value * UNIT_MAPPER[key_name] for key_name,key_value in MAPPER.items()}
            
            match sum([np.isnan(value) for value in values.values()]):
                case 1:
                    missing_variable = [name for name, value in values.items() if np.isnan(value)][0]
                    missing_variable = symbols(missing_variable)
                    logger.info('Hello')
                case 0:
                    continue
                case _ if sum([np.isnan(value) for value in values.values()]) > 1:
                    continue
                case _:
                    raise ValueError('An error that should have never appeared has appeared.')
            
            Equation_subs = Equation.subs(
                {symbols(var_name): var_value for var_name, var_value in values.items() if not(np.isnan(var_value))}
                )
            solution = solve(Equation_subs, missing_variable)
        
        
        
        ...
        


    def _calculate_impact_parameter(self):
        
        b, a, i, R_s = symbols('b, a i R_s')
        
        Equation = Eq(b,
                a*cos(i)/R_s
                )
        
        MAPPER= {
            'b': 'Planet.ImpactParameter',
            'a': 'Planet.SemiMajorAxis',
            'i': 'Planet.Inclination',
            'R_s': 'Star.Radius'
            }
        
        UNIT_MAPPER = {
            'b': 1,
            'a': u.au.to(u.m),
            'i': u.deg.to(u.rad),
            'R_s': u.R_s.to(u.R_s)
        }
        
        self._solve_equation(
            Equation=Equation,
            symbols=(b, a, i, R_s),
            MAPPER= MAPPER,
            UNIT_MAPPER=UNIT_MAPPER,
        )
        
        ...
        

    def _calculate_transit_length(self):
        
        T_14, P, R_s, R_p, a, b, i, e, omega = symbols('T_14, P R_s R_p a b i e omega')
        
        # sigma_P, sigma_R_s, sigma_R_p, sigma_a, sigma_b, sigma_i = symbols('sigma_P sigma_R_s sigma_R_p sigma_a sigma_b sigma_i')

        Equation = Eq(T_14,
                      (P / pi) * asin(
                          (R_s / a) * sqrt((1 + (R_p / R_s))**2 - b**2) / sin(i*pi/180)
                          ) #* (sqrt(1-e**2) / ((1 + e)*sin(omega)))
        )
        
        
        for ind, row in self.table.iterrows():
            values = {
                'T_14': row['Planet.TransitDuration'] * u.h.to(u.s),
                'P': row['Planet.Period'] * u.d.to(u.s),
                'R_s': row['Star.Radius'] * u.R_sun.to(u.m),
                'R_p': row['Planet.RadiusJupiter'] * u.R_jup.to(u.m),
                'a': row['Planet.SemiMajorAxis'] * u.au.to(u.m),
                'b': row['Planet.ImpactParameter'],
                'i': row['Planet.Inclination'],
                # 'e': row['Planet.Eccentricity'],
                # 'omega': row['Planet.ArgumentOfPeriastron']
            }
            
            match sum([np.isnan(value) for value in values.values()]):
                case 1:
                    missing_variable = [name for name, value in values.items() if np.isnan(value)][0]
                    missing_variable = symbols(missing_variable)
                    logger.info('Hello')
                case 0:
                    continue
                case _ if sum([np.isnan(value) for value in values.values()]) > 1:
                    continue
                case _:
                    raise ValueError('An error that should have never appeared has appeared.')
        

            Equation_subs = Equation.subs(
                {symbols(var_name): var_value for var_name, var_value in values.items() if not(np.isnan(var_value))}
                )
            solution = solve(Equation_subs, missing_variable)
        
        
        
        
        
        
        # dT_dP = T_14_expr.diff(P)
        # dT_dR_s = T_14_expr.diff(R_s)
        # dT_dR_p = T_14_expr.diff(R_p)
        # dT_da = T_14_expr.diff(a)
        # dT_db = T_14_expr.diff(b)
        # dT_di = T_14_expr.diff(i)

        # # General error propagation formula for sigma_T_14
        # sigma_T_14_expr = sqrt(
        #     (dT_dP * sigma_P)**2 +
        #     (dT_dR_s * sigma_R_s)**2 +
        #     (dT_dR_p * sigma_R_p)**2 +
        #     (dT_da * sigma_a)**2 +
        #     (dT_db * sigma_b)**2 +
        #     (dT_di * sigma_i)**2
        # )
        return None
    
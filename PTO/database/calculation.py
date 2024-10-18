import numpy as np
from dataclasses import dataclass
import astropy.units as u
import astropy.nddata as nddata
import astropy.constants as con
from sympy import symbols, Eq, solve, sin, pi, sqrt, asin, cos, acos
import logging
import sympy as smp

from ..utils.utilities import logger_default, time_function

logger = logging.getLogger(__name__)
logger = logger_default(logger) 


class CalculationUtilities():

    def __drop_errors(self,
                      key: str) -> None:
        """
        Drops values where the key has no errorbar.

        Parameters
        ----------
        key : str
            Key which to check for NaNs to drop.
        """
        condition = (self.table[key].notna() & (
            self.table[f'{key}.Error.Lower'].isna() |
            self.table[f'{key}.Error.Upper'].isna()
        ))
        indices = self.table[condition].index

        self.table.loc[indices, key] = np.nan
        self.table.loc[indices, f'{key}.Error.Lower'] = np.nan
        self.table.loc[indices, f'{key}.Error.Upper'] = np.nan

    def __replace_errors(self,
                         key: str) -> None:
        """
        Replace NaN errors with 0 for given key.

        Parameters
        ----------
        key : str
            Key which to check for NaNs to replace.
        """
        condition = (self.table[key].notna() & (
            self.table[f'{key}.Error.Lower'].isna() |
            self.table[f'{key}.Error.Upper'].isna()
        ))
        indices = self.table[condition].index

        self.table.loc[indices, f'{key}.Error.Lower'] = 0
        self.table.loc[indices, f'{key}.Error.Upper'] = 0

    def _handle_keys_without_errors(self,
                                    mode: str = 'drop') -> None:
        """
        Handles keys which don't hold error values.

        Parameters
        ----------
        mode : str, optional
            How to deal with NaNs errors, by default 'drop'. If 'drop', the values are dropped, if "replace", the errorbars are replaced with 0.

        Raises
        ------
        ValueError
            Raising error when invalid mode is provided.
        """
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
        """
        Converts units between x and y units for given keys, if not present.

        Parameters
        ----------
        key_original : str
            Original key to look to.
        unit_original : u.Unit
            Unit of the original key.
        key_result : str
            Resulting key to pass the converted value.
        unit_result : u.Unit
            Resulting unit to pass into the converted value
        """

        condition = (self.table[key_original].notna() &
                     self.table[key_result].isna()
                     )
        indices = self.table[condition].index

        if len(indices) == 0:
            return

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


    def _add_Earth_and_Jupiter_units(self) -> None:
        """
        Handles keys which are related by simple conversion, like units between Earth and Jupiter.
        """
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
    
    def __filter_degree_solution(self,
                                 solution: list,
                                 missing_variable: str,
                                 UNIT_MAPPER: dict) -> list:
        """
        Filters solutions for degree, when multiple solutions are given.
        
        The boundary is defined by np.pi, or 180 degrees.

        Parameters
        ----------
        solution : list
            Solution with more than single result. Typical for inclination values.
        missing_variable : str
            The name of the missing variable.
        UNIT_MAPPER : dict
            Mapper for units. This ensures correct units are used

        Returns
        -------
        solution : list
            Filtered solution to just single value.

        Raises
        ------
        ValueError
            For invalid units this raises error. Shouldn't happen.
        """
        
        if UNIT_MAPPER[missing_variable].unit == u.rad:
            return [value for value in solution if value < np.pi]
        elif UNIT_MAPPER[missing_variable].unit == u.deg:
            return [value for value in solution if value < 180]
        else:
            raise ValueError('Invalid unit.')
    
    def _filter_solutions(self,
                          solution: list,
                          missing_variable: str,
                          UNIT_MAPPER: dict) -> list:
        """
        Filter solutions with multiple results.

        Parameters
        ----------
        solution : list
            Solution to filter results in.
        missing_variable : str
            Which variable has been calculated
        UNIT_MAPPER : dict
            Unit mapper to double-check units

        Returns
        -------
        solution : list
            Filtered list with single solution.

        Raises
        ------
        NotImplementedError
            Not all variables are implemented by default. If triggered, it needs to be added to the match case syntax.
        """
        match missing_variable:
            case 'i':
                return self.__filter_degree_solution(solution, missing_variable, UNIT_MAPPER)
            case _:
                raise NotImplementedError('This variable is not implemented. FIXME')
    
    def _solve_equation(self,
                        Equation: Eq,
                        MAPPER: dict,
                        UNIT_MAPPER: dict,
                        ) -> None:
        """
        Solves the equation given a mapper and unit mapper.

        Parameters
        ----------
        Equation : Eq
            Equation to solve for.
        MAPPER : dict
            Mapper matching the symbols and the column names in the table.
        UNIT_MAPPER : dict
            Unit mapper to match symbols and the units for the symbol.

        Raises
        ------
        ValueError
            If multiple solutions detected but unhandled, check what happened.
        """
        for variable in MAPPER:
            other_variables = {item:MAPPER[item] for item in MAPPER if item !=variable}
            condition = (self.table[MAPPER[variable]].isna() &
                         [self.table[MAPPER[variable_other]].notna() for variable_other in other_variables]
                         )
            indices = self.table[condition].index

        if len(indices) == 0:
            return
            test = self.table[].iterrows()
        
        
        
        for ind, row in self.table.iterrows():
            values = {key_name: row[key_value] * UNIT_MAPPER[key_name].value for key_name,key_value in MAPPER.items()}
            
            match sum([np.isnan(value) for value in values.values()]):
                case 1:
                    missing_variable = [name for name, value in values.items() if np.isnan(value)][0]
                    other_variables = [name for name, value in values.items() if not(np.isnan(value))]
                    other_uncertainties = {
                        key: np.max([row[f"{MAPPER[key]}.Error.Lower"],
                                    row[f"{MAPPER[key]}.Error.Upper"]]) for key in other_variables
                    }
                    missing_variable = symbols(missing_variable)
                case 0:
                    continue
                case _ if sum([np.isnan(value) for value in values.values()]) > 1:
                    continue
            
            solution = solve(Equation, missing_variable)[0]
            
            subbed_solution = solution.subs(
                {symbols(var_name): var_value for var_name, var_value in values.items() if not(np.isnan(var_value))}
                )
            
            uncertainties = {
                variable: np.max([row[f'{MAPPER[variable]}.Error.Lower'], row[f'{MAPPER[variable]}.Error.Upper']]
                                 ) * UNIT_MAPPER[variable].value for variable in other_variables
            }
            
            total_uncertainty = 0
            for var in other_variables:
                partial_derivative = smp.diff(solution, symbols(var))
                total_uncertainty += ((partial_derivative * uncertainties[var])**2)
            total_uncertainty = smp.sqrt(total_uncertainty)
            
            uncertainty_result = total_uncertainty.subs(
                {symbols(var_name): var_value for var_name, var_value in values.items() if not(np.isnan(var_value))}
                )
            uncertainty_result = uncertainty_result / UNIT_MAPPER[str(missing_variable)].value
            
            if subbed_solution.is_finite:
                self.table.loc[ind, MAPPER[str(missing_variable)]] = float(subbed_solution)
                self.table.loc[ind, f"{MAPPER[str(missing_variable)]}.Error.Lower"] = float(uncertainty_result)
                self.table.loc[ind, f"{MAPPER[str(missing_variable)]}.Error.Upper"] = float(uncertainty_result)
        return

    def _calculate_impact_parameter(self) -> None:
        """
        Calculates the impact parameter related parameters.
        """
        
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
            'b': 1 * u.dimensionless_unscaled,
            'a': (1 * u.au).to(u.m),
            'i': (1 * u.deg).to(u.rad),
            'R_s': (1 * u.R_sun).to(u.m)
        }
        
        import time
        t1= time.time()
        self._solve_equation(
            Equation=Equation,
            MAPPER= MAPPER,
            UNIT_MAPPER=UNIT_MAPPER,
        )
        print(time.time() -t1)
        
    def __check_eccentricity(self) -> None:
        """
        Checks for eccentricity values.
        
        If a value of eccentricity is NaN or 0, it will be converted to 0. This ensures the equation for T14 is unaffected.
        Furthermore, Argument of periastron (omega) is set to 90 degrees, which gives 1 when used in 
        """
        
        condition = (self.table['Planet.Eccentricity'].isna() |
                     self.table['Planet.Eccentricity'] == 0
                     )
        indices = self.table[condition].index

        if len(indices) == 0:
            return
        
        self.table.loc[indices, 'Planet.Eccentricity'] = 0
        self.table.loc[indices, 'Planet.Eccentricity.Error.Lower'] = 0
        self.table.loc[indices, 'Planet.Eccentricity.Error.Upper'] = 0
        self.table.loc[indices, 'Planet.ArgumentOfPeriastron'] = 90
        self.table.loc[indices, 'Planet.ArgumentOfPeriastron.Error.Lower'] = 0
        self.table.loc[indices, 'Planet.ArgumentOfPeriastron.Error.Upper'] = 0
        return 
    
    
    def _calculate_transit_length(self) -> None:
        """
        Calculates transit length, if not provided.
        
        It also reformats the values for eccentricity/ argument of periastron to ensure it will not impact the equation.
        """
        
        self.__check_eccentricity()
        
        T_14, P, R_s, R_p, a, b, i, e, omega = symbols('T_14, P R_s R_p a b i e omega')
        
        Equation = Eq(T_14,
                      (P / pi) * asin(
                          (R_s / a) * sqrt((1 + (R_p / R_s))**2 - b**2) / sin(i*pi/180)
                          ) * (sqrt(1-e**2) / ((1 + e)*sin(omega)))
        )
        
        MAPPER= {
            'T_14': 'Planet.TransitDuration',
            'P': 'Planet.Period',
            'R_s': 'Star.Radius',
            'R_p': 'Planet.RadiusJupiter',
            'a': 'Planet.SemiMajorAxis',
            'b': 'Planet.ImpactParameter',
            'i': 'Planet.Inclination',
            'e': 'Planet.Eccentricity',
            'omega': 'Planet.ArgumentOfPeriastron'
            }
        
        UNIT_MAPPER = {
            'T_14': (1*u.hour).to(u.s),
            'P': (1*u.day).to(u.s),
            'R_s': (1 * u.R_sun).to(u.m),
            'R_p': (1 * u.R_jup).to(u.m),
            'a': (1 * u.au).to(u.m),
            'b': 1 * u.dimensionless_unscaled,
            'i': (1 * u.deg).to(u.rad),
            'e': 1 * u.dimensionless_unscaled,
            'omega': (1 * u.deg).to(u.rad),
        }

        return
    
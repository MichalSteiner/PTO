# About PTO (Planner of Transit Observations):
`PTO` is a code that allows the user to calculate transit windows and plot their visibility for given location. Currently, main focus of `PTO` is handling high-resolution observations of transits, but further extensions can be made.

The main capabilities of `PTO` are:

1. Load NASA Exoplanet archive table using TAP service
2. Filter targets based on common characteristics
3. Calculate positions of transit windows
4. Calculate the visibility of each transit windows

Please read the documentation at [pto.readthedocs.io](pto.readthedocs.io).

For any bug/ feature requests, please open an issue on the GitLab page [https://gitlab.unige.ch/spice_dune/pto](https://gitlab.unige.ch/spice_dune/pto).

To get started, please read the notebook at [https://pto.readthedocs.io/en/latest/Get_started.html](https://pto.readthedocs.io/en/latest/Get_started.html).
The tutorial will be merged in a single example notebook soon.

Acknowledgements: This code is redesigned by Michal Steiner based on legacy transit planner written by: Romain Allart, Mara Attia, Daniel Bayliss, Vincent Bourrier.
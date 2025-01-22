Planner for Transit Observations documentation
==============================================

Welcome to the Planner for Transit Observations (PTO) documentation. PTO is a code for calculation of transit windows, in particular with the intent of observing from ground based telescopes. It has fully implemented API download to NASA Exoplanet Archive, getting updated system parameters for all exoplanets.

This version of the code is in developement.

The main steps are as follow:

1. Selection of database and loading it

2. Filtering down the sample based on criteria

3. Calculation of transit windows

4. Plotting observability of calculated windows based on location

In case of feedback, please contact Michal Steiner (Michal.Steiner@unige.ch) or raise an issue on the Github page (https://github.com/MichalSteiner/PTO).

Furthermore, the code will allow user to simulate a spectral dataset from given spectrograph, allowing to calculate the signal to noise ratio of the observation. This will be implemented in the future.

Setup

.. toctree::
   :maxdepth: 1
   :caption: Setup:

   Installation

Notebooks:

.. toctree::
   :maxdepth: 1
   :caption: Notebooks:

   Get_started
   Catalogs
   Filtering
   Custom_ephemeris
   Telescopes

API:

.. toctree::
   :maxdepth: 2
   :caption: PTO package API:

   modules

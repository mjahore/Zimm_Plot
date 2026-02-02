# Zimm Plot: Analyzes light scattering measurements with a Zimm plot to calculate Mw, Rg, and A2.
###### Written and Developed by Michael J. A. Hore (hore@case.edu), Case Western Reserve University

## Description
This python script reads in light scattering measurements and constructs a Zimm plot. The input data should be a .csv file containing concentration, angle, and the Rayleigh ratio (in that order). Assumed units for these values are (mg/mL), degrees, and (cm^-1), respectively. The code will perform the necessary conversions. Configuration of wavelength, refractive index, and (dn/dc) are within the python script. It is assumed that the data is tabulated as blocks of varying concentration with increasing angle with each subsequent block. See `example_data.dat` for an example, using data from Hiemenz and Lodge (2nd Edition).

## Execution
Run with `./zimm_plot.py [data file] [number of concentrations] [number of angles]`.

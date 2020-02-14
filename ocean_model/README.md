# Ocean Model
The standard version of Speedy has a vertical slab ocean model, whereby the water columns cannot interact laterally with each other. The changes to the model in this folder make it possible for heat to diffuse horizontally through the ocean from the equator to the pole. This means that surplus incoming radiation can be transported to the poles by both the atmosphere and the ocean.


These changes were made by Frank Selten of the KNMI (www.knmi.nl) for a lecture series ‘Water in the Atmosphere’, in a collaboration with Pier Sibesma and Stephan de Roode from the Delft University of Technology (www.tudelft.nl).


To use this model, copy the content of this folder to the ‘update’ folder before running a calculation. The run_exp.s script will then use the contents of the ‘update’ folder for run.


These changes may be used freely and are provided ‘as is’, without warranty of any kind.

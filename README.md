# ExoSPEEDY
## The ICTP atmospheric general circulation model SPEEDY, adapted for exoplanets

The software in this repository is adapted from the ICTP AGCM software, nicknamed SPEEDY, see:
https://www.ictp.it/research/esp/models/speedy.aspx

The SPEEDY software was developed at ICTP by Franco Molteni and Fred Kucharski and is described in the publications referenced here:
http://users.ictp.it/~kucharsk/speedy-net.html

Please cite these publications when publishing your own work based on this software.

ExoSPEEDY adapts the SPEEDY software for use in climate modeling of exoplanets.
This work was done at Delft University of Technology under the guidance of Daphne Stam.

The first step in adapting the code for exoplanets is to identify the flow of earth parameters.
To this purpose, each subroutine is annotated with comment lines starting with the string 'C--IO x' where x is one of:
- r - reads from a file
- w - writes to a file
- h - gets input from a .h file
- m - sets or modifies value in a .h file
- s - sets a value in the subroutine

# GalBar

Potential with a Galactic bar as described in Thomas et al. 2023 (xxxxxxxxxxxxxxxxx).

## Updated March 18th 2023 by Guillaume Thomas *guillaume.thomas.astro .at. gmail.com*
-------------------------------------------------------------------------------------------------------------------------------



### REQUIREMENTS

You need Nemo (https://teuben.github.io/nemo/) installed on your computer and activated when you install it.



### INSTALLATION
To do the installation, just run:

> ./install.sh

The process can take several minute as it has to recompile gyrfalcON



### EXAMPLE
An example of the type of potfile used by GalBar and of the file where are defined the amplitude is presented in the example repository. The (background) axisymmetric potential is similar to the one used in Thomas et al. 2023, and the bar included in that file has a pattern speed of 39 km/s/kpc and an angle with the solar axis (alpha_0) of 28 degrees.

To call GalBar in GyrfalcON and in compatible tools, you should write accname=GalBar.

Note here that to use the mkorbit tools from NEMO with the GalBar potential, you should enter potname=GalBar, and you should imperativelly use the leapfrog mode.

# celes-ellipsoid-opt
 An extension of an early version of the CELES simulation toolbox to be compatible with ellipsoidal scatterers. In addition, includes adjoint optimization-based routines for spherical and ellipsoidal particles. Ellipsoidal T-matrices are calculated using the extended boundary condition method and are not reliable for very high aspect ratios or for very densely packed scatterers.
 
 Two examples can be found, celes-forward.m, and celes-optimize.m which are small examples showing how the code works to simulate and optimize ellipsoid scatterer configurations.

The original CELES codebase that this code is heavily based off can be found here:
https://disordered-photonics.github.io/celes/

Work on this project was done by Alan Zhan, Maksym Zhelyeznyakov, Taylor Fryett, and Shane Colburn in Professor Arka Majumdar's group at the University of Washington, Seattle.

Associated Publications:

https://www.osapublishing.org/ao/abstract.cfm?uri=ao-57-6-1437

https://advances.sciencemag.org/content/5/10/eaax4769

https://www.osapublishing.org/osac/abstract.cfm?uri=osac-3-1-89

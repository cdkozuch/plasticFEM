# plasticFEM
Finite Element Method (FEM) for plastic material implemented in MATLAB (it is also compatible with Octave 7.2). 

## Note
The current version only accepts applied displacements (i.e. no applied forces).

## References
* Kim, N. H. (2014). Introduction to nonlinear finite element analysis. Springer Science & Business Media.
* Simo, J. C., & Hughes, T. J. (2006). Computational inelasticity (Vol. 7). Springer Science & Business Media.

## Description
The code is composed of seven files: five functions and two scripts. The functions are: bmatHex8.m, globElemInd.m, globTanStiff.m, matStiffTen3D.m, and returnMap.m. Full descriptions of the purpose, input, and output of each function code file is given at the top of each file. The scripts are runFEM.m and runReturnMap.m. Brief descriptions of each file are given below:

* bmatHex8.m: computes the strain-displacement matrix for 8-node hexahedral
* globElemInd.m: determines the global degrees of freedom for an element
* globTanStiff.m: computes the global tangent stiffness matrix, global internal force vector, strain, stress, and history variables
* matStiffTen3D.m: computes 3D material stiffness tensor for isotropic material
* returnMap.m: performs return mapping
* runFEM.m: perform finite element analysis
* runReturnMap.m: performs return mapping iteratively

The code relies on two data files, nodes.dat and elements.dat, which contain the nodal coordinates and element connectivity, respectively.

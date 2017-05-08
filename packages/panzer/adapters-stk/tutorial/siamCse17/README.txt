The files in this directory solve the following problem:

  -\Delta u(x, y) + (1 - 8 \pi^2) u(x, y) = \sin(2 \pi x) \sin(2 \pi y),
    (x,y) \in \Omega = (0,1) \times (0,1),
  u(x, y) = 0, (x,y) \in \partial\Omega.

This example was originally presented at SIAM CSE 2017.  The corresponding
presentation can be found through a link at
https://trilinos.org/packages/panzer.

The following is an explanation of all the files in this directory.

CMakeLists.txt:
  This file tell CMake how to build this example with the rest of Trilinos.

input.xml:
  An input file in which we specify the mesh, physics blocks, closure models,
  boundary conditions, and solver options.

main.cpp:
  This runs the example.  There's an awful lot here that new users aren't
  intended to understand -- instead focus on the files below.

myBCStrategy*:
  All the files needed to implement the zero Dirichlet boundary condition.

myClosureModelFactory*, mySourceTerm*:
  All the files needed to account for the source term,
  \sin(2 \pi x) \sin(2 \pi y).

myEquationSet*:
  All the files needed to construct our equation set, or rather, the residual
  corresponding to our partial differential equation above.

README.txt:
  This file, obviously.

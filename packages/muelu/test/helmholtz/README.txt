Tests for MueLu::ShiftedLaplacian and Galeri-Helmholtz problems

These tests exercise the functions in Galeri which produce finite difference and
finite element discretizations of the Helmholtz equation on uniform quad/hex meshes.
The model problem is an infinite domain homogeneous medium; the boundary conditions
used to simulate the infinite domain are perfectly matched layers (on all sides).
The matrices and right hand sides are fed to MueLu::ShiftedLaplacian, which solves
the systems using MueLu as a preconditioner to Belos. The algebraic multigrid is set
up in a few different ways, and for each method, certain matrices need to be set.

In the following, it is assumed that for a stiffness matrix K, mass matrix M,
and angular frequency omega, the Helmholtz operator is A = K - \omega^2 M, and
the Shifted Laplacian operator is S = K - (1+ i \beta) (\omega^2) M.



1. For higher frequency problems, we suggest using Shifted Laplacians where the prolongation
and restriction matrices are produced using the stiffness matrix; The reason the stiffness
matrix is used is because it produces prolongation/restriction matrices for the Laplacian,
i.e. the smoothed-aggregation basis functions look similar to geometric multigrid. This method
can be done with either setupFastRAP() or setupSlowRAP(). In setupFastRAP(), the complex shift
added to the Helmholtz operator for damping is the same at every level of the multigrid
hierarchy, while in setupSlowRAP(), the complex shift can be varied. The initialization is the
same for both methods:

SLSolver = Teuchos::rcp( new MueLu::ShiftedLaplacian() );
SLSolver -> setstiff(...);
SLSolver -> setParameters(...);
SLSolver -> initialize();

This will construct the prolongation/restriction matrices with the stiffness matrix
(which you put in through the setstiff() function). Now, once the Helmholtz operator
and Shifted Laplacian operator are constructed, do the following for setupFastRAP():

SLSolver -> setProblemMatrix(...);
SLSolver -> setPreconditioningMatrix(...);
SLSolver -> setupFastRAP();

This will setup the level smoothers for multigrid, and also setup the Belos iterative
solver with MueLu as the preconditioner. For setupSlowRAP(), the mass matrix and level
complex shifts need to be defined; thus, the calls for setupSlowRAP() should go

SLSolver -> setProblemMatrix(...);
SLSolver -> setPreconditioningMatrix(...);
SLSolver -> setmass(...);
SLSolver -> setLevelShifts(...);
SLSolver -> setupSlowRAP();



2. For lower frequency problems, regular smoothed-aggregation AMG for the Shifted Laplacian
should suffice. This can be done with setupNormalRAP(). The calls for this mode are

SLSolver = Teuchos::rcp( new MueLu::ShiftedLaplacian() );
SLSolver -> setPreconditioningMatrix(Pmat);
SLSolver -> initialize();
SLSolver -> setProblemMatrix(Amat);
SLSolver -> setupNormalRAP();



3. Once you've called setupSlowRAP(), setupFastRAP(), or setupNormalRAP(), you can solve for
any RHS multivector using the solve function:

SLSolver -> solve(B,X);

Where B is the right hand side vector, and X is the solution vector.



4. For multiple frequencies, a sequence of linear systems needs to be solved. In this case,
setupFastRAP() or setupSlowRAP() modes are best, since they keep the prolongation and restriction
matrices around; this is because the stiffness matrix does not change even if the frequency
changes. Thus, you should only call the initialize() function once for a multiple frequency
problem; the calls to the Shifted Laplacian solver and loop over the frequencies should look
like this:

SLSolver = Teuchos::rcp( new MueLu::ShiftedLaplacian() );
SLSolver -> setstiff(...);
SLSolver -> setParameters(...);
SLSolver -> initialize();

for( int i=0; i < num_omegas; i++) {

  ...(construct Helmholtz and Shifted Laplacian matrices)...
  SLSolver -> setProblemMatrix(...);
  SLSolver -> setPreconditioningMatrix(...);
  SLSolver -> setupFastRAP();
  ...(compute right hand side)...  
  SLSolver -> solve(B,X);

}



5. A summary of the tests:
Helmholtz1D.cpp    -> finite differences in 1D, regular AMG, one frequency
Helmholtz2D.cpp    -> finite differences in 2D, setupSlowRAP(), solves for multiple frequencies
Helmholtz3D.cpp    -> finite differences in 3D, setupNormalRAP(), one frequency
HelmholtzFEM2D.cpp -> finite elements in 2D, setupFastRAP(), one frequency
HelmholtzFEM3D.cpp -> finite elements in 3D, setupNormalRAP(), one frequency

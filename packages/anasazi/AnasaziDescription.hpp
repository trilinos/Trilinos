/*!
   \mainpage %Anasazi: A Block Eigensolver Package

   \section intro Introduction

   %Anasazi is an extensible and interoperable framework for large-scale eigenvalue algorithms.
The motivation for this framework is to provide a generic interface to a collection of algorithms for 
solving large-scale eigenvalue problems of DOE interest.
Anasazi is interoperable because both the matrix and vectors (defining the
eigenspace) are considered to be opaque objects---only knowledge of the matrix and
vectors via elementary operations is necessary. An implementation of Anasazi
is accomplished via the use of interfaces. Current interfaces available include
epetra and so any libraries that understand epetra matrices and vectors (such
as AztecOO) may also be used in conjunction with Anasazi.

One of the goals of Anasazi is to allow the user the flexibility to specify the
data representation for the matrix and vectors and so leverage any existing software
investment. The algorithms that will be initially available through Anasazi are
block implicitly restarted Arnoldi and Lanczos methods and preconditioned eigensolvers.
These include a locally optimal block preconditioned congugate gradient iteration (LOBPCG) for 
symmetric positive definite generalized eigenvalue problems, and a restarted preconditioned eigensolver 
for nonsymmetric eigenvalue problems. 

   \section contributors Anasazi Contributors

   The following people have contributed to the development of %Anasazi:

   <ul>
 	<li> Rich Lehoucq, Sandia National Labs, rblehou@sandia.gov
	<li> Heidi Thornquist, Sandia National Labs, hkthorn@sandia.gov
   </ul>

*/

## GMRES Kokkos-based solver:

GMRES (the Generalized Minimum RESidual method) is an iterative solver for sparse linear systems Ax=b.  It can be used for symmetric or non-symmetric systems.  (For full details, see "Iterative Methods for Sparse Linear Systems," Saad, 2003.)

### Command-Line parameters for ex\_real\_A:
Current solver parameters for the real-valued example are as follows:

"**--filename**   :  The name of a matrix market (.mtx) file for matrix A (Default bcsstk09.mtx)."
"**--max-subsp**   :  The maximum size of the Kyrlov subspace before restarting (Default 50)."
"**--max-restarts:**  Maximum number of GMRES restarts (Default 50)."
"**--tol        :**  Convergence tolerance.  (Default 1e-10)."
"**--ortho       :**  Type of orthogonalization. Use 'CGS2' or 'MGS'. (Default 'CGS2')"
"**--rand\_rhs**    :  Generate a random right-hand side b.  (Without this option, the solver default generates b = vector of ones.)"

### Solver input parameters:
The gmres function takes the following input paramters:
**A:** A Kokkos::CrsMatrix for the linar system Ax=b.
**B:** A Kokkos::View that is the system right-hand side. Must have B.extent(1)=1. (Currently only one right-hand side is supported.)
**X:** A Kokkos::View that is used as both the initial vector for the GMRES iteration and the output for the solution vector.  (Must have X.extent(1)=1.)
**M:** A pointer to a KokkosSparse::Experimental::Preconditioner. Only right preconditioning is supported at this time.

### Handle input parameters:
The solver has a GMRESHandle struct to pass in solver options.  Available options are:
**tol:** The convergence tolerance for GMRES.  Based upon the relative residual. The solver will terminate when norm(b-Ax)/norm(b) <= tol. (Default: 1e-8)
**m:** The restart length (maximum subspace size) for GMRES.  (Default: 50)
**maxRestart:** The maximum number of restarts (or 'cycles') that GMRES is to perform. (Default: 50)
**ortho:** The orthogonalization type.  Can be "CGS2" (Default) or "MGS".  (Two iterations of Classical Gram-Schmidt, or one iteration of Modified Gram-Schmidt.)
**verbose:** Tells solve to print more information

### Solver Output:
The GMRESHandle struct is also used to pass back solver statistics. These include:
**numIters**: The number of iterations the GMRES solver took before terminating.
**endRelRes**: The ending relative residual norm attained in the GMRES solve.
**convFlagVal**: An enum FLAG value which will return "Conv" if the solver converged, "NoConv" if the solver did not converge, or "LOA" if the solver converged with respect to the internally computed residual norm but lost accuracy from the true residual norm.

### Test files:

You can download the matrices for the real example and the complex test from the [SuiteSparse](https://sparse.tamu.edu/) matrix collection.  The two matrices needed are:
* **young1c** (for complex-valued test)
* **bcsstk09** (for real-valued example)

The real-valued test uses a matrix generated directly by Kokkos Kernels.

### Test Measurements:
These measurements were taken on 7/23/21, running on an NVIDIA V100 GPU on Weaver7.
(Timings based upon the GMRES::TotalTime profiling region.)

**ex\_real\_A:** Converges in 2271 iterations and 0.9629 seconds.

(The two following timings total the time for the CGS2 and MGS tests.)
**test\_real\_A:** Converges in 30 iterations (with a restart size of 15) and 0.2536 seconds.

**test\_cmplx\_A:** Converges in 652 iterations (to a tolerance of 1e-5) in 2.822 seconds.

### Concerns, enhancements, or bug reporting:
If you wish to suggest an enhancement or make a bug report for this solver code, please post an issue at https://github.com/kokkos/kokkos-kernels/issues or email jloe@sandia.gov.

SAND2021-8676 O

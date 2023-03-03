## KokkosSparse Preconditioner Interface:

The `KokkosSparse_Preconditioner` class provides an abstract base class interface to use Kokkos-based preconditioners with iterative linear solvers. In particular, this class is designed to work with the Kokkos-based GMRES implementation in `examples/gmres`. It may also be useful for integrating Kokkos-based preconditioners into other solver codes and packages. (For Trilinos users: This class is loosely based upon the IfPack2::Preconditioner class.)  The Preconditioner class and derived classes sit in the `KokkosSparse::Experimental` namespace.

### Input parameters:

##### Preconditioner template paramters:
The KokkosSparse Preconditioner has the following template parameters:
**CRS:** The type of the compressed row sparse matrix. All key types should be derivable from this.

### Preconditioner Base Class Functions (See code for full details.):

**constructor**: Empty in the base class.
**apply**`( Kokkos::View<ScalarType*> &X, Kokkos::View<ScalarType*> &Y, const char transM[] = "N", ScalarType alpha = 1.0, ScalarType beta =  0):`
Returns `y = beta * y + alpha * M * x` where `M` is the preconditioning operator.  (May apply `M` transposed if `transM = 'T' or 'C'` for 'transpose' and 'conjugate transpose' respectively.)
**setParameters():** Used to set preconditioner parameters.
**initialize():** Set up the graph structure of this preconditioner. If the graph structure of the constructor's input matrix has changed, or if you have not yet called `initialize()`, you must call `initialize()` before you may call `compute()` or `apply()`. Thus, `initialize()` corresponds to the "symbolic factorization" step of a sparse factorization, whether or not the specific preconditioner actually does a sparse factorization.
**compute():** Set up the numerical values in this preconditioner. If the values of the constructor's input matrix have changed, or if you have not yet called `compute()`, you must call `compute()` before you may call `apply()`. Thus, `compute()` corresponds to the "numeric factorization" step of a sparse factorization, whether or not the specific preconditioner actually does a sparse factorization.
**isInitialized()** and **isComputed()** return whether the preconditioner has been initialized or computed, respectively.
**hasTransposeApply()** Returns true if the transposed (or conjugate transposed) operator apply is available for this preconditioner.  Base class function returns `false`.

### Derived Preconditioner Classes:

#### Matrix Prec:
A simple class that takes a `KokkosSparse::CrsMatrix` to apply as the preconditioner `M`.  This matrix is given to the preconditioner constructor function.  There are no parameters to set.  The functions `initialize` and `compute` do not perform any operations.

### Testing

##### Command-Line parameters for test\_prec:
Current solver parameters for the GMRES preconditioning test are as follows:

"**--mat-size**   :   The size of the n x n test matrix. (Default: n=1000.)
"**--max-subsp**   :  The maximum size of the Kyrlov subspace before restarting (Default 50)."
"**--max-restarts:**  Maximum number of GMRES restarts (Default 50)."
"**--tol        :**  Convergence tolerance.  (Default 1e-10)."
"**--ortho       :**  Type of orthogonalization. Use 'CGS2' or 'MGS'. (Default 'CGS2')"
"**--rand\_rhs**    :  Generate a random right-hand side b.  (Without this option, the solver default generates b = vector of ones.)"
"**--help**         : Display a help message."

#### Test Measurements:
**test\_prec:** Tests the GMRES solver with a diagonal matrix and uses its inverse as the preconditioner.  Should always converge in 1 iteration, regardless of problem size.

### Concerns, enhancements, or bug reporting:
If you wish to suggest an enhancement or make a bug report for this preconditioner code, please post an issue at https://github.com/kokkos/kokkos-kernels/issues or email jloe@sandia.gov.

SAND2021-9744 O


// File: index.xml

// File: classAztecOO.xml
%feature("docstring") AztecOO "

AztecOO: An object-oriented wrapper for Aztec.

Currently it accepts a Petra matrix, initial guess and RHS as separate
arguments, or alternatively, accepts a Epetra_LinearProblem. If
constructed using a Epetra_LinearProblem, AztecOO will infer some
solver/preconditioner, etc., options and parameters. Users may
override these choices and manually choose from among the full set of
Aztec options using the SetAztecOption() and SetAztecParam()
functions.

AztecOO will solve a linear systems of equations: $ AX=B $, using
Epetra objects and the Aztec solver library, where $A$ is an
Epetra_Operator or Epetra_RowMatrix (note that the Epetra_Operator
class is a base class for Epetra_RowMatrix so that Epetra_RowMatrix
isa Epetra_Operator.) $X$ and $B$ are Epetra_MultiVector objects.

WARNING:   AztecOO does not presently support solution of more than
one simultaneous right-hand-side.

C++ includes: AztecOO.h ";

/*  Constructors/destructors.  */

%feature("docstring")  AztecOO::AztecOO "AztecOO::AztecOO(Epetra_Operator *A, Epetra_MultiVector *X,
Epetra_MultiVector *B)

AztecOO Constructor.

Creates a AztecOO instance, passing in already-defined objects for the
linear operator (as an Epetra_Operator), left-hand-side and right-
hand-side.

Note: Use of this constructor may prohibit use of native AztecOO
preconditioners, since an Epetra_Operator is not necessarily an
Epetra_RowMatrix and all AztecOO incomplete factorization
preconditioners are based on having explicit access to matrix
coefficients. Polynomial preconditioners are available if the
Epetra_Operator passed in here has a non-trivial definition of the
NormInf() method and HasNormInf() returns true. ";

%feature("docstring")  AztecOO::AztecOO "AztecOO::AztecOO(Epetra_RowMatrix *A, Epetra_MultiVector *X,
Epetra_MultiVector *B)

AztecOO Constructor.

Creates a AztecOO instance, passing in already-defined objects for the
linear operator (as an Epetra_RowMatrix), left-hand-side and right-
hand-side.

Note: Use of this constructor allows full access to native AztecOO
preconditioners, using the Epetra_RowMatrix A passed in here as the
basis for computing the preconditioner. All AztecOO incomplete
factorization preconditioners are based on having explicit access to
matrix coefficients. Polynomial preconditioners are also available. It
is possible to change the matrix used for computing incomplete
factorization by calling the SetPrecMatrix() method. It is also
possible to provide a user-supplied preconditioner via
SetPrecOperator(). ";

%feature("docstring")  AztecOO::AztecOO "AztecOO::AztecOO(const
Epetra_LinearProblem &LinearProblem)

AztecOO Constructor.

Creates a AztecOO instance, using a Epetra_LinearProblem, passing in
an already-defined Epetra_LinearProblem object. The
Epetra_LinearProblem class is the preferred method for passing in the
linear problem to AztecOO because this class provides scaling
capabilities and self-consistency checks that are not available when
using other constructors.

Note: If the Epetra_LinearProblem passed in here has a non-trivial
pointer to an Epetra_Matrix then use of this constructor allows full
access to native AztecOO preconditioners, using the Epetra_RowMatrix A
passed in here as the basis for computing the preconditioner. All
AztecOO incomplete factorization preconditioners are based on having
explicit access to matrix coefficients. Polynomial preconditioners are
also available. It is possible to change the matrix used for computing
incomplete factorization by calling the SetPrecMatrix() method. It is
also possible to provide a user-supplied preconditioner by call
SetPrecOperator().

If the Epetra_LinearProblems passed in here has only an
Epetra_Operator, then use of this constructor may prohibit use of
native AztecOO preconditioners, since an Epetra_Operator is not
necessarily an Epetra_RowMatrix and all AztecOO incomplete
factorization preconditioners are based on having explicit access to
matrix coefficients. Polynomial preconditioners are available if the
Epetra_Operator passed in here has a non-trivial definition of the
NormInf() method and HasNormInf() returns true. ";

%feature("docstring")  AztecOO::AztecOO "AztecOO::AztecOO()

AztecOO Default constructor. ";

%feature("docstring")  AztecOO::AztecOO "AztecOO::AztecOO(const
AztecOO &Solver)

AztecOO Copy Constructor.

Makes copy of an existing AztecOO instance. ";

%feature("docstring")  AztecOO::~AztecOO "AztecOO::~AztecOO(void)

AztecOO Destructor.

Completely deletes a AztecOO object. ";

/*  Post-construction setup methods.  */

%feature("docstring")  AztecOO::SetProblem "int
AztecOO::SetProblem(const Epetra_LinearProblem &prob, bool
call_SetPrecMatrix=false)

AztecOO Epetra_LinearProblem Set.

Associates an already defined Epetra_LinearProblem as the problem that
will be solved during iterations. This method allows the user to
change which problem is being solved by an existing AztecOO object.

Internally calls SetUserMatrix() if the Epetra_LinearProblem's
operator can be cast to Epetra_RowMatrix, otherwise calls
SetUserOperator().

IMPORTANT WARNING *** This method calls SetUserMatrix(), which also
sets the *preconditioner* matrix to the matrix passed in, by
internally calling SetPrecMatrix(), but *ONLY* if SetPrecMatrix()
hasn't previously been called. If the user wants to make sure that any
pre-existing preconditioner is replaced, they must set the optional
bool argument 'call_SetPrecMatrix' to true, which will force this
function to call SetPrecMatrix().

WARNING:  The first time this method is called, the default options
and parameters are set. Therefore, this method should be called before
setting any individual options or parameters.

If a preconditioner has been pre-built and associated with this
AztecOO object, the Epetra_LinearProblem being passed in to this
method must have compatible domain and range maps. ";

%feature("docstring")  AztecOO::SetUserOperator "int
AztecOO::SetUserOperator(Epetra_Operator *UserOperator)

AztecOO User Operator Set.

Associates an already defined Epetra_Operator as the linear operator
for the linear system system that will be solved during iterations.
This set method allows the user to pass any type of linear operator to
AztecOO, as long as the operator implements the Epetra_Operator pure
virtual class, and has proper domain and range map dimensions.
Epetra_CrsMatrix and Epetra_VbrMatrix objects can be passed in through
this method. ";

%feature("docstring")  AztecOO::SetUserMatrix "int
AztecOO::SetUserMatrix(Epetra_RowMatrix *UserMatrix, bool
call_SetPrecMatrix=false)

AztecOO User Matrix Set.

Associates an already defined Epetra_Matrix as the matrix that will be
used by AztecOO as the linear operator when solving the linear system.
Epetra_CrsMatrix and Epetra_VbrMatrix objects can be passed in through
this method.

IMPORTANT WARNING *** This method sets the preconditioner matrix to
the matrix passed in here, by internally calling SetPrecMatrix(), but
*ONLY* if SetPrecMatrix() hasn't previously been called. If the user
wants to make sure that any pre-existing preconditioner is replaced,
they must set the optional bool argument 'call_SetPrecMatrix' to true,
which will force this function to call SetPrecMatrix(). ";

%feature("docstring")  AztecOO::SetLHS "int
AztecOO::SetLHS(Epetra_MultiVector *X)

AztecOO LHS Set.

Associates an already defined Epetra_MultiVector (or Epetra_Vector) as
the initial guess and location where the solution will be return. ";

%feature("docstring")  AztecOO::SetRHS "int
AztecOO::SetRHS(Epetra_MultiVector *B)

AztecOO RHS Set.

Associates an already defined Epetra_MultiVector (or Epetra_Vector) as
the right-hand-side of the linear system. ";

%feature("docstring")  AztecOO::UnsetLHSRHS "int
AztecOO::UnsetLHSRHS()

AztecOO unset LHS and RHS.

Sets to null the previously objects.. ";

%feature("docstring")  AztecOO::SetPrecMatrix "int
AztecOO::SetPrecMatrix(Epetra_RowMatrix *PrecMatrix)

AztecOO Preconditioner Matrix Set.

Associates an already defined Epetra_Matrix as the matrix that will be
used by AztecOO when constructing native AztecOO preconditioners. By
default, if AztecOO native preconditioners are used, the original
operator matrix will be used as the source for deriving the
preconditioner. However, there are instances where a user would like
to have the preconditioner be defined using a different matrix than
the original operator matrix. Another common situation is where the
user may not have the operator in matrix form but has a matrix that
approximates the operator and can be used as the basis for an
incomplete factorization. This set method allows the user to pass any
Epetra_RowMatrix to AztecOO for use in constructing an AztecOO native
preconditioner, as long as the matrix implements the Epetra_RowMatrix
pure virtual class, and has proper domain and range map dimensions.
Epetra_CrsMatrix and Epetra_VbrMatrix objects can be passed in through
this method. ";

%feature("docstring")  AztecOO::SetPrecOperator "int
AztecOO::SetPrecOperator(Epetra_Operator *PrecOperator)

AztecOO External Preconditioner Set.

Associates an already defined Epetra_Operator as the preconditioner
that will be called during iterations. This set method allows the user
to pass any type of preconditioner to AztecOO, as long as the
preconditioner implements the Epetra_Operator pure virtual class, and
has proper domain and range map dimensions. Ifpack preconditioners can
be passed in through this method. ";

%feature("docstring")  AztecOO::SetStatusTest "int
AztecOO::SetStatusTest(AztecOO_StatusTest *StatusTest)

AztecOO External Convergence/Status Test Set.

Assigns an already defined AztecOO_StatusTest object as the class that
will determine when iterations should stop, either because convergence
was reached or the iteration failed. This method allows a large
variety of convergence tests to be used with AztecOO. The
AztecOO_StatusTest class is a pure virtual class, so any class that
implements its interface can be passed in to this set method. A number
of pre-defined AztecOO_StatusTest derived classes are already
available, including AztecOO_StatusTestCombo, a class that allows
logical combinations of other status test objects for sophisticated
convergence testing. ";

%feature("docstring")  AztecOO::SetOutputStream "void
AztecOO::SetOutputStream(std::ostream &ostrm)

Set std::ostream for Aztec's screen output.

This sets the destination for output that Aztec would normally send to
stdout. ";

%feature("docstring")  AztecOO::SetErrorStream "void
AztecOO::SetErrorStream(std::ostream &errstrm)

Set std::ostream for Aztec's error output.

This sets the destination for output that Aztec would normally send to
stderr. ";

/*  Explicit preconditioner construction/assessment/destruction
methods.  */

%feature("docstring")  AztecOO::ConstructPreconditioner "int
AztecOO::ConstructPreconditioner(double &condest)

Forces explicit construction and retention of an AztecOO native
preconditioner.

AztecOO typically constructs the preconditioner on the first call to
the solve function. However, there are situations where we would like
to compute the preconditioner ahead of time. One particular case is
when we want to confirm that the preconditioner well-conditioned. This
method allows us to precompute the preconditioner. It also provides a
estimate of the condition number of the preconditioner. If  condest is
large, e.g., > 1.0e+14, it is likely the preconditioner will fail. In
this case, using threshold values (available in the incomplete
factorizations) can be used to reduce the condition number.

Note: This method does not work for user-defined preconditioners
(defined via calls to SetPrecOperator(). It will return with an error
code of -1 for this case. ";

%feature("docstring")  AztecOO::DestroyPreconditioner "int
AztecOO::DestroyPreconditioner()

Destroys a preconditioner computed using ConstructPreconditioner().

The ConstructPreconditioner() method creates a persistent
preconditioner. In other words the preconditioner will be used by all
calls to the Iterate() method. DestroyPreconditioner() deletes the
current preconditioner and restores AztecOO to a state where the
preconditioner will computed on first use of the preconditioner solve.
";

%feature("docstring")  AztecOO::Condest "double AztecOO::Condest()
const

Returns the condition number estimate for the current preconditioner,
if one exists, returns -1.0 if no estimate. ";

/*  Check/Attribute Access Methods.  */

%feature("docstring")  AztecOO::CheckInput "int AztecOO::CheckInput()
const

Prints a summary of solver parameters, performs simple sanity checks.
";

%feature("docstring")  AztecOO::GetProblem "Epetra_LinearProblem*
AztecOO::GetProblem() const

Get a pointer to the Linear Problem used to construct this solver;
returns zero if not available. ";

%feature("docstring")  AztecOO::GetUserOperator "Epetra_Operator*
AztecOO::GetUserOperator() const

Get a pointer to the user operator A. ";

%feature("docstring")  AztecOO::GetUserMatrix "Epetra_RowMatrix*
AztecOO::GetUserMatrix() const

Get a pointer to the user matrix A. ";

%feature("docstring")  AztecOO::GetPrecOperator "Epetra_Operator*
AztecOO::GetPrecOperator() const

Get a pointer to the preconditioner operator. ";

%feature("docstring")  AztecOO::GetPrecMatrix "Epetra_RowMatrix*
AztecOO::GetPrecMatrix() const

Get a pointer to the matrix used to construct the preconditioner. ";

%feature("docstring")  AztecOO::GetLHS "Epetra_MultiVector*
AztecOO::GetLHS() const

Get a pointer to the left-hand-side X. ";

%feature("docstring")  AztecOO::GetRHS "Epetra_MultiVector*
AztecOO::GetRHS() const

Get a pointer to the right-hand-side B. ";

%feature("docstring")  AztecOO::PrintLinearSystem "void
AztecOO::PrintLinearSystem(const char *name)

Print linear-system to files.

Parameters:
-----------

name:  Print the matrix to the file A_'name', and print the solution
and rhs vectors to files X_'name' and B_'name', respectively. Will
only produce a matrix file if the run-time-type of the matrix is
either Epetra_CrsMatrix or Epetra_VbrMatrix. ";

/*  Standard AztecOO option and parameter setting methods.  */

%feature("docstring")  AztecOO::SetAztecDefaults "int
AztecOO::SetAztecDefaults()

AztecOO function to restore default options/parameter settings.

This function is called automatically within AztecOO's constructor,
but if constructed using a Epetra_LinearProblem object, some options
are reset based on the ProblemDifficultyLevel associated with the
Epetra_LinearProblem.

See the Aztec 2.1 User Guide for a complete list of these options.

WARNING:  In AztecOO, the default value of options[AZ_poly_ord] is set
to 1. This is different than Aztec 2.1, but the preferred value since
Jacobi preconditioning is used much more often than polynomial
preconditioning and one step of Jacobi is far more effective than 3
steps. ";

%feature("docstring")  AztecOO::SetAztecOption "int
AztecOO::SetAztecOption(int option, int value)

AztecOO option setting function.

Set a specific Aztec option value. Example:
problem.SetAztecOption(AZ_precond, AZ_Jacobi)

See the Aztec 2.1 User Guide for a complete list of these options. ";

%feature("docstring")  AztecOO::GetAztecOption "int
AztecOO::GetAztecOption(int option)

AztecOO option getting function.

Get a specific Aztec optioin value. Example:
problem.GetAztecOption(AZ_precond)

See the Aztec 2.1 User Guide for a complete list of these options. ";

%feature("docstring")  AztecOO::SetAztecParam "int
AztecOO::SetAztecParam(int param, double value)

AztecOO param setting function.

Set a specific Aztec parameter value. Example:
problem.SetAztecParam(AZ_drop, 1.0E-6)

See the Aztec 2.1 User Guide for a complete list of these parameters.
";

%feature("docstring")  AztecOO::GetAllAztecOptions "const int*
AztecOO::GetAllAztecOptions() const

AztecOO option setting function.

Return a pointer to an array (size AZ_OPTIONS_SIZE) of all of the
currently set aztec options. ";

%feature("docstring")  AztecOO::GetAllAztecParams "const double*
AztecOO::GetAllAztecParams() const

AztecOO param setting function.

Return a pointer to an array (size AZ_PARAMS_SIZE) of all of the
currently set aztec parameters. ";

%feature("docstring")  AztecOO::SetAllAztecOptions "int
AztecOO::SetAllAztecOptions(const int *options)

AztecOO option setting function.

Set all Aztec option values using an existing Aztec options array. ";

%feature("docstring")  AztecOO::SetAllAztecParams "int
AztecOO::SetAllAztecParams(const double *params)

AztecOO param setting function.

Set all Aztec parameter values using an existing Aztec params array.
";

/*  Standard AztecOO solve methods.  */

%feature("docstring")  AztecOO::Iterate "int AztecOO::Iterate(int
MaxIters, double Tolerance)

AztecOO iteration function.

Iterates on the current problem until MaxIters or Tolerance is
reached. ";

%feature("docstring")  AztecOO::Iterate "int
AztecOO::Iterate(Epetra_RowMatrix *A, Epetra_MultiVector *X,
Epetra_MultiVector *B, int MaxIters, double Tolerance)

AztecOO iteration function.

Iterates on the specified matrix and vectors until MaxIters or
Tolerance is reached.. ";

/*  Specialist AztecOO solve method.  */

%feature("docstring")  AztecOO::recursiveIterate "int
AztecOO::recursiveIterate(int MaxIters, double Tolerance)

AztecOO iteration functions.

Iterates on the current problem until MaxIters or Tolerance is
reached.. This one should be suitable for recursive invocations of
Aztec. ";

%feature("docstring")  AztecOO::GetAztecStatus "const double*
AztecOO::GetAztecStatus() const

Return the Aztec status after iterating.

Returns pointer to the underlying Aztec Status array (of length
AZ_STATUS_SIZE). See the Aztec documenation. ";

/*  Adaptive Solve methods.  */

%feature("docstring")  AztecOO::SetUseAdaptiveDefaultsTrue "int
AztecOO::SetUseAdaptiveDefaultsTrue()

Force the AdaptiveIterate() method to use default adaptive strategy.
";

%feature("docstring")  AztecOO::SetAdaptiveParams "int
AztecOO::SetAdaptiveParams(int NumTrials, double *athresholds, double
*rthresholds, double condestThreshold, double maxFill, int maxKspace)

Set the parameter that control the AdaptiveIterate() method.

The AdaptiveIterate() method attempts to solve a given problem using
multiple preconditioner and iterative method tuning parameters. There
are defaults that are coded into AdaptiveIterate() method, but the
defaults can be over-ridden by the use of the SetAdaptiveParams()
method. Details of condition number management follow:

Parameters:
-----------

NumTrials:  In The number of Athresh and Rthresh pairs that should be
tried when attempting to stabilize the preconditioner.

athresholds:  In The list of absolute threshold values that should be
tried when attempting to stabilize the preconditioner.

rthresholds:  In The list of relative threshold values that should be
tried when attempting to stabilize the preconditioner.

condestThreshold:  In If the condition number estimate of the
preconditioner is above this number, no attempt will be made to try
iterations. Instead a new preconditioner will be computed using the
next threshold pair.

maxFill:  In In addition to managing the condest, the
AdaptiveIterate() method will also try to increase the preconditioner
fill if it is determined that this might help. maxFill specifies the
maximum fill allowed.

maxKspace:  In In addition to managing the condest, the
AdaptiveIterate() method will also try to increase the Krylov subspace
size if GMRES is being used and it is determined that this might help.
maxKspace specifies the maximum Krylov subspace allowed. ";

%feature("docstring")  AztecOO::AdaptiveIterate "int
AztecOO::AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double
Tolerance)

Attempts to solve the given linear problem using an adaptive strategy.
";

/*  Post-solve access functions  */

%feature("docstring")  AztecOO::NumIters "int AztecOO::NumIters()
const

Returns the total number of iterations performed on this problem. ";

%feature("docstring")  AztecOO::TrueResidual "double
AztecOO::TrueResidual() const

Returns the true unscaled residual for this problem. ";

%feature("docstring")  AztecOO::ScaledResidual "double
AztecOO::ScaledResidual() const

Returns the true scaled residual for this problem. ";

%feature("docstring")  AztecOO::RecursiveResidual "double
AztecOO::RecursiveResidual() const

Returns the recursive residual for this problem. ";

%feature("docstring")  AztecOO::SolveTime "double
AztecOO::SolveTime() const

Returns the solve time. ";

%feature("docstring")  AztecOO::GetAllAztecStatus "int
AztecOO::GetAllAztecStatus(double *status)

AztecOO status extraction function.

Extract Aztec status array into user-provided array. The array must be
of length AZ_STATUS_SIZE as defined in the az_aztec.h header file. ";

%feature("docstring")  AztecOO::SetPreconditioner "int
AztecOO::SetPreconditioner(AZ_PRECOND *Prec)

AztecOO External Preconditioner Set (object)

Associates an already defined Aztec preconditioner with this solve. ";

%feature("docstring")  AztecOO::SetPreconditioner "int
AztecOO::SetPreconditioner(AZ_PREC_FUN prec_function, void *prec_data)

AztecOO External Preconditioner Set (function and data)

Associates an external function and data pointer with preconditioner
";

%feature("docstring")  AztecOO::SetScaling "int
AztecOO::SetScaling(struct AZ_SCALING *Scaling)

AztecOO External Scaling Set.

Associates an already defined Aztec scaling object with this solve. ";

%feature("docstring")  AztecOO::SetMatrixName "int
AztecOO::SetMatrixName(int label)

AztecOO Label Matrix for Aztec.

This is used to label individual matrices within Aztec. This might be
useful if several Aztec invocations are involved corresponding to
different matrices. ";

%feature("docstring")  AztecOO::SetLabel "void
AztecOO::SetLabel(const char *const Label)

Set Label this AztecOO object.

Defines the label used to describe the this object. ";

%feature("docstring")  AztecOO::GetLabel "const char *
AztecOO::GetLabel() const

Get the label describing this AztecOO object.

Returns the string used to define this object. ";


// File: classAztecOO__Operator.xml
%feature("docstring") AztecOO_Operator "

AztecOO_Operator: An implementation of the Epetra_Operator class.

The AztecOO_Operator class implements Epetra_Operator using a pre-
constructed AztecOO solver object. Once constructed, an
AztecOO_Operator can be used as a preconditioner within another
AztecOO solver object.

C++ includes: AztecOO_Operator.h ";

%feature("docstring")  AztecOO_Operator::AztecOO_Operator "AztecOO_Operator::AztecOO_Operator(AztecOO *solver, int NumIters,
double Tol=0.0)

Uses an AztecOO instance to implement the Epetra_Operator interface.

Facilitates the use of an AztecOO solver instance as an operator. This
is particularly designed for using AztecOO as a preconditioner within
another AztecOO instance.

Parameters:
-----------

In:  - A fully-constructed AztecOO object.

In:  - Number of iterations that should be performed. Exactly this
many iterations will be done if Tol = 0.0.

In:  - Tolerance used for each application of AztecOO solver. ";

%feature("docstring")  AztecOO_Operator::~AztecOO_Operator "AztecOO_Operator::~AztecOO_Operator()

Destructor. ";

%feature("docstring")  AztecOO_Operator::SetUseTranspose "int
AztecOO_Operator::SetUseTranspose(bool use_transpose)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  use_transpose - If true, multiply by the transpose of operator,
otherwise just use operator.

WARNING:  - This returns -1 if use_transpose is true, because tranpse
is not supported. ";

%feature("docstring")  AztecOO_Operator::Apply "int
AztecOO_Operator::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a AztecOO_Operator applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

WARNING:  - This method has no effect and returns -1 as error code. ";

%feature("docstring")  AztecOO_Operator::ApplyInverse "int
AztecOO_Operator::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a AztecOO_Operator inverse applied to an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  AztecOO_Operator::NormInf "double
AztecOO_Operator::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  AztecOO_Operator::Label "const char*
AztecOO_Operator::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  AztecOO_Operator::Solver "AztecOO*
AztecOO_Operator::Solver() const

Returns a pointer to the AztecOO solver object that was used to create
this AztecOO_Operator object. ";

%feature("docstring")  AztecOO_Operator::NumIters "int
AztecOO_Operator::NumIters() const

Returns the number of iterations that will be performed with the
AztecOO solver. ";

%feature("docstring")  AztecOO_Operator::Tol "double
AztecOO_Operator::Tol() const

Returns the tolerance this will be used by the AztecOO solver. ";

%feature("docstring")  AztecOO_Operator::UseTranspose "bool
AztecOO_Operator::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  AztecOO_Operator::HasNormInf "bool
AztecOO_Operator::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  AztecOO_Operator::Comm "const Epetra_Comm&
AztecOO_Operator::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  AztecOO_Operator::OperatorDomainMap "const
Epetra_Map& AztecOO_Operator::OperatorDomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  AztecOO_Operator::OperatorRangeMap "const
Epetra_Map& AztecOO_Operator::OperatorRangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";


// File: classAztecOO__StatusTest.xml
%feature("docstring") AztecOO_StatusTest "

AztecOO_StatusTest: A pure virtual class for extending the status
testing capabilities of AztecOO.

C++ includes: AztecOO_StatusTest.h ";

%feature("docstring")  AztecOO_StatusTest::AztecOO_StatusTest "AztecOO_StatusTest::AztecOO_StatusTest()

Constructor. ";

%feature("docstring")  AztecOO_StatusTest::~AztecOO_StatusTest "virtual AztecOO_StatusTest::~AztecOO_StatusTest()

Destructor. ";

%feature("docstring")  AztecOO_StatusTest::ResidualVectorRequired "virtual bool AztecOO_StatusTest::ResidualVectorRequired() const =0

Indicates if residual vector is required by this convergence test.

If this method returns true, then the ResidualVector argument to the
Converged() method will defined. If this method returns false, then
the ResidualVector may not be defined when Converged() is called. Some
iterative methods do not explicitly construct the residual vector at
each iteration. Thus, if this vector is not required, this vector will
not need to be constructed if this method returns false. ";

%feature("docstring")  AztecOO_StatusTest::CheckStatus "virtual
AztecOO_StatusType AztecOO_StatusTest::CheckStatus(int CurrentIter,
Epetra_MultiVector *CurrentResVector, double CurrentResNormEst, bool
SolutionUpdated)=0

Check convergence status: Unconverged, Converged, Failed.

This method checks to see if the convergence criteria are met. The
calling routine may pass in the current native residual vector (the
one naturally produced as part of the iterative method) or a pre-
computed estimate of the two-norm of the current residual, or both or
neither. The calling routine should also indicate if the solution of
the linear problem has been updated to be compatible with the
residual. Some methods, such as GMRES do update the solution at each
iteration.

Parameters:
-----------

CurrentIter:  (In) Current iteration of iterative method.

CurrentResVector:  (In) The current residuals of the iterative
process. These values are assumed to be the residuals that are a
\"natural\" by-product of the iterative process. Typically this means
they are not explicitly generated and are therefore subject to round-
off error. These values will also reflect the influence of any Any
rigorous use of stopping criteria should not rely solely on results
using this vector. Instead, tests can be performed using this vector
and, after convergence is reached with this vector.

CurrentResNormEst:  (In) If the iterative method can cheaply provide
this estimate, as an alternative or in addition to providing the
CurrentResVector, this value will contain that estimate. The value
will be set to -1.0 if no estimate is available.

SolutionUpdated:  (In) If this argument is true, then the solution
vector that is part of the Epetra_LinearProblem object being solved is
consistent with the residual. Some iterative methods do not keep the
solution vector updated with the residual at each iteration. For
example GMRES does not generate the solution until the least- square
problem is solved at the end of the Arnoldi process.

AztecOO_StatusType: Unconverged, Converged or Failed. ";

%feature("docstring")  AztecOO_StatusTest::GetStatus "virtual
AztecOO_StatusType AztecOO_StatusTest::GetStatus() const =0

Return the result of the most recent checkStatus call. ";

%feature("docstring")  AztecOO_StatusTest::Print "virtual ostream&
AztecOO_StatusTest::Print(ostream &stream, int indent=0) const =0

Output formatted description of stopping test to output stream. ";

%feature("docstring")  AztecOO_StatusTest::PrintStatus "virtual void
AztecOO_StatusTest::PrintStatus(ostream &os, AztecOO_StatusType type)
const ";


// File: classAztecOO__StatusTestCombo.xml
%feature("docstring") AztecOO_StatusTestCombo "

AztecOO_StatusTestCombo: A class for extending the status testing
capabilities of AztecOO via logical combinations.

AztecOO_StatusTestCombo is an interface that can be implemented to
extend the convergence testing capabilities of AztecOO. This class
supports composite tests. In this situation, two or more existing
AztecOO_StatusTestCombo objects test1 and test2 can be used to create
a new test. For all combinations, if any tests returns Failed or
returns not-a-number (NaN) status, then the combination test returns
Failed. There are three possible combinations: OR combination: If an
OR combination is selected, the status returns Converged if any one of
the subtest returns as Converged.

AND combination: If an AND combination is selected, the status returns
Converged only when all subtests return as Converged.

SEQ combination: SEQ is a form of AND that will perform subtests in
sequence. If the first test returns Unconverged, Failed or NaN, no
other subtests are done, and the status is returned as Unconverged if
the first test was Unconverged, or as Failed if the first test was
Failed or NaN. If the first test returns Converged, the second test is
checked in the same fashion as the first. If the second test is
Converged, the third one is tested, and so on.

The purpose of the SEQ combination is to allow the addition of
expensive but more rigorous convergence tests. For example, we could
define a test that used the implicit residual vector (the one produced
by the iterative method) as the first subtest and define a second test
using the explicitly computed residual vector. Explicitly computing
the residual requires a matrix multiplication with the original matrix
operator, an expensive operation. By using the SEQ combination, we can
avoid the matrix multiplication associated with the explicit residual
calculation until the implicit residual is small.

WARNING:  Presently it is not valid to associate one status test
instance with two different AztecOO objects.

C++ includes: AztecOO_StatusTestCombo.h ";

%feature("docstring")
AztecOO_StatusTestCombo::AztecOO_StatusTestCombo "AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t)

Constructor. ";

%feature("docstring")
AztecOO_StatusTestCombo::AztecOO_StatusTestCombo "AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t,
AztecOO_StatusTest &a)

Constructor with a single test. ";

%feature("docstring")
AztecOO_StatusTestCombo::AztecOO_StatusTestCombo "AztecOO_StatusTestCombo::AztecOO_StatusTestCombo(ComboType t,
AztecOO_StatusTest &a, AztecOO_StatusTest &b)

Constructor with two tests. ";

%feature("docstring")  AztecOO_StatusTestCombo::AddStatusTest "AztecOO_StatusTestCombo &
AztecOO_StatusTestCombo::AddStatusTest(AztecOO_StatusTest &a)

Add another test to this combination. ";

%feature("docstring")
AztecOO_StatusTestCombo::~AztecOO_StatusTestCombo "virtual
AztecOO_StatusTestCombo::~AztecOO_StatusTestCombo()

Destructor. ";

%feature("docstring")  AztecOO_StatusTestCombo::ResidualVectorRequired
"bool AztecOO_StatusTestCombo::ResidualVectorRequired() const

Indicates if residual vector is required by this convergence test.

If this method returns true, then one or more of the
AztecOO_StatusTest objects that make up this combined test requires
the Residual Vector to perform its test. ";

%feature("docstring")  AztecOO_StatusTestCombo::CheckStatus "AztecOO_StatusType AztecOO_StatusTestCombo::CheckStatus(int
CurrentIter, Epetra_MultiVector *CurrentResVector, double
CurrentResNormEst, bool SolutionUpdated)

Check convergence status: Unconverged, Converged, Failed.

This method checks to see if the convergence criteria are met.
Depending on how the combined test is constructed this method will
return the appropriate status type using common logic principals.
However, if any subtest returns with a Failed status type, the
combined test will return a status type of Failed.

Parameters:
-----------

CurrentIter:  (In) Current iteration of iterative method.

CurrentResVector:  (In) The current residuals of the iterative
process.

CurrentResNormEst:  (In) Estimate of the two-norm of the residual. The
value will be set to -1.0 if no estimate is available.

SolutionUpdated:  (In) If this argument is true, then the solution
vector that is part of the Epetra_LinearProblem object being solved is
consistent with the residual.

AztecOO_StatusType: Unconverged, Converged or Failed. ";

%feature("docstring")  AztecOO_StatusTestCombo::GetStatus "AztecOO_StatusType AztecOO_StatusTestCombo::GetStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  AztecOO_StatusTestCombo::Print "ostream &
AztecOO_StatusTestCombo::Print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  AztecOO_StatusTestCombo::GetComboType "ComboType AztecOO_StatusTestCombo::GetComboType() const

Returns the maximum number of iterations set in the constructor. ";


// File: classAztecOO__StatusTestMaxIters.xml
%feature("docstring") AztecOO_StatusTestMaxIters "

AztecOO_StatusTestMaxIters: An AztecOO_StatusTest class specifying a
maximum number of iterations.

C++ includes: AztecOO_StatusTestMaxIters.h ";

%feature("docstring")
AztecOO_StatusTestMaxIters::AztecOO_StatusTestMaxIters "AztecOO_StatusTestMaxIters::AztecOO_StatusTestMaxIters(int MaxIters)

Constructor. ";

%feature("docstring")
AztecOO_StatusTestMaxIters::~AztecOO_StatusTestMaxIters "virtual
AztecOO_StatusTestMaxIters::~AztecOO_StatusTestMaxIters()

Destructor. ";

%feature("docstring")
AztecOO_StatusTestMaxIters::ResidualVectorRequired "bool
AztecOO_StatusTestMaxIters::ResidualVectorRequired() const

Indicates if residual vector is required by this convergence test:
returns false for this class. ";

%feature("docstring")  AztecOO_StatusTestMaxIters::CheckStatus "AztecOO_StatusType AztecOO_StatusTestMaxIters::CheckStatus(int
CurrentIter, Epetra_MultiVector *CurrentResVector, double
CurrentResNormEst, bool SolutionUpdated)

Check convergence status: Unconverged, Converged, Failed.

This method checks to see if the convergence criteria are met..

Parameters:
-----------

CurrentIter:  (In) Current iteration of iterative method. Compared
against MaxIters value passed in at construction. If CurrentIter <
MaxIters, we return with StatusType = Unconverged. Otherwise,
StatusType will be set to Failed.

CurrentResVector:  (In) Ignored by this class.

CurrentResNormEst:  (In) Ignored by this class.

SolutionUpdated:  (In) Ignored by this class.

StatusType Unconverged if CurrentIter<MaxIters, Failed if
CurrentIters>=MaxIters. ";

%feature("docstring")  AztecOO_StatusTestMaxIters::GetStatus "AztecOO_StatusType AztecOO_StatusTestMaxIters::GetStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  AztecOO_StatusTestMaxIters::Print "ostream &
AztecOO_StatusTestMaxIters::Print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  AztecOO_StatusTestMaxIters::GetMaxIters "int
AztecOO_StatusTestMaxIters::GetMaxIters() const

Returns the maximum number of iterations set in the constructor. ";

%feature("docstring")  AztecOO_StatusTestMaxIters::GetNumIters "int
AztecOO_StatusTestMaxIters::GetNumIters() const

Returns the current number of iterations from the most recent
StatusTest call. ";


// File: classAztecOO__StatusTestResNorm.xml
%feature("docstring") AztecOO_StatusTestResNorm "

AztecOO_StatusTestResNorm: An implementation of AztecOO_StatusTest
using a family of residual norms.

AztecOO_StatusTestResNorm is an implementation of AztecOO_StatusTest
that allows a user to construct one of a family of residual tests for
use as a status/convergence test for AztecOO. The form of the test is
\\\\[ \\\\frac{\\\\|r\\\\|}{\\\\sigma} \\\\le \\\\tau \\\\] where  $r$
is the residual vector, implicitly or explicitly computed (determined
by enum ResType),

$\\\\|r\\\\|$ is the residual norm determined by the enum NormType
(1-norm, 2-norm or inf-norm),

$\\\\sigma$ is the scale factor that can be passed in as a precomputed
double precision number, or can be selected from by the enum ScaleType
(norm of RHS, norm of initial residual).

$\\\\tau$ is the tolerance that is passed in as a double precision
number to the constructor. The value of $\\\\tau$ can be reset using
the ResetTolerance() method.

WARNING:  Presently it is not valid to associate one status test
instance with two different AztecOO objects.

C++ includes: AztecOO_StatusTestResNorm.h ";

%feature("docstring")
AztecOO_StatusTestResNorm::AztecOO_StatusTestResNorm "AztecOO_StatusTestResNorm::AztecOO_StatusTestResNorm(const
Epetra_Operator &Operator, const Epetra_Vector &LHS, const
Epetra_Vector &RHS, double Tolerance)

Constructor.

The constructor takes a single argument specifying the tolerance (
$\\\\tau$). If none of the form definition methods are called, we use
$\\\\|r\\\\|_2/\\\\|r^{(0)}\\\\|_2 \\\\le \\\\tau$ as the stopping
criterion, where $\\\\|r\\\\|_2$ uses the least costly form of the
2-norm of residual available from the iterative method and
$\\\\|r^{(0)}\\\\|_2$ is the corresponding norm of the initial
residual. The least costly form of the 2-norm depends on the chosen
iterative method. Most Krylov methods produce the preconditioned
residual vector in a form that would be exact in infinite precision
arithmetic. This vector may be different from the true residual either
because left scaling or preconditioning was used, or because round-off
error has introduced significant error, or both.

Parameters:
-----------

Operator:  (In) The original linear operator that was passed in to the
AztecOO solver object.

LHS:  (In) The original left hand side vector that was passed in to
the AztecOO solver object. NOTE: AztecOO accepts multivector objects,
but AztecOO_StatusTestResNorm does not handle residual tests for
multiple vectors. Most AztecOO users tend to have Epetra_Vectors, even
though they are passing them in as Epetra_MultiVectors, so this should
not be an issue. If you are truly using Epetra_MultiVector objects,
remember that for a multivector object mv, the ith vector can be
extracted as mv(i).

RHS:  (In) The original right hand side vector that was passed in to
the AztecOO solver object. See note for LHS.

Tolerance:  (In) A double value that is used to test the convergence
criterion. Can be reset using ResetTolerance(). ";

%feature("docstring")
AztecOO_StatusTestResNorm::~AztecOO_StatusTestResNorm "AztecOO_StatusTestResNorm::~AztecOO_StatusTestResNorm()

Destructor. ";

%feature("docstring")  AztecOO_StatusTestResNorm::DefineResForm "int
AztecOO_StatusTestResNorm::DefineResForm(ResType TypeOfResidual,
NormType TypeOfNorm, Epetra_Vector *Weights=0)

Define form of the residual, its norm and optional weighting vector.

This method defines the form of $\\\\|r\\\\|$. We specify: Whether the
residual vector should be explicitly computed, or taken from the
iterative method.

The norm to be used on the residual (this may be different than the
norm used in DefineScaleForm()).

A weighting vector that will be multiplied, element by element, with
the residual vector prior to computing the norm. This argument
defaults to a zero vector and no weighting will be done. ";

%feature("docstring")  AztecOO_StatusTestResNorm::DefineScaleForm "int AztecOO_StatusTestResNorm::DefineScaleForm(ScaleType
TypeOfScaling, NormType TypeOfNorm, Epetra_Vector *Weights=0, double
ScaleValue=1.0)

Define form of the scaling, its norm, its optional weighting vector,
or, alternatively, define an explicit value.

This method defines the form of how the residual is scaled (if at
all). It operates in two modes: User-provided scaling value: Set
argument TypeOfScaling to UserProvided.

Set ScaleValue to a non-zero value that the residual norm will be
divided by.

TypeOfNorm and Weights arguments will be ignored.

Sample use: Define ScaleValue = $\\\\|A\\\\|_{\\\\infty}$ where $ A $
is the matrix of the linear problem.

Use a supported Scaling Form: Define TypeOfScaling to be the norm of
the right hand side, the initial residual vector, or to none.

Define norm to be used on the scaling vector (this may be different
than the norm used in DefineResForm()).

Define a weighting vector that will be multiplied, element by element,
with the residual vector prior to computing the norm. This argument
defaults to a zero vector and no weighting will be done. ";

%feature("docstring")  AztecOO_StatusTestResNorm::ResetTolerance "int
AztecOO_StatusTestResNorm::ResetTolerance(double Tolerance)

Reset the value of the tolerance.

We allow the tolerance to be reset for cases where, in the process of
testing the residual, we find that the initial tolerance was too tight
or too lax. ";

%feature("docstring")
AztecOO_StatusTestResNorm::SetMaxNumExtraIterations "int
AztecOO_StatusTestResNorm::SetMaxNumExtraIterations(int
maxNumExtraIterations)

Set the maximum number of extra iterations that will be performed in
case the implicit residual succeeds that the true residual fails
(default is 0).

In some instance,especially with GMRES, the implicitly computed
residual is an optimistic estimate of the true residual. In these
cases, especially when the tolerance is set very small, the iterative
solver can never satisfy the tolerance with the explicit residual, so
we allow the user to limit the number of extra iterations that will be
performed. If the implicit residual is satisfied, then the value set
here will determine exactly how many extra iterations will be allowed
to try and converge the explicit residual. This value has no impact
status tests that are based on explicit residuals.

Parameters:
-----------

maxNumExtraIterations:  (In) Maximum number of extra iterations that
will be performed in case the implicit residual succeeds that the true
residual fails. ";

%feature("docstring")
AztecOO_StatusTestResNorm::GetMaxNumExtraIterations "int
AztecOO_StatusTestResNorm::GetMaxNumExtraIterations()

Return the maximum number of extra iterations that are performed if
the implicit residual succeeds while the true residual fails. ";

%feature("docstring")
AztecOO_StatusTestResNorm::ResidualVectorRequired "bool
AztecOO_StatusTestResNorm::ResidualVectorRequired() const

Indicates if residual vector is required by this convergence test.

The value returned by this method will depend on several factors. Once
an AztecOO_StatusTestResNorm object is constructed and the
DefineResForm and DefineScaleForm methods are optionally called, this
method can tested. For most Krylov solvers, there is no extra cost to
providing the residual vector. However, GMRES and Transpose-free QMR
will need to explicitly compute this vector if
ResidualVectorRequired() returns true, so this is an extra cost for
these two iterative methods. ";

%feature("docstring")  AztecOO_StatusTestResNorm::CheckStatus "AztecOO_StatusType AztecOO_StatusTestResNorm::CheckStatus(int
CurrentIter, Epetra_MultiVector *CurrentResVector, double
CurrentResNormEst, bool SolutionUpdated)

Check convergence status: Unconverged, Converged, Failed.

This method checks to see if the convergence criteria are met.
Depending on how the residual test is constructed this method will
return the appropriate status type.

Parameters:
-----------

CurrentIter:  (In) Current iteration of iterative method.

CurrentResVector:  (In) The current residuals of the iterative
process.

CurrentResNormEst:  (In) Estimate of the two-norm of the residual. The
value will be set to -1.0 if no estimate is available.

SolutionUpdated:  (In) If this argument is true, then the solution
vector that is part of the Epetra_LinearProblem object being solved is
consistent with the residual.

AztecOO_StatusType: Unconverged, Converged or Failed. ";

%feature("docstring")  AztecOO_StatusTestResNorm::GetStatus "AztecOO_StatusType AztecOO_StatusTestResNorm::GetStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  AztecOO_StatusTestResNorm::Print "ostream &
AztecOO_StatusTestResNorm::Print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  AztecOO_StatusTestResNorm::ResetStatus "void
AztecOO_StatusTestResNorm::ResetStatus()

Reset state of status test object. ";

%feature("docstring")  AztecOO_StatusTestResNorm::GetTolerance "double AztecOO_StatusTestResNorm::GetTolerance() const

Returns the value of the tolerance, $ \\\\tau $, set in the
constructor. ";

%feature("docstring")  AztecOO_StatusTestResNorm::GetTestValue "double AztecOO_StatusTestResNorm::GetTestValue() const

Returns the test value, $ \\\\frac{\\\\|r\\\\|}{\\\\sigma} $, computed
in most recent call to CheckStatus. ";

%feature("docstring")  AztecOO_StatusTestResNorm::GetResNormValue "double AztecOO_StatusTestResNorm::GetResNormValue() const

Returns the residual norm value, $ \\\\|r\\\\| $, computed in most
recent call to CheckStatus. ";

%feature("docstring")  AztecOO_StatusTestResNorm::GetScaledNormValue "double AztecOO_StatusTestResNorm::GetScaledNormValue() const

Returns the scaled norm value, $ \\\\sigma $. ";


// File: classAztecOOConditionNumber.xml
%feature("docstring") AztecOOConditionNumber "

Condition number estimator using AztecOO.

This object will estimate the condition number of an Epetra_Operator.

C++ includes: AztecOO_ConditionNumber.h ";

%feature("docstring")  AztecOOConditionNumber::AztecOOConditionNumber
"AztecOOConditionNumber::AztecOOConditionNumber()

Constructor. ";

%feature("docstring")  AztecOOConditionNumber::~AztecOOConditionNumber
"AztecOOConditionNumber::~AztecOOConditionNumber()

Destructor. ";

%feature("docstring")  AztecOOConditionNumber::initialize "void
AztecOOConditionNumber::initialize(const Epetra_Operator &op,
SolverType solverType=GMRES_, int krylovSubspaceSize=100, bool
printSolve=false)

Initialization. ";

%feature("docstring")  AztecOOConditionNumber::computeConditionNumber
"int AztecOOConditionNumber::computeConditionNumber(int maxIters,
double tol)

Estimates the condition number. ";

%feature("docstring")  AztecOOConditionNumber::getConditionNumber "double AztecOOConditionNumber::getConditionNumber()

Return condition number computed by last call to
computeConditionNumber. ";


// File: structAztecOO_1_1MatrixData.xml
%feature("docstring") AztecOO::MatrixData "";

%feature("docstring")  AztecOO::MatrixData::MatrixData "AztecOO::MatrixData::MatrixData(Epetra_RowMatrix *inA=0, Epetra_Vector
*inX=0, Epetra_Vector *inY=0, Epetra_Vector *inSourceVec=0,
Epetra_Vector *inTargetVec=0) ";

%feature("docstring")  AztecOO::MatrixData::~MatrixData "AztecOO::MatrixData::~MatrixData() ";


// File: structAztecOO_1_1OperatorData.xml
%feature("docstring") AztecOO::OperatorData "";

%feature("docstring")  AztecOO::OperatorData::OperatorData "AztecOO::OperatorData::OperatorData(Epetra_Operator *inA=0,
Epetra_Vector *inX=0, Epetra_Vector *inY=0) ";

%feature("docstring")  AztecOO::OperatorData::~OperatorData "AztecOO::OperatorData::~OperatorData() ";


// File: AztecOO_8cpp.xml
%feature("docstring")  AztecOO_uppercase "string
AztecOO_uppercase(const string &s) ";

%feature("docstring")  Epetra_Aztec_matnorminf "double
Epetra_Aztec_matnorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_matvec "void
Epetra_Aztec_matvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_operatornorminf "double
Epetra_Aztec_operatornorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_operatorvec "void
Epetra_Aztec_operatorvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_precond "void
Epetra_Aztec_precond(double x[], int input_options[], int
proc_config[], double input_params[], AZ_MATRIX *Amat, AZ_PRECOND
*prec) ";

%feature("docstring")  Epetra_Aztec_getrow "int
Epetra_Aztec_getrow(int columns[], double values[], int row_lengths[],
AZ_MATRIX *Amat, int N_requested_rows, int requested_rows[], int
allocated_space) ";

%feature("docstring")  Epetra_Aztec_comm_wrapper "int
Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) ";

%feature("docstring")  AztecOO_StatusTest_wrapper "void
AztecOO_StatusTest_wrapper(void *conv_test_obj, void *res_vector_obj,
int iteration, double *res_vector, int print_info, int sol_updated,
int *converged, int *isnan, double *rnorm, int *r_avail) ";


// File: AztecOO_8h.xml
%feature("docstring")  Epetra_Aztec_matvec "void
Epetra_Aztec_matvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_matnorminf "double
Epetra_Aztec_matnorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_operatorvec "void
Epetra_Aztec_operatorvec(double x[], double y[], AZ_MATRIX *Amat, int
proc_config[]) ";

%feature("docstring")  Epetra_Aztec_operatornorminf "double
Epetra_Aztec_operatornorminf(AZ_MATRIX *Amat) ";

%feature("docstring")  Epetra_Aztec_precond "void
Epetra_Aztec_precond(double x[], int input_options[], int
proc_config[], double input_params[], AZ_MATRIX *Amat, AZ_PRECOND
*prec) ";

%feature("docstring")  Epetra_Aztec_getrow "int
Epetra_Aztec_getrow(int columns[], double values[], int row_lengths[],
AZ_MATRIX *Amat, int N_requested_rows, int requested_rows[], int
allocated_space) ";

%feature("docstring")  Epetra_Aztec_comm_wrapper "int
Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) ";

%feature("docstring")  AztecOO_StatusTest_wrapper "void
AztecOO_StatusTest_wrapper(void *conv_test_obj, void *res_vector_obj,
int iteration, double *res_vector, int print_info, int sol_updated,
int *converged, int *isnan, double *rnorm, int *r_avail) ";


// File: AztecOO__ConditionNumber_8cpp.xml


// File: AztecOO__ConditionNumber_8h.xml


// File: AztecOO__ConfigDefs_8h.xml


// File: AztecOO__Operator_8cpp.xml


// File: AztecOO__Operator_8h.xml


// File: AztecOO__Scaling_8cpp.xml
%feature("docstring")  AZOO_Scale "int AZOO_Scale(int action,
Epetra_RowMatrix *A, double b[], double x[], int options[], AZ_SCALING
*scaling) ";

%feature("docstring")  AZOO_create_scaling_vector "Epetra_Vector *
AZOO_create_scaling_vector(Epetra_RowMatrix *A, int scaling_type) ";

%feature("docstring")  AztecOO_scale_epetra "int
AztecOO_scale_epetra(int action, AZ_MATRIX *Amat, int options[],
double b[], double x[], int proc_config[], AZ_SCALING *scaling)

A scaling function that can be assigned to the AZ_SCALING.scale
function-pointer, to provide a call-back that Aztec can use to scale
Epetra matrices passed in by AztecOO. ";


// File: AztecOO__Scaling_8h.xml
%feature("docstring")  AztecOO_scale_epetra "int
AztecOO_scale_epetra(int action, AZ_MATRIX *Amat, int options[],
double b[], double x[], int proc_config[], AZ_SCALING *scaling)

A scaling function that can be assigned to the AZ_SCALING.scale
function-pointer, to provide a call-back that Aztec can use to scale
Epetra matrices passed in by AztecOO. ";


// File: AztecOO__StatusTest_8h.xml


// File: AztecOO__StatusTestCombo_8cpp.xml


// File: AztecOO__StatusTestCombo_8h.xml


// File: AztecOO__StatusTestMaxIters_8cpp.xml


// File: AztecOO__StatusTestMaxIters_8h.xml


// File: AztecOO__StatusTestResNorm_8cpp.xml


// File: AztecOO__StatusTestResNorm_8h.xml


// File: AztecOO__StatusType_8h.xml


// File: AztecOO__string__maps_8cpp.xml


// File: AztecOO__string__maps_8h.xml


// File: AztecOO__Version_8h.xml
%feature("docstring")  AztecOO_Version "string AztecOO_Version() ";


// File: dir_7894d9566ccc3a34f56a76d817e62719.xml


// File: dir_d2a30b9de7594b3f8f120c191973da80.xml


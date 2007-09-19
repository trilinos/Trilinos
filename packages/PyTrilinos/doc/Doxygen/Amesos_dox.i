
// File: index.xml

// File: classAmesos__BaseSolver.xml
%feature("docstring") Amesos_BaseSolver "

Amesos_BaseSolver: A pure virtual class for direct solution of real-
valued double- precision operators.

The Amesos_BaseSolver class is a pure virtual class (that is, it
specifies interface only) that enables the use of real-valued double-
precision direct sparse solvers. Every Amesos class named Amesos_
SolverName derives from Amesos_BaseSolver.

Usage Examples

Basic calling sequence

The basic calling sequence solves A x = b or AT x = b without
specifying how A has changed between each call to Solve().

Re-using the symbolic factorization

The following calling sequence performs multiple solves of A x = b or
AT x = b in cases where the non-zero structure of A remains unchanged
between each call to Solve().

Re-using the numeric factorization

The following calling sequence performs multiple solves of A x = b or
AT x = b provided that A remains unchanged between each call to
Solve().

Constructor requirements

Every Amesos_SolverName class should accept an Epetra_LinearProblem

Mathematical methods

Four mathematical methods are defined in the base class
Amesos_BaseSolver: SymbolicFactorization(), NumericFactorization(),
and Solve().

Switching concrete classes

Different concrete classes, each based on a different third party
solver, will have different performance characteristics and will
accept different parameters.

Changing the values of the underlying matrix operator.

Any changes to the values of a matrix must be accompanied by a call to
NumericFactorization() before the next call to Solve() or the behavior
of Solve() is undefined. Any changes to the numerical structure of the
matrix must be followed by a call to SymbolicFactorization() and
NumericalFactorization() before the next call to Solve().

Once SymbolicFactorization() has been called, classes implementing
this interface may assume that any change made to the non-zero
structure of the underlying matrix will be accompanied by a call to
SymbolicFactorization() prior to a subsequent call to
NumericFactorization or Solve().

Named Parameters

Parameters can be changed or added at any time by calling
SetParameters(ParamList) with the new parameters specified in
ParamList.

It is left to the user to be sure that changes made to the parameters
are appropriate for the concrete class that they are using.

Examples of appropriate changes in parameters include:  Changing
iterative refinement rules between calls to Solve()

Changing drop tolerance rules between calls to NumericFactorization()

Examples of inappropriate changes in parameters include:  Changing
drop tolerance rules between solve steps.
Solver.NumericFactorization();
Solver.getList()->set(\"DropTolerance\",.001); Solver.Solve();
Results of making inappropriate changes in parameters is unpredictable
and could include an error return, a bogus result or ignoring the
parameter change.

Transpose solve

Any class implementing Amesos_BaseSolver should handle calls to
SetUseTranspose() at any point. However, the result of a call to
SetUseTranspose() which is not followed by a call to
SymbolicFactorization() and NumericFactorization() is implementation
dependent. Some third party libraries are able to solve AT x = b and
Ax = b using the same factorization. Others will require a new
factorization anytime that a call to SetUseTranspose() changes the
intended solve from AT x = b to Ax = b or vice-versa.

Performance expectations

The following is a list of performance guidelines that classes which
implement the Amesos_BaseSolver class are expected to maintain.

Memory usage:

For serial codes, no more than one extra copy of the original matrix
should be required. Except that some codes require matrix transpostion
which requires additional copies of the input matrix.

For distributed memory codes, no serial copies of the original matrix
should be required.

Robustness requirements

Failures should be caught by AMESOS_CHK_ERR(). The following error
codes should be used: 1: Singular matrix

2: Non-symmetric matrix

3: Matrix is not positive definite

4: Insufficient memory

Because we do not check to see if a matrix has changed between the
call to SymbolicFactorization() and the call to
NumericFactorization(), it is possible that a change to the matrix
will cause a potentially catastrophic error.

C++ includes: Amesos_BaseSolver.h ";


// File: classAmesos__Btf.xml
%feature("docstring") Amesos_Btf "

Amesos_Btf: Factors and solves a matrix after converting it to block
triangular form.

Amesos_Btf: Compute an ordering which reduces the matrix to block
upper triangular form.

Determine a task partitioning (static load balancing)

Redistribute the data based on the owner computes rule The process(es)
to which a given diagonal block is assigned owns all rows in that
diagonal block.

Uses an Amesos solver on each of the diagonal blocks

Uses Trilinos operations to handle the off-diagonal blocks.

C++ includes: Amesos_BTF.h ";


// File: classAmesos__Component.xml
%feature("docstring") Amesos_Component "

Amesos_Component: A pure virtual class for direct solvers to be used
within Amesos_Merikos to form a parallel direct solver.

The Amesos_Component interface specifies what Amesos_Merikos needs.
Any Amesos class that implements Amesos_Component can be used by
Amesos_Merikos to perform partial solves on subblocks of the matrix.

Member functions added by Amesos_Component.

PartialFactorization() PartialFactorization performs factors at most
the first SubMatrixSize_ rows and columns.

PartialFactorization delays the factorization of any columns which
generate unstable (i.e. too small) pivots.

PartialFactorization computes and returns the schur complement.

PartialFactorization does not need a symbolic factorization phase. It
uses the permutation given by SetRowPermutation.

Lsolve performs a raw partial solve, treating the unfactored rows and
columns as the identity without row or column permutation.

Usolve performs a raw partial solve, treating the unfactored rows and
columns as the identity without row or column permutation.

SetRowPermutation - sets the row permutation

GetRowPermutation - gets the row permutation

SetColumnPermutation - sets the column permutation

GetColumnPermutation - gets the column permutation

SetSubMatrixSize - Sets the maximum number of rows (and columns) to
factor.

GetSubMatrixSize - Returns the number of rows (and columns) actually
factored.

SchurComplement - Returns the Schur complement, i.e.
L21(SubMatrixSize+1:MatrixSize,1:SubMatrixSize) *
U12(1:SubMatrixSize,SubMatrixSize+1:MatrixSize)

Usage Examples

Basic calling sequence

Epetra_LinearProblem Problem(A,X,B);     Amesos_SolverName
Solver(Problem);

Solver.PartialFactorization() ;        ... Ancestor factorization
Solver.Lsolve() ;        ... Ancestor solves     Solver.Usolve() ;

Preconditions:  An ordering  Postconditions: Constructor requirements

Every Amesos_SolverName class should accept an Epetra_LinearProblem

C++ includes: Amesos_Component.h ";


// File: classAmesos__ComponentBaseSolver.xml
%feature("docstring") Amesos_ComponentBaseSolver "

Amesos_ComponentBaseSolver: A pure virtual class for direct solvers to
be used within Amesos_Merikos to form a parallel direct solver.

The Amesos_ComponentBaseSolver interface specifies what Amesos_Merikos
needs. Any Amesos class that implements Amesos_ComponentBaseSolver can
be used by Amesos_Merikos to perform partial solves on subblocks of
the matrix.

Member functions added by Amesos_ComponentBaseSolver.

PartialFactorization() PartialFactorization performs factors at most
the first SubMatrixSize_ rows and columns.

PartialFactorization delays the factorization of any columns which
generate unstable (i.e. too small) pivots.

PartialFactorization computes and returns the schur complement.

PartialFactorization does not need a symbolic factorization phase. It
uses the permutation given by SetRowPermutation.

Lsolve performs a raw partial solve, treating the unfactored rows and
columns as the identity without row or column permutation.

Usolve performs a raw partial solve, treating the unfactored rows and
columns as the identity without row or column permutation.

SetRowPermutation - sets the row permutation

GetRowPermutation - gets the row permutation

SetColumnPermutation - sets the column permutation

GetColumnPermutation - gets the column permutation

SetSubMatrixSize - Sets the maximum number of rows (and columns) to
factor.

GetSubMatrixSize - Returns the number of rows (and columns) actually
factored.

SchurComplement - Returns the Schur complement, i.e.
L21(SubMatrixSize+1:MatrixSize,1:SubMatrixSize) *
U12(1:SubMatrixSize,SubMatrixSize+1:MatrixSize)

Usage Examples

Basic calling sequence

Epetra_LinearProblem Problem(A,X,B);     Amesos_SolverName
Solver(Problem);

Solver.PartialFactorization() ;        ... Ancestor factorization
Solver.Lsolve() ;        ... Ancestor solves     Solver.Usolve() ;

Preconditions:  An ordering  Postconditions: Constructor requirements

Every Amesos_SolverName class should accept an Epetra_LinearProblem

C++ includes: Amesos_ComponentBaseSolver.h ";


// File: classAmesos__Control.xml
%feature("docstring") Amesos_Control "

Amesos_Control: Container for some control variables.

Marzio Sala, SNL 9214

C++ includes: Amesos_Control.h ";

%feature("docstring")  Amesos_Control::Amesos_Control "Amesos_Control::Amesos_Control()

Default constructor. ";

%feature("docstring")  Amesos_Control::~Amesos_Control "Amesos_Control::~Amesos_Control()

Default destructor. ";

%feature("docstring")  Amesos_Control::SetControlParameters "void
Amesos_Control::SetControlParameters(const Teuchos::ParameterList
&ParameterList) ";


// File: classAmesos__Dscpack.xml
%feature("docstring") Amesos_Dscpack "

Amesos_Dscpack: An object-oriented wrapper for Dscpack.

Amesos_Dscpack will solve a linear systems of equations: A X = B using
Epetra objects and the Dscpack solver library, where A is an
Epetra_RowMatrix and X and B are Epetra_MultiVector objects.

C++ includes: Amesos_Dscpack.h ";


// File: classAmesos__Dscpack__Pimpl.xml
%feature("docstring") Amesos_Dscpack_Pimpl "";


// File: classAmesos__Klu.xml
%feature("docstring") Amesos_Klu "

Interface to KLU internal solver.

C++ includes: Amesos_Umfpack.h ";


// File: classAmesos__Klu__Pimpl.xml
%feature("docstring") Amesos_Klu_Pimpl "";


// File: classAmesos__Lapack.xml
%feature("docstring") Amesos_Lapack "

Amesos_Lapack: an interface to LAPACK.

Class Amesos_Lapack enables the solution of the distributed linear
system, defined by an Epetra_LinearProblem, using LAPACK.

Amesos_Lapack stores the lineaar system matrix as an
Epetra_SerialDensMatrix. The linear problem is an
Epetra_SerialDenseProblem. Amesos_Lapack factorizes the matrix using
DGETRF().

Marzio Sala, 9214.

C++ includes: Amesos_Lapack.h ";


// File: classAmesos__MC64.xml
%feature("docstring") Amesos_MC64 "

Interface to MC64, reordering and scaling algorithm.

Marzio Sala, ETHZ.

C++ includes: Amesos_MC64.h ";


// File: classAmesos__Merikos.xml
%feature("docstring") Amesos_Merikos "

Amesos_Merikos: A parallel divide and conquer solver.

Merikos partitions the rows of a matrix into two or more disjoint
submatrices. i.e. if rows i and j are in different submatrices, A[i,j]
== 0 == A[j,i]. Rows/columns not in any of the submatrices, i.e. the
rows/columsn of the separator, are permuted to the bottom right.

Merikos factors each of the disjoint submatrices in parallel,
(potentially by calling Amesos_Merikos() recursively), updating the
rows and columns of the separator which belong to it and forming the
schur complement of those rows and columns of the separator.

Merikos updates the trailing block of the matrix and then factors it.

Merikos is a Greek word for partial, reflecting the fact that
Amesos_Merikos uses a series of partial LU factorizations, performed
in parallel, to piece together the full LU decomposition.

C++ includes: Amesos_Merikos.h ";


// File: classAmesos__Mumps.xml
%feature("docstring") Amesos_Mumps "

Amesos_Mumps: An object-oriented wrapper for the double precision
version of MUMPS.

Amesos_Mumps is an interface to the the double precision version of
the sparse parallel direct solver MUMPS. Given an Epetra_RowMatrix A,
and two Epetra_MultiVectors X and B, the solution with Amesos_Mumps
reads as follows:

Epetra_LinearProblem Problem; Amesos_BaseSolver * Solver; Amesos
Amesos_Factory;

Solver = Amesos_Factory.Create(\"Amesos_Mumps\", Problem);

if( Solver == 0 ) cerr << \"library not available\" << endl;

Problem.SetMatrix(&A);

Solver-> SymbolicFactorization();

Solver-> NumericFactorization();

Problem.SetLHS(&X);

Problem.SetLHS(&B);

Solver-> Solve();

A number of parameters is available to tune the performances of MUMPS.
We refer to the Amesos Reference Guide for a detailed overview of
these parameters. Here, we just recall that it is possible to solve
the linear system on a subset of the processes contained in the Comm
object of the Epetra_LinearProblem.

Amesos_Mumps accepts any Epetra_RowMatrix derived class. However,
special functions are available for Epetra_CrsMatrix and
Epetra_VbrMatrix objects.

As Amesos is based on Epetra, and Epetra is only double-precision, we
still require an Epetra_LinearProblem composed by a double-precision
matrix, and two double-precision vectors. The solution vector is
casted to double after solution. Single precision may be of interest
if Amesos is used with ML, to solve the coarse problem (for which
single-precision can be enough in term of numerical error, and usually
save memory and CPU time).

Amesos_Mumps is based on Amesos_EpetraBaseSolver, that is derived from
Amesos_BaseSolver. The main redistribution utilities, as well as a
getrow function, is obtained by EpetraBaseSolver.

WARNING:  This interface is compatible with MUMPS 4.5.4.

Marzio Sala, ETHZ.

C++ includes: Amesos_Mumps.h ";

%feature("docstring")  Amesos_Mumps::MatrixShapeOK "bool
Amesos_Mumps::MatrixShapeOK() const

Returns true if the solver can handle this matrix shape.

Returns true if the matrix shape is one that the underlying sparse
direct solver can handle. Classes that work only on square matrices
should return false for rectangular matrices. Classes that work only
on symmetric matrices whould return false for non-symmetric matrices.
";

%feature("docstring")  Amesos_Mumps::Comm "const Epetra_Comm&
Amesos_Mumps::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

%feature("docstring")  Amesos_Mumps::GetProblem "const
Epetra_LinearProblem* Amesos_Mumps::GetProblem() const

Gets a pointer to the Epetra_LinearProblem. ";


// File: classAmesos__NoCopiable.xml
%feature("docstring") Amesos_NoCopiable "

Amesos_NoCopiable: Simple class to prevent the usage of copy
constructor and operator =.

Marzio Sala, SNL 9214

C++ includes: Amesos_NoCopiable.h ";

%feature("docstring")  Amesos_NoCopiable::Amesos_NoCopiable "Amesos_NoCopiable::Amesos_NoCopiable()

Default constructor. ";

%feature("docstring")  Amesos_NoCopiable::~Amesos_NoCopiable "Amesos_NoCopiable::~Amesos_NoCopiable()

Default destructor. ";


// File: classAmesos__Paraklete.xml
%feature("docstring") Amesos_Paraklete "

Interface to PARAKLETE internal solver.Interface to PARAKLETE internal
solver.

C++ includes: Amesos_Paraklete.h ";


// File: classAmesos__Paraklete__Pimpl.xml
%feature("docstring") Amesos_Paraklete_Pimpl "";


// File: classAmesos__Pardiso.xml
%feature("docstring") Amesos_Pardiso "

Amesos_Pardiso: Interface to the PARDISO package.

Marzio Sala, SNL 9214

C++ includes: Amesos_Pardiso.h ";


// File: classAmesos__Reordering.xml
%feature("docstring") Amesos_Reordering "

Amesos_Reordering: base class for reordering procedures.

Marzio Sala, ETHZ.

C++ includes: Amesos_Reordering.h ";

%feature("docstring")  Amesos_Reordering::~Amesos_Reordering "Amesos_Reordering::~Amesos_Reordering() ";

%feature("docstring")  Amesos_Reordering::GetRowPerm "virtual int*
Amesos_Reordering::GetRowPerm()=0

Returns the row permutation vector, or 0 if not computed. ";

%feature("docstring")  Amesos_Reordering::GetColPerm "virtual int*
Amesos_Reordering::GetColPerm()=0

Returns the column permutation vector, or 0 if not computed. ";


// File: classAmesos__Scalapack.xml
%feature("docstring") Amesos_Scalapack "

Amesos_Scalapack: A serial and parallel dense solver. For now, we
implement only the unsymmetric ScaLAPACK solver.

Amesos_Scalapack, an object-oriented wrapper for LAPACK and ScaLAPACK,
will solve a linear systems of equations: A X = B using Epetra objects
and the ScaLAPACK library, where A is an Epetra_RowMatrix and X and B
are Epetra_MultiVector objects.

Amesos_Scalapack can be competitive for matrices that are not
particularly sparse. ScaLAPACK solves matrices for which the fill-in
is roughly 10% to 20% of the matrix size in time comparable to that
achieve by other Amesos classes. Amesos_Scalapack scales well and
hence its performance advantage will be largest when large number of
processes are involved.

Amesos_Scalapack uses the ScaLAPACK functions PDGETRF and PDGETRS if
more than one process is used. If only one process is used,
Amesos_ScaLAPACK uses the LAPACK function PDGETRF and PDGETRS.

AmesosScaLAPACK uses full partial pivoting and will therefore provide
answers that are at least as accurate as any direct sparse solver.

AmesosScalapack makes sense under the following circumstances: There
is sufficient memory to store the entrie dense matrix. 8*n^2/p bytes
will be required on each process. -AND- one of the following The
matrix is relatively small and dense. Amesos_Scalapack will solve
matrices less than 100 by 100 faster than other Amesos classes unless
the matrices are very sparse.

The matrix is relatively dense and many processes are available. If a
thousand processes are available, Amesos_Scalapack should be
competetive with other sparse direct solvers even for matrices whose L
and U factors contain only 5% non-zeros.

The matrix is quite dense. Amesos_Scalapack will be well on any matrix
whose L and U factors contain 20% or more non-zeros.

Execution time is less important than robustness. Amesos_Scalapack is
among the most robust parallel direct solvers.

Common control parameters :

Amesos_Scalapack supports the following parameters which are common to
across multiple Amesos solvers: ParamList.set(\"MaxProcs\", int
MaximumProcessesToUse ); By default, this is set to -1, which causes
Amesos_Scalapack to use a heuristic to determine how many processes to
use. If set to a postive value, MaximumProcessesToUse,
Amesos_Scalapack will use MaximumProcessesToUse provided that there
are that many processes available. Testing should be performed with
MaximumProcessesToUse set to some value larger than one to force
parallel execution.

ParamList.set(\"PrintTiming\", bool );

ParamList.set(\"PrintStatus\", bool );

ParamList.set(\"ComputeVectorNorms\", bool );

ParamList.set(\"ComputeTrueResidual\", bool );

ParamList.set(\"OutputLevel\", int );

ParamList.set(\"DebugLevel\", int );

ParamList.set(\"ComputeTrueResidual\", bool );  Amesos_Scalapack
supports the following parameters specific to Amesos_Scalapack.
Teuchos::ParameterList ScalapackParams =
ParameterList.sublist(\"Scalapack\") ; ScalapackParams.set(\"2D
distribution\", bool ); By default this is set \"true\". In general,
because a two dimensional data distribution generally produces faster
results. However, in some cases, a one dimensional data distribution
may provide faster execution time. The code for the one dimensional
data distribution uses a different data redistribution algorithm and
uses the transpose of the matrix internally (all of which is
transparent to the user).

ScalapackParams.set(\"grid_nb\", bool ); By default this is set to 32.
On some machines, it may be possible to improve performance by up to
10% by changing the value of grid_nb. (16,24,48,64 or 128) are
reasonable values to try. For testing on small matrices, small values
of grid_nb will (if \"MaxProcs\" is set to a value greater than 1)
force the code to execute in parallel. Limitations:

None of the following limitations would be particularly difficult to
remove.

The present implementation limits the number of right hand sides to
the number of rows assigned to each process. i.e. nrhs < n/p.

The present implementation does not take advantage of symmetric or
symmetric positive definite matrices, although ScaLAPACK has separate
routines to take advantages of such matrices.

C++ includes: Amesos_Scalapack.h ";


// File: classAmesos__Scaling.xml
%feature("docstring") Amesos_Scaling "

Amesos_Scaling: base class for scaling procedures.

Marzio Sala, ETHZ.

C++ includes: Amesos_Scaling.h ";

%feature("docstring")  Amesos_Scaling::~Amesos_Scaling "Amesos_Scaling::~Amesos_Scaling() ";

%feature("docstring")  Amesos_Scaling::GetRowScaling "virtual double*
Amesos_Scaling::GetRowScaling()=0

Returns the row scaling vector, or 0 if not computed. ";

%feature("docstring")  Amesos_Scaling::GetColScaling "virtual double*
Amesos_Scaling::GetColScaling()=0

Returns the column scaling vector, or 0 if not computed. ";


// File: classAmesos__StandardIndex.xml
%feature("docstring") Amesos_StandardIndex "";

%feature("docstring")  Amesos_StandardIndex::Amesos_StandardIndex "Amesos_StandardIndex::Amesos_StandardIndex(const Epetra_Map
&OriginalMap)

Default constructor. ";

%feature("docstring")  Amesos_StandardIndex::~Amesos_StandardIndex "Amesos_StandardIndex::~Amesos_StandardIndex()

Default destructor. ";


// File: classAmesos__Status.xml
%feature("docstring") Amesos_Status "

Amesos_Status: Container for some status variables.

Marzio Sala, SNL 9214

C++ includes: Amesos_Status.h ";

%feature("docstring")  Amesos_Status::Amesos_Status "Amesos_Status::Amesos_Status()

Default constructor. ";

%feature("docstring")  Amesos_Status::~Amesos_Status "Amesos_Status::~Amesos_Status()

Default destructor. ";

%feature("docstring")  Amesos_Status::SetStatusParameters "void
Amesos_Status::SetStatusParameters(const Teuchos::ParameterList
&ParameterList) ";


// File: classAmesos__Superlu.xml
%feature("docstring") Amesos_Superlu "

Amesos_Superlu: Amesos interface to Xioye Li's SuperLU 3.0 serial
code.

Class Amesos_Superlu solves the linear systems of equations A X = B,
where A is defined as an Epetra_RowMatrix, and X and B are two
Epetra_MultiVector's.

C++ includes: Amesos_Superlu.h ";


// File: classAmesos__Superlu__Pimpl.xml
%feature("docstring") Amesos_Superlu_Pimpl "";


// File: classAmesos__Superludist.xml
%feature("docstring") Amesos_Superludist "

Amesos_Superludist: An object-oriented wrapper for Superludist.

Amesos_Superludist will solve a linear systems of equations: A X = B
using Epetra objects and the Superludist solver library, where A is an
Epetra_RowMatrix and X and B are Epetra_MultiVector objects.

C++ includes: Amesos_Superludist.h ";

%feature("docstring")  Amesos_Superludist::SetParameters "int
Amesos_Superludist::SetParameters(Teuchos::ParameterList
&ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superludist::NumSymbolicFact "int
Amesos_Superludist::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Superludist::NumNumericFact "int
Amesos_Superludist::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Superludist::NumSolve "int
Amesos_Superludist::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Superludist::PrintTiming "void
Amesos_Superludist::PrintTiming() const

Print various timig. ";

%feature("docstring")  Amesos_Superludist::PrintStatus "void
Amesos_Superludist::PrintStatus() const

Print various information about the parameters used by Superludist. ";

%feature("docstring")  Amesos_Superludist::GetTiming "void
Amesos_Superludist::GetTiming(Teuchos::ParameterList
&TimingParameterList) const

Extracts timing information from the current solver and places it in
the parameter list. ";


// File: classAmesos__Support.xml
%feature("docstring") Amesos_Support "

Amesos_Support: Collection of utilities not included in Amesos.h.

Ken Stanley

C++ includes: Amesos_Support.h ";


// File: classAmesos__Taucs.xml
%feature("docstring") Amesos_Taucs "

Interface to TAUCS.

C++ includes: Amesos_Taucs.h ";


// File: classAmesos__Taucs__Pimpl.xml
%feature("docstring") Amesos_Taucs_Pimpl "";


// File: classAmesos__TestRowMatrix.xml
%feature("docstring") Amesos_TestRowMatrix "

Amesos_TestRowMatrix: a class to test Epetra_RowMatrix based codes.

Class Amesos_TestRowMatrix enables the creation of a Epetra_RowMatrix
derived class for testing purposed. This class requires another
Epetra_RowMatrix as input, and minimic the behavior of this matrix.
However, as it this object is not derived from Epetra_CrsMatrix or
Epetra_VbrMatrix, a dynamic_cast will not result in any
Epetra_CrsMatrix or Epetra_VrbMatrix object.

Marzio Sala, SNL 9214

C++ includes: Amesos_TestRowMatrix.h ";

%feature("docstring")  Amesos_TestRowMatrix::SetUseTranspose "int
Amesos_TestRowMatrix::SetUseTranspose(bool UseTranspose)

Sets use transpose (not implemented). ";

%feature("docstring")  Amesos_TestRowMatrix::UseTranspose "bool
Amesos_TestRowMatrix::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_TestRowMatrix::HasNormInf "bool
Amesos_TestRowMatrix::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Amesos_TestRowMatrix::Comm "const Epetra_Comm&
Amesos_TestRowMatrix::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_TestRowMatrix::OperatorDomainMap "const
Epetra_Map& Amesos_TestRowMatrix::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Amesos_TestRowMatrix::OperatorRangeMap "const
Epetra_Map& Amesos_TestRowMatrix::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Amesos_TestRowMatrix::Map "const
Epetra_BlockMap& Amesos_TestRowMatrix::Map() const ";

%feature("docstring")  Amesos_TestRowMatrix::Label "const char*
Amesos_TestRowMatrix::Label() const ";


// File: classAmesos__Time.xml
%feature("docstring") Amesos_Time "

Amesos_Time: Container for timing information.

Marzio Sala, SNL 9214

C++ includes: Amesos_Time.h ";

%feature("docstring")  Amesos_Time::Amesos_Time "Amesos_Time::Amesos_Time()

Default constructor to create size timers. ";

%feature("docstring")  Amesos_Time::~Amesos_Time "virtual
Amesos_Time::~Amesos_Time()

Default destructor. ";

%feature("docstring")  Amesos_Time::CreateTimer "void
Amesos_Time::CreateTimer(const Epetra_Comm &Comm, int size=1)

Initializes the Time object. ";

%feature("docstring")  Amesos_Time::ResetTimer "void
Amesos_Time::ResetTimer(const int timerID=0)

Resets the internally stored time object. ";

%feature("docstring")  Amesos_Time::AddTime "int
Amesos_Time::AddTime(const string what, int dataID, const int
timerID=0)

Adds to field what the time elapsed since last call to ResetTimer().
";

%feature("docstring")  Amesos_Time::GetTime "double
Amesos_Time::GetTime(const string what) const

Gets the cumulative time using the string. ";

%feature("docstring")  Amesos_Time::GetTime "double
Amesos_Time::GetTime(const int dataID) const

Gets the cumulative time using the dataID. ";

%feature("docstring")  Amesos_Time::GetTiming "void
Amesos_Time::GetTiming(Teuchos::ParameterList &list) const

Load up the current timing information into the parameter list. ";


// File: structAmesos__Time__Data.xml
%feature("docstring") Amesos_Time_Data "

Amesos_Time_Data: Simple struct for storing associated data for
Amesos_Time.

Heidi Thornquist, SNL 1437

C++ includes: Amesos_Time.h ";

%feature("docstring")  Amesos_Time_Data::Amesos_Time_Data "Amesos_Time_Data::Amesos_Time_Data(string timeName, double timeVal)

Constructor. ";

%feature("docstring")  Amesos_Time_Data::~Amesos_Time_Data "virtual
Amesos_Time_Data::~Amesos_Time_Data()

Destructor. ";


// File: classAmesos__Umfpack.xml
%feature("docstring") Amesos_Umfpack "

Class Amesos_Umfpack: An object-oriented wrapper for UMFPACK.

Amesos_Umfpack will solve a linear systems of equations: A X = B using
Epetra objects and the UMFPACK solver library, where A is an
Epetra_RowMatrix and X and B are Epetra_MultiVector objects.

C++ includes: Amesos_Umfpack.h ";


// File: classAmesos__Utils.xml
%feature("docstring") Amesos_Utils "

Amesos_Utils: Collections of basic utilities.

Marzio Sala, SNL 9214

C++ includes: Amesos_Utils.h ";

%feature("docstring")  Amesos_Utils::Amesos_Utils "Amesos_Utils::Amesos_Utils()

Default constructor. ";

%feature("docstring")  Amesos_Utils::~Amesos_Utils "Amesos_Utils::~Amesos_Utils()

Default destructor. ";

%feature("docstring")  Amesos_Utils::ComputeTrueResidual "void
Amesos_Utils::ComputeTrueResidual(const Epetra_RowMatrix &Matrix,
const Epetra_MultiVector &X, const Epetra_MultiVector &B, const bool
UseTranspose, const string prefix) const

Computes the true residual, B - Matrix * X, and prints the results. ";

%feature("docstring")  Amesos_Utils::ComputeVectorNorms "void
Amesos_Utils::ComputeVectorNorms(const Epetra_MultiVector &X, const
Epetra_MultiVector &B, const string prefix) const

Computes the norms of X and B and print the results. ";

%feature("docstring")  Amesos_Utils::PrintLine "void
Amesos_Utils::PrintLine() const

Prints line on cout. ";

%feature("docstring")  Amesos_Utils::SetMaxProcesses "void
Amesos_Utils::SetMaxProcesses(int &MaxProcesses, const
Epetra_RowMatrix &A) ";


// File: structSLUData.xml
%feature("docstring") SLUData "";


// File: namespace@0.xml


// File: namespaceSLU.xml


// File: namespacestd.xml


// File: namespaceTeuchos.xml


// File: Amesos__BaseSolver_8h.xml


// File: Amesos__BTF_8h.xml


// File: Amesos__Component_8h.xml


// File: Amesos__ComponentBaseSolver_8h.xml


// File: Amesos__ConfigDefs_8h.xml


// File: Amesos__Control_8cpp.xml


// File: Amesos__Control_8h.xml


// File: Amesos__Dscpack_8cpp.xml


// File: Amesos__Dscpack_8h.xml


// File: Amesos__Klu_8cpp.xml
%feature("docstring")  deallocFunctorDeleteWithCommon "DeallocFunctorDeleteWithCommon<T,DeleteFunctor>
@0::deallocFunctorDeleteWithCommon(const RCP< klu_common > &common,
DeleteFunctor deleteFunctor) ";


// File: Amesos__Klu_8h.xml


// File: Amesos__Lapack_8cpp.xml


// File: Amesos__Lapack_8h.xml


// File: Amesos__MC64_8cpp.xml
%feature("docstring")  std::F77_FUNC "void F77_FUNC(mc64id,
MC64ID)(int *) ";

%feature("docstring")  std::F77_FUNC "void F77_FUNC(mc64ad,
MC64AD)(int * ";


// File: Amesos__MC64_8h.xml


// File: Amesos__Merikos_8h.xml


// File: Amesos__Mumps_8cpp.xml


// File: Amesos__Mumps_8h.xml


// File: Amesos__NoCopiable_8h.xml


// File: Amesos__Paraklete_8cpp.xml
%feature("docstring")  deallocFunctorDeleteWithCommon "DeallocFunctorDeleteWithCommon<T,DeleteFunctor>
@0::deallocFunctorDeleteWithCommon(const RCP< paraklete_common >
&common, DeleteFunctor deleteFunctor) ";

%feature("docstring")  my_handler "void my_handler(int status, char
*file, int line, char *msg) ";


// File: Amesos__Paraklete_8h.xml


// File: Amesos__Pardiso_8cpp.xml
%feature("docstring")  F77_PARDISOINIT "int F77_PARDISOINIT(void *,
int *, int *) ";

%feature("docstring")  F77_PARDISO "int F77_PARDISO(void *, int *,
int *, int *, int *, int *, double *, int *, int *, int *, int *, int
*, int *, double *, double *, int *) ";


// File: Amesos__Pardiso_8h.xml


// File: Amesos__Pastix_8cpp.xml


// File: Amesos__Pastix_8h.xml


// File: Amesos__Reordering_8h.xml


// File: Amesos__Scalapack_8cpp.xml
%feature("docstring")  pcolnum "int pcolnum(int j, int nb, int npcol)
";


// File: Amesos__Scalapack_8h.xml


// File: Amesos__SCALAPACK__wrappers_8h.xml
%feature("docstring")  SL_INIT_F77 "void PREFIX SL_INIT_F77(int
*blacs_context, const int *nprow, const int *npcol) ";

%feature("docstring")  DESCINIT_F77 "void PREFIX DESCINIT_F77(int
*DescA, const int *m, const int *n, const int *mblock, const int
*nblock, const int *rsrc, const int *csrc, const int *blacs_context,
const int *Lda, int *ierr) ";

%feature("docstring")  BLACS_GRIDINFO_F77 "void PREFIX
BLACS_GRIDINFO_F77(int *blacs_context, const int *nprow, const int
*npcol, const int *myrow, const int *mycol) ";

%feature("docstring")  PDGETRF_F77 "void PREFIX PDGETRF_F77(const int
*m, const int *n, double *A, const int *Ai, const int *Aj, const int
*DescA, int *ipiv, int *info) ";

%feature("docstring")  PDGETRS_F77 "void PREFIX
PDGETRS_F77(Epetra_fcd, const int *n, const int *nrhs, const double
*A, const int *Ai, const int *Aj, const int *DescA, const int *ipiv,
double *X, const int *Xi, const int *Xj, const int *DescX, int *info)
";


// File: Amesos__Scaling_8h.xml


// File: Amesos__Status_8cpp.xml


// File: Amesos__Status_8h.xml


// File: Amesos__Superlu_8cpp.xml


// File: Amesos__Superlu_8h.xml


// File: Amesos__Superludist_8cpp.xml
%feature("docstring")  Superludist_NumProcRows "int
Superludist_NumProcRows(int NumProcs) ";

%feature("docstring")  SetNPRowAndCol "int SetNPRowAndCol(const int
MaxProcesses, int &nprow, int &npcol) ";


// File: Amesos__Superludist_8h.xml


// File: Amesos__Support_8cpp.xml


// File: Amesos__Support_8h.xml


// File: Amesos__Taucs_8cpp.xml
%feature("docstring")  taucs_dccs_free_ptr "void
taucs_dccs_free_ptr(taucs_ccs_matrix **taucs_ccs_matrix_in) ";

%feature("docstring")  taucs_supernodal_factor_free_ptr "void
taucs_supernodal_factor_free_ptr(taucs_ccs_matrix
**taucs_ccs_matrix_in) ";


// File: Amesos__Taucs_8h.xml


// File: Amesos__TestRowMatrix_8h.xml


// File: Amesos__Time_8h.xml


// File: Amesos__Umfpack_8cpp.xml


// File: Amesos__Umfpack_8h.xml


// File: Amesos__Utils_8h.xml


// File: dir_3b09264d69d076c0a3b5a8b401bfd1f8.xml


// File: dir_e562a79e1d7f66fb5a29a2640e5c2e71.xml


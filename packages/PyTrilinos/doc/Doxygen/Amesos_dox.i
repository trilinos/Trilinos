
// File: index.xml

// File: classAmesos__BaseSolver.xml
%feature("docstring") Amesos_BaseSolver "

Amesos_BaseSolver: A pure virtual class for direct solution of real-
valued double- precision operators.

Pure virtual class for all Amesos concrete implementions.

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

%feature("docstring")  Amesos_BaseSolver::~Amesos_BaseSolver "virtual
Amesos_BaseSolver::~Amesos_BaseSolver()

Destructor. ";

%feature("docstring")  Amesos_BaseSolver::SymbolicFactorization "virtual int Amesos_BaseSolver::SymbolicFactorization()=0

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_BaseSolver::NumericFactorization "virtual int Amesos_BaseSolver::NumericFactorization()=0

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_BaseSolver::Solve "virtual int
Amesos_BaseSolver::Solve()=0

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_BaseSolver::SetUseTranspose "virtual
int Amesos_BaseSolver::SetUseTranspose(bool UseTranspose)=0

If set true, X will be set to the solution of AT X = B (not A X = B)

If the implementation of this interface does not support transpose
use, this method should return a value of -1.

<br >Preconditions:  SetUseTranspose() should be called prior to the
call to SymbolicFactorization() If NumericFactorization() or Solve()
is called after SetUseTranspose() without an intervening call to
SymbolicFactorization() the result is implementation dependent.

<br >Postconditions: The next factorization and solve will be
performed with the new value of UseTranspose.

Parameters:
-----------

UseTranspose:  -- (In) If true, solve AT X = B, otherwise solve A X =
B.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Amesos_BaseSolver::UseTranspose "virtual bool
Amesos_BaseSolver::UseTranspose() const =0

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_BaseSolver::SetParameters "virtual int
Amesos_BaseSolver::SetParameters(Teuchos::ParameterList
&ParameterList)=0

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_BaseSolver::GetProblem "virtual const
Epetra_LinearProblem* Amesos_BaseSolver::GetProblem() const =0

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_BaseSolver::MatrixShapeOK "virtual bool
Amesos_BaseSolver::MatrixShapeOK() const =0

Returns true if the solver can handle this matrix shape.

Returns true if the matrix shape is one that the underlying sparse
direct solver can handle. Classes that work only on square matrices
should return false for rectangular matrices. Classes that work only
on symmetric matrices whould return false for non-symmetric matrices.
";

%feature("docstring")  Amesos_BaseSolver::Comm "virtual const
Epetra_Comm& Amesos_BaseSolver::Comm() const =0

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_BaseSolver::NumSymbolicFact "virtual
int Amesos_BaseSolver::NumSymbolicFact() const =0

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_BaseSolver::NumNumericFact "virtual int
Amesos_BaseSolver::NumNumericFact() const =0

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_BaseSolver::NumSolve "virtual int
Amesos_BaseSolver::NumSolve() const =0

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_BaseSolver::PrintStatus "virtual void
Amesos_BaseSolver::PrintStatus() const =0

Prints status information about the current solver. ";

%feature("docstring")  Amesos_BaseSolver::PrintTiming "virtual void
Amesos_BaseSolver::PrintTiming() const =0

Prints timing information about the current solver. ";

%feature("docstring")  Amesos_BaseSolver::setParameterList "virtual
void Amesos_BaseSolver::setParameterList(Teuchos::RCP<
Teuchos::ParameterList > const &paramList)

Redefined from Teuchos::ParameterListAcceptor. ";

%feature("docstring")  Amesos_BaseSolver::getNonconstParameterList "virtual Teuchos::RCP<Teuchos::ParameterList>
Amesos_BaseSolver::getNonconstParameterList()

This is an empty stub. ";

%feature("docstring")  Amesos_BaseSolver::unsetParameterList "virtual
Teuchos::RCP<Teuchos::ParameterList>
Amesos_BaseSolver::unsetParameterList()

This is an empty stub. ";

%feature("docstring")  Amesos_BaseSolver::GetTiming "virtual void
Amesos_BaseSolver::GetTiming(Teuchos::ParameterList
&TimingParameterList) const

Extracts timing information from the current solver and places it in
the parameter list. ";


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

%feature("docstring")  Amesos_Btf::Amesos_Btf "Amesos_Btf::Amesos_Btf(const Epetra_LinearProblem &LinearProblem)

Amesos_Btf Constructor.

Creates an Amesos_Btf instance, using an Epetra_LinearProblem, passing
in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Btf::~Amesos_Btf "Amesos_Btf::~Amesos_Btf(void)

Amesos_Btf Destructor.

Completely deletes an Amesos_Btf object. ";

%feature("docstring")  Amesos_Btf::SymbolicFactorization "int
Amesos_Btf::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

Compute an ordering which reduces the matrix to block upper triangular
form.

Determine a task partitioning (static load balancing)

Redistribute the data based on the owner computes rule

Instantiates an Amesos solver for each of the diagonal blocks

Calls SymbolicFactorization() on each of the diagonal blocks

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Btf::NumericFactorization "int
Amesos_Btf::NumericFactorization()

Performs NumericFactorization on the matrix A.

Calls NumericFactorization() on each of the diagonal blocks

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Btf::Solve "int Amesos_Btf::Solve()

Solves A X = B (or AT X = B)

Foreach block i:    For each block j      Compute x_i -= A_{i,j} x_j
Call Solve(x_i,b_i)     Broadcast x_i

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Btf::GetProblem "const
Epetra_LinearProblem* Amesos_Btf::GetProblem() const

Get a pointer to the Problem. ";

%feature("docstring")  Amesos_Btf::MatrixShapeOK "bool
Amesos_Btf::MatrixShapeOK() const

Returns true if BTF can handle this matrix shape.

Returns true if the matrix shape is one that BTF can handle. BTF only
works with square matrices. ";

%feature("docstring")  Amesos_Btf::SetUseTranspose "int
Amesos_Btf::SetUseTranspose(bool UseTranspose)

SetUseTranpose(true) causes Solve() To compute A^T X = B.

If SetUseTranspose() is set to true,

AT X = B is computed

else

A X = B is computed ";

%feature("docstring")  Amesos_Btf::UseTranspose "bool
Amesos_Btf::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Btf::Comm "const Epetra_Comm&
Amesos_Btf::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

%feature("docstring")  Amesos_Btf::SetParameters "int
Amesos_Btf::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subsequent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Btf::NumSymbolicFact "int
Amesos_Btf::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Btf::NumNumericFact "int
Amesos_Btf::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Btf::NumSolve "int
Amesos_Btf::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Btf::PrintTiming "void
Amesos_Btf::PrintTiming()

Print timing information. ";

%feature("docstring")  Amesos_Btf::PrintStatus "void
Amesos_Btf::PrintStatus()

Print information about the factorization and solution phases. ";


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

Preconditions:  An ordering  Postconditions:  Constructor requirements

Every Amesos_SolverName class should accept an Epetra_LinearProblem

C++ includes: Amesos_Component.h ";

%feature("docstring")  Amesos_Component::~Amesos_Component "virtual
Amesos_Component::~Amesos_Component()

Destructor. ";

%feature("docstring")  Amesos_Component::PartialFactorization "virtual int Amesos_Component::PartialFactorization()=0

Performs partial factorization on the matrix A.

Partial Factorization perfom

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Component::Lsolve "virtual int
Amesos_Component::Lsolve()=0

Solves L X = B (or LT x = B)

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Component::Usolve "virtual int
Amesos_Component::Usolve()=0

Solves L X = B (or LT x = B)

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Component::SetRowPermutation "virtual
int Amesos_Component::SetRowPermutation(int *RowPermutation)=0

SetRowPermutation. ";

%feature("docstring")  Amesos_Component::SetColumnPermutation "virtual int Amesos_Component::SetColumnPermutation(int
*ColumnPermutation)=0

SetColumnPermutation. ";

%feature("docstring")  Amesos_Component::SetSubMatrixSize "virtual
int Amesos_Component::SetSubMatrixSize(int SubMatrixSize)=0

SetSubMatrixSize. ";

%feature("docstring")  Amesos_Component::GetRowPermutation "virtual
int Amesos_Component::GetRowPermutation(int **RowPermutation)=0

GetRowPermutation.

RowPermutation reflects any row permutations performed by
PartialFactorization(). Note: It is not yet clear whether this row
permutation includes the RowPermuation upon input or whether it
returns only the row permuations performed by the most recent call to
PartialFactorization(). In other words, in the absence of pivoting,
RowPermutation might be identical to that given by SetRowPermutation()
or it might be the identity permutation. ";

%feature("docstring")  Amesos_Component::GetColumnPermutation "virtual int Amesos_Component::GetColumnPermutation(int
**ColumnPermutation)=0

GetColumnPermutation.

ColumnPermutation reflects any row permutations performed by
PartialFactorization(). Note: It is not yet clear whether this row
permutation includes the ColumnPermuation upon input or whether it
returns only the row permuations performed by the most recent call to
PartialFactorization(). In other words, in the absence of pivoting,
ColumnPermutation might be identical to that given by
SetColumnPermutation() or it might be the identity permutation. ";

%feature("docstring")  Amesos_Component::GetSubMatrixSize "virtual
int Amesos_Component::GetSubMatrixSize(int *SubMatrixSize)=0

GetSubMatrixSize. ";

%feature("docstring")  Amesos_Component::GetSchurComplement "virtual
int Amesos_Component::GetSchurComplement(Epetra_CrsMatrix
*SchurComplement)=0

GetSchurComplement. ";


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

Preconditions:  An ordering  Postconditions:  Constructor requirements

Every Amesos_SolverName class should accept an Epetra_LinearProblem

C++ includes: Amesos_ComponentBaseSolver.h ";

%feature("docstring")
Amesos_ComponentBaseSolver::~Amesos_ComponentBaseSolver "virtual
Amesos_ComponentBaseSolver::~Amesos_ComponentBaseSolver()

Destructor. ";

%feature("docstring")
Amesos_ComponentBaseSolver::PartialFactorization "virtual int
Amesos_ComponentBaseSolver::PartialFactorization()=0

Performs partial factorization on the matrix A.

Partial Factorization perfom

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_ComponentBaseSolver::Lsolve "virtual
int Amesos_ComponentBaseSolver::Lsolve()=0

Solves L X = B (or LT x = B)

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_ComponentBaseSolver::LsolveStart "*
virtual int Amesos_ComponentBaseSolver::LsolveStart()=0

Solves the triangular part of L X1 = B (or LT x = B)

Integer error code, set to 0 if successful, -1 if unimplimented. ";

%feature("docstring")  Amesos_ComponentBaseSolver::LsolvePart "virtual int Amesos_ComponentBaseSolver::LsolvePart(int begin, int
end)=0

Computes L[begin..end,:] X1.

Integer error code, set to 0 if successful, -1 if unimplimented. ";

%feature("docstring")  Amesos_ComponentBaseSolver::Usolve "virtual
int Amesos_ComponentBaseSolver::Usolve()=0

Solves U X = B (or UT x = B)

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_ComponentBaseSolver::UsolveStart "*
virtual int Amesos_ComponentBaseSolver::UsolveStart()=0

Solves the triangular part of U X1 = B (or LT x = B)

Integer error code, set to 0 if successful, -1 if unimplimented. ";

%feature("docstring")  Amesos_ComponentBaseSolver::UsolvePart "virtual int Amesos_ComponentBaseSolver::UsolvePart(int begin, int
end)=0

Computes U[:,begin..end] X1.

Integer error code, set to 0 if successful, -1 if unimplimented. ";

%feature("docstring")  Amesos_ComponentBaseSolver::SetRowPermutation "virtual int Amesos_ComponentBaseSolver::SetRowPermutation(int
*RowPermutation)=0

Solves U X = B (or UT x = B)

Integer error code, set to 0 if successful. SetRowPermutation ";

%feature("docstring")
Amesos_ComponentBaseSolver::SetColumnPermutation "virtual int
Amesos_ComponentBaseSolver::SetColumnPermutation(int
*ColumnPermutation)=0

SetColumnPermutation. ";

%feature("docstring")  Amesos_ComponentBaseSolver::SetSubMatrixSize "virtual int Amesos_ComponentBaseSolver::SetSubMatrixSize(int
SubMatrixSize)=0

SetSubMatrixSize. ";

%feature("docstring")  Amesos_ComponentBaseSolver::GetRowPermutation "virtual int Amesos_ComponentBaseSolver::GetRowPermutation(int
**RowPermutation)=0

GetRowPermutation.

RowPermutation reflects any row permutations performed by
PartialFactorization(). Note: It is not yet clear whether this row
permutation includes the RowPermuation upon input or whether it
returns only the row permuations performed by the most recent call to
PartialFactorization(). In other words, in the absence of pivoting,
RowPermutation might be identical to that given by SetRowPermutation()
or it might be the identity permutation. ";

%feature("docstring")
Amesos_ComponentBaseSolver::GetColumnPermutation "virtual int
Amesos_ComponentBaseSolver::GetColumnPermutation(int
**ColumnPermutation)=0

GetColumnPermutation.

ColumnPermutation reflects any row permutations performed by
PartialFactorization(). Note: It is not yet clear whether this row
permutation includes the ColumnPermuation upon input or whether it
returns only the row permuations performed by the most recent call to
PartialFactorization(). In other words, in the absence of pivoting,
ColumnPermutation might be identical to that given by
SetColumnPermutation() or it might be the identity permutation. ";

%feature("docstring")  Amesos_ComponentBaseSolver::GetSubMatrixSize "virtual int Amesos_ComponentBaseSolver::GetSubMatrixSize(int
*SubMatrixSize)=0

GetSubMatrixSize. ";

%feature("docstring")  Amesos_ComponentBaseSolver::GetSchurComplement
"virtual int
Amesos_ComponentBaseSolver::GetSchurComplement(Epetra_CrsMatrix
*SchurComplement)=0

GetSchurComplement. ";


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

%feature("docstring")  Amesos_Dscpack::Amesos_Dscpack "Amesos_Dscpack::Amesos_Dscpack(const Epetra_LinearProblem
&LinearProblem)

Amesos_Dscpack Constructor.

Creates an Amesos_Dscpack instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Dscpack::~Amesos_Dscpack "Amesos_Dscpack::~Amesos_Dscpack(void)

Amesos_Dscpack Destructor.

Completely deletes an Amesos_Dscpack object. ";

%feature("docstring")  Amesos_Dscpack::SymbolicFactorization "int
Amesos_Dscpack::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Dscpack::NumericFactorization "int
Amesos_Dscpack::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Dscpack::Solve "int
Amesos_Dscpack::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Dscpack::GetProblem "const
Epetra_LinearProblem* Amesos_Dscpack::GetProblem() const

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_Dscpack::MatrixShapeOK "bool
Amesos_Dscpack::MatrixShapeOK() const

Returns true if DSCPACK can handle this matrix shape.

Returns true if the matrix shape is one that DSCPACK can handle.
DSCPACK only works with symetric matrices. ";

%feature("docstring")  Amesos_Dscpack::SetUseTranspose "int
Amesos_Dscpack::SetUseTranspose(bool UseTranspose)

If set true, X will be set to the solution of AT X = B (not A X = B)

If the implementation of this interface does not support transpose
use, this method should return a value of -1.

<br >Preconditions:  SetUseTranspose() should be called prior to the
call to SymbolicFactorization() If NumericFactorization() or Solve()
is called after SetUseTranspose() without an intervening call to
SymbolicFactorization() the result is implementation dependent.

<br >Postconditions: The next factorization and solve will be
performed with the new value of UseTranspose.

Parameters:
-----------

UseTranspose:  -- (In) If true, solve AT X = B, otherwise solve A X =
B.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Amesos_Dscpack::UseTranspose "bool
Amesos_Dscpack::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Dscpack::Comm "const Epetra_Comm&
Amesos_Dscpack::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Dscpack::SetParameters "int
Amesos_Dscpack::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Dscpack::NumSymbolicFact "int
Amesos_Dscpack::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Dscpack::NumNumericFact "int
Amesos_Dscpack::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Dscpack::NumSolve "int
Amesos_Dscpack::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Dscpack::PrintTiming "void
Amesos_Dscpack::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Dscpack::PrintStatus "void
Amesos_Dscpack::PrintStatus() const

Prints information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Dscpack::GetTiming "void
Amesos_Dscpack::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";


// File: classAmesos__Dscpack__Pimpl.xml
%feature("docstring") Amesos_Dscpack_Pimpl "";


// File: classAmesos__Klu.xml
%feature("docstring") Amesos_Klu "

Interface to KLU internal solver.

Interface to UMFPACK.

C++ includes: Amesos_Umfpack.h ";

%feature("docstring")  Amesos_Klu::Amesos_Klu "Amesos_Klu::Amesos_Klu(const Epetra_LinearProblem &LinearProblem)

Amesos_Klu Constructor.

Creates an Amesos_Klu instance, using an Epetra_LinearProblem, passing
in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Klu::~Amesos_Klu "Amesos_Klu::~Amesos_Klu(void)

Amesos_Klu Destructor. ";

%feature("docstring")  Amesos_Klu::SymbolicFactorization "int
Amesos_Klu::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Klu::NumericFactorization "int
Amesos_Klu::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Klu::Solve "int Amesos_Klu::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Klu::GetProblem "const
Epetra_LinearProblem* Amesos_Klu::GetProblem() const

Get a pointer to the Problem. ";

%feature("docstring")  Amesos_Klu::MatrixShapeOK "bool
Amesos_Klu::MatrixShapeOK() const

Returns true if KLU can handle this matrix shape.

Returns true if the matrix shape is one that KLU can handle. KLU only
works with square matrices. ";

%feature("docstring")  Amesos_Klu::SetUseTranspose "int
Amesos_Klu::SetUseTranspose(bool UseTranspose_in)

SetUseTranpose(true) is more efficient in Amesos_Klu.

If SetUseTranspose() is set to true, $A^T X = B$ is computed. ";

%feature("docstring")  Amesos_Klu::UseTranspose "bool
Amesos_Klu::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Klu::Comm "const Epetra_Comm&
Amesos_Klu::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Klu::SetParameters "int
Amesos_Klu::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Klu::NumSymbolicFact "int
Amesos_Klu::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Klu::NumNumericFact "int
Amesos_Klu::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Klu::NumSolve "int
Amesos_Klu::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Klu::PrintTiming "void
Amesos_Klu::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Klu::PrintStatus "void
Amesos_Klu::PrintStatus() const

Prints information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Klu::GetTiming "void
Amesos_Klu::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information and places in parameter list. ";


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

%feature("docstring")  Amesos_Lapack::Amesos_Lapack "Amesos_Lapack::Amesos_Lapack(const Epetra_LinearProblem
&LinearProblem)

Amesos_Lapack Constructor.

Creates an Amesos_Lapack instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Lapack::~Amesos_Lapack "Amesos_Lapack::~Amesos_Lapack(void)

Amesos_Lapack Destructor.

Completely deletes an Amesos_Lapack object. ";

%feature("docstring")  Amesos_Lapack::SymbolicFactorization "int
Amesos_Lapack::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Lapack::NumericFactorization "int
Amesos_Lapack::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Lapack::Solve "int
Amesos_Lapack::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Lapack::GetProblem "const
Epetra_LinearProblem* Amesos_Lapack::GetProblem() const

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_Lapack::MatrixShapeOK "bool
Amesos_Lapack::MatrixShapeOK() const

Returns true if the solver can handle this matrix shape.

Returns true if the matrix shape is one that the underlying sparse
direct solver can handle. Classes that work only on square matrices
should return false for rectangular matrices. Classes that work only
on symmetric matrices whould return false for non-symmetric matrices.
";

%feature("docstring")  Amesos_Lapack::SetUseTranspose "int
Amesos_Lapack::SetUseTranspose(bool UseTranspose_in)

If set true, X will be set to the solution of AT X = B (not A X = B)

If the implementation of this interface does not support transpose
use, this method should return a value of -1.

<br >Preconditions:  SetUseTranspose() should be called prior to the
call to SymbolicFactorization() If NumericFactorization() or Solve()
is called after SetUseTranspose() without an intervening call to
SymbolicFactorization() the result is implementation dependent.

<br >Postconditions: The next factorization and solve will be
performed with the new value of UseTranspose.

Parameters:
-----------

UseTranspose:  -- (In) If true, solve AT X = B, otherwise solve A X =
B.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Amesos_Lapack::UseTranspose "bool
Amesos_Lapack::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Lapack::Comm "const Epetra_Comm&
Amesos_Lapack::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Lapack::setParameterList "void
Amesos_Lapack::setParameterList(Teuchos::RCP< Teuchos::ParameterList >
const &paramList)

Use this parameter list to read values from.

Redefined from Teuchos::ParameterListAcceptor ";

%feature("docstring")  Amesos_Lapack::unsetParameterList "Teuchos::RCP< Teuchos::ParameterList >
Amesos_Lapack::unsetParameterList()

This is an empty stub. ";

%feature("docstring")  Amesos_Lapack::SetParameters "int
Amesos_Lapack::SetParameters(Teuchos::ParameterList &ParameterList)

Deprecated - Sets parameters. ";

%feature("docstring")  Amesos_Lapack::GEEV "int
Amesos_Lapack::GEEV(Epetra_Vector &Er, Epetra_Vector &Ei)

Computes the eigenvalues of the linear system matrix using DGEEV.

Parameters:
-----------

Er:  - (Out) On processor zero only, it will contain the real
component of the eigenvalues.

Ei:  - (Out) On processor zero only, it will contain the imaginary
component of the eigenvalues.

Er and Ei must have been allocated so that the local length on
processor 0 equals the global size of the matrix. ";

%feature("docstring")  Amesos_Lapack::NumSymbolicFact "int
Amesos_Lapack::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Lapack::NumNumericFact "int
Amesos_Lapack::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Lapack::NumSolve "int
Amesos_Lapack::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Lapack::PrintTiming "void
Amesos_Lapack::PrintTiming() const

Print timing information. ";

%feature("docstring")  Amesos_Lapack::PrintStatus "void
Amesos_Lapack::PrintStatus() const

Print information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Lapack::GetTiming "void
Amesos_Lapack::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";


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

%feature("docstring")  Amesos_Merikos::Amesos_Merikos "Amesos_Merikos::Amesos_Merikos(const Epetra_LinearProblem
&LinearProblem)

Amesos_Merikos Constructor.

Creates an Amesos_Merikos instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object. ";

%feature("docstring")  Amesos_Merikos::~Amesos_Merikos "Amesos_Merikos::~Amesos_Merikos(void)

Amesos_Merikos Destructor.

Completely deletes an Amesos_Merikos object. ";

%feature("docstring")  Amesos_Merikos::RedistributeA "int
Amesos_Merikos::RedistributeA()

Performs SymbolicFactorization on the matrix A.

SymbolicFactorization() takes no action in Amesos_Merikos().

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::ConvertToScalapack "int
Amesos_Merikos::ConvertToScalapack() ";

%feature("docstring")  Amesos_Merikos::PerformNumericFactorization "int Amesos_Merikos::PerformNumericFactorization() ";

%feature("docstring")  Amesos_Merikos::SymbolicFactorization "int
Amesos_Merikos::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::NumericFactorization "int
Amesos_Merikos::NumericFactorization()

Performs NumericFactorization on the matrix A.

Static pivoting (i.e. scale and permute the matrix to produce a zero-
free diagonal and to minimize the need for pivoting later).
Partition the matrix      Redistribute the matrix to match the
partitioning     Foreach subblock of the matrix do: Note:  this will
happen in parallel       Create an instance of an Amesos solver object
(must           support the Amesos_Component interface)       Call
PartialFactorization        Add the Schur Complement into the trailing
block of the matrix.     Endfor Create an Amesos instance for the
trailing block of the matrix. Call SymbolicFactorization on the
trailing block      Call NumericFactorization on the trailing block

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::LSolve "int
Amesos_Merikos::LSolve()

Solves L X = B.

| L11   0   0  |  X1      B1      | L21 L22   0  |  X2   =  B2 | L31
L32 L33  |  X3   =  B3

Foreach subblock of the matrix do:       Note:  this will happen in
parallel       Lsolve()          i.e. L11.Solve(X1, B1) and
L22.Solve(X2, B2)        Update the elements of B corresponding to the
seperator,         i.e. B3 = B3 - L31 X1 - L32 X2      Endfor Perform
a solve on the trailing matrix:       i.e. L33.LSolve(X3,B3)

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::USolve "int
Amesos_Merikos::USolve()

Solves U X = B.

| U11 U12 U13  |  X1      B1      |   0 U22 U23  |  X2   =  B2 |   0
0 U33  |  X3   =  B3

Perform a solve on the trailing matrix:       i.e. U33.USolve(X3,B3)
Foreach subblock of the matrix do:       Note: this will happen in
parallel       Update the elements of B corresponding to this block
i.e. B2 = B2 - U23 X3 ; B1 = B1 - U13 X3        Usolve()          i.e.
U11.Solve(X1, B1) and U22.Solve(X2, B2)      Endfor

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::Solve "int
Amesos_Merikos::Solve()

Solves A X = B.

| L11     U12     U13  |  X1      B1      | L21     L22     U23 |  X2
=  B2      | L31     L32         A33  |  X3   =  B3

Foreach subblock of the matrix do:       Note:  this will happen in
parallel       Lsolve()          i.e. L11.Solve(X1, B1) and
L22.Solve(X2, B2)        Update the elements of B corresponding to the
seperator,         i.e. B3 = B3 - L31 X1 - L32 X2      Endfor Perform
a solve on the trailing matrix:       i.e. A33.Solve(X3,B3)

B = X ;     Foreach subblock of the matrix do:       Note:  this will
happen in parallel       Update the elements of B corresponding to
this block         i.e. B2 = B2 - U23 X3 ; B1 = B1 - U13 X3 Usolve()
i.e. U11.Solve(X1, B1) and U22.Solve(X2, B2) Endfor

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::GetProblem "const
Epetra_LinearProblem* Amesos_Merikos::GetProblem() const

Get a pointer to the Problem. ";

%feature("docstring")  Amesos_Merikos::MatrixShapeOK "bool
Amesos_Merikos::MatrixShapeOK() const

Returns true if MERIKOS can handle this matrix shape.

Returns true if the matrix shape is one that MERIKOS can handle.
MERIKOS only works with square matrices. ";

%feature("docstring")  Amesos_Merikos::SetUseTranspose "int
Amesos_Merikos::SetUseTranspose(bool UseTranspose)

SetUseTranpose() controls whether to compute AX=B or ATX = B. ";

%feature("docstring")  Amesos_Merikos::UseTranspose "bool
Amesos_Merikos::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Merikos::Comm "const Epetra_Comm&
Amesos_Merikos::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

%feature("docstring")  Amesos_Merikos::SetParameters "int
Amesos_Merikos::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Merikos::NumSymbolicFact "int
Amesos_Merikos::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Merikos::NumNumericFact "int
Amesos_Merikos::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Merikos::NumSolve "int
Amesos_Merikos::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Merikos::PrintTiming "void
Amesos_Merikos::PrintTiming()

Print timing information. ";

%feature("docstring")  Amesos_Merikos::PrintStatus "void
Amesos_Merikos::PrintStatus()

Print information about the factorization and solution phases. ";


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

if( Solver == 0 ) std::cerr << \"library not available\" << std::endl;

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

%feature("docstring")  Amesos_Mumps::Amesos_Mumps "Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &LinearProblem)

Amesos_Mumps Constructor.

Creates an Amesos_Mumps instance, using an Epetra_LinearProblem, ";

%feature("docstring")  Amesos_Mumps::~Amesos_Mumps "Amesos_Mumps::~Amesos_Mumps(void)

Amesos_Mumps Destructor.

Deletes an Amesos_Mumps object. ";

%feature("docstring")  Amesos_Mumps::SymbolicFactorization "int
Amesos_Mumps::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Mumps::NumericFactorization "int
Amesos_Mumps::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Mumps::Solve "int Amesos_Mumps::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Mumps::Destroy "void
Amesos_Mumps::Destroy()

Destroys all data associated with  this object. ";

%feature("docstring")  Amesos_Mumps::SetUseTranspose "int
Amesos_Mumps::SetUseTranspose(bool UseTranspose_in)

If set true, X will be set to the solution of AT X = B (not A X = B)

If the implementation of this interface does not support transpose
use, this method should return a value of -1.

<br >Preconditions:  SetUseTranspose() should be called prior to the
call to SymbolicFactorization() If NumericFactorization() or Solve()
is called after SetUseTranspose() without an intervening call to
SymbolicFactorization() the result is implementation dependent.

<br >Postconditions: The next factorization and solve will be
performed with the new value of UseTranspose.

Parameters:
-----------

UseTranspose:  -- (In) If true, solve AT X = B, otherwise solve A X =
B.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Amesos_Mumps::UseTranspose "bool
Amesos_Mumps::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Mumps::SetParameters "int
Amesos_Mumps::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Mumps::NumSymbolicFact "int
Amesos_Mumps::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Mumps::NumNumericFact "int
Amesos_Mumps::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Mumps::NumSolve "int
Amesos_Mumps::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Mumps::PrintTiming "void
Amesos_Mumps::PrintTiming() const

Prints timing information.

In the destruction phase, prints out detailed information about the
various phases: symbolic and numeric factorization, solution,
gather/scatter for vectors and matrices. ";

%feature("docstring")  Amesos_Mumps::PrintStatus "void
Amesos_Mumps::PrintStatus() const

Prints information about the factorization and solution phases.

In the destruction phase, prints out some information furnished by
MUMPS, like the amount of required memory, the MFLOPS. ";

%feature("docstring")  Amesos_Mumps::GetTiming "void
Amesos_Mumps::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";

%feature("docstring")  Amesos_Mumps::SetPrecscaling "int
Amesos_Mumps::SetPrecscaling(double *ColSca, double *RowSca)

Set prescaling.

Use double precision vectors of size N (global dimension of the
matrix) as scaling for columns and rows. ColSca and RowSca must be
defined on the host only, and allocated by the user, if the user sets
ICNTL(8) = -1.

Both input vectors are float with --enable-amesos-smumps, double
otherwise. ";

%feature("docstring")  Amesos_Mumps::SetRowScaling "int
Amesos_Mumps::SetRowScaling(double *RowSca)

Set row scaling.

Use double precision vectors of size N (global dimension of the
matrix) for row scaling.

Parameters:
-----------

RowSca:  (In) - float pointer with --enable-amesos-smumps, double
pointer otherwise. ";

%feature("docstring")  Amesos_Mumps::SetColScaling "int
Amesos_Mumps::SetColScaling(double *ColSca)

Set column scaling.

Use double precision vectors of size N (global dimension of the
matrix) for column scaling.

Parameters:
-----------

ColSca:  (In) - float pointer with --enable-amesos-smumps, double
pointer otherwise. ";

%feature("docstring")  Amesos_Mumps::SetOrdering "int
Amesos_Mumps::SetOrdering(int *PermIn)

Sets ordering.

Use integer vectors of size N (global dimension of the matrix) as
given ordering. PermIn must be defined on the host only, and allocated
by the user, if the user sets ICNTL(7) = 1. ";

%feature("docstring")  Amesos_Mumps::GetRINFO "double *
Amesos_Mumps::GetRINFO()

Gets the pointer to the RINFO array (defined on all processes).

Gets the pointer to the internally stored RINFO array, of type float
if option --enable-amesos-smumps is enabled, double otherwise. ";

%feature("docstring")  Amesos_Mumps::GetINFO "int *
Amesos_Mumps::GetINFO()

Gets the pointer to the INFO array (defined on all processes).

Gets the pointer to the internally stored INFO array, of type int. ";

%feature("docstring")  Amesos_Mumps::GetRINFOG "double *
Amesos_Mumps::GetRINFOG()

Gets the pointer to the RINFOG array (defined on host only).

Gets the pointer to the internally stored RINFOG array (defined on the
host process only), of type float if option --enable-amesos-smumps is
enabled, double otherwise. ";

%feature("docstring")  Amesos_Mumps::GetINFOG "int *
Amesos_Mumps::GetINFOG()

Get the pointer to the INFOG array (defined on host only).

Gets the pointer to the internally stored INFOG (defined on the host
process only) array, of type int. ";

%feature("docstring")  Amesos_Mumps::SetICNTL "void
Amesos_Mumps::SetICNTL(int pos, int value)

Set ICNTL[pos] to value. pos is expressed in FORTRAN style (starting
from 1). ";

%feature("docstring")  Amesos_Mumps::SetCNTL "void
Amesos_Mumps::SetCNTL(int pos, double value)

Set CNTL[pos] to value. pos is expressed in FORTRAN style (starting
from 1). ";

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

%feature("docstring")  Amesos_Paraklete::Amesos_Paraklete "Amesos_Paraklete::Amesos_Paraklete(const Epetra_LinearProblem
&LinearProblem)

Amesos_Paraklete Constructor.

Creates an Amesos_Paraklete instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Paraklete::~Amesos_Paraklete "Amesos_Paraklete::~Amesos_Paraklete(void)

Amesos_Paraklete Destructor. ";

%feature("docstring")  Amesos_Paraklete::SymbolicFactorization "int
Amesos_Paraklete::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Paraklete::NumericFactorization "int
Amesos_Paraklete::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Paraklete::Solve "int
Amesos_Paraklete::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Paraklete::GetProblem "const
Epetra_LinearProblem* Amesos_Paraklete::GetProblem() const

Get a pointer to the Problem. ";

%feature("docstring")  Amesos_Paraklete::MatrixShapeOK "bool
Amesos_Paraklete::MatrixShapeOK() const

Returns true if PARAKLETE can handle this matrix shape.

Returns true if the matrix shape is one that PARAKLETE can handle.
PARAKLETE only works with square matrices. ";

%feature("docstring")  Amesos_Paraklete::SetUseTranspose "int
Amesos_Paraklete::SetUseTranspose(bool UseTranspose_in)

SetUseTranpose()

If SetUseTranspose() is set to true, $A^T X = B$ is computed. ";

%feature("docstring")  Amesos_Paraklete::UseTranspose "bool
Amesos_Paraklete::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Paraklete::Comm "const Epetra_Comm&
Amesos_Paraklete::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Paraklete::SetParameters "int
Amesos_Paraklete::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Paraklete::NumSymbolicFact "int
Amesos_Paraklete::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Paraklete::NumNumericFact "int
Amesos_Paraklete::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Paraklete::NumSolve "int
Amesos_Paraklete::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Paraklete::PrintTiming "void
Amesos_Paraklete::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Paraklete::PrintStatus "void
Amesos_Paraklete::PrintStatus() const

Prints information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Paraklete::GetTiming "void
Amesos_Paraklete::GetTiming(Teuchos::ParameterList
&TimingParameterList) const

Extracts timing information from the current solver and places it in
the parameter list. ";


// File: classAmesos__Paraklete__Pimpl.xml
%feature("docstring") Amesos_Paraklete_Pimpl "";


// File: classAmesos__Pardiso.xml
%feature("docstring") Amesos_Pardiso "

Amesos_Pardiso: Interface to the PARDISO package.

Marzio Sala, SNL 9214

C++ includes: Amesos_Pardiso.h ";

%feature("docstring")  Amesos_Pardiso::Amesos_Pardiso "Amesos_Pardiso::Amesos_Pardiso(const Epetra_LinearProblem
&LinearProblem)

Constructor. ";

%feature("docstring")  Amesos_Pardiso::~Amesos_Pardiso "Amesos_Pardiso::~Amesos_Pardiso()

Destructor. ";

%feature("docstring")  Amesos_Pardiso::SymbolicFactorization "int
Amesos_Pardiso::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A. ";

%feature("docstring")  Amesos_Pardiso::NumericFactorization "int
Amesos_Pardiso::NumericFactorization()

Performs NumericFactorization on the matrix A. ";

%feature("docstring")  Amesos_Pardiso::Solve "int
Amesos_Pardiso::Solve()

Solves A X = B (or AT X = B) ";

%feature("docstring")  Amesos_Pardiso::GetProblem "const
Epetra_LinearProblem* Amesos_Pardiso::GetProblem() const

Get a pointer to the Problem. ";

%feature("docstring")  Amesos_Pardiso::MatrixShapeOK "bool
Amesos_Pardiso::MatrixShapeOK() const

Returns true if PARDISO can handle this matrix shape.

Returns true if the matrix shape is one that PARDISO can handle.
PARDISO only works with square matrices. ";

%feature("docstring")  Amesos_Pardiso::SetUseTranspose "int
Amesos_Pardiso::SetUseTranspose(bool UseTranspose)

SetUseTranpose()

If SetUseTranspose() is set to true, $A^T X = B$ is computed. ";

%feature("docstring")  Amesos_Pardiso::UseTranspose "bool
Amesos_Pardiso::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Pardiso::Comm "const Epetra_Comm&
Amesos_Pardiso::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

%feature("docstring")  Amesos_Pardiso::SetParameters "int
Amesos_Pardiso::SetParameters(Teuchos::ParameterList &ParameterList)

Set parameters from the input parameters list, returns 0 if
successful. ";

%feature("docstring")  Amesos_Pardiso::NumSymbolicFact "int
Amesos_Pardiso::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Pardiso::NumNumericFact "int
Amesos_Pardiso::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Pardiso::NumSolve "int
Amesos_Pardiso::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Pardiso::PrintTiming "void
Amesos_Pardiso::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Pardiso::PrintStatus "void
Amesos_Pardiso::PrintStatus() const

Prints information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Pardiso::GetTiming "void
Amesos_Pardiso::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";


// File: classAmesos__Reordering.xml
%feature("docstring") Amesos_Reordering "

Amesos_Reordering: base class for reordering procedures.

Base class for reordering procedures.

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

%feature("docstring")  Amesos_Scalapack::Amesos_Scalapack "Amesos_Scalapack::Amesos_Scalapack(const Epetra_LinearProblem
&LinearProblem)

Amesos_Scalapack Constructor.

Creates an Amesos_Scalapack instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Scalapack::~Amesos_Scalapack "Amesos_Scalapack::~Amesos_Scalapack(void)

Amesos_Scalapack Destructor.

Completely deletes an Amesos_Scalapack object. ";

%feature("docstring")  Amesos_Scalapack::SymbolicFactorization "int
Amesos_Scalapack::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

There is no symbolic factorization phase in ScaLAPACK, as it operates
only on dense matrices. Hence,
Amesos_Scalapack::SymbolicFactorization() takes no action.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Scalapack::NumericFactorization "int
Amesos_Scalapack::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6) NOT
IMPLEMENTED

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). Irrelevant for Amesos_Scalapack.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization(). Irrelevant for Amesos_Scalapack.

postconditions: nprow_, npcol_, DescA_

DenseA will be factored

Ipiv_ contains the pivots

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Scalapack::Solve "int
Amesos_Scalapack::Solve()

Solves A X = B (or AT X = B)

preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6) NOT
IMPLEMENTED

X and B must have the same shape (NOT CHECKED)

X and B must have fewer than nb right hand sides. EPETRA_CHK_ERR(-2)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The matrix should not have changed since the last call to
NumericFactorization().

postconditions: X will be set such that A X = B (or AT X = B), within
the limits of the accuracy of the the Scalapack solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Scalapack::GetProblem "const
Epetra_LinearProblem* Amesos_Scalapack::GetProblem() const

Get a pointer to the Problem. ";

%feature("docstring")  Amesos_Scalapack::MatrixShapeOK "bool
Amesos_Scalapack::MatrixShapeOK() const

Returns true if SCALAPACK can handle this matrix shape.

Returns true if the matrix shape is one that SCALAPACK can handle.
SCALAPACK only works with square matrices. ";

%feature("docstring")  Amesos_Scalapack::SetUseTranspose "int
Amesos_Scalapack::SetUseTranspose(bool UseTranspose)

SetUseTranpose(true) is more efficient in Amesos_Scalapack.

If SetUseTranspose() is set to true,

AT X = B is computed

else

A X = B is computed ";

%feature("docstring")  Amesos_Scalapack::UseTranspose "bool
Amesos_Scalapack::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Scalapack::Comm "const Epetra_Comm&
Amesos_Scalapack::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

%feature("docstring")  Amesos_Scalapack::SetParameters "int
Amesos_Scalapack::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subsequent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Scalapack::NumSymbolicFact "int
Amesos_Scalapack::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Scalapack::NumNumericFact "int
Amesos_Scalapack::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Scalapack::NumSolve "int
Amesos_Scalapack::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Scalapack::PrintTiming "void
Amesos_Scalapack::PrintTiming() const

Print timing information. ";

%feature("docstring")  Amesos_Scalapack::PrintStatus "void
Amesos_Scalapack::PrintStatus() const

Print information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Scalapack::GetTiming "void
Amesos_Scalapack::GetTiming(Teuchos::ParameterList
&TimingParameterList) const

Extracts timing information from the current solver and places it in
the parameter list. ";


// File: classAmesos__Scaling.xml
%feature("docstring") Amesos_Scaling "

Amesos_Scaling: base class for scaling procedures.

Base class for scaling procedures.

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

%feature("docstring")  Amesos_Superlu::Amesos_Superlu "Amesos_Superlu::Amesos_Superlu(const Epetra_LinearProblem
&LinearProblem)

Amesos_Superlu Constructor.

Creates an Amesos_Superlu instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Superlu::~Amesos_Superlu "Amesos_Superlu::~Amesos_Superlu()

Amesos_Superlu Destructor. ";

%feature("docstring")  Amesos_Superlu::SymbolicFactorization "int
Amesos_Superlu::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superlu::NumericFactorization "int
Amesos_Superlu::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superlu::Solve "int
Amesos_Superlu::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superlu::GetProblem "const
Epetra_LinearProblem* Amesos_Superlu::GetProblem() const

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_Superlu::MatrixShapeOK "bool
Amesos_Superlu::MatrixShapeOK() const

Returns true if the solver can handle this matrix shape.

Returns true if the matrix shape is one that the underlying sparse
direct solver can handle. Classes that work only on square matrices
should return false for rectangular matrices. Classes that work only
on symmetric matrices whould return false for non-symmetric matrices.
";

%feature("docstring")  Amesos_Superlu::SetUseTranspose "int
Amesos_Superlu::SetUseTranspose(bool UseTranspose)

If set true, X will be set to the solution of AT X = B (not A X = B)

If the implementation of this interface does not support transpose
use, this method should return a value of -1.

<br >Preconditions:  SetUseTranspose() should be called prior to the
call to SymbolicFactorization() If NumericFactorization() or Solve()
is called after SetUseTranspose() without an intervening call to
SymbolicFactorization() the result is implementation dependent.

<br >Postconditions: The next factorization and solve will be
performed with the new value of UseTranspose.

Parameters:
-----------

UseTranspose:  -- (In) If true, solve AT X = B, otherwise solve A X =
B.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Amesos_Superlu::UseTranspose "bool
Amesos_Superlu::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Superlu::Comm "const Epetra_Comm&
Amesos_Superlu::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Superlu::SetParameters "int
Amesos_Superlu::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superlu::NumSymbolicFact "int
Amesos_Superlu::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Superlu::NumNumericFact "int
Amesos_Superlu::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Superlu::NumSolve "int
Amesos_Superlu::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Superlu::PrintTiming "void
Amesos_Superlu::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Superlu::PrintStatus "void
Amesos_Superlu::PrintStatus() const

Prints status information. ";

%feature("docstring")  Amesos_Superlu::GetTiming "void
Amesos_Superlu::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";


// File: classAmesos__Superlu__Pimpl.xml
%feature("docstring") Amesos_Superlu_Pimpl "";


// File: classAmesos__Superludist.xml
%feature("docstring") Amesos_Superludist "

Amesos_Superludist: An object-oriented wrapper for Superludist.

Amesos_Superludist will solve a linear systems of equations: A X = B
using Epetra objects and the Superludist solver library, where A is an
Epetra_RowMatrix and X and B are Epetra_MultiVector objects.

C++ includes: Amesos_Superludist.h ";

%feature("docstring")  Amesos_Superludist::Amesos_Superludist "Amesos_Superludist::Amesos_Superludist(const Epetra_LinearProblem
&LinearProblem)

Amesos_Superludist Constructor.

Creates an Amesos_Superludist instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Superludist::~Amesos_Superludist "Amesos_Superludist::~Amesos_Superludist(void)

Amesos_Superludist Destructor.

Completely deletes an Amesos_Superludist object. ";

%feature("docstring")  Amesos_Superludist::SymbolicFactorization "int
Amesos_Superludist::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superludist::NumericFactorization "int
Amesos_Superludist::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superludist::Solve "int
Amesos_Superludist::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Superludist::SetUseTranspose "int
Amesos_Superludist::SetUseTranspose(bool UseTranspose)

Amesos_Superludist does not support transpose at this time.

returns 0 if UseTranspose is set to false, else 1 (failure) ";

%feature("docstring")  Amesos_Superludist::GetProblem "const
Epetra_LinearProblem* Amesos_Superludist::GetProblem() const

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_Superludist::MatrixShapeOK "bool
Amesos_Superludist::MatrixShapeOK() const

Returns true if SUPERLUDIST can handle this matrix shape.

Returns true if the matrix shape is one that SUPERLUDIST can handle.
SUPERLUDIST only works with square matrices. ";

%feature("docstring")  Amesos_Superludist::UseTranspose "bool
Amesos_Superludist::UseTranspose() const

Always returns true. ";

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

%feature("docstring")  Amesos_Taucs::Amesos_Taucs "Amesos_Taucs::Amesos_Taucs(const Epetra_LinearProblem &LinearProblem)

Default constructor. ";

%feature("docstring")  Amesos_Taucs::~Amesos_Taucs "Amesos_Taucs::~Amesos_Taucs(void)

Default destructor. ";

%feature("docstring")  Amesos_Taucs::SymbolicFactorization "int
Amesos_Taucs::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Taucs::NumericFactorization "int
Amesos_Taucs::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Taucs::Solve "int Amesos_Taucs::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Taucs::GetProblem "const
Epetra_LinearProblem* Amesos_Taucs::GetProblem() const

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_Taucs::MatrixShapeOK "bool
Amesos_Taucs::MatrixShapeOK() const

Returns true if the solver can handle this matrix shape.

Returns true if the matrix shape is one that the underlying sparse
direct solver can handle. Classes that work only on square matrices
should return false for rectangular matrices. Classes that work only
on symmetric matrices whould return false for non-symmetric matrices.
";

%feature("docstring")  Amesos_Taucs::SetUseTranspose "int
Amesos_Taucs::SetUseTranspose(bool UseTranspose)

Amesos_Taucs supports only symmetric matrices, hence transpose is
irrelevant, but harmless. ";

%feature("docstring")  Amesos_Taucs::UseTranspose "bool
Amesos_Taucs::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Taucs::Comm "const Epetra_Comm&
Amesos_Taucs::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Taucs::SetParameters "int
Amesos_Taucs::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Taucs::NumSymbolicFact "int
Amesos_Taucs::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Taucs::NumNumericFact "int
Amesos_Taucs::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Taucs::NumSolve "int
Amesos_Taucs::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Taucs::PrintTiming "void
Amesos_Taucs::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Taucs::PrintStatus "void
Amesos_Taucs::PrintStatus() const

Prints status information. ";

%feature("docstring")  Amesos_Taucs::GetTiming "void
Amesos_Taucs::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";


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

%feature("docstring")  Amesos_TestRowMatrix::Amesos_TestRowMatrix "Amesos_TestRowMatrix::Amesos_TestRowMatrix(Epetra_RowMatrix
*Matrix_in)

Constructor. ";

%feature("docstring")  Amesos_TestRowMatrix::~Amesos_TestRowMatrix "virtual Amesos_TestRowMatrix::~Amesos_TestRowMatrix()

Destructor. ";

%feature("docstring")  Amesos_TestRowMatrix::NumMyRowEntries "virtual
int Amesos_TestRowMatrix::NumMyRowEntries(int MyRow, int &NumEntries)
const

Returns the number of nonzero entries in MyRow.

Parameters:
-----------

MyRow:  - (In) Local row.

NumEntries:  - (Out) Number of nonzero values present.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_TestRowMatrix::MaxNumEntries "virtual
int Amesos_TestRowMatrix::MaxNumEntries() const

Returns the maximum of NumMyRowEntries() over all rows. ";

%feature("docstring")  Amesos_TestRowMatrix::ExtractMyRowCopy "virtual int Amesos_TestRowMatrix::ExtractMyRowCopy(int MyRow, int
Length, int &NumEntries, double *Values, int *Indices) const

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

MyRow:  - (In) Local row to extract.

Length:  - (In) Length of Values and Indices.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Indices:  - (Out) Extracted global column indices for the
corresponding values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_TestRowMatrix::ExtractDiagonalCopy "virtual int Amesos_TestRowMatrix::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Diagonal:  - (Out) Extracted main diagonal.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_TestRowMatrix::Multiply "virtual int
Amesos_TestRowMatrix::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_RowMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

TransA:  -(In) If true, multiply by the transpose of matrix, otherwise
just use matrix.

X:  - (In) A Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  -(Out) A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_TestRowMatrix::Solve "virtual int
Amesos_TestRowMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns result of a local-only solve using a triangular
Epetra_RowMatrix with Epetra_MultiVectors X and Y (NOT IMPLEMENTED).
";

%feature("docstring")  Amesos_TestRowMatrix::Apply "virtual int
Amesos_TestRowMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Amesos_TestRowMatrix::ApplyInverse "virtual
int Amesos_TestRowMatrix::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Amesos_TestRowMatrix::InvRowSums "virtual int
Amesos_TestRowMatrix::InvRowSums(Epetra_Vector &x) const

Computes the sum of absolute values of the rows of the
Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED). ";

%feature("docstring")  Amesos_TestRowMatrix::LeftScale "virtual int
Amesos_TestRowMatrix::LeftScale(const Epetra_Vector &x)

Scales the Epetra_RowMatrix on the left with a Epetra_Vector x (NOT
IMPLEMENTED). ";

%feature("docstring")  Amesos_TestRowMatrix::InvColSums "virtual int
Amesos_TestRowMatrix::InvColSums(Epetra_Vector &x) const

Computes the sum of absolute values of the columns of the
Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED). ";

%feature("docstring")  Amesos_TestRowMatrix::RightScale "virtual int
Amesos_TestRowMatrix::RightScale(const Epetra_Vector &x)

Scales the Epetra_RowMatrix on the right with a Epetra_Vector x (NOT
IMPLEMENTED). ";

%feature("docstring")  Amesos_TestRowMatrix::Filled "virtual bool
Amesos_TestRowMatrix::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  Amesos_TestRowMatrix::NormInf "virtual double
Amesos_TestRowMatrix::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Amesos_TestRowMatrix::NormOne "virtual double
Amesos_TestRowMatrix::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")  Amesos_TestRowMatrix::NumGlobalNonzeros "virtual int Amesos_TestRowMatrix::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Amesos_TestRowMatrix::NumGlobalRows "virtual
int Amesos_TestRowMatrix::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Amesos_TestRowMatrix::NumGlobalCols "virtual
int Amesos_TestRowMatrix::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Amesos_TestRowMatrix::NumGlobalDiagonals "virtual int Amesos_TestRowMatrix::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Amesos_TestRowMatrix::NumMyNonzeros "virtual
int Amesos_TestRowMatrix::NumMyNonzeros() const

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Amesos_TestRowMatrix::NumMyRows "virtual int
Amesos_TestRowMatrix::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Amesos_TestRowMatrix::NumMyCols "virtual int
Amesos_TestRowMatrix::NumMyCols() const

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")  Amesos_TestRowMatrix::NumMyDiagonals "virtual
int Amesos_TestRowMatrix::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Amesos_TestRowMatrix::LowerTriangular "virtual
bool Amesos_TestRowMatrix::LowerTriangular() const

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Amesos_TestRowMatrix::UpperTriangular "virtual
bool Amesos_TestRowMatrix::UpperTriangular() const

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Amesos_TestRowMatrix::RowMatrixRowMap "virtual
const Epetra_Map& Amesos_TestRowMatrix::RowMatrixRowMap() const

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")  Amesos_TestRowMatrix::RowMatrixColMap "virtual
const Epetra_Map& Amesos_TestRowMatrix::RowMatrixColMap() const

Returns the Epetra_Map object associated with the columns of this
matrix. ";

%feature("docstring")  Amesos_TestRowMatrix::RowMatrixImporter "virtual const Epetra_Import* Amesos_TestRowMatrix::RowMatrixImporter()
const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")  Amesos_TestRowMatrix::SetUseTranspose "int
Amesos_TestRowMatrix::SetUseTranspose(bool UseTranspose_in)

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
Amesos_Time::AddTime(const std::string what, int dataID, const int
timerID=0)

Adds to field what the time elapsed since last call to ResetTimer().
";

%feature("docstring")  Amesos_Time::GetTime "double
Amesos_Time::GetTime(const std::string what) const

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

%feature("docstring")  Amesos_Time_Data::Amesos_Time_Data "Amesos_Time_Data::Amesos_Time_Data(std::string timeName, double
timeVal)

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

%feature("docstring")  Amesos_Umfpack::Amesos_Umfpack "Amesos_Umfpack::Amesos_Umfpack(const Epetra_LinearProblem
&LinearProblem)

Amesos_Umfpack Constructor.

Creates an Amesos_Umfpack instance, using an Epetra_LinearProblem,
passing in an already- defined Epetra_LinearProblem object.

Note: The operator in LinearProblem must be an Epetra_RowMatrix. ";

%feature("docstring")  Amesos_Umfpack::~Amesos_Umfpack "Amesos_Umfpack::~Amesos_Umfpack(void)

Amesos_Umfpack Destructor.

Completely deletes an Amesos_Umfpack object. ";

%feature("docstring")  Amesos_Umfpack::SymbolicFactorization "int
Amesos_Umfpack::SymbolicFactorization()

Performs SymbolicFactorization on the matrix A.

In addition to performing symbolic factorization on the matrix A, the
call to SymbolicFactorization() implies that no change will be made to
the non-zero structure of the underlying matrix without a subsequent
call to SymbolicFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

<br >Postconditions: Symbolic Factorization will be performed (or
marked to be performed) allowing NumericFactorization() and Solve() to
be called.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Umfpack::NumericFactorization "int
Amesos_Umfpack::NumericFactorization()

Performs NumericFactorization on the matrix A.

In addition to performing numeric factorization on the matrix A, the
call to NumericFactorization() implies that no change will be made to
the underlying matrix without a subsequent call to
NumericFactorization().

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization(). (return -2 if the number of non-
zeros changes) Other changes can have arbitrary consequences.

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should be indexed from 0 to n-1, unless the parameter
\"Reindex\" was set to \"true\" prior to the call to
SymbolicFactorization(). (return -3 - if caught)

The paremeter \"Reindex\" should not be set to \"true\" except on
CrsMatrices. (return -4)

The paremeter \"Reindex\" should not be set to \"true\" unless Amesos
was built with EpetraExt, i.e. with --enable-epetraext on the
configure line. (return -4)

Internal errors retur -5.

<br >Postconditions: Numeric Factorization will be performed (or
marked to be performed) allowing Solve() to be performed correctly
despite a potential change in in the matrix values (though not in the
non-zero structure).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Umfpack::Solve "int
Amesos_Umfpack::Solve()

Solves A X = B (or AT x = B)

<br >Preconditions:  GetProblem().GetOperator() != 0 (return -1)

MatrixShapeOk( GetProblem().GetOperator()) == true (return -6)

GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for
return values)

The non-zero structure of the matrix should not have changed since the
last call to SymbolicFactorization().

The distribution of the matrix should not have changed since the last
call to SymbolicFactorization()

The matrix should not have changed since the last call to
NumericFactorization().

<br >Postconditions: X will be set such that A X = B (or AT X = B),
within the limits of the accuracy of the underlying solver.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Umfpack::GetProblem "const
Epetra_LinearProblem* Amesos_Umfpack::GetProblem() const

Returns the Epetra_LinearProblem.

Warning! Do not call return->SetOperator(...) to attempt to change the
Epetra_Operator object (even if the new matrix has the same
structure). This new operator matrix will be ignored! ";

%feature("docstring")  Amesos_Umfpack::MatrixShapeOK "bool
Amesos_Umfpack::MatrixShapeOK() const

Returns true if UMFPACK can handle this matrix shape.

Returns true if the matrix shape is one that UMFPACK can handle.
UMFPACK only works with square matrices. ";

%feature("docstring")  Amesos_Umfpack::SetUseTranspose "int
Amesos_Umfpack::SetUseTranspose(bool UseTranspose_in)

If set true, X will be set to the solution of AT X = B (not A X = B)

If the implementation of this interface does not support transpose
use, this method should return a value of -1.

<br >Preconditions:  SetUseTranspose() should be called prior to the
call to SymbolicFactorization() If NumericFactorization() or Solve()
is called after SetUseTranspose() without an intervening call to
SymbolicFactorization() the result is implementation dependent.

<br >Postconditions: The next factorization and solve will be
performed with the new value of UseTranspose.

Parameters:
-----------

UseTranspose:  -- (In) If true, solve AT X = B, otherwise solve A X =
B.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Amesos_Umfpack::UseTranspose "bool
Amesos_Umfpack::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Amesos_Umfpack::Comm "const Epetra_Comm&
Amesos_Umfpack::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Amesos_Umfpack::GetRcond "double
Amesos_Umfpack::GetRcond() const

Returns an estimate of the reciprocal of the condition number. ";

%feature("docstring")  Amesos_Umfpack::SetParameters "int
Amesos_Umfpack::SetParameters(Teuchos::ParameterList &ParameterList)

Updates internal variables.

<br >Preconditions: None.

<br >Postconditions: Internal variables controlling the factorization
and solve will be updated and take effect on all subseuent calls to
NumericFactorization() and Solve().

All parameters whose value are to differ from the default values must
be included in ParameterList. Parameters not specified in
ParameterList revert to their default values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Amesos_Umfpack::NumSymbolicFact "int
Amesos_Umfpack::NumSymbolicFact() const

Returns the number of symbolic factorizations performed by this
object. ";

%feature("docstring")  Amesos_Umfpack::NumNumericFact "int
Amesos_Umfpack::NumNumericFact() const

Returns the number of numeric factorizations performed by this object.
";

%feature("docstring")  Amesos_Umfpack::NumSolve "int
Amesos_Umfpack::NumSolve() const

Returns the number of solves performed by this object. ";

%feature("docstring")  Amesos_Umfpack::PrintTiming "void
Amesos_Umfpack::PrintTiming() const

Prints timing information. ";

%feature("docstring")  Amesos_Umfpack::PrintStatus "void
Amesos_Umfpack::PrintStatus() const

Prints information about the factorization and solution phases. ";

%feature("docstring")  Amesos_Umfpack::GetTiming "void
Amesos_Umfpack::GetTiming(Teuchos::ParameterList &TimingParameterList)
const

Extracts timing information from the current solver and places it in
the parameter list. ";


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
UseTranspose, const std::string prefix) const

Computes the true residual, B - Matrix * X, and prints the results. ";

%feature("docstring")  Amesos_Utils::ComputeVectorNorms "void
Amesos_Utils::ComputeVectorNorms(const Epetra_MultiVector &X, const
Epetra_MultiVector &B, const std::string prefix) const

Computes the norms of X and B and print the results. ";

%feature("docstring")  Amesos_Utils::PrintLine "void
Amesos_Utils::PrintLine() const

Prints line on std::cout. ";

%feature("docstring")  Amesos_Utils::SetMaxProcesses "void
Amesos_Utils::SetMaxProcesses(int &MaxProcesses, const
Epetra_RowMatrix &A) ";


// File: structSLUData.xml
%feature("docstring") SLUData "";


// File: namespace@19.xml


// File: namespace@9.xml


// File: namespaceSLU.xml


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


// File: Amesos__Klu_8h.xml


// File: Amesos__Lapack_8cpp.xml


// File: Amesos__Lapack_8h.xml


// File: Amesos__MC64_8cpp.xml
%feature("docstring")  F77_FUNC "void F77_FUNC(mc64id, MC64ID)(int *)
";

%feature("docstring")  F77_FUNC "void F77_FUNC(mc64ad, MC64AD)(int *
";


// File: Amesos__MC64_8h.xml


// File: Amesos__Merikos_8h.xml


// File: Amesos__Mumps_8cpp.xml


// File: Amesos__Mumps_8h.xml


// File: Amesos__NoCopiable_8h.xml


// File: Amesos__Paraklete_8cpp.xml
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


// File: dir_25619facb131199a9ac01607470cb94c.xml


// File: dir_0b7ccc221d1448fd2902bb281f2b4eff.xml


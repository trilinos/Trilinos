
// File: index.xml

// File: classNOX_1_1LineSearch_1_1Backtrack.xml
%feature("docstring") NOX::LineSearch::Backtrack "

Generic backtracking line search.

This line search starts with the step length defined by \"Default
Step\". It checks to see if the norm of the right hand side (RHS) has
been reduced. If so, it exits successfully. Otherwise, it reduces the
step length by the reduction factor (defaults to one-half). It
continues to repeat this procedure until it either finds a reduction
in the norm of the RHS or the step is less than that specified by
\"Minimum Step\". In the later case, the line search has failed, and
we take the step defined by \"Recovery Step\".

This line search can be called via NOX::LineSearch::Manager.

The following parameters can be specified for this line search in the
\"Backtrack\" sublist of the \"Line Search\" sublist.

\"Default Step\" - starting step length (defaults to 1.0)

\"Minimum Step\" - minimum acceptable step length (defaults to
1.0e-12)

\"Recovery Step\" - step to take when the line search fails (defaults
to value for \"Default Step\")

\"Max Iters\" - maximum number of iterations (i.e., RHS computations)

\"Reduction Factor\" - A multiplier between zero and one that reduces
the step size between line search iterations

C++ includes: NOX_LineSearch_Backtrack.H ";

%feature("docstring")  NOX::LineSearch::Backtrack::Backtrack "NOX::LineSearch::Backtrack::Backtrack(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::LineSearch::Backtrack::~Backtrack "NOX::LineSearch::Backtrack::~Backtrack()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Backtrack::reset "bool
NOX::LineSearch::Backtrack::reset(const Teuchos::RCP< NOX::GlobalData
> &gd, Teuchos::ParameterList &params) ";

%feature("docstring")  NOX::LineSearch::Backtrack::compute "bool
NOX::LineSearch::Backtrack::compute(NOX::Abstract::Group &newgrp,
double &step, const NOX::Abstract::Vector &dir, const
NOX::Solver::Generic &s)

Perform a line search.

On input:

Parameters:
-----------

grp:  The initial solution vector, $x_{\\\\rm old}$.

dir:  A vector of directions to be used in the line search, $d$.

s:  The nonlinear solver.

On output:

Parameters:
-----------

step:  The distance the direction was scaled, $ \\\\lambda $.

grp:  The grp is updated with a new solution, $ x_{\\\\rm new} $,
resulting from the linesearch. Normally, for a single direction line
search, this is computed as:

\\\\[ x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d. \\\\]

Ideally, $ \\\\|F(x_{\\\\rm new})\\\\| < \\\\|F(x_{\\\\rm old})\\\\| $
(e.g the final direction is a descent direction).

Note that the dir object is a std::vector. For typical line searches
as described in the above equation, this vector is of size one. We
have used a std::vector to allow for special cases of multi-
directional line searches such as the Bader/Schnabel curvillinear line
search.

Return value is true for a successful line search computation. ";


// File: classNOX_1_1Direction_1_1Broyden.xml
%feature("docstring") NOX::Direction::Broyden "

Broyden direction

We will calculate a limited-memory Broyden direction of the form

$d_k = -B_k^{-1} F_k.$

Here $B_k$ is a limited-memory Broyden approximation to the Jacobian
of $F$ at $x_k$, and $F_k = F(x_k)$. It is based on apply Broyden
updates to the Jacobian from some previous step.

The Broyden direction can only be used with
NOX::Solver::LineSearchBased. It cannot be used with any other solver,
include NOX::Solver::TrustRegionBased.  References

C. T. Kelley, Iterative Methods for Linear and Nonlinear Equations,
SIAM, 1995.

Parameters

To use this direction, specify that the \"Method\" is \"Broyden\" in
the \"Direction\" sublist of the parameters that are passed to the
solver (see NOX::Direction::Manager for more information on choosing
the search direction).

In \"Direction\"/\"Broyden\":

\"Restart Frequency\" - How often the Jacobian should be refreshed. A
value of 5, for example, means that the Jacobian should be updated
every 5 iterations. Defaults to 10.

\"Max Convergence Rate\" - Maximum convergence rate allowed when
reusing the Jacobian. The Jacobian will be refreshed if the
convergence rate, $ \\\\alpha $, is larger than this value. The
convergence rate is calculated by $ \\\\alpha = \\\\frac{\\\\| F_k
\\\\| }{\\\\| F_{k-1} \\\\|} $ where F is the nonlinear residual and $
k $ is the nonlinear iteration. Defaults to 1.0.

\"Memory\" - The maximum number of past updates that can be saved in
memory. Defaults to the value of \"Restart Frequency\".

\"Linear Solver\" - optional SUBLIST of linear solver parameters.

\"Linear Solver\"/\"Tolerance\" - Desired tolerance for linear solve.
Defaults to 1.0e-4. The tolerance can be computed using adaptive
forcing terms. See NOX::Direction::Utils::InexactNewton for additional
options.

C++ includes: NOX_Direction_Broyden.H ";

%feature("docstring")  NOX::Direction::Broyden::Broyden "NOX::Direction::Broyden::Broyden(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::Direction::Broyden::~Broyden "NOX::Direction::Broyden::~Broyden()

Destructor. ";

%feature("docstring")  NOX::Direction::Broyden::reset "bool
NOX::Direction::Broyden::reset(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params)

Reset direction based on possibly new parameters. ";

%feature("docstring")  NOX::Direction::Broyden::compute "bool
NOX::Direction::Broyden::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::Generic &solver)

Not supported for this direction - only works for line search based
solver. ";

%feature("docstring")  NOX::Direction::Broyden::compute "bool
NOX::Direction::Broyden::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::LineSearchBased &solver)

Same as compute( NOX::Abstract::Vector&, NOX::Abstract::Group&, const
NOX::Solver::Generic&).

Enables direct support for line search based solvers for the purpose
of efficiency since the LineSearchBased object has a getStep()
function that some directions require.

If it is not redefined in the derived class, it will just call the
compute with the NOX::Solver::Generic argument. ";


// File: classNOX_1_1Direction_1_1Broyden_1_1BroydenMemory.xml
%feature("docstring") NOX::Direction::Broyden::BroydenMemory "

Utility class for NOX::Direction::Broyden method to manage the
information stored in \"limited\" memory.

Store up to $m$ MemoryUnit objects where $m$ is passed to reset().
Every time push() is called, a new MemoryUnit is added. If there are
already $m$ MemoryUnit's, the oldest is bumped off the list. The zero-
th entry is always the oldest.

In order to avoid allocating and deallocating memory, we reuse the
BroydenMemoryUnit objects rather than destroying and re-constructing
them. However, this detail should be transparent to users of this
class.

C++ includes: NOX_Direction_Broyden.H ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemory::BroydenMemory "NOX::Direction::Broyden::BroydenMemory::BroydenMemory()

Constructor.

Does nothing. ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemory::~BroydenMemory "NOX::Direction::Broyden::BroydenMemory::~BroydenMemory()

Destructor.

Does nothing. ";

%feature("docstring")  NOX::Direction::Broyden::BroydenMemory::reset "void NOX::Direction::Broyden::BroydenMemory::reset(int m)

Reset the memory.

Sets mMax to the value of m. Sets the size of the memory vector to be
at least mMax. Sets the capacity of the index vector to be at least
mMax. Sets the size of the index vector to be zero. ";

%feature("docstring")  NOX::Direction::Broyden::BroydenMemory::reset "void NOX::Direction::Broyden::BroydenMemory::reset()

Reset the memory.

Sets the size of the index vector to be zero. ";

%feature("docstring")  NOX::Direction::Broyden::BroydenMemory::push "void NOX::Direction::Broyden::BroydenMemory::push(const
NOX::Abstract::Vector &d)

Add new information to the memory.

We need to calculate where the new udpate should be stored in memory
and update the information in index.

Let k denote the index of where the new update should be stored. If
there are current m items stored in memory and m < mMax, then we set k
= m. Otherwise, we set k equal to the location of the oldest update.
The oldest update is deleted to make room for the new update. In both
cases, index must be updated appropriately so that the first (zero)
entry points to the oldest update and the last entry points to the
newest update. ";

%feature("docstring")  NOX::Direction::Broyden::BroydenMemory::empty "bool NOX::Direction::Broyden::BroydenMemory::empty() const

Returns true if the memory is empty. ";

%feature("docstring")  NOX::Direction::Broyden::BroydenMemory::size "int NOX::Direction::Broyden::BroydenMemory::size() const

Number of items in memory. ";


// File: classNOX_1_1Direction_1_1Broyden_1_1BroydenMemoryUnit.xml
%feature("docstring") NOX::Direction::Broyden::BroydenMemoryUnit "

Utility class for NOX::Direction::Broyden::BroydenMemory.

Stores an $s$-vector and associated information for the limited-memory
Broyden update.

Throughout the docuementation of this class, we make the following
associations. sptr stores the vector $s$

lambda stores the value of $\\\\lambda$

snormsqr stores the values of $\\\\|s\\\\|_2^2$

C++ includes: NOX_Direction_Broyden.H ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::BroydenMemoryUnit "NOX::Direction::Broyden::BroydenMemoryUnit::BroydenMemoryUnit()

Constructor.

Sets sptr to NULL, and lambda and snormsqr to zero. ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::~BroydenMemoryUnit "NOX::Direction::Broyden::BroydenMemoryUnit::~BroydenMemoryUnit()

Destuctor.

Deletes sptr. ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::reset "void
NOX::Direction::Broyden::BroydenMemoryUnit::reset(const
NOX::Abstract::Vector &d)

Reset this memory unit with a new update vector.

Let the vector $d$ represent the input argument. Then we set $s=d$.
Also recalculates $\\\\|s\\\\|_2^2$ and sets $\\\\lambda = 0$.

If sptr is NULL, space is allocated by cloning the input vector (see
NOX::Abstract::Vector::clone). ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::setStep "void
NOX::Direction::Broyden::BroydenMemoryUnit::setStep(double step)

Update the step length.

Set $ \\\\lambda $ to the input argument. If $ \\\\lambda \\\\neq 1 $,
then reset $ s = \\\\lambda s $ and $ \\\\|s\\\\|_2^2 = \\\\lambda^2
\\\\|s\\\\|_2^2 $. ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::sPtr "Teuchos::RCP< const
NOX::Abstract::Vector >
NOX::Direction::Broyden::BroydenMemoryUnit::sPtr() const

Get pointer to $s$. ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::step "double
NOX::Direction::Broyden::BroydenMemoryUnit::step() const

Get the step, $\\\\lambda$. ";

%feature("docstring")
NOX::Direction::Broyden::BroydenMemoryUnit::sNormSqr "double
NOX::Direction::Broyden::BroydenMemoryUnit::sNormSqr() const

Get $\\\\|s\\\\|_2^2 $. ";


// File: classNOX_1_1Epetra_1_1BroydenOperator.xml
%feature("docstring") NOX::Epetra::BroydenOperator "

A concrete implementation of a Broyden-type operator for NOX.

This operator is intended to allow cheap updates to an existing
Jacobian or preconditioning matrix that would otherwise be difficult
or impossible to obtain by other means. It computes updates using
secant approximations emobdied in truncated Broyden updates that
preserve matrix sparsity.

This class derives from NOX::Abstract::PrePostOperator in order to
perform a Broyden-type update on an existing matrix that it holds but
does not own. This update is performed after each nonlinear iteration
within method runPostIterate(...) according to the recursive formula:

\\\\[ \\\\tilde{B}_{k+1} = \\\\tilde{B}_k + \\\\frac{({y_k -
\\\\tilde{B}_k s_k})s_k^T}{s^T s} \\\\]

where \\\\[ y_k = F_{k+1} - F_k \\\\] and \\\\[ s_k = x_{k+1} - x_k
\\\\]

The tilde on the matrices $ B $ indicates that the updates are
constrained so that the nonzero structure of the original matrix
passed into the constructor is preserved. Inasmuch as unconstrained
Broyden updates produce dense matrices, these constrained updates lead
to a loss of Broyden-matrix properties, e.g.

\\\\[ \\\\tilde{B}_{k+1} s_k \\\\ne \\\\tilde{B}_k + s_k \\\\]

\\\\[ \\\\tilde{B}_{k+1} q \\\\ne \\\\tilde{B}_k q \\\\quad \\\\forall
q : s_k^T q = 0 \\\\]

One could recover these properties by passing into the constructor a
dense Epetra_CrsMatrix, though the cost of typical use of this matrix,
e.g. applying ILU to it, would be significant. Additionally,
\"better\" values obtained from another Jacobian or preconditioning
matrix can be used to replace corresponding values in the updated
Broyden matrix by passing the Jacobian or preconditioning matrix and
its associated interface to the constructor. The structure of the
Jacobain or preconditioning matrix typically represents a subset of
the Broyden matrix, e.g. a block diagonal matrix.

C++ includes: NOX_Epetra_BroydenOperator.H ";

/*  IsValid flags  */

/* True if objects are current with respect to the currect stepVec.

*/

/*  "Is" functions  */

/* Checks to see if various objects have been computed. Returns true
if the corresponding \"compute\" function has been called since the
last update to the solution vector (via instantiation or computeX).

*/

%feature("docstring")  NOX::Epetra::BroydenOperator::isStep "bool
BroydenOperator::isStep() const ";

%feature("docstring")  NOX::Epetra::BroydenOperator::isYield "bool
BroydenOperator::isYield() const ";

%feature("docstring")  NOX::Epetra::BroydenOperator::isBroyden "bool
BroydenOperator::isBroyden() const ";

%feature("docstring")  NOX::Epetra::BroydenOperator::BroydenOperator "BroydenOperator::BroydenOperator(Teuchos::ParameterList &nlParams,
const Teuchos::RCP< NOX::Utils > &utils, Epetra_Vector &solnVec, const
Teuchos::RCP< Epetra_CrsMatrix > &broydMat0, bool verbose=false)

Constructor taking an initial matrix to be updated. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::BroydenOperator "BroydenOperator::BroydenOperator(const BroydenOperator &)

Copy Constructor. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::~BroydenOperator
"BroydenOperator::~BroydenOperator()

Destructor. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Label "const
char * BroydenOperator::Label() const

Returns a character string describing the name of the operator. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::SetUseTranspose "int BroydenOperator::SetUseTranspose(bool UseTranspose)

If set true, the transpose of this operator will be applied. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Apply "int
BroydenOperator::Apply(const Epetra_MultiVector &X, Epetra_MultiVector
&Y) const

Return the result on an Epetra_Operator applied to an
Epetra_MultiVector X in Y. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::ApplyInverse "int BroydenOperator::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Return the result on an Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::UseTranspose "bool BroydenOperator::UseTranspose() const

Returns the current use transpose setting. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::HasNormInf "bool
BroydenOperator::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::OperatorDomainMap
"const Epetra_Map & BroydenOperator::OperatorDomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::OperatorRangeMap
"const Epetra_Map & BroydenOperator::OperatorRangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Filled "bool
BroydenOperator::Filled() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumMyRowEntries "int BroydenOperator::NumMyRowEntries(int MyRow, int &NumEntries) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::MaxNumEntries "int BroydenOperator::MaxNumEntries() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::ExtractMyRowCopy
"int BroydenOperator::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")
NOX::Epetra::BroydenOperator::ExtractDiagonalCopy "int
BroydenOperator::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Multiply "int
BroydenOperator::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Solve "int
BroydenOperator::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::InvRowSums "int
BroydenOperator::InvRowSums(Epetra_Vector &x) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::LeftScale "int
BroydenOperator::LeftScale(const Epetra_Vector &x)

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::InvColSums "int
BroydenOperator::InvColSums(Epetra_Vector &x) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::RightScale "int
BroydenOperator::RightScale(const Epetra_Vector &x)

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NormInf "double
BroydenOperator::NormInf() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NormOne "double
BroydenOperator::NormOne() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumGlobalNonzeros
"int BroydenOperator::NumGlobalNonzeros() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumGlobalRows "int BroydenOperator::NumGlobalRows() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumGlobalCols "int BroydenOperator::NumGlobalCols() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")
NOX::Epetra::BroydenOperator::NumGlobalDiagonals "int
BroydenOperator::NumGlobalDiagonals() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumMyNonzeros "int BroydenOperator::NumMyNonzeros() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumMyRows "int
BroydenOperator::NumMyRows() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumMyCols "int
BroydenOperator::NumMyCols() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::NumMyDiagonals "int BroydenOperator::NumMyDiagonals() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::LowerTriangular "bool BroydenOperator::LowerTriangular() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::UpperTriangular "bool BroydenOperator::UpperTriangular() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Comm "const
Epetra_Comm & BroydenOperator::Comm() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::RowMatrixRowMap "const Epetra_Map & BroydenOperator::RowMatrixRowMap() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::RowMatrixColMap "const Epetra_Map & BroydenOperator::RowMatrixColMap() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::RowMatrixImporter
"const Epetra_Import * BroydenOperator::RowMatrixImporter() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::Map "const
Epetra_BlockMap & BroydenOperator::Map() const

See Epetra_SrcDistObj documentation. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::setStepVector "void BroydenOperator::setStepVector(Epetra_Vector &vec)

Set the current step vector, \\\\[ y_k = x_{k+1} - x_k \\\\]. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::setStepVector "void BroydenOperator::setStepVector(NOX::Epetra::Vector &vec)

Set the current step vector, \\\\[ y_k = x_{k+1} - x_k \\\\]. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::setYieldVector "void BroydenOperator::setYieldVector(Epetra_Vector &vec)

Set the current yield vector, \\\\[ y_k = F_{k+1} - F_k \\\\]. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::setYieldVector "void BroydenOperator::setYieldVector(NOX::Epetra::Vector &vec)

Set the current yield vector, \\\\[ y_k = F_{k+1} - F_k \\\\]. ";

%feature("docstring")
NOX::Epetra::BroydenOperator::computeSparseBroydenUpdate "bool
BroydenOperator::computeSparseBroydenUpdate()

Compute the sparse Broyden update. ";

%feature("docstring")
NOX::Epetra::BroydenOperator::removeEntriesFromBroydenUpdate "void
BroydenOperator::removeEntriesFromBroydenUpdate(const Epetra_CrsGraph
&graph)

Remove entries from being involved in Broyden updates. ";

%feature("docstring")  NOX::Epetra::BroydenOperator::getBroydenMatrix
"const Epetra_CrsMatrix&
NOX::Epetra::BroydenOperator::getBroydenMatrix()

Return a const reference to the Broyden matrix. The matrix is not
owned but is obtained from the client at construction. ";

%feature("docstring")
NOX::Epetra::BroydenOperator::resetBroydenMatrix "void
BroydenOperator::resetBroydenMatrix(const Epetra_CrsMatrix &mat)

Reset the values of our matrix. ";

%feature("docstring")
NOX::Epetra::BroydenOperator::addReplacementInterface "void
NOX::Epetra::BroydenOperator::addReplacementInterface(ReplacementInterface
*i)

Register replacement interface. ";


// File: classNOX_1_1Epetra_1_1BroydenOperator_1_1ReplacementInterface.xml
%feature("docstring")
NOX::Epetra::BroydenOperator::ReplacementInterface "";

%feature("docstring")
NOX::Epetra::BroydenOperator::ReplacementInterface::~ReplacementInterface
"virtual
NOX::Epetra::BroydenOperator::ReplacementInterface::~ReplacementInterface()
";

%feature("docstring")
NOX::Epetra::BroydenOperator::ReplacementInterface::getReplacementValuesMatrix
"virtual Teuchos::RCP<const Epetra_CrsMatrix>
NOX::Epetra::BroydenOperator::ReplacementInterface::getReplacementValuesMatrix(const
Epetra_Vector &x, FILL_TYPE)=0 ";


// File: classNOX_1_1StatusTest_1_1Combo.xml
%feature("docstring") NOX::StatusTest::Combo "

Arbitrary combination of status tests.

In the AND (see NOX::StatusTest::Combo::ComboType) combination, the
result is Unconverged (see NOX::StatusTest::StatusType) if any of the
tests is Unconverged. Otherwise, the result is equal to the result of
the first test in the list that is either Converged or Failed. It is
not recommended to mix Converged and Failed tests in an AND
combination.

In the OR combination, the result is Unconverged if all of the tests
are Unconverged. Otherwise, it is the result of the first test in the
list that is either Converged or Failed. Therefore, it will generally
make sense to put the Failed -type tests at the end of the OR list.

We call checkStatus on every convergence test, though some may be
called with the NOX::StatusTest::None option.

Tammy Kolda (SNL 8950) and Roger Pawlowski (SNL 1416)

C++ includes: NOX_StatusTest_Combo.H ";

%feature("docstring")  NOX::StatusTest::Combo::Combo "NOX::StatusTest::Combo::Combo(ComboType t, const NOX::Utils *u=NULL)

Constructor. Optional argument is the error stream for output. ";

%feature("docstring")  NOX::StatusTest::Combo::Combo "NOX::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP< Generic
> &a, const NOX::Utils *u=NULL)

Constructor with a single test. ";

%feature("docstring")  NOX::StatusTest::Combo::Combo "NOX::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP< Generic
> &a, const Teuchos::RCP< Generic > &b, const NOX::Utils *u=NULL)

Constructor with two tests. ";

%feature("docstring")  NOX::StatusTest::Combo::addStatusTest "NOX::StatusTest::Combo & NOX::StatusTest::Combo::addStatusTest(const
Teuchos::RCP< Generic > &a)

Add another test to this combination.

Calls isSafe() to determine if it is safe to add a to the combination.
";

%feature("docstring")  NOX::StatusTest::Combo::~Combo "NOX::StatusTest::Combo::~Combo()

Destructor. ";

%feature("docstring")  NOX::StatusTest::Combo::checkStatus "NOX::StatusTest::StatusType NOX::StatusTest::Combo::checkStatus(const
NOX::Solver::Generic &problem, NOX::StatusTest::CheckType checkType)

Tests stopping criterion.

See addOp() and orOp() for details. ";

%feature("docstring")  NOX::StatusTest::Combo::getStatus "NOX::StatusTest::StatusType NOX::StatusTest::Combo::getStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::Combo::print "ostream &
NOX::StatusTest::Combo::print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classNOX_1_1LineSearch_1_1Utils_1_1Counters.xml
%feature("docstring") NOX::LineSearch::Utils::Counters "

Common counters that all line search algorithms should report.

Output list keys

Line searches have a common set of data that should be tracked and
saved in the parameter list for the users. This class provides a set
of common data objects, accessors, and a routine to print them. A
sublist for output parameters will be created called \"Output\" in the
\"Line Search\" parameter sublist. Valid output keys for the parameter
list are:

\"Total Number of Line Search Calls\" - Total number of calls to the
compute() method of this line search.

\"Total Number of Non-trivial Line Searches\" - The total number of
steps that could not directly take a full step and meet the required
\"Convergence Criteria\" (i.e. The line search had to reduce the step
length using inner iteration calculations over iterate $ k $).

\"Total Number of Failed Line Searches\" - total number of line
searches that failed and used a recovery step.

\"Total Number of Line Search Inner Iterations\" - total number of
inner iterations $ k $ performed by this object.

C++ includes: NOX_LineSearch_Utils_Counters.H ";

/*  Common Line Search Data Members  */

/* All line searches in general should report the following
information. We save a pointer to the parameter list associated with
the line search and set an output sublist with each of the following
parameters.

*/

/*  Increment Methods  */

%feature("docstring")
NOX::LineSearch::Utils::Counters::incrementNumLineSearches "void
NOX::LineSearch::Utils::Counters::incrementNumLineSearches(int n=1)

Increment the counter for the total number of line search calls by n.
";

%feature("docstring")
NOX::LineSearch::Utils::Counters::incrementNumNonTrivialLineSearches "void
NOX::LineSearch::Utils::Counters::incrementNumNonTrivialLineSearches(int
n=1)

Increment the counter for the total number of non-trivial line search
calls by n. ";

%feature("docstring")
NOX::LineSearch::Utils::Counters::incrementNumFailedLineSearches "void
NOX::LineSearch::Utils::Counters::incrementNumFailedLineSearches(int
n=1)

Increment the counter for the total number of failed line search calls
by n. ";

%feature("docstring")
NOX::LineSearch::Utils::Counters::incrementNumIterations "void
NOX::LineSearch::Utils::Counters::incrementNumIterations(int n=1)

Increment the counter for the total number of line search inner
iterations by n. ";

/*  Accessor Methods  */

/* Returns the current counter value

*/

%feature("docstring")
NOX::LineSearch::Utils::Counters::getNumLineSearches "int
NOX::LineSearch::Utils::Counters::getNumLineSearches() const

Return the counter for the total number of line search calls. ";

%feature("docstring")
NOX::LineSearch::Utils::Counters::getNumNonTrivialLineSearches "int
NOX::LineSearch::Utils::Counters::getNumNonTrivialLineSearches() const

Return the counter for the total number of non-trivial line search
calls. ";

%feature("docstring")
NOX::LineSearch::Utils::Counters::getNumFailedLineSearches "int
NOX::LineSearch::Utils::Counters::getNumFailedLineSearches() const

Return the counter for the total number of failed line search calls.
";

%feature("docstring")
NOX::LineSearch::Utils::Counters::getNumIterations "int
NOX::LineSearch::Utils::Counters::getNumIterations() const

Return the counter for the total number of line search inner
iterations. ";

%feature("docstring")  NOX::LineSearch::Utils::Counters::Counters "NOX::LineSearch::Utils::Counters::Counters()

Default constructor. ";

%feature("docstring")  NOX::LineSearch::Utils::Counters::~Counters "NOX::LineSearch::Utils::Counters::~Counters()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Utils::Counters::reset "void
NOX::LineSearch::Utils::Counters::reset()

Reset the counters . ";

%feature("docstring")  NOX::LineSearch::Utils::Counters::setValues "bool
NOX::LineSearch::Utils::Counters::setValues(Teuchos::ParameterList
&lineSearchParams)

Sets the common line search data in an \"Output\" sublist of the
lineSearchParams list that is supplied in the method call. ";


// File: classNOX_1_1StatusTest_1_1Divergence.xml
%feature("docstring") NOX::StatusTest::Divergence "

Failure test based on a threshold value of the norm of F.

This status test returns NOX::StatusTest::Failed if the norm of F
exceeds a threshold value for n consecutive iterations, i.e.

\\\\[ \\\\| F \\\\| > {\\\\rm threshold}\\\\]

for n consecutive iterations, the status is set to
NOX::StatusTest::Failed and returned. Otherwise, the status is set to
NOX::StatusTest::Uncoverged and returned. Both n and the threshold are
specified in the constructor, by n and thresh, respectively. While a
value for thresh must be provided, a default value of n = 1 is
assumed.

C++ includes: NOX_StatusTest_Divergence.H ";

%feature("docstring")  NOX::StatusTest::Divergence::Divergence "NOX::StatusTest::Divergence::Divergence(double thresh, int n=1)

Constructor.

Parameters:
-----------

thresh:  - Threshold for divergence test

n:  - Number of consecutive nonlinear iterations ";

%feature("docstring")  NOX::StatusTest::Divergence::~Divergence "NOX::StatusTest::Divergence::~Divergence()

Destructor. ";

%feature("docstring")  NOX::StatusTest::Divergence::checkStatus "NOX::StatusTest::StatusType
NOX::StatusTest::Divergence::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)

Tests the stopping criterion.

The nature of this test is such that it must be executed at every
nonlinear iteration, so we don't use the checkType argument. ";

%feature("docstring")  NOX::StatusTest::Divergence::getStatus "NOX::StatusTest::StatusType NOX::StatusTest::Divergence::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::Divergence::print "ostream &
NOX::StatusTest::Divergence::print(ostream &stream, int indent=0)
const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  NOX::StatusTest::Divergence::getMaxNumSteps "int NOX::StatusTest::Divergence::getMaxNumSteps() const

Returns the user-specified number of steps that can consecutively fail
the threshold test before the test returns a failed status. ";

%feature("docstring")  NOX::StatusTest::Divergence::getCurrentNumSteps
"int NOX::StatusTest::Divergence::getCurrentNumSteps() const

Returns the current number of steps that have consecutively failed the
threshold test. ";

%feature("docstring")  NOX::StatusTest::Divergence::getThreshold "double NOX::StatusTest::Divergence::getThreshold() const

Returns the user-specified threshold. ";


// File: classNOX_1_1Direction_1_1Factory.xml
%feature("docstring") NOX::Direction::Factory "

Factory to build direction objects derived from
NOX::Direction::Generic.

Parameters

\"Method\" <std::string> Name of the direction. Valid choices are:
\"Newton\" NOX::Direction::Newton

\"Steepest Descent\" NOX::Direction::SteepestDescent

\"NonlinearCG\" NOX::Direction::NonlinearCG

\"Broyden\" NOX::Direction::Broyden

\"Tensor\" PRERELEASE ONLY! NOX::Direction::Tensor

\"Modified-Newton\" PRERELEASE ONLY! NOX::Direction::ModifiedNewton

\"Quasi-Newton\" PRERELEASE ONLY! NOX::Direction::QuasiNewton

\"User Defined\" - see below

\"User Defined Constructor\" - see below

\"Newton\" <sublist> Parameters to build a NOX::Direction::Newton
object.

\"Steepest Descent\" <sublist> Parameters to build a
NOX::Direction::SteepestDescent object.

\"NonlinearCG\" <sublist> Parameters to build a
NOX::Direction::NonlinearCG object.

\"Broyden\" <sublist> Parameters to build a NOX::Direction::Broyden
object.

\"Tensor\" <sublist> Parameters to build a NOX::Direction::Tensor
object.

\"Modified-Newton\" <sublist> Parameters to build a
NOX::Direction::ModifiedNewton object.

\"Quasi-Newton\" <sublist> Parameters to build a
NOX::Direction::QuasiNewton object.

\"User Defined Direction Factory\" <
RCP<NOX::Direction::UserDefinedFactory> > RCP to a
NOX::Direction::UserDefinedFactory derived object. This factory object
is used to build user defined direction objects.

Using a User-Defined Direction

The user has the option of passing in a user-defined direction. First,
they must implement their own direction, deriving from the base class
interface NOX::Direction::Generic:

Next they must write a factory to build their object, deriving from
the NOX::Direction::UserDefinedFactory base class interface:

Then under the \"Direction\" parameter sublist, they need to set the
method to \"User Defined\" and register the factory:

It is critical that the user defined factory be set in the parameter
list as a base class type object: NOX::Direction::UserDefinedFactory.

C++ includes: NOX_Direction_Factory.H ";

%feature("docstring")  NOX::Direction::Factory::Factory "NOX::Direction::Factory::Factory()

Constructor. ";

%feature("docstring")  NOX::Direction::Factory::~Factory "NOX::Direction::Factory::~Factory()

Destructor. ";

%feature("docstring")  NOX::Direction::Factory::buildDirection "Teuchos::RCP< NOX::Direction::Generic > buildDirection(const
Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Factory to build a direction object.

Parameters:
-----------

gd:  A global data pointer that contains the top level parameter list.
Without storing this inside the direction object, there is no
guarantee that the second parameter params will still exist. It can be
deleted by the top level RCP.

params:  Sublist with direction construction parameters.

Nonmember function to build a direction object. ";


// File: classNOX_1_1LineSearch_1_1Factory.xml
%feature("docstring") NOX::LineSearch::Factory "

Factory to build line search objects derived from
NOX::LineSearch::Generic.

Parameters

\"Method\" <std::string> Name of the line search. Valid choices are:
\"Full Step\" ( NOX::LineSearch::FullStep)

\"Backtrack\" ( NOX::LineSearch::Backtrack)

\"Polynomial\" ( NOX::LineSearch::Polynomial)

\"More'-Thuente\" ( NOX::LineSearch::MoreThuente)

\"User Defined\" - see below

\"User Defined Constructor\" - see below

\"Full Step\" <sublist> Parameters to build a
NOX::LineSearch::FullStep sublist.

\"Backtrack\" <sublist> Parameters to build a
NOX::LineSearch::Backtrack sublist.

\"Polynomial\" <sublist> Parameters to build a
NOX::LineSearch::Polynomial sublist.

\"More'-Thuente\" <sublist> Parameters to build a
NOX::LineSearch::MoreThuente sublist.

\"User Defined Line Search Factory\" <
RCP<NOX::LineSearch::UserDefinedFactory> > RCP to a
NOX::LineSearch::UserDefinedFactory derived object. This factory
object is used to build user defined line search objects.

Using a User-Defined Line Search

The user has the option of passing in a user-defined line search.
First, they must implement their own line search, deriving from the
base class interface NOX::LineSearch::Generic:

Next they must write a factory to build their object, deriving from
the NOX::LineSearch::UserDefinedFactory base class interface:

Then under the \"Line Search\" parameter sublist, they need to set the
method to \"User Defined\" and register the factory:

It is critical that the user defined factory be set in the parameter
list as a base class type object: NOX::LineSearch::UserDefinedFactory.

C++ includes: NOX_LineSearch_Factory.H ";

%feature("docstring")  NOX::LineSearch::Factory::Factory "NOX::LineSearch::Factory::Factory()

Constructor. ";

%feature("docstring")  NOX::LineSearch::Factory::~Factory "NOX::LineSearch::Factory::~Factory()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Factory::buildLineSearch "Teuchos::RCP< NOX::LineSearch::Generic > buildLineSearch(const
Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Factory to build a line search object.

Parameters:
-----------

gd:  A global data pointer that contains the top level parameter list.
Without storing this inside the line searchobject, there is no
guarantee that the second parameter params will still exist. It can be
deleted by the top level RCP.

params:  Sublist with line search construction parameters.

Nonmember function to build a line search object. ";


// File: classNOX_1_1Solver_1_1Factory.xml
%feature("docstring") NOX::Solver::Factory "

Factory class to control the creation of solvers derived from the
NOX::Solver::Generic object.

Parameters

The following entries may be specified in the parameter list.

\"Nonlinear Solver\" <std::string> Name of the solver method. Valid
choices are \"Line Search Based\" or \"Newton\" (
NOX::Solver::LineSearchBased) [Default]

\"Trust Region Based\" ( NOX::Solver::TrustRegionBased)

\"Inexact Trust Region Based\" ( NOX::Solver::InexactTrustRegionBased)

\"Tensor Based\" ( NOX::Solver::TensorBased)

\"Line Search Based\" <Teuchos::ParameterList> Sublist for the
NOX::Solver::LineSearchBased solver.

\"Trust Region Based\" <Teuchos::ParameterList> Sublist for the
NOX::Solver::TrustRegionBased solver.

\"Inexact Trust Region Based\" <Teuchos::ParameterList> Sublist for
the NOX::Solver::InexactTrustRegionBased solver.

\"Tensor Based\" <Teuchos::ParameterList> Sublist for the
NOX::Solver::TensorBased solver.

\"Tensor-Krylov Based\" <Teuchos::ParameterList> Sublist for the
NOX::Solver::TensorBasedTest solver. (Prerelease only)

Solvers can be constructed using a nonmember function
NOX::Solver::buildSolver instead of using this object directly.

Roger Pawlowski (SNL 1416)

C++ includes: NOX_Solver_Factory.H ";

%feature("docstring")  NOX::Solver::Factory::Factory "NOX::Solver::Factory::Factory()

Constructor. ";

%feature("docstring")  NOX::Solver::Factory::~Factory "NOX::Solver::Factory::~Factory()

Destructor. ";

%feature("docstring")  NOX::Solver::Factory::buildSolver "Teuchos::RCP< NOX::Solver::Generic > buildSolver(const Teuchos::RCP<
NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)

Nonmember helper function for the NOX::Solver::Factory. ";


// File: classNOX_1_1StatusTest_1_1Factory.xml
%feature("docstring") NOX::StatusTest::Factory "

Factory to build a set of status tests from a parameter list.

This object takes either an xml file name or a Teuchos::ParameterList
and generates an entire set (a tree) of status tests for use in a
NOX::Solver derived object.

The tagged_tests field in the constructors allows users to store tests
from the tree in a flat list in case they want to change the tolerance
values during a run. The tagged_tests flag is optional.

Please use the related nonmember functions instead of calling the
factory directly (See example below).

Valid parameters are as follows:

\"Test Type\" <std::string> Type of test this list contains. Valid
tests include: \"Combo\" - NOX::StatusTest::Combo

\"NormF\" - NOX::StatusTest::NormF

\"NormUpdate\" - NOX::StatusTest::NormUpdate

\"NomrWRMS\" - NOX::StatusTest::NormWRMS

\"MaxIters\" - NOX::StatusTest::MaxIters

\"FiniteValue\" - NOX::StatusTest::FiniteValue

\"Divergence\" - NOX::StatusTest::Divergence

\"Stagnation\" - NOX::StatusTest::Stagnation

\"User Defined\" - A user constructed test, derived from
NOX::StatusTest::Generic.

\"Tag\" <std::string> A unique identifier that will place the test in
the map for tagged_tests. This allows users to access individual tests
to change tolerances on the fly or query values while still using the
factory to build objects.

Additional parameters valid for a Combo test (
NOX::StatusTest::Combo): \"Combo Type\" <std:string> Type of combo to
use. Valid options are: \"AND\"

\"OR\"

\"Number of Tests\" <int> Number of sublists that contain tests to be
added to this combo test. The sublists must be named \"Test X\" where
\"X\" represents the test number starting with 0 and preceeding to
\"Number of Tests - 1\".

\"Test X\" <Teuchos::ParameterList> A sublist containing a test to add
to the current combo test. The \"X\" represents the number of the
test. the numbering starts with 0 and is valid through \"Number of
Tests - 1\" tests.

Additional parameters valid for a Norm F test (
NOX::StatusTest::NormF): \"Tolerance\" <double> required tolerance for
the test to return a converged status. (default = 1.0e-8)

\"Norm Type\" <std::string> Type of norm to use. Valid options are:
\"Two Norm\" (default)

\"One Norm\"

\"Max Norm\"

\"Scale Type\" <std::string> Type of scaling to use. Valid options
are: \"Unscaled\" (default)

\"Scaled\"

\"Initial Guess\" < Teuchos::RCP<NOX::Abstract::Group> > If present, a
relative tolerance will be used where the group passed in will be used
to compute $F_0$.

Additional parameters valid for a Norm Update test (
NOX::StatusTest::NormUpdate): \"Tolerance\" <double> required
tolerance for the test to return a converged status. (default =
1.0e-3)

\"Norm Type\" <std::string> Type of norm to use. Valid options are:
\"Two Norm\" (default)

\"One Norm\"

\"Max Norm\"

\"Scale Type\" <std::string> Type of scaling to use. Valid options
are: \"Unscaled\" (default)

\"Scaled\"

Additional parameters valid for a Norm WRMS test (
NOX::StatusTest::NormWRMS): \"Tolerance\" <double> (default = 1.0)

\"Relative Tolerance\" <double> (default = 1.0e-5)

\"Absolute Tolerance\" <double> or < Teuchos::RCP<const
NOX::Abstract::Vector> > (default = 1.0e-8)

\"BDF Multiplier\" <double> (default = 1.0)

\"Alpha\" <double> (default = 1.0)

\"Beta\" <double> (default = 0.5)

Additional parameters valid for a Maximum Iterations test (
NOX::StatusTest::MaxIters): \"Maximum Iterations\" <int>

Additional parameters valid for a Finite Value test (
NOX::StatusTest::FiniteValue): \"Vector Type\" <std::string> Type of
vector to check. Valid options are: \"F Vector\" (default)

\"Solution Vector\"

\"Norm Type\" <std::string> Type of norm to use. Valid options are:
\"Two Norm\" (default)

\"One Norm\"

\"Max Norm\"

Additional parameters valid for a Divergence test (
NOX::StatusTest::Divergence): \"Tolerance\" <double> (default =
1.0e+12)

\"Consecutive Iterations\" <int> (default = 1)

Additional parameters valid for a Stagnation test (
NOX::StatusTest::Stagnation): \"Tolerance\" <double> (default =
1.0e+12)

\"Consecutive Iterations\" <int> (default = 1)

Additional parameters valid for a \"User Defined\" test: \"User Status
Test\" < Teuchos::RCP<NOX::StatusTest::Generic> > A status test
suppied by the user. It is very important that when registering this
status test, that the user set it as a \"Generic\" object since there
is no implicit casting on the ParameterList's get method. See the
example below.

Example usage:

Roger Pawlowski (SNL 1416)

C++ includes: NOX_StatusTest_Factory.H ";

%feature("docstring")  NOX::StatusTest::Factory::Factory "NOX::StatusTest::Factory::Factory()

Constructor. ";

%feature("docstring")  NOX::StatusTest::Factory::~Factory "NOX::StatusTest::Factory::~Factory()

Destructor. ";

%feature("docstring")  NOX::StatusTest::Factory::buildStatusTests "Teuchos::RCP< NOX::StatusTest::Generic >
NOX::StatusTest::Factory::buildStatusTests(const std::string
&file_name, const NOX::Utils &utils, std::map< std::string,
Teuchos::RCP< NOX::StatusTest::Generic > > *tagged_tests=0) const

Returns a status test set from a parameter list xml file. ";

%feature("docstring")  NOX::StatusTest::Factory::buildStatusTests "Teuchos::RCP< NOX::StatusTest::Generic >
NOX::StatusTest::Factory::buildStatusTests(Teuchos::ParameterList &p,
const NOX::Utils &utils, std::map< std::string, Teuchos::RCP<
NOX::StatusTest::Generic > > *tagged_tests=0) const

Returns a status test set from a parameter list. ";


// File: classNOX_1_1Epetra_1_1FiniteDifference.xml
%feature("docstring") NOX::Epetra::FiniteDifference "

Concrete implementation for creating an Epetra_RowMatrix Jacobian via
finite differencing of the residual.

The Jacobian entries are calculated via 1st order finite differencing.
This requires $ N + 1 $ calls to computeF() where $ N $ is the number
of unknowns in the problem.

\\\\[ J_{ij} = \\\\frac{\\\\partial F_i}{\\\\partial x_j} =
\\\\frac{F_i(x+\\\\delta\\\\mathbf{e}_j) - F_i(x)}{\\\\delta} \\\\]

where $J$ is the Jacobian, $F$ is the function evaluation, $x$ is the
solution vector, and $\\\\delta$ is a small perturbation to the $x_j$
entry.

The perturbation, $ \\\\delta $, is calculated based on one of the
following equations:

\\\\[ \\\\delta = \\\\alpha * | x_j | + \\\\beta \\\\] \\\\[ \\\\delta
= \\\\alpha * | x_j | + \\\\beta_j \\\\]

where $ \\\\alpha $ is a scalar value (defaults to 1.0e-4) and $
\\\\beta $ can be either a scalar or a vector (defaults to a scalar
value of 1.0e-6). The choice is defined by the type of constructor
used. All parameters are supplied in the constructor. In addition to
the forward difference derivative approximation, backward or centered
differences can be used via the setDifferenceMethod function. Note
that centered difference provides second order spatial accuracy but at
the cost of twice as many function evaluations.

Since this inherits from the Epetra_RowMatrix class, it can be used as
the preconditioning matrix for AztecOO preconditioners. This method is
very inefficient when computing the Jacobian and is not recommended
for large-scale systems but only for debugging purposes.

C++ includes: NOX_Epetra_FiniteDifference.H ";

%feature("docstring")  NOX::Epetra::FiniteDifference::FiniteDifference
"FiniteDifference::FiniteDifference(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< NOX::Epetra::Interface::Required
> &i, const NOX::Epetra::Vector &initialGuess, double beta=1.0e-6,
double alpha=1.0e-4)

Constructor with scalar beta. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::FiniteDifference
"FiniteDifference::FiniteDifference(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< NOX::Epetra::Interface::Required
> &i, const NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
const Epetra_Vector > &beta, double alpha=1.0e-4)

Constructor with vector beta. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::FiniteDifference
"FiniteDifference::FiniteDifference(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< NOX::Epetra::Interface::Required
> &i, const NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
Epetra_CrsGraph > &g, double beta=1.0e-6, double alpha=1.0e-4)

Constructor that takes a pre-constructed Epetra_CrsGraph so it does
not have to determine the non-zero entries in the matrix. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::FiniteDifference
"FiniteDifference::FiniteDifference(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< NOX::Epetra::Interface::Required
> &i, const NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
Epetra_CrsGraph > &g, const Teuchos::RCP< const Epetra_Vector > &beta,
double alpha=1.0e-4)

Constructor with output control that takes a pre-constructed
Epetra_CrsGraph so it does not have to determine the non-zero entries
in the matrix. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::~FiniteDifference "FiniteDifference::~FiniteDifference()

Pure virtual destructor. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Label "const
char * FiniteDifference::Label() const

Returns a character string describing the name of the operator. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::SetUseTranspose
"int FiniteDifference::SetUseTranspose(bool UseTranspose)

If set true, the transpose of this operator will be applied. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Apply "int
FiniteDifference::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Return the result on an Epetra_Operator applied to an
Epetra_MultiVector X in Y. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::ApplyInverse "int FiniteDifference::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Return the result on an Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::UseTranspose "bool FiniteDifference::UseTranspose() const

Returns the current use transpose setting. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::HasNormInf "bool FiniteDifference::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::OperatorDomainMap "const Epetra_Map &
FiniteDifference::OperatorDomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::OperatorRangeMap
"const Epetra_Map & FiniteDifference::OperatorRangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Filled "bool
FiniteDifference::Filled() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumMyRowEntries
"int FiniteDifference::NumMyRowEntries(int MyRow, int &NumEntries)
const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::MaxNumEntries "int FiniteDifference::MaxNumEntries() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::ExtractMyRowCopy
"int FiniteDifference::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::ExtractDiagonalCopy "int
FiniteDifference::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Multiply "int
FiniteDifference::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Solve "int
FiniteDifference::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::InvRowSums "int
FiniteDifference::InvRowSums(Epetra_Vector &x) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::LeftScale "int
FiniteDifference::LeftScale(const Epetra_Vector &x)

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::InvColSums "int
FiniteDifference::InvColSums(Epetra_Vector &x) const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::RightScale "int
FiniteDifference::RightScale(const Epetra_Vector &x)

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NormInf "double
FiniteDifference::NormInf() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NormOne "double
FiniteDifference::NormOne() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::NumGlobalNonzeros "int
FiniteDifference::NumGlobalNonzeros() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumGlobalRows "int FiniteDifference::NumGlobalRows() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumGlobalCols "int FiniteDifference::NumGlobalCols() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::NumGlobalDiagonals "int
FiniteDifference::NumGlobalDiagonals() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumMyNonzeros "int FiniteDifference::NumMyNonzeros() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumMyRows "int
FiniteDifference::NumMyRows() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumMyCols "int
FiniteDifference::NumMyCols() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::NumMyDiagonals "int FiniteDifference::NumMyDiagonals() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::LowerTriangular
"bool FiniteDifference::LowerTriangular() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::UpperTriangular
"bool FiniteDifference::UpperTriangular() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Comm "const
Epetra_Comm & FiniteDifference::Comm() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::RowMatrixRowMap
"const Epetra_Map & FiniteDifference::RowMatrixRowMap() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::RowMatrixColMap
"const Epetra_Map & FiniteDifference::RowMatrixColMap() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::RowMatrixImporter "const Epetra_Import
* FiniteDifference::RowMatrixImporter() const

See Epetra_RowMatrix documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Map "const
Epetra_BlockMap & FiniteDifference::Map() const

See Epetra_SrcDistObj documentation. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::computeJacobian
"bool FiniteDifference::computeJacobian(const Epetra_Vector &x,
Epetra_Operator &Jac)

Compute Jacobian given the specified input vector, x. Returns true if
computation was successful. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::computeJacobian
"bool FiniteDifference::computeJacobian(const Epetra_Vector &x)

Compute Jacobian given the specified input vector, x. Returns true if
computation was successful. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::computePreconditioner "bool
FiniteDifference::computePreconditioner(const Epetra_Vector &x,
Epetra_Operator &Prec, Teuchos::ParameterList *precParams=0)

Compute an Epetra_RowMatrix to be used by Aztec preconditioners given
the specified input vector, x. Returns true if computation was
successful. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::setDifferenceMethod "void
FiniteDifference::setDifferenceMethod(DifferenceType type)

Set the type of perturbation method used (default is Forward). ";

%feature("docstring")
NOX::Epetra::FiniteDifference::getUnderlyingMatrix "Epetra_CrsMatrix
& FiniteDifference::getUnderlyingMatrix() const

An accessor method for the underlying Epetra_CrsMatrix. ";

%feature("docstring")  NOX::Epetra::FiniteDifference::Print "void
FiniteDifference::Print(ostream &) const

Output the underlying matrix. ";

%feature("docstring")
NOX::Epetra::FiniteDifference::setGroupForComputeF "void
FiniteDifference::setGroupForComputeF(NOX::Abstract::Group &group)

Register a NOX::Abstract::Group derived object and use the computeF()
method of that group for the perturbation instead of the
NOX::Epetra::Interface::Required::computeF() method. This is required
for LOCA to get the operators correct during homotopy. ";


// File: classNOX_1_1Epetra_1_1FiniteDifferenceColoring.xml
%feature("docstring") NOX::Epetra::FiniteDifferenceColoring "

Concrete implementation for creating an Epetra_RowMatrix Jacobian via
finite differencing of the residual using coloring.

The Jacobian entries are calculated via 1st or 2nd order finite
differencing. This requires $ N + 1 $ or $ 2N + 1 $ calls to
computeF(), respectively, where $ N $ is the number of colors.

\\\\[ J_{ij} = \\\\frac{\\\\partial F_i}{\\\\partial x_j} =
\\\\frac{F_i(x+\\\\delta\\\\mathbf{e}_j) - F_i(x)}{\\\\delta} \\\\]

where $J$ is the Jacobian, $F$ is the function evaluation, $x$ is the
solution vector, and $\\\\delta$ is a small perturbation to the $x_j$
entry.

Instead of perturbing each $ N_{dof} $ problem degrees of freedom
sequentially and then evaluating all $ N_{dof} $ functions for each
perturbation, coloring allows several degrees of freedom (all
belonging to the same color) to be perturbed at the same time. This
reduces the total number of function evaluations needed to compute
$\\\\mathbf{J}$ from $ N_{dof}^2 $ as is required using
FiniteDifference to $ N\\\\cdot N_{dof} $, often representing
substantial computational savings.

Coloring is based on a user-supplied color map generated using an
appropriate algorithm, eg greedy-algorithm - Y. Saad, \"Iterative
Methods for Sparse Linear Systems, 2nd ed.,\" chp. 3, SIAM, 2003.. Use
can be made of the coloring algorithm provided by the EpetraExt
package in Trilinos. The 1Dfem_nonlinearColoring and Brusselator
example problems located in the nox/epetra-examples subdirectory
demonstrate use of the EpetraExt package, and the
1Dfem_nonlinearColoring directory also contains a stand-alone coloring
algorithm very similar to that in EpetraExt.

The perturbation, $ \\\\delta $, is calculated using the following
equation:

\\\\[ \\\\delta = \\\\alpha * | x_j | + \\\\beta \\\\]

where $ \\\\alpha $ is a scalar value (defaults to 1.0e-4) and $
\\\\beta $ is another scalar (defaults to 1.0e-6).

Since both FiniteDifferenceColoring and FiniteDifference inherit from
the Epetra_RowMatrix class, they can be used as preconditioning
matrices for AztecOO preconditioners.

As for FiniteDifference, 1st order accurate Forward and Backward
differences as well as 2nd order accurate Centered difference can be
specified using setDifferenceMethod with the appropriate enumerated
type passed as the argument.

Using FiniteDifferenceColoring in Parallel

Two ways of using this class in a distributed parallel environment are
currently supported. From an application standpoint, the two
approaches differ only in the status of the solution iterate used in
the residual fill. If an object of this class is contructed with
parallelColoring = true the solution iterate will be passe back in a
non-ghosted form. On the contrary, setting this parameter to false in
the constructor will cause the solution iterate to be in a ghosted
form when calling back for a residual fill. When using the second
approach, the user should be aware that the perturbed vector used to
compute residuals has already been scattered to a form consistent with
the column space of the Epetra_CrsGraph. In practice, this means that
the perturbed vector used by computeF() has already been scattered to
a ghosted or overlapped state. The application should then not perform
this step but rather simply use the vector provided with the possible
exception of requiring a local index reordering to bring the column-
space based vector in sync with a potentially different ghosted index
ordering. See the Brusselator and 1Dfem_nonlinearColoring example
problems for details.

Special Case for Approximate Jacobian Construction

Provision is made for a simplified and cheaper use of coloring that
currently provides only for the diagonal of the Jacobian to be
computed. This is based on using a first-neighbors coloring of the
original Jacobian graph using the Epetra_Ext MapColoring class with
the distance1 argument set to true. This same argument should also be
set to true in the constructor to this class. The result will be a
diagonal Jacobian filled in a much more efficient manner.

C++ includes: NOX_Epetra_FiniteDifferenceColoring.H ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoring::FiniteDifferenceColoring "FiniteDifferenceColoring::FiniteDifferenceColoring(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
Epetra_MapColoring > &colorMap, const Teuchos::RCP< vector<
Epetra_IntVector > > &columns, bool parallelColoring=false, bool
distance1=false, double beta=1.0e-6, double alpha=1.0e-4)

Constructor with output control. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoring::FiniteDifferenceColoring "FiniteDifferenceColoring::FiniteDifferenceColoring(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP< Epetra_CrsGraph
> &rawGraph, const Teuchos::RCP< Epetra_MapColoring > &colorMap, const
Teuchos::RCP< vector< Epetra_IntVector > > &columns, bool
parallelColoring=false, bool distance1=false, double beta=1.0e-6,
double alpha=1.0e-4)

Constructor with output control. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoring::~FiniteDifferenceColoring "FiniteDifferenceColoring::~FiniteDifferenceColoring()

Pure virtual destructor. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoring::computeJacobian "bool
FiniteDifferenceColoring::computeJacobian(const Epetra_Vector &x,
Epetra_Operator &Jac)

Compute Jacobian given the specified input vector, x. Returns true if
computation was successful. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoring::computeJacobian "bool
FiniteDifferenceColoring::computeJacobian(const Epetra_Vector &x)

Compute Jacobian given the specified input vector, x. Returns true if
computation was successful. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoring::createColorContainers "void
FiniteDifferenceColoring::createColorContainers()

Output the coloring map, index map and underlying matrix.

Create containers for using color and index maps in parallel coloring
";


// File: classNOX_1_1Epetra_1_1FiniteDifferenceColoringWithUpdate.xml
%feature("docstring") NOX::Epetra::FiniteDifferenceColoringWithUpdate
"

Concrete implementation for creating an Epetra_RowMatrix Jacobian via
finite differencing of the residual using coloring. This method
assumes the existence of a valid parallel coloring of the columns of
the Jacobian (aka from Isorropia).

Unlike the class NOX::FiniteDifferenceColoring, this class allows for
\"update\" colorings, for use in situations where part of the Jacobian
changes from iteration to iteration, but part does not. The first time
(or any time after the forceJacobianRecompute method is called) the
method uses the complete coloring. Afterwards, it uses the \"update\"
coloring and only changes the entries that can change.

WARNING: The \"update\" coloring assumes that rows AND columns
corresponding to uncolored (aka color 0) nodes do not change from call
to call. If either the row or the column corresponding to a given node
change then you must make sure it gets colored.

WARNING: Centered Difference Coloring is NOT supported as of yet.

C++ includes: NOX_Epetra_FiniteDifferenceColoringWithUpdate.H ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate
"FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
Epetra_MapColoring > &colorMap, double beta=1.0e-6, double
alpha=1.0e-4)

Constructor with no frills. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate
"FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP< Epetra_CrsGraph
> &rawGraph, const Teuchos::RCP< Epetra_MapColoring > &colorMap,
double beta=1.0e-6, double alpha=1.0e-4)

Constructor with graph. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate
"FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
Epetra_MapColoring > &colorMap, const Teuchos::RCP< Epetra_MapColoring
> &updateColorMap, double beta=1.0e-6, double alpha=1.0e-4)

Constructor with update map. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate
"FiniteDifferenceColoringWithUpdate::FiniteDifferenceColoringWithUpdate(Teuchos::ParameterList
&printingParams, const Teuchos::RCP< Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP< Epetra_CrsGraph
> &rawGraph, const Teuchos::RCP< Epetra_MapColoring > &colorMap, const
Teuchos::RCP< Epetra_MapColoring > &updatecolorMap, double
beta=1.0e-6, double alpha=1.0e-4)

Constructor with update map and graph. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::forceJacobianRecompute
"virtual void
NOX::Epetra::FiniteDifferenceColoringWithUpdate::forceJacobianRecompute()
";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::~FiniteDifferenceColoringWithUpdate
"FiniteDifferenceColoringWithUpdate::~FiniteDifferenceColoringWithUpdate()

Pure virtual destructor. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::computeJacobian "bool FiniteDifferenceColoringWithUpdate::computeJacobian(const
Epetra_Vector &x, Epetra_Operator &Jac)

Computes (or updates) the Jacobian given the specified input vector,
x. Returns true if computation was successful. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::computeJacobian "bool FiniteDifferenceColoringWithUpdate::computeJacobian(const
Epetra_Vector &x)

Computes (or updates) Jacobian given the specified input vector, x.
Returns true if computation was successful. ";

%feature("docstring")
NOX::Epetra::FiniteDifferenceColoringWithUpdate::SetProbingDiagnostics
"void
NOX::Epetra::FiniteDifferenceColoringWithUpdate::SetProbingDiagnostics(bool
activate)

Disable/Enable the (low computational cost) probing diagnostics. ";


// File: classNOX_1_1StatusTest_1_1FiniteValue.xml
%feature("docstring") NOX::StatusTest::FiniteValue "

Failure test based on whether the norm of a vector has a finite value.

This test returns NOX::StatusTest::Failed if the norm of a vector is
calssified as a NaN or Inf. Otherwise, it returns
NOX::StatusTest::Unconverged. The user can specify whether to use the
F vector or the solution vector from the current solution group for
the check. NOX does not have access to vector entries so the number
used in the NaN/Inf check is based on the norm of a vector.

If checkStatus is called with the type set to NOX::StatusTest::None,
then the status is set to NOX::Status::Unevaluated and returned.

C++ includes: NOX_StatusTest_FiniteValue.H ";

%feature("docstring")  NOX::StatusTest::FiniteValue::FiniteValue "NOX::StatusTest::FiniteValue::FiniteValue(VectorType v=FVector,
NOX::Abstract::Vector::NormType n=NOX::Abstract::Vector::TwoNorm)

Constructor. Specify which vector to check and with what norm to use.
";

%feature("docstring")  NOX::StatusTest::FiniteValue::~FiniteValue "NOX::StatusTest::FiniteValue::~FiniteValue()

Destructor. ";

%feature("docstring")  NOX::StatusTest::FiniteValue::checkStatus "NOX::StatusTest::StatusType
NOX::StatusTest::FiniteValue::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
NOX::StatusType::None. If the test is skipped, then the status should
be set to NOX::StatusTest::Unevaluated. ";

%feature("docstring")  NOX::StatusTest::FiniteValue::getStatus "NOX::StatusTest::StatusType NOX::StatusTest::FiniteValue::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::FiniteValue::print "ostream &
NOX::StatusTest::FiniteValue::print(ostream &stream, int indent=0)
const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  NOX::StatusTest::FiniteValue::finiteNumberTest
"int NOX::StatusTest::FiniteValue::finiteNumberTest(double x) const

The finite number test algorithm.

Autoconf will test to see if the compiler implements the isnan() and
isinf() functions in the cmath or math.h headers. If so, we will use
these. If not, we supply a default implementation. The default
implementation is only guaranteed to work if the code is IEEE 748/754
compliant. The checks for isnan and isinf are separate because
compilers like the old sgi platforms support one but not the other.
See bug 2019 for more details.

This method is public so that other objects (solvers, line searches,
and directions) can use this test on their own values.

Return Values: 0 = Finite Number

-1 = NaN

-2 = Inf ";


// File: classNOX_1_1Multiphysics_1_1Solver_1_1FixedPointBased.xml
%feature("docstring") NOX::Multiphysics::Solver::FixedPointBased "

Nonlinear solver based on a line search (i.e., damping).

Solves $F(x)=0$ using an iterative line-search-based method.

Each iteration, the solver does the following.

Compute a search direction $d$ via a NOX::Direction method

Compute a step length $\\\\lambda$ and update $x$ as $x_{\\\\rm new} =
x_{\\\\rm old} + \\\\lambda d$ via a NOX::LineSearch method.

The iterations progress until the status tests (see NOX::StatusTest)
determine either failure or convergence.

To support several line searches and status tests, this version of the
solver has a getStepSize() function that returns $\\\\lambda$.  Input
Parameters

The following parameter list entries are valid for this solver:

\"Line Search\" - Sublist of the line search parameters, passed to the
NOX::LineSearch::Manager constructor. Defaults to an empty list.

\"Direction\" - Sublist of the direction parameters, passed to the
NOX::Direction::Manager constructor. Defaults to an empty list.

\"Solver Options\" - Sublist of general solver options. \"User Defined
Pre/Post Operator\" is supported. See NOX::Parameter::PrePostOperator
for more details.

Output Parameters

Every time solve() is called, a sublist for output parameters called
\"Output\" will be created and contain the following parameters.

\"Output\":

\"Nonlinear Iterations\" - Number of nonlinear iterations

\"2-Norm of Residual\" - Two-norm of final residual

Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)

C++ includes: NOX_Multiphysics_Solver_FixedPointBased.H ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::FixedPointBased "NOX::Multiphysics::Solver::FixedPointBased::FixedPointBased(const
Teuchos::RCP< vector< Teuchos::RCP< NOX::Solver::Generic > > >
&solvers, const Teuchos::RCP<
NOX::Multiphysics::DataExchange::Interface > &interface, const
Teuchos::RCP< NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)

Constructor.

See reset(NOX::Abstract::Group&, NOX::StatusTest::Generic&,
Teuchos::ParameterList&) for description ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::~FixedPointBased "NOX::Multiphysics::Solver::FixedPointBased::~FixedPointBased()

Destructor. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::reset "bool
NOX::Multiphysics::Solver::FixedPointBased::reset(const Teuchos::RCP<
vector< Teuchos::RCP< NOX::Solver::Generic > > > &solvers, const
Teuchos::RCP< NOX::Multiphysics::DataExchange::Interface > &interface,
const Teuchos::RCP< NOX::StatusTest::Generic > &tests, const
Teuchos::RCP< Teuchos::ParameterList > &params)

Reset the nonlinear solver for a new solve.

Parameters:
-----------

tests:  Status tests to check for convergence or failure. These tests
will be modified by the solver.

params:  List of parameters. These parameters will be modified by the
solver.

All the objects passed to reset() will be modified.

The group object will be cloned via NOX::Abstract::Group::clone(), and
the vectors within will also be individually cloned via
NOX::Abstract::Vector::clone().

WARNING:  If the contents of grp, tests, or params are modified by the
calling program after calling reset(), then the behavior of iterate()
and solve() are completely undefined. To remedy this, call reset()
again with the modified objects. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::reset "void
NOX::Multiphysics::Solver::FixedPointBased::reset(const
NOX::Abstract::Vector &initialGuess, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests) ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::reset "void
NOX::Multiphysics::Solver::FixedPointBased::reset(const
NOX::Abstract::Vector &initialGuess)

reset methods inherited from NOX::Solver::Generic and needed here to
avoid hiding this overloaded virtual method ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getStatus "NOX::StatusTest::StatusType
NOX::Multiphysics::Solver::FixedPointBased::getStatus()

Check current convergence and failure status. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::step "NOX::StatusTest::StatusType
NOX::Multiphysics::Solver::FixedPointBased::step()

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::solve "NOX::StatusTest::StatusType
NOX::Multiphysics::Solver::FixedPointBased::solve()

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getSolutionGroup "const
NOX::Abstract::Group &
NOX::Multiphysics::Solver::FixedPointBased::getSolutionGroup() const

Return a reference to the current solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getPreviousSolutionGroup "const NOX::Abstract::Group &
NOX::Multiphysics::Solver::FixedPointBased::getPreviousSolutionGroup()
const

Return a reference to the previous solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getNumIterations "int
NOX::Multiphysics::Solver::FixedPointBased::getNumIterations() const

Get number of iterations. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getList "const
Teuchos::ParameterList &
NOX::Multiphysics::Solver::FixedPointBased::getList() const

Return a reference to the solver parameters. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getSolutionGroupPtr "virtual Teuchos::RCP< const NOX::Abstract::Group >
NOX::Multiphysics::Solver::FixedPointBased::getSolutionGroupPtr()
const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getPreviousSolutionGroupPtr
"Teuchos::RCP< const NOX::Abstract::Group >
NOX::Multiphysics::Solver::FixedPointBased::getPreviousSolutionGroupPtr()
const

Return a RCP to the previous solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::FixedPointBased::getListPtr "virtual
Teuchos::RCP< const Teuchos::ParameterList >
NOX::Multiphysics::Solver::FixedPointBased::getListPtr() const

Return a RCP to the solver parameters. ";


// File: classNOX_1_1LineSearch_1_1FullStep.xml
%feature("docstring") NOX::LineSearch::FullStep "

Simplest line search - always take the full step.

This line search can be called via NOX::LineSearch::Manager.

The following parameters can be specified in the \"Full Step\" sublist
of the \"Line Search\" sublist:

\"Full Step\" - length of a full step (defaults to 1.0)

C++ includes: NOX_LineSearch_FullStep.H ";

%feature("docstring")  NOX::LineSearch::FullStep::FullStep "FullStep::FullStep(const Teuchos::RCP< NOX::GlobalData > &gd,
Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::LineSearch::FullStep::~FullStep "FullStep::~FullStep()

Destructor. ";

%feature("docstring")  NOX::LineSearch::FullStep::reset "bool
FullStep::reset(const Teuchos::RCP< NOX::GlobalData > &gd,
Teuchos::ParameterList &params) ";

%feature("docstring")  NOX::LineSearch::FullStep::compute "bool
FullStep::compute(NOX::Abstract::Group &newgrp, double &step, const
NOX::Abstract::Vector &dir, const NOX::Solver::Generic &s)

Perform a line search.

On input:

Parameters:
-----------

grp:  The initial solution vector, $x_{\\\\rm old}$.

dir:  A vector of directions to be used in the line search, $d$.

s:  The nonlinear solver.

On output:

Parameters:
-----------

step:  The distance the direction was scaled, $ \\\\lambda $.

grp:  The grp is updated with a new solution, $ x_{\\\\rm new} $,
resulting from the linesearch. Normally, for a single direction line
search, this is computed as:

\\\\[ x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d. \\\\]

Ideally, $ \\\\|F(x_{\\\\rm new})\\\\| < \\\\|F(x_{\\\\rm old})\\\\| $
(e.g the final direction is a descent direction).

Note that the dir object is a std::vector. For typical line searches
as described in the above equation, this vector is of size one. We
have used a std::vector to allow for special cases of multi-
directional line searches such as the Bader/Schnabel curvillinear line
search.

Return value is true for a successful line search computation. ";


// File: classNOX_1_1Direction_1_1Generic.xml
%feature("docstring") NOX::Direction::Generic "

Generic direction interface

Generic interface for calculating a search direction, $d$, to be used
in updating the iterate.

C++ includes: NOX_Direction_Generic.H ";

%feature("docstring")  NOX::Direction::Generic::Generic "NOX::Direction::Generic::Generic()

Constructor.

Constructors of derived objects should look like reset(). ";

%feature("docstring")  NOX::Direction::Generic::~Generic "virtual
NOX::Direction::Generic::~Generic()

Destructor. ";

%feature("docstring")  NOX::Direction::Generic::reset "virtual bool
NOX::Direction::Generic::reset(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params)=0

Reset direction based on possibly new parameters. ";

%feature("docstring")  NOX::Direction::Generic::compute "virtual bool
NOX::Direction::Generic::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::Generic &solver)=0

Compute the direction vector, dir, for a specific method given the
current group, grp.

The grp is not const so that we can compute the F vector, the Jacobian
matrix, the Newton vector, and so on.

Const access to the solver is used for getting additional information
such as the past solution, the iteration number, and so on. ";

%feature("docstring")  NOX::Direction::Generic::compute "bool
NOX::Direction::Generic::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::LineSearchBased &solver)

Same as compute( NOX::Abstract::Vector&, NOX::Abstract::Group&, const
NOX::Solver::Generic&).

Enables direct support for line search based solvers for the purpose
of efficiency since the LineSearchBased object has a getStep()
function that some directions require.

If it is not redefined in the derived class, it will just call the
compute with the NOX::Solver::Generic argument. ";


// File: classNOX_1_1Multiphysics_1_1Solver_1_1Generic.xml
%feature("docstring") NOX::Multiphysics::Solver::Generic "

Abstract nonlinear solver method interface.

Defines the type of access methods into the iterative nonlinear
solvers.

Instantiate or reset() the solver.

Find the solution via solve() or perform a single iterations via
iterate().

Get information about the current solver state via getSolutionGroup(),
getPreviousSolutionGroup(), getNumIterations(), and getParameterList()
--- particularily useful for NOX::StatusTest methods.

Get the current status of the solver via getStatus().

C++ includes: NOX_Multiphysics_Solver_Generic.H ";

%feature("docstring")  NOX::Multiphysics::Solver::Generic::Generic "NOX::Multiphysics::Solver::Generic::Generic()

Constructor (does nothing). ";

%feature("docstring")  NOX::Multiphysics::Solver::Generic::~Generic "virtual NOX::Multiphysics::Solver::Generic::~Generic()

Destructor (does nothing). ";

%feature("docstring")  NOX::Multiphysics::Solver::Generic::reset "virtual bool NOX::Multiphysics::Solver::Generic::reset(const
Teuchos::RCP< vector< Teuchos::RCP< NOX::Solver::Generic > > >
&solvers, const Teuchos::RCP<
NOX::Multiphysics::DataExchange::Interface > &interface, const
Teuchos::RCP< NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)=0

Reset the nonlinear solver for a new solve.

Parameters:
-----------

tests:  Status tests to check for convergence or failure. These tests
will be modified by the solver.

params:  List of parameters. These parameters will be modified by the
solver.

All the objects passed to reset() will be modified.

The group object will be cloned via NOX::Abstract::Group::clone(), and
the vectors within will also be individually cloned via
NOX::Abstract::Vector::clone().

WARNING:  If the contents of grp, tests, or params are modified by the
calling program after calling reset(), then the behavior of iterate()
and solve() are completely undefined. To remedy this, call reset()
again with the modified objects. ";

%feature("docstring")  NOX::Multiphysics::Solver::Generic::reset "virtual void NOX::Multiphysics::Solver::Generic::reset(const
NOX::Abstract::Vector &initialGuess)=0

reset methods inherited from NOX::Solver::Generic and needed here to
avoid hiding this overloaded virtual method ";

%feature("docstring")  NOX::Multiphysics::Solver::Generic::reset "virtual void NOX::Multiphysics::Solver::Generic::reset(const
NOX::Abstract::Vector &initialGuess, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests)=0 ";


// File: classNOX_1_1Solver_1_1Generic.xml
%feature("docstring") NOX::Solver::Generic "

Abstract nonlinear solver method interface.

Defines the type of access methods into the iterative nonlinear
solvers.

Instantiate or reset() the solver.

Find the solution via solve() or perform a single iterations via
iterate().

Get information about the current solver state via getSolutionGroup(),
getPreviousSolutionGroup(), getNumIterations(), and getList() ---
particularily useful for NOX::StatusTest methods.

Get the current status of the solver via getStatus().

C++ includes: NOX_Solver_Generic.H ";

%feature("docstring")  NOX::Solver::Generic::getSolutionGroupPtr "virtual Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::Generic::getSolutionGroupPtr() const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Solver::Generic::getPreviousSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::Generic::getPreviousSolutionGroupPtr() const

Return a RCP to the previous solution group. ";

%feature("docstring")  NOX::Solver::Generic::getListPtr "virtual
Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::Generic::getListPtr() const

Return a RCP to the solver parameters. ";

%feature("docstring")  NOX::Solver::Generic::Generic "NOX::Solver::Generic::Generic()

Constructor (does nothing). ";

%feature("docstring")  NOX::Solver::Generic::~Generic "virtual
NOX::Solver::Generic::~Generic()

Destructor (does nothing). ";

%feature("docstring")  NOX::Solver::Generic::reset "virtual void
NOX::Solver::Generic::reset(const NOX::Abstract::Vector
&initial_guess)=0

Resets the solver and sets a new initial guess. ";

%feature("docstring")  NOX::Solver::Generic::reset "virtual void
NOX::Solver::Generic::reset(const NOX::Abstract::Vector
&initial_guess, const Teuchos::RCP< NOX::StatusTest::Generic >
&test)=0

Resets the solver, sets a new status test, and sets a new initial
guess. ";

%feature("docstring")  NOX::Solver::Generic::getStatus "virtual
NOX::StatusTest::StatusType NOX::Solver::Generic::getStatus()=0

Check current convergence and failure status. ";

%feature("docstring")  NOX::Solver::Generic::step "virtual
NOX::StatusTest::StatusType NOX::Solver::Generic::step()=0

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")  NOX::Solver::Generic::solve "virtual
NOX::StatusTest::StatusType NOX::Solver::Generic::solve()=0

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")  NOX::Solver::Generic::getSolutionGroup "virtual const NOX::Abstract::Group&
NOX::Solver::Generic::getSolutionGroup() const =0

Return a reference to the current solution group. ";

%feature("docstring")  NOX::Solver::Generic::getPreviousSolutionGroup
"virtual const NOX::Abstract::Group&
NOX::Solver::Generic::getPreviousSolutionGroup() const =0

Return a reference to the previous solution group. ";

%feature("docstring")  NOX::Solver::Generic::getNumIterations "virtual int NOX::Solver::Generic::getNumIterations() const =0

Get number of iterations. ";

%feature("docstring")  NOX::Solver::Generic::getList "virtual const
Teuchos::ParameterList& NOX::Solver::Generic::getList() const =0

Return a reference to the solver parameters. ";


// File: classNOX_1_1LineSearch_1_1Generic.xml
%feature("docstring") NOX::LineSearch::Generic "

Base class line search interface.

Every line search should respect the following Parameter:

\"Max Iters\" - maximum number of iterations (i.e., RHS computations)

C++ includes: NOX_LineSearch_Generic.H ";

%feature("docstring")  NOX::LineSearch::Generic::Generic "NOX::LineSearch::Generic::Generic()

Default constructor. ";

%feature("docstring")  NOX::LineSearch::Generic::~Generic "virtual
NOX::LineSearch::Generic::~Generic()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Generic::compute "virtual
bool NOX::LineSearch::Generic::compute(NOX::Abstract::Group &grp,
double &step, const NOX::Abstract::Vector &dir, const
NOX::Solver::Generic &s)=0

Perform a line search.

On input:

Parameters:
-----------

grp:  The initial solution vector, $x_{\\\\rm old}$.

dir:  A vector of directions to be used in the line search, $d$.

s:  The nonlinear solver.

On output:

Parameters:
-----------

step:  The distance the direction was scaled, $ \\\\lambda $.

grp:  The grp is updated with a new solution, $ x_{\\\\rm new} $,
resulting from the linesearch. Normally, for a single direction line
search, this is computed as:

\\\\[ x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d. \\\\]

Ideally, $ \\\\|F(x_{\\\\rm new})\\\\| < \\\\|F(x_{\\\\rm old})\\\\| $
(e.g the final direction is a descent direction).

Note that the dir object is a std::vector. For typical line searches
as described in the above equation, this vector is of size one. We
have used a std::vector to allow for special cases of multi-
directional line searches such as the Bader/Schnabel curvillinear line
search.

Return value is true for a successful line search computation. ";


// File: classNOX_1_1StatusTest_1_1Generic.xml
%feature("docstring") NOX::StatusTest::Generic "

Generic status test to check for convergence or failure of the
nonlinear solver.

C++ includes: NOX_StatusTest_Generic.H ";

%feature("docstring")  NOX::StatusTest::Generic::Generic "NOX::StatusTest::Generic::Generic()

Constructor. ";

%feature("docstring")  NOX::StatusTest::Generic::~Generic "virtual
NOX::StatusTest::Generic::~Generic()

Destructor. ";

%feature("docstring")  NOX::StatusTest::Generic::checkStatus "virtual
NOX::StatusTest::StatusType
NOX::StatusTest::Generic::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)=0

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
NOX::StatusType::None. If the test is skipped, then the status should
be set to NOX::StatusTest::Unevaluated. ";

%feature("docstring")  NOX::StatusTest::Generic::getStatus "virtual
NOX::StatusTest::StatusType NOX::StatusTest::Generic::getStatus()
const =0

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::Generic::print "virtual
ostream& NOX::StatusTest::Generic::print(ostream &stream, int
indent=0) const =0

Output formatted description of stopping test to output stream. ";


// File: classNOX_1_1MeritFunction_1_1Generic.xml
%feature("docstring") NOX::MeritFunction::Generic "

Base class to support a user defined merit function that can be passed
to line searches and directions through the parameter list.

This class allows the user to define their own merit function for use
in a line search. Each line search type will specify in it's input
parameter list if it supports this functionality.

To create and use a user defined merit function:

Create a merit function that derives from
NOX::Parameter::MeritFunction. For example, the merit function Foo
might be defined as shown below.

Create the appropriate entries in the parameter list, as follows.

C++ includes: NOX_MeritFunction_Generic.H ";

%feature("docstring")  NOX::MeritFunction::Generic::Generic "NOX::MeritFunction::Generic::Generic()

Default Constructor. ";

%feature("docstring")  NOX::MeritFunction::Generic::~Generic "virtual
NOX::MeritFunction::Generic::~Generic()

Destructor. ";

%feature("docstring")  NOX::MeritFunction::Generic::computef "virtual
double NOX::MeritFunction::Generic::computef(const
NOX::Abstract::Group &grp) const =0

Computes the merit function, $ f(x) $. ";

%feature("docstring")  NOX::MeritFunction::Generic::computeGradient "virtual void NOX::MeritFunction::Generic::computeGradient(const
NOX::Abstract::Group &group, NOX::Abstract::Vector &result) const =0

Computes the gradient of the merit function, $ \\\\nabla f $, and
returns the result in the result vector. ";

%feature("docstring")  NOX::MeritFunction::Generic::computeSlope "virtual double NOX::MeritFunction::Generic::computeSlope(const
NOX::Abstract::Vector &dir, const NOX::Abstract::Group &grp) const =0

Computes the inner product of the given direction and the gradient
associated with the merit function. Returns the steepest descent
direction in the result vector.

Calculates and returns $ \\\\zeta $: \\\\[ \\\\zeta = \\\\nabla f(x)^T
d \\\\]

Here $d$ represents the input parameter dir and $\\\\nabla f(x)$ is
the gradient of the merit function. ";

%feature("docstring")
NOX::MeritFunction::Generic::computeQuadraticModel "virtual double
NOX::MeritFunction::Generic::computeQuadraticModel(const
NOX::Abstract::Vector &dir, const NOX::Abstract::Group &grp) const =0

Compute the quadratic model, $ m(d) $, for the given merit function.

Computes and returns $ m(d) $: \\\\[ m(d) = f(x) + \\\\nabla f(x)^T d
+ d^T \\\\nabla^2 f(x) d + d^T \\\\mathbf{B} d \\\\]

Here $d$ represents the input parameter dir. $ B $ is the Hessian of
the merit function, $\\\\nabla^2 f(x)$, but can be approximated with
the restriction that it is a symmetric and has uniform boundedness in
the iterate sequence (see J. Nocedal and S. J. Wright, \"Numerical
Optimization\", Springer, 1999. Chapters 4 and 6). ";

%feature("docstring")
NOX::MeritFunction::Generic::computeQuadraticMinimizer "virtual void
NOX::MeritFunction::Generic::computeQuadraticMinimizer(const
NOX::Abstract::Group &grp, NOX::Abstract::Vector &result) const =0

Computes the vector in the steepest descent direction that minimizes
the quadratic model.

The quadratic model is defined as: \\\\[ m(d) = f(x) + \\\\nabla
f(x)^T d + d^T \\\\nabla^2 f(x) d + d^T \\\\mathbf{B} d \\\\]

where $ B $ is ideally the Hessian of the merit function, $\\\\nabla^2
f(x)$, but can be approximated with the restriction that it is a
symmetric and has uniform boundedness in the iterate sequence (see J.
Nocedal and S. J. Wright, \"Numerical Optimization\", Springer, 1999.
Chapters 4 and 6).

The result vector should be computed as: \\\\[ result =
-\\\\frac{\\\\nabla f^T \\\\nabla f}{\\\\nabla f^T B \\\\nabla f}
\\\\nabla f \\\\] ";

%feature("docstring")  NOX::MeritFunction::Generic::name "virtual
const string& NOX::MeritFunction::Generic::name() const =0

Returns the name of the merit function. ";


// File: classNOX_1_1GlobalData.xml
%feature("docstring") NOX::GlobalData "

Container class to hold \"global\" NOX objects.

GlobalData is a container class that holds ref-count pointers to
\"global\" objects, i.e., objects that nearly every NOX object will
need access to. By putting them all in one container class, the
container class can be stored in each NOX object, and if a new global
object is needed, it can be added here without modifying the rest of
the code. This is an alternative to true global or static objects
which are note safe in many contexts (threading). In particular, this
approach allows multiple NOX \"invocations\" to be in memory at the
same time.

C++ includes: NOX_GlobalData.H ";

%feature("docstring")  NOX::GlobalData::GlobalData "NOX::GlobalData::GlobalData(const Teuchos::RCP< Teuchos::ParameterList
> &noxParams)

Consturctor using the top level NOX parameter list. ";

%feature("docstring")  NOX::GlobalData::GlobalData "NOX::GlobalData::GlobalData(const Teuchos::RCP< NOX::Utils > &utils,
const Teuchos::RCP< NOX::MeritFunction::Generic > &mf)

Constructor taking a ref-count pointer to each global object. ";

%feature("docstring")  NOX::GlobalData::~GlobalData "NOX::GlobalData::~GlobalData()

Destructor. ";

%feature("docstring")  NOX::GlobalData::getUtils "Teuchos::RCP<
NOX::Utils > NOX::GlobalData::getUtils() const

Returns the print utils object. ";

%feature("docstring")  NOX::GlobalData::getMeritFunction "Teuchos::RCP< NOX::MeritFunction::Generic >
NOX::GlobalData::getMeritFunction() const

Returns the merit function object. ";

%feature("docstring")  NOX::GlobalData::getNoxParameterList "Teuchos::RCP< Teuchos::ParameterList >
NOX::GlobalData::getNoxParameterList() const

Returns the top-level nox parameter list input by the user.

This list is kept in global data so that any sublists of the main
parameters list that objects may keep a refernece to is still valid.
The line searches and directions usually store data in an output
sublist for the users to query. These sublists are NOT wrapped in
reference counted smart pointers, so if the base list is deleted, the
references to the sublist will no longer be valid. To remedy this, any
object that stores a reference to a sublist should also store the
global data object. ";


// File: classNOX_1_1Epetra_1_1Group.xml
%feature("docstring") NOX::Epetra::Group "

Concrete implementation of NOX::Abstract::Group for Trilinos/Epetra.

This group is set up to use the linear algebra services provided
through the Trilinos/Epetra package with AztecOO for the linear
solver.

C++ includes: NOX_Epetra_Group.H ";

/*  Vectors  */

/* */

/*  IsValid flags  */

/* True if the current solution is up-to-date with respect to the
currect xVector.

*/

/*  Shared Operators  */

/* */

/*  "Compute" functions.  */

/* */

%feature("docstring")  NOX::Epetra::Group::setX "void
Group::setX(const NOX::Epetra::Vector &y) ";

%feature("docstring")  NOX::Epetra::Group::setX "void
Group::setX(const NOX::Abstract::Vector &y)

Set the solution vector x to y.

This should invalidate the function value, Jacobian, gradient, and
Newton direction.

Throw an error if the copy fails.

Reference to this object ";

%feature("docstring")  NOX::Epetra::Group::computeX "void
Group::computeX(const Group &grp, const NOX::Epetra::Vector &d, double
step) ";

%feature("docstring")  NOX::Epetra::Group::computeX "void
Group::computeX(const NOX::Abstract::Group &grp, const
NOX::Abstract::Vector &d, double step)

Compute x = grp.x + step * d.

Let $x$ denote this group's solution vector. Let $\\\\hat x$ denote
the result of grp.getX(). Then set \\\\[ x = \\\\hat x +
\\\\mbox{step} \\\\; d. \\\\]

This should invalidate the function value, Jacobian, gradient, and
Newton direction.

Throw an error if the copy fails.

Reference to this object ";

%feature("docstring")  NOX::Epetra::Group::computeF "Abstract::Group::ReturnType Group::computeF()

Compute and store F(x).

It's generally useful to also compute and store the 2-norm of F(x) at
this point for later access by the getNormF() function.

NOX::Abstract::Group::Failed - If the computation fails in any way

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Epetra::Group::computeJacobian "Abstract::Group::ReturnType Group::computeJacobian()

Compute and store Jacobian.

Recall that \\\\[ F(x) = \\\\left[ \\\\begin{array}{c} F_1(x) \\\\\\\\
F_2(x) \\\\\\\\ \\\\vdots \\\\\\\\ F_n(x) \\\\\\\\ \\\\end{array}
\\\\right]. \\\\]

The Jacobian is denoted by $J$ and defined by \\\\[ J_{ij} =
\\\\frac{\\\\partial F_i}{\\\\partial x_j} (x). \\\\]

If this is a shared object, this group should taken ownership of the
Jacobian before it computes it.

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::Failed - If the computation fails in any other
way

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Epetra::Group::computeGradient "Abstract::Group::ReturnType Group::computeGradient()

Compute and store gradient.

We can pose the nonlinear equation problem $F(x) = 0$ as an
optimization problem as follows: \\\\[ \\\\min f(x) \\\\equiv
\\\\frac{1}{2} \\\\|F(x)\\\\|_2^2. \\\\]

In that case, the gradient (of $f$) is defined as \\\\[ g \\\\equiv
J^T F. \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If either $F$ or $J$ has not
been computed

NOX::Abstract::Group::Failed - If the computation fails in any other
way

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Epetra::Group::computeNewton "Abstract::Group::ReturnType
Group::computeNewton(Teuchos::ParameterList &params)

Compute the Newton direction, using parameters for the linear solve.

The Newton direction is the solution, s, of \\\\[ J s = -F. \\\\]

The parameters are from the \"Linear %Solver\" sublist of the
\"Direction\" sublist that is passed to solver during construction.

The \"Tolerance\" parameter may be added/modified in the sublist of
\"Linear Solver\" parameters that is passed into this function. The
solution should be such that \\\\[ \\\\frac{\\\\| J s - (-F)
\\\\|_2}{\\\\max \\\\{ 1, \\\\|F\\\\|_2\\\\} } < \\\\mbox{Tolerance}
\\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If either $F$ or $J$ has not
been computed

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails in any other
way

NOX::Abstract::Group::Ok - Otherwise ";

/*  Jacobian operations.  */

/* Operations using the Jacobian matrix. These may not be defined in
matrix-free scenarios.

*/

%feature("docstring")  NOX::Epetra::Group::applyJacobian "Abstract::Group::ReturnType Group::applyJacobian(const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const ";

%feature("docstring")  NOX::Epetra::Group::applyJacobian "Abstract::Group::ReturnType Group::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian to the given input vector and puts the answer in the
result.

Computes \\\\[ v = J u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector.

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If the Jacobian $J$ has not been
computed

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Epetra::Group::applyJacobianTranspose "Abstract::Group::ReturnType Group::applyJacobianTranspose(const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const ";

%feature("docstring")  NOX::Epetra::Group::applyJacobianTranspose "Abstract::Group::ReturnType Group::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian-Transpose to the given input vector and puts the
answer in the result.

Computes \\\\[ v = J^T u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector.

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If $J$ has not been computed

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Epetra::Group::applyJacobianInverse "Abstract::Group::ReturnType
Group::applyJacobianInverse(Teuchos::ParameterList &params, const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const

Applies the inverse of the Jacobian matrix to the given input vector
and puts the answer in result.

Computes \\\\[ v = J^{-1} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The \"Tolerance\" parameter specifies that the solution should be such
that \\\\[ \\\\frac{\\\\| J v - u \\\\|_2}{\\\\max \\\\{ 1,
\\\\|u\\\\|_2\\\\} } < \\\\mbox{Tolerance} \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If $J$ has not been computed

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise

The parameter \"Tolerance\" may be added/modified in the list of
parameters - this is the ideal solution tolerance for an iterative
linear solve.

The parameter \"Reuse Preconditioner\" is a boolean that tells the
group to turn off control of preconditioner recalculation. This is a
dangerous flag but can really speed the computations if the user knows
what they are doing. Toggling this flag is left to the user (ideally
it should be done through a status test). Defaults to false. ";

%feature("docstring")  NOX::Epetra::Group::applyJacobianInverse "Abstract::Group::ReturnType
Group::applyJacobianInverse(Teuchos::ParameterList &params, const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies the inverse of the Jacobian matrix to the given input vector
and puts the answer in result.

Computes \\\\[ v = J^{-1} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The \"Tolerance\" parameter specifies that the solution should be such
that \\\\[ \\\\frac{\\\\| J v - u \\\\|_2}{\\\\max \\\\{ 1,
\\\\|u\\\\|_2\\\\} } < \\\\mbox{Tolerance} \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If $J$ has not been computed

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise

The parameter \"Tolerance\" may be added/modified in the list of
parameters - this is the ideal solution tolerance for an iterative
linear solve. ";

%feature("docstring")  NOX::Epetra::Group::applyRightPreconditioning "Abstract::Group::ReturnType Group::applyRightPreconditioning(bool
useTranspose, Teuchos::ParameterList &params, const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const ";

%feature("docstring")  NOX::Epetra::Group::applyRightPreconditioning "Abstract::Group::ReturnType Group::applyRightPreconditioning(bool
useTranspose, Teuchos::ParameterList &params, const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Apply right preconditiong to the given input vector.

Let $M$ be a right preconditioner for the Jacobian $J$; in other
words, $M$ is a matrix such that \\\\[ JM \\\\approx I. \\\\]

Compute \\\\[ u = M^{-1} v, \\\\] where $u$ is the input vector and
$v$ is the result vector.

If useTranspose is true, then the transpose of the preconditioner is
applied: \\\\[ u = {M^{-1}}^T v, \\\\] The transpose preconditioner is
currently only required for Tensor methods.

The \"Tolerance\" parameter specifies that the solution should be such
that \\\\[ \\\\frac{\\\\| M v - u \\\\|_2}{\\\\max \\\\{ 1,
\\\\|u\\\\|_2\\\\} } < \\\\mbox{Tolerance} \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise

The parameters are from the \"Linear %Solver\" sublist of the
\"Direction\" sublist that is passed to solver during construction. ";

/*  "Is" functions  */

/* Checks to see if various objects have been computed. Returns true
if the corresponding \"compute\" function has been called since the
last update to the solution vector (via instantiation or computeX).

*/

%feature("docstring")  NOX::Epetra::Group::isF "bool Group::isF()
const

Return true if F is valid. ";

%feature("docstring")  NOX::Epetra::Group::isJacobian "bool
Group::isJacobian() const

Return true if the Jacobian is valid.

Default implementation in NOX::Abstract::Group returns false. ";

%feature("docstring")  NOX::Epetra::Group::isGradient "bool
Group::isGradient() const

Return true if the gradient is valid.

Default implementation in NOX::Abstract::Group returns false. ";

%feature("docstring")  NOX::Epetra::Group::isNewton "bool
Group::isNewton() const

Return true if the Newton direction is valid.

Default implementation in NOX::Abstract::Group returns false. ";

%feature("docstring")  NOX::Epetra::Group::isNormNewtonSolveResidual "bool Group::isNormNewtonSolveResidual() const

Returns true if the value of the Norm of the linear model for a full
Newton step ||Js + f|| is valid with respect to the current solution
vector. ";

%feature("docstring")  NOX::Epetra::Group::isPreconditioner "bool
Group::isPreconditioner() const

Returns true if an explicitly constructed preconditioner exists (i.e.
one that is computed and saved for further use in multiple calls to
applyRightPreconditioner). ";

%feature("docstring")  NOX::Epetra::Group::isConditionNumber "bool
Group::isConditionNumber() const

Returns true if the condition number has been computed. ";

/*  "Get" functions  */

/* Note that these function do not check whether or not the vectors
are valid. Must use the \"Is\" functions for that purpose.

*/

%feature("docstring")  NOX::Epetra::Group::getX "const
Abstract::Vector & Group::getX() const

Return solution vector. ";

%feature("docstring")  NOX::Epetra::Group::getF "const
Abstract::Vector & Group::getF() const

Return F(x). ";

%feature("docstring")  NOX::Epetra::Group::getNormF "double
Group::getNormF() const

Return 2-norm of F(x).

In other words, \\\\[ \\\\sqrt{\\\\sum_{i=1}^n F_i^2} \\\\] ";

%feature("docstring")  NOX::Epetra::Group::getGradient "const
Abstract::Vector & Group::getGradient() const

Return gradient. ";

%feature("docstring")  NOX::Epetra::Group::getNewton "const
Abstract::Vector & Group::getNewton() const

Return Newton direction. ";

%feature("docstring")  NOX::Epetra::Group::getXPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Epetra::Group::getXPtr() const

Return RCP to solution vector. ";

%feature("docstring")  NOX::Epetra::Group::getFPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Epetra::Group::getFPtr() const

Return RCP to F(x). ";

%feature("docstring")  NOX::Epetra::Group::getGradientPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Epetra::Group::getGradientPtr() const

Return RCP to gradient. ";

%feature("docstring")  NOX::Epetra::Group::getNewtonPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Epetra::Group::getNewtonPtr() const

Return RCP to Newton direction. ";

%feature("docstring")
NOX::Epetra::Group::getNormLastLinearSolveResidual "Abstract::Group::ReturnType
NOX::Epetra::Group::getNormLastLinearSolveResidual(double &residual)
const

Returns the 2-norm of the residual of the linear model used in the
Newton solve computation, ||Js+f||. This does not account for line
search adjustments to the step length! ";

%feature("docstring")  NOX::Epetra::Group::Group "Group::Group(Teuchos::ParameterList &printingParams, const
Teuchos::RCP< NOX::Epetra::Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess)

Constructor with NO linear system (VERY LIMITED).

WARNING: If this constructor is used, then methods that require a
Jacobian or preconditioning will not be available. You will be limited
to simple algorithms like nonlinear-CG with no preconditioning. ";

%feature("docstring")  NOX::Epetra::Group::Group "Group::Group(Teuchos::ParameterList &printingParams, const
Teuchos::RCP< NOX::Epetra::Interface::Required > &i, const
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linSys)

Standard Constructor. ";

%feature("docstring")  NOX::Epetra::Group::Group "Group::Group(const
NOX::Epetra::Group &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. If type is DeepCopy, takes ownership of valid shared
linear system. ";

%feature("docstring")  NOX::Epetra::Group::~Group "Group::~Group()

Destructor. ";

%feature("docstring")  NOX::Epetra::Group::clone "Teuchos::RCP<
NOX::Abstract::Group > Group::clone(CopyType type=DeepCopy) const

Create a new Group of the same derived type as this one by cloning
this one, and return a ref count pointer to the new group.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\" (only the memory is allocated, the values are
not copied into the vectors and Jacobian). Returns NULL if clone is
not supported.

Any shared data should have its ownership transfered to this group
from the source for a NOX::DeepCopy. ";

%feature("docstring")  NOX::Epetra::Group::getRequiredInterface "Teuchos::RCP< NOX::Epetra::Interface::Required >
Group::getRequiredInterface()

Return the userInterface. ";

%feature("docstring")  NOX::Epetra::Group::getLinearSystem "Teuchos::RCP< const NOX::Epetra::LinearSystem >
Group::getLinearSystem() const

Return the Linear System. ";

%feature("docstring")  NOX::Epetra::Group::getLinearSystem "Teuchos::RCP< NOX::Epetra::LinearSystem > Group::getLinearSystem()

Return the Linear System. ";

%feature("docstring")
NOX::Epetra::Group::computeJacobianConditionNumber "Abstract::Group::ReturnType
NOX::Epetra::Group::computeJacobianConditionNumber(int maxIters,
double tolerance, int krylovSubspaceSize=100, bool printOutput=false)
";

%feature("docstring")  NOX::Epetra::Group::getJacobianConditionNumber
"double NOX::Epetra::Group::getJacobianConditionNumber() const

Returns the condition number of the Jacobian matrix. ";

%feature("docstring")
NOX::Epetra::Group::disableLinearResidualComputation "void
NOX::Epetra::Group::disableLinearResidualComputation(const bool
disableChoice)

Sets option to disable linear resid computation. If disabled, this
saves on a MatVec per Newton but disallows inexact Newton methods. ";


// File: classNOX_1_1Multiphysics_1_1Group.xml
%feature("docstring") NOX::Multiphysics::Group "

NOX pure abstract interface to a \"group\"; i.e., a solution vector
and the corresponding F-vector, Jacobian matrix, gradient vector, and
Newton vector.

This class is a member of the namespace NOX::Abstract.

The user should implement their own concrete implementation of this
class or use one of the implementations provided by us. Typically the
implementation is also tied to a particular NOX::Abstract::Vector
implementation.

The group may be implemented so that multiple groups can share
underlying memory space. This is particularly important when it comes
to the Jacobian, which is often to big to be replicated for every
group. Thus, we have included instructions on how shared data should
be treated for the operator=() and clone() functions.

C++ includes: NOX_Multiphysics_Group.H ";

/*  IsValid flags  */

/* True if the current solution is up-to-date with respect to the
currect xVector.

*/

%feature("docstring")  NOX::Multiphysics::Group::setX "void
NOX::Multiphysics::Group::setX(const NOX::Abstract::Vector &y)

Set the solution vector x to y.

This should invalidate the function value, Jacobian, gradient, and
Newton direction.

Throw an error if the copy fails.

Reference to this object ";

%feature("docstring")  NOX::Multiphysics::Group::computeX "void
NOX::Multiphysics::Group::computeX(const NOX::Abstract::Group &grp,
const NOX::Abstract::Vector &d, double step)

Compute x = grp.x + step * d.

Let $x$ denote this group's solution vector. Let $\\\\hat x$ denote
the result of grp.getX(). Then set \\\\[ x = \\\\hat x +
\\\\mbox{step} \\\\; d. \\\\]

This should invalidate the function value, Jacobian, gradient, and
Newton direction.

Throw an error if the copy fails.

Reference to this object ";

%feature("docstring")  NOX::Multiphysics::Group::computeF "NOX::Abstract::Group::ReturnType NOX::Multiphysics::Group::computeF()

Compute and store F(x).

It's generally useful to also compute and store the 2-norm of F(x) at
this point for later access by the getNormF() function.

NOX::Abstract::Group::Failed - If the computation fails in any way

NOX::Abstract::Group::Ok - Otherwise ";

/*  "Is" functions.  */

/* Checks to see if various objects have been computed. Returns true
if the corresponding \"compute\" function has been called since the
last change to the solution vector.

*/

%feature("docstring")  NOX::Multiphysics::Group::isF "bool
NOX::Multiphysics::Group::isF() const

Return true if F is valid. ";

/*  "Get" functions.  */

/* Note that these function do not check whether or not the vectors
are valid. Must use the \"Is\" functions for that purpose.

*/

%feature("docstring")  NOX::Multiphysics::Group::getX "const
NOX::Abstract::Vector & NOX::Multiphysics::Group::getX() const

Return solution vector. ";

%feature("docstring")  NOX::Multiphysics::Group::getF "const
NOX::Abstract::Vector & NOX::Multiphysics::Group::getF() const

Return F(x). ";

%feature("docstring")  NOX::Multiphysics::Group::getNormF "double
NOX::Multiphysics::Group::getNormF() const

Return 2-norm of F(x).

In other words, \\\\[ \\\\sqrt{\\\\sum_{i=1}^n F_i^2} \\\\] ";

%feature("docstring")  NOX::Multiphysics::Group::getGradient "const
NOX::Abstract::Vector & NOX::Multiphysics::Group::getGradient() const

Return gradient. ";

%feature("docstring")  NOX::Multiphysics::Group::getNewton "const
NOX::Abstract::Vector & NOX::Multiphysics::Group::getNewton() const

Return Newton direction. ";

%feature("docstring")  NOX::Multiphysics::Group::getXPtr "Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getXPtr() const

Return RCP to solution vector. ";

%feature("docstring")  NOX::Multiphysics::Group::getFPtr "Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getFPtr() const

Return RCP to F(x). ";

%feature("docstring")  NOX::Multiphysics::Group::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getGradientPtr() const

Return RCP to gradient. ";

%feature("docstring")  NOX::Multiphysics::Group::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getNewtonPtr() const

Return RCP to Newton direction. ";

%feature("docstring")  NOX::Multiphysics::Group::clone "Teuchos::RCP<
NOX::Abstract::Group > NOX::Multiphysics::Group::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new Group of the same derived type as this one by cloning
this one, and return a ref count pointer to the new group.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\" (only the memory is allocated, the values are
not copied into the vectors and Jacobian). Returns NULL if clone is
not supported.

Any shared data should have its ownership transfered to this group
from the source for a NOX::DeepCopy. ";

%feature("docstring")  NOX::Multiphysics::Group::Group "NOX::Multiphysics::Group::Group(const Teuchos::RCP< vector<
Teuchos::RCP< NOX::Solver::Generic > > > &solvers, const Teuchos::RCP<
NOX::StatusTest::Generic > &t, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Constructor.

Constructors for any derived object should always define a default
x-value so that getX() is always defined. ";

%feature("docstring")  NOX::Multiphysics::Group::Group "NOX::Multiphysics::Group::Group(const Group &grp, NOX::CopyType typ)

Constructor. ";

%feature("docstring")  NOX::Multiphysics::Group::~Group "NOX::Multiphysics::Group::~Group()

Destructor. ";


// File: classNOX_1_1Abstract_1_1Group.xml
%feature("docstring") NOX::Abstract::Group "

NOX pure abstract interface to a \"group\"; i.e., a solution vector
and the corresponding F-vector, Jacobian matrix, gradient vector, and
Newton vector.

This class is a member of the namespace NOX::Abstract.

The user should implement their own concrete implementation of this
class or use one of the implementations provided by us. Typically the
implementation is also tied to a particular NOX::Abstract::Vector
implementation.

The group may be implemented so that multiple groups can share
underlying memory space. This is particularly important when it comes
to the Jacobian, which is often to big to be replicated for every
group. Thus, we have included instructions on how shared data should
be treated for the operator=() and clone() functions.

C++ includes: NOX_Abstract_Group.H ";

%feature("docstring")  NOX::Abstract::Group::setX "virtual void
NOX::Abstract::Group::setX(const NOX::Abstract::Vector &y)=0

Set the solution vector x to y.

This should invalidate the function value, Jacobian, gradient, and
Newton direction.

Throw an error if the copy fails.

Reference to this object ";

%feature("docstring")  NOX::Abstract::Group::computeX "virtual void
NOX::Abstract::Group::computeX(const NOX::Abstract::Group &grp, const
NOX::Abstract::Vector &d, double step)=0

Compute x = grp.x + step * d.

Let $x$ denote this group's solution vector. Let $\\\\hat x$ denote
the result of grp.getX(). Then set \\\\[ x = \\\\hat x +
\\\\mbox{step} \\\\; d. \\\\]

This should invalidate the function value, Jacobian, gradient, and
Newton direction.

Throw an error if the copy fails.

Reference to this object ";

%feature("docstring")  NOX::Abstract::Group::computeF "virtual
NOX::Abstract::Group::ReturnType NOX::Abstract::Group::computeF()=0

Compute and store F(x).

It's generally useful to also compute and store the 2-norm of F(x) at
this point for later access by the getNormF() function.

NOX::Abstract::Group::Failed - If the computation fails in any way

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Abstract::Group::computeJacobian "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeJacobian()

Compute and store Jacobian.

Recall that \\\\[ F(x) = \\\\left[ \\\\begin{array}{c} F_1(x) \\\\\\\\
F_2(x) \\\\\\\\ \\\\vdots \\\\\\\\ F_n(x) \\\\\\\\ \\\\end{array}
\\\\right]. \\\\]

The Jacobian is denoted by $J$ and defined by \\\\[ J_{ij} =
\\\\frac{\\\\partial F_i}{\\\\partial x_j} (x). \\\\]

If this is a shared object, this group should taken ownership of the
Jacobian before it computes it.

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::Failed - If the computation fails in any other
way

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Abstract::Group::computeGradient "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeGradient()

Compute and store gradient.

We can pose the nonlinear equation problem $F(x) = 0$ as an
optimization problem as follows: \\\\[ \\\\min f(x) \\\\equiv
\\\\frac{1}{2} \\\\|F(x)\\\\|_2^2. \\\\]

In that case, the gradient (of $f$) is defined as \\\\[ g \\\\equiv
J^T F. \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If either $F$ or $J$ has not
been computed

NOX::Abstract::Group::Failed - If the computation fails in any other
way

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Abstract::Group::computeNewton "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::computeNewton(Teuchos::ParameterList &params)

Compute the Newton direction, using parameters for the linear solve.

The Newton direction is the solution, s, of \\\\[ J s = -F. \\\\]

The parameters are from the \"Linear %Solver\" sublist of the
\"Direction\" sublist that is passed to solver during construction.

The \"Tolerance\" parameter may be added/modified in the sublist of
\"Linear Solver\" parameters that is passed into this function. The
solution should be such that \\\\[ \\\\frac{\\\\| J s - (-F)
\\\\|_2}{\\\\max \\\\{ 1, \\\\|F\\\\|_2\\\\} } < \\\\mbox{Tolerance}
\\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If either $F$ or $J$ has not
been computed

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails in any other
way

NOX::Abstract::Group::Ok - Otherwise ";

/*  Jacobian operations.  */

/* Operations using the Jacobian matrix.

*/

%feature("docstring")  NOX::Abstract::Group::applyJacobian "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobian(const NOX::Abstract::Vector
&input, NOX::Abstract::Vector &result) const

Applies Jacobian to the given input vector and puts the answer in the
result.

Computes \\\\[ v = J u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector.

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If the Jacobian $J$ has not been
computed

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Abstract::Group::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian-Transpose to the given input vector and puts the
answer in the result.

Computes \\\\[ v = J^T u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector.

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If $J$ has not been computed

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise ";

%feature("docstring")  NOX::Abstract::Group::applyJacobianInverse "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies the inverse of the Jacobian matrix to the given input vector
and puts the answer in result.

Computes \\\\[ v = J^{-1} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The \"Tolerance\" parameter specifies that the solution should be such
that \\\\[ \\\\frac{\\\\| J v - u \\\\|_2}{\\\\max \\\\{ 1,
\\\\|u\\\\|_2\\\\} } < \\\\mbox{Tolerance} \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If $J$ has not been computed

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise

The parameter \"Tolerance\" may be added/modified in the list of
parameters - this is the ideal solution tolerance for an iterative
linear solve. ";

%feature("docstring")  NOX::Abstract::Group::applyRightPreconditioning
"NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioning(bool useTranspose,
Teuchos::ParameterList &params, const NOX::Abstract::Vector &input,
NOX::Abstract::Vector &result) const

Apply right preconditiong to the given input vector.

Let $M$ be a right preconditioner for the Jacobian $J$; in other
words, $M$ is a matrix such that \\\\[ JM \\\\approx I. \\\\]

Compute \\\\[ u = M^{-1} v, \\\\] where $u$ is the input vector and
$v$ is the result vector.

If useTranspose is true, then the transpose of the preconditioner is
applied: \\\\[ u = {M^{-1}}^T v, \\\\] The transpose preconditioner is
currently only required for Tensor methods.

The \"Tolerance\" parameter specifies that the solution should be such
that \\\\[ \\\\frac{\\\\| M v - u \\\\|_2}{\\\\max \\\\{ 1,
\\\\|u\\\\|_2\\\\} } < \\\\mbox{Tolerance} \\\\]

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::NotConverged - If the linear solve fails to
satisfy the \"Tolerance\" specified in params

NOX::Abstract::Group::Failed - If the computation fails

NOX::Abstract::Group::Ok - Otherwise

The parameters are from the \"Linear %Solver\" sublist of the
\"Direction\" sublist that is passed to solver during construction. ";

/*  Block Jacobian operations.  */

/* Operations using the Jacobian matrix.

*/

%feature("docstring")  NOX::Abstract::Group::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

applyJacobian for multiple right-hand sides

The default implementation here calls applyJacobian() for each right
hand side serially but should be overloaded if a block method is
available. ";

%feature("docstring")
NOX::Abstract::Group::applyJacobianTransposeMultiVector "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

applyJacobianTranspose for multiple right-hand sides

The default implementation here calls applyJacobianTranspose() for
each right hand side serially but should be overloaded if a block
method is available. ";

%feature("docstring")
NOX::Abstract::Group::applyJacobianInverseMultiVector "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

applyJacobianInverse for multiple right-hand sides

The default implementation here calls applyJacobianInverse() for each
right hand side serially but should be overloaded if a block solver is
available. ";

%feature("docstring")
NOX::Abstract::Group::applyRightPreconditioningMultiVector "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioningMultiVector(bool
useTranspose, Teuchos::ParameterList &params, const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

applyRightPreconditioning for multiple right-hand sides

The default implementation here calls applyRightPreconditioning() for
each right hand side serially but should be overloaded if a block
method is available. ";

/*  "Is" functions.  */

/* Checks to see if various objects have been computed. Returns true
if the corresponding \"compute\" function has been called since the
last change to the solution vector.

*/

%feature("docstring")  NOX::Abstract::Group::isF "virtual bool
NOX::Abstract::Group::isF() const =0

Return true if F is valid. ";

%feature("docstring")  NOX::Abstract::Group::isJacobian "bool
NOX::Abstract::Group::isJacobian() const

Return true if the Jacobian is valid.

Default implementation in NOX::Abstract::Group returns false. ";

%feature("docstring")  NOX::Abstract::Group::isGradient "bool
NOX::Abstract::Group::isGradient() const

Return true if the gradient is valid.

Default implementation in NOX::Abstract::Group returns false. ";

%feature("docstring")  NOX::Abstract::Group::isNewton "bool
NOX::Abstract::Group::isNewton() const

Return true if the Newton direction is valid.

Default implementation in NOX::Abstract::Group returns false. ";

/*  "Get" functions.  */

/* Note that these function do not check whether or not the vectors
are valid. Must use the \"Is\" functions for that purpose.

*/

%feature("docstring")  NOX::Abstract::Group::getX "virtual const
NOX::Abstract::Vector& NOX::Abstract::Group::getX() const =0

Return solution vector. ";

%feature("docstring")  NOX::Abstract::Group::getF "virtual const
NOX::Abstract::Vector& NOX::Abstract::Group::getF() const =0

Return F(x). ";

%feature("docstring")  NOX::Abstract::Group::getNormF "virtual double
NOX::Abstract::Group::getNormF() const =0

Return 2-norm of F(x).

In other words, \\\\[ \\\\sqrt{\\\\sum_{i=1}^n F_i^2} \\\\] ";

%feature("docstring")  NOX::Abstract::Group::getGradient "virtual
const NOX::Abstract::Vector& NOX::Abstract::Group::getGradient() const
=0

Return gradient. ";

%feature("docstring")  NOX::Abstract::Group::getNewton "virtual const
NOX::Abstract::Vector& NOX::Abstract::Group::getNewton() const =0

Return Newton direction. ";

%feature("docstring")  NOX::Abstract::Group::getXPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Abstract::Group::getXPtr() const

Return RCP to solution vector. ";

%feature("docstring")  NOX::Abstract::Group::getFPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Abstract::Group::getFPtr() const

Return RCP to F(x). ";

%feature("docstring")  NOX::Abstract::Group::getGradientPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Abstract::Group::getGradientPtr() const

Return RCP to gradient. ";

%feature("docstring")  NOX::Abstract::Group::getNewtonPtr "virtual
Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Abstract::Group::getNewtonPtr() const

Return RCP to Newton direction. ";

%feature("docstring")  NOX::Abstract::Group::clone "virtual
Teuchos::RCP<NOX::Abstract::Group>
NOX::Abstract::Group::clone(NOX::CopyType type=NOX::DeepCopy) const =0

Create a new Group of the same derived type as this one by cloning
this one, and return a ref count pointer to the new group.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\" (only the memory is allocated, the values are
not copied into the vectors and Jacobian). Returns NULL if clone is
not supported.

Any shared data should have its ownership transfered to this group
from the source for a NOX::DeepCopy. ";

%feature("docstring")  NOX::Abstract::Group::Group "NOX::Abstract::Group::Group()

Constructor.

Constructors for any derived object should always define a default
x-value so that getX() is always defined. ";

%feature("docstring")  NOX::Abstract::Group::~Group "virtual
NOX::Abstract::Group::~Group()

Destructor. ";

%feature("docstring")
NOX::Abstract::Group::getNormLastLinearSolveResidual "NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::getNormLastLinearSolveResidual(double &residual)
const

Return the norm of the last linear solve residual as the result of
either a call to computeNewton() or applyJacobianInverse().

NOX::Abstract::Group::NotDefined - Returned by default implementation
in NOX::Abstract::Group

NOX::Abstract::Group::BadDependency - If no linear solve has been
calculated

NOX::Abstract::Group::Failed - Any other type of failure

NOX::Abstract::Group::Ok - Otherwise ";


// File: classNOX_1_1Direction_1_1Utils_1_1InexactNewton.xml
%feature("docstring") NOX::Direction::Utils::InexactNewton "

Inexact Newton Utilities

If we use an iterative linear solver for a Newton-based solve, then
this is called an inexact Newton method. The tolerance used to
terminate the linear solve is called the forcing term. The forcing
term may be constant, or it may be adjustable. In either case, at
iteration $k$ we require, \\\\[ \\\\frac{\\\\|J_k d_k -
(-F_k)\\\\|}{\\\\|F_k\\\\|} \\\\leq \\\\eta_k. \\\\] Here $\\\\eta_k$
is the forcing term for iteration $k$.

This solution tolerance is to be enforced by the user's implementation
of NOX::Abstract::Group::computeNewton; it is passed in as the
\"Tolerance\" in the parameter list for that function.  Adjustable
forcing terms were introduced by Eisenstat and Walker (1982); here
they are implemented as described in Pernice and Walker (1998). We
have two choices for adjustable forcing terms:

Type 1

\\\\[ \\\\eta_k = \\\\left\\\\vert \\\\frac{\\\\| F_k \\\\| -
\\\\|J_{k-1} d_{k-1} - (-F_{k-1}) \\\\| } {\\\\|F_{k-1}\\\\|}
\\\\right\\\\vert \\\\]

With the following safeguards imposed: \\\\[
\\\\max\\\\{\\\\eta_{k-1}^{\\\\frac{1 + \\\\sqrt{5}}{2}},
\\\\eta_{\\\\min} \\\\} \\\\leq \\\\eta_k \\\\leq \\\\eta_{\\\\max}
\\\\]

Type 2

\\\\[ \\\\eta_k = \\\\gamma \\\\left(
\\\\frac{\\\\|F_k\\\\|}{\\\\|F_{k-1}\\\\|} \\\\right)^\\\\alpha \\\\]

With the following safeguards imposed: \\\\[ \\\\max\\\\{\\\\gamma
\\\\eta_{k-1}^{\\\\alpha}, \\\\eta_{\\\\min} \\\\} \\\\leq \\\\eta_k
\\\\leq \\\\eta_{\\\\max} \\\\]

Parameters

\"Forcing Term Method\" - Method to compute the forcing term, i.e.,
the tolerance for the linear solver. Choices are: \"Constant\"
[default]

\"Type 1\"

\"Type 2\"

\"Forcing Term Initial Tolerance\" - $\\\\eta_0$ (initial linear
solver tolerance). Defaults to 0.1.

\"Forcing Term Minimum Tolerance\" - $\\\\eta_{\\\\min}$. Defaults to
1.0e-6.

\"Forcing Term Maximum Tolerance\" - $\\\\eta_{\\\\max}$. Defaults to
0.01.

\"Forcing Term Alpha\" - $\\\\alpha$ (used only by \"Type 2\").
Defaults to 1.5.

\"Forcing Term Gamma\" - $\\\\gamma$ (used only by \"Type 2\").
Defaults to 0.9.

\"Forcing Term User Defined Norm\" (NOX::Parameter::UserNorm derived
object) - If using a Type 1 or Type 2 adjustable forcing term, the
norms used to calculate $ \\\\eta_k $ should be based on the same norm
used in the convergence test of the linear solver. Essentially this
means that the norm must account for LEFT scaling of any kind that is
applied to the linear system during a solve. If set, the computation
of $ \\\\eta_k $ will be done using a used defined function that is
passed in through a NOX::Parameter::Arbitrary derived object. It will
take the arbitrary object and cast it to a NOX::Parameter::UserNorm
object and use it for norm computations. If this parameter is not set,
this method uses the L-2 Norm for any norm computations of $ \\\\eta_k
$.

\"Set Tolerance in Parameter List\" When calling computeForcingTerm()
the value of the forcing term will be set in the parmaeter list
pointed to by the InexactNewton object. It will be set under
paramsPtr->sublist(<directionMethod>).sublist(\"Linear
Solver\").set(\"Tolerance\", eta_k). Defaults to true.

When using a forcing term, it's critically important the the residual
of the original system is used in the comparison. This can be an issue
if scaling or left preconditioning is applied to the linear system.
References

Michael Pernice and Homer F. Walker, NITSOL: A Newton Iterative Solver
for Nonlinear Systems, SISC 19(Jan 1998):302-318.

S. C. Eisenstat and H. F. Walker, Globally convergent inexact Newton
methods, SINUM 19(1982):400-408

C++ includes: NOX_Direction_Utils_InexactNewton.H ";

%feature("docstring")
NOX::Direction::Utils::InexactNewton::InexactNewton "NOX::Direction::Utils::InexactNewton::InexactNewton(const
Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList
&directionSublist)

Constructor. ";

%feature("docstring")
NOX::Direction::Utils::InexactNewton::~InexactNewton "NOX::Direction::Utils::InexactNewton::~InexactNewton()

Destructor. ";

%feature("docstring")  NOX::Direction::Utils::InexactNewton::reset "bool NOX::Direction::Utils::InexactNewton::reset(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &directionSublist)

Reset the utilities. ";

%feature("docstring")
NOX::Direction::Utils::InexactNewton::computeForcingTerm "double
NOX::Direction::Utils::InexactNewton::computeForcingTerm(const
NOX::Abstract::Group &soln, const NOX::Abstract::Group &oldSoln, int
niter, const NOX::Solver::Generic &solver, double eta_last=-1.0)

Called each iteration to reset the forcing term (ie, the convergence
tolerance for the linear solver).

if the user supplied eta_last then it will use this value for eta_km1
instead of looking for it in the \"Linear Solver\" sublist. ";


// File: classNOX_1_1Solver_1_1InexactTrustRegionBased.xml
%feature("docstring") NOX::Solver::InexactTrustRegionBased "

Newton-like solver using a trust region.

Our goal is to solve: $ F(x) = 0, $ where $ F:\\\\Re^n \\\\rightarrow
\\\\Re^n $. Alternatively, we might say that we wish to solve

$ \\\\min f(x) \\\\equiv \\\\frac{1}{2} \\\\|F(x)\\\\|^2_2. $

The trust region subproblem (TRSP) at iteration $k$ is given by

$ \\\\min \\\\; m_k(s) \\\\equiv f_k + g_k^T d + \\\\frac{1}{2} d^T
B_k d, \\\\mbox{ s.t. } \\\\|d\\\\| \\\\leq \\\\Delta_k \\\\quad
\\\\mbox{(TRSP)} $

where

$ f_k = f(x_k) = \\\\frac{1}{2} \\\\|F(x_k)\\\\|^2_2 $,

$ g_k = \\\\nabla f(x_k) = J(x_k)^T F(x_k) $,

$ B_k = J(x_k)^T J(x_k) \\\\approx \\\\nabla^2 f(x_k) $,

$ J(x_k)$ is the Jacobian of $F$ at $x_k$, and

$ \\\\Delta_k $ is the trust region radius.

The \"improvement ratio\" for a given step $ s $ is defined as

$ \\\\rho = \\\\displaystyle\\\\frac{ f(x_k) - f(x_k + d) } { m_k(0) -
m_k(d) } $

An iteration consists of the following steps.

Compute Newton-like direction: $n$

Compute Cauchy-like direction: $c$

If this is the first iteration, initialize $\\\\Delta$ as follows: If
$\\\\|n\\\\|_2 < \\\\Delta_{\\\\min}$, then $\\\\Delta = 2
\\\\Delta_{\\\\min}$; else, $\\\\Delta = \\\\|n\\\\|_2$.

Initialize $\\\\rho = -1$

While $\\\\rho < \\\\rho_{\\\\min}$ and $\\\\Delta >
\\\\Delta_{\\\\min}$, do the following.

Compute the direction $d$ as follows:

If $\\\\|n\\\\|_2 < \\\\Delta$, then take a Newton step by setting $d
= n$

Otherwise if $\\\\|c\\\\|_2 > \\\\Delta$, then take a Cauchy step by
setting $d = \\\\displaystyle\\\\frac{\\\\Delta}{\\\\|c\\\\|_2} c$

Otherwise, take a Dog Leg step by setting $ d = (1-\\\\gamma) c +
\\\\gamma n $ where $ \\\\gamma = \\\\displaystyle\\\\frac {-c^T a +
\\\\sqrt{ (c^Ta)^2 - (c^Tc - \\\\Delta^2) a^Ta}}{a^Ta} $ with $a =
n-c$.

Set $x_{\\\\rm new} = x + d$ and calculate $f_{\\\\rm new}$

If $f_{\\\\rm new} \\\\geq f$, then $\\\\rho = -1$ Otherwise $ \\\\rho
= \\\\displaystyle \\\\frac {f - f_{\\\\rm new}} {| d^T J F +
\\\\frac{1}{2} (J d)^T (J d)|} $

Update the solution: $x = x_{\\\\rm new}$

Update trust region:

If $\\\\rho < \\\\rho_{\\\\rm s}$ and $\\\\|n\\\\|_2 < \\\\Delta$,
then shrink the trust region to the size of the Newton step:
$\\\\Delta = \\\\|n\\\\|_2$.

Otherwise if $\\\\rho < \\\\rho_{\\\\rm s}$, then shrink the trust
region: $\\\\Delta = \\\\max \\\\{ \\\\beta_{\\\\rm s} \\\\Delta,
\\\\Delta_{\\\\min} \\\\} $.

Otherwise if $\\\\rho > \\\\rho_{\\\\rm e}$ and $\\\\|d\\\\|_2 =
\\\\Delta$, then expand the trust region: $\\\\Delta = \\\\min \\\\{
\\\\beta_{\\\\rm e} \\\\Delta, \\\\Delta_{\\\\rm max} \\\\} $.

Input Paramters

The following parameters should be specified in the \"Trust Region\"
sublist based to the solver.

\"Inner Iteration Method\" - Choice of trust region algorithm to use.
Choices are: \"Standard Trust Region\"

\"Inexact Trust Region\"

\"Direction\" - Sublist of the direction parameters for the Newton
point, passed to the NOX::Direction::Manager constructor. If this
sublist does not exist, it is created by default. Furthermore, if
\"Method\" is not specified in this sublist, it is added with a value
of \"Newton\".

\"Cauchy %Direction\" - Sublist of the direction parameters for the
Cauchy point, passed to the NOX::Direction::Manager constructor. If
this sublist does not exist, it is created by default. Furthermore, if
\"Method\" is not specified in this sublist, it is added with a value
of \"Steepest Descent\" Finally, if the sub-sublist \"Steepest
Descent\" does not exist, it is created and the parameter \"Scaling
Type\" is added and set to \"Quadratic Min Model\".

\"Minimum Trust Region Radius\" ( $\\\\Delta_{\\\\min}$) - Minimum
allowable trust region radius. Defaults to 1.0e-6.

\"Maximum Trust Region Radius\" ( $\\\\Delta_{\\\\max}$) - Minimum
allowable trust region radius. Defaults to 1.0e+10.

\"Minimum Improvement Ratio\" ( $\\\\rho_{\\\\min}$) - Minimum
improvement ratio to accept the step. Defaults to 1.0e-4.

\"Contraction Trigger Ratio\" ( $\\\\rho_{\\\\rm s}$) - If the
improvement ratio is less than this value, then the trust region is
contracted by the amount specified by the \"Contraction Factor\". Must
be larger than \"Minimum Improvement Ratio\". Defaults to 0.1.

\"Contraction Factor\" ( $\\\\beta_{\\\\rm s}$) - See above. Defaults
to 0.25.

\"Expansion Trigger Ratio\" ( $\\\\rho_{\\\\rm e}$) - If the
improvement ratio is greater than this value, then the trust region is
contracted by the amount specified by the \"Expansion Factor\".
Defaults to 0.75.

\"Expansion Factor\" ( $\\\\beta_{\\\\rm e}$) - See above. Defaults to
4.0.

\"Recovery Step\" - Defaults to 1.0.

\"Use Ared/Pred Ratio Calculation\" (boolean) - Defaults to false. If
set to true, this option replaces the algorithm used to compute the
improvement ratio, $ \\\\rho $, as described above. The improvement
ratio is replaced by an \"Ared/Pred\" sufficient decrease criteria
similar to that used in line search algorithms (see Eisenstat and
Walker, SIAM Journal on Optimization V4 no. 2 (1994) pp 393-422):
$\\\\rho = \\\\frac{\\\\|F(x) \\\\| - \\\\| F(x + d) \\\\| } {\\\\|
F(x) \\\\| - \\\\| F(x) + Jd \\\\| } $

\"Use Cauchy in Newton Direction\" - Boolean. Used only by the
\"Inexact Trust Region\" algorithm. If set to true, the initial guess
for the Newton direction computation will use the Cauchy direction as
the initial guess. Defaults to false.

\"Use Dogleg Segment Minimization\" - Boolean. Used only by the
\"Inexact Trust Region\" algorithm. If set to true, the $ \\\\tau $
parameter is minimized over the dogleg line segments instead of being
computed at the trust regioin radius. Used only by the \"Inexact Trust
Region\" algorithm. Defaults to false.

\"Use Counters\" - Boolean. If set to true, solver statistics will be
stored. Defaults to true.

\"Write Output Parameters\" - Boolean. If set to true, the solver
statistics will be written to the relevant \"Output\" sublists (see
Output Parameters). Defaults to true.

\"Solver Options\" - Sublist of general solver options. \"User Defined
Pre/Post Operator\" is supported. See NOX::Parameter::PrePostOperator
for more details.

Output Paramters

A sublist called \"Output\" will be created at the top level of the
parameter list and contain the following general solver parameters:

\"Nonlinear Iterations\" - Number of nonlinear iterations

\"2-Norm or Residual\" - Two-norm of final residual

A sublist called \"Output\" will be created in the \"Trust Region\"
sublist and contain the following trust region specific output
parameters:

\"Number of Cauchy Steps\" - Number of cauchy steps taken during the
solve.

\"Number of Newton Steps\" - Number of Newton steps taken during the
solve.

\"Number of Dogleg Steps\" - Number of Dogleg steps taken during the
solve.

\"Number of Trust Region Inner Iterations\" - Number of inner
iterations required to adjust the trust region radius.

\"Dogleg Steps: Average Fraction of Newton Step Length\" - Average
value of the fraction a dogleg step took compared to the full Newton
step. The fractional value is computed as $ \\\\mbox{frac} =
\\\\frac{\\\\| d \\\\|}{\\\\| n\\\\|} $.

\"Dogleg Steps: Average Fraction Between Cauchy and Newton Direction\"
- Average value of the fraction a dogleg step took between the Cauchy
and Newton directions. This is the $ \\\\gamma $ variable in the
standard dogleg algorithm and the $ \\\\tau $ parameter in the inexact
dogleg algorithm. A value of 0.0 is a full step in the Cauchy
direction and a value of 1.0 is a full step in the Newton direction.

Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)

C++ includes: NOX_Solver_InexactTrustRegionBased.H ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::InexactTrustRegionBased "NOX::Solver::InexactTrustRegionBased::InexactTrustRegionBased(const
Teuchos::RCP< NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)

Constructor.

See reset() for description. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::~InexactTrustRegionBased "NOX::Solver::InexactTrustRegionBased::~InexactTrustRegionBased()

Destructor. ";

%feature("docstring")  NOX::Solver::InexactTrustRegionBased::reset "void NOX::Solver::InexactTrustRegionBased::reset(const
NOX::Abstract::Vector &initialGuess, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests)

Resets the solver, sets a new status test, and sets a new initial
guess. ";

%feature("docstring")  NOX::Solver::InexactTrustRegionBased::reset "void NOX::Solver::InexactTrustRegionBased::reset(const
NOX::Abstract::Vector &initialGuess)

Resets the solver and sets a new initial guess. ";

%feature("docstring")  NOX::Solver::InexactTrustRegionBased::getStatus
"NOX::StatusTest::StatusType
NOX::Solver::InexactTrustRegionBased::getStatus()

Check current convergence and failure status. ";

%feature("docstring")  NOX::Solver::InexactTrustRegionBased::step "NOX::StatusTest::StatusType
NOX::Solver::InexactTrustRegionBased::step()

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")  NOX::Solver::InexactTrustRegionBased::solve "NOX::StatusTest::StatusType
NOX::Solver::InexactTrustRegionBased::solve()

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::getSolutionGroup "const
Abstract::Group &
NOX::Solver::InexactTrustRegionBased::getSolutionGroup() const

Return a reference to the current solution group. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::getPreviousSolutionGroup "const
Abstract::Group &
NOX::Solver::InexactTrustRegionBased::getPreviousSolutionGroup() const

Return a reference to the previous solution group. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::getNumIterations "int
NOX::Solver::InexactTrustRegionBased::getNumIterations() const

Get number of iterations. ";

%feature("docstring")  NOX::Solver::InexactTrustRegionBased::getList "const Teuchos::ParameterList &
NOX::Solver::InexactTrustRegionBased::getList() const

Return a reference to the solver parameters. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::getSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::InexactTrustRegionBased::getSolutionGroupPtr() const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::getPreviousSolutionGroupPtr "virtual Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::InexactTrustRegionBased::getPreviousSolutionGroupPtr()
const

Return a RCP to the previous solution group. ";

%feature("docstring")
NOX::Solver::InexactTrustRegionBased::getListPtr "virtual
Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::InexactTrustRegionBased::getListPtr() const

Return a RCP to the solver parameters. ";


// File: classNOX_1_1Multiphysics_1_1DataExchange_1_1Interface.xml
%feature("docstring") NOX::Multiphysics::DataExchange::Interface "

Provides a set of interfaces for users to provide information about
exchanging data between registered NOX solvers.

C++ includes: NOX_Multiphysics_DataExchange_Interface.H ";

%feature("docstring")
NOX::Multiphysics::DataExchange::Interface::Interface "NOX::Multiphysics::DataExchange::Interface::Interface()

Constructor. ";

%feature("docstring")
NOX::Multiphysics::DataExchange::Interface::~Interface "virtual
NOX::Multiphysics::DataExchange::Interface::~Interface()

Destructor. ";

%feature("docstring")
NOX::Multiphysics::DataExchange::Interface::exchangeAllData "virtual
bool NOX::Multiphysics::DataExchange::Interface::exchangeAllData()=0

Exchange data for all registered problems. ";

%feature("docstring")
NOX::Multiphysics::DataExchange::Interface::exchangeDataTo "virtual
bool NOX::Multiphysics::DataExchange::Interface::exchangeDataTo(int
solverId)=0

Exchange data for a specified problem - brings needed data from others
to this problem. ";


// File: classNOX_1_1Epetra_1_1Interface_1_1Jacobian.xml
%feature("docstring") NOX::Epetra::Interface::Jacobian "

Used by NOX::Epetra to provide a link to the external code for
Jacobian fills.

This is only required if the user wishes to supply their own Jacobian
operator.

C++ includes: NOX_Epetra_Interface_Jacobian.H ";

%feature("docstring")  NOX::Epetra::Interface::Jacobian::Jacobian "NOX::Epetra::Interface::Jacobian::Jacobian()

Constructor. ";

%feature("docstring")  NOX::Epetra::Interface::Jacobian::~Jacobian "virtual NOX::Epetra::Interface::Jacobian::~Jacobian()

Destructor. ";

%feature("docstring")
NOX::Epetra::Interface::Jacobian::computeJacobian "virtual bool
NOX::Epetra::Interface::Jacobian::computeJacobian(const Epetra_Vector
&x, Epetra_Operator &Jac)=0

Compute Jacobian given the specified input vector x. Returns true if
computation was successful. ";


// File: classNOX_1_1Epetra_1_1LinearSystem.xml
%feature("docstring") NOX::Epetra::LinearSystem "

Pure virtual class interface for allowing different linear solvers to
be used by the NOX::Epetra::Group.

C++ includes: NOX_Epetra_LinearSystem.H ";

%feature("docstring")  NOX::Epetra::LinearSystem::LinearSystem "NOX::Epetra::LinearSystem::LinearSystem()

Constructor. ";

%feature("docstring")  NOX::Epetra::LinearSystem::~LinearSystem "virtual NOX::Epetra::LinearSystem::~LinearSystem()

Destructor. ";

%feature("docstring")  NOX::Epetra::LinearSystem::applyJacobian "virtual bool NOX::Epetra::LinearSystem::applyJacobian(const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const =0

Applies Jacobian to the given input vector and puts the answer in the
result.

Computes \\\\[ v = J u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector. Returns true if
successful. ";

%feature("docstring")
NOX::Epetra::LinearSystem::applyJacobianTranspose "virtual bool
NOX::Epetra::LinearSystem::applyJacobianTranspose(const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const =0

Applies Jacobian-Transpose to the given input vector and puts the
answer in the result.

Computes \\\\[ v = J^T u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector. Returns true if
successful. ";

%feature("docstring")  NOX::Epetra::LinearSystem::applyJacobianInverse
"virtual bool
NOX::Epetra::LinearSystem::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector
&result)=0

Applies the inverse of the Jacobian matrix to the given input vector
and puts the answer in result.

Computes \\\\[ v = J^{-1} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The parameter list contains the linear solver options. ";

%feature("docstring")
NOX::Epetra::LinearSystem::applyRightPreconditioning "virtual bool
NOX::Epetra::LinearSystem::applyRightPreconditioning(bool
useTranspose, Teuchos::ParameterList &params, const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const =0

Apply right preconditiong to the given input vector.

Let $M$ be a right preconditioner for the Jacobian $J$; in other
words, $M$ is a matrix such that \\\\[ JM \\\\approx I. \\\\]

Compute \\\\[ u = M^{-1} v, \\\\] where $u$ is the input vector and
$v$ is the result vector.

If useTranspose is true, then the transpose of the preconditioner is
applied: \\\\[ u = {M^{-1}}^T v, \\\\] The transpose preconditioner is
currently only required for Tensor methods.

The parameter list contains the linear solver options. ";

%feature("docstring")  NOX::Epetra::LinearSystem::getScaling "virtual
Teuchos::RCP<NOX::Epetra::Scaling>
NOX::Epetra::LinearSystem::getScaling()=0

Get the scaling object. ";

%feature("docstring")  NOX::Epetra::LinearSystem::resetScaling "virtual void NOX::Epetra::LinearSystem::resetScaling(const
Teuchos::RCP< NOX::Epetra::Scaling > &s)=0

Sets the diagonal scaling vector(s) used in scaling the linear system.

See NOX::Epetra::Scaling for details on how to specify scaling of the
linear system. ";

%feature("docstring")  NOX::Epetra::LinearSystem::computeJacobian "virtual bool NOX::Epetra::LinearSystem::computeJacobian(const
NOX::Epetra::Vector &x)=0

Evaluates the Jacobian based on the solution vector x. ";

%feature("docstring")  NOX::Epetra::LinearSystem::createPreconditioner
"virtual bool NOX::Epetra::LinearSystem::createPreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &p, bool
recomputeGraph) const =0

Explicitly constructs a preconditioner based on the solution vector x
and the parameter list p.

The user has the option of recomputing the graph when a new
preconditioner is created. The NOX::Epetra::Group controls the isValid
flag for the preconditioner and will control when to call this. ";

%feature("docstring")
NOX::Epetra::LinearSystem::destroyPreconditioner "virtual bool
NOX::Epetra::LinearSystem::destroyPreconditioner() const =0

Deletes the preconditioner.

The NOX::Epetra::Group controls the isValid flag for the
preconditioner and will control when to call this. ";

%feature("docstring")
NOX::Epetra::LinearSystem::recomputePreconditioner "virtual bool
NOX::Epetra::LinearSystem::recomputePreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &linearSolverParams)
const =0

Recalculates the preconditioner using an already allocated graph.

Use this to compute a new preconditioner while using the same graph
for the preconditioner. This avoids deleting and reallocating the
memory required for the preconditioner and results in a big speed-up
for large-scale jobs. ";

%feature("docstring")
NOX::Epetra::LinearSystem::getPreconditionerPolicy "virtual
PreconditionerReusePolicyType
NOX::Epetra::LinearSystem::getPreconditionerPolicy(bool
advanceReuseCounter=true)=0

Evaluates the preconditioner policy at the current state.

NOTE: This can change values between nonlienar iterations. It is not a
static value. ";

%feature("docstring")
NOX::Epetra::LinearSystem::isPreconditionerConstructed "virtual bool
NOX::Epetra::LinearSystem::isPreconditionerConstructed() const =0

Indicates whether a preconditioner has been constructed. ";

%feature("docstring")  NOX::Epetra::LinearSystem::hasPreconditioner "virtual bool NOX::Epetra::LinearSystem::hasPreconditioner() const =0

Indicates whether the linear system has a preconditioner. ";

%feature("docstring")  NOX::Epetra::LinearSystem::getJacobianOperator
"virtual Teuchos::RCP<const Epetra_Operator>
NOX::Epetra::LinearSystem::getJacobianOperator() const =0

Return Jacobian operator. ";

%feature("docstring")  NOX::Epetra::LinearSystem::getJacobianOperator
"virtual Teuchos::RCP<Epetra_Operator>
NOX::Epetra::LinearSystem::getJacobianOperator()=0

Return Jacobian operator. ";

%feature("docstring")
NOX::Epetra::LinearSystem::getGeneratedPrecOperator "virtual
Teuchos::RCP<const Epetra_Operator>
NOX::Epetra::LinearSystem::getGeneratedPrecOperator() const =0

Return preconditioner operator.

Note: This should only be called if hasPreconditioner() returns true.
";

%feature("docstring")
NOX::Epetra::LinearSystem::getGeneratedPrecOperator "virtual
Teuchos::RCP<Epetra_Operator>
NOX::Epetra::LinearSystem::getGeneratedPrecOperator()=0

Return preconditioner operator. ";

%feature("docstring")
NOX::Epetra::LinearSystem::setJacobianOperatorForSolve "virtual void
NOX::Epetra::LinearSystem::setJacobianOperatorForSolve(const
Teuchos::RCP< const Epetra_Operator > &solveJacOp)=0

Set Jacobian operator for solve. ";

%feature("docstring")
NOX::Epetra::LinearSystem::setPrecOperatorForSolve "virtual void
NOX::Epetra::LinearSystem::setPrecOperatorForSolve(const Teuchos::RCP<
const Epetra_Operator > &solvePrecOp)=0

Set preconditioner operator for solve.

Note: This should only be called if hasPreconditioner() returns true.
";


// File: classNOX_1_1Epetra_1_1LinearSystemAztecOO.xml
%feature("docstring") NOX::Epetra::LinearSystemAztecOO "

Concrete implementation of NOX::Epetra::LinearSolver for AztecOO.

This solver provides the linear algebra services provided through the
AztecOO parallel iterative linear solver.

The NOX::Epetra::LinearSystemAztecOO object provides a flexible and
efficient way to interface an Epetra based application code to the
Aztec linear solver. This class handles construction of both the
preconditioners and AztecOO solver. All options are determined through
parameter lists and the basic constructors.

Constructing a Linear System

There are four different constructors that can be used. The difference
between constructors is based on whether the user supplies a Jacobian,
a preconditioner, neither or both.

If a Jacobian is not supplied then this object can create an
internally constructed Jacobian based on a Finite Difference or
Matrif-Free object. The user can specify which type of object to use
by setting the parameter \"Jacobian Operator\" in the parameter list.
The choices are \"Matrix-Free\" or \"Finite Difference\".

The user can supply their own preconditioner as an Epetra_Operator, or
they can supply their own matrix (an Epetra_RowMatrix derived object)
that can be used by one of the internal preconditioner libraries
(currently aztecoo or ifpack). If they supply their own preconditioner
the object must implement the Epetra_Operator::ApplyInverse method.
This is the method called during the linear solve to introduce
preconditoning into aztecoo. If the user supplies a matrix to be used
with an internal preconditioner, it must be derived from the
Epetra_RowMatrix class and must implement all functionality in the
Epetra_RowMatrix. If a Preconditioner is not supplied, then this
object can create an internal preconditioner matrix by finite
differencing or it can use the Jacobian operator if the Jacobian
derives from the Epetra_RowMatrix class. The user can specify which
type of object to use by setting the parameter \"Preconditioner
Operator\" in the parameter list. The choices are \"Use Jacobian\" or
\"Finite Difference\".

The Jacobian and preconditioner each require an interface to update
the state of the operator with respect to the solution vector and any
other parameters. There are three interfaces that can be implemented,
NOX::Epetra::Interface::Required, NOX::Epetra::Interface::Jacobian,
and NOX::Epetra::Interface::Preconditioner.

NOX::Epetra::Interface::Required supplies the computeF() function so
codes can tell NOX what the nonlinear equations are. This is the
minimum requirement to run nox through the epetra interface.
LinearSolverAztecOO requires this in some constructors so that if a
Jacobian or preconditoner is not supplied, it will use computeF from
the Required interface to estimate the Jacobian or preconditioner via
finite differences or directional derivatives.

NOX::Epetra::Interface::Jacobian is used for updating a user supplied
Jacobian opertor with respect to the solution vector and any other
parameters. It is required only in constructors in which a user
supplies a Jacobian operator.

NOX::Epetra::Interface::Preconditioner is used for updating a user
supplied preconditioner opertor/matrix with respect to the solution
vector and any other parameters. It is required only in constructors
in which a user supplies a preconditioner operator.

\"Linear Solver\" sublist parameters

A Teuchos::ParameterList called linearSolverParams is required in the
various constructors and during some method calls such as
applyJacobianInverse() and applyRightPreconditioning(). Typically,
this list is the \"Linear Solver\" sublist found in the nox parameter
list. The following parameters can be set in the linear solver sublist
and are vaild for the NOX::Epetra::LinearSolverAztecOO object: \"Aztec
Solver\" - Determine the iterative technique used in the solve. The
following options are valid:

\"GMRES\" - Restarted generalized minimal residual (default).

\"CG\" - Conjugate gradient.

\"CGS\" - Conjugate gradient squared.

\"TFQMR\" - Transpose-free quasi-minimal reasidual.

\"BiCGStab\" - Bi-conjugate gradient with stabilization.

\"LU\" - Sparse direct solve (single processor only).

\"Size of Krylov Subspace\" - When using restarted GMRES this sets the
maximum size of the Krylov subspace (defaults to 300).

\"Orthogonalization\" - The orthogonalization routine used for the
Gram-Schmidt orthogonalization procedure in Aztec. The following
options are valid:

\"Classical\" - (default).

\"Modified\"

\"Convergence Test\" - Algorithm used to calculate the residual that
is used for determining the convergence of the linear solver. See the
Aztec 2.1 manual for more information. The following options are
valid:

\"r0\" - (default)

\"rhs\"

\"norm\"

\"no scaling\"

\"sol\"

\"Tolerance\" - Tolerance used by AztecOO to determine if an iterative
linear solve has converged.

\"Ill-Conditioning Threshold\" - If the upper hessenberg matrix during
GMRES generates a condition number greater than this parameter value,
aztec will exit the linear solve returning the it's current solution.
The default is 1.0e11.

\"Preconditioner Iterations\" - Number of iterations an
AztecOO_Operator should take when solving the preconditioner. This is
only used if an AztecOO preconditioner is used and the solver makes a
call to NOX::Epetra::Group::applyRightPreconditioning(). This is NOT a
recomended approach.

\"Max Iterations\" - maximum number of iterations in the linear solve.
Default is 400.

\"Zero Initial Guess\" - boolean. Zero out the initial guess for
linear solves performed through applyJacobianInverse calls (i.e. zero
out the result vector before the linear solve). Defaults to false.

\"Throw Error on Prec Failure\" - boolean. If set to true, an
exception will be thrown if the preconditioner fails to initialize or
recompute/refactor. If set to false, a warning will br printed if the
NOX::Utils::Warning is enabled in the printing utilities (
NOX::Utils). Defaults to true.

\"Output Frequency\" - number of linear solve iterations between
output of the linear solve residual. Takes an integer, or one of the
AztecOO flags: AZ_none, AZ_last, or AZ_all as a value. Defaults to
AZ_last.

\"Jacobian Operator\" - When a constructor does not require a Jacobian
operator, the linear system will create a default operator using:

\"Matrix-Free\" (default)

\"Finite Difference\"

\"Preconditioner\" - Sets the choice of the preconditioner to use
during linear solves. The validity of the choice of preconditioner
will depend on the types of operators that are available for the
Jacobian and preconditioner. NOTE: This flag will override any
constructor details. For example, if you supply a preconditioner
operator in the constructor, it will not be used if this flag is set
to \"None\". If you supply an Epetra_Operator for the preconditioner
but the \"Preconditioner\" flag is set to \"AztecOO\" (this requires
an Epetra_RowMatrix for the preconditioner operator), this object will
exit with a failure. The valid options and any requirements on the
operator type are listed below:

\"None\" - No preconditioning. (default)

\"AztecOO\" - AztecOO internal preconditioner. This requires a
preconditioner operator that derives from the Epetra_RowMatrix class.

\"Ifpack\" - Ifpack internal preconditioner. This requires a
preconditioner object that derives from the Epetra_RowMatrix class or
it can use a Jacobian if the Jacobian derives from an
Epetra_RowMatrix. This option is deprecated. Please use \"New
Ifpack\".

\"New Ifpack\" - Ifpack internal preconditioner. This requires a
preconditioner object that derives from the Epetra_RowMatrix class or
it can use a Jacobian if the Jacobian derives from an
Epetra_RowMatrix.

\"User Defined\" - The user supplies an Epetra_Operator derived class.
Users must implement at a minimum the ApplyInverse() function of the
Epetra_Operator class since preconditioning of vectors is accomplished
through calls to this method.

\"Jacobian Operator\" - If a constructor is used that does not supply
a Jacobian operator, nox will create an internal Jacobian operator.
This flag is ONLY valid in such cases. This will determine which
Operator is used: \"Matrix-Free\" - Create a NOX::Epetra::MatrixFree
object.

\"Finite Difference\" - Create a NOX::Epetra::FiniteDifference object.

\"Preconditioner Operator\" - If a constructor is used that does not
supply a preconditioner operator, nox will create an internal
preconditioner operator. This flag is ONLY valid in such cases. This
will determine which Operator is used: \"Use Jacobian\" - Use the
Jacobian Operator (it must be an Epetra_RowMatrix derived object).

\"Finite Difference\" - Create a NOX::Epetra::FiniteDifference object.

\"Aztec Preconditioner\" - If the \"Preconditioner\" flag is set to
\"AztecOO\" then the specific AztecOO preconditioner is specified with
this flag. Currently supported preconditioners and their corresponding
parameters that can be set are shown below (See the Aztec 2.1 manual
for more information):

\"ilu\" - ilu preconditioning. This choice allows the following
additional parameters to be specified: \"Overlap\" - defaults to 0

\"Graph Fill\" - defaults to 0

\"ilut\" - ilut preconditioning. This choice allows the following
additional parameters to be specified: \"Overlap\" - defaults to 0

\"Fill Factor\" - defaults to 1.0

\"Drop Tolerance\" - defaults to 1.0e-12

\"Jacobi\" - k step Jacobi where k is set by the \"Steps\" flag:
\"Steps\" - defaults to 3.

\"Symmetric Gauss-Siedel\" - Non-overlapping domain decomposition k
step symmetric Gauss-Siedel where k is set by the \"Steps\" flag:

\"Steps\" - defaults to 3.

\"Polynomial\" - Neumann polynomial with order set by the parameter:
\"Polynomial Order\" - defaults to 3.

\"Least-squares Polynomial\" - Least-squares polynomial with order set
by the parameter: \"Polynomial Order\" - defaults to 3.

\"Ifpack\" - If the \"Preconditioner\" flag is set to \"New Ifpack\"
then any of the options supported by the Ifpack Create factory can be
specified using a Teuchos::ParameterList containing the Ifpack options
and then setting this as a parameter named \"Ifpack\" in the \"Linear
Solver\" sublist.

\"ML\" - If the \"Preconditioner\" flag is set to \"ML\" then any of
the options supported by the ML factory can be specified using a
Teuchos::ParameterList containing the ML options and then setting this
as a parameter named \"ML\" in the \"Linear Solver\" sublist.

\"Preconditioner Reuse Policy\" - (string) Allows the user to set how
and when the preconditioner should be computed. This flag supports
native Aztec, Ifpack and ML preconditioners. There are three options:
\"Rebuild\" - The \"Rebuild\" option always completely destroys and
then rebuilds the preconditioner each time a linear solve is
requested.

\"Reuse\" - The group/linear solver will not recompute the
preconditioner even if the group's solution vector changes. It just
blindly reuses what has been constructed. This turns off control of
preconditioner recalculation. This is a dangerous condition but can
really speed up the computations if the user knows what they are
doing. We don't recommend users trying this.

\"Recompute\" - Recomputes the preconditioner, but will try to
efficiently reuse any objects that don't need to be destroyed. How
efficient the \"Recompute\" option is depends on the type of
preconditioner. For example if we are using ILU from the Ifpack
library, we would like to not destroy and reallocate the graph each
solve. With this option, we tell Ifpack to reuse the graph from last
time - e.g the sparisty pattern has not changed between applications
of the preconditioner.

\"Max Age Of Prec\" - (int) If the \"Preconditioner Reuse Policy\" is
set to \"Reuse\", this integer tells the linear system how many times
to reuse the preconditioner before rebuilding it. Defaults to 1.

\"RCM Reordering\" - Enables RCM reordering in conjunction with domain
decomp incomplete factorization preconditioning. The following options
are valid:

\"Disabled\" - (default).

\"Enabled\"

\"Use Adaptive Linear Solve\" - Enables the use of AztecOO's
AdaptiveIterate() method instead of calling the Iterate() method. This
causes the preconditioning matrix to be modified to make the linear
solves easier. AztecOO will attempt to solve the linear system
multiple times now and if the solves are failing it will modify the
preconditioner and try again. Boolean value, defaults to false. NOTE:
This only works for internal Aztec preconditioners! The
\"Preconditioning\" parameter must be set to \"AztecOO: Jacobian
Matrix\" or \"AztecOO: User RowMatrix\". (NOTE: This parameter is
currently NOT supported)

\"Max Adaptive Solve Iterations\" - (integer) Maximum number of
attempts that the linear solver will make when trying to solve a
linear system. Defaults to 5. (NOTE: This parameter is currently NOT
supported)

\"Compute Scaling Manually\" - (boolean) The linear system can be
scaled if a NOX::Epetra::Scaling object is supplied to
LinearSystemAztecOO. When to compute the scaling can be handled either
manually by the user, or this object can automatically compute the
scaling prior to a linear solve. By setting this flag to true, the
user will call NOX::Epetra::Scaling::computeScaling() manually - on
their own! Setting this to false means the LinearSystemAztecOO object
will call the computeScaling function right before it applies the
scaling to the matrix in the applyJacobianInverse function. Default is
true (user will call compute scaling).

\"Output Solver Details\" - (boolean) Write the output sublist below
to the parameter list after each linear solve. default is true.

\"Write Linear System\" - (boolean) If set to true, the linear system
(Epetra_Map, Jacobian, LHS and RHS) is printed to a set of files in
matrix market format. This option requires building nox with the flag
--enable-nox-debug and building the EpetraExt library.

\"Write Linear System File Prefix\" - (string) If writing of the
linear system is enabled (see above parameter) users can change the
name of the output file prefix. The default is \"NOX_LinSys\". This
option requires building nox with the flag --enable-nox-debug and
building the EpetraExt library.

\"Output\" sublist

The parameter list passed in during calls to ApplyJacobianInverse()
will have an \"Output\" sublist created that contains the following
parameters if the flag \"Output Solver Details\" is set to true:

\"Acheived Tolerance\" - Actual tolerance achieved by the linear
solver computed via the convergence test requested.

\"Number of Linear Iterations\" - Number of iterations used by the
linear solver in the last call to applyJacobianInverse

\"Total Number of Linear Iterations\" - Total number of linear solve
iterations performed by groups that have used this input list

C++ includes: NOX_Epetra_LinearSystem_AztecOO.H ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO "NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO(Teuchos::ParameterList
&printingParams, Teuchos::ParameterList &linearSolverParams, const
Teuchos::RCP< NOX::Epetra::Interface::Required > &iReq, const
NOX::Epetra::Vector &cloneVector, const Teuchos::RCP<
NOX::Epetra::Scaling > scalingObject=Teuchos::null)

Constructor with no Operators.

Jacobian Operator will be constructed internally based on the
parameter \"Jacobian Operator\". Defaults to using a
NOX::Epetra::MatrixFree object. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO "NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO(Teuchos::ParameterList
&printingParams, Teuchos::ParameterList &linearSolverParams, const
Teuchos::RCP< NOX::Epetra::Interface::Required > &iReq, const
Teuchos::RCP< NOX::Epetra::Interface::Jacobian > &iJac, const
Teuchos::RCP< Epetra_Operator > &J, const NOX::Epetra::Vector
&cloneVector, const Teuchos::RCP< NOX::Epetra::Scaling >
scalingObject=Teuchos::null)

Constructor with a user supplied Jacobian Operator only.

Either there is no preconditioning or the preconditioner will be
used/created internally. The Jacobian (if derived from an
Epetra_RowMatrix class can be used with an internal preconditioner.
See the parameter key \"Preconditioner Operator\" for more details. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO "NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO(Teuchos::ParameterList
&printingParams, Teuchos::ParameterList &linearSolverParams, const
Teuchos::RCP< NOX::Epetra::Interface::Required > &i, const
Teuchos::RCP< NOX::Epetra::Interface::Preconditioner > &iPrec, const
Teuchos::RCP< Epetra_Operator > &M, const NOX::Epetra::Vector
&cloneVector, const Teuchos::RCP< NOX::Epetra::Scaling >
scalingObject=Teuchos::null)

Constructor with a user supplied Preconditioner Operator only.

Jacobian operator will be constructed internally based on the
parameter \"Jacobian Operator\" in the parameter list. See the
parameter key \"Jacobian Operator\" for more details. Defaults to
using a NOX::Epetra::MatrixFree object. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO "NOX::Epetra::LinearSystemAztecOO::LinearSystemAztecOO(Teuchos::ParameterList
&printingParams, Teuchos::ParameterList &linearSolverParams, const
Teuchos::RCP< NOX::Epetra::Interface::Jacobian > &iJac, const
Teuchos::RCP< Epetra_Operator > &J, const Teuchos::RCP<
NOX::Epetra::Interface::Preconditioner > &iPrec, const Teuchos::RCP<
Epetra_Operator > &M, const NOX::Epetra::Vector &cloneVector, const
Teuchos::RCP< NOX::Epetra::Scaling > scalingObject=Teuchos::null)

Constructor with user supplied separate objects for the Jacobian (J)
and Preconditioner (M). linearSolverParams is the \"Linear Solver\"
sublist of parameter list. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::~LinearSystemAztecOO "NOX::Epetra::LinearSystemAztecOO::~LinearSystemAztecOO()

Destructor. ";

%feature("docstring")  NOX::Epetra::LinearSystemAztecOO::applyJacobian
"bool NOX::Epetra::LinearSystemAztecOO::applyJacobian(const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const

Applies Jacobian to the given input vector and puts the answer in the
result.

Computes \\\\[ v = J u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector. Returns true if
successful. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::applyJacobianTranspose "bool
NOX::Epetra::LinearSystemAztecOO::applyJacobianTranspose(const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const

Applies Jacobian-Transpose to the given input vector and puts the
answer in the result.

Computes \\\\[ v = J^T u, \\\\] where $J$ is the Jacobian, $u$ is the
input vector, and $v$ is the result vector. Returns true if
successful. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::applyJacobianInverse "bool
NOX::Epetra::LinearSystemAztecOO::applyJacobianInverse(Teuchos::ParameterList
&linearSolverParams, const NOX::Epetra::Vector &input,
NOX::Epetra::Vector &result)

Applies the inverse of the Jacobian matrix to the given input vector
and puts the answer in result.

Computes \\\\[ v = J^{-1} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The parameter list contains the linear solver options. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::applyRightPreconditioning "bool
NOX::Epetra::LinearSystemAztecOO::applyRightPreconditioning(bool
useTranspose, Teuchos::ParameterList &linearSolverParams, const
NOX::Epetra::Vector &input, NOX::Epetra::Vector &result) const

Apply right preconditiong to the given input vector.

Let $M$ be a right preconditioner for the Jacobian $J$; in other
words, $M$ is a matrix such that \\\\[ JM \\\\approx I. \\\\]

Compute \\\\[ u = M^{-1} v, \\\\] where $u$ is the input vector and
$v$ is the result vector.

If useTranspose is true, then the transpose of the preconditioner is
applied: \\\\[ u = {M^{-1}}^T v, \\\\] The transpose preconditioner is
currently only required for Tensor methods.

The parameter list contains the linear solver options. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::createPreconditioner "bool
NOX::Epetra::LinearSystemAztecOO::createPreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &linearSolverParams,
bool recomputeGraph) const

Explicitly constructs a preconditioner based on the solution vector x
and the parameter list p.

The user has the option of recomputing the graph when a new
preconditioner is created. The NOX::Epetra::Group controls the isValid
flag for the preconditioner and will control when to call this. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::destroyPreconditioner "bool
NOX::Epetra::LinearSystemAztecOO::destroyPreconditioner() const

Deletes all objects associated with the chosen preconditioner. This is
called during linear solves and when the solution vector changes to
reset the preconditioner. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::recomputePreconditioner "bool
NOX::Epetra::LinearSystemAztecOO::recomputePreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &linearSolverParams)
const

Recalculates the preconditioner using an already allocated graph.

Use this to compute a new preconditioner while using the same graph
for the preconditioner. This avoids deleting and reallocating the
memory required for the preconditioner and results in a big speed-up
for large-scale jobs. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getPreconditionerPolicy "NOX::Epetra::LinearSystem::PreconditionerReusePolicyType
NOX::Epetra::LinearSystemAztecOO::getPreconditionerPolicy(bool
advanceReuseCounter=true)

Evaluates the preconditioner policy at the current state.

NOTE: This can change values between nonlienar iterations. It is not a
static value. ";

%feature("docstring")  NOX::Epetra::LinearSystemAztecOO::reset "void
NOX::Epetra::LinearSystemAztecOO::reset(Teuchos::ParameterList
&linearSolverParams)

Reset the linear solver parameters. ";

%feature("docstring")  NOX::Epetra::LinearSystemAztecOO::getScaling "Teuchos::RCP< NOX::Epetra::Scaling >
NOX::Epetra::LinearSystemAztecOO::getScaling()

Get the scaling object. ";

%feature("docstring")  NOX::Epetra::LinearSystemAztecOO::resetScaling
"void NOX::Epetra::LinearSystemAztecOO::resetScaling(const
Teuchos::RCP< NOX::Epetra::Scaling > &s)

Sets the diagonal scaling vector(s) used in scaling the linear system.
See NOX::Epetra::Scaling for details on how to specify scaling of the
linear system. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::computeJacobian "bool
NOX::Epetra::LinearSystemAztecOO::computeJacobian(const
NOX::Epetra::Vector &x)

Compute the Jacobian. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getJacobianInterface "Teuchos::RCP<
const NOX::Epetra::Interface::Jacobian >
NOX::Epetra::LinearSystemAztecOO::getJacobianInterface() const

NOX::Interface::Jacobian accessor. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getPrecInterface "Teuchos::RCP<
const NOX::Epetra::Interface::Preconditioner >
NOX::Epetra::LinearSystemAztecOO::getPrecInterface() const

NOX::Interface::Preconditioiner accessor. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::isPreconditionerConstructed "bool
NOX::Epetra::LinearSystemAztecOO::isPreconditionerConstructed() const

Indicates whether a preconditioner has been constructed. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::hasPreconditioner "bool
NOX::Epetra::LinearSystemAztecOO::hasPreconditioner() const

Indicates whether the linear system has a preconditioner. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getJacobianOperator "Teuchos::RCP<
const Epetra_Operator >
NOX::Epetra::LinearSystemAztecOO::getJacobianOperator() const

Jacobian Epetra_Operator accessor. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getJacobianOperator "Teuchos::RCP<
Epetra_Operator >
NOX::Epetra::LinearSystemAztecOO::getJacobianOperator()

Jacobian Epetra_Operator accessor. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getPrecOperator "Teuchos::RCP<
const Epetra_Operator >
NOX::Epetra::LinearSystemAztecOO::getPrecOperator() const

Preconditioner Epetra_Operator accessor (only the base matrix if using
an internal preconditioner - aztecoo or ifpack). ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getGeneratedPrecOperator "Teuchos::RCP< const Epetra_Operator >
NOX::Epetra::LinearSystemAztecOO::getGeneratedPrecOperator() const

Return preconditioner operator generated and stored in AztecOO.

Note: This should only be called if hasPreconditioner() returns true.
";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getGeneratedPrecOperator "Teuchos::RCP< Epetra_Operator >
NOX::Epetra::LinearSystemAztecOO::getGeneratedPrecOperator()

Return preconditioner operator generated and stored in AztecOO. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getTimeCreatePreconditioner "double
NOX::Epetra::LinearSystemAztecOO::getTimeCreatePreconditioner() const

Returns the total time (sec.) spent in createPreconditioner(). ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::getTimeApplyJacobianInverse "double
NOX::Epetra::LinearSystemAztecOO::getTimeApplyJacobianInverse() const

Returns the total time (sec.) spent in applyJacobianInverse(). ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::setJacobianOperatorForSolve "void
NOX::Epetra::LinearSystemAztecOO::setJacobianOperatorForSolve(const
Teuchos::RCP< const Epetra_Operator > &solveJacOp)

Set Jacobian operator for solve. ";

%feature("docstring")
NOX::Epetra::LinearSystemAztecOO::setPrecOperatorForSolve "void
NOX::Epetra::LinearSystemAztecOO::setPrecOperatorForSolve(const
Teuchos::RCP< const Epetra_Operator > &solvePrecOp)

Set preconditioner operator for solve.

Note: This should only be called if hasPreconditioner() returns true.
";


// File: classNOX_1_1Solver_1_1LineSearchBased.xml
%feature("docstring") NOX::Solver::LineSearchBased "

Nonlinear solver based on a line search (i.e., damping).

Solves $F(x)=0$ using an iterative line-search-based method.

Each iteration, the solver does the following.

Compute a search direction $d$ via a NOX::Direction method

Compute a step length $\\\\lambda$ and update $x$ as $x_{\\\\rm new} =
x_{\\\\rm old} + \\\\lambda d$ via a NOX::LineSearch method.

The iterations progress until the status tests (see NOX::StatusTest)
determine either failure or convergence.

To support several line searches and status tests, this version of the
solver has a getStepSize() function that returns $\\\\lambda$.  Input
Parameters

The following parameter list entries are valid for this solver:

\"Line Search\" - Sublist of the line search parameters, passed to the
NOX::LineSearch::Manager constructor. Defaults to an empty list.

\"Direction\" - Sublist of the direction parameters, passed to the
NOX::Direction::Manager constructor. Defaults to an empty list.

\"Solver Options\" - Sublist of general solver options. \"User Defined
Pre/Post Operator\" is supported. See NOX::Parameter::PrePostOperator
for more details.

Output Parameters

Every time solve() is called, a sublist for output parameters called
\"Output\" will be created and contain the following parameters.

\"Output\":

\"Nonlinear Iterations\" - Number of nonlinear iterations

\"2-Norm of Residual\" - Two-norm of final residual

Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)

C++ includes: NOX_Solver_LineSearchBased.H ";

%feature("docstring")  NOX::Solver::LineSearchBased::LineSearchBased "NOX::Solver::LineSearchBased::LineSearchBased(const Teuchos::RCP<
NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)

Constructor.

See reset(NOX::Abstract::Group&, NOX::StatusTest::Generic&,
Teuchos::ParameterList&) for description ";

%feature("docstring")  NOX::Solver::LineSearchBased::~LineSearchBased
"NOX::Solver::LineSearchBased::~LineSearchBased()

Destructor. ";

%feature("docstring")  NOX::Solver::LineSearchBased::reset "void
NOX::Solver::LineSearchBased::reset(const NOX::Abstract::Vector
&initialGuess, const Teuchos::RCP< NOX::StatusTest::Generic > &tests)

Resets the solver, sets a new status test, and sets a new initial
guess. ";

%feature("docstring")  NOX::Solver::LineSearchBased::reset "void
NOX::Solver::LineSearchBased::reset(const NOX::Abstract::Vector
&initialGuess)

Resets the solver and sets a new initial guess. ";

%feature("docstring")  NOX::Solver::LineSearchBased::getStatus "NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::getStatus()

Check current convergence and failure status. ";

%feature("docstring")  NOX::Solver::LineSearchBased::step "NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::step()

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")  NOX::Solver::LineSearchBased::solve "NOX::StatusTest::StatusType NOX::Solver::LineSearchBased::solve()

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")  NOX::Solver::LineSearchBased::getSolutionGroup
"const NOX::Abstract::Group &
NOX::Solver::LineSearchBased::getSolutionGroup() const

Return a reference to the current solution group. ";

%feature("docstring")
NOX::Solver::LineSearchBased::getPreviousSolutionGroup "const
NOX::Abstract::Group &
NOX::Solver::LineSearchBased::getPreviousSolutionGroup() const

Return a reference to the previous solution group. ";

%feature("docstring")  NOX::Solver::LineSearchBased::getNumIterations
"int NOX::Solver::LineSearchBased::getNumIterations() const

Get number of iterations. ";

%feature("docstring")  NOX::Solver::LineSearchBased::getList "const
Teuchos::ParameterList & NOX::Solver::LineSearchBased::getList() const

Return a reference to the solver parameters. ";

%feature("docstring")  NOX::Solver::LineSearchBased::getStepSize "double NOX::Solver::LineSearchBased::getStepSize() const ";

%feature("docstring")
NOX::Solver::LineSearchBased::getSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::LineSearchBased::getSolutionGroupPtr() const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Solver::LineSearchBased::getPreviousSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::LineSearchBased::getPreviousSolutionGroupPtr() const

Return a RCP to the previous solution group. ";

%feature("docstring")  NOX::Solver::LineSearchBased::getListPtr "virtual Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::LineSearchBased::getListPtr() const

Return a RCP to the solver parameters. ";


// File: classNOX_1_1Multiphysics_1_1Solver_1_1Manager.xml
%feature("docstring") NOX::Multiphysics::Solver::Manager "

Manager class to control the instantiation of the objects derived from
the NOX::Solver::Generic object.

Parameters

The following entries may be specified in the parameter list.

\"Nonlinear %Solver\" - Name of the solver method. Valid choices are
\"Line Search Based\" ( NOX::Solver::LineSearchBased) [Default]

\"Trust Region Based\" ( NOX::Solver::TrustRegionBased)

Deprecated The \"Nonlinear %Solver\" choices \"Newton\" and \"Line
Search\" are deprecated and revert to \"Line Search Based\". Likewise,
the choice \"Trust Region\" is deprecated and reverts to \"Trust
Region Based\".

Russell Hooper (SNL 1416)

C++ includes: NOX_Multiphysics_Solver_Manager.H ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::Manager "NOX::Multiphysics::Solver::Manager::Manager()

Empty constructor - reset called later to really construct it. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::Manager "NOX::Multiphysics::Solver::Manager::Manager(const Teuchos::RCP<
vector< Teuchos::RCP< NOX::Solver::Generic > > > &solvers, const
Teuchos::RCP< NOX::Multiphysics::DataExchange::Interface > &i, const
Teuchos::RCP< NOX::StatusTest::Generic > &t, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Constructor.

See reset() for a full description. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::Manager "NOX::Multiphysics::Solver::Manager::Manager(const Teuchos::RCP<
NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &t, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Constructor.

See reset() for a full description. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::~Manager "NOX::Multiphysics::Solver::Manager::~Manager()

Destructor. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::reset "bool NOX::Multiphysics::Solver::Manager::reset(const Teuchos::RCP<
vector< Teuchos::RCP< NOX::Solver::Generic > > > &solvers, const
Teuchos::RCP< NOX::Multiphysics::DataExchange::Interface > &i, const
Teuchos::RCP< NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params) ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::reset "void NOX::Multiphysics::Solver::Manager::reset(const
NOX::Abstract::Vector &initialGuess, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests) ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::reset "void NOX::Multiphysics::Solver::Manager::reset(const
NOX::Abstract::Vector &initialGuess)

Resets the solver and sets a new initial guess. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::getStatus "NOX::StatusTest::StatusType
NOX::Multiphysics::Solver::Manager::getStatus()

Check current convergence and failure status. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::step "NOX::StatusTest::StatusType NOX::Multiphysics::Solver::Manager::step()

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::solve "NOX::StatusTest::StatusType
NOX::Multiphysics::Solver::Manager::solve()

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")
NOX::Multiphysics::Solver::Manager::getSolutionGroup "const
NOX::Abstract::Group &
NOX::Multiphysics::Solver::Manager::getSolutionGroup() const

Return a reference to the current solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroup "const
NOX::Abstract::Group &
NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroup() const

Return a reference to the previous solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::Manager::getNumIterations "int
NOX::Multiphysics::Solver::Manager::getNumIterations() const

Get number of iterations. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::getList "const Teuchos::ParameterList &
NOX::Multiphysics::Solver::Manager::getList() const

Return a reference to the solver parameters. ";

%feature("docstring")
NOX::Multiphysics::Solver::Manager::getSolutionGroupPtr "Teuchos::RCP< const NOX::Abstract::Group >
NOX::Multiphysics::Solver::Manager::getSolutionGroupPtr() const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroupPtr "Teuchos::RCP< const NOX::Abstract::Group >
NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroupPtr()
const

Return a RCP to the previous solution group. ";

%feature("docstring")  NOX::Multiphysics::Solver::Manager::getListPtr
"Teuchos::RCP< const Teuchos::ParameterList >
NOX::Multiphysics::Solver::Manager::getListPtr() const

Return a RCP to the solver parameters. ";


// File: classNOX_1_1Epetra_1_1MatrixFree.xml
%feature("docstring") NOX::Epetra::MatrixFree "

Concrete implementation for creating an Epetra_Operator Jacobian based
on the Matrix-Free Newton-Krylov method.

Matrix-Free Newton-Krylov is a method that takes advantage of the fact
the Newton Krylov solvers do not require an explicit Jacobian matrix.
Newton-Krylov solvers only require the matrix-vector product $Jy$ in
the iteration sequence. This product can approximated by the
following:

\\\\[ Jy = \\\\frac{F(x + \\\\delta y) - F(x)}{\\\\delta} \\\\]

where $J$ is the Jacobian, $F$ is the function evaluation, $x$ is the
solution vector, $y$ is the vector to be operated on, and $\\\\delta$
is a scalar perturbation calculated by:

\\\\[ \\\\delta = \\\\lambda * (\\\\lambda + \\\\frac{\\\\|
x\\\\|}{\\\\| y\\\\|} ) \\\\]

where $ \\\\lambda = 1.0e-6 $.

C++ includes: NOX_Epetra_MatrixFree.H ";

%feature("docstring")  NOX::Epetra::MatrixFree::MatrixFree "MatrixFree::MatrixFree(Teuchos::ParameterList &printParams, const
Teuchos::RCP< NOX::Epetra::Interface::Required > &i, const
NOX::Epetra::Vector &cloneVector, bool useNewPerturbation=false)

Constructor.

The vector x is used to clone the solution vector. ";

%feature("docstring")  NOX::Epetra::MatrixFree::~MatrixFree "MatrixFree::~MatrixFree()

Pure virtual destructor. ";

%feature("docstring")  NOX::Epetra::MatrixFree::SetUseTranspose "int
MatrixFree::SetUseTranspose(bool UseTranspose)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

UseTranspose:  -If true, multiply by the transpose of operator,
otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  NOX::Epetra::MatrixFree::Apply "int
MatrixFree::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Y:  -A Epetra_MultiVector of dimension NumVectors containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  NOX::Epetra::MatrixFree::ApplyInverse "int
MatrixFree::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  - A Epetra_MultiVector of dimension NumVectors to solve for.

Y:  -A Epetra_MultiVector of dimension NumVectors containing result.

Integer error code, set to 0 if successful.

WARNING:  In order to work with AztecOO, any implementation of this
method must support the case where X and Y are the same object. ";

%feature("docstring")  NOX::Epetra::MatrixFree::NormInf "double
MatrixFree::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  NOX::Epetra::MatrixFree::Label "const char *
MatrixFree::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  NOX::Epetra::MatrixFree::UseTranspose "bool
MatrixFree::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  NOX::Epetra::MatrixFree::HasNormInf "bool
MatrixFree::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  NOX::Epetra::MatrixFree::Comm "const
Epetra_Comm & MatrixFree::Comm() const

Returns a reference to the Epetra_Comm communicator associated with
this operator. ";

%feature("docstring")  NOX::Epetra::MatrixFree::OperatorDomainMap "const Epetra_Map & MatrixFree::OperatorDomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  NOX::Epetra::MatrixFree::OperatorRangeMap "const Epetra_Map & MatrixFree::OperatorRangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  NOX::Epetra::MatrixFree::computeJacobian "bool
MatrixFree::computeJacobian(const Epetra_Vector &x, Epetra_Operator
&Jac)

Compute Jacobian given the specified input vector, x. Returns true if
computation was successful. ";

%feature("docstring")  NOX::Epetra::MatrixFree::setDifferenceMethod "void MatrixFree::setDifferenceMethod(DifferenceType type)

Set the type of perturbation method used (default is Forward). ";

%feature("docstring")  NOX::Epetra::MatrixFree::setLambda "void
MatrixFree::setLambda(double lambda_)

Allows the user to change the value of $ \\\\lambda $ in the
perturbation calculation. ";

%feature("docstring")  NOX::Epetra::MatrixFree::setComputePerturbation
"void MatrixFree::setComputePerturbation(bool bVal)

Flag that toggles whether MatrixFree should compute the perturbation
parameter $ \\\\eta $ or use a value supplied by the user through
setPerturbation(). ";

%feature("docstring")  NOX::Epetra::MatrixFree::setPerturbation "void
MatrixFree::setPerturbation(double eta_)

Set the perturbation parameter $ \\\\eta $. ";

%feature("docstring")  NOX::Epetra::MatrixFree::getPerturbation "double MatrixFree::getPerturbation() const

Returns the most recently used value of the perturbation parameter $
\\\\eta $. ";

%feature("docstring")  NOX::Epetra::MatrixFree::setGroupForComputeF "void MatrixFree::setGroupForComputeF(const NOX::Abstract::Group
&group)

Clone a NOX::Abstract::Group derived object and use the computeF()
method of that group for the perturbation instead of the
NOX::Epetra::Interface::Required::computeF() method. This is required
for LOCA to get the operators correct during homotopy. ";

%feature("docstring")
NOX::Epetra::MatrixFree::setSolverForComputeJacobian "void
MatrixFree::setSolverForComputeJacobian(const Teuchos::RCP<
NOX::Solver::Generic > &slvr)

Save a RCP to a solver, and use the Solver's current Group's
computeF() in the computeJacobian call, which can save a function call
by respecting the isValid flag. ";


// File: classNOX_1_1StatusTest_1_1MaxIters.xml
%feature("docstring") NOX::StatusTest::MaxIters "

Failure test based on the maximum number of nonlinear solver
iterations.

Let $k$ denote the current number of iterations (accessed via
NOX::Solver::getNumIterations) and $k_{\\\\max}$ denote the tolerance
set in the constructor of this status test. This test returns
NOX::StatusTest::Failed if $ k \\\\geq k_{\\\\rm max}. $ Otherwise, it
returns NOX::StatusTest::Unconverged.

If checkStatus is called with the type set to NOX::StatusTest::None,
it then the status is set to to NOX::Status::Unevaluated and returned.
(Also niters is set to -1.)

C++ includes: NOX_StatusTest_MaxIters.H ";

%feature("docstring")  NOX::StatusTest::MaxIters::MaxIters "NOX::StatusTest::MaxIters::MaxIters(int maxIterations, const
NOX::Utils *u=NULL)

Constructor. Specify the maximum number of nonlinear solver
iterations, $k_{\\\\max}$ ands optinally an error stream for printing
errors. ";

%feature("docstring")  NOX::StatusTest::MaxIters::~MaxIters "NOX::StatusTest::MaxIters::~MaxIters()

Destructor. ";

%feature("docstring")  NOX::StatusTest::MaxIters::checkStatus "NOX::StatusTest::StatusType
NOX::StatusTest::MaxIters::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
NOX::StatusType::None. If the test is skipped, then the status should
be set to NOX::StatusTest::Unevaluated. ";

%feature("docstring")  NOX::StatusTest::MaxIters::getStatus "NOX::StatusTest::StatusType NOX::StatusTest::MaxIters::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::MaxIters::print "ostream &
NOX::StatusTest::MaxIters::print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  NOX::StatusTest::MaxIters::getMaxIters "int
NOX::StatusTest::MaxIters::getMaxIters() const

Returns the Maximum number of iterations set in the constructor. ";

%feature("docstring")  NOX::StatusTest::MaxIters::getNumIters "int
NOX::StatusTest::MaxIters::getNumIters() const

Returns the current number of iterations taken by the solver.

Returns -1 if the status of this test is NOX::StatusTest::Unevaluated.
";


// File: classNOX_1_1LineSearch_1_1MoreThuente.xml
%feature("docstring") NOX::LineSearch::MoreThuente "

More'-Thuente Line Search. Original code by Dianne O'Leary, modfified
by Tammy Kolda and Roger Pawlowski for the NOX project. This version
has been slightly optimized and also supports Homer Walker's work on
adaptive forcing terms and Ared/Pred conditions. It also allows for
arbitrary merit functions and norms to be supplied by the user.

This code is based on the More'-Thuente line search from the 1983
MINPACK Project. More specifically, this code is based on Dianne
O'Leary's 1991 Matlab-implementation of the More'-Thuente line search.
The original comments are preserved in the descriptions of the
individual subroutines. What follows is an updated summary.

The merit function we are minimizing is given by

\\\\[ f(x) = 0.5 \\\\|F(x)\\\\|^2 \\\\]

(alternatively the user can define this)

The purpose of the More'-Thuente line search is to find a step which
satisfies a sufficient decrease condition and a curvature condition.
At each stage the subroutine updates an interval of uncertainty with
endpoints stx and sty. The interval of uncertainty is initially chosen
so that it contains a minimizer of the modified function

\\\\[ f(x+{\\\\rm stp} \\\\; s) - f(x) - {\\\\rm ftol} \\\\; {\\\\rm
stp} \\\\; (\\\\nabla f(x)^T s). \\\\]

If a step is obtained for which the modified function has a
nonpositive function value and nonnegative derivative, then the
interval of uncertainty is chosen so that it contains a minimizer of
$f(x+{\\\\rm stp}\\\\;s)$.

The algorithm is designed to find a step which satisfies one of two
sufficient decrease conditions:

(1) Armijo-Goldstein Condition \\\\[ f(x + {\\\\rm stp} \\\\; s)
\\\\le f(x) + {\\\\rm ftol} \\\\; {\\\\rm stp} \\\\; (\\\\nabla f(x)^T
s), \\\\]

or

(2) Ared/Pred Condtition \\\\[ F(x_{n-1}+ \\\\lambda s) \\\\le
F(x_{n-1})(1-\\\\alpha(1-\\\\eta)) \\\\]

and the curvature condition

\\\\[ \\\\vert \\\\nabla f(x + {\\\\rm stp} \\\\; s)^T s) \\\\vert
\\\\le {\\\\rm gtol} \\\\; |\\\\nabla f(x)^T s| \\\\]

If ftol is less than gtol and if, for example, the function is bounded
below, then there is always a step which satisfies both conditions. If
no step can be found which satisfies both conditions, then the
algorithm usually stops when rounding errors prevent further progress.
In this case stp only satisfies the sufficient decrease condition.

Modifications from NOX::LineSearch::MoreThuente  1. Added the option
to use Ared/Pred conditions as describe in Homer Walker's papers. 2.
Added support to use an adjustable forcing term as describe in Homer
Walker's papers. 3. Added the option to use directional derivatives in
computing the slope instead of explicitly computing the Jacobian. This
eliminates the need to recompute the Jacobian at each inner iteration
of the More'-Thuente algorithm. 4. Added the ability to use the
NOX::Parameter::UserNorm and NOX::Parameter::MeritFunction objects to
supply user defined norms and merit functions to the line search.

Implementation  This line search can be called via
NOX::LineSearch::Manager.

This line search is used if \"More'-Thuente2\" is the \"Method\" in
the \"Line Search\" sublist. (See NOX::LineSearch::Manager for
details.)

The following parameters can be specified for this line search in the
\"More'-Thuente2\" sublist of the \"Line Search\" sublist:

\"Sufficient Decrease Condition\" - Choice to use for the sufficient
decrease condition. Options are \"Ared/Pred\" or \"Armijo-Goldstein\"
(defaults to \"Armijo-Goldstein\").  1. \"Armijo-Goldstein\"
conditions: $ f(x_{n-1}+ \\\\lambda s) \\\\le f(x_{n-1}) +\\\\alpha
\\\\lambda f'(x_{n-1}) $  2. \"Ared/Pred\" conditions: $ \\\\|
F(x_{n-1}+ \\\\lambda s) \\\\| \\\\le \\\\| F(x_{n-1}) \\\\|
(1-\\\\alpha(1-\\\\eta)) $ where $ \\\\eta $ is the linear solve
tolerance in the inexact Newton method.

\"Sufficient Decrease\" - The ftol in the sufficient decrease
condition (defaults to 1.0e-4)

\"Curvature Condition\" - The gtol in the curvature condition
(defaults to 0.9999)

\"Optimize Slope Calculation\" - Boolean value. If set to true the
value of $ s^TJ^TF $ is estimated using a directional derivative in a
call to NOX::LineSearch::Common::computeSlopeWithOutJac. If false the
slope computation is computed with the
NOX::LineSearch::Common::computeSlope method. Setting this to true
eliminates having to compute the Jacobian at each inner iteration of
the More'-Thuente line search (defaults to false).

\"User Defined Norm\" - The user can redefine the norm that is used in
the Ared/Pred sufficient decrease condition by supplying a
NOX::Parameter::UserNorm derived object in the parameter list with
this key.

\"Merit Function\" - The user can supply their own merit function to
the line search by supplying a NOX::Parameter::MeritFunction derived
object with this key.

\"Interval Width\" - The maximum width of the interval containing the
minimum of the modified function (defaults to 1.0e-15)

\"Maximum Step\" - maximum allowable step length (defaults to 1.0e6)

\"Minimum Step\" - minimum allowable step length (defaults to 1.0e-12)

\"Max Iters\" - maximum number of right-hand-side and corresponding
Jacobian evaluations (defaults to 20)

\"Default Step\" - starting step length (defaults to 1.0)

\"Recovery Step Type\" - Determines the step size to take when the
line search fails. Choices are:

\"Constant\" [default] - Uses a constant value set in \"Recovery
Step\".

\"Last Computed Step\" - Uses the last value computed by the line
search algorithm.

\"Recovery Step\" - The value of the step to take when the line search
fails. Only used if the \"Recovery Step Type\" is set to \"Constant\".
Defaults to value for \"Default Step\".

Output Parameters  A sublist for output parameters will be created
called \"Output\" in the parameter list used to instantiate or reset
the class. Valid output parameters are:

\"Total Number of Line Search Calls\" - Total number of calls to the
compute() method of this line search.

\"Total Number of Non-trivial Line Searches\" - The total number of
steps that could not directly take a full step and meet the required
\"Convergence Criteria\" (i.e. The line search had to reduce the step
length using inner iteration calculations over iterate $ k $).

\"Total Number of Failed Line Searches\" - total number of line
searches that failed and used a recovery step.

\"Total Number of Line Search Inner Iterations\" - total number of
inner iterations $ k $ performed by this object.

C++ includes: NOX_LineSearch_MoreThuente.H ";

%feature("docstring")  NOX::LineSearch::MoreThuente::MoreThuente "NOX::LineSearch::MoreThuente::MoreThuente(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::LineSearch::MoreThuente::~MoreThuente "NOX::LineSearch::MoreThuente::~MoreThuente()

Destructor. ";

%feature("docstring")  NOX::LineSearch::MoreThuente::reset "bool
NOX::LineSearch::MoreThuente::reset(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params) ";

%feature("docstring")  NOX::LineSearch::MoreThuente::compute "bool
NOX::LineSearch::MoreThuente::compute(NOX::Abstract::Group &newgrp,
double &step, const NOX::Abstract::Vector &dir, const
NOX::Solver::Generic &s)

Perform a line search.

On input:

Parameters:
-----------

grp:  The initial solution vector, $x_{\\\\rm old}$.

dir:  A vector of directions to be used in the line search, $d$.

s:  The nonlinear solver.

On output:

Parameters:
-----------

step:  The distance the direction was scaled, $ \\\\lambda $.

grp:  The grp is updated with a new solution, $ x_{\\\\rm new} $,
resulting from the linesearch. Normally, for a single direction line
search, this is computed as:

\\\\[ x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d. \\\\]

Ideally, $ \\\\|F(x_{\\\\rm new})\\\\| < \\\\|F(x_{\\\\rm old})\\\\| $
(e.g the final direction is a descent direction).

Note that the dir object is a std::vector. For typical line searches
as described in the above equation, this vector is of size one. We
have used a std::vector to allow for special cases of multi-
directional line searches such as the Bader/Schnabel curvillinear line
search.

Return value is true for a successful line search computation. ";


// File: classNOX_1_1Epetra_1_1MultiVector.xml
%feature("docstring") NOX::Epetra::MultiVector "

Implementation of NOX::Abstract::MultiVector for Epetra multi-vectors.

C++ includes: NOX_Epetra_MultiVector.H ";

%feature("docstring")  NOX::Epetra::MultiVector::getEpetraMultiVector
"Epetra_MultiVector &
NOX::Epetra::MultiVector::getEpetraMultiVector()

Get reference to underlying Epetra vector. ";

%feature("docstring")  NOX::Epetra::MultiVector::getEpetraMultiVector
"const Epetra_MultiVector &
NOX::Epetra::MultiVector::getEpetraMultiVector() const

Get const reference to underlying Epetra vector. ";

%feature("docstring")  NOX::Epetra::MultiVector::init "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::init(double
value)

Initialize every element of this multi-vector with gamma. ";

%feature("docstring")  NOX::Epetra::MultiVector::random "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::random(bool
useSeed=false, int seed=1)

Initialize each element of this multi-vector with a random value. ";

%feature("docstring")  NOX::Epetra::MultiVector::setBlock "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::setBlock(const
NOX::Abstract::MultiVector &source, const vector< int > &index)

Copy the vectors in source to a set of vectors in *this. The
index.size() vectors in source are copied to a subset of vectors in
*this indicated by the indices given in index. ";

%feature("docstring")  NOX::Epetra::MultiVector::setBlock "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::setBlock(const
NOX::Epetra::MultiVector &source, const vector< int > &index) ";

%feature("docstring")  NOX::Epetra::MultiVector::augment "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::augment(const
NOX::Abstract::MultiVector &source)

Append the vectors in source to *this. ";

%feature("docstring")  NOX::Epetra::MultiVector::augment "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::augment(const
NOX::Epetra::MultiVector &source) ";

%feature("docstring")  NOX::Epetra::MultiVector::scale "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::scale(double
gamma)

Scale each element of this multivector by gamma. ";

%feature("docstring")  NOX::Epetra::MultiVector::update "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::update(double
alpha, const NOX::Abstract::MultiVector &a, double gamma=0.0)

Compute x = (alpha * a) + (gamma * x) where a is a multi-vector and x
= *this. ";

%feature("docstring")  NOX::Epetra::MultiVector::update "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::update(double
alpha, const NOX::Epetra::MultiVector &a, double gamma=0.0) ";

%feature("docstring")  NOX::Epetra::MultiVector::update "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::update(double
alpha, const NOX::Abstract::MultiVector &a, double beta, const
NOX::Abstract::MultiVector &b, double gamma=0.0)

Compute x = (alpha * a) + (beta * b) + (gamma * x) where a and b are
multi-vectors and x = *this. ";

%feature("docstring")  NOX::Epetra::MultiVector::update "NOX::Abstract::MultiVector & NOX::Epetra::MultiVector::update(double
alpha, const NOX::Epetra::MultiVector &a, double beta, const
NOX::Epetra::MultiVector &b, double gamma=0.0) ";

%feature("docstring")  NOX::Epetra::MultiVector::update "NOX::Abstract::MultiVector &
NOX::Epetra::MultiVector::update(Teuchos::ETransp transb, double
alpha, const NOX::Abstract::MultiVector &a, const
NOX::Abstract::MultiVector::DenseMatrix &b, double gamma=0.0)

Compute x = (alpha * a * b) + (gamma * x) where a is a multivector, b
is a dense matrix, x = *this, and op(b) = b if transb =
Teuchos::NO_TRANS and op(b) is b transpose if transb = Teuchos::TRANS.
";

%feature("docstring")  NOX::Epetra::MultiVector::update "NOX::Abstract::MultiVector &
NOX::Epetra::MultiVector::update(Teuchos::ETransp transb, double
alpha, const NOX::Epetra::MultiVector &a, const
NOX::Abstract::MultiVector::DenseMatrix &b, double gamma=0.0) ";

%feature("docstring")  NOX::Epetra::MultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector > NOX::Epetra::MultiVector::clone(CopyType
type=DeepCopy) const

Create a new Vector of the same underlying type by cloning \"this\",
and return a pointer to the new vector.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\". Note that there is no assumption that a vector
created by ShapeCopy is initialized to zeros.

Pointer to newly created vector or NULL if clone is not supported. ";

%feature("docstring")  NOX::Epetra::MultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector > NOX::Epetra::MultiVector::clone(int
numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")  NOX::Epetra::MultiVector::subCopy "Teuchos::RCP< NOX::Abstract::MultiVector >
NOX::Epetra::MultiVector::subCopy(const vector< int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")  NOX::Epetra::MultiVector::subView "Teuchos::RCP< NOX::Abstract::MultiVector >
NOX::Epetra::MultiVector::subView(const vector< int > &index) const

Creates a new multi-vector with ndex.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")  NOX::Epetra::MultiVector::norm "void
NOX::Epetra::MultiVector::norm(vector< double > &result,
NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm)
const

Norm. ";

%feature("docstring")  NOX::Epetra::MultiVector::multiply "void
NOX::Epetra::MultiVector::multiply(double alpha, const
NOX::Abstract::MultiVector &y, NOX::Abstract::MultiVector::DenseMatrix
&b) const

Computes the matrix-matrix product $\\\\alpha * y^T * (*this)$. ";

%feature("docstring")  NOX::Epetra::MultiVector::multiply "void
NOX::Epetra::MultiVector::multiply(double alpha, const
NOX::Epetra::MultiVector &y, NOX::Abstract::MultiVector::DenseMatrix
&b) const ";

%feature("docstring")  NOX::Epetra::MultiVector::MultiVector "NOX::Epetra::MultiVector::MultiVector(const Teuchos::RCP<
Epetra_MultiVector > &source, NOX::CopyType type=NOX::DeepCopy,
NOX::Epetra::MultiVector::MemoryType
memoryType=NOX::Epetra::MultiVector::CreateCopy)

Constructor that creates a COPY or VIEW of the Epetra_MultiVector.

NOTE: This ctor should just always create a view. It should be
implicit from the fact that a RCP object is being passed in that a
persisting relationship is present. However, since this could cause
confusion, the default is to make a copy and if a user wants a view,
they must pass in an explicit flag.

A VIEW of a vector uses the same underlying memory. WARNING: A View
can be dangerous since multiple objects can access the same memory
locations. ";

%feature("docstring")  NOX::Epetra::MultiVector::MultiVector "NOX::Epetra::MultiVector::MultiVector(const Epetra_MultiVector
&source, NOX::CopyType type=NOX::DeepCopy)

Construct by copying map and/or elements of an Epetra_MultiVector. ";

%feature("docstring")  NOX::Epetra::MultiVector::MultiVector "NOX::Epetra::MultiVector::MultiVector(const NOX::Epetra::MultiVector
&source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  NOX::Epetra::MultiVector::~MultiVector "NOX::Epetra::MultiVector::~MultiVector()

Destruct MultiVector. ";

%feature("docstring")  NOX::Epetra::MultiVector::length "int
NOX::Epetra::MultiVector::length() const

Return the length of multi-vector. ";

%feature("docstring")  NOX::Epetra::MultiVector::numVectors "int
NOX::Epetra::MultiVector::numVectors() const

Return the number of vectors in the multi-vector. ";

%feature("docstring")  NOX::Epetra::MultiVector::print "void
NOX::Epetra::MultiVector::print(std::ostream &stream) const

Print the vector. This is meant for debugging purposes only. ";


// File: classNOX_1_1MultiVector.xml
%feature("docstring") NOX::MultiVector "

Default implementation for NOX::Abstract::MultiVector using an array
of NOX::Abstract::MultiVector's.

C++ includes: NOX_MultiVector.H ";

%feature("docstring")  NOX::MultiVector::init "NOX::Abstract::MultiVector & NOX::MultiVector::init(double gamma)

Initialize every element of this multi-vector with gamma. ";

%feature("docstring")  NOX::MultiVector::random "NOX::Abstract::MultiVector & NOX::MultiVector::random(bool
useSeed=false, int seed=1)

Initialize each element of this multi-vector with a random value. ";

%feature("docstring")  NOX::MultiVector::setBlock "NOX::Abstract::MultiVector & NOX::MultiVector::setBlock(const
NOX::Abstract::MultiVector &source, const vector< int > &index)

Copy the vectors in source to a set of vectors in *this. The
index.size() vectors in source are copied to a subset of vectors in
*this indicated by the indices given in index. ";

%feature("docstring")  NOX::MultiVector::setBlock "NOX::Abstract::MultiVector & NOX::MultiVector::setBlock(const
NOX::MultiVector &source, const vector< int > &index) ";

%feature("docstring")  NOX::MultiVector::augment "NOX::Abstract::MultiVector & NOX::MultiVector::augment(const
NOX::Abstract::MultiVector &source)

Append the vectors in source to *this. ";

%feature("docstring")  NOX::MultiVector::augment "NOX::Abstract::MultiVector & NOX::MultiVector::augment(const
NOX::MultiVector &source) ";

%feature("docstring")  NOX::MultiVector::scale "NOX::Abstract::MultiVector & NOX::MultiVector::scale(double gamma)

Scale each element of this multivector by gamma. ";

%feature("docstring")  NOX::MultiVector::update "NOX::Abstract::MultiVector & NOX::MultiVector::update(double alpha,
const NOX::Abstract::MultiVector &a, double gamma=0.0)

Compute x = (alpha * a) + (gamma * x) where a is a multi-vector and x
= *this. ";

%feature("docstring")  NOX::MultiVector::update "NOX::Abstract::MultiVector & NOX::MultiVector::update(double alpha,
const NOX::MultiVector &a, double gamma=0.0) ";

%feature("docstring")  NOX::MultiVector::update "NOX::Abstract::MultiVector & NOX::MultiVector::update(double alpha,
const NOX::Abstract::MultiVector &a, double beta, const
NOX::Abstract::MultiVector &b, double gamma=0.0)

Compute x = (alpha * a) + (beta * b) + (gamma * x) where a and b are
multi-vectors and x = *this. ";

%feature("docstring")  NOX::MultiVector::update "NOX::Abstract::MultiVector & NOX::MultiVector::update(double alpha,
const NOX::MultiVector &a, double beta, const NOX::MultiVector &b,
double gamma=0.0) ";

%feature("docstring")  NOX::MultiVector::update "NOX::Abstract::MultiVector & NOX::MultiVector::update(Teuchos::ETransp
transb, double alpha, const NOX::Abstract::MultiVector &a, const
NOX::Abstract::MultiVector::DenseMatrix &b, double gamma=0.0)

Compute x = (alpha * a * b) + (gamma * x) where a is a multivector, b
is a dense matrix, x = *this, and op(b) = b if transb =
Teuchos::NO_TRANS and op(b) is b transpose if transb = Teuchos::TRANS.
";

%feature("docstring")  NOX::MultiVector::update "NOX::Abstract::MultiVector & NOX::MultiVector::update(Teuchos::ETransp
transb, double alpha, const NOX::MultiVector &a, const
NOX::Abstract::MultiVector::DenseMatrix &b, double gamma=0.0) ";

%feature("docstring")  NOX::MultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector > NOX::MultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new Vector of the same underlying type by cloning \"this\",
and return a pointer to the new vector.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\". Note that there is no assumption that a vector
created by ShapeCopy is initialized to zeros.

Pointer to newly created vector or NULL if clone is not supported. ";

%feature("docstring")  NOX::MultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector > NOX::MultiVector::clone(int numvecs)
const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")  NOX::MultiVector::subCopy "Teuchos::RCP<
NOX::Abstract::MultiVector > NOX::MultiVector::subCopy(const vector<
int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")  NOX::MultiVector::subView "Teuchos::RCP<
NOX::Abstract::MultiVector > NOX::MultiVector::subView(const vector<
int > &index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")  NOX::MultiVector::norm "void
NOX::MultiVector::norm(vector< double > &result,
NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm)
const

Norm. ";

%feature("docstring")  NOX::MultiVector::multiply "void
NOX::MultiVector::multiply(double alpha, const
NOX::Abstract::MultiVector &y, NOX::Abstract::MultiVector::DenseMatrix
&b) const

Computes the matrix-matrix product $\\\\alpha * y^T * (*this)$. ";

%feature("docstring")  NOX::MultiVector::multiply "void
NOX::MultiVector::multiply(double alpha, const NOX::MultiVector &y,
NOX::Abstract::MultiVector::DenseMatrix &b) const ";

%feature("docstring")  NOX::MultiVector::MultiVector "NOX::MultiVector::MultiVector(const NOX::Abstract::Vector &v, int
numVecs=1, NOX::CopyType type=NOX::DeepCopy)

Create MultiVector with numVecs columns out of a single
NOX::Abstract::Vector. ";

%feature("docstring")  NOX::MultiVector::MultiVector "NOX::MultiVector::MultiVector(const NOX::Abstract::Vector *const *vs,
int numVecs, NOX::CopyType type=NOX::DeepCopy)

Create MultiVector out of array of NOX::Abstract::Vector's. ";

%feature("docstring")  NOX::MultiVector::MultiVector "NOX::MultiVector::MultiVector(const MultiVector &source, NOX::CopyType
type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  NOX::MultiVector::~MultiVector "NOX::MultiVector::~MultiVector()

Destructor. ";

%feature("docstring")  NOX::MultiVector::length "int
NOX::MultiVector::length() const

Return the length of multi-vector. ";

%feature("docstring")  NOX::MultiVector::numVectors "int
NOX::MultiVector::numVectors() const

Return the number of vectors in the multi-vector. ";

%feature("docstring")  NOX::MultiVector::print "void
NOX::MultiVector::print(std::ostream &stream) const

Print the vector. This is meant for debugging purposes only. ";


// File: classNOX_1_1Abstract_1_1MultiVector.xml
%feature("docstring") NOX::Abstract::MultiVector "

Abstract interface for multi-vectors used by NOX.

C++ includes: NOX_Abstract_MultiVector.H ";

%feature("docstring")  NOX::Abstract::MultiVector::init "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::init(double
gamma)=0

Initialize every element of this multi-vector with gamma. ";

%feature("docstring")  NOX::Abstract::MultiVector::random "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::random(bool
useSeed=false, int seed=1)=0

Initialize each element of this multi-vector with a random value. ";

%feature("docstring")  NOX::Abstract::MultiVector::setBlock "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::setBlock(const
NOX::Abstract::MultiVector &source, const vector< int > &index)=0

Copy the vectors in source to a set of vectors in *this. The
index.size() vectors in source are copied to a subset of vectors in
*this indicated by the indices given in index. ";

%feature("docstring")  NOX::Abstract::MultiVector::augment "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::augment(const
NOX::Abstract::MultiVector &source)=0

Append the vectors in source to *this. ";

%feature("docstring")  NOX::Abstract::MultiVector::scale "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::scale(double
gamma)=0

Scale each element of this multivector by gamma. ";

%feature("docstring")  NOX::Abstract::MultiVector::update "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::update(double
alpha, const NOX::Abstract::MultiVector &a, double gamma=0.0)=0

Compute x = (alpha * a) + (gamma * x) where a is a multi-vector and x
= *this. ";

%feature("docstring")  NOX::Abstract::MultiVector::update "virtual
NOX::Abstract::MultiVector& NOX::Abstract::MultiVector::update(double
alpha, const NOX::Abstract::MultiVector &a, double beta, const
NOX::Abstract::MultiVector &b, double gamma=0.0)=0

Compute x = (alpha * a) + (beta * b) + (gamma * x) where a and b are
multi-vectors and x = *this. ";

%feature("docstring")  NOX::Abstract::MultiVector::update "virtual
NOX::Abstract::MultiVector&
NOX::Abstract::MultiVector::update(Teuchos::ETransp transb, double
alpha, const NOX::Abstract::MultiVector &a, const DenseMatrix &b,
double gamma=0.0)=0

Compute x = (alpha * a * op(b)) + (gamma * x) where a is a
multivector, b is a dense matrix, x = *this, and op(b) = b if transb =
Teuchos::NO_TRANS and op(b) is b transpose if transb = Teuchos::TRANS.
";

%feature("docstring")  NOX::Abstract::MultiVector::clone "virtual
Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Abstract::MultiVector::clone(NOX::CopyType type=NOX::DeepCopy)
const =0

Create a new Vector of the same underlying type by cloning \"this\",
and return a pointer to the new vector.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\". Note that there is no assumption that a vector
created by ShapeCopy is initialized to zeros.

Pointer to newly created vector or NULL if clone is not supported. ";

%feature("docstring")  NOX::Abstract::MultiVector::clone "virtual
Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Abstract::MultiVector::clone(int numvecs) const =0

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")  NOX::Abstract::MultiVector::subCopy "virtual
Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Abstract::MultiVector::subCopy(const vector< int > &index) const
=0

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")  NOX::Abstract::MultiVector::subView "virtual
Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Abstract::MultiVector::subView(const vector< int > &index) const
=0

Creates a new multi-vector with ndex.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")  NOX::Abstract::MultiVector::norm "virtual void
NOX::Abstract::MultiVector::norm(vector< double > &result,
NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm)
const =0

Norm. ";

%feature("docstring")  NOX::Abstract::MultiVector::multiply "virtual
void NOX::Abstract::MultiVector::multiply(double alpha, const
NOX::Abstract::MultiVector &y, DenseMatrix &b) const =0

Computes the matrix-matrix product $\\\\alpha * y^T * (*this)$. ";

%feature("docstring")  NOX::Abstract::MultiVector::MultiVector "NOX::Abstract::MultiVector::MultiVector()

Default constructor. Does nothing. ";

%feature("docstring")  NOX::Abstract::MultiVector::~MultiVector "virtual NOX::Abstract::MultiVector::~MultiVector()

Destructor. Does nothing. ";

%feature("docstring")  NOX::Abstract::MultiVector::length "virtual
int NOX::Abstract::MultiVector::length() const =0

Return the length of multi-vector. ";

%feature("docstring")  NOX::Abstract::MultiVector::numVectors "virtual int NOX::Abstract::MultiVector::numVectors() const =0

Return the number of vectors in the multi-vector. ";

%feature("docstring")  NOX::Abstract::MultiVector::print "virtual
void NOX::Abstract::MultiVector::print(std::ostream &stream) const =0

Print the vector. This is meant for debugging purposes only. ";


// File: classNOX_1_1Direction_1_1Newton.xml
%feature("docstring") NOX::Direction::Newton "

Newton direction computation

Computes the Newton direction by solving the Newton system. \\\\[ Jd =
-F \\\\]

Here $J$ is the n x n Jacobian matrix at the current iterate, $F$ is
the n-vector representing the nonlinear function at the current
iterate, and $d$ is the n-vector that we are solving for.

If we use an iterative linear solver for the Newton system, then this
is called an inexact Newton method. The tolerance used to terminate
the linear solve is called the forcing term. The forcing term may be
constant, or it may be adjustable. In either case, at iteration $k$ we
require, \\\\[ \\\\frac{\\\\|J_k d_k - (-F_k)\\\\|}{\\\\|F_k\\\\|}
\\\\leq \\\\eta_k. \\\\] Here $\\\\eta_k$ is the forcing term for
iteration $k$.

This solution tolerance is to be enforced by the user's implementation
of NOX::Abstract::Group::computeNewton; it is passed in as the
\"Tolerance\" in the parameter list for that function.  Adjustable
forcing terms were introduced by Eisenstat and Walker (1982); here
they are implemented as described in Pernice and Walker (1998). We
have two choices for adjustable forcing terms:

Type 1

\\\\[ \\\\eta_k = \\\\left\\\\vert \\\\frac{\\\\| F_k \\\\| -
\\\\|J_{k-1} d_{k-1} - (-F_{k-1}) \\\\| } {\\\\|F_{k-1}\\\\|}
\\\\right\\\\vert \\\\]

With the following safeguards imposed: \\\\[
\\\\max\\\\{\\\\eta_{k-1}^{\\\\frac{1 + \\\\sqrt{5}}{2}},
\\\\eta_{\\\\min} \\\\} \\\\leq \\\\eta_k \\\\leq \\\\eta_{\\\\max}
\\\\]

Type 2

\\\\[ \\\\eta_k = \\\\gamma \\\\left(
\\\\frac{\\\\|F_k\\\\|}{\\\\|F_{k-1}\\\\|} \\\\right)^\\\\alpha \\\\]

With the following safeguards imposed: \\\\[ \\\\max\\\\{\\\\gamma
\\\\eta_{k-1}^{\\\\alpha}, \\\\eta_{\\\\min} \\\\} \\\\leq \\\\eta_k
\\\\leq \\\\eta_{\\\\max} \\\\]

Parameters

\"Direction\": \"Method\" = \"Newton\" [required]

\"Direction\"/\"Newton\":

\"Forcing Term Method\" - Method to compute the forcing term, i.e.,
the tolerance for the linear solver. Choices are: \"Constant\"
[default]

\"Type 1\"

\"Type 2\"

\"Forcing Term Initial Tolerance\" - $\\\\eta_0$ (initial linear
solver tolerance). Defaults to 0.1.

\"Forcing Term Minimum Tolerance\" - $\\\\eta_{\\\\min}$. Defaults to
1.0e-6.

\"Forcing Term Maximum Tolerance\" - $\\\\eta_{\\\\max}$. Defaults to
0.01.

\"Forcing Term Alpha\" - $\\\\alpha$ (used only by \"Type 2\").
Defaults to 1.5.

\"Forcing Term Gamma\" - $\\\\gamma$ (used only by \"Type 2\").
Defaults to 0.9.

\"Rescue Bad %Newton Solve\" (Boolean) - If set to true, we will use
the computed direction even if the linear solve does not achieve the
tolerance specified by the forcing term. Defaults to true.

\"Direction\"/\"Newton\"/\"Linear Solver\":

\"Tolerance\" - Tolerance for the linear solve. This may be adjusted
automatically by the forcing calculation. Defaults to 1.0e-10. Will be
adjusted automatically by NOX if the \"Forcing Term Method\" is \"Type
1\" or \"Type 2\".

When using a forcing term, it's critically important the the residual
of the original system is used in the comparison. This can be an issue
if scaling or left preconditioning is applied to the linear system.
References

Michael Pernice and Homer F. Walker, NITSOL: A Newton Iterative Solver
for Nonlinear Systems, SISC 19(Jan 1998):302-318.

S. C. Eisenstat and H. F. Walker, Globally convergent inexact Newton
methods, SINUM 19(1982):400-408

C++ includes: NOX_Direction_Newton.H ";

%feature("docstring")  NOX::Direction::Newton::Newton "NOX::Direction::Newton::Newton(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::Direction::Newton::~Newton "NOX::Direction::Newton::~Newton()

Destructor. ";

%feature("docstring")  NOX::Direction::Newton::reset "bool
NOX::Direction::Newton::reset(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params)

Reset direction based on possibly new parameters. ";

%feature("docstring")  NOX::Direction::Newton::compute "bool
NOX::Direction::Newton::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::Generic &solver)

Compute the direction vector, dir, for a specific method given the
current group, grp.

The grp is not const so that we can compute the F vector, the Jacobian
matrix, the Newton vector, and so on.

Const access to the solver is used for getting additional information
such as the past solution, the iteration number, and so on. ";

%feature("docstring")  NOX::Direction::Newton::compute "bool
NOX::Direction::Newton::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::LineSearchBased &solver)

Same as compute( NOX::Abstract::Vector&, NOX::Abstract::Group&, const
NOX::Solver::Generic&).

Enables direct support for line search based solvers for the purpose
of efficiency since the LineSearchBased object has a getStep()
function that some directions require.

If it is not redefined in the derived class, it will just call the
compute with the NOX::Solver::Generic argument. ";


// File: classNOX_1_1LineSearch_1_1NonlinearCG.xml
%feature("docstring") NOX::LineSearch::NonlinearCG "

Use NonlinearCG linesearch.

This is a simple linesearch intended to be used with
NOX::Direction::NonlinearCG, which provides search direction $ d $, in
computing an update to the current solution vector $ x_{new} = x_{old}
+ \\\\lambda d $. It is designed to compute a step length $ \\\\lambda
$ consistent with the exact linesearch of Linear CG for linear
problems, and it avoids use of matrices by employing a directional
derivative (details below). The step length, $ \\\\lambda $ is
computed from a single evaluation of,

\\\\[ \\\\lambda = - \\\\frac{F(x_{old})^T d}{d^T J(x_{old})d} \\\\]

where $ J $ is the n x n Jacobian matrix. Explicit construction of $ J
$ is avoided by performing the product $ Jd $ using a directional
derivative (cf NOX::Epetra::MatrixFree):

\\\\[ J(x_{old})d \\\\approx \\\\frac{F(x_{old} + \\\\delta d) -
F(x_{old})}{\\\\delta} \\\\]

where $ \\\\delta = 10^{-6} (10^{-6} + ||x_{old}|| / ||d||) $ .

Derivation / Theory:

This linesearch is derived by attempting to achieve in a single step,
the following minimization:

\\\\[ \\\\min_\\\\lambda \\\\phi(\\\\lambda)\\\\equiv\\\\phi (x_{old}+
\\\\lambda d) \\\\]

where $ \\\\phi $ is a merit function chosen (but never explicitly
given) so that an equivalence to Linear CG holds, ie $
\\\\nabla\\\\phi(x) \\\\leftrightarrow F(x) $. The minimization above
can now be cast as an equation:

\\\\[ \\\\phi ' (\\\\lambda) = \\\\nabla\\\\phi (x_{old}+ \\\\lambda
d)^T d = F(x_{old}+ \\\\lambda d)^T d = 0~~. \\\\]

An approximate solution to this equation can be obtained from a
second-order expansion of \\\\[ \\\\phi(\\\\lambda) \\\\],

\\\\[ \\\\phi(\\\\lambda)\\\\approx\\\\phi (0) + \\\\phi '
(0)\\\\lambda + \\\\phi '' (0) \\\\frac{\\\\lambda^2}{2} \\\\]

from which it immediately follows

\\\\[ \\\\%lambda_{min} \\\\approx - \\\\frac{\\\\phi ' (0)}{\\\\phi
'' (0)} = - \\\\frac{F(x_{old})^T d}{d^T J(x_{old})d} \\\\]

Input Parameters

The NonlinearCG linesearch is selected using:

\"Line Search\"

\"Method\" = \"%NonlinearCG\" [required]

Currently, no adjustable parameters exist for this linesarch.

References

linesearch is adapted from ideas presented in Section 14.2 of:

Jonathan Richard Shewchuk,\"An Introduction to the Conjugate Gradient
Method Without the Agonizing Pain,\" 1994. Though presented within the
context of nonlinear optimization, the connection to solving nonlinear
equation systems is made via the equivalence $ f'(x)
\\\\leftrightarrow F(x) $.

Russ Hooper, Org. 9233, Sandia National Labs

C++ includes: NOX_LineSearch_NonlinearCG.H ";

%feature("docstring")  NOX::LineSearch::NonlinearCG::NonlinearCG "NOX::LineSearch::NonlinearCG::NonlinearCG(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::LineSearch::NonlinearCG::~NonlinearCG "NOX::LineSearch::NonlinearCG::~NonlinearCG()

Destructor. ";

%feature("docstring")  NOX::LineSearch::NonlinearCG::reset "bool
NOX::LineSearch::NonlinearCG::reset(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params) ";

%feature("docstring")  NOX::LineSearch::NonlinearCG::compute "bool
NOX::LineSearch::NonlinearCG::compute(NOX::Abstract::Group &newgrp,
double &step, const NOX::Abstract::Vector &dir, const
NOX::Solver::Generic &s)

Perform a line search.

On input:

Parameters:
-----------

grp:  The initial solution vector, $x_{\\\\rm old}$.

dir:  A vector of directions to be used in the line search, $d$.

s:  The nonlinear solver.

On output:

Parameters:
-----------

step:  The distance the direction was scaled, $ \\\\lambda $.

grp:  The grp is updated with a new solution, $ x_{\\\\rm new} $,
resulting from the linesearch. Normally, for a single direction line
search, this is computed as:

\\\\[ x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d. \\\\]

Ideally, $ \\\\|F(x_{\\\\rm new})\\\\| < \\\\|F(x_{\\\\rm old})\\\\| $
(e.g the final direction is a descent direction).

Note that the dir object is a std::vector. For typical line searches
as described in the above equation, this vector is of size one. We
have used a std::vector to allow for special cases of multi-
directional line searches such as the Bader/Schnabel curvillinear line
search.

Return value is true for a successful line search computation. ";


// File: classNOX_1_1Direction_1_1NonlinearCG.xml
%feature("docstring") NOX::Direction::NonlinearCG "

Calculates a search direction using the Nonlinear Conjugate Gradient
method.

Calculates the direction \\\\[ d = - M^{-1}(x) F(x) + \\\\beta
d_{prev} \\\\]

where $ M $ is a preconditioner and $ \\\\beta $ is an
orthogonalization parameter which can be computed in various ways (see
below), and $ d_{prev} $ is the search direction from the previous
nonlinear iteration.

This method provides a generalization of Linear CG to nonlinear
problems. It does this by computing a search direction using an
expression analogous to that of Linear CG. The negative of the current
residual vector, $ F(x) $ is taken, allowed to be preconditioned, and
then orthogonalized against the previous search direction. This
direction can sometimes be used successfully with the various choices
provided in NOX::Linesearch but is intended to be used with
NOX::Linesearch::NonlinearCG. In fact, the expected convergence
behavior of linear problems can only be achieved in this way.

To use this direction, specify that the \"Method\" is \"NonlinearCG\"
in the \"Direction\" sublist of the parameters that are passed to the
solver (see NOX::Direction::Manager for more information on choosing
the search direction).

The following options may be specified in the \"Nonlinear CG\" sublist
of the \"Direction\" sublist of the solver parameters.

\"Orthogonalize\" can be either of:

\"Fletcher-Reeves\" [default] - $ \\\\beta = \\\\frac{F(x)^T M^{-1}(x)
F(x)}{F(x_{prev})^T M^{-1}(x_{prev}) F(x_{prev})}$

\"Polak-Ribiere\" - $ \\\\beta = max \\\\left\\\\{ \\\\beta^{PR}, 0
\\\\right\\\\} $ ,

where $ \\\\beta^{PR} = \\\\frac{F(x)^T \\\\left[M^{-1}(x) F(x) -
M^{-1}(x_{prev}) F(x_{prev})\\\\right]}{F(x_{prev})^T M^{-1}(x_{prev})
F(x_{prev})}$

These comprise the two most popular choices for orthogonalization.
Both reduce to the linear result for linear problems. \"Polak-
Ribiere\" provides an implied restart capability by setting $ \\\\beta
= 0 $ anytime the computed value is less than zero.

\"Precondition\" can be either \"On\" or \"Off\" [default]: determines
whether or not to compute and apply preconditioner $ M $. If \"Off\"
is selected, no preconditioner is computed and the behavior is
equivalent to $ M = I $ where $ I $ is the identity matrix. If \"On\",
$ M $ is computed and applied as determined by the underlying
implementation of the \"applyRightPreconditioning\" method in the
Group.

\"Restart Frequency\" - An integer specification of the number of
nonlinear iterations between restarts [default = 10]. Restart
corresponds to setting $\\\\beta = 0$. A good heuristic is to limit
this value to the number of problem degrees of freedom. Setting this
value to 1 forces $ \\\\beta = 0 $ for every nonlinear iteration which
corresponds to suppressing orthogonalization against the previous
search direction.

References

information about both linear and nonlinear conjugate gradient methods
can be found in Chapter 5 of:

Nocedal & Wright, \"Numerical Optimization\", Springer-Verlag, New
York, 1999. Though presented within the context of nonlinear
optimization, the connection to nonlinear systems of equations is made
by the correspondence $ \\\\nabla f(x) \\\\leftrightarrow F(x) $ (cf
Algorithm 5.4).

Another useful useful reference is:

Jonathan Richard Shewchuk,\"An Introduction to the Conjugate Gradient
Method Without the Agonizing Pain,\" 1994. Chapter 14 provides a
summary of issues in generalizing linear CG to the nonlinear case.
Correspondence to NOX notation is made by the equivalence $ r
\\\\leftrightarrow f' \\\\leftrightarrow F(x) $ (cd Section 14.1).

C++ includes: NOX_Direction_NonlinearCG.H ";

%feature("docstring")  NOX::Direction::NonlinearCG::NonlinearCG "NonlinearCG::NonlinearCG(const Teuchos::RCP< NOX::GlobalData > &gd,
Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::Direction::NonlinearCG::~NonlinearCG "NonlinearCG::~NonlinearCG()

Destructor. ";

%feature("docstring")  NOX::Direction::NonlinearCG::reset "bool
NonlinearCG::reset(const Teuchos::RCP< NOX::GlobalData > &gd,
Teuchos::ParameterList &p)

derived ";

%feature("docstring")  NOX::Direction::NonlinearCG::compute "bool
NonlinearCG::compute(Abstract::Vector &dir, Abstract::Group &grp,
const Solver::Generic &solver)

derived ";

%feature("docstring")  NOX::Direction::NonlinearCG::compute "bool
NonlinearCG::compute(NOX::Abstract::Vector &dir, NOX::Abstract::Group
&grp, const NOX::Solver::LineSearchBased &solver)

Same as compute( NOX::Abstract::Vector&, NOX::Abstract::Group&, const
NOX::Solver::Generic&).

Enables direct support for line search based solvers for the purpose
of efficiency since the LineSearchBased object has a getStep()
function that some directions require.

If it is not redefined in the derived class, it will just call the
compute with the NOX::Solver::Generic argument. ";


// File: classNOX_1_1StatusTest_1_1NormF.xml
%feature("docstring") NOX::StatusTest::NormF "

Various convergence tests based on the norm of the residual.

Use the constructor to define the test based on the type of scaling
(see ScaleType) and the type of Tolerance (see Tolerance).

If checkStatus is called with the type set to NOX::StatusTest::None,
then the status is set to NOX::StatusTest::Unevaluated and returned.
(Also normF is set to 0.0.)

If checkStatus is called on a problem where the solution group does
not have F evaluated (i.e., problem.getSolutionGroup().isF() is
false), then the status is set to NOX::StatusTest::Unconverged and
returned. (Also normF is set to -1.0.)

Finally, we return NOX::StatusTest::Converged if $\\\\alpha <
\\\\beta$, and NOX::StatusTest::Unconverged otherwise. Here
$\\\\alpha$ represents the norm of $F(x)$ and $\\\\beta$ represents
the tolerance, as described below.

Let $\\\\gamma$ denote an optional scale factor defined as

$\\\\gamma = \\\\frac{1}{n}$ if sType in the constructor is
NOX::NormF::Scaled, and

Then $\\\\alpha$ is defined as follows:

If nType in the constructor is Abstract::Vector::TWO, then \\\\[
\\\\alpha = \\\\sqrt{ \\\\gamma \\\\sum_{i=1}^n F_i^2 } \\\\]

If nType in the constructor is Abstract::Vector::ONE, then \\\\[
\\\\alpha = \\\\gamma \\\\sum_{i=1}^n | F_i | \\\\]

If nType in the constructor is Abstract::Vector::INF, then \\\\[
\\\\alpha = \\\\gamma \\\\max_{i} | F_i | \\\\]

We set $\\\\beta$ as follows, based on the value of tolerance in the
constructor.

If an initial guess is provided, we use a relative tolerance defined
by \\\\[ \\\\beta = \\\\alpha_0 * \\\\mbox{tolerance} \\\\] Here
$\\\\alpha_0$ is the $\\\\alpha$ (as defined above) associated with
the initial guess.

Otherwise, we use an absolute tolerance defined by \\\\[ \\\\beta =
\\\\mbox{tolerance} \\\\]

C++ includes: NOX_StatusTest_NormF.H ";

%feature("docstring")  NOX::StatusTest::NormF::reset "void
NOX::StatusTest::NormF::reset(double tolerance)

Resets the user specified absolute or relative tolerance. ";

%feature("docstring")  NOX::StatusTest::NormF::reset "void
NOX::StatusTest::NormF::reset(NOX::Abstract::Group &initialGuess,
double tolerance)

Resets the user specified relative tolerance. ";

%feature("docstring")  NOX::StatusTest::NormF::getNormF "double
NOX::StatusTest::NormF::getNormF() const

Returns the value of the F-norm computed in the last call to
checkStatus. ";

%feature("docstring")  NOX::StatusTest::NormF::getTrueTolerance "double NOX::StatusTest::NormF::getTrueTolerance() const

Returns the true tolerance. ";

%feature("docstring")  NOX::StatusTest::NormF::getSpecifiedTolerance "double NOX::StatusTest::NormF::getSpecifiedTolerance() const

Returns the specified tolerance set in the constructor. ";

%feature("docstring")  NOX::StatusTest::NormF::getInitialTolerance "double NOX::StatusTest::NormF::getInitialTolerance() const

Returns the initial tolerance. ";

%feature("docstring")  NOX::StatusTest::NormF::NormF "NOX::StatusTest::NormF::NormF(double tolerance,
NOX::Abstract::Vector::NormType ntype, ScaleType stype=Scaled, const
NOX::Utils *u=NULL)

Constructor for absolute norm.

This constructor defaults to the Absolute tolerance type. ";

%feature("docstring")  NOX::StatusTest::NormF::NormF "NOX::StatusTest::NormF::NormF(double tolerance, ScaleType
stype=Scaled, const NOX::Utils *u=NULL)

Constructor for absolute norm.

This constructor defaults to the Absolute ToleranceType and TWO
NormType. ";

%feature("docstring")  NOX::StatusTest::NormF::NormF "NOX::StatusTest::NormF::NormF(NOX::Abstract::Group &initialGuess,
double tolerance, NOX::Abstract::Vector::NormType ntype, ScaleType
stype=Scaled, const NOX::Utils *u=NULL)

Constructor with initial guess (for relative norms).

This constructor defaults to the Relative tolerance type. ";

%feature("docstring")  NOX::StatusTest::NormF::NormF "NOX::StatusTest::NormF::NormF(NOX::Abstract::Group &initialGuess,
double tolerance, ScaleType stype=Scaled, const NOX::Utils *u=NULL)

Constructor with initial guess (for relative norms).

This constructor defaults to the Relative ToleranceType and TWO
NormType. ";

%feature("docstring")  NOX::StatusTest::NormF::~NormF "NOX::StatusTest::NormF::~NormF()

Destructor. ";

%feature("docstring")  NOX::StatusTest::NormF::checkStatus "NOX::StatusTest::StatusType NOX::StatusTest::NormF::checkStatus(const
NOX::Solver::Generic &problem, NOX::StatusTest::CheckType checkType)

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
NOX::StatusType::None. If the test is skipped, then the status should
be set to NOX::StatusTest::Unevaluated. ";

%feature("docstring")  NOX::StatusTest::NormF::getStatus "NOX::StatusTest::StatusType NOX::StatusTest::NormF::getStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::NormF::print "ostream &
NOX::StatusTest::NormF::print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classNOX_1_1StatusTest_1_1NormUpdate.xml
%feature("docstring") NOX::StatusTest::NormUpdate "

Various convergence tests based on the norm of the change in the
solution vector, $ x $, between outer iterations.

If checkStatusEfficiently is called with the type set to
NOX::StatusTest::None, then the status is set to
NOX::StatusTest::Unevaluated and returned. (Also normUpdate is set to
-1.0.)

If checkStatusEfficiently is called on the first iteration, then the
status is set to NOX::StatusTest::Unconverged and returned. (Also
normUpdate is set to -1.0.)

If checkStatusEfficiently is called on a problem where the solution
group does not have F evaluated (i.e.,
problem.getSolutionGroup().isF() is false), then the status is set to
NOX::StatusTest::Unconverged and returned. (Also normUpdate is set to
-1.0.)

Finally, we return NOX::StatusTest::Converged if $\\\\alpha <
\\\\beta$ and NOX::StatusTest::Uncoverged otherwise. Here $\\\\alpha$
represents the norm of $ \\\\Delta x $ and $\\\\beta$ represents the
tolerance. We define:

\\\\[ \\\\Delta x = x_k - x_{k-1} \\\\]

where $ x_k $ is the solution vector of the $ k $-th nonlinear
iterate.

Let $\\\\gamma$ denote an optional scale factor defined as

$\\\\gamma = \\\\frac{1}{n}$ if sType in the constructor is
NOX::NormF::Scaled, and

$\\\\gamma = 1$ if sType in the constructor is NOX::NormF::Unscaled.

Then $\\\\alpha$ is defined as follows:

If nType in the constructor is Abstract::Vector::TWO, then \\\\[
\\\\alpha = \\\\sqrt{ \\\\gamma \\\\sum_{i=1}^n \\\\Delta x_i^2 }
\\\\]

If nType in the constructor is Abstract::Vector::ONE, then \\\\[
\\\\alpha = \\\\gamma \\\\sum_{i=1}^n | \\\\Delta x_i | \\\\]

If nType in the constructor is Abstract::Vector::INF, then \\\\[
\\\\alpha = \\\\gamma \\\\max_{i} | \\\\Delta x_i | \\\\]

Finally, $\\\\beta$ is set to the tolerance in the constructor, i.e.,

\\\\[ \\\\beta = \\\\mbox{tolerance} \\\\]

C++ includes: NOX_StatusTest_NormUpdate.H ";

%feature("docstring")  NOX::StatusTest::NormUpdate::getNormUpdate "double NOX::StatusTest::NormUpdate::getNormUpdate() const

Returns the value of the Update-norm computed in the last call to
checkStatus. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::getTolerance "double NOX::StatusTest::NormUpdate::getTolerance() const

Returns the true tolerance. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::NormUpdate "NormUpdate::NormUpdate(double tolerance,
NOX::Abstract::Vector::NormType ntype, ScaleType stype=Scaled)

Constructor for absolute norm.

This constructor defaults to the Absolute tolerance type. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::NormUpdate "NormUpdate::NormUpdate(double tolerance, ScaleType stype=Scaled)

Constructor for absolute norm.

This constructor defaults to the Absolute ToleranceType and TWO
NormType. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::~NormUpdate "NormUpdate::~NormUpdate()

Destructor. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::checkStatus "StatusType NormUpdate::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
NOX::StatusType::None. If the test is skipped, then the status should
be set to NOX::StatusTest::Unevaluated. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::getStatus "StatusType NormUpdate::getStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::NormUpdate::print "ostream &
NormUpdate::print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classNOX_1_1StatusTest_1_1NormWRMS.xml
%feature("docstring") NOX::StatusTest::NormWRMS "

Convergence test based on the weighted root mean square norm fo the
solution update between iterations.

` If the number of iterations is zero, then the status is set to
NOX::StatusTest::Unconverged and returned. (Also, value is set to
1.0e+12.)

Otherwise, returns NOX::StatusTest::Converged if the three criteria
listed below are satisfied. Note that use of Criteria #2 and #3 depend
on the options set in the solver.

Weigthed root mean square norm is less than a specified tolerance:

\\\\[ ||\\\\delta x^k||_{wrms} < \\\\mbox{tolerance} \\\\]

where

\\\\[ ||\\\\delta x^k||_{wrms} \\\\equiv C \\\\sqrt{ \\\\frac{1}{N}
\\\\sum_{i=1}^N \\\\left( \\\\frac {(x^k_i-x^{k-1}_i)}{RTOL
|x^{k-1}_i| + ATOL_i} \\\\right) ^2 } \\\\]

Here:

$x_i^k$ denotes component $i$ of nonlinear iterate $k$.

$N$ denotes the number of unknowns

$RTOL$ denotes the relative error tolerance, specified via rtol in the
constructor

$ATOL$ denotes the absolution error tolerance, specified via atol in
the constructor. This can be a vector or a scalar.

$C$ denotes a weight, specified via the parameter BDFMultiplier in the
constructor.

If a line search based solver is used, the line search step size, $
\\\\lambda $, must be greater than a specified step size value, $
\\\\alpha $:

\\\\[ \\\\lambda > \\\\alpha \\\\]

The motivation for this test is to avoid detecting stagnation when in
fact the true problem is that the step size is just small.

The value of $\\\\alpha$ is set in the constructor via the argument
alpha. Setting $\\\\alpha$ to zero effectively eliminates this part of
the test.

The achieved linear solver tolerance, $ \\\\eta^k $ for nonlinear
iteration $ k $, must be less than a specified tolerance value, $
\\\\beta $; i.e.,

\\\\[ \\\\eta^k < \\\\beta \\\\]

The motivation for this test is to avoid detecting stagnation when in
fact the true problem is that the linear solve tolerance was not
accurate enough.

The value of $\\\\beta$ is set in the constructor via the argument
beta. Setting $\\\\beta$ to 1.0 effectively eliminates this part of
the test.

This criteria will only be used if the \"Achieved Tolerance\"
parameter (the value of $ \\\\eta^k $) is set by the linear solver in
the \"Newton\"/\"Linear Solver\"/\"Output\" sublist. The checkStatus()
method will search for this parameter.

References:

K. E. Brennam, S. L. Cambell, L. R. Petzold, Numerical Solution of
Initial-Value Problems in Differential-Algebraic Equations, Classics
in Applied Mathematics 14, SIAM 1996.

G. D. Byrne and A. C. Hindmarch, PVODE, an ODE Solver for Parallel
Computers, Technical Report UCRL-JC-132361, Rev. 1, Center for Applied
Scientific Computing (CASC), Lawrence Livermore National Lab, May
1999.

C++ includes: NOX_StatusTest_NormWRMS.H ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getNormWRMS "double
NormWRMS::getNormWRMS() const

Returns the value of WRMS norm. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getTolerance "double NormWRMS::getTolerance() const

Returns the requested tolerance set in the constructor. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getRTOL "double
NormWRMS::getRTOL() const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getATOL "double
NormWRMS::getATOL() const

Returns the absolute tolerance set in the constructor. If ATOL is a
vector, this will return a value of -1.0. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getBDFMultiplier "double NormWRMS::getBDFMultiplier() const

Returns the value of the BDFMultiplier set in the constructor. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getAlpha "double
NormWRMS::getAlpha() const

Returns the value of 'alpha' set in the constructor. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getBeta "double
NormWRMS::getBeta() const

Returns the value of 'beta' set in the constructor. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::NormWRMS "NormWRMS::NormWRMS(double rtol, double atol, double BDFMultiplier=1.0,
double tolerance=1.0, double alpha=1.0, double beta=0.5)

Constructor where ATOL is a scalar. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::NormWRMS "NormWRMS::NormWRMS(double rtol, const Teuchos::RCP< const
NOX::Abstract::Vector > &atol, double BDFMultiplier=1.0, double
tolerance=1.0, double alpha=1.0, double beta=0.5)

Constructor where ATOL is a vector. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::~NormWRMS "NormWRMS::~NormWRMS()

Destructor. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::checkStatus "StatusType NormWRMS::checkStatus(const NOX::Solver::Generic &problem,
NOX::StatusTest::CheckType checkType)

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
NOX::StatusType::None. If the test is skipped, then the status should
be set to NOX::StatusTest::Unevaluated. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::getStatus "StatusType NormWRMS::getStatus() const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::NormWRMS::print "ostream &
NormWRMS::print(ostream &stream, int indent=0) const

Output formatted description of stopping test to output stream. ";


// File: classNOX_1_1Epetra_1_1Observer.xml
%feature("docstring") NOX::Epetra::Observer "";

%feature("docstring")  NOX::Epetra::Observer::Observer "NOX::Epetra::Observer::Observer() ";

%feature("docstring")  NOX::Epetra::Observer::~Observer "virtual
NOX::Epetra::Observer::~Observer() ";

%feature("docstring")  NOX::Epetra::Observer::observeSolution "virtual void NOX::Epetra::Observer::observeSolution(const
Epetra_Vector &soln)=0

Method called by Piro NOXSolver as a hook for postprocessing. ";

%feature("docstring")  NOX::Epetra::Observer::observeSolution "virtual void NOX::Epetra::Observer::observeSolution(const
Epetra_Vector &soln, double time_or_param_val)

LOCA calls this version, which will discard param info in this default
implementation. ";


// File: classNOX_1_1LineSearch_1_1Polynomial.xml
%feature("docstring") NOX::LineSearch::Polynomial "

A polynomial line search, either quadratic or cubic.

This line search can be called via NOX::LineSearch::Manager.

The goal of the line search is to find a step $\\\\lambda$ for the
calculation $x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d$, given
$x_{\\\\rm old}$ and $ d $. To accomplish this goal, we minimize a
merit function $ \\\\phi(\\\\lambda) $ that measures the \"goodness\"
of the step $\\\\lambda$. The standard merit function is

\\\\[ \\\\phi(\\\\lambda) \\\\equiv \\\\frac{1}{2}||F (x_{\\\\rm old}
+ \\\\lambda s)||^2, \\\\]

but a user defined merit function can be used instead (see
computePhi() for details). Our first attempt is always to try a step $
\\\\lambda_0 $, and then check the stopping criteria. The value of $
\\\\lambda_0 $ is the specified by the \"Default Step\" parameter. If
the first try doesn't work, then we successively minimize polynomial
approximations, $ p_k(\\\\lambda) \\\\approx \\\\phi(\\\\lambda) $.

Stopping Criteria

The inner iterations continue until:

The sufficient decrease condition is met; see checkConvergence() for
details.

The maximum iterations are reached; see parameter \"Max Iters\". This
is considered a failure and the recovery step is used; see parameter
\"Recovery Step\".

The minimum step length is reached; see parameter \"Minimum Step\".
This is considered a line search failure and the recovery step is
used; see parameter \"Recovery Step\".

Polynomial Models of the Merit Function

We compute $ p_k(\\\\lambda) $ by interpolating $f$. For the quadratic
fit, we interpolate $ \\\\phi(0) $, $ \\\\phi'(0) $, and $
\\\\phi(\\\\lambda_{k-1}) $ where $ \\\\lambda_{k-1} $ is the $ k-1 $
approximation to the step. For the cubit fit, we additionally include
$\\\\phi(\\\\lambda{k-2})$.

The steps are calculated iteratively as follows, depending on the
choice of \"Interpolation Type\".

\"Quadratic\" - We construct a quadratic model of $\\\\phi$, and solve
for $\\\\lambda$ to get

\\\\[ \\\\lambda_{k} = \\\\frac{-\\\\phi'(0) \\\\lambda_{k-1}^2 }{2
\\\\left[ \\\\phi(\\\\lambda_{k-1}) - \\\\phi(0) -\\\\phi'(0)
\\\\lambda_{k-1} \\\\right]} \\\\]

\"Cubic\" - We construct a cubic model of $\\\\phi$, and solve for
$\\\\lambda$. We use the quadratic model to solve for $\\\\lambda_1$;
otherwise,

\\\\[ \\\\lambda_k = \\\\frac{-b+\\\\sqrt{b^2-3a \\\\phi'(0)}}{3a}
\\\\]

where

\\\\[ \\\\left[ \\\\begin{array}{c} a \\\\\\\\ \\\\\\\\ b
\\\\end{array} \\\\right] = \\\\frac{1}{\\\\lambda_{k-1} -
\\\\lambda_{k-2}} \\\\left[ \\\\begin{array}{cc} \\\\lambda_{k-1}^{-2}
& -\\\\lambda_{k-2}^{-2} \\\\\\\\ & \\\\\\\\
-\\\\lambda_{k-2}\\\\lambda_{k-1}^{-2} &
\\\\lambda_{k-1}\\\\lambda_{k-2}^{-2} \\\\end{array} \\\\right]
\\\\left[ \\\\begin{array}{c} \\\\phi(\\\\lambda_{k-1}) - \\\\phi(0) -
\\\\phi'(0)\\\\lambda_{k-1} \\\\\\\\ \\\\\\\\
\\\\phi(\\\\lambda_{k-2}) - \\\\phi(0) - \\\\phi'(0)\\\\lambda_{k-2}
\\\\end{array} \\\\right] \\\\]

\"Quadratic3\" - We construct a quadratic model of $\\\\phi$ using
$\\\\phi(0)$, $ \\\\phi(\\\\lambda_{k-1}) $ , and
$\\\\phi(\\\\lambda_{k-2})$. No derivative information for $\\\\phi$
is used. We let $\\\\lambda_1 = \\\\frac{1}{2} \\\\lambda_0$, and
otherwise

\\\\[ \\\\lambda_k = - \\\\frac{1}{2} \\\\frac{\\\\lambda_{k-1}^2
[\\\\phi(\\\\lambda_{k-2}) -\\\\phi(0)] - \\\\lambda_{k-2}^2
[\\\\phi(\\\\lambda_{k-1}) -\\\\phi(0)]} {\\\\lambda_{k-2}
[\\\\phi(\\\\lambda_{k-1}) -\\\\phi(0)] - \\\\lambda_{k-1}
[\\\\phi(\\\\lambda_{k-2}) -\\\\phi(0)]} \\\\]

For \"Quadratic\" and \"Cubic\", see computeSlope() for details on how
$ \\\\phi'(0) $ is calculated.

Bounds on the step length

We do not allow the step to grow or shrink too quickly by enforcing
the following bounds:

\\\\[ \\\\gamma_{min} \\\\; \\\\lambda_{k-1} \\\\leq \\\\lambda_k
\\\\le \\\\gamma_{max} \\\\; \\\\lambda_{k-1} \\\\]

Here $ \\\\gamma_{min} $ and $ \\\\gamma_{max} $ are defined by
parameters \"Min Bounds Factor\" and \"Max Bounds Factor\".

Input Parameters

\"Line Search\":

\"Method\" = \"Polynomial\" [required]

\"Line Search\"/\"Polynomial\":

\"Default Step\" - Starting step length, i.e., $\\\\lambda_0$.
Defaults to 1.0.

\"Max Iters\" - Maximum number of line search iterations. The search
fails if the number of iterations exceeds this value. Defaults to 100.

\"Minimum Step\" - Minimum acceptable step length. The search fails if
the computed $\\\\lambda_k$ is less than this value. Defaults to
1.0e-12.

\"Recovery Step Type\" - Determines the step size to take when the
line search fails. Choices are:

\"Constant\" [default] - Uses a constant value set in \"Recovery
Step\".

\"Last Computed Step\" - Uses the last value computed by the line
search algorithm.

\"Recovery Step\" - The value of the step to take when the line search
fails. Only used if the \"Recovery Step Type\" is set to \"Constant\".
Defaults to value for \"Default Step\".

\"Interpolation Type\" - Type of interpolation that should be used.
Choices are

\"Cubic\" [default]

\"Quadratic\"

\"Quadratic3\"

\"Min Bounds Factor\" - Choice for $ \\\\gamma_{min} $, i.e., the
factor that limits the minimum size of the new step based on the
previous step. Defaults to 0.1.

\"Max Bounds Factor\" - Choice for $ \\\\gamma_{max} $, i.e., the
factor that limits the maximum size of the new step based on the
previous step. Defaults to 0.5.

\"Sufficient Decrease Condition\" - See checkConvergence() for
details. Choices are:

\"Armijo-Goldstein\" [default]

\"Ared/Pred\"

\"None\"

\"Alpha Factor\" - Parameter choice for sufficient decrease condition.
See checkConvergence() for details. Defaults to 1.0e-4.

\"Force Interpolation\" (boolean) - Set to true if at least one
interpolation step should be used. The default is false which means
that the line search will stop if the default step length satisfies
the convergence criteria. Defaults to false.

\"Use Counters\" (boolean) - Set to true if we should use counters and
then output the result to the paramter list as described in Output
Parameters. Defaults to true.

\"Maximum Iteration for Increase\" - Maximum index of the nonlinear
iteration for which we allow a relative increase. See
checkConvergence() for further details. Defaults to 0 (zero).

\"Allowed Relative Increase\" - See checkConvergence() for details.
Defaults to 100.

\"User Defined Merit Function\" - The user can redefine the merit
function used; see computePhi() and NOX::Parameter::MeritFunction for
details.

\"User Defined Norm\" - The user can redefine the norm that is used in
the Ared/Pred sufficient decrease condition; see computeValue() and
NOX::Parameter::UserNorm for details.

Output Parameters

If the \"Use Counters\" parameter is set to true, then a sublist for
output parameters called \"Output\" will be created in the parameter
list used to instantiate or reset the class. Valid output parameters
are:

\"Total Number of Line Search Calls\" - Total number of calls to the
compute() method of this line search.

\"Total Number of Non-trivial Line Searches\" - Total number of steps
that could not directly take a full step and meet the required
\"Sufficient Decrease Condition\" (i.e., the line search had to do at
least one interpolation).

\"Total Number of Failed Line Searches\" - Total number of line
searches that failed and used a recovery step.

\"Total Number of Line Search Inner Iterations\" - Total number of
inner iterations of all calls to compute().

References

This line search is based on materials in the following:

Section 8.3.1 in C.T. Kelley, \"Iterative Methods for Linear and
Nonlinear Equations\", SIAM, 1995.

Section 6.3.2 and Algorithm 6.3.1 of J. E. Dennis Jr. and Robert B.
Schnabel, \"Numerical Methods for Unconstrained Optimization and
Nonlinear Equations,\" Prentice Hall, 1983.

Section 3.4 of Jorge Nocedal and Stephen J. Wright, \"Numerical
Optimization,\"Springer, 1999.

\"An Inexact Newton Method for Fully Coupled Solution of the Navier-
Stokes Equations with Heat and Mass Transfer\", Shadid, J. N.,
Tuminaro, R. S., and Walker, H. F., Journal of Computational Physics,
137, 155-185 (1997)

Russ Hooper, Roger Pawlowski, Tammy Kolda

C++ includes: NOX_LineSearch_Polynomial.H ";

%feature("docstring")  NOX::LineSearch::Polynomial::Polynomial "NOX::LineSearch::Polynomial::Polynomial(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  NOX::LineSearch::Polynomial::~Polynomial "NOX::LineSearch::Polynomial::~Polynomial()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Polynomial::reset "bool
NOX::LineSearch::Polynomial::reset(const Teuchos::RCP< NOX::GlobalData
> &gd, Teuchos::ParameterList &params) ";

%feature("docstring")  NOX::LineSearch::Polynomial::compute "bool
NOX::LineSearch::Polynomial::compute(NOX::Abstract::Group &newgrp,
double &step, const NOX::Abstract::Vector &dir, const
NOX::Solver::Generic &s)

Perform a line search.

On input:

Parameters:
-----------

grp:  The initial solution vector, $x_{\\\\rm old}$.

dir:  A vector of directions to be used in the line search, $d$.

s:  The nonlinear solver.

On output:

Parameters:
-----------

step:  The distance the direction was scaled, $ \\\\lambda $.

grp:  The grp is updated with a new solution, $ x_{\\\\rm new} $,
resulting from the linesearch. Normally, for a single direction line
search, this is computed as:

\\\\[ x_{\\\\rm new} = x_{\\\\rm old} + \\\\lambda d. \\\\]

Ideally, $ \\\\|F(x_{\\\\rm new})\\\\| < \\\\|F(x_{\\\\rm old})\\\\| $
(e.g the final direction is a descent direction).

Note that the dir object is a std::vector. For typical line searches
as described in the above equation, this vector is of size one. We
have used a std::vector to allow for special cases of multi-
directional line searches such as the Bader/Schnabel curvillinear line
search.

Return value is true for a successful line search computation. ";


// File: classNOX_1_1Epetra_1_1Interface_1_1Preconditioner.xml
%feature("docstring") NOX::Epetra::Interface::Preconditioner "

Used by NOX::Epetra to provide a link to the external code for
Precondtioner fills.

This is only required if the user wishes to supply their own
preconditioner operator.

C++ includes: NOX_Epetra_Interface_Preconditioner.H ";

%feature("docstring")
NOX::Epetra::Interface::Preconditioner::Preconditioner "NOX::Epetra::Interface::Preconditioner::Preconditioner()

Constructor. ";

%feature("docstring")
NOX::Epetra::Interface::Preconditioner::~Preconditioner "virtual
NOX::Epetra::Interface::Preconditioner::~Preconditioner()

Destructor. ";

%feature("docstring")
NOX::Epetra::Interface::Preconditioner::computePreconditioner "virtual bool
NOX::Epetra::Interface::Preconditioner::computePreconditioner(const
Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList
*precParams=0)=0

Computes a user defined preconditioner. ";


// File: classNOX_1_1Solver_1_1PrePostOperator.xml
%feature("docstring") NOX::Solver::PrePostOperator "

Functor to process the pre/post operator object in the parameter list.

This is a wrapper class for a user derived
NOX::Abstract::PrePostOperator (ppo) object. All solvers use this
class so we don't have to repeat all parsing code in each NOX::Solver.
This class searches the \"Solver Options\" parameter list passed into
the ctor and if a ppo is found will wrap the object.

For instructions on how to implement a PrePostOperator, see
NOX::Abstract::PrePostOperator.

C++ includes: NOX_Solver_PrePostOperator.H ";

%feature("docstring")  NOX::Solver::PrePostOperator::PrePostOperator "NOX::Solver::PrePostOperator::PrePostOperator(const Teuchos::RCP<
NOX::Utils > &utils, Teuchos::ParameterList &solverOptionsSubList)

Ctor. ";

%feature("docstring")  NOX::Solver::PrePostOperator::~PrePostOperator
"NOX::Solver::PrePostOperator::~PrePostOperator()

Destructor. ";

%feature("docstring")  NOX::Solver::PrePostOperator::reset "void
NOX::Solver::PrePostOperator::reset(const Teuchos::RCP< NOX::Utils >
&utils, Teuchos::ParameterList &solverOptionsSublist)

Resets the pre/post operator. ";

%feature("docstring")  NOX::Solver::PrePostOperator::runPreIterate "void NOX::Solver::PrePostOperator::runPreIterate(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the start of a call to
NOX::Solver::Generic::iterate(). ";

%feature("docstring")  NOX::Solver::PrePostOperator::runPostIterate "void NOX::Solver::PrePostOperator::runPostIterate(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the end of a call to
NOX::Solver::Generic::iterate(). ";

%feature("docstring")  NOX::Solver::PrePostOperator::runPreSolve "void NOX::Solver::PrePostOperator::runPreSolve(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the start of a call to
NOX::Solver::Generic::solve(). ";

%feature("docstring")  NOX::Solver::PrePostOperator::runPostSolve "void NOX::Solver::PrePostOperator::runPostSolve(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the end of a call to
NOX::Solver::Generic::solve(). ";


// File: classNOX_1_1Abstract_1_1PrePostOperator.xml
%feature("docstring") NOX::Abstract::PrePostOperator "

NOX's pure virtual class to allow users to insert pre and post
operations into nox's solvers (before and after the
NOX::Solver::Generic::iterate() and NOX::Solver::Generic::solve()
methods).

The user should implement their own concrete implementation of this
class and register it as a
Teuchos::RCP<NOX::Abstract::PrePostoperator> in the \"Solver Options\"
sublist.

To create and use a user defined pre/post operators:

Create a pre/post operator that derives from
NOX::Abstract::PrePostOperator. For example, the pre/post operator Foo
might be defined as shown below.

Create the appropriate entries in the parameter list, as follows.

Roger Pawlowski (SNL 9233)

C++ includes: NOX_Abstract_PrePostOperator.H ";

%feature("docstring")  NOX::Abstract::PrePostOperator::PrePostOperator
"NOX::Abstract::PrePostOperator::PrePostOperator()

Abstract Vector constructor (does nothing) ";

%feature("docstring")  NOX::Abstract::PrePostOperator::PrePostOperator
"NOX::Abstract::PrePostOperator::PrePostOperator(const
NOX::Abstract::PrePostOperator &source)

Copy constructor (doesnothing). ";

%feature("docstring")
NOX::Abstract::PrePostOperator::~PrePostOperator "virtual
NOX::Abstract::PrePostOperator::~PrePostOperator()

Abstract Vector destructor (does nothing) ";

%feature("docstring")  NOX::Abstract::PrePostOperator::runPreIterate "void NOX::Abstract::PrePostOperator::runPreIterate(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the start of a call to
NOX::Solver::Generic::iterate(). ";

%feature("docstring")  NOX::Abstract::PrePostOperator::runPostIterate
"void NOX::Abstract::PrePostOperator::runPostIterate(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the end of a call to
NOX::Solver::Generic::iterate(). ";

%feature("docstring")  NOX::Abstract::PrePostOperator::runPreSolve "void NOX::Abstract::PrePostOperator::runPreSolve(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the start of a call to
NOX::Solver::Generic::solve(). ";

%feature("docstring")  NOX::Abstract::PrePostOperator::runPostSolve "void NOX::Abstract::PrePostOperator::runPostSolve(const
NOX::Solver::Generic &solver)

User defined method that will be executed at the end of a call to
NOX::Solver::Generic::solve(). ";


// File: classNOX_1_1LineSearch_1_1Utils_1_1Printing.xml
%feature("docstring") NOX::LineSearch::Utils::Printing "

Common line search utilites for printing line search information to
the screen.

All line searches should print output results in a similar format.
This utility provides common output routines.

C++ includes: NOX_LineSearch_Utils_Printing.H ";

%feature("docstring")  NOX::LineSearch::Utils::Printing::Printing "NOX::LineSearch::Utils::Printing::Printing(const Teuchos::RCP<
NOX::Utils > &u)

Default constructor. ";

%feature("docstring")  NOX::LineSearch::Utils::Printing::~Printing "NOX::LineSearch::Utils::Printing::~Printing()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Utils::Printing::reset "void
NOX::LineSearch::Utils::Printing::reset(const Teuchos::RCP< NOX::Utils
> &u) ";

%feature("docstring")
NOX::LineSearch::Utils::Printing::printOpeningRemarks "void
NOX::LineSearch::Utils::Printing::printOpeningRemarks(const string
&lineSearchName) const

Prints the opening information. ";

%feature("docstring")  NOX::LineSearch::Utils::Printing::printStep "void NOX::LineSearch::Utils::Printing::printStep(int n, double step,
double oldf, double newf, const string s=\"\", bool unscaleF=true)
const

Print out step information for the inner iterations of a line search
algorithm.

Example of output from the inner iterations of a Polynomial line
search:

************************************************************************
-- Polynomial Line Search --    1: step = 1.000e+00 oldf = 2.403e+00
newf = 1.076e+03   2: step = 1.000e-01 oldf = 2.403e+00 newf =
4.440e+00   3: step = 1.000e-02 oldf = 2.403e+00 newf = 2.394e+00
(STEP ACCEPTED!)
************************************************************************

Parameters:
-----------

unscaleF:  - If this is true (the default), than the values printed
are $ \\\\sqrt{2 * {\\\\rm oldf}} $ and $ \\\\sqrt{2 * {\\\\rm newf}}
$. This is to accomodate the standard merit function, $ \\\\phi(x) =
\\\\frac{1}{2} \\\\|F(x)\\\\|^2 $. ";


// File: classNOX_1_1Random.xml
%feature("docstring") NOX::Random "

A class to compute uniformly distributed random numbers in (-1,1).

The Random class computes pseudo-random (double precision) numbers
uniformly distributed between -1.0 and 1.0 using a multiplicative
congruential generator with modulus 2^31-1 (a Lehmer generator). For a
numerical and mathematical treatment of the algorithm, see \"Random
number generators:  good ones are hard to find\" by Stephen K. Park
and Keith W. Miller, Communications of the ACM, Vol. 31 No. 10 (1988).

C++ includes: NOX_Random.H ";

%feature("docstring")  NOX::Random::Random "NOX::Random::Random()

Initialize random number generator with a random seed.

The random seed is computed using the POSIX rand() function. ";

%feature("docstring")  NOX::Random::Random "NOX::Random::Random(int
s)

Initialize random number generator with the given seed.

The seed should be an integer between 1 and 2147483646 = 2^32-2
(inclusive). If the supplied seed is invalid, an error message is
printed and the seed is replaced by 1. ";


// File: classNOX_1_1Epetra_1_1Interface_1_1Required.xml
%feature("docstring") NOX::Epetra::Interface::Required "

Supplies NOX with the set nonlinear equations.

This is the minimum required information to solve a nonlinear problem
using the NOX::Epetra objects for the linear algebra implementation.
Used by NOX::Epetra::Group to provide a link to the external code for
residual fills.

C++ includes: NOX_Epetra_Interface_Required.H ";

%feature("docstring")  NOX::Epetra::Interface::Required::Required "NOX::Epetra::Interface::Required::Required()

Constructor. ";

%feature("docstring")  NOX::Epetra::Interface::Required::~Required "virtual NOX::Epetra::Interface::Required::~Required()

Destructor. ";

%feature("docstring")  NOX::Epetra::Interface::Required::computeF "virtual bool NOX::Epetra::Interface::Required::computeF(const
Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)=0

Compute the function, F, given the specified input vector x. Returns
true if computation was successful. ";


// File: classNOX_1_1Epetra_1_1Scaling.xml
%feature("docstring") NOX::Epetra::Scaling "

Object to control scaling of vectors and linear systems.

Currently this assumes a diagonal scaling only! Once epetra can do
matrix-matrix multiplies we will expand this class.

C++ includes: NOX_Epetra_Scaling.H ";

%feature("docstring")  NOX::Epetra::Scaling::Scaling "NOX::Epetra::Scaling::Scaling()

Constructor. ";

%feature("docstring")  NOX::Epetra::Scaling::~Scaling "NOX::Epetra::Scaling::~Scaling()

Virtual destructor. ";

%feature("docstring")  NOX::Epetra::Scaling::addUserScaling "void
NOX::Epetra::Scaling::addUserScaling(ScaleType type, const
Teuchos::RCP< Epetra_Vector > &D)

Add a user supplied diagonal scale vector to the scaling object. ";

%feature("docstring")  NOX::Epetra::Scaling::addRowSumScaling "void
NOX::Epetra::Scaling::addRowSumScaling(ScaleType type, const
Teuchos::RCP< Epetra_Vector > &D)

Add \"Row Sum\" scaling to the scaling object. The supplied vector is
used to store the current row sum vector. ";

%feature("docstring")  NOX::Epetra::Scaling::addColSumScaling "void
NOX::Epetra::Scaling::addColSumScaling(ScaleType type, const
Teuchos::RCP< Epetra_Vector > &D)

Add \"Col Sum\" scaling to the scaling object. The supplied vector is
used to store the current column sum vector. ";

%feature("docstring")  NOX::Epetra::Scaling::computeScaling "void
NOX::Epetra::Scaling::computeScaling(const Epetra_LinearProblem
&problem)

Computes Row Sum scaling diagonal vectors. Only needs to be called if
a row or column sum scaling has been requested. ";

%feature("docstring")  NOX::Epetra::Scaling::scaleLinearSystem "void
NOX::Epetra::Scaling::scaleLinearSystem(Epetra_LinearProblem &problem)

Scales the linear system. ";

%feature("docstring")  NOX::Epetra::Scaling::unscaleLinearSystem "void NOX::Epetra::Scaling::unscaleLinearSystem(Epetra_LinearProblem
&problem)

Remove the scaling from the linear system. ";

%feature("docstring")  NOX::Epetra::Scaling::applyRightScaling "void
NOX::Epetra::Scaling::applyRightScaling(const Epetra_Vector &input,
Epetra_Vector &result)

Applies any RIGHT scaling vectors to an input vector. ";

%feature("docstring")  NOX::Epetra::Scaling::applyLeftScaling "void
NOX::Epetra::Scaling::applyLeftScaling(const Epetra_Vector &input,
Epetra_Vector &result)

Applies any LEFT scaling vectors to an input vector. ";

%feature("docstring")  NOX::Epetra::Scaling::print "void
NOX::Epetra::Scaling::print(ostream &os)

Printing. ";


// File: classNOX_1_1SharedObject.xml
%feature("docstring") NOX::SharedObject "

Holder for objects that are shared between NOX::Abstract::Groups.

Due to the size of certain operators (i.e. the Jacobian and
Preconditioning Matrices), we cannot afford to have multiple copies.
Instead we implement a shared object class that that all groups use.

C++ includes: NOX_SharedObjectTemplate.H ";

%feature("docstring")  NOX::SharedObject::SharedObject "NOX::SharedObject< Object, Owner >::SharedObject(const Teuchos::RCP<
Object > &newObject)

Constructor. ";

%feature("docstring")  NOX::SharedObject::~SharedObject "NOX::SharedObject< Object, Owner >::~SharedObject()

Destructor. ";

%feature("docstring")  NOX::SharedObject::getObject "Teuchos::RCP<Object> NOX::SharedObject< Object, Owner
>::getObject(const Owner *newOwner)

Get a non-const reference to the underlying object. ";

%feature("docstring")  NOX::SharedObject::getObject "Teuchos::RCP<const Object> NOX::SharedObject< Object, Owner
>::getObject() const

Return a const reference to the underlying object. ";

%feature("docstring")  NOX::SharedObject::isOwner "bool
NOX::SharedObject< Object, Owner >::isOwner(const Owner *checkOwner)
const

Return true if testOwner is the owner of the Jacobian. ";


// File: classNOX_1_1LineSearch_1_1Utils_1_1Slope.xml
%feature("docstring") NOX::LineSearch::Utils::Slope "

Common line search utilites for computing the slope of a function.

This class provides routines for computing the slope of a give
function. There are two methods, one that uses a Jacobian and the
other that estimates the action of the Jacobian by directional
derivatives.

C++ includes: NOX_LineSearch_Utils_Slope.H ";

%feature("docstring")  NOX::LineSearch::Utils::Slope::Slope "NOX::LineSearch::Utils::Slope::Slope(const Teuchos::RCP<
NOX::GlobalData > &gd)

Default constructor. ";

%feature("docstring")  NOX::LineSearch::Utils::Slope::~Slope "NOX::LineSearch::Utils::Slope::~Slope()

Destructor. ";

%feature("docstring")  NOX::LineSearch::Utils::Slope::reset "void
NOX::LineSearch::Utils::Slope::reset(const Teuchos::RCP<
NOX::GlobalData > &gd)

Reset method. ";

%feature("docstring")  NOX::LineSearch::Utils::Slope::computeSlope "double NOX::LineSearch::Utils::Slope::computeSlope(const
NOX::Abstract::Vector &dir, const NOX::Abstract::Group &grp)

Compute the inner product of the given direction and the gradient
associated with the given group.

Calculates and returns \\\\[ \\\\zeta = d^T \\\\nabla f(x). \\\\]

Here $d$ represents the input parameter dir and $\\\\nabla f(x)$ is
the gradient associated with the given group. ";

%feature("docstring")
NOX::LineSearch::Utils::Slope::computeSlopeWithOutJac "double
NOX::LineSearch::Utils::Slope::computeSlopeWithOutJac(const
NOX::Abstract::Vector &dir, const NOX::Abstract::Group &grp)

This is a variant of the computeSlope() method above optimized to work
with out having to compute an explicit Jacobian.

Calculates and returns \\\\[ \\\\zeta = d^T \\\\nabla f(x) = d^TJ^TF
\\\\]

Here $d$ represents the input parameter dir $\\\\nabla f(x)$ is the
gradient associated with the given group (for nonlinear solves this
equates to $ J^TF $ where $ J $ is the Jacobian and $ F $ is the
original nonlinear function).

We can rewrite this equation as:

\\\\[ d^TJ^TF = F^TJd \\\\]

which allows us to use directional derivatives to estimate $ J^TF $:

\\\\[ F^TJd = F^T \\\\frac{F(x + \\\\eta d) - F(x)}{\\\\eta} \\\\]

This may allow for faster computations of the slope if the Jacobian is
expensive to evaluate.

where $\\\\eta$ is a scalar perturbation calculated by:

\\\\[ \\\\eta = \\\\lambda * (\\\\lambda + \\\\frac{\\\\|
x\\\\|}{\\\\| d\\\\|} ) \\\\]

$ \\\\lambda $ is a constant fixed at 1.0e-6. ";


// File: classNOX_1_1StatusTest_1_1Stagnation.xml
%feature("docstring") NOX::StatusTest::Stagnation "

Failure test based on the convergence rate between nonlinear
iterations.

This status test returns NOX::StatusTest::Failed if we fail to reduce
the norm of $F$ by a specified tolerance for n consecutive iterations.
In other words, if

\\\\[ \\\\frac{\\\\| F_k \\\\|}{\\\\| F_{k-1} \\\\|} \\\\geq {\\\\rm
tolerance}\\\\]

for n consecutive iterations, the status is set to
NOX::StatusTest::Failed and returned. Otherwise, the status is set to
NOX::StatusTest::Uncoverged and returned. Both n and the tolerance are
specified in the constructor, by n and tol, respectively.

Based on experience the following values are recommended:

For Newton solves: n = 50, tolerance = 1.0

For Newton solves with a line search: n = 15, tolerance = 0.99

C++ includes: NOX_StatusTest_Stagnation.H ";

%feature("docstring")  NOX::StatusTest::Stagnation::Stagnation "NOX::StatusTest::Stagnation::Stagnation(int n=50, double tol=1.0)

Constructor.

Parameters:
-----------

n:  - Number of consecutive nonlinear iterations

tol:  - Tolerance for stagnation test ";

%feature("docstring")  NOX::StatusTest::Stagnation::~Stagnation "NOX::StatusTest::Stagnation::~Stagnation()

Destructor. ";

%feature("docstring")  NOX::StatusTest::Stagnation::checkStatus "NOX::StatusTest::StatusType
NOX::StatusTest::Stagnation::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)

Tests the stopping criterion.

The nature of this test is such that it must be executed at every
nonlinear iteration, so we don't use the checkType argument. ";

%feature("docstring")  NOX::StatusTest::Stagnation::getStatus "NOX::StatusTest::StatusType NOX::StatusTest::Stagnation::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  NOX::StatusTest::Stagnation::print "ostream &
NOX::StatusTest::Stagnation::print(ostream &stream, int indent=0)
const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  NOX::StatusTest::Stagnation::getMaxNumSteps "int NOX::StatusTest::Stagnation::getMaxNumSteps() const

Returns the used specified number of steps that can consecutively fail
the tolerance test before the test returns a failed status. ";

%feature("docstring")  NOX::StatusTest::Stagnation::getCurrentNumSteps
"int NOX::StatusTest::Stagnation::getCurrentNumSteps() const

Returns the current number of steps that have consecutively failed the
tolerance test. ";

%feature("docstring")  NOX::StatusTest::Stagnation::getTolerance "double NOX::StatusTest::Stagnation::getTolerance() const

Returns the user specified tolerance. ";

%feature("docstring")  NOX::StatusTest::Stagnation::getConvRate "double NOX::StatusTest::Stagnation::getConvRate() const

Returns the current convergence rate. ";


// File: classNOX_1_1Direction_1_1SteepestDescent.xml
%feature("docstring") NOX::Direction::SteepestDescent "

Calculates the steepest descent direction.

Calculates the direction \\\\[ d = - \\\\gamma \\\\nabla f(x) =
-\\\\gamma J(x)^T F(x) \\\\]

This is the (scaled) gradient of the function $f(x) = \\\\frac{1}{2}
\\\\|F(x)\\\\|^2$. The valued of $\\\\gamma$ depends on the choice of
\"Scaling Type\" below.

Parameters

\"Direction\": \"Method\" = \"Steepest Descent\"

\"Direction\"/\"Steepest Descent\":

\"Scaling Type\" can be either of:

\"2-Norm\" - $ \\\\gamma = \\\\displaystyle\\\\frac{1}{\\\\|d\\\\|_2}
$ [default]

\"Quadratic Model Min\" - $ \\\\gamma =
\\\\displaystyle\\\\frac{\\\\|d\\\\|_2^2}{d^T J^T J d} $

\"F 2-Norm\" - $ \\\\gamma =
\\\\displaystyle\\\\frac{1}{\\\\|F(x)\\\\|_2} $

\"None\" - No scaling

C++ includes: NOX_Direction_SteepestDescent.H ";

%feature("docstring")
NOX::Direction::SteepestDescent::SteepestDescent "NOX::Direction::SteepestDescent::SteepestDescent(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")
NOX::Direction::SteepestDescent::~SteepestDescent "NOX::Direction::SteepestDescent::~SteepestDescent()

Destructor. ";

%feature("docstring")  NOX::Direction::SteepestDescent::reset "bool
NOX::Direction::SteepestDescent::reset(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params)

Reset direction based on possibly new parameters. ";

%feature("docstring")  NOX::Direction::SteepestDescent::compute "bool
NOX::Direction::SteepestDescent::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::Generic &solver)

Compute the direction vector, dir, for a specific method given the
current group, grp.

The grp is not const so that we can compute the F vector, the Jacobian
matrix, the Newton vector, and so on.

Const access to the solver is used for getting additional information
such as the past solution, the iteration number, and so on. ";

%feature("docstring")  NOX::Direction::SteepestDescent::compute "bool
NOX::Direction::SteepestDescent::compute(NOX::Abstract::Vector &dir,
NOX::Abstract::Group &grp, const NOX::Solver::LineSearchBased &solver)

Same as compute( NOX::Abstract::Vector&, NOX::Abstract::Group&, const
NOX::Solver::Generic&).

Enables direct support for line search based solvers for the purpose
of efficiency since the LineSearchBased object has a getStep()
function that some directions require.

If it is not redefined in the derived class, it will just call the
compute with the NOX::Solver::Generic argument. ";


// File: classNOX_1_1MeritFunction_1_1SumOfSquares.xml
%feature("docstring") NOX::MeritFunction::SumOfSquares "

Sum of squares merit function.

A basic merit function used in many nonlinear equation solvers:

\\\\[ f = \\\\frac{1}{2} \\\\| F(x) \\\\| ^2 \\\\]

Where the norm is the 2-Norm using the NOX::Abstract::Vector's inner
product.

This is the default merit function used in nox.

This merit function is taken from: J. E. Dennis Jr. and Robert B.
Schnabel, \"Numerical Methods for Unconstrained Optimization and
Nonlinear Equations,\" Prentice Hall, 1983

C++ includes: NOX_MeritFunction_SumOfSquares.H ";

%feature("docstring")  NOX::MeritFunction::SumOfSquares::SumOfSquares
"NOX::MeritFunction::SumOfSquares::SumOfSquares(const Teuchos::RCP<
NOX::Utils > &u)

Constructor. ";

%feature("docstring")  NOX::MeritFunction::SumOfSquares::~SumOfSquares
"NOX::MeritFunction::SumOfSquares::~SumOfSquares()

Destructor. ";

%feature("docstring")  NOX::MeritFunction::SumOfSquares::computef "double NOX::MeritFunction::SumOfSquares::computef(const
NOX::Abstract::Group &grp) const

Computes the merit function, $ f(x) = \\\\frac{1}{2}\\\\| F(x) \\\\|^2
$. ";

%feature("docstring")
NOX::MeritFunction::SumOfSquares::computeGradient "void
NOX::MeritFunction::SumOfSquares::computeGradient(const
NOX::Abstract::Group &group, NOX::Abstract::Vector &result) const

Computes the gradient, $ g = \\\\nabla f(x) = J(x)^T F(x) $. ";

%feature("docstring")  NOX::MeritFunction::SumOfSquares::computeSlope
"double NOX::MeritFunction::SumOfSquares::computeSlope(const
NOX::Abstract::Vector &dir, const NOX::Abstract::Group &grp) const

Computes the slope, $ s(x,d) = d^T \\\\nabla f(x) = d^T J(x)^T F(x) $.

If the Jacobian is not computed in the grp object, then the slope can
be approximated using directional derivatives. More information can be
found in the method computeSlopeWithoutJac. ";

%feature("docstring")
NOX::MeritFunction::SumOfSquares::computeQuadraticModel "double
NOX::MeritFunction::SumOfSquares::computeQuadraticModel(const
NOX::Abstract::Vector &dir, const NOX::Abstract::Group &grp) const

Computes the quadratic model, $ m(x,d) = f(x) + \\\\nabla f(x)^T d +
d^T \\\\nabla^2 f(x) d $.

We approximate $ \\\\nabla^2f(x) \\\\approx J^TJ $:

\\\\[ m(d) = f(x) + (J(x)^T F)^T d + \\\\frac{1}{2} d^T B d \\\\] ";

%feature("docstring")
NOX::MeritFunction::SumOfSquares::computeQuadraticMinimizer "void
NOX::MeritFunction::SumOfSquares::computeQuadraticMinimizer(const
NOX::Abstract::Group &grp, NOX::Abstract::Vector &result) const

Computes the vector in the steepest descent direction that minimizes,
the quadratic model.

Computes the vector result: \\\\[ result = \\\\frac{\\\\nabla f^T
\\\\nabla f}{\\\\nabla f^T B \\\\nabla f} \\\\nabla f = -\\\\frac{(J^T
F)^T (J^T F)}{(J J^T F)^T (J J^T F)} J^T F \\\\] ";

%feature("docstring")  NOX::MeritFunction::SumOfSquares::name "const
string & NOX::MeritFunction::SumOfSquares::name() const

Returns the name of the merit function. ";


// File: classNOX_1_1Solver_1_1TensorBased.xml
%feature("docstring") NOX::Solver::TensorBased "

Nonlinear solver based on a rank-1 tensor method.

Solves $F(x)=0$ using a rank-1 tensor method and a linesearch
globalization.

At the kth nonlinear iteration, the solver does the following:

Computes the tensor direction $ d_T $ by finding the root or smallest
magnitude minimizer of the local model \\\\[ M_T(x_k+d) = F_k + J_kd +
a_k(s_k^Td)^2, \\\\] where \\\\[ a_k = 2(F_{k-1} - F_k - J_ks_k) /
(s_k^Ts_k)^2 \\\\] and \\\\[ s_k = s_{k-1} - s_k. \\\\]

Modifies the step according to a global strategy and updates $x$ as
$x_{k+1} = x_k + d(\\\\lambda) $ via a linesearch method, where $
d(\\\\lambda) $ is some function of $ \\\\lambda $. For instance, the
curvilinear step $ d_{\\\\lambda T} $ is a function of the linesearch
parameter $ \\\\lambda $ and is a parametric step that spans the
directions of the tensor step and the Newton step. At $ \\\\lambda=1
$, the curvilinear step equals the full tensor step, and as $
\\\\lambda $ nears 0, the curvilinear step approaches the Newton
direction. This step provides a monotonic decrease in the norm of the
local tensor model as $ \\\\lambda $ varies from 0 to 1.

The solver iterates until the status tests (see NOX::StatusTest)
determine either failure or convergence.

Input Parameters

To use this solver, set the \"Nonlinear Solver\" parameter to be
\"Tensor Based\". Then, specify the following sublists with the
appropriate parameters as indicated below.

\"Direction\" - Sublist of the direction parameters, passed to the
NOX::Direction::Factory constructor. Defaults to an empty list.

\"Method\" - Name of the direction to be computed in this solver.
\"Tensor\" and \"Newton\" are the only two valid choices. A sublist by
this name specifies all of the parameters to be passed to the linear
solver. See below under \"Linear Solver\".

\"Rescue Bad Newton Solve\" (Boolean) - If the linear solve does not
meet the tolerance specified by the forcing term, then use the step
anyway. Defaults to true.

\"Linear Solver\" - Sublist for the specific linear solver parameters
that are passed to NOX::Abstract::Group::computeNewton() and
NOX::Abstract::Group::applyJacobianInverse(). \"Linear Solver\" is
itself a sublist of the list specified in \"Method\" above (i.e.,
\"Tensor\" or \"Newton\"). Below is a partial list of standard
parameters usually available in common linear solvers. Check with the
specific linear solver being used for other parameters.

\"Max Iterations\" - Maximum number of Arnoldi iterations (also max
Krylov space dimension)

\"Tolerance\" - Relative tolerance for solving local model [default =
1e-4]

\"Output Frequency\" - Print output at every number of iterations
[default = 20]

\"Line Search\" - Sublist of the line search parameters. Because the
tensor step is not guaranteed to be a descent direction on the
function, not all \"basic\" line search approaches would be
appropriate. Thus, the LineSearch classes available to Newton's method
(e.g., Polynomial, More-Thuente) are not used here. Instead, this
solver class approriately handles technical considerations for tensor
methods with its own set of global strategies. The following
parameters specify the specific options for this line search:

\"Method\" - Name of the line search available to tensor methods Valid
choices are:

\"Curvilinear\" - Backtrack along the \"curvilinear\" path that spans
the tensor direction and the Newton direction and that maintains
monotonicity on the tensor model. Recommended because it tends to be
more robust and efficient than the other choices. [Default]

\"Standard\" - Backtrack along tensor direction unless it is not a
descent direction, in which case backtrack along Newton direction.

\"Dual\" - Backtrack along both the Newton and tensor directions and
choose the better of the two.

\"Full Step\" - Only use the full step and do not backtrack along both
the Newton and tensor directions and choose the better of the two.

\"Lambda selection\" - Flag for how to calculate the next linesearch
parameter lambda. Valid choices are \"Quadratic\" and \"Halving\"
(default). Quadratic constructs a quadratic interpolating polynomial
from the last trial point and uses the minimum of this function as the
next trial lambda (bounded by 0.1). Halving divides the linesearch
parameter by 2 before each trial, which is simpler but tends to
generate longer steps than quadratic.

\"Default Step\" - Starting value of the linesearch parameter
(defaults to 1.0)

\"Minimum Step\" - Minimum acceptable linesearch parameter before the
linesearch terminates (defaults to 1.0e-12). If there are many
linesearch failures, then lowering this value is one thing to try.

\"Recovery Step Type\" - Determines the step size to take when the
line search fails. Choices are:

\"Constant\" [default] - Uses a constant value set in \"Recovery
Step\".

\"Last Computed Step\" - Uses the last value computed by the line
search algorithm.

\"Recovery Step\" - Step parameter to take when the line search fails
(defaults to value for \"Default Step\")

\"Max Iters\" - Maximum number of iterations (i.e., backtracks)

\"Solver Options\" - Sublist of general solver options.

\"User Defined Pre/Post Operator\" is supported. See
NOX::Parameter::PrePostOperator for more details.

Output Parameters

Every time solve() is called, a sublist for output parameters called
\"Output\" will be created and will contain the following parameters:

\"Nonlinear Iterations\" - Number of nonlinear iterations

\"2-Norm of Residual\" - L-2 norm of the final residual $ F(x_k) $.

References

B. W. Bader, Tensor-Krylov methods for solving large-scale systems of
nonlinear equations, Ph.D. Thesis, 2003, University of Colorado,
Boulder, Colorado.

B. W. Bader, Tensor-Krylov methods for solving large-scale systems of
nonlinear equations, submitted to SIAM J. Numer. Anal.

B. W. Bader and R. B. Schnabel, Curvilinear linesearch for tensor
methods, SISC, 25(2):604-622.

R. B. Schnabel and P. D. Frank, Tensor methods for nonlinear
equations, SIAM J. Numer. Anal., 21(5):815-843.

Brett Bader (SNL 9233)

C++ includes: NOX_Solver_TensorBased.H ";

%feature("docstring")  NOX::Solver::TensorBased::TensorBased "NOX::Solver::TensorBased::TensorBased(const Teuchos::RCP<
NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)

Constructor.

See reset() for description. ";

%feature("docstring")  NOX::Solver::TensorBased::~TensorBased "NOX::Solver::TensorBased::~TensorBased()

Destructor. ";

%feature("docstring")  NOX::Solver::TensorBased::reset "void
NOX::Solver::TensorBased::reset(const NOX::Abstract::Vector
&initialGuess, const Teuchos::RCP< NOX::StatusTest::Generic > &tests)

Resets the solver, sets a new status test, and sets a new initial
guess. ";

%feature("docstring")  NOX::Solver::TensorBased::reset "void
NOX::Solver::TensorBased::reset(const NOX::Abstract::Vector
&initialGuess)

Resets the solver and sets a new initial guess. ";

%feature("docstring")  NOX::Solver::TensorBased::getStatus "NOX::StatusTest::StatusType NOX::Solver::TensorBased::getStatus()

Check current convergence and failure status. ";

%feature("docstring")  NOX::Solver::TensorBased::step "NOX::StatusTest::StatusType NOX::Solver::TensorBased::step()

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")  NOX::Solver::TensorBased::solve "NOX::StatusTest::StatusType NOX::Solver::TensorBased::solve()

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")  NOX::Solver::TensorBased::getSolutionGroup "const NOX::Abstract::Group &
NOX::Solver::TensorBased::getSolutionGroup() const

Return a reference to the current solution group. ";

%feature("docstring")
NOX::Solver::TensorBased::getPreviousSolutionGroup "const
NOX::Abstract::Group &
NOX::Solver::TensorBased::getPreviousSolutionGroup() const

Return a reference to the previous solution group. ";

%feature("docstring")  NOX::Solver::TensorBased::getNumIterations "int NOX::Solver::TensorBased::getNumIterations() const

Get number of iterations. ";

%feature("docstring")  NOX::Solver::TensorBased::getList "const
Teuchos::ParameterList & NOX::Solver::TensorBased::getList() const

Return a reference to the solver parameters. ";

%feature("docstring")  NOX::Solver::TensorBased::getSolutionGroupPtr "virtual Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::TensorBased::getSolutionGroupPtr() const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Solver::TensorBased::getPreviousSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::TensorBased::getPreviousSolutionGroupPtr() const

Return a RCP to the previous solution group. ";

%feature("docstring")  NOX::Solver::TensorBased::getListPtr "virtual
Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::TensorBased::getListPtr() const

Return a RCP to the solver parameters. ";


// File: classNOX_1_1Solver_1_1TrustRegionBased.xml
%feature("docstring") NOX::Solver::TrustRegionBased "

Newton-like solver using a trust region.

Our goal is to solve: $ F(x) = 0, $ where $ F:\\\\Re^n \\\\rightarrow
\\\\Re^n $. Alternatively, we might say that we wish to solve

$ \\\\min f(x) \\\\equiv \\\\frac{1}{2} \\\\|F(x)\\\\|^2_2. $

The trust region subproblem (TRSP) at iteration $k$ is given by

$ \\\\min \\\\; m_k(s) \\\\equiv f_k + g_k^T d + \\\\frac{1}{2} d^T
B_k d, \\\\mbox{ s.t. } \\\\|d\\\\| \\\\leq \\\\Delta_k \\\\quad
\\\\mbox{(TRSP)} $

where

$ f_k = f(x_k) = \\\\frac{1}{2} \\\\|F(x_k)\\\\|^2_2 $,

$ g_k = \\\\nabla f(x_k) = J(x_k)^T F(x_k) $,

$ B_k = J(x_k)^T J(x_k) \\\\approx \\\\nabla^2 f(x_k) $,

$ J(x_k)$ is the Jacobian of $F$ at $x_k$, and

$ \\\\Delta_k $ is the trust region radius.

The \"improvement ratio\" for a given step $ s $ is defined as

$ \\\\rho = \\\\displaystyle\\\\frac{ f(x_k) - f(x_k + d) } { m_k(0) -
m_k(d) } $

An iteration consists of the following steps.

Compute Newton-like direction: $n$

Compute Cauchy-like direction: $c$

If this is the first iteration, initialize $\\\\Delta$ as follows: If
$\\\\|n\\\\|_2 < \\\\Delta_{\\\\min}$, then $\\\\Delta = 2
\\\\Delta_{\\\\min}$; else, $\\\\Delta = \\\\|n\\\\|_2$.

Initialize $\\\\rho = -1$

While $\\\\rho < \\\\rho_{\\\\min}$ and $\\\\Delta >
\\\\Delta_{\\\\min}$, do the following.

Compute the direction $d$ as follows:

If $\\\\|n\\\\|_2 < \\\\Delta$, then take a Newton step by setting $d
= n$

Otherwise if $\\\\|c\\\\|_2 > \\\\Delta$, then take a Cauchy step by
setting $d = \\\\displaystyle\\\\frac{\\\\Delta}{\\\\|c\\\\|_2} c$

Otherwise, take a Dog Leg step by setting $ d = (1-\\\\gamma) c +
\\\\gamma n $ where $ \\\\gamma = \\\\displaystyle\\\\frac {-c^T a +
\\\\sqrt{ (c^Ta)^2 - (c^Tc - \\\\Delta^2) a^Ta}}{a^Ta} $ with $a =
n-c$.

Set $x_{\\\\rm new} = x + d$ and calculate $f_{\\\\rm new}$

If $f_{\\\\rm new} \\\\geq f$, then $\\\\rho = -1$ Otherwise $ \\\\rho
= \\\\displaystyle \\\\frac {f - f_{\\\\rm new}} {| d^T J F +
\\\\frac{1}{2} (J d)^T (J d)|} $

Update the solution: $x = x_{\\\\rm new}$

Update trust region:

If $\\\\rho < \\\\rho_{\\\\rm s}$ and $\\\\|n\\\\|_2 < \\\\Delta$,
then shrink the trust region to the size of the Newton step:
$\\\\Delta = \\\\|n\\\\|_2$.

Otherwise if $\\\\rho < \\\\rho_{\\\\rm s}$, then shrink the trust
region: $\\\\Delta = \\\\max \\\\{ \\\\beta_{\\\\rm s} \\\\Delta,
\\\\Delta_{\\\\min} \\\\} $.

Otherwise if $\\\\rho > \\\\rho_{\\\\rm e}$ and $\\\\|d\\\\|_2 =
\\\\Delta$, then expand the trust region: $\\\\Delta = \\\\min \\\\{
\\\\beta_{\\\\rm e} \\\\Delta, \\\\Delta_{\\\\rm max} \\\\} $.

Input Paramters

The following parameters should be specified in the \"Trust Region\"
sublist based to the solver.

\"Direction\" - Sublist of the direction parameters for the Newton
point, passed to the NOX::Direction::Manager constructor. If this
sublist does not exist, it is created by default. Furthermore, if
\"Method\" is not specified in this sublist, it is added with a value
of \"Newton\".

\"Cauchy %Direction\" - Sublist of the direction parameters for the
Cauchy point, passed to the NOX::Direction::Manager constructor. If
this sublist does not exist, it is created by default. Furthremore, if
\"Method\" is not specified in this sublist, it is added with a value
of \"Steepest Descent\" Finally, if the sub-sublist \"Steepest
Descent\" does not exist, it is created and the parameter \"Scaling
Type\" is added and set to \"Quadratic\".

\"Minimum Trust Region Radius\" ( $\\\\Delta_{\\\\min}$) - Minimum
allowable trust region radius. Defaults to 1.0e-6.

\"Maximum Trust Region Radius\" ( $\\\\Delta_{\\\\max}$) - Maximum
allowable trust region radius. Defaults to 1.0e+10.

\"Minimum Improvement Ratio\" ( $\\\\rho_{\\\\min}$) - Minimum
improvement ratio to accept the step. Defaults to 1.0e-4.

\"Contraction Trigger Ratio\" ( $\\\\rho_{\\\\rm s}$) - If the
improvement ratio is less than this value, then the trust region is
contracted by the amount specified by the \"Contraction Factor\". Must
be larger than \"Minimum Improvement Ratio\". Defaults to 0.1.

\"Contraction Factor\" ( $\\\\beta_{\\\\rm s}$) - See above. Defaults
to 0.25.

\"Expansion Trigger Ratio\" ( $\\\\rho_{\\\\rm e}$) - If the
improvement ratio is greater than this value, then the trust region is
contracted by the amount specified by the \"Expansion Factor\".
Defaults to 0.75.

\"Expansion Factor\" ( $\\\\beta_{\\\\rm e}$) - See above. Defaults to
4.0.

\"Recovery Step\" - Defaults to 1.0.

\"Use Ared/Pred Ratio Calculation\" (boolean) - Defaults to false. If
set to true, this option replaces the algorithm used to compute the
improvement ratio, $ \\\\rho $, as described above. The improvement
ratio is replaced by an \"Ared/Pred\" sufficient decrease criteria
similar to that used in line search algorithms (see Eisenstat and
Walker, SIAM Journal on Optimization V4 no. 2 (1994) pp 393-422):
$\\\\rho = \\\\frac{\\\\|F(x) \\\\| - \\\\| F(x + d) \\\\| } {\\\\|
F(x) \\\\| - \\\\| F(x) + Jd \\\\| } $

\"Solver Options\" - Sublist of general solver options. \"User Defined
Pre/Post Operator\" is supported. See NOX::Parameter::PrePostOperator
for more details.

Output Paramters

A sublist for output parameters called \"Output\" will be created and
contain the following parameters:

\"Nonlinear Iterations\" - Number of nonlinear iterations

\"2-Norm or Residual\" - Two-norm of final residual

Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)

C++ includes: NOX_Solver_TrustRegionBased.H ";

%feature("docstring")  NOX::Solver::TrustRegionBased::TrustRegionBased
"TrustRegionBased::TrustRegionBased(const Teuchos::RCP<
NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params)

Constructor.

See reset() for description. ";

%feature("docstring")
NOX::Solver::TrustRegionBased::~TrustRegionBased "TrustRegionBased::~TrustRegionBased()

Destructor. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::reset "void
TrustRegionBased::reset(const NOX::Abstract::Vector &initialGuess,
const Teuchos::RCP< NOX::StatusTest::Generic > &tests)

Resets the solver, sets a new status test, and sets a new initial
guess. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::reset "void
TrustRegionBased::reset(const NOX::Abstract::Vector &initialGuess)

Resets the solver and sets a new initial guess. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::getStatus "NOX::StatusTest::StatusType TrustRegionBased::getStatus()

Check current convergence and failure status. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::step "NOX::StatusTest::StatusType TrustRegionBased::step()

Do one nonlinear step in the iteration sequence and return status. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::solve "NOX::StatusTest::StatusType TrustRegionBased::solve()

Solve the nonlinear problem and return final status.

By \"solve\", we call iterate() until the NOX::StatusTest value is
either NOX::StatusTest::Converged or NOX::StatusTest::Failed. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::getSolutionGroup
"const Abstract::Group & TrustRegionBased::getSolutionGroup() const

Return a reference to the current solution group. ";

%feature("docstring")
NOX::Solver::TrustRegionBased::getPreviousSolutionGroup "const
Abstract::Group & TrustRegionBased::getPreviousSolutionGroup() const

Return a reference to the previous solution group. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::getNumIterations
"int TrustRegionBased::getNumIterations() const

Get number of iterations. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::getList "const
Teuchos::ParameterList & TrustRegionBased::getList() const

Return a reference to the solver parameters. ";

%feature("docstring")
NOX::Solver::TrustRegionBased::getSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::TrustRegionBased::getSolutionGroupPtr() const

Return a RCP to the solution group. ";

%feature("docstring")
NOX::Solver::TrustRegionBased::getPreviousSolutionGroupPtr "virtual
Teuchos::RCP< const NOX::Abstract::Group >
NOX::Solver::TrustRegionBased::getPreviousSolutionGroupPtr() const

Return a RCP to the previous solution group. ";

%feature("docstring")  NOX::Solver::TrustRegionBased::getListPtr "virtual Teuchos::RCP< const Teuchos::ParameterList >
NOX::Solver::TrustRegionBased::getListPtr() const

Return a RCP to the solver parameters. ";


// File: classNOX_1_1Direction_1_1UserDefinedFactory.xml
%feature("docstring") NOX::Direction::UserDefinedFactory "

Pure virtual interface for users to supply their own direction
objects.

C++ includes: NOX_Direction_UserDefinedFactory.H ";

%feature("docstring")
NOX::Direction::UserDefinedFactory::UserDefinedFactory "NOX::Direction::UserDefinedFactory::UserDefinedFactory()

Constructor. ";

%feature("docstring")
NOX::Direction::UserDefinedFactory::~UserDefinedFactory "virtual
NOX::Direction::UserDefinedFactory::~UserDefinedFactory()

Destructor. ";

%feature("docstring")
NOX::Direction::UserDefinedFactory::buildDirection "virtual
Teuchos::RCP<NOX::Direction::Generic>
NOX::Direction::UserDefinedFactory::buildDirection(const Teuchos::RCP<
NOX::GlobalData > &gd, Teuchos::ParameterList &params) const =0

Builds a user defined direction object.

Parameters:
-----------

gd:  A global data pointer that contains the top level parameter list.
Without storing this inside the direction object, there is no
guarantee that the second parameter params will still exist. It can be
deleted by the top level RCP.

params:  Sublist with direction construction parameters. ";


// File: classNOX_1_1LineSearch_1_1UserDefinedFactory.xml
%feature("docstring") NOX::LineSearch::UserDefinedFactory "

Pure virtual interface for users to supply their own line search
objects.

C++ includes: NOX_LineSearch_UserDefinedFactory.H ";

%feature("docstring")
NOX::LineSearch::UserDefinedFactory::UserDefinedFactory "NOX::LineSearch::UserDefinedFactory::UserDefinedFactory()

Constructor. ";

%feature("docstring")
NOX::LineSearch::UserDefinedFactory::~UserDefinedFactory "virtual
NOX::LineSearch::UserDefinedFactory::~UserDefinedFactory()

Destructor. ";

%feature("docstring")
NOX::LineSearch::UserDefinedFactory::buildLineSearch "virtual
Teuchos::RCP<NOX::LineSearch::Generic>
NOX::LineSearch::UserDefinedFactory::buildLineSearch(const
Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)
const =0

Builds a user defined line search object.

Parameters:
-----------

gd:  A global data pointer that contains the top level parameter list.
Without storing this inside the line searchobject, there is no
guarantee that the second parameter params will still exist. It can be
deleted by the top level RCP.

params:  Sublist with line search construction parameters. ";


// File: classNOX_1_1LineSearch_1_1UserDefinedFactoryT.xml
%feature("docstring") NOX::LineSearch::UserDefinedFactoryT "

Concrete instantiation of a NOX::LineSearch::UserDefinedFactory object
that uses the base objects only for constuction.

If the user writes their own line search and that object has the same
constructor arguments as the nox line searches (the gd and params as
in the buildDirection method), then users can use this object instead
of having to write their own factory.

For example, if a user writes their own line search object:

They can build that object using this factory and do not have to write
their own factory

It is critical that the user defined factory be set in the parameter
list as a base class type object: NOX::LineSearch::UserDefinedFactory.

C++ includes: NOX_LineSearch_UserDefinedFactoryT.H ";

%feature("docstring")
NOX::LineSearch::UserDefinedFactoryT::UserDefinedFactoryT "NOX::LineSearch::UserDefinedFactoryT< T >::UserDefinedFactoryT()

Constructor. ";

%feature("docstring")
NOX::LineSearch::UserDefinedFactoryT::~UserDefinedFactoryT "NOX::LineSearch::UserDefinedFactoryT< T >::~UserDefinedFactoryT()

Destructor. ";

%feature("docstring")
NOX::LineSearch::UserDefinedFactoryT::buildLineSearch "Teuchos::RCP<NOX::LineSearch::Generic>
NOX::LineSearch::UserDefinedFactoryT< T >::buildLineSearch(const
Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)
const

Builds a user defined line search object.

Parameters:
-----------

gd:  A global data pointer that contains the top level parameter list.
Without storing this inside the line searchobject, there is no
guarantee that the second parameter params will still exist. It can be
deleted by the top level RCP.

params:  Sublist with line search construction parameters. ";


// File: classNOX_1_1Direction_1_1UserDefinedFactoryT.xml
%feature("docstring") NOX::Direction::UserDefinedFactoryT "

Concrete instantiation of a NOX::Direction::UserDefinedFactory object
that uses the base objects only for constuction.

If the user writes their own direction and that object has the same
constructor arguments as the nox directions (the gd and params as in
the buildDirection method), then users can use this object instead of
having to write their own factory.

For example, if a user writes their own direction object:

They can build that object using this factory and do not have to write
their own factory

It is critical that the user defined factory be set in the parameter
list as a base class type object: NOX::Direction::UserDefinedFactory.

C++ includes: NOX_Direction_UserDefinedFactoryT.H ";

%feature("docstring")
NOX::Direction::UserDefinedFactoryT::UserDefinedFactoryT "NOX::Direction::UserDefinedFactoryT< T >::UserDefinedFactoryT()

Constructor. ";

%feature("docstring")
NOX::Direction::UserDefinedFactoryT::~UserDefinedFactoryT "NOX::Direction::UserDefinedFactoryT< T >::~UserDefinedFactoryT()

Destructor. ";

%feature("docstring")
NOX::Direction::UserDefinedFactoryT::buildDirection "Teuchos::RCP<NOX::Direction::Generic>
NOX::Direction::UserDefinedFactoryT< T >::buildDirection(const
Teuchos::RCP< NOX::GlobalData > &gd, Teuchos::ParameterList &params)
const

Builds a user defined direction object.

Parameters:
-----------

gd:  A global data pointer that contains the top level parameter list.
Without storing this inside the direction object, there is no
guarantee that the second parameter params will still exist. It can be
deleted by the top level RCP.

params:  Sublist with direction construction parameters. ";


// File: classNOX_1_1Utils.xml
%feature("docstring") NOX::Utils "

Provides printing utilities.

The following parameters are valid for this class and should be
defined in the \"Printing\" sublist of the solver parameter list.

\"Output Information\" - Can be a sublist or an integer. If it is an
integer, this value is a sum of MsgType's to specify how much
information to show. Defaults to NOX::Utils::Warning +
NOX::Utils::OuterIteration + NOX::Utils::InnerIteration +
NOX::Utils::Parameters = 0xf = 15. If this is a sublist, the following
booleans are valid (set to true to enable print option): \"Error\"

\"Warning\"

\"Outer Iteration\"

\"Inner Iteration\"

\"Parameters\"

\"Details\"

\"Outer Iteration StatusTest\"

\"Linear Solver Details\"

\"Test Details\"

\"Stepper Iteration\"

\"Stepper Details\"

\"Stepper Parameters\"

\"Debug\"

\"Output Processor\" - Specifies the designated print process.
Defaults to 0.

\"MyPID\" - Specifies this process's ID. Defaults to 0.

\"Output Precision\" - Specifis the default number of decimal places
to be used when printing floating point numbers. The default is 4.

\"Output Stream\" - A Teuchos::RCP<std::ostream> object that will be
used for standard output.

\"Error Stream\" - A Teuchos::RCP<std::ostream> object that will be
used for output of error information.

The public variables should never be modified directly.

C++ includes: NOX_Utils.H ";

%feature("docstring")  NOX::Utils::Utils "NOX::Utils::Utils(int
outputInformation=0xf, int MyPID=0, int outputProcess=0, int
outputPrecision=3, const Teuchos::RCP< std::ostream >
&outputStream=Teuchos::null, const Teuchos::RCP< std::ostream >
&errorStream=Teuchos::null)

Constructor.

The final two arguments are a reference counted pointers to ostreams.
This defaults to std::cout and std::cerr if not supplied. Users should
only supply this argument if they want output directed to a
std::ostream other than the defaults. If so, they must wrap the
ostream in a reference counted pointer for safe memory management. ";

%feature("docstring")  NOX::Utils::Utils "NOX::Utils::Utils(Teuchos::ParameterList &p)

Constructor via a parameter list. ";

%feature("docstring")  NOX::Utils::Utils "NOX::Utils::Utils(const
NOX::Utils &u)

Copy constructor. ";

%feature("docstring")  NOX::Utils::~Utils "NOX::Utils::~Utils()

Destructor. ";

%feature("docstring")  NOX::Utils::reset "void
NOX::Utils::reset(Teuchos::ParameterList &p)

Reset the printing parameters. ";

%feature("docstring")  NOX::Utils::isPrintType "bool
NOX::Utils::isPrintType(NOX::Utils::MsgType type) const

Returns true if this is a valid print type. ";

%feature("docstring")  NOX::Utils::out "std::ostream &
NOX::Utils::out() const

Returns the ostream for printing if this proces is the print process.
Returns a Teuchos::oblackholestream otherwise. ";

%feature("docstring")  NOX::Utils::out "std::ostream &
NOX::Utils::out(NOX::Utils::MsgType type) const

Returns the ostream for printing if this process is the print process
and the print type is valid. Returns a Teuchos::oblackholestream
otherwise. ";

%feature("docstring")  NOX::Utils::pout "std::ostream &
NOX::Utils::pout() const

Returns the ostream for printing regardless of the print processor.
Only use this call if you want all processes to print to the ostream.
";

%feature("docstring")  NOX::Utils::pout "std::ostream &
NOX::Utils::pout(NOX::Utils::MsgType type) const

Returns the ostream for printing if the print type matches. Only use
this call if you want all processes to print to the ostream for the
print type. ";

%feature("docstring")  NOX::Utils::err "std::ostream &
NOX::Utils::err() const

Returns the error stream for printing if this is the print process. ";

%feature("docstring")  NOX::Utils::perr "std::ostream &
NOX::Utils::perr() const

Returns the error stream for printing to all processors. Only use this
call if you want all processes to print to the error stream. ";

%feature("docstring")  NOX::Utils::print "void
NOX::Utils::print(ostream &os) const

Print this object. ";

%feature("docstring")  NOX::Utils::sciformat "NOX::Utils::Sci
NOX::Utils::sciformat(double dval) const

Creates a Sci object which can be used in an output stream for
printing a double precision number in scientific format with an
arbitrary precision. The precision is that specificed by the Utils
object.

For example,

This is modeled after the Form and Bound_form objects in Stroustrup,
C++ Programming Langauge, 3rd ed., Chapter 21.4. ";


// File: classNOX_1_1Utils_1_1Fill.xml
%feature("docstring") NOX::Utils::Fill "

Fill object - used to print the given character the number of times
specified.

C++ includes: NOX_Utils.H ";

%feature("docstring")  NOX::Utils::Fill::Fill "NOX::Utils::Fill::Fill(int ntimes, char ch)

Constructor. ";

%feature("docstring")  NOX::Utils::Fill::~Fill "NOX::Utils::Fill::~Fill()

Destructor. ";


// File: classNOX_1_1Utils_1_1Sci.xml
%feature("docstring") NOX::Utils::Sci "

Sci object - used to print the given value with the specified
precision.

C++ includes: NOX_Utils.H ";

%feature("docstring")  NOX::Utils::Sci::Sci "NOX::Utils::Sci::Sci(double val, int precision=-1)

Constructor. ";

%feature("docstring")  NOX::Utils::Sci::~Sci "NOX::Utils::Sci::~Sci()

Destructor. ";


// File: classNOX_1_1Epetra_1_1Vector.xml
%feature("docstring") NOX::Epetra::Vector "

Implementation of NOX::Abstract::Vector for Epetra vectors.

C++ includes: NOX_Epetra_Vector.H ";

%feature("docstring")  NOX::Epetra::Vector::getEpetraVector "Epetra_Vector & NOX::Epetra::Vector::getEpetraVector()

Get reference to underlying Epetra vector. ";

%feature("docstring")  NOX::Epetra::Vector::getEpetraVector "const
Epetra_Vector & NOX::Epetra::Vector::getEpetraVector() const

Get const reference to underlying Epetra vector. ";

%feature("docstring")  NOX::Epetra::Vector::init "NOX::Abstract::Vector & NOX::Epetra::Vector::init(double gamma)

Initialize every element of this vector with gamma.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\gamma \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\] Reference to
this object ";

%feature("docstring")  NOX::Epetra::Vector::random "NOX::Abstract::Vector & NOX::Epetra::Vector::random(bool
useSeed=false, int seed=1)

Initialize each element of this vector with a random value.

If useSeed is true, uses the value of seed to seed the random number
generator before filling the entries of this vector. So, if two calls
are made where useSeed is true and seed is the same, then the vectors
returned should be the same.

Default implementation throw an error. Only referenced by LOCA
methods.

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::abs "NOX::Abstract::Vector & NOX::Epetra::Vector::abs(const
NOX::Epetra::Vector &y) ";

%feature("docstring")  NOX::Epetra::Vector::abs "NOX::Abstract::Vector & NOX::Epetra::Vector::abs(const
NOX::Abstract::Vector &y)

Put element-wise absolute values of source vector y into this vector.

Here x represents this vector, and we update it as \\\\[ x_i = | y_i |
\\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::reciprocal "NOX::Abstract::Vector & NOX::Epetra::Vector::reciprocal(const
NOX::Epetra::Vector &y) ";

%feature("docstring")  NOX::Epetra::Vector::reciprocal "NOX::Abstract::Vector & NOX::Epetra::Vector::reciprocal(const
NOX::Abstract::Vector &y)

Put element-wise reciprocal of source vector y into this vector.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\frac{1}{y_i} \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::scale "NOX::Abstract::Vector & NOX::Epetra::Vector::scale(double gamma)

Scale each element of this vector by gamma.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\gamma x_i \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::scale "NOX::Abstract::Vector & NOX::Epetra::Vector::scale(const
NOX::Epetra::Vector &a) ";

%feature("docstring")  NOX::Epetra::Vector::scale "NOX::Abstract::Vector & NOX::Epetra::Vector::scale(const
NOX::Abstract::Vector &a)

Scale this vector element-by-element by the vector a.

Here x represents this vector, and we update it as \\\\[ x_i = x_i
\\\\cdot a_i \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::update "NOX::Abstract::Vector & NOX::Epetra::Vector::update(double alpha,
const NOX::Epetra::Vector &a, double gamma=0.0) ";

%feature("docstring")  NOX::Epetra::Vector::update "NOX::Abstract::Vector & NOX::Epetra::Vector::update(double alpha,
const NOX::Abstract::Vector &a, double gamma=0.0)

Compute x = (alpha * a) + (gamma * x) where x is this vector.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\alpha \\\\; a_i + \\\\gamma \\\\; x_i \\\\quad \\\\mbox{for }
i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::update "NOX::Abstract::Vector & NOX::Epetra::Vector::update(double alpha,
const NOX::Epetra::Vector &a, double beta, const NOX::Epetra::Vector
&b, double gamma=0.0) ";

%feature("docstring")  NOX::Epetra::Vector::update "NOX::Abstract::Vector & NOX::Epetra::Vector::update(double alpha,
const NOX::Abstract::Vector &a, double beta, const
NOX::Abstract::Vector &b, double gamma=0.0)

Compute x = (alpha * a) + (beta * b) + (gamma * x) where x is this
vector.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\alpha \\\\; a_i + \\\\beta \\\\; b_i + \\\\gamma \\\\; x_i
\\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Epetra::Vector::clone "Teuchos::RCP<
NOX::Abstract::Vector > NOX::Epetra::Vector::clone(CopyType
type=DeepCopy) const

Create a new Vector of the same underlying type by cloning \"this\",
and return a pointer to the new vector.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\" (the memory is allocated for the objects, but
the current values are not copied into the vector). Note that there is
no assumption that a vector created by ShapeCopy is initialized to
zeros.

Pointer to newly created vector or NULL if clone is not supported. ";

%feature("docstring")  NOX::Epetra::Vector::createMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
NOX::Epetra::Vector::createMultiVector(const NOX::Abstract::Vector
*const *vecs, int numVecs, NOX::CopyType type=NOX::DeepCopy) const

Create a MultiVector with numVecs+1 columns out of an array of
Vectors. The vector stored under this will be the first column with
the remaining numVecs columns given by vecs.

The implementation here creates a NOX::Epetra::MultiVector with either
Shape or Deep copies of the supplied vectors. ";

%feature("docstring")  NOX::Epetra::Vector::createMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
NOX::Epetra::Vector::createMultiVector(int numVecs, NOX::CopyType
type=NOX::DeepCopy) const

Create a MultiVector with numVecs columns.

The implementation here creates a NOX::Epetra::MultiVector with either
Shape or Deep copies of the supplied vector. ";

%feature("docstring")  NOX::Epetra::Vector::norm "double
NOX::Epetra::Vector::norm(NOX::Abstract::Vector::NormType
type=TwoNorm) const

Norm.

Here x represents this vector, and we compute its norm as follows: for
each NOX::Abstract::Vector::NormType:  NOX::Abstract::Vector::TwoNorm
\\\\[ \\\\|x\\\\| = \\\\sqrt{\\\\sum_{i=1}^{n} x_i^2} \\\\]

NOX::Abstract::Vector::OneNorm \\\\[ \\\\|x\\\\| = \\\\sum_{i=1}^{n}
|x_i| \\\\]

NOX::Abstract::Vector::MaxNorm \\\\[ \\\\|x\\\\| = \\\\max_{i} |x_i|
\\\\]

$\\\\|x\\\\|$ ";

%feature("docstring")  NOX::Epetra::Vector::norm "double
NOX::Epetra::Vector::norm(const NOX::Epetra::Vector &weights) const ";

%feature("docstring")  NOX::Epetra::Vector::norm "double
NOX::Epetra::Vector::norm(const NOX::Abstract::Vector &weights) const

Weighted 2-Norm.

Here x represents this vector, and we compute its weighted norm as
follows: \\\\[ \\\\|x\\\\|_w = \\\\sqrt{\\\\sum_{i=1}^{n} w_i \\\\;
x_i^2} \\\\]  $ \\\\|x\\\\|_w $ ";

%feature("docstring")  NOX::Epetra::Vector::innerProduct "double
NOX::Epetra::Vector::innerProduct(const NOX::Epetra::Vector &y) const
";

%feature("docstring")  NOX::Epetra::Vector::innerProduct "double
NOX::Epetra::Vector::innerProduct(const NOX::Abstract::Vector &y)
const

Inner product with y.

Here x represents this vector, and we compute its inner product with y
as follows: \\\\[ \\\\langle x,y \\\\rangle = \\\\sum_{i=1}^n x_i y_i
\\\\]  $\\\\langle x,y \\\\rangle$ ";

%feature("docstring")  NOX::Epetra::Vector::Vector "NOX::Epetra::Vector::Vector(const Teuchos::RCP< Epetra_Vector >
&source, NOX::Epetra::Vector::MemoryType
memoryType=NOX::Epetra::Vector::CreateCopy, NOX::CopyType
type=NOX::DeepCopy, Teuchos::RCP< NOX::Epetra::VectorSpace >
vs=Teuchos::null)

Constructor that creates a COPY or VIEW of the Epetra_Vector.

NOTE: This ctor should just always create a view. It should be
implicit from the fact that a RCP object is being passed in that a
persisting relationship is present. However, since this could cause
confusion, the default is to make a copy and if a user wants a view,
they must pass in an explicit flag.

A VIEW of a vector uses the same underlying memory. WARNING: A View
can be dangerous since multiple objects can access the same memory
locations. ";

%feature("docstring")  NOX::Epetra::Vector::Vector "NOX::Epetra::Vector::Vector(const Epetra_Vector &source, NOX::CopyType
type=NOX::DeepCopy, Teuchos::RCP< NOX::Epetra::VectorSpace >
vs=Teuchos::null)

Construct by copying map and/or elements of an Epetra_Vector.

Allocates an entirely new vector. Does NOT allow for a view. ";

%feature("docstring")  NOX::Epetra::Vector::Vector "NOX::Epetra::Vector::Vector(const NOX::Epetra::Vector &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  NOX::Epetra::Vector::~Vector "NOX::Epetra::Vector::~Vector()

Destruct Vector. ";

%feature("docstring")  NOX::Epetra::Vector::length "int
NOX::Epetra::Vector::length() const

Return the length of vector.

The length of this vector

Even if the vector is distributed across processors, this should
return the  global length of the vector. ";

%feature("docstring")  NOX::Epetra::Vector::print "void
NOX::Epetra::Vector::print(std::ostream &stream) const

Print the vector. To be used for debugging only. ";

%feature("docstring")  NOX::Epetra::Vector::getVectorSpace "Teuchos::RCP< NOX::Epetra::VectorSpace >
NOX::Epetra::Vector::getVectorSpace() const

Returns the NOX::Epetra::VectorSpace associated with this vector. ";


// File: classNOX_1_1Abstract_1_1Vector.xml
%feature("docstring") NOX::Abstract::Vector "

NOX's pure abstract vector interface for vectors that are used by the
nonlinear solver.

This class is a member of the namespace NOX::Abstract.

The user should implement their own concrete implementation of this
class or use one of the implementations provided by us.

Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)

C++ includes: NOX_Abstract_Vector.H ";

%feature("docstring")  NOX::Abstract::Vector::init "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::init(double gamma)=0

Initialize every element of this vector with gamma.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\gamma \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\] Reference to
this object ";

%feature("docstring")  NOX::Abstract::Vector::random "NOX::Abstract::Vector & NOX::Abstract::Vector::random(bool
useSeed=false, int seed=1)

Initialize each element of this vector with a random value.

If useSeed is true, uses the value of seed to seed the random number
generator before filling the entries of this vector. So, if two calls
are made where useSeed is true and seed is the same, then the vectors
returned should be the same.

Default implementation throw an error. Only referenced by LOCA
methods.

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::abs "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::abs(const
NOX::Abstract::Vector &y)=0

Put element-wise absolute values of source vector y into this vector.

Here x represents this vector, and we update it as \\\\[ x_i = | y_i |
\\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::reciprocal "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::reciprocal(const
NOX::Abstract::Vector &y)=0

Put element-wise reciprocal of source vector y into this vector.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\frac{1}{y_i} \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::scale "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::scale(double gamma)=0

Scale each element of this vector by gamma.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\gamma x_i \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::scale "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::scale(const
NOX::Abstract::Vector &a)=0

Scale this vector element-by-element by the vector a.

Here x represents this vector, and we update it as \\\\[ x_i = x_i
\\\\cdot a_i \\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::update "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::update(double alpha,
const NOX::Abstract::Vector &a, double gamma=0.0)=0

Compute x = (alpha * a) + (gamma * x) where x is this vector.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\alpha \\\\; a_i + \\\\gamma \\\\; x_i \\\\quad \\\\mbox{for }
i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::update "virtual
NOX::Abstract::Vector& NOX::Abstract::Vector::update(double alpha,
const NOX::Abstract::Vector &a, double beta, const
NOX::Abstract::Vector &b, double gamma=0.0)=0

Compute x = (alpha * a) + (beta * b) + (gamma * x) where x is this
vector.

Here x represents this vector, and we update it as \\\\[ x_i =
\\\\alpha \\\\; a_i + \\\\beta \\\\; b_i + \\\\gamma \\\\; x_i
\\\\quad \\\\mbox{for } i=1,\\\\dots,n \\\\]

Reference to this object ";

%feature("docstring")  NOX::Abstract::Vector::clone "virtual
Teuchos::RCP<NOX::Abstract::Vector>
NOX::Abstract::Vector::clone(NOX::CopyType type=NOX::DeepCopy) const
=0

Create a new Vector of the same underlying type by cloning \"this\",
and return a pointer to the new vector.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\" (the memory is allocated for the objects, but
the current values are not copied into the vector). Note that there is
no assumption that a vector created by ShapeCopy is initialized to
zeros.

Pointer to newly created vector or NULL if clone is not supported. ";

%feature("docstring")  NOX::Abstract::Vector::createMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
NOX::Abstract::Vector::createMultiVector(const NOX::Abstract::Vector
*const *vecs, int numVecs, NOX::CopyType type=NOX::DeepCopy) const

Create a MultiVector with numVecs+1 columns out of an array of
Vectors. The vector stored under this will be the first column with
the remaining numVecs columns given by vecs.

The default implementation creates a generic NOX::MultiVector with
either Shape or Deep copies of the supplied vectors. ";

%feature("docstring")  NOX::Abstract::Vector::createMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
NOX::Abstract::Vector::createMultiVector(int numVecs, NOX::CopyType
type=NOX::DeepCopy) const

Create a MultiVector with numVecs columns.

The default implementation creates a generic NOX::MultiVector with
either Shape or Deep copies of the supplied vector. ";

%feature("docstring")  NOX::Abstract::Vector::norm "virtual double
NOX::Abstract::Vector::norm(NOX::Abstract::Vector::NormType
type=NOX::Abstract::Vector::TwoNorm) const =0

Norm.

Here x represents this vector, and we compute its norm as follows: for
each NOX::Abstract::Vector::NormType:  NOX::Abstract::Vector::TwoNorm
\\\\[ \\\\|x\\\\| = \\\\sqrt{\\\\sum_{i=1}^{n} x_i^2} \\\\]

NOX::Abstract::Vector::OneNorm \\\\[ \\\\|x\\\\| = \\\\sum_{i=1}^{n}
|x_i| \\\\]

NOX::Abstract::Vector::MaxNorm \\\\[ \\\\|x\\\\| = \\\\max_{i} |x_i|
\\\\]

$\\\\|x\\\\|$ ";

%feature("docstring")  NOX::Abstract::Vector::norm "virtual double
NOX::Abstract::Vector::norm(const NOX::Abstract::Vector &weights)
const =0

Weighted 2-Norm.

Here x represents this vector, and we compute its weighted norm as
follows: \\\\[ \\\\|x\\\\|_w = \\\\sqrt{\\\\sum_{i=1}^{n} w_i \\\\;
x_i^2} \\\\]  $ \\\\|x\\\\|_w $ ";

%feature("docstring")  NOX::Abstract::Vector::innerProduct "virtual
double NOX::Abstract::Vector::innerProduct(const NOX::Abstract::Vector
&y) const =0

Inner product with y.

Here x represents this vector, and we compute its inner product with y
as follows: \\\\[ \\\\langle x,y \\\\rangle = \\\\sum_{i=1}^n x_i y_i
\\\\]  $\\\\langle x,y \\\\rangle$ ";

%feature("docstring")  NOX::Abstract::Vector::Vector "NOX::Abstract::Vector::Vector()

Abstract Vector constructor (does nothing) ";

%feature("docstring")  NOX::Abstract::Vector::~Vector "virtual
NOX::Abstract::Vector::~Vector()

Abstract Vector destructor (does nothing) ";

%feature("docstring")  NOX::Abstract::Vector::length "virtual int
NOX::Abstract::Vector::length() const =0

Return the length of vector.

The length of this vector

Even if the vector is distributed across processors, this should
return the  global length of the vector. ";

%feature("docstring")  NOX::Abstract::Vector::print "void
NOX::Abstract::Vector::print(std::ostream &stream) const

Print the vector. To be used for debugging only. ";


// File: classNOX_1_1Epetra_1_1VectorSpace.xml
%feature("docstring") NOX::Epetra::VectorSpace "

Pure virtual base class for the vector space used by
NOX::Epetra::Vectors.

This class allows users to override the inner product and norm used by
the NOX::Epetra::Vector class. The most frequent use of this class is
for introducing a weighted norm throughout NOX.

C++ includes: NOX_Epetra_VectorSpace.H ";

%feature("docstring")  NOX::Epetra::VectorSpace::VectorSpace "NOX::Epetra::VectorSpace::VectorSpace()

Constructor. ";

%feature("docstring")  NOX::Epetra::VectorSpace::~VectorSpace "virtual NOX::Epetra::VectorSpace::~VectorSpace()

Destructor. ";

%feature("docstring")  NOX::Epetra::VectorSpace::innerProduct "virtual double NOX::Epetra::VectorSpace::innerProduct(const
Epetra_Vector &a, const Epetra_Vector &b) const =0

Computes the inner product: <a,b>. ";

%feature("docstring")  NOX::Epetra::VectorSpace::norm "virtual double
NOX::Epetra::VectorSpace::norm(const Epetra_Vector &a,
NOX::Abstract::Vector::NormType=NOX::Abstract::Vector::TwoNorm) const
=0

Computes the norm.

For an L2 norm, the computation is: sqrt( <a,a> ). ";


// File: classNOX_1_1Epetra_1_1VectorSpaceL2.xml
%feature("docstring") NOX::Epetra::VectorSpaceL2 "

Concrete class for an L2 vector space.

This class allows users to override the inner product and norm used by
the NOX::Epetra::Vector class. The most frequent use of this class is
for introducing a weighted norm throughout NOX.

C++ includes: NOX_Epetra_VectorSpace_L2.H ";

%feature("docstring")  NOX::Epetra::VectorSpaceL2::VectorSpaceL2 "NOX::Epetra::VectorSpaceL2::VectorSpaceL2()

Constructor. ";

%feature("docstring")  NOX::Epetra::VectorSpaceL2::~VectorSpaceL2 "NOX::Epetra::VectorSpaceL2::~VectorSpaceL2()

Destructor. ";

%feature("docstring")  NOX::Epetra::VectorSpaceL2::innerProduct "double NOX::Epetra::VectorSpaceL2::innerProduct(const Epetra_Vector
&a, const Epetra_Vector &b) const

Computes the inner product: <a,b>. ";

%feature("docstring")  NOX::Epetra::VectorSpaceL2::norm "double
NOX::Epetra::VectorSpaceL2::norm(const Epetra_Vector &a,
NOX::Abstract::Vector::NormType=NOX::Abstract::Vector::TwoNorm) const

Computes the norm.

For an L2 norm, the computation is: sqrt( <a,a> ). ";


// File: classNOX_1_1Epetra_1_1VectorSpaceScaledL2.xml
%feature("docstring") NOX::Epetra::VectorSpaceScaledL2 "

Concrete class for a weighted L2 vector space.

This class allows users to override the inner product and norm used by
the NOX::Epetra::Vector class. The most frequent use of this class is
for introducing a weighted norm throughout NOX.

C++ includes: NOX_Epetra_VectorSpace_ScaledL2.H ";

%feature("docstring")
NOX::Epetra::VectorSpaceScaledL2::VectorSpaceScaledL2 "NOX::Epetra::VectorSpaceScaledL2::VectorSpaceScaledL2(const
Teuchos::RCP< NOX::Epetra::Scaling > &s,
NOX::Epetra::Scaling::ScaleType st=NOX::Epetra::Scaling::Left)

Constructor. ";

%feature("docstring")
NOX::Epetra::VectorSpaceScaledL2::~VectorSpaceScaledL2 "NOX::Epetra::VectorSpaceScaledL2::~VectorSpaceScaledL2()

Destructor. ";

%feature("docstring")  NOX::Epetra::VectorSpaceScaledL2::innerProduct
"double NOX::Epetra::VectorSpaceScaledL2::innerProduct(const
Epetra_Vector &a, const Epetra_Vector &b) const

Computes a scaled inner product.

Computes a scaled inner product: $ <Da, Db> $ where $D$ is the set of
scaling vectors associated with either left of right scaling. ";

%feature("docstring")  NOX::Epetra::VectorSpaceScaledL2::norm "double
NOX::Epetra::VectorSpaceScaledL2::norm(const Epetra_Vector &a,
NOX::Abstract::Vector::NormType=NOX::Abstract::Vector::TwoNorm) const

Computes the scaled norm.

Computes the scaled norm using $ Da $ where $D$ is the set of scaling
vectors associated with either left of right scaling. ";


// File: namespaceNOX.xml
%feature("docstring")  NOX::Abstract::version "string NOX::version()

Returns a string with the current version number of the NOX code. ";


// File: namespaceNOX_1_1Abstract.xml


// File: namespaceNOX_1_1Direction.xml
%feature("docstring")  NOX::Direction::Utils::buildDirection "Teuchos::RCP<NOX::Direction::Generic>
NOX::Direction::buildDirection(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params) ";


// File: namespaceNOX_1_1Direction_1_1Utils.xml


// File: namespaceNOX_1_1Epetra.xml


// File: namespaceNOX_1_1Epetra_1_1Interface.xml


// File: namespaceNOX_1_1LineSearch.xml
%feature("docstring")  NOX::LineSearch::Utils::buildLineSearch "Teuchos::RCP<NOX::LineSearch::Generic>
NOX::LineSearch::buildLineSearch(const Teuchos::RCP< NOX::GlobalData >
&gd, Teuchos::ParameterList &params) ";


// File: namespaceNOX_1_1LineSearch_1_1Utils.xml


// File: namespaceNOX_1_1MeritFunction.xml


// File: namespaceNOX_1_1Multiphysics.xml


// File: namespaceNOX_1_1Multiphysics_1_1DataExchange.xml


// File: namespaceNOX_1_1Multiphysics_1_1Solver.xml


// File: namespaceNOX_1_1Parameter.xml


// File: namespaceNOX_1_1Solver.xml
%feature("docstring")  NOX::Solver::buildSolver "Teuchos::RCP<NOX::Solver::Generic> NOX::Solver::buildSolver(const
Teuchos::RCP< NOX::Abstract::Group > &grp, const Teuchos::RCP<
NOX::StatusTest::Generic > &tests, const Teuchos::RCP<
Teuchos::ParameterList > &params) ";

%feature("docstring")  NOX::Solver::parseStatusTestCheckType "NOX::StatusTest::CheckType
NOX::Solver::parseStatusTestCheckType(Teuchos::ParameterList
&solver_options_list)

Nonmember method that returns the status test check type.

This object parses the \"Solver Options\" parameter list for a key
\"Status Test Check Type\" of type <std::string> with possible values:

\"Complete\"

\"Minimal\" (default)

\"None\"

These options correspond to the NOX::StatusTest::CheckType. Please
follow the link for this object for more information. ";


// File: namespaceNOX_1_1StatusTest.xml
%feature("docstring")  NOX::StatusTest::buildStatusTests "Teuchos::RCP<NOX::StatusTest::Generic>
NOX::StatusTest::buildStatusTests(const std::string &file_name, const
NOX::Utils &utils, std::map< std::string, Teuchos::RCP<
NOX::StatusTest::Generic > > *tagged_tests) ";

%feature("docstring")  NOX::StatusTest::buildStatusTests "Teuchos::RCP<NOX::StatusTest::Generic>
NOX::StatusTest::buildStatusTests(Teuchos::ParameterList &p, const
NOX::Utils &utils, std::map< std::string, Teuchos::RCP<
NOX::StatusTest::Generic > > *tagged_tests) ";


// File: namespacestd.xml


// File: namespaceTeuchos.xml


// File: NOX_8H.xml


// File: NOX__Abstract__Group_8C.xml


// File: NOX__Abstract__Group_8H.xml


// File: NOX__Abstract__MultiVector_8H.xml


// File: NOX__Abstract__PrePostOperator_8H.xml


// File: NOX__Abstract__Vector_8C.xml


// File: NOX__Abstract__Vector_8H.xml


// File: NOX__Common_8H.xml


// File: NOX__Description_8H.xml


// File: NOX__Direction__Broyden_8C.xml


// File: NOX__Direction__Broyden_8H.xml


// File: NOX__Direction__Factory_8C.xml


// File: NOX__Direction__Factory_8H.xml


// File: NOX__Direction__Generic_8C.xml


// File: NOX__Direction__Generic_8H.xml


// File: NOX__Direction__ModifiedNewton_8C.xml


// File: NOX__Direction__ModifiedNewton_8H.xml


// File: NOX__Direction__Newton_8C.xml


// File: NOX__Direction__Newton_8H.xml


// File: NOX__Direction__NonlinearCG_8C.xml


// File: NOX__Direction__NonlinearCG_8H.xml


// File: NOX__Direction__QuasiNewton_8C.xml


// File: NOX__Direction__QuasiNewton_8H.xml


// File: NOX__Direction__SteepestDescent_8C.xml


// File: NOX__Direction__SteepestDescent_8H.xml


// File: NOX__Direction__Tensor_8C.xml


// File: NOX__Direction__Tensor_8H.xml


// File: NOX__Direction__UserDefinedFactory_8H.xml


// File: NOX__Direction__UserDefinedFactoryT_8H.xml


// File: NOX__Direction__Utils__InexactNewton_8C.xml


// File: NOX__Direction__Utils__InexactNewton_8H.xml


// File: NOX__Epetra_8H.xml


// File: NOX__Epetra__BroydenOperator_8C.xml


// File: NOX__Epetra__BroydenOperator_8H.xml


// File: NOX__Epetra__FiniteDifference_8C.xml


// File: NOX__Epetra__FiniteDifference_8H.xml


// File: NOX__Epetra__FiniteDifferenceColoring_8C.xml


// File: NOX__Epetra__FiniteDifferenceColoring_8H.xml


// File: NOX__Epetra__FiniteDifferenceColoringWithUpdate_8C.xml


// File: NOX__Epetra__FiniteDifferenceColoringWithUpdate_8H.xml


// File: NOX__Epetra__Group_8C.xml


// File: NOX__Epetra__Group_8H.xml


// File: NOX__Epetra__Interface__Jacobian_8H.xml


// File: NOX__Epetra__Interface__Preconditioner_8H.xml


// File: NOX__Epetra__Interface__Required_8H.xml


// File: NOX__Epetra__LinearSystem_8H.xml


// File: NOX__Epetra__LinearSystem__Amesos_8C.xml


// File: NOX__Epetra__LinearSystem__Amesos_8H.xml


// File: NOX__Epetra__LinearSystem__AztecOO_8C.xml


// File: NOX__Epetra__LinearSystem__AztecOO_8H.xml


// File: NOX__Epetra__LinearSystem__Stratimikos_8C.xml


// File: NOX__Epetra__LinearSystem__Stratimikos_8H.xml


// File: NOX__Epetra__MatrixFree_8C.xml


// File: NOX__Epetra__MatrixFree_8H.xml


// File: NOX__Epetra__ModelEvaluatorInterface_8C.xml


// File: NOX__Epetra__ModelEvaluatorInterface_8H.xml


// File: NOX__Epetra__MultiVector_8C.xml


// File: NOX__Epetra__MultiVector_8H.xml


// File: NOX__Epetra__Observer_8H.xml


// File: NOX__Epetra__Scaling_8C.xml


// File: NOX__Epetra__Scaling_8H.xml


// File: NOX__Epetra__Vector_8C.xml


// File: NOX__Epetra__Vector_8H.xml


// File: NOX__Epetra__VectorSpace_8H.xml


// File: NOX__Epetra__VectorSpace__L2_8C.xml


// File: NOX__Epetra__VectorSpace__L2_8H.xml


// File: NOX__Epetra__VectorSpace__ScaledL2_8C.xml


// File: NOX__Epetra__VectorSpace__ScaledL2_8H.xml


// File: NOX__GlobalData_8C.xml


// File: NOX__GlobalData_8H.xml


// File: NOX__LineSearch__Backtrack_8C.xml


// File: NOX__LineSearch__Backtrack_8H.xml


// File: NOX__LineSearch__Factory_8C.xml


// File: NOX__LineSearch__Factory_8H.xml


// File: NOX__LineSearch__FullStep_8C.xml


// File: NOX__LineSearch__FullStep_8H.xml


// File: NOX__LineSearch__Generic_8H.xml


// File: NOX__LineSearch__MoreThuente_8C.xml


// File: NOX__LineSearch__MoreThuente_8H.xml


// File: NOX__LineSearch__NonlinearCG_8C.xml


// File: NOX__LineSearch__NonlinearCG_8H.xml


// File: NOX__LineSearch__Polynomial_8C.xml


// File: NOX__LineSearch__Polynomial_8H.xml


// File: NOX__LineSearch__Tensor_8C.xml


// File: NOX__LineSearch__Tensor_8H.xml


// File: NOX__LineSearch__UserDefinedFactory_8H.xml


// File: NOX__LineSearch__UserDefinedFactoryT_8H.xml


// File: NOX__LineSearch__Utils__Counters_8C.xml


// File: NOX__LineSearch__Utils__Counters_8H.xml


// File: NOX__LineSearch__Utils__Printing_8C.xml


// File: NOX__LineSearch__Utils__Printing_8H.xml


// File: NOX__LineSearch__Utils__Slope_8C.xml


// File: NOX__LineSearch__Utils__Slope_8H.xml


// File: NOX__MeritFunction__Generic_8H.xml


// File: NOX__MeritFunction__SumOfSquares_8C.xml


// File: NOX__MeritFunction__SumOfSquares_8H.xml


// File: NOX__Multiphysics__DataExchange__Interface_8H.xml


// File: NOX__Multiphysics__Group_8C.xml


// File: NOX__Multiphysics__Group_8H.xml


// File: NOX__Multiphysics__Solver__FixedPointBased_8C.xml


// File: NOX__Multiphysics__Solver__FixedPointBased_8H.xml


// File: NOX__Multiphysics__Solver__Generic_8H.xml


// File: NOX__Multiphysics__Solver__Manager_8C.xml


// File: NOX__Multiphysics__Solver__Manager_8H.xml


// File: NOX__MultiVector_8C.xml


// File: NOX__MultiVector_8H.xml


// File: NOX__Random_8C.xml


// File: NOX__Random_8H.xml


// File: NOX__SharedObjectTemplate_8H.xml


// File: NOX__Solver__Factory_8C.xml


// File: NOX__Solver__Factory_8H.xml


// File: NOX__Solver__Generic_8H.xml


// File: NOX__Solver__InexactTrustRegionBased_8C.xml


// File: NOX__Solver__InexactTrustRegionBased_8H.xml


// File: NOX__Solver__LineSearchBased_8C.xml


// File: NOX__Solver__LineSearchBased_8H.xml


// File: NOX__Solver__PrePostOperator_8C.xml


// File: NOX__Solver__PrePostOperator_8H.xml


// File: NOX__Solver__SolverUtils_8C.xml


// File: NOX__Solver__SolverUtils_8H.xml


// File: NOX__Solver__TensorBased_8C.xml


// File: NOX__Solver__TensorBased_8H.xml


// File: NOX__Solver__TensorBasedTest_8C.xml


// File: NOX__Solver__TensorBasedTest_8H.xml


// File: NOX__Solver__TrustRegionBased_8C.xml


// File: NOX__Solver__TrustRegionBased_8H.xml


// File: NOX__StatusTest__Combo_8C.xml


// File: NOX__StatusTest__Combo_8H.xml


// File: NOX__StatusTest__Divergence_8C.xml


// File: NOX__StatusTest__Divergence_8H.xml


// File: NOX__StatusTest__Factory_8C.xml


// File: NOX__StatusTest__Factory_8H.xml


// File: NOX__StatusTest__FiniteValue_8C.xml


// File: NOX__StatusTest__FiniteValue_8H.xml


// File: NOX__StatusTest__Generic_8C.xml


// File: NOX__StatusTest__Generic_8H.xml


// File: NOX__StatusTest__MaxIters_8C.xml


// File: NOX__StatusTest__MaxIters_8H.xml


// File: NOX__StatusTest__NormF_8C.xml


// File: NOX__StatusTest__NormF_8H.xml


// File: NOX__StatusTest__NormUpdate_8C.xml


// File: NOX__StatusTest__NormUpdate_8H.xml


// File: NOX__StatusTest__NormWRMS_8C.xml


// File: NOX__StatusTest__NormWRMS_8H.xml


// File: NOX__StatusTest__Stagnation_8C.xml


// File: NOX__StatusTest__Stagnation_8H.xml


// File: NOX__Utils_8C.xml


// File: NOX__Utils_8H.xml


// File: NOX__Version_8C.xml


// File: NOX__Version_8H.xml


// File: nox_user_information.xml


// File: step1.xml


// File: step2.xml


// File: step3.xml


// File: step4.xml


// File: nox_developer_information.xml


// File: nox_release_information.xml


// File: nox_configuration_options.xml


// File: prerelease.xml


// File: cvsrepos.xml


// File: coding.xml


// File: nox_class_overview.xml


// File: parameters.xml


// File: epetra_interface.xml


// File: thyra_interface.xml


// File: petsc_interface.xml


// File: nox_epetra_tutorial.xml


// File: portability_issues.xml


// File: deprecated.xml


// File: dir_b91ee76c155b82a02721c1190b219a22.xml


// File: dir_2d461d9e99e31007a9f2d6e54cc9f896.xml


// File: dir_508aeab7ec1c3def73e9e6be4a4635da.xml


// File: indexpage.xml


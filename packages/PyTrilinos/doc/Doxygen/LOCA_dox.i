
// File: index.xml

// File: classLOCA_1_1StatusTest_1_1Abstract.xml
%feature("docstring") LOCA::StatusTest::Abstract "";

%feature("docstring")  LOCA::StatusTest::Abstract::Abstract "LOCA::StatusTest::Abstract::Abstract()

Constructor. ";

%feature("docstring")  LOCA::StatusTest::Abstract::~Abstract "virtual
LOCA::StatusTest::Abstract::~Abstract()

Destructor. ";

%feature("docstring")  LOCA::StatusTest::Abstract::checkStatus "virtual LOCA::StatusTest::StatusType
LOCA::StatusTest::Abstract::checkStatus(const LOCA::Abstract::Iterator
&stepper, LOCA::StatusTest::CheckType checkType)=0

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
LOCA::StatusType::None. If the test is skipped, then the status should
be set to LOCA::StatusTest::Unevaluated. ";

%feature("docstring")  LOCA::StatusTest::Abstract::getStatus "virtual
LOCA::StatusTest::StatusType LOCA::StatusTest::Abstract::getStatus()
const =0

Return the result of the most recent checkStatus call. ";

%feature("docstring")  LOCA::StatusTest::Abstract::print "virtual
std::ostream& LOCA::StatusTest::Abstract::print(std::ostream &stream,
int indent=0) const =0

Output formatted description of stopping test to output stream. ";


// File: classLOCA_1_1Parameter_1_1AbstractEntry.xml
%feature("docstring") LOCA::Parameter::AbstractEntry "

Abstract interface for all entries in LOCA::Parameter::Library.

This class doesn't have much of an interface and really only serves
the purpose of having a common parent class for parameter entries of
all value types.

C++ includes: LOCA_Parameter_Entry.H ";

%feature("docstring")  LOCA::Parameter::AbstractEntry::AbstractEntry "LOCA::Parameter::AbstractEntry::AbstractEntry()

Default contructor. ";

%feature("docstring")  LOCA::Parameter::AbstractEntry::~AbstractEntry
"virtual LOCA::Parameter::AbstractEntry::~AbstractEntry()

Destructor. ";


// File: classLOCA_1_1BorderedSystem_1_1AbstractGroup.xml
%feature("docstring") LOCA::BorderedSystem::AbstractGroup "

An interface for groups that are bordered systems.

This class provides an interface for groups whose Jacobian is of the
form \\\\[ \\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\] where $A$ and $B$ are multivectors and $C$ is a dense matrix. It
provides methods for determining the width of the bordered
rows/columns and for extracting these components. It is intended to be
a recusive interface, in that if the group representing $J$ also has
this form, the extracted rows/columns should be the combined
rows/columns of this group and the underlying group (and so on).

C++ includes: LOCA_BorderedSystem_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::getBorderedWidth "virtual int
LOCA::BorderedSystem::AbstractGroup::getBorderedWidth() const =0

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::getUnborderedGroup "virtual
Teuchos::RCP<const NOX::Abstract::Group>
LOCA::BorderedSystem::AbstractGroup::getUnborderedGroup() const =0

Get bottom-level unbordered group. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::isCombinedAZero "virtual bool
LOCA::BorderedSystem::AbstractGroup::isCombinedAZero() const =0

Indicates whether combined A block is zero. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::isCombinedBZero "virtual bool
LOCA::BorderedSystem::AbstractGroup::isCombinedBZero() const =0

Indicates whether combined B block is zero. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::isCombinedCZero "virtual bool
LOCA::BorderedSystem::AbstractGroup::isCombinedCZero() const =0

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::extractSolutionComponent "virtual void
LOCA::BorderedSystem::AbstractGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const
=0

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::extractParameterComponent "virtual void
LOCA::BorderedSystem::AbstractGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const =0

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::loadNestedComponents "virtual
void LOCA::BorderedSystem::AbstractGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const =0

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")  LOCA::BorderedSystem::AbstractGroup::fillA "virtual void
LOCA::BorderedSystem::AbstractGroup::fillA(NOX::Abstract::MultiVector
&A) const =0

Fill the combined A block as described above. ";

%feature("docstring")  LOCA::BorderedSystem::AbstractGroup::fillB "virtual void
LOCA::BorderedSystem::AbstractGroup::fillB(NOX::Abstract::MultiVector
&B) const =0

Fill the combined B block as described above. ";

%feature("docstring")  LOCA::BorderedSystem::AbstractGroup::fillC "virtual void
LOCA::BorderedSystem::AbstractGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const =0

Fill the combined C block as described above. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::AbstractGroup "LOCA::BorderedSystem::AbstractGroup::AbstractGroup()

Constructor. ";

%feature("docstring")
LOCA::BorderedSystem::AbstractGroup::~AbstractGroup "virtual
LOCA::BorderedSystem::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1TimeDependent_1_1AbstractGroup.xml
%feature("docstring") LOCA::TimeDependent::AbstractGroup "

Interface to underlying groups for time dependent systems.

This abstract class provides an interface for time dependent problems,
i.e., problems with a mass matrix (typically used in eignvalue or Hopf
calculations). It provides pure virtual methods for computing and
manipulating the shifted matrix $\\\\alpha J + \\\\beta M$ where $J$
is the Jacobian matrix and $M$ is the mass matrix.

C++ includes: LOCA_TimeDependent_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::computeShiftedMatrix "virtual
NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::computeShiftedMatrix(double alpha,
double beta)=0

Compute the shifted matrix. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::applyShiftedMatrix "virtual
NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::applyShiftedMatrix(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const =0

Multiply the shifted matrix by a vector. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::applyShiftedMatrixMultiVector "virtual NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::applyShiftedMatrixMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const =0

Multiply the shifted matrix by a multi-vector. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::applyShiftedMatrixInverseMultiVector
"virtual NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::applyShiftedMatrixInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const =0

Apply the inverse of the shifted matrix by a multi-vector, as needed
by the shift-and-invert and generalized Cayley transformations. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::AbstractGroup "LOCA::TimeDependent::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::~AbstractGroup "virtual
LOCA::TimeDependent::AbstractGroup::~AbstractGroup()

Destructor. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::computeSecondShiftedMatrix "virtual NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::computeSecondShiftedMatrix(double
alpha, double beta)=0

Compute the second shifted matrix. Can avoid recomputing if two are
stored.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::applySecondShiftedMatrix "virtual
NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::applySecondShiftedMatrix(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const =0

Multiply the shifted matrix by a vector.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::TimeDependent::AbstractGroup::applySecondShiftedMatrixMultiVector
"virtual NOX::Abstract::Group::ReturnType
LOCA::TimeDependent::AbstractGroup::applySecondShiftedMatrixMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const =0

Multiply the shifted matrix by a multi-vector.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";


// File: classLOCA_1_1TurningPoint_1_1MinimallyAugmented_1_1AbstractGroup.xml
%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup "

Interface to underlying groups for turning point calculations using
the minimally augmented formulation.

This abstract class provides the required interface for underlying
groups to locate turning points using the minimally augmented turning
point formulation (see
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup for a
description of the governing equations).

This class is derived from the
LOCA::TurningPoint::MooreSpence::AbstractGroup and declares several
pure virtual methods compute various derivatives of $w^TJn$ for a
given $w$ and $n$. Default implementations for the derivatives using
finite differencing are implemented in the
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup.

C++ includes: LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::computeDwtJnDp
"virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::computeDwtJnDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &w, const
NOX::Abstract::Vector &nullVector,
NOX::Abstract::MultiVector::DenseMatrix &result, bool isValid)=0

Computes the derivative $\\\\partial w^TJn/\\\\partial p$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::computeDwtJDp "virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::computeDwtJDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &w,
NOX::Abstract::MultiVector &result, bool isValid)=0

Computes the derivative $\\\\partial w^TJ/\\\\partial p$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::computeDwtJnDx
"virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::computeDwtJnDx(const
NOX::Abstract::Vector &w, const NOX::Abstract::Vector &nullVector,
NOX::Abstract::Vector &result)=0

Computes the derivative $\\\\frac{\\\\partial w^TJn}{\\\\partial x}$.
";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::AbstractGroup "LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::~AbstractGroup
"virtual
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1AbstractGroup.xml
%feature("docstring") LOCA::TurningPoint::MooreSpence::AbstractGroup "

Interface to underlying groups for turning point calculations using
the Moore-Spence formulation.

This abstract class provides the required interface for underlying
groups to locate turning points using the bordering algorithm for the
Moore-Spence turning point formulation (see
LOCA::TurningPoint::MooreSpence::ExtendedGroup for a description of
the governing equations).

This class is derived from the LOCA::MultiContinuation::AbstractGroup
and declares several pure virtual methods compute various derivatives
of $Jn$ for a given $n$. Default implementations for the derivatives
using finite differencing are implemented in the
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup.

C++ includes: LOCA_TurningPoint_MooreSpence_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDJnDpMulti "virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDJnDpMulti(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &nullVector,
NOX::Abstract::MultiVector &result, bool isValid)=0

Computes the derivative $\\\\partial Jn/\\\\partial p$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDJnDxaMulti "virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDJnDxaMulti(const
NOX::Abstract::Vector &nullVector, const NOX::Abstract::MultiVector
&aVector, NOX::Abstract::MultiVector &result)=0

Computes the directional derivative $\\\\frac{\\\\partial
Jn}{\\\\partial x} a$ for the given direction $a$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDJnDxaMulti "virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDJnDxaMulti(const
NOX::Abstract::Vector &nullVector, const NOX::Abstract::Vector
&JnVector, const NOX::Abstract::MultiVector &aVector,
NOX::Abstract::MultiVector &result)=0

Computes the directional derivative $\\\\frac{\\\\partial
Jn}{\\\\partial x} a$ for the given direction $a$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDwtJnDxMulti "virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::AbstractGroup::computeDwtJnDxMulti(const
NOX::Abstract::MultiVector &w, const NOX::Abstract::Vector
&nullVector, NOX::Abstract::MultiVector &result)=0

Computes the derivative $\\\\frac{\\\\partial w^TJn}{\\\\partial x}$.
";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::AbstractGroup::AbstractGroup "LOCA::TurningPoint::MooreSpence::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::AbstractGroup::~AbstractGroup "virtual
LOCA::TurningPoint::MooreSpence::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1MultiContinuation_1_1AbstractGroup.xml
%feature("docstring") LOCA::MultiContinuation::AbstractGroup "

LOCA abstract interface for continuation, derived from the
NOX::Abstract::Group. This abstract class provides the interface
necessary to perform continuation, i.e., compute families of solutions
to $ F(x,p) = 0 $.

Concrete implemenations of this interface must provide implementations
of all of the methods in the NOX::Abstract::Group interface as well as
the additional interface defined here.

C++ includes: LOCA_MultiContinuation_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")  LOCA::MultiContinuation::AbstractGroup::copy "virtual void LOCA::MultiContinuation::AbstractGroup::copy(const
NOX::Abstract::Group &source)=0

Copy the group (replaces operator = ) ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::setParamsMulti "virtual void
LOCA::MultiContinuation::AbstractGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)=0

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::setParams "virtual void
LOCA::MultiContinuation::AbstractGroup::setParams(const
LOCA::ParameterVector &p)=0

Set the parameter vector in the group to p (pVector = p). ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::setParam "virtual void
LOCA::MultiContinuation::AbstractGroup::setParam(int paramID, double
val)=0

Set parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::setParam "virtual void
LOCA::MultiContinuation::AbstractGroup::setParam(std::string paramID,
double val)=0

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::getParams "virtual const
LOCA::ParameterVector&
LOCA::MultiContinuation::AbstractGroup::getParams() const =0

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::getParam "virtual double
LOCA::MultiContinuation::AbstractGroup::getParam(int paramID) const =0

Return copy of parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::getParam "virtual double
LOCA::MultiContinuation::AbstractGroup::getParam(std::string paramID)
const =0

Return copy of parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::computeDfDpMulti "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::AbstractGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)=0

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

/*  Virtual methods with default implementations  */

/* These methods should be overloaded in a concrete implementation if
more appropriate/efficient approaches are available.

*/

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::preProcessContinuationStep "void
LOCA::MultiContinuation::AbstractGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. The default implementation to empty. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::postProcessContinuationStep "void
LOCA::MultiContinuation::AbstractGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. The
default implementation to empty. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::projectToDraw "void
LOCA::MultiContinuation::AbstractGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation.

This method is called every time a solution is saved by the
multiparameter continuation code MF for later visualization and should
project the solution vector down to a few scalars. The array px will
be preallocated to the proper length given by
projectToDrawDimension().

The default implementation is the max norm of the vector. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::projectToDrawDimension "int
LOCA::MultiContinuation::AbstractGroup::projectToDrawDimension() const

Returns the dimension of the project to draw array.

The default implementation is to return 1 since the default projection
is the max norm of the vector (a scalar). ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::computeScaledDotProduct "double
LOCA::MultiContinuation::AbstractGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute a scaled dot product.

The default implementation here just computes a.dot(b) but should be
overloaded for any problem that his difficult scaling. ";

/*  Virtual methods with empty or trivial implementations  */

/* These methods should be overloaded in a concrete implementation but
their implementation is not critical to the rest of LOCA and therefore
have empty or trivial implementations.

*/

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::printSolution "virtual void
LOCA::MultiContinuation::AbstractGroup::printSolution(const double
conParam) const

Function to print out solution and parameter after successful step.

Empty default definition. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::printSolution "virtual void
LOCA::MultiContinuation::AbstractGroup::printSolution(const
NOX::Abstract::Vector &x_, const double conParam) const

Function to print out a vector and parameter after successful step.

Empty default definition. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::scaleVector "void
LOCA::MultiContinuation::AbstractGroup::scaleVector(NOX::Abstract::Vector
&x) const

Scales a vector using scaling vector.

The default definition here is to do nothing, i.e., no scaling ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::AbstractGroup "LOCA::MultiContinuation::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractGroup::~AbstractGroup "virtual
LOCA::MultiContinuation::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1PhaseTransition_1_1AbstractGroup.xml
%feature("docstring") LOCA::PhaseTransition::AbstractGroup "

Interface to underlying groups for phase transition calculations.

This abstract class provides the required interface for underlying
groups to locate phase transitions using the bordering algorithm from
the Salinger&Frink (2003) paper.

This class is derived from the LOCA::MultiContinuation::AbstractGroup
and declares a pure virtual method for computing the free energy.

C++ includes: LOCA_PhaseTransition_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::PhaseTransition::AbstractGroup::computeFreeEnergy "virtual
double LOCA::PhaseTransition::AbstractGroup::computeFreeEnergy()=0

Computes the free energy at the current solution and parameter values.
";

%feature("docstring")
LOCA::PhaseTransition::AbstractGroup::AbstractGroup "LOCA::PhaseTransition::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::PhaseTransition::AbstractGroup::~AbstractGroup "virtual
LOCA::PhaseTransition::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1Pitchfork_1_1MinimallyAugmented_1_1AbstractGroup.xml
%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup "

Interface to underlying groups for pitchfork calculations using the
minimally augmented formulation.

This abstract class provides the required interface for underlying
groups to locate pitchforks using the minimally augmented pitchfork
formulation (see LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup
for a description of the governing equations).

This class is derived from the
LOCA::Pitchfork::MooreSpence::AbstractGroup and
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup and does not
declare any new virtual methods.

C++ includes: LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup::AbstractGroup "LOCA::Pitchfork::MinimallyAugmented::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup::~AbstractGroup "virtual
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1AbstractGroup.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::AbstractGroup "

Interface to underlying groups for pitchfork calculations using the
Moore-Spence formulation.

This abstract class provides the required interface for underlying
groups to locate pitchforks using the bordering algorithms for the
Moore-Spence pitchfork formulation (see
LOCA::Pitchfork::MooreSpence::ExtendedGroup for a description of the
governing equations).

This class is derived from the
LOCA::TurningPoint::MooreSpence::AbstractGroup and declares a single
virtual method, innerProduct(), to compute the inner product of a
vector with the asymmetry vector. It has a default implementation
given by the dot product, but should be overloaded for any problem
that has a different definition for the inner product.

C++ includes: LOCA_Pitchfork_MooreSpence_AbstractGroup.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::AbstractGroup::AbstractGroup "LOCA::Pitchfork::MooreSpence::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::AbstractGroup::~AbstractGroup "virtual
LOCA::Pitchfork::MooreSpence::AbstractGroup::~AbstractGroup()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::AbstractGroup::innerProduct "virtual
double LOCA::Pitchfork::MooreSpence::AbstractGroup::innerProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute the inner product of a and b.

The default implementation is given by the dot product of a and b. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::AbstractGroup::innerProduct "virtual
void LOCA::Pitchfork::MooreSpence::AbstractGroup::innerProduct(const
NOX::Abstract::MultiVector &a, const NOX::Abstract::MultiVector &b,
NOX::Abstract::MultiVector::DenseMatrix &c) const

Compute the inner product of a and b.

The default implementation is given by the dot product of a and b. ";


// File: classLOCA_1_1Homotopy_1_1AbstractGroup.xml
%feature("docstring") LOCA::Homotopy::AbstractGroup "

Interface to underlying groups for homotopy calculations.

This abstract class provides an interface for a homotopy technique for
solving nonlinear equations. See LOCA::Homotopy::Group for a
description of the technique used. This class provides a single pure
virtual method, augmentJacobianForHomotopy(), which scales the
diagonal of the Jacobian by a constant times the identity matrix.

C++ includes: LOCA_Homotopy_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::Homotopy::AbstractGroup::augmentJacobianForHomotopy "virtual
NOX::Abstract::Group::ReturnType
LOCA::Homotopy::AbstractGroup::augmentJacobianForHomotopy(double a,
double b)=0

Replace Jacobian $J$ by $aJ+bI$ where $I$ is the identity matrix. ";

%feature("docstring")  LOCA::Homotopy::AbstractGroup::AbstractGroup "LOCA::Homotopy::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")  LOCA::Homotopy::AbstractGroup::~AbstractGroup "virtual LOCA::Homotopy::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1Hopf_1_1MinimallyAugmented_1_1AbstractGroup.xml
%feature("docstring") LOCA::Hopf::MinimallyAugmented::AbstractGroup "

Interface to underlying groups for Hopf calculations using the
minimally augmented formulation.

This abstract class provides the required interface for underlying
groups to locate Hopfs using the minimally augmented Hopf (see
LOCA::Hopf::MinimallyAugmented::ExtendedGroup for a description of the
governing equations).

This class is derived from the LOCA::Hopf::MooreSpence::AbstractGroup
and declares several pure virtual methods to compute various
derivatives of $w^H C e$ for a given $w = w_1 + i w_2$ and $e = y + i
z$ where $C = J + i \\\\omega M$. Default implementations for the
derivatives using finite differencing are implemented in the
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup.

C++ includes: LOCA_Hopf_MinimallyAugmented_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::applyComplexTranspose "virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::AbstractGroup::applyComplexTranspose(const
NOX::Abstract::Vector &input_real, const NOX::Abstract::Vector
&input_imag, NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag) const =0

Computes conjugate-tranpose matrix vector product $ (J+i\\\\omega M)^H
(x + iy) $. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::applyComplexTransposeMultiVector
"virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::AbstractGroup::applyComplexTransposeMultiVector(const
NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const =0

Computes conjugate-tranpose matrix vector product $ (J+i\\\\omega M)^H
(x + iy) $. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::applyComplexTransposeInverseMultiVector
"virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::AbstractGroup::applyComplexTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const =0

Solve $(J+i\\\\omega M)^H (x + iy) = a+ib$. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::computeDwtCeDp "virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::AbstractGroup::computeDwtCeDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &w1, const
NOX::Abstract::Vector &w2, const NOX::Abstract::Vector &y, const
NOX::Abstract::Vector &z, double omega,
NOX::Abstract::MultiVector::DenseMatrix &result_real,
NOX::Abstract::MultiVector::DenseMatrix &result_imag, bool isValid)=0

Computes the derivative $\\\\partial w^H C e/\\\\partial p$, $w = w_1
+ i w_2$, $e = y + i z$, $C = J + i \\\\omega M$. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::computeDwtCeDx "virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::AbstractGroup::computeDwtCeDx(const
NOX::Abstract::Vector &w1, const NOX::Abstract::Vector &w2, const
NOX::Abstract::Vector &y, const NOX::Abstract::Vector &z, double
omega, NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag)=0

Computes the derivative $\\\\frac{\\\\partial w^H C e}{\\\\partial
x}$, $w = w_1 + i w_2$, $e = y + i z$, $C = J + i \\\\omega M$. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::AbstractGroup "LOCA::Hopf::MinimallyAugmented::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::AbstractGroup::~AbstractGroup "virtual
LOCA::Hopf::MinimallyAugmented::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1AbstractGroup.xml
%feature("docstring") LOCA::Hopf::MooreSpence::AbstractGroup "

Interface to underlying groups for Hopf point calculations using the
Moore-Spence formulation.

This abstract class provides the required interface for underlying
groups to locate Hopf bifurcations using the bordering algorithm for
the Moore-Spence Hopf fomulation (see
LOCA::Hopf::MooreSpence::ExtendedGroup for a description of the
governing equations).

This class is derived from the
LOCA::TurningPoint::MooreSpence::AbstractGroup and
LOCa::TimeDependent::AbstractGroup and declares several pure virtual
methods for applying, solving, and computing derivatives of the
complex matrix $J+i\\\\omega B$ where $J$ is the Jacobian matrix, $B$
is the mass matrix, and $\\\\omega$ is a (real) scalar.

C++ includes: LOCA_Hopf_MooreSpence_AbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::isComplex "virtual bool
LOCA::Hopf::MooreSpence::AbstractGroup::isComplex() const =0

Is $J+i\\\\omega B$ valid. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::computeComplex "virtual
NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::computeComplex(double
frequency)=0

Compute $J+i\\\\omega B$.

The argument frequency stores $\\\\omega$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::applyComplex "virtual
NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::applyComplex(const
NOX::Abstract::Vector &input_real, const NOX::Abstract::Vector
&input_imag, NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag) const =0

Compute $(J+i\\\\omega B)(y+iz)$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::applyComplexMultiVector "virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::applyComplexMultiVector(const
NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const =0

Compute $(J+i\\\\omega B)(y+iz)$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::applyComplexInverseMultiVector
"virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::applyComplexInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const =0

Solve $(J+i\\\\omega B)(y+iz) = a+ib$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::computeDCeDp "virtual
NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::computeDCeDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &yVector,
const NOX::Abstract::Vector &zVector, double w,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag, bool isValid)=0

Computes the derivative $\\\\frac{\\\\partial (J+i\\\\omega
B)(y+iz)}{\\\\partial p}$ where $p$ is the parameter indexed by
paramIDs. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::computeDCeDxa "virtual
NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::computeDCeDxa(const
NOX::Abstract::Vector &yVector, const NOX::Abstract::Vector &zVector,
double w, const NOX::Abstract::MultiVector &aVector,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag)=0

Computes the directional derivative $\\\\frac{\\\\partial
(J+i\\\\omega B)(y+iz)}{\\\\partial x} a$ for the given direction $a$.
";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::computeDCeDxa "virtual
NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::AbstractGroup::computeDCeDxa(const
NOX::Abstract::Vector &yVector, const NOX::Abstract::Vector &zVector,
double w, const NOX::Abstract::MultiVector &aVector, const
NOX::Abstract::Vector &Ce_real, const NOX::Abstract::Vector &Ce_imag,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag)=0

Computes the directional derivative $\\\\frac{\\\\partial
(J+i\\\\omega B)(y+iz)}{\\\\partial x} a$ for the given direction $a$.
The arguments Ce_real and Ce_imag hold the real and imaginary
components of $(J+i\\\\omega B)(y+iz)$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::AbstractGroup "LOCA::Hopf::MooreSpence::AbstractGroup::AbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::AbstractGroup::~AbstractGroup "virtual
LOCA::Hopf::MooreSpence::AbstractGroup::~AbstractGroup()

Destructor. ";


// File: classLOCA_1_1BorderedSolver_1_1AbstractOperator.xml
%feature("docstring") LOCA::BorderedSolver::AbstractOperator "

Abstract interface class representing an operator for solving bordered
sets of linear equations.

C++ includes: LOCA_BorderedSolver_AbstractOperator.H ";

%feature("docstring")
LOCA::BorderedSolver::AbstractOperator::AbstractOperator "LOCA::BorderedSolver::AbstractOperator::AbstractOperator()

Constructor. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractOperator::~AbstractOperator "virtual
LOCA::BorderedSolver::AbstractOperator::~AbstractOperator()

Destructor. ";

%feature("docstring")  LOCA::BorderedSolver::AbstractOperator::apply "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractOperator::apply(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const =0

Apply the operator. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractOperator::applyTranspose "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractOperator::applyTranspose(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const =0

Apply transpose of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractOperator::applyInverse "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractOperator::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const =0

Apply inverse of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractOperator::applyInverseTranspose "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractOperator::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const =0

Apply inverse transpose of the operator. ";


// File: classLOCA_1_1BorderedSolver_1_1AbstractStrategy.xml
%feature("docstring") LOCA::BorderedSolver::AbstractStrategy "

Abstract interface class for solving bordered sets of linear
equations.

Abstract interface for solving systems of equations of the form \\\\[
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
F \\\\\\\\ G \\\\end{bmatrix} \\\\] where $J$ is an $n\\\\times n$
matrix, $A$ and $B$ are $n\\\\times m$, $C$ is $m\\\\times m$, $X$ and
$F$ are $n\\\\times p$ and $Y$ and $G$ are $m\\\\times p$. The action
of $J$ and its inverse are represnted by a
LOCA::BorderedSolver::AbstractOperator while $A$ is a
NOX::Abstract::MultiVector and $B$, $C$ are represtend by the solution
and parameter components of the derivative of a constraint contained
in LOCA::MultiContinuation::ConstraintInterface. All classes that
implement a method for computing solutions to this system of equations
should be derived from this class. Constructors for derived classes
should be of the form:

where global_data is the LOCA global data object, topParams is the
parsed top-level parameter list, and solverParams is a parameter list
of bordered-solver parameters.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_BorderedSolver_AbstractStrategy.H ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::AbstractStrategy "LOCA::BorderedSolver::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::~AbstractStrategy "virtual
LOCA::BorderedSolver::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::setMatrixBlocks "virtual void
LOCA::BorderedSolver::AbstractStrategy::setMatrixBlocks(const
Teuchos::RCP< const LOCA::BorderedSolver::AbstractOperator > &op,
const Teuchos::RCP< const NOX::Abstract::MultiVector > &blockA, const
Teuchos::RCP< const LOCA::MultiContinuation::ConstraintInterface >
&blockB, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > &blockC)=0

Set blocks.

The blockA or blockC pointer may be null if either is zero. Whether
block B is zero will be determined by querying blockB via
ConstraintInterface::isConstraintDerivativesXZero. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::setMatrixBlocksMultiVecConstraint
"void
LOCA::BorderedSolver::AbstractStrategy::setMatrixBlocksMultiVecConstraint(const
Teuchos::RCP< const LOCA::BorderedSolver::AbstractOperator > &op,
const Teuchos::RCP< const NOX::Abstract::MultiVector > &blockA, const
Teuchos::RCP< const NOX::Abstract::MultiVector > &blockB, const
Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix > &blockC)

Set blocks with multivector constraint.

This is a version of setMatrixBlocks that takes a multivector for
blockB. This method has a default implementation to generate a
LOCA::MultiContinuation::MultiVecConstraint from blockB which is then
passed to the setMatrixBlocks() method. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::initForSolve "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractStrategy::initForSolve()=0

Intialize solver for a solve.

This should be called after setMatrixBlocks(), but before
applyInverse(). ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::initForTransposeSolve "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractStrategy::initForTransposeSolve()=0

Intialize solver for a transpose solve.

This should be called after setMatrixBlocks(), but before
applyInverseTranspose(). ";

%feature("docstring")  LOCA::BorderedSolver::AbstractStrategy::apply "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractStrategy::apply(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const =0

Computed extended matrix-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} \\\\] where $U$ is
$n\\\\times p$, $V$ is $m\\\\times p$ and the other blocks are as
defined above. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::applyTranspose "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractStrategy::applyTranspose(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const =0

Computed extended matrix transpose-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J^T & B \\\\\\\\ A^T & C^T \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} \\\\] where $U$ is
$n\\\\times p$, $V$ is $m\\\\times p$ and the other blocks are as
defined above. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::applyInverse "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractStrategy::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const =0

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";

%feature("docstring")
LOCA::BorderedSolver::AbstractStrategy::applyInverseTranspose "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::AbstractStrategy::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const =0

Solves the transpose of the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1StepSize_1_1AbstractStrategy.xml
%feature("docstring") LOCA::StepSize::AbstractStrategy "

Abstract interface class for step size control strategies.

AbstractStrategy defines an abstract interface for step size control
strategies. It is used by LOCA::Stepper to the step size for each
continuation step.

The interface currently defines three pure virtual methods,
computeStepSize() to compute the step size, getPrevSteSize() to get
the step size from the previous step, and getStartStepSize() to get
the initial step size. Derived classes should implement this method
for a particular strategy. Constructors for derived classes should be
of the form:

where global_data is the LOCA global data object, topParams is the
parsed top-level parameter list, and stepsizeParams is a parameter
list of step size control parameters.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_StepSize_AbstractStrategy.H ";

%feature("docstring")
LOCA::StepSize::AbstractStrategy::AbstractStrategy "LOCA::StepSize::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::StepSize::AbstractStrategy::~AbstractStrategy "virtual
LOCA::StepSize::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")
LOCA::StepSize::AbstractStrategy::computeStepSize "virtual
NOX::Abstract::Group::ReturnType
LOCA::StepSize::AbstractStrategy::computeStepSize(LOCA::MultiContinuation::AbstractStrategy
&curGroup, const LOCA::MultiContinuation::ExtendedVector &predictor,
const NOX::Solver::Generic &solver, const
LOCA::Abstract::Iterator::StepStatus &stepStatus, const
LOCA::Abstract::Iterator &stepper, double &stepSize)=0

Compute step size.

Parameters:
-----------

curGroup:  [in] Current continuation group

predictor:  [in] Current predictor direction

solver:  [in] Solver from previous step

stepStatus:  [in] Status of previous step

stepper:  [in] Stepper

stepSize:  [out] Computed step size

ReturnType code indicating success or failure ";

%feature("docstring")
LOCA::StepSize::AbstractStrategy::getPrevStepSize "virtual double
LOCA::StepSize::AbstractStrategy::getPrevStepSize() const =0

Return the previous step size. ";

%feature("docstring")
LOCA::StepSize::AbstractStrategy::getStartStepSize "virtual double
LOCA::StepSize::AbstractStrategy::getStartStepSize() const =0

Return the initial step size. ";


// File: classLOCA_1_1Eigensolver_1_1AbstractStrategy.xml
%feature("docstring") LOCA::Eigensolver::AbstractStrategy "

Abstract interface class for Eigensolver strategies.

AbstractStrategy defines an abstract interface for eigensolver
strategies. It is used by LOCA::Stepper to compute eigenvalues of the
steady-state solution after each continuation step.

The interface currently defines one pure virtual method,
computeEigenvalues(), to compute the eigenvalues. Derived classes
should implement this method for a particular eigensolver strategy.
Constructors for derived classes should be of the form:

where global_data is the LOCA global data object, topParams is the
parsed top-level parameter list, and eigenParams is a parameter list
of eigensolver parameters.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_Eigensolver_AbstractStrategy.H ";

%feature("docstring")
LOCA::Eigensolver::AbstractStrategy::AbstractStrategy "LOCA::Eigensolver::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::Eigensolver::AbstractStrategy::~AbstractStrategy "virtual
LOCA::Eigensolver::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")
LOCA::Eigensolver::AbstractStrategy::computeEigenvalues "virtual
NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::AbstractStrategy::computeEigenvalues(NOX::Abstract::Group
&group, Teuchos::RCP< std::vector< double > > &evals_r, Teuchos::RCP<
std::vector< double > > &evals_i, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_r, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_i)=0

Compute eigenvalues/eigenvectors in group group.

Parameters:
-----------

group:  [in] NOX Group to compute eigenvalues of

evals_r:  [out] Real eigenvalues

evals_i:  [out] Imaginary eigenvalues

evecs_r:  [out] Real eigenvectors

evecs_i:  [out] Imaginary eigenvectors

ReturnType code indicating success or failure ";


// File: classLOCA_1_1AnasaziOperator_1_1AbstractStrategy.xml
%feature("docstring") LOCA::AnasaziOperator::AbstractStrategy "

Abstract interface class for Anasazi operator strategies.

AbstractStrategy defines an abstract interface for anasazi operators.
It is used by LOCA::Eigensolver::AnasaziStrategy to compute different
kinds of eigenvalues of the steady-state solution after each
continuation step.

The interface currently defines several pure virtual methods, apply(),
to apply the operator, transformEigenvalues() to transform the
computed eigenvalues back to eigenvalues of untransformed state,
rayleighQuotient to compute the rayleighQuotient for the operator, and
label() to return the name of the operator. Derived classes should
implement these method for a particular operator. Constructors for
derived classes should be of the form:

where global_data is the LOCA global data object, topParams is the
parsed top-level parameter list, eigenParams is a parameter list of
eigensolver parameters, solverParams is a parameter list of linear
solver parameters, and grp is the group representing the Jacobian and
mass matrices.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_AnasaziOperator_AbstractStrategy.H ";

%feature("docstring")
LOCA::AnasaziOperator::AbstractStrategy::AbstractStrategy "LOCA::AnasaziOperator::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::AnasaziOperator::AbstractStrategy::~AbstractStrategy "virtual
LOCA::AnasaziOperator::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")  LOCA::AnasaziOperator::AbstractStrategy::label
"virtual const std::string&
LOCA::AnasaziOperator::AbstractStrategy::label() const =0

Return name of this operator. ";

%feature("docstring")  LOCA::AnasaziOperator::AbstractStrategy::apply
"virtual void LOCA::AnasaziOperator::AbstractStrategy::apply(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &output)
const =0

Apply the operator to input with the result in output. ";

%feature("docstring")
LOCA::AnasaziOperator::AbstractStrategy::preProcessSeedVector "virtual void
LOCA::AnasaziOperator::AbstractStrategy::preProcessSeedVector(NOX::Abstract::MultiVector
&ivec)

Give strategy an opportunit to massage the random seed vector. ";

%feature("docstring")
LOCA::AnasaziOperator::AbstractStrategy::beginPostProcessing "virtual
void LOCA::AnasaziOperator::AbstractStrategy::beginPostProcessing()

Hook to precompute info for subsequent repeated calls to
tranformEigenvalue and rayleighQuotient. ";

%feature("docstring")
LOCA::AnasaziOperator::AbstractStrategy::transformEigenvalue "virtual
void
LOCA::AnasaziOperator::AbstractStrategy::transformEigenvalue(double
&ev_r, double &ev_i) const =0

Transform eigenvalue in place. ";

%feature("docstring")
LOCA::AnasaziOperator::AbstractStrategy::rayleighQuotient "virtual
NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::AbstractStrategy::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const =0

Compute Rayleigh quotient. ";


// File: classLOCA_1_1MultiContinuation_1_1AbstractStrategy.xml
%feature("docstring") LOCA::MultiContinuation::AbstractStrategy "

Abstract interface class for continuation strategies.

AbstractStrategy defines an abstract interface for continuation
strategies. This interface is used by the LOCA::Stepper to manipulate
continuation groups in a consistent manner. It defines a number pure
virtual methods that all continuation groups must implement.

C++ includes: LOCA_MultiContinuation_AbstractStrategy.H ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::AbstractStrategy "LOCA::MultiContinuation::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::~AbstractStrategy "virtual
LOCA::MultiContinuation::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")  LOCA::MultiContinuation::AbstractStrategy::copy
"virtual void LOCA::MultiContinuation::AbstractStrategy::copy(const
NOX::Abstract::Group &source)=0

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getNumParams "virtual int
LOCA::MultiContinuation::AbstractStrategy::getNumParams() const =0

Returns number of parameters. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::preProcessContinuationStep
"virtual void
LOCA::MultiContinuation::AbstractStrategy::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)=0

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::postProcessContinuationStep
"virtual void
LOCA::MultiContinuation::AbstractStrategy::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)=0

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::computePredictor "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::AbstractStrategy::computePredictor()=0

Compute predictor directions. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::isPredictor "virtual bool
LOCA::MultiContinuation::AbstractStrategy::isPredictor() const =0

Is Predictor valid. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::scaleTangent "virtual void
LOCA::MultiContinuation::AbstractStrategy::scaleTangent()=0

Scales tangent to predictor. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::setPredictorTangentDirection
"virtual void
LOCA::MultiContinuation::AbstractStrategy::setPredictorTangentDirection(const
LOCA::MultiContinuation::ExtendedVector &v, int i)=0

Sets tangent to predictor.

This is required by MF which takes the tangent space, orthogonalizes
it, and then sets it back in the group. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getPredictorTangent "virtual const LOCA::MultiContinuation::ExtendedMultiVector&
LOCA::MultiContinuation::AbstractStrategy::getPredictorTangent() const
=0

Returns tangent to predictor. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getScaledPredictorTangent "virtual const LOCA::MultiContinuation::ExtendedMultiVector&
LOCA::MultiContinuation::AbstractStrategy::getScaledPredictorTangent()
const =0

Returns scaled tangent to predictor. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::setPrevX "virtual void
LOCA::MultiContinuation::AbstractStrategy::setPrevX(const
NOX::Abstract::Vector &y)=0

Set the previous solution vector y. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getPrevX "virtual const
LOCA::MultiContinuation::ExtendedVector&
LOCA::MultiContinuation::AbstractStrategy::getPrevX() const =0

Gets the previous solution vector. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::setStepSize "virtual void
LOCA::MultiContinuation::AbstractStrategy::setStepSize(double deltaS,
int i=0)=0

Set step size for continuation constraint equation i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getStepSize "virtual
double LOCA::MultiContinuation::AbstractStrategy::getStepSize(int i=0)
const =0

Get step size for continuation constraint equation i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::setContinuationParameter "virtual void
LOCA::MultiContinuation::AbstractStrategy::setContinuationParameter(double
val, int i=0)=0

Sets the value for continuation parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameter "virtual double
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameter(int
i=0) const =0

Returns the value for continuation parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameterID
"virtual int
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameterID(int
i=0) const =0

Get the continuation parameter id for parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameterIDs
"virtual const std::vector<int>&
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameterIDs()
const =0

Get the continuation parameter ids. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameterName
"virtual std::string
LOCA::MultiContinuation::AbstractStrategy::getContinuationParameterName(int
i=0) const =0

Get the continuation parameter id for parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::getStepSizeScaleFactor "virtual double
LOCA::MultiContinuation::AbstractStrategy::getStepSizeScaleFactor(int
i=0) const =0

Returns step size scale factor for constraint equation i. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::printSolution "virtual
void LOCA::MultiContinuation::AbstractStrategy::printSolution() const
=0

Prints the group. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::computeScaledDotProduct "virtual double
LOCA::MultiContinuation::AbstractStrategy::computeScaledDotProduct(const
NOX::Abstract::Vector &x, const NOX::Abstract::Vector &y) const =0

Computes a scaled dot product between two continuation vectors. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::projectToDrawDimension "virtual int
LOCA::MultiContinuation::AbstractStrategy::projectToDrawDimension()
const =0

Returns dimension of project to draw array. ";

%feature("docstring")
LOCA::MultiContinuation::AbstractStrategy::projectToDraw "virtual
void LOCA::MultiContinuation::AbstractStrategy::projectToDraw(const
LOCA::MultiContinuation::ExtendedVector &x, double *px) const =0

Fills the project to draw array. ";


// File: classLOCA_1_1EigenvalueSort_1_1AbstractStrategy.xml
%feature("docstring") LOCA::EigenvalueSort::AbstractStrategy "

Abstract interface for eigenvalue sorting strategies.

AbstractStrategy defines an abstract interface for sorting
eigenvalues. It is used by LOCA::Eigensolver strategies to ensure the
desired eigenvalues are printed/saved.

The interface defines two pure virtual methods, sort(), to sort the
eigenvalues, optionally returning a permutation vector determining how
th eigenvalues were sorted. There is one version for real-only
eigenvalues and one version for complex. Derived classes should
implement these methods for a particular sorting strategy.
Constructors for derived classes should be of the form:

where global_data is the LOCA global data object and eigenParams is
the eigensolver parameter list (see
LOCA::Eigensolver::AbstractStrategy). In addition to any parameters
for the chosen sorting method, this list should contain the parameter
\"Sorting Order\" giving the name of the sorting strategy.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::AbstractStrategy::AbstractStrategy "LOCA::EigenvalueSort::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::EigenvalueSort::AbstractStrategy::~AbstractStrategy "virtual
LOCA::EigenvalueSort::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::AbstractStrategy::sort "virtual NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::AbstractStrategy::sort(int n, double *evals,
std::vector< int > *perm=NULL) const =0

Sort real eigenvalues optionally returning a permutation vector.

Parameters:
-----------

n:  [in] Number of eigenvalues

evals:  [in/out] Array of length n containing the eigenvalues to be
sorted.

perm:  [out] Vector of length n to store the permutation (optional).

Returns the status of the sorting routine ";

%feature("docstring")  LOCA::EigenvalueSort::AbstractStrategy::sort "virtual NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::AbstractStrategy::sort(int n, double *r_evals,
double *i_evals=NULL, std::vector< int > *perm=NULL) const =0

Sort complex eigenvalues optionally returning a permutation vector.

Parameters:
-----------

n:  [in] Number of eigenvalues

r_evals:  [in/out] Array of length n containing the real part of the
eigenvalues to be sorted.

i_evals:  [in/out] Array of length n containing the imaginary part of
the eigenvalues to be sorted.

perm:  [out] Vector of length n to store the permutation (optional).

Returns the status of the sorting routine ";


// File: classLOCA_1_1Epetra_1_1TransposeLinearSystem_1_1AbstractStrategy.xml
%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy "

A pure virtual interface for solving the transpose of a linear system.

C++ includes: LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::AbstractStrategy
"LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::~AbstractStrategy
"virtual
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::applyJacobianTransposeInverse
"virtual bool
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector
&result)=0

Applies the inverse of the Jacobian matrix transpose to the given
input vector and puts the answer in result.

Computes \\\\[ v = J^{-T} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The parameter list contains the linear solver options. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::createJacobianTranspose
"virtual bool
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::createJacobianTranspose()=0

Evaluates the Jacobian-transpose based on the solution vector x.

Note: For flexibility, this method does not compute the original
Jacobian matrix. It uses whatever is currently stored in the linear
system. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::createTransposePreconditioner
"virtual bool
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::createTransposePreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &p)=0

Explicitly constructs a preconditioner based on the solution vector x
and the parameter list p.

Note: x is only needed for user-supplied preconditioners. When using a
built- in preconditioner (e.g., Ifpack), x will note be used. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::getJacobianTransposeOperator
"virtual Teuchos::RCP<Epetra_Operator>
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::getJacobianTransposeOperator()=0

Get Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::getTransposePreconditioner
"virtual Teuchos::RCP<Epetra_Operator>
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::getTransposePreconditioner()=0

Get transpose-preconditioner. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::setJacobianTransposeOperator
"virtual void
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::setJacobianTransposeOperator(const
Teuchos::RCP< Epetra_Operator > &new_jac_trans)=0

Set Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::setTransposePreconditioner
"virtual void
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy::setTransposePreconditioner(const
Teuchos::RCP< Epetra_Operator > &new_prec_trans)=0

Set transpose-preconditioner. ";


// File: classLOCA_1_1MultiPredictor_1_1AbstractStrategy.xml
%feature("docstring") LOCA::MultiPredictor::AbstractStrategy "

Abstract interface class for predictor strategies.

AbstractStrategy defines an abstract interface for predictor
strategies. It is used by the LOCA::Stepper and
LOCA::MultiContinuation groups to compute predictors and tangent
approximations to continuation curves and surfaces.

The interface defines several pure virtual methods that derived
classes should implement for a particular predictor strategy. Note
that predictor strategies are assumed to have a state, and therefore
need to define copy constructors, assignment operators, and a clone
function. Constructors of derived classes should be of the form:

where global_data is the LOCA global data object, topParams is the
parsed top-level parameter list, and predictorParams is a parameter
list of predictor parameters.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_MultiPredictor_AbstractStrategy.H ";

%feature("docstring")
LOCA::MultiPredictor::AbstractStrategy::AbstractStrategy "LOCA::MultiPredictor::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::MultiPredictor::AbstractStrategy::~AbstractStrategy "virtual
LOCA::MultiPredictor::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::AbstractStrategy::clone "virtual Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::AbstractStrategy::clone(NOX::CopyType
type=NOX::DeepCopy) const =0

Clone function. ";

%feature("docstring")  LOCA::MultiPredictor::AbstractStrategy::compute
"virtual NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::AbstractStrategy::compute(bool baseOnSecant,
const std::vector< double > &stepSize,
LOCA::MultiContinuation::ExtendedGroup &grp, const
LOCA::MultiContinuation::ExtendedVector &prevXVec, const
LOCA::MultiContinuation::ExtendedVector &xVec)=0

Compute the predictor given the current and previous solution vectors.
Set baseOnSecant to false if the predictor orientation should not be
based on the secant vector (first or last steps of a continuation
run).

As an example for a first-order predictor, this method should compute
the approximate tangent to the continuation curve. ";

%feature("docstring")
LOCA::MultiPredictor::AbstractStrategy::evaluate "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::AbstractStrategy::evaluate(const std::vector<
double > &stepSize, const LOCA::MultiContinuation::ExtendedVector
&xVec, LOCA::MultiContinuation::ExtendedMultiVector &result) const =0

Evaluate predictor with step size stepSize.

For a first-order predictor, this method should compute result[i] =
xVec[i] + stepSize[i] * v[i] for each i, where v[i] is the ith
predictor direction. ";

%feature("docstring")
LOCA::MultiPredictor::AbstractStrategy::computeTangent "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::AbstractStrategy::computeTangent(LOCA::MultiContinuation::ExtendedMultiVector
&tangent)=0

Compute tangent to predictor and store in tangent.

For a first-order predictor, this is the predictor direction itself.
";

%feature("docstring")
LOCA::MultiPredictor::AbstractStrategy::isTangentScalable "virtual
bool LOCA::MultiPredictor::AbstractStrategy::isTangentScalable() const
=0

Is the tangent vector for this predictor scalable.

This method determines whether the approximate tangent computed by
this strategy is appropriate for scaling. ";


// File: classLOCA_1_1SaveEigenData_1_1AbstractStrategy.xml
%feature("docstring") LOCA::SaveEigenData::AbstractStrategy "

Abstract interface class strategies to save eigenvector/value data.

AbstractStrategy defines an abstract interface for saving eigenvectors
and eigenvalues that are computed at each continuation step. This is
important because this data is often useful for restarting
continuations near bifurcation points and gives the user flexibility
in how this data is stored.

The interface currently defines one pure virtual method, save(), to
save any eigenvectors or values as specified by the user. Derived
classes should implement this method for a particular strategy to save
this data, which is usually highly application code dependent.
Constructors for derived classes should be of the form:

where global_data is the LOCA global data object, topParams is the
parsed top-level parameter list, and eigenParams is a parameter list
of eigensolver parameters. This list should also specify which and how
many eigenvectors/values to save as defined by the strategy.

This class and its children follow the Strategy pattern as defined in
Erich Gamma, et al. \"Design Patterns:  Elements of Reusable   Object-
Oriented Software.\" Addison Wesley, Boston, MA, 1995.

C++ includes: LOCA_SaveEigenData_AbstractStrategy.H ";

%feature("docstring")
LOCA::SaveEigenData::AbstractStrategy::AbstractStrategy "LOCA::SaveEigenData::AbstractStrategy::AbstractStrategy()

Constructor. ";

%feature("docstring")
LOCA::SaveEigenData::AbstractStrategy::~AbstractStrategy "virtual
LOCA::SaveEigenData::AbstractStrategy::~AbstractStrategy()

Destructor. ";

%feature("docstring")  LOCA::SaveEigenData::AbstractStrategy::save "virtual NOX::Abstract::Group::ReturnType
LOCA::SaveEigenData::AbstractStrategy::save(Teuchos::RCP< std::vector<
double > > &evals_r, Teuchos::RCP< std::vector< double > > &evals_i,
Teuchos::RCP< NOX::Abstract::MultiVector > &evecs_r, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_i)=0

Save eigenvalues/eigenvectors.

Parameters:
-----------

evals_r:  [out] Real eigenvalues

evals_i:  [out] Imaginary eigenvalues

evecs_r:  [out] Real eigenvectors

evecs_i:  [out] Imaginary eigenvectors

ReturnType code indicating success or failure ";


// File: classLOCA_1_1StepSize_1_1Adaptive.xml
%feature("docstring") LOCA::StepSize::Adaptive "

Adaptive step size control strategy

This class implements an adaptive step size control strategy derived
from the strategy implemented in the LOCA::StepSize::Constant class.
If the previous step was unsucessful, the step size is cut in half as
in the constant strategy, but if the step was sucessful this strategy
increases the step size based on the number of nonlinear solver
iterations required in the previous step. In particular, the new step
size $\\\\Delta s_{new}$ is given by \\\\[ \\\\Delta s_{new} =
\\\\Delta s_{old}\\\\left(1 + a\\\\left(\\\\frac{N_{max} -
N}{N_{max}}\\\\right)^2\\\\right) \\\\] where $a\\\\in[0,1]$ is an
aggressiveness factor, $N$ is the number of nonlinear solver
iterations in the previous step, and $N_{max}$ is the maximum number
of nonlinear solver iterations.

The parameters used by this class supplied in the constructor are the
same as used by the Constant class in addition to: \"Aggressiveness\"
- Aggressiveness factor $a$ (Default 0.5)

C++ includes: LOCA_StepSize_Adaptive.H ";

%feature("docstring")  LOCA::StepSize::Adaptive::Adaptive "LOCA::StepSize::Adaptive::Adaptive(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepsizeParams)

Constructor. ";

%feature("docstring")  LOCA::StepSize::Adaptive::~Adaptive "LOCA::StepSize::Adaptive::~Adaptive()

Destructor. ";

%feature("docstring")  LOCA::StepSize::Adaptive::computeStepSize "NOX::Abstract::Group::ReturnType
LOCA::StepSize::Adaptive::computeStepSize(LOCA::MultiContinuation::AbstractStrategy
&curGroup, const LOCA::MultiContinuation::ExtendedVector &predictor,
const NOX::Solver::Generic &solver, const
LOCA::Abstract::Iterator::StepStatus &stepStatus, const
LOCA::Abstract::Iterator &stepper, double &stepSize)

Compute the step size as described above.

Parameters:
-----------

curGroup:  [in] Current continuation group

predictor:  [in] Current predictor direction

solver:  [in] Solver from previous step

stepStatus:  [in] Status of previous step

stepper:  [in] Stepper

stepSize:  [out] Computed step size

ReturnType code indicating success or failure Returns
NOX::Abstract::Group::Failed if the computed step size is smaller than
the minimum step size ";


// File: classLOCA_1_1Epetra_1_1AdaptiveSolutionManager.xml
%feature("docstring") LOCA::Epetra::AdaptiveSolutionManager "";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager "LOCA::Epetra::AdaptiveSolutionManager::AdaptiveSolutionManager(const
Teuchos::RCP< const Epetra_Map > &map_, const Teuchos::RCP< const
Epetra_Map > &overlapMap_, const Teuchos::RCP< const Epetra_CrsGraph >
&overlapJacGraph_) ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::~AdaptiveSolutionManager "virtual
LOCA::Epetra::AdaptiveSolutionManager::~AdaptiveSolutionManager() ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::updateSolution "Teuchos::RCP<
const Epetra_Vector >
LOCA::Epetra::AdaptiveSolutionManager::updateSolution() ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::buildSolutionGroup "virtual
Teuchos::RCP<LOCA::Epetra::Group>
LOCA::Epetra::AdaptiveSolutionManager::buildSolutionGroup()=0

Build the LOCA solution group. ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::projectCurrentSolution "void
LOCA::Epetra::AdaptiveSolutionManager::projectCurrentSolution()

Remap \"old\" solution into new data structures. ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::getConvergenceData "void
LOCA::Epetra::AdaptiveSolutionManager::getConvergenceData(int
&KrylovIters, int &lastSolveKrylocIters, int &linSolves, double
&tolAchieved) const

Accessor functions. ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::getAdaptParamsNonConst "virtual Teuchos::RCP<Teuchos::ParameterList>
LOCA::Epetra::AdaptiveSolutionManager::getAdaptParamsNonConst() ";

%feature("docstring")
LOCA::Epetra::AdaptiveSolutionManager::getAdaptParams "virtual
Teuchos::RCP<const Teuchos::ParameterList>
LOCA::Epetra::AdaptiveSolutionManager::getAdaptParams() const ";


// File: classLOCA_1_1Epetra_1_1AdaptiveStepper.xml
%feature("docstring") LOCA::Epetra::AdaptiveStepper "

Implementation of LOCA::Abstract::Iterator for computing points along
a continuation curve.

The AdaptiveStepper class implements the pure virtual methods of the
LOCA::Abstract::Iterator for iteratively computing points along a
continuation curve.

C++ includes: LOCA_Epetra_AdaptiveStepper.H ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::AdaptiveStepper
"LOCA::Epetra::AdaptiveStepper::AdaptiveStepper(const Teuchos::RCP<
Teuchos::ParameterList > &pList, const Teuchos::RCP<
LOCA::Epetra::AdaptiveSolutionManager > &mgr, const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
NOX::StatusTest::Generic > &nt) ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::~AdaptiveStepper
"LOCA::Epetra::AdaptiveStepper::~AdaptiveStepper()

Destructor. ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::eigensolverReset
"bool LOCA::Epetra::AdaptiveStepper::eigensolverReset(Teuchos::RCP<
Teuchos::ParameterList > &newEigensolverList)

Replaces the eigensolver parameter list. ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::getSolutionGroup
"Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Epetra::AdaptiveStepper::getSolutionGroup() const

Return the current solution group. ";

%feature("docstring")
LOCA::Epetra::AdaptiveStepper::getBifurcationGroup "Teuchos::RCP<
const LOCA::MultiContinuation::AbstractGroup >
LOCA::Epetra::AdaptiveStepper::getBifurcationGroup() const

Return the current bifurcation group.

If the current bifurcation method is \"None\", then the returned group
is the same as getSolutionGroup(), otherwise this method returns the
current bifurcation group (e.g., a turning point group). ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::getList "Teuchos::RCP< const Teuchos::ParameterList >
LOCA::Epetra::AdaptiveStepper::getList() const

Return the output parameters from the stepper algorithm. ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::getSolver "Teuchos::RCP< const NOX::Solver::Generic >
LOCA::Epetra::AdaptiveStepper::getSolver() const

Return the current nonlinear solver pointer.

Will throw an error if the solver does not exist yet. ";

%feature("docstring")
LOCA::Epetra::AdaptiveStepper::getContinuationParameter "double
LOCA::Epetra::AdaptiveStepper::getContinuationParameter() const

Return the current continuation parameter from the underlying
LOCA::MultiContinuation::AbstractStrategy. ";

%feature("docstring")  LOCA::Epetra::AdaptiveStepper::run "LOCA::Abstract::Iterator::IteratorStatus
LOCA::Epetra::AdaptiveStepper::run()

Run the iterator. ";


// File: classLOCA_1_1Eigensolver_1_1AnasaziStrategy.xml
%feature("docstring") LOCA::Eigensolver::AnasaziStrategy "

Anasazi eigensolver strategy.

This class implements an eigensolver strategy using the generic
Trilinos eigensolver package Anasazi. In particular, this strategy
uses the Anasazi::BlockKrylovSchur solver. Since Anasazi is a generic
solver, this strategy will work with any group implementation. This
strategy references the following parameters passed through the
eigenParams argument to the constructor (this list is passed directly
to the Anasazi::BlockKrylovSchulSolMgr solver manager): \"Operator\"
-- [string] (default: \"Jacobian inverse\") Operator to compute
eigenvalues of.

\"Block Size\" -- [int] (default: 1) Block size

\"Num Blocks\" -- [int] (default: 30) Maximum number of blocks (equals
the maximum length of the Arnoldi factorization

\"Num Eigenvalues\" -- [int] (default: 4) Number of requested
eigenvalues

\"Convergence Tolerance\" -- [double] (default: 1.0e-7) Tolerance for
the converged eigenvalues

\"Step Size\" -- [int] (default: 1) Checks convergence every so many
steps

\"Maximum Restarts\" -- [int] (default: 1) Number of restarts allowed

\"Symmetric\" -- [bool] (default: false) Is the operator symmetric

\"Verbosity\" -- [Anasazi::MsgType] (default:
Anasazi::Errors+Anasazi::Warnings+Anasazi::FinalSummary) Verbosity
level

\"Sorting Order\" -- [string\" (default: \"LM\") Sorting order of
printed eigenvalues

C++ includes: LOCA_Eigensolver_AnasaziStrategy.H ";

%feature("docstring")
LOCA::Eigensolver::AnasaziStrategy::AnasaziStrategy "LOCA::Eigensolver::AnasaziStrategy::AnasaziStrategy(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list. Used to obtain
linear-solver parameters and to instantiate sorting strategy.

eigenParams:  [in] Eigensolver parameters as described above. solver.
";

%feature("docstring")
LOCA::Eigensolver::AnasaziStrategy::~AnasaziStrategy "LOCA::Eigensolver::AnasaziStrategy::~AnasaziStrategy()

Destructor. ";

%feature("docstring")
LOCA::Eigensolver::AnasaziStrategy::computeEigenvalues "NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::AnasaziStrategy::computeEigenvalues(NOX::Abstract::Group
&group, Teuchos::RCP< std::vector< double > > &evals_r, Teuchos::RCP<
std::vector< double > > &evals_i, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_r, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_i)

Compute eigenvalues/eigenvectors.

The implementation here the sets up and calls the Anasazi
BlockKrylovSchur solver for computing eigenvalues. ";


// File: classLOCA_1_1MultiContinuation_1_1ArcLengthConstraint.xml
%feature("docstring") LOCA::MultiContinuation::ArcLengthConstraint "

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
arclength continuation.

This class implements the arclength constraint equation for pseudo-
arclength continuation: \\\\[
g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)= (x-x^\\\\ast)^Tv_x +
(p-p^\\\\ast) v_p - \\\\Delta s \\\\] where $v_x$, $v_p$ are the
solution and parameter components of the predictor direction $v$
respectively. Since the derivative of $g$ with respect to $x$ is just
$v$, the predictor tangent, this class implements the MVDX version of
the constraint interface.

C++ includes: LOCA_MultiContinuation_ArcLengthConstraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::copy "void
LOCA::MultiContinuation::ArcLengthConstraint::copy(const
ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::clone "Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface >
LOCA::MultiContinuation::ArcLengthConstraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::numConstraints "int
LOCA::MultiContinuation::ArcLengthConstraint::numConstraints() const

Return number of constraints. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::setX "void
LOCA::MultiContinuation::ArcLengthConstraint::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::setParam "void
LOCA::MultiContinuation::ArcLengthConstraint::setParam(int paramID,
double val)

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::setParams "void
LOCA::MultiContinuation::ArcLengthConstraint::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::computeConstraints "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ArcLengthConstraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::isConstraints "bool
LOCA::MultiContinuation::ArcLengthConstraint::isConstraints() const

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::isDX "bool
LOCA::MultiContinuation::ArcLengthConstraint::isDX() const

Return true if derivatives of constraints w.r.t. x are valid. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::getConstraints "const
NOX::Abstract::MultiVector::DenseMatrix &
LOCA::MultiContinuation::ArcLengthConstraint::getConstraints() const

Return constraint residuals. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::getDX "const
NOX::Abstract::MultiVector *
LOCA::MultiContinuation::ArcLengthConstraint::getDX() const

Return solution component of constraint derivatives. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::isDXZero "bool
LOCA::MultiContinuation::ArcLengthConstraint::isDXZero() const

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint "LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::MultiContinuation::ArcLengthGroup > &grp)

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint "LOCA::MultiContinuation::ArcLengthConstraint::ArcLengthConstraint(const
ArcLengthConstraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::~ArcLengthConstraint "LOCA::MultiContinuation::ArcLengthConstraint::~ArcLengthConstraint()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthConstraint::setArcLengthGroup "void
LOCA::MultiContinuation::ArcLengthConstraint::setArcLengthGroup(const
Teuchos::RCP< LOCA::MultiContinuation::ArcLengthGroup > &grp)

Set pointer to arclength group. ";


// File: classLOCA_1_1MultiContinuation_1_1ArcLengthGroup.xml
%feature("docstring") LOCA::MultiContinuation::ArcLengthGroup "

Specialization of LOCA::MultiContinuation::ExtendedGroup to pseudo-
arclength continuation.

Pseudo arc-length continuation corresponds to a continuation equation
$g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)=0$ with $g$ given by
\\\\[ g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)=
(x-x^\\\\ast)^Tv_x + (p-p^\\\\ast) v_p - \\\\Delta s \\\\] where
$v_x$, $v_p$ are the solution and parameter components of the
predictor direction $v$ respectively. This corresponds geometrically
to constraining the nonlinear solver steps used in calculating
$F(x,p)=0$ to be orthogonal to the predictor direction $v$. The
arclength constraint $g$ is represented by a
LOCA::MultiContinuation::ArcLengthConstraint object.

This class also reimplements the scaleTangent() and
computeScaledDotProduct() methods to implement a scaling method that
tries to ensure the solution and parameter contributions to the arc-
length equation are of the same order. Specifically, the arc-length
equation is replaced by \\\\[ (x-x^\\\\ast)^Tv_x +
\\\\theta^2(p-p^\\\\ast) v_p - \\\\Delta s = 0 \\\\] where $\\\\theta$
is chosen so that $\\\\theta^2 v_p$ is equal to a target value, 0.5 by
default. Parameters for this scaling method are passed through the
continuationParams argument to the constructor and are: \"Enable Arc
Length Scaling\" -- [bool] (default: true)

\"Initial Scale Factor\" -- [double] (default: 1.0)

\"Goal Arc Length Parameter Contribution\" -- [double] (default: 0.5)

\"Max Arc Length Parameter Contribution\" -- [double] (default: 0.8)

\"Min Scale Factor\" -- [double] (default: 1.0e-3)  Whether this
scaling method is used is determined by the \"Enable Arc Length
Scaling\", and the initial value for $\\\\theta$ is given by \"Initial
Scale Factor\". A new value of $\\\\theta$ is chosen only if
$\\\\theta^2 v_p$ is larger than the value given by \"Max Arc Length
Parameter Contribution\" and \"Min Scale Factor\" provides a minimum
value for $\\\\theta$.

C++ includes: LOCA_MultiContinuation_ArcLengthGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")  LOCA::MultiContinuation::ArcLengthGroup::clone
"Teuchos::RCP< NOX::Abstract::Group >
LOCA::MultiContinuation::ArcLengthGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. ";

/*  Implementation of LOCA::MultiContinuation::AbstractStrategy
virtual methods  */

%feature("docstring")  LOCA::MultiContinuation::ArcLengthGroup::copy "void LOCA::MultiContinuation::ArcLengthGroup::copy(const
NOX::Abstract::Group &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthGroup::scaleTangent "void
LOCA::MultiContinuation::ArcLengthGroup::scaleTangent()

Scales predictor. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthGroup::computeScaledDotProduct "double
LOCA::MultiContinuation::ArcLengthGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &x, const NOX::Abstract::Vector &y) const

Computes a scaled dot product between two continuation vectors. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthGroup::recalculateScaleFactor "void
LOCA::MultiContinuation::ArcLengthGroup::recalculateScaleFactor(double
dpds, double thetaOld, double &thetaNew)

Calculates scale factors. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup "LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &continuationParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, const Teuchos::RCP<
LOCA::MultiPredictor::AbstractStrategy > &pred, const std::vector< int
> &paramIDs)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list.

continuationParams:  [in] Continuation parameters as described above.

grp:  [in] Group representing $F$.

pred:  [in] Predictor strategy.

paramIDs:  [in] Parameter IDs of continuation parameters. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup "LOCA::MultiContinuation::ArcLengthGroup::ArcLengthGroup(const
ArcLengthGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ArcLengthGroup::~ArcLengthGroup "LOCA::MultiContinuation::ArcLengthGroup::~ArcLengthGroup()

Destructor. ";


// File: classLOCA_1_1Epetra_1_1AugmentedOp.xml
%feature("docstring") LOCA::Epetra::AugmentedOp "

Epetra operator representing a $n+m$ bordered matrix.

The LOCA::Epetra::AugmentedOp is an Epetra_Operator representing the
$n+m$ bordered matrix \\\\[ \\\\begin{bmatrix} J & A \\\\\\\\ B^T & C
\\\\end{bmatrix} \\\\] where $J$ is an Epetra_Operator representing an
$n\\\\times n$ matrix, and $A$ and $B$ are length $n$
Epetra_MultiVector's with $m$ columns, and $C$ is an $m\\\\times m$
dense matrix. It is assumed the Epetra_Map's for $A$, $B$, and $J$ are
the same and the corresponding map for the bordered matrix is
constructed from this map by storing the additional components on
processor 0. The buildEpetraAugmentedMultiVec() method can be used to
construct an Epetra_MultiVector using this map, a supplied length $n$
Epetra_MultiVector and an $m\\\\times m$ matrix, while
setEpetraAugmentedMultiVec() splits an extended multivector into its
length $n$ and $m$ components. The Apply() method performs the
$n+m\\\\times n+m$ matrix multiplication while ApplyInverse() uses a
block-elimination algorithm to compute the inverse using the
ApplyInverse() method of the underlying operator $J$. In this way,
linear systems of the form \\\\[ \\\\begin{bmatrix} J & A \\\\\\\\ B^T
& C \\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix}
= \\\\begin{bmatrix} F \\\\\\\\ G \\\\end{bmatrix} \\\\] can be solved
in a matrix-free mode using the Apply() method. This operator can also
represent a preconditioner of the form \\\\[ \\\\begin{bmatrix} M & A
\\\\\\\\ B^T & C \\\\end{bmatrix} \\\\] using the ApplyInvese()
method, where $M$ is a preconditioner for $J$. Note that if $J$ is
nearly singular, the preconditioner should not be too good because
otherwise the preconditining operation represented by ApplyInverse()
becomes unstable.

C++ includes: LOCA_Epetra_AugmentedOp.H ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::AugmentedOp "LOCA::Epetra::AugmentedOp::AugmentedOp(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP< Epetra_Operator >
&jac, const Teuchos::RCP< const Epetra_MultiVector > &a, const
Teuchos::RCP< const Epetra_MultiVector > &b, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > c)

Constructor.

Builds the bordered operator using the supplied operator jac and
Epetra_Vector's a and b. It is assumed a, b, and jac all have the same
map. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::~AugmentedOp "LOCA::Epetra::AugmentedOp::~AugmentedOp()

Destructor. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::SetUseTranspose "int LOCA::Epetra::AugmentedOp::SetUseTranspose(bool UseTranspose)

If set true, transpose of this operator will be applied.

Note that is only valid if the underlying operator $J$ supports a
transpose. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::Apply "int
LOCA::Epetra::AugmentedOp::Apply(const Epetra_MultiVector &Input,
Epetra_MultiVector &Result) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector Input in Result.

Computes the extended matrix-vector product \\\\[ \\\\begin{bmatrix} J
& A \\\\\\\\ B^T & C \\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y
\\\\end{bmatrix} = \\\\begin{bmatrix} JX + AY \\\\\\\\ B^T X + CY
\\\\end{bmatrix} \\\\] or its transpose if UseTranpose() is true. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::ApplyInverse "int
LOCA::Epetra::AugmentedOp::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector Input in Result.

Solves the extended system \\\\[ \\\\begin{bmatrix} J & A \\\\\\\\ B^T
& C \\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix}
= \\\\begin{bmatrix} F \\\\\\\\ G \\\\end{bmatrix} \\\\] using the
following block-elimination algorithm: \\\\[ \\\\begin{split} X_1 &=
J^{-1} F, \\\\\\\\ X_2 &= J^{-1} A, \\\\\\\\ Y &= (C-B^T
X_2)^{-1}(G-B^T X_1), \\\\\\\\ X &= X_1 - X_2 Y \\\\end{split} \\\\]
If UseTranpose() is true, the tranpose of the system is solved. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::NormInf "double
LOCA::Epetra::AugmentedOp::NormInf() const

Returns the infinity norm of the bordered matrix.

This is defined only if NormInf() of the underlying operator $J$ is
defined and is given by
$\\\\|J\\\\|_\\\\infty+\\\\|A\\\\|_\\\\infty+\\\\|B\\\\|_\\\\infty$.
";

%feature("docstring")  LOCA::Epetra::AugmentedOp::Label "const char *
LOCA::Epetra::AugmentedOp::Label() const

Returns a character std::string describing the operator. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::UseTranspose "bool
LOCA::Epetra::AugmentedOp::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::HasNormInf "bool
LOCA::Epetra::AugmentedOp::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::Comm "const
Epetra_Comm & LOCA::Epetra::AugmentedOp::Comm() const

Returns a reference to the Epetra_Comm communicator associated with
this operator. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::OperatorDomainMap "const Epetra_Map & LOCA::Epetra::AugmentedOp::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::OperatorRangeMap "const Epetra_Map & LOCA::Epetra::AugmentedOp::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";

%feature("docstring")  LOCA::Epetra::AugmentedOp::init "void
LOCA::Epetra::AugmentedOp::init(const Epetra_MultiVector &x)

Initialiazes operator for a solve. ";

%feature("docstring")
LOCA::Epetra::AugmentedOp::buildEpetraAugmentedMultiVec "Teuchos::RCP< Epetra_MultiVector >
LOCA::Epetra::AugmentedOp::buildEpetraAugmentedMultiVec(const
Epetra_MultiVector &x, const NOX::Abstract::MultiVector::DenseMatrix
*y, bool doCopy) const

Builds an extended vector from components.

Builds an extended vector using the map representing the bordered
matrix. If doCopy is true, the contents of x are copied into the
extended vector, otherwise only space for the extended vector is
created. ";

%feature("docstring")
LOCA::Epetra::AugmentedOp::setEpetraAugmentedMultiVec "void
LOCA::Epetra::AugmentedOp::setEpetraAugmentedMultiVec(Epetra_MultiVector
&x, NOX::Abstract::MultiVector::DenseMatrix &y, const
Epetra_MultiVector &augMultiVec) const

Sets components from extended vector.

Splits the extended vector augMultiVec into components x and y by
copying values out of extVec. ";


// File: classLOCA_1_1BorderedSolver_1_1BorderedOperator.xml
%feature("docstring") LOCA::BorderedSolver::BorderedOperator "

Bordered solver operator representing as bordered Jacobian as operator
as implemented in the NOX::Abstract::Group.

C++ includes: LOCA_BorderedSolver_BorderedOperator.H ";

%feature("docstring")
LOCA::BorderedSolver::BorderedOperator::BorderedOperator "LOCA::BorderedSolver::BorderedOperator::BorderedOperator(const
Teuchos::RCP< const LOCA::BorderedSystem::AbstractGroup > &grp)

Constructor. ";

%feature("docstring")
LOCA::BorderedSolver::BorderedOperator::~BorderedOperator "virtual
LOCA::BorderedSolver::BorderedOperator::~BorderedOperator()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::BorderedOperator::getBorderedGroup "virtual
Teuchos::RCP<const LOCA::BorderedSystem::AbstractGroup>
LOCA::BorderedSolver::BorderedOperator::getBorderedGroup() const ";

%feature("docstring")  LOCA::BorderedSolver::BorderedOperator::apply "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::BorderedOperator::apply(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const

Apply the operator. ";

%feature("docstring")
LOCA::BorderedSolver::BorderedOperator::applyTranspose "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::BorderedOperator::applyTranspose(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const

Apply transpose of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::BorderedOperator::applyInverse "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::BorderedOperator::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const

Apply inverse of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::BorderedOperator::applyInverseTranspose "virtual NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::BorderedOperator::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const

Apply inverse transpose of the operator. ";


// File: classLOCA_1_1BorderedSolver_1_1Bordering.xml
%feature("docstring") LOCA::BorderedSolver::Bordering "

Bordered system solver strategy based on bordering.

This class solves the extended system of equations \\\\[
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
F \\\\\\\\ G \\\\end{bmatrix} \\\\] via bordering (block elimination):
\\\\[ \\\\begin{aligned} X_1 &= J^{-1} F \\\\\\\\ X_2 &= J^{-1} A
\\\\\\\\ Y &= (C-B^T X_2)^{-1}(G-B^T X_1) \\\\\\\\ X &= X_1 - X_2 Y
\\\\end{aligned} \\\\] It takes advantage of any of the matrix blocks
being zero and concatenates $F$ and $A$ into a contiguous multivector
to compute $X_1$ and $X_2$ in one block solve.

To solve the transpose of the system, a similar bordering algorithm is
implemented. Note however that for the transpose, the constraint
object representing $B$ must implement the
LOCA::MultiContinuation::ConstraintInterfaceMVDX since $B$ appears on
the right-hand-side of a linear system.

C++ includes: LOCA_BorderedSolver_Bordering.H ";

%feature("docstring")  LOCA::BorderedSolver::Bordering::Bordering "LOCA::BorderedSolver::Bordering::Bordering(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Currently none are
referenced. ";

%feature("docstring")  LOCA::BorderedSolver::Bordering::~Bordering "LOCA::BorderedSolver::Bordering::~Bordering()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::Bordering::setMatrixBlocks "void
LOCA::BorderedSolver::Bordering::setMatrixBlocks(const Teuchos::RCP<
const LOCA::BorderedSolver::AbstractOperator > &op, const
Teuchos::RCP< const NOX::Abstract::MultiVector > &blockA, const
Teuchos::RCP< const LOCA::MultiContinuation::ConstraintInterface >
&blockB, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > &blockC)

Set blocks.

The blockA or blockC pointer may be null if either is zero. Whether
block B is zero will be determined by querying blockB via
ConstraintInterface::isConstraintDerivativesXZero. ";

%feature("docstring")  LOCA::BorderedSolver::Bordering::initForSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Bordering::initForSolve()

Intialize solver for a solve.

This should be called after setMatrixBlocks(), but before
applyInverse(). ";

%feature("docstring")
LOCA::BorderedSolver::Bordering::initForTransposeSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Bordering::initForTransposeSolve()

Intialize solver for a transpose solve.

This should be called after setMatrixBlocks(), but before
applyInverseTranspose(). ";

%feature("docstring")  LOCA::BorderedSolver::Bordering::apply "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Bordering::apply(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J*X + A*Y \\\\\\\\ B^T*X + C*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")  LOCA::BorderedSolver::Bordering::applyTranspose
"NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Bordering::applyTranspose(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix transpose-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J^T & B \\\\\\\\ A^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J^T*X + B*Y \\\\\\\\ A^T*X + C^T*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")  LOCA::BorderedSolver::Bordering::applyInverse "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Bordering::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system as defined above using bordering.

The params argument is the linear solver parameters. If isZeroF or
isZeroG is true, than the corresponding F or G pointers may be NULL.
";

%feature("docstring")
LOCA::BorderedSolver::Bordering::applyInverseTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Bordering::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the transpose of the extended system as defined above using
bordering.

The params argument is the linear solver parameters. If isZeroF or
isZeroG is true, than the corresponding F or G pointers may be NULL.
Note that for the transpose solve B must be of type
LOCA::MultiContinuation::ConstraintInterfaceMVDX. ";


// File: classLOCA_1_1AnasaziOperator_1_1Cayley.xml
%feature("docstring") LOCA::AnasaziOperator::Cayley "

Anasazi operator for computing generalized eigenvalues using Cayley
transformations.

This class implements the LOCA::AnasaziOperator::AbstractStrategy
interface for computing generalized eigenvalues $\\\\lambda$ and
eigenvectors $z$ of the system \\\\[ J z = \\\\lambda M z *\\\\] where
$J$ is the Jacobian matrix and $M$ is the mass matrix. The eigenvalues
are computed using a Cayley transformation, i.e. solving \\\\[ (J -
\\\\sigma M) z = (J - \\\\mu M) r \\\\] where $\\\\sigma$ is the
Cayley pole and $\\\\mu$ is the Cayley zero.

The parameters used by this class supplied in the constructor are:
\"Cayley Pole\" - $\\\\sigma$ as defined above (Default 0.0)

\"Cayley Zero\" - $\\\\mu$ as defined above (Default 0.0)

Also the grp argument to the constructor must be a child of
LOCA::TimeDependent::AbstractGroup for the shift-invert operations.

C++ includes: LOCA_AnasaziOperator_Cayley.H ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley::Cayley "LOCA::AnasaziOperator::Cayley::Cayley(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
LOCA::TimeDependent::AbstractGroup > &grp)

Constructor.

Argument grp must be of type LOCA::TimeDependent::AbstractGroup. See
class description for a list of eigenParams. ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley::~Cayley "LOCA::AnasaziOperator::Cayley::~Cayley()

Destructor. ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley::label "const
std::string & LOCA::AnasaziOperator::Cayley::label() const

Return name of this operator. ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley::apply "void
LOCA::AnasaziOperator::Cayley::apply(const NOX::Abstract::MultiVector
&input, NOX::Abstract::MultiVector &output) const

Apply the operator.

Applies the inverse of the shifted operator, i.e., solves \\\\[
(J-\\\\omega I)z = M r \\\\] for $z$, where $r = \\\\mbox{input}$ and
$z = \\\\mbox{output}$. ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley::preProcessSeedVector "void
LOCA::AnasaziOperator::Cayley::preProcessSeedVector(NOX::Abstract::MultiVector
&ivec)

PreProcess the random seed vector.

Performs one backward Euler iteration on the random initial seed
vector, to satisfy contraints ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley::transformEigenvalue "void
LOCA::AnasaziOperator::Cayley::transformEigenvalue(double &ev_r,
double &ev_i) const

Transform eigenvalue.

Transforms the given eigenvalue to the eigenvalue of the Jacobian-mass
matrix system by shifting and inverting it. ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley::rayleighQuotient
"NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::Cayley::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const

Compute Rayleigh quotient.

Computes the Rayleigh quotient $z^T J z / z^T M z$ for the eigenvector
$z$. ";


// File: classLOCA_1_1AnasaziOperator_1_1Cayley2Matrix.xml
%feature("docstring") LOCA::AnasaziOperator::Cayley2Matrix "

Anasazi operator for computing generalized eigenvalues using Cayley
transformations.

This class implements the LOCA::AnasaziOperator::AbstractStrategy
interface for computing generalized eigenvalues $\\\\lambda$ and
eigenvectors $z$ of the system \\\\[ J z = \\\\lambda M z *\\\\] where
$J$ is the Jacobian matrix and $M$ is the mass matrix. The eigenvalues
are computed using a Cayley transformation, i.e. solving \\\\[ (J -
\\\\sigma M) z = (J - \\\\mu M) r \\\\] where $\\\\sigma$ is the
Cayley pole and $\\\\mu$ is the Cayley zero.

The parameters used by this class supplied in the constructor are:
\"Cayley Pole\" - $\\\\sigma$ as defined above (Default 0.0)

\"Cayley Zero\" - $\\\\mu$ as defined above (Default 0.0)

Also the grp argument to the constructor must be a child of
LOCA::TimeDependent::AbstractGroup for the shift-invert operations.

C++ includes: LOCA_AnasaziOperator_Cayley2Matrix.H ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley2Matrix::Cayley2Matrix "LOCA::AnasaziOperator::Cayley2Matrix::Cayley2Matrix(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
LOCA::TimeDependent::AbstractGroup > &grp)

Constructor.

Argument grp must be of type LOCA::TimeDependent::AbstractGroup. See
class description for a list of eigenParams. ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley2Matrix::~Cayley2Matrix "LOCA::AnasaziOperator::Cayley2Matrix::~Cayley2Matrix()

Destructor. ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley2Matrix::label "const std::string & LOCA::AnasaziOperator::Cayley2Matrix::label()
const

Return name of this operator. ";

%feature("docstring")  LOCA::AnasaziOperator::Cayley2Matrix::apply "void LOCA::AnasaziOperator::Cayley2Matrix::apply(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &output)
const

Apply the operator.

Applies the inverse of the shifted operator, i.e., solves \\\\[
(J-\\\\omega I)z = M r \\\\] for $z$, where $r = \\\\mbox{input}$ and
$z = \\\\mbox{output}$. ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley2Matrix::preProcessSeedVector "void
LOCA::AnasaziOperator::Cayley2Matrix::preProcessSeedVector(NOX::Abstract::MultiVector
&ivec)

PreProcess the random seed vector.

Performs one backward Euler iteration on the random initial seed
vector, to satisfy contraints ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley2Matrix::beginPostProcessing "void
LOCA::AnasaziOperator::Cayley2Matrix::beginPostProcessing()

Begin PostProcessing of eigenvalues.

Compute Jacobian and mass matrix once, for use in subsequent repeated
calls to rayleighQuotient ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley2Matrix::transformEigenvalue "void
LOCA::AnasaziOperator::Cayley2Matrix::transformEigenvalue(double
&ev_r, double &ev_i) const

Transform eigenvalue.

Transforms the given eigenvalue to the eigenvalue of the Jacobian-mass
matrix system by shifting and inverting it. ";

%feature("docstring")
LOCA::AnasaziOperator::Cayley2Matrix::rayleighQuotient "NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::Cayley2Matrix::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const

Compute Rayleigh quotient.

Computes the Rayleigh quotient $z^T J z / z^T M z$ for the eigenvector
$z$. ";


// File: classLOCA_1_1StatusTest_1_1Combo.xml
%feature("docstring") LOCA::StatusTest::Combo "

Arbitrary combination of status tests.

In the AND (see LOCA::StatusTest::Combo::ComboType) combination, the
result is Unconverged (see LOCA::StatusTest::StatusType) if any of the
tests is Unconverged. Otherwise, the result is equal to the result of
the first test in the list that is either Converged or Failed. It is
not recommended to mix Converged and Failed tests in an AND
combination.

In the OR combination, the result is Unconverged if all of the tests
are Unconverged. Otherwise, it is the result of the first test in the
list that is either Converged or Failed. Therefore, it will generally
make sense to put the Failed -type tests at the end of the OR list.

We call checkStatus on every convergence test, though some may be
called with the LOCA::StatusTest::None option.

C++ includes: LOCA_StatusTest_Combo.H ";

%feature("docstring")  LOCA::StatusTest::Combo::Combo "LOCA::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP< const
LOCA::GlobalData > globalDataPtr=Teuchos::null)

Constructor. Optional argument is the error stream for output. ";

%feature("docstring")  LOCA::StatusTest::Combo::Combo "LOCA::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP<
Abstract > &a, const Teuchos::RCP< const LOCA::GlobalData >
globalDataPtr=Teuchos::null)

Constructor with a single test. ";

%feature("docstring")  LOCA::StatusTest::Combo::Combo "LOCA::StatusTest::Combo::Combo(ComboType t, const Teuchos::RCP<
Abstract > &a, const Teuchos::RCP< Abstract > &b, const Teuchos::RCP<
const LOCA::GlobalData > globalDataPtr=Teuchos::null)

Constructor with two tests. ";

%feature("docstring")  LOCA::StatusTest::Combo::addStatusTest "LOCA::StatusTest::Combo & LOCA::StatusTest::Combo::addStatusTest(const
Teuchos::RCP< Abstract > &a)

Add another test to this combination.

Calls isSafe() to determine if it is safe to add a to the combination.
";

%feature("docstring")  LOCA::StatusTest::Combo::~Combo "LOCA::StatusTest::Combo::~Combo()

Destructor. ";

%feature("docstring")  LOCA::StatusTest::Combo::checkStatus "virtual
LOCA::StatusTest::StatusType
LOCA::StatusTest::Combo::checkStatus(const LOCA::Abstract::Iterator
&stepper, LOCA::StatusTest::CheckType checkType)

Tests stopping criterion.

See addOp() and orOp() for details. ";

%feature("docstring")  LOCA::StatusTest::Combo::getStatus "LOCA::StatusTest::StatusType LOCA::StatusTest::Combo::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  LOCA::StatusTest::Combo::print "std::ostream &
LOCA::StatusTest::Combo::print(std::ostream &stream, int indent=0)
const

Output formatted description of stopping test to output stream. ";


// File: classLOCA_1_1Epetra_1_1CompactWYOp.xml
%feature("docstring") LOCA::Epetra::CompactWYOp "

An Epetra operator for solving extended sets of equations using
Householder transformations.

This class implements the $P$ operator as described in the
LOCA::BorderedSolver::EpetraHouseholder documentation for solving an
extended set of equations. It uses the $Q$ factor from a QR
factorization using the compact WY representation.

C++ includes: LOCA_Epetra_CompactWYOp.H ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::CompactWYOp "LOCA::Epetra::CompactWYOp::CompactWYOp(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP< const
Epetra_Operator > &jacOperator, const Teuchos::RCP< const
Epetra_MultiVector > &A_multiVec, const Teuchos::RCP< const
Epetra_MultiVector > &Y_x_multiVec, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > &Y_p_matrix, const
Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
&T_matrix)

Constructor.

Parameters:
-----------

global_data:  [in] The global data object

jacOperator:  [in] Jacobian operator J

A_multiVec:  [in] Multivector representing A

Y_x_multiVec:  [in] Multivector representing the solution component of
the Y matrix in the compact WY representation

Y_p_matrix:  [in] Matrix representing the parameter component of the Y
matrix in the compact WY representation

T_matrix:  [in] Matrix representing the T matrix in the compact WY
representation. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::~CompactWYOp "LOCA::Epetra::CompactWYOp::~CompactWYOp()

Destructor. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::SetUseTranspose "int LOCA::Epetra::CompactWYOp::SetUseTranspose(bool UseTranspose)

The operator currently does not support a transpose.

Setting this to true throws an error. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::Apply "int
LOCA::Epetra::CompactWYOp::Apply(const Epetra_MultiVector &Input,
Epetra_MultiVector &Result) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector Input in Result as described above. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::ApplyInverse "int
LOCA::Epetra::CompactWYOp::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

This method does nothing. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::NormInf "double
LOCA::Epetra::CompactWYOp::NormInf() const

Returns an approximate infinity norm of the operator matrix.

This is defined only if NormInf() of the underlying operator $J$ is
defined and is given by $\\\\|J\\\\|_\\\\infty+\\\\|A\\\\|_\\\\infty$.
";

%feature("docstring")  LOCA::Epetra::CompactWYOp::Label "const char *
LOCA::Epetra::CompactWYOp::Label() const

Returns a character std::string describing the operator. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::UseTranspose "bool
LOCA::Epetra::CompactWYOp::UseTranspose() const

Returns the current UseTranspose setting. Always returns false. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::HasNormInf "bool
LOCA::Epetra::CompactWYOp::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::Comm "const
Epetra_Comm & LOCA::Epetra::CompactWYOp::Comm() const

Returns a reference to the Epetra_Comm communicator associated with
this operator. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::OperatorDomainMap "const Epetra_Map & LOCA::Epetra::CompactWYOp::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::OperatorRangeMap "const Epetra_Map & LOCA::Epetra::CompactWYOp::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::init "void
LOCA::Epetra::CompactWYOp::init(const Epetra_MultiVector &x)

Initialize operator. Call this before starting a linear solve. The
Epetra_MultiVector argument x must be of the same size and
distribution as arguments to Apply(). ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::finish "void
LOCA::Epetra::CompactWYOp::finish()

Finish up solve. Call this after a linear solve is finished to inform
the operator that the solve is completed. ";

%feature("docstring")  LOCA::Epetra::CompactWYOp::applyCompactWY "void LOCA::Epetra::CompactWYOp::applyCompactWY(const
Epetra_MultiVector &x, Epetra_MultiVector &result_x,
Epetra_MultiVector &result_p) const

Applies the operator Q with a zero parameter component on input. ";


// File: classLOCA_1_1Hopf_1_1ComplexMultiVector.xml
%feature("docstring") LOCA::Hopf::ComplexMultiVector "

Multi-vector class to hold two multi-vectors to represent a complex
multi-vector.

This is not a true complex multi-vector. Operations like dot() and
multiply() are not correct for complex vectors. This class exists to
make some aspects of the real-equivalent formulation of complex linear
algebra simpler to implement.

C++ includes: LOCA_Hopf_ComplexMultiVector.H ";

%feature("docstring")
LOCA::Hopf::ComplexMultiVector::ComplexMultiVector "LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const NOX::Abstract::Vector
&cloneVec, int nColumns)

Constructor.

Generates a multivector with nColumns columns from cloneVec ";

%feature("docstring")
LOCA::Hopf::ComplexMultiVector::ComplexMultiVector "LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const NOX::Abstract::MultiVector
&realVec, const NOX::Abstract::MultiVector &imagVec)

Constructor. ";

%feature("docstring")
LOCA::Hopf::ComplexMultiVector::ComplexMultiVector "LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(const
ComplexMultiVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::ComplexMultiVector::ComplexMultiVector "LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(const
ComplexMultiVector &source, int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")
LOCA::Hopf::ComplexMultiVector::ComplexMultiVector "LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(const
ComplexMultiVector &source, const std::vector< int > &index, bool
view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")
LOCA::Hopf::ComplexMultiVector::~ComplexMultiVector "LOCA::Hopf::ComplexMultiVector::~ComplexMultiVector()

Destructor. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new multi-vector of the same underlying type by cloning
\"this\", and return a pointer to the new vector. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::clone(int numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::subCopy "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::subCopy(const std::vector< int >
&index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::subView "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::subView(const std::vector< int >
&index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::getRealMultiVec
"Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::getRealMultiVec() const

Returns the real component of extended multivector. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::getRealMultiVec
"Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::getRealMultiVec()

Returns the real component of extended multivector. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::getImagMultiVec
"Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::getImagMultiVec() const

Returns the imaginary component of extended multivector. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::getImagMultiVec
"Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::ComplexMultiVector::getImagMultiVec()

Returns the imaginary component of extended multivector. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::getColumn "Teuchos::RCP< LOCA::Hopf::ComplexVector >
LOCA::Hopf::ComplexMultiVector::getColumn(int i)

Returns ith column as an extended vector. ";

%feature("docstring")  LOCA::Hopf::ComplexMultiVector::getColumn "Teuchos::RCP< const LOCA::Hopf::ComplexVector >
LOCA::Hopf::ComplexMultiVector::getColumn(int i) const

Returns ith column as an extended vector. ";


// File: classLOCA_1_1BorderedSolver_1_1ComplexOperator.xml
%feature("docstring") LOCA::BorderedSolver::ComplexOperator "

Bordered solver operator representing the $J + i\\\\omega M$ as
implemented in the LOCA::Hopf::MooreSpence::AbstractGroup.

C++ includes: LOCA_BorderedSolver_ComplexOperator.H ";

%feature("docstring")
LOCA::BorderedSolver::ComplexOperator::ComplexOperator "LOCA::BorderedSolver::ComplexOperator::ComplexOperator(const
Teuchos::RCP< const LOCA::Hopf::MooreSpence::AbstractGroup > &grp,
double omega)

Constructor. ";

%feature("docstring")
LOCA::BorderedSolver::ComplexOperator::~ComplexOperator "LOCA::BorderedSolver::ComplexOperator::~ComplexOperator()

Destructor. ";

%feature("docstring")  LOCA::BorderedSolver::ComplexOperator::getGroup
"Teuchos::RCP< const NOX::Abstract::Group >
LOCA::BorderedSolver::ComplexOperator::getGroup() const

Get group pointer. ";

%feature("docstring")
LOCA::BorderedSolver::ComplexOperator::getFrequency "double
LOCA::BorderedSolver::ComplexOperator::getFrequency() const

Get frequency. ";

%feature("docstring")  LOCA::BorderedSolver::ComplexOperator::apply "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::apply(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const

Apply the operator. ";

%feature("docstring")
LOCA::BorderedSolver::ComplexOperator::applyTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::applyTranspose(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const

Apply transpose of the operator.

Group must be of type LOCA::Hopf::MinimallyAugmented::AbstractGroup
for this method to be implemented. ";

%feature("docstring")
LOCA::BorderedSolver::ComplexOperator::applyInverse "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const

Apply inverse of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::ComplexOperator::applyInverseTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const

Apply inverse transpose of the operator.

Group must be of type LOCA::Hopf::MinimallyAugmented::AbstractGroup
for this method to be implemented. ";


// File: classLOCA_1_1Hopf_1_1ComplexVector.xml
%feature("docstring") LOCA::Hopf::ComplexVector "

Vector class to hold two vectors to represent a complex vector.

This is not a true complex vector. Operations like innerProduct() are
not correct for complex vectors. This class exists to make some
aspects of the real-equivalent formulation of complex linear algebra
simpler to implement.

C++ includes: LOCA_Hopf_ComplexVector.H ";

%feature("docstring")  LOCA::Hopf::ComplexVector::ComplexVector "LOCA::Hopf::ComplexVector::ComplexVector(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const NOX::Abstract::Vector &realVec,
const NOX::Abstract::Vector &imagVec)

Constructor. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::ComplexVector "LOCA::Hopf::ComplexVector::ComplexVector(const ComplexVector &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::~ComplexVector "LOCA::Hopf::ComplexVector::~ComplexVector()

Destructor. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::clone "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Hopf::ComplexVector::clone(NOX::CopyType type=NOX::DeepCopy)
const

Cloning function. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::setVec "void
LOCA::Hopf::ComplexVector::setVec(const NOX::Abstract::Vector
&realVec, const NOX::Abstract::Vector &imagVec)

Sets the vector by setting its components. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::getRealVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::ComplexVector::getRealVec() const

Returns the real component of extended vector. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::getImagVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::ComplexVector::getImagVec() const

Returns the imaginary component of extended vector. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::getRealVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Hopf::ComplexVector::getRealVec()

Returns the real component of extended vector. ";

%feature("docstring")  LOCA::Hopf::ComplexVector::getImagVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Hopf::ComplexVector::getImagVec()

Returns the imaginary component of extended vector. ";


// File: classLOCA_1_1MultiContinuation_1_1CompositeConstraint.xml
%feature("docstring") LOCA::MultiContinuation::CompositeConstraint "

Implementation of LOCA::MultiContinuation::ConstraintInterface for
composite constraints, i.e., a constraint comprised of multiple,
separate constraints.

C++ includes: LOCA_MultiContinuation_CompositeConstraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::copy "void
LOCA::MultiContinuation::CompositeConstraint::copy(const
ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::clone "Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface >
LOCA::MultiContinuation::CompositeConstraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::numConstraints "int
LOCA::MultiContinuation::CompositeConstraint::numConstraints() const

Return number of constraints. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::setX "void
LOCA::MultiContinuation::CompositeConstraint::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::setParam "void
LOCA::MultiContinuation::CompositeConstraint::setParam(int paramID,
double val)

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::setParams "void
LOCA::MultiContinuation::CompositeConstraint::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::computeConstraints "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::isConstraints "bool
LOCA::MultiContinuation::CompositeConstraint::isConstraints() const

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::isDX "bool
LOCA::MultiContinuation::CompositeConstraint::isDX() const

Return true if derivatives of constraints w.r.t. x are valid. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::getConstraints "const
NOX::Abstract::MultiVector::DenseMatrix &
LOCA::MultiContinuation::CompositeConstraint::getConstraints() const

Return constraint residuals. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::multiplyDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::multiplyDX(double alpha,
const NOX::Abstract::MultiVector &input_x,
NOX::Abstract::MultiVector::DenseMatrix &result_p) const

Compute result_p = alpha * dg/dx * input_x.

Note that if there are n constraints and input_x has m columns,
result_p should be a n by m matrix and is equivalent to ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::addDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraint::addDX(Teuchos::ETransp
transb, double alpha, const NOX::Abstract::MultiVector::DenseMatrix
&b, double beta, NOX::Abstract::MultiVector &result_x) const

Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x.

Note that this should be equivalent to ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::isDXZero "bool
LOCA::MultiContinuation::CompositeConstraint::isDXZero() const

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::preProcessContinuationStep
"void
LOCA::MultiContinuation::CompositeConstraint::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. The default implementation is empty. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::postProcessContinuationStep
"void
LOCA::MultiContinuation::CompositeConstraint::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. The
default implementation is empty. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint "LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const std::vector<
Teuchos::RCP< LOCA::MultiContinuation::ConstraintInterface > >
&constraintObjects)

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint "LOCA::MultiContinuation::CompositeConstraint::CompositeConstraint(const
CompositeConstraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraint::~CompositeConstraint "LOCA::MultiContinuation::CompositeConstraint::~CompositeConstraint()

Destructor. ";


// File: classLOCA_1_1MultiContinuation_1_1CompositeConstraintMVDX.xml
%feature("docstring") LOCA::MultiContinuation::CompositeConstraintMVDX
"

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
composite constraints, i.e., a constraint comprised of multiple,
separate constraints.

C++ includes: LOCA_MultiContinuation_CompositeConstraintMVDX.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::copy "void
LOCA::MultiContinuation::CompositeConstraintMVDX::copy(const
ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::clone "Teuchos::RCP< LOCA::MultiContinuation::ConstraintInterface >
LOCA::MultiContinuation::CompositeConstraintMVDX::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::computeDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraintMVDX::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::multiplyDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraintMVDX::multiplyDX(double
alpha, const NOX::Abstract::MultiVector &input_x,
NOX::Abstract::MultiVector::DenseMatrix &result_p) const

Compute result_p = alpha * dg/dx * input_x.

Note that if there are n constraints and input_x has m columns,
result_p should be a n by m matrix and is equivalent to ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::addDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraintMVDX::addDX(Teuchos::ETransp
transb, double alpha, const NOX::Abstract::MultiVector::DenseMatrix
&b, double beta, NOX::Abstract::MultiVector &result_x) const

Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x.

Note that this should be equivalent to ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX
*/

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::getDX "const
NOX::Abstract::MultiVector *
LOCA::MultiContinuation::CompositeConstraintMVDX::getDX() const

Return solution component of constraint derivatives. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::CompositeConstraintMVDX
"LOCA::MultiContinuation::CompositeConstraintMVDX::CompositeConstraintMVDX(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const std::vector<
Teuchos::RCP< LOCA::MultiContinuation::ConstraintInterfaceMVDX > >
&constraintObjects)

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::CompositeConstraintMVDX
"LOCA::MultiContinuation::CompositeConstraintMVDX::CompositeConstraintMVDX(const
CompositeConstraintMVDX &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::CompositeConstraintMVDX::~CompositeConstraintMVDX
"LOCA::MultiContinuation::CompositeConstraintMVDX::~CompositeConstraintMVDX()

Destructor. ";


// File: classLOCA_1_1StepSize_1_1Constant.xml
%feature("docstring") LOCA::StepSize::Constant "

Constant step size control strategy

This class implements a roughly constant step size control strategy.
If the previous step was sucessful, the new step size is set equal to
the old, otherwise the step size is cut by a supplied factor. Once a
sucessful step is made, the step size is increased by a supplied
factor until the initial step size is reached.

This class also incorporates rescaling of the continuation parameter
when calculating a step size (common in arc-length continuation). For
the first continuation step, the step size is chosen so that step size
times the parameter component of the predictor is equal to the initial
step size. From then on, the step size is multiplied by the step size
scale factor (see ( LOCA::MultiContinuation::ArcLengthGroup) which
incorporates rescaling of the continuation parameter.

The parameters used by this class supplied in the constructor are:
\"Max Step Size\" - Largest valid step size (Default 1.0e+12)

\"Min Step Size\" - Smallest valid step size (Default 1.0e-12)

\"Initial Step Size\" - Initial step size (Default 1.0)

\"Failed Step Reduction Factor\" - Factor by which step size is
reduced after a failed step. (Default 0.5)

\"Successful Step Increase Factor\" - Factor by which step size is
increased after a successful step when the step size is smaller than
the initial step size. (Default 1.26 = $2^{1/3}$)

C++ includes: LOCA_StepSize_Constant.H ";

%feature("docstring")  LOCA::StepSize::Constant::Constant "LOCA::StepSize::Constant::Constant(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepsizeParams)

Constructor. ";

%feature("docstring")  LOCA::StepSize::Constant::~Constant "LOCA::StepSize::Constant::~Constant()

Destructor. ";

%feature("docstring")  LOCA::StepSize::Constant::computeStepSize "NOX::Abstract::Group::ReturnType
LOCA::StepSize::Constant::computeStepSize(LOCA::MultiContinuation::AbstractStrategy
&curGroup, const LOCA::MultiContinuation::ExtendedVector &predictor,
const NOX::Solver::Generic &solver, const
LOCA::Abstract::Iterator::StepStatus &stepStatus, const
LOCA::Abstract::Iterator &stepper, double &stepSize)

Compute the step size as described above.

Parameters:
-----------

curGroup:  [in] Current continuation group

predictor:  [in] Current predictor direction

solver:  [in] Solver from previous step

stepStatus:  [in] Status of previous step

stepper:  [in] Stepper

stepSize:  [out] Computed step size

ReturnType code indicating success or failure Returns
NOX::Abstract::Group::Failed if the computed step size is smaller than
the minimum step size ";

%feature("docstring")  LOCA::StepSize::Constant::getPrevStepSize "double LOCA::StepSize::Constant::getPrevStepSize() const

Returns previous step size. ";

%feature("docstring")  LOCA::StepSize::Constant::getStartStepSize "double LOCA::StepSize::Constant::getStartStepSize() const

Returns initial step size. ";


// File: classLOCA_1_1MultiPredictor_1_1Constant.xml
%feature("docstring") LOCA::MultiPredictor::Constant "

Constant predictor strategy

This class computes the predictor direction given by a vector of zeros
for the solution vector component and 1 for the parameter component.
When used with natural continuation, this corresponds to what is
commonly referred to as zero'th order continuation.

C++ includes: LOCA_MultiPredictor_Constant.H ";

%feature("docstring")  LOCA::MultiPredictor::Constant::Constant "LOCA::MultiPredictor::Constant::Constant(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &predParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

predParams:  [in] Predictor parameters. None are currently referenced.
";

%feature("docstring")  LOCA::MultiPredictor::Constant::~Constant "LOCA::MultiPredictor::Constant::~Constant()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::Constant::Constant "LOCA::MultiPredictor::Constant::Constant(const Constant &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::MultiPredictor::Constant::clone "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::MultiPredictor::Constant::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. ";

%feature("docstring")  LOCA::MultiPredictor::Constant::compute "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Constant::compute(bool baseOnSecant, const
std::vector< double > &stepSize,
LOCA::MultiContinuation::ExtendedGroup &grp, const
LOCA::MultiContinuation::ExtendedVector &prevXVec, const
LOCA::MultiContinuation::ExtendedVector &xVec)

Compute the predictor given the current and previous solution vectors.
Set baseOnSecant to false if the predictor orientation should not be
based on the secant vector (first or last steps of a continuation
run).

This method actually implements the predictor computation described
above ";

%feature("docstring")  LOCA::MultiPredictor::Constant::evaluate "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Constant::evaluate(const std::vector< double >
&stepSize, const LOCA::MultiContinuation::ExtendedVector &xVec,
LOCA::MultiContinuation::ExtendedMultiVector &result) const

Evaluate predictor with step size stepSize.

This method computes result[i] = xVec[i] + stepSize[i] * v[i] for each
i, where v[i] is the ith predictor direction. ";

%feature("docstring")  LOCA::MultiPredictor::Constant::computeTangent
"NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Constant::computeTangent(LOCA::MultiContinuation::ExtendedMultiVector
&tangent)

Compute tangent to predictor and store in tangent. ";

%feature("docstring")
LOCA::MultiPredictor::Constant::isTangentScalable "bool
LOCA::MultiPredictor::Constant::isTangentScalable() const

Is the tangent vector for this predictor scalable.

For the constant predictor, this always returns false. ";


// File: classLOCA_1_1MultiContinuation_1_1ConstrainedGroup.xml
%feature("docstring") LOCA::MultiContinuation::ConstrainedGroup "

Extended group representing a constrained nonlinear problem.

This class represents a constrained system of nonlinear equations:
\\\\[ \\\\begin{split} f(x,p) &= 0 \\\\\\\\ g(x,p) &= 0 \\\\end{split}
\\\\] where $x\\\\in\\\\Re^n$ is the solution vector,
$p\\\\in\\\\Re^m$ is a set of constraint parameters,
$f(x,p)\\\\in\\\\Re^n$ is represented by some
LOCA::MultiContinuation::AbstractGroup, and $g(x,p)\\\\in\\\\Re^m$ is
a constraint represented by a
LOCA::MultiContinuation::ConstraintInterface object. Newton steps for
this system are computed via some
LOCA::BorderedSolver::AbstractStrategy which is specified via the
constraintParams argument to the constructor.

C++ includes: LOCA_MultiContinuation_ConstrainedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::clone "Teuchos::RCP<
NOX::Abstract::Group >
LOCA::MultiContinuation::ConstrainedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. ";

%feature("docstring")  LOCA::MultiContinuation::ConstrainedGroup::setX
"void LOCA::MultiContinuation::ConstrainedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeX "void
LOCA::MultiContinuation::ConstrainedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute and return solution vector, x, where this.x = grp.x + step *
d. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeF()

Compute extended continuation equations. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeJacobian()

Compute extended continuation jacobian. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeGradient()

Gradient is not defined for this system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction for extended continuation system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianMultiVector "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")  LOCA::MultiContinuation::ConstrainedGroup::isF
"bool LOCA::MultiContinuation::ConstrainedGroup::isF() const

Return true if extended residual is valid. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::isJacobian "bool
LOCA::MultiContinuation::ConstrainedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::isGradient "bool
LOCA::MultiContinuation::ConstrainedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::isNewton "bool
LOCA::MultiContinuation::ConstrainedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")  LOCA::MultiContinuation::ConstrainedGroup::getX
"const NOX::Abstract::Vector &
LOCA::MultiContinuation::ConstrainedGroup::getX() const

Return extended solution vector. ";

%feature("docstring")  LOCA::MultiContinuation::ConstrainedGroup::getF
"const NOX::Abstract::Vector &
LOCA::MultiContinuation::ConstrainedGroup::getF() const

Return extended residual. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getNormF "double
LOCA::MultiContinuation::ConstrainedGroup::getNormF() const

Return 2-norm of extended residual. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::MultiContinuation::ConstrainedGroup::getGradient() const

Gradient is never valid. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::MultiContinuation::ConstrainedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getXPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getXPtr() const

Return RCP to extended solution vector. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getFPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getFPtr() const

Return RCP to extended residual. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getGradientPtr() const

Gradient is never valid. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ConstrainedGroup::getNewtonPtr() const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getNormNewtonSolveResidual
"double
LOCA::MultiContinuation::ConstrainedGroup::getNormNewtonSolveResidual()
const

Returns 2-norm of extended Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup() const

Return underlying group. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::MultiContinuation::ConstrainedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")  LOCA::MultiContinuation::ConstrainedGroup::copy
"void LOCA::MultiContinuation::ConstrainedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::setParamsMulti "void
LOCA::MultiContinuation::ConstrainedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::setParams "void
LOCA::MultiContinuation::ConstrainedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p (pVector = p). ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::setParam "void
LOCA::MultiContinuation::ConstrainedGroup::setParam(int paramID,
double val)

Set parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::setParam "void
LOCA::MultiContinuation::ConstrainedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getParams "const
LOCA::ParameterVector &
LOCA::MultiContinuation::ConstrainedGroup::getParams() const

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getParam "double
LOCA::MultiContinuation::ConstrainedGroup::getParam(int paramID) const

Return copy of parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getParam "double
LOCA::MultiContinuation::ConstrainedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::preProcessContinuationStep
"void
LOCA::MultiContinuation::ConstrainedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::postProcessContinuationStep
"void
LOCA::MultiContinuation::ConstrainedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::projectToDraw "void
LOCA::MultiContinuation::ConstrainedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::projectToDrawDimension "int
LOCA::MultiContinuation::ConstrainedGroup::projectToDrawDimension()
const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::computeScaledDotProduct "double
LOCA::MultiContinuation::ConstrainedGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute a scaled dot product. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::printSolution "void
LOCA::MultiContinuation::ConstrainedGroup::printSolution(const double
conParam) const

Function to print out solution and parameter after successful step. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::printSolution "void
LOCA::MultiContinuation::ConstrainedGroup::printSolution(const
NOX::Abstract::Vector &x, const double conParam) const

Function to print out a vector and parameter after successful step. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::scaleVector "void
LOCA::MultiContinuation::ConstrainedGroup::scaleVector(NOX::Abstract::Vector
&x) const

Scales a vector using scaling vector. ";

/*  Implementation of  */

/*  LOCA::BorderedSystem::AbstractGroup virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getBorderedWidth "int
LOCA::MultiContinuation::ConstrainedGroup::getBorderedWidth() const

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getUnborderedGroup "Teuchos::RCP< const NOX::Abstract::Group >
LOCA::MultiContinuation::ConstrainedGroup::getUnborderedGroup() const

Get bottom-level unbordered group. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::isCombinedAZero "bool
LOCA::MultiContinuation::ConstrainedGroup::isCombinedAZero() const

Indicates whether combined A block is zero. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::isCombinedBZero "bool
LOCA::MultiContinuation::ConstrainedGroup::isCombinedBZero() const

Indicates whether combined B block is zero. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::isCombinedCZero "bool
LOCA::MultiContinuation::ConstrainedGroup::isCombinedCZero() const

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::extractSolutionComponent "void
LOCA::MultiContinuation::ConstrainedGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::extractParameterComponent "void
LOCA::MultiContinuation::ConstrainedGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::loadNestedComponents "void
LOCA::MultiContinuation::ConstrainedGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::fillA "void
LOCA::MultiContinuation::ConstrainedGroup::fillA(NOX::Abstract::MultiVector
&A) const

Fill the combined A block as described above. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::fillB "void
LOCA::MultiContinuation::ConstrainedGroup::fillB(NOX::Abstract::MultiVector
&B) const

Fill the combined B block as described above. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::fillC "void
LOCA::MultiContinuation::ConstrainedGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const

Fill the combined C block as described above. ";

/*  Implementation of LOCA::Abstract::TransposeSolveGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverse
"NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Solve Jacobian-tranpose system. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstrainedGroup::applyJacobianTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Solve Jacobian-tranpose system with multiple right-hand sides. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup "LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &constraintParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, const Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface > &constraints, const
std::vector< int > &paramIDs, bool skip_dfdp=false)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list.

constraintParams:  [in] Parameter list determining the bordered solver
method.

grp:  [in] Group representing $f$.

constraints:  [in] Constraint object representing $g$.

paramIDs:  [in] Parameter IDs of the constraint parameters

skip_dfdp:  [in] Whether to skip computation of df/dp when computing
extended Jacobian. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup "LOCA::MultiContinuation::ConstrainedGroup::ConstrainedGroup(const
ConstrainedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::~ConstrainedGroup "LOCA::MultiContinuation::ConstrainedGroup::~ConstrainedGroup()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::setConstraintParameter "void
LOCA::MultiContinuation::ConstrainedGroup::setConstraintParameter(int
i, double val)

Set constraint parameter i to value val. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParameter "double
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParameter(int
i) const

Get constraint parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParamIDs "const std::vector< int > &
LOCA::MultiContinuation::ConstrainedGroup::getConstraintParamIDs()
const

Get constraint parameter IDs. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getGroup "Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup >
LOCA::MultiContinuation::ConstrainedGroup::getGroup()

Get group. ";

%feature("docstring")
LOCA::MultiContinuation::ConstrainedGroup::getConstraints "Teuchos::RCP< LOCA::MultiContinuation::ConstraintInterface >
LOCA::MultiContinuation::ConstrainedGroup::getConstraints()

Get constraints. ";


// File: classLOCA_1_1TurningPoint_1_1MinimallyAugmented_1_1Constraint.xml
%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint "

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
computing turning points for the minimally augmented turning point
formulation.

This class implements the turning point constraint equation
$\\\\sigma(x,p) = 0$ for the minimally augmented turning point
formulation where $\\\\sigma$ is defined via \\\\[ \\\\begin{bmatrix}
J & a \\\\\\\\ b^T & 0 \\\\end{bmatrix} \\\\begin{bmatrix} v \\\\\\\\
\\\\sigma_1 \\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n
\\\\end{bmatrix}, \\\\] \\\\[ \\\\begin{bmatrix} J^T & b \\\\\\\\ a^T
& 0 \\\\end{bmatrix} \\\\begin{bmatrix} w \\\\\\\\ \\\\sigma_2
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\sigma = -w^T J v/n \\\\] for any vectors $a$ and $b$
in $\\\\Re^n$. Using these relationships, it is easy to show \\\\[
\\\\begin{split} \\\\sigma_x &= -(w^T J v)_x/n = -w^T J_x v/n \\\\\\\\
\\\\sigma_p &= -(w^T J v)_p/n = -w^T J_p v/n \\\\end{split} \\\\]

The class is intialized via the tpParams parameter list argument to
the constructor. The parameters this class recognizes are: \"Symmetric
Jacobian\" -- [bool] (default: false) - Flag indicating whether
Jacobian matrix $J$ is symmetric, in which case we force $a = b$ and
therefore the second tranpose solve for $w$ is unnecessary

\"Initial Null Vector Compuation\" -- [string] (default: \"User
Provided\") - Method to compute initial $a$ and $b$ vectors. Valid
choices are: \"User Provided\" - Initial vectors are provided in the
parameter list, in which case the following parameters are relevant:
\"Initial A Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be
supplied) - Vector storing initial value for $a$ vector

\"Initial B Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be
supplied for nonsymmetric Jacobians) - Vector storing initial value
for $b$ vector

\"Solve df/dp\" - Compute $a = J^{-T}df/dp$ and $b = J^{-1} df/dp$
where $p$ is the bifurcation parameter.

\"Constant\" - Entries of $a$ and $b$ are set to 1.0

\"Null Vector Scaling\" -- [string] (default: \"Order N\") - Method to
scale $a$ and $b$. This determines the norm of these vectors and the
scaling of $\\\\sigma$. Valid choices are: \"None\" -- Use initial
scaling

\"Order 1\" -- Scale to unit norm

\"Order N\" -- Use vector length scaling

\"Update Null Vectors Every Continuation Step\" -- [bool] (default:
true) - Flag indicating whether to update $a$ and $b$ vectors via $a =
w$ and $b = v$ every continuation step

\"Update Null Vectors Every Nonlinear Iteration\" -- [bool] (default:
false) - Flag indicating whether to update $a$ and $b$ vectors via $a
= w$ and $b = v$ every nonlinear iteration

\"Multiply Null Vectors by Mass Matrix\" -- [bool] (default: false) -
Flag indicating whether to multiply $a$ and $b$ vectors by the mass
matrix $M = \\\\partial f/\\\\partial\\\\dot{x}$ at the strart of a
turning point calculation, and each time $a$ and $b$ are updated. This
can improve the scaling of these vectors, and may orthogonalize them
against structural null spaces (i.e., pressure null space for
incompressible Navier-Stokes).

C++ includes: LOCA_TurningPoint_MinimallyAugmented_Constraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::copy "void
LOCA::TurningPoint::MinimallyAugmented::Constraint::copy(const
LOCA::MultiContinuation::ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::clone "Teuchos::RCP< LOCA::MultiContinuation::ConstraintInterface >
LOCA::TurningPoint::MinimallyAugmented::Constraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::numConstraints "int
LOCA::TurningPoint::MinimallyAugmented::Constraint::numConstraints()
const

Return number of constraints. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::setX "void
LOCA::TurningPoint::MinimallyAugmented::Constraint::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::setParam "void
LOCA::TurningPoint::MinimallyAugmented::Constraint::setParam(int
paramID, double val)

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::setParams "void
LOCA::TurningPoint::MinimallyAugmented::Constraint::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::computeConstraints
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::isConstraints "bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::isConstraints()
const

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::isDX "bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::isDX() const

Return true if derivatives of constraints w.r.t. x are valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getConstraints "const NOX::Abstract::MultiVector::DenseMatrix &
LOCA::TurningPoint::MinimallyAugmented::Constraint::getConstraints()
const

Return constraint residuals. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getDX "const
NOX::Abstract::MultiVector *
LOCA::TurningPoint::MinimallyAugmented::Constraint::getDX() const

Return solution component of constraint derivatives. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::isDXZero "bool
LOCA::TurningPoint::MinimallyAugmented::Constraint::isDXZero() const

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::postProcessContinuationStep
"void
LOCA::TurningPoint::MinimallyAugmented::Constraint::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful.
Here we update the $a$ and $b$ vectors to $w$ and $v$ respectively if
requested. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::Constraint "LOCA::TurningPoint::MinimallyAugmented::Constraint::Constraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &tpParams, const Teuchos::RCP<
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup > &g, int
bif_param)

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::Constraint "LOCA::TurningPoint::MinimallyAugmented::Constraint::Constraint(const
Constraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::~Constraint "LOCA::TurningPoint::MinimallyAugmented::Constraint::~Constraint()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::setGroup "void
LOCA::TurningPoint::MinimallyAugmented::Constraint::setGroup(const
Teuchos::RCP< LOCA::TurningPoint::MinimallyAugmented::AbstractGroup >
&g)

Set the group pointer.

This method should be called when ever the constrained group is
copied, since we don't explicitly copy the underlying group here. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getLeftNullVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::Constraint::getLeftNullVec()
const

Returns left null vector w. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getRightNullVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::Constraint::getRightNullVec()
const

Returns right null vector v. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getAVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::Constraint::getAVec() const

Returns a vector. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getBVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::Constraint::getBVec() const

Returns b vector. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::Constraint::getSigma "double
LOCA::TurningPoint::MinimallyAugmented::Constraint::getSigma() const

Returns sigma. ";


// File: classLOCA_1_1Pitchfork_1_1MinimallyAugmented_1_1Constraint.xml
%feature("docstring") LOCA::Pitchfork::MinimallyAugmented::Constraint
"

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
computing pitchforks for the minimally augmented pitchfork
formulation.

This class implements the pitchfork constraint equations
$\\\\sigma(x,p) = 0$, $\\\\langle \\\\psi,x \\\\rangle = 0$ for the
minimally augmented pitchfork formulation where $\\\\sigma$ is defined
via \\\\[ \\\\begin{bmatrix} J & a \\\\\\\\ b^T & 0 \\\\end{bmatrix}
\\\\begin{bmatrix} v \\\\\\\\ \\\\sigma_1 \\\\end{bmatrix} =
\\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix}, \\\\] \\\\[
\\\\begin{bmatrix} J^T & b \\\\\\\\ a^T & 0 \\\\end{bmatrix}
\\\\begin{bmatrix} w \\\\\\\\ \\\\sigma_2 \\\\end{bmatrix} =
\\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix}, \\\\] \\\\[
\\\\sigma = -w^T J v/n \\\\] for any vectors $a$ and $b$ in
$\\\\Re^n$. Using these relationships, it is easy to show \\\\[
\\\\begin{split} \\\\sigma_x &= -(w^T J v)_x/n = -w^T J_x v/n \\\\\\\\
\\\\sigma_p &= -(w^T J v)_p/n = -w^T J_p v/n \\\\end{split} \\\\]

The class is derived from
LOCA::TurningPoint::MinimallyAugmented::Constraint. See this class for
a description of available parameters.

C++ includes: LOCA_Pitchfork_MinimallyAugmented_Constraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::copy "void
LOCA::Pitchfork::MinimallyAugmented::Constraint::copy(const
LOCA::MultiContinuation::ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::clone "Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface >
LOCA::Pitchfork::MinimallyAugmented::Constraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::numConstraints "int
LOCA::Pitchfork::MinimallyAugmented::Constraint::numConstraints()
const

Return number of constraints. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::computeConstraints "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::Constraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::Constraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::Constraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::getConstraints "const NOX::Abstract::MultiVector::DenseMatrix &
LOCA::Pitchfork::MinimallyAugmented::Constraint::getConstraints()
const

Return constraint residuals. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::getDX "const
NOX::Abstract::MultiVector *
LOCA::Pitchfork::MinimallyAugmented::Constraint::getDX() const

Return solution component of constraint derivatives. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::Constraint "LOCA::Pitchfork::MinimallyAugmented::Constraint::Constraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &pfParams, const Teuchos::RCP<
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup > &g, const
Teuchos::RCP< const NOX::Abstract::Vector > &psi, int bif_param)

Constructor. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::Constraint "LOCA::Pitchfork::MinimallyAugmented::Constraint::Constraint(const
Constraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::~Constraint "LOCA::Pitchfork::MinimallyAugmented::Constraint::~Constraint()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::Constraint::setGroup "void
LOCA::Pitchfork::MinimallyAugmented::Constraint::setGroup(const
Teuchos::RCP< LOCA::TurningPoint::MinimallyAugmented::AbstractGroup >
&g)

Set the group pointer.

This method should be called when ever the constrained group is
copied, since we don't explicitly copy the underlying group here. ";


// File: classLOCA_1_1Hopf_1_1MinimallyAugmented_1_1Constraint.xml
%feature("docstring") LOCA::Hopf::MinimallyAugmented::Constraint "

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
computing Hopf bifurcations for the minimally augmented Hopf
formulation.

This class implements the turning point constraint equation
$\\\\sigma(x,p,\\\\omega) = 0$ for the minimally augmented Hopf
formulation where $\\\\sigma$ is defined via \\\\[ \\\\begin{bmatrix}
J+i\\\\omega M & a \\\\\\\\ b^H & 0 \\\\end{bmatrix}
\\\\begin{bmatrix} v \\\\\\\\ \\\\sigma_1 \\\\end{bmatrix} =
\\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix}, \\\\] \\\\[
\\\\begin{bmatrix} J^T-i\\\\omega M^T & b \\\\\\\\ a^H & 0
\\\\end{bmatrix} \\\\begin{bmatrix} w \\\\\\\\ \\\\sigma_2
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\sigma = w^H J+i\\\\omega M v/n \\\\] for any vectors
$a$ and $b$ in $\\\\mathbb{C}^n$. Using these relationships, it is
easy to show \\\\[ \\\\begin{split} \\\\sigma_x &= (w^H(J+i\\\\omega
M)v)_x/n = w^H(J+i\\\\omega M)_x v/n \\\\\\\\ \\\\sigma_p &=
(w^H(J+i\\\\omega M)v)_p/n = w^H(J+i\\\\omega M)_p v/n \\\\end{split}
\\\\]

The class is intialized via the hpfParams parameter list argument to
the constructor. The parameters this class recognizes are: \"Update
Null Vectors Every Continuation Step\" -- [bool] (default: true) -
Flag indicating whether to update $a$ and $b$ vectors via $a = w$ and
$b = v$ every continuation step

\"Update Null Vectors Every Nonlinear Iteration\" -- [bool] (default:
false) - Flag indicating whether to update $a$ and $b$ vectors via $a
= w$ and $b = v$ every nonlinear iteration

C++ includes: LOCA_Hopf_MinimallyAugmented_Constraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::copy "void
LOCA::Hopf::MinimallyAugmented::Constraint::copy(const
LOCA::MultiContinuation::ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::clone "Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface >
LOCA::Hopf::MinimallyAugmented::Constraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::numConstraints "int
LOCA::Hopf::MinimallyAugmented::Constraint::numConstraints() const

Return number of constraints. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::setX "void
LOCA::Hopf::MinimallyAugmented::Constraint::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::setParam "void
LOCA::Hopf::MinimallyAugmented::Constraint::setParam(int paramID,
double val)

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::setParams "void
LOCA::Hopf::MinimallyAugmented::Constraint::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::computeConstraints "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::isConstraints "bool
LOCA::Hopf::MinimallyAugmented::Constraint::isConstraints() const

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::isDX "bool
LOCA::Hopf::MinimallyAugmented::Constraint::isDX() const

Return true if derivatives of constraints w.r.t. x are valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getConstraints "const
NOX::Abstract::MultiVector::DenseMatrix &
LOCA::Hopf::MinimallyAugmented::Constraint::getConstraints() const

Return constraint residuals. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getDX "const
NOX::Abstract::MultiVector *
LOCA::Hopf::MinimallyAugmented::Constraint::getDX() const

Return solution component of constraint derivatives. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::isDXZero "bool
LOCA::Hopf::MinimallyAugmented::Constraint::isDXZero() const

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::postProcessContinuationStep
"void
LOCA::Hopf::MinimallyAugmented::Constraint::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful.
Here we update the $a$ and $b$ vectors to $w$ and $v$ respectively if
requested. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::Constraint "LOCA::Hopf::MinimallyAugmented::Constraint::Constraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &hpfParams, const Teuchos::RCP<
LOCA::Hopf::MinimallyAugmented::AbstractGroup > &g, bool is_symmetric,
const NOX::Abstract::Vector &a_real, const NOX::Abstract::Vector
&a_imag, const NOX::Abstract::Vector *b_real, const
NOX::Abstract::Vector *b_imag, int bif_param, double freq)

Constructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::Constraint "LOCA::Hopf::MinimallyAugmented::Constraint::Constraint(const
Constraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::~Constraint "LOCA::Hopf::MinimallyAugmented::Constraint::~Constraint()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::setGroup "void
LOCA::Hopf::MinimallyAugmented::Constraint::setGroup(const
Teuchos::RCP< LOCA::Hopf::MinimallyAugmented::AbstractGroup > &g)

Set the group pointer.

This method should be called when ever the constrained group is
copied, since we don't explicitly copy the underlying group here. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::setFrequency "void
LOCA::Hopf::MinimallyAugmented::Constraint::setFrequency(double freq)

Set Hopf frequency. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getLeftNullVecReal "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::Constraint::getLeftNullVecReal() const

Returns real component of left null vector w. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getLeftNullVecImag "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::Constraint::getLeftNullVecImag() const

Returns imaginary component of left null vector w. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getRightNullVecReal "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::Constraint::getRightNullVecReal()
const

Returns real component of right null vector v. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getRightNullVecImag "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::Constraint::getRightNullVecImag()
const

Returns imaginary component of right null vector v. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getSigmaReal "double
LOCA::Hopf::MinimallyAugmented::Constraint::getSigmaReal() const

Returns real component of sigma. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::getSigmaImag "double
LOCA::Hopf::MinimallyAugmented::Constraint::getSigmaImag() const

Returns imaginary component of sigma. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::Constraint::computeDOmega "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::Constraint::computeDOmega(NOX::Abstract::MultiVector::DenseMatrix
&domega)

Compute derivative of sigma w.r.t. frequency omega. ";


// File: classLOCA_1_1MultiContinuation_1_1ConstraintInterface.xml
%feature("docstring") LOCA::MultiContinuation::ConstraintInterface "

Abstract interface for the constraint portion of a constrained
nonlinear system.

This class is used in conjunction with
LOCA::MultiContinuation::ConstrainedGroup to represent a constrained
nonlinear system: \\\\[ f(x,y) = 0 g(x,y) = 0 \\\\] where $f(x,y)$ is
represented by a concrete implementation of a
LOCA::MultiContinuation::AbstractGroup and $g(x,y)$ (the constraint)
is represented by an implementation of this class. Here it is assumed
the resulting system is square, i.e., $x\\\\in\\\\Re^n$,
$y\\\\in\\\\Re^m$, $f(x,y)\\\\in\\\\Re^n$ and $g(x,y)\\\\in\\\\Re^m$.

This class provides an interface to evaluate $g(x,y)$, compute the
derivatives $g_x$ and $g_y$, and apply the derivative $g_x$ to
arbitrary multi-vectors (the implementation is never required to
explicitly store $g_x$ which is impractical in many situations).

C++ includes: LOCA_MultiContinuation_ConstraintInterface.H ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::ConstraintInterface "LOCA::MultiContinuation::ConstraintInterface::ConstraintInterface()

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::~ConstraintInterface "virtual
LOCA::MultiContinuation::ConstraintInterface::~ConstraintInterface()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::copy "virtual void
LOCA::MultiContinuation::ConstraintInterface::copy(const
ConstraintInterface &source)=0

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::clone "virtual
Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::ConstraintInterface::clone(NOX::CopyType
type=NOX::DeepCopy) const =0

Cloning function. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::numConstraints "virtual
int LOCA::MultiContinuation::ConstraintInterface::numConstraints()
const =0

Return number of constraints. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::setX "virtual void
LOCA::MultiContinuation::ConstraintInterface::setX(const
NOX::Abstract::Vector &x)=0

Set the solution vector to x. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::setParam "virtual void
LOCA::MultiContinuation::ConstraintInterface::setParam(int paramID,
double val)=0

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::setParams "virtual void
LOCA::MultiContinuation::ConstraintInterface::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)=0

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::computeConstraints "virtual NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterface::computeConstraints()=0

Compute constraint residuals. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::computeDX "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterface::computeDX()=0

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::computeDP "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterface::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)=0

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::isConstraints "virtual
bool LOCA::MultiContinuation::ConstraintInterface::isConstraints()
const =0

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::isDX "virtual bool
LOCA::MultiContinuation::ConstraintInterface::isDX() const =0

Return true if derivative of constraint w.r.t. x is valid. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::getConstraints "virtual
const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::MultiContinuation::ConstraintInterface::getConstraints() const
=0

Return constraint residuals. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::multiplyDX "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterface::multiplyDX(double alpha,
const NOX::Abstract::MultiVector &input_x,
NOX::Abstract::MultiVector::DenseMatrix &result_p) const =0

Compute result_p = alpha * dg/dx * input_x.

Note that if there are n constraints and input_x has m columns,
result_p should be a n by m matrix and is equivalent to ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::addDX "virtual
NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterface::addDX(Teuchos::ETransp
transb, double alpha, const NOX::Abstract::MultiVector::DenseMatrix
&b, double beta, NOX::Abstract::MultiVector &result_x) const =0

Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x.

Note that this should be equivalent to ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::isDXZero "virtual bool
LOCA::MultiContinuation::ConstraintInterface::isDXZero() const =0

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::preProcessContinuationStep
"virtual void
LOCA::MultiContinuation::ConstraintInterface::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. The default implementation is empty. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterface::postProcessContinuationStep
"virtual void
LOCA::MultiContinuation::ConstraintInterface::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. The
default implementation is empty. ";


// File: classLOCA_1_1MultiContinuation_1_1ConstraintInterfaceMVDX.xml
%feature("docstring") LOCA::MultiContinuation::ConstraintInterfaceMVDX
"

Abstract interface for the constraint portion of a constrained
nonlinear system for constraints that support computing a solution
component derivative as a multi-vector.

This class extends the LOCA::MultiContinuation::ConstraintInterface to
support constraints that support computing the entire derivative with
respect to the solution components (x) and storing the resulting
derivative as a multivector. This interface adds one additional
method, getConstraintDerivativesX(), that returns this derivative.
Additionally, it implements the applyConstraintDerivativesX() methods
using standard multi-vector operations.

C++ includes: LOCA_MultiContinuation_ConstraintInterfaceMVDX.H ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterfaceMVDX::ConstraintInterfaceMVDX
"LOCA::MultiContinuation::ConstraintInterfaceMVDX::ConstraintInterfaceMVDX()

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterfaceMVDX::~ConstraintInterfaceMVDX
"virtual
LOCA::MultiContinuation::ConstraintInterfaceMVDX::~ConstraintInterfaceMVDX()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterfaceMVDX::getDX "virtual
const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::ConstraintInterfaceMVDX::getDX() const =0

Return solution component of constraint derivatives.

May return NULL if constraint derivative is zero ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterfaceMVDX::multiplyDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterfaceMVDX::multiplyDX(double
alpha, const NOX::Abstract::MultiVector &input_x,
NOX::Abstract::MultiVector::DenseMatrix &result_p) const

Compute result_p = alpha * dg/dx * input_x.

This method is implemented using getConstraintDerivativesX() and the
NOX::Abstract::MultiVector::multiply() method. ";

%feature("docstring")
LOCA::MultiContinuation::ConstraintInterfaceMVDX::addDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterfaceMVDX::addDX(Teuchos::ETransp
transb, double alpha, const NOX::Abstract::MultiVector::DenseMatrix
&b, double beta, NOX::Abstract::MultiVector &result_x) const

Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x.

This method is implemented using getConstraintDerivativesX() and the
NOX::Abstract::MultiVector::update() method. ";


// File: classLOCA_1_1SingularJacobianSolve_1_1Default.xml
%feature("docstring") LOCA::SingularJacobianSolve::Default "

Default singular Jacobian solve computation class

This class computes the solution to $J x = b$ using the
applyJacobianInverse method of the underlying group ignoring the null
vector data.

C++ includes: LOCA_SingularJacobianSolve_Default.H ";

%feature("docstring")  LOCA::SingularJacobianSolve::Default::Default "LOCA::SingularJacobianSolve::Default::Default(Teuchos::ParameterList
&params)

Constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Default::Default "LOCA::SingularJacobianSolve::Default::Default(const Default &source)

Copy constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Default::~Default
"LOCA::SingularJacobianSolve::Default::~Default()

Destructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Default::clone "LOCA::SingularJacobianSolve::Generic *
LOCA::SingularJacobianSolve::Default::clone() const

Clone function. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Default::reset "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::reset(Teuchos::ParameterList
&params)

Reset parameters.

There are no additional parameters for the default calculation. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Default::compute "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::compute(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector &input, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector &result)

Computes the solution as described above. ";

%feature("docstring")
LOCA::SingularJacobianSolve::Default::computeMulti "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Default::computeMulti(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector *const *inputs, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector **results, int nVecs)

Computes solution for multiple RHS using applyJacobianInverseMulti. ";


// File: classLOCA_1_1Parameter_1_1DefaultFunctor.xml
%feature("docstring") LOCA::Parameter::DefaultFunctor "

Default function object for setting a single parameter in a single
object using a data member pointer.

The constructor takes a reference to an object object of type\\\\
ObjectType and a pointer object_val_ptr to a data member of class of
ObjectType of type ValueType. The parameter is set to value via

C++ includes: LOCA_Parameter_Entry.H ";

%feature("docstring")  LOCA::Parameter::DefaultFunctor::DefaultFunctor
"LOCA::Parameter::DefaultFunctor< ObjectType, ValueType
>::DefaultFunctor(ObjectType &object, ValueType
ObjectType::*object_val_ptr)

Constructor.

object is a reference to the object to set the parameter in, and
object_val_ptr is a pointer to a data member of type ValueType of that
class. ";

%feature("docstring")
LOCA::Parameter::DefaultFunctor::~DefaultFunctor "virtual
LOCA::Parameter::DefaultFunctor< ObjectType, ValueType
>::~DefaultFunctor()

Destructor. ";

%feature("docstring")  LOCA::Parameter::DefaultFunctor::set "virtual
void LOCA::Parameter::DefaultFunctor< ObjectType, ValueType
>::set(const ValueType &value)

Set parameter using object and data member pointer. ";

%feature("docstring")  LOCA::Parameter::DefaultFunctor::get "virtual
ValueType LOCA::Parameter::DefaultFunctor< ObjectType, ValueType
>::get() const

Get parameter value this object represents. ";


// File: classLOCA_1_1Eigensolver_1_1DefaultStrategy.xml
%feature("docstring") LOCA::Eigensolver::DefaultStrategy "

Default eigensolver strategy.

This class implements a default eigensolver strategy that does not
compute any eigenvalues.

C++ includes: LOCA_Eigensolver_DefaultStrategy.H ";

%feature("docstring")
LOCA::Eigensolver::DefaultStrategy::DefaultStrategy "LOCA::Eigensolver::DefaultStrategy::DefaultStrategy(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor. ";

%feature("docstring")
LOCA::Eigensolver::DefaultStrategy::~DefaultStrategy "LOCA::Eigensolver::DefaultStrategy::~DefaultStrategy()

Destructor. ";

%feature("docstring")
LOCA::Eigensolver::DefaultStrategy::computeEigenvalues "NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::DefaultStrategy::computeEigenvalues(NOX::Abstract::Group
&group, Teuchos::RCP< std::vector< double > > &evals_r, Teuchos::RCP<
std::vector< double > > &evals_i, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_r, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_i)

Compute eigenvalues/eigenvectors.

The implementation here does nothing and always returns
NOX::Abstract::Group::Ok. Note that this implies the returned ref-
count pointers are null. ";


// File: classLOCA_1_1SaveEigenData_1_1DefaultStrategy.xml
%feature("docstring") LOCA::SaveEigenData::DefaultStrategy "

Default strategy for saving eigenvector/value data.

This class implements a default strategy for saving eigenvectors and
eigenvalues that does nothing and exists so the LOCA::Stepper always
has an object to pass eigen data to.

C++ includes: LOCA_SaveEigenData_DefaultStrategy.H ";

%feature("docstring")
LOCA::SaveEigenData::DefaultStrategy::DefaultStrategy "LOCA::SaveEigenData::DefaultStrategy::DefaultStrategy(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor. ";

%feature("docstring")
LOCA::SaveEigenData::DefaultStrategy::~DefaultStrategy "LOCA::SaveEigenData::DefaultStrategy::~DefaultStrategy()

Destructor. ";

%feature("docstring")  LOCA::SaveEigenData::DefaultStrategy::save "NOX::Abstract::Group::ReturnType
LOCA::SaveEigenData::DefaultStrategy::save(Teuchos::RCP< std::vector<
double > > &evals_r, Teuchos::RCP< std::vector< double > > &evals_i,
Teuchos::RCP< NOX::Abstract::MultiVector > &evecs_r, Teuchos::RCP<
NOX::Abstract::MultiVector > &evecs_i)

Save eigenvalues/eigenvectors.

The implementation here does nothing and always returns
NOX::Abstract::Group::Ok. ";


// File: classLOCA_1_1Homotopy_1_1DeflatedGroup.xml
%feature("docstring") LOCA::Homotopy::DeflatedGroup "

LOCA's Homotopy Algorithm.

The HomotopyGroup is a concrete implementation of the
LOCA::Continuation::AbstractGroup that modifies the set of nonlinear
equations to be solved to allow for Homotopy to be applied to the
system. This object should be used in conjunction with the
LOCA::Stepper object to drive the continuation. This algorithm solves
a system of nonlinear equations supplied by the user ( $ F(x) $)
through continuation. An artificial parameter $ \\\\lambda $ is used
to control the continuation. The idea is to solve a simple equation
starting at $ \\\\lambda $ = 0 and, using the solution from the
previous step, solve systems of equations that gets progressively
closer to the true system of interest ( at $ \\\\lambda $ = 1.0 we
recover the original equations $ F(x) $). By constraining the
definition of $ g(x, \\\\lambda) $ and using artificial parameter
contiuation, the continuation branch should be free of multiplicity
and bifurcation phenomena.

The modified system of equations, $ g(x, \\\\lambda) $, supplied by
the HomotopyGroup is defined as:

\\\\[ g(x, \\\\lambda) = \\\\lambda F(x) + (1.0 - \\\\lambda)(x -
a)(S) \\\\]

where $x$ is the solution vector, $ \\\\lambda $ is an artificial
parameter, $ F(x) $ is the set of nonlinear equations the user
supplies, $ g(x) $ is the corresponding set of homotopy equations that
LOCA will solve, $ a $ is a random vector, and $ S $ is a scaling
factor used to switch sign of the last term (typically valued 1.0 or
-1.0).

This group requires the loca Stepper for continuation from $
\\\\lambda $ = 0.0 (a simple set of equations to solve) to $
\\\\lambda $ = 1.0 (the set of equations requested by the user, $ F(x)
$). The Homotopy::Group will generate the Stepper parameter sublist in
the parameter list that is passed in to the constructor. The user is
free to modify this list (it sets default values) before passing it
into the stepper object but should NOT change the starting and
stopping values for the continuation parameter.

References:

ALGORITHM 652 HOMPACK: A Suite of Codes for Globally Convergent
Homotopy Algorithms, Watson, L.T., Billups, S.C, and Morgan, A.P., ACM
Transactions on Mathematical Software, Vol. 13, No. 3, September 1987,
pp281-310.

C++ includes: LOCA_Homotopy_DeflatedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::Homotopy::DeflatedGroup::clone(NOX::CopyType type=NOX::DeepCopy)
const

Clone function. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::setX "void
LOCA::Homotopy::DeflatedGroup::setX(const NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::computeX "void
LOCA::Homotopy::DeflatedGroup::computeX(const NOX::Abstract::Group &g,
const NOX::Abstract::Vector &d, double step)

Compute and return solution vector, x, where this.x = grp.x + step *
d. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::computeF()

Compute extended continuation equations. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::computeJacobian
"NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::computeJacobian()

Compute extended continuation jacobian. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::computeGradient
"NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::computeGradient()

Gradient is not defined for this system. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction for extended continuation system. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::applyJacobianMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::applyJacobianTransposeMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::applyJacobianInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isF "bool
LOCA::Homotopy::DeflatedGroup::isF() const

Return true if extended residual is valid. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isJacobian "bool LOCA::Homotopy::DeflatedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isGradient "bool LOCA::Homotopy::DeflatedGroup::isGradient() const

Always returns false. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isNewton "bool
LOCA::Homotopy::DeflatedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getX "const
NOX::Abstract::Vector & LOCA::Homotopy::DeflatedGroup::getX() const

Return extended solution vector. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getF "const
NOX::Abstract::Vector & LOCA::Homotopy::DeflatedGroup::getF() const

Return extended residual. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getNormF "double LOCA::Homotopy::DeflatedGroup::getNormF() const

Return 2-norm of extended residual. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getGradient "const NOX::Abstract::Vector &
LOCA::Homotopy::DeflatedGroup::getGradient() const

Gradient is never valid. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getNewton "const NOX::Abstract::Vector &
LOCA::Homotopy::DeflatedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getXPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::getXPtr() const

Return RCP to extended solution vector. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getFPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::getFPtr() const

Return RCP to extended residual. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::getGradientPtr() const

Gradient is never valid. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::DeflatedGroup::getNewtonPtr() const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::getNormNewtonSolveResidual "double
LOCA::Homotopy::DeflatedGroup::getNormNewtonSolveResidual() const

Returns 2-norm of extended Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::getUnderlyingGroup "Teuchos::RCP<
const LOCA::MultiContinuation::AbstractGroup >
LOCA::Homotopy::DeflatedGroup::getUnderlyingGroup() const

Return underlying group. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::getUnderlyingGroup "Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup >
LOCA::Homotopy::DeflatedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::copy "void
LOCA::Homotopy::DeflatedGroup::copy(const NOX::Abstract::Group
&source)

Assignment operator. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::setParamsMulti "void LOCA::Homotopy::DeflatedGroup::setParamsMulti(const std::vector<
int > &paramIDs, const NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::setParams "void
LOCA::Homotopy::DeflatedGroup::setParams(const ParameterVector &p)

Set the parameter vector in the group to p (pVector = p). ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::setParam "void
LOCA::Homotopy::DeflatedGroup::setParam(int paramID, double val)

Set parameter indexed by (integer) paramID. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::setParam "void
LOCA::Homotopy::DeflatedGroup::setParam(std::string paramID, double
val)

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getParams "const LOCA::ParameterVector &
LOCA::Homotopy::DeflatedGroup::getParams() const

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getParam "double LOCA::Homotopy::DeflatedGroup::getParam(int paramID) const

Return copy of parameter indexed by (integer) paramID. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getParam "double LOCA::Homotopy::DeflatedGroup::getParam(std::string paramID)
const

Return copy of parameter indexed by (std::string) paramID. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::computeDfDpMulti
"NOX::Abstract::Group::ReturnType
LOCA::Homotopy::DeflatedGroup::computeDfDpMulti(const std::vector< int
> &paramIDs, NOX::Abstract::MultiVector &dfdp, bool isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::preProcessContinuationStep "void
LOCA::Homotopy::DeflatedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::postProcessContinuationStep "void
LOCA::Homotopy::DeflatedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::projectToDraw "void LOCA::Homotopy::DeflatedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::projectToDrawDimension "int
LOCA::Homotopy::DeflatedGroup::projectToDrawDimension() const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::computeScaledDotProduct "double
LOCA::Homotopy::DeflatedGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute a scaled dot product. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::printSolution "void LOCA::Homotopy::DeflatedGroup::printSolution(const double
conParam) const

Function to print out solution and parameter after successful step. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::printSolution "void LOCA::Homotopy::DeflatedGroup::printSolution(const
NOX::Abstract::Vector &x, const double conParam) const

Function to print out a vector and parameter after successful step. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::scaleVector "void LOCA::Homotopy::DeflatedGroup::scaleVector(NOX::Abstract::Vector
&x) const

Scales a vector using scaling vector. ";

/*  Implementation of  */

/*  LOCA::BorderedSystem::AbstractGroup virtual methods

*/

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getBorderedWidth
"int LOCA::Homotopy::DeflatedGroup::getBorderedWidth() const

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::getUnborderedGroup "Teuchos::RCP<
const NOX::Abstract::Group >
LOCA::Homotopy::DeflatedGroup::getUnborderedGroup() const

Get bottom-level unbordered group. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isCombinedAZero
"bool LOCA::Homotopy::DeflatedGroup::isCombinedAZero() const

Indicates whether combined A block is zero. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isCombinedBZero
"bool LOCA::Homotopy::DeflatedGroup::isCombinedBZero() const

Indicates whether combined B block is zero. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::isCombinedCZero
"bool LOCA::Homotopy::DeflatedGroup::isCombinedCZero() const

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::extractSolutionComponent "void
LOCA::Homotopy::DeflatedGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::extractParameterComponent "void
LOCA::Homotopy::DeflatedGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::Homotopy::DeflatedGroup::loadNestedComponents "void
LOCA::Homotopy::DeflatedGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::fillA "void
LOCA::Homotopy::DeflatedGroup::fillA(NOX::Abstract::MultiVector &A)
const

Fill the combined A block as described above. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::fillB "void
LOCA::Homotopy::DeflatedGroup::fillB(NOX::Abstract::MultiVector &B)
const

Fill the combined B block as described above. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::fillC "void
LOCA::Homotopy::DeflatedGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const

Fill the combined C block as described above. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::DeflatedGroup "LOCA::Homotopy::DeflatedGroup::DeflatedGroup(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &hParams, const Teuchos::RCP<
LOCA::Homotopy::AbstractGroup > &grp, const Teuchos::RCP< const
NOX::Abstract::Vector > &start_vec, const std::vector< Teuchos::RCP<
const NOX::Abstract::Vector > > &prev_solns, const double
identity_sign=1.0)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list.

hParams:  [in] Homotopy parameters

grp:  [in] Group representing $f$. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::DeflatedGroup "LOCA::Homotopy::DeflatedGroup::DeflatedGroup(const DeflatedGroup
&source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::~DeflatedGroup "LOCA::Homotopy::DeflatedGroup::~DeflatedGroup()

Destructor. ";

%feature("docstring")  LOCA::Homotopy::DeflatedGroup::getHomotopyParam
"double LOCA::Homotopy::DeflatedGroup::getHomotopyParam() const

Get homotopy parameter. ";


// File: classLOCA_1_1DerivUtils.xml
%feature("docstring") LOCA::DerivUtils "

LOCA's generic derivative computation class to compute various
derivatives via finite differencing.

The DerivUtils class provides generic methods to compute the following
derivatives: \\\\[ \\\\frac{\\\\partial F}{\\\\partial
p},\\\\quad\\\\frac{\\\\partial Jn}{\\\\partial
p},\\\\quad\\\\frac{\\\\partial Jn}{\\\\partial x}a \\\\] where $J =
\\\\partial F/\\\\partial x$ and $n$, $a$ are vectors. These
quantities are calculate by finite differencing.

C++ includes: LOCA_DerivUtils.H ";

%feature("docstring")  LOCA::DerivUtils::DerivUtils "LOCA::DerivUtils::DerivUtils(const Teuchos::RCP< LOCA::GlobalData >
&global_data, double perturb=1.0e-6)

Default constructor. perturb is the relative perturbation size used in
differencing calculations. ";

%feature("docstring")  LOCA::DerivUtils::DerivUtils "LOCA::DerivUtils::DerivUtils(const DerivUtils &)

Copy constructor. ";

%feature("docstring")  LOCA::DerivUtils::~DerivUtils "LOCA::DerivUtils::~DerivUtils()

Destructor. ";

%feature("docstring")  LOCA::DerivUtils::clone "Teuchos::RCP<
LOCA::DerivUtils > LOCA::DerivUtils::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. Creates a copy of the DerivUtils object of the same
type. ";

%feature("docstring")  LOCA::DerivUtils::computeDfDp "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDfDp(LOCA::MultiContinuation::AbstractGroup
&grp, const std::vector< int > &param_ids, NOX::Abstract::MultiVector
&result, bool isValidF) const

Compute derivative of f with respect to parameter, identified by
param_id. ";

%feature("docstring")  LOCA::DerivUtils::computeDJnDp "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDJnDp(LOCA::MultiContinuation::AbstractGroup
&, const std::vector< int > &paramIDs, const NOX::Abstract::Vector
&nullVector, NOX::Abstract::MultiVector &result, bool isValid) const

Compute derivative of Jn with respect to particular parameter
param_id. ";

%feature("docstring")  LOCA::DerivUtils::computeDJnDxa "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDJnDxa(LOCA::MultiContinuation::AbstractGroup
&grp, const NOX::Abstract::Vector &nullVector, const
NOX::Abstract::MultiVector &aVector, NOX::Abstract::MultiVector
&result) const

Compute vector (Jn)_{x}a given multi-vector a. ";

%feature("docstring")  LOCA::DerivUtils::computeDJnDxa "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDJnDxa(LOCA::MultiContinuation::AbstractGroup
&grp, const NOX::Abstract::Vector &nullVector, const
NOX::Abstract::MultiVector &aVector, const NOX::Abstract::Vector
&JnVector, NOX::Abstract::MultiVector &result) const

Compute vector (Jn)_{x}a given multi-vector a, given JnVector. ";

%feature("docstring")  LOCA::DerivUtils::computeDwtJnDp "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJnDp(LOCA::MultiContinuation::AbstractGroup
&grp, const std::vector< int > &paramIDs, const NOX::Abstract::Vector
&w, const NOX::Abstract::Vector &nullVector,
NOX::Abstract::MultiVector::DenseMatrix &result, bool isValid) const

Compute derivative of w^TJn with respect to particular parameter
param_id. ";

%feature("docstring")  LOCA::DerivUtils::computeDwtJDp "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJDp(LOCA::MultiContinuation::AbstractGroup
&grp, const std::vector< int > &paramIDs, const NOX::Abstract::Vector
&w, NOX::Abstract::MultiVector &result, bool isValid) const

Compute derivative of w^TJ with respect to particular parameter
param_id. ";

%feature("docstring")  LOCA::DerivUtils::computeDwtJnDx "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJnDx(LOCA::MultiContinuation::AbstractGroup
&grp, const NOX::Abstract::Vector &w, const NOX::Abstract::Vector
&nullVector, NOX::Abstract::Vector &result) const

Compute vector (w^TJn)_{x}. ";

%feature("docstring")  LOCA::DerivUtils::computeDwtJnDx "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJnDx(LOCA::MultiContinuation::AbstractGroup
&grp, const NOX::Abstract::MultiVector &w, const NOX::Abstract::Vector
&nullVector, NOX::Abstract::MultiVector &result) const

Compute vector (w^TJn)_{x}. ";

%feature("docstring")  LOCA::DerivUtils::computeDCeDp "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDCeDp(LOCA::Hopf::MooreSpence::AbstractGroup
&grp, const std::vector< int > &paramIDs, const NOX::Abstract::Vector
&yVector, const NOX::Abstract::Vector &zVector, double w,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag, bool isValid) const

Compute derivative of (J+iwM)(y+iz) with respect to parameter,. ";

%feature("docstring")  LOCA::DerivUtils::computeDCeDxa "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDCeDxa(LOCA::Hopf::MooreSpence::AbstractGroup
&grp, const NOX::Abstract::Vector &yVector, const
NOX::Abstract::Vector &zVector, double w, const
NOX::Abstract::MultiVector &aVector, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Compute vector (J+iwM)(y+iz))_{x}a given a. ";

%feature("docstring")  LOCA::DerivUtils::computeDCeDxa "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDCeDxa(LOCA::Hopf::MooreSpence::AbstractGroup
&grp, const NOX::Abstract::Vector &yVector, const
NOX::Abstract::Vector &zVector, double w, const
NOX::Abstract::MultiVector &aVector, const NOX::Abstract::Vector
&Ce_real, const NOX::Abstract::Vector &Ce_imag,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag) const

Compute vector (J+iwM)(y+iz))_{x}a given a and (J+iwM)(y+iz) vector.
";

%feature("docstring")  LOCA::DerivUtils::computeDwtCeDp "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtCeDp(LOCA::Hopf::MinimallyAugmented::AbstractGroup
&grp, const std::vector< int > &paramIDs, const NOX::Abstract::Vector
&w1, const NOX::Abstract::Vector &w2, const NOX::Abstract::Vector
&yVector, const NOX::Abstract::Vector &zVector, double omega,
NOX::Abstract::MultiVector::DenseMatrix &result_real,
NOX::Abstract::MultiVector::DenseMatrix &result_imag, bool isValid)
const

Compute derivative of (w1+iw2)^T(J+iwM)(y+iz) w.r.t. parameter p. ";

%feature("docstring")  LOCA::DerivUtils::computeDwtCeDx "NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtCeDx(LOCA::Hopf::MinimallyAugmented::AbstractGroup
&grp, const NOX::Abstract::Vector &w1, const NOX::Abstract::Vector
&w2, const NOX::Abstract::Vector &yVector, const NOX::Abstract::Vector
&zVector, double omega, NOX::Abstract::Vector &result_real,
NOX::Abstract::Vector &result_imag) const

Compute vector (w1+iw2)^T(J+iwM)(y+iz))_{x}. ";


// File: classLOCA_1_1Parameter_1_1Entry.xml
%feature("docstring") LOCA::Parameter::Entry "

Parameter entry interface class templated on ValueType.

This class provides the interface that all parameter entry classes
should implement. It is templated on the ValueType, which is the type
that the underlying parameter is stored as.

C++ includes: LOCA_Parameter_Entry.H ";

%feature("docstring")  LOCA::Parameter::Entry::Entry "LOCA::Parameter::Entry< ValueType >::Entry()

Default constructor. ";

%feature("docstring")  LOCA::Parameter::Entry::~Entry "virtual
LOCA::Parameter::Entry< ValueType >::~Entry()

Destructor. ";

%feature("docstring")  LOCA::Parameter::Entry::setValue "virtual void
LOCA::Parameter::Entry< ValueType >::setValue(const ValueType
&value)=0

Set parameter this object represents to value. ";

%feature("docstring")  LOCA::Parameter::Entry::getValue "virtual
ValueType LOCA::Parameter::Entry< ValueType >::getValue() const =0

Get parameter value this object represents. ";

%feature("docstring")  LOCA::Parameter::Entry::setIsInLibrary "virtual void LOCA::Parameter::Entry< ValueType >::setIsInLibrary()=0

Informs entry that it is now stored in the library.

This is used primarily for informing the entry on how to delete itself
when deleting the library. ";


// File: classLOCA_1_1BorderedSolver_1_1EpetraAugmented.xml
%feature("docstring") LOCA::BorderedSolver::EpetraAugmented "

Bordered system solver strategy based on augmenting the Jacobian
operator.

This class solves the extended system of equations \\\\[
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
F \\\\\\\\ G \\\\end{bmatrix} \\\\] by forming an augmented
Epetra_Operator representing \\\\[ \\\\begin{bmatrix} J & A \\\\\\\\
B^T & C \\\\end{bmatrix} \\\\] by creating a new Epetra_Map for the
additional equations.

C++ includes: LOCA_BorderedSolver_EpetraAugmented.H ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::EpetraAugmented "LOCA::BorderedSolver::EpetraAugmented::EpetraAugmented(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Currently none are
referenced. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::~EpetraAugmented "LOCA::BorderedSolver::EpetraAugmented::~EpetraAugmented()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::setMatrixBlocks "void
LOCA::BorderedSolver::EpetraAugmented::setMatrixBlocks(const
Teuchos::RCP< const LOCA::BorderedSolver::AbstractOperator > &op,
const Teuchos::RCP< const NOX::Abstract::MultiVector > &blockA, const
Teuchos::RCP< const LOCA::MultiContinuation::ConstraintInterface >
&blockB, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > &blockC)

Set blocks.

The blockA or blockC pointer may be null if either is zero. Whether
block B is zero will be determined by querying blockB via
ConstraintInterface::isConstraintDerivativesXZero. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::initForSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::initForSolve()

Intialize solver for a solve.

This should be called after setMatrixBlocks(), but before
applyInverse(). ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::initForTransposeSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::initForTransposeSolve()

Intialize solver for a transpose solve.

This should be called after setMatrixBlocks(), but before
applyInverseTranspose(). ";

%feature("docstring")  LOCA::BorderedSolver::EpetraAugmented::apply "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::apply(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J*X + A*Y \\\\\\\\ B^T*X + C*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::applyTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::applyTranspose(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix transpose-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J^T & B \\\\\\\\ A^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J^T*X + B*Y \\\\\\\\ A^T*X + C^T*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::applyInverse "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system using the technique described above.

The params argument is the linear solver parameters. If isZeroF or
isZeroG is true, than the corresponding F or G pointers may be NULL.

Note that if either the A or B blocks are zero, the system is solved
using a simple block elimination scheme instead of the Householder
scheme. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraAugmented::applyInverseTranspose "virtual
NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraAugmented::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the transpose of the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1BorderedSolver_1_1EpetraHouseholder.xml
%feature("docstring") LOCA::BorderedSolver::EpetraHouseholder "

Bordered system solver strategy based on Householder transformations.

This class solves the extended system of equations \\\\[
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
F \\\\\\\\ G \\\\end{bmatrix} \\\\] using Householder tranformations.
The algorithm works as follows: First consider a slightly rearranged
version of the extended system of equations: \\\\[ \\\\begin{bmatrix}
C & B^T \\\\\\\\ A & J \\\\end{bmatrix} \\\\begin{bmatrix} Y \\\\\\\\
X \\\\end{bmatrix} = \\\\begin{bmatrix} G \\\\\\\\ F \\\\end{bmatrix}.
\\\\] Let \\\\[ Q^T \\\\begin{bmatrix} C^T \\\\\\\\ B \\\\end{bmatrix}
= \\\\begin{bmatrix} R \\\\\\\\ 0 \\\\end{bmatrix} \\\\] be the QR
decomposition of the constraints matrix where
$Q\\\\in\\\\Re^{n+m\\\\times n+m}$ and $R\\\\in\\\\Re^{m\\\\times m}$.
Define \\\\[ \\\\begin{bmatrix} Z_Y \\\\\\\\ Z_X \\\\end{bmatrix} =
Q^T \\\\begin{bmatrix} Y \\\\\\\\ X \\\\end{bmatrix}, \\\\] then the
extended system of equations is equivalent to \\\\[ \\\\begin{bmatrix}
R^T & 0 \\\\\\\\ [A & J] Q \\\\end{bmatrix} \\\\begin{bmatrix} Z_Y
\\\\\\\\ Z_X \\\\end{bmatrix} = \\\\begin{bmatrix} G \\\\\\\\ F
\\\\end{bmatrix} \\\\] and hence \\\\[ \\\\begin{split} Z_Y &= R^{-T}
G \\\\\\\\ [A \\\\;\\\\; J] Q \\\\begin{bmatrix} 0 \\\\\\\\ Z_X
\\\\end{bmatrix} &= F - [A \\\\;\\\\; J] Q \\\\begin{bmatrix} Z_Y
\\\\\\\\ 0 \\\\end{bmatrix}. \\\\end{split} \\\\] This last equation
equation can be written \\\\[ P Z_X = \\\\tilde{F} \\\\] where
$P\\\\in\\\\Re^{n\\\\times n}$ is given by \\\\[ P Z_X = [A \\\\;\\\\;
J] Q \\\\begin{bmatrix} 0 \\\\\\\\ Z_X \\\\end{bmatrix} \\\\] and
\\\\[ \\\\tilde{F} = F - [A \\\\;\\\\; J] Q \\\\begin{bmatrix} Z_Y
\\\\\\\\ 0 \\\\end{bmatrix}. \\\\] We then recover $X$ and $Y$ by
\\\\[ \\\\begin{bmatrix} Y \\\\\\\\ X \\\\end{bmatrix} = Q
\\\\begin{bmatrix} Z_Y \\\\\\\\ Z_X \\\\end{bmatrix}. \\\\] It can be
further shown that the $P$ operator above can be written \\\\[ P = J +
U V^T \\\\] where $U = A*Y_1 + J*Y_2$, $V = Y_2*T^T$ and $Y = [Y_1 ;
Y_2]$. The equation $P Z_X = \\\\tilde{F}$ is solved using an
iterative solver using the definition of $P Z_X$ above, in this case
AztecOO. The system is preconditioned using the preconditioner for
$J$. The operator $Q$ is generated using the standard Householder QR
algorithm (Algorithm 5.2.1, G. Golub and C. Van Loan, \"Matrix
Computations,\" 3rd Edition, Johns Hopkins, Baltimore, 1996) and is
stored using the compact WY representation: $Q = I + Y T Y^T$ (see R.
Schreiber and C. Van Loan, \"A Storage-Efficient WY Representation
for Products of Householder Transformations,\" SIAM J. Sci. Stat.
Comput., Vol. 10, No. 1, pp. 53-57, January 1989).

The operator representing $P$ is encapsulated in the class
LOCA::Epetra::LowRankUpdateRowMatrix if $J$ is an Epetra_RowMatrix and
LOCA::Epetra::LowRankUpdateOp otherwise. If the row matrix version is
available $P$ can be scaled and also used to construct a
preconditioner. If \"Include UV In Preconditioner\" is true as
discussed below, the $U$ and $V$ terms will be included when computing
this preconditioner, which can help stability when $J$ is nearly
singular.

The class is intialized via the solverParams parameter list argument
to the constructor. The parameters this class recognizes are:
\"Preconditioner Method\" -- [string] (default: \"Jacobian\") - Method
for preconditioning the $P$ operator. Choices are: \"Jacobian\"
(default) -- Use the preconditioner for $J$

\"SMW\" -- Use the Sherman-Morrison-Woodbury formula for the inverse
of $P$, replacing the inverse of $J$ with the preconditioner for $J$.

\"Scale Augmented Rows\" -- [bool] (default: true) - Scale augmented
rows to unit 2-norm before computing QR factorization.

\"Include UV In Preconditioner\" -- [bool] (default: false) - Flag
indicating whether to use the $U$ and $V$ terms in the preconditioner
for $P$ when using the \"Jacobian\" preconditioner method.

\"Use P For Preconditioner\" -- [bool] (default: false) - Flag
indicating whether to use the representation of $P$ as a
LOCA::Epetra::LowRankUpdateRowMatrix for computing the preconditioner
when using the \"Jacobian\" preconditioner method. This is valid only
for preconditioners that accept an Epetra_RowMatrix interface.

\"Transpose Solver Method\" -- [string] (default: \"Transpose
Preconditioner\") Method for preconditioning the transpose linear
system. See LOCA::Epetra::TransposeLinearSystem::Factory for available
choices.

C++ includes: LOCA_BorderedSolver_EpetraHouseholder.H ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::EpetraHouseholder "LOCA::BorderedSolver::EpetraHouseholder::EpetraHouseholder(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters as described above ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::~EpetraHouseholder "LOCA::BorderedSolver::EpetraHouseholder::~EpetraHouseholder()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::setMatrixBlocks "void
LOCA::BorderedSolver::EpetraHouseholder::setMatrixBlocks(const
Teuchos::RCP< const LOCA::BorderedSolver::AbstractOperator > &op,
const Teuchos::RCP< const NOX::Abstract::MultiVector > &blockA, const
Teuchos::RCP< const LOCA::MultiContinuation::ConstraintInterface >
&blockB, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > &blockC)

Set blocks.

The blockA or blockC pointer may be null if either is zero. Whether
block B is zero will be determined by querying blockB via
ConstraintInterface::isConstraintDerivativesXZero. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::initForSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::initForSolve()

Intialize solver for a solve.

This should be called after setMatrixBlocks(), but before
applyInverse(). ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::initForTransposeSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::initForTransposeSolve()

Intialize solver for a transpose solve.

This should be called after setMatrixBlocks(), but before
applyInverseTranspose(). ";

%feature("docstring")  LOCA::BorderedSolver::EpetraHouseholder::apply
"NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::apply(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J*X + A*Y \\\\\\\\ B^T*X + C*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::applyTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::applyTranspose(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix transpose-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J^T & B \\\\\\\\ A^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J^T*X + B*Y \\\\\\\\ A^T*X + C^T*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::applyInverse "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system using the technique described above.

The params argument is the linear solver parameters. If isZeroF or
isZeroG is true, than the corresponding F or G pointers may be NULL.

Note that if either the A or B blocks are zero, the system is solved
using a simple block elimination scheme instead of the Householder
scheme. ";

%feature("docstring")
LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::EpetraHouseholder::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the transpose of the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1ErrorCheck.xml
%feature("docstring") LOCA::ErrorCheck "

An Error checking algorithm for NOX/LOCA routines.

This object will check the return types on objects and print a warning
or throw an error if appropriate

C++ includes: LOCA_ErrorCheck.H ";

%feature("docstring")  LOCA::ErrorCheck::ErrorCheck "LOCA::ErrorCheck::ErrorCheck(const Teuchos::RCP< LOCA::GlobalData >
&global_data)

Constructor. ";

%feature("docstring")  LOCA::ErrorCheck::~ErrorCheck "LOCA::ErrorCheck::~ErrorCheck()

Destructor. ";

%feature("docstring")  LOCA::ErrorCheck::throwError "void
LOCA::ErrorCheck::throwError(const std::string
&callingFunction=\"<Unknown Method>\", const std::string
&message=\"\", const std::string &throwLabel=\"LOCA Error\")

Generic call to throw that prints info to the screen. ";

%feature("docstring")  LOCA::ErrorCheck::printWarning "void
LOCA::ErrorCheck::printWarning(const std::string
&callingFunction=\"<Unknown Method>\", const std::string
&message=\"\")

Generic printing algorithm for sending warnings to the screen. ";

%feature("docstring")  LOCA::ErrorCheck::checkReturnType "void
LOCA::ErrorCheck::checkReturnType(const
NOX::Abstract::Group::ReturnType &status, const std::string
&callingFunction=std::string(\"<Unknown Method>\"))

Checks the supplied return type and performs an appropriate action.

This routine performs the following actions depending on the value of
status NOX::Abstract::Group::Ok -- nothing

NOX::Abstract::Group::Failed -- print message and throw a std::string

NOX::Abstract::Group::NotDefined -- print message and throw a
std::string

NOX::Abstract::Group::BadDependency -- print message and throw a
std::string

NOX::Abstract::Group::NotConverged -- print a warning message ";

%feature("docstring")  LOCA::ErrorCheck::checkReturnType "void
LOCA::ErrorCheck::checkReturnType(const
NOX::Abstract::Group::ReturnType &status, const ActionType &action,
const std::string &callingFunction=std::string(\"<Unknown Method>\"),
const std::string &message=std::string(\"\"))

Checks the return type for the NOX::AbstractGroup and may throw an
error or print a warning to the screen based on the ActionType
requested. ";

%feature("docstring")  LOCA::ErrorCheck::combineReturnTypes "NOX::Abstract::Group::ReturnType
LOCA::ErrorCheck::combineReturnTypes(const
NOX::Abstract::Group::ReturnType &status1, const
NOX::Abstract::Group::ReturnType &status2)

Combines two return types.

If either return type is NOX::Abstract::Group::NotDefined, returns
NotDefined. Otherwise if either is BadDependcy, returns BadDependency,
if either is Failed, returns Failed, if either is NotConverged,
returns NotConverged, and otherwise returns Ok. ";

%feature("docstring")  LOCA::ErrorCheck::combineAndCheckReturnTypes "NOX::Abstract::Group::ReturnType
LOCA::ErrorCheck::combineAndCheckReturnTypes(const
NOX::Abstract::Group::ReturnType &status1, const
NOX::Abstract::Group::ReturnType &status2, const std::string
&callingFunction=std::string(\"<Unknown Method>\"))

Combines two return types and checks the first.

First combines status1 and status2 using combineReturnTypes() and
checks the first using checkReturnType(). ";


// File: classLOCA_1_1Epetra_1_1TransposeLinearSystem_1_1ExplicitTranspose.xml
%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose "

Method for solving the transpose of a linear system by explicitly
forming the transpose of the matrix.

C++ includes: LOCA_Epetra_TransposeLinearSystem_ExplicitTranspose.H ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::ExplicitTranspose
"LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::ExplicitTranspose(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linsys)

Constructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::~ExplicitTranspose
"LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::~ExplicitTranspose()

Destructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::applyJacobianTransposeInverse
"bool
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector
&result)

Applies the inverse of the Jacobian matrix transpose to the given
input vector and puts the answer in result.

Computes \\\\[ v = J^{-T} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The parameter list contains the linear solver options. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::createJacobianTranspose
"bool
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::createJacobianTranspose()

Evaluates the Jacobian-transpose based on the solution vector x.

Note: For flexibility, this method does not compute the original
Jacobian matrix. It uses whatever is currently stored in the linear
system. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::createTransposePreconditioner
"bool
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::createTransposePreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &p)

Explicitly constructs a preconditioner based on the solution vector x
and the parameter list p.

Note: x is only needed for user-supplied preconditioners. When using a
built- in preconditioner (e.g., Ifpack), x will note be used. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::getJacobianTransposeOperator
"Teuchos::RCP< Epetra_Operator >
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::getJacobianTransposeOperator()

Get Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::getTransposePreconditioner
"Teuchos::RCP< Epetra_Operator >
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::getTransposePreconditioner()

Get transpose-preconditioner. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::setJacobianTransposeOperator
"void
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::setJacobianTransposeOperator(const
Teuchos::RCP< Epetra_Operator > &new_jac_trans)

Set Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::setTransposePreconditioner
"void
LOCA::Epetra::TransposeLinearSystem::ExplicitTranspose::setTransposePreconditioner(const
Teuchos::RCP< Epetra_Operator > &new_prec_trans)

Set transpose-preconditioner. ";


// File: classLOCA_1_1TurningPoint_1_1MinimallyAugmented_1_1ExtendedGroup.xml
%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup "

A group representing the minimally augemented turning point equations.

The LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup is a
concrete implementation of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup and
LOCA::Extended::MultiAbstractGroup that defines the following extended
set of equations that are regular at a generic turning point: \\\\[
G(z) = \\\\left[ \\\\begin{array}{c} F(x,p) \\\\\\\\ \\\\sigma
\\\\end{array} \\\\right] = 0 \\\\] where $z = [x,
p]\\\\in\\\\Re^{n+1}$, $x$ is the solution vector, $p$ is the
bifurcation parameter, $J$ is the Jacobian of $F$, and
$\\\\sigma\\\\in\\\\Re$ is a measure of the singularity of $F$ and is
defined via \\\\[ \\\\begin{bmatrix} J & a \\\\\\\\ b^T & 0
\\\\end{bmatrix} \\\\begin{bmatrix} v \\\\\\\\ \\\\sigma_1
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\begin{bmatrix} J^T & b \\\\\\\\ a^T & 0
\\\\end{bmatrix} \\\\begin{bmatrix} w \\\\\\\\ \\\\sigma_2
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\sigma = w^T J v/n \\\\] for any vectors $a$ and $b$ in
$\\\\Re^n$. Using these relationships, it is easy to show \\\\[
\\\\begin{split} \\\\sigma_x &= (w^T J v)_x/n = w^T J_x v/n \\\\\\\\
\\\\sigma_p &= (w^T J v)_p/n = w^T J_p v/n \\\\end{split} \\\\]

The group stores an underlying group of type
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup to represent the
equations $F(x,p) = 0$ and to manipulate the underlying Jacobian $J$.
This interface defines methods for computing the derivatives $(w^T J
v)_x$ and $(w^T J v)_p$ as well.

This class implements all of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup, and
LOCA::Extended::MultiAbstractGroup methods for this extended set of
equations and therefore is a complete group which can be passed to
most NOX solvers to locate a single turning point or to the
LOCA::Stepper to compute a family of turning points in a second
parameter.

The class is intialized via the tpParams parameter list argument to
the constructor. The parameters this class recognizes are:
\"Bifurcation Parameter\" -- [string] (Must be supplied) - Name of the
bifurcation parameter $p$

\"Bordered Solver Method\" -- [string] (default \"Bordering\") Method
for solving bordered systems of equations. See
LOCA::BorderedSolver::Factory for a description.

\"Symmetric Jacobian\" -- [bool] (default: false) - Flag indicating
whether Jacobian matrix $J$ is symmetric, in which case we force $a =
b$ and therefore the second tranpose solve for $w$ is unnecessary

\"Constraint Method\" -- [string] (default: \"Default\") - Type of
constraint method to use. Valid choices are \"Default\" (
LOCA::TurningPoint::MinimallyAugmented::Constraint) Default method
described above.

\"Modified\"
(LOCA::TurningPoint::MinimallyAugmented::ModifedConstraint) A modified
method that computes updates to the null vectors every nonlinear
interation, instead of directly solving for them

\"Initial Null Vector Compuation\" -- [string] (default: \"User
Provided\") - Method to compute initial $a$ and $b$ vectors. Valid
choices are: \"User Provided\" - Initial vectors are provided in the
parameter list, in which case the following parameters are relevant:
\"Initial A Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be
supplied) - Vector storing initial value for $a$ vector

\"Initial B Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be
supplied for nonsymmetric Jacobians) - Vector storing initial value
for $b$ vector

\"Solve df/dp\" - Compute $a = J^{-T}df/dp$ and $b = J^{-1} df/dp$
where $p$ is the bifurcation parameter.

\"Constant\" - Entries of $a$ and $b$ are set to 1.0

\"Null Vector Scaling Method\" -- [string] (default: \"Order N\") -
Method to scale $a$ and $b$. This determines the norm of these vectors
and the scaling of $\\\\sigma$. Valid choices are: \"None\" -- Use
initial scaling

\"Order 1\" -- Scale to unit norm

\"Order N\" -- Use vector length scaling

\"Update Null Vectors Every Continuation Step\" -- [bool] (default:
true) - Flag indicating whether to update $a$ and $b$ vectors via $a =
w$ and $b = v$ every continuation step

\"Update Null Vectors Every Nonlinear Iteration\" -- [bool] (default:
false) - Flag indicating whether to update $a$ and $b$ vectors via $a
= w$ and $b = v$ every nonlinear iteration

\"Multiply Null Vectors by Mass Matrix\" -- [bool] (default: false) -
Flag indicating whether to multiply $a$ and $b$ vectors by the mass
matrix $M = \\\\partial f/\\\\partial\\\\dot{x}$ at the strart of a
turning point calculation, and each time $a$ and $b$ are updated. This
can improve the scaling of these vectors, and may orthogonalize them
against structural null spaces (i.e., pressure null space for
incompressible Navier-Stokes).

C++ includes: LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setX "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector, x, to y. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeX "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute this.x = grp.x + step * d. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeF()

Compute the turning point equation residual $G$.

This method fills the extended residual \\\\[ G(z) = \\\\left[
\\\\begin{array}{c} F(x,p) \\\\\\\\ \\\\sigma \\\\end{array}
\\\\right]. \\\\] ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeJacobian
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeJacobian()

Compute the blocks of the Jacobian derivative of $G$.

This method computes the $J$, $F_p$, $\\\\sigma_x$ and $\\\\sigma_p$
blocks of the extended Jacobian: \\\\[ D_z G(z) = \\\\begin{bmatrix} J
& F_p \\\\\\\\ \\\\sigma_x & \\\\sigma_p \\\\end{bmatrix} \\\\] ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeGradient
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeGradient()

Gradient computation is not defined for this group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction using applyJacobianInverse(). ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Computes the extended Jacobian vector product. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTranspose
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Computes the extended Jacobian transpose vector product. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianInverse
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies the inverse of the extended Jacobian matrix. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose for extended system. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isF "bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isF() const

Return true if the extended residual $G$ is valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isJacobian "bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isJacobian()
const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isGradient "bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isGradient()
const

Always returns false. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isNewton "bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isNewton()
const

Return true if the extended Newton direction is valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getX "const
NOX::Abstract::Vector &
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getX() const

Return extended solution vector $z$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getF "const
NOX::Abstract::Vector &
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getF() const

Return extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNormF "double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNormF()
const

Return 2-norm of $G(z)$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getGradient "const NOX::Abstract::Vector &
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getGradient()
const

Vector returned is not valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNewton "const NOX::Abstract::Vector &
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNewton()
const

Return extended Newton direction. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getXPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getXPtr() const

Return RCP to extended solution vector $z$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getFPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getFPtr() const

Return RCP to extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getGradientPtr
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getGradientPtr()
const

Vector returned is not valid. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNewtonPtr()
const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual
"double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual()
const

Return the norm of the Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup
"Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup()
const

Return underlying group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup
"Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::copy "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParamsMulti
"void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParams "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParam "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParam "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::setParam(int
paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getParams "const LOCA::ParameterVector &
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getParams()
const

Return a const reference to the paramter vector owned by the group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getParam "double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getParam(int
paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getParam "double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeDfDpMulti
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::preProcessContinuationStep
"void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::postProcessContinuationStep
"void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::projectToDraw "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::projectToDrawDimension
"int
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::projectToDrawDimension()
const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::printSolution "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::printSolution(const
double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::printSolution "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &x_, const double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

/*  Implementation of  */

/*  LOCA::BorderedSystem::AbstractGroup virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getBorderedWidth
"int
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getBorderedWidth()
const

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getUnborderedGroup
"Teuchos::RCP< const NOX::Abstract::Group >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getUnborderedGroup()
const

Get bottom-level unbordered group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isCombinedAZero
"bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isCombinedAZero()
const

Indicates whether combined A block is zero. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isCombinedBZero
"bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isCombinedBZero()
const

Indicates whether combined B block is zero. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isCombinedCZero
"bool
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::isCombinedCZero()
const

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::extractSolutionComponent
"void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::extractParameterComponent
"void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::loadNestedComponents
"void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::fillA "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::fillA(NOX::Abstract::MultiVector
&A) const

Fill the combined A block as described above. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::fillB "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::fillB(NOX::Abstract::MultiVector
&B) const

Fill the combined B block as described above. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::fillC "void
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const

Fill the combined C block as described above. ";

/*  Implementation of LOCA::Abstract::TransposeSolveGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeInverse
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Solve Jacobian-tranpose system. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Solve Jacobian-tranpose system with multiple right-hand sides. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::ExtendedGroup "LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &tpParams, const Teuchos::RCP<
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup > &g)

Constructor with initial data passed through parameter lists. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::ExtendedGroup "LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::~ExtendedGroup
"LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::~ExtendedGroup()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getBifParam "double
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getBifParam()
const

Get bifurcation parameter. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getLeftNullVec
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getLeftNullVec()
const

Returns left null vector. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getRightNullVec
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getRightNullVec()
const

Returns right null vector v. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getAVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getAVec() const

Returns \"A\" vector. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getBVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup::getBVec() const

Returns \"B\". ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1ExtendedGroup.xml
%feature("docstring") LOCA::TurningPoint::MooreSpence::ExtendedGroup "

A group representing the Moore-Spence turning point equations.

The LOCA::TurningPoint::MooreSpence::ExtendedGroup is a concrete
implementation of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup and
LOCA::Extended::MultiAbstractGroup that defines the following extended
set of equations that are regular at a generic turning point: \\\\[
G(z) = \\\\left[ \\\\begin{array}{c} F(x,p) \\\\\\\\ Jn \\\\\\\\
l^Tn-1 \\\\end{array} \\\\right] = 0 \\\\] where $z = [x, n,
p]\\\\in\\\\Re^{2n+1}$, $x$ is the solution vector, $n$ is the null
vector, $l$ is the length normalization vector and $J$ is the Jacobian
of F.

The group stores an underlying group of type
LOCA::TurningPoint::MooreSpence AbstractGroup to represent the
equations $F(x,p) = 0$ and to manipulate the underlying Jacobian $J$.
Note that the entire extended Jacobian $D_z G$ is not stored in
memory, rather a block-elimination algorithm (bordering algorithm) is
used to compute linear solves of the extended Jacobian (see
LOCA::TurningPoint::MooreSpence::SolverFactory for more details).

This class implements all of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup, and
LOCA::Extended::MultiAbstractGroup methods for this extended set of
equations and therefore is a complete group which can be passed to
most NOX solvers to locate a single turning point or to the
LOCA::Stepper to compute a family of turning points in a second
parameter.

However, Jacobian-tranpose operations and gradient calculations cannot
be implemented efficiently and therefore gradient-base nonlinear
solvers such as steepest descent and Trust region methods cannot be
used to solve the extended turning point equations.

The class is intialized via the tpParams parameter list argument to
the constructor. The parameters this class recognizes are:
\"Bifurcation Parameter\" -- [string] (Must be supplied) Name of the
bifurcation parameter $p$

\"Length Normalization Vector\" --
[Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Vector
storing length normalization vector $l$

\"Initial Null Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must
be supplied) - Vector storing initial guess for the null vector $n$

\"Perturb Initial Solution\" -- [bool] (default: false) Flag
indicating whether to perturb the initial solution

\"Relative Perturbation Size\" -- [double] (default: 1.0e-3) Relative
perturbation size if perturbing the initial solution

\"Null Vector Scaling\" -- [string] (default: \"Order N\") - Method to
scale $l$. This determines the norm of $l$ (and hence, $n$). Valid
choices are: \"None\" -- Use initial scaling

\"Order 1\" -- Scale to unit norm

\"Order N\" -- Use vector length scaling

\"Update Null Vectors Every Continuation Step\" -- [bool] (default:
false) - Flag indicating whether to update length normalization vector
$l$ via $l = n$ every continuation step

\"Multiply Null Vectors by Mass Matrix\" -- [bool] (default: false) -
Flag indicating whether to multiply length scaling vector $l$ by the
mass matrix $M = \\\\partial f/\\\\partial\\\\dot{x}$ at the start of
a turning point calculation, and each time it is updated. This can
improve the scaling of this vector, and may orthogonalize it against
structural null spaces (i.e., pressure null space for incompressible
Navier-Stokes).

C++ includes: LOCA_TurningPoint_MooreSpence_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::clone "Teuchos::RCP<
NOX::Abstract::Group >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setX "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector, x, to y. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeX "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute this.x = grp.x + step * d. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeF()

Compute the turning point equation residual $G$.

This method fills the extended residual \\\\[ G(z) = \\\\left[
\\\\begin{array}{c} F(x,p) \\\\\\\\ Jn \\\\\\\\ l^Tn-1 \\\\end{array}
\\\\right]. \\\\] The solution component residual $F(x,p)$ and the
null-vector residual $Jn$ are calculated via the computeF and
applyJacobian methods of the underlying group. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeJacobian()

Compute the blocks of the Jacobian derivative of $G$.

This method computes the $J$, $\\\\partial F/\\\\partial p$, and
$\\\\partial Jn/\\\\partial p$ blocks of the extended Jacobian: \\\\[
D_z G(z) = \\\\begin{bmatrix} J & 0 & \\\\frac{\\\\partial
F}{\\\\partial p} \\\\\\\\ \\\\frac{\\\\partial Jn}{\\\\partial x} & J
& \\\\frac{\\\\partial Jn}{\\\\partial p} \\\\\\\\ 0 & l^T & 0
\\\\end{bmatrix} \\\\] by calling the computeJacobian, computeDfDp,
and computeDJnDp methods of the underlying group. The second
derivative matrix $\\\\partial Jn/\\\\partial x$ is not calculated
since only its action on vectors is needed for linear solves using the
bordering algorithm. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeGradient()

Gradient computation is not defined for this group. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction using applyJacobianInverse(). ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Computes the extended Jacobian vector product.

This method computes the extended Jacobian vector product \\\\[
\\\\begin{bmatrix} J & 0 & \\\\frac{\\\\partial F}{\\\\partial p}
\\\\\\\\ \\\\frac{\\\\partial Jn}{\\\\partial x} & J &
\\\\frac{\\\\partial Jn}{\\\\partial p} \\\\\\\\ 0 & l^T & 0
\\\\end{bmatrix} \\\\begin{bmatrix} a \\\\\\\\ b \\\\\\\\ c
\\\\end{bmatrix} = \\\\begin{bmatrix} Ja + \\\\frac{\\\\partial
F}{\\\\partial p}c \\\\\\\\ \\\\frac{\\\\partial Jn}{\\\\partial x}a +
Jb + \\\\frac{\\\\partial Jn}{\\\\partial p}c \\\\\\\\ l^T b
\\\\end{bmatrix} \\\\] using the applyJacobian and computeDJnDxa
methods of the underlying group where $a$, $b$, and $c$ are the
solution, null-vector, and paramter components of the given vector
input. Vectors input and result must be of type
LOCA::TurningPoint::MooreSpence::ExtendedVector, otherwise an error is
thrown. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTranspose
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose product is not defined by this group. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies the inverse of the extended Jacobian matrix using the
bordering algorithm.

This method is a special case of applyJacobianInverseMultiVector() for
a single right-hand-side. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system.

Uses a LOCA::TurningPoint::MooreSpence::SolverStrategy instantiated by
the LOCA::TurningPoint::MooreSpence::SolverFactory to implement the
solve. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isF "bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isF() const

Return true if the extended residual $G$ is valid. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isJacobian "bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isGradient "bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isNewton "bool
LOCA::TurningPoint::MooreSpence::ExtendedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getX "const
NOX::Abstract::Vector &
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getX() const

Return extended solution vector $z$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getF "const
NOX::Abstract::Vector &
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getF() const

Return extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormF "double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormF() const

Return 2-norm of $G(z)$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradient() const

Vector returned is not valid. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getXPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getXPtr() const

Return extended solution vector $z$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getFPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getFPtr() const

Return extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getGradientPtr() const

Vector returned is not valid. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNewtonPtr() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual
"double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()
const

Return the norm of the Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup()
const

Return underlying group. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::copy "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParamsMulti "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParams "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam(int paramID,
double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParams "const
LOCA::ParameterVector &
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParams() const

Return a const reference to the paramter vector owned by the group. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam "double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam(int paramID)
const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam "double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::preProcessContinuationStep
"void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::postProcessContinuationStep
"void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDraw "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDrawDimension
"int
LOCA::TurningPoint::MooreSpence::ExtendedGroup::projectToDrawDimension()
const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution(const
double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &x_, const double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

/*  Implementation of LOCA::Abstract::TransposeSolveGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeInverse
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Solve Jacobian-tranpose system. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::ExtendedGroup::applyJacobianTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Solve Jacobian-tranpose system with multiple right-hand sides. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup "LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &tpParams, const Teuchos::RCP<
LOCA::TurningPoint::MooreSpence::AbstractGroup > &g)

Constructor with initial data passed through parameter lists. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup "LOCA::TurningPoint::MooreSpence::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::~ExtendedGroup "LOCA::TurningPoint::MooreSpence::ExtendedGroup::~ExtendedGroup()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getBifParam "double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getBifParam() const

Get bifurcation parameter. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm "double
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm(const
NOX::Abstract::Vector &n) const

Defines null vector normalization $l^Tn$ equation. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm "void
LOCA::TurningPoint::MooreSpence::ExtendedGroup::lTransNorm(const
NOX::Abstract::MultiVector &n, NOX::Abstract::MultiVector::DenseMatrix
&result) const

null vector normalization for multivectors

Note: result should have 1 row and n.numVectors() columns. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getLengthVector "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedGroup::getLengthVector()
const

Return length normalization vector. ";


// File: classLOCA_1_1MultiContinuation_1_1ExtendedGroup.xml
%feature("docstring") LOCA::MultiContinuation::ExtendedGroup "

Base class for all continuation groups.

Continuation is defined as computing some curve
$(x(s),p(s))\\\\in\\\\Re^{n+1}$ such that $F(x(s),p(s))=0$ for some
parameterization $s$. Given some point $(x_0,p_0)$ on the curve,
another nearby point on the curve is calculated by first computing a
predictor direction $v\\\\in\\\\Re^{n+1}$ and the approximate point
$(x^\\\\ast,p^\\\\ast) = (x_0,p_0) + v\\\\Delta s$ where $\\\\Delta s$
is the step size. Then the next point on the curve is computed by
solving the extended set of equations \\\\[ \\\\begin{array}{cc}
F(x,p) &= 0 \\\\\\\\ g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)
&= 0 \\\\end{array} \\\\] for $(x,p)$. The equation
$g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)=0$ is called the
continuation equation and different choices of $g$ yield different
continuation methods.

Mathematically, this computation amounts to repeatedly computing
solutions to a constrained nonlinear system. This class provides a
common implementation for all continuation groups in terms of the
LOCA::MultiContinuation::ConstrainedGroup using a supplied group to
represent $F$ and an implementation of
LOCA::MultiContinuation::ConstraintInterface to represent $g$.

Note that this class has no public constructor other than the copy
constructor since it is intended to only provide an implemenation of
much of the continuation work. Each derived class that implements a
specific continuation strategy should provide its own public
constructor.

C++ includes: LOCA_MultiContinuation_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::MultiContinuation::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::setX "void LOCA::MultiContinuation::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computeX "void
LOCA::MultiContinuation::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute and return solution vector, x, where this.x = grp.x + step *
d. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeF()

Compute extended continuation equations. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeJacobian()

Compute extended continuation jacobian. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeGradient()

Gradient is not defined for this system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction for extended continuation system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::applyJacobianMultiVector "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::isF "bool LOCA::MultiContinuation::ExtendedGroup::isF() const

Return true if extended residual is valid. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isJacobian "bool
LOCA::MultiContinuation::ExtendedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isGradient "bool
LOCA::MultiContinuation::ExtendedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isNewton "bool
LOCA::MultiContinuation::ExtendedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::getX "const NOX::Abstract::Vector &
LOCA::MultiContinuation::ExtendedGroup::getX() const

Return extended solution vector. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::getF "const NOX::Abstract::Vector &
LOCA::MultiContinuation::ExtendedGroup::getF() const

Return extended residual. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getNormF "double
LOCA::MultiContinuation::ExtendedGroup::getNormF() const

Return 2-norm of extended residual. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::MultiContinuation::ExtendedGroup::getGradient() const

Gradient is never valid. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::MultiContinuation::ExtendedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::getXPtr
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getXPtr() const

Return extended solution vector. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::getFPtr
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getFPtr() const

Return extended residual. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getGradientPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getGradientPtr() const

Gradient is never valid. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getNewtonPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedGroup::getNewtonPtr() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getNormNewtonSolveResidual "double
LOCA::MultiContinuation::ExtendedGroup::getNormNewtonSolveResidual()
const

Returns 2-norm of extended Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup() const

Return underlying group. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::MultiContinuation::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractStrategy  */

/* virtual methods

*/

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::copy "void LOCA::MultiContinuation::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getNumParams "int
LOCA::MultiContinuation::ExtendedGroup::getNumParams() const

Returns number of parameters. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::preProcessContinuationStep "void
LOCA::MultiContinuation::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::postProcessContinuationStep "void
LOCA::MultiContinuation::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computePredictor "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ExtendedGroup::computePredictor()

Compute predictor directions. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isPredictor "bool
LOCA::MultiContinuation::ExtendedGroup::isPredictor() const

Is Predictor valid. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::scaleTangent "void
LOCA::MultiContinuation::ExtendedGroup::scaleTangent()

Scales tangent to predictor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::setPredictorTangentDirection "void
LOCA::MultiContinuation::ExtendedGroup::setPredictorTangentDirection(const
LOCA::MultiContinuation::ExtendedVector &v, int i)

Sets tangent to predictor.

This is required by MF which takes the tangent space, orthogonalizes
it, and then sets it back in the group. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getPredictorTangent "const
LOCA::MultiContinuation::ExtendedMultiVector &
LOCA::MultiContinuation::ExtendedGroup::getPredictorTangent() const

Returns tangent to predictor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getScaledPredictorTangent "const LOCA::MultiContinuation::ExtendedMultiVector &
LOCA::MultiContinuation::ExtendedGroup::getScaledPredictorTangent()
const

Returns scaled tangent to predictor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::setPrevX "void
LOCA::MultiContinuation::ExtendedGroup::setPrevX(const
NOX::Abstract::Vector &y)

Set the previous solution vector y. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getPrevX "const
LOCA::MultiContinuation::ExtendedVector &
LOCA::MultiContinuation::ExtendedGroup::getPrevX() const

Gets the previous solution vector. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::setStepSize "void
LOCA::MultiContinuation::ExtendedGroup::setStepSize(double deltaS, int
i=0)

Set step size for continuation constraint equation i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getStepSize "double
LOCA::MultiContinuation::ExtendedGroup::getStepSize(int i=0) const

Get step size for continuation constraint equation i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::setContinuationParameter "void
LOCA::MultiContinuation::ExtendedGroup::setContinuationParameter(double
val, int i=0)

Sets the value for continuation parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameter "double
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameter(int
i=0) const

Returns the value for continuation parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterID "int
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterID(int
i=0) const

Get the continuation parameter id for parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterIDs "const std::vector< int > &
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterIDs()
const

Get the continuation parameter ids. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterName "std::string
LOCA::MultiContinuation::ExtendedGroup::getContinuationParameterName(int
i=0) const

Get the continuation parameter id for parameter i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getStepSizeScaleFactor "double
LOCA::MultiContinuation::ExtendedGroup::getStepSizeScaleFactor(int
i=0) const

Returns step size scale factor for constraint equation i. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::printSolution "void
LOCA::MultiContinuation::ExtendedGroup::printSolution() const

Prints the group. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::computeScaledDotProduct "double
LOCA::MultiContinuation::ExtendedGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &x, const NOX::Abstract::Vector &y) const

Computes a scaled dot product between two continuation vectors. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::projectToDrawDimension "int
LOCA::MultiContinuation::ExtendedGroup::projectToDrawDimension() const

Returns dimension of project to draw array. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::projectToDraw "void
LOCA::MultiContinuation::ExtendedGroup::projectToDraw(const
LOCA::MultiContinuation::ExtendedVector &x, double *px) const

Fills the project to draw array. ";

/*  Implementation of  */

/*  LOCA::BorderedSystem::AbstractGroup virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getBorderedWidth "int
LOCA::MultiContinuation::ExtendedGroup::getBorderedWidth() const

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::getUnborderedGroup "Teuchos::RCP< const NOX::Abstract::Group >
LOCA::MultiContinuation::ExtendedGroup::getUnborderedGroup() const

Get bottom-level unbordered group. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isCombinedAZero "bool
LOCA::MultiContinuation::ExtendedGroup::isCombinedAZero() const

Indicates whether combined A block is zero. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isCombinedBZero "bool
LOCA::MultiContinuation::ExtendedGroup::isCombinedBZero() const

Indicates whether combined B block is zero. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::isCombinedCZero "bool
LOCA::MultiContinuation::ExtendedGroup::isCombinedCZero() const

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::extractSolutionComponent "void
LOCA::MultiContinuation::ExtendedGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::extractParameterComponent "void
LOCA::MultiContinuation::ExtendedGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::loadNestedComponents "void
LOCA::MultiContinuation::ExtendedGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::fillA "void
LOCA::MultiContinuation::ExtendedGroup::fillA(NOX::Abstract::MultiVector
&A) const

Fill the combined A block as described above. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::fillB "void
LOCA::MultiContinuation::ExtendedGroup::fillB(NOX::Abstract::MultiVector
&B) const

Fill the combined B block as described above. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedGroup::fillC "void
LOCA::MultiContinuation::ExtendedGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const

Fill the combined C block as described above. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup "LOCA::MultiContinuation::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedGroup::~ExtendedGroup "LOCA::MultiContinuation::ExtendedGroup::~ExtendedGroup()

Destructor. ";


// File: classLOCA_1_1PhaseTransition_1_1ExtendedGroup.xml
%feature("docstring") LOCA::PhaseTransition::ExtendedGroup "";

/*  "Compute" functions.  */

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::setX "void LOCA::PhaseTransition::ExtendedGroup::setX(const
NOX::Abstract::Vector &y) ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::setX "void LOCA::PhaseTransition::ExtendedGroup::setX(const
LOCA::PhaseTransition::ExtendedVector &y)

See above. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::computeX
"void LOCA::PhaseTransition::ExtendedGroup::computeX(const
NOX::Abstract::Group &grp, const NOX::Abstract::Vector &d, double
step) ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::computeX
"void LOCA::PhaseTransition::ExtendedGroup::computeX(const
LOCA::PhaseTransition::ExtendedGroup &grp, const
LOCA::PhaseTransition::ExtendedVector &d, double step)

See above. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::computeF
"NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::computeF() ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::computeJacobian() ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params) ";

/*  Jacobian operations.  */

/* Operations using the Jacobian matrix. These may not be defined in
matrix-free scenarios.

*/

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::applyJacobian(const
LOCA::PhaseTransition::ExtendedVector &input,
LOCA::PhaseTransition::ExtendedVector &result) const ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

See above. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const LOCA::PhaseTransition::ExtendedVector &input,
LOCA::PhaseTransition::ExtendedVector &result) const ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const ";

/*  "Is" functions  */

/* Checks to see if various objects have been computed. Returns true
if the corresponding \"compute\" function has been called since the
last update to the solution vector (via instantiation or computeX).

*/

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::isF "bool LOCA::PhaseTransition::ExtendedGroup::isF() const ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::isJacobian "bool
LOCA::PhaseTransition::ExtendedGroup::isJacobian() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::isNewton
"bool LOCA::PhaseTransition::ExtendedGroup::isNewton() const ";

/*  "Get" functions  */

/* Note that these function do not check whether or not the vectors
are valid. Must use the \"Is\" functions for that purpose.

*/

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getX "const NOX::Abstract::Vector &
LOCA::PhaseTransition::ExtendedGroup::getX() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getF "const NOX::Abstract::Vector &
LOCA::PhaseTransition::ExtendedGroup::getF() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getNormF
"double LOCA::PhaseTransition::ExtendedGroup::getNormF() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getNewton
"const NOX::Abstract::Vector &
LOCA::PhaseTransition::ExtendedGroup::getNewton() const ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::PhaseTransition::ExtendedGroup::getGradient() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getXPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getXPtr() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getFPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getFPtr() const ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::getNewtonPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getNewtonPtr() const ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::getGradientPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedGroup::getGradientPtr() const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::copy "void LOCA::PhaseTransition::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Start methods for LOCA::Abstract::Group. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::setParams
"void LOCA::PhaseTransition::ExtendedGroup::setParams(const
LOCA::ParameterVector &p)

Set the parameter vector in the group to p (pVector = p). ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::setParam
"void LOCA::PhaseTransition::ExtendedGroup::setParam(int paramID,
double val)

Set parameter indexed by (integer) paramID. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::setParam
"void LOCA::PhaseTransition::ExtendedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getParams
"const LOCA::ParameterVector &
LOCA::PhaseTransition::ExtendedGroup::getParams() const

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getParam
"double LOCA::PhaseTransition::ExtendedGroup::getParam(int paramID)
const

Return copy of parameter indexed by (integer) paramID. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::getParam
"double LOCA::PhaseTransition::ExtendedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::printSolution "void
LOCA::PhaseTransition::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &solution, const double param) const

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::printSolution "void
LOCA::PhaseTransition::ExtendedGroup::printSolution(const double
param) const

Function to print out solution and parameter after successful step.

Empty default definition. ";

/*  Vectors  */

/*  IsValid flags  */

/* True if the current solution is up-to-date with respect to the
currect xVector.

*/

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::ExtendedGroup "LOCA::PhaseTransition::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > gD, const Teuchos::RCP<
Teuchos::ParameterList > &bifurcationparams_, const Teuchos::RCP<
LOCA::PhaseTransition::AbstractGroup > &grp_)

Constructor. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::ExtendedGroup "LOCA::PhaseTransition::ExtendedGroup::ExtendedGroup(const
LOCA::PhaseTransition::ExtendedGroup &source, NOX::CopyType
type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::~ExtendedGroup "LOCA::PhaseTransition::ExtendedGroup::~ExtendedGroup()

Destructor. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::PhaseTransition::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedGroup::print "void LOCA::PhaseTransition::ExtendedGroup::print() const

Print out the group. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::getUnderlyingGroup "virtual
Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::PhaseTransition::ExtendedGroup::getUnderlyingGroup() const

Return underlying group. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::getUnderlyingGroup "virtual
Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::PhaseTransition::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::setParamsMulti "void
LOCA::PhaseTransition::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::PhaseTransition::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValid_F)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";


// File: classLOCA_1_1Pitchfork_1_1MinimallyAugmented_1_1ExtendedGroup.xml
%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup "

A group representing the minimally augemented pitchfork equations.

The LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup is a concrete
implementation of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup and
LOCA::Extended::MultiAbstractGroup that defines the following extended
set of equations that are regular at a generic pitchfork: \\\\[ G(z) =
\\\\left[ \\\\begin{array}{c} F(x,p)+s\\\\psi\\\\\\\\ \\\\sigma
\\\\\\\\ \\\\langle x,\\\\psi \\\\rangle \\\\end{array} \\\\right] = 0
\\\\] where $z = [x, p, s]\\\\in\\\\Re^{n+2}$, $x$ is the solution
vector, $s$ is the slack variable representing the asymmetry, $p$ is
the bifurcation parameter, $\\\\psi$ is the asymmetric vector, and
$\\\\sigma\\\\in\\\\Re$ is a measure of the singularity of $F$ and is
defined via \\\\[ \\\\begin{bmatrix} J & a \\\\\\\\ b^T & 0
\\\\end{bmatrix} \\\\begin{bmatrix} v \\\\\\\\ \\\\sigma_1
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\begin{bmatrix} J^T & b \\\\\\\\ a^T & 0
\\\\end{bmatrix} \\\\begin{bmatrix} w \\\\\\\\ \\\\sigma_2
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\sigma = w^T J v/n \\\\] for any vectors $a$ and $b$ in
$\\\\Re^n$. Using these relationships, it is easy to show \\\\[
\\\\begin{split} \\\\sigma_x &= (w^T J v)_x/n = w^T J_x v/n \\\\\\\\
\\\\sigma_p &= (w^T J v)_p/n = w^T J_p v/n \\\\end{split} \\\\]

The group stores an underlying group of type
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup to represent the
equations $F(x,p) = 0$ and to manipulate the underlying Jacobian $J$.
This interface defines methods for computing the derivatives $(w^T J
v)_x$ and $(w^T J v)_p$ and computing the inner product $\\\\langle
\\\\psi,x \\\\rangle $ as well.

This class implements all of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup, and
LOCA::Extended::MultiAbstractGroup methods for this extended set of
equations and therefore is a complete group which can be passed to
most NOX solvers to locate a single pitchfork or to the LOCA::Stepper
to compute a family of pitchforks in a second parameter.

The class is intialized via the pfParams parameter list argument to
the constructor. The parameters this class recognizes are:
\"Bifurcation Parameter\" -- [string] (Must be supplied) - Name of the
bifurcation parameter $p$

\"Antisymmetric Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>]
(Must be supplied) Vector storing antisymmetric vector $\\\\psi$

\"Bordered Solver Method\" -- [string] (default \"Bordering\") Method
for solving bordered systems of equations. See
LOCA::BorderedSolver::Factory for a description.

\"Symmetric Jacobian\" -- [bool] (default: false) - Flag indicating
whether Jacobian matrix $J$ is symmetric, in which case we force $a =
b$ and therefore the second tranpose solve for $w$ is unnecessary

\"Initial Null Vector Compuation\" -- [string] (default: \"User
Provided\") - Method to compute initial $a$ and $b$ vectors. Valid
choices are: \"User Provided\" - Initial vectors are provided in the
parameter list, in which case the following parameters are relevant:
\"Initial A Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be
supplied) - Vector storing initial value for $a$ vector

\"Initial B Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must be
supplied for nonsymmetric Jacobians) - Vector storing initial value
for $b$ vector

\"Solve df/dp\" - Compute $a = J^{-T}df/dp$ and $b = J^{-1} df/dp$
where $p$ is the bifurcation parameter.

\"Update Null Vectors Every Continuation Step\" -- [bool] (default:
true) - Flag indicating whether to update $a$ and $b$ vectors via $a =
w$ and $b = v$ every continuation step

\"Update Null Vectors Every Nonlinear Iteration\" -- [bool] (default:
false) - Flag indicating whether to update $a$ and $b$ vectors via $a
= w$ and $b = v$ every nonlinear iteration

C++ includes: LOCA_Pitchfork_MinimallyAugmented_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setX "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeX "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute and return solution vector, x, where this.x = grp.x + step *
d. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeF()

Compute extended continuation equations. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeJacobian()

Compute extended continuation jacobian. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeGradient()

Gradient is not defined for this system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction for extended continuation system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianTranspose
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianInverse
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isF "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isF() const

Return true if extended residual is valid. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isJacobian "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isGradient "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isNewton "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getX "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getX() const

Return extended solution vector. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getF "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getF() const

Return extended residual. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNormF "double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNormF() const

Return 2-norm of extended residual. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getGradient "const NOX::Abstract::Vector &
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getGradient()
const

Gradient is never valid. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getXPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getXPtr() const

Return RCP to extended solution vector. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getFPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getFPtr() const

Return RCP to extended residual. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getGradientPtr()
const

Gradient is never valid. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNewtonPtr()
const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual
"double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual()
const

Returns 2-norm of extended Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup
"Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup()
const

Return underlying group. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup
"Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::copy "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParamsMulti "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParams "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p (pVector = p). ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParam "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParam(int
paramID, double val)

Set parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParam "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getParams "const
LOCA::ParameterVector &
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getParams() const

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getParam "double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getParam(int
paramID) const

Return copy of parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getParam "double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::preProcessContinuationStep
"void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::postProcessContinuationStep
"void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::projectToDraw "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::projectToDrawDimension
"int
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::projectToDrawDimension()
const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeScaledDotProduct
"double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute a scaled dot product. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::printSolution "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::printSolution(const
double conParam) const

Function to print out solution and parameter after successful step. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::printSolution "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &x, const double conParam) const

Function to print out a vector and parameter after successful step. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::scaleVector "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::scaleVector(NOX::Abstract::Vector
&x) const

Scales a vector using scaling vector. ";

/*  Implementation of  */

/*  LOCA::BorderedSystem::AbstractGroup virtual methods

*/

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getBorderedWidth "int
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getBorderedWidth()
const

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getUnborderedGroup
"Teuchos::RCP< const NOX::Abstract::Group >
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getUnborderedGroup()
const

Get bottom-level unbordered group. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isCombinedAZero "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isCombinedAZero()
const

Indicates whether combined A block is zero. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isCombinedBZero "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isCombinedBZero()
const

Indicates whether combined B block is zero. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isCombinedCZero "bool
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::isCombinedCZero()
const

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::extractSolutionComponent
"void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::extractParameterComponent
"void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::loadNestedComponents
"void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillA "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillA(NOX::Abstract::MultiVector
&A) const

Fill the combined A block as described above. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillB "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillB(NOX::Abstract::MultiVector
&B) const

Fill the combined B block as described above. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillC "void
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const

Fill the combined C block as described above. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::ExtendedGroup "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &pfParams, const Teuchos::RCP<
LOCA::Pitchfork::MinimallyAugmented::AbstractGroup > &grp)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list.

pfParams:  [in] Parameter list determining the bordered solver method.

grp:  [in] Group representing $f$. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::ExtendedGroup "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::~ExtendedGroup "LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::~ExtendedGroup()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getBifParam "double
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup::getBifParam()
const

Get bifurcation parameter. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1ExtendedGroup.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::ExtendedGroup "

A group representing the Moore-Spence pitchfork equations.

The LOCA::Pitchfork::MooreSpence::ExtendedGroup is a concrete
implementation of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup and
LOCA::Extended::MultiAbstractGroup that defines the following extended
set of equations that are regular at a generic pitchfork bifurcation:
\\\\[ G(z) = \\\\left[ \\\\begin{array}{c} F(x,p)+\\\\sigma\\\\psi
\\\\\\\\ Jn \\\\\\\\ \\\\langle x, \\\\psi\\\\rangle \\\\\\\\ l^Tn-1
\\\\end{array} \\\\right] = 0 \\\\] where $z = [x, n, \\\\sigma,
p]\\\\in\\\\Re^{2n+2}$, $x$ is the solution vector, $n$ is the null
vector, $\\\\psi$ is the antisymmetric vector, $\\\\sigma$ is the
slack variables representing the asymmetry in the problem, $l$ is the
length normalization vector and $J$ is the Jacobian of F.

The group stores an underlying group of type
LOCA::Pitchfork::MooreSpence AbstractGroup to represent the equations
$F(x,p) = 0$ and to manipulate the underlying Jacobian $J$. Note that
the entire extended Jacobian $D_z G$ is not stored in memory, rather a
block-elimination algorithm (bordering algorithm) is used to compute
linear solves of the extended Jacobian (see
LOCA::Pitchfork::MooreSpence::SolverFactory for more details).

This class implements all of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup, and
LOCA::Extended::MultiAbstractGroup methods for this extended set of
equations and therefore is a complete group which can be passed to
most NOX solvers to locate a single pitchfork point or to the
LOCA::Stepper to compute a family of pitchforks in a second parameter.

However, Jacobian-tranpose operations and gradient calculations cannot
be implemented efficiently and therefore gradient-base nonlinear
solvers such as steepest descent and Trust region methods cannot be
used to solve the extended pitchfork equations.

The class is intialized via the pfParams parameter list argument to
the constructor. The parameters this class recognizes are:
\"Bifurcation Parameter\" -- [string] (Must be supplied) Name of the
bifurcation parameter $p$

\"Length Normalization Vector\" --
[Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Vector
storing length normalization vector $l$

\"Initial Null Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>] (Must
be supplied) Vector storing initial guess for the null vector $n$

\"Antisymmetric Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>]
(Must be supplied) Vector storing antisymmetric vector $\\\\psi$

\"Perturb Initial Solution\" -- [bool] (default: false) Flag
indicating whether to perturb the initial solution

\"Relative Perturbation Size\" -- [double] (default: 1.0e-3) Relative
perturbation size if perturbing the initial solution

C++ includes: LOCA_Pitchfork_MooreSpence_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::clone "Teuchos::RCP<
NOX::Abstract::Group >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setX "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector, x, to y. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeX "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute this.x = grp.x + step * d. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeF()

Compute the pitchfork equation residual $G$.

This method fills the extended residual \\\\[ G(z) = \\\\left[
\\\\begin{array}{c} F(x,p)+\\\\sigma\\\\psi \\\\\\\\ Jn \\\\\\\\
\\\\langle x, \\\\psi\\\\rangle \\\\\\\\ l^Tn-1 \\\\end{array}
\\\\right] = 0 \\\\] The solution component residual $F(x,p)$ and the
null-vector residual $Jn$ are calculated via the computeF and
applyJacobian methods of the underlying group. The symmetry residual,
$\\\\langle x, \\\\psi\\\\rangle$ is computed via the innerProduct()
method of the underlying group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeJacobian()

Compute the blocks of the Jacobian derivative of $G$.

This method computes the $J$, $\\\\partial F/\\\\partial p$, and
$\\\\partial Jn/\\\\partial p$ blocks of the extended Jacobian: \\\\[
D_z G(z) = \\\\begin{bmatrix} J & 0 & \\\\psi & \\\\frac{\\\\partial
F}{\\\\partial p} \\\\\\\\ \\\\frac{\\\\partial Jn}{\\\\partial x} & J
& 0 & \\\\frac{\\\\partial Jn}{\\\\partial p} \\\\\\\\ \\\\psi^T & 0 &
0 & 0 \\\\\\\\ 0 & l^T & 0 & 0 \\\\end{bmatrix} \\\\] by calling the
computeJacobian, computeDfDp, and computeDJnDp methods of the
underlying group. The second derivative matrix $\\\\partial
Jn/\\\\partial x$ is not calculated since only its action on vectors
is needed for linear solves using the bordering algorithms. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeGradient()

Gradient computation is not defined for this group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction using applyJacobianInverse(). ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Computes the extended Jacobian vector product.

This method computes the extended Jacobian vector product \\\\[
\\\\begin{bmatrix} J & 0 & \\\\psi & \\\\frac{\\\\partial
F}{\\\\partial p} \\\\\\\\ \\\\frac{\\\\partial Jn}{\\\\partial x} & J
& 0 & \\\\frac{\\\\partial Jn}{\\\\partial p} \\\\\\\\ \\\\psi^T & 0 &
0 & 0 \\\\\\\\ 0 & l^T & 0 & 0 \\\\end{bmatrix} \\\\begin{bmatrix} a
\\\\\\\\ b \\\\\\\\ c \\\\\\\\ d \\\\end{bmatrix} = \\\\begin{bmatrix}
Ja + \\\\psi c + \\\\frac{\\\\partial F}{\\\\partial p}d \\\\\\\\
\\\\frac{\\\\partial Jn}{\\\\partial x}a + Jb + \\\\frac{\\\\partial
Jn}{\\\\partial p}d \\\\\\\\ \\\\langle a, \\\\psi\\\\rangle \\\\\\\\
l^T b \\\\end{bmatrix} \\\\] using the applyJacobian, computeDJnDxa,
and innerProdut methods of the underlying group where $a$, $b$, $c$,
$d$ are the solution, null-vector, slack, and paramter components of
the given vector input. Vectors input and result must be of type
LOCA::Pitchfork::MooreSpence::ExtendedVector, otherwise an error is
thrown. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose product is not defined by this group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies the inverse of the extended Jacobian matrix using the
bordering algorithm.

This method is a special case of applyJacobianInverseMultiVector() for
a single right-hand-side. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system.

Uses a LOCA::Pitchfork::MooreSpence::SolverStrategy instantiated by
the LOCA::Pitchfork::MooreSpence::SolverFactory to implement the
solve. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isF "bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isF() const

Return true if the extended residual $G$ is valid. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isJacobian "bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isGradient "bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isNewton "bool
LOCA::Pitchfork::MooreSpence::ExtendedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getX "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getX() const

Return extended solution vector $z$. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getF "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getF() const

Return extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormF "double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormF() const

Return 2-norm of $G(z)$. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradient() const

Vector returned is not valid. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getXPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getXPtr() const

Return RCP to extended solution vector $z$. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getFPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getFPtr() const

Return RCP to extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getGradientPtr() const

Vector returned is not valid. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNewtonPtr() const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual
"double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()
const

Return the norm of the Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getUnderlyingGroup()
const

Return underlying group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::copy "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParamsMulti "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParams "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParam "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParam(int paramID,
double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParam "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParams "const
LOCA::ParameterVector &
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParams() const

Return a const reference to the paramter vector owned by the group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParam "double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParam(int paramID)
const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParam "double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::preProcessContinuationStep
"void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::postProcessContinuationStep
"void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::projectToDraw "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::projectToDrawDimension "int
LOCA::Pitchfork::MooreSpence::ExtendedGroup::projectToDrawDimension()
const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution(const
double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &x_, const double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getBifParam "double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::getBifParam() const

Get bifurcation parameter. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::lTransNorm "double
LOCA::Pitchfork::MooreSpence::ExtendedGroup::lTransNorm(const
NOX::Abstract::Vector &n) const

Defines null vector normalization $l^Tn$ equation. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::lTransNorm "void
LOCA::Pitchfork::MooreSpence::ExtendedGroup::lTransNorm(const
NOX::Abstract::MultiVector &n, NOX::Abstract::MultiVector::DenseMatrix
&result) const

null vector normalization for multivectors

Note: result should have 1 row and n.numVectors() columns. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::ExtendedGroup "LOCA::Pitchfork::MooreSpence::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &pfParams, const Teuchos::RCP<
LOCA::Pitchfork::MooreSpence::AbstractGroup > &g)

Constructor with initial data passed through parameter lists. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::ExtendedGroup "LOCA::Pitchfork::MooreSpence::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedGroup::~ExtendedGroup "LOCA::Pitchfork::MooreSpence::ExtendedGroup::~ExtendedGroup()

Destructor. ";


// File: classLOCA_1_1Hopf_1_1MinimallyAugmented_1_1ExtendedGroup.xml
%feature("docstring") LOCA::Hopf::MinimallyAugmented::ExtendedGroup "

A group representing the minimally augemented Hopf equations.

The LOCA::Hopf::MinimallyAugmented::ExtendedGroup is a concrete
implementation of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup and
LOCA::Extended::MultiAbstractGroup that defines the following extended
set of equations that are regular at a generic Hopf: \\\\[ G(z) =
\\\\left[ \\\\begin{array}{c} F(x,p)\\\\\\\\ \\\\sigma \\\\end{array}
\\\\right] = 0 \\\\] where $z = [x, p, \\\\omega]\\\\in\\\\Re^{n+2}$,
$x$ is the solution vector, $p$ is the bifurcation parameter,
$\\\\omega$ is the Hopf frequency and $\\\\sigma\\\\in\\\\Re$ is a
measure of the singularity of $J+i\\\\omega M$ and is defined via
\\\\[ \\\\begin{bmatrix} J+i\\\\omega M & a \\\\\\\\ b^H & 0
\\\\end{bmatrix} \\\\begin{bmatrix} v \\\\\\\\ \\\\sigma_1
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\begin{bmatrix} J^T-i\\\\omega M^T & b \\\\\\\\ a^H & 0
\\\\end{bmatrix} \\\\begin{bmatrix} w \\\\\\\\ \\\\sigma_2
\\\\end{bmatrix} = \\\\begin{bmatrix} 0 \\\\\\\\ n \\\\end{bmatrix},
\\\\] \\\\[ \\\\sigma = w^H J+i\\\\omega M v/n \\\\] for any vectors
$a$ and $b$ in $\\\\mathbb{C}^n$. Using these relationships, it is
easy to show \\\\[ \\\\begin{split} \\\\sigma_x &= (w^H(J+i\\\\omega
M)v)_x/n = w^H(J+i\\\\omega M)_x v/n \\\\\\\\ \\\\sigma_p &=
(w^H(J+i\\\\omega M)v)_p/n = w^H(J+i\\\\omega M)_p v/n \\\\end{split}
\\\\]

The group stores an underlying group of type
LOCA::Hopf::MinimallyAugmented::AbstractGroup to represent the
equations $F(x,p) = 0$ and to manipulate the underlying complex matrix
$C = J+i\\\\omega M$. This interface defines methods for computing the
derivatives $(w^H C v)_x$ and $(w^H C v)_p$. Since LOCA is not able to
deal with complex vectors and matrices directly, the real-equivalent
formulation is used for all complex calculations.

This class implements all of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup, and
LOCA::Extended::MultiAbstractGroup methods for this extended set of
equations and therefore is a complete group which can be passed to
most NOX solvers to locate a single pitchfork or to the LOCA::Stepper
to compute a family of pitchforks in a second parameter.

The class is intialized via the hpfParams parameter list argument to
the constructor. The parameters this class recognizes are:
\"Bifurcation Parameter\" -- [string] (Must be supplied) - Name of the
bifurcation parameter $p$

\"Initial Frequency\" -- [double] (Must be supplied) Initial guess for
the Hopf frequency $\\\\omega$.

\"Bordered Solver Method\" -- [string] (default \"Bordering\") Method
for solving bordered systems of equations. See
LOCA::BorderedSolver::Factory for a description.

\"Symmetric Jacobian\" -- [bool] (default: false) - Flag indicating
whether Jacobian matrix $J$ is symmetric, in which case we force $a =
b$ and therefore the second tranpose solve for $w$ is unnecessary

\"Initial Null Vector Compuation\" -- [string] (default: \"User
Provided\") - Method to compute initial $a$ and $b$ vectors. Valid
choices are: \"User Provided\" - Initial vectors are provided in the
parameter list, in which case the following parameters are relevant:
\"Initial Real A Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>]
(Must be supplied) - Vector storing initial value for the real
component of the $a$ vector

\"Initial Imaginary A Vector\" --
[Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) - Vector
storing initial value for the imaginary component of the $a$ vector

\"Initial Real B Vector\" -- [Teuchos::RCP<NOX::Abstract::Vector>]
(Must be supplied for nonsymmetric Jacobians) - Vector storing initial
value for the real component of the $b$ vector

\"Initial Imaginary B Vector\" --
[Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) - Vector
storing initial value for the imaginary component of the $b$ vector

\"Solve df/dp\" - Compute $a = J^{-T}df/dp$ and $b = J^{-1} df/dp$
where $p$ is the bifurcation parameter.

\"Update Null Vectors Every Continuation Step\" -- [bool] (default:
true) - Flag indicating whether to update $a$ and $b$ vectors via $a =
w$ and $b = v$ every continuation step

\"Update Null Vectors Every Nonlinear Iteration\" -- [bool] (default:
false) - Flag indicating whether to update $a$ and $b$ vectors via $a
= w$ and $b = v$ every nonlinear iteration

C++ includes: LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::clone "Teuchos::RCP<
NOX::Abstract::Group >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setX "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeX "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute and return solution vector, x, where this.x = grp.x + step *
d. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeF()

Compute extended continuation equations. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeJacobian()

Compute extended continuation jacobian. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeGradient()

Gradient is not defined for this system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction for extended continuation system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianTranspose
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isF "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isF() const

Return true if extended residual is valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isJacobian "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isGradient "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isNewton "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getX "const
NOX::Abstract::Vector &
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getX() const

Return extended solution vector. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getF "const
NOX::Abstract::Vector &
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getF() const

Return extended residual. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNormF "double
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNormF() const

Return 2-norm of extended residual. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getGradient() const

Gradient is never valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getXPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getXPtr() const

Return RCP to extended solution vector. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getFPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getFPtr() const

Return RCP to extended residual. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getGradientPtr() const

Gradient is never valid. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNewtonPtr() const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual
"double
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getNormNewtonSolveResidual()
const

Returns 2-norm of extended Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup()
const

Return underlying group. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::copy "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParamsMulti "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParams "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p (pVector = p). ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParam "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParam(int paramID,
double val)

Set parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParam "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::setParam(std::string
paramID, double val)

Set parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getParams "const
LOCA::ParameterVector &
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getParams() const

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getParam "double
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getParam(int paramID)
const

Return copy of parameter indexed by (integer) paramID. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getParam "double
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getParam(std::string
paramID) const

Return copy of parameter indexed by (std::string) paramID. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::preProcessContinuationStep
"void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::postProcessContinuationStep
"void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::projectToDraw "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::projectToDrawDimension
"int
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::projectToDrawDimension()
const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeScaledDotProduct
"double
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::computeScaledDotProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute a scaled dot product. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::printSolution "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::printSolution(const
double conParam) const

Function to print out solution and parameter after successful step. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::printSolution "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &x, const double conParam) const

Function to print out a vector and parameter after successful step. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::scaleVector "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::scaleVector(NOX::Abstract::Vector
&x) const

Scales a vector using scaling vector. ";

/*  Implementation of  */

/*  LOCA::BorderedSystem::AbstractGroup virtual methods

*/

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getBorderedWidth "int
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getBorderedWidth()
const

Return the total width of the bordered rows/columns. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getUnborderedGroup "Teuchos::RCP< const NOX::Abstract::Group >
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getUnborderedGroup()
const

Get bottom-level unbordered group. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isCombinedAZero "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isCombinedAZero() const

Indicates whether combined A block is zero. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isCombinedBZero "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isCombinedBZero() const

Indicates whether combined B block is zero. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isCombinedCZero "bool
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::isCombinedCZero() const

Indicates whether combined C block is zero. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::extractSolutionComponent
"void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::extractSolutionComponent(const
NOX::Abstract::MultiVector &v, NOX::Abstract::MultiVector &v_x) const

Given the vector v, extract the underlying solution component
corresponding to the unbordered group. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::extractParameterComponent
"void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::extractParameterComponent(bool
use_transpose, const NOX::Abstract::MultiVector &v,
NOX::Abstract::MultiVector::DenseMatrix &v_p) const

Given the vector v, extract the parameter components of all of the
nested subvectors in v down to the solution component for the
unbordered group. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::loadNestedComponents "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::loadNestedComponents(const
NOX::Abstract::MultiVector &v_x, const
NOX::Abstract::MultiVector::DenseMatrix &v_p,
NOX::Abstract::MultiVector &v) const

Given the solution component v_x and combined parameter components
v_p, distribute these components through the nested sub-vectors in v.
";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::fillA "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::fillA(NOX::Abstract::MultiVector
&A) const

Fill the combined A block as described above. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::fillB "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::fillB(NOX::Abstract::MultiVector
&B) const

Fill the combined B block as described above. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::fillC "void
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::fillC(NOX::Abstract::MultiVector::DenseMatrix
&C) const

Fill the combined C block as described above. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::ExtendedGroup "LOCA::Hopf::MinimallyAugmented::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &hpfParams, const Teuchos::RCP<
LOCA::Hopf::MinimallyAugmented::AbstractGroup > &grp)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list.

hpfParams:  [in] Parameter list determining the bordered solver
method.

grp:  [in] Group representing $f$. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::ExtendedGroup "LOCA::Hopf::MinimallyAugmented::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::~ExtendedGroup "LOCA::Hopf::MinimallyAugmented::ExtendedGroup::~ExtendedGroup()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getBifParam "double
LOCA::Hopf::MinimallyAugmented::ExtendedGroup::getBifParam() const

Get bifurcation parameter. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1ExtendedGroup.xml
%feature("docstring") LOCA::Hopf::MooreSpence::ExtendedGroup "

A group representing the Moore-Spence Hopf equations.

The LOCA::Hopf::MooreSpence::ExtendedGroup is a concrete
implementation of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup and
LOCA::Extended::MultiAbstractGroup that defines the following extended
set of equations that are regular at a generic Hopf point: \\\\[ G(z)
= \\\\left[ \\\\begin{array}{c} F(x,p) \\\\\\\\ Jy-\\\\omega Bz
\\\\\\\\ Jz+\\\\omega By \\\\\\\\ l^Ty-1 \\\\\\\\ l^Tz \\\\end{array}
\\\\right] = 0 \\\\] where $z = [x, y, z, \\\\omega,
p]\\\\in\\\\Re^{3n+2}$, $x$ is the solution vector, $y+i\\\\omega z$
is the complex eigenvector of $J$ with corresponding eigenvalues
$\\\\pm i\\\\omega$, $l$ is the length normalization vector and $J$ is
the Jacobian of F w.r.t $x$.

The group stores an underlying group of type LOCA::Hopf::MooreSpence
AbstractGroup to represent the equations $F(x,p) = 0$ and to
manipulate the underlying Jacobian $J$. Note that the entire extended
Jacobian $D_z G$ is not stored in memory, rather a block-elimination
algorithm (bordering algorithm) is used to compute linear solves of
the extended Jacobian (see LOCA::Hopf::MooreSpence::SolverFactory()
for more details).

This class implements all of the NOX::Abstract::Group,
LOCA::MultiContinuation::AbstractGroup, and
LOCA::Extended::MultiAbstractGroup methods for this extended set of
equations and therefore is a complete group which can be passed to
most NOX solvers to locate a single Hopf point or to the LOCA::Stepper
to compute a family of Hopf points in a second parameter.

However, Jacobian-tranpose operations and gradient calculations cannot
be implemented efficiently and therefore gradient-base nonlinear
solvers such as steepest descent and Trust region methods cannot be
used to solve the extended Hopf point equations.

The class is intialized via the hopfParams parameter list argument to
the constructor. The parameters this class recognizes are:
\"Bifurcation Parameter\" -- [string] (Must be supplied) Name of the
bifurcation parameter $p$

\"Length Normalization Vector\" --
[Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Vector
storing length normalization vector $l$

\"Initial Real Eigenvector\" -- [Teuchos::RCP<NOX::Abstract::Vector>]
(Must be supplied) Vector storing initial guess for the real component
of the eigenvector $y$

\"Initial Imaginary Eigenvector\" --
[Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) Vector
storing initial guess for the imaginary component of the eigenvector
$z$

\"Initial Frequency\" -- [double] (Must be supplied) Initial guess for
the Hopf frequency $\\\\omega$.

\"Perturb Initial Solution\" -- [bool] (default: false) Flag
indicating whether to perturb the initial solution

\"Relative Perturbation Size\" -- [double] (default: 1.0e-3) Relative
perturbation size if perturbing the initial solution

C++ includes: LOCA_Hopf_MooreSpence_ExtendedGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::Hopf::MooreSpence::ExtendedGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::setX "void LOCA::Hopf::MooreSpence::ExtendedGroup::setX(const
NOX::Abstract::Vector &y)

Set the solution vector, x, to y. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::computeX "void
LOCA::Hopf::MooreSpence::ExtendedGroup::computeX(const
NOX::Abstract::Group &g, const NOX::Abstract::Vector &d, double step)

Compute this.x = grp.x + step * d. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::computeF "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeF()

Compute the Hopf point equation residual $G$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeJacobian()

Compute the blocks of the Jacobian derivative of $G$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeGradient()

Gradient computation is not defined for this group. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeNewton(Teuchos::ParameterList
&params)

Compute Newton direction using applyJacobianInverse(). ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobian(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Computes the extended Jacobian vector product. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Jacobian transpose product is not defined by this group. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies the inverse of the extended Jacobian matrix using the
bordering algorithm.

This method is a special case of applyJacobianInverseMultiVector() for
a single right-hand-side. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for extended system. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Jacobian transpose not defined for this system. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for extended system.

Uses a LOCA::Hopf::MooreSpence::SolverStrategy instantiated by the
LOCA::Hopf::MooreSpence::SolverFactory to implement the solve. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::isF "bool LOCA::Hopf::MooreSpence::ExtendedGroup::isF() const

Return true if the extended residual $G$ is valid. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::isJacobian "bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isJacobian() const

Return true if the extended Jacobian is valid. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::isGradient "bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isGradient() const

Always returns false. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::isNewton "bool
LOCA::Hopf::MooreSpence::ExtendedGroup::isNewton() const

Return true if the extended Newton direction is valid. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::getX "const NOX::Abstract::Vector &
LOCA::Hopf::MooreSpence::ExtendedGroup::getX() const

Return extended solution vector $z$. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::getF "const NOX::Abstract::Vector &
LOCA::Hopf::MooreSpence::ExtendedGroup::getF() const

Return extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getNormF "double
LOCA::Hopf::MooreSpence::ExtendedGroup::getNormF() const

Return 2-norm of $G(z)$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getGradient "const
NOX::Abstract::Vector &
LOCA::Hopf::MooreSpence::ExtendedGroup::getGradient() const

Vector returned is not valid. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getNewton "const
NOX::Abstract::Vector &
LOCA::Hopf::MooreSpence::ExtendedGroup::getNewton() const

Return extended Newton direction. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::getXPtr
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getXPtr() const

Return RCP to extended solution vector $z$. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::getFPtr
"Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getFPtr() const

Return RCP to extended equation residual $G(z)$. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getGradientPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getGradientPtr() const

Vector returned is not valid. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getNewtonPtr "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedGroup::getNewtonPtr() const

Return RCP to extended Newton direction. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual "double
LOCA::Hopf::MooreSpence::ExtendedGroup::getNormNewtonSolveResidual()
const

Return the norm of the Newton solve residual. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Hopf::MooreSpence::ExtendedGroup::getUnderlyingGroup() const

Return underlying group. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Hopf::MooreSpence::ExtendedGroup::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedGroup::copy "void LOCA::Hopf::MooreSpence::ExtendedGroup::copy(const
NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::setParamsMulti "void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParamsMulti(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::setParams "void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParams(const
ParameterVector &p)

Set the parameter vector in the group to p. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::setParam "void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParam(int paramID, double
val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::setParam "void
LOCA::Hopf::MooreSpence::ExtendedGroup::setParam(std::string paramID,
double val)

Set parameter indexed by paramID. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getParams "const
LOCA::ParameterVector &
LOCA::Hopf::MooreSpence::ExtendedGroup::getParams() const

Return a const reference to the paramter vector owned by the group. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getParam "double
LOCA::Hopf::MooreSpence::ExtendedGroup::getParam(int paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getParam "double
LOCA::Hopf::MooreSpence::ExtendedGroup::getParam(std::string paramID)
const

Return copy of parameter indexed by paramID. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::ExtendedGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::preProcessContinuationStep "void
LOCA::Hopf::MooreSpence::ExtendedGroup::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::postProcessContinuationStep "void
LOCA::Hopf::MooreSpence::ExtendedGroup::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::projectToDraw "void
LOCA::Hopf::MooreSpence::ExtendedGroup::projectToDraw(const
NOX::Abstract::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::projectToDrawDimension "int
LOCA::Hopf::MooreSpence::ExtendedGroup::projectToDrawDimension() const

Returns the dimension of the project to draw array. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution "void
LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution(const double
conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution "void
LOCA::Hopf::MooreSpence::ExtendedGroup::printSolution(const
NOX::Abstract::Vector &x_, const double conParam) const

Function to print out extended solution and continuation parameter
after successful continuation step.

This method prints the solution, null-vector, and parameter components
of the extended solution vector using the printSolution method of the
underlying group. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getBifParam "double
LOCA::Hopf::MooreSpence::ExtendedGroup::getBifParam() const

Get bifurcation parameter. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::getFrequency "double
LOCA::Hopf::MooreSpence::ExtendedGroup::getFrequency() const

Get bifurcation frequency. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::lTransNorm "double
LOCA::Hopf::MooreSpence::ExtendedGroup::lTransNorm(const
NOX::Abstract::Vector &z) const

Defines null vector normalization $l^Tz$ equation. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::lTransNorm "void
LOCA::Hopf::MooreSpence::ExtendedGroup::lTransNorm(const
NOX::Abstract::MultiVector &z, NOX::Abstract::MultiVector::DenseMatrix
&result) const

null vector normalization for multivectors

Note: result should have 1 row and z.numVectors() columns. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::ExtendedGroup "LOCA::Hopf::MooreSpence::ExtendedGroup::ExtendedGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &hopfParams, const Teuchos::RCP<
LOCA::Hopf::MooreSpence::AbstractGroup > &g)

Constructor with initial data passed through parameter lists. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::ExtendedGroup "LOCA::Hopf::MooreSpence::ExtendedGroup::ExtendedGroup(const
ExtendedGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedGroup::~ExtendedGroup "LOCA::Hopf::MooreSpence::ExtendedGroup::~ExtendedGroup()

Destructor. ";


// File: classLOCA_1_1MultiContinuation_1_1ExtendedMultiVector.xml
%feature("docstring") LOCA::MultiContinuation::ExtendedMultiVector "

MultiVector class to hold solution vectors, Newton vectors, etc. for
continuation equations.

This class uses the LOCA::Extended::MultiVector implementation to
store the solution and parameter components of the continuation vector
and merely provides an interface for naming which components of the
multivector these quantities correspond to.

C++ includes: LOCA_MultiContinuation_ExtendedMultiVector.H ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &xVec, int nColumns, int nScalarRows,
NOX::CopyType type=NOX::DeepCopy)

Constructor.

Generates a multivector with nColumns from xVec amd nScalarRows of
zeros. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::MultiVector &xVec, int nScalarRows)

Constructor.

Initializes the scalar matrix to nScalarRows rows and
xVec.numVectors() columns of zeros ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::MultiVector &xVec, const
NOX::Abstract::MultiVector::DenseMatrix &params)

Constructor.

Sets the scalar matrix explicitly ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, const std::vector< int > &index, bool
view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::~ExtendedMultiVector "LOCA::MultiContinuation::ExtendedMultiVector::~ExtendedMultiVector()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::MultiContinuation::ExtendedMultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new multi-vector of the same underlying type by cloning
\"this\", and return a pointer to the new vector. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::MultiContinuation::ExtendedMultiVector::clone(int numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::subCopy "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::MultiContinuation::ExtendedMultiVector::subCopy(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::subView "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::MultiContinuation::ExtendedMultiVector::subView(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec() const

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec()

Returns the solution vector component of extended multivector. ";


// File: classLOCA_1_1PhaseTransition_1_1ExtendedMultiVector.xml
%feature("docstring") LOCA::PhaseTransition::ExtendedMultiVector "

MultiVector class to hold solution vectors, Newton vectors, etc. for
the phase transition tracking algorithm.

This class uses the LOCA::Extended::MultiVector implementation to
store the solution1, solution2, and parameter components of the phase
transition multivector and merely provides an interface for naming
which components of the multivector these quantities correspond to.

C++ includes: LOCA_PhaseTransition_ExtendedMultiVector.H ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector "LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &cloneVec, int nColumns)

Constructor.

Generates a multivector with nColumns columns from cloneVec ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector "LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::MultiVector &xVec, const NOX::Abstract::MultiVector
&nullVec, const NOX::Abstract::MultiVector::DenseMatrix &bifParams)

Constructor.

Construct the multivector from xVec, nullVec, and bifParams ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector "LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector "LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector "LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, const std::vector< int > &index, bool
view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::~ExtendedMultiVector "LOCA::PhaseTransition::ExtendedMultiVector::~ExtendedMultiVector()

Destructor. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new multi-vector of the same underlying type by cloning
\"this\", and return a pointer to the new vector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::clone(int numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::subCopy "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::subCopy(const std::vector<
int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::subView "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::subView(const std::vector<
int > &index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::getXMultiVec() const

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::getXMultiVec()

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::getNullMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::getNullMultiVec() const

Returns the null vector component of extended multivector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::getNullMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::PhaseTransition::ExtendedMultiVector::getNullMultiVec()

Returns the null vector component of extended multivector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::getColumn "Teuchos::RCP<
LOCA::PhaseTransition::ExtendedVector >
LOCA::PhaseTransition::ExtendedMultiVector::getColumn(int i)

Returns ith column as an extended vector. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedMultiVector::getColumn "Teuchos::RCP<
const LOCA::PhaseTransition::ExtendedVector >
LOCA::PhaseTransition::ExtendedMultiVector::getColumn(int i) const

Returns ith column as an extended vector. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1ExtendedMultiVector.xml
%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector "

MultiVector class to hold solution vectors, Newton vectors, etc.for
the Moore-Spence turning point formulation.

This class uses the LOCA::Extended::MultiVector implementation to
store the solution, null, and parameter components of the turning
point multivector and merely provides an interface for naming which
components of the multivector these quantities correspond to.

C++ includes: LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &cloneVec, int nColumns)

Constructor.

Generates a multivector with nColumns columns from cloneVec ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::MultiVector &xVec, const NOX::Abstract::MultiVector
&nullVec, const NOX::Abstract::MultiVector::DenseMatrix &bifParams)

Constructor.

Construct the multivector from xVec, nullVec, and bifParams ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, const std::vector< int > &index, bool
view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector
"LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new multi-vector of the same underlying type by cloning
\"this\", and return a pointer to the new vector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone(int
numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subCopy "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subCopy(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subView "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subView(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec()
const

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec()

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec
"Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec()
const

Returns the null vector component of extended multivector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec
"Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec()

Returns the null vector component of extended multivector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn "Teuchos::RCP< LOCA::TurningPoint::MooreSpence::ExtendedVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn(int i)

Returns ith column as an extended vector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn "Teuchos::RCP< const LOCA::TurningPoint::MooreSpence::ExtendedVector >
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn(int i)
const

Returns ith column as an extended vector. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1ExtendedMultiVector.xml
%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector "

MultiVector class to hold solution vectors, Newton vectors, etc.for
the Moore-Spence pitchfork formulation.

This class uses the LOCA::Extended::MultiVector implementation to
store the solution, null, parameter, and slack components of the
pitchfork multivector and merely provides an interface for naming
which components of the multivector these quantities correspond to.

C++ includes: LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &cloneVec, int nColumns)

Constructor.

Generates a multivector with nColumns columns from cloneVec ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::MultiVector &xVec, const NOX::Abstract::MultiVector
&nullVec, const NOX::Abstract::MultiVector::DenseMatrix &slacks, const
NOX::Abstract::MultiVector::DenseMatrix &bifParams)

Constructor.

Construct the multivector from xVec, nullVec, and bifParams ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector
"LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, const std::vector< int > &index, bool
view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector
"LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new multi-vector of the same underlying type by cloning
\"this\", and return a pointer to the new vector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::clone(int numvecs)
const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::subCopy "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::subCopy(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::subView "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::subView(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getXMultiVec()
const

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getXMultiVec()

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getNullMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getNullMultiVec()
const

Returns the null vector component of extended multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getNullMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getNullMultiVec()

Returns the null vector component of extended multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getSlacks "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getSlacks() const

Returns slack component of the extended multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getSlacks "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getSlacks()

Returns slack component of the extended multivector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getBifParams "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getBifParams()
const

Returns bifurcation parameter component of the extended multivector.
";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getBifParams "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getBifParams()

Returns bifurcation parameter component of the extended multivector.
";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getColumn "Teuchos::RCP< LOCA::Pitchfork::MooreSpence::ExtendedVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getColumn(int i)

Returns ith column as an extended vector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getColumn "Teuchos::RCP< const LOCA::Pitchfork::MooreSpence::ExtendedVector >
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getColumn(int i)
const

Returns ith column as an extended vector. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1ExtendedMultiVector.xml
%feature("docstring") LOCA::Hopf::MooreSpence::ExtendedMultiVector "

Multi-vector class to hold solution vectors, Newton vectors, etc.for
the Moore-Spence Hopf eqautions.

This class uses the LOCA::Extended::MultiVector implementation to
store the solution, real and imaginary eigenvector, frequency and
parameter components of the Hopf multi vector and merely provides an
interface for naming which components of the multivector these
quantities correspond to.

C++ includes: LOCA_Hopf_MooreSpence_ExtendedMultiVector.H ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector "LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &cloneVec, int nColumns)

Constructor.

Generates a multivector with nColumns columns from cloneVec ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector "LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::MultiVector &xVec, const NOX::Abstract::MultiVector
&realEigenVec, const NOX::Abstract::MultiVector &imagEigenVec, const
NOX::Abstract::MultiVector::DenseMatrix &freqs, const
NOX::Abstract::MultiVector::DenseMatrix &bifParams)

Constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector "LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector "LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector "LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(const
ExtendedMultiVector &source, const std::vector< int > &index, bool
view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector "LOCA::Hopf::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Create a new multi-vector of the same underlying type by cloning
\"this\", and return a pointer to the new vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::clone "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::clone(int numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::subCopy "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::subCopy(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::subView "Teuchos::RCP<
NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::subView(const
std::vector< int > &index) const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getXMultiVec() const

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getXMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getXMultiVec()

Returns the solution vector component of extended multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getRealEigenMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getRealEigenMultiVec()
const

Returns the real part of the eigenvector component of extended
multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getRealEigenMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getRealEigenMultiVec()

Returns the real part of the eigenvector component of extended
multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getImagEigenMultiVec "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getImagEigenMultiVec()
const

Returns the imaginary part of the eigenvector component of extended
multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getImagEigenMultiVec "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getImagEigenMultiVec()

Returns the imaginary part of the eigenvector component of extended
multivector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getFrequencies "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getFrequencies() const

Returns frequency component of extended multi vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getFrequencies "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getFrequencies()

Returns frequency component of extended multi vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getBifParams "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getBifParams() const

Returns bifurcation parameter component of extended multi vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getBifParams "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getBifParams()

Returns bifurcation parameter component of extended multi vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getColumn "Teuchos::RCP< LOCA::Hopf::MooreSpence::ExtendedVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getColumn(int i)

Returns ith column as an extended vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getColumn "Teuchos::RCP< const LOCA::Hopf::MooreSpence::ExtendedVector >
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getColumn(int i) const

Returns ith column as an extended vector. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1ExtendedVector.xml
%feature("docstring") LOCA::Hopf::MooreSpence::ExtendedVector "

Vector class to hold solution vectors, Newton vectors, etc. for Moore-
Spence Hopf equations.

This class uses the LOCA::Extended::Vector implementation to store the
solution, real and imaginary eigenvector, frequency and parameter
components of the Hopf vector and merely provides an interface for
naming which components of the multivector these quantities correspond
to.

C++ includes: LOCA_Hopf_MooreSpence_ExtendedVector.H ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector "LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector
&realEigenVec, const NOX::Abstract::Vector &imagEigenVec, double
frequency, double bifParam)

Constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector "LOCA::Hopf::MooreSpence::ExtendedVector::ExtendedVector(const
ExtendedVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::~ExtendedVector "LOCA::Hopf::MooreSpence::ExtendedVector::~ExtendedVector()

Destructor. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedVector::clone
"Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::ExtendedVector::setVec
"void LOCA::Hopf::MooreSpence::ExtendedVector::setVec(const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector
&realEigenVec, const NOX::Abstract::Vector &imagEigenVec, double
frequency, double bifPar)

Sets the Hopf vector by setting its five components. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec "Teuchos::RCP< const
NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec() const

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec() const

Returns the real part of the eigenvector component of extended vector.
";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec() const

Returns the imaginary part of the eigenvector component of extended
vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency "double
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency() const

Returns the frequency component of the extended vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam "double
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam() const

Get Bifurcation parameter. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec "Teuchos::RCP<
NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::getXVec()

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::getRealEigenVec()

Returns the real part of the eigenvector component of extended vector.
";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Hopf::MooreSpence::ExtendedVector::getImagEigenVec()

Returns the imaginary part of the eigenvector component of extended
vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency "double &
LOCA::Hopf::MooreSpence::ExtendedVector::getFrequency()

Returns the frequency component of the extended vector. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam "double &
LOCA::Hopf::MooreSpence::ExtendedVector::getBifParam()

Get Bifurcation parameter. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1ExtendedVector.xml
%feature("docstring") LOCA::TurningPoint::MooreSpence::ExtendedVector
"

Vector class to hold solution vectors, Newton vectors, etc. for the
Moore-Spence turning point formulation.

This class uses the LOCA::Extended::Vector implementation to store the
solution, null, and parameter components of the turning point vector
and merely provides an interface for naming which components of the
vector these quantities correspond to.

C++ includes: LOCA_TurningPoint_MooreSpence_ExtendedVector.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector "LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector &nullVec,
double bifParam)

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector "LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(const
ExtendedVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::~ExtendedVector "LOCA::TurningPoint::MooreSpence::ExtendedVector::~ExtendedVector()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::clone "Teuchos::RCP<
NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::setVec "void
LOCA::TurningPoint::MooreSpence::ExtendedVector::setVec(const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector &nullVec,
double bifPar)

Sets the Vector by setting its three components. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec() const

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec() const

Returns the null vector component of extended vector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam "double
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam() const

Get Bifurcation parameter. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec()

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec()

Returns the null vector component of extended vector. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam "double
& LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam()

Get Bifurcation parameter. ";


// File: classLOCA_1_1MultiContinuation_1_1ExtendedVector.xml
%feature("docstring") LOCA::MultiContinuation::ExtendedVector "

Vector class to hold solution vectors, Newton vectors, etc. for
continuation equations.

This class uses the LOCA::Extended::Vector implementation to store the
solution and parameter components of the continuation vector and
merely provides an interface for naming which components of the
multivector these quantities correspond to.

C++ includes: LOCA_MultiContinuation_ExtendedVector.H ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedVector::ExtendedVector "LOCA::MultiContinuation::ExtendedVector::ExtendedVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &xVec, int nScalars)

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedVector::ExtendedVector "LOCA::MultiContinuation::ExtendedVector::ExtendedVector(const
ExtendedVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedVector::~ExtendedVector "LOCA::MultiContinuation::ExtendedVector::~ExtendedVector()

Destructor. ";

%feature("docstring")  LOCA::MultiContinuation::ExtendedVector::clone
"Teuchos::RCP< NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Assignment operator. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedVector::getXVec "Teuchos::RCP< const
NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedVector::getXVec() const

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::MultiContinuation::ExtendedVector::getXVec "Teuchos::RCP<
NOX::Abstract::Vector >
LOCA::MultiContinuation::ExtendedVector::getXVec()

Returns the solution vector component of extended vector. ";


// File: classLOCA_1_1PhaseTransition_1_1ExtendedVector.xml
%feature("docstring") LOCA::PhaseTransition::ExtendedVector "

Vector class to hold solution vectors, Newton vectors, etc. for the
Phase Transition tracking formulation.

This class uses the LOCA::Extended::Vector implementation to store the
solution1, solution2, and parameter components of the phase transition
vector and merely provides an interface for naming which components of
the vector these quantities correspond to.

C++ includes: LOCA_PhaseTransition_ExtendedVector.H ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedVector::ExtendedVector "LOCA::PhaseTransition::ExtendedVector::ExtendedVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &x1Vec, const NOX::Abstract::Vector &x2Vec,
double ptp)

Constructor. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedVector::ExtendedVector "LOCA::PhaseTransition::ExtendedVector::ExtendedVector(const
ExtendedVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::PhaseTransition::ExtendedVector::~ExtendedVector "LOCA::PhaseTransition::ExtendedVector::~ExtendedVector()

Destructor. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::clone "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::setVec "void LOCA::PhaseTransition::ExtendedVector::setVec(const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector &nullVec,
double bifPar)

Sets the Vector by setting its three components. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::X1 "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedVector::X1() const

Returns the solution1 vector component of extended vector. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::X2 "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedVector::X2() const

Returns the solution2 vector component of extended vector. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::PTP "double LOCA::PhaseTransition::ExtendedVector::PTP() const

Get Bifurcation parameter. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::X1 "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedVector::X1()

Returns the solution vector component of extended vector. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::X2 "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::PhaseTransition::ExtendedVector::X2()

Returns the null vector component of extended vector. ";

%feature("docstring")  LOCA::PhaseTransition::ExtendedVector::PTP "double & LOCA::PhaseTransition::ExtendedVector::PTP()

Get Bifurcation parameter. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1ExtendedVector.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::ExtendedVector "

Vector class to hold solution vectors, Newton vectors, etc. for the
Moore-Spence turning point formulation.

This class uses the LOCA::Extended::Vector implementation to store the
solution, null, and parameter components of the turning point vector
and merely provides an interface for naming which components of the
vector these quantities correspond to.

C++ includes: LOCA_Pitchfork_MooreSpence_ExtendedVector.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector "LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector &nullVec,
double slack, double bifParam)

Constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector "LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector(const
ExtendedVector &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::~ExtendedVector "LOCA::Pitchfork::MooreSpence::ExtendedVector::~ExtendedVector()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::clone "Teuchos::RCP<
NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedVector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::setVec "void
LOCA::Pitchfork::MooreSpence::ExtendedVector::setVec(const
NOX::Abstract::Vector &xVec, const NOX::Abstract::Vector &nullVec,
double slack, double bifPar)

Sets the Vector by setting its three components. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getXVec "Teuchos::RCP<
const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedVector::getXVec() const

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getNullVec "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedVector::getNullVec() const

Returns the null vector component of extended vector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getSlack "double
LOCA::Pitchfork::MooreSpence::ExtendedVector::getSlack() const

Get slack component. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getBifParam "double
LOCA::Pitchfork::MooreSpence::ExtendedVector::getBifParam() const

Get Bifurcation parameter. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getXVec "Teuchos::RCP<
NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedVector::getXVec()

Returns the solution vector component of extended vector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getNullVec "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Pitchfork::MooreSpence::ExtendedVector::getNullVec()

Returns the null vector component of extended vector. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getSlack "double &
LOCA::Pitchfork::MooreSpence::ExtendedVector::getSlack()

Get slack component. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::ExtendedVector::getBifParam "double &
LOCA::Pitchfork::MooreSpence::ExtendedVector::getBifParam()

Get Bifurcation parameter. ";


// File: classLOCA_1_1AnasaziOperator_1_1Factory.xml
%feature("docstring") LOCA::AnasaziOperator::Factory "

Factory for creating Anasazi operator strategy objects.

The parameters passed to the create() through the eigenParams argument
method should specify the \"Operator\" as described below, as well as
any additional parameters for the particular strategy. \"Operator\" -
Name of the Anasazi operator. Valid choices are \"Jacobian Inverse\" (
LOCA::AnasaziOperator::JacobianInverse) [Default]

\"Shift-Invert\" ( LOCA::AnasaziOperator::ShiftInvert)

\"Cayley\" ( LOCA::AnasaziOperator::Cayley)

There is also an Epetra specific strategy that can be instantiated by
the LOCA::Epetra::Factory. See LOCA::Epetra::AnasaziOperator::Floquet.

C++ includes: LOCA_AnasaziOperator_Factory.H ";

%feature("docstring")  LOCA::AnasaziOperator::Factory::Factory "LOCA::AnasaziOperator::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::AnasaziOperator::Factory::~Factory "LOCA::AnasaziOperator::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::AnasaziOperator::Factory::create "Teuchos::RCP< LOCA::AnasaziOperator::AbstractStrategy >
LOCA::AnasaziOperator::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Abstract::Group > &grp)

Create Anasazi operator strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

eigenParams:  [in] Eigensolver parameters as described above

solverParams:  [in] Linear solver parameters

grp:  [in] Group representing Jacobian/mass matrices ";

%feature("docstring")  LOCA::AnasaziOperator::Factory::strategyName "const std::string &
LOCA::AnasaziOperator::Factory::strategyName(Teuchos::ParameterList
&eigenParams) const

Return strategy name given by eigenParams. ";


// File: classLOCA_1_1Epetra_1_1TransposeLinearSystem_1_1Factory.xml
%feature("docstring") LOCA::Epetra::TransposeLinearSystem::Factory "

Factory for creating transpose linear system strategy objects.

The parameters passed to the create() through the solverParams
argument method should specify the \"Transpose Solver Method\" as
described below, as well as any additional parameters for the
particular strategy. \"Transpose Solver Method\" - Name of the method.
Valid choices are \"Tranpose Preconditioner\"
(NOX::Epetra::TransposeLinearSystem::TransposePreconditioner)
[Default]

\"Explicit Transpose\"
(NOX::Epetra::TransposeLinearSystem::ExplicitTranspose)

\"Left Preconditioning\"
(NOX::Epetra::TransposeLinearSystem::LeftPreconditioning)

C++ includes: LOCA_Epetra_TransposeLinearSystem_Factory.H ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::Factory::Factory "LOCA::Epetra::TransposeLinearSystem::Factory::Factory(const
Teuchos::RCP< LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::Factory::~Factory "LOCA::Epetra::TransposeLinearSystem::Factory::~Factory()

Destructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::Factory::create "Teuchos::RCP<
LOCA::Epetra::TransposeLinearSystem::AbstractStrategy >
LOCA::Epetra::TransposeLinearSystem::Factory::create(const
Teuchos::RCP< Teuchos::ParameterList > &solverParams, const
Teuchos::RCP< NOX::Epetra::LinearSystem > &linsys)

Create transpose solver strategy.

Parameters:
-----------

solverParams:  [in] Solver parameters as described above

linsys:  [in] Linear system solver ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::Factory::strategyName "const
std::string &
LOCA::Epetra::TransposeLinearSystem::Factory::strategyName(Teuchos::ParameterList
&solverParams) const

Return strategy name given by solverParams. ";


// File: classLOCA_1_1StepSize_1_1Factory.xml
%feature("docstring") LOCA::StepSize::Factory "

Factory for creating step size control strategy objects.

The parameters passed to the create() through the stepsizeParams
argument method should specify the \"Method\" as described below, as
well as any additional parameters for the particular strategy.
\"Method\" - Name of the step size control method. Valid choices are
\"Constant\" ( LOCA::StepSize::Constant)

\"Adaptive\" LOCA::StepSize::Adaptive) [Default]

C++ includes: LOCA_StepSize_Factory.H ";

%feature("docstring")  LOCA::StepSize::Factory::Factory "LOCA::StepSize::Factory::Factory(const Teuchos::RCP< LOCA::GlobalData
> &global_data)

Constructor. ";

%feature("docstring")  LOCA::StepSize::Factory::~Factory "LOCA::StepSize::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::StepSize::Factory::create "Teuchos::RCP<
LOCA::StepSize::AbstractStrategy >
LOCA::StepSize::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepsizeParams)

Create step size control strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

stepsizeParams:  [in] Step size parameters as described above ";

%feature("docstring")  LOCA::StepSize::Factory::strategyName "const
std::string &
LOCA::StepSize::Factory::strategyName(Teuchos::ParameterList
&stepsizeParams) const

Return strategy name given by stepsizeParams. ";


// File: classLOCA_1_1Eigensolver_1_1Factory.xml
%feature("docstring") LOCA::Eigensolver::Factory "

Factory for creating Eigensolver strategy objects.

The parameters passed to the create() through the eigenParams argument
method should specify the \"Method\" as described below, as well as
any additional parameters for the particular strategy. \"Method\" -
Name of the eigensolver method. Valid choices are \"Default\" (
LOCA::Eigensolver::DefaultStrategy) [Default]

\"Anasazi\" ( LOCA::Eigensolver::AnasaziStrategy)

There is also a LAPACK specific strategy that can be instantiated by
the LOCA::LAPACK::Factory. See LOCA::Eigensolver::DGGEVStrategy.

C++ includes: LOCA_Eigensolver_Factory.H ";

%feature("docstring")  LOCA::Eigensolver::Factory::Factory "LOCA::Eigensolver::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::Eigensolver::Factory::~Factory "LOCA::Eigensolver::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::Eigensolver::Factory::create "Teuchos::RCP< LOCA::Eigensolver::AbstractStrategy >
LOCA::Eigensolver::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Create eigensolver strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

eigenParams:  [in] Eigensolver parameters as described above ";

%feature("docstring")  LOCA::Eigensolver::Factory::strategyName "const std::string &
LOCA::Eigensolver::Factory::strategyName(Teuchos::ParameterList
&eigenParams) const

Return strategy name given by eigenParams. ";


// File: classLOCA_1_1EigenvalueSort_1_1Factory.xml
%feature("docstring") LOCA::EigenvalueSort::Factory "

Factory for creating EigenvalueSort strategy objects.

The parameters passed to the create() through the eigenParams argument
method should specify the sorting method \"Sorting Method\" as
described below, as well as any additional parameters for the
particular strategy. \"Sorting Order\" - Name of the sorting method.
Valid choices are \"LM\" ( LOCA::EigenvalueSort::LargestMagnitude)
[Default]

\"LR\" ( LOCA::EigenvalueSort::LargestReal)

\"LI\" ( LOCA::EigenvalueSort::LargestImaginary)

\"SM\" ( LOCA::EigenvalueSort::SmallestMagnitude)

\"SR\" ( LOCA::EigenvalueSort::SmallestReal)

\"SI\" ( LOCA::EigenvalueSort::SmallestImaginary)

\"CA\" ( LOCA::EigenvalueSort::LargestRealInverseCayley)

C++ includes: LOCA_EigenvalueSort_Factory.H ";

%feature("docstring")  LOCA::EigenvalueSort::Factory::Factory "LOCA::EigenvalueSort::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::EigenvalueSort::Factory::~Factory "LOCA::EigenvalueSort::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::Factory::create "Teuchos::RCP< LOCA::EigenvalueSort::AbstractStrategy >
LOCA::EigenvalueSort::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Create sorting strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

eigenParams:  [in] Eigensolver parameters as described above ";

%feature("docstring")  LOCA::EigenvalueSort::Factory::strategyName "const std::string &
LOCA::EigenvalueSort::Factory::strategyName(Teuchos::ParameterList
&eigenParams) const

Return strategy name given by eigenParams. ";


// File: classLOCA_1_1Abstract_1_1Factory.xml
%feature("docstring") LOCA::Abstract::Factory "

Abstract interface for providing a user-defined factory

LOCA::Abstract::Factory provides an abstract interface for providing
user-defined factories to the LOCA::Factory. The LOCA::Factory
provides a mechanism for instantiating different strategies based on
parameter list choices. This class allows additional strategies to be
instantiated by the factory without modifying the factory itself. This
is done by deriving a user-defined factory from this interface,
implementing any of the create methods for the user-defined
strategies, and passing an instance of the derived factory to the
LOCA::Factory object. Any derived class must implement the  init()
method to set the global data object which the factory can then pass
to any instantiated strategies.

C++ includes: LOCA_Abstract_Factory.H ";

/*  Strategy create methods  */

%feature("docstring")
LOCA::Abstract::Factory::createPredictorStrategy "bool
LOCA::Abstract::Factory::createPredictorStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList >
&predictorParams, Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy
> &strategy)

Create predictor strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createContinuationStrategy "bool
LOCA::Abstract::Factory::createContinuationStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList >
&stepperParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, const Teuchos::RCP<
LOCA::MultiPredictor::AbstractStrategy > &pred, const std::vector< int
> &paramIDs, Teuchos::RCP< LOCA::MultiContinuation::AbstractStrategy >
&strategy)

Create continuation strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createBifurcationStrategy "bool
LOCA::Abstract::Factory::createBifurcationStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList >
&bifurcationParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &strategy)

Create bifurcation strategy. ";

%feature("docstring")  LOCA::Abstract::Factory::createStepSizeStrategy
"bool LOCA::Abstract::Factory::createStepSizeStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepsizeParams, Teuchos::RCP<
LOCA::StepSize::AbstractStrategy > &strategy)

Create step size strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createBorderedSolverStrategy "bool
LOCA::Abstract::Factory::createBorderedSolverStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, Teuchos::RCP<
LOCA::BorderedSolver::AbstractStrategy > &strategy)

Create bordered system solver strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createEigensolverStrategy "bool
LOCA::Abstract::Factory::createEigensolverStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList > &eigenParams,
Teuchos::RCP< LOCA::Eigensolver::AbstractStrategy > &strategy)

Create eigensolver strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createEigenvalueSortStrategy "bool
LOCA::Abstract::Factory::createEigenvalueSortStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, Teuchos::RCP<
LOCA::EigenvalueSort::AbstractStrategy > &strategy)

Create eigenvalue sorting strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createSaveEigenDataStrategy "bool
LOCA::Abstract::Factory::createSaveEigenDataStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList > &eigenParams,
Teuchos::RCP< LOCA::SaveEigenData::AbstractStrategy > &strategy)

Create strategy to save eigenvector/value data. ";

%feature("docstring")
LOCA::Abstract::Factory::createAnasaziOperatorStrategy "bool
LOCA::Abstract::Factory::createAnasaziOperatorStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Abstract::Group > &grp, Teuchos::RCP<
LOCA::AnasaziOperator::AbstractStrategy > &strategy)

Create Anasazi operator. ";

%feature("docstring")
LOCA::Abstract::Factory::createMooreSpenceTurningPointSolverStrategy "bool
LOCA::Abstract::Factory::createMooreSpenceTurningPointSolverStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, Teuchos::RCP<
LOCA::TurningPoint::MooreSpence::SolverStrategy > &strategy)

Create Moore-Spence turning point solver strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createMooreSpencePitchforkSolverStrategy "bool
LOCA::Abstract::Factory::createMooreSpencePitchforkSolverStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, Teuchos::RCP<
LOCA::Pitchfork::MooreSpence::SolverStrategy > &strategy)

Create Moore-Spence pitchfork solver strategy. ";

%feature("docstring")
LOCA::Abstract::Factory::createMooreSpenceHopfSolverStrategy "bool
LOCA::Abstract::Factory::createMooreSpenceHopfSolverStrategy(const
std::string &strategyName, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, Teuchos::RCP<
LOCA::Hopf::MooreSpence::SolverStrategy > &strategy)

Create Moore-Spence Hopf solver strategy. ";

%feature("docstring")  LOCA::Abstract::Factory::Factory "LOCA::Abstract::Factory::Factory()

Constructor. ";

%feature("docstring")  LOCA::Abstract::Factory::~Factory "virtual
LOCA::Abstract::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::Abstract::Factory::init "virtual void
LOCA::Abstract::Factory::init(const Teuchos::RCP< LOCA::GlobalData >
&global_data)=0

Initialize factory.

The LOCA::Factory will call this method to initialize the user
provided factory. The user-provided factory should perform any needed
initialization here that cannot occur at construction. ";


// File: classLOCA_1_1MultiContinuation_1_1Factory.xml
%feature("docstring") LOCA::MultiContinuation::Factory "

Factory for creating continuation strategy objects.

The parameters passed to the create() through the stepperParams
argument method should specify the \"Continuation Method\" as
described below, as well as any additional parameters for the
particular strategy. \"Continuation Method\" - Name of the
continuation method. Valid choices are \"Arc Length\" (
LOCA::MultiContinuation::ArcLengthGroup) [Default]

\"Natural\" ( LOCA::MultiContinuation::NaturalGroup)

C++ includes: LOCA_MultiContinuation_Factory.H ";

%feature("docstring")  LOCA::MultiContinuation::Factory::Factory "LOCA::MultiContinuation::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::MultiContinuation::Factory::~Factory "LOCA::MultiContinuation::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::MultiContinuation::Factory::create "Teuchos::RCP< LOCA::MultiContinuation::AbstractStrategy >
LOCA::MultiContinuation::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepperParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, const Teuchos::RCP<
LOCA::MultiPredictor::AbstractStrategy > &pred, const std::vector< int
> &paramIDs)

Create continuation strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

stepperParams:  [in] Stepper parameters as described above

grp:  [in] Underlying group

pred:  [in] Predictor strategy

paramIDs:  [in] Indicies of continuation parameters ";

%feature("docstring")  LOCA::MultiContinuation::Factory::strategyName
"const std::string &
LOCA::MultiContinuation::Factory::strategyName(Teuchos::ParameterList
&stepperParams) const

Return strategy name given by stepperParams. ";


// File: classLOCA_1_1Epetra_1_1Factory.xml
%feature("docstring") LOCA::Epetra::Factory "

Implementation of the LOCA::Abstract::Factory for Epetra groups.

C++ includes: LOCA_Epetra_Factory.H ";

/*  Strategy create methods  */

%feature("docstring")
LOCA::Epetra::Factory::createBorderedSolverStrategy "bool
LOCA::Epetra::Factory::createBorderedSolverStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList >
&solverParams, Teuchos::RCP< LOCA::BorderedSolver::AbstractStrategy >
&strategy)

Create bordered system solver strategy. ";

%feature("docstring")
LOCA::Epetra::Factory::createAnasaziOperatorStrategy "bool
LOCA::Epetra::Factory::createAnasaziOperatorStrategy(const std::string
&strategyName, const Teuchos::RCP< LOCA::Parameter::SublistParser >
&topParams, const Teuchos::RCP< Teuchos::ParameterList > &eigenParams,
const Teuchos::RCP< Teuchos::ParameterList > &solverParams, const
Teuchos::RCP< NOX::Abstract::Group > &grp, Teuchos::RCP<
LOCA::AnasaziOperator::AbstractStrategy > &strategy)

Create Anasazi operator strategy for Floquet option. ";

%feature("docstring")  LOCA::Epetra::Factory::Factory "LOCA::Epetra::Factory::Factory()

Constructor. ";

%feature("docstring")  LOCA::Epetra::Factory::~Factory "LOCA::Epetra::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::Epetra::Factory::init "void
LOCA::Epetra::Factory::init(const Teuchos::RCP< LOCA::GlobalData >
&global_data)

Initialize factory. ";


// File: classLOCA_1_1MultiPredictor_1_1Factory.xml
%feature("docstring") LOCA::MultiPredictor::Factory "

Factory for creating Predictor strategy objects.

The parameters passed to the create() through the predictorParams
argument method should specify the \"Method\" as described below, as
well as any additional parameters for the particular strategy.
\"Method\" - Name of the predictor method. Valid choices are
\"Constant\" ( LOCA::MultiPredictor::Constant)

\"Tangent\" ( LOCA::MultiPredictor::Tangent)

\"Secant\" ( LOCA::MultiPredictor::Secant) [Default]

\"Random\" ( LOCA::MultiPredictor::Random)

\"Restart\" ( LOCA::MultiPredictor::Restart)

C++ includes: LOCA_MultiPredictor_Factory.H ";

%feature("docstring")  LOCA::MultiPredictor::Factory::Factory "LOCA::MultiPredictor::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::MultiPredictor::Factory::~Factory "LOCA::MultiPredictor::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::Factory::create "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::MultiPredictor::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &predictorParams)

Create predictor strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

predictorParams:  [in] Predictor parameters as described above ";

%feature("docstring")  LOCA::MultiPredictor::Factory::strategyName "const std::string &
LOCA::MultiPredictor::Factory::strategyName(Teuchos::ParameterList
&predictorParams) const

Return strategy name given by predictorParams. ";


// File: classLOCA_1_1BorderedSolver_1_1Factory.xml
%feature("docstring") LOCA::BorderedSolver::Factory "

Factory for creating BorderedSolver strategy objects.

The parameters passed to the create() through the solverParams
argument method should specify the \"Bordered Solver Method\" as
described below, as well as any additional parameters for the
particular strategy. \"Bordered Solver Method\" - Name of the method.
Valid choices are \"Bordering\" ( LOCA::BorderedSolver::Bordering)
[Default]

\"Nested\" ( LOCA::BorderedSolver::Nested)

There are also Epetra and LAPACK specific strategies that can be
instantiated by the LOCA::Epetra::Factory and LOCA::LAPACK::Factory.
See LOCA::BorderedSolver::LAPACKDirectSolve,
LOCA::BorderedSolver::EpetraHouseholder and
LOCA::BorderedSolver::Epetra::Augmented.

C++ includes: LOCA_BorderedSolver_Factory.H ";

%feature("docstring")  LOCA::BorderedSolver::Factory::Factory "LOCA::BorderedSolver::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::BorderedSolver::Factory::~Factory "LOCA::BorderedSolver::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::BorderedSolver::Factory::create "Teuchos::RCP< LOCA::BorderedSolver::AbstractStrategy >
LOCA::BorderedSolver::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Create bordered system solver strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

solverParams:  [in] Solver parameters as described above ";

%feature("docstring")  LOCA::BorderedSolver::Factory::strategyName "const std::string &
LOCA::BorderedSolver::Factory::strategyName(Teuchos::ParameterList
&solverParams) const

Return strategy name given by solverParams. ";


// File: classLOCA_1_1Factory.xml
%feature("docstring") LOCA::Factory "

Factory class for creating strategies

The Factory class provides a single location for instantiating various
strategies based on parameter list choices. It provides a create()
method for each type of strategy which instantiates strategy objects
for that type. Each create method takes as arguments a ref-count
pointer to a LOCA::Parameter::SublistParser and a parameter list. The
parameter list determines which strategy to choose and also should
provide any parameters the strategy requires. The sublist parser
provides a parsed version of the top-level parameter list and allows
strategies to easily obtain other sublists from the top-level list. A
user-supplied factory may also be provided for instantiating user-
defined strategies. If a user-defined factory is supplied, each create
method will first attempt to instantiate the strategy using it, and
then instantiate strategies itself if necessary.

C++ includes: LOCA_Factory.H ";

/*  Strategy create methods  */

%feature("docstring")  LOCA::Factory::createPredictorStrategy "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::Factory::createPredictorStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &predictorParams)

Create predictor strategy.

Instantiates a predictor strategy based on the \"Method\" parameter of
the \"Predictor\" sublist. See LOCA::MultiPredictor::Factory for a
description of available strategies. ";

%feature("docstring")  LOCA::Factory::createContinuationStrategy "Teuchos::RCP< LOCA::MultiContinuation::AbstractStrategy >
LOCA::Factory::createContinuationStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepperParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, const Teuchos::RCP<
LOCA::MultiPredictor::AbstractStrategy > &pred, const std::vector< int
> &paramIDs)

Create continuation strategy.

Instantiates a continuation strategy based on the \"Continuation
Method\" parameter of the \"Stepper\" sublist. See
LOCA::MultiContinuation::Factory for a description of available
strategies. ";

%feature("docstring")  LOCA::Factory::createBifurcationStrategy "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Factory::createBifurcationStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &bifurcationParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp)

Create bifurcation strategy.

Instantiates a bifurcation strategy based on the \"Method\" parameter
of the \"Bifurcation\" sublist. See LOCA::Bifurcation::Factory for a
description of available strategies. ";

%feature("docstring")  LOCA::Factory::createStepSizeStrategy "Teuchos::RCP< LOCA::StepSize::AbstractStrategy >
LOCA::Factory::createStepSizeStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &stepsizeParams)

Create step size control strategy.

Instantiates a step size control strategy based on the \"Method\"
parameter of the \"Step Size\" sublist. See LOCA::StepSize::Factory
for a description of available strategies. ";

%feature("docstring")  LOCA::Factory::createBorderedSolverStrategy "Teuchos::RCP< LOCA::BorderedSolver::AbstractStrategy >
LOCA::Factory::createBorderedSolverStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Create bordered system solver strategy.

Instantiates an bordered system solver strategy based on the
\"Bordered Solver Method\" parameter of the \"Linear Solver\" sublist.
See LOCA::BorderedSolver::Factory for a description of available
strategies. ";

%feature("docstring")  LOCA::Factory::createEigensolverStrategy "Teuchos::RCP< LOCA::Eigensolver::AbstractStrategy >
LOCA::Factory::createEigensolverStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Create eigensolver strategy.

Instantiates an eigensolver strategy based on the \"Method\" parameter
of the \"Eigensolver\" sublist. See LOCA::Eigensolver::Factory for a
description of available strategies. ";

%feature("docstring")  LOCA::Factory::createEigenvalueSortStrategy "Teuchos::RCP< LOCA::EigenvalueSort::AbstractStrategy >
LOCA::Factory::createEigenvalueSortStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Create eigenvalue sort strategy.

Instantiates an eigenvalue sorting strategy based on the \"Sorting
Method\" parameter of the \"Eigensolver\" sublist. See
LOCA::EigenvalueSort::Factory for a description of available
strategies. ";

%feature("docstring")  LOCA::Factory::createSaveEigenDataStrategy "Teuchos::RCP< LOCA::SaveEigenData::AbstractStrategy >
LOCA::Factory::createSaveEigenDataStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Create strategy to save eigenvector/value data.

Instantiates a strategy to save eigenvector/value data based on the
\"Save Eigen Data Method\" parameter of the \"Eigensolver\" sublist.
See LOCA::SaveEigenData::Factory for a description of available
strategies. ";

%feature("docstring")  LOCA::Factory::createAnasaziOperatorStrategy "Teuchos::RCP< LOCA::AnasaziOperator::AbstractStrategy >
LOCA::Factory::createAnasaziOperatorStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Abstract::Group > &grp)

Create Anasazi operator strategy.

Instantiates an Anasazi operator strategy based on the \"Operator\"
parameter of the \"Eigensolver\" sublist. See
LOCA::AnasaziOperator::Factory for a description of available
strategies. ";

%feature("docstring")
LOCA::Factory::createMooreSpenceTurningPointSolverStrategy "Teuchos::RCP< LOCA::TurningPoint::MooreSpence::SolverStrategy >
LOCA::Factory::createMooreSpenceTurningPointSolverStrategy(const
Teuchos::RCP< LOCA::Parameter::SublistParser > &topParams, const
Teuchos::RCP< Teuchos::ParameterList > &solverParams)

Create Moore-Spence turning point solver strategy.

Instantiates a solver strategy based on the \"Solver Method\"
parameter of the \"Bifurcation\" sublist. See
LOCA::TurningPoint::MooreSpence::SolverFactory for a description of
available strategies. ";

%feature("docstring")
LOCA::Factory::createMooreSpencePitchforkSolverStrategy "Teuchos::RCP< LOCA::Pitchfork::MooreSpence::SolverStrategy >
LOCA::Factory::createMooreSpencePitchforkSolverStrategy(const
Teuchos::RCP< LOCA::Parameter::SublistParser > &topParams, const
Teuchos::RCP< Teuchos::ParameterList > &solverParams)

Create Moore-Spence pitchfork solver strategy.

Instantiates a solver strategy based on the \"Solver Method\"
parameter of the \"Bifurcation\" sublist. See
LOCA::Pitchfork::MooreSpence::SolverFactory for a description of
available strategies. ";

%feature("docstring")
LOCA::Factory::createMooreSpenceHopfSolverStrategy "Teuchos::RCP<
LOCA::Hopf::MooreSpence::SolverStrategy >
LOCA::Factory::createMooreSpenceHopfSolverStrategy(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Create Moore-Spence Hopf solver strategy.

Instantiates a solver strategy based on the \"Solver Method\"
parameter of the \"Bifurcation\" sublist. See
LOCA::Hopf::MooreSpence::SolverFactory for a description of available
strategies. ";

%feature("docstring")  LOCA::Factory::Factory "LOCA::Factory::Factory(const Teuchos::RCP< LOCA::GlobalData >
&global_data)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object. The constructor sets the
factory member of the global data to this. ";

%feature("docstring")  LOCA::Factory::Factory "LOCA::Factory::Factory(const Teuchos::RCP< LOCA::GlobalData >
&global_data, const Teuchos::RCP< LOCA::Abstract::Factory >
&userFactory)

Constructor with user-supplied factory.

Parameters:
-----------

global_data:  [in] Global data object. The constructor sets the
factory member of the global data to this.

userFactory:  [in] A user-supplied factory for instantiating user-
defined strategies. ";

%feature("docstring")  LOCA::Factory::~Factory "LOCA::Factory::~Factory()

Destructor. ";


// File: classLOCA_1_1Bifurcation_1_1Factory.xml
%feature("docstring") LOCA::Bifurcation::Factory "

Factory for creating bifurcation strategy objects.

The parameters passed to the create() through the bifurcationParams
argument method should specify the \"Type\" and \"Formulation\" as
described below, as well as any additional parameters for the
particular strategy. \"Type\" - Name of the bifurcation type. Valid
choices are \"None\" - No bifurcation [Default]

\"Turning Point\" Turning point bifurations

\"Pitchfork\" Pitchfork bifurcations

\"Hopf\" Hopf bifurcations

\"Phase Transition\" Phase Transition points

\"User-Defined\" User defined bifurcation. Set \"User-Defined Name\"
to be the parameter list name containing user-defined strategy, which
must be of type Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>.

\"Formulation\" - Name of the bifurcation formulation. Valid choices
are For turning point bifurcations: \"Moore-Spence\" [Default] (
LOCA::TurningPoint::MooreSpence::ExtendedGroup)

\"Minimally Augmented\" (
LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup)

For pitchfork bifurcations: \"Moore-Spence\" [Default] (
LOCA::Pitchfork::MooreSpence::ExtendedGroup)

\"Minimally Augmented\" (
LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup)

For Hopf bifurcations: \"Moore-Spence\" [Default] (
LOCA::Hopf::MooreSpence::ExtendedGroup)

\"Minimally Augmented\" (
LOCA::Hopf::MinimallyAugmented::ExtendedGroup)

C++ includes: LOCA_Bifurcation_Factory.H ";

%feature("docstring")  LOCA::Bifurcation::Factory::Factory "LOCA::Bifurcation::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::Bifurcation::Factory::~Factory "LOCA::Bifurcation::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::Bifurcation::Factory::create "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Bifurcation::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &bifurcationParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp)

Create bifurcation strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

bifurcationParams:  [in] Bifurcation parameters as described above

grp:  [in] Underlying group ";

%feature("docstring")  LOCA::Bifurcation::Factory::strategyName "std::string
LOCA::Bifurcation::Factory::strategyName(Teuchos::ParameterList
&bifurcationParams) const

Return strategy name given by bifurcationParams. ";


// File: classLOCA_1_1SaveEigenData_1_1Factory.xml
%feature("docstring") LOCA::SaveEigenData::Factory "

Factory for creating strategy objects to save eigenvectors/values.

The parameters passed to the create() through the eigenParams argument
method should specify the \"Method\" as described below, as well as
any additional parameters for the particular strategy. \"Save Eigen
Data Method\" - Name of the method. Valid choices are \"Default\" (
LOCA::SaveEigenData::DefaultStrategy) [Default]

\"User-Defined\" - User defined strategy

User-defined strategies are defined by supplying the parameter \"User-
Defined Save Eigen Data Name\" which is the std::string name of the
strategy, and then a parameter with this name that is of the type
Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>. This allows the
user to supply a user-defined strategy without providing a factory to
instantiate it. By supplying the name of the parameter storing the
strategy, the user can provide multiple strategies in the parameter
list and select among them by setting \"User-Defined Save Eigen Data
Name\" to be the name of the strategy.

C++ includes: LOCA_SaveEigenData_Factory.H ";

%feature("docstring")  LOCA::SaveEigenData::Factory::Factory "LOCA::SaveEigenData::Factory::Factory(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::SaveEigenData::Factory::~Factory "LOCA::SaveEigenData::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::SaveEigenData::Factory::create "Teuchos::RCP< LOCA::SaveEigenData::AbstractStrategy >
LOCA::SaveEigenData::Factory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Create strategy to save eigenvectors/values.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

eigenParams:  [in] Eigensolver parameters as described above ";

%feature("docstring")  LOCA::SaveEigenData::Factory::strategyName "const std::string &
LOCA::SaveEigenData::Factory::strategyName(Teuchos::ParameterList
&eigenParams) const

Return strategy name given by eigenParams. ";


// File: classLOCA_1_1StatusTest_1_1Factory.xml
%feature("docstring") LOCA::StatusTest::Factory "

Factory to build a set of status tests from a parameter list.

This object takes either an XML file name or a Teuchos::ParameterList
and generates an entire set (a tree) of status tests for use in a
LOCA::Stepper derived object.

The tagged_tests field in the constructors allows users to store tests
from the tree in a flat list in case they want to change the tolerance
values during a run. The tagged_tests flag is optional.

Please use the related nonmember functions instead of calling the
factory directly (See example below).

Valid parameters are as follows:

\"Test Type\" <std::string> Type of test this list contains. Valid
tests include: \"Combo\" - NOX::StatusTest::Combo

\"MaxIters\" - LOCA::StatusTest::MaxIters

\"User Defined\" - A user constructed test, derived from
NOX::StatusTest::Generic.

\"Tag\" <std::string> A unique identifier that will place the test in
the map for tagged_tests. This allows users to access individual tests
to change tolerances on the fly or query values while still using the
factory to build objects.

Additional parameters valid for a Combo test (
LOCA::StatusTest::Combo): \"Combo Type\" <std:string> Type of combo to
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

Additional parameters valid for a Maximum Iterations test
(NOX::StatusTest::MaxIters): \"Maximum Iterations\" <int>

Additional parameters valid for a \"User Defined\" test: \"User Status
Test\" < Teuchos::RCP<LOCA::StatusTest::Abstract> > A status test
suppied by the user. It is very important that when registering this
status test, that the user set it as a \"Generic\" object since there
is no implicit casting on the ParameterList's get method. See the
example below.

Example usage:

Nico Schloemer

C++ includes: LOCA_StatusTest_Factory.H ";

%feature("docstring")  LOCA::StatusTest::Factory::Factory "LOCA::StatusTest::Factory::Factory()

Constructor. ";

%feature("docstring")  LOCA::StatusTest::Factory::~Factory "LOCA::StatusTest::Factory::~Factory()

Destructor. ";

%feature("docstring")  LOCA::StatusTest::Factory::buildStatusTests "Teuchos::RCP< LOCA::StatusTest::Abstract >
LOCA::StatusTest::Factory::buildStatusTests(const std::string
&file_name, const Teuchos::RCP< const LOCA::GlobalData > &globalData,
std::map< std::string, Teuchos::RCP< LOCA::StatusTest::Abstract > >
*tagged_tests=0) const

Returns a status test set from a parameter list xml file. ";

%feature("docstring")  LOCA::StatusTest::Factory::buildStatusTests "Teuchos::RCP< LOCA::StatusTest::Abstract >
LOCA::StatusTest::Factory::buildStatusTests(Teuchos::ParameterList &p,
const Teuchos::RCP< const LOCA::GlobalData > &globalData, std::map<
std::string, Teuchos::RCP< LOCA::StatusTest::Abstract > >
*tagged_tests=0) const

Returns a status test set from a parameter list. ";


// File: classLOCA_1_1TurningPoint_1_1MinimallyAugmented_1_1FiniteDifferenceGroup.xml
%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup "

Concrete class that provides concrete implementations of the
derivative computation methods of the
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup using first-
order finite differencing.

The finite-differencing calculations are actually implemented by the
LOCA::DerivUtils class, and a custom DerivUtils object can be passed
through the constructor of this class. However, in the future the
calculations encapsulated in the DerivUtils class may be incorporated
directly into this class and other finite- differencing child classes.

C++ includes:
LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup()

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup(const
FiniteDifferenceGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::~FiniteDifferenceGroup
"LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::~FiniteDifferenceGroup()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::computeDwtJnDp
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::computeDwtJnDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &w, const
NOX::Abstract::Vector &nullVector,
NOX::Abstract::MultiVector::DenseMatrix &result, bool isValid)

Computes the derivative $\\\\partial w^TJn/\\\\partial p$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDwtJnDp() method of the passed
LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::computeDwtJDp
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::computeDwtJDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &w,
NOX::Abstract::MultiVector &result, bool isValid)

Computes the derivative $\\\\partial w^TJ/\\\\partial p$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDwtJDp() method of the passed
LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::computeDwtJnDx
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::FiniteDifferenceGroup::computeDwtJnDx(const
NOX::Abstract::Vector &w, const NOX::Abstract::Vector &nullVector,
NOX::Abstract::Vector &result)

Computes the derivative $\\\\frac{\\\\partial w^TJn}{\\\\partial x}$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDwtJnDx() method of the passed
LOCA::DerivUtils object. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1FiniteDifferenceGroup.xml
%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup "

Concrete class that provides concrete implementations of the
derivative computation methods of the
LOCA::TurningPoint::MooreSpence::AbstractGroup using first-order
finite differencing.

The finite-differencing calculations are actually implemented by the
LOCA::DerivUtils class, and a custom DerivUtils object can be passed
through the constructor of this class. However, in the future the
calculations encapsulated in the DerivUtils class may be incorporated
directly into this class and other finite- differencing child classes.

C++ includes: LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup()

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup(const
FiniteDifferenceGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup
"LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDpMulti
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDpMulti(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &nullVector,
NOX::Abstract::MultiVector &result, bool isValid)

Computes the derivative $\\\\partial Jn/\\\\partial p$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDJnDp() method of the passed LOCA::DerivUtils
object. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDxaMulti
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDxaMulti(const
NOX::Abstract::Vector &nullVector, const NOX::Abstract::MultiVector
&aVector, NOX::Abstract::MultiVector &result)

Computes the directional derivative $\\\\frac{\\\\partial
Jn}{\\\\partial x} a$ for the given direction $a$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDJnDxa() method of the passed
LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDxaMulti
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDJnDxaMulti(const
NOX::Abstract::Vector &nullVector, const NOX::Abstract::Vector
&JnVector, const NOX::Abstract::MultiVector &aVector,
NOX::Abstract::MultiVector &result)

Computes the directional derivative $\\\\frac{\\\\partial
Jn}{\\\\partial x} a$ for the given direction $a$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDJnDxa() method of the passed
LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDwtJnDxMulti
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup::computeDwtJnDxMulti(const
NOX::Abstract::MultiVector &w, const NOX::Abstract::Vector
&nullVector, NOX::Abstract::MultiVector &result)

Computes the derivative $\\\\frac{\\\\partial w^TJn}{\\\\partial x}$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDwtJnDx() method of the passed
LOCA::DerivUtils object. ";


// File: classLOCA_1_1MultiContinuation_1_1FiniteDifferenceGroup.xml
%feature("docstring") LOCA::MultiContinuation::FiniteDifferenceGroup "

Concrete class that provides a concrete implementation of the
computeDfDp() method of the LOCA::Continuation::AbstractGroup using
first-order finite differencing.

The finite-differencing calculations are actually implemented by the
LOCA::DerivUtils class, and a custom DerivUtils object can be set by
the setDerivUtils() method. However, in the future the calculations
encapsulated in the DerivUtils class may be incorporated directly into
this class and other finite- differencing child classes.

C++ includes: LOCA_MultiContinuation_FiniteDifferenceGroup.H ";

%feature("docstring")
LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup()

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup(const
FiniteDifferenceGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::FiniteDifferenceGroup::~FiniteDifferenceGroup
"LOCA::MultiContinuation::FiniteDifferenceGroup::~FiniteDifferenceGroup()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::FiniteDifferenceGroup::copy "void
LOCA::MultiContinuation::FiniteDifferenceGroup::copy(const
NOX::Abstract::Group &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::FiniteDifferenceGroup::setDerivUtils "void
LOCA::MultiContinuation::FiniteDifferenceGroup::setDerivUtils(const
Teuchos::RCP< LOCA::DerivUtils > &deriv)

Set the LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::MultiContinuation::FiniteDifferenceGroup::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::FiniteDifferenceGroup::computeDfDpMulti(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector &dfdp, bool
isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDfDp() method of the passed LOCA::DerivUtils
object. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1FiniteDifferenceGroup.xml
%feature("docstring") LOCA::Hopf::MooreSpence::FiniteDifferenceGroup "

Concrete class that provides concrete implementations of the
derivative computation methods of the
LOCA::Hopf::MooreSpence::AbstractGroup using first-order finite
differencing.

The finite-differencing calculations are actually implemented by the
LOCA::DerivUtils class, and a custom DerivUtils object can be passed
through the constructor of this class. However, in the future the
calculations encapsulated in the DerivUtils class may be incorporated
directly into this class and other finite- differencing child classes.

C++ includes: LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H ";

%feature("docstring")
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup()

Constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup(const
FiniteDifferenceGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup
"LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDp "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &yVector,
const NOX::Abstract::Vector &zVector, double w,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag, bool isValid)

Computes the derivative $\\\\frac{\\\\partial (J+i\\\\omega
B)(y+iz)}{\\\\partial p}$ where $p$ is the parameter indexed by
paramIDs.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDCeDp() method of the passed LOCA::DerivUtils
object. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa(const
NOX::Abstract::Vector &yVector, const NOX::Abstract::Vector &zVector,
double w, const NOX::Abstract::MultiVector &aVector,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag)

Computes the directional derivative $\\\\frac{\\\\partial
(J+i\\\\omega B)(y+iz)}{\\\\partial x} a$ for the given direction $a$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDCeDxa() method of the passed
LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa(const
NOX::Abstract::Vector &yVector, const NOX::Abstract::Vector &zVector,
double w, const NOX::Abstract::MultiVector &aVector, const
NOX::Abstract::Vector &Ce_real, const NOX::Abstract::Vector &Ce_imag,
NOX::Abstract::MultiVector &result_real, NOX::Abstract::MultiVector
&result_imag)

Computes the directional derivative $\\\\frac{\\\\partial
(J+i\\\\omega B)(y+iz)}{\\\\partial x} a$ for the given direction $a$.
The arguments Ce_real and Ce_imag hold the real and imaginary
components of $(J+i\\\\omega B)(y+iz)$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDCeDxa() method of the passed
LOCA::DerivUtils object. ";


// File: classLOCA_1_1Hopf_1_1MinimallyAugmented_1_1FiniteDifferenceGroup.xml
%feature("docstring")
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup "

Concrete class that provides concrete implementations of the
derivative computation methods of the
LOCA::Hopf::MinimallyAugmented::AbstractGroup using first-order finite
differencing.

The finite-differencing calculations are actually implemented by the
LOCA::DerivUtils class, and a custom DerivUtils object can be passed
through the constructor of this class. However, in the future the
calculations encapsulated in the DerivUtils class may be incorporated
directly into this class and other finite- differencing child classes.

C++ includes: LOCA_Hopf_MinimallyAugmented_FiniteDifferenceGroup.H ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup()

Constructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup
"LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::FiniteDifferenceGroup(const
FiniteDifferenceGroup &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::~FiniteDifferenceGroup
"LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::~FiniteDifferenceGroup()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::computeDwtCeDp
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::computeDwtCeDp(const
std::vector< int > &paramIDs, const NOX::Abstract::Vector &w1, const
NOX::Abstract::Vector &w2, const NOX::Abstract::Vector &y, const
NOX::Abstract::Vector &x, double omega,
NOX::Abstract::MultiVector::DenseMatrix &result_real,
NOX::Abstract::MultiVector::DenseMatrix &result_imag, bool isValid)

Computes the derivative $\\\\partial w^TCe/\\\\partial p$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDwtCeDp() method of the passed
LOCA::DerivUtils object. ";

%feature("docstring")
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::computeDwtCeDx
"NOX::Abstract::Group::ReturnType
LOCA::Hopf::MinimallyAugmented::FiniteDifferenceGroup::computeDwtCeDx(const
NOX::Abstract::Vector &w1, const NOX::Abstract::Vector &w2, const
NOX::Abstract::Vector &y, const NOX::Abstract::Vector &z, double
omega, NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag)

Computes the derivative $\\\\frac{\\\\partial w^TCe}{\\\\partial x}$.

The calculation is implemented by calling the corresponding
LOCA::DerivUtils::computeDwtCeDx() method of the passed
LOCA::DerivUtils object. ";


// File: classLOCA_1_1Epetra_1_1AnasaziOperator_1_1Floquet.xml
%feature("docstring") LOCA::Epetra::AnasaziOperator::Floquet "

Anasazi operator for computing generalized eigenvalues using Cayley
transformations.

This class implements the LOCA::AnasaziOperator::AbstractStrategy
interface for computing generalized eigenvalues $\\\\lambda$ and
eigenvectors $z$ of the system \\\\[ J z = \\\\lambda M z *\\\\] where
$J$ is the Jacobian matrix and $M$ is the mass matrix. The eigenvalues
are computed using a Cayley transformation, i.e. solving \\\\[ (J -
\\\\sigma M) z = (J - \\\\mu M) r \\\\] where $\\\\sigma$ is the
Cayley pole and $\\\\mu$ is the Cayley zero.

The parameters used by this class supplied in the constructor are:
\"Cayley Pole\" - $\\\\sigma$ as defined above (Default 0.0)

\"Cayley Zero\" - $\\\\mu$ as defined above (Default 0.0)

Also the grp argument to the constructor must be a child of
LOCA::TimeDependent::AbstractGroup for the shift-invert operations.

C++ includes: LOCA_Epetra_AnasaziOperator_Floquet.H ";

%feature("docstring")  LOCA::Epetra::AnasaziOperator::Floquet::Floquet
"LOCA::Epetra::AnasaziOperator::Floquet::Floquet(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Abstract::Group > &grp)

Constructor.

Argument grp must be of type LOCA::TimeDependent::AbstractGroup. See
class description for a list of eigenParams. ";

%feature("docstring")
LOCA::Epetra::AnasaziOperator::Floquet::~Floquet "LOCA::Epetra::AnasaziOperator::Floquet::~Floquet()

Destructor. ";

%feature("docstring")  LOCA::Epetra::AnasaziOperator::Floquet::label "const std::string & LOCA::Epetra::AnasaziOperator::Floquet::label()
const

Return name of this operator. ";

%feature("docstring")  LOCA::Epetra::AnasaziOperator::Floquet::apply "void LOCA::Epetra::AnasaziOperator::Floquet::apply(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &output)
const

Apply the operator.

Applies the inverse of the shifted operator, i.e., solves \\\\[
(J-\\\\omega I)z = M r \\\\] for $z$, where $r = \\\\mbox{input}$ and
$z = \\\\mbox{output}$. ";

%feature("docstring")
LOCA::Epetra::AnasaziOperator::Floquet::transformEigenvalue "void
LOCA::Epetra::AnasaziOperator::Floquet::transformEigenvalue(double
&ev_r, double &ev_i) const

Transform eigenvalue.

Transforms the given eigenvalue to the eigenvalue of the Jacobian-mass
matrix system by shifting and inverting it. ";

%feature("docstring")
LOCA::Epetra::AnasaziOperator::Floquet::rayleighQuotient "NOX::Abstract::Group::ReturnType
LOCA::Epetra::AnasaziOperator::Floquet::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const

Compute Rayleigh quotient.

Computes the Rayleigh quotient $z^T J z / z^T M z$ for the eigenvector
$z$. ";


// File: classLOCA_1_1Epetra_1_1Interface_1_1FreeEnergy.xml
%feature("docstring") LOCA::Epetra::Interface::FreeEnergy "

Used by LOCA::Epetra::Group to provide a link to the external code for
computing the free energy.

C++ includes: LOCA_Epetra_Interface_FreeEnergy.H ";

%feature("docstring")  LOCA::Epetra::Interface::FreeEnergy::FreeEnergy
"LOCA::Epetra::Interface::FreeEnergy::FreeEnergy()

Constructor. ";

%feature("docstring")
LOCA::Epetra::Interface::FreeEnergy::~FreeEnergy "virtual
LOCA::Epetra::Interface::FreeEnergy::~FreeEnergy()

Destructor. ";

%feature("docstring")
LOCA::Epetra::Interface::FreeEnergy::computeFreeEnergy "virtual
double LOCA::Epetra::Interface::FreeEnergy::computeFreeEnergy(const
Epetra_Vector &x)=0

Call user routine for computing the FreeEnergy of a system, for use in
Phase Transition Tracking alg. ";


// File: classLOCA_1_1SingularJacobianSolve_1_1Generic.xml
%feature("docstring") LOCA::SingularJacobianSolve::Generic "

Generic singular jacobian solve interface.

Generic interface for solving $Jx=b$ when $J$ is (nearly) singular.
All classes the implement a method for computing solutions to nearly
singular systems should be derived from this class.

C++ includes: LOCA_SingularJacobianSolve_Generic.H ";

%feature("docstring")  LOCA::SingularJacobianSolve::Generic::Generic "LOCA::SingularJacobianSolve::Generic::Generic()

Constructor.

Constructors of derived objects should look like reset. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Generic::Generic "LOCA::SingularJacobianSolve::Generic::Generic(const Generic &source)

Copy constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Generic::~Generic
"virtual LOCA::SingularJacobianSolve::Generic::~Generic()

Destructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Generic::clone "virtual Generic* LOCA::SingularJacobianSolve::Generic::clone() const
=0

Clone function. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Generic::reset "virtual NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Generic::reset(Teuchos::ParameterList
&params)=0

Reset parameters. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Generic::compute "virtual NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Generic::compute(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector &input, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector &result)=0

Compute solution to singular system.

The passed parameters are assumed be the (nonsingular) linear solver
parameters. ";

%feature("docstring")
LOCA::SingularJacobianSolve::Generic::computeMulti "virtual
NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Generic::computeMulti(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector *const *inputs, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector **results, int nVecs)=0

Compute solution to singular system with multiple RHS.

The passed parameters are assumed be the (nonsingular) linear solver
parameters. ";


// File: classLOCA_1_1GlobalData.xml
%feature("docstring") LOCA::GlobalData "

Container class to hold \"global\" LOCA objects.

GlobalData is a container class that holds ref-count pointers to
\"global\" objects, i.e., objects that nearly every LOCA object will
need access to. By putting them all in one container class, the
container class can be stored in each LOCA object, and if a new global
object is needed, it can be added here without modifying the rest of
the code. This is an alternative to true global or static objects
which are note safe in many contexts. In particular, this approach
allows multiple LOCA \"invocations\" to be in memory at the same time.
Note that all data members are declared public.

C++ includes: LOCA_GlobalData.H ";

%feature("docstring")  LOCA::GlobalData::GlobalData "LOCA::GlobalData::GlobalData(const Teuchos::RCP< NOX::Utils >
&loca_utils, const Teuchos::RCP< LOCA::ErrorCheck > &loca_error_check,
const Teuchos::RCP< LOCA::Factory > &loca_factory)

Constructor taking a ref-count pointer to each global object. ";

%feature("docstring")  LOCA::GlobalData::~GlobalData "LOCA::GlobalData::~GlobalData()

Destructor. ";


// File: classLOCA_1_1Abstract_1_1Group.xml
%feature("docstring") LOCA::Abstract::Group "

Compatiblity class for AbstractGroup hierarchy.

This class is derived from all LOCA AbstractGroup abstract base
classes as well as all FiniteDifference groups and any other groups
that provided default implementations for AbstractGroup pure virtual
methods. This class provides definitions for all needed assignment
operators and provides definitions for some pure virtual methods by
printing error messages. This class exists primarily for compatiblity
to an older class hierarchy and will most likely be removed in the
future.

C++ includes: LOCA_Abstract_Group.H ";

/*  Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
*/

%feature("docstring")
LOCA::Abstract::Group::augmentJacobianForHomotopy "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::augmentJacobianForHomotopy(double a, double b)

Replace Jacobian $J$ by $aJ+bI$ where $I$ is the identity matrix and
$p$ is a scalar.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

/*  Implementation of LOCA::TimeDependent::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Abstract::Group::computeShiftedMatrix "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeShiftedMatrix(double alpha, double beta)

Compute the shifted matrix.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")  LOCA::Abstract::Group::applyShiftedMatrix "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrix(const NOX::Abstract::Vector
&input, NOX::Abstract::Vector &result) const

Multiply the shifted matrix by a vector.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::Abstract::Group::applyShiftedMatrixMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrixMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Multiply the shifted matrix by a multi-vector.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::Abstract::Group::applyShiftedMatrixInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrixInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Apply the inverse of the shifted matrix by a multi-vector, as needed
by the shift-and-invert and generalized Cayley transformations.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::Abstract::Group::computeSecondShiftedMatrix "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeSecondShiftedMatrix(double alpha, double
beta)

Compute the second shifted matrix. Can avoid recomputing if two are
stored.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")  LOCA::Abstract::Group::applySecondShiftedMatrix
"NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applySecondShiftedMatrix(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Multiply the shifted matrix by a vector.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::Abstract::Group::applySecondShiftedMatrixMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applySecondShiftedMatrixMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Multiply the shifted matrix by a multi-vector.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

/*  Implementation of LOCA::Hopf::Moorespence::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Abstract::Group::isComplex "bool
LOCA::Abstract::Group::isComplex() const

Is $J+i\\\\omega B$ valid.

The implementation here always returns false. ";

%feature("docstring")  LOCA::Abstract::Group::computeComplex "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeComplex(double frequency)

Compute $J+i\\\\omega B$.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")  LOCA::Abstract::Group::applyComplex "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplex(const NOX::Abstract::Vector
&input_real, const NOX::Abstract::Vector &input_imag,
NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag) const

Compute $(J+i\\\\omega B)(y+iz)$.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")  LOCA::Abstract::Group::applyComplexMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexMultiVector(const
NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Compute $(J+i\\\\omega B)(y+iz)$.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::Abstract::Group::applyComplexInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Solve $(J+i\\\\omega B)(y+iz) = a+ib$.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

/*  Implementation of LOCA::Hopf::MinimallyAugmented::AbstractGroup
virtual methods.  */

%feature("docstring")  LOCA::Abstract::Group::applyComplexTranspose "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexTranspose(const
NOX::Abstract::Vector &input_real, const NOX::Abstract::Vector
&input_imag, NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag) const

Computes conjugate-tranpose matrix vector product $ (J+i\\\\omega B)^H
(x + iy) $.

Implementation here prints an error message and returns
NOX::Abstract::Group::NotDefined. ";

%feature("docstring")
LOCA::Abstract::Group::applyComplexTransposeMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexTransposeMultiVector(const
NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Computes conjugate-tranpose matrix vector product $ (J+i\\\\omega B)^H
(x + iy) $. ";

%feature("docstring")
LOCA::Abstract::Group::applyComplexTransposeInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Solve $(J+i\\\\omega B)^H (x + iy) = a+ib$. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Abstract::Group::copy "void
LOCA::Abstract::Group::copy(const NOX::Abstract::Group &source)

Assignment operator. ";

%feature("docstring")  LOCA::Abstract::Group::setParamsMulti "void
LOCA::Abstract::Group::setParamsMulti(const std::vector< int >
&paramIDs, const NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")  LOCA::Abstract::Group::notifyCompletedStep "void LOCA::Abstract::Group::notifyCompletedStep()

Notify group that the continuation step is completed The default
implementation here is to do nothing. ";

/*  Implementation of NOX::Abstract::Group virtual methods.  */

/*  Implementation of LOCA::PhaseTransition::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Abstract::Group::computeFreeEnergy "double LOCA::Abstract::Group::computeFreeEnergy()

Computes the free energy at the current solution and parameter values.
";

%feature("docstring")  LOCA::Abstract::Group::Group "LOCA::Abstract::Group::Group(const Teuchos::RCP< LOCA::GlobalData >
&global_data)

Constructor. ";

%feature("docstring")  LOCA::Abstract::Group::Group "LOCA::Abstract::Group::Group(const Teuchos::RCP< LOCA::GlobalData >
&global_data, const Teuchos::RCP< LOCA::DerivUtils > &deriv)

Constructor. ";

%feature("docstring")  LOCA::Abstract::Group::Group "LOCA::Abstract::Group::Group(const Group &source, NOX::CopyType
type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::Abstract::Group::~Group "LOCA::Abstract::Group::~Group()

Destructor. ";


// File: classLOCA_1_1Epetra_1_1Group.xml
%feature("docstring") LOCA::Epetra::Group "

Extension of the NOX::Epetra::Group to LOCA.

This class extends the NOX::Epetra::Group to LOCA enabling
continuation and bifurcation capabilities using Epetra. It is derived
from the NOX::Epetra::Group (basic Epetra support), the
LOCA::Abstract::Group (brings in all LOCA abstract base classes), and
the LOCA::Abstract::TransposeSolveGroup (for
applyJacobianTransposeInverse() methods). It stores a parameter vector
for setting/retrieving parameter values and overloads the computeF()
and computeJacobian() methods of the NOX::Epetra::Group parent class
to set the entire contents of the parameter vector in the problem
interface before calling the NOX::Epetra::Group computeF() and
computeJacobian().

Since it is derived from the LOCA::Abstract::Group (which is in turn
derived from all FiniteDifference groups), it uses the finite-
difference implementations for all parameter derivatives and second
derivatives. However this behavior can be modified by calling the
setDerivUtils() method of the
LOCA::MultiContinuation::FiniteDifferenceGroup parent class.

This class provides complete support for all continuation and
bifurcation methods including shift-invert and Cayley methods for
computing eigenvalues and Hopf bifurcations. However this support is
only enabled by calling the appropriate constructor described below.

C++ includes: LOCA_Epetra_Group.H ";

/*  Overloaded NOX::Epetra::Group  methods.  */

%feature("docstring")  LOCA::Epetra::Group::clone "Teuchos::RCP<
NOX::Abstract::Group > LOCA::Epetra::Group::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")  LOCA::Epetra::Group::computeF "NOX::Abstract::Group::ReturnType LOCA::Epetra::Group::computeF()

Overloaded computeF()

Calls LOCA::Epetra::Interface::setParams before evalulating F. ";

%feature("docstring")  LOCA::Epetra::Group::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeJacobian()

Overloaded computeJacobian()

Calls LOCA::Epetra::Interface::setParams before evalulating J. ";

/*  Implementation of LOCA::Abstract::TransposeSolveGroup methods.  */

%feature("docstring")
LOCA::Epetra::Group::applyJacobianTransposeInverse "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Solve Jacobian-tranpose system.

In addition to all regular linear solver parameters, this method
references the following additional parameters: \"Transpose Solver
Method\" -- [string] (default: \"Transpose Preconditioner\") Method
for preconditioning the transpose linear system (
LOCA::Epetra::TransposeLinearSystem::Factory). Available choices are:
\"Transpose Preconditioner\" -- Use the transpose of the
preconditioner for the original system.

\"Left Preconditioning\" -- Use the transpose of the preconditioner,
and apply using left preconditioning.

\"Explicit Transpose\" -- Form the transpose of the matrix and compute
the preconditioner. This method is available only if Trilinos is
configured with EpetraExt support (--enable-epetraext). ";

%feature("docstring")
LOCA::Epetra::Group::applyJacobianTransposeInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyJacobianTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Solve Jacobian-tranpose system with multiple right-hand sides.

In addition to all regular linear solver parameters, this method
references the following additional parameters: \"Transpose Solver
Method\" -- [string] (default: \"Transpose Preconditioner\") Method
for preconditioning the transpose linear system (
LOCA::Epetra::TransposeLinearSystem::Factory). Available choices are:
\"Transpose Preconditioner\" -- Use the transpose of the
preconditioner for the original system.

\"Left Preconditioning\" -- Use the transpose of the preconditioner,
and apply using left preconditioning.

\"Explicit Transpose\" -- Form the transpose of the matrix and compute
the preconditioner. This method is available only if Trilinos is
configured with EpetraExt support (--enable-epetraext). ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Epetra::Group::copy "void
LOCA::Epetra::Group::copy(const NOX::Abstract::Group &source)

Copy. ";

%feature("docstring")  LOCA::Epetra::Group::setParams "void
LOCA::Epetra::Group::setParams(const ParameterVector &p)

Set the parameters. ";

%feature("docstring")  LOCA::Epetra::Group::setParam "void
LOCA::Epetra::Group::setParam(int paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")  LOCA::Epetra::Group::setParam "void
LOCA::Epetra::Group::setParam(std::string paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")  LOCA::Epetra::Group::getParams "const
LOCA::ParameterVector & LOCA::Epetra::Group::getParams() const

Return a const reference to the ParameterVector owned by the group. ";

%feature("docstring")  LOCA::Epetra::Group::getParam "double
LOCA::Epetra::Group::getParam(int paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")  LOCA::Epetra::Group::getParam "double
LOCA::Epetra::Group::getParam(std::string paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")  LOCA::Epetra::Group::preProcessContinuationStep
"void
LOCA::Epetra::Group::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. The implementation here is to call the corresponding
method in the interface. ";

%feature("docstring")
LOCA::Epetra::Group::postProcessContinuationStep "void
LOCA::Epetra::Group::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. The
implementation here is to call the corresponding method in the
interface. ";

%feature("docstring")  LOCA::Epetra::Group::projectToDraw "void
LOCA::Epetra::Group::projectToDraw(const NOX::Abstract::Vector &x,
double *px) const

Projects solution to a few scalars for multiparameter continuation.

This method is called every time a solution is saved by the
multiparameter continuation code MF for later visualization and should
project the solution vector down to a few scalars. The array px will
be preallocated to the proper length given by
projectToDrawDimension().

The implementation here is to call the corresponding method in the
interface. ";

%feature("docstring")  LOCA::Epetra::Group::projectToDrawDimension "int LOCA::Epetra::Group::projectToDrawDimension() const

Returns the dimension of the project to draw array.

The implementation here is to call the corresponding method in the
interface. ";

%feature("docstring")  LOCA::Epetra::Group::computeScaledDotProduct "double LOCA::Epetra::Group::computeScaledDotProduct(const
NOX::Abstract::Vector &a, const NOX::Abstract::Vector &b) const

Compute a scaled dot product.

The implementation here uses the scaling vector $s$ if one is
supplied: \\\\[ \\\\sum_{i=1}^n a_i*b_i*s_i*s_i. \\\\] If the scaling
vector is not provided, the standard dot product is used. ";

%feature("docstring")  LOCA::Epetra::Group::printSolution "void
LOCA::Epetra::Group::printSolution(const double conParam) const

Call the user interface print() routine, solution vector. ";

%feature("docstring")  LOCA::Epetra::Group::printSolution "void
LOCA::Epetra::Group::printSolution(const NOX::Abstract::Vector &x,
const double conParam) const

Call the user interface print() routine, any vector. ";

%feature("docstring")  LOCA::Epetra::Group::scaleVector "void
LOCA::Epetra::Group::scaleVector(NOX::Abstract::Vector &x) const

Scales a vector using scaling vector.

The implementation here uses the scaling vector $s$ if one is
supplied: \\\\[ x_i = a_i*s_i. \\\\] If the scaling vector is not
provided, the vector is rescaled by the square root of its length. ";

/*  Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
*/

%feature("docstring")  LOCA::Epetra::Group::augmentJacobianForHomotopy
"NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::augmentJacobianForHomotopy(double a, double b)

Replace Jacobian $J$ by $aJ+bI$ where $I$ is the identity matrix. ";

/*  Implementation of LOCA::TimeDependent::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Epetra::Group::computeShiftedMatrix "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeShiftedMatrix(double alpha, double beta)

Compute the shifted matrix. ";

%feature("docstring")  LOCA::Epetra::Group::applyShiftedMatrix "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrix(const NOX::Abstract::Vector
&input, NOX::Abstract::Vector &result) const

Multiply the shifted matrix by a vector. ";

%feature("docstring")
LOCA::Epetra::Group::applyShiftedMatrixMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrixMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Multiply the shifted matrix by a multi-vector. ";

%feature("docstring")
LOCA::Epetra::Group::applyShiftedMatrixInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyShiftedMatrixInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Apply the inverse of the shifted matrix by a multi-vector, as needed
by the shift-and-invert and generalized Cayley transformations. ";

%feature("docstring")  LOCA::Epetra::Group::computeSecondShiftedMatrix
"NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeSecondShiftedMatrix(double alpha, double
beta)

Compute the second shifted matrix (uses different memory then Shifted
matrix) ";

%feature("docstring")  LOCA::Epetra::Group::applySecondShiftedMatrix "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applySecondShiftedMatrix(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Multiply the second shifted matrix by a vector. ";

%feature("docstring")
LOCA::Epetra::Group::applySecondShiftedMatrixMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applySecondShiftedMatrixMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Multiply the second shifted matrix by a multi-vector. ";

/*  Implementation of LOCA::Hopf::MooreSpence::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Epetra::Group::isComplex "bool
LOCA::Epetra::Group::isComplex() const

Is $J+i\\\\omega B$ valid. ";

%feature("docstring")  LOCA::Epetra::Group::computeComplex "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeComplex(double frequency)

Compute $J+i\\\\omega B$.

The argument frequency stores $\\\\omega$. ";

%feature("docstring")  LOCA::Epetra::Group::applyComplex "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplex(const NOX::Abstract::Vector
&input_real, const NOX::Abstract::Vector &input_imag,
NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag) const

Compute $(J+i\\\\omega B)(y+iz)$. ";

%feature("docstring")  LOCA::Epetra::Group::applyComplexMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexMultiVector(const
NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Compute $(J+i\\\\omega B)(y+iz)$. ";

%feature("docstring")
LOCA::Epetra::Group::applyComplexInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Solve $(J+i\\\\omega B)(y+iz) = a+ib$. ";

/*  Implementation of LOCA::Hopf::MinimallyAugmented::AbstractGroup
virtual methods.  */

%feature("docstring")  LOCA::Epetra::Group::applyComplexTranspose "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexTranspose(const NOX::Abstract::Vector
&input_real, const NOX::Abstract::Vector &input_imag,
NOX::Abstract::Vector &result_real, NOX::Abstract::Vector
&result_imag) const

Computes conjugate-tranpose matrix vector product $ (J+i\\\\omega B)^H
(x + iy) $. ";

%feature("docstring")
LOCA::Epetra::Group::applyComplexTransposeMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexTransposeMultiVector(const
NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Computes conjugate-tranpose matrix vector product $ (J+i\\\\omega B)^H
(x + iy) $. ";

%feature("docstring")
LOCA::Epetra::Group::applyComplexTransposeInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::applyComplexTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input_real, const
NOX::Abstract::MultiVector &input_imag, NOX::Abstract::MultiVector
&result_real, NOX::Abstract::MultiVector &result_imag) const

Solve $(J+i\\\\omega B)^H (x + iy) = a+ib$. ";

/*  Implementation of LOCA::PhseTransition::AbstractGroup virtual
methods.  */

%feature("docstring")  LOCA::Epetra::Group::computeFreeEnergy "double
LOCA::Epetra::Group::computeFreeEnergy()

Computes the free energy at the current solution and parameter values.
";

%feature("docstring")  LOCA::Epetra::Group::Group "LOCA::Epetra::Group::Group(const Teuchos::RCP< LOCA::GlobalData >
&global_data, Teuchos::ParameterList &printingParams, const
Teuchos::RCP< LOCA::Epetra::Interface::Required > &i,
NOX::Epetra::Vector &initialGuess, const LOCA::ParameterVector &p)

Constructor with NO linear system (VERY LIMITED).

WARNING: If this constructor is used, then methods that require a
Jacobian or preconditioning will not be available. You will be limited
to simple algorithms like nonlinear-CG with no preconditioning. ";

%feature("docstring")  LOCA::Epetra::Group::Group "LOCA::Epetra::Group::Group(const Teuchos::RCP< LOCA::GlobalData >
&global_data, Teuchos::ParameterList &printingParams, const
Teuchos::RCP< LOCA::Epetra::Interface::Required > &i,
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linSys, const LOCA::ParameterVector &p)

Standard Constructor enabling most LOCA support.

This is the most commonly used constructor and provides support for
all LOCA algorithms except shift-invert and Cayley transformations and
Hopf bifurcations. ";

%feature("docstring")  LOCA::Epetra::Group::Group "LOCA::Epetra::Group::Group(const Teuchos::RCP< LOCA::GlobalData >
&global_data, Teuchos::ParameterList &printingParams, const
Teuchos::RCP< LOCA::Epetra::Interface::TimeDependent > &i,
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linSys, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &shiftedLinSys, const
LOCA::ParameterVector &p)

Constructor with time-dependent interface and shifted linear system.

Use this constructor to enable shift-invert and Cayley transformations
or Hopf bifurcations. It requires another interface to compute the
shifted matrix $\\\\alpha J + \\\\beta M$ where $J$ is the Jacobian
matrix and $M$ is the mass matrix, and a linear system object to solve
this system. Setting linSys = shiftedLinSys is a valid option for
passing the shifted solver, but this will cause the shifted matrix to
overwrite the Jacobian possibly resulting in more matrix fills. See
declareSeparateMatrixMemory() method below to assert separate memory.
";

%feature("docstring")  LOCA::Epetra::Group::Group "LOCA::Epetra::Group::Group(const Teuchos::RCP< LOCA::GlobalData >
&global_data, Teuchos::ParameterList &printingParams, const
Teuchos::RCP< LOCA::Epetra::Interface::TimeDependentMatrixFree > &i,
NOX::Epetra::Vector &initialGuess, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linSys, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &shiftedLinSys, const
LOCA::ParameterVector &p)

Constructor with time-dependent matrix-free interface and shifted
linear system.

This constructor may also be used for shift-invert and Cayley
transformations, but should be only be used for a matrix-free method
for solving the shifted system. ";

%feature("docstring")  LOCA::Epetra::Group::Group "LOCA::Epetra::Group::Group(const Group &source, NOX::CopyType
type=NOX::DeepCopy)

Copy constructor. If type is DeepCopy, takes ownership of valid shared
Jacobian and shared preconditioning matrix. ";

%feature("docstring")  LOCA::Epetra::Group::~Group "LOCA::Epetra::Group::~Group()

Destructor. ";

%feature("docstring")  LOCA::Epetra::Group::setFreeEnergyInterface "void LOCA::Epetra::Group::setFreeEnergyInterface(const Teuchos::RCP<
LOCA::Epetra::Interface::FreeEnergy > &iFE)

Method to inject an interface for calucatiuong the free energy. ";

%feature("docstring")
LOCA::Epetra::Group::declareSeparateMatrixMemory "void
LOCA::Epetra::Group::declareSeparateMatrixMemory(bool
separateMem=true)

Method for calling code to guarantee to LOCA that separate matrix. ";

%feature("docstring")  LOCA::Epetra::Group::getUserInterface "Teuchos::RCP< NOX::Epetra::Interface::Required >
LOCA::Epetra::Group::getUserInterface()

Return the userInterface. ";

%feature("docstring")  LOCA::Epetra::Group::printSolution "void
LOCA::Epetra::Group::printSolution(const NOX::Epetra::Vector &x, const
double conParam) const

Call the user interface print() routine, any vector. ";

%feature("docstring")  LOCA::Epetra::Group::setScaleVector "void
LOCA::Epetra::Group::setScaleVector(const NOX::Abstract::Vector &s)

Sets the scale vector. ";

%feature("docstring")
LOCA::Epetra::Group::setJacobianOperatorForSolve "void
LOCA::Epetra::Group::setJacobianOperatorForSolve(const Teuchos::RCP<
const Epetra_Operator > &op) const

Sets the Jacobian operator. ";

%feature("docstring")  LOCA::Epetra::Group::getComplexLinearSystem "Teuchos::RCP< const NOX::Epetra::LinearSystem >
LOCA::Epetra::Group::getComplexLinearSystem() const

Return the Linear System. ";

%feature("docstring")  LOCA::Epetra::Group::getComplexLinearSystem "Teuchos::RCP< NOX::Epetra::LinearSystem >
LOCA::Epetra::Group::getComplexLinearSystem()

Return the Linear System. ";

%feature("docstring")  LOCA::Epetra::Group::getComplexMaps "void
LOCA::Epetra::Group::getComplexMaps(Teuchos::RCP< const
Epetra_BlockMap > &baseMap, Teuchos::RCP< const Epetra_BlockMap >
&globalMap) const ";


// File: classLOCA_1_1Homotopy_1_1Group.xml
%feature("docstring") LOCA::Homotopy::Group "

LOCA's Homotopy Algorithm.

The HomotopyGroup is a concrete implementation of the
LOCA::Continuation::AbstractGroup that modifies the set of nonlinear
equations to be solved to allow for Homotopy to be applied to the
system. This object should be used in conjunction with the
LOCA::Stepper object to drive the continuation. This algorithm solves
a system of nonlinear equations supplied by the user ( $ F(x) $)
through continuation. An artificial parameter $ \\\\lambda $ is used
to control the continuation. The idea is to solve a simple equation
starting at $ \\\\lambda $ = 0 and, using the solution from the
previous step, solve systems of equations that gets progressively
closer to the true system of interest ( at $ \\\\lambda $ = 1.0 we
recover the original equations $ F(x) $). By constraining the
definition of $ g(x, \\\\lambda) $ and using artificial parameter
contiuation, the continuation branch should be free of multiplicity
and bifurcation phenomena.

The modified system of equations, $ g(x, \\\\lambda) $, supplied by
the HomotopyGroup is defined as:

\\\\[ g(x, \\\\lambda) = \\\\lambda F(x) + (1.0 - \\\\lambda)(x - a)
\\\\]

where $x$ is the solution vector, $ \\\\lambda $ is an artificial
parameter, $ F(x) $ is the set of nonlinear equations the user
supplies, $ g(x) $ is the corresponding set of homotopy equations that
LOCA will solve, and $ a $ is a random vector.

This group requires the loca Stepper for continuation from $
\\\\lambda $ = 0.0 (a simple set of equations to solve) to $
\\\\lambda $ = 1.0 (the set of equations requested by the user, $ F(x)
$). The Homotopy::Group will generate the Stepper parameter sublist in
the parameter list that is passed in to the constructor. The user is
free to modify this list (it sets default values) before passing it
into the stepper object but should NOT change the starting and
stopping values for the continuation parameter.

References:

ALGORITHM 652 HOMPACK: A Suite of Codes for Globally Convergent
Homotopy Algorithms, Watson, L.T., Billups, S.C, and Morgan, A.P., ACM
Transactions on Mathematical Software, Vol. 13, No. 3, September 1987,
pp281-310.

C++ includes: LOCA_Homotopy_Group.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")  LOCA::Homotopy::Group::clone "Teuchos::RCP<
NOX::Abstract::Group > LOCA::Homotopy::Group::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")  LOCA::Homotopy::Group::setX "void
LOCA::Homotopy::Group::setX(const NOX::Abstract::Vector &y)

Set the solution vector, x, to y. ";

%feature("docstring")  LOCA::Homotopy::Group::computeX "void
LOCA::Homotopy::Group::computeX(const NOX::Abstract::Group &g, const
NOX::Abstract::Vector &d, double step)

Compute this.x = grp.x + step * d. ";

%feature("docstring")  LOCA::Homotopy::Group::computeF "NOX::Abstract::Group::ReturnType LOCA::Homotopy::Group::computeF()

Compute the homotopy residual $g$. ";

%feature("docstring")  LOCA::Homotopy::Group::computeJacobian "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeJacobian()

Compute the Jacobian derivative of the homotopy residual $g$. ";

%feature("docstring")  LOCA::Homotopy::Group::computeGradient "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeGradient()

Compute gradient of homotopy residual $g$. ";

%feature("docstring")  LOCA::Homotopy::Group::computeNewton "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeNewton(Teuchos::ParameterList &params)

Compute Newton direction using applyJacobianInverse. ";

%feature("docstring")  LOCA::Homotopy::Group::applyJacobian "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobian(const NOX::Abstract::Vector
&input, NOX::Abstract::Vector &result) const

Computes the homotopy Jacobian vector product. ";

%feature("docstring")  LOCA::Homotopy::Group::applyJacobianTranspose "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianTranspose(const
NOX::Abstract::Vector &input, NOX::Abstract::Vector &result) const

Computes the homotopy Jacobian-transpose vector product. ";

%feature("docstring")  LOCA::Homotopy::Group::applyJacobianInverse "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const

Applies the inverse of the homotopy Jacobian matrix. ";

%feature("docstring")  LOCA::Homotopy::Group::applyJacobianMultiVector
"NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian for homotopy system. ";

%feature("docstring")
LOCA::Homotopy::Group::applyJacobianTransposeMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianTransposeMultiVector(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &result)
const

Applies Jacobian-transpose for homotopy system. ";

%feature("docstring")
LOCA::Homotopy::Group::applyJacobianInverseMultiVector "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::applyJacobianInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const

Applies Jacobian inverse for homotopy system. ";

%feature("docstring")  LOCA::Homotopy::Group::isF "bool
LOCA::Homotopy::Group::isF() const

Return true if the homotopy residual $g$ is valid. ";

%feature("docstring")  LOCA::Homotopy::Group::isJacobian "bool
LOCA::Homotopy::Group::isJacobian() const

Return true if the homotopy Jacobian is valid. ";

%feature("docstring")  LOCA::Homotopy::Group::isGradient "bool
LOCA::Homotopy::Group::isGradient() const

Return true if the homotopy gradient is valid. ";

%feature("docstring")  LOCA::Homotopy::Group::isNewton "bool
LOCA::Homotopy::Group::isNewton() const

Return true if the homotopy Newton direction is valid. ";

%feature("docstring")  LOCA::Homotopy::Group::getX "const
NOX::Abstract::Vector & LOCA::Homotopy::Group::getX() const

Return homotopy solution vector $x$. ";

%feature("docstring")  LOCA::Homotopy::Group::getF "const
NOX::Abstract::Vector & LOCA::Homotopy::Group::getF() const

Return homotopy residual $g$. ";

%feature("docstring")  LOCA::Homotopy::Group::getNormF "double
LOCA::Homotopy::Group::getNormF() const

Return 2-norm of $g$. ";

%feature("docstring")  LOCA::Homotopy::Group::getGradient "const
NOX::Abstract::Vector & LOCA::Homotopy::Group::getGradient() const

Return homotopy gradient. ";

%feature("docstring")  LOCA::Homotopy::Group::getNewton "const
NOX::Abstract::Vector & LOCA::Homotopy::Group::getNewton() const

Return homotopy Newton direction. ";

%feature("docstring")  LOCA::Homotopy::Group::getXPtr "Teuchos::RCP<
const NOX::Abstract::Vector > LOCA::Homotopy::Group::getXPtr() const

Return homotopy solution vector $x$. ";

%feature("docstring")  LOCA::Homotopy::Group::getFPtr "Teuchos::RCP<
const NOX::Abstract::Vector > LOCA::Homotopy::Group::getFPtr() const

Return homotopy residual $g$. ";

%feature("docstring")  LOCA::Homotopy::Group::getGradientPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::Group::getGradientPtr() const

Return homotopy gradient. ";

%feature("docstring")  LOCA::Homotopy::Group::getNewtonPtr "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Homotopy::Group::getNewtonPtr() const

Return homotopy Newton direction. ";

/*  Implementation of LOCA::Extended::MultiAbstractGroup  */

/* virtual methods

*/

%feature("docstring")  LOCA::Homotopy::Group::getUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Homotopy::Group::getUnderlyingGroup() const

Return underlying group. ";

%feature("docstring")  LOCA::Homotopy::Group::getUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Homotopy::Group::getUnderlyingGroup()

Return underlying group. ";

/*  Implementation of LOCA::MultiContinuation::AbstractGroup  */

/* virtual methods

*/

%feature("docstring")  LOCA::Homotopy::Group::copy "void
LOCA::Homotopy::Group::copy(const NOX::Abstract::Group &source)

Assignment. ";

%feature("docstring")  LOCA::Homotopy::Group::setParamsMulti "void
LOCA::Homotopy::Group::setParamsMulti(const std::vector< int >
&paramIDs, const NOX::Abstract::MultiVector::DenseMatrix &vals)

Set parameters indexed by (integer) paramIDs. ";

%feature("docstring")  LOCA::Homotopy::Group::setParams "void
LOCA::Homotopy::Group::setParams(const ParameterVector &p)

Set the parameter vector in the group to p. ";

%feature("docstring")  LOCA::Homotopy::Group::setParam "void
LOCA::Homotopy::Group::setParam(int paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")  LOCA::Homotopy::Group::setParam "void
LOCA::Homotopy::Group::setParam(std::string paramID, double val)

Set parameter indexed by paramID. ";

%feature("docstring")  LOCA::Homotopy::Group::getParams "const
LOCA::ParameterVector & LOCA::Homotopy::Group::getParams() const

Return a const reference to the paramter vector owned by the group. ";

%feature("docstring")  LOCA::Homotopy::Group::getParam "double
LOCA::Homotopy::Group::getParam(int paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")  LOCA::Homotopy::Group::getParam "double
LOCA::Homotopy::Group::getParam(std::string paramID) const

Return copy of parameter indexed by paramID. ";

%feature("docstring")  LOCA::Homotopy::Group::computeDfDpMulti "NOX::Abstract::Group::ReturnType
LOCA::Homotopy::Group::computeDfDpMulti(const std::vector< int >
&paramIDs, NOX::Abstract::MultiVector &dfdp, bool isValidF)

Compute $\\\\partial F/\\\\partial p$ for each parameter $p$ indexed
by paramIDs. The first column of dfdp holds F, which is valid if
isValidF is true. Otherwise F must be computed. ";

%feature("docstring")
LOCA::Homotopy::Group::preProcessContinuationStep "void
LOCA::Homotopy::Group::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. ";

%feature("docstring")
LOCA::Homotopy::Group::postProcessContinuationStep "void
LOCA::Homotopy::Group::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";

%feature("docstring")  LOCA::Homotopy::Group::projectToDraw "void
LOCA::Homotopy::Group::projectToDraw(const NOX::Abstract::Vector &x,
double *px) const

Projects solution to a few scalars for multiparameter continuation. ";

%feature("docstring")  LOCA::Homotopy::Group::projectToDrawDimension "int LOCA::Homotopy::Group::projectToDrawDimension() const

Returns the dimension of the project to draw array. ";

%feature("docstring")  LOCA::Homotopy::Group::printSolution "void
LOCA::Homotopy::Group::printSolution(const double conParam) const

Function to print out solution and continuation parameter after
successful continuation step. ";

%feature("docstring")  LOCA::Homotopy::Group::printSolution "void
LOCA::Homotopy::Group::printSolution(const NOX::Abstract::Vector &x_,
const double conParam) const

Function to print out solution and continuation parameter after
successful continuation step. ";

%feature("docstring")  LOCA::Homotopy::Group::Group "LOCA::Homotopy::Group::Group(Teuchos::ParameterList &locaSublist,
const Teuchos::RCP< LOCA::GlobalData > &global_data, const
Teuchos::RCP< LOCA::Homotopy::AbstractGroup > &g, double
scaleRandom=1.0, double scaleInitialGuess=0.0)

Constructor to set the base group and generate the \"%Stepper\"
sublist for homotopy continuation.

The locaSublist variable is the \"LOCA\" sublist (of type
Teuchos::ParameterList) that will be used in loca continuation runs.

The variables scalarRandomVector and scalarInitialGuess are used to
give some control over the generation of the random vector. In certain
instances we have seen the random vector force the solution to a set
of variables that are unphysical and could break the function
evaluations (cause them to return nan). For example, in heat transfer
problems, the temperature could be the dependent variable. If the
solution vector has an unphysical temperature ( the random vector
could force the temperature to negative or near zero values for the
solution at $ \\\\lambda = 0$) then property evaluations could break.
The random vector can be modified to keep the values near the initial
guess based on values supplied to the constructor of the
HomotopyGroup:

\\\\[ a = abs(r) * \\\\mbox{scalarRandom} + x_o *
\\\\mbox{scalarInitialGuess} \\\\]

where $ r $ is the random vector generated by a call to
NOX::Abstract::Vector::random(), $ \\\\mbox{scalarRandom} $ is a
scalar value, $ x_o $ is the initial guess to the solution vector, and
$ \\\\mbox{scalarInitialGuess} $ is a scalar value. The defualt values
force the random vector to be calculated as:

\\\\[ a = abs(r) \\\\]

IMPORTANT: For homotopy to work correctly you should not change the
starting and stopping parameter values (0.0 and 1.0 respectively) set
in the \"%Stepper\" sublist. ";

%feature("docstring")  LOCA::Homotopy::Group::Group "LOCA::Homotopy::Group::Group(Teuchos::ParameterList &locaSublist,
const Teuchos::RCP< LOCA::GlobalData > &global_data, const
Teuchos::RCP< LOCA::Homotopy::AbstractGroup > &g, const
NOX::Abstract::Vector &randomVector)

Constructor with a user supplied random vector. ";

%feature("docstring")  LOCA::Homotopy::Group::Group "LOCA::Homotopy::Group::Group(const Group &source, NOX::CopyType
type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::Homotopy::Group::~Group "LOCA::Homotopy::Group::~Group()

Destructor. ";


// File: classLOCA_1_1BorderedSolver_1_1HouseholderQR.xml
%feature("docstring") LOCA::BorderedSolver::HouseholderQR "

A convenience class to compute the QR factorization of a an extended
multi-vector.

This class computes the QR factorization \\\\[ Q^T \\\\begin{bmatrix}
op(C) \\\\\\\\ B \\\\end{bmatrix} = \\\\begin{bmatrix} R \\\\\\\\ 0
\\\\end{bmatrix} \\\\] where $C$ is an $m\\\\times m$ matrix, $B$ is
an $n\\\\times m$ matrix, $Q$ is an $n+m\\\\times n+m$ matrix, $R$ is
an $m\\\\times m$ matrix, and $op()$ represents either the identity
operation or the transpose. The matrix $C$ is represented by a
NOX::Abstract::MultiVector::DenseMatrix while $B$ is a
NOX::Abstract::MultiVector. Given $B$ and $C$, this class computes $Q$
and $R$ with $R$ returned as NOX::Abstract::MultiVector::DenseMatrix.
The operator $Q$ is generated using the standard Householder QR
algorithm (Algorithm 5.2.1, G. Golub and C. Van Loan, \"Matrix
Computations,\" 3rd Edition, Johns Hopkins, Baltimore, 1996) and is
stored using the compact WY representation: $Q = I + Y^T T Y$ (see R.
Schreiver and C. Van Loan, \"A Storage-Efficient WY Represntation  for
Products of Householder Transformations,\" SIAM J. Sci. Stat. Comput.,
Vol. 10, No. 1, pp. 53-57, January 1989).

C++ includes: LOCA_BorderedSolver_HouseholderQR.H ";

%feature("docstring")
LOCA::BorderedSolver::HouseholderQR::HouseholderQR "LOCA::BorderedSolver::HouseholderQR::HouseholderQR(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object ";

%feature("docstring")
LOCA::BorderedSolver::HouseholderQR::~HouseholderQR "LOCA::BorderedSolver::HouseholderQR::~HouseholderQR()

Destructor. ";

%feature("docstring")  LOCA::BorderedSolver::HouseholderQR::computeQR
"void LOCA::BorderedSolver::HouseholderQR::computeQR(const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector &B, bool use_c_transpose,
NOX::Abstract::MultiVector::DenseMatrix &Y1,
NOX::Abstract::MultiVector &Y2,
NOX::Abstract::MultiVector::DenseMatrix &T,
NOX::Abstract::MultiVector::DenseMatrix &R)

Compute QR factorization as described above.

Set use_c_transpose to true if the transpose of $C$ is required. ";

%feature("docstring")
LOCA::BorderedSolver::HouseholderQR::applyCompactWY "void
LOCA::BorderedSolver::HouseholderQR::applyCompactWY(const
NOX::Abstract::MultiVector::DenseMatrix &Y1, const
NOX::Abstract::MultiVector &Y2, const
NOX::Abstract::MultiVector::DenseMatrix &T,
NOX::Abstract::MultiVector::DenseMatrix &X1,
NOX::Abstract::MultiVector &X2, bool isZeroX1, bool isZeroX2, bool
useTranspose) const

Applies the operator Q as described above overwriting x and y. If
either of x or y are zero on input, set the corresponding isZeroX or
isZeroY flags. Set\\\\ useTranspose to true to instead apply the
transpose of Q. ";

%feature("docstring")
LOCA::BorderedSolver::HouseholderQR::applyCompactWY "void
LOCA::BorderedSolver::HouseholderQR::applyCompactWY(const
NOX::Abstract::MultiVector::DenseMatrix &Y1, const
NOX::Abstract::MultiVector &Y2, const
NOX::Abstract::MultiVector::DenseMatrix &T, const
NOX::Abstract::MultiVector::DenseMatrix *input1, const
NOX::Abstract::MultiVector *input2,
NOX::Abstract::MultiVector::DenseMatrix &result1,
NOX::Abstract::MultiVector &result2, bool useTranspose) const

Another version of applyCompactWY() that does not overwrite its
inputs. If either input is zero, set the corresponding pointer to
NULL. ";


// File: classLOCA_1_1Epetra_1_1IdentityOp.xml
%feature("docstring") LOCA::Epetra::IdentityOp "

An Epetra operator representing the identity matrix.

C++ includes: LOCA_Epetra_IdentityOp.H ";

%feature("docstring")  LOCA::Epetra::IdentityOp::IdentityOp "LOCA::Epetra::IdentityOp::IdentityOp(const Teuchos::RCP< const
Epetra_Comm > &comm, const Teuchos::RCP< const Epetra_Map > &map)

Constructor.

Parameters:
-----------

comm:  [in] Comm object

map:  [in] Map object ";

%feature("docstring")  LOCA::Epetra::IdentityOp::~IdentityOp "LOCA::Epetra::IdentityOp::~IdentityOp()

Destructor. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::SetUseTranspose "int
LOCA::Epetra::IdentityOp::SetUseTranspose(bool UseTranspose)

Set to true if the transpose of the operator is requested. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::Apply "int
LOCA::Epetra::IdentityOp::Apply(const Epetra_MultiVector &Input,
Epetra_MultiVector &Result) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector Input in Result as described above. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::ApplyInverse "int
LOCA::Epetra::IdentityOp::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of the inverse of the operator applied to a
Epetra_MultiVector Input in Result as described above. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::NormInf "double
LOCA::Epetra::IdentityOp::NormInf() const

Returns an approximate infinity norm of the operator matrix. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::Label "const char *
LOCA::Epetra::IdentityOp::Label() const

Returns a character std::string describing the operator. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::UseTranspose "bool
LOCA::Epetra::IdentityOp::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::HasNormInf "bool
LOCA::Epetra::IdentityOp::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::Comm "const
Epetra_Comm & LOCA::Epetra::IdentityOp::Comm() const

Returns a reference to the Epetra_Comm communicator associated with
this operator. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::OperatorDomainMap "const Epetra_Map & LOCA::Epetra::IdentityOp::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")  LOCA::Epetra::IdentityOp::OperatorRangeMap "const Epetra_Map & LOCA::Epetra::IdentityOp::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";


// File: classLOCA_1_1Abstract_1_1Iterator.xml
%feature("docstring") LOCA::Abstract::Iterator "

An abstract interface for implementing iteration.

The LOCA::Abstract::Iterator defines an interface for implementing
many kinds of iterative processes. In LOCA, this is used to implement
the Stepper which computes points along a continuation curve.

Many iterative processes can be abstracted in the following manner:

Initialize iteration (start)

Compute iteration (iterate) while iterator is not finished preprocess
step (preprocess)

compute step (compute)

postprocess step (posprocess)

check iterator status (stop)

Finalize iteration (finish)

The run method of the iterator implements this iterative process with
start, finish, preprocess, compute and postprocess left as pure
virtual methods to be implemented for the specific iterative process.

The iterator has one parameter, \"Max Steps\" (default 100) giving the
maximum number of steps the iterator should take. The default
implementation of stop only stops the iterator when this maximum
number of steps has been reached.

C++ includes: LOCA_Abstract_Iterator.H ";

%feature("docstring")  LOCA::Abstract::Iterator::Iterator "LOCA::Abstract::Iterator::Iterator(Teuchos::ParameterList &p)

Constructor. ";

%feature("docstring")  LOCA::Abstract::Iterator::Iterator "LOCA::Abstract::Iterator::Iterator(const Iterator &it)

Copy Constructor. ";

%feature("docstring")  LOCA::Abstract::Iterator::~Iterator "LOCA::Abstract::Iterator::~Iterator()

Destructor. ";

%feature("docstring")  LOCA::Abstract::Iterator::resetIterator "bool
LOCA::Abstract::Iterator::resetIterator(Teuchos::ParameterList &p)

Reset the iterator to start a new iteration. ";

%feature("docstring")  LOCA::Abstract::Iterator::getIteratorStatus "LOCA::Abstract::Iterator::IteratorStatus
LOCA::Abstract::Iterator::getIteratorStatus() const

Return the status of the iterator. ";

%feature("docstring")  LOCA::Abstract::Iterator::getStepNumber "int
LOCA::Abstract::Iterator::getStepNumber() const

Returns the number of accepted steps. ";

%feature("docstring")  LOCA::Abstract::Iterator::getNumFailedSteps "int LOCA::Abstract::Iterator::getNumFailedSteps() const

Returns the number of failed steps. ";

%feature("docstring")  LOCA::Abstract::Iterator::getNumTotalSteps "int LOCA::Abstract::Iterator::getNumTotalSteps() const

Returns the total number of steps attempted. ";

%feature("docstring")  LOCA::Abstract::Iterator::run "LOCA::Abstract::Iterator::IteratorStatus
LOCA::Abstract::Iterator::run()

Run the iterator. ";


// File: classLOCA_1_1SingularJacobianSolve_1_1ItRef.xml
%feature("docstring") LOCA::SingularJacobianSolve::ItRef "

This class computes the solution to $J x = b$ using one step of
iterative refinement.

This singular solve method uses one step of iterative refinement to
improve the accuracy of the solution to the linear system $J x = b$.
In particular, the algorithm used here is \\\\[ \\\\begin{aligned}
&\\\\text{Solve}\\\\; Jx_1 = b \\\\\\\\ &r = b - Jx_1 \\\\\\\\
&\\\\text{Solve}\\\\; Jx_2 = r \\\\\\\\ &x = x_1 + x_2
\\\\end{aligned} \\\\] Both solves use the underlying group's
applyJacobianInverse() method and therefore this is a generic
technique for computing solutions to nearly singular system since it
uses any supplied linear solver.

This algorithm is selected by setting the \"Method\" parameter of the
\"Singular Solve\" sublist of the NOX linear solver parameter list to
\"Iterative Refinement\".

C++ includes: LOCA_SingularJacobianSolve_ItRef.H ";

%feature("docstring")  LOCA::SingularJacobianSolve::ItRef::ItRef "LOCA::SingularJacobianSolve::ItRef::ItRef(Teuchos::ParameterList
&params)

Constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::ItRef::ItRef "LOCA::SingularJacobianSolve::ItRef::ItRef(const ItRef &source)

Copy constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::ItRef::~ItRef "LOCA::SingularJacobianSolve::ItRef::~ItRef()

Destructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::ItRef::clone "LOCA::SingularJacobianSolve::Generic *
LOCA::SingularJacobianSolve::ItRef::clone() const

Clone function. ";

%feature("docstring")  LOCA::SingularJacobianSolve::ItRef::reset "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::ItRef::reset(Teuchos::ParameterList
&params)

Reset parameters.

There are no additional parameters for the Nic calculation. ";

%feature("docstring")  LOCA::SingularJacobianSolve::ItRef::compute "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::ItRef::compute(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector &input, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector &result)

Computes the solution as described above. ";

%feature("docstring")
LOCA::SingularJacobianSolve::ItRef::computeMulti "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::ItRef::computeMulti(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector *const *inputs, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector **results, int nVecs)

Computes solution for multiple RHS. ";


// File: classLOCA_1_1AnasaziOperator_1_1JacobianInverse.xml
%feature("docstring") LOCA::AnasaziOperator::JacobianInverse "

Anasazi operator for computing eigenvalues of the inverse-Jacobian.

This class implements the LOCA::AnasaziOperator::AbstractStrategy
interface for computing eigenvalues of the inverse-Jacobian.

C++ includes: LOCA_AnasaziOperator_JacobianInverse.H ";

%feature("docstring")
LOCA::AnasaziOperator::JacobianInverse::JacobianInverse "LOCA::AnasaziOperator::JacobianInverse::JacobianInverse(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Abstract::Group > &grp)

Constructor. ";

%feature("docstring")
LOCA::AnasaziOperator::JacobianInverse::~JacobianInverse "LOCA::AnasaziOperator::JacobianInverse::~JacobianInverse()

Destructor. ";

%feature("docstring")  LOCA::AnasaziOperator::JacobianInverse::label "const std::string & LOCA::AnasaziOperator::JacobianInverse::label()
const

Return name of this operator. ";

%feature("docstring")  LOCA::AnasaziOperator::JacobianInverse::apply "void LOCA::AnasaziOperator::JacobianInverse::apply(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &output)
const

Apply the operator.

Computes $\\\\mbox{output} = J^{-1}\\\\mbox{input}$. ";

%feature("docstring")
LOCA::AnasaziOperator::JacobianInverse::beginPostProcessing "void
LOCA::AnasaziOperator::JacobianInverse::beginPostProcessing()

Begin PostProcessing of eigenvalues.

Compute Jacobian matrix once, for use in subsequent repeated calls to
rayleighQuotient ";

%feature("docstring")
LOCA::AnasaziOperator::JacobianInverse::transformEigenvalue "void
LOCA::AnasaziOperator::JacobianInverse::transformEigenvalue(double
&ev_r, double &ev_i) const

Transform eigenvalue.

Transforms the given eigenvalue to the eigenvalue of the Jacobian by
inverting it. ";

%feature("docstring")
LOCA::AnasaziOperator::JacobianInverse::rayleighQuotient "NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::JacobianInverse::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const

Compute Rayleigh quotient.

Computes the Rayleigh quotient $z^T J z$ for the eigenvector $z$. ";


// File: classLOCA_1_1BorderedSolver_1_1JacobianOperator.xml
%feature("docstring") LOCA::BorderedSolver::JacobianOperator "

Bordered solver operator representing the Jacobian as implemented in
the NOX::Abstract::Group.

C++ includes: LOCA_BorderedSolver_JacobianOperator.H ";

%feature("docstring")
LOCA::BorderedSolver::JacobianOperator::JacobianOperator "LOCA::BorderedSolver::JacobianOperator::JacobianOperator(const
Teuchos::RCP< const NOX::Abstract::Group > &grp)

Constructor. ";

%feature("docstring")
LOCA::BorderedSolver::JacobianOperator::~JacobianOperator "LOCA::BorderedSolver::JacobianOperator::~JacobianOperator()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::JacobianOperator::getGroup "Teuchos::RCP< const
NOX::Abstract::Group >
LOCA::BorderedSolver::JacobianOperator::getGroup() const

Get group pointer. ";

%feature("docstring")  LOCA::BorderedSolver::JacobianOperator::apply "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::apply(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const

Apply the operator. ";

%feature("docstring")
LOCA::BorderedSolver::JacobianOperator::applyTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::applyTranspose(const
NOX::Abstract::MultiVector &X, NOX::Abstract::MultiVector &Y) const

Apply transpose of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::JacobianOperator::applyInverse "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const

Apply inverse of the operator. ";

%feature("docstring")
LOCA::BorderedSolver::JacobianOperator::applyInverseTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::JacobianOperator::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &B,
NOX::Abstract::MultiVector &X) const

Apply inverse transpose of the operator.

Group must be of type LOCA::Abstract::TransposeSolveGroup for this
method to be defined. ";


// File: classLOCA_1_1EigenvalueSort_1_1LargestImaginary.xml
%feature("docstring") LOCA::EigenvalueSort::LargestImaginary "

Largest-imaginary sorting strategy.

Sorts eigenvalues in decreasing order according to their imaginary
part. This method requires no parameters in the eigenParams argument
to the constructor

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::LargestImaginary::LargestImaginary "LOCA::EigenvalueSort::LargestImaginary::LargestImaginary(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")
LOCA::EigenvalueSort::LargestImaginary::~LargestImaginary "LOCA::EigenvalueSort::LargestImaginary::~LargestImaginary()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestImaginary::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestImaginary::sort(int n, double *evals,
std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestImaginary::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestImaginary::sort(int n, double *r_evals,
double *i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1EigenvalueSort_1_1LargestMagnitude.xml
%feature("docstring") LOCA::EigenvalueSort::LargestMagnitude "

Largest-magnitude sorting strategy.

Sorts eigenvalues in decreasing order according to their magnitude.
This method requires no parameters in the eigenParams argument to the
constructor

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::LargestMagnitude::LargestMagnitude "LOCA::EigenvalueSort::LargestMagnitude::LargestMagnitude(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")
LOCA::EigenvalueSort::LargestMagnitude::~LargestMagnitude "LOCA::EigenvalueSort::LargestMagnitude::~LargestMagnitude()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestMagnitude::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestMagnitude::sort(int n, double *evals,
std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestMagnitude::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestMagnitude::sort(int n, double *r_evals,
double *i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1EigenvalueSort_1_1LargestReal.xml
%feature("docstring") LOCA::EigenvalueSort::LargestReal "

Largest-real sorting strategy.

Sorts eigenvalues in decreasing order according to their real part.
This method requires no parameters in the eigenParams argument to the
constructor

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")  LOCA::EigenvalueSort::LargestReal::LargestReal
"LOCA::EigenvalueSort::LargestReal::LargestReal(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestReal::~LargestReal
"LOCA::EigenvalueSort::LargestReal::~LargestReal()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestReal::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestReal::sort(int n, double *evals,
std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")  LOCA::EigenvalueSort::LargestReal::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestReal::sort(int n, double *r_evals, double
*i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1EigenvalueSort_1_1LargestRealInverseCayley.xml
%feature("docstring") LOCA::EigenvalueSort::LargestRealInverseCayley "

Largest-Real Cayley sorting strategy.

Sorts eigenvalues in decreasing order according to the real part of
their inverse-Cayley transformation. This method references the
\"CayleyPole\" and \"CayleyZero\" parameters in the eigensolver
parameter list.

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::LargestRealInverseCayley::LargestRealInverseCayley
"LOCA::EigenvalueSort::LargestRealInverseCayley::LargestRealInverseCayley(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")
LOCA::EigenvalueSort::LargestRealInverseCayley::~LargestRealInverseCayley
"LOCA::EigenvalueSort::LargestRealInverseCayley::~LargestRealInverseCayley()

Destructor. ";

%feature("docstring")
LOCA::EigenvalueSort::LargestRealInverseCayley::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestRealInverseCayley::sort(int n, double
*evals, std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")
LOCA::EigenvalueSort::LargestRealInverseCayley::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestRealInverseCayley::sort(int n, double
*r_evals, double *i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1Epetra_1_1LeftPreconditionedOp.xml
%feature("docstring") LOCA::Epetra::LeftPreconditionedOp "

An Epetra operator for implementing the operator $P = M^{-1}J$.

This class implements the Epetra_Operator interface for $P = M^{-1}J$
where $J$ and $M$ are Epetra_Operator's.

C++ includes: LOCA_Epetra_LeftPreconditionedOp.H ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::LeftPreconditionedOp "LOCA::Epetra::LeftPreconditionedOp::LeftPreconditionedOp(const
Teuchos::RCP< Epetra_Operator > &jacOperator, const Teuchos::RCP<
Epetra_Operator > &precOperator)

Constructor.

Parameters:
-----------

jacOperator:  [in] Jacobian operator J

precOperator:  [in] Preconditioner operator M ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::~LeftPreconditionedOp "LOCA::Epetra::LeftPreconditionedOp::~LeftPreconditionedOp()

Destructor. ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::SetUseTranspose "int
LOCA::Epetra::LeftPreconditionedOp::SetUseTranspose(bool UseTranspose)

Set to true if the transpose of the operator is requested. ";

%feature("docstring")  LOCA::Epetra::LeftPreconditionedOp::Apply "int
LOCA::Epetra::LeftPreconditionedOp::Apply(const Epetra_MultiVector
&Input, Epetra_MultiVector &Result) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector Input in Result as described above. ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::ApplyInverse "int
LOCA::Epetra::LeftPreconditionedOp::ApplyInverse(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of the inverse of the operator applied to a
Epetra_MultiVector Input in Result as described above. ";

%feature("docstring")  LOCA::Epetra::LeftPreconditionedOp::NormInf "double LOCA::Epetra::LeftPreconditionedOp::NormInf() const

Returns an approximate infinity norm of the operator matrix.

This is defined only if NormInf() of the underlying operators $J$ and
$M$ is defined and is given by
$\\\\|J\\\\|_\\\\infty+\\\\|M\\\\|_\\\\infty$. ";

%feature("docstring")  LOCA::Epetra::LeftPreconditionedOp::Label "const char * LOCA::Epetra::LeftPreconditionedOp::Label() const

Returns a character std::string describing the operator. ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::UseTranspose "bool
LOCA::Epetra::LeftPreconditionedOp::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  LOCA::Epetra::LeftPreconditionedOp::HasNormInf
"bool LOCA::Epetra::LeftPreconditionedOp::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  LOCA::Epetra::LeftPreconditionedOp::Comm "const Epetra_Comm & LOCA::Epetra::LeftPreconditionedOp::Comm() const

Returns a reference to the Epetra_Comm communicator associated with
this operator. ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::OperatorDomainMap "const
Epetra_Map & LOCA::Epetra::LeftPreconditionedOp::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")
LOCA::Epetra::LeftPreconditionedOp::OperatorRangeMap "const
Epetra_Map & LOCA::Epetra::LeftPreconditionedOp::OperatorRangeMap()
const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";


// File: classLOCA_1_1Epetra_1_1TransposeLinearSystem_1_1LeftPreconditioning.xml
%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning "

Method for solving the transpose of a linear system by transposing the
preconditioner and switching to left preconditioning.

C++ includes: LOCA_Epetra_TransposeLinearSystem_LeftPreconditioning.H
";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::LeftPreconditioning
"LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::LeftPreconditioning(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linsys)

Constructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::~LeftPreconditioning
"LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::~LeftPreconditioning()

Destructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::applyJacobianTransposeInverse
"bool
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector
&result)

Applies the inverse of the Jacobian matrix transpose to the given
input vector and puts the answer in result.

Computes \\\\[ v = J^{-T} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The parameter list contains the linear solver options. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::createJacobianTranspose
"bool
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::createJacobianTranspose()

Evaluates the Jacobian-transpose based on the solution vector x.

Note: For flexibility, this method does not compute the original
Jacobian matrix. It uses whatever is currently stored in the linear
system. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::createTransposePreconditioner
"bool
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::createTransposePreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &p)

Explicitly constructs a preconditioner based on the solution vector x
and the parameter list p.

Note: x is only needed for user-supplied preconditioners. When using a
built- in preconditioner (e.g., Ifpack), x will note be used. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::getJacobianTransposeOperator
"Teuchos::RCP< Epetra_Operator >
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::getJacobianTransposeOperator()

Get Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::getTransposePreconditioner
"Teuchos::RCP< Epetra_Operator >
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::getTransposePreconditioner()

Get transpose-preconditioner. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::setJacobianTransposeOperator
"void
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::setJacobianTransposeOperator(const
Teuchos::RCP< Epetra_Operator > &new_jac_trans)

Set Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::setTransposePreconditioner
"void
LOCA::Epetra::TransposeLinearSystem::LeftPreconditioning::setTransposePreconditioner(const
Teuchos::RCP< Epetra_Operator > &new_prec_trans)

Set transpose-preconditioner. ";


// File: classLOCA_1_1Parameter_1_1Library.xml
%feature("docstring") LOCA::Parameter::Library "

Class to provide a centralized library for setting/retrieving
numerical parameter values in application codes.

This class provides a mechanism for setting and retrieving arbitrary
numerical parameter values throughout an application code. Parameters
can be material properties, coefficients in source functions, etc. The
purpose of this class is to allow external libraries to set and
retrieve parameters values to perform, for example, numerical
continuation and optimization.

This class in currently under development and is far from complete.

C++ includes: LOCA_Parameter_Library.H ";

%feature("docstring")  LOCA::Parameter::Library::Library "LOCA::Parameter::Library::Library()

Default constructor. ";

%feature("docstring")  LOCA::Parameter::Library::~Library "LOCA::Parameter::Library::~Library()

Destructor. ";

%feature("docstring")  LOCA::Parameter::Library::setValue "void
LOCA::Parameter::Library::setValue(const std::string &name, const
ValueType &value)

Set parameter given by name to value value. ";

%feature("docstring")  LOCA::Parameter::Library::getValue "ValueType
LOCA::Parameter::Library::getValue(const std::string &name) const

Get parameter given by name. ";

%feature("docstring")  LOCA::Parameter::Library::addParameterEntry "bool LOCA::Parameter::Library::addParameterEntry(const std::string
&name, ObjectType &object, ValueType ObjectType::*object_val_ptr)

Add a new parameter to library using the default setting mechanism.

Returns true if successful in adding entry to library, false
otherwise. ";

%feature("docstring")  LOCA::Parameter::Library::addParameterEntry "bool LOCA::Parameter::Library::addParameterEntry(const std::string
&name, FunctorType *fctr)

Add a new parameter to library using functor setting mechanism.

Returns true if successful in adding entry to library, false
otherwise. ";

%feature("docstring")  LOCA::Parameter::Library::addParameterEntry "bool LOCA::Parameter::Library::addParameterEntry(const std::string
&name, Entry< ValueType > *entry)

Add a new parameter using custom entry.

Returns true if successful in adding entry to library, false
otherwise. ";


// File: classLOCA_1_1BorderedSolver_1_1LowerTriangularBlockElimination.xml
%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination "

Block elimination strategy for solving a block lower-triangular
system.

This class solves the extended system of equations \\\\[
\\\\begin{bmatrix} op(J) & 0 \\\\\\\\ B^T & op(C) \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
F \\\\\\\\ G \\\\end{bmatrix} \\\\] via block elimination: \\\\[
\\\\begin{aligned} X &= op(J)^{-1} F \\\\\\\\ Y &= op(C)^{-1}(G-B^T X)
\\\\end{aligned} \\\\] where $op$ represents either the identity
operation or the transpose. $C$ must be nonzero, while $B$, $F$ or $G$
may be zero. $B$ may be specified either as a
NOX::Abstract::MultiVector or a
LOCA::MultiContinuation::ConstraintInterface object. The solve for the
non-transposed system is implemented by the solve() method, while the
solve for the transposed system is implemented by the solveTranspose()
method.

C++ includes: LOCA_BorderedSolver_LowerTriangularBlockElimination.H ";

%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination::LowerTriangularBlockElimination
"LOCA::BorderedSolver::LowerTriangularBlockElimination::LowerTriangularBlockElimination(const
Teuchos::RCP< LOCA::GlobalData > &global_data)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object ";

%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination::~LowerTriangularBlockElimination
"LOCA::BorderedSolver::LowerTriangularBlockElimination::~LowerTriangularBlockElimination()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination::solve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::solve(Teuchos::ParameterList
&params, const LOCA::BorderedSolver::AbstractOperator &op, const
LOCA::MultiContinuation::ConstraintInterface &B, const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system as described above with B specified as a
LOCA::MultiContinuation::ConstraintInterface object.

Either F, or G may be zero by passing NULL. ";

%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination::solve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::solve(Teuchos::ParameterList
&params, const LOCA::BorderedSolver::AbstractOperator &op, const
NOX::Abstract::MultiVector &B, const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system as described above with B specified as a
NOX::Abstract::MultiVector.

Either F, or G may be zero by passing NULL. ";

%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination::solveTranspose
"NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::solveTranspose(Teuchos::ParameterList
&params, const LOCA::BorderedSolver::AbstractOperator &op, const
LOCA::MultiContinuation::ConstraintInterface &B, const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system using the tranpose of J and C as described
above with B specified as a
LOCA::MultiContinuation::ConstraintInterface object.

Either F, or G may be zero by passing NULL. ";

%feature("docstring")
LOCA::BorderedSolver::LowerTriangularBlockElimination::solveTranspose
"NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LowerTriangularBlockElimination::solveTranspose(Teuchos::ParameterList
&params, const LOCA::BorderedSolver::AbstractOperator &op, const
NOX::Abstract::MultiVector &B, const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system using the tranpose of J and C as described
above with B specified as a NOX::Abstract::MultiVector object.

Either F, or G may be zero by passing NULL. ";


// File: classLOCA_1_1Epetra_1_1LowRankUpdateOp.xml
%feature("docstring") LOCA::Epetra::LowRankUpdateOp "

An Epetra operator for implementing the operator $P = J + U V^T$.

This class implements the Epetra_Operator interface for $P = J + U
V^T$ where $J$ is an Epetra_Operator and $U$ and $V$ are
Epetra_MultiVectors.

C++ includes: LOCA_Epetra_LowRankUpdateOp.H ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::LowRankUpdateOp
"LOCA::Epetra::LowRankUpdateOp::LowRankUpdateOp(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP< Epetra_Operator >
&jacOperator, const Teuchos::RCP< const Epetra_MultiVector >
&U_multiVec, const Teuchos::RCP< const Epetra_MultiVector >
&V_multiVec, bool setup_for_solve)

Constructor.

Parameters:
-----------

global_data:  [in] The global data object

jacOperator:  [in] Jacobian operator J

U_multiVec:  [in] Multivector representing U

V_multiVec:  [in] Multivector representing V ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::~LowRankUpdateOp
"LOCA::Epetra::LowRankUpdateOp::~LowRankUpdateOp()

Destructor. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::SetUseTranspose
"int LOCA::Epetra::LowRankUpdateOp::SetUseTranspose(bool
UseTranspose)

Set to true if the transpose of the operator is requested. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::Apply "int
LOCA::Epetra::LowRankUpdateOp::Apply(const Epetra_MultiVector &Input,
Epetra_MultiVector &Result) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector Input in Result as described above. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::ApplyInverse "int LOCA::Epetra::LowRankUpdateOp::ApplyInverse(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const

This method does nothing. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::NormInf "double
LOCA::Epetra::LowRankUpdateOp::NormInf() const

Returns an approximate infinity norm of the operator matrix.

This is defined only if NormInf() of the underlying operator $J$ is
defined and is given by
$\\\\|J\\\\|_\\\\infty+\\\\|U\\\\|_\\\\infty\\\\|V\\\\|_\\\\infty$. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::Label "const
char * LOCA::Epetra::LowRankUpdateOp::Label() const

Returns a character std::string describing the operator. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::UseTranspose "bool LOCA::Epetra::LowRankUpdateOp::UseTranspose() const

Returns the current UseTranspose setting. Always returns false. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::HasNormInf "bool LOCA::Epetra::LowRankUpdateOp::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::Comm "const
Epetra_Comm & LOCA::Epetra::LowRankUpdateOp::Comm() const

Returns a reference to the Epetra_Comm communicator associated with
this operator. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateOp::OperatorDomainMap "const Epetra_Map &
LOCA::Epetra::LowRankUpdateOp::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateOp::OperatorRangeMap
"const Epetra_Map & LOCA::Epetra::LowRankUpdateOp::OperatorRangeMap()
const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";


// File: classLOCA_1_1Epetra_1_1LowRankUpdateRowMatrix.xml
%feature("docstring") LOCA::Epetra::LowRankUpdateRowMatrix "

An Epetra row matrix for implementing the operator $P = J + U V^T$.

This class implements the Epetra_RowMatrix interface for $P = J + U
V^T$ where $J$ is an Epetra_RowMatrix and $U$ and $V$ are
Epetra_MultiVectors. It is derived from LOCA::Epetra::LowRankUpdateOp
to implement the Epetra_Operator interface. The interface here
implements the Epetra_RowMatrix interface when the matrix $J$ is
itself a row matrix. This allows preconditioners to be computed and
scaling in linear systems to be performed when using this operator.
The implementation here merely adds the corresponding entries for $U
V^T$ to the rows of $J$. Note however this is only an approximation to
the true matrix $J + U V^T$.

This class assumes $U$ and $V$ have the same distribution as the rows
of $J$.

C++ includes: LOCA_Epetra_LowRankUpdateRowMatrix.H ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumMyRowEntries "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumMyRowEntries(int MyRow, int
&NumEntries) const

Returns the number of nonzero entries in MyRow. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::MaxNumEntries "int
LOCA::Epetra::LowRankUpdateRowMatrix::MaxNumEntries() const

Returns the maximum of NumMyRowEntries() over all rows. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::ExtractMyRowCopy "int
LOCA::Epetra::LowRankUpdateRowMatrix::ExtractMyRowCopy(int MyRow, int
Length, int &NumEntries, double *Values, int *Indices) const

Returns a copy of the specified local row in user-provided arrays. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::ExtractDiagonalCopy "int
LOCA::Epetra::LowRankUpdateRowMatrix::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const

Returns a copy of the main diagonal in a user-provided vector. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::Multiply
"int LOCA::Epetra::LowRankUpdateRowMatrix::Multiply(bool TransA,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_RowMatrix multiplied by a
Epetra_MultiVector X in Y. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::Solve "int LOCA::Epetra::LowRankUpdateRowMatrix::Solve(bool Upper, bool
Trans, bool UnitDiagonal, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns result of a local-only solve using a triangular
Epetra_RowMatrix with Epetra_MultiVectors X and Y. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::InvRowSums "int
LOCA::Epetra::LowRankUpdateRowMatrix::InvRowSums(Epetra_Vector &x)
const

Computes the sum of absolute values of the rows of the
Epetra_RowMatrix, results returned in x. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::LeftScale
"int LOCA::Epetra::LowRankUpdateRowMatrix::LeftScale(const
Epetra_Vector &x)

Scales the Epetra_RowMatrix on the left with a Epetra_Vector x. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::InvColSums "int
LOCA::Epetra::LowRankUpdateRowMatrix::InvColSums(Epetra_Vector &x)
const

Computes the sum of absolute values of the columns of the
Epetra_RowMatrix, results returned in x. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::RightScale "int
LOCA::Epetra::LowRankUpdateRowMatrix::RightScale(const Epetra_Vector
&x)

Scales the Epetra_RowMatrix on the right with a Epetra_Vector x. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::Filled "bool LOCA::Epetra::LowRankUpdateRowMatrix::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::NormInf "double LOCA::Epetra::LowRankUpdateRowMatrix::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::NormOne "double LOCA::Epetra::LowRankUpdateRowMatrix::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalNonzeros "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalNonzeros64 "long long
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalNonzeros64() const ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalRows "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalRows64 "long long
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalRows64() const ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalCols "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalCols64 "long long
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalCols64() const ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalDiagonals "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalDiagonals64 "long long
LOCA::Epetra::LowRankUpdateRowMatrix::NumGlobalDiagonals64() const ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumMyNonzeros "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumMyNonzeros() const

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::NumMyRows
"int LOCA::Epetra::LowRankUpdateRowMatrix::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::NumMyCols
"int LOCA::Epetra::LowRankUpdateRowMatrix::NumMyCols() const

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::NumMyDiagonals "int
LOCA::Epetra::LowRankUpdateRowMatrix::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::LowerTriangular "bool
LOCA::Epetra::LowRankUpdateRowMatrix::LowerTriangular() const

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::UpperTriangular "bool
LOCA::Epetra::LowRankUpdateRowMatrix::UpperTriangular() const

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::RowMatrixRowMap "const
Epetra_Map & LOCA::Epetra::LowRankUpdateRowMatrix::RowMatrixRowMap()
const

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::RowMatrixColMap "const
Epetra_Map & LOCA::Epetra::LowRankUpdateRowMatrix::RowMatrixColMap()
const

Returns the Epetra_Map object associated with the columns of this
matrix. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::RowMatrixImporter "const
Epetra_Import *
LOCA::Epetra::LowRankUpdateRowMatrix::RowMatrixImporter() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::LowRankUpdateRowMatrix "LOCA::Epetra::LowRankUpdateRowMatrix::LowRankUpdateRowMatrix(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Epetra_RowMatrix > &jacRowMatrix, const Teuchos::RCP<
Epetra_MultiVector > &U_multiVec, const Teuchos::RCP<
Epetra_MultiVector > &V_multiVec, bool setup_for_solve, bool
include_UV_terms)

Constructor.

Parameters:
-----------

global_data:  [in] The global data object

jacRowMatrix:  [in] Jacobian operator J as a row matrix

U_multiVec:  [in] Multivector representing U

V_multiVec:  [in] Multivector representing V

setup_for_solve:  [in] Setup data structures for ApplyInverse()

include_UV_terms:  [in] Include $U V^T$ terms in RowMatrix routines
ExtractRowCopy(), ExtactDiagonalCopy(), InvRowSums(), InvColSums(),
NormInf() and NormOne(). ";

%feature("docstring")
LOCA::Epetra::LowRankUpdateRowMatrix::~LowRankUpdateRowMatrix "LOCA::Epetra::LowRankUpdateRowMatrix::~LowRankUpdateRowMatrix()

Destructor. ";

%feature("docstring")  LOCA::Epetra::LowRankUpdateRowMatrix::Map "const Epetra_BlockMap & LOCA::Epetra::LowRankUpdateRowMatrix::Map()
const

Returns a reference to the Epetra_BlockMap for this object. ";


// File: classLOCA_1_1SingularJacobianSolve_1_1Manager.xml
%feature("docstring") LOCA::SingularJacobianSolve::Manager "

Manager for all singular Jacobian solve computations

The parameters passed to the constructor or reset should specify the
\"Method\", as described below, as well as any additional parameters
for that particular method.

\"Method\" - Name of the singular jacobian solve method. Valid choices
are \"Default\" ( LOCA::SingularJacobianSolve::Default) [ Default]

\"Nic\" ( LOCA::SingularJacobianSolve::Nic)

\"Nic-Day\" ( LOCA::SingularJacobianSolve::NicDay)

\"Iterative Refinement\" ( LOCA::SingularJacobianSolve::ItRef)

C++ includes: LOCA_SingularJacobianSolve_Manager.H ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::Manager "LOCA::SingularJacobianSolve::Manager::Manager(Teuchos::ParameterList
&params)

Constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::Manager "LOCA::SingularJacobianSolve::Manager::Manager(const
Teuchos::ParameterList &params=Teuchos::ParameterList())

Constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::Manager "LOCA::SingularJacobianSolve::Manager::Manager(const Manager &source)

Copy constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::~Manager
"LOCA::SingularJacobianSolve::Manager::~Manager()

Destructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::clone "LOCA::SingularJacobianSolve::Generic *
LOCA::SingularJacobianSolve::Manager::clone() const

Clone function. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::reset "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Manager::reset(Teuchos::ParameterList
&params)

Reset parameters. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::compute "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Manager::compute(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector &input, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector &result)

Computes solution based on method parameter. ";

%feature("docstring")
LOCA::SingularJacobianSolve::Manager::computeMulti "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Manager::computeMulti(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector *const *inputs, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector **results, int nVecs)

Computes solution based on method parameter for multiple RHS. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Manager::getMethod
"const std::string &
LOCA::SingularJacobianSolve::Manager::getMethod() const

Returns current method. ";


// File: classLOCA_1_1Epetra_1_1Interface_1_1MassMatrix.xml
%feature("docstring") LOCA::Epetra::Interface::MassMatrix "

Used by LOCA::Epetra::Group to provide a link to the external code for
the MassMatrix (coefficients of time dependent terms).

This is used for Hopf bifurcation tracking, linear stability analysis,
and space-time solutions (xyzt).

C++ includes: LOCA_Epetra_Interface_MassMatrix.H ";

%feature("docstring")  LOCA::Epetra::Interface::MassMatrix::MassMatrix
"LOCA::Epetra::Interface::MassMatrix::MassMatrix()

Constructor. ";

%feature("docstring")
LOCA::Epetra::Interface::MassMatrix::~MassMatrix "virtual
LOCA::Epetra::Interface::MassMatrix::~MassMatrix()

Destructor. ";

%feature("docstring")
LOCA::Epetra::Interface::MassMatrix::computeMassMatrix "virtual bool
LOCA::Epetra::Interface::MassMatrix::computeMassMatrix(const
Epetra_Vector &x)=0

Compute MassMatrix given the specified input vector x. Returns true if
computation was successful. ";

%feature("docstring")
LOCA::Epetra::Interface::MassMatrix::setOldSolution "virtual void
LOCA::Epetra::Interface::MassMatrix::setOldSolution(const
Epetra_Vector &x, const int timeStep)

Routines used in XYZT to set the old solution, the one from the
previous time step.

There is a different routine for first step, where the old solution is
set by the user and not part of the solution vector. These routines
are used by space-time (xyzt) problems, where the residual vector is a
function of the previous solution, which is also being solved for, and
where the MassMatrix is calculated as a function of a different
solution vector then the Jacobian (that is, the previous time step).
The timeStep argument is sent so the use can set the global time, for
cases when computeF, computeJacobian, computeMassMatrix fills are
functions of time (nonautonomous systems). ";

%feature("docstring")
LOCA::Epetra::Interface::MassMatrix::setOldSolutionFirstStep "virtual
void LOCA::Epetra::Interface::MassMatrix::setOldSolutionFirstStep()

See setOldSolution description. ";

%feature("docstring")
LOCA::Epetra::Interface::MassMatrix::dataForPrintSolution "virtual
void LOCA::Epetra::Interface::MassMatrix::dataForPrintSolution(const
int conStep, const int timeStep, const int totalTimeSteps)

Provides data to application for output files.

This routine is called from Interface::xyzt::printSolution() just
before the call to Interface::Required::printSolution(x,param), and
gives the application some indices that can be used for creating a
unique name/index for the output files. ";


// File: classLOCA_1_1StatusTest_1_1MaxIters.xml
%feature("docstring") LOCA::StatusTest::MaxIters "

Failure test based on the maximum number of continuation steps.

Let $k$ denote the current number of iterations (accessed via
LOCA::Stepper::getNumTotalSteps) and $k_{\\\\max}$ denote the
tolerance set in the constructor of this status test. This test
returns LOCA::StatusTest::Failed if $ k \\\\geq k_{\\\\rm max}. $
Otherwise, it returns LOCA::StatusTest::NotFinished.

If checkStatus is called with the type set to LOCA::StatusTest::None,
it then the status is set to to LOCA::Status::Unevaluated and
returned. (Also #niters is set to -1.)

C++ includes: LOCA_StatusTest_MaxIters.H ";

%feature("docstring")  LOCA::StatusTest::MaxIters::MaxIters "LOCA::StatusTest::MaxIters::MaxIters(int maxIterations, bool
return_failed_on_max_steps=true, const Teuchos::RCP< const
LOCA::GlobalData > globalData=Teuchos::null)

Constructor. Specify the maximum number of nonlinear solver
iterations, $k_{\\\\max}$ ands optinally an error stream for printing
errors. ";

%feature("docstring")  LOCA::StatusTest::MaxIters::~MaxIters "LOCA::StatusTest::MaxIters::~MaxIters()

Destructor. ";

%feature("docstring")  LOCA::StatusTest::MaxIters::checkStatus "virtual LOCA::StatusTest::StatusType
LOCA::StatusTest::MaxIters::checkStatus(const LOCA::Abstract::Iterator
&stepper, LOCA::StatusTest::CheckType status)

Test the stopping criterion

The test can (and should, if possible) be skipped if checkType is
LOCA::StatusType::None. If the test is skipped, then the status should
be set to LOCA::StatusTest::Unevaluated. ";

%feature("docstring")  LOCA::StatusTest::MaxIters::getStatus "LOCA::StatusTest::StatusType LOCA::StatusTest::MaxIters::getStatus()
const

Return the result of the most recent checkStatus call. ";

%feature("docstring")  LOCA::StatusTest::MaxIters::print "std::ostream & LOCA::StatusTest::MaxIters::print(std::ostream &stream,
int indent=0) const

Output formatted description of stopping test to output stream. ";

%feature("docstring")  LOCA::StatusTest::MaxIters::getMaxIters "int
LOCA::StatusTest::MaxIters::getMaxIters() const

Returns the Maximum number of iterations set in the constructor. ";

%feature("docstring")  LOCA::StatusTest::MaxIters::getNumIters "int
LOCA::StatusTest::MaxIters::getNumIters() const

Returns the current number of iterations taken by the solver.

Returns -1 if the status of this test is NOX::StatusTest::Unevaluated.
";


// File: classLOCA_1_1Epetra_1_1ModelEvaluatorInterface.xml
%feature("docstring") LOCA::Epetra::ModelEvaluatorInterface "

Wrapper for an EpetraExt::ModelEvaluator.

If an application interfaces their code to solvers using the
EpetraExt::ModelEvaluator, this class provides a wrapper so that the
model evaluator can be used instead of having the user write concrete
versions of the LOCA::Epetra::Interface objects.

C++ includes: LOCA_Epetra_ModelEvaluatorInterface.H ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::ModelEvaluatorInterface "LOCA::Epetra::ModelEvaluatorInterface::ModelEvaluatorInterface(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const
Teuchos::RefCountPtr< EpetraExt::ModelEvaluator > &m, double
perturb=1.0e-6)

Constructor. ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::~ModelEvaluatorInterface "LOCA::Epetra::ModelEvaluatorInterface::~ModelEvaluatorInterface()

Destructor. ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::getLOCAParameterVector "const
LOCA::ParameterVector &
LOCA::Epetra::ModelEvaluatorInterface::getLOCAParameterVector() const

Return LOCA parameter vector. ";

%feature("docstring")  LOCA::Epetra::ModelEvaluatorInterface::computeF
"bool LOCA::Epetra::ModelEvaluatorInterface::computeF(const
Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag) ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::computeJacobian "bool
LOCA::Epetra::ModelEvaluatorInterface::computeJacobian(const
Epetra_Vector &x, Epetra_Operator &Jac) ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::computePreconditioner "bool
LOCA::Epetra::ModelEvaluatorInterface::computePreconditioner(const
Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList
*precParams=0) ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::setParameters "void
LOCA::Epetra::ModelEvaluatorInterface::setParameters(const
ParameterVector &p)

Set parameters in the user's application.

Should be called prior to calling one of the compute functions. ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::computeShiftedMatrix "bool
LOCA::Epetra::ModelEvaluatorInterface::computeShiftedMatrix(double
alpha, double beta, const Epetra_Vector &x, Epetra_Operator &A)

Call user routine for computing the shifted matrix $\\\\alpha J +
\\\\beta M$ where $J$ is the Jacobian matrix and $M$ is the mass
matrix. ";

%feature("docstring")  LOCA::Epetra::ModelEvaluatorInterface::setXdot
"void LOCA::Epetra::ModelEvaluatorInterface::setXdot(const
Epetra_Vector &xdot, const double time)

Routine used in XYZT to set x_dot and time in the interface.

The computeF() routine for XYZT problems needs to be a function of
x_dot, but th NOX/LOCA computeF() does not take x_dot as an argument.
This is used to set x_dot in the application interface so the
subsequent call to computeF has the correct x_dot value. The timeStep
argument is sent so the use can set the global time, for cases when
computeF, computeJacobian, computeMassMatrix fills are functions of
time (nonautonomous systems). ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::printSolution "void
LOCA::Epetra::ModelEvaluatorInterface::printSolution(const
Epetra_Vector &x_, double conParam)

Call user's own print routine for vector-parameter pair. ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::setObserver "void
LOCA::Epetra::ModelEvaluatorInterface::setObserver(const Teuchos::RCP<
NOX::Epetra::Observer > &observer_) ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::ModelEvaluatorInterface "LOCA::Epetra::ModelEvaluatorInterface::ModelEvaluatorInterface(const
ModelEvaluatorInterface &)

Copy constructor. ";

%feature("docstring")  LOCA::Epetra::ModelEvaluatorInterface::clone "Teuchos::RCP< LOCA::DerivUtils >
LOCA::Epetra::ModelEvaluatorInterface::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone. ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::computeDfDp "NOX::Abstract::Group::ReturnType
LOCA::Epetra::ModelEvaluatorInterface::computeDfDp(LOCA::MultiContinuation::AbstractGroup
&grp, const std::vector< int > &param_ids, NOX::Abstract::MultiVector
&result, bool isValidF) const

Compute derivative of f with respect to parameter, identified by
param_id. ";

%feature("docstring")
LOCA::Epetra::ModelEvaluatorInterface::postProcessContinuationStep "void
LOCA::Epetra::ModelEvaluatorInterface::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus, LOCA::Epetra::Group &group)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. ";


// File: classLOCA_1_1TurningPoint_1_1MinimallyAugmented_1_1ModifiedConstraint.xml
%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint "

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
computing turning points for the minimally augmented turning point
formulation.

This class is a modification of
LOCA::TurningPoint::MinimallyAugmented::Constraint where updates are
computed to the left and right null vectors $w$ and $v$ every
nonlinear iteration instead of solving for them directly: \\\\[
\\\\begin{bmatrix} J & a \\\\\\\\ b^T & 0 \\\\end{bmatrix}
\\\\begin{bmatrix} \\\\Delta v \\\\\\\\ \\\\Delta \\\\sigma_1
\\\\end{bmatrix} = - \\\\begin{bmatrix} J v + \\\\sigma_1 a +
(Jv)_x\\\\Delta x + (Jv)_p \\\\Delta p\\\\\\\\ b^T a - n
\\\\end{bmatrix}, \\\\] \\\\[ \\\\begin{bmatrix} J^T & b \\\\\\\\ a^T
& 0 \\\\end{bmatrix} \\\\begin{bmatrix} \\\\Delta w \\\\\\\\ \\\\Delta
\\\\sigma_2 \\\\end{bmatrix} = - \\\\begin{bmatrix} J^T w +
\\\\sigma_2 b + (J^T w)_x \\\\Delta x + (J^T w)_p \\\\Delta p \\\\\\\\
a^T w - n \\\\end{bmatrix}, \\\\]

The class is intialized via the tpParams parameter list argument to
the constructor. This class recognizes all paramters for
LOCA::TurningPoint::MinimallyAugmented::Constraint plus the following:
\"Include Newton Terms\" -- [bool] (default: false) - Flag indicating
whether to include the $\\\\Delta x$ and $\\\\Delta p$ terms above
when computing the null vector updates

C++ includes:
LOCA_TurningPoint_MinimallyAugmented_ModifiedConstraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::copy "void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::copy(const
LOCA::MultiContinuation::ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::clone "Teuchos::RCP< LOCA::MultiContinuation::ConstraintInterface >
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::computeConstraints
"NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::preProcessContinuationStep
"void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. Here we set up the constraint class to solve for $w$ and
$v$ for the first nonlinear iteration. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::postProcessContinuationStep
"void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful.
Here we set up the constraint class to solve for $w$ and $v$ for the
first nonlinear iteration. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::ModifiedConstraint
"LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::ModifiedConstraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &tpParams, const Teuchos::RCP<
LOCA::TurningPoint::MinimallyAugmented::AbstractGroup > &g, int
bif_param)

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::ModifiedConstraint
"LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::ModifiedConstraint(const
ModifiedConstraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::~ModifiedConstraint
"LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::~ModifiedConstraint()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::setNewtonUpdates
"void
LOCA::TurningPoint::MinimallyAugmented::ModifiedConstraint::setNewtonUpdates(const
NOX::Abstract::Vector &dx, double dp, double step)

Set the newton update for x and p. ";


// File: classLOCA_1_1Extended_1_1MultiAbstractGroup.xml
%feature("docstring") LOCA::Extended::MultiAbstractGroup "

LOCA abstract interface for extended groups, derived from the
NOX::Abstract::Group, i.e., an abstract interface for \"super\" groups
that have an underlying group component.

Concrete implemenations of this interface must provide implementations
of all of the methods in the NOX::Abstract::Group interface as well as
the additional interface defined here.

C++ includes: LOCA_Extended_MultiAbstractGroup.H ";

/*  Pure virtual methods  */

/* These methods must be defined by any concrete implementation

*/

%feature("docstring")
LOCA::Extended::MultiAbstractGroup::getUnderlyingGroup "virtual
Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getUnderlyingGroup() const =0

Return underlying group.

This method should the underlying group data member. ";

%feature("docstring")
LOCA::Extended::MultiAbstractGroup::getUnderlyingGroup "virtual
Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getUnderlyingGroup()=0

Return underlying group.

This method should the underlying group data member. ";

/*  Virtual methods with default implementations  */

/* These methods should be overloaded in a concrete implementation if
more appropriate/efficient approaches are available.

*/

%feature("docstring")
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup()
const

Return base-level underlying group.

This method is intended for composite groups (such as extended
bifurcation groups) which have an underlying group as a data member.
This method is supposed to return the base level group and has a
default recursive implementation that should work in most cases. ";

%feature("docstring")
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup "Teuchos::RCP< LOCA::MultiContinuation::AbstractGroup >
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup()

Return base-level underlying group.

This method is intended for composite groups (such as extended
bifurcation groups) which have an underlying group as a data member.
This method is supposed to return the base level group and has a
default recursive implementation that should work in most cases. ";

%feature("docstring")
LOCA::Extended::MultiAbstractGroup::MultiAbstractGroup "LOCA::Extended::MultiAbstractGroup::MultiAbstractGroup()

Default constructor. ";

%feature("docstring")
LOCA::Extended::MultiAbstractGroup::~MultiAbstractGroup "virtual
LOCA::Extended::MultiAbstractGroup::~MultiAbstractGroup()

Destructor. ";


// File: classLOCA_1_1MultiContinuation_1_1MultiVecConstraint.xml
%feature("docstring") LOCA::MultiContinuation::MultiVecConstraint "

Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for
a simple linear multivector constraint.

C++ includes: LOCA_MultiContinuation_MultiVecConstraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX
*/

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::copy "void
LOCA::MultiContinuation::MultiVecConstraint::copy(const
ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::clone "Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface >
LOCA::MultiContinuation::MultiVecConstraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::numConstraints "int
LOCA::MultiContinuation::MultiVecConstraint::numConstraints() const

Return number of constraints. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::setX "void
LOCA::MultiContinuation::MultiVecConstraint::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::setParam "void
LOCA::MultiContinuation::MultiVecConstraint::setParam(int paramID,
double val)

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::setParams "void
LOCA::MultiContinuation::MultiVecConstraint::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::computeConstraints "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::MultiVecConstraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::MultiVecConstraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::MultiVecConstraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::isConstraints "bool
LOCA::MultiContinuation::MultiVecConstraint::isConstraints() const

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::isDX "bool
LOCA::MultiContinuation::MultiVecConstraint::isDX() const

Return true if derivatives of constraints w.r.t. x are valid. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::getConstraints "const
NOX::Abstract::MultiVector::DenseMatrix &
LOCA::MultiContinuation::MultiVecConstraint::getConstraints() const

Return constraint residuals. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::getDX "const
NOX::Abstract::MultiVector *
LOCA::MultiContinuation::MultiVecConstraint::getDX() const

Return solution component of constraint derivatives. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::isDXZero "bool
LOCA::MultiContinuation::MultiVecConstraint::isDXZero() const

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::notifyCompletedStep "void
LOCA::MultiContinuation::MultiVecConstraint::notifyCompletedStep()

Notify constraint that the continuation step is completed.

Here we do nothing ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::MultiVecConstraint "LOCA::MultiContinuation::MultiVecConstraint::MultiVecConstraint(const
Teuchos::RCP< const NOX::Abstract::MultiVector > &dx)

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::MultiVecConstraint "LOCA::MultiContinuation::MultiVecConstraint::MultiVecConstraint(const
MultiVecConstraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::~MultiVecConstraint "LOCA::MultiContinuation::MultiVecConstraint::~MultiVecConstraint()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::MultiVecConstraint::setDx "void
LOCA::MultiContinuation::MultiVecConstraint::setDx(const Teuchos::RCP<
const NOX::Abstract::MultiVector > &dx)

Set constraint vector. ";


// File: classLOCA_1_1Extended_1_1MultiVector.xml
%feature("docstring") LOCA::Extended::MultiVector "

Implemenatation of the NOX::Abstract::MultiVector class for extended
multi-vectors comprised of an arbitrary number of multi-vectors and
scalars.

C++ includes: LOCA_Extended_MultiVector.H ";

%feature("docstring")  LOCA::Extended::MultiVector::init "NOX::Abstract::MultiVector & LOCA::Extended::MultiVector::init(double
gamma)

Initialize every element of this multi-vector with gamma. ";

%feature("docstring")  LOCA::Extended::MultiVector::random "NOX::Abstract::MultiVector & LOCA::Extended::MultiVector::random(bool
useSeed=false, int seed=1)

Initialize each element of this multi-vector with a random value. ";

%feature("docstring")  LOCA::Extended::MultiVector::setBlock "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::setBlock(const NOX::Abstract::MultiVector
&source, const std::vector< int > &index)

Copy the vectors in source to a set of vectors in *this. The
index.size() vectors in source are copied to a subset of vectors in
*this indicated by the indices given in index. ";

%feature("docstring")  LOCA::Extended::MultiVector::setBlock "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::setBlock(const MultiVector &source, const
std::vector< int > &index) ";

%feature("docstring")  LOCA::Extended::MultiVector::augment "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::augment(const NOX::Abstract::MultiVector
&source)

Append the vectors in source to *this. ";

%feature("docstring")  LOCA::Extended::MultiVector::augment "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::augment(const MultiVector &source) ";

%feature("docstring")  LOCA::Extended::MultiVector::scale "NOX::Abstract::MultiVector & LOCA::Extended::MultiVector::scale(double
gamma)

Scale each element of this multivector by gamma. ";

%feature("docstring")  LOCA::Extended::MultiVector::update "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::update(double alpha, const
NOX::Abstract::MultiVector &a, double gamma=0.0)

Compute x = (alpha * a) + (gamma * x) where a is a multi-vector and x
= *this. ";

%feature("docstring")  LOCA::Extended::MultiVector::update "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::update(double alpha, const MultiVector
&a, double gamma=0.0) ";

%feature("docstring")  LOCA::Extended::MultiVector::update "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::update(double alpha, const
NOX::Abstract::MultiVector &a, double beta, const
NOX::Abstract::MultiVector &b, double gamma=0.0)

Compute x = (alpha * a) + (beta * b) + (gamma * x) where a and b are
multi-vectors and x = *this. ";

%feature("docstring")  LOCA::Extended::MultiVector::update "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::update(double alpha, const MultiVector
&a, double beta, const MultiVector &b, double gamma=0.0) ";

%feature("docstring")  LOCA::Extended::MultiVector::update "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::update(Teuchos::ETransp transb, double
alpha, const NOX::Abstract::MultiVector &a, const
NOX::Abstract::MultiVector::DenseMatrix &b, double gamma=0.0)

Compute x = (alpha * a * b) + (gamma * x) where a is a multivector, b
is a dense matrix, x = *this, and op(b) = b if transb =
Teuchos::NO_TRANS and op(b) is b transpose if transb = Teuchos::TRANS.
";

%feature("docstring")  LOCA::Extended::MultiVector::update "NOX::Abstract::MultiVector &
LOCA::Extended::MultiVector::update(Teuchos::ETransp transb, double
alpha, const MultiVector &a, const
NOX::Abstract::MultiVector::DenseMatrix &b, double gamma=0.0) ";

%feature("docstring")  LOCA::Extended::MultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::MultiVector::clone(NOX::CopyType type=NOX::DeepCopy)
const

Create a new MultiVector of the same underlying type by cloning
\"this\", and return a pointer to the new vector.

If type is NOX::DeepCopy, then we need to create an exact replica of
\"this\". Otherwise, if type is NOX::ShapeCopy, we need only replicate
the shape of \"this\". Note that there is no assumption that a vector
created by ShapeCopy is initialized to zeros.

Pointer to newly created vector or NULL if clone is not supported. ";

%feature("docstring")  LOCA::Extended::MultiVector::clone "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::MultiVector::clone(int numvecs) const

Creates a new multi-vector with numvecs columns. ";

%feature("docstring")  LOCA::Extended::MultiVector::subCopy "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::MultiVector::subCopy(const std::vector< int > &index)
const

Creates a new multi-vector with index.size() columns whose columns are
copies of the columns of *this given by index. ";

%feature("docstring")  LOCA::Extended::MultiVector::subView "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::MultiVector::subView(const std::vector< int > &index)
const

Creates a new multi-vector with index.size() columns that shares the
columns of *this given by index. ";

%feature("docstring")  LOCA::Extended::MultiVector::norm "void
LOCA::Extended::MultiVector::norm(std::vector< double > &result,
NOX::Abstract::Vector::NormType type=NOX::Abstract::Vector::TwoNorm)
const

Norm. ";

%feature("docstring")  LOCA::Extended::MultiVector::multiply "void
LOCA::Extended::MultiVector::multiply(double alpha, const
NOX::Abstract::MultiVector &y, NOX::Abstract::MultiVector::DenseMatrix
&b) const

Computes the matrix-matrix product $\\\\alpha * y^T * (*this)$. ";

%feature("docstring")  LOCA::Extended::MultiVector::multiply "void
LOCA::Extended::MultiVector::multiply(double alpha, const MultiVector
&y, NOX::Abstract::MultiVector::DenseMatrix &b) const ";

%feature("docstring")  LOCA::Extended::MultiVector::length "int
LOCA::Extended::MultiVector::length() const

Return the length of multi-vector. ";

%feature("docstring")  LOCA::Extended::MultiVector::numVectors "int
LOCA::Extended::MultiVector::numVectors() const

Return the number of vectors in the multi-vector. ";

%feature("docstring")  LOCA::Extended::MultiVector::print "void
LOCA::Extended::MultiVector::print(std::ostream &stream) const

Print the vector. This is meant for debugging purposes only. ";

%feature("docstring")  LOCA::Extended::MultiVector::getMultiVector "Teuchos::RCP< const NOX::Abstract::MultiVector >
LOCA::Extended::MultiVector::getMultiVector(int i) const

Returns const ref-count pointer to the ith multi-vector. ";

%feature("docstring")  LOCA::Extended::MultiVector::getMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::MultiVector::getMultiVector(int i)

Returns ref-count pointer to the ith multi-vector. ";

%feature("docstring")  LOCA::Extended::MultiVector::getScalars "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Extended::MultiVector::getScalars() const

Returns const ref-count pointer to scalar matrix. ";

%feature("docstring")  LOCA::Extended::MultiVector::getScalars "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Extended::MultiVector::getScalars()

Returns ref-count pointer to scalar matrix. ";

%feature("docstring")  LOCA::Extended::MultiVector::getScalarRows "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Extended::MultiVector::getScalarRows(int num_rows, int row)
const

Returns const ref-count pointer to num_rows rows of scalar matrix
starting at row row. ";

%feature("docstring")  LOCA::Extended::MultiVector::getScalarRows "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Extended::MultiVector::getScalarRows(int num_rows, int row)

Returns ref-count pointer to num_rows rows of scalar matrix starting
at row row. ";

%feature("docstring")  LOCA::Extended::MultiVector::getScalar "const
double & LOCA::Extended::MultiVector::getScalar(int i, int j) const

Returns const reference to the scalar for row i, column j. ";

%feature("docstring")  LOCA::Extended::MultiVector::getScalar "double
& LOCA::Extended::MultiVector::getScalar(int i, int j)

Returns reference to the scalar for row i, column j. ";

%feature("docstring")  LOCA::Extended::MultiVector::getVector "Teuchos::RCP< LOCA::Extended::Vector >
LOCA::Extended::MultiVector::getVector(int i)

Return a ref-count pointer to the i-th column of the multivector as an
abstract vector. ";

%feature("docstring")  LOCA::Extended::MultiVector::getVector "Teuchos::RCP< const LOCA::Extended::Vector >
LOCA::Extended::MultiVector::getVector(int i) const

Return a const ref-count pointer to the i-th column of the multivector
as an abstract vector. ";

%feature("docstring")  LOCA::Extended::MultiVector::getNumScalarRows "int LOCA::Extended::MultiVector::getNumScalarRows() const

Returns number of scalars rows. ";

%feature("docstring")  LOCA::Extended::MultiVector::getNumMultiVectors
"int LOCA::Extended::MultiVector::getNumMultiVectors() const

Returns number of multi vectors. ";

%feature("docstring")  LOCA::Extended::MultiVector::MultiVector "LOCA::Extended::MultiVector::MultiVector(const MultiVector &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::Extended::MultiVector::MultiVector "LOCA::Extended::MultiVector::MultiVector(const MultiVector &source,
int nColumns)

Copy constructor that creates a new multivector with nColumns columns.
";

%feature("docstring")  LOCA::Extended::MultiVector::MultiVector "LOCA::Extended::MultiVector::MultiVector(const MultiVector &source,
const std::vector< int > &index, bool view)

Copy constructor that creates a sub copy or view of the given
multivector. ";

%feature("docstring")  LOCA::Extended::MultiVector::~MultiVector "LOCA::Extended::MultiVector::~MultiVector()

Vector destructor. ";


// File: classLOCA_1_1MultiContinuation_1_1NaturalConstraint.xml
%feature("docstring") LOCA::MultiContinuation::NaturalConstraint "

Implementation of LOCA::MultiContinuation::ConstraintInterface for
natural continuation.

This class implements the natural constraint equation for natural
continuation: \\\\[ g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)=
p-p_0-v_p \\\\Delta s \\\\] where $v_p$ is the parameter component of
the predictor direction $v$.

C++ includes: LOCA_MultiContinuation_NaturalConstraint.H ";

/*  Implementation of LOCA::MultiContinuation::ConstraintInterface  */

/* virtual methods

*/

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::copy "void
LOCA::MultiContinuation::NaturalConstraint::copy(const
ConstraintInterface &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::clone "Teuchos::RCP<
LOCA::MultiContinuation::ConstraintInterface >
LOCA::MultiContinuation::NaturalConstraint::clone(NOX::CopyType
type=NOX::DeepCopy) const

Cloning function. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::numConstraints "int
LOCA::MultiContinuation::NaturalConstraint::numConstraints() const

Return number of constraints. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::setX "void
LOCA::MultiContinuation::NaturalConstraint::setX(const
NOX::Abstract::Vector &y)

Set the solution vector to y. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::setParam "void
LOCA::MultiContinuation::NaturalConstraint::setParam(int paramID,
double val)

Sets parameter indexed by paramID. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::setParams "void
LOCA::MultiContinuation::NaturalConstraint::setParams(const
std::vector< int > &paramIDs, const
NOX::Abstract::MultiVector::DenseMatrix &vals)

Sets parameters indexed by paramIDs. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::computeConstraints "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeConstraints()

Compute continuation constraint equations. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::computeDX "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeDX()

Compute derivative of constraints w.r.t. solution vector x. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::computeDP "NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::NaturalConstraint::computeDP(const
std::vector< int > &paramIDs, NOX::Abstract::MultiVector::DenseMatrix
&dgdp, bool isValidG)

Compute derivative of constraints w.r.t. supplied parameters.

The first column of dgdp should be filled with the constraint
residuals $g$ if isValidG is false. If isValidG is true, then the dgdp
contains $g$ on input. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::isConstraints "bool
LOCA::MultiContinuation::NaturalConstraint::isConstraints() const

Return true if constraint residuals are valid. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::isDX "bool
LOCA::MultiContinuation::NaturalConstraint::isDX() const

Return true if derivatives of constraints w.r.t. x are valid. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::getConstraints "const
NOX::Abstract::MultiVector::DenseMatrix &
LOCA::MultiContinuation::NaturalConstraint::getConstraints() const

Return constraint residuals. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::getDX "const
NOX::Abstract::MultiVector *
LOCA::MultiContinuation::NaturalConstraint::getDX() const

Return solution component of constraint derivatives.

Since the solution component of the derivative is always zero, this
always returns NULL. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::isDXZero "bool
LOCA::MultiContinuation::NaturalConstraint::isDXZero() const

Return true if solution component of constraint derivatives is zero.
";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint "LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::MultiContinuation::NaturalGroup > &grp)

Constructor. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint "LOCA::MultiContinuation::NaturalConstraint::NaturalConstraint(const
NaturalConstraint &source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::~NaturalConstraint "LOCA::MultiContinuation::NaturalConstraint::~NaturalConstraint()

Destructor. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalConstraint::setNaturalGroup "void
LOCA::MultiContinuation::NaturalConstraint::setNaturalGroup(const
Teuchos::RCP< LOCA::MultiContinuation::NaturalGroup > &grp)

Set pointer to natural group. ";


// File: classLOCA_1_1MultiContinuation_1_1NaturalGroup.xml
%feature("docstring") LOCA::MultiContinuation::NaturalGroup "

Specialization of LOCA::MultiContinuation::ExtendedGroup to natural
continuation.

Natural continuation corresponds to a continuation equation
$g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)=0$ with $g$ given by
\\\\[ g(x,p,x_0,p_0,x^\\\\ast,p^\\\\ast,v,\\\\Delta s)= p-p_0-v_p
\\\\Delta s \\\\] where $v_p$ is the parameter component of the
predictor direction $v$. This corresponds geometrically to
constraining the nonlinear solver steps used in calculating $F(x,p)=0$
to be orthogonal to the parameter axis. The natural constraint $g$ is
represented by a LOCA::MultiContinuation::NaturalConstraint object.

C++ includes: LOCA_MultiContinuation_NaturalGroup.H ";

/*  Implementation of NOX::Abstract::Group virtual methods  */

%feature("docstring")  LOCA::MultiContinuation::NaturalGroup::clone "Teuchos::RCP< NOX::Abstract::Group >
LOCA::MultiContinuation::NaturalGroup::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. ";

/*  Implementation of LOCA::MultiContinuation::AbstractStrategy
virtual methods  */

%feature("docstring")  LOCA::MultiContinuation::NaturalGroup::copy "void LOCA::MultiContinuation::NaturalGroup::copy(const
NOX::Abstract::Group &source)

Copy. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalGroup::NaturalGroup "LOCA::MultiContinuation::NaturalGroup::NaturalGroup(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &continuationParams, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &grp, const Teuchos::RCP<
LOCA::MultiPredictor::AbstractStrategy > &pred, const std::vector< int
> &paramIDs)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list.

continuationParams:  [in] Continuation parameters.

grp:  [in] Group representing $F$.

pred:  [in] Predictor strategy.

paramIDs:  [in] Parameter IDs of continuation parameters. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalGroup::NaturalGroup "LOCA::MultiContinuation::NaturalGroup::NaturalGroup(const NaturalGroup
&source, NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")
LOCA::MultiContinuation::NaturalGroup::~NaturalGroup "LOCA::MultiContinuation::NaturalGroup::~NaturalGroup()

Destructor. ";


// File: classLOCA_1_1BorderedSolver_1_1Nested.xml
%feature("docstring") LOCA::BorderedSolver::Nested "

Bordered system solver strategy for nested bordered systems.

This class implements a bordered solver strategy for the bordered
system \\\\[ \\\\begin{bmatrix} J & A \\\\\\\\ B^T & C
\\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} =
\\\\begin{bmatrix} F \\\\\\\\ G \\\\end{bmatrix} \\\\] when $J$ itself
has this block form. It combines the blocks for $A$, $B$, and $C$ and
then instantiates a solver as specified by the \"Nested Bordered
Solver\" sublist of the solverParams pass through the constructor.
This sublist should specify the \"Bordered Solver Method\" for the
solver as well as any other parameters for that method, and any method
that can be instantiated through the LOCA::Factory is available.

Note that the operator representing $J$ must implement the
LOCA::BorderedSolver::BorderedOperator interface, and the constraint
object representing $B$ must be of type
LOCA::MultiContinuation::ConstraintInterfaceMVDX.

C++ includes: LOCA_BorderedSolver_Nested.H ";

%feature("docstring")  LOCA::BorderedSolver::Nested::Nested "LOCA::BorderedSolver::Nested::Nested(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters as described above ";

%feature("docstring")  LOCA::BorderedSolver::Nested::~Nested "LOCA::BorderedSolver::Nested::~Nested()

Destructor. ";

%feature("docstring")  LOCA::BorderedSolver::Nested::setMatrixBlocks "void LOCA::BorderedSolver::Nested::setMatrixBlocks(const Teuchos::RCP<
const LOCA::BorderedSolver::AbstractOperator > &op, const
Teuchos::RCP< const NOX::Abstract::MultiVector > &blockA, const
Teuchos::RCP< const LOCA::MultiContinuation::ConstraintInterface >
&blockB, const Teuchos::RCP< const
NOX::Abstract::MultiVector::DenseMatrix > &blockC)

Set blocks.

The blockA or blockC pointer may be null if either is zero. Whether
block B is zero will be determined by querying blockB via
ConstraintInterface::isConstraintDerivativesXZero. ";

%feature("docstring")  LOCA::BorderedSolver::Nested::initForSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::initForSolve()

Intialize solver for a solve.

This should be called after setMatrixBlocks(), but before
applyInverse(). ";

%feature("docstring")
LOCA::BorderedSolver::Nested::initForTransposeSolve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::initForTransposeSolve()

Intialize solver for a transpose solve.

This should be called after setMatrixBlocks(), but before
applyInverseTranspose(). ";

%feature("docstring")  LOCA::BorderedSolver::Nested::apply "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::apply(const NOX::Abstract::MultiVector
&X, const NOX::Abstract::MultiVector::DenseMatrix &Y,
NOX::Abstract::MultiVector &U, NOX::Abstract::MultiVector::DenseMatrix
&V) const

Computed extended matrix-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J & A \\\\\\\\ B^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J*X + A*Y \\\\\\\\ B^T*X + C*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")  LOCA::BorderedSolver::Nested::applyTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::applyTranspose(const
NOX::Abstract::MultiVector &X, const
NOX::Abstract::MultiVector::DenseMatrix &Y, NOX::Abstract::MultiVector
&U, NOX::Abstract::MultiVector::DenseMatrix &V) const

Computed extended matrix transpose-multivector product.

Computes \\\\[ \\\\begin{bmatrix} U \\\\\\\\ V \\\\end{bmatrix} =
\\\\begin{bmatrix} J^T & B \\\\\\\\ A^T & C \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
J^T*X + B*Y \\\\\\\\ A^T*X + C^T*Y \\\\end{bmatrix}. \\\\] ";

%feature("docstring")  LOCA::BorderedSolver::Nested::applyInverse "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::applyInverse(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system as defined above using bordering.

The params argument is the linear solver parameters. If isZeroF or
isZeroG is true, than the corresponding F or G pointers may be NULL.
";

%feature("docstring")
LOCA::BorderedSolver::Nested::applyInverseTranspose "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::Nested::applyInverseTranspose(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the transpose of the extended system as defined above using
bordering.

The params argument is the linear solver parameters. If isZeroF or
isZeroG is true, than the corresponding F or G pointers may be NULL.
";


// File: classLOCA_1_1SingularJacobianSolve_1_1Nic.xml
%feature("docstring") LOCA::SingularJacobianSolve::Nic "

This class computes the solution to $J x = b$ using the Nic method.

The idea here is to use deflation of the right hand side to improve
the conditioning of the linear system. Typically a solution to $J x =
b$ when $J$ is nearly singular will have a large component in the
direction of the null vector $v$. The idea then is to deflate $Jv$ out
of the right hand side $b$. The complete algorithm used here is: \\\\[
\\\\begin{aligned} &\\\\tilde{b} = b - \\\\frac{v^T b}{v^T J v} Jv
\\\\\\\\ &\\\\text{Solve}\\\\; J\\\\tilde{x} = \\\\tilde{b} \\\\\\\\
&x = \\\\tilde{x} + \\\\frac{v^T b}{v^T J v} v \\\\end{aligned} \\\\]
The solve $J\\\\tilde{x} = \\\\tilde{b}$ uses the underlying group's
applyJacobianInverse() method and therefore this is a generic
technique for computing solutions to nearly singular system since it
uses any supplied linear solver.

This algorithm is selected by setting the \"Method\" parameter of the
\"Singular Solve\" sublist of the NOX linear solver parameter list to
\"Nic\". The idea for this algorithm is taken from: R. A. Nicolaides,
\"Deflation of Conjugate   Gradients With Applications to Boundary
Value Problems,\" SIAM J. Numer. Anal., 24(2), 1987.

C++ includes: LOCA_SingularJacobianSolve_Nic.H ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::Nic "LOCA::SingularJacobianSolve::Nic::Nic(Teuchos::ParameterList &params)

Constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::Nic "LOCA::SingularJacobianSolve::Nic::Nic(const Nic &source)

Copy constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::~Nic "LOCA::SingularJacobianSolve::Nic::~Nic()

Destructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::clone "LOCA::SingularJacobianSolve::Generic *
LOCA::SingularJacobianSolve::Nic::clone() const

Clone function. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::reset "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Nic::reset(Teuchos::ParameterList
&params)

Reset parameters.

There are no additional parameters for the Nic calculation. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::compute "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Nic::compute(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector &input, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector &result)

Computes the solution as described above. ";

%feature("docstring")  LOCA::SingularJacobianSolve::Nic::computeMulti
"NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::Nic::computeMulti(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector *const *inputs, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector **results, int nVecs)

Computes solution for multiple RHS. ";


// File: classLOCA_1_1SingularJacobianSolve_1_1NicDay.xml
%feature("docstring") LOCA::SingularJacobianSolve::NicDay "

This class computes the solution to $J x = b$ using the Nic-Day
method.

This singular solve method is a modification of the deflation idea
implemented in LOCA::SingularJacobianSolve::Nic where deflation of the
right hand side is used to improve the conditioning of the linear
system. Typically a solution to $J x = b$ when $J$ is nearly singular
will have a large component in the direction of the null vector $v$.
The idea then is to deflate $Jv$ out of the right hand side $b$. The
complete algorithm used here is: \\\\[ \\\\begin{aligned}
&\\\\tilde{b} = b - \\\\frac{b^T J v}{v^T J^T J v} Jv \\\\\\\\
&\\\\text{Solve}\\\\; J\\\\tilde{x} = \\\\tilde{b} \\\\\\\\ &x =
\\\\tilde{x} + \\\\frac{b^T J v}{v^T J^T J v} v \\\\end{aligned} \\\\]
The solve $J\\\\tilde{x} = \\\\tilde{b}$ uses the underlying group's
applyJacobianInverse() method and therefore this is a generic
technique for computing solutions to nearly singular system since it
uses any supplied linear solver.

This algorithm is selected by setting the \"Method\" parameter of the
\"Singular Solve\" sublist of the NOX linear solver parameter list to
\"Nic-Day\". The idea for this algorithm is taken from: R. A.
Nicolaides, \"Deflation of Conjugate   Gradients With Applications to
Boundary Value Problems,\" SIAM J. Numer. Anal., 24(2), 1987.

C++ includes: LOCA_SingularJacobianSolve_NicDay.H ";

%feature("docstring")  LOCA::SingularJacobianSolve::NicDay::NicDay "LOCA::SingularJacobianSolve::NicDay::NicDay(Teuchos::ParameterList
&params)

Constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::NicDay::NicDay "LOCA::SingularJacobianSolve::NicDay::NicDay(const NicDay &source)

Copy constructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::NicDay::~NicDay "LOCA::SingularJacobianSolve::NicDay::~NicDay()

Destructor. ";

%feature("docstring")  LOCA::SingularJacobianSolve::NicDay::clone "LOCA::SingularJacobianSolve::Generic *
LOCA::SingularJacobianSolve::NicDay::clone() const

Clone function. ";

%feature("docstring")  LOCA::SingularJacobianSolve::NicDay::reset "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::NicDay::reset(Teuchos::ParameterList
&params)

Reset parameters.

There are no additional parameters for the NicDay calculation. ";

%feature("docstring")  LOCA::SingularJacobianSolve::NicDay::compute "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::NicDay::compute(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector &input, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector &result)

Computes the solution as described above. ";

%feature("docstring")
LOCA::SingularJacobianSolve::NicDay::computeMulti "NOX::Abstract::Group::ReturnType
LOCA::SingularJacobianSolve::NicDay::computeMulti(Teuchos::ParameterList
&params, LOCA::Continuation::AbstractGroup &grp, const
NOX::Abstract::Vector *const *inputs, const NOX::Abstract::Vector
&approxNullVec, const NOX::Abstract::Vector &jacApproxNullVec,
NOX::Abstract::Vector **results, int nVecs)

Computes solution for multiple RHS. ";


// File: classLOCA_1_1Bifurcation_1_1TPBord_1_1StatusTest_1_1NullVectorNormWRMS.xml
%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS "

A convergence test based on the weighted root-mean-square norm of the
update to the null vector component for turning point location.

Let $n$ be the approximation to the null vector for turning point
tracking (see LOCA::Bifurcation::TPBord::ExtendedGroup). This
convergence test defines convergence for the null vector when the
following is true \\\\[
\\\\sqrt{\\\\frac{1}{N}\\\\sum_{i=1}^N\\\\left(\\\\frac{|n_i-(n_0)_i|}{\\\\epsilon_r|n_i|
+ \\\\epsilon_a}\\\\right)} < \\\\tau \\\\] where $n_0$ is the
previous approximation to the null vector, $N$ is the length of $n$,
$\\\\epsilon_r$ is the relative tolerance, $\\\\epsilon_a$ is the
absolute tolerance, and $\\\\tau$ is an overall scale factor
(typically $\\\\tau = 1$).

Note that this status test deals only with the null vector component
of the turning point equations. This status test should be combined
with other status tests for the solution and parameter components
(using NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build
a composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a turning point group
(i.e., not derived from LOCA::Bifurcation::TPBord::ExtendedGroup),
checkStatus() returns NOX::StatusTest::Converged. This allows the
status test to be used in situations other than turning point
tracking, e.g., steady- state solves, without raising error
conditions.

C++ includes: LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H
";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getNullVectorNormWRMS
"double
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getNullVectorNormWRMS()
const

Returns the value of weighted parameter update norm. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getRTOL "double
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getRTOL()
const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getATOL "double
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getATOL()
const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getTOL "double
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getTOL()
const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::NullVectorNormWRMS
"LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::NullVectorNormWRMS(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::~NullVectorNormWRMS
"LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::~NullVectorNormWRMS()

Destructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::checkStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getStatus "NOX::StatusTest::StatusType
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getStatus()
const

Returns status as defined above. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::print "virtual std::ostream&
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1Bifurcation_1_1PitchforkBord_1_1StatusTest_1_1NullVectorNormWRMS.xml
%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS "

A convergence test based on the weighted root-mean-square norm of the
update to the null vector component for pitchfork location.

Let $n$ be the approximation to the null vector for pitchfork tracking
(see LOCA::Bifurcation::PitchforkBord::ExtendedGroup). This
convergence test defines convergence for the null vector when the
following is true \\\\[
\\\\sqrt{\\\\frac{1}{N}\\\\sum_{i=1}^N\\\\left(\\\\frac{|n_i-(n_0)_i|}{\\\\epsilon_r|n_i|
+ \\\\epsilon_a}\\\\right)} < \\\\tau \\\\] where $n_0$ is the
previous approximation to the null vector, $N$ is the length of $n$,
$\\\\epsilon_r$ is the relative tolerance, $\\\\epsilon_a$ is the
absolute tolerance, and $\\\\tau$ is an overall scale factor
(typically $\\\\tau = 1$).

Note that this status test deals only with the null vector component
of the pitchfork equations. This status test should be combined with
other status tests for the solution and parameter components (using
NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build a
composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a pitchfork group
(i.e., not derived from
LOCA::Bifurcation::PitchforkBord::ExtendedGroup), checkStatus()
returns NOX::StatusTest::Converged. This allows the status test to be
used in situations other than pitchfork tracking, e.g., steady-state
solves, without raising error conditions.

C++ includes: LOCA_Bifurcation_PitchforkBord_NullVectorNormWRMS.H ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getNullVectorNormWRMS
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getNullVectorNormWRMS()
const

Returns the value of weighted parameter update norm. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getRTOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getRTOL()
const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getATOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getATOL()
const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getTOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getTOL()
const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::NullVectorNormWRMS
"LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::NullVectorNormWRMS(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::~NullVectorNormWRMS
"LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::~NullVectorNormWRMS()

Destructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::checkStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::getStatus()
const

Returns status as defined above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::print
"virtual std::ostream&
LOCA::Bifurcation::PitchforkBord::StatusTest::NullVectorNormWRMS::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1Continuation_1_1StatusTest_1_1ParameterResidualNorm.xml
%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm "

A convergence test based on the parameter component of the residual
for continuation.

Consider a continuation method with parameter equation $g = 0$ (see
LOCA::Continuation::ExtendedGroup). This convergence test defines
convergence of the parameter equation when the following is true \\\\[
\\\\frac{|g|}{\\\\epsilon_r|\\\\Delta s| + \\\\epsilon_a} < \\\\tau
\\\\] where $\\\\Delta s$ is the current step size, $\\\\epsilon_r$ is
the relative tolerance, $\\\\epsilon_a$ is the absolute tolerance, and
$\\\\tau$ is an overall scale factor (typically $\\\\tau = 1$).

Note that this status test deals only with the parameter component of
the continuation equations. This status test should be combined with
other status tests for the solution component (using
NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build a
composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a continuation group
(i.e., not derived from LOCA::Continuation::ExtendedGroup),
checkStatus() returns NOX::StatusTest::Converged. This allows the
status test to be used in situations other than continuation, e.g.,
steady-state solves, without raising error conditions.

C++ includes: LOCA_Continuation_StatusTest_ParameterResidualNorm.H ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::getResidualNorm
"double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getResidualNorm()
const

Returns the value of scaled parameter residual norm. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::getRTOL "double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getRTOL() const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::getATOL "double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getATOL() const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::getTOL "double
LOCA::Continuation::StatusTest::ParameterResidualNorm::getTOL() const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::ParameterResidualNorm
"LOCA::Continuation::StatusTest::ParameterResidualNorm::ParameterResidualNorm(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::~ParameterResidualNorm
"LOCA::Continuation::StatusTest::ParameterResidualNorm::~ParameterResidualNorm()

Destructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::checkStatus "NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterResidualNorm::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::getStatus "NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterResidualNorm::getStatus()
const

Returns status as defined above. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterResidualNorm::print "ostream
& LOCA::Continuation::StatusTest::ParameterResidualNorm::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1Bifurcation_1_1TPBord_1_1StatusTest_1_1ParameterUpdateNorm.xml
%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm "

A convergence test based on the update of the parameter component for
turning point location.

Let $p$ be the turning point parameter (see
LOCA::Bifurcation::TPBord::ExtendedGroup). This convergence test
defines convergence for the parameter when the following is true \\\\[
\\\\frac{|p-p_0|}{\\\\epsilon_r|p| + \\\\epsilon_a} < \\\\tau \\\\]
where $p_0$ is the previous parameter value, $\\\\epsilon_r$ is the
relative tolerance, $\\\\epsilon_a$ is the absolute tolerance, and
$\\\\tau$ is an overall scale factor (typically $\\\\tau = 1$).

Note that this status test deals only with the parameter component of
the turning point equations. This status test should be combined with
other status tests for the solution and null vector components (using
NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build a
composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a turning point group
(i.e., not derived from LOCA::Bifurcation::TPBord::ExtendedGroup),
checkStatus() returns NOX::StatusTest::Converged. This allows the
status test to be used in situations other than turning point
tracking, e.g., steady- state solves, without raising error
conditions.

C++ includes: LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H
";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getUpdateNorm
"double
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getUpdateNorm()
const

Returns the value of weighted parameter update norm. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getRTOL "double
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getRTOL()
const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getATOL "double
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getATOL()
const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getTOL "double
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getTOL()
const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm
"LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm
"LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm()

Destructor. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::checkStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getStatus()
const

Returns status as defined above. ";

%feature("docstring")
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::print "virtual std::ostream&
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1Bifurcation_1_1PitchforkBord_1_1StatusTest_1_1ParameterUpdateNorm.xml
%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm "

A convergence test based on the update of the parameter component for
pitchfork location.

Let $p$ be the pitchfork parameter (see
LOCA::Bifurcation::PitchforkBord::ExtendedGroup). This convergence
test defines convergence for the parameter when the following is true
\\\\[ \\\\frac{|p-p_0|}{\\\\epsilon_r|p| + \\\\epsilon_a} < \\\\tau
\\\\] where $p_0$ is the previous parameter value, $\\\\epsilon_r$ is
the relative tolerance, $\\\\epsilon_a$ is the absolute tolerance, and
$\\\\tau$ is an overall scale factor (typically $\\\\tau = 1$).

Note that this status test deals only with the parameter component of
the pitchfork equations. This status test should be combined with
other status tests for the solution and null vector components (using
NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build a
composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a pitchfork group
(i.e., not derived from
LOCA::Bifurcation::Pitchforkbord::ExtendedGroup), checkStatus()
returns NOX::StatusTest::Converged. This allows the status test to be
used in situations other than pitchfork tracking, e.g., steady-state
solves, without raising error conditions.

C++ includes: LOCA_Bifurcation_PitchforkBord_ParameterUpdateNorm.H ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getUpdateNorm
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getUpdateNorm()
const

Returns the value of weighted parameter update norm. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getRTOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getRTOL()
const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getATOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getATOL()
const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getTOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getTOL()
const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm
"LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm
"LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm()

Destructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::checkStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::getStatus()
const

Returns status as defined above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::print
"virtual std::ostream&
LOCA::Bifurcation::PitchforkBord::StatusTest::ParameterUpdateNorm::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1Continuation_1_1StatusTest_1_1ParameterUpdateNorm.xml
%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm "

A convergence test based on the update of the parameter component for
continuation.

Consider a continuation method with parameter equation $g = 0$ (see
LOCA::Continuation::ExtendedGroup). This convergence test defines
convergence of the parameter equation when the following is true \\\\[
\\\\frac{|p-p_0|}{\\\\epsilon_r|p_0| + \\\\epsilon_a} < \\\\tau \\\\]
where $p$ is the current paramter value, $p_0$ is the previous
parameter value, $\\\\epsilon_r$ is the relative tolerance,
$\\\\epsilon_a$ is the absolute tolerance, and $\\\\tau$ is an overall
scale factor (typically $\\\\tau = 1$).

Note that this status test deals only with the parameter component of
the continuation equations. This status test should be combined with
other status tests for the solution component (using
NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build a
composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a continuation group
(i.e., not derived from LOCA::Continuation::ExtendedGroup),
checkStatus() returns NOX::StatusTest::Converged. This allows the
status test to be used in situations other than continuation, e.g.,
steady-state solves, without raising error conditions.

C++ includes: LOCA_Continuation_StatusTest_ParameterUpdateNorm.H ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getUpdateNorm "double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getUpdateNorm()
const

Returns the value of weighted parameter update norm. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getRTOL "double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getRTOL() const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getATOL "double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getATOL() const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getTOL "double
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getTOL() const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm
"LOCA::Continuation::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm
"LOCA::Continuation::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm()

Destructor. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::checkStatus "NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterUpdateNorm::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getStatus "NOX::StatusTest::StatusType
LOCA::Continuation::StatusTest::ParameterUpdateNorm::getStatus() const

Returns status as defined above. ";

%feature("docstring")
LOCA::Continuation::StatusTest::ParameterUpdateNorm::print "ostream &
LOCA::Continuation::StatusTest::ParameterUpdateNorm::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1ParameterVector.xml
%feature("docstring") LOCA::ParameterVector "

LOCA's container for holding a set of parameters that are used by the
LOCA continuation routines.

Roger Pawlowski (SNL 9233)

C++ includes: LOCA_Parameter_Vector.H ";

%feature("docstring")  LOCA::ParameterVector::ParameterVector "LOCA::ParameterVector::ParameterVector()

Constructor. ";

%feature("docstring")  LOCA::ParameterVector::ParameterVector "LOCA::ParameterVector::ParameterVector(const ParameterVector &source)

Copy constructor. ";

%feature("docstring")  LOCA::ParameterVector::clone "LOCA::ParameterVector * LOCA::ParameterVector::clone() const

Clone. ";

%feature("docstring")  LOCA::ParameterVector::~ParameterVector "LOCA::ParameterVector::~ParameterVector()

Destructor. ";

%feature("docstring")  LOCA::ParameterVector::addParameter "int
LOCA::ParameterVector::addParameter(std::string label, double
value=0.0)

Adds a parameter to the list. Returns the index value assigned to the
parameter. ";

%feature("docstring")  LOCA::ParameterVector::init "bool
LOCA::ParameterVector::init(double value)

Initialize the vector. Returns true if successful. ";

%feature("docstring")  LOCA::ParameterVector::scale "bool
LOCA::ParameterVector::scale(double value)

Scales the entire vector by value. Returns true if successful. ";

%feature("docstring")  LOCA::ParameterVector::scale "bool
LOCA::ParameterVector::scale(const ParameterVector &p)

Scales the vactor with another vector (element-wise multiply). Returns
true if successful. ";

%feature("docstring")  LOCA::ParameterVector::update "bool
LOCA::ParameterVector::update(double alpha, const ParameterVector
&alphaVector, double b)

Updates the parameter vector: this = alpha * alphaVector + b * this.
Returns true if successful. ";

%feature("docstring")  LOCA::ParameterVector::setValue "void
LOCA::ParameterVector::setValue(unsigned int i, double value)

Set the value of the parameter with index i. Will throw an error if
index is out of range. ";

%feature("docstring")  LOCA::ParameterVector::setValue "void
LOCA::ParameterVector::setValue(std::string label, double value)

Set the value of the parameter with the corresponding label. Will
throw an error if \"label\" is not valid. ";

%feature("docstring")  LOCA::ParameterVector::getValue "double
LOCA::ParameterVector::getValue(unsigned int i) const

Returns the value of the parameter with index i. Will throw an error
if index is out of range. ";

%feature("docstring")  LOCA::ParameterVector::getValue "double
LOCA::ParameterVector::getValue(std::string label) const

Returns the value of the parameter with the corresponding label. Will
throw an error if \"label\" is not valid. ";

%feature("docstring")  LOCA::ParameterVector::getIndex "int
LOCA::ParameterVector::getIndex(std::string label) const

Returns the index of the parameter with the corresponding label.
Returns a -1 if \"label\" is not found. ";

%feature("docstring")  LOCA::ParameterVector::getDoubleArrayPointer "double * LOCA::ParameterVector::getDoubleArrayPointer()

Returns a pointer to a C-style array of the parameter values. ";

%feature("docstring")  LOCA::ParameterVector::isParameter "bool
LOCA::ParameterVector::isParameter(std::string label) const

Returns true if the parameter std::string \"label\" corresponds to a
parameter label in the object. ";

%feature("docstring")  LOCA::ParameterVector::getLabel "std::string
LOCA::ParameterVector::getLabel(unsigned int i) const

Returns the label of the parameter with index i. ";

%feature("docstring")  LOCA::ParameterVector::length "int
LOCA::ParameterVector::length() const

Returns the length of parameter vector. ";

%feature("docstring")  LOCA::ParameterVector::print "void
LOCA::ParameterVector::print(std::ostream &stream) const

Prints the vector to cout. ";

%feature("docstring")  LOCA::ParameterVector::getValuesVector "const
std::vector< double > & LOCA::ParameterVector::getValuesVector() const

Accessor to get the underlying stl vector with all parameter values.
";

%feature("docstring")  LOCA::ParameterVector::getNamesVector "const
std::vector< std::string > & LOCA::ParameterVector::getNamesVector()
const

Accessor to get the underlying stl vector with all parameter names. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1PhippsBordering.xml
%feature("docstring") LOCA::TurningPoint::MooreSpence::PhippsBordering
"

Moore-Spence turning point solver strategy based on \"Phipps\"
bordering which is the 5-solve modified turning point bordering
algorithm that uses bordered linear solves.

This class solves the Moore-Spence turning point Newton equations:
\\\\[ \\\\begin{bmatrix} J & 0 & f_p \\\\\\\\ (Jv)_x & J & (Jv)_p
\\\\\\\\ 0 & \\\\phi^T & 0 \\\\end{bmatrix} \\\\begin{bmatrix} X
\\\\\\\\ Y \\\\\\\\ z \\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\
G \\\\\\\\ h \\\\end{bmatrix} \\\\] via the following modified block
elimination scheme: \\\\[ \\\\begin{split} \\\\begin{bmatrix} J & u
\\\\\\\\ v^T & 0 \\\\end{bmatrix} \\\\begin{bmatrix} A & B \\\\\\\\ a
& b \\\\end{bmatrix} &= \\\\begin{bmatrix} F & f_p \\\\\\\\ 0 & 0
\\\\end{bmatrix} \\\\\\\\ \\\\begin{bmatrix} J & u \\\\\\\\ v^T & 0
\\\\end{bmatrix} \\\\begin{bmatrix} C & D & E\\\\\\\\ c & d & e
\\\\end{bmatrix} &= \\\\begin{bmatrix} (Jv)_x A - G & (Jv)_x B -
(Jv)_p & (Jv)_x v \\\\\\\\ 0 & 0 & 0 \\\\end{bmatrix} \\\\\\\\
\\\\begin{bmatrix} \\\\sigma & 0 & b \\\\\\\\ e & \\\\sigma & -d
\\\\\\\\ -\\\\phi^T E & \\\\phi^T v & \\\\phi^T D \\\\end{bmatrix}
\\\\begin{bmatrix} \\\\alpha \\\\\\\\ \\\\beta \\\\\\\\ z
\\\\end{bmatrix} &= \\\\begin{bmatrix} a \\\\\\\\ c \\\\\\\\ h +
\\\\phi^T C \\\\end{bmatrix} \\\\\\\\ X &= A - B z + v \\\\alpha
\\\\\\\\ Y &= -C + d z - E \\\\alpha + v \\\\beta \\\\end{split} \\\\]
where $s = \\\\|J v\\\\|$ and $u = J v/s$. Each bordered solve is
implemented by a LOCA::BorderedSolver::AbstractStrategy strategy
object.

C++ includes: LOCA_TurningPoint_MooreSpence_PhippsBordering.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::PhippsBordering::PhippsBordering "LOCA::TurningPoint::MooreSpence::PhippsBordering::PhippsBordering(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Instantiates a
bordered solver for solving the bordeded systems described above. See
LOCA::BorderedSolver::Factory for a description of available solvers.
";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::PhippsBordering::~PhippsBordering "LOCA::TurningPoint::MooreSpence::PhippsBordering::~PhippsBordering()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::PhippsBordering::setBlocks "void
LOCA::TurningPoint::MooreSpence::PhippsBordering::setBlocks(const
Teuchos::RCP< LOCA::TurningPoint::MooreSpence::AbstractGroup > &group,
const Teuchos::RCP< LOCA::TurningPoint::MooreSpence::ExtendedGroup >
&tpGroup, const Teuchos::RCP< const NOX::Abstract::Vector >
&nullVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&JnVector, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&dfdp, const Teuchos::RCP< const NOX::Abstract::MultiVector > &dJndp)

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

tpGroup:  [in] Turning point group representing the turning point
equations.

nullVector:  [in] Vector representing v

JnVector:  [in] Vector representing Jv

dfdp:  [in] Vector representing df/dp

dJndp:  [in] Vector representing d(Jv)/dp ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::PhippsBordering::solve "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::PhippsBordering::solve(Teuchos::ParameterList
&params, const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector
&input, LOCA::TurningPoint::MooreSpence::ExtendedMultiVector &result)
const

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::PhippsBordering::solveTranspose "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::PhippsBordering::solveTranspose(Teuchos::ParameterList
&params, const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector
&input, LOCA::TurningPoint::MooreSpence::ExtendedMultiVector &result)
const

Solves the transpose of the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1PhippsBordering.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::PhippsBordering "

Moore-Spence pitchfork solver strategy based on \"Phipps\" bordering
which is the 7-solve modified pitchfork bordering algorithm that uses
bordered linear solves.

This class solves the Moore-Spence pitchfork Newton equations: \\\\[
\\\\begin{bmatrix} J & 0 & \\\\psi & f_p \\\\\\\\ (Jv)_x & J & 0 &
(Jv)_p \\\\\\\\ \\\\psi^T & 0 & 0 & 0 \\\\\\\\ 0 & \\\\phi^T & 0 & 0
\\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\\\\\ w \\\\\\\\ z
\\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\ G \\\\\\\\ s \\\\\\\\
h \\\\end{bmatrix} \\\\] via the following modified block elimination
scheme: \\\\[ \\\\begin{split} \\\\begin{bmatrix} J & u \\\\\\\\ v^T &
0 \\\\end{bmatrix} \\\\begin{bmatrix} A & B & C \\\\\\\\ a & b & c
\\\\end{bmatrix} &= \\\\begin{bmatrix} F & f_p & \\\\psi \\\\\\\\ 0 &
0 & 0 \\\\end{bmatrix} \\\\\\\\ \\\\begin{bmatrix} J & u \\\\\\\\ v^T
& 0 \\\\end{bmatrix} \\\\begin{bmatrix} D & E & K & L \\\\\\\\ d & e &
k & l \\\\end{bmatrix} &= \\\\begin{bmatrix} G - (Jv)_x A & (Jv)_p -
(Jv)_x B & -(Jv)_x C & -(Jv)_x v \\\\\\\\ 0 & 0 & 0 & 0
\\\\end{bmatrix} \\\\\\\\ \\\\begin{bmatrix} \\\\sigma & 0 & b & c
\\\\\\\\ -l & \\\\sigma & e & k \\\\\\\\ \\\\langle
v,\\\\psi\\\\rangle & 0 & -\\\\langle B,\\\\psi\\\\rangle &
-\\\\langle C,\\\\psi\\\\rangle \\\\\\\\ \\\\phi^T L & \\\\phi^T v &
-\\\\phi^T E & -\\\\phi^T K \\\\end{bmatrix} \\\\begin{bmatrix}
\\\\alpha \\\\\\\\ \\\\beta \\\\\\\\ z \\\\\\\\ w \\\\end{bmatrix} &=
\\\\begin{bmatrix} a \\\\\\\\ d \\\\\\\\ s - \\\\langle
A,\\\\psi\\\\rangle h - \\\\phi^T D \\\\end{bmatrix} \\\\\\\\ X &= A -
B z - C w + v \\\\alpha \\\\\\\\ Y &= D - E z - K w + L \\\\alpha + v
\\\\beta \\\\end{split} \\\\] where $\\\\sigma = \\\\|J v\\\\|$ and $u
= J v/\\\\sigma$. Each bordered solve is implemented by a
LOCA::BorderedSolver::AbstractStrategy strategy object.

C++ includes: LOCA_Pitchfork_MooreSpence_PhippsBordering.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::PhippsBordering::PhippsBordering "LOCA::Pitchfork::MooreSpence::PhippsBordering::PhippsBordering(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Instantiates a
bordered solver for solving the bordeded systems described above. See
LOCA::BorderedSolver::Factory for a description of available solvers.
";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::PhippsBordering::~PhippsBordering "LOCA::Pitchfork::MooreSpence::PhippsBordering::~PhippsBordering()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::PhippsBordering::setBlocks "void
LOCA::Pitchfork::MooreSpence::PhippsBordering::setBlocks(const
Teuchos::RCP< LOCA::Pitchfork::MooreSpence::AbstractGroup > &group,
const Teuchos::RCP< LOCA::Pitchfork::MooreSpence::ExtendedGroup >
&pfGroup, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&asymMultiVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&nullVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&JnVector, const Teuchos::RCP< const NOX::Abstract::Vector > &dfdp,
const Teuchos::RCP< const NOX::Abstract::Vector > &dJndp)

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

pfGroup:  [in] Pitchfork group representing the pitchfork equations.

asymMultiVector:  [in] Multivector representing the asymmetric vector

nullVector:  [in] Vector representing v

JnVector:  [in] Vector representing Jv

dfdp:  [in] Vector representing df/dp

dJndp:  [in] Vector representing d(Jv)/dp ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::PhippsBordering::solve "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::PhippsBordering::solve(Teuchos::ParameterList
&params, const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector
&input, LOCA::Pitchfork::MooreSpence::ExtendedMultiVector &result)
const

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1MultiPredictor_1_1Random.xml
%feature("docstring") LOCA::MultiPredictor::Random "

Random predictor strategy

This class computes the predictor direction where the solution
component is filled with random values and the parameter component
equal to 1. Each componenet of the solution vector $v_i$ of the
predictor is given by $v_i = \\\\epsilon r_i x_i$ where $r_i$ is a
random value between -1 and 1, $x_i$ is the corresponding component of
the solution vector, and $\\\\epsilon$ is a parameter.

The parameters used by this class supplied in the constructor are:
\"Epsilon\" - $\\\\epsilon$ as defined above (Default 1.0e-3)

C++ includes: LOCA_MultiPredictor_Random.H ";

%feature("docstring")  LOCA::MultiPredictor::Random::Random "LOCA::MultiPredictor::Random::Random(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &predParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

predParams:  [in] Predictor parameters as described above. ";

%feature("docstring")  LOCA::MultiPredictor::Random::~Random "LOCA::MultiPredictor::Random::~Random()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::Random::Random "LOCA::MultiPredictor::Random::Random(const Random &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::MultiPredictor::Random::clone "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::MultiPredictor::Random::clone(NOX::CopyType type=NOX::DeepCopy)
const

Clone function. ";

%feature("docstring")  LOCA::MultiPredictor::Random::compute "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Random::compute(bool baseOnSecant, const
std::vector< double > &stepSize,
LOCA::MultiContinuation::ExtendedGroup &grp, const
LOCA::MultiContinuation::ExtendedVector &prevXVec, const
LOCA::MultiContinuation::ExtendedVector &xVec)

Compute the predictor given the current and previous solution vectors.
Set baseOnSecant to false if the predictor orientation should not be
based on the secant vector (first or last steps of a continuation
run).

This method actually implements the predictor computation described
above ";

%feature("docstring")  LOCA::MultiPredictor::Random::evaluate "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Random::evaluate(const std::vector< double >
&stepSize, const LOCA::MultiContinuation::ExtendedVector &xVec,
LOCA::MultiContinuation::ExtendedMultiVector &result) const

Evaluate predictor with step size stepSize.

This method computes result[i] = xVec[i] + stepSize[i] * v[i] for each
i, where v[i] is the ith predictor direction. ";

%feature("docstring")  LOCA::MultiPredictor::Random::computeTangent "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Random::computeTangent(LOCA::MultiContinuation::ExtendedMultiVector
&tangent)

Compute tangent to predictor and store in tangent. ";

%feature("docstring")  LOCA::MultiPredictor::Random::isTangentScalable
"bool LOCA::MultiPredictor::Random::isTangentScalable() const

Is the tangent vector for this predictor scalable.

For the random predictor, this always returns false. ";


// File: classLOCA_1_1Epetra_1_1Interface_1_1Required.xml
%feature("docstring") LOCA::Epetra::Interface::Required "

Used by LOCA::Epetra::Group to provide a link to the external code for
setting problem parameters.

This interface is derived from the NOX::Epetra::Interface::Required
and additionally provides a method for setting problem parameters.

C++ includes: LOCA_Epetra_Interface_Required.H ";

%feature("docstring")  LOCA::Epetra::Interface::Required::Required "LOCA::Epetra::Interface::Required::Required()

Constructor. ";

%feature("docstring")  LOCA::Epetra::Interface::Required::~Required "virtual LOCA::Epetra::Interface::Required::~Required()

Destructor. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::setParameters "virtual void
LOCA::Epetra::Interface::Required::setParameters(const ParameterVector
&p)=0

Set parameters in the user's application.

Should be called prior to calling one of the compute functions. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::printSolution "virtual void
LOCA::Epetra::Interface::Required::printSolution(const Epetra_Vector
&x_, double conParam)

Call user's own print routine for vector-parameter pair. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::dataForPrintSolution "virtual void
LOCA::Epetra::Interface::Required::dataForPrintSolution(const int
conStep, const int timeStep, const int totalTimeSteps)

Provides data to application for output files.

This routine is called from Interface::xyzt::printSolution() just
before the call to Interface::Required::printSolution(x,param), and
gives the application some indices that can be used for creating a
unique name/index for the output files. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::setMultiPointParameter "virtual
void LOCA::Epetra::Interface::Required::setMultiPointParameter(const
int stepNum)

Set multipoint parameter in the user's application.

Should be called prior to calling one of the compute functions. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::preProcessContinuationStep "virtual void
LOCA::Epetra::Interface::Required::preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus, LOCA::Epetra::Group &group)

Perform any preprocessing before a continuation step starts.

The stepStatus argument indicates whether the previous step was
successful. The default implementation here is empty. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::postProcessContinuationStep "virtual void
LOCA::Epetra::Interface::Required::postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus
stepStatus, LOCA::Epetra::Group &group)

Perform any postprocessing after a continuation step finishes.

The stepStatus argument indicates whether the step was successful. The
default implementation here is empty. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::projectToDraw "virtual void
LOCA::Epetra::Interface::Required::projectToDraw(const
NOX::Epetra::Vector &x, double *px) const

Projects solution to a few scalars for multiparameter continuation.

Default implementation is the max norm. ";

%feature("docstring")
LOCA::Epetra::Interface::Required::projectToDrawDimension "virtual
int LOCA::Epetra::Interface::Required::projectToDrawDimension() const

Returns the dimension of the projection to draw array. ";


// File: classLOCA_1_1MultiPredictor_1_1Restart.xml
%feature("docstring") LOCA::MultiPredictor::Restart "

Restart predictor strategy

This class implements a predictor that is restarted from a previous
computation. In other words, this class takes a predictor vector that
would be computed previously and uses it as the predictor.

The parameters used by this class supplied in the constructor are:
\"Restart Vector\" - Teuchos::RCP to a
LOCA::MultiContinuation::ExtendedVector or ExtendedMultiVector storing
the predictor direction.

C++ includes: LOCA_MultiPredictor_Restart.H ";

%feature("docstring")  LOCA::MultiPredictor::Restart::Restart "LOCA::MultiPredictor::Restart::Restart(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &predParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

predParams:  [in] Predictor parameters as described above. ";

%feature("docstring")  LOCA::MultiPredictor::Restart::~Restart "LOCA::MultiPredictor::Restart::~Restart()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::Restart::Restart "LOCA::MultiPredictor::Restart::Restart(const Restart &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::MultiPredictor::Restart::clone "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::MultiPredictor::Restart::clone(NOX::CopyType type=NOX::DeepCopy)
const

Clone function. ";

%feature("docstring")  LOCA::MultiPredictor::Restart::compute "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Restart::compute(bool baseOnSecant, const
std::vector< double > &stepSize,
LOCA::MultiContinuation::ExtendedGroup &grp, const
LOCA::MultiContinuation::ExtendedVector &prevXVec, const
LOCA::MultiContinuation::ExtendedVector &xVec)

Compute the predictor given the current and previous solution vectors.
Set baseOnSecant to false if the predictor orientation should not be
based on the secant vector (first or last steps of a continuation
run).

This method actually implements the predictor computation described
above ";

%feature("docstring")  LOCA::MultiPredictor::Restart::evaluate "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Restart::evaluate(const std::vector< double >
&stepSize, const LOCA::MultiContinuation::ExtendedVector &xVec,
LOCA::MultiContinuation::ExtendedMultiVector &result) const

Evaluate predictor with step size stepSize.

This method computes result[i] = xVec[i] + stepSize[i] * v[i] for each
i, where v[i] is the ith predictor direction. ";

%feature("docstring")  LOCA::MultiPredictor::Restart::computeTangent "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Restart::computeTangent(LOCA::MultiContinuation::ExtendedMultiVector
&tangent)

Compute tangent to predictor and store in tangent. ";

%feature("docstring")
LOCA::MultiPredictor::Restart::isTangentScalable "bool
LOCA::MultiPredictor::Restart::isTangentScalable() const

Is the tangent vector for this predictor scalable.

For the restart predictor, this always returns false. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1SalingerBordering.xml
%feature("docstring")
LOCA::TurningPoint::MooreSpence::SalingerBordering "

Moore-Spence turning point solver strategy based on \"Salinger\"
bordering. This is the classic 4-solve bordering method.

This class solves the Moore-Spence turning point Newton equations:
\\\\[ \\\\begin{bmatrix} J & 0 & f_p \\\\\\\\ (Jv)_x & J & (Jv)_p
\\\\\\\\ 0 & \\\\phi^T & 0 \\\\end{bmatrix} \\\\begin{bmatrix} X
\\\\\\\\ Y \\\\\\\\ z \\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\
G \\\\\\\\ h \\\\end{bmatrix} \\\\] via the following block
elimination scheme: \\\\[ \\\\begin{split} J [A \\\\; b] &= [F \\\\;
f_p] \\\\\\\\ J [C \\\\; d] &= (Jv)_x[A \\\\; b] - [G \\\\; (Jv)_p]
\\\\\\\\ z &= (h + \\\\phi^T C) / \\\\phi^T d \\\\\\\\ X &= A - b z
\\\\\\\\ Y &= -C + d z \\\\end{split} \\\\]

C++ includes: LOCA_TurningPoint_MooreSpence_SalingerBordering.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SalingerBordering::SalingerBordering
"LOCA::TurningPoint::MooreSpence::SalingerBordering::SalingerBordering(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Currently none are
referenced. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SalingerBordering::~SalingerBordering
"LOCA::TurningPoint::MooreSpence::SalingerBordering::~SalingerBordering()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SalingerBordering::setBlocks "void
LOCA::TurningPoint::MooreSpence::SalingerBordering::setBlocks(const
Teuchos::RCP< LOCA::TurningPoint::MooreSpence::AbstractGroup > &group,
const Teuchos::RCP< LOCA::TurningPoint::MooreSpence::ExtendedGroup >
&tpGroup, const Teuchos::RCP< const NOX::Abstract::Vector >
&nullVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&JnVector, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&dfdp, const Teuchos::RCP< const NOX::Abstract::MultiVector > &dJndp)

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

tpGroup:  [in] Turning point group representing the turning point
equations.

nullVector:  [in] Vector representing v

JnVector:  [in] Vector representing Jv

dfdp:  [in] Vector representing df/dp

dJndp:  [in] Vector representing d(Jv)/dp ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SalingerBordering::solve "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::SalingerBordering::solve(Teuchos::ParameterList
&params, const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector
&input, LOCA::TurningPoint::MooreSpence::ExtendedMultiVector &result)
const

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SalingerBordering::solveTranspose "NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::SalingerBordering::solveTranspose(Teuchos::ParameterList
&params, const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector
&input, LOCA::TurningPoint::MooreSpence::ExtendedMultiVector &result)
const

Solves the transpose of the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1SalingerBordering.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::SalingerBordering
"

Moore-Spence pitchfork solver strategy based on \"Salinger\"
bordering. This is the classic 6-solve bordering method.

This class solves the Moore-Spence pitchfork Newton equations: \\\\[
\\\\begin{bmatrix} J & 0 & \\\\psi & f_p \\\\\\\\ (Jv)_x & J & 0 &
(Jv)_p \\\\\\\\ \\\\psi^T & 0 & 0 & 0 \\\\\\\\ 0 & \\\\phi^T & 0 & 0
\\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\\\\\ w \\\\\\\\ z
\\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\ G \\\\\\\\ s \\\\\\\\
h \\\\end{bmatrix} \\\\] via the following block elimination scheme:
\\\\[ \\\\begin{split} J [A \\\\; b \\\\; c] &= [F \\\\; f_p \\\\;
\\\\psi] \\\\\\\\ J [D \\\\; e \\\\; g] &= [G \\\\; (Jv)_p \\\\; 0] -
(Jv)_x[A \\\\; b \\\\; c] \\\\\\\\ w &= \\\\frac{\\\\phi^T
e(\\\\langle A, \\\\psi\\\\rangle - s) - \\\\langle b,
\\\\psi\\\\rangle(\\\\phi^T D - h)} {\\\\phi^T e \\\\langle c,
\\\\psi\\\\rangle - \\\\phi^T g \\\\langle b, \\\\psi\\\\rangle}
\\\\\\\\ z &= \\\\frac{\\\\phi^T D - h - (\\\\phi^T g) w}{\\\\phi^T e}
\\\\\\\\ X &= A - b z - c w \\\\\\\\ Y &= D - e z - g w \\\\end{split}
\\\\]

C++ includes: LOCA_Pitchfork_MooreSpence_SalingerBordering.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SalingerBordering::SalingerBordering "LOCA::Pitchfork::MooreSpence::SalingerBordering::SalingerBordering(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Currently none are
referenced. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SalingerBordering::~SalingerBordering "LOCA::Pitchfork::MooreSpence::SalingerBordering::~SalingerBordering()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SalingerBordering::setBlocks "void
LOCA::Pitchfork::MooreSpence::SalingerBordering::setBlocks(const
Teuchos::RCP< LOCA::Pitchfork::MooreSpence::AbstractGroup > &group,
const Teuchos::RCP< LOCA::Pitchfork::MooreSpence::ExtendedGroup >
&pfGroup, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&asymMultiVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&nullVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&JnVector, const Teuchos::RCP< const NOX::Abstract::Vector > &dfdp,
const Teuchos::RCP< const NOX::Abstract::Vector > &dJndp)

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

pfGroup:  [in] Pitchfork group representing the pitchfork equations.

asymMultiVector:  [in] Multivector representing the asymmetric vector

nullVector:  [in] Vector representing v

JnVector:  [in] Vector representing Jv

dfdp:  [in] Vector representing df/dp

dJndp:  [in] Vector representing d(Jv)/dp ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SalingerBordering::solve "NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::SalingerBordering::solve(Teuchos::ParameterList
&params, const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector
&input, LOCA::Pitchfork::MooreSpence::ExtendedMultiVector &result)
const

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1SalingerBordering.xml
%feature("docstring") LOCA::Hopf::MooreSpence::SalingerBordering "

Moore-Spence Hopf solver strategy based on \"Salinger\" bordering.
This is the classic 5-solve Hopf bordering method.

This class solves the Moore-Spence Hopf Newton equations: \\\\[
\\\\begin{bmatrix} J & 0 & 0 & 0 & f_p \\\\\\\\ (Jy-wBz)_x & J & -wB &
-Bz & (Jy-wBz)_p \\\\\\\\ (Jz+wBy)_x & wB & J & By & (Jz+wBy)_p
\\\\\\\\ 0 & \\\\phi^T & 0 & 0 & 0 \\\\\\\\ 0 & 0 & \\\\phi^T & 0 & 0
\\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\\\\\ Z \\\\\\\\
\\\\omega \\\\\\\\ \\\\lambda \\\\end{bmatrix} = \\\\begin{bmatrix} F
\\\\\\\\ G \\\\\\\\ H \\\\\\\\ u \\\\\\\\ v \\\\end{bmatrix}. \\\\]
via the following block elimination scheme: \\\\[ \\\\begin{split} J
[A \\\\; b] &= [F \\\\; f_p] \\\\\\\\ \\\\begin{bmatrix} J & -wB
\\\\\\\\ wB & J \\\\end{bmatrix} \\\\begin{bmatrix} C & e & g \\\\\\\\
D & f & h \\\\end{bmatrix} &= \\\\begin{bmatrix} G - (Jy-wBz)_x A &
(Jy-wBz)_p - (Jy-wBz)_x b & -Bz \\\\\\\\ H - (Jz+wBy)_x A & (Jz+wBy)_p
- (Jy-wBz)_x b & By \\\\end{bmatrix} \\\\\\\\ \\\\lambda &=
\\\\frac{(\\\\phi^T h)(\\\\phi^T C-u)-(\\\\phi^T g)(\\\\phi^T D-v)}
{(\\\\phi^T h)(\\\\phi^T e)-(\\\\phi^T g)(\\\\phi^T f)} \\\\\\\\
\\\\omega &= \\\\frac{\\\\phi^T D - v - (\\\\phi^T
f)\\\\lambda}{\\\\phi^T h} \\\\\\\\ X &= A - b \\\\lambda \\\\\\\\
\\\\begin{bmatrix} Y \\\\\\\\ Z \\\\end{bmatrix} &= \\\\begin{bmatrix}
C \\\\\\\\ D \\\\end{bmatrix} - \\\\begin{bmatrix} e \\\\\\\\ f
\\\\end{bmatrix}\\\\lambda - \\\\begin{bmatrix} g \\\\\\\\ h
\\\\end{bmatrix}\\\\omega \\\\end{split} \\\\]

C++ includes: LOCA_Hopf_MooreSpence_SalingerBordering.H ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SalingerBordering::SalingerBordering "LOCA::Hopf::MooreSpence::SalingerBordering::SalingerBordering(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

topParams:  [in] Parsed top-level parameter list

solverParams:  [in] Bordered solver parameters. Currently none are
referenced. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SalingerBordering::~SalingerBordering "LOCA::Hopf::MooreSpence::SalingerBordering::~SalingerBordering()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SalingerBordering::setBlocks "void
LOCA::Hopf::MooreSpence::SalingerBordering::setBlocks(const
Teuchos::RCP< LOCA::Hopf::MooreSpence::AbstractGroup > &group, const
Teuchos::RCP< LOCA::Hopf::MooreSpence::ExtendedGroup > &hopfGroup,
const Teuchos::RCP< const NOX::Abstract::Vector > &yVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &zVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &CeRealVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &CeImagVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &dfdp, const Teuchos::RCP<
const NOX::Abstract::Vector > &dCedpReal, const Teuchos::RCP< const
NOX::Abstract::Vector > &dCedpImag, const Teuchos::RCP< const
NOX::Abstract::Vector > &ByVector, const Teuchos::RCP< const
NOX::Abstract::Vector > &mBzVector, double w)

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

hopfGroup:  [in] Hopf group representing the Hopf equations.

yVector:  [in] Vector representing y

zVector:  [in] Vector representing z

CeRealVector:  [in] Vector representing Jy-wBz

CeImagVector:  [in] Vector representing Jz+wBy

dfdp:  [in] Vector representing df/dp

dCedpReal:  [in] Vector representing d(Jy-wBz)/dp

dCedpImag:  [in] Vector representing d(Jz+wBy)/dp

ByVector:  [in] Vector representing By

mBzVector:  [in] Vector representing -Bz

w:  [in] Bifurcation frequency w ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SalingerBordering::solve "NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::SalingerBordering::solve(Teuchos::ParameterList
&params, const LOCA::Hopf::MooreSpence::ExtendedMultiVector &input,
LOCA::Hopf::MooreSpence::ExtendedMultiVector &result) const

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1MultiPredictor_1_1Secant.xml
%feature("docstring") LOCA::MultiPredictor::Secant "

Secant predictor strategy

This class implements a predictor strategy based on computing the
secant vector $v$ to the continuation curve given by \\\\[ v = x - x_o
\\\\] where $x$ is the current solution vector and $x_o$ is the
previous solution vector. Note that for multi-parameter continuation,
the solution component for each secant direction is given as above,
with the parameter components given by the identity matrix.

For the first step of a continuation run, $x_o$ is not defined, and so
a different predictor is used for this step. This predictor is
specified via the \"First Step Predictor\" sublist of the
\"Predictor\" sublist. This predictor is instantiated using the
LOCA::Factory as usual.

The parameters used by this class supplied in the constructor are:
\"First Step Predictor\" Predictor sublist for first step

C++ includes: LOCA_MultiPredictor_Secant.H ";

%feature("docstring")  LOCA::MultiPredictor::Secant::Secant "LOCA::MultiPredictor::Secant::Secant(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &predParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object used for LOCA::Factory

topParams:  [in] Parsed top-level parameter list used when creating
first step predictor

predParams:  [in] Predictor parameters used to obtain \"First Step
Predictor\" as described above. ";

%feature("docstring")  LOCA::MultiPredictor::Secant::~Secant "LOCA::MultiPredictor::Secant::~Secant()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::Secant::Secant "LOCA::MultiPredictor::Secant::Secant(const Secant &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::MultiPredictor::Secant::clone "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::MultiPredictor::Secant::clone(NOX::CopyType type=NOX::DeepCopy)
const

Clone function. ";

%feature("docstring")  LOCA::MultiPredictor::Secant::compute "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Secant::compute(bool baseOnSecant, const
std::vector< double > &stepSize,
LOCA::MultiContinuation::ExtendedGroup &grp, const
LOCA::MultiContinuation::ExtendedVector &prevXVec, const
LOCA::MultiContinuation::ExtendedVector &xVec)

Compute the predictor given the current and previous solution vectors.
Set baseOnSecant to false if the predictor orientation should not be
based on the secant vector (first or last steps of a continuation
run).

This method actually implements the secant calculation described above
";

%feature("docstring")  LOCA::MultiPredictor::Secant::evaluate "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Secant::evaluate(const std::vector< double >
&stepSize, const LOCA::MultiContinuation::ExtendedVector &xVec,
LOCA::MultiContinuation::ExtendedMultiVector &result) const

Evaluate predictor with step size stepSize.

This method computes result[i] = xVec[i] + stepSize[i] * v[i] for each
i, where v[i] is the ith predictor direction. ";

%feature("docstring")  LOCA::MultiPredictor::Secant::computeTangent "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Secant::computeTangent(LOCA::MultiContinuation::ExtendedMultiVector
&tangent)

Compute tangent to predictor and store in tangent. ";

%feature("docstring")  LOCA::MultiPredictor::Secant::isTangentScalable
"bool LOCA::MultiPredictor::Secant::isTangentScalable() const

Is the tangent vector for this predictor scalable.

For the secant predictor, this always returns true. ";


// File: classLOCA_1_1AnasaziOperator_1_1ShiftInvert.xml
%feature("docstring") LOCA::AnasaziOperator::ShiftInvert "

Anasazi operator for computing generalized eigenvalues using shift-
invert.

This class implements the LOCA::AnasaziOperator::AbstractStrategy
interface for computing generalized eigenvalues $\\\\lambda$ and
eigenvectors $z$ of the system \\\\[ J z = \\\\lambda M z *\\\\] where
$J$ is the Jacobian matrix and $M$ is the mass matrix. The right-most
eigenvalues are computed using shift-invert, i.e. solving \\\\[ (J -
\\\\omega M) z = \\\\lambda M z \\\\] where $\\\\omega$ is a real
scalar. The resulting eigenvalue is $\\\\lambda + \\\\omega$.

The parameters used by this class supplied in the constructor are:
\"Shift\" - $\\\\omega$ as defined above (Default 0.0)

Also the grp argument to the constructor must be a child of
LOCA::TimeDependent::AbstractGroup for the shift-invert operations.

C++ includes: LOCA_AnasaziOperator_ShiftInvert.H ";

%feature("docstring")  LOCA::AnasaziOperator::ShiftInvert::ShiftInvert
"LOCA::AnasaziOperator::ShiftInvert::ShiftInvert(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
LOCA::TimeDependent::AbstractGroup > &grp)

Constructor.

Argument grp must be of type LOCA::TimeDependent::AbstractGroup. See
class description for a list of eigenParams. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert::~ShiftInvert "LOCA::AnasaziOperator::ShiftInvert::~ShiftInvert()

Destructor. ";

%feature("docstring")  LOCA::AnasaziOperator::ShiftInvert::label "const std::string & LOCA::AnasaziOperator::ShiftInvert::label() const

Return name of this operator. ";

%feature("docstring")  LOCA::AnasaziOperator::ShiftInvert::apply "void LOCA::AnasaziOperator::ShiftInvert::apply(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &output)
const

Apply the operator.

Applies the inverse of the shifted operator, i.e., solves \\\\[
(J-\\\\omega I)z = M r \\\\] for $z$, where $r = \\\\mbox{input}$ and
$z = \\\\mbox{output}$. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert::transformEigenvalue "void
LOCA::AnasaziOperator::ShiftInvert::transformEigenvalue(double &ev_r,
double &ev_i) const

Transform eigenvalue.

Transforms the given eigenvalue to the eigenvalue of the Jacobian-mass
matrix system by shifting and inverting it. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert::rayleighQuotient "NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::ShiftInvert::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const

Compute Rayleigh quotient.

Computes the Rayleigh quotient $z^T J z / z^T M z$ for the eigenvector
$z$. ";


// File: classLOCA_1_1AnasaziOperator_1_1ShiftInvert2Matrix.xml
%feature("docstring") LOCA::AnasaziOperator::ShiftInvert2Matrix "

Anasazi operator for computing generalized eigenvalues using shift-
invert.

This class implements the LOCA::AnasaziOperator::AbstractStrategy
interface for computing generalized eigenvalues $\\\\lambda$ and
eigenvectors $z$ of the system \\\\[ J z = \\\\lambda M z *\\\\] where
$J$ is the Jacobian matrix and $M$ is the mass matrix. The right-most
eigenvalues are computed using shift-invert, i.e. solving \\\\[ (J -
\\\\omega M) z = \\\\lambda M z \\\\] where $\\\\omega$ is a real
scalar. The resulting eigenvalue is $\\\\lambda + \\\\omega$.

The parameters used by this class supplied in the constructor are:
\"Shift\" - $\\\\omega$ as defined above (Default 0.0)

Also the grp argument to the constructor must be a child of
LOCA::TimeDependent::AbstractGroup for the shift-invert operations.

C++ includes: LOCA_AnasaziOperator_ShiftInvert2Matrix.H ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::ShiftInvert2Matrix "LOCA::AnasaziOperator::ShiftInvert2Matrix::ShiftInvert2Matrix(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
LOCA::TimeDependent::AbstractGroup > &grp)

Constructor.

Argument grp must be of type LOCA::TimeDependent::AbstractGroup. See
class description for a list of eigenParams. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::~ShiftInvert2Matrix "LOCA::AnasaziOperator::ShiftInvert2Matrix::~ShiftInvert2Matrix()

Destructor. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::label "const std::string &
LOCA::AnasaziOperator::ShiftInvert2Matrix::label() const

Return name of this operator. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::apply "void
LOCA::AnasaziOperator::ShiftInvert2Matrix::apply(const
NOX::Abstract::MultiVector &input, NOX::Abstract::MultiVector &output)
const

Apply the operator.

Applies the inverse of the shifted operator, i.e., solves \\\\[
(J-\\\\omega I)z = M r \\\\] for $z$, where $r = \\\\mbox{input}$ and
$z = \\\\mbox{output}$. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::beginPostProcessing "void
LOCA::AnasaziOperator::ShiftInvert2Matrix::beginPostProcessing()

Begin PostProcessing of eigenvalues.

Compute Jacobian and mass matrix once, for use in subsequent repeated
calls to rayleighQuotient ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::transformEigenvalue "void
LOCA::AnasaziOperator::ShiftInvert2Matrix::transformEigenvalue(double
&ev_r, double &ev_i) const

Transform eigenvalue.

Transforms the given eigenvalue to the eigenvalue of the Jacobian-mass
matrix system by shifting and inverting it. ";

%feature("docstring")
LOCA::AnasaziOperator::ShiftInvert2Matrix::rayleighQuotient "NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::ShiftInvert2Matrix::rayleighQuotient(NOX::Abstract::Vector
&evec_r, NOX::Abstract::Vector &evec_i, double &rq_r, double &rq_i)
const

Compute Rayleigh quotient.

Computes the Rayleigh quotient $z^T J z / z^T M z$ for the eigenvector
$z$. ";


// File: classLOCA_1_1Epetra_1_1ShiftInvertInterface.xml
%feature("docstring") LOCA::Epetra::ShiftInvertInterface "

Interface for LOCA::Epetra::ShifterInvertOperator.

C++ includes: LOCA_Epetra_ShiftInvertOperator.H ";

%feature("docstring")
LOCA::Epetra::ShiftInvertInterface::ShiftInvertInterface "LOCA::Epetra::ShiftInvertInterface::ShiftInvertInterface()

Constructor. ";

%feature("docstring")
LOCA::Epetra::ShiftInvertInterface::~ShiftInvertInterface "LOCA::Epetra::ShiftInvertInterface::~ShiftInvertInterface()

Destructor. ";

%feature("docstring")
LOCA::Epetra::ShiftInvertInterface::computeJacobian "bool
LOCA::Epetra::ShiftInvertInterface::computeJacobian(const
Epetra_Vector &x, Epetra_Operator &Jac)

Compute Jacobian $J$. ";


// File: classLOCA_1_1Epetra_1_1ShiftInvertOperator.xml
%feature("docstring") LOCA::Epetra::ShiftInvertOperator "

Epetra operator for $(J-\\\\sigma M)^{-1}$.

C++ includes: LOCA_Epetra_ShiftInvertOperator.H ";

%feature("docstring")
LOCA::Epetra::ShiftInvertOperator::ShiftInvertOperator "LOCA::Epetra::ShiftInvertOperator::ShiftInvertOperator(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
LOCA::Epetra::Group > &grp, const Teuchos::RCP< const Epetra_Operator
> &jac, double shift)

Constructor. ";

%feature("docstring")
LOCA::Epetra::ShiftInvertOperator::~ShiftInvertOperator "LOCA::Epetra::ShiftInvertOperator::~ShiftInvertOperator()

Destructor. ";

%feature("docstring")
LOCA::Epetra::ShiftInvertOperator::SetUseTranspose "int
LOCA::Epetra::ShiftInvertOperator::SetUseTranspose(bool UseTranspose)

Set transpose. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::Apply "int
LOCA::Epetra::ShiftInvertOperator::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply shifted operator. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::ApplyInverse
"int LOCA::Epetra::ShiftInvertOperator::ApplyInverse(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Apply shifted operator inverse. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::NormInf "double LOCA::Epetra::ShiftInvertOperator::NormInf() const

Computing infinity norm. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::Label "const char * LOCA::Epetra::ShiftInvertOperator::Label() const

Label. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::UseTranspose
"bool LOCA::Epetra::ShiftInvertOperator::UseTranspose() const

Transpose. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::HasNormInf "bool LOCA::Epetra::ShiftInvertOperator::HasNormInf() const

Have norm-inf. ";

%feature("docstring")  LOCA::Epetra::ShiftInvertOperator::Comm "const
Epetra_Comm & LOCA::Epetra::ShiftInvertOperator::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")
LOCA::Epetra::ShiftInvertOperator::OperatorDomainMap "const
Epetra_Map & LOCA::Epetra::ShiftInvertOperator::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")
LOCA::Epetra::ShiftInvertOperator::OperatorRangeMap "const Epetra_Map
& LOCA::Epetra::ShiftInvertOperator::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";


// File: classLOCA_1_1Bifurcation_1_1PitchforkBord_1_1StatusTest_1_1SlackUpdateNorm.xml
%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm "

A convergence test based on the update of the slack variable component
for pitchfork location.

Let $\\\\sigma$ be the pitchfork slack variable (see
LOCA::Bifurcation::PitchforkBord::ExtendedGroup). This convergence
test defines convergence for the slack variable when the following is
true \\\\[ \\\\frac{|\\\\sigma-\\\\sigma_0|}{\\\\epsilon_r|\\\\sigma|
+ \\\\epsilon_a} < \\\\tau \\\\] where $\\\\sigma_0$ is the previous
parameter value, $\\\\epsilon_r$ is the relative tolerance,
$\\\\epsilon_a$ is the absolute tolerance, and $\\\\tau$ is an overall
scale factor (typically $\\\\tau = 1$).

Note that this status test deals only with the slack component of the
pitchfork equations. This status test should be combined with other
status tests for the solution and null vector components (using
NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build a
composite status test for the entire system.

Also note that if the group returned by the getSolutionGroup() method
of the solver supplied in checkStatus() is not a pitchfork group
(i.e., not derived from
LOCA::Bifurcation::Pitchforkbord::ExtendedGroup), checkStatus()
returns NOX::StatusTest::Converged. This allows the status test to be
used in situations other than pitchfork tracking, e.g., steady-state
solves, without raising error conditions.

C++ includes: LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getSlackUpdateNorm
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getSlackUpdateNorm()
const

Returns the value of weighted parameter update norm. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getRTOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getRTOL()
const

Returns the realative tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getATOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getATOL()
const

Returns the absolute tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getTOL
"double
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getTOL()
const

Returns the tolerance set in the constructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::SlackUpdateNorm
"LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::SlackUpdateNorm(double
rtol, double atol, double tol)

Constructor.

rtol is the relative tolerance $\\\\epsilon_r$, atol is the absolute
tolerance $\\\\epsilon_a$, and tol is the overall scale factor
$\\\\tau$ defined above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::~SlackUpdateNorm
"LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::~SlackUpdateNorm()

Destructor. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::checkStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::checkStatus(const
NOX::Solver::Generic &problem)

Evaluates convergence criteria specified above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getStatus
"NOX::StatusTest::StatusType
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getStatus()
const

Returns status as defined above. ";

%feature("docstring")
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::print "virtual std::ostream&
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::print(ostream
&stream, int indent=0) const

Prints current status. ";


// File: classLOCA_1_1EigenvalueSort_1_1SmallestImaginary.xml
%feature("docstring") LOCA::EigenvalueSort::SmallestImaginary "

Smallest-imaginary sorting strategy.

Sorts eigenvalues in increasing order according to their imaginary
part. This method requires no parameters in the eigenParams argument
to the constructor

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::SmallestImaginary::SmallestImaginary "LOCA::EigenvalueSort::SmallestImaginary::SmallestImaginary(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")
LOCA::EigenvalueSort::SmallestImaginary::~SmallestImaginary "LOCA::EigenvalueSort::SmallestImaginary::~SmallestImaginary()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::SmallestImaginary::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestImaginary::sort(int n, double *evals,
std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")  LOCA::EigenvalueSort::SmallestImaginary::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestImaginary::sort(int n, double *r_evals,
double *i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1EigenvalueSort_1_1SmallestMagnitude.xml
%feature("docstring") LOCA::EigenvalueSort::SmallestMagnitude "

Smallest-magnitude sorting strategy.

Sorts eigenvalues in increasing order according to their magnitude.
This method requires no parameters in the eigenParams argument to the
constructor

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::SmallestMagnitude::SmallestMagnitude "LOCA::EigenvalueSort::SmallestMagnitude::SmallestMagnitude(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")
LOCA::EigenvalueSort::SmallestMagnitude::~SmallestMagnitude "LOCA::EigenvalueSort::SmallestMagnitude::~SmallestMagnitude()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::SmallestMagnitude::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestMagnitude::sort(int n, double *evals,
std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")  LOCA::EigenvalueSort::SmallestMagnitude::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestMagnitude::sort(int n, double *r_evals,
double *i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1EigenvalueSort_1_1SmallestReal.xml
%feature("docstring") LOCA::EigenvalueSort::SmallestReal "

Smallest-real sorting strategy.

Sorts eigenvalues in increasing order according to their real part.
This method requires no parameters in the eigenParams argument to the
constructor

C++ includes: LOCA_EigenvalueSort_Strategies.H ";

%feature("docstring")
LOCA::EigenvalueSort::SmallestReal::SmallestReal "LOCA::EigenvalueSort::SmallestReal::SmallestReal(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &eigenParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

eigenParams:  [in] Eigensolver parameters. ";

%feature("docstring")
LOCA::EigenvalueSort::SmallestReal::~SmallestReal "LOCA::EigenvalueSort::SmallestReal::~SmallestReal()

Destructor. ";

%feature("docstring")  LOCA::EigenvalueSort::SmallestReal::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestReal::sort(int n, double *evals,
std::vector< int > *perm=NULL) const

Sort real eigenvalues. ";

%feature("docstring")  LOCA::EigenvalueSort::SmallestReal::sort "NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestReal::sort(int n, double *r_evals,
double *i_evals, std::vector< int > *perm=NULL) const

Sort complex eigenvalues. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1SolverFactory.xml
%feature("docstring") LOCA::TurningPoint::MooreSpence::SolverFactory "

Factory for creating solver objects for solving Moore-Spence turning
point equations.

The parameters passed to the create() through the solverParams
argument method should specify the \"Solver Method\" as described
below, as well as any additional parameters for the particular
strategy. \"Solver Method\" - Name of the method. Valid choices are
\"Salinger Bordering\" (
LOCA::TurningPoint::MooreSpence::SalingerBordering) [Default]

\"Phipps Bordering\" (
LOCA::TurningPoint::MooreSpence::PhippsBordering)

C++ includes: LOCA_TurningPoint_MooreSpence_SolverFactory.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverFactory::SolverFactory "LOCA::TurningPoint::MooreSpence::SolverFactory::SolverFactory(const
Teuchos::RCP< LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverFactory::~SolverFactory "LOCA::TurningPoint::MooreSpence::SolverFactory::~SolverFactory()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverFactory::create "Teuchos::RCP<
LOCA::TurningPoint::MooreSpence::SolverStrategy >
LOCA::TurningPoint::MooreSpence::SolverFactory::create(const
Teuchos::RCP< LOCA::Parameter::SublistParser > &topParams, const
Teuchos::RCP< Teuchos::ParameterList > &solverParams)

Create solver strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

solverParams:  [in] Solver parameters as described above ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverFactory::strategyName "const
std::string &
LOCA::TurningPoint::MooreSpence::SolverFactory::strategyName(Teuchos::ParameterList
&solverParams) const

Return strategy name given by solverParams. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1SolverFactory.xml
%feature("docstring") LOCA::Hopf::MooreSpence::SolverFactory "

Factory for creating solver objects for solving Moore-Spence Hopf
equations.

The parameters passed to the create() through the solverParams
argument method should specify the \"Solver Method\" as described
below, as well as any additional parameters for the particular
strategy. \"Solver Method\" - Name of the method. Valid choices are
\"Salinger Bordering\" ( LOCA::Hopf::MooreSpence::SalingerBordering)
[Default]

C++ includes: LOCA_Hopf_MooreSpence_SolverFactory.H ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SolverFactory::SolverFactory "LOCA::Hopf::MooreSpence::SolverFactory::SolverFactory(const
Teuchos::RCP< LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SolverFactory::~SolverFactory "LOCA::Hopf::MooreSpence::SolverFactory::~SolverFactory()

Destructor. ";

%feature("docstring")  LOCA::Hopf::MooreSpence::SolverFactory::create
"Teuchos::RCP< LOCA::Hopf::MooreSpence::SolverStrategy >
LOCA::Hopf::MooreSpence::SolverFactory::create(const Teuchos::RCP<
LOCA::Parameter::SublistParser > &topParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Create solver strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

solverParams:  [in] Solver parameters as described above ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SolverFactory::strategyName "const
std::string &
LOCA::Hopf::MooreSpence::SolverFactory::strategyName(Teuchos::ParameterList
&solverParams) const

Return strategy name given by solverParams. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1SolverFactory.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::SolverFactory "

Factory for creating solver objects for solving Moore-Spence pitchfork
equations.

The parameters passed to the create() through the solverParams
argument method should specify the \"Solver Method\" as described
below, as well as any additional parameters for the particular
strategy. \"Solver Method\" - Name of the method. Valid choices are
\"Salinger Bordering\" (
LOCA::Pitchfork::MooreSpence::SalingerBordering) [Default]

\"Phipps Bordering\" ( LOCA::Pitchfork::MooreSpence::PhippsBordering)

C++ includes: LOCA_Pitchfork_MooreSpence_SolverFactory.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverFactory::SolverFactory "LOCA::Pitchfork::MooreSpence::SolverFactory::SolverFactory(const
Teuchos::RCP< LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverFactory::~SolverFactory "LOCA::Pitchfork::MooreSpence::SolverFactory::~SolverFactory()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverFactory::create "Teuchos::RCP<
LOCA::Pitchfork::MooreSpence::SolverStrategy >
LOCA::Pitchfork::MooreSpence::SolverFactory::create(const
Teuchos::RCP< LOCA::Parameter::SublistParser > &topParams, const
Teuchos::RCP< Teuchos::ParameterList > &solverParams)

Create solver strategy.

Parameters:
-----------

topParams:  [in] Parsed top-level parameter list.

solverParams:  [in] Solver parameters as described above ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverFactory::strategyName "const
std::string &
LOCA::Pitchfork::MooreSpence::SolverFactory::strategyName(Teuchos::ParameterList
&solverParams) const

Return strategy name given by solverParams. ";


// File: classLOCA_1_1TurningPoint_1_1MooreSpence_1_1SolverStrategy.xml
%feature("docstring") LOCA::TurningPoint::MooreSpence::SolverStrategy
"

Abstract strategy for solving the Moore-Spence turning point
equations.

This class provides an abstract interface for solver strategies to
solve the Moore-Spence turning point Newton system: \\\\[
\\\\begin{bmatrix} J & 0 & f_p \\\\\\\\ (Jv)_x & J & (Jv)_p \\\\\\\\ 0
& \\\\phi^T & 0 \\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y
\\\\\\\\ z \\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\ G \\\\\\\\
h \\\\end{bmatrix}. \\\\] After instantiating a solver (via
LOCA::TurningPoint::MooreSpence::SolverFactory), the linear system is
set up by setBlocks() and can then be solved by solve().

C++ includes: LOCA_TurningPoint_MooreSpence_SolverStrategy.H ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverStrategy::SolverStrategy "LOCA::TurningPoint::MooreSpence::SolverStrategy::SolverStrategy()

Constructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverStrategy::~SolverStrategy "virtual
LOCA::TurningPoint::MooreSpence::SolverStrategy::~SolverStrategy()

Destructor. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverStrategy::setBlocks "virtual
void LOCA::TurningPoint::MooreSpence::SolverStrategy::setBlocks(const
Teuchos::RCP< LOCA::TurningPoint::MooreSpence::AbstractGroup > &group,
const Teuchos::RCP< LOCA::TurningPoint::MooreSpence::ExtendedGroup >
&tpGroup, const Teuchos::RCP< const NOX::Abstract::Vector >
&nullVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&JnVector, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&dfdp, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&dJndp)=0

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

tpGroup:  [in] Turning point group representing the turning point
equations.

nullVector:  [in] Vector representing v

JnVector:  [in] Vector representing Jv

dfdp:  [in] Vector representing df/dp

dJndp:  [in] Vector representing d(Jv)/dp ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverStrategy::solve "virtual
NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::SolverStrategy::solve(Teuchos::ParameterList
&params, const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector
&input, LOCA::TurningPoint::MooreSpence::ExtendedMultiVector &result)
const =0

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";

%feature("docstring")
LOCA::TurningPoint::MooreSpence::SolverStrategy::solveTranspose "virtual NOX::Abstract::Group::ReturnType
LOCA::TurningPoint::MooreSpence::SolverStrategy::solveTranspose(Teuchos::ParameterList
&params, const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector
&input, LOCA::TurningPoint::MooreSpence::ExtendedMultiVector &result)
const

Solves the transpose of the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1Pitchfork_1_1MooreSpence_1_1SolverStrategy.xml
%feature("docstring") LOCA::Pitchfork::MooreSpence::SolverStrategy "

Abstract strategy for solving the Moore-Spence pitchfork equations.

This class provides an abstract interface for solver strategies to
solve the Moore-Spence pitchfork Newton system: \\\\[
\\\\begin{bmatrix} J & 0 & \\\\psi & f_p \\\\\\\\ (Jv)_x & J & 0 &
(Jv)_p \\\\\\\\ \\\\psi^T & 0 & 0 & 0 \\\\\\\\ 0 & \\\\phi^T & 0 & 0
\\\\end{bmatrix} \\\\begin{bmatrix} X \\\\\\\\ Y \\\\\\\\ w \\\\\\\\ z
\\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\ G \\\\\\\\ s \\\\\\\\
h \\\\end{bmatrix}. \\\\] After instantiating a solver Solvers (via
LOCA::Pitchfork::MooreSpence::SolverFactory), the linear system is set
up by setBlocks() and can then be solved by solve().

C++ includes: LOCA_Pitchfork_MooreSpence_SolverStrategy.H ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverStrategy::SolverStrategy "LOCA::Pitchfork::MooreSpence::SolverStrategy::SolverStrategy()

Constructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverStrategy::~SolverStrategy "virtual
LOCA::Pitchfork::MooreSpence::SolverStrategy::~SolverStrategy()

Destructor. ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverStrategy::setBlocks "virtual void
LOCA::Pitchfork::MooreSpence::SolverStrategy::setBlocks(const
Teuchos::RCP< LOCA::Pitchfork::MooreSpence::AbstractGroup > &group,
const Teuchos::RCP< LOCA::Pitchfork::MooreSpence::ExtendedGroup >
&pfGroup, const Teuchos::RCP< const NOX::Abstract::MultiVector >
&asymMultiVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&nullVector, const Teuchos::RCP< const NOX::Abstract::Vector >
&JnVector, const Teuchos::RCP< const NOX::Abstract::Vector > &dfdp,
const Teuchos::RCP< const NOX::Abstract::Vector > &dJndp)=0

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

pfGroup:  [in] Pitchfork group representing the pitchfork equations.

asymMultiVector:  [in] Multivector representing the asymmetric vector

nullVector:  [in] Vector representing v

JnVector:  [in] Vector representing Jv

dfdp:  [in] Vector representing df/dp

dJndp:  [in] Vector representing d(Jv)/dp ";

%feature("docstring")
LOCA::Pitchfork::MooreSpence::SolverStrategy::solve "virtual
NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MooreSpence::SolverStrategy::solve(Teuchos::ParameterList
&params, const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector
&input, LOCA::Pitchfork::MooreSpence::ExtendedMultiVector &result)
const =0

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1Hopf_1_1MooreSpence_1_1SolverStrategy.xml
%feature("docstring") LOCA::Hopf::MooreSpence::SolverStrategy "

Abstract strategy for solving the Moore-Spence Hopf equations.

This class provides an abstract interface for solver strategies to
solve the Moore-Spence Hopf Newton system: \\\\[ \\\\begin{bmatrix} J
& 0 & 0 & 0 & f_p \\\\\\\\ (Jy-wBz)_x & J & -wB & -Bz & (Jy-wBz)_p
\\\\\\\\ (Jz+wBy)_x & wB & J & By & (Jz+wBy)_p \\\\\\\\ 0 & \\\\phi^T
& 0 & 0 & 0 \\\\\\\\ 0 & 0 & \\\\phi^T & 0 & 0 \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\\\\\ Z \\\\\\\\ \\\\Omega^T
\\\\\\\\ \\\\Lambda^T \\\\end{bmatrix} = \\\\begin{bmatrix} F \\\\\\\\
G \\\\\\\\ H \\\\\\\\ U^T \\\\\\\\ V^T \\\\end{bmatrix}. \\\\] After
instantiating a solver (via LOCA::Hopf::MooreSpence::SolverFactory),
the linear system is set up by setBlocks() and can then be solved by
solve().

C++ includes: LOCA_Hopf_MooreSpence_SolverStrategy.H ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SolverStrategy::SolverStrategy "LOCA::Hopf::MooreSpence::SolverStrategy::SolverStrategy()

Constructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SolverStrategy::~SolverStrategy "virtual
LOCA::Hopf::MooreSpence::SolverStrategy::~SolverStrategy()

Destructor. ";

%feature("docstring")
LOCA::Hopf::MooreSpence::SolverStrategy::setBlocks "virtual void
LOCA::Hopf::MooreSpence::SolverStrategy::setBlocks(const Teuchos::RCP<
LOCA::Hopf::MooreSpence::AbstractGroup > &group, const Teuchos::RCP<
LOCA::Hopf::MooreSpence::ExtendedGroup > &hopfGroup, const
Teuchos::RCP< const NOX::Abstract::Vector > &yVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &zVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &CeRealVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &CeImagVector, const
Teuchos::RCP< const NOX::Abstract::Vector > &dfdp, const Teuchos::RCP<
const NOX::Abstract::Vector > &dCedpReal, const Teuchos::RCP< const
NOX::Abstract::Vector > &dCedpImag, const Teuchos::RCP< const
NOX::Abstract::Vector > &ByVector, const Teuchos::RCP< const
NOX::Abstract::Vector > &mBzVector, double w)=0

Set blocks in extended linear system.

Parameters:
-----------

group:  [in] Underlying group representing J

hopfGroup:  [in] Hopf group representing the Hopf equations.

yVector:  [in] Vector representing y

zVector:  [in] Vector representing z

CeRealVector:  [in] Vector representing Jy-wBz

CeImagVector:  [in] Vector representing Jz+wBy

dfdp:  [in] Vector representing df/dp

dCedpReal:  [in] Vector representing d(Jy-wBz)/dp

dCedpImag:  [in] Vector representing d(Jz+wBy)/dp

ByVector:  [in] Vector representing By

mBzVector:  [in] Vector representing -Bz

w:  [in] Bifurcation frequency w ";

%feature("docstring")  LOCA::Hopf::MooreSpence::SolverStrategy::solve
"virtual NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::SolverStrategy::solve(Teuchos::ParameterList
&params, const LOCA::Hopf::MooreSpence::ExtendedMultiVector &input,
LOCA::Hopf::MooreSpence::ExtendedMultiVector &result) const =0

Solves the extended system as defined above.

The params argument is the linear solver parameters. ";


// File: classLOCA_1_1Parameter_1_1StandardEntry.xml
%feature("docstring") LOCA::Parameter::StandardEntry "

Standard parameter entry class using a function object.

This is the standard parameter entry class that uses a function object
to actually set/retrieve parameter values. The nice thing about using
a function object is it allows one to set parameters that don't
actually exist in the code, for example, setting a dimensionless group
value by modifiying a number of physical parameters. By supplying an
appropriate function object, this class should suffice for
setting/retrieving parameter values in nearly all cases.

The constructor takes a pointer to the supplied function object. It is
assumed that this class then owns that pointer, and in particular,
calls delete in the destructor if the entry is successfully added to
the library. It does not delete the function object otherwise.

C++ includes: LOCA_Parameter_Entry.H ";

%feature("docstring")  LOCA::Parameter::StandardEntry::StandardEntry "LOCA::Parameter::StandardEntry< FunctorType, ValueType
>::StandardEntry(FunctorType *fctr)

Constructor. ";

%feature("docstring")  LOCA::Parameter::StandardEntry::~StandardEntry
"virtual LOCA::Parameter::StandardEntry< FunctorType, ValueType
>::~StandardEntry()

Destructor. ";

%feature("docstring")  LOCA::Parameter::StandardEntry::setValue "virtual void LOCA::Parameter::StandardEntry< FunctorType, ValueType
>::setValue(const ValueType &value)

Set parameter this object represents to value. ";

%feature("docstring")  LOCA::Parameter::StandardEntry::getValue "virtual ValueType LOCA::Parameter::StandardEntry< FunctorType,
ValueType >::getValue() const

Get parameter value this object represents. ";

%feature("docstring")  LOCA::Parameter::StandardEntry::setIsInLibrary
"virtual void LOCA::Parameter::StandardEntry< FunctorType, ValueType
>::setIsInLibrary()

Informs entry that it is now stored in the library.

This is used primarily for informing the entry on how to delete itself
when deleting the library. ";


// File: classLOCA_1_1Stepper.xml
%feature("docstring") LOCA::Stepper "

Implementation of LOCA::Abstract::Iterator for computing points along
a continuation curve.

The Stepper class implements the pure virtual methods of the
LOCA::Abstract::Iterator for iteratively computing points along a
continuation curve.

C++ includes: LOCA_Stepper.H ";

%feature("docstring")  LOCA::Stepper::Stepper "LOCA::Stepper::Stepper(const Teuchos::RCP< LOCA::GlobalData >
&global_data, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &initialGuess, const
Teuchos::RCP< LOCA::StatusTest::Abstract > &lt, const Teuchos::RCP<
NOX::StatusTest::Generic > &nt, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Constructor with LOCA::StatusTest. ";

%feature("docstring")  LOCA::Stepper::Stepper "LOCA::Stepper::Stepper(const Teuchos::RCP< LOCA::GlobalData >
&global_data, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &initialGuess, const
Teuchos::RCP< NOX::StatusTest::Generic > &nt, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Obsolete constructor without LOCA::StatusTest.Deprecated Use the
constructor with LOCA::StatusTest instead. ";

%feature("docstring")  LOCA::Stepper::~Stepper "LOCA::Stepper::~Stepper()

Destructor. ";

%feature("docstring")  LOCA::Stepper::reset "bool
LOCA::Stepper::reset(const Teuchos::RCP< LOCA::GlobalData >
&global_data, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &initialGuess, const
Teuchos::RCP< LOCA::StatusTest::Abstract > &lt, const Teuchos::RCP<
NOX::StatusTest::Generic > &nt, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Reset the Stepper to start a new continuation run. Version with
LOCA::StatusTest ";

%feature("docstring")  LOCA::Stepper::reset "bool
LOCA::Stepper::reset(const Teuchos::RCP< LOCA::GlobalData >
&global_data, const Teuchos::RCP<
LOCA::MultiContinuation::AbstractGroup > &initialGuess, const
Teuchos::RCP< NOX::StatusTest::Generic > &nt, const Teuchos::RCP<
Teuchos::ParameterList > &p)

Reset the Stepper to start a new continuation run. Obsolete version
without LOCA::StatusTest.Deprecated Use reset() with LOCA::StatusTest
instead. ";

%feature("docstring")  LOCA::Stepper::eigensolverReset "bool
LOCA::Stepper::eigensolverReset(Teuchos::RCP< Teuchos::ParameterList >
&newEigensolverList)

Replaces the eigensolver parameter list. ";

%feature("docstring")  LOCA::Stepper::getSolutionGroup "Teuchos::RCP<
const LOCA::MultiContinuation::AbstractGroup >
LOCA::Stepper::getSolutionGroup() const

Return the current solution group. ";

%feature("docstring")  LOCA::Stepper::getBifurcationGroup "Teuchos::RCP< const LOCA::MultiContinuation::AbstractGroup >
LOCA::Stepper::getBifurcationGroup() const

Return the current bifurcation group.

If the current bifurcation method is \"None\", then the returned group
is the same as getSolutionGroup(), otherwise this method returns the
current bifurcation group (e.g., a turning point group). ";

%feature("docstring")  LOCA::Stepper::getList "Teuchos::RCP< const
Teuchos::ParameterList > LOCA::Stepper::getList() const

Return the output parameters from the stepper algorithm. ";

%feature("docstring")  LOCA::Stepper::getSolver "Teuchos::RCP< const
NOX::Solver::Generic > LOCA::Stepper::getSolver() const

Return the current nonlinear solver pointer.

Will throw an error if the solver does not exist yet. ";

%feature("docstring")  LOCA::Stepper::getContinuationParameter "double LOCA::Stepper::getContinuationParameter() const

Return the current continuation parameter from the underlying
LOCA::MultiContinuation::AbstractStrategy. ";


// File: classLOCA_1_1Parameter_1_1SublistParser.xml
%feature("docstring") LOCA::Parameter::SublistParser "

Class to parse a parameter list for sublists.

This class parses a supplied parameter list and locates various
sublists. This saves the code from having to traverse the parameter
list to find sublists itself, and puts in one location the hard-coded
structure of the parameter list.

C++ includes: LOCA_Parameter_SublistParser.H ";

%feature("docstring")  LOCA::Parameter::SublistParser::SublistParser "LOCA::Parameter::SublistParser::SublistParser(const Teuchos::RCP<
LOCA::GlobalData > &global_data)

Constructor. ";

%feature("docstring")  LOCA::Parameter::SublistParser::~SublistParser
"LOCA::Parameter::SublistParser::~SublistParser()

Destructor. ";

%feature("docstring")  LOCA::Parameter::SublistParser::parseSublists "void LOCA::Parameter::SublistParser::parseSublists(const Teuchos::RCP<
Teuchos::ParameterList > &topLevelParams)

Parse parameter list to find sublists. ";

%feature("docstring")  LOCA::Parameter::SublistParser::getSublist "Teuchos::RCP< Teuchos::ParameterList >
LOCA::Parameter::SublistParser::getSublist(const std::string &name)

Return sublist of name name. ";


// File: classLOCA_1_1MultiPredictor_1_1Tangent.xml
%feature("docstring") LOCA::MultiPredictor::Tangent "

Tangent predictor strategy

This class implements a predictor strategy based on computing the
tangent to the continuation curve. If $p$ is the vector of
continuation parameters, then the solution component of the tangent
vectors $v_x$ are computed by solving \\\\[ J v_x = -
\\\\frac{\\\\partial f}{\\\\partial p}. \\\\] The parameter component
$v_p$ is set to the identity matrix.

C++ includes: LOCA_MultiPredictor_Tangent.H ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::Tangent "LOCA::MultiPredictor::Tangent::Tangent(const Teuchos::RCP<
LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &predParams, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object

predParams:  [in] Predictor parameters. Currently no parameters are
used by the Tangent predictor.

solverParams:  [in] Linear solver parameters used in linear solve to
compute tangent vectors $v$. ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::~Tangent "LOCA::MultiPredictor::Tangent::~Tangent()

Destructor. ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::Tangent "LOCA::MultiPredictor::Tangent::Tangent(const Tangent &source,
NOX::CopyType type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::clone "Teuchos::RCP< LOCA::MultiPredictor::AbstractStrategy >
LOCA::MultiPredictor::Tangent::clone(NOX::CopyType type=NOX::DeepCopy)
const

Clone function. ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::compute "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Tangent::compute(bool baseOnSecant, const
std::vector< double > &stepSize,
LOCA::MultiContinuation::ExtendedGroup &grp, const
LOCA::MultiContinuation::ExtendedVector &prevXVec, const
LOCA::MultiContinuation::ExtendedVector &xVec)

Compute the predictor given the current and previous solution vectors.
Set baseOnSecant to false if the predictor orientation should not be
based on the secant vector (first or last steps of a continuation
run).

This method actually implements the predictor solve described above ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::evaluate "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Tangent::evaluate(const std::vector< double >
&stepSize, const LOCA::MultiContinuation::ExtendedVector &xVec,
LOCA::MultiContinuation::ExtendedMultiVector &result) const

Evaluate predictor with step size stepSize.

This method computes result[i] = xVec[i] + stepSize[i] * v[i] for each
i, where v[i] is the ith predictor direction. ";

%feature("docstring")  LOCA::MultiPredictor::Tangent::computeTangent "NOX::Abstract::Group::ReturnType
LOCA::MultiPredictor::Tangent::computeTangent(LOCA::MultiContinuation::ExtendedMultiVector
&tangent)

Compute tangent to predictor and store in tangent. ";

%feature("docstring")
LOCA::MultiPredictor::Tangent::isTangentScalable "bool
LOCA::MultiPredictor::Tangent::isTangentScalable() const

Is the tangent vector for this predictor scalable.

For the tangent predictor, this always returns true. ";


// File: classLOCA_1_1Epetra_1_1Interface_1_1TimeDependent.xml
%feature("docstring") LOCA::Epetra::Interface::TimeDependent "

Used by LOCA::Epetra::Group to provide a link to the external code for
computing the shifted matrix.

This interface is derived from the NOX::Epetra::Interface::Jacobian
and additionally provides a method for computing the shifted matrix
$\\\\alpha J + \\\\beta M$. This is needed for linear stability
analysis and Hopf tracking.

C++ includes: LOCA_Epetra_Interface_TimeDependent.H ";

%feature("docstring")
LOCA::Epetra::Interface::TimeDependent::TimeDependent "LOCA::Epetra::Interface::TimeDependent::TimeDependent()

Constructor. ";

%feature("docstring")
LOCA::Epetra::Interface::TimeDependent::~TimeDependent "virtual
LOCA::Epetra::Interface::TimeDependent::~TimeDependent()

Destructor. ";

%feature("docstring")
LOCA::Epetra::Interface::TimeDependent::computeShiftedMatrix "virtual
bool
LOCA::Epetra::Interface::TimeDependent::computeShiftedMatrix(double
alpha, double beta, const Epetra_Vector &x, Epetra_Operator &A)=0

Call user routine for computing the shifted matrix $\\\\alpha J +
\\\\beta M$ where $J$ is the Jacobian matrix and $M$ is the mass
matrix. ";

%feature("docstring")  LOCA::Epetra::Interface::TimeDependent::setXdot
"virtual void LOCA::Epetra::Interface::TimeDependent::setXdot(const
Epetra_Vector &xdot, const double time)

Routine used in XYZT to set x_dot and time in the interface.

The computeF() routine for XYZT problems needs to be a function of
x_dot, but th NOX/LOCA computeF() does not take x_dot as an argument.
This is used to set x_dot in the application interface so the
subsequent call to computeF has the correct x_dot value. The timeStep
argument is sent so the use can set the global time, for cases when
computeF, computeJacobian, computeMassMatrix fills are functions of
time (nonautonomous systems). ";


// File: classLOCA_1_1Epetra_1_1Interface_1_1TimeDependentMatrixFree.xml
%feature("docstring") LOCA::Epetra::Interface::TimeDependentMatrixFree
"

Used by LOCA::Epetra::Group to provide a link to the external code for
applying the shifted matrix in a matrix-free setting.

This interface is derived from the NOX::Epetra::Interface::Required
and additionally provides a method for applying the shifted matrix
$\\\\alpha J + \\\\beta M$. This is needed for linear stability
analysis and Hopf tracking.

C++ includes: LOCA_Epetra_Interface_TimeDependentMatrixFree.H ";

%feature("docstring")
LOCA::Epetra::Interface::TimeDependentMatrixFree::TimeDependentMatrixFree
"LOCA::Epetra::Interface::TimeDependentMatrixFree::TimeDependentMatrixFree()

Constructor. ";

%feature("docstring")
LOCA::Epetra::Interface::TimeDependentMatrixFree::~TimeDependentMatrixFree
"virtual
LOCA::Epetra::Interface::TimeDependentMatrixFree::~TimeDependentMatrixFree()

Destructor. ";

%feature("docstring")
LOCA::Epetra::Interface::TimeDependentMatrixFree::applyShiftedMatrix "virtual bool
LOCA::Epetra::Interface::TimeDependentMatrixFree::applyShiftedMatrix(double
alpha, double beta, const NOX::Epetra::Vector &input,
NOX::Epetra::Vector &result) const =0

Call user routine for applying the shifted matrix $\\\\alpha J +
\\\\beta M$ where $J$ is the Jacobian matrix and $M$ is the mass
matrix. ";


// File: classLOCA_1_1Epetra_1_1TransposeLinearSystem_1_1TransposePreconditioner.xml
%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner "

Method for solving the transpose of a linear system by using the
transpose of the preconditioner.

C++ includes:
LOCA_Epetra_TransposeLinearSystem_TransposePreconditioner.H ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::TransposePreconditioner
"LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::TransposePreconditioner(const
Teuchos::RCP< LOCA::GlobalData > &global_data, const Teuchos::RCP<
Teuchos::ParameterList > &solverParams, const Teuchos::RCP<
NOX::Epetra::LinearSystem > &linsys)

Constructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::~TransposePreconditioner
"LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::~TransposePreconditioner()

Destructor. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::applyJacobianTransposeInverse
"bool
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Epetra::Vector &input, NOX::Epetra::Vector
&result)

Applies the inverse of the Jacobian matrix transpose to the given
input vector and puts the answer in result.

Computes \\\\[ v = J^{-T} u, \\\\] where $J$ is the Jacobian, $u$ is
the input vector, and $v$ is the result vector.

The parameter list contains the linear solver options. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::createJacobianTranspose
"bool
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::createJacobianTranspose()

Evaluates the Jacobian-transpose based on the solution vector x.

Note: For flexibility, this method does not compute the original
Jacobian matrix. It uses whatever is currently stored in the linear
system. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::createTransposePreconditioner
"bool
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::createTransposePreconditioner(const
NOX::Epetra::Vector &x, Teuchos::ParameterList &p)

Explicitly constructs a preconditioner based on the solution vector x
and the parameter list p.

Note: x is only needed for user-supplied preconditioners. When using a
built- in preconditioner (e.g., Ifpack), x will note be used. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::getJacobianTransposeOperator
"Teuchos::RCP< Epetra_Operator >
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::getJacobianTransposeOperator()

Get Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::getTransposePreconditioner
"Teuchos::RCP< Epetra_Operator >
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::getTransposePreconditioner()

Get transpose-preconditioner. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::setJacobianTransposeOperator
"void
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::setJacobianTransposeOperator(const
Teuchos::RCP< Epetra_Operator > &new_jac_trans)

Set Jacobian-transpose operator. ";

%feature("docstring")
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::setTransposePreconditioner
"void
LOCA::Epetra::TransposeLinearSystem::TransposePreconditioner::setTransposePreconditioner(const
Teuchos::RCP< Epetra_Operator > &new_prec_trans)

Set transpose-preconditioner. ";


// File: classLOCA_1_1Abstract_1_1TransposeSolveGroup.xml
%feature("docstring") LOCA::Abstract::TransposeSolveGroup "

Abstract group interface class for solving the transpose of the
Jacobian.

This interface, derived from NOX::Abstract::Group, provides the
additional interface for solving the transpose of the Jacobian.

C++ includes: LOCA_Abstract_TransposeSolveGroup.H ";

%feature("docstring")
LOCA::Abstract::TransposeSolveGroup::TransposeSolveGroup "LOCA::Abstract::TransposeSolveGroup::TransposeSolveGroup()

Constructor. ";

%feature("docstring")
LOCA::Abstract::TransposeSolveGroup::~TransposeSolveGroup "virtual
LOCA::Abstract::TransposeSolveGroup::~TransposeSolveGroup()

Destructor. ";

%feature("docstring")
LOCA::Abstract::TransposeSolveGroup::applyJacobianTransposeInverse "virtual NOX::Abstract::Group::ReturnType
LOCA::Abstract::TransposeSolveGroup::applyJacobianTransposeInverse(Teuchos::ParameterList
&params, const NOX::Abstract::Vector &input, NOX::Abstract::Vector
&result) const =0

Solve Jacobian-tranpose system. ";

%feature("docstring")
LOCA::Abstract::TransposeSolveGroup::applyJacobianTransposeInverseMultiVector
"virtual NOX::Abstract::Group::ReturnType
LOCA::Abstract::TransposeSolveGroup::applyJacobianTransposeInverseMultiVector(Teuchos::ParameterList
&params, const NOX::Abstract::MultiVector &input,
NOX::Abstract::MultiVector &result) const =0

Solve Jacobian-tranpose system with multiple right-hand sides. ";


// File: classLOCA_1_1BorderedSolver_1_1UpperTriangularBlockElimination.xml
%feature("docstring")
LOCA::BorderedSolver::UpperTriangularBlockElimination "

Block elimination strategy for solving a block upper-triangular
system.

This class solves the extended system of equations \\\\[
\\\\begin{bmatrix} op(J) & A \\\\\\\\ 0 & op(C) \\\\end{bmatrix}
\\\\begin{bmatrix} X \\\\\\\\ Y \\\\end{bmatrix} = \\\\begin{bmatrix}
F \\\\\\\\ G \\\\end{bmatrix} \\\\] via block elimination: \\\\[
\\\\begin{aligned} Y &= op(C)^{-1}G \\\\\\\\ X &= op(J)^{-1}(F - A Y)
\\\\end{aligned} \\\\] where $op$ represents either the identity
operation or the transpose. $C$ must be nonzero, while $A$, $F$ or $G$
may be zero. The solve for the non-transposed system is implemented by
the solve() method, while the solve for the transposed system is
implemented by the solveTranspose() method.

C++ includes: LOCA_BorderedSolver_UpperTriangularBlockElimination.H ";

%feature("docstring")
LOCA::BorderedSolver::UpperTriangularBlockElimination::UpperTriangularBlockElimination
"LOCA::BorderedSolver::UpperTriangularBlockElimination::UpperTriangularBlockElimination(const
Teuchos::RCP< LOCA::GlobalData > &global_data)

Constructor.

Parameters:
-----------

global_data:  [in] Global data object ";

%feature("docstring")
LOCA::BorderedSolver::UpperTriangularBlockElimination::~UpperTriangularBlockElimination
"LOCA::BorderedSolver::UpperTriangularBlockElimination::~UpperTriangularBlockElimination()

Destructor. ";

%feature("docstring")
LOCA::BorderedSolver::UpperTriangularBlockElimination::solve "NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::UpperTriangularBlockElimination::solve(Teuchos::ParameterList
&params, const LOCA::BorderedSolver::AbstractOperator &op, const
NOX::Abstract::MultiVector *A, const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system as described above.

Either A, F, or G may be zero by passing NULL. ";

%feature("docstring")
LOCA::BorderedSolver::UpperTriangularBlockElimination::solveTranspose
"NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::UpperTriangularBlockElimination::solveTranspose(Teuchos::ParameterList
&params, const LOCA::BorderedSolver::AbstractOperator &op, const
NOX::Abstract::MultiVector *A, const
NOX::Abstract::MultiVector::DenseMatrix &C, const
NOX::Abstract::MultiVector *F, const
NOX::Abstract::MultiVector::DenseMatrix *G, NOX::Abstract::MultiVector
&X, NOX::Abstract::MultiVector::DenseMatrix &Y) const

Solves the extended system using the tranpose of J and C as described
above.

Either A, F, or G may be zero by passing NULL. ";


// File: classLOCA_1_1Extended_1_1Vector.xml
%feature("docstring") LOCA::Extended::Vector "

Implemenatation of the NOX::Abstract::Vector class for extended
vectors comprised of an arbitrary number of vectors and scalars.

Many continuation and bifurcation calculations can be viewed as the
solution to an extended set of equations. For example, calculating a
turning point can be viewed as computing a solution to $G(z) = 0$
where $z = [x, n, p]\\\\in\\\\Re^{2n+1}$ and \\\\[ G(z) = \\\\left[
\\\\begin{array}{c} F(x,p) \\\\\\\\ Jn \\\\\\\\ n^Tn-1 \\\\end{array}
\\\\right] \\\\] The extended vector $z$ is comprised of the two
vectors $x$ and $n$ as well as the scalar $p$. This class provides an
implementation of the NOX::Abstract::Vector interface for such
extended vectors. It stores an array of pointers to
NOX::Abstract::Vector's as well as an array of scalars using the STL
vector class.

The implementations of the NOX::Abstract::Vector methods are defined
in terms of the implementations of each stored abstract vector.

C++ includes: LOCA_Extended_Vector.H ";

%feature("docstring")  LOCA::Extended::Vector::Vector "LOCA::Extended::Vector::Vector(const Vector &source, NOX::CopyType
type=NOX::DeepCopy)

Copy constructor. ";

%feature("docstring")  LOCA::Extended::Vector::~Vector "LOCA::Extended::Vector::~Vector()

Vector destructor. ";

%feature("docstring")  LOCA::Extended::Vector::clone "Teuchos::RCP<
NOX::Abstract::Vector > LOCA::Extended::Vector::clone(NOX::CopyType
type=NOX::DeepCopy) const

Clone function. Applies clone to each stored vector. ";

%feature("docstring")  LOCA::Extended::Vector::createMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::Vector::createMultiVector(const NOX::Abstract::Vector
*const *vecs, int numVecs, NOX::CopyType type=NOX::DeepCopy) const

Create a MultiVector with numVecs+1 columns out of an array of
Vectors. The vector stored under this will be the first column with
the remaining numVecs columns given by vecs. ";

%feature("docstring")  LOCA::Extended::Vector::createMultiVector "Teuchos::RCP< NOX::Abstract::MultiVector >
LOCA::Extended::Vector::createMultiVector(int numVecs, NOX::CopyType
type=NOX::DeepCopy) const

Create a MultiVector with numVecs columns. ";

%feature("docstring")  LOCA::Extended::Vector::init "NOX::Abstract::Vector & LOCA::Extended::Vector::init(double gamma)

NOX::Abstract::Vector init function. Initializes each stored vector
and scalar. ";

%feature("docstring")  LOCA::Extended::Vector::random "NOX::Abstract::Vector & LOCA::Extended::Vector::random(bool
useSeed=false, int seed=1)

Initialize every element of this vector with random values. ";

%feature("docstring")  LOCA::Extended::Vector::abs "NOX::Abstract::Vector & LOCA::Extended::Vector::abs(const
NOX::Abstract::Vector &y)

NOX::Abstract::Vector abs function. Compues absolute value of each
stored vector and scalar. ";

%feature("docstring")  LOCA::Extended::Vector::reciprocal "NOX::Abstract::Vector & LOCA::Extended::Vector::reciprocal(const
NOX::Abstract::Vector &y)

NOX::Abstract::Vector reciprocal function. Computes reciprocal of each
stored vector and scalar. ";

%feature("docstring")  LOCA::Extended::Vector::scale "NOX::Abstract::Vector & LOCA::Extended::Vector::scale(double gamma)

NOX::Abstract::Vector scale function. Scales each stored vector and
scalar. ";

%feature("docstring")  LOCA::Extended::Vector::scale "NOX::Abstract::Vector & LOCA::Extended::Vector::scale(const
NOX::Abstract::Vector &a)

NOX::Abstract::Vector scale function. Scales each stored vector and
scalar. ";

%feature("docstring")  LOCA::Extended::Vector::update "NOX::Abstract::Vector & LOCA::Extended::Vector::update(double alpha,
const NOX::Abstract::Vector &a, double gamma=0.0)

NOX::Abstract::Vector update function. Applies vector update to each
stored vector and scalar. ";

%feature("docstring")  LOCA::Extended::Vector::update "NOX::Abstract::Vector & LOCA::Extended::Vector::update(double alpha,
const NOX::Abstract::Vector &a, double beta, const
NOX::Abstract::Vector &b, double gamma=0.0)

NOX::Abstract::Vector update function. Applies vector update to each
stored vector and scalar. ";

%feature("docstring")  LOCA::Extended::Vector::norm "double
LOCA::Extended::Vector::norm(NormType type=TwoNorm) const

NOX::Abstract::Vector norm function. Computes norm of each stored
vector and combines to compute appropriate norm. ";

%feature("docstring")  LOCA::Extended::Vector::norm "double
LOCA::Extended::Vector::norm(const NOX::Abstract::Vector &weights)
const

NOX::Abstract::Vector weighted norm function. Computes weighted norm
of each stored vector and combines to compute appropriate norm. ";

%feature("docstring")  LOCA::Extended::Vector::innerProduct "double
LOCA::Extended::Vector::innerProduct(const NOX::Abstract::Vector &y)
const

NOX::Abstract::Vector innerProduct function. Computes inner product *
of each stored vector and combines to compute inner product. ";

%feature("docstring")  LOCA::Extended::Vector::length "int
LOCA::Extended::Vector::length() const

NOX::Abstract::Vector length function. Computes sum of lengths of
stored vectors plus number of scalars. ";

%feature("docstring")  LOCA::Extended::Vector::print "void
LOCA::Extended::Vector::print(std::ostream &stream) const

NOX::Abstract::Vector print function. For debugging purposes. ";

%feature("docstring")  LOCA::Extended::Vector::setVector "void
LOCA::Extended::Vector::setVector(int i, const NOX::Abstract::Vector
&v)

Sets the ith vector. ";

%feature("docstring")  LOCA::Extended::Vector::setVectorView "void
LOCA::Extended::Vector::setVectorView(int i, const Teuchos::RCP<
NOX::Abstract::Vector > &v)

Sets the ith vector as a view. ";

%feature("docstring")  LOCA::Extended::Vector::setScalar "void
LOCA::Extended::Vector::setScalar(int i, double s)

Sets the ith scalar. ";

%feature("docstring")  LOCA::Extended::Vector::setScalarArray "void
LOCA::Extended::Vector::setScalarArray(double *sv)

Sets the scalar array. ";

%feature("docstring")  LOCA::Extended::Vector::getVector "Teuchos::RCP< const NOX::Abstract::Vector >
LOCA::Extended::Vector::getVector(int i) const

Returns const ref-count pointer to the ith vector. ";

%feature("docstring")  LOCA::Extended::Vector::getVector "Teuchos::RCP< NOX::Abstract::Vector >
LOCA::Extended::Vector::getVector(int i)

Returns ref-count pointer to the ith vector. ";

%feature("docstring")  LOCA::Extended::Vector::getScalar "double
LOCA::Extended::Vector::getScalar(int i) const

Returns copy of the ith scalar. ";

%feature("docstring")  LOCA::Extended::Vector::getScalar "double &
LOCA::Extended::Vector::getScalar(int i)

Returns reference to the ith scalar. ";

%feature("docstring")  LOCA::Extended::Vector::getScalars "Teuchos::RCP< const NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Extended::Vector::getScalars() const

Returns array of scalars. ";

%feature("docstring")  LOCA::Extended::Vector::getScalars "Teuchos::RCP< NOX::Abstract::MultiVector::DenseMatrix >
LOCA::Extended::Vector::getScalars()

Returns array of scalars. ";

%feature("docstring")  LOCA::Extended::Vector::getNumScalars "int
LOCA::Extended::Vector::getNumScalars() const

Returns number of scalars. ";

%feature("docstring")  LOCA::Extended::Vector::getNumVectors "int
LOCA::Extended::Vector::getNumVectors() const

Returns number of vectors. ";


// File: classLOCA_1_1Solver_1_1Wrapper.xml
%feature("docstring") LOCA::Solver::Wrapper "

A wrapper class for wrapping a NOX solver.

The LOCA::Solver::Wrapper class provides a wrapper for NOX solvers to
change the group data returned by getSolutionGroup() and
getPreviousSolutionGroup() so that status tests can operate correctly.
LOCA continuation and bifurcation algorithms are implemented by
extended groups which augment the nonlinear equations defining
equilibrium solutions with appropriate continuation/bifurcation
equations. Therefore the groups returned by getSolutionGroup() and
getPreviousSolutionGroup() will be groups corresponding to the
continuation/bifurcation equations, not the underlying group. Status
tests that are designed to use concrete data from the original
underlying group (e.g., NOX::StatusTest::NormWRMS) will then fail
(usually via a segmentation fault or raised exception) when applied to
the extended continuation or bifurcation groups.

This solver wrapper class fixes this problem by reimplementing the
solution group accessor methods to return the underlying groups of the
solution groups (via
LOCA::Extended::MultiAbstractGroup::getUnderlyingGroup()) if they are
extended groups (derived from LOCA::Extended::MultiAbstractGroup). If
the groups are not extended groups, the original solution groups are
returned. All other NOX::Solver::Generic methods are passed to the
wrapped solver.

The LOCA::StatusTest::Wrapper class uses this wrapper class to wrap
the solver supplied via the checkStatus method which is then forwarded
to the original status test. When used properly, the group \"seen\" by
the status test is of the appropriate type for the status test.

C++ includes: LOCA_Solver_Wrapper.H ";

%feature("docstring")  LOCA::Solver::Wrapper::Wrapper "LOCA::Solver::Wrapper::Wrapper(const Teuchos::RCP<
NOX::Solver::Generic > &solver)

Constructor with a non-const ref-count pointer to a NOX solver.

The constructor calls resetWrapper() to grab the proper solution
groups. ";

%feature("docstring")  LOCA::Solver::Wrapper::Wrapper "LOCA::Solver::Wrapper::Wrapper(const Teuchos::RCP< const
NOX::Solver::Generic > &solver)

Constructor with a const ref-count pointer to a NOX solver.

The constructor calls resetWrapper() to grab the proper solution
groups. ";

%feature("docstring")  LOCA::Solver::Wrapper::~Wrapper "LOCA::Solver::Wrapper::~Wrapper()

Destructor. ";

%feature("docstring")  LOCA::Solver::Wrapper::reset "void
LOCA::Solver::Wrapper::reset(const NOX::Abstract::Vector
&initialGuess)

Implementation of reset method (forwarded to wrapped solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::reset "void
LOCA::Solver::Wrapper::reset(const NOX::Abstract::Vector
&initialGuess, const Teuchos::RCP< NOX::StatusTest::Generic > &tests)

Implementation of reset method (forwarded to wrapped solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::getStatus "NOX::StatusTest::StatusType LOCA::Solver::Wrapper::getStatus()

Implementation of getStatus method (forwarded to wrapped solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::step "NOX::StatusTest::StatusType LOCA::Solver::Wrapper::step()

Implementation of step method (forwarded to wrapped solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::solve "NOX::StatusTest::StatusType LOCA::Solver::Wrapper::solve()

Implementation of solve method (forwarded to wrapped solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::getSolutionGroup "const
NOX::Abstract::Group & LOCA::Solver::Wrapper::getSolutionGroup() const

Returns underlying group if solution group is extended, solution group
otherwise. ";

%feature("docstring")  LOCA::Solver::Wrapper::getPreviousSolutionGroup
"const NOX::Abstract::Group &
LOCA::Solver::Wrapper::getPreviousSolutionGroup() const

Returns underlying group if previous solution group is extended,
previous solution group otherwise. ";

%feature("docstring")  LOCA::Solver::Wrapper::getNumIterations "int
LOCA::Solver::Wrapper::getNumIterations() const

Implementation of getNumIterations method (forwarded to wrapped
solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::getList "const
Teuchos::ParameterList & LOCA::Solver::Wrapper::getList() const

Implementation of getList method (forwarded to wrapped solver) ";

%feature("docstring")  LOCA::Solver::Wrapper::getSolutionGroupPtr "Teuchos::RCP< const NOX::Abstract::Group >
LOCA::Solver::Wrapper::getSolutionGroupPtr() const

Returns underlying group if solution group is extended, solution group
otherwise. ";

%feature("docstring")
LOCA::Solver::Wrapper::getPreviousSolutionGroupPtr "Teuchos::RCP<
const NOX::Abstract::Group >
LOCA::Solver::Wrapper::getPreviousSolutionGroupPtr() const

Returns underlying group if previous solution group is extended,
previous solution group otherwise. ";

%feature("docstring")  LOCA::Solver::Wrapper::getListPtr "Teuchos::RCP< const Teuchos::ParameterList >
LOCA::Solver::Wrapper::getListPtr() const

Implementation of getListPtr method (forwarded to wrapped solver) ";


// File: classLOCA_1_1StatusTest_1_1Wrapper.xml
%feature("docstring") LOCA::StatusTest::Wrapper "

A wrapper class for wrapping a NOX status test.

The LOCA::StatusTest::Wrapper class provides a wrapper for NOX status
tests to change the solver passed to the wrapped status test. The
solver passed through the checkStatus() method is wrapped via the
LOCA::Solver::Wrapper class and then forwarded to the checkStatus()
method of the wrapped status test. The purpose of this is to allow
status tests that use concrete group data to function correctly when
the group is stored in an extended continuation or bifurcation group.
(See LOCA::Solver::Wrapper for more details or the LOCA status tests
page for examples on how to effectively use this class.)

C++ includes: LOCA_StatusTest_Wrapper.H ";

%feature("docstring")  LOCA::StatusTest::Wrapper::Wrapper "LOCA::StatusTest::Wrapper::Wrapper(const Teuchos::RCP<
NOX::StatusTest::Generic > &s)

Constructor. ";

%feature("docstring")  LOCA::StatusTest::Wrapper::~Wrapper "LOCA::StatusTest::Wrapper::~Wrapper()

Destructor. ";

%feature("docstring")  LOCA::StatusTest::Wrapper::checkStatus "NOX::StatusTest::StatusType
LOCA::StatusTest::Wrapper::checkStatus(const NOX::Solver::Generic
&problem, NOX::StatusTest::CheckType checkType)

Calls checkStatus of underlying status test. ";

%feature("docstring")  LOCA::StatusTest::Wrapper::getStatus "NOX::StatusTest::StatusType LOCA::StatusTest::Wrapper::getStatus()
const

Calls getStatus of underlying status test. ";

%feature("docstring")  LOCA::StatusTest::Wrapper::print "std::ostream
& LOCA::StatusTest::Wrapper::print(std::ostream &stream, int indent=0)
const

Calls print of underlying status test. ";

%feature("docstring")
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest "Teuchos::RCP<
NOX::StatusTest::Generic >
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest()

Returns underlying status test. ";

%feature("docstring")
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest "Teuchos::RCP<
const NOX::StatusTest::Generic >
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest() const

Returns underlying status test. ";


// File: classLOCA_1_1Epetra_1_1xyztPrec.xml
%feature("docstring") LOCA::Epetra::xyztPrec "

Preconditioner operator class for solving space-time (XYZT) systems.

Implements right preconditioning operators for use in global XYZT
Jacobian matrix solves.

Global - applies a right preconditioner to the global XYZT Jacobian
matrix

Sequential - applies single block right preconditioning sequentially
in time. This preconditioner is intended as an efficient competitor to
the Global preconditioner by preconditioning using only the nonzero
blocks.

Parallel - simultaneously applies sequential right preconditioning
across the decoupled time domains. This means there is no
communication of solutions between time doamins.

BlockDiagonal - similar to the Parallel preconditioner, simultaneously
applies sequential right preconditioning across the decoupled time
domains, but only using the diagonal blocks of the Jacobian matrix.
Note that the BlockDiagonal and Parallel preconditioners are
equivalent when each time domain contains only one time step.

Parareal - two pass right preconditioning applying Sequential
preconditioner over a coarse time grid (first time steps of each time
domain) and then Parallel preconditioning across the decoupled time
domains. This can be thought of as a linearized parareal strategy for
acheiving parallelism in time. The benefit over Parallel
preconditioning alone is that an estimate of the solution from the
time step on the previous domain is computed to help accelerate
convergence.

BDSDT (block diagonal in space, diagonal in time)

C++ includes: LOCA_Epetra_xyztPrec.H ";

%feature("docstring")  LOCA::Epetra::xyztPrec::xyztPrec "LOCA::Epetra::xyztPrec::xyztPrec(EpetraExt::BlockCrsMatrix &jacobian,
Epetra_CrsMatrix &splitJac, EpetraExt::BlockVector &solution,
EpetraExt::BlockVector &solutionOverlap, Epetra_Import
&overlapImporter, Teuchos::ParameterList &precPrintParams,
Teuchos::ParameterList &precLSParams, const Teuchos::RCP<
EpetraExt::MultiComm > globalComm_)

Constructor.

Builds a preconditioner operator for a full XYZT Jacobian matrix
jacobian. Right preconditioner applies are controlled using the
parameters in precLSParams. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::~xyztPrec "LOCA::Epetra::xyztPrec::~xyztPrec()

Destructor. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::SetUseTranspose "int
LOCA::Epetra::xyztPrec::SetUseTranspose(bool UseTranspose)

Set transpose. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::Apply "int
LOCA::Epetra::xyztPrec::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply XYZT preconditioner operator. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::ApplyInverse "int
LOCA::Epetra::xyztPrec::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Apply XYZT preconditioner operator inverse. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::NormInf "double
LOCA::Epetra::xyztPrec::NormInf() const

Computing infinity norm. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::Label "const char *
LOCA::Epetra::xyztPrec::Label() const

Label. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::UseTranspose "bool
LOCA::Epetra::xyztPrec::UseTranspose() const

Transpose. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::HasNormInf "bool
LOCA::Epetra::xyztPrec::HasNormInf() const

Have norm-inf. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::Comm "const
Epetra_Comm & LOCA::Epetra::xyztPrec::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::OperatorDomainMap "const Epetra_Map & LOCA::Epetra::xyztPrec::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::OperatorRangeMap "const Epetra_Map & LOCA::Epetra::xyztPrec::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::computeF "bool
LOCA::Epetra::xyztPrec::computeF(const Epetra_Vector &, Epetra_Vector
&, const NOX::Epetra::Interface::Required::FillType)

Compute residual $F$. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::computeJacobian "bool
LOCA::Epetra::xyztPrec::computeJacobian(const Epetra_Vector &,
Epetra_Operator &)

Compute Jacobian $J$. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::computePreconditioner "bool LOCA::Epetra::xyztPrec::computePreconditioner(const Epetra_Vector
&x, Epetra_Operator &Prec, Teuchos::ParameterList *p=0)

Compute preconditioner $M$. ";

%feature("docstring")  LOCA::Epetra::xyztPrec::throwError "void
LOCA::Epetra::xyztPrec::throwError(const std::string &functionName,
const std::string &errorMsg) const

Exception handler for the XYZT preconditioner class. ";


// File: namespaceEpetraExt.xml


// File: namespaceLOCA.xml
%feature("docstring")  LOCA::Abstract::createGlobalData "Teuchos::RCP< LOCA::GlobalData > LOCA::createGlobalData(const
Teuchos::RCP< Teuchos::ParameterList > &paramList, const Teuchos::RCP<
LOCA::Abstract::Factory > &userFactory=Teuchos::null)

Creates and initializes a LOCA::GlobalData object. ";

%feature("docstring")  LOCA::Abstract::destroyGlobalData "void
LOCA::destroyGlobalData(const Teuchos::RCP< LOCA::GlobalData >
&globalData)

De-initializes a LOCA::GlobalData object for destruction.

Sets the data members to Teuchos::null to remove circular references
";


// File: namespaceLOCA_1_1Abstract.xml


// File: namespaceLOCA_1_1AnasaziOperator.xml


// File: namespaceLOCA_1_1Bifurcation.xml


// File: namespaceLOCA_1_1Bifurcation_1_1PitchforkBord.xml


// File: namespaceLOCA_1_1Bifurcation_1_1PitchforkBord_1_1StatusTest.xml


// File: namespaceLOCA_1_1Bifurcation_1_1TPBord.xml


// File: namespaceLOCA_1_1Bifurcation_1_1TPBord_1_1StatusTest.xml


// File: namespaceLOCA_1_1BorderedSolver.xml


// File: namespaceLOCA_1_1BorderedSystem.xml


// File: namespaceLOCA_1_1Continuation.xml


// File: namespaceLOCA_1_1Continuation_1_1StatusTest.xml


// File: namespaceLOCA_1_1Eigensolver.xml


// File: namespaceLOCA_1_1EigenvalueSort.xml


// File: namespaceLOCA_1_1Epetra.xml


// File: namespaceLOCA_1_1Epetra_1_1AnasaziOperator.xml


// File: namespaceLOCA_1_1Epetra_1_1Interface.xml


// File: namespaceLOCA_1_1Epetra_1_1TransposeLinearSystem.xml


// File: namespaceLOCA_1_1Extended.xml


// File: namespaceLOCA_1_1Homotopy.xml


// File: namespaceLOCA_1_1Hopf.xml


// File: namespaceLOCA_1_1Hopf_1_1MinimallyAugmented.xml


// File: namespaceLOCA_1_1Hopf_1_1MooreSpence.xml


// File: namespaceLOCA_1_1MultiContinuation.xml


// File: namespaceLOCA_1_1MultiPredictor.xml


// File: namespaceLOCA_1_1Parameter.xml


// File: namespaceLOCA_1_1PhaseTransition.xml


// File: namespaceLOCA_1_1Pitchfork.xml


// File: namespaceLOCA_1_1Pitchfork_1_1MinimallyAugmented.xml


// File: namespaceLOCA_1_1Pitchfork_1_1MooreSpence.xml


// File: namespaceLOCA_1_1SaveEigenData.xml


// File: namespaceLOCA_1_1SingularJacobianSolve.xml


// File: namespaceLOCA_1_1Solver.xml


// File: namespaceLOCA_1_1StatusTest.xml


// File: namespaceLOCA_1_1StepSize.xml


// File: namespaceLOCA_1_1TimeDependent.xml


// File: namespaceLOCA_1_1TurningPoint.xml


// File: namespaceLOCA_1_1TurningPoint_1_1MinimallyAugmented.xml


// File: namespaceLOCA_1_1TurningPoint_1_1MooreSpence.xml


// File: namespaceNOX.xml


// File: namespaceNOX_1_1Epetra.xml


// File: namespaceNOX_1_1Parameter.xml


// File: namespaceNOX_1_1Solver.xml


// File: namespaceTeuchos.xml


// File: LOCA_8H.xml


// File: LOCA__Abstract__Factory_8C.xml


// File: LOCA__Abstract__Factory_8H.xml


// File: LOCA__Abstract__Group_8C.xml


// File: LOCA__Abstract__Group_8H.xml


// File: LOCA__Abstract__Iterator_8C.xml


// File: LOCA__Abstract__Iterator_8H.xml


// File: LOCA__Abstract__TransposeSolveGroup_8H.xml


// File: LOCA__AnasaziOperator__AbstractStrategy_8H.xml


// File: LOCA__AnasaziOperator__Cayley_8C.xml


// File: LOCA__AnasaziOperator__Cayley_8H.xml


// File: LOCA__AnasaziOperator__Cayley2Matrix_8C.xml


// File: LOCA__AnasaziOperator__Cayley2Matrix_8H.xml


// File: LOCA__AnasaziOperator__Factory_8C.xml


// File: LOCA__AnasaziOperator__Factory_8H.xml


// File: LOCA__AnasaziOperator__JacobianInverse_8C.xml


// File: LOCA__AnasaziOperator__JacobianInverse_8H.xml


// File: LOCA__AnasaziOperator__ShiftInvert_8C.xml


// File: LOCA__AnasaziOperator__ShiftInvert_8H.xml


// File: LOCA__AnasaziOperator__ShiftInvert2Matrix_8C.xml


// File: LOCA__AnasaziOperator__ShiftInvert2Matrix_8H.xml


// File: LOCA__Bifurcation__Factory_8C.xml


// File: LOCA__Bifurcation__Factory_8H.xml


// File: LOCA__Bifurcation__PitchforkBord__NullVectorNormWRMS_8C.xml


// File: LOCA__Bifurcation__PitchforkBord__NullVectorNormWRMS_8H.xml


// File: LOCA__Bifurcation__PitchforkBord__ParameterUpdateNorm_8C.xml


// File: LOCA__Bifurcation__PitchforkBord__ParameterUpdateNorm_8H.xml


// File: LOCA__Bifurcation__PitchforkBord__SlackUpdateNorm_8C.xml


// File: LOCA__Bifurcation__PitchforkBord__SlackUpdateNorm_8H.xml


// File: LOCA__Bifurcation__TPBord__StatusTest__NullVectorNormWRMS_8C.xml


// File: LOCA__Bifurcation__TPBord__StatusTest__NullVectorNormWRMS_8H.xml


// File: LOCA__Bifurcation__TPBord__StatusTest__ParameterUpdateNorm_8C.xml


// File: LOCA__Bifurcation__TPBord__StatusTest__ParameterUpdateNorm_8H.xml


// File: LOCA__BorderedSolver__AbstractOperator_8H.xml


// File: LOCA__BorderedSolver__AbstractStrategy_8C.xml


// File: LOCA__BorderedSolver__AbstractStrategy_8H.xml


// File: LOCA__BorderedSolver__BorderedOperator_8H.xml


// File: LOCA__BorderedSolver__Bordering_8C.xml


// File: LOCA__BorderedSolver__Bordering_8H.xml


// File: LOCA__BorderedSolver__ComplexOperator_8C.xml


// File: LOCA__BorderedSolver__ComplexOperator_8H.xml


// File: LOCA__BorderedSolver__EpetraAugmented_8C.xml


// File: LOCA__BorderedSolver__EpetraAugmented_8H.xml


// File: LOCA__BorderedSolver__EpetraHouseholder_8C.xml


// File: LOCA__BorderedSolver__EpetraHouseholder_8H.xml


// File: LOCA__BorderedSolver__Factory_8C.xml


// File: LOCA__BorderedSolver__Factory_8H.xml


// File: LOCA__BorderedSolver__HouseholderQR_8C.xml


// File: LOCA__BorderedSolver__HouseholderQR_8H.xml


// File: LOCA__BorderedSolver__JacobianOperator_8C.xml


// File: LOCA__BorderedSolver__JacobianOperator_8H.xml


// File: LOCA__BorderedSolver__LowerTriangularBlockElimination_8C.xml


// File: LOCA__BorderedSolver__LowerTriangularBlockElimination_8H.xml


// File: LOCA__BorderedSolver__Nested_8C.xml


// File: LOCA__BorderedSolver__Nested_8H.xml


// File: LOCA__BorderedSolver__UpperTriangularBlockElimination_8C.xml


// File: LOCA__BorderedSolver__UpperTriangularBlockElimination_8H.xml


// File: LOCA__BorderedSystem__AbstractGroup_8H.xml


// File: LOCA__Continuation__StatusTest__ParameterResidualNorm_8C.xml


// File: LOCA__Continuation__StatusTest__ParameterResidualNorm_8H.xml


// File: LOCA__Continuation__StatusTest__ParameterUpdateNorm_8C.xml


// File: LOCA__Continuation__StatusTest__ParameterUpdateNorm_8H.xml


// File: LOCA__DerivUtils_8C.xml


// File: LOCA__DerivUtils_8H.xml


// File: LOCA__Description_8H.xml


// File: LOCA__Eigensolver__AbstractStrategy_8H.xml


// File: LOCA__Eigensolver__AnasaziStrategy_8C.xml


// File: LOCA__Eigensolver__AnasaziStrategy_8H.xml


// File: LOCA__Eigensolver__DefaultStrategy_8C.xml


// File: LOCA__Eigensolver__DefaultStrategy_8H.xml


// File: LOCA__Eigensolver__Factory_8C.xml


// File: LOCA__Eigensolver__Factory_8H.xml


// File: LOCA__EigenvalueSort__Factory_8C.xml


// File: LOCA__EigenvalueSort__Factory_8H.xml


// File: LOCA__EigenvalueSort__Strategies_8C.xml


// File: LOCA__EigenvalueSort__Strategies_8H.xml


// File: LOCA__Epetra_8H.xml


// File: LOCA__Epetra__AdaptiveSolutionManager_8C.xml


// File: LOCA__Epetra__AdaptiveSolutionManager_8H.xml


// File: LOCA__Epetra__AdaptiveStepper_8C.xml


// File: LOCA__Epetra__AdaptiveStepper_8H.xml


// File: LOCA__Epetra__AnasaziOperator__Floquet_8C.xml


// File: LOCA__Epetra__AnasaziOperator__Floquet_8H.xml


// File: LOCA__Epetra__AugmentedOp_8C.xml


// File: LOCA__Epetra__AugmentedOp_8H.xml


// File: LOCA__Epetra__CompactWYOp_8C.xml


// File: LOCA__Epetra__CompactWYOp_8H.xml


// File: LOCA__Epetra__Factory_8C.xml


// File: LOCA__Epetra__Factory_8H.xml


// File: LOCA__Epetra__Group_8C.xml


// File: LOCA__Epetra__Group_8H.xml


// File: LOCA__Epetra__IdentityOp_8C.xml


// File: LOCA__Epetra__IdentityOp_8H.xml


// File: LOCA__Epetra__Interface__FreeEnergy_8H.xml


// File: LOCA__Epetra__Interface__MassMatrix_8H.xml


// File: LOCA__Epetra__Interface__MultiPoint_8C.xml


// File: LOCA__Epetra__Interface__MultiPoint_8H.xml


// File: LOCA__Epetra__Interface__Required_8H.xml


// File: LOCA__Epetra__Interface__TimeDependent_8H.xml


// File: LOCA__Epetra__Interface__TimeDependentMatrixFree_8H.xml


// File: LOCA__Epetra__Interface__xyzt_8C.xml


// File: LOCA__Epetra__Interface__xyzt_8H.xml


// File: LOCA__Epetra__LeftPreconditionedOp_8C.xml


// File: LOCA__Epetra__LeftPreconditionedOp_8H.xml


// File: LOCA__Epetra__LowRankUpdateOp_8C.xml


// File: LOCA__Epetra__LowRankUpdateOp_8H.xml


// File: LOCA__Epetra__LowRankUpdateRowMatrix_8C.xml


// File: LOCA__Epetra__LowRankUpdateRowMatrix_8H.xml


// File: LOCA__Epetra__ModelEvaluatorInterface_8C.xml


// File: LOCA__Epetra__ModelEvaluatorInterface_8H.xml


// File: LOCA__Epetra__ShiftInvertOperator_8C.xml


// File: LOCA__Epetra__ShiftInvertOperator_8H.xml


// File: LOCA__Epetra__TransposeLinearSystem__AbstractStrategy_8H.xml


// File: LOCA__Epetra__TransposeLinearSystem__ExplicitTranspose_8C.xml


// File: LOCA__Epetra__TransposeLinearSystem__ExplicitTranspose_8H.xml


// File: LOCA__Epetra__TransposeLinearSystem__Factory_8C.xml


// File: LOCA__Epetra__TransposeLinearSystem__Factory_8H.xml


// File: LOCA__Epetra__TransposeLinearSystem__LeftPreconditioning_8C.xml


// File: LOCA__Epetra__TransposeLinearSystem__LeftPreconditioning_8H.xml


// File: LOCA__Epetra__TransposeLinearSystem__TransposePreconditioner_8C.xml


// File: LOCA__Epetra__TransposeLinearSystem__TransposePreconditioner_8H.xml


// File: LOCA__Epetra__xyztPrec_8C.xml


// File: LOCA__Epetra__xyztPrec_8H.xml


// File: LOCA__ErrorCheck_8C.xml


// File: LOCA__ErrorCheck_8H.xml


// File: LOCA__Extended__MultiAbstractGroup_8C.xml


// File: LOCA__Extended__MultiAbstractGroup_8H.xml


// File: LOCA__Extended__MultiVector_8C.xml


// File: LOCA__Extended__MultiVector_8H.xml


// File: LOCA__Extended__Vector_8C.xml


// File: LOCA__Extended__Vector_8H.xml


// File: LOCA__Factory_8C.xml


// File: LOCA__Factory_8H.xml


// File: LOCA__GlobalData_8C.xml


// File: LOCA__GlobalData_8H.xml


// File: LOCA__Homotopy__AbstractGroup_8H.xml


// File: LOCA__Homotopy__DeflatedGroup_8C.xml


// File: LOCA__Homotopy__DeflatedGroup_8H.xml


// File: LOCA__Homotopy__Group_8C.xml


// File: LOCA__Homotopy__Group_8H.xml


// File: LOCA__Hopf__ComplexMultiVector_8C.xml


// File: LOCA__Hopf__ComplexMultiVector_8H.xml


// File: LOCA__Hopf__ComplexVector_8C.xml


// File: LOCA__Hopf__ComplexVector_8H.xml


// File: LOCA__Hopf__MinimallyAugmented__AbstractGroup_8H.xml


// File: LOCA__Hopf__MinimallyAugmented__Constraint_8C.xml


// File: LOCA__Hopf__MinimallyAugmented__Constraint_8H.xml


// File: LOCA__Hopf__MinimallyAugmented__ExtendedGroup_8C.xml


// File: LOCA__Hopf__MinimallyAugmented__ExtendedGroup_8H.xml


// File: LOCA__Hopf__MinimallyAugmented__FiniteDifferenceGroup_8C.xml


// File: LOCA__Hopf__MinimallyAugmented__FiniteDifferenceGroup_8H.xml


// File: LOCA__Hopf__MooreSpence__AbstractGroup_8H.xml


// File: LOCA__Hopf__MooreSpence__ExtendedGroup_8C.xml


// File: LOCA__Hopf__MooreSpence__ExtendedGroup_8H.xml


// File: LOCA__Hopf__MooreSpence__ExtendedMultiVector_8C.xml


// File: LOCA__Hopf__MooreSpence__ExtendedMultiVector_8H.xml


// File: LOCA__Hopf__MooreSpence__ExtendedVector_8C.xml


// File: LOCA__Hopf__MooreSpence__ExtendedVector_8H.xml


// File: LOCA__Hopf__MooreSpence__FiniteDifferenceGroup_8C.xml


// File: LOCA__Hopf__MooreSpence__FiniteDifferenceGroup_8H.xml


// File: LOCA__Hopf__MooreSpence__SalingerBordering_8C.xml


// File: LOCA__Hopf__MooreSpence__SalingerBordering_8H.xml


// File: LOCA__Hopf__MooreSpence__SolverFactory_8C.xml


// File: LOCA__Hopf__MooreSpence__SolverFactory_8H.xml


// File: LOCA__Hopf__MooreSpence__SolverStrategy_8H.xml


// File: LOCA__MultiContinuation__AbstractGroup_8C.xml


// File: LOCA__MultiContinuation__AbstractGroup_8H.xml


// File: LOCA__MultiContinuation__AbstractStrategy_8H.xml


// File: LOCA__MultiContinuation__ArcLengthConstraint_8C.xml


// File: LOCA__MultiContinuation__ArcLengthConstraint_8H.xml


// File: LOCA__MultiContinuation__ArcLengthGroup_8C.xml


// File: LOCA__MultiContinuation__ArcLengthGroup_8H.xml


// File: LOCA__MultiContinuation__CompositeConstraint_8C.xml


// File: LOCA__MultiContinuation__CompositeConstraint_8H.xml


// File: LOCA__MultiContinuation__CompositeConstraintMVDX_8C.xml


// File: LOCA__MultiContinuation__CompositeConstraintMVDX_8H.xml


// File: LOCA__MultiContinuation__ConstrainedGroup_8C.xml


// File: LOCA__MultiContinuation__ConstrainedGroup_8H.xml


// File: LOCA__MultiContinuation__ConstraintInterface_8H.xml


// File: LOCA__MultiContinuation__ConstraintInterfaceMVDX_8C.xml


// File: LOCA__MultiContinuation__ConstraintInterfaceMVDX_8H.xml


// File: LOCA__MultiContinuation__ExtendedGroup_8C.xml


// File: LOCA__MultiContinuation__ExtendedGroup_8H.xml


// File: LOCA__MultiContinuation__ExtendedMultiVector_8C.xml


// File: LOCA__MultiContinuation__ExtendedMultiVector_8H.xml


// File: LOCA__MultiContinuation__ExtendedVector_8C.xml


// File: LOCA__MultiContinuation__ExtendedVector_8H.xml


// File: LOCA__MultiContinuation__Factory_8C.xml


// File: LOCA__MultiContinuation__Factory_8H.xml


// File: LOCA__MultiContinuation__FiniteDifferenceGroup_8C.xml


// File: LOCA__MultiContinuation__FiniteDifferenceGroup_8H.xml


// File: LOCA__MultiContinuation__MultiVecConstraint_8C.xml


// File: LOCA__MultiContinuation__MultiVecConstraint_8H.xml


// File: LOCA__MultiContinuation__NaturalConstraint_8C.xml


// File: LOCA__MultiContinuation__NaturalConstraint_8H.xml


// File: LOCA__MultiContinuation__NaturalGroup_8C.xml


// File: LOCA__MultiContinuation__NaturalGroup_8H.xml


// File: LOCA__MultiPredictor__AbstractStrategy_8C.xml


// File: LOCA__MultiPredictor__AbstractStrategy_8H.xml


// File: LOCA__MultiPredictor__Constant_8C.xml


// File: LOCA__MultiPredictor__Constant_8H.xml


// File: LOCA__MultiPredictor__Factory_8C.xml


// File: LOCA__MultiPredictor__Factory_8H.xml


// File: LOCA__MultiPredictor__Random_8C.xml


// File: LOCA__MultiPredictor__Random_8H.xml


// File: LOCA__MultiPredictor__Restart_8C.xml


// File: LOCA__MultiPredictor__Restart_8H.xml


// File: LOCA__MultiPredictor__Secant_8C.xml


// File: LOCA__MultiPredictor__Secant_8H.xml


// File: LOCA__MultiPredictor__Tangent_8C.xml


// File: LOCA__MultiPredictor__Tangent_8H.xml


// File: LOCA__Parameter__Entry_8H.xml


// File: LOCA__Parameter__Library_8C.xml


// File: LOCA__Parameter__Library_8H.xml


// File: LOCA__Parameter__LibraryT_8H.xml


// File: LOCA__Parameter__SublistParser_8C.xml


// File: LOCA__Parameter__SublistParser_8H.xml


// File: LOCA__Parameter__Vector_8C.xml


// File: LOCA__Parameter__Vector_8H.xml


// File: LOCA__PhaseTransition__AbstractGroup_8H.xml


// File: LOCA__PhaseTransition__ExtendedGroup_8C.xml
%feature("docstring")  post_process "void post_process(double **,
char *, int *, double *, int, int) ";


// File: LOCA__PhaseTransition__ExtendedGroup_8H.xml


// File: LOCA__PhaseTransition__ExtendedMultiVector_8C.xml


// File: LOCA__PhaseTransition__ExtendedMultiVector_8H.xml


// File: LOCA__PhaseTransition__ExtendedVector_8C.xml


// File: LOCA__PhaseTransition__ExtendedVector_8H.xml


// File: LOCA__Pitchfork__MinimallyAugmented__AbstractGroup_8H.xml


// File: LOCA__Pitchfork__MinimallyAugmented__Constraint_8C.xml


// File: LOCA__Pitchfork__MinimallyAugmented__Constraint_8H.xml


// File: LOCA__Pitchfork__MinimallyAugmented__ExtendedGroup_8C.xml


// File: LOCA__Pitchfork__MinimallyAugmented__ExtendedGroup_8H.xml


// File: LOCA__Pitchfork__MooreSpence__AbstractGroup_8H.xml


// File: LOCA__Pitchfork__MooreSpence__ExtendedGroup_8C.xml


// File: LOCA__Pitchfork__MooreSpence__ExtendedGroup_8H.xml


// File: LOCA__Pitchfork__MooreSpence__ExtendedMultiVector_8C.xml


// File: LOCA__Pitchfork__MooreSpence__ExtendedMultiVector_8H.xml


// File: LOCA__Pitchfork__MooreSpence__ExtendedVector_8C.xml


// File: LOCA__Pitchfork__MooreSpence__ExtendedVector_8H.xml


// File: LOCA__Pitchfork__MooreSpence__PhippsBordering_8C.xml


// File: LOCA__Pitchfork__MooreSpence__PhippsBordering_8H.xml


// File: LOCA__Pitchfork__MooreSpence__SalingerBordering_8C.xml


// File: LOCA__Pitchfork__MooreSpence__SalingerBordering_8H.xml


// File: LOCA__Pitchfork__MooreSpence__SolverFactory_8C.xml


// File: LOCA__Pitchfork__MooreSpence__SolverFactory_8H.xml


// File: LOCA__Pitchfork__MooreSpence__SolverStrategy_8H.xml


// File: LOCA__SaveEigenData__AbstractStrategy_8H.xml


// File: LOCA__SaveEigenData__DefaultStrategy_8C.xml


// File: LOCA__SaveEigenData__DefaultStrategy_8H.xml


// File: LOCA__SaveEigenData__Factory_8C.xml


// File: LOCA__SaveEigenData__Factory_8H.xml


// File: LOCA__SingularJacobianSolve__Default_8C.xml


// File: LOCA__SingularJacobianSolve__Default_8H.xml


// File: LOCA__SingularJacobianSolve__Generic_8H.xml


// File: LOCA__SingularJacobianSolve__ItRef_8C.xml


// File: LOCA__SingularJacobianSolve__ItRef_8H.xml


// File: LOCA__SingularJacobianSolve__Manager_8C.xml


// File: LOCA__SingularJacobianSolve__Manager_8H.xml


// File: LOCA__SingularJacobianSolve__Nic_8C.xml


// File: LOCA__SingularJacobianSolve__Nic_8H.xml


// File: LOCA__SingularJacobianSolve__NicDay_8C.xml


// File: LOCA__SingularJacobianSolve__NicDay_8H.xml


// File: LOCA__Solver__Wrapper_8C.xml


// File: LOCA__Solver__Wrapper_8H.xml


// File: LOCA__StatusTest__Abstract_8C.xml


// File: LOCA__StatusTest__Abstract_8H.xml


// File: LOCA__StatusTest__Combo_8C.xml
%feature("docstring")  checkStatus "LOCA::StatusTest::StatusType
LOCA::StatusTest::Combo:: checkStatus(const LOCA::Abstract::Iterator
&stepper, LOCA::StatusTest::CheckType checkType) ";


// File: LOCA__StatusTest__Combo_8H.xml


// File: LOCA__StatusTest__Factory_8C.xml


// File: LOCA__StatusTest__Factory_8H.xml


// File: LOCA__StatusTest__MaxIters_8C.xml
%feature("docstring")  checkStatus "LOCA::StatusTest::StatusType
LOCA::StatusTest::MaxIters:: checkStatus(const
LOCA::Abstract::Iterator &stepper, LOCA::StatusTest::CheckType
checkType) ";


// File: LOCA__StatusTest__MaxIters_8H.xml


// File: LOCA__StatusTest__Wrapper_8C.xml


// File: LOCA__StatusTest__Wrapper_8H.xml


// File: LOCA__Stepper_8C.xml


// File: LOCA__Stepper_8H.xml


// File: LOCA__StepSize__AbstractStrategy_8H.xml


// File: LOCA__StepSize__Adaptive_8C.xml


// File: LOCA__StepSize__Adaptive_8H.xml


// File: LOCA__StepSize__Constant_8C.xml


// File: LOCA__StepSize__Constant_8H.xml


// File: LOCA__StepSize__Factory_8C.xml


// File: LOCA__StepSize__Factory_8H.xml


// File: LOCA__TimeDependent__AbstractGroup_8H.xml


// File: LOCA__TurningPoint__MinimallyAugmented__AbstractGroup_8H.xml


// File: LOCA__TurningPoint__MinimallyAugmented__Constraint_8C.xml


// File: LOCA__TurningPoint__MinimallyAugmented__Constraint_8H.xml


// File: LOCA__TurningPoint__MinimallyAugmented__ExtendedGroup_8C.xml


// File: LOCA__TurningPoint__MinimallyAugmented__ExtendedGroup_8H.xml


// File: LOCA__TurningPoint__MinimallyAugmented__FiniteDifferenceGroup_8C.xml


// File: LOCA__TurningPoint__MinimallyAugmented__FiniteDifferenceGroup_8H.xml


// File: LOCA__TurningPoint__MinimallyAugmented__ModifiedConstraint_8C.xml


// File: LOCA__TurningPoint__MinimallyAugmented__ModifiedConstraint_8H.xml


// File: LOCA__TurningPoint__MooreSpence__AbstractGroup_8H.xml


// File: LOCA__TurningPoint__MooreSpence__ExtendedGroup_8C.xml


// File: LOCA__TurningPoint__MooreSpence__ExtendedGroup_8H.xml


// File: LOCA__TurningPoint__MooreSpence__ExtendedMultiVector_8C.xml


// File: LOCA__TurningPoint__MooreSpence__ExtendedMultiVector_8H.xml


// File: LOCA__TurningPoint__MooreSpence__ExtendedVector_8C.xml


// File: LOCA__TurningPoint__MooreSpence__ExtendedVector_8H.xml


// File: LOCA__TurningPoint__MooreSpence__FiniteDifferenceGroup_8C.xml


// File: LOCA__TurningPoint__MooreSpence__FiniteDifferenceGroup_8H.xml


// File: LOCA__TurningPoint__MooreSpence__PhippsBordering_8C.xml


// File: LOCA__TurningPoint__MooreSpence__PhippsBordering_8H.xml


// File: LOCA__TurningPoint__MooreSpence__SalingerBordering_8C.xml


// File: LOCA__TurningPoint__MooreSpence__SalingerBordering_8H.xml


// File: LOCA__TurningPoint__MooreSpence__SolverFactory_8C.xml


// File: LOCA__TurningPoint__MooreSpence__SolverFactory_8H.xml


// File: LOCA__TurningPoint__MooreSpence__SolverStrategy_8H.xml


// File: loca_overview.xml


// File: loca_class_overview.xml


// File: loca_parameters.xml


// File: loca_user_info.xml


// File: loca_continuation_tutorial.xml


// File: loca_tp_continuation_tutorial.xml


// File: new_loca_framework.xml


// File: deprecated.xml


// File: dir_c9905962d2938ddd59412d843d73cc24.xml


// File: dir_bfeef42343a264a4e64fa04f2f91485c.xml


// File: dir_b94a64db05dcca28bd6c00c353b931ee.xml


// File: dir_2b7c7ba56e6f30a4947b8163393ea120.xml


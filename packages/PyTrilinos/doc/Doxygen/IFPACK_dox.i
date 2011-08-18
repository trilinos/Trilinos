
// File: index.xml

// File: structcs__symbolic.xml
%feature("docstring") cs_symbolic "";


// File: structcsr__dmperm__results.xml
%feature("docstring") csr_dmperm_results "";


// File: structcsr__numeric.xml
%feature("docstring") csr_numeric "";


// File: classIfpack.xml
%feature("docstring") Ifpack "

Ifpack: a function class to define Ifpack preconditioners.

Class Ifpack is a function class, that contains just one method:
Create(). Using Create(), users can easily define a variety of IFPACK
preconditioners.

Create requires 3 arguments: a string, indicating the preconditioner
to be built;

a pointer to an Epetra_RowMatrix, representing the matrix to be used
to define the preconditioner;

an interger (defaulted to 0), that specifies the amount of overlap
among the processes.

The first argument can assume the following values:  \"point
relaxation\" : returns an instance of
Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>

\"point relaxation stand-alone\" : returns an instance of
Ifpack_PointRelaxation (value of overlap is ignored).

\"block relaxation\" : returns an instance of
Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation>

\"block relaxation stand-alone)\" : returns an instance of
Ifpack_BlockRelaxation.

\"Amesos\" : returns an instance of
Ifpack_AdditiveSchwarz<Ifpack_Amesos>.

\"Amesos stand-alone\" : returns an instance of Ifpack_Amesos.

\"IC\" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_IC>.

\"IC stand-alone\" : returns an instance of
Ifpack_AdditiveSchwarz<Ifpack_IC>.

\"ICT\" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ICT>.

\"ICT stand-alone\" : returns an instance of Ifpack_ICT.

\"ILU\" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILU>.

\"ILU stand-alone\" : returns an instance of Ifpack_ILU.

\"ILUT\" : returns an instance of Ifpack_AdditiveSchwarz<Ifpack_ILUT>.

\"ILUT stand-alone\" : returns an instance of Ifpack_ILUT.

otherwise, Create() returns 0.

Objects in stand-alone mode cannot use reordering, variable overlap,
and singleton filters. However, their construction can be slightly
faster than the non stand-alone counterpart.

The following fragment of code shows the basic usage of this class.

Marzio Sala, (formally) SNL org. 1414

C++ includes: Ifpack.h ";

%feature("docstring")  Ifpack::Create "Ifpack_Preconditioner *
Ifpack::Create(const string PrecType, Epetra_RowMatrix *Matrix, const
int overlap=0)

Creates an instance of Ifpack_Preconditioner given the string name of
the preconditioner type (can fail with bad input).

Parameters:
-----------

PrecType:  (In) - String name of preconditioner type to be created.

Matrix:  (In) - Matrix used to define the preconditioner

overlap:  (In) - specified overlap, defaulted to 0.

Returns 0 if the preconditioner with that input name does not exist.
Otherwise, return a newly created preconditioner object. Note that the
client is responsible for calling delete on the returned object once
it is finished using it! ";

%feature("docstring")  Ifpack::SetParameters "int
Ifpack::SetParameters(int argc, char *argv[], Teuchos::ParameterList
&List, string &PrecType, int &Overlap)

Sets the options in List from the command line.

Note: If you want full support for all parameters, consider reading in
a parameter list from an XML file as supported by the Teuchos helper
function Teuchos::updateParametersFromXmlFile() or
Teuchos::updateParametersFromXmlStream(). ";


// File: classIfpack__AbsComp.xml
%feature("docstring") Ifpack_AbsComp "";


// File: classIfpack__AdditiveSchwarz.xml
%feature("docstring") Ifpack_AdditiveSchwarz "

Ifpack_AdditiveSchwarz: a class to define Additive Schwarz
preconditioners of Epetra_RowMatrix's.

Class Ifpack_AdditiveSchwarz enables the construction of Additive
Schwarz (one-level overlapping domain decomposition) preconditioners,
for a given Epetra_RowMatrix. Ifpack_AdditiveSchwarz is derived from
Ifpack_Preconditioner, itself derived from Epetra_Operator. An
application of the Additive Schwarz preconditioner can be obtained by
calling method ApplyInverse().

One-level overlapping domain decomposition preconditioners use local
solvers, of Dirichlet type. This means that the inverse of the local
matrix (with minimal or wider overlap) is applied to the residual to
be preconditioned.

The preconditioner can be written as: \\\\[ P_{AS}^{-1} =
\\\\sum_{i=1}^M P_i A_i^{-1} R_i , \\\\] where $M$ is the number of
subdomains (that is, the number of processors in the computation),
$R_i$ is an operator that restricts the global vector to the vector
lying on subdomain $i$, $P_i$ is the prolongator operator, and \\\\[
A_i = R_i A P_i. \\\\]

The construction of Schwarz preconditioners is mainly composed by two
steps: definition of the restriction and prolongation operator $R_i$
and $R_i^T$. If minimal overlap is chosen, their implementation is
trivial, $R_i$ will return all the local components. For wider
overlaps, instead, Epetra_Import and Epetra_Export will be used to
import/export data. The user must provide both the matrix to be
preconditioned (which is suppose to have minimal-overlap) and the
matrix with wider overlap.

definition of a technique to apply the inverse of $A_i$. To solve on
each subdomain, the user can adopt any class, derived from
Ifpack_Preconditioner. This can be easily accomplished, as
Ifpack_AdditiveSchwarz is templated with the solver for each
subdomain.

The local matrix $A_i$ can be filtered, to eliminate singletons, and
reordered. At the present time, RCM and METIS can be used to reorder
the local matrix.

The complete list of supported parameters is reported in page
ifp_params.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_AdditiveSchwarz.h ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Ifpack_AdditiveSchwarz
"Ifpack_AdditiveSchwarz< T >::Ifpack_AdditiveSchwarz(Epetra_RowMatrix
*Matrix_in, int OverlapLevel_in=0)

Ifpack_AdditiveSchwarz constructor with given Epetra_RowMatrix.

Creates an Ifpack_AdditiveSchwarz preconditioner with overlap. To use
minimal-overlap, OverlappingMatrix is omitted (as defaulted to 0).

Parameters:
-----------

Matrix:  - (In) Pointer to matrix to be preconditioned

OverlappingMatrix:  - (In) Pointer to the matrix extended with the
desired level of overlap. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::~Ifpack_AdditiveSchwarz
"virtual Ifpack_AdditiveSchwarz< T >::~Ifpack_AdditiveSchwarz()

Destructor. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::SetUseTranspose "int
Ifpack_AdditiveSchwarz< T >::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied (not
implemented).

This flag allows the transpose of the given operator to be used
implicitly.

Parameters:
-----------

UseTranspose_in:  - (In) If true, multiply by the transpose of
operator, otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Apply "int
Ifpack_AdditiveSchwarz< T >::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the matrix to X, returns the result in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  -(Out) A Epetra_MultiVector of dimension NumVectors containing the
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::ApplyInverse "int
Ifpack_AdditiveSchwarz< T >::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the preconditioner to X, returns the result in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to be
preconditioned.

Y:  -(Out) A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

WARNING:  In order to work with AztecOO, any implementation of this
method must support the case where X and Y are the same object. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::NormInf "double
Ifpack_AdditiveSchwarz< T >::NormInf() const

Returns the infinity norm of the global matrix (not implemented) ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Label "const char *
Ifpack_AdditiveSchwarz< T >::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::UseTranspose "bool
Ifpack_AdditiveSchwarz< T >::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::HasNormInf "bool
Ifpack_AdditiveSchwarz< T >::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Comm "const
Epetra_Comm & Ifpack_AdditiveSchwarz< T >::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::OperatorDomainMap "const Epetra_Map & Ifpack_AdditiveSchwarz< T >::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::OperatorRangeMap "const Epetra_Map & Ifpack_AdditiveSchwarz< T >::OperatorRangeMap()
const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Initialize "int
Ifpack_AdditiveSchwarz< T >::Initialize()

Initialized the preconditioner. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Compute "int
Ifpack_AdditiveSchwarz< T >::Compute()

Computes the preconditioner. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Condest "double
Ifpack_AdditiveSchwarz< T >::Condest(const Ifpack_CondestType
CT=Ifpack_Cheap, const int MaxIters=1550, const double Tol=1e-9,
Epetra_RowMatrix *Matrix_in=0)

Computes the estimated condition number and returns its value. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Condest "virtual
double Ifpack_AdditiveSchwarz< T >::Condest() const

Returns the estimated condition number, or -1.0 if not computed. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Matrix "virtual const
Epetra_RowMatrix& Ifpack_AdditiveSchwarz< T >::Matrix() const

Returns a refernence to the internally stored matrix. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::IsOverlapping "virtual
bool Ifpack_AdditiveSchwarz< T >::IsOverlapping() const

Returns true is an overlapping matrix is present. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Print "std::ostream &
Ifpack_AdditiveSchwarz< T >::Print(std::ostream &) const

Prints major information about this preconditioner. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::Inverse "virtual const
T* Ifpack_AdditiveSchwarz< T >::Inverse() const ";

%feature("docstring")  Ifpack_AdditiveSchwarz::NumInitialize "virtual
int Ifpack_AdditiveSchwarz< T >::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_AdditiveSchwarz::NumCompute "virtual
int Ifpack_AdditiveSchwarz< T >::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_AdditiveSchwarz::NumApplyInverse "virtual int Ifpack_AdditiveSchwarz< T >::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_AdditiveSchwarz::InitializeTime "virtual double Ifpack_AdditiveSchwarz< T >::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_AdditiveSchwarz::ComputeTime "virtual
double Ifpack_AdditiveSchwarz< T >::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_AdditiveSchwarz::ApplyInverseTime "virtual double Ifpack_AdditiveSchwarz< T >::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_AdditiveSchwarz::InitializeFlops "virtual double Ifpack_AdditiveSchwarz< T >::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::ComputeFlops "virtual
double Ifpack_AdditiveSchwarz< T >::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::ApplyInverseFlops "virtual double Ifpack_AdditiveSchwarz< T >::ApplyInverseFlops() const

Returns the number of flops in the application of the preconditioner.
";

%feature("docstring")  Ifpack_AdditiveSchwarz::OverlapLevel "virtual
int Ifpack_AdditiveSchwarz< T >::OverlapLevel() const

Returns the level of overlap. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::List "virtual const
Teuchos::ParameterList& Ifpack_AdditiveSchwarz< T >::List() const

Returns a reference to the internally stored list. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::IsInitialized "virtual
bool Ifpack_AdditiveSchwarz< T >::IsInitialized() const

Returns true if the preconditioner has been successfully initialized.
";

%feature("docstring")  Ifpack_AdditiveSchwarz::IsComputed "virtual
bool Ifpack_AdditiveSchwarz< T >::IsComputed() const

Returns true if the preconditioner has been successfully computed. ";

%feature("docstring")  Ifpack_AdditiveSchwarz::SetParameters "int
Ifpack_AdditiveSchwarz< T >::SetParameters(Teuchos::ParameterList
&List)

Sets the parameters.

Sets the parameter for the additive Schwarz preconditioner, as well as
for all the preconditioners that may need to be defined on each
subblock. Parameters accepted by List are:  \"schwarz: combine mode\"
: It must be an Epetra_CombineMode. Default: Zero. It Can be assume of
the following values: Add: Components on the receiving processor will
be added together;

Zero: Off-processor components will be ignored;

Insert: Off-processor components will be inserted into locations on
receiving processor replacing existing values.

Average: Off-processor components will be averaged with existing;

AbsMax: Magnitudes of Off-processor components will be maxed with
magnitudes of existing components on the receiving processor.

\"schwarz: compute condest\" : if true,  Compute() will estimate the
condition number of the preconditioner. Default: true. ";


// File: structIfpack__AIJMatrix.xml
%feature("docstring") Ifpack_AIJMatrix "";


// File: classIfpack__AMDReordering.xml
%feature("docstring") Ifpack_AMDReordering "

Ifpack_AMDReordering: approximate minimum degree reordering.

C++ includes: Ifpack_AMDReordering.h ";

%feature("docstring")  Ifpack_AMDReordering::Ifpack_AMDReordering "Ifpack_AMDReordering::Ifpack_AMDReordering()

Constructor for Ifpack_Graph's. ";

%feature("docstring")  Ifpack_AMDReordering::Ifpack_AMDReordering "Ifpack_AMDReordering::Ifpack_AMDReordering(const Ifpack_AMDReordering
&RHS)

Copy Constructor. ";

%feature("docstring")  Ifpack_AMDReordering::~Ifpack_AMDReordering "virtual Ifpack_AMDReordering::~Ifpack_AMDReordering()

Destructor. ";

%feature("docstring")  Ifpack_AMDReordering::SetParameter "int
Ifpack_AMDReordering::SetParameter(const string Name, const int Value)

Sets integer parameters `Name'. ";

%feature("docstring")  Ifpack_AMDReordering::SetParameter "int
Ifpack_AMDReordering::SetParameter(const string Name, const double
Value)

Sets double parameters `Name'. ";

%feature("docstring")  Ifpack_AMDReordering::SetParameters "int
Ifpack_AMDReordering::SetParameters(Teuchos::ParameterList &List)

Sets all parameters. ";

%feature("docstring")  Ifpack_AMDReordering::Compute "int
Ifpack_AMDReordering::Compute(const Ifpack_Graph &Graph)

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_AMDReordering::Compute "int
Ifpack_AMDReordering::Compute(const Epetra_RowMatrix &Matrix)

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_AMDReordering::IsComputed "bool
Ifpack_AMDReordering::IsComputed() const

Returns true is the reordering object has been successfully
initialized, false otherwise. ";

%feature("docstring")  Ifpack_AMDReordering::Reorder "int
Ifpack_AMDReordering::Reorder(const int i) const

Returns the reordered index of row i. ";

%feature("docstring")  Ifpack_AMDReordering::InvReorder "int
Ifpack_AMDReordering::InvReorder(const int i) const

Returns the inverse reordered index of row i. ";

%feature("docstring")  Ifpack_AMDReordering::P "int
Ifpack_AMDReordering::P(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &Xreord) const

Applies reordering to multivector X, whose local length equals the
number of local rows. ";

%feature("docstring")  Ifpack_AMDReordering::Pinv "int
Ifpack_AMDReordering::Pinv(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &Xinvreord) const

Applies inverse reordering to multivector X, whose local length equals
the number of local rows. ";

%feature("docstring")  Ifpack_AMDReordering::Print "ostream &
Ifpack_AMDReordering::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_AMDReordering::NumMyRows "int
Ifpack_AMDReordering::NumMyRows() const

Returns the number of local rows. ";


// File: classIfpack__Amesos.xml
%feature("docstring") Ifpack_Amesos "

Ifpack_Amesos: a class to use Amesos' factorizations as
preconditioners.

Class Ifpack_Amesos enables the use of Amesos' factorizations as
Ifpack_Preconditioners.

Ifpack_Amesos is just a bare-bone wrap to Amesos. Currently, the only
parameter required recognized by SetParameters() is \"amesos: solver
type\" (defaulted to \"Amesos_Klu\"), which defined the Amesos solver.
The Teuchos list in input to SetParameters() is copied, then the
copied list is used to set the parameters of the Amesos object.

This class works with matrices whose communicator contains only one
process, that is, either serial matrices, or Ifpack_LocalFilter'd
matrices.

WARNING:  The number of flops is NOT updated.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_Amesos.h ";

%feature("docstring")  Ifpack_Amesos::Ifpack_Amesos "Ifpack_Amesos::Ifpack_Amesos(Epetra_RowMatrix *Matrix)

Constructor. ";

%feature("docstring")  Ifpack_Amesos::Ifpack_Amesos "Ifpack_Amesos::Ifpack_Amesos(const Ifpack_Amesos &rhs)

Copy constructor. ";

%feature("docstring")  Ifpack_Amesos::~Ifpack_Amesos "virtual
Ifpack_Amesos::~Ifpack_Amesos()

Destructor. ";

%feature("docstring")  Ifpack_Amesos::SetUseTranspose "int
Ifpack_Amesos::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied (not
implemented).

This flag allows the transpose of the given operator to be used
implicitly.

Parameters:
-----------

UseTranspose_in:  - (In) If true, multiply by the transpose of
operator, otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Ifpack_Amesos::Apply "int
Ifpack_Amesos::Apply(const Epetra_MultiVector &X, Epetra_MultiVector
&Y) const

Applies the matrix to an Epetra_MultiVector.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  - (Out) A Epetra_MultiVector of dimension NumVectors containing
the result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_Amesos::ApplyInverse "int
Ifpack_Amesos::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the preconditioner to X, returns the result in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to be
preconditioned.

Y:  - (Out) A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

WARNING:  In order to work with AztecOO, any implementation of this
method must support the case where X and Y are the same object. ";

%feature("docstring")  Ifpack_Amesos::NormInf "double
Ifpack_Amesos::NormInf() const

Returns the infinity norm of the global matrix (not implemented) ";

%feature("docstring")  Ifpack_Amesos::Label "const char *
Ifpack_Amesos::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_Amesos::UseTranspose "bool
Ifpack_Amesos::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_Amesos::HasNormInf "bool
Ifpack_Amesos::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_Amesos::Comm "const Epetra_Comm &
Ifpack_Amesos::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_Amesos::OperatorDomainMap "const
Epetra_Map & Ifpack_Amesos::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_Amesos::OperatorRangeMap "const
Epetra_Map & Ifpack_Amesos::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_Amesos::IsInitialized "virtual bool
Ifpack_Amesos::IsInitialized() const

Returns true is the preconditioner has been successfully initialized.
";

%feature("docstring")  Ifpack_Amesos::Initialize "int
Ifpack_Amesos::Initialize()

Initializes the preconditioners.

0 if successful, 1 if problems occurred. ";

%feature("docstring")  Ifpack_Amesos::IsComputed "virtual bool
Ifpack_Amesos::IsComputed() const

Returns true if the preconditioner has been successfully computed. ";

%feature("docstring")  Ifpack_Amesos::Compute "int
Ifpack_Amesos::Compute()

Computes the preconditioners.

0 if successful, 1 if problems occurred. ";

%feature("docstring")  Ifpack_Amesos::SetParameters "int
Ifpack_Amesos::SetParameters(Teuchos::ParameterList &List)

Sets all the parameters for the preconditioner.

Parameters currently supported:  \"amesos: solver type\" : Specifies
the solver type for Amesos. Default: Amesos_Klu.

The input list will be copied, then passed to the Amesos object
through Amesos::SetParameters(). ";

%feature("docstring")  Ifpack_Amesos::Matrix "virtual const
Epetra_RowMatrix& Ifpack_Amesos::Matrix() const

Returns a const reference to the internally stored matrix. ";

%feature("docstring")  Ifpack_Amesos::Condest "double
Ifpack_Amesos::Condest(const Ifpack_CondestType CT=Ifpack_Cheap, const
int MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix
*Matrix_in=0)

Returns the estimated condition number, computes it if necessary. ";

%feature("docstring")  Ifpack_Amesos::Condest "virtual double
Ifpack_Amesos::Condest() const

Returns the estimated condition number, never computes it. ";

%feature("docstring")  Ifpack_Amesos::NumInitialize "virtual int
Ifpack_Amesos::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_Amesos::NumCompute "virtual int
Ifpack_Amesos::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_Amesos::NumApplyInverse "virtual int
Ifpack_Amesos::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_Amesos::InitializeTime "virtual double
Ifpack_Amesos::InitializeTime() const

Returns the total time spent in Initialize(). ";

%feature("docstring")  Ifpack_Amesos::ComputeTime "virtual double
Ifpack_Amesos::ComputeTime() const

Returns the total time spent in Compute(). ";

%feature("docstring")  Ifpack_Amesos::ApplyInverseTime "virtual
double Ifpack_Amesos::ApplyInverseTime() const

Returns the total time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_Amesos::InitializeFlops "virtual double
Ifpack_Amesos::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_Amesos::ComputeFlops "virtual double
Ifpack_Amesos::ComputeFlops() const

Returns the total number of flops to computate the preconditioner. ";

%feature("docstring")  Ifpack_Amesos::ApplyInverseFlops "virtual
double Ifpack_Amesos::ApplyInverseFlops() const

Returns the total number of flops to apply the preconditioner. ";

%feature("docstring")  Ifpack_Amesos::List "virtual const
Teuchos::ParameterList& Ifpack_Amesos::List() const ";

%feature("docstring")  Ifpack_Amesos::Print "std::ostream &
Ifpack_Amesos::Print(std::ostream &os) const

Prints on ostream basic information about this object. ";


// File: classIfpack__BlockRelaxation.xml
%feature("docstring") Ifpack_BlockRelaxation "

Ifpack_BlockRelaxation: a class to define block relaxation
preconditioners of Epetra_RowMatrix's.

The Ifpack_BlockRelaxation class enables the construction of block
relaxation preconditioners of an Epetra_RowMatrix.
Ifpack_PointRelaxation is derived from the Ifpack_Preconditioner
class, which is derived from Epetra_Operator. Therefore this object
can be used as preconditioner everywhere an ApplyInverse() method is
required in the preconditioning step.

The class currently support: block Jacobi;

block Gauss-Seidel;

symmetric block Gauss-Seidel.

The idea of block relaxation method is to extend their point
relaxation counterpart (implemented in Ifpack_PointRelaxation), by
working on a group of equation simulteneously. Generally, larger
blocks result in better convergence and increased robusteness.

The user can decide: the number of blocks (say, NumBlocks). If
NumBlocks is equal to the number of rows, then the resulting scheme is
equivalent to a point relaxation scheme;

how to apply the inverse of each diagonal block, by choosing a dense
container or a sparse container. The implementation of block
relaxation schemes requires the application of the inverse of each
diagonal block. This can be done using LAPACK (dense container), or
any Ifpack_Preconditioner derived class (sparse container);

blocks can be defined using a linear decomposition, by a simple greedy
algorithm, or by resorting to METIS.

The following is an example of usage of this preconditioner with dense
containers. First, we include the header files:

Then, we declare the preconditioner. Note that this is done through
the class Ifpack_AdditiveSchwarz (see note below in this section).

The complete list of supported parameters is reported in page
ifp_params. For a presentation of basic relaxation schemes, please
refer to page Ifpack_PointRelaxation.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_BlockRelaxation.h ";

%feature("docstring")  Ifpack_BlockRelaxation::Ifpack_BlockRelaxation
"Ifpack_BlockRelaxation< T >::Ifpack_BlockRelaxation(const
Epetra_RowMatrix *Matrix)

Ifpack_BlockRelaxation constructor with given Epetra_RowMatrix.

Creates an Ifpack_Preconditioner preconditioner.

Parameters:
-----------

In:  Matrix - Pointer to matrix to be preconditioned. ";

%feature("docstring")  Ifpack_BlockRelaxation::~Ifpack_BlockRelaxation
"Ifpack_BlockRelaxation< T >::~Ifpack_BlockRelaxation() ";

%feature("docstring")  Ifpack_BlockRelaxation::Apply "int
Ifpack_BlockRelaxation< T >::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the matrix to an Epetra_MultiVector.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing the
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_BlockRelaxation::ApplyInverse "int
Ifpack_BlockRelaxation< T >::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the block Jacobi preconditioner to X, returns the result in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to be
preconditioned.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_BlockRelaxation::NormInf "virtual
double Ifpack_BlockRelaxation< T >::NormInf() const

Returns the infinity norm of the global matrix (not implemented) ";

%feature("docstring")  Ifpack_BlockRelaxation::SetUseTranspose "virtual int Ifpack_BlockRelaxation< T >::SetUseTranspose(bool
UseTranspose_in) ";

%feature("docstring")  Ifpack_BlockRelaxation::Label "const char *
Ifpack_BlockRelaxation< T >::Label() const ";

%feature("docstring")  Ifpack_BlockRelaxation::UseTranspose "virtual
bool Ifpack_BlockRelaxation< T >::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_BlockRelaxation::HasNormInf "virtual
bool Ifpack_BlockRelaxation< T >::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_BlockRelaxation::Comm "const
Epetra_Comm & Ifpack_BlockRelaxation< T >::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_BlockRelaxation::OperatorDomainMap "const Epetra_Map & Ifpack_BlockRelaxation< T >::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_BlockRelaxation::OperatorRangeMap "const Epetra_Map & Ifpack_BlockRelaxation< T >::OperatorRangeMap()
const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_BlockRelaxation::NumLocalBlocks "int
Ifpack_BlockRelaxation< T >::NumLocalBlocks() const

Returns the number local blocks. ";

%feature("docstring")  Ifpack_BlockRelaxation::IsInitialized "virtual
bool Ifpack_BlockRelaxation< T >::IsInitialized() const

Returns true if the preconditioner has been successfully computed. ";

%feature("docstring")  Ifpack_BlockRelaxation::IsComputed "virtual
bool Ifpack_BlockRelaxation< T >::IsComputed() const

Returns true if the preconditioner has been successfully computed. ";

%feature("docstring")  Ifpack_BlockRelaxation::SetParameters "int
Ifpack_BlockRelaxation< T >::SetParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the preconditioner. ";

%feature("docstring")  Ifpack_BlockRelaxation::Initialize "int
Ifpack_BlockRelaxation< T >::Initialize()

Initializes the preconditioner. ";

%feature("docstring")  Ifpack_BlockRelaxation::Compute "int
Ifpack_BlockRelaxation< T >::Compute()

Computes the preconditioner. ";

%feature("docstring")  Ifpack_BlockRelaxation::Matrix "virtual const
Epetra_RowMatrix& Ifpack_BlockRelaxation< T >::Matrix() const

Returns a pointer to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_BlockRelaxation::Condest "virtual
double Ifpack_BlockRelaxation< T >::Condest(const Ifpack_CondestType
CT=Ifpack_Cheap, const int MaxIters=1550, const double Tol=1e-9,
Epetra_RowMatrix *Matrix_in=0)

Computes the condition number estimate, returns its value. ";

%feature("docstring")  Ifpack_BlockRelaxation::Condest "virtual
double Ifpack_BlockRelaxation< T >::Condest() const

Returns the computed condition number estimate, or -1.0 if not
computed. ";

%feature("docstring")  Ifpack_BlockRelaxation::Print "std::ostream&
Ifpack_BlockRelaxation< T >::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_BlockRelaxation::NumInitialize "virtual
int Ifpack_BlockRelaxation< T >::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_BlockRelaxation::NumCompute "virtual
int Ifpack_BlockRelaxation< T >::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_BlockRelaxation::NumApplyInverse "virtual int Ifpack_BlockRelaxation< T >::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_BlockRelaxation::InitializeTime "virtual double Ifpack_BlockRelaxation< T >::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_BlockRelaxation::ComputeTime "virtual
double Ifpack_BlockRelaxation< T >::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_BlockRelaxation::ApplyInverseTime "virtual double Ifpack_BlockRelaxation< T >::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_BlockRelaxation::InitializeFlops "virtual double Ifpack_BlockRelaxation< T >::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_BlockRelaxation::ComputeFlops "virtual
double Ifpack_BlockRelaxation< T >::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_BlockRelaxation::ApplyInverseFlops "virtual double Ifpack_BlockRelaxation< T >::ApplyInverseFlops() const

Returns the number of flops in the application of the preconditioner.
";


// File: classIfpack__Chebyshev.xml
%feature("docstring") Ifpack_Chebyshev "

Ifpack_Chebyshev: class for preconditioning with Chebyshev polynomials
in Ifpack.

The Ifpack_Chebyshev class enables the construction of preconditioners
based on Chebyshev polynomials for an Epetra_RowMatrix.
Ifpack_Chebyshev is derived from the Ifpack_Preconditioner class,
which is itself derived from Epetra_Operator. Therefore this object
can be used as preconditioner everywhere an ApplyInverse() method is
required in the preconditioning step.

The class is an adaptation of the routine ML_Cheby in
Smoother/ml_smoother.h

(04/04/06) Flops are not counted in the routine ApplyInverse()

(04/04/06) The switch to use the transpose matrix is not used in
ApplyInverse()

The list of parameters is EigRatio_ = List.get(\"chebyshev: ratio
eigenvalue\", EigRatio_); this is the ratio to define the lower bound
on the spectrum; lambda^* = LambdaMax_ / EigRatio_; a typical value
used in ML is 30.0 (30.0 is the default value).

LambdaMin_ = List.get(\"chebyshev: min eigenvalue\", LambdaMin_); this
is the smallest eigenvalue; this parameter is optional and is only
accessed to check whether the input matrix is equal to identity.

LambdaMax_ = List.get(\"chebyshev: max eigenvalue\", LambdaMax_); this
is the largest eigenvalue of the matrix.

PolyDegree_ = List.get(\"chebyshev: degree\",PolyDegree_); this is the
polynomial degree.

MinDiagonalValue_ = List.get(\"chebyshev: min diagonal value\",
MinDiagonalValue_); this defines the threshold for diagonal values
under which they are not inverted

ZeroStartingSolution_ = List.get(\"chebyshev: zero starting
solution\", ZeroStartingSolution_); this flag allows to set a non-zero
initial guess.

Ulrich Hetmaniuk. SNL 1414.

C++ includes: Ifpack_Chebyshev.h ";

%feature("docstring")  Ifpack_Chebyshev::Ifpack_Chebyshev "Ifpack_Chebyshev::Ifpack_Chebyshev(const Epetra_Operator *Matrix)

Ifpack_Chebyshev constructor with given
Epetra_Operator/Epetra_RowMatrix.

Creates an instance of Ifpack_Chebyshev class.

Parameters:
-----------

Matrix:  - (In) Pointer to the operator to precondition. ";

%feature("docstring")  Ifpack_Chebyshev::Ifpack_Chebyshev "Ifpack_Chebyshev::Ifpack_Chebyshev(const Epetra_RowMatrix *Matrix)

Ifpack_Chebyshev constructor with given
Epetra_Operator/Epetra_RowMatrix.

Creates an instance of Ifpack_Chebyshev class.

Parameters:
-----------

Matrix:  - (In) Pointer to the matrix to precondition. ";

%feature("docstring")  Ifpack_Chebyshev::~Ifpack_Chebyshev "virtual
Ifpack_Chebyshev::~Ifpack_Chebyshev()

Destructor. ";

%feature("docstring")  Ifpack_Chebyshev::Apply "int
Ifpack_Chebyshev::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the matrix to an Epetra_MultiVector.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  - (Out) A Epetra_MultiVector of dimension NumVectors containing
the result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_Chebyshev::ApplyInverse "int
Ifpack_Chebyshev::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the preconditioner to X, returns the result in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to be
preconditioned.

Y:  - (InOut) A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

WARNING:  This routine is NOT AztecOO complaint. ";

%feature("docstring")  Ifpack_Chebyshev::NormInf "virtual double
Ifpack_Chebyshev::NormInf() const

Returns the infinity norm of the global matrix (not implemented) ";

%feature("docstring")  Ifpack_Chebyshev::Label "virtual const char*
Ifpack_Chebyshev::Label() const ";

%feature("docstring")  Ifpack_Chebyshev::UseTranspose "virtual bool
Ifpack_Chebyshev::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_Chebyshev::HasNormInf "virtual bool
Ifpack_Chebyshev::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_Chebyshev::Comm "const Epetra_Comm &
Ifpack_Chebyshev::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_Chebyshev::OperatorDomainMap "const
Epetra_Map & Ifpack_Chebyshev::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_Chebyshev::OperatorRangeMap "const
Epetra_Map & Ifpack_Chebyshev::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_Chebyshev::Initialize "int
Ifpack_Chebyshev::Initialize()

Computes all it is necessary to initialize the preconditioner. ";

%feature("docstring")  Ifpack_Chebyshev::IsInitialized "virtual bool
Ifpack_Chebyshev::IsInitialized() const

Returns true if the preconditioner has been successfully initialized,
false otherwise. ";

%feature("docstring")  Ifpack_Chebyshev::IsComputed "virtual bool
Ifpack_Chebyshev::IsComputed() const

Returns true if the preconditioner has been successfully computed. ";

%feature("docstring")  Ifpack_Chebyshev::Compute "int
Ifpack_Chebyshev::Compute()

Computes the preconditioners. ";

%feature("docstring")  Ifpack_Chebyshev::GetLambdaMax "virtual double
Ifpack_Chebyshev::GetLambdaMax()

Returns an approximation to the largest eigenvalue. ";

%feature("docstring")  Ifpack_Chebyshev::GetLambdaMin "virtual double
Ifpack_Chebyshev::GetLambdaMin()

Contains an approximation to the smallest eigenvalue. ";

%feature("docstring")  Ifpack_Chebyshev::Matrix "virtual const
Epetra_RowMatrix& Ifpack_Chebyshev::Matrix() const

Returns a pointer to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_Chebyshev::Condest "double
Ifpack_Chebyshev::Condest(const Ifpack_CondestType CT=Ifpack_Cheap,
const int MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix
*Matrix_in=0)

Computes the condition number estimates and returns the value. ";

%feature("docstring")  Ifpack_Chebyshev::Condest "virtual double
Ifpack_Chebyshev::Condest() const

Returns the condition number estimate, or -1.0 if not computed. ";

%feature("docstring")  Ifpack_Chebyshev::SetParameters "int
Ifpack_Chebyshev::SetParameters(Teuchos::ParameterList &List)

Sets all the parameters for the preconditioner. ";

%feature("docstring")  Ifpack_Chebyshev::Print "ostream &
Ifpack_Chebyshev::Print(ostream &os) const

Prints object to an output stream. ";

%feature("docstring")  Ifpack_Chebyshev::NumInitialize "virtual int
Ifpack_Chebyshev::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_Chebyshev::NumCompute "virtual int
Ifpack_Chebyshev::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_Chebyshev::NumApplyInverse "virtual int
Ifpack_Chebyshev::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_Chebyshev::InitializeTime "virtual
double Ifpack_Chebyshev::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_Chebyshev::ComputeTime "virtual double
Ifpack_Chebyshev::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_Chebyshev::ApplyInverseTime "virtual
double Ifpack_Chebyshev::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_Chebyshev::InitializeFlops "virtual
double Ifpack_Chebyshev::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_Chebyshev::ComputeFlops "virtual double
Ifpack_Chebyshev::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_Chebyshev::ApplyInverseFlops "virtual
double Ifpack_Chebyshev::ApplyInverseFlops() const

Returns the number of flops for the application of the preconditioner.
";

%feature("docstring")  Ifpack_Chebyshev::PowerMethod "int
Ifpack_Chebyshev::PowerMethod(const Epetra_Operator &Operator, const
Epetra_Vector &InvPointDiagonal, const int MaximumIterations, double
&LambdaMax)

Simple power method to compute lambda_max. ";

%feature("docstring")  Ifpack_Chebyshev::CG "int
Ifpack_Chebyshev::CG(const Epetra_Operator &Operator, const
Epetra_Vector &InvPointDiagonal, const int MaximumIterations, double
&lambda_min, double &lambda_max)

Uses AztecOO's CG to estimate lambda_min and lambda_max. ";

%feature("docstring")  Ifpack_Chebyshev::SetUseTranspose "virtual int
Ifpack_Chebyshev::SetUseTranspose(bool UseTranspose_in)

This flag can be used to apply the preconditioner to the transpose of
the input operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";


// File: classIfpack__Container.xml
%feature("docstring") Ifpack_Container "

Ifpack_Container: a pure virtual class for creating and solving local
linear problems.

Class Ifpack_Container provides the abstract interfaces for
containers. A \"container\" is an object that hosts all it is
necessary to create, populate, and solve local linear problems. The
local linear problem matrix, B, is a submatrix of the local components
of a distributed matrix, A. The idea of container is to specify the
rows of A that are contained in B, then extract B from A, and compute
all it is necessary to solve a linear system in B. Then, set starting
solution (if necessary) and right-hand side for B, and solve the
linear system in B.

A container should be used in the following manner: Create an
container object, specifying the number of rows of B.

If necessary, set parameters for the solution using SetParameters().

Initialize the container by calling Initialize().

Specify the ID of the local rows of A that are contained in B, using
ID().

Prepare the linear system solver using Compute().

set LHS and/or RHS elements using LHS() and RHS().

Solve the linear system using ApplyInverse().

Get the componenets of the computed solution using LHS().

The number of vectors can be set using SetNumVectors(), and it is
defaulted to 1.

Containers are currently used by class Ifpack_BlockRelaxation.

Ifpack_Container is a pure virtual class. Two concrete implementations
are provided in classes Ifpack_SparseContainer (that stores matrices
in sparse the format Epetra_CrsMatrix) and Ifpack_DenseContainer (for
relatively small matrices, as matrices are stored as
Epetra_SerialDenseMatrix's).

Still to do: Flops count has to be tested.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_Container.h ";

%feature("docstring")  Ifpack_Container::~Ifpack_Container "virtual
Ifpack_Container::~Ifpack_Container()

Destructor. ";

%feature("docstring")  Ifpack_Container::NumRows "virtual int
Ifpack_Container::NumRows() const =0

Returns the number of rows of the matrix and LHS/RHS. ";

%feature("docstring")  Ifpack_Container::NumVectors "virtual int
Ifpack_Container::NumVectors() const =0

Returns the number of vectors in LHS/RHS. ";

%feature("docstring")  Ifpack_Container::SetNumVectors "virtual int
Ifpack_Container::SetNumVectors(const int i)=0

Sets the number of vectors for LHS/RHS. ";

%feature("docstring")  Ifpack_Container::LHS "virtual double&
Ifpack_Container::LHS(const int i, const int Vector=0)=0

Returns the i-th component of the vector Vector of LHS. ";

%feature("docstring")  Ifpack_Container::RHS "virtual double&
Ifpack_Container::RHS(const int i, const int Vector=0)=0

Returns the i-th component of the vector Vector of RHS. ";

%feature("docstring")  Ifpack_Container::ID "virtual int&
Ifpack_Container::ID(const int i)=0

Returns the ID associated to local row i.

The set of (local) rows assigned to this container is defined by
calling ID(i) = j, where i (from 0 to NumRows()) indicates the
container-row, and j indicates the local row in the calling process.

This is usually used to recorder the local row ID (on calling process)
of the i-th row in the container. ";

%feature("docstring")  Ifpack_Container::SetMatrixElement "virtual
int Ifpack_Container::SetMatrixElement(const int row, const int col,
const double value)=0

Set the matrix element (row,col) to value. ";

%feature("docstring")  Ifpack_Container::Initialize "virtual int
Ifpack_Container::Initialize()=0

Initializes the container, by performing all operations that only
require matrix structure. ";

%feature("docstring")  Ifpack_Container::Compute "virtual int
Ifpack_Container::Compute(const Epetra_RowMatrix &A)=0

Finalizes the linear system matrix and prepares for the application of
the inverse. ";

%feature("docstring")  Ifpack_Container::SetParameters "virtual int
Ifpack_Container::SetParameters(Teuchos::ParameterList &List)=0

Sets all necessary parameters. ";

%feature("docstring")  Ifpack_Container::IsInitialized "virtual bool
Ifpack_Container::IsInitialized() const =0

Returns true is the container has been successfully initialized. ";

%feature("docstring")  Ifpack_Container::IsComputed "virtual bool
Ifpack_Container::IsComputed() const =0

Returns true is the container has been successfully computed. ";

%feature("docstring")  Ifpack_Container::Apply "virtual int
Ifpack_Container::Apply()=0

Apply the matrix to RHS, results are stored in LHS. ";

%feature("docstring")  Ifpack_Container::ApplyInverse "virtual int
Ifpack_Container::ApplyInverse()=0

Apply the inverse of the matrix to RHS, results are stored in LHS. ";

%feature("docstring")  Ifpack_Container::Label "virtual const char*
Ifpack_Container::Label() const =0

Returns the label of this container. ";

%feature("docstring")  Ifpack_Container::InitializeFlops "virtual
double Ifpack_Container::InitializeFlops() const =0

Returns the flops in Initialize(). ";

%feature("docstring")  Ifpack_Container::ComputeFlops "virtual double
Ifpack_Container::ComputeFlops() const =0

Returns the flops in Compute(). ";

%feature("docstring")  Ifpack_Container::ApplyFlops "virtual double
Ifpack_Container::ApplyFlops() const =0

Returns the flops in Apply(). ";

%feature("docstring")  Ifpack_Container::ApplyInverseFlops "virtual
double Ifpack_Container::ApplyInverseFlops() const =0

Returns the flops in ApplyInverse(). ";

%feature("docstring")  Ifpack_Container::Print "virtual ostream&
Ifpack_Container::Print(std::ostream &os) const =0

Prints out basic information about the container. ";


// File: classIfpack__CrsGraph.xml
%feature("docstring") Ifpack_CrsGraph "";

%feature("docstring")  Ifpack_CrsGraph::~Ifpack_CrsGraph "virtual
Ifpack_CrsGraph::~Ifpack_CrsGraph() ";

%feature("docstring")  Ifpack_CrsGraph::NumRows "virtual int
Ifpack_CrsGraph::NumRows() const =0 ";

%feature("docstring")  Ifpack_CrsGraph::NumCols "virtual int
Ifpack_CrsGraph::NumCols() const =0 ";

%feature("docstring")  Ifpack_CrsGraph::IndexBase "virtual int
Ifpack_CrsGraph::IndexBase() const =0 ";

%feature("docstring")  Ifpack_CrsGraph::NumIndices "virtual int
Ifpack_CrsGraph::NumIndices(int Row) const =0 ";

%feature("docstring")  Ifpack_CrsGraph::ExtractRowCopy "virtual int*
Ifpack_CrsGraph::ExtractRowCopy(int Row, int LenOfIndices, int
&NumIndices, int *&Indices) const =0 ";


// File: classIfpack__CrsIct.xml
%feature("docstring") Ifpack_CrsIct "

Ifpack_CrsIct: A class for constructing and using an incomplete
Cholesky factorization of a given Epetra_CrsMatrix.

The Ifpack_CrsIct class computes a threshold based incomplete LDL^T
factorization of a given Epetra_CrsMatrix. The factorization that is
produced is a function of several parameters: Maximum number of
entries per row/column in factor - The factorization will contain at
most this number of nonzero terms in each row/column of the
factorization.

Diagonal perturbation - Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing. If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A. Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only. Below we discuss the
details of diagonal perturbations. The absolute and relative threshold
values are set by calling SetAbsoluteThreshold() and
SetRelativeThreshold(), respectively.

Estimating Preconditioner Condition Numbers

For ill-conditioned matrices, we often have difficulty computing
usable incomplete factorizations. The most common source of problems
is that the factorization may encounter a small or zero pivot, in
which case the factorization can fail, or even if the factorization
succeeds, the factors may be so poorly conditioned that use of them in
the iterative phase produces meaningless results. Before we can fix
this problem, we must be able to detect it. To this end, we use a
simple but effective condition number estimate for $(LU)^{-1}$.

The condition of a matrix $B$, called $cond_p(B)$, is defined as
$cond_p(B) = \\\\|B\\\\|_p\\\\|B^{-1}\\\\|_p$ in some appropriate norm
$p$. $cond_p(B)$ gives some indication of how many accurate floating
point digits can be expected from operations involving the matrix and
its inverse. A condition number approaching the accuracy of a given
floating point number system, about 15 decimal digits in IEEE double
precision, means that any results involving $B$ or $B^{-1}$ may be
meaningless.

The $\\\\infty$-norm of a vector $y$ is defined as the maximum of the
absolute values of the vector entries, and the $\\\\infty$-norm of a
matrix C is defined as $\\\\|C\\\\|_\\\\infty =
\\\\max_{\\\\|y\\\\|_\\\\infty = 1} \\\\|Cy\\\\|_\\\\infty$. A crude
lower bound for the $cond_\\\\infty(C)$ is
$\\\\|C^{-1}e\\\\|_\\\\infty$ where $e = (1, 1, \\\\ldots, 1)^T$. It
is a lower bound because $cond_\\\\infty(C) =
\\\\|C\\\\|_\\\\infty\\\\|C^{-1}\\\\|_\\\\infty \\\\ge
\\\\|C^{-1}\\\\|_\\\\infty \\\\ge |C^{-1}e\\\\|_\\\\infty$.

For our purposes, we want to estimate $cond_\\\\infty(LU)$, where $L$
and $U$ are our incomplete factors. Edmond in his Ph.D. thesis
demonstrates that $\\\\|(LU)^{-1}e\\\\|_\\\\infty$ provides an
effective estimate for $cond_\\\\infty(LU)$. Furthermore, since
finding $z$ such that $LUz = y$ is a basic kernel for applying the
preconditioner, computing this estimate of $cond_\\\\infty(LU)$ is
performed by setting $y = e$, calling the solve kernel to compute $z$
and then computing $\\\\|z\\\\|_\\\\infty$.

A priori Diagonal Perturbations

Given the above method to estimate the conditioning of the incomplete
factors, if we detect that our factorization is too ill-conditioned we
can improve the conditioning by perturbing the matrix diagonal and
restarting the factorization using this more diagonally dominant
matrix. In order to apply perturbation, prior to starting the
factorization, we compute a diagonal perturbation of our matrix $A$
and perform the factorization on this perturbed matrix. The overhead
cost of perturbing the diagonal is minimal since the first step in
computing the incomplete factors is to copy the matrix $A$ into the
memory space for the incomplete factors. We simply compute the
perturbed diagonal at this point.

The actual perturbation values we use are the diagonal values $(d_1,
d_2, \\\\ldots, d_n)$ with $d_i = sgn(d_i)\\\\alpha + d_i\\\\rho$,
$i=1, 2, \\\\ldots, n$, where $n$ is the matrix dimension and
$sgn(d_i)$ returns the sign of the diagonal entry. This has the effect
of forcing the diagonal values to have minimal magnitude of
$\\\\alpha$ and to increase each by an amount proportional to
$\\\\rho$, and still keep the sign of the original diagonal entry.

Constructing Ifpack_CrsIct objects

Constructing Ifpack_CrsIct objects is a multi-step process. The basic
steps are as follows: Create Ifpack_CrsIct instance, including
storage, via constructor.

Enter values via one or more Put or SumInto functions.

Complete construction via FillComplete call.

Note that, even after a matrix is constructed, it is possible to
update existing matrix entries. It is not possible to create new
entries.

Counting Floating Point Operations

Each Ifpack_CrsIct object keep track of the number of serial floating
point operations performed using the specified object as the this
argument to the function. The Flops() function returns this number as
a double precision number. Using this information, in conjunction with
the Epetra_Time class, one can get accurate parallel performance
numbers. The ResetFlops() function resets the floating point counter.

WARNING:  A Epetra_Map is required for the Ifpack_CrsIct constructor.

C++ includes: Ifpack_CrsIct.h ";

%feature("docstring")  Ifpack_CrsIct::Label "const char*
Ifpack_CrsIct::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_CrsIct::SetUseTranspose "int
Ifpack_CrsIct::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose_in -If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_CrsIct::Apply "int
Ifpack_CrsIct::Apply(const Epetra_MultiVector &X, Epetra_MultiVector
&Y) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Note that this implementation of Apply does NOT perform a forward back
solve with the LDU factorization. Instead it applies these operators
via multiplication with U, D and L respectively. The ApplyInverse()
method performs a solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsIct::ApplyInverse "int
Ifpack_CrsIct::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

In this implementation, we use several existing attributes to
determine how virtual method ApplyInverse() should call the concrete
method Solve(). We pass in the UpperTriangular(), the
Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most
notable warning is that if a matrix has no diagonal values we assume
that there is an implicit unit diagonal that should be accounted for
when doing a triangular solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsIct::NormInf "double
Ifpack_CrsIct::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_CrsIct::HasNormInf "bool
Ifpack_CrsIct::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_CrsIct::UseTranspose "bool
Ifpack_CrsIct::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_CrsIct::OperatorDomainMap "const
Epetra_Map& Ifpack_CrsIct::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_CrsIct::OperatorRangeMap "const
Epetra_Map& Ifpack_CrsIct::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_CrsIct::Comm "const Epetra_Comm&
Ifpack_CrsIct::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_CrsIct::Ifpack_CrsIct "Ifpack_CrsIct::Ifpack_CrsIct(const Epetra_CrsMatrix &A, double
Droptol=1.0E-4, int Lfil=20)

Ifpack_CrsIct constuctor with variable number of indices per row.

Creates a Ifpack_CrsIct object and allocates storage.

Parameters:
-----------

In:  A - User matrix to be factored.

In:  Graph - Graph generated by Ifpack_IlukGraph. ";

%feature("docstring")  Ifpack_CrsIct::Ifpack_CrsIct "Ifpack_CrsIct::Ifpack_CrsIct(const Ifpack_CrsIct &IctOperator)

Copy constructor. ";

%feature("docstring")  Ifpack_CrsIct::~Ifpack_CrsIct "Ifpack_CrsIct::~Ifpack_CrsIct()

Ifpack_CrsIct Destructor. ";

%feature("docstring")  Ifpack_CrsIct::SetAbsoluteThreshold "void
Ifpack_CrsIct::SetAbsoluteThreshold(double Athresh)

Set absolute threshold value. ";

%feature("docstring")  Ifpack_CrsIct::SetRelativeThreshold "void
Ifpack_CrsIct::SetRelativeThreshold(double Rthresh)

Set relative threshold value. ";

%feature("docstring")  Ifpack_CrsIct::SetOverlapMode "void
Ifpack_CrsIct::SetOverlapMode(Epetra_CombineMode OverlapMode)

Set overlap mode type. ";

%feature("docstring")  Ifpack_CrsIct::SetParameters "int
Ifpack_CrsIct::SetParameters(const Teuchos::ParameterList
&parameterlist, bool cerr_warning_if_unused=false)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_CrsIct::InitValues "int
Ifpack_CrsIct::InitValues(const Epetra_CrsMatrix &A)

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_CrsIct::ValuesInitialized "bool
Ifpack_CrsIct::ValuesInitialized() const

If values have been initialized, this query returns true, otherwise it
returns false. ";

%feature("docstring")  Ifpack_CrsIct::Factor "int
Ifpack_CrsIct::Factor()

Compute IC factor U using the specified graph, diagonal perturbation
thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_CrsIct::Factored "bool
Ifpack_CrsIct::Factored() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_CrsIct::Solve "int
Ifpack_CrsIct::Solve(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_CrsIct forward/back solve on a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsIct::Multiply "int
Ifpack_CrsIct::Multiply(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of multiplying U, D and U^T in that order on an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, multiply by L^T, D and U^T in that order.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsIct::Condest "int
Ifpack_CrsIct::Condest(bool Trans, double &ConditionNumberEstimate)
const

Returns the maximum over all the condition number estimate for each
local ILU set of factors.

This functions computes a local condition number estimate on each
processor and return the maximum over all processor of the estimate.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

Out:  ConditionNumberEstimate - The maximum across all processors of
the infinity-norm estimate of the condition number of the inverse of
LDU. ";

%feature("docstring")  Ifpack_CrsIct::GetAbsoluteThreshold "double
Ifpack_CrsIct::GetAbsoluteThreshold()

Get absolute threshold value. ";

%feature("docstring")  Ifpack_CrsIct::GetRelativeThreshold "double
Ifpack_CrsIct::GetRelativeThreshold()

Get relative threshold value. ";

%feature("docstring")  Ifpack_CrsIct::GetOverlapMode "Epetra_CombineMode Ifpack_CrsIct::GetOverlapMode()

Get overlap mode type. ";

%feature("docstring")  Ifpack_CrsIct::NumGlobalNonzeros "int
Ifpack_CrsIct::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_CrsIct::NumMyNonzeros "int
Ifpack_CrsIct::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_CrsIct::D "const Epetra_Vector&
Ifpack_CrsIct::D() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_CrsIct::U "const Epetra_CrsMatrix&
Ifpack_CrsIct::U() const

Returns the address of the U factor associated with this factored
matrix. ";


// File: classIfpack__CrsIlut.xml
%feature("docstring") Ifpack_CrsIlut "

Ifpack_CrsIlut: ILUT preconditioner of a given Epetra_RowMatrix.

C++ includes: Ifpack_CrsIlut.h ";

%feature("docstring")  Ifpack_CrsIlut::Ifpack_CrsIlut "Ifpack_CrsIlut::Ifpack_CrsIlut(const Ifpack_OverlapGraph
*OverlapGraph, double DropTol=1.0E-4, double FillTol=1.0)

Constructor using Ifpack_OverlapGraph.

Creates an object from the overlap graph.

Parameters:
-----------

OverlapGraph:  (In) - Graph describing the graph that should be used
for the factors.

DropTol:  (In/Default) - Drop tolerance used by ILUT algorithm.

FillTol:  (In/Default) - Fill tolerance used by ILUT algorithm. ";

%feature("docstring")  Ifpack_CrsIlut::Ifpack_CrsIlut "Ifpack_CrsIlut::Ifpack_CrsIlut(const Epetra_RowMatrix *UserMatrix,
double DropTol=1.0E-4, double FillTol=1.0)

Constructor using Epetra_RowMatrix.

Creates an Ifpack_Graph object from the user graph implicitly defined
by the Epetra_RowMatrix interface.

Parameters:
-----------

RowMatrix:  (In) - An object that has implemented the Epetra_RowMatrix
interface.

DropTol:  (In/Default) - Drop tolerance used by ILUT algorithm.

FillTol:  (In/Default) - Fill tolerance used by ILUT algorithm. ";

%feature("docstring")  Ifpack_CrsIlut::Ifpack_CrsIlut "Ifpack_CrsIlut::Ifpack_CrsIlut(const Ifpack_CrsIlut &Source)

Copy constructor. ";

%feature("docstring")  Ifpack_CrsIlut::~Ifpack_CrsIlut "virtual
Ifpack_CrsIlut::~Ifpack_CrsIlut()

Ifpack_CrsIlut Destructor. ";

%feature("docstring")  Ifpack_CrsIlut::SetDropTol "int
Ifpack_CrsIlut::SetDropTol(double DropTol)

Set Drop tolerance value as defined by the ILUT algorithm. ";

%feature("docstring")  Ifpack_CrsIlut::SetFillTol "int
Ifpack_CrsIlut::SetFillTol(double FillTol)

Set fill tolerance value as defined by the ILUT algorithm. ";

%feature("docstring")  Ifpack_CrsIlut::SetParameters "int
Ifpack_CrsIlut::SetParameters(const Teuchos::ParameterList
&parameterlist, bool cerr_warning_if_unused=false)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_CrsIlut::DropTol "double
Ifpack_CrsIlut::DropTol() const

Set Drop tolerance value as defined by the ILUT algorithm. ";

%feature("docstring")  Ifpack_CrsIlut::FillTol "double
Ifpack_CrsIlut::FillTol() const

Set fill tolerance value as defined by the ILUT algorithm. ";


// File: classIfpack__CrsRick.xml
%feature("docstring") Ifpack_CrsRick "

Ifpack_CrsRick: A class for constructing and using an incomplete
lower/upper (ILU) factorization of a given Epetra_CrsMatrix.

The Ifpack_CrsRick class computes a \"Relaxed\" ILU factorization with
level k fill of a given Epetra_CrsMatrix. The factorization that is
produced is a function of several parameters: The pattern of the
matrix - All fill is derived from the original matrix nonzero
structure. Level zero fill is defined as the original matrix pattern
(nonzero structure), even if the matrix value at an entry is stored as
a zero. (Thus it is possible to add entries to the ILU factors by
adding zero entries the original matrix.)

Level of fill - Starting with the original matrix pattern as level
fill of zero, the next level of fill is determined by analyzing the
graph of the previous level and determining nonzero fill that is a
result of combining entries that were from previous level only (not
the current level). This rule limits fill to entries that are direct
decendents from the previous level graph. Fill for level k is
determined by applying this rule recursively. For sufficiently large
values of k, the fill would eventually be complete and an exact LU
factorization would be computed. Level of fill is defined during the
construction of the Ifpack_IlukGraph object.

Level of overlap - All Ifpack preconditioners work on parallel
distributed memory computers by using the row partitioning the user
input matrix to determine the partitioning for local ILU factors. If
the level of overlap is set to zero, the rows of the user matrix that
are stored on a given processor are treated as a self-contained local
matrix and all column entries that reach to off-processor entries are
ignored. Setting the level of overlap to one tells Ifpack to increase
the size of the local matrix by adding rows that are reached to by
rows owned by this processor. Increasing levels of overlap are defined
recursively in the same way. For sufficiently large levels of overlap,
the entire matrix would be part of each processor's local ILU
factorization process. Level of overlap is defined during the
construction of the Ifpack_IlukGraph object.

Once the factorization is computed, applying the factorization
\\\\(LUy = x\\\\) results in redundant approximations for any elements
of y that correspond to rows that are part of more than one local ILU
factor. The OverlapMode (changed by calling SetOverlapMode()) defines
how these redundancies are handled using the Epetra_CombineMode enum.
The default is to zero out all values of y for rows that were not part
of the original matrix row distribution.

Fraction of relaxation - Ifpack_CrsRick computes the ILU factorization
row-by-row. As entries at a given row are computed, some number of
them will be dropped because they do match the prescribed sparsity
pattern. The relaxation factor determines how these dropped values
will be handled. If the RelaxValue (changed by calling
SetRelaxValue()) is zero, then these extra entries will by dropped.
This is a classical ILU approach. If the RelaxValue is 1, then the sum
of the extra entries will be added to the diagonal. This is a
classical Modified ILU (MILU) approach. If RelaxValue is between 0 and
1, then RelaxValue times the sum of extra entries will be added to the
diagonal.

For most situations, RelaxValue should be set to zero. For certain
kinds of problems, e.g., reservoir modeling, there is a conservation
principle involved such that any operator should obey a zero row-sum
property. MILU was designed for these cases and you should set the
RelaxValue to 1. For other situations, setting RelaxValue to some
nonzero value may improve the stability of factorization, and can be
used if the computed ILU factors are poorly conditioned.

Diagonal perturbation - Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing. If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A. Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only. Below we discuss the
details of diagonal perturbations. The absolute and relative threshold
values are set by calling SetAbsoluteThreshold() and
SetRelativeThreshold(), respectively.

Estimating Preconditioner Condition Numbers

For ill-conditioned matrices, we often have difficulty computing
usable incomplete factorizations. The most common source of problems
is that the factorization may encounter a small or zero pivot, in
which case the factorization can fail, or even if the factorization
succeeds, the factors may be so poorly conditioned that use of them in
the iterative phase produces meaningless results. Before we can fix
this problem, we must be able to detect it. To this end, we use a
simple but effective condition number estimate for $(LU)^{-1}$.

The condition of a matrix $B$, called $cond_p(B)$, is defined as
$cond_p(B) = \\\\|B\\\\|_p\\\\|B^{-1}\\\\|_p$ in some appropriate norm
$p$. $cond_p(B)$ gives some indication of how many accurate floating
point digits can be expected from operations involving the matrix and
its inverse. A condition number approaching the accuracy of a given
floating point number system, about 15 decimal digits in IEEE double
precision, means that any results involving $B$ or $B^{-1}$ may be
meaningless.

The $\\\\infty$-norm of a vector $y$ is defined as the maximum of the
absolute values of the vector entries, and the $\\\\infty$-norm of a
matrix C is defined as $\\\\|C\\\\|_\\\\infty =
\\\\max_{\\\\|y\\\\|_\\\\infty = 1} \\\\|Cy\\\\|_\\\\infty$. A crude
lower bound for the $cond_\\\\infty(C)$ is
$\\\\|C^{-1}e\\\\|_\\\\infty$ where $e = (1, 1, \\\\ldots, 1)^T$. It
is a lower bound because $cond_\\\\infty(C) =
\\\\|C\\\\|_\\\\infty\\\\|C^{-1}\\\\|_\\\\infty \\\\ge
\\\\|C^{-1}\\\\|_\\\\infty \\\\ge |C^{-1}e\\\\|_\\\\infty$.

For our purposes, we want to estimate $cond_\\\\infty(LU)$, where $L$
and $U$ are our incomplete factors. Edmond in his Ph.D. thesis
demonstrates that $\\\\|(LU)^{-1}e\\\\|_\\\\infty$ provides an
effective estimate for $cond_\\\\infty(LU)$. Furthermore, since
finding $z$ such that $LUz = y$ is a basic kernel for applying the
preconditioner, computing this estimate of $cond_\\\\infty(LU)$ is
performed by setting $y = e$, calling the solve kernel to compute $z$
and then computing $\\\\|z\\\\|_\\\\infty$.

A priori Diagonal Perturbations

Given the above method to estimate the conditioning of the incomplete
factors, if we detect that our factorization is too ill-conditioned we
can improve the conditioning by perturbing the matrix diagonal and
restarting the factorization using this more diagonally dominant
matrix. In order to apply perturbation, prior to starting the
factorization, we compute a diagonal perturbation of our matrix $A$
and perform the factorization on this perturbed matrix. The overhead
cost of perturbing the diagonal is minimal since the first step in
computing the incomplete factors is to copy the matrix $A$ into the
memory space for the incomplete factors. We simply compute the
perturbed diagonal at this point.

The actual perturbation values we use are the diagonal values $(d_1,
d_2, \\\\ldots, d_n)$ with $d_i = sgn(d_i)\\\\alpha + d_i\\\\rho$,
$i=1, 2, \\\\ldots, n$, where $n$ is the matrix dimension and
$sgn(d_i)$ returns the sign of the diagonal entry. This has the effect
of forcing the diagonal values to have minimal magnitude of
$\\\\alpha$ and to increase each by an amount proportional to
$\\\\rho$, and still keep the sign of the original diagonal entry.

Constructing Ifpack_CrsRick objects

Constructing Ifpack_CrsRick objects is a multi-step process. The basic
steps are as follows: Create Ifpack_CrsRick instance, including
storage, via constructor.

Enter values via one or more Put or SumInto functions.

Complete construction via FillComplete call.

Note that, even after a matrix is constructed, it is possible to
update existing matrix entries. It is not possible to create new
entries.

Counting Floating Point Operations

Each Ifpack_CrsRick object keep track of the number of serial floating
point operations performed using the specified object as the this
argument to the function. The Flops() function returns this number as
a double precision number. Using this information, in conjunction with
the Epetra_Time class, one can get accurate parallel performance
numbers. The ResetFlops() function resets the floating point counter.

WARNING:  A Epetra_Map is required for the Ifpack_CrsRick constructor.

C++ includes: Ifpack_CrsRick.h ";

%feature("docstring")  Ifpack_CrsRick::Label "char*
Ifpack_CrsRick::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_CrsRick::SetUseTranspose "int
Ifpack_CrsRick::SetUseTranspose(bool UseTranspose)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose -If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_CrsRick::Apply "int
Ifpack_CrsRick::Apply(const Epetra_MultiVector &X, Epetra_MultiVector
&Y) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Note that this implementation of Apply does NOT perform a forward back
solve with the LDU factorization. Instead it applies these operators
via multiplication with U, D and L respectively. The ApplyInverse()
method performs a solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRick::ApplyInverse "int
Ifpack_CrsRick::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

In this implementation, we use several existing attributes to
determine how virtual method ApplyInverse() should call the concrete
method Solve(). We pass in the UpperTriangular(), the
Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most
notable warning is that if a matrix has no diagonal values we assume
that there is an implicit unit diagonal that should be accounted for
when doing a triangular solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRick::NormInf "double
Ifpack_CrsRick::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_CrsRick::HasNormInf "bool
Ifpack_CrsRick::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_CrsRick::UseTranspose "bool
Ifpack_CrsRick::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_CrsRick::OperatorDomainMap "const
Epetra_Map& Ifpack_CrsRick::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_CrsRick::OperatorRangeMap "const
Epetra_Map& Ifpack_CrsRick::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_CrsRick::Ifpack_CrsRick "Ifpack_CrsRick::Ifpack_CrsRick(const Epetra_CrsMatrix &A, const
Ifpack_IlukGraph &Graph)

Ifpack_CrsRick constuctor with variable number of indices per row.

Creates a Ifpack_CrsRick object and allocates storage.

Parameters:
-----------

In:  A - User matrix to be factored.

In:  Graph - Graph generated by Ifpack_IlukGraph. ";

%feature("docstring")  Ifpack_CrsRick::Ifpack_CrsRick "Ifpack_CrsRick::Ifpack_CrsRick(const Ifpack_CrsRick &Matrix)

Copy constructor. ";

%feature("docstring")  Ifpack_CrsRick::~Ifpack_CrsRick "Ifpack_CrsRick::~Ifpack_CrsRick()

Ifpack_CrsRick Destructor. ";

%feature("docstring")  Ifpack_CrsRick::InitValues "int
Ifpack_CrsRick::InitValues()

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U. ";

%feature("docstring")  Ifpack_CrsRick::ValuesInitialized "bool
Ifpack_CrsRick::ValuesInitialized() const

If values have been initialized, this query returns true, otherwise it
returns false. ";

%feature("docstring")  Ifpack_CrsRick::SetRelaxValue "void
Ifpack_CrsRick::SetRelaxValue(double RelaxValue)

Set RILU(k) relaxation parameter. ";

%feature("docstring")  Ifpack_CrsRick::SetAbsoluteThreshold "void
Ifpack_CrsRick::SetAbsoluteThreshold(double Athresh)

Set absolute threshold value. ";

%feature("docstring")  Ifpack_CrsRick::SetRelativeThreshold "void
Ifpack_CrsRick::SetRelativeThreshold(double Rthresh)

Set relative threshold value. ";

%feature("docstring")  Ifpack_CrsRick::SetOverlapMode "void
Ifpack_CrsRick::SetOverlapMode(Epetra_CombineMode OverlapMode)

Set overlap mode type. ";

%feature("docstring")  Ifpack_CrsRick::SetParameters "int
Ifpack_CrsRick::SetParameters(const Teuchos::ParameterList
&parameterlist, bool cerr_warning_if_unused=false)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_CrsRick::Factor "int
Ifpack_CrsRick::Factor()

Compute ILU factors L and U using the specified graph, diagonal
perturbation thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_CrsRick::Factored "bool
Ifpack_CrsRick::Factored() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_CrsRick::Solve "int
Ifpack_CrsRick::Solve(bool Trans, const Epetra_Vector &x,
Epetra_Vector &y) const

Returns the result of a Ifpack_CrsRick forward/back solve on a
Epetra_Vector x in y.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  x -A Epetra_Vector to solve for.

Out:  y -A Epetra_Vector containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRick::Solve "int
Ifpack_CrsRick::Solve(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_CrsRick forward/back solve on a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRick::Multiply "int
Ifpack_CrsRick::Multiply(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of multiplying U, D and U^T in that order on an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, multiply by L^T, D and U^T in that order.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRick::Condest "int
Ifpack_CrsRick::Condest(bool Trans, double &ConditionNumberEstimate)
const

Returns the maximum over all the condition number estimate for each
local ILU set of factors.

This functions computes a local condition number estimate on each
processor and return the maximum over all processor of the estimate.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

Out:  ConditionNumberEstimate - The maximum across all processors of
the infinity-norm estimate of the condition number of the inverse of
LDU. ";

%feature("docstring")  Ifpack_CrsRick::GetRelaxValue "double
Ifpack_CrsRick::GetRelaxValue()

Get RILU(k) relaxation parameter. ";

%feature("docstring")  Ifpack_CrsRick::GetAbsoluteThreshold "double
Ifpack_CrsRick::GetAbsoluteThreshold()

Get absolute threshold value. ";

%feature("docstring")  Ifpack_CrsRick::GetRelativeThreshold "double
Ifpack_CrsRick::GetRelativeThreshold()

Get relative threshold value. ";

%feature("docstring")  Ifpack_CrsRick::GetOverlapMode "Epetra_CombineMode Ifpack_CrsRick::GetOverlapMode()

Get overlap mode type. ";

%feature("docstring")  Ifpack_CrsRick::NumGlobalRows "int
Ifpack_CrsRick::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Ifpack_CrsRick::NumGlobalCols "int
Ifpack_CrsRick::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Ifpack_CrsRick::NumGlobalNonzeros "int
Ifpack_CrsRick::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_CrsRick::NumGlobalDiagonals "virtual
int Ifpack_CrsRick::NumGlobalDiagonals() const

Returns the number of diagonal entries found in the global input
graph. ";

%feature("docstring")  Ifpack_CrsRick::NumMyRows "int
Ifpack_CrsRick::NumMyRows() const

Returns the number of local matrix rows. ";

%feature("docstring")  Ifpack_CrsRick::NumMyCols "int
Ifpack_CrsRick::NumMyCols() const

Returns the number of local matrix columns. ";

%feature("docstring")  Ifpack_CrsRick::NumMyNonzeros "int
Ifpack_CrsRick::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_CrsRick::NumMyDiagonals "virtual int
Ifpack_CrsRick::NumMyDiagonals() const

Returns the number of diagonal entries found in the local input graph.
";

%feature("docstring")  Ifpack_CrsRick::IndexBase "int
Ifpack_CrsRick::IndexBase() const

Returns the index base for row and column indices for this graph. ";

%feature("docstring")  Ifpack_CrsRick::Graph "const Ifpack_IlukGraph&
Ifpack_CrsRick::Graph() const

Returns the address of the Ifpack_IlukGraph associated with this
factored matrix. ";

%feature("docstring")  Ifpack_CrsRick::D "const Epetra_Vector&
Ifpack_CrsRick::D() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_CrsRick::U "const Epetra_CrsMatrix&
Ifpack_CrsRick::U() const

Returns the address of the U factor associated with this factored
matrix. ";


// File: classIfpack__CrsRiluk.xml
%feature("docstring") Ifpack_CrsRiluk "

Ifpack_CrsRiluk: A class for constructing and using an incomplete
lower/upper (ILU) factorization of a given Epetra_RowMatrix.

The Ifpack_CrsRiluk class computes a \"Relaxed\" ILU factorization
with level k fill of a given Epetra_CrsMatrix. The factorization that
is produced is a function of several parameters: The pattern of the
matrix - All fill is derived from the original matrix nonzero
structure. Level zero fill is defined as the original matrix pattern
(nonzero structure), even if the matrix value at an entry is stored as
a zero. (Thus it is possible to add entries to the ILU factors by
adding zero entries the original matrix.)

Level of fill - Starting with the original matrix pattern as level
fill of zero, the next level of fill is determined by analyzing the
graph of the previous level and determining nonzero fill that is a
result of combining entries that were from previous level only (not
the current level). This rule limits fill to entries that are direct
decendents from the previous level graph. Fill for level k is
determined by applying this rule recursively. For sufficiently large
values of k, the fill would eventually be complete and an exact LU
factorization would be computed. Level of fill is defined during the
construction of the Ifpack_IlukGraph object.

Level of overlap - All Ifpack preconditioners work on parallel
distributed memory computers by using the row partitioning the user
input matrix to determine the partitioning for local ILU factors. If
the level of overlap is set to zero, the rows of the user matrix that
are stored on a given processor are treated as a self-contained local
matrix and all column entries that reach to off-processor entries are
ignored. Setting the level of overlap to one tells Ifpack to increase
the size of the local matrix by adding rows that are reached to by
rows owned by this processor. Increasing levels of overlap are defined
recursively in the same way. For sufficiently large levels of overlap,
the entire matrix would be part of each processor's local ILU
factorization process. Level of overlap is defined during the
construction of the Ifpack_IlukGraph object.

Once the factorization is computed, applying the factorization
\\\\(LUy = x\\\\) results in redundant approximations for any elements
of y that correspond to rows that are part of more than one local ILU
factor. The OverlapMode (changed by calling SetOverlapMode()) defines
how these redundancies are handled using the Epetra_CombineMode enum.
The default is to zero out all values of y for rows that were not part
of the original matrix row distribution.

Fraction of relaxation - Ifpack_CrsRiluk computes the ILU
factorization row-by-row. As entries at a given row are computed, some
number of them will be dropped because they do match the prescribed
sparsity pattern. The relaxation factor determines how these dropped
values will be handled. If the RelaxValue (changed by calling
SetRelaxValue()) is zero, then these extra entries will by dropped.
This is a classical ILU approach. If the RelaxValue is 1, then the sum
of the extra entries will be added to the diagonal. This is a
classical Modified ILU (MILU) approach. If RelaxValue is between 0 and
1, then RelaxValue times the sum of extra entries will be added to the
diagonal.

For most situations, RelaxValue should be set to zero. For certain
kinds of problems, e.g., reservoir modeling, there is a conservation
principle involved such that any operator should obey a zero row-sum
property. MILU was designed for these cases and you should set the
RelaxValue to 1. For other situations, setting RelaxValue to some
nonzero value may improve the stability of factorization, and can be
used if the computed ILU factors are poorly conditioned.

Diagonal perturbation - Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing. If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A. Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only. Below we discuss the
details of diagonal perturbations. The absolute and relative threshold
values are set by calling SetAbsoluteThreshold() and
SetRelativeThreshold(), respectively.

Estimating Preconditioner Condition Numbers

For ill-conditioned matrices, we often have difficulty computing
usable incomplete factorizations. The most common source of problems
is that the factorization may encounter a small or zero pivot, in
which case the factorization can fail, or even if the factorization
succeeds, the factors may be so poorly conditioned that use of them in
the iterative phase produces meaningless results. Before we can fix
this problem, we must be able to detect it. To this end, we use a
simple but effective condition number estimate for $(LU)^{-1}$.

The condition of a matrix $B$, called $cond_p(B)$, is defined as
$cond_p(B) = \\\\|B\\\\|_p\\\\|B^{-1}\\\\|_p$ in some appropriate norm
$p$. $cond_p(B)$ gives some indication of how many accurate floating
point digits can be expected from operations involving the matrix and
its inverse. A condition number approaching the accuracy of a given
floating point number system, about 15 decimal digits in IEEE double
precision, means that any results involving $B$ or $B^{-1}$ may be
meaningless.

The $\\\\infty$-norm of a vector $y$ is defined as the maximum of the
absolute values of the vector entries, and the $\\\\infty$-norm of a
matrix C is defined as $\\\\|C\\\\|_\\\\infty =
\\\\max_{\\\\|y\\\\|_\\\\infty = 1} \\\\|Cy\\\\|_\\\\infty$. A crude
lower bound for the $cond_\\\\infty(C)$ is
$\\\\|C^{-1}e\\\\|_\\\\infty$ where $e = (1, 1, \\\\ldots, 1)^T$. It
is a lower bound because $cond_\\\\infty(C) =
\\\\|C\\\\|_\\\\infty\\\\|C^{-1}\\\\|_\\\\infty \\\\ge
\\\\|C^{-1}\\\\|_\\\\infty \\\\ge |C^{-1}e\\\\|_\\\\infty$.

For our purposes, we want to estimate $cond_\\\\infty(LU)$, where $L$
and $U$ are our incomplete factors. Edmond in his Ph.D. thesis
demonstrates that $\\\\|(LU)^{-1}e\\\\|_\\\\infty$ provides an
effective estimate for $cond_\\\\infty(LU)$. Furthermore, since
finding $z$ such that $LUz = y$ is a basic kernel for applying the
preconditioner, computing this estimate of $cond_\\\\infty(LU)$ is
performed by setting $y = e$, calling the solve kernel to compute $z$
and then computing $\\\\|z\\\\|_\\\\infty$.

A priori Diagonal Perturbations

Given the above method to estimate the conditioning of the incomplete
factors, if we detect that our factorization is too ill-conditioned we
can improve the conditioning by perturbing the matrix diagonal and
restarting the factorization using this more diagonally dominant
matrix. In order to apply perturbation, prior to starting the
factorization, we compute a diagonal perturbation of our matrix $A$
and perform the factorization on this perturbed matrix. The overhead
cost of perturbing the diagonal is minimal since the first step in
computing the incomplete factors is to copy the matrix $A$ into the
memory space for the incomplete factors. We simply compute the
perturbed diagonal at this point.

The actual perturbation values we use are the diagonal values $(d_1,
d_2, \\\\ldots, d_n)$ with $d_i = sgn(d_i)\\\\alpha + d_i\\\\rho$,
$i=1, 2, \\\\ldots, n$, where $n$ is the matrix dimension and
$sgn(d_i)$ returns the sign of the diagonal entry. This has the effect
of forcing the diagonal values to have minimal magnitude of
$\\\\alpha$ and to increase each by an amount proportional to
$\\\\rho$, and still keep the sign of the original diagonal entry.

Constructing Ifpack_CrsRiluk objects

Constructing Ifpack_CrsRiluk objects is a multi-step process. The
basic steps are as follows: Create Ifpack_CrsRiluk instance, including
storage, via constructor.

Enter values via one or more Put or SumInto functions.

Complete construction via FillComplete call.

Note that, even after a matrix is constructed, it is possible to
update existing matrix entries. It is not possible to create new
entries.

Counting Floating Point Operations

Each Ifpack_CrsRiluk object keep track of the number of serial
floating point operations performed using the specified object as the
this argument to the function. The Flops() function returns this
number as a double precision number. Using this information, in
conjunction with the Epetra_Time class, one can get accurate parallel
performance numbers. The ResetFlops() function resets the floating
point counter.

WARNING:  A Epetra_Map is required for the Ifpack_CrsRiluk
constructor.

C++ includes: Ifpack_CrsRiluk.h ";

%feature("docstring")  Ifpack_CrsRiluk::Label "const char*
Ifpack_CrsRiluk::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_CrsRiluk::SetUseTranspose "int
Ifpack_CrsRiluk::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose_in -If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_CrsRiluk::Apply "int
Ifpack_CrsRiluk::Apply(const Epetra_MultiVector &X, Epetra_MultiVector
&Y) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Note that this implementation of Apply does NOT perform a forward back
solve with the LDU factorization. Instead it applies these operators
via multiplication with U, D and L respectively. The ApplyInverse()
method performs a solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRiluk::ApplyInverse "int
Ifpack_CrsRiluk::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

In this implementation, we use several existing attributes to
determine how virtual method ApplyInverse() should call the concrete
method Solve(). We pass in the UpperTriangular(), the
Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most
notable warning is that if a matrix has no diagonal values we assume
that there is an implicit unit diagonal that should be accounted for
when doing a triangular solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRiluk::NormInf "double
Ifpack_CrsRiluk::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_CrsRiluk::HasNormInf "bool
Ifpack_CrsRiluk::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_CrsRiluk::UseTranspose "bool
Ifpack_CrsRiluk::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_CrsRiluk::OperatorDomainMap "const
Epetra_Map& Ifpack_CrsRiluk::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_CrsRiluk::OperatorRangeMap "const
Epetra_Map& Ifpack_CrsRiluk::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_CrsRiluk::Comm "const Epetra_Comm&
Ifpack_CrsRiluk::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_CrsRiluk::Ifpack_CrsRiluk "Ifpack_CrsRiluk::Ifpack_CrsRiluk(const Ifpack_IlukGraph &Graph_in)

Ifpack_CrsRiluk constuctor with variable number of indices per row.

Creates a Ifpack_CrsRiluk object and allocates storage.

Parameters:
-----------

In:  Graph_in - Graph generated by Ifpack_IlukGraph. ";

%feature("docstring")  Ifpack_CrsRiluk::Ifpack_CrsRiluk "Ifpack_CrsRiluk::Ifpack_CrsRiluk(const Ifpack_CrsRiluk &Matrix)

Copy constructor. ";

%feature("docstring")  Ifpack_CrsRiluk::~Ifpack_CrsRiluk "Ifpack_CrsRiluk::~Ifpack_CrsRiluk()

Ifpack_CrsRiluk Destructor. ";

%feature("docstring")  Ifpack_CrsRiluk::InitValues "int
Ifpack_CrsRiluk::InitValues(const Epetra_CrsMatrix &A)

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_CrsRiluk::InitValues "int
Ifpack_CrsRiluk::InitValues(const Epetra_VbrMatrix &A)

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_CrsRiluk::ValuesInitialized "bool
Ifpack_CrsRiluk::ValuesInitialized() const

If values have been initialized, this query returns true, otherwise it
returns false. ";

%feature("docstring")  Ifpack_CrsRiluk::SetRelaxValue "void
Ifpack_CrsRiluk::SetRelaxValue(double RelaxValue)

Set RILU(k) relaxation parameter. ";

%feature("docstring")  Ifpack_CrsRiluk::SetAbsoluteThreshold "void
Ifpack_CrsRiluk::SetAbsoluteThreshold(double Athresh)

Set absolute threshold value. ";

%feature("docstring")  Ifpack_CrsRiluk::SetRelativeThreshold "void
Ifpack_CrsRiluk::SetRelativeThreshold(double Rthresh)

Set relative threshold value. ";

%feature("docstring")  Ifpack_CrsRiluk::SetOverlapMode "void
Ifpack_CrsRiluk::SetOverlapMode(Epetra_CombineMode OverlapMode)

Set overlap mode type. ";

%feature("docstring")  Ifpack_CrsRiluk::SetParameters "int
Ifpack_CrsRiluk::SetParameters(const Teuchos::ParameterList
&parameterlist, bool cerr_warning_if_unused=false)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_CrsRiluk::Factor "int
Ifpack_CrsRiluk::Factor()

Compute ILU factors L and U using the specified graph, diagonal
perturbation thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_CrsRiluk::Factored "bool
Ifpack_CrsRiluk::Factored() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_CrsRiluk::Solve "int
Ifpack_CrsRiluk::Solve(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_CrsRiluk forward/back solve on a
Epetra_MultiVector X in Y (works for Epetra_Vectors also).

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRiluk::Multiply "int
Ifpack_CrsRiluk::Multiply(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of multiplying U, D and L in that order on an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, multiply by L^T, D and U^T in that order.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_CrsRiluk::Condest "int
Ifpack_CrsRiluk::Condest(bool Trans, double &ConditionNumberEstimate)
const

Returns the maximum over all the condition number estimate for each
local ILU set of factors.

This functions computes a local condition number estimate on each
processor and return the maximum over all processor of the estimate.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

Out:  ConditionNumberEstimate - The maximum across all processors of
the infinity-norm estimate of the condition number of the inverse of
LDU. ";

%feature("docstring")  Ifpack_CrsRiluk::GetRelaxValue "double
Ifpack_CrsRiluk::GetRelaxValue()

Get RILU(k) relaxation parameter. ";

%feature("docstring")  Ifpack_CrsRiluk::GetAbsoluteThreshold "double
Ifpack_CrsRiluk::GetAbsoluteThreshold()

Get absolute threshold value. ";

%feature("docstring")  Ifpack_CrsRiluk::GetRelativeThreshold "double
Ifpack_CrsRiluk::GetRelativeThreshold()

Get relative threshold value. ";

%feature("docstring")  Ifpack_CrsRiluk::GetOverlapMode "Epetra_CombineMode Ifpack_CrsRiluk::GetOverlapMode()

Get overlap mode type. ";

%feature("docstring")  Ifpack_CrsRiluk::NumGlobalRows "int
Ifpack_CrsRiluk::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Ifpack_CrsRiluk::NumGlobalCols "int
Ifpack_CrsRiluk::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Ifpack_CrsRiluk::NumGlobalNonzeros "int
Ifpack_CrsRiluk::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_CrsRiluk::NumGlobalBlockDiagonals "virtual int Ifpack_CrsRiluk::NumGlobalBlockDiagonals() const

Returns the number of diagonal entries found in the global input
graph. ";

%feature("docstring")  Ifpack_CrsRiluk::NumMyRows "int
Ifpack_CrsRiluk::NumMyRows() const

Returns the number of local matrix rows. ";

%feature("docstring")  Ifpack_CrsRiluk::NumMyCols "int
Ifpack_CrsRiluk::NumMyCols() const

Returns the number of local matrix columns. ";

%feature("docstring")  Ifpack_CrsRiluk::NumMyNonzeros "int
Ifpack_CrsRiluk::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_CrsRiluk::NumMyBlockDiagonals "virtual
int Ifpack_CrsRiluk::NumMyBlockDiagonals() const

Returns the number of diagonal entries found in the local input graph.
";

%feature("docstring")  Ifpack_CrsRiluk::NumMyDiagonals "virtual int
Ifpack_CrsRiluk::NumMyDiagonals() const

Returns the number of nonzero diagonal values found in matrix. ";

%feature("docstring")  Ifpack_CrsRiluk::IndexBase "int
Ifpack_CrsRiluk::IndexBase() const

Returns the index base for row and column indices for this graph. ";

%feature("docstring")  Ifpack_CrsRiluk::Graph "const
Ifpack_IlukGraph& Ifpack_CrsRiluk::Graph() const

Returns the address of the Ifpack_IlukGraph associated with this
factored matrix. ";

%feature("docstring")  Ifpack_CrsRiluk::L "const Epetra_CrsMatrix&
Ifpack_CrsRiluk::L() const

Returns the address of the L factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_CrsRiluk::D "const Epetra_Vector&
Ifpack_CrsRiluk::D() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_CrsRiluk::U "const Epetra_CrsMatrix&
Ifpack_CrsRiluk::U() const

Returns the address of the L factor associated with this factored
matrix. ";


// File: classIfpack__DenseContainer.xml
%feature("docstring") Ifpack_DenseContainer "

Ifpack_DenseContainer: a class to define containers for dense
matrices.

To understand what an IFPACK container is, please refer to the
documentation of the pure virtual class Ifpack_Container. Currently,
containers are used by class Ifpack_BlockRelaxation.

Using block methods, one needs to store all diagonal blocks and to be
also to apply the inverse of each diagonal block. Using class
Ifpack_DenseContainer, one can store the blocks as dense matrices,
which can be advantageous when the blocks are small. Otherwise, class
Ifpack_SparseContainer is probably more appropriate.

A typical use of a container is as follows:

A call to Compute() computes the LU factorization of the linear system
matrix, using LAPACK (more precisely, by calling the corresponding
routines in Epetra_SerialDenseSolver). The default behavior is to
store the matrix factors by overwriting the linear system matrix
itself. This way, method Apply() fails, as the original matrix does no
longer exists. An alternative is to call KeepNonFactoredMatrix(true),
which forces Ifpack_DenseContainer to maintain in memory a copy of the
non-factored matrix.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_DenseContainer.h ";

%feature("docstring")  Ifpack_DenseContainer::Ifpack_DenseContainer "Ifpack_DenseContainer::Ifpack_DenseContainer(const int NumRows_in,
const int NumVectors_in=1)

Default constructor. ";

%feature("docstring")  Ifpack_DenseContainer::Ifpack_DenseContainer "Ifpack_DenseContainer::Ifpack_DenseContainer(const
Ifpack_DenseContainer &rhs)

Copy constructor. ";

%feature("docstring")  Ifpack_DenseContainer::~Ifpack_DenseContainer "virtual Ifpack_DenseContainer::~Ifpack_DenseContainer()

Destructor. ";

%feature("docstring")  Ifpack_DenseContainer::NumRows "int
Ifpack_DenseContainer::NumRows() const

Returns the number of rows of the matrix and LHS/RHS. ";

%feature("docstring")  Ifpack_DenseContainer::NumVectors "virtual int
Ifpack_DenseContainer::NumVectors() const

Returns the number of vectors in LHS/RHS. ";

%feature("docstring")  Ifpack_DenseContainer::SetNumVectors "virtual
int Ifpack_DenseContainer::SetNumVectors(const int NumVectors_in)

Sets the number of vectors for LHS/RHS. ";

%feature("docstring")  Ifpack_DenseContainer::LHS "double &
Ifpack_DenseContainer::LHS(const int i, const int Vector=0)

Returns the i-th component of the vector Vector of LHS. ";

%feature("docstring")  Ifpack_DenseContainer::RHS "double &
Ifpack_DenseContainer::RHS(const int i, const int Vector=0)

Returns the i-th component of the vector Vector of RHS. ";

%feature("docstring")  Ifpack_DenseContainer::ID "int &
Ifpack_DenseContainer::ID(const int i)

Returns the ID associated to local row i.

The set of (local) rows assigned to this container is defined by
calling ID(i) = j, where i (from 0 to NumRows()) indicates the
container-row, and j indicates the local row in the calling process.

This is usually used to recorder the local row ID (on calling process)
of the i-th row in the container. ";

%feature("docstring")  Ifpack_DenseContainer::SetMatrixElement "int
Ifpack_DenseContainer::SetMatrixElement(const int row, const int col,
const double value)

Set the matrix element (row,col) to value. ";

%feature("docstring")  Ifpack_DenseContainer::SetParameters "virtual
int Ifpack_DenseContainer::SetParameters(Teuchos::ParameterList &List)

Sets all necessary parameters. ";

%feature("docstring")  Ifpack_DenseContainer::IsInitialized "virtual
bool Ifpack_DenseContainer::IsInitialized() const

Returns true is the container has been successfully initialized. ";

%feature("docstring")  Ifpack_DenseContainer::IsComputed "virtual
bool Ifpack_DenseContainer::IsComputed() const

Returns true is the container has been successfully computed. ";

%feature("docstring")  Ifpack_DenseContainer::Label "virtual const
char* Ifpack_DenseContainer::Label() const

Returns the label of this container. ";

%feature("docstring")  Ifpack_DenseContainer::SetKeepNonFactoredMatrix
"virtual int Ifpack_DenseContainer::SetKeepNonFactoredMatrix(const
bool flag)

If flag is true, keeps a copy of the non-factored matrix. ";

%feature("docstring")  Ifpack_DenseContainer::KeepNonFactoredMatrix "virtual bool Ifpack_DenseContainer::KeepNonFactoredMatrix() const

Returns KeepNonFactoredMatrix_. ";

%feature("docstring")  Ifpack_DenseContainer::LHS "virtual const
Epetra_SerialDenseMatrix& Ifpack_DenseContainer::LHS() const

Returns the dense vector containing the LHS. ";

%feature("docstring")  Ifpack_DenseContainer::RHS "virtual const
Epetra_SerialDenseMatrix& Ifpack_DenseContainer::RHS() const

Returns the dense vector containing the RHS. ";

%feature("docstring")  Ifpack_DenseContainer::Matrix "virtual const
Epetra_SerialDenseMatrix& Ifpack_DenseContainer::Matrix() const

Returns the dense matrix or its factors. ";

%feature("docstring")  Ifpack_DenseContainer::NonFactoredMatrix "virtual const Epetra_SerialDenseMatrix&
Ifpack_DenseContainer::NonFactoredMatrix() const

Returns the non-factored dense matrix (only if stored). ";

%feature("docstring")  Ifpack_DenseContainer::ID "virtual const
Epetra_IntSerialDenseVector& Ifpack_DenseContainer::ID() const

Returns the integer dense vector of IDs. ";

%feature("docstring")  Ifpack_DenseContainer::Initialize "int
Ifpack_DenseContainer::Initialize()

Initialize the container. ";

%feature("docstring")  Ifpack_DenseContainer::Compute "int
Ifpack_DenseContainer::Compute(const Epetra_RowMatrix &Matrix_in)

Finalizes the linear system matrix and prepares for the application of
the inverse. ";

%feature("docstring")  Ifpack_DenseContainer::Apply "int
Ifpack_DenseContainer::Apply()

Apply the matrix to RHS, results are stored in LHS. ";

%feature("docstring")  Ifpack_DenseContainer::ApplyInverse "int
Ifpack_DenseContainer::ApplyInverse()

Apply the inverse of the matrix to RHS, results are stored in LHS. ";

%feature("docstring")  Ifpack_DenseContainer::InitializeFlops "virtual double Ifpack_DenseContainer::InitializeFlops() const

Returns the flops in Initialize(). ";

%feature("docstring")  Ifpack_DenseContainer::ComputeFlops "virtual
double Ifpack_DenseContainer::ComputeFlops() const

Returns the flops in Compute(). ";

%feature("docstring")  Ifpack_DenseContainer::ApplyFlops "virtual
double Ifpack_DenseContainer::ApplyFlops() const

Returns the flops in Apply(). ";

%feature("docstring")  Ifpack_DenseContainer::ApplyInverseFlops "virtual double Ifpack_DenseContainer::ApplyInverseFlops() const

Returns the flops in ApplyInverse(). ";

%feature("docstring")  Ifpack_DenseContainer::Print "virtual ostream&
Ifpack_DenseContainer::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";


// File: classIfpack__DiagonalFilter.xml
%feature("docstring") Ifpack_DiagonalFilter "

Ifpack_DiagonalFilter: Filter to modify the diagonal entries of a
given Epetra_RowMatrix.

Ifpack_DiagonalFilter modifies the elements on the diagonal.

A typical use is as follows:

Marzio Sala, SNL 9214.  Last modified on 24-Jan-05.

C++ includes: Ifpack_DiagonalFilter.h ";

%feature("docstring")  Ifpack_DiagonalFilter::Ifpack_DiagonalFilter "Ifpack_DiagonalFilter::Ifpack_DiagonalFilter(const
Teuchos::RefCountPtr< Epetra_RowMatrix > &Matrix, double
AbsoluteThreshold, double RelativeThreshold)

Constructor. ";

%feature("docstring")  Ifpack_DiagonalFilter::~Ifpack_DiagonalFilter "virtual Ifpack_DiagonalFilter::~Ifpack_DiagonalFilter()

Destructor. ";

%feature("docstring")  Ifpack_DiagonalFilter::NumMyRowEntries "virtual int Ifpack_DiagonalFilter::NumMyRowEntries(int MyRow, int
&NumEntries) const

Returns the number of entries in MyRow. ";

%feature("docstring")  Ifpack_DiagonalFilter::MaxNumEntries "virtual
int Ifpack_DiagonalFilter::MaxNumEntries() const

Returns the maximum number of entries. ";

%feature("docstring")  Ifpack_DiagonalFilter::ExtractMyRowCopy "int
Ifpack_DiagonalFilter::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  Ifpack_DiagonalFilter::ExtractDiagonalCopy "virtual int Ifpack_DiagonalFilter::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const ";

%feature("docstring")  Ifpack_DiagonalFilter::Multiply "int
Ifpack_DiagonalFilter::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DiagonalFilter::Solve "virtual int
Ifpack_DiagonalFilter::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  Ifpack_DiagonalFilter::Apply "virtual int
Ifpack_DiagonalFilter::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DiagonalFilter::ApplyInverse "virtual
int Ifpack_DiagonalFilter::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DiagonalFilter::InvRowSums "virtual int
Ifpack_DiagonalFilter::InvRowSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_DiagonalFilter::LeftScale "virtual int
Ifpack_DiagonalFilter::LeftScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_DiagonalFilter::InvColSums "virtual int
Ifpack_DiagonalFilter::InvColSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_DiagonalFilter::RightScale "virtual int
Ifpack_DiagonalFilter::RightScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_DiagonalFilter::Filled "virtual bool
Ifpack_DiagonalFilter::Filled() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NormInf "virtual double
Ifpack_DiagonalFilter::NormInf() const

Not implemented for efficiency reasons. ";

%feature("docstring")  Ifpack_DiagonalFilter::NormOne "virtual double
Ifpack_DiagonalFilter::NormOne() const

Not implemented for efficiency reasons. ";

%feature("docstring")  Ifpack_DiagonalFilter::NumGlobalNonzeros "virtual int Ifpack_DiagonalFilter::NumGlobalNonzeros() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumGlobalRows "virtual
int Ifpack_DiagonalFilter::NumGlobalRows() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumGlobalCols "virtual
int Ifpack_DiagonalFilter::NumGlobalCols() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumGlobalDiagonals "virtual int Ifpack_DiagonalFilter::NumGlobalDiagonals() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumMyNonzeros "virtual
int Ifpack_DiagonalFilter::NumMyNonzeros() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumMyRows "virtual int
Ifpack_DiagonalFilter::NumMyRows() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumMyCols "virtual int
Ifpack_DiagonalFilter::NumMyCols() const ";

%feature("docstring")  Ifpack_DiagonalFilter::NumMyDiagonals "virtual
int Ifpack_DiagonalFilter::NumMyDiagonals() const ";

%feature("docstring")  Ifpack_DiagonalFilter::LowerTriangular "virtual bool Ifpack_DiagonalFilter::LowerTriangular() const ";

%feature("docstring")  Ifpack_DiagonalFilter::UpperTriangular "virtual bool Ifpack_DiagonalFilter::UpperTriangular() const ";

%feature("docstring")  Ifpack_DiagonalFilter::RowMatrixRowMap "virtual const Epetra_Map& Ifpack_DiagonalFilter::RowMatrixRowMap()
const ";

%feature("docstring")  Ifpack_DiagonalFilter::RowMatrixColMap "virtual const Epetra_Map& Ifpack_DiagonalFilter::RowMatrixColMap()
const ";

%feature("docstring")  Ifpack_DiagonalFilter::RowMatrixImporter "virtual const Epetra_Import*
Ifpack_DiagonalFilter::RowMatrixImporter() const ";

%feature("docstring")  Ifpack_DiagonalFilter::SetUseTranspose "int
Ifpack_DiagonalFilter::SetUseTranspose(bool UseTranspose_in) ";

%feature("docstring")  Ifpack_DiagonalFilter::UseTranspose "bool
Ifpack_DiagonalFilter::UseTranspose() const ";

%feature("docstring")  Ifpack_DiagonalFilter::HasNormInf "bool
Ifpack_DiagonalFilter::HasNormInf() const

Not implemented for efficiency reasons. ";

%feature("docstring")  Ifpack_DiagonalFilter::Comm "const
Epetra_Comm& Ifpack_DiagonalFilter::Comm() const ";

%feature("docstring")  Ifpack_DiagonalFilter::OperatorDomainMap "const Epetra_Map& Ifpack_DiagonalFilter::OperatorDomainMap() const ";

%feature("docstring")  Ifpack_DiagonalFilter::OperatorRangeMap "const
Epetra_Map& Ifpack_DiagonalFilter::OperatorRangeMap() const ";

%feature("docstring")  Ifpack_DiagonalFilter::Map "const
Epetra_BlockMap& Ifpack_DiagonalFilter::Map() const ";

%feature("docstring")  Ifpack_DiagonalFilter::Label "const char*
Ifpack_DiagonalFilter::Label() const ";


// File: classIfpack__DiagPreconditioner.xml
%feature("docstring") Ifpack_DiagPreconditioner "

Ifpack_DiagPreconditioner: a class for diagonal preconditioning.

C++ includes: Ifpack_DiagPreconditioner.h ";

%feature("docstring")
Ifpack_DiagPreconditioner::Ifpack_DiagPreconditioner "Ifpack_DiagPreconditioner::Ifpack_DiagPreconditioner(const Epetra_Map
&DomainMap, const Epetra_Map &RangeMap, const Epetra_Vector &diag)

ctor ";

%feature("docstring")
Ifpack_DiagPreconditioner::~Ifpack_DiagPreconditioner "Ifpack_DiagPreconditioner::~Ifpack_DiagPreconditioner()

dtor ";

%feature("docstring")  Ifpack_DiagPreconditioner::SetUseTranspose "int Ifpack_DiagPreconditioner::SetUseTranspose(bool UseTranspose_in)
";

%feature("docstring")  Ifpack_DiagPreconditioner::Apply "int
Ifpack_DiagPreconditioner::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DiagPreconditioner::ApplyInverse "int
Ifpack_DiagPreconditioner::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Y.ReciprocalMultiply(1.0, diag_, X, 0.0); ";

%feature("docstring")  Ifpack_DiagPreconditioner::NormInf "double
Ifpack_DiagPreconditioner::NormInf() const ";

%feature("docstring")  Ifpack_DiagPreconditioner::Label "const char*
Ifpack_DiagPreconditioner::Label() const ";

%feature("docstring")  Ifpack_DiagPreconditioner::UseTranspose "bool
Ifpack_DiagPreconditioner::UseTranspose() const ";

%feature("docstring")  Ifpack_DiagPreconditioner::HasNormInf "bool
Ifpack_DiagPreconditioner::HasNormInf() const ";

%feature("docstring")  Ifpack_DiagPreconditioner::Comm "const
Epetra_Comm& Ifpack_DiagPreconditioner::Comm() const ";

%feature("docstring")  Ifpack_DiagPreconditioner::OperatorDomainMap "const Epetra_Map& Ifpack_DiagPreconditioner::OperatorDomainMap() const
";

%feature("docstring")  Ifpack_DiagPreconditioner::OperatorRangeMap "const Epetra_Map& Ifpack_DiagPreconditioner::OperatorRangeMap() const
";

%feature("docstring")  Ifpack_DiagPreconditioner::Map "const
Epetra_BlockMap& Ifpack_DiagPreconditioner::Map() const ";


// File: classIfpack__DropFilter.xml
%feature("docstring") Ifpack_DropFilter "

Ifpack_DropFilter: Filter based on matrix entries.

Ifpack_DropFilter enables the dropping of all elements whose absolute
value is below a specified threshold.

A typical use is as follows:

It is supposed that Ifpack_DropFilter is used on localized matrices.

Marzio Sala, SNL 9214.  Last modified: Oct-04.

C++ includes: Ifpack_DropFilter.h ";

%feature("docstring")  Ifpack_DropFilter::Ifpack_DropFilter "Ifpack_DropFilter::Ifpack_DropFilter(const Teuchos::RefCountPtr<
Epetra_RowMatrix > &Matrix, double DropTol)

Constructor. ";

%feature("docstring")  Ifpack_DropFilter::~Ifpack_DropFilter "virtual
Ifpack_DropFilter::~Ifpack_DropFilter()

Destructor. ";

%feature("docstring")  Ifpack_DropFilter::NumMyRowEntries "virtual
int Ifpack_DropFilter::NumMyRowEntries(int MyRow, int &NumEntries)
const

Returns the number of entries in MyRow. ";

%feature("docstring")  Ifpack_DropFilter::MaxNumEntries "virtual int
Ifpack_DropFilter::MaxNumEntries() const

Returns the maximum number of entries. ";

%feature("docstring")  Ifpack_DropFilter::ExtractMyRowCopy "int
Ifpack_DropFilter::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  Ifpack_DropFilter::ExtractDiagonalCopy "int
Ifpack_DropFilter::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const
";

%feature("docstring")  Ifpack_DropFilter::Multiply "int
Ifpack_DropFilter::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DropFilter::Solve "int
Ifpack_DropFilter::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DropFilter::Apply "int
Ifpack_DropFilter::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DropFilter::ApplyInverse "int
Ifpack_DropFilter::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_DropFilter::InvRowSums "int
Ifpack_DropFilter::InvRowSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_DropFilter::LeftScale "virtual int
Ifpack_DropFilter::LeftScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_DropFilter::InvColSums "int
Ifpack_DropFilter::InvColSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_DropFilter::RightScale "virtual int
Ifpack_DropFilter::RightScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_DropFilter::Filled "virtual bool
Ifpack_DropFilter::Filled() const ";

%feature("docstring")  Ifpack_DropFilter::NormInf "virtual double
Ifpack_DropFilter::NormInf() const ";

%feature("docstring")  Ifpack_DropFilter::NormOne "virtual double
Ifpack_DropFilter::NormOne() const ";

%feature("docstring")  Ifpack_DropFilter::NumGlobalNonzeros "virtual
int Ifpack_DropFilter::NumGlobalNonzeros() const ";

%feature("docstring")  Ifpack_DropFilter::NumGlobalRows "virtual int
Ifpack_DropFilter::NumGlobalRows() const ";

%feature("docstring")  Ifpack_DropFilter::NumGlobalCols "virtual int
Ifpack_DropFilter::NumGlobalCols() const ";

%feature("docstring")  Ifpack_DropFilter::NumGlobalDiagonals "virtual
int Ifpack_DropFilter::NumGlobalDiagonals() const ";

%feature("docstring")  Ifpack_DropFilter::NumMyNonzeros "virtual int
Ifpack_DropFilter::NumMyNonzeros() const ";

%feature("docstring")  Ifpack_DropFilter::NumMyRows "virtual int
Ifpack_DropFilter::NumMyRows() const ";

%feature("docstring")  Ifpack_DropFilter::NumMyCols "virtual int
Ifpack_DropFilter::NumMyCols() const ";

%feature("docstring")  Ifpack_DropFilter::NumMyDiagonals "virtual int
Ifpack_DropFilter::NumMyDiagonals() const ";

%feature("docstring")  Ifpack_DropFilter::LowerTriangular "virtual
bool Ifpack_DropFilter::LowerTriangular() const ";

%feature("docstring")  Ifpack_DropFilter::UpperTriangular "virtual
bool Ifpack_DropFilter::UpperTriangular() const ";

%feature("docstring")  Ifpack_DropFilter::RowMatrixRowMap "virtual
const Epetra_Map& Ifpack_DropFilter::RowMatrixRowMap() const ";

%feature("docstring")  Ifpack_DropFilter::RowMatrixColMap "virtual
const Epetra_Map& Ifpack_DropFilter::RowMatrixColMap() const ";

%feature("docstring")  Ifpack_DropFilter::RowMatrixImporter "virtual
const Epetra_Import* Ifpack_DropFilter::RowMatrixImporter() const ";

%feature("docstring")  Ifpack_DropFilter::SetUseTranspose "int
Ifpack_DropFilter::SetUseTranspose(bool UseTranspose) ";

%feature("docstring")  Ifpack_DropFilter::UseTranspose "bool
Ifpack_DropFilter::UseTranspose() const ";

%feature("docstring")  Ifpack_DropFilter::HasNormInf "bool
Ifpack_DropFilter::HasNormInf() const ";

%feature("docstring")  Ifpack_DropFilter::Comm "const Epetra_Comm&
Ifpack_DropFilter::Comm() const ";

%feature("docstring")  Ifpack_DropFilter::OperatorDomainMap "const
Epetra_Map& Ifpack_DropFilter::OperatorDomainMap() const ";

%feature("docstring")  Ifpack_DropFilter::OperatorRangeMap "const
Epetra_Map& Ifpack_DropFilter::OperatorRangeMap() const ";

%feature("docstring")  Ifpack_DropFilter::Map "const Epetra_BlockMap&
Ifpack_DropFilter::Map() const ";

%feature("docstring")  Ifpack_DropFilter::Label "const char*
Ifpack_DropFilter::Label() const ";


// File: classIfpack__Element.xml
%feature("docstring") Ifpack_Element "";

%feature("docstring")  Ifpack_Element::Ifpack_Element "Ifpack_Element::Ifpack_Element() ";

%feature("docstring")  Ifpack_Element::Ifpack_Element "Ifpack_Element::Ifpack_Element(const Ifpack_Element &rhs) ";

%feature("docstring")  Ifpack_Element::Index "int
Ifpack_Element::Index() const ";

%feature("docstring")  Ifpack_Element::Value "double
Ifpack_Element::Value() const ";

%feature("docstring")  Ifpack_Element::AbsValue "double
Ifpack_Element::AbsValue() const ";

%feature("docstring")  Ifpack_Element::SetIndex "void
Ifpack_Element::SetIndex(const int i) ";

%feature("docstring")  Ifpack_Element::SetValue "void
Ifpack_Element::SetValue(const double val) ";


// File: classIfpack__EquationPartitioner.xml
%feature("docstring") Ifpack_EquationPartitioner "

Ifpack_EquationPartitioner: A class to decompose an Ifpack_Graph so
that each block will contain all the rows for a different equation.

Ifpack_EquationPartitioner enables a decomposition into blocks of
equations. Suppose that the input Ifpack_Graph is based on an
Epetra_RowMatrix, whose rows represent (U_i,V_i,P_i) for each grid
node i. This partitioner will decompose the graph into three
subgraphs, each of them containing the rows of U, then V, than P.

The number of equations is set as the number of local partitions.

It is required that NumRows % NumLocalParts() = 0.

C++ includes: Ifpack_EquationPartitioner.h ";

%feature("docstring")
Ifpack_EquationPartitioner::Ifpack_EquationPartitioner "Ifpack_EquationPartitioner::Ifpack_EquationPartitioner(const
Ifpack_Graph *Graph)

Constructor. ";

%feature("docstring")
Ifpack_EquationPartitioner::~Ifpack_EquationPartitioner "virtual
Ifpack_EquationPartitioner::~Ifpack_EquationPartitioner()

Destructor. ";

%feature("docstring")
Ifpack_EquationPartitioner::SetPartitionParameters "int
Ifpack_EquationPartitioner::SetPartitionParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner. ";

%feature("docstring")  Ifpack_EquationPartitioner::ComputePartitions "int Ifpack_EquationPartitioner::ComputePartitions()

Computes the partitions. Returns 0 if successful. ";


// File: classIfpack__Graph.xml
%feature("docstring") Ifpack_Graph "

Ifpack_Graph: a pure virtual class that defines graphs for IFPACK.

Class Ifpack_Graph defines the abstract interface to use graphs in
IFPACK. This class contains all the functions that are required by
IFPACK classes.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_Graph.h ";

%feature("docstring")  Ifpack_Graph::~Ifpack_Graph "virtual
Ifpack_Graph::~Ifpack_Graph()

Destructor. ";

%feature("docstring")  Ifpack_Graph::NumMyRows "virtual int
Ifpack_Graph::NumMyRows() const =0

Returns the number of local rows. ";

%feature("docstring")  Ifpack_Graph::NumMyCols "virtual int
Ifpack_Graph::NumMyCols() const =0

Returns the number of local columns. ";

%feature("docstring")  Ifpack_Graph::NumGlobalRows "virtual int
Ifpack_Graph::NumGlobalRows() const =0

Returns the number of global rows. ";

%feature("docstring")  Ifpack_Graph::NumGlobalCols "virtual int
Ifpack_Graph::NumGlobalCols() const =0

Returns the number of global columns. ";

%feature("docstring")  Ifpack_Graph::MaxMyNumEntries "virtual int
Ifpack_Graph::MaxMyNumEntries() const =0

Returns the maximun number of entries for row. ";

%feature("docstring")  Ifpack_Graph::NumMyNonzeros "virtual int
Ifpack_Graph::NumMyNonzeros() const =0

Returns the number of local nonzero entries. ";

%feature("docstring")  Ifpack_Graph::Filled "virtual bool
Ifpack_Graph::Filled() const =0

Returns true is graph is filled. ";

%feature("docstring")  Ifpack_Graph::GRID "virtual int
Ifpack_Graph::GRID(int) const =0

Returns the global row ID of input local row. ";

%feature("docstring")  Ifpack_Graph::GCID "virtual int
Ifpack_Graph::GCID(int) const =0

Returns the global column ID of input local column. ";

%feature("docstring")  Ifpack_Graph::LRID "virtual int
Ifpack_Graph::LRID(int) const =0

Returns the local row ID of input global row. ";

%feature("docstring")  Ifpack_Graph::LCID "virtual int
Ifpack_Graph::LCID(int) const =0

Returns the local column ID of input global column. ";

%feature("docstring")  Ifpack_Graph::ExtractMyRowCopy "virtual int
Ifpack_Graph::ExtractMyRowCopy(int MyRow, int LenOfIndices, int
&NumIndices, int *Indices) const =0

Extracts a copy of input local row. ";

%feature("docstring")  Ifpack_Graph::Comm "virtual const Epetra_Comm&
Ifpack_Graph::Comm() const =0

Returns the communicator object of the graph. ";

%feature("docstring")  Ifpack_Graph::Print "virtual ostream&
Ifpack_Graph::Print(std::ostream &os) const =0

Prints basic information about the graph object. ";


// File: classIfpack__Graph__Epetra__CrsGraph.xml
%feature("docstring") Ifpack_Graph_Epetra_CrsGraph "

Ifpack_Graph_Epetra_CrsGraph: a class to define Ifpack_Graph as a
light-weight conversion of Epetra_CrsGraph's.

Class Ifpack_Graph_Epetra_CrsGraph enables the construction of an
Ifpack_Graph based on the input Epetra_CrsGraph. Note that data are
not copied to this object; instead, wrappers are furnished.

C++ includes: Ifpack_Graph_Epetra_CrsGraph.h ";

%feature("docstring")
Ifpack_Graph_Epetra_CrsGraph::Ifpack_Graph_Epetra_CrsGraph "Ifpack_Graph_Epetra_CrsGraph::Ifpack_Graph_Epetra_CrsGraph(const
Teuchos::RefCountPtr< const Epetra_CrsGraph > &CrsGraph)

Constructor. ";

%feature("docstring")
Ifpack_Graph_Epetra_CrsGraph::~Ifpack_Graph_Epetra_CrsGraph "virtual
Ifpack_Graph_Epetra_CrsGraph::~Ifpack_Graph_Epetra_CrsGraph()

Destructor. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::NumMyRows "int
Ifpack_Graph_Epetra_CrsGraph::NumMyRows() const

Returns the number of local rows. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::NumMyCols "int
Ifpack_Graph_Epetra_CrsGraph::NumMyCols() const

Returns the number of local columns. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::NumGlobalRows "int Ifpack_Graph_Epetra_CrsGraph::NumGlobalRows() const

Returns the number of global rows. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::NumGlobalCols "int Ifpack_Graph_Epetra_CrsGraph::NumGlobalCols() const

Returns the number of global columns. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::MaxMyNumEntries "int Ifpack_Graph_Epetra_CrsGraph::MaxMyNumEntries() const

Returns the maximun number of entries for row. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::NumMyNonzeros "int Ifpack_Graph_Epetra_CrsGraph::NumMyNonzeros() const

Returns the number of local nonzero entries. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::Filled "bool
Ifpack_Graph_Epetra_CrsGraph::Filled() const

Returns true is graph is filled. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::GRID "int
Ifpack_Graph_Epetra_CrsGraph::GRID(int) const

Returns the global row ID of input local row. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::GCID "int
Ifpack_Graph_Epetra_CrsGraph::GCID(int) const

Returns the global column ID of input local column. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::LRID "int
Ifpack_Graph_Epetra_CrsGraph::LRID(int) const

Returns the local row ID of input global row. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::LCID "int
Ifpack_Graph_Epetra_CrsGraph::LCID(int) const

Returns the local column ID of input global column. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::ExtractMyRowCopy
"int Ifpack_Graph_Epetra_CrsGraph::ExtractMyRowCopy(int GlobalRow,
int LenOfIndices, int &NumIndices, int *Indices) const

Extracts a copy of input local row. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::Comm "const
Epetra_Comm & Ifpack_Graph_Epetra_CrsGraph::Comm() const

Returns the communicator object of the graph. ";

%feature("docstring")  Ifpack_Graph_Epetra_CrsGraph::Print "ostream &
Ifpack_Graph_Epetra_CrsGraph::Print(std::ostream &os) const

Prints basic information about the graph object. ";


// File: classIfpack__Graph__Epetra__RowMatrix.xml
%feature("docstring") Ifpack_Graph_Epetra_RowMatrix "

Ifpack_Graph_Epetra_RowMatrix: a class to define Ifpack_Graph as a
light-weight conversion of Epetra_RowMatrix's.

Class Ifpack_Graph_Epetra_RowMatrix enables the construction of an
Ifpack_Graph based on the input Epetra_RowMatrix. Note that data are
not copied to this object; instead, wrappers are furnished.

Marzio Sala, SNL 9214

C++ includes: Ifpack_Graph_Epetra_RowMatrix.h ";

%feature("docstring")
Ifpack_Graph_Epetra_RowMatrix::Ifpack_Graph_Epetra_RowMatrix "Ifpack_Graph_Epetra_RowMatrix::Ifpack_Graph_Epetra_RowMatrix(const
Teuchos::RefCountPtr< const Epetra_RowMatrix > &RowMatrix)

Constructor. ";

%feature("docstring")
Ifpack_Graph_Epetra_RowMatrix::~Ifpack_Graph_Epetra_RowMatrix "virtual
Ifpack_Graph_Epetra_RowMatrix::~Ifpack_Graph_Epetra_RowMatrix()

Destructor. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::NumMyRows "int
Ifpack_Graph_Epetra_RowMatrix::NumMyRows() const

Returns the number of local rows. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::NumMyCols "int
Ifpack_Graph_Epetra_RowMatrix::NumMyCols() const

Returns the number of local columns. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::NumGlobalRows "int Ifpack_Graph_Epetra_RowMatrix::NumGlobalRows() const

Returns the number of global rows. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::NumGlobalCols "int Ifpack_Graph_Epetra_RowMatrix::NumGlobalCols() const

Returns the number of global columns. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::MaxMyNumEntries
"int Ifpack_Graph_Epetra_RowMatrix::MaxMyNumEntries() const

Returns the maximun number of entries for row. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::NumMyNonzeros "int Ifpack_Graph_Epetra_RowMatrix::NumMyNonzeros() const

Returns the number of local nonzero entries. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::Filled "bool
Ifpack_Graph_Epetra_RowMatrix::Filled() const

Returns true is graph is filled. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::GRID "int
Ifpack_Graph_Epetra_RowMatrix::GRID(int) const

Returns the global row ID of input local row. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::GCID "int
Ifpack_Graph_Epetra_RowMatrix::GCID(int) const

Returns the global column ID of input local column. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::LRID "int
Ifpack_Graph_Epetra_RowMatrix::LRID(int) const

Returns the local row ID of input global row. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::LCID "int
Ifpack_Graph_Epetra_RowMatrix::LCID(int) const

Returns the local column ID of input global column. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::ExtractMyRowCopy
"int Ifpack_Graph_Epetra_RowMatrix::ExtractMyRowCopy(int GlobalRow,
int LenOfIndices, int &NumIndices, int *Indices) const

Extracts a copy of input local row. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::Comm "const
Epetra_Comm & Ifpack_Graph_Epetra_RowMatrix::Comm() const

Returns the communicator object of the graph. ";

%feature("docstring")  Ifpack_Graph_Epetra_RowMatrix::Print "ostream
& Ifpack_Graph_Epetra_RowMatrix::Print(std::ostream &os) const

Prints basic information abobut the graph object. ";


// File: classIfpack__GreedyPartitioner.xml
%feature("docstring") Ifpack_GreedyPartitioner "

Ifpack_GreedyPartitioner: A class to decompose Ifpack_Graph's using a
simple greedy algorithm.

C++ includes: Ifpack_GreedyPartitioner.h ";

%feature("docstring")
Ifpack_GreedyPartitioner::Ifpack_GreedyPartitioner "Ifpack_GreedyPartitioner::Ifpack_GreedyPartitioner(const Ifpack_Graph
*Graph)

Constructor. ";

%feature("docstring")
Ifpack_GreedyPartitioner::~Ifpack_GreedyPartitioner "virtual
Ifpack_GreedyPartitioner::~Ifpack_GreedyPartitioner()

Destructor. ";

%feature("docstring")
Ifpack_GreedyPartitioner::SetPartitionParameters "int
Ifpack_GreedyPartitioner::SetPartitionParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner (root node). ";

%feature("docstring")  Ifpack_GreedyPartitioner::ComputePartitions "int Ifpack_GreedyPartitioner::ComputePartitions()

Computes the partitions. Returns 0 if successful. ";


// File: classIfpack__HashTable.xml
%feature("docstring") Ifpack_HashTable "";

%feature("docstring")  Ifpack_HashTable::Ifpack_HashTable "Ifpack_HashTable::Ifpack_HashTable(const int n_keys=1031, const int
n_sets=1)

constructor. ";

%feature("docstring")  Ifpack_HashTable::get "double
Ifpack_HashTable::get(const int key)

Returns an element from the hash table, or 0.0 if not found. ";

%feature("docstring")  Ifpack_HashTable::set "void
Ifpack_HashTable::set(const int key, const double value, const bool
addToValue=false)

Sets an element in the hash table. ";

%feature("docstring")  Ifpack_HashTable::reset "void
Ifpack_HashTable::reset()

Resets the entries of the already allocated memory. This method can be
used to clean an array, to be reused without additional memory
allocation/deallocation. ";

%feature("docstring")  Ifpack_HashTable::getNumEntries "int
Ifpack_HashTable::getNumEntries() const

Returns the number of stored entries. ";

%feature("docstring")  Ifpack_HashTable::arrayify "void
Ifpack_HashTable::arrayify(int *key_array, double *val_array)

Converts the contents in array format for both keys and values. ";

%feature("docstring")  Ifpack_HashTable::print "void
Ifpack_HashTable::print()

Basic printing routine. ";

%feature("docstring")  Ifpack_HashTable::getRecommendedHashSize "int
Ifpack_HashTable::getRecommendedHashSize(int n) ";


// File: classIfpack__IC.xml
%feature("docstring") Ifpack_IC "

Ifpack_IC: A class for constructing and using an incomplete Cholesky
factorization of a given Epetra_RowMatrix.

The Ifpack_IC class computes a threshold based incomplete LDL^T
factorization of a given Epetra_RowMatrix. The factorization that is
produced is a function of several parameters: Maximum number of
entries per row/column in factor - The factorization will contain at
most this number of nonzero terms in each row/column of the
factorization.

Diagonal perturbation - Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing. If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A. Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only. Details can be found
in ifp_diag_pert.

C++ includes: Ifpack_IC.h ";

%feature("docstring")  Ifpack_IC::SetUseTranspose "int
Ifpack_IC::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose_in -If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_IC::NormInf "double
Ifpack_IC::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_IC::HasNormInf "bool
Ifpack_IC::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_IC::UseTranspose "bool
Ifpack_IC::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_IC::OperatorDomainMap "const
Epetra_Map& Ifpack_IC::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_IC::OperatorRangeMap "const Epetra_Map&
Ifpack_IC::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_IC::Comm "const Epetra_Comm&
Ifpack_IC::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_IC::Ifpack_IC "Ifpack_IC::Ifpack_IC(Epetra_RowMatrix *A)

Ifpack_IC constuctor with variable number of indices per row.

Creates a Ifpack_IC object and allocates storage.

Parameters:
-----------

In:  A - User matrix to be factored.

In:  Graph - Graph generated by Ifpack_IlukGraph. ";

%feature("docstring")  Ifpack_IC::~Ifpack_IC "Ifpack_IC::~Ifpack_IC()

Ifpack_IC Destructor. ";

%feature("docstring")  Ifpack_IC::SetAbsoluteThreshold "void
Ifpack_IC::SetAbsoluteThreshold(double Athresh)

Set absolute threshold value. ";

%feature("docstring")  Ifpack_IC::SetRelativeThreshold "void
Ifpack_IC::SetRelativeThreshold(double Rthresh)

Set relative threshold value. ";

%feature("docstring")  Ifpack_IC::SetParameters "int
Ifpack_IC::SetParameters(Teuchos::ParameterList &parameterlis)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_IC::SetParameter "int
Ifpack_IC::SetParameter(const string Name, const int Value) ";

%feature("docstring")  Ifpack_IC::SetParameter "int
Ifpack_IC::SetParameter(const string Name, const double Value) ";

%feature("docstring")  Ifpack_IC::Matrix "const Epetra_RowMatrix&
Ifpack_IC::Matrix() const

Returns a pointer to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_IC::Matrix "Epetra_RowMatrix&
Ifpack_IC::Matrix() ";

%feature("docstring")  Ifpack_IC::IsInitialized "bool
Ifpack_IC::IsInitialized() const

Returns true if the preconditioner has been successfully initialized,
false otherwise. ";

%feature("docstring")  Ifpack_IC::Initialize "int
Ifpack_IC::Initialize()

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_IC::Compute "int Ifpack_IC::Compute()

Compute IC factor U using the specified graph, diagonal perturbation
thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_IC::ComputeSetup "int
Ifpack_IC::ComputeSetup() ";

%feature("docstring")  Ifpack_IC::IsComputed "bool
Ifpack_IC::IsComputed() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_IC::ApplyInverse "int
Ifpack_IC::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_IC forward/back solve on a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_IC::Apply "int Ifpack_IC::Apply(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_IC::Condest "double
Ifpack_IC::Condest(const Ifpack_CondestType CT=Ifpack_Cheap, const int
MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix *Matrix_in=0)

Returns the maximum over all the condition number estimate for each
local ILU set of factors.

This functions computes a local condition number estimate on each
processor and return the maximum over all processor of the estimate.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

Out:  ConditionNumberEstimate - The maximum across all processors of
the infinity-norm estimate of the condition number of the inverse of
LDU. ";

%feature("docstring")  Ifpack_IC::Condest "double
Ifpack_IC::Condest() const

Returns the computed condition number estimate, or -1.0 if not
computed. ";

%feature("docstring")  Ifpack_IC::GetAbsoluteThreshold "double
Ifpack_IC::GetAbsoluteThreshold()

Get absolute threshold value. ";

%feature("docstring")  Ifpack_IC::GetRelativeThreshold "double
Ifpack_IC::GetRelativeThreshold()

Get relative threshold value. ";

%feature("docstring")  Ifpack_IC::NumGlobalNonzeros "int
Ifpack_IC::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_IC::NumMyNonzeros "int
Ifpack_IC::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_IC::D "const Epetra_Vector&
Ifpack_IC::D() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_IC::U "const Epetra_CrsMatrix&
Ifpack_IC::U() const

Returns the address of the U factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_IC::Label "const char*
Ifpack_IC::Label() const ";

%feature("docstring")  Ifpack_IC::SetLabel "int
Ifpack_IC::SetLabel(const char *Label_in) ";

%feature("docstring")  Ifpack_IC::Print "std::ostream &
Ifpack_IC::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_IC::NumInitialize "virtual int
Ifpack_IC::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_IC::NumCompute "virtual int
Ifpack_IC::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_IC::NumApplyInverse "virtual int
Ifpack_IC::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_IC::InitializeTime "virtual double
Ifpack_IC::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_IC::ComputeTime "virtual double
Ifpack_IC::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_IC::ApplyInverseTime "virtual double
Ifpack_IC::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_IC::InitializeFlops "virtual double
Ifpack_IC::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_IC::ComputeFlops "virtual double
Ifpack_IC::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_IC::ApplyInverseFlops "virtual double
Ifpack_IC::ApplyInverseFlops() const

Returns the number of flops in the application of the preconditioner.
";


// File: classIfpack__ICT.xml
%feature("docstring") Ifpack_ICT "

Ifpack_ICT: A class for constructing and using an incomplete Cholesky
factorization of a given Epetra_RowMatrix.

The Ifpack_ICT class computes a threshold based incomplete LDL^T
factorization of a given Epetra_RowMatrix. The factorization that is
produced is a function of several parameters: Maximum number of
entries per row/column in factor - The factorization will contain at
most this number of nonzero terms in each row/column of the
factorization.

Diagonal perturbation - Prior to computing the factorization, it is
possible to modify the diagonal entries of the matrix for which the
factorization will be computing. If the absolute and relative
perturbation values are zero and one, respectively, the factorization
will be compute for the original user matrix A. Otherwise, the
factorization will computed for a matrix that differs from the
original user matrix in the diagonal values only. Details can be found
in ifp_diag_pert.

C++ includes: Ifpack_ICT.h ";

%feature("docstring")  Ifpack_ICT::SetUseTranspose "int
Ifpack_ICT::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose_in -If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_ICT::NormInf "double
Ifpack_ICT::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_ICT::HasNormInf "bool
Ifpack_ICT::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_ICT::UseTranspose "bool
Ifpack_ICT::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_ICT::OperatorDomainMap "const
Epetra_Map& Ifpack_ICT::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_ICT::OperatorRangeMap "const
Epetra_Map& Ifpack_ICT::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_ICT::Comm "const Epetra_Comm&
Ifpack_ICT::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_ICT::Ifpack_ICT "Ifpack_ICT::Ifpack_ICT(const Epetra_RowMatrix *A)

Ifpack_ICT constuctor with variable number of indices per row.

Creates a Ifpack_ICT object and allocates storage.

Parameters:
-----------

In:  A - User matrix to be factored.

In:  Graph - Graph generated by Ifpack_IlukGraph. ";

%feature("docstring")  Ifpack_ICT::~Ifpack_ICT "Ifpack_ICT::~Ifpack_ICT()

Ifpack_ICT Destructor. ";

%feature("docstring")  Ifpack_ICT::SetParameters "int
Ifpack_ICT::SetParameters(Teuchos::ParameterList &parameterlis)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_ICT::Matrix "const Epetra_RowMatrix&
Ifpack_ICT::Matrix() const

Returns a reference to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_ICT::IsInitialized "bool
Ifpack_ICT::IsInitialized() const

Returns true is the preconditioner has been successfully initialized.
";

%feature("docstring")  Ifpack_ICT::Initialize "int
Ifpack_ICT::Initialize()

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_ICT::Compute "int Ifpack_ICT::Compute()

Compute IC factor U using the specified graph, diagonal perturbation
thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_ICT::IsComputed "bool
Ifpack_ICT::IsComputed() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_ICT::ApplyInverse "int
Ifpack_ICT::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_ICT forward/back solve on a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_ICT::Apply "int Ifpack_ICT::Apply(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_ICT::Condest "double
Ifpack_ICT::Condest(const Ifpack_CondestType CT=Ifpack_Cheap, const
int MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix
*Matrix_in=0)

Returns the maximum over all the condition number estimate for each
local ILU set of factors.

This functions computes a local condition number estimate on each
processor and return the maximum over all processor of the estimate.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

Out:  ConditionNumberEstimate - The maximum across all processors of
the infinity-norm estimate of the condition number of the inverse of
LDU. ";

%feature("docstring")  Ifpack_ICT::Condest "double
Ifpack_ICT::Condest() const

Returns the computed condition number estimate, or -1.0 if not
computed. ";

%feature("docstring")  Ifpack_ICT::NumGlobalNonzeros "int
Ifpack_ICT::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_ICT::NumMyNonzeros "int
Ifpack_ICT::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_ICT::H "const Epetra_CrsMatrix&
Ifpack_ICT::H() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_ICT::Label "const char*
Ifpack_ICT::Label() const ";

%feature("docstring")  Ifpack_ICT::SetLabel "int
Ifpack_ICT::SetLabel(const char *Label_in) ";

%feature("docstring")  Ifpack_ICT::Print "std::ostream &
Ifpack_ICT::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_ICT::NumInitialize "virtual int
Ifpack_ICT::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_ICT::NumCompute "virtual int
Ifpack_ICT::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_ICT::NumApplyInverse "virtual int
Ifpack_ICT::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_ICT::InitializeTime "virtual double
Ifpack_ICT::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_ICT::ComputeTime "virtual double
Ifpack_ICT::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_ICT::ApplyInverseTime "virtual double
Ifpack_ICT::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_ICT::InitializeFlops "virtual double
Ifpack_ICT::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_ICT::ComputeFlops "virtual double
Ifpack_ICT::ComputeFlops() const

Returns the number of flops in all applications of Compute(). ";

%feature("docstring")  Ifpack_ICT::ApplyInverseFlops "virtual double
Ifpack_ICT::ApplyInverseFlops() const

Returns the number of flops in all applications of ApplyInverse(). ";

%feature("docstring")  Ifpack_ICT::LevelOfFill "double
Ifpack_ICT::LevelOfFill() const

Returns the level-of-fill.

: if 1.0, then the factored matrix contains approximatively the same
number of elements of A. ";

%feature("docstring")  Ifpack_ICT::AbsoluteThreshold "double
Ifpack_ICT::AbsoluteThreshold() const

Returns the absolute threshold. ";

%feature("docstring")  Ifpack_ICT::RelativeThreshold "double
Ifpack_ICT::RelativeThreshold() const

Returns the relative threshold. ";

%feature("docstring")  Ifpack_ICT::RelaxValue "double
Ifpack_ICT::RelaxValue() const

Returns the relaxation value. ";

%feature("docstring")  Ifpack_ICT::DropTolerance "double
Ifpack_ICT::DropTolerance() const

Returns the drop threshold. ";


// File: classIfpack__IKLU.xml
%feature("docstring") Ifpack_IKLU "

Ifpack_IKLU: A class for constructing and using an incomplete LU
factorization of a given Epetra_RowMatrix.

The Ifpack_IKLU class computes a \"Relaxed\" IKLU factorization with
level k fill of a given Epetra_RowMatrix.

Please refer to ifp_ilu for a general description of the ILU
algorithm.

The complete list of supported parameters is reported in page
ifp_params.

Heidi Thornquist, Org. 1437

C++ includes: Ifpack_IKLU.h ";

%feature("docstring")  Ifpack_IKLU::Ifpack_IKLU "Ifpack_IKLU::Ifpack_IKLU(const Epetra_RowMatrix *A)

Ifpack_IKLU constuctor with variable number of indices per row. ";

%feature("docstring")  Ifpack_IKLU::~Ifpack_IKLU "Ifpack_IKLU::~Ifpack_IKLU()

Ifpack_IKLU Destructor. ";

%feature("docstring")  Ifpack_IKLU::SetParameters "int
Ifpack_IKLU::SetParameters(Teuchos::ParameterList &parameterlis)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_IKLU::Initialize "int
Ifpack_IKLU::Initialize()

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_IKLU::IsInitialized "bool
Ifpack_IKLU::IsInitialized() const

Returns true if the preconditioner has been successfully initialized.
";

%feature("docstring")  Ifpack_IKLU::Compute "int
Ifpack_IKLU::Compute()

Compute IC factor U using the specified graph, diagonal perturbation
thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_IKLU::IsComputed "bool
Ifpack_IKLU::IsComputed() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_IKLU::ApplyInverse "int
Ifpack_IKLU::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_IKLU forward/back solve on a
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to solve for.

Y:  - (Out) A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_IKLU::Apply "int
Ifpack_IKLU::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  Ifpack_IKLU::Condest "double
Ifpack_IKLU::Condest(const Ifpack_CondestType CT=Ifpack_Cheap, const
int MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix
*Matrix_in=0)

Computed the estimated condition number and returns the value. ";

%feature("docstring")  Ifpack_IKLU::Condest "double
Ifpack_IKLU::Condest() const

Returns the computed estimated condition number, or -1.0 if no
computed. ";

%feature("docstring")  Ifpack_IKLU::SetUseTranspose "int
Ifpack_IKLU::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

UseTranspose_in:  - (In) If true, multiply by the transpose of
operator, otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_IKLU::NormInf "double
Ifpack_IKLU::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_IKLU::HasNormInf "bool
Ifpack_IKLU::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_IKLU::UseTranspose "bool
Ifpack_IKLU::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_IKLU::OperatorDomainMap "const
Epetra_Map& Ifpack_IKLU::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_IKLU::OperatorRangeMap "const
Epetra_Map& Ifpack_IKLU::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_IKLU::Comm "const Epetra_Comm&
Ifpack_IKLU::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_IKLU::Matrix "const Epetra_RowMatrix&
Ifpack_IKLU::Matrix() const

Returns a reference to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_IKLU::L "const Epetra_CrsMatrix&
Ifpack_IKLU::L() const

Returns a reference to the L factor. ";

%feature("docstring")  Ifpack_IKLU::U "const Epetra_CrsMatrix&
Ifpack_IKLU::U() const

Returns a reference to the U factor. ";

%feature("docstring")  Ifpack_IKLU::Label "const char*
Ifpack_IKLU::Label() const

Returns the label of this object. ";

%feature("docstring")  Ifpack_IKLU::SetLabel "int
Ifpack_IKLU::SetLabel(const char *Label_in)

Sets the label for this object. ";

%feature("docstring")  Ifpack_IKLU::Print "std::ostream &
Ifpack_IKLU::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_IKLU::NumInitialize "virtual int
Ifpack_IKLU::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_IKLU::NumCompute "virtual int
Ifpack_IKLU::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_IKLU::NumApplyInverse "virtual int
Ifpack_IKLU::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_IKLU::InitializeTime "virtual double
Ifpack_IKLU::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_IKLU::ComputeTime "virtual double
Ifpack_IKLU::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_IKLU::ApplyInverseTime "virtual double
Ifpack_IKLU::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_IKLU::InitializeFlops "virtual double
Ifpack_IKLU::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_IKLU::ComputeFlops "virtual double
Ifpack_IKLU::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_IKLU::ApplyInverseFlops "virtual double
Ifpack_IKLU::ApplyInverseFlops() const

Returns the number of flops in the application of the preconditioner.
";

%feature("docstring")  Ifpack_IKLU::LevelOfFill "double
Ifpack_IKLU::LevelOfFill() const ";

%feature("docstring")  Ifpack_IKLU::RelaxValue "double
Ifpack_IKLU::RelaxValue() const

Set relative threshold value. ";

%feature("docstring")  Ifpack_IKLU::AbsoluteThreshold "double
Ifpack_IKLU::AbsoluteThreshold() const

Get absolute threshold value. ";

%feature("docstring")  Ifpack_IKLU::RelativeThreshold "double
Ifpack_IKLU::RelativeThreshold() const

Get relative threshold value. ";

%feature("docstring")  Ifpack_IKLU::DropTolerance "double
Ifpack_IKLU::DropTolerance() const

Gets the dropping tolerance. ";

%feature("docstring")  Ifpack_IKLU::NumGlobalNonzeros "int
Ifpack_IKLU::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_IKLU::NumMyNonzeros "int
Ifpack_IKLU::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";


// File: classIfpack__ILU.xml
%feature("docstring") Ifpack_ILU "

Ifpack_ILU: A class for constructing and using an incomplete
lower/upper (ILU) factorization of a given Epetra_RowMatrix.

The Ifpack_ILU class computes a \"Relaxed\" ILU factorization with
level k fill of a given Epetra_RowMatrix.

Please refer to ifp_ilu for a general description of the ILU
algorithm.

The complete list of supported parameters is reported in page
ifp_params.

Mike Heroux, Marzio Sala, SNL 9214.

C++ includes: Ifpack_ILU.h ";

%feature("docstring")  Ifpack_ILU::Ifpack_ILU "Ifpack_ILU::Ifpack_ILU(Epetra_RowMatrix *A)

Constructor. ";

%feature("docstring")  Ifpack_ILU::~Ifpack_ILU "Ifpack_ILU::~Ifpack_ILU()

Destructor. ";

%feature("docstring")  Ifpack_ILU::Initialize "int
Ifpack_ILU::Initialize()

Initialize the preconditioner, does not touch matrix values. ";

%feature("docstring")  Ifpack_ILU::IsInitialized "bool
Ifpack_ILU::IsInitialized() const

Returns true if the preconditioner has been successfully initialized.
";

%feature("docstring")  Ifpack_ILU::Compute "int Ifpack_ILU::Compute()

Compute ILU factors L and U using the specified graph, diagonal
perturbation thresholds and relaxation parameters.

This function computes the ILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the ILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_ILU::IsComputed "bool
Ifpack_ILU::IsComputed() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_ILU::SetParameters "int
Ifpack_ILU::SetParameters(Teuchos::ParameterList &parameterlist)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_ILU::SetUseTranspose "int
Ifpack_ILU::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

UseTranspose_in:  - (In) If true, multiply by the transpose of
operator, otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_ILU::Apply "int Ifpack_ILU::Apply(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_ILU::Multiply "int
Ifpack_ILU::Multiply(bool Trans, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_ILU::ApplyInverse "int
Ifpack_ILU::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

In this implementation, we use several existing attributes to
determine how virtual method ApplyInverse() should call the concrete
method Solve(). We pass in the UpperTriangular(), the
Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most
notable warning is that if a matrix has no diagonal values we assume
that there is an implicit unit diagonal that should be accounted for
when doing a triangular solve.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y - (Out) A Epetra_MultiVector of dimension NumVectors
containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_ILU::Condest "double
Ifpack_ILU::Condest(const Ifpack_CondestType CT=Ifpack_Cheap, const
int MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix
*Matrix_in=0)

Computes the estimated condition number and returns the value. ";

%feature("docstring")  Ifpack_ILU::Condest "double
Ifpack_ILU::Condest() const

Returns the computed estimated condition number, or -1.0 if not
computed. ";

%feature("docstring")  Ifpack_ILU::L "const Epetra_CrsMatrix&
Ifpack_ILU::L() const

Returns the address of the L factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_ILU::D "const Epetra_Vector&
Ifpack_ILU::D() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_ILU::U "const Epetra_CrsMatrix&
Ifpack_ILU::U() const

Returns the address of the L factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_ILU::Label "const char*
Ifpack_ILU::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_ILU::SetLabel "int
Ifpack_ILU::SetLabel(const char *Label_in)

Sets label for this object. ";

%feature("docstring")  Ifpack_ILU::NormInf "double
Ifpack_ILU::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_ILU::HasNormInf "bool
Ifpack_ILU::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_ILU::UseTranspose "bool
Ifpack_ILU::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_ILU::OperatorDomainMap "const
Epetra_Map& Ifpack_ILU::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_ILU::OperatorRangeMap "const
Epetra_Map& Ifpack_ILU::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_ILU::Comm "const Epetra_Comm&
Ifpack_ILU::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_ILU::Matrix "const Epetra_RowMatrix&
Ifpack_ILU::Matrix() const

Returns a reference to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_ILU::Print "virtual ostream&
Ifpack_ILU::Print(ostream &os) const

Prints on stream basic information about this object. ";

%feature("docstring")  Ifpack_ILU::NumInitialize "virtual int
Ifpack_ILU::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_ILU::NumCompute "virtual int
Ifpack_ILU::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_ILU::NumApplyInverse "virtual int
Ifpack_ILU::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_ILU::InitializeTime "virtual double
Ifpack_ILU::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_ILU::ComputeTime "virtual double
Ifpack_ILU::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_ILU::ApplyInverseTime "virtual double
Ifpack_ILU::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_ILU::InitializeFlops "virtual double
Ifpack_ILU::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_ILU::ComputeFlops "virtual double
Ifpack_ILU::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_ILU::ApplyInverseFlops "virtual double
Ifpack_ILU::ApplyInverseFlops() const

Returns the number of flops in the application of the preconditioner.
";


// File: classIfpack__IlukGraph.xml
%feature("docstring") Ifpack_IlukGraph "

Ifpack_IlukGraph: A class for constructing level filled graphs for use
with ILU(k) class preconditioners.

The Ifpack_IlukGraph class enable the construction matrix graphs using
level-fill algorithms. The only function required for construction is
an ExtractRowView capability, i.e., the matrix that is passed in to
the constructor must implement the Ifpack_CrsGraph interface defined
in Ifpack_CrsMatrix.h

Constructing Ifpack_IlukGraph objects

Constructing Ifpack_IlukGraph objects is usually a two step process of
passing in a Ifpack_CrsGraph object and an integer indicating the
desired level of fill and then calling the ConstructFilledGraph
function to complete the process. This allows warning error codes to
be returned to the calling routine.

It is worth noting that an Ifpack_IlukGraph object has two
Epetra_CrsGraph objects containing L and U, the graphs for the lower
and upper triangular parts of the ILU(k) graph. Thus, it is possible
to manually insert and delete graph entries in L and U via the
Epetra_CrsGraph InsertIndices and RemoveIndices functions. However, in
this case FillComplete must be called before the graph is used for
subsequent operations.

C++ includes: Ifpack_IlukGraph.h ";

%feature("docstring")  Ifpack_IlukGraph::Ifpack_IlukGraph "Ifpack_IlukGraph::Ifpack_IlukGraph(const Epetra_CrsGraph &Graph_in,
int LevelFill_in, int LevelOverlap_in)

Ifpack_IlukGraph constuctor.

Creates a Ifpack_IlukGraph object using the input graph and specified
level of fill.

Parameters:
-----------

In:  Graph_in - An existing Ifpack_CrsGraph. This object must
implement the Ifpack_CrsGraph functions that provide graph dimension
and pattern information.

In:  LevelFill_in - The level of fill to compute via ILU(k) algorithm.

In:  LevelOverlap_in - The level of between subdomains.

WARNING:  Actual construction occurs in ConstructFilledGraph. This
allows error codes to be passed back to the user. ";

%feature("docstring")  Ifpack_IlukGraph::Ifpack_IlukGraph "Ifpack_IlukGraph::Ifpack_IlukGraph(const Ifpack_IlukGraph &Graph_in)

Copy constructor. ";

%feature("docstring")  Ifpack_IlukGraph::~Ifpack_IlukGraph "Ifpack_IlukGraph::~Ifpack_IlukGraph()

Ifpack_IlukGraph Destructor. ";

%feature("docstring")  Ifpack_IlukGraph::SetParameters "int
Ifpack_IlukGraph::SetParameters(const Teuchos::ParameterList
&parameterlist, bool cerr_warning_if_unused=false)

Set parameters using Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_IlukGraph::ConstructFilledGraph "int
Ifpack_IlukGraph::ConstructFilledGraph()

Does the actual construction of the graph. ";

%feature("docstring")  Ifpack_IlukGraph::ConstructOverlapGraph "int
Ifpack_IlukGraph::ConstructOverlapGraph()

Does the actual construction of the overlap matrix graph. ";

%feature("docstring")  Ifpack_IlukGraph::LevelFill "virtual int
Ifpack_IlukGraph::LevelFill() const

Returns the level of fill used to construct this graph. ";

%feature("docstring")  Ifpack_IlukGraph::LevelOverlap "virtual int
Ifpack_IlukGraph::LevelOverlap() const

Returns the level of overlap used to construct this graph. ";

%feature("docstring")  Ifpack_IlukGraph::NumGlobalBlockRows "int
Ifpack_IlukGraph::NumGlobalBlockRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Ifpack_IlukGraph::NumGlobalBlockCols "int
Ifpack_IlukGraph::NumGlobalBlockCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Ifpack_IlukGraph::NumGlobalRows "int
Ifpack_IlukGraph::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Ifpack_IlukGraph::NumGlobalCols "int
Ifpack_IlukGraph::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Ifpack_IlukGraph::NumGlobalNonzeros "int
Ifpack_IlukGraph::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_IlukGraph::NumGlobalBlockDiagonals "virtual int Ifpack_IlukGraph::NumGlobalBlockDiagonals() const

Returns the number of diagonal entries found in the global input
graph. ";

%feature("docstring")  Ifpack_IlukGraph::NumMyBlockRows "int
Ifpack_IlukGraph::NumMyBlockRows() const

Returns the number of local matrix rows. ";

%feature("docstring")  Ifpack_IlukGraph::NumMyBlockCols "int
Ifpack_IlukGraph::NumMyBlockCols() const

Returns the number of local matrix columns. ";

%feature("docstring")  Ifpack_IlukGraph::NumMyRows "int
Ifpack_IlukGraph::NumMyRows() const

Returns the number of local matrix rows. ";

%feature("docstring")  Ifpack_IlukGraph::NumMyCols "int
Ifpack_IlukGraph::NumMyCols() const

Returns the number of local matrix columns. ";

%feature("docstring")  Ifpack_IlukGraph::NumMyNonzeros "int
Ifpack_IlukGraph::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_IlukGraph::NumMyBlockDiagonals "virtual
int Ifpack_IlukGraph::NumMyBlockDiagonals() const

Returns the number of diagonal entries found in the local input graph.
";

%feature("docstring")  Ifpack_IlukGraph::IndexBase "int
Ifpack_IlukGraph::IndexBase() const

Returns the index base for row and column indices for this graph. ";

%feature("docstring")  Ifpack_IlukGraph::L_Graph "virtual
Epetra_CrsGraph& Ifpack_IlukGraph::L_Graph()

Returns the graph of lower triangle of the ILU(k) graph as a
Epetra_CrsGraph. ";

%feature("docstring")  Ifpack_IlukGraph::U_Graph "virtual
Epetra_CrsGraph& Ifpack_IlukGraph::U_Graph()

Returns the graph of upper triangle of the ILU(k) graph as a
Epetra_CrsGraph. ";

%feature("docstring")  Ifpack_IlukGraph::L_Graph "virtual
Epetra_CrsGraph& Ifpack_IlukGraph::L_Graph() const

Returns the graph of lower triangle of the ILU(k) graph as a
Epetra_CrsGraph. ";

%feature("docstring")  Ifpack_IlukGraph::U_Graph "virtual
Epetra_CrsGraph& Ifpack_IlukGraph::U_Graph() const

Returns the graph of upper triangle of the ILU(k) graph as a
Epetra_CrsGraph. ";

%feature("docstring")  Ifpack_IlukGraph::OverlapImporter "virtual
Epetra_Import* Ifpack_IlukGraph::OverlapImporter() const

Returns the importer used to create the overlapped graph. ";

%feature("docstring")  Ifpack_IlukGraph::OverlapGraph "virtual
Epetra_CrsGraph* Ifpack_IlukGraph::OverlapGraph() const

Returns the the overlapped graph. ";

%feature("docstring")  Ifpack_IlukGraph::DomainMap "virtual const
Epetra_BlockMap& Ifpack_IlukGraph::DomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  Ifpack_IlukGraph::RangeMap "virtual const
Epetra_BlockMap& Ifpack_IlukGraph::RangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_IlukGraph::Comm "virtual const
Epetra_Comm& Ifpack_IlukGraph::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";


// File: classIfpack__ILUT.xml
%feature("docstring") Ifpack_ILUT "

Ifpack_ILUT: A class for constructing and using an incomplete LU
factorization of a given Epetra_RowMatrix.

The Ifpack_ILUT class computes a \"Relaxed\" ILUT factorization with
dual threshold dropping of small elements of a given Epetra_RowMatrix.

This implementation does not use the algorithm that is described in
ifp_ilu. The algorithm drops entries in a row (i) of matrix A that are
smaller than drop_tolerance even before the factorization of row i
then computes the factorization for that row. This is different than
the usual algorithm where the drop tolerance is applied to the
factored rows.

The complete list of supported parameters is reported in page
ifp_params.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_ILUT.h ";

%feature("docstring")  Ifpack_ILUT::Ifpack_ILUT "Ifpack_ILUT::Ifpack_ILUT(const Epetra_RowMatrix *A)

Ifpack_ILUT constuctor with variable number of indices per row. ";

%feature("docstring")  Ifpack_ILUT::~Ifpack_ILUT "Ifpack_ILUT::~Ifpack_ILUT()

Ifpack_ILUT Destructor. ";

%feature("docstring")  Ifpack_ILUT::SetParameters "int
Ifpack_ILUT::SetParameters(Teuchos::ParameterList &parameterlis)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_ILUT::Initialize "int
Ifpack_ILUT::Initialize()

Initialize L and U with values from user matrix A.

Copies values from the user's matrix into the nonzero pattern of L and
U.

Parameters:
-----------

In:  A - User matrix to be factored.

WARNING:  The graph of A must be identical to the graph passed in to
Ifpack_IlukGraph constructor. ";

%feature("docstring")  Ifpack_ILUT::IsInitialized "bool
Ifpack_ILUT::IsInitialized() const

Returns true if the preconditioner has been successfully initialized.
";

%feature("docstring")  Ifpack_ILUT::Compute "int
Ifpack_ILUT::Compute()

Compute IC factor U using the specified graph, diagonal perturbation
thresholds and relaxation parameters.

This function computes the RILU(k) factors L and U using the current:
Ifpack_IlukGraph specifying the structure of L and U.

Value for the RILU(k) relaxation parameter.

Value for the a priori diagonal threshold values.  InitValues() must
be called before the factorization can proceed. ";

%feature("docstring")  Ifpack_ILUT::IsComputed "bool
Ifpack_ILUT::IsComputed() const

If factor is completed, this query returns true, otherwise it returns
false. ";

%feature("docstring")  Ifpack_ILUT::ApplyInverse "int
Ifpack_ILUT::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Ifpack_ILUT forward/back solve on a
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to solve for.

Y:  - (Out) A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_ILUT::Apply "int
Ifpack_ILUT::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  Ifpack_ILUT::Condest "double
Ifpack_ILUT::Condest(const Ifpack_CondestType CT=Ifpack_Cheap, const
int MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix
*Matrix_in=0)

Computed the estimated condition number and returns the value. ";

%feature("docstring")  Ifpack_ILUT::Condest "double
Ifpack_ILUT::Condest() const

Returns the computed estimated condition number, or -1.0 if no
computed. ";

%feature("docstring")  Ifpack_ILUT::SetUseTranspose "int
Ifpack_ILUT::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

UseTranspose_in:  - (In) If true, multiply by the transpose of
operator, otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_ILUT::NormInf "double
Ifpack_ILUT::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_ILUT::HasNormInf "bool
Ifpack_ILUT::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_ILUT::UseTranspose "bool
Ifpack_ILUT::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_ILUT::OperatorDomainMap "const
Epetra_Map& Ifpack_ILUT::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_ILUT::OperatorRangeMap "const
Epetra_Map& Ifpack_ILUT::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_ILUT::Comm "const Epetra_Comm&
Ifpack_ILUT::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Ifpack_ILUT::Matrix "const Epetra_RowMatrix&
Ifpack_ILUT::Matrix() const

Returns a reference to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_ILUT::L "const Epetra_CrsMatrix&
Ifpack_ILUT::L() const

Returns a reference to the L factor. ";

%feature("docstring")  Ifpack_ILUT::U "const Epetra_CrsMatrix&
Ifpack_ILUT::U() const

Returns a reference to the U factor. ";

%feature("docstring")  Ifpack_ILUT::Label "const char*
Ifpack_ILUT::Label() const

Returns the label of this object. ";

%feature("docstring")  Ifpack_ILUT::SetLabel "int
Ifpack_ILUT::SetLabel(const char *Label_in)

Sets the label for this object. ";

%feature("docstring")  Ifpack_ILUT::Print "std::ostream &
Ifpack_ILUT::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_ILUT::NumInitialize "virtual int
Ifpack_ILUT::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_ILUT::NumCompute "virtual int
Ifpack_ILUT::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_ILUT::NumApplyInverse "virtual int
Ifpack_ILUT::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_ILUT::InitializeTime "virtual double
Ifpack_ILUT::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_ILUT::ComputeTime "virtual double
Ifpack_ILUT::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_ILUT::ApplyInverseTime "virtual double
Ifpack_ILUT::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_ILUT::InitializeFlops "virtual double
Ifpack_ILUT::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_ILUT::ComputeFlops "virtual double
Ifpack_ILUT::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_ILUT::ApplyInverseFlops "virtual double
Ifpack_ILUT::ApplyInverseFlops() const

Returns the number of flops in the application of the preconditioner.
";

%feature("docstring")  Ifpack_ILUT::LevelOfFill "double
Ifpack_ILUT::LevelOfFill() const ";

%feature("docstring")  Ifpack_ILUT::RelaxValue "double
Ifpack_ILUT::RelaxValue() const

Set relative threshold value. ";

%feature("docstring")  Ifpack_ILUT::AbsoluteThreshold "double
Ifpack_ILUT::AbsoluteThreshold() const

Get absolute threshold value. ";

%feature("docstring")  Ifpack_ILUT::RelativeThreshold "double
Ifpack_ILUT::RelativeThreshold() const

Get relative threshold value. ";

%feature("docstring")  Ifpack_ILUT::DropTolerance "double
Ifpack_ILUT::DropTolerance() const

Gets the dropping tolerance. ";

%feature("docstring")  Ifpack_ILUT::NumGlobalNonzeros "int
Ifpack_ILUT::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_ILUT::NumMyNonzeros "int
Ifpack_ILUT::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";


// File: classIfpack__LinearPartitioner.xml
%feature("docstring") Ifpack_LinearPartitioner "

Ifpack_LinearPartitioner: A class to define linear partitions.

C++ includes: Ifpack_LinearPartitioner.h ";

%feature("docstring")
Ifpack_LinearPartitioner::Ifpack_LinearPartitioner "Ifpack_LinearPartitioner::Ifpack_LinearPartitioner(const Ifpack_Graph
*Graph)

Constructor. ";

%feature("docstring")
Ifpack_LinearPartitioner::~Ifpack_LinearPartitioner "virtual
Ifpack_LinearPartitioner::~Ifpack_LinearPartitioner()

Destructor. ";

%feature("docstring")
Ifpack_LinearPartitioner::SetPartitionParameters "int
Ifpack_LinearPartitioner::SetPartitionParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner (none for linear
partioning). ";

%feature("docstring")  Ifpack_LinearPartitioner::ComputePartitions "int Ifpack_LinearPartitioner::ComputePartitions()

Computes the partitions. Returns 0 if successful. ";


// File: classIfpack__LocalFilter.xml
%feature("docstring") Ifpack_LocalFilter "

Ifpack_LocalFilter a class for light-weight extraction of the
submatrix corresponding to local rows and columns.

Class Ifpack_LocalFilter enables a light-weight contruction of an
Epetra_RowMatrix-derived object, containing only the elements of the
original, distributed matrix with local row and column ID. The local
submatrix is based on a communicator containing the local process
only. Each process will have its local object, corresponding to the
local submatrix. Submatrices may or may not overlap.

The following instructions can be used to create \"localized\"
matrices:

Once created, LocalA defined, on each process, the submatrix
corresponding to local rows and columns only. The creation and use of
LocalA is \"cheap\", as the elements of the local matrix are obtained
through calls to ExtractMyRowCopy on the original, distributed matrix,
say A. This means that A must remain in scope every time LocalA is
accessed.

A very convenient use of this class is to use Ifpack solvers to
compute the LU factorizations of local blocks. If applied to a
localized matrix, several Ifpack objects can operator in the same
phase in a safe way, without non- required data exchange.

Marzio Sala, SNL 9214

C++ includes: Ifpack_LocalFilter.h ";

%feature("docstring")  Ifpack_LocalFilter::Ifpack_LocalFilter "Ifpack_LocalFilter::Ifpack_LocalFilter(const Teuchos::RefCountPtr<
const Epetra_RowMatrix > &Matrix)

Constructor. ";

%feature("docstring")  Ifpack_LocalFilter::~Ifpack_LocalFilter "virtual Ifpack_LocalFilter::~Ifpack_LocalFilter()

Destructor. ";

%feature("docstring")  Ifpack_LocalFilter::NumMyRowEntries "virtual
int Ifpack_LocalFilter::NumMyRowEntries(int MyRow, int &NumEntries)
const

Returns the number of nonzero entries in MyRow.

Parameters:
-----------

MyRow:  - (In) Local row.

NumEntries:  - (Out) Number of nonzero values present.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_LocalFilter::MaxNumEntries "virtual int
Ifpack_LocalFilter::MaxNumEntries() const

Returns the maximum of NumMyRowEntries() over all rows. ";

%feature("docstring")  Ifpack_LocalFilter::ExtractMyRowCopy "int
Ifpack_LocalFilter::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

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

%feature("docstring")  Ifpack_LocalFilter::ExtractDiagonalCopy "int
Ifpack_LocalFilter::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Diagonal:  - (Out) Extracted main diagonal.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_LocalFilter::Multiply "virtual int
Ifpack_LocalFilter::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

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

%feature("docstring")  Ifpack_LocalFilter::Solve "virtual int
Ifpack_LocalFilter::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns result of a local-only solve using a triangular
Epetra_RowMatrix with Epetra_MultiVectors X and Y (NOT IMPLEMENTED).
";

%feature("docstring")  Ifpack_LocalFilter::Apply "int
Ifpack_LocalFilter::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_LocalFilter::ApplyInverse "int
Ifpack_LocalFilter::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_LocalFilter::InvRowSums "virtual int
Ifpack_LocalFilter::InvRowSums(Epetra_Vector &x) const

Computes the sum of absolute values of the rows of the
Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED). ";

%feature("docstring")  Ifpack_LocalFilter::LeftScale "virtual int
Ifpack_LocalFilter::LeftScale(const Epetra_Vector &x)

Scales the Epetra_RowMatrix on the left with a Epetra_Vector x (NOT
IMPLEMENTED). ";

%feature("docstring")  Ifpack_LocalFilter::InvColSums "virtual int
Ifpack_LocalFilter::InvColSums(Epetra_Vector &x) const

Computes the sum of absolute values of the columns of the
Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED). ";

%feature("docstring")  Ifpack_LocalFilter::RightScale "virtual int
Ifpack_LocalFilter::RightScale(const Epetra_Vector &x)

Scales the Epetra_RowMatrix on the right with a Epetra_Vector x (NOT
IMPLEMENTED). ";

%feature("docstring")  Ifpack_LocalFilter::Filled "virtual bool
Ifpack_LocalFilter::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  Ifpack_LocalFilter::NormInf "virtual double
Ifpack_LocalFilter::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Ifpack_LocalFilter::NormOne "virtual double
Ifpack_LocalFilter::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")  Ifpack_LocalFilter::NumGlobalNonzeros "virtual
int Ifpack_LocalFilter::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Ifpack_LocalFilter::NumGlobalRows "virtual int
Ifpack_LocalFilter::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Ifpack_LocalFilter::NumGlobalCols "virtual int
Ifpack_LocalFilter::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Ifpack_LocalFilter::NumGlobalDiagonals "virtual int Ifpack_LocalFilter::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Ifpack_LocalFilter::NumMyNonzeros "virtual int
Ifpack_LocalFilter::NumMyNonzeros() const

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Ifpack_LocalFilter::NumMyRows "virtual int
Ifpack_LocalFilter::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Ifpack_LocalFilter::NumMyCols "virtual int
Ifpack_LocalFilter::NumMyCols() const

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")  Ifpack_LocalFilter::NumMyDiagonals "virtual
int Ifpack_LocalFilter::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Ifpack_LocalFilter::LowerTriangular "virtual
bool Ifpack_LocalFilter::LowerTriangular() const

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Ifpack_LocalFilter::UpperTriangular "virtual
bool Ifpack_LocalFilter::UpperTriangular() const

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Ifpack_LocalFilter::RowMatrixRowMap "virtual
const Epetra_Map& Ifpack_LocalFilter::RowMatrixRowMap() const

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")  Ifpack_LocalFilter::RowMatrixColMap "virtual
const Epetra_Map& Ifpack_LocalFilter::RowMatrixColMap() const

Returns the Epetra_Map object associated with the columns of this
matrix. ";

%feature("docstring")  Ifpack_LocalFilter::RowMatrixImporter "virtual
const Epetra_Import* Ifpack_LocalFilter::RowMatrixImporter() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")  Ifpack_LocalFilter::SetOwnership "int
Ifpack_LocalFilter::SetOwnership(bool ownership)

Sets ownership. ";

%feature("docstring")  Ifpack_LocalFilter::SetUseTranspose "int
Ifpack_LocalFilter::SetUseTranspose(bool UseTranspose_in)

Sets use transpose (not implemented). ";

%feature("docstring")  Ifpack_LocalFilter::UseTranspose "bool
Ifpack_LocalFilter::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_LocalFilter::HasNormInf "bool
Ifpack_LocalFilter::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_LocalFilter::Comm "const Epetra_Comm&
Ifpack_LocalFilter::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_LocalFilter::OperatorDomainMap "const
Epetra_Map& Ifpack_LocalFilter::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_LocalFilter::OperatorRangeMap "const
Epetra_Map& Ifpack_LocalFilter::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_LocalFilter::Map "const Epetra_BlockMap
& Ifpack_LocalFilter::Map() const ";

%feature("docstring")  Ifpack_LocalFilter::Label "const char*
Ifpack_LocalFilter::Label() const ";


// File: classIfpack__METISPartitioner.xml
%feature("docstring") Ifpack_METISPartitioner "

Ifpack_METISPartitioner: A class to decompose Ifpack_Graph's using
METIS.

Class Ifpack_METISPartitioner enables the decomposition of the local
Ifpack_Graph's using METIS. In order to work properly, this class
requires IFPACK to be configured with option --enable-ifpack-metis.
Otherwise, this class will always create one partition.

C++ includes: Ifpack_METISPartitioner.h ";

%feature("docstring")
Ifpack_METISPartitioner::Ifpack_METISPartitioner "Ifpack_METISPartitioner::Ifpack_METISPartitioner(const Ifpack_Graph
*Graph)

Constructor. ";

%feature("docstring")
Ifpack_METISPartitioner::~Ifpack_METISPartitioner "virtual
Ifpack_METISPartitioner::~Ifpack_METISPartitioner()

Destructor. ";

%feature("docstring")  Ifpack_METISPartitioner::SetPartitionParameters
"int
Ifpack_METISPartitioner::SetPartitionParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner (none at moment). ";

%feature("docstring")  Ifpack_METISPartitioner::ComputePartitions "int Ifpack_METISPartitioner::ComputePartitions()

Computes the partitions. Returns 0 if successful. ";


// File: classIfpack__METISReordering.xml
%feature("docstring") Ifpack_METISReordering "

Ifpack_METISReordering: A class to reorder a graph using METIS.

C++ includes: Ifpack_METISReordering.h ";

%feature("docstring")  Ifpack_METISReordering::Ifpack_METISReordering
"Ifpack_METISReordering::Ifpack_METISReordering()

Constructor. ";

%feature("docstring")  Ifpack_METISReordering::~Ifpack_METISReordering
"virtual Ifpack_METISReordering::~Ifpack_METISReordering()

Destructor. ";

%feature("docstring")  Ifpack_METISReordering::SetParameter "virtual
int Ifpack_METISReordering::SetParameter(const string Name, const int
Value)

Sets integer parameters `Name'. ";

%feature("docstring")  Ifpack_METISReordering::SetParameter "virtual
int Ifpack_METISReordering::SetParameter(const string Name, const
double Value)

Sets double parameters `Name'. ";

%feature("docstring")  Ifpack_METISReordering::SetParameters "virtual
int Ifpack_METISReordering::SetParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner (none at moment). ";

%feature("docstring")  Ifpack_METISReordering::Compute "int
Ifpack_METISReordering::Compute(const Ifpack_Graph &Graph)

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_METISReordering::Compute "int
Ifpack_METISReordering::Compute(const Epetra_RowMatrix &Matrix)

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_METISReordering::IsComputed "virtual
bool Ifpack_METISReordering::IsComputed() const

Returns true is the reordering object has been successfully
initialized, false otherwise. ";

%feature("docstring")  Ifpack_METISReordering::Reorder "int
Ifpack_METISReordering::Reorder(const int i) const

Returns the reordered index of row i. ";

%feature("docstring")  Ifpack_METISReordering::InvReorder "int
Ifpack_METISReordering::InvReorder(const int i) const

Returns the inverse reordered index of row i. ";

%feature("docstring")  Ifpack_METISReordering::P "int
Ifpack_METISReordering::P(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &X) const

Applies reordering to multivector Xorig, whose local length equals the
number of local rows, stores result in X. ";

%feature("docstring")  Ifpack_METISReordering::Pinv "int
Ifpack_METISReordering::Pinv(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &X) const

Applies inverse reordering to multivector Xorig, whose local length
equals the number of local rows, stores result in X. ";

%feature("docstring")  Ifpack_METISReordering::Print "ostream &
Ifpack_METISReordering::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";


// File: classIfpack__OverlapFactorObject.xml
%feature("docstring") Ifpack_OverlapFactorObject "

Ifpack_OverlapFactorObject: Supports functionality common to Ifpack
overlap factorization classes.

C++ includes: Ifpack_OverlapFactorObject.h ";

%feature("docstring")
Ifpack_OverlapFactorObject::Ifpack_OverlapFactorObject "Ifpack_OverlapFactorObject::Ifpack_OverlapFactorObject(const
Ifpack_OverlapGraph *OverlapGraph)

Constructor using Ifpack_OverlapGraph.

Creates an object from the overlap graph.

Parameters:
-----------

In:  OverlapGraph - Graph describing the graph that should be used for
the factors. ";

%feature("docstring")
Ifpack_OverlapFactorObject::Ifpack_OverlapFactorObject "Ifpack_OverlapFactorObject::Ifpack_OverlapFactorObject(const
Epetra_RowMatrix *UserMatrix)

Constructor using Epetra_RowMatrix.

Creates an Ifpack_Graph object from the user graph implicitly defined
by the Epetra_RowMatrix interface.

Parameters:
-----------

In:  RowMatrix - An object that has implemented the Epetra_RowMatrix
interface. ";

%feature("docstring")
Ifpack_OverlapFactorObject::Ifpack_OverlapFactorObject "Ifpack_OverlapFactorObject::Ifpack_OverlapFactorObject(const
Ifpack_OverlapFactorObject &Source)

Copy constructor. ";

%feature("docstring")
Ifpack_OverlapFactorObject::~Ifpack_OverlapFactorObject "virtual
Ifpack_OverlapFactorObject::~Ifpack_OverlapFactorObject()

Ifpack_OverlapFactorObject Destructor. ";

%feature("docstring")  Ifpack_OverlapFactorObject::InitValues "virtual int Ifpack_OverlapFactorObject::InitValues(const
Epetra_RowMatrix *UserMatrix)

Initialize values from user matrix A, can be called repeatedly as
matrix values change.

Processes matrix values, primarily handling overlap if any has been
requested. This method then calls ProcessOverlapMatrix(), a virtual
method that must be implemented by any class that derives from this
class.

Parameters:
-----------

In:  UserMatrix - User matrix to be processed. ";

%feature("docstring")  Ifpack_OverlapFactorObject::Factor "virtual
int Ifpack_OverlapFactorObject::Factor()

Compute factors.

This function computes factors using the method DerivedFactor() that
is implemented by the derived class. InitValues() must be called
before the factorization can proceed. ";

%feature("docstring")  Ifpack_OverlapFactorObject::Allocated "bool
Ifpack_OverlapFactorObject::Allocated() const

If storage has been allocated, this query returns true, otherwise it
returns false. ";

%feature("docstring")  Ifpack_OverlapFactorObject::ValuesInitialized "bool Ifpack_OverlapFactorObject::ValuesInitialized() const

If values have been initialized, this query returns true, otherwise it
returns false. ";

%feature("docstring")  Ifpack_OverlapFactorObject::Factored "bool
Ifpack_OverlapFactorObject::Factored() const

If factor is completed, this query returns true, otherwise it returns
false. ";


// File: classIfpack__OverlapGraph.xml
%feature("docstring") Ifpack_OverlapGraph "

Ifpack_OverlapGraph: Constructs a graph for use with Ifpack
preconditioners.

C++ includes: Ifpack_OverlapGraph.h ";

%feature("docstring")  Ifpack_OverlapGraph::Ifpack_OverlapGraph "Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Teuchos::RefCountPtr<
const Epetra_CrsGraph > &UserMatrixGraph_in, int OverlapLevel_in)

Constructor using Epetra_CrsGraph.

Creates an Ifpack_OverlapGraph object from the user graph.

Parameters:
-----------

In:  UserMatrixGraph_in - Graph from user matrix. ";

%feature("docstring")  Ifpack_OverlapGraph::Ifpack_OverlapGraph "Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Teuchos::RefCountPtr<
const Epetra_RowMatrix > &UserMatrix_in, int OverlapLevel_in)

Constructor using Epetra_RowMatrix.

Creates an Ifpack_OverlapGraph object from the user graph implicitly
defined by the Epetra_RowMatrix interface.

Parameters:
-----------

In:  RowMatrix - An object that has implemented the Epetra_RowMatrix
interface. ";

%feature("docstring")  Ifpack_OverlapGraph::Ifpack_OverlapGraph "Ifpack_OverlapGraph::Ifpack_OverlapGraph(const Ifpack_OverlapGraph
&Source)

Copy constructor. ";

%feature("docstring")  Ifpack_OverlapGraph::~Ifpack_OverlapGraph "virtual Ifpack_OverlapGraph::~Ifpack_OverlapGraph()

Ifpack_CrsIlut Destructor. ";

%feature("docstring")  Ifpack_OverlapGraph::SetParameters "int
Ifpack_OverlapGraph::SetParameters(const Teuchos::ParameterList
&parameterlist, bool cerr_warning_if_unused=false)

Set parameters using a Teuchos::ParameterList object. ";

%feature("docstring")  Ifpack_OverlapGraph::OverlapGraph "const
Epetra_CrsGraph& Ifpack_OverlapGraph::OverlapGraph() const

Returns the overlap graph object. ";

%feature("docstring")  Ifpack_OverlapGraph::OverlapRowMap "const
Epetra_BlockMap& Ifpack_OverlapGraph::OverlapRowMap() const

Returns the RowMap associated with the overlap graph. ";

%feature("docstring")  Ifpack_OverlapGraph::OverlapImporter "const
Epetra_Import& Ifpack_OverlapGraph::OverlapImporter() const

Returns the overlap graph object. ";

%feature("docstring")  Ifpack_OverlapGraph::OverlapLevel "int
Ifpack_OverlapGraph::OverlapLevel() const

Returns the level of overlap used to create this graph.

The graph created by this class uses a recursive definition 0f
overlap. Level one overlap is created by copying all off-processor
rows that are reached to be at least one column of the rows that are
on processor. Level two overlap is the same process used on the level
one graph. ";

%feature("docstring")  Ifpack_OverlapGraph::Print "void
Ifpack_OverlapGraph::Print(ostream &os) const ";


// File: classIfpack__OverlappingPartitioner.xml
%feature("docstring") Ifpack_OverlappingPartitioner "";

%feature("docstring")
Ifpack_OverlappingPartitioner::Ifpack_OverlappingPartitioner "Ifpack_OverlappingPartitioner::Ifpack_OverlappingPartitioner(const
Ifpack_Graph *Graph)

Constructor. ";

%feature("docstring")
Ifpack_OverlappingPartitioner::~Ifpack_OverlappingPartitioner "Ifpack_OverlappingPartitioner::~Ifpack_OverlappingPartitioner()

Destructor. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::NumLocalParts "int Ifpack_OverlappingPartitioner::NumLocalParts() const

Returns the number of computed local partitions. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::OverlappingLevel
"int Ifpack_OverlappingPartitioner::OverlappingLevel() const

Returns the overlapping level. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::NumRowsInPart "int Ifpack_OverlappingPartitioner::NumRowsInPart(const int Part) const

Returns the number of rows contained in specified partition. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::RowsInPart "int
Ifpack_OverlappingPartitioner::RowsInPart(const int Part, int *List)
const

Copies into List the rows in the (overlapping) partition Part. ";

%feature("docstring")
Ifpack_OverlappingPartitioner::NonOverlappingPartition "const int*
Ifpack_OverlappingPartitioner::NonOverlappingPartition() const

Returns a pointer to the integer vector containing the non-overlapping
partition ID of each local row. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::SetParameters "int
Ifpack_OverlappingPartitioner::SetParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner.

The supported parameters are:  \"partitioner: overlap\" (int, default
= 0).

\"partitioner: local parts\" (int, default = 1).

\"partitioner: print level\" (int, default = 0). ";

%feature("docstring")
Ifpack_OverlappingPartitioner::SetPartitionParameters "virtual int
Ifpack_OverlappingPartitioner::SetPartitionParameters(Teuchos::ParameterList
&List)=0

Sets all the parameters for the partitioner.

This function is used by derived classes to set their own parameters.
These classes should not derive SetParameters(), so that common
parameters can be set just once. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::Compute "int
Ifpack_OverlappingPartitioner::Compute()

Computes the partitions. Returns 0 if successful. ";

%feature("docstring")
Ifpack_OverlappingPartitioner::ComputePartitions "virtual int
Ifpack_OverlappingPartitioner::ComputePartitions()=0

Computes the partitions. Returns 0 if successful. ";

%feature("docstring")
Ifpack_OverlappingPartitioner::ComputeOverlappingPartitions "int
Ifpack_OverlappingPartitioner::ComputeOverlappingPartitions()

Computes the partitions. Returns 0 if successful. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::IsComputed "bool Ifpack_OverlappingPartitioner::IsComputed()

Returns true if partitions have been computed successfully. ";

%feature("docstring")  Ifpack_OverlappingPartitioner::Print "virtual
ostream& Ifpack_OverlappingPartitioner::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";


// File: classIfpack__OverlappingRowMatrix.xml
%feature("docstring") Ifpack_OverlappingRowMatrix "

Ifpack_OverlappingRowMatrix: matrix with ghost rows, based on
Epetra_RowMatrix.

C++ includes: Ifpack_OverlappingRowMatrix.h ";

%feature("docstring")
Ifpack_OverlappingRowMatrix::Ifpack_OverlappingRowMatrix "Ifpack_OverlappingRowMatrix::Ifpack_OverlappingRowMatrix(const
Teuchos::RefCountPtr< const Epetra_RowMatrix > &Matrix_in, int
OverlapLevel_in) ";

%feature("docstring")
Ifpack_OverlappingRowMatrix::~Ifpack_OverlappingRowMatrix "Ifpack_OverlappingRowMatrix::~Ifpack_OverlappingRowMatrix() ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumMyRowEntries "int Ifpack_OverlappingRowMatrix::NumMyRowEntries(int MyRow, int
&NumEntries) const

Returns the number of nonzero entries in MyRow.

Parameters:
-----------

MyRow:  - (In) Local row.

NumEntries:  - (Out) Number of nonzero values present.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::MaxNumEntries "virtual int Ifpack_OverlappingRowMatrix::MaxNumEntries() const

Returns the maximum of NumMyRowEntries() over all rows. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::ExtractMyRowCopy "int Ifpack_OverlappingRowMatrix::ExtractMyRowCopy(int MyRow, int
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

%feature("docstring")
Ifpack_OverlappingRowMatrix::ExtractDiagonalCopy "int
Ifpack_OverlappingRowMatrix::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Diagonal:  - (Out) Extracted main diagonal.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::Multiply "int
Ifpack_OverlappingRowMatrix::Multiply(bool TransA, const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const

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

%feature("docstring")  Ifpack_OverlappingRowMatrix::Solve "virtual
int Ifpack_OverlappingRowMatrix::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const

Returns result of a local-only solve using a triangular
Epetra_RowMatrix with Epetra_MultiVectors X and Y (NOT IMPLEMENTED).
";

%feature("docstring")  Ifpack_OverlappingRowMatrix::Apply "int
Ifpack_OverlappingRowMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::ApplyInverse "int
Ifpack_OverlappingRowMatrix::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::InvRowSums "virtual int Ifpack_OverlappingRowMatrix::InvRowSums(Epetra_Vector &x)
const

Computes the sum of absolute values of the rows of the
Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED). ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::LeftScale "virtual int Ifpack_OverlappingRowMatrix::LeftScale(const Epetra_Vector
&x)

Scales the Epetra_RowMatrix on the left with a Epetra_Vector x (NOT
IMPLEMENTED). ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::InvColSums "virtual int Ifpack_OverlappingRowMatrix::InvColSums(Epetra_Vector &x)
const

Computes the sum of absolute values of the columns of the
Epetra_RowMatrix, results returned in x (NOT IMPLEMENTED). ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::RightScale "virtual int Ifpack_OverlappingRowMatrix::RightScale(const
Epetra_Vector &x)

Scales the Epetra_RowMatrix on the right with a Epetra_Vector x (NOT
IMPLEMENTED). ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::Filled "virtual
bool Ifpack_OverlappingRowMatrix::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NormInf "virtual
double Ifpack_OverlappingRowMatrix::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NormOne "virtual
double Ifpack_OverlappingRowMatrix::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumGlobalNonzeros
"virtual int Ifpack_OverlappingRowMatrix::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumGlobalRows "virtual int Ifpack_OverlappingRowMatrix::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumGlobalCols "virtual int Ifpack_OverlappingRowMatrix::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumGlobalDiagonals
"virtual int Ifpack_OverlappingRowMatrix::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumMyNonzeros "virtual int Ifpack_OverlappingRowMatrix::NumMyNonzeros() const

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumMyRows "virtual int Ifpack_OverlappingRowMatrix::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumMyCols "virtual int Ifpack_OverlappingRowMatrix::NumMyCols() const

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")  Ifpack_OverlappingRowMatrix::NumMyDiagonals "virtual int Ifpack_OverlappingRowMatrix::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::LowerTriangular "virtual bool Ifpack_OverlappingRowMatrix::LowerTriangular() const

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::UpperTriangular "virtual bool Ifpack_OverlappingRowMatrix::UpperTriangular() const

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::RowMatrixRowMap "virtual const Epetra_Map&
Ifpack_OverlappingRowMatrix::RowMatrixRowMap() const

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")  Ifpack_OverlappingRowMatrix::RowMatrixColMap "virtual const Epetra_Map&
Ifpack_OverlappingRowMatrix::RowMatrixColMap() const

Returns the Epetra_Map object associated with the columns of this
matrix. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::RowMatrixImporter
"virtual const Epetra_Import*
Ifpack_OverlappingRowMatrix::RowMatrixImporter() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::SetOwnership "int
Ifpack_OverlappingRowMatrix::SetOwnership(bool ownership)

Sets ownership. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::SetUseTranspose "int Ifpack_OverlappingRowMatrix::SetUseTranspose(bool UseTranspose_in)

Sets use transpose (not implemented). ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::UseTranspose "bool Ifpack_OverlappingRowMatrix::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::HasNormInf "bool
Ifpack_OverlappingRowMatrix::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::Comm "const
Epetra_Comm& Ifpack_OverlappingRowMatrix::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::OperatorDomainMap
"const Epetra_Map& Ifpack_OverlappingRowMatrix::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::OperatorRangeMap "const Epetra_Map& Ifpack_OverlappingRowMatrix::OperatorRangeMap()
const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::Map "const
Epetra_BlockMap & Ifpack_OverlappingRowMatrix::Map() const ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::Label "const
char* Ifpack_OverlappingRowMatrix::Label() const ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::OverlapLevel "int
Ifpack_OverlappingRowMatrix::OverlapLevel() const ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::ImportMultiVector
"int Ifpack_OverlappingRowMatrix::ImportMultiVector(const
Epetra_MultiVector &X, Epetra_MultiVector &OvX, Epetra_CombineMode
CM=Insert) ";

%feature("docstring")  Ifpack_OverlappingRowMatrix::ExportMultiVector
"int Ifpack_OverlappingRowMatrix::ExportMultiVector(const
Epetra_MultiVector &OvX, Epetra_MultiVector &X, Epetra_CombineMode
CM=Add) ";


// File: classIfpack__OverlapSolveObject.xml
%feature("docstring") Ifpack_OverlapSolveObject "

Ifpack_OverlapSolveObject: Provides Overlapped Forward/back solve
services for Ifpack.

C++ includes: Ifpack_OverlapSolveObject.h ";

%feature("docstring")
Ifpack_OverlapSolveObject::Ifpack_OverlapSolveObject "Ifpack_OverlapSolveObject::Ifpack_OverlapSolveObject(char *Label,
const Epetra_Comm &Comm)

Constructor. ";

%feature("docstring")
Ifpack_OverlapSolveObject::Ifpack_OverlapSolveObject "Ifpack_OverlapSolveObject::Ifpack_OverlapSolveObject(const
Ifpack_OverlapSolveObject &Source)

Copy constructor. ";

%feature("docstring")
Ifpack_OverlapSolveObject::~Ifpack_OverlapSolveObject "Ifpack_OverlapSolveObject::~Ifpack_OverlapSolveObject()

Ifpack_OverlapSolveObject Destructor. ";

%feature("docstring")  Ifpack_OverlapSolveObject::SetOverlapMode "void Ifpack_OverlapSolveObject::SetOverlapMode(Epetra_CombineMode
OverlapMode)

Generate Ifpack_OverlapGraph object using current settings. ";

%feature("docstring")  Ifpack_OverlapSolveObject::SetLowerOperator "int Ifpack_OverlapSolveObject::SetLowerOperator(Epetra_CrsMatrix *L,
bool UseLTrans)

Define the operator to be used for the lower triangle. ";

%feature("docstring")  Ifpack_OverlapSolveObject::SetDiagonal "int
Ifpack_OverlapSolveObject::SetDiagonal(Epetra_Vector *D, bool UseDInv)

Define the vector to be used for the diagonal. ";

%feature("docstring")  Ifpack_OverlapSolveObject::SetUpperOperator "int Ifpack_OverlapSolveObject::SetUpperOperator(Epetra_CrsMatrix *U,
bool UseUTrans)

Define the operator to be used for the upper triangle. ";

%feature("docstring")  Ifpack_OverlapSolveObject::Solve "int
Ifpack_OverlapSolveObject::Solve(bool Trans, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const

Returns the result of a Ifpack_CrsIlut forward/back solve on a
Epetra_MultiVector X in Y (works for Epetra_Vectors also).

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_OverlapSolveObject::Multiply "int
Ifpack_OverlapSolveObject::Multiply(bool Trans, const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of multiplying U, D and L in that order on an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Trans -If true, multiply by L^T, D and U^T in that order.

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_OverlapSolveObject::Condest "int
Ifpack_OverlapSolveObject::Condest(bool Trans, double
&ConditionNumberEstimate) const

Returns the maximum over all the condition number estimate for each
local ILU set of factors.

This functions computes a local condition number estimate on each
processor and return the maximum over all processor of the estimate.

Parameters:
-----------

In:  Trans -If true, solve transpose problem.

Out:  ConditionNumberEstimate - The maximum across all processors of
the infinity-norm estimate of the condition number of the inverse of
LDU. ";

%feature("docstring")  Ifpack_OverlapSolveObject::OverlapMode "Epetra_CombineMode Ifpack_OverlapSolveObject::OverlapMode() const

Returns the overlap mode used to combine terms that are redundantly
computed.

Since rows of the graph, and any related matrices are multiply owned,
some values in the subdomain solves will be computed on multiple
processors. The overlap mode is used to determine how the redundant
values that come in from other processors will be handled. ";

%feature("docstring")  Ifpack_OverlapSolveObject::NumGlobalNonzeros "int Ifpack_OverlapSolveObject::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global graph. ";

%feature("docstring")  Ifpack_OverlapSolveObject::NumMyNonzeros "int
Ifpack_OverlapSolveObject::NumMyNonzeros() const

Returns the number of nonzero entries in the local graph. ";

%feature("docstring")  Ifpack_OverlapSolveObject::L "const
Epetra_CrsMatrix& Ifpack_OverlapSolveObject::L() const

Returns the address of the L factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_OverlapSolveObject::D "const
Epetra_Vector& Ifpack_OverlapSolveObject::D() const

Returns the address of the D factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_OverlapSolveObject::U "const
Epetra_CrsMatrix& Ifpack_OverlapSolveObject::U() const

Returns the address of the L factor associated with this factored
matrix. ";

%feature("docstring")  Ifpack_OverlapSolveObject::Label "char*
Ifpack_OverlapSolveObject::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Ifpack_OverlapSolveObject::SetUseTranspose "int Ifpack_OverlapSolveObject::SetUseTranspose(bool UseTranspose)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose -If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Ifpack_OverlapSolveObject::Apply "int
Ifpack_OverlapSolveObject::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Note that this implementation of Apply does NOT perform a forward back
solve with the LDU factorization. Instead it applies these operators
via multiplication with U, D and L respectively. The ApplyInverse()
method performs a solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_OverlapSolveObject::ApplyInverse "int
Ifpack_OverlapSolveObject::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

In this implementation, we use several existing attributes to
determine how virtual method ApplyInverse() should call the concrete
method Solve(). We pass in the UpperTriangular(), the
Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most
notable warning is that if a matrix has no diagonal values we assume
that there is an implicit unit diagonal that should be accounted for
when doing a triangular solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_OverlapSolveObject::NormInf "double
Ifpack_OverlapSolveObject::NormInf() const

Returns 0.0 because this class cannot compute Inf-norm. ";

%feature("docstring")  Ifpack_OverlapSolveObject::HasNormInf "bool
Ifpack_OverlapSolveObject::HasNormInf() const

Returns false because this class cannot compute an Inf-norm. ";

%feature("docstring")  Ifpack_OverlapSolveObject::UseTranspose "bool
Ifpack_OverlapSolveObject::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_OverlapSolveObject::OperatorDomainMap "const Epetra_Map& Ifpack_OverlapSolveObject::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_OverlapSolveObject::OperatorRangeMap "const Epetra_Map& Ifpack_OverlapSolveObject::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_OverlapSolveObject::Comm "const
Epetra_Comm& Ifpack_OverlapSolveObject::Comm() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";


// File: classIfpack__Partitioner.xml
%feature("docstring") Ifpack_Partitioner "

Ifpack_Partitioner: A class to decompose local Ifpack_Graph's.

Class Ifpack_Partitioner enables the decomposition of a local
Ifpack_Graph's. It is supposed that the graph refers to a localized
matrix (that is, a matrix that has been filtered through
Ifpack_LocalFilter).

The overloaded operator (int i) can be used to extract the local
partition ID of local row i.

The partitions created by Ifpack_Partitioner derived clased are non-
overlapping in graph sense. This means that each row (or, more
approriately, vertex) of G is assigned to exactly one partition.

Partitioner can be extended using the functionalities of class
Ifpack_OverlappingPartitioner (itself derived from Ifpack_Partitioner.
This class extends the non-overlapping partitions by the required
amount of overlap, considering local nodes only (that is, this overlap
do not modify the overlap among the processes).

Ifpack_Partitioner is a pure virtual class. Concrete implementations
are:  Ifpack_LinearPartitioner, which allows the decomposition of the
rows of the graph in simple consecutive chunks;

Ifpack_METISPartitioner, which calls METIS to decompose the graph
(this requires the configuration option --enable-ifpack-metis);

Ifpack_GreedyPartitioner, a simple greedy algorith;

Ifpack_EquationPartitioner, which creates NumPDEEqns parts (where
NumPDEEqns is the number of equations in the linear system). It is
supposed that all the equations referring to the same grid node are
ordered consecutively. Besides, the number of equations per node must
be constant in the domain.

Generically, a constructor requires an Ifpack_Graph object.
Ifpack_Graph is a pure virtual class. Concrete implentations are:
Ifpack_Graph_Epetra_CrsGraph, a light-weight class to wrap
Epetra_CrsGraph objects as Ifpack_Graph objects;

Ifpack_Graph_Epetra_RowMatrix, a light-weight class to wrap
Epetra_RowMatrix objects as Ifpack_Graph objects.

An example of use is an Ifpack_Partitioner derived class is as
follows:

When overlapping partitiones are created, the user can get the row ID
contained in each partition as follows:

Ifpack_Partitioner is used to create the subblocks in
Ifpack_BlockJacobi, Ifpack_BlockGaussSeidel, and
Ifpack_BlockSymGaussSeidel.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_Partitioner.h ";

%feature("docstring")  Ifpack_Partitioner::~Ifpack_Partitioner "virtual Ifpack_Partitioner::~Ifpack_Partitioner()

Destructor. ";

%feature("docstring")  Ifpack_Partitioner::NumLocalParts "virtual int
Ifpack_Partitioner::NumLocalParts() const =0

Returns the number of computed local partitions. ";

%feature("docstring")  Ifpack_Partitioner::OverlappingLevel "virtual
int Ifpack_Partitioner::OverlappingLevel() const =0

Returns the overlapping level. ";

%feature("docstring")  Ifpack_Partitioner::NumRowsInPart "virtual int
Ifpack_Partitioner::NumRowsInPart(const int Part) const =0

Returns the number of rows contained in specified partition. ";

%feature("docstring")  Ifpack_Partitioner::RowsInPart "virtual int
Ifpack_Partitioner::RowsInPart(const int Part, int *List) const =0

Copies into List the rows in the (overlapping) partition Part. ";

%feature("docstring")  Ifpack_Partitioner::NonOverlappingPartition "virtual const int* Ifpack_Partitioner::NonOverlappingPartition() const
=0

Returns a pointer to the integer vector containing the non-overlapping
partition ID of each local row. ";

%feature("docstring")  Ifpack_Partitioner::SetParameters "virtual int
Ifpack_Partitioner::SetParameters(Teuchos::ParameterList &List)=0

Sets all the parameters for the partitioner. ";

%feature("docstring")  Ifpack_Partitioner::Compute "virtual int
Ifpack_Partitioner::Compute()=0

Computes the partitions. Returns 0 if successful. ";

%feature("docstring")  Ifpack_Partitioner::IsComputed "virtual bool
Ifpack_Partitioner::IsComputed()=0

Returns true if partitions have been computed successfully. ";

%feature("docstring")  Ifpack_Partitioner::Print "virtual ostream&
Ifpack_Partitioner::Print(std::ostream &os) const =0

Prints basic information about the partitioning object. ";


// File: classIfpack__PointRelaxation.xml
%feature("docstring") Ifpack_PointRelaxation "

Ifpack_PointRelaxation: a class to define point relaxation
preconditioners of for Epetra_RowMatrix's.

The Ifpack_PointRelaxation class enables the construction of point
relaxation preconditioners of an Epetra_RowMatrix.
Ifpack_PointRelaxation is derived from the Ifpack_Preconditioner
class, which is itself derived from Epetra_Operator. Therefore this
object can be used as preconditioner everywhere an ApplyInverse()
method is required in the preconditioning step.

This class enables the construction of the following simple
preconditioners: Jacobi;

Gauss-Seidel;

symmetric Gauss-Seidel.

We now briefly describe the main features of the above
preconditioners. Consider a linear system of type \\\\[ A x = b, \\\\]
where $A$ is a square, real matrix, and $x, b$ are two real vectors.
We begin with the decomposition \\\\[ A = D - E - F \\\\] where $D$ is
the diagonal of A, $-E$ is the strict lower part, and $-F$ is the
strict upper part. It is assumed that the diagonal entries of $A$ are
different from zero.

Given an starting solution $x_0$, an iteration of the (damped) Jacobi
method can be written in matrix form as follows: \\\\[ x_{k+1} =
\\\\omega D^{-1}(E + F) x_k + D_{-1}b, \\\\] for $k < k_{max}$, and
$\\\\omega $ a damping parameter.

Using Ifpack_Jacobi, the user can apply the specified number of sweeps
( $k_{max}$), and the damping parameter. If only one sweep is used,
then the class simply applies the inverse of the diagonal of A to the
input vector.

Given an starting solution $x_0$, an iteration of the (damped)
GaussSeidel method can be written in matrix form as follows: \\\\[ (D
- E) x_{k+1} = \\\\omega F x_k + b, \\\\] for $k < k_{max}$, and
$\\\\omega $ a damping parameter. Equivalently, the Gauss-Seidel
preconditioner can be defined as \\\\[ P_{GS}^{-1} = (D - E)^{-1}.
\\\\] Clearly, the role of E and F can be interchanged. However,
Ifpack_GaussSeidel does not consider backward Gauss-Seidel methods.

For a list of supported parameters, please refer to page ifp_params.

The complete list of supported parameters is reported in page
ifp_params. For a presentation of basic relaxation schemes, please
refer to page Ifpack_PointRelaxation.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_PointRelaxation.h ";

%feature("docstring")  Ifpack_PointRelaxation::Ifpack_PointRelaxation
"Ifpack_PointRelaxation::Ifpack_PointRelaxation(const
Epetra_RowMatrix *Matrix)

Ifpack_PointRelaxation constructor with given Epetra_RowMatrix.

Creates an instance of Ifpack_PointRelaxation class.

Parameters:
-----------

Matrix:  - (In) Pointer to matrix to precondition. ";

%feature("docstring")  Ifpack_PointRelaxation::~Ifpack_PointRelaxation
"virtual Ifpack_PointRelaxation::~Ifpack_PointRelaxation()

Destructor. ";

%feature("docstring")  Ifpack_PointRelaxation::Apply "virtual int
Ifpack_PointRelaxation::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the matrix to an Epetra_MultiVector.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  - (Out) A Epetra_MultiVector of dimension NumVectors containing
the result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Ifpack_PointRelaxation::ApplyInverse "int
Ifpack_PointRelaxation::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the preconditioner to X, returns the result in Y.

Parameters:
-----------

X:  - (In) A Epetra_MultiVector of dimension NumVectors to be
preconditioned.

Y:  - (InOut) A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

WARNING:  This routine is NOT AztecOO complaint. ";

%feature("docstring")  Ifpack_PointRelaxation::NormInf "virtual
double Ifpack_PointRelaxation::NormInf() const

Returns the infinity norm of the global matrix (not implemented) ";

%feature("docstring")  Ifpack_PointRelaxation::Label "virtual const
char* Ifpack_PointRelaxation::Label() const ";

%feature("docstring")  Ifpack_PointRelaxation::UseTranspose "virtual
bool Ifpack_PointRelaxation::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Ifpack_PointRelaxation::HasNormInf "virtual
bool Ifpack_PointRelaxation::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Ifpack_PointRelaxation::Comm "const
Epetra_Comm & Ifpack_PointRelaxation::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Ifpack_PointRelaxation::OperatorDomainMap "const Epetra_Map & Ifpack_PointRelaxation::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Ifpack_PointRelaxation::OperatorRangeMap "const Epetra_Map & Ifpack_PointRelaxation::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this
operator. ";

%feature("docstring")  Ifpack_PointRelaxation::Initialize "int
Ifpack_PointRelaxation::Initialize()

Computes all it is necessary to initialize the preconditioner. ";

%feature("docstring")  Ifpack_PointRelaxation::IsInitialized "virtual
bool Ifpack_PointRelaxation::IsInitialized() const

Returns true if the preconditioner has been successfully initialized,
false otherwise. ";

%feature("docstring")  Ifpack_PointRelaxation::IsComputed "virtual
bool Ifpack_PointRelaxation::IsComputed() const

Returns true if the preconditioner has been successfully computed. ";

%feature("docstring")  Ifpack_PointRelaxation::Compute "int
Ifpack_PointRelaxation::Compute()

Computes the preconditioners. ";

%feature("docstring")  Ifpack_PointRelaxation::Matrix "virtual const
Epetra_RowMatrix& Ifpack_PointRelaxation::Matrix() const

Returns a pointer to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_PointRelaxation::Condest "double
Ifpack_PointRelaxation::Condest(const Ifpack_CondestType
CT=Ifpack_Cheap, const int MaxIters=1550, const double Tol=1e-9,
Epetra_RowMatrix *Matrix=0)

Computes the condition number estimates and returns the value. ";

%feature("docstring")  Ifpack_PointRelaxation::Condest "virtual
double Ifpack_PointRelaxation::Condest() const

Returns the condition number estimate, or -1.0 if not computed. ";

%feature("docstring")  Ifpack_PointRelaxation::SetParameters "int
Ifpack_PointRelaxation::SetParameters(Teuchos::ParameterList &List)

Sets all the parameters for the preconditioner. ";

%feature("docstring")  Ifpack_PointRelaxation::Print "ostream &
Ifpack_PointRelaxation::Print(ostream &os) const

Prints object to an output stream. ";

%feature("docstring")  Ifpack_PointRelaxation::NumInitialize "virtual
int Ifpack_PointRelaxation::NumInitialize() const

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_PointRelaxation::NumCompute "virtual
int Ifpack_PointRelaxation::NumCompute() const

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_PointRelaxation::NumApplyInverse "virtual int Ifpack_PointRelaxation::NumApplyInverse() const

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_PointRelaxation::InitializeTime "virtual double Ifpack_PointRelaxation::InitializeTime() const

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_PointRelaxation::ComputeTime "virtual
double Ifpack_PointRelaxation::ComputeTime() const

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_PointRelaxation::ApplyInverseTime "virtual double Ifpack_PointRelaxation::ApplyInverseTime() const

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_PointRelaxation::InitializeFlops "virtual double Ifpack_PointRelaxation::InitializeFlops() const

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_PointRelaxation::ComputeFlops "virtual
double Ifpack_PointRelaxation::ComputeFlops() const

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_PointRelaxation::ApplyInverseFlops "virtual double Ifpack_PointRelaxation::ApplyInverseFlops() const

Returns the number of flops for the application of the preconditioner.
";

%feature("docstring")  Ifpack_PointRelaxation::SetUseTranspose "virtual int Ifpack_PointRelaxation::SetUseTranspose(bool
UseTranspose_in)

This flag can be used to apply the preconditioner to the transpose of
the input operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";


// File: classIfpack__Preconditioner.xml
%feature("docstring") Ifpack_Preconditioner "

Ifpack_Preconditioner: basic class for preconditioning in Ifpack.

Class Ifpack_Preconditioner is a pure virtual class, and it defines
the structure of all Ifpack preconditioners.

This class is a simple extension to Epetra_Operator. It provides the
following additional methods:  Initialize() performs all operations
based on the graph of the matrix (without considering the numerical
values);

IsInitialized() returns true if the preconditioner has been
successfully initialized;

Compute() computes all is required to apply the preconditioner, using
matrix values (and assuming that the sparsity of the matrix has not
been changed);

IsComputed() should return true if the preconditioner has been
successfully computed, false otherwise.

Condest() returns an estimation of the condition number, or -1.0 if
not available

Matrix() returns a reference to the matrix to be preconditioned.

It is required that Compute() call Initialize() if IsInitialized()
returns false. The preconditioner is applied by ApplyInverse() (which
returns if IsComputed() is false). Every time that Initialize() is
called, the object destroys all the previously allocated information,
and re-initialize the preconditioner. Every time Compute() is called,
the object re-computed the actual values of the preconditioner.

Estimating Preconditioner Condition Numbers

The condition of a matrix $B$, called $cond_p(B)$, is defined as
$cond_p(B) = \\\\|B\\\\|_p\\\\|B^{-1}\\\\|_p$ in some appropriate norm
$p$. $cond_p(B)$ gives some indication of how many accurate floating
point digits can be expected from operations involving the matrix and
its inverse. A condition number approaching the accuracy of a given
floating point number system, about 15 decimal digits in IEEE double
precision, means that any results involving $B$ or $B^{-1}$ may be
meaningless.

Method Compute() can be use to estimate of the condition number.
Compute() requires one parameter, of type Ifpack_CondestType (default
value is Ifpack_Cheap; other valid choices are Ifpack_CG and
Ifpack_GMRES).

While Ifpack_CG and Ifpack_GMRES construct and AztecOO solver, and use
methods AZ_cg_condnum and AZ_gmres_condnum to evaluate an accurate
(but very expensive) estimate of the condition number, Ifpack_Cheap
computes $\\\\|(P)^{-1}e\\\\|_\\\\infty$, which is only a very crude
estimation of the actual condition number. Note that this estimated
number can be less than 1.0. However, this approach has the following
advantages: since finding $z$ such that $P z = y$ is a basic kernel
for applying the preconditioner, computing this estimate of
$cond_\\\\infty(P^{-1})$ is performed by setting $y = e$, calling the
solve kernel to compute $z$ and then computing
$\\\\|z\\\\|_\\\\infty$;

the only cost is one application of the preconditioner.

If this estimate is very large, the application of the computed
preconditioner may generate large numerical errors. Hence, the user
may check this number, and decide to recompute the preconditioner is
the computed estimate is larger than a given threshold. This is
particularly useful in ICT and RILUK factorizations, as for ill-
conditioned matrices, we often have difficulty computing usable
incomplete factorizations. The most common source of problems is that
the factorization may encounter a small or zero pivot, in which case
the factorization can fail, or even if the factorization succeeds, the
factors may be so poorly conditioned that use of them in the iterative
phase produces meaningless results. Before we can fix this problem, we
must be able to detect it.

If IFPACK is configured with Teuchos support, method SetParameters()
should be adopted. Otherwise, users can set parameters (one
at-a-time), using methods SetParameter(), for integers and doubles.
Ifpack_Preconditioner objects overload the << operator. Derived
classes should specify a Print() method, that will be used in operator
<<.

C++ includes: Ifpack_Preconditioner.h ";

%feature("docstring")  Ifpack_Preconditioner::SetParameters "virtual
int Ifpack_Preconditioner::SetParameters(Teuchos::ParameterList
&List)=0

Sets all parameters for the preconditioner. ";

%feature("docstring")  Ifpack_Preconditioner::Initialize "virtual int
Ifpack_Preconditioner::Initialize()=0

Computes all it is necessary to initialize the preconditioner. ";

%feature("docstring")  Ifpack_Preconditioner::IsInitialized "virtual
bool Ifpack_Preconditioner::IsInitialized() const =0

Returns true if the preconditioner has been successfully initialized,
false otherwise. ";

%feature("docstring")  Ifpack_Preconditioner::Compute "virtual int
Ifpack_Preconditioner::Compute()=0

Computes all it is necessary to apply the preconditioner. ";

%feature("docstring")  Ifpack_Preconditioner::IsComputed "virtual
bool Ifpack_Preconditioner::IsComputed() const =0

Returns true if the preconditioner has been successfully computed,
false otherwise. ";

%feature("docstring")  Ifpack_Preconditioner::Condest "virtual double
Ifpack_Preconditioner::Condest(const Ifpack_CondestType
CT=Ifpack_Cheap, const int MaxIters=1550, const double Tol=1e-9,
Epetra_RowMatrix *Matrix=0)=0

Computes the condition number estimate, returns its value. ";

%feature("docstring")  Ifpack_Preconditioner::Condest "virtual double
Ifpack_Preconditioner::Condest() const =0

Returns the computed condition number estimate, or -1.0 if not
computed. ";

%feature("docstring")  Ifpack_Preconditioner::ApplyInverse "virtual
int Ifpack_Preconditioner::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const =0

Applies the preconditioner to vector X, returns the result in Y. ";

%feature("docstring")  Ifpack_Preconditioner::Matrix "virtual const
Epetra_RowMatrix& Ifpack_Preconditioner::Matrix() const =0

Returns a pointer to the matrix to be preconditioned. ";

%feature("docstring")  Ifpack_Preconditioner::NumInitialize "virtual
int Ifpack_Preconditioner::NumInitialize() const =0

Returns the number of calls to Initialize(). ";

%feature("docstring")  Ifpack_Preconditioner::NumCompute "virtual int
Ifpack_Preconditioner::NumCompute() const =0

Returns the number of calls to Compute(). ";

%feature("docstring")  Ifpack_Preconditioner::NumApplyInverse "virtual int Ifpack_Preconditioner::NumApplyInverse() const =0

Returns the number of calls to ApplyInverse(). ";

%feature("docstring")  Ifpack_Preconditioner::InitializeTime "virtual
double Ifpack_Preconditioner::InitializeTime() const =0

Returns the time spent in Initialize(). ";

%feature("docstring")  Ifpack_Preconditioner::ComputeTime "virtual
double Ifpack_Preconditioner::ComputeTime() const =0

Returns the time spent in Compute(). ";

%feature("docstring")  Ifpack_Preconditioner::ApplyInverseTime "virtual double Ifpack_Preconditioner::ApplyInverseTime() const =0

Returns the time spent in ApplyInverse(). ";

%feature("docstring")  Ifpack_Preconditioner::InitializeFlops "virtual double Ifpack_Preconditioner::InitializeFlops() const =0

Returns the number of flops in the initialization phase. ";

%feature("docstring")  Ifpack_Preconditioner::ComputeFlops "virtual
double Ifpack_Preconditioner::ComputeFlops() const =0

Returns the number of flops in the computation phase. ";

%feature("docstring")  Ifpack_Preconditioner::ApplyInverseFlops "virtual double Ifpack_Preconditioner::ApplyInverseFlops() const =0

Returns the number of flops in the application of the preconditioner.
";

%feature("docstring")  Ifpack_Preconditioner::Print "virtual ostream&
Ifpack_Preconditioner::Print(std::ostream &os) const =0

Prints basic information on iostream. This function is used by
operator<<. ";


// File: classIfpack__RCMReordering.xml
%feature("docstring") Ifpack_RCMReordering "

Ifpack_RCMReordering: reverse Cuthill-McKee reordering.

C++ includes: Ifpack_RCMReordering.h ";

%feature("docstring")  Ifpack_RCMReordering::Ifpack_RCMReordering "Ifpack_RCMReordering::Ifpack_RCMReordering()

Constructor for Ifpack_Graph's. ";

%feature("docstring")  Ifpack_RCMReordering::Ifpack_RCMReordering "Ifpack_RCMReordering::Ifpack_RCMReordering(const Ifpack_RCMReordering
&RHS)

Copy Constructor. ";

%feature("docstring")  Ifpack_RCMReordering::~Ifpack_RCMReordering "virtual Ifpack_RCMReordering::~Ifpack_RCMReordering()

Destructor. ";

%feature("docstring")  Ifpack_RCMReordering::SetParameter "int
Ifpack_RCMReordering::SetParameter(const string Name, const int Value)

Sets integer parameters `Name'. ";

%feature("docstring")  Ifpack_RCMReordering::SetParameter "int
Ifpack_RCMReordering::SetParameter(const string Name, const double
Value)

Sets double parameters `Name'. ";

%feature("docstring")  Ifpack_RCMReordering::SetParameters "int
Ifpack_RCMReordering::SetParameters(Teuchos::ParameterList &List)

Sets all parameters. ";

%feature("docstring")  Ifpack_RCMReordering::Compute "int
Ifpack_RCMReordering::Compute(const Ifpack_Graph &Graph)

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_RCMReordering::Compute "int
Ifpack_RCMReordering::Compute(const Epetra_RowMatrix &Matrix)

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_RCMReordering::IsComputed "virtual bool
Ifpack_RCMReordering::IsComputed() const

Returns true is the reordering object has been successfully
initialized, false otherwise. ";

%feature("docstring")  Ifpack_RCMReordering::Reorder "int
Ifpack_RCMReordering::Reorder(const int i) const

Returns the reordered index of row i. ";

%feature("docstring")  Ifpack_RCMReordering::InvReorder "int
Ifpack_RCMReordering::InvReorder(const int i) const

Returns the inverse reordered index of row i. ";

%feature("docstring")  Ifpack_RCMReordering::P "int
Ifpack_RCMReordering::P(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &Xreord) const

Applies reordering to multivector X, whose local length equals the
number of local rows. ";

%feature("docstring")  Ifpack_RCMReordering::Pinv "int
Ifpack_RCMReordering::Pinv(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &Xinvreord) const

Applies inverse reordering to multivector X, whose local length equals
the number of local rows. ";

%feature("docstring")  Ifpack_RCMReordering::Print "ostream &
Ifpack_RCMReordering::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";

%feature("docstring")  Ifpack_RCMReordering::NumMyRows "virtual int
Ifpack_RCMReordering::NumMyRows() const

Returns the number of local rows. ";

%feature("docstring")  Ifpack_RCMReordering::RootNode "virtual int
Ifpack_RCMReordering::RootNode() const

Returns the root node. ";


// File: classIfpack__ReorderFilter.xml
%feature("docstring") Ifpack_ReorderFilter "

Ifpack_ReorderFilter: a class for light-weight reorder of local rows
and columns of an Epetra_RowMatrix.

Class Ifpack_ReorderFilter enables a light-weight construction of
reordered matrices.

This class is used in Ifpack_AdditiveSchwarz to reorder (if required
by the user) the localized matrix. As the localized matrix is defined
on a serial communicator only, all maps are trivial (as all elements
reside on the same process). This class does not attemp to define
properly reordered maps, hence it should not be used for distributed
matrices.

To improve the performances of Ifpack_AdditiveSchwarz, some operations
are not performed in the construction phase (like for instance the
computation of the 1-norm and infinite-norm, of check whether the
reordered matrix is lower/upper triangular or not).

Marzio Sala, SNL 9214.

C++ includes: Ifpack_ReorderFilter.h ";

%feature("docstring")  Ifpack_ReorderFilter::Ifpack_ReorderFilter "Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Teuchos::RefCountPtr<
Epetra_RowMatrix > &Matrix_in, const Teuchos::RefCountPtr<
Ifpack_Reordering > &Reordering_in) ";

%feature("docstring")  Ifpack_ReorderFilter::Ifpack_ReorderFilter "Ifpack_ReorderFilter::Ifpack_ReorderFilter(const Ifpack_ReorderFilter
&RHS)

Copy constructor. ";

%feature("docstring")  Ifpack_ReorderFilter::~Ifpack_ReorderFilter "virtual Ifpack_ReorderFilter::~Ifpack_ReorderFilter()

Destructor. ";

%feature("docstring")  Ifpack_ReorderFilter::NumMyRowEntries "virtual
int Ifpack_ReorderFilter::NumMyRowEntries(int MyRow, int &NumEntries)
const

Returns the number of local row entries. ";

%feature("docstring")  Ifpack_ReorderFilter::MaxNumEntries "virtual
int Ifpack_ReorderFilter::MaxNumEntries() const

Returns maximum num entries. ";

%feature("docstring")  Ifpack_ReorderFilter::ExtractMyRowCopy "int
Ifpack_ReorderFilter::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  Ifpack_ReorderFilter::ExtractDiagonalCopy "int
Ifpack_ReorderFilter::ExtractDiagonalCopy(Epetra_Vector &Diagonal)
const

Extracts a copy of the diagonal of the reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::Multiply "int
Ifpack_ReorderFilter::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const

Multiplies multi-vector X with the reordered matrix, returns result in
Y. ";

%feature("docstring")  Ifpack_ReorderFilter::Solve "int
Ifpack_ReorderFilter::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Solve, not implemented. ";

%feature("docstring")  Ifpack_ReorderFilter::Apply "int
Ifpack_ReorderFilter::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the reordered matrix to multi-vector X, returns the result in
Y. ";

%feature("docstring")  Ifpack_ReorderFilter::ApplyInverse "virtual
int Ifpack_ReorderFilter::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Applies the inverse of this operator (not implemented). ";

%feature("docstring")  Ifpack_ReorderFilter::InvRowSums "virtual int
Ifpack_ReorderFilter::InvRowSums(Epetra_Vector &x) const

Inverse of row sums (not implemented). ";

%feature("docstring")  Ifpack_ReorderFilter::LeftScale "virtual int
Ifpack_ReorderFilter::LeftScale(const Epetra_Vector &x)

Left scale of the matrix (not implemented). ";

%feature("docstring")  Ifpack_ReorderFilter::InvColSums "virtual int
Ifpack_ReorderFilter::InvColSums(Epetra_Vector &x) const

Inverse of column sums (not implemented). ";

%feature("docstring")  Ifpack_ReorderFilter::RightScale "virtual int
Ifpack_ReorderFilter::RightScale(const Epetra_Vector &x)

Right scale of the matrix (not implemented). ";

%feature("docstring")  Ifpack_ReorderFilter::Filled "virtual bool
Ifpack_ReorderFilter::Filled() const

Returns true is the matrix called FillComplete(). ";

%feature("docstring")  Ifpack_ReorderFilter::NormInf "virtual double
Ifpack_ReorderFilter::NormInf() const

Returns the infinite-norm. ";

%feature("docstring")  Ifpack_ReorderFilter::NormOne "virtual double
Ifpack_ReorderFilter::NormOne() const

Returns the 1-norm. ";

%feature("docstring")  Ifpack_ReorderFilter::NumGlobalNonzeros "virtual int Ifpack_ReorderFilter::NumGlobalNonzeros() const

Returns the number of global nonzero elements. ";

%feature("docstring")  Ifpack_ReorderFilter::NumGlobalRows "virtual
int Ifpack_ReorderFilter::NumGlobalRows() const

Returns the number of global rows. ";

%feature("docstring")  Ifpack_ReorderFilter::NumGlobalCols "virtual
int Ifpack_ReorderFilter::NumGlobalCols() const

Returns the number of global columns. ";

%feature("docstring")  Ifpack_ReorderFilter::NumGlobalDiagonals "virtual int Ifpack_ReorderFilter::NumGlobalDiagonals() const

Returns the number of global diagonals. ";

%feature("docstring")  Ifpack_ReorderFilter::NumMyNonzeros "virtual
int Ifpack_ReorderFilter::NumMyNonzeros() const

Returns the number of local nonzero elements. ";

%feature("docstring")  Ifpack_ReorderFilter::NumMyRows "virtual int
Ifpack_ReorderFilter::NumMyRows() const

Returns the number of local rows. ";

%feature("docstring")  Ifpack_ReorderFilter::NumMyCols "virtual int
Ifpack_ReorderFilter::NumMyCols() const

Returns the number of local columns. ";

%feature("docstring")  Ifpack_ReorderFilter::NumMyDiagonals "virtual
int Ifpack_ReorderFilter::NumMyDiagonals() const

Returns the number of local diagonals. ";

%feature("docstring")  Ifpack_ReorderFilter::LowerTriangular "virtual
bool Ifpack_ReorderFilter::LowerTriangular() const

Returns true is the reordered matrix is lower triangular. ";

%feature("docstring")  Ifpack_ReorderFilter::UpperTriangular "virtual
bool Ifpack_ReorderFilter::UpperTriangular() const

Returns true is the reordered matrix is upper triangular. ";

%feature("docstring")  Ifpack_ReorderFilter::RowMatrixRowMap "virtual
const Epetra_Map& Ifpack_ReorderFilter::RowMatrixRowMap() const

Returns the row matrix of the non-reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::RowMatrixColMap "virtual
const Epetra_Map& Ifpack_ReorderFilter::RowMatrixColMap() const

Returns the column matrix of the non-reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::RowMatrixImporter "virtual const Epetra_Import* Ifpack_ReorderFilter::RowMatrixImporter()
const

Returns the importer of the non-reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::SetUseTranspose "int
Ifpack_ReorderFilter::SetUseTranspose(bool UseTranspose_in)

Sets the use of the transpose. ";

%feature("docstring")  Ifpack_ReorderFilter::UseTranspose "bool
Ifpack_ReorderFilter::UseTranspose() const

Returns true if the transpose of this matrix is used. ";

%feature("docstring")  Ifpack_ReorderFilter::HasNormInf "bool
Ifpack_ReorderFilter::HasNormInf() const

Returns true if this matrix has the infinite norm. ";

%feature("docstring")  Ifpack_ReorderFilter::Comm "const Epetra_Comm&
Ifpack_ReorderFilter::Comm() const

Returns the communicator. ";

%feature("docstring")  Ifpack_ReorderFilter::OperatorDomainMap "const
Epetra_Map& Ifpack_ReorderFilter::OperatorDomainMap() const

Returns the operator domain map of the non-reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::OperatorRangeMap "const
Epetra_Map& Ifpack_ReorderFilter::OperatorRangeMap() const

Returns the operator domain range of the non-reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::Map "const
Epetra_BlockMap& Ifpack_ReorderFilter::Map() const

Returns the map of the non-reordered matrix. ";

%feature("docstring")  Ifpack_ReorderFilter::Label "const char*
Ifpack_ReorderFilter::Label() const

Returns the label of this object. ";

%feature("docstring")  Ifpack_ReorderFilter::Matrix "Teuchos::RefCountPtr<Epetra_RowMatrix> Ifpack_ReorderFilter::Matrix()
const

Returns a reference-counted pointer to the internally stored pointer
to Epetra_RowMatrix. ";

%feature("docstring")  Ifpack_ReorderFilter::Reordering "Teuchos::RefCountPtr<Ifpack_Reordering>
Ifpack_ReorderFilter::Reordering() const

Returns a reference-counted pointer to the internally stored pointer
to Ifpack_Reordering.. ";


// File: classIfpack__Reordering.xml
%feature("docstring") Ifpack_Reordering "

Ifpack_Reordering: basic class for reordering for a Ifpack_Graph
object.

Class Ifpack_Reordering is a pure virtual class that defines the
structure of all Ifpack reordering.

The Ifpack_Graph object is used only by method Compute().

A typical code reads as follows (using for instance RCM reordering):

An Ifpack_Reordering object is a tool used by class
Ifpack_Preconditioner, to reorder the localized matrix (with or
without overlap). As its basic usage is for localized matrices, this
class takes care of reordering the local rows only. It is also
supposed that the input graph contains no singletons. This is not a
limitation, as class Ifpack_AdditiveSchwarz will filter the graph
using Ifpack_SingletonFilter before using reordering.

If IFPACK is configure with Teuchos support, method SetParameters()
should be adopted. Otherwise, users can set parameters (one
at-a-time), using methods SetParameter(), for integers and doubles.

Ifpack_Preconditioner objects overload the << operator. Derived
classes should specify a Print() method, that will be used in operator
<<.

Marzio Sala, SNL 9214.

C++ includes: Ifpack_Reordering.h ";

%feature("docstring")  Ifpack_Reordering::~Ifpack_Reordering "virtual
Ifpack_Reordering::~Ifpack_Reordering()

Destructor. ";

%feature("docstring")  Ifpack_Reordering::SetParameter "virtual int
Ifpack_Reordering::SetParameter(const string Name, const int Value)=0

Sets integer parameters `Name'. ";

%feature("docstring")  Ifpack_Reordering::SetParameter "virtual int
Ifpack_Reordering::SetParameter(const string Name, const double
Value)=0

Sets double parameters `Name'. ";

%feature("docstring")  Ifpack_Reordering::SetParameters "virtual int
Ifpack_Reordering::SetParameters(Teuchos::ParameterList &List)=0

Sets all parameters. ";

%feature("docstring")  Ifpack_Reordering::Compute "virtual int
Ifpack_Reordering::Compute(const Ifpack_Graph &Graph)=0

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_Reordering::Compute "virtual int
Ifpack_Reordering::Compute(const Epetra_RowMatrix &Matrix)=0

Computes all it is necessary to initialize the reordering object. ";

%feature("docstring")  Ifpack_Reordering::IsComputed "virtual bool
Ifpack_Reordering::IsComputed() const =0

Returns true is the reordering object has been successfully
initialized, false otherwise. ";

%feature("docstring")  Ifpack_Reordering::Reorder "virtual int
Ifpack_Reordering::Reorder(const int i) const =0

Returns the reordered index of row i. ";

%feature("docstring")  Ifpack_Reordering::InvReorder "virtual int
Ifpack_Reordering::InvReorder(const int i) const =0

Returns the inverse reordered index of row i. ";

%feature("docstring")  Ifpack_Reordering::P "virtual int
Ifpack_Reordering::P(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &X) const =0

Applies reordering to multivector Xorig, whose local length equals the
number of local rows, stores reordered vector in X. ";

%feature("docstring")  Ifpack_Reordering::Pinv "virtual int
Ifpack_Reordering::Pinv(const Epetra_MultiVector &Xorig,
Epetra_MultiVector &X) const =0

Applies inverse reordering to multivector Xorig, whose local length
equals the number of local rows, stores inverse reordered vector in X.
";

%feature("docstring")  Ifpack_Reordering::Print "virtual ostream&
Ifpack_Reordering::Print(std::ostream &os) const =0

Prints basic information on iostream. This function is used by
operator<<. ";


// File: classIfpack__SingletonFilter.xml
%feature("docstring") Ifpack_SingletonFilter "

Ifpack_SingletonFilter: Filter based on matrix entries.

C++ includes: Ifpack_SingletonFilter.h ";

%feature("docstring")  Ifpack_SingletonFilter::Ifpack_SingletonFilter
"Ifpack_SingletonFilter::Ifpack_SingletonFilter(const
Teuchos::RefCountPtr< Epetra_RowMatrix > &Matrix)

Constructor. ";

%feature("docstring")  Ifpack_SingletonFilter::~Ifpack_SingletonFilter
"virtual Ifpack_SingletonFilter::~Ifpack_SingletonFilter()

Destructor. ";

%feature("docstring")  Ifpack_SingletonFilter::NumMyRowEntries "virtual int Ifpack_SingletonFilter::NumMyRowEntries(int MyRow, int
&NumEntries) const

Returns the number of entries in MyRow. ";

%feature("docstring")  Ifpack_SingletonFilter::MaxNumEntries "virtual
int Ifpack_SingletonFilter::MaxNumEntries() const

Returns the maximum number of entries. ";

%feature("docstring")  Ifpack_SingletonFilter::ExtractMyRowCopy "int
Ifpack_SingletonFilter::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  Ifpack_SingletonFilter::ExtractDiagonalCopy "int Ifpack_SingletonFilter::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const ";

%feature("docstring")  Ifpack_SingletonFilter::Multiply "int
Ifpack_SingletonFilter::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_SingletonFilter::Solve "int
Ifpack_SingletonFilter::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  Ifpack_SingletonFilter::Apply "int
Ifpack_SingletonFilter::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_SingletonFilter::ApplyInverse "int
Ifpack_SingletonFilter::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_SingletonFilter::InvRowSums "virtual
int Ifpack_SingletonFilter::InvRowSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_SingletonFilter::LeftScale "virtual int
Ifpack_SingletonFilter::LeftScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_SingletonFilter::InvColSums "virtual
int Ifpack_SingletonFilter::InvColSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_SingletonFilter::RightScale "virtual
int Ifpack_SingletonFilter::RightScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_SingletonFilter::Filled "virtual bool
Ifpack_SingletonFilter::Filled() const ";

%feature("docstring")  Ifpack_SingletonFilter::NormInf "virtual
double Ifpack_SingletonFilter::NormInf() const ";

%feature("docstring")  Ifpack_SingletonFilter::NormOne "virtual
double Ifpack_SingletonFilter::NormOne() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumGlobalNonzeros "virtual int Ifpack_SingletonFilter::NumGlobalNonzeros() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumGlobalRows "virtual
int Ifpack_SingletonFilter::NumGlobalRows() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumGlobalCols "virtual
int Ifpack_SingletonFilter::NumGlobalCols() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumGlobalDiagonals "virtual int Ifpack_SingletonFilter::NumGlobalDiagonals() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumMyNonzeros "virtual
int Ifpack_SingletonFilter::NumMyNonzeros() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumMyRows "virtual int
Ifpack_SingletonFilter::NumMyRows() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumMyCols "virtual int
Ifpack_SingletonFilter::NumMyCols() const ";

%feature("docstring")  Ifpack_SingletonFilter::NumMyDiagonals "virtual int Ifpack_SingletonFilter::NumMyDiagonals() const ";

%feature("docstring")  Ifpack_SingletonFilter::LowerTriangular "virtual bool Ifpack_SingletonFilter::LowerTriangular() const ";

%feature("docstring")  Ifpack_SingletonFilter::UpperTriangular "virtual bool Ifpack_SingletonFilter::UpperTriangular() const ";

%feature("docstring")  Ifpack_SingletonFilter::RowMatrixRowMap "virtual const Epetra_Map& Ifpack_SingletonFilter::RowMatrixRowMap()
const ";

%feature("docstring")  Ifpack_SingletonFilter::RowMatrixColMap "virtual const Epetra_Map& Ifpack_SingletonFilter::RowMatrixColMap()
const ";

%feature("docstring")  Ifpack_SingletonFilter::RowMatrixImporter "virtual const Epetra_Import*
Ifpack_SingletonFilter::RowMatrixImporter() const ";

%feature("docstring")  Ifpack_SingletonFilter::SetUseTranspose "int
Ifpack_SingletonFilter::SetUseTranspose(bool UseTranspose_in) ";

%feature("docstring")  Ifpack_SingletonFilter::UseTranspose "bool
Ifpack_SingletonFilter::UseTranspose() const ";

%feature("docstring")  Ifpack_SingletonFilter::HasNormInf "bool
Ifpack_SingletonFilter::HasNormInf() const ";

%feature("docstring")  Ifpack_SingletonFilter::Comm "const
Epetra_Comm& Ifpack_SingletonFilter::Comm() const ";

%feature("docstring")  Ifpack_SingletonFilter::OperatorDomainMap "const Epetra_Map& Ifpack_SingletonFilter::OperatorDomainMap() const ";

%feature("docstring")  Ifpack_SingletonFilter::OperatorRangeMap "const Epetra_Map& Ifpack_SingletonFilter::OperatorRangeMap() const ";

%feature("docstring")  Ifpack_SingletonFilter::Map "const
Epetra_BlockMap& Ifpack_SingletonFilter::Map() const ";

%feature("docstring")  Ifpack_SingletonFilter::Label "const char*
Ifpack_SingletonFilter::Label() const ";

%feature("docstring")  Ifpack_SingletonFilter::SolveSingletons "int
Ifpack_SingletonFilter::SolveSingletons(const Epetra_MultiVector &RHS,
Epetra_MultiVector &LHS) ";

%feature("docstring")  Ifpack_SingletonFilter::CreateReducedRHS "int
Ifpack_SingletonFilter::CreateReducedRHS(const Epetra_MultiVector
&LHS, const Epetra_MultiVector &RHS, Epetra_MultiVector &ReducedRHS)
";

%feature("docstring")  Ifpack_SingletonFilter::UpdateLHS "int
Ifpack_SingletonFilter::UpdateLHS(const Epetra_MultiVector
&ReducedLHS, Epetra_MultiVector &LHS) ";


// File: classIfpack__SparseContainer.xml
%feature("docstring") Ifpack_SparseContainer "

Ifpack_SparseContainer: a class for storing and solving linear systems
using sparse matrices.

To understand what an IFPACK container is, please refer to the
documentation of the pure virtual class Ifpack_Container. Currently,
containers are used by class Ifpack_BlockRelaxation.

Using block methods, one needs to store all diagonal blocks and to be
also to apply the inverse of each diagonal block. Using class
Ifpack_DenseContainer, one can store the blocks as sparse matrices
(Epetra_CrsMatrix), which can be advantageous when the blocks are
large. Otherwise, class Ifpack_DenseContainer is probably more
appropriate.

Sparse containers are templated with a type T, which represent the
class to use in the application of the inverse. (T is not used in
Ifpack_DenseContainer). In SparseContainer, T must be an
Ifpack_Preconditioner derived class. The container will allocate a T
object, use SetParameters() and Compute(), then use T every time the
linear system as to be solved (using the ApplyInverse() method of T).

Marzio Sala, SNL 9214.

C++ includes: Ifpack_SparseContainer.h ";

%feature("docstring")  Ifpack_SparseContainer::Ifpack_SparseContainer
"Ifpack_SparseContainer< T >::Ifpack_SparseContainer(const int
NumRows, const int NumVectors=1)

Constructor. ";

%feature("docstring")  Ifpack_SparseContainer::Ifpack_SparseContainer
"Ifpack_SparseContainer< T >::Ifpack_SparseContainer(const
Ifpack_SparseContainer< T > &rhs)

Copy constructor. ";

%feature("docstring")  Ifpack_SparseContainer::~Ifpack_SparseContainer
"Ifpack_SparseContainer< T >::~Ifpack_SparseContainer()

Destructor. ";

%feature("docstring")  Ifpack_SparseContainer::NumRows "int
Ifpack_SparseContainer< T >::NumRows() const

Returns the number of rows of the matrix and LHS/RHS. ";

%feature("docstring")  Ifpack_SparseContainer::NumVectors "virtual
int Ifpack_SparseContainer< T >::NumVectors() const

Returns the number of vectors in LHS/RHS. ";

%feature("docstring")  Ifpack_SparseContainer::SetNumVectors "virtual
int Ifpack_SparseContainer< T >::SetNumVectors(const int
NumVectors_in)

Sets the number of vectors for LHS/RHS. ";

%feature("docstring")  Ifpack_SparseContainer::LHS "double &
Ifpack_SparseContainer< T >::LHS(const int i, const int Vector=0)

Returns the i-th component of the vector Vector of LHS. ";

%feature("docstring")  Ifpack_SparseContainer::RHS "double &
Ifpack_SparseContainer< T >::RHS(const int i, const int Vector=0)

Returns the i-th component of the vector Vector of RHS. ";

%feature("docstring")  Ifpack_SparseContainer::ID "int &
Ifpack_SparseContainer< T >::ID(const int i)

Returns the ID associated to local row i.

The set of (local) rows assigned to this container is defined by
calling ID(i) = j, where i (from 0 to NumRows()) indicates the
container-row, and j indicates the local row in the calling process.

This is usually used to recorder the local row ID (on calling process)
of the i-th row in the container. ";

%feature("docstring")  Ifpack_SparseContainer::SetMatrixElement "int
Ifpack_SparseContainer< T >::SetMatrixElement(const int row, const int
col, const double value)

Set the matrix element (row,col) to value. ";

%feature("docstring")  Ifpack_SparseContainer::IsInitialized "virtual
bool Ifpack_SparseContainer< T >::IsInitialized() const

Returns true is the container has been successfully initialized. ";

%feature("docstring")  Ifpack_SparseContainer::IsComputed "virtual
bool Ifpack_SparseContainer< T >::IsComputed() const

Returns true is the container has been successfully computed. ";

%feature("docstring")  Ifpack_SparseContainer::SetParameters "int
Ifpack_SparseContainer< T >::SetParameters(Teuchos::ParameterList
&List)

Sets all necessary parameters. ";

%feature("docstring")  Ifpack_SparseContainer::Label "virtual const
char* Ifpack_SparseContainer< T >::Label() const

Returns the label of this container. ";

%feature("docstring")  Ifpack_SparseContainer::Map "Teuchos::RCP<const Epetra_Map> Ifpack_SparseContainer< T >::Map()
const

Returns a pointer to the internally stored map. ";

%feature("docstring")  Ifpack_SparseContainer::LHS "Teuchos::RCP<const Epetra_MultiVector> Ifpack_SparseContainer< T
>::LHS() const

Returns a pointer to the internally stored solution multi-vector. ";

%feature("docstring")  Ifpack_SparseContainer::RHS "Teuchos::RCP<const Epetra_MultiVector> Ifpack_SparseContainer< T
>::RHS() const

Returns a pointer to the internally stored rhs multi-vector. ";

%feature("docstring")  Ifpack_SparseContainer::Matrix "Teuchos::RCP<const Epetra_CrsMatrix> Ifpack_SparseContainer< T
>::Matrix() const

Returns a pointer to the internally stored matrix. ";

%feature("docstring")  Ifpack_SparseContainer::ID "const
Epetra_IntSerialDenseVector& Ifpack_SparseContainer< T >::ID() const

Returns a pointer to the internally stored ID's. ";

%feature("docstring")  Ifpack_SparseContainer::Inverse "Teuchos::RCP<const T> Ifpack_SparseContainer< T >::Inverse() const

Returns a pointer to the internally stored inverse operator. ";

%feature("docstring")  Ifpack_SparseContainer::Initialize "int
Ifpack_SparseContainer< T >::Initialize()

Initializes the container, by completing all the operations based on
matrix structure.

After a call to Initialize(), no new matrix entries can be added. ";

%feature("docstring")  Ifpack_SparseContainer::Compute "int
Ifpack_SparseContainer< T >::Compute(const Epetra_RowMatrix
&Matrix_in)

Finalizes the linear system matrix and prepares for the application of
the inverse. ";

%feature("docstring")  Ifpack_SparseContainer::Apply "int
Ifpack_SparseContainer< T >::Apply()

Apply the matrix to RHS, result is stored in LHS. ";

%feature("docstring")  Ifpack_SparseContainer::ApplyInverse "int
Ifpack_SparseContainer< T >::ApplyInverse()

Apply the inverse of the matrix to RHS, result is stored in LHS. ";

%feature("docstring")  Ifpack_SparseContainer::Destroy "int
Ifpack_SparseContainer< T >::Destroy()

Destroys all data. ";

%feature("docstring")  Ifpack_SparseContainer::InitializeFlops "virtual double Ifpack_SparseContainer< T >::InitializeFlops() const

Returns the flops in Compute(). ";

%feature("docstring")  Ifpack_SparseContainer::ComputeFlops "virtual
double Ifpack_SparseContainer< T >::ComputeFlops() const

Returns the flops in Compute(). ";

%feature("docstring")  Ifpack_SparseContainer::ApplyFlops "virtual
double Ifpack_SparseContainer< T >::ApplyFlops() const

Returns the flops in Apply(). ";

%feature("docstring")  Ifpack_SparseContainer::ApplyInverseFlops "virtual double Ifpack_SparseContainer< T >::ApplyInverseFlops() const

Returns the flops in ApplyInverse(). ";

%feature("docstring")  Ifpack_SparseContainer::Print "ostream &
Ifpack_SparseContainer< T >::Print(std::ostream &os) const

Prints basic information on iostream. This function is used by
operator<<. ";


// File: classIfpack__SparsityFilter.xml
%feature("docstring") Ifpack_SparsityFilter "

Ifpack_SparsityFilter: a class to drop based on sparsity.

C++ includes: Ifpack_SparsityFilter.h ";

%feature("docstring")  Ifpack_SparsityFilter::Ifpack_SparsityFilter "Ifpack_SparsityFilter::Ifpack_SparsityFilter(const
Teuchos::RefCountPtr< Epetra_RowMatrix > &Matrix, int
AllowedNumEntries, int AllowedBandwidth=-1) ";

%feature("docstring")  Ifpack_SparsityFilter::~Ifpack_SparsityFilter "virtual Ifpack_SparsityFilter::~Ifpack_SparsityFilter() ";

%feature("docstring")  Ifpack_SparsityFilter::NumMyRowEntries "virtual int Ifpack_SparsityFilter::NumMyRowEntries(int MyRow, int
&NumEntries) const ";

%feature("docstring")  Ifpack_SparsityFilter::MaxNumEntries "virtual
int Ifpack_SparsityFilter::MaxNumEntries() const ";

%feature("docstring")  Ifpack_SparsityFilter::ExtractMyRowCopy "int
Ifpack_SparsityFilter::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  Ifpack_SparsityFilter::ExtractDiagonalCopy "int Ifpack_SparsityFilter::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const ";

%feature("docstring")  Ifpack_SparsityFilter::Multiply "int
Ifpack_SparsityFilter::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_SparsityFilter::Solve "int
Ifpack_SparsityFilter::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  Ifpack_SparsityFilter::Apply "int
Ifpack_SparsityFilter::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_SparsityFilter::ApplyInverse "int
Ifpack_SparsityFilter::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Ifpack_SparsityFilter::InvRowSums "virtual int
Ifpack_SparsityFilter::InvRowSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_SparsityFilter::LeftScale "virtual int
Ifpack_SparsityFilter::LeftScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_SparsityFilter::InvColSums "virtual int
Ifpack_SparsityFilter::InvColSums(Epetra_Vector &x) const ";

%feature("docstring")  Ifpack_SparsityFilter::RightScale "virtual int
Ifpack_SparsityFilter::RightScale(const Epetra_Vector &x) ";

%feature("docstring")  Ifpack_SparsityFilter::Filled "virtual bool
Ifpack_SparsityFilter::Filled() const ";

%feature("docstring")  Ifpack_SparsityFilter::NormInf "virtual double
Ifpack_SparsityFilter::NormInf() const ";

%feature("docstring")  Ifpack_SparsityFilter::NormOne "virtual double
Ifpack_SparsityFilter::NormOne() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumGlobalNonzeros "virtual int Ifpack_SparsityFilter::NumGlobalNonzeros() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumGlobalRows "virtual
int Ifpack_SparsityFilter::NumGlobalRows() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumGlobalCols "virtual
int Ifpack_SparsityFilter::NumGlobalCols() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumGlobalDiagonals "virtual int Ifpack_SparsityFilter::NumGlobalDiagonals() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumMyNonzeros "virtual
int Ifpack_SparsityFilter::NumMyNonzeros() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumMyRows "virtual int
Ifpack_SparsityFilter::NumMyRows() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumMyCols "virtual int
Ifpack_SparsityFilter::NumMyCols() const ";

%feature("docstring")  Ifpack_SparsityFilter::NumMyDiagonals "virtual
int Ifpack_SparsityFilter::NumMyDiagonals() const ";

%feature("docstring")  Ifpack_SparsityFilter::LowerTriangular "virtual bool Ifpack_SparsityFilter::LowerTriangular() const ";

%feature("docstring")  Ifpack_SparsityFilter::UpperTriangular "virtual bool Ifpack_SparsityFilter::UpperTriangular() const ";

%feature("docstring")  Ifpack_SparsityFilter::RowMatrixRowMap "virtual const Epetra_Map& Ifpack_SparsityFilter::RowMatrixRowMap()
const ";

%feature("docstring")  Ifpack_SparsityFilter::RowMatrixColMap "virtual const Epetra_Map& Ifpack_SparsityFilter::RowMatrixColMap()
const ";

%feature("docstring")  Ifpack_SparsityFilter::RowMatrixImporter "virtual const Epetra_Import*
Ifpack_SparsityFilter::RowMatrixImporter() const ";

%feature("docstring")  Ifpack_SparsityFilter::SetUseTranspose "int
Ifpack_SparsityFilter::SetUseTranspose(bool UseTranspose) ";

%feature("docstring")  Ifpack_SparsityFilter::UseTranspose "bool
Ifpack_SparsityFilter::UseTranspose() const ";

%feature("docstring")  Ifpack_SparsityFilter::HasNormInf "bool
Ifpack_SparsityFilter::HasNormInf() const ";

%feature("docstring")  Ifpack_SparsityFilter::Comm "const
Epetra_Comm& Ifpack_SparsityFilter::Comm() const ";

%feature("docstring")  Ifpack_SparsityFilter::OperatorDomainMap "const Epetra_Map& Ifpack_SparsityFilter::OperatorDomainMap() const ";

%feature("docstring")  Ifpack_SparsityFilter::OperatorRangeMap "const
Epetra_Map& Ifpack_SparsityFilter::OperatorRangeMap() const ";

%feature("docstring")  Ifpack_SparsityFilter::Map "const
Epetra_BlockMap& Ifpack_SparsityFilter::Map() const ";

%feature("docstring")  Ifpack_SparsityFilter::Label "const char*
Ifpack_SparsityFilter::Label() const ";


// File: classIfpack__UserPartitioner.xml
%feature("docstring") Ifpack_UserPartitioner "

Ifpack_UserPartitioner: A class to define linear partitions.

C++ includes: Ifpack_UserPartitioner.h ";

%feature("docstring")  Ifpack_UserPartitioner::Ifpack_UserPartitioner
"Ifpack_UserPartitioner::Ifpack_UserPartitioner(const Ifpack_Graph
*Graph)

Constructor. ";

%feature("docstring")  Ifpack_UserPartitioner::~Ifpack_UserPartitioner
"virtual Ifpack_UserPartitioner::~Ifpack_UserPartitioner()

Destructor. ";

%feature("docstring")  Ifpack_UserPartitioner::SetPartitionParameters
"int
Ifpack_UserPartitioner::SetPartitionParameters(Teuchos::ParameterList
&List)

Sets all the parameters for the partitioner (none for linear
partioning). ";

%feature("docstring")  Ifpack_UserPartitioner::ComputePartitions "int
Ifpack_UserPartitioner::ComputePartitions()

Computes the partitions. Returns 0 if successful. ";


// File: structrow__matrix.xml
%feature("docstring") row_matrix "";


// File: namespace@0.xml


// File: namespacestd.xml


// File: namespaceTeuchos.xml


// File: Ifpack_8cpp.xml


// File: Ifpack_8h.xml


// File: Ifpack__AdditiveSchwarz_8h.xml


// File: Ifpack__AMDReordering_8cpp.xml


// File: Ifpack__AMDReordering_8h.xml


// File: Ifpack__Amesos_8cpp.xml


// File: Ifpack__Amesos_8h.xml


// File: Ifpack__BlockRelaxation_8h.xml


// File: Ifpack__Chebyshev_8cpp.xml
%feature("docstring")  Apply_Transpose "void
Apply_Transpose(Teuchos::RCP< const Epetra_Operator > Operator_, const
Epetra_MultiVector &X, Epetra_MultiVector &Y) ";


// File: Ifpack__Chebyshev_8h.xml


// File: Ifpack__Condest_8cpp.xml
%feature("docstring")  Ifpack_Condest "double Ifpack_Condest(const
Ifpack_Preconditioner &IFP, const Ifpack_CondestType CT, const int
MaxIters, const double Tol, Epetra_RowMatrix *Matrix) ";


// File: Ifpack__Condest_8h.xml
%feature("docstring")  Ifpack_Condest "double Ifpack_Condest(const
Ifpack_Preconditioner &IFP, const Ifpack_CondestType CT, const int
MaxIters=1550, const double Tol=1e-9, Epetra_RowMatrix *Matrix=0) ";


// File: Ifpack__CondestType_8h.xml


// File: Ifpack__ConfigDefs_8h.xml


// File: Ifpack__Container_8h.xml


// File: Ifpack__CrsGraph_8h.xml


// File: Ifpack__CrsIct_8cpp.xml


// File: Ifpack__CrsIct_8h.xml


// File: Ifpack__CrsIlut_8cpp.xml


// File: Ifpack__CrsIlut_8h.xml


// File: Ifpack__CrsRick_8cpp.xml


// File: Ifpack__CrsRick_8h.xml


// File: Ifpack__CrsRiluk_8cpp.xml


// File: Ifpack__CrsRiluk_8h.xml


// File: Ifpack__DenseContainer_8cpp.xml


// File: Ifpack__DenseContainer_8h.xml


// File: Ifpack__DiagonalFilter_8cpp.xml


// File: Ifpack__DiagonalFilter_8h.xml


// File: Ifpack__DiagPreconditioner_8cpp.xml


// File: Ifpack__DiagPreconditioner_8h.xml


// File: Ifpack__DropFilter_8cpp.xml


// File: Ifpack__DropFilter_8h.xml


// File: Ifpack__EquationPartitioner_8cpp.xml


// File: Ifpack__EquationPartitioner_8h.xml


// File: Ifpack__Euclid_8cpp.xml


// File: Ifpack__Euclid_8h.xml


// File: Ifpack__Graph_8h.xml


// File: Ifpack__Graph__Epetra__CrsGraph_8cpp.xml


// File: Ifpack__Graph__Epetra__CrsGraph_8h.xml


// File: Ifpack__Graph__Epetra__RowMatrix_8cpp.xml


// File: Ifpack__Graph__Epetra__RowMatrix_8h.xml


// File: Ifpack__GreedyPartitioner_8cpp.xml


// File: Ifpack__GreedyPartitioner_8h.xml


// File: Ifpack__HashTable_8cpp.xml


// File: Ifpack__HashTable_8h.xml


// File: Ifpack__HIPS_8cpp.xml


// File: Ifpack__HIPS_8h.xml


// File: Ifpack__Hypre_8cpp.xml


// File: Ifpack__Hypre_8h.xml


// File: Ifpack__IC_8cpp.xml


// File: Ifpack__IC_8h.xml


// File: Ifpack__IC__Utils_8cpp.xml
%feature("docstring")  Ifpack_AIJMatrix_alloc "void
Ifpack_AIJMatrix_alloc(Ifpack_AIJMatrix *a, int n, int nnz) ";

%feature("docstring")  Ifpack_AIJMatrix_dealloc "void
Ifpack_AIJMatrix_dealloc(Ifpack_AIJMatrix *a) ";

%feature("docstring")  qsplit "static void qsplit(double *a, int
*ind, int n, int ncut) ";

%feature("docstring")  update_column "static void update_column(int
k, const int *ia, const int *ja, const double *a, const int *ifirst,
const int *ifirst2, const int *list2, const double *multipliers, const
double *d, int *marker, double *ta, int *itcol, int *ptalen) ";

%feature("docstring")  update_lists "static void update_lists(int k,
const int *ia, const int *ja, int *ifirst, int *list) ";

%feature("docstring")  update_lists_newcol "static void
update_lists_newcol(int k, int isk, int iptr, int *ifirst, int *list)
";

%feature("docstring")  crout_ict "void crout_ict(int n, const
Ifpack_AIJMatrix *AL, const double *Adiag, double droptol, int lfil,
Ifpack_AIJMatrix *L, double **pdiag) ";


// File: Ifpack__IC__Utils_8h.xml
%feature("docstring")  ifpack_quicksort "void ifpack_quicksort(int
*const pbase, double *const daux, int total_elems) ";

%feature("docstring")  Ifpack_AIJMatrix_dealloc "void
Ifpack_AIJMatrix_dealloc(Ifpack_AIJMatrix *a) ";

%feature("docstring")  crout_ict "void crout_ict(int n, const
Ifpack_AIJMatrix *AL, const double *Adiag, double droptol, int lfil,
Ifpack_AIJMatrix *L, double **pdiag) ";


// File: Ifpack__ICT_8cpp.xml


// File: Ifpack__ICT_8h.xml


// File: Ifpack__IHSS_8cpp.xml


// File: Ifpack__IHSS_8h.xml


// File: Ifpack__IKLU_8cpp.xml


// File: Ifpack__IKLU_8h.xml


// File: Ifpack__IKLU__Utils_8cpp.xml
%feature("docstring")  csr_spalloc "csr* csr_spalloc(int m, int n,
int nzmax, int values, int triplet) ";

%feature("docstring")  csr_sprealloc "int csr_sprealloc(csr *A, int
nzmax) ";

%feature("docstring")  csr_realloc "void* csr_realloc(void *p, int n,
size_t size, int *ok) ";

%feature("docstring")  csr_spfree "csr* csr_spfree(csr *A) ";

%feature("docstring")  csr_sfree "css* csr_sfree(css *S) ";

%feature("docstring")  csr_nfree "csrn* csr_nfree(csrn *N) ";

%feature("docstring")  csr_done "csr* csr_done(csr *C, void *w, void
*x, int ok) ";

%feature("docstring")  csr_ndone "csrn* csr_ndone(csrn *N, csr *C,
void *w, void *x, int ok) ";

%feature("docstring")  csr_idone "int* csr_idone(int *p, csr *C, void
*w, int ok) ";

%feature("docstring")  csr_cumsum "double csr_cumsum(int *p, int *c,
int n) ";

%feature("docstring")  csr_scatter "int csr_scatter(const csr *B, int
i, double alpha, int *w, double *x, int mark, csr *C, int nz) ";

%feature("docstring")  csr_add "csr* csr_add(const csr *A, const csr
*B, double alpha, double beta) ";

%feature("docstring")  csr_transpose "csr* csr_transpose(const csr
*A, int values) ";

%feature("docstring")  csr_multiply "csr* csr_multiply(const csr *A,
const csr *B) ";

%feature("docstring")  csr_sqr "css* csr_sqr(int order, const csr *A)
";

%feature("docstring")  csr_reach "int csr_reach(csr *G, const csr *B,
int k, int *xi, const int *pinv) ";

%feature("docstring")  csr_dfs "int csr_dfs(int j, csr *G, int top,
int *xi, int *pstack, const int *pinv) ";

%feature("docstring")  csr_tdfs "int csr_tdfs(int j, int k, int
*head, const int *next, int *post, int *stack) ";

%feature("docstring")  csr_lu "csrn* csr_lu(const csr *A, const css
*S, double tol) ";

%feature("docstring")  csr_spsolve "int csr_spsolve(csr *G, const csr
*B, int k, int *xi, double *x, const int *pinv, int up) ";

%feature("docstring")  csr_wclear "static int csr_wclear(int mark,
int lemax, int *w, int n) ";

%feature("docstring")  csr_diag "static int csr_diag(int i, int j,
double aij, void *other) ";

%feature("docstring")  csr_amd "int* csr_amd(int order, const csr *A)
";

%feature("docstring")  csr_print "int csr_print(const csr *A, int
brief) ";

%feature("docstring")  csr_norm "double csr_norm(const csr *A) ";

%feature("docstring")  csr_fkeep "int csr_fkeep(csr *A,
int(*fkeep)(int, int, double, void *), void *other) ";


// File: Ifpack__IKLU__Utils_8h.xml
%feature("docstring")  csr_add "csr* csr_add(const csr *A, const csr
*B, double alpha, double beta) ";

%feature("docstring")  csr_multiply "csr* csr_multiply(const csr *A,
const csr *B) ";

%feature("docstring")  csr_norm "double csr_norm(const csr *A) ";

%feature("docstring")  csr_print "int csr_print(const csr *A, int
brief) ";

%feature("docstring")  csr_transpose "csr* csr_transpose(const csr
*A, int values) ";

%feature("docstring")  csr_realloc "void* csr_realloc(void *p, int n,
size_t size, int *ok) ";

%feature("docstring")  csr_spalloc "csr* csr_spalloc(int m, int n,
int nzmax, int values, int triplet) ";

%feature("docstring")  csr_spfree "csr* csr_spfree(csr *A) ";

%feature("docstring")  csr_sprealloc "int csr_sprealloc(csr *A, int
nzmax) ";

%feature("docstring")  csr_amd "int* csr_amd(int order, const csr *A)
";

%feature("docstring")  csr_droptol "int csr_droptol(csr *A, double
tol) ";

%feature("docstring")  csr_dropzeros "int csr_dropzeros(csr *A) ";

%feature("docstring")  csr_lsolve "int csr_lsolve(const csr *L,
double *x) ";

%feature("docstring")  csr_lu "csrn* csr_lu(const csr *A, const css
*S, double tol) ";

%feature("docstring")  csr_permute "csr* csr_permute(const csr *A,
const int *pinv, const int *q, int values) ";

%feature("docstring")  csr_sqr "css* csr_sqr(int order, const csr *A)
";

%feature("docstring")  csr_usolve "int csr_usolve(const csr *U,
double *x) ";

%feature("docstring")  csr_sfree "css* csr_sfree(css *S) ";

%feature("docstring")  csr_dfree "csrd* csr_dfree(csrd *D) ";

%feature("docstring")  csr_nfree "csrn* csr_nfree(csrn *N) ";

%feature("docstring")  csr_cumsum "double csr_cumsum(int *p, int *c,
int n) ";

%feature("docstring")  csr_dfs "int csr_dfs(int j, csr *G, int top,
int *xi, int *pstack, const int *pinv) ";

%feature("docstring")  csr_reach "int csr_reach(csr *G, const csr *B,
int k, int *xi, const int *pinv) ";

%feature("docstring")  csr_scatter "int csr_scatter(const csr *A, int
j, double beta, int *w, double *x, int mark, csr *C, int nz) ";

%feature("docstring")  csr_scc "csrd* csr_scc(csr *A) ";

%feature("docstring")  csr_spsolve "int csr_spsolve(csr *G, const csr
*B, int k, int *xi, double *x, const int *pinv, int up) ";

%feature("docstring")  csr_tdfs "int csr_tdfs(int j, int k, int
*head, const int *next, int *post, int *stack) ";

%feature("docstring")  csr_dalloc "csrd* csr_dalloc(int m, int n) ";

%feature("docstring")  csr_ddone "csrd* csr_ddone(csrd *D, csr *C,
void *w, int ok) ";

%feature("docstring")  csr_done "csr* csr_done(csr *C, void *w, void
*x, int ok) ";

%feature("docstring")  csr_idone "int* csr_idone(int *p, csr *C, void
*w, int ok) ";

%feature("docstring")  csr_ndone "csrn* csr_ndone(csrn *N, csr *C,
void *w, void *x, int ok) ";

%feature("docstring")  csr_fkeep "int csr_fkeep(csr *A,
int(*fkeep)(int, int, double, void *), void *other) ";


// File: Ifpack__ILU_8cpp.xml


// File: Ifpack__ILU_8h.xml


// File: Ifpack__IlukGraph_8cpp.xml


// File: Ifpack__IlukGraph_8h.xml


// File: Ifpack__ILUT_8cpp.xml


// File: Ifpack__ILUT_8h.xml


// File: Ifpack__LinearPartitioner_8cpp.xml


// File: Ifpack__LinearPartitioner_8h.xml


// File: Ifpack__LocalFilter_8cpp.xml


// File: Ifpack__LocalFilter_8h.xml


// File: Ifpack__METISPartitioner_8cpp.xml


// File: Ifpack__METISPartitioner_8h.xml


// File: Ifpack__METISReordering_8cpp.xml


// File: Ifpack__METISReordering_8h.xml


// File: Ifpack__NodeFilter_8cpp.xml


// File: Ifpack__NodeFilter_8h.xml


// File: Ifpack__OverlapFactor_8cpp.xml


// File: Ifpack__OverlapFactorObject_8h.xml


// File: Ifpack__OverlapGraph_8cpp.xml


// File: Ifpack__OverlapGraph_8h.xml


// File: Ifpack__OverlappingPartitioner_8cpp.xml


// File: Ifpack__OverlappingPartitioner_8h.xml


// File: Ifpack__OverlappingRowMatrix_8cpp.xml


// File: Ifpack__OverlappingRowMatrix_8h.xml


// File: Ifpack__OverlapSolveObject_8cpp.xml


// File: Ifpack__OverlapSolveObject_8h.xml


// File: Ifpack__Partitioner_8h.xml


// File: Ifpack__PerturbedMatrix_8h.xml


// File: Ifpack__PointRelaxation_8cpp.xml


// File: Ifpack__PointRelaxation_8h.xml


// File: Ifpack__Preconditioner_8h.xml


// File: Ifpack__RCMReordering_8cpp.xml


// File: Ifpack__RCMReordering_8h.xml


// File: Ifpack__ReorderFilter_8cpp.xml


// File: Ifpack__ReorderFilter_8h.xml


// File: Ifpack__Reordering_8h.xml


// File: Ifpack__ScalingType_8h.xml


// File: Ifpack__SILU_8cpp.xml


// File: Ifpack__SILU_8h.xml


// File: Ifpack__SingletonFilter_8cpp.xml


// File: Ifpack__SingletonFilter_8h.xml


// File: Ifpack__SORa_8cpp.xml


// File: Ifpack__SORa_8h.xml


// File: Ifpack__SparseContainer_8h.xml


// File: Ifpack__SparsityFilter_8cpp.xml


// File: Ifpack__SparsityFilter_8h.xml


// File: Ifpack__SPARSKIT_8cpp.xml


// File: Ifpack__SPARSKIT_8h.xml


// File: Ifpack__UserPartitioner_8cpp.xml


// File: Ifpack__UserPartitioner_8h.xml


// File: Ifpack__Utils_8cpp.xml
%feature("docstring")  Ifpack_PrintLine "void Ifpack_PrintLine()

Prints a line of `=' on cout. ";

%feature("docstring")  Ifpack_BreakForDebugger "void
Ifpack_BreakForDebugger(Epetra_Comm &Comm)

Stops the execution of code, so that a debugger can be attached. ";

%feature("docstring")  Ifpack_CreateOverlappingCrsMatrix "Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(const
Epetra_RowMatrix *Matrix, const int OverlappingLevel)

Creates an overlapping Epetra_CrsMatrix. Returns 0 if OverlappingLevel
is 0. ";

%feature("docstring")  Ifpack_CreateOverlappingCrsMatrix "Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(const
Epetra_CrsGraph *Graph, const int OverlappingLevel)

Creates an overlapping Epetra_CrsGraph. Returns 0 if OverlappingLevel
is 0. ";

%feature("docstring")  Ifpack_toString "string Ifpack_toString(const
int &x)

Convertes an integer to string. ";

%feature("docstring")  Ifpack_toString "string Ifpack_toString(const
double &x)

Convertes a double to string. ";

%feature("docstring")  Ifpack_PrintResidual "int
Ifpack_PrintResidual(char *Label, const Epetra_RowMatrix &A, const
Epetra_MultiVector &X, const Epetra_MultiVector &Y)

Prints on cout the true residual. ";

%feature("docstring")  Ifpack_PrintResidual "int
Ifpack_PrintResidual(const int iter, const Epetra_RowMatrix &A, const
Epetra_MultiVector &X, const Epetra_MultiVector &Y) ";

%feature("docstring")  Ifpack_PrintSparsity_Simple "void
Ifpack_PrintSparsity_Simple(const Epetra_RowMatrix &A) ";

%feature("docstring")  Ifpack_FrobeniusNorm "double
Ifpack_FrobeniusNorm(const Epetra_RowMatrix &A) ";

%feature("docstring")  print "static void print() ";

%feature("docstring")  print "static void print(const char str[], T
val) ";

%feature("docstring")  print "static void print(const char str[], T
val, double percentage) ";

%feature("docstring")  print "static void print(const char str[], T
one, T two, T three, bool equal=true) ";

%feature("docstring")  Ifpack_Analyze "int Ifpack_Analyze(const
Epetra_RowMatrix &A, const bool Cheap, const int NumPDEEqns)

Analyzes the basic properties of the input matrix A; see ifp_analyze.
";

%feature("docstring")  Ifpack_AnalyzeVectorElements "int
Ifpack_AnalyzeVectorElements(const Epetra_Vector &Diagonal, const bool
abs, const int steps)

Analyzes the distribution of values of the input vector Diagonal.

Parameters:
-----------

Diagonal:  - (In) Vector to be analyzed.

abs:  - (In) if true, the function will analyze vector B, whose
elements are defined as $ B_{i} = | D_{i}| $.

steps:  - (In) number of intervals for the analysis.

An example of output is reported ifp_vector. ";

%feature("docstring")  Ifpack_AnalyzeMatrixElements "int
Ifpack_AnalyzeMatrixElements(const Epetra_RowMatrix &A, const bool
abs, const int steps)

Analyzes the distribution of values of the input matrix A.

Parameters:
-----------

A:  - (In) matrix to be analyzed.

abs:  - (In) if true, the function will analyze matrix B, whose
elements are defined as $ B_{i,i} = | A_{i,i}| $.

steps:  - (In) number of intervals for the analysis.

An example of output is reported ifp_matrix. ";

%feature("docstring")  Ifpack_PrintSparsity "int
Ifpack_PrintSparsity(const Epetra_RowMatrix &A, const char
*InputFileName, const int NumPDEEqns) ";


// File: Ifpack__Utils_8h.xml
/*  Largely inspired from Yousef Saad's SPARSKIT plot function.  */

/* Plots the sparsity pattern of an Epetra_RowMatrix into a PS file.

Parameters:
-----------

A:  (In) - Epetra_RowMatrix whose sparsity pattern will be plotted.

FileName:  (In) - char string containing the filename. If 0, then the
matrix label is used as file name, after appending .ps.

NumPDEEqns:  (In) - number of PDE equations. The function will plot
the block structure of the matrix if NumPDEEqns > 1

*/

%feature("docstring")  Ifpack_PrintSparsity "int
Ifpack_PrintSparsity(const Epetra_RowMatrix &A, const char
*FileName=0, const int NumPDEEqns=1) ";

%feature("docstring")  Ifpack_PrintLine "void Ifpack_PrintLine()

Prints a line of `=' on cout. ";

%feature("docstring")  Ifpack_BreakForDebugger "void
Ifpack_BreakForDebugger(Epetra_Comm &Comm)

Stops the execution of code, so that a debugger can be attached. ";

%feature("docstring")  Ifpack_CreateOverlappingCrsMatrix "Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(const
Epetra_RowMatrix *Matrix, const int OverlappingLevel)

Creates an overlapping Epetra_CrsMatrix. Returns 0 if OverlappingLevel
is 0. ";

%feature("docstring")  Ifpack_CreateOverlappingCrsMatrix "Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(const
Epetra_CrsGraph *Graph, const int OverlappingLevel)

Creates an overlapping Epetra_CrsGraph. Returns 0 if OverlappingLevel
is 0. ";

%feature("docstring")  Ifpack_toString "string Ifpack_toString(const
int &x)

Convertes an integer to string. ";

%feature("docstring")  Ifpack_toString "string Ifpack_toString(const
double &x)

Convertes a double to string. ";

%feature("docstring")  Ifpack_PrintResidual "int
Ifpack_PrintResidual(char *Label, const Epetra_RowMatrix &A, const
Epetra_MultiVector &X, const Epetra_MultiVector &Y)

Prints on cout the true residual. ";

%feature("docstring")  Ifpack_PrintResidual "int
Ifpack_PrintResidual(const int iter, const Epetra_RowMatrix &A, const
Epetra_MultiVector &X, const Epetra_MultiVector &Y) ";

%feature("docstring")  Ifpack_PrintSparsity_Simple "void
Ifpack_PrintSparsity_Simple(const Epetra_RowMatrix &A) ";

%feature("docstring")  Ifpack_Analyze "int Ifpack_Analyze(const
Epetra_RowMatrix &A, const bool Cheap=false, const int NumPDEEqns=1)

Analyzes the basic properties of the input matrix A; see ifp_analyze.
";

%feature("docstring")  Ifpack_AnalyzeMatrixElements "int
Ifpack_AnalyzeMatrixElements(const Epetra_RowMatrix &A, const bool
abs=false, const int steps=10)

Analyzes the distribution of values of the input matrix A.

Parameters:
-----------

A:  - (In) matrix to be analyzed.

abs:  - (In) if true, the function will analyze matrix B, whose
elements are defined as $ B_{i,i} = | A_{i,i}| $.

steps:  - (In) number of intervals for the analysis.

An example of output is reported ifp_matrix. ";

%feature("docstring")  Ifpack_AnalyzeVectorElements "int
Ifpack_AnalyzeVectorElements(const Epetra_Vector &Diagonal, const bool
abs=false, const int steps=10)

Analyzes the distribution of values of the input vector Diagonal.

Parameters:
-----------

Diagonal:  - (In) Vector to be analyzed.

abs:  - (In) if true, the function will analyze vector B, whose
elements are defined as $ B_{i} = | D_{i}| $.

steps:  - (In) number of intervals for the analysis.

An example of output is reported ifp_vector. ";


// File: Ifpack__ValidParameters_8cpp.xml
%feature("docstring")  Ifpack_GetValidParameters "Teuchos::ParameterList Ifpack_GetValidParameters()

Returns a list which contains all the parameters possibly used by
IFPACK. ";


// File: Ifpack__ValidParameters_8h.xml
%feature("docstring")  Ifpack_GetValidParameters "Teuchos::ParameterList Ifpack_GetValidParameters()

Returns a list which contains all the parameters possibly used by
IFPACK. ";


// File: Ifpack__Version_8h.xml
%feature("docstring")  Ifpack_Version "string Ifpack_Version() ";


// File: dir_9b150e231a3088b15fbe04166109c688.xml


// File: dir_0edeec2379d245cfc73672c515d4a02d.xml


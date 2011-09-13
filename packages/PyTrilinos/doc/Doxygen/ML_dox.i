
// File: index.xml

// File: classMLAPI_1_1BaseLinearCombination.xml
%feature("docstring") MLAPI::BaseLinearCombination "";

%feature("docstring")
MLAPI::BaseLinearCombination::~BaseLinearCombination "virtual
MLAPI::BaseLinearCombination::~BaseLinearCombination() ";

%feature("docstring")  MLAPI::BaseLinearCombination::GetVectorSpace "virtual const Space MLAPI::BaseLinearCombination::GetVectorSpace()
const =0

Returns the vector space of the underlying object. ";

%feature("docstring")  MLAPI::BaseLinearCombination::Update "virtual
void MLAPI::BaseLinearCombination::Update(MultiVector &v) const =0 ";

%feature("docstring")  MLAPI::BaseLinearCombination::Set "virtual
void MLAPI::BaseLinearCombination::Set(MultiVector &v) const =0 ";


// File: classMLAPI_1_1BaseObject.xml
%feature("docstring") MLAPI::BaseObject "

Basic class for MLAPI objects.

BaseObject is the basic class for all MLAPI objects. Currently, it
contains the label of the object and method Print().

Marzio Sala, SNL 9214

C++ includes: MLAPI_BaseObject.h ";

%feature("docstring")  MLAPI::BaseObject::BaseObject "MLAPI::BaseObject::BaseObject()

Constructor with empty label. ";

%feature("docstring")  MLAPI::BaseObject::BaseObject "MLAPI::BaseObject::BaseObject(const string &Label)

Constructor with given Label. ";

%feature("docstring")  MLAPI::BaseObject::~BaseObject "virtual
MLAPI::BaseObject::~BaseObject()

Destructor. ";

%feature("docstring")  MLAPI::BaseObject::SetLabel "void
MLAPI::BaseObject::SetLabel(const string &Label)

Sets the Label of this object to Label. ";

%feature("docstring")  MLAPI::BaseObject::GetLabel "const string&
MLAPI::BaseObject::GetLabel() const

Returns the Label of this object. ";

%feature("docstring")  MLAPI::BaseObject::Print "virtual
std::ostream& MLAPI::BaseObject::Print(std::ostream &os, const bool
Verbose=true) const =0

Prints information on stream. ";


// File: classMLAPI_1_1BaseOperator.xml
%feature("docstring") MLAPI::BaseOperator "

Base class for all MLAPI objects.

Marzio Sala, SNL 9214.

C++ includes: MLAPI_BaseOperator.h ";

%feature("docstring")  MLAPI::BaseOperator::~BaseOperator "virtual
MLAPI::BaseOperator::~BaseOperator()

Virtual destructor. ";

%feature("docstring")  MLAPI::BaseOperator::Apply "virtual int
MLAPI::BaseOperator::Apply(const MultiVector &LHS, MultiVector &RHS)
const =0

Applies the operator to X, using Y as starting solution. Returns the
solution in Y. ";

%feature("docstring")  MLAPI::BaseOperator::GetOperatorDomainSpace "virtual const Space MLAPI::BaseOperator::GetOperatorDomainSpace()
const =0

Returns a copy of the domain space of this object. ";

%feature("docstring")  MLAPI::BaseOperator::GetOperatorRangeSpace "virtual const Space MLAPI::BaseOperator::GetOperatorRangeSpace() const
=0

Returns a copy of the range space of this object. ";


// File: classMLAPI_1_1BaseOperatorTimesMultiVector.xml
%feature("docstring") MLAPI::BaseOperatorTimesMultiVector "";

%feature("docstring")
MLAPI::BaseOperatorTimesMultiVector::BaseOperatorTimesMultiVector "MLAPI::BaseOperatorTimesMultiVector::BaseOperatorTimesMultiVector(const
BaseOperator &A, const MultiVector &x) ";

%feature("docstring")
MLAPI::BaseOperatorTimesMultiVector::GetVectorSpace "const Space
MLAPI::BaseOperatorTimesMultiVector::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")
MLAPI::BaseOperatorTimesMultiVector::GetBaseOperator "const
BaseOperator& MLAPI::BaseOperatorTimesMultiVector::GetBaseOperator()
const ";

%feature("docstring")
MLAPI::BaseOperatorTimesMultiVector::GetMultiVector "const
MultiVector& MLAPI::BaseOperatorTimesMultiVector::GetMultiVector()
const ";

%feature("docstring")  MLAPI::BaseOperatorTimesMultiVector::Update "void MLAPI::BaseOperatorTimesMultiVector::Update(MultiVector &v) const
";

%feature("docstring")  MLAPI::BaseOperatorTimesMultiVector::Set "void
MLAPI::BaseOperatorTimesMultiVector::Set(MultiVector &v) const ";


// File: classMLAPI_1_1CompObject.xml
%feature("docstring") MLAPI::CompObject "

Class to count flops.

Marzio Sala, SNL 9214

C++ includes: MLAPI_CompObject.h ";

%feature("docstring")  MLAPI::CompObject::CompObject "MLAPI::CompObject::CompObject()

Constructor, set counter to 0.0. ";

%feature("docstring")  MLAPI::CompObject::~CompObject "MLAPI::CompObject::~CompObject()

Destructor. ";

%feature("docstring")  MLAPI::CompObject::GetFlops "double
MLAPI::CompObject::GetFlops() const

Returns the internal counter of flops. ";

%feature("docstring")  MLAPI::CompObject::SetFlops "void
MLAPI::CompObject::SetFlops(double Flops) const

Sets internal counter to Flops. ";

%feature("docstring")  MLAPI::CompObject::UpdateFlops "void
MLAPI::CompObject::UpdateFlops(double Flops) const

Updates internal counter by summing Flops. ";


// File: classMLAPI_1_1DistributedMatrix.xml
%feature("docstring") MLAPI::DistributedMatrix "";

%feature("docstring")  MLAPI::DistributedMatrix::DistributedMatrix "MLAPI::DistributedMatrix::DistributedMatrix(const Space &RowSpace,
const Space &ColSpace) ";

%feature("docstring")  MLAPI::DistributedMatrix::NumMyRowEntries "virtual int MLAPI::DistributedMatrix::NumMyRowEntries(int MyRow, int
&NumEntries) const ";

%feature("docstring")  MLAPI::DistributedMatrix::MaxNumEntries "virtual int MLAPI::DistributedMatrix::MaxNumEntries() const ";

%feature("docstring")  MLAPI::DistributedMatrix::ExtractMyRowCopy "virtual int MLAPI::DistributedMatrix::ExtractMyRowCopy(int MyRow, int
Length, int &NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  MLAPI::DistributedMatrix::ExtractDiagonalCopy "virtual int
MLAPI::DistributedMatrix::ExtractDiagonalCopy(Epetra_Vector &Diagonal)
const ";

%feature("docstring")  MLAPI::DistributedMatrix::Multiply "virtual
int MLAPI::DistributedMatrix::Multiply(bool TransA, const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  MLAPI::DistributedMatrix::Solve "virtual int
MLAPI::DistributedMatrix::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  MLAPI::DistributedMatrix::InvRowSums "virtual
int MLAPI::DistributedMatrix::InvRowSums(Epetra_Vector &x) const ";

%feature("docstring")  MLAPI::DistributedMatrix::LeftScale "virtual
int MLAPI::DistributedMatrix::LeftScale(const Epetra_Vector &x) ";

%feature("docstring")  MLAPI::DistributedMatrix::InvColSums "virtual
int MLAPI::DistributedMatrix::InvColSums(Epetra_Vector &x) const ";

%feature("docstring")  MLAPI::DistributedMatrix::RightScale "virtual
int MLAPI::DistributedMatrix::RightScale(const Epetra_Vector &x) ";

%feature("docstring")  MLAPI::DistributedMatrix::Filled "virtual bool
MLAPI::DistributedMatrix::Filled() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NormInf "virtual
double MLAPI::DistributedMatrix::NormInf() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NormOne "virtual
double MLAPI::DistributedMatrix::NormOne() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumGlobalNonzeros "virtual int MLAPI::DistributedMatrix::NumGlobalNonzeros() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumGlobalRows "virtual int MLAPI::DistributedMatrix::NumGlobalRows() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumGlobalCols "virtual int MLAPI::DistributedMatrix::NumGlobalCols() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumGlobalDiagonals "virtual int MLAPI::DistributedMatrix::NumGlobalDiagonals() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumMyNonzeros "virtual int MLAPI::DistributedMatrix::NumMyNonzeros() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumMyRows "virtual
int MLAPI::DistributedMatrix::NumMyRows() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumMyCols "virtual
int MLAPI::DistributedMatrix::NumMyCols() const ";

%feature("docstring")  MLAPI::DistributedMatrix::NumMyDiagonals "virtual int MLAPI::DistributedMatrix::NumMyDiagonals() const ";

%feature("docstring")  MLAPI::DistributedMatrix::LowerTriangular "virtual bool MLAPI::DistributedMatrix::LowerTriangular() const ";

%feature("docstring")  MLAPI::DistributedMatrix::UpperTriangular "virtual bool MLAPI::DistributedMatrix::UpperTriangular() const ";

%feature("docstring")  MLAPI::DistributedMatrix::RowMatrixRowMap "virtual const Epetra_Map& MLAPI::DistributedMatrix::RowMatrixRowMap()
const ";

%feature("docstring")  MLAPI::DistributedMatrix::RowMatrixColMap "virtual const Epetra_Map& MLAPI::DistributedMatrix::RowMatrixColMap()
const ";

%feature("docstring")  MLAPI::DistributedMatrix::RowMatrixImporter "virtual const Epetra_Import*
MLAPI::DistributedMatrix::RowMatrixImporter() const ";

%feature("docstring")  MLAPI::DistributedMatrix::OperatorDomainMap "virtual const Epetra_Map&
MLAPI::DistributedMatrix::OperatorDomainMap() const ";

%feature("docstring")  MLAPI::DistributedMatrix::OperatorRangeMap "virtual const Epetra_Map& MLAPI::DistributedMatrix::OperatorRangeMap()
const ";

%feature("docstring")  MLAPI::DistributedMatrix::Map "virtual const
Epetra_Map& MLAPI::DistributedMatrix::Map() const ";

%feature("docstring")  MLAPI::DistributedMatrix::SetUseTranspose "virtual int MLAPI::DistributedMatrix::SetUseTranspose(bool what) ";

%feature("docstring")  MLAPI::DistributedMatrix::Apply "int
MLAPI::DistributedMatrix::Apply(const MultiVector &X, MultiVector &Y)
const

Applies this operator to LHS, returns the result in RHS. ";

%feature("docstring")  MLAPI::DistributedMatrix::Apply "virtual int
MLAPI::DistributedMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  MLAPI::DistributedMatrix::ApplyInverse "virtual int MLAPI::DistributedMatrix::ApplyInverse(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  MLAPI::DistributedMatrix::Label "virtual const
char* MLAPI::DistributedMatrix::Label() const ";

%feature("docstring")  MLAPI::DistributedMatrix::UseTranspose "virtual bool MLAPI::DistributedMatrix::UseTranspose() const ";

%feature("docstring")  MLAPI::DistributedMatrix::HasNormInf "virtual
bool MLAPI::DistributedMatrix::HasNormInf() const ";

%feature("docstring")  MLAPI::DistributedMatrix::Comm "virtual const
Epetra_Comm& MLAPI::DistributedMatrix::Comm() const ";

%feature("docstring")  MLAPI::DistributedMatrix::Print "std::ostream&
MLAPI::DistributedMatrix::Print(std::ostream &os, const bool
verbose=true) const

Prints basic information about this object. ";

%feature("docstring")  MLAPI::DistributedMatrix::GetDomainSpace "Space MLAPI::DistributedMatrix::GetDomainSpace() const

Returns a reference to the internally stored domain space. ";

%feature("docstring")  MLAPI::DistributedMatrix::GetRangeSpace "Space
MLAPI::DistributedMatrix::GetRangeSpace() const

Returns a reference to the internally stored range space. ";

%feature("docstring")  MLAPI::DistributedMatrix::ReplaceElement "void
MLAPI::DistributedMatrix::ReplaceElement(const int GRID, const int
GCID, const double value) ";

%feature("docstring")  MLAPI::DistributedMatrix::FillComplete "void
MLAPI::DistributedMatrix::FillComplete() ";

%feature("docstring")  MLAPI::DistributedMatrix::IsFillCompleted "bool MLAPI::DistributedMatrix::IsFillCompleted() const ";


// File: classMLAPI_1_1DoubleVector.xml
%feature("docstring") MLAPI::DoubleVector "";

%feature("docstring")  MLAPI::DoubleVector::DoubleVector "MLAPI::DoubleVector::DoubleVector(const int size) ";

%feature("docstring")  MLAPI::DoubleVector::DoubleVector "MLAPI::DoubleVector::DoubleVector(double *ptr) ";

%feature("docstring")  MLAPI::DoubleVector::~DoubleVector "MLAPI::DoubleVector::~DoubleVector() ";

%feature("docstring")  MLAPI::DoubleVector::Values "double*
MLAPI::DoubleVector::Values() ";

%feature("docstring")  MLAPI::DoubleVector::Values "const double*
MLAPI::DoubleVector::Values() const ";


// File: classMLAPI_1_1Epetra__SerialMatrix.xml
%feature("docstring") MLAPI::Epetra_SerialMatrix "";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Epetra_SerialMatrix
"MLAPI::Epetra_SerialMatrix::Epetra_SerialMatrix(const Space
&RowSpace, const Space &ColSpace) ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumMyRowEntries "virtual int MLAPI::Epetra_SerialMatrix::NumMyRowEntries(int MyRow, int
&NumEntries) const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::MaxNumEntries "virtual int MLAPI::Epetra_SerialMatrix::MaxNumEntries() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::ExtractMyRowCopy "virtual int MLAPI::Epetra_SerialMatrix::ExtractMyRowCopy(int MyRow,
int Length, int &NumEntries, double *Values, int *Indices) const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::ExtractDiagonalCopy
"virtual int
MLAPI::Epetra_SerialMatrix::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Multiply "virtual
int MLAPI::Epetra_SerialMatrix::Multiply(bool TransA, const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Solve "virtual int
MLAPI::Epetra_SerialMatrix::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::InvRowSums "virtual int MLAPI::Epetra_SerialMatrix::InvRowSums(Epetra_Vector &x)
const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::LeftScale "virtual
int MLAPI::Epetra_SerialMatrix::LeftScale(const Epetra_Vector &x) ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::InvColSums "virtual int MLAPI::Epetra_SerialMatrix::InvColSums(Epetra_Vector &x)
const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::RightScale "virtual int MLAPI::Epetra_SerialMatrix::RightScale(const Epetra_Vector
&x) ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Filled "virtual
bool MLAPI::Epetra_SerialMatrix::Filled() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NormInf "virtual
double MLAPI::Epetra_SerialMatrix::NormInf() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NormOne "virtual
double MLAPI::Epetra_SerialMatrix::NormOne() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumGlobalNonzeros "virtual int MLAPI::Epetra_SerialMatrix::NumGlobalNonzeros() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumGlobalRows "virtual int MLAPI::Epetra_SerialMatrix::NumGlobalRows() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumGlobalCols "virtual int MLAPI::Epetra_SerialMatrix::NumGlobalCols() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumGlobalDiagonals
"virtual int MLAPI::Epetra_SerialMatrix::NumGlobalDiagonals() const
";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumMyNonzeros "virtual int MLAPI::Epetra_SerialMatrix::NumMyNonzeros() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumMyRows "virtual
int MLAPI::Epetra_SerialMatrix::NumMyRows() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumMyCols "virtual
int MLAPI::Epetra_SerialMatrix::NumMyCols() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::NumMyDiagonals "virtual int MLAPI::Epetra_SerialMatrix::NumMyDiagonals() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::LowerTriangular "virtual bool MLAPI::Epetra_SerialMatrix::LowerTriangular() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::UpperTriangular "virtual bool MLAPI::Epetra_SerialMatrix::UpperTriangular() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::RowMatrixRowMap "virtual const Epetra_Map&
MLAPI::Epetra_SerialMatrix::RowMatrixRowMap() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::RowMatrixColMap "virtual const Epetra_Map&
MLAPI::Epetra_SerialMatrix::RowMatrixColMap() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::RowMatrixImporter "virtual const Epetra_Import*
MLAPI::Epetra_SerialMatrix::RowMatrixImporter() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::OperatorDomainMap "virtual const Epetra_Map&
MLAPI::Epetra_SerialMatrix::OperatorDomainMap() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::OperatorRangeMap "virtual const Epetra_Map&
MLAPI::Epetra_SerialMatrix::OperatorRangeMap() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Map "virtual const
Epetra_Map& MLAPI::Epetra_SerialMatrix::Map() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::SetUseTranspose "virtual int MLAPI::Epetra_SerialMatrix::SetUseTranspose(bool) ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Apply "virtual int
MLAPI::Epetra_SerialMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::ApplyInverse "virtual int MLAPI::Epetra_SerialMatrix::ApplyInverse(const
Epetra_MultiVector &X, Epetra_MultiVector &Y) const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Label "virtual
const char* MLAPI::Epetra_SerialMatrix::Label() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::UseTranspose "virtual bool MLAPI::Epetra_SerialMatrix::UseTranspose() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::HasNormInf "virtual bool MLAPI::Epetra_SerialMatrix::HasNormInf() const ";

%feature("docstring")  MLAPI::Epetra_SerialMatrix::Comm "virtual
const Epetra_Comm& MLAPI::Epetra_SerialMatrix::Comm() const ";


// File: classMLAPI_1_1EpetraBaseOperator.xml
%feature("docstring") MLAPI::EpetraBaseOperator "

Basic class to wrap MLAPI::InverseOperator into Epetra_Operator.

Marzio Sala, SNL 9214.

C++ includes: MLAPI_EpetraBaseOperator.h ";

%feature("docstring")  MLAPI::EpetraBaseOperator::EpetraBaseOperator "MLAPI::EpetraBaseOperator::EpetraBaseOperator(const Epetra_Map &Map,
const BaseOperator &Op)

Constructor. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::~EpetraBaseOperator
"virtual MLAPI::EpetraBaseOperator::~EpetraBaseOperator()

Destructor. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::ApplyInverse "int
MLAPI::EpetraBaseOperator::ApplyInverse(const Epetra_MultiVector
&X_Epetra, Epetra_MultiVector &Y_Epetra) const

Applies the operator to X, returns the results in Y.

Apply() and ApplyInverse() are the SAME function! ";

%feature("docstring")  MLAPI::EpetraBaseOperator::SetUseTranspose "virtual int MLAPI::EpetraBaseOperator::SetUseTranspose(bool
UseTranspose)

Sets the use of tranpose (NOT IMPLEMENTED). ";

%feature("docstring")  MLAPI::EpetraBaseOperator::Apply "virtual int
MLAPI::EpetraBaseOperator::Apply(const Epetra_MultiVector &X_Epetra,
Epetra_MultiVector &Y_Epetra) const

Applies the operator to X, returns the results in Y. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::NormInf "virtual
double MLAPI::EpetraBaseOperator::NormInf() const

NOT IMPLEMENTED. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::Label "virtual
const char* MLAPI::EpetraBaseOperator::Label() const

Returns the label of this object. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::UseTranspose "virtual bool MLAPI::EpetraBaseOperator::UseTranspose() const

Returns false. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::HasNormInf "virtual
bool MLAPI::EpetraBaseOperator::HasNormInf() const

NOT IMPLEMENTED. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::Comm "virtual const
Epetra_Comm& MLAPI::EpetraBaseOperator::Comm() const

Returns a reference to the communicator object. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::OperatorDomainMap "virtual const Epetra_Map&
MLAPI::EpetraBaseOperator::OperatorDomainMap() const

Returns a reference to the OperatorDomainMap. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::OperatorRangeMap "virtual const Epetra_Map&
MLAPI::EpetraBaseOperator::OperatorRangeMap() const

Returns a reference to the OperatorRangeMap. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::Map "virtual const
Epetra_Map& MLAPI::EpetraBaseOperator::Map() const

Returns a reference to the Map of this object. ";

%feature("docstring")  MLAPI::EpetraBaseOperator::GetBaseOperator "const BaseOperator& MLAPI::EpetraBaseOperator::GetBaseOperator() const
";


// File: classMLAPI_1_1InverseOperator.xml
%feature("docstring") MLAPI::InverseOperator "

InverseOperator: basic class to define smoother and coarse solvers.

Marzio Sala, D-INFK/ETHZ.

C++ includes: MLAPI_LoadBalanceInverseOperator.h ";

%feature("docstring")  MLAPI::InverseOperator::InverseOperator "MLAPI::InverseOperator::InverseOperator()

Empty constructor. ";

%feature("docstring")  MLAPI::InverseOperator::InverseOperator "MLAPI::InverseOperator::InverseOperator(const Operator &Op, const
string Type)

Constructor for a given Operator and type, and default parameters. ";

%feature("docstring")  MLAPI::InverseOperator::InverseOperator "MLAPI::InverseOperator::InverseOperator(const Operator &Op, const
string Type, Teuchos::ParameterList &List)

Constructor for a given Operator, type and parameters. ";

%feature("docstring")  MLAPI::InverseOperator::InverseOperator "MLAPI::InverseOperator::InverseOperator(const InverseOperator &RHS)

Copy constructor. ";

%feature("docstring")  MLAPI::InverseOperator::~InverseOperator "MLAPI::InverseOperator::~InverseOperator()

Destructor. ";

%feature("docstring")  MLAPI::InverseOperator::Reshape "void
MLAPI::InverseOperator::Reshape()

Resets this object. ";

%feature("docstring")  MLAPI::InverseOperator::Reshape "void
MLAPI::InverseOperator::Reshape(const Operator &Op, const string Type)

Reshapes the object with default values. ";

%feature("docstring")  MLAPI::InverseOperator::Reshape "void
MLAPI::InverseOperator::Reshape(const Operator &Op, const string Type,
Teuchos::ParameterList &List, Teuchos::ParameterList *pushlist=NULL)

Reshapes the object by setting the Operator and the specified type. ";

%feature("docstring")  MLAPI::InverseOperator::Reshape "void
MLAPI::InverseOperator::Reshape(Ifpack_Preconditioner *prec, const
Operator &Op, const bool ownership)

Reshape with preconstructed smoother as Ifpack_Preconditioner. ";

%feature("docstring")  MLAPI::InverseOperator::GetOperatorRangeSpace "const Space MLAPI::InverseOperator::GetOperatorRangeSpace() const

Returns a reference to the range space of this object. ";

%feature("docstring")  MLAPI::InverseOperator::GetOperatorDomainSpace
"const Space MLAPI::InverseOperator::GetOperatorDomainSpace() const

Returns a reference to the domain space of this object. ";

%feature("docstring")  MLAPI::InverseOperator::GetRangeSpace "const
Space MLAPI::InverseOperator::GetRangeSpace() const

Returns a reference to the range space of this object. ";

%feature("docstring")  MLAPI::InverseOperator::GetDomainSpace "const
Space MLAPI::InverseOperator::GetDomainSpace() const

Returns a reference to the domain space of this object. ";

%feature("docstring")  MLAPI::InverseOperator::RCPRowMatrix "const
Teuchos::RefCountPtr<Epetra_RowMatrix>
MLAPI::InverseOperator::RCPRowMatrix() const

Returns pointer of the internally stored ML_Epetra::RowMatrix object.
";

%feature("docstring")  MLAPI::InverseOperator::RowMatrix "Epetra_RowMatrix* MLAPI::InverseOperator::RowMatrix() const

Returns pointer of the internally stored ML_Epetra::RowMatrix object.
";

%feature("docstring")  MLAPI::InverseOperator::GetOperator "const
Operator& MLAPI::InverseOperator::GetOperator() const

Returns a reference to the Operator of which this object defines the
inverse. ";

%feature("docstring")  MLAPI::InverseOperator::GetRCPData "Teuchos::RefCountPtr<Ifpack_Preconditioner>&
MLAPI::InverseOperator::GetRCPData()

Returns a pointer to the internally stored IFPACK preconditioner. ";

%feature("docstring")  MLAPI::InverseOperator::GetRCPMLPrec "Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>&
MLAPI::InverseOperator::GetRCPMLPrec()

Returns a pointer to the internally stored IFPACK preconditioner. ";

%feature("docstring")  MLAPI::InverseOperator::GetRCPData "const
Teuchos::RefCountPtr<Ifpack_Preconditioner>&
MLAPI::InverseOperator::GetRCPData() const

Returns a pointer to the internally stored IFPACK preconditioner. ";

%feature("docstring")  MLAPI::InverseOperator::GetRCPMLPrec "const
Teuchos::RefCountPtr<ML_Epetra::MultiLevelPreconditioner>&
MLAPI::InverseOperator::GetRCPMLPrec() const

Returns a pointer to the internally stored ML preconditioner. ";

%feature("docstring")  MLAPI::InverseOperator::Apply "int
MLAPI::InverseOperator::Apply(const MultiVector &x, MultiVector &y)
const

Applies this object to vector lhs, returns values in rhs. ";

%feature("docstring")  MLAPI::InverseOperator::Print "ostream&
MLAPI::InverseOperator::Print(std::ostream &os, const bool
verbose=true) const

Prints out basic information about this object. ";


// File: classMLAPI_1_1LinearCombinationAdd.xml
%feature("docstring") MLAPI::LinearCombinationAdd "";

%feature("docstring")
MLAPI::LinearCombinationAdd::LinearCombinationAdd "MLAPI::LinearCombinationAdd::LinearCombinationAdd(const
BaseLinearCombination &left, const BaseLinearCombination &right) ";

%feature("docstring")  MLAPI::LinearCombinationAdd::GetVectorSpace "const Space MLAPI::LinearCombinationAdd::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")  MLAPI::LinearCombinationAdd::Update "void
MLAPI::LinearCombinationAdd::Update(MultiVector &v) const ";

%feature("docstring")  MLAPI::LinearCombinationAdd::Set "void
MLAPI::LinearCombinationAdd::Set(MultiVector &v) const ";


// File: classMLAPI_1_1LinearCombinationMixed.xml
%feature("docstring") MLAPI::LinearCombinationMixed "";

%feature("docstring")
MLAPI::LinearCombinationMixed::LinearCombinationMixed "MLAPI::LinearCombinationMixed::LinearCombinationMixed(const
BaseLinearCombination &left, const MultiVector &right, double alpha)
";

%feature("docstring")  MLAPI::LinearCombinationMixed::GetVectorSpace "const Space MLAPI::LinearCombinationMixed::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")  MLAPI::LinearCombinationMixed::Update "void
MLAPI::LinearCombinationMixed::Update(MultiVector &v) const ";

%feature("docstring")  MLAPI::LinearCombinationMixed::Set "void
MLAPI::LinearCombinationMixed::Set(MultiVector &v) const ";


// File: classMLAPI_1_1LinearCombinationScaled.xml
%feature("docstring") MLAPI::LinearCombinationScaled "";

%feature("docstring")
MLAPI::LinearCombinationScaled::LinearCombinationScaled "MLAPI::LinearCombinationScaled::LinearCombinationScaled(const
BaseLinearCombination &left, double scalar) ";

%feature("docstring")  MLAPI::LinearCombinationScaled::GetVectorSpace
"const Space MLAPI::LinearCombinationScaled::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")  MLAPI::LinearCombinationScaled::Set "void
MLAPI::LinearCombinationScaled::Set(MultiVector &v) const ";

%feature("docstring")  MLAPI::LinearCombinationScaled::Update "void
MLAPI::LinearCombinationScaled::Update(MultiVector &v) const ";


// File: classMLAPI_1_1LoadBalanceInverseOperator.xml
%feature("docstring") MLAPI::LoadBalanceInverseOperator "";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::LoadBalanceInverseOperator "MLAPI::LoadBalanceInverseOperator::LoadBalanceInverseOperator()

Empty constructor. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::LoadBalanceInverseOperator "MLAPI::LoadBalanceInverseOperator::LoadBalanceInverseOperator(const
LoadBalanceInverseOperator &RHS)

Copy constructor. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::~LoadBalanceInverseOperator "virtual
MLAPI::LoadBalanceInverseOperator::~LoadBalanceInverseOperator()

Destructor. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::Reshape "void MLAPI::LoadBalanceInverseOperator::Reshape()

Resets this object. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::Reshape "void MLAPI::LoadBalanceInverseOperator::Reshape(Ifpack_Preconditioner
*prec, const LoadBalanceOperator &Op, const bool ownership)

Reshape with preconstructed smoother as Ifpack_Preconditioner. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::GetParticipation "virtual bool
MLAPI::LoadBalanceInverseOperator::GetParticipation() const

Returns a bool indicating whether this proc participates in the
operator application. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::GetOperatorRangeSpace "const Space
MLAPI::LoadBalanceInverseOperator::GetOperatorRangeSpace() const

Returns a reference to the range space of this object. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::GetOperatorDomainSpace "const
Space MLAPI::LoadBalanceInverseOperator::GetOperatorDomainSpace()
const

Returns a reference to the domain space of this object. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::GetRangeSpace "const Space
MLAPI::LoadBalanceInverseOperator::GetRangeSpace() const

Returns a reference to the range space of this object. ";

%feature("docstring")
MLAPI::LoadBalanceInverseOperator::GetDomainSpace "const Space
MLAPI::LoadBalanceInverseOperator::GetDomainSpace() const

Returns a reference to the domain space of this object. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::RCPRowMatrix
"const Teuchos::RCP<Epetra_RowMatrix>
MLAPI::LoadBalanceInverseOperator::RCPRowMatrix() const

Returns pointer of the internally stored ML_Epetra::RowMatrix object.
";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::RowMatrix "Epetra_RowMatrix* MLAPI::LoadBalanceInverseOperator::RowMatrix() const

Returns pointer of the internally stored ML_Epetra::RowMatrix object.
";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::GetOperator
"const LoadBalanceOperator&
MLAPI::LoadBalanceInverseOperator::GetOperator() const

Returns a reference to the Operator of which this object defines the
inverse. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::GetRCPData "Teuchos::RCP<Ifpack_Preconditioner>&
MLAPI::LoadBalanceInverseOperator::GetRCPData()

Returns a pointer to the internally stored IFPACK preconditioner. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::GetRCPData "const Teuchos::RCP<Ifpack_Preconditioner>&
MLAPI::LoadBalanceInverseOperator::GetRCPData() const

Returns a pointer to the internally stored IFPACK preconditioner. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::Apply "int
MLAPI::LoadBalanceInverseOperator::Apply(const MultiVector &x,
MultiVector &y) const

Applies this object to vector lhs, returns values in rhs. ";

%feature("docstring")  MLAPI::LoadBalanceInverseOperator::Print "ostream& MLAPI::LoadBalanceInverseOperator::Print(std::ostream &os,
const bool verbose=true) const

Prints out basic information about this object. ";


// File: classMLAPI_1_1LoadBalanceOperator.xml
%feature("docstring") MLAPI::LoadBalanceOperator "";

%feature("docstring")  MLAPI::LoadBalanceOperator::LoadBalanceOperator
"MLAPI::LoadBalanceOperator::LoadBalanceOperator()

Default constructor. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::LoadBalanceOperator
"MLAPI::LoadBalanceOperator::LoadBalanceOperator(const Space
&DomainSpace, const Space &RangeSpace, ML_Operator *Op, bool
Ownership=true, Teuchos::RefCountPtr< ML_Operator_Box >
AuxOp=Teuchos::null)

Constructor with given already computed ML_Operator pointer. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::LoadBalanceOperator
"MLAPI::LoadBalanceOperator::LoadBalanceOperator(const Space
&DomainSpace, const Space &RangeSpace, Epetra_RowMatrix *Matrix, bool
Ownership=true, Teuchos::RefCountPtr< ML_Operator_Box >
AuxOp=Teuchos::null)

Constructor with given already FillComplete()'d object. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::LoadBalanceOperator
"MLAPI::LoadBalanceOperator::LoadBalanceOperator(const
LoadBalanceOperator &RHS)

Copy constructor. ";

%feature("docstring")
MLAPI::LoadBalanceOperator::~LoadBalanceOperator "MLAPI::LoadBalanceOperator::~LoadBalanceOperator()

Destructor. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::Reshape "void
MLAPI::LoadBalanceOperator::Reshape()

Resets this object. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::Reshape "void
MLAPI::LoadBalanceOperator::Reshape(const Space &DomainSpace, const
Space &RangeSpace, ML_Operator *Op, bool Ownership=true,
Teuchos::RefCountPtr< ML_Operator_Box > AuxOp=Teuchos::null)

Reshape with given already computed ML_Operator pointer. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::Reshape "void
MLAPI::LoadBalanceOperator::Reshape(const Space &DomainSpace, const
Space &RangeSpace, Epetra_RowMatrix *Matrix, bool Ownership=true,
Teuchos::RCP< ML_Operator_Box > AuxOp=Teuchos::null)

Reshape with given already FillComplete()'d object. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetParticipation "virtual bool MLAPI::LoadBalanceOperator::GetParticipation() const

Returns a bool indicating whether this proc participates in the
operator application. ";

%feature("docstring")
MLAPI::LoadBalanceOperator::GetOperatorDomainSpace "const Space
MLAPI::LoadBalanceOperator::GetOperatorDomainSpace() const

Returns a reference to the internally stored domain space. ";

%feature("docstring")
MLAPI::LoadBalanceOperator::GetOperatorRangeSpace "const Space
MLAPI::LoadBalanceOperator::GetOperatorRangeSpace() const

Returns a reference to the internally stored range space. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetDomainSpace "const Space MLAPI::LoadBalanceOperator::GetDomainSpace() const

Returns a reference to the internally stored domain space. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetRangeSpace "const Space MLAPI::LoadBalanceOperator::GetRangeSpace() const

Returns a reference to the internally stored range space. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetColumnSpace "const Space MLAPI::LoadBalanceOperator::GetColumnSpace() const

Returns a reference to the internally stored column space. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetNumGlobalRows "int MLAPI::LoadBalanceOperator::GetNumGlobalRows() const

Returns the number of global rows. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetNumMyRows "int
MLAPI::LoadBalanceOperator::GetNumMyRows() const

Returns the number of local rows. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetNumGlobalCols "int MLAPI::LoadBalanceOperator::GetNumGlobalCols() const

Returns the number of global columns. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetNumMyCols "int
MLAPI::LoadBalanceOperator::GetNumMyCols() const

Returns the number of local columns. ";

%feature("docstring")
MLAPI::LoadBalanceOperator::GetNumGlobalNonzeros "int
MLAPI::LoadBalanceOperator::GetNumGlobalNonzeros() const

Returns the global number of nonzeros. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetNumMyNonzeros "int MLAPI::LoadBalanceOperator::GetNumMyNonzeros() const

Returns the local number of nonzeros. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetRowMatrix "const Epetra_RowMatrix* MLAPI::LoadBalanceOperator::GetRowMatrix()
const

Returns the RefCountPtr of OperatorBox_. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetML_Operator "ML_Operator* MLAPI::LoadBalanceOperator::GetML_Operator() const

Returns the RefCountPtr of OperatorBox_. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetRCPOperatorBox "const Teuchos::RCP<ML_Operator_Box>&
MLAPI::LoadBalanceOperator::GetRCPOperatorBox() const

Returns the RefCountPtr of OperatorBox_. ";

%feature("docstring")
MLAPI::LoadBalanceOperator::GetRCPAuxOperatorBox "const
Teuchos::RCP<ML_Operator_Box>&
MLAPI::LoadBalanceOperator::GetRCPAuxOperatorBox() const

Returns the RefCountPtr of AuxOperatorBox_. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetRCPRowMatrix "const Teuchos::RCP<Epetra_RowMatrix>&
MLAPI::LoadBalanceOperator::GetRCPRowMatrix() const

Returns the RefCountPtr of RowMatrix_. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetGRID "int
MLAPI::LoadBalanceOperator::GetGRID(const int LRID) const

Returns the global ID of local row ID LRID. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::GetGCID "int
MLAPI::LoadBalanceOperator::GetGCID(const int LCID) const

Returns the global ID of local column ID LCID. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::Apply "int
MLAPI::LoadBalanceOperator::Apply(const MultiVector &X, MultiVector
&Y) const

Applies this operator to LHS, returns the result in RHS. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::Print "ostream&
MLAPI::LoadBalanceOperator::Print(std::ostream &os, const bool
verbose=true) const

Prints basic information about this object. ";

%feature("docstring")  MLAPI::LoadBalanceOperator::BuildColumnSpace "void MLAPI::LoadBalanceOperator::BuildColumnSpace()

Build the column space, by computing the GID of all local columns. ";


// File: classMLAPI_1_1MATLABStream.xml
%feature("docstring") MLAPI::MATLABStream "

Basic stream to save in a MATLAB-compatible file MLAPI objects.

For an example of usage, see ml_blackboard_cpp

Marzio Sala, SNL 9214

C++ includes: MLAPI_MATLABStream.h ";

%feature("docstring")  MLAPI::MATLABStream::MATLABStream "MLAPI::MATLABStream::MATLABStream(const string &FileName, bool
UseSparse=true)

Opens the specified file for writing. ";

%feature("docstring")  MLAPI::MATLABStream::~MATLABStream "MLAPI::MATLABStream::~MATLABStream()

Finally closes the output file. ";

%feature("docstring")  MLAPI::MATLABStream::GetUseSparse "bool
MLAPI::MATLABStream::GetUseSparse() const

Returns true if the stream uses sparse MATLAB format. ";

%feature("docstring")  MLAPI::MATLABStream::SetUseSparse "void
MLAPI::MATLABStream::SetUseSparse(const bool UseSparse)

Toggles the use of sparse MATLAB formats. ";

%feature("docstring")  MLAPI::MATLABStream::GetFileName "string
MLAPI::MATLABStream::GetFileName() const

Returns the name of the output file. ";


// File: classMLAPI_1_1ML__Operator__Box.xml
%feature("docstring") MLAPI::ML_Operator_Box "

Simple wrapper for ML_Operator struct.

Marzio Sala, SNL 9214.

C++ includes: MLAPI_Operator_Box.h ";

%feature("docstring")  MLAPI::ML_Operator_Box::ML_Operator_Box "MLAPI::ML_Operator_Box::ML_Operator_Box(ML_Operator *Op, bool
Ownership=true)

Constructor. ";

%feature("docstring")  MLAPI::ML_Operator_Box::~ML_Operator_Box "MLAPI::ML_Operator_Box::~ML_Operator_Box()

Destructor. ";

%feature("docstring")  MLAPI::ML_Operator_Box::GetData "ML_Operator*
MLAPI::ML_Operator_Box::GetData() const

Returns a pointer to the internally stored ML_Operator. ";


// File: classMLAPI_1_1MultiLevelAdaptiveSA.xml
%feature("docstring") MLAPI::MultiLevelAdaptiveSA "

Black-box multilevel adaptive smoothed aggregation preconditioner.

This class implements an adaptive smoothed aggregation preconditioner.
An example of usage is reported in file ml_adaptivesa. We note that
the usage of this class is slightly different from that of
MultiLevelSA.

An instance of this class can be created as follows:

Important methods of this class: The number of PDE equations on the
finest level can be queried using GetInputNumPDEEqns().

GetNumPDEEqns() returns the number of PDE equations on the current
level. This value can be set via SetNumPDEEqns().

GetNullSpace() returns a reference to the internally stored null
space; the null space is set using SetNullSpace().

GetMaxLevels() returns the number of levels. If called before
Compute(), GetMaxLevels() returns the maximum number of levels used in
the constructor, otherwise returns the actual number of levels.

GetSmootherType() and GetCoarseType() return the smoother and coarse
type.

The number of application of the cycle in IncrementNullSpace() is
given by GetNumItersCoarse() and GetNumItersFine().

Methods A(level), P(level), R(level) and S(level) return a reference
to the internally stored operators.

Method SetList() can be used at any time to reset the internally
stored list.

The general usage is: Specify the null space using SetNullSpace(NS),
where NS is a MultiVector, then compute the hierarchy using Compute(),
or

Compute the first component of the null space using
SetupInitialNullSpace(). This will define a single-vector null space,
and store it using SetNullSpace(NS).

When a non-empty null space is provided, the user can increment by one
the dimension of the null space by calling IncrementNullSpace().

Method AdaptCompute() performs all these operations.

Marzio Sala, Ray Tuminaro, Jonathan Hu, Michael Gee, Marian Brezina.

C++ includes: MLAPI_MultiLevelAdaptiveSA.h ";

%feature("docstring")
MLAPI::MultiLevelAdaptiveSA::MultiLevelAdaptiveSA "MLAPI::MultiLevelAdaptiveSA::MultiLevelAdaptiveSA(const Operator
FineMatrix, Teuchos::ParameterList &List, const int NumPDEEqns, const
int MaxLevels=20)

Constructs the hierarchy for given Operator and parameters. ";

%feature("docstring")
MLAPI::MultiLevelAdaptiveSA::~MultiLevelAdaptiveSA "virtual
MLAPI::MultiLevelAdaptiveSA::~MultiLevelAdaptiveSA()

Destructor. ";

%feature("docstring")
MLAPI::MultiLevelAdaptiveSA::GetOperatorDomainSpace "const Space
MLAPI::MultiLevelAdaptiveSA::GetOperatorDomainSpace() const

Returns a copy of the internally stored domain space. ";

%feature("docstring")
MLAPI::MultiLevelAdaptiveSA::GetOperatorRangeSpace "const Space
MLAPI::MultiLevelAdaptiveSA::GetOperatorRangeSpace() const

Returns a copy of the internally stored range space. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetDomainSpace "const Space MLAPI::MultiLevelAdaptiveSA::GetDomainSpace() const

Returns a copy of the internally stored domain space. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetRangeSpace "const Space MLAPI::MultiLevelAdaptiveSA::GetRangeSpace() const

Returns a copy of the internally stored range space. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::R "Operator&
MLAPI::MultiLevelAdaptiveSA::R(const int i)

Returns a reference to the restriction operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::R "const
Operator& MLAPI::MultiLevelAdaptiveSA::R(const int i) const

Returns a reference to the restriction operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::A "Operator&
MLAPI::MultiLevelAdaptiveSA::A(const int i)

Returns a reference to the operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::A "const
Operator& MLAPI::MultiLevelAdaptiveSA::A(const int i) const

Returns a reference to the operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::P "Operator&
MLAPI::MultiLevelAdaptiveSA::P(const int i)

Returns a reference to the prolongator operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::P "const
Operator& MLAPI::MultiLevelAdaptiveSA::P(const int i) const

Returns a reference to the prolongator operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::S "InverseOperator& MLAPI::MultiLevelAdaptiveSA::S(const int i)

Returns a reference to the inverse operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::S "const
InverseOperator& MLAPI::MultiLevelAdaptiveSA::S(const int i) const

Returns a reference to the inverse operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetMaxLevels "int
MLAPI::MultiLevelAdaptiveSA::GetMaxLevels() const

Returns the actual number of levels. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::SetMaxLevels "void MLAPI::MultiLevelAdaptiveSA::SetMaxLevels(const int MaxLevels)

Returns the actual number of levels. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetNullSpace "const MultiVector MLAPI::MultiLevelAdaptiveSA::GetNullSpace() const

Gets a reference to the internally stored null space. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::SetNullSpace "void MLAPI::MultiLevelAdaptiveSA::SetNullSpace(MultiVector &NullSpace)

Sets the null space multi-vector to NullSpace. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::IsComputed "bool
MLAPI::MultiLevelAdaptiveSA::IsComputed() const

Returns true if the hierarchy has been successfully computed. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::SetList "void
MLAPI::MultiLevelAdaptiveSA::SetList(Teuchos::ParameterList &List)

Sets the internally stored list to List. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetSmootherType "string MLAPI::MultiLevelAdaptiveSA::GetSmootherType()

Returns the smoother solver type. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetCoarseType "string MLAPI::MultiLevelAdaptiveSA::GetCoarseType()

Returns the coarse solver type. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::SetInputNumPDEEqns
"void MLAPI::MultiLevelAdaptiveSA::SetInputNumPDEEqns(const int n)

Returns the number of PDE equations on the finest level. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetInputNumPDEEqns
"int MLAPI::MultiLevelAdaptiveSA::GetInputNumPDEEqns()

Returns the number of PDE equations on the current level. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetNumPDEEqns "int MLAPI::MultiLevelAdaptiveSA::GetNumPDEEqns()

Sets the number of PDE equations on the current level. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::SetNumPDEEqns "void MLAPI::MultiLevelAdaptiveSA::SetNumPDEEqns(const int NumPDEEqns)
";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetMaxCoarseSize "int MLAPI::MultiLevelAdaptiveSA::GetMaxCoarseSize()

Returns the maximum allowed coarse size. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetMaxReduction "double MLAPI::MultiLevelAdaptiveSA::GetMaxReduction()

Returns the maximum allowed reduction. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetNumItersCoarse
"int MLAPI::MultiLevelAdaptiveSA::GetNumItersCoarse()

Returns the maximum number of applications on the coarser levels. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetNumItersFine "int MLAPI::MultiLevelAdaptiveSA::GetNumItersFine()

Returns the maximum number of applications on the finest level. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::GetComplexity "double MLAPI::MultiLevelAdaptiveSA::GetComplexity()

Returns the multigrid preconditioner operator complexity. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::Compute "void
MLAPI::MultiLevelAdaptiveSA::Compute()

Creates an hierarchy using the provided or default null space. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::AdaptCompute "void MLAPI::MultiLevelAdaptiveSA::AdaptCompute(const bool
UseDefaultOrSpecified, int AdditionalCandidates)

Setup the adaptive multilevel hierarchy. ";

%feature("docstring")
MLAPI::MultiLevelAdaptiveSA::SetupInitialNullSpace "void
MLAPI::MultiLevelAdaptiveSA::SetupInitialNullSpace()

Computes the first component of the null space. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::IncrementNullSpace
"bool MLAPI::MultiLevelAdaptiveSA::IncrementNullSpace()

Increments the null space dimension by one. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::Apply "int
MLAPI::MultiLevelAdaptiveSA::Apply(const MultiVector &b_f, MultiVector
&x_f) const

Applies the preconditioner to b_f, returns the result in x_f. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::SolveMultiLevelSA
"int MLAPI::MultiLevelAdaptiveSA::SolveMultiLevelSA(const MultiVector
&b_f, MultiVector &x_f, int level) const

Recursively called core of the multi level preconditioner. ";

%feature("docstring")  MLAPI::MultiLevelAdaptiveSA::Print "std::ostream& MLAPI::MultiLevelAdaptiveSA::Print(std::ostream &os,
const bool verbose=true) const

Prints basic information about this preconditioner. ";


// File: classMultiLevelPreconditioner.xml
%feature("docstring") MultiLevelPreconditioner "

ML black-box preconditioner for Epetra_RowMatrix derived classes.

C++ includes: ml_MultiLevelPreconditioner.h ";


// File: classMLAPI_1_1MultiLevelSA.xml
%feature("docstring") MLAPI::MultiLevelSA "

Black-box multilevel smoothed aggregation preconditioner.

Marzio Sala, SNL 9214

C++ includes: MLAPI_MultiLevelSA.h ";

%feature("docstring")  MLAPI::MultiLevelSA::MultiLevelSA "MLAPI::MultiLevelSA::MultiLevelSA(const Operator FineMatrix,
Teuchos::ParameterList &List, const bool ConstructNow=true)

Constructs the hierarchy for given Operator and parameters. ";

%feature("docstring")  MLAPI::MultiLevelSA::~MultiLevelSA "virtual
MLAPI::MultiLevelSA::~MultiLevelSA()

Destructor. ";

%feature("docstring")  MLAPI::MultiLevelSA::GetOperatorDomainSpace "const Space MLAPI::MultiLevelSA::GetOperatorDomainSpace() const

Returns a copy of the internally stored domain space. ";

%feature("docstring")  MLAPI::MultiLevelSA::GetOperatorRangeSpace "const Space MLAPI::MultiLevelSA::GetOperatorRangeSpace() const

Returns a copy of the internally stored range space. ";

%feature("docstring")  MLAPI::MultiLevelSA::GetDomainSpace "const
Space MLAPI::MultiLevelSA::GetDomainSpace() const

Returns a copy of the internally stored domain space. ";

%feature("docstring")  MLAPI::MultiLevelSA::GetRangeSpace "const
Space MLAPI::MultiLevelSA::GetRangeSpace() const

Returns a copy of the internally stored range space. ";

%feature("docstring")  MLAPI::MultiLevelSA::R "const Operator&
MLAPI::MultiLevelSA::R(const int i) const

Returns a reference to the restriction operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelSA::A "const Operator&
MLAPI::MultiLevelSA::A(const int i) const

Returns a reference to the operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelSA::P "const Operator&
MLAPI::MultiLevelSA::P(const int i) const

Returns a reference to the prolongator operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelSA::S "const InverseOperator&
MLAPI::MultiLevelSA::S(const int i) const

Returns a reference to the inverse operator of level i. ";

%feature("docstring")  MLAPI::MultiLevelSA::GetMaxLevels "int
MLAPI::MultiLevelSA::GetMaxLevels() const

Returns the actual number of levels. ";

%feature("docstring")  MLAPI::MultiLevelSA::IsComputed "bool
MLAPI::MultiLevelSA::IsComputed() const

Returns true if the hierarchy has been successfully computed, false
otherwise. ";

%feature("docstring")  MLAPI::MultiLevelSA::Compute "void
MLAPI::MultiLevelSA::Compute()

Computes the hierarchy. ";

%feature("docstring")  MLAPI::MultiLevelSA::Apply "int
MLAPI::MultiLevelSA::Apply(const MultiVector &b_f, MultiVector &x_f)
const

Applies the preconditioner to b_f, returns the result in x_f. ";

%feature("docstring")  MLAPI::MultiLevelSA::SolveMultiLevelSA "int
MLAPI::MultiLevelSA::SolveMultiLevelSA(const MultiVector &b_f,
MultiVector &x_f, int level) const

Recursively called core of the multi level preconditioner. ";

%feature("docstring")  MLAPI::MultiLevelSA::Print "std::ostream&
MLAPI::MultiLevelSA::Print(std::ostream &os, const bool verbose=true)
const

Prints basic information about this preconditioner. ";


// File: classMLAPI_1_1MultiVector.xml
%feature("docstring") MLAPI::MultiVector "

Basic class for distributed double-precision vectors.

Marzio Sala, SNL 9214.

C++ includes: MLAPI_MultiVector.h ";

%feature("docstring")  MLAPI::MultiVector::MultiVector "MLAPI::MultiVector::MultiVector()

Default constructor. ";

%feature("docstring")  MLAPI::MultiVector::MultiVector "MLAPI::MultiVector::MultiVector(const Space &VectorSpace, const int
NumVectors=1, bool SetToZero=true)

Constructor for a given Space. ";

%feature("docstring")  MLAPI::MultiVector::MultiVector "MLAPI::MultiVector::MultiVector(const Space &VectorSpace, double
**Values, const int NumVectors=1)

Constructor with a given Space, and user-provided array of values. ";

%feature("docstring")  MLAPI::MultiVector::MultiVector "MLAPI::MultiVector::MultiVector(const Space &VectorSpace,
Teuchos::RefCountPtr< DoubleVector > RCPValues)

Constructor with a given Space, and user-provided RefCountPtr. ";

%feature("docstring")  MLAPI::MultiVector::MultiVector "MLAPI::MultiVector::MultiVector(const Space &VectorSpace, std::vector<
Teuchos::RefCountPtr< DoubleVector > > RCPValues)

Constructor with a given Space, and user-provided array of values. ";

%feature("docstring")  MLAPI::MultiVector::MultiVector "MLAPI::MultiVector::MultiVector(const MultiVector &rhs)

Copy constructor. ";

%feature("docstring")  MLAPI::MultiVector::~MultiVector "MLAPI::MultiVector::~MultiVector()

Destructor. ";

%feature("docstring")  MLAPI::MultiVector::Reshape "void
MLAPI::MultiVector::Reshape()

Resets this object. ";

%feature("docstring")  MLAPI::MultiVector::Reshape "void
MLAPI::MultiVector::Reshape(const Space &S, const int NumVectors=1,
const bool SetToZero=true)

Sets the space of this vector. ";

%feature("docstring")  MLAPI::MultiVector::Append "void
MLAPI::MultiVector::Append(const int NumVectors=1, const bool
SetToZero=true)

Appends a new vector. ";

%feature("docstring")  MLAPI::MultiVector::Append "void
MLAPI::MultiVector::Append(MultiVector rhs)

Appends a new vector. ";

%feature("docstring")  MLAPI::MultiVector::Delete "void
MLAPI::MultiVector::Delete(const int v)

Deletes the last vector. ";

%feature("docstring")  MLAPI::MultiVector::GetVectorSpace "const
Space& MLAPI::MultiVector::GetVectorSpace() const

Returns the Space on which this vector is defined. ";

%feature("docstring")  MLAPI::MultiVector::GetVectorSpace "Space&
MLAPI::MultiVector::GetVectorSpace()

Returns the Space on which this vector is defined (non-const) ";

%feature("docstring")  MLAPI::MultiVector::GetNumVectors "int
MLAPI::MultiVector::GetNumVectors() const

Returns the number of vectors. ";

%feature("docstring")  MLAPI::MultiVector::GetMyLength "int
MLAPI::MultiVector::GetMyLength() const

Returns the local length of each vector. ";

%feature("docstring")  MLAPI::MultiVector::GetGlobalLength "int
MLAPI::MultiVector::GetGlobalLength() const

Returns the global length of each vector. ";

%feature("docstring")  MLAPI::MultiVector::GetValues "double*
MLAPI::MultiVector::GetValues(const int v)

Returns a pointer to the double array (non-const version) ";

%feature("docstring")  MLAPI::MultiVector::GetValues "const double*
MLAPI::MultiVector::GetValues(const int v) const

Returns a pointer to the double array (const version) ";

%feature("docstring")  MLAPI::MultiVector::GetRCPValues "Teuchos::RefCountPtr<DoubleVector>&
MLAPI::MultiVector::GetRCPValues(const int v)

Returns a pointer to the double array (non-const version) ";

%feature("docstring")  MLAPI::MultiVector::GetRCPValues "const
Teuchos::RefCountPtr<DoubleVector>&
MLAPI::MultiVector::GetRCPValues(const int v) const

Returns a pointer to the double array (const version) ";

%feature("docstring")  MLAPI::MultiVector::SetRCPValues "void
MLAPI::MultiVector::SetRCPValues(const Teuchos::RefCountPtr<
DoubleVector > &RCPValues, const int v)

Sets the RefCountPtr<Values_> ";

%feature("docstring")  MLAPI::MultiVector::Update "void
MLAPI::MultiVector::Update(const double alpha, int v=-1)

Sets this(v) = rhs. ";

%feature("docstring")  MLAPI::MultiVector::Update "void
MLAPI::MultiVector::Update(const MultiVector &rhs)

Sets this = rhs. ";

%feature("docstring")  MLAPI::MultiVector::Update "void
MLAPI::MultiVector::Update(double alpha, const MultiVector &rhs)

Sets this = alpha * rhs. ";

%feature("docstring")  MLAPI::MultiVector::Update "void
MLAPI::MultiVector::Update(double alpha, const MultiVector &x, double
beta, const MultiVector &y)

Sets this = alpha * x + beta * y. ";

%feature("docstring")  MLAPI::MultiVector::Update "void
MLAPI::MultiVector::Update(double alpha, const MultiVector &rhs,
double beta)

Sets this = alpha * rhs + beta * this. ";

%feature("docstring")  MLAPI::MultiVector::DotProduct "double
MLAPI::MultiVector::DotProduct(const MultiVector &rhs, int v=-1) const

Computes the dot product between this vector and rhs. ";

%feature("docstring")  MLAPI::MultiVector::Norm2 "double
MLAPI::MultiVector::Norm2(int v=-1) const

Computes the 2-norm of this vector. ";

%feature("docstring")  MLAPI::MultiVector::NormInf "double
MLAPI::MultiVector::NormInf(int v=-1) const

Computes the infinite norm of this vector. ";

%feature("docstring")  MLAPI::MultiVector::NormOne "double
MLAPI::MultiVector::NormOne(int v=-1) const

Computes the one norm of this vector. ";

%feature("docstring")  MLAPI::MultiVector::Reciprocal "void
MLAPI::MultiVector::Reciprocal(int v=-1)

Replaces each element of the vector with its reciprocal. ";

%feature("docstring")  MLAPI::MultiVector::Scale "void
MLAPI::MultiVector::Scale(const double Factor, int v=-1)

Scales each element by the specified factor. ";

%feature("docstring")  MLAPI::MultiVector::Random "void
MLAPI::MultiVector::Random(int v=-1)

Populates the vector with random elements. ";

%feature("docstring")  MLAPI::MultiVector::Sort "void
MLAPI::MultiVector::Sort(int v=-1, const bool IsIncreasing=false)

Sorts the component of the vector. ";

%feature("docstring")  MLAPI::MultiVector::Print "virtual
std::ostream& MLAPI::MultiVector::Print(std::ostream &os, const bool
verbose=true) const

Prints basic information about this object on ostream. ";

%feature("docstring")  MLAPI::MultiVector::IsAlias "bool
MLAPI::MultiVector::IsAlias(const MultiVector &rhs) const ";


// File: classMLAPI_1_1MultiVectorCombination.xml
%feature("docstring") MLAPI::MultiVectorCombination "";

%feature("docstring")
MLAPI::MultiVectorCombination::MultiVectorCombination "MLAPI::MultiVectorCombination::MultiVectorCombination(const double
alpha, const MultiVector x, const double beta, const MultiVector y) ";

%feature("docstring")  MLAPI::MultiVectorCombination::GetVectorSpace "const Space MLAPI::MultiVectorCombination::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")
MLAPI::MultiVectorCombination::GetLeftMultiVector "const MultiVector
MLAPI::MultiVectorCombination::GetLeftMultiVector() const ";

%feature("docstring")  MLAPI::MultiVectorCombination::GetLeftScalar "double MLAPI::MultiVectorCombination::GetLeftScalar() const ";

%feature("docstring")
MLAPI::MultiVectorCombination::GetRightMultiVector "const MultiVector
MLAPI::MultiVectorCombination::GetRightMultiVector() const ";

%feature("docstring")  MLAPI::MultiVectorCombination::GetRightScalar "double MLAPI::MultiVectorCombination::GetRightScalar() const ";

%feature("docstring")  MLAPI::MultiVectorCombination::Update "void
MLAPI::MultiVectorCombination::Update(MultiVector &v) const ";

%feature("docstring")  MLAPI::MultiVectorCombination::Set "void
MLAPI::MultiVectorCombination::Set(MultiVector &v) const ";


// File: classMLAPI_1_1MultiVectorScaled.xml
%feature("docstring") MLAPI::MultiVectorScaled "";

%feature("docstring")  MLAPI::MultiVectorScaled::MultiVectorScaled "MLAPI::MultiVectorScaled::MultiVectorScaled(const MultiVector &vector,
const double alpha) ";

%feature("docstring")  MLAPI::MultiVectorScaled::GetVectorSpace "const Space MLAPI::MultiVectorScaled::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")  MLAPI::MultiVectorScaled::GetMultiVector "const MultiVector& MLAPI::MultiVectorScaled::GetMultiVector() const ";

%feature("docstring")  MLAPI::MultiVectorScaled::GetScalar "double
MLAPI::MultiVectorScaled::GetScalar() const ";

%feature("docstring")  MLAPI::MultiVectorScaled::Update "void
MLAPI::MultiVectorScaled::Update(MultiVector &v) const ";

%feature("docstring")  MLAPI::MultiVectorScaled::Set "void
MLAPI::MultiVectorScaled::Set(MultiVector &v) const ";


// File: classMLAPI_1_1Operator.xml
%feature("docstring") MLAPI::Operator "

Operator: basic class to define operators within MLAPI.

Michael Gee, TU Munich.

Marzio Sala, SNL 9214

C++ includes: MLAPI_Operator.h ";

%feature("docstring")  MLAPI::Operator::Operator "MLAPI::Operator::Operator()

Default constructor. ";

%feature("docstring")  MLAPI::Operator::Operator "MLAPI::Operator::Operator(const Space &DomainSpace, const Space
&RangeSpace, ML_Operator *Op, bool Ownership=true,
Teuchos::RefCountPtr< ML_Operator_Box > AuxOp=Teuchos::null)

Constructor with given already computed ML_Operator pointer. ";

%feature("docstring")  MLAPI::Operator::Operator "MLAPI::Operator::Operator(const Space &DomainSpace, const Space
&RangeSpace, Epetra_RowMatrix *Matrix, bool Ownership=true,
Teuchos::RefCountPtr< ML_Operator_Box > AuxOp=Teuchos::null)

Constructor with given already FillComplete()'d object. ";

%feature("docstring")  MLAPI::Operator::Operator "MLAPI::Operator::Operator(const Operator &RHS)

Copy constructor. ";

%feature("docstring")  MLAPI::Operator::~Operator "MLAPI::Operator::~Operator()

Destructor. ";

%feature("docstring")  MLAPI::Operator::Reshape "void
MLAPI::Operator::Reshape()

Resets this object. ";

%feature("docstring")  MLAPI::Operator::Reshape "void
MLAPI::Operator::Reshape(const Space &DomainSpace, const Space
&RangeSpace, ML_Operator *Op, bool Ownership=true,
Teuchos::RefCountPtr< ML_Operator_Box > AuxOp=Teuchos::null)

Reshape with given already computed ML_Operator pointer. ";

%feature("docstring")  MLAPI::Operator::Reshape "void
MLAPI::Operator::Reshape(const Space &DomainSpace, const Space
&RangeSpace, Epetra_RowMatrix *Matrix, bool Ownership=true,
Teuchos::RefCountPtr< ML_Operator_Box > AuxOp=Teuchos::null)

Reshape with given already FillComplete()'d object. ";

%feature("docstring")  MLAPI::Operator::GetOperatorDomainSpace "const
Space MLAPI::Operator::GetOperatorDomainSpace() const

Returns a reference to the internally stored domain space. ";

%feature("docstring")  MLAPI::Operator::GetOperatorRangeSpace "const
Space MLAPI::Operator::GetOperatorRangeSpace() const

Returns a reference to the internally stored range space. ";

%feature("docstring")  MLAPI::Operator::GetDomainSpace "const Space
MLAPI::Operator::GetDomainSpace() const

Returns a reference to the internally stored domain space. ";

%feature("docstring")  MLAPI::Operator::GetRangeSpace "const Space
MLAPI::Operator::GetRangeSpace() const

Returns a reference to the internally stored range space. ";

%feature("docstring")  MLAPI::Operator::GetColumnSpace "const Space
MLAPI::Operator::GetColumnSpace() const

Returns a reference to the internally stored column space. ";

%feature("docstring")  MLAPI::Operator::GetNumGlobalRows "int
MLAPI::Operator::GetNumGlobalRows() const

Returns the number of global rows. ";

%feature("docstring")  MLAPI::Operator::GetNumMyRows "int
MLAPI::Operator::GetNumMyRows() const

Returns the number of local rows. ";

%feature("docstring")  MLAPI::Operator::GetNumGlobalCols "int
MLAPI::Operator::GetNumGlobalCols() const

Returns the number of global columns. ";

%feature("docstring")  MLAPI::Operator::GetNumMyCols "int
MLAPI::Operator::GetNumMyCols() const

Returns the number of local columns. ";

%feature("docstring")  MLAPI::Operator::GetNumGlobalNonzeros "int
MLAPI::Operator::GetNumGlobalNonzeros() const

Returns the global number of nonzeros. ";

%feature("docstring")  MLAPI::Operator::GetNumMyNonzeros "int
MLAPI::Operator::GetNumMyNonzeros() const

Returns the local number of nonzeros. ";

%feature("docstring")  MLAPI::Operator::GetRowMatrix "const
Epetra_RowMatrix* MLAPI::Operator::GetRowMatrix() const

Returns the RefCountPtr of OperatorBox_. ";

%feature("docstring")  MLAPI::Operator::GetML_Operator "ML_Operator*
MLAPI::Operator::GetML_Operator() const

Returns the RefCountPtr of OperatorBox_. ";

%feature("docstring")  MLAPI::Operator::GetRCPOperatorBox "const
Teuchos::RefCountPtr<ML_Operator_Box>&
MLAPI::Operator::GetRCPOperatorBox() const

Returns the RefCountPtr of OperatorBox_. ";

%feature("docstring")  MLAPI::Operator::GetRCPAuxOperatorBox "const
Teuchos::RefCountPtr<ML_Operator_Box>&
MLAPI::Operator::GetRCPAuxOperatorBox() const

Returns the RefCountPtr of AuxOperatorBox_. ";

%feature("docstring")  MLAPI::Operator::GetRCPRowMatrix "const
Teuchos::RefCountPtr<Epetra_RowMatrix>&
MLAPI::Operator::GetRCPRowMatrix() const

Returns the RefCountPtr of RowMatrix_. ";

%feature("docstring")  MLAPI::Operator::GetGRID "int
MLAPI::Operator::GetGRID(const int LRID) const

Returns the global ID of local row ID LRID. ";

%feature("docstring")  MLAPI::Operator::GetGCID "int
MLAPI::Operator::GetGCID(const int LCID) const

Returns the global ID of local column ID LCID. ";

%feature("docstring")  MLAPI::Operator::Apply "int
MLAPI::Operator::Apply(const MultiVector &X, MultiVector &Y) const

Applies this operator to LHS, returns the result in RHS. ";

%feature("docstring")  MLAPI::Operator::Print "ostream&
MLAPI::Operator::Print(std::ostream &os, const bool verbose=true)
const

Prints basic information about this object. ";

%feature("docstring")  MLAPI::Operator::BuildColumnSpace "void
MLAPI::Operator::BuildColumnSpace()

Build the column space, by computing the GID of all local columns. ";


// File: classMLAPI_1_1Residual.xml
%feature("docstring") MLAPI::Residual "";

%feature("docstring")  MLAPI::Residual::Residual "MLAPI::Residual::Residual(double alpha, const MultiVector &b, double
beta, const BaseOperator &A, const MultiVector &x) ";

%feature("docstring")  MLAPI::Residual::GetVectorSpace "const Space
MLAPI::Residual::GetVectorSpace() const

Returns the vector space of the underlying object. ";

%feature("docstring")  MLAPI::Residual::Update "void
MLAPI::Residual::Update(MultiVector &v) const ";

%feature("docstring")  MLAPI::Residual::Set "void
MLAPI::Residual::Set(MultiVector &v) const ";


// File: classMLAPI_1_1SerialMatrix.xml
%feature("docstring") MLAPI::SerialMatrix "";

%feature("docstring")  MLAPI::SerialMatrix::SerialMatrix "MLAPI::SerialMatrix::SerialMatrix() ";

%feature("docstring")  MLAPI::SerialMatrix::SerialMatrix "MLAPI::SerialMatrix::SerialMatrix(const Space &RowSpace, const Space
&ColSpace) ";

%feature("docstring")  MLAPI::SerialMatrix::Print "std::ostream&
MLAPI::SerialMatrix::Print(std::ostream &os, const bool verbose=true)
const

Prints basic information about this object. ";


// File: classMLAPI_1_1Space.xml
%feature("docstring") MLAPI::Space "

Specifies the number and distribution among processes of elements.

Marzio Sala, SNL 9214

C++ includes: MLAPI_Space.h ";

%feature("docstring")  MLAPI::Space::Space "MLAPI::Space::Space()

Default constructor, defines an empty space. ";

%feature("docstring")  MLAPI::Space::Space "MLAPI::Space::Space(const
int NumGlobalElements, const int NumMyElements=-1)

Constructor with specified number of global and local elements.

Constructs a space with linear distribution.

Parameters:
-----------

NumGlobalElements:  - (In) number of global elements.

NumMyElements:  - (In) number of local elements. If different from -1,
then eithere NumGlobalElements == -1, or the sum across of processors
of NumMyElements equals NumGlobalElements ";

%feature("docstring")  MLAPI::Space::Space "MLAPI::Space::Space(const
Epetra_Map &Map)

Constructor with specified Epetra_Map. ";

%feature("docstring")  MLAPI::Space::Space "MLAPI::Space::Space(const
int NumGlobalElements, const int NumMyElements, const int
*MyGlobalElements)

Constructor for non-linear distributions.

Parameters:
-----------

NumGlobalElements:  - (In) number of global elements. Set to -1 to
compute it automatically.

NumMyElements:  - (In) number of local elements. Cannot be set to -1.

MyGlobalElements:  - (In) contains the global ID of each local node.

Global ID always starts from 0. ";

%feature("docstring")  MLAPI::Space::Space "MLAPI::Space::Space(const
Space &RHS)

Copy constructor. ";

%feature("docstring")  MLAPI::Space::~Space "MLAPI::Space::~Space()

Destructor. ";

%feature("docstring")  MLAPI::Space::Reshape "void
MLAPI::Space::Reshape()

Resets this object. ";

%feature("docstring")  MLAPI::Space::Reshape "void
MLAPI::Space::Reshape(const int NumGlobalElements, const int
NumMyElements=-1)

Resets the dimension of the space by specifying the local number of
elements. ";

%feature("docstring")  MLAPI::Space::Reshape "void
MLAPI::Space::Reshape(const int NumGlobalElements, const int
NumMyElements, const int *MyGlobalElements)

Reset the dimension of the space by specifying the local number of
elements and their global numbering (starting from 0). ";

%feature("docstring")  MLAPI::Space::GetNumMyElements "int
MLAPI::Space::GetNumMyElements() const

Returns the local number of elements on the calling process. ";

%feature("docstring")  MLAPI::Space::GetNumGlobalElements "int
MLAPI::Space::GetNumGlobalElements() const

Returns the global number of elements. ";

%feature("docstring")  MLAPI::Space::GetOffset "int
MLAPI::Space::GetOffset() const

Returns the global ID of the first element on the calling process (for
linear distributions only). ";

%feature("docstring")  MLAPI::Space::IsLinear "bool
MLAPI::Space::IsLinear() const

Returns true if the decomposition among processors is linear. ";

%feature("docstring")  MLAPI::Space::GetRCPMyGlobalElements "const
Teuchos::RefCountPtr<Epetra_IntSerialDenseVector>
MLAPI::Space::GetRCPMyGlobalElements() const

Returns a pointer to the list of global nodes. ";

%feature("docstring")  MLAPI::Space::Print "std::ostream&
MLAPI::Space::Print(std::ostream &os, const bool verbose=true) const

Prints on ostream basic information about this object. ";


// File: structMLAPI_1_1StackEntry.xml
%feature("docstring") MLAPI::StackEntry "";


// File: classMLAPI_1_1TimeObject.xml
%feature("docstring") MLAPI::TimeObject "

Class to track time spent in an object.

Marzio Sala, SNL 9214

C++ includes: MLAPI_TimeObject.h ";

%feature("docstring")  MLAPI::TimeObject::TimeObject "MLAPI::TimeObject::TimeObject()

Constructor, set counter to 0.0. ";

%feature("docstring")  MLAPI::TimeObject::~TimeObject "MLAPI::TimeObject::~TimeObject()

Destructor. ";

%feature("docstring")  MLAPI::TimeObject::ResetTimer "void
MLAPI::TimeObject::ResetTimer() const

Resets the internal timer. ";

%feature("docstring")  MLAPI::TimeObject::UpdateTime "void
MLAPI::TimeObject::UpdateTime() const

Updates the internal timer with the time spent since the last call to
ResetTimer(). ";

%feature("docstring")  MLAPI::TimeObject::UpdateTime "void
MLAPI::TimeObject::UpdateTime(double t) const

Updates the internal timer with input value t. ";

%feature("docstring")  MLAPI::TimeObject::GetTime "double
MLAPI::TimeObject::GetTime() const

Returns the internally stored counter. ";


// File: namespaceMLAPI.xml
%feature("docstring")  MLAPI::GetPtent "void MLAPI::GetPtent(const
Operator &A, Teuchos::ParameterList &List, const MultiVector &ThisNS,
Operator &Ptent, MultiVector &NextNS)

Builds the tentative prolongator using aggregation. ";

%feature("docstring")  MLAPI::GetPtent "void MLAPI::GetPtent(const
Operator &A, Teuchos::ParameterList &List, Operator &Ptent)

Builds the tentative prolongator with default null space. ";

%feature("docstring")  MLAPI::GetPtent "void MLAPI::GetPtent(const
Epetra_RowMatrix &A, Teuchos::ParameterList &List, double *thisns,
Teuchos::RCP< Epetra_CrsMatrix > &Ptent, Teuchos::RCP<
Epetra_MultiVector > &NextNS, const int domainoffset=0)

Builds the tentative prolongator using aggregation.

Build Ptent and NextNS as usual but as Epetra objects.

Parameters:
-----------

A:  (in): Matrix to be aggregated on

List:  (in): ParameterList containing ML options

thisns:  (in): nullspace in format ML accepts

Ptent(out):  ::  Matrix containing tentative prolongator

NextNS:  (out): MultiVector containing next level nullspace.

domainoffset:  (in,optional): give an offset such that the domainmap
of Ptent starts global numbering from domainoffset instead from zero.
This is useful to create block matrices.

Michael Gee (gee@lnm.mw.tum.de) ";

%feature("docstring")  MLAPI::GetPtent "void MLAPI::GetPtent(const
Epetra_RowMatrix &A, Teuchos::ParameterList &List, double *thisns,
Teuchos::RCP< Epetra_CrsMatrix > &Ptent, const int domainoffset=0)

Builds the tentative prolongator using aggregation.

Build Ptent and NextNS as usual but as Epetra objects.

Parameters:
-----------

A:  (in): Matrix to be aggregated on

List:  (in): ParameterList containing ML options

thisns:  (in): nullspace in format ML accepts

Ptent(out):  ::  Matrix containing tentative prolongator

domainoffset:  (in,optional): give an offset such that the domainmap
of Ptent starts global numbering from domainoffset instead from zero.
This is useful to create block matrices.

Michael Gee (gee@lnm.mw.tum.de) ";

%feature("docstring")  MLAPI::GetAggregates "int
MLAPI::GetAggregates(Epetra_RowMatrix &A, Teuchos::ParameterList
&List, double *thisns, Epetra_IntVector &aggrinfo)

Call ML aggregation on A according to parameters supplied in List.
Return aggregates in aggrinfo.

On input, map of aggrinfo has to map row map of A. On output,
aggrinfo[i] contains number of aggregate the row belongs to, where
aggregates are numbered starting from 0. Return value is the
processor-local number of aggregates build. If aggrinfo[i] >= return-
value, then i is a processor local row of a row that ML has detected
to be on a Dirichlet BC.

Parameters:
-----------

A:  (in): Matrix to be aggregated on

List:  (in): ParameterList containing ML options

thisns:  (in): nullspace

aggrinfo(out):  ::  vector containing aggregation information

Map of aggrinfo has to match rowmap of A on input.

returns processor-local number of aggregates

Michael Gee (gee@lnm.mw.tum.de) ";

%feature("docstring")  MLAPI::GetGlobalAggregates "int
MLAPI::GetGlobalAggregates(Epetra_RowMatrix &A, Teuchos::ParameterList
&List, double *thisns, Epetra_IntVector &aggrinfo)

Call ML aggregation on A according to parameters supplied in List.
Return aggregates in aggrinfo.

On input, map of aggrinfo has to map row map of A. On output,
aggrinfo[i] contains number of global aggregate the row belongs to,
where aggregates are numbered starting from 0 globally. Return value
is the processor-local number of aggregates build. If aggrinfo[i] < 0,
then i is a processor local row that ML has detected to be on a
Dirichlet BC. if aggrinfo[i] >= 0, then i is a processor local row and
aggrinfo[i] is a global aggregate id.

Parameters:
-----------

A:  (in): Matrix to be aggregated on

List:  (in): ParameterList containing ML options

thisns:  (in): nullspace

aggrinfo(out):  ::  vector containing aggregation information in
global numbering

Map of aggrinfo has to match rowmap of A on input.

returns processor-local number of aggregates

Michael Gee (gee@lnm.mw.tum.de) ";

%feature("docstring")  MLAPI::SetDefaults "void
MLAPI::SetDefaults(Teuchos::ParameterList &List)

Sets default values in input List. ";

%feature("docstring")  MLAPI::MaxEigAnorm "double
MLAPI::MaxEigAnorm(const Operator &Op, const bool
DiagonalScaling=false)

Computes the maximum eigenvalue of Op using the A-norm of the
operator. ";

%feature("docstring")  MLAPI::MaxEigCG "double MLAPI::MaxEigCG(const
Operator &Op, const bool DiagonalScaling=false)

Computes the maximum eigenvalue of Op using the CG method. ";

%feature("docstring")  MLAPI::MaxEigPowerMethod "double
MLAPI::MaxEigPowerMethod(const Operator &Op, const bool
DiagonalScaling=false)

Computes the maximum eigenvalue of Op using the power method. ";

%feature("docstring")  MLAPI::MaxEigAnasazi "double
MLAPI::MaxEigAnasazi(const Operator &Op, const bool
DiagonalScaling=false)

Computes the maximum eigenvalue of Op using Anasazi. ";

%feature("docstring")  MLAPI::Eig "void MLAPI::Eig(const Operator
&Op, MultiVector &ER, MultiVector &EI)

Computes eigenvalues and eigenvectors using LAPACK (w/ one process
only). ";

%feature("docstring")  MLAPI::Eigs "void MLAPI::Eigs(const Operator
&A, int NumEigenvalues, MultiVector &ER, MultiVector &EI) ";

%feature("docstring")  MLAPI::StackPop "void MLAPI::StackPop() ";

%feature("docstring")  MLAPI::StackPrint "void MLAPI::StackPrint() ";

%feature("docstring")  MLAPI::Gallery "Operator MLAPI::Gallery(const
string ProblemType, const Space &MySpace)

Creates a matrix using the TRIUTILS gallery. ";

%feature("docstring")  MLAPI::GetShiftedLaplacian1D "Operator
MLAPI::GetShiftedLaplacian1D(const int NX, const double Factor=0.99)

Creates a 1D shifted Laplacian. ";

%feature("docstring")  MLAPI::GetShiftedLaplacian2D "Operator
MLAPI::GetShiftedLaplacian2D(const int NX, const int NY, const double
Factor=0.99, const bool RandomScale=false)

Creates a 2D shifted Laplacian. ";

%feature("docstring")  MLAPI::ReadMatrix "Operator
MLAPI::ReadMatrix(const char *FileName)

Reads a matrix in MATLAB format. ";

%feature("docstring")  MLAPI::GetRecirc2D "Operator
MLAPI::GetRecirc2D(const int NX, const int NY, const double conv,
const double diff)

Creates a recirculation problem in 2D. ";

%feature("docstring")  MLAPI::ReadParameterList "Teuchos::ParameterList MLAPI::ReadParameterList(const char *FileName)

Populates a list from specified file. ";

%feature("docstring")  MLAPI::Krylov "void MLAPI::Krylov(const
Operator &A, const MultiVector &LHS, const MultiVector &RHS, const
BaseOperator &Prec, Teuchos::ParameterList &List) ";

%feature("docstring")  MLAPI::Duplicate "MultiVector
MLAPI::Duplicate(const MultiVector &y)

Creates a new vector, x, such that x = y. ";

%feature("docstring")  MLAPI::Duplicate "MultiVector
MLAPI::Duplicate(const MultiVector &y, const int v)

Creates a new vector, x, such that x = y(:,v) ";

%feature("docstring")  MLAPI::Extract "MultiVector
MLAPI::Extract(const MultiVector &y, const int v)

Extracts a component from a vector. ";

%feature("docstring")  MLAPI::Redistribute "MultiVector
MLAPI::Redistribute(const MultiVector &y, const int NumEquations)

Redistributes the entry of a vector as a multivector. ";

%feature("docstring")  MLAPI::GetRAP "Operator MLAPI::GetRAP(const
Operator &R, const Operator &A, const Operator &P)

Performs a triple matrix-matrix product, res = R * A *P. ";

%feature("docstring")  MLAPI::GetTranspose "Operator
MLAPI::GetTranspose(const Operator &A, const bool byrow=true)

Returns a newly created transpose of A. ";

%feature("docstring")  MLAPI::GetIdentity "Operator
MLAPI::GetIdentity(const Space &DomainSpace, const Space &RangeSpace)

Returns the identity matrix. ";

%feature("docstring")  MLAPI::GetDiagonal "MultiVector
MLAPI::GetDiagonal(const Operator &A)

Returns a vector containing the diagonal elements of A. ";

%feature("docstring")  MLAPI::GetDiagonal "MultiVector
MLAPI::GetDiagonal(const Operator &A, const int offset)

Returns a vector containing the diagonal elements of A. ";

%feature("docstring")  MLAPI::GetDiagonal "Operator
MLAPI::GetDiagonal(const MultiVector &D)

Returns a newly created operator, containing D on the diagonal. ";

%feature("docstring")  MLAPI::GetJacobiIterationOperator "Operator
MLAPI::GetJacobiIterationOperator(const Operator &Amat, double
Damping)

Returns an operator defined as (I - Damping A). ";

%feature("docstring")  MLAPI::GetPtent1D "Operator
MLAPI::GetPtent1D(const MultiVector &D, const int offset=0)

Returns a newly created operator, containing D on the diagonal. ";

%feature("docstring")  MLAPI::ML_Operator_Add2 "int
MLAPI::ML_Operator_Add2(ML_Operator *A, ML_Operator *B, ML_Operator
*C, int matrix_type, double scalarA, double scalarB) ";

%feature("docstring")  MLAPI::AnalyzeCheap "void
MLAPI::AnalyzeCheap(const Operator &A)

Performs a cheap analysis of the properties of the input operator. ";

%feature("docstring")  MLAPI::PrintSparsity "void
MLAPI::PrintSparsity(const Operator &A, int NumPDEEquations=1)

Prints on file the sparsity structure of input operator. ";

%feature("docstring")  MLAPI::GetScaledOperator "Operator
MLAPI::GetScaledOperator(const Operator &A, const double alpha)

Multiply A by a double value, alpha. ";

%feature("docstring")  MLAPI::Duplicate "Operator
MLAPI::Duplicate(const Operator &A)

Duplicates a given operator. ";

%feature("docstring")  MLAPI::ReadSAMISMatrix "void
MLAPI::ReadSAMISMatrix(const char *filen, Operator &A, int
&NumPDEEqns)

Reads symmetric matrix from SAMIS binary format. ";

%feature("docstring")  MLAPI::ReadSAMISKernel "void
MLAPI::ReadSAMISKernel(const char *myKerFileName, MultiVector &A,
const int limKer=-1)

Reads null space vectors from SAMIS binary format. ";

%feature("docstring")  MLAPI::GetML_Comm "ML_Comm*
MLAPI::GetML_Comm()

Returns a pointer to the ML_Comm object defined on MPI_COMM_WORLD. ";

%feature("docstring")  MLAPI::GetEpetra_Comm "Epetra_Comm&
MLAPI::GetEpetra_Comm()

Returns a reference to the Epetra_Comm object defined on
MPI_COMM_WORLD. ";

%feature("docstring")  MLAPI::Barrier "void MLAPI::Barrier()

Calls Mpi_Barrier() if MPI is enabled. ";

%feature("docstring")  MLAPI::GetMyPID "int MLAPI::GetMyPID()

Returns the ID of the calling process. ";

%feature("docstring")  MLAPI::GetNumProcs "int MLAPI::GetNumProcs()

Returns the total number of processes in the computation. ";

%feature("docstring")  MLAPI::GetPrintLevel "int
MLAPI::GetPrintLevel()

Retutns the level of output (always 0 if MyPID() != 0). ";

%feature("docstring")  MLAPI::SetPrintLevel "void
MLAPI::SetPrintLevel(int Level)

Sets the level of output (from 0 to 10, 0 being verbose). ";

%feature("docstring")  MLAPI::Init "void MLAPI::Init()

Initialize the MLAPI workspace. ";

%feature("docstring")  MLAPI::Finalize "void MLAPI::Finalize()

Destroys the MLAPI workspace. ";

%feature("docstring")  MLAPI::GetString "string
MLAPI::GetString(const int &x) ";

%feature("docstring")  MLAPI::GetString "string
MLAPI::GetString(const double &x) ";

%feature("docstring")  MLAPI::GetMatrixType "int
MLAPI::GetMatrixType() ";


// File: namespaceTeuchos.xml


// File: ml__MultiLevelPreconditioner_8cpp.xml


// File: ml__MultiLevelPreconditioner_8h.xml


// File: MLAPI_8h.xml


// File: MLAPI__Aggregation_8cpp.xml


// File: MLAPI__Aggregation_8h.xml


// File: MLAPI__BaseLinearCombination_8h.xml


// File: MLAPI__BaseObject_8cpp.xml


// File: MLAPI__BaseObject_8h.xml


// File: MLAPI__BaseOperator_8h.xml


// File: MLAPI__CompObject_8h.xml


// File: MLAPI__Defaults_8cpp.xml


// File: MLAPI__Defaults_8h.xml


// File: MLAPI__DistributedMatrix_8cpp.xml


// File: MLAPI__DistributedMatrix_8h.xml


// File: MLAPI__Eig_8cpp.xml


// File: MLAPI__Eig_8h.xml


// File: MLAPI__EpetraBaseOperator_8h.xml


// File: MLAPI__Error_8cpp.xml


// File: MLAPI__Error_8h.xml


// File: MLAPI__Expressions_8cpp.xml


// File: MLAPI__Expressions_8h.xml


// File: MLAPI__Gallery_8cpp.xml


// File: MLAPI__Gallery_8h.xml


// File: MLAPI__InverseOperator_8cpp.xml


// File: MLAPI__InverseOperator_8h.xml


// File: MLAPI__Krylov_8cpp.xml


// File: MLAPI__Krylov_8h.xml


// File: MLAPI__LinearCombinations_8cpp.xml


// File: MLAPI__LinearCombinations_8h.xml


// File: MLAPI__LoadBalanceInverseOperator_8cpp.xml


// File: MLAPI__LoadBalanceInverseOperator_8h.xml


// File: MLAPI__LoadBalanceOperator_8h.xml


// File: MLAPI__MATLABStream_8h.xml


// File: MLAPI__MultiLevelAdaptiveSA_8h.xml


// File: MLAPI__MultiLevelSA_8h.xml


// File: MLAPI__MultiVector_8h.xml


// File: MLAPI__MultiVector__Utils_8cpp.xml


// File: MLAPI__MultiVector__Utils_8h.xml


// File: MLAPI__Operator_8h.xml


// File: MLAPI__Operator__Box_8h.xml


// File: MLAPI__Operator__Utils_8cpp.xml


// File: MLAPI__Operator__Utils_8h.xml


// File: MLAPI__SAMIS_8cpp.xml


// File: MLAPI__SAMIS_8h.xml


// File: MLAPI__SerialMatrix_8h.xml


// File: MLAPI__Space_8h.xml


// File: MLAPI__TimeObject_8h.xml


// File: MLAPI__Workspace_8cpp.xml


// File: MLAPI__Workspace_8h.xml


// File: dir_155fcb0d4cda658caac02c78710535e4.xml


// File: dir_cc1e09988ecec1a328546f0af31f254d.xml


// File: dir_b81d57d7b52d20c51be34aee0631351e.xml


// File: dir_3e2c3fa1e52511a70276f9374dc35ced.xml


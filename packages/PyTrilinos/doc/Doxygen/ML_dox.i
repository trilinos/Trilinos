
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
MLAPI::BaseOperatorTimesMultiVector::GetVectorSpace() const ";

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
const ";

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
verbose=true) const ";

%feature("docstring")  MLAPI::DistributedMatrix::GetDomainSpace "Space MLAPI::DistributedMatrix::GetDomainSpace() const ";

%feature("docstring")  MLAPI::DistributedMatrix::GetRangeSpace "Space
MLAPI::DistributedMatrix::GetRangeSpace() const ";

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

C++ includes: MLAPI_InverseOperator.h ";


// File: classMLAPI_1_1LinearCombinationAdd.xml
%feature("docstring") MLAPI::LinearCombinationAdd "";

%feature("docstring")
MLAPI::LinearCombinationAdd::LinearCombinationAdd "MLAPI::LinearCombinationAdd::LinearCombinationAdd(const
BaseLinearCombination &left, const BaseLinearCombination &right) ";

%feature("docstring")  MLAPI::LinearCombinationAdd::GetVectorSpace "const Space MLAPI::LinearCombinationAdd::GetVectorSpace() const ";

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

%feature("docstring")  MLAPI::LinearCombinationMixed::GetVectorSpace "const Space MLAPI::LinearCombinationMixed::GetVectorSpace() const ";

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
";

%feature("docstring")  MLAPI::LinearCombinationScaled::Set "void
MLAPI::LinearCombinationScaled::Set(MultiVector &v) const ";

%feature("docstring")  MLAPI::LinearCombinationScaled::Update "void
MLAPI::LinearCombinationScaled::Update(MultiVector &v) const ";


// File: classMLAPI_1_1MATLABStream.xml
%feature("docstring") MLAPI::MATLABStream "

Basic stream to save in a MATLAB-compatible file MLAPI objects.

For an example of usage, see ml_blackboard_cpp

Marzio Sala, SNL 9214

C++ includes: MLAPI_MATLABStream.h ";


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
An example of usage is reported in file ml_adaptivesa . We note that
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


// File: classMultiLevelPreconditioner.xml
%feature("docstring") MultiLevelPreconditioner "

ML black-box preconditioner for Epetra_RowMatrix derived classes.

C++ includes: ml_MultiLevelPreconditioner.h ";


// File: classMLAPI_1_1MultiLevelSA.xml
%feature("docstring") MLAPI::MultiLevelSA "

Black-box multilevel smoothed aggregation preconditioner.

Marzio Sala, SNL 9214

C++ includes: MLAPI_MultiLevelSA.h ";


// File: classMLAPI_1_1MultiVector.xml
%feature("docstring") MLAPI::MultiVector "

Basic class for distributed double-precision vectors.

Marzio Sala, SNL 9214.

C++ includes: MLAPI_MultiVector.h ";

%feature("docstring")  MLAPI::MultiVector::IsAlias "bool
MLAPI::MultiVector::IsAlias(const MultiVector &rhs) const ";


// File: classMLAPI_1_1MultiVectorCombination.xml
%feature("docstring") MLAPI::MultiVectorCombination "";

%feature("docstring")
MLAPI::MultiVectorCombination::MultiVectorCombination "MLAPI::MultiVectorCombination::MultiVectorCombination(const double
alpha, const MultiVector x, const double beta, const MultiVector y) ";

%feature("docstring")  MLAPI::MultiVectorCombination::GetVectorSpace "const Space MLAPI::MultiVectorCombination::GetVectorSpace() const ";

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

%feature("docstring")  MLAPI::MultiVectorScaled::GetVectorSpace "const Space MLAPI::MultiVectorScaled::GetVectorSpace() const ";

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

Marzio Sala, SNL 9214

C++ includes: MLAPI_Operator.h ";


// File: classMLAPI_1_1Residual.xml
%feature("docstring") MLAPI::Residual "";

%feature("docstring")  MLAPI::Residual::Residual "MLAPI::Residual::Residual(double alpha, const MultiVector &b, double
beta, const BaseOperator &A, const MultiVector &x) ";

%feature("docstring")  MLAPI::Residual::GetVectorSpace "const Space
MLAPI::Residual::GetVectorSpace() const ";

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
const ";


// File: classMLAPI_1_1Space.xml
%feature("docstring") MLAPI::Space "

Specifies the number and distribution among processes of elements.

Marzio Sala, SNL 9214

C++ includes: MLAPI_Space.h ";


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

Creates a new vector, x, such that x = y(:,v). ";

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


// File: namespacestd.xml


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


// File: dir_675942c6029ac094066b3b01798a20e5.xml


// File: dir_cbc5cef1c09d94196b66e1045b0d879a.xml


// File: dir_5c2a07a4854ec895e04db044a77b08e2.xml


// File: dir_c42f661d3164c34551290526b5f2c443.xml


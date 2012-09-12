
// File: index.xml

// File: structEpetra__CrsGraphData_1_1EntriesInOneRow.xml


// File: classEpetra__BasicDirectory.xml
%feature("docstring") Epetra_BasicDirectory "

Epetra_BasicDirectory: This class allows Epetra_Map objects to
reference non-local elements.

For Epetra_BlockMap objects, a Epetra_Directory object must be created
to allow referencing of non-local elements. The Epetra_BasicDirectory
produces and contains a uniform linear Epetra_BlockMap and a ProcList_
allowing blocks of non-local elements to be accessed by dereferencing
through the Epetra_BasicDirectory.

This class currently has one constructor, taking a Epetra_BlockMap
object.

C++ includes: Epetra_BasicDirectory.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_BasicDirectory::Epetra_BasicDirectory "Epetra_BasicDirectory::Epetra_BasicDirectory(const Epetra_BlockMap
&Map)

Epetra_BasicDirectory constructor. ";

%feature("docstring")  Epetra_BasicDirectory::Epetra_BasicDirectory "Epetra_BasicDirectory::Epetra_BasicDirectory(const
Epetra_BasicDirectory &Directory)

Epetra_BasicDirectory copy constructor. ";

%feature("docstring")  Epetra_BasicDirectory::~Epetra_BasicDirectory "Epetra_BasicDirectory::~Epetra_BasicDirectory(void)

Epetra_BasicDirectory destructor. ";

/*  Query method  */

%feature("docstring")  Epetra_BasicDirectory::GetDirectoryEntries "int Epetra_BasicDirectory::GetDirectoryEntries(const Epetra_BlockMap
&Map, const int NumEntries, const int *GlobalEntries, int *Procs, int
*LocalEntries, int *EntrySizes, bool high_rank_sharing_procs=false)
const

GetDirectoryEntries : Returns proc and local id info for non-local map
entries.

Given a list of Global Entry IDs, this function returns the list of
processor IDs and local IDs on the owning processor that correspond to
the list of entries. If LocalEntries is 0, then local IDs are not
returned. If EntrySizes is nonzero, it will contain a list of
corresponding element sizes for the requested global entries.

Parameters:
-----------

In:  NumEntries - Number of Global IDs being passed in.

In:  GlobalEntries - List of Global IDs being passed in.

InOut:  Procs - User allocated array of length at least NumEntries. On
return contains list of processors owning the Global IDs in question.
If any of the GIDs is shared by more than one processor, then the
lowest- numbered processor is listed in this array, unless the
optional argument 'high_rank_sharing_procs' is given as true.

InOut:  LocalEntries - User allocated array of length at least
NumEntries. On return contains the local ID of the global on the
owning processor. If LocalEntries is zero, no local ID information is
returned.

InOut:  EntrySizes - User allocated array of length at least
NumEntries. On return contains the size of the object associated with
this global ID. If LocalEntries is zero, no size information is
returned.

In:  high_rank_sharing_procs Optional argument, defaults to true. If
any GIDs appear on multiple processors (referred to as \"sharing
procs\"), this specifies whether the lowest-rank proc or the highest-
rank proc is chosen as the \"owner\".

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicDirectory::GetDirectoryEntries "int Epetra_BasicDirectory::GetDirectoryEntries(const Epetra_BlockMap
&Map, const int NumEntries, const long long *GlobalEntries, int
*Procs, int *LocalEntries, int *EntrySizes, bool
high_rank_sharing_procs=false) const ";

%feature("docstring")  Epetra_BasicDirectory::GIDsAllUniquelyOwned "bool Epetra_BasicDirectory::GIDsAllUniquelyOwned() const

GIDsAllUniquelyOwned: returns true if all GIDs appear on just one
processor.

If any GIDs are owned by multiple processors, returns false. ";

/*  I/O Methods  */

%feature("docstring")  Epetra_BasicDirectory::Print "void
Epetra_BasicDirectory::Print(ostream &os) const

Print method. ";


// File: classEpetra__BasicRowMatrix.xml
%feature("docstring") Epetra_BasicRowMatrix "

Epetra_BasicRowMatrix: A class for simplifying the development of
Epetra_RowMatrix adapters.

The Epetra_BasicRowMatrix is an adapter class for Epetra_RowMatrix
that implements most of the Epetra_RowMatrix methods using reasonable
default implementations. The Epetra_RowMatrix class has 39 pure
virtual methods, requiring the adapter class to implement all of them.
Epetra_BasicRowMatrix has only 4 pure virtual methods that must be
implemented (See Epetra_JadMatrix for an example): ExtractMyRowCopy:
Provide a row of values and indices for a specified local row.

ExtractMyEntryView (const and non-const versions): Provide the memory
address of the ith nonzero term stored on the calling processor, along
with its corresponding local row and column index, where i goes from 0
to the NumMyNonzeros()-1. The order in which the nonzeros are
traversed is not specified and is up to the adapter implementation.

NumMyRowEntries: Provide the number of entries for a specified local
row.

An alternative is possible if you do not want to provide a non-trivial
implementation of the ExtraMyEntryView methods (See
Epetra_VbrRowMatrix for an example): Implement ExtractMyRowCopy and
NumMyRowEntries as above.

Implement ExtractMyEntryView (both versions) returning a -1 integer
code with no other executable code.

Implement the RightScale and LeftScale methods non-trivially.

In addition, most adapters will probably re-implement the Multiply()
method and perhaps the Solve() method, although one or the other may
be implemented to return -1, signaling that there is no valid
implementation. By default, the Multiply() method is implemented using
ExtractMyRowCopy, which can usual be improved upon. By default Solve()
and ApplyInverse() are implemented to return -1 (not implemented).

All other implemented methods in Epetra_BasicRowMatrix should not
exhibit a signficant performance degradation, either because they are
relatively small and fast, or because they are not a significant
portion of the runtime for most codes. All methods are virtual, so
they can be re-implemented by the adapter.

In addition to implementing the above methods, an adapter must inherit
the Epetra_BasicRowMatrix interface and call the Epetra_BasicRowMatrix
constructor as part of the adapter constructor. There are two
constructors. The first requires the user to pass in the RowMap and
ColMap, both of which are Epetra_Map objects. On each processor the
RowMap (ColMap) must contain the global IDs (GIDs) of the rows
(columns) that the processor cares about. The first constructor
requires only these two maps, assuming that the RowMap will also serve
as the DomainMap and RangeMap. In this case, the RowMap must be
1-to-1, meaning that if a global ID appears on one processor, it
appears only once on that processor and does not appear on any other
processor. For many sparse matrix data structures, it is the case that
a given row is completely owned by one processor and that the global
matrix is square. The first constructor is for this situation.

The second constructor allows the caller to specify all four maps. In
this case the DomainMap, the layout of multivectors/vectors that are
in the domain of the matrix (the x vector if computing y = A*x), must
be 1-to-1. Also, the RangeMap, the layout of y must be 1-to-1. The
RowMap and ColMap do not need to be 1-to-1, but the GIDs must be found
in the RangeMap and DomainMap, respectively.

Note that Epetra_Operator is a base class for Epetra_RowMatrix, so any
adapter for Epetra_BasicRowMatrix (or Epetra_RowMatrix) is also an
adapter for Epetra_Operator.

An example of how to provide an adapter for Epetra_BasicRowMatrix can
be found by looking at Epetra_JadMatrix.

C++ includes: Epetra_BasicRowMatrix.h ";

/*  Constructor/Destructor  */

%feature("docstring")  Epetra_BasicRowMatrix::Epetra_BasicRowMatrix "Epetra_BasicRowMatrix::Epetra_BasicRowMatrix(const Epetra_Comm &Comm)

Epetra_BasicRowMatrix constuctor. ";

%feature("docstring")  Epetra_BasicRowMatrix::~Epetra_BasicRowMatrix "Epetra_BasicRowMatrix::~Epetra_BasicRowMatrix()

Epetra_BasicRowMatrix Destructor. ";

/*  Setup functions  */

%feature("docstring")  Epetra_BasicRowMatrix::SetMaps "void
Epetra_BasicRowMatrix::SetMaps(const Epetra_Map &RowMap, const
Epetra_Map &ColMap)

Set maps (Version 1); call this function or the next, but not both. ";

%feature("docstring")  Epetra_BasicRowMatrix::SetMaps "void
Epetra_BasicRowMatrix::SetMaps(const Epetra_Map &RowMap, const
Epetra_Map &ColMap, const Epetra_Map &DomainMap, const Epetra_Map
&RangeMap)

Set maps (Version 2); call this function or the previous, but not
both. ";

/*  User-required implementation methods  */

%feature("docstring")  Epetra_BasicRowMatrix::ExtractMyRowCopy "virtual int Epetra_BasicRowMatrix::ExtractMyRowCopy(int MyRow, int
Length, int &NumEntries, double *Values, int *Indices) const =0

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

MyRow:  (In) - Local row to extract.

Length:  (In) - Length of Values and Indices.

NumEntries:  (Out) - Number of nonzero entries extracted.

Values:  (Out) - Extracted values for this row.

Indices:  (Out) - Extracted global column indices for the
corresponding values.

Integer error code, set to 0 if successful, set to -1 if MyRow not
valid, -2 if Length is too short (NumEntries will have required
length). ";

%feature("docstring")  Epetra_BasicRowMatrix::ExtractMyEntryView "virtual int Epetra_BasicRowMatrix::ExtractMyEntryView(int CurEntry,
double *&Value, int &RowIndex, int &ColIndex)=0

Returns a reference to the ith entry in the matrix, along with its row
and column index.

Parameters:
-----------

CurEntry:  (In) - Index of local entry (from 0 to NumMyNonzeros()-1)
to extract.

Value:  (Out) - Extracted reference to current values.

RowIndex:  (Out) - Row index for current entry.

ColIndex:  (Out) - Column index for current entry.

Integer error code, set to 0 if successful, set to -1 if CurEntry not
valid. ";

%feature("docstring")  Epetra_BasicRowMatrix::ExtractMyEntryView "virtual int Epetra_BasicRowMatrix::ExtractMyEntryView(int CurEntry,
double const *&Value, int &RowIndex, int &ColIndex) const =0

Returns a const reference to the ith entry in the matrix, along with
its row and column index.

Parameters:
-----------

CurEntry:  (In) - Index of local entry (from 0 to NumMyNonzeros()-1)
to extract.

Value:  (Out) - Extracted reference to current values.

RowIndex:  (Out) - Row index for current entry.

ColIndex:  (Out) - Column index for current entry.

Integer error code, set to 0 if successful, set to -1 if CurEntry not
valid. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumMyRowEntries "virtual int Epetra_BasicRowMatrix::NumMyRowEntries(int MyRow, int
&NumEntries) const =0

Return the current number of values stored for the specified local
row.

Similar to NumMyEntries() except NumEntries is returned as an argument
and error checking is done on the input value MyRow.

Parameters:
-----------

MyRow:  (In) - Local row.

NumEntries:  (Out) - Number of nonzero values.

Integer error code, set to 0 if successful, set to -1 if MyRow not
valid. ";

/*  Computational methods  */

%feature("docstring")  Epetra_BasicRowMatrix::Multiply "int
Epetra_BasicRowMatrix::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_BasicRowMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

TransA:  (In) - If true, multiply by the transpose of matrix,
otherwise just use matrix.

X:  (Out) - An Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  (Out) - An Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::Solve "virtual int
Epetra_BasicRowMatrix::Solve(bool Upper, bool Trans, bool
UnitDiagonal, const Epetra_MultiVector &X, Epetra_MultiVector &Y)
const

Returns the result of a Epetra_BasicRowMatrix solve with a
Epetra_MultiVector X in Y (not implemented).

Parameters:
-----------

Upper:  (In) - If true, solve Ux = y, otherwise solve Lx = y.

Trans:  (In) - If true, solve transpose problem.

UnitDiagonal:  (In) - If true, assume diagonal is unit (whether it's
stored or not).

X:  (In) - An Epetra_MultiVector of dimension NumVectors to solve for.

Y:  (Out) - An Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::ExtractDiagonalCopy "int Epetra_BasicRowMatrix::ExtractDiagonalCopy(Epetra_Vector
&Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Diagonal:  (Out) - Extracted main diagonal.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::InvRowSums "int
Epetra_BasicRowMatrix::InvRowSums(Epetra_Vector &x) const

Computes the sum of absolute values of the rows of the
Epetra_BasicRowMatrix, results returned in x.

The vector x will return such that x[i] will contain the inverse of
sum of the absolute values of the this matrix will be scaled such that
A(i,j) = x(i)*A(i,j) where i denotes the global row number of A and j
denotes the global column number of A. Using the resulting vector from
this function as input to LeftScale() will make the infinity norm of
the resulting matrix exactly 1.

Parameters:
-----------

x:  (Out) - An Epetra_Vector containing the row sums of the this
matrix.

WARNING:  It is assumed that the distribution of x is the same as the
rows of this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::LeftScale "int
Epetra_BasicRowMatrix::LeftScale(const Epetra_Vector &x)

Scales the Epetra_BasicRowMatrix on the left with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i
denotes the row number of A and j denotes the column number of A.

Parameters:
-----------

x:  (In) - An Epetra_Vector to solve for.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::InvColSums "int
Epetra_BasicRowMatrix::InvColSums(Epetra_Vector &x) const

Computes the sum of absolute values of the columns of the
Epetra_BasicRowMatrix, results returned in x.

The vector x will return such that x[j] will contain the inverse of
sum of the absolute values of the this matrix will be sca such that
A(i,j) = x(j)*A(i,j) where i denotes the global row number of A and j
denotes the global column number of A. Using the resulting vector from
this function as input to RighttScale() will make the one norm of the
resulting matrix exactly 1.

Parameters:
-----------

x:  (Out) - An Epetra_Vector containing the column sums of the this
matrix.

WARNING:  It is assumed that the distribution of x is the same as the
rows of this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::RightScale "int
Epetra_BasicRowMatrix::RightScale(const Epetra_Vector &x)

Scales the Epetra_BasicRowMatrix on the right with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i
denotes the global row number of A and j denotes the global column
number of A.

Parameters:
-----------

x:  (In) - The Epetra_Vector used for scaling this.

Integer error code, set to 0 if successful. ";

/*  Matrix Properties Query Methods  */

%feature("docstring")  Epetra_BasicRowMatrix::Filled "virtual bool
Epetra_BasicRowMatrix::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false, presently always returns true. ";

%feature("docstring")  Epetra_BasicRowMatrix::LowerTriangular "bool
Epetra_BasicRowMatrix::LowerTriangular() const

If matrix is lower triangular, this query returns true, otherwise it
returns false. ";

%feature("docstring")  Epetra_BasicRowMatrix::UpperTriangular "virtual bool Epetra_BasicRowMatrix::UpperTriangular() const

If matrix is upper triangular, this query returns true, otherwise it
returns false. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_BasicRowMatrix::NormInf "virtual double
Epetra_BasicRowMatrix::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Epetra_BasicRowMatrix::NormOne "virtual double
Epetra_BasicRowMatrix::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalNonzeros "virtual int Epetra_BasicRowMatrix::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalNonzeros64 "virtual long long Epetra_BasicRowMatrix::NumGlobalNonzeros64() const
";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalRows "virtual
int Epetra_BasicRowMatrix::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalRows64 "virtual long long Epetra_BasicRowMatrix::NumGlobalRows64() const ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalCols "virtual
int Epetra_BasicRowMatrix::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalCols64 "virtual long long Epetra_BasicRowMatrix::NumGlobalCols64() const ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalDiagonals "virtual int Epetra_BasicRowMatrix::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumGlobalDiagonals64 "virtual long long Epetra_BasicRowMatrix::NumGlobalDiagonals64() const
";

%feature("docstring")  Epetra_BasicRowMatrix::NumMyNonzeros "virtual
int Epetra_BasicRowMatrix::NumMyNonzeros() const

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumMyRows "virtual int
Epetra_BasicRowMatrix::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Epetra_BasicRowMatrix::NumMyCols "virtual int
Epetra_BasicRowMatrix::NumMyCols() const

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")  Epetra_BasicRowMatrix::NumMyDiagonals "virtual
int Epetra_BasicRowMatrix::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries. ";

%feature("docstring")  Epetra_BasicRowMatrix::MaxNumEntries "virtual
int Epetra_BasicRowMatrix::MaxNumEntries() const

Returns the maximum number of nonzero entries across all rows on this
processor. ";

%feature("docstring")  Epetra_BasicRowMatrix::OperatorDomainMap "virtual const Epetra_Map& Epetra_BasicRowMatrix::OperatorDomainMap()
const

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Epetra_BasicRowMatrix::OperatorRangeMap "virtual const Epetra_Map& Epetra_BasicRowMatrix::OperatorRangeMap()
const

Returns the Epetra_Map object associated with the range of this
operator (same as domain). ";

%feature("docstring")  Epetra_BasicRowMatrix::Map "virtual const
Epetra_BlockMap& Epetra_BasicRowMatrix::Map() const

Implement the Epetra_SrcDistObjec::Map() function. ";

%feature("docstring")  Epetra_BasicRowMatrix::RowMatrixRowMap "virtual const Epetra_Map& Epetra_BasicRowMatrix::RowMatrixRowMap()
const

Returns the Row Map object needed for implementing Epetra_RowMatrix.
";

%feature("docstring")  Epetra_BasicRowMatrix::RowMatrixColMap "virtual const Epetra_Map& Epetra_BasicRowMatrix::RowMatrixColMap()
const

Returns the Column Map object needed for implementing
Epetra_RowMatrix. ";

%feature("docstring")  Epetra_BasicRowMatrix::RowMatrixImporter "virtual const Epetra_Import*
Epetra_BasicRowMatrix::RowMatrixImporter() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")  Epetra_BasicRowMatrix::Comm "virtual const
Epetra_Comm& Epetra_BasicRowMatrix::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

/*  I/O Methods  */

%feature("docstring")  Epetra_BasicRowMatrix::Print "void
Epetra_BasicRowMatrix::Print(ostream &os) const

Print method. ";

/*  Additional methods required to support the Epetra_RowMatrix
interface  */

%feature("docstring")  Epetra_BasicRowMatrix::SetUseTranspose "virtual int Epetra_BasicRowMatrix::SetUseTranspose(bool use_transpose)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

use_transpose:  (In) - If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Epetra_BasicRowMatrix::Label "virtual const
char* Epetra_BasicRowMatrix::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_BasicRowMatrix::Apply "virtual int
Epetra_BasicRowMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_RowMatrix applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  (In) - A Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  (Out) - A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_BasicRowMatrix::ApplyInverse "virtual
int Epetra_BasicRowMatrix::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_RowMatrix inverse applied to an
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  (In) - A Epetra_MultiVector of dimension NumVectors to solve for.

Y:  (Out) - A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code = -1.

WARNING:  This method is NOT supported. ";

%feature("docstring")  Epetra_BasicRowMatrix::HasNormInf "bool
Epetra_BasicRowMatrix::HasNormInf() const

Returns true because this class can compute an Inf-norm. ";

%feature("docstring")  Epetra_BasicRowMatrix::UseTranspose "virtual
bool Epetra_BasicRowMatrix::UseTranspose() const

Returns the current UseTranspose setting. ";

/*  Additional accessor methods  */

%feature("docstring")  Epetra_BasicRowMatrix::Importer "virtual const
Epetra_Import* Epetra_BasicRowMatrix::Importer() const

Returns the Epetra_Import object that contains the import operations
for distributed operations, returns zero if none.

If RowMatrixColMap!=OperatorDomainMap, then this method returns a
pointer to an Epetra_Import object that imports objects from an
OperatorDomainMap layout to a RowMatrixColMap layout. This operation
is needed for sparse matrix- vector multiplication, y = Ax, to gather
x elements for local multiplication operations.

If RowMatrixColMap==OperatorDomainMap, then the pointer will be
returned as 0.

Raw pointer to importer. This importer will be valid as long as the
Epetra_RowMatrix object is valid. ";

%feature("docstring")  Epetra_BasicRowMatrix::Exporter "virtual const
Epetra_Export* Epetra_BasicRowMatrix::Exporter() const

Returns the Epetra_Export object that contains the export operations
for distributed operations, returns zero if none.

If RowMatrixRowMap!=OperatorRangeMap, then this method returns a
pointer to an Epetra_Export object that exports objects from an
RowMatrixRowMap layout to a OperatorRangeMap layout. This operation is
needed for sparse matrix- vector multiplication, y = Ax, to scatter-
add y elements generated during local multiplication operations.

If RowMatrixRowMap==OperatorRangeMap, then the pointer will be
returned as 0. For a typical Epetra_RowMatrix object, this pointer
will be zero since it is often the case that
RowMatrixRowMap==OperatorRangeMap.

Raw pointer to exporter. This exporter will be valid as long as the
Epetra_RowMatrix object is valid. ";

/*  Post-construction modifications  */


// File: classEpetra__BLAS.xml
%feature("docstring") Epetra_BLAS "

Epetra_BLAS: The Epetra BLAS Wrapper Class.

The Epetra_BLAS class is a wrapper that encapsulates the BLAS (Basic
Linear Algebra Subprograms). The BLAS provide portable, high-
performance implementations of kernels such as dense vectoer
multiplication, dot products, dense matrix-vector multiplication and
dense matrix-matrix multiplication.

The standard BLAS interface is Fortran-specific. Unfortunately, the
interface between C++ and Fortran is not standard across all computer
platforms. The Epetra_BLAS class provides C++ wrappers for the BLAS
kernels in order to insulate the rest of Epetra from the details of
C++ to Fortran translation. A Epetra_BLAS object is essentially
nothing, but allows access to the BLAS wrapper functions.

Epetra_BLAS is a serial interface only. This is appropriate since the
standard BLAS are only specified for serial execution (or shared
memory parallel).

C++ includes: Epetra_BLAS.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_BLAS::Epetra_BLAS "Epetra_BLAS::Epetra_BLAS(void)

Epetra_BLAS Constructor.

Builds an instance of a serial BLAS object. ";

%feature("docstring")  Epetra_BLAS::Epetra_BLAS "Epetra_BLAS::Epetra_BLAS(const Epetra_BLAS &BLAS)

Epetra_BLAS Copy Constructor.

Makes an exact copy of an existing Epetra_BLAS instance. ";

%feature("docstring")  Epetra_BLAS::~Epetra_BLAS "Epetra_BLAS::~Epetra_BLAS(void)

Epetra_BLAS Destructor. ";

/*  Level 1 BLAS  */

%feature("docstring")  Epetra_BLAS::ASUM "float
Epetra_BLAS::ASUM(const int N, const float *X, const int INCX=1) const

Epetra_BLAS one norm function (SASUM). ";

%feature("docstring")  Epetra_BLAS::ASUM "double
Epetra_BLAS::ASUM(const int N, const double *X, const int INCX=1)
const

Epetra_BLAS one norm function (DASUM). ";

%feature("docstring")  Epetra_BLAS::DOT "float Epetra_BLAS::DOT(const
int N, const float *X, const float *Y, const int INCX=1, const int
INCY=1) const

Epetra_BLAS dot product function (SDOT). ";

%feature("docstring")  Epetra_BLAS::DOT "double
Epetra_BLAS::DOT(const int N, const double *X, const double *Y, const
int INCX=1, const int INCY=1) const

Epetra_BLAS dot product function (DDOT). ";

%feature("docstring")  Epetra_BLAS::NRM2 "float
Epetra_BLAS::NRM2(const int N, const float *X, const int INCX=1) const

Epetra_BLAS norm function (SNRM2). ";

%feature("docstring")  Epetra_BLAS::NRM2 "double
Epetra_BLAS::NRM2(const int N, const double *X, const int INCX=1)
const

Epetra_BLAS norm function (DNRM2). ";

%feature("docstring")  Epetra_BLAS::SCAL "void
Epetra_BLAS::SCAL(const int N, const float ALPHA, float *X, const int
INCX=1) const

Epetra_BLAS vector scale function (SSCAL) ";

%feature("docstring")  Epetra_BLAS::SCAL "void
Epetra_BLAS::SCAL(const int N, const double ALPHA, double *X, const
int INCX=1) const

Epetra_BLAS vector scale function (DSCAL) ";

%feature("docstring")  Epetra_BLAS::COPY "void
Epetra_BLAS::COPY(const int N, const float *X, float *Y, const int
INCX=1, const int INCY=1) const

Epetra_BLAS vector copy function (SCOPY) ";

%feature("docstring")  Epetra_BLAS::COPY "void
Epetra_BLAS::COPY(const int N, const double *X, double *Y, const int
INCX=1, const int INCY=1) const

Epetra_BLAS vector scale function (DCOPY) ";

%feature("docstring")  Epetra_BLAS::IAMAX "int
Epetra_BLAS::IAMAX(const int N, const float *X, const int INCX=1)
const

Epetra_BLAS arg maximum of absolute value function (ISAMAX) ";

%feature("docstring")  Epetra_BLAS::IAMAX "int
Epetra_BLAS::IAMAX(const int N, const double *X, const int INCX=1)
const

Epetra_BLAS arg maximum of absolute value function (IDAMAX) ";

%feature("docstring")  Epetra_BLAS::AXPY "void
Epetra_BLAS::AXPY(const int N, const float ALPHA, const float *X,
float *Y, const int INCX=1, const int INCY=1) const

Epetra_BLAS vector update function (SAXPY) ";

%feature("docstring")  Epetra_BLAS::AXPY "void
Epetra_BLAS::AXPY(const int N, const double ALPHA, const double *X,
double *Y, const int INCX=1, const int INCY=1) const

Epetra_BLAS vector update function (DAXPY) ";

/*  Level 2 BLAS  */

%feature("docstring")  Epetra_BLAS::GEMV "void
Epetra_BLAS::GEMV(const char TRANS, const int M, const int N, const
float ALPHA, const float *A, const int LDA, const float *X, const
float BETA, float *Y, const int INCX=1, const int INCY=1) const

Epetra_BLAS matrix-vector multiply function (SGEMV) ";

%feature("docstring")  Epetra_BLAS::GEMV "void
Epetra_BLAS::GEMV(const char TRANS, const int M, const int N, const
double ALPHA, const double *A, const int LDA, const double *X, const
double BETA, double *Y, const int INCX=1, const int INCY=1) const

Epetra_BLAS matrix-vector multiply function (DGEMV) ";

/*  Level 3 BLAS  */

%feature("docstring")  Epetra_BLAS::GEMM "void
Epetra_BLAS::GEMM(const char TRANSA, const char TRANSB, const int M,
const int N, const int K, const float ALPHA, const float *A, const int
LDA, const float *B, const int LDB, const float BETA, float *C, const
int LDC) const

Epetra_BLAS matrix-matrix multiply function (SGEMM) ";

%feature("docstring")  Epetra_BLAS::GEMM "void
Epetra_BLAS::GEMM(const char TRANSA, const char TRANSB, const int M,
const int N, const int K, const double ALPHA, const double *A, const
int LDA, const double *B, const int LDB, const double BETA, double *C,
const int LDC) const

Epetra_BLAS matrix-matrix multiply function (DGEMM) ";

%feature("docstring")  Epetra_BLAS::SYMM "void
Epetra_BLAS::SYMM(const char SIDE, const char UPLO, const int M, const
int N, const float ALPHA, const float *A, const int LDA, const float
*B, const int LDB, const float BETA, float *C, const int LDC) const

Epetra_BLAS symmetric matrix-matrix multiply function (SSYMM) ";

%feature("docstring")  Epetra_BLAS::SYMM "void
Epetra_BLAS::SYMM(const char SIDE, const char UPLO, const int M, const
int N, const double ALPHA, const double *A, const int LDA, const
double *B, const int LDB, const double BETA, double *C, const int LDC)
const

Epetra_BLAS matrix-matrix multiply function (DSYMM) ";

%feature("docstring")  Epetra_BLAS::TRMM "void
Epetra_BLAS::TRMM(const char SIDE, const char UPLO, const char TRANSA,
const char DIAG, const int M, const int N, const float ALPHA, const
float *A, const int LDA, float *B, const int LDB) const

Epetra_BLAS triangular matrix-matrix multiply function (STRMM) ";

%feature("docstring")  Epetra_BLAS::TRMM "void
Epetra_BLAS::TRMM(const char SIDE, const char UPLO, const char TRANSA,
const char DIAG, const int M, const int N, const double ALPHA, const
double *A, const int LDA, double *B, const int LDB) const

Epetra_BLAS triangular matrix-matrix multiply function (DTRMM) ";

%feature("docstring")  Epetra_BLAS::SYRK "void
Epetra_BLAS::SYRK(const char UPLO, const char TRANS, const int N,
const int K, const float ALPHA, const float *A, const int LDA, const
float BETA, float *C, const int LDC) const

Eperta_BLAS symetric rank k funtion (ssyrk) ";

%feature("docstring")  Epetra_BLAS::SYRK "void
Epetra_BLAS::SYRK(const char UPLO, const char TRANS, const int N,
const int K, const double ALPHA, const double *A, const int LDA, const
double BETA, double *C, const int LDC) const

Eperta_BLAS symetric rank k funtion (dsyrk) ";


// File: classEpetra__BlockMap.xml
%feature("docstring") Epetra_BlockMap "

Epetra_BlockMap: A class for partitioning block element vectors and
matrices.

It is often the case that multiple matrix and vector objects have an
identical distribution of elements on a parallel machine. The
Epetra_BlockMap class keeps information that describes this
distribution for matrices and vectors that have block elements. The
definition of an element can vary depending on the situation. For
vectors (and multi-vectors), an element is a span of one or more
contiguous entries. For matrices, it is a span of one or more matrix
rows. More generally, an element in the BlockMap class is an ordered
list of points. (NOTE: Points do not have global ID's.) Two additional
definitions useful in understanding the BlockMap class follow:
BlockMap - A distributed ordered list of elements.

First Point - First ordered point in an element

This class has a variety of constructors that can be separated into
two categories: Fixed element size constructors: All map elements have
an identical size. This corresponds to a block partitioning of
matrices and vectors where the element size is the same for all
elements. A common example is multiple degrees of freedom per mesh
node in finite element computations where the number of degrees of
freedom is the same for all nodes.

Variable element size constructor: Map element sizes may vary and are
individually defined via a list of element sizes. This is the most
general case and corresponds to a variable block partitioning of the
matrices and vectors. A common example is multiple degrees of freedom
per mesh node in finite element computations where the number of
degrees of freedom varies. This happens, for example, if regions have
differing material types or there are chemical reactions in the
simulation.

Epetra_BlockMap allows the storage and retrieval of the following
information. Depending on the constructor that is used, some of the
information is defined by the user and some is determined by the
constructor. Once an Epetra_BlockMap is constructed any of the
following can be obtained by calling a query function that has the
same name as the attribute, e.g. to get the value of
NumGlobalElements, you can call a function NumGlobalElements(). For
attributes that are lists, the query functions return the list values
in a user allocated array.

NumGlobalElements - The total number of elements across all
processors. If this parameter and NumMyElements are both passed in to
the constructor, one of the three cases will apply: If
NumGlobalElements = NumMyElements (and not equal to zero) the map is
defined to be a local replicated map. In this case, objects
constructed using this map will be identically replicated across all
processors in the communicator.

If NumGlobalElements = -1 and NumMyElements is passed in then
NumGlobalElements will be computed as the sum of NumMyElements across
all processors.

If neither of the above is true, NumGlobalElements will be checked
against the sum of NumMyElements across all processors. An error is
issued if the comparison is not equal.

NumMyElements - The number of elements owned by the calling processor.

MyGlobalElements - A list of length NumMyElements that contains the
global element IDs of the elements owned by the calling processor.

ElementSize - The size of elements if the size of all elements is the
same. This will be the case if the query function
ConstantElementSize() returns true. Otherwise this value will be set
to zero.

ElementSizeList - A list of the element sizes for elements owned by
the calling processor. This list is always accessible, even if the
element sizes are all one or of constant value. However, in these
cases, the ElementSizeList will not be generated unless a query for
the list is called.

IndexBase - The base integer value for indexed array references.
Typically this is 0 for C/C++ and 1 for Fortran, but it can be set to
any integer value.

Comm - The Epetra_Comm communicator. This communicator can in turn be
queried for processor rank and size information.

In addition to the information above that is passed in to or created
by the Epetra_BlockMap constructor, the following attributes are
computed and available via query to the user using the same scheme as
above, e.g., use NumGlobalPoints() to get the value of
NumGlobalPoints.

NumGlobalPoints - The total number of points across all processors.

NumMyPoints - The number of points on the calling processor.

MinAllGID - The minimum global index value across all processors.

MaxAllGID - The maximum global index value across all processors.

MinMyGID - The minimum global index value on the calling processor.

MaxMyGID - The maximum global index value on the calling processor.

MinLID - The minimum local index value on the calling processor.

MaxLID - The maximum local index value on the calling processor.

MinElementSize - The minimum element size across all processors.

MaxElementSize - The maximum element size across all processors.

The following functions allow boolean tests for certain properties.

ConstantElementSize() - Returns true if the element size for this map
is the same for all elements.

LinearMap() - Returns true if the elements are distributed linear
across processors, i.e., processor 0 gets the first n/p elements,
processor 1 gets the next n/p elements, etc. where n is the number of
elements and p is the number of processors.

DistributedGlobal() - Returns true if the element space of the map
spans more than one processor. This will be true in most cases, but
will be false on in serial and for objects that are created via the
derived Epetra_LocalMap class.

WARNING:  A Epetra_Comm object is required for all Epetra_BlockMap
constructors.  {error handling}

Most methods in Epetra_BlockMap return an integer error code. If the
error code is 0, then no error occurred. If > 0 then a warning error
occurred. If < 0 then a fatal error occurred.

Epetra_BlockMap constructors will throw an exception of an error
occurrs. These exceptions will alway be negative integer values as
follows: -1 NumGlobalElements < -1. Should be >= -1 (Should be >= 0
for first BlockMap constructor).

-2 NumMyElements < 0. Should be >= 0.

-3 ElementSize <= 0. Should be > 0.

-4 Invalid NumGlobalElements. Should equal sum of MyGlobalElements, or
set to -1 to compute automatically.

-5 Minimum global element index is less than index base.

-99 Internal Epetra_BlockMap error. Contact developer.

For robust code, Epetra_BlockMap constructor calls should be caught
using the try {...} catch {...} mechanism. For example:

try {      Epetra_BlockMap * map = new
Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm);   }
catch (int Error) {     if (Error==-1) { // handle error }     if
(Error==-2) ...

{ In the current implementation, Epetra_BlockMap is the base class
for:  Epetra_Map.

Epetra_LocalBlockMap.  }

C++ includes: Epetra_BlockMap.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int
ElementSize, int IndexBase, const Epetra_Comm &Comm)

Epetra_BlockMap constructor for a Epetra-defined uniform linear
distribution of constant size elements.

Creates a map that distributes NumGlobalElements elements evenly
across all processors in the Epetra_Comm communicator. If
NumGlobalElements does not divide exactly into the number of
processors, the first processors in the communicator get one extra
element until the remainder is gone.

The elements are defined to have a constant fixed size specified by
ElementSize.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute.

In:  ElementSize - Number of points or vector entries per element.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_BlockMap object. ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(unsigned int NumGlobalElements, int
ElementSize, int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(long long NumGlobalElements, int
ElementSize, int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(unsigned long long NumGlobalElements,
int ElementSize, int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int
NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm
&Comm)

Epetra_BlockMap constructor for a user-defined linear distribution of
constant size elements.

Creates a map that puts NumMyElements on the calling processor. If
NumGlobalElements=-1, the number of global elements will be the
computed sum of NumMyElements across all processors in the Epetra_Comm
communicator.

The elements are defined to have a constant fixed size specified by
ElementSize.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute. Must be
either -1 or equal to the computed sum of NumMyElements across all
processors in the Epetra_Comm communicator.

In:  NumMyElements - Number of elements owned by the calling
processor.

In:  ElementSize - Number of points or vector entries per element.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_BlockMap object. ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(unsigned int NumGlobalElements, int
NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm
&Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(long long NumGlobalElements, int
NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm
&Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(unsigned long long NumGlobalElements,
int NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm
&Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int
NumMyElements, const int *MyGlobalElements, int ElementSize, int
IndexBase, const Epetra_Comm &Comm)

Epetra_BlockMap constructor for a user-defined arbitrary distribution
of constant size elements.

Creates a map that puts NumMyElements on the calling processor. The
indices of the elements are determined from the list MyGlobalElements.
If NumGlobalElements=-1, the number of global elements will be the
computed sum of NumMyElements across all processors in the Epetra_Comm
communicator.

The elements are defined to have a constant fixed size specified by
ElementSize.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute. Must be
either -1 or equal to the computed sum of NumMyElements across all
processors in the Epetra_Comm communicator.

In:  NumMyElements - Number of elements owned by the calling
processor.

In:  MyGlobalElements - Integer array of length NumMyElements. The ith
entry contains the global index value of the ith element on this
processor. Index values are not required to be contiguous on a
processor, or to be within the range of 0 to NumGlobalElements. As
long as the index values are consistently defined and used, any set of
NumGlobalElements distinct integer values is acceptable.

In:  ElementSize - Number of points or vector entries per element.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_BlockMap object. ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(long long NumGlobalElements, int
NumMyElements, const long long *MyGlobalElements, int ElementSize, int
IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(int NumGlobalElements, int
NumMyElements, const int *MyGlobalElements, const int
*ElementSizeList, int IndexBase, const Epetra_Comm &Comm)

Epetra_BlockMap constructor for a user-defined arbitrary distribution
of variable size elements.

Creates a map that puts NumMyElements on the calling processor. If
NumGlobalElements=-1, the number of global elements will be the
computed sum of NumMyElements across all processors in the Epetra_Comm
communicator.

The elements are defined to have a variable size defined by
ElementSizeList.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute. Must be
either -1 or equal to the computed sum of NumMyElements across all
processors in the Epetra_Comm communicator.

In:  NumMyElements - Number of elements owned by the calling
processor.

In:  MyGlobalElements - Integer array of length NumMyElements. The ith
entry contains the global index value of the ith element on this
processor. Index values are not required to be contiguous on a
processor, or to be within the range of 0 to NumGlobalElements. As
long as the index values are consistently defined and used, any set of
NumGlobalElements distinct integer values is acceptable.

In:  ElementSizeList - A list of the element sizes for elements owned
by the calling processor. The ith entry contains the element size of
the ith element on this processor.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_BlockMap object. ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(long long NumGlobalElements, int
NumMyElements, const long long *MyGlobalElements, const int
*ElementSizeList, int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_BlockMap::Epetra_BlockMap "Epetra_BlockMap::Epetra_BlockMap(const Epetra_BlockMap &map)

Epetra_BlockMap copy constructor. ";

%feature("docstring")  Epetra_BlockMap::~Epetra_BlockMap "Epetra_BlockMap::~Epetra_BlockMap(void)

Epetra_BlockMap destructor. ";

/*  Local/Global ID accessor methods  */

%feature("docstring")  Epetra_BlockMap::RemoteIDList "int
Epetra_BlockMap::RemoteIDList(int NumIDs, const int *GIDList, int
*PIDList, int *LIDList) const

Returns the processor IDs and corresponding local index value for a
given list of global indices.

For each element (GID) of a given list of global element numbers
(stored in GIDList) of length NumIDs, this function returns (in
PIDList) the ID (rank) of the processor that owns the GID for this map
and returns the local index (in LIDList) of the GID on that processor.

If a GID is present on more than one processor, the lowest rank
processor ID is used, as is the LID for that processor. If a GID is
not present on any processor, the corresponding PID will return as -1.
";

%feature("docstring")  Epetra_BlockMap::RemoteIDList "int
Epetra_BlockMap::RemoteIDList(int NumIDs, const long long *GIDList,
int *PIDList, int *LIDList) const ";

%feature("docstring")  Epetra_BlockMap::RemoteIDList "int
Epetra_BlockMap::RemoteIDList(int NumIDs, const int *GIDList, int
*PIDList, int *LIDList, int *SizeList) const

Returns the processor IDs, corresponding local index value, and
element size for a given list of global indices.

For each element (GID) of a given a list of global element numbers
(stored in GIDList) of length NumIDs, this function returns (in
PIDList) the with processor that owns the GID for this map and returns
the local index (in LIDList) of the GID on that processor. Finally it
returns the element sizes in SizeList. ";

%feature("docstring")  Epetra_BlockMap::RemoteIDList "int
Epetra_BlockMap::RemoteIDList(int NumIDs, const long long *GIDList,
int *PIDList, int *LIDList, int *SizeList) const ";

%feature("docstring")  Epetra_BlockMap::LID "int
Epetra_BlockMap::LID(int GID) const

Returns local ID of global ID, return -1 if not found on this
processor. ";

%feature("docstring")  Epetra_BlockMap::LID "int
Epetra_BlockMap::LID(unsigned int GID) const ";

%feature("docstring")  Epetra_BlockMap::LID "int
Epetra_BlockMap::LID(long long GID) const ";

%feature("docstring")  Epetra_BlockMap::LID "int
Epetra_BlockMap::LID(unsigned long long GID) const ";

%feature("docstring")  Epetra_BlockMap::GID "int
Epetra_BlockMap::GID(int LID) const

Returns global ID of local ID, return IndexBase-1 if not found on this
processor. ";

%feature("docstring")  Epetra_BlockMap::GID64 "long long
Epetra_BlockMap::GID64(int LID) const ";

%feature("docstring")  Epetra_BlockMap::FindLocalElementID "int
Epetra_BlockMap::FindLocalElementID(int PointID, int &ElementID, int
&ElementOffset) const

Returns the LID of the element that contains the given local PointID,
and the Offset of the point in that element. ";

%feature("docstring")  Epetra_BlockMap::MyGID "bool
Epetra_BlockMap::MyGID(int GID_in) const

Returns true if the GID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_BlockMap::MyGID "bool
Epetra_BlockMap::MyGID(unsigned int GID_in) const ";

%feature("docstring")  Epetra_BlockMap::MyGID "bool
Epetra_BlockMap::MyGID(long long GID_in) const ";

%feature("docstring")  Epetra_BlockMap::MyLID "bool
Epetra_BlockMap::MyLID(int LID_in) const

Returns true if the LID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_BlockMap::MinAllGID "int
Epetra_BlockMap::MinAllGID() const

Returns the minimum global ID across the entire map. ";

%feature("docstring")  Epetra_BlockMap::MinAllGID64 "long long
Epetra_BlockMap::MinAllGID64() const ";

%feature("docstring")  Epetra_BlockMap::MaxAllGID "int
Epetra_BlockMap::MaxAllGID() const

Returns the maximum global ID across the entire map. ";

%feature("docstring")  Epetra_BlockMap::MaxAllGID64 "long long
Epetra_BlockMap::MaxAllGID64() const ";

%feature("docstring")  Epetra_BlockMap::MinMyGID "int
Epetra_BlockMap::MinMyGID() const

Returns the maximum global ID owned by this processor. ";

%feature("docstring")  Epetra_BlockMap::MinMyGID64 "long long
Epetra_BlockMap::MinMyGID64() const ";

%feature("docstring")  Epetra_BlockMap::MaxMyGID "int
Epetra_BlockMap::MaxMyGID() const

Returns the maximum global ID owned by this processor. ";

%feature("docstring")  Epetra_BlockMap::MaxMyGID64 "long long
Epetra_BlockMap::MaxMyGID64() const ";

%feature("docstring")  Epetra_BlockMap::MinLID "int
Epetra_BlockMap::MinLID() const

The minimum local index value on the calling processor. ";

%feature("docstring")  Epetra_BlockMap::MaxLID "int
Epetra_BlockMap::MaxLID() const

The maximum local index value on the calling processor. ";

/*  Size and dimension accessor functions  */

%feature("docstring")  Epetra_BlockMap::NumGlobalElements "int
Epetra_BlockMap::NumGlobalElements() const

Number of elements across all processors. ";

%feature("docstring")  Epetra_BlockMap::NumGlobalElements64 "long
long Epetra_BlockMap::NumGlobalElements64() const ";

%feature("docstring")  Epetra_BlockMap::NumMyElements "int
Epetra_BlockMap::NumMyElements() const

Number of elements on the calling processor. ";

%feature("docstring")  Epetra_BlockMap::MyGlobalElements "int
Epetra_BlockMap::MyGlobalElements(int *MyGlobalElementList) const

Puts list of global elements on this processor into the user-provided
array. ";

%feature("docstring")  Epetra_BlockMap::MyGlobalElements "int
Epetra_BlockMap::MyGlobalElements(long long *MyGlobalElementList)
const ";

%feature("docstring")  Epetra_BlockMap::MyGlobalElementsPtr "int
Epetra_BlockMap::MyGlobalElementsPtr(int *&MyGlobalElementList) const
";

%feature("docstring")  Epetra_BlockMap::MyGlobalElementsPtr "int
Epetra_BlockMap::MyGlobalElementsPtr(long long *&MyGlobalElementList)
const ";

%feature("docstring")  Epetra_BlockMap::ElementSize "int
Epetra_BlockMap::ElementSize() const

Returns the size of elements in the map; only valid if map has
constant element size. ";

%feature("docstring")  Epetra_BlockMap::ElementSize "int
Epetra_BlockMap::ElementSize(int LID) const

Size of element for specified LID. ";

%feature("docstring")  Epetra_BlockMap::FirstPointInElement "int
Epetra_BlockMap::FirstPointInElement(int LID) const

Returns the requested entry in the FirstPointInElementList; see
FirstPointInElementList() for details.

This function provides similar functionality to
FirstPointInElementList(), but for simple maps may avoid the explicit
construction of the FirstPointInElementList array. Returns -1 if LID
is out-of-range. ";

%feature("docstring")  Epetra_BlockMap::IndexBase "int
Epetra_BlockMap::IndexBase() const

Index base for this map. ";

%feature("docstring")  Epetra_BlockMap::NumGlobalPoints "int
Epetra_BlockMap::NumGlobalPoints() const

Number of global points for this map; equals the sum of all element
sizes across all processors. ";

%feature("docstring")  Epetra_BlockMap::NumGlobalPoints64 "long long
Epetra_BlockMap::NumGlobalPoints64() const ";

%feature("docstring")  Epetra_BlockMap::NumMyPoints "int
Epetra_BlockMap::NumMyPoints() const

Number of local points for this map; equals the sum of all element
sizes on the calling processor. ";

%feature("docstring")  Epetra_BlockMap::MinMyElementSize "int
Epetra_BlockMap::MinMyElementSize() const

Minimum element size on the calling processor. ";

%feature("docstring")  Epetra_BlockMap::MaxMyElementSize "int
Epetra_BlockMap::MaxMyElementSize() const

Maximum element size on the calling processor. ";

%feature("docstring")  Epetra_BlockMap::MinElementSize "int
Epetra_BlockMap::MinElementSize() const

Minimum element size across all processors. ";

%feature("docstring")  Epetra_BlockMap::MaxElementSize "int
Epetra_BlockMap::MaxElementSize() const

Maximum element size across all processors. ";

/*  Miscellaneous boolean tests  */

%feature("docstring")  Epetra_BlockMap::UniqueGIDs "bool
Epetra_BlockMap::UniqueGIDs() const

Returns true if map GIDs are 1-to-1.

Certain operations involving Epetra_BlockMap and Epetra_Map objects
are well-defined only if the map GIDs are uniquely present in the map.
In other words, if a GID occurs in the map, it occurs only once on a
single processor and nowhere else. This boolean test returns true if
this property is true, otherwise it returns false. ";

%feature("docstring")  Epetra_BlockMap::GlobalIndicesInt "bool
Epetra_BlockMap::GlobalIndicesInt() const

Returns true if map create with int NumGlobalElements. ";

%feature("docstring")  Epetra_BlockMap::GlobalIndicesLongLong "bool
Epetra_BlockMap::GlobalIndicesLongLong() const

Returns true if map create with long long NumGlobalElements. ";

%feature("docstring")  Epetra_BlockMap::GlobalIndicesIsType "bool
Epetra_BlockMap::GlobalIndicesIsType< long long >() const ";

%feature("docstring")  Epetra_BlockMap::GlobalIndicesTypeValid "bool
Epetra_BlockMap::GlobalIndicesTypeValid() const ";

%feature("docstring")  Epetra_BlockMap::GlobalIndicesTypeMatch "bool
Epetra_BlockMap::GlobalIndicesTypeMatch(const Epetra_BlockMap &other)
const ";

%feature("docstring")  Epetra_BlockMap::ConstantElementSize "bool
Epetra_BlockMap::ConstantElementSize() const

Returns true if map has constant element size. ";

%feature("docstring")  Epetra_BlockMap::SameAs "bool
Epetra_BlockMap::SameAs(const Epetra_BlockMap &Map) const

Returns true if this and Map are identical maps. ";

%feature("docstring")  Epetra_BlockMap::PointSameAs "bool
Epetra_BlockMap::PointSameAs(const Epetra_BlockMap &Map) const

Returns true if this and Map have identical point-wise structure.

If both maps have the same number of global points and the same point
distribution across processors then this method returns true. ";

%feature("docstring")  Epetra_BlockMap::LinearMap "bool
Epetra_BlockMap::LinearMap() const

Returns true if the global ID space is contiguously divided (but not
necessarily uniformly) across all processors. ";

%feature("docstring")  Epetra_BlockMap::DistributedGlobal "bool
Epetra_BlockMap::DistributedGlobal() const

Returns true if map is defined across more than one processor. ";

/*  Array accessor functions  */

%feature("docstring")  Epetra_BlockMap::MyGlobalElements "int *
Epetra_BlockMap::MyGlobalElements() const

Pointer to internal array containing list of global IDs assigned to
the calling processor. ";

%feature("docstring")  Epetra_BlockMap::MyGlobalElements64 "long long
* Epetra_BlockMap::MyGlobalElements64() const ";

%feature("docstring")  Epetra_BlockMap::FirstPointInElementList "int
* Epetra_BlockMap::FirstPointInElementList() const

Pointer to internal array containing a mapping between the local
elements and the first local point number in each element.

This array is a scan sum of the ElementSizeList such that the ith
entry in FirstPointInElementList is the sum of the first i-1 entries
of ElementSizeList(). ";

%feature("docstring")  Epetra_BlockMap::ElementSizeList "int *
Epetra_BlockMap::ElementSizeList() const

List of the element sizes corresponding to the array
MyGlobalElements(). ";

%feature("docstring")  Epetra_BlockMap::PointToElementList "int *
Epetra_BlockMap::PointToElementList() const

For each local point, indicates the local element ID that the point
belongs to. ";

%feature("docstring")  Epetra_BlockMap::ElementSizeList "int
Epetra_BlockMap::ElementSizeList(int *ElementSizeList) const

Same as ElementSizeList() except it fills the user array that is
passed in. ";

%feature("docstring")  Epetra_BlockMap::FirstPointInElementList "int
Epetra_BlockMap::FirstPointInElementList(int *FirstPointInElementList)
const

Same as FirstPointInElementList() except it fills the user array that
is passed in. ";

%feature("docstring")  Epetra_BlockMap::PointToElementList "int
Epetra_BlockMap::PointToElementList(int *PointToElementList) const

Same as PointToElementList() except it fills the user array that is
passed in. ";

/*  Miscellaneous  */

%feature("docstring")  Epetra_BlockMap::Print "void
Epetra_BlockMap::Print(ostream &os) const

Print object to an output stream. ";

%feature("docstring")  Epetra_BlockMap::Comm "const Epetra_Comm&
Epetra_BlockMap::Comm() const

Access function for Epetra_Comm communicator. ";

%feature("docstring")  Epetra_BlockMap::IsOneToOne "bool
Epetra_BlockMap::IsOneToOne() const ";

/*  Expert Users and Developers Only  */

%feature("docstring")  Epetra_BlockMap::ReferenceCount "int
Epetra_BlockMap::ReferenceCount() const

Returns the reference count of BlockMapData.

(Intended for testing purposes.) ";

%feature("docstring")  Epetra_BlockMap::DataPtr "const
Epetra_BlockMapData* Epetra_BlockMap::DataPtr() const

Returns a pointer to the BlockMapData instance this BlockMap uses.

(Intended for developer use only for testing purposes.) ";


// File: classEpetra__BlockMapData.xml
%feature("docstring") Epetra_BlockMapData "

Epetra_BlockMapData: The Epetra BlockMap Data Class.

The Epetra_BlockMapData class is an implementation detail of
Epetra_BlockMap. It is reference-counted, and can be shared by
multiple Epetra_BlockMap instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_BlockMapData.h ";

/*  Constructor/Destructor Methods  */


// File: classEpetra__Comm.xml
%feature("docstring") Epetra_Comm "

Epetra_Comm: The Epetra Communication Abstract Base Class.

The Epetra_Comm class is an interface that encapsulates the general
information and services needed for other Epetra classes to run on a
parallel computer. An Epetra_Comm object is required for building all
Epetra Map objects, which in turn are required for all other Epetra
classes.

Epetra_Comm has default implementations, via Epetra_SerialComm and
Epetra_MpiComm, for both serial execution and MPI distributed memory
execution. It is meant to insulate the user from the specifics of
communication that are not required for normal manipulation of linear
algebra objects. Most Epetra_Comm interfaces are similar to MPI
interfaces, except that the type of data is not required as an
argument since C++ can bind to the appropriate interface based on
argument typing.

Any implementation of the Epetra_Comm interface is also responsible
for generating an Epetra_Distributor and Epetra_Directory object.

C++ includes: Epetra_Comm.h ";

/*  Constructor / Destructor  */

%feature("docstring")  Epetra_Comm::Clone "virtual Epetra_Comm*
Epetra_Comm::Clone() const =0

Epetra_Comm clone constructor.

The clone function will return a new heap-allocated Comm instance. It
is the responsibility of the caller to ensure that this new instance
is properly destroyed. ";

%feature("docstring")  Epetra_Comm::~Epetra_Comm "virtual
Epetra_Comm::~Epetra_Comm()

Epetra_Comm Destructor. ";

/*  Barrier Methods  */

%feature("docstring")  Epetra_Comm::Barrier "virtual void
Epetra_Comm::Barrier() const =0

Epetra_Comm Barrier function.

Each processor must wait at the point the barrier is called until all
processors have arrived. ";

/*  Broadcast Methods  */

%feature("docstring")  Epetra_Comm::Broadcast "virtual int
Epetra_Comm::Broadcast(double *MyVals, int Count, int Root) const =0

Epetra_Comm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_Comm::Broadcast "virtual int
Epetra_Comm::Broadcast(int *MyVals, int Count, int Root) const =0

Epetra_Comm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_Comm::Broadcast "virtual int
Epetra_Comm::Broadcast(long *MyVals, int Count, int Root) const =0

Epetra_Comm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_Comm::Broadcast "virtual int
Epetra_Comm::Broadcast(char *MyVals, int Count, int Root) const =0

Epetra_Comm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

/*  Gather Methods  */

%feature("docstring")  Epetra_Comm::GatherAll "virtual int
Epetra_Comm::GatherAll(double *MyVals, double *AllVals, int Count)
const =0

Epetra_Comm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must be of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_Comm::GatherAll "virtual int
Epetra_Comm::GatherAll(int *MyVals, int *AllVals, int Count) const =0

Epetra_Comm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must be of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_Comm::GatherAll "virtual int
Epetra_Comm::GatherAll(long *MyVals, long *AllVals, int Count) const
=0

Epetra_Comm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must be of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_Comm::GatherAll "virtual int
Epetra_Comm::GatherAll(long long *MyVals, long long *AllVals, int
Count) const =0

Epetra_Comm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must be of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

/*  Sum Methods  */

%feature("docstring")  Epetra_Comm::SumAll "virtual int
Epetra_Comm::SumAll(double *PartialSums, double *GlobalSums, int
Count) const =0

Epetra_Comm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::SumAll "virtual int
Epetra_Comm::SumAll(int *PartialSums, int *GlobalSums, int Count)
const =0

Epetra_Comm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::SumAll "virtual int
Epetra_Comm::SumAll(long *PartialSums, long *GlobalSums, int Count)
const =0

Epetra_Comm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::SumAll "virtual int
Epetra_Comm::SumAll(long long *PartialSums, long long *GlobalSums, int
Count) const =0

Epetra_Comm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Max/Min Methods  */

%feature("docstring")  Epetra_Comm::MaxAll "virtual int
Epetra_Comm::MaxAll(double *PartialMaxs, double *GlobalMaxs, int
Count) const =0

Epetra_Comm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MaxAll "virtual int
Epetra_Comm::MaxAll(int *PartialMaxs, int *GlobalMaxs, int Count)
const =0

Epetra_Comm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MaxAll "virtual int
Epetra_Comm::MaxAll(long *PartialMaxs, long *GlobalMaxs, int Count)
const =0

Epetra_Comm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MaxAll "virtual int
Epetra_Comm::MaxAll(long long *PartialMaxs, long long *GlobalMaxs, int
Count) const =0

Epetra_Comm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MinAll "virtual int
Epetra_Comm::MinAll(double *PartialMins, double *GlobalMins, int
Count) const =0

Epetra_Comm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MinAll "virtual int
Epetra_Comm::MinAll(int *PartialMins, int *GlobalMins, int Count)
const =0

Epetra_Comm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MinAll "virtual int
Epetra_Comm::MinAll(long *PartialMins, long *GlobalMins, int Count)
const =0

Epetra_Comm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::MinAll "virtual int
Epetra_Comm::MinAll(long long *PartialMins, long long *GlobalMins, int
Count) const =0

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Parallel Prefix Methods  */

%feature("docstring")  Epetra_Comm::ScanSum "virtual int
Epetra_Comm::ScanSum(double *MyVals, double *ScanSums, int Count)
const =0

Epetra_Comm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::ScanSum "virtual int
Epetra_Comm::ScanSum(int *MyVals, int *ScanSums, int Count) const =0

Epetra_Comm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::ScanSum "virtual int
Epetra_Comm::ScanSum(long *MyVals, long *ScanSums, int Count) const =0

Epetra_Comm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_Comm::ScanSum "virtual int
Epetra_Comm::ScanSum(long long *MyVals, long long *ScanSums, int
Count) const =0

Epetra_Comm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

/*  Attribute Accessor Methods  */

%feature("docstring")  Epetra_Comm::MyPID "virtual int
Epetra_Comm::MyPID() const =0

Return my process ID.

In MPI mode returns the rank of the calling process. In serial mode
returns 0. ";

%feature("docstring")  Epetra_Comm::NumProc "virtual int
Epetra_Comm::NumProc() const =0

Returns total number of processes.

In MPI mode returns the size of the MPI communicator. In serial mode
returns 1. ";

/*  Gather/Scatter and Directory Constructors  */

%feature("docstring")  Epetra_Comm::CreateDistributor "virtual
Epetra_Distributor* Epetra_Comm::CreateDistributor() const =0

Create a distributor object. ";

%feature("docstring")  Epetra_Comm::CreateDirectory "virtual
Epetra_Directory* Epetra_Comm::CreateDirectory(const Epetra_BlockMap
&Map) const =0

Create a directory object for the given Epetra_BlockMap. ";

/*  I/O methods  */

%feature("docstring")  Epetra_Comm::PrintInfo "virtual void
Epetra_Comm::PrintInfo(ostream &os) const =0

Print object to an output stream. ";


// File: classEpetra__CompObject.xml
%feature("docstring") Epetra_CompObject "

Epetra_CompObject: Functionality and data that is common to all
computational classes.

The Epetra_CompObject is a base class for all Epetra computational
objects. It provides the basic mechanisms and interface specifications
for floating point operations using Epetra_Flops objects.

C++ includes: Epetra_CompObject.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_CompObject::Epetra_CompObject "Epetra_CompObject::Epetra_CompObject()

Basic Epetra_CompObject constuctor. ";

%feature("docstring")  Epetra_CompObject::Epetra_CompObject "Epetra_CompObject::Epetra_CompObject(const Epetra_CompObject &Source)

Epetra_CompObject copy constructor. ";

%feature("docstring")  Epetra_CompObject::~Epetra_CompObject "Epetra_CompObject::~Epetra_CompObject()

Epetra_CompObject destructor. ";

/*  Set/Get counter method  */

%feature("docstring")  Epetra_CompObject::SetFlopCounter "void
Epetra_CompObject::SetFlopCounter(const Epetra_Flops &FlopCounter_in)

Set the internal Epetra_Flops() pointer. ";

%feature("docstring")  Epetra_CompObject::SetFlopCounter "void
Epetra_CompObject::SetFlopCounter(const Epetra_CompObject &CompObject)

Set the internal Epetra_Flops() pointer to the flop counter of another
Epetra_CompObject. ";

%feature("docstring")  Epetra_CompObject::UnsetFlopCounter "void
Epetra_CompObject::UnsetFlopCounter()

Set the internal Epetra_Flops() pointer to 0 (no flops counted). ";

%feature("docstring")  Epetra_CompObject::GetFlopCounter "Epetra_Flops* Epetra_CompObject::GetFlopCounter() const

Get the pointer to the Epetra_Flops() object associated with this
object, returns 0 if none. ";

/*  Set flop count methods  */

%feature("docstring")  Epetra_CompObject::ResetFlops "void
Epetra_CompObject::ResetFlops() const

Resets the number of floating point operations to zero for this multi-
vector. ";

%feature("docstring")  Epetra_CompObject::Flops "double
Epetra_CompObject::Flops() const

Returns the number of floating point operations with this multi-
vector. ";

/*  Update flop count methods  */

%feature("docstring")  Epetra_CompObject::UpdateFlops "void
Epetra_CompObject::UpdateFlops(int Flops_in) const

Increment Flop count for this object. ";

%feature("docstring")  Epetra_CompObject::UpdateFlops "void
Epetra_CompObject::UpdateFlops(long int Flops_in) const

Increment Flop count for this object. ";

%feature("docstring")  Epetra_CompObject::UpdateFlops "void
Epetra_CompObject::UpdateFlops(long long Flops_in) const

Increment Flop count for this object. ";

%feature("docstring")  Epetra_CompObject::UpdateFlops "void
Epetra_CompObject::UpdateFlops(double Flops_in) const

Increment Flop count for this object. ";

%feature("docstring")  Epetra_CompObject::UpdateFlops "void
Epetra_CompObject::UpdateFlops(float Flops_in) const

Increment Flop count for this object. ";


// File: classEpetra__CrsGraph.xml
%feature("docstring") Epetra_CrsGraph "

Epetra_CrsGraph: A class for constructing and using sparse compressed
row graphs.

Epetra_CrsGraph enables the piecewise construction and use of sparse
matrix graphs (the integer structure without values) where entries are
intended for row access.

Epetra_CrsGraph is an attribute of all Epetra row-based matrix
classes, defining their nonzero structure and also holding their
Epetra_Map attributes.

Constructing Epetra_CrsGraph objects

Constructing Epetra_CrsGraph objects is a multi-step process. The
basic steps are as follows: Create Epetra_CrsGraph instance, including
some initial storage, via constructor. In addition to the copy
constructor, Epetra_CrsGraph has four different constructors. All four
of these constructors have an argument, StaticProfile, which by
default is set to false. If it is set to true, then the profile (the
number of indices per row as defined by NumIndicesPerRow) will be
rigidly enforced. Although this takes away flexibility, it allows a
single array to be allocated for all indices. This decreases memory
fragmentation and improves performance across many operations. A more
detailed discussion of the StaticProfile option is found below. User-
provided row map, variable nonzero profile: This constructor is used
to define the row distribution of the graph and specify a varying
number of nonzero entries per row. It is best to use this constructor
when the user will be inserting entries using global index values and
wants every column index to be included in the graph. Note that in
this case, the column map will be built for the user when
FillComplete() is called. This constructor is also appropriate for
when there is a large variation in the number of indices per row. If
this is not the case, the next constructor may be more convenient to
use.

User-provided row map, fixed nonzero profile: This constructor is used
to define the row distribution of the graph and specify a fixed number
of nonzero entries per row. It is best to use this constructor when
the user will be inserting entries using global index values and wants
every column index to be included in the graph. Note that in this
case, the column map will be built for the user when FillComplete() is
called. This constructor is also appropriate for when there is little
or no variation in the number of indices per row.

User-provided row map, user-provided column map and variable nonzero
profile: This constructor is used to define the row and column
distribution of the graph, and specify a varying number of nonzero
entries per row. It is best to use this constructor when the user will
be inserting entries and already knows which columns of the matrix
should be included on each processor. Note that in this case, the
column map will not be built for the user when FillComplete() is
called. Also, if the user attempts to insert a column index whose GID
is not part of the column map on that process, the index will be
discarded. This property can be used to \"filter out\" column entries
that should be ignored. This constructor is also appropriate for when
there is a large variation in the number of indices per row. If this
is not the case, the next constructor may be more convenient to use.

User-provided row map, user-provided column map and fixed nonzero
profile: This constructor is used to define the row and column
distribution of the graph, and specify a fixed number of nonzero
entries per row. It is best to use this constructor when the user will
be inserting entries and already knows which columns of the matrix
should be included on each processor. Note that in this case, the
column map will not be built for the user when FillComplete() is
called. Also, if the user attempts to insert a column index whose GID
is not part of the column map on that process, the index will be
discarded. This constructor is also appropriate for when there is
little or no variation in the number of indices per row.

Enter row and column entry information via calls to the
InsertGlobalIndices method.

Complete construction via FillComplete call, which performs the
following tasks: Transforms indices to local index space (after this,
IndicesAreLocal()==true)

Sorts column-indices within each row

Compresses out any redundant indices within rows

Computes global data such as num-nonzeros, maximum row-lengths, etc.

(Optional) Optimize the graph storage via a call to OptimizeStorage.

Performance Enhancement Issues

The Epetra_CrsGraph class attempts to address four basic types of
situations, depending on the user's primary concern:

Simple, flexible construction over minimal memory use or control of
column indices: In this case the user wants to provide only a row
distribution of the graph and insert indices without worrying about
memory allocation performance. This type of user is best served by the
constructor that requires only a row map, and a fixed number of
indices per row. In fact, setting NumIndicesPerRow=0 is probably the
best option.

Stronger control over memory allocation performance and use over
flexibility and simplicity: In this case the user explicitly set
StaticProfile to true and will provide values, either a single global
int or an array of int's, for NumIndicesPerRow, such that the actual
number of indices submitted to the graph will not exceed the
estimates. Because we know that NumIndicesPerRow will not be exceeded,
we can pre-allocate all of the storage for the graph as a single
array. This is typically much more efficient.

Explicit control over column indices: In this case the user prescribes
the column map. Given the column map, any index that is submitted for
entry into the graph will be included only if they are present in the
list of GIDs for the column map on the processor that submits the
index. This feature allows the user to define a filter such that only
certain columns will be kept. The user also prescribes the local
ordering via this technique, since the ordering of GIDs in the column
map imposes the local ordering.

Construction using local indices only: In some situations, users may
want to build a graph using local index values only. In this case, the
user must explicitly assign GIDs. This is done by prescribing the
column map, in the same way as the previous situation.

Notes: In all but the most advanced uses, users will typically not
specify the column map. In other words, graph entries will be
submitted using GIDs not LIDs and all entries that are submitted are
intended to be inserted into the graph.

If a user is not particularly worried about performance, or really
needs the flexibility associated with the first situation, then there
is no need to explicitly manage the NumIndicesPerRow values or set
StaticProfile to true. In this case, it is best to set
NumIndicesPerRow to zero.

Users who are concerned about performance should carefully manage
NumIndicesPerRow and set StaticProfile to true. This will give the
best performance and use the least amount of memory.

A compromise approach would be to not set StaticProfile to true,
giving the user flexibility, but then calling OptimizeStorage() once
FillComplete() has been called. This approach requires additional
temporary memory because the graph will be copied into an efficient
data structure and the old memory deleted. However, once the copy has
been made, the resulting data structure is as efficient as when
StaticProfile is used.

Epetra_Map attributes

Epetra_CrsGraph objects have four Epetra_Map attributes.

The Epetra_Map attributes can be obtained via these accessor methods:
RowMap() Describes the numbering and distribution of the rows of the
graph. The row-map exists and is valid for the entire life of the
graph, having been passed in as a constructor argument. The set of
graph rows is defined by the row-map and may not be changed. Rows may
not be inserted or deleted by the user. The only change that may be
made is that the user can replace the row-map with a compatible row-
map (which is the same except for re-numbering) by calling the
ReplaceRowMap() method.

ColMap() Describes the set of column-indices that appear in the rows
in each processor's portion of the graph. Unless provided by the user
at construction time, a valid column-map doesn't exist until
FillComplete() is called.

RangeMap() Describes the range of the matrix operator. e.g., for a
matrix-vector product operation, the result vector's map must be
compatible with the range-map of the matrix operator. The range-map is
usually the same as the row-map. The range-map is set equal to the
row-map at graph creation time, but may be specified by the user when
FillComplete() is called.

DomainMap() Describes the domain of the matrix operator. The domain-
map can be specified by the user when FillComplete() is called. Until
then, it is set equal to the row-map.

It is important to note that while the row-map and the range-map are
often the same, the column-map and the domain-map are almost never the
same. The set of entries in a distributed column-map almost always
form overlapping sets, with entries being associated with more than
one processor. A domain-map, on the other hand, must be a 1-to-1 map,
with entries being associated with only a single processor.

Global versus Local indices

After creation and before FillComplete() has been called, the column-
indices of the graph are in the global space as received from the
user. One of the tasks performed by FillComplete() is to transform the
indices to a local index space. The query methods IndicesAreGlobal()
and IndicesAreLocal() return true or false depending on whether this
transformation has been performed or not.

Note the behavior of several graph methods:  InsertGlobalIndices()
returns an error if IndicesAreLocal()==true or
StorageOptimized()==true

InsertMyIndices() returns an error if IndicesAreGlobal()==true or
StorageOptimized()==true

RemoveGlobalIndices() returns an error if IndicesAreLocal()==true or
if graph was constructed in View mode

RemoveMyIndices() returns an error if IndicesAreGlobal()==true or if
graph was constructed in View mode

ExtractGlobalRowCopy() works regardless of state of indices

ExtractMyRowCopy() returns an error if IndicesAreGlobal()==true

ExtractGlobalRowView() returns an error if IndicesAreLocal()==true

ExtractMyRowView() returns an error if IndicesAreGlobal()==true

Note that even after a graph is constructed, it is possible to add or
remove entries. However, FillComplete must then be called again to
restore the graph to a consistent state.

C++ includes: Epetra_CrsGraph.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_CrsGraph::Epetra_CrsGraph "Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const int *NumIndicesPerRow, bool
StaticProfile=false)

Epetra_CrsGraph constuctor with variable number of indices per row.

Creates a Epetra_CrsGraph object and allocates storage.

Parameters:
-----------

CV:  - (In) A Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap)
listing the rows that this processor will contribute to.In

NumIndicesPerRow:  - (In) An integer array of length NumMyRows such
that NumIndicesPerRow[i] indicates the (approximate if
StaticProfile=false) number of entries in the ith row.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsGraph::Epetra_CrsGraph "Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int NumIndicesPerRow, bool
StaticProfile=false)

Epetra_CrsGraph constuctor with fixed number of indices per row.

Creates a Epetra_CrsGraph object and allocates storage.

Parameters:
-----------

CV:  - (In) A Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap)
listing the rows that this processor will contribute to.

NumIndicesPerRow:  - (In) An integer that indicates the (approximate
if StaticProfile=false) number of entries in the each row. Note that
it is possible to use 0 for this value and let fill occur during the
insertion phase.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsGraph::Epetra_CrsGraph "Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, const int
*NumIndicesPerRow, bool StaticProfile=false)

Epetra_CrsGraph constuctor with variable number of indices per row.

Creates a Epetra_CrsGraph object and allocates storage.

Parameters:
-----------

CV:  - (In) A Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap)
listing the rows that this processor will contribute to.

ColMap:  - (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap)
listing the columns that this processor will contribute to.

NumIndicesPerRow:  - (In) An integer array of length NumMyRows such
that NumIndicesPerRow[i] indicates the (approximate if
StaticProfile=false) number of entries in the ith row.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsGraph::Epetra_CrsGraph "Epetra_CrsGraph::Epetra_CrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
NumIndicesPerRow, bool StaticProfile=false)

Epetra_CrsGraph constuctor with fixed number of indices per row.

Creates a Epetra_CrsGraph object and allocates storage.

Parameters:
-----------

CV:  - (In) A Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap)
listing the rows that this processor will contribute to.

ColMap:  - (In) An Epetra_BlockMap (or Epetra_Map or Epetra_LocalMap)
listing the columns that this processor will contribute to.

In:  NumIndicesPerRow - An integer that indicates the (approximate if
StaticProfile=false) number of entries in the each row. Note that it
is possible to use 0 for this value and let fill occur during the
insertion phase.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsGraph::Epetra_CrsGraph "Epetra_CrsGraph::Epetra_CrsGraph(const Epetra_CrsGraph &Graph)

Copy constructor.

This will create a Level 1 deep copy. This Graph will share ownership
of the CrsGraphData object with the right hand side Graph. ";

%feature("docstring")  Epetra_CrsGraph::~Epetra_CrsGraph "Epetra_CrsGraph::~Epetra_CrsGraph()

Epetra_CrsGraph Destructor. ";

/*  Insertion/Removal methods  */

%feature("docstring")  Epetra_CrsGraph::InsertGlobalIndices "int
Epetra_CrsGraph::InsertGlobalIndices(int GlobalRow, int NumIndices,
int *Indices)

Enter a list of elements in a specified global row of the graph.

Parameters:
-----------

Row:  - (In) Global row number of indices.

NumIndices:  - (In) Number of Indices.

Indices:  - (In) Global column indices to insert.

Integer error code, set to 0 if successful. If the insertion requires
that additional memory be allocated for the row, a positive error code
of 1 is returned. If the graph is a 'View' mode graph, then a positive
warning code of 2 will be returned if the specified row already
exists. Returns 1 if underlying graph data is shared by multiple graph
instances.

IndicesAreGlobal()==true, StorageOptimized()==false ";

%feature("docstring")  Epetra_CrsGraph::InsertGlobalIndices "int
Epetra_CrsGraph::InsertGlobalIndices(long long GlobalRow, int
NumIndices, long long *Indices) ";

%feature("docstring")  Epetra_CrsGraph::RemoveGlobalIndices "int
Epetra_CrsGraph::RemoveGlobalIndices(int GlobalRow, int NumIndices,
int *Indices)

Remove a list of elements from a specified global row of the graph.

Parameters:
-----------

Row:  - (In) Global row number of indices.

NumIndices:  - (In) Number of Indices.

Indices:  - (In) Global column indices to remove.

Integer error code, set to 0 if successful. Returns 1 if data is
shared.

IndicesAreGlobal()==true, StorageOptimized()==false ";

%feature("docstring")  Epetra_CrsGraph::RemoveGlobalIndices "int
Epetra_CrsGraph::RemoveGlobalIndices(long long GlobalRow, int
NumIndices, long long *Indices) ";

%feature("docstring")  Epetra_CrsGraph::RemoveGlobalIndices "int
Epetra_CrsGraph::RemoveGlobalIndices(long long Row)

Remove all indices from a specified global row of the graph.

Parameters:
-----------

Row:  - (In) Global row number of indices.

Integer error code, set to 0 if successful. Returns 1 if data is
shared.

IndicesAreGlobal()==true, StorageOptimized()==false ";

%feature("docstring")  Epetra_CrsGraph::InsertMyIndices "int
Epetra_CrsGraph::InsertMyIndices(int LocalRow, int NumIndices, int
*Indices)

Enter a list of elements in a specified local row of the graph.

Parameters:
-----------

Row:  - (In) Local row number of indices.

NumIndices:  - (In) Number of Indices.

Indices:  - (In) Local column indices to insert.

Integer error code, set to 0 if successful. If the insertion requires
that additional memory be allocated for the row, a positive error code
of 1 is returned. If one or more of the indices is ignored (due to not
being contained in the column-map), then a positive warning code of 2
is returned. If the graph is a 'View' mode graph, then a positive
warning code of 3 will be returned if the specified row already
exists. Returns 1 if underlying graph data is shared by multiple graph
instances.

IndicesAreLocal()==true, StorageOptimized()==false ";

%feature("docstring")  Epetra_CrsGraph::RemoveMyIndices "int
Epetra_CrsGraph::RemoveMyIndices(int LocalRow, int NumIndices, int
*Indices)

Remove a list of elements from a specified local row of the graph.

Parameters:
-----------

Row:  - (In) Local row number of indices.

NumIndices:  - (In) Number of Indices.

Indices:  - (In) Local column indices to remove.

Integer error code, set to 0 if successful. Returns 1 if data is
shared.

IndicesAreLocal()==true, StorageOptimized()==false ";

%feature("docstring")  Epetra_CrsGraph::RemoveMyIndices "int
Epetra_CrsGraph::RemoveMyIndices(int Row)

Remove all indices from a specified local row of the graph.

Parameters:
-----------

Row:  - (In) Local row number of indices.

Integer error code, set to 0 if successful. Returns 1 if data is
shared.

IndicesAreLocal()==true, StorageOptimized()==false ";

/*  Transformation methods  */

%feature("docstring")  Epetra_CrsGraph::FillComplete "int
Epetra_CrsGraph::FillComplete()

Tranform to local index space. Perform other operations to allow
optimal matrix operations.

This overloading of the FillComplete method assumes that the domain-
map and range-map both equal the row-map, and simply calls
FillComplete( RowMap(), RowMap()). Integer error code, set to 0 if
successful. Returns 1 if data is shared (i.e., if the underlying
graph-data object has a reference- count greater than 1).

IndicesAreLocal()==true, Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::FillComplete "int
Epetra_CrsGraph::FillComplete(const Epetra_BlockMap &DomainMap, const
Epetra_BlockMap &RangeMap)

Transform to local index space using specified Domain/Range maps.
Perform other operations to allow optimal matrix operations.

Performs this sequence of operations: Transform indices to local index
space

Sort column-indices within each row

Compress out any redundant indices within rows

Compute global data such as num-nonzeros, maximum row-lengths, etc.

Integer error code, set to 0 if successful. Returns 1 if data is
shared (i.e., if the underlying graph-data object has a reference-
count greater than 1).

IndicesAreLocal()==true, Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::OptimizeStorage "int
Epetra_CrsGraph::OptimizeStorage()

Make consecutive row index sections contiguous, minimize internal
storage used for constructing graph.

After construction and during initialization (when indices are being
added via InsertGlobalIndices() etc.), the column- indices for each
row are held in a separate piece of allocated memory. This method
moves the column-indices for all rows into one large contiguous array
and eliminates internal storage that is not needed after graph
construction. Calling this method can have a significant impact on
memory costs and machine performance.

If this object was constructed in View mode then this method can't
make non-contiguous indices contiguous and will return a warning code
of 1 if the viewed data isn't already contiguous. Integer error code,
set to 0 if successful.

Filled()==true.

If CV=View when the graph was constructed, then this method will be
effective  if the indices of the graph were already contiguous. In
this case, the indices are left untouched and internal storage for the
graph is minimized.

StorageOptimized()==true, if successful ";

/*  Extraction methods  */

%feature("docstring")  Epetra_CrsGraph::ExtractGlobalRowCopy "int
Epetra_CrsGraph::ExtractGlobalRowCopy(int GlobalRow, int LenOfIndices,
int &NumIndices, int *Indices) const

Extract a list of elements in a specified global row of the graph. Put
into storage allocated by calling routine.

Parameters:
-----------

Row:  - (In) Global row number to get indices.

LenOfIndices:  - (In) Length of Indices array.

NumIndices:  - (Out) Number of Indices.

Indices:  - (Out) Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_CrsGraph::ExtractGlobalRowCopy "int
Epetra_CrsGraph::ExtractGlobalRowCopy(long long GlobalRow, int
LenOfIndices, int &NumIndices, long long *Indices) const ";

%feature("docstring")  Epetra_CrsGraph::ExtractMyRowCopy "int
Epetra_CrsGraph::ExtractMyRowCopy(int LocalRow, int LenOfIndices, int
&NumIndices, int *Indices) const

Extract a list of elements in a specified local row of the graph. Put
into storage allocated by calling routine.

Parameters:
-----------

Row:  - (In) Local row number to get indices.

LenOfIndices:  - (In) Length of Indices array.

NumIndices:  - (Out) Number of Indices.

Indices:  - (Out) Local column indices corresponding to values.

Integer error code, set to 0 if successful.

IndicesAreLocal()==true ";

%feature("docstring")  Epetra_CrsGraph::ExtractGlobalRowView "int
Epetra_CrsGraph::ExtractGlobalRowView(int GlobalRow, int &NumIndices,
int *&Indices) const

Get a view of the elements in a specified global row of the graph.

This function requires that the graph not be completed (
FillComplete() was not called).

Parameters:
-----------

Row:  - (In) Global row number to get indices.

NumIndices:  - (Out) Number of Indices.

Indices:  - (Out) Global column indices corresponding to values.

Integer error code, set to 0 if successful. Returns -1 if invalid row.
Returns -2 if graph is completed.

IndicesAreLocal()==false ";

%feature("docstring")  Epetra_CrsGraph::ExtractGlobalRowView "int
Epetra_CrsGraph::ExtractGlobalRowView(long long GlobalRow, int
&NumIndices, long long *&Indices) const ";

%feature("docstring")  Epetra_CrsGraph::ExtractMyRowView "int
Epetra_CrsGraph::ExtractMyRowView(int LocalRow, int &NumIndices, int
*&Indices) const

Get a view of the elements in a specified local row of the graph.

This function requires that the graph be completed FillComplete() was
called).

Parameters:
-----------

Row:  - (In) Local row number to get indices.

NumIndices:  - (Out) Number of Indices.

Indices:  - (Out) Column indices corresponding to values.

Integer error code, set to 0 if successful. Returns -1 if invalid row.
Returns -2 if graph is not completed.

IndicesAreLocal()==true ";

/*  Graph Properties Query Methods  */

%feature("docstring")  Epetra_CrsGraph::Filled "bool
Epetra_CrsGraph::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  Epetra_CrsGraph::StorageOptimized "bool
Epetra_CrsGraph::StorageOptimized() const

If OptimizeStorage() has been called, this query returns true,
otherwise it returns false. ";

%feature("docstring")  Epetra_CrsGraph::IndicesAreGlobal "bool
Epetra_CrsGraph::IndicesAreGlobal() const

If column indices are in global range, this query returns true,
otherwise it returns false. ";

%feature("docstring")  Epetra_CrsGraph::IndicesAreLocal "bool
Epetra_CrsGraph::IndicesAreLocal() const

If column indices are in local range, this query returns true,
otherwise it returns false. ";

%feature("docstring")  Epetra_CrsGraph::LowerTriangular "bool
Epetra_CrsGraph::LowerTriangular() const

If graph is lower triangular in local index space, this query returns
true, otherwise it returns false.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::UpperTriangular "bool
Epetra_CrsGraph::UpperTriangular() const

If graph is upper triangular in local index space, this query returns
true, otherwise it returns false.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NoDiagonal "bool
Epetra_CrsGraph::NoDiagonal() const

If graph has no diagonal entries in global index space, this query
returns true, otherwise it returns false.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::MyGlobalRow "bool
Epetra_CrsGraph::MyGlobalRow(int GID) const

Returns true of GID is owned by the calling processor, otherwise it
returns false. ";

%feature("docstring")  Epetra_CrsGraph::MyGlobalRow "bool
Epetra_CrsGraph::MyGlobalRow(long long GID) const ";

%feature("docstring")  Epetra_CrsGraph::HaveColMap "bool
Epetra_CrsGraph::HaveColMap() const

Returns true if we have a well-defined ColMap, and returns false
otherwise.

We have a well-defined ColMap if a) a ColMap was passed in at
construction, or b) the MakeColMap function has been called. (Calling
either of the FillComplete functions will result in MakeColMap being
called.) ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_CrsGraph::NumMyRows "int
Epetra_CrsGraph::NumMyRows() const

Returns the number of matrix rows on this processor. ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalRows "int
Epetra_CrsGraph::NumGlobalRows() const

Returns the number of matrix rows in global matrix. ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalRows64 "long long
Epetra_CrsGraph::NumGlobalRows64() const ";

%feature("docstring")  Epetra_CrsGraph::NumMyCols "int
Epetra_CrsGraph::NumMyCols() const

Returns the number of entries in the set of column-indices that appear
on this processor.

The set of column-indices that appear on this processor is the union
of column-indices that appear in all local rows. The size of this set
isn't available until FillComplete() has been called.  Filled()==true
";

%feature("docstring")  Epetra_CrsGraph::NumGlobalCols "int
Epetra_CrsGraph::NumGlobalCols() const

Returns the number of matrix columns in global matrix.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalCols64 "long long
Epetra_CrsGraph::NumGlobalCols64() const ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalNonzeros "int
Epetra_CrsGraph::NumGlobalNonzeros() const

Returns the number of indices in the global graph.

Note that if the graph's maps are defined such that some nonzeros
appear on more than one processor, then those nonzeros will be counted
more than once. If the user wishes to assemble a graph from
overlapping data, they can use Epetra_FECrsGraph.  Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalNonzeros64 "long
long Epetra_CrsGraph::NumGlobalNonzeros64() const ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalDiagonals "int
Epetra_CrsGraph::NumGlobalDiagonals() const

Returns the number of diagonal entries in the global graph, based on
global row/column index comparisons.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalDiagonals64 "long
long Epetra_CrsGraph::NumGlobalDiagonals64() const ";

%feature("docstring")  Epetra_CrsGraph::NumMyDiagonals "int
Epetra_CrsGraph::NumMyDiagonals() const

Returns the number of diagonal entries in the local graph, based on
global row/column index comparisons.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumMyBlockRows "int
Epetra_CrsGraph::NumMyBlockRows() const

Returns the number of block matrix rows on this processor. ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalBlockRows "int
Epetra_CrsGraph::NumGlobalBlockRows() const

Returns the number of Block matrix rows in global matrix. ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalBlockRows64 "long
long Epetra_CrsGraph::NumGlobalBlockRows64() const ";

%feature("docstring")  Epetra_CrsGraph::NumMyBlockCols "int
Epetra_CrsGraph::NumMyBlockCols() const

Returns the number of Block matrix columns on this processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalBlockCols "int
Epetra_CrsGraph::NumGlobalBlockCols() const

Returns the number of Block matrix columns in global matrix.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalBlockCols64 "long
long Epetra_CrsGraph::NumGlobalBlockCols64() const ";

%feature("docstring")  Epetra_CrsGraph::NumMyBlockDiagonals "int
Epetra_CrsGraph::NumMyBlockDiagonals() const

Returns the number of Block diagonal entries in the local graph, based
on global row/column index comparisons.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalBlockDiagonals "int
Epetra_CrsGraph::NumGlobalBlockDiagonals() const

Returns the number of Block diagonal entries in the global graph,
based on global row/column index comparisons.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalBlockDiagonals64 "long long Epetra_CrsGraph::NumGlobalBlockDiagonals64() const ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalEntries "int
Epetra_CrsGraph::NumGlobalEntries() const

Returns the number of entries in the global graph.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalEntries64 "long long
Epetra_CrsGraph::NumGlobalEntries64() const ";

%feature("docstring")  Epetra_CrsGraph::NumMyEntries "int
Epetra_CrsGraph::NumMyEntries() const

Returns the number of entries on this processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::MaxRowDim "int
Epetra_CrsGraph::MaxRowDim() const

Returns the max row dimension of block entries on the processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::GlobalMaxRowDim "int
Epetra_CrsGraph::GlobalMaxRowDim() const

Returns the max row dimension of block entries across all processors.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::MaxColDim "int
Epetra_CrsGraph::MaxColDim() const

Returns the max column dimension of block entries on the processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::GlobalMaxColDim "int
Epetra_CrsGraph::GlobalMaxColDim() const

Returns the max column dimension of block entries across all
processors.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumMyNonzeros "int
Epetra_CrsGraph::NumMyNonzeros() const

Returns the number of indices in the local graph.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumGlobalIndices "int
Epetra_CrsGraph::NumGlobalIndices(long long Row) const

Returns the current number of nonzero entries in specified global row
on this processor. ";

%feature("docstring")  Epetra_CrsGraph::NumAllocatedGlobalIndices "int Epetra_CrsGraph::NumAllocatedGlobalIndices(long long Row) const

Returns the allocated number of nonzero entries in specified global
row on this processor. ";

%feature("docstring")  Epetra_CrsGraph::MaxNumIndices "int
Epetra_CrsGraph::MaxNumIndices() const

Returns the maximum number of nonzero entries across all rows on this
processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::GlobalMaxNumIndices "int
Epetra_CrsGraph::GlobalMaxNumIndices() const

Returns the maximun number of nonzero entries across all rows across
all processors.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::MaxNumNonzeros "int
Epetra_CrsGraph::MaxNumNonzeros() const

Returns the maximum number of nonzero points across all rows on this
processor.

For each entry in the graph, let i = the GRID of the entry and j = the
CGID of the entry. Then the entry size is the product of the rowmap
elementsize of i and the colmap elementsize of i. Let ki = sum of all
entry sizes for the entries in the ith row. For example, if the ith
block row had 5 block entries and the element size of each entry was
4-by-4, ki would be 80. Then this function returns the max over all ki
for all row on this processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::GlobalMaxNumNonzeros "int
Epetra_CrsGraph::GlobalMaxNumNonzeros() const

Returns the maximun number of nonzero points across all rows across
all processors.

This function returns the max over all processor of MaxNumNonzeros().

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::NumMyIndices "int
Epetra_CrsGraph::NumMyIndices(int Row) const

Returns the current number of nonzero entries in specified local row
on this processor. ";

%feature("docstring")  Epetra_CrsGraph::NumAllocatedMyIndices "int
Epetra_CrsGraph::NumAllocatedMyIndices(int Row) const

Returns the allocated number of nonzero entries in specified local row
on this processor. ";

%feature("docstring")  Epetra_CrsGraph::IndexBase "int
Epetra_CrsGraph::IndexBase() const

Returns the index base for row and column indices for this graph. ";

%feature("docstring")  Epetra_CrsGraph::RowMap "const
Epetra_BlockMap& Epetra_CrsGraph::RowMap() const

Returns the RowMap associated with this graph. ";

%feature("docstring")  Epetra_CrsGraph::ReplaceRowMap "int
Epetra_CrsGraph::ReplaceRowMap(const Epetra_BlockMap &newmap)

Replaces the current RowMap with the user-specified map object, but
only if currentmap->PointSameAs(newmap) is true. This is a collective
function. Returns 0 if map is replaced, -1 if not.

RowMap().PointSameAs(newmap)==true ";

%feature("docstring")  Epetra_CrsGraph::ReplaceColMap "int
Epetra_CrsGraph::ReplaceColMap(const Epetra_BlockMap &newmap)

Replaces the current ColMap with the user-specified map object, but
only if no entries have been inserted into the graph yet (both
IndicesAreLocal() and IndicesAreGlobal() are false) or
currentmap->PointSameAs(newmap) is true. This is a collective
function. Returns 0 if map is replaced, -1 if not.

( IndicesAreLocal()==false && IndicesAreGlobal()==false) ||
ColMap().PointSameAs(newmap)==true ";

%feature("docstring")  Epetra_CrsGraph::ColMap "const
Epetra_BlockMap& Epetra_CrsGraph::ColMap() const

Returns the Column Map associated with this graph.

HaveColMap()==true ";

%feature("docstring")  Epetra_CrsGraph::DomainMap "const
Epetra_BlockMap& Epetra_CrsGraph::DomainMap() const

Returns the DomainMap associated with this graph.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::RangeMap "const
Epetra_BlockMap& Epetra_CrsGraph::RangeMap() const

Returns the RangeMap associated with this graph.

Filled()==true ";

%feature("docstring")  Epetra_CrsGraph::Importer "const
Epetra_Import* Epetra_CrsGraph::Importer() const

Returns the Importer associated with this graph. ";

%feature("docstring")  Epetra_CrsGraph::Exporter "const
Epetra_Export* Epetra_CrsGraph::Exporter() const

Returns the Exporter associated with this graph. ";

%feature("docstring")  Epetra_CrsGraph::Comm "const Epetra_Comm&
Epetra_CrsGraph::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
graph. ";

/*  Local/Global ID methods  */

%feature("docstring")  Epetra_CrsGraph::LRID "int
Epetra_CrsGraph::LRID(int GRID_in) const

Returns the local row index for given global row index, returns -1 if
no local row for this global row. ";

%feature("docstring")  Epetra_CrsGraph::LRID "int
Epetra_CrsGraph::LRID(long long GRID_in) const ";

%feature("docstring")  Epetra_CrsGraph::GRID "int
Epetra_CrsGraph::GRID(int LRID_in) const

Returns the global row index for give local row index, returns
IndexBase-1 if we don't have this local row. ";

%feature("docstring")  Epetra_CrsGraph::GRID64 "long long
Epetra_CrsGraph::GRID64(int LRID_in) const ";

%feature("docstring")  Epetra_CrsGraph::LCID "int
Epetra_CrsGraph::LCID(int GCID_in) const

Returns the local column index for given global column index, returns
-1 if no local column for this global column.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsGraph::LCID "int
Epetra_CrsGraph::LCID(long long GCID_in) const ";

%feature("docstring")  Epetra_CrsGraph::GCID "int
Epetra_CrsGraph::GCID(int LCID_in) const

Returns the global column index for give local column index, returns
IndexBase-1 if we don't have this local column.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsGraph::GCID64 "long long
Epetra_CrsGraph::GCID64(int LCID_in) const ";

%feature("docstring")  Epetra_CrsGraph::MyGRID "bool
Epetra_CrsGraph::MyGRID(int GRID_in) const

Returns true if the GRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_CrsGraph::MyGRID "bool
Epetra_CrsGraph::MyGRID(long long GRID_in) const ";

%feature("docstring")  Epetra_CrsGraph::MyLRID "bool
Epetra_CrsGraph::MyLRID(int LRID_in) const

Returns true if the LRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_CrsGraph::MyGCID "bool
Epetra_CrsGraph::MyGCID(int GCID_in) const

Returns true if the GCID passed in belongs to the calling processor in
this map, otherwise returns false.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsGraph::MyGCID "bool
Epetra_CrsGraph::MyGCID(long long GCID_in) const ";

%feature("docstring")  Epetra_CrsGraph::MyLCID "bool
Epetra_CrsGraph::MyLCID(int LCID_in) const

Returns true if the LRID passed in belongs to the calling processor in
this map, otherwise returns false.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

/*  Inlined Operator Methods  */

/*  I/O Methods  */

%feature("docstring")  Epetra_CrsGraph::Print "void
Epetra_CrsGraph::Print(ostream &os) const

Print method. ";

%feature("docstring")  Epetra_CrsGraph::PrintGraphData "void
Epetra_CrsGraph::PrintGraphData(ostream &os) const ";

%feature("docstring")  Epetra_CrsGraph::PrintGraphData "void
Epetra_CrsGraph::PrintGraphData(ostream &os, int level) const ";

/*  Deprecated methods:  These methods still work, but will be removed
in a future version  */

%feature("docstring")  Epetra_CrsGraph::ImportMap "const
Epetra_BlockMap& Epetra_CrsGraph::ImportMap() const

Use ColMap() instead. ";

%feature("docstring")  Epetra_CrsGraph::TransformToLocal "int
Epetra_CrsGraph::TransformToLocal()

Use FillComplete() instead. ";

%feature("docstring")  Epetra_CrsGraph::TransformToLocal "int
Epetra_CrsGraph::TransformToLocal(const Epetra_BlockMap *DomainMap,
const Epetra_BlockMap *RangeMap)

Use FillComplete(const Epetra_BlockMap& DomainMap, const
Epetra_BlockMap& RangeMap) instead. ";

/*  Expert Users and Developers Only  */

%feature("docstring")  Epetra_CrsGraph::ReferenceCount "int
Epetra_CrsGraph::ReferenceCount() const

Returns the reference count of CrsGraphData.

(Intended for testing purposes.) ";

%feature("docstring")  Epetra_CrsGraph::DataPtr "const
Epetra_CrsGraphData* Epetra_CrsGraph::DataPtr() const

Returns a pointer to the CrsGraphData instance this CrsGraph uses.

(Intended for developer use only for testing purposes.) ";

%feature("docstring")
Epetra_CrsGraph::SortGhostsAssociatedWithEachProcessor "void
Epetra_CrsGraph::SortGhostsAssociatedWithEachProcessor(bool Flag)

Forces FillComplete() to locally order ghostnodes associated with each
remote processor in ascending order.

To be compliant with AztecOO, FillComplete() already locally orders
ghostnodes such that information received from processor k has a lower
local numbering than information received from processor j if k is
less than j. SortGhostsAssociatedWithEachProcessor(True) further
forces FillComplete() to locally number all ghostnodes received from
processor k in ascending order. That is, the local numbering of b is
less than c if the global numbering of b is less than c and if both b
and c are owned by the same processor. This is done to be compliant
with some limited block features within ML. In particular, some ML
features require that a block structure of the matrix be maintained
even within the ghost variables. ";


// File: classEpetra__CrsGraphData.xml
%feature("docstring") Epetra_CrsGraphData "

Epetra_CrsGraphData: The Epetra CrsGraph Data Class.

The Epetra_CrsGraphData class is an implementation detail of
Epetra_CrsGraph. It is reference-counted, and can be shared by
multiple Epetra_CrsGraph instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_CrsGraphData.h ";

/*  Constructor/Destructor Methods  */

/*  Helper methods called in CrsGraph. Mainly memory allocations and
deallocations.  */


// File: classEpetra__CrsMatrix.xml
%feature("docstring") Epetra_CrsMatrix "

Epetra_CrsMatrix: A class for constructing and using real-valued
double-precision sparse compressed row matrices.

The Epetra_CrsMatrix class is a sparse compressed row matrix object.
This matrix can be used in a parallel setting, with data distribution
described by Epetra_Map attributes. The structure or graph of the
matrix is defined by an Epetra_CrsGraph attribute.

In addition to coefficient access, the primary operations provided by
Epetra_CrsMatrix are matrix times vector and matrix times multi-vector
multiplication.

Epetra_CrsMatrix matrices can be square or rectangular.

Creating and filling Epetra_CrsMatrix objects

Constructing Epetra_CrsMatrix objects is a multi-step process. The
basic steps are as follows: Create Epetra_CrsMatrix instance,
including storage, via one of the constructors: Constructor that
accepts one Epetra_Map object, a row-map defining the distribution of
matrix rows.

Constructor that accepts two Epetra_Map objects. (The second map is a
column-map, and describes the set of column-indices that appear in
each processor's portion of the matrix. Generally these are
overlapping sets -- column-indices may appear on more than one
processor.)

Constructor that accepts an Epetra_CrsGraph object, defining the non-
zero structure of the matrix.  Note that the constructors which accept
Epetra_Map arguments also accept an argument that gives an estimate of
the number of nonzeros per row. This allows storage to be pre-
allocated and can improve the performance of the data input methods.
The estimate need not be accurate, as additional storage is allocated
automatically when needed. However, a more accurate estimate helps
performance by reducing the amount of extra memory allocation.

Enter values via one or more Insert/Replace/SumInto functions.

Complete construction by calling FillComplete.

Note that, even after a matrix is constructed (FillComplete has been
called), it is possible to update existing matrix entries. It is not
possible to create new entries.

Epetra_Map attributes

Epetra_CrsMatrix objects have four Epetra_Map attributes, which are
held by the Epetra_CrsGraph attribute.

The Epetra_Map attributes can be obtained via these accessor methods:
RowMap() Describes the numbering and distribution of the rows of the
matrix. The row-map exists and is valid for the entire life of the
matrix. The set of matrix rows is defined by the row-map and may not
be changed. Rows may not be inserted or deleted by the user. The only
change that may be made is that the user can replace the row-map with
a compatible row-map (which is the same except for re-numbering) by
calling the ReplaceRowMap() method.

ColMap() Describes the set of column-indices that appear in the rows
in each processor's portion of the matrix. Unless provided by the user
at construction time, a valid column-map doesn't exist until
FillComplete() is called.

RangeMap() Describes the range of the matrix operator. e.g., for a
matrix-vector product operation, the result vector's map must be
compatible with the range-map of this matrix. The range-map is usually
the same as the row-map. The range-map is set equal to the row-map at
matrix creation time, but may be specified by the user when
FillComplete() is called.

DomainMap() Describes the domain of the matrix operator. The domain-
map can be specified by the user when FillComplete() is called. Until
then, it is set equal to the row-map.

It is important to note that while the row-map and the range-map are
often the same, the column-map and the domain-map are almost never the
same. The set of entries in a distributed column-map almost always
form overlapping sets, with entries being associated with more than
one processor. A domain-map, on the other hand, must be a 1-to-1 map,
with entries being associated with only a single processor.

Local versus Global Indices

Epetra_CrsMatrix has query functions IndicesAreLocal() and
IndicesAreGlobal(), which are used to determine whether the underlying
Epetra_CrsGraph attribute's column-indices have been transformed into
a local index space or not. (This transformation occurs when the
method Epetra_CrsGraph::FillComplete() is called, which happens when
the method Epetra_CrsMatrix::FillComplete() is called.) The state of
the indices in the graph determines the behavior of many
Epetra_CrsMatrix methods. If an Epetra_CrsMatrix instance is
constructed using one of the constructors that does not accept a pre-
existing Epetra_CrsGraph object, then an Epetra_CrsGraph attribute is
created internally and its indices remain untransformed (
IndicesAreGlobal()==true) until Epetra_CrsMatrix::FillComplete() is
called. The query function Epetra_CrsMatrix::Filled() returns true if
Epetra_CrsMatrix::FillComplete() has been called.

Note the following method characteristics:

InsertGlobalValues() may only be used to insert new nonzeros in the
matrix if indices are global.

SumIntoGlobalValues() may be used regardless of whether indices are
global or local, but can only be used to update matrix locations that
already exist; it can never be used to establish new nonzero
locations.

ReplaceGlobalValues() may also be used only to update matrix locations
that already exist, and works regardless of whether indices are local
or global.

SumIntoMyValues() and ReplaceMyValues() may only be used if indices
are local.

Multiply() may only be used after FillComplete() has been called.

Most methods have preconditions documented, check documentation for
specific methods not mentioned here.

Counting Floating Point Operations

Each Epetra_CrsMatrix object keeps track of the number of serial
floating point operations performed using the specified object as the
this argument to the function. The Flops() function returns this
number as a double precision number. Using this information, in
conjunction with the Epetra_Time class, one can get accurate parallel
performance numbers. The ResetFlops() function resets the floating
point counter.

WARNING:  A Epetra_Map is required for the Epetra_CrsMatrix
constructor.

C++ includes: Epetra_CrsMatrix.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_CrsMatrix::Epetra_CrsMatrix "Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, const int *NumEntriesPerRow, bool
StaticProfile=false)

Epetra_CrsMatrix constructor with variable number of indices per row.

Creates a Epetra_CrsMatrix object and allocates storage.

Parameters:
-----------

CV:  - (In) An Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_Map defining the numbering and distribution
of matrix rows.

NumEntriesPerRow:  - (In) An integer array of length NumRows such that
NumEntriesPerRow[i] indicates the (approximate if StaticProfile=false)
number of entries in the ith row.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsMatrix::Epetra_CrsMatrix "Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, int NumEntriesPerRow, bool StaticProfile=false)

Epetra_CrsMatrix constructor with fixed number of indices per row.

Creates a Epetra_CrsMatrix object and allocates storage.

Parameters:
-----------

CV:  - (In) An Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_Map defining the numbering and distribution
of matrix rows.

NumEntriesPerRow:  - (In) An integer that indicates the (approximate)
number of entries in the each row. Note that it is possible to use 0
for this value and let fill occur during the insertion phase.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsMatrix::Epetra_CrsMatrix "Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, const Epetra_Map &ColMap, const int
*NumEntriesPerRow, bool StaticProfile=false)

Epetra_CrsMatrix constructor with variable number of indices per row.

Creates a Epetra_CrsMatrix object and allocates storage.

Parameters:
-----------

CV:  - (In) An Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_Map defining the numbering and distribution
of matrix rows.

ColMap:  - (In) An Epetra_Map defining the set of column-indices that
appear in each processor's locally owned matrix rows.

NumEntriesPerRow:  - (In) An integer array of length NumRows such that
NumEntriesPerRow[i] indicates the (approximate if StaticProfile=false)
number of entries in the ith row.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsMatrix::Epetra_CrsMatrix "Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, const Epetra_Map &ColMap, int NumEntriesPerRow,
bool StaticProfile=false)

Epetra_CrsMatrix constuctor with fixed number of indices per row.

Creates a Epetra_CrsMatrix object and allocates storage.

Parameters:
-----------

CV:  - (In) An Epetra_DataAccess enumerated type set to Copy or View.

RowMap:  - (In) An Epetra_Map defining the numbering and distribution
of matrix rows.

ColMap:  - (In) An Epetra_Map defining the set of column-indices that
appear in each processor's locally owned matrix rows.

NumEntriesPerRow:  - (In) An integer that indicates the (approximate
if StaticProfile=false) number of entries in the each row. Note that
it is possible to use 0 for this value and let fill occur during the
insertion phase.

StaticProfile:  - (In) Optional argument that indicates whether or not
NumIndicesPerRow should be interpreted as an exact count of nonzeros,
or should be used as an approximation. By default this value is false,
allowing the profile to be determined dynamically. If the user sets it
to true, then the memory allocation for the Epetra_CrsGraph object
will be done in one large block, saving on memory fragmentation and
generally improving the performance of matrix multiplication and solve
kernels. ";

%feature("docstring")  Epetra_CrsMatrix::Epetra_CrsMatrix "Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const
Epetra_CrsGraph &Graph)

Construct a matrix using an existing Epetra_CrsGraph object.

Allows the nonzero structure from another matrix, or a structure that
was constructed independently, to be used for this matrix.

Parameters:
-----------

CV:  - (In) An Epetra_DataAccess enumerated type set to Copy or View.

Graph:  - (In) A Epetra_CrsGraph object, constructed directly or
extracted from another Epetra matrix object. ";

%feature("docstring")  Epetra_CrsMatrix::Epetra_CrsMatrix "Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix &Matrix)

Copy constructor. ";

%feature("docstring")  Epetra_CrsMatrix::~Epetra_CrsMatrix "Epetra_CrsMatrix::~Epetra_CrsMatrix()

Epetra_CrsMatrix Destructor. ";

/*  Insertion/Replace/SumInto methods  */

%feature("docstring")  Epetra_CrsMatrix::PutScalar "int
Epetra_CrsMatrix::PutScalar(double ScalarConstant)

Initialize all values in the matrix with constant value.

Parameters:
-----------

ScalarConstant:  - (In) Value to use.

Integer error code, set to 0 if successful.

None.

All values in this set to ScalarConstant. ";

%feature("docstring")  Epetra_CrsMatrix::Scale "int
Epetra_CrsMatrix::Scale(double ScalarConstant)

Multiply all values in the matrix by a constant value (in place: A <-
ScalarConstant * A).

Parameters:
-----------

ScalarConstant:  - (In) Value to use.

Integer error code, set to 0 if successful.

None.

All values of this have been multiplied by ScalarConstant. ";

%feature("docstring")  Epetra_CrsMatrix::InsertGlobalValues "int
Epetra_CrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
const double *Values, const int *Indices)

Insert a list of elements in a given global row of the matrix.

This method is used to construct a matrix for the first time. It
cannot be used if the matrix structure has already been fixed (via a
call to FillComplete()). If multiple values are inserted for the same
matrix entry, the values are initially stored separately, so memory
use will grow as a result. However, when FillComplete is called the
values will be summed together and the additional memory will be
released.

For example, if the values 2.0, 3.0 and 4.0 are all inserted in Row 1,
Column 2, extra storage is used to store each of the three values
separately. In this way, the insert process does not require any
searching and can be faster. However, when FillComplete() is called,
the values will be summed together to equal 9.0 and only a single
entry will remain in the matrix for Row 1, Column 2.

Parameters:
-----------

GlobalRow:  - (In) Row number (in global coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Global column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, a positive warning code will be
returned.

WARNING:  This method may not be called once FillComplete() has been
called.

IndicesAreLocal()==false && IndicesAreContiguous()==false ";

%feature("docstring")  Epetra_CrsMatrix::InsertGlobalValues "int
Epetra_CrsMatrix::InsertGlobalValues(long long GlobalRow, int
NumEntries, const double *Values, const long long *Indices) ";

%feature("docstring")  Epetra_CrsMatrix::InsertGlobalValues "int
Epetra_CrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
double *Values, int *Indices)

Insert a list of elements in a given global row of the matrix.

This method is used to construct a matrix for the first time. It
cannot be used if the matrix structure has already been fixed (via a
call to FillComplete()). If multiple values are inserted for the same
matrix entry, the values are initially stored separately, so memory
use will grow as a result. However, when FillComplete is called the
values will be summed together and the additional memory will be
released.

For example, if the values 2.0, 3.0 and 4.0 are all inserted in Row 1,
Column 2, extra storage is used to store each of the three values
separately. In this way, the insert process does not require any
searching and can be faster. However, when FillComplete() is called,
the values will be summed together to equal 9.0 and only a single
entry will remain in the matrix for Row 1, Column 2.

Parameters:
-----------

GlobalRow:  - (In) Row number (in global coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Global column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, a positive warning code will be
returned.

WARNING:  This method may not be called once FillComplete() has been
called.

IndicesAreLocal()==false && IndicesAreContiguous()==false ";

%feature("docstring")  Epetra_CrsMatrix::InsertGlobalValues "int
Epetra_CrsMatrix::InsertGlobalValues(long long GlobalRow, int
NumEntries, double *Values, long long *Indices) ";

%feature("docstring")  Epetra_CrsMatrix::ReplaceGlobalValues "int
Epetra_CrsMatrix::ReplaceGlobalValues(int GlobalRow, int NumEntries,
const double *Values, const int *Indices)

Replace specified existing values with this list of entries for a
given global row of the matrix.

Parameters:
-----------

GlobalRow:  - (In) Row number (in global coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Global column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if a value is
not already present for the specified location in the matrix, the
input value will be ignored and a positive warning code will be
returned.

IndicesAreLocal()==false && IndicesAreContiguous()==false ";

%feature("docstring")  Epetra_CrsMatrix::ReplaceGlobalValues "int
Epetra_CrsMatrix::ReplaceGlobalValues(long long GlobalRow, int
NumEntries, const double *Values, const long long *Indices) ";

%feature("docstring")  Epetra_CrsMatrix::SumIntoGlobalValues "int
Epetra_CrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries,
const double *Values, const int *Indices)

Add this list of entries to existing values for a given global row of
the matrix.

Parameters:
-----------

GlobalRow:  - (In) Row number (in global coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Global column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if a value is
not already present for the specified location in the matrix, the
input value will be ignored and a positive warning code will be
returned.

IndicesAreLocal()==false && IndicesAreContiguous()==false ";

%feature("docstring")  Epetra_CrsMatrix::SumIntoGlobalValues "int
Epetra_CrsMatrix::SumIntoGlobalValues(long long GlobalRow, int
NumEntries, const double *Values, const long long *Indices) ";

%feature("docstring")  Epetra_CrsMatrix::InsertMyValues "int
Epetra_CrsMatrix::InsertMyValues(int MyRow, int NumEntries, const
double *Values, const int *Indices)

Insert a list of elements in a given local row of the matrix.

Parameters:
-----------

MyRow:  - (In) Row number (in local coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Local column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, a positive warning code will be
returned.

IndicesAreGlobal()==false && ( IndicesAreContiguous()==false ||
CV_==View)

The given local row of the matrix has been updated as described above.
";

%feature("docstring")  Epetra_CrsMatrix::InsertMyValues "int
Epetra_CrsMatrix::InsertMyValues(int MyRow, int NumEntries, double
*Values, int *Indices)

Insert a list of elements in a given local row of the matrix.

Parameters:
-----------

MyRow:  - (In) Row number (in local coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Local column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, a positive warning code will be
returned.

IndicesAreGlobal()==false && ( IndicesAreContiguous()==false ||
CV_==View)

The given local row of the matrix has been updated as described above.
";

%feature("docstring")  Epetra_CrsMatrix::ReplaceMyValues "int
Epetra_CrsMatrix::ReplaceMyValues(int MyRow, int NumEntries, const
double *Values, const int *Indices)

Replace current values with this list of entries for a given local row
of the matrix.

Parameters:
-----------

MyRow:  - (In) Row number (in local coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Local column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if a value is
not already present for the specified location in the matrix, the
input value will be ignored and a positive warning code will be
returned.

IndicesAreLocal()==true

MyRow contains the given list of Values at the given Indices. ";

%feature("docstring")  Epetra_CrsMatrix::SumIntoMyValues "int
Epetra_CrsMatrix::SumIntoMyValues(int MyRow, int NumEntries, const
double *Values, const int *Indices)

Add this list of entries to existing values for a given local row of
the matrix.

Parameters:
-----------

MyRow:  - (In) Row number (in local coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Local column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, a positive warning code will be
returned.

IndicesAreLocal()==true

The given Values at the given Indices have been summed into the
entries of MyRow. ";

%feature("docstring")  Epetra_CrsMatrix::ReplaceDiagonalValues "int
Epetra_CrsMatrix::ReplaceDiagonalValues(const Epetra_Vector &Diagonal)

Replaces diagonal values of the matrix with those in the user-provided
vector.

This routine is meant to allow replacement of { existing} diagonal
values. If a diagonal value does not exist for a given row, the
corresponding value in the input Epetra_Vector will be ignored and the
return code will be set to 1.

The Epetra_Map associated with the input Epetra_Vector must be
compatible with the RowMap of the matrix.

Parameters:
-----------

Diagonal:  - (In) New values to be placed in the main diagonal.

Integer error code, set to 0 if successful, set to 1 on the calling
processor if one or more diagonal entries not present in matrix.

Filled()==true

Diagonal values have been replaced with the values of Diagonal. ";

/*  Transformation methods  */

%feature("docstring")  Epetra_CrsMatrix::FillComplete "int
Epetra_CrsMatrix::FillComplete(bool OptimizeDataStorage=true)

Signal that data entry is complete. Perform transformations to local
index space. ";

%feature("docstring")  Epetra_CrsMatrix::FillComplete "int
Epetra_CrsMatrix::FillComplete(const Epetra_Map &DomainMap, const
Epetra_Map &RangeMap, bool OptimizeDataStorage=true)

Signal that data entry is complete. Perform transformations to local
index space. ";

%feature("docstring")  Epetra_CrsMatrix::OptimizeStorage "int
Epetra_CrsMatrix::OptimizeStorage()

Make consecutive row index sections contiguous, minimize internal
storage used for constructing graph.

After construction and during initialization (when values are being
added), the matrix coefficients for each row are managed as separate
segments of memory. This method moves the values for all rows into one
large contiguous array and eliminates internal storage that is not
needed after matrix construction. Calling this method can have a
significant impact on memory costs and machine performance.

If this object was constructed in View mode then this method can't
make non-contiguous values contiguous and will return a warning code
of 1 if the viewed data isn't already contiguous.

A call to this method will also call the OptimizeStorage method for
the associated Epetra_CrsGraph object. If the storage for this graph
has already been optimized this additional call will have no effect.

Integer error code, set to 0 if successful.

Filled()==true.

If CV=View when the graph was constructed, then this method will be
effective  if the indices of the graph were already contiguous. In
this case, the indices are left untouched and internal storage for the
graph is minimized.

StorageOptimized()==true, if successful.

Graph(). StorageOptimized()==true, if successful. ";

%feature("docstring")  Epetra_CrsMatrix::MakeDataContiguous "int
Epetra_CrsMatrix::MakeDataContiguous()

Eliminates memory that is used for construction. Make consecutive row
index sections contiguous. ";

/*  Extraction methods  */

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowCopy "int
Epetra_CrsMatrix::ExtractGlobalRowCopy(int GlobalRow, int Length, int
&NumEntries, double *Values, int *Indices) const

Returns a copy of the specified global row in user-provided arrays.

Parameters:
-----------

GlobalRow:  - (In) Global row to extract.

ILength:  - (In) Length of Values and Indices.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Indices:  - (Out) Extracted global column indices for the
corresponding values.

Integer error code, set to 0 if successful, non-zero if global row is
not owned by calling process or if the number of entries in this row
exceed the Length parameter. ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowCopy "int
Epetra_CrsMatrix::ExtractGlobalRowCopy(long long GlobalRow, int
Length, int &NumEntries, double *Values, long long *Indices) const ";

%feature("docstring")  Epetra_CrsMatrix::ExtractMyRowCopy "int
Epetra_CrsMatrix::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

MyRow:  - (In) Local row to extract.

Length:  - (In) Length of Values and Indices.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Indices:  - (Out) Extracted local column indices for the corresponding
values.

Integer error code, set to 0 if successful.

IndicesAreLocal()==true ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowCopy "int
Epetra_CrsMatrix::ExtractGlobalRowCopy(int GlobalRow, int Length, int
&NumEntries, double *Values) const

Returns a copy of the specified global row values in user-provided
array.

Parameters:
-----------

GlobalRow:  - (In) Global row to extract.

Length:  - (In) Length of Values.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowCopy "int
Epetra_CrsMatrix::ExtractGlobalRowCopy(long long GlobalRow, int
Length, int &NumEntries, double *Values) const ";

%feature("docstring")  Epetra_CrsMatrix::ExtractMyRowCopy "int
Epetra_CrsMatrix::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values) const

Returns a copy of the specified local row values in user-provided
array.

Parameters:
-----------

MyRow:  - (In) Local row to extract.

Length:  - (In) Length of Values.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_CrsMatrix::ExtractDiagonalCopy "int
Epetra_CrsMatrix::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Diagonal:  - (Out) Extracted main diagonal.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowView "int
Epetra_CrsMatrix::ExtractGlobalRowView(int GlobalRow, int &NumEntries,
double *&Values, int *&Indices) const

Returns a view of the specified global row values via pointers to
internal data.

Parameters:
-----------

GlobalRow:  - (In) Global row to view.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Indices:  - (Out) Extracted global column indices for the
corresponding values.

Integer error code, set to 0 if successful. Returns -1 of row not on
this processor. Returns -2 if matrix is not in global form (if
FillComplete() has already been called).

IndicesAreGlobal()==true ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowView "int
Epetra_CrsMatrix::ExtractGlobalRowView(long long GlobalRow, int
&NumEntries, double *&Values, long long *&Indices) const ";

%feature("docstring")  Epetra_CrsMatrix::ExtractMyRowView "int
Epetra_CrsMatrix::ExtractMyRowView(int MyRow, int &NumEntries, double
*&Values, int *&Indices) const

Returns a view of the specified local row values via pointers to
internal data.

Parameters:
-----------

MyRow:  - (In) Local row to view.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Indices:  - (Out) Extracted local column indices for the corresponding
values.

Integer error code, set to 0 if successful. Returns -1 of row not on
this processor. Returns -2 if matrix is not in local form (if
FillComplete() has not been called).

IndicesAreLocal()==true ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowView "int
Epetra_CrsMatrix::ExtractGlobalRowView(int GlobalRow, int &NumEntries,
double *&Values) const

Returns a view of the specified global row values via pointers to
internal data.

Parameters:
-----------

GlobalRow:  - (In) Global row to extract.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_CrsMatrix::ExtractGlobalRowView "int
Epetra_CrsMatrix::ExtractGlobalRowView(long long GlobalRow, int
&NumEntries, double *&Values) const ";

%feature("docstring")  Epetra_CrsMatrix::ExtractMyRowView "int
Epetra_CrsMatrix::ExtractMyRowView(int MyRow, int &NumEntries, double
*&Values) const

Returns a view of the specified local row values via pointers to
internal data.

Parameters:
-----------

MyRow:  - (In) Local row to extract.

NumEntries:  - (Out) Number of nonzero entries extracted.

Values:  - (Out) Extracted values for this row.

Integer error code, set to 0 if successful. ";

/*  Computational methods  */

%feature("docstring")  Epetra_CrsMatrix::Multiply "int
Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y) const

Returns the result of a Epetra_CrsMatrix multiplied by a Epetra_Vector
x in y.

Parameters:
-----------

TransA:  - (In) If true, multiply by the transpose of matrix,
otherwise just use matrix.

x:  - (In) An Epetra_Vector to multiply by.

y:  - (Out) An Epetra_Vector containing result.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::Multiply1 "int
Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y) const ";

%feature("docstring")  Epetra_CrsMatrix::Multiply "int
Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_CrsMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

TransA:  - (In) If true, multiply by the transpose of matrix,
otherwise just use matrix.

X:  - (In) An Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  - (Out) An Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::Multiply1 "int
Epetra_CrsMatrix::Multiply1(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  Epetra_CrsMatrix::Solve "int
Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_Vector &x, Epetra_Vector &y) const

Returns the result of a local solve using the Epetra_CrsMatrix on a
Epetra_Vector x in y.

This method solves a triangular system of equations asynchronously on
each processor.

Parameters:
-----------

Upper:  - (In) If true, solve Uy = x, otherwise solve Ly = x.

Trans:  - (In) If true, solve transpose problem.

UnitDiagonal:  - (In) If true, assume diagonal is unit (whether it's
stored or not).

x:  - (In) An Epetra_Vector to solve for.

y:  - (Out) An Epetra_Vector containing result.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::Solve "int
Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of a local solve using the Epetra_CrsMatrix a
Epetra_MultiVector X in Y.

This method solves a triangular system of equations asynchronously on
each processor.

Parameters:
-----------

Upper:  - (In) If true, solve Uy = x, otherwise solve Ly = x.

Trans:  - (In) If true, solve transpose problem.

UnitDiagonal:  - (In) If true, assume diagonal is unit (whether it's
stored or not).

X:  - (In) An Epetra_MultiVector of dimension NumVectors to solve for.

Y:  - (Out) An Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::InvRowSums "int
Epetra_CrsMatrix::InvRowSums(Epetra_Vector &x) const

Computes the inverse of the sum of absolute values of the rows of the
Epetra_CrsMatrix, results returned in x.

The vector x will return such that x[i] will contain the inverse of
the sum of the absolute values of the entries in the ith row of the
this matrix. Using the resulting vector from this function as input to
LeftScale() will make the infinity norm of the resulting matrix
exactly 1. WARNING:  The NormInf() method will not properly calculate
the infinity norm for a matrix that has entries that are replicated on
multiple processors. In this case, if the rows are fully replicated,
NormInf() will return a value equal to the maximum number of
processors that any individual row of the matrix is replicated on.

Parameters:
-----------

x:  - (Out) An Epetra_Vector containing the inverse of the row sums of
the this matrix.

WARNING:  When rows are fully replicated on multiple processors, it is
assumed that the distribution of x is the same as the rows (
RowMap())of this. When multiple processors contain partial sums for
individual entries, the distribution of x is assumed to be the same as
the RangeMap() of this. When each row of this is uniquely owned, the
distribution of x can be that of the RowMap() or the RangeMap().

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::InvRowMaxs "int
Epetra_CrsMatrix::InvRowMaxs(Epetra_Vector &x) const

Computes the inverse of the max of absolute values of the rows of the
Epetra_CrsMatrix, results returned in x.

The vector x will return such that x[i] will contain the inverse of
max of the absolute values of the entries in the ith row of the this
matrix. WARNING:  This method will not work when multiple processors
contain partial sums for individual entries.

Parameters:
-----------

x:  - (Out) An Epetra_Vector containing the inverse of the row maxs of
the this matrix.

WARNING:  When rows are fully replicated on multiple processors, it is
assumed that the distribution of x is the same as the rows (
RowMap())of this. When each row of this is uniquely owned, the
distribution of x can be that of the RowMap() or the RangeMap().

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::LeftScale "int
Epetra_CrsMatrix::LeftScale(const Epetra_Vector &x)

Scales the Epetra_CrsMatrix on the left with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i
denotes the row number of A and j denotes the column number of A.

Parameters:
-----------

x:  - (In) An Epetra_Vector to scale with.

Integer error code, set to 0 if successful.

Filled()==true

The matrix will be scaled as described above. ";

%feature("docstring")  Epetra_CrsMatrix::InvColSums "int
Epetra_CrsMatrix::InvColSums(Epetra_Vector &x) const

Computes the inverse of the sum of absolute values of the columns of
the Epetra_CrsMatrix, results returned in x.

The vector x will return such that x[j] will contain the inverse of
the sum of the absolute values of the entries in the jth column of the
this matrix. Using the resulting vector from this function as input to
RightScale() will make the one norm of the resulting matrix exactly 1.
WARNING:  The NormOne() method will not properly calculate the one
norm for a matrix that has entries that are replicated on multiple
processors. In this case, if the columns are fully replicated,
NormOne() will return a value equal to the maximum number of
processors that any individual column of the matrix is repliated on.

Parameters:
-----------

x:  - (Out) An Epetra_Vector containing the column sums of the this
matrix.

WARNING:  When columns are fully replicated on multiple processors, it
is assumed that the distribution of x is the same as the columns (
ColMap()) of this. When multiple processors contain partial sums for
entries, the distribution of x is assumed to be the same as the
DomainMap() of this. When each column of this is uniquely owned, the
distribution of x can be that of the ColMap() or the DomainMap().

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::InvColMaxs "int
Epetra_CrsMatrix::InvColMaxs(Epetra_Vector &x) const

Computes the max of absolute values of the columns of the
Epetra_CrsMatrix, results returned in x.

The vector x will return such that x[j] will contain the inverse of
max of the absolute values of the entries in the jth row of the this
matrix. WARNING:  This method will not work when multiple processors
contain partial sums for individual entries.

Parameters:
-----------

x:  - (Out) An Epetra_Vector containing the column maxs of the this
matrix.

WARNING:  When columns are fully replicated on multiple processors, it
is assumed that the distribution of x is the same as the columns (
ColMap()) of this. When each column of this is uniquely owned, the
distribution of x can be that of the ColMap() or the DomainMap().

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::RightScale "int
Epetra_CrsMatrix::RightScale(const Epetra_Vector &x)

Scales the Epetra_CrsMatrix on the right with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i
denotes the global row number of A and j denotes the global column
number of A.

Parameters:
-----------

x:  - (In) The Epetra_Vector used for scaling this.

Integer error code, set to 0 if successful.

Filled()==true

The matrix will be scaled as described above. ";

/*  Matrix Properties Query Methods  */

%feature("docstring")  Epetra_CrsMatrix::Filled "bool
Epetra_CrsMatrix::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  Epetra_CrsMatrix::StorageOptimized "bool
Epetra_CrsMatrix::StorageOptimized() const

If OptimizeStorage() has been called, this query returns true,
otherwise it returns false. ";

%feature("docstring")  Epetra_CrsMatrix::IndicesAreGlobal "bool
Epetra_CrsMatrix::IndicesAreGlobal() const

If matrix indices has not been transformed to local, this query
returns true, otherwise it returns false. ";

%feature("docstring")  Epetra_CrsMatrix::IndicesAreLocal "bool
Epetra_CrsMatrix::IndicesAreLocal() const

If matrix indices has been transformed to local, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_CrsMatrix::IndicesAreContiguous "bool
Epetra_CrsMatrix::IndicesAreContiguous() const

If matrix indices are packed into single array (done in
OptimizeStorage()) return true, otherwise false. ";

%feature("docstring")  Epetra_CrsMatrix::LowerTriangular "bool
Epetra_CrsMatrix::LowerTriangular() const

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_CrsMatrix::UpperTriangular "bool
Epetra_CrsMatrix::UpperTriangular() const

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_CrsMatrix::NoDiagonal "bool
Epetra_CrsMatrix::NoDiagonal() const

If matrix has no diagonal entries in global index space, this query
returns true, otherwise it returns false. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_CrsMatrix::NormInf "double
Epetra_CrsMatrix::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Epetra_CrsMatrix::NormOne "double
Epetra_CrsMatrix::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")  Epetra_CrsMatrix::NormFrobenius "double
Epetra_CrsMatrix::NormFrobenius() const

Returns the frobenius norm of the global matrix. ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalNonzeros "int
Epetra_CrsMatrix::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalNonzeros64 "long
long Epetra_CrsMatrix::NumGlobalNonzeros64() const ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalRows "int
Epetra_CrsMatrix::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalRows64 "long long
Epetra_CrsMatrix::NumGlobalRows64() const ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalCols "int
Epetra_CrsMatrix::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalCols64 "long long
Epetra_CrsMatrix::NumGlobalCols64() const ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalDiagonals "int
Epetra_CrsMatrix::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalDiagonals64 "long
long Epetra_CrsMatrix::NumGlobalDiagonals64() const ";

%feature("docstring")  Epetra_CrsMatrix::NumMyNonzeros "int
Epetra_CrsMatrix::NumMyNonzeros() const

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Epetra_CrsMatrix::NumMyRows "int
Epetra_CrsMatrix::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Epetra_CrsMatrix::NumMyCols "int
Epetra_CrsMatrix::NumMyCols() const

Returns the number of entries in the set of column-indices that appear
on this processor.

The set of column-indices that appear on this processor is the union
of column-indices that appear in all local rows. The size of this set
isn't available until FillComplete() has been called.  Filled()==true
";

%feature("docstring")  Epetra_CrsMatrix::NumMyDiagonals "int
Epetra_CrsMatrix::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons.

Filled()==true ";

%feature("docstring")  Epetra_CrsMatrix::NumGlobalEntries "int
Epetra_CrsMatrix::NumGlobalEntries(long long Row) const

Returns the current number of nonzero entries in specified global row
on this processor. ";

%feature("docstring")  Epetra_CrsMatrix::NumAllocatedGlobalEntries "int Epetra_CrsMatrix::NumAllocatedGlobalEntries(int Row) const

Returns the allocated number of nonzero entries in specified global
row on this processor. ";

%feature("docstring")  Epetra_CrsMatrix::MaxNumEntries "int
Epetra_CrsMatrix::MaxNumEntries() const

Returns the maximum number of nonzero entries across all rows on this
processor.

Filled()==true ";

%feature("docstring")  Epetra_CrsMatrix::GlobalMaxNumEntries "int
Epetra_CrsMatrix::GlobalMaxNumEntries() const

Returns the maximum number of nonzero entries across all rows on all
processors.

Filled()==true ";

%feature("docstring")  Epetra_CrsMatrix::NumMyEntries "int
Epetra_CrsMatrix::NumMyEntries(int Row) const

Returns the current number of nonzero entries in specified local row
on this processor. ";

%feature("docstring")  Epetra_CrsMatrix::NumAllocatedMyEntries "int
Epetra_CrsMatrix::NumAllocatedMyEntries(int Row) const

Returns the allocated number of nonzero entries in specified local row
on this processor. ";

%feature("docstring")  Epetra_CrsMatrix::IndexBase "int
Epetra_CrsMatrix::IndexBase() const

Returns the index base for row and column indices for this graph. ";

%feature("docstring")  Epetra_CrsMatrix::StaticGraph "bool
Epetra_CrsMatrix::StaticGraph()

Returns true if the graph associated with this matrix was pre-
constructed and therefore not changeable. ";

%feature("docstring")  Epetra_CrsMatrix::Graph "const
Epetra_CrsGraph& Epetra_CrsMatrix::Graph() const

Returns a reference to the Epetra_CrsGraph object associated with this
matrix. ";

%feature("docstring")  Epetra_CrsMatrix::RowMap "const Epetra_Map&
Epetra_CrsMatrix::RowMap() const

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")  Epetra_CrsMatrix::ReplaceRowMap "int
Epetra_CrsMatrix::ReplaceRowMap(const Epetra_BlockMap &newmap)

Replaces the current RowMap with the user-specified map object.

Replaces the current RowMap with the user-specified map object, but
only if currentmap->PointSameAs(newmap) is true. This is a collective
function. Returns 0 if map is replaced, -1 if not.

RowMap().PointSameAs(newmap)==true ";

%feature("docstring")  Epetra_CrsMatrix::HaveColMap "bool
Epetra_CrsMatrix::HaveColMap() const

Returns true if we have a well-defined ColMap, and returns false
otherwise.

We have a well-defined ColMap if a) a ColMap was passed in at
construction, or b) the MakeColMap function has been called. (Calling
either of the FillComplete functions will result in MakeColMap being
called.) ";

%feature("docstring")  Epetra_CrsMatrix::ReplaceColMap "int
Epetra_CrsMatrix::ReplaceColMap(const Epetra_BlockMap &newmap)

Replaces the current ColMap with the user-specified map object.

Replaces the current ColMap with the user-specified map object, but
only if no entries have been inserted into the matrix (both
IndicesAreLocal() and IndicesAreGlobal() are false) or
currentmap->PointSameAs(newmap) is true. This is a collective
function. Returns 0 if map is replaced, -1 if not.

( IndicesAreLocal()==false && IndicesAreGlobal()==false) ||
ColMap().PointSameAs(newmap)==true ";

%feature("docstring")  Epetra_CrsMatrix::ColMap "const Epetra_Map&
Epetra_CrsMatrix::ColMap() const

Returns the Epetra_Map object that describes the set of column-indices
that appear in each processor's locally owned matrix rows.

Note that if the matrix was constructed with only a row-map, then
until FillComplete() is called, this method returns a column-map that
is a copy of the row-map. That 'initial' column-map is replaced with a
computed column- map (that contains the set of column-indices
appearing in each processor's local portion of the matrix) when
FillComplete() is called.

HaveColMap()==true ";

%feature("docstring")  Epetra_CrsMatrix::DomainMap "const Epetra_Map&
Epetra_CrsMatrix::DomainMap() const

Returns the Epetra_Map object associated with the domain of this
matrix operator.

Filled()==true ";

%feature("docstring")  Epetra_CrsMatrix::RangeMap "const Epetra_Map&
Epetra_CrsMatrix::RangeMap() const

Returns the Epetra_Map object associated with the range of this matrix
operator.

Filled()==true ";

%feature("docstring")  Epetra_CrsMatrix::Importer "const
Epetra_Import* Epetra_CrsMatrix::Importer() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")  Epetra_CrsMatrix::Exporter "const
Epetra_Export* Epetra_CrsMatrix::Exporter() const

Returns the Epetra_Export object that contains the export operations
for distributed operations. ";

%feature("docstring")  Epetra_CrsMatrix::Comm "const Epetra_Comm&
Epetra_CrsMatrix::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

/*  Local/Global ID methods  */

%feature("docstring")  Epetra_CrsMatrix::LRID "int
Epetra_CrsMatrix::LRID(int GRID_in) const

Returns the local row index for given global row index, returns -1 if
no local row for this global row. ";

%feature("docstring")  Epetra_CrsMatrix::LRID "int
Epetra_CrsMatrix::LRID(long long GRID_in) const ";

%feature("docstring")  Epetra_CrsMatrix::GRID "int
Epetra_CrsMatrix::GRID(int LRID_in) const

Returns the global row index for give local row index, returns
IndexBase-1 if we don't have this local row. ";

%feature("docstring")  Epetra_CrsMatrix::GRID64 "long long
Epetra_CrsMatrix::GRID64(int LRID_in) const ";

%feature("docstring")  Epetra_CrsMatrix::LCID "int
Epetra_CrsMatrix::LCID(int GCID_in) const

Returns the local column index for given global column index, returns
-1 if no local column for this global column.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsMatrix::LCID "int
Epetra_CrsMatrix::LCID(long long GCID_in) const ";

%feature("docstring")  Epetra_CrsMatrix::GCID "int
Epetra_CrsMatrix::GCID(int LCID_in) const

Returns the global column index for give local column index, returns
IndexBase-1 if we don't have this local column.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsMatrix::GCID64 "long long
Epetra_CrsMatrix::GCID64(int LCID_in) const ";

%feature("docstring")  Epetra_CrsMatrix::MyGRID "bool
Epetra_CrsMatrix::MyGRID(int GRID_in) const

Returns true if the GRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_CrsMatrix::MyGRID "bool
Epetra_CrsMatrix::MyGRID(long long GRID_in) const ";

%feature("docstring")  Epetra_CrsMatrix::MyLRID "bool
Epetra_CrsMatrix::MyLRID(int LRID_in) const

Returns true if the LRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_CrsMatrix::MyGCID "bool
Epetra_CrsMatrix::MyGCID(int GCID_in) const

Returns true if the GCID passed in belongs to the calling processor in
this map, otherwise returns false.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsMatrix::MyGCID "bool
Epetra_CrsMatrix::MyGCID(long long GCID_in) const ";

%feature("docstring")  Epetra_CrsMatrix::MyLCID "bool
Epetra_CrsMatrix::MyLCID(int LCID_in) const

Returns true if the LRID passed in belongs to the calling processor in
this map, otherwise returns false.

HaveColMap()==true (If HaveColMap()==false, returns -1) ";

%feature("docstring")  Epetra_CrsMatrix::MyGlobalRow "bool
Epetra_CrsMatrix::MyGlobalRow(int GID) const

Returns true of GID is owned by the calling processor, otherwise it
returns false. ";

%feature("docstring")  Epetra_CrsMatrix::MyGlobalRow "bool
Epetra_CrsMatrix::MyGlobalRow(long long GID) const ";

/*  I/O Methods  */

%feature("docstring")  Epetra_CrsMatrix::Print "void
Epetra_CrsMatrix::Print(ostream &os) const

Print method. ";

/*  Additional methods required to support the Epetra_Operator
interface  */

%feature("docstring")  Epetra_CrsMatrix::Label "const char*
Epetra_CrsMatrix::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_CrsMatrix::SetUseTranspose "int
Epetra_CrsMatrix::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

UseTranspose:  - (In) If true, multiply by the transpose of operator,
otherwise just use operator.

Always returns 0. ";

%feature("docstring")  Epetra_CrsMatrix::Apply "int
Epetra_CrsMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

X:  - (In) An Epetra_MultiVector of dimension NumVectors to multiply
with matrix.

Y:  -(Out) An Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::ApplyInverse "int
Epetra_CrsMatrix::ApplyInverse(const Epetra_MultiVector &X,
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

X:  - (In) An Epetra_MultiVector of dimension NumVectors to solve for.

Y:  - (Out) An Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::HasNormInf "bool
Epetra_CrsMatrix::HasNormInf() const

Returns true because this class can compute an Inf-norm. ";

%feature("docstring")  Epetra_CrsMatrix::UseTranspose "bool
Epetra_CrsMatrix::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_CrsMatrix::OperatorDomainMap "const
Epetra_Map& Epetra_CrsMatrix::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")  Epetra_CrsMatrix::OperatorRangeMap "const
Epetra_Map& Epetra_CrsMatrix::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";

/*  Additional methods required to implement Epetra_RowMatrix
interface  */

%feature("docstring")  Epetra_CrsMatrix::NumMyRowEntries "int
Epetra_CrsMatrix::NumMyRowEntries(int MyRow, int &NumEntries) const

Return the current number of values stored for the specified local
row.

Similar to NumMyEntries() except NumEntries is returned as an argument
and error checking is done on the input value MyRow.

Parameters:
-----------

MyRow:  - (In) Local row.

NumEntries:  - (Out) Number of nonzero values.

Integer error code, set to 0 if successful.

None.

Unchanged. ";

%feature("docstring")  Epetra_CrsMatrix::Map "const Epetra_BlockMap&
Epetra_CrsMatrix::Map() const

Map() method inherited from Epetra_DistObject. ";

%feature("docstring")  Epetra_CrsMatrix::RowMatrixRowMap "const
Epetra_Map& Epetra_CrsMatrix::RowMatrixRowMap() const

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")  Epetra_CrsMatrix::RowMatrixColMap "const
Epetra_Map& Epetra_CrsMatrix::RowMatrixColMap() const

Returns the Epetra_Map object associated with columns of this matrix.
";

%feature("docstring")  Epetra_CrsMatrix::RowMatrixImporter "const
Epetra_Import* Epetra_CrsMatrix::RowMatrixImporter() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

/*  Inlined Operator Methods  */

/*  Expert-only methods:  These methods are intended for experts only
and have some risk of changing in the future, since they rely on
underlying data structure assumptions  */

%feature("docstring")  Epetra_CrsMatrix::ExtractCrsDataPointers "int
Epetra_CrsMatrix::ExtractCrsDataPointers(int *&IndexOffset, int
*&Indices, double *&Values_in) const

Returns internal data pointers associated with Crs matrix format.

Returns data pointers to facilitate optimized code within external
packages.

Parameters:
-----------

IndexOffset:  - (Out) Extracted array of indices into Values[] and
Indices[]. Local row k is stored in
Values[IndexOffset[k]:IndexOffset[k+1]-1] and
Indices[IndexOffset[k]:IndexOffset[k+1]-1].

Values:  - (Out) Extracted values for all local rows.

Indices:  - (Out) Extracted local column indices for the corresponding
values.

Integer error code, set to 0 if successful. Returns -1 if FillComplete
has not been performed or Storage has not been Optimized.

WARNING:  This method is intended for expert only, its use may require
user code modifications in future versions of Epetra. ";

%feature("docstring")
Epetra_CrsMatrix::SortGhostsAssociatedWithEachProcessor "int
Epetra_CrsMatrix::SortGhostsAssociatedWithEachProcessor(bool Flag)

Forces FillComplete() to locally order ghostnodes associated with each
remote processor in ascending order.

To be compliant with AztecOO, FillComplete() already locally orders
ghostnodes such that information received from processor k has a lower
local numbering than information received from processor j if k is
less than j. SortGhostsAssociatedWithEachProcessor(True) further
forces FillComplete() to locally number all ghostnodes received from
processor k in ascending order. That is, the local numbering of b is
less than c if the global numbering of b is less than c and if both b
and c are owned by the same processor. This is done to be compliant
with some limited block features within ML. In particular, some ML
features require that a block structure of the matrix be maintained
even within the ghost variables. Always returns 0. ";

/*  Deprecated methods:  These methods still work, but will be removed
in a future version  */

%feature("docstring")  Epetra_CrsMatrix::ImportMap "const Epetra_Map&
Epetra_CrsMatrix::ImportMap() const

Use ColMap() instead. ";

%feature("docstring")  Epetra_CrsMatrix::TransformToLocal "int
Epetra_CrsMatrix::TransformToLocal()

Use FillComplete() instead. ";

%feature("docstring")  Epetra_CrsMatrix::TransformToLocal "int
Epetra_CrsMatrix::TransformToLocal(const Epetra_Map *DomainMap, const
Epetra_Map *RangeMap)

Use FillComplete(const Epetra_Map& DomainMap, const Epetra_Map&
RangeMap) instead. ";


// File: classEpetra__CrsSingletonFilter.xml
%feature("docstring") Epetra_CrsSingletonFilter "

Epetra_CrsSingletonFilter: A class for explicitly eliminating matrix
rows and columns.

The Epetra_CrsSingletonFilter class takes an existing
Epetra_LinearProblem object, analyzes it structure and explicitly
eliminates singleton rows and columns from the matrix and
appropriately modifies the RHS and LHS of the linear problem. The
result of this process is a reduced system of equations that is itself
an Epetra_LinearProblem object. The reduced system can then be solved
using any solver that is understands an Epetra_LinearProblem. The
solution for the full system is obtained by calling
ComputeFullSolution().

Singleton rows are defined to be rows that have a single nonzero entry
in the matrix. The equation associated with this row can be explicitly
eliminated because it involved only one variable. For example if row i
has a single nonzero value in column j, call it A(i,j), we can
explicitly solve for x(j) = b(i)/A(i,j), where b(i) is the ith entry
of the RHS and x(j) is the jth entry of the LHS.

Singleton columns are defined to be columns that have a single nonzero
entry in the matrix. The variable associated with this column is fully
dependent, meaning that the solution for all other variables does not
depend on it. If this entry is A(i,j) then the ith row and jth column
can be removed from the system and x(j) can be solved after the
solution for all other variables is determined.

By removing singleton rows and columns, we can often produce a reduced
system that is smaller and far less dense, and in general having
better numerical properties.

The basic procedure for using this class is as follows: Construct full
problem: Construct and Epetra_LinearProblem containing the \"full\"
matrix, RHS and LHS. This is done outside of Epetra_CrsSingletonFilter
class. Presumably, you have some reason to believe that this system
may contain singletons.

Construct an Epetra_CrsSingletonFilter instance: Constructor needs no
arguments.

Analyze matrix: Invoke the Analyze() method, passing in the
Epetra_RowMatrix object from your full linear problem mentioned in the
first step above.

Go/No Go decision to construct reduced problem: Query the results of
the Analyze method using the SingletonsDetected() method. This method
returns \"true\" if there were singletons found in the matrix. You can
also query any of the other methods in the Filter Statistics section
to determine if you want to proceed with the construction of the
reduced system.

Construct reduced problem: If, in the previous step, you determine
that you want to proceed with the construction of the reduced problem,
you should next call the ConstructReducedProblem() method, passing in
the full linear problem object from the first step. This method will
use the information from the Analyze() method to construct a reduce
problem that has explicitly eliminated the singleton rows, solved for
the corresponding LHS values and updated the RHS. This step will also
remove singleton columns from the reduced system. Once the solution of
the reduced problem is is computed (via any solver that understands an
Epetra_LinearProblem), you should call the ComputeFullSolution()
method to compute the LHS values assocaited with the singleton
columns.

Solve reduced problem: Obtain a pointer to the reduced problem using
the ReducedProblem() method. Using the solver of your choice, solve
the reduced system.

Compute solution to full problem: Once the solution the reduced
problem is determined, the ComputeFullSolution() method will place the
reduced solution values into the appropriate locations of the full
solution LHS and then compute the values associated with column
singletons. At this point, you have a complete solution to the
original full problem.

Solve a subsequent full problem that differs from the original problem
only in values: It is often the case that the structure of a problem
will be the same for a sequence of linear problems. In this case, the
UpdateReducedProblem() method can be useful. After going through the
above process one time, if you have a linear problem that is
structural identical to the previous problem, you can minimize memory
and time costs by using the UpdateReducedProblem() method, passing in
the subsequent problem. Once you have called the
UpdateReducedProblem() method, you can then solve the reduce problem
problem as you wish, and then compute the full solution as before. The
pointer generated by ReducedProblem() will not change when
UpdateReducedProblem() is called.

C++ includes: Epetra_CrsSingletonFilter.h ";

/*  Constructors/Destructor  */

%feature("docstring")
Epetra_CrsSingletonFilter::Epetra_CrsSingletonFilter "Epetra_CrsSingletonFilter::Epetra_CrsSingletonFilter()

Epetra_CrsSingletonFilter default constructor. ";

%feature("docstring")
Epetra_CrsSingletonFilter::~Epetra_CrsSingletonFilter "Epetra_CrsSingletonFilter::~Epetra_CrsSingletonFilter()

Epetra_CrsSingletonFilter Destructor. ";

/*  Analyze methods  */

%feature("docstring")  Epetra_CrsSingletonFilter::Analyze "int
Epetra_CrsSingletonFilter::Analyze(Epetra_RowMatrix *FullMatrix)

Analyze the input matrix, removing row/column pairs that have
singletons.

Analyzes the user's input matrix to determine rows and columns that
should be explicitly eliminated to create the reduced system. Look for
rows and columns that have single entries. These rows/columns can
easily be removed from the problem. The results of calling this method
are two MapColoring objects accessible via RowMapColors() and
ColMapColors() accessor methods. All rows/columns that would be
eliminated in the reduced system have a color of 1 in the
corresponding RowMapColors/ColMapColors object. All kept rows/cols
have a color of 0. ";

%feature("docstring")  Epetra_CrsSingletonFilter::SingletonsDetected "bool Epetra_CrsSingletonFilter::SingletonsDetected() const

Returns true if singletons were detected in this matrix (must be
called after Analyze() to be effective). ";

/*  Reduce methods  */

%feature("docstring")
Epetra_CrsSingletonFilter::ConstructReducedProblem "int
Epetra_CrsSingletonFilter::ConstructReducedProblem(Epetra_LinearProblem
*Problem)

Return a reduced linear problem based on results of Analyze().

Creates a new Epetra_LinearProblem object based on the results of the
Analyze phase. A pointer to the reduced problem is obtained via a call
to ReducedProblem().

Error code, set to 0 if no error. ";

%feature("docstring")  Epetra_CrsSingletonFilter::UpdateReducedProblem
"int
Epetra_CrsSingletonFilter::UpdateReducedProblem(Epetra_LinearProblem
*Problem)

Update a reduced linear problem using new values.

Updates an existing Epetra_LinearProblem object using new matrix, LHS
and RHS values. The matrix structure must be identical to the matrix
that was used to construct the original reduced problem.

Error code, set to 0 if no error. ";

/*  Methods to construct Full System Solution  */

%feature("docstring")  Epetra_CrsSingletonFilter::ComputeFullSolution
"int Epetra_CrsSingletonFilter::ComputeFullSolution()

Compute a solution for the full problem using the solution of the
reduced problem, put in LHS of FullProblem().

After solving the reduced linear system, this method can be called to
compute the solution to the original problem, assuming the solution
for the reduced system is valid. The solution of the unreduced,
original problem will be in the LHS of the original
Epetra_LinearProblem. ";

/*  Filter Statistics  */

%feature("docstring")  Epetra_CrsSingletonFilter::NumRowSingletons "int Epetra_CrsSingletonFilter::NumRowSingletons() const

Return number of rows that contain a single entry, returns -1 if
Analysis not performed yet. ";

%feature("docstring")  Epetra_CrsSingletonFilter::NumColSingletons "int Epetra_CrsSingletonFilter::NumColSingletons() const

Return number of columns that contain a single entry that are not
associated with singleton row, returns -1 if Analysis not performed
yet. ";

%feature("docstring")  Epetra_CrsSingletonFilter::NumSingletons "int
Epetra_CrsSingletonFilter::NumSingletons() const

Return total number of singletons detected, returns -1 if Analysis not
performed yet.

Return total number of singletons detected across all processors. This
method will not return a valid result until after the Analyze() method
is called. The dimension of the reduced system can be computed by
subtracting this number from dimension of full system. WARNING:  This
method returns -1 if Analyze() method has not been called. ";

%feature("docstring")  Epetra_CrsSingletonFilter::RatioOfDimensions "double Epetra_CrsSingletonFilter::RatioOfDimensions() const

Returns ratio of reduced system to full system dimensions, returns
-1.0 if reduced problem not constructed. ";

%feature("docstring")  Epetra_CrsSingletonFilter::RatioOfNonzeros "double Epetra_CrsSingletonFilter::RatioOfNonzeros() const

Returns ratio of reduced system to full system nonzero count, returns
-1.0 if reduced problem not constructed. ";

/*  Attribute Access Methods  */

%feature("docstring")  Epetra_CrsSingletonFilter::FullProblem "Epetra_LinearProblem* Epetra_CrsSingletonFilter::FullProblem() const

Returns pointer to the original unreduced Epetra_LinearProblem. ";

%feature("docstring")  Epetra_CrsSingletonFilter::ReducedProblem "Epetra_LinearProblem* Epetra_CrsSingletonFilter::ReducedProblem()
const

Returns pointer to the derived reduced Epetra_LinearProblem. ";

%feature("docstring")  Epetra_CrsSingletonFilter::FullMatrix "Epetra_RowMatrix* Epetra_CrsSingletonFilter::FullMatrix() const

Returns pointer to Epetra_CrsMatrix from full problem. ";

%feature("docstring")  Epetra_CrsSingletonFilter::ReducedMatrix "Epetra_CrsMatrix* Epetra_CrsSingletonFilter::ReducedMatrix() const

Returns pointer to Epetra_CrsMatrix from full problem. ";

%feature("docstring")  Epetra_CrsSingletonFilter::RowMapColors "Epetra_MapColoring* Epetra_CrsSingletonFilter::RowMapColors() const

Returns pointer to Epetra_MapColoring object: color 0 rows are part of
reduced system. ";

%feature("docstring")  Epetra_CrsSingletonFilter::ColMapColors "Epetra_MapColoring* Epetra_CrsSingletonFilter::ColMapColors() const

Returns pointer to Epetra_MapColoring object: color 0 columns are part
of reduced system. ";

%feature("docstring")  Epetra_CrsSingletonFilter::ReducedMatrixRowMap
"Epetra_Map* Epetra_CrsSingletonFilter::ReducedMatrixRowMap() const

Returns pointer to Epetra_Map describing the reduced system row
distribution. ";

%feature("docstring")  Epetra_CrsSingletonFilter::ReducedMatrixColMap
"Epetra_Map* Epetra_CrsSingletonFilter::ReducedMatrixColMap() const

Returns pointer to Epetra_Map describing the reduced system column
distribution. ";

%feature("docstring")
Epetra_CrsSingletonFilter::ReducedMatrixDomainMap "Epetra_Map*
Epetra_CrsSingletonFilter::ReducedMatrixDomainMap() const

Returns pointer to Epetra_Map describing the domain map for the
reduced system. ";

%feature("docstring")
Epetra_CrsSingletonFilter::ReducedMatrixRangeMap "Epetra_Map*
Epetra_CrsSingletonFilter::ReducedMatrixRangeMap() const

Returns pointer to Epetra_Map describing the range map for the reduced
system. ";


// File: classEpetra__Data.xml
%feature("docstring") Epetra_Data "

Epetra_Data: The Epetra Base Data Class.

The Epetra_Data class is a base class for all Epetra Data Classes. It
provides a mechanism so that one data object can be shared by multiple
class instances. However, it is meant only to be used internally by
another Epetra class. It does not provide smart pointer like
capabilities. Incrementing and decrementing the reference count, and
deleting the data class instance (if necessary), are duties of the
Epetra class utilizing Epetra_Data.

All of Epetra_Data's methods are protected. This is because
Epetra_Data should never be used directly. Rather, a class that
derives from Epetra_Data should be used instead. For example,
Epetra_MpiCommData or Epetra_BlockMapData.

DEVELOPER NOTES: (1) Any class that inherits from Epetra_Data may need
to define an assignment operator, if it adds pointers. Epetra_Data
doesn't have any, and so the default (compiler-generated) assignment
operator is good enough. (2) The behavior of a derived class is left
up to the implementer(s) of that class. As such, it cannot be assumed
that just because a class inherits from Epetra_Data, that it supports
copy construction or assignment, or that it will perform as expected.

C++ includes: Epetra_Data.h ";

/*  Constructor/Destructor Methods  */

/*  Reference-Counting Methods  */


// File: classEpetra__Directory.xml
%feature("docstring") Epetra_Directory "

Epetra_Directory: This class is a pure virtual class whose interface
allows Epetra_Map and Epetr_BlockMap objects to reference non-local
elements.

For Epetra_BlockMap objects, a Epetra_Directory object must be created
by a call to the Epetra_Comm CreateDirectory method. The Directory is
needed to allow referencing of non-local elements.

C++ includes: Epetra_Directory.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_Directory::~Epetra_Directory "virtual
Epetra_Directory::~Epetra_Directory()

Epetra_Directory destructor. ";

/*  Query method  */

%feature("docstring")  Epetra_Directory::GetDirectoryEntries "virtual
int Epetra_Directory::GetDirectoryEntries(const Epetra_BlockMap &Map,
const int NumEntries, const int *GlobalEntries, int *Procs, int
*LocalEntries, int *EntrySizes, bool high_rank_sharing_procs=false)
const =0

GetDirectoryEntries : Returns proc and local id info for non-local map
entries.

Given a list of Global Entry IDs, this function returns the list of
processor IDs and local IDs on the owning processor that correspond to
the list of entries. If LocalEntries is 0, then local IDs are not
returned. If EntrySizes is nonzero, it will contain a list of
corresponding element sizes for the requested global entries.

Parameters:
-----------

In:  NumEntries - Number of Global IDs being passed in.

In:  GlobalEntries - List of Global IDs being passed in.

InOut:  Procs - User allocated array of length at least NumEntries. On
return contains list of processors owning the Global IDs in question.

InOut:  LocalEntries - User allocated array of length at least
NumEntries. On return contains the local ID of the global on the
owning processor. If LocalEntries is zero, no local ID information is
returned.

InOut:  EntrySizes - User allocated array of length at least
NumEntries. On return contains the size of the object associated with
this global ID. If LocalEntries is zero, no size information is
returned.

In:  high_rank_sharing_procs Optional argument, defaults to true. If
any GIDs appear on multiple processors (referred to as \"sharing
procs\"), this specifies whether the lowest-rank proc or the highest-
rank proc is chosen as the \"owner\".

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_Directory::GetDirectoryEntries "virtual
int Epetra_Directory::GetDirectoryEntries(const Epetra_BlockMap &Map,
const int NumEntries, const long long *GlobalEntries, int *Procs, int
*LocalEntries, int *EntrySizes, bool high_rank_sharing_procs=false)
const =0 ";

%feature("docstring")  Epetra_Directory::GIDsAllUniquelyOwned "virtual bool Epetra_Directory::GIDsAllUniquelyOwned() const =0

GIDsAllUniquelyOwned: returns true if all GIDs appear on just one
processor.

If any GIDs are owned by multiple processors, returns false. ";


// File: classEpetra__DistObject.xml
%feature("docstring") Epetra_DistObject "

Epetra_DistObject: A class for constructing and using dense multi-
vectors, vectors and matrices in parallel.

The Epetra_DistObject is a base class for all Epetra distributed
global objects. It provides the basic mechanisms and interface
specifications for importing and exporting operations using
Epetra_Import and Epetra_Export objects.

Distributed Global vs. Replicated Local.

Distributed Global objects - In most instances, a distributed object
will be partitioned across multiple memory images associated with
multiple processors. In this case, there is a unique copy of each
element and elements are spread across all processors specified by the
Epetra_Comm communicator.

Replicated Local Objects - Some algorithms use objects that are too
small to be distributed across all processors, the Hessenberg matrix
in a GMRES computation. In other cases, such as with block iterative
methods, block dot product functions produce small dense matrices that
are required by all processors. Replicated local objectss handle these
types of situation.

C++ includes: Epetra_DistObject.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_DistObject::Epetra_DistObject "Epetra_DistObject::Epetra_DistObject(const Epetra_BlockMap &Map)

Basic Epetra_DistObject constuctor.

Creates a Epetra_DistObject object.

Parameters:
-----------

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

WARNING:  Note that, because Epetra_LocalMap derives from Epetra_Map
and Epetra_Map derives from Epetra_BlockMap, this constructor works
for all three types of Epetra map classes.

Pointer to a Epetra_DistObject. ";

%feature("docstring")  Epetra_DistObject::Epetra_DistObject "Epetra_DistObject::Epetra_DistObject(const Epetra_BlockMap &Map, const
char *const Label)

Creates a Epetra_DistObject object.

Parameters:
-----------

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

WARNING:  Note that, because Epetra_LocalMap derives from Epetra_Map
and Epetra_Map derives from Epetra_BlockMap, this constructor works
for all three types of Epetra map classes.

Parameters:
-----------

In:  Label - An identifier for this object. By default, set to the
name of the object class.

Pointer to a Epetra_DistObject. ";

%feature("docstring")  Epetra_DistObject::Epetra_DistObject "Epetra_DistObject::Epetra_DistObject(const Epetra_DistObject &Source)

Epetra_DistObject copy constructor. ";

%feature("docstring")  Epetra_DistObject::~Epetra_DistObject "Epetra_DistObject::~Epetra_DistObject()

Epetra_DistObject destructor. ";

/*  Import/Export Methods  */

%feature("docstring")  Epetra_DistObject::Import "int
Epetra_DistObject::Import(const Epetra_SrcDistObject &A, const
Epetra_Import &Importer, Epetra_CombineMode CombineMode, const
Epetra_OffsetIndex *Indexor=0)

Imports an Epetra_DistObject using the Epetra_Import object.

Parameters:
-----------

In:  Source - Distributed object that will be imported into the
\"\\\\e this\" object.

In:  Importer - A Epetra_Import object specifying the communication
required.

In:  CombineMode - A Epetra_CombineMode enumerated type specifying how
results should be combined on the receiving processor.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_DistObject::Import "int
Epetra_DistObject::Import(const Epetra_SrcDistObject &A, const
Epetra_Export &Exporter, Epetra_CombineMode CombineMode, const
Epetra_OffsetIndex *Indexor=0)

Imports an Epetra_DistObject using the Epetra_Export object.

Parameters:
-----------

In:  Source - Distributed object that will be imported into the
\"\\\\e this\" object.

In:  Exporter - A Epetra_Export object specifying the communication
required.

In:  CombineMode - A Epetra_CombineMode enumerated type specifying how
results should be combined on the receiving processor.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_DistObject::Export "int
Epetra_DistObject::Export(const Epetra_SrcDistObject &A, const
Epetra_Import &Importer, Epetra_CombineMode CombineMode, const
Epetra_OffsetIndex *Indexor=0)

Exports an Epetra_DistObject using the Epetra_Import object.

Parameters:
-----------

In:  Source - Distributed object that will be exported to the \"\\\\e
this\" object.

In:  Importer - A Epetra_Import object specifying the communication
required.

In:  CombineMode - A Epetra_CombineMode enumerated type specifying how
results should be combined on the receiving processor.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_DistObject::Export "int
Epetra_DistObject::Export(const Epetra_SrcDistObject &A, const
Epetra_Export &Exporter, Epetra_CombineMode CombineMode, const
Epetra_OffsetIndex *Indexor=0)

Exports an Epetra_DistObject using the Epetra_Export object.

Parameters:
-----------

In:  Source - Distributed object that will be exported to the \"\\\\e
this\" multivector.

In:  Exporter - A Epetra_Export object specifying the communication
required.

In:  CombineMode - A Epetra_CombineMode enumerated type specifying how
results should be combined on the receiving processor.

Integer error code, set to 0 if successful. ";

/*  Attribute accessor methods  */

%feature("docstring")  Epetra_DistObject::Map "const Epetra_BlockMap&
Epetra_DistObject::Map() const

Returns the address of the Epetra_BlockMap for this multi-vector. ";

%feature("docstring")  Epetra_DistObject::Comm "const Epetra_Comm&
Epetra_DistObject::Comm() const

Returns the address of the Epetra_Comm for this multi-vector. ";

%feature("docstring")  Epetra_DistObject::DistributedGlobal "bool
Epetra_DistObject::DistributedGlobal() const

Returns true if this multi-vector is distributed global, i.e., not
local replicated. ";

/*  Miscellaneous  */

%feature("docstring")  Epetra_DistObject::Print "void
Epetra_DistObject::Print(ostream &os) const

Print method. ";

/*  Internal utilities  */

/*  Virtual methods to be implemented by derived class  */


// File: classEpetra__Distributor.xml
%feature("docstring") Epetra_Distributor "

Epetra_Distributor: The Epetra Gather/Scatter Setup Base Class.

The Epetra_Distributor class is an interface that encapsulates the
general information and services needed for other Epetra classes to
perform gather/scatter operations on a parallel computer. An
Epetra_Distributor object is actually produced by calling a method in
the Epetra_Comm class.

Epetra_Distributor has default implementations, via
Epetra_SerialDistributor and Epetra_MpiDistributor, for both serial
execution and MPI distributed memory execution. It is meant to
insulate the user from the specifics of communication that are not
required for normal manipulation of linear algebra objects..

C++ includes: Epetra_Distributor.h ";

/*  Constructor and Destructor  */

%feature("docstring")  Epetra_Distributor::Clone "virtual
Epetra_Distributor* Epetra_Distributor::Clone()=0

Epetra_Distributor clone constructor. ";

%feature("docstring")  Epetra_Distributor::~Epetra_Distributor "virtual Epetra_Distributor::~Epetra_Distributor()

Epetra_Distributor Destructor. ";

/*  Gather/Scatter Constructors  */

%feature("docstring")  Epetra_Distributor::CreateFromSends "virtual
int Epetra_Distributor::CreateFromSends(const int &NumExportIDs, const
int *ExportPIDs, bool Deterministic, int &NumRemoteIDs)=0

Create Distributor object using list of process IDs to which we
export.

Take a list of Process IDs and construct a plan for efficiently
scattering to these processes. Return the number of IDs being sent to
me.

Parameters:
-----------

NumExportIDs:  In Number of IDs that need to be sent from this
processor.

ExportPIDs:  In List of processors that will get the exported IDs.

Deterministic:  In No op.

NumRemoteIDs:  Out Number of IDs this processor will be receiving. ";

%feature("docstring")  Epetra_Distributor::CreateFromRecvs "virtual
int Epetra_Distributor::CreateFromRecvs(const int &NumRemoteIDs, const
int *RemoteGIDs, const int *RemotePIDs, bool Deterministic, int
&NumExportIDs, int *&ExportGIDs, int *&ExportPIDs)=0

Create Distributor object using list of Remote global IDs and
corresponding PIDs.

Take a list of global IDs and construct a plan for efficiently
scattering to these processes. Return the number and list of IDs being
sent by me.

Parameters:
-----------

NumRemoteIDs:  In Number of IDs this processor will be receiving.

RemoteGIDs:  In List of IDs that this processor wants.

RemotePIDs:  In List of processors that will send the remote IDs.

Deterministic:  In No op.

NumExportIDs:  Out Number of IDs that need to be sent from this
processor.

ExportPIDs:  Out List of processors that will get the exported IDs. ";

%feature("docstring")  Epetra_Distributor::CreateFromRecvs "virtual
int Epetra_Distributor::CreateFromRecvs(const int &NumRemoteIDs, const
long long *RemoteGIDs, const int *RemotePIDs, bool Deterministic, int
&NumExportIDs, long long *&ExportGIDs, int *&ExportPIDs)=0 ";

/*  Execute Gather/Scatter Operations (Constant size objects)  */

%feature("docstring")  Epetra_Distributor::Do "virtual int
Epetra_Distributor::Do(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)=0

Execute plan on buffer of export objects in a single step. ";

%feature("docstring")  Epetra_Distributor::DoReverse "virtual int
Epetra_Distributor::DoReverse(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)=0

Execute reverse of plan on buffer of export objects in a single step.
";

%feature("docstring")  Epetra_Distributor::DoPosts "virtual int
Epetra_Distributor::DoPosts(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)=0

Post buffer of export objects (can do other local work before
executing Waits) ";

%feature("docstring")  Epetra_Distributor::DoWaits "virtual int
Epetra_Distributor::DoWaits()=0

Wait on a set of posts. ";

%feature("docstring")  Epetra_Distributor::DoReversePosts "virtual
int Epetra_Distributor::DoReversePosts(char *export_objs, int
obj_size, int &len_import_objs, char *&import_objs)=0

Do reverse post of buffer of export objects (can do other local work
before executing Waits) ";

%feature("docstring")  Epetra_Distributor::DoReverseWaits "virtual
int Epetra_Distributor::DoReverseWaits()=0

Wait on a reverse set of posts. ";

/*  Execute Gather/Scatter Operations (Non-constant size objects)  */

%feature("docstring")  Epetra_Distributor::Do "virtual int
Epetra_Distributor::Do(char *export_objs, int obj_size, int *&sizes,
int &len_import_objs, char *&import_objs)=0

Execute plan on buffer of export objects in a single step (object size
may vary) ";

%feature("docstring")  Epetra_Distributor::DoReverse "virtual int
Epetra_Distributor::DoReverse(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)=0

Execute reverse of plan on buffer of export objects in a single step
(object size may vary) ";

%feature("docstring")  Epetra_Distributor::DoPosts "virtual int
Epetra_Distributor::DoPosts(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)=0

Post buffer of export objects (can do other local work before
executing Waits) ";

%feature("docstring")  Epetra_Distributor::DoReversePosts "virtual
int Epetra_Distributor::DoReversePosts(char *export_objs, int
obj_size, int *&sizes, int &len_import_objs, char *&import_objs)=0

Do reverse post of buffer of export objects (can do other local work
before executing Waits) ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_Distributor::Print "virtual void
Epetra_Distributor::Print(ostream &os) const =0 ";


// File: classEpetra__Export.xml
%feature("docstring") Epetra_Export "

Epetra_Export: This class builds an export object for efficient
exporting of off- processor elements.

Epetra_Export is used to construct a communication plan that can be
called repeatedly by computational classes such the Epetra matrix,
vector and multivector classes to efficiently send data to a target
processor.

This class currently has one constructor, taking two Epetra_Map or
Epetra_BlockMap objects. The first map specifies the global IDs that
are owned by the calling processor. The second map specifies the
global IDs of elements that we want to export to later.

C++ includes: Epetra_Export.h ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_Export::Print "void
Epetra_Export::Print(ostream &os) const

Print object to an output stream Print method ";

%feature("docstring")  Epetra_Export::Epetra_Export "Epetra_Export::Epetra_Export(const Epetra_BlockMap &SourceMap, const
Epetra_BlockMap &TargetMap)

Constructs a Epetra_Export object from the source and target maps.

This constructor builds an Epetra_Export object by comparing the GID
lists of the source and target maps.

Parameters:
-----------

SourceMap:  (In) Map containing the GIDs from which data should be
exported from each processor to the target map whenever an export
operation is performed using this exporter.

TargetMap:  (In) Map containing the GIDs that should be used for
exporting data.

WARNING:  Note that the TargetMap must have GIDs uniquely owned, each
GID of the target map can occur only once.  Builds an export object
that will transfer objects built with SourceMap to objects built with
TargetMap.

A Epetra_Export object categorizes the elements of the target map into
three sets as follows: All elements in the target map that have the
same GID as the corresponding element of the source map, starting with
the first element in the target map, going up to the first element
that is different from the source map. The number of these IDs is
returned by NumSameIDs().

All elements that are local to the processor, but are not part of the
first set of elements. These elements have GIDs that are owned by the
calling processor, but at least the first element of this list is
permuted. Even if subsequent elements are not permuted, they are
included in this list. The number of permuted elements is returned by
NumPermutedIDs(). The list of elements (local IDs) in the source map
that are permuted can be found in the list PermuteFromLIDs(). The list
of elements (local IDs) in the target map that are the new locations
of the source elements can be found in the list PermuteToLIDs().

All remaining elements of the target map correspond to global IDs that
are owned by remote processors. The number of these elements is
returned by NumRemoteIDs() and the list of these is returned by
RemoteLIDs().

Given the above information, the Epetra_Export constructor builds a
list of elements that must be communicated to other processors as a
result of export requests. The number of exported elements (where
multiple sends of the same element to different processors is counted)
is returned by NumExportIDs(). The local IDs to be sent are returned
by the list ExportLIDs(). The processors to which each of the elements
will be sent in returned in a list of the same length by ExportPIDs().

The total number of elements that will be sent by the calling
processor is returned by NumSend(). The total number of elements that
will be received is returned by NumRecv().

The following example illustrates the basic concepts.

Assume we have 3 processors and 9 global elements with each processor
owning 3 elements as follows  PE 0 Elements |  PE 1 Elements  |  PE 2
Elements     0  1  2 3  4  5           6  7  8

The above layout essentially defines the target map argument of the
export object.

This could correspond to a 9-entry forcing vector with the first three
entries on PE 0, and so on. Suppose that the entries of this forcing
vector are computed by integrating over linear \"hat\" functions:

^  ^  ^  ^  ^  ^  ^  ^  ^  \\\\/ \\\\/ \\\\/ \\\\/ \\\\/ \\\\/ \\\\/
\\\\/   /\\\\ /\\\\ /\\\\ /\\\\ /\\\\ /\\\\ /\\\\ /\\\\
+--+--+--+--+--+--+--+--+ 0  1  2  3  4  5  6  7  8

In this case, PE 0 will make contributions to entries 0 through 3, PE
1 will make contributions to entries 2 through 6 and PE 2 will make
contributions to entries 5 through 8. A convenient way to compute
these contributions is to create a forcing vector with replicated
entries for the shared contributions. Specifically the following
SourceMap works for this scenario:

PE 0 Elements    |  PE 1 Elements    |  PE 2 Elements      0  1 2  3
2  3  4  5  6        5  6  7  8

A vector constructed using this SourceMap can be used to collect each
processor's contributions to the forcing vector. Note that the
ordering of the elements on each processor is not unique, but has been
chosen for illustration.

With these two maps passed into the Epetra_Export constructor, we get
the following attribute definitions:

On PE 0:

NumSameIDs      = 3  NumPermuteIDs   = 0 PermuteToLIDs   = 0
PermuteFromLIDs = 0  NumRemoteIDs    = 1 RemoteLIDs      = [2]
NumExportIDs    = 1 ExportLIDs      = [3] ExportPIDs      = [1]
NumSend         = 1 NumRecv         = 1

On PE 1:

NumSameIDs      = 0  NumPermuteIDs   = 3 PermuteToLIDs   = [0, 1, 2]
PermuteFromLIDs = [1, 2, 3]  NumRemoteIDs    = 2 RemoteLIDs      = [0,
2]  NumExportIDs    = 2 ExportLIDs      = [0, 4] ExportPIDs      = [0,
2]  NumSend         = 2 NumRecv         = 2

On PE 2:

NumSameIDs      = 0  NumPermuteIDs   = 3 PermuteToLIDs   = [0, 1, 2]
PermuteFromLIDs = [1, 2, 3]  NumRemoteIDs    = 1 RemoteLIDs      = [0]
NumExportIDs    = 1 ExportLIDs      = [0] ExportPIDs      = [1]
NumSend         = 1 NumRecv         = 1

Using Epetra_Export Objects

Once a Epetra_Export object has been constructed, it can be used by
any of the Epetra classes that support distributed global objects,
namely Epetra_Vector, Epetra_MultiVector, Epetra_CrsGraph,
Epetra_CrsMatrix and Epetra_VbrMatrix. All of these classes have
Export and Export methods that will fill new objects whose
distribution is described by the target map, taking elements from the
source object whose distribution is described by the source map.
Details of usage for each class is given in the appropriate class
documentation.

In the above example, if x_integrate is constructed using the
SourceMap and then filled with local contributions, and x_force is
constructed using the target map, the following operation will fill
x_force with the combined results of x_integrate:
x_force.Export(x_integrate, exporter, Add); The third argument above
tells the export operation to add results that come from multiple
processors for the same GID.

Epetra_Export objects can also be used by Import operations to perform
the reverse operation. For example, if x_force in the above example
had boundary conditions that should be sent to processors that share a
boundary element, the following operation would send replicated values
to x_integrate: x_integrate.Import(x_force, exporter, Insert); At the
end of this operation, x_integrate would have replicated values from
x_force of entries 2 and 3 on PEs 0 and 1, and entries 5 and 6 on PEs
1 and 2. ";

%feature("docstring")  Epetra_Export::Epetra_Export "Epetra_Export::Epetra_Export(const Epetra_Export &Exporter)

Epetra_Export copy constructor. ";

%feature("docstring")  Epetra_Export::~Epetra_Export "Epetra_Export::~Epetra_Export(void)

Epetra_Export destructor. ";

%feature("docstring")  Epetra_Export::NumSameIDs "int
Epetra_Export::NumSameIDs() const

Returns the number of elements that are identical between the source
and target maps, up to the first different ID. ";

%feature("docstring")  Epetra_Export::NumPermuteIDs "int
Epetra_Export::NumPermuteIDs() const

Returns the number of elements that are local to the calling
processor, but not part of the first NumSameIDs() elements. ";

%feature("docstring")  Epetra_Export::PermuteFromLIDs "int*
Epetra_Export::PermuteFromLIDs() const

List of elements in the source map that are permuted. ";

%feature("docstring")  Epetra_Export::PermuteToLIDs "int*
Epetra_Export::PermuteToLIDs() const

List of elements in the target map that are permuted. ";

%feature("docstring")  Epetra_Export::NumRemoteIDs "int
Epetra_Export::NumRemoteIDs() const

Returns the number of elements that are not on the calling processor.
";

%feature("docstring")  Epetra_Export::RemoteLIDs "int*
Epetra_Export::RemoteLIDs() const

List of elements in the target map that are coming from other
processors. ";

%feature("docstring")  Epetra_Export::NumExportIDs "int
Epetra_Export::NumExportIDs() const

Returns the number of elements that must be sent by the calling
processor to other processors. ";

%feature("docstring")  Epetra_Export::ExportLIDs "int*
Epetra_Export::ExportLIDs() const

List of elements that will be sent to other processors. ";

%feature("docstring")  Epetra_Export::ExportPIDs "int*
Epetra_Export::ExportPIDs() const

List of processors to which elements will be sent, ExportLIDs() [i]
will be sent to processor ExportPIDs() [i]. ";

%feature("docstring")  Epetra_Export::NumSend "int
Epetra_Export::NumSend() const

Total number of elements to be sent. ";

%feature("docstring")  Epetra_Export::NumRecv "int
Epetra_Export::NumRecv() const

Total number of elements to be received. ";

%feature("docstring")  Epetra_Export::SourceMap "const
Epetra_BlockMap& Epetra_Export::SourceMap() const

Returns the SourceMap used to construct this exporter. ";

%feature("docstring")  Epetra_Export::TargetMap "const
Epetra_BlockMap& Epetra_Export::TargetMap() const

Returns the TargetMap used to construct this exporter. ";

%feature("docstring")  Epetra_Export::Distributor "Epetra_Distributor& Epetra_Export::Distributor() const ";


// File: classEpetra__FECrsGraph.xml
%feature("docstring") Epetra_FECrsGraph "

Epetra Finite-Element CrsGraph. This class provides the ability to
insert indices into a matrix-graph, where the indices represent dense
submatrices such as element-stiffnesses that might arise from a
finite-element application.

In a parallel setting, indices may be submitted on the local processor
for rows that do not reside in the local portion of the row-map. After
all indices have been submitted, the GlobalAssemble method gathers all
non-local graph rows to the appropriate 'owning' processors (an owning
processor is a processor which has the row in its row-map).

C++ includes: Epetra_FECrsGraph.h ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int *NumIndicesPerRow, bool
ignoreNonLocalEntries=false, bool buildNonlocalGraph=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int NumIndicesPerRow, bool
ignoreNonLocalEntries=false, bool buildNonlocalGraph=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
*NumIndicesPerRow, bool ignoreNonLocalEntries=false, bool
buildNonlocalGraph=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
NumIndicesPerRow, bool ignoreNonLocalEntries=false, bool
buildNonlocalGraph=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(const Epetra_FECrsGraph &Graph)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::~Epetra_FECrsGraph "Epetra_FECrsGraph::~Epetra_FECrsGraph()

Destructor ";

%feature("docstring")  Epetra_FECrsGraph::InsertGlobalIndices "int
Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const int *rows,
int numCols, const int *cols)

Insert a rectangular, dense 'submatrix' of entries (matrix nonzero
positions) into the graph.

Parameters:
-----------

numRows:  Number of rows in the submatrix.

rows:  List of row-numbers for the submatrix.

numCols:  Number of columns in the submatrix.

cols:  List of column-indices that will be used for each row in the
'rows' list. ";

%feature("docstring")  Epetra_FECrsGraph::InsertGlobalIndices "int
Epetra_FECrsGraph::InsertGlobalIndices(int numRows, const long long
*rows, int numCols, const long long *cols) ";

%feature("docstring")  Epetra_FECrsGraph::GlobalAssemble "int
Epetra_FECrsGraph::GlobalAssemble(bool callFillComplete=true)

Gather any overlapping/shared data into the non-overlapping
partitioning defined by the Map that was passed to this matrix at
construction time. Data imported from other processors is stored on
the owning processor with a \"sumInto\" or accumulate operation. This
is a collective method -- every processor must enter it before any
will complete it.

NOTE***: When GlobalAssemble() calls FillComplete(), it passes the
arguments 'DomainMap()' and 'RangeMap()', which are the map attributes
held by the base-class CrsMatrix and its graph. If a rectangular
matrix is being assembled, the domain-map and range-map must be
specified by calling the other overloading of this method. Otherwise,
GlobalAssemble() has no way of knowing what these maps should really
be.

Parameters:
-----------

callFillComplete:  option argument, defaults to true. Determines
whether GlobalAssemble() internally calls the FillComplete() method on
this matrix.

error-code 0 if successful, non-zero if some error occurs ";

%feature("docstring")  Epetra_FECrsGraph::GlobalAssemble "int
Epetra_FECrsGraph::GlobalAssemble(const Epetra_Map &domain_map, const
Epetra_Map &range_map, bool callFillComplete=true)

Gather any overlapping/shared data into the non-overlapping
partitioning defined by the Map that was passed to this matrix at
construction time. Data imported from other processors is stored on
the owning processor with a \"sumInto\" or accumulate operation. This
is a collective method -- every processor must enter it before any
will complete it.

NOTE***: When GlobalAssemble() (the other overloading of this method)
calls FillComplete(), it passes the arguments 'DomainMap()' and
'RangeMap()', which are the map attributes already held by the base-
class CrsMatrix and its graph. If a rectangular matrix is being
assembled, the domain-map and range-map must be specified. Otherwise,
GlobalAssemble() has no way of knowing what these maps should really
be.

Parameters:
-----------

domain_map:  user-supplied domain map for this matrix

range_map:  user-supplied range map for this matrix

callFillComplete:  option argument, defaults to true. Determines
whether GlobalAssemble() internally calls the FillComplete() method on
this matrix.

error-code 0 if successful, non-zero if some error occurs ";

%feature("docstring")  Epetra_FECrsGraph::UseNonlocalGraph "bool
Epetra_FECrsGraph::UseNonlocalGraph() const ";


// File: classEpetra__FECrsMatrix.xml
%feature("docstring") Epetra_FECrsMatrix "

Epetra Finite-Element CrsMatrix. This class provides the ability to
input finite-element style sub-matrix data, including sub-matrices
with non-local rows (which could correspond to shared finite-element
nodes for example). This class inherits Epetra_CrsMatrix, and so all
Epetra_CrsMatrix functionality is also available.

It is intended that this class will be used as follows: Construct with
either a map or graph that describes a (non- overlapping) data
distribution.

Input data, including non-local data, using the methods
InsertGlobalValues(), SumIntoGlobalValues() and/or
ReplaceGlobalValues().

Call the method GlobalAssemble(), which gathers all non-local data
onto the owning processors as determined by the map provided at
construction. Users should note that the GlobalAssemble() method has
an optional argument which determines whether GlobalAssemble() in turn
calls FillComplete() after the data-exchange has occurred. If not
explicitly supplied, this argument defaults to true. NOTE***: When
GlobalAssemble() calls FillComplete(), it passes the arguments
'DomainMap()' and 'RangeMap()', which are the map attributes held by
the base-class CrsMatrix and its graph. If a rectangular matrix is
being assembled, the correct domain-map and range-map must be passed
to GlobalAssemble (there are two overloadings of this method) --
otherwise, it has no way of knowing what these maps should really be.

Sub-matrix data, which is assumed to be a rectangular 'table' of
coefficients accompanied by 'scatter-indices', can be provided in
three forms: Fortran-style packed 1-D array.

C-style double-pointer, or list-of-rows.

Epetra_SerialDenseMatrix object.  In all cases, a \"format\" parameter
specifies whether the data is laid out in row-major or column-major
order (i.e., whether coefficients for a row lie contiguously or
whether coefficients for a column lie contiguously). See the
documentation for the methods SumIntoGlobalValues() and
ReplaceGlobalValues().

Important notes: Since Epetra_FECrsMatrix inherits Epetra_CrsMatrix,
the semantics of the Insert/SumInto/Replace methods are the same as
they are on Epetra_CrsMatrix, which is:  InsertGlobalValues() inserts
values into the matrix only if the graph has not yet been finalized (
FillComplete() has not yet been called). For non-local values, the
call to InsertGlobalValues() may succeed but the GlobalAssemble()
method may then fail because the non-local data is not actually
inserted in the underlying matrix until GlobalAssemble() is called.

SumIntoGlobalValues() and ReplaceGlobalValues() only work for values
that already exist in the matrix. In other words, these methods can
not be used to put new values into the matrix.

C++ includes: Epetra_FECrsMatrix.h ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, int *NumEntriesPerRow, bool
ignoreNonLocalEntries=false)

Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, int NumEntriesPerRow, bool
ignoreNonLocalEntries=false)

Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, const Epetra_Map &ColMap, int *NumEntriesPerRow,
bool ignoreNonLocalEntries=false)

Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const
Epetra_Map &RowMap, const Epetra_Map &ColMap, int NumEntriesPerRow,
bool ignoreNonLocalEntries=false)

Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const
Epetra_CrsGraph &Graph, bool ignoreNonLocalEntries=false)

Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(Epetra_DataAccess CV, const
Epetra_FECrsGraph &Graph, bool ignoreNonLocalEntries=false)

Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(const Epetra_FECrsMatrix &src)

Copy Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::~Epetra_FECrsMatrix "Epetra_FECrsMatrix::~Epetra_FECrsMatrix()

Destructor. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int GlobalRow, int NumEntries,
const double *Values, const int *Indices)

override base-class Epetra_CrsMatrix::SumIntoGlobalValues method ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(long long GlobalRow, int
NumEntries, const double *Values, const long long *Indices) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
const double *Values, const int *Indices)

override base-class Epetra_CrsMatrix::InsertGlobalValues method ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(long long GlobalRow, int
NumEntries, const double *Values, const long long *Indices) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int GlobalRow, int NumEntries,
double *Values, int *Indices)

override base-class Epetra_CrsMatrix::InsertGlobalValues method ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(long long GlobalRow, int
NumEntries, double *Values, long long *Indices) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int GlobalRow, int NumEntries,
const double *Values, const int *Indices)

override base-class Epetra_CrsMatrix::ReplaceGlobalValues method ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(long long GlobalRow, int
NumEntries, const double *Values, const long long *Indices) ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int
*indices, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Sum a Fortran-style table (single-dimensional packed-list) of
coefficients into the matrix, adding them to any coefficients that may
already exist at the specified row/column locations.

Parameters:
-----------

numIndices:  Number of rows (and columns) in the sub-matrix.

indices:  List of scatter-indices (rows and columns) for the sub-
matrix.

values:  List, length numIndices*numIndices. Square sub-matrix of
coefficients, packed in a 1-D array. Data is packed either
contiguously by row or by column, specified by the final parameter
'format'.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const long
long *indices, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int *rows,
int numCols, const int *cols, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Sum a Fortran-style table (single-dimensional packed-list) of
coefficients into the matrix, adding them to any coefficients that may
already exist at the specified row/column locations.

Parameters:
-----------

numRows:  Number of rows in the sub-matrix.

rows:  List of row-numbers (scatter-indices) for the sub-matrix.

numCols:  Number of columns in the sub-matrix.

cols:  List of column-numbers (scatter-indices) for the sub-matrix.

values:  List, length numRows*numCols. Rectangular sub-matrix of
coefficients, packed in a 1-D array. Data is packed either
contiguously by row or by column, specified by the final parameter
'format'.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const long long
*rows, int numCols, const long long *cols, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const int
*indices, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR)

Sum C-style table (double-pointer, or list of lists) of coefficients
into the matrix, adding them to any coefficients that may already
exist at the specified row/column locations.

Parameters:
-----------

numIndices:  Number of rows (and columns) in the sub-matrix.

indices:  List of scatter-indices (rows and columns) for the sub-
matrix.

values:  Square sub-matrix of coefficients, provided in a 2-D array,
or double- pointer.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is ROW_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numIndices, const long
long *indices, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const int *rows,
int numCols, const int *cols, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR)

Sum C-style table (double-pointer, or list of lists) of coefficients
into the matrix, adding them to any coefficients that may already
exist at the specified row/column locations.

Parameters:
-----------

numRows:  Number of rows in the sub-matrix.

rows:  List of row-numbers (scatter-indices) for the sub-matrix.

numCols:  Number of columns in the sub-matrix.

cols:  List of column-numbers (scatter-indices) for the sub-matrix.

values:  Rectangular sub-matrix of coefficients, provided in a 2-D
array, or double-pointer.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is ROW_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(int numRows, const long long
*rows, int numCols, const long long *cols, const double *const
*values, int format=Epetra_FECrsMatrix::ROW_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int
*indices, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Insert a Fortran-style table (single-dimensional packed-list) of
coefficients into the matrix.

Parameters:
-----------

numIndices:  Number of rows (and columns) in the sub-matrix.

indices:  List of scatter-indices (rows and columns) for the sub-
matrix.

values:  List, length numIndices*numIndices. Square sub-matrix of
coefficients, packed in a 1-D array. Data is packed either
contiguously by row or by column, specified by the final parameter
'format'.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const long long
*indices, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int *rows,
int numCols, const int *cols, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Insert a Fortran-style table (single-dimensional packed-list) of
coefficients into the matrix.

Parameters:
-----------

numRows:  Number of rows in the sub-matrix.

rows:  List of row-numbers (scatter-indices) for the sub-matrix.

numCols:  Number of columns in the sub-matrix.

cols:  List of column-numbers (scatter-indices) for the sub-matrix.

values:  List, length numRows*numCols. Rectangular sub-matrix of
coefficients, packed in a 1-D array. Data is packed either
contiguously by row or by column, specified by the final parameter
'format'.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const long long
*rows, int numCols, const long long *cols, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const int
*indices, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR)

Insert a C-style table (double-pointer, or list of lists) of
coefficients into the matrix.

Parameters:
-----------

numIndices:  Number of rows (and columns) in the sub-matrix.

indices:  List of scatter-indices (rows and columns) for the sub-
matrix.

values:  Square sub-matrix of coefficients, provided in a 2-D array,
or double- pointer.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is ROW_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numIndices, const long long
*indices, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const int *rows,
int numCols, const int *cols, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR)

Insert a C-style table (double-pointer, or list of lists) of
coefficients into the matrix.

Parameters:
-----------

numRows:  Number of rows in the sub-matrix.

rows:  List of row-numbers (scatter-indices) for the sub-matrix.

numCols:  Number of columns in the sub-matrix.

cols:  List of column-numbers (scatter-indices) for the sub-matrix.

values:  Rectangular sub-matrix of coefficients, provided in a 2-D
array, or double-pointer.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is ROW_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(int numRows, const long long
*rows, int numCols, const long long *cols, const double *const
*values, int format=Epetra_FECrsMatrix::ROW_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int
*indices, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Copy a Fortran-style table (single-dimensional packed-list) of
coefficients into the matrix, replacing any coefficients that may
already exist at the specified row/column locations.

Parameters:
-----------

numIndices:  Number of rows (and columns) in the sub-matrix.

indices:  List of scatter-indices (rows and columns) for the sub-
matrix.

values:  List, length numIndices*numIndices. Square sub-matrix of
coefficients, packed in a 1-D array. Data is packed either
contiguously by row or by column, specified by the final parameter
'format'.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const long
long *indices, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int *rows,
int numCols, const int *cols, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Copy Fortran-style table (single-dimensional packed-list) of
coefficients into the matrix, replacing any coefficients that may
already exist at the specified row/column locations.

Parameters:
-----------

numRows:  Number of rows in the sub-matrix.

rows:  List of row-numbers (scatter-indices) for the sub-matrix.

numCols:  Number of columns in the sub-matrix.

cols:  List, of column-numbers (scatter-indices) for the sub-matrix.

values:  List, length numRows*numCols. Rectangular sub-matrix of
coefficients, packed in a 1-D array. Data is packed either
contiguously by row or by column, specified by the final parameter
'format'.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const long long
*rows, int numCols, const long long *cols, const double *values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const int
*indices, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR)

Copy C-style table (double-pointer, or list of lists) of coefficients
into the matrix, replacing any coefficients that may already exist at
the specified row/column locations.

Parameters:
-----------

numIndices:  Number of rows (and columns) in the sub-matrix.

indices:  List of scatter-indices (rows and columns) for the sub-
matrix.

values:  Square sub-matrix of coefficients, provided in a 2-D array,
or double- pointer.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is ROW_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numIndices, const long
long *indices, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const int *rows,
int numCols, const int *cols, const double *const *values, int
format=Epetra_FECrsMatrix::ROW_MAJOR)

Copy C-style table (double-pointer, or list of lists) of coefficients
into the matrix, replacing any coefficients that may already exist at
the specified row/column locations.

Parameters:
-----------

numRows:  Number of rows in the sub-matrix.

rows:  List of row-numbers (scatter-indices) for the sub-matrix.

numCols:  Number of columns in the sub-matrix.

cols:  List of column-numbers (scatter-indices) for the sub-matrix.

values:  Rectangular sub-matrix of coefficients, provided in a 2-D
array, or double-pointer.

format:  Specifies whether the data in 'values' is packed in column-
major or row-major order. Valid values are
Epetra_FECrsMatrix::ROW_MAJOR or Epetra_FECrsMatrix::COLUMN_MAJOR.
This is an optional parameter, default value is ROW_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(int numRows, const long long
*rows, int numCols, const long long *cols, const double *const
*values, int format=Epetra_FECrsMatrix::ROW_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(const
Epetra_IntSerialDenseVector &indices, const Epetra_SerialDenseMatrix
&values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Sum a square structurally-symmetric sub-matrix into the global matrix.
For non-square sub-matrices, see the other overloading of this method.

Parameters:
-----------

indices:  List of scatter-indices. indices.Length() must be the same
as values.M() and values.N().

values:  Sub-matrix of coefficients. Must be square.

format:  Optional format specifier, defaults to COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(const
Epetra_LongLongSerialDenseVector &indices, const
Epetra_SerialDenseMatrix &values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(const
Epetra_IntSerialDenseVector &rows, const Epetra_IntSerialDenseVector
&cols, const Epetra_SerialDenseMatrix &values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Sum a general sub-matrix into the global matrix. For square
structurally-symmetric sub-matrices, see the other overloading of this
method.

Parameters:
-----------

rows:  List of row-indices. rows.Length() must be the same as
values.M().

cols:  List of column-indices. cols.Length() must be the same as
values.N().

values:  Sub-matrix of coefficients.

format:  Optional format specifier, defaults to COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::SumIntoGlobalValues "int
Epetra_FECrsMatrix::SumIntoGlobalValues(const
Epetra_LongLongSerialDenseVector &rows, const
Epetra_LongLongSerialDenseVector &cols, const Epetra_SerialDenseMatrix
&values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(const
Epetra_IntSerialDenseVector &indices, const Epetra_SerialDenseMatrix
&values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Insert a square structurally-symmetric sub-matrix into the global
matrix. For non-square sub-matrices, see the other overloading of this
method.

Parameters:
-----------

indices:  List of scatter-indices. indices.Length() must be the same
as values.M() and values.N().

values:  Sub-matrix of coefficients. Must be square.

format:  Optional format specifier, defaults to COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(const
Epetra_LongLongSerialDenseVector &indices, const
Epetra_SerialDenseMatrix &values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(const
Epetra_IntSerialDenseVector &rows, const Epetra_IntSerialDenseVector
&cols, const Epetra_SerialDenseMatrix &values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Insert a general sub-matrix into the global matrix. For square
structurally-symmetric sub-matrices, see the other overloading of this
method.

Parameters:
-----------

rows:  List of row-indices. rows.Length() must be the same as
values.M().

cols:  List of column-indices. cols.Length() must be the same as
values.N().

values:  Sub-matrix of coefficients.

format:  Optional format specifier, defaults to COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::InsertGlobalValues "int
Epetra_FECrsMatrix::InsertGlobalValues(const
Epetra_LongLongSerialDenseVector &rows, const
Epetra_LongLongSerialDenseVector &cols, const Epetra_SerialDenseMatrix
&values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(const
Epetra_IntSerialDenseVector &indices, const Epetra_SerialDenseMatrix
&values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Use a square structurally-symmetric sub-matrix to replace existing
values in the global matrix. For non-square sub-matrices, see the
other overloading of this method.

Parameters:
-----------

indices:  List of scatter-indices. indices.Length() must be the same
as values.M() and values.N().

values:  Sub-matrix of coefficients. Must be square.

format:  Optional format specifier, defaults to COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(const
Epetra_LongLongSerialDenseVector &indices, const
Epetra_SerialDenseMatrix &values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(const
Epetra_IntSerialDenseVector &rows, const Epetra_IntSerialDenseVector
&cols, const Epetra_SerialDenseMatrix &values, int
format=Epetra_FECrsMatrix::COLUMN_MAJOR)

Use a general sub-matrix to replace existing values. For square
structurally-symmetric sub-matrices, see the other overloading of this
method.

Parameters:
-----------

rows:  List of row-indices. rows.Length() must be the same as
values.M().

cols:  List of column-indices. cols.Length() must be the same as
values.N().

values:  Sub-matrix of coefficients.

format:  Optional format specifier, defaults to COLUMN_MAJOR. ";

%feature("docstring")  Epetra_FECrsMatrix::ReplaceGlobalValues "int
Epetra_FECrsMatrix::ReplaceGlobalValues(const
Epetra_LongLongSerialDenseVector &rows, const
Epetra_LongLongSerialDenseVector &cols, const Epetra_SerialDenseMatrix
&values, int format=Epetra_FECrsMatrix::COLUMN_MAJOR) ";

%feature("docstring")  Epetra_FECrsMatrix::GlobalAssemble "int
Epetra_FECrsMatrix::GlobalAssemble(bool callFillComplete=true,
Epetra_CombineMode combineMode=Add, bool
save_off_and_reuse_map_exporter=false)

Gather any overlapping/shared data into the non-overlapping
partitioning defined by the Map that was passed to this matrix at
construction time. Data imported from other processors is stored on
the owning processor with a \"sumInto\" or accumulate operation. This
is a collective method -- every processor must enter it before any
will complete it.

NOTE***: When GlobalAssemble() calls FillComplete(), it passes the
arguments 'DomainMap()' and 'RangeMap()', which are the map attributes
held by the base-class CrsMatrix and its graph. If a rectangular
matrix is being assembled, the domain-map and range-map must be
specified by calling the other overloading of this method. Otherwise,
GlobalAssemble() has no way of knowing what these maps should really
be.

Parameters:
-----------

callFillComplete:  option argument, defaults to true. Determines
whether GlobalAssemble() internally calls the FillComplete() method on
this matrix.

error-code 0 if successful, non-zero if some error occurs ";

%feature("docstring")  Epetra_FECrsMatrix::GlobalAssemble "int
Epetra_FECrsMatrix::GlobalAssemble(const Epetra_Map &domain_map, const
Epetra_Map &range_map, bool callFillComplete=true, Epetra_CombineMode
combineMode=Add, bool save_off_and_reuse_map_exporter=false)

Gather any overlapping/shared data into the non-overlapping
partitioning defined by the Map that was passed to this matrix at
construction time. Data imported from other processors is stored on
the owning processor with a \"sumInto\" or accumulate operation. This
is a collective method -- every processor must enter it before any
will complete it.

NOTE***: When GlobalAssemble() (the other overloading of this method)
calls FillComplete(), it passes the arguments 'DomainMap()' and
'RangeMap()', which are the map attributes already held by the base-
class CrsMatrix and its graph. If a rectangular matrix is being
assembled, the domain-map and range-map must be specified. Otherwise,
GlobalAssemble() has no way of knowing what these maps should really
be.

Parameters:
-----------

domain_map:  user-supplied domain map for this matrix

range_map:  user-supplied range map for this matrix

callFillComplete:  option argument, defaults to true. Determines
whether GlobalAssemble() internally calls the FillComplete() method on
this matrix.

error-code 0 if successful, non-zero if some error occurs ";

%feature("docstring")  Epetra_FECrsMatrix::setIgnoreNonLocalEntries "void Epetra_FECrsMatrix::setIgnoreNonLocalEntries(bool flag)

Set whether or not non-local data values should be ignored. By
default, non-local data values are NOT ignored. ";


// File: classEpetra__FEVbrMatrix.xml
%feature("docstring") Epetra_FEVbrMatrix "

Epetra Finite-Element VbrMatrix. This class provides the ability to
input finite-element style sub-matrix data, including sub-matrices
with non-local rows (which could correspond to shared finite-element
nodes for example). This class inherits Epetra_VbrMatrix, and so all
Epetra_VbrMatrix functionality is also available.

C++ includes: Epetra_FEVbrMatrix.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_FEVbrMatrix::Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int *NumBlockEntriesPerRow, bool
ignoreNonLocalEntries=false)

Epetra_FEVbrMatrix constuctor with variable number of indices per row.

Creates a Epetra_FEVbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - A Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  NumBlockEntriesPerRow - An integer array of length NumRows such
that NumBlockEntriesPerRow[i] indicates the (approximate) number of
Block entries in the ith row. ";

%feature("docstring")  Epetra_FEVbrMatrix::Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int NumBlockEntriesPerRow, bool
ignoreNonLocalEntries=false)

Epetra_FEVbrMatrix constuctor with fixed number of indices per row.

Creates a Epetra_FEVbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - An Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  NumBlockEntriesPerRow - An integer that indicates the
(approximate) number of Block entries in the each Block row. Note that
it is possible to use 0 for this value and let fill occur during the
insertion phase. ";

%feature("docstring")  Epetra_FEVbrMatrix::Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
*NumBlockEntriesPerRow, bool ignoreNonLocalEntries=false)

Epetra_FEVbrMatrix constuctor with variable number of indices per row.

Creates a Epetra_FEVbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - A Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  ColMap - A Epetra_BlockMap listing the block columns to be
contained on this processor.

In:  NumBlockEntriesPerRow - An integer array of length NumRows such
that NumBlockEntriesPerRow[i] indicates the (approximate) number of
Block entries in the ith row. ";

%feature("docstring")  Epetra_FEVbrMatrix::Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
NumBlockEntriesPerRow, bool ignoreNonLocalEntries=false)

Epetra_FEVbrMatrix constuctor with fixed number of indices per row.

Creates a Epetra_FEVbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - An Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  ColMap - An Epetra_BlockMap listing the block columns to be
contained on this processor.

In:  NumBlockEntriesPerRow - An integer that indicates the
(approximate) number of Block entries in the each Block row. Note that
it is possible to use 0 for this value and let fill occur during the
insertion phase. ";

%feature("docstring")  Epetra_FEVbrMatrix::Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(Epetra_DataAccess CV, const
Epetra_CrsGraph &Graph, bool ignoreNonLocalEntries=false)

Constructor with pre-constructed Graph. ";

%feature("docstring")  Epetra_FEVbrMatrix::Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::Epetra_FEVbrMatrix(const Epetra_FEVbrMatrix &src)

Copy Constructor. ";

%feature("docstring")  Epetra_FEVbrMatrix::~Epetra_FEVbrMatrix "Epetra_FEVbrMatrix::~Epetra_FEVbrMatrix()

Epetra_VbrMatrix Destructor. ";

/*  Insertion/Replace/SumInto methods  */

%feature("docstring")  Epetra_FEVbrMatrix::PutScalar "int
Epetra_FEVbrMatrix::PutScalar(double ScalarConstant)

Initialize all values in graph of the matrix with constant value.

Parameters:
-----------

In:  ScalarConstant - Value to use.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_FEVbrMatrix::BeginInsertGlobalValues "int Epetra_FEVbrMatrix::BeginInsertGlobalValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate insertion of a list of elements in a given global row of the
matrix, values are inserted via SubmitEntry().

Parameters:
-----------

In:  BlockRow - Block Row number (in global coordinates) to put
elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_FEVbrMatrix::BeginReplaceGlobalValues "int Epetra_FEVbrMatrix::BeginReplaceGlobalValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate replacement of current values with this list of entries for a
given global row of the matrix, values are replaced via SubmitEntry()

Parameters:
-----------

In:  Row - Block Row number (in global coordinates) to put elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_FEVbrMatrix::BeginSumIntoGlobalValues "int Epetra_FEVbrMatrix::BeginSumIntoGlobalValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate summing into current values with this list of entries for a
given global row of the matrix, values are replaced via SubmitEntry()

Parameters:
-----------

In:  Row - Block Row number (in global coordinates) to put elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_FEVbrMatrix::SubmitBlockEntry "int
Epetra_FEVbrMatrix::SubmitBlockEntry(double *Values, int LDA, int
NumRows, int NumCols)

Submit a block entry to the indicated block row and column specified
in the Begin routine. ";

%feature("docstring")  Epetra_FEVbrMatrix::EndSubmitEntries "int
Epetra_FEVbrMatrix::EndSubmitEntries()

Completes processing of all data passed in for the current block row.

This function completes the processing of all block entries submitted
via SubmitBlockEntry(). It also checks to make sure that
SubmitBlockEntry was called the correct number of times as specified
by the Begin routine that initiated the entry process. ";

%feature("docstring")  Epetra_FEVbrMatrix::GlobalAssemble "int
Epetra_FEVbrMatrix::GlobalAssemble(bool callFillComplete=true) ";


// File: classEpetra__FEVector.xml
%feature("docstring") Epetra_FEVector "

Epetra Finite-Element Vector. This class inherits Epetra_MultiVector
and thus provides all Epetra_MultiVector functionality.

The added functionality provided by Epetra_FEVector is the ability to
perform finite-element style vector assembly. It accepts sub-vector
contributions, such as those that would come from element-load
vectors, etc. These sub-vectors need not be owned by the local
processor. In other words, the user can assemble overlapping data
(e.g., corresponding to shared finite-element nodes). When the user is
finished assembling their vector data, they then call the method
Epetra_FEVector::GlobalAssemble() which gathers the overlapping data
(all non-local data that was input on each processor) into the data-
distribution specified by the map with which the Epetra_FEVector was
constructed.

C++ includes: Epetra_FEVector.h ";

%feature("docstring")  Epetra_FEVector::Epetra_FEVector "Epetra_FEVector::Epetra_FEVector(const Epetra_BlockMap &Map, int
numVectors=1, bool ignoreNonLocalEntries=false)

Constructor that requires a map specifying a non-overlapping data
layout.

Parameters:
-----------

Map:  Map describing a non-overlapping distribution for the underlying
Epetra_MultiVector into which this Epetra_FEVector will funnel data.

numVectors:  Optional argument, default value is 1. (See the
documentation for Epetra_MultiVector for the meaning of this argument.

ignoreNonLocalEntries:  Optional argument, default value is false.
Under certain special circumstances it is desirable to have non-local
contributions ignored rather than saving them for the GlobalAssemble
step. ";

%feature("docstring")  Epetra_FEVector::Epetra_FEVector "Epetra_FEVector::Epetra_FEVector(Epetra_DataAccess CV, const
Epetra_BlockMap &Map, double *A, int MyLDA, int NumVectors, bool
ignoreNonLocalEntries=false)

Set multi-vector values from two-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  A - Pointer to an array of double precision numbers. The first
vector starts at A. The second vector starts at A+MyLDA, the third at
A+2*MyLDA, and so on.

In:  MyLDA - The \"Leading Dimension\", or stride between vectors in
memory.

WARNING:  This value refers to the stride on the calling processor.
Thus it is a local quantity, not a global quantity.

Parameters:
-----------

In:  NumVectors - Number of vectors in multi-vector.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_FEVector::Epetra_FEVector "Epetra_FEVector::Epetra_FEVector(Epetra_DataAccess CV, const
Epetra_BlockMap &Map, double **ArrayOfPointers, int NumVectors, bool
ignoreNonLocalEntries=false)

Set multi-vector values from array of pointers.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  ArrayOfPointers - An array of pointers such that
ArrayOfPointers[i] points to the memory location containing ith vector
to be copied.

In:  NumVectors - Number of vectors in multi-vector.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_FEVector::Epetra_FEVector "Epetra_FEVector::Epetra_FEVector(const Epetra_FEVector &source)

Copy constructor. ";

%feature("docstring")  Epetra_FEVector::~Epetra_FEVector "Epetra_FEVector::~Epetra_FEVector()

Destructor ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int *GIDs,
const double *values, int vectorIndex=0)

Accumulate values into the vector, adding them to any values that
already exist for the specified indices. ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(int numIDs, const long long
*GIDs, const double *values, int vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(const Epetra_IntSerialDenseVector
&GIDs, const Epetra_SerialDenseVector &values, int vectorIndex=0)

Accumulate values into the vector, adding them to any values that
already exist for the specified GIDs.

Parameters:
-----------

GIDs:  List of global ids. Must be the same length as the accompanying
list of values.

values:  List of coefficient values. Must be the same length as the
accompanying list of GIDs. ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(const
Epetra_LongLongSerialDenseVector &GIDs, const Epetra_SerialDenseVector
&values, int vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int *GIDs,
const double *values, int vectorIndex=0)

Copy values into the vector overwriting any values that already exist
for the specified indices. ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(int numIDs, const long long
*GIDs, const double *values, int vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(const Epetra_IntSerialDenseVector
&GIDs, const Epetra_SerialDenseVector &values, int vectorIndex=0)

Copy values into the vector, replacing any values that already exist
for the specified GIDs.

Parameters:
-----------

GIDs:  List of global ids. Must be the same length as the accompanying
list of values.

values:  List of coefficient values. Must be the same length as the
accompanying list of GIDs. ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(const
Epetra_LongLongSerialDenseVector &GIDs, const Epetra_SerialDenseVector
&values, int vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int *GIDs,
const int *numValuesPerID, const double *values, int vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(int numIDs, const long long
*GIDs, const int *numValuesPerID, const double *values, int
vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int *GIDs,
const int *numValuesPerID, const double *values, int vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(int numIDs, const long long
*GIDs, const int *numValuesPerID, const double *values, int
vectorIndex=0) ";

%feature("docstring")  Epetra_FEVector::GlobalAssemble "int
Epetra_FEVector::GlobalAssemble(Epetra_CombineMode mode=Add, bool
reuse_map_and_exporter=false)

Gather any overlapping/shared data into the non-overlapping
partitioning defined by the Map that was passed to this vector at
construction time. Data imported from other processors is stored on
the owning processor with a \"sumInto\" or accumulate operation. This
is a collective method -- every processor must enter it before any
will complete it.

Optimization for power-users: The optional parameter
'reuse_map_and_exporter' defaults to false. By default, a map that
describes the non-local data is re-created at each call to
GlobalAssemble, along with an exporter used to do the communication.
This is expensive. If you know that the layout of your nonlocal data
has not changed since your previous call to GlobalAssemble, you can
set this flag to true and it will reuse the previously created map and
exporter rather than creating new ones. ";

%feature("docstring")  Epetra_FEVector::setIgnoreNonLocalEntries "void Epetra_FEVector::setIgnoreNonLocalEntries(bool flag)

Set whether or not non-local data values should be ignored. ";


// File: classEpetra__Flops.xml
%feature("docstring") Epetra_Flops "

Epetra_Flops: The Epetra Floating Point Operations Class.

The Epetra_Flops class provides basic support and consistent
interfaces for counting and reporting floating point operations
performed in the Epetra computational classes. All classes based on
the Epetra_CompObject can count flops by the user creating an
Epetra_Flops object and calling the SetFlopCounter() method for an
Epetra_CompObject.

C++ includes: Epetra_Flops.h ";

%feature("docstring")  Epetra_Flops::Epetra_Flops "Epetra_Flops::Epetra_Flops(void)

Epetra_Flops Constructor.

Creates a Epetra_Flops instance. This instance can be queried for the
number of floating point operations performed for the associated this
object. ";

%feature("docstring")  Epetra_Flops::Epetra_Flops "Epetra_Flops::Epetra_Flops(const Epetra_Flops &Flops_in)

Epetra_Flops Copy Constructor.

Makes an exact copy of an existing Epetra_Flops instance. ";

%feature("docstring")  Epetra_Flops::Flops "double
Epetra_Flops::Flops() const

Returns the number of floating point operations with this object and
resets the count. ";

%feature("docstring")  Epetra_Flops::ResetFlops "void
Epetra_Flops::ResetFlops()

Resets the number of floating point operations to zero for this multi-
vector. ";

%feature("docstring")  Epetra_Flops::~Epetra_Flops "Epetra_Flops::~Epetra_Flops(void)

Epetra_Flops Destructor.

Completely deletes a Epetra_Flops object. ";


// File: classEpetra__HashTable.xml
%feature("docstring") Epetra_HashTable "";

%feature("docstring")  Epetra_HashTable::Epetra_HashTable "Epetra_HashTable< value_type >::Epetra_HashTable(const int size, const
unsigned int seed=(2654435761U)) ";

%feature("docstring")  Epetra_HashTable::Epetra_HashTable "Epetra_HashTable< value_type >::Epetra_HashTable(const
Epetra_HashTable &obj) ";

%feature("docstring")  Epetra_HashTable::~Epetra_HashTable "Epetra_HashTable< value_type >::~Epetra_HashTable() ";

%feature("docstring")  Epetra_HashTable::Add "void Epetra_HashTable<
value_type >::Add(const long long key, const value_type value) ";

%feature("docstring")  Epetra_HashTable::Get "value_type
Epetra_HashTable< value_type >::Get(const long long key) ";


// File: classEpetra__Import.xml
%feature("docstring") Epetra_Import "

Epetra_Import: This class builds an import object for efficient
importing of off- processor elements.

Epetra_Import is used to construct a communication plan that can be
called repeatedly by computational classes such the Epetra matrix,
vector and multivector classes to efficiently obtain off-processor
elements.

This class currently has one constructor, taking two Epetra_Map or
Epetra_BlockMap objects. The first map specifies the global IDs of
elements that we want to import later. The second map specifies the
global IDs that are owned by the calling processor.

C++ includes: Epetra_Import.h ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_Import::Print "void
Epetra_Import::Print(ostream &os) const

Print object to an output stream Print method ";

%feature("docstring")  Epetra_Import::Epetra_Import "Epetra_Import::Epetra_Import(const Epetra_BlockMap &TargetMap, const
Epetra_BlockMap &SourceMap)

Constructs a Epetra_Import object from the source and target maps.

This constructor builds an Epetra_Import object by comparing the GID
lists of the source and target maps.

Parameters:
-----------

TargetMap:  (In) Map containing the GIDs from which data should be
imported to each processor from the source map whenever an import
operation is performed using this importer.

SourceMap:  (In) Map containing the GIDs that should be used for
importing data.

WARNING:  Note that the SourceMap must have GIDs uniquely owned, each
GID of the source map can occur only once.  Builds an import object
that will transfer objects built with SourceMap to objects built with
TargetMap.

A Epetra_Import object categorizes the elements of the target map into
three sets as follows: All elements in the target map that have the
same GID as the corresponding element of the source map, starting with
the first element in the target map, going up to the first element
that is different from the source map. The number of these IDs is
returned by NumSameIDs().

All elements that are local to the processor, but are not part of the
first set of elements. These elements have GIDs that are owned by the
calling processor, but at least the first element of this list is
permuted. Even if subsequent elements are not permuted, they are
included in this list. The number of permuted elements is returned by
NumPermutedIDs(). The list of elements (local IDs) in the source map
that are permuted can be found in the list PermuteFromLIDs(). The list
of elements (local IDs) in the target map that are the new locations
of the source elements can be found in the list PermuteToLIDs().

All remaining elements of the target map correspond to global IDs that
are owned by remote processors. The number of these elements is
returned by NumRemoteIDs() and the list of these is returned by
RemoteLIDs().

Given the above information, the Epetra_Import constructor builds a
list of elements that must be communicated to other processors as a
result of import requests. The number of exported elements (where
multiple sends of the same element to different processors is counted)
is returned by NumExportIDs(). The local IDs to be sent are returned
by the list ExportLIDs(). The processors to which each of the elements
will be sent in returned in a list of the same length by ExportPIDs().

The total number of elements that will be sent by the calling
processor is returned by NumSend(). The total number of elements that
will be received is returned by NumRecv().

The following example illustrates the basic concepts.

Assume we have 3 processors and 9 global elements with each processor
owning 3 elements as follows  PE 0 Elements |  PE 1 Elements  |  PE 2
Elements     0  1  2 3  4  5           6  7  8

The above layout essentially defines the source map argument of the
import object.

This could correspond to a 9 by 9 matrix with the first three rows on
PE 0, and so on. Suppose that this matrix is periodic tridiagonal
having the following sparsity pattern:

PE 0 Rows:    X  X  0  0  0  0  0  0  X   X  X  X  0  0  0  0  0  0 0
X  X  X  0  0  0  0  0  PE 1 Rows:    0  0  X  X  X  0  0  0  0   0 0
0  X  X  X  0  0  0   0  0  0  0  X  X  X  0  0  PE 2 Rows:    0  0 0
0  0  X  X  X  0   0  0  0  0  0  0  X  X  X   X  0  0  0  0  0  0 X
X

To perform a matrix vector multiplication operation y = A*x (assuming
that x has the same distribution as the rows of the matrix A) each
processor will need to import elements of x that are not local. To do
this, we build a target map on each processor as follows:     PE 0
Elements    |  PE 1 Elements    |  PE 2 Elements     0  1  2 3  8
2  3  4  5  6       0  5  6  7  8

The above list is the elements that will be needed to perform the
matrix vector multiplication locally on each processor. Note that the
ordering of the elements on each processor is not unique, but has been
chosen for illustration.

With these two maps passed into the Epetra_Import constructor, we get
the following attribute definitions:

On PE 0:

NumSameIDs      = 3  NumPermuteIDs   = 0 PermuteToLIDs   = 0
PermuteFromLIDs = 0  NumRemoteIDs    = 2 RemoteLIDs      = [3, 4]
NumExportIDs    = 2 ExportLIDs      = [0, 2] ExportPIDs      = [1, 2]
NumSend         = 2 NumRecv         = 2

On PE 1:

NumSameIDs      = 0  NumPermuteIDs   = 3 PermuteFromLIDs = [0, 1, 2]
PermuteToLIDs   = [1, 2, 3]  NumRemoteIDs    = 2 RemoteLIDs      = [0,
4]  NumExportIDs    = 2 ExportLIDs      = [0, 2] ExportPIDs      = [0,
2]  NumSend         = 2 NumRecv         = 2

On PE 2:

NumSameIDs      = 0  NumPermuteIDs   = 3 PermuteFromLIDs = [0, 1, 2]
PermuteToLIDs   = [2, 3, 4]  NumRemoteIDs    = 2 RemoteLIDs      = [0,
1]  NumExportIDs    = 2 ExportLIDs      = [0, 2] ExportPIDs      = [0,
1]  NumSend         = 2 NumRecv         = 2

Using Epetra_Import Objects

Once a Epetra_Import object has been constructed, it can be used by
any of the Epetra classes that support distributed global objects,
namely Epetra_Vector, Epetra_MultiVector, Epetra_CrsGraph,
Epetra_CrsMatrix and Epetra_VbrMatrix. All of these classes have
Import and Export methods that will fill new objects whose
distribution is described by the target map, taking elements from the
source object whose distribution is described by the source map.
Details of usage for each class is given in the appropriate class
documentation.

Note that the reverse operation, an export, using this importer is
also possible and appropriate in some instances. For example, if we
compute y = A^Tx, the transpose matrix-multiplication operation, then
we can use the importer we constructed in the above example to do an
export operation to y, adding the contributions that come from
multiple processors. ";

%feature("docstring")  Epetra_Import::Epetra_Import "Epetra_Import::Epetra_Import(const Epetra_Import &Importer)

Epetra_Import copy constructor. ";

%feature("docstring")  Epetra_Import::~Epetra_Import "Epetra_Import::~Epetra_Import(void)

Epetra_Import destructor. ";

%feature("docstring")  Epetra_Import::NumSameIDs "int
Epetra_Import::NumSameIDs() const

Returns the number of elements that are identical between the source
and target maps, up to the first different ID. ";

%feature("docstring")  Epetra_Import::NumPermuteIDs "int
Epetra_Import::NumPermuteIDs() const

Returns the number of elements that are local to the calling
processor, but not part of the first NumSameIDs() elements. ";

%feature("docstring")  Epetra_Import::PermuteFromLIDs "int*
Epetra_Import::PermuteFromLIDs() const

List of elements in the source map that are permuted. ";

%feature("docstring")  Epetra_Import::PermuteToLIDs "int*
Epetra_Import::PermuteToLIDs() const

List of elements in the target map that are permuted. ";

%feature("docstring")  Epetra_Import::NumRemoteIDs "int
Epetra_Import::NumRemoteIDs() const

Returns the number of elements that are not on the calling processor.
";

%feature("docstring")  Epetra_Import::RemoteLIDs "int*
Epetra_Import::RemoteLIDs() const

List of elements in the target map that are coming from other
processors. ";

%feature("docstring")  Epetra_Import::NumExportIDs "int
Epetra_Import::NumExportIDs() const

Returns the number of elements that must be sent by the calling
processor to other processors. ";

%feature("docstring")  Epetra_Import::ExportLIDs "int*
Epetra_Import::ExportLIDs() const

List of elements that will be sent to other processors. ";

%feature("docstring")  Epetra_Import::ExportPIDs "int*
Epetra_Import::ExportPIDs() const

List of processors to which elements will be sent, ExportLIDs() [i]
will be sent to processor ExportPIDs() [i]. ";

%feature("docstring")  Epetra_Import::NumSend "int
Epetra_Import::NumSend() const

Total number of elements to be sent. ";

%feature("docstring")  Epetra_Import::NumRecv "int
Epetra_Import::NumRecv() const

Total number of elements to be received. ";

%feature("docstring")  Epetra_Import::SourceMap "const
Epetra_BlockMap& Epetra_Import::SourceMap() const

Returns the SourceMap used to construct this importer. ";

%feature("docstring")  Epetra_Import::TargetMap "const
Epetra_BlockMap& Epetra_Import::TargetMap() const

Returns the TargetMap used to construct this importer. ";

%feature("docstring")  Epetra_Import::Distributor "Epetra_Distributor& Epetra_Import::Distributor() const ";


// File: classEpetra__IntSerialDenseMatrix.xml
%feature("docstring") Epetra_IntSerialDenseMatrix "

Epetra_IntSerialDenseMatrix: A class for constructing and using
general dense integer matrices.

The Epetra_IntSerialDenseMatrix class enables the construction and use
of integer-valued, general dense matrices.

The Epetra_IntSerialDenseMatrix class is intended to provide very
basic support for dense rectangular matrices.

Constructing Epetra_IntSerialDenseMatrix Objects

There are four Epetra_IntSerialDenseMatrix constructors. The first
constructs a zero-sized object which should be made to appropriate
length using the Shape() or Reshape() functions and then filled with
the [] or () operators. The second constructs an object sized to the
dimensions specified, which should be filled with the [] or ()
operators. The third is a constructor that accepts user data as a 2D
array, and the fourth is a copy constructor. The third constructor has
two data access modes (specified by the Epetra_DataAccess argument):
Copy mode - Allocates memory and makes a copy of the user-provided
data. In this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Epetra_IntSerialDenseMatrix constructors will
throw an exception if an error occurrs. These exceptions will alway be
negative integer values as follows: -1 Invalid dimension specified.

-2 Shape returned non-zero.

-3 Null pointer specified for user's data.

-99 Internal Epetra_IntSerialDenseMatrix error. Contact developer.

Other Epetra_IntSerialDenseMatrix functions that do not return an
integer error code (such as operators () and [] ) will throw an
exception if an error occurrs. These exceptions will be integer values
as follows: -1 Invalid row specified.

-2 Invalid column specified.

-5 Invalid assignment (type mismatch).

-99 Internal Epetra_IntSerialDenseMatrix error. Contact developer.

b Extracting Data from Epetra_IntSerialDenseMatrix Objects

Once a Epetra_IntSerialDenseMatrix is constructed, it is possible to
view the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.  Vector and Utility Functions

Once a Epetra_IntSerialDenseMatrix is constructed, several
mathematical functions can be applied to the object. Specifically:
Multiplication.

Norms.

C++ includes: Epetra_IntSerialDenseMatrix.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix "Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix()

Default constructor; defines a zero size object.

Epetra_IntSerialDenseMatrix objects defined by the default constructor
should be sized with the Shape() or Reshape functions. Values should
be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix "Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(int NumRows,
int NumCols)

Shaped constructor; defines a variable-sized object.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Epetra_SerialDenseMatrix objects defined by the shaped constructor are
already shaped to the dimensions given as a parameters. All values are
initialized to 0. Calling this constructor is equivalent to using the
default constructor, and then calling the Shape function on it. Values
should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix "Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(Epetra_DataAccess
CV, int *A, int LDA, int NumRows, int NumCols)

Set object values from two-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  A - Pointer to an array of integer numbers. The first vector
starts at A. The second vector starts at A+LDA, the third at A+2*LDA,
and so on.

In:  LDA - The \"Leading Dimension\", or stride between vectors in
memory.

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

See Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix "Epetra_IntSerialDenseMatrix::Epetra_IntSerialDenseMatrix(const
Epetra_IntSerialDenseMatrix &Source)

Epetra_IntSerialDenseMatrix copy constructor.

This matrix will take on the data access mode of the Source matrix. ";

%feature("docstring")
Epetra_IntSerialDenseMatrix::~Epetra_IntSerialDenseMatrix "Epetra_IntSerialDenseMatrix::~Epetra_IntSerialDenseMatrix()

Epetra_IntSerialDenseMatrix destructor. ";

/*  Shaping/sizing Methods  */

%feature("docstring")  Epetra_IntSerialDenseMatrix::Shape "int
Epetra_IntSerialDenseMatrix::Shape(int NumRows, int NumCols)

Set dimensions of a Epetra_IntSerialDenseMatrix object; init values to
zero.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Allows user to define the dimensions of a Epetra_IntSerialDenseMatrix
at any point. This function can be called at any point after
construction. Any values that were previously in this object are
destroyed and the resized matrix starts off with all zero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::Reshape "int
Epetra_IntSerialDenseMatrix::Reshape(int NumRows, int NumCols)

Reshape a Epetra_IntSerialDenseMatrix object.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Allows user to define the dimensions of a Epetra_IntSerialDenseMatrix
at any point. This function can be called at any point after
construction. Any values that were previously in this object are
copied into the new shape. If the new shape is smaller than the
original, the upper left portion of the original matrix (the principal
submatrix) is copied to the new matrix.

Integer error code, set to 0 if successful. ";

/*  Data Accessor methods  */

%feature("docstring")  Epetra_IntSerialDenseMatrix::OneNorm "int
Epetra_IntSerialDenseMatrix::OneNorm()

Computes the 1-Norm of the this matrix.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::InfNorm "int
Epetra_IntSerialDenseMatrix::InfNorm()

Computes the Infinity-Norm of the this matrix. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::Random "int
Epetra_IntSerialDenseMatrix::Random()

Set matrix values to random numbers.

IntSerialDenseMatrix uses the random number generator provided by
Epetra_Util. The matrix values will be set to random values on the
interval (0, 2^31 - 1).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::M "int
Epetra_IntSerialDenseMatrix::M() const

Returns row dimension of system. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::N "int
Epetra_IntSerialDenseMatrix::N() const

Returns column dimension of system. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::A "const int*
Epetra_IntSerialDenseMatrix::A() const

Returns const pointer to the this matrix. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::A "int*
Epetra_IntSerialDenseMatrix::A()

Returns pointer to the this matrix. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::LDA "int
Epetra_IntSerialDenseMatrix::LDA() const

Returns the leading dimension of the this matrix. ";

%feature("docstring")  Epetra_IntSerialDenseMatrix::CV "Epetra_DataAccess Epetra_IntSerialDenseMatrix::CV() const

Returns the data access mode of the this matrix. ";

/*  I/O methods  */

%feature("docstring")  Epetra_IntSerialDenseMatrix::Print "void
Epetra_IntSerialDenseMatrix::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";

/*  Expert-only unsupported methods  */

%feature("docstring")  Epetra_IntSerialDenseMatrix::MakeViewOf "int
Epetra_IntSerialDenseMatrix::MakeViewOf(const
Epetra_IntSerialDenseMatrix &Source)

Reset an existing IntSerialDenseMatrix to point to another Matrix.

Allows an existing IntSerialDenseMatrix to become a View of another
matrix's data, regardless of the DataAccess mode of the Source matrix.
It is assumed that the Source matrix is an independent matrix, and no
checking is done to verify this.

This is used by Epetra_CrsGraph in the OptimizeStorage method. It is
used so that an existing (Copy) matrix can be converted to a View.
This frees up memory that CrsGraph no longer needs.

Parameters:
-----------

Source:  The IntSerialDenseMatrix this will become a view of.

Integer error code, set to 0 if successful, and set to -1 if a type
mismatch occured.

WARNING:  This method is extremely dangerous and should only be used
by experts. ";


// File: classEpetra__IntSerialDenseVector.xml
%feature("docstring") Epetra_IntSerialDenseVector "

Epetra_IntSerialDenseVector: A class for constructing and using dense
vectors.

The Epetra_IntSerialDenseVector class enables the construction and use
of integer-valued, dense vectors. It derives from the
Epetra_IntSerialDenseMatrix class.

The Epetra_IntSerialDenseVector class is intended to provide
convenient vector notation but derives all signficant functionality
from Epetra_IntSerialDenseMatrix.

Constructing Epetra_IntSerialDenseVector Objects

There are three Epetra_IntSerialDenseVector constructors. The first
constructs a zero-length object which should be made to appropriate
length using the Size() or Resize() functions and then filled with the
[] or () operators. The second constructs an object sized to the
dimension specified, which should be filled with the [] or ()
operators. The third is a constructor that accepts user data as a 1D
array, and the fourth is a copy constructor. The third constructor has
two data access modes (specified by the Epetra_DataAccess argument):
Copy mode - Allocates memory and makes a copy of the user-provided
data. In this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Extracting Data from Epetra_IntSerialDenseVector
Objects

Once a Epetra_IntSerialDenseVector is constructed, it is possible to
view the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.

C++ includes: Epetra_IntSerialDenseVector.h ";

/*  I/O methods  */

%feature("docstring")  Epetra_IntSerialDenseVector::Print "void
Epetra_IntSerialDenseVector::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";

/*  Expert-only unsupported methods  */

%feature("docstring")  Epetra_IntSerialDenseVector::MakeViewOf "int
Epetra_IntSerialDenseVector::MakeViewOf(const
Epetra_IntSerialDenseVector &Source)

Reset an existing IntSerialDenseVector to point to another Vector.

Allows an existing IntSerialDenseVector to become a View of another
vector's data, regardless of the DataAccess mode of the Source vector.
It is assumed that the Source vector is an independent vector, and no
checking is done to verify this.

This is used by Epetra_CrsGraph in the OptimizeStorage method. It is
used so that an existing (Copy) vector can be converted to a View.
This frees up memory that CrsGraph no longer needs.

Parameters:
-----------

Source:  The IntSerialDenseVector this will become a view of.

Integer error code, set to 0 if successful.

WARNING:  This method is extremely dangerous and should only be used
by experts. ";

%feature("docstring")
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector()

Default constructor; defines a zero size object.

Epetra_IntSerialDenseVector objects defined by the default constructor
should be sized with the Size() or Resize functions. Values should be
defined by using the [] or () operators. ";

%feature("docstring")
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(int
Length_in)

Sized constructor; defines a variable-sized object.

Parameters:
-----------

In:  Length - Length of vector.

Epetra_IntSerialDenseVector objects defined by the sized constructor
are already sized to the dimension given as a parameter. All values
are initialized to 0. Calling this constructor is equivalent to using
the default constructor, and then calling the Size function on it.
Values should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(Epetra_DataAccess
CV_in, int *Values_in, int Length_in)

Set object values from one-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Values - Pointer to an array of integer numbers containing the
values.

In:  Length - Length of vector.

See Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(const
Epetra_IntSerialDenseVector &Source)

Epetra_IntSerialDenseVector copy constructor. ";

%feature("docstring")  Epetra_IntSerialDenseVector::Size "int
Epetra_IntSerialDenseVector::Size(int Length_in)

Set length of a Epetra_IntSerialDenseVector object; init values to
zero.

Parameters:
-----------

In:  Length - Length of vector object.

Allows user to define the dimension of a Epetra_IntSerialDenseVector.
This function can be called at any point after construction. Any
values that were previously in this object are destroyed and the
resized vector starts off with all zero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_IntSerialDenseVector::Resize "int
Epetra_IntSerialDenseVector::Resize(int Length_in)

Resize a Epetra_IntSerialDenseVector object.

Parameters:
-----------

In:  Length - Length of vector object.

Allows user to define the dimension of a Epetra_IntSerialDenseVector.
This function can be called at any point after construction. Any
values that were previously in this object are copied into the new
size. If the new shape is smaller than the original, the first Length
values are copied to the new vector.

Integer error code, set to 0 if successful. ";

%feature("docstring")
Epetra_IntSerialDenseVector::~Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::~Epetra_IntSerialDenseVector()

Epetra_IntSerialDenseVector destructor. ";

%feature("docstring")  Epetra_IntSerialDenseVector::Random "int
Epetra_IntSerialDenseVector::Random()

Set vector values to random numbers.

IntSerialDenseVector uses the random number generator provided by
Epetra_Util. The vector values will be set to random values on the
interval (0, 2^31 - 1).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_IntSerialDenseVector::Length "int
Epetra_IntSerialDenseVector::Length() const

Returns length of vector. ";

%feature("docstring")  Epetra_IntSerialDenseVector::Values "int*
Epetra_IntSerialDenseVector::Values()

Returns pointer to the values in vector. ";

%feature("docstring")  Epetra_IntSerialDenseVector::Values "const
int* Epetra_IntSerialDenseVector::Values() const

Returns const pointer to the values in vector. ";

%feature("docstring")  Epetra_IntSerialDenseVector::CV "Epetra_DataAccess Epetra_IntSerialDenseVector::CV() const

Returns the data access mode of the this vector. ";


// File: classEpetra__IntVector.xml
%feature("docstring") Epetra_IntVector "

Epetra_IntVector: A class for constructing and using dense integer
vectors on a parallel computer.

The Epetra_IntVector class enables the construction and use of integer
dense vectors in a distributed memory environment. The distribution of
the dense vector is determined in part by a Epetra_Comm object and a
Epetra_Map (or Epetra_LocalMap or Epetra_BlockMap).

Distributed Global vs. Replicated Local Distributed Global Vectors -
In most instances, a multi-vector will be partitioned across multiple
memory images associated with multiple processors. In this case, there
is a unique copy of each element and elements are spread across all
processors specified by the Epetra_Comm communicator.

Replicated Local Vectors - Some algorithms use vectors that are too
small to be distributed across all processors. Replicated local
vectors handle these types of situation.

Constructing Epetra_IntVectors

There are four Epetra_IntVector constructors. The first is a basic
constructor that allocates space and sets all values to zero, the
second is a copy constructor. The third and fourth constructors work
with user data. These constructors have two data access modes: Copy
mode - Allocates memory and makes a copy of the user-provided data. In
this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the vector.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  All Epetra_IntVector constructors require a map
argument that describes the layout of elements on the parallel
machine. Specifically, map is a Epetra_Map, Epetra_LocalMap or
Epetra_BlockMap object describing the desired memory layout for the
vector.

There are four different Epetra_IntVector constructors: Basic - All
values are zero.

Copy - Copy an existing vector.

Copy from or make view of user int array.

Extracting Data from Epetra_IntVectors

Once a Epetra_IntVector is constructed, it is possible to extract a
copy of the values or create a view of them.

WARNING:  ExtractView functions are extremely dangerous from a data
hiding perspective. For both ExtractView fuctions, there is a
corresponding ExtractCopy function. We strongly encourage users to
develop code using ExtractCopy functions first and only use the
ExtractView functions in a secondary optimization phase.  There are
two Extract functions: ExtractCopy - Copy values into a user-provided
array.

ExtractView - Set user-provided array to point to Epetra_IntVector
data.

WARNING:  A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is
required for all Epetra_IntVector constructors.

C++ includes: Epetra_IntVector.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_IntVector::Epetra_IntVector "Epetra_IntVector::Epetra_IntVector(const Epetra_BlockMap &Map, bool
zeroOut=true)

Basic Epetra_IntVector constuctor.

Creates a Epetra_IntVector object and, by default, fills with zero
values.

Parameters:
-----------

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

WARNING:  Note that, because Epetra_LocalMap derives from Epetra_Map
and Epetra_Map derives from Epetra_BlockMap, this constructor works
for all three types of Epetra map classes.

Parameters:
-----------

In:  zeroOut - If true then the allocated memory will be zeroed out
initialy. If false then this memory will not be touched which can be
significantly faster.

Pointer to a Epetra_IntVector. ";

%feature("docstring")  Epetra_IntVector::Epetra_IntVector "Epetra_IntVector::Epetra_IntVector(const Epetra_IntVector &Source)

Epetra_IntVector copy constructor. ";

%feature("docstring")  Epetra_IntVector::Epetra_IntVector "Epetra_IntVector::Epetra_IntVector(Epetra_DataAccess CV, const
Epetra_BlockMap &Map, int *V)

Set vector values from user array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  V - Pointer to an array of integer numbers..

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_IntVector::~Epetra_IntVector "Epetra_IntVector::~Epetra_IntVector()

Epetra_IntVector destructor. ";

/*  Post-construction modification methods  */

%feature("docstring")  Epetra_IntVector::PutValue "int
Epetra_IntVector::PutValue(int Value)

Set all elements of the vector to Value. ";

/*  Extraction methods  */

%feature("docstring")  Epetra_IntVector::ExtractCopy "int
Epetra_IntVector::ExtractCopy(int *V) const

Put vector values into user-provided array.

Parameters:
-----------

Out:  V - Pointer to memory space that will contain the vector values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_IntVector::ExtractView "int
Epetra_IntVector::ExtractView(int **V) const

Set user-provided address of V.

Parameters:
-----------

Out:  V - Address of a pointer to that will be set to point to the
values of the vector.

Integer error code, set to 0 if successful. ";

/*  Mathematical methods  */

%feature("docstring")  Epetra_IntVector::MaxValue "int
Epetra_IntVector::MaxValue()

Find maximum value.

Maximum value across all processors. ";

%feature("docstring")  Epetra_IntVector::MinValue "int
Epetra_IntVector::MinValue()

Find minimum value.

Minimum value across all processors. ";

/*  Overloaded operators  */

/*  Attribute access functions  */

%feature("docstring")  Epetra_IntVector::Values "int*
Epetra_IntVector::Values() const

Returns a pointer to an array containing the values of this vector. ";

%feature("docstring")  Epetra_IntVector::MyLength "int
Epetra_IntVector::MyLength() const

Returns the local vector length on the calling processor of vectors in
the multi-vector. ";

%feature("docstring")  Epetra_IntVector::GlobalLength "int
Epetra_IntVector::GlobalLength() const

Returns the global vector length of vectors in the multi-vector. ";

%feature("docstring")  Epetra_IntVector::GlobalLength64 "long long
Epetra_IntVector::GlobalLength64() const ";

/*  I/O methods  */

%feature("docstring")  Epetra_IntVector::Print "void
Epetra_IntVector::Print(ostream &os) const

Print method. ";


// File: classEpetra__InvOperator.xml
%feature("docstring") Epetra_InvOperator "

Epetra_InvOperator: An implementation of the Epetra_Operator class
that reverses the role of Apply() and ApplyInverse() methods.

The Epetra_InvOperator class implements Epetra_Operator using another
pre-constructed Epetra_Operator object. Once constructed, an
Epetra_InvOperator can be used as the inverse of the input operator
object as long as the appropriate Apply and ApplyInverse methods are
implemented in the original Epetra_Operator object.

C++ includes: Epetra_InvOperator.h ";

/*  Constructor  */

%feature("docstring")  Epetra_InvOperator::Epetra_InvOperator "Epetra_InvOperator::Epetra_InvOperator(Epetra_Operator *operatorIn)

Uses an Epetra_Operator instance to implement the Epetra_Operator
interface.

Facilitates the use of an Epetra_Operator instance as an inverse
operator.

Parameters:
-----------

In:  - A fully-constructed Epetra_Operator object. ";

%feature("docstring")  Epetra_InvOperator::~Epetra_InvOperator "Epetra_InvOperator::~Epetra_InvOperator()

Destructor. ";

/*  Attribute set methods  */

%feature("docstring")  Epetra_InvOperator::SetUseTranspose "int
Epetra_InvOperator::SetUseTranspose(bool UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose_in - If true, multiply by the transpose of operator,
otherwise just use operator.

WARNING:  - This method has no effect and returns -1 as error code. ";

/*  Mathematical functions  */

%feature("docstring")  Epetra_InvOperator::Apply "int
Epetra_InvOperator::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_InvOperator applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

WARNING:  - This method has no effect and returns -1 as error code. ";

%feature("docstring")  Epetra_InvOperator::ApplyInverse "int
Epetra_InvOperator::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_InvOperator inverse applied to an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_InvOperator::NormInf "double
Epetra_InvOperator::NormInf() const

Returns the infinity norm of the global matrix. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_InvOperator::Label "const char*
Epetra_InvOperator::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_InvOperator::Operator "Epetra_Operator*
Epetra_InvOperator::Operator() const

Returns a pointer to the Epetra_Operator operator object that was used
to create this Epetra_InvOperator object. ";

%feature("docstring")  Epetra_InvOperator::UseTranspose "bool
Epetra_InvOperator::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_InvOperator::HasNormInf "bool
Epetra_InvOperator::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Epetra_InvOperator::Comm "const Epetra_Comm&
Epetra_InvOperator::Comm() const

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Epetra_InvOperator::OperatorDomainMap "const
Epetra_Map& Epetra_InvOperator::OperatorDomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  Epetra_InvOperator::OperatorRangeMap "const
Epetra_Map& Epetra_InvOperator::OperatorRangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";


// File: classEpetra__JadMatrix.xml
%feature("docstring") Epetra_JadMatrix "

Epetra_JadMatrix: A class for constructing matrix objects optimized
for common kernels.

The Epetra_JadMatrix class takes an existing Epetra_RowMatrix ojbect,
analyzes it and builds a jagged diagonal equivalent of it. Once
constructed, it is also possible to update the values of the matrix
with values from another Epetra_RowMatrix that has the identical
structure.

C++ includes: Epetra_JadMatrix.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_JadMatrix::Epetra_JadMatrix "Epetra_JadMatrix::Epetra_JadMatrix(const Epetra_RowMatrix &Matrix)

Epetra_JadMatrix constuctor. ";

%feature("docstring")  Epetra_JadMatrix::~Epetra_JadMatrix "Epetra_JadMatrix::~Epetra_JadMatrix()

Epetra_JadMatrix Destructor. ";

/*  Post-construction modifications  */

%feature("docstring")  Epetra_JadMatrix::UpdateValues "int
Epetra_JadMatrix::UpdateValues(const Epetra_RowMatrix &Matrix, bool
CheckStructure=false)

Update values using a matrix with identical structure. ";

/*  Methods required for implementing Epetra_BasicRowMatrix  */

%feature("docstring")  Epetra_JadMatrix::ExtractMyRowCopy "int
Epetra_JadMatrix::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

MyRow:  (In) - Local row to extract.

Length:  (In) - Length of Values and Indices.

NumEntries:  (Out) - Number of nonzero entries extracted.

Values:  (Out) - Extracted values for this row.

Indices:  (Out) - Extracted global column indices for the
corresponding values.

Integer error code, set to 0 if successful, set to -1 if MyRow not
valid, -2 if Length is too short (NumEntries will have required
length). ";

%feature("docstring")  Epetra_JadMatrix::ExtractMyEntryView "int
Epetra_JadMatrix::ExtractMyEntryView(int CurEntry, double *&Value, int
&RowIndex, int &ColIndex)

Returns a reference to the ith entry in the matrix, along with its row
and column index.

Parameters:
-----------

CurEntry:  (In) - Local entry to extract.

Value:  (Out) - Extracted reference to current values.

RowIndex:  (Out) - Row index for current entry.

ColIndex:  (Out) - Column index for current entry.

Integer error code, set to 0 if successful, set to -1 if CurEntry not
valid. ";

%feature("docstring")  Epetra_JadMatrix::ExtractMyEntryView "int
Epetra_JadMatrix::ExtractMyEntryView(int CurEntry, double const
*&Value, int &RowIndex, int &ColIndex) const

Returns a const reference to the ith entry in the matrix, along with
its row and column index.

Parameters:
-----------

CurEntry:  (In) - Local entry to extract.

Value:  (Out) - Extracted reference to current values.

RowIndex:  (Out) - Row index for current entry.

ColIndex:  (Out) - Column index for current entry.

Integer error code, set to 0 if successful, set to -1 if CurEntry not
valid. ";

%feature("docstring")  Epetra_JadMatrix::NumMyRowEntries "int
Epetra_JadMatrix::NumMyRowEntries(int MyRow, int &NumEntries) const

Return the current number of values stored for the specified local
row.

Similar to NumMyEntries() except NumEntries is returned as an argument
and error checking is done on the input value MyRow.

Parameters:
-----------

MyRow:  - (In) Local row.

NumEntries:  - (Out) Number of nonzero values.

Integer error code, set to 0 if successful, set to -1 if MyRow not
valid.

None.

Unchanged. ";

/*  Computational methods  */

%feature("docstring")  Epetra_JadMatrix::Multiply "int
Epetra_JadMatrix::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_JadMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  TransA -If true, multiply by the transpose of matrix, otherwise
just use matrix.

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_JadMatrix::Solve "int
Epetra_JadMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_JadMatrix solve with a
Epetra_MultiVector X in Y (not implemented).

Parameters:
-----------

In:  Upper -If true, solve Ux = y, otherwise solve Lx = y.

In:  Trans -If true, solve transpose problem.

In:  UnitDiagonal -If true, assume diagonal is unit (whether it's
stored or not).

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";


// File: classEpetra__LAPACK.xml
%feature("docstring") Epetra_LAPACK "

Epetra_LAPACK: The Epetra LAPACK Wrapper Class.

The Epetra_LAPACK class is a wrapper that encapsulates LAPACK (Linear
Algebra Package). LAPACK provides portable, high- performance
implementations of linear, eigen, SVD, etc solvers.

The standard LAPACK interface is Fortran-specific. Unfortunately, the
interface between C++ and Fortran is not standard across all computer
platforms. The Epetra_LAPACK class provides C++ wrappers for the
LAPACK kernels in order to insulate the rest of Epetra from the
details of C++ to Fortran translation. A Epetra_LAPACK object is
essentially nothing, but allows access to the LAPACK wrapper
functions.

Epetra_LAPACK is a serial interface only. This is appropriate since
the standard LAPACK are only specified for serial execution (or shared
memory parallel).

C++ includes: Epetra_LAPACK.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_LAPACK::Epetra_LAPACK "Epetra_LAPACK::Epetra_LAPACK(void)

Epetra_LAPACK Constructor.

Builds an instance of a serial LAPACK object. ";

%feature("docstring")  Epetra_LAPACK::Epetra_LAPACK "Epetra_LAPACK::Epetra_LAPACK(const Epetra_LAPACK &LAPACK)

Epetra_LAPACK Copy Constructor.

Makes an exact copy of an existing Epetra_LAPACK instance. ";

%feature("docstring")  Epetra_LAPACK::~Epetra_LAPACK "Epetra_LAPACK::~Epetra_LAPACK(void)

Epetra_LAPACK Destructor. ";

/*  Symmetric Positive Definite linear system routines  */

%feature("docstring")  Epetra_LAPACK::POTRF "void
Epetra_LAPACK::POTRF(const char UPLO, const int N, float *A, const int
LDA, int *INFO) const

Epetra_LAPACK factorization for positive definite matrix (SPOTRF) ";

%feature("docstring")  Epetra_LAPACK::POTRF "void
Epetra_LAPACK::POTRF(const char UPLO, const int N, double *A, const
int LDA, int *INFO) const

Epetra_LAPACK factorization for positive definite matrix (DPOTRF) ";

%feature("docstring")  Epetra_LAPACK::POTRS "void
Epetra_LAPACK::POTRS(const char UPLO, const int N, const int NRHS,
const float *A, const int LDA, float *X, const int LDX, int *INFO)
const

Epetra_LAPACK solve (after factorization) for positive definite matrix
(SPOTRS) ";

%feature("docstring")  Epetra_LAPACK::POTRS "void
Epetra_LAPACK::POTRS(const char UPLO, const int N, const int NRHS,
const double *A, const int LDA, double *X, const int LDX, int *INFO)
const

Epetra_LAPACK solve (after factorization) for positive definite matrix
(DPOTRS) ";

%feature("docstring")  Epetra_LAPACK::POTRI "void
Epetra_LAPACK::POTRI(const char UPLO, const int N, float *A, const int
LDA, int *INFO) const

Epetra_LAPACK inversion for positive definite matrix (SPOTRI) ";

%feature("docstring")  Epetra_LAPACK::POTRI "void
Epetra_LAPACK::POTRI(const char UPLO, const int N, double *A, const
int LDA, int *INFO) const

Epetra_LAPACK inversion for positive definite matrix (DPOTRI) ";

%feature("docstring")  Epetra_LAPACK::POCON "void
Epetra_LAPACK::POCON(const char UPLO, const int N, const float *A,
const int LDA, const float ANORM, float *RCOND, float *WORK, int
*IWORK, int *INFO) const

Epetra_LAPACK condition number estimator for positive definite matrix
(SPOCON) ";

%feature("docstring")  Epetra_LAPACK::POCON "void
Epetra_LAPACK::POCON(const char UPLO, const int N, const double *A,
const int LDA, const double ANORM, double *RCOND, double *WORK, int
*IWORK, int *INFO) const

Epetra_LAPACK condition number estimator for positive definite matrix
(DPOCON) ";

%feature("docstring")  Epetra_LAPACK::POSV "void
Epetra_LAPACK::POSV(const char UPLO, const int N, const int NRHS,
float *A, const int LDA, float *X, const int LDX, int *INFO) const

Epetra_LAPACK factor and solve for positive definite matrix (SPOSV) ";

%feature("docstring")  Epetra_LAPACK::POSV "void
Epetra_LAPACK::POSV(const char UPLO, const int N, const int NRHS,
double *A, const int LDA, double *X, const int LDX, int *INFO) const

Epetra_LAPACK factor and solve for positive definite matrix (DPOSV) ";

%feature("docstring")  Epetra_LAPACK::POEQU "void
Epetra_LAPACK::POEQU(const int N, const float *A, const int LDA, float
*S, float *SCOND, float *AMAX, int *INFO) const

Epetra_LAPACK equilibration for positive definite matrix (SPOEQU) ";

%feature("docstring")  Epetra_LAPACK::POEQU "void
Epetra_LAPACK::POEQU(const int N, const double *A, const int LDA,
double *S, double *SCOND, double *AMAX, int *INFO) const

Epetra_LAPACK equilibration for positive definite matrix (DPOEQU) ";

%feature("docstring")  Epetra_LAPACK::PORFS "void
Epetra_LAPACK::PORFS(const char UPLO, const int N, const int NRHS,
const float *A, const int LDA, const float *AF, const int LDAF, const
float *B, const int LDB, float *X, const int LDX, float *FERR, float
*BERR, float *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK solve driver for positive definite matrix (SPOSVX) ";

%feature("docstring")  Epetra_LAPACK::PORFS "void
Epetra_LAPACK::PORFS(const char UPLO, const int N, const int NRHS,
const double *A, const int LDA, const double *AF, const int LDAF,
const double *B, const int LDB, double *X, const int LDX, double
*FERR, double *BERR, double *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK solve driver for positive definite matrix (DPOSVX) ";

%feature("docstring")  Epetra_LAPACK::POSVX "void
Epetra_LAPACK::POSVX(const char FACT, const char UPLO, const int N,
const int NRHS, float *A, const int LDA, float *AF, const int LDAF,
const char EQUED, float *S, float *B, const int LDB, float *X, const
int LDX, float *RCOND, float *FERR, float *BERR, float *WORK, int
*IWORK, int *INFO) const

Epetra_LAPACK solve driver for positive definite matrix (SPOSVX) ";

%feature("docstring")  Epetra_LAPACK::POSVX "void
Epetra_LAPACK::POSVX(const char FACT, const char UPLO, const int N,
const int NRHS, double *A, const int LDA, double *AF, const int LDAF,
const char EQUED, double *S, double *B, const int LDB, double *X,
const int LDX, double *RCOND, double *FERR, double *BERR, double
*WORK, int *IWORK, int *INFO) const

Epetra_LAPACK solve driver for positive definite matrix (DPOSVX) ";

/*  General linear system routines  */

%feature("docstring")  Epetra_LAPACK::GELS "void
Epetra_LAPACK::GELS(const char TRANS, const int M, const int N, const
int NRHS, double *A, const int LDA, double *B, const int LDB, double
*WORK, const int LWORK, int *INFO) const

Epetra_LAPACK simple driver to solve least-squares systems. ";

%feature("docstring")  Epetra_LAPACK::GETRF "void
Epetra_LAPACK::GETRF(const int M, const int N, float *A, const int
LDA, int *IPIV, int *INFO) const

Epetra_LAPACK factorization for general matrix (SGETRF) ";

%feature("docstring")  Epetra_LAPACK::GETRF "void
Epetra_LAPACK::GETRF(const int M, const int N, double *A, const int
LDA, int *IPIV, int *INFO) const

Epetra_LAPACK factorization for general matrix (DGETRF) ";

%feature("docstring")  Epetra_LAPACK::GEQRF "void
Epetra_LAPACK::GEQRF(const int M, const int N, float *A, const int
LDA, float *TAU, float *WORK, const int lwork, int *INFO) const

Epetra_LAPACK QR factorization for general matrix (SGEQRF) ";

%feature("docstring")  Epetra_LAPACK::GEQRF "void
Epetra_LAPACK::GEQRF(const int M, const int N, double *A, const int
LDA, double *TAU, double *WORK, const int lwork, int *INFO) const

Epetra_LAPACK factorization for general matrix (DGEQRF) ";

%feature("docstring")  Epetra_LAPACK::GETRS "void
Epetra_LAPACK::GETRS(const char TRANS, const int N, const int NRHS,
const float *A, const int LDA, const int *IPIV, float *X, const int
LDX, int *INFO) const

Epetra_LAPACK solve (after factorization) for general matrix (SGETRS)
";

%feature("docstring")  Epetra_LAPACK::GETRS "void
Epetra_LAPACK::GETRS(const char TRANS, const int N, const int NRHS,
const double *A, const int LDA, const int *IPIV, double *X, const int
LDX, int *INFO) const

Epetra_LAPACK solve (after factorization) for general matrix (DGETRS)
";

%feature("docstring")  Epetra_LAPACK::GETRI "void
Epetra_LAPACK::GETRI(const int N, float *A, const int LDA, int *IPIV,
float *WORK, const int *LWORK, int *INFO) const

Epetra_LAPACK inversion for general matrix (SGETRI) ";

%feature("docstring")  Epetra_LAPACK::GETRI "void
Epetra_LAPACK::GETRI(const int N, double *A, const int LDA, int *IPIV,
double *WORK, const int *LWORK, int *INFO) const

Epetra_LAPACK inversion for general matrix (DGETRI) ";

%feature("docstring")  Epetra_LAPACK::GECON "void
Epetra_LAPACK::GECON(const char NORM, const int N, const float *A,
const int LDA, const float ANORM, float *RCOND, float *WORK, int
*IWORK, int *INFO) const

Epetra_LAPACK condition number estimator for general matrix (SGECON)
";

%feature("docstring")  Epetra_LAPACK::GECON "void
Epetra_LAPACK::GECON(const char NORM, const int N, const double *A,
const int LDA, const double ANORM, double *RCOND, double *WORK, int
*IWORK, int *INFO) const

Epetra_LAPACK condition number estimator for general matrix (DGECON)
";

%feature("docstring")  Epetra_LAPACK::GESV "void
Epetra_LAPACK::GESV(const int N, const int NRHS, float *A, const int
LDA, int *IPIV, float *X, const int LDX, int *INFO) const

Epetra_LAPACK factor and solve for general matrix (SGESV) ";

%feature("docstring")  Epetra_LAPACK::GESV "void
Epetra_LAPACK::GESV(const int N, const int NRHS, double *A, const int
LDA, int *IPIV, double *X, const int LDX, int *INFO) const

Epetra_LAPACK factor and solve for general matrix (DGESV) ";

%feature("docstring")  Epetra_LAPACK::GEEQU "void
Epetra_LAPACK::GEEQU(const int M, const int N, const float *A, const
int LDA, float *R, float *C, float *ROWCND, float *COLCND, float
*AMAX, int *INFO) const

Epetra_LAPACK equilibration for general matrix (SGEEQU) ";

%feature("docstring")  Epetra_LAPACK::GEEQU "void
Epetra_LAPACK::GEEQU(const int M, const int N, const double *A, const
int LDA, double *R, double *C, double *ROWCND, double *COLCND, double
*AMAX, int *INFO) const

Epetra_LAPACK equilibration for general matrix (DGEEQU) ";

%feature("docstring")  Epetra_LAPACK::GERFS "void
Epetra_LAPACK::GERFS(const char TRANS, const int N, const int NRHS,
const float *A, const int LDA, const float *AF, const int LDAF, const
int *IPIV, const float *B, const int LDB, float *X, const int LDX,
float *FERR, float *BERR, float *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK Refine solution (GERFS) ";

%feature("docstring")  Epetra_LAPACK::GERFS "void
Epetra_LAPACK::GERFS(const char TRANS, const int N, const int NRHS,
const double *A, const int LDA, const double *AF, const int LDAF,
const int *IPIV, const double *B, const int LDB, double *X, const int
LDX, double *FERR, double *BERR, double *WORK, int *IWORK, int *INFO)
const

Epetra_LAPACK Refine solution (GERFS) ";

%feature("docstring")  Epetra_LAPACK::GESVX "void
Epetra_LAPACK::GESVX(const char FACT, const char TRANS, const int N,
const int NRHS, float *A, const int LDA, float *AF, const int LDAF,
int *IPIV, const char EQUED, float *R, float *C, float *B, const int
LDB, float *X, const int LDX, float *RCOND, float *FERR, float *BERR,
float *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK solve driver for general matrix (SGESVX) ";

%feature("docstring")  Epetra_LAPACK::GESVX "void
Epetra_LAPACK::GESVX(const char FACT, const char TRANS, const int N,
const int NRHS, double *A, const int LDA, double *AF, const int LDAF,
int *IPIV, const char EQUED, double *R, double *C, double *B, const
int LDB, double *X, const int LDX, double *RCOND, double *FERR, double
*BERR, double *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK solve driver for general matrix (DGESVX) ";

%feature("docstring")  Epetra_LAPACK::GEHRD "void
Epetra_LAPACK::GEHRD(const int N, const int ILO, const int IHI, float
*A, const int LDA, float *TAU, float *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper for reduction to Hessenberg form (SGEHRD) ";

%feature("docstring")  Epetra_LAPACK::GEHRD "void
Epetra_LAPACK::GEHRD(const int N, const int ILO, const int IHI, double
*A, const int LDA, double *TAU, double *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper for reduction to Hessenberg form (DGEHRD) ";

/*  Hessenberg routines  */

%feature("docstring")  Epetra_LAPACK::HSEQR "void
Epetra_LAPACK::HSEQR(const char JOB, const char COMPZ, const int N,
const int ILO, const int IHI, float *H, const int LDH, float *WR,
float *WI, float *Z, const int LDZ, float *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper for computing the eigenvalues of a real upper
Hessenberg matrix (SHSEQR) ";

%feature("docstring")  Epetra_LAPACK::HSEQR "void
Epetra_LAPACK::HSEQR(const char JOB, const char COMPZ, const int N,
const int ILO, const int IHI, double *H, const int LDH, double *WR,
double *WI, double *Z, const int LDZ, double *WORK, const int LWORK,
int *INFO) const

Epetra_LAPACK wrapper for computing the eigenvalues of a real upper
Hessenberg matrix (DHSEQR) ";

/*  Orthogonal matrix routines  */

%feature("docstring")  Epetra_LAPACK::ORGQR "void
Epetra_LAPACK::ORGQR(const int M, const int N, const int K, float *A,
const int LDA, float *TAU, float *WORK, const int LWORK, int *INFO)
const

Epetra_LAPACK wrapper for generating a m x n real matrix Q with
orthonormal columns, defined as the product of k elementary
reflectors. (SORGQR) ";

%feature("docstring")  Epetra_LAPACK::ORGQR "void
Epetra_LAPACK::ORGQR(const int M, const int N, const int K, double *A,
const int LDA, double *TAU, double *WORK, const int LWORK, int *INFO)
const

Epetra_LAPACK wrapper for generating a m x n real matrix Q with
orthonormal columns, defined as the product of k elementary
reflectors. (DORGQR) ";

%feature("docstring")  Epetra_LAPACK::ORGHR "void
Epetra_LAPACK::ORGHR(const int N, const int ILO, const int IHI, float
*A, const int LDA, float *TAU, float *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper for generating a real orthogonal matrix Q
defined by elementary reflectors. (SORGHR) ";

%feature("docstring")  Epetra_LAPACK::ORGHR "void
Epetra_LAPACK::ORGHR(const int N, const int ILO, const int IHI, double
*A, const int LDA, double *TAU, double *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper for generating a real orthogonal matrix Q
defined by elementary reflectors. (DORGHR) ";

%feature("docstring")  Epetra_LAPACK::ORMHR "void
Epetra_LAPACK::ORMHR(const char SIDE, const char TRANS, const int M,
const int N, const int ILO, const int IHI, const float *A, const int
LDA, const float *TAU, float *C, const int LDC, float *WORK, const int
LWORK, int *INFO) const

Epetra_LAPACK wrapper for applying an orthogonal matrix in-place
(SORMHR) ";

%feature("docstring")  Epetra_LAPACK::ORMHR "void
Epetra_LAPACK::ORMHR(const char SIDE, const char TRANS, const int M,
const int N, const int ILO, const int IHI, const double *A, const int
LDA, const double *TAU, double *C, const int LDC, double *WORK, const
int LWORK, int *INFO) const

Epetra_LAPACK wrapper for applying an orthogonal matrix in-place
(DORMHR) ";

%feature("docstring")  Epetra_LAPACK::LARFT "void
Epetra_LAPACK::LARFT(const char DIRECT, const char STOREV, const int
N, const int K, double *V, const int LDV, double *TAU, double *T,
const int LDT) const

Epetra_LAPACK for forming the triangular factor of a product of
elementary Householder reflectors (SLARFT). ";

%feature("docstring")  Epetra_LAPACK::LARFT "void
Epetra_LAPACK::LARFT(const char DIRECT, const char STOREV, const int
N, const int K, float *V, const int LDV, float *TAU, float *T, const
int LDT) const

Epetra_LAPACK for forming the triangular factor of a product of
elementary Householder reflectors (DLARFT). ";

/*  Triangular matrix routines  */

%feature("docstring")  Epetra_LAPACK::TREVC "void
Epetra_LAPACK::TREVC(const char SIDE, const char HOWMNY, int *SELECT,
const int N, const float *T, const int LDT, float *VL, const int LDVL,
float *VR, const int LDVR, const int MM, int *M, float *WORK, int
*INFO) const

Epetra_LAPACK wrapper for computing eigenvectors of a quasi-
triangular/triagnular matrix (STREVC)

WARNING:  HOWMNY = 'S\" is not supported. ";

%feature("docstring")  Epetra_LAPACK::TREVC "void
Epetra_LAPACK::TREVC(const char SIDE, const char HOWMNY, int *SELECT,
const int N, const double *T, const int LDT, double *VL, const int
LDVL, double *VR, const int LDVR, const int MM, int *M, double *WORK,
int *INFO) const

Epetra_LAPACK wrapper for computing eigenvectors of a quasi-
triangular/triagnular matrix (DTREVC)

WARNING:  HOWMNY = 'S\" is not supported. ";

%feature("docstring")  Epetra_LAPACK::TREXC "void
Epetra_LAPACK::TREXC(const char COMPQ, const int N, float *T, const
int LDT, float *Q, const int LDQ, int IFST, int ILST, float *WORK, int
*INFO) const

Epetra_LAPACK wrapper for reordering the real-Schur/Schur
factorization of a matrix (STREXC) ";

%feature("docstring")  Epetra_LAPACK::TREXC "void
Epetra_LAPACK::TREXC(const char COMPQ, const int N, double *T, const
int LDT, double *Q, const int LDQ, int IFST, int ILST, double *WORK,
int *INFO) const

Epetra_LAPACK wrapper for reordering the real-Schur/Schur
factorization of a matrix (DTREXC) ";

/*  Singular Value Decomposition matrix routines  */

%feature("docstring")  Epetra_LAPACK::GESVD "void
Epetra_LAPACK::GESVD(const char JOBU, const char JOBVT, const int M,
const int N, float *A, const int LDA, float *S, float *U, const int
LDU, float *VT, const int LDVT, float *WORK, const int *LWORK, int
*INFO) const

Epetra_LAPACK wrapper for computing the singular value decomposition
(SGESVD) ";

%feature("docstring")  Epetra_LAPACK::GESVD "void
Epetra_LAPACK::GESVD(const char JOBU, const char JOBVT, const int M,
const int N, double *A, const int LDA, double *S, double *U, const int
LDU, double *VT, const int LDVT, double *WORK, const int *LWORK, int
*INFO) const

Epetra_LAPACK wrapper for computing the singular value decomposition
(DGESVD) ";

%feature("docstring")  Epetra_LAPACK::GGSVD "void
Epetra_LAPACK::GGSVD(const char JOBU, const char JOBV, const char
JOBQ, const int M, const int N, const int P, int *K, int *L, double
*A, const int LDA, double *B, const int LDB, double *ALPHA, double
*BETA, double *U, const int LDU, double *V, const int LDV, double *Q,
const int LDQ, double *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK wrapper to compute the generalized singular value
decomposition (GSVD) of an M-by-N real matrix A and P-by-N real matrix
B. ";

%feature("docstring")  Epetra_LAPACK::GGSVD "void
Epetra_LAPACK::GGSVD(const char JOBU, const char JOBV, const char
JOBQ, const int M, const int N, const int P, int *K, int *L, float *A,
const int LDA, float *B, const int LDB, float *ALPHA, float *BETA,
float *U, const int LDU, float *V, const int LDV, float *Q, const int
LDQ, float *WORK, int *IWORK, int *INFO) const

Epetra_LAPACK wrapper to compute the generalized singular value
decomposition (GSVD) of an M-by-N real matrix A and P-by-N real matrix
B. ";

/*  Eigenvalue/Eigenvector routines  */

%feature("docstring")  Epetra_LAPACK::GEEV "void
Epetra_LAPACK::GEEV(const char JOBVL, const char JOBVR, const int N,
double *A, const int LDA, double *WR, double *WI, double *VL, const
int LDVL, double *VR, const int LDVR, double *WORK, const int LWORK,
int *INFO) const

Epetra_LAPACK wrapper to compute for an N-by-N real nonsymmetric
matrix A, the eigenvalues and, optionally, the left and/or right
eigenvectors. ";

%feature("docstring")  Epetra_LAPACK::GEEV "void
Epetra_LAPACK::GEEV(const char JOBVL, const char JOBVR, const int N,
float *A, const int LDA, float *WR, float *WI, float *VL, const int
LDVL, float *VR, const int LDVR, float *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper to compute for an N-by-N real nonsymmetric
matrix A, the eigenvalues and, optionally, the left and/or right
eigenvectors. ";

%feature("docstring")  Epetra_LAPACK::SPEV "void
Epetra_LAPACK::SPEV(const char JOBZ, const char UPLO, const int N,
double *AP, double *W, double *Z, int LDZ, double *WORK, int *INFO)
const

Epetra_LAPACK wrapper to compute all the eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A in packed storage. ";

%feature("docstring")  Epetra_LAPACK::SPEV "void
Epetra_LAPACK::SPEV(const char JOBZ, const char UPLO, const int N,
float *AP, float *W, float *Z, int LDZ, float *WORK, int *INFO) const

Epetra_LAPACK wrapper to compute all the eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A in packed storage. ";

%feature("docstring")  Epetra_LAPACK::SPGV "void
Epetra_LAPACK::SPGV(const int ITYPE, const char JOBZ, const char UPLO,
const int N, double *AP, double *BP, double *W, double *Z, const int
LDZ, double *WORK, int *INFO) const

Epetra_LAPACK wrapper to compute all the eigenvalues and, optionally,
the eigenvectors of a real generalized symmetric-definite
eigenproblem, of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
B*A*x=(lambda)*x. ";

%feature("docstring")  Epetra_LAPACK::SPGV "void
Epetra_LAPACK::SPGV(const int ITYPE, const char JOBZ, const char UPLO,
const int N, float *AP, float *BP, float *W, float *Z, const int LDZ,
float *WORK, int *INFO) const

Epetra_LAPACK wrapper to compute all the eigenvalues and, optionally,
the eigenvectors of a real generalized symmetric-definite
eigenproblem, of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
B*A*x=(lambda)*x. ";

%feature("docstring")  Epetra_LAPACK::SYEV "void
Epetra_LAPACK::SYEV(const char JOBZ, const char UPLO, const int N,
double *A, const int LDA, double *W, double *WORK, const int LWORK,
int *INFO) const

Epetra_LAPACK wrapper to compute all eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A. ";

%feature("docstring")  Epetra_LAPACK::SYEV "void
Epetra_LAPACK::SYEV(const char JOBZ, const char UPLO, const int N,
float *A, const int LDA, float *W, float *WORK, const int LWORK, int
*INFO) const

Epetra_LAPACK wrapper to compute all eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A. ";

%feature("docstring")  Epetra_LAPACK::SYEVD "void
Epetra_LAPACK::SYEVD(const char JOBZ, const char UPLO, const int N,
double *A, const int LDA, double *W, double *WORK, const int LWORK,
int *IWORK, const int LIWORK, int *INFO) const

Epetra_LAPACK wrapper to compute all eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A. ";

%feature("docstring")  Epetra_LAPACK::SYEVD "void
Epetra_LAPACK::SYEVD(const char JOBZ, const char UPLO, const int N,
float *A, const int LDA, float *W, float *WORK, const int LWORK, int
*IWORK, const int LIWORK, int *INFO) const

Epetra_LAPACK wrapper to compute all eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A. ";

%feature("docstring")  Epetra_LAPACK::SYEVX "void
Epetra_LAPACK::SYEVX(const char JOBZ, const char RANGE, const char
UPLO, const int N, double *A, const int LDA, const double *VL, const
double *VU, const int *IL, const int *IU, const double ABSTOL, int *M,
double *W, double *Z, const int LDZ, double *WORK, const int LWORK,
int *IWORK, int *IFAIL, int *INFO) const

Epetra_LAPACK wrapper to compute selected eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A. ";

%feature("docstring")  Epetra_LAPACK::SYEVX "void
Epetra_LAPACK::SYEVX(const char JOBZ, const char RANGE, const char
UPLO, const int N, float *A, const int LDA, const float *VL, const
float *VU, const int *IL, const int *IU, const float ABSTOL, int *M,
float *W, float *Z, const int LDZ, float *WORK, const int LWORK, int
*IWORK, int *IFAIL, int *INFO) const

Epetra_LAPACK wrapper to compute selected eigenvalues and, optionally,
eigenvectors of a real symmetric matrix A. ";

%feature("docstring")  Epetra_LAPACK::SYGV "void
Epetra_LAPACK::SYGV(const int ITYPE, const char JOBZ, const char UPLO,
const int N, double *A, const int LDA, double *B, const int LDB,
double *W, double *WORK, const int LWORK, int *INFO) const

Epetra_LAPACK wrapper to compute all the eigenvalues, and optionally,
the eigenvectors of a real generalized symmetric-definite
eigenproblem, of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
B*A*x=(lambda)*x. ";

%feature("docstring")  Epetra_LAPACK::SYGV "void
Epetra_LAPACK::SYGV(const int ITYPE, const char JOBZ, const char UPLO,
const int N, float *A, const int LDA, float *B, const int LDB, float
*W, float *WORK, const int LWORK, int *INFO) const

Epetra_LAPACK wrapper to compute all the eigenvalues, and optionally,
the eigenvectors of a real generalized symmetric-definite
eigenproblem, of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
B*A*x=(lambda)*x. ";

%feature("docstring")  Epetra_LAPACK::SYGVX "void
Epetra_LAPACK::SYGVX(const int ITYPE, const char JOBZ, const char
RANGE, const char UPLO, const int N, double *A, const int LDA, double
*B, const int LDB, const double *VL, const double *VU, const int *IL,
const int *IU, const double ABSTOL, int *M, double *W, double *Z,
const int LDZ, double *WORK, const int LWORK, int *IWORK, int *IFAIL,
int *INFO) const

Epetra_LAPACK wrapper to compute selected eigenvalues, and optionally,
eigenvectors of a real generalized symmetric-definite eigenproblem, of
the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. ";

%feature("docstring")  Epetra_LAPACK::SYGVX "void
Epetra_LAPACK::SYGVX(const int ITYPE, const char JOBZ, const char
RANGE, const char UPLO, const int N, float *A, const int LDA, float
*B, const int LDB, const float *VL, const float *VU, const int *IL,
const int *IU, const float ABSTOL, int *M, float *W, float *Z, const
int LDZ, float *WORK, const int LWORK, int *IWORK, int *IFAIL, int
*INFO) const

Epetra_LAPACK wrapper to compute selected eigenvalues, and optionally,
eigenvectors of a real generalized symmetric-definite eigenproblem, of
the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x. ";

%feature("docstring")  Epetra_LAPACK::SYEVR "void
Epetra_LAPACK::SYEVR(const char JOBZ, const char RANGE, const char
UPLO, const int N, double *A, const int LDA, const double *VL, const
double *VU, const int *IL, const int *IU, const double ABSTOL, int *M,
double *W, double *Z, const int LDZ, int *ISUPPZ, double *WORK, const
int LWORK, int *IWORK, const int LIWORK, int *INFO) const

Epetra_LAPACK wrapper to compute selected eigenvalues and, optionally,
eigenvectors of a real symmetric matrix T. ";

%feature("docstring")  Epetra_LAPACK::SYEVR "void
Epetra_LAPACK::SYEVR(const char JOBZ, const char RANGE, const char
UPLO, const int N, float *A, const int LDA, const float *VL, const
float *VU, const int *IL, const int *IU, const float ABSTOL, int *M,
float *W, float *Z, const int LDZ, int *ISUPPZ, float *WORK, const int
LWORK, int *IWORK, const int LIWORK, int *INFO) const

Epetra_LAPACK wrapper to compute selected eigenvalues and, optionally,
eigenvectors of a real symmetric matrix T. ";

%feature("docstring")  Epetra_LAPACK::GEEVX "void
Epetra_LAPACK::GEEVX(const char BALANC, const char JOBVL, const char
JOBVR, const char SENSE, const int N, double *A, const int LDA, double
*WR, double *WI, double *VL, const int LDVL, double *VR, const int
LDVR, int *ILO, int *IHI, double *SCALE, double *ABNRM, double
*RCONDE, double *RCONDV, double *WORK, const int LWORK, int *IWORK,
int *INFO) const

Epetra_LAPACK wrapper to compute for an N-by-N real nonsymmetric
matrix A, the eigenvalues and, optionally, the left and/or right
eigenvectors. ";

%feature("docstring")  Epetra_LAPACK::GEEVX "void
Epetra_LAPACK::GEEVX(const char BALANC, const char JOBVL, const char
JOBVR, const char SENSE, const int N, float *A, const int LDA, float
*WR, float *WI, float *VL, const int LDVL, float *VR, const int LDVR,
int *ILO, int *IHI, float *SCALE, float *ABNRM, float *RCONDE, float
*RCONDV, float *WORK, const int LWORK, int *IWORK, int *INFO) const

Epetra_LAPACK wrapper to compute for an N-by-N real nonsymmetric
matrix A, the eigenvalues and, optionally, the left and/or right
eigenvectors. ";

%feature("docstring")  Epetra_LAPACK::GESDD "void
Epetra_LAPACK::GESDD(const char JOBZ, const int M, const int N, double
*A, const int LDA, double *S, double *U, const int LDU, double *VT,
const int LDVT, double *WORK, const int LWORK, int *IWORK, int *INFO)
const

Epetra_LAPACK wrapper to compute the singular value decomposition
(SVD) of a real M-by-N matrix A, optionally computing the left and
right singular vectors. ";

%feature("docstring")  Epetra_LAPACK::GESDD "void
Epetra_LAPACK::GESDD(const char JOBZ, const int M, const int N, float
*A, const int LDA, float *S, float *U, const int LDU, float *VT, const
int LDVT, float *WORK, const int LWORK, int *IWORK, int *INFO) const

Epetra_LAPACK wrapper to. ";

%feature("docstring")  Epetra_LAPACK::GGEV "void
Epetra_LAPACK::GGEV(const char JOBVL, const char JOBVR, const int N,
double *A, const int LDA, double *B, const int LDB, double *ALPHAR,
double *ALPHAI, double *BETA, double *VL, const int LDVL, double *VR,
const int LDVR, double *WORK, const int LWORK, int *INFO) const

Epetra_LAPACK wrapper to compute for a pair of N-by-N real
nonsymmetric matrices (A,B) the generalized eigenvalues, and
optionally, the left and/or right generalized eigenvectors. ";

%feature("docstring")  Epetra_LAPACK::GGEV "void
Epetra_LAPACK::GGEV(const char JOBVL, const char JOBVR, const int N,
float *A, const int LDA, float *B, const int LDB, float *ALPHAR, float
*ALPHAI, float *BETA, float *VL, const int LDVL, float *VR, const int
LDVR, float *WORK, const int LWORK, int *INFO) const

Epetra_LAPACK wrapper to compute for a pair of N-by-N real
nonsymmetric matrices (A,B) the generalized eigenvalues, and
optionally, the left and/or right generalized eigenvectors. ";

/*  Linear Least Squares  */

%feature("docstring")  Epetra_LAPACK::GGLSE "void
Epetra_LAPACK::GGLSE(const int M, const int N, const int P, double *A,
const int LDA, double *B, const int LDB, double *C, double *D, double
*X, double *WORK, const int LWORK, int *INFO) const

Epetra_LAPACK wrapper to solve the linear equality-constrained least
squares (LSE) problem. ";

%feature("docstring")  Epetra_LAPACK::GGLSE "void
Epetra_LAPACK::GGLSE(const int M, const int N, const int P, float *A,
const int LDA, float *B, const int LDB, float *C, float *D, float *X,
float *WORK, const int LWORK, int *INFO) const

Epetra_LAPACK wrapper to solve the linear equality-constrained least
squares (LSE) problem. ";

/*  Machine characteristics routines  */

%feature("docstring")  Epetra_LAPACK::LAMCH "void
Epetra_LAPACK::LAMCH(const char CMACH, float &T) const

Epetra_LAPACK wrapper for DLAMCH routine. On out, T holds machine
double precision floating point characteristics. This information is
returned by the Lapack routine. ";

%feature("docstring")  Epetra_LAPACK::LAMCH "void
Epetra_LAPACK::LAMCH(const char CMACH, double &T) const

Epetra_LAPACK wrapper for SLAMCH routine. On out, T holds machine
single precision floating point characteristics. This information is
returned by the Lapack routine. ";

/*  Triangular solve  */

%feature("docstring")  Epetra_LAPACK::TRTRS "void
Epetra_LAPACK::TRTRS(const char UPLO, const char TRANS, const char
DIAG, const int N, const int NRHS, const float *A, const int LDA,
float *B, const int LDB, int *INFO) const

Epetra_LAPACK wrapper for TRTRS routine. ";

%feature("docstring")  Epetra_LAPACK::TRTRS "void
Epetra_LAPACK::TRTRS(const char UPLO, const char TRANS, const char
DIAG, const int N, const int NRHS, const double *A, const int LDA,
double *B, const int LDB, int *INFO) const

Epetra_LAPACK wrapper for TRTRS routine. ";


// File: classEpetra__LinearProblem.xml
%feature("docstring") Epetra_LinearProblem "

Epetra_LinearProblem: The Epetra Linear Problem Class.

The Epetra_LinearProblem class is a wrapper that encapsulates the
general information needed for solving a linear system of equations.
Currently it accepts a Epetra matrix, initial guess and RHS and
returns the solution. the elapsed time for each calling processor.

C++ includes: Epetra_LinearProblem.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_LinearProblem::Epetra_LinearProblem "Epetra_LinearProblem::Epetra_LinearProblem(void)

Epetra_LinearProblem Default Constructor.

Creates an empty Epetra_LinearProblem instance. The operator A, left-
hand-side X and right-hand-side B must be set use the SetOperator(),
SetLHS() and SetRHS() methods respectively. ";

%feature("docstring")  Epetra_LinearProblem::Epetra_LinearProblem "Epetra_LinearProblem::Epetra_LinearProblem(Epetra_RowMatrix *A,
Epetra_MultiVector *X, Epetra_MultiVector *B)

Epetra_LinearProblem Constructor to pass in an operator as a matrix.

Creates a Epetra_LinearProblem instance where the operator is passed
in as a matrix. ";

%feature("docstring")  Epetra_LinearProblem::Epetra_LinearProblem "Epetra_LinearProblem::Epetra_LinearProblem(Epetra_Operator *A,
Epetra_MultiVector *X, Epetra_MultiVector *B)

Epetra_LinearProblem Constructor to pass in a basic Epetra_Operator.

Creates a Epetra_LinearProblem instance for the case where an operator
is not necessarily a matrix. ";

%feature("docstring")  Epetra_LinearProblem::Epetra_LinearProblem "Epetra_LinearProblem::Epetra_LinearProblem(const Epetra_LinearProblem
&Problem)

Epetra_LinearProblem Copy Constructor.

Makes copy of an existing Epetra_LinearProblem instance. ";

%feature("docstring")  Epetra_LinearProblem::~Epetra_LinearProblem "Epetra_LinearProblem::~Epetra_LinearProblem(void)

Epetra_LinearProblem Destructor.

Completely deletes a Epetra_LinearProblem object. ";

/*  Integrity check method  */

%feature("docstring")  Epetra_LinearProblem::CheckInput "int
Epetra_LinearProblem::CheckInput() const

Check input parameters for existence and size consistency.

Returns 0 if all input parameters are valid. Returns +1 if operator is
not a matrix. This is not necessarily an error, but no scaling can be
done if the user passes in an Epetra_Operator that is not an
Epetra_Matrix ";

/*  Set methods  */

%feature("docstring")  Epetra_LinearProblem::AssertSymmetric "void
Epetra_LinearProblem::AssertSymmetric() ";

%feature("docstring")  Epetra_LinearProblem::SetPDL "void
Epetra_LinearProblem::SetPDL(ProblemDifficultyLevel PDL)

Set problem difficulty level.

Sets Aztec options and parameters based on a definition of easy
moderate or hard problem. Relieves the user from explicitly setting a
large number of individual parameter values. This function can be used
in conjunction with the SetOptions() and SetParams() functions. ";

%feature("docstring")  Epetra_LinearProblem::SetOperator "void
Epetra_LinearProblem::SetOperator(Epetra_RowMatrix *A)

Set Operator A of linear problem AX = B using an Epetra_RowMatrix.

Sets a pointer to a Epetra_RowMatrix. No copy of the operator is made.
";

%feature("docstring")  Epetra_LinearProblem::SetOperator "void
Epetra_LinearProblem::SetOperator(Epetra_Operator *A)

Set Operator A of linear problem AX = B using an Epetra_Operator.

Sets a pointer to a Epetra_Operator. No copy of the operator is made.
";

%feature("docstring")  Epetra_LinearProblem::SetLHS "void
Epetra_LinearProblem::SetLHS(Epetra_MultiVector *X)

Set left-hand-side X of linear problem AX = B.

Sets a pointer to a Epetra_MultiVector. No copy of the object is made.
";

%feature("docstring")  Epetra_LinearProblem::SetRHS "void
Epetra_LinearProblem::SetRHS(Epetra_MultiVector *B)

Set right-hand-side B of linear problem AX = B.

Sets a pointer to a Epetra_MultiVector. No copy of the object is made.
";

/*  Computational methods  */

%feature("docstring")  Epetra_LinearProblem::LeftScale "int
Epetra_LinearProblem::LeftScale(const Epetra_Vector &D)

Perform left scaling of a linear problem.

Applies the scaling vector D to the left side of the matrix A() and to
the right hand side B(). Note that the operator must be an
Epetra_RowMatrix, not just an Epetra_Operator (the base class of
Epetra_RowMatrix).

Parameters:
-----------

In:  D - Vector containing scaling values. D[i] will be applied to the
ith row of A() and B().

Integer error code, set to 0 if successful. Return -1 if operator is
not a matrix. ";

%feature("docstring")  Epetra_LinearProblem::RightScale "int
Epetra_LinearProblem::RightScale(const Epetra_Vector &D)

Perform right scaling of a linear problem.

Applies the scaling vector D to the right side of the matrix A().
Apply the inverse of D to the initial guess. Note that the operator
must be an Epetra_RowMatrix, not just an Epetra_Operator (the base
class of Epetra_RowMatrix).

Parameters:
-----------

In:  D - Vector containing scaling values. D[i] will be applied to the
ith row of A(). 1/D[i] will be applied to the ith row of B().

Integer error code, set to 0 if successful. Return -1 if operator is
not a matrix. ";

/*  Accessor methods  */

%feature("docstring")  Epetra_LinearProblem::GetOperator "Epetra_Operator* Epetra_LinearProblem::GetOperator() const

Get a pointer to the operator A. ";

%feature("docstring")  Epetra_LinearProblem::GetMatrix "Epetra_RowMatrix* Epetra_LinearProblem::GetMatrix() const

Get a pointer to the matrix A. ";

%feature("docstring")  Epetra_LinearProblem::GetLHS "Epetra_MultiVector* Epetra_LinearProblem::GetLHS() const

Get a pointer to the left-hand-side X. ";

%feature("docstring")  Epetra_LinearProblem::GetRHS "Epetra_MultiVector* Epetra_LinearProblem::GetRHS() const

Get a pointer to the right-hand-side B. ";

%feature("docstring")  Epetra_LinearProblem::GetPDL "ProblemDifficultyLevel Epetra_LinearProblem::GetPDL() const

Get problem difficulty level. ";

%feature("docstring")  Epetra_LinearProblem::IsOperatorSymmetric "bool Epetra_LinearProblem::IsOperatorSymmetric() const

Get operator symmetry bool. ";


// File: classEpetra__LinearProblemRedistor.xml
%feature("docstring") Epetra_LinearProblemRedistor "

Epetra_LinearProblemRedistor: A class for redistributing an
Epetra_LinearProblem object.

This class provides capabilities to redistribute an existing
Epetra_LinearProblem object across a parallel distributed memory
machine. All or part of a linear problem object can be redistributed.
Reverse distributions, value updates and matrix transposition are also
supported. Specification of the redistribution can be done by
providing a target Epetra_Map object describing the new distribution,
or

by specifying the number of processors to use and stating whether or
not the problem should be completely replicated on all processors.

C++ includes: Epetra_LinearProblemRedistor.h ";

/*  Constructors/destructors  */

%feature("docstring")
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor "Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(Epetra_LinearProblem
*OrigProblem, const Epetra_Map &RedistMap)

Epetra_LinearProblemRedistor constructor using pre-defined layout.

Parameters:
-----------

Problem:  (In) An existing Epetra_LinearProblem object. The
Epetra_RowMatrix, the LHS and RHS pointers do not need to be defined
before this constructor is called.

RedistMap:  (In) An Epetra_Map describing the target layout of the
redistribution.

Pointer to a Epetra_LinearProblemRedistor object. ";

%feature("docstring")
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor "Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(Epetra_LinearProblem
*OrigProblem, int NumProc, bool Replicate)

Epetra_LinearProblemRedistor constructor specifying number of
processor and replication bool.

Parameters:
-----------

Problem:  (In) An existing Epetra_LinearProblem object. The
Epetra_RowMatrix, the LHS and RHS pointers do not need to be defined
before this constructor is called.

NumProc:  (In) Number of processors to use when redistributing the
problem. Must be between 1 and the number of processor on the parallel
machine.

Replicate:  (In) A bool that indicates if the linear problem should be
fully replicated on all processors. If true, then a complete copy of
the linear problem will be made on each processor. If false, then the
problem will be roughly evenly spread across the total number of
processors.

Pointer to a Epetra_LinearProblemRedistor object. ";

%feature("docstring")
Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor "Epetra_LinearProblemRedistor::Epetra_LinearProblemRedistor(const
Epetra_LinearProblemRedistor &Source)

Epetra_LinearProblemRedistor copy constructor. ";

%feature("docstring")
Epetra_LinearProblemRedistor::~Epetra_LinearProblemRedistor "Epetra_LinearProblemRedistor::~Epetra_LinearProblemRedistor()

Epetra_LinearProblemRedistor destructor. ";

/*  Forward transformation methods  */

%feature("docstring")
Epetra_LinearProblemRedistor::CreateRedistProblem "int
Epetra_LinearProblemRedistor::CreateRedistProblem(const bool
ConstructTranspose, const bool MakeDataContiguous,
Epetra_LinearProblem *&RedistProblem)

Generate a new Epetra_LinearProblem as a redistribution of the one
passed into the constructor.

Constructs a new Epetra_LinearProblem that is a copy of the one passed
in to the constructor. The new problem will have redistributed copies
of the RowMatrix, LHS and RHS from the original problem. If any of
these three objects are 0 pointers, then the corresponding pointer
will be zero in the redistributed object.

The redistributed matrix will constructed as an Epetra_CrsMatrix. The
LHS and RHS will be Epetra_MultiVector objects.

Two bools can be set when calling this method. The first,
ConstructTranspose, will cause the Redistribute method to construct
the transpose of the original row matrix. The second,
MakeDataContiguous, forces the memory layout of the output matrix, RHS
and LHS to be stored so that it is compatible with Fortran. In
particular, the Epetra_CrsMatrix is stored so that value from row to
row are contiguous, as are the indices. This is compatible with the
Harwell-Boeing compressed row and compressed column format. The RHS
and LHS are created so that there is a constant stride between the
columns of the multivector. This is compatible with Fortran 2D array
storage.

Parameters:
-----------

ConstructTranspose:  (In) Causes the output matrix to be transposed.
This feature can be used to support solvers that need the matrix to be
stored in column format. This option has no impact on the LHS and RHS
of the output problem.

MakeDataContiguous:  (In) Causes the output matrix, LHS and RHS to be
stored in a form compatible with Fortran-style solvers. The output
matrix will be compatible with the Harwell-Boeing compressed column
format. The RHS and LHS will be stored such that the last value in
column j of the multivector is stored next to the first value in
column j+1.

RedistProblem:  (Out) The redistributed Epetra_LinearProblem. The
RowMatrix, LHS and RHS that are generated as part of this problem will
be destroyed when the Epetra_LinearProblemRedistor object is
destroyed.

Integer error code, 0 if no errors, positive value if one or more of
the LHS or RHS pointers were 0. Negative if some other fatal error
occured. ";

%feature("docstring")
Epetra_LinearProblemRedistor::UpdateRedistProblemValues "int
Epetra_LinearProblemRedistor::UpdateRedistProblemValues(Epetra_LinearProblem
*ProblemWithNewValues)

Update the values of an already-redistributed problem.

Updates the values of an already-redistributed problem. This method
allows updating the redistributed problem without allocating new
storage. All three objects in the RedistProblem will be updated,
namely the Matrix, LHS and RHS. If the LHS or RHS are 0 pointers, they
will be ignored.

Parameters:
-----------

ProblemWithNewValues:  (In) The values from ProblemWithNewValues will
be copied into the RedistProblem. The ProblemWithNewValues object must
be identical in structure to the Epetra_LinearProblem object used to
create this instance of Epetra_LinearProblemRedistor.

Integer error code, 0 if no errors, positive value if one or more of
the input LHS or RHS pointers were 0. Negative if some other fatal
error occured. ";

%feature("docstring")  Epetra_LinearProblemRedistor::UpdateRedistRHS "int Epetra_LinearProblemRedistor::UpdateRedistRHS(Epetra_MultiVector
*RHSWithNewValues)

Update the values of an already-redistributed RHS.

Updates the values of an already-redistributed RHS. This method allows
updating the redistributed RHS without allocating new storage. This
method updates only the RHS, and no other part of the
RedistLinearProblem.

Parameters:
-----------

RHSWithNewValues:  (In) The values from RHSWithNewValues will be
copied into the RHS of the RedistProblem. The RHSWithNewValues object
must be identical in structure to the Epetra_MultiVector object used
to create this instance of Epetra_LinearProblemRedistor.

Integer error code, 0 if no errors. ";

/*  Reverse transformation methods  */

%feature("docstring")  Epetra_LinearProblemRedistor::UpdateOriginalLHS
"int
Epetra_LinearProblemRedistor::UpdateOriginalLHS(Epetra_MultiVector
*LHS)

Update LHS of original Linear Problem object.

Copies the values from the LHS of the RedistProblem Object into the
LHS passed in to the method. If the RedistProblem is replicated, the
LHS will be computed as an average from all processor.

Parameters:
-----------

LHS:  (Out) On exit, the values in LHS will contain the values from
the current RedistLinearProblem LHS. If the GIDs of the RedistMap are
not one-to-one, e.g., if the map is replicated, the output values for
each GID will be an average of all values at that GID. If the
RedistProblem is being solved redundantly in any fashion, this
approach to computing the values of LHS should produce a valid answer
no matter how the RedistProblem LHS is distributed.

Error code, returns 0 if no error. ";

/*  Attribute accessor methods  */

%feature("docstring")  Epetra_LinearProblemRedistor::RedistMap "const
Epetra_Map& Epetra_LinearProblemRedistor::RedistMap() const

Returns const reference to the Epetra_Map that describes the layout of
the RedistLinearProblem.

The RedistMap object can be used to construct other Epetra_DistObject
objects whose maps are compatible with the redistributed linear
problem map. WARNING:  This method must not be called until after
CreateRedistProblem() is called. ";

%feature("docstring")  Epetra_LinearProblemRedistor::RedistExporter "const Epetra_Export& Epetra_LinearProblemRedistor::RedistExporter()
const

Returns const reference to the Epetra_Export object used to
redistribute the original linear problem.

The RedistExporter object can be used to redistribute other
Epetra_DistObject objects whose maps are compatible with the original
linear problem map, or a reverse distribution for objects compatible
with the RedistMap(). WARNING:  This method must not be called until
after CreateRedistProblem() is called. ";

/*  Utility methods  */

%feature("docstring")  Epetra_LinearProblemRedistor::ExtractHbData "int Epetra_LinearProblemRedistor::ExtractHbData(int &M, int &N, int
&nz, int *&ptr, int *&ind, double *&val, int &Nrhs, double *&rhs, int
&ldrhs, double *&lhs, int &ldlhs) const

Extract the redistributed problem data in a form usable for other
codes that require Harwell-Boeing format.

This method extract data from the linear problem for use with other
packages, such as SuperLU, that require the matrix, rhs and lhs in
Harwell-Boeing format. Note that the arrays returned by this method
are owned by the Epetra_LinearProblemRedistor class, and they will be
deleted when the owning class is destroyed. ";

/*  I/O methods  */

%feature("docstring")  Epetra_LinearProblemRedistor::Print "virtual
void Epetra_LinearProblemRedistor::Print(ostream &os) const

Print method. ";


// File: classEpetra__LocalMap.xml
%feature("docstring") Epetra_LocalMap "

Epetra_LocalMap: A class for replicating vectors and matrices across
multiple processors.

Small matrix and vector objects are often replicated on distributed
memory parallel machines. The Epetra_LocalMap class allows
construction of these replicated local objects and keeps information
that describes this distribution.

Epetra_LocalMap allows the storage and retrieval of the following
information. Once a Epetra_Map is constructed any of the following
attributes can be obtained by calling a query function that has the
name as the attribute, e.g. to get the value of NumGlobalPoints, you
can call a function NumGlobalElements(). For attributes that are
lists, the query functions return the list values in a user allocated
array.

NumMyElements - The number of elements owned by the calling processor.

IndexBase - The base integer value for indexed array references.
Typically this is 0 for C/C++ and 1 for Fortran, but it can be set to
any integer value.

Comm - The Epetra_Comm communicator. This communicator can in turn be
queried for processor rank and size information.

The Epetra_LocalMap class is actually a derived class of Epetra_Map.
Epetra_Map is in turn derived from Epetra_BlockMap. As such,
Epetra_LocalMap has full access to all the functions in these other
map classes.

In particular, the following function allows a boolean test:

DistributedGlobal() - Returns false for a Epetra_LocalMap object.

WARNING:  A Epetra_Comm object is required for all Epetra_LocalMap
constructors.

C++ includes: Epetra_LocalMap.h ";

%feature("docstring")  Epetra_LocalMap::Epetra_LocalMap "Epetra_LocalMap::Epetra_LocalMap(int NumMyElements, int IndexBase,
const Epetra_Comm &Comm)

Epetra_LocalMap constructor for a user-defined replicate distribution
of elements.

Creates a map that puts NumMyElements on the calling processor. Each
processor should pass in the same value for NumMyElements.

Parameters:
-----------

In:  NumMyElements - Number of elements owned by the calling
processor.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_Map object. ";

%feature("docstring")  Epetra_LocalMap::Epetra_LocalMap "Epetra_LocalMap::Epetra_LocalMap(unsigned int NumMyElements, int
IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_LocalMap::Epetra_LocalMap "Epetra_LocalMap::Epetra_LocalMap(long long NumMyElements, int
IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_LocalMap::Epetra_LocalMap "Epetra_LocalMap::Epetra_LocalMap(const Epetra_LocalMap &map)

Epetra_LocalMap copy constructor. ";

%feature("docstring")  Epetra_LocalMap::~Epetra_LocalMap "Epetra_LocalMap::~Epetra_LocalMap()

Epetra_LocalMap destructor. ";


// File: classEpetra__LongLongSerialDenseMatrix.xml
%feature("docstring") Epetra_LongLongSerialDenseMatrix "

Epetra_LongLongSerialDenseMatrix: A class for constructing and using
general dense integer matrices.

The Epetra_LongLongSerialDenseMatrix class enables the construction
and use of integer-valued, general dense matrices.

The Epetra_LongLongSerialDenseMatrix class is intended to provide very
basic support for dense rectangular matrices.

Constructing Epetra_LongLongSerialDenseMatrix Objects

There are four Epetra_LongLongSerialDenseMatrix constructors. The
first constructs a zero-sized object which should be made to
appropriate length using the Shape() or Reshape() functions and then
filled with the [] or () operators. The second constructs an object
sized to the dimensions specified, which should be filled with the []
or () operators. The third is a constructor that accepts user data as
a 2D array, and the fourth is a copy constructor. The third
constructor has two data access modes (specified by the
Epetra_DataAccess argument): Copy mode - Allocates memory and makes a
copy of the user-provided data. In this case, the user data is not
needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Epetra_LongLongSerialDenseMatrix constructors
will throw an exception if an error occurrs. These exceptions will
alway be negative integer values as follows: -1 Invalid dimension
specified.

-2 Shape returned non-zero.

-3 Null pointer specified for user's data.

-99 Internal Epetra_LongLongSerialDenseMatrix error. Contact
developer.

Other Epetra_LongLongSerialDenseMatrix functions that do not return an
integer error code (such as operators () and [] ) will throw an
exception if an error occurrs. These exceptions will be integer values
as follows: -1 Invalid row specified.

-2 Invalid column specified.

-5 Invalid assignment (type mismatch).

-99 Internal Epetra_LongLongSerialDenseMatrix error. Contact
developer.

b Extracting Data from Epetra_LongLongSerialDenseMatrix Objects

Once a Epetra_LongLongSerialDenseMatrix is constructed, it is possible
to view the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.  Vector and Utility Functions

Once a Epetra_LongLongSerialDenseMatrix is constructed, several
mathematical functions can be applied to the object. Specifically:
Multiplication.

Norms.

C++ includes: Epetra_LongLongSerialDenseMatrix.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")
Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix "Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix()

Default constructor; defines a zero size object.

Epetra_LongLongSerialDenseMatrix objects defined by the default
constructor should be sized with the Shape() or Reshape functions.
Values should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix "Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix(int
NumRows, int NumCols)

Shaped constructor; defines a variable-sized object.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Epetra_SerialDenseMatrix objects defined by the shaped constructor are
already shaped to the dimensions given as a parameters. All values are
initialized to 0. Calling this constructor is equivalent to using the
default constructor, and then calling the Shape function on it. Values
should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix "Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix(Epetra_DataAccess
CV, long long *A, int LDA, int NumRows, int NumCols)

Set object values from two-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  A - Pointer to an array of integer numbers. The first vector
starts at A. The second vector starts at A+LDA, the third at A+2*LDA,
and so on.

In:  LDA - The \"Leading Dimension\", or stride between vectors in
memory.

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

See Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix "Epetra_LongLongSerialDenseMatrix::Epetra_LongLongSerialDenseMatrix(const
Epetra_LongLongSerialDenseMatrix &Source)

Epetra_LongLongSerialDenseMatrix copy constructor.

This matrix will take on the data access mode of the Source matrix. ";

%feature("docstring")
Epetra_LongLongSerialDenseMatrix::~Epetra_LongLongSerialDenseMatrix "Epetra_LongLongSerialDenseMatrix::~Epetra_LongLongSerialDenseMatrix()

Epetra_LongLongSerialDenseMatrix destructor. ";

/*  Shaping/sizing Methods  */

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::Shape "int
Epetra_LongLongSerialDenseMatrix::Shape(int NumRows, int NumCols)

Set dimensions of a Epetra_LongLongSerialDenseMatrix object; init
values to zero.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Allows user to define the dimensions of a
Epetra_LongLongSerialDenseMatrix at any point. This function can be
called at any point after construction. Any values that were
previously in this object are destroyed and the resized matrix starts
off with all zero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::Reshape "int
Epetra_LongLongSerialDenseMatrix::Reshape(int NumRows, int NumCols)

Reshape a Epetra_LongLongSerialDenseMatrix object.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Allows user to define the dimensions of a
Epetra_LongLongSerialDenseMatrix at any point. This function can be
called at any point after construction. Any values that were
previously in this object are copied into the new shape. If the new
shape is smaller than the original, the upper left portion of the
original matrix (the principal submatrix) is copied to the new matrix.

Integer error code, set to 0 if successful. ";

/*  Data Accessor methods  */

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::OneNorm "long long Epetra_LongLongSerialDenseMatrix::OneNorm()

Computes the 1-Norm of the this matrix.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::InfNorm "long long Epetra_LongLongSerialDenseMatrix::InfNorm()

Computes the Infinity-Norm of the this matrix. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::Random "int
Epetra_LongLongSerialDenseMatrix::Random()

Set matrix values to random numbers.

LongLongSerialDenseMatrix uses the random number generator provided by
Epetra_Util. The matrix values will be set to random values on the
interval (0, 2^31 - 1).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::M "int
Epetra_LongLongSerialDenseMatrix::M() const

Returns row dimension of system. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::N "int
Epetra_LongLongSerialDenseMatrix::N() const

Returns column dimension of system. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::A "const
long long* Epetra_LongLongSerialDenseMatrix::A() const

Returns const pointer to the this matrix. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::A "long
long* Epetra_LongLongSerialDenseMatrix::A()

Returns pointer to the this matrix. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::LDA "int
Epetra_LongLongSerialDenseMatrix::LDA() const

Returns the leading dimension of the this matrix. ";

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::CV "Epetra_DataAccess Epetra_LongLongSerialDenseMatrix::CV() const

Returns the data access mode of the this matrix. ";

/*  I/O methods  */

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::Print "void
Epetra_LongLongSerialDenseMatrix::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";

/*  Expert-only unsupported methods  */

%feature("docstring")  Epetra_LongLongSerialDenseMatrix::MakeViewOf "int Epetra_LongLongSerialDenseMatrix::MakeViewOf(const
Epetra_LongLongSerialDenseMatrix &Source)

Reset an existing LongLongSerialDenseMatrix to point to another
Matrix.

Allows an existing LongLongSerialDenseMatrix to become a View of
another matrix's data, regardless of the DataAccess mode of the Source
matrix. It is assumed that the Source matrix is an independent matrix,
and no checking is done to verify this.

This is used by Epetra_CrsGraph in the OptimizeStorage method. It is
used so that an existing (Copy) matrix can be converted to a View.
This frees up memory that CrsGraph no longer needs.

Parameters:
-----------

Source:  The LongLongSerialDenseMatrix this will become a view of.

Integer error code, set to 0 if successful, and set to -1 if a type
mismatch occured.

WARNING:  This method is extremely dangerous and should only be used
by experts. ";


// File: classEpetra__LongLongSerialDenseVector.xml
%feature("docstring") Epetra_LongLongSerialDenseVector "

Epetra_LongLongSerialDenseVector: A class for constructing and using
dense vectors.

The Epetra_LongLongSerialDenseVector class enables the construction
and use of integer-valued, dense vectors. It derives from the
Epetra_LongLongSerialDenseMatrix class.

The Epetra_LongLongSerialDenseVector class is intended to provide
convenient vector notation but derives all signficant functionality
from Epetra_LongLongSerialDenseMatrix.

Constructing Epetra_LongLongSerialDenseVector Objects

There are three Epetra_LongLongSerialDenseVector constructors. The
first constructs a zero-length object which should be made to
appropriate length using the Size() or Resize() functions and then
filled with the [] or () operators. The second constructs an object
sized to the dimension specified, which should be filled with the []
or () operators. The third is a constructor that accepts user data as
a 1D array, and the fourth is a copy constructor. The third
constructor has two data access modes (specified by the
Epetra_DataAccess argument): Copy mode - Allocates memory and makes a
copy of the user-provided data. In this case, the user data is not
needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Extracting Data from
Epetra_LongLongSerialDenseVector Objects

Once a Epetra_LongLongSerialDenseVector is constructed, it is possible
to view the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.

C++ includes: Epetra_LongLongSerialDenseVector.h ";

/*  I/O methods  */

%feature("docstring")  Epetra_LongLongSerialDenseVector::Print "void
Epetra_LongLongSerialDenseVector::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";

/*  Expert-only unsupported methods  */

%feature("docstring")  Epetra_LongLongSerialDenseVector::MakeViewOf "int Epetra_LongLongSerialDenseVector::MakeViewOf(const
Epetra_LongLongSerialDenseVector &Source)

Reset an existing LongLongSerialDenseVector to point to another
Vector.

Allows an existing LongLongSerialDenseVector to become a View of
another vector's data, regardless of the DataAccess mode of the Source
vector. It is assumed that the Source vector is an independent vector,
and no checking is done to verify this.

This is used by Epetra_CrsGraph in the OptimizeStorage method. It is
used so that an existing (Copy) vector can be converted to a View.
This frees up memory that CrsGraph no longer needs.

Parameters:
-----------

Source:  The LongLongSerialDenseVector this will become a view of.

Integer error code, set to 0 if successful.

WARNING:  This method is extremely dangerous and should only be used
by experts. ";

%feature("docstring")
Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector "Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector()

Default constructor; defines a zero size object.

Epetra_LongLongSerialDenseVector objects defined by the default
constructor should be sized with the Size() or Resize functions.
Values should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector "Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector(int
Length_in)

Sized constructor; defines a variable-sized object.

Parameters:
-----------

In:  Length - Length of vector.

Epetra_LongLongSerialDenseVector objects defined by the sized
constructor are already sized to the dimension given as a parameter.
All values are initialized to 0. Calling this constructor is
equivalent to using the default constructor, and then calling the Size
function on it. Values should be defined by using the [] or ()
operators. ";

%feature("docstring")
Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector "Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector(Epetra_DataAccess
CV_in, long long *Values_in, int Length_in)

Set object values from one-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Values - Pointer to an array of integer numbers containing the
values.

In:  Length - Length of vector.

See Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector "Epetra_LongLongSerialDenseVector::Epetra_LongLongSerialDenseVector(const
Epetra_LongLongSerialDenseVector &Source)

Epetra_LongLongSerialDenseVector copy constructor. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::Size "int
Epetra_LongLongSerialDenseVector::Size(int Length_in)

Set length of a Epetra_LongLongSerialDenseVector object; init values
to zero.

Parameters:
-----------

In:  Length - Length of vector object.

Allows user to define the dimension of a
Epetra_LongLongSerialDenseVector. This function can be called at any
point after construction. Any values that were previously in this
object are destroyed and the resized vector starts off with all zero
values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::Resize "int
Epetra_LongLongSerialDenseVector::Resize(int Length_in)

Resize a Epetra_LongLongSerialDenseVector object.

Parameters:
-----------

In:  Length - Length of vector object.

Allows user to define the dimension of a
Epetra_LongLongSerialDenseVector. This function can be called at any
point after construction. Any values that were previously in this
object are copied into the new size. If the new shape is smaller than
the original, the first Length values are copied to the new vector.

Integer error code, set to 0 if successful. ";

%feature("docstring")
Epetra_LongLongSerialDenseVector::~Epetra_LongLongSerialDenseVector "Epetra_LongLongSerialDenseVector::~Epetra_LongLongSerialDenseVector()

Epetra_LongLongSerialDenseVector destructor. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::Random "int
Epetra_LongLongSerialDenseVector::Random()

Set vector values to random numbers.

LongLongSerialDenseVector uses the random number generator provided by
Epetra_Util. The vector values will be set to random values on the
interval (0, 2^31 - 1).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::Length "int
Epetra_LongLongSerialDenseVector::Length() const

Returns length of vector. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::Values "long
long* Epetra_LongLongSerialDenseVector::Values()

Returns pointer to the values in vector. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::Values "const long long* Epetra_LongLongSerialDenseVector::Values() const

Returns const pointer to the values in vector. ";

%feature("docstring")  Epetra_LongLongSerialDenseVector::CV "Epetra_DataAccess Epetra_LongLongSerialDenseVector::CV() const

Returns the data access mode of the this vector. ";


// File: classEpetra__LongLongVector.xml
%feature("docstring") Epetra_LongLongVector "

Epetra_LongLongVector: A class for constructing and using dense
integer vectors on a parallel computer.

The Epetra_LongLongVector class enables the construction and use of
integer dense vectors in a distributed memory environment. The
distribution of the dense vector is determined in part by a
Epetra_Comm object and a Epetra_Map (or Epetra_LocalMap or
Epetra_BlockMap).

Distributed Global vs. Replicated Local Distributed Global Vectors -
In most instances, a multi-vector will be partitioned across multiple
memory images associated with multiple processors. In this case, there
is a unique copy of each element and elements are spread across all
processors specified by the Epetra_Comm communicator.

Replicated Local Vectors - Some algorithms use vectors that are too
small to be distributed across all processors. Replicated local
vectors handle these types of situation.

Constructing Epetra_LongLongVectors

There are four Epetra_LongLongVector constructors. The first is a
basic constructor that allocates space and sets all values to zero,
the second is a copy constructor. The third and fourth constructors
work with user data. These constructors have two data access modes:
Copy mode - Allocates memory and makes a copy of the user-provided
data. In this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the vector.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  All Epetra_LongLongVector constructors require a
map argument that describes the layout of elements on the parallel
machine. Specifically, map is a Epetra_Map, Epetra_LocalMap or
Epetra_BlockMap object describing the desired memory layout for the
vector.

There are four different Epetra_LongLongVector constructors: Basic -
All values are zero.

Copy - Copy an existing vector.

Copy from or make view of user int array.

Extracting Data from Epetra_LongLongVectors

Once a Epetra_LongLongVector is constructed, it is possible to extract
a copy of the values or create a view of them.

WARNING:  ExtractView functions are extremely dangerous from a data
hiding perspective. For both ExtractView fuctions, there is a
corresponding ExtractCopy function. We strongly encourage users to
develop code using ExtractCopy functions first and only use the
ExtractView functions in a secondary optimization phase.  There are
two Extract functions: ExtractCopy - Copy values into a user-provided
array.

ExtractView - Set user-provided array to point to
Epetra_LongLongVector data.

WARNING:  A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is
required for all Epetra_LongLongVector constructors.

C++ includes: Epetra_LongLongVector.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_LongLongVector::Epetra_LongLongVector "Epetra_LongLongVector::Epetra_LongLongVector(const Epetra_BlockMap
&Map, bool zeroOut=true)

Basic Epetra_LongLongVector constuctor.

Creates a Epetra_LongLongVector object and, by default, fills with
zero values.

Parameters:
-----------

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

WARNING:  Note that, because Epetra_LocalMap derives from Epetra_Map
and Epetra_Map derives from Epetra_BlockMap, this constructor works
for all three types of Epetra map classes.

Parameters:
-----------

In:  zeroOut - If true then the allocated memory will be zeroed out
initialy. If false then this memory will not be touched which can be
significantly faster.

Pointer to a Epetra_LongLongVector. ";

%feature("docstring")  Epetra_LongLongVector::Epetra_LongLongVector "Epetra_LongLongVector::Epetra_LongLongVector(const
Epetra_LongLongVector &Source)

Epetra_LongLongVector copy constructor. ";

%feature("docstring")  Epetra_LongLongVector::Epetra_LongLongVector "Epetra_LongLongVector::Epetra_LongLongVector(Epetra_DataAccess CV,
const Epetra_BlockMap &Map, long long *V)

Set vector values from user array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  V - Pointer to an array of long long numbers..

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_LongLongVector::~Epetra_LongLongVector "Epetra_LongLongVector::~Epetra_LongLongVector()

Epetra_LongLongVector destructor. ";

/*  Post-construction modification methods  */

%feature("docstring")  Epetra_LongLongVector::PutValue "int
Epetra_LongLongVector::PutValue(long long Value)

Set all elements of the vector to Value. ";

/*  Extraction methods  */

%feature("docstring")  Epetra_LongLongVector::ExtractCopy "int
Epetra_LongLongVector::ExtractCopy(long long *V) const

Put vector values into user-provided array.

Parameters:
-----------

Out:  V - Pointer to memory space that will contain the vector values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_LongLongVector::ExtractView "int
Epetra_LongLongVector::ExtractView(long long **V) const

Set user-provided address of V.

Parameters:
-----------

Out:  V - Address of a pointer to that will be set to point to the
values of the vector.

Integer error code, set to 0 if successful. ";

/*  Mathematical methods  */

%feature("docstring")  Epetra_LongLongVector::MaxValue "long long
Epetra_LongLongVector::MaxValue()

Find maximum value.

Maximum value across all processors. ";

%feature("docstring")  Epetra_LongLongVector::MinValue "long long
Epetra_LongLongVector::MinValue()

Find minimum value.

Minimum value across all processors. ";

/*  Overloaded operators  */

/*  Attribute access functions  */

%feature("docstring")  Epetra_LongLongVector::Values "long long*
Epetra_LongLongVector::Values() const

Returns a pointer to an array containing the values of this vector. ";

%feature("docstring")  Epetra_LongLongVector::MyLength "int
Epetra_LongLongVector::MyLength() const

Returns the local vector length on the calling processor of vectors in
the multi-vector. ";

%feature("docstring")  Epetra_LongLongVector::GlobalLength "long long
Epetra_LongLongVector::GlobalLength() const

Returns the global vector length of vectors in the multi-vector. ";

/*  I/O methods  */

%feature("docstring")  Epetra_LongLongVector::Print "void
Epetra_LongLongVector::Print(ostream &os) const

Print method. ";


// File: classEpetra__Map.xml
%feature("docstring") Epetra_Map "

Epetra_Map: A class for partitioning vectors and matrices.

It is often the case that multiple matrix and vector objects have an
identical distribution of elements on a parallel machine. The
Epetra_Map class keep information that describes this distribution for
matrices and vectors.

Epetra_Map allows the storage and retrieval of the following
information. Depending on the constructor that is used, some of the
information is defined by the user and some is determined by the
constructor. Once a Epetra_Map is constructed any of the following
attributes can be obtained by calling a query function that has the
name as the attribute, e.g. to get the value of NumGlobalElements, you
can call a function NumGlobalElements(). For attributes that are
lists, the query functions return the list values in a user allocated
array.

NumGlobalElements - The total number of elements across all
processors. If this parameter and NumMyElements are both passed into
the constructor, one of the three cases will apply: If
NumGlobalElements = NumMyElements (and not equal to zero) the map is
defined to be a local replicated map. In this case, objects
constructed using this map will be identically replicated across all
processors in the communicator.

If NumGlobalElements = -1 and NumMyElements is passed in then
NumGlobalElements will be computed as the sum of NumMyElements across
all processors.

If neither of the above is true, NumGlobalElements will be checked
against the sum of NumMyElements across all processors. An error is
issued if the comparison is not equal.

NumMyElements - The number of elements owned by the calling processor.

MyGlobalElements - A list of length NumMyElements that contains the
global element IDs of the elements owned by the calling processor.

IndexBase - The base integer value for indexed array references.
Typically this is 0 for C/C++ and 1 for Fortran, but it can be set to
any integer value.

Comm - The Epetra_Comm communicator. This communicator can in turn be
queried for processor rank and size information.

In addition to the information above that is passed in to or created
by the Epetra_Map constructor, the following attributes are computed
and available via query to the user using the same scheme as above,
e.g., use NumGlobalPoints() to get the value of NumGlobalPoints.

NumGlobalPoints - The total number of points across all processors.

NumMyPoints - The number of points on the calling processor.

MinAllGID - The minimum global index value across all processors.

MaxAllGID - The maximum global index value across all processors.

MinMyGID - The minimum global index value on the calling processor.

MaxMyGID - The maximum global index value on the calling processor.

MinLID - The minimum local index value on the calling processor.

MaxLID - The maximum local index value on the calling processor.

The following functions allow boolean tests for certain properties.

LinearMap() - Returns true if the elements are distributed linear
across processors, i.e., processor 0 gets the first n/p elements,
processor 1 gets the next n/p elements, etc. where n is the number of
elements and p is the number of processors.

DistributedGlobal() - Returns true if the element space of the map
spans more than one processor. This will be true in most cases, but
will be false in serial cases and for objects that are created via the
derived Epetra_LocalMap class.

WARNING:  An Epetra_Comm object is required for all Epetra_Map
constructors.

In the current implementation, Epetra_BlockMap is the base class for
Epetra_Map.

C++ includes: Epetra_Map.h ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(int NumGlobalElements, int IndexBase, const
Epetra_Comm &Comm)

Epetra_Map constructor for a Epetra-defined uniform linear
distribution of elements.

Creates a map that distributes NumGlobalElements elements evenly
across all processors in the Epetra_Comm communicator. If
NumGlobalElements does not divide exactly into the number of
processors, the first processors in the communicator get one extra
element until the remainder is gone.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_Map object. ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(unsigned int NumGlobalElements, int IndexBase,
const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(long long NumGlobalElements, int IndexBase,
const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(unsigned long long NumGlobalElements, int
IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements, int
IndexBase, const Epetra_Comm &Comm)

Epetra_Map constructor for a user-defined linear distribution of
elements.

Creates a map that puts NumMyElements on the calling processor. If
NumGlobalElements=-1, the number of global elements will be the
computed sum of NumMyElements across all processors in the Epetra_Comm
communicator.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute. Must be
either -1 or equal to the computed sum of NumMyElements across all
processors in the Epetra_Comm communicator.

In:  NumMyElements - Number of elements owned by the calling
processor.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_Map object. ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(unsigned int NumGlobalElements, int
NumMyElements, int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(long long NumGlobalElements, int NumMyElements,
int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(unsigned long long NumGlobalElements, int
NumMyElements, int IndexBase, const Epetra_Comm &Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(int NumGlobalElements, int NumMyElements, const
int *MyGlobalElements, int IndexBase, const Epetra_Comm &Comm)

Epetra_Map constructor for a user-defined arbitrary distribution of
elements.

Creates a map that puts NumMyElements on the calling processor. The
indices of the elements are determined from the list MyGlobalElements.
If NumGlobalElements=-1, the number of global elements will be the
computed sum of NumMyElements across all processors in the Epetra_Comm
communicator.

Parameters:
-----------

In:  NumGlobalElements - Number of elements to distribute. Must be
either -1 or equal to the computed sum of NumMyElements across all
processors in the Epetra_Comm communicator.

In:  NumMyElements - Number of elements owned by the calling
processor.

In:  MyGlobalElements - Integer array of length NumMyElements. The ith
entry contains the global index value of the ith element on this
processor. Index values are not required to be contiguous on a
processor, or to be within the range of 0 to NumGlobalElements. As
long as the index values are consistently defined and used, any set of
NumGlobalElements distinct integer values is acceptable.

In:  IndexBase - Minimum index value used for arrays that use this
map. Typically 0 for C/C++ and 1 for Fortran.

In:  Comm - Epetra_Comm communicator containing information on the
number of processors.

Pointer to a Epetra_Map object. ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(long long NumGlobalElements, int NumMyElements,
const long long *MyGlobalElements, int IndexBase, const Epetra_Comm
&Comm) ";

%feature("docstring")  Epetra_Map::Epetra_Map "Epetra_Map::Epetra_Map(const Epetra_Map &map)

Epetra_Map copy constructor. ";

%feature("docstring")  Epetra_Map::~Epetra_Map "Epetra_Map::~Epetra_Map(void)

Epetra_Map destructor. ";


// File: classEpetra__MapColoring.xml
%feature("docstring") Epetra_MapColoring "

Epetra_MapColoring: A class for coloring Epetra_Map and
Epetra_BlockMap objects.

This class allows the user to associate an integer value, i.e., a
color, to each element of an existing Epetra_Map or Epetra_BlockMap
object. Colors may be assigned at construction, or via set methods.
Any elements that are not explicitly assigned a color are assigned the
color 0 (integer zero).

This class has the following features:

A color (arbitrary integer label) can be associated locally with each
element of a map. Color assignment can be done all-at-once via the
constructor, or

via operator[] (using LIDs) one-at-a-time

operator() (using GIDs) one-at-a-time

or some combination of the above.

Any element that is not explicitly colored takes on the default color.
The default color is implicitly zero, unless specified differently at
the time of construction.

Color information may be accessed in the following ways: By local
element ID (LID) - Returns the color of a specified LID, where the LID
is associated with the Epetra_Map or BlockMap that was passed in to
the Epetra_MapColoring constructor.

By global element ID (GID) - Returns the color of the specified GID.
There two methods for accessing GIDs, one assumes the request is for
GIDs owned by the calling processor, the second allows arbitrary
requested for GIDs, as long as the GID is defined on some processor
for the Epetra_Map or Epetra_BlockMap.

By color groups - Elements are grouped by color so that all elements
of a given color can be accessed.

Epetra_Map/Epetra_BlockMap pointers for a specified color - This
facilitates use of coloring with Epetra distributed objects that are
distributed via the map that was colored. For example, if users want
to work with all rows of a matrix that have a certain color, they can
create a map for that color and use it to access only those rows.

The Epetra_MapColoring class implements the Epetra_DistObject
interface. Therefore, a map coloring can be computed for a map with a
given distribution and then redistributed across the parallel machine.
For example, it would be possible to compute a map coloring on a
single processor (perhaps because the algorithm for computing the
color assignment is too difficult to implement in parallel or because
it is cheap to run and not worth parallelizing), and then re-
distribute the coloring using an Epetra_Export or Epetra_Import
object.

C++ includes: Epetra_MapColoring.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_MapColoring::Epetra_MapColoring "Epetra_MapColoring::Epetra_MapColoring(const Epetra_BlockMap &Map,
const int DefaultColor=0)

Epetra_MapColoring basic constructor.

Parameters:
-----------

In:  Map - An Epetra_Map or Epetra_BlockMap (Note: Epetra_BlockMap is
a base class of Epetra_Map, so either can be passed in to this
constructor.

In:  DefaultColor - The integer value to use as the default color for
this map. This constructor will initially define the color of all map
elements to the default color.

Pointer to a Epetra_MapColoring object. ";

%feature("docstring")  Epetra_MapColoring::Epetra_MapColoring "Epetra_MapColoring::Epetra_MapColoring(const Epetra_BlockMap &Map, int
*ElementColors, const int DefaultColor=0)

Epetra_MapColoring constructor.

Parameters:
-----------

In:  Map - An Epetra_Map or Epetra_BlockMap (Note: Epetra_BlockMap is
a base class of Epetra_Map, so either can be passed in to this
constructor.

In:  ElementColors - Array of dimension Map.NumMyElements() containing
the list of colors that should be assigned the map elements on this
processor. If this argument is set to 0 (zero), all elements will
initially be assigned color 0 (zero). Element colors can be modified
by using methods described below.

In:  DefaultColor - The color that will be assigned by default when no
other value is specified. This value has no meaning for this
constructor, but is used by certain methods now and in the future.

Pointer to a Epetra_MapColoring object. ";

%feature("docstring")  Epetra_MapColoring::Epetra_MapColoring "Epetra_MapColoring::Epetra_MapColoring(const Epetra_MapColoring
&Source)

Epetra_MapColoring copy constructor. ";

%feature("docstring")  Epetra_MapColoring::~Epetra_MapColoring "Epetra_MapColoring::~Epetra_MapColoring()

Epetra_MapColoring destructor. ";

/*  Set Color methods  */

/*  Local/Global color accessor methods  */

/*  Color Information Access Methods  */

%feature("docstring")  Epetra_MapColoring::NumColors "int
Epetra_MapColoring::NumColors() const

Returns number of colors on the calling processor. ";

%feature("docstring")  Epetra_MapColoring::MaxNumColors "int
Epetra_MapColoring::MaxNumColors() const

Returns maximum over all processors of the number of colors. ";

%feature("docstring")  Epetra_MapColoring::ListOfColors "int*
Epetra_MapColoring::ListOfColors() const

Array of length NumColors() containing List of color values used in
this coloring.

Color values can be arbitrary integer values. As a result, a user of a
previously constructed MapColoring object may need to know exactly
which color values are present. This array contains that information
as a sorted list of integer values. ";

%feature("docstring")  Epetra_MapColoring::DefaultColor "int
Epetra_MapColoring::DefaultColor() const

Returns default color. ";

%feature("docstring")  Epetra_MapColoring::NumElementsWithColor "int
Epetra_MapColoring::NumElementsWithColor(int Color) const

Returns number of map elements on calling processor having specified
Color. ";

%feature("docstring")  Epetra_MapColoring::ColorLIDList "int *
Epetra_MapColoring::ColorLIDList(int Color) const

Returns pointer to array of Map LIDs associated with the specified
color.

Returns a pointer to a list of Map LIDs associated with the specified
color. This is a purely local list with no information about other
processors. If there are no LIDs associated with the specified color,
the pointer is set to zero. ";

%feature("docstring")  Epetra_MapColoring::ElementColors "int*
Epetra_MapColoring::ElementColors() const

Returns pointer to array of the colors associated with the LIDs on the
calling processor.

Returns a pointer to the list of colors associated with the elements
on this processor such that ElementColor[LID] is the color assigned to
that LID. ";

/*  Epetra_Map and Epetra_BlockMap generators  */

%feature("docstring")  Epetra_MapColoring::GenerateMap "Epetra_Map *
Epetra_MapColoring::GenerateMap(int Color) const

Generates an Epetra_Map of the GIDs associated with the specified
color.

This method will create an Epetra_Map such that on each processor the
GIDs associated with the specified color will be part of the map on
that processor. Note that this method always generates an Epetra_Map,
not an Epetra_BlockMap, even if the map associated with this map
coloring is a block map. Once the map is generated, the user is
responsible for deleting it. ";

%feature("docstring")  Epetra_MapColoring::GenerateBlockMap "Epetra_BlockMap * Epetra_MapColoring::GenerateBlockMap(int Color)
const

Generates an Epetra_BlockMap of the GIDs associated with the specified
color.

This method will create an Epetra_BlockMap such that on each processor
the GIDs associated with the specified color will be part of the map
on that processor. Note that this method will generate an
Epetra_BlockMap such that each element as the same element size as the
corresponding element of map associated with the map coloring. Once
the map is generated, the user is responsible for deleting it. ";

/*  I/O methods  */

%feature("docstring")  Epetra_MapColoring::Print "void
Epetra_MapColoring::Print(ostream &os) const

Print method. ";


// File: classEpetra__MpiComm.xml
%feature("docstring") Epetra_MpiComm "

Epetra_MpiComm: The Epetra MPI Communication Class.

The Epetra_MpiComm class is an implementation of Epetra_Comm that
encapsulates the general information and services needed for other
Epetra classes to run on a parallel computer using MPI.

C++ includes: Epetra_MpiComm.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")  Epetra_MpiComm::Epetra_MpiComm "Epetra_MpiComm::Epetra_MpiComm(MPI_Comm comm)

Epetra_MpiComm MPI Constructor.

Creates a Epetra_MpiComm instance for use with MPI. If no specialized
MPI communicator is needed, this constuctor can be called with the
argument MPI_COMM_WORLD. ";

%feature("docstring")  Epetra_MpiComm::Epetra_MpiComm "Epetra_MpiComm::Epetra_MpiComm(const Epetra_MpiComm &Comm)

Epetra_MpiComm Copy Constructor.

Makes an exact copy of an existing Epetra_MpiComm instance. ";

%feature("docstring")  Epetra_MpiComm::Clone "Epetra_Comm*
Epetra_MpiComm::Clone() const

Clone method. ";

%feature("docstring")  Epetra_MpiComm::~Epetra_MpiComm "Epetra_MpiComm::~Epetra_MpiComm()

Epetra_MpiComm Destructor.

Completely deletes a Epetra_MpiComm object. WARNING:  Note: All
objects that depend on a Epetra_MpiComm instance should be destroyed
prior to calling this function. ";

/*  Barrier Methods  */

%feature("docstring")  Epetra_MpiComm::Barrier "void
Epetra_MpiComm::Barrier() const

Epetra_MpiComm Barrier function.

Causes each processor in the communicator to wait until all processors
have arrived. ";

/*  Broadcast Methods  */

%feature("docstring")  Epetra_MpiComm::Broadcast "int
Epetra_MpiComm::Broadcast(double *MyVals, int Count, int Root) const

Epetra_MpiComm Broadcast function.

Takes list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_MpiComm::Broadcast "int
Epetra_MpiComm::Broadcast(int *MyVals, int Count, int Root) const

Epetra_MpiComm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_MpiComm::Broadcast "int
Epetra_MpiComm::Broadcast(long *MyVals, int Count, int Root) const

Epetra_MpiComm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_MpiComm::Broadcast "int
Epetra_MpiComm::Broadcast(char *MyVals, int Count, int Root) const

Epetra_MpiComm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

/*  Gather Methods  */

%feature("docstring")  Epetra_MpiComm::GatherAll "int
Epetra_MpiComm::GatherAll(double *MyVals, double *AllVals, int Count)
const

Epetra_MpiComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_MpiComm::GatherAll "int
Epetra_MpiComm::GatherAll(int *MyVals, int *AllVals, int Count) const

Epetra_MpiComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_MpiComm::GatherAll "int
Epetra_MpiComm::GatherAll(long *MyVals, long *AllVals, int Count)
const

Epetra_MpiComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_MpiComm::GatherAll "int
Epetra_MpiComm::GatherAll(long long *MyVals, long long *AllVals, int
Count) const

Epetra_MpiComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

/*  Sum Methods  */

%feature("docstring")  Epetra_MpiComm::SumAll "int
Epetra_MpiComm::SumAll(double *PartialSums, double *GlobalSums, int
Count) const

Epetra_MpiComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::SumAll "int
Epetra_MpiComm::SumAll(int *PartialSums, int *GlobalSums, int Count)
const

Epetra_MpiComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::SumAll "int
Epetra_MpiComm::SumAll(long *PartialSums, long *GlobalSums, int Count)
const

Epetra_MpiComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::SumAll "int
Epetra_MpiComm::SumAll(long long *PartialSums, long long *GlobalSums,
int Count) const

Epetra_MpiComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Max/Min Methods  */

%feature("docstring")  Epetra_MpiComm::MaxAll "int
Epetra_MpiComm::MaxAll(double *PartialMaxs, double *GlobalMaxs, int
Count) const

Epetra_MpiComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MaxAll "int
Epetra_MpiComm::MaxAll(int *PartialMaxs, int *GlobalMaxs, int Count)
const

Epetra_MpiComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MaxAll "int
Epetra_MpiComm::MaxAll(long *PartialMaxs, long *GlobalMaxs, int Count)
const

Epetra_MpiComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MaxAll "int
Epetra_MpiComm::MaxAll(long long *PartialMaxs, long long *GlobalMaxs,
int Count) const

Epetra_MpiComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MinAll "int
Epetra_MpiComm::MinAll(double *PartialMins, double *GlobalMins, int
Count) const

Epetra_MpiComm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MinAll "int
Epetra_MpiComm::MinAll(int *PartialMins, int *GlobalMins, int Count)
const

Epetra_MpiComm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MinAll "int
Epetra_MpiComm::MinAll(long *PartialMins, long *GlobalMins, int Count)
const

Epetra_MpiComm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::MinAll "int
Epetra_MpiComm::MinAll(long long *PartialMins, long long *GlobalMins,
int Count) const

Epetra_MpiComm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Parallel Prefix Methods  */

%feature("docstring")  Epetra_MpiComm::ScanSum "int
Epetra_MpiComm::ScanSum(double *MyVals, double *ScanSums, int Count)
const

Epetra_MpiComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::ScanSum "int
Epetra_MpiComm::ScanSum(int *MyVals, int *ScanSums, int Count) const

Epetra_MpiComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::ScanSum "int
Epetra_MpiComm::ScanSum(long *MyVals, long *ScanSums, int Count) const

Epetra_MpiComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiComm::ScanSum "int
Epetra_MpiComm::ScanSum(long long *MyVals, long long *ScanSums, int
Count) const

Epetra_MpiComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

/*  Attribute Accessor Methods  */

%feature("docstring")  Epetra_MpiComm::Comm "MPI_Comm
Epetra_MpiComm::Comm() const

Extract MPI Communicator from a Epetra_MpiComm object. ";

%feature("docstring")  Epetra_MpiComm::MyPID "int
Epetra_MpiComm::MyPID() const

Return my process ID.

In MPI mode returns the rank of the calling process. In serial mode
returns 0. ";

%feature("docstring")  Epetra_MpiComm::NumProc "int
Epetra_MpiComm::NumProc() const

Returns total number of processes.

In MPI mode returns the size of the MPI communicator. In serial mode
returns 1. ";

/*  Gather/Scatter and Directory Constructors  */

%feature("docstring")  Epetra_MpiComm::CreateDistributor "Epetra_Distributor * Epetra_MpiComm::CreateDistributor() const

Create a distributor object. ";

%feature("docstring")  Epetra_MpiComm::CreateDirectory "Epetra_Directory * Epetra_MpiComm::CreateDirectory(const
Epetra_BlockMap &Map) const

Create a directory object for the given Epetra_BlockMap. ";

/*  MPI-specific Methods  */

%feature("docstring")  Epetra_MpiComm::GetMpiTag "int
Epetra_MpiComm::GetMpiTag() const

Acquire an MPI tag from the Epetra range of 24050-24099, increment
tag. ";

%feature("docstring")  Epetra_MpiComm::GetMpiComm "MPI_Comm
Epetra_MpiComm::GetMpiComm() const

Get the MPI Communicator (identical to Comm() method; used when we
know we are MPI. ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_MpiComm::Print "void
Epetra_MpiComm::Print(ostream &os) const

Print method that implements Epetra_Object virtual Print method. ";

%feature("docstring")  Epetra_MpiComm::PrintInfo "void
Epetra_MpiComm::PrintInfo(ostream &os) const

Print method that implements Epetra_Comm virtual PrintInfo method. ";

/*  Expert Users and Developers Only  */

%feature("docstring")  Epetra_MpiComm::ReferenceCount "int
Epetra_MpiComm::ReferenceCount() const

Returns the reference count of MpiCommData.

(Intended for testing purposes.) ";

%feature("docstring")  Epetra_MpiComm::DataPtr "const
Epetra_MpiCommData* Epetra_MpiComm::DataPtr() const

Returns a pointer to the MpiCommData instance this MpiComm uses.

(Intended for developer use only for testing purposes.) ";


// File: classEpetra__MpiCommData.xml
%feature("docstring") Epetra_MpiCommData "

Epetra_MpiCommData: The Epetra Mpi Communication Data Class.

The Epetra_MpiCommData class is an implementation detail of
Epetra_MpiComm. It is reference-counted, and can be shared by multiple
Epetra_MpiComm instances. It derives from Epetra_Data, and inherits
reference-counting from it.

C++ includes: Epetra_MpiCommData.h ";

/*  Constructor/Destructor Methods  */


// File: classEpetra__MpiDistributor.xml
%feature("docstring") Epetra_MpiDistributor "

MPI implementation of Epetra_Distributor.

This class is an MPI implementation of  Epetra_Distributor. It
encapsulates the general information and services needed for other
Epetra classes to perform gather/scatter operations on a parallel
computer. An Epetra_MpiDistributor instance is actually produced by
calling a method in the Epetra_MpiComm class.

C++ includes: Epetra_MpiDistributor.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_MpiDistributor::Epetra_MpiDistributor "Epetra_MpiDistributor::Epetra_MpiDistributor(const Epetra_MpiComm
&Comm)

Default constructor. ";

%feature("docstring")  Epetra_MpiDistributor::Epetra_MpiDistributor "Epetra_MpiDistributor::Epetra_MpiDistributor(const
Epetra_MpiDistributor &Distributor)

Copy constructor. ";

%feature("docstring")  Epetra_MpiDistributor::Clone "Epetra_Distributor* Epetra_MpiDistributor::Clone()

Clone method. ";

%feature("docstring")  Epetra_MpiDistributor::~Epetra_MpiDistributor "Epetra_MpiDistributor::~Epetra_MpiDistributor()

Destructor (declared virtual for memory safety). ";

/*  Gather/Scatter Constructors  */

%feature("docstring")  Epetra_MpiDistributor::CreateFromSends "int
Epetra_MpiDistributor::CreateFromSends(const int &NumExportIDs, const
int *ExportPIDs, bool Deterministic, int &NumRemoteIDs)

Create a communication plan from send list.

Given a list of process IDs to which to send the given number of data
IDs, construct a communication plan for efficiently scattering data to
these processes.

The number of data IDs being sent to me.

Parameters:
-----------

NumExportIDs:  [in] Number of data IDs that need to be sent from the
calling process.

ExportPIDs:  [in] List of process IDs that will get the exported data
IDs.

Deterministic:  [in] Currently has no effect.

NumRemoteIDs:  [out] Number of data IDs the calling process will be
receiving. ";

%feature("docstring")  Epetra_MpiDistributor::CreateFromRecvs "int
Epetra_MpiDistributor::CreateFromRecvs(const int &NumRemoteIDs, const
int *RemoteGIDs, const int *RemotePIDs, bool Deterministic, int
&NumExportIDs, int *&ExportGIDs, int *&ExportPIDs)

Create a communication plan from receive list.

Given a list of remote data IDs and corresponding process IDs from
which to receive data, construct a communication plan for efficiently
scattering data to these processes.

The number and list of data IDs being sent by me.

Parameters:
-----------

NumRemoteIDs:  [in] Number of data IDs the calling process will be
receiving.

RemoteGIDs:  [in] List of data IDs that the calling process wants to
receive.

RemotePIDs:  [in] List of IDs of the processes that will send the
remote data IDs to the calling process.

Deterministic:  [in] Currently has no effect.

NumExportIDs:  [out] Number of data IDs that need to be sent from the
calling process.

ExportGIDs:  [out] List of data IDs that the calling process will send
out.

ExportPIDs:  [out] List of IDs of the processes that will receive the
data IDs sent by the calling process.

This method allocates the output arrays using new. The caller is
responsible for deallocating them after use. ";

%feature("docstring")  Epetra_MpiDistributor::CreateFromRecvs "int
Epetra_MpiDistributor::CreateFromRecvs(const int &NumRemoteIDs, const
long long *RemoteGIDs, const int *RemotePIDs, bool Deterministic, int
&NumExportIDs, long long *&ExportGIDs, int *&ExportPIDs) ";

/*  Execute Gather/Scatter Operations  */

%feature("docstring")  Epetra_MpiDistributor::Do "int
Epetra_MpiDistributor::Do(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)

Execute plan on buffer of export objects in a single step. ";

%feature("docstring")  Epetra_MpiDistributor::DoReverse "int
Epetra_MpiDistributor::DoReverse(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)

Execute reverse of plan on buffer of export objects in a single step.
";

%feature("docstring")  Epetra_MpiDistributor::DoPosts "int
Epetra_MpiDistributor::DoPosts(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)

Post buffer of export objects (can do other local work before
executing Waits) ";

%feature("docstring")  Epetra_MpiDistributor::DoWaits "int
Epetra_MpiDistributor::DoWaits()

Wait on a set of posts. ";

%feature("docstring")  Epetra_MpiDistributor::DoReversePosts "int
Epetra_MpiDistributor::DoReversePosts(char *export_objs, int obj_size,
int &len_import_objs, char *&import_objs)

Do reverse post of buffer of export objects (can do other local work
before executing Waits) ";

%feature("docstring")  Epetra_MpiDistributor::DoReverseWaits "int
Epetra_MpiDistributor::DoReverseWaits()

Wait on a reverse set of posts. ";

/*  Execute Gather/Scatter Operations (Non-constant size objects)  */

%feature("docstring")  Epetra_MpiDistributor::Do "int
Epetra_MpiDistributor::Do(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Execute plan on buffer of export objects in a single step (object size
may vary) ";

%feature("docstring")  Epetra_MpiDistributor::DoReverse "int
Epetra_MpiDistributor::DoReverse(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Execute reverse of plan on buffer of export objects in a single step
(object size may vary) ";

%feature("docstring")  Epetra_MpiDistributor::DoPosts "int
Epetra_MpiDistributor::DoPosts(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Post buffer of export objects (can do other local work before
executing Waits) ";

%feature("docstring")  Epetra_MpiDistributor::DoReversePosts "int
Epetra_MpiDistributor::DoReversePosts(char *export_objs, int obj_size,
int *&sizes, int &len_import_objs, char *&import_objs)

Do reverse post of buffer of export objects (can do other local work
before executing Waits) ";

/*  Attribute Accessor Methods  */

%feature("docstring")  Epetra_MpiDistributor::NumReceives "int
Epetra_MpiDistributor::NumReceives() const

The number of procs from which we will receive data. ";

%feature("docstring")  Epetra_MpiDistributor::NumSends "int
Epetra_MpiDistributor::NumSends() const

The number of procs to which we will send data. ";

%feature("docstring")  Epetra_MpiDistributor::MaxSendLength "int
Epetra_MpiDistributor::MaxSendLength() const

Maximum number of values that this proc is sending to another single
proc. ";

%feature("docstring")  Epetra_MpiDistributor::TotalReceiveLength "int
Epetra_MpiDistributor::TotalReceiveLength() const

Total number of values that this proc is receiving from other procs.
";

%feature("docstring")  Epetra_MpiDistributor::ProcsFrom "const int*
Epetra_MpiDistributor::ProcsFrom() const

A list of procs sending values to this proc. ";

%feature("docstring")  Epetra_MpiDistributor::ProcsTo "const int*
Epetra_MpiDistributor::ProcsTo() const

A list of procs to which this proc is sending values. ";

%feature("docstring")  Epetra_MpiDistributor::LengthsFrom "const int*
Epetra_MpiDistributor::LengthsFrom() const

Number of values we're receiving from each proc.

We will receive LengthsFrom[i] values from proc ProcsFrom[i]. ";

%feature("docstring")  Epetra_MpiDistributor::LengthsTo "const int*
Epetra_MpiDistributor::LengthsTo() const

Number of values we're sending to each proc.

We will send LengthsTo[i] values to procs ProcsTo[i]. ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_MpiDistributor::Print "void
Epetra_MpiDistributor::Print(ostream &os) const ";


// File: classEpetra__MpiSmpComm.xml
%feature("docstring") Epetra_MpiSmpComm "

Epetra_MpiSmpComm: The Epetra MPI Shared Memory Parallel Communication
Class.

The Epetra_MpiSmpComm class is an implementation of Epetra_Comm that
encapsulates the general information and services needed for other
Epetra classes to run on a parallel computer using MPI and shared
memory threads. WARNING:  This is an experimental class that
marginally supported nested share memory parallelism within MPI
processes.

This class has been deprecated.

C++ includes: Epetra_MpiSmpComm.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")  Epetra_MpiSmpComm::Epetra_MpiSmpComm "Epetra_MpiSmpComm::Epetra_MpiSmpComm(MPI_Comm comm)

Epetra_MpiSmpComm MPI Constructor.

Creates a Epetra_MpiSmpComm instance for use with MPI. If no
specialized MPI communicator is needed, this constuctor can be called
with the argument MPI_COMM_WORLD. ";

%feature("docstring")  Epetra_MpiSmpComm::Epetra_MpiSmpComm "Epetra_MpiSmpComm::Epetra_MpiSmpComm(const Epetra_MpiSmpComm &Comm)

Epetra_MpiSmpComm Copy Constructor.

Makes an exact copy of an existing Epetra_MpiSmpComm instance. ";

%feature("docstring")  Epetra_MpiSmpComm::Clone "Epetra_Comm*
Epetra_MpiSmpComm::Clone() const

Clone method. ";

%feature("docstring")  Epetra_MpiSmpComm::~Epetra_MpiSmpComm "Epetra_MpiSmpComm::~Epetra_MpiSmpComm()

Epetra_MpiSmpComm Destructor.

Completely deletes a Epetra_MpiSmpComm object. WARNING:  Note: All
objects that depend on a Epetra_MpiSmpComm instance should be
destroyed prior to calling this function. ";

/*  Barrier Methods  */

%feature("docstring")  Epetra_MpiSmpComm::Barrier "void
Epetra_MpiSmpComm::Barrier() const

Epetra_MpiSmpComm Barrier function.

Causes each processor in the communicator to wait until all processors
have arrived. ";

/*  Broadcast Methods  */

%feature("docstring")  Epetra_MpiSmpComm::Broadcast "int
Epetra_MpiSmpComm::Broadcast(double *MyVals, int Count, int Root)
const

Epetra_MpiSmpComm Broadcast function.

Takes list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_MpiSmpComm::Broadcast "int
Epetra_MpiSmpComm::Broadcast(int *MyVals, int Count, int Root) const

Epetra_MpiSmpComm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_MpiSmpComm::Broadcast "int
Epetra_MpiSmpComm::Broadcast(long *MyVals, int Count, int Root) const

Epetra_MpiSmpComm Broadcast function.

Take list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

%feature("docstring")  Epetra_MpiSmpComm::Broadcast "int
Epetra_MpiSmpComm::Broadcast(char *MyVals, int Count, int Root) const

Epetra_MpiSmpComm Broadcast function.

Takes list of input values from the root processor and sends to all
other processors.

Parameters:
-----------

Values:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of Values.

Root:  In On entry, contains the processor from which all processors
will receive a copy of Values. ";

/*  Gather Methods  */

%feature("docstring")  Epetra_MpiSmpComm::GatherAll "int
Epetra_MpiSmpComm::GatherAll(double *MyVals, double *AllVals, int
Count) const

Epetra_MpiSmpComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_MpiSmpComm::GatherAll "int
Epetra_MpiSmpComm::GatherAll(int *MyVals, int *AllVals, int Count)
const

Epetra_MpiSmpComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_MpiSmpComm::GatherAll "int
Epetra_MpiSmpComm::GatherAll(long *MyVals, long *AllVals, int Count)
const

Epetra_MpiSmpComm All Gather function.

Take list of input values from all processors in the communicator and
creates an ordered contiguous list of those values on each processor.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

/*  Sum Methods  */

%feature("docstring")  Epetra_MpiSmpComm::SumAll "int
Epetra_MpiSmpComm::SumAll(double *PartialSums, double *GlobalSums, int
Count) const

Epetra_MpiSmpComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::SumAll "int
Epetra_MpiSmpComm::SumAll(int *PartialSums, int *GlobalSums, int
Count) const

Epetra_MpiSmpComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::SumAll "int
Epetra_MpiSmpComm::SumAll(long *PartialSums, long *GlobalSums, int
Count) const

Epetra_MpiSmpComm Global Sum function.

Take list of input values from all processors in the communicator,
computes the sum and returns the sum to all processors.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Max/Min Methods  */

%feature("docstring")  Epetra_MpiSmpComm::MaxAll "int
Epetra_MpiSmpComm::MaxAll(double *PartialMaxs, double *GlobalMaxs, int
Count) const

Epetra_MpiSmpComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::MaxAll "int
Epetra_MpiSmpComm::MaxAll(int *PartialMaxs, int *GlobalMaxs, int
Count) const

Epetra_MpiSmpComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::MaxAll "int
Epetra_MpiSmpComm::MaxAll(long *PartialMaxs, long *GlobalMaxs, int
Count) const

Epetra_MpiSmpComm Global Max function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::MinAll "int
Epetra_MpiSmpComm::MinAll(double *PartialMins, double *GlobalMins, int
Count) const

Epetra_MpiSmpComm Global Min function.

Take list of input values from all processors in the communicator,
computes the min and returns the min to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::MinAll "int
Epetra_MpiSmpComm::MinAll(int *PartialMins, int *GlobalMins, int
Count) const

Epetra_MpiSmpComm Global Min function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::MinAll "int
Epetra_MpiSmpComm::MinAll(long *PartialMins, long *GlobalMins, int
Count) const

Epetra_MpiSmpComm Global Min function.

Take list of input values from all processors in the communicator,
computes the max and returns the max to all processors.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Parallel Prefix Methods  */

%feature("docstring")  Epetra_MpiSmpComm::ScanSum "int
Epetra_MpiSmpComm::ScanSum(double *MyVals, double *ScanSums, int
Count) const

Epetra_MpiSmpComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::ScanSum "int
Epetra_MpiSmpComm::ScanSum(int *MyVals, int *ScanSums, int Count)
const

Epetra_MpiSmpComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_MpiSmpComm::ScanSum "int
Epetra_MpiSmpComm::ScanSum(long *MyVals, long *ScanSums, int Count)
const

Epetra_MpiSmpComm Scan Sum function.

Take list of input values from all processors in the communicator,
computes the scan sum and returns it to all processors such that
processor i contains the sum of values from processor 0 up to and
including processor i.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

/*  Attribute Accessor Methods  */

%feature("docstring")  Epetra_MpiSmpComm::Comm "MPI_Comm
Epetra_MpiSmpComm::Comm() const

Extract MPI Communicator from a Epetra_MpiSmpComm object. ";

%feature("docstring")  Epetra_MpiSmpComm::MyPID "int
Epetra_MpiSmpComm::MyPID() const

Return my process ID.

In MPI mode returns the rank of the calling process. In serial mode
returns 0. ";

%feature("docstring")  Epetra_MpiSmpComm::NumProc "int
Epetra_MpiSmpComm::NumProc() const

Returns total number of processes.

In MPI mode returns the size of the MPI communicator. In serial mode
returns 1. ";

/*  Gather/Scatter and Directory Constructors  */

%feature("docstring")  Epetra_MpiSmpComm::CreateDistributor "Epetra_Distributor * Epetra_MpiSmpComm::CreateDistributor() const

Create a distributor object. ";

%feature("docstring")  Epetra_MpiSmpComm::CreateDirectory "Epetra_Directory* Epetra_MpiSmpComm::CreateDirectory(const
Epetra_BlockMap &Map) const

Create a directory object for the given Epetra_BlockMap. ";

/*  MPI-specific Methods  */

%feature("docstring")  Epetra_MpiSmpComm::GetMpiTag "int
Epetra_MpiSmpComm::GetMpiTag() const

Acquire an MPI tag from the Epetra range of 24050-24099, increment
tag. ";

%feature("docstring")  Epetra_MpiSmpComm::GetMpiComm "MPI_Comm
Epetra_MpiSmpComm::GetMpiComm() const

Acquire an MPI tag from the Epetra range of 24050-24099, increment
tag. ";

/*  Experimental SMP cluster methods (not rigorously implemented)  */

%feature("docstring")  Epetra_MpiSmpComm::NodeBarrier "void
Epetra_MpiSmpComm::NodeBarrier() const

Epetra_MpiSmpComm Node Barrier function.

A no-op for a serial communicator. For MPI, it causes each process on
a given node in the communicator to wait until all processes on that
node have arrived.

This function can be used to select a subset of MPI processes that are
associated with a group of threaded processes and synchronize only
with this subset. ";

%feature("docstring")  Epetra_MpiSmpComm::MyThreadID "int
Epetra_MpiSmpComm::MyThreadID() const

Return my thread ID.

If SetMyThreadID was called to set a thread value, this function
returns the thread ID of the calling process. Otherwise returns 0. ";

%feature("docstring")  Epetra_MpiSmpComm::MyNodeID "int
Epetra_MpiSmpComm::MyNodeID() const

Return my node ID.

If SetMyNodeD was called to set a node value, this function returns
the thread ID of the calling process. Otherwise returns the same value
as MyPID(). ";

%feature("docstring")  Epetra_MpiSmpComm::SetNumThreads "int
Epetra_MpiSmpComm::SetNumThreads(int NumThreads)

Set number of threads on this node.

Sets the number of threads on the node that owns the calling process.
By default the number of threads is 1. ";

%feature("docstring")  Epetra_MpiSmpComm::NumThreads "int
Epetra_MpiSmpComm::NumThreads() const

Get number of threads on this node.

Sets the number of threads on the node that owns the calling process.
By default the number of threads is 1. ";

%feature("docstring")  Epetra_MpiSmpComm::SetMyThreadID "int
Epetra_MpiSmpComm::SetMyThreadID(int ThreadID)

Set my thread ID.

Sets the thread ID for the calling process. Can be used to facilitate
threaded programming across an MPI application by allowing multiple
MPI processes to be considered threads of a virtual shared memory
process. Threads and nodes should be used together. By default the
thread ID is zero. ";

%feature("docstring")  Epetra_MpiSmpComm::SetMyNodeID "int
Epetra_MpiSmpComm::SetMyNodeID(int NodeID)

Set my node ID.

Sets the node ID for the calling process. Can be used to facilitate
threaded programming across an MPI application by associating several
MPI processes with a single node. By default, each MPI process is
associated with a single node with the same ID. ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_MpiSmpComm::Print "void
Epetra_MpiSmpComm::Print(ostream &os) const

Print method that implements Epetra_Object virtual Print method. ";

%feature("docstring")  Epetra_MpiSmpComm::PrintInfo "void
Epetra_MpiSmpComm::PrintInfo(ostream &os) const

Print method that implements Epetra_Comm virtual PrintInfo method. ";


// File: classEpetra__MpiSmpCommData.xml
%feature("docstring") Epetra_MpiSmpCommData "

Epetra_MpiSmpCommData: The Epetra Mpi Shared Memory
ParallelCommunication Data Class.

The Epetra_MpiSmpCommData class is an implementation detail of
Epetra_MpiSmpComm. It is reference-counted, and can be shared by
multiple Epetra_MpiSmpComm instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_MpiSmpCommData.h ";

/*  Constructor/Destructor Methods  */


// File: classEpetra__MultiVector.xml
%feature("docstring") Epetra_MultiVector "

Epetra_MultiVector: A class for constructing and using dense multi-
vectors, vectors and matrices in parallel.

The Epetra_MultiVector class enables the construction and use of real-
valued, double- precision dense vectors, multi-vectors, and matrices
in a distributed memory environment. The dimensions and distribution
of the dense multi-vectors is determined in part by a Epetra_Comm
object, a Epetra_Map (or Epetra_LocalMap or Epetra_BlockMap) and the
number of vectors passed to the constructors described below.

There are several concepts that important for understanding the
Epetra_MultiVector class:

Multi-vectors, Vectors and Matrices. Vector - A list of real-valued,
double-precision numbers. Also a multi-vector with one vector.

Multi-Vector - A collection of one or more vectors, all having the
same length and distribution.

(Dense) Matrix - A special form of multi-vector such that stride in
memory between any two consecutive vectors in the multi-vector is the
same for all vectors. This is identical to a two-dimensional array in
Fortran and plays an important part in high performance computations.

Distributed Global vs. Replicated Local. Distributed Global Multi-
vectors - In most instances, a multi-vector will be partitioned across
multiple memory images associated with multiple processors. In this
case, there is a unique copy of each element and elements are spread
across all processors specified by the Epetra_Comm communicator.

Replicated Local Multi-vectors - Some algorithms use multi-vectors
that are too small to be distributed across all processors, the
Hessenberg matrix in a GMRES computation. In other cases, such as with
block iterative methods, block dot product functions produce small
dense matrices that are required by all processors. Replicated local
multi-vectors handle these types of situation.

Multi-vector Functions vs. Dense Matrix Functions. Multi-vector
functions - These functions operate simultaneously but independently
on each vector in the multi-vector and produce individual results for
each vector.

Dense matrix functions - These functions operate on the multi-vector
as a matrix, providing access to selected dense BLAS and LAPACK
operations.

Constructing Epetra_MultiVectors

Except for the basic constructor and copy constructor,
Epetra_MultiVector constructors have two data access modes: Copy mode
- Allocates memory and makes a copy of the user-provided data. In this
case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the multi-
vector.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  All Epetra_MultiVector constructors require a map
argument that describes the layout of elements on the parallel
machine. Specifically, map is a Epetra_Map, Epetra_LocalMap or
Epetra_BlockMap object describing the desired memory layout for the
multi-vector.

There are six different Epetra_MultiVector constructors: Basic - All
values are zero.

Copy - Copy an existing multi-vector.

Copy from or make view of two-dimensional Fortran style array.

Copy from or make view of an array of pointers.

Copy or make view of a list of vectors from another Epetra_MultiVector
object.

Copy or make view of a range of vectors from another
Epetra_MultiVector object.

Extracting Data from Epetra_MultiVectors

Once a Epetra_MultiVector is constructed, it is possible to extract a
copy of the values or create a view of them.

WARNING:  ExtractView functions are extremely dangerous from a data
hiding perspective. For both ExtractView fuctions, there is a
corresponding ExtractCopy function. We strongly encourage users to
develop code using ExtractCopy functions first and only use the
ExtractView functions in a secondary optimization phase.  There are
four Extract functions: ExtractCopy - Copy values into a user-provided
two-dimensional array.

ExtractCopy - Copy values into a user-provided array of pointers.

ExtractView - Set user-provided two-dimensional array parameters to
point to Epetra_MultiVector data.

ExtractView - Set user-provided array of pointer parameters to point
to Epetra_MultiVector data.

Vector, Matrix and Utility Functions

Once a Epetra_MultiVector is constructed, a variety of mathematical
functions can be applied to the individual vectors. Specifically: Dot
Products.

Vector Updates.

p Norms.

Weighted Norms.

Minimum, Maximum and Average Values.

In addition, a matrix-matrix multiply function supports a variety of
operations on any viable combination of global distributed and local
replicated multi-vectors using calls to DGEMM, a high performance
kernel for matrix operations. In the near future we will add support
for calls to other selected BLAS and LAPACK functions.

Counting Floating Point Operations

Each Epetra_MultiVector object keep track of the number of serial
floating point operations performed using the specified object as the
this argument to the function. The Flops() function returns this
number as a double precision number. Using this information, in
conjunction with the Epetra_Time class, one can get accurate parallel
performance numbers. The ResetFlops() function resets the floating
point counter.

WARNING:  A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is
required for all Epetra_MultiVector constructors.

C++ includes: Epetra_MultiVector.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_MultiVector::Epetra_MultiVector "Epetra_MultiVector::Epetra_MultiVector(const Epetra_BlockMap &Map, int
NumVectors, bool zeroOut=true)

Basic Epetra_MultiVector constuctor.

Creates a Epetra_MultiVector object and, by default, fills with zero
values.

Parameters:
-----------

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

WARNING:  Note that, because Epetra_LocalMap derives from Epetra_Map
and Epetra_Map derives from Epetra_BlockMap, this constructor works
for all three types of Epetra map classes.

Parameters:
-----------

In:  NumVectors - Number of vectors in multi-vector.

In:  zeroOut - If true then the allocated memory will be zeroed out
initialy. If false then this memory will not be touched which can be
significantly faster.

Pointer to a Epetra_MultiVector. ";

%feature("docstring")  Epetra_MultiVector::Epetra_MultiVector "Epetra_MultiVector::Epetra_MultiVector(const Epetra_MultiVector
&Source)

Epetra_MultiVector copy constructor. ";

%feature("docstring")  Epetra_MultiVector::Epetra_MultiVector "Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const
Epetra_BlockMap &Map, double *A, int MyLDA, int NumVectors)

Set multi-vector values from two-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  A - Pointer to an array of double precision numbers. The first
vector starts at A. The second vector starts at A+MyLDA, the third at
A+2*MyLDA, and so on.

In:  MyLDA - The \"Leading Dimension\", or stride between vectors in
memory.

WARNING:  This value refers to the stride on the calling processor.
Thus it is a local quantity, not a global quantity.

Parameters:
-----------

In:  NumVectors - Number of vectors in multi-vector.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::Epetra_MultiVector "Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const
Epetra_BlockMap &Map, double **ArrayOfPointers, int NumVectors)

Set multi-vector values from array of pointers.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  ArrayOfPointers - An array of pointers such that
ArrayOfPointers[i] points to the memory location containing ith vector
to be copied.

In:  NumVectors - Number of vectors in multi-vector.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::Epetra_MultiVector "Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const
Epetra_MultiVector &Source, int *Indices, int NumVectors)

Set multi-vector values from list of vectors in an existing
Epetra_MultiVector.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Source - An existing fully constructed Epetra_MultiVector.

In:  Indices - Integer list of the vectors to copy.

In:  NumVectors - Number of vectors in multi-vector.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::Epetra_MultiVector "Epetra_MultiVector::Epetra_MultiVector(Epetra_DataAccess CV, const
Epetra_MultiVector &Source, int StartIndex, int NumVectors)

Set multi-vector values from range of vectors in an existing
Epetra_MultiVector.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Source - An existing fully constructed Epetra_MultiVector.

In:  StartIndex - First of the vectors to copy.

In:  NumVectors - Number of vectors in multi-vector.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::~Epetra_MultiVector "Epetra_MultiVector::~Epetra_MultiVector()

Epetra_MultiVector destructor. ";

/*  Post-construction modification routines  */

%feature("docstring")  Epetra_MultiVector::ReplaceGlobalValue "int
Epetra_MultiVector::ReplaceGlobalValue(int GlobalRow, int VectorIndex,
double ScalarValue)

Replace current value at the specified (GlobalRow, VectorIndex)
location with ScalarValue.

Replaces the existing value for a single entry in the multivector. The
specified global row must correspond to a GID owned by the map of the
multivector on the calling processor. In other words, this method does
not perform cross-processor communication.

If the map associated with this multivector is an Epetra_BlockMap,
only the first point entry associated with the global row will be
modified. To modify a different point entry, use the other version of
this method

Parameters:
-----------

In:  GlobalRow - Row of Multivector to modify in global index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if GlobalRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(). ";

%feature("docstring")  Epetra_MultiVector::ReplaceGlobalValue "int
Epetra_MultiVector::ReplaceGlobalValue(long long GlobalRow, int
VectorIndex, double ScalarValue) ";

%feature("docstring")  Epetra_MultiVector::ReplaceGlobalValue "int
Epetra_MultiVector::ReplaceGlobalValue(int GlobalBlockRow, int
BlockRowOffset, int VectorIndex, double ScalarValue)

Replace current value at the specified (GlobalBlockRow,
BlockRowOffset, VectorIndex) location with ScalarValue.

Replaces the existing value for a single entry in the multivector. The
specified global block row and block row offset must correspond to a
GID owned by the map of the multivector on the calling processor. In
other words, this method does not perform cross-processor
communication.

Parameters:
-----------

In:  GlobalBlockRow - BlockRow of Multivector to modify in global
index space.

In:  BlockRowOffset - Offset into BlockRow of Multivector to modify in
global index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if GlobalRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(), set to -2 if BlockRowOffset is out-of-range. ";

%feature("docstring")  Epetra_MultiVector::ReplaceGlobalValue "int
Epetra_MultiVector::ReplaceGlobalValue(long long GlobalBlockRow, int
BlockRowOffset, int VectorIndex, double ScalarValue) ";

%feature("docstring")  Epetra_MultiVector::SumIntoGlobalValue "int
Epetra_MultiVector::SumIntoGlobalValue(int GlobalRow, int VectorIndex,
double ScalarValue)

Adds ScalarValue to existing value at the specified (GlobalRow,
VectorIndex) location.

Sums the given value into the existing value for a single entry in the
multivector. The specified global row must correspond to a GID owned
by the map of the multivector on the calling processor. In other
words, this method does not perform cross-processor communication.

If the map associated with this multivector is an Epetra_BlockMap,
only the first point entry associated with the global row will be
modified. To modify a different point entry, use the other version of
this method

Parameters:
-----------

In:  GlobalRow - Row of Multivector to modify in global index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if GlobalRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(). ";

%feature("docstring")  Epetra_MultiVector::SumIntoGlobalValue "int
Epetra_MultiVector::SumIntoGlobalValue(long long GlobalRow, int
VectorIndex, double ScalarValue) ";

%feature("docstring")  Epetra_MultiVector::SumIntoGlobalValue "int
Epetra_MultiVector::SumIntoGlobalValue(int GlobalBlockRow, int
BlockRowOffset, int VectorIndex, double ScalarValue)

Adds ScalarValue to existing value at the specified (GlobalBlockRow,
BlockRowOffset, VectorIndex) location.

Sums the given value into the existing value for a single entry in the
multivector. The specified global block row and block row offset must
correspond to a GID owned by the map of the multivector on the calling
processor. In other words, this method does not perform cross-
processor communication.

Parameters:
-----------

In:  GlobalBlockRow - BlockRow of Multivector to modify in global
index space.

In:  BlockRowOffset - Offset into BlockRow of Multivector to modify in
global index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if GlobalRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(), set to -2 if BlockRowOffset is out-of-range. ";

%feature("docstring")  Epetra_MultiVector::SumIntoGlobalValue "int
Epetra_MultiVector::SumIntoGlobalValue(long long GlobalBlockRow, int
BlockRowOffset, int VectorIndex, double ScalarValue) ";

%feature("docstring")  Epetra_MultiVector::ReplaceMyValue "int
Epetra_MultiVector::ReplaceMyValue(int MyRow, int VectorIndex, double
ScalarValue)

Replace current value at the specified (MyRow, VectorIndex) location
with ScalarValue.

Replaces the existing value for a single entry in the multivector. The
specified local row must correspond to a GID owned by the map of the
multivector on the calling processor. In other words, this method does
not perform cross-processor communication.

This method is intended for use with vectors based on an Epetra_Map.
If used on a vector based on a non-trivial Epetra_BlockMap, this will
update only block row 0, i.e.

Epetra_MultiVector::ReplaceMyValue ( MyRow, VectorIndex, ScalarValue )
is equivalent to: Epetra_MultiVector::ReplaceMyValue ( 0, MyRow,
VectorIndex, ScalarValue )

Parameters:
-----------

In:  MyRow - Row of Multivector to modify in local index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if MyRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(). ";

%feature("docstring")  Epetra_MultiVector::ReplaceMyValue "int
Epetra_MultiVector::ReplaceMyValue(int MyBlockRow, int BlockRowOffset,
int VectorIndex, double ScalarValue)

Replace current value at the specified (MyBlockRow, BlockRowOffset,
VectorIndex) location with ScalarValue.

Replaces the existing value for a single entry in the multivector. The
specified local block row and block row offset must correspond to a
GID owned by the map of the multivector on the calling processor. In
other words, this method does not perform cross-processor
communication.

Parameters:
-----------

In:  MyBlockRow - BlockRow of Multivector to modify in local index
space.

In:  BlockRowOffset - Offset into BlockRow of Multivector to modify in
local index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if MyRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(), set to -2 if BlockRowOffset is out-of-range. ";

%feature("docstring")  Epetra_MultiVector::SumIntoMyValue "int
Epetra_MultiVector::SumIntoMyValue(int MyRow, int VectorIndex, double
ScalarValue)

Adds ScalarValue to existing value at the specified (MyRow,
VectorIndex) location.

Sums the given value into the existing value for a single entry in the
multivector. The specified local row must correspond to a GID owned by
the map of the multivector on the calling processor. In other words,
this method does not perform cross-processor communication.

If the map associated with this multivector is an Epetra_BlockMap,
only the first point entry associated with the local row will be
modified. To modify a different point entry, use the other version of
this method

Parameters:
-----------

In:  MyRow - Row of Multivector to modify in local index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if MyRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(). ";

%feature("docstring")  Epetra_MultiVector::SumIntoMyValue "int
Epetra_MultiVector::SumIntoMyValue(int MyBlockRow, int BlockRowOffset,
int VectorIndex, double ScalarValue)

Adds ScalarValue to existing value at the specified (MyBlockRow,
BlockRowOffset, VectorIndex) location.

Sums the given value into the existing value for a single entry in the
multivector. The specified local block row and block row offset must
correspond to a GID owned by the map of the multivector on the calling
processor. In other words, this method does not perform cross-
processor communication.

Parameters:
-----------

In:  MyBlockRow - BlockRow of Multivector to modify in local index
space.

In:  BlockRowOffset - Offset into BlockRow of Multivector to modify in
local index space.

In:  VectorIndex - Vector within MultiVector that should to modify.

In:  ScalarValue - Value to add to existing value.

Integer error code, set to 0 if successful, set to 1 if MyRow not
associated with calling processor set to -1 if VectorIndex >=
NumVectors(), set to -2 if BlockRowOffset is out-of-range. ";

%feature("docstring")  Epetra_MultiVector::PutScalar "int
Epetra_MultiVector::PutScalar(double ScalarConstant)

Initialize all values in a multi-vector with constant value.

Parameters:
-----------

In:  ScalarConstant - Value to use.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Random "int
Epetra_MultiVector::Random()

Set multi-vector values to random numbers.

MultiVector uses the random number generator provided by Epetra_Util.
The multi-vector values will be set to random values on the interval
(-1.0, 1.0).

Integer error code, set to 0 if successful. ";

/*  Extraction methods  */

%feature("docstring")  Epetra_MultiVector::ExtractCopy "int
Epetra_MultiVector::ExtractCopy(double *A, int MyLDA) const

Put multi-vector values into user-provided two-dimensional array.

Parameters:
-----------

Out:  A - Pointer to memory space that will contain the multi-vector
values. The first vector will be copied to the memory pointed to by A.
The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so
on.

In:  MyLDA - The \"Leading Dimension\", or stride between vectors in
memory.

WARNING:  This value refers to the stride on the calling processor.
Thus it is a local quantity, not a global quantity.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::ExtractCopy "int
Epetra_MultiVector::ExtractCopy(double **ArrayOfPointers) const

Put multi-vector values into user-provided array of pointers.

Parameters:
-----------

Out:  ArrayOfPointers - An array of pointers to memory space that will
contain the multi-vector values, such that ArrayOfPointers[i] points
to the memory location where the ith vector to be copied.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::ExtractView "int
Epetra_MultiVector::ExtractView(double **A, int *MyLDA) const

Set user-provided addresses of A and MyLDA.

Parameters:
-----------

A:  (Out) - Address of a pointer to that will be set to point to the
values of the multi-vector. The first vector will be at the memory
pointed to by A. The second vector starts at A+MyLDA, the third at
A+2*MyLDA, and so on.

MyLDA:  (Out) - Address of the \"Leading Dimension\", or stride
between vectors in memory.

WARNING:  This value refers to the stride on the calling processor.
Thus it is a local quantity, not a global quantity.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_MultiVector::ExtractView "int
Epetra_MultiVector::ExtractView(double ***ArrayOfPointers) const

Set user-provided addresses of ArrayOfPointers.

Parameters:
-----------

ArrayOfPointers:  (Out) - Address of array of pointers to memory space
that will set to the multi-vector array of pointers, such that
ArrayOfPointers[i] points to the memory location where the ith vector
is located.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

/*  Mathematical methods  */

%feature("docstring")  Epetra_MultiVector::Dot "int
Epetra_MultiVector::Dot(const Epetra_MultiVector &A, double *Result)
const

Computes dot product of each corresponding pair of vectors.

Parameters:
-----------

In:  A - Multi-vector to be used with the \"\\\\e this\" multivector.

Out:  Result - Result[i] will contain the ith dot product result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Abs "int
Epetra_MultiVector::Abs(const Epetra_MultiVector &A)

Puts element-wise absolute values of input Multi-vector in target.

Parameters:
-----------

In:  A - Input Multi-vector.

Out:   this will contain the absolute values of the entries of A.

Integer error code, set to 0 if successful.  Note: It is possible to
use the same argument for A and this. ";

%feature("docstring")  Epetra_MultiVector::Reciprocal "int
Epetra_MultiVector::Reciprocal(const Epetra_MultiVector &A)

Puts element-wise reciprocal values of input Multi-vector in target.

Parameters:
-----------

In:  A - Input Multi-vector.

Out:   this will contain the element-wise reciprocal values of the
entries of A.

Integer error code, set to 0 if successful. Returns 2 if some entry is
too small, but not zero. Returns 1 if some entry is zero.  Note: It is
possible to use the same argument for A and this. Also, if a given
value of A is smaller than Epetra_DoubleMin (defined in
Epetra_Epetra.h), but nonzero, then the return code is 2. If an entry
is zero, the return code is 1. However, in all cases the reciprocal
value is still used, even if a NaN is the result. ";

%feature("docstring")  Epetra_MultiVector::Scale "int
Epetra_MultiVector::Scale(double ScalarValue)

Scale the current values of a multi-vector, this = ScalarValue* this.

Parameters:
-----------

In:  ScalarValue - Scale value.

Out:   This - Multi-vector with scaled values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Scale "int
Epetra_MultiVector::Scale(double ScalarA, const Epetra_MultiVector &A)

Replace multi-vector values with scaled values of A, this = ScalarA*A.

Parameters:
-----------

In:  ScalarA - Scale value.

In:  A - Multi-vector to copy.

Out:   This - Multi-vector with values overwritten by scaled values of
A.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Update "int
Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector
&A, double ScalarThis)

Update multi-vector values with scaled values of A, this = ScalarThis*
this + ScalarA*A.

Parameters:
-----------

In:  ScalarA - Scale value for A.

In:  A - Multi-vector to add.

In:  ScalarThis - Scale value for this.

Out:   This - Multi-vector with updatede values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Update "int
Epetra_MultiVector::Update(double ScalarA, const Epetra_MultiVector
&A, double ScalarB, const Epetra_MultiVector &B, double ScalarThis)

Update multi-vector with scaled values of A and B, this = ScalarThis*
this + ScalarA*A + ScalarB*B.

Parameters:
-----------

In:  ScalarA - Scale value for A.

In:  A - Multi-vector to add.

In:  ScalarB - Scale value for B.

In:  B - Multi-vector to add.

In:  ScalarThis - Scale value for this.

Out:   This - Multi-vector with updatede values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Norm1 "int
Epetra_MultiVector::Norm1(double *Result) const

Compute 1-norm of each vector in multi-vector.

Parameters:
-----------

Out:  Result - Result[i] contains 1-norm of ith vector.

WARNING:  Map of the this multivector must have unique GIDs
(UniqueGIDs() must return true).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Norm2 "int
Epetra_MultiVector::Norm2(double *Result) const

Compute 2-norm of each vector in multi-vector.

Parameters:
-----------

Out:  Result - Result[i] contains 2-norm of ith vector.

WARNING:  Map of the this multivector must have unique GIDs
(UniqueGIDs() must return true).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::NormInf "int
Epetra_MultiVector::NormInf(double *Result) const

Compute Inf-norm of each vector in multi-vector.

Parameters:
-----------

Out:  Result - Result[i] contains Inf-norm of ith vector.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::NormWeighted "int
Epetra_MultiVector::NormWeighted(const Epetra_MultiVector &Weights,
double *Result) const

Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.

Parameters:
-----------

In:  Weights - Multi-vector of weights. If Weights contains a single
vector, that vector will be used as the weights for all vectors of
this. Otherwise, Weights should have the same number of vectors as
this.

Out:  Result - Result[i] contains the weighted 2-norm of ith vector.
Specifically if we denote the ith vector in the multivector by $x$,
and the ith weight vector by $w$ and let j represent the jth entry of
each vector, on return Result[i] will contain the following result:
\\\\[\\\\sqrt{(1/n)\\\\sum_{j=1}^n(x_j/w_j)^2}\\\\], where $n$ is the
global length of the vectors.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::MinValue "int
Epetra_MultiVector::MinValue(double *Result) const

Compute minimum value of each vector in multi-vector.

Note that the vector contents must be already initialized for this
function to compute a well-defined result. The length of the vector
need not be greater than zero on all processors. If length is greater
than zero on any processor then a valid result will be computed.

Parameters:
-----------

Out:  Result - Result[i] contains minimum value of ith vector.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::MaxValue "int
Epetra_MultiVector::MaxValue(double *Result) const

Compute maximum value of each vector in multi-vector.

Note that the vector contents must be already initialized for this
function to compute a well-defined result. The length of the vector
need not be greater than zero on all processors. If length is greater
than zero on any processor then a valid result will be computed.

Parameters:
-----------

Out:  Result - Result[i] contains maximum value of ith vector.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::MeanValue "int
Epetra_MultiVector::MeanValue(double *Result) const

Compute mean (average) value of each vector in multi-vector.

Parameters:
-----------

Out:  Result - Result[i] contains mean value of ith vector.

WARNING:  Map of the this multivector must have unique GIDs
(UniqueGIDs() must return true).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Multiply "int
Epetra_MultiVector::Multiply(char TransA, char TransB, double
ScalarAB, const Epetra_MultiVector &A, const Epetra_MultiVector &B,
double ScalarThis)

Matrix-Matrix multiplication, this = ScalarThis* this + ScalarAB*A*B.

This function performs a variety of matrix-matrix multiply operations,
interpreting the Epetra_MultiVectors ( this-aka C , A and B) as 2D
matrices. Variations are due to the fact that A, B and C can be local
replicated or global distributed Epetra_MultiVectors and that we may
or may not operate with the transpose of A and B. Possible cases are:
Total of 32 case (2^5).     Num     OPERATIONS case  Notes     1)
C(local) = A^X(local) * B^X(local)  4 (X=Transpose or Not, No comm
needed)      2) C(local) = A^T(distr) * B (distr)  1   (2D dot
product, replicate C)     3) C(distr) = A (distr) * B^X(local)  2
(2D vector update, no comm needed)      Note that the following
operations are not meaningful for      1D distributions:      1)
C(local) = A^T(distr) * B^T(distr)  1     2) C(local) = A  (distr) *
B^X(distr)  2     3) C(distr) = A^X(local) * B^X(local)  4     4)
C(distr) = A^X(local) * B^X(distr)  4     5) C(distr) = A^T(distr) *
B^X(local)  2     6) C(local) = A^X(distr) * B^X(local)  4     7)
C(distr) = A^X(distr) * B^X(local)  4     8) C(local) = A^X(local) *
B^X(distr)  4

Parameters:
-----------

In:  TransA - Operate with the transpose of A if = 'T', else no
transpose if = 'N'.

In:  TransB - Operate with the transpose of B if = 'T', else no
transpose if = 'N'.

In:  ScalarAB - Scalar to multiply with A*B.

In:  A - Multi-vector.

In:  B - Multi-vector.

In:  ScalarThis - Scalar to multiply with this.

WARNING:  Map of the distributed multivectors must have unique GIDs
(UniqueGIDs() must return true).

Integer error code, set to 0 if successful.

WARNING:  {Each multi-vector A, B and this is checked if it has
constant stride using the ConstantStride() query function. If it does
not have constant stride, a temporary copy is made and used for the
computation. This activity is transparent to the user, except that
there is memory and computation overhead. All temporary space is
deleted prior to exit.}

A.ConstantStride() || B.ConstantStride() ) EPETRA_CHK_ERR(-1); //
Return error ";

%feature("docstring")  Epetra_MultiVector::Multiply "int
Epetra_MultiVector::Multiply(double ScalarAB, const Epetra_MultiVector
&A, const Epetra_MultiVector &B, double ScalarThis)

Multiply a Epetra_MultiVector with another, element-by-element.

This function supports diagonal matrix multiply. A is usually a single
vector while B and this may have one or more columns. Note that B and
this must have the same shape. A can be one vector or have the same
shape as B. The actual computation is this = ScalarThis * this +
ScalarAB * B @ A where @ denotes element-wise multiplication. ";

%feature("docstring")  Epetra_MultiVector::ReciprocalMultiply "int
Epetra_MultiVector::ReciprocalMultiply(double ScalarAB, const
Epetra_MultiVector &A, const Epetra_MultiVector &B, double ScalarThis)

Multiply a Epetra_MultiVector by the reciprocal of another, element-
by-element.

This function supports diagonal matrix scaling. A is usually a single
vector while B and this may have one or more columns. Note that B and
this must have the same shape. A can be one vector or have the same
shape as B. The actual computation is this = ScalarThis * this +
ScalarAB * B @ A where @ denotes element-wise division. ";

/*  Random number utilities  */

%feature("docstring")  Epetra_MultiVector::SetSeed "int
Epetra_MultiVector::SetSeed(unsigned int Seed_in)

Set seed for Random function.

Parameters:
-----------

In:  Seed - Should be an integer on the interval (0, 2^31-1).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_MultiVector::Seed "unsigned int
Epetra_MultiVector::Seed()

Get seed from Random function.

Current random number seed. ";

/*  Overloaded operators  */

/*  Attribute access functions  */

%feature("docstring")  Epetra_MultiVector::NumVectors "int
Epetra_MultiVector::NumVectors() const

Returns the number of vectors in the multi-vector. ";

%feature("docstring")  Epetra_MultiVector::MyLength "int
Epetra_MultiVector::MyLength() const

Returns the local vector length on the calling processor of vectors in
the multi-vector. ";

%feature("docstring")  Epetra_MultiVector::GlobalLength "int
Epetra_MultiVector::GlobalLength() const

Returns the global vector length of vectors in the multi-vector. ";

%feature("docstring")  Epetra_MultiVector::GlobalLength64 "long long
Epetra_MultiVector::GlobalLength64() const ";

%feature("docstring")  Epetra_MultiVector::Stride "int
Epetra_MultiVector::Stride() const

Returns the stride between vectors in the multi-vector (only
meaningful if ConstantStride() is true). ";

%feature("docstring")  Epetra_MultiVector::ConstantStride "bool
Epetra_MultiVector::ConstantStride() const

Returns true if this multi-vector has constant stride between vectors.
";

/*  I/O methods  */

%feature("docstring")  Epetra_MultiVector::Print "void
Epetra_MultiVector::Print(ostream &os) const

Print method. ";

/*  Expert-only unsupported methods  */

%feature("docstring")  Epetra_MultiVector::ResetView "int
Epetra_MultiVector::ResetView(double **ArrayOfPointers)

Reset the view of an existing multivector to point to new user data.

Allows the (very) light-weight replacement of multivector values for
an existing multivector that was constructed using an
Epetra_DataAccess mode of View. No checking is performed to see if the
array of values passed in contains valid data. It is assumed that the
user has verified the integrity of data before calling this method.
This method is useful for situations where a multivector is needed for
use with an Epetra operator or matrix and the user is not passing in a
multivector, or the multivector is being passed in with another map
that is not exactly compatible with the operator, but has the correct
number of entries.

This method is used by AztecOO and Ifpack in the matvec, and solve
methods to improve performance and reduce repeated calls to
constructors and destructors.

Parameters:
-----------

ArrayOfPointers:  Contains the array of pointers containing the
multivector data.

Integer error code, set to 0 if successful, -1 if the multivector was
not created as a View.

WARNING:  This method is extremely dangerous and should only be used
by experts. ";

%feature("docstring")  Epetra_MultiVector::Values "double*
Epetra_MultiVector::Values() const

Get pointer to MultiVector values. ";

%feature("docstring")  Epetra_MultiVector::Pointers "double**
Epetra_MultiVector::Pointers() const

Get pointer to individual vector pointers. ";

%feature("docstring")  Epetra_MultiVector::ReplaceMap "int
Epetra_MultiVector::ReplaceMap(const Epetra_BlockMap &map)

Replace map, only if new map has same point-structure as current map.
return 0 if map is replaced, -1 if not. ";

%feature("docstring")  Epetra_MultiVector::Reduce "int
Epetra_MultiVector::Reduce() ";


// File: classEpetra__Object.xml
%feature("docstring") Epetra_Object "

Epetra_Object: The base Epetra class.

The Epetra_Object class provides capabilities common to all Epetra
objects, such as a label that identifies an object instance, constant
definitions, enum types.

C++ includes: Epetra_Object.h ";

/*  Constructors/destructor  */

%feature("docstring")  Epetra_Object::Epetra_Object "Epetra_Object::Epetra_Object(int TracebackModeIn=-1, bool
set_label=true)

Epetra_Object Constructor.

Epetra_Object is the primary base class in Epetra. All Epetra class
are derived from it, directly or indirectly. This class is seldom used
explictly. ";

%feature("docstring")  Epetra_Object::Epetra_Object "Epetra_Object::Epetra_Object(const char *const Label, int
TracebackModeIn=-1)

Epetra_Object Constructor.

Creates a Epetra_Object with the given label. ";

%feature("docstring")  Epetra_Object::Epetra_Object "Epetra_Object::Epetra_Object(const Epetra_Object &Object)

Epetra_Object Copy Constructor.

Makes an exact copy of an existing Epetra_Object instance. ";

%feature("docstring")  Epetra_Object::~Epetra_Object "Epetra_Object::~Epetra_Object()

Epetra_Object Destructor.

Completely deletes a Epetra_Object object. ";

/*  Attribute set/get methods  */

%feature("docstring")  Epetra_Object::SetLabel "void
Epetra_Object::SetLabel(const char *const Label)

Epetra_Object Label definition using char *.

Defines the label used to describe the this object. ";

%feature("docstring")  Epetra_Object::Label "const char *
Epetra_Object::Label() const

Epetra_Object Label access funtion.

Returns the string used to define this object. ";

%feature("docstring")  Epetra_Object::SetTracebackMode "void
Epetra_Object::SetTracebackMode(int TracebackModeValue)

Set the value of the Epetra_Object error traceback report mode.

Sets the integer error traceback behavior. TracebackMode controls
whether or not traceback information is printed when run time integer
errors are detected:

<= 0 - No information report

= 1 - Fatal (negative) values are reported

>= 2 - All values (except zero) reported.

Default is set to 1. ";

%feature("docstring")  Epetra_Object::GetTracebackMode "int
Epetra_Object::GetTracebackMode()

Get the value of the Epetra_Object error report mode. ";

%feature("docstring")  Epetra_Object::GetTracebackStream "std::ostream & Epetra_Object::GetTracebackStream()

Get the output stream for error reporting. ";

/*  Miscellaneous  */

%feature("docstring")  Epetra_Object::Print "void
Epetra_Object::Print(ostream &os) const

Print object to an output stream Print method ";

%feature("docstring")  Epetra_Object::ReportError "int
Epetra_Object::ReportError(const string Message, int ErrorCode) const

Error reporting method. ";


// File: classEpetra__OffsetIndex.xml
%feature("docstring") Epetra_OffsetIndex "

Epetra_OffsetIndex: This class builds index for efficient mapping of
data from one Epetra_CrsGraph based object to another.

Epetra_OffsetIndex generates and index of offsets allowing direct
access to data for Import/Export operations on Epetra_CrsGraph based
objects such as Epetra_CrsMatrix.

C++ includes: Epetra_OffsetIndex.h ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_OffsetIndex::Print "void
Epetra_OffsetIndex::Print(ostream &os) const

Print object to an output stream Print method ";

%feature("docstring")  Epetra_OffsetIndex::Epetra_OffsetIndex "Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_CrsGraph
&SourceGraph, const Epetra_CrsGraph &TargetGraph, Epetra_Import
&Importer)

Constructs a Epetra_OffsetIndex object from the graphs and an
importer. ";

%feature("docstring")  Epetra_OffsetIndex::Epetra_OffsetIndex "Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_CrsGraph
&SourceGraph, const Epetra_CrsGraph &TargetGraph, Epetra_Export
&Exporter)

Constructs a Epetra_OffsetIndex object from the graphs and an
exporter. ";

%feature("docstring")  Epetra_OffsetIndex::Epetra_OffsetIndex "Epetra_OffsetIndex::Epetra_OffsetIndex(const Epetra_OffsetIndex
&Indexor)

Epetra_OffsetIndex copy constructor. ";

%feature("docstring")  Epetra_OffsetIndex::~Epetra_OffsetIndex "Epetra_OffsetIndex::~Epetra_OffsetIndex(void)

Epetra_OffsetIndex destructor. ";

%feature("docstring")  Epetra_OffsetIndex::SameOffsets "int**
Epetra_OffsetIndex::SameOffsets() const

Accessor. ";

%feature("docstring")  Epetra_OffsetIndex::PermuteOffsets "int**
Epetra_OffsetIndex::PermuteOffsets() const

Accessor. ";

%feature("docstring")  Epetra_OffsetIndex::RemoteOffsets "int**
Epetra_OffsetIndex::RemoteOffsets() const

Accessor. ";


// File: classEpetra__Operator.xml
%feature("docstring") Epetra_Operator "

Epetra_Operator: A pure virtual class for using real-valued double-
precision operators.

The Epetra_Operator class is a pure virtual class (specifies interface
only) that enable the use of real-valued double-precision operators.
It is currently implemented by both the Epetra_CrsMatrix and
Epetra_VbrMatrix classes and the Ifpack_CrsRiluk preconditioner class.

C++ includes: Epetra_Operator.h ";

/*  Destructor  */

%feature("docstring")  Epetra_Operator::~Epetra_Operator "virtual
Epetra_Operator::~Epetra_Operator()

Destructor. ";

/*  Attribute set methods  */

%feature("docstring")  Epetra_Operator::SetUseTranspose "virtual int
Epetra_Operator::SetUseTranspose(bool UseTranspose)=0

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose -If true, multiply by the transpose of operator,
otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

/*  Mathematical functions  */

%feature("docstring")  Epetra_Operator::Apply "virtual int
Epetra_Operator::Apply(const Epetra_MultiVector &X, Epetra_MultiVector
&Y) const =0

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_Operator::ApplyInverse "virtual int
Epetra_Operator::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const =0

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful.

WARNING:  In order to work with AztecOO, any implementation of this
method must support the case where X and Y are the same object. ";

%feature("docstring")  Epetra_Operator::NormInf "virtual double
Epetra_Operator::NormInf() const =0

Returns the infinity norm of the global matrix. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_Operator::Label "virtual const char*
Epetra_Operator::Label() const =0

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_Operator::UseTranspose "virtual bool
Epetra_Operator::UseTranspose() const =0

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_Operator::HasNormInf "virtual bool
Epetra_Operator::HasNormInf() const =0

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Epetra_Operator::Comm "virtual const
Epetra_Comm& Epetra_Operator::Comm() const =0

Returns a pointer to the Epetra_Comm communicator associated with this
operator. ";

%feature("docstring")  Epetra_Operator::OperatorDomainMap "virtual
const Epetra_Map& Epetra_Operator::OperatorDomainMap() const =0

Returns the Epetra_Map object associated with the domain of this
operator. ";

%feature("docstring")  Epetra_Operator::OperatorRangeMap "virtual
const Epetra_Map& Epetra_Operator::OperatorRangeMap() const =0

Returns the Epetra_Map object associated with the range of this
operator. ";


// File: classEpetra__OskiError.xml
%feature("docstring") Epetra_OskiError "

Epetra_OskiError: The Epetra OSKI Class to provide access to get and
set error handling routines in OSKI.

C++ includes: Epetra_OskiError.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_OskiError::Epetra_OskiError "Epetra_OskiError::Epetra_OskiError()

Default Constructor. ";

%feature("docstring")  Epetra_OskiError::~Epetra_OskiError "virtual
Epetra_OskiError::~Epetra_OskiError()

Destructor. ";

/*  Set/Get  */

%feature("docstring")  Epetra_OskiError::OskiGetErrorHandler "Epetra_OskiError Epetra_OskiError::OskiGetErrorHandler()

Gets a pointer to the current error handler routine being used by
OSKI. ";

%feature("docstring")  Epetra_OskiError::OskiSetErrorHandler "void
Epetra_OskiError::OskiSetErrorHandler(Epetra_OskiError
&NewErrorHandler)

Sets the error handling routine to be used by OSKI to NewErrorHandler.
";


// File: classEpetra__OskiMatrix.xml
%feature("docstring") Epetra_OskiMatrix "

Epetra_OskiMatrix: A class for constructing and using OSKI Matrices
within Epetra. For information on known issues with OSKI see the
detailed description.

OSKI is a high-performance sparse matrix kernel package written by the
UC Berkeley Benchmarking and Optimization Group. The Epetra_OskiMatrix
class is a lightweight interface to allow Epetra users access to OSKI
functionality. This interface includes lightweight conversion of
Epetra matrices to OSKI matrices, runtime tuning based on parameter
list hints, and access to high-performance computational kernels to
perform matrix-vector/matrix multi-vector calculations.

The Epetra_OskiMatrix class provides access to the entire OSKI
interface. However, the following features are not fully implemented
in OSKI version oski-1.0.1h, and therefore are unsupported in the
Epetra_OskiMatrix class:

OSKI does not provide stock composed kernels. Hence, the tune function
must be called in order to see performance gains when using the
MatTransMatMultiply and MultiplyAndMatTransMultiply kernels.

The MatPowMultiply kernel does not work.

Optimized multivector kernels are not created by default when
installing OSKI.

The tune function cannot transform a (nearly) symmetric matrix to be
stored as such.

In order to use the $A^TA$ OSKI kernel (MatTransMatMultiply), in
oski/src/MBCSR/ata.c you must replace the lines with

OSKI does not convert between CSR and CSC when it could be profitable,
such as when performing $AA^T$ on a CSR matrix.

OSKI may be incompatible with the following architectures: Barcelona
(quad-core Opteron): errors during \"make install\" (confirmed with
OSKI developers)

single core Xeon: OSKI installs, but never transforms matrices. This
includes cases where other machines will transform the same matrices,
and where one would expect the matrix to be transformed, based on OSKI
tuning data.

C++ includes: Epetra_OskiMatrix.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_OskiMatrix::Epetra_OskiMatrix "Epetra_OskiMatrix::Epetra_OskiMatrix(const Epetra_OskiMatrix &Source)

Copy constructor. ";

%feature("docstring")  Epetra_OskiMatrix::Epetra_OskiMatrix "Epetra_OskiMatrix::Epetra_OskiMatrix(const Epetra_CrsMatrix &Source,
const Teuchos::ParameterList &List)

Constructor creates an Epetra_OskiMatrix from an Epetra_CrsMatrix.

Parameters:
-----------

Source:  (In) An Epetra_CrsMatrix that is to be wrapped as an
Epetra_OskiMatrix.

List:  (In) Any options or data wanted or needed for the conversion.

Pointer to an Epetra_OskiMatrix.  Options that can be passed to the
List are presented below. They are: \"<type> <option name> <default
value>: <description of purpose>\"

bool autotune false: If true, Epetra tries to set as many hints as
possible based on its knowledge of the matrix.

string matrixtype general: Other types that can be taken are:
uppertri, lowertri, uppersymm, lowersymm, fullsymm, upperherm,
lowerherm and fullherm.

bool diagstored false: If true, the diagonal entries are not stored in
the matrix and are all assumed to be 1.

bool zerobased false: If true, the array is zero based, as in C.
Otherwise, it is 1 based, as in Fortran.

bool sorted false: If true, all elements in the passed in array are
sorted.

bool unique false: If true, a value in a column only appears once in
each row.

bool deepcopy false: If true, when the OSKI matrix is created it will
be a deepcopy of the data in the function. ";

%feature("docstring")  Epetra_OskiMatrix::~Epetra_OskiMatrix "virtual
Epetra_OskiMatrix::~Epetra_OskiMatrix()

Destructor. ";

/*  Extract/Replace Values  */

%feature("docstring")  Epetra_OskiMatrix::ReplaceMyValues "int
Epetra_OskiMatrix::ReplaceMyValues(int MyRow, int NumEntries, double
*Values, int *Indices)

Replace current values with this list of entries for a given local row
of the matrix. Warning this could be expensive.

The reason this function could be expensive is its underlying
implementation. Both the OSKI and Epetra versions of the matrix must
be changed when the matrix has been permuted. When this is the case, a
call must be made to the Epetra ReplaceMyValues, and NumEntries calls
must be made to a function that changes the OSKI matrix's values one
at a time.

Parameters:
-----------

MyRow:  (In) Row number (in local coordinates) to put elements.

NumEntries:  (In) Number of entries.

Values:  (In) Values to enter.

Indices:  (In) Local column indices corresponding to the values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, Oski will fail. A positive
warning code may be returned, but this should be treated as a fatal
error; part of the data will be changed, and OSKI cannot support
adding in new data values.

IndicesAreLocal()==true

The given Values at the given Indices have been summed into the
entries of MyRow. ";

%feature("docstring")  Epetra_OskiMatrix::SumIntoMyValues "int
Epetra_OskiMatrix::SumIntoMyValues(int MyRow, int NumEntries, double
*Values, int *Indices)

Add this list of entries to existing values for a given local row of
the matrix. WARNING: this could be expensive.

The reason this function could be expensive is its underlying
implementation. Both the OSKI and Epetra versions of the Matrix must
be changed when the matrix has been permuted. When this is the case, a
call must be made to the Epetra SumIntoMyValues, and NumEntries calls
must be made to a function that changes the OSKI matrix's values one
at a time.

Parameters:
-----------

MyRow:  - (In) Row number (in local coordinates) to put elements.

NumEntries:  - (In) Number of entries.

Values:  - (In) Values to enter.

Indices:  - (In) Local column indices corresponding to values.

Integer error code, set to 0 if successful. Note that if the allocated
length of the row has to be expanded, a positive warning code may be
returned. This should be treated as a fatal error, as part of the data
will be changed, and OSKI cannot support adding in new data values.

IndicesAreLocal()==true

The given Values at the given Indices have been summed into the
entries of MyRow. ";

%feature("docstring")  Epetra_OskiMatrix::ExtractDiagonalCopy "int
Epetra_OskiMatrix::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Diagonal:  - (Out) Extracted main diagonal.

Integer error code, set to 0 if successful and non-zero if not.

Filled()==true

Unchanged. ";

%feature("docstring")  Epetra_OskiMatrix::ReplaceDiagonalValues "int
Epetra_OskiMatrix::ReplaceDiagonalValues(const Epetra_OskiVector
&Diagonal)

Replaces diagonal values of the matrix with those in the user-provided
vector.

This routine is meant to allow replacement of { existing} diagonal
values. If a diagonal value does not exist for a given row, the
corresponding value in the input Epetra_OskiVector will be ignored,
and the return code will be set to 1.

Parameters:
-----------

Diagonal:  - (In) New values to be placed in the main diagonal.

Integer error code, set to 0 if successful, set to 1 on the calling
processor if one or more diagonal entries not present in matrix. Other
error codes can be returned as well, indicating improperly constructed
matrices or vectors.

Filled()==true

Diagonal values have been replaced with the values of Diagonal. ";

/*  Computational methods  */

%feature("docstring")  Epetra_OskiMatrix::Multiply "int
Epetra_OskiMatrix::Multiply(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y) const

Performs a matrix vector multiply of y = this^TransA*x.

The vectors x and y can be either Epetra_Vectors or
Epetra_OskiVectors.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
computing the product.

x:  (In) The vector the matrix is multiplied by.

y:  (Out) The vector where the calculation result is stored.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Multiply "int
Epetra_OskiMatrix::Multiply(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y, double Alpha, double Beta=0.0) const

Performs a matrix vector multiply of y = Alpha*this^TransA*x + Beta*y.

The vectors x and y can be either Epetra_Vectors or
Epetra_OskiVectors.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
computing the product.

x:  (In) The vector the matrix is multiplied by.

y:  (In/Out) The vector where the calculation result is stored.

Alpha:  (In) A scalar constant used to scale x.

Beta:  (In) A scalar constant used to scale y.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Multiply "int
Epetra_OskiMatrix::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Performs a matrix multi-vector multiply of Y = this^TransA*X.

The multi-vectors X and Y can be either Epetra_MultiVectors or
Epetra_OskiMultiVectors.

Parameters:
-----------

TransA:  (In) If Trans = TRUE then use the transpose of the matrix in
computing the product.

X:  (In) The multi-vector the matrix is multiplied by.

Y:  (Out) The multi-vector where the calculation result is stored.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Multiply "int
Epetra_OskiMatrix::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y, double Alpha, double Beta=0.0) const

Performs a matrix multi-vector multiply of Y = Alpha*this^TransA*X +
Beta*Y.

The multi-vectors X and Y can be either Epetra_MultiVectors or
Epetra_OskiMultiVectors.

Parameters:
-----------

TransA:  (In) If Trans = TRUE then use the transpose of the matrix in
computing the product.

X:  (In) The multi-vector the matrix is multiplied by.

Y:  (In/Out) The multi-vector where the calculation result is stored.

Alpha:  (In) A scalar constant used to scale X.

Beta:  (In) A scalar constant used to scale Y.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Solve "int
Epetra_OskiMatrix::Solve(bool Upper, bool TransA, bool UnitDiagonal,
const Epetra_Vector &x, Epetra_Vector &y) const

Performs a triangular solve of y = (this^TransA)^-1*x where this is a
triangular matrix.

The vector x can be either be an Epetra_Vector or Epetra_OskiVector.
The OskiMatrix must already be triangular, and the UnitDiagonal
setting associated with it will be used.

Parameters:
-----------

Upper:  (In) This parameter is ignored, and is here only to match the
Epetra_CrsMatrix::Solve syntax.

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
solving the equations.

UnitDiagonal:  (In) This parameter is ignored only here to match the
Epetra_CrsMatrix::Solve syntax.

x:  (In) The vector solved against.

y:  (Out) The solution vector.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Solve "int
Epetra_OskiMatrix::Solve(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y, double Alpha=1.0) const

Performs a triangular solve of y = Alpha*(this^TransA)^-1*x where this
is a triangular matrix.

The vector x can be either be an Epetra_Vector or Epetra_OskiVector.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
solving the equations.

x:  (In) The vector solved against.

y:  (Out) The solution vector.

Alpha:  (In) A scalar constant used to scale x.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Solve "int
Epetra_OskiMatrix::Solve(bool Upper, bool TransA, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Performs a triangular solve of Y = (this^TransA)^-1*X where this is a
triangular matrix.

The vector X can be either be an Epetra_MultiVector or
Epetra_OskiMultiVector. The OskiMatrix must already be triangular, and
the UnitDiagonal setting associated with it will be used.

Parameters:
-----------

Upper:  (In) This parameter is ignored only here to match the
Epetra_CrsMatrix::Solve syntax.

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
solving the equations.

UnitDiagonal:  (In) This parameter is ignored only here to match the
Epetra_CrsMatrix::Solve syntax.

X:  (In) The multi-vector solved against.

Y:  (Out) The solution multi-vector.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::Solve "int
Epetra_OskiMatrix::Solve(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y, double Alpha=1.0) const

Performs a triangular solve of Y = Alpha*(this^TransA)^-1*X where this
is a triangular matrix.

The vector X can be either be an Epetra_MultiVector or
Epetra_OskiMultiVector.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
solving the equations.

X:  (In) The multi-vector solved against.

Y:  (Out) The solution multi-vector.

Alpha:  (In) A scalar constant used to scale X.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::MatTransMatMultiply "int
Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, const Epetra_Vector
&x, Epetra_Vector &y, Epetra_Vector *t, double Alpha=1.0, double
Beta=0.0) const

Performs two matrix vector multiplies of y = Alpha*this^TransA*this*x
+ Beta*y or y = Alpha*this*this^TransA*x + Beta*y.

The vectors x, y and t can be either Epetra_Vectors or
Epetra_OskiVectors. This composed routine is most commonly used in
linear least squares and bidiagonalization methods. The parallel
version of y = Alpha*this*this^TransA*x + Beta*y uses calls to the
Multiply routine under the hood, as it is not possible to perform both
multiplies automatically.

Parameters:
-----------

ATA:  (In) If TransA = TRUE then compute this^T*this*x otherwise
compute this*this^T*x.

x:  (In) The vector the matrix is multiplied by.

y:  (In/Out) The vector where the calculation result is stored.

t:  (Out) The vector where the result of the this*x is stored if
TransA = true and this^T*x is stored otherwise.

Alpha:  (In) A scalar constant used to scale x.

Beta:  (In) A scalar constant used to scale y.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::MatTransMatMultiply "int
Epetra_OskiMatrix::MatTransMatMultiply(bool ATA, const
Epetra_MultiVector &X, Epetra_MultiVector &Y, Epetra_MultiVector *T,
double Alpha=1.0, double Beta=0.0) const

Performs two matrix multi-vector multiplies of Y =
Alpha*this^TransA*this*X + Beta*Y or Y = Alpha*this*this^TransA*X +
Beta*Y.

The multi-vectors X, Y and T can be either Epetra_MultiVectors or
Epetra_OskiMultiVectors. This composed routine is most commonly used
in linear least squares and bidiagonalization methods. The parallel
version of Y = Alpha*this*this^TransA*X + Beta*Y uses calls to the
Multiply routine under the hood, as it is not possible to perform both
multiplies automatically.

Parameters:
-----------

ATA:  (In) If TransA = TRUE then compute this^T*this*X otherwise
compute this*this^T*X.

X:  (In) The vector the matrix is multiplied by.

Y:  (In/Out) The vector where the calculation result is stored.

T:  (Out) The multi-vector where the result of the this*X is stored if
TransA = true and this^T*X is stored otherwise.

Alpha:  (In) A scalar constant used to scale X.

Beta:  (In) A scalar constant used to scale Y.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::MultiplyAndMatTransMultiply
"int Epetra_OskiMatrix::MultiplyAndMatTransMultiply(bool TransA,
const Epetra_Vector &x, Epetra_Vector &y, const Epetra_Vector &w,
Epetra_Vector &z, double Alpha=1.0, double Beta=0.0, double Omega=1.0,
double Zeta=0.0) const

Performs the two matrix vector multiplies of y = Alpha*this*x + Beta*y
and z = Omega*this^TransA*w + Zeta*z.

The vectors x, y, w and z can be either Epetra_Vectors or
Epetra_OskiVectors. This composed routine is most commonly used in bi-
conjugate gradient calculations.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
computing the second product.

x:  (In) A vector the matrix is multiplied by.

y:  (In/Out) A vector where the calculation result of the first
multiply is stored.

w:  (In) A vector the matrix is multiplied by.

z:  (In/Out) A vector where the calculation result of the second
multiply is stored.

Alpha:  (In) A scalar constant used to scale x.

Beta:  (In) A scalar constant used to scale y.

Omega:  (In) A scalar constant used to scale w.

Zeta:  (In) A scalar constant used to scale z.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::MultiplyAndMatTransMultiply
"int Epetra_OskiMatrix::MultiplyAndMatTransMultiply(bool TransA,
const Epetra_MultiVector &X, Epetra_MultiVector &Y, const
Epetra_MultiVector &W, Epetra_MultiVector &Z, double Alpha=1.0, double
Beta=0.0, double Omega=1.0, double Zeta=0.0) const

Performs the two matrix multi-vector multiplies of Y = Alpha*this*X +
Beta*Y and Z = Omega*this^TransA*W + Zeta*Z.

The multi-vectors X, Y, W and Z can be either Epetra_MultiVectors or
Epetra_OskiMultiVectors. This composed routine is most commonly used
in bi-conjugate gradient calculations.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
computing the second product.

X:  (In) A multi-vector the matrix is multiplied by.

Y:  (In/Out) A multi-vector where the calculation result of the first
multiply is stored.

W:  (In) A multi-vector the matrix is multiplied by.

Z:  (In/Out) A multi-vector where the calculation result of the second
multiply is stored.

Alpha:  (In) A scalar constant used to scale X.

Beta:  (In) A scalar constant used to scale Y.

Omega:  (In) A scalar constant used to scale W.

Zeta:  (In) A scalar constant used to scale Z.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::MatPowMultiply "int
Epetra_OskiMatrix::MatPowMultiply(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y, Epetra_MultiVector &T, int Power=2, double
Alpha=1.0, double Beta=0.0) const

Performs a matrix vector multiply of y = Alpha*(this^TransA)^Power*x +
Beta*y. This is not implemented as described in the detailed
description.

The vectors x and y can be either Epetra_Vectors or
Epetra_OskiVectors. The vector T can be either an Epetra_MultiVector
or Epetra_OskiMultiVector. This composed routine is used in power and
S-step methods. This routine is not implemented due a bug in the
oski-1.01h kernel that makes testing of correctness impossible.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
computing the product.

x:  (In) The vector the matrix is multiplied by.

y:  (In/Out) The vector where the calculation result is stored.

T:  (Out) The multi-vector where the result of each subsequent
multiplication this*x ... this^(Power-1)*x is stored.

Power:  (In) The power to raise the matrix to in the calculation.

Alpha:  (In) A scalar constant used to scale x.

Beta:  (In) A scalar constant used to scale y.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

%feature("docstring")  Epetra_OskiMatrix::MatPowMultiply "int
Epetra_OskiMatrix::MatPowMultiply(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y, int Power=2, double Alpha=1.0, double Beta=0.0)
const

Performs a matrix vector multiply of y = Alpha*(this^TransA)^Power*x +
Beta*y. This is not implemented as described in the detailed
description.

The vectors x and y can be either Epetra_Vectors or
Epetra_OskiVectors. This composed routine is used in power and S-step
methods. This routine is not implemented due a bug in the oski-1.01h
kernel that makes testing of correctness impossible.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the matrix in
computing the product.

x:  (In) The vector the matrix is multiplied by.

y:  (In/Out) The vector where the calculation result is stored.

Power:  (In) The power to raise the matrix to in the calculation.

Alpha:  (In) A scalar constant used to scale x.

Beta:  (In) A scalar constant used to scale y.

Integer error code, set to 0 if successful.

Filled()==true

Unchanged ";

/*  Tuning  */

%feature("docstring")  Epetra_OskiMatrix::SetHint "int
Epetra_OskiMatrix::SetHint(const Teuchos::ParameterList &List)

Stores the hints in List in the matrix structure.

Parameters:
-----------

List:  (In) A list of hints and options to register along with the
matrix used for tuning purposes. The full list is given below. It may
be moved to the user guide in the future.

On successful storage of the hint 0 is returned. On failure an error
code is returned.  Options that can be passed to the List are
presented below. They are: \"<type> <option name>: <description of
purpose>\". The available hints are grouped by section, and only one
hint from each section can be true for a given matrix. For options
where multiple arguments can be passed in at once the interface only
supports up to 5. This means only 5 block sizes or 5 diaganols can be
passed in at once. If you have more changing the code to support your
needs should be simple, but every case you use must be enumerated. Of
course you can just say there are diagonals and blocks and not pass in
specific sizes as wells.

bool noblocks: If true, the matrix has no block structure

bool singleblocksize: If true, the matrix structure is dominated by
blocks of the size of the next two parameters. int row: The number of
rows in each block.

int col: The number of columns in each block.

bool multipleblocksize: If true, the matrix consists of multiple block
sizes. The next 3 parameters describe these and are optional. int
blocks: The number of block sizes in the matrix.

int row<x>: Where x is the block number, and x goes from 1 to blocks.
This is the number of rows in block x.

int col<x>: Where x is the block number, and x goes from 1 to blocks.
This is the number of cols in block x.

bool alignedblocks: If true, all blocks are aligned to a grid.

bool unalignedblocks: If true, blocks are not aligned to a grid.

bool symmetricpattern: If true, the matrix is either symmetric or
nearly symmetric.

bool nonsymmetricpattern: If true, the matrix has a very unsymmetric
pattern.

bool randompattern: If true, the matrix's non-zeros are distributed in
a random pattern.

bool correlatedpattern: If true, the row and column indices for non-
zeros are highly correlated.

bool nodiags : If true, the matrix has little or no diagonal
structure.

bool diags: If true, the matrix consists of diagonal structure
described the next two optional parameters. int numdiags: The number
of diagonal sizes known to be present others not listed could be
present.

int diag<x>: Where x is the diagonal number, and x goes from 1 to
numdiags. This is the size of the diagonal. ";

%feature("docstring")  Epetra_OskiMatrix::SetHintMultiply "int
Epetra_OskiMatrix::SetHintMultiply(bool TransA, double Alpha, const
Epetra_OskiMultiVector &InVec, double Beta, const
Epetra_OskiMultiVector &OutVec, int NumCalls, const
Teuchos::ParameterList &List)

Workload hints for computing a matrix-vector multiply used by
OskiTuneMat to optimize the data structure storage, and the routine to
compute the calculation.

In parallel the routine uses symbolic vectors. This is done for two
reasons. Doing this saves on data allocation and potentially
communication overhead. For a matrix-vector routine there should be no
advantage to having the actual vector, as its size must be the same as
a matrix dimension. For a matrix-multivector routine there could be
gains from knowing the number of vectors in the multi-vector. However,
OSKI does not perform multi-vector optimizations, so there is no need
to add the overhead.

Parameters:
-----------

Trans:  (In) If Trans = true then the transpose of the matrix will be
used in computing the product.

Alpha:  (In) A scalar constant used to scale InVec.

InVec:  (In) The vector the matrix is multiplied by or whether it is a
single vector or multi-vector.

Beta:  (In) A scalar constant used to scale OutVec.

OutVec:  (In) The vector where the calculation result is stored or
whether it is a single vector or multi-vector.

NumCalls:  (In) The number of times the operation is called or the
tuning level wanted.

List:  (In) Used for denoting the use of symbolic vectors for both
InVec and OutVec, as well as for level of aggressive tuning if either
NumCalls not known or to be overridden. Options are shown below. It
should be noted that by using these options the associated vector or
NumCalls becomes invalid.

Stores the workload hint in the matrix if the operation is valid. If
the operation is not valid an error code is returned.  Options that
can be passed to the List are presented below. They are: \"<type>
<option name>: <description of purpose>\". The available hints are
grouped by section, and only one hint from each section can be true
for a given matrix.

These replace InVec. bool symminvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symminmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace OutVec. bool symmoutvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symmoutmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

bool tune: If true, have OSKI tune moderately rather than using the
number of calls passed in.

bool tuneaggressive: If true, have OSKI tune aggressively rather than
using the number of calls passed in. ";

%feature("docstring")  Epetra_OskiMatrix::SetHintSolve "int
Epetra_OskiMatrix::SetHintSolve(bool TransA, double Alpha, const
Epetra_OskiMultiVector &Vector, int NumCalls, const
Teuchos::ParameterList &List)

Workload hints for computing a triangular solve used by OskiTuneMat to
optimize the data structure storage, and the routine to compute the
calculation.

In parallel the routine uses symbolic vectors. This is done for two
reasons. Doing this saves on data allocation and potentially
communication overhead. For a matrix-vector routine there should be no
advantage to having the actual vector, as its size must be the same as
a matrix dimension. For a matrix-multivector routine there could be
gains from knowing the number of vectors in the multi-vector. However,
OSKI does not perform multi-vector optimizations, so there is no need
to add the overhead.

Parameters:
-----------

Trans:  (In) If Trans = true then the transpose of the matrix will be
used in computing the product.

Alpha:  (In) A scalar constant used to scale InVec.

Vector:  (In) The vector being used in the solve and to store the
solution.

NumCalls:  (In) The number of times the operation is called or the
tuning level wanted.

List:  (In) Used for denoting the use of a symbolic vectors, as well
as for level of aggressive tuning if either NumCalls not known or to
be overridden. Options are shown below. It should be noted that by
using these options the associated vector or NumCalls becomes invalid.

Stores the workload hint in the matrix if the operation is valid. If
the operation is not valid an error code is returned.  Options that
can be passed to the List are presented below. They are: \"<type>
<option name>: <description of purpose>\". The available hints are
grouped by section, and only one hint from each section can be true
for a given matrix.

These replace Vector. bool symmvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symmmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

bool tune: If true, have OSKI tune moderately rather than using the
number of calls passed in.

bool tuneaggressive: If true, have OSKI tune aggressively rather than
using the number of calls passed in. ";

%feature("docstring")  Epetra_OskiMatrix::SetHintMatTransMatMultiply "int Epetra_OskiMatrix::SetHintMatTransMatMultiply(bool ATA, double
Alpha, const Epetra_OskiMultiVector &InVec, double Beta, const
Epetra_OskiMultiVector &OutVec, const Epetra_OskiMultiVector
&Intermediate, int NumCalls, const Teuchos::ParameterList &List)

Workload hints for computing a two matrix-vector multiplies that are
composed used by OskiTuneMat to optimize the data structure storage,
and the routine to compute the calculation.

In parallel the routine uses symbolic vectors. This is done for two
reasons. Doing this saves on data allocation and potentially
communication overhead. For a matrix-vector routine there should be no
advantage to having the actual vector, as its size must be the same as
a matrix dimension. For a matrix-multivector routine there could be
gains from knowing the number of vectors in the multi-vector. However,
OSKI does not perform multi-vector optimizations, so there is no need
to add the overhead.

Parameters:
-----------

ATA:  (In) If ATA = true then this^T*this*x will be computed otherwise
this*this^T*x will be.

Alpha:  (In) A scalar constant used to scale InVec.

InVec:  (In) The vector the matrix is multiplied by or whether it is a
single vector or multi-vector.

Beta:  (In) A scalar constant used to scale OutVec.

OutVec:  (In) The vector where the calculation result is stored or
whether it is a single vector or multi-vector.

Intermediate:  (In) The vector where result of the first product can
be stored or whether it is a single vector or multi-vector. If this
quantity is NULL then the intermediate product is not stored.

NumCalls:  (In) The number of times the operation is called or the
tuning level wanted.

List:  (In) Used for denoting the use of symbolic vectors for InVec,
OutVec and Intermediate, along with the level of aggressive tuning if
either NumCalls not known or to be overridden. Options are shown
below. It should be noted that by using these options the associated
vector or NumCalls becomes invalid.

Stores the workload hint in the matrix if the operation is valid. If
the operation is not valid an error code is returned.  Options that
can be passed to the List are presented below. They are: \"<type>
<option name>: <description of purpose>\". The available hints are
grouped by section, and only one hint from each section can be true
for a given matrix.

These replace InVec. bool symminvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symminmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace OutVec. bool symmoutvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symmoutmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace Intermediate. bool symmintervec: If true, use a symbolic
vector rather than the vector passed in for tuning purposes.

bool symmintermultivec: If true, use a symbolic multi-vector rather
than the multi-vector passed in for tuning purposes.

bool tune: If true, have OSKI tune moderately rather than using the
number of calls passed in.

bool tuneaggressive: If true, have OSKI tune aggressively rather than
using the number of calls passed in. ";

%feature("docstring")
Epetra_OskiMatrix::SetHintMultiplyAndMatTransMultiply "int
Epetra_OskiMatrix::SetHintMultiplyAndMatTransMultiply(bool TransA,
double Alpha, const Epetra_OskiMultiVector &InVec, double Beta, const
Epetra_OskiMultiVector &OutVec, double Omega, const
Epetra_OskiMultiVector &InVec2, double Zeta, const
Epetra_OskiMultiVector &OutVec2, int NumCalls, const
Teuchos::ParameterList &List)

Workload hints for computing two matrix-vector multiplies used by
OskiTuneMat to optimize the data structure storage, and the routine to
compute the calculation.

In parallel the routine uses symbolic vectors. This is done for two
reasons. Doing this saves on data allocation and potentially
communication overhead. For a matrix-vector routine there should be no
advantage to having the actual vector, as its size must be the same as
a matrix dimension. For a matrix-multivector routine there could be
gains from knowing the number of vectors in the multi-vector. However,
OSKI does not perform multi-vector optimizations, so there is no need
to add the overhead.

Parameters:
-----------

Trans:  (In) If Trans = true then the transpose of the matrix will be
used in computing the product.

Alpha:  (In) A scalar constant used to scale InVec.

InVec:  (In) The vector the matrix is multiplied by or whether it is a
single vector or multi-vector.

Beta:  (In) A scalar constant used to scale OutVec.

OutVec:  (In) The vector where the calculation result is stored or
whether it is a single vector or multi-vector.

Omega:  (In) A scalar constant used to scale InVec2.

InVec2:  (In) The vector the matrix is multiplied by or whether it is
a single vector or multi-vector.

Zeta:  (In) A scalar constant used to scale OutVec2.

OutVec2:  (In) The vector where the calculation result is stored or
whether it is a single vector or multi-vector.

NumCalls:  (In) The number of times the operation is called or the
tuning level wanted.

List:  (In) Used for denoting the use of symbolic vectors for both
InVec and OutVec, as well as for level of aggressive tuning if either
NumCalls not known or to be overridden. Options are shown below. It
should be noted that by using these options the associated vector or
NumCalls becomes invalid.

Stores the workload hint in the matrix if the operation is valid. If
the operation is not valid an error code is returned.  Options that
can be passed to the List are presented below. They are: \"<type>
<option name>: <description of purpose>\". The available hints are
grouped by section, and only one hint from each section can be true
for a given matrix.

These replace InVec. bool symminvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symminmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace OutVec. bool symmoutvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symmoutmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace InVec2. bool symminvec2: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symminmultivec2: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace OutVec2. bool symmoutvec2: If true, use a symbolic
vector rather than the vector passed in for tuning purposes.

bool symmoutmultivec2: If true, use a symbolic multi-vector rather
than the multi-vector passed in for tuning purposes.

bool tune: If true, have OSKI tune moderately rather than using the
number of calls passed in.

bool tuneaggressive: If true, have OSKI tune aggressively rather than
using the number of calls passed in. ";

%feature("docstring")  Epetra_OskiMatrix::SetHintPowMultiply "int
Epetra_OskiMatrix::SetHintPowMultiply(bool TransA, double Alpha, const
Epetra_OskiMultiVector &InVec, double Beta, const
Epetra_OskiMultiVector &OutVec, const Epetra_OskiMultiVector
&Intermediate, int Power, int NumCalls, const Teuchos::ParameterList
&List)

Workload hints for computing a matrix-vector multiply performed Power
times used by OskiTuneMat to optimize the data structure storage and
the routine to compute the calculation.

In parallel the routine uses symbolic vectors. This is done for two
reasons. Doing this saves on data allocation and potentially
communication overhead. For a matrix-vector routine there should be no
advantage to having the actual vector, as its size must be the same as
a matrix dimension. For a matrix-multivector routine there could be
gains from knowing the number of vectors in the multi-vector. However,
OSKI does not perform multi-vector optimizations, so there is no need
to add the overhead.

Parameters:
-----------

Trans:  (In) If Trans = true then the transpose of the matrix will be
used in computing the product.

Alpha:  (In) A scalar constant used to scale InVec.

InVec:  (In) The vector the matrix is multiplied by or whether it is a
single vector or multi-vector.

Beta:  (In) A scalar constant used to scale OutVec.

OutVec:  (In) The vector where the calculation result is stored or
whether it is a single vector or multi-vector.

Intermediate:  (In) The multi-vector where result of the first product
can be stored or whether it is a single vector or multi-vector. If
this quantity is NULL then the intermediate product is not stored.

Power:  (In) The power to raise the matrix to in the calculation.

NumCalls:  (In) The number of times the operation is called or the
tuning level wanted.

List:  (In) Used for denoting the use of symbolic vectors for both
InVec and OutVec, as well as for level of aggressive tuning if either
NumCalls not known or to be overridden. Options are shown below. It
should be noted that by using these options the associated vector or
NumCalls becomes invalid.

Stores the workload hint in the matrix if the operation is valid. If
the operation is not valid an error code is returned.  Options that
can be passed to the List are presented below. They are: \"<type>
<option name>: <description of purpose>\". The available hints are
grouped by section, and only one hint from each section can be true
for a given matrix.

These replace InVec. bool symminvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symminmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

These replace OutVec. bool symmoutvec: If true, use a symbolic vector
rather than the vector passed in for tuning purposes.

bool symmoutmultivec: If true, use a symbolic multi-vector rather than
the multi-vector passed in for tuning purposes.

This replaces Intermediate. bool symmintermultivec: If true, use a
symbolic multi-vector rather than the multi-vector passed in for
tuning purposes.

bool tune: If true, have OSKI tune moderately rather than using the
number of calls passed in.

bool tuneaggressive: If true, have OSKI tune aggressively rather than
using the number of calls passed in. ";

%feature("docstring")  Epetra_OskiMatrix::TuneMatrix "int
Epetra_OskiMatrix::TuneMatrix()

Tunes the matrix multiply if its deemed profitable.

The routine tunes based upon user provided hints if given. If hints
are not given the tuning is performed based on expected future
workload for the calculation. On success returns a non-negative status
code of the transformations performed. On failure an error code is
returned. ";

/*  Data Structure Transformation Methods  */

%feature("docstring")  Epetra_OskiMatrix::IsMatrixTransformed "int
Epetra_OskiMatrix::IsMatrixTransformed() const

Returns 1 if the matrix has been reordered by tuning, and 0 if it has
not been. ";

%feature("docstring")  Epetra_OskiMatrix::ViewTransformedMat "const
Epetra_OskiMatrix& Epetra_OskiMatrix::ViewTransformedMat() const

Returns the transformed version of InMat if InMat has been
transformed. If InMat has not been transformed then the return will
equal InMat. ";

%feature("docstring")  Epetra_OskiMatrix::ViewRowPermutation "const
Epetra_OskiPermutation& Epetra_OskiMatrix::ViewRowPermutation() const

Returns a read only row/left permutation of the Matrix. ";

%feature("docstring")  Epetra_OskiMatrix::ViewColumnPermutation "const Epetra_OskiPermutation&
Epetra_OskiMatrix::ViewColumnPermutation() const

Returns a read only column/right permutation of the Matrix. ";

%feature("docstring")  Epetra_OskiMatrix::GetMatrixTransforms "char*
Epetra_OskiMatrix::GetMatrixTransforms() const

Returns a string holding the transformations performed on the matrix
when it was tuned.

Upon success returns a newly-allocated string that stores the
transformations applied to the matrix during tuning. NULL is returned
upon an error. It is the users responsibility to deallocate the
returned string. ";

%feature("docstring")  Epetra_OskiMatrix::ApplyMatrixTransforms "int
Epetra_OskiMatrix::ApplyMatrixTransforms(const char *Transforms)

Replaces the current data structure of the matrix with the one
specified in Transforms.

If a previously tuned copy of the Matrix existed it is now replaced by
one specified in Transforms.

Parameters:
-----------

Transforms:  (In) A string that holds the transformations to be
applied to the matrix. If Transforms is NULL or the empty string then
no changes to the data structure are made.

If the transformation was successful 0 is returned. Otherwise an error
code is returned. ";


// File: classEpetra__OskiMultiVector.xml
%feature("docstring") Epetra_OskiMultiVector "

Epetra_OskiMultiVector: A class for constructing and using dense Oski
multi-vectors on a single processor or a single core of a multi-
processor.

The Epetra_OskiMultiVector class enables the construction and use of
real-valued, double- precision dense vectors and multi-vectors, in a
serial environment. The dimensions of the dense multi-vectors comes
from the inherited Epetra_MultiVector object. All values and data
layouts are kept the same and the multi- vector is wrapped for use
with OSKI.

C++ includes: Epetra_OskiMultiVector.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_OskiMultiVector::Epetra_OskiMultiVector
"Epetra_OskiMultiVector::Epetra_OskiMultiVector(const
Epetra_OskiMultiVector &Source)

Copy constructor. ";

%feature("docstring")  Epetra_OskiMultiVector::Epetra_OskiMultiVector
"Epetra_OskiMultiVector::Epetra_OskiMultiVector(const
Epetra_MultiVector &Source)

Constructor creates and Epetra_OskiMultiVector from an
Epetra_MultiVector.

Parameters:
-----------

Source:  (In) An Epetra_MultiVector that is wrapped as an
Epetra_OskiMultiVector.

Pointer to an Epetra_OskiMultiVector.

If the Epetra_MultiVector is not stored contigously according to the
BLAS standard then a deep copy is made. ";

%feature("docstring")  Epetra_OskiMultiVector::~Epetra_OskiMultiVector
"virtual Epetra_OskiMultiVector::~Epetra_OskiMultiVector()

Destructor. ";

/*  Extraction Methods  */

%feature("docstring")  Epetra_OskiMultiVector::Copy_Created "bool
Epetra_OskiMultiVector::Copy_Created() const

Returns true if a deep copy of the multi-vector was created by the
constructor. ";

%feature("docstring")  Epetra_OskiMultiVector::Oski_View "oski_vecview_t Epetra_OskiMultiVector::Oski_View() const

Returns the Oski portion of the Multi-Vector. ";

%feature("docstring")  Epetra_OskiMultiVector::Epetra_View "const
Epetra_MultiVector* Epetra_OskiMultiVector::Epetra_View() const

Returns the Epetra portion of the Multi-Vector. ";

/*  Operators  */


// File: classEpetra__OskiPermutation.xml
%feature("docstring") Epetra_OskiPermutation "

Epetra_OskiPermutation: A class for storing the permutation performed
on a Epetra_OskiMatrix.

The Epetra_OskiPermutation is one of the possible transformations that
OSKI can perform on a matrix. The permutation is stored with the
matrix in OSKI. Using this class a Epetra_OskiPermutation can be
applied to an Epetra_OskiMultiVector.

C++ includes: Epetra_OskiPermutation.h ";

/*  Constructor/Destructor  */

%feature("docstring")  Epetra_OskiPermutation::Epetra_OskiPermutation
"Epetra_OskiPermutation::Epetra_OskiPermutation()

Default Constructor. ";

%feature("docstring")  Epetra_OskiPermutation::Epetra_OskiPermutation
"Epetra_OskiPermutation::Epetra_OskiPermutation(const
Epetra_OskiPermutation &Source)

Copy Constructor. ";

%feature("docstring")  Epetra_OskiPermutation::Epetra_OskiPermutation
"Epetra_OskiPermutation::Epetra_OskiPermutation(bool RowPerm, const
Epetra_OskiMatrix &Source)

Constructor creates an Epetra_OskiPermutation from an
Epetra_OskiMatrix.

Acquires the permutation from the passed in matrix and stores it
within the object. If RowPerm is true this is a row permutation and if
RowPerm is false this is a column permutation. ";

%feature("docstring")  Epetra_OskiPermutation::~Epetra_OskiPermutation
"virtual Epetra_OskiPermutation::~Epetra_OskiPermutation()

Destructor. ";

/*  Replace Method  */

%feature("docstring")  Epetra_OskiPermutation::ReplacePermutation "void Epetra_OskiPermutation::ReplacePermutation(const oski_perm_t
&InPerm)

Stores a permutation in the data structure. ";

/*  Apply  */

%feature("docstring")  Epetra_OskiPermutation::PermuteVector "int
Epetra_OskiPermutation::PermuteVector(const bool TransA,
Epetra_OskiMultiVector &Vector) const

Permutes Vector according to the Permutation. If a transpose is
desired it performs that operation.

The Vector passed into this function is a view. It is the underlying
object that is permuted by the function.

Parameters:
-----------

TransA:  (In) If TransA = TRUE then use the transpose of the
permutation and apply it to Vector.

Vector:  (In/Out) The vector that is permuted by Permutation^Trans.

When successful returns 0. On error Vector is not permuted and a error
code is returned. ";


// File: classEpetra__OskiUtils.xml
%feature("docstring") Epetra_OskiUtils "

Epetra_OskiUtils: The Epetra OSKI Class to handle all operations that
do not involve the use of a matrix, vector, error or permutation
object.

The Epetra_OskiUtils class is a helper class used to call OSKI
functions that do not use matrix, vector, error or permutation
objects. It provides an interface to access the initialization and
finalize routines of OSKI.

All functions are public to allow access to methods needed by programs
using OSKI. There are no data members of the class as all data is kept
in the matrix, vector, multi-vector, error and permutation classes.

C++ includes: Epetra_OskiUtils.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_OskiUtils::Epetra_OskiUtils "Epetra_OskiUtils::Epetra_OskiUtils()

Default Constructor. ";

%feature("docstring")  Epetra_OskiUtils::~Epetra_OskiUtils "virtual
Epetra_OskiUtils::~Epetra_OskiUtils()

Destructor. ";

/*  Start/End  */

%feature("docstring")  Epetra_OskiUtils::Init "void
Epetra_OskiUtils::Init()

Initializes OSKI.

Calls the OSKI routine to initialize the use of OSKI. This routine is
required before OSKI can be used. ";

%feature("docstring")  Epetra_OskiUtils::Close "void
Epetra_OskiUtils::Close()

Finalizes the use of OSKI.

When done using OSKI this routine performs cleanup operations. While
not strictly required it is highly recommended to be called when OSKI
is no longer being used. ";


// File: classEpetra__OskiVector.xml
%feature("docstring") Epetra_OskiVector "

Epetra_OskiVector: A class for constructing and using dense OSKI
vectors on a single processor or a single core of a multi-processor.

The Epetra_OskiVector class enables the construction and use of real-
valued, double- precision dense vectors in a serial environment. The
dimensions of the dense vectors comes from the inherited Epetra_Vector
object. All values and data layouts are kept the same and the vector
is wrapped for use with OSKI.

C++ includes: Epetra_OskiVector.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_OskiVector::Epetra_OskiVector "Epetra_OskiVector::Epetra_OskiVector(const Epetra_OskiVector &Source)

Copy constructor. ";

%feature("docstring")  Epetra_OskiVector::Epetra_OskiVector "Epetra_OskiVector::Epetra_OskiVector(const Epetra_Vector &Source)

Constructor creates and Epetra_OskiVector from an Epetra_Vector.

Parameters:
-----------

Source:  (In) An Epetra_Vector that is to be wrapped as an
Epetra_OskiVector.

Pointer to an Epetra_OskiVector. ";

%feature("docstring")  Epetra_OskiVector::~Epetra_OskiVector "virtual
Epetra_OskiVector::~Epetra_OskiVector()

Destructor. ";

/*  Extraction Method  */

%feature("docstring")  Epetra_OskiVector::Epetra_View "const
Epetra_Vector* Epetra_OskiVector::Epetra_View() const

Returns a view to the Epetra Object. ";

/*  Operators  */


// File: classEpetra__RowMatrix.xml
%feature("docstring") Epetra_RowMatrix "

Epetra_RowMatrix: A pure virtual class for using real-valued double-
precision row matrices.

The Epetra_RowMatrix class is a pure virtual class (specifies
interface only) that enable the use of real-valued double-precision
sparse matrices where matrix entries are intended for row access. It
is currently implemented by both the Epetra_CrsMatrix and
Epetra_VbrMatrix classes.

C++ includes: Epetra_RowMatrix.h ";

/*  Destructor  */

%feature("docstring")  Epetra_RowMatrix::~Epetra_RowMatrix "virtual
Epetra_RowMatrix::~Epetra_RowMatrix()

Destructor. ";

/*  Matrix data extraction routines  */

%feature("docstring")  Epetra_RowMatrix::NumMyRowEntries "virtual int
Epetra_RowMatrix::NumMyRowEntries(int MyRow, int &NumEntries) const =0

Returns the number of nonzero entries in MyRow.

Parameters:
-----------

In:  MyRow - Local row.

Out:  NumEntries - Number of nonzero values present.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::MaxNumEntries "virtual int
Epetra_RowMatrix::MaxNumEntries() const =0

Returns the maximum of NumMyRowEntries() over all rows. ";

%feature("docstring")  Epetra_RowMatrix::ExtractMyRowCopy "virtual
int Epetra_RowMatrix::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const =0

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

In:  MyRow - Local row to extract.

In:  Length - Length of Values and Indices.

Out:  NumEntries - Number of nonzero entries extracted.

Out:  Values - Extracted values for this row.

Out:  Indices - Extracted local column indices for the corresponding
values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::ExtractDiagonalCopy "virtual
int Epetra_RowMatrix::ExtractDiagonalCopy(Epetra_Vector &Diagonal)
const =0

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Out:  Diagonal - Extracted main diagonal.

Integer error code, set to 0 if successful. ";

/*  Mathematical functions  */

%feature("docstring")  Epetra_RowMatrix::Multiply "virtual int
Epetra_RowMatrix::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const =0

Returns the result of a Epetra_RowMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  TransA -If true, multiply by the transpose of matrix, otherwise
just use matrix.

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::Solve "virtual int
Epetra_RowMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const =0

Returns result of a local-only solve using a triangular
Epetra_RowMatrix with Epetra_MultiVectors X and Y.

This method will perform a triangular solve independently on each
processor of the parallel machine. No communication is performed.

Parameters:
-----------

In:  Upper -If true, solve Ux = y, otherwise solve Lx = y.

In:  Trans -If true, solve transpose problem.

In:  UnitDiagonal -If true, assume diagonal is unit (whether it's
stored or not).

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::InvRowSums "virtual int
Epetra_RowMatrix::InvRowSums(Epetra_Vector &x) const =0

Computes the sum of absolute values of the rows of the
Epetra_RowMatrix, results returned in x.

The vector x will return such that x[i] will contain the inverse of
sum of the absolute values of the this matrix will be scaled such that
A(i,j) = x(i)*A(i,j) where i denotes the global row number of A and j
denotes the global column number of A. Using the resulting vector from
this function as input to LeftScale() will make the infinity norm of
the resulting matrix exactly 1.

Parameters:
-----------

Out:  x -A Epetra_Vector containing the row sums of the this matrix.

WARNING:  It is assumed that the distribution of x is the same as the
rows of this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::LeftScale "virtual int
Epetra_RowMatrix::LeftScale(const Epetra_Vector &x)=0

Scales the Epetra_RowMatrix on the left with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i
denotes the row number of A and j denotes the column number of A.

Parameters:
-----------

In:  x -A Epetra_Vector to solve for.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::InvColSums "virtual int
Epetra_RowMatrix::InvColSums(Epetra_Vector &x) const =0

Computes the sum of absolute values of the columns of the
Epetra_RowMatrix, results returned in x.

The vector x will return such that x[j] will contain the inverse of
sum of the absolute values of the this matrix will be sca such that
A(i,j) = x(j)*A(i,j) where i denotes the global row number of A and j
denotes the global column number of A. Using the resulting vector from
this function as input to RighttScale() will make the one norm of the
resulting matrix exactly 1.

Parameters:
-----------

Out:  x -A Epetra_Vector containing the column sums of the this
matrix.

WARNING:  It is assumed that the distribution of x is the same as the
rows of this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_RowMatrix::RightScale "virtual int
Epetra_RowMatrix::RightScale(const Epetra_Vector &x)=0

Scales the Epetra_RowMatrix on the right with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i
denotes the global row number of A and j denotes the global column
number of A.

Parameters:
-----------

In:  x -The Epetra_Vector used for scaling this.

Integer error code, set to 0 if successful. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_RowMatrix::Filled "virtual bool
Epetra_RowMatrix::Filled() const =0

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

%feature("docstring")  Epetra_RowMatrix::NormInf "virtual double
Epetra_RowMatrix::NormInf() const =0

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Epetra_RowMatrix::NormOne "virtual double
Epetra_RowMatrix::NormOne() const =0

Returns the one norm of the global matrix. ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalNonzeros "virtual
int Epetra_RowMatrix::NumGlobalNonzeros() const =0

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalNonzeros64 "virtual
long long Epetra_RowMatrix::NumGlobalNonzeros64() const =0 ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalRows "virtual int
Epetra_RowMatrix::NumGlobalRows() const =0

Returns the number of global matrix rows. ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalRows64 "virtual
long long Epetra_RowMatrix::NumGlobalRows64() const =0 ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalCols "virtual int
Epetra_RowMatrix::NumGlobalCols() const =0

Returns the number of global matrix columns. ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalCols64 "virtual
long long Epetra_RowMatrix::NumGlobalCols64() const =0 ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalDiagonals "virtual
int Epetra_RowMatrix::NumGlobalDiagonals() const =0

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Epetra_RowMatrix::NumGlobalDiagonals64 "virtual long long Epetra_RowMatrix::NumGlobalDiagonals64() const =0 ";

%feature("docstring")  Epetra_RowMatrix::NumMyNonzeros "virtual int
Epetra_RowMatrix::NumMyNonzeros() const =0

Returns the number of nonzero entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Epetra_RowMatrix::NumMyRows "virtual int
Epetra_RowMatrix::NumMyRows() const =0

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Epetra_RowMatrix::NumMyCols "virtual int
Epetra_RowMatrix::NumMyCols() const =0

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")  Epetra_RowMatrix::NumMyDiagonals "virtual int
Epetra_RowMatrix::NumMyDiagonals() const =0

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Epetra_RowMatrix::LowerTriangular "virtual
bool Epetra_RowMatrix::LowerTriangular() const =0

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_RowMatrix::UpperTriangular "virtual
bool Epetra_RowMatrix::UpperTriangular() const =0

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_RowMatrix::RowMatrixRowMap "virtual
const Epetra_Map& Epetra_RowMatrix::RowMatrixRowMap() const =0

Returns the Epetra_Map object associated with the rows of this matrix.
";

%feature("docstring")  Epetra_RowMatrix::RowMatrixColMap "virtual
const Epetra_Map& Epetra_RowMatrix::RowMatrixColMap() const =0

Returns the Epetra_Map object associated with the columns of this
matrix. ";

%feature("docstring")  Epetra_RowMatrix::RowMatrixImporter "virtual
const Epetra_Import* Epetra_RowMatrix::RowMatrixImporter() const =0

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";


// File: classEpetra__RowMatrixTransposer.xml
%feature("docstring") Epetra_RowMatrixTransposer "

Epetra_RowMatrixTransposer: A class for transposing an
Epetra_RowMatrix object.

This class provides capabilities to construct a transpose matrix of an
existing Epetra_RowMatrix object and (optionally) redistribute it
across a parallel distributed memory machine.

C++ includes: Epetra_RowMatrixTransposer.h ";

/*  Constructors/destructors  */

%feature("docstring")
Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer "Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer(Epetra_RowMatrix
*OrigMatrix)

Primary Epetra_RowMatrixTransposer constructor.

Parameters:
-----------

Matrix:  (In) An existing Epetra_RowMatrix object. The
Epetra_RowMatrix, the LHS and RHS pointers do not need to be defined
before this constructor is called.

Pointer to a Epetra_RowMatrixTransposer object. ";

%feature("docstring")
Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer "Epetra_RowMatrixTransposer::Epetra_RowMatrixTransposer(const
Epetra_RowMatrixTransposer &Source)

Epetra_RowMatrixTransposer copy constructor. ";

%feature("docstring")
Epetra_RowMatrixTransposer::~Epetra_RowMatrixTransposer "Epetra_RowMatrixTransposer::~Epetra_RowMatrixTransposer()

Epetra_RowMatrixTransposer destructor. ";

/*  Forward transformation methods  */

%feature("docstring")  Epetra_RowMatrixTransposer::CreateTranspose "int Epetra_RowMatrixTransposer::CreateTranspose(const bool
MakeDataContiguous, Epetra_CrsMatrix *&TransposeMatrix, Epetra_Map
*TransposeRowMap=0)

Generate a new Epetra_CrsMatrix as the transpose of an
Epetra_RowMatrix passed into the constructor.

Constructs a new Epetra_CrsMatrix that is a copy of the
Epetra_RowMatrix passed in to the constructor.

Parameters:
-----------

MakeDataContiguous:  (In) Causes the output matrix, LHS and RHS to be
stored in a form compatible with Fortran-style solvers. The output
matrix will be compatible with the Harwell-Boeing compressed column
format. The RHS and LHS will be stored such that the last value in
column j of the multivector is stored next to the first value in
column j+1.

TransposeRowMap:  (Optional/In) If this argument is defined, the
transpose matrix will be distributed using this map as the row map for
the transpose. If it is set to zero, the transpose matrix will use the
OrigMatrix->RowMatrixDomainMap as the row map.

Integer error code, 0 if no errors. Negative if some fatal error
occured. ";

%feature("docstring")
Epetra_RowMatrixTransposer::UpdateTransposeValues "int
Epetra_RowMatrixTransposer::UpdateTransposeValues(Epetra_RowMatrix
*MatrixWithNewValues)

Update the values of an already-redistributed problem.

Updates the values of an already-redistributed problem. This method
allows updating the redistributed problem without allocating new
storage.

Parameters:
-----------

MatrixWithNewValues:  (In) The values from MatrixWithNewValues will be
copied into the TransposeMatrix. The MatrixWithNewValues object must
be identical in structure to the original matrix object used to create
this instance of Epetra_RowMatrixTransposer.

Integer error code, 0 if no errors. Negative if some fatal error
occured. ";

/*  Reverse transformation methods  */

%feature("docstring")
Epetra_RowMatrixTransposer::UpdateOriginalMatrixValues "int
Epetra_RowMatrixTransposer::UpdateOriginalMatrixValues()

Update values of original matrix (Not implemented and not sure if we
will implement this). ";

/*  Attribute accessor methods  */

%feature("docstring")  Epetra_RowMatrixTransposer::TransposeRowMap "const Epetra_Map& Epetra_RowMatrixTransposer::TransposeRowMap() const

Returns const reference to the Epetra_Map object describing the row
distribution of the transpose matrix.

The RedistExporter object can be used to redistribute other
Epetra_DistObject objects whose maps are compatible with the original
linear problem map, or with the RedistMap(). WARNING:  Must not be
called before CreateTranspose()is called. ";

%feature("docstring")  Epetra_RowMatrixTransposer::TransposeExporter "const Epetra_Export& Epetra_RowMatrixTransposer::TransposeExporter()
const

Returns const reference to the Epetra_Export object used to
redistribute the original matrix.

The TransposeExporter object can be used to redistribute other
Epetra_DistObject objects whose maps are compatible with the original
matrix. WARNING:  Must not be called before CreateTranspose() is
called. ";


// File: classEpetra__SerialComm.xml
%feature("docstring") Epetra_SerialComm "

Epetra_SerialComm: The Epetra Serial Communication Class.

The Epetra_SerialComm class is an implementation of Epetra_Comm,
providing the general information and services needed for other Epetra
classes to run on a serial computer.

C++ includes: Epetra_SerialComm.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")  Epetra_SerialComm::Epetra_SerialComm "Epetra_SerialComm::Epetra_SerialComm()

Epetra_SerialComm Serial Constructor.

Builds an instance of a serial communicator. Even if the application
is running in parallel via MPI, this communicator will execute in
serial. The access functions return the number of processors to be 1
and the processor ID to be 0. ";

%feature("docstring")  Epetra_SerialComm::Epetra_SerialComm "Epetra_SerialComm::Epetra_SerialComm(const Epetra_SerialComm &Comm)

Epetra_SerialComm Copy Constructor.

Makes an exact copy of an existing Epetra_SerialComm instance. ";

%feature("docstring")  Epetra_SerialComm::Clone "Epetra_Comm*
Epetra_SerialComm::Clone() const

Clone method. ";

%feature("docstring")  Epetra_SerialComm::~Epetra_SerialComm "Epetra_SerialComm::~Epetra_SerialComm()

Epetra_SerialComm Destructor.

Completely deletes a Epetra_SerialComm object. WARNING:  Note: All
objects that depend on a Epetra_SerialComm instance should be
destroyed prior to calling this function. ";

/*  Barrier Methods  */

%feature("docstring")  Epetra_SerialComm::Barrier "void
Epetra_SerialComm::Barrier() const

Epetra_SerialComm Barrier function.

A no-op for a serial communicator. ";

/*  Broadcast Methods  */

%feature("docstring")  Epetra_SerialComm::Broadcast "int
Epetra_SerialComm::Broadcast(double *MyVals, int Count, int Root)
const

Epetra_SerialComm Broadcast function.

A no-op for a serial communicator.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of MyVals.

Root:  In On entry, contains the processor from which all processors
will receive a copy of MyVals. ";

%feature("docstring")  Epetra_SerialComm::Broadcast "int
Epetra_SerialComm::Broadcast(int *MyVals, int Count, int Root) const

Epetra_SerialComm Broadcast function.

A no-op for a serial communicator.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of MyVals.

Root:  In On entry, contains the processor from which all processors
will receive a copy of MyVals. ";

%feature("docstring")  Epetra_SerialComm::Broadcast "int
Epetra_SerialComm::Broadcast(long *MyVals, int Count, int Root) const

Epetra_SerialComm Broadcast function.

A no-op for a serial communicator.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of MyVals.

Root:  In On entry, contains the processor from which all processors
will receive a copy of MyVals. ";

%feature("docstring")  Epetra_SerialComm::Broadcast "int
Epetra_SerialComm::Broadcast(char *MyVals, int Count, int Root) const

Epetra_SerialComm Broadcast function.

A no-op for a serial communicator.

Parameters:
-----------

MyVals:  InOut On entry, the root processor contains the list of
values. On exit, all processors will have the same list of values.
Note that values must be allocated on all processor before the
broadcast.

Count:  In On entry, contains the length of the list of MyVals.

Root:  In On entry, contains the processor from which all processors
will receive a copy of MyVals. ";

/*  Gather Methods  */

%feature("docstring")  Epetra_SerialComm::GatherAll "int
Epetra_SerialComm::GatherAll(double *MyVals, double *AllVals, int
Count) const

Epetra_SerialComm All Gather function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_SerialComm::GatherAll "int
Epetra_SerialComm::GatherAll(int *MyVals, int *AllVals, int Count)
const

Epetra_SerialComm All Gather function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_SerialComm::GatherAll "int
Epetra_SerialComm::GatherAll(long *MyVals, long *AllVals, int Count)
const

Epetra_SerialComm All Gather function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

%feature("docstring")  Epetra_SerialComm::GatherAll "int
Epetra_SerialComm::GatherAll(long long *MyVals, long long *AllVals,
int Count) const

Epetra_SerialComm All Gather function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values, to be sent to all
processors.

AllVals:  Out On exit, contains the list of values from all
processors. Must by of size NumProc*Count.

Count:  In On entry, contains the length of the list of MyVals. ";

/*  Sum Methods  */

%feature("docstring")  Epetra_SerialComm::SumAll "int
Epetra_SerialComm::SumAll(double *PartialSums, double *GlobalSums, int
Count) const

Epetra_SerialComm Global Sum function.

A copy for a serial communicator.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::SumAll "int
Epetra_SerialComm::SumAll(int *PartialSums, int *GlobalSums, int
Count) const

Epetra_SerialComm Global Sum function.

A copy for a serial communicator.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::SumAll "int
Epetra_SerialComm::SumAll(long *PartialSums, long *GlobalSums, int
Count) const

Epetra_SerialComm Global Sum function.

A copy for a serial communicator.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::SumAll "int
Epetra_SerialComm::SumAll(long long *PartialSums, long long
*GlobalSums, int Count) const

Epetra_SerialComm Global Sum function.

A copy for a serial communicator.

Parameters:
-----------

PartialSums:  In On entry, contains the list of values, usually
partial sums computed locally, to be summed across all processors.

GlobalSums:  Out On exit, contains the list of values summed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Max/Min Methods  */

%feature("docstring")  Epetra_SerialComm::MaxAll "int
Epetra_SerialComm::MaxAll(double *PartialMaxs, double *GlobalMaxs, int
Count) const

Epetra_SerialComm Global Max function.

A copy for a serial communicator.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally, using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MaxAll "int
Epetra_SerialComm::MaxAll(int *PartialMaxs, int *GlobalMaxs, int
Count) const

Epetra_SerialComm Global Max function.

A copy for a serial communicator.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MaxAll "int
Epetra_SerialComm::MaxAll(long *PartialMaxs, long *GlobalMaxs, int
Count) const

Epetra_SerialComm Global Max function.

A copy for a serial communicator.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MaxAll "int
Epetra_SerialComm::MaxAll(long long *PartialMaxs, long long
*GlobalMaxs, int Count) const

Epetra_SerialComm Global Max function.

A copy for a serial communicator.

Parameters:
-----------

PartialMaxs:  In On entry, contains the list of values, usually
partial maxs computed locally; using these Partial Maxs, the max
across all processors will be computed.

GlobalMaxs:  Out On exit, contains the list of maxs computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MinAll "int
Epetra_SerialComm::MinAll(double *PartialMins, double *GlobalMins, int
Count) const

Epetra_SerialComm Global Min function.

A copy for a serial communicator.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MinAll "int
Epetra_SerialComm::MinAll(int *PartialMins, int *GlobalMins, int
Count) const

Epetra_SerialComm Global Min function.

A copy for a serial communicator.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MinAll "int
Epetra_SerialComm::MinAll(long *PartialMins, long *GlobalMins, int
Count) const

Epetra_SerialComm Global Min function.

A copy for a serial communicator.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::MinAll "int
Epetra_SerialComm::MinAll(long long *PartialMins, long long
*GlobalMins, int Count) const

Epetra_SerialComm Global Min function.

A copy for a serial communicator.

Parameters:
-----------

PartialMins:  In On entry, contains the list of values, usually
partial mins computed locally; using these Partial Mins, the min
across all processors will be computed.

GlobalMins:  Out On exit, contains the list of mins computed across
all processors.

Count:  In On entry, contains the length of the list of values. ";

/*  Parallel Prefix Methods  */

%feature("docstring")  Epetra_SerialComm::ScanSum "int
Epetra_SerialComm::ScanSum(double *MyVals, double *ScanSums, int
Count) const

Epetra_SerialComm Scan Sum function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::ScanSum "int
Epetra_SerialComm::ScanSum(int *MyVals, int *ScanSums, int Count)
const

Epetra_SerialComm Scan Sum function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::ScanSum "int
Epetra_SerialComm::ScanSum(long *MyVals, long *ScanSums, int Count)
const

Epetra_SerialComm Scan Sum function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

%feature("docstring")  Epetra_SerialComm::ScanSum "int
Epetra_SerialComm::ScanSum(long long *MyVals, long long *ScanSums, int
Count) const

Epetra_SerialComm Scan Sum function.

A copy for a serial communicator.

Parameters:
-----------

MyVals:  In On entry, contains the list of values to be summed across
all processors.

ScanSums:  Out On exit, contains the list of values summed across
processors 0 through i.

Count:  In On entry, contains the length of the list of values. ";

/*  Attribute Accessor Methods  */

%feature("docstring")  Epetra_SerialComm::MyPID "int
Epetra_SerialComm::MyPID() const

Return my process ID.

In MPI mode returns the rank of the calling process. In serial mode
returns 0. ";

%feature("docstring")  Epetra_SerialComm::NumProc "int
Epetra_SerialComm::NumProc() const

Returns total number of processes (always returns 1 for SerialComm).
";

/*  Gather/Scatter and Directory Constructors  */

%feature("docstring")  Epetra_SerialComm::CreateDistributor "Epetra_Distributor * Epetra_SerialComm::CreateDistributor() const

Create a distributor object. ";

%feature("docstring")  Epetra_SerialComm::CreateDirectory "Epetra_Directory * Epetra_SerialComm::CreateDirectory(const
Epetra_BlockMap &Map) const

Create a directory object for the given Epetra_BlockMap. ";

/*  Print object to an output stream  */

%feature("docstring")  Epetra_SerialComm::Print "void
Epetra_SerialComm::Print(ostream &os) const

Print method that implements Epetra_Object virtual Print method. ";

%feature("docstring")  Epetra_SerialComm::PrintInfo "void
Epetra_SerialComm::PrintInfo(ostream &os) const

Print method that implements Epetra_Comm virtual PrintInfo method. ";

/*  Expert Users and Developers Only  */

%feature("docstring")  Epetra_SerialComm::ReferenceCount "int
Epetra_SerialComm::ReferenceCount() const

Returns the reference count of SerialCommData.

(Intended for testing purposes.) ";

%feature("docstring")  Epetra_SerialComm::DataPtr "const
Epetra_SerialCommData* Epetra_SerialComm::DataPtr() const

Returns a pointer to the SerialCommData instance this SerialComm uses.

(Intended for developer use only for testing purposes.) ";


// File: classEpetra__SerialCommData.xml
%feature("docstring") Epetra_SerialCommData "

Epetra_SerialCommData: The Epetra Serial Communication Data Class.

The Epetra_SerialCommData class is an implementation detail of
Epetra_SerialComm. It is reference-counted, and can be shared by
multiple Epetra_SerialComm instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_SerialCommData.h ";

/*  Constructor/Destructor Methods  */


// File: classEpetra__SerialDenseMatrix.xml
%feature("docstring") Epetra_SerialDenseMatrix "

Epetra_SerialDenseMatrix: A class for constructing and using real
double precision general dense matrices.

The Epetra_SerialDenseMatrix class enables the construction and use of
real-valued, general, double-precision dense matrices. It is built on
the BLAS, and derives from the Epetra_BLAS.

The Epetra_SerialDenseMatrix class is intended to provide very basic
support for dense rectangular matrices.

Constructing Epetra_SerialDenseMatrix Objects

There are four Epetra_SerialDenseMatrix constructors. The first
constructs a zero-sized object which should be made to appropriate
length using the Shape() or Reshape() functions and then filled with
the [] or () operators. The second constructs an object sized to the
dimensions specified, which should be filled with the [] or ()
operators. The third is a constructor that accepts user data as a 2D
array, and the fourth is a copy constructor. The third constructor has
two data access modes (specified by the Epetra_DataAccess argument):
Copy mode - Allocates memory and makes a copy of the user-provided
data. In this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Extracting Data from Epetra_SerialDenseMatrix
Objects

Once a Epetra_SerialDenseMatrix is constructed, it is possible to view
the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.  Vector and Utility Functions

Once a Epetra_SerialDenseMatrix is constructed, several mathematical
functions can be applied to the object. Specifically: Multiplication.

Norms.

Counting floating point operations The Epetra_SerialDenseMatrix class
has Epetra_CompObject as a base class. Thus, floating point operations
are counted and accumulated in the Epetra_Flop object (if any) that
was set using the SetFlopCounter() method in the Epetra_CompObject
base class.

C++ includes: Epetra_SerialDenseMatrix.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix "Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(bool
set_object_label=true)

Default constructor; defines a zero size object.

Epetra_SerialDenseMatrix objects defined by the default constructor
should be sized with the Shape() or Reshape functions. Values should
be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix "Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(int NumRows, int
NumCols, bool set_object_label=true)

Shaped constructor; defines a variable-sized object.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Epetra_SerialDenseMatrix objects defined by the shaped constructor are
already shaped to the dimensions given as a parameters. All values are
initialized to 0. Calling this constructor is equivalent to using the
default constructor, and then calling the Shape function on it. Values
should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix "Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(Epetra_DataAccess
CV, double *A_in, int LDA_in, int NumRows, int NumCols, bool
set_object_label=true)

Set object values from two-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  A - Pointer to an array of double precision numbers. The first
vector starts at A. The second vector starts at A+LDA, the third at
A+2*LDA, and so on.

In:  LDA - The \"Leading Dimension\", or stride between vectors in
memory.

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

See Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix "Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(const
Epetra_SerialDenseMatrix &Source)

Epetra_SerialDenseMatrix copy constructor. ";

%feature("docstring")
Epetra_SerialDenseMatrix::~Epetra_SerialDenseMatrix "Epetra_SerialDenseMatrix::~Epetra_SerialDenseMatrix()

Epetra_SerialDenseMatrix destructor. ";

/*  Shaping/sizing Methods  */

%feature("docstring")  Epetra_SerialDenseMatrix::Shape "int
Epetra_SerialDenseMatrix::Shape(int NumRows, int NumCols)

Set dimensions of a Epetra_SerialDenseMatrix object; init values to
zero.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Allows user to define the dimensions of a Epetra_SerialDenseMatrix at
any point. This function can be called at any point after
construction. Any values that were previously in this object are
destroyed and the resized matrix starts off with all zero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::Reshape "int
Epetra_SerialDenseMatrix::Reshape(int NumRows, int NumCols)

Reshape a Epetra_SerialDenseMatrix object.

Parameters:
-----------

In:  NumRows - Number of rows in object.

In:  NumCols - Number of columns in object.

Allows user to define the dimensions of a Epetra_SerialDenseMatrix at
any point. This function can be called at any point after
construction. Any values that were previously in this object are
copied into the new shape. If the new shape is smaller than the
original, the upper left portion of the original matrix (the principal
submatrix) is copied to the new matrix.

Integer error code, set to 0 if successful. ";

/*  Mathematical methods  */

%feature("docstring")  Epetra_SerialDenseMatrix::Multiply "int
Epetra_SerialDenseMatrix::Multiply(char TransA, char TransB, double
ScalarAB, const Epetra_SerialDenseMatrix &A, const
Epetra_SerialDenseMatrix &B, double ScalarThis)

Matrix-Matrix multiplication, this = ScalarThis* this + ScalarAB*A*B.

This function performs a variety of matrix-matrix multiply operations.

Parameters:
-----------

In:  TransA - Operate with the transpose of A if = 'T', else no
transpose if = 'N'.

In:  TransB - Operate with the transpose of B if = 'T', else no
transpose if = 'N'.

In:  ScalarAB - Scalar to multiply with A*B.

In:  A - Dense Matrix.

In:  B - Dense Matrix.

In:  ScalarThis - Scalar to multiply with this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::Multiply "int
Epetra_SerialDenseMatrix::Multiply(bool transA, const
Epetra_SerialDenseMatrix &x, Epetra_SerialDenseMatrix &y)

Matrix-Vector multiplication, y = A*x, where 'this' == A. ";

%feature("docstring")  Epetra_SerialDenseMatrix::Multiply "int
Epetra_SerialDenseMatrix::Multiply(char SideA, double ScalarAB, const
Epetra_SerialSymDenseMatrix &A, const Epetra_SerialDenseMatrix &B,
double ScalarThis)

Matrix-Matrix multiplication with a symmetric matrix A.

If SideA = 'L', compute this = ScalarThis* this + ScalarAB*A*B. If
SideA = 'R', compute this = ScalarThis* this + ScalarAB*B*A.

This function performs a variety of matrix-matrix multiply operations.

Parameters:
-----------

In:  SideA - Specifies order of A relative to B.

In:  ScalarAB - Scalar to multiply with A*B.

In:  A - Symmetric Dense Matrix, either upper or lower triangle will
be used depending on value of A.Upper().

In:  B - Dense Matrix.

In:  ScalarThis - Scalar to multiply with this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::Scale "int
Epetra_SerialDenseMatrix::Scale(double ScalarA)

Inplace scalar-matrix product A = a A.

Scale a matrix, entry-by-entry using the value ScalarA.

Parameters:
-----------

ScalarA:  (In) Scalar to multiply with A.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::NormOne "double
Epetra_SerialDenseMatrix::NormOne() const

Computes the 1-Norm of the this matrix.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::NormInf "double
Epetra_SerialDenseMatrix::NormInf() const

Computes the Infinity-Norm of the this matrix. ";

/*  Data Accessor methods  */

%feature("docstring")  Epetra_SerialDenseMatrix::Random "int
Epetra_SerialDenseMatrix::Random()

Set matrix values to random numbers.

SerialDenseMatrix uses the random number generator provided by
Epetra_Util. The matrix values will be set to random values on the
interval (-1.0, 1.0).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::M "int
Epetra_SerialDenseMatrix::M() const

Returns row dimension of system. ";

%feature("docstring")  Epetra_SerialDenseMatrix::N "int
Epetra_SerialDenseMatrix::N() const

Returns column dimension of system. ";

%feature("docstring")  Epetra_SerialDenseMatrix::A "double*
Epetra_SerialDenseMatrix::A() const

Returns pointer to the this matrix. ";

%feature("docstring")  Epetra_SerialDenseMatrix::A "double*
Epetra_SerialDenseMatrix::A()

Returns pointer to the this matrix. ";

%feature("docstring")  Epetra_SerialDenseMatrix::LDA "int
Epetra_SerialDenseMatrix::LDA() const

Returns the leading dimension of the this matrix. ";

%feature("docstring")  Epetra_SerialDenseMatrix::CV "Epetra_DataAccess Epetra_SerialDenseMatrix::CV() const

Returns the data access mode of the this matrix. ";

/*  I/O methods  */

%feature("docstring")  Epetra_SerialDenseMatrix::Print "void
Epetra_SerialDenseMatrix::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";

/*  Deprecated methods (will be removed in later versions of this
class)  */

%feature("docstring")  Epetra_SerialDenseMatrix::OneNorm "virtual
double Epetra_SerialDenseMatrix::OneNorm() const

Computes the 1-Norm of the this matrix (identical to NormOne()
method).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::InfNorm "virtual
double Epetra_SerialDenseMatrix::InfNorm() const

Computes the Infinity-Norm of the this matrix (identical to NormInf()
method). ";

/*  Additional methods to support Epetra_SerialDenseOperator interface
*/

%feature("docstring")  Epetra_SerialDenseMatrix::SetUseTranspose "virtual int Epetra_SerialDenseMatrix::SetUseTranspose(bool
UseTranspose_in)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose -If true, multiply by the transpose of operator,
otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Epetra_SerialDenseMatrix::Apply "int
Epetra_SerialDenseMatrix::Apply(const Epetra_SerialDenseMatrix &X,
Epetra_SerialDenseMatrix &Y)

Returns the result of a Epetra_SerialDenseOperator applied to a
Epetra_SerialDenseMatrix X in Y.

Parameters:
-----------

In:  X - A Epetra_SerialDenseMatrix to multiply with operator.

Out:  Y -A Epetra_SerialDenseMatrix containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::ApplyInverse "virtual int Epetra_SerialDenseMatrix::ApplyInverse(const
Epetra_SerialDenseMatrix &X, Epetra_SerialDenseMatrix &Y)

Returns the result of a Epetra_SerialDenseOperator inverse applied to
an Epetra_SerialDenseMatrix X in Y.

Parameters:
-----------

In:  X - A Epetra_SerialDenseMatrix to solve for.

Out:  Y -A Epetra_SerialDenseMatrix containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseMatrix::Label "virtual const
char* Epetra_SerialDenseMatrix::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_SerialDenseMatrix::UseTranspose "virtual bool Epetra_SerialDenseMatrix::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_SerialDenseMatrix::HasNormInf "virtual
bool Epetra_SerialDenseMatrix::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Epetra_SerialDenseMatrix::RowDim "virtual int
Epetra_SerialDenseMatrix::RowDim() const

Returns the row dimension of operator. ";

%feature("docstring")  Epetra_SerialDenseMatrix::ColDim "virtual int
Epetra_SerialDenseMatrix::ColDim() const

Returns the column dimension of operator. ";


// File: classEpetra__SerialDenseOperator.xml
%feature("docstring") Epetra_SerialDenseOperator "

Epetra_SerialDenseOperator: A pure virtual class for using real-valued
double-precision operators.

The Epetra_SerialDenseOperator class is a pure virtual class
(specifies interface only) that enable the use of real-valued double-
precision operators. It is currently implemented by the
Epetra_SerialDenseMatrix, Epetra_SerialDenseSolver and
Epetra_SerialDenseSVD classes.

C++ includes: Epetra_SerialDenseOperator.h ";

/*  Destructor  */

%feature("docstring")
Epetra_SerialDenseOperator::~Epetra_SerialDenseOperator "virtual
Epetra_SerialDenseOperator::~Epetra_SerialDenseOperator()

Destructor. ";

/*  Attribute set methods  */

%feature("docstring")  Epetra_SerialDenseOperator::SetUseTranspose "virtual int Epetra_SerialDenseOperator::SetUseTranspose(bool
UseTranspose)=0

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  UseTranspose -If true, multiply by the transpose of operator,
otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

/*  Mathematical functions  */

%feature("docstring")  Epetra_SerialDenseOperator::Apply "virtual int
Epetra_SerialDenseOperator::Apply(const Epetra_SerialDenseMatrix &X,
Epetra_SerialDenseMatrix &Y)=0

Returns the result of a Epetra_SerialDenseOperator applied to a
Epetra_SerialDenseMatrix X in Y.

Parameters:
-----------

In:  X - A Epetra_SerialDenseMatrix to multiply with operator.

Out:  Y -A Epetra_SerialDenseMatrix containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseOperator::ApplyInverse "virtual int Epetra_SerialDenseOperator::ApplyInverse(const
Epetra_SerialDenseMatrix &X, Epetra_SerialDenseMatrix &Y)=0

Returns the result of a Epetra_SerialDenseOperator inverse applied to
an Epetra_SerialDenseMatrix X in Y.

Parameters:
-----------

In:  X - A Epetra_SerialDenseMatrix to solve for.

Out:  Y -A Epetra_SerialDenseMatrix containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseOperator::NormInf "virtual
double Epetra_SerialDenseOperator::NormInf() const =0

Returns the infinity norm of the global matrix. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_SerialDenseOperator::Label "virtual
const char* Epetra_SerialDenseOperator::Label() const =0

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_SerialDenseOperator::UseTranspose "virtual bool Epetra_SerialDenseOperator::UseTranspose() const =0

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_SerialDenseOperator::HasNormInf "virtual bool Epetra_SerialDenseOperator::HasNormInf() const =0

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Epetra_SerialDenseOperator::RowDim "virtual
int Epetra_SerialDenseOperator::RowDim() const =0

Returns the row dimension of operator. ";

%feature("docstring")  Epetra_SerialDenseOperator::ColDim "virtual
int Epetra_SerialDenseOperator::ColDim() const =0

Returns the column dimension of operator. ";


// File: classEpetra__SerialDenseSolver.xml
%feature("docstring") Epetra_SerialDenseSolver "

Epetra_SerialDenseSolver: A class for solving dense linear problems.

The Epetra_SerialDenseSolver class enables the definition, in terms of
Epetra_SerialDenseMatrix and Epetra_SerialDenseVector objects, of a
dense linear problem, followed by the solution of that problem via the
most sophisticated techniques available in LAPACK.

The Epetra_SerialDenseSolver class is intended to provide full-
featured support for solving linear problems for general dense
rectangular (or square) matrices. It is written on top of BLAS and
LAPACK and thus has excellent performance and numerical capabilities.
Using this class, one can either perform simple factorizations and
solves or apply all the tricks available in LAPACK to get the best
possible solution for very ill-conditioned problems.

Epetra_SerialDenseSolver vs. Epetra_LAPACK

The Epetra_LAPACK class provides access to most of the same
functionality as Epetra_SerialDenseSolver. The primary difference is
that Epetra_LAPACK is a \"thin\" layer on top of LAPACK and
Epetra_SerialDenseSolver attempts to provide easy access to the more
sophisticated aspects of solving dense linear and eigensystems. When
you should use Epetra_LAPACK: If you are simply looking for a
convenient wrapper around the Fortran LAPACK routines and you have a
well-conditioned problem, you should probably use Epetra_LAPACK
directly.

When you should use Epetra_SerialDenseSolver: If you want to (or
potentially want to) solve ill-conditioned problems or want to work
with a more object-oriented interface, you should probably use
Epetra_SerialDenseSolver.

Constructing Epetra_SerialDenseSolver Objects

There is a single Epetra_SerialDenseSolver constructor. However, the
matrix, right hand side and solution vectors must be set prior to
executing most methods in this class.

Setting vectors used for linear solves

The matrix A, the left hand side X and the right hand side B (when
solving AX = B, for X), can be set by appropriate set methods. Each of
these three objects must be an Epetra_SerialDenseMatrix or and
Epetra_SerialDenseVector object. The set methods are as follows:
SetMatrix() - Sets the matrix.

SetVectors() - Sets the left and right hand side vector(s).

Vector and Utility Functions

Once a Epetra_SerialDenseSolver is constructed, several mathematical
functions can be applied to the object. Specifically: Factorizations.

Solves.

Condition estimates.

Equilibration.

Norms.

Counting floating point operations The Epetra_SerialDenseSolver class
has Epetra_CompObject as a base class. Thus, floating point operations
are counted and accumulated in the Epetra_Flop object (if any) that
was set using the SetFlopCounter() method in the Epetra_CompObject
base class.

Strategies for Solving Linear Systems In many cases, linear systems
can be accurately solved by simply computing the LU factorization of
the matrix and then performing a forward back solve with a given set
of right hand side vectors. However, in some instances, the
factorization may be very poorly conditioned and this simple approach
may not work. In these situations, equilibration and iterative
refinement may improve the accuracy, or prevent a breakdown in the
factorization.

Epetra_SerialDenseSolver will use equilibration with the factorization
if, once the object is constructed and before it is factored, you call
the function FactorWithEquilibration(true) to force equilibration to
be used. If you are uncertain if equilibration should be used, you may
call the function ShouldEquilibrate() which will return true if
equilibration could possibly help. ShouldEquilibrate() uses guidelines
specified in the LAPACK User Guide, namely if SCOND < 0.1 and AMAX <
Underflow or AMAX > Overflow, to determine if equilibration might be
useful.

Epetra_SerialDenseSolver will use iterative refinement after a
forward/back solve if you call SolveToRefinedSolution(true). It will
also compute forward and backward error estimates if you call
EstimateSolutionErrors(true). Access to the forward (back) error
estimates is available via FERR() ( BERR()).

Examples using Epetra_SerialDenseSolver can be found in the Epetra
test directories.

C++ includes: Epetra_SerialDenseSolver.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")
Epetra_SerialDenseSolver::Epetra_SerialDenseSolver "Epetra_SerialDenseSolver::Epetra_SerialDenseSolver()

Default constructor; matrix should be set using SetMatrix(), LHS and
RHS set with SetVectors(). ";

%feature("docstring")
Epetra_SerialDenseSolver::~Epetra_SerialDenseSolver "Epetra_SerialDenseSolver::~Epetra_SerialDenseSolver()

Epetra_SerialDenseSolver destructor. ";

/*  Set Methods  */

%feature("docstring")  Epetra_SerialDenseSolver::SetMatrix "int
Epetra_SerialDenseSolver::SetMatrix(Epetra_SerialDenseMatrix &A)

Sets the pointers for coefficient matrix. ";

%feature("docstring")  Epetra_SerialDenseSolver::SetVectors "int
Epetra_SerialDenseSolver::SetVectors(Epetra_SerialDenseMatrix &X,
Epetra_SerialDenseMatrix &B)

Sets the pointers for left and right hand side vector(s).

Row dimension of X must match column dimension of matrix A, row
dimension of B must match row dimension of A. X and B must have the
same dimensions. ";

/*  Strategy modifying Methods  */

%feature("docstring")
Epetra_SerialDenseSolver::FactorWithEquilibration "void
Epetra_SerialDenseSolver::FactorWithEquilibration(bool Flag)

Causes equilibration to be called just before the matrix factorization
as part of the call to Factor.

This function must be called before the factorization is performed. ";

%feature("docstring")  Epetra_SerialDenseSolver::SolveWithTranspose "void Epetra_SerialDenseSolver::SolveWithTranspose(bool Flag)

If Flag is true, causes all subsequent function calls to work with the
transpose of this matrix, otherwise not. ";

%feature("docstring")
Epetra_SerialDenseSolver::SolveToRefinedSolution "void
Epetra_SerialDenseSolver::SolveToRefinedSolution(bool Flag)

Causes all solves to compute solution to best ability using iterative
refinement. ";

%feature("docstring")
Epetra_SerialDenseSolver::EstimateSolutionErrors "void
Epetra_SerialDenseSolver::EstimateSolutionErrors(bool Flag)

Causes all solves to estimate the forward and backward solution error.

Error estimates will be in the arrays FERR and BERR, resp, after the
solve step is complete. These arrays are accessible via the FERR() and
BERR() access functions. ";

/*  Factor/Solve/Invert Methods  */

%feature("docstring")  Epetra_SerialDenseSolver::Factor "int
Epetra_SerialDenseSolver::Factor(void)

Computes the in-place LU factorization of the matrix using the LAPACK
routine DGETRF.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseSolver::Solve "int
Epetra_SerialDenseSolver::Solve(void)

Computes the solution X to AX = B for the this matrix and the B
provided to SetVectors()..

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseSolver::Invert "int
Epetra_SerialDenseSolver::Invert(void)

Inverts the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")
Epetra_SerialDenseSolver::ComputeEquilibrateScaling "int
Epetra_SerialDenseSolver::ComputeEquilibrateScaling(void)

Computes the scaling vector S(i) = 1/sqrt(A(i,i)) of the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialDenseSolver::EquilibrateMatrix "int Epetra_SerialDenseSolver::EquilibrateMatrix(void)

Equilibrates the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialDenseSolver::EquilibrateRHS "int
Epetra_SerialDenseSolver::EquilibrateRHS(void)

Equilibrates the current RHS.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialDenseSolver::ApplyRefinement "int
Epetra_SerialDenseSolver::ApplyRefinement(void)

Apply Iterative Refinement.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialDenseSolver::UnequilibrateLHS "int Epetra_SerialDenseSolver::UnequilibrateLHS(void)

Unscales the solution vectors if equilibration was used to solve the
system.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")
Epetra_SerialDenseSolver::ReciprocalConditionEstimate "int
Epetra_SerialDenseSolver::ReciprocalConditionEstimate(double &Value)

Returns the reciprocal of the 1-norm condition number of the this
matrix.

Parameters:
-----------

Value:  Out On return contains the reciprocal of the 1-norm condition
number of the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

/*  Query methods  */

%feature("docstring")  Epetra_SerialDenseSolver::Transpose "bool
Epetra_SerialDenseSolver::Transpose()

Returns true if transpose of this matrix has and will be used. ";

%feature("docstring")  Epetra_SerialDenseSolver::Factored "bool
Epetra_SerialDenseSolver::Factored()

Returns true if matrix is factored (factor available via AF() and
LDAF()). ";

%feature("docstring")  Epetra_SerialDenseSolver::A_Equilibrated "bool
Epetra_SerialDenseSolver::A_Equilibrated()

Returns true if factor is equilibrated (factor available via AF() and
LDAF()). ";

%feature("docstring")  Epetra_SerialDenseSolver::B_Equilibrated "bool
Epetra_SerialDenseSolver::B_Equilibrated()

Returns true if RHS is equilibrated (RHS available via B() and LDB()).
";

%feature("docstring")  Epetra_SerialDenseSolver::ShouldEquilibrate "virtual bool Epetra_SerialDenseSolver::ShouldEquilibrate()

Returns true if the LAPACK general rules for equilibration suggest you
should equilibrate the system. ";

%feature("docstring")
Epetra_SerialDenseSolver::SolutionErrorsEstimated "bool
Epetra_SerialDenseSolver::SolutionErrorsEstimated()

Returns true if forward and backward error estimated have been
computed (available via FERR() and BERR()). ";

%feature("docstring")  Epetra_SerialDenseSolver::Inverted "bool
Epetra_SerialDenseSolver::Inverted()

Returns true if matrix inverse has been computed (inverse available
via AF() and LDAF()). ";

%feature("docstring")
Epetra_SerialDenseSolver::ReciprocalConditionEstimated "bool
Epetra_SerialDenseSolver::ReciprocalConditionEstimated()

Returns true if the condition number of the this matrix has been
computed (value available via ReciprocalConditionEstimate()). ";

%feature("docstring")  Epetra_SerialDenseSolver::Solved "bool
Epetra_SerialDenseSolver::Solved()

Returns true if the current set of vectors has been solved. ";

%feature("docstring")  Epetra_SerialDenseSolver::SolutionRefined "bool Epetra_SerialDenseSolver::SolutionRefined()

Returns true if the current set of vectors has been refined. ";

/*  Data Accessor methods  */

%feature("docstring")  Epetra_SerialDenseSolver::Matrix "Epetra_SerialDenseMatrix* Epetra_SerialDenseSolver::Matrix() const

Returns pointer to current matrix. ";

%feature("docstring")  Epetra_SerialDenseSolver::FactoredMatrix "Epetra_SerialDenseMatrix* Epetra_SerialDenseSolver::FactoredMatrix()
const

Returns pointer to factored matrix (assuming factorization has been
performed). ";

%feature("docstring")  Epetra_SerialDenseSolver::LHS "Epetra_SerialDenseMatrix* Epetra_SerialDenseSolver::LHS() const

Returns pointer to current LHS. ";

%feature("docstring")  Epetra_SerialDenseSolver::RHS "Epetra_SerialDenseMatrix* Epetra_SerialDenseSolver::RHS() const

Returns pointer to current RHS. ";

%feature("docstring")  Epetra_SerialDenseSolver::M "int
Epetra_SerialDenseSolver::M() const

Returns row dimension of system. ";

%feature("docstring")  Epetra_SerialDenseSolver::N "int
Epetra_SerialDenseSolver::N() const

Returns column dimension of system. ";

%feature("docstring")  Epetra_SerialDenseSolver::A "double*
Epetra_SerialDenseSolver::A() const

Returns pointer to the this matrix. ";

%feature("docstring")  Epetra_SerialDenseSolver::LDA "int
Epetra_SerialDenseSolver::LDA() const

Returns the leading dimension of the this matrix. ";

%feature("docstring")  Epetra_SerialDenseSolver::B "double*
Epetra_SerialDenseSolver::B() const

Returns pointer to current RHS. ";

%feature("docstring")  Epetra_SerialDenseSolver::LDB "int
Epetra_SerialDenseSolver::LDB() const

Returns the leading dimension of the RHS. ";

%feature("docstring")  Epetra_SerialDenseSolver::NRHS "int
Epetra_SerialDenseSolver::NRHS() const

Returns the number of current right hand sides and solution vectors.
";

%feature("docstring")  Epetra_SerialDenseSolver::X "double*
Epetra_SerialDenseSolver::X() const

Returns pointer to current solution. ";

%feature("docstring")  Epetra_SerialDenseSolver::LDX "int
Epetra_SerialDenseSolver::LDX() const

Returns the leading dimension of the solution. ";

%feature("docstring")  Epetra_SerialDenseSolver::AF "double*
Epetra_SerialDenseSolver::AF() const

Returns pointer to the factored matrix (may be the same as A() if
factorization done in place). ";

%feature("docstring")  Epetra_SerialDenseSolver::LDAF "int
Epetra_SerialDenseSolver::LDAF() const

Returns the leading dimension of the factored matrix. ";

%feature("docstring")  Epetra_SerialDenseSolver::IPIV "int*
Epetra_SerialDenseSolver::IPIV() const

Returns pointer to pivot vector (if factorization has been computed),
zero otherwise. ";

%feature("docstring")  Epetra_SerialDenseSolver::ANORM "double
Epetra_SerialDenseSolver::ANORM() const

Returns the 1-Norm of the this matrix (returns -1 if not yet
computed). ";

%feature("docstring")  Epetra_SerialDenseSolver::RCOND "double
Epetra_SerialDenseSolver::RCOND() const

Returns the reciprocal of the condition number of the this matrix
(returns -1 if not yet computed). ";

%feature("docstring")  Epetra_SerialDenseSolver::ROWCND "double
Epetra_SerialDenseSolver::ROWCND() const

Ratio of smallest to largest row scale factors for the this matrix
(returns -1 if not yet computed).

If ROWCND() is >= 0.1 and AMAX() is not close to overflow or
underflow, then equilibration is not needed. ";

%feature("docstring")  Epetra_SerialDenseSolver::COLCND "double
Epetra_SerialDenseSolver::COLCND() const

Ratio of smallest to largest column scale factors for the this matrix
(returns -1 if not yet computed).

If COLCND() is >= 0.1 then equilibration is not needed. ";

%feature("docstring")  Epetra_SerialDenseSolver::AMAX "double
Epetra_SerialDenseSolver::AMAX() const

Returns the absolute value of the largest entry of the this matrix
(returns -1 if not yet computed). ";

%feature("docstring")  Epetra_SerialDenseSolver::FERR "double*
Epetra_SerialDenseSolver::FERR() const

Returns a pointer to the forward error estimates computed by LAPACK.
";

%feature("docstring")  Epetra_SerialDenseSolver::BERR "double*
Epetra_SerialDenseSolver::BERR() const

Returns a pointer to the backward error estimates computed by LAPACK.
";

%feature("docstring")  Epetra_SerialDenseSolver::R "double*
Epetra_SerialDenseSolver::R() const

Returns a pointer to the row scaling vector used for equilibration. ";

%feature("docstring")  Epetra_SerialDenseSolver::C "double*
Epetra_SerialDenseSolver::C() const

Returns a pointer to the column scale vector used for equilibration.
";

/*  I/O methods  */

%feature("docstring")  Epetra_SerialDenseSolver::Print "void
Epetra_SerialDenseSolver::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";


// File: classEpetra__SerialDenseSVD.xml
%feature("docstring") Epetra_SerialDenseSVD "

Epetra_SerialDenseSVD: A class for SVDing dense linear problems.

The Epetra_SerialDenseSVD class enables the definition, in terms of
Epetra_SerialDenseMatrix and Epetra_SerialDenseVector objects, of a
dense linear problem, followed by the solution of that problem via the
most sophisticated techniques available in LAPACK.

The Epetra_SerialDenseSVD class is intended to provide full-featured
support for solving linear problems for general dense rectangular (or
square) matrices. It is written on top of BLAS and LAPACK and thus has
excellent performance and numerical capabilities. Using this class,
one can either perform simple factorizations and solves or apply all
the tricks available in LAPACK to get the best possible solution for
very ill-conditioned problems.

Epetra_SerialDenseSVD vs. Epetra_LAPACK

The Epetra_LAPACK class provides access to most of the same
functionality as Epetra_SerialDenseSolver. The primary difference is
that Epetra_LAPACK is a \"thin\" layer on top of LAPACK and
Epetra_SerialDenseSolver attempts to provide easy access to the more
sophisticated aspects of solving dense linear and eigensystems. When
you should use Epetra_LAPACK: If you are simply looking for a
convenient wrapper around the Fortran LAPACK routines and you have a
well-conditioned problem, you should probably use Epetra_LAPACK
directly.

When you should use Epetra_SerialDenseSolver: If you want to (or
potentially want to) solve ill-conditioned problems or want to work
with a more object-oriented interface, you should probably use
Epetra_SerialDenseSolver.

Constructing Epetra_SerialDenseSVD Objects

There is a single Epetra_SerialDenseSVD constructor. However, the
matrix, right hand side and solution vectors must be set prior to
executing most methods in this class.

Setting vectors used for linear solves

The matrix A, the left hand side X and the right hand side B (when
solving AX = B, for X), can be set by appropriate set methods. Each of
these three objects must be an Epetra_SerialDenseMatrix or and
Epetra_SerialDenseVector object. The set methods are as follows:
SetMatrix() - Sets the matrix.

SetVectors() - Sets the left and right hand side vector(s).

Vector and Utility Functions

Once a Epetra_SerialDenseSVD is constructed, several mathematical
functions can be applied to the object. Specifically: Factorizations.

Solves.

Condition estimates.

Norms.

Counting floating point operations The Epetra_SerialDenseSVD class has
Epetra_CompObject as a base class. Thus, floating point operations are
counted and accumulated in the Epetra_Flop object (if any) that was
set using the SetFlopCounter() method in the Epetra_CompObject base
class.

Examples using Epetra_SerialDenseSVD can be found in the Epetra test
directories.

C++ includes: Epetra_SerialDenseSVD.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")  Epetra_SerialDenseSVD::Epetra_SerialDenseSVD "Epetra_SerialDenseSVD::Epetra_SerialDenseSVD()

Default constructor; matrix should be set using SetMatrix(), LHS and
RHS set with SetVectors(). ";

%feature("docstring")  Epetra_SerialDenseSVD::~Epetra_SerialDenseSVD "Epetra_SerialDenseSVD::~Epetra_SerialDenseSVD()

Epetra_SerialDenseSVD destructor. ";

/*  Set Methods  */

%feature("docstring")  Epetra_SerialDenseSVD::SetMatrix "int
Epetra_SerialDenseSVD::SetMatrix(Epetra_SerialDenseMatrix &A)

Sets the pointers for coefficient matrix. ";

%feature("docstring")  Epetra_SerialDenseSVD::SetVectors "int
Epetra_SerialDenseSVD::SetVectors(Epetra_SerialDenseMatrix &X,
Epetra_SerialDenseMatrix &B)

Sets the pointers for left and right hand side vector(s).

Row dimension of X must match column dimension of matrix A, row
dimension of B must match row dimension of A. X and B must have the
same dimensions. ";

/*  Strategy modifying Methods  */

%feature("docstring")  Epetra_SerialDenseSVD::SolveWithTranspose "void Epetra_SerialDenseSVD::SolveWithTranspose(bool Flag)

Causes equilibration to be called just before the matrix factorization
as part of the call to Factor.

This function must be called before the factorization is performed. If
Flag is true, causes all subsequent function calls to work with the
transpose of this matrix, otherwise not. ";

/*  Factor/Solve/Invert Methods  */

/* Causes all solves to compute solution to best ability using
iterative refinement.

*/

%feature("docstring")  Epetra_SerialDenseSVD::Factor "int
Epetra_SerialDenseSVD::Factor(void) ";

%feature("docstring")  Epetra_SerialDenseSVD::Solve "int
Epetra_SerialDenseSVD::Solve(void)

Computes the solution X to AX = B for the this matrix and the B
provided to SetVectors()..

Inverse of Matrix must be formed Integer error code, set to 0 if
successful. ";

%feature("docstring")  Epetra_SerialDenseSVD::Invert "int
Epetra_SerialDenseSVD::Invert(double rthresh=0.0, double athresh=0.0)

Inverts the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

/*  Query methods  */

%feature("docstring")  Epetra_SerialDenseSVD::Transpose "bool
Epetra_SerialDenseSVD::Transpose()

Returns true if transpose of this matrix has and will be used. ";

%feature("docstring")  Epetra_SerialDenseSVD::Factored "bool
Epetra_SerialDenseSVD::Factored()

Returns true if matrix is factored (factor available via AF() and
LDAF()). ";

%feature("docstring")  Epetra_SerialDenseSVD::Inverted "bool
Epetra_SerialDenseSVD::Inverted()

Returns true if matrix inverse has been computed (inverse available
via AF() and LDAF()). ";

%feature("docstring")  Epetra_SerialDenseSVD::Solved "bool
Epetra_SerialDenseSVD::Solved()

Returns true if the current set of vectors has been solved. ";

/*  Data Accessor methods  */

%feature("docstring")  Epetra_SerialDenseSVD::Matrix "Epetra_SerialDenseMatrix* Epetra_SerialDenseSVD::Matrix() const

Returns pointer to current matrix. ";

%feature("docstring")  Epetra_SerialDenseSVD::InvertedMatrix "Epetra_SerialDenseMatrix* Epetra_SerialDenseSVD::InvertedMatrix()
const

Returns pointer to inverted matrix (assuming inverse has been
performed). ";

%feature("docstring")  Epetra_SerialDenseSVD::LHS "Epetra_SerialDenseMatrix* Epetra_SerialDenseSVD::LHS() const

Returns pointer to current LHS. ";

%feature("docstring")  Epetra_SerialDenseSVD::RHS "Epetra_SerialDenseMatrix* Epetra_SerialDenseSVD::RHS() const

Returns pointer to current RHS. ";

%feature("docstring")  Epetra_SerialDenseSVD::M "int
Epetra_SerialDenseSVD::M() const

Returns row dimension of system. ";

%feature("docstring")  Epetra_SerialDenseSVD::N "int
Epetra_SerialDenseSVD::N() const

Returns column dimension of system. ";

%feature("docstring")  Epetra_SerialDenseSVD::A "double*
Epetra_SerialDenseSVD::A() const

Returns pointer to the this matrix. ";

%feature("docstring")  Epetra_SerialDenseSVD::LDA "int
Epetra_SerialDenseSVD::LDA() const

Returns the leading dimension of the this matrix. ";

%feature("docstring")  Epetra_SerialDenseSVD::B "double*
Epetra_SerialDenseSVD::B() const

Returns pointer to current RHS. ";

%feature("docstring")  Epetra_SerialDenseSVD::LDB "int
Epetra_SerialDenseSVD::LDB() const

Returns the leading dimension of the RHS. ";

%feature("docstring")  Epetra_SerialDenseSVD::NRHS "int
Epetra_SerialDenseSVD::NRHS() const

Returns the number of current right hand sides and solution vectors.
";

%feature("docstring")  Epetra_SerialDenseSVD::X "double*
Epetra_SerialDenseSVD::X() const

Returns pointer to current solution. ";

%feature("docstring")  Epetra_SerialDenseSVD::LDX "int
Epetra_SerialDenseSVD::LDX() const

Returns the leading dimension of the solution. ";

%feature("docstring")  Epetra_SerialDenseSVD::S "double*
Epetra_SerialDenseSVD::S() const ";

%feature("docstring")  Epetra_SerialDenseSVD::AI "double*
Epetra_SerialDenseSVD::AI() const

Returns pointer to the inverted matrix (may be the same as A() if
factorization done in place). ";

%feature("docstring")  Epetra_SerialDenseSVD::LDAI "int
Epetra_SerialDenseSVD::LDAI() const

Returns the leading dimension of the inverted matrix. ";

%feature("docstring")  Epetra_SerialDenseSVD::ANORM "double
Epetra_SerialDenseSVD::ANORM() const

Returns the 1-Norm of the this matrix (returns -1 if not yet
computed). ";

/*  I/O methods  */

%feature("docstring")  Epetra_SerialDenseSVD::Print "void
Epetra_SerialDenseSVD::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";

/*  Additional methods for support of Epetra_SerialDenseOperator
interface  */

%feature("docstring")  Epetra_SerialDenseSVD::SetUseTranspose "virtual int Epetra_SerialDenseSVD::SetUseTranspose(bool use_transpose)

If set true, transpose of this operator will be applied.

This flag allows the transpose of the given operator to be used
implicitly. Setting this flag affects only the Apply() and
ApplyInverse() methods. If the implementation of this interface does
not support transpose use, this method should return a value of -1.

Parameters:
-----------

In:  use_transpose -If true, multiply by the transpose of operator,
otherwise just use operator.

Integer error code, set to 0 if successful. Set to -1 if this
implementation does not support transpose. ";

%feature("docstring")  Epetra_SerialDenseSVD::Apply "virtual int
Epetra_SerialDenseSVD::Apply(const Epetra_SerialDenseMatrix &Xmat,
Epetra_SerialDenseMatrix &Ymat)

Returns the result of a Epetra_SerialDenseOperator applied to a
Epetra_SerialDenseMatrix X in Y.

Parameters:
-----------

In:  X - A Epetra_SerialDenseMatrix to multiply with operator.

Out:  Y -A Epetra_SerialDenseMatrix containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseSVD::ApplyInverse "virtual
int Epetra_SerialDenseSVD::ApplyInverse(const Epetra_SerialDenseMatrix
&Xmat, Epetra_SerialDenseMatrix &Ymat)

Returns the result of a Epetra_SerialDenseOperator inverse applied to
an Epetra_SerialDenseMatrix X in Y.

Parameters:
-----------

In:  X - A Epetra_SerialDenseMatrix to solve for.

Out:  Y -A Epetra_SerialDenseMatrix containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseSVD::NormInf "virtual double
Epetra_SerialDenseSVD::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Epetra_SerialDenseSVD::Label "virtual const
char* Epetra_SerialDenseSVD::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_SerialDenseSVD::UseTranspose "virtual
bool Epetra_SerialDenseSVD::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_SerialDenseSVD::HasNormInf "virtual
bool Epetra_SerialDenseSVD::HasNormInf() const

Returns true if the this object can provide an approximate Inf-norm,
false otherwise. ";

%feature("docstring")  Epetra_SerialDenseSVD::RowDim "virtual int
Epetra_SerialDenseSVD::RowDim() const

Returns the row dimension of operator. ";

%feature("docstring")  Epetra_SerialDenseSVD::ColDim "virtual int
Epetra_SerialDenseSVD::ColDim() const

Returns the column dimension of operator. ";

%feature("docstring")  Epetra_SerialDenseSVD::AllocateWORK "void
Epetra_SerialDenseSVD::AllocateWORK() ";

%feature("docstring")  Epetra_SerialDenseSVD::AllocateIWORK "void
Epetra_SerialDenseSVD::AllocateIWORK() ";

%feature("docstring")  Epetra_SerialDenseSVD::InitPointers "void
Epetra_SerialDenseSVD::InitPointers() ";

%feature("docstring")  Epetra_SerialDenseSVD::DeleteArrays "void
Epetra_SerialDenseSVD::DeleteArrays() ";

%feature("docstring")  Epetra_SerialDenseSVD::ResetMatrix "void
Epetra_SerialDenseSVD::ResetMatrix() ";

%feature("docstring")  Epetra_SerialDenseSVD::ResetVectors "void
Epetra_SerialDenseSVD::ResetVectors() ";


// File: classEpetra__SerialDenseVector.xml
%feature("docstring") Epetra_SerialDenseVector "

Epetra_SerialDenseVector: A class for constructing and using dense
vectors.

The Epetra_SerialDenseVector class enables the construction and use of
real-valued, double- precision dense vectors. It is built on the BLAS
and LAPACK and derives from the Epetra_SerialDenseMatrix class.

The Epetra_SerialDenseVector class is intended to provide convenient
vector notation but derives all signficant functionality from
Epetra_SerialDenseMatrix.

Constructing Epetra_SerialDenseVector Objects

There are four Epetra_SerialDenseVector constructors. The first
constructs a zero-length object which should be made to appropriate
length using the Size() or Resize() functions and then filled with the
[] or () operators. The second constructs an object sized to the
dimension specified, which should be filled with the [] or ()
operators. The third is a constructor that accepts user data as a 1D
array, and the fourth is a copy constructor. The third constructor has
two data access modes (specified by the Epetra_DataAccess argument):
Copy mode - Allocates memory and makes a copy of the user-provided
data. In this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Extracting Data from Epetra_SerialDenseVector
Objects

Once a Epetra_SerialDenseVector is constructed, it is possible to view
the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.  The final useful function is Flops().
Each Epetra_SerialDenseVector object keep track of the number of
serial floating point operations performed using the specified object
as the this argument to the function. The Flops() function returns
this number as a double precision number. Using this information, in
conjunction with the Epetra_Time class, one can get accurate parallel
performance numbers.

C++ includes: Epetra_SerialDenseVector.h ";

/*  Constructors/destructors  */

%feature("docstring")
Epetra_SerialDenseVector::Epetra_SerialDenseVector "Epetra_SerialDenseVector::Epetra_SerialDenseVector()

Default constructor; defines a zero size object.

Epetra_SerialDenseVector objects defined by the default constructor
should be sized with the Size() or Resize functions. Values should be
defined by using the [] or () operators. ";

%feature("docstring")
Epetra_SerialDenseVector::Epetra_SerialDenseVector "Epetra_SerialDenseVector::Epetra_SerialDenseVector(int Length)

Sized constructor; defines a variable-sized object.

Parameters:
-----------

In:  Length - Length of vector.

Epetra_SerialDenseVector objects defined by the sized constructor are
already sized to the dimension given as a parameter. All values are
initialized to 0. Calling this constructor is equivalent to using the
default constructor, and then calling the Size function on it. Values
should be defined by using the [] or () operators. ";

%feature("docstring")
Epetra_SerialDenseVector::Epetra_SerialDenseVector "Epetra_SerialDenseVector::Epetra_SerialDenseVector(Epetra_DataAccess
CV, double *Values, int Length)

Set object values from one-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Values - Pointer to an array of double precision numbers
containing the values.

In:  Length - Length of vector.

See Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_SerialDenseVector::Epetra_SerialDenseVector "Epetra_SerialDenseVector::Epetra_SerialDenseVector(const
Epetra_SerialDenseVector &Source)

Epetra_SerialDenseVector copy constructor. ";

%feature("docstring")
Epetra_SerialDenseVector::~Epetra_SerialDenseVector "Epetra_SerialDenseVector::~Epetra_SerialDenseVector()

Epetra_SerialDenseVector destructor. ";

/*  Post-construction modification routines  */

%feature("docstring")  Epetra_SerialDenseVector::Size "int
Epetra_SerialDenseVector::Size(int Length_in)

Set length of a Epetra_SerialDenseVector object; init values to zero.

Parameters:
-----------

In:  Length - Length of vector object.

Allows user to define the dimension of a Epetra_SerialDenseVector.
This function can be called at any point after construction. Any
values that were previously in this object are destroyed and the
resized vector starts off with all zero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseVector::Resize "int
Epetra_SerialDenseVector::Resize(int Length_in)

Resize a Epetra_SerialDenseVector object.

Parameters:
-----------

In:  Length - Length of vector object.

Allows user to define the dimension of a Epetra_SerialDenseVector.
This function can be called at any point after construction. Any
values that were previously in this object are copied into the new
size. If the new shape is smaller than the original, the first Length
values are copied to the new vector.

Integer error code, set to 0 if successful. ";

/*  Element access methods  */

/*  Mathematical methods  */

%feature("docstring")  Epetra_SerialDenseVector::Random "int
Epetra_SerialDenseVector::Random()

Set vector values to random numbers.

SerialDenseVector uses the random number generator provided by
Epetra_Util. The vector values will be set to random values on the
interval (-1.0, 1.0).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialDenseVector::Dot "double
Epetra_SerialDenseVector::Dot(const Epetra_SerialDenseVector &x) const

Compute 1-norm of each vector in multi-vector.

Parameters:
-----------

x:  (In) Input vector x.

Dot-product of the this vector and x. ";

%feature("docstring")  Epetra_SerialDenseVector::Norm1 "double
Epetra_SerialDenseVector::Norm1() const

Compute 1-norm of each vector in multi-vector.

1-norm of the vector. ";

%feature("docstring")  Epetra_SerialDenseVector::Norm2 "double
Epetra_SerialDenseVector::Norm2() const

Compute 2-norm of each vector in multi-vector.

Parameters:
-----------

Out:

2-norm of the vector. ";

%feature("docstring")  Epetra_SerialDenseVector::NormInf "double
Epetra_SerialDenseVector::NormInf() const

Compute Inf-norm of each vector in multi-vector.

Infinity-norm of the vector. ";

/*  Attribute access methods  */

%feature("docstring")  Epetra_SerialDenseVector::Length "int
Epetra_SerialDenseVector::Length() const

Returns length of vector. ";

%feature("docstring")  Epetra_SerialDenseVector::Values "double*
Epetra_SerialDenseVector::Values() const

Returns pointer to the values in vector. ";

%feature("docstring")  Epetra_SerialDenseVector::CV "Epetra_DataAccess Epetra_SerialDenseVector::CV() const

Returns the data access mode of the this vector. ";

/*  I/O methods  */

%feature("docstring")  Epetra_SerialDenseVector::Print "void
Epetra_SerialDenseVector::Print(ostream &os) const

Print service methods; defines behavior of ostream << operator. ";


// File: classEpetra__SerialDistributor.xml
%feature("docstring") Epetra_SerialDistributor "

Epetra_SerialDistributor: The Epetra Serial implementation of the
Epetra_Distributor Gather/Scatter Setup Class.

The Epetra_SerialDistributor class is an Serial implement of
Epetra_Distributor that is essentially a trivial class since a serial
machine is a trivial parallel machine. An Epetra_SerialDistributor
object is actually produced by calling a method in the
Epetra_SerialComm class.

C++ includes: Epetra_SerialDistributor.h ";

/*  Constructor/Destructor  */

%feature("docstring")
Epetra_SerialDistributor::Epetra_SerialDistributor "Epetra_SerialDistributor::Epetra_SerialDistributor(const
Epetra_SerialComm &Comm)

Constructor. ";

%feature("docstring")
Epetra_SerialDistributor::Epetra_SerialDistributor "Epetra_SerialDistributor::Epetra_SerialDistributor(const
Epetra_SerialDistributor &Plan)

Epetra_SerialDistributor Copy Constructor. ";

%feature("docstring")  Epetra_SerialDistributor::Clone "Epetra_Distributor* Epetra_SerialDistributor::Clone()

Clone method. ";

%feature("docstring")
Epetra_SerialDistributor::~Epetra_SerialDistributor "Epetra_SerialDistributor::~Epetra_SerialDistributor()

Epetra_Comm Destructor. ";

%feature("docstring")  Epetra_SerialDistributor::CreateFromSends "int
Epetra_SerialDistributor::CreateFromSends(const int &NumExportIDs,
const int *ExportPIDs, bool Deterministic, int &NumRemoteIDs)

Create Distributor object using list of process IDs to which we
export.

Take a list of Process IDs and construct a plan for efficiently
scattering to these processes. Return the number of IDs being sent to
me.

Parameters:
-----------

NumExportIDs:  In Number of IDs that need to be sent from this
processor.

ExportPIDs:  In List of processors that will get the exported IDs.

Deterministic:  In No op.

NumRemoteIDs:  Out Number of IDs this processor will be receiving. ";

%feature("docstring")  Epetra_SerialDistributor::CreateFromRecvs "int
Epetra_SerialDistributor::CreateFromRecvs(const int &NumRemoteIDs,
const int *RemoteGIDs, const int *RemotePIDs, bool Deterministic, int
&NumExportIDs, int *&ExportGIDs, int *&ExportPIDs)

Create Distributor object using list of Remote global IDs and
corresponding PIDs.

Take a list of global IDs and construct a plan for efficiently
scattering to these processes. Return the number and list of IDs being
sent by me.

Parameters:
-----------

NumRemoteIDs:  In Number of IDs this processor will be receiving.

RemoteGIDs:  In List of IDs that this processor wants.

RemotePIDs:  In List of processors that will send the remote IDs.

Deterministic:  In No op.

NumExportIDs:  Out Number of IDs that need to be sent from this
processor.

ExportPIDs:  Out List of processors that will get the exported IDs. ";

%feature("docstring")  Epetra_SerialDistributor::CreateFromRecvs "int
Epetra_SerialDistributor::CreateFromRecvs(const int &NumRemoteIDs,
const long long *RemoteGIDs, const int *RemotePIDs, bool
Deterministic, int &NumExportIDs, long long *&ExportGIDs, int
*&ExportPIDs) ";

%feature("docstring")  Epetra_SerialDistributor::Do "int
Epetra_SerialDistributor::Do(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)

Execute plan on buffer of export objects in a single step. ";

%feature("docstring")  Epetra_SerialDistributor::DoReverse "int
Epetra_SerialDistributor::DoReverse(char *export_objs, int obj_size,
int &len_import_objs, char *&import_objs)

Execute reverse of plan on buffer of export objects in a single step.
";

%feature("docstring")  Epetra_SerialDistributor::DoPosts "int
Epetra_SerialDistributor::DoPosts(char *export_objs, int obj_size, int
&len_import_objs, char *&import_objs)

Post buffer of export objects (can do other local work before
executing Waits) ";

%feature("docstring")  Epetra_SerialDistributor::DoWaits "int
Epetra_SerialDistributor::DoWaits()

Wait on a set of posts. ";

%feature("docstring")  Epetra_SerialDistributor::DoReversePosts "int
Epetra_SerialDistributor::DoReversePosts(char *export_objs, int
obj_size, int &len_import_objs, char *&import_objs)

Do reverse post of buffer of export objects (can do other local work
before executing Waits) ";

%feature("docstring")  Epetra_SerialDistributor::DoReverseWaits "int
Epetra_SerialDistributor::DoReverseWaits()

Wait on a reverse set of posts. ";

%feature("docstring")  Epetra_SerialDistributor::Do "int
Epetra_SerialDistributor::Do(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Execute plan on buffer of export objects in a single step (object size
may vary) ";

%feature("docstring")  Epetra_SerialDistributor::DoReverse "int
Epetra_SerialDistributor::DoReverse(char *export_objs, int obj_size,
int *&sizes, int &len_import_objs, char *&import_objs)

Execute reverse of plan on buffer of export objects in a single step
(object size may vary) ";

%feature("docstring")  Epetra_SerialDistributor::DoPosts "int
Epetra_SerialDistributor::DoPosts(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Post buffer of export objects (can do other local work before
executing Waits) ";

%feature("docstring")  Epetra_SerialDistributor::DoReversePosts "int
Epetra_SerialDistributor::DoReversePosts(char *export_objs, int
obj_size, int *&sizes, int &len_import_objs, char *&import_objs)

Do reverse post of buffer of export objects (can do other local work
before executing Waits) ";

%feature("docstring")  Epetra_SerialDistributor::Print "void
Epetra_SerialDistributor::Print(ostream &os) const ";


// File: classEpetra__SerialSpdDenseSolver.xml
%feature("docstring") Epetra_SerialSpdDenseSolver "

Epetra_SerialSpdDenseSolver: A class for constructing and using
symmetric positive definite dense matrices.

The Epetra_SerialSpdDenseSolver class enables the construction and use
of real-valued, symmetric positive definite, double-precision dense
matrices. It is built on the Epetra_DenseMatrix class which in turn is
built on the BLAS and LAPACK via the Epetra_BLAS and Epetra_LAPACK
classes.

The Epetra_SerialSpdDenseSolver class is intended to provide full-
featured support for solving linear and eigen system problems for
symmetric positive definite matrices. It is written on top of BLAS and
LAPACK and thus has excellent performance and numerical capabilities.
Using this class, one can either perform simple factorizations and
solves or apply all the tricks available in LAPACK to get the best
possible solution for very ill-conditioned problems.

Epetra_SerialSpdDenseSolver vs. Epetra_LAPACK

The Epetra_LAPACK class provides access to most of the same
functionality as Epetra_SerialSpdDenseSolver. The primary difference
is that Epetra_LAPACK is a \"thin\" layer on top of LAPACK and
Epetra_SerialSpdDenseSolver attempts to provide easy access to the
more sophisticated aspects of solving dense linear and eigensystems.
When you should use Epetra_LAPACK: If you are simply looking for a
convenient wrapper around the Fortran LAPACK routines and you have a
well-conditioned problem, you should probably use Epetra_LAPACK
directly.

When you should use Epetra_SerialSpdDenseSolver: If you want to (or
potentially want to) solve ill-conditioned problems or want to work
with a more object-oriented interface, you should probably use
Epetra_SerialSpdDenseSolver.

Constructing Epetra_SerialSpdDenseSolver Objects

There are three Epetra_DenseMatrix constructors. The first constructs
a zero-sized object which should be made to appropriate length using
the Shape() or Reshape() functions and then filled with the [] or ()
operators. The second is a constructor that accepts user data as a 2D
array, the third is a copy constructor. The second constructor has two
data access modes (specified by the Epetra_DataAccess argument): Copy
mode - Allocates memory and makes a copy of the user-provided data. In
this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Setting vectors used for linear solves

Setting the X and B vectors (which are Epetra_DenseMatrix objects)
used for solving linear systems is done separately from the
constructor. This allows a single matrix factor to be used for
multiple solves. Similar to the constructor, the vectors X and B can
be copied or viewed using the Epetra_DataAccess argument.

Extracting Data from Epetra_SerialSpdDenseSolver Objects

Once a Epetra_SerialSpdDenseSolver is constructed, it is possible to
view the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.  Vector and Utility Functions

Once a Epetra_SerialSpdDenseSolver is constructed, several
mathematical functions can be applied to the object. Specifically:
Factorizations.

Solves.

Condition estimates.

Equilibration.

Norms.

The final useful function is Flops(). Each Epetra_SerialSpdDenseSolver
object keep track of the number of serial floating point operations
performed using the specified object as the this argument to the
function. The Flops() function returns this number as a double
precision number. Using this information, in conjunction with the
Epetra_Time class, one can get accurate parallel performance numbers.

Strategies for Solving Linear Systems In many cases, linear systems
can be accurately solved by simply computing the Cholesky
factorization of the matrix and then performing a forward back solve
with a given set of right hand side vectors. However, in some
instances, the factorization may be very poorly conditioned and the
simple approach may not work. In these situations, equilibration and
iterative refinement may improve the accuracy, or prevent a breakdown
in the factorization.

Epetra_SerialSpdDenseSolver will use equilibration with the
factorization if, once the object is constructed and before it is
factored, you call the function FactorWithEquilibration(true) to force
equilibration to be used. If you are uncertain if equilibration should
be used, you may call the function ShouldEquilibrate() which will
return true if equilibration could possibly help. ShouldEquilibrate()
uses guidelines specified in the LAPACK User Guide, namely if SCOND <
0.1 and AMAX < Underflow or AMAX > Overflow, to determine if
equilibration might be useful.

Epetra_SerialSpdDenseSolver will use iterative refinement after a
forward/back solve if you call SolveToRefinedSolution(true). It will
also compute forward and backward error estimates if you call
EstimateSolutionErrors(true). Access to the forward (back) error
estimates is available via FERR() ( BERR()).

Examples using Epetra_SerialSpdDenseSolver can be found in the Epetra
test directories.

C++ includes: Epetra_SerialSpdDenseSolver.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")
Epetra_SerialSpdDenseSolver::Epetra_SerialSpdDenseSolver "Epetra_SerialSpdDenseSolver::Epetra_SerialSpdDenseSolver()

Default constructor; matrix should be set using SetMatrix(), LHS and
RHS set with SetVectors(). ";

%feature("docstring")
Epetra_SerialSpdDenseSolver::~Epetra_SerialSpdDenseSolver "Epetra_SerialSpdDenseSolver::~Epetra_SerialSpdDenseSolver()

Epetra_SerialDenseSolver destructor. ";

/*  Set Methods  */

%feature("docstring")  Epetra_SerialSpdDenseSolver::SetMatrix "int
Epetra_SerialSpdDenseSolver::SetMatrix(Epetra_SerialSymDenseMatrix
&A_in)

Sets the pointers for coefficient matrix; special version for
symmetric matrices. ";

/*  Factor/Solve/Invert Methods  */

%feature("docstring")  Epetra_SerialSpdDenseSolver::Factor "int
Epetra_SerialSpdDenseSolver::Factor(void)

Computes the in-place Cholesky factorization of the matrix using the
LAPACK routine DPOTRF.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::Solve "int
Epetra_SerialSpdDenseSolver::Solve(void)

Computes the solution X to AX = B for the this matrix and the B
provided to SetVectors()..

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::Invert "int
Epetra_SerialSpdDenseSolver::Invert(void)

Inverts the this matrix.

Note: This function works a little differently that DPOTRI in that it
fills the entire matrix with the inverse, independent of the UPLO
specification.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")
Epetra_SerialSpdDenseSolver::ComputeEquilibrateScaling "int
Epetra_SerialSpdDenseSolver::ComputeEquilibrateScaling(void)

Computes the scaling vector S(i) = 1/sqrt(A(i,i) of the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::EquilibrateMatrix
"int Epetra_SerialSpdDenseSolver::EquilibrateMatrix(void)

Equilibrates the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::EquilibrateRHS "int Epetra_SerialSpdDenseSolver::EquilibrateRHS(void)

Equilibrates the current RHS.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::ApplyRefinement "int Epetra_SerialSpdDenseSolver::ApplyRefinement(void)

Apply Iterative Refinement.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::UnequilibrateLHS "int Epetra_SerialSpdDenseSolver::UnequilibrateLHS(void)

Unscales the solution vectors if equilibration was used to solve the
system.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

%feature("docstring")
Epetra_SerialSpdDenseSolver::ReciprocalConditionEstimate "int
Epetra_SerialSpdDenseSolver::ReciprocalConditionEstimate(double
&Value)

Returns the reciprocal of the 1-norm condition number of the this
matrix.

Parameters:
-----------

Value:  Out On return contains the reciprocal of the 1-norm condition
number of the this matrix.

Integer error code, set to 0 if successful. Otherwise returns the
LAPACK error code INFO. ";

/*  Query methods  */

%feature("docstring")  Epetra_SerialSpdDenseSolver::ShouldEquilibrate
"bool Epetra_SerialSpdDenseSolver::ShouldEquilibrate()

Returns true if the LAPACK general rules for equilibration suggest you
should equilibrate the system. ";

/*  Data Accessor methods  */

%feature("docstring")  Epetra_SerialSpdDenseSolver::SymMatrix "Epetra_SerialSymDenseMatrix* Epetra_SerialSpdDenseSolver::SymMatrix()
const

Returns pointer to current matrix. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::SymFactoredMatrix
"Epetra_SerialSymDenseMatrix*
Epetra_SerialSpdDenseSolver::SymFactoredMatrix() const

Returns pointer to factored matrix (assuming factorization has been
performed). ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::SCOND "double
Epetra_SerialSpdDenseSolver::SCOND()

Ratio of smallest to largest equilibration scale factors for the this
matrix (returns -1 if not yet computed).

If SCOND() is >= 0.1 and AMAX() is not close to overflow or underflow,
then equilibration is not needed. ";

%feature("docstring")  Epetra_SerialSpdDenseSolver::AMAX "double
Epetra_SerialSpdDenseSolver::AMAX()

Returns the absolute value of the largest entry of the this matrix
(returns -1 if not yet computed). ";


// File: classEpetra__SerialSymDenseMatrix.xml
%feature("docstring") Epetra_SerialSymDenseMatrix "

Epetra_SerialSymDenseMatrix: A class for constructing and using
symmetric positive definite dense matrices.

The Epetra_SerialSymDenseMatrix class enables the construction and use
of real-valued, symmetric positive definite, double-precision dense
matrices. It is built on the Epetra_SerialDenseMatrix class which in
turn is built on the BLAS via the Epetra_BLAS class.

The Epetra_SerialSymDenseMatrix class is intended to provide full-
featured support for solving linear and eigen system problems for
symmetric positive definite matrices. It is written on top of BLAS and
LAPACK and thus has excellent performance and numerical capabilities.
Using this class, one can either perform simple factorizations and
solves or apply all the tricks available in LAPACK to get the best
possible solution for very ill-conditioned problems.

Epetra_SerialSymDenseMatrix vs. Epetra_LAPACK

The Epetra_LAPACK class provides access to most of the same
functionality as Epetra_SerialSymDenseMatrix. The primary difference
is that Epetra_LAPACK is a \"thin\" layer on top of LAPACK and
Epetra_SerialSymDenseMatrix attempts to provide easy access to the
more sophisticated aspects of solving dense linear and eigensystems.
When you should use Epetra_LAPACK: If you are simply looking for a
convenient wrapper around the Fortran LAPACK routines and you have a
well-conditioned problem, you should probably use Epetra_LAPACK
directly.

When you should use Epetra_SerialSymDenseMatrix: If you want to (or
potentially want to) solve ill-conditioned problems or want to work
with a more object-oriented interface, you should probably use
Epetra_SerialSymDenseMatrix.

Constructing Epetra_SerialSymDenseMatrix Objects

There are three Epetra_DenseMatrix constructors. The first constructs
a zero-sized object which should be made to appropriate length using
the Shape() or Reshape() functions and then filled with the [] or ()
operators. The second is a constructor that accepts user data as a 2D
array, the third is a copy constructor. The second constructor has two
data access modes (specified by the Epetra_DataAccess argument): Copy
mode - Allocates memory and makes a copy of the user-provided data. In
this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the object.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  Extracting Data from Epetra_SerialSymDenseMatrix
Objects

Once a Epetra_SerialSymDenseMatrix is constructed, it is possible to
view the data via access functions.

WARNING:  Use of these access functions cam be extremely dangerous
from a data hiding perspective.  Vector and Utility Functions

Once a Epetra_SerialSymDenseMatrix is constructed, several
mathematical functions can be applied to the object. Specifically:
Multiplication.

Norms.

Counting floating point operations The Epetra_SerialSymDenseMatrix
class has Epetra_CompObject as a base class. Thus, floating point
operations are counted and accumulated in the Epetra_Flop object (if
any) that was set using the SetFlopCounter() method in the
Epetra_CompObject base class.

C++ includes: Epetra_SerialSymDenseMatrix.h ";

/*  Constructor/Destructor Methods  */

%feature("docstring")
Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix "Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix(void)

Default constructor; defines a zero size object.

Epetra_SerialSymDenseMatrix objects defined by the default constructor
should be sized with the Shape() or Reshape() functions. Values should
be defined by using the [] or ()operators.

Note: By default the active part of the matrix is assumed to be in the
lower triangle. To set the upper part as active, call SetUpper(). See
Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix "Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix(Epetra_DataAccess
CV, double *A, int LDA, int NumRowsCols)

Set object values from two-dimensional array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  A - Pointer to an array of double precision numbers. The first
vector starts at A. The second vector starts at A+LDA, the third at
A+2*LDA, and so on.

In:  LDA - The \"Leading Dimension\", or stride between vectors in
memory.

In:  NumRowsCols - Number of rows and columns in object.

Note: By default the active part of the matrix is assumed to be in the
lower triangle. To set the upper part as active, call SetUpper(). See
Detailed Description section for further discussion. ";

%feature("docstring")
Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix "Epetra_SerialSymDenseMatrix::Epetra_SerialSymDenseMatrix(const
Epetra_SerialSymDenseMatrix &Source)

Epetra_SerialSymDenseMatrix copy constructor. ";

%feature("docstring")
Epetra_SerialSymDenseMatrix::~Epetra_SerialSymDenseMatrix "Epetra_SerialSymDenseMatrix::~Epetra_SerialSymDenseMatrix()

Epetra_SerialSymDenseMatrix destructor. ";

/*  Set Methods  */

%feature("docstring")  Epetra_SerialSymDenseMatrix::Shape "int
Epetra_SerialSymDenseMatrix::Shape(int NumRowsCols)

Set dimensions of a Epetra_SerialSymDenseMatrix object; init values to
zero.

Parameters:
-----------

In:  NumRowsCols - Number of rows and columns in object.

Allows user to define the dimensions of a Epetra_DenseMatrix at any
point. This function can be called at any point after construction.
Any values that were previously in this object are destroyed and the
resized matrix starts off with all zero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::Reshape "int
Epetra_SerialSymDenseMatrix::Reshape(int NumRowsCols)

Reshape a Epetra_SerialSymDenseMatrix object.

Parameters:
-----------

In:  NumRowsCols - Number of rows and columns in object.

Allows user to define the dimensions of a Epetra_SerialSymDenseMatrix
at any point. This function can be called at any point after
construction. Any values that were previously in this object are
copied into the new shape. If the new shape is smaller than the
original, the upper left portion of the original matrix (the principal
submatrix) is copied to the new matrix.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::SetLower "void
Epetra_SerialSymDenseMatrix::SetLower()

Specify that the lower triangle of the this matrix should be used. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::SetUpper "void
Epetra_SerialSymDenseMatrix::SetUpper()

Specify that the upper triangle of the this matrix should be used. ";

/*  Query methods  */

%feature("docstring")  Epetra_SerialSymDenseMatrix::Upper "bool
Epetra_SerialSymDenseMatrix::Upper() const

Returns true if upper triangle of this matrix has and will be used. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::UPLO "char
Epetra_SerialSymDenseMatrix::UPLO() const

Returns character value of UPLO used by LAPACK routines. ";

/*  Mathematical Methods  */

%feature("docstring")  Epetra_SerialSymDenseMatrix::Scale "int
Epetra_SerialSymDenseMatrix::Scale(double ScalarA)

Inplace scalar-matrix product A = a A.

Scale a matrix, entry-by-entry using the value ScalarA. This method is
sensitive to the UPLO() parameter.

Parameters:
-----------

ScalarA:  (In) Scalar to multiply with A.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::NormOne "double
Epetra_SerialSymDenseMatrix::NormOne() const

Computes the 1-Norm of the this matrix.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::NormInf "double
Epetra_SerialSymDenseMatrix::NormInf() const

Computes the Infinity-Norm of the this matrix. ";

/*  Deprecated methods (will be removed in later versions of this
class)  */

%feature("docstring")  Epetra_SerialSymDenseMatrix::OneNorm "double
Epetra_SerialSymDenseMatrix::OneNorm() const

Computes the 1-Norm of the this matrix (identical to NormOne()
method).

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::InfNorm "double
Epetra_SerialSymDenseMatrix::InfNorm() const

Computes the Infinity-Norm of the this matrix (identical to NormInf()
method). ";

%feature("docstring")  Epetra_SerialSymDenseMatrix::CopyUPLOMat "void
Epetra_SerialSymDenseMatrix::CopyUPLOMat(bool Upper, double *A, int
LDA, int NumRows) ";


// File: classEpetra__SrcDistObject.xml
%feature("docstring") Epetra_SrcDistObject "

Epetra_SrcDistObject: A class for supporting flexible source
distributed objects for import/export operations.

The Epetra_SrcDistObject is a base class for all Epetra distributed
global objects that are potential source objects for the general
Epetra_DistObject class. It provides a way to send a very general
distributed object as the potential source object for an import or
export object. For example, it is possible to pass an Epetra_RowMatrix
object as the source object for an import/export where the target is
an Epetra_CrsMatrix, or an Epetra_CrsGraph (where the RowMatrix values
will be ignored).

C++ includes: Epetra_SrcDistObject.h ";

/*  Destructor  */

%feature("docstring")  Epetra_SrcDistObject::~Epetra_SrcDistObject "virtual Epetra_SrcDistObject::~Epetra_SrcDistObject()

Epetra_SrcDistObject destructor. ";

/*  Attribute accessor methods  */

%feature("docstring")  Epetra_SrcDistObject::Map "virtual const
Epetra_BlockMap& Epetra_SrcDistObject::Map() const =0

Returns a reference to the Epetra_BlockMap for this object. ";


// File: classEpetra__Time.xml
%feature("docstring") Epetra_Time "

Epetra_Time: The Epetra Timing Class.

The Epetra_Time class is a wrapper that encapsulates the general
information needed getting timing information. Currently it return the
elapsed time for each calling processor.. A Epetra_Comm object is
required for building all Epetra_Time objects.

Epetra_Time support both serial execution and (via MPI) parallel
distributed memory execution. It is meant to insulate the user from
the specifics of timing across a variety of platforms.

C++ includes: Epetra_Time.h ";

%feature("docstring")  Epetra_Time::Epetra_Time "Epetra_Time::Epetra_Time(const Epetra_Comm &Comm)

Epetra_Time Constructor.

Creates a Epetra_Time instance. This instance can be queried for
elapsed time on the calling processor. StartTime is also set for use
with the ElapsedTime function. ";

%feature("docstring")  Epetra_Time::Epetra_Time "Epetra_Time::Epetra_Time(const Epetra_Time &Time)

Epetra_Time Copy Constructor.

Makes an exact copy of an existing Epetra_Time instance. ";

%feature("docstring")  Epetra_Time::WallTime "double
Epetra_Time::WallTime(void) const

Epetra_Time wall-clock time function.

Returns the wall-clock time in seconds. A code section can be timed by
putting it between two calls to WallTime and taking the difference of
the times. ";

%feature("docstring")  Epetra_Time::ResetStartTime "void
Epetra_Time::ResetStartTime(void)

Epetra_Time function to reset the start time for a timer object.

Resets the start time for the timer object to the current time A code
section can be timed by putting it between a call to ResetStartTime
and ElapsedTime. ";

%feature("docstring")  Epetra_Time::ElapsedTime "double
Epetra_Time::ElapsedTime(void) const

Epetra_Time elapsed time function.

Returns the elapsed time in seconds since the timer object was
constructed, or since the ResetStartTime function was called. A code
section can be timed by putting it between the Epetra_Time constructor
and a call to ElapsedTime, or between a call to ResetStartTime and
ElapsedTime. ";

%feature("docstring")  Epetra_Time::~Epetra_Time "Epetra_Time::~Epetra_Time(void)

Epetra_Time Destructor.

Completely deletes a Epetra_Time object. ";


// File: classEpetra__Util.xml
%feature("docstring") Epetra_Util "

Epetra_Util: The Epetra Util Wrapper Class.

The Epetra_Util class is a collection of useful functions that cut
across a broad set of other classes. A random number generator is
provided, along with methods to set and retrieve the random-number
seed.

The random number generator is a multiplicative linear congruential
generator, with multiplier 16807 and modulus 2^31 - 1. It is based on
the algorithm described in \"Random Number Generators: Good Ones Are
Hard To Find\", S. K. Park and K. W. Miller, Communications of the
ACM, vol. 31, no. 10, pp. 1192-1201.

Sorting is provided by a static function on this class (i.e., it is
not necessary to construct an instance of this class to use the Sort
function).

A static function is provided for creating a new Epetra_Map object
with 1-to-1 ownership of entries from an existing map which may have
entries that appear on multiple processors.

Epetra_Util is a serial interface only. This is appropriate since the
standard utilities are only specified for serial execution (or shared
memory parallel).

C++ includes: Epetra_Util.h ";

/*  Random number utilities  */

%feature("docstring")  Epetra_Util::RandomInt "unsigned int
Epetra_Util::RandomInt()

Returns a random integer on the interval (0, 2^31-1) ";

%feature("docstring")  Epetra_Util::RandomDouble "double
Epetra_Util::RandomDouble()

Returns a random double on the interval (-1.0,1.0) ";

%feature("docstring")  Epetra_Util::Seed "unsigned int
Epetra_Util::Seed() const

Get seed from Random function.

Current random number seed. ";

%feature("docstring")  Epetra_Util::SetSeed "int
Epetra_Util::SetSeed(unsigned int Seed_in)

Set seed for Random function.

Parameters:
-----------

In:  Seed - An integer on the interval [1, 2^31-2]

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_Util::Epetra_Util "Epetra_Util::Epetra_Util()

Epetra_Util Constructor.

Builds an instance of a serial Util object. ";

%feature("docstring")  Epetra_Util::Epetra_Util "Epetra_Util::Epetra_Util(const Epetra_Util &Util)

Epetra_Util Copy Constructor.

Makes an exact copy of an existing Epetra_Util instance. ";

%feature("docstring")  Epetra_Util::~Epetra_Util "Epetra_Util::~Epetra_Util()

Epetra_Util Destructor. ";


// File: classEpetra__VbrMatrix.xml
%feature("docstring") Epetra_VbrMatrix "

Epetra_VbrMatrix: A class for the construction and use of real-valued
double-precision variable block-row sparse matrices.

The Epetra_VbrMatrix class is a sparse variable block row matrix
object. This matrix can be used in a parallel setting, with data
distribution described by Epetra_Map attributes. The structure or
graph of the matrix is defined by an Epetra_CrsGraph attribute.

In addition to coefficient access, the primary operations provided by
Epetra_VbrMatrix are matrix times vector and matrix times multi-vector
multiplication.

Creating and filling Epetra_VbrMatrix objects

Constructing Epetra_VbrMatrix objects is a multi-step process. The
basic steps are as follows: Create Epetra_VbrMatrix instance via one
of the constructors: Constructor that accepts one Epetra_Map object, a
row-map defining the distribution of matrix rows.

Constructor that accepts two Epetra_Map objects. (The second map is a
column-map, and describes the set of column-indices that appear in
each processor's portion of the matrix. Generally these are
overlapping sets -- column-indices may appear on more than one
processor.)

Constructor that accepts an Epetra_CrsGraph object, defining the non-
zero structure of the matrix.

Input coefficient values (more detail on this below).

Complete construction by calling FillComplete.

Note that even after FillComplete() has been called, it is possible to
update existing matrix entries but it is not possible to create new
entries.

Using Epetra_VbrMatrix as an Epetra_RowMatrix

Although Epetra_VbrMatrix does inherit from Epetra_RowMatrix, a design
flaw in the inheritance structure of Epetra prohibits the use of an
Epetra_VbrMatrix object as an Epetra_RowMatrix in some important
situations. Therefore we recommend the use of the Epetra_VbrRowMatrix
class to wrap an Epetra_VbrMatrix object for use as an
Epetra_RowMatrix. The Epetra_VbrRowMatrix object does not duplicate
data in the Epetra_VbrMatrix object, but uses it to satisfy the
Epetra_RowMatrix interface.

Epetra_Map attributes

Epetra_VbrMatrix objects have four Epetra_Map attributes, which are
held by the Epetra_CrsGraph attribute.

The Epetra_Map attributes can be obtained via these accessor methods:
RowMap() Describes the numbering and distribution of the rows of the
matrix. The row-map exists and is valid for the entire life of the
matrix. The set of matrix rows is defined by the row-map and may not
be changed. Rows may not be inserted or deleted by the user. The only
change that may be made is that the user can replace the row-map with
a compatible row-map (which is the same except for re-numbering) by
calling the ReplaceRowMap() method.

ColMap() Describes the set of column-indices that appear in the rows
in each processor's portion of the matrix. Unless provided by the user
at construction time, a valid column-map doesn't exist until
FillComplete() is called.

RangeMap() Describes the range of the matrix operator. e.g., for a
matrix-vector product operation, the result vector's map must be
compatible with the range-map of this matrix. The range-map is usually
the same as the row-map. The range-map is set equal to the row-map at
matrix creation time, but may be specified by the user when
FillComplete() is called.

DomainMap() Describes the domain of the matrix operator. The domain-
map can be specified by the user when FillComplete() is called. Until
then, it is set equal to the row-map.

It is important to note that while the row-map and the range-map are
often the same, the column-map and the domain-map are almost never the
same. The set of entries in a distributed column-map almost always
form overlapping sets, with entries being associated with more than
one processor. A domain-map, on the other hand, must be a 1-to-1 map,
with entries being associated with only a single processor.

Local versus Global Indices

Epetra_VbrMatrix has query functions IndicesAreLocal() and
IndicesAreGlobal(), which are used to determine whether the underlying
Epetra_CrsGraph attribute's column-indices have been transformed into
a local index space or not. (This transformation occurs when the
method Epetra_CrsGraph::FillComplete() is called, which happens when
the method Epetra_VbrMatrix::FillComplete() is called.) The state of
the indices in the graph determines the behavior of many
Epetra_VbrMatrix methods. If an Epetra_VbrMatrix instance is
constructed using one of the constructors that does not accept a pre-
existing Epetra_CrsGraph object, then an Epetra_CrsGraph attribute is
created internally and its indices remain untransformed (
IndicesAreGlobal()==true) until Epetra_VbrMatrix::FillComplete() is
called. The query function Epetra_VbrMatrix::Filled() returns true if
Epetra_VbrMatrix::FillComplete() has been called.

Inputting coefficient values

The process for inputting block-entry coefficients is as follows:
Indicate that values for a specified row are about to be provided by
calling one of these methods which specify a block-row and a list of
block-column-indices:  BeginInsertGlobalValues()

BeginInsertMyValues()

BeginReplaceGlobalValues()

BeginReplaceMyValues()

BeginSumIntoGlobalValues()

BeginSumIntoMyValues()

Loop over the list of block-column-indices and pass each block-entry
to the matrix using the method SubmitBlockEntry().

Complete the process for the specified block-row by calling the method
EndSubmitEntries().

Note that the 'GlobalValues' methods have the precondition that
IndicesAreGlobal() must be true, and the 'MyValues' methods have the
precondition that IndicesAreLocal() must be true. Furthermore, the
'SumInto' and 'Replace' methods may only be used to update matrix
entries which already exist, and the 'Insert' methods may only be used
if IndicesAreContiguous() is false.

Counting Floating Point Operations

Each Epetra_VbrMatrix object keeps track of the number of serial
floating point operations performed using the specified object as the
this argument to the function. The Flops() function returns this
number as a double precision number. Using this information, in
conjunction with the Epetra_Time class, one can get accurate parallel
performance numbers. The ResetFlops() function resets the floating
point counter.

C++ includes: Epetra_VbrMatrix.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_VbrMatrix::Epetra_VbrMatrix "Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int *NumBlockEntriesPerRow)

Epetra_VbrMatrix constuctor with variable number of indices per row.

Creates a Epetra_VbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - A Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  NumBlockEntriesPerRow - An integer array of length NumRows such
that NumBlockEntriesPerRow[i] indicates the (approximate) number of
Block entries in the ith row. ";

%feature("docstring")  Epetra_VbrMatrix::Epetra_VbrMatrix "Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int NumBlockEntriesPerRow)

Epetra_VbrMatrix constuctor with fixed number of indices per row.

Creates a Epetra_VbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - An Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  NumBlockEntriesPerRow - An integer that indicates the
(approximate) number of Block entries in the each Block row. Note that
it is possible to use 0 for this value and let fill occur during the
insertion phase. ";

%feature("docstring")  Epetra_VbrMatrix::Epetra_VbrMatrix "Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
*NumBlockEntriesPerRow)

Epetra_VbrMatrix constuctor with variable number of indices per row.

Creates a Epetra_VbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - A Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  ColMap - A Epetra_BlockMap.

In:  NumBlockEntriesPerRow - An integer array of length NumRows such
that NumBlockEntriesPerRow[i] indicates the (approximate) number of
Block entries in the ith row. ";

%feature("docstring")  Epetra_VbrMatrix::Epetra_VbrMatrix "Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
NumBlockEntriesPerRow)

Epetra_VbrMatrix constuctor with fixed number of indices per row.

Creates a Epetra_VbrMatrix object and allocates storage.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  RowMap - A Epetra_BlockMap listing the block rows that this
processor will contribute to.

In:  ColMap - An Epetra_BlockMap listing the block columns that this
processor will contribute to.

In:  NumBlockEntriesPerRow - An integer that indicates the
(approximate) number of Block entries in the each Block row. Note that
it is possible to use 0 for this value and let fill occur during the
insertion phase. ";

%feature("docstring")  Epetra_VbrMatrix::Epetra_VbrMatrix "Epetra_VbrMatrix::Epetra_VbrMatrix(Epetra_DataAccess CV, const
Epetra_CrsGraph &Graph)

Construct a matrix using an existing Epetra_CrsGraph object.

Allows the nonzero structure from another matrix, or a structure that
was constructed independently, to be used for this matrix.

Parameters:
-----------

In:  CV - A Epetra_DataAccess enumerated type set to Copy or View.

In:  Graph - A Epetra_CrsGraph object, extracted from another Epetra
matrix object or constructed directly from using the Epetra_CrsGraph
constructors. ";

%feature("docstring")  Epetra_VbrMatrix::Epetra_VbrMatrix "Epetra_VbrMatrix::Epetra_VbrMatrix(const Epetra_VbrMatrix &Matrix)

Copy constructor. ";

%feature("docstring")  Epetra_VbrMatrix::~Epetra_VbrMatrix "Epetra_VbrMatrix::~Epetra_VbrMatrix()

Epetra_VbrMatrix Destructor. ";

/*  Insertion/Replace/SumInto methods  */

%feature("docstring")  Epetra_VbrMatrix::PutScalar "int
Epetra_VbrMatrix::PutScalar(double ScalarConstant)

Initialize all values in graph of the matrix with constant value.

Parameters:
-----------

In:  ScalarConstant - Value to use.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::Scale "int
Epetra_VbrMatrix::Scale(double ScalarConstant)

Multiply all values in the matrix by a constant value (in place: A <-
ScalarConstant * A).

Parameters:
-----------

In:  ScalarConstant - Value to use.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::DirectSubmitBlockEntry "int
Epetra_VbrMatrix::DirectSubmitBlockEntry(int GlobalBlockRow, int
GlobalBlockCol, const double *values, int LDA, int NumRows, int
NumCols, bool sum_into)

Submit a block-entry directly into the matrix (without using a
begin/end sequence)

Experimental method which allows submitting a block-entry without
first calling BeginInsertGlobalValues. This method copies the input
data directly into the matrix storage. The block-entry is specified by
global block-row and block-col indices. ";

%feature("docstring")  Epetra_VbrMatrix::BeginInsertGlobalValues "int
Epetra_VbrMatrix::BeginInsertGlobalValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate insertion of a list of elements in a given global row of the
matrix, values are inserted via SubmitEntry().

Parameters:
-----------

In:  BlockRow - Block Row number (in global coordinates) to put
elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginInsertMyValues "int
Epetra_VbrMatrix::BeginInsertMyValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate insertion of a list of elements in a given local row of the
matrix, values are inserted via SubmitEntry().

Parameters:
-----------

In:  BlockRow - Block Row number (in local coordinates) to put
elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Local column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginReplaceGlobalValues "int Epetra_VbrMatrix::BeginReplaceGlobalValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate replacement of current values with this list of entries for a
given global row of the matrix, values are replaced via SubmitEntry()

Parameters:
-----------

In:  Row - Block Row number (in global coordinates) to put elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginReplaceMyValues "int
Epetra_VbrMatrix::BeginReplaceMyValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate replacement of current values with this list of entries for a
given local row of the matrix, values are replaced via SubmitEntry()

Parameters:
-----------

In:  Row - Block Row number (in local coordinates) to put elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Local column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginSumIntoGlobalValues "int Epetra_VbrMatrix::BeginSumIntoGlobalValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate summing into current values with this list of entries for a
given global row of the matrix, values are replaced via SubmitEntry()

Parameters:
-----------

In:  Row - Block Row number (in global coordinates) to put elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Global column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginSumIntoMyValues "int
Epetra_VbrMatrix::BeginSumIntoMyValues(int BlockRow, int
NumBlockEntries, int *BlockIndices)

Initiate summing into current values with this list of entries for a
given local row of the matrix, values are replaced via SubmitEntry()

Parameters:
-----------

In:  Row - Block Row number (in local coordinates) to put elements.

In:  NumBlockEntries - Number of entries.

In:  Indices - Local column indices corresponding to values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::SubmitBlockEntry "int
Epetra_VbrMatrix::SubmitBlockEntry(double *Values, int LDA, int
NumRows, int NumCols)

Submit a block entry to the indicated block row and column specified
in the Begin routine. ";

%feature("docstring")  Epetra_VbrMatrix::SubmitBlockEntry "int
Epetra_VbrMatrix::SubmitBlockEntry(Epetra_SerialDenseMatrix &Mat)

Submit a block entry to the indicated block row and column specified
in the Begin routine. ";

%feature("docstring")  Epetra_VbrMatrix::EndSubmitEntries "int
Epetra_VbrMatrix::EndSubmitEntries()

Completes processing of all data passed in for the current block row.

This function completes the processing of all block entries submitted
via SubmitBlockEntry(). It also checks to make sure that
SubmitBlockEntry was called the correct number of times as specified
by the Begin routine that initiated the entry process. ";

%feature("docstring")  Epetra_VbrMatrix::ReplaceDiagonalValues "int
Epetra_VbrMatrix::ReplaceDiagonalValues(const Epetra_Vector &Diagonal)

Replaces diagonal values of the with those in the user-provided
vector.

This routine is meant to allow replacement of { existing} diagonal
values. If a diagonal value does not exist for a given row, the
corresponding value in the input Epetra_Vector will be ignored and the
return code will be set to 1.

The Epetra_Map associated with the input Epetra_Vector must be
compatible with the RowMap of the matrix.

Parameters:
-----------

Diagonal:  (In) - New values to be placed in the main diagonal.

Integer error code, set to 0 if successful, 1 of one or more diagonal
entries not present in matrix. ";

%feature("docstring")  Epetra_VbrMatrix::FillComplete "int
Epetra_VbrMatrix::FillComplete()

Signal that data entry is complete, perform transformations to local
index space. ";

%feature("docstring")  Epetra_VbrMatrix::FillComplete "int
Epetra_VbrMatrix::FillComplete(const Epetra_BlockMap &DomainMap, const
Epetra_BlockMap &RangeMap)

Signal that data entry is complete, perform transformations to local
index space. ";

%feature("docstring")  Epetra_VbrMatrix::Filled "bool
Epetra_VbrMatrix::Filled() const

If FillComplete() has been called, this query returns true, otherwise
it returns false. ";

/*  Extraction methods  */

%feature("docstring")  Epetra_VbrMatrix::ExtractGlobalBlockRowPointers
"int Epetra_VbrMatrix::ExtractGlobalBlockRowPointers(int BlockRow,
int MaxNumBlockEntries, int &RowDim, int &NumBlockEntries, int
*BlockIndices, Epetra_SerialDenseMatrix **&Values) const

Copy the block indices into user-provided array, set pointers for rest
of data for specified global block row.

This function provides the lightest weight approach to accessing a
global block row when the matrix may be be stored in local or global
index space. In other words, this function will always work because
the block indices are returned in user-provided space. All other array
arguments are independent of whether or not indices are local or
global. Other than the BlockIndices array, all other array argument
are returned as pointers to internal data.

Parameters:
-----------

In:  BlockRow - Global block row to extract.

In:  MaxNumBlockEntries - Length of user-provided BlockIndices array.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries actually extracted.

Out:  BlockIndices - Extracted global column indices for the
corresponding block entries.

Out:  Values - Pointer to list of pointers to block entries. Note that
the actual values are not copied.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractMyBlockRowPointers "int Epetra_VbrMatrix::ExtractMyBlockRowPointers(int BlockRow, int
MaxNumBlockEntries, int &RowDim, int &NumBlockEntries, int
*BlockIndices, Epetra_SerialDenseMatrix **&Values) const

Copy the block indices into user-provided array, set pointers for rest
of data for specified local block row.

This function provides the lightest weight approach to accessing a
local block row when the matrix may be be stored in local or global
index space. In other words, this function will always work because
the block indices are returned in user-provided space. All other array
arguments are independent of whether or not indices are local or
global. Other than the BlockIndices array, all other array argument
are returned as pointers to internal data.

Parameters:
-----------

In:  BlockRow - Local block row to extract.

In:  MaxNumBlockEntries - Length of user-provided BlockIndices array.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries actually extracted.

Out:  BlockIndices - Extracted local column indices for the
corresponding block entries.

Out:  Values - Pointer to list of pointers to block entries. Note that
the actual values are not copied.

Integer error code, set to 0 if successful. ";

%feature("docstring")
Epetra_VbrMatrix::BeginExtractGlobalBlockRowCopy "int
Epetra_VbrMatrix::BeginExtractGlobalBlockRowCopy(int BlockRow, int
MaxNumBlockEntries, int &RowDim, int &NumBlockEntries, int
*BlockIndices, int *ColDims) const

Initiates a copy of the specified global row in user-provided arrays.

Parameters:
-----------

In:  BlockRow - Global block row to extract.

In:  MaxNumBlockEntries - Length of user-provided BlockIndices,
ColDims, and LDAs arrays.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries actually extracted.

Out:  BlockIndices - Extracted global column indices for the
corresponding block entries.

Out:  ColDim - List of column dimensions for each corresponding block
entry that will be extracted.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginExtractMyBlockRowCopy "int Epetra_VbrMatrix::BeginExtractMyBlockRowCopy(int BlockRow, int
MaxNumBlockEntries, int &RowDim, int &NumBlockEntries, int
*BlockIndices, int *ColDims) const

Initiates a copy of the specified local row in user-provided arrays.

Parameters:
-----------

In:  BlockRow - Local block row to extract.

In:  MaxNumBlockEntries - Length of user-provided BlockIndices,
ColDims, and LDAs arrays.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries actually extracted.

Out:  BlockIndices - Extracted local column indices for the
corresponding block entries.

Out:  ColDim - List of column dimensions for each corresponding block
entry that will be extracted.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractEntryCopy "int
Epetra_VbrMatrix::ExtractEntryCopy(int SizeOfValues, double *Values,
int LDA, bool SumInto) const

Extract a copy of an entry from the block row specified by one of the
BeginExtract routines.

Once BeginExtractGlobalBlockRowCopy() or BeginExtractMyBlockRowCopy()
is called, you can extract the block entries of specified block row
one-entry-at-a-time. The entries will be extracted in an order
corresponding to the BlockIndices list that was returned by the
BeginExtract routine.

Parameters:
-----------

In:  SizeOfValues - Amount of memory associated with Values. This must
be at least as big as LDA*NumCol, where NumCol is the column dimension
of the block entry being copied

InOut:  Values - Starting location where the block entry will be
copied.

In:  LDA - Specifies the stride that will be used when copying columns
into Values.

In:  SumInto - If set to true, the block entry values will be summed
into existing values. ";

%feature("docstring")
Epetra_VbrMatrix::BeginExtractGlobalBlockRowView "int
Epetra_VbrMatrix::BeginExtractGlobalBlockRowView(int BlockRow, int
&RowDim, int &NumBlockEntries, int *&BlockIndices) const

Initiates a view of the specified global row, only works if matrix
indices are in global mode.

Parameters:
-----------

In:  BlockRow - Global block row to view.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries to be viewed.

Out:  BlockIndices - Pointer to global column indices for the
corresponding block entries.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginExtractMyBlockRowView "int Epetra_VbrMatrix::BeginExtractMyBlockRowView(int BlockRow, int
&RowDim, int &NumBlockEntries, int *&BlockIndices) const

Initiates a view of the specified local row, only works if matrix
indices are in local mode.

Parameters:
-----------

In:  BlockRow - Local block row to view.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries to be viewed.

Out:  BlockIndices - Pointer to local column indices for the
corresponding block entries.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractEntryView "int
Epetra_VbrMatrix::ExtractEntryView(Epetra_SerialDenseMatrix *&entry)
const

Returns a pointer to the current block entry.

After a call to BeginExtractGlobal() or
BlockRowViewBeginExtractMyBlockRowView(), ExtractEntryView() can be
called up to NumBlockEntries times to get each block entry in the
specified block row.

Parameters:
-----------

InOut:  entry - A pointer that will be set to the current block entry.
";

%feature("docstring")  Epetra_VbrMatrix::ExtractGlobalBlockRowView "int Epetra_VbrMatrix::ExtractGlobalBlockRowView(int BlockRow, int
&RowDim, int &NumBlockEntries, int *&BlockIndices,
Epetra_SerialDenseMatrix **&Values) const

Initiates a view of the specified global row, only works if matrix
indices are in global mode.

Parameters:
-----------

In:  BlockRow - Global block row to view.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries to be viewed.

Out:  BlockIndices - Pointer to global column indices for the
corresponding block entries.

Out:  Values - Pointer to an array of pointers to the block entries in
the specified block row.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractMyBlockRowView "int
Epetra_VbrMatrix::ExtractMyBlockRowView(int BlockRow, int &RowDim, int
&NumBlockEntries, int *&BlockIndices, Epetra_SerialDenseMatrix
**&Values) const

Initiates a view of the specified local row, only works if matrix
indices are in local mode.

Parameters:
-----------

In:  BlockRow - Local block row to view.

Out:  RowDim - Number of equations in the requested block row.

Out:  NumBlockEntries - Number of nonzero entries to be viewed.

Out:  BlockIndices - Pointer to local column indices for the
corresponding block entries.

Out:  Values - Pointer to an array of pointers to the block entries in
the specified block row.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractDiagonalCopy "int
Epetra_VbrMatrix::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const

Returns a copy of the main diagonal in a user-provided vector.

Parameters:
-----------

Out:  Diagonal - Extracted main diagonal.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::BeginExtractBlockDiagonalCopy
"int Epetra_VbrMatrix::BeginExtractBlockDiagonalCopy(int
MaxNumBlockDiagonalEntries, int &NumBlockDiagonalEntries, int
*RowColDims) const

Initiates a copy of the block diagonal entries to user-provided
arrays.

Parameters:
-----------

In:  MaxNumBlockDiagonalEntries - Length of user-provided RowColDims
array.

Out:  NumBlockDiagonalEntries - Number of block diagonal entries that
can actually be extracted.

Out:  RowColDim - List of row and column dimension for corresponding
block diagonal entries.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractBlockDiagonalEntryCopy
"int Epetra_VbrMatrix::ExtractBlockDiagonalEntryCopy(int
SizeOfValues, double *Values, int LDA, bool SumInto) const

Extract a copy of a block diagonal entry from the matrix.

Once BeginExtractBlockDiagonalCopy() is called, you can extract the
block diagonal entries one-entry- at-a-time. The entries will be
extracted in ascending order.

Parameters:
-----------

In:  SizeOfValues - Amount of memory associated with Values. This must
be at least as big as LDA*NumCol, where NumCol is the column dimension
of the block entry being copied

InOut:  Values - Starting location where the block entry will be
copied.

In:  LDA - Specifies the stride that will be used when copying columns
into Values.

In:  SumInto - If set to true, the block entry values will be summed
into existing values. ";

%feature("docstring")  Epetra_VbrMatrix::BeginExtractBlockDiagonalView
"int Epetra_VbrMatrix::BeginExtractBlockDiagonalView(int
&NumBlockDiagonalEntries, int *&RowColDims) const

Initiates a view of the block diagonal entries.

Parameters:
-----------

Out:  NumBlockDiagonalEntries - Number of block diagonal entries that
can be viewed.

Out:  RowColDim - Pointer to list of row and column dimension for
corresponding block diagonal entries.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractBlockDiagonalEntryView
"int Epetra_VbrMatrix::ExtractBlockDiagonalEntryView(double *&Values,
int &LDA) const

Extract a view of a block diagonal entry from the matrix.

Once BeginExtractBlockDiagonalView() is called, you can extract a view
of the block diagonal entries one- entry-at-a-time. The views will be
extracted in ascending order.

Parameters:
-----------

Out:  Values - Pointer to internal copy of block entry.

Out:  LDA - Column stride of Values. ";

/*  Computational methods  */

%feature("docstring")  Epetra_VbrMatrix::Multiply1 "int
Epetra_VbrMatrix::Multiply1(bool TransA, const Epetra_Vector &x,
Epetra_Vector &y) const

Returns the result of a Epetra_VbrMatrix multiplied by a Epetra_Vector
x in y.

Parameters:
-----------

In:  TransA - If true, multiply by the transpose of matrix, otherwise
just use matrix.

In:  x - A Epetra_Vector to multiply by.

Out:  y - A Epetra_Vector containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::Multiply "int
Epetra_VbrMatrix::Multiply(bool TransA, const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_VbrMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  TransA -If true, multiply by the transpose of matrix, otherwise
just use matrix.

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::Solve "int
Epetra_VbrMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_Vector &x, Epetra_Vector &y) const

Returns the result of a solve using the Epetra_VbrMatrix on a
Epetra_Vector x in y.

Parameters:
-----------

In:  Upper -If true, solve Ux = y, otherwise solve Lx = y.

In:  Trans -If true, solve transpose problem.

In:  UnitDiagonal -If true, assume diagonal is unit (whether it's
stored or not).

In:  x -A Epetra_Vector to solve for.

Out:  y -A Epetra_Vector containing result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::Solve "int
Epetra_VbrMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_VbrMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  Upper -If true, solve Ux = y, otherwise solve Lx = y.

In:  Trans -If true, solve transpose problem.

In:  UnitDiagonal -If true, assume diagonal is unit (whether it's
stored or not).

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::InvRowSums "int
Epetra_VbrMatrix::InvRowSums(Epetra_Vector &x) const

Computes the sum of absolute values of the rows of the
Epetra_VbrMatrix, results returned in x.

The vector x will return such that x[i] will contain the inverse of
sum of the absolute values of the this matrix will be scaled such that
A(i,j) = x(i)*A(i,j) where i denotes the global row number of A and j
denotes the global column number of A. Using the resulting vector from
this function as input to LeftScale() will make the infinity norm of
the resulting matrix exactly 1.

Parameters:
-----------

Out:  x -A Epetra_Vector containing the row sums of the this matrix.

WARNING:  It is assumed that the distribution of x is the same as the
rows of this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::LeftScale "int
Epetra_VbrMatrix::LeftScale(const Epetra_Vector &x)

Scales the Epetra_VbrMatrix on the left with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i
denotes the row number of A and j denotes the column number of A.

Parameters:
-----------

In:  x -A Epetra_Vector to solve for.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::InvColSums "int
Epetra_VbrMatrix::InvColSums(Epetra_Vector &x) const

Computes the sum of absolute values of the columns of the
Epetra_VbrMatrix, results returned in x.

The vector x will return such that x[j] will contain the inverse of
sum of the absolute values of the this matrix will be sca such that
A(i,j) = x(j)*A(i,j) where i denotes the global row number of A and j
denotes the global column number of A. Using the resulting vector from
this function as input to RighttScale() will make the one norm of the
resulting matrix exactly 1.

Parameters:
-----------

Out:  x -A Epetra_Vector containing the column sums of the this
matrix.

WARNING:  It is assumed that the distribution of x is the same as the
rows of this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::RightScale "int
Epetra_VbrMatrix::RightScale(const Epetra_Vector &x)

Scales the Epetra_VbrMatrix on the right with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i
denotes the global row number of A and j denotes the global column
number of A.

Parameters:
-----------

In:  x -The Epetra_Vector used for scaling this.

Integer error code, set to 0 if successful. ";

/*  Matrix Properties Query Methods  */

%feature("docstring")  Epetra_VbrMatrix::OptimizeStorage "int
Epetra_VbrMatrix::OptimizeStorage()

Eliminates memory that is used for construction. Make consecutive row
index sections contiguous. ";

%feature("docstring")  Epetra_VbrMatrix::StorageOptimized "bool
Epetra_VbrMatrix::StorageOptimized() const

If OptimizeStorage() has been called, this query returns true,
otherwise it returns false. ";

%feature("docstring")  Epetra_VbrMatrix::IndicesAreGlobal "bool
Epetra_VbrMatrix::IndicesAreGlobal() const

If matrix indices has not been transformed to local, this query
returns true, otherwise it returns false. ";

%feature("docstring")  Epetra_VbrMatrix::IndicesAreLocal "bool
Epetra_VbrMatrix::IndicesAreLocal() const

If matrix indices has been transformed to local, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_VbrMatrix::IndicesAreContiguous "bool
Epetra_VbrMatrix::IndicesAreContiguous() const

If matrix indices are packed into single array (done in
OptimizeStorage()) return true, otherwise false. ";

%feature("docstring")  Epetra_VbrMatrix::LowerTriangular "bool
Epetra_VbrMatrix::LowerTriangular() const

If matrix is lower triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_VbrMatrix::UpperTriangular "bool
Epetra_VbrMatrix::UpperTriangular() const

If matrix is upper triangular in local index space, this query returns
true, otherwise it returns false. ";

%feature("docstring")  Epetra_VbrMatrix::NoDiagonal "bool
Epetra_VbrMatrix::NoDiagonal() const

If matrix has no diagonal entries based on global row/column index
comparisons, this query returns true, otherwise it returns false. ";

/*  Attribute access functions  */

%feature("docstring")  Epetra_VbrMatrix::NormInf "double
Epetra_VbrMatrix::NormInf() const

Returns the infinity norm of the global matrix. ";

%feature("docstring")  Epetra_VbrMatrix::NormOne "double
Epetra_VbrMatrix::NormOne() const

Returns the one norm of the global matrix. ";

%feature("docstring")  Epetra_VbrMatrix::NormFrobenius "double
Epetra_VbrMatrix::NormFrobenius() const

Returns the frobenius norm of the global matrix. ";

%feature("docstring")  Epetra_VbrMatrix::MaxRowDim "int
Epetra_VbrMatrix::MaxRowDim() const

Returns the maximum row dimension of all block entries on this
processor. ";

%feature("docstring")  Epetra_VbrMatrix::MaxColDim "int
Epetra_VbrMatrix::MaxColDim() const

Returns the maximum column dimension of all block entries on this
processor. ";

%feature("docstring")  Epetra_VbrMatrix::GlobalMaxRowDim "int
Epetra_VbrMatrix::GlobalMaxRowDim() const

Returns the maximum row dimension of all block entries across all
processors. ";

%feature("docstring")  Epetra_VbrMatrix::GlobalMaxColDim "int
Epetra_VbrMatrix::GlobalMaxColDim() const

Returns the maximum column dimension of all block entries across all
processors. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyRows "int
Epetra_VbrMatrix::NumMyRows() const

Returns the number of matrix rows owned by the calling processor. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyCols "int
Epetra_VbrMatrix::NumMyCols() const

Returns the number of matrix columns owned by the calling processor.
";

%feature("docstring")  Epetra_VbrMatrix::NumMyNonzeros "int
Epetra_VbrMatrix::NumMyNonzeros() const

Returns the number of nonzero entriesowned by the calling processor .
";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalRows "int
Epetra_VbrMatrix::NumGlobalRows() const

Returns the number of global matrix rows. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalRows64 "long long
Epetra_VbrMatrix::NumGlobalRows64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalCols "int
Epetra_VbrMatrix::NumGlobalCols() const

Returns the number of global matrix columns. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalCols64 "long long
Epetra_VbrMatrix::NumGlobalCols64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalNonzeros "int
Epetra_VbrMatrix::NumGlobalNonzeros() const

Returns the number of nonzero entries in the global matrix. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalNonzeros64 "long
long Epetra_VbrMatrix::NumGlobalNonzeros64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumMyBlockRows "int
Epetra_VbrMatrix::NumMyBlockRows() const

Returns the number of Block matrix rows owned by the calling
processor. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyBlockCols "int
Epetra_VbrMatrix::NumMyBlockCols() const

Returns the number of Block matrix columns owned by the calling
processor. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyBlockEntries "int
Epetra_VbrMatrix::NumMyBlockEntries() const

Returns the number of nonzero block entries in the calling processor's
portion of the matrix. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyBlockDiagonals "int
Epetra_VbrMatrix::NumMyBlockDiagonals() const

Returns the number of local nonzero block diagonal entries, based on
global row/column index comparisons. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyDiagonals "int
Epetra_VbrMatrix::NumMyDiagonals() const

Returns the number of local nonzero diagonal entries, based on global
row/column index comparisons. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockRows "int
Epetra_VbrMatrix::NumGlobalBlockRows() const

Returns the number of global Block matrix rows. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockRows64 "long
long Epetra_VbrMatrix::NumGlobalBlockRows64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockCols "int
Epetra_VbrMatrix::NumGlobalBlockCols() const

Returns the number of global Block matrix columns. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockCols64 "long
long Epetra_VbrMatrix::NumGlobalBlockCols64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockEntries "int
Epetra_VbrMatrix::NumGlobalBlockEntries() const

Returns the number of nonzero block entries in the global matrix. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockEntries64 "long long Epetra_VbrMatrix::NumGlobalBlockEntries64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockDiagonals "int
Epetra_VbrMatrix::NumGlobalBlockDiagonals() const

Returns the number of global nonzero block diagonal entries, based on
global row/column index comparisions. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockDiagonals64 "long long Epetra_VbrMatrix::NumGlobalBlockDiagonals64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalDiagonals "int
Epetra_VbrMatrix::NumGlobalDiagonals() const

Returns the number of global nonzero diagonal entries, based on global
row/column index comparisions. ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalDiagonals64 "long
long Epetra_VbrMatrix::NumGlobalDiagonals64() const ";

%feature("docstring")  Epetra_VbrMatrix::NumGlobalBlockEntries "int
Epetra_VbrMatrix::NumGlobalBlockEntries(int Row) const

Returns the current number of nonzero Block entries in specified
global row on this processor. ";

%feature("docstring")
Epetra_VbrMatrix::NumAllocatedGlobalBlockEntries "int
Epetra_VbrMatrix::NumAllocatedGlobalBlockEntries(int Row) const

Returns the allocated number of nonzero Block entries in specified
global row on this processor. ";

%feature("docstring")  Epetra_VbrMatrix::MaxNumBlockEntries "int
Epetra_VbrMatrix::MaxNumBlockEntries() const

Returns the maximum number of nonzero entries across all rows on this
processor. ";

%feature("docstring")  Epetra_VbrMatrix::GlobalMaxNumBlockEntries "int Epetra_VbrMatrix::GlobalMaxNumBlockEntries() const

Returns the maximum number of nonzero entries across all rows on this
processor. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyBlockEntries "int
Epetra_VbrMatrix::NumMyBlockEntries(int Row) const

Returns the current number of nonzero Block entries in specified local
row on this processor. ";

%feature("docstring")  Epetra_VbrMatrix::NumAllocatedMyBlockEntries "int Epetra_VbrMatrix::NumAllocatedMyBlockEntries(int Row) const

Returns the allocated number of nonzero Block entries in specified
local row on this processor. ";

%feature("docstring")  Epetra_VbrMatrix::MaxNumNonzeros "int
Epetra_VbrMatrix::MaxNumNonzeros() const

Returns the maximum number of nonzero entries across all block rows on
this processor.

Let ki = the number of nonzero values in the ith block row of the
VbrMatrix object. For example, if the ith block row had 5 block
entries and the size of each entry was 4-by-4, ki would be 80. Then
this function return the max over all ki for all row on this
processor. ";

%feature("docstring")  Epetra_VbrMatrix::GlobalMaxNumNonzeros "int
Epetra_VbrMatrix::GlobalMaxNumNonzeros() const

Returns the maximum number of nonzero entries across all block rows on
all processors.

This function returns the max over all processor of MaxNumNonzeros().
";

%feature("docstring")  Epetra_VbrMatrix::IndexBase "int
Epetra_VbrMatrix::IndexBase() const

Returns the index base for row and column indices for this graph. ";

%feature("docstring")  Epetra_VbrMatrix::Graph "const
Epetra_CrsGraph& Epetra_VbrMatrix::Graph() const

Returns a pointer to the Epetra_CrsGraph object associated with this
matrix. ";

%feature("docstring")  Epetra_VbrMatrix::Importer "const
Epetra_Import* Epetra_VbrMatrix::Importer() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

%feature("docstring")  Epetra_VbrMatrix::Exporter "const
Epetra_Export* Epetra_VbrMatrix::Exporter() const

Returns the Epetra_Export object that contains the export operations
for distributed operations. ";

%feature("docstring")  Epetra_VbrMatrix::DomainMap "const
Epetra_BlockMap& Epetra_VbrMatrix::DomainMap() const

Returns the Epetra_BlockMap object associated with the domain of this
matrix operator. ";

%feature("docstring")  Epetra_VbrMatrix::RangeMap "const
Epetra_BlockMap& Epetra_VbrMatrix::RangeMap() const

Returns the Epetra_BlockMap object associated with the range of this
matrix operator. ";

%feature("docstring")  Epetra_VbrMatrix::RowMap "const
Epetra_BlockMap& Epetra_VbrMatrix::RowMap() const

Returns the RowMap object as an Epetra_BlockMap (the Epetra_Map base
class) needed for implementing Epetra_RowMatrix. ";

%feature("docstring")  Epetra_VbrMatrix::ColMap "const
Epetra_BlockMap& Epetra_VbrMatrix::ColMap() const

Returns the ColMap as an Epetra_BlockMap (the Epetra_Map base class)
needed for implementing Epetra_RowMatrix. ";

%feature("docstring")  Epetra_VbrMatrix::Comm "const Epetra_Comm&
Epetra_VbrMatrix::Comm() const

Fills a matrix with rows from a source matrix based on the specified
importer.

Returns a pointer to the Epetra_Comm communicator associated with this
matrix. ";

/*  Local/Global ID methods  */

%feature("docstring")  Epetra_VbrMatrix::LRID "int
Epetra_VbrMatrix::LRID(int GRID_in) const

Returns the local row index for given global row index, returns -1 if
no local row for this global row. ";

%feature("docstring")  Epetra_VbrMatrix::LRID "int
Epetra_VbrMatrix::LRID(long long GRID_in) const ";

%feature("docstring")  Epetra_VbrMatrix::GRID "int
Epetra_VbrMatrix::GRID(int LRID_in) const

Returns the global row index for give local row index, returns
IndexBase-1 if we don't have this local row. ";

%feature("docstring")  Epetra_VbrMatrix::GRID64 "long long
Epetra_VbrMatrix::GRID64(int LRID_in) const ";

%feature("docstring")  Epetra_VbrMatrix::LCID "int
Epetra_VbrMatrix::LCID(int GCID_in) const

Returns the local column index for given global column index, returns
-1 if no local column for this global column. ";

%feature("docstring")  Epetra_VbrMatrix::LCID "int
Epetra_VbrMatrix::LCID(long long GCID_in) const ";

%feature("docstring")  Epetra_VbrMatrix::GCID "int
Epetra_VbrMatrix::GCID(int LCID_in) const

Returns the global column index for give local column index, returns
IndexBase-1 if we don't have this local column. ";

%feature("docstring")  Epetra_VbrMatrix::GCID64 "long long
Epetra_VbrMatrix::GCID64(int LCID_in) const ";

%feature("docstring")  Epetra_VbrMatrix::MyGRID "bool
Epetra_VbrMatrix::MyGRID(int GRID_in) const

Returns true if the GRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_VbrMatrix::MyGRID "bool
Epetra_VbrMatrix::MyGRID(long long GRID_in) const ";

%feature("docstring")  Epetra_VbrMatrix::MyLRID "bool
Epetra_VbrMatrix::MyLRID(int LRID_in) const

Returns true if the LRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_VbrMatrix::MyGCID "bool
Epetra_VbrMatrix::MyGCID(int GCID_in) const

Returns true if the GCID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_VbrMatrix::MyGCID "bool
Epetra_VbrMatrix::MyGCID(long long GCID_in) const ";

%feature("docstring")  Epetra_VbrMatrix::MyLCID "bool
Epetra_VbrMatrix::MyLCID(int LCID_in) const

Returns true if the LRID passed in belongs to the calling processor in
this map, otherwise returns false. ";

%feature("docstring")  Epetra_VbrMatrix::MyGlobalBlockRow "bool
Epetra_VbrMatrix::MyGlobalBlockRow(int GID) const

Returns true of GID is owned by the calling processor, otherwise it
returns false. ";

%feature("docstring")  Epetra_VbrMatrix::MyGlobalBlockRow "bool
Epetra_VbrMatrix::MyGlobalBlockRow(long long GID) const ";

/*  I/O Methods  */

%feature("docstring")  Epetra_VbrMatrix::Print "void
Epetra_VbrMatrix::Print(ostream &os) const

Print method. ";

/*  Additional methods required to support the Epetra_Operator
interface  */

%feature("docstring")  Epetra_VbrMatrix::Label "const char*
Epetra_VbrMatrix::Label() const

Returns a character string describing the operator. ";

%feature("docstring")  Epetra_VbrMatrix::SetUseTranspose "int
Epetra_VbrMatrix::SetUseTranspose(bool UseTranspose_in)

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

%feature("docstring")  Epetra_VbrMatrix::Apply "int
Epetra_VbrMatrix::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator applied to a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ApplyInverse "int
Epetra_VbrMatrix::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const

Returns the result of a Epetra_Operator inverse applied to an
Epetra_MultiVector X in Y.

In this implementation, we use several existing attributes to
determine how virtual method ApplyInverse() should call the concrete
method Solve(). We pass in the UpperTriangular(), the
Epetra_VbrMatrix::UseTranspose(), and NoDiagonal() methods. The most
notable warning is that if a matrix has no diagonal values we assume
that there is an implicit unit diagonal that should be accounted for
when doing a triangular solve.

Parameters:
-----------

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::HasNormInf "bool
Epetra_VbrMatrix::HasNormInf() const

Returns true because this class can compute an Inf-norm. ";

%feature("docstring")  Epetra_VbrMatrix::UseTranspose "bool
Epetra_VbrMatrix::UseTranspose() const

Returns the current UseTranspose setting. ";

%feature("docstring")  Epetra_VbrMatrix::OperatorDomainMap "const
Epetra_Map& Epetra_VbrMatrix::OperatorDomainMap() const

Returns the Epetra_Map object associated with the domain of this
matrix operator. ";

%feature("docstring")  Epetra_VbrMatrix::OperatorRangeMap "const
Epetra_Map& Epetra_VbrMatrix::OperatorRangeMap() const

Returns the Epetra_Map object associated with the range of this matrix
operator. ";

/*  Additional methods required to implement RowMatrix interface  */

%feature("docstring")  Epetra_VbrMatrix::ExtractGlobalRowCopy "int
Epetra_VbrMatrix::ExtractGlobalRowCopy(int GlobalRow, int Length, int
&NumEntries, double *Values, int *Indices) const

Returns a copy of the specified global row in user-provided arrays.

Parameters:
-----------

In:  GlobalRow - Global row to extract.

In:  Length - Length of Values and Indices.

Out:  NumEntries - Number of nonzero entries extracted.

Out:  Values - Extracted values for this row.

Out:  Indices - Extracted global column indices for the corresponding
values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::ExtractMyRowCopy "int
Epetra_VbrMatrix::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

In:  MyRow - Local row to extract.

In:  Length - Length of Values and Indices.

Out:  NumEntries - Number of nonzero entries extracted.

Out:  Values - Extracted values for this row.

Out:  Indices - Extracted local column indices for the corresponding
values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::NumMyRowEntries "int
Epetra_VbrMatrix::NumMyRowEntries(int MyRow, int &NumEntries) const

Return the current number of values stored for the specified local
row.

Parameters:
-----------

In:  MyRow - Local row.

Out:  NumEntries - Number of nonzero values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrMatrix::MaxNumEntries "int
Epetra_VbrMatrix::MaxNumEntries() const

Returns the maximum of NumMyRowEntries() over all rows. ";

%feature("docstring")  Epetra_VbrMatrix::Map "const Epetra_BlockMap&
Epetra_VbrMatrix::Map() const

Map() method inherited from Epetra_DistObject. ";

%feature("docstring")  Epetra_VbrMatrix::RowMatrixRowMap "const
Epetra_Map& Epetra_VbrMatrix::RowMatrixRowMap() const

Returns the EpetraMap object associated with the rows of this matrix.
";

%feature("docstring")  Epetra_VbrMatrix::RowMatrixColMap "const
Epetra_Map& Epetra_VbrMatrix::RowMatrixColMap() const

Returns the Epetra_Map object associated with columns of this matrix.
";

%feature("docstring")  Epetra_VbrMatrix::RowMatrixImporter "const
Epetra_Import* Epetra_VbrMatrix::RowMatrixImporter() const

Returns the Epetra_Import object that contains the import operations
for distributed operations. ";

/*  Deprecated methods:  These methods still work, but will be removed
in a future version  */

%feature("docstring")  Epetra_VbrMatrix::BlockImportMap "const
Epetra_BlockMap& Epetra_VbrMatrix::BlockImportMap() const

Use BlockColMap() instead. ";

%feature("docstring")  Epetra_VbrMatrix::TransformToLocal "int
Epetra_VbrMatrix::TransformToLocal()

Use FillComplete() instead. ";

%feature("docstring")  Epetra_VbrMatrix::TransformToLocal "int
Epetra_VbrMatrix::TransformToLocal(const Epetra_BlockMap *DomainMap,
const Epetra_BlockMap *RangeMap)

Use FillComplete(const Epetra_BlockMap& DomainMap, const
Epetra_BlockMap& RangeMap) instead. ";


// File: classEpetra__VbrRowMatrix.xml
%feature("docstring") Epetra_VbrRowMatrix "

Epetra_VbrRowMatrix: A class for using an existing Epetra_VbrMatrix
object as an Epetra_RowMatrix object.

The Epetra_VbrRowMatrix class takes an existing Epetra_VbrMatrix
object and allows its use as an Epetra_RowMatrix without allocating
additional storage. Although the Epetra_VbrMatrix itself inherits from
Epetra_RowMatrix, a design flaw in the inheritance structure of Epetra
prohibits the use of an Epetra_VbrMatrix object as an Epetra_RowMatrix
in some important situations. Therefore we recommend the use of this
class to wrap an Epetra_VbrMatrix object.

WARNING:  This class takes a pointer to an existing Epetra_VbrMatrix
object. It is assumed that the user will pass in a pointer to a valid
Epetra_VbrMatrix object, and will retain it throughout the life of the
Epetra_VbrRowMatrix object.

C++ includes: Epetra_VbrRowMatrix.h ";

/*  Constructors/Destructor  */

%feature("docstring")  Epetra_VbrRowMatrix::Epetra_VbrRowMatrix "Epetra_VbrRowMatrix::Epetra_VbrRowMatrix(Epetra_VbrMatrix *Matrix)

Epetra_VbrRowMatrix constuctor. ";

%feature("docstring")  Epetra_VbrRowMatrix::~Epetra_VbrRowMatrix "virtual Epetra_VbrRowMatrix::~Epetra_VbrRowMatrix()

Epetra_VbrRowMatrix Destructor. ";

/*  Post-construction modifications  */

%feature("docstring")  Epetra_VbrRowMatrix::UpdateMatrix "int
Epetra_VbrRowMatrix::UpdateMatrix(Epetra_VbrMatrix *Matrix)

Update the matrix to which this object points. ";

/*  Methods required for implementing Epetra_BasicRowMatrix  */

%feature("docstring")  Epetra_VbrRowMatrix::ExtractMyRowCopy "int
Epetra_VbrRowMatrix::ExtractMyRowCopy(int MyRow, int Length, int
&NumEntries, double *Values, int *Indices) const

Returns a copy of the specified local row in user-provided arrays.

Parameters:
-----------

MyRow:  (In) - Local row to extract.

Length:  (In) - Length of Values and Indices.

NumEntries:  (Out) - Number of nonzero entries extracted.

Values:  (Out) - Extracted values for this row.

Indices:  (Out) - Extracted global column indices for the
corresponding values.

Integer error code, set to 0 if successful, set to -1 if MyRow not
valid, -2 if Length is too short (NumEntries will have required
length). ";

%feature("docstring")  Epetra_VbrRowMatrix::ExtractMyEntryView "int
Epetra_VbrRowMatrix::ExtractMyEntryView(int CurEntry, double *&Value,
int &RowIndex, int &ColIndex)

Returns a reference to the ith entry in the matrix, along with its row
and column index.

Parameters:
-----------

CurEntry:  (In) - Local entry to extract.

Value:  (Out) - Extracted reference to current values.

RowIndex:  (Out) - Row index for current entry.

ColIndex:  (Out) - Column index for current entry.

Integer error code, set to 0 if successful, set to -1 if CurEntry not
valid. ";

%feature("docstring")  Epetra_VbrRowMatrix::ExtractMyEntryView "int
Epetra_VbrRowMatrix::ExtractMyEntryView(int CurEntry, double const
*&Value, int &RowIndex, int &ColIndex) const

Returns a const reference to the ith entry in the matrix, along with
its row and column index.

Parameters:
-----------

CurEntry:  (In) - Local entry to extract.

Value:  (Out) - Extracted reference to current values.

RowIndex:  (Out) - Row index for current entry.

ColIndex:  (Out) - Column index for current entry.

Integer error code, set to 0 if successful, set to -1 if CurEntry not
valid. ";

%feature("docstring")  Epetra_VbrRowMatrix::NumMyRowEntries "int
Epetra_VbrRowMatrix::NumMyRowEntries(int MyRow, int &NumEntries) const

Return the current number of values stored for the specified local
row.

Similar to NumMyEntries() except NumEntries is returned as an argument
and error checking is done on the input value MyRow.

Parameters:
-----------

MyRow:  - (In) Local row.

NumEntries:  - (Out) Number of nonzero values.

Integer error code, set to 0 if successful, set to -1 if MyRow not
valid.

None.

Unchanged. ";

/*  Computational methods  */

%feature("docstring")  Epetra_VbrRowMatrix::RightScale "int
Epetra_VbrRowMatrix::RightScale(const Epetra_Vector &x)

Scales the Epetra_VbrMatrix on the right with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i
denotes the global row number of A and j denotes the global column
number of A.

Parameters:
-----------

In:  x -The Epetra_Vector used for scaling this.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrRowMatrix::LeftScale "int
Epetra_VbrRowMatrix::LeftScale(const Epetra_Vector &x)

Scales the Epetra_VbrMatrix on the left with a Epetra_Vector x.

The this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i
denotes the row number of A and j denotes the column number of A.

Parameters:
-----------

In:  x -A Epetra_Vector to solve for.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrRowMatrix::Multiply "int
Epetra_VbrRowMatrix::Multiply(bool TransA, const Epetra_MultiVector
&X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_VbrRowMatrix multiplied by a
Epetra_MultiVector X in Y.

Parameters:
-----------

In:  TransA -If true, multiply by the transpose of matrix, otherwise
just use matrix.

In:  X - A Epetra_MultiVector of dimension NumVectors to multiply with
matrix.

Out:  Y -A Epetra_MultiVector of dimension NumVectorscontaining
result.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_VbrRowMatrix::Solve "int
Epetra_VbrRowMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal,
const Epetra_MultiVector &X, Epetra_MultiVector &Y) const

Returns the result of a Epetra_VbrRowMatrix solve with a
Epetra_MultiVector X in Y (not implemented).

Parameters:
-----------

In:  Upper -If true, solve Ux = y, otherwise solve Lx = y.

In:  Trans -If true, solve transpose problem.

In:  UnitDiagonal -If true, assume diagonal is unit (whether it's
stored or not).

In:  X - A Epetra_MultiVector of dimension NumVectors to solve for.

Out:  Y -A Epetra_MultiVector of dimension NumVectors containing
result.

Integer error code, set to 0 if successful. ";


// File: classEpetra__Vector.xml
%feature("docstring") Epetra_Vector "

Epetra_Vector: A class for constructing and using dense vectors on a
parallel computer.

The Epetra_Vector class enables the construction and use of real-
valued, double- precision dense vectors in a distributed memory
environment. The distribution of the dense vector is determined in
part by a Epetra_Comm object and a Epetra_Map (or Epetra_LocalMap or
Epetra_BlockMap).

This class is derived from the Epetra_MultiVector class. As such, it
has full access to all of the functionality provided in the
Epetra_MultiVector class.

Distributed Global vs. Replicated Local Distributed Global Vectors -
In most instances, a multi-vector will be partitioned across multiple
memory images associated with multiple processors. In this case, there
is a unique copy of each element and elements are spread across all
processors specified by the Epetra_Comm communicator.

Replicated Local Vectors - Some algorithms use vectors that are too
small to be distributed across all processors. Replicated local
vectors handle these types of situation.

Constructing Epetra_Vectors

There are four Epetra_Vector constructors. The first is a basic
constructor that allocates space and sets all values to zero, the
second is a copy constructor. The third and fourth constructors work
with user data. These constructors have two data access modes: Copy
mode - Allocates memory and makes a copy of the user-provided data. In
this case, the user data is not needed after construction.

View mode - Creates a \"view\" of the user data. In this case, the
user data is required to remain intact for the life of the vector.

WARNING:  View mode is extremely dangerous from a data hiding
perspective. Therefore, we strongly encourage users to develop code
using Copy mode first and only use the View mode in a secondary
optimization phase.  All Epetra_Vector constructors require a map
argument that describes the layout of elements on the parallel
machine. Specifically, map is a Epetra_Map, Epetra_LocalMap or
Epetra_BlockMap object describing the desired memory layout for the
vector.

There are four different Epetra_Vector constructors: Basic - All
values are zero.

Copy - Copy an existing vector.

Copy from or make view of user double array.

Copy or make view of a vector from a Epetra_MultiVector object.

Extracting Data from Epetra_Vectors

Once a Epetra_Vector is constructed, it is possible to extract a copy
of the values or create a view of them.

WARNING:  ExtractView functions are extremely dangerous from a data
hiding perspective. For both ExtractView fuctions, there is a
corresponding ExtractCopy function. We strongly encourage users to
develop code using ExtractCopy functions first and only use the
ExtractView functions in a secondary optimization phase.  There are
two Extract functions: ExtractCopy - Copy values into a user-provided
array.

ExtractView - Set user-provided array to point to Epetra_Vector data.

Vector and Utility Functions

Once a Epetra_Vector is constructed, a variety of mathematical
functions can be applied to the vector. Specifically: Dot Products.

Vector Updates.

p Norms.

Weighted Norms.

Minimum, Maximum and Average Values.

The final useful function is Flops(). Each Epetra_Vector object keep
track of the number of serial floating point operations performed
using the specified object as the this argument to the function. The
Flops() function returns this number as a double precision number.
Using this information, in conjunction with the Epetra_Time class, one
can get accurate parallel performance numbers.

WARNING:  A Epetra_Map, Epetra_LocalMap or Epetra_BlockMap object is
required for all Epetra_Vector constructors.

C++ includes: Epetra_Vector.h ";

/*  Constructors/destructors  */

%feature("docstring")  Epetra_Vector::Epetra_Vector "Epetra_Vector::Epetra_Vector(const Epetra_BlockMap &Map, bool
zeroOut=true)

Basic Epetra_Vector constuctor.

Creates a Epetra_Vector object and fills with zero values.

Parameters:
-----------

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  zeroOut - If true then the allocated memory will be zeroed out
initialy. If false then this memory will not be touched which can be
significantly faster.

WARNING:  Note that, because Epetra_LocalMap derives from Epetra_Map
and Epetra_Map derives from Epetra_BlockMap, this constructor works
for all three types of Epetra map classes.

Pointer to a Epetra_Vector. ";

%feature("docstring")  Epetra_Vector::Epetra_Vector "Epetra_Vector::Epetra_Vector(const Epetra_Vector &Source)

Epetra_Vector copy constructor. ";

%feature("docstring")  Epetra_Vector::Epetra_Vector "Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const
Epetra_BlockMap &Map, double *V)

Set vector values from user array.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  V - Pointer to an array of double precision numbers..

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_Vector::Epetra_Vector "Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const
Epetra_MultiVector &Source, int Index)

Set vector values from a vector in an existing Epetra_MultiVector.

Parameters:
-----------

In:  Epetra_DataAccess - Enumerated type set to Copy or View.

In:  Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

In:  Source - An existing fully constructed Epetra_MultiVector.

In:  Index - Index of vector to access.

Integer error code, set to 0 if successful.  See Detailed Description
section for further discussion. ";

%feature("docstring")  Epetra_Vector::~Epetra_Vector "Epetra_Vector::~Epetra_Vector()

Epetra_Vector destructor. ";

/*  Post-construction modification routines  */

%feature("docstring")  Epetra_Vector::ReplaceGlobalValues "int
Epetra_Vector::ReplaceGlobalValues(int NumEntries, const double
*Values, const int *Indices)

Replace values in a vector with a given indexed list of values,
indices are in global index space.

Replace the Indices[i] entry in the this object with Values[i], for
i=0; i<NumEntries. The indices are in global index space.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in global index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::ReplaceGlobalValues "int
Epetra_Vector::ReplaceGlobalValues(int NumEntries, const double
*Values, const long long *Indices) ";

%feature("docstring")  Epetra_Vector::ReplaceMyValues "int
Epetra_Vector::ReplaceMyValues(int NumEntries, const double *Values,
const int *Indices)

Replace values in a vector with a given indexed list of values,
indices are in local index space.

Replace the Indices[i] entry in the this object with Values[i], for
i=0; i<NumEntries. The indices are in local index space.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in local index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::SumIntoGlobalValues "int
Epetra_Vector::SumIntoGlobalValues(int NumEntries, const double
*Values, const int *Indices)

Sum values into a vector with a given indexed list of values, indices
are in global index space.

Sum Values[i] into the Indices[i] entry in the this object, for i=0;
i<NumEntries. The indices are in global index space.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in global index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::SumIntoMyValues "int
Epetra_Vector::SumIntoMyValues(int NumEntries, const double *Values,
const int *Indices)

Sum values into a vector with a given indexed list of values, indices
are in local index space.

Sum Values[i] into the Indices[i] entry in the this object, for i=0;
i<NumEntries. The indices are in local index space.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in local index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::ReplaceGlobalValues "int
Epetra_Vector::ReplaceGlobalValues(int NumEntries, int BlockOffset,
const double *Values, const int *Indices)

Replace values in a vector with a given indexed list of values at the
specified BlockOffset, indices are in global index space.

Replace the Indices[i] entry in the this object with Values[i], for
i=0; i<NumEntries. The indices are in global index space. This method
is intended for vector that are defined using block maps. In this
situation, an index value is associated with one or more vector
entries, depending on the element size of the given index. The
BlockOffset argument indicates which vector entry to modify as an
offset from the first vector entry associated with the given index.
The offset is used for each entry in the input list.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  BlockOffset - Offset from the first vector entry associated with
each of the given indices.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in global index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::ReplaceMyValues "int
Epetra_Vector::ReplaceMyValues(int NumEntries, int BlockOffset, const
double *Values, const int *Indices)

Replace values in a vector with a given indexed list of values at the
specified BlockOffset, indices are in local index space.

Replace the (Indices[i], BlockOffset) entry in the this object with
Values[i], for i=0; i<NumEntries. The indices are in local index
space. This method is intended for vector that are defined using block
maps. In this situation, an index value is associated with one or more
vector entries, depending on the element size of the given index. The
BlockOffset argument indicates which vector entry to modify as an
offset from the first vector entry associated with the given index.
The offset is used for each entry in the input list.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  BlockOffset - Offset from the first vector entry associated with
each of the given indices.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in local index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::SumIntoGlobalValues "int
Epetra_Vector::SumIntoGlobalValues(int NumEntries, int BlockOffset,
const double *Values, const int *Indices)

Sum values into a vector with a given indexed list of values at the
specified BlockOffset, indices are in global index space.

Sum Values[i] into the Indices[i] entry in the this object, for i=0;
i<NumEntries. The indices are in global index space. This method is
intended for vector that are defined using block maps. In this
situation, an index value is associated with one or more vector
entries, depending on the element size of the given index. The
BlockOffset argument indicates which vector entry to modify as an
offset from the first vector entry associated with the given index.
The offset is used for each entry in the input list.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  BlockOffset - Offset from the first vector entry associated with
each of the given indices.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in global index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

%feature("docstring")  Epetra_Vector::SumIntoMyValues "int
Epetra_Vector::SumIntoMyValues(int NumEntries, int BlockOffset, const
double *Values, const int *Indices)

Sum values into a vector with a given indexed list of values at the
specified BlockOffset, indices are in local index space.

Sum Values[i] into the Indices[i] entry in the this object, for i=0;
i<NumEntries. The indices are in local index space. This method is
intended for vector that are defined using block maps. In this
situation, an index value is associated with one or more vector
entries, depending on the element size of the given index. The
BlockOffset argument indicates which vector entry to modify as an
offset from the first vector entry associated with the given index.
The offset is used for each entry in the input list.

Parameters:
-----------

In:  NumEntries - Number of vector entries to modify.

In:  BlockOffset - Offset from the first vector entry associated with
each of the given indices.

In:  Values - Values which will replace existing values in vector, of
length NumEntries.

In:  Indices - Indices in local index space corresponding to Values.

Integer error code, set to 0 if successful, set to 1 if one or more
indices are not associated with calling processor. ";

/*  Extraction methods  */

%feature("docstring")  Epetra_Vector::ExtractCopy "int
Epetra_Vector::ExtractCopy(double *V) const

Put vector values into user-provided array.

Parameters:
-----------

Out:  V - Pointer to memory space that will contain the vector values.

Integer error code, set to 0 if successful. ";

%feature("docstring")  Epetra_Vector::ExtractView "int
Epetra_Vector::ExtractView(double **V) const

Set user-provided address of V.

Parameters:
-----------

Out:  V - Address of a pointer to that will be set to point to the
values of the vector.

Integer error code, set to 0 if successful. ";

/*  Overloaded operators  */

/*  Expert-only unsupported methods  */

%feature("docstring")  Epetra_Vector::ResetView "int
Epetra_Vector::ResetView(double *Values_in)

Reset the view of an existing vector to point to new user data.

Allows the (very) light-weight replacement of multivector values for
an existing vector that was constructed using an Epetra_DataAccess
mode of View. No checking is performed to see if the values passed in
contain valid data. It is assumed that the user has verified the
integrity of data before calling this method. This method is useful
for situations where a vector is needed for use with an Epetra
operator or matrix and the user is not passing in a multivector, or
the multivector is being passed in with another map that is not
exactly compatible with the operator, but has the correct number of
entries.

This method is used by AztecOO and Ifpack in the matvec and solve
methods to improve performance and reduce repeated calls to
constructors and destructors.

Parameters:
-----------

Values:  Vector data.

Integer error code, set to 0 if successful, -1 if the multivector was
not created as a View.

WARNING:  This method is extremely dangerous and should only be used
by experts. ";


// File: structEpetra__CrsGraphData_1_1IndexData_3_01int_01_4.xml
%feature("docstring") Epetra_CrsGraphData::IndexData< int > " ";

%feature("docstring")  Epetra_CrsGraphData::IndexData< int
>::IndexData " Epetra_CrsGraphData::IndexData< int >::IndexData(int
NumMyBlockRows, bool AllocSorted) ";

%feature("docstring")  Epetra_CrsGraphData::IndexData< int >::Allocate
" void Epetra_CrsGraphData::IndexData< int >::Allocate(int
NumMyBlockRows, bool AllocSorted) ";

%feature("docstring")  Epetra_CrsGraphData::IndexData< int
>::Deallocate " void Epetra_CrsGraphData::IndexData< int
>::Deallocate() ";


// File: structEpetra__CrsGraphData_1_1IndexData_3_01long_01long_01_4.xml
%feature("docstring") Epetra_CrsGraphData::IndexData< long long > " ";

%feature("docstring")  Epetra_CrsGraphData::IndexData< long long
>::IndexData " Epetra_CrsGraphData::IndexData< long long
>::IndexData(int NumMyBlockRows, bool AllocSorted) ";

%feature("docstring")  Epetra_CrsGraphData::IndexData< long long
>::Allocate " void Epetra_CrsGraphData::IndexData< long long
>::Allocate(int NumMyBlockRows, bool AllocSorted) ";

%feature("docstring")  Epetra_CrsGraphData::IndexData< long long
>::Deallocate " void Epetra_CrsGraphData::IndexData< long long
>::Deallocate() ";


// File: structEpetra__MapColoring_1_1ListItem.xml


// File: structEpetra__HashTable_1_1Node.xml


// File: Epetra__BasicDirectory_8cpp.xml


// File: Epetra__BasicDirectory_8h.xml


// File: Epetra__BasicRowMatrix_8cpp.xml


// File: Epetra__BasicRowMatrix_8h.xml


// File: Epetra__BLAS_8cpp.xml


// File: Epetra__BLAS_8h.xml


// File: Epetra__BLAS__wrappers_8h.xml
%feature("docstring")  DASUM_F77 "double PREFIX DASUM_F77(const int
*n, const double x[], const int *incx) ";

%feature("docstring")  DAXPY_F77 "void PREFIX DAXPY_F77(const int *n,
const double *alpha, const double x[], const int *incx, double y[],
const int *incy) ";

%feature("docstring")  DCOPY_F77 "void PREFIX DCOPY_F77(const int *n,
const double *x, const int *incx, double *y, const int *incy) ";

%feature("docstring")  DDOT_F77 "double PREFIX DDOT_F77(const int *n,
const double x[], const int *incx, const double y[], const int *incy)
";

%feature("docstring")  DNRM2_F77 "double PREFIX DNRM2_F77(const int
*n, const double x[], const int *incx) ";

%feature("docstring")  DSCAL_F77 "void PREFIX DSCAL_F77(const int *n,
const double *alpha, double *x, const int *incx) ";

%feature("docstring")  IDAMAX_F77 "int PREFIX IDAMAX_F77(const int
*n, const double *x, const int *incx) ";

%feature("docstring")  SASUM_F77 "float PREFIX SASUM_F77(const int
*n, const float x[], const int *incx) ";

%feature("docstring")  SAXPY_F77 "void PREFIX SAXPY_F77(const int *n,
const float *alpha, const float x[], const int *incx, float y[], const
int *incy) ";

%feature("docstring")  SCOPY_F77 "void PREFIX SCOPY_F77(const int *n,
const float *x, const int *incx, float *y, const int *incy) ";

%feature("docstring")  SDOT_F77 "float PREFIX SDOT_F77(const int *n,
const float x[], const int *incx, const float y[], const int *incy) ";

%feature("docstring")  SNRM2_F77 "float PREFIX SNRM2_F77(const int
*n, const float x[], const int *incx) ";

%feature("docstring")  SSCAL_F77 "void PREFIX SSCAL_F77(const int *n,
const float *alpha, float *x, const int *incx) ";

%feature("docstring")  ISAMAX_F77 "int PREFIX ISAMAX_F77(const int
*n, const float *x, const int *incx) ";

%feature("docstring")  DGEMV_F77 "void PREFIX DGEMV_F77(Epetra_fcd,
const int *m, const int *n, const double *alpha, const double A[],
const int *lda, const double x[], const int *incx, const double *beta,
double y[], const int *incy) ";

%feature("docstring")  DTRMV_F77 "void PREFIX DTRMV_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *n, const double *a, const int *lda,
double *x, const int *incx) ";

%feature("docstring")  DGER_F77 "void PREFIX DGER_F77(const int *m,
const int *n, const double *alpha, const double *x, const int *incx,
const double *y, const int *incy, double *a, const int *lda) ";

%feature("docstring")  SGEMV_F77 "void PREFIX SGEMV_F77(Epetra_fcd,
const int *m, const int *n, const float *alpha, const float A[], const
int *lda, const float x[], const int *incx, const float *beta, float
y[], const int *incy) ";

%feature("docstring")  STRMV_F77 "void PREFIX STRMV_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *n, const float *a, const int *lda,
float *x, const int *incx) ";

%feature("docstring")  SGER_F77 "void PREFIX SGER_F77(const int *m,
const int *n, const float *alpha, const float *x, const int *incx,
const float *y, const int *incy, float *a, const int *lda) ";

%feature("docstring")  DGEMM_F77 "void PREFIX DGEMM_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, const int *k, const double
*alpha, const double *a, const int *lda, const double *b, const int
*ldb, const double *beta, double *c, const int *ldc) ";

%feature("docstring")  DSYMM_F77 "void PREFIX DSYMM_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, const double *alpha, const
double *a, const int *lda, const double *b, const int *ldb, const
double *beta, double *c, const int *ldc) ";

%feature("docstring")  DTRMM_F77 "void PREFIX DTRMM_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *m, const int *n, const
double *alpha, const double *a, const int *lda, double *b, const int
*ldb) ";

%feature("docstring")  DTRSM_F77 "void PREFIX DTRSM_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *m, const int *n, const
double *alpha, const double *a, const int *lda, double *b, const int
*ldb) ";

%feature("docstring")  EPETRA_DCRSMV_F77 "void PREFIX
EPETRA_DCRSMV_F77(const int *, const int *, const int *, const double
*, const int *, const int *, double *, double *) ";

%feature("docstring")  EPETRA_DCRSMM_F77 "void PREFIX
EPETRA_DCRSMM_F77(const int *, const int *, const int *, const double
*, const int *, const int *, double *, int *, double *, int *, int *)
";

%feature("docstring")  EPETRA_DCRSSV_F77 "void PREFIX
EPETRA_DCRSSV_F77(const int *, const int *, const int *, const int *,
const int *, const int *, const double *, const int *, const int *,
double *, double *, const int *) ";

%feature("docstring")  EPETRA_DCRSSM_F77 "void PREFIX
EPETRA_DCRSSM_F77(const int *, const int *, const int *, const int *,
const int *, const int *, const double *, const int *, const int *,
double *, const int *, double *, const int *, const int *, const int
*) ";

%feature("docstring")  DSYRK_F77 "void PREFIX DSYRK_F77(Epetra_fcd
uplo, Epetra_fcd trans, const int *n, const int *k, const double
*alpha, const double *a, const int *lda, const double *beta, double
*c, const int *ldc) ";

%feature("docstring")  SGEMM_F77 "void PREFIX SGEMM_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, const int *k, const float
*alpha, const float *a, const int *lda, const float *b, const int
*ldb, const float *beta, float *c, const int *ldc) ";

%feature("docstring")  SSYMM_F77 "void PREFIX SSYMM_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, const float *alpha, const
float *a, const int *lda, const float *b, const int *ldb, const float
*beta, float *c, const int *ldc) ";

%feature("docstring")  STRMM_F77 "void PREFIX STRMM_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *m, const int *n, const
float *alpha, const float *a, const int *lda, float *b, const int
*ldb) ";

%feature("docstring")  STRSM_F77 "void PREFIX STRSM_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *m, const int *n, const
float *alpha, const float *a, const int *lda, float *b, const int
*ldb) ";

%feature("docstring")  XERBLA_F77 "void PREFIX XERBLA_F77(Epetra_fcd,
int *info) ";

%feature("docstring")  SSYRK_F77 "void PREFIX SSYRK_F77(Epetra_fcd
uplo, Epetra_fcd trans, const int *n, const int *k, const float
*alpha, const float *a, const int *lda, const float *beta, float *c,
const int *ldc) ";


// File: Epetra__BlockMap_8cpp.xml


// File: Epetra__BlockMap_8h.xml


// File: Epetra__BlockMapData_8cpp.xml


// File: Epetra__BlockMapData_8h.xml


// File: Epetra__C__wrappers_8cpp.xml
%feature("docstring")  epetra_serialcomm_create "EPETRA_OBJECT_PTR
MANGLE() epetra_serialcomm_create()

Epetra_Comm ";

%feature("docstring")  epetra_comm_mypid "int MANGLE()
epetra_comm_mypid(EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_comm_numproc "int MANGLE()
epetra_comm_numproc(EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_comm_barrier "void MANGLE()
epetra_comm_barrier(EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_comm_destroy "void MANGLE()
epetra_comm_destroy(EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create1 "EPETRA_OBJECT_PTR MANGLE()
epetra_map_create1(EPETRA_INT numGlobalElements, EPETRA_INT indexBase,
EPETRA_OBJECT_REF comm)

Epetra_Map ";

%feature("docstring")  epetra_map_create2 "EPETRA_OBJECT_PTR MANGLE()
epetra_map_create2(EPETRA_INT numGlobalElements, EPETRA_INT
numMyElements, EPETRA_INT indexBase, EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create3 "EPETRA_OBJECT_PTR MANGLE()
epetra_map_create3(EPETRA_INT numGlobalElements, EPETRA_INT
numLocalElements, int *updateList, EPETRA_INT indexBase,
EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create1_64 "EPETRA_OBJECT_PTR
MANGLE() epetra_map_create1_64(EPETRA_LONG_LONG numGlobalElements,
EPETRA_INT indexBase, EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create2_64 "EPETRA_OBJECT_PTR
MANGLE() epetra_map_create2_64(EPETRA_LONG_LONG numGlobalElements,
EPETRA_INT numMyElements, EPETRA_INT indexBase, EPETRA_OBJECT_REF
comm) ";

%feature("docstring")  epetra_map_create3_64 "EPETRA_OBJECT_PTR
MANGLE() epetra_map_create3_64(EPETRA_LONG_LONG numGlobalElements,
EPETRA_INT numLocalElements, long long *updateList, EPETRA_INT
indexBase, EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_nummyelements "int MANGLE()
epetra_map_nummyelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_numglobalelements "long long
MANGLE() epetra_map_numglobalelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_myglobalelements "int* MANGLE()
epetra_map_myglobalelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_myglobalelements_64 "long long*
MANGLE() epetra_map_myglobalelements_64(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_comm "EPETRA_OBJECT_PTR MANGLE()
epetra_map_comm(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_destroy "void MANGLE()
epetra_map_destroy(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_vector_create1 "EPETRA_OBJECT_PTR
MANGLE() epetra_vector_create1(EPETRA_OBJECT_REF map)

Epetra_Vector ";

%feature("docstring")  epetra_vector_create2 "EPETRA_OBJECT_PTR
MANGLE() epetra_vector_create2(EPETRA_INT CopyValues,
EPETRA_OBJECT_REF map, double *V) ";

%feature("docstring")  epetra_vector_putscalar "int MANGLE()
epetra_vector_putscalar(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalar) ";

%feature("docstring")  epetra_vector_norm1 "int MANGLE()
epetra_vector_norm1(EPETRA_OBJECT_REF x, double *scalar) ";

%feature("docstring")  epetra_vector_norm2 "int MANGLE()
epetra_vector_norm2(EPETRA_OBJECT_REF x, double *scalar) ";

%feature("docstring")  epetra_vector_random "int MANGLE()
epetra_vector_random(EPETRA_OBJECT_REF x) ";

%feature("docstring")  epetra_vector_update "int MANGLE()
epetra_vector_update(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalara,
EPETRA_OBJECT_REF a, EPETRA_DOUBLE scalarb, EPETRA_OBJECT_REF b,
EPETRA_DOUBLE scalarx) ";

%feature("docstring")  epetra_vector_print "void MANGLE()
epetra_vector_print(EPETRA_OBJECT_REF x) ";

%feature("docstring")  epetra_vector_destroy "void MANGLE()
epetra_vector_destroy(EPETRA_OBJECT_REF x) ";


// File: Epetra__C__wrappers_8h.xml
%feature("docstring")  epetra_serialcomm_create "EPETRA_OBJECT_PTR
MANGLE() epetra_serialcomm_create()

Epetra_Comm ";

%feature("docstring")  epetra_comm_mypid "int MANGLE()
epetra_comm_mypid(EPETRA_OBJECT_REF communicator) ";

%feature("docstring")  epetra_comm_numproc "int MANGLE()
epetra_comm_numproc(EPETRA_OBJECT_REF communicator) ";

%feature("docstring")  epetra_comm_barrier "void MANGLE()
epetra_comm_barrier(EPETRA_OBJECT_REF communicator) ";

%feature("docstring")  epetra_comm_destroy "void MANGLE()
epetra_comm_destroy(EPETRA_OBJECT_REF communicator) ";

%feature("docstring")  epetra_map_create1 "EPETRA_OBJECT_PTR MANGLE()
epetra_map_create1(EPETRA_INT numGlobalEquations, EPETRA_INT
indexBase, EPETRA_OBJECT_REF comm)

Epetra_Map ";

%feature("docstring")  epetra_map_create2 "EPETRA_OBJECT_PTR MANGLE()
epetra_map_create2(EPETRA_INT numGlobalEquations, EPETRA_INT
numMyElements, EPETRA_INT indexBase, EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create3 "EPETRA_OBJECT_PTR MANGLE()
epetra_map_create3(EPETRA_INT numGlobalEquations, EPETRA_INT
numlocalEquations, int *updateList, EPETRA_INT indexBase,
EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create1_64 "EPETRA_OBJECT_PTR
MANGLE() epetra_map_create1_64(EPETRA_LONG_LONG numGlobalEquations,
EPETRA_INT indexBase, EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_create2_64 "EPETRA_OBJECT_PTR
MANGLE() epetra_map_create2_64(EPETRA_LONG_LONG numGlobalEquations,
EPETRA_INT numMyElements, EPETRA_INT indexBase, EPETRA_OBJECT_REF
comm) ";

%feature("docstring")  epetra_map_create3_64 "EPETRA_OBJECT_PTR
MANGLE() epetra_map_create3_64(EPETRA_LONG_LONG numGlobalEquations,
EPETRA_INT numlocalEquations, long long *updateList, EPETRA_INT
indexBase, EPETRA_OBJECT_REF comm) ";

%feature("docstring")  epetra_map_nummyelements "int MANGLE()
epetra_map_nummyelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_numglobalelements "long long
MANGLE() epetra_map_numglobalelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_myglobalelements "int* MANGLE()
epetra_map_myglobalelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_myglobalelements_64 "long long*
MANGLE() epetra_map_myglobalelements_64(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_comm "EPETRA_OBJECT_PTR MANGLE()
epetra_map_comm(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_destroy "void MANGLE()
epetra_map_destroy(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_vector_create1 "EPETRA_OBJECT_PTR
MANGLE() epetra_vector_create1(EPETRA_OBJECT_REF map)

Epetra_Vector ";

%feature("docstring")  epetra_vector_create2 "EPETRA_OBJECT_PTR
MANGLE() epetra_vector_create2(EPETRA_INT Copy, EPETRA_OBJECT_REF map,
double *V) ";

%feature("docstring")  epetra_vector_putscalar "int MANGLE()
epetra_vector_putscalar(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalar) ";

%feature("docstring")  epetra_vector_update "int MANGLE()
epetra_vector_update(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalara,
EPETRA_OBJECT_REF a, EPETRA_DOUBLE scalarb, EPETRA_OBJECT_REF b,
EPETRA_DOUBLE scalarx) ";

%feature("docstring")  epetra_vector_norm1 "int MANGLE()
epetra_vector_norm1(EPETRA_OBJECT_REF x, double *result) ";

%feature("docstring")  epetra_vector_norm2 "int MANGLE()
epetra_vector_norm2(EPETRA_OBJECT_REF x, double *result) ";

%feature("docstring")  epetra_vector_random "int MANGLE()
epetra_vector_random(EPETRA_OBJECT_REF x) ";

%feature("docstring")  epetra_vector_print "void MANGLE()
epetra_vector_print(EPETRA_OBJECT_REF x) ";

%feature("docstring")  epetra_vector_destroy "void MANGLE()
epetra_vector_destroy(EPETRA_OBJECT_REF x) ";


// File: Epetra__CombineMode_8h.xml


// File: Epetra__Comm_8h.xml


// File: Epetra__CompObject_8cpp.xml


// File: Epetra__CompObject_8h.xml


// File: Epetra__ConfigDefs_8h.xml


// File: Epetra__CrsGraph_8cpp.xml
%feature("docstring")  epetra_shellsort "void epetra_shellsort(int
*list, int length) ";

%feature("docstring")  epetra_crsgraph_compress_out_duplicates "void
epetra_crsgraph_compress_out_duplicates(int len, int *list, int
&newlen)

*!*!*!

*!*!*! ";


// File: Epetra__CrsGraph_8h.xml


// File: Epetra__CrsGraphData_8cpp.xml


// File: Epetra__CrsGraphData_8h.xml


// File: Epetra__CrsMatrix_8cpp.xml


// File: Epetra__CrsMatrix_8h.xml


// File: Epetra__CrsSingletonFilter_8cpp.xml


// File: Epetra__CrsSingletonFilter_8h.xml


// File: Epetra__Data_8cpp.xml


// File: Epetra__Data_8h.xml


// File: Epetra__DataAccess_8h.xml


// File: Epetra__Directory_8h.xml


// File: Epetra__DistObject_8cpp.xml


// File: Epetra__DistObject_8h.xml


// File: Epetra__Distributor_8h.xml


// File: Epetra__Export_8cpp.xml


// File: Epetra__Export_8h.xml


// File: Epetra__FECrsGraph_8cpp.xml


// File: Epetra__FECrsGraph_8h.xml


// File: Epetra__FECrsMatrix_8cpp.xml


// File: Epetra__FECrsMatrix_8h.xml


// File: Epetra__FEVbrMatrix_8cpp.xml


// File: Epetra__FEVbrMatrix_8h.xml


// File: Epetra__FEVector_8cpp.xml


// File: Epetra__FEVector_8h.xml


// File: Epetra__Flops_8cpp.xml


// File: Epetra__Flops_8h.xml


// File: Epetra__Fortran__wrappers_8cpp.xml


// File: Epetra__Fortran__wrappers_8h.xml


// File: Epetra__HashTable_8h.xml


// File: Epetra__Import_8cpp.xml


// File: Epetra__Import_8h.xml


// File: Epetra__IntSerialDenseMatrix_8cpp.xml


// File: Epetra__IntSerialDenseMatrix_8h.xml


// File: Epetra__IntSerialDenseVector_8cpp.xml


// File: Epetra__IntSerialDenseVector_8h.xml


// File: Epetra__IntVector_8cpp.xml


// File: Epetra__IntVector_8h.xml


// File: Epetra__InvOperator_8cpp.xml


// File: Epetra__InvOperator_8h.xml


// File: Epetra__JadMatrix_8cpp.xml


// File: Epetra__JadMatrix_8h.xml


// File: Epetra__LAPACK_8cpp.xml


// File: Epetra__LAPACK_8h.xml


// File: Epetra__LAPACK__wrappers_8h.xml
%feature("docstring")  DGECON_F77 "void PREFIX DGECON_F77(Epetra_fcd
norm, const int *n, const double *a, const int *lda, const double
*anorm, double *rcond, double *work, int *iwork, int *info) ";

%feature("docstring")  DGEEQU_F77 "void PREFIX DGEEQU_F77(const int
*m, const int *n, const double *a, const int *lda, double *r, double
*c, double *rowcnd, double *colcnd, double *amax, int *info) ";

%feature("docstring")  DGEEV_F77 "void PREFIX DGEEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, double *a, const int *lda, double *wr,
double *wi, double *vl, const int *ldvl, double *vr, const int *ldvr,
double *work, const int *lwork, int *info) ";

%feature("docstring")  DGEEVX_F77 "void PREFIX DGEEVX_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *n, double *a, const int
*lda, double *wr, double *wi, double *vl, const int *ldvl, double *vr,
const int *ldvr, int *ilo, int *ihi, double *scale, double *abnrm,
double *rconde, double *rcondv, double *work, const int *lwork, int
*iwork, int *info) ";

%feature("docstring")  DGEHRD_F77 "void PREFIX DGEHRD_F77(const int
*n, const int *ilo, const int *ihi, double *A, const int *lda, double
*tau, double *work, const int *lwork, int *info) ";

%feature("docstring")  DGELS_F77 "void PREFIX DGELS_F77(Epetra_fcd
ch, const int *m, const int *n, const int *nrhs, double *a, const int
*lda, double *b, const int *ldb, double *work, const int *lwork, int
*info) ";

%feature("docstring")  DGELSS_F77 "void PREFIX DGELSS_F77(const int
*m, const int *n, const int *nrhs, double *a, const int *lda, double
*b, const int *ldb, double *s, const double *rcond, int *rank, double
*work, const int *lwork, int *info) ";

%feature("docstring")  DGEQPF_F77 "void PREFIX DGEQPF_F77(const int
*m, const int *n, double *a, const int *lda, int *jpvt, double *tau,
double *work, int *info) ";

%feature("docstring")  DGERFS_F77 "void PREFIX DGERFS_F77(Epetra_fcd,
const int *n, const int *nrhs, const double *a, const int *lda, const
double *af, const int *ldaf, const int *ipiv, const double *b, const
int *ldb, double *x, const int *ldx, double *ferr, double *berr,
double *work, int *iwork, int *info) ";

%feature("docstring")  DGESDD_F77 "void PREFIX DGESDD_F77(Epetra_fcd,
const int *m, const int *n, double *a, const int *lda, double *s,
double *u, const int *ldu, double *vt, const int *ldvt, double *work,
const int *lwork, int *iwork, int *info) ";

%feature("docstring")  DGESVD_F77 "void PREFIX DGESVD_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, double *a, const int *lda,
double *s, double *u, const int *ldu, double *vt, const int *ldvt,
double *work, const int *lwork, int *info) ";

%feature("docstring")  DGESV_F77 "void PREFIX DGESV_F77(const int *n,
const int *nrhs, double *a, const int *lda, int *ipiv, double *x,
const int *ldx, int *info) ";

%feature("docstring")  DGESVX_F77 "void PREFIX DGESVX_F77(Epetra_fcd,
Epetra_fcd, const int *n, const int *nrhs, double *a, const int *lda,
double *af, const int *ldaf, int *ipiv, Epetra_fcd, double *r, double
*c, double *b, const int *ldb, double *x, const int *ldx, double
*rcond, double *ferr, double *berr, double *work, int *iwork, int
*info) ";

%feature("docstring")  DGETRF_F77 "void PREFIX DGETRF_F77(const int
*m, const int *n, double *a, const int *lda, int *ipiv, int *info) ";

%feature("docstring")  DGEQRF_F77 "void PREFIX DGEQRF_F77(const int
*m, const int *n, double *a, const int *lda, double *tau, double
*work, const int *lwork, int *info) ";

%feature("docstring")  DGETRI_F77 "void PREFIX DGETRI_F77(const int
*n, double *a, const int *lda, int *ipiv, double *work, const int
*lwork, int *info) ";

%feature("docstring")  DGETRS_F77 "void PREFIX DGETRS_F77(Epetra_fcd,
const int *n, const int *nrhs, const double *a, const int *lda, const
int *ipiv, double *x, const int *ldx, int *info) ";

%feature("docstring")  DGGEV_F77 "void PREFIX DGGEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, double *a, const int *lda, double *b, const
int *ldb, double *alphar, double *alphai, double *beta, double *vl,
const int *ldvl, double *vr, const int *ldvr, double *work, const int
*lwork, int *info) ";

%feature("docstring")  DGGLSE_F77 "void PREFIX DGGLSE_F77(const int
*m, const int *n, const int *p, double *a, const int *lda, double *b,
const int *ldb, double *c, double *d, double *x, double *work, const
int *lwork, int *info) ";

%feature("docstring")  DGGSVD_F77 "void PREFIX DGGSVD_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *m, const int *n, const int *p, int
*k, int *l, double *a, const int *lda, double *b, const int *ldb,
double *alpha, double *beta, double *u, const int *ldu, double *v,
const int *ldv, double *q, const int *ldq, double *work, int *iwork,
int *info) ";

%feature("docstring")  DHSEQR_F77 "void PREFIX DHSEQR_F77(Epetra_fcd
job, Epetra_fcd, const int *n, const int *ilo, const int *ihi, double
*h, const int *ldh, double *wr, double *wi, double *z, const int *ldz,
double *work, const int *lwork, int *info) ";

%feature("docstring")  DLAMCH_F77 "double PREFIX
DLAMCH_F77(Epetra_fcd) ";

%feature("docstring")  DLARFT_F77 "void PREFIX DLARFT_F77(Epetra_fcd
direct, Epetra_fcd storev, const int *n, const int *k, double *v,
const int *ldv, double *tau, double *t, const int *ldt) ";

%feature("docstring")  DORGQR_F77 "void PREFIX DORGQR_F77(const int
*m, const int *n, const int *k, double *a, const int *lda, const
double *tau, double *work, const int *lwork, int *info) ";

%feature("docstring")  DORGHR_F77 "void PREFIX DORGHR_F77(const int
*n, const int *ilo, const int *ihi, double *a, const int *lda, const
double *tau, double *work, const int *lwork, int *info) ";

%feature("docstring")  DORMHR_F77 "void PREFIX DORMHR_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, const int *ilo, const int
*ihi, const double *a, const int *lda, const double *tau, double *c,
const int *ldc, double *work, const int *lwork, int *info) ";

%feature("docstring")  DPOCON_F77 "void PREFIX DPOCON_F77(Epetra_fcd,
const int *n, const double *a, const int *lda, const double *anorm,
double *rcond, double *work, int *iwork, int *info) ";

%feature("docstring")  DPOEQU_F77 "void PREFIX DPOEQU_F77(const int
*n, const double *a, const int *lda, double *s, double *scond, double
*amax, int *info) ";

%feature("docstring")  DPORFS_F77 "void PREFIX DPORFS_F77(Epetra_fcd,
const int *n, const int *nrhs, const double *a, const int *lda, const
double *af, const int *ldaf, const double *b, const int *ldb, double
*x, const int *ldx, double *ferr, double *berr, double *work, int
*iwork, int *info) ";

%feature("docstring")  DPOSV_F77 "void PREFIX DPOSV_F77(Epetra_fcd,
const int *n, const int *nrhs, const double *a, const int *lda, double
*x, const int *ldx, int *info) ";

%feature("docstring")  DPOSVX_F77 "void PREFIX DPOSVX_F77(Epetra_fcd,
Epetra_fcd, const int *n, const int *nrhs, double *a, const int *lda,
double *af, const int *ldaf, Epetra_fcd, double *s, double *b, const
int *ldb, double *x, const int *ldx, double *rcond, double *ferr,
double *berr, double *work, int *iwork, int *info) ";

%feature("docstring")  DPOTRF_F77 "void PREFIX DPOTRF_F77(Epetra_fcd,
const int *n, double *a, const int *lda, int *info) ";

%feature("docstring")  DPOTRI_F77 "void PREFIX DPOTRI_F77(Epetra_fcd,
const int *n, double *a, const int *lda, int *info) ";

%feature("docstring")  DPOTRS_F77 "void PREFIX DPOTRS_F77(Epetra_fcd,
const int *n, const int *nrhs, const double *a, const int *lda, double
*x, const int *ldx, int *info) ";

%feature("docstring")  DSPEV_F77 "void PREFIX DSPEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, double *ap, double *w, double *z, const int
*ldz, double *work, int *info) ";

%feature("docstring")  DSPGV_F77 "void PREFIX DSPGV_F77(const int
*itype, Epetra_fcd, Epetra_fcd, const int *n, double *ap, double *bp,
double *w, double *z, const int *ldz, double *work, int *info) ";

%feature("docstring")  DSTEV_F77 "void PREFIX DSTEV_F77(Epetra_fcd
jobz, const int *n, double *d, double *e, double *z, const int *ldz,
double *work, int *info) ";

%feature("docstring")  DSYEVD_F77 "void PREFIX DSYEVD_F77(Epetra_fcd,
Epetra_fcd, const int *n, double *a, const int *lda, double *w, double
*work, const int *lwork, int *iwork, const int *liwork, int *info) ";

%feature("docstring")  DSYEV_F77 "void PREFIX DSYEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, double *a, const int *lda, double *w, double
*work, const int *lwork, int *info) ";

%feature("docstring")  DSYEVR_F77 "void PREFIX DSYEVR_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *n, double *a, const int *lda, const
double *vl, const double *vu, const int *il, const int *iu, const
double *abstol, int *m, double *w, double *z, const int *ldz, int
*isuppz, double *work, const int *lwork, int *iwork, const int
*liwork, int *info) ";

%feature("docstring")  DSYEVX_F77 "void PREFIX DSYEVX_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *n, double *a, const int *lda, const
double *vl, const double *vu, const int *il, const int *iu, const
double *abstol, int *m, double *w, double *z, const int *ldz, double
*work, const int *lwork, int *iwork, int *ifail, int *info) ";

%feature("docstring")  DSYGV_F77 "void PREFIX DSYGV_F77(const int
*itype, Epetra_fcd, Epetra_fcd, const int *n, double *a, const int
*lda, double *b, const int *ldb, double *w, double *work, const int
*lwork, int *info) ";

%feature("docstring")  DSYGVX_F77 "void PREFIX DSYGVX_F77(const int
*itype, Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *n, double *a,
const int *lda, double *b, const int *ldb, const double *vl, const
double *vu, const int *il, const int *iu, const double *abstol, int
*m, double *w, double *z, const int *ldz, double *work, const int
*lwork, int *iwork, int *ifail, int *info) ";

%feature("docstring")  DTREVC_F77 "void PREFIX DTREVC_F77(Epetra_fcd,
Epetra_fcd, int *select, const int *n, const double *t, const int
*ldt, double *vl, const int *ldvl, double *vr, const int *ldvr, const
int *mm, int *m, double *work, int *info) ";

%feature("docstring")  DTREXC_F77 "void PREFIX DTREXC_F77(Epetra_fcd,
const int *n, double *t, const int *ldt, double *q, const int *ldq,
int *ifst, int *ilst, double *work, int *info) ";

%feature("docstring")  DTRTRS_F77 "void PREFIX DTRTRS_F77(Epetra_fcd
uplo, Epetra_fcd trans, Epetra_fcd diag, const int *n, const int
*nrhs, const double *a, const int *lda, double *b, const int *ldb, int
*info) ";

%feature("docstring")  SGECON_F77 "void PREFIX SGECON_F77(Epetra_fcd
norm, const int *n, const float *a, const int *lda, const float
*anorm, float *rcond, float *work, int *iwork, int *info) ";

%feature("docstring")  SGEEQU_F77 "void PREFIX SGEEQU_F77(const int
*m, const int *n, const float *a, const int *lda, float *r, float *c,
float *rowcnd, float *colcnd, float *amax, int *info) ";

%feature("docstring")  SGEEV_F77 "void PREFIX SGEEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, float *a, const int *lda, float *wr, float
*wi, float *vl, const int *ldvl, float *vr, const int *ldvr, float
*work, const int *lwork, int *info) ";

%feature("docstring")  SGEEVX_F77 "void PREFIX SGEEVX_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *n, float *a, const int
*lda, float *wr, float *wi, float *vl, const int *ldvl, float *vr,
const int *ldvr, int *ilo, int *ihi, float *scale, float *abnrm, float
*rconde, float *rcondv, float *work, const int *lwork, int *iwork, int
*info) ";

%feature("docstring")  SGEHRD_F77 "void PREFIX SGEHRD_F77(const int
*n, const int *ilo, const int *ihi, float *A, const int *lda, float
*tau, float *work, const int *lwork, int *info) ";

%feature("docstring")  SGELS_F77 "void PREFIX SGELS_F77(Epetra_fcd
ch, const int *m, const int *n, const int *nrhs, float *a, const int
*lda, float *b, const int *ldb, float *work, const int *lwork, int
*info) ";

%feature("docstring")  SGELSS_F77 "void PREFIX SGELSS_F77(const int
*m, const int *n, const int *nrhs, float *a, const int *lda, float *b,
const int *ldb, float *s, const float *rcond, int *rank, float *work,
const int *lwork, int *info) ";

%feature("docstring")  SGEQPF_F77 "void PREFIX SGEQPF_F77(const int
*m, const int *n, float *a, const int *lda, int *jpvt, float *tau,
float *work, int *info) ";

%feature("docstring")  SGERFS_F77 "void PREFIX SGERFS_F77(Epetra_fcd,
const int *n, const int *nrhs, const float *a, const int *lda, const
float *af, const int *ldaf, const int *ipiv, const float *b, const int
*ldb, float *x, const int *ldx, float *ferr, float *berr, float *work,
int *iwork, int *info) ";

%feature("docstring")  SGESDD_F77 "void PREFIX SGESDD_F77(Epetra_fcd,
const int *m, const int *n, float *a, const int *lda, float *s, float
*u, const int *ldu, float *vt, const int *ldvt, float *work, const int
*lwork, int *iwork, int *info) ";

%feature("docstring")  SGESVD_F77 "void PREFIX SGESVD_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, float *a, const int *lda,
float *s, float *u, const int *ldu, float *vt, const int *ldvt, float
*work, const int *lwork, int *info) ";

%feature("docstring")  SGESV_F77 "void PREFIX SGESV_F77(const int *n,
const int *nrhs, float *a, const int *lda, int *ipiv, float *x, const
int *ldx, int *info) ";

%feature("docstring")  SGESVX_F77 "void PREFIX SGESVX_F77(Epetra_fcd,
Epetra_fcd, const int *n, const int *nrhs, float *a, const int *lda,
float *af, const int *ldaf, int *ipiv, Epetra_fcd, float *r, float *c,
float *b, const int *ldb, float *x, const int *ldx, float *rcond,
float *ferr, float *berr, float *work, int *iwork, int *info) ";

%feature("docstring")  SGETRF_F77 "void PREFIX SGETRF_F77(const int
*m, const int *n, float *a, const int *lda, int *ipiv, int *info) ";

%feature("docstring")  SGEQRF_F77 "void PREFIX SGEQRF_F77(const int
*m, const int *n, float *a, const int *lda, float *tau, float *work,
const int *lwork, int *info) ";

%feature("docstring")  SGETRI_F77 "void PREFIX SGETRI_F77(const int
*n, float *a, const int *lda, int *ipiv, float *work, const int
*lwork, int *info) ";

%feature("docstring")  SGETRS_F77 "void PREFIX SGETRS_F77(Epetra_fcd,
const int *n, const int *nrhs, const float *a, const int *lda, const
int *ipiv, float *x, const int *ldx, int *info) ";

%feature("docstring")  SGGEV_F77 "void PREFIX SGGEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, float *a, const int *lda, float *b, const
int *ldb, float *alphar, float *alphai, float *beta, float *vl, const
int *ldvl, float *vr, const int *ldvr, float *work, const int *lwork,
int *info) ";

%feature("docstring")  SGGLSE_F77 "void PREFIX SGGLSE_F77(const int
*m, const int *n, const int *p, float *a, const int *lda, float *b,
const int *ldb, float *c, float *d, float *x, float *work, const int
*lwork, int *info) ";

%feature("docstring")  SGGSVD_F77 "void PREFIX SGGSVD_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *m, const int *n, const int *p, int
*k, int *l, float *a, const int *lda, float *b, const int *ldb, float
*alpha, float *beta, float *u, const int *ldu, float *v, const int
*ldv, float *q, const int *ldq, float *work, int *iwork, int *info) ";

%feature("docstring")  SHSEQR_F77 "void PREFIX SHSEQR_F77(Epetra_fcd
job, Epetra_fcd, const int *n, const int *ilo, const int *ihi, float
*h, const int *ldh, float *wr, float *wi, float *z, const int *ldz,
float *work, const int *lwork, int *info) ";

%feature("docstring")  SLAMCH_F77 "float PREFIX
SLAMCH_F77(Epetra_fcd) ";

%feature("docstring")  SLARFT_F77 "void PREFIX SLARFT_F77(Epetra_fcd
direct, Epetra_fcd storev, const int *n, const int *k, float *v, const
int *ldv, float *tau, float *t, const int *ldt) ";

%feature("docstring")  SORGQR_F77 "void PREFIX SORGQR_F77(const int
*m, const int *n, const int *k, float *a, const int *lda, const float
*tau, float *work, const int *lwork, int *info) ";

%feature("docstring")  SORGHR_F77 "void PREFIX SORGHR_F77(const int
*n, const int *ilo, const int *ihi, float *a, const int *lda, const
float *tau, float *work, const int *lwork, int *info) ";

%feature("docstring")  SORMHR_F77 "void PREFIX SORMHR_F77(Epetra_fcd,
Epetra_fcd, const int *m, const int *n, const int *ilo, const int
*ihi, const float *a, const int *lda, const float *tau, float *c,
const int *ldc, float *work, const int *lwork, int *info) ";

%feature("docstring")  SPOCON_F77 "void PREFIX SPOCON_F77(Epetra_fcd,
const int *n, const float *a, const int *lda, const float *anorm,
float *rcond, float *work, int *iwork, int *info) ";

%feature("docstring")  SPOEQU_F77 "void PREFIX SPOEQU_F77(const int
*n, const float *a, const int *lda, float *s, float *scond, float
*amax, int *info) ";

%feature("docstring")  SPORFS_F77 "void PREFIX SPORFS_F77(Epetra_fcd,
const int *n, const int *nrhs, const float *a, const int *lda, const
float *af, const int *ldaf, const float *b, const int *ldb, float *x,
const int *ldx, float *ferr, float *berr, float *work, int *iwork, int
*info) ";

%feature("docstring")  SPOSV_F77 "void PREFIX SPOSV_F77(Epetra_fcd,
const int *n, const int *nrhs, const float *a, const int *lda, float
*x, const int *ldx, int *info) ";

%feature("docstring")  SPOSVX_F77 "void PREFIX SPOSVX_F77(Epetra_fcd,
Epetra_fcd, const int *n, const int *nrhs, float *a, const int *lda,
float *af, const int *ldaf, Epetra_fcd, float *s, float *b, const int
*ldb, float *x, const int *ldx, float *rcond, float *ferr, float
*berr, float *work, int *iwork, int *info) ";

%feature("docstring")  SPOTRF_F77 "void PREFIX SPOTRF_F77(Epetra_fcd,
const int *n, float *a, const int *lda, int *info) ";

%feature("docstring")  SPOTRI_F77 "void PREFIX SPOTRI_F77(Epetra_fcd,
const int *n, float *a, const int *lda, int *info) ";

%feature("docstring")  SPOTRS_F77 "void PREFIX SPOTRS_F77(Epetra_fcd,
const int *n, const int *nrhs, const float *a, const int *lda, float
*x, const int *ldx, int *info) ";

%feature("docstring")  SSPEV_F77 "void PREFIX SSPEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, float *ap, float *w, float *z, const int
*ldz, float *work, int *info) ";

%feature("docstring")  SSPGV_F77 "void PREFIX SSPGV_F77(const int
*itype, Epetra_fcd, Epetra_fcd, const int *n, float *ap, float *bp,
float *w, float *z, const int *ldz, float *work, int *info) ";

%feature("docstring")  SSTEV_F77 "void PREFIX SSTEV_F77(Epetra_fcd
jobz, const int *n, float *d, float *e, float *z, const int *ldz,
float *work, int *info) ";

%feature("docstring")  SSYEVD_F77 "void PREFIX SSYEVD_F77(Epetra_fcd,
Epetra_fcd, const int *n, float *a, const int *lda, float *w, float
*work, const int *lwork, int *iwork, const int *liwork, int *info) ";

%feature("docstring")  SSYEV_F77 "void PREFIX SSYEV_F77(Epetra_fcd,
Epetra_fcd, const int *n, float *a, const int *lda, float *w, float
*work, const int *lwork, int *info) ";

%feature("docstring")  SSYEVR_F77 "void PREFIX SSYEVR_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *n, float *a, const int *lda, const
float *vl, const float *vu, const int *il, const int *iu, const float
*abstol, int *m, float *w, float *z, const int *ldz, int *isuppz,
float *work, const int *lwork, int *iwork, const int *liwork, int
*info) ";

%feature("docstring")  SSYEVX_F77 "void PREFIX SSYEVX_F77(Epetra_fcd,
Epetra_fcd, Epetra_fcd, const int *n, float *a, const int *lda, const
float *vl, const float *vu, const int *il, const int *iu, const float
*abstol, int *m, float *w, float *z, const int *ldz, float *work,
const int *lwork, int *iwork, int *ifail, int *info) ";

%feature("docstring")  SSYGV_F77 "void PREFIX SSYGV_F77(const int
*itype, Epetra_fcd, Epetra_fcd, const int *n, float *a, const int
*lda, float *b, const int *ldb, float *w, float *work, const int
*lwork, int *info) ";

%feature("docstring")  SSYGVX_F77 "void PREFIX SSYGVX_F77(const int
*itype, Epetra_fcd, Epetra_fcd, Epetra_fcd, const int *n, float *a,
const int *lda, float *b, const int *ldb, const float *vl, const float
*vu, const int *il, const int *iu, const float *abstol, int *m, float
*w, float *z, const int *ldz, float *work, const int *lwork, int
*iwork, int *ifail, int *info) ";

%feature("docstring")  STREVC_F77 "void PREFIX STREVC_F77(Epetra_fcd,
Epetra_fcd, int *select, const int *n, const float *t, const int *ldt,
float *vl, const int *ldvl, float *vr, const int *ldvr, const int *mm,
int *m, float *work, int *info) ";

%feature("docstring")  STREXC_F77 "void PREFIX STREXC_F77(Epetra_fcd,
const int *n, float *t, const int *ldt, float *q, const int *ldq, int
*ifst, int *ilst, float *work, int *info) ";

%feature("docstring")  STRTRS_F77 "void PREFIX STRTRS_F77(Epetra_fcd
uplo, Epetra_fcd trans, Epetra_fcd diag, const int *n, const int
*nrhs, const float *a, const int *lda, float *b, const int *ldb, int
*info) ";


// File: Epetra__LinearProblem_8cpp.xml


// File: Epetra__LinearProblem_8h.xml


// File: Epetra__LinearProblemRedistor_8cpp.xml


// File: Epetra__LinearProblemRedistor_8h.xml


// File: Epetra__LocalMap_8cpp.xml


// File: Epetra__LocalMap_8h.xml


// File: Epetra__LongLongSerialDenseMatrix_8cpp.xml


// File: Epetra__LongLongSerialDenseMatrix_8h.xml


// File: Epetra__LongLongSerialDenseVector_8cpp.xml


// File: Epetra__LongLongSerialDenseVector_8h.xml


// File: Epetra__LongLongVector_8cpp.xml


// File: Epetra__LongLongVector_8h.xml


// File: Epetra__Map_8cpp.xml


// File: Epetra__Map_8h.xml


// File: Epetra__MapColoring_8cpp.xml


// File: Epetra__MapColoring_8h.xml


// File: Epetra__MpiComm_8cpp.xml


// File: Epetra__MpiComm_8h.xml


// File: Epetra__MpiCommData_8cpp.xml


// File: Epetra__MpiCommData_8h.xml


// File: Epetra__MpiDistributor_8cpp.xml


// File: Epetra__MpiDistributor_8h.xml


// File: Epetra__MpiSmpComm_8cpp.xml


// File: Epetra__MpiSmpComm_8h.xml


// File: Epetra__MpiSmpCommData_8cpp.xml


// File: Epetra__MpiSmpCommData_8h.xml


// File: Epetra__MultiVector_8cpp.xml


// File: Epetra__MultiVector_8h.xml


// File: Epetra__Object_8cpp.xml


// File: Epetra__Object_8h.xml


// File: Epetra__OffsetIndex_8cpp.xml


// File: Epetra__OffsetIndex_8h.xml


// File: Epetra__Operator_8h.xml


// File: Epetra__OskiError_8cpp.xml


// File: Epetra__OskiError_8h.xml


// File: Epetra__OskiMatrix_8cpp.xml


// File: Epetra__OskiMatrix_8h.xml


// File: Epetra__OskiMultiVector_8cpp.xml


// File: Epetra__OskiMultiVector_8h.xml


// File: Epetra__OskiPermutation_8cpp.xml


// File: Epetra__OskiPermutation_8h.xml


// File: Epetra__OskiUtils_8cpp.xml


// File: Epetra__OskiUtils_8h.xml


// File: Epetra__OskiVector_8cpp.xml


// File: Epetra__OskiVector_8h.xml


// File: Epetra__RowMatrix_8h.xml


// File: Epetra__RowMatrixTransposer_8cpp.xml


// File: Epetra__RowMatrixTransposer_8h.xml


// File: Epetra__SerialComm_8cpp.xml


// File: Epetra__SerialComm_8h.xml


// File: Epetra__SerialCommData_8cpp.xml


// File: Epetra__SerialCommData_8h.xml


// File: Epetra__SerialDenseMatrix_8cpp.xml


// File: Epetra__SerialDenseMatrix_8h.xml


// File: Epetra__SerialDenseOperator_8h.xml


// File: Epetra__SerialDenseSolver_8cpp.xml


// File: Epetra__SerialDenseSolver_8h.xml


// File: Epetra__SerialDenseSVD_8cpp.xml


// File: Epetra__SerialDenseSVD_8h.xml


// File: Epetra__SerialDenseVector_8cpp.xml


// File: Epetra__SerialDenseVector_8h.xml


// File: Epetra__SerialDistributor_8cpp.xml


// File: Epetra__SerialDistributor_8h.xml


// File: Epetra__SerialSpdDenseSolver_8cpp.xml


// File: Epetra__SerialSpdDenseSolver_8h.xml


// File: Epetra__SerialSymDenseMatrix_8cpp.xml


// File: Epetra__SerialSymDenseMatrix_8h.xml


// File: Epetra__SrcDistObject_8h.xml


// File: Epetra__Time_8cpp.xml


// File: Epetra__Time_8h.xml


// File: Epetra__Util_8cpp.xml
%feature("docstring")  Epetra_Util_binary_search "int
Epetra_Util_binary_search(T item, const T *list, int len, int
&insertPoint)

Utility function to perform a binary-search on a list of data.
Important assumption: data is assumed to be sorted.

Parameters:
-----------

item:  to be searched for

list:  to be searched in

len:  Length of list

insertPoint:  Input/Output. If item is found, insertPoint is not
referenced. If item is not found, insertPoint is set to the offset at
which item should be inserted in list such that order (sortedness)
would be maintained.

offset Location in list at which item was found. -1 if not found. ";

%feature("docstring")  Epetra_Util_binary_search "int
Epetra_Util_binary_search(int item, const int *list, int len, int
&insertPoint) ";

%feature("docstring")  Epetra_Util_binary_search "int
Epetra_Util_binary_search(long long item, const long long *list, int
len, int &insertPoint) ";

%feature("docstring")  Epetra_Util_ExtractHbData "int
Epetra_Util_ExtractHbData(Epetra_CrsMatrix *A, Epetra_MultiVector
*LHS, Epetra_MultiVector *RHS, int &M, int &N, int &nz, int *&ptr, int
*&ind, double *&val, int &Nrhs, double *&rhs, int &ldrhs, double
*&lhs, int &ldlhs)

Harwell-Boeing data extraction routine.

This routine will extract data from an existing Epetra_Crs Matrix, and
optionally from related rhs and lhs objects in a form that is
compatible with software that requires the Harwell-Boeing data format.
The matrix must be passed in, but the RHS and LHS arguments may be set
to zero (either or both of them). For each of the LHS or RHS
arguments, if non-trivial and contain more than one vector, the
vectors must have strided access. If both LHS and RHS are non-trivial,
they must have the same number of vectors. If the input objects are
distributed, the returned matrices will contain the local part of the
matrix and vectors only.

Parameters:
-----------

A:  (In) Epetra_CrsMatrix.

LHS:  (In) Left hand side multivector. Set to zero if none not
available or needed.

RHS:  (In) Right hand side multivector. Set to zero if none not
available or needed.

M:  (Out) Local row dimension of matrix.

N:  (Out) Local column dimension of matrix.

nz:  (Out) Number of nonzero entries in matrix.

ptr:  (Out) Offsets into ind and val arrays pointing to start of each
row's data.

ind:  (Out) Column indices of the matrix, in compressed form.

val:  (Out) Matrix values, in compressed form corresponding to the ind
array.

Nrhs:  (Out) Number of right/left hand sides found (if any) in RHS and
LHS.

rhs:  (Out) Fortran-style 2D array of RHS values.

ldrhs:  (Out) Stride between columns of rhs.

lhs:  (Out) Fortran-style 2D array of LHS values.

ldrhs:  (Out) Stride between columns of lhs. ";


// File: Epetra__Util_8h.xml
%feature("docstring")  Epetra_Util_binary_search "int
Epetra_Util_binary_search(T item, const T *list, int len, int
&insertPoint)

Utility function to perform a binary-search on a list of data.
Important assumption: data is assumed to be sorted.

Parameters:
-----------

item:  to be searched for

list:  to be searched in

len:  Length of list

insertPoint:  Input/Output. If item is found, insertPoint is not
referenced. If item is not found, insertPoint is set to the offset at
which item should be inserted in list such that order (sortedness)
would be maintained.

offset Location in list at which item was found. -1 if not found. ";

%feature("docstring")  Epetra_Util_binary_search "EPETRA_LIB_DLL_EXPORT int Epetra_Util_binary_search(int item, const
int *list, int len, int &insertPoint) ";

%feature("docstring")  Epetra_Util_binary_search "EPETRA_LIB_DLL_EXPORT int Epetra_Util_binary_search(long long item,
const long long *list, int len, int &insertPoint) ";

%feature("docstring")  Epetra_Util_insert_empty_positions "int
Epetra_Util_insert_empty_positions(T *&array, int &usedLength, int
&allocatedLength, int insertOffset, int numPositions, int
allocChunkSize=32) ";

%feature("docstring")  Epetra_Util_insert "int Epetra_Util_insert(T
item, int offset, T *&list, int &usedLength, int &allocatedLength, int
allocChunkSize=32)

Function to insert an item in a list, at a specified offset.

Parameters:
-----------

item:  to be inserted

offset:  location at which to insert item

list:  array into which item is to be inserted. This array may be re-
allocated by this function.

usedLength:  number of items already present in list. Will be updated
to reflect the new length.

allocatedLength:  current allocated length of list. Will be updated to
reflect the new allocated-length, if applicable. Re-allocation occurs
only if usedLength==allocatedLength on entry.

allocChunkSize:  Optional argument, defaults to 32. Increment by which
the array should be expanded, if re-allocation is necessary.

error-code 0 if successful. -1 if input parameters don't make sense.
";

%feature("docstring")  Epetra_Util_ExtractHbData "int
Epetra_Util_ExtractHbData(Epetra_CrsMatrix *A, Epetra_MultiVector
*LHS, Epetra_MultiVector *RHS, int &M, int &N, int &nz, int *&ptr, int
*&ind, double *&val, int &Nrhs, double *&rhs, int &ldrhs, double
*&lhs, int &ldlhs)

Harwell-Boeing data extraction routine.

This routine will extract data from an existing Epetra_Crs Matrix, and
optionally from related rhs and lhs objects in a form that is
compatible with software that requires the Harwell-Boeing data format.
The matrix must be passed in, but the RHS and LHS arguments may be set
to zero (either or both of them). For each of the LHS or RHS
arguments, if non-trivial and contain more than one vector, the
vectors must have strided access. If both LHS and RHS are non-trivial,
they must have the same number of vectors. If the input objects are
distributed, the returned matrices will contain the local part of the
matrix and vectors only.

Parameters:
-----------

A:  (In) Epetra_CrsMatrix.

LHS:  (In) Left hand side multivector. Set to zero if none not
available or needed.

RHS:  (In) Right hand side multivector. Set to zero if none not
available or needed.

M:  (Out) Local row dimension of matrix.

N:  (Out) Local column dimension of matrix.

nz:  (Out) Number of nonzero entries in matrix.

ptr:  (Out) Offsets into ind and val arrays pointing to start of each
row's data.

ind:  (Out) Column indices of the matrix, in compressed form.

val:  (Out) Matrix values, in compressed form corresponding to the ind
array.

Nrhs:  (Out) Number of right/left hand sides found (if any) in RHS and
LHS.

rhs:  (Out) Fortran-style 2D array of RHS values.

ldrhs:  (Out) Stride between columns of rhs.

lhs:  (Out) Fortran-style 2D array of LHS values.

ldrhs:  (Out) Stride between columns of lhs. ";


// File: Epetra__VbrMatrix_8cpp.xml


// File: Epetra__VbrMatrix_8h.xml


// File: Epetra__VbrRowMatrix_8h.xml


// File: Epetra__Vector_8cpp.xml


// File: Epetra__Vector_8h.xml


// File: Epetra__Version_8h.xml
%feature("docstring")  Epetra_Version "string Epetra_Version() ";


// File: dir_233cfbe96499141f251547db95e4fda3.xml


// File: dir_29eb5d4e506afeb59b852f05b8c1a238.xml



// File: index.xml

// File: classEpetra__BasicDirectory.xml
%feature("docstring") Epetra_BasicDirectory "

Epetra_BasicDirectory: This class allows Epetra_Map objects to
reference non-local elements.

For Epetra_BlockMap objects, a Epetra_Directory object must be created
to allow referencing of non-local elements. The Epetra_BasicDirectory
produces and contains a uniform linear Epetra_BlockMap and a ProcList_
allowing blocks of non-local elements to be accessed by dereferencing
throught the Epetra_BasicDirectory.

This class currently has one constructor, taking a Epetra_BlockMap
object.

C++ includes: Epetra_BasicDirectory.h ";


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
Epetra_VbrRowMatrix for and example): Implement ExtractMyRowCopy and
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


// File: classEpetra__BlockMapData.xml
%feature("docstring") Epetra_BlockMapData "

Epetra_BlockMapData: The Epetra BlockMap Data Class.

The Epetra_BlockMapData class is an implementation detail of
Epetra_BlockMap. It is reference-counted, and can be shared by
multiple Epetra_BlockMap instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_BlockMapData.h ";


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


// File: classEpetra__CompObject.xml
%feature("docstring") Epetra_CompObject "

Epetra_CompObject: Functionality and data that is common to all
computational classes.

The Epetra_CompObject is a base class for all Epetra computational
objects. It provides the basic mechanisms and interface specifications
for floating point operations using Epetra_Flops objects.

C++ includes: Epetra_CompObject.h ";


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


// File: classEpetra__CrsGraphData.xml
%feature("docstring") Epetra_CrsGraphData "

Epetra_CrsGraphData: The Epetra CrsGraph Data Class.

The Epetra_CrsGraphData class is an implementation detail of
Epetra_CrsGraph. It is reference-counted, and can be shared by
multiple Epetra_CrsGraph instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_CrsGraphData.h ";


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


// File: classEpetra__Directory.xml
%feature("docstring") Epetra_Directory "

Epetra_Directory: This class is a pure virtual class whose interface
allows Epetra_Map and Epetr_BlockMap objects to reference non-local
elements.

For Epetra_BlockMap objects, a Epetra_Directory object must be created
by a call to the Epetra_Comm CreateDirectory method. The Directory is
needed to allow referencing of non-local elements.

C++ includes: Epetra_Directory.h ";


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
owning 3 elements as follows PE 0 Elements |  PE 1 Elements  |  PE 2
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

PE 0 Elements    |  PE 1 Elements    |  PE 2 Elements      0  1  2 3
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
x_force with the combined results of
x_integrate:x_force.Export(x_integrate, exporter, Add); The third
argument above tells the export operation to add results that come
from multiple processors for the same GID.

Epetra_Export objects can also be used by Import operations to perform
the reverse operation. For example, if x_force in the above example
had boundary conditions that should be sent to processors that share a
boundary element, the following operation would send replicated values
to x_integrate:x_integrate.Import(x_force, exporter, Insert); At the
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
ignoreNonLocalEntries=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, int NumIndicesPerRow, bool
ignoreNonLocalEntries=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
*NumIndicesPerRow, bool ignoreNonLocalEntries=false)

Constructor ";

%feature("docstring")  Epetra_FECrsGraph::Epetra_FECrsGraph "Epetra_FECrsGraph::Epetra_FECrsGraph(Epetra_DataAccess CV, const
Epetra_BlockMap &RowMap, const Epetra_BlockMap &ColMap, int
NumIndicesPerRow, bool ignoreNonLocalEntries=false)

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

%feature("docstring")  Epetra_FECrsMatrix::Epetra_FECrsMatrix "Epetra_FECrsMatrix::Epetra_FECrsMatrix(const Epetra_FECrsMatrix &src)

Copy Constructor. ";

%feature("docstring")  Epetra_FECrsMatrix::~Epetra_FECrsMatrix "Epetra_FECrsMatrix::~Epetra_FECrsMatrix()

Destructor. ";

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

%feature("docstring")  Epetra_FECrsMatrix::GlobalAssemble "int
Epetra_FECrsMatrix::GlobalAssemble(bool callFillComplete=true)

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


// File: classEpetra__FEVector.xml
%feature("docstring") Epetra_FEVector "

Epetra Finite-Element Vector. This class inherits Epetra_MultiVector
and thus provides all Epetra_MultiVector functionality, with one
restriction: currently an Epetra_FEVector only has 1 internal vector.

The added functionality provided by Epetra_FEVector is the ability to
perform finite-element style vector assembly. It accepts sub-vector
contributions, such as those that would come from element-load
vectors, etc., and these sub-vectors need not be wholly locally owned.
In other words, the user can assemble overlapping data (e.g.,
corresponding to shared finite-element nodes). When the user is
finished assembling their vector data, they then call the method
Epetra_FEVector::GlobalAssemble() which gathers the overlapping data
(all non-local data that was input on each processor) into the data-
distribution specified by the map that the Epetra_FEVector is
constructed with.

Note: At the current time (Sept 6, 2002) the methods in this
implementation assume that there is only 1 point associated with each
map element. This limitation will be removed in the near future.

C++ includes: Epetra_FEVector.h ";

%feature("docstring")  Epetra_FEVector::Epetra_FEVector "Epetra_FEVector::Epetra_FEVector(const Epetra_BlockMap &Map, bool
ignoreNonLocalEntries=false)

Constructor that requires a map specifying a non-overlapping data
layout. The methods SumIntoGlobalValues() and ReplaceGlobalValues()
will accept any global IDs, and GlobalAssemble() will move any non-
local data onto the appropriate owning processors. ";

%feature("docstring")  Epetra_FEVector::Epetra_FEVector "Epetra_FEVector::Epetra_FEVector(const Epetra_FEVector &source)

Copy constructor. ";

%feature("docstring")  Epetra_FEVector::~Epetra_FEVector "Epetra_FEVector::~Epetra_FEVector()

Destructor ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int *GIDs,
const double *values)

Accumulate values into the vector, adding them to any values that
already exist for the specified indices. ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(const Epetra_IntSerialDenseVector
&GIDs, const Epetra_SerialDenseVector &values)

Accumulate values into the vector, adding them to any values that
already exist for the specified GIDs.

Parameters:
-----------

GIDs:  List of global ids. Must be the same length as the accompanying
list of values.

values:  List of coefficient values. Must be the same length as the
accompanying list of GIDs. ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int *GIDs,
const double *values)

Copy values into the vector overwriting any values that already exist
for the specified indices. ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(const Epetra_IntSerialDenseVector
&GIDs, const Epetra_SerialDenseVector &values)

Copy values into the vector, replacing any values that already exist
for the specified GIDs.

Parameters:
-----------

GIDs:  List of global ids. Must be the same length as the accompanying
list of values.

values:  List of coefficient values. Must be the same length as the
accompanying list of GIDs. ";

%feature("docstring")  Epetra_FEVector::SumIntoGlobalValues "int
Epetra_FEVector::SumIntoGlobalValues(int numIDs, const int *GIDs,
const int *numValuesPerID, const double *values) ";

%feature("docstring")  Epetra_FEVector::ReplaceGlobalValues "int
Epetra_FEVector::ReplaceGlobalValues(int numIDs, const int *GIDs,
const int *numValuesPerID, const double *values) ";

%feature("docstring")  Epetra_FEVector::GlobalAssemble "int
Epetra_FEVector::GlobalAssemble(Epetra_CombineMode mode=Add)

Gather any overlapping/shared data into the non-overlapping
partitioning defined by the Map that was passed to this vector at
construction time. Data imported from other processors is stored on
the owning processor with a \"sumInto\" or accumulate operation. This
is a collective method -- every processor must enter it before any
will complete it. ";

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

%feature("docstring")  Epetra_Flops::Epetra_Flops "Epetra_Flops::Epetra_Flops(const Epetra_Flops &Flops)

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

%feature("docstring")  Epetra_HashTable::Epetra_HashTable "Epetra_HashTable::Epetra_HashTable(const int size, const unsigned int
seed=(2654435761U)) ";

%feature("docstring")  Epetra_HashTable::Epetra_HashTable "Epetra_HashTable::Epetra_HashTable(const Epetra_HashTable &obj) ";

%feature("docstring")  Epetra_HashTable::~Epetra_HashTable "Epetra_HashTable::~Epetra_HashTable() ";

%feature("docstring")  Epetra_HashTable::Add "void
Epetra_HashTable::Add(const int key, const int value) ";

%feature("docstring")  Epetra_HashTable::Get "int
Epetra_HashTable::Get(const int key) ";


// File: structEpetra__HashTable_1_1Node.xml


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
owning 3 elements as follows PE 0 Elements |  PE 1 Elements  |  PE 2
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
this, we build a target map on each processor as follows:    PE 0
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

NumSameIDs      = 0  NumPermuteIDs   = 3 PermuteToLIDs   = [0, 1, 2]
PermuteFromLIDs = [1, 2, 3]  NumRemoteIDs    = 2 RemoteLIDs      = [0,
4]  NumExportIDs    = 2 ExportLIDs      = [0, 2] ExportPIDs      = [0,
2]  NumSend         = 2 NumRecv         = 2

On PE 2:

NumSameIDs      = 0  NumPermuteIDs   = 3 PermuteToLIDs   = [0, 1, 2]
PermuteFromLIDs = [2, 3, 4]  NumRemoteIDs    = 2 RemoteLIDs      = [0,
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

%feature("docstring")
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector()

Default constructor; defines a zero size object.

Epetra_IntSerialDenseVector objects defined by the default constructor
should be sized with the Size() or Resize functions. Values should be
defined by using the [] or () operators. ";

%feature("docstring")
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector "Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(int Length)

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
CV, int *Values, int Length)

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
Epetra_IntSerialDenseVector::Size(int Length)

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
Epetra_IntSerialDenseVector::Resize(int Length)

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


// File: classEpetra__LinearProblem.xml
%feature("docstring") Epetra_LinearProblem "

Epetra_LinearProblem: The Epetra Linear Problem Class.

The Epetra_LinearProblem class is a wrapper that encapsulates the
general information needed for solving a linear system of equations.
Currently it accepts a Epetra matrix, initial guess and RHS and
returns the solution. the elapsed time for each calling processor.

C++ includes: Epetra_LinearProblem.h ";


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

%feature("docstring")  Epetra_LocalMap::Epetra_LocalMap "Epetra_LocalMap::Epetra_LocalMap(const Epetra_LocalMap &map)

Epetra_LocalMap copy constructor. ";

%feature("docstring")  Epetra_LocalMap::~Epetra_LocalMap "Epetra_LocalMap::~Epetra_LocalMap()

Epetra_LocalMap destructor. ";


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


// File: structEpetra__MapColoring_1_1ListItem.xml


// File: classEpetra__MpiComm.xml
%feature("docstring") Epetra_MpiComm "

Epetra_MpiComm: The Epetra MPI Communication Class.

The Epetra_MpiComm class is an implementation of Epetra_Comm that
encapsulates the general information and services needed for other
Epetra classes to run on a parallel computer using MPI.

C++ includes: Epetra_MpiComm.h ";


// File: classEpetra__MpiCommData.xml
%feature("docstring") Epetra_MpiCommData "

Epetra_MpiCommData: The Epetra Mpi Communication Data Class.

The Epetra_MpiCommData class is an implementation detail of
Epetra_MpiComm. It is reference-counted, and can be shared by multiple
Epetra_MpiComm instances. It derives from Epetra_Data, and inherits
reference-counting from it.

C++ includes: Epetra_MpiCommData.h ";


// File: classEpetra__MpiDistributor.xml
%feature("docstring") Epetra_MpiDistributor "

Epetra_MpiDistributor: The Epetra MPI implementation of the
Epetra_Distributor Gather/Scatter Setup Class.

The Epetra_MpiDistributor class is an MPI implement of
Epetra_Distributor that encapsulates the general information and
services needed for other Epetra classes to perform gather/scatter
operations on a parallel computer. An Epetra_MpiDistributor object is
actually produced by calling a method in the Epetra_MpiComm class.

C++ includes: Epetra_MpiDistributor.h ";


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

C++ includes: Epetra_MpiSmpComm.h ";


// File: classEpetra__MpiSmpCommData.xml
%feature("docstring") Epetra_MpiSmpCommData "

Epetra_MpiSmpCommData: The Epetra Mpi Shared Memory
ParallelCommunication Data Class.

The Epetra_MpiSmpCommData class is an implementation detail of
Epetra_MpiSmpComm. It is reference-counted, and can be shared by
multiple Epetra_MpiSmpComm instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_MpiSmpCommData.h ";


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


// File: classEpetra__OffsetIndex.xml
%feature("docstring") Epetra_OffsetIndex "

Epetra_OffsetIndex: This class builds index for efficient mapping of
data from one Epetra_CrsGraph based object to another.

Epetra_OffsetIndex generates and index of offsets allowing direct
access to data for Import/Export operations on Epetra_CrsGraph based
objects such as Epetra_CrsMatrix.

C++ includes: Epetra_OffsetIndex.h ";

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


// File: classEpetra__RowMatrixTransposer.xml
%feature("docstring") Epetra_RowMatrixTransposer "

Epetra_RowMatrixTransposer: A class for transposing an
Epetra_RowMatrix object.

This class provides capabilities to construct a transpose matrix of an
existing Epetra_RowMatrix object and (optionally) redistribute it
across a parallel distributed memory machine.

C++ includes: Epetra_RowMatrixTransposer.h ";


// File: classEpetra__SerialComm.xml
%feature("docstring") Epetra_SerialComm "

Epetra_SerialComm: The Epetra Serial Communication Class.

The Epetra_SerialComm class is an implementation of Epetra_Comm,
providing the general information and services needed for other Epetra
classes to run on a serial computer.

C++ includes: Epetra_SerialComm.h ";


// File: classEpetra__SerialCommData.xml
%feature("docstring") Epetra_SerialCommData "

Epetra_SerialCommData: The Epetra Serial Communication Data Class.

The Epetra_SerialCommData class is an implementation detail of
Epetra_SerialComm. It is reference-counted, and can be shared by
multiple Epetra_SerialComm instances. It derives from Epetra_Data, and
inherits reference-counting from it.

C++ includes: Epetra_SerialCommData.h ";


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
executing Waits). ";

%feature("docstring")  Epetra_SerialDistributor::DoWaits "int
Epetra_SerialDistributor::DoWaits()

Wait on a set of posts. ";

%feature("docstring")  Epetra_SerialDistributor::DoReversePosts "int
Epetra_SerialDistributor::DoReversePosts(char *export_objs, int
obj_size, int &len_import_objs, char *&import_objs)

Do reverse post of buffer of export objects (can do other local work
before executing Waits). ";

%feature("docstring")  Epetra_SerialDistributor::DoReverseWaits "int
Epetra_SerialDistributor::DoReverseWaits()

Wait on a reverse set of posts. ";

%feature("docstring")  Epetra_SerialDistributor::Do "int
Epetra_SerialDistributor::Do(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Execute plan on buffer of export objects in a single step (object size
may vary). ";

%feature("docstring")  Epetra_SerialDistributor::DoReverse "int
Epetra_SerialDistributor::DoReverse(char *export_objs, int obj_size,
int *&sizes, int &len_import_objs, char *&import_objs)

Execute reverse of plan on buffer of export objects in a single step
(object size may vary). ";

%feature("docstring")  Epetra_SerialDistributor::DoPosts "int
Epetra_SerialDistributor::DoPosts(char *export_objs, int obj_size, int
*&sizes, int &len_import_objs, char *&import_objs)

Post buffer of export objects (can do other local work before
executing Waits). ";

%feature("docstring")  Epetra_SerialDistributor::DoReversePosts "int
Epetra_SerialDistributor::DoReversePosts(char *export_objs, int
obj_size, int *&sizes, int &len_import_objs, char *&import_objs)

Do reverse post of buffer of export objects (can do other local work
before executing Waits). ";

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

%feature("docstring")  epetra_map_nummyelements "int MANGLE()
epetra_map_nummyelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_numglobalelements "int MANGLE()
epetra_map_numglobalelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_myglobalelements "int* MANGLE()
epetra_map_myglobalelements(EPETRA_OBJECT_REF map) ";

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

%feature("docstring")  epetra_map_nummyelements "int MANGLE()
epetra_map_nummyelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_numglobalelements "int MANGLE()
epetra_map_numglobalelements(EPETRA_OBJECT_REF map) ";

%feature("docstring")  epetra_map_myglobalelements "int* MANGLE()
epetra_map_myglobalelements(EPETRA_OBJECT_REF map) ";

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


// File: Epetra__LinearProblem_8cpp.xml


// File: Epetra__LinearProblem_8h.xml


// File: Epetra__LinearProblemRedistor_8cpp.xml


// File: Epetra__LinearProblemRedistor_8h.xml


// File: Epetra__LocalMap_8cpp.xml


// File: Epetra__LocalMap_8h.xml


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
Epetra_Util_binary_search(int item, const int *list, int len, int
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
Epetra_Util_binary_search(int item, const int *list, int len, int
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

%feature("docstring")  Epetra_Util_insert "int Epetra_Util_insert(T
item, int offset, T *&list, int &usedLength, int &allocatedLength, int
allocChunkSize=32)

Function to insert an item in a list, at a specified offset. error-
code 0 if successful, -1 if input parameters seem unreasonable (offset
> usedLength, offset<0, etc).

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


// File: dir_18b2d64510239fed06b88e74196cfd3f.xml


// File: dir_4368af47e412e90c65d06ecb9459c00d.xml


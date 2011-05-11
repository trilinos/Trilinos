
// File: index.xml

// File: classEpetraExt_1_1CrsGraph__SymmRCM_1_1BFT.xml


// File: classEpetraExt_1_1BlockAdjacencyGraph.xml
%feature("docstring") EpetraExt::BlockAdjacencyGraph "";

%feature("docstring")
EpetraExt::BlockAdjacencyGraph::~BlockAdjacencyGraph "EpetraExt::BlockAdjacencyGraph::~BlockAdjacencyGraph()

Destructor ";

%feature("docstring")
EpetraExt::BlockAdjacencyGraph::BlockAdjacencyGraph "EpetraExt::BlockAdjacencyGraph::BlockAdjacencyGraph()

Constructor ";

%feature("docstring")  EpetraExt::BlockAdjacencyGraph::compute "Teuchos::RCP< Epetra_CrsGraph >
EpetraExt::BlockAdjacencyGraph::compute(Epetra_CrsGraph &B, int nbrr,
std::vector< int > &r, std::vector< double > &weights, bool
verbose=false)

Constructs an adjacency graph representing the block connectivity of
the input graph, where nbrr is the number of block rows in B and r
contains the row index where each block begins. A reference-counted
pointer to an Epetra_CrsGraph that has nbrr rows is returned as well
as the vector of weights. This vector is of length nbrr returns some
weighting on the block adjacency graph that can be used to balance the
original graph B. Right now, that weighting is just the number of rows
in each block. ";


// File: classEpetraExt_1_1CrsGraph__AMD.xml
%feature("docstring") EpetraExt::CrsGraph_AMD "

EpetraExt::CrsGraph_AMD: A transform for Approximate Minimum Degree
Reordering using Tim Daley's AMD Algorithm.

C++ includes: EpetraExt_AMD_CrsGraph.h ";

%feature("docstring")  EpetraExt::CrsGraph_AMD::~CrsGraph_AMD "EpetraExt::CrsGraph_AMD::~CrsGraph_AMD()

EpetraExt::CrsGraph_AMD Destructor. ";

%feature("docstring")  EpetraExt::CrsGraph_AMD::CrsGraph_AMD "EpetraExt::CrsGraph_AMD::CrsGraph_AMD(bool verbose=false)

EpetraExt::CrsGraph_AMD Constructor.

Creates a transform for AMD reordering of a Epetra_CrsGraph

Parameters:
-----------

In:  verbose - Turns on verbosity for debugging ";


// File: classEpetraExt_1_1CrsGraph__MapColoring.xml
%feature("docstring") EpetraExt::CrsGraph_MapColoring "

Map Coloring of independent columns in a Graph

Generates a Epetra_MapColoring object for which all column indices
form independent sets.

C++ includes: EpetraExt_MapColoring.h ";

%feature("docstring")
EpetraExt::CrsGraph_MapColoring::~CrsGraph_MapColoring "EpetraExt::CrsGraph_MapColoring::~CrsGraph_MapColoring()

Destructor ";

%feature("docstring")
EpetraExt::CrsGraph_MapColoring::CrsGraph_MapColoring "EpetraExt::CrsGraph_MapColoring::CrsGraph_MapColoring(ColoringAlgorithm
algo=GREEDY, int reordering=0, bool distance1=false, int verbosity=0)

Constructor ";


// File: classEpetraExt_1_1CrsGraph__MapColoringIndex.xml
%feature("docstring") EpetraExt::CrsGraph_MapColoringIndex "

Generates a std::vector of Epetra_IntVector's to be used to map
perturbation contributions to a CrsGraph/CrsMatrix from a perturbed
vector.

C++ includes: EpetraExt_MapColoringIndex.h ";

%feature("docstring")
EpetraExt::CrsGraph_MapColoringIndex::~CrsGraph_MapColoringIndex "EpetraExt::CrsGraph_MapColoringIndex::~CrsGraph_MapColoringIndex()

Destructor ";

%feature("docstring")
EpetraExt::CrsGraph_MapColoringIndex::CrsGraph_MapColoringIndex "EpetraExt::CrsGraph_MapColoringIndex::CrsGraph_MapColoringIndex(const
Epetra_MapColoring &ColorMap)

Constructor input param ColorMap defines the perturbation coloring ";


// File: classEpetraExt_1_1CrsGraph__Overlap.xml
%feature("docstring") EpetraExt::CrsGraph_Overlap "

Given an input Epetra_CrsGraph, a \"overlapped\" Epetra_CrsGraph is
generated including rows associated with off processor contributions.

C++ includes: EpetraExt_Overlap_CrsGraph.h ";

%feature("docstring")  EpetraExt::CrsGraph_Overlap::~CrsGraph_Overlap
"EpetraExt::CrsGraph_Overlap::~CrsGraph_Overlap()

Destructor ";

%feature("docstring")  EpetraExt::CrsGraph_Overlap::CrsGraph_Overlap "EpetraExt::CrsGraph_Overlap::CrsGraph_Overlap(int overlap, bool
squareLocalBlock=false)

Constructor ";


// File: structEpetraExt_1_1CrsGraph__SymmRCM.xml
%feature("docstring") EpetraExt::CrsGraph_SymmRCM "

Generates the symmetric RCM reordered version of a Epetra_CrsGraph.

C++ includes: EpetraExt_SymmRCM_CrsGraph.h ";

%feature("docstring")  EpetraExt::CrsGraph_SymmRCM::~CrsGraph_SymmRCM
"EpetraExt::CrsGraph_SymmRCM::~CrsGraph_SymmRCM()

Destructor. ";

%feature("docstring")  EpetraExt::CrsGraph_SymmRCM::CrsGraph_SymmRCM "EpetraExt::CrsGraph_SymmRCM::CrsGraph_SymmRCM(bool BruteForce=false,
int testLeafWidth=5)

Constructor. ";


// File: classEpetraExt_1_1CrsGraph__Transpose.xml
%feature("docstring") EpetraExt::CrsGraph_Transpose "

Transform to generate the explicit transpose of a Epetra_CrsGraph.

C++ includes: EpetraExt_Transpose_CrsGraph.h ";

%feature("docstring")
EpetraExt::CrsGraph_Transpose::~CrsGraph_Transpose "EpetraExt::CrsGraph_Transpose::~CrsGraph_Transpose()

Destructor. ";

%feature("docstring")
EpetraExt::CrsGraph_Transpose::CrsGraph_Transpose "EpetraExt::CrsGraph_Transpose::CrsGraph_Transpose(bool
IgnoreNonLocalCols=false)

Constructor.

Parameters:
-----------

In:  IgnoreNonLocalCols - Whether or not to include non-local columns
in the transpose. ";


// File: classEpetraExt_1_1CrsGraph__View.xml
%feature("docstring") EpetraExt::CrsGraph_View "

Generates a sub-block view of a Epetra_CrsGraph.

C++ includes: EpetraExt_View_CrsGraph.h ";

%feature("docstring")  EpetraExt::CrsGraph_View::~CrsGraph_View "EpetraExt::CrsGraph_View::~CrsGraph_View()

Destructor. ";

%feature("docstring")  EpetraExt::CrsGraph_View::CrsGraph_View "EpetraExt::CrsGraph_View::CrsGraph_View(const Epetra_BlockMap
*new_row_map, const Epetra_BlockMap *new_col_map=0)

Constructor. ";


// File: classEpetraExt_1_1CrsMatrix__Dirichlet.xml
%feature("docstring") EpetraExt::CrsMatrix_Dirichlet "

Given an input Epetra_LinearProblem, apply given dirichlet conditions

C++ includes: EpetraExt_Dirichlet_CrsMatrix.h ";

%feature("docstring")
EpetraExt::CrsMatrix_Dirichlet::~CrsMatrix_Dirichlet "EpetraExt::CrsMatrix_Dirichlet::~CrsMatrix_Dirichlet()

Destructor ";

%feature("docstring")
EpetraExt::CrsMatrix_Dirichlet::CrsMatrix_Dirichlet "EpetraExt::CrsMatrix_Dirichlet::CrsMatrix_Dirichlet(const
Epetra_IntVector &Locations, bool Symmetric=false)

Constructor

Parameters:
-----------

Locations:  Integer vector containing 1's for Dirichlet BC rows and
0's otherwise

Symmetric:  Boolean flag indicating whether to enforce symmetry by
zeroing out columns, false by default ";

%feature("docstring")  EpetraExt::CrsMatrix_Dirichlet::fwd "bool
EpetraExt::CrsMatrix_Dirichlet::fwd()

Applies Dirichlet BC's ";

%feature("docstring")  EpetraExt::CrsMatrix_Dirichlet::rvs "bool
EpetraExt::CrsMatrix_Dirichlet::rvs()

NoOp ";


// File: classEpetraExt_1_1CrsMatrix__Reindex.xml
%feature("docstring") EpetraExt::CrsMatrix_Reindex "

Given an Epetra_CrsMatrix, a \"reindexed\" version is returned based
on the new row map. The row map must be conformal to the original. The
Matrix data will be shared by the new Matrix using the new indexing

C++ includes: EpetraExt_Reindex_CrsMatrix.h ";

%feature("docstring")
EpetraExt::CrsMatrix_Reindex::~CrsMatrix_Reindex "EpetraExt::CrsMatrix_Reindex::~CrsMatrix_Reindex()

Destructor ";

%feature("docstring")  EpetraExt::CrsMatrix_Reindex::CrsMatrix_Reindex
"EpetraExt::CrsMatrix_Reindex::CrsMatrix_Reindex(const Epetra_Map
&new_row_map)

Constructor ";


// File: classEpetraExt_1_1CrsMatrix__SolverMap.xml
%feature("docstring") EpetraExt::CrsMatrix_SolverMap "

Given an input Epetra_CrsMatrix, the column map is checked for missing
indices associated with the local rows. If found, a view of the
Epetra_CrsMatrix is formed using a new Epetra_CrsGraph with a fixed
column mapping including all local row indices.

C++ includes: EpetraExt_SolverMap_CrsMatrix.h ";

%feature("docstring")
EpetraExt::CrsMatrix_SolverMap::~CrsMatrix_SolverMap "EpetraExt::CrsMatrix_SolverMap::~CrsMatrix_SolverMap()

Destructor ";

%feature("docstring")
EpetraExt::CrsMatrix_SolverMap::CrsMatrix_SolverMap "EpetraExt::CrsMatrix_SolverMap::CrsMatrix_SolverMap()

Constructor ";


// File: classEpetraExt_1_1CrsMatrix__SubCopy.xml
%feature("docstring") EpetraExt::CrsMatrix_SubCopy "

Generates a sub-block view of a Epetra_CrsMatrix.

C++ includes: EpetraExt_SubCopy_CrsMatrix.h ";

%feature("docstring")
EpetraExt::CrsMatrix_SubCopy::~CrsMatrix_SubCopy "EpetraExt::CrsMatrix_SubCopy::~CrsMatrix_SubCopy()

Destructor. ";

%feature("docstring")  EpetraExt::CrsMatrix_SubCopy::CrsMatrix_SubCopy
"EpetraExt::CrsMatrix_SubCopy::CrsMatrix_SubCopy(const Epetra_Map
&newMap)

Constructor. ";

%feature("docstring")  EpetraExt::CrsMatrix_SubCopy::CrsMatrix_SubCopy
"EpetraExt::CrsMatrix_SubCopy::CrsMatrix_SubCopy(const Epetra_Map
&newRangeAndRowMap, const Epetra_Map &newDomainMap)

Constructor. ";

%feature("docstring")  EpetraExt::CrsMatrix_SubCopy::fwd "bool
EpetraExt::CrsMatrix_SubCopy::fwd()

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful.

Preconditions: ";

%feature("docstring")  EpetraExt::CrsMatrix_SubCopy::rvs "bool
EpetraExt::CrsMatrix_SubCopy::rvs()

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful.

Preconditions: ";


// File: classEpetraExt_1_1CrsMatrix__View.xml
%feature("docstring") EpetraExt::CrsMatrix_View "

Generates a sub-block view of a Epetra_CrsMatrix.

C++ includes: EpetraExt_View_CrsMatrix.h ";

%feature("docstring")  EpetraExt::CrsMatrix_View::~CrsMatrix_View "EpetraExt::CrsMatrix_View::~CrsMatrix_View()

Destructor. ";

%feature("docstring")  EpetraExt::CrsMatrix_View::CrsMatrix_View "EpetraExt::CrsMatrix_View::CrsMatrix_View(const Epetra_CrsGraph
&orig_graph, const Epetra_CrsGraph &new_graph)

Constructor. ";


// File: classEpetraExt_1_1CrsMatrixStruct.xml
%feature("docstring") EpetraExt::CrsMatrixStruct "";

%feature("docstring")  EpetraExt::CrsMatrixStruct::CrsMatrixStruct "EpetraExt::CrsMatrixStruct::CrsMatrixStruct() ";

%feature("docstring")  EpetraExt::CrsMatrixStruct::~CrsMatrixStruct "EpetraExt::CrsMatrixStruct::~CrsMatrixStruct() ";

%feature("docstring")  EpetraExt::CrsMatrixStruct::deleteContents "void EpetraExt::CrsMatrixStruct::deleteContents() ";


// File: classEpetraExt_1_1CrsWrapper.xml
%feature("docstring") EpetraExt::CrsWrapper "";

%feature("docstring")  EpetraExt::CrsWrapper::~CrsWrapper "virtual
EpetraExt::CrsWrapper::~CrsWrapper() ";

%feature("docstring")  EpetraExt::CrsWrapper::RowMap "virtual const
Epetra_Map& EpetraExt::CrsWrapper::RowMap() const =0 ";

%feature("docstring")  EpetraExt::CrsWrapper::Filled "virtual bool
EpetraExt::CrsWrapper::Filled()=0 ";

%feature("docstring")  EpetraExt::CrsWrapper::InsertGlobalValues "virtual int EpetraExt::CrsWrapper::InsertGlobalValues(int GlobalRow,
int NumEntries, double *Values, int *Indices)=0 ";

%feature("docstring")  EpetraExt::CrsWrapper::SumIntoGlobalValues "virtual int EpetraExt::CrsWrapper::SumIntoGlobalValues(int GlobalRow,
int NumEntries, double *Values, int *Indices)=0 ";


// File: classEpetraExt_1_1CrsWrapper__Epetra__CrsMatrix.xml
%feature("docstring") EpetraExt::CrsWrapper_Epetra_CrsMatrix "";

%feature("docstring")
EpetraExt::CrsWrapper_Epetra_CrsMatrix::CrsWrapper_Epetra_CrsMatrix "EpetraExt::CrsWrapper_Epetra_CrsMatrix::CrsWrapper_Epetra_CrsMatrix(Epetra_CrsMatrix
&epetracrsmatrix) ";

%feature("docstring")
EpetraExt::CrsWrapper_Epetra_CrsMatrix::~CrsWrapper_Epetra_CrsMatrix "EpetraExt::CrsWrapper_Epetra_CrsMatrix::~CrsWrapper_Epetra_CrsMatrix()
";

%feature("docstring")  EpetraExt::CrsWrapper_Epetra_CrsMatrix::RowMap
"const Epetra_Map & EpetraExt::CrsWrapper_Epetra_CrsMatrix::RowMap()
const ";

%feature("docstring")  EpetraExt::CrsWrapper_Epetra_CrsMatrix::Filled
"bool EpetraExt::CrsWrapper_Epetra_CrsMatrix::Filled() ";

%feature("docstring")
EpetraExt::CrsWrapper_Epetra_CrsMatrix::InsertGlobalValues "int
EpetraExt::CrsWrapper_Epetra_CrsMatrix::InsertGlobalValues(int
GlobalRow, int NumEntries, double *Values, int *Indices) ";

%feature("docstring")
EpetraExt::CrsWrapper_Epetra_CrsMatrix::SumIntoGlobalValues "int
EpetraExt::CrsWrapper_Epetra_CrsMatrix::SumIntoGlobalValues(int
GlobalRow, int NumEntries, double *Values, int *Indices) ";


// File: classEpetraExt_1_1CrsWrapper__GraphBuilder.xml
%feature("docstring") EpetraExt::CrsWrapper_GraphBuilder "";

%feature("docstring")
EpetraExt::CrsWrapper_GraphBuilder::CrsWrapper_GraphBuilder "EpetraExt::CrsWrapper_GraphBuilder::CrsWrapper_GraphBuilder(const
Epetra_Map &emap) ";

%feature("docstring")
EpetraExt::CrsWrapper_GraphBuilder::~CrsWrapper_GraphBuilder "EpetraExt::CrsWrapper_GraphBuilder::~CrsWrapper_GraphBuilder() ";

%feature("docstring")  EpetraExt::CrsWrapper_GraphBuilder::RowMap "const Epetra_Map& EpetraExt::CrsWrapper_GraphBuilder::RowMap() const
";

%feature("docstring")  EpetraExt::CrsWrapper_GraphBuilder::Filled "bool EpetraExt::CrsWrapper_GraphBuilder::Filled() ";

%feature("docstring")
EpetraExt::CrsWrapper_GraphBuilder::InsertGlobalValues "int
EpetraExt::CrsWrapper_GraphBuilder::InsertGlobalValues(int GlobalRow,
int NumEntries, double *Values, int *Indices) ";

%feature("docstring")
EpetraExt::CrsWrapper_GraphBuilder::SumIntoGlobalValues "int
EpetraExt::CrsWrapper_GraphBuilder::SumIntoGlobalValues(int GlobalRow,
int NumEntries, double *Values, int *Indices) ";

%feature("docstring")  EpetraExt::CrsWrapper_GraphBuilder::get_graph "std::map< int, std::set< int > * > &
EpetraExt::CrsWrapper_GraphBuilder::get_graph() ";

%feature("docstring")
EpetraExt::CrsWrapper_GraphBuilder::get_max_row_length "int
EpetraExt::CrsWrapper_GraphBuilder::get_max_row_length() ";


// File: classEpetraExt_1_1DistArray.xml
%feature("docstring") EpetraExt::DistArray "

DistArray<T>: A class to store row-oriented multivectors of type T.

Class DistArray allows the construction and usage of multivectors.
These vectors contain element of type T, and the storage is row-
oriented, and not column-oriented as in class Epetra_MultiVector. As
such, this class should be used as a container for data, on which no
BLAS-like operations are performed.

DistArray objects are indentified by an Epetra_Map and a RowSize. The
map specifies the distribution of the elements across the processors
and therefore the number of local elements, while the RowSize gives
the total number of data assigned to each node. RowSize is constant
for all elements.

DistArray is derived from Epetra_DistObject, and it can therefore be
redistributed using Import/Export instructions.

The typical usage of this class is as follows:

Marzio Sala, ETHZ/D-INFK.

C++ includes: EpetraExt_DistArray.h ";

%feature("docstring")  EpetraExt::DistArray::DistArray "EpetraExt::DistArray< T >::DistArray(const Epetra_Map &Map, const int
RowSize)

Constructor for a given Map and RowSize. ";

%feature("docstring")  EpetraExt::DistArray::MyLength "int
EpetraExt::DistArray< T >::MyLength() const

Returns the length of the locally owned array. ";

%feature("docstring")  EpetraExt::DistArray::GlobalLength "int
EpetraExt::DistArray< T >::GlobalLength() const

Returns the global length of the array. ";

%feature("docstring")  EpetraExt::DistArray::RowSize "int
EpetraExt::DistArray< T >::RowSize() const

Returns the row size, that is, the amount of data associated with each
element. ";

%feature("docstring")  EpetraExt::DistArray::Print "void
EpetraExt::DistArray< T >::Print(std::ostream &os) const

Prints the array on the specified stream. ";

%feature("docstring")  EpetraExt::DistArray::NextGID "int
EpetraExt::DistArray< T >::NextGID() ";

%feature("docstring")  EpetraExt::DistArray::FirstGID "int
EpetraExt::DistArray< T >::FirstGID() ";

%feature("docstring")  EpetraExt::DistArray::ExtractView "const
std::vector<T>& EpetraExt::DistArray< T >::ExtractView() const

Extracts a view of the array. ";

%feature("docstring")  EpetraExt::DistArray::Values "T*
EpetraExt::DistArray< T >::Values()

Returns a pointer to the internally stored data (non-const version).
";

%feature("docstring")  EpetraExt::DistArray::Values "const T*
EpetraExt::DistArray< T >::Values() const

Returns a pointer to the internally stored data (const version). ";


// File: classEpetraExt_1_1Exception.xml
%feature("docstring") EpetraExt::Exception "";

%feature("docstring")  EpetraExt::Exception::Exception "EpetraExt::Exception::Exception(const string FileName, const int
LineNumber, const string Line1, const string Line2=\"\", const string
Line3=\"\", const string Line4=\"\", const string Line5=\"\", const
string Line6=\"\") ";

%feature("docstring")  EpetraExt::Exception::Print "void
EpetraExt::Exception::Print() ";


// File: classEpetraExt_1_1InPlaceTransform.xml
%feature("docstring") EpetraExt::InPlaceTransform "";

%feature("docstring")  EpetraExt::InPlaceTransform::~InPlaceTransform
"virtual EpetraExt::InPlaceTransform< T >::~InPlaceTransform() ";


// File: classEpetraExt_1_1LinearProblem__CrsSingletonFilter.xml
%feature("docstring") EpetraExt::LinearProblem_CrsSingletonFilter "

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

C++ includes: EpetraExt_CrsSingletonFilter_LinearProblem.h ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::LinearProblem_CrsSingletonFilter
"EpetraExt::LinearProblem_CrsSingletonFilter::LinearProblem_CrsSingletonFilter(bool
verbose=false)

Epetra_CrsSingletonFilter default constructor. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::~LinearProblem_CrsSingletonFilter
"EpetraExt::LinearProblem_CrsSingletonFilter::~LinearProblem_CrsSingletonFilter()

Epetra_CrsSingletonFilter Destructor. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::Analyze "int
EpetraExt::LinearProblem_CrsSingletonFilter::Analyze(Epetra_RowMatrix
*FullMatrix)

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

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::SingletonsDetected "bool
EpetraExt::LinearProblem_CrsSingletonFilter::SingletonsDetected()
const

Returns true if singletons were detected in this matrix (must be
called after Analyze() to be effective). ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ConstructReducedProblem "int
EpetraExt::LinearProblem_CrsSingletonFilter::ConstructReducedProblem(Epetra_LinearProblem
*Problem)

Return a reduced linear problem based on results of Analyze().

Creates a new Epetra_LinearProblem object based on the results of the
Analyze phase. A pointer to the reduced problem is obtained via a call
to ReducedProblem().

Error code, set to 0 if no error. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::UpdateReducedProblem "int
EpetraExt::LinearProblem_CrsSingletonFilter::UpdateReducedProblem(Epetra_LinearProblem
*Problem)

Update a reduced linear problem using new values.

Updates an existing Epetra_LinearProblem object using new matrix, LHS
and RHS values. The matrix structure must be identical to the matrix
that was used to construct the original reduced problem.

Error code, set to 0 if no error. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ComputeFullSolution "int
EpetraExt::LinearProblem_CrsSingletonFilter::ComputeFullSolution()

Compute a solution for the full problem using the solution of the
reduced problem, put in LHS of FullProblem().

After solving the reduced linear system, this method can be called to
compute the solution to the original problem, assuming the solution
for the reduced system is valid. The solution of the unreduced,
original problem will be in the LHS of the original
Epetra_LinearProblem. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::NumRowSingletons "int
EpetraExt::LinearProblem_CrsSingletonFilter::NumRowSingletons() const

Return number of rows that contain a single entry, returns -1 if
Analysis not performed yet. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::NumColSingletons "int
EpetraExt::LinearProblem_CrsSingletonFilter::NumColSingletons() const

Return number of columns that contain a single entry that are not
associated with singleton row, returns -1 if Analysis not performed
yet. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::NumSingletons "int
EpetraExt::LinearProblem_CrsSingletonFilter::NumSingletons() const

Return total number of singletons detected, returns -1 if Analysis not
performed yet.

Return total number of singletons detected across all processors. This
method will not return a valid result until after the Analyze() method
is called. The dimension of the reduced system can be computed by
subtracting this number from dimension of full system. WARNING:  This
method returns -1 if Analyze() method has not been called. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::RatioOfDimensions "double
EpetraExt::LinearProblem_CrsSingletonFilter::RatioOfDimensions() const

Returns ratio of reduced system to full system dimensions, returns
-1.0 if reduced problem not constructed. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::RatioOfNonzeros "double
EpetraExt::LinearProblem_CrsSingletonFilter::RatioOfNonzeros() const

Returns ratio of reduced system to full system nonzero count, returns
-1.0 if reduced problem not constructed. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::FullProblem "Epetra_LinearProblem*
EpetraExt::LinearProblem_CrsSingletonFilter::FullProblem() const

Returns pointer to the original unreduced Epetra_LinearProblem. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedProblem "Epetra_LinearProblem*
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedProblem() const

Returns pointer to the derived reduced Epetra_LinearProblem. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::FullMatrix "Epetra_RowMatrix*
EpetraExt::LinearProblem_CrsSingletonFilter::FullMatrix() const

Returns pointer to Epetra_CrsMatrix from full problem. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrix "Epetra_CrsMatrix*
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrix() const

Returns pointer to Epetra_CrsMatrix from full problem. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::RowMapColors "Epetra_MapColoring*
EpetraExt::LinearProblem_CrsSingletonFilter::RowMapColors() const

Returns pointer to Epetra_MapColoring object: color 0 rows are part of
reduced system. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ColMapColors "Epetra_MapColoring*
EpetraExt::LinearProblem_CrsSingletonFilter::ColMapColors() const

Returns pointer to Epetra_MapColoring object: color 0 columns are part
of reduced system. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixRowMap "Epetra_Map*
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixRowMap()
const

Returns pointer to Epetra_Map describing the reduced system row
distribution. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixColMap "Epetra_Map*
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixColMap()
const

Returns pointer to Epetra_Map describing the reduced system column
distribution. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixDomainMap "Epetra_Map*
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixDomainMap()
const

Returns pointer to Epetra_Map describing the domain map for the
reduced system. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixRangeMap "Epetra_Map*
EpetraExt::LinearProblem_CrsSingletonFilter::ReducedMatrixRangeMap()
const

Returns pointer to Epetra_Map describing the range map for the reduced
system. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::analyze "bool
EpetraExt::LinearProblem_CrsSingletonFilter::analyze(OriginalTypeRef
orig)

Initial analysis phase of transform. Returns true if the transform is
possible allowing methods  construct(),  fwd() and  rvs() to be
successfully utilized.

Preconditions: default implementation calls method operator() and
stores the resulting object in an internal attribute newObj_. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::construct "LinearProblem_CrsSingletonFilter::NewTypeRef
EpetraExt::LinearProblem_CrsSingletonFilter::construct()

Construction of new object as a result of the transform.

Preconditions: default implementation returns internal attribute
newObj_. ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::fwd "bool
EpetraExt::LinearProblem_CrsSingletonFilter::fwd()

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful.

Preconditions: ";

%feature("docstring")
EpetraExt::LinearProblem_CrsSingletonFilter::rvs "bool
EpetraExt::LinearProblem_CrsSingletonFilter::rvs()

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful.

Preconditions: ";


// File: classEpetraExt_1_1LinearProblem__GraphTrans.xml
%feature("docstring") EpetraExt::LinearProblem_GraphTrans "

EpetraExt::LinearProblem_GraphTrans: Adaptation of a Epetra_CrsGraph
Transform to a Epetra_LinearProblem Transform.

C++ includes: EpetraExt_LPTrans_From_GraphTrans.h ";

%feature("docstring")
EpetraExt::LinearProblem_GraphTrans::~LinearProblem_GraphTrans "EpetraExt::LinearProblem_GraphTrans::~LinearProblem_GraphTrans()

EpetraExt::LinearProblem_GraphTrans Destructor. ";

%feature("docstring")
EpetraExt::LinearProblem_GraphTrans::LinearProblem_GraphTrans "EpetraExt::LinearProblem_GraphTrans::LinearProblem_GraphTrans(StructuralSameTypeTransform<
Epetra_CrsGraph > &graph_trans)

EpetraExt::LinearProblem_GraphTrans Constructor.

Constructs a LinearProblem Transform based on the input CrsGraph
Transform

Parameters:
-----------

In:  graph_trans - Base Epetra_CrsGraph Transform from which a
consistent Epetra_LinearProblem Transform is generated ";

%feature("docstring")  EpetraExt::LinearProblem_GraphTrans::fwd "bool
EpetraExt::LinearProblem_GraphTrans::fwd()

Forward migration of data from original to transformed object. ";

%feature("docstring")  EpetraExt::LinearProblem_GraphTrans::rvs "bool
EpetraExt::LinearProblem_GraphTrans::rvs()

Reverse migration of data from transformed to original object. ";


// File: classEpetraExt_1_1LinearProblem__MatrixTrans.xml
%feature("docstring") EpetraExt::LinearProblem_MatrixTrans "

Adaptation of an Epetra_CrsMatrix Transform to a Epetra_LinearProblem
Transform.

C++ includes: EpetraExt_LPTrans_From_MatrixTrans.h ";

%feature("docstring")
EpetraExt::LinearProblem_MatrixTrans::~LinearProblem_MatrixTrans "EpetraExt::LinearProblem_MatrixTrans::~LinearProblem_MatrixTrans()

Destructor. ";

%feature("docstring")
EpetraExt::LinearProblem_MatrixTrans::LinearProblem_MatrixTrans "EpetraExt::LinearProblem_MatrixTrans::LinearProblem_MatrixTrans(SameTypeTransform<
Epetra_CrsMatrix > &matrix_trans)

Constructor. ";

%feature("docstring")  EpetraExt::LinearProblem_MatrixTrans::fwd "bool EpetraExt::LinearProblem_MatrixTrans::fwd()

Forward Data Migration. ";

%feature("docstring")  EpetraExt::LinearProblem_MatrixTrans::rvs "bool EpetraExt::LinearProblem_MatrixTrans::rvs()

Reverse Data Migration. ";


// File: classEpetraExt_1_1LinearProblem__Reindex.xml
%feature("docstring") EpetraExt::LinearProblem_Reindex "

Given and input Epetra_LinearProblem, a \"reindexed\" version will be
returned using the given NewRowMap. If a null map is given, a
lexigraphically indexed LP will be returned. The data in the new E_LP
is a \"reindexed\" view of the original.

C++ includes: EpetraExt_Reindex_LinearProblem.h ";

%feature("docstring")
EpetraExt::LinearProblem_Reindex::~LinearProblem_Reindex "EpetraExt::LinearProblem_Reindex::~LinearProblem_Reindex()

Destructor ";

%feature("docstring")
EpetraExt::LinearProblem_Reindex::LinearProblem_Reindex "EpetraExt::LinearProblem_Reindex::LinearProblem_Reindex(Epetra_Map
*NewRowMap)

Constructor ";


// File: classEpetraExt_1_1LinearProblem__Scale.xml
%feature("docstring") EpetraExt::LinearProblem_Scale "

Given an input Epetra_LinearProblem, recursive, left and right scaling
are performed.

C++ includes: EpetraExt_Scale_LinearProblem.h ";

%feature("docstring")
EpetraExt::LinearProblem_Scale::~LinearProblem_Scale "EpetraExt::LinearProblem_Scale::~LinearProblem_Scale()

Destructor ";

%feature("docstring")
EpetraExt::LinearProblem_Scale::LinearProblem_Scale "EpetraExt::LinearProblem_Scale::LinearProblem_Scale(ScaleType
left=Sum, ScaleType right=Sum, double exp_fac=1.0, int iterations=1)

Constructor ";

%feature("docstring")  EpetraExt::LinearProblem_Scale::fwd "bool
EpetraExt::LinearProblem_Scale::fwd()

Applies forward scaling ";

%feature("docstring")  EpetraExt::LinearProblem_Scale::rvs "bool
EpetraExt::LinearProblem_Scale::rvs()

Reverses scaling ";


// File: classEpetraExt_1_1LinearProblem__SolverMap.xml
%feature("docstring") EpetraExt::LinearProblem_SolverMap "

Constructs a LinearProblem with a \"fixed\" Column Map for the
CrsMatrix. Almost entirely a view except for the \"fixed\"
Epetra_CrsGraph.

C++ includes: EpetraExt_SolverMap_LinearProblem.h ";

%feature("docstring")
EpetraExt::LinearProblem_SolverMap::~LinearProblem_SolverMap "EpetraExt::LinearProblem_SolverMap::~LinearProblem_SolverMap()

Destructor ";

%feature("docstring")
EpetraExt::LinearProblem_SolverMap::LinearProblem_SolverMap "EpetraExt::LinearProblem_SolverMap::LinearProblem_SolverMap()

Constructor ";


// File: classEpetraExt_1_1MatrixMatrix.xml
%feature("docstring") EpetraExt::MatrixMatrix "

Collection of matrix-matrix operations. This class basically functions
as a namespace, containing only static methods. See the program
epetraext/test/MatrixMatrix/cxx_main.cpp for a usage example.

C++ includes: EpetraExt_MatrixMatrix.h ";

%feature("docstring")  EpetraExt::MatrixMatrix::~MatrixMatrix "virtual EpetraExt::MatrixMatrix::~MatrixMatrix()

destructor ";


// File: classEpetraExt_1_1MultiVector__Reindex.xml
%feature("docstring") EpetraExt::MultiVector_Reindex "

Given an input Epetra_MultiVector, a \"reindexed\" view is returned.

C++ includes: EpetraExt_Reindex_MultiVector.h ";

%feature("docstring")
EpetraExt::MultiVector_Reindex::~MultiVector_Reindex "EpetraExt::MultiVector_Reindex::~MultiVector_Reindex()

Destructor ";

%feature("docstring")
EpetraExt::MultiVector_Reindex::MultiVector_Reindex "EpetraExt::MultiVector_Reindex::MultiVector_Reindex(const Epetra_Map
&new_row_map)

Constructor ";


// File: classEpetraExt_1_1MultiVector__View.xml
%feature("docstring") EpetraExt::MultiVector_View "

Generates a sub-block view of a Epetra_MultiVector.

C++ includes: EpetraExt_View_MultiVector.h ";

%feature("docstring")  EpetraExt::MultiVector_View::~MultiVector_View
"EpetraExt::MultiVector_View::~MultiVector_View()

Destructor. ";

%feature("docstring")  EpetraExt::MultiVector_View::MultiVector_View "EpetraExt::MultiVector_View::MultiVector_View(const Epetra_BlockMap
&orig_map, const Epetra_BlockMap &new_map, const int num_vec=-1)

Constructor. ";


// File: structEpetraExt_1_1Perm__traits.xml
%feature("docstring") EpetraExt::Perm_traits "

Define some traits to make it easier to deal with template-parameters
which are objects to be permuted. Given a template parameter, we'll
want to have the following operations available: determine the type

construct an instance of it

replace its row-map

produce a column-permutation of it

First the default definition, which catches all types \"T\", followed
by some specializations for anticipated types. Any type other than the
types specifically anticipated will be handled by this default
definition, allowing the Permutation class to abort or return NULL
where appropriate.

We define these trait structs in this file rather than in a separate
file in an attempt to avoid some template-instantiation
complications...

C++ includes: EpetraExt_Permutation_impl.h ";


// File: structEpetraExt_1_1Perm__traits_3_01Epetra__CrsGraph_01_4.xml
%feature("docstring") EpetraExt::Perm_traits< Epetra_CrsGraph > "

A specialization of Perm_traits for the specific type Epetra_CrsGraph.

C++ includes: EpetraExt_Permutation_impl.h ";


// File: structEpetraExt_1_1Perm__traits_3_01Epetra__CrsMatrix_01_4.xml
%feature("docstring") EpetraExt::Perm_traits< Epetra_CrsMatrix > "

A specialization of Perm_traits for the specific type
Epetra_CrsMatrix.

C++ includes: EpetraExt_Permutation_impl.h ";


// File: structEpetraExt_1_1Perm__traits_3_01Epetra__MultiVector_01_4.xml
%feature("docstring") EpetraExt::Perm_traits< Epetra_MultiVector > "

A specialization of Perm_traits for the specific type
Epetra_MultiVector.

C++ includes: EpetraExt_Permutation_impl.h ";


// File: classEpetraExt_1_1Permutation.xml
%feature("docstring") EpetraExt::Permutation "

Permutation stores and describes a permutation matrix P. As described
in \"Matrix Computations\" (Golub and Van Loan), a permutation matrix
is the identity matrix with its rows re-ordered. The permutation is
internally stored as an integer vector p, where p[i] is the column-
index of the \"1\" in P's i-th row. Consider the example of permuting
a matrix A by applying the permutation matrix P to form the result B.
i.e., B = PA. If p[i] = j, then row j of A becomes row i of B.

This Permutation class is templated on the type of the object to be
permuted. However, not all objects are eligible to be template
parameters. Currently the following objects may be used:
Epetra_CrsMatrix, Epetra_CrsGraph and Epetra_MultiVector.

A test program which exercises this Permutation class is located in
the directory packages/epetraext/test/Permutation.

Implementation Notes: Permutation currently inherits
StructuralSameTypeTransform, which in turn     inherits Transform
through SameTypeTransform. Permutation, and its base classes,     are
templates. A couple of noteworthy consequences result from this:

1. A separate instantiation of Permutation must be created for each
type        of object to be permuted. Example:
Epetra_CrsGraph& graph = ... Epetra_CrsMatrix& A = ...
Permutation<Epetra_CrsGraph> graph_perm(...);
Permutation<Epetra_CrsMatrix> matrix_perm(...);

Epetra_CrsMatrix& PA = matrix_perm(A); Epetra_CrsGraph& Pgraph =
graph_perm(graph);

2. Following the semantics of Transform, when the Permutation class is
used        to create a new permuted copy of an object, ownership of
the new copy is        retained by Permutation. Permutation will
destroy the new object. This means        that only one object should
be permuted by a Permutation instance.

It is not clear that these are desirable behaviors for permutations.
It is     possible that Permutation will be altered to remove these
limitations, as     follows:

1. If Permutation doesn't inherit Transform, then Permutation need not
be        a template and instead we could either overload or template-
ize the        operator() method member. This would allow a single
instantiation of Permutation to be used for permuting all of the
eligible target types.

2. Allowing the caller (user) to take ownership of the newly- produced
permuted objects would allow a single Permutation instance to be used
repeatedly since it would no longer need to hold a pointer to the new
object        for later deletion.

Then, example usage could look like this: Epetra_CrsMatrix& A = ...
Epetra_MultiVector& v = ... Permutation P(...);

Epetra_CrsMatrix PA = P(A);          Epetra_MultiVector Pv = P(v);

Questions and comments about this class may be directed to Alan
Williams.

C++ includes: EpetraExt_Permutation.h ";

%feature("docstring")  EpetraExt::Permutation::Permutation "EpetraExt::Permutation< T >::Permutation(Epetra_DataAccess CV, const
Epetra_BlockMap &map, int *permutation)

Constructor

Parameters:
-----------

CV:  Set to either Copy or View.

map:  Defines the index space to be permuted.

permutation:  Array defining the permutation. The length of this array
must be 'map.NumMyElements()'. This array is the local portion of the
'p' vector described in the 'Detailed Description' section. ";

%feature("docstring")  EpetraExt::Permutation::Permutation "EpetraExt::Permutation< T >::Permutation(const Epetra_BlockMap &map)

Constructor. This constructor creates an empty permutation object. The
contents must then be set using regular Epetra_IntVector methods.

Parameters:
-----------

map:  Defines the index space to be permuted. ";

%feature("docstring")  EpetraExt::Permutation::Permutation "EpetraExt::Permutation< T >::Permutation(const Permutation< T > &src)

Copy Constructor ";

%feature("docstring")  EpetraExt::Permutation::~Permutation "EpetraExt::Permutation< T >::~Permutation()

Destructor ";


// File: classEpetraExt_1_1ProductOperator.xml
%feature("docstring") EpetraExt::ProductOperator "

Implements Epetra_Operator as a product of one or more Epetra_Operator
objects.

This class implements a product operator of the form:

M = M[0]*M[1]*...*M[num_Op-1]

and operator applications are performed one constituent operator at a
time as:

Forward Mat-vec: Y = M * X     T[k-1] = M[k]*T[k]       for k =
num_Op-1...0         where: T[num_Op-1] = X (input vector) where:
T[-1]       = Y (output vector)   Adjoint Mat-vec: Y = M' * X T[k] =
M[k]'*T[k-1]       for k = 0...num_Op-1         where: T[-1] = X
(input vector)        where: T[num_Op-1] = Y (output vector)

Likewise, the inverse can also be applied (if all of the constituent
operators support the inverse operation) as:

Forward Inverse Mat-vec: Y = inv(M) * X     T[k] = inv(M[k])*T[k-1]
for k = 0...num_Op-1       for k = 0...num_Op-1 where: T[-1]       = X
(input vector)        where: T[num_Op-1] = Y (output vector)   Adjoint
Inverse Mat-vec: Y = inv(M') * X     T[k] = inv(M[k]')*T[k-1]
for k = num_Op-1...0         where: T[num_Op-1] = X (input vector)
where: T[-1]       = Y (output vector)

Note that maps for the result of the inverse of an operator is the
same as the result of the adjoint of the operator and the map for the
result of the inverse of the adjoint is the same as for the result the
non-inverse forward opeator.

The client constructs this object with a list of Epetra_Operator
objects an how the non-transposed operator is to be viewed and if it
is to be views as its inverse or not (see  initialize()).

Note: The Epetra_Map objects returned from OperatorDomainMap() and
OperatorRangeMap() must always be with respect to the non-transposed
operator! This is very strange behavior and is totally undocumented in
the Epetra_Operator interface but it seems to be the case.

C++ includes: EpetraExt_ProductOperator.h ";

/*  Public types  */

/*  Constructors / initializers / accessors  */

%feature("docstring")  EpetraExt::ProductOperator::ProductOperator "EpetraExt::ProductOperator::ProductOperator()

Construct to uninitialized. ";

%feature("docstring")  EpetraExt::ProductOperator::ProductOperator "EpetraExt::ProductOperator::ProductOperator(const int num_Op, const
Teuchos::RefCountPtr< const Epetra_Operator > Op[], const
Teuchos::ETransp Op_trans[], const EApplyMode Op_inverse[])

Calls  initialize(). ";

%feature("docstring")  EpetraExt::ProductOperator::initialize "void
EpetraExt::ProductOperator::initialize(const int num_Op, const
Teuchos::RefCountPtr< const Epetra_Operator > Op[], const
Teuchos::ETransp Op_trans[], const EApplyMode Op_inverse[])

Setup with constituent operators.

Parameters:
-----------

num_Op:  [in] Number of constinuent operators.

Op:  [in] Array (length num_Op) of smart pointers to Epetra_Operator
objects that implement each constituent operator.

Op_trans:    [in] Array (length num_Op) that determines the transpose
mode of the constituent operator (see below).

Op_inverse:    [in] Array (length num_Op) that determines the inverse
apply mode of the constituent operator (see below).

Preconditions:   num_Op > 0

Op!=NULL

Op[k].get()!=NULL, for k=0...num_Op-1

Op_trans!=NULL

Op_inverse!=NULL

Postconditions:  this-> num_Op()==num_Op

this->Op(k).get()==Op[k].get(), for k=0...num_Op-1

this->Op_trans(k)==Op_trans[k], for k=0...num_Op-1

this->Op_inverse(k)==Op_inverse[k], for k=0...num_Op-1

The forward constituent operator T[k-1] = M[k]*T[k] described in the
main documenatation above is defined as follows:

Op[k]->SetUseTranspose( Op_trans[k]!=Teuchos::NO_TRANS ); if(
Op_inverse[k]==APPLY_MODE_APPLY )      Op[k]->Apply( T[k], T[k-1] );
else         Op[k]->ApplyInverse( T[k], T[k-1] );

The inverse constituent operator T[k] = inv(M[k])*T[k-1] described in
the main documenatation above is defined as follows:

Op[k]->SetUseTranspose( Op_trans[k]!=Teuchos::NO_TRANS ); if(
Op_inverse[k]==APPLY_MODE_APPLY )      Op[k]->ApplyInverse( T[k-1],
T[k] );          else         Op[k]->Apply( T[k-1], T[k] );

The other transposed constituent operators M[k]' and inv(M[k]') are
defined by simply changing the value of the transpose as
Op[k]->SetUseTranspose( Op_trans[k]==Teuchos::NO_TRANS );. Note, that
Op[k]->SetUseTranspose(...) is called immediately before
Op[k]->Apply(...) or Op[k]->ApplyInverse(...) is called to avoid
tricky problems that could occur with multiple uses of the same
operator. ";

%feature("docstring")  EpetraExt::ProductOperator::uninitialize "void
EpetraExt::ProductOperator::uninitialize(int *num_Op,
Teuchos::RefCountPtr< const Epetra_Operator > Op[], Teuchos::ETransp
Op_trans[], EApplyMode p_inverse[])

Set to an uninitialized state and wipe out memory.

Postconditions:  this-> num_Op()==0 ";

%feature("docstring")  EpetraExt::ProductOperator::applyConstituent "void EpetraExt::ProductOperator::applyConstituent(const int k,
Teuchos::ETransp Op_trans, EApplyMode Op_inverse, const
Epetra_MultiVector &X_k, Epetra_MultiVector *Y_k) const

Apply the kth aggregate operator M[k] correctly.

Parameters:
-----------

k:  [in] Gives the index (zero-based) of the constituent operator to
apply.

Op_trans:  [in] Determines if the transpose of the constituent
operator is to be applied.

Op_inverse:  [in] Determines if the operator or its inverse (if
supported) should be applied.

X_k:  [in] The input vector.

Y_k:  [out] The output vector.

Clients should call this function to correctly apply a constitient
operator! ";

%feature("docstring")  EpetraExt::ProductOperator::num_Op "int
EpetraExt::ProductOperator::num_Op() const

Return the number of aggregate opeators.

A return value of 0 is a flag that this is not initialized yet. ";

%feature("docstring")  EpetraExt::ProductOperator::Op "Teuchos::RefCountPtr< const Epetra_Operator >
EpetraExt::ProductOperator::Op(int k) const

Access the kth operator (zero-based).

Preconditions:  0 <= k <= this-> num_Op()-1

Warning! This is the raw opeator passed into initialize(...). In order
to apply the constituent operator M[k] you must call
ApplyConstituent(). ";

%feature("docstring")  EpetraExt::ProductOperator::Op_trans "Teuchos::ETransp EpetraExt::ProductOperator::Op_trans(int k) const

Access the transpose mode of the kth operator (zero-based).

Preconditions:  0 <= k <= this-> num_Op()-1 ";

%feature("docstring")  EpetraExt::ProductOperator::Op_inverse "ProductOperator::EApplyMode EpetraExt::ProductOperator::Op_inverse(int
k) const

Access the inverse mode of the kth operator (zero-based).

Preconditions:  0 <= k <= this-> num_Op()-1 ";

/*  Overridden from Epetra_Operator  */

%feature("docstring")  EpetraExt::ProductOperator::SetUseTranspose "int EpetraExt::ProductOperator::SetUseTranspose(bool UseTranspose) ";

%feature("docstring")  EpetraExt::ProductOperator::Apply "int
EpetraExt::ProductOperator::Apply(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  EpetraExt::ProductOperator::ApplyInverse "int
EpetraExt::ProductOperator::ApplyInverse(const Epetra_MultiVector &X,
Epetra_MultiVector &Y) const ";

%feature("docstring")  EpetraExt::ProductOperator::NormInf "double
EpetraExt::ProductOperator::NormInf() const ";

%feature("docstring")  EpetraExt::ProductOperator::Label "const char
* EpetraExt::ProductOperator::Label() const ";

%feature("docstring")  EpetraExt::ProductOperator::UseTranspose "bool
EpetraExt::ProductOperator::UseTranspose() const ";

%feature("docstring")  EpetraExt::ProductOperator::HasNormInf "bool
EpetraExt::ProductOperator::HasNormInf() const ";

%feature("docstring")  EpetraExt::ProductOperator::Comm "const
Epetra_Comm & EpetraExt::ProductOperator::Comm() const ";

%feature("docstring")  EpetraExt::ProductOperator::OperatorDomainMap "const Epetra_Map & EpetraExt::ProductOperator::OperatorDomainMap()
const ";

%feature("docstring")  EpetraExt::ProductOperator::OperatorRangeMap "const Epetra_Map & EpetraExt::ProductOperator::OperatorRangeMap()
const ";


// File: classEpetraExt_1_1RowMatrix__Transpose.xml
%feature("docstring") EpetraExt::RowMatrix_Transpose "

Transform to form the explicit transpose of a Epetra_RowMatrix.

C++ includes: EpetraExt_Transpose_RowMatrix.h ";

%feature("docstring")
EpetraExt::RowMatrix_Transpose::~RowMatrix_Transpose "EpetraExt::RowMatrix_Transpose::~RowMatrix_Transpose()

Destructor. ";

%feature("docstring")
EpetraExt::RowMatrix_Transpose::RowMatrix_Transpose "EpetraExt::RowMatrix_Transpose::RowMatrix_Transpose(bool
MakeDataContiguous=false, Epetra_Map *TransposeRowMap=0, bool
IgnoreNonLocalCols=false)

Constructor.

Parameters:
-----------

In:  MakeDataContiguous - Whether to optimize form of matrix to be
contiguous data storage.

In:  TransposeRowMap - Map to be used for row mapping of transpose
matrix

In:  IgnoreNonLocalCols - Whether to ignore non-local columns for the
transpose ";

%feature("docstring")  EpetraExt::RowMatrix_Transpose::fwd "bool
EpetraExt::RowMatrix_Transpose::fwd()

Foward Data Migration. ";

%feature("docstring")  EpetraExt::RowMatrix_Transpose::rvs "bool
EpetraExt::RowMatrix_Transpose::rvs()

Reverse Data Migration. ";


// File: classEpetraExt_1_1SameTypeTransform.xml
%feature("docstring") EpetraExt::SameTypeTransform "";

%feature("docstring")
EpetraExt::SameTypeTransform::~SameTypeTransform "virtual
EpetraExt::SameTypeTransform< T >::~SameTypeTransform() ";


// File: classEpetraExt_1_1StructuralSameTypeTransform.xml
%feature("docstring") EpetraExt::StructuralSameTypeTransform "";

%feature("docstring")  EpetraExt::StructuralSameTypeTransform::fwd "bool EpetraExt::StructuralSameTypeTransform< T >::fwd()

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful.

Preconditions: ";

%feature("docstring")  EpetraExt::StructuralSameTypeTransform::rvs "bool EpetraExt::StructuralSameTypeTransform< T >::rvs()

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful.

Preconditions: ";

%feature("docstring")
EpetraExt::StructuralSameTypeTransform::~StructuralSameTypeTransform "virtual EpetraExt::StructuralSameTypeTransform< T
>::~StructuralSameTypeTransform() ";


// File: classEpetraExt_1_1StructuralTransform.xml
%feature("docstring") EpetraExt::StructuralTransform "";

%feature("docstring")  EpetraExt::StructuralTransform::fwd "bool
EpetraExt::StructuralTransform< T, U >::fwd()

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful.

Preconditions: ";

%feature("docstring")  EpetraExt::StructuralTransform::rvs "bool
EpetraExt::StructuralTransform< T, U >::rvs()

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful.

Preconditions: ";

%feature("docstring")
EpetraExt::StructuralTransform::~StructuralTransform "virtual
EpetraExt::StructuralTransform< T, U >::~StructuralTransform() ";


// File: classEpetraExt_1_1Transform.xml
%feature("docstring") EpetraExt::Transform "

Base Class for all Epetra Transform Operators.

This is the abstract definition for all Epetra Transform Operators.
Depending on the type of Transform, several specializations are
available: Structural, SameType, InPlace, View.

C++ includes: EpetraExt_Transform.h ";

/*  Typedefs for templated classes  */

/*  Pure Virtual Methods which must be implemented by subclasses  */

%feature("docstring")  EpetraExt::Transform::fwd "virtual bool
EpetraExt::Transform< T, U >::fwd()=0

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful.

Preconditions: ";

%feature("docstring")  EpetraExt::Transform::rvs "virtual bool
EpetraExt::Transform< T, U >::rvs()=0

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful.

Preconditions: ";

/*  Virtual functions with default implements allowing for optional
*/

/* implementation by the Transform developer

*/

%feature("docstring")  EpetraExt::Transform::analyze "bool
EpetraExt::Transform< T, U >::analyze(OriginalTypeRef orig)

Initial analysis phase of transform. Returns true if the transform is
possible allowing methods  construct(),  fwd() and  rvs() to be
successfully utilized.

Preconditions: default implementation calls method operator() and
stores the resulting object in an internal attribute newObj_. ";

%feature("docstring")  EpetraExt::Transform::construct "Transform< T,
U >::NewTypeRef EpetraExt::Transform< T, U >::construct()

Construction of new object as a result of the transform.

Preconditions: default implementation returns internal attribute
newObj_. ";

%feature("docstring")  EpetraExt::Transform::isConstructed "bool
EpetraExt::Transform< T, U >::isConstructed()

Check for whether transformed object has been constructed

Preconditions: default implementation returns true if newObj_ != 0. ";

%feature("docstring")  EpetraExt::Transform::~Transform "virtual
EpetraExt::Transform< T, U >::~Transform() ";


// File: classEpetraExt_1_1Transform__Composite.xml
%feature("docstring") EpetraExt::Transform_Composite "

Composition Class for Epetra Transform SameType Operators.

This class allows SameType Transforms to be composed as a single
Transform.

C++ includes: EpetraExt_Transform_Composite.h ";

%feature("docstring")
EpetraExt::Transform_Composite::Transform_Composite "EpetraExt::Transform_Composite< T >::Transform_Composite()

EpetraExt::Transform_Composite Constructor. ";

%feature("docstring")
EpetraExt::Transform_Composite::~Transform_Composite "EpetraExt::Transform_Composite< T >::~Transform_Composite()

EpetraExt::Transform_Composite Destructor. ";

%feature("docstring")  EpetraExt::Transform_Composite::addTransform "void EpetraExt::Transform_Composite< T
>::addTransform(TransformTypePtr new_trans)

Transform Addition.

Add SameType Transform to composition. Order of Addition == Order of
Application ";

%feature("docstring")  EpetraExt::Transform_Composite::fwd "bool
EpetraExt::Transform_Composite< T >::fwd()

Forward Data Transfer.

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful. ";

%feature("docstring")  EpetraExt::Transform_Composite::rvs "bool
EpetraExt::Transform_Composite< T >::rvs()

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful. ";


// File: classEpetraExt_1_1Vector__Dirichlet.xml
%feature("docstring") EpetraExt::Vector_Dirichlet "

Given an input Epetra_Vector, apply given dirichlet conditions

C++ includes: EpetraExt_Dirichlet_Vector.h ";

%feature("docstring")  EpetraExt::Vector_Dirichlet::~Vector_Dirichlet
"EpetraExt::Vector_Dirichlet::~Vector_Dirichlet()

Destructor ";

%feature("docstring")  EpetraExt::Vector_Dirichlet::Vector_Dirichlet "EpetraExt::Vector_Dirichlet::Vector_Dirichlet(const Epetra_IntVector
&Locations, const Epetra_Vector &Values)

Constructor

Parameters:
-----------

Locations:  Integer Vector containing 1's for Dirichlet rows and 0's
if not

Values:  Vector containing values of the Dirichlet BC's ";

%feature("docstring")  EpetraExt::Vector_Dirichlet::fwd "bool
EpetraExt::Vector_Dirichlet::fwd()

Applies Dirichlet BC's ";

%feature("docstring")  EpetraExt::Vector_Dirichlet::rvs "bool
EpetraExt::Vector_Dirichlet::rvs()

NoOp ";


// File: classEpetraExt_1_1ViewTransform.xml
%feature("docstring") EpetraExt::ViewTransform "";

%feature("docstring")  EpetraExt::ViewTransform::fwd "bool
EpetraExt::ViewTransform< T >::fwd()

Forward transfer of data from orig object input in the operator()
method call to the new object created in this same call. Returns true
is operation is successful.

Preconditions: ";

%feature("docstring")  EpetraExt::ViewTransform::rvs "bool
EpetraExt::ViewTransform< T >::rvs()

Reverse transfer of data from new object created in the operator()
method call to the orig object input to this same method. Returns true
if operation is successful.

Preconditions: ";

%feature("docstring")  EpetraExt::ViewTransform::~ViewTransform "virtual EpetraExt::ViewTransform< T >::~ViewTransform() ";


// File: classEpetraExt_1_1XMLReader.xml
%feature("docstring") EpetraExt::XMLReader "

class XMLReader: A class for reading Epetra objects stored in XML
files.

Class EpetraExt::XMLReader allows to read several Trilinos objects
stored in XML files. The XML data format is specified in the
documentation of class EpetraExt::XMLWriter, which also contains a
MATLAB script. A typical usage of this class is reported in file
epetraext/example/inout/XML_IO.cpp.

This class requires Teuchos to be configured with the option --enable-
teuchos-expat.

Reading objects from a file requires the following steps. First, we
define an XMLReader object, Then, we define a set of pointers, Reading
simply goes as follows: In distributed environments,
Epetra_MultiVector, Epetra_CrsGraph and Epetra_CrsMatrix objects have
a linear distribution. Epetra_Map objects can be read only when using
the same number of processors used for writing.

WARNING:  All the created objects must be deleted from the user using
delete.

Marzio Sala, D-INFK/ETHZ

C++ includes: EpetraExt_XMLReader.h ";

%feature("docstring")  EpetraExt::XMLReader::XMLReader "EpetraExt::XMLReader::XMLReader(const Epetra_Comm &Comm, const
std::string &FileName)

ctor ";

%feature("docstring")  EpetraExt::XMLReader::~XMLReader "EpetraExt::XMLReader::~XMLReader()

dtor ";

%feature("docstring")  EpetraExt::XMLReader::Read "void
EpetraExt::XMLReader::Read(const std::string &Label, Epetra_Map *&Map)

Reads the Epetra_Map stored with label Label. ";

%feature("docstring")  EpetraExt::XMLReader::Read "void
EpetraExt::XMLReader::Read(const std::string &Label, Epetra_CrsGraph
*&Graph)

Reads the Epetra_CrsGraph stored with label Label. ";

%feature("docstring")  EpetraExt::XMLReader::Read "void
EpetraExt::XMLReader::Read(const std::string &Label, Epetra_CrsMatrix
*&Matrix)

Reads the Epetra_CrsMatrix stored with label Label. ";

%feature("docstring")  EpetraExt::XMLReader::Read "void
EpetraExt::XMLReader::Read(const std::string &Label,
Epetra_MultiVector *&MultiVector)

Reads the Epetra_MultiVector stored with label Label. ";

%feature("docstring")  EpetraExt::XMLReader::Read "void
EpetraExt::XMLReader::Read(const std::string &Label, std::vector<
std::string > &Content)

Reads a std::vector of strings with label Label. ";

%feature("docstring")  EpetraExt::XMLReader::Read "void
EpetraExt::XMLReader::Read(const std::string &Label,
Teuchos::ParameterList &List)

Reads the Teuchos::ParameterList stored with label Label. ";


// File: classEpetraExt_1_1XMLWriter.xml
%feature("docstring") EpetraExt::XMLWriter "

class XMLWriter: A class for writing Trilinos objects to XML files.

Class EpetraExt::XMLWriter writes several Trilinos objects in an XML-
compatible format. The list of supported objects contains: Epetra_Map;

Epetra_MultiVector;

Epetra_CrsGraph;

Epetra_CrsMatrix;

Epetra_RowMatrix;

Teuchos::ParameterList.

All objects can be read and written, with the std::exception of
Epetra_RowMatrix objects, that can only be written to files.

An example of usage is reported in file
epetraext/example/inout/XML_IO.cpp.

Writing objects goes as follows. Let Map, Matrix, LHS and RHS an
Epetra_Map, Epetra_CrsMatrix, and two Epetra_MultiVector's,
respectively. First, we define an XMLWriter object and we open the
file using MyProblem label: Writing objects simply goes as A
Teuchos::ParameterList (List), a std::string, and a
std::vector<std::string> can be written as Finally, we close the file
Note that only processor 0 writes the Teuchos::ParameterList,
std::string, and std::vector<std::string>.

The written file is as follows:

This class requires Teuchos to be configured with the option --enable-
teuchos-expat.

Marzio Sala, D-INFK/ETHZ

C++ includes: EpetraExt_XMLWriter.h ";

%feature("docstring")  EpetraExt::XMLWriter::XMLWriter "EpetraExt::XMLWriter::XMLWriter(const Epetra_Comm &Comm, const
std::string &FileName)

ctor ";

%feature("docstring")  EpetraExt::XMLWriter::~XMLWriter "EpetraExt::XMLWriter::~XMLWriter()

dtor ";

%feature("docstring")  EpetraExt::XMLWriter::Create "void
EpetraExt::XMLWriter::Create(const std::string &Label)

Creates the file, giving Label to the whole object. ";

%feature("docstring")  EpetraExt::XMLWriter::Close "void
EpetraExt::XMLWriter::Close()

Closes the file. No Write operations can follow. ";

%feature("docstring")  EpetraExt::XMLWriter::Write "void
EpetraExt::XMLWriter::Write(const std::string &Label, const Epetra_Map
&Map)

Writes an Epetra_Map using label Label. ";

%feature("docstring")  EpetraExt::XMLWriter::Write "void
EpetraExt::XMLWriter::Write(const std::string &Label, const
Epetra_RowMatrix &Matrix)

Writes an Epetra_RowMatrix using label Label. ";

%feature("docstring")  EpetraExt::XMLWriter::Write "void
EpetraExt::XMLWriter::Write(const std::string &Label, const
Epetra_MultiVector &MultiVector)

Writes an Epetra_MultiVector using label Label. ";

%feature("docstring")  EpetraExt::XMLWriter::Write "void
EpetraExt::XMLWriter::Write(const std::string &Label, const
std::vector< std::string > &Content)

Writes the std::vector of std::string's using label Label. ";

%feature("docstring")  EpetraExt::XMLWriter::Write "void
EpetraExt::XMLWriter::Write(const std::string &Label, const
std::string &Text)

Writes input std::string using label Label. ";

%feature("docstring")  EpetraExt::XMLWriter::Write "void
EpetraExt::XMLWriter::Write(const std::string &Label,
Teuchos::ParameterList &List)

Writes a Teuchos::ParameterList using label Label. ";


// File: namespaceEpetraExt.xml
%feature("docstring")  EpetraExt::sparsedot "double
EpetraExt::sparsedot(double *u, int *u_ind, int u_len, double *v, int
*v_ind, int v_len)

Method for internal use... sparsedot forms a dot-product between two
sparsely-populated 'vectors'. Important assumption: assumes the
indices in u_ind and v_ind are sorted. ";

%feature("docstring")  EpetraExt::mult_A_B "int
EpetraExt::mult_A_B(CrsMatrixStruct &Aview, CrsMatrixStruct &Bview,
CrsWrapper &C) ";

%feature("docstring")  EpetraExt::mult_A_Btrans "int
EpetraExt::mult_A_Btrans(CrsMatrixStruct &Aview, CrsMatrixStruct
&Bview, CrsWrapper &C) ";

%feature("docstring")  EpetraExt::mult_Atrans_B "int
EpetraExt::mult_Atrans_B(CrsMatrixStruct &Aview, CrsMatrixStruct
&Bview, CrsWrapper &C) ";

%feature("docstring")  EpetraExt::mult_Atrans_Btrans "int
EpetraExt::mult_Atrans_Btrans(CrsMatrixStruct &Aview, CrsMatrixStruct
&Bview, CrsWrapper &C) ";

%feature("docstring")  EpetraExt::import_and_extract_views "int
EpetraExt::import_and_extract_views(const Epetra_CrsMatrix &M, const
Epetra_Map &targetMap, CrsMatrixStruct &Mview) ";

%feature("docstring")  EpetraExt::distribute_list "int
EpetraExt::distribute_list(const Epetra_Comm &Comm, int lenSendList,
const int *sendList, int &maxSendLen, int *&recvList) ";

%feature("docstring")  EpetraExt::create_map_from_imported_rows "Epetra_Map* EpetraExt::create_map_from_imported_rows(const Epetra_Map
*map, int totalNumSend, int *sendRows, int numProcs, int
*numSendPerProc) ";

%feature("docstring")  EpetraExt::form_map_union "int
EpetraExt::form_map_union(const Epetra_Map *map1, const Epetra_Map
*map2, const Epetra_Map *&mapunion) ";

%feature("docstring")  EpetraExt::find_rows_containing_cols "Epetra_Map* EpetraExt::find_rows_containing_cols(const
Epetra_CrsMatrix &M, const Epetra_Map *colmap) ";

%feature("docstring")  EpetraExt::dumpCrsMatrixStruct "int
EpetraExt::dumpCrsMatrixStruct(const CrsMatrixStruct &M) ";

%feature("docstring")  EpetraExt::insert_matrix_locations "void
EpetraExt::insert_matrix_locations(CrsWrapper_GraphBuilder
&graphbuilder, Epetra_CrsMatrix &C) ";

%feature("docstring")  EpetraExt::EpetraExt_Version "std::string
EpetraExt::EpetraExt_Version() ";

%feature("docstring")  EpetraExt::MatrixMarketFileToMap "int
EpetraExt::MatrixMarketFileToMap(const char *filename, const
Epetra_Comm &comm, Epetra_Map *&map)

Constructs an Epetra_BlockMap object from a Matrix Market format file.

This function constructs an Epetra_BlockMap or Epetra_Map object by
reading a Matrix Market file. If the file was created using the
EpetraExt::BlockMapOut functions, special information was encoded in
the comment field of this map that allows for identical reproduction
of the map, including distribution across processors and element size
information. If the same of processors is being used to create the
object as were used to write it, the object will be an exact
reproduction of the original. Otherwise, a uniform distribution of the
GIDs will be created.

The first column of the input file will must be the list of GIDs in
the map. If the block map has non-uniform sizes, a second column must
contain the element sizes.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the map GIDs. A
second column may be present if the BlockMap has nonuniform sizes.

comm:  (In) An Epetra_Comm object describing the parallel machine.

map:  (Out) An Epetra_Map object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system, -2 if file
contained nontrivial Epetra_BlockMap, 1 if number of processors
differs from file creator. ";

%feature("docstring")  EpetraExt::MatrixMarketFileToBlockMap "int
EpetraExt::MatrixMarketFileToBlockMap(const char *filename, const
Epetra_Comm &comm, Epetra_BlockMap *&blockMap)

Constructs an Epetra_BlockMap object from a Matrix Market format file.

This function constructs an Epetra_BlockMap or Epetra_Map object by
reading a Matrix Market file. If the file was created using the
EpetraExt::BlockMapOut functions, special information was encoded in
the comment field of this map that allows for identical reproduction
of the map, including distribution across processors and element size
information. If the same of processors is being used to create the
object as were used to write it, the object will be an exact
reproduction of the original. Otherwise, a uniform distribution of the
GIDs will be created.

The first column of the input file will must be the list of GIDs in
the map. If the block map has non-uniform sizes, a second column must
contain the element sizes.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the map GIDs. A
second column may be present if the BlockMap has nonuniform sizes.

comm:  (In) An Epetra_Comm object describing the parallel machine.

blockMap:  (Out) An Epetra_BlockMap object constructed from file
contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system, returns 1
if number of processors differs from file creator. ";

%feature("docstring")  EpetraExt::MatrixMarketFileToRowMap "int
EpetraExt::MatrixMarketFileToRowMap(const char *filename, const
Epetra_Comm &comm, Epetra_BlockMap *&rowmap) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToBlockMaps "int
EpetraExt::MatrixMarketFileToBlockMaps(const char *filename, const
Epetra_Comm &comm, Epetra_BlockMap *&rowmap, Epetra_BlockMap *&colmap,
Epetra_BlockMap *&rangemap, Epetra_BlockMap *&domainmap)

Constructs row,col,range and domain maps from a matrix-market matrix
file. ";

%feature("docstring")  EpetraExt::BlockMapToMatrixMarketFile "int
EpetraExt::BlockMapToMatrixMarketFile(const char *filename, const
Epetra_BlockMap &blockMap, const char *mapName=0, const char
*mapDescription=0, bool writeHeader=true)

Writes an Epetra_BlockMap or Epetra_Map object to a Matrix Market
format file.

This function takes an Epetra_BlockMap or Epetra_Map object and writes
it to the specified file. The map can be distributed or serial. The
user can provide a strings containing the object name, a description,
and specify that header information should or should not be printed to
the file.

Special information is encoded in the comment field of this map that
allows for identical reproduction of the map, including distribution
across processors and element size information.

The first column of the output file will be the list of GIDs in the
map. If the block map has non-uniform sizes, a second column will be
generated containing the element sizes.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the map GIDs. A
second column may be present if the BlockMap has nonuniform sizes.

blockMap:  (In) An Epetra_BlockMap Object containing the user map to
be dumped to file.

mapName:  (In) A C-style string pointer to a name that will be stored
in the comment field of the file. This is not a required argument.
Note that it is possible to pass in the method A.Label().

mapDescription:  (In) A C-style string pointer to a map description
that will be stored in the comment field of the file.

writeHeader:  (In) If true, the header will be written, otherwise only
the map entries will be written.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::BlockMapToHandle "int
EpetraExt::BlockMapToHandle(FILE *handle, const Epetra_BlockMap &map)
";

%feature("docstring")  EpetraExt::writeBlockMap "int
EpetraExt::writeBlockMap(FILE *handle, int length, const int *v1,
const int *v2, bool doSizes) ";

%feature("docstring")  EpetraExt::BlockMapToHandle "int
EpetraExt::BlockMapToHandle(std::FILE *handle, const Epetra_BlockMap
&blockMap)

Writes an Epetra_BlockMap or Epetra_Map object to a file handle.

This function takes an Epetra_BlockMap or Epetra_Map object and writes
it to the specified file handle.

Parameters:
-----------

handle:  (In) A C-style file handle, already opened. On exit, the file
associated with this handle will have appended to it a row for each
multivector row.

blockMap:  (In) An Epetra_BlockMap object containing the user object
to be dumped to file.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::writeBlockMap "int
EpetraExt::writeBlockMap(std::FILE *handle, int length, const int *v1,
const int *v2, bool doSizes) ";

%feature("docstring")  EpetraExt::sort_three "static void
EpetraExt::sort_three(int *list, int *parlista, double *parlistb, int
start, int end) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Comm &comm, Epetra_CrsMatrix *&A) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Comm &comm, Epetra_CrsMatrix *&A, const bool transpose) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Comm &comm, Epetra_CrsMatrix *&A, const bool transpose=0, const
bool verbose=0)

Constructs an Epetra_CrsMatrix object from a Matrix Market format
file, simplest version: requires matrix to be square, distributes rows
evenly across processors.

This function constructs an Epetra_CrsMatrix object by reading a
Matrix Market file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matrix Market coordinate format.

comm:  (In) An Epetra_Comm object.

transpose:  (In) A boolean value indicating whether the reader should
transpose the matrix as it is read into matrix A. (default = 0).

verbose:  (In) A boolean value indicating whether the reader should
print diagnostic statements to stdout. (default = 0).

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system.   (See the
<a href=\"http://math.nist.gov/MatrixMarket\">Matrix Market</a> home
page for details.) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Map &rowMap, const Epetra_Map &rangeMap, const Epetra_Map
&domainMap, Epetra_CrsMatrix *&A, const bool transpose=0, const bool
verbose=0)

Constructs an Epetra_CrsMatrix object from a Matrix Market format
file, row, range and domain map specified; typically used for
rectangular matrices.

Reads an Epetra_CrsMatrix object from a matrix-market file, but uses
the specified maps for constructing and 'FillComplete()'ing the
matrix. Successfully creates rectangular matrices.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matrix Market coordinate format.

rowMap:  (In) An Epetra_Map object describing the desired row
distribution of the matrix.

rangeMap:  (In) An Epetra_Map object describing the distribution of
range vectors that will be used with this matrix, must be 1-to-1.

domainMap:  (In) An Epetra_Map object describing the distribution of
domain vectors that will be used with this matrix, must be 1-to-1.

transpose:  (In) A boolean value indicating whether the reader should
transpose the matrix as it is read into matrix A. (default = 0).

verbose:  (In) A boolean value indicating whether the reader should
print diagnostic statements to stdout. (default = 0).

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system.  (See the
<a href=\"http://math.nist.gov/MatrixMarket\">Matrix Market</a> home
page for details.) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Map &rowMap, Epetra_CrsMatrix *&A, const bool transpose=0,
const bool verbose=0)

Constructs an Epetra_CrsMatrix object from a Matrix Market format
file, only row map specified; allows user defined distribution of
matrix rows, requires square matrix.

This function constructs an Epetra_CrsMatrix object by reading a
Matrix Market file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matrix Market coordinate format.

rowMap:  (In) An Epetra_Map object describing the desired row
distribution of the matrix.

transpose:  (In) A boolean value indicating whether the reader should
transpose the matrix as it is read into matrix A. (default = 0).

verbose:  (In) A boolean value indicating whether the reader should
print diagnostic statements to stdout. (default = 0).

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system.   (See the
<a href=\"http://math.nist.gov/MatrixMarket\">Matrix Market</a> home
page for details.) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Map &rowMap, const Epetra_Map &colMap, Epetra_CrsMatrix *&A,
const bool transpose=0, const bool verbose=0)

Constructs an Epetra_CrsMatrix object from a Matrix Market format
file, both row and column map specified; this version is seldom used
unless you want explicit control over column map.

This function constructs an Epetra_CrsMatrix object by reading a
Matrix Market file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matrix Market coordinate format.

rowMap:  (In) An Epetra_Map object describing the desired row
distribution of the matrix.

colMap:  (In) An Epetra_Map object describing the desired column
distribution of the matrix.

transpose:  (In) A boolean value indicating whether the reader should
transpose the matrix as it is read into matrix A. (default = 0).

verbose:  (In) A boolean value indicating whether the reader should
print diagnostic statements to stdout. (default = 0).

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system.   (See the
<a href=\"http://math.nist.gov/MatrixMarket\">Matrix Market</a> home
page for details.) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrix "int
EpetraExt::MatrixMarketFileToCrsMatrix(const char *filename, const
Epetra_Map &rowMap, const Epetra_Map &colMap, const Epetra_Map
&rangeMap, const Epetra_Map &domainMap, Epetra_CrsMatrix *&A, const
bool transpose=0, const bool verbose=0)

Constructs an Epetra_CrsMatrix object from a Matrix Market format
file, row, column, range and domain map specified; this version is
seldom required unless you want explicit control over column map.

Reads an Epetra_CrsMatrix object from a matrix-market file, but uses
the specified maps for constructing and 'FillComplete()'ing the
matrix. Successfully creates rectangular matrices.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matrix Market coordinate format.

rowMap:  (In) An Epetra_Map object describing the desired row
distribution of the matrix.

colMap:  (In) An Epetra_Map object describing the desired column
distribution of the matrix.

rangeMap:  (In) An Epetra_Map object describing the distribution of
range vectors that will be used with this matrix, must be 1-to-1.

domainMap:  (In) An Epetra_Map object describing the distribution of
domain vectors that will be used with this matrix, must be 1-to-1.

transpose:  (In) A boolean value indicating whether the reader should
transpose the matrix as it is read into matrix A. (default = 0).

verbose:  (In) A boolean value indicating whether the reader should
print diagnostic statements to stdout. (default = 0).

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system.  (See the
<a href=\"http://math.nist.gov/MatrixMarket\">Matrix Market</a> home
page for details.) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToCrsMatrixHandle "int EpetraExt::MatrixMarketFileToCrsMatrixHandle(const char *filename,
const Epetra_Comm &comm, Epetra_CrsMatrix *&A, const Epetra_Map
*rowMap, const Epetra_Map *colMap, const Epetra_Map *rangeMap, const
Epetra_Map *domainMap, const bool transpose, const bool verbose) ";

%feature("docstring")  EpetraExt::quickpart_list_inc_int "static void
EpetraExt::quickpart_list_inc_int(int *list, int *parlista, double
*parlistb, int start, int end, int *equal, int *larger) ";

%feature("docstring")  EpetraExt::MatlabFileToCrsMatrix "int
EpetraExt::MatlabFileToCrsMatrix(const char *filename, const
Epetra_Comm &comm, Epetra_CrsMatrix *&A)

Constructs an Epetra_CrsMatrix object from a Matlab format file,
distributes rows evenly across processors.

This function constructs an Epetra_CrsMatrix object by reading a
Matlab (i,j,value) format file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matlab coordinate format.

comm:  (In) An Epetra_Comm object. The matrix will have its rows
distributed evenly by row-count across the parallel machine.

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.
The input matrix can be any dimension, square or rectangular.

WARNING:  User must delete matrix A!!.

Returns 0 if no error, -1 if any problems with file system.  Notes:
The file will be read twice: first to get the maximum row and column
dimensions. Next to insert values.

The global row and column dimensions will be determined by the maximum
row and column index, respectively, contained in the file. If some
rows or columns are empty they will still be present in the matrix.

The format expected for the input file is a list of nonzero entries
with one entry per row. Each row will have the row index, column index
and value listed with space in between each item. The number of lines
in the file should be exactly the number of entries of the matrix. For
example, consider the following matrix where only the nonzero values
are stored:

\\\\[ \\\\left[\\\\begin{array}{cccc} 5 & 7 & 0 & 0 \\\\\\\\ 3 & 2 & 0
& 1 \\\\\\\\ 0 & 0 & 0 & 4 \\\\\\\\ \\\\end{array}\\\\right]. \\\\]

A Matlab format file for this matrix would be: 1 1 5.0 1 2 7.0 2 1 3.0
2 2 2.0 2 4 1.0 4 4 4.0

Note that the entries can be listed in any order and that the matrix
does not need to be square. Values in the first and second columns
must be integer values and in those in the third column must be
floating point format.

(See the <a href=\"http://www.mathworks.com\">Matlab</a> home page for
details.) ";

%feature("docstring")  EpetraExt::HypreFileToCrsMatrix "int
EpetraExt::HypreFileToCrsMatrix(const char *filename, const
Epetra_Comm &comm, Epetra_CrsMatrix *&A)

Constructs an Epetra_CrsMatrix object from a Hypre Matrix Print
command, the row map is specified.

Reads an Epetra_CrsMatrix object from a Hypre Matrix Printout, the
matrix should be square.

Parameters:
-----------

filename:  (In) A filename not including the processor id extension,
including path if desired.

comm:  (In) An Epetra_Comm object describing the communication among
processors.

A:  (Out) An Epetra_CrsMatrix object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::mm_read_unsymmetric_sparse "int
EpetraExt::mm_read_unsymmetric_sparse(const char *fname, int *M_, int
*N_, int *nz_, double **val_, int **I_, int **J_) ";

%feature("docstring")  EpetraExt::mm_is_valid "int
EpetraExt::mm_is_valid(MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_read_banner "int
EpetraExt::mm_read_banner(FILE *f, MM_typecode *matcode) ";

%feature("docstring")  EpetraExt::mm_write_mtx_crd_size "int
EpetraExt::mm_write_mtx_crd_size(FILE *f, int M, int N, int nz) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd_size "int
EpetraExt::mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz) ";

%feature("docstring")  EpetraExt::mm_read_mtx_array_size "int
EpetraExt::mm_read_mtx_array_size(FILE *f, int *M, int *N) ";

%feature("docstring")  EpetraExt::mm_write_mtx_array_size "int
EpetraExt::mm_write_mtx_array_size(FILE *f, int M, int N) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd_data "int
EpetraExt::mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int
I[], int J[], double val[], MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd_entry "int
EpetraExt::mm_read_mtx_crd_entry(FILE *f, int *I, int *J, double
*real, double *imag, MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd "int
EpetraExt::mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int
**I, int **J, double **val, MM_typecode *matcode) ";

%feature("docstring")  EpetraExt::mm_write_banner "int
EpetraExt::mm_write_banner(FILE *f, MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_write_mtx_crd "int
EpetraExt::mm_write_mtx_crd(char fname[], int M, int N, int nz, int
I[], int J[], double val[], MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_typecode_to_str "void
EpetraExt::mm_typecode_to_str(MM_typecode matcode, char *buffer) ";

%feature("docstring")  EpetraExt::mm_read_banner "int
EpetraExt::mm_read_banner(std::FILE *f, MM_typecode *matcode) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd_size "int
EpetraExt::mm_read_mtx_crd_size(std::FILE *f, int *M, int *N, int *nz)
";

%feature("docstring")  EpetraExt::mm_read_mtx_array_size "int
EpetraExt::mm_read_mtx_array_size(std::FILE *f, int *M, int *N) ";

%feature("docstring")  EpetraExt::mm_write_banner "int
EpetraExt::mm_write_banner(std::FILE *f, MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_write_mtx_crd_size "int
EpetraExt::mm_write_mtx_crd_size(std::FILE *f, int M, int N, int nz)
";

%feature("docstring")  EpetraExt::mm_write_mtx_array_size "int
EpetraExt::mm_write_mtx_array_size(std::FILE *f, int M, int N) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd_data "int
EpetraExt::mm_read_mtx_crd_data(std::FILE *f, int M, int N, int nz,
int I[], int J[], double val[], MM_typecode matcode) ";

%feature("docstring")  EpetraExt::mm_read_mtx_crd_entry "int
EpetraExt::mm_read_mtx_crd_entry(std::FILE *f, int *I, int *J, double
*real, double *img, MM_typecode matcode) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToMultiVector "int
EpetraExt::MatrixMarketFileToMultiVector(const char *filename, const
Epetra_BlockMap &map, Epetra_MultiVector *&A)

Constructs an Epetra_MultiVector object from a Matrix Market format
file.

This function constructs an Epetra_MultiVector object by reading a
Matrix Market file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The multivector
to be read should be in this file in Matrix Market array format.

map:  (In) An Epetra_Map or Epetra_BlockMap object describing the
desired distribution of the multivector.

A:  (Out) An Epetra_MultiVector object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::MultiVectorToMatlabFile "int
EpetraExt::MultiVectorToMatlabFile(const char *filename, const
Epetra_MultiVector &A)

Writes an Epetra_MultiVector object to a file that is compatible with
Matlab.

This function takes any matrix that implements the Epetra_MultiVector
interface and writes it to the specified file. The matrix can be
distributed or serial. This function is a convenience wrapper around
MultiVectorToMatrixMarketFile. The following Matlab commands can be
used to read the resulting file and convert to it to a Matlab sparse
matrix: load filename;  For example: load A.dat;  The above produces a
dense matrix A with each vector in the multivector as a column in A.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contain a row for each row of the multivector.

A:  (In) An Epetra_MultiVector Object containing the user matrix to be
dumped to file.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::MultiVectorToMatrixMarketFile "int
EpetraExt::MultiVectorToMatrixMarketFile(const char *filename, const
Epetra_MultiVector &A, const char *matrixName=0, const char
*matrixDescription=0, bool writeHeader=true)

Writes an Epetra_MultiVector object to a Matrix Market format file.

This function takes an Epetra_MultiVector object and writes it to the
specified file. The multivector can be distributed or serial. The user
can provide a strings containing the object name, a description, and
specify that header information should or should not be printed to the
file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the matrix
coefficients. The file will contain a row for each entry. All entries
for a column are listed before going to the next column.

A:  (In) An Epetra_MultiVector Object containing the user matrix to be
dumped to file.

matrixName:  (In) A C-style string pointer to a name that will be
stored in the comment field of the file. This is not a required
argument. Note that it is possible to pass in the method A.Label().

matrixDescription:  (In) A C-style string pointer to a matrix
description that will be stored in the comment field of the file.

writeHeader:  (In) If true, the header will be written, otherwise only
the matrix entries will be written.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::MultiVectorToMatlabHandle "int
EpetraExt::MultiVectorToMatlabHandle(FILE *handle, const
Epetra_MultiVector &A) ";

%feature("docstring")  EpetraExt::MultiVectorToMatrixMarketHandle "int EpetraExt::MultiVectorToMatrixMarketHandle(FILE *handle, const
Epetra_MultiVector &A) ";

%feature("docstring")  EpetraExt::MultiVectorToHandle "int
EpetraExt::MultiVectorToHandle(FILE *handle, const Epetra_MultiVector
&A, bool mmFormat) ";

%feature("docstring")  EpetraExt::writeMultiVector "int
EpetraExt::writeMultiVector(FILE *handle, const Epetra_MultiVector &A,
bool mmFormat) ";

%feature("docstring")  EpetraExt::MultiVectorToMatrixMarketHandle "int EpetraExt::MultiVectorToMatrixMarketHandle(std::FILE *handle,
const Epetra_MultiVector &A)

Writes an Epetra_MultiVector object that is compatible with Matrix
Market array format to a file handle.

This function takes an Epetra_MultiVector and writes it to the
specified file handle.

Parameters:
-----------

handle:  (In) A C-style file handle, already opened. On exit, the file
associated with this handle will have appended to it a row for each
multivector row.

A:  (In) An Epetra_MultiVector Object containing the user object to be
dumped to file.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::MultiVectorToMatlabHandle "int
EpetraExt::MultiVectorToMatlabHandle(std::FILE *handle, const
Epetra_MultiVector &A)

Writes an Epetra_MultiVector object that is compatible with Matlab to
a file handle.

This function takes an Epetra_MultiVector and writes it to the
specified file handle.

Parameters:
-----------

handle:  (In) A C-style file handle, already opened. On exit, the file
associated with this handle will have appended to it a row for each
multivector row.

A:  (In) An Epetra_MultiVector Object containing the user object to be
dumped to file.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::MultiVectorToHandle "int
EpetraExt::MultiVectorToHandle(std::FILE *handle, const
Epetra_MultiVector &A, bool mmFormat) ";

%feature("docstring")  EpetraExt::writeMultiVector "int
EpetraExt::writeMultiVector(std::FILE *handle, const
Epetra_MultiVector &A, bool mmFormat) ";

%feature("docstring")  EpetraExt::OperatorToMatlabFile "int
EpetraExt::OperatorToMatlabFile(const char *filename, const
Epetra_Operator &A)

Writes an Epetra_Operator object to a file that is compatible with
Matlab.

This function takes any matrix that implements the Epetra_Operator
interface and writes it to the specified file. The matrix can be
distributed or serial. This function is a convenience wrapper around
OperatorToMatrixMarketFile. The following Matlab commands can be used
to read the resulting file and convert to it to a Matlab sparse
matrix: load filename;

matrix_name = spconvert(filename_root);  For example: load A.dat;

A = spconvert(filename_root);  The above produces a sparse matrix A.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contain a row for each matrix entry The first column is the global row
index, using base 1, the second column is the global column index of
the entry, the third value is the matrix value for that entry.

A:  (In) An Epetra_Operator Object containing the implicit user matrix
to be dumped to file. Any object that implements the Epetra_Operator
interface can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface, as is
AztecOO_Operator, all Ifpack and ML preconditioners.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::OperatorToMatrixMarketFile "int
EpetraExt::OperatorToMatrixMarketFile(const char *filename, const
Epetra_Operator &A, const char *matrixName=0, const char
*matrixDescription=0, bool writeHeader=true)

Writes an Epetra_Operator object to a Matrix Market format file,
forming the coefficients by applying the operator to the e_j vectors.

This function takes any linear operator that implements the
Epetra_Operator interface and writes it to the specified file. The
operator can be distributed or serial. The user can provide a strings
containing the matrix name, a matrix description, and specify that
header information should or should not be printed to the file.

The coeffients are formed by applying the operator to the canonical
vectors \\\\[ e_j = (0, \\\\ldots, 0, 1, 0, \\\\ldots, 0) \\\\] where
the value 1 appears in the the jth entry. The number of canonical
vectors used is determined by the size of the OperatorDomainMap() and
the lengths by the size of OperatorRangeMap().

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the matrix
coefficients. The file will contain a row for each matrix entry The
first column is the global row index, using base 1, the second column
is the global column index of the entry, the third value is the matrix
value for that entry.

A:  (In) An Epetra_Operator Object containing the implicit user matrix
to be dumped to file. Any object that implements the Epetra_Operator
interface can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface, as is
AztecOO_Operator, all Ifpack and ML preconditioners.

matrixName:  (In) A C-style string pointer to a name that will be
stored in the comment field of the file. This is not a required
argument. Note that it is possible to pass in the method A.Label().

matrixDescription:  (In) A C-style string pointer to a matrix
description that will be stored in the comment field of the file.

writeHeader:  (In) If true, the header will be written, otherwise only
the matrix entries will be written.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::OperatorToHandle "int
EpetraExt::OperatorToHandle(FILE *handle, const Epetra_Operator &A) ";

%feature("docstring")  EpetraExt::writeOperatorStrip "int
EpetraExt::writeOperatorStrip(FILE *handle, const Epetra_MultiVector
&y, const Epetra_Map &rootDomainMap, const Epetra_Map &rootRangeMap,
int startColumn) ";

%feature("docstring")  EpetraExt::get_nz "int EpetraExt::get_nz(const
Epetra_Operator &A, int &nz) ";

%feature("docstring")  EpetraExt::OperatorToHandle "int
EpetraExt::OperatorToHandle(std::FILE *handle, const Epetra_Operator
&A)

Writes an Epetra_Operator object to a format file that is compatible
with Matlab.

This function takes any matrix that implements the Epetra_Operator
interface and writes it to the specified file handle. The matrix can
be distributed or serial. This function is a convenience wrapper
around OperatorToMatrixMarketFile.

Parameters:
-----------

handle:  (In) A C-style file handle, already opened. On exit, the file
associated with this handle will have appended to it a row for each
matrix entry The first column is the global row index, using base 1,
the second column is the global column index of the entry, the third
value is the matrix value for that entry.

A:  (In) An Epetra_Operator Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_Operator
interface can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface, as is
AztecOO_Operator, all Ifpack and ML preconditioners.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::writeOperatorStrip "int
EpetraExt::writeOperatorStrip(std::FILE *handle, const
Epetra_MultiVector &y, const Epetra_Map &rootDomainMap, const
Epetra_Map &rootRangeMap, int startColumn) ";

%feature("docstring")  EpetraExt::readEpetraLinearSystem "void
EpetraExt::readEpetraLinearSystem(const std::string &fileName, const
Epetra_Comm &comm, Teuchos::RefCountPtr< Epetra_CrsMatrix > *A=NULL,
Teuchos::RefCountPtr< Epetra_Map > *map=NULL, Teuchos::RefCountPtr<
Epetra_Vector > *x=NULL, Teuchos::RefCountPtr< Epetra_Vector >
*b=NULL, Teuchos::RefCountPtr< Epetra_Vector > *xExact=NULL)

Read in an Epetra linear system from a file.

Parameters:
-----------

fileName:  [in] Name of the file to read in the linear system (see
file formats below).

comm:  [in] The communicator

map:  [out] The createe map. map==NULL is allowed on input in which
case this will not be returned.

A:  [out] The created matrix. A==NULL is allowed on input in which
case this will not be returned.

x:  [out] The created LHS vector. x==NULL is allowed on input in which
case this will not be returned.

b:  [out] The created RHS vector. b==NULL is allowed on input in which
case this will not be returned.

xExact:  [out] The created exact LHS vector (if known). xExact==NULL
is allowed on input in which case this will not be returned.

This function reads from a number file formats (*.triU, *.triS, *.mtx,
*.hb)

ToDo: Finish documentation!

ToDo: Put this in EpetraExt after the release is finished. ";

%feature("docstring")  EpetraExt::RowMatrixToMatlabFile "int
EpetraExt::RowMatrixToMatlabFile(const char *filename, const
Epetra_RowMatrix &A)

Writes an Epetra_RowMatrix object to a file that is compatible with
Matlab.

This function takes any matrix that implements the Epetra_RowMatrix
interface and writes it to the specified file. The matrix can be
distributed or serial. This function is a convenience wrapper around
RowMatrixToMatrixMarketFile. The following Matlab commands can be used
to read the resulting file and convert to it to a Matlab sparse
matrix: load filename;

matrix_name = spconvert(filename_root);  For example: load A.dat;

A = spconvert(filename_root);  The above produces a sparse matrix A.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contain a row for each matrix entry The first column is the global row
index, using base 1, the second column is the global column index of
the entry, the third value is the matrix value for that entry.

A:  (In) An Epetra_RowMatrix Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_RowMatrix
interface can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::RowMatrixToMatrixMarketFile "int
EpetraExt::RowMatrixToMatrixMarketFile(const char *filename, const
Epetra_RowMatrix &A, const char *matrixName=0, const char
*matrixDescription=0, bool writeHeader=true)

Writes an Epetra_RowMatrix object to a Matrix Market format file.

This function takes any matrix that implements the Epetra_RowMatrix
interface and writes it to the specified file. The matrix can be
distributed or serial. The user can provide a strings containing the
matrix name, a matrix description, and specify that header information
should or should not be printed to the file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the matrix
coefficients. The file will contain a row for each matrix entry The
first column is the global row index, using base 1, the second column
is the global column index of the entry, the third value is the matrix
value for that entry.

A:  (In) An Epetra_RowMatrix Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_RowMatrix
interface can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface.

matrixName:  (In) A C-style std::string pointer to a name that will be
stored in the comment field of the file. This is not a required
argument. Note that it is possible to pass in the method A.Label() if
the matrix is one of the four types: Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix.

matrixDescription:  (In) A C-style std::string pointer to a matrix
description that will be stored in the comment field of the file.

writeHeader:  (In) If true, the header will be written, otherwise only
the matrix entries will be written.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::RowMatrixToHandle "int
EpetraExt::RowMatrixToHandle(FILE *handle, const Epetra_RowMatrix &A)
";

%feature("docstring")  EpetraExt::writeRowMatrix "int
EpetraExt::writeRowMatrix(FILE *handle, const Epetra_RowMatrix &A) ";

%feature("docstring")  EpetraExt::RowMatrixToHandle "int
EpetraExt::RowMatrixToHandle(std::FILE *handle, const Epetra_RowMatrix
&A)

Writes an Epetra_RowMatrix object to a format file that is compatible
with Matlab.

This function takes any matrix that implements the Epetra_RowMatrix
interface and writes it to the specified file handle. The matrix can
be distributed or serial. This function is a convenience wrapper
around RowMatrixToMatrixMarketFile.

Parameters:
-----------

handle:  (In) A C-style file handle, already opened. On exit, the file
associated with this handle will have appended to it a row for each
matrix entry The first column is the global row index, using base 1,
the second column is the global column index of the entry, the third
value is the matrix value for that entry.

A:  (In) An Epetra_RowMatrix Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_RowMatrix
interface can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::writeRowMatrix "int
EpetraExt::writeRowMatrix(std::FILE *handle, const Epetra_RowMatrix
&A) ";

%feature("docstring")  EpetraExt::toString "std::string
EpetraExt::toString(const int &x) ";

%feature("docstring")  EpetraExt::toString "std::string
EpetraExt::toString(const unsigned int &x) ";

%feature("docstring")  EpetraExt::toString "std::string
EpetraExt::toString(const double &x) ";

%feature("docstring")  EpetraExt::MatrixMarketFileToVector "int
EpetraExt::MatrixMarketFileToVector(const char *filename, const
Epetra_BlockMap &map, Epetra_Vector *&A)

Constructs an Epetra_Vector object from a Matrix Market format file.

This function constructs an Epetra_Vector object by reading a Matrix
Market file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. The matrix to
be read should be in this file in Matrix Market array format.

map:  (In) An Epetra_Map or Epetra_BlockMap object describing the
desired distribution of the vector.

A:  (Out) An Epetra_Vector object constructed from file contents.

WARNING:  User must delete!!.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::VectorToMatlabFile "int
EpetraExt::VectorToMatlabFile(const char *filename, const
Epetra_Vector &A)

Writes an Epetra_Vector object to a file that is compatible with
Matlab.

This function takes any matrix that implements the Epetra_Vector
interface and writes it to the specified file. The matrix can be
distributed or serial. This function is a convenience wrapper around
VectorToMatrixMarketFile. The following Matlab commands can be used to
read the resulting file and convert to it to a Matlab sparse matrix:
load filename;

matrix_name = spconvert(filename_root);  For example: load A.dat;

A = spconvert(filename_root);  The above produces a sparse matrix A.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contain a row for each matrix entry The first column is the global row
index, using base 1, the second column is the global column index of
the entry, the third value is the matrix value for that entry.

A:  (In) An Epetra_Vector Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_Vector interface
can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::VectorToMatrixMarketFile "int
EpetraExt::VectorToMatrixMarketFile(const char *filename, const
Epetra_Vector &A, const char *matrixName=0, const char
*matrixDescription=0, bool writeHeader=true)

Writes an Epetra_Vector object to a Matrix Market format file.

This function takes any matrix that implements the Epetra_Vector
interface and writes it to the specified file. The matrix can be
distributed or serial. The user can provide a strings containing the
matrix name, a matrix description, and specify that header information
should or should not be printed to the file.

Parameters:
-----------

filename:  (In) A filename, including path if desired. If a file with
this name already exists, it will be deleted. On exit, this file will
contained any requested header information followed by the matrix
coefficients. The file will contain a row for each matrix entry The
first column is the global row index, using base 1, the second column
is the global column index of the entry, the third value is the matrix
value for that entry.

A:  (In) An Epetra_Vector Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_Vector interface
can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface.

matrixName:  (In) A C-style string pointer to a name that will be
stored in the comment field of the file. This is not a required
argument. Note that it is possible to pass in the method A.Label() if
the matrix is one of the four types: Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix.

matrixDescription:  (In) A C-style string pointer to a matrix
description that will be stored in the comment field of the file.

writeHeader:  (In) If true, the header will be written, otherwise only
the matrix entries will be written.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::VectorToHandle "int
EpetraExt::VectorToHandle(FILE *handle, const Epetra_Vector &A) ";

%feature("docstring")  EpetraExt::writeVector "int
EpetraExt::writeVector(FILE *handle, const Epetra_Vector &A) ";

%feature("docstring")  EpetraExt::VectorToHandle "int
EpetraExt::VectorToHandle(std::FILE *handle, const Epetra_Vector &A)

Writes an Epetra_Vector object to a format file that is compatible
with Matlab.

This function takes any matrix that implements the Epetra_Vector
interface and writes it to the specified file handle. The matrix can
be distributed or serial. This function is a convenience wrapper
around VectorToMatrixMarketFile.

Parameters:
-----------

handle:  (In) A C-style file handle, already opened. On exit, the file
associated with this handle will have appended to it a row for each
matrix entry The first column is the global row index, using base 1,
the second column is the global column index of the entry, the third
value is the matrix value for that entry.

A:  (In) An Epetra_Vector Object containing the user matrix to be
dumped to file. Any object that implements the Epetra_Vector interface
can be passed in. In particular, the Epetra_CrsMatrix,
Epetra_VbrMatrix, Epetra_FECrsMatrix, Epetra_FEVbrMatrix and
Epetra_MsrMatrix classes are compatible with this interface.

Returns 0 if no error, -1 if any problems with file system. ";

%feature("docstring")  EpetraExt::writeVector "int
EpetraExt::writeVector(std::FILE *handle, const Epetra_Vector &A) ";

%feature("docstring")  EpetraExt::compare_ints "int
EpetraExt::compare_ints(const void *a, const void *b)

Given an Epetra_CrsGraph that has block structure, an adjacency graph
is constructed representing the block connectivity of the original
graph.

David Day, Heidi Thornquist ";

%feature("docstring")  EpetraExt::ceil31log2 "int
EpetraExt::ceil31log2(int n) ";


// File: namespaceTeuchos.xml


// File: EpetraExt__AMD__CrsGraph_8cpp.xml


// File: EpetraExt__AMD__CrsGraph_8h.xml


// File: EpetraExt__BlockAdjacencyGraph_8cpp.xml


// File: EpetraExt__BlockAdjacencyGraph_8h.xml


// File: EpetraExt__BlockMapIn_8cpp.xml


// File: EpetraExt__BlockMapIn_8h.xml


// File: EpetraExt__BlockMapOut_8cpp.xml


// File: EpetraExt__BlockMapOut_8h.xml


// File: EpetraExt__ConfigDefs_8h.xml


// File: EpetraExt__CrsMatrixIn_8cpp.xml


// File: EpetraExt__CrsMatrixIn_8h.xml


// File: EpetraExt__CrsSingletonFilter__LinearProblem_8cpp.xml


// File: EpetraExt__CrsSingletonFilter__LinearProblem_8h.xml


// File: EpetraExt__Dirichlet__CrsMatrix_8cpp.xml


// File: EpetraExt__Dirichlet__CrsMatrix_8h.xml


// File: EpetraExt__Dirichlet__Vector_8cpp.xml


// File: EpetraExt__Dirichlet__Vector_8h.xml


// File: EpetraExt__DistArray_8h.xml


// File: EpetraExt__Exception_8h.xml


// File: EpetraExt__HDF5_8cpp.xml


// File: EpetraExt__HDF5_8h.xml


// File: EpetraExt__HDF5__DistObject_8cpp.xml


// File: EpetraExt__HDF5__Handle_8h.xml


// File: EpetraExt__LPTrans__From__GraphTrans_8cpp.xml


// File: EpetraExt__LPTrans__From__GraphTrans_8h.xml


// File: EpetraExt__LPTrans__From__MatrixTrans_8cpp.xml


// File: EpetraExt__LPTrans__From__MatrixTrans_8h.xml


// File: EpetraExt__MapColoring_8cpp.xml


// File: EpetraExt__MapColoring_8h.xml


// File: EpetraExt__MapColoringIndex_8cpp.xml


// File: EpetraExt__MapColoringIndex_8h.xml


// File: EpetraExt__MatrixMatrix_8cpp.xml


// File: EpetraExt__MatrixMatrix_8h.xml


// File: EpetraExt__MMHelpers_8cpp.xml


// File: EpetraExt__MMHelpers_8h.xml


// File: EpetraExt__mmio_8cpp.xml


// File: EpetraExt__mmio_8h.xml


// File: EpetraExt__MultiVectorIn_8cpp.xml


// File: EpetraExt__MultiVectorIn_8h.xml


// File: EpetraExt__MultiVectorOut_8cpp.xml


// File: EpetraExt__MultiVectorOut_8h.xml


// File: EpetraExt__OperatorOut_8cpp.xml


// File: EpetraExt__OperatorOut_8h.xml


// File: EpetraExt__Overlap__CrsGraph_8cpp.xml


// File: EpetraExt__Overlap__CrsGraph_8h.xml


// File: EpetraExt__Permutation_8cpp.xml


// File: EpetraExt__Permutation_8h.xml


// File: EpetraExt__Permutation__impl_8h.xml


// File: EpetraExt__ProductOperator_8cpp.xml


// File: EpetraExt__ProductOperator_8h.xml


// File: EpetraExt__readEpetraLinearSystem_8cpp.xml


// File: EpetraExt__readEpetraLinearSystem_8h.xml


// File: EpetraExt__Reindex__CrsMatrix_8cpp.xml


// File: EpetraExt__Reindex__CrsMatrix_8h.xml


// File: EpetraExt__Reindex__LinearProblem_8cpp.xml


// File: EpetraExt__Reindex__LinearProblem_8h.xml


// File: EpetraExt__Reindex__MultiVector_8cpp.xml


// File: EpetraExt__Reindex__MultiVector_8h.xml


// File: EpetraExt__RowMatrixOut_8cpp.xml


// File: EpetraExt__RowMatrixOut_8h.xml


// File: EpetraExt__Scale__LinearProblem_8cpp.xml


// File: EpetraExt__Scale__LinearProblem_8h.xml


// File: EpetraExt__SolverMap__CrsMatrix_8cpp.xml


// File: EpetraExt__SolverMap__CrsMatrix_8h.xml


// File: EpetraExt__SolverMap__LinearProblem_8cpp.xml


// File: EpetraExt__SolverMap__LinearProblem_8h.xml


// File: EpetraExt__StaticCondensation__LinearProblem_8cpp.xml


// File: EpetraExt__SubCopy__CrsMatrix_8cpp.xml


// File: EpetraExt__SubCopy__CrsMatrix_8h.xml


// File: EpetraExt__SymmRCM__CrsGraph_8cpp.xml


// File: EpetraExt__SymmRCM__CrsGraph_8h.xml


// File: EpetraExt__TimedEpetraOperator_8cpp.xml


// File: EpetraExt__Transform_8h.xml


// File: EpetraExt__Transform__Composite_8h.xml


// File: EpetraExt__Transpose__CrsGraph_8cpp.xml


// File: EpetraExt__Transpose__CrsGraph_8h.xml


// File: EpetraExt__Transpose__RowMatrix_8cpp.xml


// File: EpetraExt__Transpose__RowMatrix_8h.xml


// File: EpetraExt__Utils_8cpp.xml


// File: EpetraExt__Utils_8h.xml


// File: EpetraExt__VectorIn_8cpp.xml


// File: EpetraExt__VectorIn_8h.xml


// File: EpetraExt__VectorOut_8cpp.xml


// File: EpetraExt__VectorOut_8h.xml


// File: EpetraExt__Version_8h.xml


// File: EpetraExt__View__CrsGraph_8cpp.xml


// File: EpetraExt__View__CrsGraph_8h.xml


// File: EpetraExt__View__CrsMatrix_8cpp.xml


// File: EpetraExt__View__CrsMatrix_8h.xml


// File: EpetraExt__View__MultiVector_8cpp.xml


// File: EpetraExt__View__MultiVector_8h.xml


// File: EpetraExt__XMLReader_8cpp.xml
%feature("docstring")  Tokenize "static void Tokenize(const
std::string &str, std::vector< std::string > &tokens, const
std::string &delimiters=\" \") ";


// File: EpetraExt__XMLReader_8h.xml


// File: EpetraExt__XMLWriter_8cpp.xml


// File: EpetraExt__XMLWriter_8h.xml


// File: dir_11c281d317fac26a80734cfdf01c3276.xml


// File: dir_390ceb73d984add9550a91a5a05c10bc.xml


// File: dir_df40b4453a44ce93def4751f410362de.xml


// File: dir_1446c85883e8456091748f28f5dc619d.xml


// File: dir_7af1f200b770e582760f7ec5a7e1d7c3.xml


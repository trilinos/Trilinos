
// File: index.xml

// File: classGaleri_1_1FiniteElements_1_1AbstractGrid.xml
%feature("docstring") Galeri::FiniteElements::AbstractGrid "

Abstract interface to access finite element grids.

AbstractGrid is a pure virtual function, that specifies the interface
methods that a grid class must implement. Following an approach
similar to getrow() for matrices, there is no grid format; instead,
the user must implement a set of getelement() and similar methods.
Therefore, it is possible to define grids that are built on-the-fly, a
feature that is particularly convenient in testing phase.

This format is based on the following assumptions: All grid elements
are of the same type (for examples, all triangles). It is not possible
to have mixed grids (i.e., with some triangles and some
quadrilateral).

All grid elements are 3D. For 2D problems, the user must specify a
z-coordinate, for example 0.0.

Elements, vertices and faces are numbered locally. The local-to-global
mapping is required for vertices only, and must be defined using
Epetra_Map's.

Two Epetra_Map's must be created: The VertexMap() locally contains all
the vertices that belong to the local finite elements. A global vertex
can be replicated over more than one process.

The RowMap() globally contains all the global vertices. A global
vertex is assigned to exactly one process.

We require two maps for the following reason: In parallel runs,
vertices corresponding to internal boundary faces can are replicated
over processors. This makes the contruction of the finite element
matrix easier, since it is possible to work on local quantities only.
However, such a distribution is not compatible with AztecOO and ML
(and several other Trilinos packages), since these libraries require
each row to belong to exactly one processor. Methods
ExportToVertexMap() and ExportToRowMap() export objects from one map
to the other. Typically, ExportToRowMap() is used in the assembly
phase, whese local stiffness matrix and right-hand side are built
using VertexMap(), then exported to RowMap(). ExportToVertexMap(),
instead, is used after the solution of the linear system is computed,
to update the values of the solution on the local vertices, so that
norms can be computed, and the solution visualized. These methods are
trivial in the serial case.

Marzio Sala, SNL 9214

C++ includes: Galeri_AbstractGrid.h ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::~AbstractGrid "virtual
Galeri::FiniteElements::AbstractGrid::~AbstractGrid()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumDimensions "virtual int
Galeri::FiniteElements::AbstractGrid::NumDimensions() const =0

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumVerticesPerElement "virtual
int Galeri::FiniteElements::AbstractGrid::NumVerticesPerElement()
const =0

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::AbstractGrid::NumFacesPerElement() const =0

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::AbstractGrid::NumVerticesPerFace() const =0

Returns the number of vertices contained in each face. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ElementType "virtual string
Galeri::FiniteElements::AbstractGrid::ElementType() const =0

Returns a string containing the element type.

Returns a string containing the type of element. This string is used
in the quadrature class. Currently supported options are:
\"GALERI_TRIANGLE\"

\"GALERI_QUAD\"

\"GALERI_HEX\"

\"GALERI_TET\" ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumNeighborsPerElement "virtual
int Galeri::FiniteElements::AbstractGrid::NumNeighborsPerElement()
const =0

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumMyElements "virtual int
Galeri::FiniteElements::AbstractGrid::NumMyElements() const =0

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumGlobalElements "virtual int
Galeri::FiniteElements::AbstractGrid::NumGlobalElements() const =0

Returns the global number of finite elements. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumMyVertices "virtual int
Galeri::FiniteElements::AbstractGrid::NumMyVertices() const =0

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumGlobalVertices "virtual int
Galeri::FiniteElements::AbstractGrid::NumGlobalVertices() const =0

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::AbstractGrid::NumMyBoundaryFaces() const =0

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::NumGlobalBoundaryFaces "virtual
int Galeri::FiniteElements::AbstractGrid::NumGlobalBoundaryFaces()
const =0

Returns the global number of boundary faces. ";

%feature("docstring")  Galeri::FiniteElements::AbstractGrid::MyVolume
"virtual double Galeri::FiniteElements::AbstractGrid::MyVolume()
const =0

Returns the volume of all local elements. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::GlobalVolume "virtual double
Galeri::FiniteElements::AbstractGrid::GlobalVolume() const =0

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::VertexCoord "virtual void
Galeri::FiniteElements::AbstractGrid::VertexCoord(const int
LocalVertex, double *coord) const =0

Returns the coordinates of local vertex LocalVertex in vector coord.

Parameters:
-----------

LocalVertex:  - (In) Local ID of the vertex for whic coordinates are
required. Must be contained in the interval [0, NumMyVertices())

coord:  - (Out) double array of size 3. In output, contains the x-, y-
and z-coordinate of the specified vertex.

Parameter coord must be allocated of size 3 for both 2D and 3D
problems. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::VertexCoord "virtual void
Galeri::FiniteElements::AbstractGrid::VertexCoord(const int Length,
const int *IDs, double *x, double *y, double *z) const =0

Returns the coordinates of specified local vertices.

Parameters:
-----------

Length:  - (In) Length of array IDs.

IDs:  - (In) Contains the list of vertices of which coordinates are
required.

x:  - (Out) double array of size Length. In output, contains the
x-coordinates of the specified vertices.

y:  - (Out) double array of size Length. In output, contains the
y-coordinates of the specified vertices.

z:  - (Out) double array of size Length. In output, contains the
z-coordinates of the specified vertices.

The z array must be allocated for both 2D and 3D problems. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ElementVertices "virtual void
Galeri::FiniteElements::AbstractGrid::ElementVertices(const int
LocalElement, int *elements) const =0

Returns the local vertex IDs of the specified local finite element.

Parameters:
-----------

LocalElement:  - (In) ID of the required local element.

elements:  - (Out) array of length NumElementVertices(), in output
will contain the local ID of the vertices of the specified element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ElementNeighbors "virtual void
Galeri::FiniteElements::AbstractGrid::ElementNeighbors(const int
LocalElement, int *elements) const =0

Returns the local IDs of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::FaceVertices "virtual void
Galeri::FiniteElements::AbstractGrid::FaceVertices(const int
LocalFace, int &tag, int *IDs) const =0

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")  Galeri::FiniteElements::AbstractGrid::FacePatch
"virtual int Galeri::FiniteElements::AbstractGrid::FacePatch(const
int LocalFace) const =0

Returns the patch ID of the specified face.

Returns an integer ID that identifies the given boundary face as
belonging to a given part of the domain. It can be used by the user to
specify the value and the type of the boundary condition. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ElementMinLength "virtual
double Galeri::FiniteElements::AbstractGrid::ElementMinLength(const
int LocalElement) const =0

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ElementMaxLength "virtual
double Galeri::FiniteElements::AbstractGrid::ElementMaxLength(const
int LocalElement) const =0

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ElementVolume "virtual double
Galeri::FiniteElements::AbstractGrid::ElementVolume(const int
LocalElement) const =0

Returns the volume of the specified local finite element.

Returns the area (in 2D) or the volume (in 3D) of the specified local
element ";

%feature("docstring")  Galeri::FiniteElements::AbstractGrid::FaceArea
"virtual double Galeri::FiniteElements::AbstractGrid::FaceArea(const
int LocalFace) const =0

Returns the area of the specified local face.

Returns the length (in 2D) or the area (in 3D) of the specified
boundary face ";

%feature("docstring")  Galeri::FiniteElements::AbstractGrid::VertexMap
"virtual const Epetra_Map&
Galeri::FiniteElements::AbstractGrid::VertexMap() const =0

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")  Galeri::FiniteElements::AbstractGrid::RowMap "virtual const Epetra_Map&
Galeri::FiniteElements::AbstractGrid::RowMap() const =0

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ExportToVertexMap "virtual void
Galeri::FiniteElements::AbstractGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const
=0

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::AbstractGrid::ExportToRowMap "virtual void
Galeri::FiniteElements::AbstractGrid::ExportToRowMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const
=0

Exports distributed object from VertexMap() to RowMap(). ";

%feature("docstring")  Galeri::FiniteElements::AbstractGrid::Comm "virtual const Epetra_Comm&
Galeri::FiniteElements::AbstractGrid::Comm() const =0

Returns a reference to the communicator object. ";


// File: classGaleri_1_1FiniteElements_1_1AbstractProblem.xml
%feature("docstring") Galeri::FiniteElements::AbstractProblem "

Abstract interface to define linear problems.

AbstractProblem defines a set of abstract interfaces, used to
construct the linear system corresponding to the finite element
discretization of a scalar PDE problem. A concrete implementation will
require an AbstractGrid and an AbstractVariational object; the former
is used to query for the grid elements, the latter to integrate the
variational form over such elements. The role of AbstractProblem is to
take the elemental matrices, given by AbstractVariational, and insert
them into the global, distributed matrix (whose RowMatrixRowMap() is
given by Grid().RowMap()).

Marzio Sala, SNL 9214.

C++ includes: Galeri_AbstractProblem.h ";

%feature("docstring")
Galeri::FiniteElements::AbstractProblem::~AbstractProblem "virtual
Galeri::FiniteElements::AbstractProblem::~AbstractProblem()

Destructor. ";

%feature("docstring")  Galeri::FiniteElements::AbstractProblem::A "virtual Epetra_RowMatrix&
Galeri::FiniteElements::AbstractProblem::A()=0

Returns a reference to the linear system matrix. ";

%feature("docstring")  Galeri::FiniteElements::AbstractProblem::RHS "virtual Epetra_MultiVector&
Galeri::FiniteElements::AbstractProblem::RHS()=0

Returns a reference to the multi-vector of right-hand side. ";

%feature("docstring")  Galeri::FiniteElements::AbstractProblem::LHS "virtual Epetra_MultiVector&
Galeri::FiniteElements::AbstractProblem::LHS()=0

Returns a reference to the multi-vector of starting solution. ";

%feature("docstring")  Galeri::FiniteElements::AbstractProblem::Grid "virtual const AbstractGrid&
Galeri::FiniteElements::AbstractProblem::Grid() const =0

Returns a reference to the grid object. ";

%feature("docstring")
Galeri::FiniteElements::AbstractProblem::Variational "virtual const
AbstractVariational&
Galeri::FiniteElements::AbstractProblem::Variational() const =0

Returns a reference to the variational object. ";

%feature("docstring")
Galeri::FiniteElements::AbstractProblem::Compute "virtual void
Galeri::FiniteElements::AbstractProblem::Compute()=0

Computes the linear system matrix, LHS and RHS. ";

%feature("docstring")
Galeri::FiniteElements::AbstractProblem::ComputeNorms "virtual void
Galeri::FiniteElements::AbstractProblem::ComputeNorms(Epetra_MultiVector
&RowMatrixField, int(*ExactSolution)(double, double, double, double
*), const bool verbose=true, double *SolutionNorm=0, double
*ExactNorm=0, double *DiffNorm=0)=0

Computes the norm of computed solution, exact solution, and error.

Parameters:
-----------

RowMatrixField:  - (In) Multi-vector defined on Grid().RowMap() which
contains the numerical solution.

ExactSolution:  - (In) Function defined as in the following example:
ExactSolution(double x, double y, double x, double* sol) will contain
the value of the solution in sol[0], the x-derivative in sol[1], the
y-derivative in sol[2], and the z-derivative in sol[2].

verbose:  - (In) If true, prints out the results.

SolutionNorm:  - (Out) a double array of size 3, which will contain
the L2 norm, the semi-H1 norm and the H1-norm of the numerical
solution.

ExactNorm:  - (Out) a double array of size 3, which will contain the
L2 norm, the semi-H1 norm and the H1-norm of the exact solution.

DiffNorm:  - (Out) a double array of size 3, which will contain the L2
norm, the semi-H1 norm and the H1-norm of the error. ";


// File: classGaleri_1_1FiniteElements_1_1AbstractQuadrature.xml
%feature("docstring") Galeri::FiniteElements::AbstractQuadrature "

Interfaces for quadrature over elements.

AbstractQuadrature is a pure virtual class that defines a set of
abstract interfaces to basis and test functions (and their
derivatives), and also furnishes all the tools required to numerically
integrate over an element.

Marzio Sala, SNL 9214.

C++ includes: Galeri_AbstractQuadrature.h ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::~AbstractQuadrature "virtual
Galeri::FiniteElements::AbstractQuadrature::~AbstractQuadrature()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::NumQuadrNodes "virtual
int Galeri::FiniteElements::AbstractQuadrature::NumQuadrNodes() const
=0

Returns the number of quadrature node per element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::NumPhiFunctions "virtual
int Galeri::FiniteElements::AbstractQuadrature::NumPhiFunctions()
const =0

Returns the number of basis function on the reference element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::NumPsiFunctions "virtual
int Galeri::FiniteElements::AbstractQuadrature::NumPsiFunctions()
const =0

Returns the number of test function on the reference element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::ComputeJacobian "virtual
void Galeri::FiniteElements::AbstractQuadrature::ComputeJacobian(const
int QuadrNode, const double *x, const double *y, const double *z)
const =0

Computes the Jacobian at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::ComputeQuadrNodes "virtual void
Galeri::FiniteElements::AbstractQuadrature::ComputeQuadrNodes(const
int QuadrNode, const double *x, const double *y, const double *z,
double &xq, double &yq, double &zq) const =0

Maps the quadrature nodes from the reference element to the actual
one. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::ComputeDerivatives "virtual void
Galeri::FiniteElements::AbstractQuadrature::ComputeDerivatives(const
int QuadrNode) const =0

Computes the derivatives at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::QuadrWeight "virtual
double Galeri::FiniteElements::AbstractQuadrature::QuadrWeight(const
int QuadrNode) const =0

Computes the weight at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::DetJacobian "virtual
double Galeri::FiniteElements::AbstractQuadrature::DetJacobian(const
int QuadrNode) const =0

Computes the determinant of the Jacobian matrix at the quadrature
node. ";

%feature("docstring")  Galeri::FiniteElements::AbstractQuadrature::Phi
"virtual double Galeri::FiniteElements::AbstractQuadrature::Phi(const
int i) const =0

Returns the value of the i-th basis function on the reference element.
";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::PhiX "virtual double
Galeri::FiniteElements::AbstractQuadrature::PhiX(const int i) const =0

Returns the value of the x-derivative i-th basis function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::PhiY "virtual double
Galeri::FiniteElements::AbstractQuadrature::PhiY(const int i) const =0

Returns the value of the y-derivative i-th basis function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::PhiZ "virtual double
Galeri::FiniteElements::AbstractQuadrature::PhiZ(const int i) const =0

Returns the value of the z-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::AbstractQuadrature::Psi
"virtual double Galeri::FiniteElements::AbstractQuadrature::Psi(const
int i) const =0

Returns the value of the i-th test function on the reference element.
";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::PsiX "virtual double
Galeri::FiniteElements::AbstractQuadrature::PsiX(const int i) const =0

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::PsiY "virtual double
Galeri::FiniteElements::AbstractQuadrature::PsiY(const int i) const =0

Returns the value of the y-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractQuadrature::PsiZ "virtual double
Galeri::FiniteElements::AbstractQuadrature::PsiZ(const int i) const =0

Returns the value of the z-derivative i-th test function on the
reference element. ";


// File: classGaleri_1_1FiniteElements_1_1AbstractVariational.xml
%feature("docstring") Galeri::FiniteElements::AbstractVariational "

Pure virtual class that defines the variational form.

AbstractVariational is a pure virtual class, that specifies a set of
abstract interfaces, required to integrate the variational form and
the right-hand side over an element, and to compute the norm of
compute solution, exact solution, and error over the element. A
concrete implementation of this class also defined how boundary
conditions are resolved.

The element on which the integration is performed is specified by
providing the coordinates of the local vertices.

Marzio Sala, SNL 9214.

C++ includes: Galeri_AbstractVariational.h ";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::~AbstractVariational "virtual
Galeri::FiniteElements::AbstractVariational::~AbstractVariational()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::LHS "virtual double
Galeri::FiniteElements::AbstractVariational::LHS(const double Phi,
const double Psi, const double PhiX, const double PsiX, const double
PhiY, const double PsiY, const double PhiZ, const double PsiZ, const
double x, const double y, const double z) const =0

Evaluates the bilinear form (without integral) at point (x,y,z). ";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::RHS "virtual double
Galeri::FiniteElements::AbstractVariational::RHS(const double Psi,
const double PsiX, const double PsiY, const double PsiZ, const double
x, const double y, const double z) const =0

Returns the value of the right-hand side (without integral) at point
(x, y, z). ";

%feature("docstring")  Galeri::FiniteElements::AbstractVariational::BC
"virtual int Galeri::FiniteElements::AbstractVariational::BC(const
int PatchID) const =0

Returns an integer identifying the boundary condition assigned to the
specified patch. ";

%feature("docstring")  Galeri::FiniteElements::AbstractVariational::BC
"virtual double Galeri::FiniteElements::AbstractVariational::BC(const
double x, const double y, const double z, const int PatchID) const =0

Returns the value of the boundary condition at point (x, y, z). ";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::IntegrateOverElement "virtual int
Galeri::FiniteElements::AbstractVariational::IntegrateOverElement(const
AbstractVariational &Variational, const double *x, const double *y,
const double *z, const double *data, double *ElementMatrix, double
*ElementRHS) const =0

Integrates the bilinear form and the right-hand side over the element.
";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::ElementNorm "virtual int
Galeri::FiniteElements::AbstractVariational::ElementNorm(const double
*LocalSol, const double *x, const double *y, const double *z, double
*Norm) const =0

Computes the norm of the computed solution over the element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::ElementNorm "virtual int
Galeri::FiniteElements::AbstractVariational::ElementNorm(int(*ExactSolution)(double,
double, double, double *), const double *x, const double *y, const
double *z, double *Norm) const =0

Computed the norm of the exact solution over the element. ";

%feature("docstring")
Galeri::FiniteElements::AbstractVariational::ElementNorm "virtual int
Galeri::FiniteElements::AbstractVariational::ElementNorm(const double
*LocalSol, int(*ExactSolution)(double, double, double, double *),
const double *x, const double *y, const double *z, double *Norm) const
=0

Computed the norm of the computed and exact solution over the element.
";


// File: classGaleri_1_1Exception.xml
%feature("docstring") Galeri::Exception "";

%feature("docstring")  Galeri::Exception::Exception "Galeri::Exception::Exception(const string FileName, const int
LineNumber, const string Line1, const string Line2=\"\", const string
Line3=\"\", const string Line4=\"\", const string Line5=\"\", const
string Line6=\"\") ";

%feature("docstring")  Galeri::Exception::Print "void
Galeri::Exception::Print() ";


// File: classGaleri_1_1FiniteElements_1_1FileGrid.xml
%feature("docstring") Galeri::FiniteElements::FileGrid "";

%feature("docstring")  Galeri::FiniteElements::FileGrid::FileGrid "Galeri::FiniteElements::FileGrid::FileGrid(const Epetra_Comm &Comm,
const string FileName)

Constructor.

Parameters:
-----------

Comm:  - (In) Communicator object.

FileName:  - (In) Name of grid file. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::~FileGrid "virtual Galeri::FiniteElements::FileGrid::~FileGrid() ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::NumDimensions
"virtual int Galeri::FiniteElements::FileGrid::NumDimensions() const

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumVerticesPerElement "virtual int
Galeri::FiniteElements::FileGrid::NumVerticesPerElement() const

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::FileGrid::NumFacesPerElement() const

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::FileGrid::NumVerticesPerFace() const

Returns the number of vertices contained in each face. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::ElementType "virtual string Galeri::FiniteElements::FileGrid::ElementType() const

Returns a string containing the element type.

Returns a string containing the type of element. This string is used
in the quadrature class. Currently supported options are:
\"GALERI_TRIANGLE\"

\"GALERI_QUAD\"

\"GALERI_HEX\"

\"GALERI_TET\" ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::Comm "virtual const Epetra_Comm& Galeri::FiniteElements::FileGrid::Comm()
const

Returns a reference to the communicator object. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::NumMyElements
"virtual int Galeri::FiniteElements::FileGrid::NumMyElements() const

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumGlobalElements "virtual int
Galeri::FiniteElements::FileGrid::NumGlobalElements() const

Returns the global number of finite elements. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::NumMyVertices
"virtual int Galeri::FiniteElements::FileGrid::NumMyVertices() const

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumGlobalVertices "virtual int
Galeri::FiniteElements::FileGrid::NumGlobalVertices() const

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::FileGrid::NumMyBoundaryFaces() const

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumGlobalBoundaryFaces "virtual int
Galeri::FiniteElements::FileGrid::NumGlobalBoundaryFaces() const

Returns the global number of boundary faces. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::VertexCoord "virtual void Galeri::FiniteElements::FileGrid::VertexCoord(const int
LocalID, double *coord) const

Returns the coordinates of local vertex LocalVertex in vector coord.

Parameters:
-----------

LocalVertex:  - (In) Local ID of the vertex for whic coordinates are
required. Must be contained in the interval [0, NumMyVertices())

coord:  - (Out) double array of size 3. In output, contains the x-, y-
and z-coordinate of the specified vertex.

Parameter coord must be allocated of size 3 for both 2D and 3D
problems. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::VertexCoord "virtual void Galeri::FiniteElements::FileGrid::VertexCoord(const int
Length, const int *IDs, double *x, double *y, double *z) const

Returns the coordinates of specified local vertices.

Parameters:
-----------

Length:  - (In) Length of array IDs.

IDs:  - (In) Contains the list of vertices of which coordinates are
required.

x:  - (Out) double array of size Length. In output, contains the
x-coordinates of the specified vertices.

y:  - (Out) double array of size Length. In output, contains the
y-coordinates of the specified vertices.

z:  - (Out) double array of size Length. In output, contains the
z-coordinates of the specified vertices.

The z array must be allocated for both 2D and 3D problems. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::ElementVertices "virtual void
Galeri::FiniteElements::FileGrid::ElementVertices(const int LocalID,
int *elements) const

Returns the local vertex IDs of the specified local finite element.

Parameters:
-----------

LocalElement:  - (In) ID of the required local element.

elements:  - (Out) array of length NumElementVertices(), in output
will contain the local ID of the vertices of the specified element. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::ElementMinLength "virtual double
Galeri::FiniteElements::FileGrid::ElementMinLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::ElementMaxLength "virtual double
Galeri::FiniteElements::FileGrid::ElementMaxLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::RCPVertexMap
"virtual const RefCountPtr<Epetra_Map>
Galeri::FiniteElements::FileGrid::RCPVertexMap() const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::RCPElementMap
"virtual const RefCountPtr<Epetra_Map>
Galeri::FiniteElements::FileGrid::RCPElementMap() const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::VertexMap "virtual const Epetra_Map&
Galeri::FiniteElements::FileGrid::VertexMap() const

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")  Galeri::FiniteElements::FileGrid::ElementMap "virtual const Epetra_Map&
Galeri::FiniteElements::FileGrid::ElementMap() const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::RowMap "virtual const Epetra_Map& Galeri::FiniteElements::FileGrid::RowMap()
const

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")  Galeri::FiniteElements::FileGrid::Importer "virtual const Epetra_Import&
Galeri::FiniteElements::FileGrid::Importer() const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::ElementTag "virtual int Galeri::FiniteElements::FileGrid::ElementTag(const int
LocalID) const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::VertexTag "virtual int Galeri::FiniteElements::FileGrid::VertexTag(const int
LocalID) const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::ElementVolume
"virtual double Galeri::FiniteElements::FileGrid::ElementVolume()
const ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::FaceVertices
"virtual void Galeri::FiniteElements::FileGrid::FaceVertices(const
int LocalFace, int &tag, int *IDs) const

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::FacePatch "int Galeri::FiniteElements::FileGrid::FacePatch(const int LocalFace)
const

Returns the patch ID of the specified face.

Returns an integer ID that identifies the given boundary face as
belonging to a given part of the domain. It can be used by the user to
specify the value and the type of the boundary condition. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::ElementVolume
"virtual double Galeri::FiniteElements::FileGrid::ElementVolume(const
int LocalElement) const

Returns the volume of the specified local finite element.

Returns the area (in 2D) or the volume (in 3D) of the specified local
element ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::FaceArea "virtual double Galeri::FiniteElements::FileGrid::FaceArea(const int
LocalFace) const

Returns the area of the specified local face.

Returns the length (in 2D) or the area (in 3D) of the specified
boundary face ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::MyVolume "virtual double Galeri::FiniteElements::FileGrid::MyVolume() const

Returns the volume of all local elements. ";

%feature("docstring")  Galeri::FiniteElements::FileGrid::GlobalVolume
"virtual double Galeri::FiniteElements::FileGrid::GlobalVolume()
const

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::ExportToVertexMap "void
Galeri::FiniteElements::FileGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::ExportToRowMap "void
Galeri::FiniteElements::FileGrid::ExportToRowMap(const
Epetra_DistObject &VertexObject, Epetra_DistObject &RowObject) const

Exports distributed object from VertexMap() to RowMap(). ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::NumNeighborsPerElement "int
Galeri::FiniteElements::FileGrid::NumNeighborsPerElement() const

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::FileGrid::ElementNeighbors "void
Galeri::FiniteElements::FileGrid::ElementNeighbors(int, int *) const

Returns the local IDs of neighboring elements. ";


// File: classGaleri__FileGrid.xml
%feature("docstring") Galeri_FileGrid "

Reads a grid from file, in Galeri format.

Marzio Sala, ETHZ/COLAB.

C++ includes: Galeri_FileGrid.h ";


// File: classGaleri_1_1FiniteElements_1_1GalerkinVariational.xml
%feature("docstring") Galeri::FiniteElements::GalerkinVariational "

Defines a pure Galerkin variational form of a scalar PDE.

This class defines a pure Galerkin variational form of a second order,
symmetric scalar PDE, discretized using Lagrange finite elements. The
class is templated with an AbstractQuadrature class, which will be
used to specify the quadrature formula, and the values of test and
basis functions at the quadrature node. The constructor requires
function pointers, that specify the values of the coefficients.

Marzio Sala, SNL 9214.

C++ includes: Galeri_GalerkinVariational.h ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::GalerkinVariational "Galeri::FiniteElements::GalerkinVariational< T
>::GalerkinVariational(const int NumQuadratureNodes,
double(*diff)(const double &, const double &, const double &),
double(*source)(const double &, const double &, const double &),
double(*force)(const double &, const double &, const double &),
double(*bc)(const double &, const double &, const double &, const int
&), int(*bc_type)(const int &))

Constructor. ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::~GalerkinVariational "Galeri::FiniteElements::GalerkinVariational< T
>::~GalerkinVariational()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::diff "double
Galeri::FiniteElements::GalerkinVariational< T >::diff(const double x,
const double y, const double z) const

Evaluates the diffusion coefficient at point (x, y, z). ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::source "double
Galeri::FiniteElements::GalerkinVariational< T >::source(const double
x, const double y, const double z) const

Evaluates the source term at point (x, y, z). ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::force "double
Galeri::FiniteElements::GalerkinVariational< T >::force(const double
x, const double y, const double z) const

Evaluates the force term at point (x, y, z). ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::IntegrateOverElement "virtual int Galeri::FiniteElements::GalerkinVariational< T
>::IntegrateOverElement(const AbstractVariational &Variational, const
double *x, const double *y, const double *z, const double *data,
double *ElementMatrix, double *ElementRHS) const

Integrates the variational form and the right-hand side. ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::ElementNorm "virtual int
Galeri::FiniteElements::GalerkinVariational< T >::ElementNorm(const
double *LocalSol, const double *x, const double *y, const double *z,
double *Norm) const

Computes the norm of the numerical solution over an element. ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::ElementNorm "virtual int
Galeri::FiniteElements::GalerkinVariational< T
>::ElementNorm(int(*ExactSolution)(double, double, double, double *),
const double *x, const double *y, const double *z, double *Norm) const

Computes the norm of the exact solution over an element. ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::ElementNorm "virtual int
Galeri::FiniteElements::GalerkinVariational< T >::ElementNorm(const
double *LocalSol, int(*ExactSolution)(double, double, double, double
*), const double *x, const double *y, const double *z, double *Norm)
const

Computes the norm of the error over an element. ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::LHS "double
Galeri::FiniteElements::GalerkinVariational< T >::LHS(const double
Phi, const double Psi, const double PhiX, const double PsiX, const
double PhiY, const double PsiY, const double PhiZ, const double PsiZ,
const double x, const double y, const double z) const

Evaluates the left-hand side at point (x, y, z). ";

%feature("docstring")
Galeri::FiniteElements::GalerkinVariational::RHS "double
Galeri::FiniteElements::GalerkinVariational< T >::RHS(const double
Psi, const double PsiX, const double PsiY, const double PsiZ, const
double x, const double y, const double z) const

Evaluates the right-hand side at point (x, y, z). ";

%feature("docstring")  Galeri::FiniteElements::GalerkinVariational::BC
"int Galeri::FiniteElements::GalerkinVariational< T >::BC(const int
PatchID) const

Returns the boundary condition type of the specified patch. ";

%feature("docstring")  Galeri::FiniteElements::GalerkinVariational::BC
"double Galeri::FiniteElements::GalerkinVariational< T >::BC(const
double x, const double y, const double z, const int PatchID) const

Returns the value of the boundary condition at point (x, y, z). ";


// File: classGaleri_1_1FiniteElements_1_1HexCubeGrid.xml
%feature("docstring") Galeri::FiniteElements::HexCubeGrid "

Creates a grid composed by hexahedra in a cube.

Marzio Sala, SNL 9214.

C++ includes: Galeri_HexCubeGrid.h ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::HexCubeGrid "Galeri::FiniteElements::HexCubeGrid::HexCubeGrid(Epetra_Comm &Comm,
const int nx, const int ny, const int nz, const int mx, const int my,
const int mz, const double lx=1.0, const double ly=1.0, const double
lz=1.0)

Constructor.

Parameters:
-----------

Comm:  - (In) Communicator object.

nx:  - (In) number of elements along the X-axis.

ny:  - (In) number of elements along the Y-axis.

mx:  - (In) Number of subdomains along the X-axis.

my:  - (In) Number of subdomains along the Y-axis.

mz:  - (In) Number of subdomains along the Z-axis.

lx:  - (In) Length of the cube along the X-axis.

ly:  - (In) Length of the cube along the Y-axis.

lz:  - (In) Length of the cube along the Z-axis.

The total number of processors must equal mx * my. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::~HexCubeGrid "virtual
Galeri::FiniteElements::HexCubeGrid::~HexCubeGrid()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumDimensions "virtual int
Galeri::FiniteElements::HexCubeGrid::NumDimensions() const

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumVerticesPerElement "virtual
int Galeri::FiniteElements::HexCubeGrid::NumVerticesPerElement() const

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::HexCubeGrid::NumFacesPerElement() const

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::HexCubeGrid::NumVerticesPerFace() const

Returns the number of vertices contained in each face. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementType "virtual string
Galeri::FiniteElements::HexCubeGrid::ElementType() const

Returns GALERI_HEX. ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::Comm "virtual const Epetra_Comm& Galeri::FiniteElements::HexCubeGrid::Comm()
const

Returns a reference to the communicator object. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyElements "virtual int
Galeri::FiniteElements::HexCubeGrid::NumMyElements() const

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalElements "virtual int
Galeri::FiniteElements::HexCubeGrid::NumGlobalElements() const

Returns the global number of finite elements. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyVertices "virtual int
Galeri::FiniteElements::HexCubeGrid::NumMyVertices() const

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalVertices "virtual int
Galeri::FiniteElements::HexCubeGrid::NumGlobalVertices() const

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::HexCubeGrid::NumMyBoundaryFaces() const

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalBoundaryFaces "virtual
int Galeri::FiniteElements::HexCubeGrid::NumGlobalBoundaryFaces()
const

Returns the global number of boundary faces. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::VertexCoord "virtual void
Galeri::FiniteElements::HexCubeGrid::VertexCoord(const int LocalID,
double *coord) const

Returns the coordinates of local vertex LocalVertex in vector coord.

Parameters:
-----------

LocalVertex:  - (In) Local ID of the vertex for whic coordinates are
required. Must be contained in the interval [0, NumMyVertices())

coord:  - (Out) double array of size 3. In output, contains the x-, y-
and z-coordinate of the specified vertex.

Parameter coord must be allocated of size 3 for both 2D and 3D
problems. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::VertexCoord "virtual void
Galeri::FiniteElements::HexCubeGrid::VertexCoord(const int Length,
const int *IDs, double *x, double *y, double *z) const

Returns the coordinates of specified local vertices.

Parameters:
-----------

Length:  - (In) Length of array IDs.

IDs:  - (In) Contains the list of vertices of which coordinates are
required.

x:  - (Out) double array of size Length. In output, contains the
x-coordinates of the specified vertices.

y:  - (Out) double array of size Length. In output, contains the
y-coordinates of the specified vertices.

z:  - (Out) double array of size Length. In output, contains the
z-coordinates of the specified vertices.

The z array must be allocated for both 2D and 3D problems. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementVertices "virtual void
Galeri::FiniteElements::HexCubeGrid::ElementVertices(const int
LocalID, int *elements) const

Returns the local vertex IDs of the specified local finite element.

Parameters:
-----------

LocalElement:  - (In) ID of the required local element.

elements:  - (Out) array of length NumElementVertices(), in output
will contain the local ID of the vertices of the specified element. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementMinLength "virtual double
Galeri::FiniteElements::HexCubeGrid::ElementMinLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementMaxLength "virtual double
Galeri::FiniteElements::HexCubeGrid::ElementMaxLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::RCPVertexMap "virtual const
RefCountPtr<Epetra_Map>
Galeri::FiniteElements::HexCubeGrid::RCPVertexMap() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::RCPElementMap "virtual const
RefCountPtr<Epetra_Map>
Galeri::FiniteElements::HexCubeGrid::RCPElementMap() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::VertexMap
"virtual const Epetra_Map&
Galeri::FiniteElements::HexCubeGrid::VertexMap() const

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::ElementMap
"virtual const Epetra_Map&
Galeri::FiniteElements::HexCubeGrid::ElementMap() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::RowMap "virtual const Epetra_Map&
Galeri::FiniteElements::HexCubeGrid::RowMap() const

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::Importer "virtual const Epetra_Import&
Galeri::FiniteElements::HexCubeGrid::Importer() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::ElementTag
"virtual int Galeri::FiniteElements::HexCubeGrid::ElementTag(const
int ID) const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::VertexTag
"virtual int Galeri::FiniteElements::HexCubeGrid::VertexTag(const int
ID) const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementVolume "virtual double
Galeri::FiniteElements::HexCubeGrid::ElementVolume() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::FaceVertices "virtual void
Galeri::FiniteElements::HexCubeGrid::FaceVertices(const int LocalFace,
int &tag, int *IDs) const

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::FacePatch
"int Galeri::FiniteElements::HexCubeGrid::FacePatch(const int
LocalFace) const

Returns the patch ID of the specified face.

Returns an integer ID that identifies the given boundary face as
belonging to a given part of the domain. It can be used by the user to
specify the value and the type of the boundary condition. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyElementsX "int
Galeri::FiniteElements::HexCubeGrid::NumMyElementsX() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyElementsY "int
Galeri::FiniteElements::HexCubeGrid::NumMyElementsY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyElementsXY "int
Galeri::FiniteElements::HexCubeGrid::NumMyElementsXY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyElementsZ "int
Galeri::FiniteElements::HexCubeGrid::NumMyElementsZ() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesX "int
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesX() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesY "int
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesXY "int
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesXY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesZ "int
Galeri::FiniteElements::HexCubeGrid::NumMyVerticesZ() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsX "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsX() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsY "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsXY "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsXY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsZ "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalElementsZ() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesX "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesX() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesY "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesXY "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesXY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesZ "int
Galeri::FiniteElements::HexCubeGrid::NumGlobalVerticesZ() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::LengthX "double Galeri::FiniteElements::HexCubeGrid::LengthX() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::LengthY "double Galeri::FiniteElements::HexCubeGrid::LengthY() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::LengthZ "double Galeri::FiniteElements::HexCubeGrid::LengthZ() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::DeltaX "double Galeri::FiniteElements::HexCubeGrid::DeltaX() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::DeltaY "double Galeri::FiniteElements::HexCubeGrid::DeltaY() const ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::DeltaZ "double Galeri::FiniteElements::HexCubeGrid::DeltaZ() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementVolume "virtual double
Galeri::FiniteElements::HexCubeGrid::ElementVolume(const int
LocalElement) const

Returns the volume of the specified local finite element.

Returns the area (in 2D) or the volume (in 3D) of the specified local
element ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::FaceArea "virtual double Galeri::FiniteElements::HexCubeGrid::FaceArea(const int
LocalFace) const

Returns the area of the specified local face.

Returns the length (in 2D) or the area (in 3D) of the specified
boundary face ";

%feature("docstring")  Galeri::FiniteElements::HexCubeGrid::MyVolume "virtual double Galeri::FiniteElements::HexCubeGrid::MyVolume() const

Returns the volume of all local elements. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::GlobalVolume "virtual double
Galeri::FiniteElements::HexCubeGrid::GlobalVolume() const

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumDomainsX "int
Galeri::FiniteElements::HexCubeGrid::NumDomainsX() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumDomainsY "int
Galeri::FiniteElements::HexCubeGrid::NumDomainsY() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumDomainsZ "int
Galeri::FiniteElements::HexCubeGrid::NumDomainsZ() const ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ExportToVertexMap "void
Galeri::FiniteElements::HexCubeGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ExportToRowMap "void
Galeri::FiniteElements::HexCubeGrid::ExportToRowMap(const
Epetra_DistObject &VertexObject, Epetra_DistObject &RowObject) const

Exports distributed object from VertexMap() to RowMap(). ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::NumNeighborsPerElement "virtual
int Galeri::FiniteElements::HexCubeGrid::NumNeighborsPerElement()
const

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::HexCubeGrid::ElementNeighbors "virtual void
Galeri::FiniteElements::HexCubeGrid::ElementNeighbors(const int
LocalElement, int *elements) const

Returns the local IDs of neighboring elements. ";


// File: classGaleri_1_1FiniteElements_1_1HexQuadrature.xml
%feature("docstring") Galeri::FiniteElements::HexQuadrature "

Quadrature formula on hexahedra.

Marzio Sala, SNL 9214.

Last updated on 02-Apr-05.

C++ includes: Galeri_HexQuadrature.h ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::HexQuadrature "Galeri::FiniteElements::HexQuadrature::HexQuadrature(const int
NumQuadrNodes)

Constructor.

Parameters:
-----------

NumQuadrNodes:  - (In) Number of quadrature nodes per element. Valid
choices are: 1. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::~HexQuadrature "Galeri::FiniteElements::HexQuadrature::~HexQuadrature() ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::ComputeJacobian "void
Galeri::FiniteElements::HexQuadrature::ComputeJacobian(const int
QuadrNode, const double *x_hex, const double *y_hex, const double
*z_hex) const

Computes the Jacobian at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::ComputeQuadrNodes "void
Galeri::FiniteElements::HexQuadrature::ComputeQuadrNodes(const int ii,
const double *x, const double *y, const double *z, double &xq, double
&yq, double &zq) const

Maps the quadrature nodes from the reference element to the actual
one. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::ComputeDerivatives "void
Galeri::FiniteElements::HexQuadrature::ComputeDerivatives(const int
QuadrNode) const

Computes the derivatives at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::QuadrWeight "double
Galeri::FiniteElements::HexQuadrature::QuadrWeight(const int
QuadrNode) const

Computes the weight at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::DetJacobian "double
Galeri::FiniteElements::HexQuadrature::DetJacobian(const int
QuadrNode) const

Computes the determinant of the Jacobian matrix at the quadrature
node. ";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::Phi "double Galeri::FiniteElements::HexQuadrature::Phi(const int i) const

Returns the value of the i-th basis function on the reference element.
";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::PhiX "double Galeri::FiniteElements::HexQuadrature::PhiX(const int i) const

Returns the value of the x-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::PhiY "double Galeri::FiniteElements::HexQuadrature::PhiY(const int i) const

Returns the value of the y-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::PhiZ "double Galeri::FiniteElements::HexQuadrature::PhiZ(const int i) const

Returns the value of the z-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::Psi "double Galeri::FiniteElements::HexQuadrature::Psi(const int i) const

Returns the value of the i-th test function on the reference element.
";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::PsiX "double Galeri::FiniteElements::HexQuadrature::PsiX(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::PsiY "double Galeri::FiniteElements::HexQuadrature::PsiY(const int i) const

Returns the value of the y-derivative i-th test function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::HexQuadrature::PsiZ "double Galeri::FiniteElements::HexQuadrature::PsiZ(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::NumQuadrNodes "int
Galeri::FiniteElements::HexQuadrature::NumQuadrNodes() const

Returns the number of quadrature node per element. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::NumPhiFunctions "int
Galeri::FiniteElements::HexQuadrature::NumPhiFunctions() const

Returns the number of basis function on the reference element. ";

%feature("docstring")
Galeri::FiniteElements::HexQuadrature::NumPsiFunctions "int
Galeri::FiniteElements::HexQuadrature::NumPsiFunctions() const

Returns the number of test function on the reference element. ";


// File: classGaleri_1_1FiniteElements_1_1LinearProblem.xml
%feature("docstring") Galeri::FiniteElements::LinearProblem "

Basic implementation of scalar finite element problem.

This class fill the linea system matrix (defined as an
Epetra_CrsMatrix), the right-hand side (defined as an
Epetra_MultiVector) and the starting solution (defined as a zero
Epetra_MultiVector). In the current implementation, only one rhs is
created.

Neumann boundary conditions are still to be fixed.

Marzio Sala, SNL 9214.

C++ includes: Galeri_LinearProblem.h ";

%feature("docstring")
Galeri::FiniteElements::LinearProblem::LinearProblem "Galeri::FiniteElements::LinearProblem::LinearProblem(const
AbstractGrid &Grid, const AbstractVariational &Variational,
Epetra_CrsMatrix &A, Epetra_MultiVector &LHS, Epetra_MultiVector &RHS)

Constructor.

Parameters:
-----------

Grid:  - (In) Reference to an AbstractGrid object

Variational:  - (In) Reference to an AbstractVariational object

A:  - (In/Out) Epetra_CrsMatrix, whose Map is Grid().RowMap(), that
will contain the linear system matrix.

LHS:  - (In/Out) Epetra_MultiVector, whose Map is Grid().RowMap(),
that will contain the starting solution (zero vector).

RHS:  - (In/Out) Epetra_MultiVector, whose Map is Grid().RowMap(),
that will contain the right-hand side. ";

%feature("docstring")
Galeri::FiniteElements::LinearProblem::~LinearProblem "virtual
Galeri::FiniteElements::LinearProblem::~LinearProblem()

Destructor. ";

%feature("docstring")  Galeri::FiniteElements::LinearProblem::Compute
"void Galeri::FiniteElements::LinearProblem::Compute()

Fills the linear system matrix and the right-hand side, zeros out the
solution. ";

%feature("docstring")
Galeri::FiniteElements::LinearProblem::ComputeNorms "void
Galeri::FiniteElements::LinearProblem::ComputeNorms(Epetra_MultiVector
&RowMatrixField, int(*ExactSolution)(double, double, double, double
*), const bool verbose=true, double *Solution=0, double *Exact=0,
double *Diff=0)

Computes L2, semi-H1 and H1 norms.

double xq, yq, zq; ";

%feature("docstring")  Galeri::FiniteElements::LinearProblem::A "virtual Epetra_RowMatrix& Galeri::FiniteElements::LinearProblem::A()

Returns a reference to the linear system matrix. ";

%feature("docstring")  Galeri::FiniteElements::LinearProblem::CrsA "virtual Epetra_CrsMatrix&
Galeri::FiniteElements::LinearProblem::CrsA()

Returns a reference to the linear system matrix as Epetra_CrsMatrix.
";

%feature("docstring")  Galeri::FiniteElements::LinearProblem::RHS "virtual Epetra_MultiVector&
Galeri::FiniteElements::LinearProblem::RHS()

Returns a reference to the multi-vector of right-hand side. ";

%feature("docstring")  Galeri::FiniteElements::LinearProblem::LHS "virtual Epetra_MultiVector&
Galeri::FiniteElements::LinearProblem::LHS()

Returns a reference to the multi-vector of starting solution. ";

%feature("docstring")  Galeri::FiniteElements::LinearProblem::Grid "virtual const AbstractGrid&
Galeri::FiniteElements::LinearProblem::Grid() const

Returns a reference to the grid object. ";

%feature("docstring")
Galeri::FiniteElements::LinearProblem::Variational "virtual const
AbstractVariational&
Galeri::FiniteElements::LinearProblem::Variational() const

Returns a reference to the variational object. ";


// File: classGaleri_1_1FiniteElements_1_1MEDITInterface.xml
%feature("docstring") Galeri::FiniteElements::MEDITInterface "";

%feature("docstring")
Galeri::FiniteElements::MEDITInterface::MEDITInterface "Galeri::FiniteElements::MEDITInterface::MEDITInterface(Epetra_Comm
&Comm) ";

%feature("docstring")
Galeri::FiniteElements::MEDITInterface::~MEDITInterface "Galeri::FiniteElements::MEDITInterface::~MEDITInterface() ";

%feature("docstring")  Galeri::FiniteElements::MEDITInterface::Comm "const Epetra_Comm& Galeri::FiniteElements::MEDITInterface::Comm()
const ";

%feature("docstring")  Galeri::FiniteElements::MEDITInterface::Write "void Galeri::FiniteElements::MEDITInterface::Write(const AbstractGrid
&data, const string &BaseName, const Epetra_MultiVector &Field)

int zzz = data.NumMyVertices(); ";


// File: classGaleri_1_1FiniteElements_1_1QuadQuadrature.xml
%feature("docstring") Galeri::FiniteElements::QuadQuadrature "

Quadrature formula on quadrilaterals.

Marzio Sala, SNL 9214.

Last updated on 31-Mar-05.

C++ includes: Galeri_QuadQuadrature.h ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::QuadQuadrature "Galeri::FiniteElements::QuadQuadrature::QuadQuadrature(const int
NumQuadrNodes)

Constructor.

Parameters:
-----------

NumQuadrNodes:  - (In) Number of quadrature nodes per element. Valid
choices are: 1, 4, 9. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::~QuadQuadrature "Galeri::FiniteElements::QuadQuadrature::~QuadQuadrature()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::ComputeJacobian "void
Galeri::FiniteElements::QuadQuadrature::ComputeJacobian(const int
QuadrNode, const double *x, const double *y, const double *z) const

Computes the Jacobian at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::ComputeQuadrNodes "void
Galeri::FiniteElements::QuadQuadrature::ComputeQuadrNodes(const int
QuadrNode, const double *x, const double *y, const double *z, double
&xq, double &yq, double &zq) const

Maps the quadrature nodes from the reference element to the actual
one. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::ComputeDerivatives "void
Galeri::FiniteElements::QuadQuadrature::ComputeDerivatives(const int
QuadrNode) const

Computes the derivatives at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::QuadrWeight "double
Galeri::FiniteElements::QuadQuadrature::QuadrWeight(const int
QuadrNode) const

Computes the weight at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::DetJacobian "double
Galeri::FiniteElements::QuadQuadrature::DetJacobian(const int
QuadrNode) const

Computes the determinant of the Jacobian matrix at the quadrature
node. ";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::Phi "double Galeri::FiniteElements::QuadQuadrature::Phi(const int i) const

Returns the value of the i-th basis function on the reference element.
";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::PhiX "double Galeri::FiniteElements::QuadQuadrature::PhiX(const int i) const

Returns the value of the x-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::PhiY "double Galeri::FiniteElements::QuadQuadrature::PhiY(const int i) const

Returns the value of the y-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::PhiZ "double Galeri::FiniteElements::QuadQuadrature::PhiZ(const int i) const

Returns the value of the z-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::Psi "double Galeri::FiniteElements::QuadQuadrature::Psi(const int i) const

Returns the value of the i-th test function on the reference element.
";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::PsiX "double Galeri::FiniteElements::QuadQuadrature::PsiX(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::PsiY "double Galeri::FiniteElements::QuadQuadrature::PsiY(const int i) const

Returns the value of the y-derivative i-th test function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::QuadQuadrature::PsiZ "double Galeri::FiniteElements::QuadQuadrature::PsiZ(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::NumQuadrNodes "int
Galeri::FiniteElements::QuadQuadrature::NumQuadrNodes() const

Returns the number of quadrature node per element. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::NumPhiFunctions "int
Galeri::FiniteElements::QuadQuadrature::NumPhiFunctions() const

Returns the number of basis function on the reference element. ";

%feature("docstring")
Galeri::FiniteElements::QuadQuadrature::NumPsiFunctions "int
Galeri::FiniteElements::QuadQuadrature::NumPsiFunctions() const

Returns the number of test function on the reference element. ";


// File: classGaleri_1_1FiniteElements_1_1QuadRectangleGrid.xml
%feature("docstring") Galeri::FiniteElements::QuadRectangleGrid "

Creates a grid with quadrilaterals on a rectangle.

Marzio Sala, SNL 9214.

C++ includes: Galeri_QuadRectangleGrid.h ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::QuadRectangleGrid "Galeri::FiniteElements::QuadRectangleGrid::QuadRectangleGrid(const
Epetra_Comm &Comm, const int nx, const int ny, const int mx, const int
my, const double lx=1.0, const double ly=1.0)

Constructor.

Parameters:
-----------

Comm:  - (In) Communicator object.

nx:  - (In) number of elements along the X-axis.

ny:  - (In) number of elements along the Y-axis.

mx:  - (In) Number of subdomains along the X-axis.

my:  - (In) Number of subdomains along the Y-axis.

lx:  - (In) Length of the rectangle along the X-axis.

ly:  - (In) Length of the rectangle along the Y-axis.

The total number of processors must equal mx * my. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::~QuadRectangleGrid "virtual
Galeri::FiniteElements::QuadRectangleGrid::~QuadRectangleGrid()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumDimensions "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumDimensions() const

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumVerticesPerElement "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumVerticesPerElement()
const

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumFacesPerElement() const

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumVerticesPerFace() const

Returns the number of vertices contained in each face. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementType "virtual
string Galeri::FiniteElements::QuadRectangleGrid::ElementType() const

Returns GALERI_QUAD. ";

%feature("docstring")  Galeri::FiniteElements::QuadRectangleGrid::Comm
"virtual const Epetra_Comm&
Galeri::FiniteElements::QuadRectangleGrid::Comm() const

Returns a reference to the communicator object. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyElements "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumMyElements() const

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalElements "virtual
int Galeri::FiniteElements::QuadRectangleGrid::NumGlobalElements()
const

Returns the global number of finite elements. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyVertices "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumMyVertices() const

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalVertices "virtual
int Galeri::FiniteElements::QuadRectangleGrid::NumGlobalVertices()
const

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumMyBoundaryFaces() const

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalBoundaryFaces "virtual int
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalBoundaryFaces()
const

Returns the global number of boundary faces. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::VertexCoord "virtual void
Galeri::FiniteElements::QuadRectangleGrid::VertexCoord(const int
LocalID, double *coord) const

Returns the coordinates of local vertex LocalVertex in vector coord.

Parameters:
-----------

LocalVertex:  - (In) Local ID of the vertex for whic coordinates are
required. Must be contained in the interval [0, NumMyVertices())

coord:  - (Out) double array of size 3. In output, contains the x-, y-
and z-coordinate of the specified vertex.

Parameter coord must be allocated of size 3 for both 2D and 3D
problems. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::VertexCoord "virtual void
Galeri::FiniteElements::QuadRectangleGrid::VertexCoord(const int
Length, const int *IDs, double *x, double *y, double *z) const

Returns the coordinates of specified local vertices.

Parameters:
-----------

Length:  - (In) Length of array IDs.

IDs:  - (In) Contains the list of vertices of which coordinates are
required.

x:  - (Out) double array of size Length. In output, contains the
x-coordinates of the specified vertices.

y:  - (Out) double array of size Length. In output, contains the
y-coordinates of the specified vertices.

z:  - (Out) double array of size Length. In output, contains the
z-coordinates of the specified vertices.

The z array must be allocated for both 2D and 3D problems. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementVertices "virtual
void Galeri::FiniteElements::QuadRectangleGrid::ElementVertices(const
int LocalID, int *elements) const

Returns the local vertex IDs of the specified local finite element.

Parameters:
-----------

LocalElement:  - (In) ID of the required local element.

elements:  - (Out) array of length NumElementVertices(), in output
will contain the local ID of the vertices of the specified element. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementMinLength "virtual
double
Galeri::FiniteElements::QuadRectangleGrid::ElementMinLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementMaxLength "virtual
double
Galeri::FiniteElements::QuadRectangleGrid::ElementMaxLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::RCPVertexMap "virtual
const RefCountPtr<Epetra_Map>
Galeri::FiniteElements::QuadRectangleGrid::RCPVertexMap() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::RCPElementMap "virtual
const RefCountPtr<Epetra_Map>
Galeri::FiniteElements::QuadRectangleGrid::RCPElementMap() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::VertexMap "virtual const
Epetra_Map& Galeri::FiniteElements::QuadRectangleGrid::VertexMap()
const

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementMap "virtual const
Epetra_Map& Galeri::FiniteElements::QuadRectangleGrid::ElementMap()
const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::BoundaryFaceMap "virtual
const Epetra_Map&
Galeri::FiniteElements::QuadRectangleGrid::BoundaryFaceMap() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::RowMap "virtual const
Epetra_Map& Galeri::FiniteElements::QuadRectangleGrid::RowMap() const

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::Importer "virtual const
Epetra_Import& Galeri::FiniteElements::QuadRectangleGrid::Importer()
const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementTag "virtual int
Galeri::FiniteElements::QuadRectangleGrid::ElementTag(const int ID)
const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::VertexTag "virtual int
Galeri::FiniteElements::QuadRectangleGrid::VertexTag(const int ID)
const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementVolume "virtual
double Galeri::FiniteElements::QuadRectangleGrid::ElementVolume()
const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::FaceVertices "virtual void
Galeri::FiniteElements::QuadRectangleGrid::FaceVertices(const int
LocalFace, int &tag, int *IDs) const

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::FacePatch "int
Galeri::FiniteElements::QuadRectangleGrid::FacePatch(const int
LocalFace) const

Returns the patch ID of the specified face.

Returns an integer ID that identifies the given boundary face as
belonging to a given part of the domain. It can be used by the user to
specify the value and the type of the boundary condition. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyElementsX "int
Galeri::FiniteElements::QuadRectangleGrid::NumMyElementsX() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyElementsY "int
Galeri::FiniteElements::QuadRectangleGrid::NumMyElementsY() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyVerticesX "int
Galeri::FiniteElements::QuadRectangleGrid::NumMyVerticesX() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumMyVerticesY "int
Galeri::FiniteElements::QuadRectangleGrid::NumMyVerticesY() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalElementsX "int
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalElementsX() const
";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalElementsY "int
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalElementsY() const
";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalVerticesX "int
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalVerticesX() const
";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalVerticesY "int
Galeri::FiniteElements::QuadRectangleGrid::NumGlobalVerticesY() const
";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::LengthX "double
Galeri::FiniteElements::QuadRectangleGrid::LengthX() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::LengthY "double
Galeri::FiniteElements::QuadRectangleGrid::LengthY() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::DeltaX "double
Galeri::FiniteElements::QuadRectangleGrid::DeltaX() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::DeltaY "double
Galeri::FiniteElements::QuadRectangleGrid::DeltaY() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementVolume "virtual
double Galeri::FiniteElements::QuadRectangleGrid::ElementVolume(const
int LocalElement) const

Returns the volume of the specified local finite element.

Returns the area (in 2D) or the volume (in 3D) of the specified local
element ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::FaceArea "virtual double
Galeri::FiniteElements::QuadRectangleGrid::FaceArea(const int
LocalFace) const

Returns the area of the specified local face.

Returns the length (in 2D) or the area (in 3D) of the specified
boundary face ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::MyVolume "virtual double
Galeri::FiniteElements::QuadRectangleGrid::MyVolume() const

Returns the volume of all local elements. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::GlobalVolume "virtual
double Galeri::FiniteElements::QuadRectangleGrid::GlobalVolume() const

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumDomainsX "int
Galeri::FiniteElements::QuadRectangleGrid::NumDomainsX() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumDomainsY "int
Galeri::FiniteElements::QuadRectangleGrid::NumDomainsY() const ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ExportToVertexMap "void
Galeri::FiniteElements::QuadRectangleGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ExportToRowMap "void
Galeri::FiniteElements::QuadRectangleGrid::ExportToRowMap(const
Epetra_DistObject &VertexObject, Epetra_DistObject &RowObject) const

Exports distributed object from VertexMap() to RowMap(). ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::NumNeighborsPerElement "int
Galeri::FiniteElements::QuadRectangleGrid::NumNeighborsPerElement()
const

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::QuadRectangleGrid::ElementNeighbors "void
Galeri::FiniteElements::QuadRectangleGrid::ElementNeighbors(int, int
*) const

Returns the local IDs of neighboring elements. ";


// File: classGaleri_1_1FiniteElements_1_1SUPGVariational.xml
%feature("docstring") Galeri::FiniteElements::SUPGVariational "

SUPG discretization of an advection-diffusion PDE.

This class performs the finite element discretization of a scalar,
advection-diffusion PDE, using the SUPG stabilization and the coth
formula for the definition of tau. This class works only with
triangles and tetrahedra.

Marzio Sala, SNL 9214.

C++ includes: Galeri_SUPGVariational.h ";

%feature("docstring")
Galeri::FiniteElements::SUPGVariational::SUPGVariational "Galeri::FiniteElements::SUPGVariational< T >::SUPGVariational(const
int NumQuadratureNodes, double(*diff)(const double &, const double &,
const double &), double(*bx)(const double &, const double &, const
double &), double(*by)(const double &, const double &, const double
&), double(*bz)(const double &, const double &, const double &),
double(*source)(const double &, const double &, const double &),
double(*force)(const double &, const double &, const double &),
double(*bc)(const double &, const double &, const double &, const int
&), int(*bc_type)(const int &))

Constructor. ";

%feature("docstring")
Galeri::FiniteElements::SUPGVariational::~SUPGVariational "Galeri::FiniteElements::SUPGVariational< T >::~SUPGVariational()

Destructor. ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::diff "double Galeri::FiniteElements::SUPGVariational< T >::diff(const double
x, const double y, const double z) const

Evaluates the diffusion coefficient at point (x, y, z). ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::source
"double Galeri::FiniteElements::SUPGVariational< T >::source(const
double x, const double y, const double z) const

Evaluates the source term at point (x, y, z). ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::force
"double Galeri::FiniteElements::SUPGVariational< T >::force(const
double x, const double y, const double z) const

Evaluates the force term at point (x, y, z). ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::conv_x
"double Galeri::FiniteElements::SUPGVariational< T >::conv_x(const
double x, const double y, const double z) const

Evaluates the x-component of the convective term at point (x, y, z).
";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::conv_y
"double Galeri::FiniteElements::SUPGVariational< T >::conv_y(const
double x, const double y, const double z) const

Evaluates the y-component of the convective term at point (x, y, z).
";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::conv_z
"double Galeri::FiniteElements::SUPGVariational< T >::conv_z(const
double x, const double y, const double z) const

Evaluates the z-component of the convective term at point (x, y, z).
";

%feature("docstring")
Galeri::FiniteElements::SUPGVariational::IntegrateOverElement "virtual int Galeri::FiniteElements::SUPGVariational< T
>::IntegrateOverElement(const AbstractVariational &Variational, const
double *x, const double *y, const double *z, const double *data,
double *ElementMatrix, double *ElementRHS) const

Integrates the bilinear form and the right-hand side over the element.
";

%feature("docstring")
Galeri::FiniteElements::SUPGVariational::ElementNorm "virtual int
Galeri::FiniteElements::SUPGVariational< T >::ElementNorm(const double
*LocalSol, const double *x, const double *y, const double *z, double
*Norm) const

Computes the norm of the computed solution over the element. ";

%feature("docstring")
Galeri::FiniteElements::SUPGVariational::ElementNorm "virtual int
Galeri::FiniteElements::SUPGVariational< T
>::ElementNorm(int(*ExactSolution)(double, double, double, double *),
const double *x, const double *y, const double *z, double *Norm) const

Computed the norm of the exact solution over the element. ";

%feature("docstring")
Galeri::FiniteElements::SUPGVariational::ElementNorm "virtual int
Galeri::FiniteElements::SUPGVariational< T >::ElementNorm(const double
*LocalSol, int(*ExactSolution)(double, double, double, double *),
const double *x, const double *y, const double *z, double *Norm) const

Computed the norm of the computed and exact solution over the element.
";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::LHS "double Galeri::FiniteElements::SUPGVariational< T >::LHS(const double
Phi, const double Psi, const double PhiX, const double PsiX, const
double PhiY, const double PsiY, const double PhiZ, const double PsiZ,
const double x, const double y, const double z) const

Evaluates the bilinear form (without integral) at point (x,y,z). ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::RHS "double Galeri::FiniteElements::SUPGVariational< T >::RHS(const double
Psi, const double PsiX, const double PsiY, const double PsiZ, const
double x, const double y, const double z) const

Returns the value of the right-hand side (without integral) at point
(x, y, z). ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::BC "int Galeri::FiniteElements::SUPGVariational< T >::BC(const int
PatchID) const

Returns an integer identifying the boundary condition assigned to the
specified patch. ";

%feature("docstring")  Galeri::FiniteElements::SUPGVariational::BC "double Galeri::FiniteElements::SUPGVariational< T >::BC(const double
x, const double y, const double z, const int Patch) const

Returns the value of the boundary condition at point (x, y, z). ";


// File: classGaleri_1_1FiniteElements_1_1TetCubeGrid.xml
%feature("docstring") Galeri::FiniteElements::TetCubeGrid "

Creates a grid with tetrahedral elements in a cube.

Marzio Sala, SNL 9214.

C++ includes: Galeri_TetCubeGrid.h ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::TetCubeGrid "Galeri::FiniteElements::TetCubeGrid::TetCubeGrid(const Epetra_Comm
&Comm, const int nx, const int ny, const int nz, const int mx, const
int my, const int mz, const double lx=1.0, const double ly=1.0, const
double lz=1.0)

Constructor.

Parameters:
-----------

Comm:  - (In) Communicator object.

nx:  - (In) number of elements along the X-axis.

ny:  - (In) number of elements along the Y-axis.

mx:  - (In) Number of subdomains along the X-axis.

my:  - (In) Number of subdomains along the Y-axis.

mz:  - (In) Number of subdomains along the Z-axis.

lx:  - (In) Length of the cube along the X-axis.

ly:  - (In) Length of the cube along the Y-axis.

lz:  - (In) Length of the cube along the Z-axis.

The total number of processors must equal mx * my * mz. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::~TetCubeGrid "virtual
Galeri::FiniteElements::TetCubeGrid::~TetCubeGrid()

Destructor. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumDimensions "virtual int
Galeri::FiniteElements::TetCubeGrid::NumDimensions() const

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumVerticesPerElement "virtual
int Galeri::FiniteElements::TetCubeGrid::NumVerticesPerElement() const

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::TetCubeGrid::NumFacesPerElement() const

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::TetCubeGrid::NumVerticesPerFace() const

Returns the number of vertices contained in each face. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementType "virtual string
Galeri::FiniteElements::TetCubeGrid::ElementType() const

Returns GALERI_TET. ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::Comm "virtual const Epetra_Comm& Galeri::FiniteElements::TetCubeGrid::Comm()
const

Returns a reference to the communicator object. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyElements "virtual int
Galeri::FiniteElements::TetCubeGrid::NumMyElements() const

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalElements "virtual int
Galeri::FiniteElements::TetCubeGrid::NumGlobalElements() const

Returns the global number of finite elements. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyVertices "virtual int
Galeri::FiniteElements::TetCubeGrid::NumMyVertices() const

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalVertices "virtual int
Galeri::FiniteElements::TetCubeGrid::NumGlobalVertices() const

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::TetCubeGrid::NumMyBoundaryFaces() const

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalBoundaryFaces "virtual
int Galeri::FiniteElements::TetCubeGrid::NumGlobalBoundaryFaces()
const

Returns the global number of boundary faces. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::VertexCoord "virtual void
Galeri::FiniteElements::TetCubeGrid::VertexCoord(const int LocalID,
double *coord) const

Returns the coordinates of local vertex LocalVertex in vector coord.

Parameters:
-----------

LocalVertex:  - (In) Local ID of the vertex for whic coordinates are
required. Must be contained in the interval [0, NumMyVertices())

coord:  - (Out) double array of size 3. In output, contains the x-, y-
and z-coordinate of the specified vertex.

Parameter coord must be allocated of size 3 for both 2D and 3D
problems. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::VertexCoord "virtual void
Galeri::FiniteElements::TetCubeGrid::VertexCoord(const int Length,
const int *IDs, double *x, double *y, double *z) const

Returns the coordinates of specified local vertices.

Parameters:
-----------

Length:  - (In) Length of array IDs.

IDs:  - (In) Contains the list of vertices of which coordinates are
required.

x:  - (Out) double array of size Length. In output, contains the
x-coordinates of the specified vertices.

y:  - (Out) double array of size Length. In output, contains the
y-coordinates of the specified vertices.

z:  - (Out) double array of size Length. In output, contains the
z-coordinates of the specified vertices.

The z array must be allocated for both 2D and 3D problems. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementVertices "virtual void
Galeri::FiniteElements::TetCubeGrid::ElementVertices(const int
LocalID, int *elements) const

Returns the local vertex IDs of the specified local finite element.

Parameters:
-----------

LocalElement:  - (In) ID of the required local element.

elements:  - (Out) array of length NumElementVertices(), in output
will contain the local ID of the vertices of the specified element. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementMinLength "virtual double
Galeri::FiniteElements::TetCubeGrid::ElementMinLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementMaxLength "virtual double
Galeri::FiniteElements::TetCubeGrid::ElementMaxLength(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::RCPVertexMap "virtual const
RefCountPtr<Epetra_Map>
Galeri::FiniteElements::TetCubeGrid::RCPVertexMap() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::RCPElementMap "virtual const
RefCountPtr<Epetra_Map>
Galeri::FiniteElements::TetCubeGrid::RCPElementMap() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::VertexMap
"virtual const Epetra_Map&
Galeri::FiniteElements::TetCubeGrid::VertexMap() const

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::ElementMap
"virtual const Epetra_Map&
Galeri::FiniteElements::TetCubeGrid::ElementMap() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::RowMap "virtual const Epetra_Map&
Galeri::FiniteElements::TetCubeGrid::RowMap() const

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::Importer "virtual const Epetra_Import&
Galeri::FiniteElements::TetCubeGrid::Importer() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::ElementTag
"virtual int Galeri::FiniteElements::TetCubeGrid::ElementTag(const
int ID) const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::VertexTag
"virtual int Galeri::FiniteElements::TetCubeGrid::VertexTag(const int
ID) const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementVolume "virtual double
Galeri::FiniteElements::TetCubeGrid::ElementVolume() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::FaceVertices "virtual void
Galeri::FiniteElements::TetCubeGrid::FaceVertices(const int LocalFace,
int &tag, int *IDs) const

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::FacePatch
"int Galeri::FiniteElements::TetCubeGrid::FacePatch(const int
LocalFace) const

Returns the patch ID of the specified face.

Returns an integer ID that identifies the given boundary face as
belonging to a given part of the domain. It can be used by the user to
specify the value and the type of the boundary condition. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumNeighborsPerElement "int
Galeri::FiniteElements::TetCubeGrid::NumNeighborsPerElement() const

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementNeighbors "void
Galeri::FiniteElements::TetCubeGrid::ElementNeighbors(int, int *)
const

Returns the local IDs of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyElementsX "int
Galeri::FiniteElements::TetCubeGrid::NumMyElementsX() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyElementsY "int
Galeri::FiniteElements::TetCubeGrid::NumMyElementsY() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyElementsZ "int
Galeri::FiniteElements::TetCubeGrid::NumMyElementsZ() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyVerticesX "int
Galeri::FiniteElements::TetCubeGrid::NumMyVerticesX() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyVerticesY "int
Galeri::FiniteElements::TetCubeGrid::NumMyVerticesY() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumMyVerticesZ "int
Galeri::FiniteElements::TetCubeGrid::NumMyVerticesZ() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalElementsX "int
Galeri::FiniteElements::TetCubeGrid::NumGlobalElementsX() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalElementsY "int
Galeri::FiniteElements::TetCubeGrid::NumGlobalElementsY() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalElementsZ "int
Galeri::FiniteElements::TetCubeGrid::NumGlobalElementsZ() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalVerticesX "int
Galeri::FiniteElements::TetCubeGrid::NumGlobalVerticesX() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalVerticesY "int
Galeri::FiniteElements::TetCubeGrid::NumGlobalVerticesY() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumGlobalVerticesZ "int
Galeri::FiniteElements::TetCubeGrid::NumGlobalVerticesZ() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::LengthX "double Galeri::FiniteElements::TetCubeGrid::LengthX() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::LengthY "double Galeri::FiniteElements::TetCubeGrid::LengthY() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::LengthZ "double Galeri::FiniteElements::TetCubeGrid::LengthZ() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::DeltaX "double Galeri::FiniteElements::TetCubeGrid::DeltaX() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::DeltaY "double Galeri::FiniteElements::TetCubeGrid::DeltaY() const ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::DeltaZ "double Galeri::FiniteElements::TetCubeGrid::DeltaZ() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ElementVolume "virtual double
Galeri::FiniteElements::TetCubeGrid::ElementVolume(const int
LocalElement) const

Returns the volume of the specified local finite element.

Returns the area (in 2D) or the volume (in 3D) of the specified local
element ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::FaceArea "virtual double Galeri::FiniteElements::TetCubeGrid::FaceArea(const int
LocalFace) const

Returns the area of the specified local face.

Returns the length (in 2D) or the area (in 3D) of the specified
boundary face ";

%feature("docstring")  Galeri::FiniteElements::TetCubeGrid::MyVolume "virtual double Galeri::FiniteElements::TetCubeGrid::MyVolume() const

Returns the volume of all local elements. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::GlobalVolume "virtual double
Galeri::FiniteElements::TetCubeGrid::GlobalVolume() const

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumDomainsX "int
Galeri::FiniteElements::TetCubeGrid::NumDomainsX() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumDomainsY "int
Galeri::FiniteElements::TetCubeGrid::NumDomainsY() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::NumDomainsZ "int
Galeri::FiniteElements::TetCubeGrid::NumDomainsZ() const ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ExportToVertexMap "void
Galeri::FiniteElements::TetCubeGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::TetCubeGrid::ExportToRowMap "void
Galeri::FiniteElements::TetCubeGrid::ExportToRowMap(const
Epetra_DistObject &VertexObject, Epetra_DistObject &RowObject) const

Exports distributed object from VertexMap() to RowMap(). ";


// File: classGaleri_1_1FiniteElements_1_1TetQuadrature.xml
%feature("docstring") Galeri::FiniteElements::TetQuadrature "

Quadrature formula on tetrahedra.

Marzio Sala, SNL 9214.

Last updated on Apr-05.

C++ includes: Galeri_TetQuadrature.h ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::TetQuadrature "Galeri::FiniteElements::TetQuadrature::TetQuadrature(const int
NumQuadrNodes)

Constructor.

Parameters:
-----------

NumQuadrNodes:  - (In) Number of quadrature nodes per element. Valid
choices are: 1. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::~TetQuadrature "Galeri::FiniteElements::TetQuadrature::~TetQuadrature() ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::ComputeJacobian "void
Galeri::FiniteElements::TetQuadrature::ComputeJacobian(const int
QuadrNode, const double *x, const double *y, const double *z) const

Computes the Jacobian at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::ComputeQuadrNodes "void
Galeri::FiniteElements::TetQuadrature::ComputeQuadrNodes(const int ii,
const double *x, const double *y, const double *z, double &xq, double
&yq, double &zq) const

Maps the quadrature nodes from the reference element to the actual
one. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::ComputeDerivatives "void
Galeri::FiniteElements::TetQuadrature::ComputeDerivatives(const int
QuadrNode) const

Computes the derivatives at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::QuadrWeight "double
Galeri::FiniteElements::TetQuadrature::QuadrWeight(const int
QuadrNode) const

Computes the weight at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::DetJacobian "double
Galeri::FiniteElements::TetQuadrature::DetJacobian(const int
QuadrNode) const

Computes the determinant of the Jacobian matrix at the quadrature
node. ";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::Phi "double Galeri::FiniteElements::TetQuadrature::Phi(const int i) const

Returns the value of the i-th basis function on the reference element.
";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::PhiX "double Galeri::FiniteElements::TetQuadrature::PhiX(const int i) const

Returns the value of the x-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::PhiY "double Galeri::FiniteElements::TetQuadrature::PhiY(const int i) const

Returns the value of the y-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::PhiZ "double Galeri::FiniteElements::TetQuadrature::PhiZ(const int i) const

Returns the value of the z-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::Psi "double Galeri::FiniteElements::TetQuadrature::Psi(const int i) const

Returns the value of the i-th test function on the reference element.
";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::PsiX "double Galeri::FiniteElements::TetQuadrature::PsiX(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::PsiY "double Galeri::FiniteElements::TetQuadrature::PsiY(const int i) const

Returns the value of the y-derivative i-th test function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::TetQuadrature::PsiZ "double Galeri::FiniteElements::TetQuadrature::PsiZ(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::NumQuadrNodes "int
Galeri::FiniteElements::TetQuadrature::NumQuadrNodes() const

Returns the number of quadrature node per element. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::NumPhiFunctions "int
Galeri::FiniteElements::TetQuadrature::NumPhiFunctions() const

Returns the number of basis function on the reference element. ";

%feature("docstring")
Galeri::FiniteElements::TetQuadrature::NumPsiFunctions "int
Galeri::FiniteElements::TetQuadrature::NumPsiFunctions() const

Returns the number of test function on the reference element. ";


// File: classGaleri_1_1FiniteElements_1_1TRIANGLEGrid.xml
%feature("docstring") Galeri::FiniteElements::TRIANGLEGrid "";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::TRIANGLEGrid "Galeri::FiniteElements::TRIANGLEGrid::TRIANGLEGrid(const Epetra_Comm
&Comm, const int NumPoints, const double *x, const double *y, const
double MaxArea) ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::~TRIANGLEGrid "Galeri::FiniteElements::TRIANGLEGrid::~TRIANGLEGrid() ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumDimensions "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumDimensions() const

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumVerticesPerElement "virtual
int Galeri::FiniteElements::TRIANGLEGrid::NumVerticesPerElement()
const

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumFacesPerElement() const

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumVerticesPerFace() const

Returns the number of vertices contained in each face. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ElementType "virtual string
Galeri::FiniteElements::TRIANGLEGrid::ElementType() const

Returns a string containing the element type. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumNeighborsPerElement "virtual
int Galeri::FiniteElements::TRIANGLEGrid::NumNeighborsPerElement()
const

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumMyElements "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumMyElements() const

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumGlobalElements "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumGlobalElements() const

Returns the global number of finite elements. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumMyVertices "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumMyVertices() const

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumGlobalVertices "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumGlobalVertices() const

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::TRIANGLEGrid::NumMyBoundaryFaces() const

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::NumGlobalBoundaryFaces "virtual
int Galeri::FiniteElements::TRIANGLEGrid::NumGlobalBoundaryFaces()
const

Returns the global number of boundary faces. ";

%feature("docstring")  Galeri::FiniteElements::TRIANGLEGrid::MyVolume
"virtual double Galeri::FiniteElements::TRIANGLEGrid::MyVolume()
const

Returns the volume of all local elements. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::GlobalVolume "virtual double
Galeri::FiniteElements::TRIANGLEGrid::GlobalVolume() const

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::VertexCoord "virtual void
Galeri::FiniteElements::TRIANGLEGrid::VertexCoord(const int
LocalVertex, double *coord) const

Returns the coordinates of local vertex LocalVertex in vector coord.
";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::VertexCoord "virtual void
Galeri::FiniteElements::TRIANGLEGrid::VertexCoord(const int Length,
const int *IDs, double *x, double *y, double *z) const

Returns the coordinates of specified local vertices. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ElementVertices "virtual void
Galeri::FiniteElements::TRIANGLEGrid::ElementVertices(const int
LocalElement, int *elements) const

Returns the local vertex IDs of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::FaceVertices "virtual void
Galeri::FiniteElements::TRIANGLEGrid::FaceVertices(const int
LocalFace, int &tag, int *IDs) const

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ElementNeighbors "virtual void
Galeri::FiniteElements::TRIANGLEGrid::ElementNeighbors(const int
LocalElement, int *elements) const

Returns the local IDs of neighboring elements. ";

%feature("docstring")  Galeri::FiniteElements::TRIANGLEGrid::FacePatch
"virtual int Galeri::FiniteElements::TRIANGLEGrid::FacePatch(const
int LocalFace) const

Returns the patch ID of the specified face. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ElementMinLength "virtual
double Galeri::FiniteElements::TRIANGLEGrid::ElementMinLength(const
int LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ElementMaxLength "virtual
double Galeri::FiniteElements::TRIANGLEGrid::ElementMaxLength(const
int LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ElementVolume "virtual double
Galeri::FiniteElements::TRIANGLEGrid::ElementVolume(const int
LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")  Galeri::FiniteElements::TRIANGLEGrid::FaceArea
"virtual double Galeri::FiniteElements::TRIANGLEGrid::FaceArea(const
int LocalFace) const

Returns the area of the specified local face. ";

%feature("docstring")  Galeri::FiniteElements::TRIANGLEGrid::VertexMap
"virtual const Epetra_Map&
Galeri::FiniteElements::TRIANGLEGrid::VertexMap() const

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")  Galeri::FiniteElements::TRIANGLEGrid::RowMap "virtual const Epetra_Map&
Galeri::FiniteElements::TRIANGLEGrid::RowMap() const

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ExportToVertexMap "virtual void
Galeri::FiniteElements::TRIANGLEGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::TRIANGLEGrid::ExportToRowMap "virtual void
Galeri::FiniteElements::TRIANGLEGrid::ExportToRowMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from VertexMap() to RowMap(). ";

%feature("docstring")  Galeri::FiniteElements::TRIANGLEGrid::Comm "virtual const Epetra_Comm&
Galeri::FiniteElements::TRIANGLEGrid::Comm() const

Returns a reference to the communicator object. ";


// File: classGaleri_1_1FiniteElements_1_1TriangleQuadrature.xml
%feature("docstring") Galeri::FiniteElements::TriangleQuadrature "

Quadrature formula on triangles.

Marzio Sala, SNL 9214.

Last updated on Apr-05.

C++ includes: Galeri_TriangleQuadrature.h ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::TriangleQuadrature "Galeri::FiniteElements::TriangleQuadrature::TriangleQuadrature(const
int NumQuadrNodes)

Constructor.

Parameters:
-----------

NumQuadrNodes:  - (In) Number of quadrature nodes per element. Valid
choices are: 1, 3, 4, 7. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::~TriangleQuadrature "Galeri::FiniteElements::TriangleQuadrature::~TriangleQuadrature()

Deastructor. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::ComputeJacobian "void
Galeri::FiniteElements::TriangleQuadrature::ComputeJacobian(const int
QuadrNode, const double *x_triangle, const double *y_triangle, const
double *z_triangle) const

Computes the Jacobian at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::ComputeQuadrNodes "void
Galeri::FiniteElements::TriangleQuadrature::ComputeQuadrNodes(const
int ii, const double *x, const double *y, const double *z, double &xq,
double &yq, double &zq) const

Maps the quadrature nodes from the reference element to the actual
one. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::ComputeDerivatives "void
Galeri::FiniteElements::TriangleQuadrature::ComputeDerivatives(const
int QuadrNode) const

Computes the derivatives at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::QuadrWeight "double
Galeri::FiniteElements::TriangleQuadrature::QuadrWeight(const int
QuadrNode) const

Computes the weight at the specified quadrature node. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::DetJacobian "double
Galeri::FiniteElements::TriangleQuadrature::DetJacobian(const int
QuadrNode) const

Computes the determinant of the Jacobian matrix at the quadrature
node. ";

%feature("docstring")  Galeri::FiniteElements::TriangleQuadrature::Phi
"double Galeri::FiniteElements::TriangleQuadrature::Phi(const int i)
const

Returns the value of the i-th basis function on the reference element.
";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::PhiX "double
Galeri::FiniteElements::TriangleQuadrature::PhiX(const int i) const

Returns the value of the x-derivative i-th basis function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::PhiY "double
Galeri::FiniteElements::TriangleQuadrature::PhiY(const int i) const

Returns the value of the y-derivative i-th basis function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::PhiZ "double
Galeri::FiniteElements::TriangleQuadrature::PhiZ(const int i) const

Returns the value of the z-derivative i-th basis function on the
reference element. ";

%feature("docstring")  Galeri::FiniteElements::TriangleQuadrature::Psi
"double Galeri::FiniteElements::TriangleQuadrature::Psi(const int i)
const

Returns the value of the i-th test function on the reference element.
";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::PsiX "double
Galeri::FiniteElements::TriangleQuadrature::PsiX(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::PsiY "double
Galeri::FiniteElements::TriangleQuadrature::PsiY(const int i) const

Returns the value of the y-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::PsiZ "double
Galeri::FiniteElements::TriangleQuadrature::PsiZ(const int i) const

Returns the value of the z-derivative i-th test function on the
reference element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::NumQuadrNodes "int
Galeri::FiniteElements::TriangleQuadrature::NumQuadrNodes() const

Returns the number of quadrature node per element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::NumPhiFunctions "int
Galeri::FiniteElements::TriangleQuadrature::NumPhiFunctions() const

Returns the number of basis function on the reference element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleQuadrature::NumPsiFunctions "int
Galeri::FiniteElements::TriangleQuadrature::NumPsiFunctions() const

Returns the number of test function on the reference element. ";


// File: classGaleri_1_1FiniteElements_1_1TriangleRectangleGrid.xml
%feature("docstring") Galeri::FiniteElements::TriangleRectangleGrid "

Creates a grid composed by triangles, the domain is a rectangle.

This class defined, on-the-fly, the triangulation of a 2D rectangular
domain. The elements are all triangles. For parallel run, the
rectangle is subdivided along the X- and Y-axis, as specified by the
user.

Marzio Sala, SNL 9214.

C++ includes: Galeri_TriangleRectangleGrid.h ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::TriangleRectangleGrid "Galeri::FiniteElements::TriangleRectangleGrid::TriangleRectangleGrid(const
Epetra_Comm &Comm, const int nx, const int ny, const int mx, const int
my, const double lx=1.0, const double ly=1.0)

Constructor.

Parameters:
-----------

Comm:  - (In) Communicator object.

nx:  - (In) number of elements along the X-axis.

ny:  - (In) number of elements along the Y-axis.

mx:  - (In) Number of subdomains along the X-axis.

my:  - (In) Number of subdomains along the Y-axis.

lx:  - (In) Length of the rectangle along the X-axis.

ly:  - (In) Length of the rectangle along the Y-axis.

The total number of processors must equal mx * my. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::~TriangleRectangleGrid
"virtual
Galeri::FiniteElements::TriangleRectangleGrid::~TriangleRectangleGrid()
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumDimensions "virtual
int Galeri::FiniteElements::TriangleRectangleGrid::NumDimensions()
const

Returns the number of dimensions of the grid. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumVerticesPerElement "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumVerticesPerElement()
const

Returns the number of vertices contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumFacesPerElement "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumFacesPerElement()
const

Returns the number of faces contained in each element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumVerticesPerFace "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumVerticesPerFace()
const

Returns the number of vertices contained in each face. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementType "virtual
string Galeri::FiniteElements::TriangleRectangleGrid::ElementType()
const

Returns GALERI_TRIANGLE. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::Comm "virtual const
Epetra_Comm& Galeri::FiniteElements::TriangleRectangleGrid::Comm()
const

Returns a reference to the communicator object. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyElements "virtual
int Galeri::FiniteElements::TriangleRectangleGrid::NumMyElements()
const

Returns the number of finite elements on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalElements "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalElements()
const

Returns the global number of finite elements. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyVertices "virtual
int Galeri::FiniteElements::TriangleRectangleGrid::NumMyVertices()
const

Returns the number of vertices on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalVertices "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalVertices()
const

Returns the global number of vertices. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyBoundaryFaces "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumMyBoundaryFaces()
const

Returns the number of boundary faces on the calling process. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalBoundaryFaces
"virtual int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalBoundaryFaces()
const

Returns the global number of boundary faces. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::VertexCoord "virtual
void Galeri::FiniteElements::TriangleRectangleGrid::VertexCoord(const
int LocalID, double *coord) const

Returns the coordinates of local vertex LocalVertex in vector coord.

Parameters:
-----------

LocalVertex:  - (In) Local ID of the vertex for whic coordinates are
required. Must be contained in the interval [0, NumMyVertices())

coord:  - (Out) double array of size 3. In output, contains the x-, y-
and z-coordinate of the specified vertex.

Parameter coord must be allocated of size 3 for both 2D and 3D
problems. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::VertexCoord "virtual
void Galeri::FiniteElements::TriangleRectangleGrid::VertexCoord(const
int Length, const int *IDs, double *x, double *y, double *z) const

Returns the coordinates of specified local vertices.

Parameters:
-----------

Length:  - (In) Length of array IDs.

IDs:  - (In) Contains the list of vertices of which coordinates are
required.

x:  - (Out) double array of size Length. In output, contains the
x-coordinates of the specified vertices.

y:  - (Out) double array of size Length. In output, contains the
y-coordinates of the specified vertices.

z:  - (Out) double array of size Length. In output, contains the
z-coordinates of the specified vertices.

The z array must be allocated for both 2D and 3D problems. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementVertices "virtual void
Galeri::FiniteElements::TriangleRectangleGrid::ElementVertices(const
int LocalID, int *elements) const

Returns the local vertex IDs of the specified local finite element.

Parameters:
-----------

LocalElement:  - (In) ID of the required local element.

elements:  - (Out) array of length NumElementVertices(), in output
will contain the local ID of the vertices of the specified element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementMinLength "virtual double
Galeri::FiniteElements::TriangleRectangleGrid::ElementMinLength(const
int LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementMaxLength "virtual double
Galeri::FiniteElements::TriangleRectangleGrid::ElementMaxLength(const
int LocalElement) const

Returns the volume of the specified local finite element. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::RCPVertexMap "virtual
const RefCountPtr<Epetra_Map>
Galeri::FiniteElements::TriangleRectangleGrid::RCPVertexMap() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::RCPElementMap "virtual
const RefCountPtr<Epetra_Map>
Galeri::FiniteElements::TriangleRectangleGrid::RCPElementMap() const
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::VertexMap "virtual
const Epetra_Map&
Galeri::FiniteElements::TriangleRectangleGrid::VertexMap() const

Returns a reference to the map representing the vertex distribution.
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementMap "virtual
const Epetra_Map&
Galeri::FiniteElements::TriangleRectangleGrid::ElementMap() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::RowMap "virtual const
Epetra_Map& Galeri::FiniteElements::TriangleRectangleGrid::RowMap()
const

Returns a reference to the map representing the distribution of rows.
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::Importer "virtual
const Epetra_Import&
Galeri::FiniteElements::TriangleRectangleGrid::Importer() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementTag "virtual
int Galeri::FiniteElements::TriangleRectangleGrid::ElementTag(const
int ID) const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::VertexTag "virtual int
Galeri::FiniteElements::TriangleRectangleGrid::VertexTag(const int ID)
const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementVolume "virtual
double Galeri::FiniteElements::TriangleRectangleGrid::ElementVolume()
const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::FaceVertices "virtual
void Galeri::FiniteElements::TriangleRectangleGrid::FaceVertices(const
int LocalFace, int &tag, int *IDs) const

Returns the local vertex IDs of vertices contained in the specified
boundary face. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::FacePatch "int
Galeri::FiniteElements::TriangleRectangleGrid::FacePatch(const int
LocalFace) const

Returns the patch ID of the specified face.

Returns an integer ID that identifies the given boundary face as
belonging to a given part of the domain. It can be used by the user to
specify the value and the type of the boundary condition. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyElementsX "int
Galeri::FiniteElements::TriangleRectangleGrid::NumMyElementsX() const
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyElementsY "int
Galeri::FiniteElements::TriangleRectangleGrid::NumMyElementsY() const
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyVerticesX "int
Galeri::FiniteElements::TriangleRectangleGrid::NumMyVerticesX() const
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumMyVerticesY "int
Galeri::FiniteElements::TriangleRectangleGrid::NumMyVerticesY() const
";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalElementsX "int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalElementsX()
const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalElementsY "int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalElementsY()
const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalVerticesX "int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalVerticesX()
const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalVerticesY "int
Galeri::FiniteElements::TriangleRectangleGrid::NumGlobalVerticesY()
const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::LengthX "double
Galeri::FiniteElements::TriangleRectangleGrid::LengthX() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::LengthY "double
Galeri::FiniteElements::TriangleRectangleGrid::LengthY() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::DeltaX "double
Galeri::FiniteElements::TriangleRectangleGrid::DeltaX() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::DeltaY "double
Galeri::FiniteElements::TriangleRectangleGrid::DeltaY() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementVolume "virtual
double
Galeri::FiniteElements::TriangleRectangleGrid::ElementVolume(const int
LocalElement) const

Returns the volume of the specified local finite element.

Returns the area (in 2D) or the volume (in 3D) of the specified local
element ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::FaceArea "virtual
double Galeri::FiniteElements::TriangleRectangleGrid::FaceArea(const
int LocalFace) const

Returns the area of the specified local face.

Returns the length (in 2D) or the area (in 3D) of the specified
boundary face ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::MyVolume "virtual
double Galeri::FiniteElements::TriangleRectangleGrid::MyVolume() const

Returns the volume of all local elements. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::GlobalVolume "virtual
double Galeri::FiniteElements::TriangleRectangleGrid::GlobalVolume()
const

Returns the global volume of the grid. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumDomainsX "int
Galeri::FiniteElements::TriangleRectangleGrid::NumDomainsX() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumDomainsY "int
Galeri::FiniteElements::TriangleRectangleGrid::NumDomainsY() const ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ExportToVertexMap "void
Galeri::FiniteElements::TriangleRectangleGrid::ExportToVertexMap(const
Epetra_DistObject &RowObject, Epetra_DistObject &VertexObject) const

Exports distributed object from RowMap() to VertexMap(). ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ExportToRowMap "void
Galeri::FiniteElements::TriangleRectangleGrid::ExportToRowMap(const
Epetra_DistObject &VertexObject, Epetra_DistObject &RowObject) const

Exports distributed object from VertexMap() to RowMap(). ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::NumNeighborsPerElement
"int
Galeri::FiniteElements::TriangleRectangleGrid::NumNeighborsPerElement()
const

Returns the number of neighboring elements. ";

%feature("docstring")
Galeri::FiniteElements::TriangleRectangleGrid::ElementNeighbors "void
Galeri::FiniteElements::TriangleRectangleGrid::ElementNeighbors(int,
int *) const

Returns the local IDs of neighboring elements. ";


// File: namespaceGaleri.xml
%feature("docstring")  Galeri::FiniteElements::CreateCrsMatrix "Epetra_CrsMatrix * Galeri::CreateCrsMatrix(const string MatrixType,
const Epetra_Map *Map, Teuchos::ParameterList &List) ";

%feature("docstring")  Galeri::FiniteElements::CreateMap "Epetra_Map
* Galeri::CreateMap(string MapType, Epetra_Comm &Comm,
Teuchos::ParameterList &List) ";

%feature("docstring")  Galeri::FiniteElements::CSRCSC "int
Galeri::CSRCSC(int n, int n2, int job, int ipos, double *a, int *ja,
int *ia, double *ao, int *jao, int *iao) ";

%feature("docstring")  Galeri::FiniteElements::SSRCSR "int
Galeri::SSRCSR(int job, int value2, int nrow, double *a, int *ja, int
*ia, int nzmax, double *ao, int *jao, int *iao, int *indu, int *iwk)
";

%feature("docstring")  Galeri::FiniteElements::SCSCMV "void
Galeri::SCSCMV(int isym, int m, int n, double *val, int *indx, int
*pntr, double *x, double *y) ";

%feature("docstring")  Galeri::FiniteElements::SCSCRES "double
Galeri::SCSCRES(int isym, int m, int n, double *val, int *indx, int
*pntr, double *x, double *b) ";

%feature("docstring")  Galeri::FiniteElements::ReadHB "void
Galeri::ReadHB(const char *data_file, const Epetra_Comm &comm,
Epetra_Map *&map, Epetra_CrsMatrix *&A, Epetra_Vector *&x,
Epetra_Vector *&b, Epetra_Vector *&xexact) ";

%feature("docstring")  Galeri::FiniteElements::Solve "void
Galeri::Solve(const Epetra_LinearProblem Problem) ";

%feature("docstring")  Galeri::FiniteElements::Solve "void
Galeri::Solve(const Epetra_RowMatrix *Matrix, const Epetra_MultiVector
*LHS, const Epetra_MultiVector *RHS) ";

%feature("docstring")  Galeri::FiniteElements::ComputeNorm "double
Galeri::ComputeNorm(const Epetra_MultiVector *LHS, const
Epetra_MultiVector *RHS) ";

%feature("docstring")  Galeri::FiniteElements::ComputeNorm "double
Galeri::ComputeNorm(const Epetra_RowMatrix *A, const
Epetra_MultiVector *LHS, const Epetra_MultiVector *RHS) ";

%feature("docstring")
Galeri::FiniteElements::CreateCartesianCoordinates "Epetra_MultiVector * Galeri::CreateCartesianCoordinates(const string
CoordType, const Epetra_BlockMap *BlockMap, Teuchos::ParameterList
&List) ";

%feature("docstring")  Galeri::FiniteElements::toString "string
Galeri::toString(const int &x) ";

%feature("docstring")  Galeri::FiniteElements::toString "string
Galeri::toString(const unsigned int &x) ";

%feature("docstring")  Galeri::FiniteElements::toString "string
Galeri::toString(const double &x) ";

%feature("docstring")
Galeri::FiniteElements::GetNeighboursCartesian2d "void
Galeri::GetNeighboursCartesian2d(const int i, const int nx, const int
ny, int &left, int &right, int &lower, int &upper) ";

%feature("docstring")
Galeri::FiniteElements::GetNeighboursCartesian2d "void
Galeri::GetNeighboursCartesian2d(const int i, const int nx, const int
ny, int &left, int &right, int &lower, int &upper, int &left2, int
&right2, int &lower2, int &upper2) ";

%feature("docstring")
Galeri::FiniteElements::GetNeighboursCartesian3d "void
Galeri::GetNeighboursCartesian3d(const int i, const int nx, const int
ny, const int nz, int &left, int &right, int &lower, int &upper, int
&below, int &above) ";

%feature("docstring")  Galeri::FiniteElements::PrintStencil2D "void
Galeri::PrintStencil2D(const Epetra_CrsMatrix *Matrix, const int nx,
const int ny, int GID) ";

%feature("docstring")  Galeri::FiniteElements::CreateVbrMatrix "Epetra_VbrMatrix * Galeri::CreateVbrMatrix(const Epetra_CrsMatrix
*CrsMatrix, const int NumPDEs) ";


// File: namespaceGaleri_1_1FiniteElements.xml
%feature("docstring")  Galeri::FiniteElements::Length "double
Galeri::FiniteElements::Length(const double x1, const double y1, const
double z1, const double x2, const double y2, const double z2)

Returns the distance between two points in space. ";

%feature("docstring")  Galeri::FiniteElements::Length "double
Galeri::FiniteElements::Length(const double *x, const double *y, const
double *z)

Returns the distance between two points in space. ";

%feature("docstring")  Galeri::FiniteElements::AreaOfTriangle "double
Galeri::FiniteElements::AreaOfTriangle(const double *x, const double
*y, const double *z)

Computes the area of a triangle in space. ";

%feature("docstring")  Galeri::FiniteElements::AreaOfQuad "double
Galeri::FiniteElements::AreaOfQuad(const double *x, const double *y,
const double *z)

Computes the are of a quadrilateral in space. ";

%feature("docstring")  Galeri::FiniteElements::VolumeOfTet "double
Galeri::FiniteElements::VolumeOfTet(const double *X, const double *Y,
const double *Z)

Computes the volume of a tetrahedron. ";


// File: namespaceGaleri_1_1Maps.xml
%feature("docstring")  Galeri::Maps::Cartesian2D "Epetra_Map*
Galeri::Maps::Cartesian2D(const Epetra_Comm &Comm, const int nx, const
int ny, const int mx, const int my) ";

%feature("docstring")  Galeri::Maps::Cartesian3D "Epetra_Map*
Galeri::Maps::Cartesian3D(const Epetra_Comm &Comm, const int nx, const
int ny, const int nz, const int mx, const int my, const int mz) ";

%feature("docstring")  Galeri::Maps::Interlaced "Epetra_Map*
Galeri::Maps::Interlaced(Epetra_Comm &Comm, int NumGlobalElements) ";

%feature("docstring")  Galeri::Maps::Linear "Epetra_Map*
Galeri::Maps::Linear(Epetra_Comm &Comm, int NumGlobalElements) ";

%feature("docstring")  Galeri::Maps::NodeCartesian2D "Epetra_Map*
Galeri::Maps::NodeCartesian2D(const Epetra_Comm &Comm, const
Epetra_Comm &NodeComm, const int MyNodeID, const int nx, const int ny,
const int ndx, const int ndy, const int px, const int py) ";

%feature("docstring")  Galeri::Maps::Random "Epetra_Map*
Galeri::Maps::Random(const Epetra_Comm &Comm, const int n) ";


// File: namespaceGaleri_1_1Matrices.xml
%feature("docstring")  Galeri::Matrices::BentPipe2D "Epetra_CrsMatrix* Galeri::Matrices::BentPipe2D(const Epetra_Map *Map,
const int nx, const int ny, const double lx, const double ly, const
double conv, const double diff) ";

%feature("docstring")  Galeri::Matrices::BigCross2D "Epetra_CrsMatrix* Galeri::Matrices::BigCross2D(const Epetra_Map *Map,
const int nx, const int ny, const double a, const double b, const
double c, const double d, const double e, const double bb, const
double cc, const double dd, const double ee) ";

%feature("docstring")  Galeri::Matrices::BigStar2D "Epetra_CrsMatrix*
Galeri::Matrices::BigStar2D(const Epetra_Map *Map, const int nx, const
int ny, const double a, const double b, const double c, const double
d, const double e, const double z1, const double z2, const double z3,
const double z4, const double bb, const double cc, const double dd,
const double ee) ";

%feature("docstring")  Galeri::Matrices::Cauchy "Epetra_CrsMatrix*
Galeri::Matrices::Cauchy(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::Cross2D "Epetra_CrsMatrix*
Galeri::Matrices::Cross2D(const Epetra_Map *Map, const int nx, const
int ny, const double a, const double b, const double c, const double
d, const double e) ";

%feature("docstring")  Galeri::Matrices::Cross2D "Epetra_CrsMatrix*
Galeri::Matrices::Cross2D(const Epetra_Map *Map, const int nx, const
int ny, const Epetra_Vector &A, const Epetra_Vector &B, const
Epetra_Vector &C, const Epetra_Vector &D, const Epetra_Vector &E) ";

%feature("docstring")  Galeri::Matrices::Cross3D "Epetra_CrsMatrix*
Galeri::Matrices::Cross3D(const Epetra_Map *Map, const int nx, const
int ny, const int nz, const double a, const double b, const double c,
const double d, const double e, const double f, const double g) ";

%feature("docstring")  Galeri::Matrices::Diag "Epetra_CrsMatrix*
Galeri::Matrices::Diag(const Epetra_Map *Map, double Value) ";

%feature("docstring")  Galeri::Matrices::Diag "Epetra_CrsMatrix*
Galeri::Matrices::Diag(const Epetra_Map *Map, Epetra_Vector &Vector)
";

%feature("docstring")  Galeri::Matrices::Fielder "Epetra_CrsMatrix*
Galeri::Matrices::Fielder(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::Hanowa "Epetra_CrsMatrix*
Galeri::Matrices::Hanowa(const Epetra_Map *Map, const double value) ";

%feature("docstring")  Galeri::Matrices::Hilbert "Epetra_CrsMatrix*
Galeri::Matrices::Hilbert(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::JordanBlock "Epetra_CrsMatrix* Galeri::Matrices::JordanBlock(const Epetra_Map *Map,
const double value) ";

%feature("docstring")  Galeri::Matrices::KMS "Epetra_CrsMatrix*
Galeri::Matrices::KMS(const Epetra_Map *Map, const double value) ";

%feature("docstring")  Galeri::Matrices::Laplace1DNeumann "Epetra_CrsMatrix* Galeri::Matrices::Laplace1DNeumann(const Epetra_Map
*Map) ";

%feature("docstring")  Galeri::Matrices::Lehmer "Epetra_CrsMatrix*
Galeri::Matrices::Lehmer(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::Minij "Epetra_CrsMatrix*
Galeri::Matrices::Minij(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::Ones "Epetra_CrsMatrix*
Galeri::Matrices::Ones(const Epetra_Map *Map, const double value) ";

%feature("docstring")  Galeri::Matrices::Parter "Epetra_CrsMatrix*
Galeri::Matrices::Parter(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::Pei "Epetra_CrsMatrix*
Galeri::Matrices::Pei(const Epetra_Map *Map, const double value) ";

%feature("docstring")  Galeri::Matrices::Recirc2D "Epetra_CrsMatrix*
Galeri::Matrices::Recirc2D(const Epetra_Map *Map, const int nx, const
int ny, const double lx, const double ly, const double conv, const
double diff) ";

%feature("docstring")  Galeri::Matrices::Ris "Epetra_CrsMatrix*
Galeri::Matrices::Ris(const Epetra_Map *Map) ";

%feature("docstring")  Galeri::Matrices::Star2D "Epetra_CrsMatrix*
Galeri::Matrices::Star2D(const Epetra_Map *Map, const int nx, const
int ny, const double a, const double b, const double c, const double
d, const double e, const double z1, const double z2, const double z3,
const double z4) ";

%feature("docstring")  Galeri::Matrices::Stretched2D "Epetra_CrsMatrix* Galeri::Matrices::Stretched2D(const Epetra_Map *Map,
const int nx, const int ny, const double epsilon) ";

%feature("docstring")  Galeri::Matrices::Tridiag "Epetra_CrsMatrix*
Galeri::Matrices::Tridiag(const Epetra_Map *Map, const double a, const
double b, const double c) ";

%feature("docstring")  Galeri::Matrices::UniFlow2D "Epetra_CrsMatrix*
Galeri::Matrices::UniFlow2D(const Epetra_Map *Map, const int nx, const
int ny, const double lx, const double ly, const double conv, const
double diff, const double alpha) ";

%feature("docstring")  Galeri::Matrices::Vander "Epetra_CrsMatrix*
Galeri::Matrices::Vander(const Epetra_Map *Map, const double value) ";


// File: namespacestd.xml


// File: namespaceTeuchos.xml


// File: Galeri__AbstractGrid_8h.xml


// File: Galeri__AbstractProblem_8h.xml


// File: Galeri__AbstractQuadrature_8h.xml


// File: Galeri__AbstractVariational_8h.xml


// File: Galeri__BentPipe2D_8h.xml


// File: Galeri__BigCross2D_8h.xml


// File: Galeri__BigStar2D_8h.xml


// File: Galeri__Cartesian2D_8h.xml


// File: Galeri__Cartesian3D_8h.xml


// File: Galeri__Cauchy_8h.xml


// File: Galeri__ConfigDefs_8h.xml


// File: Galeri__Cross2D_8h.xml


// File: Galeri__Cross3D_8h.xml


// File: Galeri__CrsMatrices_8cpp.xml


// File: Galeri__CrsMatrices_8h.xml


// File: Galeri__Diag_8h.xml


// File: Galeri__Exception_8h.xml


// File: Galeri__Fielder_8h.xml


// File: Galeri__FileGrid_8h.xml


// File: Galeri__FiniteElements_8h.xml


// File: Galeri__GalerkinVariational_8h.xml


// File: Galeri__Hanowa_8h.xml


// File: Galeri__HexCubeGrid_8h.xml


// File: Galeri__HexQuadrature_8h.xml


// File: Galeri__Hilbert_8h.xml


// File: Galeri__Interlaced_8h.xml


// File: Galeri__JordanBlock_8h.xml


// File: Galeri__KMS_8h.xml


// File: Galeri__Laplace1DNeumann_8h.xml


// File: Galeri__Lehmer_8h.xml


// File: Galeri__Linear_8h.xml


// File: Galeri__LinearProblem_8h.xml


// File: Galeri__Maps_8cpp.xml


// File: Galeri__Maps_8h.xml


// File: Galeri__MEDITInterface_8h.xml


// File: Galeri__Minij_8h.xml


// File: Galeri__NodeCartesian2D_8h.xml


// File: Galeri__Ones_8h.xml


// File: Galeri__Parter_8h.xml


// File: Galeri__Pei_8h.xml


// File: Galeri__QuadQuadrature_8h.xml


// File: Galeri__QuadRectangleGrid_8h.xml


// File: Galeri__Random_8h.xml


// File: Galeri__ReadHB_8cpp.xml


// File: Galeri__ReadHB_8h.xml


// File: Galeri__Recirc2D_8h.xml


// File: Galeri__Ris_8h.xml


// File: Galeri__Star2D_8h.xml


// File: Galeri__Stretched2D_8h.xml


// File: Galeri__SUPGVariational_8h.xml


// File: Galeri__TetCubeGrid_8h.xml


// File: Galeri__TetQuadrature_8h.xml


// File: Galeri__TRIANGLEGrid_8h.xml


// File: Galeri__TriangleQuadrature_8h.xml


// File: Galeri__TriangleRectangleGrid_8h.xml


// File: Galeri__Tridiag_8h.xml


// File: Galeri__UniFlow2D_8h.xml


// File: Galeri__Utils_8cpp.xml


// File: Galeri__Utils_8h.xml


// File: Galeri__Vander_8h.xml


// File: Galeri__VbrMatrices_8cpp.xml


// File: Galeri__VbrMatrices_8h.xml


// File: Galeri__Version_8h.xml
%feature("docstring")  Galeri_Version "string Galeri_Version() ";


// File: Galeri__Workspace_8h.xml


// File: Laplacian3D.xml


// File: AdvDiff2D.xml


// File: ml_MLAPI.xml


// File: dir_83c6dce543ad695bb9704f6bf70daa7f.xml


// File: dir_da9a02cbd345ffddae58604b9dfe1707.xml


// File: dir_3e4054923dec99bc1d6b48fcf5433e95.xml


// File: dir_467b95ce87aa0ef5cd28c6f7bf50ec25.xml


// File: dir_ffaff38ab9f6978c2d165634937e41a3.xml


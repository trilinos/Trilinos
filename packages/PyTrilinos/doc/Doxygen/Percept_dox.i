
// File: index.xml

// File: classA.xml
%feature("docstring") A "";


// File: classA_3_01int_00_01U_01_4.xml
%feature("docstring") A< int, U > " ";


// File: classB.xml
%feature("docstring") B "";


// File: classstk_1_1percept_1_1IntrepidManager_1_1Bases.xml
%feature("docstring") stk::percept::IntrepidManager::Bases "

these are the \"transformed_basis_values\" at each point in each cell
in the work set ([C],[ B],[P]), or ([C],[ B],[P],[D]) for GRAD here we
assume that [ B] is equivalent to [V]

C++ includes: IntrepidManager.hpp ";

%feature("docstring")  stk::percept::IntrepidManager::Bases::Bases "stk::percept::IntrepidManager::Bases::Bases(IM &im) ";


// File: classBases.xml
%feature("docstring") Bases "";

%feature("docstring")  Bases::Bases "Bases< T, I, Topology,
ReferenceElement >::Bases(T pt[Topology::dim]) ";

%feature("docstring")  Bases::basis "T Bases< T, I, Topology,
ReferenceElement >::basis(I i, T pt[Topology::dim]) ";


// File: classBasis.xml
%feature("docstring") Basis "";

%feature("docstring")  Basis::Basis "Basis< T, I, Topology,
ReferenceElement >::Basis() ";

%feature("docstring")  Basis::Basis "Basis< T, I, Topology,
ReferenceElement >::Basis(T pt[Topology::dim]) ";


// File: classBasis_3_01T_00_01I_00_01Quadrilateral4_00_01StandardReferenceElement_01_4.xml
%feature("docstring") Basis< T, I, Quadrilateral4,
StandardReferenceElement > " ";

%feature("docstring")  Basis< T, I, Quadrilateral4,
StandardReferenceElement >::Basis " Basis< T, I, Quadrilateral4,
StandardReferenceElement >::Basis(T pt[Base::topology::dim]) ";

%feature("docstring")  Basis< T, I, Quadrilateral4,
StandardReferenceElement >::basis " T Basis< T, I, Quadrilateral4,
StandardReferenceElement >::basis(I i, T pt[Base::topology::dim]) ";


// File: classIntrepid_1_1Basis__HGRAD__HEX__C2__Serendipity__FEM.xml
%feature("docstring") Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM "

Implementation of the default H(grad)-compatible FEM basis of
INCOMPLETE degree 2 on Hexahedron cell.

Implements Lagrangian basis of incomplete degree 2 on the reference
Hexahedron cell. The basis has cardinality 20 and spans an INCOMPLETE
tri-quadratic polynomial space. Basis functions are dual to a
unisolvent set of degrees-of-freedom (DoF) defined and enumerated as
follows:

=================================================================================================
|         |           degree-of-freedom-tag table                    |
|   |   DoF
|----------------------------------------------------------|      DoF
definition       |   | ordinal |  subc dim    | subc ordinal | subc
DoF ord |subc num DoF |                           |
|=========|==============|==============|==============|=============|===========================|
|    0    |       0      |       0      |       0      |      1      |
L_0(u) = u(-1,-1,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    1    |       0      |       1      |       0      |      1      |
L_1(u) = u( 1,-1,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    2    |       0      |       2      |       0      |      1      |
L_2(u) = u( 1, 1,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    3    |       0      |       3      |       0      |      1      |
L_3(u) = u(-1, 1,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    4    |       0      |       4      |       0      |      1      |
L_4(u) = u(-1,-1, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    5    |       0      |       5      |       0      |      1      |
L_5(u) = u( 1,-1, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    6    |       0      |       6      |       0      |      1      |
L_6(u) = u( 1, 1, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    7    |       0      |       7      |       0      |      1      |
L_7(u) = u(-1, 1, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    8    |       1      |       0      |       0      |      1      |
L_8(u) = u( 0,-1,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    9    |       1      |       1      |       0      |      1      |
L_9(u) = u( 1, 0,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   10    |       1      |       2      |       0      |      1      |
L_10(u) = u( 0, 1,-1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   11    |       1      |       3      |       0      |      1      |
L_11(u) = u(-1, 0,-1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   12    |       1      |       8      |       0      |      1      |
L_12(u) = u(-1,-1, 0)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   13    |       1      |       9      |       0      |      1      |
L_13(u) = u( 1,-1, 0)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   14    |       1      |      10      |       0      |      1      |
L_14(u) = u( 1, 1, 0)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   15    |       1      |      11      |       0      |      1      |
L_15(u) = u(-1, 1, 0)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   16    |       1      |       4      |       0      |      1      |
L_16(u) = u( 0,-1, 1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   17    |       1      |       5      |       0      |      1      |
L_17(u) = u( 1, 0, 1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   18    |       1      |       6      |       0      |      1      |
L_18(u) = u( 0, 1, 1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   19    |       1      |       7      |       0      |      1      |
L_19(u) = u(-1, 0, 1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|=========|==============|==============|==============|=============|===========================|
|   MAX   |  maxScDim=2  |  maxScOrd=12 |  maxDfOrd=0  |      -      |
|
|=========|==============|==============|==============|=============|===========================|

Ordering of DoFs follows the node order in Hexahedron<20> topology.
Note that node order in this topology does not follow the natural oder
of k-subcells where the nodes are located, except for nodes 0 to 7
which coincide with the vertices of the base Hexahedrn <8> topology.
As a result, L_0 to L_7 are associated with nodes 0 to 7, but L_8 to
L_19 are not associated with edges 0 to 12 in that order.

C++ includes: Intrepid_HGRAD_HEX_C2_Serendipity_FEM.hpp ";

%feature("docstring")
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::Basis_HGRAD_HEX_C2_Serendipity_FEM
"Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM< Scalar, ArrayScalar
>::Basis_HGRAD_HEX_C2_Serendipity_FEM()

Constructor. ";

%feature("docstring")
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::getValues "void
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getValues(ArrayScalar &outputValues, const ArrayScalar
&inputPoints, const EOperator operatorType) const

Evaluation of a FEM basis on a reference Hexahedron cell.

Returns values of operatorType acting on FEM basis functions for a set
of points in the reference Hexahedron cell. For rank and dimensions of
I/O array arguments see Section basis_md_array_sec.

Parameters:
-----------

outputValues:  [out] - rank-2 or 3 array with the computed basis
values

inputPoints:  [in] - rank-2 array with dimensions (P,D) containing
reference points

operatorType:  [in] - operator applied to basis functions ";

%feature("docstring")
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::getValues "void
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getValues(ArrayScalar &outputValues, const ArrayScalar
&inputPoints, const ArrayScalar &cellVertices, const EOperator
operatorType=OPERATOR_VALUE) const

FVD basis evaluation: invocation of this method throws an exception.
";

%feature("docstring")
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM::getDofCoords "void
Intrepid::Basis_HGRAD_HEX_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getDofCoords(ArrayScalar &DofCoords) const

Returns spatial locations (coordinates) of degrees of freedom on a
reference Quadrilateral.

Parameters:
-----------

DofCoords:  [out] - array with the coordinates of degrees of freedom,
dimensioned (F,D) ";


// File: classIntrepid_1_1Basis__HGRAD__QUAD__C2__Serendipity__FEM.xml
%feature("docstring") Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM "

Implementation of the default H(grad)-compatible FEM serendipity basis
of incomplete degree 2 on Quadrilateral cell.

Implements Lagrangian serendipity basis of degree 2 on the reference
Quadrilateral cell. The basis has cardinality 8 and spans an
INCOMPLETE bi-quadratic polynomial space. Basis functions are dual to
a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as
follows:

=================================================================================================
|         |           degree-of-freedom-tag table                    |
|   |   DoF
|----------------------------------------------------------|      DoF
definition       |   | ordinal |  subc dim    | subc ordinal | subc
DoF ord |subc num DoF |                           |
|=========|==============|==============|==============|=============|===========================|
|    0    |       0      |       0      |       0      |      1      |
L_0(u) = u(-1,-1)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    1    |       0      |       1      |       0      |      1      |
L_1(u) = u( 1,-1)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    2    |       0      |       2      |       0      |      1      |
L_2(u) = u( 1, 1)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    3    |       0      |       3      |       0      |      1      |
L_3(u) = u(-1, 1)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    4    |       1      |       0      |       0      |      1      |
L_4(u) = u( 0,-1)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    5    |       1      |       1      |       0      |      1      |
L_5(u) = u( 1, 0)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    6    |       1      |       2      |       0      |      1      |
L_6(u) = u( 0, 1)       |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    7    |       1      |       3      |       0      |      1      |
L_7(u) = u(-1, 0)       |
|=========|==============|==============|==============|=============|===========================|
|   MAX   |  maxScDim=2  |  maxScOrd=3  |  maxDfOrd=0  |     -       |
|
|=========|==============|==============|==============|=============|===========================|

C++ includes: Intrepid_HGRAD_QUAD_C2_Serendipity_FEM.hpp ";

%feature("docstring")
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::Basis_HGRAD_QUAD_C2_Serendipity_FEM
"Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM< Scalar, ArrayScalar
>::Basis_HGRAD_QUAD_C2_Serendipity_FEM()

Constructor. ";

%feature("docstring")
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::getValues "void
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getValues(ArrayScalar &outputValues, const ArrayScalar
&inputPoints, const EOperator operatorType) const

FEM basis evaluation on a reference Quadrilateral cell.

Returns values of operatorType acting on FEM basis functions for a set
of points in the reference Quadrilateral cell. For rank and dimensions
of I/O array arguments see Section basis_md_array_sec .

Parameters:
-----------

outputValues:  [out] - rank-2 or 3 array with the computed basis
values

inputPoints:  [in] - rank-2 array with dimensions (P,D) containing
reference points

operatorType:  [in] - operator applied to basis functions ";

%feature("docstring")
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::getValues "void
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getValues(ArrayScalar &outputValues, const ArrayScalar
&inputPoints, const ArrayScalar &cellVertices, const EOperator
operatorType=OPERATOR_VALUE) const

FVD basis evaluation: invocation of this method throws an exception.
";

%feature("docstring")
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM::getDofCoords "void
Intrepid::Basis_HGRAD_QUAD_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getDofCoords(ArrayScalar &DofCoords) const

Returns spatial locations (coordinates) of degrees of freedom on a
reference Quadrilateral.

Parameters:
-----------

DofCoords:  [out] - array with the coordinates of degrees of freedom,
dimensioned (F,D) ";


// File: classIntrepid_1_1Basis__HGRAD__WEDGE__C2__Serendipity__FEM.xml
%feature("docstring") Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM "

Implementation of the default H(grad)-compatible FEM basis of
incomplete degree 2 on Wedge cell.

Implements Lagrangian basis of degree 2 on the reference Wedge cell.
The basis has cardinality 15 and spans a INCOMPLETE bi-quadratic
polynomial space. Basis functions are dual to a unisolvent set of
degrees-of-freedom (DoF) defined and enumerated as follows:

=================================================================================================
|         |           degree-of-freedom-tag table                    |
|   |   DoF
|----------------------------------------------------------|      DoF
definition       |   | ordinal |  subc dim    | subc ordinal | subc
DoF ord |subc num DoF |                           |
|=========|==============|==============|==============|=============|===========================|
|    0    |       0      |       0      |       0      |      1      |
L_0(u) = u( 0, 0,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    1    |       0      |       1      |       0      |      1      |
L_1(u) = u( 1, 0,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    2    |       0      |       2      |       0      |      1      |
L_2(u) = u( 0, 1,-1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    3    |       0      |       3      |       0      |      1      |
L_3(u) = u( 0, 0, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    4    |       0      |       4      |       0      |      1      |
L_4(u) = u( 1, 0, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    5    |       0      |       5      |       0      |      1      |
L_5(u) = u( 0, 1, 1)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    6    |       1      |       0      |       0      |      1      |
L_6(u) = u(1/2, 0,-1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    7    |       1      |       1      |       0      |      1      |
L_7(u) = u(1/2,1/2,-1)  |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    8    |       1      |       2      |       0      |      1      |
L_8(u) = u( 0,1/2,-1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|    9    |       1      |       6      |       0      |      1      |
L_9(u) = u( 0, 0, 0)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   10    |       1      |       7      |       0      |      1      |
L_10(u)= u( 1, 0, 0)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   11    |       1      |       8      |       0      |      1      |
L_11(u)= u( 0, 1, 0)    |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   12    |       1      |       3      |       0      |      1      |
L_12(u)= u(1/2, 0, 1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   13    |       1      |       4      |       0      |      1      |
L_13(u)= u(1/2,1/2, 1)  |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|   14    |       1      |       5      |       0      |      1      |
L_14(u)= u( 0,1/2, 1)   |
|---------|--------------|--------------|--------------|-------------|---------------------------|
|=========|==============|==============|==============|=============|===========================|
|   MAX   |  maxScDim=2  |  maxScOrd=8  |  maxDfOrd=0  |      -      |
|
|=========|==============|==============|==============|=============|===========================|

Ordering of DoFs follows the node order in Wedge<15> topology. Note
that node order in this topology does not follow the natural oder of
k-subcells where the nodes are located, except for nodes 0 to 5 which
coincide with the vertices of the base Wedge<6> topology. As a result,
L_0 to L_5 are associated with nodes 0 to 5, but L_6 to L_14 are not
associated with edges 0 to 9 in that order.

C++ includes: Intrepid_HGRAD_WEDGE_C2_Serendipity_FEM.hpp ";

%feature("docstring")
Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM::Basis_HGRAD_WEDGE_C2_Serendipity_FEM
"Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM< Scalar, ArrayScalar
>::Basis_HGRAD_WEDGE_C2_Serendipity_FEM()

Constructor. ";

%feature("docstring")
Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM::getValues "void
Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getValues(ArrayScalar &outputValues, const ArrayScalar
&inputPoints, const EOperator operatorType) const

FEM basis evaluation on a reference Wedge cell.

Returns values of operatorType acting on FEM basis functions for a set
of points in the reference Wedge cell. For rank and dimensions of I/O
array arguments see Section basis_md_array_sec .

Parameters:
-----------

outputValues:  [out] - rank-2 or 3 array with the computed basis
values

inputPoints:  [in] - rank-2 array with dimensions (P,D) containing
reference points

operatorType:  [in] - operator applied to basis functions ";

%feature("docstring")
Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM::getValues "void
Intrepid::Basis_HGRAD_WEDGE_C2_Serendipity_FEM< Scalar, ArrayScalar
>::getValues(ArrayScalar &outputValues, const ArrayScalar
&inputPoints, const ArrayScalar &cellVertices, const EOperator
operatorType=OPERATOR_VALUE) const

FVD basis evaluation: invocation of this method throws an exception.
";


// File: classBasisBase.xml
%feature("docstring") BasisBase "";

%feature("docstring")  BasisBase::BasisBase "BasisBase< T, I,
Topology, ReferenceElement >::BasisBase(T pt[Topology::dim]) ";

%feature("docstring")  BasisBase::basis "T BasisBase< T, I, Topology,
ReferenceElement >::basis(I i, T pt[Topology::dim]) ";


// File: classBasisImpl.xml
%feature("docstring") BasisImpl "";


// File: classBasisImpl_3_01T_00_01I_00_01Quadrilateral4_00_01StandardReferenceElement_01_4.xml
%feature("docstring") BasisImpl< T, I, Quadrilateral4,
StandardReferenceElement > " ";

%feature("docstring")  BasisImpl< T, I, Quadrilateral4,
StandardReferenceElement >::basis " T BasisImpl< T, I, Quadrilateral4,
StandardReferenceElement >::basis(I i, T pt[Quadrilateral4::dim]) ";


// File: classstk_1_1percept_1_1BeamFixture.xml
%feature("docstring") stk::percept::BeamFixture "

Use case with mixed element topologies and field relations to provide
fast access to node field data from an element.

copied from stk_mesh and modified

C++ includes: BeamFixture.hpp ";

%feature("docstring")  stk::percept::BeamFixture::~BeamFixture "stk::percept::BeamFixture::~BeamFixture() ";

%feature("docstring")  stk::percept::BeamFixture::BeamFixture "stk::percept::BeamFixture::BeamFixture(stk::ParallelMachine comm, bool
doCommit=true) ";

%feature("docstring")  stk::percept::BeamFixture::populate "void
stk::percept::BeamFixture::populate() ";


// File: classstk_1_1percept_1_1BucketOp.xml
%feature("docstring") stk::percept::BucketOp "";

%feature("docstring")  stk::percept::BucketOp::~BucketOp "virtual
stk::percept::BucketOp::~BucketOp() ";


// File: classstk_1_1percept_1_1BuildBoundingBoxes.xml
%feature("docstring") stk::percept::BuildBoundingBoxes "";

%feature("docstring")
stk::percept::BuildBoundingBoxes::BuildBoundingBoxes "stk::percept::BuildBoundingBoxes< SpatialDim
>::BuildBoundingBoxes(std::vector< AABoundingBox > &boxes,
VectorFieldType *coords_field) ";

%feature("docstring")  stk::percept::BuildBoundingBoxes::init "void
stk::percept::BuildBoundingBoxes< SpatialDim >::init(std::vector<
AABoundingBox > &boxes, VectorFieldType *coords_field) ";

%feature("docstring")
stk::percept::BuildBoundingBoxes::init_elementOp "void
stk::percept::BuildBoundingBoxes< SpatialDim >::init_elementOp() ";

%feature("docstring")
stk::percept::BuildBoundingBoxes::fini_elementOp "void
stk::percept::BuildBoundingBoxes< SpatialDim >::fini_elementOp() ";

%feature("docstring")
stk::percept::BuildBoundingBoxes::getBoundingBox "AABoundingBox
stk::percept::BuildBoundingBoxes< SpatialDim >::getBoundingBox(const
stk::mesh::Entity &element, const mesh::BulkData &bulkData) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1CellWorkSet.xml
%feature("docstring") stk::percept::IntrepidManager::CellWorkSet "

([C], [V], [D])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::CellWorkSet::CellWorkSet "stk::percept::IntrepidManager::CellWorkSet::CellWorkSet(IM &im)

([C], [V], [D]) ";


// File: classstk_1_1percept_1_1unit__tests_1_1CheckCoordMag.xml
%feature("docstring") stk::percept::unit_tests::CheckCoordMag "";

%feature("docstring")
stk::percept::unit_tests::CheckCoordMag::CheckCoordMag "stk::percept::unit_tests::CheckCoordMag::CheckCoordMag(std::string
name=\"\") ";


// File: classCircle.xml
%feature("docstring") Circle "";

%feature("docstring")  Circle::Circle "Circle::Circle(double r) ";

%feature("docstring")  Circle::area "virtual double
Circle::area(void) ";

%feature("docstring")  Circle::perimeter "virtual double
Circle::perimeter(void) ";


// File: classstk_1_1percept_1_1CodeVerifier.xml
%feature("docstring") stk::percept::CodeVerifier "";


// File: classstk_1_1adapt_1_1Colorer.xml
%feature("docstring") stk::adapt::Colorer "";

%feature("docstring")  stk::adapt::Colorer::Colorer "stk::adapt::Colorer::Colorer(std::vector< ColorerSetType >
&element_colors, std::vector< mesh::EntityRank > ranks) ";

%feature("docstring")  stk::adapt::Colorer::Colorer "stk::adapt::Colorer::Colorer(std::vector< mesh::EntityRank > ranks) ";

%feature("docstring")  stk::adapt::Colorer::color "void
stk::adapt::Colorer::color(percept::PerceptMesh &eMesh, unsigned
*elementType=0, mesh::PartVector *fromParts=0, mesh::FieldBase
*element_color_field=0) ";

%feature("docstring")  stk::adapt::Colorer::setNoColoring "void
stk::adapt::Colorer::setNoColoring(bool no_coloring)

Set to true to avoid the coloring step and just return elements in a
single color. ";

%feature("docstring")  stk::adapt::Colorer::getNoColoring "bool
stk::adapt::Colorer::getNoColoring() ";

%feature("docstring")  stk::adapt::Colorer::getElementColors "std::vector< ColorerSetType > &
stk::adapt::Colorer::getElementColors() ";


// File: structstk_1_1adapt_1_1CompareSDSEntityType.xml
%feature("docstring") stk::adapt::CompareSDSEntityType "";


// File: classstk_1_1percept_1_1CompositeFunction.xml
%feature("docstring") stk::percept::CompositeFunction "";

%feature("docstring")
stk::percept::CompositeFunction::CompositeFunction "stk::percept::CompositeFunction::CompositeFunction(const char *name,
Function &func_1, Function &func_2, Dimensions
domain_dimensions=Dimensions(), Dimensions
codomain_dimensions=Dimensions(), unsigned integration_order=0)

compose two functions to be able to apply in turn as in
func_2(func_1(x)), or more specifically as: func_1(domain,
codomain_temp); Note that since this is a Function also, one can make
multiple compositions e.g. h(g(f(x))) by CompositeFunction g_of_f (f,
g) CompositeFunction h_of_g_of_f (g_of_f, h); The first function in
the list is always applied first. ";


// File: classstk_1_1percept_1_1util_1_1CompositeGeneralFunction.xml
%feature("docstring") stk::percept::util::CompositeGeneralFunction "";

%feature("docstring")
stk::percept::util::CompositeGeneralFunction::CompositeGeneralFunction
"stk::percept::util::CompositeGeneralFunction< domain_f,
codomain_f_and_domain_g, codomain_g
>::CompositeGeneralFunction(GeneralFunction< codomain_f_and_domain_g,
codomain_g > &g, GeneralFunction< domain_f, codomain_f_and_domain_g >
&f, codomain_f_and_domain_g &tmp_fv=0) ";


// File: classstk_1_1percept_1_1ComputeBases.xml
%feature("docstring") stk::percept::ComputeBases "";

%feature("docstring")  stk::percept::ComputeBases::ComputeBases "stk::percept::ComputeBases::ComputeBases() ";

%feature("docstring")  stk::percept::ComputeBases::getBases "void
stk::percept::ComputeBases::getBases(const stk::mesh::Bucket &bucket,
const MDArray &parametric_coordinates, MDArray
&transformed_basis_values, int which_cell=-1)

([P],[D]) ";

%feature("docstring")  stk::percept::ComputeBases::getBases "void
stk::percept::ComputeBases::getBases(const stk::mesh::Entity &element,
const MDArray &parametric_coordinates, MDArray
&transformed_basis_values) ";


// File: classstk_1_1percept_1_1ComputeFieldValues.xml
%feature("docstring") stk::percept::ComputeFieldValues "";

%feature("docstring")
stk::percept::ComputeFieldValues::get_fieldValues "void
stk::percept::ComputeFieldValues::get_fieldValues(const
stk::mesh::Entity &element, MDArray &transformed_basis_values,
mesh::FieldBase *field, MDArray &output_field_values)

NOTE: this is needed since FunctionSpaceTools::evaluate method assumes
the output array is initialized to 0 ";


// File: classstk_1_1percept_1_1ConstantFunction.xml
%feature("docstring") stk::percept::ConstantFunction "";

%feature("docstring")
stk::percept::ConstantFunction::ConstantFunction "stk::percept::ConstantFunction::ConstantFunction(double value, const
char *name, Dimensions domain_dimensions=Dimensions(), Dimensions
codomain_dimensions=Dimensions(), unsigned integration_order=0) ";

%feature("docstring")  stk::percept::ConstantFunction::getValue "double& stk::percept::ConstantFunction::getValue() ";

%feature("docstring")  stk::percept::ConstantFunction::setValue "void
stk::percept::ConstantFunction::setValue(double &v) ";


// File: classstk_1_1percept_1_1ConstantFunctionVec.xml
%feature("docstring") stk::percept::ConstantFunctionVec "";

%feature("docstring")
stk::percept::ConstantFunctionVec::ConstantFunctionVec "stk::percept::ConstantFunctionVec::ConstantFunctionVec(std::vector<
double > &value, const char *name, Dimensions
domain_dimensions=Dimensions(), Dimensions
codomain_dimensions=Dimensions(), unsigned integration_order=0) ";

%feature("docstring")  stk::percept::ConstantFunctionVec::getValue "std::vector<double>& stk::percept::ConstantFunctionVec::getValue() ";

%feature("docstring")  stk::percept::ConstantFunctionVec::setValue "void stk::percept::ConstantFunctionVec::setValue(std::vector< double >
&v) ";


// File: classCPPArray.xml
%feature("docstring") CPPArray "";

%feature("docstring")  CPPArray::dimension "const unsigned CPPArray<
D0, D1, D2, T >::dimension(const unsigned d) const ";

%feature("docstring")  CPPArray::rank "const unsigned CPPArray< D0,
D1, D2, T >::rank() const ";

%feature("docstring")  CPPArray::dimension "const unsigned CPPArray<
D0, D1, D2, T >::dimension(const unsigned d) ";

%feature("docstring")  CPPArray::rank "const unsigned CPPArray< D0,
D1, D2, T >::rank() const ";

%feature("docstring")  CPPArray::CPPArray "CPPArray< D0, D1, D2, T
>::CPPArray() ";

%feature("docstring")  CPPArray::init "void CPPArray< D0, D1, D2, T
>::init() ";

%feature("docstring")  CPPArray::dimension "const unsigned CPPArray<
D0, D1, D2, T >::dimension(const unsigned d) const ";

%feature("docstring")  CPPArray::rank "const unsigned CPPArray< D0,
D1, D2, T >::rank() const ";

%feature("docstring")  CPPArray::CPPArray "CPPArray< D0, D1, D2, T
>::CPPArray() ";

%feature("docstring")  CPPArray::dimension "const unsigned CPPArray<
D0, D1, D2, T >::dimension(const unsigned d) const ";

%feature("docstring")  CPPArray::rank "const unsigned CPPArray< D0,
D1, D2, T >::rank() const ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1CubaturePoints.xml
%feature("docstring") stk::percept::IntrepidManager::CubaturePoints "

([P],[D])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::CubaturePoints::CubaturePoints "stk::percept::IntrepidManager::CubaturePoints::CubaturePoints(IM &im)

([P],[D]) ";

%feature("docstring")
stk::percept::IntrepidManager::CubaturePoints::copyTo "void
stk::percept::IntrepidManager::CubaturePoints::copyTo(MDArray &mda) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1CubatureWeights.xml
%feature("docstring") stk::percept::IntrepidManager::CubatureWeights "

([P])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::CubatureWeights::CubatureWeights "stk::percept::IntrepidManager::CubatureWeights::CubatureWeights(IM
&im)

([P]) ";


// File: structDim.xml
%feature("docstring") Dim "";


// File: classstk_1_1percept_1_1Dimensions.xml
%feature("docstring") stk::percept::Dimensions "";

%feature("docstring")  stk::percept::Dimensions::Dimensions "stk::percept::Dimensions::Dimensions() ";

%feature("docstring")  stk::percept::Dimensions::Dimensions "stk::percept::Dimensions::Dimensions(int i0) ";

%feature("docstring")  stk::percept::Dimensions::Dimensions "stk::percept::Dimensions::Dimensions(int i0, int i1) ";

%feature("docstring")  stk::percept::Dimensions::Dimensions "stk::percept::Dimensions::Dimensions(int i0, int i1, int i2) ";


// File: classElementFunctor.xml
%feature("docstring") ElementFunctor "";


// File: classstk_1_1percept_1_1ElementOp.xml
%feature("docstring") stk::percept::ElementOp "";

%feature("docstring")  stk::percept::ElementOp::init_elementOp "virtual void stk::percept::ElementOp::init_elementOp()=0 ";

%feature("docstring")  stk::percept::ElementOp::fini_elementOp "virtual void stk::percept::ElementOp::fini_elementOp()=0 ";

%feature("docstring")  stk::percept::ElementOp::~ElementOp "virtual
stk::percept::ElementOp::~ElementOp() ";


// File: structstk_1_1adapt_1_1ElementRefinePredicate.xml
%feature("docstring") stk::adapt::ElementRefinePredicate "";

%feature("docstring")
stk::adapt::ElementRefinePredicate::ElementRefinePredicate "stk::adapt::ElementRefinePredicate::ElementRefinePredicate(stk::mesh::Selector
*selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) ";


// File: structstk_1_1percept_1_1interface__table_1_1elemInfoType.xml
%feature("docstring") stk::percept::interface_table::elemInfoType "";


// File: classstk_1_1percept_1_1Example2FunctionWithIntrepidRequest.xml
%feature("docstring")
stk::percept::Example2FunctionWithIntrepidRequest "";

%feature("docstring")
stk::percept::Example2FunctionWithIntrepidRequest::Example2FunctionWithIntrepidRequest
"stk::percept::Example2FunctionWithIntrepidRequest::Example2FunctionWithIntrepidRequest(BulkData
&bulkData) ";


// File: classstk_1_1percept_1_1ExampleFunctionWithIntrepidRequest.xml
%feature("docstring") stk::percept::ExampleFunctionWithIntrepidRequest
"";

%feature("docstring")
stk::percept::ExampleFunctionWithIntrepidRequest::ExampleFunctionWithIntrepidRequest
"stk::percept::ExampleFunctionWithIntrepidRequest::ExampleFunctionWithIntrepidRequest(BulkData
&bulkData) ";


// File: classstk_1_1percept_1_1ExceptionWatch.xml
%feature("docstring") stk::percept::ExceptionWatch "";

%feature("docstring")  stk::percept::ExceptionWatch::ExceptionWatch "stk::percept::ExceptionWatch::ExceptionWatch(int line, char const
*pfname) ";

%feature("docstring")  stk::percept::ExceptionWatch::~ExceptionWatch "stk::percept::ExceptionWatch::~ExceptionWatch() ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1FaceNormal.xml
%feature("docstring") stk::percept::IntrepidManager::FaceNormal "

([C], [P], [D])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::FaceNormal::FaceNormal "stk::percept::IntrepidManager::FaceNormal::FaceNormal(IM &im)

([C], [P], [D]) ";


// File: structstk_1_1percept_1_1FieldCreateOrder.xml
%feature("docstring") stk::percept::FieldCreateOrder "";

%feature("docstring")
stk::percept::FieldCreateOrder::FieldCreateOrder "stk::percept::FieldCreateOrder::FieldCreateOrder() ";

%feature("docstring")
stk::percept::FieldCreateOrder::FieldCreateOrder "stk::percept::FieldCreateOrder::FieldCreateOrder(const std::string
name, const unsigned entity_rank, const std::vector< int > dimensions,
const mesh::Part *part) ";


// File: classstk_1_1percept_1_1FieldFunction.xml
%feature("docstring") stk::percept::FieldFunction "

Evaluate the function at this input point (or points) returning
value(s) in output_field_values

In the following, the arrays are dimensioned using the notation (from
Intrepid's doc):

[C] - num. integration domains (cells/elements) [F] - num. Intrepid
\"fields\" (number of bases within an element == num. nodes typically)
[P] - num. integration (or interpolation) points within the element
[D] - spatial dimension [D1], [D2] - spatial dimension

Locally, we introduce this notation:

[DOF] - number of degrees-of-freedom per node of the interpolated stk
Field. For example, a vector field in 3D has [DOF] = 3

C++ includes: FieldFunction.hpp ";

%feature("docstring")  stk::percept::FieldFunction::FieldFunction "stk::percept::FieldFunction::FieldFunction(const char *name,
mesh::FieldBase *field, PerceptMesh &mesh, int domain_dimension=3, int
codomain_dimension=1, SearchType searchType=SIMPLE_SEARCH, unsigned
integration_order=0) ";

%feature("docstring")  stk::percept::FieldFunction::get_field "mesh::FieldBase * stk::percept::FieldFunction::get_field() ";

%feature("docstring")  stk::percept::FieldFunction::interpolateFrom "void stk::percept::FieldFunction::interpolateFrom(Function &function)

for each node in the codomain, evaluate the function_to_interpolate's
function, assign to the codomain field ";

%feature("docstring")  stk::percept::FieldFunction::FieldFunction "stk::percept::FieldFunction::FieldFunction(const char *name,
mesh::FieldBase *field, mesh::BulkData *bulk, Dimensions
domain_dimensions=Dimensions(), Dimensions
codomain_dimensions=Dimensions(), SearchType searchType=SIMPLE_SEARCH,
unsigned integration_order=0) ";

%feature("docstring")  stk::percept::FieldFunction::FieldFunction "stk::percept::FieldFunction::FieldFunction(const char *name,
mesh::FieldBase *field, PerceptMesh &eMesh, Dimensions
domain_dimensions, Dimensions codomain_dimensions, SearchType
searchType=SIMPLE_SEARCH, unsigned integration_order=0) ";

%feature("docstring")  stk::percept::FieldFunction::~FieldFunction "stk::percept::FieldFunction::~FieldFunction() ";

%feature("docstring")  stk::percept::FieldFunction::derivative "virtual Teuchos::RCP<Function >
stk::percept::FieldFunction::derivative(MDArrayString &deriv_spec)

Return a function that is the derivative of this function. The
derivative is specified as a rank-2 array of strings that specify what
derivative to take and how many derivatives. For example, ";

%feature("docstring")  stk::percept::FieldFunction::gradient "virtual
Teuchos::RCP<Function > stk::percept::FieldFunction::gradient(int
spatialDim=3) ";

%feature("docstring")  stk::percept::FieldFunction::localEvaluation "void stk::percept::FieldFunction::localEvaluation(MDArray &in, MDArray
&out, double time_value_optional=0.0) ";

%feature("docstring")  stk::percept::FieldFunction::setup_searcher "void stk::percept::FieldFunction::setup_searcher(int D_) ";

%feature("docstring")  stk::percept::FieldFunction::helper "void
stk::percept::FieldFunction::helper(MDArray &input_phy_points, MDArray
&output_field_values, const BucketOrEntity &bucket_or_element, const
MDArray &parametric_coordinates, double time_value_optional)

NOTE: this is needed since FunctionSpaceTools::evaluate method assumes
the output array is initialized to 0 ";

%feature("docstring")  stk::percept::FieldFunction::get_bulk_data "mesh::BulkData * stk::percept::FieldFunction::get_bulk_data() ";

%feature("docstring")
stk::percept::FieldFunction::getFoundOnLocalOwnedPart "bool
stk::percept::FieldFunction::getFoundOnLocalOwnedPart() ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1FieldValues.xml
%feature("docstring") stk::percept::IntrepidManager::FieldValues "

([C],[P],[DOF]): evaluated field values at each integration point in
each cell:

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::FieldValues::FieldValues "stk::percept::IntrepidManager::FieldValues::FieldValues(IM &im)

([C],[P],[DOF]): evaluated field values at each integration point in
each cell: ";


// File: structstk_1_1percept_1_1unit__tests_1_1FindMapItem1.xml
%feature("docstring") stk::percept::unit_tests::FindMapItem1 "";


// File: structstk_1_1percept_1_1unit__tests_1_1FindMapItem2.xml
%feature("docstring") stk::percept::unit_tests::FindMapItem2 "";


// File: classFlux__.xml
%feature("docstring") Flux_ "";

%feature("docstring")  Flux_::Flux_ "Flux_< Mu, Mu2, Temperature
>::Flux_(Mu mu, Mu2 mu2, Temperature temperature) ";

%feature("docstring")  Flux_::Flux_ "Flux_< Mu, Mu2, Temperature
>::Flux_(Mu mu, Mu2 mu2, Temperature temperature) ";


// File: classstk_1_1percept_1_1Function.xml
%feature("docstring") stk::percept::Function "";

%feature("docstring")  stk::percept::Function::Function "stk::percept::Function::Function()

Create a function with the given name, domain dimensions and codomain
dimensions, and integration order. If domain_dimensions and
codomain_dimensions are defaulted, then they are dimensioned using the
defaults of domain_dimensions = {s_spatialDimDefault} and
codomain_dimensions = {s_codomainDimDefault} If integration_order is
not specified, it defaults to s_integration_order_default. The
defaults of the defaults is s_spatialDimDefault=3,
s_codomainDimDefault=1, and s_integration_order_default=1

[DEPRECATED START] NOTE: we assume that input arrays have either
{x,y,t} 2D+time or {x,y,z,t} for 3D+time. There is no separate
argument for time. This means arrays should be dimensioned with length
3 (2D+time), or 4 (3D+time) respectively. [DEPRECATED END]: now we
have time in the operator() arg list

A Function has a { core} dimensioning of its domain and codamain, as
specified in its constructor (or by setting with accessors
setDomainDimensions() and setCodomainDimensions(). This means that all
Function's implementations of operator()(MDArray& input, MDArray&
output) expect the rightmost dimensions of input to be equal to the
rightmost dimensions of its domain dimensions, and similar for output
and codomain_dimensions

Conventions: 1. core dimensions given by getDomainDimensions() and
getCodomainDimensions() properties 2. rightmost dimensions of input
array must match rightmost dimensions of Function's domain 3.
rightmost dimensions of output array must match rightmost dimensions
of Function's codomain 4. domain and codomain can have different ranks
and dimensions 5. usage of conventions 1-4 are in scenarios where the
input and output arrays contain multiple points that the client is
requesting the Function to evaluate. For example:

int numElements = bucket.size(); int numQuadPoints =
quadratureFactory.getNumPoints(); int spaceDim = 3; int numValuesAtQP
= 1;

MDArray quadrature_points(numElements, numQuadPoints, spaceDim);
MDArray function_value_at_quadrature_points(numElements,
numQuadPoints, numValuesAtQP);

StringFunction sf(\"<... some definition ...>\", \"<... some name
...>\", Dimensions(spaceDim), Dimensions(numValuesAtQP) );

MDArray single_point(spaceDim); MDArray
value_at_single_point(numValuesAtQP);

// note that this same sf can be evaluated on either of these pairs of
in/out arrays string_function(quadrature_points,
function_value_at_quadrature_points); string_function(single_point,
value_at_single_point); ";

%feature("docstring")  stk::percept::Function::Function "stk::percept::Function::Function(const char *name, Dimensions
domain_dimensions=Dimensions(), Dimensions
codomain_dimensions=Dimensions(), unsigned integration_order=0) ";

%feature("docstring")  stk::percept::Function::value "void
stk::percept::Function::value(MDArray &domain, MDArray &MDOutVal,
double time_value_optional=0.0)

this version uses the MDOutVal argument as a predefined array to
ensure the python returned value is properly sized ";

%feature("docstring")  stk::percept::Function::value "MDArray
stk::percept::Function::value(MDArray &domain, double
time_value_optional=0.0)

this version creates a dummy output and returns it based on the
codomain-dimensions specified at construction ";

%feature("docstring")  stk::percept::Function::derivative "virtual
Teuchos::RCP<Function >
stk::percept::Function::derivative(MDArrayString &deriv_spec)

Return a function that is the derivative of this function. The
derivative is specified as a rank-2 array of strings that specify what
derivative to take and how many derivatives. For example, ";

%feature("docstring")  stk::percept::Function::gradient "virtual
Teuchos::RCP<Function > stk::percept::Function::gradient(int
spatialDim=3) ";

%feature("docstring")  stk::percept::Function::derivativeAtPoint "void stk::percept::Function::derivativeAtPoint(MDArrayString
&deriv_spec, MDArray &domain, MDArray &codomain, double time=0.0) ";

%feature("docstring")  stk::percept::Function::add_alias "Function *
stk::percept::Function::add_alias(const char *alias)

allow this function to have one or more aliases ";

%feature("docstring")  stk::percept::Function::getName "std::string&
stk::percept::Function::getName() ";

%feature("docstring")  stk::percept::Function::setIntegrationOrder "void stk::percept::Function::setIntegrationOrder(unsigned iord) ";

%feature("docstring")  stk::percept::Function::getIntegrationOrder "unsigned stk::percept::Function::getIntegrationOrder(void) ";

%feature("docstring")  stk::percept::Function::setDomainDimensions "void stk::percept::Function::setDomainDimensions(const Dimensions
dims) ";

%feature("docstring")  stk::percept::Function::setCodomainDimensions "void stk::percept::Function::setCodomainDimensions(const Dimensions
dims) ";

%feature("docstring")  stk::percept::Function::argsAreValid "bool
stk::percept::Function::argsAreValid(const MDArray &in, const MDArray
&out)

Verify that the last dimensions of

Parameters:
-----------

in:  and

out:  are the same; this allows Functions to be invoked at multiple
points where the first M indices represent an M-d array of points to
evaluate the function, while the last N indices should match the
Functions domain and codomain dimensions. ";


// File: classstk_1_1percept_1_1FunctionOperator.xml
%feature("docstring") stk::percept::FunctionOperator "";

%feature("docstring")
stk::percept::FunctionOperator::FunctionOperator "stk::percept::FunctionOperator::FunctionOperator(mesh::BulkData
&bulkData, mesh::Part *part=0) ";

%feature("docstring")
stk::percept::FunctionOperator::FunctionOperator "stk::percept::FunctionOperator::FunctionOperator(mesh::BulkData
&bulkData, mesh::Selector *selector) ";

%feature("docstring")  stk::percept::FunctionOperator::init "void
stk::percept::FunctionOperator::init(mesh::Part *part) ";

%feature("docstring")  stk::percept::FunctionOperator::init "void
stk::percept::FunctionOperator::init(mesh::Selector *selector) ";

%feature("docstring")
stk::percept::FunctionOperator::~FunctionOperator "virtual
stk::percept::FunctionOperator::~FunctionOperator() ";


// File: classstk_1_1percept_1_1FunctionWithIntrepidRequest.xml
%feature("docstring") stk::percept::FunctionWithIntrepidRequest "";

%feature("docstring")
stk::percept::FunctionWithIntrepidRequest::getFunction "Function*
stk::percept::FunctionWithIntrepidRequest::getFunction() ";

%feature("docstring")
stk::percept::FunctionWithIntrepidRequest::FunctionWithIntrepidRequest
"stk::percept::FunctionWithIntrepidRequest::FunctionWithIntrepidRequest()
";

%feature("docstring")
stk::percept::FunctionWithIntrepidRequest::FunctionWithIntrepidRequest
"stk::percept::FunctionWithIntrepidRequest::FunctionWithIntrepidRequest(Function
*func, Request values=Request(), Request gradient=Request(), Request
higherDerivs=Request()) ";


// File: classstk_1_1percept_1_1util_1_1GeneralFunction.xml
%feature("docstring") stk::percept::util::GeneralFunction "";

%feature("docstring")
stk::percept::util::GeneralFunction::GeneralFunction "stk::percept::util::GeneralFunction< domain, codomain
>::GeneralFunction() ";


// File: classstk_1_1percept_1_1util_1_1GeneralFunctionWithGrad.xml
%feature("docstring") stk::percept::util::GeneralFunctionWithGrad "";


// File: classstk_1_1percept_1_1GenericFunction.xml
%feature("docstring") stk::percept::GenericFunction "";

%feature("docstring")  stk::percept::GenericFunction::GenericFunction
"stk::percept::GenericFunction::GenericFunction(Dimensions
domain_dimensions=Dimensions(), Dimensions
codomain_dimensions=Dimensions()) ";

%feature("docstring")  stk::percept::GenericFunction::~GenericFunction
"virtual stk::percept::GenericFunction::~GenericFunction() ";

%feature("docstring")
stk::percept::GenericFunction::isSpatialOperator "bool
stk::percept::GenericFunction::isSpatialOperator() ";

%feature("docstring")
stk::percept::GenericFunction::setIsSpatialOperator "void
stk::percept::GenericFunction::setIsSpatialOperator(bool so) ";

%feature("docstring")
stk::percept::GenericFunction::getDomainDimensions "Dimensions
stk::percept::GenericFunction::getDomainDimensions() ";

%feature("docstring")
stk::percept::GenericFunction::getCodomainDimensions "Dimensions
stk::percept::GenericFunction::getCodomainDimensions() ";

%feature("docstring")  stk::percept::GenericFunction::getNewDomain "MDArray stk::percept::GenericFunction::getNewDomain() ";

%feature("docstring")  stk::percept::GenericFunction::getNewCodomain "MDArray stk::percept::GenericFunction::getNewCodomain() ";

%feature("docstring")  stk::percept::GenericFunction::getNewCodomain "MDArray stk::percept::GenericFunction::getNewCodomain() const ";


// File: classstk_1_1percept_1_1GenericVectorOfObjectPointers.xml
%feature("docstring") stk::percept::GenericVectorOfObjectPointers "";

%feature("docstring")
stk::percept::GenericVectorOfObjectPointers::GenericVectorOfObjectPointers
"stk::percept::GenericVectorOfObjectPointers< VecType
>::GenericVectorOfObjectPointers(int n) ";

%feature("docstring")
stk::percept::GenericVectorOfObjectPointers::GenericVectorOfObjectPointers
"stk::percept::GenericVectorOfObjectPointers< VecType
>::GenericVectorOfObjectPointers(VecType *vt1=0, VecType *vt2=0,
VecType *vt3=0, VecType *vt4=0, VecType *vt5=0, VecType *vt6=0,
VecType *vt7=0, VecType *vt8=0) ";


// File: structGeometryEvaluator.xml
%feature("docstring") GeometryEvaluator "";

%feature("docstring")  GeometryEvaluator::GeometryEvaluator "GeometryEvaluator::GeometryEvaluator(Part *part) ";


// File: classGeometryFactory.xml
%feature("docstring") GeometryFactory "";

%feature("docstring")  GeometryFactory::GeometryFactory "GeometryFactory::GeometryFactory(GeometryKernel *kernel, MeshGeometry
*geometry) ";

%feature("docstring")  GeometryFactory::~GeometryFactory "GeometryFactory::~GeometryFactory() ";

%feature("docstring")  GeometryFactory::read_file "bool
GeometryFactory::read_file(const std::string &filename, PerceptMesh
*mesh) ";


// File: classGeometryKernel.xml
%feature("docstring") GeometryKernel "";

%feature("docstring")  GeometryKernel::GeometryKernel "GeometryKernel::GeometryKernel() ";

%feature("docstring")  GeometryKernel::~GeometryKernel "virtual
GeometryKernel::~GeometryKernel() ";

%feature("docstring")  GeometryKernel::read_file "virtual bool
GeometryKernel::read_file(const std::string &file_name, std::vector<
GeometryHandle > &geometry_entities)=0 ";

%feature("docstring")  GeometryKernel::get_attribute "virtual
std::string GeometryKernel::get_attribute(GeometryHandle geom)=0 ";

%feature("docstring")  GeometryKernel::snap_to "virtual void
GeometryKernel::snap_to(KernelPoint &point, GeometryHandle geom,
double *converged_tolerance=NULL, double *uvw_computed=NULL, double
*uvw_hint=NULL)=0 ";

%feature("docstring")  GeometryKernel::normal_at "virtual void
GeometryKernel::normal_at(KernelPoint &point, GeometryHandle geom,
std::vector< double > &normal)=0 ";

%feature("docstring")  GeometryKernel::is_curve "virtual bool
GeometryKernel::is_curve(GeometryHandle geom) const =0 ";

%feature("docstring")  GeometryKernel::is_surface "virtual bool
GeometryKernel::is_surface(GeometryHandle geom) const =0 ";


// File: classGeometryKernelOpenNURBS.xml
%feature("docstring") GeometryKernelOpenNURBS "";

%feature("docstring")
GeometryKernelOpenNURBS::GeometryKernelOpenNURBS "GeometryKernelOpenNURBS::GeometryKernelOpenNURBS() ";

%feature("docstring")
GeometryKernelOpenNURBS::~GeometryKernelOpenNURBS "GeometryKernelOpenNURBS::~GeometryKernelOpenNURBS() ";

%feature("docstring")  GeometryKernelOpenNURBS::read_file "bool
GeometryKernelOpenNURBS::read_file(const std::string &file_name,
std::vector< GeometryHandle > &geometry_entities) ";

%feature("docstring")  GeometryKernelOpenNURBS::get_attribute "std::string GeometryKernelOpenNURBS::get_attribute(GeometryHandle
geom) ";

%feature("docstring")  GeometryKernelOpenNURBS::snap_to "void
GeometryKernelOpenNURBS::snap_to(KernelPoint &point, GeometryHandle
geom, double *converged_tolerance=NULL, double *uvw_computed=NULL,
double *uvw_hint=NULL) ";

%feature("docstring")  GeometryKernelOpenNURBS::normal_at "void
GeometryKernelOpenNURBS::normal_at(KernelPoint &point, GeometryHandle
geom, std::vector< double > &normal) ";

%feature("docstring")  GeometryKernelOpenNURBS::is_curve "bool
GeometryKernelOpenNURBS::is_curve(GeometryHandle geom) const ";

%feature("docstring")  GeometryKernelOpenNURBS::is_surface "bool
GeometryKernelOpenNURBS::is_surface(GeometryHandle geom) const ";


// File: classGeometryKernelStupid.xml
%feature("docstring") GeometryKernelStupid "";

%feature("docstring")  GeometryKernelStupid::GeometryKernelStupid "GeometryKernelStupid::GeometryKernelStupid() ";

%feature("docstring")  GeometryKernelStupid::~GeometryKernelStupid "virtual GeometryKernelStupid::~GeometryKernelStupid() ";

%feature("docstring")  GeometryKernelStupid::read_file "virtual bool
GeometryKernelStupid::read_file(const std::string &file_name,
std::vector< GeometryHandle > &geometry_entities) ";

%feature("docstring")  GeometryKernelStupid::get_attribute "virtual
std::string GeometryKernelStupid::get_attribute(GeometryHandle geom)
";

%feature("docstring")  GeometryKernelStupid::snap_to "virtual void
GeometryKernelStupid::snap_to(KernelPoint &point, GeometryHandle geom)
";

%feature("docstring")  GeometryKernelStupid::normal_at "virtual void
GeometryKernelStupid::normal_at(KernelPoint &point, GeometryHandle
geom, std::vector< double > &normal) ";

%feature("docstring")  GeometryKernelStupid::is_curve "virtual bool
GeometryKernelStupid::is_curve(GeometryHandle geom) const ";

%feature("docstring")  GeometryKernelStupid::is_surface "virtual bool
GeometryKernelStupid::is_surface(GeometryHandle geom) const ";


// File: classstk_1_1percept_1_1GeometryVerifier.xml
%feature("docstring") stk::percept::GeometryVerifier "";

%feature("docstring")
stk::percept::GeometryVerifier::GeometryVerifier "stk::percept::GeometryVerifier::GeometryVerifier(bool dump=false,
double badJac=1.e-10) ";

%feature("docstring")  stk::percept::GeometryVerifier::isGeometryBad "bool stk::percept::GeometryVerifier::isGeometryBad(stk::mesh::BulkData
&bulk, bool printTable=false)

Check for nonpositive Jacobian

note: we're using cubature here instead of explicitly specifying some
reference points the idea is that we'll get a good estimate of the
Jacobian's sign by testing it at all the cubature points ";


// File: classstk_1_1percept_1_1GMeshSpec.xml
%feature("docstring") stk::percept::GMeshSpec "";

%feature("docstring")  stk::percept::GMeshSpec::GMeshSpec "stk::percept::GMeshSpec::GMeshSpec(const std::string &name) ";


// File: classGRAD.xml
%feature("docstring") GRAD "";


// File: classstk_1_1percept_1_1H1__NormOp.xml
%feature("docstring") stk::percept::H1_NormOp "";

%feature("docstring")  stk::percept::H1_NormOp::H1_NormOp "stk::percept::H1_NormOp::H1_NormOp(Function &integrand, int
spatialDim=3)

integrand tells what fields Intrepid should compute, etc. ";

%feature("docstring")  stk::percept::H1_NormOp::finalOp "void
stk::percept::H1_NormOp::finalOp(const std::vector< double > &vin,
std::vector< double > &vout) ";


// File: classstk_1_1percept_1_1H1Norm.xml
%feature("docstring") stk::percept::H1Norm "

compute the H1 norm or semi-norm

C++ includes: H1Norm.hpp ";

%feature("docstring")  stk::percept::H1Norm::H1Norm "stk::percept::H1Norm::H1Norm(mesh::BulkData &bulkData, std::string
partName, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false)
";

%feature("docstring")  stk::percept::H1Norm::H1Norm "stk::percept::H1Norm::H1Norm(mesh::BulkData &bulkData, MDArrayString
&partNames, TurboOption turboOpt=TURBO_NONE, bool
is_surface_norm=false) ";

%feature("docstring")  stk::percept::H1Norm::H1Norm "stk::percept::H1Norm::H1Norm(mesh::BulkData &bulkData, mesh::Part
*part=0, TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false)
";

%feature("docstring")  stk::percept::H1Norm::H1Norm "stk::percept::H1Norm::H1Norm(mesh::BulkData &bulkData, mesh::Selector
*selector, TurboOption turboOpt=TURBO_NONE, bool
is_surface_norm=false) ";


// File: classstk_1_1percept_1_1HasConstValue.xml
%feature("docstring") stk::percept::HasConstValue "";

%feature("docstring")  stk::percept::HasConstValue::getValue "virtual
ValueType& stk::percept::HasConstValue< ValueType >::getValue()=0 ";

%feature("docstring")  stk::percept::HasConstValue::~HasConstValue "virtual stk::percept::HasConstValue< ValueType >::~HasConstValue() ";


// File: classstk_1_1percept_1_1HasFinalOp.xml
%feature("docstring") stk::percept::HasFinalOp "";

%feature("docstring")  stk::percept::HasFinalOp::finalOp "virtual
void stk::percept::HasFinalOp< ValueType >::finalOp(const ValueType
&vin, ValueType &vout)=0 ";


// File: classstk_1_1percept_1_1HasValue.xml
%feature("docstring") stk::percept::HasValue "";

%feature("docstring")  stk::percept::HasValue::getValue "virtual
ValueType& stk::percept::HasValue< ValueType >::getValue()=0 ";

%feature("docstring")  stk::percept::HasValue::setValue "virtual void
stk::percept::HasValue< ValueType >::setValue(ValueType &)=0 ";

%feature("docstring")  stk::percept::HasValue::~HasValue "virtual
stk::percept::HasValue< ValueType >::~HasValue() ";


// File: classstk_1_1percept_1_1HeterogeneousFixture.xml
%feature("docstring") stk::percept::HeterogeneousFixture "

Use case with mixed element topologies and field relations to provide
fast access to node field data from an element.

copied from stk_mesh and modified

C++ includes: HeterogeneousFixture.hpp ";

%feature("docstring")
stk::percept::HeterogeneousFixture::~HeterogeneousFixture "stk::percept::HeterogeneousFixture::~HeterogeneousFixture() ";

%feature("docstring")
stk::percept::HeterogeneousFixture::HeterogeneousFixture "stk::percept::HeterogeneousFixture::HeterogeneousFixture(stk::ParallelMachine
comm, bool doCommit=true, bool do_sidesets=false) ";

%feature("docstring")  stk::percept::HeterogeneousFixture::populate "void stk::percept::HeterogeneousFixture::populate() ";


// File: classHex__.xml
%feature("docstring") Hex_ "";


// File: classstk_1_1adapt_1_1IAdapter.xml
%feature("docstring") stk::adapt::IAdapter "";

%feature("docstring")  stk::adapt::IAdapter::buildUnrefineList "virtual ElementUnrefineCollection
stk::adapt::IAdapter::buildUnrefineList()=0 ";


// File: classstk_1_1percept_1_1IdentityFunction.xml
%feature("docstring") stk::percept::IdentityFunction "";

%feature("docstring")
stk::percept::IdentityFunction::IdentityFunction "stk::percept::IdentityFunction::IdentityFunction() ";


// File: classstk_1_1adapt_1_1IEdgeAdapter.xml
%feature("docstring") stk::adapt::IEdgeAdapter "

An IEdgeAdapter is an abstract base class for derived classes that are
required to overload the mark method, which provides info such as the
element the edge belongs to, which edge ordinal it is, the nodes of
the edge and the edge coordinates.

C++ includes: IEdgeAdapter.hpp ";

%feature("docstring")  stk::adapt::IEdgeAdapter::IEdgeAdapter "stk::adapt::IEdgeAdapter::IEdgeAdapter(percept::PerceptMesh &eMesh,
UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")  stk::adapt::IEdgeAdapter::buildUnrefineList "ElementUnrefineCollection
stk::adapt::IEdgeAdapter::buildUnrefineList()

can be overriden ";


// File: structstk_1_1adapt_1_1IEdgeBasedAdapterPredicate.xml
%feature("docstring") stk::adapt::IEdgeBasedAdapterPredicate "

Signatures for predicate objects that can be used to select entities
(elements, edges, faces,...) for refinement or unrefinement. The class
must supply an operator() that takes an entity and decides on whether
it should be refined, or unrefined, or ignored.

The Selector pattern as shown below is useful for selecting only
entities that belong to particular mesh Parts, for example, or any
other definition of Selector.

We follow the unary_function pattern to enable the structs to be used
in STL algorithms that know about unary_functions.

The following are a couple of examples of what refine and unrefine
predicates might look like.

Following these examples we show the prototype for the operations that
are performed on these predicates.

C++ includes: IEdgeBasedAdapterPredicate.hpp ";


// File: classstk_1_1adapt_1_1IElementAdapter.xml
%feature("docstring") stk::adapt::IElementAdapter "

An IElementAdapter is an abstract base class for derived classes that
are required to overload the mark method, which supplies the derived
class with the element to be marked for refine, unrefine, or both (
See:  IAdapter::AdaptInstruction)

C++ includes: IElementAdapter.hpp ";

%feature("docstring")  stk::adapt::IElementAdapter::IElementAdapter "stk::adapt::IElementAdapter::IElementAdapter(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")  stk::adapt::IElementAdapter::buildUnrefineList
"ElementUnrefineCollection
stk::adapt::IElementAdapter::buildUnrefineList() ";


// File: structstk_1_1adapt_1_1IElementBasedAdapterPredicate.xml
%feature("docstring") stk::adapt::IElementBasedAdapterPredicate "

Signatures for predicate objects that can be used to select entities
(elements, edges, faces,...) for refinement or unrefinement. The class
must supply an operator() that takes an entity and decides on whether
it should be refined, or unrefined, or ignored.

The Selector pattern as shown below is useful for selecting only
entities that belong to particular mesh Parts, for example, or any
other definition of Selector.

We follow the unary_function pattern to enable the structs to be used
in STL algorithms that know about unary_functions.

The following are a couple of examples of what refine and unrefine
predicates might look like.

Following these examples we show the prototype for the operations that
are performed on these predicates.

C++ includes: IElementBasedAdapterPredicate.hpp ";


// File: structInit.xml
%feature("docstring") Init "";


// File: classstk_1_1percept_1_1IntrepidManager_1_1Integral.xml
%feature("docstring") stk::percept::IntrepidManager::Integral "

([C])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::Integral::Integral "stk::percept::IntrepidManager::Integral::Integral(IM &im)

([C]) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1IntegralDOF.xml
%feature("docstring") stk::percept::IntrepidManager::IntegralDOF "

([C], [DOF])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::IntegralDOF::IntegralDOF "stk::percept::IntrepidManager::IntegralDOF::IntegralDOF(IM &im)

([C], [DOF]) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1IntegrandValues.xml
%feature("docstring") stk::percept::IntrepidManager::IntegrandValues "

([C], [P])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::IntegrandValues::IntegrandValues "stk::percept::IntrepidManager::IntegrandValues::IntegrandValues(IM
&im)

([C], [P]) ";

%feature("docstring")
stk::percept::IntrepidManager::IntegrandValues::copyFrom "void
stk::percept::IntrepidManager::IntegrandValues::copyFrom(MDArray &mda)
";

%feature("docstring")
stk::percept::IntrepidManager::IntegrandValues::copyFrom "void
stk::percept::IntrepidManager::IntegrandValues::copyFrom(IntrepidManager
&im, MDArray &mda, int iDof) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1IntegrandValuesDOF.xml
%feature("docstring")
stk::percept::IntrepidManager::IntegrandValuesDOF "

([C], [P], [DOF])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::IntegrandValuesDOF::IntegrandValuesDOF
"stk::percept::IntrepidManager::IntegrandValuesDOF::IntegrandValuesDOF(IM
&im)

([C], [P], [DOF]) ";

%feature("docstring")
stk::percept::IntrepidManager::IntegrandValuesDOF::copyFrom "void
stk::percept::IntrepidManager::IntegrandValuesDOF::copyFrom(MDArray
&mda) ";


// File: classstk_1_1percept_1_1IntegratedOp.xml
%feature("docstring") stk::percept::IntegratedOp "";

%feature("docstring")  stk::percept::IntegratedOp::IntegratedOp "stk::percept::IntegratedOp::IntegratedOp(Function &integrand,
TurboOption turboOpt=TURBO_NONE, mesh::FieldBase *field=0) ";

%feature("docstring")  stk::percept::IntegratedOp::setAccumulationType
"void
stk::percept::IntegratedOp::setAccumulationType(AccumulationType type)
";

%feature("docstring")  stk::percept::IntegratedOp::getAccumulationType
"AccumulationType stk::percept::IntegratedOp::getAccumulationType()
";

%feature("docstring")  stk::percept::IntegratedOp::setCubDegree "void
stk::percept::IntegratedOp::setCubDegree(unsigned cubDegree) ";

%feature("docstring")  stk::percept::IntegratedOp::getCubDegree "unsigned stk::percept::IntegratedOp::getCubDegree() ";

%feature("docstring")  stk::percept::IntegratedOp::init "void
stk::percept::IntegratedOp::init() ";

%feature("docstring")  stk::percept::IntegratedOp::getValue "std::vector<double>& stk::percept::IntegratedOp::getValue(void) ";

%feature("docstring")  stk::percept::IntegratedOp::getElementCount "unsigned stk::percept::IntegratedOp::getElementCount() ";

%feature("docstring")  stk::percept::IntegratedOp::init_elementOp "void stk::percept::IntegratedOp::init_elementOp() ";

%feature("docstring")  stk::percept::IntegratedOp::fini_elementOp "void stk::percept::IntegratedOp::fini_elementOp() ";


// File: classstk_1_1percept_1_1IntrepidManager.xml
%feature("docstring") stk::percept::IntrepidManager "

|-------------------------------------------------------------------------------------------------|
| Index type | Dimension | Description |
|---------------------------|-----------|---------------------------------------------------------|
| point | [P] | number of points stored in an MD array | | vertex |
[V] | number of nodes stored in an MD aray | | field | [F] | number of
fields stored in an MD array | | basis field | [ B] | number of basis
fields stored in an MD array | | cell | [C] | number of cells stored
in an MD array | | field coordinate | [D] | space dimension | |
derivative ordinal | [K] | cardinality of the set of kth derivatives |
| | | | | dof | [DOF] | number of DOFs stored in an MD array |
|-------------------------------------------------------------------------------------------------|

Note: Intrepid really doesn't have a concept of \"DOF\" at a node.
It's either a single variable, or a vector- or tensor-valued variable.
So, no DOF-related arrays as used herein can be used with Intrepid -
you must call Intrepd one DOF at a time.

FieldContainer<double> cub_points(numCubPoints, spaceDim);
FieldContainer<double> cub_weights(numCubPoints);

FieldContainer<double> cell_nodes(numCells, numNodes, spaceDim);

FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim,
spaceDim); FieldContainer<double> jacobian_inv(numCells, numCubPoints,
spaceDim, spaceDim); FieldContainer<double> jacobian_det(numCells,
numCubPoints); FieldContainer<double> weighted_measure(numCells,
numCubPoints);

FieldContainer<double> grad_at_cub_points(numFields, numCubPoints,
spaceDim); FieldContainer<double>
transformed_grad_at_cub_points(numCells, numFields, numCubPoints,
spaceDim); FieldContainer<double>
weighted_transformed_grad_at_cub_points(numCells, numFields,
numCubPoints, spaceDim); FieldContainer<double>
stiffness_matrices(numCells, numFields, numFields);

C++ includes: IntrepidManager.hpp ";

%feature("docstring")  stk::percept::IntrepidManager::IntrepidManager
"stk::percept::IntrepidManager::IntrepidManager(Elements_Tag el,
Cub_Points_Tag ct, NodesPerElem_Tag nc, Spatial_Dim_Tag st, DOFs_Tag
dt) ";

%feature("docstring")  stk::percept::IntrepidManager::IntrepidManager
"stk::percept::IntrepidManager::IntrepidManager(Elements_Tag el,
CellTopology &cellTopo, unsigned cubDegree=2) ";

%feature("docstring")  stk::percept::IntrepidManager::setupCubature "void stk::percept::IntrepidManager::setupCubature(CellTopology
&cellTopo, unsigned cubDegree=2) ";


// File: classstk_1_1percept_1_1IsInElement.xml
%feature("docstring") stk::percept::IsInElement "";

%feature("docstring")  stk::percept::IsInElement::IsInElement "stk::percept::IsInElement::IsInElement(MDArray &input_phy_points,
MDArray &found_parametric_coordinates) ";

%feature("docstring")  stk::percept::IsInElement::init_elementOp "void stk::percept::IsInElement::init_elementOp() ";

%feature("docstring")  stk::percept::IsInElement::fini_elementOp "void stk::percept::IsInElement::fini_elementOp() ";


// File: structstk_1_1percept_1_1jacData.xml
%feature("docstring") stk::percept::jacData "";

%feature("docstring")  stk::percept::jacData::jacData "stk::percept::jacData::jacData() ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1Jacobian.xml
%feature("docstring") stk::percept::IntrepidManager::Jacobian "

([C], [P], [D], [D])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::Jacobian::Jacobian "stk::percept::IntrepidManager::Jacobian::Jacobian(IM &im)

([C], [P], [D], [D]) ";

%feature("docstring")  stk::percept::IntrepidManager::Jacobian::copyTo
"void stk::percept::IntrepidManager::Jacobian::copyTo(MDArray &mda)
";


// File: classstk_1_1percept_1_1IntrepidManager_1_1JacobianDet.xml
%feature("docstring") stk::percept::IntrepidManager::JacobianDet "

([C], [P])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::JacobianDet::JacobianDet "stk::percept::IntrepidManager::JacobianDet::JacobianDet(IM &im)

([C], [P]) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1JacobianInverse.xml
%feature("docstring") stk::percept::IntrepidManager::JacobianInverse "

([C], [P], [D], [D])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::JacobianInverse::JacobianInverse "stk::percept::IntrepidManager::JacobianInverse::JacobianInverse(IM
&im)

([C], [P], [D], [D]) ";


// File: classstk_1_1percept_1_1l2NormOpScalar.xml
%feature("docstring") stk::percept::l2NormOpScalar "";

%feature("docstring")  stk::percept::l2NormOpScalar::l2NormOpScalar "stk::percept::l2NormOpScalar::l2NormOpScalar(FieldBase *field) ";


// File: classstk_1_1percept_1_1l2NormOpScalarFunction.xml
%feature("docstring") stk::percept::l2NormOpScalarFunction "";

%feature("docstring")
stk::percept::l2NormOpScalarFunction::l2NormOpScalarFunction "stk::percept::l2NormOpScalarFunction::l2NormOpScalarFunction(Function
&f) ";


// File: classLine__.xml
%feature("docstring") Line_ "";


// File: classstk_1_1percept_1_1LN__NormOp.xml
%feature("docstring") stk::percept::LN_NormOp "

for Power = -1, compute the inf-norm

C++ includes: Norm.hpp ";

%feature("docstring")  stk::percept::LN_NormOp::LN_NormOp "stk::percept::LN_NormOp< Power >::LN_NormOp(Function &integrand)

integrand tells what fields Intrepid should compute, etc. Note: this
function is intended to be used to wrap a Function using
CompositeFunction and thus its domain and codomain are the same as the
wrapped function's codomain ";

%feature("docstring")  stk::percept::LN_NormOp::finalOp "void
stk::percept::LN_NormOp< Power >::finalOp(const std::vector< double >
&vin, std::vector< double > &vout) ";


// File: structstk_1_1percept_1_1unit__tests_1_1LocalFixture.xml
%feature("docstring") stk::percept::unit_tests::LocalFixture "";

%feature("docstring")
stk::percept::unit_tests::LocalFixture::LocalFixture "stk::percept::unit_tests::LocalFixture::LocalFixture(size_t num_xyz=4,
size_t num_y=0, size_t num_z=0, bool sidesets=false, bool commit=true)
";

%feature("docstring")  stk::percept::unit_tests::LocalFixture::init "int stk::percept::unit_tests::LocalFixture::init(size_t num_xyz,
size_t num_y_arg, size_t num_z_arg, bool sidesets=false, bool
commit=true) ";


// File: structstk_1_1percept_1_1lstr.xml
%feature("docstring") stk::percept::lstr "";


// File: structstk_1_1percept_1_1ltstr.xml
%feature("docstring") stk::percept::ltstr "";


// File: classstk_1_1percept_1_1Math.xml
%feature("docstring") stk::percept::Math "";


// File: classstk_1_1percept_1_1MaxOfNodeValues.xml
%feature("docstring") stk::percept::MaxOfNodeValues "";

%feature("docstring")  stk::percept::MaxOfNodeValues::MaxOfNodeValues
"stk::percept::MaxOfNodeValues::MaxOfNodeValues(int spatialDim,
Function &integrand) ";


// File: classstk_1_1percept_1_1MDArrayString.xml
%feature("docstring") stk::percept::MDArrayString "";

%feature("docstring")  stk::percept::MDArrayString::MDArrayString "stk::percept::MDArrayString::MDArrayString() ";

%feature("docstring")  stk::percept::MDArrayString::MDArrayString "stk::percept::MDArrayString::MDArrayString(int dim) ";

%feature("docstring")  stk::percept::MDArrayString::MDArrayString "stk::percept::MDArrayString::MDArrayString(int dim0, int dim1) ";

%feature("docstring")  stk::percept::MDArrayString::MDArrayString "stk::percept::MDArrayString::MDArrayString(const MDArrayString &mda)
";

%feature("docstring")  stk::percept::MDArrayString::resize "void
stk::percept::MDArrayString::resize(int dim) ";

%feature("docstring")  stk::percept::MDArrayString::resize "void
stk::percept::MDArrayString::resize(int dim0, int dim1) ";

%feature("docstring")  stk::percept::MDArrayString::rank "const int
stk::percept::MDArrayString::rank() const ";

%feature("docstring")  stk::percept::MDArrayString::setValues "void
stk::percept::MDArrayString::setValues(std::string *data) ";

%feature("docstring")  stk::percept::MDArrayString::dimension "int
stk::percept::MDArrayString::dimension(int i1) const ";


// File: structstk_1_1adapt_1_1regression__tests_1_1MemoryInfo.xml
%feature("docstring") stk::adapt::regression_tests::MemoryInfo "";

%feature("docstring")
stk::adapt::regression_tests::MemoryInfo::MemoryInfo "stk::adapt::regression_tests::MemoryInfo::MemoryInfo() ";

%feature("docstring")
stk::adapt::regression_tests::MemoryInfo::get_memory_usage "void
stk::adapt::regression_tests::MemoryInfo::get_memory_usage() ";

%feature("docstring")
stk::adapt::regression_tests::MemoryInfo::set_state "void
stk::adapt::regression_tests::MemoryInfo::set_state() ";

%feature("docstring")
stk::adapt::regression_tests::MemoryInfo::get_increment "void
stk::adapt::regression_tests::MemoryInfo::get_increment() ";


// File: structstk_1_1adapt_1_1MemoryMultipliers.xml
%feature("docstring") stk::adapt::MemoryMultipliers "";

%feature("docstring")
stk::adapt::MemoryMultipliers::MemoryMultipliers "stk::adapt::MemoryMultipliers::MemoryMultipliers(MemMultType
mult_hex8=1490, MemMultType mult_tet4=702, MemMultType mult_nodes=0)
";

%feature("docstring")  stk::adapt::MemoryMultipliers::read_simple "void stk::adapt::MemoryMultipliers::read_simple(std::string file_name)
";

%feature("docstring")  stk::adapt::MemoryMultipliers::estimate_memory
"MemorySizeType stk::adapt::MemoryMultipliers::estimate_memory() ";

%feature("docstring")  stk::adapt::MemoryMultipliers::estimate_memory
"MemorySizeType
stk::adapt::MemoryMultipliers::estimate_memory(std::vector<
RefinementInfoByType > &refInfo) ";


// File: classstk_1_1percept_1_1MeshDifference.xml
%feature("docstring") stk::percept::MeshDifference "";

%feature("docstring")  stk::percept::MeshDifference::MeshDifference "stk::percept::MeshDifference::MeshDifference() ";

%feature("docstring")  stk::percept::MeshDifference::run "void
stk::percept::MeshDifference::run(int argc, char **argv) ";

%feature("docstring")  stk::percept::MeshDifference::process_options "void stk::percept::MeshDifference::process_options(RunEnvironment &re)
";


// File: classMeshGeometry.xml
%feature("docstring") MeshGeometry "";

%feature("docstring")  MeshGeometry::MeshGeometry "MeshGeometry::MeshGeometry(GeometryKernel *geom, double
doCheckMovement=0.0, double doCheckCpuTime=0.0, bool
cache_bucket_selectors_is_active=false, bool doPrint=false) ";

%feature("docstring")  MeshGeometry::~MeshGeometry "MeshGeometry::~MeshGeometry() ";

%feature("docstring")  MeshGeometry::print_node_movement_summary "void MeshGeometry::print_node_movement_summary() ";

%feature("docstring")  MeshGeometry::add_evaluator "void
MeshGeometry::add_evaluator(GeometryEvaluator *evaluator) ";

%feature("docstring")  MeshGeometry::add_evaluators "void
MeshGeometry::add_evaluators(std::vector< GeometryEvaluator * >
evaluators) ";

%feature("docstring")  MeshGeometry::snap_points_to_geometry "void
MeshGeometry::snap_points_to_geometry(PerceptMesh *mesh_data) ";

%feature("docstring")  MeshGeometry::snap_points_to_geometry "void
MeshGeometry::snap_points_to_geometry(PerceptMesh *mesh_data,
std::vector< stk::mesh::Entity * > &nodes) ";

%feature("docstring")  MeshGeometry::normal_at "void
MeshGeometry::normal_at(PerceptMesh *eMesh, stk::mesh::Entity *node,
std::vector< double > &normal) ";

%feature("docstring")  MeshGeometry::classify_node "int
MeshGeometry::classify_node(const stk::mesh::Entity &node, size_t
&curveOrSurfaceEvaluator)

Return 0,1,2,3 if the node or bucket is on a geometry vertex, curve,
surface or domain. Return the found evaluators in the curveEvaluators
and surfEvaluators.

Return 0,1,2,3 if the node is on a geometry vertex, curve, surface or
domain. ";

%feature("docstring")  MeshGeometry::classify_bucket "int
MeshGeometry::classify_bucket(const stk::mesh::Bucket &bucket, size_t
&curveOrSurfaceEvaluator)

Return 0,1,2,3 if the bucket is on a geometry vertex, curve, surface
or domain. ";

%feature("docstring")  MeshGeometry::getGeomEvaluators "const
std::vector< GeometryEvaluator * > & MeshGeometry::getGeomEvaluators()
";


// File: classstk_1_1percept_1_1MeshGeometryVerifier.xml
%feature("docstring") stk::percept::MeshGeometryVerifier "";


// File: classstk_1_1adapt_1_1Elem_1_1MeshObjRefinementTopology.xml
%feature("docstring") stk::adapt::Elem::MeshObjRefinementTopology "

Class to hold refinement information for partial and heterogeneous

C++ includes: RefinementTopology.hpp ";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::MeshObjRefinementTopology
"stk::adapt::Elem::MeshObjRefinementTopology::MeshObjRefinementTopology()
";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::~MeshObjRefinementTopology
"stk::adapt::Elem::MeshObjRefinementTopology::~MeshObjRefinementTopology()
";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::num_child "UInt
stk::adapt::Elem::MeshObjRefinementTopology::num_child() const ";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::num_child_nodes "UInt
stk::adapt::Elem::MeshObjRefinementTopology::num_child_nodes() const
";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::child_cell_topology "CellTopology
stk::adapt::Elem::MeshObjRefinementTopology::child_cell_topology(UInt
child) const ";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::child_node "const UInt *
stk::adapt::Elem::MeshObjRefinementTopology::child_node(UInt child)
const ";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::homogeneous_refinement "bool
stk::adapt::Elem::MeshObjRefinementTopology::homogeneous_refinement()
const ";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::full_refinement "bool
stk::adapt::Elem::MeshObjRefinementTopology::full_refinement() const
";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::child_face "std::pair<
UInt, UInt >
stk::adapt::Elem::MeshObjRefinementTopology::child_face(const UInt
face_ordinal, const UInt face_child_ordinal, const Elem::CellTopology
&objTop, const RefinementKey &objDesiredKey) const

Mapping of parent->face->child to parent->child->face face_ordinal <
num_faces()

face_child_ordinal < face_topology(face_ordinal)->num_child() objTop
need to have objects regular topology as well

(child_ordinal, child_face_ordinal) ";

%feature("docstring")
stk::adapt::Elem::MeshObjRefinementTopology::child_edge "std::pair<
UInt, UInt >
stk::adapt::Elem::MeshObjRefinementTopology::child_edge(const UInt
edge_ordinal, const UInt edge_child_ordinal, const Elem::CellTopology
&objTop) const

Mapping of parent->edge->child to parent->child->edge edge_ordinal <
getEdgeCount()

edge_child_ordinal < edge_topology(edge_ordinal)->num_child() objTop
need to have objects regular topology as well

(child_ordinal, child_edge_ordinal) ";


// File: classstk_1_1adapt_1_1Elem_1_1MeshObjTopology.xml
%feature("docstring") stk::adapt::Elem::MeshObjTopology "";

%feature("docstring")
stk::adapt::Elem::MeshObjTopology::MeshObjTopology "stk::adapt::Elem::MeshObjTopology::MeshObjTopology(const
CellTopologyData *cell_topology_data) ";

%feature("docstring")
stk::adapt::Elem::MeshObjTopology::~MeshObjTopology "stk::adapt::Elem::MeshObjTopology::~MeshObjTopology() ";

%feature("docstring")
stk::adapt::Elem::MeshObjTopology::getCellTopology "Elem::CellTopology
stk::adapt::Elem::MeshObjTopology::getCellTopology() const ";


// File: classstk_1_1percept_1_1MeshTopologyVerifier.xml
%feature("docstring") stk::percept::MeshTopologyVerifier "";


// File: classstk_1_1percept_1_1MeshTransformer.xml
%feature("docstring") stk::percept::MeshTransformer "";

%feature("docstring")  stk::percept::MeshTransformer::MeshTransformer
"stk::percept::MeshTransformer::MeshTransformer() ";

%feature("docstring")  stk::percept::MeshTransformer::MeshTransformer
"stk::percept::MeshTransformer::MeshTransformer(Math::Matrix &m) ";


// File: classstk_1_1percept_1_1MeshUtil.xml
%feature("docstring") stk::percept::MeshUtil "";


// File: classstk_1_1percept_1_1MeshVerifier.xml
%feature("docstring") stk::percept::MeshVerifier "";


// File: structstk_1_1percept_1_1minMaxAve.xml
%feature("docstring") stk::percept::minMaxAve "";

%feature("docstring")  stk::percept::minMaxAve::minMaxAve "stk::percept::minMaxAve::minMaxAve() ";

%feature("docstring")  stk::percept::minMaxAve::registerValue "void
stk::percept::minMaxAve::registerValue(unsigned id, double val) ";

%feature("docstring")  stk::percept::minMaxAve::finish "void
stk::percept::minMaxAve::finish(stk::mesh::BulkData &mesh) ";

%feature("docstring")  stk::percept::minMaxAve::setStandardRanges "void stk::percept::minMaxAve::setStandardRanges() ";


// File: classMu2__.xml
%feature("docstring") Mu2_ "";

%feature("docstring")  Mu2_::Mu2_ "Mu2_< Temperature
>::Mu2_(Temperature temperature) ";

%feature("docstring")  Mu2_::Mu2_ "Mu2_< Temperature
>::Mu2_(Temperature temperature) ";


// File: classMu__.xml
%feature("docstring") Mu_ "";

%feature("docstring")  Mu_::Mu_ "Mu_< Temperature >::Mu_(Temperature
temperature) ";

%feature("docstring")  Mu_::Mu_ "Mu_< Temperature >::Mu_(Temperature
temperature) ";


// File: classstk_1_1percept_1_1MultipleFieldFunction.xml
%feature("docstring") stk::percept::MultipleFieldFunction "

This class compbines several Fields into a single function, e.g.,
density_field*temperature_field/pressure_field It uses Intrepid to
compute basis functions then evaluates all the required fields and
stores them in MDArray's ready for the operator() to compute the
actual function.

C++ includes: MultipleFieldFunction.hpp ";


// File: structstk_1_1adapt_1_1my__equal__to.xml
%feature("docstring") stk::adapt::my_equal_to "";


// File: structstk_1_1adapt_1_1my__equal__to__old.xml
%feature("docstring") stk::adapt::my_equal_to_old "";


// File: structstk_1_1adapt_1_1my__fast__equal__to.xml
%feature("docstring") stk::adapt::my_fast_equal_to "";


// File: structstk_1_1adapt_1_1my__fast__hash.xml
%feature("docstring") stk::adapt::my_fast_hash "";


// File: structstk_1_1adapt_1_1my__hash.xml
%feature("docstring") stk::adapt::my_hash "";


// File: structstk_1_1adapt_1_1my__hash__old.xml
%feature("docstring") stk::adapt::my_hash_old "";


// File: structstk_1_1adapt_1_1my__tuple__hash.xml
%feature("docstring") stk::adapt::my_tuple_hash "";


// File: classstk_1_1percept_1_1MyEdge.xml
%feature("docstring") stk::percept::MyEdge "";

%feature("docstring")  stk::percept::MyEdge::MyEdge "stk::percept::MyEdge< IdType >::MyEdge(IdType i0, IdType i1) ";

%feature("docstring")  stk::percept::MyEdge::getId0 "IdType
stk::percept::MyEdge< IdType >::getId0() const ";

%feature("docstring")  stk::percept::MyEdge::getId1 "IdType
stk::percept::MyEdge< IdType >::getId1() const ";


// File: structstk_1_1adapt_1_1regression__tests_1_1MyEdgeBasedRefinePredicate.xml
%feature("docstring")
stk::adapt::regression_tests::MyEdgeBasedRefinePredicate "";

%feature("docstring")
stk::adapt::regression_tests::MyEdgeBasedRefinePredicate::MyEdgeBasedRefinePredicate
"stk::adapt::regression_tests::MyEdgeBasedRefinePredicate::MyEdgeBasedRefinePredicate(stk::mesh::Selector
*selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) ";


// File: structstk_1_1adapt_1_1myVec.xml
%feature("docstring") stk::adapt::myVec "";


// File: classstk_1_1percept_1_1Math_1_1MyVector.xml
%feature("docstring") stk::percept::Math::MyVector "";

%feature("docstring")  stk::percept::Math::MyVector::MyVector "stk::percept::Math::MyVector::MyVector(double x=0.0) ";

%feature("docstring")  stk::percept::Math::MyVector::MyVector "stk::percept::Math::MyVector::MyVector(double *x) ";


// File: classstk_1_1percept_1_1Name.xml
%feature("docstring") stk::percept::Name "

Useful in other places where two strings are passed into a function or
constructor.

this is to avoid a common bug where the name of the String Function is
given instead of function_string, ie., first two args are accidentally
reversed - this is essentially a model of a \"named argument\", as
opposed to a positional one; of course, it is still only a hint
(though a strong one) to the user

C++ includes: Name.hpp ";

%feature("docstring")  stk::percept::Name::Name "stk::percept::Name::Name(const std::string name) ";

%feature("docstring")  stk::percept::Name::getName "const
std::string& stk::percept::Name::getName() const ";


// File: structstk_1_1adapt_1_1NodeIdsOnSubDimEntityType.xml
%feature("docstring") stk::adapt::NodeIdsOnSubDimEntityType "

data on a sub-dim entity (global node ids on the entity, the owning
element's id)

C++ includes: NodeRegistry.hpp ";

%feature("docstring")
stk::adapt::NodeIdsOnSubDimEntityType::NodeIdsOnSubDimEntityType "stk::adapt::NodeIdsOnSubDimEntityType::NodeIdsOnSubDimEntityType(unsigned
sz=1, NodeIdsOnSubDimEntityTypeQuantum allValues=0) ";

%feature("docstring")  stk::adapt::NodeIdsOnSubDimEntityType::resize "void stk::adapt::NodeIdsOnSubDimEntityType::resize(size_t sz) ";

%feature("docstring")  stk::adapt::NodeIdsOnSubDimEntityType::pack "void stk::adapt::NodeIdsOnSubDimEntityType::pack(CommBuffer &buff) ";

%feature("docstring")  stk::adapt::NodeIdsOnSubDimEntityType::unpack "void stk::adapt::NodeIdsOnSubDimEntityType::unpack(PerceptMesh &eMesh,
CommBuffer &buff) ";


// File: classstk_1_1adapt_1_1NodeRegistry.xml
%feature("docstring") stk::adapt::NodeRegistry "";

%feature("docstring")  stk::adapt::NodeRegistry::NodeRegistry "stk::adapt::NodeRegistry::NodeRegistry(percept::PerceptMesh &eMesh,
bool useCustomGhosting=false) ";

%feature("docstring")  stk::adapt::NodeRegistry::~NodeRegistry "stk::adapt::NodeRegistry::~NodeRegistry() ";

%feature("docstring")  stk::adapt::NodeRegistry::init_comm_all "void
stk::adapt::NodeRegistry::init_comm_all() ";

%feature("docstring")  stk::adapt::NodeRegistry::init_entity_repo "void stk::adapt::NodeRegistry::init_entity_repo() ";

%feature("docstring")  stk::adapt::NodeRegistry::clear_dangling_nodes
"void stk::adapt::NodeRegistry::clear_dangling_nodes(SetOfEntities
*nodes_to_be_deleted) ";

%feature("docstring")  stk::adapt::NodeRegistry::initialize "void
stk::adapt::NodeRegistry::initialize() ";

%feature("docstring")  stk::adapt::NodeRegistry::beginRegistration "void stk::adapt::NodeRegistry::beginRegistration() ";

%feature("docstring")  stk::adapt::NodeRegistry::endRegistration "void stk::adapt::NodeRegistry::endRegistration() ";

%feature("docstring")  stk::adapt::NodeRegistry::beginLocalMeshMods "void stk::adapt::NodeRegistry::beginLocalMeshMods() ";

%feature("docstring")  stk::adapt::NodeRegistry::endLocalMeshMods "void stk::adapt::NodeRegistry::endLocalMeshMods() ";

%feature("docstring")  stk::adapt::NodeRegistry::beginCheckForRemote "void stk::adapt::NodeRegistry::beginCheckForRemote() ";

%feature("docstring")  stk::adapt::NodeRegistry::endCheckForRemote "void stk::adapt::NodeRegistry::endCheckForRemote() ";

%feature("docstring")  stk::adapt::NodeRegistry::beginGetFromRemote "void stk::adapt::NodeRegistry::beginGetFromRemote() ";

%feature("docstring")  stk::adapt::NodeRegistry::endGetFromRemote "void stk::adapt::NodeRegistry::endGetFromRemote() ";

%feature("docstring")  stk::adapt::NodeRegistry::removePseudoEntities
"void stk::adapt::NodeRegistry::removePseudoEntities() ";

%feature("docstring")
stk::adapt::NodeRegistry::setAllReceivedNodeData "void
stk::adapt::NodeRegistry::setAllReceivedNodeData() ";

%feature("docstring")
stk::adapt::NodeRegistry::removeUnmarkedSubDimEntities "void
stk::adapt::NodeRegistry::removeUnmarkedSubDimEntities()

when a sub-dim entity is visited during node registration but is
flagged as not being marked, and thus not requiring any new nodes, we
flag it with NR_MARK_NONE, then remove it here ";

%feature("docstring")  stk::adapt::NodeRegistry::is_empty "bool
stk::adapt::NodeRegistry::is_empty(const stk::mesh::Entity &element,
stk::mesh::EntityRank needed_entity_rank, unsigned iSubDimOrd) ";

%feature("docstring")  stk::adapt::NodeRegistry::registerNeedNewNode "bool stk::adapt::NodeRegistry::registerNeedNewNode(const
stk::mesh::Entity &element, NeededEntityType &needed_entity_rank,
unsigned iSubDimOrd, bool needNodes)

Register the need for a new node on the sub-dimensional entity

Parameters:
-----------

subDimEntity:  on element

element.:  If the element is a ghost element, the entity is still
registered: the locality/ownership of the new entity can be determined
by the locality of the element (ghost or not).

once it's in, the assertion should be: owning_elementId <
non_owning_elementId && owning_elementRank >= non_owning_elementRank
";

%feature("docstring")  stk::adapt::NodeRegistry::checkForRemote "bool
stk::adapt::NodeRegistry::checkForRemote(const stk::mesh::Entity
&element, NeededEntityType &needed_entity_rank, unsigned iSubDimOrd,
bool needNodes_notUsed)

check the newly registered node from the registry, which does one of
three things, depending on what mode we are in: 1. counts buffer in
prep for sending (just does a pack) 2. packs the buffer (after buffers
are alloc'd) 3. returns the new node after all communications are done
";

%feature("docstring")  stk::adapt::NodeRegistry::getFromRemote "bool
stk::adapt::NodeRegistry::getFromRemote(const stk::mesh::Entity
&element, NeededEntityType &needed_entity_rank, unsigned iSubDimOrd,
bool needNodes_notUsed) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_using_find
"LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_using_find(stk::mesh::EntityRank
&rank, const stk::mesh::EntityId &id) const ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity(stk::mesh::BulkData &bulk,
stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_I "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_I(stk::mesh::BulkData &bulk,
stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_Ia "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_Ia(stk::mesh::BulkData &bulk,
stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_Ib "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_Ib(stk::mesh::BulkData &bulk,
stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_element "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_element(stk::mesh::BulkData
&bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_node_I "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_node_I(stk::mesh::BulkData &bulk,
stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_node_Ia "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_node_Ia(stk::mesh::BulkData
&bulk, stk::mesh::EntityRank rank, NodeIdsOnSubDimEntityType
&nodeIds_onSE, unsigned index) ";

%feature("docstring")  stk::adapt::NodeRegistry::get_entity_node_Ib "LOCAL_INLINE stk::mesh::Entity*
stk::adapt::NodeRegistry::get_entity_node_Ib(stk::mesh::BulkData
&bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id) ";

%feature("docstring")
stk::adapt::NodeRegistry::getNewNodesOnSubDimEntity "NodeIdsOnSubDimEntityType *
stk::adapt::NodeRegistry::getNewNodesOnSubDimEntity(const
stk::mesh::Entity &element, stk::mesh::EntityRank &needed_entity_rank,
unsigned iSubDimOrd) ";

%feature("docstring")  stk::adapt::NodeRegistry::makeCentroidCoords "void stk::adapt::NodeRegistry::makeCentroidCoords(const
stk::mesh::Entity &element, stk::mesh::EntityRank needed_entity_rank,
unsigned iSubDimOrd)

makes coordinates of this new node be the centroid of its sub entity
";

%feature("docstring")  stk::adapt::NodeRegistry::makeCentroidField "void stk::adapt::NodeRegistry::makeCentroidField(const
stk::mesh::Entity &element, stk::mesh::EntityRank needed_entity_rank,
unsigned iSubDimOrd, stk::mesh::FieldBase *field)

!stkmesh::Entity * node =
get_entity_node_II(*m_eMesh.get_bulk_data(),Node, nodeId); ";

%feature("docstring")  stk::adapt::NodeRegistry::makeCentroid "void
stk::adapt::NodeRegistry::makeCentroid(stk::mesh::FieldBase *field)

makes coordinates of this new node be the centroid of its sub entity -
this version does it for all new nodes

!element_p = get_entity_element(*m_eMesh.get_bulk_data(),
m_eMesh.element_rank(), elementId);

!stkmesh::Entity * node = get_entity_node_II(*m_eMesh.get_bulk_data(),
mesh::Node, nodeId); ";

%feature("docstring")  stk::adapt::NodeRegistry::interpolateFields "void stk::adapt::NodeRegistry::interpolateFields(const
stk::mesh::Entity &element, stk::mesh::EntityRank needed_entity_rank,
unsigned iSubDimOrd)

do interpolation for all fields ";

%feature("docstring")  stk::adapt::NodeRegistry::interpolateFields "void stk::adapt::NodeRegistry::interpolateFields()

do interpolation for all fields ";

%feature("docstring")  stk::adapt::NodeRegistry::addToExistingParts "void stk::adapt::NodeRegistry::addToExistingParts(const
stk::mesh::Entity &element, stk::mesh::EntityRank needed_entity_rank,
unsigned iSubDimOrd)

check for adding new nodes to existing parts based on sub-entity part
ownership

!Entity * node = get_entity_node_II(*m_eMesh.get_bulk_data(),Node,
nodeId); ";

%feature("docstring")  stk::adapt::NodeRegistry::addToExistingPartsNew
"void stk::adapt::NodeRegistry::addToExistingPartsNew()

Check for adding new nodes to existing parts based on sub-entity part
ownership. This version does it in bulk and thus avoids repeats on
shared sub-dim entities. ";

%feature("docstring")
stk::adapt::NodeRegistry::getNewNodeAndOwningElement "SubDimCellData&
stk::adapt::NodeRegistry::getNewNodeAndOwningElement(SubDimCell_SDSEntityType
&subDimEntity) ";

%feature("docstring")  stk::adapt::NodeRegistry::getFromMapPtr "SubDimCellData* stk::adapt::NodeRegistry::getFromMapPtr(const
SubDimCell_SDSEntityType &subDimEntity) const ";

%feature("docstring")  stk::adapt::NodeRegistry::getFromMap "SubDimCellData& stk::adapt::NodeRegistry::getFromMap(const
SubDimCell_SDSEntityType &subDimEntity) const ";

%feature("docstring")  stk::adapt::NodeRegistry::putInMap "void
stk::adapt::NodeRegistry::putInMap(SubDimCell_SDSEntityType
&subDimEntity, SubDimCellData &data) ";

%feature("docstring")  stk::adapt::NodeRegistry::doForAllSubEntities "void
stk::adapt::NodeRegistry::doForAllSubEntities(ElementFunctionPrototype
function, const stk::mesh::Entity &element, vector< NeededEntityType >
&needed_entity_ranks)

this is a helper method that loops over all sub-dimensional entities
whose rank matches on of those in

Parameters:
-----------

needed_entity_ranks:  and registers that sub-dimensional entity as
needing a new node.

isGhost:  should be true if this element is a ghost, in which case
this will call the appropriate method to set up for

note: at this level of granularity we can do single edge refinement,
hanging nodes, etc. ";

%feature("docstring")
stk::adapt::NodeRegistry::noInline_getSubDimEntity "void
stk::adapt::NodeRegistry::noInline_getSubDimEntity(SubDimCell_SDSEntityType
&subDimEntity, const stk::mesh::Entity &element, stk::mesh::EntityRank
needed_entity_rank, unsigned iSubDimOrd)

fill

Parameters:
-----------

subDimEntity:  with the stk::mesh::EntityId's of the ordinal

iSubDimOrd:  sub-dimensional entity of

element:  of rank

needed_entity_rank:

!subDimEntity.insert(element.identifier());

!subDimEntity.insert(elem_nodes[inodes[jnode]].entity()->identifier());
";

%feature("docstring")  stk::adapt::NodeRegistry::getSubDimEntity "void
stk::adapt::NodeRegistry::getSubDimEntity(SubDimCell_SDSEntityType
&subDimEntity, const stk::mesh::Entity &element, stk::mesh::EntityRank
needed_entity_rank, unsigned iSubDimOrd)

fill

Parameters:
-----------

subDimEntity:  with the stk::mesh::EntityId's of the ordinal

iSubDimOrd:  sub-dimensional entity of

element:  of rank

needed_entity_rank:

!subDimEntity.insert(element.identifier());

!subDimEntity.insert(elem_nodes[inodes[jnode]].entity()->identifier());
";

%feature("docstring")  stk::adapt::NodeRegistry::total_size "unsigned
stk::adapt::NodeRegistry::total_size() ";

%feature("docstring")  stk::adapt::NodeRegistry::local_size "unsigned
stk::adapt::NodeRegistry::local_size() ";

%feature("docstring")  stk::adapt::NodeRegistry::isParallelRun "bool
stk::adapt::NodeRegistry::isParallelRun(unsigned size) ";

%feature("docstring")  stk::adapt::NodeRegistry::checkDB "void
stk::adapt::NodeRegistry::checkDB(std::string msg=\"\") ";

%feature("docstring")  stk::adapt::NodeRegistry::allocateBuffers "bool stk::adapt::NodeRegistry::allocateBuffers()

allocate the send/recv buffers for all-to-all communication ";

%feature("docstring")  stk::adapt::NodeRegistry::communicate "void
stk::adapt::NodeRegistry::communicate() ";

%feature("docstring")  stk::adapt::NodeRegistry::unpack "void
stk::adapt::NodeRegistry::unpack() ";

%feature("docstring")
stk::adapt::NodeRegistry::createNewNodesInParallel "void
stk::adapt::NodeRegistry::createNewNodesInParallel()

after registering all needed nodes, this method is used to request new
nodes on this processor ";

%feature("docstring")  stk::adapt::NodeRegistry::createNodeAndConnect
"void stk::adapt::NodeRegistry::createNodeAndConnect(CommDataType
&buffer_entry, NodeIdsOnSubDimEntityType &nodeIds_onSE, unsigned
from_proc, vector< stk::mesh::EntityProc > &nodes_to_ghost)

unpacks the incoming information in

Parameters:
-----------

buffer_entry:  and adds that information to my local node registry
(i.e. the map of sub-dimensional entity to global node id is updated)
";

%feature("docstring")  stk::adapt::NodeRegistry::getMap "SubDimCellToDataMap& stk::adapt::NodeRegistry::getMap() ";

%feature("docstring")  stk::adapt::NodeRegistry::getMesh "PerceptMesh& stk::adapt::NodeRegistry::getMesh() ";

%feature("docstring")  stk::adapt::NodeRegistry::getUseCustomGhosting
"bool stk::adapt::NodeRegistry::getUseCustomGhosting() ";

%feature("docstring")  stk::adapt::NodeRegistry::cleanDeletedNodes "void stk::adapt::NodeRegistry::cleanDeletedNodes(std::set<
stk::mesh::Entity * > &deleted_nodes, bool debug=false) ";

%feature("docstring")
stk::adapt::NodeRegistry::clear_element_owner_data_phase_2 "void
stk::adapt::NodeRegistry::clear_element_owner_data_phase_2() ";

%feature("docstring")
stk::adapt::NodeRegistry::clear_element_owner_data "void
stk::adapt::NodeRegistry::clear_element_owner_data(std::set<
stk::mesh::Entity * > &elems_to_be_deleted) ";

%feature("docstring")  stk::adapt::NodeRegistry::dumpDB "void
stk::adapt::NodeRegistry::dumpDB(std::string msg=\"\") ";

%feature("docstring")  stk::adapt::NodeRegistry::get_memory_usage "unsigned stk::adapt::NodeRegistry::get_memory_usage() ";


// File: classstk_1_1percept_1_1NoMallocArray.xml
%feature("docstring") stk::percept::NoMallocArray "";

%feature("docstring")  stk::percept::NoMallocArray::NoMallocArray "stk::percept::NoMallocArray< T, N >::NoMallocArray() ";

%feature("docstring")  stk::percept::NoMallocArray::NoMallocArray "stk::percept::NoMallocArray< T, N >::NoMallocArray(size_type sz, const
T &val) ";

%feature("docstring")  stk::percept::NoMallocArray::clear "void
stk::percept::NoMallocArray< T, N >::clear() ";

%feature("docstring")  stk::percept::NoMallocArray::begin "iterator
stk::percept::NoMallocArray< T, N >::begin() ";

%feature("docstring")  stk::percept::NoMallocArray::begin "const_iterator stk::percept::NoMallocArray< T, N >::begin() const ";

%feature("docstring")  stk::percept::NoMallocArray::end "iterator
stk::percept::NoMallocArray< T, N >::end() ";

%feature("docstring")  stk::percept::NoMallocArray::end "const_iterator stk::percept::NoMallocArray< T, N >::end() const ";

%feature("docstring")  stk::percept::NoMallocArray::resize "void
stk::percept::NoMallocArray< T, N >::resize(size_type n) ";

%feature("docstring")  stk::percept::NoMallocArray::rbegin "reverse_iterator stk::percept::NoMallocArray< T, N >::rbegin() ";

%feature("docstring")  stk::percept::NoMallocArray::rbegin "const_reverse_iterator stk::percept::NoMallocArray< T, N >::rbegin()
const ";

%feature("docstring")  stk::percept::NoMallocArray::rend "reverse_iterator stk::percept::NoMallocArray< T, N >::rend() ";

%feature("docstring")  stk::percept::NoMallocArray::rend "const_reverse_iterator stk::percept::NoMallocArray< T, N >::rend()
const ";

%feature("docstring")  stk::percept::NoMallocArray::at "reference
stk::percept::NoMallocArray< T, N >::at(size_type i) ";

%feature("docstring")  stk::percept::NoMallocArray::at "const_reference stk::percept::NoMallocArray< T, N >::at(size_type i)
const ";

%feature("docstring")  stk::percept::NoMallocArray::front "reference
stk::percept::NoMallocArray< T, N >::front() ";

%feature("docstring")  stk::percept::NoMallocArray::front "const_reference stk::percept::NoMallocArray< T, N >::front() const ";

%feature("docstring")  stk::percept::NoMallocArray::back "reference
stk::percept::NoMallocArray< T, N >::back() ";

%feature("docstring")  stk::percept::NoMallocArray::back "const_reference stk::percept::NoMallocArray< T, N >::back() const ";

%feature("docstring")  stk::percept::NoMallocArray::size "size_type
stk::percept::NoMallocArray< T, N >::size() const ";

%feature("docstring")  stk::percept::NoMallocArray::insert "void
stk::percept::NoMallocArray< T, N >::insert(T val) ";

%feature("docstring")  stk::percept::NoMallocArray::empty "bool
stk::percept::NoMallocArray< T, N >::empty() ";

%feature("docstring")  stk::percept::NoMallocArray::max_size "size_type stk::percept::NoMallocArray< T, N >::max_size() ";

%feature("docstring")  stk::percept::NoMallocArray::max_capacity "size_type stk::percept::NoMallocArray< T, N >::max_capacity() ";

%feature("docstring")  stk::percept::NoMallocArray::contains "bool
stk::percept::NoMallocArray< T, N >::contains(T val) ";

%feature("docstring")  stk::percept::NoMallocArray::swap "void
stk::percept::NoMallocArray< T, N >::swap(NoMallocArray< T, N > &y) ";

%feature("docstring")  stk::percept::NoMallocArray::data "const T*
stk::percept::NoMallocArray< T, N >::data() const ";

%feature("docstring")  stk::percept::NoMallocArray::data "T*
stk::percept::NoMallocArray< T, N >::data() ";

%feature("docstring")  stk::percept::NoMallocArray::c_array "T*
stk::percept::NoMallocArray< T, N >::c_array() ";

%feature("docstring")  stk::percept::NoMallocArray::assign "void
stk::percept::NoMallocArray< T, N >::assign(const T &value) ";

%feature("docstring")  stk::percept::NoMallocArray::rangecheck "void
stk::percept::NoMallocArray< T, N >::rangecheck(size_type i) ";


// File: classstk_1_1percept_1_1Norm.xml
%feature("docstring") stk::percept::Norm "

for Power = -1, compute the inf-norm

C++ includes: Norm.hpp ";

%feature("docstring")  stk::percept::Norm::Norm "stk::percept::Norm<
Power >::Norm(mesh::BulkData &bulkData, std::string partName,
TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) ";

%feature("docstring")  stk::percept::Norm::Norm "stk::percept::Norm<
Power >::Norm(mesh::BulkData &bulkData, MDArrayString &partNames,
TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) ";

%feature("docstring")  stk::percept::Norm::Norm "stk::percept::Norm<
Power >::Norm(mesh::BulkData &bulkData, mesh::Part *part=0,
TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) ";

%feature("docstring")  stk::percept::Norm::Norm "stk::percept::Norm<
Power >::Norm(mesh::BulkData &bulkData, mesh::Selector *selector,
TurboOption turboOpt=TURBO_NONE, bool is_surface_norm=false) ";

%feature("docstring")  stk::percept::Norm::~Norm "virtual
stk::percept::Norm< Power >::~Norm() ";

%feature("docstring")  stk::percept::Norm::setCubDegree "void
stk::percept::Norm< Power >::setCubDegree(unsigned cubDegree) ";

%feature("docstring")  stk::percept::Norm::getCubDegree "unsigned
stk::percept::Norm< Power >::getCubDegree() ";

%feature("docstring")  stk::percept::Norm::set_is_surface_norm "void
stk::percept::Norm< Power >::set_is_surface_norm(bool is_surface_norm)
";

%feature("docstring")  stk::percept::Norm::get_is_surface_norm "bool
stk::percept::Norm< Power >::get_is_surface_norm() ";

%feature("docstring")  stk::percept::Norm::evaluate "double
stk::percept::Norm< Power >::evaluate(Function &integrand) ";

%feature("docstring")  stk::percept::Norm::error_check_is_surface_norm
"void stk::percept::Norm< Power >::error_check_is_surface_norm()

if a Selector is specified with part(s) that are not auto-declared,
make sure all parts are of the same rank, and that m_is_surface_norm
is set correctly (if not, warn...) ";


// File: classstk_1_1percept_1_1Observable.xml
%feature("docstring") stk::percept::Observable "";

%feature("docstring")  stk::percept::Observable::Observable "stk::percept::Observable< DATA_TYPE >::Observable() ";

%feature("docstring")  stk::percept::Observable::addObserver "void
stk::percept::Observable< DATA_TYPE >::addObserver(Observer< DATA_TYPE
> &observer) ";

%feature("docstring")  stk::percept::Observable::removeObserver "void
stk::percept::Observable< DATA_TYPE >::removeObserver(Observer<
DATA_TYPE > &observer) ";

%feature("docstring")  stk::percept::Observable::getObservers "Observers& stk::percept::Observable< DATA_TYPE >::getObservers() ";

%feature("docstring")  stk::percept::Observable::notifyObservers "void stk::percept::Observable< DATA_TYPE >::notifyObservers(DATA_TYPE
*data) ";


// File: classstk_1_1percept_1_1Observer.xml
%feature("docstring") stk::percept::Observer "";

%feature("docstring")  stk::percept::Observer::Observer "stk::percept::Observer< DATA_TYPE >::Observer(Observable< DATA_TYPE >
&observable) ";

%feature("docstring")  stk::percept::Observer::~Observer "stk::percept::Observer< DATA_TYPE >::~Observer() ";

%feature("docstring")  stk::percept::Observer::notify "virtual void
stk::percept::Observer< DATA_TYPE >::notify(DATA_TYPE *data)=0 ";


// File: classOnePointRule.xml
%feature("docstring") OnePointRule "";


// File: classstk_1_1percept_1_1ParallelMachineFinalize.xml
%feature("docstring") stk::percept::ParallelMachineFinalize "";

%feature("docstring")
stk::percept::ParallelMachineFinalize::ParallelMachineFinalize "stk::percept::ParallelMachineFinalize::ParallelMachineFinalize(bool
need_to_finalize=false) ";

%feature("docstring")
stk::percept::ParallelMachineFinalize::~ParallelMachineFinalize "stk::percept::ParallelMachineFinalize::~ParallelMachineFinalize() ";


// File: classstk_1_1percept_1_1PartOp.xml
%feature("docstring") stk::percept::PartOp "";


// File: classstk_1_1percept_1_1PerceptMesh.xml
%feature("docstring") stk::percept::PerceptMesh "";

%feature("docstring")  stk::percept::PerceptMesh::PerceptMesh "stk::percept::PerceptMesh::PerceptMesh(size_t spatialDimension=3u,
stk::ParallelMachine comm=MPI_COMM_WORLD)

high-level interface

Create a Mesh object that owns its constituent FEMMetaData and
BulkData (which are created by this object) ";

%feature("docstring")  stk::percept::PerceptMesh::PerceptMesh "stk::percept::PerceptMesh::PerceptMesh(const
stk::mesh::fem::FEMMetaData *metaData, stk::mesh::BulkData *bulkData,
bool isCommitted=true)

Create a Mesh object that doesn't own its constituent FEMMetaData and
BulkData, pointers to which are adopted by this constructor. ";

%feature("docstring")  stk::percept::PerceptMesh::open_read_only "void stk::percept::PerceptMesh::open_read_only(const std::string
&in_filename)

reads and commits mesh, editing disabled ";

%feature("docstring")  stk::percept::PerceptMesh::open "void
stk::percept::PerceptMesh::open(const std::string &in_filename)

reads but doesn't commit mesh, enabling edit ";

%feature("docstring")  stk::percept::PerceptMesh::new_mesh_read_only "void stk::percept::PerceptMesh::new_mesh_read_only(const GMeshSpec
gmesh_spec)

creates a new mesh using the GeneratedMesh fixture with spec

Parameters:
-----------

gmesh_spec:  Read:  Only mode, no edits allowed ";

%feature("docstring")  stk::percept::PerceptMesh::new_mesh "void
stk::percept::PerceptMesh::new_mesh(const GMeshSpec gmesh_spec)

creates a new mesh using the GeneratedMesh fixture with spec

Parameters:
-----------

gmesh_spec:  ";

%feature("docstring")  stk::percept::PerceptMesh::add_field "stk::mesh::FieldBase* stk::percept::PerceptMesh::add_field(const
std::string &name, const unsigned entity_rank, int vectorDimension=0,
const std::string part_name=\"universal_part\")

add a field to the mesh ";

%feature("docstring")  stk::percept::PerceptMesh::get_field "stk::mesh::FieldBase * stk::percept::PerceptMesh::get_field(const
std::string &name) ";

%feature("docstring")  stk::percept::PerceptMesh::commit "void
stk::percept::PerceptMesh::commit()

commits mesh - any operations done on a non-committed mesh, except to
add fields will throw an exception ";

%feature("docstring")  stk::percept::PerceptMesh::reopen "void
stk::percept::PerceptMesh::reopen(const std::string
temp_file_name=\"percept_tmp.e\")

reopens the mesh for editing - warning, this operation writes the mesh
to a temp file then re-reads it and thus recreates the internal
FEMMetaData and BulkData ";

%feature("docstring")  stk::percept::PerceptMesh::save_as "void
stk::percept::PerceptMesh::save_as(const std::string &out_filename)

commits mesh if not committed and saves it in new file ";

%feature("docstring")  stk::percept::PerceptMesh::close "void
stk::percept::PerceptMesh::close()

closes this mesh, deleting its data

closes this mesh to further changes ";

%feature("docstring")  stk::percept::PerceptMesh::print_info "void
stk::percept::PerceptMesh::print_info(std::ostream &stream,
std::string header=\"\", int print_level=0, bool do_endl=true)

print number of parts and fields, and info on each ";

%feature("docstring")  stk::percept::PerceptMesh::print_info "void
stk::percept::PerceptMesh::print_info(std::string header=\"\", int
print_level=0, bool do_endl=true)

print number of parts and fields, and info on each ";

%feature("docstring")  stk::percept::PerceptMesh::print_fields "void
stk::percept::PerceptMesh::print_fields(std::string header=\"\")

print the fields defined on the mesh ";

%feature("docstring")  stk::percept::PerceptMesh::get_spatial_dim "int stk::percept::PerceptMesh::get_spatial_dim() ";

%feature("docstring")  stk::percept::PerceptMesh::get_number_elements
"int stk::percept::PerceptMesh::get_number_elements() ";

%feature("docstring")  stk::percept::PerceptMesh::get_number_nodes "int stk::percept::PerceptMesh::get_number_nodes() ";

%feature("docstring")  stk::percept::PerceptMesh::get_number_edges "int stk::percept::PerceptMesh::get_number_edges() ";

%feature("docstring")
stk::percept::PerceptMesh::get_number_elements_locally_owned "int
stk::percept::PerceptMesh::get_number_elements_locally_owned() ";

%feature("docstring")  stk::percept::PerceptMesh::get_rank "unsigned
stk::percept::PerceptMesh::get_rank()

parallel rank ";

%feature("docstring")  stk::percept::PerceptMesh::get_parallel_rank "unsigned stk::percept::PerceptMesh::get_parallel_rank() ";

%feature("docstring")  stk::percept::PerceptMesh::get_parallel_size "unsigned stk::percept::PerceptMesh::get_parallel_size() ";

%feature("docstring")  stk::percept::PerceptMesh::print_entity "void
stk::percept::PerceptMesh::print_entity(const stk::mesh::Entity
&entity, stk::mesh::FieldBase *field=0)

print a node, edge, element, etc; optionally pass in a field to dump
data associated with the entity ";

%feature("docstring")  stk::percept::PerceptMesh::print_entity_compact
"std::string stk::percept::PerceptMesh::print_entity_compact(const
stk::mesh::Entity &entity, stk::mesh::FieldBase *field=0)

shorter output for print_entity ";

%feature("docstring")  stk::percept::PerceptMesh::dump_elements "void
stk::percept::PerceptMesh::dump_elements(const std::string
&partName=\"\")

print elements on the given part ";

%feature("docstring")
stk::percept::PerceptMesh::dump_elements_compact "void
stk::percept::PerceptMesh::dump_elements_compact(const std::string
&partName=\"\")

compact print of elements on the given part ";

%feature("docstring")  stk::percept::PerceptMesh::get_bulk_data "stk::mesh::BulkData * stk::percept::PerceptMesh::get_bulk_data()

get the low-level bulk data pointer from stk_mesh ";

%feature("docstring")  stk::percept::PerceptMesh::get_fem_meta_data "stk::mesh::fem::FEMMetaData *
stk::percept::PerceptMesh::get_fem_meta_data()

get the low-level meta data pointer from stk_mesh ";

%feature("docstring")  stk::percept::PerceptMesh::get_part "mesh::Part* stk::percept::PerceptMesh::get_part(const std::string
&part_name)

get a pointer to a stk_mesh Part with the given name ";

%feature("docstring")  stk::percept::PerceptMesh::get_field_data "double stk::percept::PerceptMesh::get_field_data(const
stk::mesh::FieldBase *field, const mesh::Entity *entity, unsigned
ordinal=0)

get the value of a field on the given entity; if a vector field, pass
in the index of the vector required (ordinal) ";

%feature("docstring")  stk::percept::PerceptMesh::set_field_data "void stk::percept::PerceptMesh::set_field_data(double value, const
stk::mesh::FieldBase *field, const mesh::Entity *entity, unsigned
ordinal=0)

set the value of a field on the given entity; if a vector field, pass
in the index of the vector required (ordinal) ";

%feature("docstring")  stk::percept::PerceptMesh::get_node_field_data
"double
stk::percept::PerceptMesh::get_node_field_data(stk::mesh::FieldBase
*field, const mesh::EntityId node_id, unsigned ordinal=0)

get the value of a field on the given node; if a vector field, pass in
the index of the vector required (ordinal) ";

%feature("docstring")  stk::percept::PerceptMesh::set_node_field_data
"void stk::percept::PerceptMesh::set_node_field_data(double value,
stk::mesh::FieldBase *field, const mesh::EntityId node_id, unsigned
ordinal=0)

set the value of a field on the given node; if a vector field, pass in
the index of the vector required (ordinal) ";

%feature("docstring")  stk::percept::PerceptMesh::get_node "stk::mesh::Entity* stk::percept::PerceptMesh::get_node(const
mesh::EntityId node_id)

get a pointer to a node with given id ";

%feature("docstring")  stk::percept::PerceptMesh::get_element "stk::mesh::Entity* stk::percept::PerceptMesh::get_element(const
mesh::EntityId element_id)

get a pointer to an element with given id ";

%feature("docstring")  stk::percept::PerceptMesh::get_entity "stk::mesh::Entity*
stk::percept::PerceptMesh::get_entity(mesh::EntityRank rank, const
mesh::EntityId id)

get a pointer to an entity with given id ";

%feature("docstring")  stk::percept::PerceptMesh::get_node "stk::mesh::Entity * stk::percept::PerceptMesh::get_node(double x,
double y, double z=0, double t=0)

find and return pointer to node closest to given point - in parallel,
check return for null (if null, closest node is on another proc) ";

%feature("docstring")  stk::percept::PerceptMesh::get_element "stk::mesh::Entity * stk::percept::PerceptMesh::get_element(double x,
double y, double z=0, double t=0)

find and return pointer to element that contains given point - in
parallel, check return for null (if null, element containing point is
on another proc)

find element that contains or is closest to given point ";

%feature("docstring")
stk::percept::PerceptMesh::get_coordinates_field "VectorFieldType*
stk::percept::PerceptMesh::get_coordinates_field()

return a pointer to the field containing node coordinates ";

%feature("docstring")  stk::percept::PerceptMesh::node_rank "stk::mesh::EntityRank stk::percept::PerceptMesh::node_rank() const

return the rank of a node ";

%feature("docstring")  stk::percept::PerceptMesh::edge_rank "stk::mesh::EntityRank stk::percept::PerceptMesh::edge_rank() const

Returns the edge rank which changes depending on spatial dimension. ";

%feature("docstring")  stk::percept::PerceptMesh::face_rank "stk::mesh::EntityRank stk::percept::PerceptMesh::face_rank() const

Returns the face rank which changes depending on spatial dimension. ";

%feature("docstring")  stk::percept::PerceptMesh::side_rank "stk::mesh::EntityRank stk::percept::PerceptMesh::side_rank() const

Returns the side rank which changes depending on spatial dimension. ";

%feature("docstring")  stk::percept::PerceptMesh::element_rank "stk::mesh::EntityRank stk::percept::PerceptMesh::element_rank() const

Returns the element rank which is always equal to spatial dimension.
";

%feature("docstring")
stk::percept::PerceptMesh::read_database_at_step "void
stk::percept::PerceptMesh::read_database_at_step(int step)

set the current data in fields to the given Exodus step by reading
from the database ";

%feature("docstring")
stk::percept::PerceptMesh::read_database_at_time "void
stk::percept::PerceptMesh::read_database_at_time(double time)

set the current data in fields to the given Exodus time by reading
from the database (finds the closest step to the given time (no
interpolation yet)) ";

%feature("docstring")
stk::percept::PerceptMesh::get_current_database_step "int
stk::percept::PerceptMesh::get_current_database_step()

return the current state of the Exodus database, 0 if not loaded yet
(steps are 1-based in Exodus) ";

%feature("docstring")
stk::percept::PerceptMesh::get_current_database_time "double
stk::percept::PerceptMesh::get_current_database_time()

return the current state of the Exodus database (time associated with
current step) ";

%feature("docstring")
stk::percept::PerceptMesh::get_database_step_at_time "int
stk::percept::PerceptMesh::get_database_step_at_time(double time)

return the step number closest to specified time, thus
read_database_at_time(time) is equivalent to
read_database_at_step(get_database_step_at_time(time)) ";

%feature("docstring")
stk::percept::PerceptMesh::get_database_time_at_step "double
stk::percept::PerceptMesh::get_database_time_at_step(int step)

return the state time associated with given step ";

%feature("docstring")
stk::percept::PerceptMesh::get_database_time_step_count "int
stk::percept::PerceptMesh::get_database_time_step_count()

return the number of steps in the database ";

%feature("docstring")  stk::percept::PerceptMesh::transform_mesh "void stk::percept::PerceptMesh::transform_mesh(MDArray &matrix)

transform mesh by a given 3x3 matrix ";

%feature("docstring")
stk::percept::PerceptMesh::add_coordinate_state_fields "void
stk::percept::PerceptMesh::add_coordinate_state_fields()

add coordinate-like fields needed, for example, to use smoothing of
geometry-projected refined meshes Must be called before commit() ";

%feature("docstring")  stk::percept::PerceptMesh::add_spacing_fields "void stk::percept::PerceptMesh::add_spacing_fields()

add spacing fields for having refinement obey the spacing (i.e.
putting new nodes not at midpoint) ";

%feature("docstring")  stk::percept::PerceptMesh::set_proc_rank_field
"void
stk::percept::PerceptMesh::set_proc_rank_field(stk::mesh::FieldBase
*proc_rank_field=0)

set proc_rank on each element ";

%feature("docstring")
stk::percept::PerceptMesh::has_coordinate_state_fields "bool
stk::percept::PerceptMesh::has_coordinate_state_fields()

get number of coordinate field states needed ";

%feature("docstring")  stk::percept::PerceptMesh::copy_field_state "void stk::percept::PerceptMesh::copy_field_state(stk::mesh::FieldBase
*field, unsigned dest_state, unsigned src_state)

copy field state data from one state (src_state) to another
(dest_state) ";

%feature("docstring")  stk::percept::PerceptMesh::copy_field "void
stk::percept::PerceptMesh::copy_field(stk::mesh::FieldBase
*field_dest, stk::mesh::FieldBase *field_src)

copy field data from one field (field_src) to another (field_dest) ";

%feature("docstring")
stk::percept::PerceptMesh::nodal_field_state_axpby "void
stk::percept::PerceptMesh::nodal_field_state_axpby(stk::mesh::FieldBase
*field, double alpha, unsigned x_state, double beta, unsigned y_state)

axpby calculates: y = alpha*x + beta*y ";

%feature("docstring")  stk::percept::PerceptMesh::nodal_field_axpby "void stk::percept::PerceptMesh::nodal_field_axpby(double alpha,
stk::mesh::FieldBase *field_x, double beta, stk::mesh::FieldBase
*field_y)

axpby calculates: y = alpha*x + beta*y ";

%feature("docstring")
stk::percept::PerceptMesh::nodal_field_state_axpbypgz "void
stk::percept::PerceptMesh::nodal_field_state_axpbypgz(stk::mesh::FieldBase
*field, double alpha, unsigned x_state, double beta, unsigned y_state,
double gamma, unsigned z_state)

axpbypgz calculates: z = alpha*x + beta*y + gamma*z ";

%feature("docstring")  stk::percept::PerceptMesh::nodal_field_axpbypgz
"void stk::percept::PerceptMesh::nodal_field_axpbypgz(double alpha,
stk::mesh::FieldBase *field_x, double beta, stk::mesh::FieldBase
*field_y, double gamma, stk::mesh::FieldBase *field_z)

axpbypgz calculates: z = alpha*x + beta*y + gamma*z ";

%feature("docstring")  stk::percept::PerceptMesh::nodal_field_dot "double stk::percept::PerceptMesh::nodal_field_dot(stk::mesh::FieldBase
*field_x, stk::mesh::FieldBase *field_y)

dot calculates: x.y ";

%feature("docstring")
stk::percept::PerceptMesh::nodal_field_set_value "void
stk::percept::PerceptMesh::nodal_field_set_value(stk::mesh::FieldBase
*field_x, double value=0.0)

set field to constant value ";

%feature("docstring")
stk::percept::PerceptMesh::remove_geometry_blocks_on_output "void
stk::percept::PerceptMesh::remove_geometry_blocks_on_output(std::string
geometry_file_name)

remove blocks in the mesh used solely for geometry association, during
output of the mesh to Exodus.

Parameters:
-----------

geometry_file_name:  = name of the OpenNURBS file (*.3dm) containing
the geometry info

Only available when Percept is configured with
STK_PERCEPT_HAS_GEOMETRY ";

%feature("docstring")  stk::percept::PerceptMesh::set_sync_io_regions
"void stk::percept::PerceptMesh::set_sync_io_regions(bool val) ";

%feature("docstring")  stk::percept::PerceptMesh::transform_mesh "void stk::percept::PerceptMesh::transform_mesh(Math::Matrix &matrix)

transform mesh by a given 3x3 matrix ";

%feature("docstring")  stk::percept::PerceptMesh::get_non_const_part "stk::mesh::Part * stk::percept::PerceptMesh::get_non_const_part(const
std::string &part_name) ";

%feature("docstring")  stk::percept::PerceptMesh::openEmpty "void
stk::percept::PerceptMesh::openEmpty()

opens an empty mesh, with a commit ";

%feature("docstring")  stk::percept::PerceptMesh::get_closest_node "stk::mesh::Entity * stk::percept::PerceptMesh::get_closest_node(double
x, double y, double z=0, double t=0, double *sum_min_ret=0)

find node closest to given point ";

%feature("docstring")  stk::percept::PerceptMesh::~PerceptMesh "stk::percept::PerceptMesh::~PerceptMesh() ";

%feature("docstring")  stk::percept::PerceptMesh::init "void
stk::percept::PerceptMesh::init(stk::ParallelMachine
comm=MPI_COMM_WORLD, bool no_alloc=false) ";

%feature("docstring")  stk::percept::PerceptMesh::destroy "void
stk::percept::PerceptMesh::destroy() ";

%feature("docstring")  stk::percept::PerceptMesh::getPart "const
stk::mesh::Part * stk::percept::PerceptMesh::getPart(const std::string
&part_name) ";

%feature("docstring")  stk::percept::PerceptMesh::print_entity "void
stk::percept::PerceptMesh::print_entity(std::ostream &out, const
stk::mesh::Entity &entity, stk::mesh::FieldBase *field=0) ";

%feature("docstring")  stk::percept::PerceptMesh::setSpatialDim "void
stk::percept::PerceptMesh::setSpatialDim(int sd) ";

%feature("docstring")  stk::percept::PerceptMesh::setStreamingSize "void stk::percept::PerceptMesh::setStreamingSize(int streaming_size)
";

%feature("docstring")  stk::percept::PerceptMesh::getStreamingSize "int stk::percept::PerceptMesh::getStreamingSize() ";

%feature("docstring")  stk::percept::PerceptMesh::dump "void
stk::percept::PerceptMesh::dump(const std::string &file=\"\")

reads the given file into a temporary model and prints info about it

Read in the model given by.

Parameters:
-----------

file:  and print some info about the file to stdout ";

%feature("docstring")  stk::percept::PerceptMesh::isGhostElement "bool stk::percept::PerceptMesh::isGhostElement(const stk::mesh::Entity
&element) ";

%feature("docstring")
stk::percept::PerceptMesh::check_entity_duplicate "bool
stk::percept::PerceptMesh::check_entity_duplicate(stk::mesh::Entity
&entity) ";

%feature("docstring")  stk::percept::PerceptMesh::delete_side_sets "void stk::percept::PerceptMesh::delete_side_sets() ";

%feature("docstring")
stk::percept::PerceptMesh::addParallelInfoFields "void
stk::percept::PerceptMesh::addParallelInfoFields(bool elemental, bool
nodal, std::string elemental_proc_rank_name=\"proc_rank\", std::string
nodal_fixed_flag=\"fixed\", std::string
nodal_global_id_name=\"GLOBAL_ID\", std::string
nodal_proc_id_name=\"PROCESSOR_ID\", std::string
nodal_local_id_name=\"LOCAL_ID\")

add some fields that are useful for debugging or for exporting meshes
to Mesquite - must be done before commit() ";

%feature("docstring")
stk::percept::PerceptMesh::populateParallelInfoFields "void
stk::percept::PerceptMesh::populateParallelInfoFields(bool elemental,
bool nodal, stk::mesh::Selector *fixed_node_selector=0, std::string
elemental_proc_rank_name=\"proc_rank\", std::string
nodal_fixed_flag=\"fixed\", std::string
nodal_global_id_name=\"GLOBAL_ID\", std::string
nodal_proc_id_name=\"PROCESSOR_ID\", std::string
nodal_local_id_name=\"LOCAL_ID\")

fill the fields from addParallelInfoFields with data from stk_mesh
database ";

%feature("docstring")
stk::percept::PerceptMesh::getFamilyTreeRelationIndex "unsigned
stk::percept::PerceptMesh::getFamilyTreeRelationIndex(FamiltyTreeLevel
level, const stk::mesh::Entity &element)

A family tree relation holds the parent/child relations for a refined
mesh.

Case 0: a single refinement of a parent P_0 and its children C_0_0,
C_0_1,...,C_0_N leads to a new family tree entity FT_0 that has down
relations to {P_0, C_0_0, C_0_1,...,C_0_N} The back pointers from P_0,
C_0_0, ... are initially stored as the 0'th index of their relations,
i.e.: P_0.relations(FAMILY_TREE_RANK)[0] --> FT_0,
C_0_0.relations(FAMILY_TREE_RANK)[0] --> FT_0, etc. Case 1: a
previously refined child, say C_0_1, renamed to P_0_1, gets further
refined leading to a new family tree entity, FT_1 pointing to: {P_0_1,
C_0_1_0, C_0_1_1,... } but, now the relations indexing changes
(actually, we can't predict what it will be, thus the need for this
function getFamilyTreeRelationIndex):
P_0_1.relations(FAMILY_TREE_RANK)[0] --> FT_1
P_0_1.relations(FAMILY_TREE_RANK)[1] --> FT_0 etc. So, we use this
function to look for the family tree corresponding to if we are
looking for the first level (if there's only one level, or we are
looking for the family tree associated with the element when it was a
child for the first time), orthe \"level 1\" family tree
(corresponding to Case 1 where we are looking for the family tree of
the element associated with it being a parent). ";

%feature("docstring")  stk::percept::PerceptMesh::isChildElement "bool stk::percept::PerceptMesh::isChildElement(const stk::mesh::Entity
&element, bool check_for_family_tree=true)

the element is not a parent of the 0'th family_tree relation ";

%feature("docstring")  stk::percept::PerceptMesh::isLeafElement "bool
stk::percept::PerceptMesh::isLeafElement(const stk::mesh::Entity
&element) ";

%feature("docstring")  stk::percept::PerceptMesh::isChildElementLeaf "bool stk::percept::PerceptMesh::isChildElementLeaf(const
stk::mesh::Entity &element, bool check_for_family_tree=true)

the element is not a parent of any family tree relation ";

%feature("docstring")  stk::percept::PerceptMesh::hasFamilyTree "bool
stk::percept::PerceptMesh::hasFamilyTree(const stk::mesh::Entity
&element) ";

%feature("docstring")  stk::percept::PerceptMesh::isParentElement "bool stk::percept::PerceptMesh::isParentElement(const
stk::mesh::Entity &element, bool check_for_family_tree=true)

if the element is a parent at any level, return true ";

%feature("docstring")  stk::percept::PerceptMesh::isParentElementLeaf
"bool stk::percept::PerceptMesh::isParentElementLeaf(const
stk::mesh::Entity &element, bool check_for_family_tree=true)

is element a parent at the leaf level (either there is only one level,
and it's a parent, or if more than one, the element is a child and a
parent and its children have no children) ";

%feature("docstring")
stk::percept::PerceptMesh::isParentElementLevel2 "bool
stk::percept::PerceptMesh::isParentElementLevel2(const
stk::mesh::Entity &element, bool check_for_family_tree=true)

is element a parent at level 2 (meaning that it is both a child and a
parent) ";

%feature("docstring")  stk::percept::PerceptMesh::isChildWithoutNieces
"bool stk::percept::PerceptMesh::isChildWithoutNieces(const
stk::mesh::Entity &element, bool check_for_family_tree=true)

is element a child with siblings with no nieces or nephews (siblings
with children) (alternative would be \"is child and is parent not a
grandparent\") ";

%feature("docstring")  stk::percept::PerceptMesh::getChildren "bool
stk::percept::PerceptMesh::getChildren(const stk::mesh::Entity
&element, std::vector< stk::mesh::Entity * > &children, bool
check_for_family_tree=true, bool only_if_element_is_parent_leaf=false)
";

%feature("docstring")  stk::percept::PerceptMesh::printParentChildInfo
"void stk::percept::PerceptMesh::printParentChildInfo(const
stk::mesh::Entity &element, bool check_for_family_tree=true) ";

%feature("docstring")  stk::percept::PerceptMesh::createOrGetNode "stk::mesh::Entity &
stk::percept::PerceptMesh::createOrGetNode(stk::mesh::EntityId nid,
double *x=0) ";

%feature("docstring")  stk::percept::PerceptMesh::createEntities "void stk::percept::PerceptMesh::createEntities(stk::mesh::EntityRank
entityRank, int count, std::vector< stk::mesh::Entity * >
&requested_entities) ";

%feature("docstring")  stk::percept::PerceptMesh::node_field_data "double*
stk::percept::PerceptMesh::node_field_data(stk::mesh::FieldBase
*field, const mesh::EntityId node_id) ";

%feature("docstring")  stk::percept::PerceptMesh::nodalOpLoop "void
stk::percept::PerceptMesh::nodalOpLoop(GenericFunction &nodalOp,
stk::mesh::FieldBase *field=0, stk::mesh::Selector *selector=0) ";

%feature("docstring")  stk::percept::PerceptMesh::elementOpLoop "void
stk::percept::PerceptMesh::elementOpLoop(ElementOp &elementOp,
stk::mesh::FieldBase *field=0, stk::mesh::Part *part=0)

Loop over all elements and apply.

Parameters:
-----------

elementOp:  passing in the argument

field:  to

elementOp:  ";

%feature("docstring")  stk::percept::PerceptMesh::bucketOpLoop "void
stk::percept::PerceptMesh::bucketOpLoop(BucketOp &bucketOp,
stk::mesh::FieldBase *field=0, stk::mesh::Part *part=0)

Loop over all buckets and apply.

Parameters:
-----------

bucketOp:  passing in the argument

field:  to

bucketOp:  ";

%feature("docstring")  stk::percept::PerceptMesh::elementOpLoop "void
stk::percept::PerceptMesh::elementOpLoop(ElementOp &elementOp,
stk::mesh::FieldBase *field, stk::mesh::Selector *selector, bool
is_surface_norm=false)

Loop over all elements and apply.

Parameters:
-----------

elementOp:  passing in the argument

field:  to

elementOp:  ";

%feature("docstring")  stk::percept::PerceptMesh::bucketOpLoop "void
stk::percept::PerceptMesh::bucketOpLoop(BucketOp &bucketOp,
stk::mesh::FieldBase *field, stk::mesh::Selector *selector, bool
is_surface_norm=false) ";

%feature("docstring")  stk::percept::PerceptMesh::edge_length_ave "double stk::percept::PerceptMesh::edge_length_ave(const
stk::mesh::Entity &entity, mesh::FieldBase *coord_field=0) ";

%feature("docstring")
stk::percept::PerceptMesh::adapt_parent_to_child_relations "SameRankRelation&
stk::percept::PerceptMesh::adapt_parent_to_child_relations() ";

%feature("docstring")  stk::percept::PerceptMesh::isBoundarySurface "bool stk::percept::PerceptMesh::isBoundarySurface(mesh::Part &block,
mesh::Part &surface) ";

%feature("docstring")  stk::percept::PerceptMesh::get_io_omitted_parts
"const stk::mesh::PartVector&
stk::percept::PerceptMesh::get_io_omitted_parts() ";

%feature("docstring")  stk::percept::PerceptMesh::set_io_omitted_parts
"void
stk::percept::PerceptMesh::set_io_omitted_parts(stk::mesh::PartVector
&io_omitted_parts) ";

%feature("docstring")  stk::percept::PerceptMesh::get_cell_topology "const CellTopologyData*
stk::percept::PerceptMesh::get_cell_topology(const stk::mesh::Part
&part) ";

%feature("docstring")  stk::percept::PerceptMesh::get_cell_topology "const CellTopologyData*
stk::percept::PerceptMesh::get_cell_topology(const mesh::Part &part)
";

%feature("docstring")  stk::percept::PerceptMesh::fillCellNodes "void
stk::percept::PerceptMesh::fillCellNodes(const mesh::Bucket &bucket,
mesh::FieldBase *field, ArrayType &cellNodes, unsigned dataStrideArg)
";


// File: classstk_1_1percept_1_1IntrepidManager_1_1PhysicalCoords.xml
%feature("docstring") stk::percept::IntrepidManager::PhysicalCoords "

([C], [P], [D])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::PhysicalCoords::PhysicalCoords "stk::percept::IntrepidManager::PhysicalCoords::PhysicalCoords(IM &im)

([C], [P], [D]) ";

%feature("docstring")
stk::percept::IntrepidManager::PhysicalCoords::copyTo "void
stk::percept::IntrepidManager::PhysicalCoords::copyTo(MDArray &mda) ";


// File: structstk_1_1adapt_1_1regression__tests_1_1PlaneShock.xml
%feature("docstring") stk::adapt::regression_tests::PlaneShock "";

%feature("docstring")
stk::adapt::regression_tests::PlaneShock::PlaneShock "stk::adapt::regression_tests::PlaneShock::PlaneShock() ";

%feature("docstring")
stk::adapt::regression_tests::PlaneShock::setCurrentPlanePoint "void
stk::adapt::regression_tests::PlaneShock::setCurrentPlanePoint(double
shock_displacement) ";


// File: classstk_1_1adapt_1_1PredicateBasedEdgeAdapter.xml
%feature("docstring") stk::adapt::PredicateBasedEdgeAdapter "";

%feature("docstring")
stk::adapt::PredicateBasedEdgeAdapter::PredicateBasedEdgeAdapter "stk::adapt::PredicateBasedEdgeAdapter< RefinePredicate
>::PredicateBasedEdgeAdapter(RefinePredicate &predicate_refine,
percept::PerceptMesh &eMesh, UniformRefinerPatternBase &bp,
stk::mesh::FieldBase *proc_rank_field=0) ";

%feature("docstring")
stk::adapt::PredicateBasedEdgeAdapter::getRefinePredicate "RefinePredicate& stk::adapt::PredicateBasedEdgeAdapter<
RefinePredicate >::getRefinePredicate() ";

%feature("docstring")  stk::adapt::PredicateBasedEdgeAdapter::mark "virtual int stk::adapt::PredicateBasedEdgeAdapter< RefinePredicate
>::mark(const stk::mesh::Entity &element, unsigned which_edge,
stk::mesh::Entity &node0, stk::mesh::Entity &node1, double *coord0,
double *coord1, std::vector< int > *existing_edge_marks)

DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE. ";


// File: classstk_1_1adapt_1_1PredicateBasedElementAdapter.xml
%feature("docstring") stk::adapt::PredicateBasedElementAdapter "";

%feature("docstring")
stk::adapt::PredicateBasedElementAdapter::PredicateBasedElementAdapter
"stk::adapt::PredicateBasedElementAdapter< RefinePredicate
>::PredicateBasedElementAdapter(RefinePredicate &predicate_refine,
percept::PerceptMesh &eMesh, UniformRefinerPatternBase &bp,
stk::mesh::FieldBase *proc_rank_field=0) ";

%feature("docstring")
stk::adapt::PredicateBasedElementAdapter::buildUnrefineList "virtual
ElementUnrefineCollection stk::adapt::PredicateBasedElementAdapter<
RefinePredicate >::buildUnrefineList() ";


// File: classstk_1_1percept_1_1PrintFieldOp.xml
%feature("docstring") stk::percept::PrintFieldOp "";

%feature("docstring")  stk::percept::PrintFieldOp::PrintFieldOp "stk::percept::PrintFieldOp::PrintFieldOp(std::string name, PerceptMesh
&eMesh, int dom, int codom) ";


// File: classstk_1_1percept_1_1ProgressMeter.xml
%feature("docstring") stk::percept::ProgressMeter "";

%feature("docstring")  stk::percept::ProgressMeter::ProgressMeter "stk::percept::ProgressMeter::ProgressMeter(Observable<
ProgressMeterData > &observable) ";

%feature("docstring")  stk::percept::ProgressMeter::notify "void
stk::percept::ProgressMeter::notify(ProgressMeterData *data) ";


// File: structstk_1_1percept_1_1ProgressMeterData.xml
%feature("docstring") stk::percept::ProgressMeterData "";

%feature("docstring")
stk::percept::ProgressMeterData::ProgressMeterData "stk::percept::ProgressMeterData::ProgressMeterData(STATE state, double
data, std::string stage=\"\") ";


// File: classstk_1_1percept_1_1PyramidFixture.xml
%feature("docstring") stk::percept::PyramidFixture "

Use case with mixed element topologies and field relations to provide
fast access to node field data from an element.

copied from stk_mesh and modified

C++ includes: PyramidFixture.hpp ";

%feature("docstring")  stk::percept::PyramidFixture::~PyramidFixture "stk::percept::PyramidFixture::~PyramidFixture() ";

%feature("docstring")  stk::percept::PyramidFixture::PyramidFixture "stk::percept::PyramidFixture::PyramidFixture(stk::ParallelMachine
comm, bool doCommit=true, bool do_sidesets=false) ";

%feature("docstring")  stk::percept::PyramidFixture::populate "void
stk::percept::PyramidFixture::populate() ";


// File: classstk_1_1percept_1_1QuadFixture.xml
%feature("docstring") stk::percept::QuadFixture "

Topology can also be Triangle<3>

C++ includes: QuadFixture.hpp ";

%feature("docstring")  stk::percept::QuadFixture::~QuadFixture "stk::percept::QuadFixture< Scalar, Topology >::~QuadFixture() ";

%feature("docstring")  stk::percept::QuadFixture::QuadFixture "stk::percept::QuadFixture< Scalar, Topology
>::QuadFixture(stk::ParallelMachine pm, unsigned nx, unsigned ny, bool
generate_sidesets_in, bool debug_geom_side_sets_as_blocks_in=false) ";

%feature("docstring")  stk::percept::QuadFixture::set_bounding_box "void stk::percept::QuadFixture< Scalar, Topology
>::set_bounding_box(double xmin, double xmax, double ymin, double
ymax) ";

%feature("docstring")  stk::percept::QuadFixture::generate_mesh "void
stk::percept::QuadFixture< Scalar, Topology >::generate_mesh() ";

%feature("docstring")  stk::percept::QuadFixture::generate_mesh "void
stk::percept::QuadFixture< Scalar, Topology
>::generate_mesh(std::vector< stk::mesh::EntityId >
&element_ids_on_this_processor) ";

%feature("docstring")  stk::percept::QuadFixture::node_id "stk::mesh::EntityId stk::percept::QuadFixture< Scalar, Topology
>::node_id(unsigned ix, unsigned iy) const ";

%feature("docstring")  stk::percept::QuadFixture::elem_id "stk::mesh::EntityId stk::percept::QuadFixture< Scalar, Topology
>::elem_id(unsigned ix, unsigned iy) const ";

%feature("docstring")  stk::percept::QuadFixture::node "stk::mesh::Entity* stk::percept::QuadFixture< Scalar, Topology
>::node(unsigned ix, unsigned iy) const ";

%feature("docstring")  stk::percept::QuadFixture::node_ix_iy "void
stk::percept::QuadFixture< Scalar, Topology
>::node_ix_iy(stk::mesh::EntityId entity_id, unsigned &ix, unsigned
&iy) const ";

%feature("docstring")  stk::percept::QuadFixture::elem_ix_iy "void
stk::percept::QuadFixture< Scalar, Topology
>::elem_ix_iy(stk::mesh::EntityId entity_id, unsigned &ix, unsigned
&iy) const ";

%feature("docstring")  stk::percept::QuadFixture::elem "stk::mesh::Entity* stk::percept::QuadFixture< Scalar, Topology
>::elem(unsigned ix, unsigned iy) const ";

%feature("docstring")  stk::percept::QuadFixture::generate_sides_meta
"void stk::percept::QuadFixture< Scalar, Topology
>::generate_sides_meta() ";

%feature("docstring")  stk::percept::QuadFixture::generate_sides_bulk
"void stk::percept::QuadFixture< Scalar, Topology
>::generate_sides_bulk(std::vector< stk::mesh::EntityId >
&element_ids_on_this_processor) ";


// File: structQuadrilateral4.xml
%feature("docstring") Quadrilateral4 "";


// File: structstk_1_1adapt_1_1RefinementInfoByType.xml
%feature("docstring") stk::adapt::RefinementInfoByType "";

%feature("docstring")
stk::adapt::RefinementInfoByType::RefinementInfoByType "stk::adapt::RefinementInfoByType::RefinementInfoByType() ";


// File: classstk_1_1adapt_1_1Elem_1_1RefinementKey.xml
%feature("docstring") stk::adapt::Elem::RefinementKey "";

%feature("docstring")  stk::adapt::Elem::RefinementKey::RefinementKey
"stk::adapt::Elem::RefinementKey::RefinementKey()

Default Constructer that assigns the value_ to 0x0 ";

%feature("docstring")  stk::adapt::Elem::RefinementKey::RefinementKey
"stk::adapt::Elem::RefinementKey::RefinementKey(UInt valueArg)

Constructer where value is known ";

%feature("docstring")  stk::adapt::Elem::RefinementKey::~RefinementKey
"stk::adapt::Elem::RefinementKey::~RefinementKey()

Destructor ";

%feature("docstring")  stk::adapt::Elem::RefinementKey::value "UInt
stk::adapt::Elem::RefinementKey::value()

Returns the value_ ";

%feature("docstring")  stk::adapt::Elem::RefinementKey::assign_value "void stk::adapt::Elem::RefinementKey::assign_value(UInt valueArg)

Assigns value_ ";

%feature("docstring")  stk::adapt::Elem::RefinementKey::assign_value "void stk::adapt::Elem::RefinementKey::assign_value(std::vector< UInt >
&edge_order)

function that assigns the value_ based on edges to be cut ";

%feature("docstring")
stk::adapt::Elem::RefinementKey::ordered_cut_edges "std::vector< UInt
> stk::adapt::Elem::RefinementKey::ordered_cut_edges(UInt numEdges)
const

function that returns ordered edges to be cut for this method need the
objTopologyes number of edges ";

%feature("docstring")
stk::adapt::Elem::RefinementKey::full_refinement "bool
stk::adapt::Elem::RefinementKey::full_refinement(UInt numEdges)

Function returns boolean for full refinement ";


// File: classstk_1_1adapt_1_1Elem_1_1RefinementTopology.xml
%feature("docstring") stk::adapt::Elem::RefinementTopology "";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::RefinementTopology "stk::adapt::Elem::RefinementTopology::RefinementTopology(MeshObjTopology
*mesh_obj_topology, UInt num_child, const MeshObjTopology *const
*child_topology, UInt num_child_nodes, const UInt *const *child_nodes,
const UInt num_edges, const UInt *const *edge_node, const UInt
num_faces, const UInt *const *face_node, UInt num_orientations, const
UInt *const *perm_node, const UInt *const *perm_edge, bool
homogeneous_child) ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::RefinementTopology "stk::adapt::Elem::RefinementTopology::RefinementTopology(const
CellTopology &mesh_obj_topology, UInt num_child, const MeshObjTopology
*const *child_topology, UInt num_child_nodes, const UInt *const
*child_nodes, const UInt num_edges, const UInt *const *edge_node,
const UInt num_faces, const UInt *const *face_node, UInt
num_orientations, const UInt *const *perm_node, const UInt *const
*perm_edge, bool homogeneous_child) ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::~RefinementTopology "stk::adapt::Elem::RefinementTopology::~RefinementTopology() ";

%feature("docstring")  stk::adapt::Elem::RefinementTopology::num_child
"UInt stk::adapt::Elem::RefinementTopology::num_child() const

Number of refinement child topologies. ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::homogeneous_child "bool
stk::adapt::Elem::RefinementTopology::homogeneous_child() const

Query if the refined mesh object topologies are homogeneous. ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::child_cell_topology "const
CellTopology*
stk::adapt::Elem::RefinementTopology::child_cell_topology() const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::child_cell_topology "CellTopology
stk::adapt::Elem::RefinementTopology::child_cell_topology(UInt
ordinal) const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::num_child_nodes "UInt
stk::adapt::Elem::RefinementTopology::num_child_nodes() const

Total number of unique nodes of the connected child objects. Nodes
that are shared by child objects are only counted once. Every node of
the parent is also one of the child nodes. ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::child_nodes "const UInt* const*
stk::adapt::Elem::RefinementTopology::child_nodes() const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::child_node "const UInt*
stk::adapt::Elem::RefinementTopology::child_node(UInt child) const

Map ordinals of child object topology nodes to child nodes. Array
dimension is [ child_topology(child)->num_nodes() ] child <
num_child()

node_of_child < child_topology(child)->num_nodes()

0 <= child_node(child)[i] < num_child_nodes ";

%feature("docstring")  stk::adapt::Elem::RefinementTopology::edge_node
"const UInt* stk::adapt::Elem::RefinementTopology::edge_node(UInt
edge) const ";

%feature("docstring")  stk::adapt::Elem::RefinementTopology::face_node
"const UInt* stk::adapt::Elem::RefinementTopology::face_node(UInt
face) const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::node_permutation "const UInt *
stk::adapt::Elem::RefinementTopology::node_permutation(UInt
orientation) const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::edge_permutation "const UInt *
stk::adapt::Elem::RefinementTopology::edge_permutation(UInt
orientation) const

Permutation vector from the expected orientation to the input
permutation. ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::query_refinement_topology "bool
stk::adapt::Elem::RefinementTopology::query_refinement_topology(RefinementKey
&object_key, MeshObjRefinementTopology &refTop) const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::refine_rivara_tri "bool
stk::adapt::Elem::RefinementTopology::refine_rivara_tri(RefinementKey
&, MeshObjRefinementTopology &refTop) const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::refine_rivara_tet "bool
stk::adapt::Elem::RefinementTopology::refine_rivara_tet(RefinementKey
&, MeshObjRefinementTopology &refTop) const ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::child_face "std::pair< UInt,
UInt > stk::adapt::Elem::RefinementTopology::child_face(const UInt
face_ordinal, const UInt face_child_ordinal) const

Mapping of parent->face->child to parent->child->face face_ordinal <
num_faces()

face_child_ordinal < face_topology(face_ordinal)-> num_child()

(child_ordinal, child_face_ordinal) ";

%feature("docstring")
stk::adapt::Elem::RefinementTopology::child_edge "std::pair< UInt,
UInt > stk::adapt::Elem::RefinementTopology::child_edge(const UInt
edge_ordinal, const UInt edge_child_ordinal) const

Mapping of parent->edge->child to parent->child->edge edge_ordinal <
getEdgeCount()

edge_child_ordinal < edge_topology(edge_ordinal)-> num_child()

(child_ordinal, child_edge_ordinal) ";


// File: classstk_1_1adapt_1_1Elem_1_1StdMeshObjTopologies_1_1RefinementTopologyExtra.xml
%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra "";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Line< 2 > >::refinement_topology 0 1 PARENT Linear Edge
Element Nodes (SPACE_DIM = 1!) o--------------- o

After refinement:

0 2 1 CHILD Linear Edge Element Nodes (new nodes = *)
o-------*-------o

| CHILD Linear Edge Node Maps (global node numbers!) 0 1 | o-------o |
E#1 | Element (or edge) 0: childNodeMap[0] = { 0, 2 }; | 0 1 | o
-------o | E#2 | Element (or edge) 1: childNodeMap[1] = { 2, 1 };

Refined Linear Edge (or Linear Bar element) PERMUTATION Node Maps:

Polarity = 1 { 0, 1; 2 } Polarity = 0 { 1, 0; 2 } 0 2 1 PARENT 3-Node
Line Object Nodes o-------o-------o

After refinement:

0 3 2 4 1 CHILD Objects (new nodes = *) o---*---o---*---o

| CHILD Line Node Maps (global node numbers!) 0 2 1 | o---o---o | E#1
| Object (or line) 0: childNodeMap[0] = { 0, 2, 3 }; | 0 2 1 | o---o
---o | E#2 | Object (or line) 1: childNodeMap[1] = { 2, 1, 4 };

Refined 3-Node Line Object PERMUTATION Node Maps:

Polarity = 1 { 0, 1, 2; 3, 4 } Polarity = 0 { 1, 0, 2; 4, 3 } New ref
topo info ------------------

{Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc,
param-coord}

struct RefinementTopologyExtraEntry { unsigned ordinal_of_node; //
ordinal of node in the total list of nodes - corresponds to the shards
node ordinal unsigned rank_of_subcell; // rank of the subcell this
node is associated with unsigned ordinal_of_subcell; // ordinal of the
subcell in the shards numbering (e.g. edge # 3) unsigned
ordinal_of_node_on_subcell; // ordinal of the node on the subcell
(whcih node it is on a subcell that has multiple nodes) unsigned
num_nodes_on_subcell; // how many nodes exist on the subcell double
parametric_coordinates[3]; };

Bootstrapping this file: to create this file, run the regression test
RegressionTestUniformRefiner.cpp :: generate_tables after putting in a
dummy entry in ./sierra_element/GeneratedRefinementTable.hpp. The run
will produce a local file, generated_refinement_tables.hpp which can
be checked against the gold copy of GeneratedRefinementTable.hpp, then
copied over it. Add a call below to generate the actual new table
data. ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Beam< 2 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::ShellLine< 2 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::ShellLine< 3 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Quadrilateral< 4 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Triangle< 3 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::ShellTriangle< 3 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::ShellTriangle< 6 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::ShellQuadrilateral< 4 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::ShellQuadrilateral< 8 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Tetrahedron< 4 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Hexahedron< 8 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Wedge< 6 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Wedge< 18 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Wedge< 15 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Pyramid< 5 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Pyramid< 13 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Line< 3 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Beam< 3 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Triangle< 6 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Quadrilateral< 8 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Quadrilateral< 9 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Hexahedron< 27 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Hexahedron< 20 > >::refinement_topology";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra::refinement_topology
"RefTopoX
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtra<
shards::Tetrahedron< 10 > >::refinement_topology";


// File: structstk_1_1adapt_1_1Elem_1_1StdMeshObjTopologies_1_1RefinementTopologyExtraEntry.xml
%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefinementTopologyExtraEntry "";


// File: classRefinePredicate.xml
%feature("docstring") RefinePredicate "

Predicate-based marker

The functor an operator() that returns an entry from AdaptInstruction,
either to do nothing, refine, unrefine, or both refine & unrefine
(useful for unit testing, etc.) ";


// File: classRefinePredicate.xml
%feature("docstring") RefinePredicate "

Predicate-based marker

The functor an operator() that returns an entry from AdaptInstruction,
either to do nothing, refine, unrefine, or both refine & unrefine
(useful for unit testing, etc.) ";


// File: classstk_1_1adapt_1_1Refiner.xml
%feature("docstring") stk::adapt::Refiner "

e.g. UniformRefiner<shards::Hex<8>, shards::Tet<4> >

C++ includes: Refiner.hpp ";

%feature("docstring")  stk::adapt::Refiner::Refiner "stk::adapt::Refiner::Refiner(percept::PerceptMesh &eMesh,
UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")  stk::adapt::Refiner::Refiner "stk::adapt::Refiner::Refiner(percept::PerceptMesh &eMesh, Pattern
refine_type, stk::mesh::FieldBase *proc_rank_field=0) ";

%feature("docstring")  stk::adapt::Refiner::~Refiner "stk::adapt::Refiner::~Refiner() ";

%feature("docstring")  stk::adapt::Refiner::doBreak "void
stk::adapt::Refiner::doBreak()

FIXME add part info

communicate all-to-all the new node creation information which also
updates the node registry so it can be queried locally now for any
ghost or non-ghost element

Global element ops: here's where we e.g. connect the new elements by
declaring new relations
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Global node loop operations: this is where we perform ops like adding
new nodes to the right parts, interpolating fields, etc.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
";

%feature("docstring")  stk::adapt::Refiner::setRemoveOldElements "void stk::adapt::Refiner::setRemoveOldElements(bool do_remove) ";

%feature("docstring")  stk::adapt::Refiner::getRemoveOldElements "bool stk::adapt::Refiner::getRemoveOldElements() ";

%feature("docstring")  stk::adapt::Refiner::setGeometryFile "void
stk::adapt::Refiner::setGeometryFile(std::string file_name) ";

%feature("docstring")  stk::adapt::Refiner::setSmoothGeometry "void
stk::adapt::Refiner::setSmoothGeometry(bool do_smooth) ";

%feature("docstring")  stk::adapt::Refiner::getSmoothGeometry "bool
stk::adapt::Refiner::getSmoothGeometry() ";

%feature("docstring")  stk::adapt::Refiner::setRemoveGeometryBlocks "void stk::adapt::Refiner::setRemoveGeometryBlocks(bool do_remove) ";

%feature("docstring")  stk::adapt::Refiner::getRemoveGeometryBlocks "bool stk::adapt::Refiner::getRemoveGeometryBlocks() ";

%feature("docstring")  stk::adapt::Refiner::setIgnoreSideSets "void
stk::adapt::Refiner::setIgnoreSideSets(bool ignore_sidesets) ";

%feature("docstring")  stk::adapt::Refiner::getIgnoreSideSets "bool
stk::adapt::Refiner::getIgnoreSideSets() ";

%feature("docstring")  stk::adapt::Refiner::getRefinementInfoByType "std::vector< RefinementInfoByType > &
stk::adapt::Refiner::getRefinementInfoByType() ";

%feature("docstring")  stk::adapt::Refiner::setQueryPassOnly "void
stk::adapt::Refiner::setQueryPassOnly(bool doQueryOnly) ";

%feature("docstring")  stk::adapt::Refiner::setDoProgressMeter "void
stk::adapt::Refiner::setDoProgressMeter(bool do_progress) ";

%feature("docstring")  stk::adapt::Refiner::getDoProgressMeter "bool
stk::adapt::Refiner::getDoProgressMeter() ";

%feature("docstring")  stk::adapt::Refiner::unrefineTheseElements "void
stk::adapt::Refiner::unrefineTheseElements(ElementUnrefineCollection
&elements_to_unref) ";

%feature("docstring")  stk::adapt::Refiner::unrefineAll "void
stk::adapt::Refiner::unrefineAll() ";

%feature("docstring")
stk::adapt::Refiner::setAlwaysInitializeNodeRegistry "void
stk::adapt::Refiner::setAlwaysInitializeNodeRegistry(bool do_init) ";

%feature("docstring")
stk::adapt::Refiner::getAlwaysInitializeNodeRegistry "bool
stk::adapt::Refiner::getAlwaysInitializeNodeRegistry() ";

%feature("docstring")  stk::adapt::Refiner::deleteParentElements "void stk::adapt::Refiner::deleteParentElements()

Delete all elements that aren't child elements. ";

%feature("docstring")  stk::adapt::Refiner::check_db "void
stk::adapt::Refiner::check_db(std::string msg=\"\") ";

%feature("docstring")  stk::adapt::Refiner::check_sidesets "void
stk::adapt::Refiner::check_sidesets(std::string msg=\"\") ";

%feature("docstring")  stk::adapt::Refiner::check_sidesets_1 "void
stk::adapt::Refiner::check_sidesets_1(std::string msg) ";

%feature("docstring")  stk::adapt::Refiner::check_sidesets_2 "void
stk::adapt::Refiner::check_sidesets_2(std::string msg) ";

%feature("docstring")  stk::adapt::Refiner::fix_side_sets_1 "void
stk::adapt::Refiner::fix_side_sets_1() ";

%feature("docstring")  stk::adapt::Refiner::fix_side_sets_2 "void
stk::adapt::Refiner::fix_side_sets_2() ";

%feature("docstring")  stk::adapt::Refiner::fix_side_sets_3 "void
stk::adapt::Refiner::fix_side_sets_3(bool checkParentChild,
SidePartMap &side_part_map) ";

%feature("docstring")  stk::adapt::Refiner::get_side_part_relations "void stk::adapt::Refiner::get_side_part_relations(bool
checkParentChild, SidePartMap &side_part_map)

determine side part to elem part relations ";

%feature("docstring")  stk::adapt::Refiner::connectSides "bool
stk::adapt::Refiner::connectSides(stk::mesh::Entity *element,
stk::mesh::Entity *side_elem, SidePartMap *side_part_map=0) ";

%feature("docstring")  stk::adapt::Refiner::fixElementSides2 "void
stk::adapt::Refiner::fixElementSides2() ";

%feature("docstring")  stk::adapt::Refiner::fixSides "void
stk::adapt::Refiner::fixSides(stk::mesh::Entity *parent) ";

%feature("docstring")  stk::adapt::Refiner::getNodeRegistry "NodeRegistry& stk::adapt::Refiner::getNodeRegistry() ";

%feature("docstring")  stk::adapt::Refiner::getMesh "percept::PerceptMesh& stk::adapt::Refiner::getMesh() ";


// File: classstk_1_1adapt_1_1RefinerPattern.xml
%feature("docstring") stk::adapt::RefinerPattern "";


// File: classstk_1_1adapt_1_1RefinerPattern_3_01shards_1_1Line_3_012_01_4_00_01shards_1_1Line_3_012_01_4_00-1_01_4.xml
%feature("docstring") stk::adapt::RefinerPattern< shards::Line< 2 >,
shards::Line< 2 >,-1 > " ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Line< 2 >,
shards::Line< 2 >,-1 >::RefinerPattern " stk::adapt::RefinerPattern<
shards::Line< 2 >, shards::Line< 2 >,-1
>::RefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Line< 2 >,
shards::Line< 2 >,-1 >::doBreak " virtual void
stk::adapt::RefinerPattern< shards::Line< 2 >, shards::Line< 2 >,-1
>::doBreak() ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Line< 2 >,
shards::Line< 2 >,-1 >::fillNeededEntities " void
stk::adapt::RefinerPattern< shards::Line< 2 >, shards::Line< 2 >,-1
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Line< 2 >,
shards::Line< 2 >,-1 >::getNumNewElemPerElem " virtual unsigned
stk::adapt::RefinerPattern< shards::Line< 2 >, shards::Line< 2 >,-1
>::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Line< 2 >,
shards::Line< 2 >,-1 >::createNewElements " void
stk::adapt::RefinerPattern< shards::Line< 2 >, shards::Line< 2 >,-1
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1RefinerPattern_3_01shards_1_1Tetrahedron_3_014_01_4_00_01shards_1_1Tetrahedron_3_014_01_4_00-1_01_4.xml
%feature("docstring") stk::adapt::RefinerPattern< shards::Tetrahedron<
4 >, shards::Tetrahedron< 4 >,-1 > "

general refinement pattern

C++ includes: RefinerPattern_Tet4_Tet4_N.hpp ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::RefinerPattern " stk::adapt::RefinerPattern< shards::Tetrahedron< 4
>, shards::Tetrahedron< 4 >,-1 >::RefinerPattern(percept::PerceptMesh
&eMesh, BlockNamesType block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::~RefinerPattern " stk::adapt::RefinerPattern< shards::Tetrahedron<
4 >, shards::Tetrahedron< 4 >,-1 >::~RefinerPattern() ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::setSubPatterns " void stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1 >::doBreak "
virtual void stk::adapt::RefinerPattern< shards::Tetrahedron< 4 >,
shards::Tetrahedron< 4 >,-1 >::doBreak() ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::fillNeededEntities " void stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::getNumNewElemPerElem " virtual unsigned stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::createNewElements " void stk::adapt::RefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >,-1
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, std::vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1RefinerPattern_3_01shards_1_1Triangle_3_013_01_4_00_01shards_1_1Triangle_3_013_01_4_00_012_01_4.xml
%feature("docstring") stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 > "

this is for testing only - or could be unsed in future for a
bisection-based refinement scheme

C++ includes: RefinerPattern_Tri3_Tri3_2.hpp ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::RefinerPattern "
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::RefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::~RefinerPattern "
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::~RefinerPattern() ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::setSubPatterns " void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::doBreak " virtual void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::doBreak() ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::fillNeededEntities " void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::getNumNewElemPerElem " virtual unsigned
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 2 >::createNewElements " void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>, 2 >::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1RefinerPattern_3_01shards_1_1Triangle_3_013_01_4_00_01shards_1_1Triangle_3_013_01_4_00-1_01_4.xml
%feature("docstring") stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 > "

general refinement pattern

C++ includes: RefinerPattern_Tri3_Tri3_N.hpp ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::RefinerPattern "
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::RefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::~RefinerPattern "
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::~RefinerPattern() ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::setSubPatterns " void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::doBreak " virtual void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::doBreak() ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::fillNeededEntities " void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::getNumNewElemPerElem " virtual unsigned
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::RefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >,-1 >::createNewElements " void
stk::adapt::RefinerPattern< shards::Triangle< 3 >, shards::Triangle< 3
>,-1 >::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1RefinerUtil.xml
%feature("docstring") stk::adapt::RefinerUtil "";


// File: classstk_1_1adapt_1_1Elem_1_1StdMeshObjTopologies_1_1RefTopoX1.xml
%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::RefTopoX1 "";


// File: classstk_1_1percept_1_1RunEnvironment.xml
%feature("docstring") stk::percept::RunEnvironment "";

%feature("docstring")  stk::percept::RunEnvironment::RunEnvironment "stk::percept::RunEnvironment::RunEnvironment(int *argc, char ***argv,
bool debug=false) ";

%feature("docstring")  stk::percept::RunEnvironment::RunEnvironment "stk::percept::RunEnvironment::RunEnvironment(int *argc, char ***argv,
stk::ParallelMachine comm, bool debug=false) ";

%feature("docstring")
stk::percept::RunEnvironment::processCommandLine "int
stk::percept::RunEnvironment::processCommandLine(int argc, char
**argv) ";

%feature("docstring")
stk::percept::RunEnvironment::processCommandLine "int
stk::percept::RunEnvironment::processCommandLine() ";

%feature("docstring")  stk::percept::RunEnvironment::~RunEnvironment "stk::percept::RunEnvironment::~RunEnvironment() ";

%feature("docstring")  stk::percept::RunEnvironment::printHelp "void
stk::percept::RunEnvironment::printHelp() ";

%feature("docstring")
stk::percept::RunEnvironment::build_log_description "std::string
stk::percept::RunEnvironment::build_log_description(const std::string
&working_directory, int parallel_rank, int parallel_size) ";

%feature("docstring")  stk::percept::RunEnvironment::get_argc "int
stk::percept::RunEnvironment::get_argc() ";

%feature("docstring")  stk::percept::RunEnvironment::get_argv "char**
stk::percept::RunEnvironment::get_argv() ";


// File: classstk_1_1percept_1_1Searcher.xml
%feature("docstring") stk::percept::Searcher "";

%feature("docstring")  stk::percept::Searcher::findElement "virtual
const stk::mesh::Entity* stk::percept::Searcher::findElement(MDArray
&input_phy_points, MDArray &found_parametric_coordinates, unsigned
&found_it, const mesh::Entity *hint_element)=0

Find the element containing this physical point and return if found
(also set the found_it flag to 1, else 0). If hint_element is non-
null, use it to check first if it contains the point to potentially
avoid a more costly search.

Dimensions of input_phy_points = ([P]=1, [D]) Dimensions of
found_parametric_coordinates = ([P]=1, [D]) ";

%feature("docstring")  stk::percept::Searcher::setupSearch "virtual
void stk::percept::Searcher::setupSearch() ";

%feature("docstring")  stk::percept::Searcher::tearDownSearch "virtual void stk::percept::Searcher::tearDownSearch() ";

%feature("docstring")  stk::percept::Searcher::~Searcher "virtual
stk::percept::Searcher::~Searcher() ";


// File: classSENS.xml
%feature("docstring") SENS "";


// File: classstk_1_1adapt_1_1SerializeNodeRegistry.xml
%feature("docstring") stk::adapt::SerializeNodeRegistry "";

%feature("docstring")
stk::adapt::SerializeNodeRegistry::SerializeNodeRegistry "stk::adapt::SerializeNodeRegistry::SerializeNodeRegistry(PerceptMesh
&eMesh, NodeRegistry *nodeRegistry, std::string input_mesh_name,
std::string output_mesh_name, int M, int iM, int W=1, int iW=0, int
M_0=-1, int M_1=-1) ";


// File: classstk_1_1adapt_1_1unit__tests_1_1SetRefineField.xml
%feature("docstring") stk::adapt::unit_tests::SetRefineField "";

%feature("docstring")
stk::adapt::unit_tests::SetRefineField::SetRefineField "stk::adapt::unit_tests::SetRefineField::SetRefineField(percept::PerceptMesh
&eMesh) ";

%feature("docstring")
stk::adapt::unit_tests::SetRefineField::init_elementOp "virtual void
stk::adapt::unit_tests::SetRefineField::init_elementOp() ";

%feature("docstring")
stk::adapt::unit_tests::SetRefineField::fini_elementOp "virtual void
stk::adapt::unit_tests::SetRefineField::fini_elementOp() ";


// File: classstk_1_1adapt_1_1regression__tests_1_1SetRefineField.xml
%feature("docstring") stk::adapt::regression_tests::SetRefineField "";

%feature("docstring")
stk::adapt::regression_tests::SetRefineField::SetRefineField "stk::adapt::regression_tests::SetRefineField::SetRefineField(percept::PerceptMesh
&eMesh) ";

%feature("docstring")
stk::adapt::regression_tests::SetRefineField::init_elementOp "virtual
void stk::adapt::regression_tests::SetRefineField::init_elementOp() ";

%feature("docstring")
stk::adapt::regression_tests::SetRefineField::fini_elementOp "virtual
void stk::adapt::regression_tests::SetRefineField::fini_elementOp() ";


// File: classstk_1_1adapt_1_1unit__tests_1_1SetUnrefineField.xml
%feature("docstring") stk::adapt::unit_tests::SetUnrefineField "";

%feature("docstring")
stk::adapt::unit_tests::SetUnrefineField::SetUnrefineField "stk::adapt::unit_tests::SetUnrefineField::SetUnrefineField(percept::PerceptMesh
&eMesh) ";

%feature("docstring")
stk::adapt::unit_tests::SetUnrefineField::init_elementOp "virtual
void stk::adapt::unit_tests::SetUnrefineField::init_elementOp() ";

%feature("docstring")
stk::adapt::unit_tests::SetUnrefineField::fini_elementOp "virtual
void stk::adapt::unit_tests::SetUnrefineField::fini_elementOp() ";


// File: classstk_1_1adapt_1_1regression__tests_1_1SetUnrefineField.xml
%feature("docstring") stk::adapt::regression_tests::SetUnrefineField "";

%feature("docstring")
stk::adapt::regression_tests::SetUnrefineField::SetUnrefineField "stk::adapt::regression_tests::SetUnrefineField::SetUnrefineField(percept::PerceptMesh
&eMesh) ";

%feature("docstring")
stk::adapt::regression_tests::SetUnrefineField::init_elementOp "virtual void
stk::adapt::regression_tests::SetUnrefineField::init_elementOp() ";

%feature("docstring")
stk::adapt::regression_tests::SetUnrefineField::fini_elementOp "virtual void
stk::adapt::regression_tests::SetUnrefineField::fini_elementOp() ";


// File: classShape.xml
%feature("docstring") Shape "";

%feature("docstring")  Shape::Shape "Shape::Shape() ";

%feature("docstring")  Shape::~Shape "virtual Shape::~Shape() ";

%feature("docstring")  Shape::move "void Shape::move(double dx,
double dy) ";

%feature("docstring")  Shape::area "virtual double
Shape::area(void)=0 ";

%feature("docstring")  Shape::perimeter "virtual double
Shape::perimeter(void)=0 ";


// File: classShape__.xml
%feature("docstring") Shape_ "";

%feature("docstring")  Shape_::Shape_ "Shape_< ElementType
>::Shape_(ElementType element_type) ";

%feature("docstring")  Shape_::Shape_ "Shape_< ElementType
>::Shape_(ElementType element_type) ";


// File: classstk_1_1percept_1_1ShardsInterfaceTable.xml
%feature("docstring") stk::percept::ShardsInterfaceTable "";

%feature("docstring")
stk::percept::ShardsInterfaceTable::lookupShardsId "int
stk::percept::ShardsInterfaceTable::lookupShardsId(const char *) ";


// File: structstk_1_1adapt_1_1regression__tests_1_1ShockBasedRefinePredicate.xml
%feature("docstring")
stk::adapt::regression_tests::ShockBasedRefinePredicate "";

%feature("docstring")
stk::adapt::regression_tests::ShockBasedRefinePredicate::ShockBasedRefinePredicate
"stk::adapt::regression_tests::ShockBasedRefinePredicate::ShockBasedRefinePredicate(stk::mesh::FieldBase
*nodal_refine_field, percept::PerceptMesh &eMesh, stk::mesh::Selector
*selector, stk::mesh::FieldBase *field, double tolerance, PlaneShock
shock, double shock_displacement=0, double shock_diff_criterion=0.4)
";


// File: structstk_1_1adapt_1_1SierraPort.xml
%feature("docstring") stk::adapt::SierraPort "";


// File: classstk_1_1percept_1_1SimpleSearcher.xml
%feature("docstring") stk::percept::SimpleSearcher "";

%feature("docstring")  stk::percept::SimpleSearcher::SimpleSearcher "stk::percept::SimpleSearcher::SimpleSearcher(stk::mesh::BulkData
*bulk) ";

%feature("docstring")  stk::percept::SimpleSearcher::~SimpleSearcher "virtual stk::percept::SimpleSearcher::~SimpleSearcher() ";

%feature("docstring")  stk::percept::SimpleSearcher::findElement "const stk::mesh::Entity *
stk::percept::SimpleSearcher::findElement(MDArray &input_phy_points,
MDArray &found_parametric_coordinates, unsigned &found_it, const
mesh::Entity *hint_element)

Dimensions of input_phy_points = ([P]=1, [D]) Dimensions of
found_parametric_coordinates = ([P]=1, [D]) ";


// File: classmoab_1_1SimplexTemplateRefiner.xml
%feature("docstring") moab::SimplexTemplateRefiner "

This class comes from MOAB, modified for STK Percept/Adapt

MOAB 4.1.0RC1 Released June 1, 2011

http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB

Modifications center around exposing only the refine_3_simplex method
as a static method, and changing the algorithm to avoid introducing
face nodes to disambiguate cases with equal edge lengths. Also, memory
management and output functors removed to just triangulate a given tet
with given edge marks with one level of subdivision.

Changes are mostly identify with comments starting with \"p\".

This is a concrete subclass of EntityRefiner that implements
refinement using templates applied to simplices. Entities that are not
simplices are divided into tetrahedra, triangles, or lines before
being processed. Points are passed through unchanged.

David Thompson

Philippe Pebay

C++ includes: Percept_MOAB_SimplexTemplateRefiner.hpp ";

%feature("docstring")
moab::SimplexTemplateRefiner::SimplexTemplateRefiner "moab::SimplexTemplateRefiner::SimplexTemplateRefiner() ";

%feature("docstring")
moab::SimplexTemplateRefiner::~SimplexTemplateRefiner "virtual
moab::SimplexTemplateRefiner::~SimplexTemplateRefiner() ";

%feature("docstring")  moab::SimplexTemplateRefiner::refine_3_simplex
"bool moab::SimplexTemplateRefiner::refine_3_simplex(std::vector<
TetTupleInt > &new_tets, unsigned edge_marks[6], int max_depth, double
*v0, void *t0, EntityHandle h0, double *v1, void *t1, EntityHandle h1,
double *v2, void *t2, EntityHandle h2, double *v3, void *t3,
EntityHandle h3)

Refine a tetrahedron. ";

%feature("docstring")
moab::SimplexTemplateRefiner::heap_coord_storage "double*
moab::SimplexTemplateRefiner::heap_coord_storage() ";

%feature("docstring")  moab::SimplexTemplateRefiner::heap_tag_storage
"void* moab::SimplexTemplateRefiner::heap_tag_storage() ";

%feature("docstring")  moab::SimplexTemplateRefiner::best_tets "int
moab::SimplexTemplateRefiner::best_tets(int *alternates, double *[14],
int, int) ";


// File: classstk_1_1percept_1_1SingleTetFixture.xml
%feature("docstring") stk::percept::SingleTetFixture "

Use case with mixed element topologies and field relations to provide
fast access to node field data from an element.

copied from stk_mesh and modified

C++ includes: SingleTetFixture.hpp ";

%feature("docstring")
stk::percept::SingleTetFixture::~SingleTetFixture "stk::percept::SingleTetFixture::~SingleTetFixture() ";

%feature("docstring")
stk::percept::SingleTetFixture::SingleTetFixture "stk::percept::SingleTetFixture::SingleTetFixture(stk::ParallelMachine
comm, bool doCommit=true, unsigned npts=0, Point *points=0, unsigned
ntets=0, TetIds *tetIds=0, stk::mesh::EntityId elem_id_start=0) ";

%feature("docstring")  stk::percept::SingleTetFixture::populate "void
stk::percept::SingleTetFixture::populate() ";


// File: structstk_1_1adapt_1_1Specialization.xml
%feature("docstring") stk::adapt::Specialization "";


// File: classSquare.xml
%feature("docstring") Square "";

%feature("docstring")  Square::Square "Square::Square(double w) ";

%feature("docstring")  Square::area "virtual double
Square::area(void) ";

%feature("docstring")  Square::perimeter "virtual double
Square::perimeter(void) ";


// File: structStandardReferenceElement_3_01Quadrilateral4_01_4.xml
%feature("docstring") StandardReferenceElement< Quadrilateral4 > " ";


// File: structstk_1_1adapt_1_1STK__Adapt__Auto__Part.xml
%feature("docstring") stk::adapt::STK_Adapt_Auto_Part "

signifies a part that has been defined automatically during adaptivity

C++ includes: UniformRefinerPattern.hpp ";


// File: classstk_1_1percept_1_1STKSearcher.xml
%feature("docstring") stk::percept::STKSearcher "";

%feature("docstring")  stk::percept::STKSearcher::STKSearcher "stk::percept::STKSearcher< SpatialDim
>::STKSearcher(stk::mesh::BulkData *bulk) ";

%feature("docstring")  stk::percept::STKSearcher::~STKSearcher "stk::percept::STKSearcher< SpatialDim >::~STKSearcher() ";

%feature("docstring")  stk::percept::STKSearcher::setupSearch "void
stk::percept::STKSearcher< SpatialDim >::setupSearch() ";

%feature("docstring")  stk::percept::STKSearcher::tearDownSearch "void stk::percept::STKSearcher< SpatialDim >::tearDownSearch() ";

%feature("docstring")  stk::percept::STKSearcher::findElement "const
stk::mesh::Entity * stk::percept::STKSearcher< SpatialDim
>::findElement(MDArray &input_phy_points, MDArray
&found_parametric_coordinates, unsigned &found_it, const mesh::Entity
*hint_element)

Dimensions of input_phy_points = ([P]=1, [D]) Dimensions of
found_parametric_coordinates = ([P]=1, [D]) ";


// File: classstk_1_1percept_1_1StringFunction.xml
%feature("docstring") stk::percept::StringFunction "";

%feature("docstring")  stk::percept::StringFunction::StringFunction "stk::percept::StringFunction::StringFunction(const char
*function_string, Name name=Name(\"noname\"), int domain_dimension=3,
int codomain_dimension=1, unsigned integration_order=0) ";

%feature("docstring")  stk::percept::StringFunction::getFunctionString
"std::string stk::percept::StringFunction::getFunctionString() ";

%feature("docstring")
stk::percept::StringFunction::set_gradient_strings "void
stk::percept::StringFunction::set_gradient_strings(std::string
gstring[3], int len) ";

%feature("docstring")
stk::percept::StringFunction::set_gradient_strings "void
stk::percept::StringFunction::set_gradient_strings(MDArrayString
&gstring) ";

%feature("docstring")  stk::percept::StringFunction::derivative "Teuchos::RCP< Function >
stk::percept::StringFunction::derivative(MDArrayString &deriv_spec)

Return a function that is the derivative of this function. The
derivative is specified as a rank-2 array of strings that specify what
derivative to take and how many derivatives. For example, ";

%feature("docstring")  stk::percept::StringFunction::gradient "Teuchos::RCP< Function > stk::percept::StringFunction::gradient(int
spatialDim=3) ";

%feature("docstring")  stk::percept::StringFunction::StringFunction "stk::percept::StringFunction::StringFunction(const char
*function_string, Name name, Dimensions domain_dimensions, Dimensions
codomain_dimensions, unsigned integration_order=0) ";

%feature("docstring")  stk::percept::StringFunction::StringFunction "stk::percept::StringFunction::StringFunction(const StringFunction &s)
";

%feature("docstring")  stk::percept::StringFunction::resolve "void
stk::percept::StringFunction::resolve(stk::expreval::VariableMap::iterator
&var_it) ";

%feature("docstring")  stk::percept::StringFunction::derivative_test "Teuchos::RCP< Function >
stk::percept::StringFunction::derivative_test(MDArrayString
&deriv_spec) ";

%feature("docstring")
stk::percept::StringFunction::derivative_test_fd "Teuchos::RCP<
Function >
stk::percept::StringFunction::derivative_test_fd(MDArrayString
&deriv_spec, double eps=1.e-6) ";


// File: classstk_1_1adapt_1_1SubDimCell.xml
%feature("docstring") stk::adapt::SubDimCell "

We assume we don't have any sub-dimensional entities with more than 4
nodes.

C++ includes: SubDimCell.hpp ";

%feature("docstring")  stk::adapt::SubDimCell::SubDimCell "stk::adapt::SubDimCell< T, N, CompareClass >::SubDimCell() ";

%feature("docstring")  stk::adapt::SubDimCell::SubDimCell "stk::adapt::SubDimCell< T, N, CompareClass >::SubDimCell(unsigned n)
";

%feature("docstring")  stk::adapt::SubDimCell::insert "void
stk::adapt::SubDimCell< T, N, CompareClass >::insert(T val) ";

%feature("docstring")  stk::adapt::SubDimCell::hashCode "int
stk::adapt::SubDimCell< T, N, CompareClass >::hashCode() ";

%feature("docstring")  stk::adapt::SubDimCell::getHash "unsigned
stk::adapt::SubDimCell< T, N, CompareClass >::getHash() const ";

%feature("docstring")  stk::adapt::SubDimCell::setHash "void
stk::adapt::SubDimCell< T, N, CompareClass >::setHash(std::size_t
hash) ";

%feature("docstring")  stk::adapt::SubDimCell::clear "void
stk::adapt::SubDimCell< T, N, CompareClass >::clear() ";

%feature("docstring")  stk::adapt::SubDimCell::hashCode "int
stk::adapt::SubDimCell< SDSEntityType, 4, SubDimCellCompare<
SDSEntityType > >::hashCode() ";


// File: structstk_1_1adapt_1_1SubDimCell__compare.xml
%feature("docstring") stk::adapt::SubDimCell_compare "";


// File: structstk_1_1adapt_1_1SubDimCellCompare.xml
%feature("docstring") stk::adapt::SubDimCellCompare "";


// File: classstk_1_1utils_1_1SweepMesher.xml
%feature("docstring") stk::utils::SweepMesher "

A simple utility to product tensor product (line, quad, hex) meshes by
sweeping as well as non-tensor product mesh by breaking into sub-
elements (tri, tet, wedge, pyramid).

Steve Kennon, Brian Carnes, Kevin Copps  Usage: initialize with a
simple pair of node, element arrays, such as

double coords[][3] = { {0,0,0}, {1,0,0}, {2,2,0}, {0,3,0}, {0,0,1},
{1,0,1}, {2,2,1}, {0,3,1} };

unsigned quad4Elems[] = { 0,1,2,3, 4,5,6,7 };

SweepMesher tp; tp.initNodes(coords, 8); tp.initElems(elemType, // one
of enum's defined below quad4Elems, 2);

Then use sweep to create a hex mesh (this example breaks a quad to
create two Tri's, then creates a mixed hex/wedge mesh)

boost::array< double, 3> dir = {0,0,1}; std::vector<Transform *>
xforms(1, &TransformDir( dir ) );

// break one of the quads into tris unsigned quadElemIndex = 1;
tp2.breakElem<SweepMesher::ET_Quad4,
SweepMesher::ET_Tri3>(quadElemIndex); std::cout << \"after
break\\\\n\"; tp2.dump();

// sweep to make a hex mesh boost::array< double, 3> dir1 =
{0,0,2.345}; xforms[0] = &TransformDir(dir1); tp2.sweep(
SweepMesher::ET_Quad4, SweepMesher::ET_Hex8, xforms);

C++ includes: SweepMesher.hpp ";


// File: classstk_1_1percept_1_1SweepMesher.xml
%feature("docstring") stk::percept::SweepMesher "";

%feature("docstring")  stk::percept::SweepMesher::SweepMesher "stk::percept::SweepMesher::SweepMesher(unsigned spatialDim=3)

only a few sweep types allowed so far; later could add quadratic
sweeping ";

%feature("docstring")  stk::percept::SweepMesher::~SweepMesher "stk::percept::SweepMesher::~SweepMesher() ";

%feature("docstring")  stk::percept::SweepMesher::initialize "void
stk::percept::SweepMesher::initialize() ";

%feature("docstring")  stk::percept::SweepMesher::CopyFromBasicMesh "void stk::percept::SweepMesher::CopyFromBasicMesh(SweepMesher &source)
";

%feature("docstring")  stk::percept::SweepMesher::get_bulk_data "stk::mesh::BulkData* stk::percept::SweepMesher::get_bulk_data() ";

%feature("docstring")  stk::percept::SweepMesher::getMetaData "stk::mesh::fem::FEMMetaData* stk::percept::SweepMesher::getMetaData()
";

%feature("docstring")  stk::percept::SweepMesher::initNodes "void
stk::percept::SweepMesher::initNodes(double coords[][3], unsigned
numNodes) ";

%feature("docstring")  stk::percept::SweepMesher::initNodes "void
stk::percept::SweepMesher::initNodes(Coord coords[], unsigned
numNodes) ";

%feature("docstring")  stk::percept::SweepMesher::initElems "void
stk::percept::SweepMesher::initElems(unsigned elemType, unsigned
indices[], unsigned numElem) ";

%feature("docstring")  stk::percept::SweepMesher::sweep "void
stk::percept::SweepMesher::sweep(unsigned elemType, unsigned
sweptElemType, std::vector< Transform * > xforms) ";

%feature("docstring")  stk::percept::SweepMesher::sweep "void
stk::percept::SweepMesher::sweep(VectorOfInt elemTypes, VectorOfInt
sweptElemTypes, std::vector< Transform * > xforms)

for a specified group of element types in the mesh at once ";

%feature("docstring")  stk::percept::SweepMesher::sweep "void
stk::percept::SweepMesher::sweep(std::vector< Transform * > xforms)

for all element types in the mesh at once ";

%feature("docstring")  stk::percept::SweepMesher::sweep "void
stk::percept::SweepMesher::sweep(const double path[][3], unsigned
npts)

for all element types in the mesh at once - path following ";

%feature("docstring")  stk::percept::SweepMesher::sweep "void
stk::percept::SweepMesher::sweep(const VectorOfCoord &path) ";

%feature("docstring")  stk::percept::SweepMesher::sweep "void
stk::percept::SweepMesher::sweep(const VectorOfCoord &path, const
VectorOfCoord &dir) ";

%feature("docstring")  stk::percept::SweepMesher::transform "void
stk::percept::SweepMesher::transform(Transform &xform)

apply a single transformation to all nodes' coordinates ";

%feature("docstring")  stk::percept::SweepMesher::squareMesh "void
stk::percept::SweepMesher::squareMesh(unsigned nx, unsigned ny, double
xlength, double ylength, double xorigin=0.0, double yorigin=0.0) ";

%feature("docstring")  stk::percept::SweepMesher::cubeMesh "void
stk::percept::SweepMesher::cubeMesh(unsigned nx, unsigned ny, unsigned
nz, double xlength, double ylength, double zlength, double
xorigin=0.0, double yorigin=0.0, double zorigin=0.0) ";

%feature("docstring")  stk::percept::SweepMesher::debug "void
stk::percept::SweepMesher::debug(const char *str) ";

%feature("docstring")  stk::percept::SweepMesher::debug "void
stk::percept::SweepMesher::debug(const char *str, const int i) ";

%feature("docstring")  stk::percept::SweepMesher::breakElement "void
stk::percept::SweepMesher::breakElement< shards_Hexahedron_8,
shards_Tetrahedron_4 >(unsigned elemIndex)

Break a wedge element into 3 (common case) or 8 (rare case) tets using
a constructive algorithm.

Break a hex element into 6 (common case) or 12 (rare case) tets using
a constructive algorithm. Note: the 5-tet case cannot be constructed
with this algorithm - a table lookup scheme would be more efficient
and allow for the 5-tet case

Rare case where all valences are 3 (each face is broken in same
twisting direction - this is the classic \"un-tetrahedralizable\"
configuration, the Schonhardt prism) - in this case, we have to add a
Steiner point. The point can be added along one face's diagonal
midpoint, but this introduces a need to investigate neighbors, so, we
simply choose to create more tets by using the centroid as the Steiner
point.
(cf.http://www.ams.org/journals/spmj/2005-16-04/S1061-0022-05-00872-1/S1061-0022-05-00872-1.pdf
St. Petersburg Math. J. Tom. 16 (2004), vyp. 4 Vol. 16 (2005), No. 4,
Pages 673690 S 1061-0022(05)00872-1 Article electronically published
on June 24, 2005 REGULAR TRIANGULATIONS AND STEINER POINTS M. YU.
ZVAGELSKI I, A. V. PROSKURNIKOV, AND YU. R. ROMANOVSKI I )

normal case - connect max valence node to other faces without that
node (should be 3 tets always)

Rare case - create tets by joining centroid to each face - for now,
just throw an exception to see how often this case occurs - FIXME -
take this exception out later

normal case - connect max valence node to other faces without that
node (should be 6 tets always) Note: there is a 5-tet configuration
that exists for some face diagonal configurations - FIXME - could add
this case later The 5-tet case consists of an interior tet with no
boundary faces, and 4 corner tets; the boundary faces have to each
have alternating diagonals along the 3 axis directions for this
configuration to exist ";

%feature("docstring")  stk::percept::SweepMesher::breakAllElements "void stk::percept::SweepMesher::breakAllElements() ";

%feature("docstring")  stk::percept::SweepMesher::dumpSTK "void
stk::percept::SweepMesher::dumpSTK() ";

%feature("docstring")  stk::percept::SweepMesher::dump "void
stk::percept::SweepMesher::dump(bool onOff) ";

%feature("docstring")  stk::percept::SweepMesher::dump "void
stk::percept::SweepMesher::dump() ";

%feature("docstring")  stk::percept::SweepMesher::stkMeshCreate "void
stk::percept::SweepMesher::stkMeshCreate(stk::ParallelMachine &)

create a std::mesh representation of this

based on UseCase_3 in stk_mesh/use_cases - creates nodes and elements
in stk::mesh database ";

%feature("docstring")
stk::percept::SweepMesher::stkMeshCreateMetaNoCommit "void
stk::percept::SweepMesher::stkMeshCreateMetaNoCommit(stk::ParallelMachine
&) ";

%feature("docstring")
stk::percept::SweepMesher::stkMeshCreateBulkAfterMetaCommit "void
stk::percept::SweepMesher::stkMeshCreateBulkAfterMetaCommit(stk::ParallelMachine
&) ";

%feature("docstring")  stk::percept::SweepMesher::writeSTKMesh "void
stk::percept::SweepMesher::writeSTKMesh(const char *filename) ";


// File: classTCoeff.xml
%feature("docstring") TCoeff "";

%feature("docstring")  TCoeff::TCoeff "TCoeff::TCoeff() ";

%feature("docstring")  TCoeff::TCoeff "TCoeff::TCoeff() ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1temp.xml
%feature("docstring") stk::percept::IntrepidManager::temp "";


// File: classTemperature__.xml
%feature("docstring") Temperature_ "";

%feature("docstring")  Temperature_::Temperature_ "Temperature_<
TCoeff, Shape >::Temperature_(TCoeff t, Shape shape) ";

%feature("docstring")  Temperature_::Temperature_ "Temperature_<
TCoeff, Shape >::Temperature_(TCoeff t, Shape shape) ";


// File: classstk_1_1adapt_1_1TestLocalRefiner.xml
%feature("docstring") stk::adapt::TestLocalRefiner "

A test implementation that does uniform refinement but uses non-
uniform methods

C++ includes: TestLocalRefiner.hpp ";

%feature("docstring")  stk::adapt::TestLocalRefiner::TestLocalRefiner
"stk::adapt::TestLocalRefiner::TestLocalRefiner(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTet__N__1.xml
%feature("docstring") stk::adapt::TestLocalRefinerTet_N_1 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_3

C++ includes: TestLocalRefinerTet_N_1.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_1::TestLocalRefinerTet_N_1 "stk::adapt::TestLocalRefinerTet_N_1::TestLocalRefinerTet_N_1(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTet__N__2.xml
%feature("docstring") stk::adapt::TestLocalRefinerTet_N_2 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_3

C++ includes: TestLocalRefinerTet_N_2.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_2::TestLocalRefinerTet_N_2 "stk::adapt::TestLocalRefinerTet_N_2::TestLocalRefinerTet_N_2(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0, unsigned mark_first_n_edges=1) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTet__N__2__1.xml
%feature("docstring") stk::adapt::TestLocalRefinerTet_N_2_1 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_3

C++ includes: TestLocalRefinerTet_N_2_1.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_2_1::TestLocalRefinerTet_N_2_1 "stk::adapt::TestLocalRefinerTet_N_2_1::TestLocalRefinerTet_N_2_1(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0, unsigned edge_mark_bitcode=1) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTet__N__3.xml
%feature("docstring") stk::adapt::TestLocalRefinerTet_N_3 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_3

C++ includes: TestLocalRefinerTet_N_3.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_3::TestLocalRefinerTet_N_3 "stk::adapt::TestLocalRefinerTet_N_3::TestLocalRefinerTet_N_3(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTet__N__3__1.xml
%feature("docstring") stk::adapt::TestLocalRefinerTet_N_3_1 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_3_1

C++ includes: TestLocalRefinerTet_N_3_1.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_3_1::TestLocalRefinerTet_N_3_1 "stk::adapt::TestLocalRefinerTet_N_3_1::TestLocalRefinerTet_N_3_1(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0, unsigned edge_mark_bitcode=1) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTet__N__4.xml
%feature("docstring") stk::adapt::TestLocalRefinerTet_N_4 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_4

C++ includes: TestLocalRefinerTet_N_4.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_4::TestLocalRefinerTet_N_4 "stk::adapt::TestLocalRefinerTet_N_4::TestLocalRefinerTet_N_4(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0, unsigned edge_mark_bitcode=1) ";

%feature("docstring")
stk::adapt::TestLocalRefinerTet_N_4::buildTestUnrefineList "ElementUnrefineCollection
stk::adapt::TestLocalRefinerTet_N_4::buildTestUnrefineList() ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri "

A test implementation that does uniform refinement but uses non-
uniform methods

C++ includes: TestLocalRefinerTri.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri::TestLocalRefinerTri "stk::adapt::TestLocalRefinerTri::TestLocalRefinerTri(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri1.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri1 "

A test implementation that does uniform refinement but uses non-
uniform methods

C++ includes: TestLocalRefinerTri1.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri1::TestLocalRefinerTri1 "stk::adapt::TestLocalRefinerTri1::TestLocalRefinerTri1(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri2.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri2 "

A test implementation that does uniform refinement but uses non-
uniform methods

C++ includes: TestLocalRefinerTri2.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri2::TestLocalRefinerTri2 "stk::adapt::TestLocalRefinerTri2::TestLocalRefinerTri2(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0, bool diagonals=true) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri_N "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N

C++ includes: TestLocalRefinerTri_N.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N::TestLocalRefinerTri_N "stk::adapt::TestLocalRefinerTri_N::TestLocalRefinerTri_N(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N::buildTestUnrefineList "ElementUnrefineCollection
stk::adapt::TestLocalRefinerTri_N::buildTestUnrefineList() ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__1.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri_N_1 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_1

C++ includes: TestLocalRefinerTri_N_1.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_1::TestLocalRefinerTri_N_1 "stk::adapt::TestLocalRefinerTri_N_1::TestLocalRefinerTri_N_1(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_1::buildTestUnrefineList "ElementUnrefineCollection
stk::adapt::TestLocalRefinerTri_N_1::buildTestUnrefineList() ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__2.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri_N_2 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_2

C++ includes: TestLocalRefinerTri_N_2.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_2::TestLocalRefinerTri_N_2 "stk::adapt::TestLocalRefinerTri_N_2::TestLocalRefinerTri_N_2(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_2::buildTestUnrefineList "ElementUnrefineCollection
stk::adapt::TestLocalRefinerTri_N_2::buildTestUnrefineList() ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__3.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri_N_3 "

A test implementation that marks some edges randomly to test
RefinerPattern_Tri3_Tri3_N_3

C++ includes: TestLocalRefinerTri_N_3.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3::TestLocalRefinerTri_N_3 "stk::adapt::TestLocalRefinerTri_N_3::TestLocalRefinerTri_N_3(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3::buildTestUnrefineList "ElementUnrefineCollection
stk::adapt::TestLocalRefinerTri_N_3::buildTestUnrefineList() ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__3__EdgeBasedAnisotropic.xml
%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_EdgeBasedAnisotropic "

A test implementation as a use case for EdgeBasedAnisotropic

C++ includes: TestLocalRefinerTri_N_3_EdgeBasedAnisotropic.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_EdgeBasedAnisotropic::TestLocalRefinerTri_N_3_EdgeBasedAnisotropic
"stk::adapt::TestLocalRefinerTri_N_3_EdgeBasedAnisotropic::TestLocalRefinerTri_N_3_EdgeBasedAnisotropic(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, VectorFieldType
*nodal_hessian_field, stk::mesh::FieldBase *proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__3__IEdgeAdapter.xml
%feature("docstring") stk::adapt::TestLocalRefinerTri_N_3_IEdgeAdapter
"

A test implementation as a use case for IEdgeAdapter

C++ includes: TestLocalRefinerTri_N_3_IEdgeAdapter.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_IEdgeAdapter::TestLocalRefinerTri_N_3_IEdgeAdapter
"stk::adapt::TestLocalRefinerTri_N_3_IEdgeAdapter::TestLocalRefinerTri_N_3_IEdgeAdapter(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__3__IElementAdapter.xml
%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_IElementAdapter "

A test implementation as a use case for IElementAdapter

C++ includes: TestLocalRefinerTri_N_3_IElementAdapter.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_IElementAdapter::TestLocalRefinerTri_N_3_IElementAdapter
"stk::adapt::TestLocalRefinerTri_N_3_IElementAdapter::TestLocalRefinerTri_N_3_IElementAdapter(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1TestLocalRefinerTri__N__3__MeshSizeRatio.xml
%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_MeshSizeRatio "

A test implementation as a use case for MeshSizeRatio

C++ includes: TestLocalRefinerTri_N_3_MeshSizeRatio.hpp ";

%feature("docstring")
stk::adapt::TestLocalRefinerTri_N_3_MeshSizeRatio::TestLocalRefinerTri_N_3_MeshSizeRatio
"stk::adapt::TestLocalRefinerTri_N_3_MeshSizeRatio::TestLocalRefinerTri_N_3_MeshSizeRatio(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, ScalarFieldType
*elem_ratio_field, stk::mesh::FieldBase *proc_rank_field=0) ";


// File: classstk_1_1percept_1_1TopologyVerifier.xml
%feature("docstring") stk::percept::TopologyVerifier "";

%feature("docstring")
stk::percept::TopologyVerifier::TopologyVerifier "stk::percept::TopologyVerifier::TopologyVerifier() ";

%feature("docstring")  stk::percept::TopologyVerifier::isTopologyBad "bool stk::percept::TopologyVerifier::isTopologyBad(mesh::Entity &elem)

return true if topology is bad ";

%feature("docstring")  stk::percept::TopologyVerifier::isTopologyBad "bool stk::percept::TopologyVerifier::isTopologyBad(mesh::BulkData
&mesh_bulk_data) ";


// File: classstk_1_1percept_1_1Transform.xml
%feature("docstring") stk::percept::Transform "";

%feature("docstring")  stk::percept::Transform::Transform "stk::percept::Transform::Transform() ";


// File: classstk_1_1percept_1_1TransformDir.xml
%feature("docstring") stk::percept::TransformDir "";

%feature("docstring")  stk::percept::TransformDir::TransformDir "stk::percept::TransformDir::TransformDir(Coord dir) ";


// File: classstk_1_1percept_1_1util_1_1TransformPath.xml
%feature("docstring") stk::percept::util::TransformPath "";

%feature("docstring")
stk::percept::util::TransformPath::TransformPath "stk::percept::util::TransformPath::TransformPath(const Coord &from,
const Coord &from_dir, const Coord &to, const Coord &to_dir)

Given some points on a plane (\"from\"), we want to move them to a new
plane by rotation and translation The initial plane is defined by an
origin ([in] from), and it's normal ([in] from_dir). The final plane
is defined by origin/normal: [in] to, to_dir The translation delta is
just the vector {to - from}. The algorithm rotates the points in the
initial plane to the new plane's direction, then does the translate.
Note that the points on the initial plane don't have to lie in a
plane, just easier to visualize that way. ";


// File: classstk_1_1adapt_1_1UniformRefiner.xml
%feature("docstring") stk::adapt::UniformRefiner "";

%feature("docstring")  stk::adapt::UniformRefiner::UniformRefiner "stk::adapt::UniformRefiner::UniformRefiner(percept::PerceptMesh
&eMesh, UniformRefinerPatternBase &bp, stk::mesh::FieldBase
*proc_rank_field=0) ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern.xml
%feature("docstring") stk::adapt::UniformRefinerPattern "";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Beam_3_012_01_4_00_01shards_1_1Beam_3_03519a7ae913f8910e3a3ddb9fe2df6b0.xml
%feature("docstring") stk::adapt::UniformRefinerPattern< shards::Beam<
2 >, shards::Beam< 2 >, 2, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Beam< 2 >,
shards::Beam< 2 >, 2, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Beam< 2 >, shards::Beam< 2
>, 2, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 2 >, 2, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Beam_3_012_01_4_00_01shards_1_1Beam_3_080bac1ee2fd59d19e985daebf1216a29.xml
%feature("docstring") stk::adapt::UniformRefinerPattern< shards::Beam<
2 >, shards::Beam< 3 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Beam< 2 >,
shards::Beam< 3 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort >::setSubPatterns
" void stk::adapt::UniformRefinerPattern< shards::Beam< 2 >,
shards::Beam< 3 >, 1, SierraPort >::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Beam< 2 >, shards::Beam< 3
>, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Beam< 2 >, shards::Beam< 3 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Beam_3_013_01_4_00_01shards_1_1Beam_3_059ac980035fc1205cd0415423628387f.xml
%feature("docstring") stk::adapt::UniformRefinerPattern< shards::Beam<
3 >, shards::Beam< 3 >, 2, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Beam< 3 >,
shards::Beam< 3 >, 2, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Beam< 3 >, shards::Beam< 3
>, 2, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Beam< 3 >, shards::Beam< 3 >, 2, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_0120_01_4_00_01shards_1_1H8021c51ddccde0744fb078dc28f4c799.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort > "
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 20 >,
shards::Hexahedron< 20 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 20 >, shards::Hexahedron< 20 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_0127_01_4_00_01shards_1_1H5bf9dad70361e8b239dcfc2ebea00561.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort > "
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 27 >,
shards::Hexahedron< 27 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 27 >, shards::Hexahedron< 27 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_018_01_4_00_01shards_1_1Hedb807cf313e3c9d19f93af2873214965.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort > "
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8 >,
shards::Hexahedron< 20 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 20 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_018_01_4_00_01shards_1_1Hed7669878dcd0afef94f0abc98228f9ba.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort > "
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8 >,
shards::Hexahedron< 27 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 27 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_018_01_4_00_01shards_1_1Heb1712b5fd5edebd325d5502c66a325d7.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8 >,
shards::Hexahedron< 8 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Hexahedron< 8 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_018_01_4_00_01shards_1_1Tee8bef9c857dea3743b6c73531852876f.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24 > "

Topological traits: Dimension = 3, Sides = 6, Edges = 12, Vertices =
8, and Nodes = 8, 20, or 27.

From Shards_BasicTopologies.hpp

Linear 8-Node Hexahedron node locations.

7                    6 o------------------o /|                 /| / |
/ | /  |               /  | /   |              /   | /    |
/    | /     |            /     | 4 /      |         5 /      | o
------------------o       | |       |          |       | |     3
o----------|-------o 2 |      /           |      / |     /
|     / |    /             |    / |   /              |   / |  /
|  / | /                | / |/                 |/ o------------------o
0                    1

face numbering for symmetric hex to tet break pattern |   typedef |
MakeTypeList< IndexList< 0 , 1 ,   8 > , 7
| IndexList< 1 , 2 ,   9 > , o------------------o 6
| IndexList< 2 , 3 ,  10 > , /|                 /|
| IndexList< 3 , 0 ,  11 > ,
/ |                / |                                              |
IndexList< 4 , 5 ,  16 > ,                                         / |
13          /  |                                              |
IndexList< 5 , 6 ,  17 > ,                                        / |
o         /   |                                              |
IndexList< 6 , 7 ,  18 > ,                                       / |
o10   /    |     Node #14 is at centroid of element       | IndexList<
7 , 4 ,  19 > ,                                      / |            /
|                                              | IndexList< 0 , 4 ,
12 > ,                                   4 / |         5 /      |
\"2D surface\" containing nodes            | IndexList< 1 , 5 ,  13 >
,                                    o ------------------o    9  |
0,1,5,4 has node 25 at center.... |                   IndexList< 2 , 6
,  14 > , | 11o   | 3        |   o   | |                   IndexList<
3 , 7 ,  15 > >type |       o----------|-------o 2 |
HexahedronEdgeNodeMap ; |      /           |      / | |     /   8
|     / |   typedef |    /    o        |    / |     MakeTypeList<
IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > , |   /        o12   |
/ |                   IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > ,
|  /               |  / |                   IndexList< 2, 3, 7, 6,
10, 15, 18, 14,   26 > , | /                | / |
IndexList< 0, 4, 7, 3,  12, 19, 15, 11,   23 > , |/                 |/
|                   IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > , o
------------------o |                   IndexList< 4, 5, 6, 7,  16,
17, 18, 19,   22 > >type       0                    1 |
HexahedronFaceNodeMap ; | |

C++ includes: UniformRefinerPattern_Hex8_Tet4_24.hpp ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24 >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8
>, shards::Tetrahedron< 4 >, 24 >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8 >,
shards::Tetrahedron< 4 >, 24 >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 24
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Hexahedron_3_018_01_4_00_01shards_1_1Teea6db3e46ff803882dc213da0e816d67.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6 > "

Topological traits: Dimension = 3, Sides = 6, Edges = 12, Vertices =
8, and Nodes = 8, 20, or 27.

From Shards_BasicTopologies.hpp

Linear 8-Node Hexahedron node locations.

7                    6 o------------------o /|                 /| / |
/ | /  |               /  | /   |              /   | /    |
/    | /     |            /     | 4 /      |         5 /      | o
------------------o       | |       |          |       | |     3
o----------|-------o 2 |      /           |      / |     /
|     / |    /             |    / |   /              |   / |  /
|  / | /                | / |/                 |/ o------------------o
0                    1

face numbering for symmetric hex to tet break pattern |   typedef |
MakeTypeList< IndexList< 0 , 1 ,   8 > , 7
| IndexList< 1 , 2 ,   9 > , o------------------o 6
| IndexList< 2 , 3 ,  10 > , /|                 /|
| IndexList< 3 , 0 ,  11 > ,
/ |                / |                                              |
IndexList< 4 , 5 ,  16 > ,                                         / |
13          /  |                                              |
IndexList< 5 , 6 ,  17 > ,                                        / |
o         /   |                                              |
IndexList< 6 , 7 ,  18 > ,                                       / |
o10   /    |     Node #14 is at centroid of element       | IndexList<
7 , 4 ,  19 > ,                                      / |            /
|                                              | IndexList< 0 , 4 ,
12 > ,                                   4 / |         5 /      |
\"2D surface\" containing nodes            | IndexList< 1 , 5 ,  13 >
,                                    o ------------------o    9  |
0,1,5,4 has node 25 at center.... |                   IndexList< 2 , 6
,  14 > , | 11o   | 3        |   o   | |                   IndexList<
3 , 7 ,  15 > >type |       o----------|-------o 2 |
HexahedronEdgeNodeMap ; |      /           |      / | |     /   8
|     / |   typedef |    /    o        |    / |     MakeTypeList<
IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > , |   /        o12   |
/ |                   IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > ,
|  /               |  / |                   IndexList< 2, 3, 7, 6,
10, 15, 18, 14,   26 > , | /                | / |
IndexList< 0, 4, 7, 3,  12, 19, 15, 11,   23 > , |/                 |/
|                   IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > , o
------------------o |                   IndexList< 4, 5, 6, 7,  16,
17, 18, 19,   22 > >type       0                    1 |
HexahedronFaceNodeMap ; | |

C++ includes: UniformRefinerPattern_Hex8_Tet4_6_12.hpp ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6 >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8
>, shards::Tetrahedron< 4 >, 6 >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6 >::setSubPatterns
" void stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8 >,
shards::Tetrahedron< 4 >, 6 >::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Hexahedron< 8 >,
shards::Tetrahedron< 4 >, 6 >::getNumNewElemPerElem()

NOTE: we create additional un-used elements if the Hex8 can be broken
into 6 tets. ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Hexahedron< 8 >, shards::Tetrahedron< 4 >, 6
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element

Rare case - create tets by joining centroid to each face - for now,
just throw an exception to see how often this case occurs - FIXME -
take this exception out later

normal case - connect max valence node to other faces without that
node (should be 6 tets always) Note: there is a 5-tet configuration
that exists for some face diagonal configurations - FIXME - could add
this case later The 5-tet case consists of an interior tet with no
boundary faces, and 4 corner tets; the boundary faces have to each
have alternating diagonals along the 3 axis directions for this
configuration to exist ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Line_3_012_01_4_00_01shards_1_1Line_3_0eb73084bdefe28105346a456a31628e2.xml
%feature("docstring") stk::adapt::UniformRefinerPattern< shards::Line<
2 >, shards::Line< 2 >, 2, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Line< 2 >,
shards::Line< 2 >, 2, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Line< 2 >, shards::Line< 2
>, 2, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 2 >, 2, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Line_3_012_01_4_00_01shards_1_1Line_3_066a5ad0685bedb17c2b98fdea46088dc.xml
%feature("docstring") stk::adapt::UniformRefinerPattern< shards::Line<
2 >, shards::Line< 3 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Line< 2 >,
shards::Line< 3 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort >::setSubPatterns
" void stk::adapt::UniformRefinerPattern< shards::Line< 2 >,
shards::Line< 3 >, 1, SierraPort >::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Line< 2 >, shards::Line< 3
>, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Line< 2 >, shards::Line< 3 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Line_3_013_01_4_00_01shards_1_1Line_3_0f742194bc9d2883e16f39e10b5178ec6.xml
%feature("docstring") stk::adapt::UniformRefinerPattern< shards::Line<
3 >, shards::Line< 3 >, 2, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Line< 3 >,
shards::Line< 3 >, 2, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Line< 3 >, shards::Line< 3
>, 2, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Line< 3 >, shards::Line< 3 >, 2, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Pyramid_3_0113_01_4_00_01shards_1_1Pyra4fed58b5628ca1720e1146837cd090c8.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Pyramid< 13 >,
shards::Pyramid< 13 >, 10, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 13 >, shards::Pyramid< 13 >, 10, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Pyramid_3_015_01_4_00_01shards_1_1Pyram49fbce98a5575bb93770f4a86852e8e3.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort >::doBreak
" virtual void stk::adapt::UniformRefinerPattern< shards::Pyramid< 5
>, shards::Pyramid< 13 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 13 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 13 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Pyramid_3_015_01_4_00_01shards_1_1Pyram2ea8d2598167b47d065c996a10a50c15.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort >::doBreak
" virtual void stk::adapt::UniformRefinerPattern< shards::Pyramid< 5
>, shards::Pyramid< 5 >, 10, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::getFromTypeKey " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 10, SierraPort >::getFromTypeKey() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::getFromTopoPartName " virtual std::string
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 10, SierraPort >::getFromTopoPartName()

provided by this class
------------------------------------------------------------------------------------------------------------------------
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::getToTopoPartName " virtual std::string
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 10, SierraPort >::getToTopoPartName() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::getFromTopology " virtual const CellTopologyData*
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 10, SierraPort >::getFromTopology() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::getToTopology " virtual const CellTopologyData*
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 10, SierraPort >::getToTopology() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 10, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 10, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Pyramid_3_015_01_4_00_01shards_1_1Pyramab19c3fb962f501ddef2b4aac97182a5.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 6, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Pyramid< 5 >, 6, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Pyramid< 5 >, 6, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Pyramid_3_015_01_4_00_01shards_1_1Tetraa89de1ab4e8283aa79cb4c808ccac13c.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Pyramid< 5 >,
shards::Tetrahedron< 4 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Pyramid< 5 >, shards::Tetrahedron< 4 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_f72d17577cedec1a029247d4f5161333.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4 > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4 >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Quadrilateral<
4 >, shards::Quadrilateral< 4 >, 4 >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Quadrilateral< 4 >, 4 >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_4655f30c0bc4dbca0cec571b46e5ab4c.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Quadrilateral< 4 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 4 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_f3c6b28b6c1a6a12c54dc917efa27cb3.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Quadrilateral< 8 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 8 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_f3bdf5cdf5d0507af1303134649c06c7.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Quadrilateral< 9 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::fixSurfaceAndEdgeSetNamesMap " virtual StringStringMap
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Quadrilateral< 9 >, 1, SierraPort
>::fixSurfaceAndEdgeSetNamesMap()

for i/o to work properly, supply string replacements such as for
hex-->tet breaking, you would supply \"quad\"-->\"tri\" etc. string
maps ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Quadrilateral< 9 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_005c58985ac9d612de4415eb15a0bcca.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2 > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2 >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Quadrilateral<
4 >, shards::Triangle< 3 >, 2 >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Triangle< 3 >, 2 >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::fixSurfaceAndEdgeSetNamesMap " virtual StringStringMap
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Triangle< 3 >, 2 >::fixSurfaceAndEdgeSetNamesMap()

for i/o to work properly, supply string replacements such as for
hex-->tet breaking, you would supply \"quad\"-->\"tri\" etc. string
maps ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 2
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_461207b6251d0064c8d5b56e9a0271b9.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization >
" ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Triangle< 3 >, 4, Specialization >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::fixSurfaceAndEdgeSetNamesMap " virtual StringStringMap
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Triangle< 3 >, 4, Specialization
>::fixSurfaceAndEdgeSetNamesMap()

for i/o to work properly, supply string replacements such as for
hex-->tet breaking, you would supply \"quad\"-->\"tri\" etc. string
maps ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 4, Specialization
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_014_01_4_00_01shards_1_779cc24671f5495df8217c1c52c6ee24.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6 > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6 >::setSubPatterns
" void stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Triangle< 3 >, 6 >::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6 >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Quadrilateral<
4 >, shards::Triangle< 3 >, 6 >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 4 >,
shards::Triangle< 3 >, 6 >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 4 >, shards::Triangle< 3 >, 6
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element

[above] at (p4.side 1){2}; [left] at (p4.side 2){3}; [below] at
(p4.side 3){0}; [right] at (p4.side 4){1}; ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_018_01_4_00_01shards_1_967bcbe24117f1aba5b7bd72e159f462.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 8 >,
shards::Quadrilateral< 8 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 8 >, shards::Quadrilateral< 8 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Quadrilateral_3_019_01_4_00_01shards_1_4c7cf4e80e025a0bbeb1d5e8c98220b0.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Quadrilateral< 9 >,
shards::Quadrilateral< 9 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Quadrilateral< 9 >, shards::Quadrilateral< 9 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellLine_3_012_01_4_00_01shards_1_1She25e5783b043b37be41a34bb285317d14.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellLine< 2 >,
shards::ShellLine< 2 >, 2, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 2 >, 2, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellLine_3_012_01_4_00_01shards_1_1She32f3111d4cc0c58b08dc7c6c97185ea5.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellLine< 2 >,
shards::ShellLine< 3 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 2 >, shards::ShellLine< 3 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellLine_3_013_01_4_00_01shards_1_1She282897e9bced26a4c5efb1e03b9049f0.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellLine< 3 >,
shards::ShellLine< 3 >, 2, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::ShellLine< 3 >, shards::ShellLine< 3 >, 2, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellQuadrilateral_3_014_01_4_00_01shard917e1321f6798d8f4ae5042ee82581f.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 4 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::~UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 4 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::setSubPatterns(std::vector< UniformRefinerPatternBase *
> &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::doBreak " virtual void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 4 >, 4, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::fillNeededEntities " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 4 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 4 >, 4, SierraPort
>::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 4 >, 4,
SierraPort >::createNewElements " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 4 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellQuadrilateral_3_014_01_4_00_01sharc8402000d5fe6fe4721f31f2f1d87158.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 8 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::~UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 8 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::setSubPatterns(std::vector< UniformRefinerPatternBase *
> &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::doBreak " virtual void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 8 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::fillNeededEntities " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 8 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 8 >, 1, SierraPort
>::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 8 >, 1,
SierraPort >::createNewElements " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 8 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellQuadrilateral_3_014_01_4_00_01sharcfda37310a64d3b19332d60f0d1b4d7d.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 9 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::~UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 9 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::setSubPatterns(std::vector< UniformRefinerPatternBase *
> &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::doBreak " virtual void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 9 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::fillNeededEntities " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 9 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 9 >, 1, SierraPort
>::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 4 >, shards::ShellQuadrilateral< 9 >, 1,
SierraPort >::createNewElements " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 4 >,
shards::ShellQuadrilateral< 9 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellQuadrilateral_3_018_01_4_00_01shar8eadcc16bba8d59a778995c552583530.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 8 >,
shards::ShellQuadrilateral< 8 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::~UniformRefinerPattern "
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 8 >,
shards::ShellQuadrilateral< 8 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::setSubPatterns(std::vector< UniformRefinerPatternBase *
> &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::doBreak " virtual void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 8 >,
shards::ShellQuadrilateral< 8 >, 4, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::fillNeededEntities " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 8 >,
shards::ShellQuadrilateral< 8 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 8 >,
shards::ShellQuadrilateral< 8 >, 4, SierraPort
>::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellQuadrilateral< 8 >, shards::ShellQuadrilateral< 8 >, 4,
SierraPort >::createNewElements " void
stk::adapt::UniformRefinerPattern< shards::ShellQuadrilateral< 8 >,
shards::ShellQuadrilateral< 8 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellTriangle_3_013_01_4_00_01shards_1_46e930dce8b0b36d120c7f322680e49a.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellTriangle< 3 >,
shards::ShellTriangle< 3 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 3 >, shards::ShellTriangle< 3 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1ShellTriangle_3_016_01_4_00_01shards_1_4ade42dcdd3b147165b6c76cb82723eb.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
> " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::ShellTriangle< 6 >,
shards::ShellTriangle< 6 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::ShellTriangle< 6 >, shards::ShellTriangle< 6 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Tetrahedron_3_0110_01_4_00_01shards_1_1c7c7e63ced74b2b9d2e4b874c5a303d9.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort >
" ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Tetrahedron< 10 >,
shards::Tetrahedron< 10 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 10 >, shards::Tetrahedron< 10 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Tetrahedron_3_014_01_4_00_01shards_1_1T26ea4ee662d5a303dda37a954db365c5.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort > "
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Tetrahedron< 4 >,
shards::Tetrahedron< 10 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 10 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Tetrahedron_3_014_01_4_00_01shards_1_1T6569d2f59d103085265c7c0bd5d57f66.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort > "
";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::doBreak " virtual void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Tetrahedron< 4 >,
shards::Tetrahedron< 4 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Tetrahedron< 4 >, shards::Tetrahedron< 4 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Triangle_3_013_01_4_00_01shards_1_1Tria39630e1acdd9a0b2b96c6d103a9cce86.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort >::doBreak
" virtual void stk::adapt::UniformRefinerPattern< shards::Triangle< 3
>, shards::Triangle< 3 >, 4, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Triangle< 3 >,
shards::Triangle< 3 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 3 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Triangle_3_013_01_4_00_01shards_1_1Tria34e0ed26dd6f237c4e42df4b98e5001c.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort >::doBreak
" virtual void stk::adapt::UniformRefinerPattern< shards::Triangle< 3
>, shards::Triangle< 6 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Triangle< 3 >,
shards::Triangle< 6 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 3 >, shards::Triangle< 6 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Triangle_3_016_01_4_00_01shards_1_1Tria44b116d2d6f638e5c54117a6bd1ae130.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort >::doBreak
" virtual void stk::adapt::UniformRefinerPattern< shards::Triangle< 6
>, shards::Triangle< 6 >, 4, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Triangle< 6 >,
shards::Triangle< 6 >, 4, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Triangle< 6 >, shards::Triangle< 6 >, 4, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Wedge_3_0115_01_4_00_01shards_1_1Wedge_81b3700b6a0fb0153b52d21ea1743b77.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Wedge< 15 >,
shards::Wedge< 15 >, 8, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Wedge< 15 >, shards::Wedge<
15 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 15 >, shards::Wedge< 15 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Wedge_3_0118_01_4_00_01shards_1_1Wedge_e5a958df3e802a57ee72996381872595.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Wedge< 18 >,
shards::Wedge< 18 >, 8, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Wedge< 18 >, shards::Wedge<
18 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 18 >, shards::Wedge< 18 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Wedge_3_016_01_4_00_01shards_1_1Wedge_34049bec8435403a14fd496e28fdfe807.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Wedge< 6 >,
shards::Wedge< 15 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Wedge< 6 >, shards::Wedge<
15 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 15 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Wedge_3_016_01_4_00_01shards_1_1Wedge_36c02b343d4fdb4a75fb3d6575e72d162.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Wedge< 6 >,
shards::Wedge< 18 >, 1, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Wedge< 6 >, shards::Wedge<
18 >, 1, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 18 >, 1, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPattern_3_01shards_1_1Wedge_3_016_01_4_00_01shards_1_1Wedge_309e7909fabfee73b98e6847af508960c.xml
%feature("docstring") stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort > " ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::UniformRefinerPattern(percept::PerceptMesh &eMesh, BlockNamesType
block_names=BlockNamesType()) ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::~UniformRefinerPattern " stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::~UniformRefinerPattern() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::setSubPatterns " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::setSubPatterns(std::vector< UniformRefinerPatternBase * > &bp,
percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort >::doBreak "
virtual void stk::adapt::UniformRefinerPattern< shards::Wedge< 6 >,
shards::Wedge< 6 >, 8, SierraPort >::doBreak() ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::fillNeededEntities " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::fillNeededEntities(std::vector< NeededEntityType >
&needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::getNumNewElemPerElem " virtual unsigned
stk::adapt::UniformRefinerPattern< shards::Wedge< 6 >, shards::Wedge<
6 >, 8, SierraPort >::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")  stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::createNewElements " void stk::adapt::UniformRefinerPattern<
shards::Wedge< 6 >, shards::Wedge< 6 >, 8, SierraPort
>::createNewElements(percept::PerceptMesh &eMesh, NodeRegistry
&nodeRegistry, stk::mesh::Entity &element, NewSubEntityNodesType
&new_sub_entity_nodes, vector< stk::mesh::Entity * >::iterator
&element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1UniformRefinerPatternBase.xml
%feature("docstring") stk::adapt::UniformRefinerPatternBase "

The base class for all refinement patterns
------------------------------------------------------------------------------------------------------------------------

C++ includes: UniformRefinerPattern.hpp ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::UniformRefinerPatternBase "stk::adapt::UniformRefinerPatternBase::UniformRefinerPatternBase() ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::~UniformRefinerPatternBase "virtual
stk::adapt::UniformRefinerPatternBase::~UniformRefinerPatternBase() ";

%feature("docstring")  stk::adapt::UniformRefinerPatternBase::doBreak
"virtual void stk::adapt::UniformRefinerPatternBase::doBreak()=0 ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getFromTypeKey "virtual
unsigned stk::adapt::UniformRefinerPatternBase::getFromTypeKey()=0 ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getFromTopology "virtual const
CellTopologyData*
stk::adapt::UniformRefinerPatternBase::getFromTopology()=0 ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getToTopology "virtual const
CellTopologyData*
stk::adapt::UniformRefinerPatternBase::getToTopology()=0 ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getPrimaryEntityRank "stk::mesh::EntityRank
stk::adapt::UniformRefinerPatternBase::getPrimaryEntityRank() ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::fillNeededEntities "virtual
void
stk::adapt::UniformRefinerPatternBase::fillNeededEntities(std::vector<
NeededEntityType > &needed_entities)=0

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::setNeededParts "virtual void
stk::adapt::UniformRefinerPatternBase::setNeededParts(percept::PerceptMesh
&eMesh, BlockNamesType block_names_ranks, bool sameTopology=true) ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getNumNewElemPerElem "virtual
unsigned
stk::adapt::UniformRefinerPatternBase::getNumNewElemPerElem()=0

supply the number of new elements per element during refinement ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::createNewElements "virtual
void
stk::adapt::UniformRefinerPatternBase::createNewElements(percept::PerceptMesh
&eMesh, NodeRegistry &nodeRegistry, stk::mesh::Entity &element,
NewSubEntityNodesType &new_sub_entity_nodes, vector< stk::mesh::Entity
* >::iterator &element_pool, stk::mesh::FieldBase
*proc_rank_field=0)=0

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::set_parent_child_relations "void
stk::adapt::UniformRefinerPatternBase::set_parent_child_relations(percept::PerceptMesh
&eMesh, stk::mesh::Entity &old_owning_elem, stk::mesh::Entity
&newElement, unsigned ordinal, unsigned *numChild=0)

if numChild is passed in as non-null, use that value, else use
getNumNewElemPerElem() as size of child vector ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::interpolateElementFields "void
stk::adapt::UniformRefinerPatternBase::interpolateElementFields(percept::PerceptMesh
&eMesh, stk::mesh::Entity &old_owning_elem, stk::mesh::Entity
&newElement) ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::findSideRelations "bool
stk::adapt::UniformRefinerPatternBase::findSideRelations(percept::PerceptMesh
&eMesh, stk::mesh::Entity *parent, stk::mesh::Entity *child)

given a new element (child) that is a child of an original element
(parent), look at parent's side to elem relations and from the
children of the element, choose an element to connect the new side to
(using connectSides) ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::connectSides "bool
stk::adapt::UniformRefinerPatternBase::connectSides(percept::PerceptMesh
&eMesh, stk::mesh::Entity *element, stk::mesh::Entity *side_elem) ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::setSubPatterns "virtual void
stk::adapt::UniformRefinerPatternBase::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::fixSurfaceAndEdgeSetNamesMap "virtual StringStringMap
stk::adapt::UniformRefinerPatternBase::fixSurfaceAndEdgeSetNamesMap()

for i/o to work properly, supply string replacements such as for
hex-->tet breaking, you would supply \"quad\"-->\"tri\" etc. string
maps ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getFromTopoPartName "virtual
std::string
stk::adapt::UniformRefinerPatternBase::getFromTopoPartName()=0

provided by this class
------------------------------------------------------------------------------------------------------------------------
";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getToTopoPartName "virtual
std::string
stk::adapt::UniformRefinerPatternBase::getToTopoPartName()=0 ";

%feature("docstring")  stk::adapt::UniformRefinerPatternBase::getName
"virtual std::string stk::adapt::UniformRefinerPatternBase::getName()
";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getToParts "stk::mesh::PartVector&
stk::adapt::UniformRefinerPatternBase::getToParts() ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getFromParts "stk::mesh::PartVector&
stk::adapt::UniformRefinerPatternBase::getFromParts() ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getAppendConvertString "const
std::string&
stk::adapt::UniformRefinerPatternBase::getAppendConvertString() ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getAppendOriginalString "const
std::string&
stk::adapt::UniformRefinerPatternBase::getAppendOriginalString() ";

%feature("docstring")  stk::adapt::UniformRefinerPatternBase::setToOne
"void stk::adapt::UniformRefinerPatternBase::setToOne(std::vector<
NeededEntityType > &needed_entities)

utilities --------- sets the needed number of nodes on each sub-entity
to 1 - this is just a helper - in general, edges and faces have 1 new
node for linear elements, and multiple new nodes in the case of
quadratic elements ";

%feature("docstring")  stk::adapt::UniformRefinerPatternBase::midPoint
"double* stk::adapt::UniformRefinerPatternBase::midPoint(const double
*p1, const double *p2, int spatialDim, double *x) ";

%feature("docstring")
stk::adapt::UniformRefinerPatternBase::getCentroid "double*
stk::adapt::UniformRefinerPatternBase::getCentroid(double *pts[], int
len, int spatialDim, double *x) ";


// File: classstk_1_1adapt_1_1unit__tests_1_1UnitTestSupport.xml
%feature("docstring") stk::adapt::unit_tests::UnitTestSupport "";


// File: classstk_1_1adapt_1_1URP.xml
%feature("docstring") stk::adapt::URP "";

%feature("docstring")  stk::adapt::URP::getFromTypeKey "virtual
unsigned stk::adapt::URP< FromTopology, ToTopology >::getFromTypeKey()
";

%feature("docstring")  stk::adapt::URP::getFromTopology "virtual
const CellTopologyData* stk::adapt::URP< FromTopology, ToTopology
>::getFromTopology() ";

%feature("docstring")  stk::adapt::URP::getToTopology "virtual const
CellTopologyData* stk::adapt::URP< FromTopology, ToTopology
>::getToTopology() ";

%feature("docstring")  stk::adapt::URP::getFromTopoPartName "virtual
std::string stk::adapt::URP< FromTopology, ToTopology
>::getFromTopoPartName()

provided by this class
------------------------------------------------------------------------------------------------------------------------
";

%feature("docstring")  stk::adapt::URP::getToTopoPartName "virtual
std::string stk::adapt::URP< FromTopology, ToTopology
>::getToTopoPartName() ";

%feature("docstring")  stk::adapt::URP::getName "virtual std::string
stk::adapt::URP< FromTopology, ToTopology >::getName() ";

%feature("docstring")  stk::adapt::URP::~URP "virtual
stk::adapt::URP< FromTopology, ToTopology >::~URP() ";

%feature("docstring")  stk::adapt::URP::get_effective_topo "const
CellTopologyData* stk::adapt::URP< FromTopology, ToTopology
>::get_effective_topo(mesh::Part &part) ";

%feature("docstring")  stk::adapt::URP::setNeededParts "void
stk::adapt::URP< FromTopology, ToTopology
>::setNeededParts(percept::PerceptMesh &eMesh, BlockNamesType
block_names_ranks, bool sameTopology=true) ";

%feature("docstring")  stk::adapt::URP::change_entity_parts "void
stk::adapt::URP< FromTopology, ToTopology
>::change_entity_parts(percept::PerceptMesh &eMesh, stk::mesh::Entity
&old_owning_elem, stk::mesh::Entity &newElement) ";


// File: classstk_1_1adapt_1_1URP1.xml
%feature("docstring") stk::adapt::URP1 "

Utility intermediate base class providing more support for standard
refinement operations
------------------------------------------------------------------------------------------------------------------------

C++ includes: UniformRefinerPattern.hpp ";


// File: classstk_1_1adapt_1_1URP__Heterogeneous__3D.xml
%feature("docstring") stk::adapt::URP_Heterogeneous_3D "";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::URP_Heterogeneous_3D "stk::adapt::URP_Heterogeneous_3D::URP_Heterogeneous_3D(percept::PerceptMesh
&eMesh, BlockNamesType block_names=BlockNamesType())

setNeededParts(eMesh, block_names, true); ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::~URP_Heterogeneous_3D "stk::adapt::URP_Heterogeneous_3D::~URP_Heterogeneous_3D() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::setSubPatterns "void
stk::adapt::URP_Heterogeneous_3D::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")  stk::adapt::URP_Heterogeneous_3D::doBreak "virtual void stk::adapt::URP_Heterogeneous_3D::doBreak() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::getFromTypeKey "virtual unsigned
stk::adapt::URP_Heterogeneous_3D::getFromTypeKey() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::getFromTopoPartName "virtual
std::string stk::adapt::URP_Heterogeneous_3D::getFromTopoPartName()

provided by this class
------------------------------------------------------------------------------------------------------------------------
";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::getToTopoPartName "virtual
std::string stk::adapt::URP_Heterogeneous_3D::getToTopoPartName() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::getFromTopology "virtual const
CellTopologyData* stk::adapt::URP_Heterogeneous_3D::getFromTopology()
";

%feature("docstring")  stk::adapt::URP_Heterogeneous_3D::getToTopology
"virtual const CellTopologyData*
stk::adapt::URP_Heterogeneous_3D::getToTopology() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::fillNeededEntities "void
stk::adapt::URP_Heterogeneous_3D::fillNeededEntities(std::vector<
NeededEntityType > &needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::getNumNewElemPerElem "virtual
unsigned stk::adapt::URP_Heterogeneous_3D::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_3D::createNewElements "void
stk::adapt::URP_Heterogeneous_3D::createNewElements(percept::PerceptMesh
&eMesh, NodeRegistry &nodeRegistry, stk::mesh::Entity &element,
NewSubEntityNodesType &new_sub_entity_nodes, vector< stk::mesh::Entity
* >::iterator &element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1URP__Heterogeneous__Enrich__3D.xml
%feature("docstring") stk::adapt::URP_Heterogeneous_Enrich_3D "";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::URP_Heterogeneous_Enrich_3D "stk::adapt::URP_Heterogeneous_Enrich_3D::URP_Heterogeneous_Enrich_3D(percept::PerceptMesh
&eMesh, BlockNamesType block_names=BlockNamesType()) ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::~URP_Heterogeneous_Enrich_3D
"stk::adapt::URP_Heterogeneous_Enrich_3D::~URP_Heterogeneous_Enrich_3D()
";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::setSubPatterns "void
stk::adapt::URP_Heterogeneous_Enrich_3D::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::getFromTopoPartName "virtual
std::string
stk::adapt::URP_Heterogeneous_Enrich_3D::getFromTopoPartName()

provided by this class
------------------------------------------------------------------------------------------------------------------------
";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::getToTopoPartName "virtual
std::string
stk::adapt::URP_Heterogeneous_Enrich_3D::getToTopoPartName() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::doBreak "virtual void
stk::adapt::URP_Heterogeneous_Enrich_3D::doBreak() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::getFromTypeKey "virtual
unsigned stk::adapt::URP_Heterogeneous_Enrich_3D::getFromTypeKey() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::getFromTopology "virtual
const CellTopologyData*
stk::adapt::URP_Heterogeneous_Enrich_3D::getFromTopology() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::getToTopology "virtual const
CellTopologyData*
stk::adapt::URP_Heterogeneous_Enrich_3D::getToTopology() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::fillNeededEntities "void
stk::adapt::URP_Heterogeneous_Enrich_3D::fillNeededEntities(std::vector<
NeededEntityType > &needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::getNumNewElemPerElem "virtual unsigned
stk::adapt::URP_Heterogeneous_Enrich_3D::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_Enrich_3D::createNewElements "void
stk::adapt::URP_Heterogeneous_Enrich_3D::createNewElements(percept::PerceptMesh
&eMesh, NodeRegistry &nodeRegistry, stk::mesh::Entity &element,
NewSubEntityNodesType &new_sub_entity_nodes, vector< stk::mesh::Entity
* >::iterator &element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1adapt_1_1URP__Heterogeneous__QuadraticRefine__3D.xml
%feature("docstring") stk::adapt::URP_Heterogeneous_QuadraticRefine_3D
"";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::URP_Heterogeneous_QuadraticRefine_3D
"stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::URP_Heterogeneous_QuadraticRefine_3D(percept::PerceptMesh
&eMesh, BlockNamesType block_names=BlockNamesType()) ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::~URP_Heterogeneous_QuadraticRefine_3D
"stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::~URP_Heterogeneous_QuadraticRefine_3D()
";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::setSubPatterns "void
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::setSubPatterns(std::vector<
UniformRefinerPatternBase * > &bp, percept::PerceptMesh &eMesh)

optionally overridden (must be overridden if sidesets are to work
properly) to provide info on which sub pattern should be used to
refine side sets (and edge sets)

default is only this pattern ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getFromTopoPartName
"virtual std::string
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getFromTopoPartName()

provided by this class
------------------------------------------------------------------------------------------------------------------------
";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getToTopoPartName "virtual std::string
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getToTopoPartName()
";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::doBreak "virtual
void stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::doBreak() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getFromTypeKey "virtual unsigned
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getFromTypeKey() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getFromTopology "virtual const CellTopologyData*
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getFromTopology() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getToTopology "virtual const CellTopologyData*
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getToTopology() ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::fillNeededEntities "void
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::fillNeededEntities(std::vector<
NeededEntityType > &needed_entities)

must be provided by derived classes
------------------------------------------------------------------------------------------------------------------------
supplies the ranks of the sub entities needed during refinement (eg.
m_eMesh.face_rank(), m_eMesh.edge_rank(),..) 10/02/10 and the number
of nodes needed for each sub entity ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getNumNewElemPerElem
"virtual unsigned
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::getNumNewElemPerElem()

supply the number of new elements per element during refinement ";

%feature("docstring")
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::createNewElements "void
stk::adapt::URP_Heterogeneous_QuadraticRefine_3D::createNewElements(percept::PerceptMesh
&eMesh, NodeRegistry &nodeRegistry, stk::mesh::Entity &element,
NewSubEntityNodesType &new_sub_entity_nodes, vector< stk::mesh::Entity
* >::iterator &element_pool, stk::mesh::FieldBase *proc_rank_field=0)

given the node database ( NodeRegistry), and the newly created nodes,
and an iterator for the elements in the element pool, create all new
sub-elements of the refined element ";


// File: classstk_1_1percept_1_1Util.xml
%feature("docstring") stk::percept::Util "";


// File: classVector.xml
%feature("docstring") Vector "";

%feature("docstring")  Vector::init "void Vector< T >::init(double t,
RefElem i) ";

%feature("docstring")  Vector::init "void Vector< T, I >::init(T t, I
i) ";

%feature("docstring")  Vector::value "T Vector< T >::value() ";


// File: classvector.xml
%feature("docstring") vector "";

%feature("docstring")  vector::init "void vector< T >::init(T t) ";

%feature("docstring")  vector::value "T vector< T >::value() ";

%feature("docstring")  vector::init "void vector< double
>::init(double t) ";


// File: classvector1.xml
%feature("docstring") vector1 "";

%feature("docstring")  vector1::init "void vector1< T, V >::init(T t)
";

%feature("docstring")  vector1::init "void vector1< double, vector
>::init(double d) ";


// File: classvector2.xml
%feature("docstring") vector2 "";

%feature("docstring")  vector2::init "void vector2< T, V >::init(T t)
";

%feature("docstring")  vector2::init "void vector2< double, vector
>::init(double d) ";


// File: classVectorB.xml
%feature("docstring") VectorB "";

%feature("docstring")  VectorB::init "void VectorB< T, RefElem
>::init(T t, RefElem i) ";

%feature("docstring")  VectorB::value "T VectorB< T, RefElem
>::value() ";


// File: classstk_1_1percept_1_1Verifier.xml
%feature("docstring") stk::percept::Verifier "";

%feature("docstring")  stk::percept::Verifier::Verifier "stk::percept::Verifier::Verifier() ";

%feature("docstring")  stk::percept::Verifier::verify "void
stk::percept::Verifier::verify(int argc, char **argv) ";

%feature("docstring")  stk::percept::Verifier::process_options "void
stk::percept::Verifier::process_options(RunEnvironment &re) ";


// File: classstk_1_1percept_1_1WedgeFixture.xml
%feature("docstring") stk::percept::WedgeFixture "";

%feature("docstring")  stk::percept::WedgeFixture::createMesh "mesh::BulkData*
stk::percept::WedgeFixture::createMesh(stk::ParallelMachine
parallel_machine, unsigned n_nodes_x, unsigned n_nodes_y, unsigned
n_nodes_z, double xmin, double xmax, double ymin, double ymax, double
zmin, double zmax, std::string output_filename) ";

%feature("docstring")  stk::percept::WedgeFixture::getMetaData "mesh::fem::FEMMetaData* stk::percept::WedgeFixture::getMetaData() ";

%feature("docstring")
stk::percept::WedgeFixture::createBulkAfterMetaCommit "void
stk::percept::WedgeFixture::createBulkAfterMetaCommit(stk::ParallelMachine
parallel_machine) ";

%feature("docstring")  stk::percept::WedgeFixture::createFixedSizeMesh
"void
stk::percept::WedgeFixture::createFixedSizeMesh(stk::ParallelMachine
parallel_machine, std::string output_filename) ";


// File: classstk_1_1percept_1_1IntrepidManager_1_1WeightedMeasure.xml
%feature("docstring") stk::percept::IntrepidManager::WeightedMeasure "

weights multiplied by Jacobian det at cubature points ([C], [P])

C++ includes: IntrepidManager.hpp ";

%feature("docstring")
stk::percept::IntrepidManager::WeightedMeasure::WeightedMeasure "stk::percept::IntrepidManager::WeightedMeasure::WeightedMeasure(IM
&im)

([C], [P]) ";


// File: classstk_1_1percept_1_1unit__tests_1_1XF1.xml
%feature("docstring") stk::percept::unit_tests::XF1 "";


// File: namespace@303.xml


// File: namespace@313.xml


// File: namespaceIntrepid.xml


// File: namespacemoab.xml
%feature("docstring")  moab::error "static void moab::error(int i) ";

%feature("docstring")  moab::compute_edge_length_squared "double
moab::compute_edge_length_squared(double *c0, double *c1) ";


// File: namespacepercept.xml


// File: namespaceshards.xml


// File: namespacestd.xml


// File: namespacestk.xml


// File: namespacestk_1_1adapt.xml
%feature("docstring")
stk::adapt::Elem::dummyAllocateMethodForMacBuildWarning "void
stk::adapt::dummyAllocateMethodForMacBuildWarning() ";

%feature("docstring")  stk::adapt::Elem::contains "bool
stk::adapt::contains(STD_Set &set, Key key) ";

%feature("docstring")  stk::adapt::Elem::BOOST_STATIC_ASSERT "stk::adapt::BOOST_STATIC_ASSERT(sizeof(uint64_t)==sizeof(size_t)) ";

%feature("docstring")  stk::adapt::Elem::MegaByte "static double
stk::adapt::MegaByte(MemorySizeType x) ";

%feature("docstring")  stk::adapt::Elem::memory_dump "static
MemorySizeType stk::adapt::memory_dump(int dump_level, const
stk::ParallelMachine &comm, stk::mesh::BulkData &bulkData,
NodeRegistry *node_reg, std::string msg) ";

%feature("docstring")  stk::adapt::Elem::test_memory "void
stk::adapt::test_memory(percept::PerceptMesh &eMesh, MemorySizeType
n_elements, MemorySizeType n_nodes) ";

%feature("docstring")  stk::adapt::Elem::checkInput "static void
stk::adapt::checkInput(std::string option, std::string value,
std::string allowed_values, RunEnvironment &run_environment) ";

%feature("docstring")  stk::adapt::Elem::print_simple_usage "static
void stk::adapt::print_simple_usage(int argc, char **argv) ";

%feature("docstring")  stk::adapt::Elem::adapt_main "int
stk::adapt::adapt_main(int argc, char **argv) ";

%feature("docstring")  stk::adapt::Elem::adapt_main_full_options "int
stk::adapt::adapt_main_full_options(int argc, char **argv) ";

%feature("docstring")  stk::adapt::Elem::check_for_simple_options "static int stk::adapt::check_for_simple_options(int argc, char **argv)
";

%feature("docstring")  stk::adapt::Elem::adapt_main_simple_options "int stk::adapt::adapt_main_simple_options(int argc_in, char **argv_in)
";

%feature("docstring")  stk::adapt::Elem::dump_args "static void
stk::adapt::dump_args(int argc, char **argv) ";

%feature("docstring")  stk::adapt::Elem::doPrintSizes "static void
stk::adapt::doPrintSizes() ";

%feature("docstring")  stk::adapt::Elem::getChildVectorPtr "static
const SameRankRelationValue*
stk::adapt::getChildVectorPtr(SameRankRelation &repo,
stk::mesh::Entity *parent) ";

%feature("docstring")  stk::adapt::Elem::get_random_sequence "static
std::vector<unsigned> stk::adapt::get_random_sequence(int min, int
max, int num) ";

%feature("docstring")  stk::adapt::Elem::exact_hessian "void
stk::adapt::exact_hessian(const int id, const double *xyz, double
*hess, const int spatial_dim) ";

%feature("docstring")  stk::adapt::Elem::interp_nodal_hessian "void
stk::adapt::interp_nodal_hessian(const int hess_id, PerceptMesh
&eMesh, VectorFieldType *nodal_hessian_field) ";

%feature("docstring")  stk::adapt::Elem::exact_nodal_solution "void
stk::adapt::exact_nodal_solution(const double *xyz, double *field,
const int spatial_dim) ";

%feature("docstring")  stk::adapt::Elem::triangle_area "double
stk::adapt::triangle_area(const std::vector< double > &nodal_coords)
";

%feature("docstring")  stk::adapt::Elem::compute_elem_mesh_size_ratio
"void stk::adapt::compute_elem_mesh_size_ratio(PerceptMesh &eMesh,
ScalarFieldType *elem_ratio_field, const double &global_error_tol) ";


// File: namespacestk_1_1adapt_1_1Elem.xml
%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getBasicCellTopology "CellTopology stk::adapt::Elem::getBasicCellTopology(const char *name)
";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getBasicCellTopology "CellTopology stk::adapt::Elem::getBasicCellTopology(TopologyId id) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getCellTopologyId "TopologyId
stk::adapt::Elem::getCellTopologyId(const CellTopology &cell_topology)
";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::edgeCellTopology "CellTopology stk::adapt::Elem::edgeCellTopology(const
Elem::CellTopology &cell_topology, UInt ordinal) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::faceCellTopology "CellTopology stk::adapt::Elem::faceCellTopology(const
Elem::CellTopology &cell_topology, UInt ordinal) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::sideCellTopology "CellTopology stk::adapt::Elem::sideCellTopology(const
Elem::CellTopology &cell_topology, UInt ordinal)

Query 2D edge or 3D face topologies. ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::isElement "bool
stk::adapt::Elem::isElement(const Elem::CellTopology &cell_topology,
unsigned spatial_dimension) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::isSolidElement "bool
stk::adapt::Elem::isSolidElement(const Elem::CellTopology
&cell_topology, unsigned spatial_dimension) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::isShellElement "bool
stk::adapt::Elem::isShellElement(const Elem::CellTopology
&cell_topology) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::isRodElement "bool
stk::adapt::Elem::isRodElement(const Elem::CellTopology
&cell_topology) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::isParticleElement "bool
stk::adapt::Elem::isParticleElement(const Elem::CellTopology
&cell_topology) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getNodesOfEdge "const
unsigned * stk::adapt::Elem::getNodesOfEdge(const Elem::CellTopology
&cell_topology, unsigned edge) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getNodesOfFace "const
unsigned * stk::adapt::Elem::getNodesOfFace(const Elem::CellTopology
&cell_topology, unsigned face) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getNodesOfSide "const
unsigned * stk::adapt::Elem::getNodesOfSide(const Elem::CellTopology
&cell_topology, unsigned side) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getEdgeNode "int
stk::adapt::Elem::getEdgeNode(const Elem::CellTopology &cell_topology,
unsigned edge, unsigned node_of_edge) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getFaceNode "int
stk::adapt::Elem::getFaceNode(const Elem::CellTopology &cell_topology,
unsigned face, unsigned node_of_face) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getSideNode "int
stk::adapt::Elem::getSideNode(const Elem::CellTopology &cell_topology,
unsigned side, unsigned node_of_side) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getFaceEdge "int
stk::adapt::Elem::getFaceEdge(const Elem::CellTopology &cell_topology,
unsigned face, unsigned edge_of_face) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::nodeCellTopology "CellTopology stk::adapt::Elem::nodeCellTopology(const
Elem::CellTopology &cell_topology, UInt ordinal) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::isCellTopologySubsetOf "bool
stk::adapt::Elem::isCellTopologySubsetOf(const Elem::CellTopology
&cell_topology, const Elem::CellTopology &richer) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getParametricDimension "unsigned stk::adapt::Elem::getParametricDimension(const
Elem::CellTopology &cell_topology) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::findReversePermutation "int
stk::adapt::Elem::findReversePermutation(const CellTopologyData &top,
int permutation_ord) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getCellTopology "CellTopology
stk::adapt::Elem::getCellTopology() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getRefinementTopology "const
RefinementTopology * stk::adapt::Elem::getRefinementTopology(const
Elem::CellTopology &cell_topology) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getRefinementEdgeNode "const
UInt * stk::adapt::Elem::getRefinementEdgeNode(const
Elem::CellTopology &cell_topology, UInt edge) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getRefinementFaceNode "const
UInt * stk::adapt::Elem::getRefinementFaceNode(const
Elem::CellTopology &cell_topology, UInt face) ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::getRefinementEdgePermutation "const UInt * stk::adapt::Elem::getRefinementEdgePermutation(const
Elem::CellTopology &cell_topology, UInt permutation_ordinal) ";


// File: namespacestk_1_1adapt_1_1Elem_1_1@398.xml


// File: namespacestk_1_1adapt_1_1Elem_1_1@404.xml


// File: namespacestk_1_1adapt_1_1Elem_1_1StdMeshObjTopologies.xml
%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::node "const MeshObjTopology * stk::adapt::Elem::StdMeshObjTopologies::node()
";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::point "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::point() ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::line "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::line(UInt eclass, UInt nnode)
";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::tri "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::tri(UInt eclass, UInt nnode)

2 PARENT Linear 3-Node Triangle Element Nodes o (SPACE_DIM = 2!) /
\\\\ / \\\\ (PARENT) Linear 3-Node Triangle Edge Node Map: / \\\\ /
\\\\ { {0, 1}, {1, 2}, {2, 0} }; / \\\\ / \\\\ / \\\\ o---------------
o 0 1

After refinement:

2 CHILD Linear 3-Node Triangle Element Nodes o (new nodes = *) / \\\\
/ \\\\ / \\\\ 5 *-------* 4 / \\\\ / \\\\ / \\\\ / \\\\ / \\\\ / \\\\
o-------*-------o 0 3 1

| CHILD Linear 3-Node Triangle Element Node Maps: | | | static const
UInt child_0[] = { 0, 3, 5 }; | static const UInt child_1[] = { 3, 1,
4 }; | static const UInt child_2[] = { 5, 4, 2 }; | static const UInt
child_3[] = { 4, 5, 3 }; | 2 PARENT 6-Node Triangle Object Nodes o /
\\\\ / \\\\ (PARENT) 6-Node Triangle Object Edge Node Map: / \\\\ 5 o
o 4 { {0, 1, 3}, {1, 2, 4}, {2, 0, 5} }; / \\\\ / \\\\ / \\\\ o-------
o-------o 0 3 1

After refinement:

2 CHILD 6-Node Triangle Object Nodes o (new nodes = *) / \\\\ 10 * * 9
/ 14 \\\\ 5 o---*---o 4 / \\\\ / \\\\ 11 * 12* *13 * 8 / \\\\ / \\\\
o---*---o---*---o 0 6 3 7 1

| CHILD 6-Node Triangle Object Node Maps: | | static const UInt
child_0[] = { 0, 3, 5, 6, 12, 11 }; | static const UInt child_1[] = {
3, 1, 4, 7, 8, 13 }; | static const UInt child_2[] = { 5, 4, 2, 14, 9,
10 }; | static const UInt child_3[] = { 4, 5, 3, 14, 12, 13 }; |

Refined 6-Node Triangle Object PERMUTATION Node Maps:

Rotation Polarity 0 1 { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
14 }; 0 0 { 0, 2, 1, 5, 4, 3, 11, 10, 9, 8, 7, 6, 12, 14, 13 }; 1 1 {
2, 0, 1, 5, 3, 4, 10, 11, 6, 7, 8, 9, 14, 12, 13 }; 1 0 { 2, 1, 0, 4,
3, 5, 9, 8, 7, 6, 11, 10, 14, 13, 12 }; 2 1 { 1, 2, 0, 4, 5, 3, 8, 9,
10, 11, 6, 7, 13, 14, 12 }; 2 0 { 1, 0, 2, 3, 5, 4, 7, 6, 11, 10, 9, 8
13, 12, 14 }; ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::tri4 "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::tri4(UInt eclass)

2 PARENT 4-Node Triangle Object Nodes o / \\\\ / \\\\ (PARENT) 4-Node
Triangle Object Edge Node Map: / \\\\ / \\\\ { {0, 1}, {1, 2}, {2, 0}
}; / o \\\\ / 3 \\\\ / \\\\ o---------------o 0 1

After refinement:

2 CHILD 4-Node Triangle Object Nodes o (new nodes = *) / \\\\ / 9 \\\\
/ * \\\\ 6 *-------* 5 / \\\\ o / \\\\ / 7 \\\\ 3 / 8 \\\\ / * \\\\ /
* \\\\ o-------*-------o 0 4 1

| CHILD 4-Node Triangle Object Node Maps: | | static const UInt
child_0[] = { 0, 4, 6, 7 }; | static const UInt child_1[] = { 4, 1, 5,
8 }; | static const UInt child_2[] = { 6, 5, 2, 9 }; | static const
UInt child_3[] = { 5, 6, 4, 3 }; |

4-Node Triangle Object PERMUTATION Node Maps:

Original Parent 4-Node Triangle Object: Rotation Polarity 0 1 { 0, 1,
2, 3 } 0 0 { 0, 2, 1, 3 } 1 1 { 2, 0, 1, 3 } 1 0 { 2, 1, 0, 3 } 2 1 {
1, 2, 0, 3 } 2 0 { 1, 0, 2, 3 }

After Refinement and using child node numbering: Rotation Polarity

0 1 { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }; 0 0 { 0, 2, 1, 3, 6, 5, 4, 7, 9,
8 }; 1 1 { 2, 0, 1, 3, 6, 4, 5, 9, 7, 8 }; 1 0 { 2, 1, 0, 3, 6, 5, 4,
9, 8, 7 }; 2 1 { 1, 2, 0, 3, 5, 6, 4, 8, 9, 7 }; 2 0 { 1, 0, 2, 3, 4,
6, 5, 8, 7, 9 }; ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::quad "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::quad(UInt eclass, UInt nnode)

3 2 PARENT Linear 4-Node Quadrilateral Element Nodes o---------------o
(SPACE_DIM = 2!) | | | | | | | | | | (PARENT) Linear 4-Node
Quadrilateral | | Element Edge Node Map: | | | | { {0, 1}, {1, 2}, {2,
3} {3, 0} }; | | o---------------o 0 1

After refinement:

3 6 2 CHILD Linear 4-Node Quadrilateral Element Nodes
o-------*-------o (SPACE_DIM = 2!) (new nodes = *) | | | | | | | | | |
8| | 7*-------*-------*5 | | | | | | | | | | | | o-------*-------o 0 4
1 | CHILD Linear 4-Node Quadrilateral Element Node Maps: | | | Element
0: childNodeMap[0] = { 0, 4, 8, 7 } | Element 1: childNodeMap[1] = {
4, 1, 5, 8 } | Element 2: childNodeMap[2] = { 8, 5, 2, 6 } | Element
3: childNodeMap[3] = { 7, 8, 6, 3 } |

New ref topo info Quad4 ------------------

{Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc,
param-coord} { {0, 0, 0, 0, 1, {0.0, 0.0, 0.0} }, {1, 0, 1, 0, 1,
{1.0, 0.0, 0.0} }, {2, 0, 2, 0, 1, {1.0, 1.0, 0.0} }, {3, 0, 3, 0, 1,
{0.0, 1.0, 0.0} }, {4, 1, 0, 0, 1, {0.5, 0.0, 0.0} }, {5, 1, 1, 0, 1,
{1.0, 0.5, 0.0} }, {6, 1, 2, 0, 1, {0.5, 1.0, 0.0} }, {7, 1, 3, 0, 1,
{0.0, 0.5, 0.0} }, {8, 2, 0, 0, 1, {0.5, 0.5, 0.0} } }

Refined Linear 4-Node Quadrilateral Element PERMUTATION Node Maps:

Rotation Polarity 0 1 { 0, 1, 2, 3; 4, 5, 6, 7, 8 } 0 0 { 0, 3, 2, 1;
7, 6, 5, 4, 8 } 1 1 { 3, 0, 1, 2; 7, 4, 5, 6, 8 } 1 0 { 3, 2, 1, 0; 6,
5, 4, 7, 8 } 2 1 { 2, 3, 0, 1; 6, 7, 4, 5, 8 } 2 0 { 2, 1, 0, 3; 5, 4,
7, 6, 8 } 3 1 { 1, 2, 3, 0; 5, 6, 7, 4, 8 } 3 0 { 1, 0, 3, 2; 4, 7, 6,
5, 8 } 3 6 2 PARENT 9-Node Quadrilateral Object Nodes o-------o-------
o | | | | | 8 | 7 o o o 5 (PARENT) 9-Node Quadrilateral Object's | |
Edge Node Map: | | | | { {0, 1, 4}, {1, 2, 5}, {2, 3, 6} {3, 0, 7} };
o-------o-------o 0 4 1

After refinement:

3 14 6 13 2 CHILD 9-Node Quadrilateral Object Nodes
o----*----o----*----o (new nodes = *) | | | | 24 | 23 | 15* * *19 *
*12 | | | | 8| 18 | 7 o----*----o----*----o 5 | 20 | | | | | 16* * 17*
* *11 | 21 | 22 | | | | o----*----o----*----o 0 9 4 10 1

CHILD 9-Node Quadrilateral Object Node Maps: | | | Object 0:
childNodeMap[0] = { 0, 4, 8, 7, 9, 17, 20, 16; 21 } | Object 1:
childNodeMap[1] = { 4, 1, 5, 8, 10, 11, 18, 17; 22 } | Object 2:
childNodeMap[2] = { 8, 5, 2, 6, 18, 12, 13, 19; 23 } | Object 3:
childNodeMap[3] = { 7, 8, 6, 3, 20, 19, 14, 15; 24 } | New ref topo
info Quad9 ------------------

{Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc,
param-coord} { {0, 0, 0, 0, 1, {0.0, 0.0, 0.0} }, {1, 0, 1, 0, 1,
{1.0, 0.0, 0.0} }, {2, 0, 2, 0, 1, {1.0, 1.0, 0.0} }, {3, 0, 3, 0, 1,
{0.0, 1.0, 0.0} }, {4, 1, 0, 0, 1, {0.5, 0.0, 0.0} }, {5, 1, 1, 0, 1,
{1.0, 0.5, 0.0} }, {6, 1, 2, 0, 1, {0.5, 1.0, 0.0} }, {7, 1, 3, 0, 1,
{0.0, 0.5, 0.0} }, {8, 2, 0, 8, 9, {0.5, 0.5, 0.0} },

{9, 1, 0, 1, 3, {0.25, 0.00, 0.00} }, {10, 1, 0, 2, 3, {0.75, 0.00,
0.00} }, {11, 1, 1, 1, 3, {1.00, 0.25, 0.00} }, {12, 1, 1, 2, 3,
{1.00, 0.75, 0.00} }, {13, 1, 2, 1, 3, {0.75, 1.00, 0.00} }, {14, 1,
2, 2, 3, {0.25, 1.00, 0.00} }, {15, 1, 3, 1, 3, {0.00, 0.75, 0.00} },
{16, 1, 3, 2, 3, {0.00, 0.25, 0.00} }

{17, 2, 0, 4, 9, {0.50, 0.25, 0.00} }, {18, 2, 0, 5, 9, {0.75, 0.50,
0.00} }, {19, 2, 0, 6, 9, {0.50, 0.75, 0.00} }, {20, 2, 0, 7, 9,
{0.25, 0.50, 0.00} },

{21, 2, 0, 0, 9, {0.25, 0.25, 0.00} }, {22, 2, 0, 1, 9, {0.75, 0.25,
0.00} }, {23, 2, 0, 2, 9, {0.75, 0.75, 0.00} }, {24, 2, 0, 3, 9,
{0.25, 0.75, 0.00} },

}

Refined 9-Node Quadrilateral Object PERMUTATION Node Maps:

Rotation Polarity 0 1 { 0, 1, 2, 3, 4, 5, 6, 7; 8, 9, 10, 11, 12, 13,
14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 } 0 0 { 0, 3, 2, 1, 7, 6,
5, 4; 8, 16, 15, 14, 13, 12, 11, 10, 9, 20, 19, 18, 17, 21, 24, 23, 22
} 1 1 { 3, 0, 1, 2, 7, 4, 5, 6; 8, 15, 16, 9, 10, 11, 12, 13, 14, 20,
17, 18, 19, 24, 21, 22, 23 } 1 0 { 3, 2, 1, 0, 6, 5, 4, 7; 8, 14, 13,
12, 11, 10, 9, 16, 15, 19, 18, 17, 20, 24, 23, 22, 21 } 1 1 { 1, 2, 3,
0, 5, 6, 7, 4; 8, 11, 12, 13, 14, 15, 16, 9, 10, 18, 19, 20, 17, 22,
23, 24, 21 } 1 0 { 1, 0, 3, 2, 4, 7, 6, 5; 8, 10, 9, 16, 15, 14, 13,
12, 11, 17, 20, 19, 18, 22, 21, 24, 23 } 3 1 { 2, 3, 0, 1, 6, 7, 4, 5;
8, 13, 14, 15, 16, 9, 10, 11, 12, 19, 20, 17, 18, 23, 24, 21, 22 } 3 0
{ 2, 1, 0, 3, 5, 4, 7, 6; 8, 12, 11, 10, 9, 16, 15, 14, 13, 18, 17,
20, 19, 23, 22, 21, 24 } ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::tet "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::tet(UInt nnode)

PARENT 4-Node Tetrahedron Object Nodes 3 o /|\\\\ / | \\\\ (PARENT)
4-Node Tetrahedron Object / | \\\\ Edge Node Map: / | \\\\ / | \\\\ {
{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} }; 0 o-----|-----o 2
\\\\ | / \\\\ | / \\\\ | / \\\\ | / \\\\|/ o 1

After refinement (new nodes = *):

3 o /|\\\\ / | \\\\ 7 * | * 9 / | \\\\ / 6| \\\\ 0 o----*|-----o 2
\\\\ *8 / \\\\ | / 4 * | * 5 \\\\ | / \\\\|/ o 1

CHILD 4-Node Tetrahedron 3D Object Node Maps: | | static const UInt
child_0[] = { 0, 4, 6, 7 }; // srkenno 091410 fixed (used to be {0, 4,
8, 7} ) | static const UInt child_1[] = { 4, 1, 5, 8 }; | static const
UInt child_2[] = { 6, 5, 2, 9 }; | static const UInt child_3[] = { 7,
8, 9, 3 }; | static const UInt child_4[] = { 8, 7, 6, 4 }; | static
const UInt child_5[] = { 6, 9, 8, 5 }; | static const UInt child_6[] =
{ 9, 8, 7, 6 }; | static const UInt child_7[] = { 5, 6, 4, 8 }; |
PARENT 10-Node Tetrahedron Object Nodes 3 o /|\\\\ / | \\\\ 7 o | o 9
(PARENT) 10-Node Tetrahedron Object / | \\\\ Edge Node Map: / 6| \\\\
0 o----o|-----o 2 { {0, 1, 4}, {1, 2, 5}, {2, 0, 6}, \\\\ o8 / {0, 3,
7}, {1, 3, 8}, {2, 3, 9} }; \\\\ | / 4 o | o 5 \\\\ | / \\\\|/ o 1

After refinement (new nodes = *):

3 o /|\\\\ * | * / | \\\\ 7 o | o 9 / * \\\\ * | * / 6| \\\\ 0
o---*--o|--*----o 2 \\\\ o8 / * | * \\\\ | / 4 o | o 5 \\\\ * / * | *
\\\\|/ o 1

| // Child edge node tables | | static const UInt edge_0[] = { 0, 1,
4, 10, 11 }; | static const UInt edge_1[] = { 1, 2, 5, 12, 13 }; |
static const UInt edge_2[] = { 2, 0, 6, 14, 15 }; | static const UInt
edge_3[] = { 0, 3, 7, 16, 19 }; | static const UInt edge_4[] = { 1, 3,
8, 17, 20 }; | static const UInt edge_5[] = { 2, 3, 9, 18, 21 }; | |
// Child face node (cfn) tables: | // Local Face (LF) LF0 uses [0:5],
LF1 uses [6:11], LF2 uses [12:17], LF3 uses [18:23] | | static const
UInt cfn_0[] = {0, 4, 7, 10, 27, 16, 4, 6, 7, 31, 25, 27, 0, 7, 6, 16,
25, 15, 0, 6, 4, 15, 31, 10 }; | static const UInt cfn_1[] = {4, 1, 8,
11, 17, 34, 1, 5, 8, 12, 26, 17, 4, 8, 5, 34, 26, 20, 4, 5, 1, 30, 12,
11 }; | static const UInt cfn_2[] = {6, 5, 9, 24, 28, 33, 5, 2, 9, 13,
18, 28, 6, 9, 2, 33, 18, 14, 6, 2, 5, 14, 13, 24 }; | static const
UInt cfn_3[] = {7, 8, 3, 23, 20, 19, 8, 9, 3, 32, 21, 20, 7, 3, 9, 19,
21, 29, 7, 9, 8, 29, 32, 33 }; | static const UInt cfn_4[] = {8, 7, 4,
23, 27, 34, 7, 6, 4, 25, 31, 27, 8, 4, 6, 34, 31, 22, 8, 6, 7, 22, 25,
23 }; | static const UInt cfn_5[] = {6, 9, 5, 33, 28, 24, 9, 8, 5, 32,
26, 28, 6, 5, 8, 24, 26, 22, 6, 8, 9, 22, 32, 33 }; | static const
UInt cfn_6[] = {9, 8, 6, 32, 22, 33, 8, 7, 6, 23, 25, 22, 9, 6, 7, 33,
25, 29, 9, 7, 8, 29, 23, 32 }; | static const UInt cfn_7[] = {5, 6, 8,
24, 22, 26, 6, 4, 8, 31, 34, 22, 5, 8, 4, 26, 34, 30, 5, 4, 6, 30, 31,
24 }; |

Face #0: Face #1: Face #2:

3 3 3 o o o / \\\\ / \\\\ / \\\\ 19* *20 20* *21 21* *19 / 23 \\\\ /
32 \\\\ / 29 \\\\ < 7 o---*---o 8 8 o---*---o 9 9 o---*---o 7 \\\\ /
\\\\ / \\\\ / \\\\ / \\\\ / \\\\ / \\\\ \\\\ 16* 27* *34 *17 17* 26*
*28 *18 18* 33* *25 *16 | / \\\\ / \\\\ / \\\\ / \\\\ / \\\\ / \\\\ |
o---*---o---*---o o---*---o---*---o o---*---o---*---o # 0 10 4 11 1 1
12 5 13 2 2 14 6 15 0 # # \\\\ \\\\ \\\\ \\\\ --> -->

Face #3:

1 o / \\\\ 11* *12 / 30 \\\\ 4 o---*---o 5 / \\\\ / \\\\ 10* 31* *24
*13 / \\\\ / \\\\ o---*---o---*---o 0 15 6 14 2 # \\\\ \\\\ -->

Various Interior Faces of Children:

7 9 7 25 6 8 26 5 o o o---*---o o---*---o / \\\\ / \\\\ \\\\ / \\\\ /
23* *25 32* *33 27* *31 34* *30 / 22 \\\\ / 22 \\\\ \\\\ / \\\\ /
8o---*---o10 8o---*---o6 o o \\\\ / \\\\ / 4 4 26* *24 34* *31 \\\\ /
\\\\ / o o 5 4

9 9 o o / \\\\ / \\\\ 33* *28 29* *32 / \\\\ / \\\\ o---*---o
o---*---o 6 24 5 7 23 8

Edge node tables for refined 10-Node Tetrahedrons: | | static const
UInt edge_0[] = { 0, 1, 4, 10, 11 }; | static const UInt edge_1[] = {
1, 2, 5, 12, 13 }; | static const UInt edge_2[] = { 2, 0, 6, 14, 15 };
| static const UInt edge_3[] = { 0, 3, 7, 16, 19 }; | static const
UInt edge_4[] = { 1, 3, 8, 17, 20 }; | static const UInt edge_5[] = {
2, 3, 9, 18, 21 }; |

Face node tables for refined 10-Node Tetrahedrons: | | static const
UInt face_0[] = { 0, 1, 3, 4, 8, 7, 10, 11, 17, 20, 19, 16, 27, 34, 23
}; | static const UInt face_1[] = { 1, 2, 3, 5, 9, 8, 12, 13, 18, 21,
20, 17, 26, 28, 32 }; | static const UInt face_2[] = { 0, 3, 2, 7, 9,
6, 16, 19, 21, 18, 14, 15, 25, 29, 33 }; | static const UInt face_3[]
= { 0, 2, 1, 6, 5, 4, 15, 14, 13, 12, 11, 10, 31, 24, 30 }; |

CHILD 10-Node Tetrahedron Object Node Maps:

| | static const UInt child_0[] = { 0, 4, 8, 7, 10, 31, 15, 16, 27, 25
}; | static const UInt child_1[] = { 4, 1, 5, 8, 11, 12, 30, 34, 17,
26 }; | static const UInt child_2[] = { 6, 5, 2, 9, 24, 13, 14, 33,
28, 18 }; | static const UInt child_3[] = { 7, 8, 9, 3, 23, 32, 29,
19, 20, 21 }; | static const UInt child_4[] = { 8, 7, 6, 4, 23, 25,
22, 34, 27, 31 }; | static const UInt child_5[] = { 6, 9, 8, 5, 33,
32, 22, 24, 28, 26 }; | static const UInt child_6[] = { 9, 8, 7, 6,
32, 23, 29, 33, 22, 25 }; | static const UInt child_7[] = { 5, 6, 4,
8, 24, 31, 30, 26, 22, 34 }; | PARENT Semi-Linear 8-Node Tetrahedron
Nodes 3 (SPACE_DIM = 3!) o /|\\\\ / | \\\\ / | \\\\ \"Front faces\" /
| \\\\ Node 4 on 2D surface containing nodes 0, 1, 3 / | \\\\ Node 5
on 2D surface containing nodes 1, 2, 3 0 o-----|-----o 2 \\\\ 4 | + /
(PARENT) Semi-Linear 8-Node Tetrahedron \\\\ + | 5 / 3D Element Edge
Node Map: \\\\ | / \\\\ | / { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3},
{2, 3} } \\\\|/ o 1

3 o /|\\\\ \"Back faces (note mid-face-node does not follow face
ordering!)\" / | \\\\ Node 7 on 2D surface containing nodes 0, 3, 2 /
| \\\\ Node 6 on 2D surface containing nodes 0, 2, 1 / |+ \\\\ / |7
\\\\ 0 o-----|-----o 2 \\\\ 6| / \\\\ +| / \\\\ | / \\\\ | / \\\\|/ o
1

After refinement (new nodes = *):

3 o /|\\\\ / | \\\\ 11* | *13 / | \\\\ / 10| \\\\ 0 o----*|-----o 2
\\\\ *12 / \\\\ | / 8 * | * 9 \\\\ | / \\\\|/ o 1

| // Child edge node tables | | static const UInt edge_0[] = { 0, 1, 8
}; | static const UInt edge_1[] = { 1, 2, 9 }; | static const UInt
edge_2[] = { 2, 0, 10 }; | static const UInt edge_3[] = { 0, 3, 11 };
| static const UInt edge_4[] = { 1, 3, 12 }; | static const UInt
edge_5[] = { 2, 3, 13 }; |

| | // Child face node (cfn) tables: | // Local Face (LF) | LF0 uses
[0:3], LF1 uses [4:7], LF2 uses [8:11], LF3 uses [12:15] | | static
const UInt cfn_0[] = | { 0, 8, 11, 14, 8, 10, 11, 30, 0, 11, 10, 17,
0, 10, 8, 16 }; | static const UInt cfn_1[] = | { 8, 1, 12, 18, 1, 9,
12, 15, 8, 12, 9, 31, 8, 9, 1, 24 }; | static const UInt cfn_2[] = | {
10, 9, 13, 32, 9, 2, 13, 19, 10, 13, 2, 25, 10, 2, 9, 20 }; | static
const UInt cfn_3[] = | { 11, 12, 3, 22, 12, 13, 3, 23, 11, 3, 13, 21,
11, 13, 12, 33 }; | static const UInt cfn_4[] = | { 12, 11, 8, 4, 11,
10, 8, 30, 12, 8, 10, 29, 12, 10, 11, 26 }; | static const UInt
cfn_5[] = | { 10, 13, 9, 32, 13, 12, 9, 5, 10, 9, 12, 27, 10, 12, 13,
28 }; | static const UInt cfn_6[] = | { 13, 12, 10, 28, 12, 11, 10,
26, 13, 10, 11, 7, 13, 11, 12, 33 }; | static const UInt cfn_7[] = | {
9, 10, 12, 27, 10, 8, 12, 29, 9, 12, 8, 31, 9, 8, 10, 6 }; |

Face #0: Face #1: Face #2:

3 3 3 o o o / \\\\ / \\\\ / \\\\ / 16\\\\ / 19\\\\ / 24\\\\ / * \\\\ /
* \\\\ / * \\\\ < 11*-------*12 12*-------*13 13*-------*11 \\\\ /
\\\\ o / \\\\ / \\\\ o / \\\\ / \\\\ o / \\\\ \\\\ /14 \\\\ 4 /15 \\\\
/17 \\\\ 5 /18 \\\\ /25 \\\\ 7 / 23\\\\ | / * \\\\ / * \\\\ / * \\\\ /
* \\\\ / * \\\\ / * \\\\ | o-------*-------o o-------*-------o
o-------*-------o # 0 8 1 1 9 2 2 10 0 # # \\\\ \\\\ \\\\ \\\\ --> -->

Face #3:

1 o / \\\\ /22 \\\\ / * \\\\ 8*-------*9 / \\\\ o / \\\\ /20 \\\\ 6 /
21\\\\ / * \\\\ / * \\\\ o-------*-------o 0 10 2 # \\\\ \\\\ -->

Various Interior Faces of Children:

11 13 11 10 12 9 * * *-------* *-------* / \\\\ / \\\\ \\\\ * / \\\\ *
/ /26 \\\\ /28 \\\\  /  / / * \\\\ / * \\\\ \\\\ / \\\\ /
12*-------*10 10*-------*12 * * \\\\ * / \\\\ * / 8 8  /  / \\\\ /
\\\\ / * * 9 8

13 13 * * / \\\\ / \\\\ /32 \\\\ /33 \\\\ / * \\\\ / * \\\\ *-------*
*-------* 10 9 11 12

| // tet8 face node tables | | static const UInt t8_face_0[] = { 0, 1,
3, 4, 8, 12, 11, 14, 15, 16 }; | static const UInt t8_face_1[] = { 1,
2, 3, 5, 9, 13, 12, 17, 18, 19 }; | static const UInt t8_face_2[] = {
0, 3, 2, 7, 11, 13, 10, 23, 24, 25 }; | static const UInt t8_face_3[]
= { 0, 2, 1, 6, 10, 9, 8, 20, 21, 22 }; |

CHILD 8-Node Tetrahedron Object Node Maps:

| | static const UInt child_0[] = { 0, 8, 10, 11, 14, 30, 23, 20 }; |
static const UInt child_1[] = { 8, 1, 9, 12, 15, 17, 31, 22 }; |
static const UInt child_2[] = { 10, 9, 2, 13, 32, 18, 25, 21 }; |
static const UInt child_3[] = { 11, 12, 13, 3, 16, 19, 24, 33 }; |
static const UInt child_4[] = { 12, 11, 10, 8, 4, 30, 29, 26 }; |
static const UInt child_5[] = { 10, 13, 12, 9, 32, 5, 27, 28 }; |
static const UInt child_6[] = { 13, 12, 11, 10, 28, 26, 7, 33 }; |
static const UInt child_7[] = { 9, 10, 8, 12, 27, 29, 31, 6 }; | ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::hex "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::hex(UInt nnode)

PARENT Linear 8-Node Hexahedron Nodes 7 6 (SPACE_DIM = 3!) o
------------------o /| /| / | / | / | / | / | / | / | / | / | / | 4 /
| 5 / | o------------------o | | | | | | 3 o----------|-------o 2 | /
| / | / | / | / | / | / | / | / | / | / | / |/ |/ o------------------o
0 1

(PARENT) Linear 8-Node Hexahedron 3D Element Edge Node Map:

{ {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0,
4}, {1, 5}, {2, 6}, {3, 7} };

3D Element Face Node Map:

{ {0, 1, 5, 4}, {1, 2, 6, 5}, { 2, 3, 7, 6}, { 0, 4, 7, 3}, { 0, 3, 2,
1}, { 4, 5, 6, 7} }; Shards face list info:

typedef MakeTypeList< IndexList< 0, 1, 5, 4,   8, 13, 16, 12,   25 > ,
IndexList< 1, 2, 6, 5,   9, 14, 17, 13,   24 > , IndexList< 2, 3, 7,
6,  10, 15, 18, 14,   26 > , IndexList< 0, 4, 7, 3,  12, 19, 15, 11,
23 > , IndexList< 0, 3, 2, 1,  11, 10,  9,  8,   21 > , IndexList< 4,
5, 6, 7,  16, 17, 18, 19,   22 > >type HexahedronFaceNodeMap ;

After refinement (new nodes = *):

7 18 6 o---------*--------o /| /| /| / | / | / | 19 / | 22/ | / |
*---------*--------*17 | /| 15*----/|---*---/|---*14 / | /| / | /|26/
| /| | (PARENT) Linear 8-Node Hexahedron 4 / | / |16/ | / | / | / | |
3D Element Edge Node to mid-edge quadratic node map o---------*--20
----o5 |/ | | | 23*---|-|---*- 10|---*24 | | { 8, 9, 10, 11, | /|
3o-|--/|---*|--/|---o 2 | 16, 17, 18, 19, | / | / | / | / | / | / |
12, 13, 14, 15 } |/ | / 25|/ | / |/ | / | 12*---------*--------*13 |/
| Face to mid-face quadratic node map | 11*-----|---*----|---* 9 | {
25, 24, 26, 23, 21, 22 } | / | /21 | / | 0, 1, 2, 3, 4, 5 | / | / | /
|/ |/ |/ o---------*--------o 0 8 1

CHILD Linear 8-Node Hexahedron 3D Element Node Maps: | | static const
UInt child_0[] = { 0, 8, 21, 11, 12, 25, 20, 23 }; | static const UInt
child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20 }; | static const UInt
child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26 }; | static const UInt
child_3[] = { 11, 21, 10, 3, 23, 20, 26, 15 }; | static const UInt
child_4[] = { 12, 25, 20, 23, 4, 16, 22, 19 }; | static const UInt
child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22 }; | static const UInt
child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18 }; | static const UInt
child_7[] = { 23, 20, 26, 15, 19, 22, 18, 7 }; | PARENT Quadratic
20-Node Hexahedron Nodes 7 18 6 (SPACE_DIM = 3!) o--------o---------o
/| /| / | / | / | / | 19o | 17o | / 15o / o14 / | / | 4 / | 16 / | o
---------o--------o 5 | | | 10 | | | 3 o-------o--|-------o 2 | / | /
| / | / 12o / o13 / | o11 | o9 | / | / | / | / |/ |/ o---------o
--------o 0 8 1

PARENT Quadratic 20-Node Hexahedron 3D Element Edge Node Map:

| | static const UInt edge_0[] = { 0, 1, 8 }; | static const UInt
edge_1[] = { 1, 2, 9 }; | static const UInt edge_2[] = { 2, 3, 10 }; |
static const UInt edge_3[] = { 3, 0, 11 }; | static const UInt
edge_4[] = { 4, 5, 16 }; | static const UInt edge_5[] = { 5, 6, 17 };
| static const UInt edge_6[] = { 6, 7, 18 }; | static const UInt
edge_7[] = { 7, 4, 19 }; | static const UInt edge_8[] = { 0, 4, 12 };
| static const UInt edge_9[] = { 1, 5, 13 }; | static const UInt
edge_10[] = { 2, 6, 14 }; | static const UInt edge_11[] = { 3, 7, 15
}; |

CHILD Quadratic 20-Node Hexahedron 3D Element Node Maps: | | // Child
node tables use Two groups of nodes: | // a) vertices | // b) outer
edge mid-points | | static const UInt child_0[] = { 0, 8, 21, 11, 12,
25, 20, 23, | 27, 60, 67, 34, 35, 59, 79, 74, 51, 75, 77, 58 }; | |
static const UInt child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20, | 28, 29,
68, 60, 59, 36, 69, 79, 52, 53, 78, 75 }; | | static const UInt
child_2[] = { 21, 9, 2, 10, 20, 24, 14, 26, | 68, 30, 31, 61, 79, 69,
37, 62, 78, 54, 55, 76 }; | | static const UInt child_3[] = { 11, 21,
10, 3, 23, 20, 26, 15, | 67, 61, 32, 33, 74, 79, 62, 38, 77, 76, 56,
57 }; | | static const UInt child_4[] = { 12, 25, 20, 23, 4, 16, 22,
19, | 51, 75, 77, 58, 39, 66, 80, 73, 43, 65, 72, 50 }; | | static
const UInt child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22, | 52, 53, 78,
75, 66, 40, 70, 80, 44, 45, 71, 65 }; | | static const UInt child_6[]
= { 20, 24, 14, 26, 22, 17, 6, 18, | 78, 54, 55, 76, 80, 70, 41, 63,
71, 46, 47, 64 }; | | static const UInt child_7[] = { 23, 20, 26, 15,
19, 22, 18, 7, | 77, 76, 56, 57, 73, 80, 63, 42, 72, 64, 48, 49 }; |
PARENT Quadratic 27-Node Hexahedron Nodes 7 18 6 (SPACE_DIM = 3!) o
--------o---------o /| /| / | / | / | / | 19o | 17o | / 15o / o14 / |
/ | 4 / | 16 / | o---------o--------o 5 | | | 10 | | | 3 o-------
o--|-------o 2 | / | / | / | / 12o / o13 / | o11 | o9 | / | / | / | /
|/ |/ o---------o--------o 0 8 1

x--------x---------x /| /| / | / | / | 22 / | x | o x | / x o26 / x
(Node #20 is at centroid of element) / | / | / | / | \"2D surface\"
containing nodes 0, 8, 1, 13, 5, 16, 4, 12 has x---------x--------x |
node 25 at center.... | 23o | | o24 | | x-------x--|-------x | / | / |
/ 25 | / x / o x / | x o21 | x | / | / | / | / |/ |/ x---------x
--------x

PARENT Quadratic 27-Node Hexahedron 3D Element Edge Node Map:

| | static const UInt edge_0[] = { 0, 1, 8 }; | static const UInt
edge_1[] = { 1, 2, 9 }; | static const UInt edge_2[] = { 2, 3, 10 }; |
static const UInt edge_3[] = { 3, 0, 11 }; | static const UInt
edge_4[] = { 4, 5, 16 }; | static const UInt edge_5[] = { 5, 6, 17 };
| static const UInt edge_6[] = { 6, 7, 18 }; | static const UInt
edge_7[] = { 7, 4, 19 }; | static const UInt edge_8[] = { 0, 4, 12 };
| static const UInt edge_9[] = { 1, 5, 13 }; | static const UInt
edge_10[] = { 2, 6, 14 }; | static const UInt edge_11[] = { 3, 7, 15
}; |

Refined 27-Node Hexahedron Edge node tables: | | static const UInt
edge_0[] = { 0, 1, 8, 27, 28 }; | static const UInt edge_1[] = { 1, 2,
9, 29, 30 }; | static const UInt edge_2[] = { 2, 3, 10, 31, 32 }; |
static const UInt edge_3[] = { 3, 0, 11, 33, 34 }; | static const UInt
edge_4[] = { 4, 5, 16, 43, 44 }; | static const UInt edge_5[] = { 5,
6, 17, 45, 46 }; | static const UInt edge_6[] = { 6, 7, 18, 47, 48 };
| static const UInt edge_7[] = { 7, 4, 19, 49, 50 }; | static const
UInt edge_8[] = { 0, 4, 12, 35, 39 }; | static const UInt edge_9[] = {
1, 5, 13, 36, 40 }; | static const UInt edge_10[] = { 2, 6, 14, 37, 41
}; | static const UInt edge_11[] = { 3, 7, 15, 38, 42 }; |

CHILD 27-Node Hexahedron 3D Element Node Maps: | | // Child node
tables use Four groups of nodes: | // a) vertices | // b) outer edge
mid-points | // c) centroid | // d) mid-face points | | static const
UInt child_0[] = { 0, 8, 21, 11, 12, 25, 20, 23, | 27, 60, 67, 34, 35,
59, 79, 74, 51, 75, 77, 58, | 81, 89, 117, 97, 113, 105, 121 }; | |
static const UInt child_1[] = { 8, 1, 9, 21, 25, 13, 24, 20, | 28, 29,
68, 60, 59, 36, 69, 79, 52, 53, 78, 75, | 82, 92, 118, 113, 101, 106,
122 }; | | static const UInt child_2[] = { 21, 9, 2, 10, 20, 24, 14,
26, | 68, 30, 31, 61, 79, 69, 37, 62, 78, 54, 55, 76, | 83, 91, 119,
114, 102, 122, 109 }; | | static const UInt child_3[] = { 11, 21, 10,
3, 23, 20, 26, 15, | 67, 61, 32, 33, 74, 79, 62, 38, 77, 76, 56, 57, |
84, 90, 120, 100, 114, 121, 110 }; | | static const UInt child_4[] = {
12, 25, 20, 23, 4, 16, 22, 19, | 51, 75, 77, 58, 39, 66, 80, 73, 43,
65, 72, 50, | 85, 117, 93, 98, 116, 108, 124 }; | | static const UInt
child_5[] = { 25, 13, 24, 20, 16, 5, 17, 22, | 52, 53, 78, 75, 66, 40,
70, 80, 44, 45, 71, 65, | 86, 118, 94, 116, 104, 107, 123 }; | |
static const UInt child_6[] = { 20, 24, 14, 26, 22, 17, 6, 18, | 78,
54, 55, 76, 80, 70, 41, 63, 71, 46, 47, 64, | 87, 119, 95, 115, 103,
123, 112 }; | | static const UInt child_7[] = { 23, 20, 26, 15, 19,
22, 18, 7, | 77, 76, 56, 57, 73, 80, 63, 42, 72, 64, 48, 49, | 88,
120, 96, 99, 115, 124, 111 }; | | | Refined Hexagonal Element
\"Exterior\" Faces | (Local Face node numbering for
'Hierarchical/Consistent' Hex objects) | | 3 14 6 13 2 |
o----*----o----*----o | | | | | | 24 | 23 | | 15* * *19 * *12 | | | |
| | 8| 18 | | 7 o----*----o----*----o 5 | | 20 | | | | | | | 16* * 17*
* *11 | | 21 | 22 | | | | | | o----*----o----*----o | 0 9 4 10 1 | | |
Hexagonal object face topology child-nodes: | | // Face node tables
use Six groups of nodes: | // a) vertices (Local nodes: 0-1-2-3 ) | //
b) edge mid-points (Local nodes: 4-5-6-7 ) | // c) centroid (Local
node : 8 ) | // d) edge quater points (Local nodes:
9-10-11-12-13-14-15-16) | // e) interior edge mid-points (Local nodes:
17-18-19-20 ) | // f) mid-quadrant points (Local nodes: 21-22-23-24 )
| | static const UInt face_0[] = { 0, 1, 5, 4, 8, 13, 16, 12, 25, |
27, 28, 36, 40, 44, 43, 39, 35, | 59, 52, 66, 51, 105, 106, 107, 108
}; | | static const UInt face_1[] = { 1, 2, 6, 5, 9, 14, 17, 13, 24, |
29, 30, 37, 41, 46, 45, 40, 36, | 69, 54, 70, 53, 101, 102, 103, 104
}; | | static const UInt face_2[] = { 2, 3, 7, 6, 10, 15, 18, 14, 26,
| 31, 32, 38, 42, 48, 47, 41, 37, | 62, 56, 63, 55, 109, 110, 111, 112
}; | | static const UInt face_3[] = { 0, 4, 7, 3, 12, 19, 15, 11, 23,
| 35, 39, 50, 49, 42, 38, 33, 34, | 58, 73, 57, 74 97, 98, 99, 100 };
| | static const UInt face_4[] = { 0, 3, 2, 1, 11, 10, 9, 8 21, | 34,
33, 32, 31, 30, 29, 28, 27, | 67, 61, 68, 60, 89, 90, 91, 92 }; | |
static const UInt face_5[] = { 4, 5, 6, 7, 16, 17, 18, 19, 22, | 43,
44, 45, 46, 47, 48, 49, 50, | 65, 71, 64, 72, 93, 94, 95, 96 }; ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::wedge "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::wedge(UInt nnode)

Wedge Element with 6 Nodes Pentahedral Element with 6 Nodes 3-Sided
Prism Element with 6 Nodes

PARENT Linear 6-Node Wedge Nodes 5 (SPACE_DIM = 3!) . o . / \\\\ . /
\\\\ Face_Quad_4_3D() 0-1-4-3 . / \\\\ Face_Quad_4_3D() 1-4-5-2 . /
\\\\ Face_Quad_4_3D() 0-2-5-3 2 . o---------o 4 Face_Tri_3_3D() 0-1-2
o . 3 . Face_Tri_3_3D() 3-4-5 /. . /. \\\\ . /. \\\\ . /. \\\\ . o
---------o 0 1

CHILD Wedge6 node tables; | | static const UInt child_0[] = { 0, 6, 8,
9, 15, 17 }; | static const UInt child_1[] = { 6, 1, 7, 15, 10, 16 };
| static const UInt child_2[] = { 8, 7, 2, 17, 16, 11 }; | static
const UInt child_3[] = { 7, 8, 6, 16, 17, 15 }; | static const UInt
child_4[] = { 9, 15, 17, 3, 12, 14 }; | static const UInt child_5[] =
{ 15, 10, 16, 12, 4, 13 }; | static const UInt child_6[] = { 17, 16,
11, 14, 13, 5 }; | static const UInt child_7[] = { 16, 17, 15, 13, 14,
12 };

Refined Wedge6 Edge node tables: | | static const UInt edge_0[] = { 0,
1, 6 }; | static const UInt edge_1[] = { 1, 2, 7 }; | static const
UInt edge_2[] = { 2, 0, 8 }; | static const UInt edge_3[] = { 3, 4, 12
}; | static const UInt edge_4[] = { 4, 5, 13 }; | static const UInt
edge_5[] = { 5, 3, 14 }; | static const UInt edge_6[] = { 0, 3, 9 }; |
static const UInt edge_7[] = { 1, 4, 10 }; | static const UInt
edge_8[] = { 2, 5, 11 };

Refined Wedge6 Face node tables: | | static const UInt face_0[] = { 0,
1, 4, 3, 6, 10, 12, 9, 15 }; | static const UInt face_1[] = { 1, 2, 5,
4, 7, 11, 13, 10, 16 }; | static const UInt face_2[] = { 0, 3, 5, 2,
9, 14, 11, 8, 17 }; | static const UInt face_3[] = { 0, 2, 1, 8, 7, 6,
-1, -1, -1}; | static const UInt face_4[] = { 3, 4, 5, 12, 13, 14, -1,
-1, -1}; Wedge Element with 15 Nodes Pentahedral Element with 15 Nodes
3-Sided Prism Element with 15 Nodes

PARENT Quadratic 15-Node Wedge Nodes (SPACE_DIM = 3!)

5 o / \\\\ 14 / \\\\ o o 13 / \\\\ 3 / \\\\ o-----o-----o 12 4 11 o
Face_Quad_8_3D() 0-6-1-10-4-12-3-9 . . Face_Quad_8_3D()
1-10-4-13-5-11-2-7 . . Face_Quad_8_3D() 0-8-2-11-5-14-3-9 . .
Face_Tri_6_3D() 0-6-1-7-2-8 . . Face_Tri_6_3D() 3-12-4-13-5-14 . . 9
o...........o 10

2 o / \\\\ / \\\\ 8 o o 7 / \\\\ / \\\\ o-----o-----o 0 6 1

| | After Refinement: |

Face #0 Face #1 Face #2

3 30 12 31 4 4 32 13 33 5 5 34 14 35 3 o----*----o----*----o
o----*----o----*----o o----*----o----*----o | | | | | | | | | | | | |
| | | | | 51* *52 *53 53* *54 *55 55* *56 *51 | | | | | | | | | | 15|
25 | | 16| 27 | | 17| 29 | 9 o----*----o----*----o 10 10
o----*----o----*----o 11 11 o----*----o----*----o 9 | 24 | | | 26 | |
| 28 | | | | | | | | | | | 45* 46* *47 47* 48* 49* 49* 50* *45 ^ | | |
| | | | | | | | | | | | | | | | | o----*----o----*----o
o----*----o----*----o o----*----o----*----o # 0 18 6 19 1 1 20 7 21 2
2 22 8 23 0 # # \\\\ \\\\ \\\\ \\\\ --> -->

Face #4 Face #5

2 21 7 20 1 5 33 13 32 4 o-----*-----o-----*-----o
o-----*-----o-----*-----o \\\\ / \\\\ / \\\\ / \\\\ / \\\\ / \\\\ /
\\\\ / \\\\ / 22* 36* *38 *19 34* 42* *44 *31 \\\\ / \\\\ / \\\\ /
\\\\ / \\\\ / 37 \\\\ / \\\\ / 43 \\\\ / 8 o-----*-----o 6 14
o-----*-----o 12 \\\\ / \\\\ / \\\\ / \\\\ / ^ 23* *18 35* *30 ^ \\\\
\\\\ / \\\\ / / \\\\ \\\\ / \\\\ / / # o o # 0 3

CHILD Wedge15 node tables: | | static const UInt child_0[] = { 0, 6,
8, 9, 15, 17, 18, 37, 28, 45, 46, 50, 24, 40, 29 }; | static const
UInt child_1[] = { 6, 1, 7, 15, 10, 16, 19, 20, 38, 46, 47, 48, 25,
26, 41 }; | static const UInt child_2[] = { 8, 7, 2, 17, 16, 11, 36,
21, 22, 50, 48, 49, 39, 27, 28 }; | static const UInt child_3[] = { 9,
15, 17, 3, 12, 14, 24, 40, 29, 51, 52, 56, 30, 43, 35 }; | static
const UInt child_4[] = { 15, 10, 16, 12, 4, 13, 25, 26, 41, 52, 53,
54, 31, 32, 44 }; | static const UInt child_5[] = { 17, 16, 11, 14,
13, 5, 39, 27, 28, 56, 54, 55, 42, 33, 34 }; | static const UInt
child_6[] = { 7, 8, 6, 16, 17, 15, 36, 37, 38, 48, 50, 46, 39, 40, 41
}; | static const UInt child_7[] = { 16, 17, 15, 13, 14, 12, 39, 40,
41, 54, 56, 52, 42, 43, 44 }; |

Refined Wedge15 Edge node tables: | | static const UInt edge_0[] = {
0, 1, 6, 18, 19 }; | static const UInt edge_1[] = { 1, 2, 7, 20, 21 };
| static const UInt edge_2[] = { 2, 0, 8, 22, 23 }; | static const
UInt edge_3[] = { 3, 4, 12, 30, 31 }; | static const UInt edge_4[] = {
4, 5, 13, 32, 33 }; | static const UInt edge_5[] = { 5, 3, 14, 34, 35
}; | static const UInt edge_6[] = { 0, 3, 9, 45, 51 }; | static const
UInt edge_7[] = { 1, 4, 10, 47, 53 }; | static const UInt edge_8[] = {
2, 5, 11, 49, 55 };

Refined Wedge15 Face node tables: | | static const UInt face_0[] = {0,
1, 4, 3, 6, 10, 12, 9, 15, 18, 19, 47, 53, 31, 30, 51, 45, 46, 25, 52,
24}; | static const UInt face_1[] = {1, 2, 5, 4, 7, 11, 13, 10, 16,
20, 21, 49, 55, 33, 32, 53, 47, 48, 27, 54, 26}; | static const UInt
face_2[] = {0, 3, 5, 2, 9, 14, 11, 8, 17, 45, 51, 35, 34, 55, 49, 22,
23, 29, 56, 28, 50}; | static const UInt face_3[] = {0, 2, 1, 8, 7, 6,
23, 22, 21, 20, 19, 18, 37, 36, 38 }; | static const UInt face_4[] =
{3, 4, 5, 12, 13, 14, 30, 31, 32, 33, 34, 35, 43, 44, 42 }; | ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::pyramid
"const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::pyramid(UInt nnode) ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::Node_0
"const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Node_0() ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::Edge_2
"const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Edge_2() ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::Edge_3
"const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Edge_3() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Face_Tri_3 "const
MeshObjTopology * stk::adapt::Elem::StdMeshObjTopologies::Face_Tri_3()
";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Face_Tri_4 "const
MeshObjTopology * stk::adapt::Elem::StdMeshObjTopologies::Face_Tri_4()
";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Face_Tri_6 "const
MeshObjTopology * stk::adapt::Elem::StdMeshObjTopologies::Face_Tri_6()
";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Face_Quad_4 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Face_Quad_4() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Face_Quad_8 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Face_Quad_8() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Face_Quad_9 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Face_Quad_9() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Particle_1 "const
MeshObjTopology * stk::adapt::Elem::StdMeshObjTopologies::Particle_1()
";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::Rod_2 "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Rod_2() ";

%feature("docstring")  stk::adapt::Elem::StdMeshObjTopologies::Rod_3 "const MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Rod_3() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Shell_Line_2 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Shell_Line_2() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Shell_Line_3 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Shell_Line_3() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Shell_Tri_3 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Shell_Tri_3() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Shell_Tri_6 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Shell_Tri_6() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Shell_Quad_4 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Shell_Quad_4() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Shell_Quad_9 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Shell_Quad_9() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Tet_4 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Tet_4() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Tet_8 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Tet_8() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Tet_10 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Tet_10() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Wedge_6 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Wedge_6() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Wedge_15 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Wedge_15() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Hex_8 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Hex_8() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Hex_20 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Hex_20() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Hex_27 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Hex_27() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Pyramid_5 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Pyramid_5() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::Solid_Pyramid_13 "const
MeshObjTopology *
stk::adapt::Elem::StdMeshObjTopologies::Solid_Pyramid_13() ";

%feature("docstring")
stk::adapt::Elem::StdMeshObjTopologies::bootstrap "void
stk::adapt::Elem::StdMeshObjTopologies::bootstrap() ";


// File: namespacestk_1_1adapt_1_1regression__tests.xml
%feature("docstring")  stk::adapt::regression_tests::normalize "static void stk::adapt::regression_tests::normalize(double
input_normal[3], double normal[3]) ";

%feature("docstring")  stk::adapt::regression_tests::normalize "static void stk::adapt::regression_tests::normalize(double
input_output_normal[3]) ";

%feature("docstring")  stk::adapt::regression_tests::distance "static
double stk::adapt::regression_tests::distance(double c0[3], double
c1[3]) ";

%feature("docstring")  stk::adapt::regression_tests::difference "static void stk::adapt::regression_tests::difference(double v01[3],
double c0[3], double c1[3]) ";

%feature("docstring")  stk::adapt::regression_tests::dot "static
double stk::adapt::regression_tests::dot(double c0[3], double c1[3])
";

%feature("docstring")  stk::adapt::regression_tests::plane_dot_product
"static double stk::adapt::regression_tests::plane_dot_product(double
plane_point[3], double plane_normal[3], double point[3]) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_localRefiner,
break_tet_to_tet_N_5_ElementBased) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_localRefiner,
break_tet_to_tet_N_5_EdgeBased) ";

%feature("docstring")  stk::adapt::regression_tests::shock_function "static double stk::adapt::regression_tests::shock_function(double x)
";

%feature("docstring")  stk::adapt::regression_tests::shock_diff "static double
stk::adapt::regression_tests::shock_diff(stk::mesh::FieldBase
*nodal_refine_field, percept::PerceptMesh &eMesh, stk::mesh::Entity
&node0, stk::mesh::Entity &node1, double *coord0, double *coord1,
PlaneShock &shock, double shock_displacement) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_localRefiner,
break_tet_to_tet_N_5_EdgeBased_shock) ";

%feature("docstring")
stk::adapt::regression_tests::do_moving_shock_test "static void
stk::adapt::regression_tests::do_moving_shock_test(int num_time_steps,
bool save_intermediate=false, bool delete_parents=false) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_localRefiner,
break_tet_to_tet_N_5_EdgeBased_moving_shock) ";

%feature("docstring")
stk::adapt::regression_tests::do_moving_shock_test_cyl_sidesets "static void
stk::adapt::regression_tests::do_moving_shock_test_cyl_sidesets(int
num_time_steps, bool save_intermediate=false, bool
delete_parents=false) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_localRefiner,
break_tet_to_tet_N_5_EdgeBased_moving_shock_cyl_sidesets) ";

%feature("docstring")
stk::adapt::regression_tests::do_moving_shock_test_large_test "static
void stk::adapt::regression_tests::do_moving_shock_test_large_test(int
num_time_steps, bool save_intermediate=false, bool
delete_parents=false) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_localRefiner,
break_tet_to_tet_N_5_EdgeBased_moving_shock_large_test) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(nodeRegistry_regr,
test_parallel_0) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(nodeRegistry_regr,
test_parallel_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(nodeRegistry_regr,
test_parallel_2) ";

%feature("docstring")  stk::adapt::regression_tests::MegaByte "double
stk::adapt::regression_tests::MegaByte(MemorySizeType x) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(adapt, count_memory)
";

%feature("docstring")  stk::adapt::regression_tests::output_draw "static void stk::adapt::regression_tests::output_draw(std::string
filename, std::string toFile) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
pyramid_mesh) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
beam_enrich) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
beam_refine) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
generate_tables) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad4_to_quad8_to_quad8_shell) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_quad_shell) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tri_to_tri_shell) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
draw1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
draw) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_tri_6) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_tri_4) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_quad) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_quad_sierra)

uses the Sierra-ported tables from
framework/{element,mesh_modification} ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_quad_sierra_sidesets)

uses the Sierra-ported tables from
framework/{element,mesh_modification} ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_tet4_24_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_tet4_6_12_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_tet4_6_12_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
quad4_quad4_4_test_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_quad_sierra_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad_to_quad_sierra_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad4_to_quad9) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad4_to_quad8) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad8_to_quad8) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad4_to_quad9_to_quad9_0) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_quad4_to_quad9_to_quad9) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tri_to_tri_sierra_0) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tri_to_tri_sierra_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tri3_to_tri6_sierra) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tri3_to_tri6_to_tri6_sierra) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tet4_tet4_0) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tet4_tet4_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tet4_tet10_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tet4_tet10_tet10_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex8_8_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex8_8_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex27_1_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex27_1_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex20_1_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex20_1_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex20_hex20_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex20_hex20_1_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex27_hex27_0) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
hex8_hex27_hex27_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
wedge6_2) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
wedge6_enrich_1) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
wedge6_enrich_refine) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
heterogeneous_mesh) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
heterogeneous_mesh_sidesets) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
pyramid_mesh_enrich) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
biplane_refine) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_tet_shell3_tet) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
break_hex_shell4_hex) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
heterogeneous_mesh_enrich) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
heterogeneous_quadratic_refine) ";

%feature("docstring")  stk::adapt::regression_tests::STKUNIT_UNIT_TEST
"stk::adapt::regression_tests::STKUNIT_UNIT_TEST(regr_uniformRefiner,
wedge6_wedge18_enrich) ";


// File: namespacestk_1_1adapt_1_1unit__tests.xml
%feature("docstring")  stk::adapt::unit_tests::dw "static
stk::diag::Writer& stk::adapt::unit_tests::dw() ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(mesh_colorer, test1) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(mesh_colorer, test_quad) ";

%feature("docstring")  stk::adapt::unit_tests::save_or_diff "static
void stk::adapt::unit_tests::save_or_diff(PerceptMesh &eMesh,
std::string filename, int option=0)

This function either writes the given mesh to a file in Exodus format
(option 0) or, under option 1, checks if the file already exists, and
if so, treats that file as the \"gold\" copy and does a regression
difference check. ";

%feature("docstring")  stk::adapt::unit_tests::tet_volume "static
double stk::adapt::unit_tests::tet_volume(SingleTetFixture::Point
*node_coord_data, SingleTetFixture::TetIds &tetra_node_ids, unsigned
node_id_offset=0) ";

%feature("docstring")  stk::adapt::unit_tests::totalVolume "static
double stk::adapt::unit_tests::totalVolume(PerceptMesh &eMesh) ";

%feature("docstring")
stk::adapt::unit_tests::fixture_setup_NxNxN_box_hex_and_tet_mesh "static void
stk::adapt::unit_tests::fixture_setup_NxNxN_box_hex_and_tet_mesh() ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
triangulate_tet)

check triangulate_tet ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
triangulate_tet_2)

check triangulate_tet - two tets sharing a face ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
triangulate_tet_2_rand)

check triangulate_tet - two tets sharing a face, random coords ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
triangulate_tet_64)

check triangulate_tet - all 64 cases for a single tet ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
triangulate_tet_planes)

check triangulate_tet ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri)

Refine a triangle mesh.

Create a triangle mesh using the QuadFixture with the option of
breaking the quads into triangles Refine the triangle mesh, write the
results. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_1)

Refine a triangle mesh by trying to mark only one edge per triangle,
in a random-ish way. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_2)

Refine a triangle mesh by trying to mark only one edge per triangle,
in a random-ish way. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_1) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_2) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_3_1) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_3_2) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_3_1_IEdgeAdapter) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_3_1_IElementAdapter) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_5_ElementBased) ";

%feature("docstring")  stk::adapt::unit_tests::set_node_coords "static void
stk::adapt::unit_tests::set_node_coords(percept::PerceptMesh &eMesh,
mesh::PairIterRelation &elem_nodes, double tri_coords[3][3]) ";

%feature("docstring")  stk::adapt::unit_tests::convert_tuple "std::vector<int>
stk::adapt::unit_tests::convert_tuple(tri_tuple_type_local &tuple) ";

%feature("docstring")  stk::adapt::unit_tests::in_set "static bool
stk::adapt::unit_tests::in_set(tri_tuple_type_local &expected, vector<
tri_tuple_type_local > &base, bool reverse=false) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_MeshSizeRatio) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
break_tri_to_tri_N_EdgeBasedAnisotropic) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_localRefiner,
check_triangulate_face)

Create a single triangle mesh and mark the edges, call
RefinerPattern_Tri3_Tri3_N::triangulate_face and check properties of
the result - reverse the triangle polarity and check for consistency
";

%feature("docstring")  stk::adapt::unit_tests::s_diagWriter "static
stk::diag::Writer
stk::adapt::unit_tests::s_diagWriter(std::cout.rdbuf(), dw_enabled) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(nodeRegistry,
createAddNodes_serial_and_1st_parallel) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(nodeRegistry,
test_parallel_0) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(nodeRegistry,
test_parallel_1) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(nodeRegistry,
test_parallel_1_0)

set cell topology for the part block_hex_20 ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(SubDimCell, test1) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(SubDimCell, test2) ";

%feature("docstring")  stk::adapt::unit_tests::fixture_setup_0 "static void stk::adapt::unit_tests::fixture_setup_0()

Creates meshes for use in later tests 1. Create hex mesh from a
fixture and write it in Exodus format for use later. 2. Read the hex
mesh and convert it to tet elements using stk_adapt/UniformRefiner,
write it in Exodus format ";

%feature("docstring")  stk::adapt::unit_tests::fixture_setup_1 "static void stk::adapt::unit_tests::fixture_setup_1()

Creates meshes for use in later tests - quad meshes with and without
sidesets. ";

%feature("docstring")  stk::adapt::unit_tests::fixture_setup "static
void stk::adapt::unit_tests::fixture_setup() ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
stk_fixture_workaround) ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_quad_sierra_1_test)

Refine a quad mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tri_to_tri_sierra_1_test)

Refine a triangle mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner,
break_quad_to_quad_sierra)

Refine quad elements uses the Sierra-ported tables from
framework/{element,mesh_modification}

!eMesh, \"./square_quad4_ref_sierra_out.e\"); ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner,
break_quad_to_quad_sierra_1)

Refine quad elements with beam elements for the \"side sets\". ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner,
break_tri_to_tri_sierra)

Create a triangle mesh using the QuadFixture with the option of
breaking the quads into triangles Refine the triangle mesh, write the
results. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner,
hex8_hex8_8_1)

Refine a hex8 mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner,
wedge6_1)

Create and write a wedge mesh using the WedgeFixture. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
beam_enrich)

Create a Beam mesh and enrich it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
beam_refine)

Create a beam mesh and refine it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
generate_tables)

This code generates C++ tables used by stk_adapt - tables contain node
numbering, parametric coordinates consistent with Intrepid, and
related information needed by UniformRefiner. The generated code
should be compared with and merged into
<stk_adapt/sierra_element/GeneratedRefinementTable.hpp> as appropriate
if there is a change in the tables in that package, or an additional
element type is added. This comment is from the generated code and
tells how to bootstrap this process. Bootstrapping this file: to
create this file, run the regression test
RegressionTestUniformRefiner.cpp :: generate_tables after putting in a
dummy entry in ./sierra_element/GeneratedRefinementTable.hpp. The run
will produce a local file, generated_refinement_tables.hpp which can
be checked against the gold copy of GeneratedRefinementTable.hpp, then
copied over it. Add a call below to generate the actual new table
data. ";

%feature("docstring")  stk::adapt::unit_tests::output_draw "static
void stk::adapt::unit_tests::output_draw(std::string filename,
std::string toFile)

Code to generate Dot/Graphviz files representing the topology of
element refinement based on the internal tables. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner, draw1)
";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit_uniformRefiner, draw)
";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_tri_6)

Convert a quad mesh to triangles with 6 triangles per quad. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_tri_4)

Convert a quad mesh to triangles with 4 triangles per quad. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_quad)

Refine a quad mesh using the \"standalone\" refinement pattern:
UniformRefinerPattern<shards::Quadrilateral<4>,
shards::Quadrilateral<4>, 4 > This pattern is an example (like the
convert-type patterns) showing how to write a new pattern with no
dependencies ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_quad_sierra)

Refine a quad mesh uses the Sierra-ported tables from
framework/{element,mesh_modification} ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_quad_sierra_sidesets)

Refine a quad mesh with sidesets uses the Sierra-ported tables from
framework/{element,mesh_modification} ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_tet4_24_1)

Convert a hex mesh to tets using 24 tets per hex. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_tet4_6_12_1)

Convert a hex mesh using 6 tets per hex. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_tet4_6_12_2)

Convert a hex mesh using 6 tets per hex. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
quad4_quad4_4_test_1)

Refine a quad mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_quad_sierra_1)

Refine a quad mesh; test the multiple refinement feature. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad_to_quad_sierra_2)

Refine a quad mesh; test the multiple refinement feature. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad4_to_quad9)

Enrich a quad mesh (convert linear Quad4 elements to quadratic Quad9).
";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad4_to_quad8)

Enrich a quad mesh (convert linear Quad4 elements to serendepity
Quad8). ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad8_to_quad8)

Refine a quad8/serendepity mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad4_to_quad9_to_quad9_0)

Enrich a quad4 mesh to quad9 then refine it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_quad4_to_quad9_to_quad9)

Enrich a quad4 mesh to quad9 then refine it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tri_to_tri_sierra_0)

Refine a triangle mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tri_to_tri_sierra_1)

Refine a triangle mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tri3_to_tri6_sierra)

Enrich a triangle mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tri3_to_tri6_to_tri6_sierra)

Enrich a triangle mesh then refine it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tet4_tet4_0)

Refine a linear tet mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tet4_tet4_1)

Refine a linear tet mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tet4_tet10_1)

Enrich a linear tet mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
break_tet4_tet10_tet10_1)

Enrich a linear tet mesh then refine it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex8_8_1)

Refine a linear hex mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex8_8_2)

Refine a linear hex mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex27_1_1)

Enrich a linear hex mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex27_1_2)

Enrich a linear hex mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex20_1_1)

Enrich a linear hex mesh to serendepity hex20 elements. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex20_1_2)

Enrich a linear hex mesh to serendepity hex20 elements. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex20_hex20_1)

Refine a serendepity hex20 mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex20_hex20_1_2)

Refine a serendepity hex20 mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex27_hex27_0)

Refine a quadratic hex27 mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
hex8_hex27_hex27_1)

Refine a quadratic hex27 mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
wedge6_2)

Refine a linear wedge mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
wedge6_enrich_1)

Enrich a linear wedge mesh to serendepity Wedge15. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
wedge6_enrich_refine)

Enrich a linear wedge mesh to serendepity Wedge15 then refine it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
heterogeneous_mesh)

Generate a heterogeneous mesh (tet, hex, wedge elements) then refine
it. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
heterogeneous_mesh_enrich)

Enrich a heterogeneous mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
heterogeneous_quadratic_refine)

Refine the enriched heterogeneous mesh. ";

%feature("docstring")  stk::adapt::unit_tests::STKUNIT_UNIT_TEST "stk::adapt::unit_tests::STKUNIT_UNIT_TEST(unit1_uniformRefiner,
wedge6_wedge18_enrich)

Enrich a wedge6 mesh to wedge18. ";


// File: namespacestk_1_1diag.xml


// File: namespacestk_1_1mesh.xml


// File: namespacestk_1_1percept.xml
%feature("docstring")  stk::percept::interface_table::verifyMesh "bool stk::percept::verifyMesh(const BeamFixture &mesh) ";

%feature("docstring")  stk::percept::interface_table::verifyMesh "bool stk::percept::verifyMesh(const HeterogeneousFixture &mesh) ";

%feature("docstring")  stk::percept::interface_table::verifyMesh "bool stk::percept::verifyMesh(const PyramidFixture &mesh) ";

%feature("docstring")  stk::percept::interface_table::join "std::string stk::percept::join(std::string str1, std::string str2) ";

%feature("docstring")  stk::percept::interface_table::eval "double
stk::percept::eval(double x, double y, double z, double t,
Teuchos::RCP< Function > &func) ";

%feature("docstring")  stk::percept::interface_table::eval_print "void stk::percept::eval_print(double x, double y, double z, double t,
Teuchos::RCP< Function > &func) ";

%feature("docstring")  stk::percept::interface_table::eval_print2 "void stk::percept::eval_print2(double x, double y, double t,
Teuchos::RCP< Function > &func) ";

%feature("docstring")  stk::percept::interface_table::eval_vec3 "MDArray stk::percept::eval_vec3(double x, double y, double z, double
t, Teuchos::RCP< Function > &func) ";

%feature("docstring")  stk::percept::interface_table::eval_vec3_print
"void stk::percept::eval_vec3_print(double x, double y, double z,
double t, Teuchos::RCP< Function > &func) ";

%feature("docstring")  stk::percept::interface_table::eval "double
stk::percept::eval(double x, double y, double z, double t, Function
&func) ";

%feature("docstring")  stk::percept::interface_table::eval2 "double
stk::percept::eval2(double x, double y, double t, Teuchos::RCP<
Function > &func) ";

%feature("docstring")  stk::percept::interface_table::eval2 "double
stk::percept::eval2(double x, double y, double t, Function &func) ";

%feature("docstring")  stk::percept::interface_table::eval_print "void stk::percept::eval_print(double x, double y, double z, double t,
Function &func) ";

%feature("docstring")  stk::percept::interface_table::eval_print2 "void stk::percept::eval_print2(double x, double y, double t, Function
&func) ";

%feature("docstring")  stk::percept::interface_table::eval_vec3 "MDArray stk::percept::eval_vec3(double x, double y, double z, double
t, Function &func) ";

%feature("docstring")  stk::percept::interface_table::eval_vec3_print
"void stk::percept::eval_vec3_print(double x, double y, double z,
double t, Function &func) ";

%feature("docstring")  stk::percept::interface_table::first_dimensions
"static void stk::percept::first_dimensions(MDArray &arr, int
arr_offset, int *n_points, int max_rank=3) ";

%feature("docstring")
stk::percept::interface_table::get_entity_rank_names "static
std::vector<std::string> stk::percept::get_entity_rank_names(unsigned
dim) ";

%feature("docstring")  stk::percept::interface_table::push_back "void
stk::percept::push_back(std::vector< T > &dst, const std::vector< T >
&src) ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION
"stk::percept::IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION(DOFs_Tag)
IM_SHARDS_ARRAY_DIM_TAG_IMPLEMENTATION(BasisFields_Tag) void tni(void)
";

%feature("docstring")  stk::percept::interface_table::test "void
stk::percept::test() ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION "stk::percept::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION(Elements_Tag) ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION "stk::percept::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION(Cub_Points_Tag) ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION "stk::percept::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION(NodesPerElem_Tag) ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION "stk::percept::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION(Spatial_Dim_Tag) ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION "stk::percept::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION(DOFs_Tag) ";

%feature("docstring")
stk::percept::interface_table::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION "stk::percept::IM_SHARDS_ARRAY_DIM_TAG_DECLARATION(BasisFields_Tag) ";

%feature("docstring")  stk::percept::interface_table::tni "void
stk::percept::tni(void) ";

%feature("docstring")  stk::percept::interface_table::square "double
stk::percept::square(double x)

[DEPRECATED] ";

%feature("docstring")  stk::percept::interface_table::noendl "std::ostream& stk::percept::noendl(std::ostream &os) ";

%feature("docstring")  stk::percept::interface_table::checkOmit "static void stk::percept::checkOmit(const std::vector< T * >
&collection, std::string omit_part) ";

%feature("docstring")
stk::percept::interface_table::checkForPartsToAvoidReading "static
void stk::percept::checkForPartsToAvoidReading(Ioss::Region
&in_region, std::string omit_part) ";

%feature("docstring")
stk::percept::interface_table::setup_spatialDim_metaData "static void
stk::percept::setup_spatialDim_metaData(Ioss::Region &region,
stk::mesh::fem::FEMMetaData &meta, int &spatial_dim) ";

%feature("docstring")
stk::percept::interface_table::decipher_filename "static bool
stk::percept::decipher_filename(std::string filename_in, int
&my_processor, int &processor_count) ";

%feature("docstring")
stk::percept::interface_table::percept_create_output_mesh "static
void stk::percept::percept_create_output_mesh(const std::string
&filename, stk::ParallelMachine comm, stk::mesh::BulkData &bulk_data,
stk::io::MeshData &mesh_data) ";

%feature("docstring")  stk::percept::interface_table::out "std::ostream & stk::percept::out()

Normal output stream. ";

%feature("docstring")  stk::percept::interface_table::pout "std::ostream & stk::percept::pout()

Per-processor output stream (See RuntimeDeferredx). ";

%feature("docstring")  stk::percept::interface_table::dout "std::ostream & stk::percept::dout()

Diagnostic output stream. ";

%feature("docstring")  stk::percept::interface_table::tout "std::ostream & stk::percept::tout()

Regression test textual output stream. ";

%feature("docstring")  stk::percept::interface_table::dwout "std::ostream & stk::percept::dwout()

Diagnostic writer stream. ";

%feature("docstring")  stk::percept::interface_table::dw "stk::diag::Writer & stk::percept::dw() ";

%feature("docstring")  stk::percept::interface_table::report_handler "void stk::percept::report_handler(const char *message, int type) ";

%feature("docstring")  stk::percept::interface_table::timerSet "stk::diag::TimerSet & stk::percept::timerSet() ";

%feature("docstring")  stk::percept::interface_table::timer "stk::diag::Timer & stk::percept::timer() ";

%feature("docstring")  stk::percept::interface_table::panic "static
void stk::percept::panic() ";

%feature("docstring")  stk::percept::interface_table::runCommand "static void stk::percept::runCommand(std::string command) ";

%feature("docstring")
stk::percept::interface_table::my_report_handler "void
stk::percept::my_report_handler(const char *message, int type) ";

%feature("docstring")  stk::percept::interface_table::get_memory_info
"static void stk::percept::get_memory_info(size_t &memory_usage,
size_t &faults) ";

%feature("docstring")  stk::percept::interface_table::get_heap_info "static void stk::percept::get_heap_info(size_t &heap_size, size_t
&largest_free) ";

%feature("docstring")  stk::percept::interface_table::perceptTimerSet
"stk::diag::TimerSet & stk::percept::perceptTimerSet() ";

%feature("docstring")  stk::percept::interface_table::perceptTimer "stk::diag::Timer & stk::percept::perceptTimer() ";

%feature("docstring")  stk::percept::interface_table::getLapTime "LapTimeType stk::percept::getLapTime(stk::diag::Timer &lap_timer) ";

%feature("docstring")
stk::percept::interface_table::getAccumulatedLap "LapCountType
stk::percept::getAccumulatedLap(stk::diag::Timer &timer, bool option)
";

%feature("docstring")  stk::percept::interface_table::toString "std::string stk::percept::toString(T t) ";

%feature("docstring")  stk::percept::interface_table::square "T
stk::percept::square(T t) ";

%feature("docstring")  stk::percept::interface_table::SQR "T
stk::percept::SQR(T t) ";

%feature("docstring")  stk::percept::interface_table::toInt "int
stk::percept::toInt(std::string t) ";


// File: namespacestk_1_1percept_1_1@198.xml


// File: namespacestk_1_1percept_1_1@202.xml


// File: namespacestk_1_1percept_1_1@204.xml


// File: namespacestk_1_1percept_1_1@207.xml


// File: namespacestk_1_1percept_1_1@304.xml


// File: namespacestk_1_1percept_1_1@305.xml


// File: namespacestk_1_1percept_1_1interface__table.xml


// File: namespacestk_1_1percept_1_1io__util.xml
%feature("docstring")
stk::percept::io_util::process_read_nodeblocks_meta "void
stk::percept::io_util::process_read_nodeblocks_meta(Ioss::Region
&region, stk::mesh::fem::FEMMetaData &meta, int &spatial_dim)

Declare \"coordinates\" field and put it on the universal part. This
example also defines all Ioss::Field::TRANSIENT fields that exist on
the Ioss::Nodeblock as fields on the universal part. ";

%feature("docstring")
stk::percept::io_util::process_read_elementblocks_meta "void
stk::percept::io_util::process_read_elementblocks_meta(Ioss::Region
&region, stk::mesh::fem::FEMMetaData &meta)

Declare a part for each element block on the Ioss::Region 'region'
unless the element block has the \"omitted\" property set to the value
1. The example then iterates each element block and defines any
Ioss::Field::ATTRIBUTE and Ioss::Field::TRANSIENT fields that exist on
the Ioss::ElementBlock as fields on the corresponding part. ";

%feature("docstring")
stk::percept::io_util::process_read_nodesets_meta "void
stk::percept::io_util::process_read_nodesets_meta(Ioss::Region
&region, stk::mesh::fem::FEMMetaData &meta)

Declare a part for each Ioss::NodeSet on the Ioss::Region 'region'
unless the nodeset has the \"omitted\" property set to the value 1.
The example then iterates each nodeset and defines any \"distribution
factor\" and Ioss::Field::TRANSIENT fields that exist on the
Ioss::NodeSet as fields on the corresponding part. ";

%feature("docstring")
stk::percept::io_util::process_read_sidesets_meta "void
stk::percept::io_util::process_read_sidesets_meta(Ioss::Region
&region, stk::mesh::fem::FEMMetaData &meta)

Declare a part for each Ioss::SideSet on the Ioss::Region 'region'
unless the sideset has the \"omitted\" property set to the value 1.
The example then iterates each sideset and defines any \"distribution
factor\" and Ioss::Field::TRANSIENT fields that exist on the
Ioss::SideSet as fields on the corresponding part.

Each sideblock in the active sidesets is then processed by defining a
part for each Ioss::SideBlock on the Ioss::SideSet unless the
sideblock has the \"omitted\" property set to the value 1. The example
then iterates each sideblock and defines any \"distribution factor\"
and Ioss::Field::TRANSIENT fields that exist on the Ioss::SideBlock as
fields on the corresponding part. ";

%feature("docstring")
stk::percept::io_util::process_read_nodeblocks_bulk "void
stk::percept::io_util::process_read_nodeblocks_bulk(Ioss::Region
&region, stk::mesh::BulkData &bulk)

NOTE: This must be called after the process_read_elementblocks() call
since there may be nodes that exist in the database that are not part
of the analysis mesh due to subsetting of the element blocks.

Populates the \"coordinates\" field for all active nodes in the model.
";

%feature("docstring")
stk::percept::io_util::process_read_elementblocks_bulk "void
stk::percept::io_util::process_read_elementblocks_bulk(Ioss::Region
&region, stk::mesh::BulkData &bulk)

NOTE: This should be the first function called of any of the
\"process_read_X\" type functions that take an stk::mesh::BulkData
argument, especially if the input Ioss::Region mesh is going to be
subsetted (have element blocks omitted).

This function iterates all non-omitted element blocks and declares
each element (and the corresponding nodes) in the element block. If
there are any Ioss::Field::ATTRIBUTE fields on the element block (for
example, shell thickness or particle radius), then that field data is
alse read and the corresponding stk::mesh::Field populated. ";

%feature("docstring")
stk::percept::io_util::process_read_nodesets_bulk "void
stk::percept::io_util::process_read_nodesets_bulk(Ioss::Region
&region, stk::mesh::BulkData &bulk)

Iterates each non-omitted Ioss::NodeSet and then iterates each node in
the Ioss::NodeSet. If the node exists (that is, it is connected to a
non-omitted Ioss::ElementBlock), then that node is associated with the
part corresponding to this Ioss::NodeSet. If the
\"distribution_factor\" field exists, then that data is also
associated with the field. ";

%feature("docstring")
stk::percept::io_util::process_read_sidesets_bulk "void
stk::percept::io_util::process_read_sidesets_bulk(Ioss::Region
&region, stk::mesh::BulkData &bulk)

Process each non-omitted Ioss::SideSet and the contained non-omitted
Ioss::SideBlock and associate each element-side pair with the
corresponding part if the underlying element is active. If the
\"distribution_factor\" field exists, then that data is also
associated with the corresponding field. ";

%feature("docstring")
stk::percept::io_util::process_read_input_request "void
stk::percept::io_util::process_read_input_request(Ioss::Region
&region, stk::mesh::BulkData &bulk, int step)

A minimal example function showing how field data on the Ioss::Region
entities can be periodically transferred to the corresponding field(s)
on the stk::mesh entities. This would be used to bring in initial
condition data or interpolation data or any other scenario in which
data on the mesh file needs to be transferred to the stk::mesh fields.
";

%feature("docstring")  stk::percept::io_util::process_output_request "void stk::percept::io_util::process_output_request(Ioss::Region
&region, stk::mesh::BulkData &bulk, int step)

A minimal example function showing how stk::mesh field data can
periodically be output to a results, history, heartbeat, or restart
database. The scheduling would be done either in this function or at a
higher level and is not shown here. The function iterates all parts
and if there is a corresponding Ioss part on the Ioss::Region, all
fields defined to be output are iterated and their data output to the
corresponding Ioss::Field. The function calls the
stk::io::is_valid_part_field() function to determine whether the field
should be output and then calls the stk::io::field_data_to_ioss()
function to do the actual output of the field. ";


// File: namespacestk_1_1percept_1_1regression__tests.xml
%feature("docstring")
stk::percept::regression_tests::STKUNIT_UNIT_TEST "stk::percept::regression_tests::STKUNIT_UNIT_TEST(perceptMesh,
open_new_close_PerceptMesh)

eval_print(x,y,z,time, sf_error); ";

%feature("docstring")
stk::percept::regression_tests::STKUNIT_UNIT_TEST "stk::percept::regression_tests::STKUNIT_UNIT_TEST(perceptMesh,
open_new_close_PerceptMesh_2) ";

%feature("docstring")
stk::percept::regression_tests::STKUNIT_UNIT_TEST "stk::percept::regression_tests::STKUNIT_UNIT_TEST(perceptMesh,
open_new_close_PerceptMesh_3) ";

%feature("docstring")
stk::percept::regression_tests::STKUNIT_UNIT_TEST "stk::percept::regression_tests::STKUNIT_UNIT_TEST(perceptMesh,
open_new_reopen_PerceptMesh)

reopen the mesh to allow for more fields to be added - note that this
involves a db write/read operation ";


// File: namespacestk_1_1percept_1_1unit__tests.xml
%feature("docstring")  stk::percept::unit_tests::test_shards_array "void stk::percept::unit_tests::test_shards_array() ";

%feature("docstring")  stk::percept::unit_tests::testSweepMesher "int
stk::percept::unit_tests::testSweepMesher(stk::ParallelMachine
parallel_machine) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_demo_1_0_0) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_read_print) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_demo_1) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_demo_2) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_readMesh_createField_interpolateFrom)

the coordinates field is always created by the PerceptMesh read
operation, here we just get the field

get the new field created by readModelCreateOptionalFields()

create a field function from the existing coordinates field

here we could evaluate this field function

create a field function to represent the new coordinate magnitude
field, and interpolate the string function to its nodes

check that the coordinates mag field is set correctly ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_multiplePoints)

the coordinates field is always created by the PerceptMesh read
operation, here we just get the field ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_point_eval_verify)

test evaluation of field function at a point

the coordinates field is always created by the PerceptMesh read
operation, here we just get the field

create a field function from the existing coordinates field

here we evaluate this field function ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_point_eval_deriv_verify)

test evaluation of field function at a point

the coordinates field is always created by the PerceptMesh read
operation, here we just get the field

create a field function from the existing coordinates field

here we evaluate this field function ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
fieldFunction_point_eval_timing)

test evaluation of field function at a point

the coordinates field is always created by the PerceptMesh read
operation, here we just get the field

create a field function from the existing coordinates field

! pts(0) = 0.2; pts(1) = 0.3; pts(2)= 0.4; ";

%feature("docstring")  stk::percept::unit_tests::s_diagWriter "static
stk::diag::Writer
stk::percept::unit_tests::s_diagWriter(std::cout.rdbuf(), dw_enabled)
";

%feature("docstring")  stk::percept::unit_tests::dw "static
stk::diag::Writer& stk::percept::unit_tests::dw() ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(geom, volume) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, volume)

create a field function from the existing coordinates field

the function to be integrated - here it is just the identity, and when
integrated should produce the volume

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work

get the l2 norm of identity

Note: need to create new fields each time, which requires a change to
the meta data ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, surface_area)

create a field function from the existing coordinates field

the function to be integrated - here it is just the identity, and when
integrated should produce the area of faces

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work

get the l2 norm of identity ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, string_function)

Create the operator that will do the work get the l2 norm

the function to be integrated: (Integral[ abs(x), dxdydz]) =?= (2 *
|x|^2/2 @ [0, 0.5]) ==> .25)

the function to be integrated: (Max[ x^2+y^3+z^4, dxdydz]) =?= (@
[-0.5, 0.5]^3 ) ==> .5^2+.5^3+.5^4)

indirection

the function to be integrated: sqrt(Integral[(x*y*z)^2, dxdydz]) =?=
(see unitTest1.py)

indirection

the function to be integrated (but over a rotated domain):
sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py) now rotate
the mesh ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, string_function_1)

Create the operator that will do the work get the l2 norm

the function to be integrated: sqrt(Integral[(x*y*z)^2, dxdydz]) =?=
(see unitTest1.py)

indirection ";

%feature("docstring")
stk::percept::unit_tests::TEST_norm_string_function_turbo_verify_correctness
"void
stk::percept::unit_tests::TEST_norm_string_function_turbo_verify_correctness(TurboOption
turboOpt)

This test uses a back door to the function that passes in the element
to avoid the lookup of the element when the StringFunction contains
references to FieldFunctions

create a field function from the existing coordinates field

the function to be integrated: sqrt(Integral[x^2, dxdydz]) =?=
sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work get the l2 norm

the function to be integrated: (Integral[ abs(x), dxdydz]) =?= (2 *
|x|^2/2 @ [0, 0.5]) ==> .25)

----- here the function to be integrated: sqrt(Integral[(x*y*z)^2,
dxdydz]) =?= (see unitTest1.py)

the function to be integrated (but over a rotated domain):
sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py) now rotate
the mesh ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm,
string_function_turbo_verify_correctness_element) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm,
string_function_turbo_verify_correctness_bucket) ";

%feature("docstring")
stk::percept::unit_tests::TEST_norm_string_function_turbo_timings "void
stk::percept::unit_tests::TEST_norm_string_function_turbo_timings(TurboOption
turboOpt)

This test uses a back door to the function that passes in the element
to avoid the lookup of the element when the StringFunction contains
references to FieldFunctions

create a meta data/bulk data empty pair

the coordinates field is always created by the PerceptMesh read
operation, here we just get the field

create a field function from the existing coordinates field

the function to be integrated: sqrt(Integral[x^2, dxdydz]) =?=
sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)

the function to be integrated: sqrt(Integral[x^2, dxdydz]) =?=
sqrt(x^3/3 @ [-0.5, 0.5]) ==> sqrt(0.25/3)

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work get the l2 norm

the function to be integrated: (Integral[ abs(x), dxdydz]) =?= (2 *
|x|^2/2 @ [0, 0.5]) ==> .25)

the function to be integrated: sqrt(Integral[(x*y*z)^2, dxdydz]) =?=
(see unitTest1.py)

the function to be integrated (but over a rotated domain):
sqrt(Integral[(x*y*z)^2, dxdydz]) =?= (see unitTest2.py) now rotate
the mesh ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm,
string_function_turbo_timings) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm,
string_function_turbo_timings_bucket) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, field_function)

Create the operator that will do the work get the l2 norm

the function to be integrated: (Integral[ abs(x), dxdydz]) =?= (2 *
|x|^2/2 @ [0, 0.5]) ==> .25) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, h1_volume)

create a field function from the existing coordinates field

the function to be integrated - here it is just the identity, and when
integrated should produce the volume

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work

get the l2 norm of identity ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, h1_volume_1)

create a field function from the existing coordinates field

the function to be integrated - here it is just the identity, and when
integrated should produce the volume

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work

get the l2 norm of plane ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(norm, h1_volume_2)

create a field function from the existing coordinates field

the function to be integrated - here it is just the identity, and when
integrated should produce the volume

A place to hold the result. This is a \"writable\" function (we may
want to make this explicit - StringFunctions are not writable;
FieldFunctions are since we interpolate values to them from other
functions).

Create the operator that will do the work

get the l2 norm of plane ";

%feature("docstring")  stk::percept::unit_tests::fixture_setup_0 "static void stk::percept::unit_tests::fixture_setup_0()

create a mesh of hex elements for use in other tests below ";

%feature("docstring")  stk::percept::unit_tests::fixture_setup_1 "static void stk::percept::unit_tests::fixture_setup_1()

using the QuadFixture, generate meshes with and without sidesets ";

%feature("docstring")  stk::percept::unit_tests::fixture_setup "static void stk::percept::unit_tests::fixture_setup() ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(unit_perceptMesh,
wedge6_1)

generate a mesh with wedge elements ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(perceptMesh, walk_nodes)

a tutorial on some of the innards of stk_mesh database ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(perceptMesh,
test_mesh_diff)

Test the mesh_difference capability of PerceptMesh and the interface
to stk_io 1. read (and write and read back in) meshes generated above
(quad_fixture) 2. invoke PerceptMesh::print_info(ostringstream...) to
create a string representation of the mesh 3. compare the string with
the saved, gold value of the string 4. invoke mesh_difference to
ensure it behaves as expected (two meshes are shown as identical) 5.
modify one mesh and ensure mesh_difference shows the meshes as being
different ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(perceptMesh,
create_skewed_mesh) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(perceptMesh,
create_quad_streaming_mesh) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(perceptMesh, test_states)
";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(search, test1)

dw().m(LOG_SEARCH) << \"Use case 1\" << stk::diag::push <<
stk::diag::dendl; ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_xy_basic) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_xy_basic_1) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_xy_basic_2) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_test_alias) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_vector_valued) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_arithmetic_ops) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_derivative) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_derivative_1) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_derivative_2) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_multiplePoints)

indirection ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_expressions) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(function,
stringFunction_timing) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(unit_tests_percept,
noMallocArray) ";

%feature("docstring")  stk::percept::unit_tests::setupMap "static
void stk::percept::unit_tests::setupMap(MAP &map, unsigned N) ";

%feature("docstring")  stk::percept::unit_tests::find1 "static
unsigned* stk::percept::unit_tests::find1(MAP &map, ITER &i, unsigned
key) ";

%feature("docstring")  stk::percept::unit_tests::find2 "static
unsigned* stk::percept::unit_tests::find2(MAP &map, ITER &i, unsigned
key) ";

%feature("docstring")  stk::percept::unit_tests::dot1 "static double
stk::percept::unit_tests::dot1(MAP &map, ITER &it, unsigned N,
unsigned niter, FUNC &fm) ";

%feature("docstring")  stk::percept::unit_tests::doTest "static void
stk::percept::unit_tests::doTest(MAP &map, ITER &it, unsigned N,
unsigned niter, FUNC &fm, std::string msg) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(time_maps,
compare_different_maps) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(topo, testCrossedElems) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(geom, geomPrints) ";

%feature("docstring")  stk::percept::unit_tests::STKUNIT_UNIT_TEST "stk::percept::unit_tests::STKUNIT_UNIT_TEST(geom, geomEqui) ";


// File: namespacestk_1_1percept_1_1util.xml


// File: namespacestk_1_1utils.xml


// File: namespacestk__example__io.xml
%feature("docstring")  stk_example_io::process_nodeblocks "void
stk_example_io::process_nodeblocks(Ioss::Region &region,
stk::mesh::fem::FEMMetaData &meta)

Declare \"coordinates\" field and put it on the universal part. This
example also defines all Ioss::Field::TRANSIENT fields that exist on
the Ioss::Nodeblock as fields on the universal part. ";

%feature("docstring")  stk_example_io::process_elementblocks "void
stk_example_io::process_elementblocks(Ioss::Region &region,
stk::mesh::fem::FEMMetaData &meta)

Declare a part for each element block on the Ioss::Region 'region'
unless the element block has the \"omitted\" property set to the value
1. The example then iterates each element block and defines any
Ioss::Field::ATTRIBUTE and Ioss::Field::TRANSIENT fields that exist on
the Ioss::ElementBlock as fields on the corresponding part. ";

%feature("docstring")  stk_example_io::process_nodesets "void
stk_example_io::process_nodesets(Ioss::Region &region,
stk::mesh::fem::FEMMetaData &meta)

Declare a part for each Ioss::NodeSet on the Ioss::Region 'region'
unless the nodeset has the \"omitted\" property set to the value 1.
The example then iterates each nodeset and defines any \"distribution
factor\" and Ioss::Field::TRANSIENT fields that exist on the
Ioss::NodeSet as fields on the corresponding part. ";

%feature("docstring")  stk_example_io::process_sidesets "void
stk_example_io::process_sidesets(Ioss::Region &region,
stk::mesh::fem::FEMMetaData &meta)

Declare a part for each Ioss::SideSet on the Ioss::Region 'region'
unless the sideset has the \"omitted\" property set to the value 1.
The example then iterates each sideset and defines any \"distribution
factor\" and Ioss::Field::TRANSIENT fields that exist on the
Ioss::SideSet as fields on the corresponding part.

Each sideblock in the active sidesets is then processed by defining a
part for each Ioss::SideBlock on the Ioss::SideSet unless the
sideblock has the \"omitted\" property set to the value 1. The example
then iterates each sideblock and defines any \"distribution factor\"
and Ioss::Field::TRANSIENT fields that exist on the Ioss::SideBlock as
fields on the corresponding part.

<

< ";

%feature("docstring")  stk_example_io::process_nodeblocks "void
stk_example_io::process_nodeblocks(Ioss::Region &region,
stk::mesh::BulkData &bulk)

NOTE: This must be called after the process_elementblocks() call since
there may be nodes that exist in the database that are not part of the
analysis mesh due to subsetting of the element blocks.

Populates the \"coordinates\" field for all active nodes in the model.
";

%feature("docstring")  stk_example_io::process_elementblocks "void
stk_example_io::process_elementblocks(Ioss::Region &region,
stk::mesh::BulkData &bulk)

NOTE: This should be the first function called of any of the
\"process_X\" type functions that take an stk::mesh::BulkData
argument, especially if the input Ioss::Region mesh is going to be
subsetted (have element blocks omitted).

This function iterates all non-omitted element blocks and declares
each element (and the corresponding nodes) in the element block. If
there are any Ioss::Field::ATTRIBUTE fields on the element block (for
example, shell thickness or particle radius), then that field data is
alse read and the corresponding stk::mesh::Field populated. ";

%feature("docstring")  stk_example_io::process_nodesets "void
stk_example_io::process_nodesets(Ioss::Region &region,
stk::mesh::BulkData &bulk)

Iterates each non-omitted Ioss::NodeSet and then iterates each node in
the Ioss::NodeSet. If the node exists (that is, it is connected to a
non-omitted Ioss::ElementBlock), then that node is associated with the
part corresponding to this Ioss::NodeSet. If the
\"distribution_factor\" field exists, then that data is also
associated with the field. ";

%feature("docstring")  stk_example_io::process_sidesets "void
stk_example_io::process_sidesets(Ioss::Region &region,
stk::mesh::BulkData &bulk)

Process each non-omitted Ioss::SideSet and the contained non-omitted
Ioss::SideBlock and associate each element-side pair with the
corresponding part if the underlying element is active. If the
\"distribution_factor\" field exists, then that data is also
associated with the corresponding field. ";

%feature("docstring")  stk_example_io::process_input_request "void
stk_example_io::process_input_request(Ioss::Region &region,
stk::mesh::BulkData &bulk, int step)

A minimal example function showing how field data on the Ioss::Region
entities can be periodically transferred to the corresponding field(s)
on the stk::mesh entities. This would be used to bring in initial
condition data or interpolation data or any other scenario in which
data on the mesh file needs to be transferred to the stk::mesh fields.
";

%feature("docstring")  stk_example_io::process_output_request "void
stk_example_io::process_output_request(Ioss::Region &region,
stk::mesh::BulkData &bulk, int step)

A minimal example function showing how stk::mesh field data can
periodically be output to a results, history, heartbeat, or restart
database. The scheduling would be done either in this function or at a
higher level and is not shown here. The function iterates all parts
and if there is a corresponding Ioss part on the Ioss::Region, all
fields defined to be output are iterated and their data output to the
corresponding Ioss::Field. The function calls the
stk::io::is_valid_part_field() function to determine whether the field
should be output and then calls the stk::io::field_data_to_ioss()
function to do the actual output of the field. ";

%feature("docstring")  stk_example_io::my_test "void
stk_example_io::my_test(mesh::BulkData &M, const unsigned elem_type,
const VectorFieldType &coord_field, const VectorFieldType
&elem_centroid_field) ";

%feature("docstring")  stk_example_io::io_example "void
stk_example_io::io_example(stk::ParallelMachine comm, const
std::string &in_filename, const std::string &out_filename)

The real app would also only register a subset of the stk::mesh fields
as output fields and would probably have a mapping from the internally
used name to some name picked by the user. In this example, all
Ioss::Field::TRANSIENT fields defined on the stk::mesh are output to
the results database and the internal stk::mesh field name is used as
the name on the database.... ";

%feature("docstring")  stk_example_io::process_surface_entity "void
stk_example_io::process_surface_entity(Ioss::SideSet *entity,
stk::mesh::fem::FEMMetaData &meta, stk::mesh::EntityRank entity_rank)
";

%feature("docstring")  stk_example_io::process_surface_entity "void
stk_example_io::process_surface_entity(const Ioss::SideSet *io,
stk::mesh::BulkData &bulk) ";

%feature("docstring")  stk_example_io::get_field_data "void
stk_example_io::get_field_data(stk::mesh::BulkData &bulk,
stk::mesh::Part &part, stk::mesh::EntityRank part_type,
Ioss::GroupingEntity *io_entity, Ioss::Field::RoleType filter_role) ";

%feature("docstring")  stk_example_io::put_field_data "void
stk_example_io::put_field_data(stk::mesh::BulkData &bulk,
stk::mesh::Part &part, stk::mesh::EntityRank part_type,
Ioss::GroupingEntity *io_entity, Ioss::Field::RoleType filter_role) ";


// File: namespacestk__percept__unit.xml
%feature("docstring")  stk_percept_unit::use_encr_case_1_driver "void
stk_percept_unit::use_encr_case_1_driver(MPI_Comm comm)

set cell topology for the part block_1

set cell topology for the part block_1 ";

%feature("docstring")  stk_percept_unit::myMain "int
stk_percept_unit::myMain() ";

%feature("docstring")  stk_percept_unit::STKUNIT_UNIT_TEST "stk_percept_unit::STKUNIT_UNIT_TEST(topo, test1) ";

%feature("docstring")  stk_percept_unit::use_encr_case_1_generate_mesh
"void stk_percept_unit::use_encr_case_1_generate_mesh(mesh::BulkData
&mesh, const unsigned N[], const VectorFieldType &node_coord, const
ElementNodePointerFieldType &elem_node_coord, mesh::Part &hex_block)
";


// File: AdaptMain_8cpp.xml
%feature("docstring")  stk::main "int main(int argc, char **argv) ";


// File: Allocate_8cpp.xml


// File: Basis_81_8hpp.xml


// File: Basis_82_8hpp.xml


// File: Basis_8hpp.xml


// File: BeamFixture_8cpp.xml


// File: BeamFixture_8hpp.xml


// File: BucketOp_8hpp.xml


// File: BuildBoundingBoxes_8hpp.xml


// File: BuildBoundingBoxesDef_8hpp.xml


// File: CellTopology_8cpp.xml


// File: CellTopology_8hpp.xml


// File: build_81_8nogit_2CMakeFiles_2CompilerIdCXX_2CMakeCXXCompilerId_8cpp.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";


// File: build_8dir_2CMakeFiles_2CompilerIdCXX_2CMakeCXXCompilerId_8cpp.xml
%feature("docstring")  main "int main(int argc, char *argv[]) ";


// File: Colorer_8cpp.xml


// File: Colorer_8hpp.xml


// File: CompositeFunction_8hpp.xml


// File: ComputeBases_8hpp.xml


// File: ComputeFieldValues_8hpp.xml


// File: ConstantFunction_8hpp.xml


// File: cpp-templates_81_8cpp.xml
%feature("docstring")  main "int main() ";


// File: cpp-templates_82_8cpp.xml
%feature("docstring")  main "int main() ";


// File: cpp-templates_83_8cpp.xml
%feature("docstring")  init "void init(T t, I i) ";

%feature("docstring")  init "void init(T t, int i) ";

%feature("docstring")  main "int main() ";


// File: cpp-templates_84_8cpp.xml
%feature("docstring")  basis1 "T basis1(BasisBase< T, I, Topology,
ReferenceElement > &b, I i, T pt[Topology::dim]) ";

%feature("docstring")  basis1 "double basis1(BasisBase< double, int,
Quadrilateral4, StandardReferenceElement > &b, int i, double
pt[Quadrilateral4::dim]) ";

%feature("docstring")  main "int main() ";


// File: cpp-templates_85_8cpp.xml
%feature("docstring")  main "int main() ";


// File: cpp-templates_8cpp.xml
%feature("docstring")  main "int main() ";


// File: CPPArray_8hpp.xml


// File: CPPArray__unit_81_8cpp.xml
%feature("docstring")  main "int main() ";


// File: CPPArray__unit_81_8hpp.xml
%feature("docstring")  init "static void init(double val[2], int
dim[1]) ";


// File: CPPArray__unit_8cpp.xml
%feature("docstring")  init "static void init(double val[2], int
dim[1]) ";

%feature("docstring")  main "int main() ";


// File: Dimensions_8cpp.xml


// File: Dimensions_8hpp.xml


// File: Documentation_8hpp.xml


// File: Edge_8hpp.xml


// File: ElementOp_8hpp.xml


// File: example_8cpp.xml


// File: example_8hpp.xml


// File: example__wrap_8cpp.xml


// File: ExceptionWatch_8hpp.xml


// File: FieldFunction_8cpp.xml


// File: FieldFunction_8hpp.xml


// File: FiniteElement_81_8hpp.xml


// File: FiniteElement_8hpp.xml


// File: Fixture_8cpp.xml


// File: Fixture_8hpp.xml


// File: foo_8cpp.xml
%feature("docstring")  foo "int foo() ";


// File: Function_8cpp.xml


// File: Function_8hpp.xml


// File: FunctionOperator_8cpp.xml


// File: FunctionOperator_8hpp.xml


// File: FunctionWithIntrepidRequest_8hpp.xml


// File: GeneralFunction_8hpp.xml


// File: generated__refinement__tables_8hpp.xml


// File: GeneratedRefinementTable_8hpp.xml


// File: GenericFunction_8cpp.xml


// File: GenericFunction_8hpp.xml


// File: GeometryFactory_8cpp.xml


// File: GeometryFactory_8hpp.xml


// File: GeometryKernel_8hpp.xml


// File: GeometryKernelOpenNURBS_8cpp.xml


// File: GeometryKernelOpenNURBS_8hpp.xml


// File: GeometryKernelStupid_8hpp.xml


// File: GeometryVerifier_8cpp.xml


// File: GeometryVerifier_8hpp.xml


// File: H1Norm_8hpp.xml


// File: HasValue_8hpp.xml


// File: HeterogeneousFixture_8cpp.xml


// File: HeterogeneousFixture_8hpp.xml


// File: IAdapter_8hpp.xml


// File: IEdgeAdapter_8hpp.xml


// File: IEdgeBasedAdapterPredicate_8hpp.xml


// File: IElementAdapter_8hpp.xml


// File: IElementBasedAdapterPredicate_8hpp.xml


// File: IntegratedOp_8hpp.xml


// File: Intrepid__HGRAD__HEX__C2__Serendipity__FEM_8hpp.xml


// File: Intrepid__HGRAD__HEX__C2__Serendipity__FEMDef_8hpp.xml


// File: Intrepid__HGRAD__QUAD__C2__Serendipity__FEM_8hpp.xml


// File: Intrepid__HGRAD__QUAD__C2__Serendipity__FEMDef_8hpp.xml


// File: Intrepid__HGRAD__WEDGE__C2__Serendipity__FEM_8hpp.xml


// File: Intrepid__HGRAD__WEDGE__C2__Serendipity__FEMDef_8hpp.xml


// File: IntrepidManager_8cpp.xml


// File: IntrepidManager_8hpp.xml


// File: IsInElement_8cpp.xml


// File: IsInElement_8hpp.xml


// File: JacobianUtil_8cpp.xml


// File: JacobianUtil_8hpp.xml


// File: main-sun_8cpp.xml
%feature("docstring")  main "int main() ";


// File: main_8cpp.xml
%feature("docstring")  main "int main() ";


// File: main1_8cpp.xml
%feature("docstring")  main "int main(int argc, char **argv) ";


// File: mainsh_8cpp.xml


// File: Math_8hpp.xml


// File: MDArray_8hpp.xml


// File: MeshDifference_8cpp.xml


// File: MeshDifference_8hpp.xml


// File: MeshDifferenceMain_8cpp.xml
%feature("docstring")  main "int main(int argc, char **argv) ";


// File: MeshGeometry_8cpp.xml


// File: MeshGeometry_8hpp.xml


// File: MeshObjTopology_8hpp.xml


// File: MeshUtil_8cpp.xml


// File: MeshUtil_8hpp.xml


// File: MultipleFieldFunction_8hpp.xml


// File: Name_8hpp.xml


// File: NodeRegistry_8cpp.xml


// File: NodeRegistry_8hpp.xml


// File: NodeRegistryDef_8hpp.xml


// File: NoMallocArray_8hpp.xml


// File: Norm_8hpp.xml


// File: Observable_8hpp.xml


// File: old-shape_8hpp.xml
%feature("docstring")  main "int main() ";


// File: OptionMask_8hpp.xml


// File: ParallelUtil_8hpp.xml


// File: ParallelUtilDef_8hpp.xml


// File: PartOp_8hpp.xml


// File: Percept_8hpp.xml


// File: Percept__MOAB__SimplexTemplateRefiner_8cpp.xml


// File: Percept__MOAB__SimplexTemplateRefiner_8hpp.xml


// File: PerceptBoostArray_8hpp.xml


// File: PerceptMesh_8cpp.xml


// File: PerceptMesh_8hpp.xml


// File: PerceptMeshReadWrite_8hpp.xml


// File: PerceptMesquiteMesh_8cpp.xml


// File: PerceptMesquiteMesh_8hpp.xml


// File: PerceptMesquiteMeshDomain_8cpp.xml


// File: PerceptMesquiteMeshDomain_8hpp.xml


// File: PMMLaplaceSmoother_8hpp.xml


// File: PMMLaplaceSmoother1_8cpp.xml


// File: PMMLaplaceSmoother1_8hpp.xml


// File: PMMMsqMatrix_8hpp.xml


// File: PMMParallelReferenceMeshSmoother_8cpp.xml


// File: PMMParallelReferenceMeshSmoother_8hpp.xml


// File: PMMParallelReferenceMeshSmoother1_8cpp.xml


// File: PMMParallelReferenceMeshSmoother1_8hpp.xml


// File: PMMParallelReferenceMeshSmoother2_8cpp.xml


// File: PMMParallelReferenceMeshSmoother2_8hpp.xml


// File: PMMParallelReferenceMeshSmoother3_8cpp.xml


// File: PMMParallelReferenceMeshSmoother3_8hpp.xml


// File: PMMParallelShapeImprover_8cpp.xml


// File: PMMParallelShapeImprover_8hpp.xml


// File: PMMShapeImprover_8cpp.xml


// File: PMMShapeImprover_8hpp.xml


// File: PMMShapeSizeOrientImprover_8cpp.xml


// File: PMMShapeSizeOrientImprover_8hpp.xml


// File: PMMSmootherMetric_8hpp.xml


// File: PredicateBasedEdgeAdapter_8hpp.xml


// File: PredicateBasedElementAdapter_8hpp.xml


// File: ProgressMeter_8cpp.xml


// File: ProgressMeter_8hpp.xml


// File: PyramidFixture_8cpp.xml


// File: PyramidFixture_8hpp.xml


// File: QuadFixture_8hpp.xml


// File: RefinementInfoByType_8cpp.xml


// File: RefinementInfoByType_8hpp.xml


// File: RefinementKey_8cpp.xml


// File: RefinementKey_8hpp.xml


// File: RefinementTopology_8cpp.xml


// File: RefinementTopology_8hpp.xml


// File: Refiner_8cpp.xml


// File: Refiner_8hpp.xml


// File: RefinerPattern__Line2__Line2__N_8hpp.xml


// File: RefinerPattern__Tet4__Tet4__N_8hpp.xml


// File: RefinerPattern__Tri3__Tri3__2_8hpp.xml


// File: RefinerPattern__Tri3__Tri3__N_8hpp.xml


// File: RefinerUnrefine_8cpp.xml


// File: RefinerUtil_8cpp.xml


// File: RefinerUtil_8hpp.xml


// File: RegressionTestFileLoc_8hpp.xml


// File: RegressionTestLocalRefiner_8cpp.xml


// File: stk__percept_2regression__tests_2RegressionTestMain_8cpp.xml


// File: stk__adapt_2regression__tests_2RegressionTestMain_8cpp.xml


// File: RegressionTestMeshColorer_8cpp.xml


// File: RegressionTestNodeRegistry_8cpp.xml


// File: RegressionTestPerceptMeshFieldFunction_8cpp.xml


// File: RegressionTestSTKMeshMemory_8cpp.xml


// File: RegressionTestUniformRefiner_8cpp.xml


// File: build_81_8nogit_2packages_2rtop_2src_2RTOpPack__ROpGetSubVector_8hpp.xml


// File: build_8dir_2packages_2rtop_2src_2RTOpPack__ROpGetSubVector_8hpp.xml


// File: build_81_8nogit_2packages_2rtop_2src_2RTOpPack__RTOpSubRangeDecorator_8hpp.xml


// File: build_8dir_2packages_2rtop_2src_2RTOpPack__RTOpSubRangeDecorator_8hpp.xml


// File: build_81_8nogit_2packages_2rtop_2src_2RTOpPack__RTOpT_8hpp.xml


// File: build_8dir_2packages_2rtop_2src_2RTOpPack__RTOpT_8hpp.xml


// File: build_81_8nogit_2packages_2rtop_2src_2RTOpPack__RTOpTHelpers_8hpp.xml


// File: build_8dir_2packages_2rtop_2src_2RTOpPack__RTOpTHelpers_8hpp.xml


// File: build_81_8nogit_2packages_2rtop_2src_2RTOpPack__SPMD__apply__op_8hpp.xml


// File: build_8dir_2packages_2rtop_2src_2RTOpPack__SPMD__apply__op_8hpp.xml


// File: build_81_8nogit_2packages_2rtop_2src_2RTOpPack__TOpLinearCombination_8hpp.xml


// File: build_8dir_2packages_2rtop_2src_2RTOpPack__TOpLinearCombination_8hpp.xml


// File: RunEnvironment_8cpp.xml


// File: RunEnvironment_8hpp.xml


// File: SameRankRelation_8hpp.xml


// File: Searcher_8hpp.xml


// File: SerializeNodeRegistry_8hpp.xml


// File: ShardsInterfaceTable_8cpp.xml


// File: ShardsInterfaceTable_8hpp.xml


// File: SimpleSearcher_8cpp.xml


// File: SimpleSearcher_8hpp.xml


// File: SingleTetFixture_8cpp.xml


// File: SingleTetFixture_8hpp.xml


// File: SpacingFieldUtil_8cpp.xml


// File: SpacingFieldUtil_8hpp.xml


// File: StdMeshObjTopologies_8cpp.xml


// File: StdMeshObjTopologies_8hpp.xml


// File: stk__mesh_8hpp.xml


// File: stk__percept__code__types_8hpp.xml


// File: STKSearcher_8hpp.xml


// File: STKSearcherDef_8hpp.xml


// File: StringFunction_8cpp.xml


// File: StringFunction_8hpp.xml


// File: SubDimCell_8hpp.xml


// File: SweepMesher_8cpp.xml


// File: SweepMesher_8hpp.xml


// File: TestLocalRefiner_8cpp.xml


// File: TestLocalRefiner_8hpp.xml


// File: TestLocalRefinerTet__N__1_8hpp.xml


// File: TestLocalRefinerTet__N__2_8hpp.xml


// File: TestLocalRefinerTet__N__2__1_8hpp.xml


// File: TestLocalRefinerTet__N__3_8hpp.xml


// File: TestLocalRefinerTet__N__3__1_8cpp.xml


// File: TestLocalRefinerTet__N__3__1_8hpp.xml


// File: TestLocalRefinerTet__N__4_8hpp.xml


// File: TestLocalRefinerTri_8cpp.xml


// File: TestLocalRefinerTri_8hpp.xml


// File: TestLocalRefinerTri1_8cpp.xml


// File: TestLocalRefinerTri1_8hpp.xml


// File: TestLocalRefinerTri2_8cpp.xml


// File: TestLocalRefinerTri2_8hpp.xml


// File: TestLocalRefinerTri__N_8cpp.xml


// File: TestLocalRefinerTri__N_8hpp.xml


// File: TestLocalRefinerTri__N__1_8cpp.xml


// File: TestLocalRefinerTri__N__1_8hpp.xml


// File: TestLocalRefinerTri__N__2_8cpp.xml


// File: TestLocalRefinerTri__N__2_8hpp.xml


// File: TestLocalRefinerTri__N__3_8hpp.xml


// File: TestLocalRefinerTri__N__3__EdgeBasedAnisotropic_8hpp.xml


// File: TestLocalRefinerTri__N__3__IEdgeAdapter_8hpp.xml


// File: TestLocalRefinerTri__N__3__IElementAdapter_8hpp.xml


// File: TestLocalRefinerTri__N__3__MeshSizeRatio_8hpp.xml


// File: testpgi_8cpp.xml
%feature("docstring")  mymainpgi "int mymainpgi(int argc, char
**argv) ";

%feature("docstring")  main "int main(int argc, char **argv) ";


// File: testpgi1_8cpp.xml
%feature("docstring")  mymainpgi "int mymainpgi(int argc, char
**argv) ";


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__apply__op__helper_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__apply__op__helper_8hpp.xml


// File: build_81_8nogit_2packages_2stratimikos_2adapters_2belos_2src_2Thyra__BelosLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2stratimikos_2adapters_2belos_2src_2Thyra__BelosLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2stratimikos_2adapters_2belos_2src_2Thyra__BelosLinearOpWithSolveFactory_8hpp.xml


// File: build_8dir_2packages_2stratimikos_2adapters_2belos_2src_2Thyra__BelosLinearOpWithSolveFactory_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultAddedLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultAddedLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultAdjointLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultAdjointLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultBlockedLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultBlockedLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultBlockedTriangularLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultBlockedTriangularLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultBlockedTriangularLinearOpWithSolveFactory_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultBlockedTriangularLinearOpWithSolveFactory_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultClusteredSpmdProductVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultClusteredSpmdProductVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultClusteredSpmdProductVectorSpace_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultClusteredSpmdProductVectorSpace_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultColumnwiseMultiVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultColumnwiseMultiVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultDiagonalLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultDiagonalLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultDiagonalLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultDiagonalLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultFiniteDifferenceModelEvaluator_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultFiniteDifferenceModelEvaluator_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultIdentityLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultIdentityLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultInverseLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultInverseLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultLinearOpSource_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultLinearOpSource_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultMultipliedLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultMultipliedLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultMultiVectorLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultMultiVectorLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultMultiVectorProductVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultMultiVectorProductVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultMultiVectorProductVectorSpace_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultMultiVectorProductVectorSpace_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultPreconditioner_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultPreconditioner_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultProductMultiVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultProductMultiVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultProductVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultProductVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultProductVectorSpace_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultProductVectorSpace_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultScaledAdjointLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultScaledAdjointLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultSerialDenseLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultSerialDenseLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultSerialDenseLinearOpWithSolveFactory_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultSerialDenseLinearOpWithSolveFactory_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdMultiVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdMultiVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdVector_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdVector_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdVectorSpace_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdVectorSpace_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdVectorSpaceFactory_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultSpmdVectorSpaceFactory_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DefaultZeroLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DefaultZeroLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DelayedLinearOpWithSolve_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DelayedLinearOpWithSolve_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DelayedLinearOpWithSolveFactory_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DelayedLinearOpWithSolveFactory_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__describeLinearOp_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__describeLinearOp_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2test_2nonlinear_2models_2Thyra__DiagonalQuadraticResponseOnlyModelEvaluator_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2test_2nonlinear_2models_2Thyra__DiagonalQuadraticResponseOnlyModelEvaluator_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2test_2nonlinear_2models_2Thyra__DiagonalScalarProd_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2test_2nonlinear_2models_2Thyra__DiagonalScalarProd_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__DirectionalFiniteDiffCalculator_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__DirectionalFiniteDiffCalculator_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__EuclideanScalarProd_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__EuclideanScalarProd_8hpp.xml


// File: build_81_8nogit_2packages_2stratimikos_2adapters_2belos_2src_2Thyra__GeneralSolveCriteriaBelosStatusTest_8hpp.xml


// File: build_8dir_2packages_2stratimikos_2adapters_2belos_2src_2Thyra__GeneralSolveCriteriaBelosStatusTest_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpDefaultBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpDefaultBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpScalarProd_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpScalarProd_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpTester_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpWithSolveBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpWithSolveBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpWithSolveFactoryBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpWithSolveFactoryBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__LinearOpWithSolveTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__LinearOpWithSolveTester_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__ModelEvaluatorBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__ModelEvaluatorBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__MultiVectorAdapterBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__MultiVectorAdapterBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__MultiVectorBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__MultiVectorBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__MultiVectorDefaultBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__MultiVectorDefaultBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__MultiVectorStdOps_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__MultiVectorStdOps_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__MultiVectorStdOpsTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__MultiVectorStdOpsTester_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__MultiVectorTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__MultiVectorTester_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__PreconditionerFactoryBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__PreconditionerFactoryBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__ScalarProdBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__ScalarProdBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__ScalarProdVectorSpaceBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__ScalarProdVectorSpaceBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__ScaledAdjointLinearOpBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__ScaledAdjointLinearOpBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__ScaledModelEvaluator_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__ScaledModelEvaluator_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2test_2nonlinear_2models_2Thyra__Simple2DModelEvaluator_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2test_2nonlinear_2models_2Thyra__Simple2DModelEvaluator_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__SpmdMultiVectorBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__SpmdMultiVectorBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__SpmdMultiVectorSerializer_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__SpmdMultiVectorSerializer_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__SpmdVectorBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__SpmdVectorBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__SpmdVectorSpaceBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__SpmdVectorSpaceBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__SpmdVectorSpaceDefaultBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__SpmdVectorSpaceDefaultBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorDefaultBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorDefaultBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorSpaceBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorSpaceBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorSpaceDefaultBase_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorSpaceDefaultBase_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorSpaceTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorSpaceTester_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorStdOps_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorStdOps_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorStdOpsTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorStdOpsTester_8hpp.xml


// File: build_81_8nogit_2packages_2thyra_2core_2src_2Thyra__VectorTester_8hpp.xml


// File: build_8dir_2packages_2thyra_2core_2src_2Thyra__VectorTester_8hpp.xml


// File: TopologyVerifier_8cpp.xml


// File: TopologyVerifier_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__BlockCrsGraph_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__BlockCrsGraph_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__BlockMap_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__BlockMap_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__BlockMultiVector_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__BlockMultiVector_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__CrsGraph_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__CrsGraph_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__CrsMatrix_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__CrsMatrix_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__CrsMatrixMultiplyOp_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__CrsMatrixMultiplyOp_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__CrsMatrixSolveOp_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__CrsMatrixSolveOp_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__Directory_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__Directory_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__Map_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__Map_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2inout_2Tpetra__MatrixIO_8hpp.xml


// File: build_8dir_2packages_2tpetra_2inout_2Tpetra__MatrixIO_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__MultiVector_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__MultiVector_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__RowMatrixTransposer_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__RowMatrixTransposer_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__VbrMatrix_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__VbrMatrix_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2src_2Tpetra__Vector_8hpp.xml


// File: build_8dir_2packages_2tpetra_2src_2Tpetra__Vector_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2ext_2TpetraExt__BlockExtraction_8hpp.xml


// File: build_8dir_2packages_2tpetra_2ext_2TpetraExt__BlockExtraction_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2ext_2TpetraExt__MatrixMatrix_8hpp.xml


// File: build_8dir_2packages_2tpetra_2ext_2TpetraExt__MatrixMatrix_8hpp.xml


// File: build_81_8nogit_2packages_2tpetra_2ext_2TpetraExt__MMHelpers_8hpp.xml


// File: build_8dir_2packages_2tpetra_2ext_2TpetraExt__MMHelpers_8hpp.xml


// File: TransformPath_8hpp.xml


// File: UniformRefiner_8cpp.xml


// File: UniformRefiner_8hpp.xml


// File: UniformRefinerPattern_8cpp.xml


// File: UniformRefinerPattern_8hpp.xml


// File: UniformRefinerPattern__Beam2__Beam2__2__sierra_8hpp.xml


// File: UniformRefinerPattern__Beam2__Beam3__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Beam3__Beam3__2__sierra_8hpp.xml


// File: UniformRefinerPattern__Hex20__Hex20__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Hex27__Hex27__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Hex8__Hex20__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Hex8__Hex27__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Hex8__Hex8__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Hex8__Tet4__24_8hpp.xml


// File: UniformRefinerPattern__Hex8__Tet4__6__12_8hpp.xml


// File: UniformRefinerPattern__Line2__Line2__2__sierra_8hpp.xml


// File: UniformRefinerPattern__Line2__Line3__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Line3__Line3__2__sierra_8hpp.xml


// File: UniformRefinerPattern__Pyramid13__Pyramid13__10__sierra_8hpp.xml


// File: UniformRefinerPattern__Pyramid5__Pyramid13__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Pyramid5__Pyramid5__10__sierra_8hpp.xml


// File: UniformRefinerPattern__Quad4__Quad4__4_8hpp.xml


// File: UniformRefinerPattern__Quad4__Quad4__4__sierra_8hpp.xml


// File: UniformRefinerPattern__Quad4__Quad8__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Quad4__Quad9__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Quad4__Tri3__2_8hpp.xml


// File: UniformRefinerPattern__Quad4__Tri3__4_8hpp.xml


// File: UniformRefinerPattern__Quad4__Tri3__6_8hpp.xml


// File: UniformRefinerPattern__Quad8__Quad8__4__sierra_8hpp.xml


// File: UniformRefinerPattern__Quad9__Quad9__4__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellLine2__ShellLine2__2__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellLine2__ShellLine3__1__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellLine3__ShellLine3__2__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellQuad4__ShellQuad4__4__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellQuad4__ShellQuad8__1__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellQuad4__ShellQuad9__1__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellQuad8__ShellQuad8__4__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellTri3__ShellTri3__4__sierra_8hpp.xml


// File: UniformRefinerPattern__ShellTri6__ShellTri6__4__sierra_8hpp.xml


// File: UniformRefinerPattern__Tet10__Tet10__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Tet4__Tet10__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Tet4__Tet4__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Tri3__Tri3__4__sierra_8hpp.xml


// File: UniformRefinerPattern__Tri3__Tri6__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Tri6__Tri6__4__sierra_8hpp.xml


// File: UniformRefinerPattern__Wedge15__Wedge15__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Wedge18__Wedge18__8__sierra_8hpp.xml


// File: UniformRefinerPattern__Wedge6__Wedge15__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Wedge6__Wedge18__1__sierra_8hpp.xml


// File: UniformRefinerPattern__Wedge6__Wedge6__8__sierra_8hpp.xml


// File: UnitTestFieldFunction_8cpp.xml


// File: UnitTestGeometryVerifier_8cpp.xml


// File: UnitTestLocalRefiner_8cpp.xml


// File: stk__percept_2unit__tests_2UnitTestMain_8cpp.xml


// File: stk__adapt_2unit__tests_2UnitTestMain_8cpp.xml


// File: UnitTestMeshColorer_8cpp.xml


// File: UnitTestNodeRegistry_8cpp.xml


// File: UnitTestNorm_8cpp.xml


// File: UnitTestPerceptMesh_8cpp.xml


// File: UnitTestPerceptMesquiteMesh_8cpp.xml


// File: UnitTestSearch_8cpp.xml


// File: UnitTestStringFunction_8cpp.xml


// File: UnitTestSubDimCell_8cpp.xml


// File: UnitTestSupport_8cpp.xml


// File: UnitTestSupport_8hpp.xml


// File: UnitTestTimeMaps_8cpp.xml


// File: UnitTestTopoCheck_8cpp.xml
%feature("docstring")  stk_percept_unit::isParallel "static bool
isParallel() ";


// File: UnitTestUniformRefiner_8cpp.xml


// File: UnitTestUtilities_8cpp.xml
%feature("docstring")  stk_example_io::myMain "int myMain(int argc,
char **argv) ";


// File: URP__Heterogeneous__3D_8hpp.xml


// File: URP__Heterogeneous__Enrich__3D_8hpp.xml


// File: URP__Heterogeneous__QuadraticRefine__3D_8hpp.xml


// File: Util_8cpp.xml


// File: Util_8hpp.xml


// File: Verifier_8cpp.xml


// File: Verifier_8hpp.xml


// File: WedgeFixture_8hpp.xml


// File: group__stk__io__module.xml


// File: parallel__node__registry.xml


// File: deprecated.xml


// File: dir_43b85d1b9756cd92573a05b256534e97.xml


// File: dir_39663f98425a89e72081ace7340e200d.xml


// File: dir_28fe238ec277346dd5ac9abf6486ecdc.xml


// File: dir_c4f432dcc56274f262ee25f6c7e2797d.xml


// File: dir_4cd0eb8fce7a6ee8ac0f33e12658ea31.xml


// File: dir_6d019fe9a3e912022e34420a98559810.xml


// File: dir_93acfb67a15c26e3fe5fe49fc9b7326b.xml


// File: dir_74894c83a9b783e238234574b90570d4.xml


// File: dir_ba265262617a39be59e5c0b75c304539.xml


// File: dir_165c6557e1301ecb3c9f5905eec3b00b.xml


// File: dir_88acda998ade427962a1aa8508daf788.xml


// File: dir_920c0bfc63489478efe34b800b46c97f.xml


// File: dir_0443bedd395825987c9352c6b26347a0.xml


// File: dir_c7c5b2cad158fb171fec82ab94b075cc.xml


// File: dir_33761180d013d3a426bf969aabed0981.xml


// File: dir_87355dd6c023fe3112a1eb08b5bd3d90.xml


// File: dir_97ccff9159689a66d84f5690d42f6e91.xml


// File: dir_9ef0d9f94c5b08ad32905ee42c99da0b.xml


// File: dir_656d53a56da6f6249384b0180a8a3dc9.xml


// File: dir_0fb21209bb7281c17d4777c1420e9f9d.xml


// File: dir_5469e35ff415ec73490a0d10cd6a4da4.xml


// File: dir_a0f0fc4a584be87a71e7cbd43b9a7c3b.xml


// File: dir_36ebc1658037cfbe2a3bddb93ac684c2.xml


// File: dir_cd830b685ea2bb92ad9375e95f5a1d8c.xml


// File: dir_84ad5f77c193d063050821de314d5ffb.xml


// File: dir_242febc40ebf037f72d55dec0deae29e.xml


// File: dir_dd4251f6c1cd951b1931a0d3f29c3170.xml


// File: dir_c371071c4adb1bdceb3bdd8c7b9f6296.xml


// File: dir_80e6ada6a090e8cbadc3d76932d1a232.xml


// File: dir_224a14681bb5c169504c3ecec1f071a7.xml


// File: dir_1ee449185ae9f96b294df79597402ef9.xml


// File: dir_e64ab482b65c8c2adc5f5c6863c92f30.xml


// File: dir_df3262facab5696542d8ae3effebf9e8.xml


// File: dir_15365b9e2d352a1ca0edd7cbf7a3464a.xml


// File: dir_410c1fe22e0084673b69062c6c54d879.xml


// File: dir_073c6910f180136e38f3faaac4e2a78d.xml


// File: dir_ff477d759c31d91d840e8f68d4d3307e.xml


// File: dir_702756e862c498edd6b289e650bd4f23.xml


// File: dir_d44f8c6f9fb6d7790031ce62a1900ecc.xml


// File: dir_3108127aafc8f30586aa41e578b1c88a.xml


// File: dir_45ee8ae98f67f343e71aedeb361668fe.xml


// File: dir_2165f1784ae20b5da3935b09d83f066e.xml


// File: dir_d558d26576acd1c39a83ffd0505d3402.xml


// File: dir_5ac9855eb2a152d35b8c8730350743d8.xml


// File: dir_33bc1e7d328ef6eaf2b5b327efca805e.xml


// File: dir_fc26ebcb0c645981cb8a320bc0cf8db6.xml


// File: dir_128b25ae309c0099aa78423b81222b12.xml


// File: dir_e349b3b9e712ee09db2a3cb7f9ac55ed.xml


// File: dir_b236f165722599467367314bd533ebb3.xml


// File: dir_642888a80233615f89234a856e08879d.xml


// File: dir_ed42a66720eeb3d62cedc9bc0b4314bc.xml


// File: dir_036f272561d74278959e0ca0a17c87ce.xml


// File: dir_b7b6fa893d7aa49fe64cc0c1725eb058.xml


// File: dir_8476efa39462d68fa0fcb0fd9fbc6f41.xml


// File: dir_bb76877bacd7903b45d18d5b03d8db4d.xml


// File: dir_745a945aedff9456fe32200265e844a9.xml


// File: dir_54419b89dba1c3538f7068ee95c894c4.xml


// File: dir_e39d2153e012a72ddf68a6c835aa1fb2.xml


// File: dir_3befe1e5fb3228176e604f3801f5b895.xml


// File: dir_534e19014992359ef56b8dac61596d0e.xml


// File: dir_be24f826505eaa65a3c6b4a01e254fe2.xml


// File: dir_c3265abd3a68485f5633c533094a6871.xml


// File: dir_b0b96e322f3cbdb53a38b8603baaa664.xml


// File: dir_b88e0f0f53cb1d5e809f45e5b3718b5a.xml


// File: dir_0cb005c1013ddf5c357d41a0c66cb4a8.xml


// File: dir_eec3f98673451bab17591b540eb377d4.xml


// File: dir_3a8e978858dd8ded4ec2aa714334a98d.xml


// File: dir_64589e2da2495338c450da853bf13b83.xml


// File: dir_e6076efb2092f98e23962f456d0814ea.xml


// File: dir_3f187bb3e470ae6aaf7c43d405290954.xml


// File: indexpage.xml


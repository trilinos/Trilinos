#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_BOUNDINGBOXMESH_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_BOUNDINGBOXMESH_HPP_

#include <Akri_BoundingBoxMesh.hpp>
#include <stk_topology/topology.hpp>

namespace krino {

class BoundingBoxMeshOnWorld : public BoundingBoxMesh
{
public:
  BoundingBoxMeshOnWorld(stk::topology elemTopology) : BoundingBoxMesh(elemTopology, MPI_COMM_WORLD) {} // CHECK: ALLOW MPI_COMM_WORLD
};

class BoundingBoxMeshTri3 : public BoundingBoxMeshOnWorld
{
public:
  BoundingBoxMeshTri3() : BoundingBoxMeshOnWorld(stk::topology::TRIANGLE_3_2D) {}
};

class BoundingBoxMeshQuad4 : public BoundingBoxMeshOnWorld
{
public:
  BoundingBoxMeshQuad4() : BoundingBoxMeshOnWorld(stk::topology::QUADRILATERAL_4_2D) {}
};

class BoundingBoxMeshTet4 : public BoundingBoxMeshOnWorld
{
public:
  BoundingBoxMeshTet4() : BoundingBoxMeshOnWorld(stk::topology::TETRAHEDRON_4) {}
};

class BoundingBoxMeshHex8 : public BoundingBoxMeshOnWorld
{
public:
  BoundingBoxMeshHex8() : BoundingBoxMeshOnWorld(stk::topology::HEXAHEDRON_8) {}
};

}

#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_BOUNDINGBOXMESH_HPP_ */

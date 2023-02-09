#ifndef MESH_RELATIONAL_DATA_H
#define MESH_RELATIONAL_DATA_H

#include "field.h"
#include "mesh.h"
#include "predicates/intersection_common.h"
#include "projection.h"
#include "variable_size_field.h"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

struct FakeVert
{
    int id;
    utils::Point pt;
};

inline bool operator<(const FakeVert& lhs, const FakeVert& rhs)
{
  return lhs.id < rhs.id;
}

class FakeVertGenerator
{
  public:
    FakeVert get_vert() { return {m_currId++, utils::Point()}; }

    FakeVert get_vert(const utils::Point& pt) { return {m_currId++, pt}; }

    int get_num_verts() const { return m_currId; }

  private:
    int m_currId = 0;
};

// TODO: get rid of EdgeSplitRecord in favor of VertOnEdge?
struct EdgeSplitRecord
{
    FakeVert vert;
    double xi;
};

struct VertOnEdge
{
    FakeVert vert;
    double xi;
};

using stk::middle_mesh::mesh::impl::FieldShape;

struct MeshRelationalData
{
    MeshRelationalData(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                       std::shared_ptr<mesh::Mesh> meshIn)
      : vertsInClassOnMesh1(mesh::create_field<predicates::impl::PointRecord>(meshIn, FieldShape(1, 0, 0), 1))
      , verts2ClassOnMesh1(mesh::create_field<predicates::impl::PointRecord>(mesh2, FieldShape(1, 0, 0), 1))
      , verts1ToFakeVerts(mesh::create_field<FakeVert>(mesh1, FieldShape(1, 0, 0), 1))
      , verts2ToFakeVerts(mesh::create_field<FakeVert>(mesh2, FieldShape(1, 0, 0), 1))
      , edges2ToFakeVertsIn(mesh::create_variable_size_field<VertOnEdge>(mesh2, FieldShape(0, 1, 0)))
      , mesh1EdgesToSplit(mesh::create_variable_size_field<EdgeSplitRecord>(mesh1, FieldShape(0, 1, 0)))
      , mesh1ElsToVertsIn(mesh::create_variable_size_field<mesh::MeshEntityPtr>(mesh1, FieldShape(0, 0, 1)))
      , meshInElementsToMesh1Elements(mesh::create_field<mesh::MeshEntityPtr>(meshIn, FieldShape(0, 0, 1), 1, nullptr))
      , meshInElementsToMesh2Elements(mesh::create_field<mesh::MeshEntityPtr>(meshIn, FieldShape(0, 0, 1), 1, nullptr))
      , mesh1ElementsToMesh2Elements(mesh::create_variable_size_field<mesh::MeshEntityPtr>(mesh1, FieldShape(0, 0, 1)))
    // edges1_to_in(mesh::create_variable_size_field<MeshEntityPtr>(mesh1, FieldShape(0, 1, 0)))
    {}

    mesh::FieldPtr<predicates::impl::PointRecord> vertsInClassOnMesh1; // classification of mesh_in vertices on mesh1
    mesh::FieldPtr<predicates::impl::PointRecord> verts2ClassOnMesh1;  // classification of mesh2 verts on mesh2
    mesh::FieldPtr<FakeVert> verts1ToFakeVerts; // maps mesh1 verts to placeholder verts that will be created
                                                // on mesh_in
    mesh::FieldPtr<FakeVert> verts2ToFakeVerts; // maps mesh2 verts to placeholder verts that will be created
                                                // on mesh_in
    mesh::VariableSizeFieldPtr<VertOnEdge> edges2ToFakeVertsIn; // TODO: the _in suffix is not necessary:
                                                                //       FakeVerts are always in mesh_in
    mesh::VariableSizeFieldPtr<EdgeSplitRecord>
        mesh1EdgesToSplit; // edges on mesh1 that are split by mesh2 edges or vertices

    // VariableSizeFieldPtr<MeshEntityPtr> edges1_to_in;  // edges on mesh_in that are on mesh1 edges
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
        mesh1ElsToVertsIn;                               // maps the elements of mesh1 to the vertices
                                                         // on mesh_in that are inside them (or on their
                                                         // closure)
                                                         // TODO: get rid of this field: store a seed
                                                         //       vertex (a mesh1 element corner vertex)
                                                         //       and do an search outwards until no more
                                                         //       verts classified on the element are found
    std::vector<mesh::MeshEntityPtr> fakeVertsToVertsIn; // map of fake_verts to mesh_in verts once they
                                                         // are created
    mesh::FieldPtr<mesh::MeshEntityPtr> meshInElementsToMesh1Elements; // maps mesh_in elements to mesh1 elements
    mesh::FieldPtr<mesh::MeshEntityPtr> meshInElementsToMesh2Elements; // maps mesh_in elements to mesh2 elements
    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
        mesh1ElementsToMesh2Elements; // stores the candidate set of mesh2 elements
                                      // that *may* overlap a given mesh1 element.
                                      // This set of elements is required to cover
                                      // the mesh1 elements, but may contain additional
                                      // elements that are nearby
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
#ifdef STK_BUILT_IN_SIERRA

#include "create_stk_mesh.hpp"
#include "field.hpp"
#include "constants.hpp"
#include "stk_mesh/base/FEMHelpers.hpp" //TODO: ifdef this

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

void StkMeshCreator::declare_stk_vert_field()
{
  std::string name = NAME_PREFIX + "nonconformal_interface_vert_field";
  m_metaDataPtr->enable_late_fields();
  m_stkNodeField = &(m_metaDataPtr->declare_field<FieldScalarType>(stk::topology::NODE_RANK, name));
}

void StkMeshCreator::load_mesh(const std::string& fname)
{
  stk::io::StkMeshIoBroker reader(m_bulkDataPtr->parallel());
  reader.set_bulk_data(*m_bulkDataPtr);
  reader.add_mesh_database(fname, stk::io::READ_MESH);
  reader.create_input_mesh();
  reader.populate_bulk_data();
}

MeshPart StkMeshCreator::create_mesh_from_part(const std::string& name)
{
  std::cout << "creating mesh from part " << name << std::endl;
  std::shared_ptr<mesh::Mesh> mesh = mesh::make_empty_mesh();
  MeshFieldPtr stkEls =
      mesh::create_field<MeshFieldPtr::element_type::value_type>(mesh, mesh::impl::FieldShape(0, 0, 1), 1);
  m_part = m_metaDataPtr->get_part(name);
  stk::mesh::put_field_on_mesh(*m_stkNodeField, *m_part, 1, 0);

  create_nodes(mesh);
  if (m_bulkDataPtr->does_sideset_exist(*m_part))
    create_faces_from_sideset(mesh, stkEls);
  else
    create_faces_from_shells(mesh, stkEls);

  return {mesh, stkEls};
}

// copies the coordinates from the Mesh to the STK mesh
void StkMeshCreator::write_back_coords(std::shared_ptr<mesh::Mesh> mesh, const std::string& name)
{
  m_part                = m_metaDataPtr->get_part(name);
  const auto& meshVerts = mesh->get_vertices();

  // TODO: distinguish between NODE_RANK and VERT_RANK
  const stk::mesh::FieldBase& coordField = *(m_metaDataPtr->coordinate_field());
  stk::mesh::Selector select             = *m_part;
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, select);

  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& vert : *bucket)
    {
      FieldScalarType* vertIdxs = (stk::mesh::field_data(*m_stkNodeField, vert));
      mesh::MeshEntityPtr vert2 = meshVerts[vertIdxs[0]];
      utils::Point pt           = vert2->get_point_orig(0);

      double* coordsV = reinterpret_cast<double*>(stk::mesh::field_data(coordField, vert));
      coordsV[0]      = pt.x;
      coordsV[1]      = pt.y;
      coordsV[2]      = pt.z;
    }
}

void StkMeshCreator::create_nodes(std::shared_ptr<mesh::Mesh> mesh)
{
  // TODO: distinguish between NODE_RANK and VERT_RANK
  const stk::mesh::FieldBase& coordField = *(m_metaDataPtr->coordinate_field());
  stk::mesh::Selector select             = *m_part;
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, select);

  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& vert : *bucket)
    {
      double* coordsV = reinterpret_cast<double*>(stk::mesh::field_data(coordField, vert));

      auto vert2 = mesh->create_vertex(coordsV[0], coordsV[1], coordsV[2]);

      // record relation between meshes
      FieldScalarType* nodeFieldV = stk::mesh::field_data(*m_stkNodeField, vert);
      nodeFieldV[0]               = vert2->get_id();
    }
}
/*
void StkMeshCreator::createFaces(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stk_els)
{
  const auto& mesh_verts = mesh->get_vertices();
  stk::mesh::Selector select = *m_part;
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::FACE_RANK, select);

  std::array<mesh::MeshEntityPtr, 4> verts;
  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& face : *bucket)
    {
      std::cout << "\nface " << face << std::endl;
      int i = 0;
      for (auto it = m_bulkDataPtr->begin_nodes(face); it != m_bulkDataPtr->end_nodes(face); ++it)
      {
        FieldScalarType* vert_idxs = (stk::mesh::field_data(*m_stkNodeField, *it));
        verts.at(i) = mesh_verts[vert_idxs[0]];
        ++i;
      }

      assert(i == 4);
      mesh::MeshEntityPtr el = mesh->create_quad_from_verts(verts[0], verts[1],
                                                   verts[2], verts[3]);
      std::cout << "creating quad from verts " << verts[0]->get_point_orig(0)
                << ", " << verts[1]->get_point_orig(0)
                << ", " << verts[2]->get_point_orig(0)
                << ", " << verts[3]->get_point_orig(0) << std::endl;
      (*stk_els)(el, 0, 0) = face;
    }
}
*/

void StkMeshCreator::create_faces_from_sideset(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls)
{
  stk::mesh::Entity invalidEntity;
  assert(m_bulkDataPtr->does_sideset_exist(*m_part));
  const stk::mesh::SideSet& sideset = m_bulkDataPtr->get_sideset(*m_part);
  const auto& meshVerts             = mesh->get_vertices();

  std::array<mesh::MeshEntityPtr, 4> verts;
  for (auto sIt = sideset.begin(); sIt != sideset.end(); ++sIt)
  {
    stk::mesh::SideSetEntry entry = *sIt;
    stk::mesh::Entity face        = get_side_entity_for_elem_side_pair(*m_bulkDataPtr, entry.element, entry.side);
    assert(face != invalidEntity); // check entity is valid
    int i = 0;
    for (auto it = m_bulkDataPtr->begin_nodes(face); it != m_bulkDataPtr->end_nodes(face); ++it)
    {
      FieldScalarType* vertIdxs = (stk::mesh::field_data(*m_stkNodeField, *it));
      verts.at(i)               = meshVerts[vertIdxs[0]];
      ++i;
    }

    assert(i == 4);
    mesh::MeshEntityPtr el = mesh->create_quad_from_verts(verts[0], verts[1], verts[2], verts[3]);
    (*stkEls)(el, 0, 0)    = entry;
  }
}

void StkMeshCreator::create_faces_from_shells(std::shared_ptr<mesh::Mesh> mesh, MeshFieldPtr stkEls)
{
  const auto& meshVerts = mesh->get_vertices();

  stk::mesh::Selector select             = *m_part;
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::ELEMENT_RANK, select);

  stk::topology::topology_t topo = m_part->topology();
  if (!(topo == stk::topology::SHELL_TRI_3 || topo == stk::topology::SHELL_QUAD_4))
    throw std::runtime_error(
        "part topology is not SHELL_TRI_3 or SHELL_QUAD_4.  These are the only supported topologies");

  std::array<mesh::MeshEntityPtr, 4> verts;
  for (stk::mesh::Bucket* bucket : buckets)
    for (auto& stkEl : *bucket)
    {
      size_t idx = 0;
      for (auto it = m_bulkDataPtr->begin_nodes(stkEl); it != m_bulkDataPtr->end_nodes(stkEl); ++it)
      {
        assert(idx < verts.size());
        FieldScalarType* vertIdxs = (stk::mesh::field_data(*m_stkNodeField, *it));
        verts[idx]                = meshVerts[vertIdxs[0]];
        ++idx;
      }

      mesh::MeshEntityPtr meshEl = nullptr;
      if (topo == stk::topology::SHELL_TRI_3)
        meshEl = mesh->create_triangle_from_verts(verts[0], verts[1], verts[2]);
      else if (topo == stk::topology::SHELL_QUAD_4)
        meshEl = mesh->create_quad_from_verts(verts[0], verts[1], verts[2], verts[3]);
      else
        throw std::runtime_error("found element that is not triangle or quad");

      (*stkEls)(meshEl, 0, 0) = stk::mesh::SideSetEntry(stkEl);
    }
}

} // namespace impl
} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

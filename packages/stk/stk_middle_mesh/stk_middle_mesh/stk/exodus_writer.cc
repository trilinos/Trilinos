#ifdef STK_BUILT_IN_SIERRA

#include "stk/exodus_writer.h"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

void ExodusWriter::initialize_output_mesh()
{
  using CoordFieldType       = stk::mesh::Field<double, stk::mesh::Cartesian>;
  CoordFieldType& coordField = m_metaDataOutPtr->declare_field<CoordFieldType>(stk::topology::NODE_RANK, "coord_field");
  stk::io::set_field_role(coordField, Ioss::Field::MESH);
  m_metaDataOutPtr->set_coordinate_field(&coordField);

  stk::mesh::put_field_on_mesh(coordField, m_metaDataOutPtr->universal_part(), m_metaDataOutPtr->spatial_dimension(),
                               (stk::mesh::FieldTraits<CoordFieldType>::data_type*)nullptr);
}

void ExodusWriter::create_part(const std::string& name)
{
  m_part = &(m_metaDataOutPtr->declare_part_with_topology(name, stk::topology::SHELL_TRI_3));
  stk::io::put_io_part_attribute(*m_part);
}

mesh::FieldPtr<stk::mesh::Entity> ExodusWriter::create_part_verts()
{
  std::vector<stk::mesh::Part*> partVec{m_part};
  auto stkVerts = mesh::create_field<stk::mesh::Entity>(m_mesh, mesh::impl::FieldShape(1, 0, 0), 1);

  m_bulkDataOutPtr->modification_begin();

  std::vector<size_t> requests(m_metaDataOutPtr->entity_rank_count(), 0);
  requests[stk::topology::NODE_RANK] = count_valid(m_mesh->get_vertices());
  std::vector<stk::mesh::Entity> entities;
  m_bulkDataOutPtr->generate_new_entities(requests, entities);
  for (auto& e : entities)
    m_bulkDataOutPtr->change_entity_parts(e, partVec);

  m_bulkDataOutPtr->modification_end();

  const stk::mesh::FieldBase& coordField = *(m_metaDataOutPtr->coordinate_field());
  int entityIdx                          = 0;
  for (auto& vert : m_mesh->get_vertices())
    if (vert)
    {
      double* coords  = static_cast<double*>(stk::mesh::field_data(coordField, entities[entityIdx]));
      utils::Point pt = vert->get_point_orig(0);

      coords[0]               = pt.x;
      coords[1]               = pt.y;
      coords[2]               = pt.z;
      (*stkVerts)(vert, 0, 0) = entities[entityIdx];

      entityIdx++;
    }

  return stkVerts;
}

void ExodusWriter::create_part_elements(mesh::FieldPtr<stk::mesh::Entity> stkVerts)
{
  // create elements:
  //   get Mesh element
  //   get Mesh vertices
  //   get corresponding Stk vertices
  //   create Stk face
  //   call declare_relation(stk_face, stk_vert)

  int nels = count_valid(m_mesh->get_elements());

  m_bulkDataOutPtr->modification_begin();

  std::vector<stk::mesh::Entity> entities;
  std::vector<size_t> requests(m_metaDataOutPtr->entity_rank_count(), 0);
  requests[stk::topology::ELEMENT_RANK] = nels;
  m_bulkDataOutPtr->generate_new_entities(requests, entities);

  std::vector<stk::mesh::Part*> partVec{m_part};
  int entityIdx = 0;
  for (auto& el : m_mesh->get_elements())
    if (el)
    {
      mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
      int ndown = get_downward(el, 0, verts);
      assert(ndown == 3);

      m_bulkDataOutPtr->change_entity_parts(entities[entityIdx], partVec);
      for (int j = 0; j < ndown; ++j)
      {
        auto stkVert = (*stkVerts)(verts[j], 0, 0);
        m_bulkDataOutPtr->declare_relation(entities[entityIdx], stkVert, j);
      }

      entityIdx++;
    }

  m_bulkDataOutPtr->modification_end();
}

void ExodusWriter::write_output_stk(const std::string& fname)
{
  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(m_bulkDataOutPtr);
  auto fh = stkIo.create_output_mesh(fname, stk::io::WRITE_RESULTS);
  stkIo.write_output_mesh(fh);
}

} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk
#endif

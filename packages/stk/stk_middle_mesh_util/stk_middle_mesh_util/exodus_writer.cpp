#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"

#include "stk_middle_mesh_util/field_output_adaptor.hpp"
#include "stk_middle_mesh_util/exodus_writer.hpp"
#include "stk_io/WriteMesh.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

void ExodusWriter::initialize_output_mesh()
{
  using CoordFieldType = stk::mesh::Field<double>;
  CoordFieldType& coordField = m_metaDataOutPtr->declare_field<double>(stk::topology::NODE_RANK, "coord_field");
  stk::io::set_field_role(coordField, Ioss::Field::MESH);
  m_metaDataOutPtr->set_coordinate_field(&coordField);

  stk::mesh::put_field_on_mesh(coordField, m_metaDataOutPtr->universal_part(), m_metaDataOutPtr->spatial_dimension(), nullptr);
}

void ExodusWriter::create_part(const std::string& name)
{
  m_vertPart = &(m_metaDataOutPtr->declare_part_with_topology(name + "_verts", stk::topology::NODE));
  // Setting the io attribute on the vert part causes IOSS to output a nodeset, and any vertex fields
  // will be nodeset fields rather than element fields.  This causes Paraview to display the field
  // as discrete points rather than interpolate it over elements.
  //stk::io::put_io_part_attribute(*m_vertPart);

  int numTris = mesh::count_entities_of_type(m_mesh, mesh::MeshEntityType::Triangle);
  int numQuads = mesh::count_entities_of_type(m_mesh, mesh::MeshEntityType::Quad);

  if (numTris > 0)
  {
    m_triPart  = &(m_metaDataOutPtr->declare_part_with_topology(name + "_tris", stk::topology::SHELL_TRI_3));
    stk::io::put_io_part_attribute(*m_triPart);
  }

  if (numQuads > 0)
  {
    m_quadPart = &(m_metaDataOutPtr->declare_part_with_topology(name + "_shell_quad4", stk::topology::SHELL_QUAD_4));
    stk::io::put_io_part_attribute(*m_quadPart);
  }
}

void ExodusWriter::declare_fields()
{
  m_metaDataOutPtr->enable_late_fields();

  for (auto& fieldAdaptor : m_meshFields)
  {
    StkFields stkFields;
    stk::middle_mesh::mesh::FieldShape fshape = fieldAdaptor->get_field_shape();
    if (fshape.count[0] > 0)
    {
      std::string name = fieldAdaptor->get_name() + "_verts";
      stkFields.vertField = &(m_metaDataOutPtr->declare_field<double>(stk::topology::NODE_RANK, name));
      stk::mesh::put_field_on_mesh(*(stkFields.vertField), *m_vertPart, fieldAdaptor->get_num_comp(), fshape.count[0], 0);
      stk::io::set_field_output_type((*stkFields.vertField), fieldAdaptor->get_field_output_type());
      stk::io::set_field_role(*(stkFields.vertField), Ioss::Field::TRANSIENT);
    }

    if (fshape.count[1] > 0)
    {
      throw std::runtime_error("edge fields cannot be output to Exodus");
    }

    if (fshape.count[2] > 0)
    {
      std::string name = fieldAdaptor->get_name() + "_elems";
      stkFields.elemField = &(m_metaDataOutPtr->declare_field<double>(stk::topology::ELEM_RANK, name));
      if (m_triPart)
        stk::mesh::put_field_on_mesh(*(stkFields.elemField), *m_triPart, fieldAdaptor->get_num_comp(), fshape.count[2], 0);

      if (m_quadPart)
        stk::mesh::put_field_on_mesh(*(stkFields.elemField), *m_quadPart, fieldAdaptor->get_num_comp(), fshape.count[2], 0);
      stk::io::set_field_output_type(*(stkFields.elemField), fieldAdaptor->get_field_output_type());
      stk::io::set_field_role(*(stkFields.elemField), Ioss::Field::TRANSIENT);
    }

    m_stkFields.push_back(stkFields);
  }
}

mesh::FieldPtr<stk::mesh::Entity> ExodusWriter::create_part_verts()
{
  std::vector<stk::mesh::Part*> partVec{m_vertPart};
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
  int nels = count_valid(m_mesh->get_elements());

  m_bulkDataOutPtr->modification_begin();

  std::vector<stk::mesh::Entity> entities;
  std::vector<size_t> requests(m_metaDataOutPtr->entity_rank_count(), 0);
  requests[stk::topology::ELEMENT_RANK] = nels;
  m_bulkDataOutPtr->generate_new_entities(requests, entities);

  std::vector<stk::mesh::Part*> partVec(1);
  int entityIdx = 0;
  for (auto& el : m_mesh->get_elements())
    if (el)
    {
      partVec[0] = el->get_type() == mesh::MeshEntityType::Triangle ? m_triPart : m_quadPart;
      mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
      int ndown = get_downward(el, 0, verts);
      //assert(ndown == 3);

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



void ExodusWriter::copy_fields_to_stk_mesh()
{
  for (size_t i=0; i < m_meshFields.size(); ++i)
  {
    auto fieldAdaptor = m_meshFields[i];
    StkFields stkField = m_stkFields[i];

    if (stkField.vertField)
    {
      copy_field_to_stk_mesh(PartToUse::Vert, fieldAdaptor, *(stkField.vertField));
    }

    if (stkField.elemField)
    {
      if (m_triPart)
        copy_field_to_stk_mesh(PartToUse::Triangle, fieldAdaptor, *(stkField.elemField));

      if (m_quadPart)
        copy_field_to_stk_mesh(PartToUse::Quad, fieldAdaptor, *(stkField.elemField));
    }
  }
}

void ExodusWriter::copy_field_to_stk_mesh(PartToUse partEnum, FieldOutputAdaptorPtr meshField, StkFieldType& stkField)
{
  stk::topology::rank_t entityRank;
  stk::mesh::Part* part;
  int dim = partEnum == PartToUse::Vert ? 0 : 2;
  mesh::MeshEntityType meshEntityType;

  switch (partEnum)
  {
    case PartToUse::Vert:
    {
      entityRank = stk::topology::NODE_RANK;
      part = m_vertPart;
      dim = 0;
      meshEntityType = mesh::MeshEntityType::Vertex;
      break;
    }
    case PartToUse::Triangle:
    {
      entityRank = stk::topology::ELEM_RANK;
      part = m_triPart;
      dim = 2;
      meshEntityType = mesh::MeshEntityType::Triangle;
      break;
    }  

    case PartToUse::Quad:
    {
      entityRank = stk::topology::ELEM_RANK;
      part = m_quadPart;
      dim = 2;
      meshEntityType = mesh::MeshEntityType::Quad;
      break;
    }
    default:
      throw std::runtime_error("unhandled enum");   
  }


  std::vector<double> vals(meshField->get_num_comp());
  auto& meshEntities = m_mesh->get_mesh_entities(dim);
  size_t idx = 0;
  for (stk::mesh::Bucket* bucket : m_bulkDataOutPtr->get_buckets(entityRank, *part))
    for (size_t bucketIdx = 0; bucketIdx < bucket->size(); ++bucketIdx)
    {
      while (!meshEntities[idx] || meshEntities[idx]->get_type() != meshEntityType)
        idx++;

      stk::middle_mesh::mesh::MeshEntityPtr meshEntity = meshEntities[idx];

      double* stkFieldData = stk::mesh::field_data(stkField, *bucket, bucketIdx);
      for (int i=0; i < meshField->get_field_shape().count[dim]; ++i)
      {
        meshField->get_data(meshEntity, i, vals);
        for (int j=0; j < meshField->get_num_comp(); ++j)
          stkFieldData[i * meshField->get_num_comp() + j] = vals[j];
      }

      idx++;
    }
}

void ExodusWriter::write_output_stk(const std::string& fname)
{
  stk::io::write_mesh_with_fields(fname, *m_bulkDataOutPtr, 1);
}


} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk


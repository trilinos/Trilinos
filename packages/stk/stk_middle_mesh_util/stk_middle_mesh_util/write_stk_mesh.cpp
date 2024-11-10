
#include "write_stk_mesh.hpp"
#include "stk_middle_mesh/field.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {
namespace impl {

// writes meshes defined by classes to stk Part
void StkMeshWriter::populate_part(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                                  MeshPart& pmesh1, MeshPart& pmesh2, stk::mesh::Part* part)
{
  std::vector<stk::mesh::Entity> entities;
  MNodeFieldPtr meshInVertsToStkVerts = create_part_verts(nonconformalAbstract, part);
  create_part_elements(nonconformalAbstract, meshInVertsToStkVerts, part, entities);
  write_gid_field(nonconformalAbstract, entities, pmesh1, pmesh2);
}

StkMeshWriter::MNodeFieldPtr
StkMeshWriter::create_part_verts(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                                 stk::mesh::Part* part)
{
  std::vector<stk::mesh::Part*> partVec{part};
  std::shared_ptr<mesh::Mesh> meshIn = nonconformalAbstract->get_mesh_in();
  MNodeFieldPtr stkVerts =
      mesh::create_field<stk::mesh::Entity>(meshIn, ::stk::middle_mesh::mesh::impl::FieldShape(1, 0, 0), 1);

  m_bulkDataOut.modification_begin();

  std::vector<size_t> requests(m_metaDataOut.entity_rank_count(), 0); // TODO: what is entity_rank_count
  requests[stk::topology::NODE_RANK] = count_valid(meshIn->get_vertices());
  std::vector<stk::mesh::Entity> entities;
  m_bulkDataOut.generate_new_entities(requests, entities);
  for (auto& e : entities)
    m_bulkDataOut.change_entity_parts(e, partVec);

  m_bulkDataOut.modification_end();

  // TODO: is it legal to modify a Stk Field inside a modification
  //       cycle
  const stk::mesh::FieldBase& coordField = *(m_metaDataOut.coordinate_field());
  int entityIdx                          = 0;
  for (auto& vert : meshIn->get_vertices())
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

void StkMeshWriter::create_part_elements(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                                         MNodeFieldPtr stkVerts, stk::mesh::Part* part,
                                         std::vector<stk::mesh::Entity>& entities)
{
  // create elements:
  //   get Mesh element
  //   get Mesh vertices
  //   get corresponding Stk vertices
  //   create Stk face
  //   call declare_relation(stk_face, stk_vert)

  std::shared_ptr<mesh::Mesh> meshIn = nonconformalAbstract->get_mesh_in();
  int nels                           = count_valid(meshIn->get_elements());

  m_bulkDataOut.modification_begin();

  entities.clear();
  std::vector<size_t> requests(m_metaDataOut.entity_rank_count(), 0);
  requests[stk::topology::ELEMENT_RANK] = nels;
  m_bulkDataOut.generate_new_entities(requests, entities);

  std::vector<stk::mesh::Part*> partVec{part};
  int entityIdx = 0;
  for (auto& el : meshIn->get_elements())
    if (el)
    {
      mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
      int ndown = get_downward(el, 0, verts);
      assert(ndown == 3);

      m_bulkDataOut.change_entity_parts(entities[entityIdx], partVec);
      for (int j = 0; j < ndown; ++j)
      {
        auto stkVert = (*stkVerts)(verts[j], 0, 0);
        m_bulkDataOut.declare_relation(entities[entityIdx], stkVert, j);
      }

      entityIdx++;
    }

  m_bulkDataOut.modification_end();
}

void StkMeshWriter::write_gid_field(std::shared_ptr<nonconformal::impl::NonconformalAbstract> nonconformalAbstract,
                                    std::vector<stk::mesh::Entity> entities, MeshPart& pmesh1, MeshPart& pmesh2)
{
  std::shared_ptr<mesh::Mesh> meshIn = nonconformalAbstract->get_mesh_in();
  auto& meshInElsToMesh1Els          = *(nonconformalAbstract->get_mesh1_classification());
  auto& meshInElsToMesh2Els          = *(nonconformalAbstract->get_mesh2_classification());

  int entityIdx = 0;
  for (auto& el : meshIn->get_elements())
    if (el)
    {
      auto stkEl    = entities[entityIdx];
      auto meshElL  = meshInElsToMesh1Els(el, 0, 0);
      auto meshElR  = meshInElsToMesh2Els(el, 0, 0);
      auto stkFaceL = (*pmesh1.stkEls)(meshElL, 0, 0);
      auto stkFaceR = (*pmesh2.stkEls)(meshElR, 0, 0);

      VertIdType* gidField = stk::mesh::field_data(*m_gidField, stkEl);
      // Note: stk local face ids are zero-based, but Exodus is 1-based,
      gidField[0] = m_bulkDataIn.identifier(stkFaceL.element);
      gidField[1] = stkFaceL.side + 1;
      gidField[2] = m_bulkDataIn.identifier(stkFaceR.element);
      gidField[3] = stkFaceR.side + 1;

      entityIdx++;
    }
}

} // namespace impl

} // namespace stk_interface
} // namespace middle_mesh
} // namespace stk


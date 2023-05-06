#include "element_mesh_extractor.hpp"
#include "mesh_entity.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

using predicates::impl::PointClassification;

ElementMeshData make_element_mesh_data(mesh::MeshEntityPtr el1)
{
  ElementMeshData elementMeshData;
  elementMeshData.el1                    = el1;
  elementMeshData.elementMeshIn          = mesh::make_empty_mesh(MPI_COMM_SELF);
  elementMeshData.elementMeshInVertClass = mesh::create_field<predicates::impl::PointRecord>(
      elementMeshData.elementMeshIn, mesh::FieldShape(1, 0, 0), 1);
  elementMeshData.elementMeshEntitiesToMeshInEntities =
      mesh::create_field<mesh::MeshEntityPtr>(elementMeshData.elementMeshIn, mesh::FieldShape(1, 0, 1), 1);

  return elementMeshData;
}

ElementMeshData ElementMeshExtractor::extract_element_mesh(mesh::MeshEntityPtr el1)
{
  auto& mesh1ElsToVertsIn   = *(m_relationalData->mesh1ElsToVertsIn);
  auto& vertsInClassOnMesh1 = *(m_relationalData->vertsInClassOnMesh1);

  ElementMeshData elementMeshData      = make_element_mesh_data(el1);
  auto& elementMeshVertsInClassOnMesh1 = *(elementMeshData.elementMeshInVertClass);
  auto& elementMeshVertsToMeshInVerts  = *(elementMeshData.elementMeshEntitiesToMeshInEntities);

  std::vector<EdgeVertexIdPair> edgesIn;
  get_el1_edges_in(el1, edgesIn);
  for (int i = 0; i < mesh1ElsToVertsIn.get_num_comp(el1, 0); ++i)
  {
    mesh::MeshEntityPtr vertIn = mesh1ElsToVertsIn(el1, 0, i);
    auto newVert               = elementMeshData.elementMeshIn->create_vertex(vertIn->get_point_orig(0));

    elementMeshVertsInClassOnMesh1(newVert, 0, 0) =
        convert_classification_to_element(el1, vertsInClassOnMesh1(vertIn, 0, 0));
    elementMeshVertsToMeshInVerts(newVert, 0, 0) = vertIn;
  }

  for (auto& edge : edgesIn)
  {
    mesh::MeshEntityPtr v1 = elementMeshData.elementMeshIn->get_vertices()[edge.id1];
    mesh::MeshEntityPtr v2 = elementMeshData.elementMeshIn->get_vertices()[edge.id2];

    elementMeshData.elementMeshIn->create_edge(v1, v2);
  }

  return elementMeshData;
}

void ElementMeshExtractor::write_elements_back_to_middle_grid(ElementMeshData& elementMeshData, mesh::MeshEntityPtr el1)
{
  auto& elementMeshEntitiesToMeshInEntities = *(elementMeshData.elementMeshEntitiesToMeshInEntities);
  auto& meshInElementsToMesh1Elements       = *(m_relationalData->meshInElementsToMesh1Elements);

  for (auto& el : elementMeshData.elementMeshIn->get_elements())
    if (el)
    {
      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
      [[maybe_unused]] int nverts  = get_downward(el, 0, verts.data());
      mesh::MeshEntityPtr meshInEl = nullptr;
      if (el->get_type() == mesh::MeshEntityType::Triangle)
      {
        assert(nverts == 3);
        meshInEl = m_meshIn->create_triangle_from_verts(elementMeshEntitiesToMeshInEntities(verts[0], 0, 0),
                                                        elementMeshEntitiesToMeshInEntities(verts[1], 0, 0),
                                                        elementMeshEntitiesToMeshInEntities(verts[2], 0, 0));
      } else if (el->get_type() == mesh::MeshEntityType::Quad)
      {
        assert(nverts == 4);
        meshInEl = m_meshIn->create_quad_from_verts(
            elementMeshEntitiesToMeshInEntities(verts[0], 0, 0), elementMeshEntitiesToMeshInEntities(verts[1], 0, 0),
            elementMeshEntitiesToMeshInEntities(verts[2], 0, 0), elementMeshEntitiesToMeshInEntities(verts[3], 0, 0));
      } else
        throw std::runtime_error("unrecognized element type");

      meshInElementsToMesh1Elements(meshInEl, 0, 0) = el1;
      elementMeshEntitiesToMeshInEntities(el, 0, 0) = meshInEl;
    }
}

void ElementMeshExtractor::get_el1_edges_in(mesh::MeshEntityPtr el1, std::vector<EdgeVertexIdPair>& edgesIn)
{
  auto& mesh1ElsToVertsIn = *(m_relationalData->mesh1ElsToVertsIn);
  edgesIn.clear();
  sort_mesh_in_verts(el1);

  for (int i = 0; i < mesh1ElsToVertsIn.get_num_comp(el1, 0); ++i)
  {
    mesh::MeshEntityPtr vertIn = mesh1ElsToVertsIn(el1, 0, i);
    for (int j = 0; j < vertIn->count_up(); ++j)
    {
      mesh::MeshEntityPtr edge      = vertIn->get_up(j);
      mesh::MeshEntityPtr otherVert = get_other_vert(vertIn, edge);

      int idx = get_vert_in_mesh1_el_index(el1, otherVert);
      if (idx != -1)
        edgesIn.push_back(EdgeVertexIdPair{i, idx});
    }
  }

  std::sort(edgesIn.begin(), edgesIn.end());
  auto it = std::unique(edgesIn.begin(), edgesIn.end());
  edgesIn.erase(it, edgesIn.end());
}

void ElementMeshExtractor::sort_mesh_in_verts(mesh::MeshEntityPtr el1)
{
  auto& mesh1ElsToVertsIn         = *(m_relationalData->mesh1ElsToVertsIn);
  int nverts                      = mesh1ElsToVertsIn.get_num_comp(el1, 0);
  mesh::MeshEntityPtr* vertsStart = &(mesh1ElsToVertsIn(el1, 0, 0));
  mesh::MeshEntityPtr* vertsEnd   = &(mesh1ElsToVertsIn(el1, 0, nverts - 1)) + 1;
  std::sort(vertsStart, vertsEnd, mesh::is_less);

  assert(std::unique(vertsStart, vertsEnd) == vertsEnd);
}

int ElementMeshExtractor::get_vert_in_mesh1_el_index(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr vertIn)
{
  auto& mesh1ElsToVertsIn = *(m_relationalData->mesh1ElsToVertsIn);
  int nverts              = mesh1ElsToVertsIn.get_num_comp(el1, 0);
  assert(nverts > 0);

  mesh::MeshEntityPtr* vertsStart = &(mesh1ElsToVertsIn(el1, 0, 0));
  mesh::MeshEntityPtr* vertsEnd   = &(mesh1ElsToVertsIn(el1, 0, nverts - 1)) + 1;

  auto it = std::lower_bound(vertsStart, vertsEnd, vertIn, mesh::is_less);
  if (it == vertsEnd || *it != vertIn)
    return -1;
  else
    return std::distance(vertsStart, it);
}

predicates::impl::PointRecord
ElementMeshExtractor::convert_classification_to_element(mesh::MeshEntityPtr el1,
                                                        const predicates::impl::PointRecord& record)
{
  assert(record.type != PointClassification::Exterior);

  if (record.type == PointClassification::Interior)
  {
    assert(record.el == el1);
    return record;
  } else
  {
    mesh::MeshEntityPtr entity = get_entity(record);
    int newLocalId             = predicates::impl::get_entity_id(el1, entity);
    return predicates::impl::PointRecord(record.type, newLocalId, el1);
  }
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk

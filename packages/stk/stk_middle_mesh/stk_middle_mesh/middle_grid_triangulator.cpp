#include "middle_grid_triangulator.hpp"

#include "element_mesh_classifier.hpp"
#include "element_mesh_extractor.hpp"
#include "element_mesh_triangulator.hpp"
#include "mesh_io.hpp"
#include "adjacency_search.hpp"
#include "utils.hpp"
#include "bounding_box_search.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

void MiddleGridTriangulator::triangulate()
{
  ElementMeshExtractor extractor(m_meshIn, m_relationalData);
  ElementMeshTriangulator triangulator(m_relationalData);
  ElementMeshClassifier elementClassifier(m_relationalData, m_classifier);
  std::vector<mesh::MeshEntityPtr> mesh2Els;
  mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> mesh1ElsToMesh2ElsPtr = compute_mesh1_to_mesh2_element_maps();
  auto& mesh1ElsToMesh2Els = *mesh1ElsToMesh2ElsPtr;

  for (auto el1 : m_mesh1->get_elements())
  {
    if (el1)
    {
      mesh2Els.assign(mesh1ElsToMesh2Els(el1, 0).begin(), mesh1ElsToMesh2Els(el1, 0).end());

      if (m_output)
      {
        std::cout << "\nTriangulating el1 = " << el1 << ", id = " << el1->get_id() << std::endl;
        for (int i = 0; i < el1->count_down(); ++i)
          std::cout << "el1 edge id = " << el1->get_down(i)->get_id() << std::endl;
      }

      ElementMeshData elementMeshData = extractor.extract_element_mesh(el1);
      int numConstraintEdges          = count_valid(elementMeshData.elementMeshIn->get_edges());

      if (m_output)
      {
        mesh::impl::print_vert_edges("mesh_in_element", elementMeshData.elementMeshIn);
      }

      triangulator.triangulate(elementMeshData);

      extractor.write_elements_back_to_middle_grid(elementMeshData, el1);

      elementClassifier.classify(elementMeshData, numConstraintEdges, mesh2Els);
    }
  }
}

mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
MiddleGridTriangulator::compute_mesh1_to_mesh2_element_maps()
{
  using SearchMesh = mesh::impl::SearchMeshElementBoundingBoxBase;
  using CoarseSearch = search::ElementToElementBoundingBoxSearch;
  auto mesh1ElsToMesh2Els = mesh::create_variable_size_field<mesh::MeshEntityPtr>(m_mesh1, mesh::FieldShape(0, 0, 1));


  predicates::impl::AveragedNormalField averagedNormalField(m_mesh2);
  auto searchMesh1 = std::make_shared<mesh::impl::SearchMeshElementBoundingBox>(m_mesh1, MPI_COMM_SELF);
  auto searchMesh2 = std::make_shared<mesh::impl::SearchMeshElementBoundingBoxNormal>(m_mesh2, MPI_COMM_SELF,
                                        averagedNormalField.get_field(), m_searchOpts.normalDirectionFactor);

  bool doParallelSearch = false;
  CoarseSearch coarseSearch(searchMesh1, searchMesh2, "MiddleGridTriangulator::local_coarse_search", MPI_COMM_SELF, doParallelSearch, m_searchOpts);
  coarseSearch.coarse_search();

  if (coarseSearch.get_unpaired_recv_entities().size() > 0)
  {
    throw std::runtime_error("local coarse search could not find candidate elements for some vertices");
  }

  const CoarseSearch::EntityProcRelationVec& mesh2To1Relations = coarseSearch.get_range_to_domain();
  for (const CoarseSearch::EntityProcRelation& relation : mesh2To1Relations)
  {
    SearchMesh::EntityKey mesh2ElId = relation.first.id();
    SearchMesh::EntityKey mesh1ElId = relation.second.id();
    mesh::MeshEntityPtr mesh1El = m_mesh1->get_elements()[mesh1ElId];
    mesh::MeshEntityPtr mesh2El = m_mesh2->get_elements()[mesh2ElId];
    mesh1ElsToMesh2Els->insert(mesh1El, 0, mesh2El);
  }

  return mesh1ElsToMesh2Els;
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk

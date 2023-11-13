#include "middle_grid_triangulator.hpp"

#include "element_mesh_classifier.hpp"
#include "element_mesh_extractor.hpp"
#include "element_mesh_triangulator.hpp"
#include "mesh_io.hpp"
#include "adjacency_search.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

void MiddleGridTriangulator::triangulate()
{
  ElementMeshExtractor extractor(m_meshIn, m_relationalData);
  ElementMeshTriangulator triangulator(m_relationalData);
  ElementMeshClassifier elementClassifier(m_relationalData, m_classifier);
  mesh::impl::AdjacencySearch adjSearch(m_mesh1, m_mesh2);
  mesh::MeshEntityPtr el1;
  std::vector<mesh::MeshEntityPtr> mesh2Els;
  while ( (el1 = adjSearch.get_next(mesh2Els)) )
  {
    if (el1)
    {
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

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk

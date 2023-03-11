
#include "create_mesh.hpp"
#include "field.hpp"
#include "variable_size_field.hpp"

namespace stk {
namespace middle_mesh {
namespace test_util {

void test_area_per_element(std::shared_ptr<mesh::Mesh> mesh,
                           mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inverseClassificationPtr);

void test_areas_positive(std::shared_ptr<mesh::Mesh> mesh);

double compute_mesh_area(std::shared_ptr<mesh::Mesh> mesh);

void test_total_areas_same(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                           std::shared_ptr<mesh::Mesh> meshIn);

void test_number_of_elements(std::shared_ptr<mesh::Mesh> mesh, int numel);

void test_every_element_classified(std::shared_ptr<mesh::Mesh> meshIn,
                                   mesh::FieldPtr<mesh::MeshEntityPtr> elementClassificationPtr);

void test_every_element_classified_inverse(std::shared_ptr<mesh::Mesh> mesh,
                                           mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inverseClassificationPtr);
} // namespace test_util
} // namespace middle_mesh
} // namespace stk

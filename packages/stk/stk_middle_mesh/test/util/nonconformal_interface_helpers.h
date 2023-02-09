
#include "create_mesh.h"
#include "field.h"
#include "variable_size_field.h"

namespace stk {
namespace middle_mesh {
namespace impl {

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
} // namespace impl
} // namespace middle_mesh
} // namespace stk

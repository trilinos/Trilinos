#include "util/nonconformal_interface_helpers.hpp"
#include "stk_middle_mesh/element_operations_2d.hpp"
#include "stk_middle_mesh/field.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace test_util {

void test_area_per_element(std::shared_ptr<mesh::Mesh> mesh,
                           mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inverseClassificationPtr)
{
  mesh::impl::ElementOperations2D elemOps;
  auto& inverseClassification = *inverseClassificationPtr;

  for (auto& el : mesh->get_elements())
    if (el)
    {
      double elArea         = elemOps.compute_area(el);
      double middleGridArea = 0;
      for (int i = 0; i < inverseClassification.get_num_comp(el, 0); ++i)
      {
        mesh::MeshEntityPtr elIn = inverseClassification(el, 0, i);
        middleGridArea += elemOps.compute_area(elIn);
      }

      EXPECT_NEAR(elArea, middleGridArea, 1e-13);
    }
}

void test_areas_positive(std::shared_ptr<mesh::Mesh> mesh)
{
  mesh::impl::ElementOperations2D elemOps;
  for (auto& el : mesh->get_elements())
    if (el)
    {
      EXPECT_GT(elemOps.compute_area(el), 0);
    }
}




double compute_mesh_area(std::shared_ptr<mesh::Mesh> mesh, MPI_Comm unionComm)
{
  double area = 0;
  double areaLocal = 0;

  if (mesh)
  {
    mesh::impl::ElementOperations2D elemOps;
    areaLocal = elemOps.compute_area(mesh);
  }

  MPI_Allreduce(&areaLocal, &area, 1, MPI_DOUBLE, MPI_SUM, unionComm);

  return area;
}

double compute_mesh_area(std::shared_ptr<mesh::Mesh> mesh)
{
  return compute_mesh_area(mesh, mesh->get_comm());
}

void test_total_areas_same(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                           std::shared_ptr<mesh::Mesh> meshIn, MPI_Comm unionComm)
{
  double area1 = 0, area2 = 2, areaIn = 0;
  area1  = compute_mesh_area(mesh1, unionComm);
  area2  = compute_mesh_area(mesh2, unionComm);
  areaIn = compute_mesh_area(meshIn, unionComm);


  EXPECT_NEAR(area1, area2, 1e-13);
  EXPECT_NEAR(areaIn, area1, 1e-13);
  EXPECT_NEAR(areaIn, area2, 1e-13);
}

void test_number_of_elements(std::shared_ptr<mesh::Mesh> mesh, int numel)
{
  EXPECT_EQ(count_valid(mesh->get_elements()), numel);
}

void test_every_element_classified(std::shared_ptr<mesh::Mesh> meshIn,
                                   mesh::FieldPtr<mesh::MeshEntityPtr> elementClassificationPtr)
{
  auto& elementClassification = *elementClassificationPtr;
  for (auto& el : meshIn->get_elements())
    if (el)
    {
      EXPECT_NE(elementClassification(el, 0, 0), nullptr);
    }
}

void test_every_element_classified_inverse(std::shared_ptr<mesh::Mesh> mesh,
                                           mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> inverseClassificationPtr)
{
  auto& inverseClassification = *inverseClassificationPtr;
  for (auto& el : mesh->get_elements())
    if (el)
    {
      EXPECT_GT(inverseClassification.get_num_comp(el, 0), 0);
    }
}

} // namespace test_util
} // namespace middle_mesh
} // namespace stk

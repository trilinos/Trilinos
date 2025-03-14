#include "gtest/gtest.h"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_io.hpp"
#include "stk_middle_mesh/middle_mesh_point_projection.hpp"
#include "stk_middle_mesh/nonconformal4.hpp"
#include "stk_middle_mesh/create_mesh.hpp"

namespace {

using namespace stk::middle_mesh;

class MiddleMeshPointProjectionTester
{
  public:
    MiddleMeshPointProjectionTester(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2) :
      m_mesh1(mesh1),
      m_mesh2(mesh2)
    {}

    void runtest(std::shared_ptr<XiCoordinates> xiCoords)
    {
      stk::middle_mesh::nonconformal4::impl::Nonconformal4 interface(m_mesh1, m_mesh2);
      auto meshIn          = interface.create();
      auto mesh1Class      = interface.get_mesh1_classification();
      auto mesh2Class      = interface.get_mesh2_classification();
      auto xiCoordsOnMesh1 = interface.compute_points_on_mesh1(xiCoords);
      auto xiCoordsOnMesh2 = interface.compute_points_on_mesh2(xiCoords);

      test_coordinate_interpolation(meshIn, xiCoords, m_mesh1, xiCoordsOnMesh1, mesh1Class);

      test_coordinate_interpolation(meshIn, xiCoords, m_mesh2, xiCoordsOnMesh2, mesh2Class);
    }


  private:
    void test_coordinate_interpolation(std::shared_ptr<mesh::Mesh> meshIn, std::shared_ptr<XiCoordinates> xiCoords,
                                       std::shared_ptr<mesh::Mesh> /*inputMesh*/, mesh::FieldPtr<utils::Point> xiCoordsOnInputMeshPtr,
                                       mesh::FieldPtr<mesh::MeshEntityPtr> meshInToInputMeshPtr)
    {
      auto& xiCoordsOnInputMesh = *xiCoordsOnInputMeshPtr;
      auto& meshInToInputMesh   = *meshInToInputMeshPtr;
      for (auto& elIn : meshIn->get_elements())
        if (elIn)
        {
          mesh::MeshEntityPtr elInput = meshInToInputMesh(elIn, 0, 0);
          auto& xiCoordsEl = xiCoords->get_xi_coords(elIn->get_type());
          std::pair<double, double> xiRangeIn    = xiCoords->get_xi_coord_range(elIn->get_type());
          std::pair<double, double> xiRangeInput = xiCoords->get_xi_coord_range(elInput->get_type());

          for (size_t i=0; i < xiCoordsEl.size(); ++i)
          {
            utils::Point xiCoordsIn = mesh::convert_xi_coords_from_range(xiRangeIn.first, xiRangeIn.second, xiCoordsEl[i]);
            auto middleMeshCoords = mesh::compute_coords_from_xi_3d(elIn, xiCoordsIn);

            utils::Point xiCoordsInput = mesh::convert_xi_coords_from_range(xiRangeInput.first, xiRangeInput.second, xiCoordsOnInputMesh(elIn, i, 0));
            auto inputMeshCoords  = mesh::compute_coords_from_xi_3d(elInput, xiCoordsInput);

            for (int d=0; d < 3; ++d)
              EXPECT_NEAR(middleMeshCoords[d], inputMeshCoords[d], 1e-13);
          }

        }
    }

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
};

class XiCoordinatesForTest : public XiCoordinates
{
  public:
    const std::vector<utils::Point>& get_xi_coords(mesh::MeshEntityType type) override
    {
      if (type == mesh::MeshEntityType::Quad)
        return m_quadXiCoords;
      else if (type == mesh::MeshEntityType::Triangle)
        return m_triangleXiCoords;
      else
        throw std::runtime_error("only triangles and quads supported");
    }

    std::pair<double, double> get_xi_coord_range(mesh::MeshEntityType /*type*/) override { return std::make_pair(0, 2); }

  private:
    std::vector<utils::Point> m_triangleXiCoords = {utils::Point(0.1, 0.2), utils::Point(1.5, 0.3)};
    std::vector<utils::Point> m_quadXiCoords = {utils::Point(0.1, 0.2), utils::Point(0.6, 0.7), utils::Point(1.4, 0.4)};
};


}

TEST(MiddleMeshPointProjection, Quads)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec1, spec2;

  spec1.xmin = 0;   spec2.xmin = 0;
  spec1.xmax = 1;   spec2.xmax = 1;
  spec1.ymin = 0;   spec2.ymin = 0;
  spec1.ymax = 1;   spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 3;
  spec1.numelY = 2; spec2.numelY = 3;

  auto f = [](const utils::Point& pt) { return pt; };

  auto mesh1 = mesh::impl::create_mesh(spec1, f);
  auto mesh2 = mesh::impl::create_mesh(spec2, f);
  auto xiCoords = std::make_shared<XiCoordinatesForTest>();

  MiddleMeshPointProjectionTester tester(mesh1, mesh2);
  tester.runtest(xiCoords);
}

TEST(MiddleMeshPointProjection, QuadsSkewed)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  mesh::impl::MeshSpec spec1, spec2;

  spec1.xmin = 0;   spec2.xmin = 0;
  spec1.xmax = 1;   spec2.xmax = 1;
  spec1.ymin = 0;   spec2.ymin = 0;
  spec1.ymax = 1;   spec2.ymax = 1;
  spec1.numelX = 2; spec2.numelX = 3;
  spec1.numelY = 2; spec2.numelY = 3;

  auto f = [](const utils::Point& pt) { return utils::Point(pt.x, pt.y + 0.5*pt.x*pt.y); };

  auto mesh1 = mesh::impl::create_mesh(spec1, f);
  auto mesh2 = mesh::impl::create_mesh(spec2, f);
  auto xiCoords = std::make_shared<XiCoordinatesForTest>();

  MiddleMeshPointProjectionTester tester(mesh1, mesh2);
  tester.runtest(xiCoords);
}
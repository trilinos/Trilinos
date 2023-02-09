#include "create_mesh.h"
#include "mesh.h"
#include "mesh_io.h"
#include "predicates/point_classifier_normal_wrapper.h"
#include "gtest/gtest.h"

#include "stk_search/Box.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

#if 0

class MeshScatterSpec {
  void addDestination(MeshEntityPtr entity, int dest_rank) {};

  void getDestinations(MeshEntityPtr entity, std::vector<int>& dest_ranks) {};
};

class TransferMesh {
 public:
  using Entity = MeshEntity;
  using EntityKey = int;

  using EntityProc = std::pair<EntityKey, int>;
  using EntityProcVec = std::vector<EntityProc>;
  // using Point = stk::search::Point<double>;
  using Box = stk::search::Box<double>;
  using BoundingBox = std::pair<Box, EntityProc>;

  TransferMesh(std::shared_ptr<Mesh> inputMesh, MeshEntityType inputTransferType)
    : mesh(inputMesh),
      meshEntityType(inputTransferType)
    {}

  void create_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const {};

 private:
  std::shared_ptr<Mesh> mesh;
  MeshEntityType meshEntityType;
};

TEST(ParallelSearch, CreateBoundingBox) {
  MeshSpec spec, spec2;
  spec.numel_x = 5;
  spec.numel_y = 5;
  spec.xmin = 0;
  spec.xmax = 1;
  spec.ymin = 0;
  spec.ymax = 1;

  spec2.numel_x = 7;
  spec2.numel_y = 7;
  spec2.xmin = 0;
  spec2.xmax = 1;
  spec2.ymin = 0;
  spec2.ymax = 1;

  auto func = [&](const Point& pt) { return Point(pt.x, pt.y, 0); };

  std::shared_ptr<Mesh> mesh;
  std::shared_ptr<PointClassifierNormalWrapper> point_classifier = nullptr;


  // std::vector<

  // TransferMesh meshA(mesh, MeshEntityType::Vertex);
  // meshA.create_bounding_boxes(
}

TEST(ParallelSearch, CheckBoundingBoxProp) {

}

TEST(ParallelSearch, CheckBoundingBoxInflation) {

}

TEST(ParallelSearch, CreateSendMesh) {

}

TEST(ParallelSearch, CheckSendMeshProp) {

}

TEST(ParallelSearch, CreateRecvMesh) {

}

TEST(ParallelSearch, CheckRecvMeshProp) {

}

TEST(ParallelSearch, CompareTwoEntities) {

}

TEST(ParallelSearch, FiveToSeven)
{
  int nprocs_mesh1 = 2, nprocs_mesh2 = 3;
  int myrank_world, commsize_world;
  MPI_Comm_size(MPI_COMM_WORLD, &commsize_world);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_world);
  int color = myrank_world < nprocs_mesh1 ? 0 : 1;

  if (commsize_world != (nprocs_mesh1 + nprocs_mesh2))
    GTEST_SKIP();

  MPI_Comm split_comm;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &split_comm);

  MeshSpec spec, spec2;
  spec.numel_x = 5;
  spec.numel_y = 5;
  spec.xmin = 0;
  spec.xmax = 1;
  spec.ymin = 0;
  spec.ymax = 1;

  spec2.numel_x = 7;
  spec2.numel_y = 7;
  spec2.xmin = 0;
  spec2.xmax = 1;
  spec2.ymin = 0;
  spec2.ymax = 1;

  auto func = [&](const Point& pt) { return Point(pt.x, pt.y, 0); };

  std::shared_ptr<Mesh> mesh;
  std::shared_ptr<PointClassifierNormalWrapper> point_classifier = nullptr;
  if (color == 0)
  {
    //std::cout << "creating mesh1" << std::endl;
    mesh = createMesh(spec, func, false, split_comm);
  } else
  {
    //std::cout << "\ncreating mesh2" << std::endl;
    mesh = createMesh(spec2, func, false, split_comm);
    point_classifier = std::make_shared<PointClassifierNormalWrapper>(mesh);  // point_classifier owns the normal vector field
                                                                              // used to define the bounding boxes
  }


  //MeshScatterSpec spec = /* do parallel search here: send color == 1 mesh to color == 0 mesh */


  // test that for each element on color 0 mesh, every point inside that element
  // is inside a color 1 mesh element that is sent the owning process of the color 0 element


}

#endif
} // namespace impl
} // namespace middle_mesh
} // namespace stk

#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_middle_mesh_util/create_stk_mesh.hpp>
#include <stk_middle_mesh/incremental_mesh_boundary_snapper.hpp>
#include <stk_middle_mesh/nonconformal4.hpp>
#include <stk_middle_mesh/application_interface.hpp>
#include <stk_middle_mesh/create_mesh.hpp>
#include <stk_middle_mesh/utils.hpp>
#include <stk_coupling/SplitComms.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/timer.hpp>
#include <memory>
#include <iostream>

namespace
{

stk::middle_mesh::stk_interface::MeshPart
create_surface_mesh(stk::ParallelMachine comm,
                    const std::string& fileName,
                    const std::string& surfacePartName)
{
  stk::middle_mesh::stk_interface::StkMeshCreator stkMeshCreator(fileName, "NONE", comm);
  stk::mesh::Part* surfacePart = stkMeshCreator.get_meta_data().get_part(surfacePartName);
  EXPECT_TRUE((surfacePart != nullptr));
  stk::mesh::Selector ownedSurface = stkMeshCreator.get_meta_data().locally_owned_part() & *surfacePart;
  const unsigned stkNumFaces = stk::mesh::count_entities(stkMeshCreator.get_bulk_data(), stk::topology::FACE_RANK, ownedSurface);
  stk::middle_mesh::stk_interface::MeshPart meshAndField = stkMeshCreator.create_mesh_from_part(surfacePartName);
  const int myRank = stk::parallel_machine_rank(comm);

  std::cout<<"localRank="<<myRank<<", mesh="<<fileName
           <<", stk surface numFaces: " << stkNumFaces
           <<", surface-mesh size: "<< meshAndField.mesh->get_elements().size()
           <<std::endl;
  return meshAndField;
}

TEST(MiddleMesh, StkMeshCreator)
{
  stk::ParallelMachine commWorld = MPI_COMM_WORLD;

  std::string meshFileName = stk::unit_test_util::get_command_line_option<std::string>("--mesh", "generated:2x1000x1000|sideset:x");
  std::string surfacePartName = stk::unit_test_util::get_command_line_option<std::string>("--surface", "surface_1");

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 1;
  stk::unit_test_util::BatchTimer batchTimer(commWorld);

  batchTimer.initialize_batch_timer();
  for(unsigned run=0; run<NUM_RUNS; ++run) {
    batchTimer.start_batch_timer();

    create_surface_mesh(commWorld, meshFileName, surfacePartName);

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMesh, CreateMeshUsingMeshGenerator)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  const int NUM_RUNS = 5;
  const int NUM_ITERS = 1;
  stk::middle_mesh::mesh::impl::MeshSpec spec{2000, 2000, 0, 1, 0, 1};

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();
    stk::middle_mesh::mesh::impl::create_mesh(spec, [&](stk::middle_mesh::utils::Point const& pt) { return pt; } );
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}


TEST(MiddleMesh, CreateMiddleMeshFromIdenticalMeshes)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  const int NUM_RUNS = 5;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{600, 600, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();

    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);

    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::Parallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMesh, CreateMiddleMeshFromOffByOneMeshes)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  const int NUM_RUNS = 5;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{500, 500, 0, 1, 0, 1};
  stk::middle_mesh::mesh::impl::MeshSpec spec2{501, 501, 0, 1, 0, 1};

  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();
    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);


    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::Parallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMesh, CreateMiddleMeshFromMeshesWithLargeRefinementRatio)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  const int NUM_RUNS = 5;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{40, 40, 0, 1, 0, 1};
  stk::middle_mesh::mesh::impl::MeshSpec spec2{800, 800, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();

    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);

    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::Parallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMesh, CreateMiddleMeshFromMeshesWithLargeRefinementRatio2)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  stk::unit_test_util::BatchTimer batchTimer(comm);
  const int NUM_RUNS = 5;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{800, 800, 0, 1, 0, 1};
  stk::middle_mesh::mesh::impl::MeshSpec spec2{40, 40, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();

    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);

    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::Parallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);
}


namespace {
class CommSplitter
{
  public:
    CommSplitter(MPI_Comm comm) :
      m_inputComm(comm)
    {
      split_comm();
    }

    ~CommSplitter() { MPI_Comm_free(&m_splitComm); }

    int get_color() const { return m_color; }

    MPI_Comm get_comm() const { return m_splitComm; }

  private:

    void split_comm()
    {
      int inputCommSize = stk::middle_mesh::utils::impl::comm_size(m_inputComm);
      MPI_Comm newComm;
      int color = 0;
      if (inputCommSize == 1)
      {
        MPI_Comm_dup(m_inputComm, &newComm);
      } else
      {
        int nprocs1 = inputCommSize / 2 + (inputCommSize % 2);
        //int nprocs2 = inputCommSize / 2;

        color = stk::middle_mesh::utils::impl::comm_rank(m_inputComm) < nprocs1  ? 0 : 1;

        MPI_Comm_split(m_inputComm, color, 0, &newComm);
      }

      m_color = color;
      m_splitComm = newComm;
    }

      MPI_Comm m_inputComm;
      int m_color;
      MPI_Comm m_splitComm;

};
}

TEST(MiddleMesh, CreateMiddleMeshMPMD)
{
  using namespace stk::middle_mesh;
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  int commSize = utils::impl::comm_size(comm);

  CommSplitter commSplitter(comm);
  stk::unit_test_util::BatchTimer batchTimer(comm);
  const int NUM_RUNS = 5;
  const int NUM_ITERS = 1;

  stk::middle_mesh::mesh::impl::MeshSpec spec{400, 400, 0, 1, 0, 1};
  stk::middle_mesh::mesh::impl::MeshSpec spec2{401, 401, 0, 1, 0, 1};
  auto func = [&](const utils::Point& pt) { return utils::Point(pt.x, 0, pt.y); };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();  

    std::shared_ptr<mesh::Mesh> mesh1, mesh2;

    if (commSize == 1 || commSplitter.get_color() == 0)
      mesh1 = create_mesh(spec, func, commSplitter.get_comm());

    if (commSize == 1 || commSplitter.get_color() == 1)    
      mesh2 = create_mesh(spec2, func, commSplitter.get_comm());

    double eps = 1e-12;
    BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.incrementalMeshBoundarySnapOpts.qualityImproverOpts.nlayers=5;

    MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = impl::EdgeTracerTolerances(eps);
    auto interface = application_interface_factory(ApplicationInterfaceType::Parallel, mesh1, mesh2,
                                                    MPI_COMM_WORLD, nullptr, ParallelSearchOpts(), 
                                                    VolumeSnapOpts(),
                                                    snapOpts,
                                                    middleMeshOpts);

    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(NUM_ITERS);  
}



}


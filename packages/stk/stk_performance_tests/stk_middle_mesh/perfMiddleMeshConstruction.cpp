#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_middle_mesh_util/create_stk_mesh.hpp>
#include <stk_middle_mesh/incremental_mesh_boundary_snapper.hpp>
#include <stk_middle_mesh/nonconformal4.hpp>
#include <stk_middle_mesh/application_interface.hpp>
#include <stk_middle_mesh/create_mesh.hpp>
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
  stk::middle_mesh::stk_interface::StkMeshCreator stkMeshCreator(fileName, comm);
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

TEST(MiddleMesh, creationPerformance_StkMeshCreator_MPMD)
{
  stk::ParallelMachine commWorld = MPI_COMM_WORLD;

  const int localRank = stk::parallel_machine_rank(commWorld);
  int myColor = stk::unit_test_util::simple_fields::get_command_line_option<int>("--color", 0);

  stk::coupling::SplitComms splitComms(commWorld, myColor);
  const int numColors = 1 + splitComms.get_other_colors().size();
  if (numColors == 1) { GTEST_SKIP(); }
  ASSERT_EQ(2, numColors)<<"This MPMD test only supports 2 colors";

  const bool localMeshIsMesh1 = myColor < splitComms.get_other_colors()[0];

  std::string meshFileName = stk::unit_test_util::simple_fields::get_command_line_option<std::string>("--mesh", "generated:2x8x8|sideset:x");
  std::string surfacePartName = stk::unit_test_util::simple_fields::get_command_line_option<std::string>("--surface", "surface_1");

  const unsigned NUM_RUNS = 5;
  const int NUM_ITERS = 1;
  stk::unit_test_util::BatchTimer batchTimer(commWorld);

  stk::middle_mesh::stk_interface::MeshPart meshAndField = create_surface_mesh(splitComms.get_split_comm(), meshFileName, surfacePartName);

  batchTimer.initialize_batch_timer();
  for(unsigned run=0; run<NUM_RUNS; ++run) {
    batchTimer.start_batch_timer();

    std::shared_ptr<stk::middle_mesh::mesh::Mesh> surfaceMesh1;
    std::shared_ptr<stk::middle_mesh::mesh::Mesh> surfaceMesh2;

    if (localMeshIsMesh1) {
      surfaceMesh1 = meshAndField.mesh;
    }
    else {
      surfaceMesh2 = meshAndField.mesh;
    }

    stk::middle_mesh::ApplicationInterfaceType appInterfaceType = stk::middle_mesh::ApplicationInterfaceType::FakeParallel;

    std::shared_ptr<stk::middle_mesh::ApplicationInterface> appInterface =
      stk::middle_mesh::application_interface_factory(appInterfaceType, surfaceMesh1, surfaceMesh2, commWorld);

    appInterface->create_middle_grid();

    const std::string middleMeshSizeMsg = localMeshIsMesh1 ?
      ", mesh1 middle-mesh size: "+std::to_string(appInterface->get_middle_grid_for_mesh1()->get_elements().size())
      :
      ", mesh2 middle-mesh size: "+std::to_string(appInterface->get_middle_grid_for_mesh2()->get_elements().size());

    if (stk::parallel_machine_rank(splitComms.get_split_comm()) == 0) {
      std::cout << "iter " << run <<": commWorldRank "<<localRank
                << ", color-"<<myColor<<" rank "<<stk::parallel_machine_rank(splitComms.get_split_comm())
                << middleMeshSizeMsg
                << std::endl;
    }

    batchTimer.stop_batch_timer();
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

TEST(MiddleMesh, CreateMeshUsingMeshGenerator)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 100);
  stk::middle_mesh::mesh::impl::MeshSpec spec{150, 200, 0, 1, 0, 1};

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    batchTimer.start_batch_timer();
    stk::middle_mesh::mesh::impl::create_mesh(spec, [&](stk::middle_mesh::utils::Point const& pt) { return pt; } );
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(numRuns);
}

TEST(MiddleMesh, CreateMeshFromSTKMesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 100);
  std::string meshSpec = stk::unit_test_util::simple_fields::get_command_line_option<std::string>("-s", "generated:200x50x5|sideset:x");

  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  stk::io::fill_mesh(meshSpec, *bulkDataPtr);

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    stk::middle_mesh::stk_interface::StkMeshCreator stkMeshCreator(bulkDataPtr);

    const std::string surfacePartName = "surface_1";
    stk::mesh::Part* surfacePart = stkMeshCreator.get_meta_data().get_part(surfacePartName);
    ASSERT_TRUE((surfacePart != nullptr));

    stk::mesh::Selector ownedSurface = stkMeshCreator.get_meta_data().locally_owned_part() & *surfacePart;

    batchTimer.start_batch_timer();
    auto meshAndField = stkMeshCreator.create_mesh_from_part(surfacePartName);
    batchTimer.stop_batch_timer();

    std::shared_ptr<stk::middle_mesh::mesh::Mesh> surfaceMesh = meshAndField.mesh;
    stk::middle_mesh::mesh::check_topology(surfaceMesh);
  }
  batchTimer.print_batch_timing(numRuns);
}

TEST(MiddleMesh, CreateMiddleMeshFromIdenticalMeshes)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 50);

  stk::middle_mesh::mesh::impl::MeshSpec spec{5, 5, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    batchTimer.start_batch_timer();

    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);

    stk::middle_mesh::mesh::impl::ElementOperations2D elemOps;

    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::FakeParallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(numRuns);
}

TEST(MiddleMesh, CreateMiddleMeshFromOffByNMeshes)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 50);
  int numMeshes = stk::unit_test_util::simple_fields::get_command_line_option<int>("-m", 50);

  stk::middle_mesh::mesh::impl::MeshSpec spec{5, 5, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    batchTimer.start_batch_timer();
    for (int j = 1; j < numMeshes; ++j) {
      stk::middle_mesh::mesh::impl::MeshSpec spec2{5+i, 5+i, 0, 1, 0, 1};
      auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
      auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);

      stk::middle_mesh::mesh::impl::ElementOperations2D elemOps;

      double eps = 1e-12;
      stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
      snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

      stk::middle_mesh::MiddleGridOpts middleMeshOpts;
      middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
      middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
      auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::FakeParallel,
                                                                      mesh1, mesh2, comm, nullptr,
                                                                      stk::middle_mesh::ParallelSearchOpts(),
                                                                      stk::middle_mesh::VolumeSnapOpts(),
                                                                      snapOpts, middleMeshOpts);
      interface->create_middle_grid();
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(numRuns);
}

TEST(MiddleMesh, CreateMiddleMeshFromMeshesWithLargeRefinementRatio)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 100);

  stk::middle_mesh::mesh::impl::MeshSpec spec{10, 10, 0, 1, 0, 1};
  stk::middle_mesh::mesh::impl::MeshSpec spec2{100, 100, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    batchTimer.start_batch_timer();

    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);

    stk::middle_mesh::mesh::impl::ElementOperations2D elemOps;

    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::FakeParallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(numRuns);
}

TEST(MiddleMesh, CreateMiddleMeshFromMeshesWithLargeRefinementRatio2)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if (stk::middle_mesh::utils::impl::comm_size(comm) > 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 100);

  stk::middle_mesh::mesh::impl::MeshSpec spec{100, 100, 0, 1, 0, 1};
  stk::middle_mesh::mesh::impl::MeshSpec spec2{10, 10, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    batchTimer.start_batch_timer();

    auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
    auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);

    stk::middle_mesh::mesh::impl::ElementOperations2D elemOps;

    double eps = 1e-12;
    stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
    snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

    stk::middle_mesh::MiddleGridOpts middleMeshOpts;
    middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
    middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
    auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::FakeParallel,
                                                                    mesh1, mesh2, comm, nullptr,
                                                                    stk::middle_mesh::ParallelSearchOpts(),
                                                                    stk::middle_mesh::VolumeSnapOpts(),
                                                                    snapOpts, middleMeshOpts);
    interface->create_middle_grid();

    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(numRuns);
}

TEST(MiddleMesh, CreateMiddleMeshSPMD)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  auto commSize = stk::middle_mesh::utils::impl::comm_size(comm);
  if (stk::middle_mesh::utils::impl::comm_size(comm) <= 1)
    GTEST_SKIP();

  stk::unit_test_util::BatchTimer batchTimer(comm);
  int numRuns = stk::unit_test_util::simple_fields::get_command_line_option<unsigned>("-r", 50);
  int numMeshes = stk::unit_test_util::simple_fields::get_command_line_option<int>("-m", 50);

  stk::middle_mesh::mesh::impl::MeshSpec spec{commSize*2, commSize*2, 0, 1, 0, 1};
  auto func = [&](stk::middle_mesh::utils::Point const& pt) { return pt; };

  batchTimer.initialize_batch_timer();
  for (int i = 0; i < numRuns; ++i) {
    batchTimer.start_batch_timer();
    for (int j = 1; j < numMeshes; ++j) {
      stk::middle_mesh::mesh::impl::MeshSpec spec2{commSize*2+j, commSize*2+j, 0, 1, 0, 1};
      auto mesh1 = stk::middle_mesh::mesh::impl::create_mesh(spec, func);
      auto mesh2 = stk::middle_mesh::mesh::impl::create_mesh(spec2, func);

      stk::middle_mesh::mesh::impl::ElementOperations2D elemOps;

      double eps = 1e-12;
      stk::middle_mesh::BoundarySnapAndQualityImprovementOpts snapOpts;
      snapOpts.type = stk::middle_mesh::BoundarySnapAndQualityImprovementType::None;

      stk::middle_mesh::MiddleGridOpts middleMeshOpts;
      middleMeshOpts.normalProjectionOpts.classifierTolerances = stk::middle_mesh::impl::PointClassifierNormalWrapperTolerances(eps);
      middleMeshOpts.normalProjectionOpts.edgeTracerTolerances = stk::middle_mesh::impl::EdgeTracerTolerances(eps);
      auto interface = stk::middle_mesh::application_interface_factory(stk::middle_mesh::ApplicationInterfaceType::FakeParallel,
                                                                      mesh1, mesh2, comm, nullptr,
                                                                      stk::middle_mesh::ParallelSearchOpts(),
                                                                      stk::middle_mesh::VolumeSnapOpts(),
                                                                      snapOpts, middleMeshOpts);
      interface->create_middle_grid();
    }
    batchTimer.stop_batch_timer();
  }
  batchTimer.print_batch_timing(numRuns);
}

}


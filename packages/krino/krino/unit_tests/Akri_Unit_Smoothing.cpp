#include <Akri_AuxMetaData.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_MeshFromFileFixture.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_Quality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_Smooth.hpp>
#include <Akri_Unit_LogRedirecter.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <random>

namespace krino {

void move_node(const CoordinatesFieldRef coordsField, const stk::mesh::Entity node, const stk::math::Vector3d & newNodeLocation)
{
  double * coords = field_data<double>(coordsField, node);
    for (unsigned d = 0; d < coordsField.dim(); ++d)
      coords[d] = newNodeLocation[d];
}

void move_node(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::math::Vector3d & newNodeLocation)
{
  const CoordinatesFieldRef coordsField(mesh.mesh_meta_data().coordinate_field(), mesh.mesh_meta_data().spatial_dimension());
  return move_node(coordsField, node, newNodeLocation);
}

stk::math::Vector3d get_node_location(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  const CoordinatesFieldRef coordsField(mesh.mesh_meta_data().coordinate_field(), mesh.mesh_meta_data().spatial_dimension());
  return get_vector_field(mesh, coordsField, node, coordsField.dim());
}

template <typename MESHSPEC>
class SmoothingFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
public:
  SmoothingFixture()
  {
  }

  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mBuilder;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mComm;
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::write_mesh;

  double get_element_quality(const stk::mesh::Entity elem)
  {
    return element_quality(mMesh, mMesh.mesh_meta_data().spatial_dimension(), mMesh.mesh_meta_data().coordinate_field(), elem);
  }
protected:
  MESHSPEC meshSpec;
};

void optimize_node_location_using_krino_lbfgs(const NodeObjFn &objFn, const NodeObjFnSens &gradObjFn, stk::math::Vector3d& nodeLocation)
{
  const double xTol = 1.e-4;
  krino::lbfgs(objFn, gradObjFn, nodeLocation, xTol);
}

void optimize_node_location_using_krino_steepest_descent(const NodeObjFn &objFn, const NodeObjFnSens &gradObjFn, stk::math::Vector3d& nodeLocation)
{
  const double xTol = 1.e-4;
  krino::steepest_descent(objFn, gradObjFn, nodeLocation, xTol);
}

void optimize_node_location_using_ROL_lbfgs(const NodeObjFn &objFn, const NodeObjFnSens &gradObjFn, stk::math::Vector3d& nodeLocation)
{
  const double xTol = 1.e-4;
  rol_optimize(objFn, gradObjFn, nodeLocation, xTol);
}

void optimize_mesh_node_locations_using_krino_lbfgs(const MeshNodesObjFn &objFn, const MeshNodesObjFnSens &gradObjFn, DistributedVector& nodeLocations)
{
  const double xTol = 1.e-4;
  krino::lbfgs(objFn, gradObjFn, nodeLocations, xTol);
}

void optimize_mesh_node_locations_using_krino_steepest_descent(const MeshNodesObjFn &objFn, const MeshNodesObjFnSens &gradObjFn, DistributedVector& nodeLocations)
{
  const double xTol = 1.e-4;
  krino::steepest_descent(objFn, gradObjFn, nodeLocations, xTol);
}

void optimize_mesh_node_locations_using_ROL_lbfgs(const MeshNodesObjFn &objFn, const MeshNodesObjFnSens &gradObjFn, DistributedVector& nodeLocations)
{
  const double xTol = 1.e-4;
  rol_optimize(objFn, gradObjFn, nodeLocations, xTol);
}

class CubeOf12TetsSmoothing : public SmoothingFixture<CubeOf12Tets>
{
public:
  CubeOf12TetsSmoothing()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1,1,1,1,1,1,1,1,1});
  }
  stk::mesh::Entity get_center_node() { return get_assigned_node_for_index(8); }
  stk::math::Vector3d get_original_center_node_location() { return stk::math::Vector3d::ZERO; }
  stk::math::Vector3d get_center_node_location() { return get_node_location(mMesh, get_center_node()); }
  void displace_center_node(const stk::math::Vector3d & disp)
  {
    const stk::math::Vector3d newNodeLocation = get_original_center_node_location() + disp;
    move_node(mMesh, get_center_node(), newNodeLocation);
  }
  CoordinatesFieldRef get_coordinates_field() { const CoordinatesFieldRef coordsField(mMesh.mesh_meta_data().coordinate_field(), mMesh.mesh_meta_data().spatial_dimension()); return coordsField; }

  void test_smoothing_of_center_node(const OptimizeNodeLocation & optimize_node_location)
  {
    displace_center_node(stk::math::Vector3d(0.3,-0.3,0.3));

    stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    improve_quality_by_optimized_mean_ratio_smoothing_on_interior(mMesh, get_coordinates_field(), elemSelector, optimize_node_location);

    expect_near_absolute(get_center_node_location(), get_original_center_node_location(), get_location_tolerance());
  }

  double get_mesh_size() { return 1.0; }
  double get_location_tolerance() { return 1.e-3; }
};

TEST_F(CubeOf12TetsSmoothing, perturbedNode_smoothNodeUsingODT_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    displace_center_node(stk::math::Vector3d(0.3,-0.3,0.3));

    stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    improve_quality_by_ODT_smoothing_on_interior(mMesh, get_coordinates_field(), elemSelector);

    expect_near_absolute(get_center_node_location(), get_original_center_node_location(), get_location_tolerance());
  }
}

TEST_F(CubeOf12TetsSmoothing, perturbedNode_smoothNodeUsingQuality_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    test_smoothing_of_center_node(optimize_node_location_using_krino_lbfgs);
    test_smoothing_of_center_node(optimize_node_location_using_krino_steepest_descent);
    test_smoothing_of_center_node(optimize_node_location_using_ROL_lbfgs);
  }
}

static void zero_out_components_on_faces_of_2x2x2_cube(const stk::math::Vector3d& coords, stk::math::Vector3d& dir)
{
  for (int i=0; i<3; ++i)
    if (coords[i] == -1. || coords[i] == 1.)
      dir[i] = 0.;
}

std::function<void(stk::mesh::Entity, stk::math::Vector3d&)> build_node_search_direction_filter_for_2x2x2_cube(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField)
{
  const std::function<void(stk::mesh::Entity, stk::math::Vector3d&)> fn = [coordsField](stk::mesh::Entity node, stk::math::Vector3d& dir)
  {
    const stk::math::Vector3d nodeCoords(field_data<double>(coordsField, node), coordsField.dim());
    zero_out_components_on_faces_of_2x2x2_cube(nodeCoords, dir);
  };
  return fn;
}

class CubeOf24TetsSmoothing : public SmoothingFixture<CubeOf24Tets>
{
public:
  CubeOf24TetsSmoothing()
  {
    set_valid_proc_sizes_for_test({1,2,4});
    const std::vector<unsigned> elementBlockIDs = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, {0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1});
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,1,2,3,0,1,2,3});

    originalNodeLocs = meshSpec.nodeLocs;
  }

  void displace_center_and_face_nodes()
  {
    std::uniform_real_distribution<> randomDist(-0.99, 0.99);
    for (int i=8; i<15; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const stk::math::Vector3d nodeCoords = get_node_location(mMesh, node);
        stk::math::Vector3d disp(randomDist(randomizer), randomDist(randomizer), randomDist(randomizer));
        zero_out_components_on_faces_of_2x2x2_cube(nodeCoords, disp);
        move_node(mMesh, node, nodeCoords+disp);
      }
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  void expect_restored_center_and_face_nodes()
  {
    for (int i=8; i<15; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        expect_near_absolute(originalNodeLocs[i], get_node_location(mMesh, node), get_location_tolerance());
      }
    }
  }

  CoordinatesFieldRef get_coordinates_field() { const CoordinatesFieldRef coordsField(mMesh.mesh_meta_data().coordinate_field(), mMesh.mesh_meta_data().spatial_dimension()); return coordsField; }

  double get_location_tolerance() { return 1.e-2; }

  std::function<void(stk::mesh::Entity, stk::math::Vector3d&)> build_node_search_direction_filter()
  {
    return build_node_search_direction_filter_for_2x2x2_cube(get_mesh(), get_coordinates_field());
  }

  void test_node_at_a_time_smoothing(const OptimizeNodeLocation & optimize_node_location)
  {
    displace_center_and_face_nodes();
    improve_quality_by_optimized_mean_ratio_smoothing(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), optimize_node_location, build_node_search_direction_filter());
    expect_restored_center_and_face_nodes();
  }

  void test_simultaneous_node_mesh_smoothing(const OptimizeMeshNodeLocations & optimize_mesh_node_locations)
  {
    displace_center_and_face_nodes();
    improve_quality_by_simultaneous_optimized_mean_ratio_smoothing(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), optimize_mesh_node_locations, build_node_search_direction_filter());
    expect_restored_center_and_face_nodes();
  }

  std::vector<stk::math::Vector3d> originalNodeLocs;
  std::mt19937 randomizer;
};


TEST_F(CubeOf24TetsSmoothing, perturbedNodes_simultaneousSmoothNodeUsingQuality_restoreOriginalTets)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    std::cout << "ROL L-BFGS\n";
    test_simultaneous_node_mesh_smoothing(optimize_mesh_node_locations_using_ROL_lbfgs);
  }
  else if(is_valid_proc_size_for_test())
  {
    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "Krino L-BFGS\n";
    test_simultaneous_node_mesh_smoothing(optimize_mesh_node_locations_using_krino_lbfgs);

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "Krino Steepest descent\n";
    test_simultaneous_node_mesh_smoothing(optimize_mesh_node_locations_using_krino_steepest_descent);
  }
}

TEST_F(CubeOf24TetsSmoothing, perturbedNodes_nodeAtATimeSmoothNodeUsingQuality_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    LogRedirecter log;

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "Krino L-BFGS\n";
    test_node_at_a_time_smoothing(optimize_node_location_using_krino_lbfgs);

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "Krino Steepest descent\n";
    test_node_at_a_time_smoothing(optimize_node_location_using_krino_steepest_descent);

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "ROL L-BFGS\n";
    test_node_at_a_time_smoothing(optimize_node_location_using_ROL_lbfgs);

    std::cout << log.get_log() << std::endl;
  }
}

using Smoother = MeshFromFileFixture;

TEST_F(Smoother, smoothUsingODT)
{
  if (read_mesh_if_present_and_supported("mesh.g"))
  {
    LogRedirecter log;

    improve_quality_by_ODT_smoothing_on_interior(get_mesh(), get_coordinates_field(), get_aux_meta().active_part());

    output_composed_mesh_with_fields(get_mesh(), get_mesh().mesh_meta_data().universal_part(), "meshSmoothedODT.e", 1, 0.);

    if (0 == stk::parallel_machine_rank(myComm))
      std::cout << log.get_log() << std::endl;
  }
}

TEST_F(Smoother, smoothNodeAtATimeUsingOptimizedMeanRatioOn2x2x2Cube)
{
  if (read_mesh_if_present_and_supported("mesh.g"))
  {
    LogRedirecter log;

    const auto node_search_direction_filter = build_node_search_direction_filter_for_2x2x2_cube(get_mesh(), get_coordinates_field());
    improve_quality_by_optimized_mean_ratio_smoothing(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), optimize_node_location_using_krino_lbfgs, node_search_direction_filter);

    output_composed_mesh_with_fields(get_mesh(), get_mesh().mesh_meta_data().universal_part(), "meshSmoothNodeAtATimeOptimizedMeanRatio.e", 1, 0.);

    if (0 == stk::parallel_machine_rank(myComm))
      std::cout << log.get_log() << std::endl;
  }
}

TEST_F(Smoother, smoothNodesSimultaneouslyUsingOptimizedMeanRatioOn2x2x2Cube)
{
  if (read_mesh_if_present_and_supported("mesh.g"))
  {
    LogRedirecter log;

    const auto node_search_direction_filter = build_node_search_direction_filter_for_2x2x2_cube(get_mesh(), get_coordinates_field());
    improve_quality_by_simultaneous_optimized_mean_ratio_smoothing(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), optimize_mesh_node_locations_using_krino_steepest_descent, node_search_direction_filter);

    output_composed_mesh_with_fields(get_mesh(), get_mesh().mesh_meta_data().universal_part(), "meshSmoothNodesSimultaneouslyOptimizedMeanRatio.e", 1, 0.);

    if (0 == stk::parallel_machine_rank(myComm))
      std::cout << log.get_log() << std::endl;
  }
}



}



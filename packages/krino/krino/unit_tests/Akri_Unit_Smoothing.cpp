#include <Akri_AuxMetaData.hpp>
#include <Akri_StkMeshFixture.hpp>
#include <Akri_MeshFromFileFixture.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_Objective_MeanRatioAboutNode.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothMesh.hpp>
#include <Akri_SmoothNode.hpp>
#include <Akri_LogRedirecter.hpp>
#include <Akri_Objective_KinematicMesh.hpp>
#include <Akri_Objective_ModifiedNeoHookean.hpp>
#include <Akri_Objective_SethHill.hpp>
#include <Akri_Objective_SizeShape.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <random>

namespace krino {

TEST(MeanRatioNormObjective, behaviorForPositiveZeroAndInvertedElements)
{
  const double positiveMeanQuality =  1.e-5;
  const double zeroMeanQuality = 0;
  const double invertedMeanQuality = -1.e-5;

  // default behavior
  EXPECT_TRUE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(positiveMeanQuality)));
  EXPECT_FALSE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(zeroMeanQuality)));
  EXPECT_FALSE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(invertedMeanQuality)));

  // unlimited behavior
  MeanRatioNormObjective::set_quality_epsilon_for_handling_inversion(0.);
  EXPECT_TRUE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(positiveMeanQuality)));
  EXPECT_FALSE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(zeroMeanQuality)));
  EXPECT_FALSE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(invertedMeanQuality)));

  // limited behavior
  MeanRatioNormObjective::set_quality_epsilon_for_handling_inversion(1.e-3);
  EXPECT_TRUE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(positiveMeanQuality)));
  EXPECT_TRUE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(zeroMeanQuality)));
  EXPECT_TRUE(std::isfinite(MeanRatioNormObjective::global_objective_element_contribution(invertedMeanQuality)));
}

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

static std::vector<double> get_owned_element_sizes(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::EntityId> & elemIds, const std::vector<double> elemSizes, const std::vector<stk::mesh::Entity> & ownedElems)
{
  std::vector<double> ownedElemSizes;
  ownedElemSizes.reserve(ownedElems.size());
  for (auto elem : ownedElems)
  {
    const auto iter = std::find(elemIds.begin(), elemIds.end(), mesh.identifier(elem)); // slow, but fine for test
    STK_ThrowRequire(iter != elemIds.end());
    ownedElemSizes.push_back(elemSizes[std::distance(elemIds.begin(), iter)]);
  }
  return ownedElemSizes;
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

  CoordinatesFieldRef get_coordinates_field() { const CoordinatesFieldRef coordsField(mMesh.mesh_meta_data().coordinate_field(), mMesh.mesh_meta_data().spatial_dimension()); return coordsField; }

  double compute_mean_ratio_quality()
  {
    MeanRatioQualityMetric qualityMetric;
    return compute_mesh_quality(mMesh, mMesh.mesh_meta_data().universal_part(), mMesh.mesh_meta_data().coordinate_field(), qualityMetric);
  }
protected:
  MESHSPEC meshSpec;
  MeanRatioNormObjective meanRatioObj;
  SizeShapeObjective sizeShapeObj;
  ModifiedNeoHookeanObjective modNeoHookean;
  SethHillParams sethHillParams{1.0, 1.0, 2, 2};
  SethHillSmoothingObjective sethHillObj{sethHillParams};
};

void optimize_node_location_using_krino_lbfgs(const Objective3DInterface &objFn, stk::math::Vector3d& nodeLocation)
{
  const double xTol = 1.e-4;
  krino::lbfgs(objFn, nodeLocation, xTol);
}

void optimize_node_location_using_krino_steepest_descent(const Objective3DInterface &objFn, stk::math::Vector3d& nodeLocation)
{
  const double xTol = 1.e-4;
  krino::steepest_descent(objFn, nodeLocation, xTol);
}

void optimize_node_location_using_ROL_lbfgs(const Objective3DInterface &objFn, stk::math::Vector3d& nodeLocation)
{
  const double xTol = 1.e-5;
  rol_optimize(objFn, nodeLocation, xTol);
}

void optimize_mesh_using_krino_lbfgs(const ObjectiveInterface &objFn, DistributedVector& nodeLocations)
{
  const double xTol = 1.e-4;
  const double initialEnergy = objFn.compute_value(nodeLocations);
  krino::lbfgs(objFn, nodeLocations, xTol, 0, initialEnergy);
}

void optimize_mesh_using_krino_steepest_descent(const ObjectiveInterface &objFn, DistributedVector& nodeLocations)
{
  const double xTol = 1.e-4;
  krino::steepest_descent(objFn, nodeLocations, xTol);
}

void optimize_mesh_using_ROL_lbfgs(const ObjectiveInterface &objFn, DistributedVector& nodeLocations)
{
  const double xTol = 1.e-4;
  rol_optimize(objFn, nodeLocations, xTol);
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

  void test_smoothing_of_center_node(const Optimizer3DFunction & optimizer)
  {
    displace_center_node(stk::math::Vector3d(0.3,-0.3,0.3));

    const stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    const stk::mesh::Selector interiorNodeSelector = elemSelector & !get_aux_meta().exposed_boundary_part() & !get_aux_meta().block_boundary_part();
    MeanRatioAboutNodeObjective nodeCoordObj;
    smooth_mesh_by_optimizing_interior_node_locations(mMesh, get_coordinates_field(), elemSelector, interiorNodeSelector, optimizer, nodeCoordObj);

    expect_near_absolute(get_center_node_location(), get_original_center_node_location(), get_location_tolerance());
  }

  double get_location_tolerance() { return 1.e-3; }
};

TEST_F(CubeOf12TetsSmoothing, perturbedNode_smoothNodeUsingODT_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    displace_center_node(stk::math::Vector3d(0.3,-0.3,0.3));

    const stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    const stk::mesh::Selector interiorNodeSelector = elemSelector & !get_aux_meta().exposed_boundary_part() & !get_aux_meta().block_boundary_part();
    smooth_mesh_by_ODT_of_interior_node_locations(mMesh, get_coordinates_field(), elemSelector, interiorNodeSelector);

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
  for (unsigned i=0; i<3; ++i)
    if (coords[i] == -1. || coords[i] == 1.)
      dir[i] = 0.;
}

GradientFilter build_node_search_direction_filter_for_2x2x2_cube()
{
  const GradientFilter fn = [](const stk::mesh::Entity node, const stk::math::Vector3d& nodeCoords, stk::math::Vector3d& dir)
  {
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

    mGradientFilter = build_node_search_direction_filter();
  }

  void test_center_node_location(const stk::math::Vector3d & goldLoc)
  {
    const unsigned centerNodeIndex = 14;
    const stk::mesh::Entity node = get_assigned_node_for_index(centerNodeIndex);
    if (mMesh.is_valid(node))
    {
      expect_near_absolute(goldLoc, get_node_location(mMesh, node), get_location_tolerance());
    }
  }

  std::array<stk::math::Vector3d,8> face_and_center_node_random_displacements()
  {
    std::uniform_real_distribution<> randomDist(-0.99, 0.99);
    std::array<stk::math::Vector3d,8> disp;
    for (auto & nodeDisp : disp)
      nodeDisp = stk::math::Vector3d(randomDist(randomizer), randomDist(randomizer), randomDist(randomizer));
    return disp;
  }

  std::array<stk::math::Vector3d,8> face_and_center_node_displacements_with_magnitude(const double mag)
  {
    const std::array<stk::math::Vector3d,8> disp =
        {{
            { 0.0, -0.3, -0.1},
            { 0.0,  mag,  0.4},
            { 0.1,  0.0,  0.8},
            { 0.1,  0.0,  0.0},
            { 0.2,  0.7,  0.0},
            {-mag,  mag,  0.0},
            { 0.4, -0.3,  mag}
        }};
    return disp;
  }

  std::array<stk::math::Vector3d,8> face_and_center_node_displacements_for_stiff_case_found_randomly()
  {
    return face_and_center_node_displacements_with_magnitude(0.99);
  }

  std::array<stk::math::Vector3d,8> face_and_center_node_displacements_for_less_stiff_case()
  {
    return face_and_center_node_displacements_with_magnitude(0.8);
  }

  void displace_face_and_center_nodes(const std::array<stk::math::Vector3d,8> & faceAndCenterNodeDisp)
  {
    for (int i=8; i<15; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const stk::math::Vector3d origNodeCoords = meshSpec.nodeLocs[i];
        stk::math::Vector3d disp = faceAndCenterNodeDisp[i-8];
        zero_out_components_on_faces_of_2x2x2_cube(origNodeCoords, disp);
        move_node(mMesh, node, origNodeCoords+disp);
      }
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  bool are_face_and_center_nodes_restored()
  {
    bool areNodesRestored = true;
    for (int i=8; i<15; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        if ((meshSpec.nodeLocs[i]-get_node_location(mMesh, node)).length() > get_location_tolerance())
        {
          areNodesRestored = false;
          std::cout << "Node " << i << " at " << get_node_location(mMesh, node) << " too far from gold " << meshSpec.nodeLocs[i] << std::endl;
        }
      }
    }
    return areNodesRestored;
  }

  double get_location_tolerance() { return 1.e-2; }

  GradientFilter build_node_search_direction_filter()
  {
    return build_node_search_direction_filter_for_2x2x2_cube();
  }

  void test_node_at_a_time_smoothing(const Optimizer3DFunction & optimizer, const std::array<stk::math::Vector3d,8> & faceAndCenterNodeDisp)
  {
    displace_face_and_center_nodes(faceAndCenterNodeDisp);
    MeanRatioAboutNodeObjective nodeCoordObj;
    smooth_mesh_by_optimizing_node_locations(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), optimizer, nodeCoordObj, build_node_search_direction_filter());
    EXPECT_TRUE(are_face_and_center_nodes_restored());
  }

  void test_mesh_smoothing(const OptimizerFunction & optimizer, MeshObjectiveBase & smoothingObj, const std::string & probDesc)
  {
    if (0 == stk::parallel_machine_rank(mComm)) std::cout << probDesc << std::endl;
    smoothingObj.set_gradient_filter(mGradientFilter);
    smooth_mesh_by_optimization(get_mesh(), get_coordinates_field(), optimizer, smoothingObj);
    EXPECT_TRUE(are_face_and_center_nodes_restored()) << probDesc;
  }

  void reset_coordinates()
  {
    for (int i=0; i<15; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
        move_node(mMesh, node, meshSpec.nodeLocs[i]);
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  void reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(const KinematicElementObjective & elemObj, const std::vector<double> & elemSizes, const stk::math::Vector3d & goldCenterNodeLoc)
  {
    reset_coordinates();
    KinematicMeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), elemObj);

    smoothingObj.set_element_ref_sizes(get_owned_element_sizes(mMesh, get_assigned_element_global_ids(), elemSizes, smoothingObj.get_owned_elements()));

    smoothingObj.set_gradient_filter(mGradientFilter);
    smooth_mesh_by_optimization(get_mesh(), get_coordinates_field(), optimize_mesh_using_krino_lbfgs, smoothingObj);

    test_center_node_location(goldCenterNodeLoc);
  }

  void reset_coords_and_test_mean_ratio_mesh_smoothing(const OptimizerFunction & optimizer, const std::array<stk::math::Vector3d,8> & faceAndCenterNodeDisp, const std::string & probDesc)
  {
    displace_face_and_center_nodes(faceAndCenterNodeDisp);

    MeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), meanRatioObj);
    test_mesh_smoothing(optimizer, smoothingObj, probDesc);
  }

  void reset_coords_and_test_kinematic_mesh_smoothing(const OptimizerFunction & optimizer, const KinematicElementObjective & elemObj, const std::array<stk::math::Vector3d,8> & faceAndCenterNodeDisp, const std::string & probDesc)
  {
    displace_face_and_center_nodes(faceAndCenterNodeDisp);

    KinematicMeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), elemObj);
    smoothingObj.set_constant_element_ref_size(1.0); // constant size, magnitude is not relevant
    test_mesh_smoothing(optimizer, smoothingObj, probDesc);
  }

  std::mt19937 randomizer;
  GradientFilter mGradientFilter;
};


TEST_F(CubeOf24TetsSmoothing, perturbedNodes_simultaneousSmoothNodeUsingQuality_restoreOriginalTets)
{
  LogRedirecter log;

  const auto stiff = face_and_center_node_displacements_for_stiff_case_found_randomly();
  const auto lessStiff = face_and_center_node_displacements_for_less_stiff_case();

  if(stk::parallel_machine_size(mComm) == 1)
  {
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, lessStiff, "ROL L-BFGS using mean ratio");
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, stiff, "ROL L-BFGS using mean ratio (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, lessStiff, "ROL L-BFGS using size-shape");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, stiff, "ROL L-BFGS using size-shape (stiff)");

    // only lessStiff because ROL BFGS can't handle Seth-Hill on it
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sethHillObj, lessStiff, "ROL L-BFGS using Seth-Hill");
  }
  if (is_valid_proc_size_for_test())
  {
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_lbfgs, lessStiff, "Krino L-BFGS using mean ratio");
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_lbfgs, stiff, "Krino L-BFGS using mean ratio (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, lessStiff, "Krino L-BFGS using size-shape");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, stiff, "Krino L-BFGS using size-shape (stiff)");

    // only lessStiff because krino BFGS can't handle Seth-Hill on it
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sethHillObj, lessStiff, "Krino L-BFGS using Seth-Hill");

    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, lessStiff, "Krino Steepest descent using mean ratio");
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, stiff, "Krino Steepest descent using mean ratio (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, lessStiff, "Krino Steepest descent using size-shape");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, stiff, "Krino Steepest descent using size-shape (stiff)");

    // only lessStiff because krino steepest_descent can't handle Seth-Hill on it
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sethHillObj, lessStiff, "Krino Steepest descent using Seth-Hill");
  }

  if (0 == stk::parallel_machine_rank(mComm))
      std::cout << log.get_log() << std::endl;
}

TEST_F(CubeOf24TetsSmoothing, perturbedNodes_nodeAtATimeSmoothNodeUsingQuality_restoreOriginalTets)
{
  if(is_valid_proc_size_for_test())
  {
    LogRedirecter log;
    const auto stiff = face_and_center_node_displacements_for_stiff_case_found_randomly();

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "Krino L-BFGS\n";
    test_node_at_a_time_smoothing(optimize_node_location_using_krino_lbfgs, stiff);

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "Krino Steepest descent\n";
    test_node_at_a_time_smoothing(optimize_node_location_using_krino_steepest_descent, stiff);

    if (0 == stk::parallel_machine_rank(mComm)) std::cout << "ROL L-BFGS\n";
    test_node_at_a_time_smoothing(optimize_node_location_using_ROL_lbfgs, stiff);

    std::cout << log.get_log() << std::endl;
  }
}

using Smoother = MeshFromFileFixture;

TEST_F(Smoother, smoothUsingODT)
{
  if (read_mesh_if_present_and_supported("mesh.g"))
  {
    LogRedirecter log;

    const stk::mesh::Selector elemSelector(get_aux_meta().active_part());
    const stk::mesh::Selector interiorNodeSelector = elemSelector & !get_aux_meta().exposed_boundary_part() & !get_aux_meta().block_boundary_part();
    smooth_mesh_by_ODT_of_interior_node_locations(get_mesh(), get_coordinates_field(), elemSelector, interiorNodeSelector);

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

    const auto node_search_direction_filter = build_node_search_direction_filter_for_2x2x2_cube();
    MeanRatioAboutNodeObjective nodeCoordObj;
    smooth_mesh_by_optimizing_node_locations(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), optimize_node_location_using_krino_lbfgs, nodeCoordObj, node_search_direction_filter);

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

    const auto node_search_direction_filter = build_node_search_direction_filter_for_2x2x2_cube();
    const SolutionProjector empty_solution_projector;
    MeanRatioNormObjective smoothingObj;
    smooth_mesh_by_optimization(get_mesh(),
        get_coordinates_field(),
        get_aux_meta().active_part(),
        optimize_mesh_using_krino_lbfgs,
        smoothingObj,
        node_search_direction_filter,
        empty_solution_projector);

    output_composed_mesh_with_fields(get_mesh(), get_mesh().mesh_meta_data().universal_part(), "meshSmoothNodesSimultaneouslyOptimizedMeanRatio.e", 1, 0.);

    if (0 == stk::parallel_machine_rank(myComm))
      std::cout << log.get_log() << std::endl;
  }
}

static void zero_out_normal_components_for_sphere(const stk::math::Vector3d & center, const double radius, const stk::math::Vector3d& coords, stk::math::Vector3d& dir)
{
  const stk::math::Vector3d normal = (coords-center).unit_vector();
  dir = dir - Dot(normal, dir)*normal;
}

static void project_onto_sphere(const stk::math::Vector3d & center, const double radius, stk::math::Vector3d& coords)
{
  const stk::math::Vector3d normal = (coords-center).unit_vector();
  coords = center + radius*normal;
}

GradientFilter build_node_search_direction_filter_for_sphere(const stk::math::Vector3d & center, const double radius)
{
  const GradientFilter fn = [center, radius](const stk::mesh::Entity node, const stk::math::Vector3d& nodeCoords, stk::math::Vector3d& dir)
  {
    zero_out_normal_components_for_sphere(center, radius, nodeCoords, dir);
  };
  return fn;
}

SolutionProjector build_solution_projector_for_sphere(const stk::math::Vector3d & center, const double radius)
{
  const SolutionProjector fn = [center, radius](const stk::mesh::Entity node, stk::math::Vector3d& nodeCoords)
  {
    project_onto_sphere(center, radius, nodeCoords);
  };
  return fn;
}

class RightTetOnSphereSmoothing : public SmoothingFixture<RightTet>
{
public:
  RightTetOnSphereSmoothing()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
    mGradientFilter = build_node_search_direction_filter();
    mSolutionProjector = build_solution_projector();
  }

  void expect_nodes_on_sphere(const double tol)
  {
    const stk::math::Vector3d sphereCenter = sphere_center();
    const double sphereRadius = sphere_radius();
    for (int i=0; i<4; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const double distErr = (get_node_location(mMesh, node) - sphereCenter).length() - sphereRadius;
        EXPECT_LE(std::abs(distErr), tol);
      }
    }
  }

  void reset_coordinates()
  {
    for (int i=0; i<4; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
        move_node(mMesh, node, meshSpec.nodeLocs[i]);
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  void test_mesh_smoothing(const OptimizerFunction & optimizer, MeshObjectiveBase & smoothingObj, const std::string & probDesc, const double locationTol, const double qualityTol)
  {
    if (0 == stk::parallel_machine_rank(mComm)) std::cout << probDesc << std::endl;
    smoothingObj.set_gradient_filter(mGradientFilter);
    smoothingObj.set_solution_projector(mSolutionProjector);
    smooth_mesh_by_optimization(get_mesh(), get_coordinates_field(), optimizer, smoothingObj);

    EXPECT_GT(compute_mean_ratio_quality(), qualityTol) << probDesc;
    expect_nodes_on_sphere(locationTol);
  }

  void reset_coords_and_test_mean_ratio_mesh_smoothing(const OptimizerFunction & optimizer, const std::string & probDesc, const double locationTol, const double qualityTol)
  {
    reset_coordinates();

    MeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), meanRatioObj);
    test_mesh_smoothing(optimizer, smoothingObj, probDesc, locationTol, qualityTol);
  }

  void reset_coords_and_test_kinematic_mesh_smoothing(const OptimizerFunction & optimizer,
    const KinematicElementObjective & elemObj,
    const std::string & probDesc,
    const double sizeChangeRatio,
    const double locationTol,
    const double qualityTol)
  {
    reset_coordinates();

    KinematicMeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), elemObj);
    const double edgeLenForInscribingSphere = 2.*std::sqrt(6.)/3. * sphere_radius();
    smoothingObj.set_constant_element_ref_size(sizeChangeRatio*edgeLenForInscribingSphere);
    test_mesh_smoothing(optimizer, smoothingObj, probDesc, locationTol, qualityTol);
  }

  GradientFilter build_node_search_direction_filter()
  {
    return build_node_search_direction_filter_for_sphere(sphere_center(), sphere_radius());
  }

  SolutionProjector build_solution_projector()
  {
    return build_solution_projector_for_sphere(sphere_center(), sphere_radius());
  }

  stk::math::Vector3d sphere_center() const { return stk::math::Vector3d(0.5,0.5,0.5); }
  double sphere_radius() const { return std::sqrt(3.0)/2.; }
  GradientFilter mGradientFilter;
  SolutionProjector mSolutionProjector;
};

TEST_F(RightTetOnSphereSmoothing, smoothNodesSimultaneouslyUsingOptimizedMeanRatioOnSphere)
{
  const double distFromSphereTol = 1.e-9;

  const double baseQualityTol = 0.999;
  const double qualityWithSizeChangeTol= 0.88;
  const double sizeChangeRatio = 0.8;

  LogRedirecter log;

  if(stk::parallel_machine_size(mComm) == 1)
  {
    // ROL does not implement solution projection
    const double distFromSphereNoProjTol = 0.12;
    const double distFromSphereNoProjCorrectVolTol = 0.01;
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, "ROL L-BFGS using mean ratio", distFromSphereNoProjTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, "ROL L-BFGS using size-shape (matching elem vol)", 1.0, distFromSphereNoProjCorrectVolTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, "ROL L-BFGS using size-shape (non-matching elem vol)", sizeChangeRatio, distFromSphereNoProjTol, qualityWithSizeChangeTol);
    const double qualitySethHillTol = 0.98; // Due to stiffness?
    const double distFromSphereNoProjCorrectVolSethHillTol = 0.21; // Due to stiffness? This is quite bad.
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sethHillObj, "ROL L-BFGS using Seth-Hill (matching elem vol)", 1.0, distFromSphereNoProjCorrectVolSethHillTol, qualitySethHillTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sethHillObj, "ROL L-BFGS using Seth-Hill (non-matching elem vol)", sizeChangeRatio, distFromSphereNoProjTol, qualityWithSizeChangeTol);
  }
  if(is_valid_proc_size_for_test())
  {
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_lbfgs, "Krino L-BFGS using mean ratio", distFromSphereTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, "Krino L-BFGS using size-shape (matching elem vol)", 1.0, distFromSphereTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, "Krino L-BFGS using size-shape (non-matching elem vol)", sizeChangeRatio, distFromSphereTol, qualityWithSizeChangeTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sethHillObj, "Krino L-BFGS using Seth-Hill (matching elem vol)", 1.0, distFromSphereTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sethHillObj, "Krino L-BFGS using Seth-Hill (non-matching elem vol)", sizeChangeRatio, distFromSphereTol, qualityWithSizeChangeTol);

    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, "Krino Steepest descent using mean ratio", distFromSphereTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, "Krino Steepest descent using size-shape (matching elem vol)", 1.0, distFromSphereTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, "Krino Steepest descent using size-shape (non-matching elem vol)", sizeChangeRatio, distFromSphereTol, qualityWithSizeChangeTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sethHillObj, "Krino Steepest descent using Seth-Hill (matching elem vol)", 1.0, distFromSphereTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sethHillObj, "Krino Steepest descent using Seth-Hill (non-matching elem vol)", sizeChangeRatio, distFromSphereTol, qualityWithSizeChangeTol);
  }

  if (0 == stk::parallel_machine_rank(mComm))
    std::cout << log.get_log() << std::endl;
}

class SquareOf8TrisSmoothing : public SmoothingFixture<SquareOf8Tris>
{
public:
  SquareOf8TrisSmoothing()
  {
    set_valid_proc_sizes_for_test({1,2,4});
    const std::vector<unsigned> elementBlockIDs = {1,1,1,1,1,1,1,1};
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, {0,1,0,1,0,1,0,1});
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, elementBlockIDs, {0,1,2,3,0,1,2,3});

    mGradientFilter = build_node_search_direction_filter();
  }

  stk::mesh::Entity get_center_node() { return get_assigned_node_for_index(8); }

  void test_center_node_location(const stk::math::Vector3d & goldLoc)
  {
    const unsigned centerNodeIndex = 8;
    const stk::mesh::Entity node = get_assigned_node_for_index(centerNodeIndex);
    if (mMesh.is_valid(node))
    {
      expect_near_absolute(goldLoc, get_node_location(mMesh, node), get_location_tolerance());
    }
  }

  void reset_coordinates()
  {
    for (int i=0; i<9; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const stk::math::Vector3d origNodeCoords(meshSpec.nodeLocs[i][0], meshSpec.nodeLocs[i][1], 0.);
        move_node(mMesh, node, origNodeCoords);
      }
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  void reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(const KinematicElementObjective & elemObj, const std::vector<double> & elemSizes, const stk::math::Vector3d & goldCenterNodeLoc)
  {
    reset_coordinates();
    KinematicMeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), elemObj);

    smoothingObj.set_element_ref_sizes(get_owned_element_sizes(mMesh, get_assigned_element_global_ids(), elemSizes, smoothingObj.get_owned_elements()));

    smoothingObj.set_gradient_filter(mGradientFilter);
    smooth_mesh_by_optimization(get_mesh(), get_coordinates_field(), optimize_mesh_using_krino_lbfgs, smoothingObj);

    test_center_node_location(goldCenterNodeLoc);
  }

  std::array<stk::math::Vector2d,5> face_and_side_node_displacements_with_magnitude(const double mag)
  {
    const std::array<stk::math::Vector2d,5> disp =
        {{
            {-mag,  0.0},
            { 0.0,  mag},
            { mag,  0.0},
            { 0.0,  0.0},
            {-mag,  mag},
        }};
    return disp;
  }

  void displace_side_and_center_nodes(const std::array<stk::math::Vector2d,5> & faceAndCenterNodeDisp)
  {
    for (int i=4; i<9; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const stk::math::Vector3d origNodeCoords(meshSpec.nodeLocs[i][0], meshSpec.nodeLocs[i][1], 0.);
        stk::math::Vector3d disp(faceAndCenterNodeDisp[i-4][0], faceAndCenterNodeDisp[i-4][1], 0.);
        zero_out_components_on_faces_of_2x2x2_cube(origNodeCoords, disp);
        move_node(mMesh, node, origNodeCoords+disp);
      }
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  bool are_side_and_center_nodes_restored()
  {
    bool areNodesRestored = true;
    for (int i=4; i<9; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const stk::math::Vector3d origNodeCoords(meshSpec.nodeLocs[i][0], meshSpec.nodeLocs[i][1], 0.);
        if ((origNodeCoords-get_node_location(mMesh, node)).length() > get_location_tolerance())
        {
          areNodesRestored = false;
          std::cout << "Node " << i << " at " << get_node_location(mMesh, node) << " too far from gold " << origNodeCoords << std::endl;
        }
      }
    }
    return areNodesRestored;
  }

  double get_location_tolerance() { return 1.e-2; }

  GradientFilter build_node_search_direction_filter()
  {
    return build_node_search_direction_filter_for_2x2x2_cube(); // works for square on z=0
  }

  void test_mesh_smoothing(const OptimizerFunction & optimizer, MeshObjectiveBase & smoothingObj, const std::string & probDesc)
  {
    if (0 == stk::parallel_machine_rank(mComm)) std::cout << probDesc << std::endl;
    smoothingObj.set_gradient_filter(mGradientFilter);
    smooth_mesh_by_optimization(get_mesh(), get_coordinates_field(), optimizer, smoothingObj);
    EXPECT_TRUE(are_side_and_center_nodes_restored()) << probDesc;
  }

  void reset_coords_and_test_mean_ratio_mesh_smoothing(const OptimizerFunction & optimizer, const std::array<stk::math::Vector2d,5> & sideAndCenterNodeDisp, const std::string & probDesc)
  {
    displace_side_and_center_nodes(sideAndCenterNodeDisp);

    MeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), meanRatioObj);
    test_mesh_smoothing(optimizer, smoothingObj, probDesc);
  }

  void reset_coords_and_test_kinematic_mesh_smoothing(const OptimizerFunction & optimizer, const KinematicElementObjective & elemObj, const std::array<stk::math::Vector2d,5> & sideAndCenterNodeDisp, const std::string & probDesc)
  {
    displace_side_and_center_nodes(sideAndCenterNodeDisp);

    KinematicMeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), elemObj);
    smoothingObj.set_constant_element_ref_size(1.0); // constant size, magnitude is not relevant
    test_mesh_smoothing(optimizer, smoothingObj, probDesc);
  }

  GradientFilter mGradientFilter;
};


TEST_F(SquareOf8TrisSmoothing, perturbedNodes_meshSmoothing_restoreOriginalTris)
{
  LogRedirecter log;

  const auto stiff = face_and_side_node_displacements_with_magnitude(0.98);
  const auto lessStiff = face_and_side_node_displacements_with_magnitude(0.8);

  if(stk::parallel_machine_size(mComm) == 1)
  {
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, lessStiff, "ROL L-BFGS using mean ratio");
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, stiff, "ROL L-BFGS using mean ratio (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, lessStiff, "ROL L-BFGS using size-shape");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, stiff, "ROL L-BFGS using size-shape (stiff)");
  }
  if (is_valid_proc_size_for_test())
  {
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_lbfgs, lessStiff, "Krino L-BFGS using mean ratio");
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_lbfgs, stiff, "Krino L-BFGS using mean ratio (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, lessStiff, "Krino L-BFGS using size-shape");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, stiff, "Krino L-BFGS using size-shape (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sethHillObj, lessStiff, "Krino L-BFGS using Seth-Hill");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sethHillObj, stiff, "Krino L-BFGS using Seth-Hill (stiff)");

    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, lessStiff, "Krino Steepest descent using mean ratio");
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, stiff, "Krino Steepest descent using mean ratio (stiff)");

    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, lessStiff, "Krino Steepest descent using size-shape");
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, stiff, "Krino Steepest descent using size-shape (stiff)");
  }

  if (0 == stk::parallel_machine_rank(mComm))
      std::cout << log.get_log() << std::endl;
}

class RightTriOnCircleSmoothing : public SmoothingFixture<RightTri>
{
public:
  RightTriOnCircleSmoothing()
  {
    set_valid_proc_sizes_for_test({1});
    StkMeshTriFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
    mGradientFilter = build_node_search_direction_filter();
    mSolutionProjector = build_solution_projector();
  }

  void expect_nodes_on_circle(const double tol)
  {
    const stk::math::Vector3d sphereCenter = sphere_center();
    const double circleRadius = circle_radius();
    for (int i=0; i<4; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const double distErr = (get_node_location(mMesh, node) - sphereCenter).length() - circleRadius;
        EXPECT_LE(std::abs(distErr), tol);
      }
    }
  }

  void reset_coordinates()
  {
    for (int i=0; i<4; ++i)
    {
      const stk::mesh::Entity node = get_assigned_node_for_index(i);
      if (mMesh.is_valid(node))
      {
        const stk::math::Vector3d origNodeCoords(meshSpec.nodeLocs[i][0], meshSpec.nodeLocs[i][1], 0.);
        move_node(mMesh, node, origNodeCoords);
      }
    }
    parallel_sync_fields(mMesh, {mMesh.mesh_meta_data().coordinate_field()});
  }

  void test_mesh_smoothing(const OptimizerFunction & optimizer, MeshObjectiveBase & smoothingObj, const std::string & probDesc, const double locationTol, const double qualityTol)
  {
    if (0 == stk::parallel_machine_rank(mComm)) std::cout << probDesc << std::endl;
    smoothingObj.set_gradient_filter(mGradientFilter);
    smoothingObj.set_solution_projector(mSolutionProjector);
    smooth_mesh_by_optimization(get_mesh(), get_coordinates_field(), optimizer, smoothingObj);

    EXPECT_GT(compute_mean_ratio_quality(), qualityTol) << probDesc;
    expect_nodes_on_circle(locationTol);
  }

  void reset_coords_and_test_mean_ratio_mesh_smoothing(const OptimizerFunction & optimizer, const std::string & probDesc, const double locationTol, const double qualityTol)
  {
    reset_coordinates();

    MeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), meanRatioObj);
    test_mesh_smoothing(optimizer, smoothingObj, probDesc, locationTol, qualityTol);
  }

  void reset_coords_and_test_kinematic_mesh_smoothing(const OptimizerFunction & optimizer,
    const KinematicElementObjective & elemObj,
    const std::string & probDesc,
    const double sizeChangeRatio,
    const double locationTol,
    const double qualityTol)
  {
    reset_coordinates();

    KinematicMeshObjective smoothingObj(get_mesh(), get_coordinates_field(), get_aux_meta().active_part(), elemObj);
    const double edgeLenForInscribingCircle = std::sqrt(3.0/2.0) * circle_radius();
    smoothingObj.set_constant_element_ref_size(sizeChangeRatio*edgeLenForInscribingCircle);
    test_mesh_smoothing(optimizer, smoothingObj, probDesc, locationTol, qualityTol);
  }

  GradientFilter build_node_search_direction_filter()
  {
    return build_node_search_direction_filter_for_sphere(sphere_center(), circle_radius());
  }

  SolutionProjector build_solution_projector()
  {
    return build_solution_projector_for_sphere(sphere_center(), circle_radius());
  }

  stk::math::Vector3d sphere_center() const { return stk::math::Vector3d(0.5,0.5,0.); }
  double circle_radius() const { return 1./std::sqrt(2.0); }
  GradientFilter mGradientFilter;
  SolutionProjector mSolutionProjector;
};

TEST_F(RightTriOnCircleSmoothing, smoothNodesSimultaneouslyUsingOptimizedMeanRatioOnCircle)
{
  const double distFromCircleTol = 1.e-9;

  const double baseQualityTol = 0.999;
  const double qualitySizeShapeTol = 0.93; // Due to stiffness?
  const double qualityWithSizeChangeTol= 0.66;
  const double sizeChangeRatio = 0.8;

  LogRedirecter log;

  if(stk::parallel_machine_size(mComm) == 1)
  {
    // ROL does not implement solution projection
    const double distFromCircleNoProjTol = 0.12;
    const double distFromCircleNoProjCorrectVolTol = 0.025;
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, "ROL L-BFGS using mean ratio", distFromCircleNoProjTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, "ROL L-BFGS using size-shape (matching elem vol)", 1.0, distFromCircleNoProjCorrectVolTol, qualitySizeShapeTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_ROL_lbfgs, sizeShapeObj, "ROL L-BFGS using size-shape (non-matching elem vol)", sizeChangeRatio, distFromCircleNoProjTol, qualityWithSizeChangeTol);
  }
  if(is_valid_proc_size_for_test())
  {
    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_lbfgs, "Krino L-BFGS using mean ratio", distFromCircleTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, "Krino L-BFGS using size-shape (matching elem vol)", 1.0, distFromCircleTol, qualitySizeShapeTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sethHillObj, "Krino L-BFGS using Seth-Hill (matching elem vol)", 1.0, distFromCircleTol, qualitySizeShapeTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_lbfgs, sizeShapeObj, "Krino L-BFGS using size-shape (non-matching elem vol)", sizeChangeRatio, distFromCircleTol, qualityWithSizeChangeTol);

    reset_coords_and_test_mean_ratio_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, "Krino Steepest descent using mean ratio", distFromCircleTol, baseQualityTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, "Krino Steepest descent using size-shape (matching elem vol)", 1.0, distFromCircleTol, qualitySizeShapeTol);
    reset_coords_and_test_kinematic_mesh_smoothing(optimize_mesh_using_krino_steepest_descent, sizeShapeObj, "Krino Steepest descent using size-shape (non-matching elem vol)", sizeChangeRatio, distFromCircleTol, qualityWithSizeChangeTol);
  }

  if (0 == stk::parallel_machine_rank(mComm))
    std::cout << log.get_log() << std::endl;
}

TEST_F(SquareOf8TrisSmoothing, nonUniformIsotropicSize_regressionTestForCenterNodeLocation)
{
  const std::vector<double> elemSizesGrowingAlongxy{ 0.1, 0.25, 0.5, 1.0, 1.0, 0.5, 0.25, 0.1 };

  const stk::math::Vector3d goldRegressionSethHillLocAlongxy(-0.7235881, -0.7235881, 0.); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sethHillObj, elemSizesGrowingAlongxy, goldRegressionSethHillLocAlongxy);

  const stk::math::Vector3d goldRegressionSizeShapeLocAlongxy(-0.08700252, -0.08700252, 0); // right direction, but very small magnitude. is it right?
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sizeShapeObj, elemSizesGrowingAlongxy, goldRegressionSizeShapeLocAlongxy);

  const stk::math::Vector3d goldRegressionModNeoHookeanLocAlongxy(-0.7755784, -0.7755784, 0.); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(modNeoHookean, elemSizesGrowingAlongxy, goldRegressionModNeoHookeanLocAlongxy);

  const std::vector<double> elemSizesGrowingAlongy{ 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 1.0, 0.1 };

  const stk::math::Vector3d goldRegressionSethHillLocAlongy(0,-0.89069324,0); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sethHillObj, elemSizesGrowingAlongy, goldRegressionSethHillLocAlongy);

  const stk::math::Vector3d goldRegressionSizeShapeLocAlongy(0,-0.1522262544,0); // right direction, but very small magnitude. is it right?
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sizeShapeObj, elemSizesGrowingAlongy, goldRegressionSizeShapeLocAlongy);

  const stk::math::Vector3d goldRegressionModNeoHookeanLocAlongy(0,-0.92212396,0); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(modNeoHookean, elemSizesGrowingAlongy, goldRegressionModNeoHookeanLocAlongy);
}

TEST_F(CubeOf24TetsSmoothing, nonUniformIsotropicSize_regressionTestForCenterNodeLocation)
{
  const std::vector<double> elemSizesGrowingAlongxyz{ 0.5, 0.5, 1., 1., 1., 0.5, 0.5, 1., 0.3, 0.1, 0.1, 0.3, 0.1, 0.1, 0.3, 0.3, 0.5, 1., 1., 0.5, 0.1, 0.1, 0.3, 0.3};

  const stk::math::Vector3d goldRegressionSethHillLocAlongxyz(-0.829575843, -0.829575843, -0.829575843); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sethHillObj, elemSizesGrowingAlongxyz, goldRegressionSethHillLocAlongxyz);

  const stk::math::Vector3d goldRegressionSizeShapeLocAlongxyz(0.192329888, 0.192329888, 0.192329888); // seems totally wrong! moved in the opposite direction
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sizeShapeObj, elemSizesGrowingAlongxyz, goldRegressionSizeShapeLocAlongxyz);

  const stk::math::Vector3d goldRegressionModNeoHookeanLocAlongxyz(-0.9546447, -0.9546447, -0.9546447); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(modNeoHookean, elemSizesGrowingAlongxyz, goldRegressionModNeoHookeanLocAlongxyz);

  const std::vector<double> elemSizesGrowingAlongz{ 0.5, 0.1, 0.5, 1.0, 0.5, 0.1, 0.5, 1.0, 0.5, 0.1, 0.5, 1.0, 0.5, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1};

  const stk::math::Vector3d goldRegressionSethHillLocAlongz(0,0,-0.953956); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sethHillObj, elemSizesGrowingAlongz, goldRegressionSethHillLocAlongz);

  const stk::math::Vector3d goldRegressionSizeShapeLocAlongz(0,0,0.2679499); // seems totally wrong! moved in the opposite direction
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(sizeShapeObj, elemSizesGrowingAlongz, goldRegressionSizeShapeLocAlongz);

  const stk::math::Vector3d goldRegressionModNeoHookeanLocAlongz(0,0,-0.9949982); // seems reasonable
  reset_coords_and_test_nonuniform_isotropic_smoothing_with_kinematic_objective(modNeoHookean, elemSizesGrowingAlongz, goldRegressionModNeoHookeanLocAlongz);
}

}



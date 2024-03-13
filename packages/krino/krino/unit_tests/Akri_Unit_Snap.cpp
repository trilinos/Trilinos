#include <type_traits>

#include <Akri_StkMeshFixture.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_Snap.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <Akri_SharpFeature.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_VolumePreservingSnappingLimiter.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_Unit_InterfaceGeometry.hpp>

namespace krino {

class RegularTriWithSides : public StkMeshTriFixture
{
protected:
  void create_sides_and_build_mesh(const std::vector<unsigned> &sidesetIds)
  {
    mBuilder.create_sideset_parts(sidesetIds);

    RegularTri meshSpec;
    build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
  }

  stk::mesh::Entity get_side_1() const { return mBuilder.get_side_with_nodes({get_assigned_node_for_index(0), get_assigned_node_for_index(1)}); }
  stk::mesh::Entity get_side_2() const { return mBuilder.get_side_with_nodes({get_assigned_node_for_index(0), get_assigned_node_for_index(2)}); }

  void expect_which_snaps_are_allowed(const std::vector<bool> & goldWhichSnapsAreAllowed, const std::vector<unsigned> & intersectionNodeIndices)
  {
    std::vector<stk::mesh::Entity> intersectionNodes;
    for (auto intersectionNodeIndex : intersectionNodeIndices)
      intersectionNodes.push_back(get_assigned_node_for_index(intersectionNodeIndex));

    const std::vector<bool> whichSnapsAreAllowed = which_intersection_point_nodes_are_compatible_for_snapping(mMesh, mBuilder.get_aux_meta(), mBuilder.get_phase_support(), intersectionNodes);
    EXPECT_EQ(goldWhichSnapsAreAllowed, whichSnapsAreAllowed);
  }
};

TEST_F(RegularTriWithSides, triMeshWithNoSidesets_attemptSnapToIntPointOnSide_snapsAllowed)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    create_sides_and_build_mesh({});

    expect_which_snaps_are_allowed({true,true}, {1,2});
  }
}

TEST_F(RegularTriWithSides, triMeshWithOneSidesetOnOneSide_attemptSnapToIntPointOnThirdSide_oneSnapAllowed)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    const unsigned sideset1Id = 1;
    create_sides_and_build_mesh({sideset1Id});

    mBuilder.add_sides_to_sidesets({get_side_1()}, {{sideset1Id}});

    expect_which_snaps_are_allowed({false,true}, {1,2});
  }
}

TEST_F(RegularTriWithSides, triMeshWithTwoSidesetOnTwoSides_attemptSnapToIntPointOnThirdSide_noSnapAllowed)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    const unsigned sideset1Id = 1;
    const unsigned sideset2Id = 2;
    create_sides_and_build_mesh({sideset1Id, sideset2Id});

    mBuilder.add_sides_to_sidesets({get_side_1(), get_side_2()}, {{sideset1Id},{sideset2Id}});

    expect_which_snaps_are_allowed({false,false}, {1,2});
  }
}

TEST_F(RegularTriWithSides, triMeshWithOneSidesetOnTwoSides_attemptSnapToIntPointOnThirdSide_noSnapAllowed)
{
  // This is the sideset keyhole problem.
  if(stk::parallel_machine_size(mComm) == 1)
  {
    const unsigned sideset1Id = 1;
    create_sides_and_build_mesh({sideset1Id});

    mBuilder.add_sides_to_sidesets({get_side_1(), get_side_2()}, {{sideset1Id},{sideset1Id}});

    expect_which_snaps_are_allowed({false,false}, {1,2});
  }
}

TEST_F(RegularTriWithSides, triMeshWithOneSidesetOnTwoSides_attemptSnapToIntPointOnVolume_noSnapAllowed)
{
  // This is a volume intersection point version of the keyhole problem.
  if(stk::parallel_machine_size(mComm) == 1)
  {
    const unsigned sideset1Id = 1;
    create_sides_and_build_mesh({sideset1Id});

    mBuilder.add_sides_to_sidesets({get_side_1(), get_side_2()}, {{sideset1Id},{sideset1Id}});

    expect_which_snaps_are_allowed({false,false,false}, {0,1,2});
  }
}

template <typename MESHSPEC>
class SharpFeatureFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
protected:
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;

  void find_sharp_features()
  {
    mySharpFeatureInfo.find_sharp_features(mMesh, mMesh.mesh_meta_data().coordinate_field(), mMesh.mesh_meta_data().universal_part(), myCosFeatureAngle);
  }

  void build_mesh_and_find_sharp_features()
  {
    this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    find_sharp_features();
  }

  void test_is_node_pinned(const stk::mesh::Entity node, const bool goldIsNodePinned)
  {
    const SharpFeatureConstraint * constraint = mySharpFeatureInfo.get_constraint(node);
    if (goldIsNodePinned)
    {
      EXPECT_TRUE(constraint != nullptr && constraint->is_pinned());
    }
    else
    {
      EXPECT_TRUE(constraint == nullptr || !constraint->is_pinned());
    }
  }

  bool is_node_in_assigned_nodes_for_indices(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const std::vector<unsigned> & nodeIndices)
  {
    for (auto nodeIndex : nodeIndices)
      if (this->get_assigned_node_for_index(nodeIndex) == node)
        return true;
    return false;
  }

  void test_are_nodes_pinned(const std::vector<unsigned> & goldPinnedNodeIndices)
  {
    std::vector<stk::mesh::Entity> ownedNodes;
    stk::mesh::get_selected_entities( mMesh.mesh_meta_data().locally_owned_part(), mMesh.buckets( stk::topology::NODE_RANK ), ownedNodes );

    for (auto && node : ownedNodes)
    {
      const bool goldIsNodePinned = is_node_in_assigned_nodes_for_indices(mMesh, node, goldPinnedNodeIndices);
      test_is_node_pinned(node, goldIsNodePinned);
    }
  }

  void test_is_node_constrained_on_edge(const unsigned sharpEdgeNodeIndex, const std::array<unsigned,2> & goldSharpEdgeNodeNbrIndices)
  {
    stk::mesh::Entity sharpEdgeNode = this->get_assigned_node_for_index(sharpEdgeNodeIndex);
    if (mMesh.is_valid(sharpEdgeNode) && mMesh.parallel_owner_rank(sharpEdgeNode) == mMesh.parallel_rank())
    {
      const SharpFeatureConstraint * constraint = mySharpFeatureInfo.get_constraint(sharpEdgeNode);
      ASSERT_TRUE(constraint != nullptr && constraint->is_constrained_on_edge());
      const std::array<stk::mesh::Entity,2> sharpEdgeNodes = constraint->get_sharp_edge_nodes();
      for (auto goldSharpEdgeNodeNbrIndex : goldSharpEdgeNodeNbrIndices)
      {
        stk::mesh::Entity goldSharpEdgeNodeNbr = this->get_assigned_node_for_index(goldSharpEdgeNodeNbrIndex);
        EXPECT_TRUE(sharpEdgeNodes[0] == goldSharpEdgeNodeNbr || sharpEdgeNodes[1] == goldSharpEdgeNodeNbr);
      }
    }
  }

  void test_are_all_nodes_pinned()
  {
    std::vector<stk::mesh::Entity> ownedNodes;
    stk::mesh::get_selected_entities( mMesh.mesh_meta_data().locally_owned_part(), mMesh.buckets( stk::topology::NODE_RANK ), ownedNodes );

    for (auto && node : ownedNodes)
      test_is_node_pinned(node, true);
  }

  MESHSPEC meshSpec;
  double myCosFeatureAngle{std::cos(M_PI/180.*135.0)};
  SharpFeatureInfo mySharpFeatureInfo;
};

typedef SharpFeatureFixture<RegularTet> SharpFeatureRegularTetFixture;

TEST_F(SharpFeatureRegularTetFixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<RightTet> SharpFeatureRightTetFixture;

TEST_F(SharpFeatureRightTetFixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<RegularTri> SharpFeatureRegularTriFixture;

TEST_F(SharpFeatureRegularTriFixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<Tri306090> SharpFeatureTri306090Fixture;

TEST_F(SharpFeatureTri306090Fixture, meshWithAllNodesOnCorners_allNodesArePinned)
{
  if(stk::parallel_machine_size(mComm) == 1)
  {
    build_mesh_and_find_sharp_features();

    test_are_all_nodes_pinned();
  }
}

typedef SharpFeatureFixture<TwoTri306090> SharpFeatureTwoTri306090Fixture;

TEST_F(SharpFeatureTwoTri306090Fixture, meshWithCornerNodesAndUnconstrainedNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});

    find_sharp_features();

    const std::vector<unsigned> goldPinnedNodeIndices{1,2,3};
    test_are_nodes_pinned(goldPinnedNodeIndices);
  }
}

typedef SharpFeatureFixture<FourRightTets> SharpFeatureFourRightTetsFixture;

TEST_F(SharpFeatureFourRightTetsFixture, meshWithCornerNodesAndUnconstrainedNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 4)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,0,1,1});
    else if(stk::parallel_machine_size(mComm) == 3)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,2});
    else if(stk::parallel_machine_size(mComm) == 4)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1,1,1}, {0,1,2,3});

    find_sharp_features();

    const std::vector<unsigned> goldPinnedNodeIndices{1,2,3,4,5};
    test_are_nodes_pinned(goldPinnedNodeIndices);
  }
}

typedef SharpFeatureFixture<TwoRightTets> SharpFeatureTwoRightTetsFixture;

TEST_F(SharpFeatureTwoRightTetsFixture, meshWithCornerNodesAndEdgeNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});

    find_sharp_features();

    const std::vector<unsigned> goldPinnedNodeIndices{1,2,3,4};
    test_are_nodes_pinned(goldPinnedNodeIndices);

    test_is_node_constrained_on_edge(0, {{1,3}});
  }
}

template <typename MESHSPEC>
class VolumePreservingSnappingLimiterFixture : public StkMeshFixture<MESHSPEC::TOPOLOGY>
{
protected:
  using StkMeshFixture<MESHSPEC::TOPOLOGY>::mMesh;

  stk::math::Vector3d compute_snap_location(const std::vector<unsigned> & snapNodeIndices, const std::vector<double> & snapNodeWeights)
  {
    stk::math::Vector3d snapLocation = stk::math::Vector3d::ZERO;
    for (size_t i=0; i<snapNodeIndices.size(); ++i)
    {
      stk::mesh::Entity snapNode = this->get_assigned_node_for_index(snapNodeIndices[i]);
      const stk::math::Vector3d nodeLocation(field_data<double>(*mMesh.mesh_meta_data().coordinate_field(), snapNode), mMesh.mesh_meta_data().spatial_dimension());
      snapLocation += snapNodeWeights[i] * nodeLocation;
    }
    return snapLocation;
  }

  void test_is_snap_allowed_based_on_volume_change(const bool goldIsSnapAllowed, const unsigned nodeIndex, const std::vector<unsigned> & snapNodeIndices, const std::vector<double> & snapNodeWeights)
  {
    stk::mesh::Entity node = this->get_assigned_node_for_index(nodeIndex);
    if (mMesh.is_valid(node) && mMesh.parallel_owner_rank(node) == mMesh.parallel_rank())
    {
      EXPECT_EQ(goldIsSnapAllowed, myVolumePreservingSnappingLimiter->is_snap_allowed(node, compute_snap_location(snapNodeIndices, snapNodeWeights)));
    }
  }

  VolumePreservingSnappingLimiter::ElementToBlockConverter build_element_to_block_converter()
  {
    auto converter = [](const stk::mesh::BulkData & mesh, const stk::mesh::Entity elem)
    {
      for (auto && part : mesh.bucket(elem).supersets())
        if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && !stk::mesh::is_auto_declared_part(*part))
          return part;
      stk::mesh::Part * blockPart = nullptr;
      return blockPart;
    };
    return converter;
  }

  void setup_volume_preserving_snapping_limiter()
  {
    myVolumePreservingSnappingLimiter = std::make_unique<VolumePreservingSnappingLimiter>(mMesh, *mMesh.mesh_meta_data().coordinate_field(), build_element_to_block_converter(), myVolumeConservationTol);
  }

  MESHSPEC meshSpec;
  double myVolumeConservationTol{0.05};
  std::unique_ptr<VolumePreservingSnappingLimiter> myVolumePreservingSnappingLimiter;
};

typedef VolumePreservingSnappingLimiterFixture<TwoRightTets> VolumePreservingSnappingLimiterTwoRightTetsFixture;

TEST_F(VolumePreservingSnappingLimiterTwoRightTetsFixture, meshWithOneBlockWithCornerNodesAndEdgeNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});

    setup_volume_preserving_snapping_limiter();

    test_is_snap_allowed_based_on_volume_change(true, 1, {1,2}, {0.99,0.01});
    test_is_snap_allowed_based_on_volume_change(false, 1, {1,2}, {0.5,0.5});

    test_is_snap_allowed_based_on_volume_change(true, 0, {0,1}, {0.5,0.5});
    test_is_snap_allowed_based_on_volume_change(true, 0, {0,3}, {0.5,0.5});
  }
}

TEST_F(VolumePreservingSnappingLimiterTwoRightTetsFixture, meshWithTwoBlockWithCornerNodesAndEdgeNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,0});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,1});

    setup_volume_preserving_snapping_limiter();

    test_is_snap_allowed_based_on_volume_change(false, 0, {0,1}, {0.5,0.5});
    test_is_snap_allowed_based_on_volume_change(false, 0, {0,3}, {0.5,0.5});
  }
}

typedef VolumePreservingSnappingLimiterFixture<TwoRightTris> VolumePreservingSnappingLimiterTwoRightTrisFixture;

TEST_F(VolumePreservingSnappingLimiterTwoRightTrisFixture, meshWithOneBlockWithCornerNodesAndEdgeNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,1}, {0,1});

    setup_volume_preserving_snapping_limiter();

    test_is_snap_allowed_based_on_volume_change(true, 1, {1,2}, {0.99,0.01});
    test_is_snap_allowed_based_on_volume_change(false, 1, {1,2}, {0.5,0.5});

    test_is_snap_allowed_based_on_volume_change(true, 0, {0,1}, {0.5,0.5});
    test_is_snap_allowed_based_on_volume_change(true, 0, {0,3}, {0.5,0.5});
  }
}

TEST_F(VolumePreservingSnappingLimiterTwoRightTrisFixture, meshWithTwoBlockWithCornerNodesAndEdgeNode_constraintsAreCorrect)
{
  if(stk::parallel_machine_size(mComm) <= 2)
  {
    if(stk::parallel_machine_size(mComm) == 1)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,0});
    else if(stk::parallel_machine_size(mComm) == 2)
      this->build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1,2}, {0,1});

    setup_volume_preserving_snapping_limiter();

    test_is_snap_allowed_based_on_volume_change(false, 0, {0,1}, {0.5,0.5});
    test_is_snap_allowed_based_on_volume_change(false, 0, {0,3}, {0.5,0.5});
  }
}

class SnapTrisFixture : public StkMeshTriFixture
{
public:
  SnapTrisFixture() {}
protected:
  double estimate_cut_quality_for_node(const std::vector<IntersectionPoint> & intersectionPoints, const mapFromEntityToIntPtIndexAndSnapAllowed & nodeToIntPtIndicesAndWhichSnapsAllowed, const stk::mesh::EntityId nodeId)
  {
    FieldRef coordsField = mMesh.mesh_meta_data().coordinate_field();
    const double minIntPtWeightForEstimatingCutQuality = 1.e-6;
    const bool globalIDsAreParallelConsistent = true;
    const auto domainsToNodesToQuality = determine_quality_per_node_per_domain(mMesh, mMesh.mesh_meta_data().universal_part(), coordsField, intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed, qualityMetric, minIntPtWeightForEstimatingCutQuality, globalIDsAreParallelConsistent);

    STK_ThrowRequire(1u == domainsToNodesToQuality.size());
    return domainsToNodesToQuality.begin()->second.at(nodeId);
  }

  std::vector<IntersectionPoint> get_edge_intersection_points(const std::vector<double> & nodeLSValues)
  {
    IntersectionPointFromNodalLevelsetInterfaceGeometry geometry;
    geometry.set_nodal_levelset(mMesh, get_assigned_node_global_ids(), nodeLSValues);
    NodeToCapturedDomainsMap nodesToCapturedDomains;
    return build_all_intersection_points(mMesh, mMesh.mesh_meta_data().universal_part(), geometry, nodesToCapturedDomains);
  }

  std::vector<size_t> get_allowed_snap_indices(const std::vector<IntersectionPoint> & intersectionPoints, const stk::mesh::Entity snapNode)
  {
    const SharpFeatureInfo * sharpFeatureInfo = nullptr;
    const auto nodeToIntPtIndicesAndWhichSnapsAllowed = get_node_to_intersection_point_indices_and_which_snaps_allowed(mMesh, sharpFeatureInfo, maxSnapForEdges, intersectionPoints);

    const FieldRef coordsField = mMesh.mesh_meta_data().coordinate_field();
    const double minQualityThatIsForSureNotInverted = 0.;
    const double cutQualityEstimate = estimate_cut_quality_for_node(intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed, mMesh.identifier(snapNode));
    const NodeToCapturedDomainsMap nodesToCapturedDomains;
    const stk::mesh::Selector elementSelector = mMesh.mesh_meta_data().universal_part();

    const auto nodeIntersectionPointIndicesAndWhichSnapsAllowed = nodeToIntPtIndicesAndWhichSnapsAllowed.at(snapNode);

    std::vector<size_t> allowedSnapIndices;
    for (auto && intPtIndexAndIsSnapAllowed : nodeIntersectionPointIndicesAndWhichSnapsAllowed)
    {
      const size_t intPtIndex = intPtIndexAndIsSnapAllowed.first;
      const bool isSnapAllowed = intPtIndexAndIsSnapAllowed.second;
      const IntersectionPoint & intersectionPoint = intersectionPoints[intPtIndex];

      if (isSnapAllowed && domains_already_snapped_to_node_are_also_at_intersection_point(nodesToCapturedDomains, snapNode, intersectionPoint.get_sorted_domains()))
      {
        const stk::math::Vector3d snapLocation = compute_intersection_point_location(mMesh.mesh_meta_data().spatial_dimension(), coordsField, intersectionPoint);
        const double minAcceptableQuality = std::max(minQualityThatIsForSureNotInverted, cutQualityEstimate);

        const double postSnapQuality = compute_quality_if_node_is_snapped_terminating_early_if_below_threshold(mMesh, elementSelector, coordsField, snapNode, snapLocation, qualityMetric, minAcceptableQuality);
        if (qualityMetric.is_first_quality_metric_better_than_second(postSnapQuality, minAcceptableQuality))
        {
          allowedSnapIndices.push_back(intPtIndex);
          std::cout << "Allowing snap of " << mMesh.identifier(snapNode) << " to " << snapLocation << " at " << debug_output(mMesh, intersectionPoint) << " with snap quality at or below " << postSnapQuality << " and estimated cut quality " << cutQualityEstimate << std::endl;
        }
        else
        {
          std::cout << "Skipping snap of " << mMesh.identifier(snapNode) << " to " << snapLocation << " at " << debug_output(mMesh, intersectionPoint) << " with snap quality at or below " << postSnapQuality << " and estimated cut quality " << cutQualityEstimate << std::endl;
        }
      }
    }

    return allowedSnapIndices;
  }

  bool has_matching_gold_snap(const IntersectionPoint & intPt, const unsigned snapNodeIndex, const std::vector<unsigned> & goldValidSnapEdgeNodeIndices)
  {
    const auto & intPtNodes = intPt.get_nodes();
    STK_ThrowRequire(2u == intPtNodes.size());
    const stk::mesh::Entity snapNode = get_assigned_node_for_index(snapNodeIndex);
    for (unsigned goldValidSnapEdgeNodeIndex : goldValidSnapEdgeNodeIndices)
      if ((intPtNodes[0] == snapNode && intPtNodes[1] == get_assigned_node_for_index(goldValidSnapEdgeNodeIndex)) ||
          (intPtNodes[1] == snapNode && intPtNodes[0] == get_assigned_node_for_index(goldValidSnapEdgeNodeIndex)))
        return true;
    return false;
  }

  void build_mesh_and_test_valid_edge_snaps_for_given_level_set(const unsigned snapNodeIndex, const std::vector<unsigned> & goldValidSnapEdgeNodeIndices, const std::vector<double> & nodeLSValues)
  {
    if(stk::parallel_machine_size(mComm) == 1)
    {
      PatchOfRegularTrisAroundNode meshSpec;
      build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});

      const auto intersectionPoints = get_edge_intersection_points(nodeLSValues);
      const std::vector<size_t> allowedSnapIndices = get_allowed_snap_indices(intersectionPoints, get_assigned_node_for_index(snapNodeIndex));

      EXPECT_EQ(goldValidSnapEdgeNodeIndices.size(), allowedSnapIndices.size());

      for (size_t intPtIndex : allowedSnapIndices)
      {
        EXPECT_TRUE(has_matching_gold_snap(intersectionPoints[intPtIndex], snapNodeIndex, goldValidSnapEdgeNodeIndices));
      }
    }
  }

  double maxSnapForEdges{1.0};
  const ScaledJacobianQualityMetric qualityMetric;
};


TEST_F(SnapTrisFixture, noSnapsAllowed)
{
  const std::vector<double> nodeLSValues = {1.1, -1., -1., -1., 1., 1., -1.};
  const std::vector<unsigned> goldValidSnapEdgeNodeIndices{};
  build_mesh_and_test_valid_edge_snaps_for_given_level_set(0, goldValidSnapEdgeNodeIndices, nodeLSValues);
}

TEST_F(SnapTrisFixture, allInternalEdgeSnapsAllowed)
{
  const std::vector<double> nodeLSValues = {0.9, -1., -1., -1., 1., 1., -1.};
  const std::vector<unsigned> goldValidSnapEdgeNodeIndices{ 1, 2, 3, 6 };
  build_mesh_and_test_valid_edge_snaps_for_given_level_set(0, goldValidSnapEdgeNodeIndices, nodeLSValues);
}

TEST_F(SnapTrisFixture, noSnapsAllowed_evenWithTerminatingIntersectionPointNearNode)
{
  const std::vector<double> nodeLSValues = {1.1, -1., -1., -1., 1., 1., -1.e-1};
  const std::vector<unsigned> goldValidSnapEdgeNodeIndices{};
  build_mesh_and_test_valid_edge_snaps_for_given_level_set(0, goldValidSnapEdgeNodeIndices, nodeLSValues);
}

TEST_F(SnapTrisFixture, noSnapsAllowed_interfaceJustInsideSupport)
{
  const std::vector<double> nodeLSValues = {1., 1., -2.e-2, -1.e-2, 1., 2., 2.};
  const std::vector<unsigned> goldValidSnapEdgeNodeIndices{};
  build_mesh_and_test_valid_edge_snaps_for_given_level_set(0, goldValidSnapEdgeNodeIndices, nodeLSValues);
}

TEST_F(SnapTrisFixture, becauseOfSnapLimit_onlySnapAllowedIsNearbyOne)
{
  maxSnapForEdges = 0.5;
  const std::vector<double> nodeLSValues = {1.0, -1.e-1, -1.e2, -1.e-1, 1., 1., -1.e-1};
  const std::vector<unsigned> goldValidSnapEdgeNodeIndices{2};
  build_mesh_and_test_valid_edge_snaps_for_given_level_set(0, goldValidSnapEdgeNodeIndices, nodeLSValues);
}
}

#include <type_traits>

#include <Akri_StkMeshFixture.hpp>
#include <Akri_FieldRef.hpp>
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
#include <Akri_MeshQuality.hpp>
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
      const stk::math::Vector3d nodeLocation(mMesh.mesh_meta_data().coordinate_field()->template data<double>().entity_values(snapNode), mMesh.mesh_meta_data().spatial_dimension());
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
    this->build_mesh_from_spec(meshSpec);

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
    this->build_mesh_from_spec_in_blocks(meshSpec, {1,2});

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
    this->build_mesh_from_spec(meshSpec);

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
    this->build_mesh_from_spec_in_blocks(meshSpec, {1,2});

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
    const CutQualityEstimator cutQualityEstimator(mMesh,
        mMesh.mesh_meta_data().universal_part(),
        coordsField,
        qualityMetric,
        minIntPtWeightForEstimatingCutQuality,
        globalIDsAreParallelConsistent);
    const auto domainsToNodesToQuality = cutQualityEstimator.estimate_quality_per_node_per_domain(
        intersectionPoints, nodeToIntPtIndicesAndWhichSnapsAllowed);

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
    const auto nodeToIntPtIndicesAndWhichSnapsAllowed = get_node_to_intersection_point_indices_and_which_snaps_allowed(mMesh, sharpFeatureInfo, maxSnapForEdges, minWeightAtNode, intersectionPoints);

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

        const double postSnapQuality = compute_quality_if_node_is_moved_terminating_early_if_below_threshold(mMesh, elementSelector, coordsField, snapNode, snapLocation, qualityMetric, minAcceptableQuality);
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
  double minWeightAtNode{0.0};
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

// A node may snap only to the points in its sector, where its weight exceeds (1-maxEdgeSnap)/maxEdgeSnap
// times the weight of every other node of the point.  With maxEdgeSnap 0.5, node 0 may snap to a point
// in the triangle where its weight is the largest, even where it would move more than half way to the
// opposite edge, but not to the others.
TEST_F(SnapTrisFixture, pointInsideTri_snapAllowedOnlyInSectorOfNode)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    PatchOfRegularTrisAroundNode meshSpec;
    build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});

    const std::vector<stk::mesh::Entity> triNodes{get_assigned_node_for_index(0),
        get_assigned_node_for_index(1),
        get_assigned_node_for_index(2)};
    const std::vector<IntersectionPoint> intersectionPoints{
        IntersectionPoint(true, triNodes, {0.04, 0.48, 0.48}, {0}),
        IntersectionPoint(true, triNodes, {0.3, 0.35, 0.35}, {0}),
        IntersectionPoint(true, triNodes, {0.4, 0.3, 0.3}, {0})};

    const std::vector<std::pair<size_t, bool>> allAllowed{{0, true}, {1, true}, {2, true}};
    EXPECT_EQ(allAllowed,
        get_node_to_intersection_point_indices_and_which_snaps_allowed(
            mMesh, nullptr, 1.0, 0.0, intersectionPoints)
            .at(triNodes[0]));

    const std::vector<std::pair<size_t, bool>> onlyThirdAllowed{{0, false}, {1, false}, {2, true}};
    EXPECT_EQ(onlyThirdAllowed,
        get_node_to_intersection_point_indices_and_which_snaps_allowed(
            mMesh, nullptr, 0.5, 0.0, intersectionPoints)
            .at(triNodes[0]));
  }
}

// Whatever the sector, a node may not snap to a point where its weight is at most the minimum weight
// (the edge tolerance): the point is within the tolerance of the opposite edge.
TEST_F(SnapTrisFixture, pointInsideTri_snapAllowedOnlyAboveMinimumWeightAtNode)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    PatchOfRegularTrisAroundNode meshSpec;
    build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});

    const std::vector<stk::mesh::Entity> triNodes{get_assigned_node_for_index(0),
        get_assigned_node_for_index(1),
        get_assigned_node_for_index(2)};
    const std::vector<IntersectionPoint> intersectionPoints{
        IntersectionPoint(true, triNodes, {0.04, 0.48, 0.48}, {0}),
        IntersectionPoint(true, triNodes, {0.4, 0.3, 0.3}, {0})};

    const std::vector<std::pair<size_t, bool>> onlySecondAllowed{{0, false}, {1, true}};
    EXPECT_EQ(onlySecondAllowed,
        get_node_to_intersection_point_indices_and_which_snaps_allowed(
            mMesh, nullptr, 1.0, 0.05, intersectionPoints)
            .at(triNodes[0]));
    EXPECT_EQ(onlySecondAllowed,
        get_node_to_intersection_point_indices_and_which_snaps_allowed(
            mMesh, nullptr, 0.5, 0.05, intersectionPoints)
            .at(triNodes[0]));
  }
}

// The sectors change continuously as a point moves off an edge: the interior node 0 may snap to a
// point just inside the triangle, beside a point on its edge 0-1, exactly when it may snap to that
// point.
TEST_F(SnapTrisFixture, pointBesideEdgeOfTri_sameSnapsAllowedAsOnEdge)
{
  if (stk::parallel_machine_size(mComm) == 1)
  {
    PatchOfRegularTrisAroundNode meshSpec;
    build_mesh(meshSpec.nodeLocs, {meshSpec.allElementConn});

    const std::vector<stk::mesh::Entity> triNodes{get_assigned_node_for_index(0),
        get_assigned_node_for_index(1),
        get_assigned_node_for_index(2)};
    const std::vector<stk::mesh::Entity> edgeNodes{triNodes[0], triNodes[1]};
    const std::vector<IntersectionPoint> intersectionPoints{
        IntersectionPoint(true, edgeNodes, {0.55, 0.45}, {0}),
        IntersectionPoint(true, triNodes, {0.55, 0.45 - 1.e-12, 1.e-12}, {0}),
        IntersectionPoint(true, edgeNodes, {0.45, 0.55}, {0}),
        IntersectionPoint(true, triNodes, {0.45, 0.55 - 1.e-12, 1.e-12}, {0})};

    const std::vector<std::pair<size_t, bool>> onlyFirstPairAllowed{{0, true}, {1, true}, {2, false}, {3, false}};
    EXPECT_EQ(onlyFirstPairAllowed,
        get_node_to_intersection_point_indices_and_which_snaps_allowed(
            mMesh, nullptr, 0.5, 0.0, intersectionPoints)
            .at(triNodes[0]));
  }
}
}

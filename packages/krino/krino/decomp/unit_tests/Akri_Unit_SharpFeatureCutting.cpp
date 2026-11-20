// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh.hpp>
#include <Akri_CreateInterfaceGeometry.hpp>
#include <Akri_Element.hpp>
#include <Akri_Element_Intersections.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MultipleEdgeCrossings.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SubElement.hpp>
#include <Akri_Unit_DecompositionFixture.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace krino
{

namespace
{

// procSizes are the processor counts the tests run on; empty for 1 up to the number of elements,
// each element on its own processor.
template <typename MESHSPEC, typename LS_FIELD_POLICY = LSPerInterfacePolicy, unsigned NUM_LS = 1>
class SharpFeatureCuttingFixture : public DecompositionFixture<MESHSPEC, LS_FIELD_POLICY, NUM_LS>
{
public:
  explicit SharpFeatureCuttingFixture(const std::vector<int> & procSizes)
  {
    if (!procSizes.empty()) this->set_valid_proc_sizes_for_test(procSizes);
    this->cdfem_support().set_cdfem_edge_degeneracy_handling(
        SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE);
    this->cdfem_support().register_cdfem_snap_displacements_field();
    this->build_mesh_from_spec(this->meshSpec);
    if (this->is_valid_proc_size_for_test()) this->setup_ls_fields();
  }

  // The steps of CDMesh::decompose_mesh before the sharp features are cut, except snapping, so that
  // the cuts are made on the mesh as built: the elements have their cutters, and those with an
  // interface their base subelement.
  void prepare_to_cut(const Surface & surface)
  {
    const Phase_Support & phaseSupport = Phase_Support::get(this->mMesh.mesh_meta_data());
    prepare_to_cut(std::make_unique<AnalyticSurfaceInterfaceGeometry>(
        std::vector<Surface_Identifier>{this->levelset_fields()[0].identifier},
        std::vector<const Surface *>{&surface},
        this->get_aux_meta().active_part(),
        this->cdfem_support(),
        phaseSupport));
  }

  void prepare_to_cut_level_sets()
  {
    const stk::mesh::MetaData & meta = this->mMesh.mesh_meta_data();
    prepare_to_cut(create_levelset_geometry(meta.spatial_dimension(),
        this->get_aux_meta().active_part(),
        this->cdfem_support(),
        Phase_Support::get(meta),
        this->levelset_fields()));
  }

  void prepare_to_cut(std::unique_ptr<InterfaceGeometry> geometry)
  {
    myGeometry = std::move(geometry);
    CDMesh & mesh = *this->cdmesh;
    fix_node_ownership_to_assure_selected_owned_element(this->mMesh, mesh.get_active_part());
    myGeometry->prepare_to_decompose_elements(this->mMesh, NodeToCapturedDomainsMap());
    mesh.generate_nonconformal_elements();
    mesh.set_phase_of_uncut_elements(*myGeometry);
    mesh.triangulate(*myGeometry);
  }

  // All of the sharp-feature cutting (CDMesh::cut_sharp_features): corners, rounds of sharp-edge
  // crossings and swaps, and multiple edge crossings.
  void cut_sharp_features() { this->cdmesh->cut_sharp_features(*myGeometry); }

  Mesh_Element & element(const size_t elemIndex)
  {
    Mesh_Element * elem = this->cdmesh->find_mesh_element(
        this->mMesh.identifier(this->get_assigned_element_for_index(elemIndex)));
    STK_ThrowRequire(elem);
    return *elem;
  }

  const SubElementNode * node(const size_t nodeIndex)
  {
    const stk::mesh::Entity entity = this->get_assigned_node_for_index(nodeIndex);
    for (auto && meshElem : this->cdmesh->elements)
      for (auto * elemNode : meshElem->get_nodes())
        if (elemNode->entity() == entity) return elemNode;
    return nullptr;
  }

  IntersectionPoint point(const std::vector<size_t> & nodeIndices,
      const std::vector<double> & weights,
      const std::vector<int> & sortedDomains)
  {
    std::vector<stk::mesh::Entity> nodes;
    for (size_t nodeIndex : nodeIndices)
      nodes.push_back(this->get_assigned_node_for_index(nodeIndex));
    return IntersectionPoint(true, nodes, weights, sortedDomains);
  }

  std::vector<const SubElement *> leaves(const size_t elemIndex)
  {
    std::vector<const SubElement *> elemLeaves;
    element(elemIndex).get_subelements(elemLeaves);
    return elemLeaves;
  }

  // The leaf nodes of the element.
  std::set<const SubElementNode *> leaf_nodes(const size_t elemIndex)
  {
    std::set<const SubElementNode *> nodes;
    for (auto * leaf : leaves(elemIndex))
      nodes.insert(leaf->get_nodes().begin(), leaf->get_nodes().end());
    return nodes;
  }

  // The child nodes on a face or edge of the element (all their mesh ancestors are nodes of it, but
  // not all of them) that are not nodes of its leaves.
  size_t num_hanging_nodes(const size_t elemIndex)
  {
    const NodeVec & elemNodes = element(elemIndex).get_nodes();
    const std::set<const SubElementNode *> elemNodeSet(elemNodes.begin(), elemNodes.end());
    const std::set<const SubElementNode *> nodesOfLeaves = leaf_nodes(elemIndex);
    size_t numHanging = 0;
    for (auto && meshNode : this->cdmesh->nodes)
    {
      if (meshNode->is_mesh_node()) continue;
      NodeSet ancestors;
      meshNode->get_ancestors(ancestors);
      const bool onFaceOrEdge = ancestors.size() < elemNodes.size() &&
          std::includes(elemNodeSet.begin(), elemNodeSet.end(), ancestors.begin(), ancestors.end());
      if (onFaceOrEdge && nodesOfLeaves.count(meshNode.get()) == 0) ++numHanging;
    }
    return numHanging;
  }

  // The leaf triangles of the element on its face of these mesh nodes, as sorted nodes.
  std::set<std::array<const SubElementNode *, 3>>
  leaf_triangles_on_face(const size_t elemIndex, const std::set<const SubElementNode *> & faceNodes)
  {
    std::set<std::array<const SubElementNode *, 3>> triangles;
    for (auto * leaf : leaves(elemIndex))
    {
      const NodeVec & leafNodes = leaf->get_nodes();
      for (int iOmit = 0; iOmit < 4; ++iOmit)
      {
        std::array<const SubElementNode *, 3> triangle;
        bool onFace = true;
        for (int i = 0, j = 0; i < 4; ++i)
        {
          if (i == iOmit) continue;
          NodeSet ancestors;
          leafNodes[i]->get_ancestors(ancestors);
          onFace = onFace &&
              std::includes(faceNodes.begin(), faceNodes.end(), ancestors.begin(), ancestors.end());
          triangle[j++] = leafNodes[i];
        }
        std::sort(triangle.begin(), triangle.end());
        if (onFace) triangles.insert(triangle);
      }
    }
    return triangles;
  }

  // The centroids of the leaf triangles of the element on its face of these mesh nodes, sorted, as
  // one vector (3 coordinates each); empty if the element is not owned here.
  std::vector<double>
  face_split_centroids(const size_t elemIndex, const std::vector<size_t> & faceNodeIndices)
  {
    std::vector<double> centroids;
    if (!this->is_owned_element(elemIndex)) return centroids;
    std::set<const SubElementNode *> faceNodes;
    for (size_t nodeIndex : faceNodeIndices)
      faceNodes.insert(node(nodeIndex));
    std::vector<std::array<double, 3>> sortedCentroids;
    for (auto & triangle : leaf_triangles_on_face(elemIndex, faceNodes))
    {
      const stk::math::Vector3d centroid =
          (triangle[0]->coordinates() + triangle[1]->coordinates() + triangle[2]->coordinates()) /
          3.;
      sortedCentroids.push_back({{centroid[0], centroid[1], centroid[2]}});
    }
    std::sort(sortedCentroids.begin(), sortedCentroids.end());
    for (auto & centroid : sortedCentroids)
      centroids.insert(centroids.end(), centroid.begin(), centroid.end());
    return centroids;
  }

  // Collective.  The face split of the element, as on the processor that owns it.
  std::vector<double>
  global_face_split_centroids(const size_t elemIndex, const std::vector<size_t> & faceNodeIndices)
  {
    std::vector<double> local = face_split_centroids(elemIndex, faceNodeIndices);
    const size_t localSize = local.size();
    size_t globalSize = 0;
    stk::all_reduce_max(this->mMesh.parallel(), &localSize, &globalSize, 1);
    local.resize(globalSize, 0.);
    std::vector<double> global(globalSize, 0.);
    stk::all_reduce_sum(this->mMesh.parallel(), local.data(), global.data(), globalSize);
    return global;
  }

  // Collective.  The two elements, wherever they are owned, split their face of these nodes into
  // the same leaf triangles.
  void expect_face_split_alike(
      const size_t elemIndex0, const size_t elemIndex1, const std::vector<size_t> & faceNodeIndices)
  {
    const std::vector<double> split0 = global_face_split_centroids(elemIndex0, faceNodeIndices);
    const std::vector<double> split1 = global_face_split_centroids(elemIndex1, faceNodeIndices);
    ASSERT_EQ(split0.size(), split1.size());
    EXPECT_FALSE(split0.empty());

    // Each centroid of one split is one of the other, to roundoff: centroids equal to roundoff can
    // be sorted differently, as the tets sum the nodes of a triangle in different orders.
    constexpr double tol = 1.e-12;
    const size_t numTriangles = split0.size() / 3;
    std::vector<bool> matched(numTriangles, false);
    for (size_t i = 0; i < numTriangles; ++i)
    {
      bool found = false;
      for (size_t j = 0; j < numTriangles && !found; ++j)
      {
        if (matched[j]) continue;
        found = std::abs(split0[3 * i] - split1[3 * j]) <= tol &&
            std::abs(split0[3 * i + 1] - split1[3 * j + 1]) <= tol &&
            std::abs(split0[3 * i + 2] - split1[3 * j + 2]) <= tol;
        if (found) matched[j] = true;
      }
      EXPECT_TRUE(found) << "triangle with centroid (" << split0[3 * i] << ", " << split0[3 * i + 1]
                         << ", " << split0[3 * i + 2] << ") of element " << elemIndex0
                         << " is not one of element " << elemIndex1;
    }
  }

  int spatial_dimension() const { return this->mMesh.mesh_meta_data().spatial_dimension(); }

  void cut_corners() { this->cdmesh->cut_sharp_feature_points(spatial_dimension()); }

  // In 2D the corners are all the sharp features.
  void cut_corners_then_edges()
  {
    this->cdmesh->cut_sharp_feature_points(spatial_dimension());
    if (spatial_dimension() == 3) this->cdmesh->cut_sharp_feature_points(2);
  }

  void cut_multiple_edge_crossings() { this->cdmesh->cut_multiple_edge_crossings_in_rounds(); }

  // The child nodes here with domains on the mesh edge of these nodes, within this distance of the
  // point.
  size_t num_child_nodes_with_domains_on_mesh_edge_near(const SubElementNode * edgeNode0,
      const SubElementNode * edgeNode1,
      const stk::math::Vector3d & point,
      const double distance)
  {
    size_t numNear = 0;
    for (auto && childNode : this->cdmesh->nodes)
    {
      if (childNode->is_mesh_node() || childNode->get_sorted_node_domains().empty()) continue;
      NodeSet ancestors;
      childNode->get_ancestors(ancestors);
      if (ancestors == NodeSet{edgeNode0, edgeNode1} &&
          (childNode->coordinates() - point).length() < distance)
        ++numNear;
    }
    return numNear;
  }

  // The child nodes here with domains within this distance of the node.
  size_t num_child_nodes_with_domains_near(const SubElementNode * meshNode, const double distance)
  {
    size_t numNear = 0;
    for (auto && childNode : this->cdmesh->nodes)
      if (!childNode->is_mesh_node() && !childNode->get_sorted_node_domains().empty() &&
          (childNode->coordinates() - meshNode->coordinates()).length() < distance)
        ++numNear;
    return numNear;
  }

  // Collective.  The number of multiple edge crossings the tets here would cut, over all processors
  // (each tet once, where it is owned).
  size_t global_num_multiple_edge_crossing_points()
  {
    size_t numPoints = 0;
    for (size_t elemIndex = 0; elemIndex < this->meshSpec.allElementConn.size(); ++elemIndex)
      if (this->is_owned_element(elemIndex))
        numPoints += element(elemIndex).get_multiple_edge_crossing_points(*this->cdmesh).size();
    size_t globalNumPoints = 0;
    stk::all_reduce_sum(this->mMesh.parallel(), &numPoints, &globalNumPoints, 1);
    return globalNumPoints;
  }

  // The child nodes here with sharp-feature domains.
  std::vector<const SubElementNode *> sharp_feature_child_nodes()
  {
    std::vector<const SubElementNode *> childNodes;
    for (auto && meshNode : this->cdmesh->nodes)
      if (!meshNode->is_mesh_node() && !meshNode->get_sorted_node_domains().empty())
        childNodes.push_back(meshNode.get());
    return childNodes;
  }

  // Level set i is -scales[i] at the nodes in phase i, and 0 at the others.
  void set_level_sets_by_node_phase(
      const std::vector<unsigned> & nodePhases, const std::vector<double> & scales)
  {
    for (size_t nodeIndex = 0; nodeIndex < nodePhases.size(); ++nodeIndex)
    {
      const stk::mesh::Entity meshNode = this->get_assigned_node_for_index(nodeIndex);
      if (!this->mMesh.is_valid(meshNode)) continue;
      for (unsigned ls = 0; ls < NUM_LS; ++ls)
        this->node_ls_field(this->get_ls_field(ls), meshNode) =
            (nodePhases[nodeIndex] == ls) ? -scales[ls] : 0.;
    }
  }

  // The internal faces of the element's leaves (each shared by two of them) that a sharp edge
  // crosses, uncaptured by the nodes of the face.
  size_t num_internal_leaf_faces_with_uncaptured_sharp_edge_crossings(const size_t elemIndex)
  {
    std::map<std::array<const SubElementNode *, 3>, NodeVec> facesNodes;
    std::map<std::array<const SubElementNode *, 3>, int> numLeavesWithFace;
    for (auto * leaf : leaves(elemIndex))
    {
      const NodeVec & leafNodes = leaf->get_nodes();
      for (int iOmit = 0; iOmit < 4; ++iOmit)
      {
        NodeVec faceNodes;
        for (int i = 0; i < 4; ++i)
          if (i != iOmit) faceNodes.push_back(leafNodes[i]);
        std::array<const SubElementNode *, 3> sortedFace{
            {faceNodes[0], faceNodes[1], faceNodes[2]}};
        std::sort(sortedFace.begin(), sortedFace.end());
        facesNodes[sortedFace] = faceNodes;
        ++numLeavesWithFace[sortedFace];
      }
    }
    size_t numCrossed = 0;
    for (auto & [sortedFace, faceNodes] : facesNodes)
      if (numLeavesWithFace[sortedFace] == 2 &&
          !element(elemIndex).get_uncaptured_sharp_edge_crossings_of_tet_face(faceNodes, this->cdfem_support().get_snapper().get_edge_tolerance()).empty())
        ++numCrossed;
    return numCrossed;
  }

  // The sum of the volumes of the element's leaves, and the smallest.
  std::pair<double, double> leaf_volumes(const size_t elemIndex)
  {
    double total = 0.;
    double smallest = std::numeric_limits<double>::max();
    for (auto * leaf : leaves(elemIndex))
    {
      const NodeVec & n = leaf->get_nodes();
      const double volume = Dot(Cross(n[1]->coordinates() - n[0]->coordinates(),
                                    n[2]->coordinates() - n[0]->coordinates()),
                                n[3]->coordinates() - n[0]->coordinates()) /
          6.;
      total += volume;
      smallest = std::min(smallest, volume);
    }
    return {total, smallest};
  }

  void expect_no_hanging_nodes()
  {
    for (size_t elemIndex = 0; elemIndex < this->meshSpec.allElementConn.size(); ++elemIndex)
    {
      if (this->is_owned_element(elemIndex))
      {
        EXPECT_EQ(0u, num_hanging_nodes(elemIndex)) << "element " << elemIndex;
      }
    }
  }

protected:
  std::unique_ptr<InterfaceGeometry> myGeometry;
};

// Inserting points into one element is checked in serial; cutting a rank, on as many
// processors as elements.
template <typename MESHSPEC>
class SerialSharpFeatureCuttingFixture : public SharpFeatureCuttingFixture<MESHSPEC>
{
public:
  SerialSharpFeatureCuttingFixture() : SharpFeatureCuttingFixture<MESHSPEC>({1}) {}
};

template <typename MESHSPEC, typename LS_FIELD_POLICY = LSPerInterfacePolicy, unsigned NUM_LS = 1>
class ParallelSharpFeatureCuttingFixture
    : public SharpFeatureCuttingFixture<MESHSPEC, LS_FIELD_POLICY, NUM_LS>
{
public:
  ParallelSharpFeatureCuttingFixture()
      : SharpFeatureCuttingFixture<MESHSPEC, LS_FIELD_POLICY, NUM_LS>({})
  {
  }
};

typedef SerialSharpFeatureCuttingFixture<RightTet> RightTetSharpFeatureCutting;
typedef SerialSharpFeatureCuttingFixture<TwoRightTets> TwoRightTetsSharpFeatureCutting;
typedef ParallelSharpFeatureCuttingFixture<TwoRightTets> TwoRightTetsSharpFeatureCorners;
typedef ParallelSharpFeatureCuttingFixture<FourRightTets> FourRightTetsSharpFeatureCorners;
typedef ParallelSharpFeatureCuttingFixture<FourRightTets> FourRightTetsSharpFeatureSharing;
typedef ParallelSharpFeatureCuttingFixture<TwoRightTets> TwoRightTetsSharpFeatureRepeats;
typedef SerialSharpFeatureCuttingFixture<RightTri> RightTriSharpFeatureCutting;
typedef ParallelSharpFeatureCuttingFixture<TwoRightTris> TwoRightTrisSharpFeatureCutting;
typedef ParallelSharpFeatureCuttingFixture<RightTri> RightTriSharpFeatureCorners;
typedef ParallelSharpFeatureCuttingFixture<RightTri> RightTriMultipleEdgeCrossings;
typedef ParallelSharpFeatureCuttingFixture<TwoRightTris, LSPerPhasePolicy, 3>
    TwoRightTris3LSPerPhaseSharpFeatureCorners;
typedef ParallelSharpFeatureCuttingFixture<RightTet> RightTetSharpFeatureEdges;
typedef ParallelSharpFeatureCuttingFixture<TwoRightTets> TwoRightTetsSharpFeatureEdges;
typedef ParallelSharpFeatureCuttingFixture<FourRightTets> FourRightTetsSharpFeatureEdges;
typedef ParallelSharpFeatureCuttingFixture<TwoRightTets, LSPerPhasePolicy, 3>
    TwoRightTets3LSPerPhaseSharpFeatureEdges;
typedef ParallelSharpFeatureCuttingFixture<RightTet, LSPerPhasePolicy, 4>
    RightTet4LSPerPhaseSharpFeatureEdges;

// A cuboid whose corner is at (0.2,0.2,0.2), inside the right tet and away from its nodes, so that
// the tet has an interface.
const Cuboid cuboidWithCornerInsideRightTet(
    stk::math::Vector3d(0.7, 0.7, 0.7), stk::math::Vector3d(1., 1., 1.));

const std::vector<int> cornerDomains{0, 0, 0};

// In 2D: a rectangle whose corner is at (0.2,0.2), inside the right tri and away from its nodes, so
// that the tri has an interface.  A corner of a curve in 2D has the domains of a sharp edge in 3D.
const Cuboid rectangleWithCornerInsideRightTri(
    stk::math::Vector3d(0.7, 0.7, 0.), stk::math::Vector3d(1., 1., 10.));

const std::vector<int> cornerDomains2D{0, 0};

} // namespace

TEST_F(RightTetSharpFeatureCutting, insertInteriorPoint_fourLeavesJoinedAtNewNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {0.25, 0.25, 0.25, 0.25}, cornerDomains));

    const std::vector<const SubElement *> elemLeaves = leaves(0);
    ASSERT_EQ(4u, elemLeaves.size());
    const SubElementNode * newNode =
        SubElementNode::common_child({node(0), node(1), node(2), node(3)});
    ASSERT_TRUE(newNode);
    EXPECT_EQ(cornerDomains, newNode->get_sorted_node_domains());
    for (auto * leaf : elemLeaves)
    {
      EXPECT_EQ(1, std::count(leaf->get_nodes().begin(), leaf->get_nodes().end(), newNode));
    }
  }
}

TEST_F(RightTetSharpFeatureCutting, insertFacePoint_threeLeavesAndFaceNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    element(0).insert_intersection_point(*cdmesh, point({1, 2, 3}, {0.2, 0.3, 0.5}, cornerDomains));

    EXPECT_EQ(3u, leaves(0).size());
    const SubElementNode * faceNode = SubElementNode::common_child({node(1), node(2), node(3)});
    ASSERT_TRUE(faceNode);
    EXPECT_EQ(cornerDomains, faceNode->get_sorted_node_domains());
    EXPECT_EQ(1u, leaf_nodes(0).count(faceNode));
    EXPECT_EQ(0u, num_hanging_nodes(0));
  }
}

TEST_F(RightTetSharpFeatureCutting, insertEdgePoint_twoLeavesAndEdgeNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    element(0).insert_intersection_point(*cdmesh, point({1, 2}, {0.4, 0.6}, cornerDomains));

    EXPECT_EQ(2u, leaves(0).size());
    const SubElementNode * edgeNode = SubElementNode::common_child({node(1), node(2)});
    ASSERT_TRUE(edgeNode);
    EXPECT_EQ(cornerDomains, edgeNode->get_sorted_node_domains());
    EXPECT_EQ(0u, num_hanging_nodes(0));
  }
}

TEST_F(RightTetSharpFeatureCutting, insertPointAtNode_noCutAndDomainsAddedCountingRepeats)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    element(0).insert_intersection_point(*cdmesh, point({0}, {1.}, {0}));
    element(0).insert_intersection_point(*cdmesh, point({0}, {1.}, cornerDomains));
    element(0).insert_intersection_point(*cdmesh, point({0}, {1.}, {0, 0}));

    EXPECT_EQ(1u, leaves(0).size());
    EXPECT_EQ(cornerDomains, node(0)->get_sorted_node_domains());
  }
}

// The second point is on the face, inside the tet, between the first point's node and nodes 0
// and 1: both leaves with that face are split at it.
TEST_F(RightTetSharpFeatureCutting, insertPointOnFaceBetweenLeaves_bothLeavesSplit)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {0.25, 0.25, 0.25, 0.25}, cornerDomains));
    const SubElementNode * interiorNode =
        SubElementNode::common_child({node(0), node(1), node(2), node(3)});
    ASSERT_TRUE(interiorNode);

    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {0.375, 0.375, 0.125, 0.125}, cornerDomains));

    EXPECT_EQ(8u, leaves(0).size());
    const SubElementNode * faceNode =
        SubElementNode::common_child({interiorNode, node(0), node(1)});
    ASSERT_TRUE(faceNode);
    EXPECT_EQ(1u, leaf_nodes(0).count(faceNode));
  }
}

// A point on the face shared by the tets, inserted by the first: the second does not have it yet.
TEST_F(TwoRightTetsSharpFeatureCutting, insertPointOnSharedFace_neighbourHasHangingNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    const size_t numNeighbourLeaves = leaves(1).size();
    element(0).insert_intersection_point(*cdmesh, point({0, 2, 4}, {0.2, 0.3, 0.5}, cornerDomains));

    EXPECT_EQ(3u, leaves(0).size());
    EXPECT_EQ(0u, num_hanging_nodes(0));
    EXPECT_EQ(numNeighbourLeaves, leaves(1).size());
    EXPECT_EQ(1u, num_hanging_nodes(1));
  }
}

// A corner at (0.0002, 0.3, 0.3), inside the first tet within the tolerance of the face it shares
// with the second (which the cuboid does not reach, so it has no subelements yet): a node on that
// face, and both tets split at it alike.
TEST_F(TwoRightTetsSharpFeatureCorners, cornerNearSharedFace_bothTetsSplitAtFaceNode)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0.5002, 0.8, 0.8), stk::math::Vector3d(1., 1., 1.));
    prepare_to_cut(cuboid);
    if (is_owned_element(1))
    {
      EXPECT_FALSE(element(1).have_subelements());
    }
    cut_corners();

    for (size_t elemIndex : {0u, 1u})
    {
      if (is_owned_element(elemIndex))
      {
        const SubElementNode * faceNode = SubElementNode::common_child({node(0), node(2), node(4)});
        EXPECT_TRUE(faceNode && cornerDomains == faceNode->get_sorted_node_domains());
        EXPECT_EQ(3u, leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// A corner at (0.2, 0.2, 0.2), inside the first tet: only it is split.
TEST_F(TwoRightTetsSharpFeatureCorners, cornerInsideTet_onlyThatTetSplit)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    cut_corners();

    if (is_owned_element(0))
    {
      EXPECT_EQ(4u, leaves(0).size());
    }
    if (is_owned_element(1))
    {
      EXPECT_FALSE(element(1).have_subelements());
    }
    expect_no_hanging_nodes();
  }
}

// A corner at (0.0002, 0.0002, 0.0002), within the tolerance of node 0: no cut, and node 0 is the
// corner, on every processor that has it.
TEST_F(TwoRightTetsSharpFeatureCorners, cornerNearNode_nodeCapturesCorner)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(
        stk::math::Vector3d(0.5002, 0.5002, 0.5002), stk::math::Vector3d(1., 1., 1.));
    prepare_to_cut(cuboid);
    cut_corners();

    if (is_owned_element(0))
    {
      EXPECT_EQ(1u, leaves(0).size());
    }
    if (is_owned_element(1))
    {
      EXPECT_FALSE(element(1).have_subelements());
    }
    ASSERT_TRUE(node(0));
    EXPECT_EQ(cornerDomains, node(0)->get_sorted_node_domains());
  }
}

// Two corners of the cuboid, at (0.0002, 0.1, 0.1) and (0.0002, 0.5, 0.1), both within the
// tolerance of the shared face: the second is on a face of a leaf made by the first, so its node
// has the first's as a parent, and the second tet is split at both in that order.
TEST_F(TwoRightTetsSharpFeatureCorners, twoCornersNearSharedFace_nestedNodesBothTetsSplitAlike)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0.5002, 0.3, 0.6), stk::math::Vector3d(1., 0.4, 1.));
    prepare_to_cut(cuboid);
    cut_corners();

    for (size_t elemIndex : {0u, 1u})
    {
      if (is_owned_element(elemIndex))
      {
        EXPECT_TRUE(SubElementNode::common_child({node(0), node(2), node(4)}));
        EXPECT_EQ(5u, leaves(elemIndex).size()) << "element " << elemIndex;
        EXPECT_EQ(5u, leaf_triangles_on_face(elemIndex, {node(0), node(2), node(4)}).size())
            << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// A corner at (0.0002, 0.0002, 0.4), inside the first tet within the tolerance of the edge from
// node 0 to node 5 shared by all four tets: a node on that edge, and every tet bisected at it.
TEST_F(FourRightTetsSharpFeatureCorners, cornerNearSharedEdge_allTetsBisected)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0.5002, 0.5002, 0.9), stk::math::Vector3d(1., 1., 1.));
    prepare_to_cut(cuboid);
    cut_corners();

    for (size_t elemIndex = 0; elemIndex < 4; ++elemIndex)
    {
      if (is_owned_element(elemIndex))
      {
        const SubElementNode * edgeNode = SubElementNode::common_child({node(0), node(5)});
        EXPECT_TRUE(edgeNode && cornerDomains == edgeNode->get_sorted_node_domains());
        EXPECT_EQ(2u, leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 5});
    expect_face_split_alike(1, 2, {0, 3, 5});
    expect_face_split_alike(2, 3, {0, 4, 5});
    expect_face_split_alike(3, 0, {0, 1, 5});
  }
}

// A sharp edge along x at y=z=0.3 crosses the face x=0 shared by the two tets, and an outer face of
// each: one node on the shared face, and both tets split at it alike.
TEST_F(TwoRightTetsSharpFeatureEdges, sharpEdgeCrossingSharedFace_bothTetsSplitAtFaceNode)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(0., 0.8, 0.8), stk::math::Vector3d(4., 1., 1.));
    prepare_to_cut(cuboid);
    cut_corners_then_edges();

    for (size_t elemIndex : {0u, 1u})
    {
      if (is_owned_element(elemIndex))
      {
        const SubElementNode * faceNode = SubElementNode::common_child({node(0), node(2), node(4)});
        EXPECT_TRUE(faceNode && (std::vector<int>{0, 0} == faceNode->get_sorted_node_domains()));
        EXPECT_EQ(5u, leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// A sharp edge along y at x=0.0003, z=0.4 crosses the face y=0 shared by the first and last tets
// within the tolerance of their edge from node 0 to node 5, and their outer faces through node 1
// within the tolerance of the edges from node 5 to nodes 2 and 4: three edge nodes, and every tet
// with one of those edges bisected at it, in the same order as the tet that made it.
TEST_F(FourRightTetsSharpFeatureEdges, sharpEdgeCrossingFacesNearTetEdges_edgeNodesThroughTheRing)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(1.00015, 0., 0.9), stk::math::Vector3d(1.9997, 4., 1.));
    prepare_to_cut(cuboid);
    cut_corners_then_edges();

    for (size_t elemIndex = 0; elemIndex < 4; ++elemIndex)
    {
      if (is_owned_element(elemIndex))
      {
        EXPECT_TRUE(SubElementNode::common_child({node(0), node(5)})) << "element " << elemIndex;
        EXPECT_EQ(3u, leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 5});
    expect_face_split_alike(1, 2, {0, 3, 5});
    expect_face_split_alike(2, 3, {0, 4, 5});
    expect_face_split_alike(3, 0, {0, 1, 5});
  }
}

// A cuboid whose corner is eta inside the face x+y+z=1 of the right tet (the face of nodes 1, 2,
// 3): the corner's weight for node 0 is eta, and its three edges cross that face eta from it.
static Cuboid cuboid_with_corner_inside_face_of_right_tet(const double eta)
{
  const double a = (1. - eta) / 3.;
  return Cuboid(stk::math::Vector3d(a + 0.5, a + 0.5, a + 0.5), stk::math::Vector3d(1., 1., 1.));
}

// Within the tolerance of the face, the corner is a node on it, which captures its edges' crossings
// of the face: no other node.
TEST_F(
    RightTetSharpFeatureEdges, cornerJustInsideFaceWithinTolerance_cornerNodeCapturesItsCrossings)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid = cuboid_with_corner_inside_face_of_right_tet(1.e-4);
    prepare_to_cut(cuboid);
    cut_corners_then_edges();

    const std::vector<const SubElementNode *> sharpNodes = sharp_feature_child_nodes();
    ASSERT_EQ(1u, sharpNodes.size());
    EXPECT_EQ(sharpNodes[0], SubElementNode::common_child({node(1), node(2), node(3)}));
    EXPECT_EQ(cornerDomains, sharpNodes[0]->get_sorted_node_domains());
    EXPECT_EQ(3u, leaves(0).size());
  }
}

// Beyond the tolerance, the corner is a node inside the tet, and its edges' crossings of the face
// are three more nodes on it.
TEST_F(
    RightTetSharpFeatureEdges, cornerJustInsideFaceBeyondTolerance_cornerNodeAndThreeCrossingNodes)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid = cuboid_with_corner_inside_face_of_right_tet(1.e-2);
    prepare_to_cut(cuboid);
    cut_corners_then_edges();

    const std::vector<const SubElementNode *> sharpNodes = sharp_feature_child_nodes();
    EXPECT_EQ(4u, sharpNodes.size());
    EXPECT_EQ(1,
        std::count_if(sharpNodes.begin(),
            sharpNodes.end(),
            [](const SubElementNode * sharpNode)
            { return sharpNode->get_sorted_node_domains().size() == 3; }));
    EXPECT_EQ(10u, leaves(0).size());
  }
}

// The kink case: a corner at (0.4, 0.5988, 6e-4), with the cuboid extending in +x, +y and -z, near
// the edge of nodes 1 and 2 on both of its faces (its weights for nodes 0 and 3 are both 6e-4). The
// corner is a node on that edge, and its node captures the crossings of its three edges, including
// the -z one 0.0012 from the edge on the face z=0: no other node, so no kink.
TEST_F(RightTetSharpFeatureEdges, cornerNearTetEdge_cornerNodeOnEdgeCapturesItsCrossings)
{
  if (is_valid_proc_size_for_test())
  {
    const double eta = 6.e-4;
    const double cz = 6.e-4;
    const double cx = 0.4;
    const double cy = 1. - eta - cz - cx;
    const Cuboid cuboid(
        stk::math::Vector3d(cx + 0.5, cy + 0.5, cz - 0.5), stk::math::Vector3d(1., 1., 1.));
    prepare_to_cut(cuboid);
    cut_corners_then_edges();

    const std::vector<const SubElementNode *> sharpNodes = sharp_feature_child_nodes();
    ASSERT_EQ(1u, sharpNodes.size());
    EXPECT_EQ(sharpNodes[0], SubElementNode::common_child({node(1), node(2)}));
    EXPECT_EQ(cornerDomains, sharpNodes[0]->get_sorted_node_domains());
    EXPECT_EQ(2u, leaves(0).size());
  }
}

// Level sets per phase: the three phases of the nodes of the shared face meet at its centroid, a
// node on it, with both tets split at it alike.
TEST_F(TwoRightTets3LSPerPhaseSharpFeatureEdges, triplePointOnSharedFace_bothTetsSplitAtFaceNode)
{
  if (is_valid_proc_size_for_test())
  {
    set_level_sets_by_node_phase({0, 0, 1, 0, 2}, {1., 1., 1.});
    prepare_to_cut_level_sets();
    cut_corners_then_edges();

    for (size_t elemIndex : {0u, 1u})
    {
      if (is_owned_element(elemIndex))
      {
        const SubElementNode * faceNode = SubElementNode::common_child({node(0), node(2), node(4)});
        EXPECT_TRUE(faceNode && (std::vector<int>{0, 1, 2} == faceNode->get_sorted_node_domains()));
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// Each node of the tet in its own phase, phase 0 scaled so that the quadruple point has weight 1e-4
// at node 0: it is a node on the face opposite node 0, at the point there where phases 1, 2 and 3
// meet (both at the centroid), so that point is no other node.
TEST_F(RightTet4LSPerPhaseSharpFeatureEdges,
    quadruplePointWithinToleranceOfFace_nodeOnFaceAtItsTriplePoint)
{
  if (is_valid_proc_size_for_test())
  {
    const double eta = 1.e-4;
    set_level_sets_by_node_phase({0, 1, 2, 3}, {(1. / eta - 1.) / 3., 1., 1., 1.});
    prepare_to_cut_level_sets();
    cut_corners_then_edges();

    const SubElementNode * quadrupleNode =
        SubElementNode::common_child({node(1), node(2), node(3)});
    ASSERT_TRUE(quadrupleNode);
    EXPECT_EQ((std::vector<int>{0, 1, 2, 3}), quadrupleNode->get_sorted_node_domains());
    for (auto * sharpNode : sharp_feature_child_nodes())
    {
      EXPECT_NE((std::vector<int>{1, 2, 3}), sharpNode->get_sorted_node_domains());
    }
  }
}

// The cdfem_hidden_interface tet, mapped onto the right tet: four planar level sets, each a phase
// where it is the smallest.  The line where phases 1, 2 and 3 meet runs straight through the tet, and
// after the first round of cuts it crosses an internal leaf edge exactly, within roundoff of the
// faces around that edge, whichever of them finds it.  The swap of that edge makes the line one
// leaf edge, between the nodes where it leaves the tet, with no node on it in between.
TEST_F(RightTet4LSPerPhaseSharpFeatureEdges, tripleLineCrossingInternalLeafEdge_recoveredBySwapOfTheEdge)
{
  if (is_valid_proc_size_for_test())
  {
    const std::vector<std::array<double, 4>> nodeLevelSets{
        {{2., -0.4, 0., 0.}}, {{2., -0.4, 1., -1.}}, {{2., -0.4, -1., 1.}}, {{-1., 2.6, 0., 0.}}};
    for (size_t nodeIndex = 0; nodeIndex < nodeLevelSets.size(); ++nodeIndex)
    {
      const stk::mesh::Entity meshNode = get_assigned_node_for_index(nodeIndex);
      if (!mMesh.is_valid(meshNode)) continue;
      for (unsigned ls = 0; ls < 4; ++ls)
        node_ls_field(get_ls_field(ls), meshNode) = nodeLevelSets[nodeIndex][ls];
    }
    prepare_to_cut_level_sets();
    cut_sharp_features();

    const std::vector<int> tripleLineDomains{1, 2, 3};
    std::vector<const SubElementNode *> tripleLineNodes;
    for (auto * leafNode : leaf_nodes(0))
      if (leafNode->get_sorted_node_domains() == tripleLineDomains) tripleLineNodes.push_back(leafNode);
    ASSERT_EQ(2u, tripleLineNodes.size());
    bool isLeafEdge = false;
    for (auto * leaf : leaves(0))
    {
      const NodeVec & leafNodes = leaf->get_nodes();
      if (std::count(leafNodes.begin(), leafNodes.end(), tripleLineNodes[0]) == 1 &&
          std::count(leafNodes.begin(), leafNodes.end(), tripleLineNodes[1]) == 1)
        isLeafEdge = true;
    }
    EXPECT_TRUE(isLeafEdge);
    const auto [volume, smallestVolume] = leaf_volumes(0);
    EXPECT_NEAR(1. / 6., volume, 1.e-12);
    EXPECT_LT(0., smallestVolume);
  }
}

// The first tet finds a point on the edge from node 0 to node 5, shared by all four tets, and the
// second a point on the edge from node 2 to node 5, which it shares with the first: every tet with
// either edge gets the points on it, wherever they were found, and the first two split their shared
// face of nodes 0, 2 and 5 at both in the same order.
TEST_F(FourRightTetsSharpFeatureSharing, pointsOnSharedEdgesFoundByDifferentTets_tetsSplitAlike)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    std::vector<IntersectionPoint> foundPoints;
    if (is_owned_element(0)) foundPoints.push_back(point({0, 5}, {0.6, 0.4}, {0, 0}));
    if (is_owned_element(1)) foundPoints.push_back(point({2, 5}, {0.7, 0.3}, {0, 0}));
    cdmesh->share_and_insert_intersection_points(foundPoints, 1.e-12);

    const std::array<size_t, 4> numLeaves{{3, 3, 2, 2}};
    for (size_t elemIndex = 0; elemIndex < 4; ++elemIndex)
    {
      if (is_owned_element(elemIndex))
      {
        EXPECT_EQ(numLeaves[elemIndex], leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 5});
    expect_face_split_alike(1, 2, {0, 3, 5});
    expect_face_split_alike(2, 3, {0, 4, 5});
    expect_face_split_alike(3, 0, {0, 1, 5});
  }
}

// A point at node 0 found by the first tet: node 0 gets its domains wherever it is.
TEST_F(FourRightTetsSharpFeatureSharing, pointAtSharedNode_domainsAddedEverywhere)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    std::vector<IntersectionPoint> foundPoints;
    if (is_owned_element(0)) foundPoints.push_back(point({0}, {1.}, cornerDomains));
    cdmesh->share_and_insert_intersection_points(foundPoints, 1.e-12);

    ASSERT_TRUE(node(0));
    EXPECT_EQ(cornerDomains, node(0)->get_sorted_node_domains());
  }
}

// A cuboid whose corner is outside the right tet, at (0.15, 0.45, 0.65), with the cuboid extending
// in -x, -y and -z: cutting the crossings of its edges leaves a segment of one crossing an internal
// leaf face between the leaves of the nodes at its ends, which a 2-3 swap makes an edge.
TEST_F(RightTetSharpFeatureEdges, segmentAcrossInternalLeafFace_recoveredBySwap)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid cuboid(stk::math::Vector3d(-0.35, -0.05, 0.15), stk::math::Vector3d(1., 1., 1.));
    prepare_to_cut(cuboid);
    cut_corners_then_edges();
    ASSERT_EQ(1u, num_internal_leaf_faces_with_uncaptured_sharp_edge_crossings(0));
    const size_t numLeavesBefore = leaves(0).size();

    element(0).recover_sharp_edges_in_tet_by_swaps(*cdmesh);
    EXPECT_EQ(0u, num_internal_leaf_faces_with_uncaptured_sharp_edge_crossings(0));
    EXPECT_EQ(numLeavesBefore + 1, leaves(0).size());
    const auto [volume, smallestVolume] = leaf_volumes(0);
    EXPECT_NEAR(1. / 6., volume, 1.e-12);
    EXPECT_LT(0., smallestVolume);
  }
}

// A sharp edge parallel to z through (0.3, 0.2), with the nodes at its crossings of the faces
// z = 0 and x+y+z = 1 inserted after a node inside the tet at (0.05, 0.35, 0.05), which is not on
// it: the segment between them crosses two internal leaf faces, which no 2-3 swap recovers, but an
// edge swap does, making it an edge.
class RightTetSegmentAroundInteriorNode : public RightTetSharpFeatureEdges
{
public:
  void insert_nodes()
  {
    prepare_to_cut(cuboid);
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {0.55, 0.05, 0.35, 0.05}, {}));
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {0.5, 0.3, 0.2, 0.}, sharpEdgeDomains));
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {0., 0.3, 0.2, 0.5}, sharpEdgeDomains));
  }
  const SubElementNode * bottom_node()
  {
    return SubElementNode::common_child({node(0), node(1), node(2)});
  }
  const SubElementNode * top_node()
  {
    return SubElementNode::common_child({node(1), node(2), node(3)});
  }
  bool is_leaf_edge(const SubElementNode * node0, const SubElementNode * node1)
  {
    for (auto * leaf : leaves(0))
    {
      const NodeVec & leafNodes = leaf->get_nodes();
      if (std::count(leafNodes.begin(), leafNodes.end(), node0) == 1 &&
          std::count(leafNodes.begin(), leafNodes.end(), node1) == 1)
        return true;
    }
    return false;
  }

  // The cutter refers to the surface.
  const Cuboid cuboid{stk::math::Vector3d(-0.2, -0.3, 0.), stk::math::Vector3d(1., 1., 3.)};
  const std::vector<int> sharpEdgeDomains{0, 0};
};

TEST_F(RightTetSegmentAroundInteriorNode, segmentAcrossTwoInternalLeafFaces_recoveredByEdgeSwap)
{
  if (is_valid_proc_size_for_test())
  {
    insert_nodes();
    ASSERT_TRUE(bottom_node() && top_node());
    ASSERT_FALSE(is_leaf_edge(bottom_node(), top_node()));
    ASSERT_EQ(2u, num_internal_leaf_faces_with_uncaptured_sharp_edge_crossings(0));
    ASSERT_EQ(8u, leaves(0).size());

    element(0).recover_sharp_edges_in_tet_by_swaps(*cdmesh);
    EXPECT_EQ(0u, num_internal_leaf_faces_with_uncaptured_sharp_edge_crossings(0));
    EXPECT_TRUE(is_leaf_edge(bottom_node(), top_node()));
    EXPECT_EQ(9u, leaves(0).size());
    const auto [volume, smallestVolume] = leaf_volumes(0);
    EXPECT_NEAR(1. / 6., volume, 1.e-12);
    EXPECT_LT(0., smallestVolume);
  }
}

// A node inside the tet whose leaf edge to node 3, at (0,0,1), passes 1e-4 beside the sharp edge
// parallel to z through (0.3, 0.2): within the edge tolerance of the leaf edge, but beyond roundoff.
class RightTetSharpEdgeBesideInternalLeafEdge : public RightTetSharpFeatureEdges
{
public:
  void insert_interior_node()
  {
    prepare_to_cut(cuboid);
    const stk::math::Vector3d top(0., 0., 1.);
    const stk::math::Vector3d besideEdge = crossingOfLeafEdge + distanceFromSharpEdge *
        stk::math::Vector3d(-0.2, 0.3, 0.).unit_vector();
    const stk::math::Vector3d x = top + 1.2 * (besideEdge - top);
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2, 3}, {1. - x[0] - x[1] - x[2], x[0], x[1], x[2]}, {}));
  }
  const SubElementNode * interior_node() { return SubElementNode::common_child({node(0), node(1), node(2), node(3)}); }

  // The child nodes with the sharp edge's domains within this distance of where it passes the leaf edge.
  std::vector<const SubElementNode *> sharp_edge_nodes_near_leaf_edge(const double distance)
  {
    std::vector<const SubElementNode *> nearNodes;
    for (auto * leafNode : leaf_nodes(0))
      if (leafNode->get_sorted_node_domains() == sharpEdgeDomains &&
          (leafNode->coordinates() - crossingOfLeafEdge).length() < distance)
        nearNodes.push_back(leafNode);
    return nearNodes;
  }

  const Cuboid cuboid{stk::math::Vector3d(-0.2, -0.3, 0.), stk::math::Vector3d(1., 1., 3.)};
  const std::vector<int> sharpEdgeDomains{0, 0};
  const stk::math::Vector3d crossingOfLeafEdge{0.3, 0.2, 0.25};
  const double distanceFromSharpEdge = 1.e-4;
};

// The sharp edge passes within the edge tolerance of an internal leaf edge that is not a chord of
// it: the faces around the leaf edge find it as crossings near each other, which are one point
// (cut_sharp_feature_points), so one node is cut, on the sharp edge.  It is a node of the face it
// was found on, beside the leaf edge, which makes slivers of the leaves at it; moving it onto the
// leaf edge instead moves it off the sharp edge, and made the sphere-and-binder decomposition
// slightly worse, so it is left on the sharp edge.
TEST_F(RightTetSharpEdgeBesideInternalLeafEdge, sharpEdgeBesideInternalLeafEdge_oneNodeOnTheSharpEdge)
{
  if (is_valid_proc_size_for_test())
  {
    insert_interior_node();
    ASSERT_TRUE(interior_node());
    cdmesh->cut_sharp_feature_points(2);

    const std::vector<const SubElementNode *> nearNodes = sharp_edge_nodes_near_leaf_edge(0.01);
    ASSERT_EQ(1u, nearNodes.size());
    EXPECT_NEAR(crossingOfLeafEdge[0], nearNodes[0]->coordinates()[0], 1.e-12);
    EXPECT_NEAR(crossingOfLeafEdge[1], nearNodes[0]->coordinates()[1], 1.e-12);
    const auto [volume, smallestVolume] = leaf_volumes(0);
    EXPECT_NEAR(1. / 6., volume, 1.e-12);
    EXPECT_LT(0., smallestVolume);
  }
}

// The sharp edge passes through the node at its crossing of the face z = 0, which has its domains:
// the faces of the leaves at that node have no crossing that it does not capture.
TEST_F(RightTetSegmentAroundInteriorNode, sharpEdgeThroughNodeWithItsDomains_capturedByTheNode)
{
  if (is_valid_proc_size_for_test())
  {
    insert_nodes();
    const SubElementNode * bottom = bottom_node();
    ASSERT_TRUE(bottom);
    size_t numFacesAtNode = 0;
    for (auto * leaf : leaves(0))
    {
      const NodeVec & leafNodes = leaf->get_nodes();
      for (int iOmit = 0; iOmit < 4; ++iOmit)
      {
        NodeVec faceNodes;
        for (int i = 0; i < 4; ++i)
          if (i != iOmit) faceNodes.push_back(leafNodes[i]);
        if (std::count(faceNodes.begin(), faceNodes.end(), bottom) == 0) continue;
        ++numFacesAtNode;
        EXPECT_TRUE(element(0).get_uncaptured_sharp_edge_crossings_of_tet_face(faceNodes, cdfem_support().get_snapper().get_edge_tolerance()).empty());
      }
    }
    EXPECT_LT(0u, numFacesAtNode);
  }
}

// A slab between the planes x = xLo and x = xHi (normal to x, and with its sharp edges far away):
// an edge of the mesh across it, with both ends outside, crosses it twice, which the signs at its
// ends do not show.
static Cuboid slab_normal_to_x(const double xLo, const double xHi)
{
  return Cuboid(
      stk::math::Vector3d(0.5 * (xLo + xHi), 0., 0.), stk::math::Vector3d(xHi - xLo, 10., 10.));
}

const std::vector<int> surfaceDomains{0};

// The edges of the right tet from node 1, to nodes 0, 2 and 3, cross the slab 0.4 <= x <= 0.6
// twice each, and no others do: six points, each on its mesh edge.
TEST_F(RightTetSharpFeatureEdges, edgesAcrossSlab_twoCrossingsOnEach)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid slab = slab_normal_to_x(0.4, 0.6);
    prepare_to_cut(slab);

    const std::vector<IntersectionPoint> points =
        element(0).get_multiple_edge_crossing_points(*cdmesh);
    ASSERT_EQ(6u, points.size());
    const stk::mesh::Entity node1 = get_assigned_node_for_index(1);
    for (auto & point : points)
    {
      EXPECT_EQ(surfaceDomains, point.get_sorted_domains());
      ASSERT_EQ(2u, point.get_nodes().size());
      const auto iter = std::find(point.get_nodes().begin(), point.get_nodes().end(), node1);
      ASSERT_TRUE(iter != point.get_nodes().end());
      const double weightOfNode1 =
          point.get_weights()[std::distance(point.get_nodes().begin(), iter)];
      EXPECT_TRUE(std::abs(weightOfNode1 - 0.4) < 1.e-12 || std::abs(weightOfNode1 - 0.6) < 1.e-12)
          << weightOfNode1;
    }
  }
}

// The slab is 1e-4 from node 0, within the tolerance of it along its edges: node 0 is on the
// surface, so it gets the surface's domains, and no node is made beside it.
TEST_F(RightTetSharpFeatureEdges, slabWithinToleranceOfNode_nodeGetsTheDomains)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid slab = slab_normal_to_x(1.e-4, 0.6);
    prepare_to_cut(slab);
    cut_multiple_edge_crossings();

    EXPECT_EQ(surfaceDomains, node(0)->get_sorted_node_domains());
    EXPECT_EQ(0u, num_child_nodes_with_domains_near(node(0), 0.1));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// The slab 0.2 from node 0 is far beyond the edge tolerance from it: both its faces get nodes on
// the edge 0-1, with snapping or only cutting alike.
static void expect_slab_beside_node_cut_on_edge01(RightTetSharpFeatureEdges & fixture)
{
  const Cuboid slab = slab_normal_to_x(0.2, 0.6);
  fixture.prepare_to_cut(slab);
  fixture.cut_multiple_edge_crossings();

  EXPECT_TRUE(fixture.node(0)->get_sorted_node_domains().empty());
  for (const double position : {0.2, 0.6})
  {
    const stk::math::Vector3d slabFaceOnEdge01 =
        (1. - position) * fixture.node(0)->coordinates() + position * fixture.node(1)->coordinates();
    EXPECT_EQ(1u,
        fixture.num_child_nodes_with_domains_on_mesh_edge_near(
            fixture.node(0), fixture.node(1), slabFaceOnEdge01, 0.05));
  }
  EXPECT_EQ(0u, fixture.global_num_multiple_edge_crossing_points());
}

TEST_F(RightTetSharpFeatureEdges, snappingSlabBesideNode_bothFacesCut)
{
  if (is_valid_proc_size_for_test()) expect_slab_beside_node_cut_on_edge01(*this);
}

// With the sharp features captured only during snapping, the slab's multiple crossings are cut with
// the end tolerance of a surface without sharp features (0.33): its face 0.2 from node 0 is at the
// node, which gets the domains, and its far face is cut.
TEST_F(RightTetSharpFeatureEdges, slabBesideNodeCapturingSharpFeaturesOnlyDuringSnapping_nearFaceIsTheNode)
{
  if (is_valid_proc_size_for_test())
  {
    cdfem_support().set_capture_sharp_features_only_during_snapping(true);
    const Cuboid slab = slab_normal_to_x(0.2, 0.6);
    prepare_to_cut(slab);
    cut_sharp_features();

    EXPECT_EQ(surfaceDomains, node(0)->get_sorted_node_domains());
    const stk::math::Vector3d nearFaceOnEdge01 = 0.8 * node(0)->coordinates() + 0.2 * node(1)->coordinates();
    const stk::math::Vector3d farFaceOnEdge01 = 0.4 * node(0)->coordinates() + 0.6 * node(1)->coordinates();
    EXPECT_EQ(0u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), nearFaceOnEdge01, 0.05));
    EXPECT_EQ(1u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), farFaceOnEdge01, 0.05));
  }
}

// With the sharp features captured only during snapping, a corner of a cuboid inside the tet is not
// cut: the tet gets no node.
TEST_F(RightTetSharpFeatureEdges, cornerInsideTetCapturingSharpFeaturesOnlyDuringSnapping_notCut)
{
  if (is_valid_proc_size_for_test())
  {
    cdfem_support().set_capture_sharp_features_only_during_snapping(true);
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    cut_sharp_features();
    EXPECT_TRUE(sharp_feature_child_nodes().empty());
  }
}

TEST_F(RightTetSharpFeatureEdges, cuttingOnlySlabBesideNode_bothFacesCut)
{
  if (is_valid_proc_size_for_test())
  {
    cdfem_support().set_cdfem_edge_degeneracy_handling(SNAP_TO_NODE);
    expect_slab_beside_node_cut_on_edge01(*this);
  }
}

// Node 0 is on the surface, and the slab's near face crosses the edge 0-1 beside it, beyond the edge
// tolerance: that crossing is cut, as is the far face.  The edge may run close to the surface between
// node 0 and the crossing, as where a node snapped onto one face of a cuboid has an edge running
// almost along that face to the adjacent face; not cutting the crossing would leave the elements with
// the edge in the wrong phase up to it.
TEST_F(RightTetSharpFeatureEdges, crossingBesideEndOnTheSurface_bothFacesCut)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid slab = slab_normal_to_x(3. * cdfem_support().get_snapper().get_edge_tolerance(), 0.6);
    prepare_to_cut(slab);
    node(0)->set_node_domains(surfaceDomains);
    cut_multiple_edge_crossings();

    EXPECT_EQ(1u, num_child_nodes_with_domains_near(node(0), 0.1));
    const stk::math::Vector3d farFaceOnEdge01 = 0.4 * node(0)->coordinates() + 0.6 * node(1)->coordinates();
    EXPECT_EQ(1u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), farFaceOnEdge01, 0.05));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// A slab 0.2 thick, tilted (rotated 0.25 about z) across the shared face, crosses edges of it
// twice each.  The tets, on one processor or two, cut the crossings on it in the same order, so
// they split it alike: cut tet by tet, each in its own order, they split it differently.
TEST_F(TwoRightTetsSharpFeatureEdges, multipleCrossingsOnEdgesOfSharedFace_bothTetsSplitAlike)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid slab(stk::math::Vector3d(0.1, 0.3, 0.2),
        stk::math::Vector3d(10., 0.2, 10.),
        stk::math::Vector3d(0., 0., 0.25));
    prepare_to_cut(slab);
    cut_multiple_edge_crossings();

    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// The slab 0.4 <= z <= 0.6 crosses the shared edge 0-5 of the four tets, and the edges from node 5
// on their shared faces, twice each.
TEST_F(FourRightTetsSharpFeatureEdges, multipleCrossingsOnSharedEdgeAndFaces_allTetsSplitAlike)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid slab(stk::math::Vector3d(0., 0., 0.5), stk::math::Vector3d(10., 10., 0.2));
    prepare_to_cut(slab);
    cut_multiple_edge_crossings();

    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 5});
    expect_face_split_alike(1, 2, {0, 3, 5});
    expect_face_split_alike(2, 3, {0, 4, 5});
    expect_face_split_alike(3, 0, {0, 1, 5});
  }
}

// Points of one feature found more than once, by different tets, on the same entity: the edge 0-4,
// shared by both tets.  Within the edge tolerance of each other they are one point, the first in
// canonical order: one node, the same in both tets, wherever they are.
TEST_F(TwoRightTetsSharpFeatureRepeats, pointsOnSharedEdgeWithinTolerance_oneNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    const double edgeTol = cdfem_support().get_snapper().get_edge_tolerance();
    const double offset = 0.5 * edgeTol;
    std::vector<IntersectionPoint> foundPoints;
    if (is_owned_element(0)) foundPoints.push_back(point({0, 4}, {0.6, 0.4}, {0, 0}));
    if (is_owned_element(1))
      foundPoints.push_back(point({0, 4}, {0.6 + offset, 0.4 - offset}, {0, 0}));
    cdmesh->share_and_insert_intersection_points(foundPoints, edgeTol);

    const std::vector<const SubElementNode *> sharpNodes = sharp_feature_child_nodes();
    ASSERT_EQ(1u, sharpNodes.size());
    const stk::math::Vector3d expectedLocation =
        0.6 * node(0)->coordinates() + 0.4 * node(4)->coordinates();
    EXPECT_NEAR(0., (sharpNodes[0]->coordinates() - expectedLocation).length(), 1.e-12);
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// Beyond the edge tolerance of each other, they are two points.
TEST_F(TwoRightTetsSharpFeatureRepeats, pointsOnSharedEdgeBeyondTolerance_twoNodes)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    const double edgeTol = cdfem_support().get_snapper().get_edge_tolerance();
    const double offset = 2. * edgeTol;
    std::vector<IntersectionPoint> foundPoints;
    if (is_owned_element(0)) foundPoints.push_back(point({0, 4}, {0.6, 0.4}, {0, 0}));
    if (is_owned_element(1))
      foundPoints.push_back(point({0, 4}, {0.6 + offset, 0.4 - offset}, {0, 0}));
    cdmesh->share_and_insert_intersection_points(foundPoints, edgeTol);

    EXPECT_EQ(2u, sharp_feature_child_nodes().size());
    expect_no_hanging_nodes();
    expect_face_split_alike(0, 1, {0, 2, 4});
  }
}

// A corner captures a crossing of one of its sharp edges within the edge tolerance of it on the
// same entity (its domains contain the crossing's), wherever each was found: only the corner's
// node.
TEST_F(TwoRightTetsSharpFeatureRepeats, crossingWithinToleranceOfCornerOnSharedEdge_cornerNodeOnly)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(cuboidWithCornerInsideRightTet);
    const double edgeTol = cdfem_support().get_snapper().get_edge_tolerance();
    const double offset = 0.5 * edgeTol;
    std::vector<IntersectionPoint> foundPoints;
    if (is_owned_element(1))
      foundPoints.push_back(point({0, 4}, {0.6 + offset, 0.4 - offset}, {0, 0}));
    if (is_owned_element(0)) foundPoints.push_back(point({0, 4}, {0.6, 0.4}, cornerDomains));
    cdmesh->share_and_insert_intersection_points(foundPoints, edgeTol);

    const std::vector<const SubElementNode *> sharpNodes = sharp_feature_child_nodes();
    ASSERT_EQ(1u, sharpNodes.size());
    EXPECT_EQ(cornerDomains, sharpNodes[0]->get_sorted_node_domains());
  }
}

// A box between the planes x = xLo and x = xHi, reaching far beyond the tets in y and z, as a
// connected faceted surface without sharp features (no sharp feature angle).  The facets are all
// owned by processor 0.
static std::unique_ptr<MeshFacetedSurface> faceted_slab_normal_to_x_without_sharp_features(
    const double xLo, const double xHi)
{
  constexpr double w = 5.;
  const std::vector<stk::math::Vector3d> vertices{{xLo, -w, -w},
      {xHi, -w, -w},
      {xHi, w, -w},
      {xLo, w, -w},
      {xLo, -w, w},
      {xHi, -w, w},
      {xHi, w, w},
      {xLo, w, w}};
  std::vector<std::array<unsigned, 3>> facets;
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
    facets = {{{0, 4, 7}},
        {{0, 7, 3}},
        {{1, 2, 6}},
        {{1, 6, 5}},
        {{0, 1, 5}},
        {{0, 5, 4}},
        {{2, 3, 7}},
        {{2, 7, 6}},
        {{0, 3, 2}},
        {{0, 2, 1}},
        {{4, 5, 6}},
        {{4, 6, 7}}};
  return std::make_unique<MeshFacetedSurface>(MPI_COMM_WORLD, 3, vertices, facets);
}

// The multiple crossings of a faceted surface without sharp features are cut too: both faces of a
// slab across the edge 0-1, away from its ends, get nodes on it.
TEST_F(RightTetSharpFeatureEdges, facetedSlabWithoutSharpFeaturesAwayFromNodes_bothFacesCut)
{
  if (is_valid_proc_size_for_test())
  {
    const auto slab = faceted_slab_normal_to_x_without_sharp_features(0.4, 0.6);
    ASSERT_FALSE(slab->has_sharp_features());
    prepare_to_cut(*slab);
    cut_multiple_edge_crossings();

    const stk::math::Vector3d slabOnEdge01 = 0.5 * node(0)->coordinates() + 0.5 * node(1)->coordinates();
    EXPECT_EQ(2u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), slabOnEdge01, 0.15));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// For a surface without sharp features, a crossing within 0.33 of the edge from an end is that end:
// the slab's near face, 0.2 from node 0, gives node 0 the surface's domains, and its far face gets a
// node on the edge 0-1.
TEST_F(RightTetSharpFeatureEdges, facetedSlabWithoutSharpFeaturesBesideNode_nearFaceIsTheNode)
{
  if (is_valid_proc_size_for_test())
  {
    const auto slab = faceted_slab_normal_to_x_without_sharp_features(0.2, 0.6);
    prepare_to_cut(*slab);
    cut_multiple_edge_crossings();

    EXPECT_EQ(surfaceDomains, node(0)->get_sorted_node_domains());
    EXPECT_EQ(0u, num_child_nodes_with_domains_near(node(0), 0.3));
    const stk::math::Vector3d farFaceOnEdge01 = 0.4 * node(0)->coordinates() + 0.6 * node(1)->coordinates();
    EXPECT_EQ(1u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), farFaceOnEdge01, 0.05));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// What decomposition_has_changed checks against: the edge crossings of a faceted surface without
// sharp features include those that the signs at an edge's ends do not show, the slab's two faces
// on the edge 0-1.
TEST_F(RightTetSharpFeatureEdges, facetedSlabWithoutSharpFeatures_allEdgeIntersectionPointsIncludeMultipleCrossings)
{
  if (is_valid_proc_size_for_test())
  {
    const auto slab = faceted_slab_normal_to_x_without_sharp_features(0.4, 0.6);
    const AnalyticSurfaceInterfaceGeometry geometry({levelset_fields()[0].identifier},
        {slab.get()},
        get_aux_meta().active_part(),
        cdfem_support(),
        Phase_Support::get(mMesh.mesh_meta_data()));
    const std::vector<IntersectionPoint> points = geometry.get_all_edge_intersection_points(mMesh);

    const std::vector<stk::mesh::Entity> edge01{get_assigned_node_for_index(0), get_assigned_node_for_index(1)};
    size_t numOnEdge01 = 0;
    for (auto & point : points)
    {
      std::vector<stk::mesh::Entity> pointNodes = point.get_nodes();
      std::vector<stk::mesh::Entity> edgeNodes = edge01;
      std::sort(pointNodes.begin(), pointNodes.end());
      std::sort(edgeNodes.begin(), edgeNodes.end());
      if (pointNodes == edgeNodes) ++numOnEdge01;
    }
    EXPECT_EQ(2u, numOnEdge01);
  }
}

// The faces of a slab thinner than the edge tolerance cross the edge 0-1 within it of each other:
// they are close crossings, where the edge crosses the slab and back, so neither is cut, and no node
// is made for them.
TEST_F(RightTetSharpFeatureEdges, slabThinnerThanTolerance_ignored)
{
  if (is_valid_proc_size_for_test())
  {
    const double thickness = 0.5 * cdfem_support().get_snapper().get_edge_tolerance();
    const Cuboid slab = slab_normal_to_x(0.4, 0.4 + thickness);
    prepare_to_cut(slab);
    cut_multiple_edge_crossings();

    const stk::math::Vector3d slabOnEdge01 =
        0.6 * node(0)->coordinates() + 0.4 * node(1)->coordinates();
    EXPECT_EQ(
        0u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), slabOnEdge01, 0.1));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// Twice as thick as the tolerance, both faces of the slab get nodes on the edge 0-1.
TEST_F(RightTetSharpFeatureEdges, slabThickerThanTolerance_bothFacesCut)
{
  if (is_valid_proc_size_for_test())
  {
    const double thickness = 2. * cdfem_support().get_snapper().get_edge_tolerance();
    const Cuboid slab = slab_normal_to_x(0.4, 0.4 + thickness);
    prepare_to_cut(slab);
    cut_multiple_edge_crossings();

    const stk::math::Vector3d slabOnEdge01 =
        0.6 * node(0)->coordinates() + 0.4 * node(1)->coordinates();
    EXPECT_EQ(2u,
        num_child_nodes_with_domains_on_mesh_edge_near(
            node(0), node(1), slabOnEdge01, 2. * thickness));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// With both ends on the surface (their domains capture it), the edge 0-1 is a chord of it: the
// faces of the slab crossing it are not cut, as cutting a chord of a curved surface where the
// surface crosses it only makes chords that cross it again, without end.
TEST_F(RightTetSharpFeatureEdges, slabCrossingChordOfItsSurface_notCut)
{
  if (is_valid_proc_size_for_test())
  {
    const double thickness = 2. * cdfem_support().get_snapper().get_edge_tolerance();
    const Cuboid slab = slab_normal_to_x(0.4, 0.4 + thickness);
    prepare_to_cut(slab);
    node(0)->set_node_domains(surfaceDomains);
    node(1)->set_node_domains(surfaceDomains);
    cut_multiple_edge_crossings();

    const stk::math::Vector3d slabOnEdge01 =
        0.6 * node(0)->coordinates() + 0.4 * node(1)->coordinates();
    EXPECT_EQ(0u,
        num_child_nodes_with_domains_on_mesh_edge_near(
            node(0), node(1), slabOnEdge01, 2. * thickness));
  }
}

TEST_F(RightTriSharpFeatureCutting, insertInteriorPoint_threeLeavesJoinedAtNewNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(rectangleWithCornerInsideRightTri);
    element(0).insert_intersection_point(
        *cdmesh, point({0, 1, 2}, {0.4, 0.3, 0.3}, cornerDomains2D));

    const std::vector<const SubElement *> elemLeaves = leaves(0);
    ASSERT_EQ(3u, elemLeaves.size());
    const SubElementNode * newNode = SubElementNode::common_child({node(0), node(1), node(2)});
    ASSERT_TRUE(newNode);
    EXPECT_EQ(cornerDomains2D, newNode->get_sorted_node_domains());
    for (auto * leaf : elemLeaves)
    {
      EXPECT_EQ(1, std::count(leaf->get_nodes().begin(), leaf->get_nodes().end(), newNode));
    }
  }
}

TEST_F(RightTriSharpFeatureCutting, insertEdgePoint_twoLeavesAndEdgeNode)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(rectangleWithCornerInsideRightTri);
    element(0).insert_intersection_point(*cdmesh, point({0, 1}, {0.6, 0.4}, cornerDomains2D));

    EXPECT_EQ(2u, leaves(0).size());
    const SubElementNode * edgeNode = SubElementNode::common_child({node(0), node(1)});
    ASSERT_TRUE(edgeNode);
    EXPECT_EQ(cornerDomains2D, edgeNode->get_sorted_node_domains());
  }
}

TEST_F(RightTriSharpFeatureCutting, insertPointAtNode_noCutAndDomainsAdded)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(rectangleWithCornerInsideRightTri);
    element(0).insert_intersection_point(*cdmesh, point({0}, {1.}, cornerDomains2D));

    EXPECT_EQ(1u, leaves(0).size());
    EXPECT_EQ(cornerDomains2D, node(0)->get_sorted_node_domains());
  }
}

// A point on the edge 0-2 shared by the two tris, found by the first: both are bisected at it,
// wherever each is.
TEST_F(TwoRightTrisSharpFeatureCutting, pointOnSharedEdge_bothTrisBisected)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(rectangleWithCornerInsideRightTri);
    std::vector<IntersectionPoint> foundPoints;
    if (is_owned_element(0)) foundPoints.push_back(point({0, 2}, {0.5, 0.5}, cornerDomains2D));
    cdmesh->share_and_insert_intersection_points(foundPoints, 1.e-12);

    for (size_t elemIndex = 0; elemIndex < 2; ++elemIndex)
    {
      if (is_owned_element(elemIndex))
      {
        EXPECT_EQ(2u, leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
  }
}

// A rectangle whose corner is at (x, y), extending in +x and +y.
static Cuboid rectangle_with_corner_at(const double x, const double y)
{
  return Cuboid(stk::math::Vector3d(x + 0.5, y + 0.5, 0.), stk::math::Vector3d(1., 1., 10.));
}

TEST_F(RightTriSharpFeatureCorners, cornerInsideTri_nodeInsideAndThreeLeaves)
{
  if (is_valid_proc_size_for_test())
  {
    prepare_to_cut(rectangleWithCornerInsideRightTri);
    cut_corners();

    const std::vector<const SubElementNode *> sharpNodes = sharp_feature_child_nodes();
    ASSERT_EQ(1u, sharpNodes.size());
    EXPECT_EQ(sharpNodes[0], SubElementNode::common_child({node(0), node(1), node(2)}));
    EXPECT_EQ(cornerDomains2D, sharpNodes[0]->get_sorted_node_domains());
    EXPECT_EQ(3u, leaves(0).size());
  }
}

// Within the tolerance of node 0, the corner is at node 0, which gets its domains.
TEST_F(RightTriSharpFeatureCorners, cornerWithinToleranceOfNode_nodeGetsTheDomains)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid rectangle = rectangle_with_corner_at(1.e-4, 1.e-4);
    prepare_to_cut(rectangle);
    cut_corners();

    EXPECT_TRUE(sharp_feature_child_nodes().empty());
    EXPECT_EQ(cornerDomains2D, node(0)->get_sorted_node_domains());
  }
}

// Within the tolerance of the edge 0-2 that the tris share, the corner is a node on it, which both
// tris are split at, wherever each is.
TEST_F(TwoRightTrisSharpFeatureCutting, cornerWithinToleranceOfSharedEdge_bothTrisSplitAtEdgeNode)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid rectangle = rectangle_with_corner_at(1.e-4, 0.4);
    prepare_to_cut(rectangle);
    cut_corners();

    const SubElementNode * edgeNode = SubElementNode::common_child({node(0), node(2)});
    ASSERT_TRUE(edgeNode);
    EXPECT_EQ(cornerDomains2D, edgeNode->get_sorted_node_domains());
    for (size_t elemIndex = 0; elemIndex < 2; ++elemIndex)
    {
      if (is_owned_element(elemIndex))
      {
        EXPECT_EQ(2u, leaves(elemIndex).size()) << "element " << elemIndex;
      }
    }
    expect_no_hanging_nodes();
  }
}

// Beyond the tolerance, the corner is a node inside the first tri; the second, without the
// interface, is not split.
TEST_F(TwoRightTrisSharpFeatureCutting, cornerBeyondToleranceOfSharedEdge_nodeInsideFirstTri)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid rectangle = rectangle_with_corner_at(1.e-2, 0.4);
    prepare_to_cut(rectangle);
    cut_corners();

    EXPECT_EQ(nullptr, SubElementNode::common_child({node(0), node(2)}));
    if (is_owned_element(0))
    {
      EXPECT_EQ(3u, leaves(0).size());
    }
    if (is_owned_element(1))
    {
      EXPECT_FALSE(element(1).have_subelements());
    }
    expect_no_hanging_nodes();
  }
}

// Level sets per phase: nodes 0, 2 and 1 in phases 0, 1 and 2, with phase 2 scaled so that the
// point where the three meet has weight 1e-4 at node 1: within the tolerance of the shared edge
// 0-2, it is a node on it, at its middle, which both tris are split at.
TEST_F(
    TwoRightTris3LSPerPhaseSharpFeatureCorners, triplePointNearSharedEdge_bothTrisSplitAtEdgeNode)
{
  if (is_valid_proc_size_for_test())
  {
    const double eta = 1.e-4;
    set_level_sets_by_node_phase({0, 2, 1, 0}, {1., 1., 0.5 * (1. - eta) / eta});
    prepare_to_cut_level_sets();
    cut_corners();

    const SubElementNode * edgeNode = SubElementNode::common_child({node(0), node(2)});
    ASSERT_TRUE(edgeNode);
    EXPECT_EQ((std::vector<int>{0, 1, 2}), edgeNode->get_sorted_node_domains());
    EXPECT_NEAR(0.,
        (edgeNode->coordinates() - 0.5 * (node(0)->coordinates() + node(2)->coordinates()))
            .length(),
        1.e-12);
    expect_no_hanging_nodes();
  }
}

// In 2D: the edges of the right tri from node 1, to nodes 0 and 2, cross the band 0.4 <= x <= 0.6
// twice each, and the third does not: four points, each on its mesh edge.
TEST_F(RightTriMultipleEdgeCrossings, edgesAcrossBand_twoCrossingsOnEach)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid band = slab_normal_to_x(0.4, 0.6);
    prepare_to_cut(band);

    const std::vector<IntersectionPoint> points =
        element(0).get_multiple_edge_crossing_points(*cdmesh);
    ASSERT_EQ(4u, points.size());
    const stk::mesh::Entity node1 = get_assigned_node_for_index(1);
    for (auto & point : points)
    {
      EXPECT_EQ(surfaceDomains, point.get_sorted_domains());
      ASSERT_EQ(2u, point.get_nodes().size());
      const auto iter = std::find(point.get_nodes().begin(), point.get_nodes().end(), node1);
      ASSERT_TRUE(iter != point.get_nodes().end());
      const double weightOfNode1 =
          point.get_weights()[std::distance(point.get_nodes().begin(), iter)];
      EXPECT_TRUE(std::abs(weightOfNode1 - 0.4) < 1.e-12 || std::abs(weightOfNode1 - 0.6) < 1.e-12)
          << weightOfNode1;
    }
  }
}

// The band is 1e-4 from node 0, within the tolerance of it along its edge: node 0 gets the
// surface's domains, and no node is made beside it.
TEST_F(RightTriMultipleEdgeCrossings, bandWithinToleranceOfNode_nodeGetsTheDomains)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid band = slab_normal_to_x(1.e-4, 0.6);
    prepare_to_cut(band);
    cut_multiple_edge_crossings();

    EXPECT_EQ(surfaceDomains, node(0)->get_sorted_node_domains());
    EXPECT_EQ(0u, num_child_nodes_with_domains_near(node(0), 0.1));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// The faces of a band thinner than the edge tolerance cross the edge 0-1 within it of each other:
// they are close crossings, where the edge crosses the band and back, so neither is cut, and no node
// is made for them.
TEST_F(RightTriMultipleEdgeCrossings, bandThinnerThanTolerance_ignored)
{
  if (is_valid_proc_size_for_test())
  {
    const double thickness = 0.5 * cdfem_support().get_snapper().get_edge_tolerance();
    const Cuboid band = slab_normal_to_x(0.4, 0.4 + thickness);
    prepare_to_cut(band);
    cut_multiple_edge_crossings();

    const stk::math::Vector3d bandOnEdge01 =
        0.6 * node(0)->coordinates() + 0.4 * node(1)->coordinates();
    EXPECT_EQ(
        0u, num_child_nodes_with_domains_on_mesh_edge_near(node(0), node(1), bandOnEdge01, 0.1));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// Twice as thick as the tolerance, both faces of the band get nodes on the edge 0-1.
TEST_F(RightTriMultipleEdgeCrossings, bandThickerThanTolerance_bothFacesCut)
{
  if (is_valid_proc_size_for_test())
  {
    const double thickness = 2. * cdfem_support().get_snapper().get_edge_tolerance();
    const Cuboid band = slab_normal_to_x(0.4, 0.4 + thickness);
    prepare_to_cut(band);
    cut_multiple_edge_crossings();

    const stk::math::Vector3d bandOnEdge01 =
        0.6 * node(0)->coordinates() + 0.4 * node(1)->coordinates();
    EXPECT_EQ(2u,
        num_child_nodes_with_domains_on_mesh_edge_near(
            node(0), node(1), bandOnEdge01, 2. * thickness));
    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
  }
}

// The band crosses the edge 0-2 that the two tris share twice: both tris get the nodes, wherever
// each is.
TEST_F(TwoRightTrisSharpFeatureCutting, multipleCrossingsOnSharedEdge_bothTrisSplitAtThem)
{
  if (is_valid_proc_size_for_test())
  {
    const Cuboid band(stk::math::Vector3d(0., 0.5, 0.), stk::math::Vector3d(10., 0.2, 10.));
    prepare_to_cut(band);
    cut_multiple_edge_crossings();

    EXPECT_EQ(0u, global_num_multiple_edge_crossing_points());
    expect_no_hanging_nodes();
  }
}

// Seven crossings beside an end, with a tolerance of 0.1 of the edge: pairing them from the end at
// the origin, (0.11, 0.12), (0.24, 0.25) and (0.26, 0.35) are close crossings, and the edge is cut
// once, at 0.13.  Walking the edge from its other end gives the same cut.
TEST(MultipleCrossingPositionsToCut, closeCrossings_sameCutWalkingTheEdgeEitherWay)
{
  const stk::math::Vector3d xOrigin(0., 0., 0.);
  const stk::math::Vector3d xOther(1., 0., 0.);
  const std::vector<double> positionsFromOrigin{0.11, 0.12, 0.13, 0.24, 0.25, 0.26, 0.35};
  std::vector<double> positionsFromOther;
  for (auto it = positionsFromOrigin.rbegin(); it != positionsFromOrigin.rend(); ++it)
    positionsFromOther.push_back(1. - *it);

  const std::vector<double> cutFromOrigin =
      multiple_crossing_positions_to_cut(positionsFromOrigin, xOrigin, xOther, 0.1, 0.1);
  const std::vector<double> cutFromOther =
      multiple_crossing_positions_to_cut(positionsFromOther, xOther, xOrigin, 0.1, 0.1);

  ASSERT_EQ(1u, cutFromOrigin.size());
  EXPECT_NEAR(0.13, cutFromOrigin[0], 1.e-12);
  ASSERT_EQ(1u, cutFromOther.size());
  EXPECT_NEAR(1. - 0.13, cutFromOther[0], 1.e-12);
}

// A crossing within the end tolerance of an end is that end, from either end.
TEST(MultipleCrossingPositionsToCut, crossingWithinEndTolerance_isTheEnd)
{
  const stk::math::Vector3d x0(0., 0., 0.);
  const stk::math::Vector3d x1(2., 0., 0.);
  EXPECT_EQ(std::vector<double>({0., 0.5}), multiple_crossing_positions_to_cut({0.04, 0.5}, x0, x1, 0.01, 0.1));
  EXPECT_EQ(std::vector<double>({0.5, 1.}), multiple_crossing_positions_to_cut({0.5, 0.96}, x1, x0, 0.01, 0.1));
}

// An edge whose ends both have the interface's domains is a chord of it; one end is not.
TEST(MultipleCrossingPositionsToInsert, edgeWithBothEndsOnInterface_isChord)
{
  const std::vector<int> interfaceDomains{0, 1};
  const std::vector<int> noDomains;
  EXPECT_TRUE(edge_is_chord_of_interface({0, 1}, {0, 1, 2}, interfaceDomains));
  EXPECT_FALSE(edge_is_chord_of_interface({0, 1}, noDomains, interfaceDomains));
  EXPECT_FALSE(edge_is_chord_of_interface({0, 2}, {0, 1}, interfaceDomains));
}

// A crossing within the end tolerance of an end is not inserted where that end's node already has
// the interface's domains: it is that node.
TEST(MultipleCrossingPositionsToInsert, crossingAtEndWithTheDomains_isTheNode)
{
  const stk::math::Vector3d x0(0., 0., 0.);
  const stk::math::Vector3d x1(2., 0., 0.);
  const std::vector<int> interfaceDomains{0, 1};
  const std::vector<int> noDomains;
  EXPECT_EQ(std::vector<double>({0., 0.5}),
      multiple_crossing_positions_to_insert({0.04, 0.5}, x0, x1, noDomains, noDomains, interfaceDomains, 0.01, 0.1));
  EXPECT_EQ(std::vector<double>({0.5}),
      multiple_crossing_positions_to_insert({0.04, 0.5}, x0, x1, interfaceDomains, noDomains, interfaceDomains, 0.01, 0.1));
  EXPECT_EQ(std::vector<double>({0., 0.5}),
      multiple_crossing_positions_to_insert({0.04, 0.5}, x0, x1, noDomains, interfaceDomains, interfaceDomains, 0.01, 0.1));
}

} // namespace krino

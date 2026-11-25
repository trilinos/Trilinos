// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AuxMetaData.hpp>
#include <Akri_KinematicMeshUtils.hpp>
#include <Akri_KinematicObjectives.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Objective_KinematicMesh.hpp>
#include <Akri_Objective_MeanRatioAboutNode.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_SecondOrderMeshSpecs.hpp>
#include <Akri_SmoothMesh.hpp>
#include <Akri_SmoothNode.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_StkMeshFixture.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

namespace krino
{

namespace
{

stk::math::Vector3d node_location(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  return get_vector_field(mesh, coordsField, node, mesh.mesh_meta_data().spatial_dimension());
}

void move_node_to(
    stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::math::Vector3d & loc)
{
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  auto coordsData = coordsField.data<double, stk::mesh::ReadWrite>();
  auto nodeCoords = coordsData.entity_values(node);
  for (unsigned j = 0; j < mesh.mesh_meta_data().spatial_dimension(); ++j)
    nodeCoords(stk::mesh::ComponentIdx(j)) = loc[j];
}

CoordinatesFieldRef coordinates_field(const stk::mesh::BulkData & mesh)
{
  return CoordinatesFieldRef(
      mesh.mesh_meta_data().coordinate_field(), mesh.mesh_meta_data().spatial_dimension());
}

// No gradient normal to the sides of the [-1,1]^dim box, so its boundary nodes stay on it.
GradientFilter boundary_of_box_filter()
{
  return [](const stk::mesh::Entity, const stk::math::Vector3d & x, stk::math::Vector3d & dir)
  {
    for (unsigned i = 0; i < 3; ++i)
      if (x[i] == -1. || x[i] == 1.) dir[i] = 0.;
  };
}

void krino_lbfgs(const ObjectiveInterface & objFn, DistributedVector & x)
{
  const double xTol = 1.e-4;
  lbfgs(objFn, x, xTol, 0, objFn.compute_value(x));
}

// A mesh spec built twice, with linear (tet4, tri3) and quadratic (tet10, tri6) elements: the
// vertex nodes of the quadratic mesh are the nodes of the linear one, with the same ids, followed
// by a mid-edge node at the middle of each edge.
template <typename LINEAR_SPEC, stk::topology::topology_t QUADRATIC> struct LinearAndQuadraticMeshes
{
  static constexpr stk::topology::topology_t LINEAR = LINEAR_SPEC::TOPOLOGY;

  LinearAndQuadraticMeshes(
      const typename StkMeshAndBuilder<LINEAR>::ElementOwnersByProcCount & elementOwnersByProcCount = {})
  {
    for (stk::mesh::BulkData * mesh : {&linear.get_mesh(), &quadratic.get_mesh()})
    {
      stk::mesh::MetaData & meta = mesh->mesh_meta_data();
      const unsigned dim = meta.spatial_dimension();
      auto & field = meta.declare_field<double>(stk::topology::NODE_RANK, "nodal_metric");
      stk::mesh::put_field_on_mesh(field, meta.universal_part(), dim * dim, nullptr);
    }
    linear.build_mesh_from_spec(linearSpec, elementOwnersByProcCount);
    quadratic.build_mesh_from_spec(quadraticSpec, elementOwnersByProcCount);
  }

  // A diagonal metric that varies with position, at every node (vertex or mid-edge) of the mesh.
  static FieldRef set_varying_nodal_metric(stk::mesh::BulkData & mesh)
  {
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    const FieldRef field(mesh.mesh_meta_data().get_field(stk::topology::NODE_RANK, "nodal_metric"));
    auto data = field.data<double, stk::mesh::ReadWrite>();
    for (auto * bucket : mesh.buckets(stk::topology::NODE_RANK))
      for (auto node : *bucket)
      {
        const stk::math::Vector3d x = node_location(mesh, node);
        auto metric = data.entity_values(node);
        for (unsigned i = 0; i < dim; ++i)
          for (unsigned j = 0; j < dim; ++j)
            metric(stk::mesh::ComponentIdx(i * dim + j)) =
                (i != j) ? 0. : 1. / ((1.5 + 0.4 * x[i]) * (1.5 + 0.4 * x[i]));
      }
    return field;
  }

  size_t num_vertex_nodes() const { return linearSpec.nodeLocs.size(); }

  stk::math::Vector3d spec_location(const size_t i) const
  {
    return stk::math::Vector3d(
        linearSpec.nodeLocs[i].data(), linear.get_mesh().mesh_meta_data().spatial_dimension());
  }

  // Moves vertex node i of both meshes by disp, as computed by both from the spec.
  void displace_vertex_node(const size_t i, const stk::math::Vector3d & disp)
  {
    const stk::math::Vector3d loc = spec_location(i);
    for (stk::mesh::BulkData * mesh : {&linear.get_mesh(), &quadratic.get_mesh()})
    {
      const stk::mesh::Entity node = (mesh == &linear.get_mesh()) ?
          linear.get_assigned_node_for_index(i) :
          quadratic.get_assigned_node_for_index(i);
      if (mesh->is_valid(node)) move_node_to(*mesh, node, loc + disp);
    }
  }

  // Moves every mid-edge node of the quadratic mesh off its edge midpoint, by an offset that varies
  // from node to node: their locations are to be ignored by smoothing, and then reset on the edges
  // with a moved vertex node.
  void displace_mid_edge_nodes_off_edge_midpoints()
  {
    const unsigned dim = quadratic.get_mesh().mesh_meta_data().spatial_dimension();
    for (size_t i = num_vertex_nodes(); i < quadraticSpec.nodeLocs.size(); ++i)
    {
      const stk::mesh::Entity node = quadratic.get_assigned_node_for_index(i);
      if (quadratic.get_mesh().is_valid(node))
      {
        const double t = 0.05 * (1 + i % 3);
        const stk::math::Vector3d offset =
            (dim == 3) ? stk::math::Vector3d(t, -t, 0.5 * t) : stk::math::Vector3d(t, -t, 0.);
        move_node_to(quadratic.get_mesh(),
            node,
            stk::math::Vector3d(quadraticSpec.nodeLocs[i].data(), dim) + offset);
      }
    }
  }

  void sync_coordinates()
  {
    for (stk::mesh::BulkData * mesh : {&linear.get_mesh(), &quadratic.get_mesh()})
      stk::mesh::communicate_field_data(*mesh, {mesh->mesh_meta_data().coordinate_field()});
  }

  // Collective.  Gives one mid-edge node owned by processor 0, that no other processor shares, of an
  // edge with a vertex node of index at least firstMovedVertexIndex, to processor 1, so that its owner
  // owns no element with its edge, as can happen after the mesh is modified in parallel.  Returns the
  // number of mid-edge nodes given, 1 if there was one to give.
  int give_unshared_mid_edge_node_to_proc_1(const size_t firstMovedVertexIndex)
  {
    const auto is_moved_vertex = [&](const stk::mesh::Entity node)
    {
      for (size_t i = firstMovedVertexIndex; i < num_vertex_nodes(); ++i)
        if (quadratic.get_assigned_node_for_index(i) == node) return true;
      return false;
    };
    stk::mesh::BulkData & mesh = quadratic.get_mesh();
    stk::mesh::EntityProcVec midNodeToGive;
    if (mesh.parallel_rank() == 0)
    {
      std::array<unsigned, 3> edgeNodeOrdinals;
      for (auto elem : quadratic.get_owned_elements())
      {
        const stk::topology topology = mesh.bucket(elem).topology();
        const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
        for (unsigned e = 0; e < topology.num_edges() && midNodeToGive.empty(); ++e)
        {
          topology.edge_node_ordinals(e, edgeNodeOrdinals.data());
          const stk::mesh::Entity midNode = elemNodes[edgeNodeOrdinals[2]];
          const bool isEdgeOfMovedVertex = is_moved_vertex(elemNodes[edgeNodeOrdinals[0]]) ||
              is_moved_vertex(elemNodes[edgeNodeOrdinals[1]]);
          if (isEdgeOfMovedVertex && mesh.bucket(midNode).owned() && !mesh.bucket(midNode).shared())
            midNodeToGive.emplace_back(midNode, 1);
        }
        if (!midNodeToGive.empty()) break;
      }
    }
    mesh.change_entity_owner(midNodeToGive);
    int numGiven = midNodeToGive.size();
    int globalNumGiven = 0;
    stk::all_reduce_sum(mesh.parallel(), &numGiven, &globalNumGiven, 1);
    return globalNumGiven;
  }

  // The locations of the nodes of the quadratic mesh on this proc, once displaced and before
  // smoothing.
  void record_quadratic_node_locations_before_smoothing()
  {
    const stk::mesh::BulkData & mesh = quadratic.get_mesh();
    for (auto * bucket : mesh.buckets(stk::topology::NODE_RANK))
      for (auto node : *bucket)
        quadraticNodeLocationsBeforeSmoothing[node] = node_location(mesh, node);
  }

  bool is_quadratic_node_moved_by_smoothing(const stk::mesh::Entity node) const
  {
    return node_location(quadratic.get_mesh(), node) !=
        quadraticNodeLocationsBeforeSmoothing.at(node);
  }

  // Every vertex node of the quadratic mesh is exactly where the same node of the linear mesh is.
  // Every mid-edge node of an edge with a moved vertex node is exactly at the middle of its edge,
  // and every mid-edge node of an edge whose vertex nodes did not move is where it was, so that
  // curved edges that smoothing leaves alone stay curved.
  void expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints()
  {
    for (size_t i = 0; i < num_vertex_nodes(); ++i)
    {
      const stk::mesh::Entity linearNode = linear.get_assigned_node_for_index(i);
      const stk::mesh::Entity quadraticNode = quadratic.get_assigned_node_for_index(i);
      if (linear.get_mesh().is_valid(linearNode))
      {
        EXPECT_EQ(node_location(linear.get_mesh(), linearNode),
            node_location(quadratic.get_mesh(), quadraticNode))
            << "vertex node " << i;
      }
    }

    const stk::mesh::BulkData & mesh = quadratic.get_mesh();
    for_each_owned_edge(
        [&](const stk::mesh::Entity elem,
            const unsigned e,
            const stk::mesh::Entity node0,
            const stk::mesh::Entity node1,
            const stk::mesh::Entity midNode)
        {
          if (is_quadratic_edge_moved_by_smoothing(node0, node1))
          {
            const stk::math::Vector3d midpoint =
                0.5 * (node_location(mesh, node0) + node_location(mesh, node1));
            EXPECT_EQ(midpoint, node_location(mesh, midNode))
                << "element " << mesh.identifier(elem) << " edge " << e;
          }
          else
          {
            EXPECT_FALSE(is_quadratic_node_moved_by_smoothing(midNode))
                << "element " << mesh.identifier(elem) << " edge " << e;
          }
        });
  }

  // Collective.  Some edge of the quadratic mesh has vertex nodes that smoothing did not move, so
  // that the mid-edge nodes left where they were are checked.
  void expect_some_edges_with_vertex_nodes_not_moved_by_smoothing()
  {
    int numUnmovedEdges = 0;
    for_each_owned_edge(
        [&](const stk::mesh::Entity,
            const unsigned,
            const stk::mesh::Entity node0,
            const stk::mesh::Entity node1,
            const stk::mesh::Entity)
        {
          if (!is_quadratic_edge_moved_by_smoothing(node0, node1)) ++numUnmovedEdges;
        });
    int globalNumUnmovedEdges = 0;
    stk::all_reduce_sum(quadratic.get_mesh().parallel(), &numUnmovedEdges, &globalNumUnmovedEdges, 1);
    EXPECT_GT(globalNumUnmovedEdges, 0);
  }

  template <typename EDGE_FUNCTION> void for_each_owned_edge(const EDGE_FUNCTION & edgeFunction) const
  {
    const stk::mesh::BulkData & mesh = quadratic.get_mesh();
    std::array<unsigned, 3> edgeNodeOrdinals;
    for (auto elem : quadratic.get_owned_elements())
    {
      const stk::topology topology = mesh.bucket(elem).topology();
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      for (unsigned e = 0; e < topology.num_edges(); ++e)
      {
        topology.edge_node_ordinals(e, edgeNodeOrdinals.data());
        edgeFunction(elem,
            e,
            elemNodes[edgeNodeOrdinals[0]],
            elemNodes[edgeNodeOrdinals[1]],
            elemNodes[edgeNodeOrdinals[2]]);
      }
    }
  }

  bool is_quadratic_edge_moved_by_smoothing(
      const stk::mesh::Entity node0, const stk::mesh::Entity node1) const
  {
    return is_quadratic_node_moved_by_smoothing(node0) || is_quadratic_node_moved_by_smoothing(node1);
  }

  // Collective.  Whether vertex node i of the linear mesh is within tol of where the spec has it,
  // as when smoothing restores the regular mesh after the node was moved.
  bool is_linear_vertex_node_restored(const size_t i, const double tol)
  {
    const stk::mesh::Entity node = linear.get_assigned_node_for_index(i);
    int notRestored = 0;
    if (linear.get_mesh().is_valid(node) && linear.get_mesh().bucket(node).owned())
    {
      notRestored = (node_location(linear.get_mesh(), node) - spec_location(i)).length() > tol;
    }
    int anyNotRestored = 0;
    stk::all_reduce_max(linear.get_mesh().parallel(), &notRestored, &anyNotRestored, 1);
    return anyNotRestored == 0;
  }

  LINEAR_SPEC linearSpec;
  HigherOrderSimplexMeshSpec<LINEAR_SPEC, QUADRATIC> quadraticSpec;
  StkMeshAndBuilder<LINEAR> linear;
  StkMeshAndBuilder<QUADRATIC> quadratic;
  std::unordered_map<stk::mesh::Entity, stk::math::Vector3d> quadraticNodeLocationsBeforeSmoothing;
};

using CubeOf24Tet4AndTet10 = LinearAndQuadraticMeshes<CubeOf24Tets, stk::topology::TETRAHEDRON_10>;
using SquareOf8Tri3AndTri6 = LinearAndQuadraticMeshes<SquareOf8Tris, stk::topology::TRIANGLE_6_2D>;

// The cube's face and center nodes (8 to 14) moved within their faces, as in the CubeOf24Tets
// smoothing tests.
void displace_cube_face_and_center_nodes(CubeOf24Tet4AndTet10 & meshes)
{
  const std::array<stk::math::Vector3d, 7> disp{{{0.0, -0.3, -0.1},
      {0.0, 0.8, 0.4},
      {0.1, 0.0, 0.8},
      {0.1, 0.0, 0.0},
      {0.2, 0.7, 0.0},
      {-0.8, 0.8, 0.0},
      {0.4, -0.3, 0.8}}};
  for (size_t i = 8; i < 15; ++i)
  {
    stk::math::Vector3d nodeDisp = disp[i - 8];
    const stk::math::Vector3d specLoc = meshes.spec_location(i);
    for (unsigned j = 0; j < 3; ++j)
      if (specLoc[j] == -1. || specLoc[j] == 1.) nodeDisp[j] = 0.;
    meshes.displace_vertex_node(i, nodeDisp);
  }
  meshes.displace_mid_edge_nodes_off_edge_midpoints();
  meshes.sync_coordinates();
  meshes.record_quadratic_node_locations_before_smoothing();
}

template <typename MESHES> void smooth_both_with_kinematic_objective(MESHES & meshes)
{
  const std::shared_ptr<KinematicElementObjective> elemObjPtr =
      make_modified_neo_hookean_objective();
  const KinematicElementObjective & elemObj = *elemObjPtr;
  for (auto * mesh : {&meshes.linear.get_mesh(), &meshes.quadratic.get_mesh()})
  {
    const stk::mesh::Selector elemSelector = mesh->mesh_meta_data().universal_part();
    KinematicConstantIsotropicMeshObjective smoothingObj(
        *mesh, coordinates_field(*mesh), elemSelector, elemObj, 1.0);
    smoothingObj.set_gradient_filter(boundary_of_box_filter());
    smooth_mesh_by_optimization(*mesh, coordinates_field(*mesh), krino_lbfgs, smoothingObj);
  }
}

template <typename MESHES> void smooth_both_with_mean_ratio_objective(MESHES & meshes)
{
  const MeanRatioNormObjective elemObj;
  for (auto * mesh : {&meshes.linear.get_mesh(), &meshes.quadratic.get_mesh()})
    smooth_mesh_by_optimization(*mesh,
        coordinates_field(*mesh),
        mesh->mesh_meta_data().universal_part(),
        krino_lbfgs,
        elemObj,
        boundary_of_box_filter(),
        SolutionProjector{});
}

void krino_lbfgs_3d(const Objective3DInterface & objFn, stk::math::Vector3d & x)
{
  const double xTol = 1.e-4;
  lbfgs(objFn, x, xTol);
}

void smooth_both_node_at_a_time_with_mean_ratio(CubeOf24Tet4AndTet10 & meshes)
{
  const MeanRatioAboutNodeObjective nodeObj;
  for (auto * mesh : {&meshes.linear.get_mesh(), &meshes.quadratic.get_mesh()})
    smooth_mesh_by_optimizing_node_locations(*mesh,
        coordinates_field(*mesh),
        mesh->mesh_meta_data().universal_part(),
        krino_lbfgs_3d,
        nodeObj,
        boundary_of_box_filter());
}

template <typename MESH_AND_BUILDER>
void smooth_interior_nodes_by_ODT(MESH_AND_BUILDER & meshAndBuilder)
{
  stk::mesh::BulkData & mesh = meshAndBuilder.get_mesh();
  const stk::mesh::Selector elemSelector(meshAndBuilder.get_aux_meta().active_part());
  const stk::mesh::Selector interiorNodeSelector = elemSelector &
      !meshAndBuilder.get_aux_meta().exposed_boundary_part() &
      !meshAndBuilder.get_aux_meta().block_boundary_part();
  smooth_mesh_by_ODT_of_interior_node_locations(
      mesh, coordinates_field(mesh), elemSelector, interiorNodeSelector);
}

} // namespace

// The mesh objective of the tet10 mesh has only its vertex nodes, and the same value and gradient
// as that of the tet4 mesh.
TEST(QuadraticSmoothing, tet10MeshObjective_onlyVertexNodesAndSameValueAsTet4)
{
  CubeOf24Tet4AndTet10 meshes;
  displace_cube_face_and_center_nodes(meshes);

  const std::shared_ptr<KinematicElementObjective> elemObjPtr =
      make_modified_neo_hookean_objective();
  const KinematicElementObjective & elemObj = *elemObjPtr;
  const stk::mesh::BulkData & linearMesh = meshes.linear.get_mesh();
  const stk::mesh::BulkData & quadraticMesh = meshes.quadratic.get_mesh();
  KinematicConstantIsotropicMeshObjective linearObj(linearMesh,
      coordinates_field(linearMesh),
      linearMesh.mesh_meta_data().universal_part(),
      elemObj,
      1.0);
  KinematicConstantIsotropicMeshObjective quadraticObj(quadraticMesh,
      coordinates_field(quadraticMesh),
      quadraticMesh.mesh_meta_data().universal_part(),
      elemObj,
      1.0);

  ASSERT_EQ(linearObj.get_nodes().size(), quadraticObj.get_nodes().size());
  for (size_t n = 0; n < linearObj.get_nodes().size(); ++n)
  {
    EXPECT_EQ(linearMesh.identifier(linearObj.get_nodes()[n]),
        quadraticMesh.identifier(quadraticObj.get_nodes()[n]));
  }

  const DistributedVector linearX =
      gather_node_coordinates(linearMesh, coordinates_field(linearMesh), linearObj.get_nodes());
  const DistributedVector quadraticX = gather_node_coordinates(
      quadraticMesh, coordinates_field(quadraticMesh), quadraticObj.get_nodes());
  EXPECT_EQ(linearObj.compute_value(linearX), quadraticObj.compute_value(quadraticX));
  DistributedVector linearGrad, quadraticGrad;
  linearObj.fill_gradient(linearX, linearGrad);
  quadraticObj.fill_gradient(quadraticX, quadraticGrad);
  ASSERT_EQ(linearGrad.size(), quadraticGrad.size());
  for (size_t i = 0; i < linearGrad.size(); ++i)
  {
    EXPECT_EQ(linearGrad[i], quadraticGrad[i]) << i;
  }
}

TEST(
    QuadraticSmoothing, tet10WithKinematicObjective_vertexNodesAsTet4AndMidEdgeNodesOfMovedEdgesAtEdgeMidpoints)
{
  CubeOf24Tet4AndTet10 meshes;
  displace_cube_face_and_center_nodes(meshes);
  smooth_both_with_kinematic_objective(meshes);
  EXPECT_TRUE(meshes.is_linear_vertex_node_restored(14, 1.e-2)); // the center node, which was moved
  meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  meshes.expect_some_edges_with_vertex_nodes_not_moved_by_smoothing();
}

TEST(
    QuadraticSmoothing, tet10WithMeanRatioObjective_vertexNodesAsTet4AndMidEdgeNodesOfMovedEdgesAtEdgeMidpoints)
{
  CubeOf24Tet4AndTet10 meshes;
  displace_cube_face_and_center_nodes(meshes);
  smooth_both_with_mean_ratio_objective(meshes);
  meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  meshes.expect_some_edges_with_vertex_nodes_not_moved_by_smoothing();
}

TEST(QuadraticSmoothing,
    tri6WithKinematicAndMeanRatioObjectives_vertexNodesAsTri3AndMidEdgeNodesOfMovedEdgesAtEdgeMidpoints)
{
  const size_t centerNode = 8;
  {
    SquareOf8Tri3AndTri6 meshes;
    meshes.displace_vertex_node(centerNode, stk::math::Vector3d(0.3, -0.4, 0.));
    meshes.displace_mid_edge_nodes_off_edge_midpoints();
    meshes.sync_coordinates();
    meshes.record_quadratic_node_locations_before_smoothing();
    smooth_both_with_kinematic_objective(meshes);
    EXPECT_TRUE(meshes.is_linear_vertex_node_restored(
        centerNode, 1.e-2)); // the center node, which was moved
    meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  }
  {
    SquareOf8Tri3AndTri6 meshes;
    meshes.displace_vertex_node(centerNode, stk::math::Vector3d(0.3, -0.4, 0.));
    meshes.displace_mid_edge_nodes_off_edge_midpoints();
    meshes.sync_coordinates();
    meshes.record_quadratic_node_locations_before_smoothing();
    smooth_both_with_mean_ratio_objective(meshes);
    meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  }
}

// The element reference sizes and the metric mappings from a nodal metric field depend only on the
// vertex nodes, so they are the same for the tet10 mesh as for the tet4 mesh.
TEST(QuadraticSmoothing, tet10RefSizesAndNodalMetricMappings_sameAsTet4)
{
  CubeOf24Tet4AndTet10 meshes;
  displace_cube_face_and_center_nodes(meshes);

  const auto & linearElems = meshes.linear.get_owned_elements();
  const auto & quadraticElems = meshes.quadratic.get_owned_elements();
  const FieldRef linearCoords(meshes.linear.get_mesh().mesh_meta_data().coordinate_field());
  const FieldRef quadraticCoords(meshes.quadratic.get_mesh().mesh_meta_data().coordinate_field());
  EXPECT_EQ(compute_element_ref_sizes_from_reference(
                meshes.linear.get_mesh(), linearCoords, linearElems, SmoothReference::EqualVolume),
      compute_element_ref_sizes_from_reference(meshes.quadratic.get_mesh(),
          quadraticCoords,
          quadraticElems,
          SmoothReference::EqualVolume));

  const FieldRef linearMetric =
      CubeOf24Tet4AndTet10::set_varying_nodal_metric(meshes.linear.get_mesh());
  const FieldRef quadraticMetric =
      CubeOf24Tet4AndTet10::set_varying_nodal_metric(meshes.quadratic.get_mesh());
  EXPECT_EQ(compute_element_metric_mappings_from_nodal_metric_field(
                meshes.linear.get_mesh(), linearMetric, linearElems),
      compute_element_metric_mappings_from_nodal_metric_field(
          meshes.quadratic.get_mesh(), quadraticMetric, quadraticElems));
}

// Node-at-a-time smoothing moves only the vertex nodes of the tet10 mesh, as for the tet4 mesh, and
// then moves its mid-edge nodes to the midpoints of their edges.
TEST(QuadraticSmoothing,
    tet10NodeAtATimeWithMeanRatio_vertexNodesAsTet4AndMidEdgeNodesOfMovedEdgesAtEdgeMidpoints)
{
  CubeOf24Tet4AndTet10 meshes;
  displace_cube_face_and_center_nodes(meshes);
  smooth_both_node_at_a_time_with_mean_ratio(meshes);
  EXPECT_TRUE(meshes.is_linear_vertex_node_restored(14, 1.e-2)); // the center node, which was moved
  meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  meshes.expect_some_edges_with_vertex_nodes_not_moved_by_smoothing();
}

// Each mid-edge node is moved by its owner, which need not own an element with its edge.
TEST(QuadraticSmoothing,
    tet10WithMidEdgeNodeOwnedByProcWithoutElementOnItsEdge_vertexNodesAsTet4AndMidEdgeNodesOfMovedEdgesAtEdgeMidpoints)
{
  const std::vector<int> halfOfElementsOnEachProc{
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  CubeOf24Tet4AndTet10 meshes({{2, halfOfElementsOnEachProc}});
  if (meshes.quadratic.parallel_size() == 2)
  {
    const size_t firstFaceOrCenterNodeIndex = 8;
    ASSERT_EQ(1, meshes.give_unshared_mid_edge_node_to_proc_1(firstFaceOrCenterNodeIndex));
    displace_cube_face_and_center_nodes(meshes);
    smooth_both_with_mean_ratio_objective(meshes);
    meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  }
}

TEST(QuadraticSmoothing, tet10NodeAtATimeWithODT_vertexNodesAsTet4AndMidEdgeNodesOfMovedEdgesAtEdgeMidpoints)
{
  CubeOf24Tet4AndTet10 meshes;
  meshes.displace_vertex_node(
      14, stk::math::Vector3d(0.3, -0.3, 0.3)); // the center node, the one interior node
  meshes.displace_mid_edge_nodes_off_edge_midpoints();
  meshes.sync_coordinates();
  meshes.record_quadratic_node_locations_before_smoothing();
  smooth_interior_nodes_by_ODT(meshes.linear);
  smooth_interior_nodes_by_ODT(meshes.quadratic);
  EXPECT_TRUE(meshes.is_linear_vertex_node_restored(14, 1.e-2));
  meshes.expect_quadratic_mesh_smoothed_as_linear_with_mid_edge_nodes_of_moved_edges_at_edge_midpoints();
  meshes.expect_some_edges_with_vertex_nodes_not_moved_by_smoothing();
}

// Quadratic elements reach the kinematic objectives as their vertices, so an objective given all
// the nodes of a tet10 or tri6 rejects them rather than reading only some.
TEST(QuadraticSmoothing, kinematicObjectiveGivenAllQuadraticNodes_throws)
{
  const std::shared_ptr<KinematicElementObjective> elemObj = make_modified_neo_hookean_objective();
  const std::array<double, 9> identity3d{{1., 0., 0., 0., 1., 0., 0., 0., 1.}};
  const std::array<double, 4> identity2d{{1., 0., 0., 1.}};
  const std::array<double, 3> x{{0., 0., 0.}};
  const std::array<const double *, 10> nodeCoords{{x.data(),
      x.data(),
      x.data(),
      x.data(),
      x.data(),
      x.data(),
      x.data(),
      x.data(),
      x.data(),
      x.data()}};
  std::vector<double> grad;
  EXPECT_ANY_THROW(elemObj->compute_element_objective(3, 10, identity3d.data(), nodeCoords.data()));
  EXPECT_ANY_THROW(
      elemObj->fill_element_sensitivity(3, 10, identity3d.data(), nodeCoords.data(), grad));
  EXPECT_ANY_THROW(elemObj->compute_element_objective(2, 6, identity2d.data(), nodeCoords.data()));
  EXPECT_ANY_THROW(
      elemObj->fill_element_sensitivity(2, 6, identity2d.data(), nodeCoords.data(), grad));
}

} // namespace krino

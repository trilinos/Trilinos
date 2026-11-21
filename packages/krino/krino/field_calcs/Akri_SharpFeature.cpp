#include <Akri_SharpFeature.hpp>

#include <array>
#include <type_traits>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Edge.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_SharpFeatureCalcs.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

static std::array<Edge,3> get_tri_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element)
{
  StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  return { edge_from_edge_nodes(mesh, elementNodes[0],elementNodes[1]),
           edge_from_edge_nodes(mesh, elementNodes[1],elementNodes[2]),
           edge_from_edge_nodes(mesh, elementNodes[2],elementNodes[0]) };
}

static void fill_tet_face_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity face, std::vector<Edge> & sideEdges)
{
  const std::array<Edge,3> triEdges = get_tri_edges(mesh, face);
  sideEdges.assign(triEdges.begin(), triEdges.end());
}

static bool does_entity_have_selected_element(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, const stk::mesh::Selector & elementSelector)
{
  for (auto && elem : StkMeshEntities{mesh.begin_elements(entity), mesh.end_elements(entity)})
    if (elementSelector(mesh.bucket(elem)))
      return true;
  return false;
}

static bool is_node_to_constrain(
    const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const bool includeSharedNodes)
{
  return mesh.parallel_rank() == mesh.parallel_owner_rank(node) ||
      (includeSharedNodes && mesh.bucket(node).shared());
}

static bool edge_has_node_to_constrain(
    const stk::mesh::BulkData & mesh, const Edge edge, const bool includeSharedNodes)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  return is_node_to_constrain(mesh, edgeNodes[0], includeSharedNodes) ||
      is_node_to_constrain(mesh, edgeNodes[1], includeSharedNodes);
}

static std::vector<Edge>
get_edges_with_nodes_to_constrain_of_selected_faces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const bool includeSharedNodes)
{
  std::vector<Edge> edges;
  std::vector<Edge> sideEdges;
  for(const auto & bucketPtr : mesh.buckets(stk::topology::FACE_RANK))
  {
    if (sideSelector(*bucketPtr))
    {
      for(const auto & side : *bucketPtr)
      {
        if (does_entity_have_selected_element(mesh, side, elementSelector))
        {
          fill_tet_face_edges(mesh, side, sideEdges);
          for (auto edge : sideEdges)
            if (edge_has_node_to_constrain(mesh, edge, includeSharedNodes)) edges.push_back(edge);
        }
      }
    }
  }

  stk::util::sort_and_unique(edges);

  return edges;
}

static std::vector<stk::mesh::Entity>
get_nodes_to_constrain_of_edges_with_selected_sides(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const bool includeSharedNodes)
{
  std::vector<stk::mesh::Entity> edgeNodes;
  for(const auto & bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
    if ((bucketPtr->owned() || (includeSharedNodes && bucketPtr->shared())) &&
        sideSelector(*bucketPtr))
      for(const auto & node : *bucketPtr)
        if (does_entity_have_selected_element(mesh, node, elementSelector))
          edgeNodes.push_back(node);
  return edgeNodes;
}

bool is_intersection_point_node_compatible_for_snapping_based_on_sharp_features(const SharpFeatureInfo & sharpFeatureInfo, const stk::mesh::Entity intPtNode, const std::vector<stk::mesh::Entity> & intPtNodes)
{
  const SharpFeatureConstraint * constraint = sharpFeatureInfo.get_constraint(intPtNode);

  if (constraint == nullptr) return true;
  if (intPtNodes.size() != 2 || constraint->is_pinned()) return false;

  const std::array<stk::mesh::Entity,2> sharpEdgeNodes = constraint->get_sharp_edge_nodes();
  for (auto && sharpEdgeNode : sharpEdgeNodes)
    if (intPtNodes[0] == sharpEdgeNode || intPtNodes[1] == sharpEdgeNode)
      return true;
  return false;
}

stk::math::Vector3d SharpFeatureInfo::NodeCoordinates::operator()(
    const stk::mesh::Entity node, const unsigned dim) const
{
  const stk::math::Vector3d x(coordsField.data<double>().entity_values(node), dim);
  if (!displacements.valid()) return x;
  return x - stk::math::Vector3d(displacements.data<double>().entity_values(node), dim);
}

void SharpFeatureInfo::find_sharp_features(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const double cosFeatureAngle,
    const bool includeSharedNodes,
    const FieldRef coordsDisplacements)
{
  const NodeCoordinates coords{coordsField, coordsDisplacements};

  if (mesh.mesh_meta_data().spatial_dimension() == 2)
    find_sharp_features_2D(
        mesh, coords, elementSelector, sideSelector, cosFeatureAngle, includeSharedNodes);
  find_sharp_features_3D(
      mesh, coords, elementSelector, sideSelector, cosFeatureAngle, includeSharedNodes);

}

void SharpFeatureInfo::find_sharp_features_3D(const stk::mesh::BulkData & mesh,
    const NodeCoordinates & coords,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const double cosFeatureAngle,
    const bool includeSharedNodes)
{
  std::map<stk::mesh::Entity,std::vector<stk::mesh::Entity>> nodeToSharpEdgeNeighbors;

  const std::vector<Edge> edgesWithNodesToConstrain =
      get_edges_with_nodes_to_constrain_of_selected_faces(
          mesh, elementSelector, sideSelector, includeSharedNodes);
  for (auto edge : edgesWithNodesToConstrain)
  {
    if (edge_has_sharp_feature_3D(
            mesh, coords, elementSelector, sideSelector, cosFeatureAngle, edge))
    {
      const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
      if (is_node_to_constrain(mesh, edgeNodes[0], includeSharedNodes))
        nodeToSharpEdgeNeighbors[edgeNodes[0]].push_back(edgeNodes[1]);
      if (is_node_to_constrain(mesh, edgeNodes[1], includeSharedNodes))
        nodeToSharpEdgeNeighbors[edgeNodes[1]].push_back(edgeNodes[0]);
    }
  }

  for (auto && entry : nodeToSharpEdgeNeighbors)
    if (entry.second.size() == 2)
      myNodeToConstrainedNeighbors.insert({entry.first, SharpFeatureConstraint::edge_constraint(entry.second[0], entry.second[1])});
    else if (entry.second.size() > 2)
      myNodeToConstrainedNeighbors.insert({entry.first, SharpFeatureConstraint::pinned_constraint()});
}

void SharpFeatureInfo::find_sharp_features_2D(const stk::mesh::BulkData & mesh,
    const NodeCoordinates & coords,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const double cosFeatureAngle,
    const bool includeSharedNodes)
{
  std::map<stk::mesh::Entity,std::vector<stk::mesh::Entity>> nodeToSharpEdgeNeighbors;

  const std::vector<stk::mesh::Entity> sideNodes =
      get_nodes_to_constrain_of_edges_with_selected_sides(
          mesh, elementSelector, sideSelector, includeSharedNodes);
  for (auto node : sideNodes)
    if (node_has_sharp_feature_2D(
            mesh, coords, elementSelector, sideSelector, cosFeatureAngle, node))
      myNodeToConstrainedNeighbors.insert({node, SharpFeatureConstraint::pinned_constraint()});
}

static void filter_sides_based_on_attached_element_and_side_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, std::vector<stk::mesh::Entity> & sides)
{
  size_t numRetainedSides = 0;
  for (auto && side : sides)
    if (sideSelector(mesh.bucket(side)) && does_entity_have_selected_element(mesh, side, elementSelector))
      sides[numRetainedSides++] = side;
  sides.resize(numRetainedSides);
}

bool SharpFeatureInfo::angle_is_sharp_between_any_two_sides_2D(const stk::mesh::BulkData & mesh,
    const NodeCoordinates & coords,
    const double cosFeatureAngle,
    const stk::mesh::Entity node,
    const std::vector<stk::mesh::Entity> & sidesOfEdge)
{
  if (sidesOfEdge.size() > 1)
  {
    const stk::math::Vector3d nodeCoords = coords(node, 2);

    std::vector<stk::math::Vector3d> otherSideNodeCoords;
    otherSideNodeCoords.reserve(sidesOfEdge.size());
    for (auto && side : sidesOfEdge)
    {
      StkMeshEntities sideNodes{mesh.begin_nodes(side), mesh.end_nodes(side)};
      STK_ThrowAssertMsg(sideNodes[0] == node || sideNodes[1] == node, "Did not find side node for segment.");
      const stk::mesh::Entity sideNode = (sideNodes[1] == node) ? sideNodes[0] : sideNodes[1];
      otherSideNodeCoords.push_back(coords(sideNode, 2));
    }

    for (size_t i = 0; i < otherSideNodeCoords.size(); ++i)
      for (size_t j = i + 1; j < otherSideNodeCoords.size(); ++j)
        if (is_angle_sharper_than_feature_angle(
                cosine_of_angle_between_segments_at_vertex(
                    nodeCoords, otherSideNodeCoords[i], otherSideNodeCoords[j]),
                cosFeatureAngle))
          return true;
  }
  return false;
}

static stk::mesh::Entity get_face_node_not_on_edge(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,2> & edgeNodes, const stk::mesh::Entity sideOfEdge)
{
  StkMeshEntities faceNodes{mesh.begin_nodes(sideOfEdge), mesh.end_nodes(sideOfEdge)};
  for (auto && faceNode : faceNodes)
    if (faceNode != edgeNodes[0] && faceNode != edgeNodes[1])
      return faceNode;
  STK_ThrowErrorMsg("Did not find face node not on edge.");
  return stk::mesh::Entity(); // unreachable; STK_ThrowErrorMsg does not return
}

bool is_sharp_between_tri_sides(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const std::array<stk::mesh::Entity, 2> & edgeNodes,
    const stk::mesh::Entity side0,
    const stk::mesh::Entity side1,
    const double cosFeatureAngle)
{
  const SharpFeatureInfo::NodeCoordinates coords{coordsField, FieldRef()};
  return is_angle_sharper_than_feature_angle(
      cosine_of_angle_between_triangles_at_edge(coords(edgeNodes[0], 3),
          coords(edgeNodes[1], 3),
          coords(get_face_node_not_on_edge(mesh, edgeNodes, side0), 3),
          coords(get_face_node_not_on_edge(mesh, edgeNodes, side1), 3)),
      cosFeatureAngle);
}

static stk::mesh::Entity get_other_node_of_line_side(
    const stk::mesh::BulkData & mesh, const stk::mesh::Entity node, const stk::mesh::Entity side)
{
  StkMeshEntities sideNodes{mesh.begin_nodes(side), mesh.end_nodes(side)};
  STK_ThrowAssertMsg(
      sideNodes[0] == node || sideNodes[1] == node, "Did not find side node for segment.");
  return (sideNodes[1] == node) ? sideNodes[0] : sideNodes[1];
}

bool is_sharp_between_line_sides(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Entity node,
    const stk::mesh::Entity side0,
    const stk::mesh::Entity side1,
    const double cosFeatureAngle)
{
  const SharpFeatureInfo::NodeCoordinates coords{coordsField, FieldRef()};
  return is_angle_sharper_than_feature_angle(
      cosine_of_angle_between_segments_at_vertex(coords(node, 2),
          coords(get_other_node_of_line_side(mesh, node, side0), 2),
          coords(get_other_node_of_line_side(mesh, node, side1), 2)),
      cosFeatureAngle);
}

bool SharpFeatureInfo::angle_is_sharp_between_any_two_sides_3D(const stk::mesh::BulkData & mesh,
    const NodeCoordinates & coords,
    const double cosFeatureAngle,
    const std::array<stk::mesh::Entity, 2> & edgeNodes,
    const std::vector<stk::mesh::Entity> & sidesOfEdge)
{
  if (sidesOfEdge.size() > 1)
  {
    const stk::math::Vector3d edgeNodeCoords0 = coords(edgeNodes[0], 3);
    const stk::math::Vector3d edgeNodeCoords1 = coords(edgeNodes[1], 3);

    std::vector<stk::math::Vector3d> coordsOfNonEdgeNodeOfSides;
    coordsOfNonEdgeNodeOfSides.reserve(sidesOfEdge.size());
    for (auto && side : sidesOfEdge)
      coordsOfNonEdgeNodeOfSides.push_back(
          coords(get_face_node_not_on_edge(mesh, edgeNodes, side), 3));

    for (size_t i = 0; i < coordsOfNonEdgeNodeOfSides.size(); ++i)
      for (size_t j = i + 1; j < coordsOfNonEdgeNodeOfSides.size(); ++j)
        if (is_angle_sharper_than_feature_angle(
                cosine_of_angle_between_triangles_at_edge(edgeNodeCoords0,
                    edgeNodeCoords1,
                    coordsOfNonEdgeNodeOfSides[i],
                    coordsOfNonEdgeNodeOfSides[j]),
                cosFeatureAngle))
          return true;
  }
  return false;
}

bool SharpFeatureInfo::edge_has_sharp_feature_3D(const stk::mesh::BulkData & mesh,
    const NodeCoordinates & coords,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const double cosFeatureAngle,
    const Edge edge)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  std::vector<stk::mesh::Entity> sidesOfEdge;
  stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edgeNodes[0], edgeNodes[1]}, stk::topology::FACE_RANK, sidesOfEdge);
  if (sidesOfEdge.size() > 1)
    filter_sides_based_on_attached_element_and_side_parts(mesh, elementSelector, sideSelector, sidesOfEdge);
  return angle_is_sharp_between_any_two_sides_3D(
      mesh, coords, cosFeatureAngle, edgeNodes, sidesOfEdge);
}

bool SharpFeatureInfo::node_has_sharp_feature_2D(const stk::mesh::BulkData & mesh,
    const NodeCoordinates & coords,
    const stk::mesh::Selector & elementSelector,
    const stk::mesh::Selector & sideSelector,
    const double cosFeatureAngle,
    const stk::mesh::Entity node)
{
  std::vector<stk::mesh::Entity> sidesOfEdge(mesh.begin_edges(node), mesh.end_edges(node));
  if (sidesOfEdge.size() > 1)
    filter_sides_based_on_attached_element_and_side_parts(mesh, elementSelector, sideSelector, sidesOfEdge);
  return angle_is_sharp_between_any_two_sides_2D(mesh, coords, cosFeatureAngle, node, sidesOfEdge);
}

const SharpFeatureConstraint * SharpFeatureInfo::get_constraint(const stk::mesh::Entity node) const
{
  const auto iter = myNodeToConstrainedNeighbors.find(node);
  if (iter != myNodeToConstrainedNeighbors.end())
    return &(iter->second);
  return nullptr;
}

}

#include "Akri_SharpFeature.hpp"

#include <array>
#include <type_traits>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_MeshHelpers.hpp>
#include "Akri_AuxMetaData.hpp"
#include "Akri_CDMesh_Utils.hpp"
#include "Akri_Edge.hpp"
#include "Akri_Phase_Support.hpp"
#include <stk_math/StkVector.hpp>

namespace krino {

std::array<Edge,6> get_tet_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element)
{
  StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  return { edge_from_edge_nodes(mesh, elementNodes[0],elementNodes[1]),
           edge_from_edge_nodes(mesh, elementNodes[1],elementNodes[2]),
           edge_from_edge_nodes(mesh, elementNodes[2],elementNodes[0]),
           edge_from_edge_nodes(mesh, elementNodes[3],elementNodes[0]),
           edge_from_edge_nodes(mesh, elementNodes[3],elementNodes[1]),
           edge_from_edge_nodes(mesh, elementNodes[3],elementNodes[2]) };
}

std::array<Edge,3> get_tri_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element)
{
  StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  return { edge_from_edge_nodes(mesh, elementNodes[0],elementNodes[1]),
           edge_from_edge_nodes(mesh, elementNodes[1],elementNodes[2]),
           edge_from_edge_nodes(mesh, elementNodes[2],elementNodes[0]) };
}

Edge get_segment_edge(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element)
{
  StkMeshEntities elementNodes{mesh.begin_nodes(element), mesh.end_nodes(element)};
  return edge_from_edge_nodes(mesh, elementNodes[0],elementNodes[1]);
}

void fill_element_edges(const stk::mesh::BulkData & mesh, const unsigned dim, const stk::mesh::Entity element, std::vector<Edge> & elementEdges)
{
  if (dim == 2)
  {
    const std::array<Edge,3> triEdges = get_tri_edges(mesh, element);
    elementEdges.assign(triEdges.begin(), triEdges.end());
    return;
  }

  const std::array<Edge,6> tetEdges = get_tet_edges(mesh, element);
  elementEdges.assign(tetEdges.begin(), tetEdges.end());
}

void fill_tet_face_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity face, std::vector<Edge> & sideEdges)
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

static stk::mesh::Selector
build_side_selector(const stk::mesh::BulkData & mesh)
{
  const AuxMetaData & auxMeta = AuxMetaData::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const stk::mesh::EntityRank sideRank = mesh.mesh_meta_data().side_rank();

  stk::mesh::PartVector sideParts;
  for (auto && part : mesh.mesh_meta_data().get_parts())
    if (is_part_to_check_for_snapping_compatibility(phaseSupport, auxMeta, sideRank, *part))
      sideParts.push_back(part);

  return stk::mesh::selectUnion(sideParts);
}

bool edge_has_owned_node(const stk::mesh::BulkData & mesh, const Edge edge)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  return mesh.parallel_rank() == mesh.parallel_owner_rank(edgeNodes[0]) || mesh.parallel_rank() == mesh.parallel_owner_rank(edgeNodes[1]);
}

std::vector<Edge> get_edges_with_owned_nodes_of_selected_faces(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector)
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
            if (edge_has_owned_node(mesh, edge))
              edges.push_back(edge);
        }
      }
    }
  }

  stk::util::sort_and_unique(edges);

  return edges;
}

std::vector<stk::mesh::Entity> get_owned_nodes_of_edges_with_selected_sides(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector)
{
  std::vector<stk::mesh::Entity> edgeNodes;
  for(const auto & bucketPtr : mesh.buckets(stk::topology::NODE_RANK))
    if (bucketPtr->owned() && sideSelector(*bucketPtr))
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

void SharpFeatureInfo::find_sharp_features(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const double cosFeatureAngle)
{
  const stk::mesh::Selector sideSelector = build_side_selector(mesh);

  if (mesh.mesh_meta_data().spatial_dimension() == 2)
    find_sharp_features_2D(mesh, coordsField, elementSelector, sideSelector, cosFeatureAngle);
  find_sharp_features_3D(mesh, coordsField, elementSelector, sideSelector, cosFeatureAngle);

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    for (auto && entry : myNodeToConstrainedNeighbors)
    {
      stk::mesh::Entity node = entry.first;
      const SharpFeatureConstraint & constraint = entry.second;
      krinolog << "Node " << mesh.identifier(node) << " is ";
      if (constraint.is_pinned())
      {
        krinolog << "pinned." << stk::diag::dendl;
      }
      else
      {
        const std::array<stk::mesh::Entity,2> sharpEdgeNbrs = constraint.get_sharp_edge_nodes();
        krinolog << "constrained to move along edge between nodes " << mesh.identifier(sharpEdgeNbrs[0]) << " and " << mesh.identifier(sharpEdgeNbrs[1]) << "." << std::endl;
      }
    }
  }
}

void SharpFeatureInfo::find_sharp_features_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle)
{
  std::map<stk::mesh::Entity,std::vector<stk::mesh::Entity>> nodeToSharpEdgeNeighbors;
  const int parallelRank = mesh.parallel_rank();

  const std::vector<Edge> edgesWithOwnedNodes = get_edges_with_owned_nodes_of_selected_faces(mesh, elementSelector, sideSelector);
  for (auto edge : edgesWithOwnedNodes)
  {
    if (edge_has_sharp_feature_3D(mesh, coordsField, elementSelector, sideSelector, cosFeatureAngle, edge))
    {
      const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
      if (parallelRank == mesh.parallel_owner_rank(edgeNodes[0]))
        nodeToSharpEdgeNeighbors[edgeNodes[0]].push_back(edgeNodes[1]);
      if (parallelRank == mesh.parallel_owner_rank(edgeNodes[1]))
        nodeToSharpEdgeNeighbors[edgeNodes[1]].push_back(edgeNodes[0]);
    }
  }

  for (auto && entry : nodeToSharpEdgeNeighbors)
    if (entry.second.size() == 2)
      myNodeToConstrainedNeighbors.insert({entry.first, SharpFeatureConstraint::edge_constraint(entry.second[0], entry.second[1])});
    else if (entry.second.size() > 2)
      myNodeToConstrainedNeighbors.insert({entry.first, SharpFeatureConstraint::pinned_constraint()});
}

void SharpFeatureInfo::find_sharp_features_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle)
{
  std::map<stk::mesh::Entity,std::vector<stk::mesh::Entity>> nodeToSharpEdgeNeighbors;

  const std::vector<stk::mesh::Entity> ownedSideNodes = get_owned_nodes_of_edges_with_selected_sides(mesh, elementSelector, sideSelector);
  for (auto node : ownedSideNodes)
    if (node_has_sharp_feature_2D(mesh, coordsField, elementSelector, sideSelector, cosFeatureAngle, node))
      myNodeToConstrainedNeighbors.insert({node, SharpFeatureConstraint::pinned_constraint()});
}

void filter_sides_based_on_attached_element_and_side_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, std::vector<stk::mesh::Entity> & sides)
{
  size_t numRetainedSides = 0;
  for (auto && side : sides)
    if (sideSelector(mesh.bucket(side)) && does_entity_have_selected_element(mesh, side, elementSelector))
      sides[numRetainedSides++] = side;
  sides.resize(numRetainedSides);
}

bool SharpFeatureInfo::angle_is_sharp_between_any_two_sides_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double cosFeatureAngle, const stk::mesh::Entity node, const std::vector<stk::mesh::Entity> & sidesOfEdge)
{
  if (sidesOfEdge.size() > 1)
  {
    const stk::math::Vector3d nodeCoords(field_data<double>(coordsField, node),2);

    std::vector<stk::math::Vector3d> sideVec;
    sideVec.reserve(sidesOfEdge.size());
    for (auto && side : sidesOfEdge)
    {
      StkMeshEntities sideNodes{mesh.begin_nodes(side), mesh.end_nodes(side)};
      STK_ThrowAssertMsg(sideNodes[0] == node || sideNodes[1] == node, "Did not find side node for segment.");
      const stk::mesh::Entity sideNode = (sideNodes[1] == node) ? sideNodes[0] : sideNodes[1];
      const stk::math::Vector3d coordsOfSideNode(field_data<double>(coordsField, sideNode),2);
      sideVec.push_back((coordsOfSideNode - nodeCoords).unit_vector());
    }

    for (size_t i=0; i<sideVec.size(); ++i)
      for (size_t j=i+1; j<sideVec.size(); ++j)
        if (Dot(sideVec[i], sideVec[j]) > cosFeatureAngle)
          return true;
  }
  return false;
}

double cosine_of_dihedral_angle_3D(const stk::math::Vector3d & edgeVec, const stk::math::Vector3d & faceTangent0, const stk::math::Vector3d & faceTangent1)
{
  // https://en.wikipedia.org/wiki/Dihedral_angle
  const stk::math::Vector3d crossEdgeFace0 = Cross(edgeVec, faceTangent0);
  const stk::math::Vector3d crossEdgeFace1 = Cross(edgeVec, faceTangent1);

  return Dot(crossEdgeFace0,crossEdgeFace1) / (crossEdgeFace0.length()*crossEdgeFace1.length());
}

stk::mesh::Entity get_face_node_not_on_edge(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,2> & edgeNodes, const stk::mesh::Entity sideOfEdge)
{
  StkMeshEntities faceNodes{mesh.begin_nodes(sideOfEdge), mesh.end_nodes(sideOfEdge)};
  for (auto && faceNode : faceNodes)
    if (faceNode != edgeNodes[0] && faceNode != edgeNodes[1])
      return faceNode;
  ThrowRuntimeError("Did not find face node not on edge.");
}

bool SharpFeatureInfo::angle_is_sharp_between_any_two_sides_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const double cosFeatureAngle, const std::array<stk::mesh::Entity,2> & edgeNodes, const std::vector<stk::mesh::Entity> & sidesOfEdge)
{
  if (sidesOfEdge.size() > 1)
  {
    const stk::math::Vector3d edgeNodeCoords0(field_data<double>(coordsField, edgeNodes[0]));
    const stk::math::Vector3d edgeNodeCoords1(field_data<double>(coordsField, edgeNodes[1]));
    const stk::math::Vector3d edgeVec = edgeNodeCoords1 - edgeNodeCoords0;

    std::vector<stk::math::Vector3d> faceTangent;
    faceTangent.reserve(sidesOfEdge.size());
    for (auto && side : sidesOfEdge)
    {
      const stk::math::Vector3d coordsOfNonEdgeNodeOfSide(field_data<double>(coordsField, get_face_node_not_on_edge(mesh, edgeNodes, side)));
      faceTangent.push_back(coordsOfNonEdgeNodeOfSide - edgeNodeCoords0);
    }

    for (size_t i=0; i<faceTangent.size(); ++i)
      for (size_t j=i+1; j<faceTangent.size(); ++j)
        if (cosine_of_dihedral_angle_3D(edgeVec, faceTangent[i], faceTangent[j]) > cosFeatureAngle)
          return true;
  }
  return false;
}

bool SharpFeatureInfo::edge_has_sharp_feature_3D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle, const Edge edge)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  std::vector<stk::mesh::Entity> sidesOfEdge;
  stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edgeNodes[0], edgeNodes[1]}, stk::topology::FACE_RANK, sidesOfEdge);
  if (sidesOfEdge.size() > 1)
    filter_sides_based_on_attached_element_and_side_parts(mesh, elementSelector, sideSelector, sidesOfEdge);
  return angle_is_sharp_between_any_two_sides_3D(mesh, coordsField, cosFeatureAngle, edgeNodes, sidesOfEdge);
}

bool SharpFeatureInfo::node_has_sharp_feature_2D(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elementSelector, const stk::mesh::Selector & sideSelector, const double cosFeatureAngle, const stk::mesh::Entity node)
{
  std::vector<stk::mesh::Entity> sidesOfEdge(mesh.begin_edges(node), mesh.end_edges(node));
  if (sidesOfEdge.size() > 1)
    filter_sides_based_on_attached_element_and_side_parts(mesh, elementSelector, sideSelector, sidesOfEdge);
  return angle_is_sharp_between_any_two_sides_2D(mesh, coordsField, cosFeatureAngle, node, sidesOfEdge);
}

const SharpFeatureConstraint * SharpFeatureInfo::get_constraint(const stk::mesh::Entity node) const
{
  const auto iter = myNodeToConstrainedNeighbors.find(node);
  if (iter != myNodeToConstrainedNeighbors.end())
    return &(iter->second);
  return nullptr;
}

}

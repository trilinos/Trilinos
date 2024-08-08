#include <Akri_HexRefiner.hpp>
#include <Akri_Refinement.hpp>

#include <numeric>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "Akri_ParallelErrorMessage.hpp"
#include "Akri_ChildNodeCreator.hpp"
#include "Akri_DiagWriter.hpp"
#include "Akri_Edge.hpp"
#include "Akri_EntityIdPool.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_MOAB_TetRefiner.hpp"
#include "Akri_NodeRefiner.hpp"
#include "Akri_TransitionElementEdgeMarker.hpp"
#include <Akri_QuadRefiner.hpp>
#include "Akri_TriRefiner.hpp"
#include "Akri_NodeRefiner.hpp"

namespace krino {

void Refinement::declare_refinement_parts()
{
  myParentPart = &myMeta.declare_part("Refinement_Parent", stk::topology::ELEMENT_RANK, true);
  myChildPart = &myMeta.declare_part("Refinement_Child", stk::topology::ELEMENT_RANK, true);
  myRefinedEdgeNodePart = &myMeta.declare_part_with_topology("Refinement_Edge_Node", stk::topology::NODE);
  myRefinedQuadFaceNodePart  = &myMeta.declare_part_with_topology("Refinement_QuadFace_Node", stk::topology::NODE);
}

void Refinement::declare_refinement_fields()
{
  if (3 == myMeta.spatial_dimension())
  {
    myChildElementIds8Field = &myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_CHILD_ELEMENT_IDS_8");

    const stk::topology tet4Topology = stk::topology::TETRAHEDRON_4;
    const stk::mesh::Part & tet4TopologyPart = myMeta.get_topology_root_part(tet4Topology);

    const stk::topology hex8Topology = stk::topology::HEXAHEDRON_8;
    const stk::mesh::Part & hex8TopologyPart = myMeta.get_topology_root_part(hex8Topology);

    stk::mesh::put_field_on_mesh(*myChildElementIds8Field, parent_part() & tet4TopologyPart, get_num_children_when_fully_refined(tet4Topology), nullptr);
    stk::mesh::put_field_on_mesh(*myChildElementIds8Field, parent_part() & hex8TopologyPart, get_num_children_when_fully_refined(hex8Topology), nullptr);
  }
  else if (2 == myMeta.spatial_dimension())
  {
    myChildElementIds4Field = &myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_CHILD_ELEMENT_IDS_4");
    const stk::topology tri3Topology = stk::topology::TRIANGLE_3_2D;
    const stk::mesh::Part & tri3TopologyPart = myMeta.get_topology_root_part(tri3Topology);

    const stk::topology quad4Topology = stk::topology::QUADRILATERAL_4_2D;
    const stk::mesh::Part & quad4TopologyPart = myMeta.get_topology_root_part(quad4Topology);

    stk::mesh::put_field_on_mesh(*myChildElementIds4Field, parent_part() & tri3TopologyPart, get_num_children_when_fully_refined(tri3Topology), nullptr);
    stk::mesh::put_field_on_mesh(*myChildElementIds4Field, parent_part() & quad4TopologyPart, get_num_children_when_fully_refined(quad4Topology), nullptr);
  }

  {
    myChildElementIds2Field = &myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_CHILD_ELEMENT_IDS_2");

    const stk::topology beam2Topology = stk::topology::BEAM_2;
    const stk::mesh::Part & beam2TopologyPart = myMeta.get_topology_root_part(beam2Topology);

    stk::mesh::put_field_on_mesh(*myChildElementIds2Field, parent_part() & beam2TopologyPart, get_num_children_when_fully_refined(beam2Topology), nullptr);
  }

  myRefinementLevelField = &myMeta.declare_field<int>(stk::topology::ELEMENT_RANK, "REFINEMENT_LEVEL");
  stk::mesh::put_field_on_mesh(*myRefinementLevelField, myMeta.universal_part(), 1, nullptr); // needed everywhere for restart, otherwise could be child_part

  myParentElementIdField = &myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_PARENT_ELEMENT_ID");
  stk::mesh::put_field_on_mesh(*myParentElementIdField, myMeta.universal_part(), 1, &stk::mesh::InvalidEntityId); // needed everywhere for restart, otherwise could be child_part

  myRefinedEdgeNodeParentIdsField = &myMeta.declare_field<uint64_t>(stk::topology::NODE_RANK, "REFINEMENT_REFINED_EDGE_NODE_PARENTS_IDS");
  stk::mesh::put_field_on_mesh(*myRefinedEdgeNodeParentIdsField, refined_edge_node_part(), 2, nullptr);

  myRefinedQuadFaceNodeParentIdsField = &myMeta.declare_field<uint64_t>(stk::topology::NODE_RANK, "REFINEMENT_REFINED_QUAD_FACE_NODE_PARENTS_IDS");
  stk::mesh::put_field_on_mesh(*myRefinedQuadFaceNodeParentIdsField, refined_quad_face_node_part(), 4, nullptr);

  myOriginatingProcForParentElementField = &myMeta.declare_field<int>(stk::topology::ELEMENT_RANK, "ORIGINATING_PROC_FOR_PARENT_ELEMENT");
  stk::mesh::put_field_on_mesh(*myOriginatingProcForParentElementField, myMeta.universal_part(), 1, nullptr); // needed everywhere for restart, otherwise could be parent_part
}

Refinement::Refinement(stk::mesh::MetaData & meta,
    stk::mesh::Part * activePart,
    const bool force64Bit,
    const bool assert32Bit,
    stk::diag::Timer & parentTimer)
    : myMeta(meta),
      myForce64Bit(force64Bit),
      myAssert32Bit(assert32Bit),
      myNodeRefiner(force64Bit, assert32Bit),
      myEntityIdPool(meta),
      myActivePart(activePart),
      refineTimer(parentTimer),
      unrefineTimer(parentTimer),
      myFixPartsandOwnersTimer("Fix Parts and Owners", parentTimer)
{
  myCoordsField = static_cast<const stk::mesh::Field<double>*>(myMeta.coordinate_field());
  declare_refinement_parts();
  declare_refinement_fields();
}

Refinement::Refinement(
    stk::mesh::MetaData & meta, stk::mesh::Part * activePart, stk::diag::Timer & parentTimer)
    : Refinement(meta, activePart, false, false, parentTimer)
{
}

Refinement::Refinement(stk::mesh::MetaData & meta, stk::diag::Timer & parentTimer)
  : Refinement(meta, nullptr, false, false, parentTimer)
{
}

bool Refinement::is_parent(const stk::mesh::Bucket & bucket) const
{
  return bucket.member(parent_part());
}

bool Refinement::is_parent(const stk::mesh::Entity elem) const
{
  return myMeta.mesh_bulk_data().bucket(elem).member(parent_part());
}

bool Refinement::is_child(const stk::mesh::Bucket & bucket) const
{
  return bucket.member(child_part());
}

bool Refinement::is_child(const stk::mesh::Entity elem) const
{
  return myMeta.mesh_bulk_data().bucket(elem).member(child_part());
}

int Refinement::refinement_level(const stk::mesh::Entity elem) const
{
  if (!is_child(elem))
    return 0;
  STK_ThrowAssertMsg(myRefinementLevelField, "Refinement level field is not defined.");
  const auto * refineLevel = stk::mesh::field_data(*myRefinementLevelField, elem);
  STK_ThrowAssertMsg(refineLevel != nullptr, "Refinement level field missing on " << myMeta.mesh_bulk_data().entity_key(elem));
  return *refineLevel;
}

void Refinement::set_refinement_level(const stk::mesh::Entity elem, const int refinementLevel) const
{
  STK_ThrowAssertMsg(myRefinementLevelField, "Refinement level field is not defined.");
  auto * refineLevel = stk::mesh::field_data(*myRefinementLevelField, elem);
  STK_ThrowAssertMsg(is_child(elem) && refineLevel != nullptr, "Refinement level field missing on " << myMeta.mesh_bulk_data().entity_key(elem));
  *refineLevel = refinementLevel;
}


int Refinement::get_originating_processor_for_parent_element(const stk::mesh::Entity elem) const
{
  STK_ThrowAssertMsg(is_parent(elem), "Call to get_originating_processor_for_parent_element() for non-parent element " << myMeta.mesh_bulk_data().entity_key(elem));
  const auto * originatingProc = stk::mesh::field_data(*myOriginatingProcForParentElementField, elem);
  STK_ThrowAssertMsg(originatingProc != nullptr, "ORIGINATING_PROC_FOR_PARENT_ELEMENT field missing on " << myMeta.mesh_bulk_data().entity_key(elem));
  return *originatingProc;
}

void Refinement::set_originating_processor_for_parent_element(const stk::mesh::Entity elem, const int originatingProc) const
{
  STK_ThrowAssertMsg(is_parent(elem), "Call to set_originating_processor_for_parent_element() for non-parent element " << myMeta.mesh_bulk_data().entity_key(elem));
  auto * originatingProcData = stk::mesh::field_data(*myOriginatingProcForParentElementField, elem);
  STK_ThrowAssertMsg(originatingProcData != nullptr, "ORIGINATING_PROC_FOR_PARENT_ELEMENT field missing on " << myMeta.mesh_bulk_data().entity_key(elem));
  *originatingProcData = originatingProc;
}

unsigned Refinement::get_num_children_when_fully_refined(const stk::topology elementTopology)
{
  switch(elementTopology)
  {
  case stk::topology::BEAM_2:
      return 2;
  case stk::topology::TRIANGLE_3_2D:
  case stk::topology::QUADRILATERAL_4_2D:
      return 4;
  case stk::topology::TETRAHEDRON_4:
    case stk::topology::HEXAHEDRON_8:
      return 8;
  default:
      ThrowRuntimeError("Element topology not found in get_num_children_when_fully_refined: " << elementTopology.name());
      break;
  }
}

unsigned Refinement::get_num_children_when_fully_refined(const stk::mesh::Entity elem) const
{
  return get_num_children_when_fully_refined(myMeta.mesh_bulk_data().bucket(elem).topology());
}

std::array<stk::mesh::EntityId,2> Refinement::get_edge_parent_node_ids(const stk::mesh::Entity edgeNode) const
{
  STK_ThrowAssertMsg(myRefinedEdgeNodeParentIdsField, "Edge Node Ids field is not defined.");
  auto * edgeNodeIds = stk::mesh::field_data(*myRefinedEdgeNodeParentIdsField, edgeNode);
  STK_ThrowAssertMsg(edgeNodeIds != nullptr, "Edge Node Ids field missing on node " << myMeta.mesh_bulk_data().identifier(edgeNode));
  return {{edgeNodeIds[0], edgeNodeIds[1]}};
}

std::array<stk::mesh::Entity,2> Refinement::get_edge_parent_nodes(const stk::mesh::Entity edgeNode) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  STK_ThrowAssertMsg(myRefinedEdgeNodeParentIdsField, "Edge Node Ids field is not defined.");
  auto * edgeNodeIds = stk::mesh::field_data(*myRefinedEdgeNodeParentIdsField, edgeNode);
  STK_ThrowAssertMsg(edgeNodeIds != nullptr, "Edge Node Ids field missing on node " << mesh.identifier(edgeNode));
  return {{mesh.get_entity(stk::topology::NODE_RANK, edgeNodeIds[0]), mesh.get_entity(stk::topology::NODE_RANK, edgeNodeIds[1])}};
}

static unsigned num_children_with_valid_id(const uint64_t * childElemIdsData, const unsigned numChildWhenFullyRefined)
{
  unsigned numChild = 0;
  for (unsigned iChild=0; iChild<numChildWhenFullyRefined; ++iChild)
    if (stk::mesh::InvalidEntityId != childElemIdsData[iChild])
      ++numChild;
  return numChild;
}

static bool has_child_with_invalid_id(const uint64_t * childElemIdsData, const unsigned numChildWhenFullyRefined)
{
  return (stk::mesh::InvalidEntityId == childElemIdsData[numChildWhenFullyRefined-1]);
}

bool Refinement::is_this_parent_element_partially_refined(const stk::mesh::Entity parentElem) const
{
  STK_ThrowAssert(is_parent(parentElem));
  //const auto & [childElemIdsData, numChildWhenFullyRefined] = get_child_ids_and_num_children_when_fully_refined(parentElem);
  const auto & childElemIdsDataAndnumChildWhenFullyRefined = get_child_ids_and_num_children_when_fully_refined(parentElem);
  const auto & childElemIdsData = std::get<0>(childElemIdsDataAndnumChildWhenFullyRefined);
  const auto & numChildWhenFullyRefined = std::get<1>(childElemIdsDataAndnumChildWhenFullyRefined);
  return (has_child_with_invalid_id(childElemIdsData, numChildWhenFullyRefined));
}

unsigned Refinement::get_num_children(const stk::mesh::Entity elem) const
{
  if (is_parent(elem))
  {
    //const auto & [childElemIdsData, numChildWhenFullyRefined] = get_child_ids_and_num_children_when_fully_refined(elem);
    const auto & childElemIdsDataAndnumChildWhenFullyRefined = get_child_ids_and_num_children_when_fully_refined(elem);
    const auto & childElemIdsData = std::get<0>(childElemIdsDataAndnumChildWhenFullyRefined);
    const auto & numChildWhenFullyRefined = std::get<1>(childElemIdsDataAndnumChildWhenFullyRefined);
    return num_children_with_valid_id(childElemIdsData, numChildWhenFullyRefined);
  }
  return 0;
}

std::tuple<const uint64_t *,unsigned> Refinement::get_child_ids_and_num_children_when_fully_refined(const stk::mesh::Entity elem) const
{
  const unsigned numChildWhenFullyRefined = get_num_children_when_fully_refined(elem);
  const stk::mesh::Field<uint64_t> & childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
  const auto * childElemIdsData = stk::mesh::field_data(childElementIdsField, elem);
  STK_ThrowAssertMsg(nullptr != childElemIdsData, "Element is parent, but does not have " << childElementIdsField.name() << " defined.");
  return std::make_tuple(childElemIdsData, numChildWhenFullyRefined);
}

void Refinement::fill_child_element_ids(const stk::mesh::Entity elem, std::vector<stk::mesh::EntityId> & childElemIds) const
{
  childElemIds.clear();
  if (is_parent(elem))
  {
    //const auto & [childElemIdsData, numChildWhenFullyRefined] = get_child_ids_and_num_children_when_fully_refined(elem);
    const auto & childElemIdsDataAndnumChildWhenFullyRefined = get_child_ids_and_num_children_when_fully_refined(elem);
    const auto & childElemIdsData = std::get<0>(childElemIdsDataAndnumChildWhenFullyRefined);
    const auto & numChildWhenFullyRefined = std::get<1>(childElemIdsDataAndnumChildWhenFullyRefined);
    for (unsigned iChild=0; iChild<numChildWhenFullyRefined; ++iChild)
    {
      if ((stk::mesh::InvalidEntityId == childElemIdsData[iChild]))
        return;
      childElemIds.push_back(childElemIdsData[iChild]);
    }
  }
}

void Refinement::fill_children(const stk::mesh::Entity elem, std::vector<stk::mesh::Entity> & children) const
{
  children.clear();
  if (is_parent(elem))
  {
    //const auto & [childElemIdsData, numChildWhenFullyRefined] = get_child_ids_and_num_children_when_fully_refined(elem);
    const auto & childElemIdsDataAndnumChildWhenFullyRefined = get_child_ids_and_num_children_when_fully_refined(elem);
    const auto & childElemIdsData = std::get<0>(childElemIdsDataAndnumChildWhenFullyRefined);
    const auto & numChildWhenFullyRefined = std::get<1>(childElemIdsDataAndnumChildWhenFullyRefined);
    const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
    for (unsigned iChild=0; iChild<numChildWhenFullyRefined; ++iChild)
    {
      if ((stk::mesh::InvalidEntityId == childElemIdsData[iChild]))
        return;
      children.push_back(mesh.get_entity(stk::topology::ELEMENT_RANK, childElemIdsData[iChild]));
    }
  }
}

std::vector<stk::mesh::Entity> Refinement::get_children(const stk::mesh::Entity elem) const
{
   std::vector<stk::mesh::Entity> children;
   fill_children(elem, children);
   return children;
}

stk::mesh::Entity Refinement::get_parent(const stk::mesh::Entity elem) const
{
  stk::mesh::Entity parent;
  if (is_child(elem))
    parent = myMeta.mesh_bulk_data().get_entity(stk::topology::ELEMENT_RANK, get_parent_id(elem));
  return parent;
}

bool Refinement::is_refined_edge_node(const stk::mesh::Entity node) const
{
  return myMeta.mesh_bulk_data().bucket(node).member(refined_edge_node_part());
}

stk::mesh::Part & Refinement::child_part() const
{
  STK_ThrowAssert(myChildPart);
  return *myChildPart;
}

stk::mesh::Part & Refinement::parent_part() const
{
  STK_ThrowAssert(myParentPart);
  return *myParentPart;
}

stk::mesh::Part & Refinement::refined_edge_node_part() const
{
  STK_ThrowAssert(myRefinedEdgeNodePart);
  return *myRefinedEdgeNodePart;
}

stk::mesh::Part & Refinement::refined_quad_face_node_part() const
{
  STK_ThrowAssert(myRefinedQuadFaceNodePart);
  return *myRefinedQuadFaceNodePart;
}

stk::math::Vector3d Refinement::get_coordinates(const stk::mesh::Entity node, const int dim) const
{
  STK_ThrowAssertMsg(myCoordsField, "Coordinates field is not defined.");
  const double * coordsData = stk::mesh::field_data(*myCoordsField, node);
  STK_ThrowAssertMsg(nullptr != coordsData, "Node does not have " << myCoordsField->name() << " defined.");
  return stk::math::Vector3d(coordsData, dim);
}

template<size_t NUMREFINEDPARENTNODES, size_t NUMCHILDNODES>
static stk::mesh::Entity declare_child_element(stk::mesh::BulkData & mesh,
    EntityIdPool & entityIdPool,
    const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity elem,
    const std::array<stk::mesh::Entity,NUMREFINEDPARENTNODES> & parentNodes,
    const std::array<int,NUMCHILDNODES> & childNodeIndices)
{
  const stk::mesh::EntityId elemId = entityIdPool.get_EntityId(stk::topology::ELEMENT_RANK);
  stk::mesh::Entity childElem = mesh.declare_element(elemId, childParts);
  for (unsigned n=0; n<NUMCHILDNODES; ++n)
  {
    stk::mesh::Entity node = parentNodes[childNodeIndices[n]];
    mesh.declare_relation(childElem, node, n);
  }
  return childElem;
}

std::vector<std::pair<stk::mesh::Entity,stk::mesh::PartVector>> get_parent_sides_and_parts(const stk::mesh::BulkData & mesh, const stk::mesh::Entity parentElem)
{
  const unsigned numParentSides = mesh.bucket(parentElem).topology().num_sides();
  std::vector<std::pair<stk::mesh::Entity,stk::mesh::PartVector>> parentSidesAndParts;
  parentSidesAndParts.reserve(numParentSides);

  for (unsigned side=0; side<numParentSides; ++side)
  {
    const stk::mesh::Entity parentSide = find_entity_by_ordinal(mesh, parentElem, mesh.mesh_meta_data().side_rank(), side);
    if (mesh.is_valid(parentSide))
    {
      stk::mesh::PartVector childParts = get_removable_parts(mesh, mesh.bucket(parentSide));
      parentSidesAndParts.emplace_back(parentSide, childParts);
    }
    else
    {
      static stk::mesh::PartVector emptyParts;
      parentSidesAndParts.emplace_back(stk::mesh::Entity::InvalidEntity, emptyParts);
    }
  }

  return parentSidesAndParts;
}

template<class CHILDSIDECONTAINER>
static void attach_child_to_existing_sides_and_append_missing_sides_to_sides_to_create(stk::mesh::BulkData & mesh,
    const stk::mesh::Entity parentElem,
    const stk::mesh::Entity childElem,
    const CHILDSIDECONTAINER & childSideIndices,
    std::vector<SideDescription> & sideRequests)
{
  const stk::topology topology = mesh.bucket(childElem).topology();
  const stk::mesh::Entity * childElemNodes = mesh.begin_nodes(childElem);

  const auto parentSidesAndParts = get_parent_sides_and_parts(mesh, parentElem);

  for (unsigned s=0; s<childSideIndices.size(); ++s)
  {
    const stk::topology sideTopology = topology.side_topology(s);
    std::vector<stk::mesh::Entity> childSideNodes(sideTopology.num_nodes());
    topology.side_nodes(childElemNodes, s, childSideNodes.data());

    std::vector<stk::mesh::Entity> sides;
    stk::mesh::get_entities_through_relations(mesh, childSideNodes, mesh.mesh_meta_data().side_rank(), sides);

    if (!sides.empty())
    {
      STK_ThrowRequire(sides.size() == 1);

      attach_entity_to_element(mesh, mesh.mesh_meta_data().side_rank(), sides[0], childElem);
    }
    else if (childSideIndices[s] >= 0)
    {
      //const auto & [parentSide, parentSideParts] = parentSidesAndParts[childSideIndices[s]];
      const auto & parentSideAndparentSideParts = parentSidesAndParts[childSideIndices[s]];
      const auto & parentSide = parentSideAndparentSideParts.first;
      const auto & parentSideParts = parentSideAndparentSideParts.second;
      if (mesh.is_valid(parentSide))
      {
        sideRequests.emplace_back(childElem, s, parentSideParts);
      }
    }
  }
}

template<class CHILDDESCRIPTION, size_t NUMREFINEDPARENTNODES>
static void declare_child_elements_and_append_sides_to_create(stk::mesh::BulkData & mesh,
    EntityIdPool & entityIdPool,
    const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity parentElem,
    const std::array<stk::mesh::Entity,NUMREFINEDPARENTNODES> & parentNodes,
    const std::vector<CHILDDESCRIPTION> & childElemDescs,
    std::vector<stk::mesh::Entity> & childElements,
    std::vector<SideDescription> & sideRequests)
{
  childElements.clear();
  childElements.reserve(childElemDescs.size());
  for (auto && childElemDesc : childElemDescs)
  {
    stk::mesh::Entity childElement = declare_child_element(mesh, entityIdPool, childParts, parentElem, parentNodes, childElemDesc.nodeIds);
    childElements.push_back(childElement);
  }

  for(size_t iChild=0; iChild<childElemDescs.size(); ++iChild)
  {
    attach_child_to_existing_sides_and_append_missing_sides_to_sides_to_create(mesh, parentElem, childElements[iChild], childElemDescs[iChild].sideIds, sideRequests);
  }
}

stk::mesh::Field<uint64_t> & Refinement::get_child_element_ids_field(const unsigned numChildWhenFullyRefined) const
{
  switch(numChildWhenFullyRefined)
  {
  case 2:
      STK_ThrowAssert(myChildElementIds2Field);
      return *myChildElementIds2Field;
  case 4:
      STK_ThrowAssert(myChildElementIds4Field);
      return *myChildElementIds4Field;
  case 8:
      STK_ThrowAssert(myChildElementIds8Field);
      return *myChildElementIds8Field;
  default:
      ThrowRuntimeError("Unsupported number of children when fully refined in get_child_element_ids_field: " << numChildWhenFullyRefined);
      break;
  }
}

stk::mesh::EntityId * Refinement::get_child_element_ids(const unsigned numChildWhenFullyRefined, const stk::mesh::Entity parent) const
{
  stk::mesh::Field<uint64_t> & childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
  stk::mesh::EntityId * childElemIdsData = stk::mesh::field_data(childElementIdsField, parent);
  STK_ThrowAssertMsg(nullptr != childElemIdsData, "Element is does not have " << childElementIdsField.name() << " defined.");
  return childElemIdsData;
}

stk::mesh::EntityId * Refinement::get_child_element_ids(const unsigned numChildWhenFullyRefined, const stk::mesh::Bucket & bucket) const
{
  stk::mesh::Field<uint64_t> & childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
  stk::mesh::EntityId * childElemIdsData = stk::mesh::field_data(childElementIdsField, bucket);
  STK_ThrowAssertMsg(nullptr != childElemIdsData, "Bucket does not have " << childElementIdsField.name() << " defined.");
  return childElemIdsData;
}

template <size_t SIZE>
std::array<int,SIZE> get_rank_of_nodes_based_on_coordinates(const std::array<stk::math::Vector3d,SIZE> & nodeCoords)
{
  // initialize original index locations
  std::array<size_t,SIZE> index;
  std::iota(index.begin(), index.end(), 0);

  std::stable_sort(index.begin(), index.end(),
       [&nodeCoords](size_t i1, size_t i2) {return is_less_than_in_x_then_y_then_z(nodeCoords[i1], nodeCoords[i2]);});

  // invert indices to get node rank by coordinates
  std::array<int,SIZE> rank;
  for (size_t i=0; i < SIZE; ++i)
    rank[index[i]] = i;
  return rank;
}

void Refinement::set_parent_parts_and_parent_child_relation_fields(const stk::mesh::Entity parentElement, const std::vector<stk::mesh::Entity> & childElements, const unsigned numChildWhenFullyRefined)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  const bool isNewParent = !mesh.bucket(parentElement).member(*myParentPart);

  if (isNewParent)
  {
    stk::mesh::ConstPartVector removeParts;
    if (myActivePart)
      removeParts.push_back(myActivePart);
    mesh.change_entity_parts(parentElement, stk::mesh::ConstPartVector{myParentPart}, removeParts);
  }

  const stk::mesh::EntityId parentElementId = mesh.identifier(parentElement);
  const int parentRefinementLevel = refinement_level(parentElement);
  auto * childElemIds = get_child_element_ids(numChildWhenFullyRefined, parentElement);
  std::fill(childElemIds, childElemIds+numChildWhenFullyRefined, stk::mesh::InvalidEntityId);
  for (size_t iChild=0; iChild<childElements.size(); ++iChild)
  {
    childElemIds[iChild] = mesh.identifier(childElements[iChild]);
    set_parent_id(childElements[iChild], parentElementId);
    set_refinement_level(childElements[iChild], parentRefinementLevel+1);
  }

  if (isNewParent)
  {
    STK_ThrowAssert(mesh.bucket(parentElement).owned());
    set_originating_processor_for_parent_element(parentElement, mesh.parallel_rank());
  }
}

static void prolong_element_fields(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElems)
{
  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();
  for ( auto && stkField : allFields )
  {
    const FieldRef field(stkField);
    if( field.entity_rank() == stk::topology::ELEM_RANK && field.type_is<double>() )
    {

      const auto * parentElemData = field_data<double>(field, parentElem);
      if (nullptr != parentElemData)
      {
        for(const auto & childElem : childElems)
        {
          auto * childElemData = field_data<double>(field, childElem);

          if(!childElemData) continue;

          const unsigned fieldLength = field.length();
          for (unsigned i = 0; i < fieldLength; ++i)
          {
            childElemData[i] = parentElemData[i];
          }
        }
      }
    }
  }
}

static void restrict_element_fields(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & childElems)
{
  const stk::mesh::FieldVector & allFields = mesh.mesh_meta_data().get_fields();
  for ( auto && stkField : allFields )
  {
    const FieldRef field(stkField);
    if( field.entity_rank() == stk::topology::ELEM_RANK && field.type_is<double>() )
    {

      auto * parentElemData = field_data<double>(field, parentElem);
      if (nullptr != parentElemData)
      {
        const unsigned fieldLength = field.length();
        std::vector<double> averagedChildElemData(fieldLength,0.);
        double numActiveChildren = 0.0;
        for(const auto & childElem : childElems)
        {
          auto * childElemData = field_data<double>(field, childElem);

          if(!childElemData) continue;
          numActiveChildren += 1.;
          for (unsigned i = 0; i < fieldLength; ++i)
          {
            averagedChildElemData[i] += childElemData[i];
          }
        }
        if (numActiveChildren > 0.)
        {
          for (unsigned i = 0; i < fieldLength; ++i)
          {
            parentElemData[i] = averagedChildElemData[i]/numActiveChildren;
          }
        }
      }
    }
  }
}

void Refinement::refine_beam_2_and_append_sides_to_create(const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    const int caseId,
    std::vector<SideDescription> & sideRequests)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Entity * parentNodes = mesh.begin_nodes(parentElem);
  const std::array<stk::mesh::Entity,3> parentElemNodes{ parentNodes[0], parentNodes[1], elemChildEdgeNodes[0] };

  struct BeamDescription
  {
    std::array<int, 2> nodeIds;
    std::array<int, 0> sideIds;
  };

  // Do we need to handle sides of a beam?
  const std::vector<BeamDescription> newElems = { BeamDescription{ {{0,2}}, {{}} }, BeamDescription{ {{2,1}}, {{}} } };

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newElems, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 2);

  prolong_element_fields(mesh, parentElem, childElements);
}

void Refinement::refine_tri_3_and_append_sides_to_create(const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    const int caseId,
    std::vector<SideDescription> & sideRequests)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Entity * parentNodes = mesh.begin_nodes(parentElem);
  const std::array<stk::mesh::Entity,6> parentElemNodes{ parentNodes[0], parentNodes[1], parentNodes[2], elemChildEdgeNodes[0], elemChildEdgeNodes[1], elemChildEdgeNodes[2] };

  const std::array<stk::math::Vector3d,3> parentNodeCoords{{ get_coordinates(parentElemNodes[0], 2), get_coordinates(parentElemNodes[1], 2), get_coordinates(parentElemNodes[2], 2) }};

  const std::array<int,3> parentNodeRank = get_rank_of_nodes_based_on_coordinates(parentNodeCoords);

  const std::vector<TriRefiner::TriDescription> newTris = TriRefiner::refinement_child_nodes_and_sides_tri3(caseId, parentNodeCoords, parentNodeRank);

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newTris, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 4);

  prolong_element_fields(mesh, parentElem, childElements);
}

void Refinement::refine_quad_4_and_append_sides_to_create(const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    const int caseId,
    std::vector<SideDescription> & sideRequests)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Entity * parentNodes = mesh.begin_nodes(parentElem);
  const stk::mesh::Entity elemCentroidNode = myNodeRefiner.get_element_centroid_child_node(parentElem);
  const std::array<stk::mesh::Entity,9> parentElemNodes{ parentNodes[0], parentNodes[1], parentNodes[2], parentNodes[3], elemChildEdgeNodes[0], elemChildEdgeNodes[1], elemChildEdgeNodes[2], elemChildEdgeNodes[3], elemCentroidNode };

  const std::vector<QuadRefiner::QuadDescription> newElems = QuadRefiner::refinement_child_nodes_and_sides_quad4(caseId);

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newElems, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 4);

  prolong_element_fields(mesh, parentElem, childElements);
}

void Refinement::refine_tet_4_and_append_sides_to_create(const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    const int caseId,
    std::vector<SideDescription> & sideRequests)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Entity * parentNodes = mesh.begin_nodes(parentElem);
  const std::array<stk::mesh::Entity,10> parentElemNodes{ parentNodes[0], parentNodes[1], parentNodes[2], parentNodes[3], elemChildEdgeNodes[0], elemChildEdgeNodes[1], elemChildEdgeNodes[2], elemChildEdgeNodes[3], elemChildEdgeNodes[4], elemChildEdgeNodes[5] };

  const std::array<stk::math::Vector3d,4> parentNodeCoords{{
    get_coordinates(parentNodes[0]),
    get_coordinates(parentNodes[1]),
    get_coordinates(parentNodes[2]),
    get_coordinates(parentNodes[3]) }};

  const std::array<int,4> parentNodeRank= get_rank_of_nodes_based_on_coordinates(parentNodeCoords);

  const double needSides = mesh.num_sides(parentElem) > 0;
  const std::vector<moab::SimplexTemplateRefiner::TetDescription> newTets = moab::SimplexTemplateRefiner::refinement_child_nodes_and_sides_tet4(caseId, parentNodeCoords, parentNodeRank, needSides);

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newTets, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 8);
  //prolong element fields
  prolong_element_fields(mesh, parentElem, childElements);
}

void Refinement::refine_hex_8_and_append_sides_to_create(const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity parentElem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    const int caseId,
    std::vector<SideDescription> & sideRequests)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Entity * parentNodes = mesh.begin_nodes(parentElem);
  const stk::mesh::Entity elemCentroidNode = myNodeRefiner.get_element_centroid_child_node(parentElem);
  const std::vector<stk::mesh::Entity> elemChildFaceNodes = myNodeRefiner.get_element_child_face_nodes(mesh, parentElem);
  // This is not exactly obvious, but comes out of topology_data.hpp
  const std::array<stk::mesh::Entity,27> parentElemNodes{
    parentNodes[0], parentNodes[1], parentNodes[2], parentNodes[3], parentNodes[4], parentNodes[5], parentNodes[6], parentNodes[7],
    elemChildEdgeNodes[0], elemChildEdgeNodes[1], elemChildEdgeNodes[2], elemChildEdgeNodes[3],
    elemChildEdgeNodes[8], elemChildEdgeNodes[9], elemChildEdgeNodes[10], elemChildEdgeNodes[11],
    elemChildEdgeNodes[4], elemChildEdgeNodes[5], elemChildEdgeNodes[6], elemChildEdgeNodes[7],
    elemCentroidNode,
    elemChildFaceNodes[4], elemChildFaceNodes[5], elemChildFaceNodes[3], elemChildFaceNodes[1], elemChildFaceNodes[0], elemChildFaceNodes[2]};

  const std::vector<HexRefiner::HexDescription> newElems = HexRefiner::refinement_child_nodes_and_sides_hex8(caseId);

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newElems, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 8);

  prolong_element_fields(mesh, parentElem, childElements);
}

static unsigned num_new_child_elements_for_case_id(const stk::topology & elemTopology, const int caseId)
{
  switch(elemTopology())
    {
        case stk::topology::BEAM_2:
        {
          STK_ThrowRequire(caseId == 1);
          return 2;
        }
    case stk::topology::TRI_3:
    case stk::topology::TRI_3_2D:
        return TriRefiner::num_new_child_elements_tri3(caseId);
    case stk::topology::QUAD_4:
    case stk::topology::QUAD_4_2D:
        return QuadRefiner::num_new_child_elements_quad4(caseId);
    case stk::topology::TETRAHEDRON_4:
        return moab::SimplexTemplateRefiner::num_new_child_elements_tet4(caseId);
    case stk::topology::HEXAHEDRON_8:
        return HexRefiner::num_new_child_elements_hex8(caseId);
    default:
        ThrowRuntimeError("Element topology not found in refine_element: " << elemTopology.name());
    }
}

void Refinement::adapt_element_and_append_sides_to_create(const stk::topology & elemTopology,
    const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity elem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    const int preCaseId,
    const int postCaseId,
    std::vector<SideDescription> & sideRequests,
    std::vector<stk::mesh::Entity> & elementsToDelete,
    std::vector<stk::mesh::Entity> & elementsThatAreNoLongerParents)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  if (0 != preCaseId)
  {
    const std::vector<stk::mesh::Entity> existingChildrenToDelete = get_children(elem);
    restrict_element_fields(mesh, elem, existingChildrenToDelete);
    elementsToDelete.insert(elementsToDelete.end(), existingChildrenToDelete.begin(), existingChildrenToDelete.end());
  }

  if (0 == postCaseId)
  {
    elementsThatAreNoLongerParents.push_back(elem);
    return;
  }

  switch(elemTopology())
    {
    case stk::topology::BEAM_2:
        refine_beam_2_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, postCaseId, sideRequests);
        break;
    case stk::topology::TRI_3:
    case stk::topology::TRI_3_2D:
        refine_tri_3_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, postCaseId, sideRequests);
        break;
    case stk::topology::QUAD_4:
    case stk::topology::QUAD_4_2D:
        refine_quad_4_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, postCaseId, sideRequests);
        break;
    case stk::topology::TETRAHEDRON_4:
        refine_tet_4_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, postCaseId, sideRequests);
        break;
    case stk::topology::HEXAHEDRON_8:
        refine_hex_8_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, postCaseId, sideRequests);
        break;
    default:
        ThrowRuntimeError("Element topology not found in refine_element: " << elemTopology.name());
        break;
    }


}

size_t Refinement::count_new_child_elements(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData, const bool doingRefinement) const
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();

  std::vector<stk::mesh::Entity> elemChildEdgeNodes;
  ElementEdgeCaseIds elementEdgeCaseIds;

  size_t numNewChildElems = 0;
  for(const auto & bucketData : bucketsData)
  {
    for(const auto & elem : std::get<2>(bucketData))
    {
      const stk::topology bucketTopology = std::get<0>(bucketData);
      edgeMarker.fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(myNodeRefiner, elem, bucketTopology, doingRefinement, elementEdgeCaseIds, elemChildEdgeNodes);
      if (elementEdgeCaseIds.has_changed())
        numNewChildElems += num_new_child_elements_for_case_id(bucketTopology, elementEdgeCaseIds.post_adapt_case_id());
    }
  }

  return numNewChildElems;
}

stk::mesh::PartVector Refinement::get_parts_for_child_elements(const stk::mesh::Bucket & parentBucket) const
{
  const stk::mesh::PartVector& parentParts = parentBucket.supersets();

  stk::mesh::PartVector childParts;
  childParts.reserve(parentParts.size());

  for ( auto&& part : parentParts )
  {
    stk::mesh::EntityRank part_rank = part->primary_entity_rank();
    if ((part_rank == stk::topology::INVALID_RANK || part_rank == parentBucket.entity_rank()) &&
        part != myParentPart &&
        (!stk::mesh::is_auto_declared_part(*part) || stk::mesh::is_topology_root_part(*part)))
    {
      childParts.push_back(part);
    }
  }

  childParts.push_back(myChildPart);
  if (myActivePart)
    childParts.push_back(myActivePart);

  return childParts;
}

std::vector<Refinement::BucketData> Refinement::get_buckets_data_for_candidate_elements_to_adapt(const EdgeMarkerInterface & edgeMarker, const bool doingRefinement) const
{
  // Cache off bucket data to avoid looping over buckets while modifying elements
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();

  stk::mesh::EntityVector bucketCandidateElements;
  std::vector<std::tuple<stk::topology,stk::mesh::PartVector,stk::mesh::EntityVector>> bucketsData;

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    bucketCandidateElements.clear();

    for (auto && elem : *bucketPtr)
      if (edgeMarker.is_element_a_candidate_for_adaptation(elem, doingRefinement))
        bucketCandidateElements.push_back(elem);

    if (!bucketCandidateElements.empty())
      bucketsData.emplace_back(bucketPtr->topology(), get_parts_for_child_elements(*bucketPtr), bucketCandidateElements);
  }
  return bucketsData;
}

static int case_id_for_fully_refined(const stk::topology elementTopology)
{
  return (1<<elementTopology.num_edges())-1;
}

static bool element_going_from_partially_refined_to_fully_refined(const stk::topology elementTopology, const bool doingRefinement, const ElementEdgeCaseIds & elementEdgeCaseIds)
{
  return doingRefinement && elementEdgeCaseIds.pre_adapt_case_id() > 0 && elementEdgeCaseIds.post_adapt_case_id() == case_id_for_fully_refined(elementTopology);
}

void Refinement::adapt_elements_and_store_sides_to_create(const EdgeMarkerInterface & edgeMarker,
    const std::vector<BucketData> & bucketsData,
    const bool doingRefinement,
    std::vector<SideDescription> & sideRequests,
    std::vector<stk::mesh::Entity> & elementsToDelete,
    std::vector<stk::mesh::Entity> & elementsThatAreNoLongerParents,
    std::vector<BucketData> & bucketDataForNewChildElementsThatMightNeedToBeRefined)
{
  const size_t numNewElements = count_new_child_elements(edgeMarker, bucketsData, doingRefinement);
  myEntityIdPool.reserve(stk::topology::ELEMENT_RANK, numNewElements, myAssert32Bit, myForce64Bit);

  bucketDataForNewChildElementsThatMightNeedToBeRefined.clear();
  std::vector<stk::mesh::Entity> bucketChildElementsThatMightNeedToBeRefined;
  std::vector<stk::mesh::Entity> elemChildEdgeNodes;
  ElementEdgeCaseIds elementEdgeCaseIds;
  for(const auto & [bucketTopology, bucketChildParts, bucketElements] : bucketsData)
  {
    bucketChildElementsThatMightNeedToBeRefined.clear();
    for(const auto & elem : bucketElements)
    {
      edgeMarker.fill_adaptation_caseIds_and_refined_edge_nodes_if_changed(myNodeRefiner, elem, bucketTopology, doingRefinement, elementEdgeCaseIds, elemChildEdgeNodes);
      if (elementEdgeCaseIds.has_changed())
      {
        adapt_element_and_append_sides_to_create(bucketTopology, bucketChildParts, elem, elemChildEdgeNodes, elementEdgeCaseIds.pre_adapt_case_id(), elementEdgeCaseIds.post_adapt_case_id(), sideRequests, elementsToDelete, elementsThatAreNoLongerParents);

        if (element_going_from_partially_refined_to_fully_refined(bucketTopology, doingRefinement, elementEdgeCaseIds)) // A second level of refinement possibly needed to remove hanging nodes
        {
          auto childElements = get_children(elem);
          bucketChildElementsThatMightNeedToBeRefined.insert(bucketChildElementsThatMightNeedToBeRefined.end(), childElements.begin(), childElements.end());
        }
      }
    }
    if (!bucketChildElementsThatMightNeedToBeRefined.empty())
      bucketDataForNewChildElementsThatMightNeedToBeRefined.emplace_back(bucketTopology, bucketChildParts, bucketChildElementsThatMightNeedToBeRefined);
  }
}

stk::mesh::PartVector Refinement::get_parts_for_new_refined_edge_nodes() const
{
  stk::mesh::PartVector refinedEdgeNodeParts = { &myMeta.get_topology_root_part(stk::topology::NODE), &refined_edge_node_part() };
  if (myActivePart)
    refinedEdgeNodeParts.push_back(myActivePart);
  return refinedEdgeNodeParts;
}

stk::mesh::PartVector Refinement::get_parts_for_new_refined_element_centroid_nodes() const
{
  stk::mesh::PartVector refinedElemCentroidNodeParts = { &myMeta.get_topology_root_part(stk::topology::NODE) };
  if (myActivePart)
    refinedElemCentroidNodeParts.push_back(myActivePart);
  return refinedElemCentroidNodeParts;
}

stk::mesh::PartVector Refinement::get_parts_for_new_refined_quad_face_nodes() const
{
  stk::mesh::PartVector refinedQuadFaceNodeParts = { &myMeta.get_topology_root_part(stk::topology::NODE), &refined_quad_face_node_part() };
  if (myActivePart)
    refinedQuadFaceNodeParts.push_back(myActivePart);
  return refinedQuadFaceNodeParts;
}

bool Refinement::locally_have_any_hanging_refined_nodes() const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part() & !parent_part();

  std::array<stk::mesh::Entity,3> edgeNodes;
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const stk::topology elemTopology = bucketPtr->topology();
    const unsigned numEdges = elemTopology.num_edges();

    for (auto && elem : *bucketPtr)
    {
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);

      for (unsigned iEdge = 0; iEdge < numEdges; ++iEdge)
      {
        elemTopology.edge_nodes(elemNodes, iEdge, edgeNodes.data());
        const Edge edge = edge_from_edge_nodes(mesh, edgeNodes[0], edgeNodes[1]);
        if (myNodeRefiner.is_edge_marked_for_refinement(edge))
        {
          krinolog << "Found element with hanging node on edge " << debug_edge(mesh, edge) << " of element " << debug_entity_1line(mesh, elem) << stk::diag::dendl;
          return true;
        }
      }
    }
  }
  return false;
}

bool Refinement::have_any_hanging_refined_nodes() const
{
  return stk::is_true_on_any_proc(myMeta.mesh_bulk_data().parallel(), locally_have_any_hanging_refined_nodes());
}

std::string Refinement::locally_check_leaf_children_have_parents_on_same_proc() const
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  std::ostringstream localErrorMsg;
  auto buckets = mesh.get_buckets(stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part());
  for(auto bucket : buckets)
  {
    for(auto && elem : *bucket)
    {
      if(is_child(elem) && !is_parent(elem))
      {
        auto parent = get_parent(elem);
        if(!mesh.is_valid(parent) || !mesh.bucket(parent).owned())
        {
          localErrorMsg << debug_entity_1line(mesh, elem) << "\n";
          return localErrorMsg.str();
        }
      }
    }
  }
  return localErrorMsg.str();
}

void Refinement::check_leaf_children_have_parents_on_same_proc() const
{
  RequireEmptyErrorMsg(myMeta.mesh_bulk_data().parallel(), locally_check_leaf_children_have_parents_on_same_proc(), "Leaf child without parent owned on same proc.");
}

unsigned Refinement::rebalance_element_count_incorporating_parallel_owner_constraints(const stk::mesh::Entity elem) const
{
  if(is_parent(elem))
  {
    std::vector<stk::mesh::Entity> elemDependents;
    fill_child_elements_that_must_stay_on_same_proc_as_parent(elem, elemDependents);
    return elemDependents.size(); // child cost
  }

  if(is_child(elem)) // if not a parent but is a child, must be leaf element -> cost included with parent
    return 0;
  return 1; // if not a parent or child, must be completed unadapted element. -> self cost
}


bool Refinement::has_parallel_owner_rebalance_constraint(const stk::mesh::Entity entity) const
{
  if(is_parent(entity))
  {
    //if a parent, check if any children are parents or invalid elements (already moved)
    //if so, constrained. If not, not constrained
    const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
    std::vector<stk::mesh::Entity> children;
    fill_children(entity, children);
    for(auto && child : children)
      if (!mesh.is_valid(child) || is_parent(child))
        return true;
    return false;
  }

  if(is_child(entity)) //if not a parent but is a child, must be leaf element, constrained
    return true;
  return false; //if not a parent or child, must be completed unadapted element. No constraint
}

void Refinement::fill_child_elements_that_must_stay_on_same_proc_as_parent(const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & dependents) const
{
  // For non-parents, this will correctly just produce an empty dependents vector
  fill_children(parent, dependents);
  if (!dependents.empty())
  {
    const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
    size_t iKeep = 0;
    for(size_t i=0; i<dependents.size(); ++i)
    {
      const stk::mesh::Entity & child = dependents[i];
      if (mesh.is_valid(child) && !is_parent(child))
        dependents[iKeep++] = child;
    }
    dependents.resize(iKeep);
  }
}

static void adjust_parent_and_child_rebalance_weights(stk::mesh::Field<double> & elemWtField, const stk::mesh::Entity parent, std::vector<stk::mesh::Entity> & children)
{
  if (children.empty())
    return;

  double childWtSum = 0.;
  for (auto && child : children)
  {
    double & childWt = *stk::mesh::field_data(elemWtField, child);
    childWtSum += childWt;
    childWt = 0.;
  }
  double & parentWt = *stk::mesh::field_data(elemWtField, parent);
  parentWt += childWtSum;
}

void Refinement::update_element_rebalance_weights_incorporating_parallel_owner_constraints(stk::mesh::Field<double> & elemWtField) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  std::vector<stk::mesh::Entity> elemDependents;
  for (auto && bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, stk::mesh::selectField(elemWtField) & myMeta.locally_owned_part()))
  {
    for (auto && elem : *bucketPtr)
    {
      fill_child_elements_that_must_stay_on_same_proc_as_parent(elem, elemDependents);
      adjust_parent_and_child_rebalance_weights(elemWtField, elem, elemDependents);
    }
  }
}

void Refinement::adapt_elements_and_sides(const EdgeMarkerInterface & edgeMarker, const bool doingRefinement)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  const std::vector<BucketData> bucketsData = get_buckets_data_for_candidate_elements_to_adapt(edgeMarker, doingRefinement);
  std::vector<BucketData> bucketDataForNewChildElementsThatMightNeedToBeRefined;

  if(stk::is_true_on_any_proc(mesh.parallel(), !bucketsData.empty()))
  {
    std::vector<SideDescription> sideRequests;

    mesh.modification_begin();

    destroy_custom_ghostings();

    if (doingRefinement)
    {
      myNodeRefiner.create_refined_element_centroid_nodes(mesh, get_parts_for_new_refined_element_centroid_nodes());

      myNodeRefiner.create_refined_quad_face_nodes(mesh, get_parts_for_new_refined_quad_face_nodes(), myRefinedQuadFaceNodeParentIdsField);

      myNodeRefiner.create_refined_edge_nodes(mesh, get_parts_for_new_refined_edge_nodes(), myRefinedEdgeNodeParentIdsField);
    }

    std::vector<stk::mesh::Entity> elementsToDelete;
    std::vector<stk::mesh::Entity> elementsThatAreNoLongerParents;

    adapt_elements_and_store_sides_to_create(edgeMarker, bucketsData, doingRefinement, sideRequests, elementsToDelete, elementsThatAreNoLongerParents, bucketDataForNewChildElementsThatMightNeedToBeRefined);
    stk::mesh::destroy_elements_no_mod_cycle(mesh, elementsToDelete, mesh.mesh_meta_data().universal_part());
    remove_parent_parts(elementsThatAreNoLongerParents);

    mesh.modification_end();

    if(stk::is_true_on_any_proc(mesh.parallel(), !sideRequests.empty()))
    {
      batch_create_sides(mesh, sideRequests);
    }
  }

  if(stk::is_true_on_any_proc(mesh.parallel(), !bucketDataForNewChildElementsThatMightNeedToBeRefined.empty()))
  {
    std::vector<SideDescription> sideRequests;

    mesh.modification_begin();

    std::vector<stk::mesh::Entity> shouldBeEmpty_ElementsToDelete;
    std::vector<stk::mesh::Entity> shouldBeEmpty_ElementsThatAreNoLongerParents;
    std::vector<BucketData> shouldBeEmpty_bucketDataForNextRound;
    adapt_elements_and_store_sides_to_create(edgeMarker, bucketDataForNewChildElementsThatMightNeedToBeRefined, doingRefinement, sideRequests, shouldBeEmpty_ElementsToDelete, shouldBeEmpty_ElementsThatAreNoLongerParents, shouldBeEmpty_bucketDataForNextRound);

    const bool logicError = !shouldBeEmpty_ElementsToDelete.empty() || !shouldBeEmpty_ElementsThatAreNoLongerParents.empty() || !shouldBeEmpty_bucketDataForNextRound.empty();

    mesh.modification_end();

    STK_ThrowRequireMsg(stk::is_true_on_all_procs(mesh.parallel(), !logicError), "Unexpected error in adapt_elements_and_sides");

    if(stk::is_true_on_any_proc(mesh.parallel(), !sideRequests.empty()))
    {
      batch_create_sides(mesh, sideRequests);
    }
  }
}

void Refinement::remove_parent_parts(const std::vector<stk::mesh::Entity> & elements)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  stk::mesh::ConstPartVector removeParts{myParentPart};
  stk::mesh::ConstPartVector addParts;
  if (myActivePart)
    addParts.push_back(myActivePart);

  for (auto && element : elements)
    mesh.change_entity_parts(element, addParts, removeParts);
}

void Refinement::respect_originating_proc_for_parents_modified_by_unrefinement(const std::vector<stk::mesh::Entity> & parentsModifiedByUnrefinement, const std::vector<int> & originatingProcForParentsModifiedByUnrefinement)
{
  STK_ThrowAssert(parentsModifiedByUnrefinement.size() == originatingProcForParentsModifiedByUnrefinement.size());
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  if (1 == mesh.parallel_size())
    return;

  std::vector<stk::mesh::EntityProc> entitiesToMove;
  for (size_t i=0; i<parentsModifiedByUnrefinement.size(); ++i)
  {
    if (is_parent(parentsModifiedByUnrefinement[i]))
      set_originating_processor_for_parent_element(parentsModifiedByUnrefinement[i], originatingProcForParentsModifiedByUnrefinement[i]);
    else
      entitiesToMove.emplace_back(parentsModifiedByUnrefinement[i], originatingProcForParentsModifiedByUnrefinement[i]);
  }

  if(stk::is_true_on_any_proc(mesh.parallel(), !entitiesToMove.empty()))
  {
    mesh.change_entity_owner(entitiesToMove);
    fix_face_and_edge_ownership(mesh);  // If this becomes a performance hotspot, then we should consider merging the edge/face movement with the element movement
  }
}

std::vector<int> Refinement::get_originating_procs_for_elements(const std::vector<stk::mesh::Entity> & elements) const
{
  std::vector<int> originatingProcsForElems;
  if (1 == myMeta.mesh_bulk_data().parallel_size())
  {
    originatingProcsForElems.assign(elements.size(), 0);
  }
  else
  {
    originatingProcsForElems.reserve(elements.size());
    for (auto elem : elements)
      originatingProcsForElems.push_back(get_originating_processor_for_parent_element(elem));
  }
  return originatingProcsForElems;
}

void Refinement::destroy_custom_ghostings()
{
  // stk::mesh::destroy_elements_no_mod_cycle() currently can create parallel inconsistencies in
  // the mesh when there are custom ghostings (e.g. from contact or transfers) present. See
  // COMPSIMHD-15326 for example. Since those are going to have to be rebuilt after the adaptivity
  // anyway we clear them out here as a workaround.
  krino::destroy_custom_ghostings(myMeta.mesh_bulk_data());
}

void Refinement::finalize()
{
  stk::diag::TimeBlock timer_(myFixPartsandOwnersTimer);
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  activate_selected_entities_touching_active_elements(
      mesh, myMeta.side_rank(), myMeta.universal_part(), *myActivePart);
  fix_node_owners_to_assure_active_owned_element_for_node(mesh, *myActivePart);
}

bool Refinement::refine_elements(const EdgeMarkerInterface & edgeMarker)
{
  stk::diag::TimeBlock timer_(refineTimer.rootTimer);

  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  {
    stk::diag::TimeBlock timer_2(refineTimer.checkLeafChildren);
    check_leaf_children_have_parents_on_same_proc();
  }

  {
    stk::diag::TimeBlock timer_2(refineTimer.findEdgesToRefine);
    find_edges_to_refine(edgeMarker);
  }
  bool didMakeAnyChanges = false;

  if (stk::is_true_on_any_proc(mesh.parallel(), locally_have_edges_to_refine()))
  {
    didMakeAnyChanges = true;

    {
      stk::diag::TimeBlock timer_2(refineTimer.doRefinement);
      adapt_elements_and_sides(edgeMarker, true);
    }

    {
      stk::diag::TimeBlock timer_2(refineTimer.prolongNodes);
      myNodeRefiner.prolong_refined_nodes(mesh);
    }
  }

  return didMakeAnyChanges;
}

bool Refinement::unrefine_elements(const EdgeMarkerInterface & edgeMarker)
{
  stk::diag::TimeBlock timer_(unrefineTimer.rootTimer);

  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  {
    stk::diag::TimeBlock timer_2(unrefineTimer.checkLeafChildren);
    check_leaf_children_have_parents_on_same_proc();
  }

  bool didMakeAnyChanges = false;
  if(stk::is_true_on_any_proc(mesh.parallel(), edgeMarker.locally_have_elements_to_unrefine()))
  {
    std::vector<stk::mesh::Entity> ownedParentElementsModifiedByUnrefinement;
    std::vector<int> originatingProcForParentsBeingModified;
    {
      stk::diag::TimeBlock timer_2(unrefineTimer.findEdgesToUnrefine);
      edgeMarker.mark_edges_to_be_unrefined(myNodeRefiner);
      ownedParentElementsModifiedByUnrefinement = edgeMarker.get_parent_elements_that_will_be_modified_by_unrefinement(myNodeRefiner);
      originatingProcForParentsBeingModified = get_originating_procs_for_elements(ownedParentElementsModifiedByUnrefinement);
    }

    if(stk::is_true_on_any_proc(mesh.parallel(), !ownedParentElementsModifiedByUnrefinement.empty()))
    {
      didMakeAnyChanges = true;

      {
        stk::diag::TimeBlock timer_2(unrefineTimer.doUnrefinement);
        adapt_elements_and_sides(edgeMarker, false);
      }

      {
        stk::diag::TimeBlock timer_2(unrefineTimer.fixFaceEdgeOwnership);

        respect_originating_proc_for_parents_modified_by_unrefinement(
            ownedParentElementsModifiedByUnrefinement, originatingProcForParentsBeingModified);
      }
    }
  }

  {
    stk::diag::TimeBlock timer_2(unrefineTimer.checkLeafChildren);
    check_leaf_children_have_parents_on_same_proc();
  }

  return didMakeAnyChanges;
}

bool Refinement::do_refinement(const EdgeMarkerInterface & edgeMarker)
{
  bool didMakeAnyChanges = false;

  didMakeAnyChanges |= refine_elements(edgeMarker);
  didMakeAnyChanges |= unrefine_elements(edgeMarker);

  if (didMakeAnyChanges && myActivePart)
  {
    finalize();
  }

  STK_ThrowAssertMsg(!have_any_hanging_refined_nodes(), "Mesh has hanging refined node.");

  check_leaf_children_have_parents_on_same_proc();

  return didMakeAnyChanges;
}

bool Refinement::do_uniform_refinement(const int numUniformRefinementLevels)
{
  UniformEdgeMarker uniformMarker(myMeta.mesh_bulk_data(), *this);
  bool didMakeAnyChanges = false;
  for (int i=0; i<numUniformRefinementLevels; ++i)
    didMakeAnyChanges |= do_refinement(uniformMarker);
  return didMakeAnyChanges;
}

void Refinement::fully_unrefine_mesh()
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  std::vector<stk::mesh::Entity> allChildElems;
  stk::mesh::get_selected_entities( child_part(), mesh.buckets(stk::topology::ELEMENT_RANK), allChildElems );

  const stk::mesh::Selector locallyOwnedRootParents = mesh.mesh_meta_data().locally_owned_part() & parent_part() & !child_part();
  std::vector<stk::mesh::Entity> ownedParentElems;
  stk::mesh::get_selected_entities( locallyOwnedRootParents, mesh.buckets(stk::topology::ELEMENT_RANK), ownedParentElems );

  stk::mesh::ConstPartVector addParts;
  if (myActivePart)
    addParts.push_back(myActivePart);

  mesh.modification_begin();
  destroy_custom_ghostings();
  stk::mesh::destroy_elements_no_mod_cycle(mesh, allChildElems, mesh.mesh_meta_data().universal_part());
  mesh.change_entity_parts(ownedParentElems, addParts, stk::mesh::ConstPartVector{myParentPart});
  mesh.modification_end();
}

void Refinement::delete_parent_elements()
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  std::vector<stk::mesh::Entity> allParentElems;
  stk::mesh::get_selected_entities( parent_part(), mesh.buckets(stk::topology::ELEMENT_RANK), allParentElems );

  std::vector<stk::mesh::Entity> ownedChildElems;
  stk::mesh::get_selected_entities(  mesh.mesh_meta_data().locally_owned_part() & child_part() & !parent_part(), mesh.buckets(stk::topology::ELEMENT_RANK), ownedChildElems );

  mesh.modification_begin();
  destroy_custom_ghostings();
  stk::mesh::destroy_elements_no_mod_cycle(mesh, allParentElems, mesh.mesh_meta_data().universal_part());
  stk::mesh::ConstPartVector addParts;
  mesh.change_entity_parts(ownedChildElems, addParts, stk::mesh::ConstPartVector{myChildPart});
  mesh.modification_end();
}

void Refinement::find_edges_to_refine(const EdgeMarkerInterface & edgeMarker)
{
  edgeMarker.mark_edges_to_be_refined(myNodeRefiner);
}

stk::mesh::EntityId Refinement::get_parent_id(const stk::mesh::Entity elem) const
{
  STK_ThrowAssertMsg(myParentElementIdField, "Parent ids field is not defined.");
  const auto * parentElemIdData = stk::mesh::field_data(*myParentElementIdField, elem);
  STK_ThrowAssertMsg(nullptr != parentElemIdData, "Element is does not have " << myParentElementIdField->name() << " defined.");
  return *parentElemIdData;
}

std::pair<stk::mesh::EntityId,int> Refinement::get_parent_id_and_parallel_owner_rank(const stk::mesh::Entity child) const
{
  const int parentParallelOwnerRank = is_parent(child) ? get_originating_processor_for_parent_element(child) : myMeta.mesh_bulk_data().parallel_owner_rank(child);
  return {get_parent_id(child), parentParallelOwnerRank};
}

void Refinement::set_parent_id(const stk::mesh::Entity elem, const stk::mesh::EntityId parentElemId) const
{
  STK_ThrowAssertMsg(myParentElementIdField, "Parent ids field is not defined.");
  auto * parentElemIdData = stk::mesh::field_data(*myParentElementIdField, elem);
  STK_ThrowAssertMsg(nullptr != parentElemIdData, "Element is does not have " << myParentElementIdField->name() << " defined.");
  *parentElemIdData = parentElemId;
}

namespace {

struct OriginatingProcParentIdChildId {
  OriginatingProcParentIdChildId(const int inOriginatingProc, const stk::mesh::EntityId inParentId, const stk::mesh::EntityId inChildId) : originatingProc(inOriginatingProc), parentId(inParentId), childId(inChildId)  {}

  int originatingProc;
  stk::mesh::EntityId parentId;
  stk::mesh::EntityId childId;
};

}

static
void pack_parent_and_child_ids_for_originating_proc(const stk::mesh::BulkData & mesh,
    const std::vector<OriginatingProcParentIdChildId> & originatingProcsParentIdChildId,
    stk::CommSparse &commSparse)
{
  std::vector<int> elemCommProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto & originatingProcParentIdChildId : originatingProcsParentIdChildId)
    {
      commSparse.send_buffer(originatingProcParentIdChildId.originatingProc).pack(originatingProcParentIdChildId.parentId);
      commSparse.send_buffer(originatingProcParentIdChildId.originatingProc).pack(originatingProcParentIdChildId.childId);
    }
  });
}

static
void unpack_parent_and_child_id_on_originating_proc(const stk::mesh::BulkData & mesh,
    std::vector<std::pair<stk::mesh::Entity, stk::mesh::EntityId>> & parentsAndOffProcChildId,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId parentId;
      commSparse.recv_buffer(procId).unpack(parentId);
      stk::mesh::EntityId childId;
      commSparse.recv_buffer(procId).unpack(childId);
      stk::mesh::Entity parent = mesh.get_entity(stk::topology::ELEMENT_RANK, parentId);
      STK_ThrowRequire(mesh.is_valid(parent));
      parentsAndOffProcChildId.emplace_back(parent, childId);
    }
  });
}

static void communicate_to_get_parents_with_off_proc_child(const stk::mesh::BulkData & mesh, const std::vector<OriginatingProcParentIdChildId> & originatingProcsParentIdChildId, std::vector<std::pair<stk::mesh::Entity, stk::mesh::EntityId>> & parentsAndOffProcChildId)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_parent_and_child_ids_for_originating_proc(mesh, originatingProcsParentIdChildId, commSparse);
  unpack_parent_and_child_id_on_originating_proc(mesh, parentsAndOffProcChildId, commSparse);
}

void Refinement::fill_parents_and_children_and_parents_with_off_proc_child(std::vector<stk::mesh::Entity> & parents, std::vector<stk::mesh::Entity> & children, std::vector<ParentAndChildId> & parentsAndOffProcChildId) const
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  parents.clear();
  children.clear();
  parentsAndOffProcChildId.clear();
  std::vector<OriginatingProcParentIdChildId> originatingProcsParentIdChildId;

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, myMeta.locally_owned_part()))
  {
    for (auto && elem : *bucket)
    {
      const stk::mesh::EntityId parentId = get_parent_id(elem);
      if (stk::mesh::InvalidEntityId != parentId)
      {
        const stk::mesh::Entity parent = myMeta.mesh_bulk_data().get_entity(stk::topology::ELEMENT_RANK, parentId);
        if (mesh.is_valid(parent) && mesh.bucket(parent).owned())
        {
          children.push_back(elem);
          parents.push_back(parent);
        }
        else
        {
          STK_ThrowAssertMsg(is_parent(elem), "In fill_parents_and_children_and_parents_with_off_proc_child(), found non-parent element " << mesh.identifier(elem) << " without local parent.");
          children.push_back(elem);
          originatingProcsParentIdChildId.emplace_back(get_originating_processor_for_parent_element(elem), parentId, mesh.identifier(elem));
        }
      }
    }
  }

  communicate_to_get_parents_with_off_proc_child(mesh, originatingProcsParentIdChildId, parentsAndOffProcChildId);

  for (auto && parentAndChildId : parentsAndOffProcChildId)
    parents.push_back(parentAndChildId.first);

  stk::util::sort_and_unique(parents);
}

void Refinement::restore_parent_and_child_element_parts(const std::vector<stk::mesh::Entity> & parents, const std::vector<stk::mesh::Entity> & children)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  const stk::mesh::PartVector childAddParts{&child_part()};
  const stk::mesh::PartVector childRemoveParts;
  mesh.batch_change_entity_parts(children, childAddParts, childRemoveParts);

  const stk::mesh::PartVector parentAddParts{&parent_part()};
  stk::mesh::PartVector parentRemoveParts;
  if (myActivePart)
    parentRemoveParts.push_back(myActivePart);

  mesh.batch_change_entity_parts(parents, parentAddParts, parentRemoveParts);
}

void Refinement::add_child_to_parent(const stk::mesh::EntityId childId, const stk::mesh::Entity parent)
{
  const unsigned numChildWhenFullyRefined = get_num_children_when_fully_refined(parent);
  auto * childElemIdsData = get_child_element_ids(numChildWhenFullyRefined, parent);

  unsigned iChild = 0;
  while(iChild<numChildWhenFullyRefined && childElemIdsData[iChild] != stk::mesh::InvalidEntityId)
    ++iChild;

  STK_ThrowRequireMsg(iChild < numChildWhenFullyRefined, "Logic error when assigning child id for parent element.");

  childElemIdsData[iChild] = childId;
}

void Refinement::parallel_sync_child_element_ids_fields()
{
  std::vector< const stk::mesh::FieldBase *> syncFields;
  if (myChildElementIds4Field)
    syncFields.push_back(myChildElementIds4Field);
  if (myChildElementIds8Field)
    syncFields.push_back(myChildElementIds8Field);
  stk::mesh::communicate_field_data(myMeta.mesh_bulk_data(), syncFields);
}

void Refinement::restore_child_element_ids_field(const std::vector<ParentAndChildId> & parentsAndOffProcChildId)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, myMeta.locally_owned_part() & parent_part()))
  {
    const unsigned numChildWhenFullyRefined = get_num_children_when_fully_refined(bucket->topology());
    auto * childElemIds = get_child_element_ids(numChildWhenFullyRefined, *bucket);
    std::fill(childElemIds, childElemIds+bucket->size()*numChildWhenFullyRefined, stk::mesh::InvalidEntityId);
  }

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, myMeta.locally_owned_part() & child_part()))
  {
    for (auto && child : *bucket)
    {
      stk::mesh::Entity parent = get_parent(child);
      if (mesh.is_valid(parent) && mesh.bucket(parent).owned())
      {
        add_child_to_parent(mesh.identifier(child), parent);
      }
    }
  }

  for (auto & parentAndOffProcChildId : parentsAndOffProcChildId)
    add_child_to_parent(parentAndOffProcChildId.second, parentAndOffProcChildId.first);

  parallel_sync_child_element_ids_fields();
}

void Refinement::restore_after_restart()
{
  std::vector<stk::mesh::Entity> parents;
  std::vector<stk::mesh::Entity> children;
  std::vector<ParentAndChildId> parentsAndOffProcChildId;

  fill_parents_and_children_and_parents_with_off_proc_child(parents, children, parentsAndOffProcChildId);

  restore_parent_and_child_element_parts(parents, children);

  restore_child_element_ids_field(parentsAndOffProcChildId);
}

}


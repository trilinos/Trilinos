#include <Akri_Refinement.hpp>

#include <numeric>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include "Akri_ChildNodeCreator.hpp"
#include "Akri_DiagWriter.hpp"
#include "Akri_Edge.hpp"
#include "Akri_EntityIdPool.hpp"
#include "Akri_MeshHelpers.hpp"
#include "Akri_MOAB_TetRefiner.hpp"
#include "Akri_NodeRefiner.hpp"
#include "Akri_QualityMetric.hpp"
#include "Akri_TransitionElementEdgeMarker.hpp"
#include "Akri_TriRefiner.hpp"
#include "Akri_NodeRefiner.hpp"

namespace krino {

void Refinement::declare_refinement_parts()
{
  myParentPart = &myMeta.declare_part("Refinement_Parent", stk::topology::ELEMENT_RANK, true);
  myChildPart = &myMeta.declare_part("Refinement_Child", stk::topology::ELEMENT_RANK, true);
  myRefinedEdgeNodePart = &myMeta.declare_part_with_topology("Refinement_Edge_Node", stk::topology::NODE);
}

void Refinement::declare_refinement_fields()
{
  if (3 == myMeta.spatial_dimension())
  {
    myChildElementIds8Field = myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_CHILD_ELEMENT_IDS_8");
    const stk::topology elemTopology = stk::topology::TETRAHEDRON_4;
    const stk::mesh::Part & tet4TopologyPart = myMeta.get_topology_root_part(elemTopology);
    stk::mesh::put_field_on_mesh(myChildElementIds8Field.field(), parent_part() & tet4TopologyPart, get_num_children_when_fully_refined(elemTopology), nullptr);
  }
  else if (2 == myMeta.spatial_dimension())
  {
    myChildElementIds4Field = myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_CHILD_ELEMENT_IDS_4");
    const stk::topology elemTopology = stk::topology::TRIANGLE_3_2D;
    const stk::mesh::Part & tri3TopologyPart = myMeta.get_topology_root_part(elemTopology);
    stk::mesh::put_field_on_mesh(myChildElementIds4Field.field(), parent_part() & tri3TopologyPart, get_num_children_when_fully_refined(elemTopology), nullptr);
  }

  myRefinementLevelField = myMeta.declare_field<int>(stk::topology::ELEMENT_RANK, "REFINEMENT_LEVEL");
  stk::mesh::put_field_on_mesh(myRefinementLevelField.field(), myMeta.universal_part(), 1, nullptr); // needed everywhere for restart, otherwise could be child_part

  myParentElementIdField = myMeta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "REFINEMENT_PARENT_ELEMENT_ID");
  stk::mesh::put_field_on_mesh(myParentElementIdField.field(), myMeta.universal_part(), 1, nullptr); // needed everywhere for restart, otherwise could be child_part

  myRefinedEdgeNodeParentIdsField = myMeta.declare_field<uint64_t>(stk::topology::NODE_RANK, "REFINEMENT_REFINED_EDGE_NODE_PARENTS_IDS");
  stk::mesh::put_field_on_mesh(myRefinedEdgeNodeParentIdsField.field(), refined_edge_node_part(), 2, nullptr);
}

Refinement::Refinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart, const bool force64Bit, const bool assert32Bit)
  : myMeta(meta),
    myForce64Bit(force64Bit),
    myAssert32Bit(assert32Bit),
    myNodeRefiner(force64Bit, assert32Bit),
    myEntityIdPool(meta),
    myActivePart(activePart)
{
  myCoordsField = myMeta.coordinate_field();
  declare_refinement_parts();
  declare_refinement_fields();
}

Refinement::Refinement(stk::mesh::MetaData & meta, stk::mesh::Part * activePart)
  : Refinement(meta, activePart, false, false)
{
}

Refinement::Refinement(stk::mesh::MetaData & meta)
  : Refinement(meta, nullptr, false, false)
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
  const auto * refineLevel = field_data<int>(myRefinementLevelField, elem);
  ThrowAssertMsg(refineLevel != nullptr, "Refinement level field missing on " << myMeta.mesh_bulk_data().entity_key(elem));
  return *refineLevel;
}

void Refinement::set_refinement_level(const stk::mesh::Entity elem, const int refinementLevel) const
{
  auto * refineLevel = field_data<int>(myRefinementLevelField, elem);
  ThrowAssertMsg(is_child(elem) && refineLevel != nullptr, "Refinement level field missing on " << myMeta.mesh_bulk_data().entity_key(elem));
  *refineLevel = refinementLevel;
}

unsigned Refinement::get_num_children_when_fully_refined(const stk::topology elementTopology)
{
  switch(elementTopology)
  {
  case stk::topology::TRIANGLE_3_2D:
      return 4;
  case stk::topology::TETRAHEDRON_4:
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

std::array<stk::mesh::Entity,2> Refinement::get_edge_parent_nodes(const stk::mesh::Entity edgeNode) const
{
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  auto * edgeNodeIds = field_data<uint64_t>(myRefinedEdgeNodeParentIdsField, edgeNode);
  ThrowAssertMsg(edgeNodeIds != nullptr, "Edge Node Ids field missing on node " << mesh.identifier(edgeNode));
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
  ThrowAssert(is_parent(parentElem));
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
  FieldRef childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
  const auto * childElemIdsData = field_data<uint64_t>(childElementIdsField, elem);
  ThrowAssertMsg(nullptr != childElemIdsData, "Element is parent, but does not have " << childElementIdsField.name() << " defined.");
  return std::make_tuple(childElemIdsData, numChildWhenFullyRefined);
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
  {
    const auto * parentElemIdData = field_data<uint64_t>(myParentElementIdField, elem);
    ThrowAssertMsg(nullptr != parentElemIdData, "Element is child, but does not have " << myParentElementIdField.name() << " defined.");
    parent = myMeta.mesh_bulk_data().get_entity(stk::topology::ELEMENT_RANK, *parentElemIdData);
  }
  return parent;
}

bool Refinement::is_refined_edge_node(const stk::mesh::Entity node) const
{
  return myMeta.mesh_bulk_data().bucket(node).member(refined_edge_node_part());
}

stk::mesh::Part & Refinement::child_part() const
{
  ThrowAssert(myChildPart);
  return *myChildPart;
}

stk::mesh::Part & Refinement::parent_part() const
{
  ThrowAssert(myParentPart);
  return *myParentPart;
}

stk::mesh::Part & Refinement::refined_edge_node_part() const
{
  ThrowAssert(myRefinedEdgeNodePart);
  return *myRefinedEdgeNodePart;
}

static int get_edge_refinement_case_id(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & elemChildEdgeNodes)
{
  int caseId = 0;
  for (size_t i=0; i<elemChildEdgeNodes.size(); ++i)
    if (mesh.is_valid(elemChildEdgeNodes[i]))
      caseId += 1<<i;
  return caseId;
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
    if (childSideIndices[s] >= 0)
    {
      const stk::topology sideTopology = topology.side_topology(s);
      std::vector<stk::mesh::Entity> childSideNodes(sideTopology.num_nodes());
      topology.side_nodes(childElemNodes, s, childSideNodes.data());

      std::vector<stk::mesh::Entity> sides;
      stk::mesh::get_entities_through_relations(mesh, childSideNodes, mesh.mesh_meta_data().side_rank(), sides);

      if (sides.empty())
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
      else
      {
        ThrowRequire(sides.size() == 1);
        attach_entity_to_element(mesh, mesh.mesh_meta_data().side_rank(), sides[0], childElem);
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

FieldRef Refinement::get_child_element_ids_field(const unsigned numChildWhenFullyRefined) const
{
  switch(numChildWhenFullyRefined)
  {
  case 4:
      return myChildElementIds4Field;
  case 8:
      return myChildElementIds8Field;
  default:
      ThrowRuntimeError("Unsupported number of children when fully refined in get_child_element_ids_field: " << numChildWhenFullyRefined);
      break;
  }
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

  stk::mesh::ConstPartVector removeParts;
  if (myActivePart)
    removeParts.push_back(myActivePart);
  mesh.change_entity_parts(parentElement, stk::mesh::ConstPartVector{myParentPart}, removeParts);

  const stk::mesh::EntityId parentElementId = mesh.identifier(parentElement);
  const int parentRefinementLevel = refinement_level(parentElement);
  FieldRef childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
  auto * childElemIds = field_data<uint64_t>(childElementIdsField, parentElement);
  std::fill(childElemIds, childElemIds+numChildWhenFullyRefined, stk::mesh::InvalidEntityId);
  for (size_t iChild=0; iChild<childElements.size(); ++iChild)
  {
    childElemIds[iChild] = mesh.identifier(childElements[iChild]);
    auto * parentElemId = field_data<uint64_t>(myParentElementIdField, childElements[iChild]);
    *parentElemId = parentElementId;
    set_refinement_level(childElements[iChild], parentRefinementLevel+1);
  }
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

  const std::array<stk::math::Vector3d,3> parentNodeCoords{{
    stk::math::Vector3d(field_data<double>(myCoordsField, parentElemNodes[0]), 2),
    stk::math::Vector3d(field_data<double>(myCoordsField, parentElemNodes[1]), 2),
    stk::math::Vector3d(field_data<double>(myCoordsField, parentElemNodes[2]), 2) }};

  const std::array<int,3> parentNodeRank = get_rank_of_nodes_based_on_coordinates(parentNodeCoords);

  const std::vector<TriRefiner::TriDescription> newTris = TriRefiner::refinement_child_nodes_and_sides_tri3(caseId, parentNodeCoords, parentNodeRank);

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newTris, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 4);
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
    stk::math::Vector3d(field_data<double>(myCoordsField, parentNodes[0])),
    stk::math::Vector3d(field_data<double>(myCoordsField, parentNodes[1])),
    stk::math::Vector3d(field_data<double>(myCoordsField, parentNodes[2])),
    stk::math::Vector3d(field_data<double>(myCoordsField, parentNodes[3])) }};

  const std::array<int,4> parentNodeRank= get_rank_of_nodes_based_on_coordinates(parentNodeCoords);

  const double needSides = mesh.num_sides(parentElem) > 0;
  const std::vector<moab::SimplexTemplateRefiner::TetDescription> newTets = moab::SimplexTemplateRefiner::refinement_child_nodes_and_sides_tet4(caseId, parentNodeCoords, parentNodeRank, needSides);

  std::vector<stk::mesh::Entity> childElements;
  declare_child_elements_and_append_sides_to_create(mesh, myEntityIdPool, childParts, parentElem, parentElemNodes, newTets, childElements, sideRequests);

  set_parent_parts_and_parent_child_relation_fields(parentElem, childElements, 8);
}

static unsigned num_new_child_elements_for_case_id(const stk::topology & elemTopology, const int caseId)
{
  switch(elemTopology())
    {
    case stk::topology::TRI_3:
    case stk::topology::TRI_3_2D:
        return TriRefiner::num_new_child_elements_tri3(caseId);
    case stk::topology::TETRAHEDRON_4:
        return moab::SimplexTemplateRefiner::num_new_child_elements_tet4(caseId);
    default:
        ThrowRuntimeError("Element topology not found in refine_element: " << elemTopology.name());
    }
}

void Refinement::refine_element_if_it_has_refined_edges_and_append_sides_to_create(const stk::topology & elemTopology,
    const stk::mesh::PartVector & childParts,
    const stk::mesh::Entity elem,
    const std::vector<stk::mesh::Entity> & elemChildEdgeNodes,
    std::vector<SideDescription> & sideRequests,
    std::vector<stk::mesh::Entity> & elementsToDelete)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const int caseId = get_edge_refinement_case_id(mesh, elemChildEdgeNodes);
  if (0 == caseId)
    return;

  const std::vector<stk::mesh::Entity> existingChildrenToDelete = get_children(elem);

  switch(elemTopology())
    {
    case stk::topology::TRI_3:
    case stk::topology::TRI_3_2D:
        refine_tri_3_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, caseId, sideRequests);
        break;
    case stk::topology::TETRAHEDRON_4:
        refine_tet_4_and_append_sides_to_create(childParts, elem, elemChildEdgeNodes, caseId, sideRequests);
        break;
    default:
        ThrowRuntimeError("Element topology not found in refine_element: " << elemTopology.name());
        break;
    }

  elementsToDelete.insert(elementsToDelete.end(), existingChildrenToDelete.begin(), existingChildrenToDelete.end());
}

size_t Refinement::count_new_child_elements(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData) const
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();

  std::vector<stk::mesh::Entity> elemChildEdgeNodes;

  size_t numNewChildElems = 0;
  for(const auto & bucketData : bucketsData)
  {
    for(const auto & elem : std::get<2>(bucketData))
    {
      const stk::topology bucketTopology = std::get<0>(bucketData);
      edgeMarker.fill_element_refined_edge_nodes(myNodeRefiner, elem, bucketTopology, elemChildEdgeNodes);
      const int caseId = get_edge_refinement_case_id(mesh, elemChildEdgeNodes);
      if (0 != caseId)
        numNewChildElems += num_new_child_elements_for_case_id(bucketTopology, caseId);
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

std::vector<Refinement::BucketData> Refinement::get_buckets_data_for_candidate_elements_to_refine(const EdgeMarkerInterface & edgeMarker) const
{
  // Cache off bucket data to avoid looping over buckets while modifying elements
  const stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part();

  const stk::mesh::EntityVector emptyVector;
  std::vector<std::tuple<stk::topology,stk::mesh::PartVector,stk::mesh::EntityVector>> bucketsData;

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, selector))
  {
    const stk::mesh::PartVector childParts = get_parts_for_child_elements(*bucketPtr);
    bucketsData.emplace_back(bucketPtr->topology(), childParts, emptyVector);

    stk::mesh::EntityVector & bucketElements = std::get<2>(bucketsData.back());
    bucketElements.reserve(bucketPtr->size());
    for (auto && elem : *bucketPtr)
      if (edgeMarker.is_element_a_candidate_for_refinement(elem))
        bucketElements.push_back(elem);
  }
  return bucketsData;
}

void Refinement::refine_elements_with_refined_edges_and_store_sides_to_create(const EdgeMarkerInterface & edgeMarker, const std::vector<BucketData> & bucketsData, std::vector<SideDescription> & sideRequests, std::vector<stk::mesh::Entity> & elementsToDelete)
{
  std::vector<stk::mesh::Entity> elemChildEdgeNodes;
  //for(const auto & [bucketTopology, bucketChildParts, bucketElements] : bucketsData)
  for(const auto & bucketData : bucketsData)
  {
    const auto & bucketTopology = std::get<0>(bucketData);
    const auto & bucketChildParts = std::get<1>(bucketData);
    const auto & bucketElements = std::get<2>(bucketData);
    for(const auto & elem : bucketElements)
    {
      edgeMarker.fill_element_refined_edge_nodes(myNodeRefiner, elem, bucketTopology, elemChildEdgeNodes);
      refine_element_if_it_has_refined_edges_and_append_sides_to_create(bucketTopology, bucketChildParts, elem, elemChildEdgeNodes, sideRequests, elementsToDelete);
    }
  }
}

stk::mesh::PartVector Refinement::get_parts_for_new_refined_edge_nodes() const
{
  stk::mesh::PartVector refinedEdgeNodeParts = { &myMeta.get_topology_root_part(stk::topology::NODE), &refined_edge_node_part() };
  if (myActivePart)
    refinedEdgeNodeParts.push_back(myActivePart);
  return refinedEdgeNodeParts;
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
          krinolog << "Found element with hanging node on edge " << debug_edge(mesh, edge) << " of element " << debug_entity_1line(mesh, elem) << stk::diag::dendl;;
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

void Refinement::create_refined_nodes_elements_and_sides(const EdgeMarkerInterface & edgeMarker)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  mesh.modification_begin();

  myNodeRefiner.create_refined_edge_nodes(mesh, get_parts_for_new_refined_edge_nodes(), myRefinedEdgeNodeParentIdsField);

  const std::vector<BucketData> bucketsData = get_buckets_data_for_candidate_elements_to_refine(edgeMarker);
  const size_t numNewElements = count_new_child_elements(edgeMarker, bucketsData);
  myEntityIdPool.reserve(stk::topology::ELEMENT_RANK, numNewElements, myAssert32Bit, myForce64Bit);

  std::vector<stk::mesh::Entity> elementsToDelete;
  std::vector<SideDescription> sideRequests;
  refine_elements_with_refined_edges_and_store_sides_to_create(edgeMarker, bucketsData, sideRequests, elementsToDelete);
  delete_mesh_entities(mesh, elementsToDelete);
  mesh.modification_end();

  if(stk::is_true_on_any_proc(mesh.parallel(), !sideRequests.empty()))
    batch_create_sides(mesh, sideRequests);

  myNodeRefiner.prolong_refined_edge_nodes(mesh);
}

void Refinement::create_another_layer_of_refined_elements_and_sides_to_eliminate_hanging_nodes(const EdgeMarkerInterface & edgeMarker)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();
  const std::vector<BucketData> bucketsData = get_buckets_data_for_candidate_elements_to_refine(edgeMarker);
  const size_t numNewElements = count_new_child_elements(edgeMarker, bucketsData);

  myEntityIdPool.reserve(stk::topology::ELEMENT_RANK, numNewElements, myAssert32Bit, myForce64Bit);

  std::vector<stk::mesh::Entity> elementsToDelete;
  std::vector<SideDescription> sideRequests;

  if(stk::is_true_on_any_proc(mesh.parallel(), numNewElements > 0))
  {
    mesh.modification_begin();
    refine_elements_with_refined_edges_and_store_sides_to_create(edgeMarker, bucketsData, sideRequests, elementsToDelete);
    delete_mesh_entities(mesh, elementsToDelete);
    mesh.modification_end();

    if(stk::is_true_on_any_proc(mesh.parallel(), !sideRequests.empty()))
      batch_create_sides(mesh, sideRequests);
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

void Refinement::mark_already_refined_edges()
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  myNodeRefiner.clear_edges_to_refine();

  const stk::mesh::Selector selector = refined_edge_node_part();

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, selector))
  {
    for (auto && existingRefinedNode : *bucketPtr)
    {
      const auto edgeNodeParents = get_edge_parent_nodes(existingRefinedNode);
      if (mesh.is_valid(edgeNodeParents[0]) && mesh.is_valid(edgeNodeParents[1]))
      {
        myNodeRefiner.mark_already_refined_edge(edge_from_edge_nodes(mesh, edgeNodeParents[0], edgeNodeParents[1]), existingRefinedNode);
      }
    }
  }
}

void Refinement::do_unrefinement(const EdgeMarkerInterface & edgeMarker)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  if(stk::is_true_on_any_proc(mesh.parallel(), edgeMarker.locally_have_elements_to_unrefine()))
  {
    std::vector<stk::mesh::Entity> childElementsToDeleteForUnrefinement;
    std::vector<stk::mesh::Entity> elementsWithoutChildrenAfterUnrefinement;
    edgeMarker.fill_elements_modified_by_unrefinement(elementsWithoutChildrenAfterUnrefinement, childElementsToDeleteForUnrefinement);

    mesh.modification_begin();
    remove_parent_parts(elementsWithoutChildrenAfterUnrefinement);
    delete_mesh_entities(mesh, childElementsToDeleteForUnrefinement);
    mesh.modification_end();

    mark_already_refined_edges();

    create_another_layer_of_refined_elements_and_sides_to_eliminate_hanging_nodes(edgeMarker);
  }
}

void Refinement::do_refinement(const EdgeMarkerInterface & edgeMarker)
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  find_edges_to_refine(edgeMarker);

  const bool haveEdgesToRefineLocally = get_num_edges_to_refine() > 0;
  if(stk::is_true_on_any_proc(mesh.parallel(), haveEdgesToRefineLocally))
  {
    create_refined_nodes_elements_and_sides(edgeMarker);

    create_another_layer_of_refined_elements_and_sides_to_eliminate_hanging_nodes(edgeMarker);

    if (myActivePart)
      activate_selected_sides_touching_active_elements(mesh, myMeta.universal_part(), *myActivePart);

    myNodeRefiner.prolong_refined_edge_nodes(mesh);
  }

  do_unrefinement(edgeMarker);

  ThrowAssertMsg(!have_any_hanging_refined_nodes(), "Mesh has hanging refined node.");
}

void Refinement::do_uniform_refinement(const int numUniformRefinementLevels)
{
  UniformEdgeMarker uniformMarker(myMeta.mesh_bulk_data(), *this);
  for (int i=0; i<numUniformRefinementLevels; ++i)
    do_refinement(uniformMarker);
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
  delete_mesh_entities(mesh, allChildElems);
  mesh.change_entity_parts(ownedParentElems, addParts, stk::mesh::ConstPartVector{myParentPart});
  mesh.modification_end();
}

void Refinement::find_edges_to_refine(const EdgeMarkerInterface & edgeMarker)
{
  edgeMarker.mark_edges_to_be_refined(myNodeRefiner);
}

void Refinement::restore_parent_and_child_element_parts_after_restart()
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  std::vector<stk::mesh::Entity> children;
  std::vector<stk::mesh::Entity> parents;
  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, myMeta.locally_owned_part()))
  {
    for (auto && elem : *bucket)
    {
      const auto * parentElemIdData = field_data<uint64_t>(myParentElementIdField, elem);
      ThrowAssertMsg(nullptr != parentElemIdData, "Element is does not have " << myParentElementIdField.name() << " defined.");
      const stk::mesh::Entity parent = myMeta.mesh_bulk_data().get_entity(stk::topology::ELEMENT_RANK, *parentElemIdData);
      if (mesh.is_valid(parent))
      {
        children.push_back(elem);
        parents.push_back(parent);
      }
    }
  }

  stk::util::sort_and_unique(parents);

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
  FieldRef childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
  auto * childElemIdsData = field_data<uint64_t>(childElementIdsField, parent);

  unsigned iChild = 0;
  while(iChild<numChildWhenFullyRefined && childElemIdsData[iChild] != stk::mesh::InvalidEntityId)
    ++iChild;

  ThrowRequireMsg(iChild < numChildWhenFullyRefined, "Logic error when assigning child id for parent element.");

  childElemIdsData[iChild] = childId;
}

void Refinement::parallel_sync_child_element_ids_fields()
{
  std::vector< const stk::mesh::FieldBase *> syncFields;
  if (myChildElementIds4Field.valid())
    syncFields.push_back(&myChildElementIds4Field.field());
  if (myChildElementIds8Field.valid())
    syncFields.push_back(&myChildElementIds8Field.field());
  stk::mesh::communicate_field_data(myMeta.mesh_bulk_data(), syncFields);
}

void Refinement::restore_child_element_ids_field_after_restart()
{
  stk::mesh::BulkData & mesh = myMeta.mesh_bulk_data();

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, myMeta.locally_owned_part() & parent_part()))
  {
    const unsigned numChildWhenFullyRefined = get_num_children_when_fully_refined(bucket->topology());
    FieldRef childElementIdsField = get_child_element_ids_field(numChildWhenFullyRefined);
    auto * childElemIds = field_data<uint64_t>(childElementIdsField, *bucket);
    std::fill(childElemIds, childElemIds+bucket->size()*numChildWhenFullyRefined, stk::mesh::InvalidEntityId);
  }

  for (auto && bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, myMeta.locally_owned_part() & child_part()))
  {
    for (auto && child : *bucket)
    {
      stk::mesh::Entity parent = get_parent(child);
      if (mesh.is_valid(parent))
      {
        add_child_to_parent(mesh.identifier(child), parent);
      }
    }
  }

  parallel_sync_child_element_ids_fields();
}

void Refinement::restore_after_restart()
{
  restore_parent_and_child_element_parts_after_restart();

  restore_child_element_ids_field_after_restart();
}

}


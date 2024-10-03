#include <Akri_DetermineElementSign.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Edge.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Surface.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino {

typedef std::pair<stk::mesh::Entity, int8_t> EntityAndSign;

static ElementToSignsMap initialize_element_signs(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector, const size_t numSurfaces)
{
  const std::vector<int8_t> initVals(numSurfaces, -2);
  ElementToSignsMap elementsToSigns;
  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector))
    for (const auto & elem : *bucketPtr)
      elementsToSigns[elem] = initVals;
  return elementsToSigns;
}

static void assign_element_sign_if_not_already_assigned(const stk::mesh::BulkData & mesh,
    const unsigned surfIndex,
    const stk::mesh::Entity element,
    const int8_t sign,
    ElementToSignsMap & elementsToSigns,
    bool & isNewAssignment)
{
  int8_t & elementSign = elementsToSigns.at(element)[surfIndex];
//  const int8_t prevSign = elementSign; //debug
  isNewAssignment = (elementSign == -2);

//  krinolog << "  Attempting to assign element " << mesh.identifier(element) << " " << int(prevSign) << " + " << int(sign);

  if (elementSign != 0)
  {
//    if (!isNewAssignment && elementSign != sign)
//      krinolog << "\nInconceivable " << debug_entity(mesh, element) << stk::diag::dendl;
    STK_ThrowRequire(isNewAssignment || elementSign == sign);
    elementSign = sign;
  }

//  krinolog << " -> " << int(elementSign) << stk::diag::dendl;
}

static void assign_interface_element_sign(const unsigned surfIndex,
    const stk::mesh::Entity element,
    ElementToSignsMap & elementsToSigns)
{
  int8_t & elementSign = elementsToSigns.at(element)[surfIndex];
  STK_ThrowRequire(elementSign == -2 || elementSign == 0);
  elementSign = 0;
}

static void fill_edge_elements(const stk::mesh::BulkData & mesh, const Edge & edge, std::vector<stk::mesh::Entity> & edgeElements)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  stk::mesh::get_entities_through_relations(mesh, stk::mesh::EntityVector{edgeNodes[0], edgeNodes[1]}, stk::topology::ELEMENT_RANK, edgeElements);
}

static void check_edge_intersections_to_assign_crossed_elements_and_find_nodes_on_either_side_of_surface(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const unsigned surfIndex,
    const Surface & surface,
    ElementToSignsMap & elementsToSigns,
    std::set<EntityAndSign> & nodesAndSigns)
{
  const double largeEdgeTolToKeepCostLow = 0.1;
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<stk::mesh::Entity> edgeElements;
  nodesAndSigns.clear();

  const std::vector<Edge> edges = get_edges_of_selected_elements(mesh, elementSelector);
  for (auto & edge : edges)
  {
    const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
    //krinolog << "Checking  edge " << mesh.identifier(edgeNodes[0]) << " " << mesh.identifier(edgeNodes[1]) << ": " << stk::diag::dendl;
    const auto [crossingSign, position] = surface.compute_intersection_with_segment(get_vector_field(mesh, coordsField, edgeNodes[0], dim), get_vector_field(mesh, coordsField, edgeNodes[1], dim), largeEdgeTolToKeepCostLow);
    if (crossingSign != 0)
    {
      nodesAndSigns.emplace(edgeNodes[0], -crossingSign);
      nodesAndSigns.emplace(edgeNodes[1], crossingSign);

      fill_edge_elements(mesh, edge, edgeElements);
      for (auto edgeElem : edgeElements)
        if (elementSelector(mesh.bucket(edgeElem)))
          assign_interface_element_sign(surfIndex, edgeElem, elementsToSigns);
    }
  }
}

static void assign_element_sign_and_append_nodes_if_new(const stk::mesh::BulkData & mesh, const unsigned surfIndex, const stk::mesh::Entity element, const int8_t sign, std::set<EntityAndSign> & nextIterNodesAndSigns, ElementToSignsMap & elementsToSigns)
{
  bool isNewAssignment = false;
  assign_element_sign_if_not_already_assigned(mesh, surfIndex, element, sign, elementsToSigns, isNewAssignment);

  if (isNewAssignment)
    for (auto elementNode : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
      nextIterNodesAndSigns.emplace(elementNode,  sign);
}

static
void pack_unowned_elements_and_signs_for_owners(const stk::mesh::BulkData & mesh,
    const std::set<EntityAndSign> & elemsToCommunicate,
    stk::CommSparse &commSparse)
{
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto & [elem, sign] : elemsToCommunicate)
    {
      if (!mesh.bucket(elem).owned())
      {
        const int owner = mesh.parallel_owner_rank(elem);
        commSparse.send_buffer(owner).pack(mesh.identifier(elem));
        commSparse.send_buffer(owner).pack(sign);
      }
    }
  });
}

static
void unpack_elements_and_signs(const stk::mesh::BulkData & mesh,
    std::set<EntityAndSign> & communicatedElementsAndSigns,
    stk::CommSparse &commSparse)
{
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId elemId;
      commSparse.recv_buffer(procId).unpack(elemId);
      stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEMENT_RANK, elemId);
      STK_ThrowRequire(mesh.is_valid(elem));
      int8_t sign = 0;
      commSparse.recv_buffer(procId).unpack(sign);
      communicatedElementsAndSigns.emplace(elem, sign);
    }
  });
}

static
void pack_owned_elements_and_signs_for_ghosting_procs(const stk::mesh::BulkData & mesh,
    const std::set<EntityAndSign> & ownedElemsAndSigns,
    stk::CommSparse &commSparse)
{
  std::vector<int> elemCommProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto & [elem, sign] : ownedElemsAndSigns)
    {
      if (mesh.bucket(elem).owned())
      {
        mesh.comm_procs(elem, elemCommProcs);
        for (int procId : elemCommProcs)
        {
          if (procId != commSparse.parallel_rank())
          {
            commSparse.send_buffer(procId).pack(mesh.identifier(elem));
            commSparse.send_buffer(procId).pack(sign);
          }
        }
      }
    }
  });
}

std::set<EntityAndSign> communicate_elements_and_signs(const stk::mesh::BulkData & mesh, const std::set<EntityAndSign> & elemsToCommunicate)
{
  std::set<EntityAndSign> communicatedElementsAndSigns;
  for (auto & [elem, sign] : elemsToCommunicate)
      if (mesh.bucket(elem).owned())
        communicatedElementsAndSigns.emplace(elem, sign);

  {
    stk::CommSparse commSparse(mesh.parallel());
    pack_unowned_elements_and_signs_for_owners(mesh, elemsToCommunicate, commSparse);
    unpack_elements_and_signs(mesh, communicatedElementsAndSigns, commSparse);
  }

  {
    stk::CommSparse commSparse(mesh.parallel());
    pack_owned_elements_and_signs_for_ghosting_procs(mesh, communicatedElementsAndSigns, commSparse);
    unpack_elements_and_signs(mesh, communicatedElementsAndSigns, commSparse);
  }

  return communicatedElementsAndSigns;
}

void assign_elements_on_either_side_of_surface(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const unsigned surfIndex,
    ElementToSignsMap & elementsToSigns,
    std::set<EntityAndSign> & nodesAndSigns)
{
  while (stk::is_true_on_any_proc(mesh.parallel(), !nodesAndSigns.empty()))
  {
    std::set<EntityAndSign> nextIterNodesAndSigns;
    std::set<EntityAndSign> elemsToCommunicate;

    for (auto & [node, sign] : nodesAndSigns)
    {
      //krinolog << "Processing node " << mesh.identifier(node) << " with sign " << int(sign) << stk::diag::dendl;
      for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
      {
        const stk::mesh::Bucket & elemBucket = mesh.bucket(elem);
        if (elementSelector(elemBucket))
        {
          if (!elemBucket.owned() || mesh.in_send_ghost(elem))
            elemsToCommunicate.emplace(elem, sign);

          assign_element_sign_and_append_nodes_if_new(mesh, surfIndex, elem, sign, nextIterNodesAndSigns, elementsToSigns);
        }
      }
    }

    const std::set<EntityAndSign> elemsFromOtherProcs = communicate_elements_and_signs(mesh, elemsToCommunicate);
    for (auto & [element, sign] : elemsFromOtherProcs)
    {
      assign_element_sign_and_append_nodes_if_new(mesh, surfIndex, element, sign, nextIterNodesAndSigns, elementsToSigns);
    }

    nodesAndSigns.swap(nextIterNodesAndSigns);
  }
}

ElementToSignsMap determine_element_signs(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Selector> & perSurfaceElementSelector,
    const std::vector<const Surface*> & surfaces)
{
  STK_ThrowAssert(perSurfaceElementSelector.size() == surfaces.size());
  ElementToSignsMap elementsToSigns = initialize_element_signs(mesh, elementSelector, surfaces.size());

  for (unsigned surfIndex=0; surfIndex<surfaces.size(); ++surfIndex)
  {
    stk::mesh::Selector elemSurfSelector = elementSelector & perSurfaceElementSelector[surfIndex];
    std::set<EntityAndSign> nodesAndSigns;
    check_edge_intersections_to_assign_crossed_elements_and_find_nodes_on_either_side_of_surface(mesh, coordsField, elemSurfSelector, surfIndex, *surfaces[surfIndex], elementsToSigns, nodesAndSigns);
    assign_elements_on_either_side_of_surface(mesh, elemSurfSelector, surfIndex, elementsToSigns, nodesAndSigns);
  }

  return elementsToSigns;
}

}


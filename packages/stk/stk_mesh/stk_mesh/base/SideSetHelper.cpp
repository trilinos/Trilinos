#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SideSetHelper.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include <stk_util/parallel/ParallelReduce.hpp>
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemGraphCoincidentElems.hpp"

namespace stk
{
namespace mesh
{

void SideSetHelper::remove_element_entries_from_sidesets(SideSetVector& sidesets, const Entity entity, std::set<const Part*> *touchedSidesetParts)
{
  if (mesh.entity_rank(entity) == stk::topology::ELEMENT_RANK && mesh.num_sides(entity) > 0) {
    for (SideSet* sideset : sidesets) {
      std::vector<SideSetEntry>::iterator lowerBound = std::lower_bound(sideset->begin(), sideset->end(), SideSetEntry(entity, 0));
      std::vector<SideSetEntry>::iterator upperBound = std::upper_bound(sideset->begin(), sideset->end(), SideSetEntry(entity, INVALID_CONNECTIVITY_ORDINAL));

      sideset->erase(lowerBound, upperBound);

      if(nullptr != touchedSidesetParts) {
        const Part* part = sideset->get_part();

        if(nullptr != part) {
          touchedSidesetParts->insert(part);
        }
      }
    }
  }
}

void SideSetHelper::remove_element_entries_from_sidesets(const Entity entity, std::set<const Part*> *touchedSidesetParts)
{
  SideSetVector sidesets = mesh.get_sidesets();
  remove_element_entries_from_sidesets(sidesets, entity, touchedSidesetParts);
}

void SideSetHelper::remove_side_entry_from_sideset(SideSet* sideset, const SideSetEntry& entry,
                                                   std::set<const Part*> *touchedSidesetParts)
{
  STK_ThrowAssert(entry.side != INVALID_CONNECTIVITY_ORDINAL && mesh.entity_rank(entry.element) == stk::topology::ELEM_RANK);

  std::vector<SideSetEntry>::iterator lowerBound = std::lower_bound(sideset->begin(), sideset->end(), entry);

  if(lowerBound != sideset->end() && *lowerBound == entry) {
    sideset->erase(lowerBound);
  }

  if(nullptr != touchedSidesetParts) {
    const Part* part = sideset->get_part();

    if(nullptr != part) {
      touchedSidesetParts->insert(part);
    }
  }
}

void SideSetHelper::remove_side_entry_from_sidesets(SideSetVector& sidesets, Entity elem, ConnectivityOrdinal ordinal,
                                                    std::set<const Part*> *touchedSidesetParts)
{
  if(ordinal != INVALID_CONNECTIVITY_ORDINAL && mesh.entity_rank(elem) == stk::topology::ELEM_RANK) {
    SideSetEntry entry(elem, ordinal);
    for (SideSet* sideset : sidesets) {
      remove_side_entry_from_sideset(sideset, entry, touchedSidesetParts);
    }
  }
}

void SideSetHelper::remove_side_entry_from_sidesets(Entity elem, ConnectivityOrdinal ordinal, std::set<const Part*> *touchedSidesetParts)
{
    SideSetVector sidesets = mesh.get_sidesets();
    remove_side_entry_from_sidesets(sidesets, elem, ordinal, touchedSidesetParts);
}

void SideSetHelper::remove_side_entries_from_sidesets(SideSetVector& sidesets, const Entity entity, std::set<const Part*> *touchedSidesetParts)
{
  if (sidesets.size() > 0 && mesh.entity_rank(entity) == mesh.mesh_meta_data().side_rank()) {
    unsigned numElements = mesh.num_elements(entity);
    const Entity* elems = mesh.begin_elements(entity);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(entity);

    for(unsigned i=0; i<numElements; ++i) {
      remove_side_entry_from_sidesets(sidesets, elems[i],  ordinals[i], touchedSidesetParts);
    }
  }
}

void SideSetHelper::remove_side_entries_from_sidesets(const Entity entity, std::set<const Part*> *touchedSidesetParts)
{
  SideSetVector sidesets = mesh.get_sidesets();
  remove_side_entries_from_sidesets(sidesets, entity, touchedSidesetParts);
}

void SideSetHelper::remove_entity_entries_from_sidesets(SideSetVector& sidesets, const Entity entity, std::set<const Part*> *touchedSidesetParts)
{
  EntityRank rank = mesh.entity_rank(entity);

  if (rank == stk::topology::ELEMENT_RANK){
    remove_element_entries_from_sidesets(sidesets, entity, touchedSidesetParts);
  }
  else if(rank == mesh.mesh_meta_data().side_rank()) {
    remove_side_entries_from_sidesets(sidesets, entity, touchedSidesetParts);
  }
}

void SideSetHelper::remove_entity_entries_from_sidesets(const Entity entity, std::set<const Part*> *touchedSidesetParts)
{
  SideSetVector sidesets = mesh.get_sidesets();
  remove_entity_entries_from_sidesets(sidesets, entity, touchedSidesetParts);
}

Entity get_side(const BulkData& bulk, Entity elem, ConnectivityOrdinal ordinal)
{
  const ConnectivityOrdinal* ordinals = bulk.begin_ordinals(elem, bulk.mesh_meta_data().side_rank());
  const Entity* sides = bulk.begin(elem, bulk.mesh_meta_data().side_rank());
  unsigned numSides = bulk.num_connectivity(elem, bulk.mesh_meta_data().side_rank());
  Entity side;

  for(unsigned i=0; i<numSides; ++i) {
    if(ordinals[i] == ordinal) {
      side = sides[i];
    }
  }

  return side;
}

void SideSetHelper::add_element_side_entry_to_sideset(SideSetSelector& sidesetSelector, Entity elem, ConnectivityOrdinal ordinal)
{
  if(ordinal != INVALID_CONNECTIVITY_ORDINAL && mesh.entity_rank(elem) == stk::topology::ELEM_RANK) {
    SideSet* sideset = sidesetSelector.sideset();
    const Selector* selector = sidesetSelector.selector();

    STK_ThrowAssertMsg(nullptr != selector, "NULL Element block selector for sideset: " << sideset->get_name());

    if(mesh.bucket(elem).owned() && (*selector)(mesh.bucket(elem))) {
      sideset->add(elem, ordinal);
    }
  }
}

void SideSetHelper::add_element_side_entry_to_sidesets(SideSetSelectorVector& sidesets, Entity elem, ConnectivityOrdinal ordinal)
{
  for (SideSetSelector& sidesetSelector : sidesets) {
    add_element_side_entry_to_sideset(sidesetSelector, elem, ordinal);
  }
}

void SideSetHelper::fill_coincident_sideset_entries_for_side(const Entity side, std::vector<SideSetEntry>& entries)
{
  if(mesh.has_face_adjacent_element_graph()) {
    fill_coincident_sideset_entries_for_side_using_elem_elem_graph(side, entries);
  } else {
    fill_coincident_sideset_entries_for_side_using_connectivity(side, entries);
  }
}

void SideSetHelper::fill_coincident_sideset_entries_for_side_using_elem_elem_graph(const Entity side, std::vector<SideSetEntry>& entries)
{
  entries.clear();

  if(mesh.entity_rank(side) == mesh.mesh_meta_data().side_rank()) {

    unsigned numElements = mesh.num_elements(side);
    const Entity* elems = mesh.begin_elements(side);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(side);

    for(unsigned i=0; i<numElements; ++i) {
      if(element_side_has_coincidence_using_elem_elem_graph(side, elems[i], ordinals[i])) {
          stk::util::insert_keep_sorted_and_unique(SideSetEntry(elems[i], ordinals[i]), entries);
      }
    }
  }
}

void SideSetHelper::fill_coincident_sideset_entries_for_side_using_connectivity(const Entity side, std::vector<SideSetEntry>& entries)
{
  entries.clear();

  if(mesh.entity_rank(side) == mesh.mesh_meta_data().side_rank()) {
    unsigned numElements = mesh.num_elements(side);
    const Entity* elems = mesh.begin_elements(side);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(side);

    if(numElements == 0) {
      return;
    }

    if(numElements == 1) {
      entries.emplace_back(elems[0], ordinals[0]);
      return;
    }

    unsigned numPossibleEntries = numElements*(numElements-1)/2;
    std::vector<std::pair<SideSetEntry, SideSetEntry>> coincidentEdges(numPossibleEntries);

    EntityVector sideNodes_i, sideNodes_j;
    for(unsigned i=0; i<numElements; ++i) {
      impl::fill_element_side_nodes_from_topology(mesh.bucket(elems[i]).topology(), mesh.begin_nodes(elems[i]), ordinals[i], sideNodes_i);

      for(unsigned j=i+1; j<numElements; ++j) {
        impl::fill_element_side_nodes_from_topology(mesh.bucket(elems[j]).topology(), mesh.begin_nodes(elems[j]), ordinals[j], sideNodes_j);

        if(impl::is_coincident_connection(mesh, elems[i], sideNodes_i, ordinals[i], mesh.bucket(elems[j]).topology(), sideNodes_j, ordinals[j])) {
          coincidentEdges.emplace_back(std::make_pair(SideSetEntry(elems[i], ordinals[i]), SideSetEntry(elems[j], ordinals[j])));
        }
      }
    }

    for(const std::pair<SideSetEntry, SideSetEntry>& coincidentEdge : coincidentEdges) {
      stk::util::insert_keep_sorted_and_unique(coincidentEdge.first, entries);
      stk::util::insert_keep_sorted_and_unique(coincidentEdge.second, entries);
    }
  }
}


bool SideSetHelper::graph_edge_can_be_distinguished(const ElemElemGraph& eeGraph, const GraphEdge& graphEdge, const Selector& selector, bool selectorValue)
{
  bool isParallelEdge = !impl::is_local_element(graphEdge.elem2());

  if(!isParallelEdge) {
    // Local connectivity
    Entity otherElement = eeGraph.get_entity(graphEdge.elem2());
    bool elemSelectorValue = selector(mesh.bucket(otherElement)) && activeSelector(mesh.bucket(otherElement));
    if(selectorValue != elemSelectorValue) {
      return true;
    }
  } else {
    // Remote connectivity
    const auto iter = parallelPartInfo.find(graphEdge.elem2());
    if(iter == parallelPartInfo.end()) { return false; }

    const std::vector<stk::mesh::PartOrdinal>& otherElementPartOrdinals = iter->second.elementPartOrdinals;

    stk::mesh::EntityId otherElemId = eeGraph.convert_negative_local_id_to_global_id(graphEdge.elem2());
    stk::mesh::Part* otherElementBlockPart = get_element_block_part(mesh, otherElementPartOrdinals, otherElemId);

    bool elemSelectorValue = selector(otherElementBlockPart) && activeSelector(otherElementBlockPart);
    if(selectorValue != elemSelectorValue) {
      return true;
    }
  }

  return false;
}

bool SideSetHelper::element_side_can_be_distinguished_using_elem_elem_graph(const Entity /*side*/, const Entity element,
                                                                            const ConnectivityOrdinal ordinal, const SideSetSelector& sideset)
{
  if(mesh.entity_rank(element) == stk::topology::ELEM_RANK && mesh.bucket(element).owned()) {
    const ElemElemGraph& eeGraph = mesh.get_face_adjacent_element_graph();
    impl::LocalId elemLocalId = eeGraph.get_local_element_id(element);

    const Selector& selector = *sideset.selector();
    bool selectorValue = selector(mesh.bucket(element)) && activeSelector(mesh.bucket(element));

    GraphEdgesForElement graphEdges = eeGraph.get_edges_for_element(elemLocalId);
    for(size_t i = 0; i < graphEdges.size(); ++i)
    {
      const GraphEdge& graphEdge =  eeGraph.get_graph().get_edge_for_element(elemLocalId, i);
      Ordinal sideOrdinal = graphEdge.side1();

      if(ordinal == sideOrdinal) {
        if(graph_edge_can_be_distinguished(eeGraph, graphEdge, selector, selectorValue)) {
          return true;
        }
      }
    }

    internalSidesetOrdinals.insert(sideset.part().mesh_meta_data_ordinal());
  }

  return false;
}

bool SideSetHelper::element_side_can_be_distinguished_using_connectivity(const Entity side, const Entity element,
                                                                         const ConnectivityOrdinal ordinal, const SideSetSelector& sideset)
{
  if(mesh.entity_rank(element) == stk::topology::ELEM_RANK) {
    unsigned numElements = mesh.num_elements(side);
    const Entity* elems = mesh.begin_elements(side);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(side);

    if(numElements == 0) {
      return true;
    }

    if(numElements == 1 && elems[0] == element && ordinals[0] == ordinal) {
      return true;
    }

    const Selector& selector = *sideset.selector();
    bool selectorValue = selector(mesh.bucket(element)) && activeSelector(mesh.bucket(element));

    for(unsigned i=0; i<numElements; ++i) {
      if(elems[i] != element) {
        bool elemSelectorValue = selector(mesh.bucket(elems[i])) && activeSelector(mesh.bucket(elems[i]));
        if(selectorValue != elemSelectorValue) {
          return true;
        }
      }
    }

    internalSidesetOrdinals.insert(sideset.part().mesh_meta_data_ordinal());
  }
  return false;
}

bool SideSetHelper::element_side_has_local_coincidence_using_elem_elem_graph(Entity element,  ConnectivityOrdinal ordinal)
{
  bool hasLocalCoincidence = false;

  if(mesh.entity_rank(element) == stk::topology::ELEM_RANK && mesh.bucket(element).owned()) {
    const stk::mesh::ElemElemGraph& eeGraph = mesh.get_face_adjacent_element_graph();

    stk::mesh::impl::LocalId elemLocalId = eeGraph.get_local_element_id(element);
    std::vector<stk::mesh::impl::LocalId> localCoincidents = eeGraph.get_coincident_elements_on_side(elemLocalId, ordinal);

    hasLocalCoincidence = localCoincidents.size() > 0;
  }

  return hasLocalCoincidence;
}

bool SideSetHelper::element_side_has_remote_coincidence_using_elem_elem_graph(Entity element,  ConnectivityOrdinal ordinal)
{
  bool hasRemoteCoincidence = false;

  if(mesh.entity_rank(element) == stk::topology::ELEM_RANK && mesh.bucket(element).owned()) {
    const stk::mesh::ElemElemGraph& eeGraph = mesh.get_face_adjacent_element_graph();

    stk::mesh::impl::LocalId elemLocalId = eeGraph.get_local_element_id(element);

    stk::mesh::EntityRank sideRank = mesh.mesh_meta_data().side_rank();
    EntityVector sideNodes;
    impl::fill_element_side_nodes_from_topology(mesh.bucket(element).topology(), mesh.begin_nodes(element), ordinal, sideNodes);
    stk::mesh::OrdinalAndPermutation connectedOrdAndPerm = stk::mesh::get_ordinal_and_permutation(mesh, element, sideRank, sideNodes);
    STK_ThrowRequire(connectedOrdAndPerm.second != INVALID_PERMUTATION);
    bool localPolarity = mesh.bucket(element).topology().sub_topology(sideRank, connectedOrdAndPerm.first).is_positive_polarity(connectedOrdAndPerm.second);

    for(const stk::mesh::GraphEdge & graphEdge : eeGraph.get_edges_for_element(elemLocalId)) {
      bool isParallelEdge = !impl::is_local_element(graphEdge.elem2());

      if (isParallelEdge) {
        const impl::ParallelInfo &parallelInfo = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
        bool remotePolarity = parallelInfo.m_remote_element_topology.sub_topology(sideRank, graphEdge.side2()).is_positive_polarity(parallelInfo.m_permutation);

        hasRemoteCoincidence = hasRemoteCoincidence || (remotePolarity == localPolarity);
      }
    }
  }

  return hasRemoteCoincidence;
}

bool SideSetHelper::element_side_has_coincidence_using_elem_elem_graph( Entity /*side*/,  Entity element,  ConnectivityOrdinal ordinal)
{
  bool hasLocalCoincidence  = element_side_has_local_coincidence_using_elem_elem_graph(element, ordinal);
  bool hasRemoteCoincidence = element_side_has_remote_coincidence_using_elem_elem_graph(element, ordinal);

  bool hasCoincidence = hasLocalCoincidence || hasRemoteCoincidence;
  return hasCoincidence;
}

bool SideSetHelper::element_side_has_coincidence_using_connectivity(const Entity side, const Entity element, const ConnectivityOrdinal ordinal)
{
  if(mesh.entity_rank(element) == stk::topology::ELEM_RANK) {

    unsigned numElements = mesh.num_elements(side);
    const Entity* elems = mesh.begin_elements(side);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(side);

    if(numElements == 1 && elems[0] == element && ordinals[0] == ordinal) {
      return true;
    }

    EntityVector inputSideNodes, sideNodes;
    impl::fill_element_side_nodes_from_topology(mesh.bucket(element).topology(), mesh.begin_nodes(element), ordinal, inputSideNodes);

    for(unsigned i=0; i<numElements; ++i) {
      if(elems[i] != element) {
        impl::fill_element_side_nodes_from_topology(mesh.bucket(elems[i]).topology(), mesh.begin_nodes(elems[i]), ordinals[i], sideNodes);

        if(impl::is_coincident_connection(mesh, element, inputSideNodes, ordinal, mesh.bucket(elems[i]).topology(), sideNodes, ordinals[i])) {
          return true;
        }
      }
    }
  }

  return false;
}

bool SideSetHelper::element_side_has_coincidence(SideSetSelector& sideset, const Entity side,
                                                 const Entity element, const ConnectivityOrdinal ordinal)
{
  if(mesh.bucket(element).owned() && mesh.has_face_adjacent_element_graph()) {
    return element_side_can_be_distinguished_using_elem_elem_graph(side, element, ordinal, sideset) ||
           element_side_has_coincidence_using_elem_elem_graph(side, element, ordinal);
  }

  return element_side_can_be_distinguished_using_connectivity(side, element, ordinal, sideset) ||
         element_side_has_coincidence_using_connectivity(side, element, ordinal);
}

void SideSetHelper::add_coincident_element_side_entry_to_sideset(SideSetSelector& sideset, const Entity side,
                                                                 const Entity element, const ConnectivityOrdinal ordinal)
{
  if (mesh.entity_rank(side) == mesh.mesh_meta_data().side_rank()) {
    if(element_side_has_coincidence(sideset, side, element, ordinal)) {
      add_element_side_entry_to_sideset(sideset, element,  ordinal);
    }
  }
}

void SideSetHelper::add_coincident_side_entries_to_sideset(SideSetSelector& sideset, const Entity entity)
{
  if (mesh.entity_rank(entity) == mesh.mesh_meta_data().side_rank()) {
    std::vector<SideSetEntry> entries;
    fill_coincident_sideset_entries_for_side(entity, entries);

    for(const SideSetEntry& entry : entries) {
      add_element_side_entry_to_sideset(sideset, entry.element,  entry.side);
    }
  }
}

void SideSetHelper::add_side_entries_to_sideset(SideSetSelector& sideset, const Entity entity)
{
  if (mesh.entity_rank(entity) == mesh.mesh_meta_data().side_rank()) {
    unsigned numElements = mesh.num_elements(entity);
    const Entity* elems = mesh.begin_elements(entity);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(entity);

    for(unsigned i=0; i<numElements; ++i) {
      add_element_side_entry_to_sideset(sideset, elems[i],  ordinals[i]);
    }
  }
}

void SideSetHelper::add_side_entries_to_sidesets(SideSetSelectorVector& sidesets, const Entity entity)
{
  if (sidesets.size() > 0 && mesh.entity_rank(entity) == mesh.mesh_meta_data().side_rank()) {
    unsigned numElements = mesh.num_elements(entity);
    const Entity* elems = mesh.begin_elements(entity);
    const ConnectivityOrdinal * ordinals = mesh.begin_element_ordinals(entity);

    for(unsigned i=0; i<numElements; ++i) {
      add_element_side_entry_to_sidesets(sidesets, elems[i],  ordinals[i]);
    }
  }
}

void SideSetHelper::add_sideset_entry_for_element_selected_by_sidesets(Entity entity, SideSetSelectorVector& sidesetsAndSelectors)
{
  if(mesh.bucket_ptr(entity) == nullptr) { return; }

  const unsigned numSides = stk::mesh::num_sides(mesh, entity);

  if(sidesetsAndSelectors.size() > 0 && mesh.entity_rank(entity) == stk::topology::ELEM_RANK && numSides > 0) {
    const std::vector<stk::mesh::ConnectivityOrdinal> ordinals = stk::mesh::get_side_ordinals(mesh, entity);
    const stk::mesh::EntityVector sides = stk::mesh::get_sides(mesh, entity);

    stk::mesh::SideSetEntry entry(entity);

    for(unsigned i=0; i<numSides; i++) {
      entry.side = ordinals[i];

     for(SideSetSelector& sidesetSelector : sidesetsAndSelectors) {
       const Selector& elemSelector = *sidesetSelector.selector();
       SideSet* sideset = sidesetSelector.sideset();

       if(elemSelector(mesh.bucket(entity)) && mesh.bucket(sides[i]).member(sidesetSelector.part())) {
         std::vector<stk::mesh::SideSetEntry>::iterator lowerBound = std::lower_bound(sideset->begin(), sideset->end(), entry);

         if(lowerBound == sideset->end() || *lowerBound != entry) {
           add_element_side_entry_to_sideset(sidesetSelector, entity, ordinals[i]);
         }
       }
     }
    }
  }
}

void SideSetHelper::warn_about_internal_sidesets()
{
  if(warnAboutInternalSidesets == false) { return; }

  bool hasInternalSidesets = internalSidesetOrdinals.size() > 0;

  if(mesh.parallel_size() > 1) {
    hasInternalSidesets =  stk::is_true_on_any_proc(mesh.parallel(), hasInternalSidesets);
  }

  if(!hasInternalSidesets) { return; }

  std::vector<int> localValues;

  std::vector<const Part *> surfacesInMap = mesh.mesh_meta_data().get_surfaces_in_surface_to_block_map();

  localValues.assign(surfacesInMap.size(), 0);

  for(size_t i = 0; i<surfacesInMap.size(); ++i) {
    const Part * surface = surfacesInMap[i];
    if(internalSidesetOrdinals.find(surface->mesh_meta_data_ordinal()) != internalSidesetOrdinals.end()) {
      localValues[i] = 1;
    }
  }

  std::vector<int> globalValues(localValues);
  if(mesh.parallel_size() > 1) {
    stk::all_reduce_max(mesh.parallel(), localValues.data(), globalValues.data(), globalValues.size());
  }

  for(size_t i=0; i<globalValues.size(); ++i) {
    if (!internalSidesetWarningHasBeenIssued && globalValues[i] != 0 && mesh.parallel_rank()==0 && outputStream != nullptr)
    {
      issue_internal_sideset_warning(surfacesInMap[i]->name(), *outputStream);
    }
    internalSidesetWarningHasBeenIssued = true;
  }
}

void SideSetHelper::reset_internal_sideset_detection(bool rankedPartsChangedOnElemsOrSides)
{
  internalSidesetOrdinals.clear();

  if(mesh.parallel_size() > 1 && mesh.has_face_adjacent_element_graph()) {
    const ElemElemGraph &graph = mesh.get_face_adjacent_element_graph();

    bool needToUpdateParallelPartInfo = (graph.mod_cycle_when_graph_modified() > m_modCycleWhenParallelPartInfoUpdated);
    bool globalRankedPartsChangedOnElemsOrSides = stk::is_true_on_any_proc(mesh.parallel(), rankedPartsChangedOnElemsOrSides);

    if (needToUpdateParallelPartInfo || globalRankedPartsChangedOnElemsOrSides) {
      parallelPartInfo.clear();

      impl::populate_part_ordinals_for_remote_edges(mesh, graph, parallelPartInfo);
      m_modCycleWhenParallelPartInfoUpdated = mesh.synchronized_count();
    }
  }
}

void SideSetHelper::warn_internal_sideset_detection()
{
  warn_about_internal_sidesets();

  internalSidesetOrdinals.clear();
}

void SideSetHelper::set_warn_about_internal_sideset(bool flag)
{
  warnAboutInternalSidesets = flag;
}

}
}

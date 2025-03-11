#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/base/SidesetUpdater.hpp"
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include "stk_mesh/baseImpl/SideSetUtilImpl.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/diag/StringUtil.hpp"           // for Type, etc
#include "stk_util/util/string_case_compare.hpp"
#include <stk_util/environment/RuntimeWarning.hpp>
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include <algorithm>
#include <utility>

namespace stk {
namespace mesh {
namespace impl {

PolarityDataCache::PolarityDataCache(const BulkData &bulk_, const Selector& activeSelector_, const ParallelPartInfo& parallelPartInfo_)
  : bulk(bulk_),
    meta(bulk.mesh_meta_data()),
    activeSelector(activeSelector_),
    parallelPartInfo(parallelPartInfo_),
    sideRank(meta.side_rank())
{

}

bool is_remote_permutation_positive(const ElemElemGraph& eeGraph, const GraphEdge& graphEdge,
                                    const Permutation localSidePermutation,
                                    const EntityKey localSideKey,
                                    const Permutation remoteSidePermutation,
                                    const EntityKey remoteSideKey)
{
  const impl::ParallelInfo &parallelInfo = eeGraph.get_parallel_info_for_graph_edge(graphEdge);
  const BulkData& bulk = eeGraph.get_mesh();
  const EntityRank sideRank = bulk.mesh_meta_data().side_rank();

  topology remoteSideTopology = parallelInfo.m_remote_element_topology.sub_topology(sideRank, graphEdge.side2());

  Entity localElement = eeGraph.get_entity(graphEdge.elem1());
  topology localSideTopology = bulk.bucket(localElement).topology().sub_topology(sideRank, graphEdge.side1());

  if(localSideTopology != remoteSideTopology) {
    stk::RuntimeWarningP0() <<
      "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" <<
      "!                                                                                             !\n" <<
      "! WARNING/ERROR: LOCAL SIDE TOPOLOGY <" << localSideTopology.name() <<
      "> DOES NOT MATCH REMOTE SIDE TOPLOGY <" << remoteSideTopology.name() << ">. !\n" <<
      "!                THIS SHOULD BE INTERPRETED AS AN ERROR.  CONTACT sierra-help@sandia.gov FOR  !\n" <<
      "!                A RESOLUTION BEFORE CONTINUING.                                              !\n" <<
      "!                                                                                             !\n" <<
      "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }

  if(localSideKey != remoteSideKey) {
    stk::RuntimeWarningP0() <<
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" <<
        "!                                                                                             !\n" <<
        "! WARNING/ERROR: LOCAL FACE: <" << localSideKey << "> for local element: <" <<
        eeGraph.get_mesh().identifier(localElement) << "> on side: <" << graphEdge.side1() <<
        "> does not match  !\n" <<
        "!                remote face: <" << remoteSideKey << "> for remote element: <" << -graphEdge.elem2() <<
        "> on side: <" << graphEdge.side2() << "> from P<" << parallelInfo.get_proc_rank_of_neighbor() <<
        ">  !\n" <<
        "!                THIS SHOULD BE INTERPRETED AS AN ERROR.  CONTACT sierra-help@sandia.gov FOR  !\n" <<
        "!                A RESOLUTION BEFORE CONTINUING.                                              !\n" <<
        "!                                                                                             !\n" <<
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }

  const bool  localSidePolarity =  localSideTopology.is_positive_polarity( localSidePermutation);

  if(remoteSidePermutation >= remoteSideTopology.num_permutations()) {
    return !localSidePolarity;
  }

  const bool remoteSidePolarity = remoteSideTopology.is_positive_polarity(remoteSidePermutation);

  if(localSidePolarity == remoteSidePolarity) {
    stk::RuntimeWarningP0() <<
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" <<
        "!                                                                                             !\n" <<
        "! WARNING/ERROR: NON-OPPOSING POLARITIES FOR FACE BETWEEN ELEMENTS " <<
        eeGraph.get_mesh().identifier(localElement) << " AND " << -graphEdge.elem2() << ".    !\n" <<
        "                 LOCAL PERMUTATION: " << (int)localSidePermutation <<
        " REMOTE PERMUTATION: " << (int)remoteSidePermutation <<
        "                               !\n" <<
        "!                THIS SHOULD BE INTERPRETED AS AN ERROR.  CONTACT sierra-help@sandia.gov FOR  !\n" <<
        "!                A RESOLUTION BEFORE CONTINUING.                                              !\n" <<
        "!                                                                                             !\n" <<
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  }

  return remoteSidePolarity;
}

void add_boundary_face_connectivity_info(const ElemElemGraph& eeGraph, const GraphEdge& graphEdge,
                                         const Entity face,
                                         const Permutation localSidePermutation, const SideSet& sset,
                                         PolarityDataCache& cache)
{
  const auto iter = cache.parallelPartInfo.find(graphEdge.elem2());
  bool isPhysicalBoundary = (iter == cache.parallelPartInfo.end());
  if(isPhysicalBoundary) {
    return;
  }

  const stk::mesh::impl::ParallelPartInfoValue& partInfoValue = iter->second;
  const std::vector<PartOrdinal>& otherElementPartOrdinals = partInfoValue.elementPartOrdinals;

  const stk::mesh::impl::ParallelSideInfoValue* sideInfoValue = partInfoValue.get_parallel_side_info(graphEdge);

  const std::vector<PartOrdinal>& otherSidesetPartOrdinals = sideInfoValue->sidesetPartOrdinals;
  const Permutation remoteSidePermutation = sideInfoValue->sidePermutation;
  const EntityKey remoteSideKey = sideInfoValue->sideKey;

  stk::mesh::EntityId otherElemId = eeGraph.convert_negative_local_id_to_global_id(graphEdge.elem2());
  Part* remoteElementBlockPart = get_element_block_part(cache.bulk, otherElementPartOrdinals, otherElemId);

  const bool isRemoteActive = cache.activeSelector(remoteElementBlockPart);
  const bool isRemoteInSideset = std::binary_search(otherSidesetPartOrdinals.begin(), otherSidesetPartOrdinals.end(),
                                                    sset.get_part()->mesh_meta_data_ordinal());
  const bool isRemotePermutationPositive = is_remote_permutation_positive(eeGraph, graphEdge,
                                                                          localSidePermutation, cache.bulk.entity_key(face),
                                                                          remoteSidePermutation, remoteSideKey);

  Entity localElement = eeGraph.get_entity(graphEdge.elem1());
  bool isLocalActive = cache.activeSelector(cache.bulk.bucket(localElement));
  bool isLocalInSideset = isLocalActive && sset.contains(localElement, graphEdge.side1());

  bool isInternalSidesetOnProcessorBoundary = (isLocalInSideset == isRemoteInSideset) && (isLocalActive == isRemoteActive);
  if(isInternalSidesetOnProcessorBoundary) {
    return;
  }

  cache.isConnectedElementActive.push_back(isRemoteActive);
  cache.isConnectedElementInSideset.push_back(isRemoteInSideset);
  cache.isConnectedElementPermutationPositive.push_back(isRemotePermutationPositive);
}

void insert_face_connectivity_info(const SideSet& sset,
                                   const Entity element,
                                   const ConnectivityOrdinal ordinal,
                                   const Permutation permutation,
                                   PolarityDataCache& cache)
{
  const Bucket& bucket = cache.bulk.bucket(element);
  const topology sideTopology = bucket.topology().sub_topology(cache.sideRank, ordinal);

  const bool isActive = cache.activeSelector(bucket);
  const bool isInSideset = isActive && sset.contains(element, ordinal);
  const bool isPermutationPositive = sideTopology.is_positive_polarity(permutation);

  cache.isConnectedElementActive.push_back(isActive);
  cache.isConnectedElementInSideset.push_back(isInSideset);
  cache.isConnectedElementPermutationPositive.push_back(isPermutationPositive);
}

void add_interior_face_connectivity_info(const ElemElemGraph& eeGraph, const GraphEdge& graphEdge,
                                         const SideSet& sset,
                                         const unsigned numSideElements, const Entity* sideElements,
                                         const Permutation * sidePermutations,
                                         PolarityDataCache& cache)
{
  Entity localElement2 = eeGraph.get_entity(graphEdge.elem2());
  ConnectivityOrdinal localOrdinal2 = graphEdge.side2();
  Permutation localPermutation2 = INVALID_PERMUTATION;

  for(unsigned j = 0; j < numSideElements; j++) {
    if(localElement2 == sideElements[j]) {
      localPermutation2 = sidePermutations[j];
      break;
    }
  }

  if(localPermutation2 == INVALID_PERMUTATION) {
    return;
  }

  STK_ThrowRequireMsg(localPermutation2 != INVALID_PERMUTATION,
                      "Could not find permutation for other element: " <<
                      cache.bulk.entity_key(localElement2) <<
                      " on graph edge: " << graphEdge);

  insert_face_connectivity_info(sset, localElement2, localOrdinal2, localPermutation2, cache);
}

void add_non_local_face_connectivity_info(const Entity /*face*/, const SideSet& sset,
                                          const unsigned numSideElements, const Entity* sideElements,
                                          const ConnectivityOrdinal * sideOrdinals,
                                          const Permutation * sidePermutations,
                                          PolarityDataCache& cache)
{
   for(unsigned i=0; i<numSideElements; ++i) {
     insert_face_connectivity_info(sset, sideElements[i], sideOrdinals[i], sidePermutations[i], cache);
   }
}

void fill_face_connectivity_info(const Entity face, const SideSet& sset, PolarityDataCache& cache)
{
  // Update with graph using first locally owned connected element as reference
  const unsigned numSideElements = cache.bulk.num_elements(face);
  const Entity* sideElements = cache.bulk.begin_elements(face);
  const ConnectivityOrdinal * sideOrdinals = cache.bulk.begin_element_ordinals(face);
  const Permutation * sidePermutations = cache.bulk.begin_element_permutations(face);

  Entity localElement;
  ConnectivityOrdinal localOrdinal = INVALID_CONNECTIVITY_ORDINAL;
  Permutation localPermutation = INVALID_PERMUTATION;

  for(unsigned i = 0; i < numSideElements; i++) {
    if(cache.bulk.bucket(sideElements[i]).owned()) {
      localElement = sideElements[i];
      localOrdinal = sideOrdinals[i];
      localPermutation = sidePermutations[i];
      break;
    }
  }

  if(INVALID_CONNECTIVITY_ORDINAL == localOrdinal) {
    // None of the connected elements are local .. cannot use graph
    add_non_local_face_connectivity_info(face, sset, numSideElements,
                                         sideElements, sideOrdinals, sidePermutations, cache);
    return;
  }

  const ElemElemGraph& eeGraph = cache.bulk.get_face_adjacent_element_graph();

  impl::LocalId elemLocalId = eeGraph.get_local_element_id(localElement);
  GraphEdgesForElement graphEdges = eeGraph.get_edges_for_element(elemLocalId);

  for(size_t i = 0; i < graphEdges.size(); ++i) {
    const GraphEdge& graphEdge = eeGraph.get_graph().get_edge_for_element(elemLocalId, i);

    if(localOrdinal == graphEdge.side1()) {
      // Add info based on the local element
      insert_face_connectivity_info(sset, localElement, localOrdinal, localPermutation, cache);

      bool isInteriorSide = impl::is_local_element(graphEdge.elem2());
      if(isInteriorSide) {
        // Add info based on the other local element on the graph edge
        add_interior_face_connectivity_info(eeGraph, graphEdge, sset, numSideElements, sideElements,
                                            sidePermutations, cache);
      } else {
        add_boundary_face_connectivity_info(eeGraph, graphEdge, face, localPermutation, sset, cache);
      }

      break;
    }
  }
}

const SideSet* get_sideset_pointer(const BulkData& bulk, const Part& sideSetPart, const SideSet* inputSidesetPtr)
{
  const SideSet* ssetPtr = inputSidesetPtr;
  if (ssetPtr == nullptr) {
    const Part &parentPart = get_sideset_parent(sideSetPart);
    if(bulk.does_sideset_exist(parentPart)) {
      ssetPtr = &bulk.get_sideset(parentPart);
    }
  }

  return ssetPtr;
}

std::pair<bool,bool> internal_is_positive_sideset_polarity(const Part& sideSetPart,
                                                           const Entity face,
                                                           PolarityDataCache& cache,
                                                           const SideSet* inputSidesetPtr)
{
  const MeshIndex& meshIndex = cache.bulk.mesh_index(face);
  const Bucket& bucket = *meshIndex.bucket;

  STK_ThrowRequire(bucket.entity_rank() == cache.bulk.mesh_meta_data().side_rank());
  STK_ThrowAssert(bucket.member(sideSetPart));

  const SideSet* ssetPtr = get_sideset_pointer(cache.bulk, sideSetPart, inputSidesetPtr);
  if (ssetPtr == nullptr) {
    return std::make_pair(false,false);
  }

  const SideSet& sset = *ssetPtr;

  cache.isConnectedElementActive.clear();
  cache.isConnectedElementInSideset.clear();
  cache.isConnectedElementPermutationPositive.clear();

  fill_face_connectivity_info(face, sset, cache);

  const unsigned numSideElements = cache.isConnectedElementActive.size();

  if (numSideElements == 1) {
    //hopefully the most common case, fast/early return:
    const bool isActive = cache.isConnectedElementActive[0];
    const bool inSideset = isActive && cache.isConnectedElementInSideset[0];
    const bool positivePermutation = cache.isConnectedElementPermutationPositive[0];
    const bool positivePolarity = (inSideset) ? positivePermutation : !positivePermutation;
    return std::make_pair(inSideset, positivePolarity);
  }

  int numFound = 0;
  bool foundElemWithPermutationPositive = false;

  for(unsigned i=0; i<numSideElements; ++i) {
    const bool isActive = cache.isConnectedElementActive[i];
    const bool inSideset = isActive && cache.isConnectedElementInSideset[i];

    if (!isActive) { continue; }
    if (inSideset) {
      numFound++;

      if (cache.isConnectedElementPermutationPositive[i]) {
        foundElemWithPermutationPositive = true;
      }
    }
  }

  const bool status = (numFound > 0);
  const bool hasPositivePolarity = (numFound == 1 ? foundElemWithPermutationPositive : false);
  auto returnValue = std::make_pair(status, hasPositivePolarity);

  return returnValue;
}


std::pair<bool,bool> internal_is_positive_sideset_face_polarity(const Entity face, PolarityDataCache& cache)
{
  std::pair<bool,bool> returnValue(false,false);

  std::vector<const stk::mesh::Part*> sidesetParts = get_sideset_io_parts(cache.bulk, face);

  if(sidesetParts.size() == 0) {
    return returnValue;
  }

  auto it = std::find_if(sidesetParts.begin(), sidesetParts.end(),
                         [&](const stk::mesh::Part* p){
                           returnValue = internal_is_positive_sideset_polarity(*p, face, cache, nullptr);
                           return returnValue.first;
                         });

  if (sidesetParts.size() > 1) {
    for(auto& ssetIter = it; ssetIter<sidesetParts.end(); ++ssetIter)
    {
      const stk::mesh::Part& sidesetPart = *(*ssetIter);
      std::pair<bool,bool> partPolarity = internal_is_positive_sideset_polarity(sidesetPart, face, cache, nullptr);

      if (partPolarity.first && partPolarity != returnValue) {
        stk::RuntimeWarningP0() << "Polarity for face: " << cache.bulk.identifier(face) << " on sideset: "
                                << (*it)->name() << " has polarity = " << returnValue.first << ", " << returnValue.second
                                << " which does not match that on sideset: " << sidesetPart.name()
                                << " and has polarity = " << partPolarity.first << ", " << partPolarity.second;
      } else if (returnValue.first == false && partPolarity.first == true) {
        returnValue = partPolarity;
      }
    }
  }

  return returnValue;
}

}}}


#include <algorithm>

#include "ParallelInfoForGraph.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/SideSetUtil.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/EquivalentEntityBlocks.hpp"
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{
void pack_data_for_part_ordinals(stk::CommSparse &comm, const ElemElemGraph& graph, const stk::mesh::BulkData& bulkData);
void pack_edge(stk::CommSparse &comm, const ElemElemGraph& graph, const stk::mesh::BulkData& bulkData, const stk::mesh::GraphEdge& edge, int other_proc);
void unpack_and_update_part_ordinals(stk::CommSparse &comm, const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph);
void unpack_and_update_part_ordinals(stk::CommSparse &comm, const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, ParallelPartInfo &parallelPartInfo);
stk::mesh::GraphEdge unpack_edge(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, int proc_id);

const stk::mesh::impl::ParallelSideInfoValue*
ParallelPartInfoValue::get_parallel_side_info(const GraphEdge& edge) const
{
  auto iter = std::lower_bound(sideInfo.begin(), sideInfo.end(), edge.side2());
  if((iter != sideInfo.end()) && (iter->sideOrdinal == edge.side2())) {
    return &(*iter);
  }

  return nullptr;
}

stk::mesh::impl::ParallelSideInfoValue*
ParallelPartInfoValue::get_parallel_side_info(const GraphEdge& edge)
{
  auto iter = std::lower_bound(sideInfo.begin(), sideInfo.end(), edge.side2());
  if((iter != sideInfo.end()) && (iter->sideOrdinal == edge.side2())) {
    return &(*iter);
  }

  return nullptr;
}

void populate_part_ordinals_for_remote_edges(const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, ParallelPartInfo &parallelPartInfo)
{
    parallelPartInfo.clear();
    stk::CommSparse comm(bulkData.parallel());
    pack_data_for_part_ordinals(comm, graph, bulkData);
    comm.allocate_buffers();
    pack_data_for_part_ordinals(comm, graph, bulkData);
    comm.communicate();
    unpack_and_update_part_ordinals(comm, bulkData, graph, parallelPartInfo);
}

std::pair<stk::mesh::Entity, stk::mesh::Permutation>
get_face_and_permutation_on_element_side_ordinal(const stk::mesh::BulkData& bulkData,
                                                 const stk::mesh::Entity localElement,
                                                 const stk::mesh::ConnectivityOrdinal localOrdinal)
{
  stk::mesh::Entity face;
  stk::mesh::Permutation permutation = INVALID_PERMUTATION;

  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  const stk::mesh::EntityRank sideRank = metaData.side_rank();

  const stk::mesh::MeshIndex& meshIndex = bulkData.mesh_index(localElement);
  const stk::mesh::Bucket& bucket = *meshIndex.bucket;

  const unsigned bucketOrdinal = meshIndex.bucket_ordinal;

  const unsigned numSides = bucket.num_connectivity(bucketOrdinal, sideRank);
  const stk::mesh::Entity* sides = bucket.begin(bucketOrdinal, sideRank);
  const stk::mesh::ConnectivityOrdinal * sideOrdinals = bucket.begin_ordinals(bucketOrdinal, sideRank);
  const stk::mesh::Permutation * sidePermutations= bucket.begin_permutations(bucketOrdinal, sideRank);

  for(unsigned i=0; i<numSides; i++) {
    if(localOrdinal == sideOrdinals[i]) {
      face = sides[i];
      permutation = sidePermutations[i];
      break;
    }
  }

  return std::make_pair(face, permutation);
}

void get_sideset_part_ordinals_and_permutation(const stk::mesh::Entity localElement,
                                               const stk::mesh::ConnectivityOrdinal localOrdinal,
                                               const stk::mesh::BulkData& bulkData,
                                               std::vector<stk::mesh::PartOrdinal>& sidesetPartOrdinals,
                                               stk::mesh::Permutation& sidePermutation,
                                               stk::mesh::EntityKey& sideKey)
{
  sidesetPartOrdinals.clear();
  const auto faceAndPermutation = get_face_and_permutation_on_element_side_ordinal(bulkData, localElement, localOrdinal);
  stk::mesh::Entity face = faceAndPermutation.first;
  sidePermutation = faceAndPermutation.second;
  sideKey = stk::mesh::EntityKey();
  if(bulkData.is_valid(face)) {
    sideKey = bulkData.entity_key(face);
    std::vector<const stk::mesh::Part*> sidesetParts = stk::mesh::get_sideset_io_parts(bulkData, face);
    for(const stk::mesh::Part *sidesetPart : sidesetParts) {
      const stk::mesh::Part &parentPart = stk::mesh::get_sideset_parent(*sidesetPart);
      if(bulkData.does_sideset_exist(parentPart)) {
        const stk::mesh::SideSet& sset = bulkData.get_sideset(parentPart);

        if(sset.contains(localElement, localOrdinal)) {
          sidesetPartOrdinals.push_back(parentPart.mesh_meta_data_ordinal());
        }
      }
    }

    stk::util::sort_and_unique(sidesetPartOrdinals);
  }
}

void pack_remote_side_info(stk::CommSparse &comm,
                           const stk::mesh::BulkData& bulkData,
                           const ElemElemGraph& graph, int proc,
                           const stk::mesh::GraphEdge& edge)
{
  std::vector<stk::mesh::PartOrdinal> sidesetPartOrdinals;
  stk::mesh::Entity localElement = graph.get_entity(edge.elem1());
  stk::mesh::ConnectivityOrdinal localOrdinal = edge.side1();

  stk::mesh::Permutation sidePermutation;
  stk::mesh::EntityKey sideKey;
  get_sideset_part_ordinals_and_permutation(localElement, localOrdinal, bulkData,
                                            sidesetPartOrdinals, sidePermutation, sideKey);

  bool hasFace = true;
  if(sideKey == stk::mesh::EntityKey()) {
    hasFace = false;
    comm.send_buffer(proc).pack<bool>(hasFace);
    return;
  }

  comm.send_buffer(proc).pack<bool>(hasFace);
  unsigned localSidePermutation = sidePermutation;
  comm.send_buffer(proc).pack<unsigned>(localSidePermutation);
  comm.send_buffer(proc).pack<stk::mesh::EntityKey>(sideKey);
  comm.send_buffer(proc).pack<stk::mesh::PartOrdinal>(sidesetPartOrdinals);
}

void pack_data_for_part_ordinals(stk::CommSparse &comm, const ElemElemGraph& graph, const stk::mesh::BulkData& bulkData)
{
    const stk::mesh::impl::ParallelGraphInfo& parallel_info = graph.get_parallel_graph().get_parallel_graph_info();
    std::vector<stk::mesh::PartOrdinal> elemPartOrdinals;
    const stk::mesh::Bucket* lastBucketPtr = nullptr;
    for(const auto& item : parallel_info)
    {
        const stk::mesh::GraphEdge &edge = item.first;
        const stk::mesh::impl::ParallelInfo &pinfo = item.second;
        stk::mesh::Entity localElement = graph.get_entity(edge.elem1());
        STK_ThrowRequireWithSierraHelpMsg(bulkData.is_valid(localElement));
        const stk::mesh::Bucket* currentBucketPtr = bulkData.bucket_ptr(localElement);
        if (lastBucketPtr != currentBucketPtr) {
          lastBucketPtr = currentBucketPtr;
          stk::mesh::impl::get_element_block_part_ordinals(localElement, bulkData, elemPartOrdinals);
        }

        int other_proc = pinfo.get_proc_rank_of_neighbor();
        pack_edge(comm, graph, bulkData, edge, other_proc);
        comm.send_buffer(other_proc).pack<stk::mesh::PartOrdinal>(elemPartOrdinals);

        pack_remote_side_info(comm, bulkData, graph, other_proc, edge);
    }
}

void pack_edge(stk::CommSparse &comm, const ElemElemGraph& graph, const stk::mesh::BulkData& bulkData, const stk::mesh::GraphEdge& edge, int other_proc)
{
    stk::mesh::EntityId id1 = bulkData.identifier(graph.get_entity(edge.elem1()));
    unsigned side1 = edge.side1();
    stk::mesh::EntityId id2 = -edge.elem2();
    unsigned side2 = edge.side2();
    comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(id1);
    comm.send_buffer(other_proc).pack<unsigned>(side1);
    comm.send_buffer(other_proc).pack<stk::mesh::EntityId>(id2);
    comm.send_buffer(other_proc).pack<unsigned>(side2);
}

void unpack_remote_side_info_value(stk::CommSparse &comm,
                                   const stk::mesh::BulkData& bulkData,
                                   const ElemElemGraph& graph, int proc,
                                   const stk::mesh::GraphEdge& edge,
                                   stk::mesh::impl::ParallelSideInfoValue& sideInfoValue)
{
  sideInfoValue.sideOrdinal = edge.side2();

  bool hasFace;
  comm.recv_buffer(proc).unpack<bool>(hasFace);

  if(!hasFace) {
    sideInfoValue.sidePermutation = INVALID_PERMUTATION;
    sideInfoValue.sideKey = stk::mesh::EntityKey();
    sideInfoValue.sidesetPartOrdinals.clear();
    return;
  }

  unsigned sidePermutation;
  comm.recv_buffer(proc).unpack<unsigned>(sidePermutation);
  sideInfoValue.sidePermutation = static_cast<stk::mesh::Permutation>(sidePermutation);
  comm.recv_buffer(proc).unpack<stk::mesh::EntityKey>(sideInfoValue.sideKey);

  std::vector<stk::mesh::PartOrdinal>& sidesetPartOrdinals = sideInfoValue.sidesetPartOrdinals;
  comm.recv_buffer(proc).unpack<stk::mesh::PartOrdinal>(sidesetPartOrdinals);
}

void unpack_remote_side_info(stk::CommSparse &comm,
                             const stk::mesh::BulkData& bulkData,
                             const ElemElemGraph& graph, int proc,
                             const stk::mesh::GraphEdge& edge,
                             std::vector<stk::mesh::impl::ParallelSideInfoValue>& sideInfo)
{
  auto iter = std::lower_bound(sideInfo.begin(), sideInfo.end(), edge.side2());
  if((iter != sideInfo.end()) && (iter->sideOrdinal == edge.side2())) {
    unpack_remote_side_info_value(comm, bulkData, graph, proc, edge, *iter);
  } else {
    stk::mesh::impl::ParallelSideInfoValue sideInfoValue;
    unpack_remote_side_info_value(comm, bulkData, graph, proc, edge, sideInfoValue);
    stk::util::insert_keep_sorted_and_unique(sideInfoValue, sideInfo);
  }
}

void unpack_and_update_part_ordinals(stk::CommSparse &comm, const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, ParallelPartInfo &parallelPartInfo)
{
    for(int i=0;i<bulkData.parallel_size();++i)
    {
        while(comm.recv_buffer(i).remaining())
        {
            stk::mesh::GraphEdge edge = unpack_edge(comm, bulkData, graph, i);

            stk::mesh::impl::ParallelPartInfoValue& infoValue = parallelPartInfo[edge.elem2()];
            std::vector<stk::mesh::PartOrdinal>& elemPartOrdinals = infoValue.elementPartOrdinals;
            comm.recv_buffer(i).unpack<stk::mesh::PartOrdinal>(elemPartOrdinals);

            unpack_remote_side_info(comm, bulkData, graph, i, edge, infoValue.sideInfo);
        }
    }
}

stk::mesh::GraphEdge unpack_edge(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, const ElemElemGraph& graph, int proc_id)
{
    stk::mesh::EntityId id1 = 0, id2 = 0;
    unsigned side1 = 0, side2 = 0;
    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(id1);
    comm.recv_buffer(proc_id).unpack<unsigned>(side1);
    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(id2);
    comm.recv_buffer(proc_id).unpack<unsigned>(side2);

    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, id2);
    STK_ThrowRequireWithSierraHelpMsg(bulkData.is_valid(element));

    stk::mesh::impl::LocalId localId2 = graph.get_local_element_id(element);
    stk::mesh::GraphEdge edge(localId2, side2, -id1, side1);
    return edge;
}

void pack_selected_value_for_par_info(stk::CommSparse &comm, int procRank, const stk::mesh::BulkData& bulkData, stk::mesh::Entity local_element, const stk::mesh::Selector& sel)
{
    if(sel(bulkData.bucket(local_element)))
        comm.send_buffer(procRank).pack<int64_t>(bulkData.identifier(local_element));
}

void pack_data_for_selector(stk::CommSparse &comm, const ElemElemGraph& graph, const stk::mesh::BulkData& bulkData, const stk::mesh::Selector& sel)
{
    for(const auto& item : graph.get_parallel_graph().get_parallel_graph_info())
    {
        const stk::mesh::GraphEdge &edge = item.first;
        const stk::mesh::impl::ParallelInfo &parInfo = item.second;
        pack_selected_value_for_par_info(comm, parInfo.get_proc_rank_of_neighbor(), bulkData, graph.get_entity(edge.elem1()), sel);
    }
}

void pack_and_communicate_selector(const stk::mesh::BulkData& bulkData,
                                   stk::CommSparse& comm,
                                   const ElemElemGraph& graph,
                                   const stk::mesh::Selector& selector)
{
    stk::pack_and_communicate(comm, [&comm, &graph, &bulkData, &selector]()
    {
        pack_data_for_selector(comm, graph, bulkData, selector);
    });
}

class RemoteSelectedValue
{
public:
    void set_id_as_selected(int64_t id)
    {
        idsSelected.insert(id);
    }
    bool is_id_selected(int64_t id) const
    {
        return idsSelected.find(id) != idsSelected.end();
    }
private:
    std::set<int64_t> idsSelected;
};

void unpack_selector_value(stk::CommSparse& comm, int rank, RemoteSelectedValue &remoteSelectedValue)
{
    int64_t id;
    comm.recv_buffer(rank).unpack<int64_t>(id);
    remoteSelectedValue.set_id_as_selected(id);
}

void unpack_selected_states(const stk::mesh::BulkData& bulkData,
                            stk::CommSparse& comm,
                            RemoteSelectedValue &remoteSelectedValue)
{
    stk::unpack_communications(comm, [&comm, &remoteSelectedValue](int rank)
    {
        unpack_selector_value(comm, rank, remoteSelectedValue);
    });
}

void update_selected_values(const ElemElemGraph& graph,
                            const RemoteSelectedValue &remoteSelectedValue,
                            ParallelSelectedInfo &selInfo)
{
    const stk::mesh::impl::ParallelGraphInfo& parallel_info = graph.get_parallel_graph().get_parallel_graph_info();
    for(const stk::mesh::impl::ParallelGraphInfo::value_type& edgeAndParInfo : parallel_info)
    {
        auto graphEdge_elem2 = edgeAndParInfo.first.elem2();
        selInfo[graphEdge_elem2] = remoteSelectedValue.is_id_selected(-graphEdge_elem2);
    }
}

void unpack_and_update_selector_value(stk::CommSparse &comm,
                                      const stk::mesh::BulkData& bulkData,
                                      const ElemElemGraph& graph,
                                      ParallelSelectedInfo &selInfo)
{
    RemoteSelectedValue remoteSelectedValue;
    unpack_selected_states(bulkData, comm, remoteSelectedValue);
    update_selected_values(graph, remoteSelectedValue, selInfo);
}

void populate_selected_value_for_remote_elements(const stk::mesh::BulkData& bulkData,
                                                 const ElemElemGraph& graph,
                                                 const stk::mesh::Selector& selector,
                                                 ParallelSelectedInfo &selInfo)
{
    selInfo.clear();
    stk::CommSparse comm(bulkData.parallel());
    pack_and_communicate_selector(bulkData, comm, graph, selector);
    unpack_and_update_selector_value(comm, bulkData, graph, selInfo);
}

}
}
}

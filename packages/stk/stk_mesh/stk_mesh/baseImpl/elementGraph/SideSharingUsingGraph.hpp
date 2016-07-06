#ifndef SIDESHARINGUSINGGRAPH_HPP_
#define SIDESHARINGUSINGGRAPH_HPP_
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>
#include <vector>
#include <stk_util/parallel/CommSparse.hpp>

namespace stk { namespace mesh { class BulkData; } }

namespace stk { namespace mesh {

struct SideSharingData
{
    stk::mesh::impl::IdViaSidePair elementAndSide;
    stk::mesh::Entity side;
    int sharingProc;
    int owningProc;
    stk::mesh::EntityId chosenSideId;
    std::vector<stk::mesh::EntityId> sideNodes;
    stk::mesh::OrdinalVector partOrdinals;
    SideSharingData() : elementAndSide({0,-1}), side(stk::mesh::Entity()), sharingProc(-1), owningProc(-1), chosenSideId(stk::mesh::InvalidEntityId) {}
    SideSharingData(const stk::mesh::impl::IdViaSidePair& sidePair, stk::mesh::Entity sideIn, int sharing_proc, int owning_proc, stk::mesh::EntityId chosen_id)
    : elementAndSide(sidePair), side(sideIn), sharingProc(sharing_proc), owningProc(owning_proc), chosenSideId(chosen_id) {}
};

stk::mesh::EntityVector fill_shared_entities_that_need_fixing(const stk::mesh::BulkData& bulkData);

void fill_sharing_data(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph &graph, const stk::mesh::EntityVector& sidesThatNeedFixing, std::vector<SideSharingData>& sideSharingDataThisProc, std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides);

void allocate_and_send(stk::CommSparse& comm, const std::vector<SideSharingData>& sideSharingDataThisProc, const std::vector<stk::mesh::impl::IdViaSidePair>& idAndSides);

void unpack_data(stk::CommSparse& comm, int my_proc_id, int num_procs, std::vector<SideSharingData>& sideSharingDataThisProc);

}}



#endif /* SIDESHARINGUSINGGRAPH_HPP_ */

#include "SideIdPool.hpp"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>


namespace stk
{
namespace mesh
{

SideIdPool::SideIdPool(stk::mesh::BulkData &bulk) :
        mBulkData(bulk),
        mSuggestedIds(),
        mNumIdsUsed(0)
{
}

void SideIdPool::generate_initial_ids(unsigned numIdsNeeded)
{
    mNumIdsUsed = 0;
    mBulkData.generate_new_ids(mBulkData.mesh_meta_data().side_rank(), 2*numIdsNeeded, mSuggestedIds);
}

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
STK_DEPRECATED
stk::mesh::EntityId SideIdPool::get_available_id()
{
    STK_ThrowRequireWithSierraHelpMsg(mNumIdsUsed < mSuggestedIds.size());
    stk::mesh::EntityId availId = *(mSuggestedIds.begin() + mNumIdsUsed);
    mNumIdsUsed++;
    return availId;
}
#endif

void SideIdPool::generate_additional_ids_collective(size_t num_additional_ids_needed)
{
    size_t num_available = mSuggestedIds.size() - mNumIdsUsed;
    
    bool need_additional_ids = num_additional_ids_needed > num_available;
    if(stk::is_true_on_any_proc(mBulkData.parallel(), need_additional_ids))
    {
        size_t num_additional_to_request = num_additional_ids_needed + mSuggestedIds.size();
        stk::mesh::EntityIdVector new_ids;
        mBulkData.generate_new_ids_given_reserved_ids(mBulkData.mesh_meta_data().side_rank(),
                                                      num_additional_to_request,
                                                      mSuggestedIds,
                                                      new_ids);
        mSuggestedIds.insert(mSuggestedIds.end(), new_ids.begin(), new_ids.end());
    }
}

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
STK_DEPRECATED
void SideIdPool::reset_suggested_side_id_iter(size_t numIdsNotReallyUsed)
{
    mNumIdsUsed -= numIdsNotReallyUsed;
}
#endif

}
}

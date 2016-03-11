#include "SideIdPool.hpp"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

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
    mBulkData.generate_new_ids(mBulkData.mesh_meta_data().side_rank(), numIdsNeeded, mSuggestedIds);
}

stk::mesh::EntityId SideIdPool::get_available_id()
{
    ThrowRequireWithSierraHelpMsg(mNumIdsUsed < mSuggestedIds.size());
    stk::mesh::EntityId availId = *(mSuggestedIds.begin() + mNumIdsUsed);
    mNumIdsUsed++;
    return availId;
}

void SideIdPool::generate_additional_ids_collective(size_t num_additional_ids_needed)
{
    stk::mesh::EntityIdVector new_ids;
    mBulkData.generate_new_ids_given_reserved_ids(mBulkData.mesh_meta_data().side_rank(),
                                                  num_additional_ids_needed,
                                                  mSuggestedIds,
                                                  new_ids);
    mSuggestedIds.insert(mSuggestedIds.end(), new_ids.begin(), new_ids.end());
}

void SideIdPool::reset_suggested_side_id_iter(size_t numIdsNotReallyUsed)
{
    mNumIdsUsed -= numIdsNotReallyUsed;
}

}
}

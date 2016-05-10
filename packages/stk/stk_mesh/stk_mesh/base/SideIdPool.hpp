#ifndef STK_SIDE_ID_POOL_HPP
#define STK_SIDE_ID_POOL_HPP

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }

namespace stk
{
namespace mesh
{

class SideIdPool
{
public:
    SideIdPool(stk::mesh::BulkData &bulk);

    void generate_initial_ids(unsigned numIdsNeeded);

    stk::mesh::EntityId get_available_id();

    void generate_additional_ids_collective(size_t num_additional_ids_needed);

    void reset_suggested_side_id_iter(size_t numIdsNotReallyUsed);

    void clear()
    {
        mNumIdsUsed = 0;
        mNumIdsRequested = 0;
        mSuggestedIds.clear();
    }

protected:
    const stk::mesh::EntityIdVector &get_all_ids()
    {
        return mSuggestedIds;
    }

private:
    stk::mesh::BulkData &mBulkData;
    stk::mesh::EntityIdVector mSuggestedIds;
    size_t mNumIdsUsed;
    size_t mNumIdsRequested;
};

}}

#endif

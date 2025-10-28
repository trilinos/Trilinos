// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_SIDE_ID_POOL_HPP
#define STK_SIDE_ID_POOL_HPP

#include <stk_util/stk_config.h>
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

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
    STK_DEPRECATED stk::mesh::EntityId get_available_id();
#endif

    void generate_additional_ids_collective(size_t num_additional_ids_needed);

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
    STK_DEPRECATED void reset_suggested_side_id_iter(size_t numIdsNotReallyUsed);
#endif

    void clear()
    {
        mNumIdsUsed = 0;
        mSuggestedIds.clear();
    }

    const stk::mesh::EntityIdVector &get_all_ids() const
    {
        return mSuggestedIds;
    }

private:
    stk::mesh::BulkData &mBulkData;
    stk::mesh::EntityIdVector mSuggestedIds;
    size_t mNumIdsUsed;
};

}}

#endif

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
// 
#ifndef _stk_mesh_unit_test_utils_BucketTester_hpp_
#define _stk_mesh_unit_test_utils_BucketTester_hpp_

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk { namespace mesh { namespace unit_test {

class BucketTester : public stk::mesh::Bucket
{
public:

    virtual ~BucketTester()
    {
    }

    void my_change_connected_nodes(unsigned bucket_ordinal, Entity* new_nodes)
    {
        this->change_existing_connectivity(bucket_ordinal, new_nodes);
    }

    void my_change_existing_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity,
            ConnectivityOrdinal elem_connectivity_ordinal, Permutation permut)
    {
        this->change_existing_permutation_for_connected_element(bucket_ordinal_of_lower_ranked_entity,
                elem_connectivity_ordinal, permut);
    }

    void my_change_existing_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, ConnectivityOrdinal edge_connectivity_ordinal,
            Permutation permut)
    {
        this->change_existing_permutation_for_connected_edge(bucket_ordinal_of_higher_ranked_entity, edge_connectivity_ordinal, permut);
    }

    void my_change_existing_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, ConnectivityOrdinal face_connectivity_ordinal,
            Permutation permut)
    {
        this->change_existing_permutation_for_connected_face(bucket_ordinal_of_higher_ranked_entity, face_connectivity_ordinal, permut);
    }
};

} } } // namespace stk mesh unit_test

#endif

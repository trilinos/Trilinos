// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#ifndef _BucketTesterr_hpp_
#define _BucketTesterr_hpp_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>

namespace stk { namespace mesh { namespace unit_test {

class BucketTester : public stk::mesh::Bucket
{
public:

//    BucketTester(stk::mesh::Bucket& bucket)
//    {
//    }

    virtual ~BucketTester()
    {
    }

    void my_change_exisiting_connectivity(unsigned bucket_ordinal, stk::mesh::Entity* new_nodes)
    {
        this->change_existing_connectivity(bucket_ordinal, new_nodes);
    }

    void my_change_exisiting_permutation_for_connected_element(unsigned bucket_ordinal_of_lower_ranked_entity,
            unsigned elem_connectivity_ordinal, stk::mesh::Permutation permut)
    {
        this->change_existing_permutation_for_connected_element(bucket_ordinal_of_lower_ranked_entity,
                elem_connectivity_ordinal, permut);
    }

    void my_change_exisiting_permutation_for_connected_edge(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned edge_connectivity_ordinal,
            stk::mesh::Permutation permut)
    {
        this->change_existing_permutation_for_connected_edge(bucket_ordinal_of_higher_ranked_entity, edge_connectivity_ordinal, permut);
    }

    void my_change_exisiting_permutation_for_connected_face(unsigned bucket_ordinal_of_higher_ranked_entity, unsigned face_connectivity_ordinal,
            stk::mesh::Permutation permut)
    {
        this->change_existing_permutation_for_connected_face(bucket_ordinal_of_higher_ranked_entity, face_connectivity_ordinal, permut);
    }
};

} } } // namespace stk mesh unit_test

#endif

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

#include <stk_mesh/base/GetBuckets.hpp>
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for sort, unique
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <utility>                      // for pair
#include <vector>                       // for vector, etc

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

void copy_ids( std::vector<unsigned> & v , const PartVector & p )
{
  {
    const size_t n = p.size();
    v.resize( n );
    for ( size_t k = 0 ; k < n ; ++k ) {
      v[k] = p[k]->mesh_meta_data_ordinal();
    }
  }

  {
    std::vector<unsigned>::iterator i = v.begin() , j = v.end();
    std::sort( i , j );
    i = std::unique( i , j );
    v.erase( i , j );
  }
}

void get_involved_parts(
    const PartVector & union_parts,
    const Bucket & candidate,
    PartVector & involved_parts
    )
{
  involved_parts.clear();
  if (union_parts.empty()) {
    return;
  }

  // Used to convert part ordinals to part pointers:
  MetaData & meta_data = MetaData::get( * union_parts[0]);
  const PartVector & all_parts = meta_data.get_parts();

  const std::pair<const unsigned *,const unsigned *>
    bucket_part_begin_end_iterators = candidate.superset_part_ordinals(); // sorted and unique

  std::vector<unsigned> union_parts_ids;
  copy_ids( union_parts_ids , union_parts ); // sorted and unique
  std::vector<unsigned>::const_iterator union_part_id_it = union_parts_ids.begin();
  const unsigned * bucket_part_id_it = bucket_part_begin_end_iterators.first ;

  while ( union_part_id_it != union_parts_ids.end() &&
          bucket_part_id_it != bucket_part_begin_end_iterators.second )
  {
    if      ( *union_part_id_it  < *bucket_part_id_it ) {
      ++union_part_id_it ;
    }
    else if ( *bucket_part_id_it < *union_part_id_it )  {
      ++bucket_part_id_it ;
    }
    else {
      // Find every match:
      Part * const part = all_parts[ *union_part_id_it ];
      involved_parts.push_back( part );
      ++union_part_id_it;
      ++bucket_part_id_it;
    }
  }

}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk



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

#include <stk_mesh/base/GetEntities.hpp>
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for sort
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for  BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityLess.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Selector.hpp"   // for Selector


namespace stk {
namespace mesh {

//----------------------------------------------------------------------

void get_entities( const BulkData & mesh , EntityRank entity_rank ,
                   std::vector< Entity> & entities,
                   bool sortByGlobalId )
{
  const BucketVector & ks = mesh.buckets( entity_rank );
  entities.clear();

  size_t count = 0;

  const BucketVector::const_iterator ie = ks.end();
        BucketVector::const_iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) { count += (*ik)->size(); }

  entities.reserve(count);

  ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    const Bucket & k = **ik ;
    size_t n = k.size();
    for(size_t i = 0; i < n; ++i) {
      entities.push_back(k[i]);
    }
  }

  if (sortByGlobalId) {
    std::sort(entities.begin(), entities.end(), EntityLess(mesh));
  }
}

std::vector<Entity> get_entities(const BulkData & mesh, EntityRank entity_rank, bool sortByGlobalId)
{
  std::vector<Entity> entities;
  get_entities(mesh, entity_rank, entities, sortByGlobalId);
  return entities;
}

unsigned count_selected_entities(
  const Selector & selector ,
  const BucketVector & input_buckets )
{
  size_t count = 0;

  const BucketVector::const_iterator ie = input_buckets.end();
        BucketVector::const_iterator ik = input_buckets.begin();

  for ( ; ik != ie ; ++ik ) {
    const Bucket & k = ** ik ;
    if ( selector( k ) ) { count += k.size(); }
  }

  return count ;
}

unsigned count_entities( const BulkData& bulk,
                         const EntityRank rank,
                         const Selector & selector )
{
  size_t count = 0;

  const BucketVector& buckets = bulk.get_buckets(rank, selector);

  for(const Bucket* bptr : buckets) {
    count += bptr->size();
  }

  return count ;
}


void get_selected_entities( const Selector & selector ,
                            const BucketVector & input_buckets ,
                            std::vector< Entity> & entities ,
                            bool sortByGlobalId )
{
  size_t count = count_selected_entities(selector,input_buckets);

  entities.resize(count);

  const BucketVector::const_iterator ie = input_buckets.end();
        BucketVector::const_iterator ik = input_buckets.begin();

  for ( size_t j = 0 ; ik != ie ; ++ik ) {
    const Bucket & k = ** ik ;
    if ( selector( k ) ) {
      const size_t n = k.size();
      for ( size_t i = 0; i < n; ++i, ++j ) {
        entities[j] = k[i] ;
      }
    }
  }

  if (input_buckets.size() > 0 && sortByGlobalId) {
    std::sort(entities.begin(), entities.end(), EntityLess(input_buckets[0]->mesh()));
  }
}

void get_entities( const BulkData& bulk,
                   const EntityRank rank,
                   const Selector & selector ,
                   std::vector< Entity> & entities ,
                   bool sortByGlobalId )
{
  const BucketVector& buckets = bulk.get_buckets(rank, selector);

  size_t count = 0;
  for (const Bucket* bptr : buckets) {
    count += bptr->size();
  }

  entities.clear();
  entities.reserve(count);

  for (const Bucket* bptr : buckets) {
    entities.insert(entities.end(), bptr->begin(), bptr->end());
  }

  if (entities.size() > 0 && sortByGlobalId) {
    std::sort(entities.begin(), entities.end(), EntityLess(bulk));
  }
}

std::vector<Entity> get_entities(const BulkData& bulk,
                                 const EntityRank rank,
                                 const Selector & selector,
                                 bool sortByGlobalId)
{
  std::vector<Entity> entities;
  get_entities(bulk, rank, selector, entities, sortByGlobalId);
  return entities;
}

//----------------------------------------------------------------------

void count_entities(
  const Selector & selector ,
  const BulkData & mesh ,
  std::vector<size_t> & count )
{
  const size_t nranks = mesh.mesh_meta_data().entity_rank_count();

  count.resize( nranks );

  for ( size_t i = 0 ; i < nranks ; ++i ) {
    count[i] = 0 ;

    const BucketVector & ks = mesh.buckets( static_cast<EntityRank>(i) );

    BucketVector::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( selector(**ik) ) {
        count[i] += (*ik)->size();
      }
    }
  }
}

unsigned get_num_entities(const stk::mesh::BulkData &bulk)
{
    unsigned numEntities = 0;
    std::vector<size_t> countPerRank;
    stk::mesh::count_entities(bulk.mesh_meta_data().universal_part(), bulk, countPerRank);
    for(unsigned count : countPerRank)
        numEntities += count;
    return numEntities;
}

} // namespace mesh
} // namespace stk

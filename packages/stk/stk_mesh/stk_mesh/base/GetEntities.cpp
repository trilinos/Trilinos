/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/GetEntities.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

void get_entities( const BulkData & mesh , EntityRank type ,
                   std::vector< Entity*> & entities )
{
  const std::vector<Bucket*> & ks = mesh.buckets( type );
  entities.clear();

  size_t count = 0;

  const std::vector<Bucket*>::const_iterator ie = ks.end();
        std::vector<Bucket*>::const_iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) { count += (*ik)->size(); }

  entities.reserve(count);

  ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    const Bucket & k = **ik ;
    size_t n = k.size();
    for(size_t i = 0; i < n; ++i) {
      entities.push_back(&k[i]);
    }
  }

  std::sort(entities.begin(), entities.end(), EntityLess());
}

unsigned count_selected_entities(
  const Selector & selector ,
  const std::vector< Bucket * > & input_buckets )
{
  size_t count = 0;

  const std::vector<Bucket*>::const_iterator ie = input_buckets.end();
        std::vector<Bucket*>::const_iterator ik = input_buckets.begin();

  for ( ; ik != ie ; ++ik ) {
    const Bucket & k = ** ik ;
    if ( selector( k ) ) { count += k.size(); }
  }

  return count ;
}


void get_selected_entities( const Selector & selector ,
                            const std::vector< Bucket * > & input_buckets ,
                            std::vector< Entity * > & entities )
{
  size_t count = count_selected_entities(selector,input_buckets);

  entities.resize(count);

  const std::vector<Bucket*>::const_iterator ie = input_buckets.end();
        std::vector<Bucket*>::const_iterator ik = input_buckets.begin();

  for ( size_t j = 0 ; ik != ie ; ++ik ) {
    const Bucket & k = ** ik ;
    if ( selector( k ) ) {
      const size_t n = k.size();
      for ( size_t i = 0; i < n; ++i, ++j ) {
        entities[j] = &k[i] ;
      }
    }
  }

  std::sort(entities.begin(), entities.end(), EntityLess());
}


//----------------------------------------------------------------------

void count_entities(
  const Selector & selector ,
  const BulkData & mesh ,
  std::vector< EntityRank > & count )
{
  const size_t ntype = MetaData::get(mesh).entity_rank_count();

  count.resize( ntype );

  for ( size_t i = 0 ; i < ntype ; ++i ) {
    count[i] = 0 ;

    const std::vector<Bucket*> & ks = mesh.buckets( i );

    std::vector<Bucket*>::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( selector(**ik) ) {
        count[i] += (*ik)->size();
      }
    }
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk


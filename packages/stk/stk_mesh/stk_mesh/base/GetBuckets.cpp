/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//----------------------------------------------------------------------

#include <algorithm>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

void get_buckets( const Selector & selector ,
                  const std::vector< Bucket * > & input ,
                        std::vector< Bucket * > & output )
{
  output.clear();
  for ( std::vector< Bucket * >::const_iterator
        i = input.begin() ; i != input.end() ; ++i ) {
    Bucket * const b = *i ;
    if ( selector( *b ) ) { output.push_back( b ); }
  }
}

AllSelectedBucketsRange get_buckets( const Selector & selector, const BulkData& mesh )
{
  AllBucketsRange all_buckets = mesh.get_bucket_range();
  return get_selected_bucket_range(all_buckets, selector);
}

AllBucketsRange get_buckets( const BulkData& mesh )
{
  return mesh.get_bucket_range();
}

AllBucketsRange get_buckets( EntityRank entity_rank, const BulkData& mesh )
{
  return mesh.get_bucket_range(entity_rank);
}

AllSelectedBucketsRange get_buckets( const Selector & selector, const AllBucketsRange& range)
{
  return get_selected_bucket_range(range, selector);
}

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
  if (union_parts.size() == 0) {
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



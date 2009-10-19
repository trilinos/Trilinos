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

void get_buckets( const SelectorInterface & selector ,
                  const std::vector< Bucket * > & input ,
                        std::vector< Bucket * > & output )
{
  for ( std::vector< Bucket * >::const_iterator
        i = input.begin() ; i != input.end() ; ++i ) {
    Bucket * const b = *i ;
    if ( selector( *b ) ) { output.push_back( b ); }
  }
}

void get_buckets( const SelectorInterface & selector ,
                  const std::vector< Bucket * > & input ,
                        std::vector< BucketAndParts > & output )
{
  for ( std::vector< Bucket * >::const_iterator
        i = input.begin() ; i != input.end() ; ++i ) {
    BucketAndParts tmp( *i );
    if ( selector( *tmp.bucket , tmp.parts ) ) { output.push_back( tmp ); }
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk



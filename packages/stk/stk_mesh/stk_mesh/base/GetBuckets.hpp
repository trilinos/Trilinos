#ifndef stk_mesh_GetBucket_hpp
#define stk_mesh_GetBucket_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------

/* \brief  Select buckets from the input to the output. */
void get_buckets( const SelectorInterface & selector ,
                  const std::vector< Bucket * > & input ,
                        std::vector< Bucket * > & output );

/* \brief  Select buckets and parts from the input to the output. */
void get_buckets( const SelectorInterface       & selector ,
                  const std::vector< Bucket * >       & input ,
                        std::vector< BucketAndParts > & output );

//----------------------------------------------------------------------
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif


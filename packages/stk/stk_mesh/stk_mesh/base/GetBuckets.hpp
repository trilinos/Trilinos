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
void get_buckets( const Selector & selector ,
                  const std::vector< Bucket * > & input ,
                        std::vector< Bucket * > & output );


/* \brief  Get the parts from the union part vector that the bucket is
 *         contained in. 
 */
void get_involved_parts( 
    const PartVector & union_parts,
    const Bucket & candidate,
    PartVector & involved_parts
    );

//----------------------------------------------------------------------
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif


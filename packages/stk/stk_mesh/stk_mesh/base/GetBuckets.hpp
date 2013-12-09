/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_GetBucket_hpp
#define stk_mesh_GetBucket_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Iterators.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------

/* \brief  Select buckets from the input to the output. Buckets in the input
 *         vector will be placed in the output vector if the bucket is
 *         selected by the selector argument.
 *         On entry, the output vector is cleared before being filled with
 *         selected buckets.
 */
template<class SELECTOR>
inline
void get_buckets( const SELECTOR & selector ,
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

/* \brief  Get the parts from the union part vector that the bucket is
 *         contained in.
 */
void get_involved_parts( const PartVector & union_parts,
                         const Bucket & candidate,
                         PartVector & involved_parts);

//----------------------------------------------------------------------
/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif


/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Comm_hpp
#define stk_mesh_Comm_hpp

#include <vector>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk_classic {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_bulk_data_parallel
 *  \{
 */

/** \brief Global counts for a mesh's entities. */
bool comm_mesh_counts( BulkData & ,
                       std::vector<size_t> & counts ,
                       bool = false );

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk_classic

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif


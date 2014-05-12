/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Comm_hpp
#define stk_mesh_Comm_hpp

#include <stddef.h>                     // for size_t
#include <vector>                       // for vector
namespace stk { namespace mesh { class BulkData; } }

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

class Selector;

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_bulk_data_parallel
 *  \{
 */

/** \brief Global counts for a mesh's entities. */
void comm_mesh_counts( const BulkData &bulk_data ,
                       std::vector<size_t> & counts ,
                       const Selector *select = 0);

void comm_mesh_counts( const BulkData & M ,
                       std::vector<size_t> & counts,
                       std::vector<size_t> & min_counts,
                       std::vector<size_t> & max_counts,
                       const Selector *selector = 0);

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif


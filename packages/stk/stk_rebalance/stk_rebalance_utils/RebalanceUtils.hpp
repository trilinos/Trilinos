/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

//----------------------------------------------------------------------

namespace stk {
  namespace rebalance {

    bool verify_dependent_ownership(const stk::mesh::EntityRank & parent_rank,
                                    stk::mesh::EntityVector & entities, 
                                    stk::mesh::DefaultFEM & fem);

  } // namepsace rebalance
} // namepsace stk

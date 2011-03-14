/*--------------------------------------------------------------------*/
/*    Copyright 2010 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001,2002 Sandia Corporation, Albuquerque, NM.

#ifndef stk_rebalance_Rebalance_hpp
#define stk_rebalance_Rebalance_hpp

#include <string>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_rebalance/Partition.hpp>

/** \file Rebalance.hpp
 *   \ingroup stk_rebalance_module
 *
 * \brief Static functions for dynamic load balancing.
 *
 *  The rebalance namespace is intended to provide an application
 *  the top level functions to perform a mesh redistribution.
 *  The action is controlled by instances of the Partition 
 *  class that is passed into the rebalance function.
 */

namespace stk {
namespace rebalance {

/** \brief Rebalance with a Partition object.
 *
 * \param bulk_data      BulkData must be in a parallel consistent state.
 *
 * \param selector       Used to select a subset of mesh entities to compute measure.
 *
 * \param coord_ref      The field containing the nodal coordinates. For the default
 *                       ZoltanPartition class in stk::reblance, this should be non-NULL.
 *
 * \param elem_weight_ref This field will be used by the \a Partition class and 
 *                        can be NULL.
 *
 * \param Partition       The base class of a derived class that is used to 
 *                        determine the new partition.  See the \a ZoltanPartition
 *                        class for an example.
 *
 * \param rank            Rank of the entities \a elem_weight_ref is defined on.
 *
 * This \a rebalance function will use the \a Partition object passed
 * to perform the rebalancing.  It will be necessary to use one of the 
 * pre-defined derived classes in stk::rebalance, like \a ZoltanPartition,
 * or to define your own.
 */
bool rebalance(mesh::BulkData & bulk_data ,
               const mesh::Selector & selector ,
               const VectorField * coord_ref ,
               const ScalarField * elem_weight_ref,
               Partition & partition,
               const stk::mesh::EntityRank rank = stk::mesh::InvalidEntityRank);

/** \} */

}
} // namespace stk

#endif

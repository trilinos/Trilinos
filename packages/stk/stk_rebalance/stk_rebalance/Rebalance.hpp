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

/** @file Rebalance.h
 * @brief Static functions for dynamic load balancing.
 *
 *  The rebalance namespace is intended to provide an application
 *  the top level functions to perform a mesh redistribution.
 *  The action is controlled by instances of the Partition 
 *  class that is passed into the rebalance function.
 */

namespace stk {
namespace rebalance {

/** \brief Determine if rebalancing is needed.
 *
 * \param bulk_data      BulkData must be in a parallel consistent state.
 *
 * \param load_measure   Field defined on mesh objects of rank "rank". 
 *                       Can be a NULL pointer.
 *
 * \param imbalance_threshold  rebalance needed if MAX divided by average load
 *                             measure exceeds this value.
 *
 * \param rank                 rank of mesh entities to define load measure.
 *
 * \param selector             used to select a subset of mesh objects to compute measure.
 *
 * This function calculates the total weight of the load on each processor by summing
 * the load_measure field over the selector objects of rank rank.  If selector is not
 * specified, all ejects of rank are summed.  If load_balance is not specified, it is 
 * assumed to be 1 for each object and the weight per processor is just the number
 * of objects on each processor.  After a processor weight is defined the MAX over the
 * processing grid is divided by the average to get a global imbalance which is
 * compared to imbalance_threshold.  True is returned if the global imbalance is
 * greater than imbalance_threshold.
 */

bool rebalance_needed(mesh::BulkData            & bulk_data,
                      const mesh::Field<double> * load_measure,
                      double                    & imbalance_threshold,
                      const stk::mesh::EntityRank rank = stk::mesh::InvalidEntityRank,
                      const mesh::Selector      * selector=NULL);

/** \brief Rebalance with a Partition object.
 * \param bulk_data      BulkData must be in a parallel consistent state.
 *
 * \param selector       Used to select a subset of mesh objects to compute measure.
 *
 * \param coord_ref      The field containing the nodal coordinates. For the default
 *                       ZoltanPartition class in stk::reblance, this should be non-NULL.
 *
 * \param elem_weight_ref This field will be used by the Partition class and is the  
 *                        same as that used by the rebalance_needed function.
 *                        can be NULL.
 *
 * \param Partition       The base class of a derived class that is used to 
 *                        determine the new partition.  See the ZoltanPartition
 *                        class for an example.
 *
 * \param rank            Rank of the entities elem_weight_ref is defined on.
 *
 * This rebalance function will use the Partition object passed
 * to perform the rebalancing.
 */
bool rebalance(mesh::BulkData & bulk_data ,
               const mesh::Selector & selector ,
               const VectorField * coord_ref ,
               const ScalarField * elem_weight_ref,
               Partition & partition,
               const stk::mesh::EntityRank rank = stk::mesh::InvalidEntityRank);

}
} // namespace stk

#endif

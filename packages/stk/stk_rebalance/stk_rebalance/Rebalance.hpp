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

/*---------------------------------------------------------------*/
//: Rebalance class is intended to provide an application
//: interface, API, for Engineering Science application codes to
//: permit dynamic load balancing (e.g., using Zoltan). The class
//: has no constructors, and only contains a handful of static
//: method functions.
/*---------------------------------------------------------------*/

namespace stk {
namespace rebalance {

/** Determine if rebalancing is needed.
 */

// comm can also be obtained from Bulk data (as parallel_machine).
bool rebalance_needed(mesh::BulkData &    bulk_data,
                      const mesh::Field<double> * load_measure,
                      ParallelMachine    comm,
                      double & imbalance_threshold);

/** Rebalance with a Partition object.
 * This rebalance function will use the Partition object passed
 * to perform the rebalancing.
 */
bool rebalance(mesh::BulkData & bulk_data ,
               const mesh::Selector & selector ,
               const VectorField * coord_ref ,
               const ScalarField * elem_weight_ref,
               Partition & partition);

}
} // namespace stk

#endif

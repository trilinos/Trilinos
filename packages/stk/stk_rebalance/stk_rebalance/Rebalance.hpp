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
bool rebalance_needed(mesh::BulkData            & rebal_bulk_data ,
                      const Real            imbalance_threshold);

/** Rebalance with a Partition object.
 * This rebalance function will use the Partition object passed
 * to perform the rebalancing.
 */
bool full_rebalance(mesh::BulkData & rebal_bulk_data ,
                    const mesh::Selector & rebal_selector ,
                    const stk::mesh::Field<double> * rebal_coord_ref ,
                    const stk::mesh::Field<double> * rebal_elem_weight_ref,
                    rebalance::Partition & partition);

}
} // namespace stk

#endif

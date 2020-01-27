// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <percept/PerceptMesh.hpp>


//----------------------------------------------------------------------

namespace stk {
  namespace rebalance {

/** \addtogroup stk_rebalance_module
 *  \{
 */

/** \brief Determine if rebalancing is needed.
 *
 * \param bulk_data      BulkData must be in a parallel consistent state.
 *
 * \param load_measure   Field defined on mesh entities of rank \a rank.
 *                       Can be a NULL pointer.
 *
 * \param imbalance_threshold  Rebalance needed if MAX divided by average load
 *                             measure exceeds this value.
 *
 * \param rank                 Rank of mesh entities to define load measure.
 *
 * \param selector             Used to select a subset of mesh entities to compute measure.
 *
 * This function calculates the total weight of the load on each processor by summing
 * the \a load_measure field over the \a selector entities of rank \a rank.  If \a selector is not
 * specified, all ejects of rank are summed.  If \a load_balance is not specified, it is
 * assumed to be 1 for each entity and the weight per processor is just the number
 * of entities on each processor.  After a processor weight is defined the MAX over the
 * processing grid is divided by the average to get a global imbalance which is
 * compared to \a imbalance_threshold.  True is returned if the global imbalance is
 * greater than \a imbalance_threshold.
 */

    void check_ownership(mesh::BulkData& bulk_data,
                         stk::mesh::EntityVector & entities,
                         const std::string& msg
                         );

    double check_balance(mesh::BulkData      & bulk_data,
                      const mesh::Field<double> * load_measure,
                      const stk::mesh::EntityRank rank,
                         const mesh::Selector      * selector=NULL,
                         stk::mesh::EntityVector *provided_elements=0, bool debug_print=true);

    bool verify_dependent_ownership(mesh::BulkData& bulk_data, const stk::mesh::EntityRank & parent_rank,
                                    stk::mesh::EntityVector & entities );

  } // namepsace rebalance
} // namepsace stk

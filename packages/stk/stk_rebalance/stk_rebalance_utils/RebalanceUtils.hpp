/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/base/Selector.hpp>


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

    double check_balance(mesh::BulkData      & bulk_data,
                      const mesh::Field<double> * load_measure,
                      const stk::mesh::EntityRank rank,
                      const mesh::Selector      * selector=NULL);

    bool verify_dependent_ownership(const stk::mesh::EntityRank & parent_rank,
                                    stk::mesh::EntityVector & entities );

  } // namepsace rebalance
} // namepsace stk

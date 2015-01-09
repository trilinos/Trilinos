/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Rebalance_Use_Cases_UseCase_9_hpp
#define Stk_Rebalance_Use_Cases_UseCase_9_hpp

#include <stk_mesh/fixtures/HexFixture.hpp>

namespace stk_classic {
namespace rebalance {
namespace use_cases {

  bool test_greedy_sideset ( stk_classic::ParallelMachine comm );

} //namespace use_cases
} //namespace rebalance
} //namespace stk_classic

#endif // Stk_Rebalance_Use_Cases_UseCase_1_hpp

/// \page stk_rebalance_use_case_4
///  \ingroup stk_rebalance_use_case_module
///
/// \section stk_rebalance_use_case_4_description Use Case 4: User-customization following default rebalance
///
/// This use case demonstrates additional user customization following a
/// default rebalance in order to enforce constraints for new partitions.  In
/// this case, the constraint is that two quad4 elements sharing edge #7
/// be collocated on the same proc following rebalance.  This is enforced using
/// a greedy sideset class which inherits the determine_new_partition method.
///
/// The following quad4 mesh is used in this use case:
///
///  Global node and element numbering
///  <pre>
///
///     13      14      15      16
///     +-------+-------+-------+
///     |       |       |       |
///     |  e7   |  e8   |  e9   |
///     |       |       |       |
///   9 +-------+-------+-------+ 12    Y
///     |       |       |       |       |
///     |  e4   |  e5   |  e6   |       |
///     |       |       |       |       *--> X
///   5 +-------+-------+-------+ 8
///     |       |       |       | 
///     |  e1   |  e2   |  e3   |
///     |       |       |       | 
///     +-------+-------+-------+ 
///     1       2       3       4
/// </pre>
///
///  Local node numbering
///  <pre>
///           3       4
///           +-------+
///           |       |
///           |  e1   |
///           |       |
///           +-------+
///           1       2
/// </pre>
/// 
/// See \ref UseCase_Rebal_4.cpp for the complete source listing.
/// 
/// 
/// 
/// 

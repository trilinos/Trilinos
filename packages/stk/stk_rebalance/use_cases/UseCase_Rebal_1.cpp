/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Rebal_1.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace rebalance {
namespace use_cases {

UseCase_1_Rebalance::UseCase_1_Rebalance( stk::ParallelMachine comm )
{
}

UseCase_1_Rebalance::~UseCase_1_Rebalance()
{ }



} //namespace use_cases
} //namespace rebalance
} //namespace stk



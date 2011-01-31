/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>

#include <use_cases/UseCase_Rebal_1.hpp>
#include <use_cases/UseCase_Rebal_2.hpp>
#include <use_cases/UseCase_Rebal_3.hpp>
#include <use_cases/UseCase_Rebal_4.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_rebalance/Rebalance.hpp>

void printStatus(bool status)
{
  if (status) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }
}

int main ( int argc, char * argv[] )
{
  stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);
  //const bool single_process =
  //  stk::parallel_machine_size( parallel_machine ) <= 1 ;

  bool status = true;
  {
    std::cout << "Use Case 1, unequal element weights ... ";
    bool local_status = stk::rebalance::use_cases::test_unequal_weights(parallel_machine);
    stk::all_reduce(parallel_machine, stk::ReduceMin<1>(&local_status));
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case 2, heavy entities ... ";
    bool local_status = stk::rebalance::use_cases::test_heavy_nodes(parallel_machine);
    stk::all_reduce(parallel_machine, stk::ReduceMin<1>(&local_status));
    printStatus(local_status);
    status = status && local_status;
  }

  {
    std::cout << "Use Case 3, contact surfaces ... ";
    bool local_status = stk::rebalance::use_cases::test_contact_surfaces(parallel_machine);
    stk::all_reduce(parallel_machine, stk::ReduceMin<1>(&local_status));
    printStatus(local_status);
    status = status && local_status;
  }
  {
    std::cout << "Use Case 4, greedy sideset ... ";
    bool local_status = stk::rebalance::use_cases::test_greedy_sideset(parallel_machine);
    stk::all_reduce(parallel_machine, stk::ReduceMin<1>(&local_status));
    printStatus(local_status);
    status = status && local_status;
  }

  int return_code = -1;
  if (status) {
    return_code = 0;
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  else {
    std::cout << "End Result: TEST FAILED" << std::endl;
  }

  stk::parallel_machine_finalize();

  return return_code;
}

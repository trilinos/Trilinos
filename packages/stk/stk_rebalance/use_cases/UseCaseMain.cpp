/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>

#include <use_cases/UseCase_Rebal_1.hpp>

#include <stk_util/parallel/Parallel.hpp>
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

   //if ( single_process ) {
    std::cout << "Use Case 1 ... ";
    bool local_status = true ;
    try {
      local_status = stk::rebalance::use_cases::test_unequal_weights(parallel_machine);
      printStatus(local_status);
    }
    catch ( const std::exception & x ) {
      local_status = false ;
      printStatus(local_status);
      std::cout << x.what();
    }
    status = status && local_status;
   //}

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

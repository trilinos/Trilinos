/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/parallel/Parallel.hpp>

/*--------------------------------------------------------------------*/
/* Parallel operations */

#if defined( STK_HAS_MPI )

namespace stk {

unsigned parallel_machine_size( ParallelMachine parallel_machine )
{
  int value = 0 ;
  if (parallel_machine != MPI_COMM_NULL) {
    if ( MPI_SUCCESS != MPI_Comm_size( parallel_machine , &value ) ) {
      value = 0 ;
    }
  }
  return value ;
}

unsigned parallel_machine_rank( ParallelMachine parallel_machine )
{
  int value = 0 ;
  if (parallel_machine != MPI_COMM_NULL) {
    if ( MPI_SUCCESS != MPI_Comm_rank( parallel_machine , &value ) ) {
      value = 0 ;
    }
  }
  return value ;
}

void parallel_machine_barrier( ParallelMachine parallel_machine )
{
  if (parallel_machine != MPI_COMM_NULL) {
    MPI_Barrier( parallel_machine );
  }
}

}

#else

namespace stk {

unsigned parallel_machine_size( ParallelMachine parallel_machine) { return 1 ; }

unsigned parallel_machine_rank( ParallelMachine parallel_machine) { return 0 ; }

void parallel_machine_barrier( ParallelMachine parallel_machine) {}

}

#endif

/*--------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/RunEnvironment.hpp>

#include <stk_percept/pyencore.h>

#if !PY_PERCEPT 
STKUNIT_MAIN(argc, argv)
#else
  int main() {return 0;}
#endif

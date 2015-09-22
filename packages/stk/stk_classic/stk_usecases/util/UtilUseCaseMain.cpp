/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

int use_case_message(use_case::UseCaseEnvironment &use_case_environment);
int use_case_timer(use_case::UseCaseEnvironment &use_case_environment);

int
main(
  int           argc,
  char **        argv)
{
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  use_case_message(use_case_environment);
  use_case_timer(use_case_environment);
}

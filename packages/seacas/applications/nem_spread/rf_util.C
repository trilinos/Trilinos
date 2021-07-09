/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <cstdio>  // for stderr
#include <cstdlib> // for exit
#include <fmt/ostream.h>

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
void check_exodus_error(int error, const char *function_name)
{

  if (error == -1) {
    fmt::print(stderr, "ERROR returned from {}!\n", function_name);
    exit(1);
  }

} /* check_exodus_error */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_line(const char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes; i++) {
    fmt::print("{}", *charstr);
  }
  fmt::print("\n");
}

/*****************************************************************************/
/*                END OF FILE rf_util.c                                      */
/*****************************************************************************/

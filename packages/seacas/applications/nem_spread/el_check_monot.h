/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <cstdlib>

template <typename INT> int check_monot(INT *vector, size_t length)
{
  /*
   * This function checks to see if an integer vector is in monotonically
   * increasing order. It returns TRUE (i.e., 1), if this is so. It returns
   * false (i.e., 0), if this is not so.
   */
  size_t i;
  for (i = 1; i < length; i++) {
    if (vector[i] < vector[i - 1]) {
      return (0);
    }
  }
  return (1);
}

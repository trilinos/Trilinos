// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

int main(int narg, char **arg)
{
  if (narg > 1)
    std::cout << arg[1] << std::endl;
  else
    std::cout << "BUMMER" << std::endl;
  return 0;
}

// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

int main(int argc, char **argv)
{
  Kokkos::initialize(argc,argv);

  testing::InitGoogleTest(&argc, argv);

  int returnVal = RUN_ALL_TESTS();

  Kokkos::finalize();

  return returnVal;
}

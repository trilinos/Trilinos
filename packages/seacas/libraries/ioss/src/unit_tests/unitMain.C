// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <gtest/gtest.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  testing::InitGoogleTest(&argc, argv);
  int errorCode = RUN_ALL_TESTS();

  MPI_Finalize();

  return errorCode;
}

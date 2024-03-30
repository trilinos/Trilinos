// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "iocatalyst_export.h"
#include <gtest/gtest.h>
#include <SEACASIoss_config.h>
#include <Ionit_Initializer.h>

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

int IOCATALYST_EXPORT main(int argc, char **argv)
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  Ioss::Init::Initializer io;
  testing::InitGoogleTest(&argc, argv);
  int errorCode = RUN_ALL_TESTS();

#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif

  return errorCode;
}

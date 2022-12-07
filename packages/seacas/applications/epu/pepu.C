/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

// concatenates EXODUS/GENESIS output from parallel processors to a single file

#include <exodusII.h>
#ifdef SEACAS_HAVE_MPI
#ifndef DISABLE_PARALLEL_EPU
#define ENABLE_PARALLEL_EPU 1
#endif
#endif

#if ENABLE_PARALLEL_EPU
#include <mpi.h>

#include "epu.C"
#else
#error "PEPU: Parallel capability not enabled. This should not happen."
#endif

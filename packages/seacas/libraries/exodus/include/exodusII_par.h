/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 * NOTE: Contents of this include file have been moved to exodusII.h
 * Retained here just for backward compatibility
 */
#pragma once

#include "exodusII.h"

#if !defined(PARALLEL_AWARE_EXODUS)
#error "Parallel-aware exodusII_par.h included in non-parallel context"
#endif

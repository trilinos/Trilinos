/*
 * NOTE: Contents of this include file have been moved to exodusII.h
 * Retained here just for backward compatibility
 */

#ifndef EXODUS_II_PAR_HDR
#define EXODUS_II_PAR_HDR

#include "exodusII.h"

#if !defined(PARALLEL_AWARE_EXODUS)
#error "Parallel-aware exodusII_par.h included in non-parallel context"
#endif

#endif

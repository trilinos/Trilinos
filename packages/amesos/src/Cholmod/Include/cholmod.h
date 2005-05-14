/* ========================================================================== */
/* === Include/cholmod.h ==================================================== */
/* ========================================================================== */

/*
 * CHOLMOD version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis,
 * William W. Hager, and Sivasankaran Rajamanickam.  Note that each of
 * CHOLMOD's modules are licensed separately.  The GNU LGPL license applies to
 * the Core, Check, and Partition modules.  The Cholesky, MatrixOps, Modify,
 * and Supernodal modules are not yet released; the GNU LGPL license will NOT
 * directly apply to them.
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD include file, for inclusion user programs. */

#ifndef CHOLMOD_H
#define CHOLMOD_H

#include "cholmod_config.h"

/* CHOLMOD always includes the Core module. */
#include "cholmod_core.h"

#ifndef NCHECK
#include "cholmod_check.h"
#endif

#ifndef NCHOLESKY
#include "cholmod_cholesky.h"
#endif

#ifndef NMATRIXOPS
#include "cholmod_matrixops.h"
#endif

#ifndef NMODIFY
#include "cholmod_modify.h"
#endif

#ifndef NPARTITION
#include "cholmod_partition.h"
#endif

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

#endif

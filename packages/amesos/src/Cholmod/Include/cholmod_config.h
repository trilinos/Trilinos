/* ========================================================================== */
/* === Include/cholmod_config.h ============================================= */
/* ========================================================================== */

/*
 * CHOLMOD version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis,
 * William W. Hager, and Sivasankaran Rajamanickam.  Note that each of
 * CHOLMOD's modules are licensed separately.  The GNU LGPL license applies to
 * the Core, Check, Cholesky, and Partition modules.  The MatrixOps, Modify,
 * and Supernodal modules are not yet released; the GNU LGPL license will NOT
 * directly apply to them.  Since this file is required by all modules,
 * the GNU LGPL does apply to this file.
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD configuration file, for inclusion in user programs.
 *
 * You do not have to edit any CHOLMOD files to compile and install CHOLMOD.
 * However, if you do not use all of CHOLMOD's modules, you need to compile
 * with the appropriate flag, or edit this file to add the appropriate #define.
 *
 * Only two are required.  If you do not have the Partition module, then to use
 * the Cholesky module you must compile with -DNPARTITIION.  Likewise, to use
 * the Cholesky module without the Supernodal modules, you must compile with
 * -DNSUPERNODAL.
 *
 * Compiler flags for CHOLMOD:
 *
 * -DNCHECK	    do not include the Check module.
 * -DNCHOLESKY	    do not include the Cholesky module.
 * -DNMATRIXOPS	    do not include the MatrixOps module.
 * -DNMODIFY	    do not include the Modify module.
 * -DNPARTITION	    do not include the Partition module.
 * -DNSUPERNODAL    do not include the Supernodal module.
 *
 * -DNPRINT	    do not print anything; do not include <stdio.h>
 *		    (Note that METIS includes <stdio.h>, which is used by the
 *		    Partition module).
 */

#ifndef CHOLMOD_CONFIG_H
#define CHOLMOD_CONFIG_H

/* Use the compiler flag, or uncomment the definition(s), if you want to use
 * one or more non-default installation options:

#define NCHECK
#define NCHOLESKY
#define NPARTITION
#define NPRINT
 */

/* This version is a subset of CHOLMOD, used for PARAKLETE version 0.1 */

#define NMATRIXOPS
#define NMODIFY
#define NSUPERNODAL

#endif

/* ========================================================================== */
/* === Include/cholmod_config.h ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_config.h.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod_config.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* CHOLMOD configuration file, for inclusion in user programs.
 *
 * You do not have to edit any CHOLMOD files to compile and install CHOLMOD.
 * However, if you do not use all of CHOLMOD's modules, you need to compile
 * with the appropriate flag, or edit this file to add the appropriate #define.
 *
 * Compiler flags for CHOLMOD:
 *
 * -DNCHECK	    do not include the Check module.        License: GNU LGPL
 * -DNCHOLESKY	    do not include the Cholesky module.     License: GNU LGPL
 * -DNPARTITION	    do not include the Partition module.    License: GNU LGPL
 *
 * -DNPRINT	    do not print anything
 *
 * -D'LONGBLAS=long' or -DLONGBLAS='long long' defines the integers used by
 *		    LAPACK and the BLAS.  Use LONGBLAS=long on Solaris to use
 *		    the 64-bit Sun Performance BLAS in cholmod_l_* routines.
 *		    You may need to use -D'LONGBLAS=long long' on the SGI
 *		    (this is not tested).
 *
 * -DNSUNPERF	    for Solaris only.  If defined, do not use the Sun
 *		    Performance Library.  The default is to use SunPerf.
 *		    You must compile CHOLMOD with -xlic_lib=sunperf.
 *
 * The Core Module (License GNU LGPL) is always included in the CHOLMOD library.
 */

#ifndef AMESOS_CHOLMOD_CONFIG_H
#define AMESOS_CHOLMOD_CONFIG_H

/* Use the compiler flag, or uncomment the definition(s), if you want to use
 * one or more non-default installation options: */

/*
#define NCHECK
#define NCHOLESKY
#define NPARTITION

#define NPRINT

#define LONGBLAS long
#define LONGBLAS long long
#define NSUNPERF
*/

/* Turning off all code that uses the GPL'ed modules */
#define NMATRIXOPS
#define NMODIFY
#define NSUPERNODAL

#endif

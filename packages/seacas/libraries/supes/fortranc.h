/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 * Define the Fortran/C interface for the system in question
 */
#ifndef _FORTRANC_
#define _FORTRANC_

/* Default system */
#if Build64
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#define FTNREAL double
#define FTNINT  long long int
#else
#define FTNREAL double
#define FTNINT  long int
#endif
#else
#define FTNREAL float
#define FTNINT  int
#endif

#endif /* _FORTRANC_ */

/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
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
#define FTNREAL double
#define FTNINT long int
#else
#define FTNREAL float
#define FTNINT int
#endif

#endif /* _FORTRANC_ */

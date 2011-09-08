/*
 * Define the Fortran/C interface for the system in question
 */
#ifndef _FORTRANC_
#define _FORTRANC_

/* Default system */
#if Build64
#    define FTNREAL double
#    define FTNINT  long int
#else
#    define FTNREAL float
#    define FTNINT  int
#endif

#endif /* _FORTRANC_ */


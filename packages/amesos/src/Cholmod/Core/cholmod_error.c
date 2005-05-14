/* ========================================================================== */
/* === Core/cholmod_error =================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD error-handling routine.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

#ifndef NPRINT
#ifdef MATLAB_MEX_FILE
/* MATLAB mexFunction: print error message on Command window */
#include "mex.h"
#define PRINT_MESSAGE(s,msg) mexPrintf ("CHOLMOD %s: %s\n", s, msg) ;
#else
/* ANSI C: print error message on Common->file */
#include <stdio.h>
#define FP ((Common->file == NULL) ? stdout : (Common->file))
#define PRINT_MESSAGE(s,msg) fprintf (FP, "CHOLMOD %s: %s\n", s, msg) ;
#endif
#else
/* no printing */
#define PRINT_MESSAGE(s,msg)
#endif


/* ========================================================================== */
/* ==== cholmod_error ======================================================= */
/* ========================================================================== */

/* An error has occurred.  Set the status, optionally print an error message,
 * and call the user error-handling routine (if it exists).  If
 * Common->try_catch is TRUE, then CHOLMOD is inside a try/catch block.
 * The status is set, but no message is printed and the user error handler
 * is not called.  This is not (yet) an error, since CHOLMOD may recover.
 *
 * In CHOLMOD v1.0, this try/catch mechanism is used only in cholmod_analyze,
 * which tries multiple ordering methods and picks the best one.  If one or
 * more ordering method fails, it keeps going.  Only one ordering needs to
 * succeed for cholmod_analyze to succed.
 */

int cholmod_error
(
    int status,
    char *msg,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;

    Common->status = status ;

    if (!(Common->try_catch))
    {

	/* print a warning or error message */
	if (status > 0 && Common->print > 1)
	{
	    PRINT_MESSAGE ("warning", msg) ;
	}
	else if (Common->print > 0)
	{
	    PRINT_MESSAGE ("ERROR", msg) ;
	}

	/* call the user error handler, if it exists */
	if (Common->error_handler != NULL)
	{
	    Common->error_handler (status, msg) ;
	}
    }

    return (TRUE) ;
}

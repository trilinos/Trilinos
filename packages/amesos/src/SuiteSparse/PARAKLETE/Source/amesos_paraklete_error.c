/* ========================================================================== */
/* === paraklete_error ====================================================== */
/* ========================================================================== */

#include "amesos_paraklete_decl.h"

/* Error handling for PARAKLETE.  Currently, this function aborts the entire
 * set of MPI processes.  It does this because coordinating a clean exit among
 * the processes for all possible failure modes is very difficult.  Future work:
 * remove the call to MPI_Abort from this function, and handle the errors more
 * gracefully without terminating the user's parallel application.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

void amesos_paraklete_error (Int status, char *filename, Int line, char *message,
    paraklete_common *Common)
{
    Common->status = status ;
    printf ("Error in proc "ID", status "ID" in file %s, line "ID" : %s\n",
        Common->myid, status, filename, line, message) ;
    CHOLMOD (error) ((int) status, filename, (int) line, message,
        &(Common->cm)) ;
    MPI (MPI_Abort (MPI_COMM_WORLD, status)) ;
    /* in case this version is not compiled with MPI */
    abort ( ) ;
}

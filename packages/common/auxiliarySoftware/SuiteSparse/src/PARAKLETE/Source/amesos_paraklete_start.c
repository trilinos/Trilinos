/* ========================================================================== */
/* === paraklete_start ====================================================== */
/* ========================================================================== */

#include "amesos_paraklete_decl.h"

/* paraklete_start (nproc, myid, Common) must be called by all processes prior
 * to using any PARAKLETE or CHOLMOD functions.
 *
 * Normally, the nproc and myid inputs to this function are from here:
 *
 *   MPI_Comm_size (MPI_COMM_WORLD, &nproc) ;
 *   MPI_Comm_rank (MPI_COMM_WORLD, &myid) ;
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

void amesos_paraklete_start (Int nproc, Int myid, paraklete_common *Common)
{
    cholmod_common *cm ;
    cm = &(Common->cm) ;
    CHOLMOD (start) (cm) ;
    PK_DEBUG_INIT ("pk", cm) ;
    Common->nproc = nproc ;
    Common->nleaves = 32 ;
    Common->myid = myid ;
    Common->file = NULL ;
    cm->print = 1 ;
    cm->precise = TRUE ;
    Common->tol_diag = 0.01 ;
    Common->tol_offdiag = 1.0 ;
    Common->growth = 2. ;
    Common->dump = 0 ;
    TRILINOS_KLU_defaults (&(Common->km)) ;
    Common->km.btf = 0 ;
}


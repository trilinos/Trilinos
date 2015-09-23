/* ========================================================================== */
/* === paraklete_finish ===================================================== */
/* ========================================================================== */

#include "amesos_paraklete_decl.h"

/* paraklete_finish (Common) must be called by all processes after
 * using any PARAKLETE and/or CHOLMOD functions.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

void amesos_paraklete_finish (paraklete_common *Common)
{
    cholmod_common *cm ;
    cm = &(Common->cm) ;
    CHOLMOD (finish) (cm) ;
    if (cm->malloc_count != 0 || cm->memory_inuse != 0)
    {
	printf ("%% proc "ID" finish: "ID" "ID"\n", Common->myid,
            (Int) (cm->malloc_count), (Int) (cm->memory_inuse)) ;
    }
    ASSERT (cm->malloc_count == 0 && cm->memory_inuse == 0) ;
    if (Common->file != NULL)
    {
        fclose (Common->file) ;
    }
    Common->file = NULL ;
}


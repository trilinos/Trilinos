/* ========================================================================== */
/* === klu_btf_defaults ===================================================== */
/* ========================================================================== */

#include "klu_btf_internal.h"

void klu_btf_defaults
(
    klu_control *control
)
{
    if (control != (klu_control *) NULL)
    {
	control->tol = 0.001 ;		/* pivot tolerance for diagonal */
	control->growth = 1.5 ;		/* realloc growth size */
	control->initmem_amd = 1.2 ;	/* init. mem with AMD:  c*nnz(L) + n */
	control->initmem = 10 ;		/* init. mem otherwise: c*nnz(A) + n */
	control->btf = TRUE ;		/* use BTF pre-ordering, or not */
	control->ordering = 0 ;		/* 0: AMD, 1: COLAMD, 2: user P and Q */
	control->scale = 1 ;		/* scale: 0: none, 1: sum, 2: max */
    }
}

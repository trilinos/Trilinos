/* ========================================================================== */
/* === klu_btf_free_numeric ================================================= */
/* ========================================================================== */

/* Free the KLU Numeric object. */

#include "klu_btf_internal.h"

void klu_btf_free_numeric (klu_numeric **NumericHandle)
{
    klu_numeric *Numeric ;
    int block, nblocks, **Lbp, **Lbi, **Ubp, **Ubi ;
    double **Lbx, **Ubx ;

    if (NumericHandle == (klu_numeric **) NULL)
    {
	/* nothing to do */
	return ;
    }

    Numeric = *NumericHandle ;
    nblocks = Numeric->nblocks ;

    FREE (Numeric->Pnum, int) ;
    FREE (Numeric->Offp, int) ;
    FREE (Numeric->Offi, int) ;
    FREE (Numeric->Offx, double) ;
    FREE (Numeric->Singleton, double) ;
    FREE (Numeric->Rs, double) ;
    FREE (Numeric->Pinv, int) ;
    FREE (Numeric->X, double) ;

    Lbp = Numeric->Lbp ;
    Lbi = Numeric->Lbi ;
    Lbx = Numeric->Lbx ;

    Ubp = Numeric->Ubp ;
    Ubi = Numeric->Ubi ;
    Ubx = Numeric->Ubx ;

    for (block = 0 ; block < nblocks ; block++)
    {
	if (Lbp != (int    **) NULL) FREE (Lbp [block], int) ;
	if (Lbi != (int    **) NULL) FREE (Lbi [block], int) ;
	if (Lbx != (double **) NULL) FREE (Lbx [block], double) ;
	if (Ubp != (int    **) NULL) FREE (Ubp [block], int) ;
	if (Ubi != (int    **) NULL) FREE (Ubi [block], int) ;
	if (Ubx != (double **) NULL) FREE (Ubx [block], double) ;
    }

    FREE (Numeric->Lbp, int *) ;
    FREE (Numeric->Lbi, int *) ;
    FREE (Numeric->Lbx, double *) ;

    FREE (Numeric->Ubp, int *) ;
    FREE (Numeric->Ubi, int *) ;
    FREE (Numeric->Ubx, double *) ;

    FREE (Numeric, klu_numeric) ;
}

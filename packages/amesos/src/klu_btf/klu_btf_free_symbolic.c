#include "klu_btf_internal.h"

void klu_btf_free_symbolic (klu_symbolic **SymbolicHandle)
{
    klu_symbolic *Symbolic ;

    if (SymbolicHandle == (klu_symbolic **) NULL)
    {
	/* nothing to do */
	return ;
    }

    Symbolic = *SymbolicHandle ;

    FREE (Symbolic->P, int) ;
    FREE (Symbolic->Q, int) ;
    FREE (Symbolic->R, int) ;
    FREE (Symbolic->Lnz, double) ;
    FREE (Symbolic, klu_symbolic) ;
}

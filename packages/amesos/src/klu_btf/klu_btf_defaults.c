/* ========================================================================== */
/* === klu_btf_defaults ===================================================== */
/* ========================================================================== */

#include "klu_btf_internal.h"

void klu_btf_defaults
(
    double Control [KLU_BTF_CONTROL]
)
{
    int i ;

    if (!Control)
    {
	/* silently return if no Control array */
	return ;
    }

    for (i = 0 ; i < KLU_BTF_CONTROL ; i++)
    {
	Control [i] = 0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* default control settings: can be modified at run-time */
    /* ---------------------------------------------------------------------- */

    /* TODO used by printing routines */
    Control [KLU_BTF_CONTROL_PRL] = KLU_BTF_CONTROL_PRL_DEFAULT ;

    /* used in klu_btf_analyze */
    Control [KLU_BTF_CONTROL_BTF] = KLU_BTF_CONTROL_BTF_DEFAULT ;
    Control [KLU_BTF_CONTROL_AMD_DENSE] = KLU_BTF_CONTROL_AMD_DENSE_DEFAULT ;
    Control [KLU_BTF_CONTROL_ORDERING] = KLU_BTF_CONTROL_ORDERING_DEFAULT ;

    /* used in klu_btf_factor */
    Control [KLU_BTF_CONTROL_TOL] = KLU_BTF_CONTROL_TOL_DEFAULT ;
    Control [KLU_BTF_CONTROL_GROWTH] = KLU_BTF_CONTROL_GROWTH_DEFAULT ;
    Control [KLU_BTF_CONTROL_INITMEM_AMD] = KLU_BTF_CONTROL_INITMEM_AMD_DEFAULT;
    Control [KLU_BTF_CONTROL_INITMEM_COLAMD] =
	KLU_BTF_CONTROL_INITMEM_COLAMD_DEFAULT ;
}

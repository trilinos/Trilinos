/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include "superlu_ddefs.h"

void pxerbla(char *srname, gridinfo_t *grid, int_t info)
{
    printf("{%4d,%4d}: On entry to %6s, parameter number %2d had an illegal value\n",
	   MYROW(grid->iam, grid), MYCOL(grid->iam, grid), srname, info);

}

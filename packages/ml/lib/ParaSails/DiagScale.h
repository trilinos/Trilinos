/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * DiagScale.h header file.
 *
 *****************************************************************************/

#include <stdio.h>
#include "Hash.h"
#include "Matrix.h"
#include "Numbering.h"

#ifndef _DIAGSCALE_H
#define _DIAGSCALE_H

typedef struct
{
    int     offset;      /* number of on-processor entries */
    double *local_diags; /* on-processor entries */
    double *ext_diags;   /* off-processor entries */
}
DiagScale;

DiagScale *DiagScaleCreate(Matrix *A, Numbering *numb);
void DiagScaleDestroy(DiagScale *p);
double DiagScaleGet(DiagScale *p, int index);

#endif /* _DIAGSCALE_H */

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
 * RowPatt.h header file.
 *
 *****************************************************************************/

#include <stdio.h>

#ifndef _ROWPATT_H
#define _ROWPATT_H

typedef struct
{
    int  maxlen;
    int  len;
    int  prev_len;
    int *ind;
    int *mark;
    int *buffer; /* buffer used for outputting indices */
    int  buflen; /* length of this buffer */
}
RowPatt;

RowPatt *RowPattCreate(int maxlen);
void RowPattDestroy(RowPatt *p);
void RowPattReset(RowPatt *p);
void RowPattMerge(RowPatt *p, int len, int *ind);
void RowPattMergeExt(RowPatt *p, int len, int *ind, int num_loc);
void RowPattGet(RowPatt *p, int *lenp, int **indp);
void RowPattPrevLevel(RowPatt *p, int *lenp, int **indp);

#endif /* _ROWPATT_H */

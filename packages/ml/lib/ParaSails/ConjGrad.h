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
 * ConjGrad.h header file.
 *
 *****************************************************************************/

#ifndef _CONJGRAD_H
#define _CONJGRAD_H

void PCG_ParaSails(Matrix *mat, ParaSails *ps, double *b, double *x,
   double tol, int max_iter);
void FGMRES_ParaSails(Matrix *mat, ParaSails *ps, double *b, double *x,
   int dim, double tol, int max_iter);

#endif /* _CONJGRAD_H */

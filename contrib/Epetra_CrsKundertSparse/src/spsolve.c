/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_file_id = "$Id$";
#endif

/*
 *  MATRIX SOLVE MODULE
 *
 *  Author:                     Advising professor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains the forward and backward substitution routines for
 *  the sparse matrix routines.
 *
 *  >>> User accessible functions contained in this file:
 *  spSolve
 *  spSolveTransposed
 *
 *  >>> Other functions contained in this file:
 *  SolveComplexMatrix
 *  SolveComplexTransposedMatrix
 */


/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88,89,90
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and
 *  its documentation for any purpose and without fee is hereby granted,
 *  provided that the copyright notices appear in all copies and
 *  supporting documentation and that the authors and the University of
 *  California are properly credited.  The authors and the University of
 *  California make no representations as to the suitability of this
 *  software for any purpose.  It is provided `as is', without express
 *  or implied warranty.
 */

#ifdef notdef
static char copyright[] =
    "Sparse1.3: Copyright (c) 1985,86,87,88,89,90 by Kenneth S. Kundert";
static char RCSid[] =
    "@(#)$Header$";
#endif

#include "spice.h"

#ifdef SHARED_MEM
#include "shared_mem.h"
#include <stdlib.h>
#endif

/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spConfig.h
 *     Macros that customize the sparse matrix routines.
 *  spMatrix.h
 *     Macros and declarations to be imported by the user.
 *  spDefs.h
 *     Matrix type and macro definitions for the sparse matrix routines.
 */

#define spINSIDE_SPARSE
#include "spconfig.h"
#include "spmatrix.h"
#include "spdefs.h"

double total_solve_time;


/*
 * Function declarations
 */

#ifdef __STDC__
#if spSEPARATED_COMPLEX_VECTORS
static void SolveComplexMatrix( MatrixPtr,
                        RealVector, RealVector, RealVector, RealVector );
static void SolveComplexTransposedMatrix( MatrixPtr,
                        RealVector, RealVector, RealVector, RealVector );
#else
static void SolveComplexMatrix( MatrixPtr, RealVector, RealVector );
static void SolveComplexTransposedMatrix( MatrixPtr, RealVector, RealVector );
#endif
#else /* __STDC__ */
static void SolveComplexMatrix();
static void SolveComplexTransposedMatrix();
#endif /* __STDC__ */
/*
#define WRITE_MAT
*/


/*
 *  SOLVE MATRIX EQUATION
 *
 *  Performs forward elimination and back substitution to find the
 *  unknown vector from the RHS vector and factored matrix.  This
 *  routine assumes that the pivots are associated with the lower
 *  triangular (L) matrix and that the diagonal of the upper triangular
 *  (U) matrix consists of ones.  This routine arranges the computation
 *  in different way than is traditionally used in order to exploit the
 *  sparsity of the right-hand side.  See the reference in spRevision.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  RHS  <input>  (RealVector)
 *      RHS is the input data array, the right hand side. This data is
 *      undisturbed and may be reused for other solves.
 *  Solution  <output>  (RealVector)
 *      Solution is the output data array. This routine is constructed such that
 *      RHS and Solution can be the same array.
 *  iRHS  <input>  (RealVector)
 *      iRHS is the imaginary portion of the input data array, the right
 *      hand side. This data is undisturbed and may be reused for other solves.
 *      This argument is only necessary if matrix is complex and if
 *      spSEPARATED_COMPLEX_VECTOR is set true.
 *  iSolution  <output>  (RealVector)
 *      iSolution is the imaginary portion of the output data array. This
 *      routine is constructed such that iRHS and iSolution can be
 *      the same array.  This argument is only necessary if matrix is complex
 *      and if spSEPARATED_COMPLEX_VECTOR is set true.
 *
 *  >>> Local variables:
 *  Intermediate  (RealVector)
 *      Temporary storage for use in forward elimination and backward
 *      substitution.  Commonly referred to as c, when the LU factorization
 *      equations are given as  Ax = b, Lc = b, Ux = c Local version of
 *      Matrix->Intermediate, which was created during the initial
 *      factorization in function spcCreateInternalVectors() in the matrix
 *      factorization module.
 *  pElement  (ElementPtr)
 *      Pointer used to address elements in both the lower and upper triangle
 *      matrices.
 *  pExtOrder  (int *)
 *      Pointer used to sequentially access each entry in IntToExtRowMap
 *      and IntToExtColMap arrays.  Used to quickly scramble and unscramble
 *      RHS and Solution to account for row and column interchanges.
 *  pPivot  (ElementPtr)
 *      Pointer that points to current pivot or diagonal element.
 *  Size  (int)
 *      Size of matrix. Made local to reduce indirection.
 *  Temp  (RealNumber)
 *      Temporary storage for entries in arrays.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

/*VARARGS3*/

/* #define FINGERPRINT */

#ifdef FINGERPRINT
void FingerPrint (double *RHS, int *RowMap, int *ColMap, int Size)
{
    int i;
    double Temp, finger;
    unsigned long Row, Col, RCmask;

    RCmask = 1;
    for (i=0 ; i<8*(sizeof(long)-2)-2 ; i++) {
      RCmask <<= 1;
      RCmask += 1;
    }
    finger = 0;
    Row = 0;
    Col = 0;
    for (i=1 ; i<=Size ; i++) {
      Row = Row*99991 + RowMap[i];
      Row &= RCmask;
      Col = Col*99991 + ColMap[i];
      Col &= RCmask;
      if (RHS) {
        Temp = RHS[i];
        finger += Temp;
/*      printf ("RHS[%d] = %g\n",i,RHS[i]); */
      }
    }
    printf ("RHS Fingerprint: %.12g,   Row: %ld,   Col: %ld\n",finger, Row, Col);
    return;
}
#endif

extern int *col_start;

void spExpandFormat(MatrixPtr);

int
spSolve( eMatrix, RHS, Solution IMAG_VECTORS )

char *eMatrix;
RealVector  RHS, Solution IMAG_VECTORS;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement;
register  RealVector  Intermediate;
register  RealNumber  Temp;
register  int  I, j, *pExtOrder, Size;
ElementPtr  pPivot;
void SolveComplexMatrix();
static int numCall, nextCheck, nextInterval;
#ifdef WRITE_MAT
FILE *fxb;
double *orig_rhs;
#endif

/* Begin `spSolve'. */
    ASSERT( IS_VALID(Matrix));
    spExpandFormat(Matrix);
    ASSERT( IS_FACTORED(Matrix) );
#ifndef SHARED_MEM
#endif

#if spCOMPLEX
    if (Matrix->Complex)
    {   SolveComplexMatrix( Matrix, RHS, Solution IMAG_VECTORS );
        return 0;
    }
#endif

#if REAL
    Intermediate = Matrix->Intermediate;
    Size = Matrix->Size;

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
    --RHS;
    --Solution;
#endif

#ifdef WRITE_MAT
    orig_rhs = (double *) tmalloc((1+Size)*sizeof(double));
    memcpy (orig_rhs, RHS, (1+Size)*sizeof(double));
#endif

/* Initialize Intermediate vector. */
    pExtOrder = &Matrix->IntToExtRowMap[Size];
    for (I = Size; I > 0; I--) {
#if REORDER_SCALING
      if (Matrix->has_scale_factors)
        Intermediate[I] = RHS[*pExtOrder]*Matrix->row_scale_factors[*pExtOrder];
      else
        Intermediate[I] = RHS[*pExtOrder];
#else
      Intermediate[I] = RHS[*pExtOrder];
#endif
      pExtOrder--;
    }

/* Forward elimination. Solves Lc = b.*/
    for (I = 1; I <= Size; I++) {   
      if ((Temp = Intermediate[I]) != 0.0) {
        pPivot = Matrix->Diag[I];
        Intermediate[I] = (Temp *= pPivot->Real);

        pElement = pPivot->NextInCol;
        while (pElement != NULL) {
            Intermediate[pElement->Row] -= Temp * pElement->Real;
            pElement = pElement->NextInCol;
        }
      }
    }
    if (Matrix->Error == spZERO_DIAG  || Matrix->Error == spSINGULAR)
      return (Matrix->Error);

/* Backward Substitution. Solves Ux = c.*/
    pExtOrder = &Matrix->IntToExtColMap[Size];
    for (I = Size; I > 0; I--)
    {   Temp = Intermediate[I];
        pElement = Matrix->Diag[I]->NextInRow;
        j = 0;
        while (pElement != NULL) {
          Matrix->Intermediate2[j++] = pElement->Real * Intermediate[pElement->Col];
          pElement = pElement->NextInRow;
        }
        while (j>0) {
          Temp -= Matrix->Intermediate2[--j];
        }
        Intermediate[I] = Temp;
        Solution[*(pExtOrder--)] = Intermediate[I];
    }
#if REORDER_SCALING
    if (Matrix->has_scale_factors) {
      for (I=1 ; I<=Size ; I++) {
/*
        printf ("Col scale factor[%d] = %g\n", I, Matrix->col_scale_factors[I]);
*/
        Solution[I] *= Matrix->col_scale_factors[I];
      }
    }
#endif
#ifdef WRITE_MAT
    fxb = fopen ("xb.vec","w");
    for (I = 1; I<=Size; I++) {
      fprintf (fxb,"%.12g %.12g\n",orig_rhs[I],Solution[I]);
    }
    fclose(fxb);
    tfree(orig_rhs);
#endif

/*
    for (I=1 ; I<=Size ; I++) {
      printf ("%d (%10s) :: %g %g\n",I,CKTnodName(Matrix->Ckt,I), RHS[I],Solution[I]);
    }
*/

#ifdef FINGERPRINT
    FingerPrint(RHS, Matrix->IntToExtRowMap, Matrix->IntToExtColMap, Size);
#endif

#ifdef SHARED_MEM
    numCall++;
#endif
    return 0;
#endif /* REAL */
}

#ifdef SELECTIVE_SOLVE
int
spSolveBack( eMatrix, RHS, Solution IMAG_VECTORS )

char *eMatrix;
RealVector  RHS, Solution IMAG_VECTORS;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement;
register  RealVector  Intermediate;
register  RealNumber  Temp;
register  int  I, *pExtOrder, Size;
ElementPtr  pPivot;
void SolveComplexMatrix();

/* Begin `spSolve'. */
    ASSERT( IS_VALID(Matrix));

#if spCOMPLEX
    if (Matrix->Complex)
    {   SolveComplexMatrix( Matrix, RHS, Solution IMAG_VECTORS );
        return 0;
    }
#endif

#if REAL
    Intermediate = Matrix->Intermediate;
    Size = Matrix->Size;

/*
    for (I=0 ; I<Size ; I++) {
      printf ("%.10g  ",RHS[I]);
    }
    printf ("\n");
*/
/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
    --RHS;
    --Solution;
#endif

/* Initialize Intermediate vector. */
    pExtOrder = &Matrix->IntToExtRowMap[Size];
    for (I = Size; I > 0; I--)
        Intermediate[I] = RHS[*(pExtOrder--)];

/* Forward elimination. Solves Lc = b.*/
    for (I = 1; I <= Size; I++)
    {   
/* This step of the elimination is skipped if Temp equals zero. */
        if ((Temp = Intermediate[I]) != 0.0)
        {   pPivot = Matrix->Diag[I];
            Intermediate[I] = (Temp *= pPivot->Back);

            pElement = pPivot->NextInCol;
            while (pElement != NULL)
            {   Intermediate[pElement->Row] -= Temp * pElement->Back;
                pElement = pElement->NextInCol;
            }
        }
    }
    if (Matrix->Error == spZERO_DIAG  || Matrix->Error == spSINGULAR)
      return (Matrix->Error);
#ifdef SHARED_MEM
    timer_start(20);
#endif

/* Backward Substitution. Solves Ux = c.*/
    pExtOrder = &Matrix->IntToExtColMap[Size];
    for (I = Size; I > 0; I--)
    {   Temp = Intermediate[I];
        pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   Temp -= pElement->Back * Intermediate[pElement->Col];
            pElement = pElement->NextInRow;
        }
        Intermediate[I] = Temp;
        Solution[*(pExtOrder--)] = Intermediate[I];
    }

/* Unscramble Intermediate vector while placing data in to Solution vector. */
/*
    pExtOrder = &Matrix->IntToExtColMap[Size];
    for (I = Size; I > 0; I--)
        Solution[*(pExtOrder--)] = Intermediate[I];
*/

#ifdef SHARED_MEM
    timer_stop(20);
    total_solve_time = timer_get(20);
#endif
    return 0;
#endif /* REAL */
}
#endif










#if spCOMPLEX
/*
 *  SOLVE COMPLEX MATRIX EQUATION
 *
 *  Performs forward elimination and back substitution to find the
 *  unknown vector from the RHS vector and factored matrix.  This
 *  routine assumes that the pivots are associated with the lower
 *  triangular (L) matrix and that the diagonal of the upper triangular
 *  (U) matrix consists of ones.  This routine arranges the computation
 *  in different way than is traditionally used in order to exploit the
 *  sparsity of the right-hand side.  See the reference in spRevision.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  RHS  <input>  (RealVector)
 *      RHS is the real portion of the input data array, the right hand
 *      side. This data is undisturbed and may be reused for other solves.
 *  Solution  <output>  (RealVector)
 *      Solution is the real portion of the output data array. This routine
 *      is constructed such that RHS and Solution can be the same
 *      array.
 *  iRHS  <input>  (RealVector)
 *      iRHS is the imaginary portion of the input data array, the right
 *      hand side. This data is undisturbed and may be reused for other solves.
 *      If spSEPARATED_COMPLEX_VECTOR is set false, there is no need to
 *      supply this array.
 *  iSolution  <output>  (RealVector)
 *      iSolution is the imaginary portion of the output data array. This
 *      routine is constructed such that iRHS and iSolution can be
 *      the same array.  If spSEPARATED_COMPLEX_VECTOR is set false, there is no
 *      need to supply this array.
 *
 *  >>> Local variables:
 *  Intermediate  (ComplexVector)
 *      Temporary storage for use in forward elimination and backward
 *      substitution.  Commonly referred to as c, when the LU factorization
 *      equations are given as  Ax = b, Lc = b, Ux = c.
 *      Local version of Matrix->Intermediate, which was created during
 *      the initial factorization in function spcCreateInternalVectors() in the
 *      matrix factorization module.
 *  pElement  (ElementPtr)
 *      Pointer used to address elements in both the lower and upper triangle
 *      matrices.
 *  pExtOrder  (int *)
 *      Pointer used to sequentially access each entry in IntToExtRowMap
 *      and IntToExtColMap arrays.  Used to quickly scramble and unscramble
 *      RHS and Solution to account for row and column interchanges.
 *  pPivot  (ElementPtr)
 *      Pointer that points to current pivot or diagonal element.
 *  Size  (int)
 *      Size of matrix. Made local to reduce indirection.
 *  Temp  (ComplexNumber)
 *      Temporary storage for entries in arrays.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

static void
SolveComplexMatrix( Matrix, RHS, Solution IMAG_VECTORS )

MatrixPtr  Matrix;
RealVector  RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  ComplexVector  Intermediate;
register  int  I, *pExtOrder, Size;
ElementPtr  pPivot;
register ComplexVector  ExtVector;
ComplexNumber  Temp;

/* Begin `SolveComplexMatrix'. */

    Size = Matrix->Size;
    Intermediate = (ComplexVector)Matrix->Intermediate;

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
    --RHS;      --iRHS;
    --Solution; --iSolution;
#else
    RHS -= 2; Solution -= 2;
#endif
#endif

/* Initialize Intermediate vector. */
    pExtOrder = &Matrix->IntToExtRowMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Size; I > 0; I--)
    {   Intermediate[I].Real = RHS[*(pExtOrder)];
        Intermediate[I].Imag = iRHS[*(pExtOrder--)];
    }
#else
    ExtVector = (ComplexVector)RHS;
    for (I = Size; I > 0; I--)
        Intermediate[I] = ExtVector[*(pExtOrder--)];
#endif

/* Forward substitution. Solves Lc = b.*/
    for (I = 1; I <= Size; I++)
    {   Temp = Intermediate[I];

/* This step of the substitution is skipped if Temp equals zero. */
        if ((Temp.Real != 0.0) OR (Temp.Imag != 0.0))
        {   pPivot = Matrix->Diag[I];
/* Cmplx expr: Temp *= (1.0 / Pivot). */
            CMPLX_MULT_ASSIGN(Temp, *pPivot);
            Intermediate[I] = Temp;
            pElement = pPivot->NextInCol;
            while (pElement != NULL)
            {
/* Cmplx expr: Intermediate[Element->Row] -= Temp * *Element. */
                CMPLX_MULT_SUBT_ASSIGN(Intermediate[pElement->Row],
                                       Temp, *pElement);
                pElement = pElement->NextInCol;
            }
        }
    }

/* Backward Substitution. Solves Ux = c.*/
    for (I = Size; I > 0; I--)
    {   Temp = Intermediate[I];
        pElement = Matrix->Diag[I]->NextInRow;

        while (pElement != NULL)
        {
/* Cmplx expr: Temp -= *Element * Intermediate[Element->Col]. */
            CMPLX_MULT_SUBT_ASSIGN(Temp, *pElement,Intermediate[pElement->Col]);
            pElement = pElement->NextInRow;
        }
        Intermediate[I] = Temp;
    }

/* Unscramble Intermediate vector while placing data in to Solution vector. */
    pExtOrder = &Matrix->IntToExtColMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Size; I > 0; I--)
    {   Solution[*(pExtOrder)] = Intermediate[I].Real;
        iSolution[*(pExtOrder--)] = Intermediate[I].Imag;
    }
#else
    ExtVector = (ComplexVector)Solution;
    for (I = Size; I > 0; I--)
        ExtVector[*(pExtOrder--)] = Intermediate[I];
#endif

    return;
}
#endif /* spCOMPLEX */














#if TRANSPOSE
/*
 *  SOLVE TRANSPOSED MATRIX EQUATION
 *
 *  Performs forward elimination and back substitution to find the
 *  unknown vector from the RHS vector and transposed factored
 *  matrix. This routine is useful when performing sensitivity analysis
 *  on a circuit using the adjoint method.  This routine assumes that
 *  the pivots are associated with the untransposed lower triangular
 *  (L) matrix and that the diagonal of the untransposed upper
 *  triangular (U) matrix consists of ones.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  RHS  <input>  (RealVector)
 *      RHS is the input data array, the right hand side. This data is
 *      undisturbed and may be reused for other solves.
 *  Solution  <output>  (RealVector)
 *      Solution is the output data array. This routine is constructed such that
 *      RHS and Solution can be the same array.
 *  iRHS  <input>  (RealVector)
 *      iRHS is the imaginary portion of the input data array, the right
 *      hand side. This data is undisturbed and may be reused for other solves.
 *      If spSEPARATED_COMPLEX_VECTOR is set false, or if matrix is real, there
 *      is no need to supply this array.
 *  iSolution  <output>  (RealVector)
 *      iSolution is the imaginary portion of the output data array. This
 *      routine is constructed such that iRHS and iSolution can be
 *      the same array.  If spSEPARATED_COMPLEX_VECTOR is set false, or if
 *      matrix is real, there is no need to supply this array.
 *
 *  >>> Local variables:
 *  Intermediate  (RealVector)
 *      Temporary storage for use in forward elimination and backward
 *      substitution.  Commonly referred to as c, when the LU factorization
 *      equations are given as  Ax = b, Lc = b, Ux = c.  Local version of
 *      Matrix->Intermediate, which was created during the initial
 *      factorization in function spcCreateInternalVectors() in the matrix
 *      factorization module.
 *  pElement  (ElementPtr)
 *      Pointer used to address elements in both the lower and upper triangle
 *      matrices.
 *  pExtOrder  (int *)
 *      Pointer used to sequentially access each entry in IntToExtRowMap
 *      and IntToExtRowMap arrays.  Used to quickly scramble and unscramble
 *      RHS and Solution to account for row and column interchanges.
 *  pPivot  (ElementPtr)
 *      Pointer that points to current pivot or diagonal element.
 *  Size  (int)
 *      Size of matrix. Made local to reduce indirection.
 *  Temp  (RealNumber)
 *      Temporary storage for entries in arrays.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

/*VARARGS3*/

void
spSolveTransposed( eMatrix, RHS, Solution IMAG_VECTORS )

char *eMatrix;
RealVector  RHS, Solution IMAG_VECTORS;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement;
register  RealVector  Intermediate;
register  int  I, *pExtOrder, Size;
ElementPtr  pPivot;
RealNumber  Temp;
void SolveComplexTransposedMatrix();

/* Begin `spSolveTransposed'. */
    spExpandFormat(Matrix);
    ASSERT( IS_VALID(Matrix) AND IS_FACTORED(Matrix) );

#if spCOMPLEX
    if (Matrix->Complex)
    {   SolveComplexTransposedMatrix( Matrix, RHS, Solution IMAG_VECTORS );
        return;
    }
#endif

#if REAL
    Size = Matrix->Size;
    Intermediate = Matrix->Intermediate;

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
    --RHS;
    --Solution;
#endif

/* Initialize Intermediate vector. */
    pExtOrder = &Matrix->IntToExtColMap[Size];
    for (I = Size; I > 0; I--)
        Intermediate[I] = RHS[*(pExtOrder--)];

/* Forward elimination. */
    for (I = 1; I <= Size; I++)
    {   
/* This step of the elimination is skipped if Temp equals zero. */
        if ((Temp = Intermediate[I]) != 0.0)
        {   pElement = Matrix->Diag[I]->NextInRow;
            while (pElement != NULL)
            {   Intermediate[pElement->Col] -= Temp * pElement->Real;
                pElement = pElement->NextInRow;
            }

        }
    }

/* Backward Substitution. */
    for (I = Size; I > 0; I--)
    {   pPivot = Matrix->Diag[I];
        Temp = Intermediate[I];
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   Temp -= pElement->Real * Intermediate[pElement->Row];
            pElement = pElement->NextInCol;
        }
        Intermediate[I] = Temp * pPivot->Real;
    }

/* Unscramble Intermediate vector while placing data in to Solution vector. */
    pExtOrder = &Matrix->IntToExtRowMap[Size];
    for (I = Size; I > 0; I--)
        Solution[*(pExtOrder--)] = Intermediate[I];

    return;
#endif /* REAL */
}
#endif /* TRANSPOSE */










#if TRANSPOSE AND spCOMPLEX
/*
 *  SOLVE COMPLEX TRANSPOSED MATRIX EQUATION
 *
 *  Performs forward elimination and back substitution to find the
 *  unknown vector from the RHS vector and transposed factored
 *  matrix. This routine is useful when performing sensitivity analysis
 *  on a circuit using the adjoint method.  This routine assumes that
 *  the pivots are associated with the untransposed lower triangular
 *  (L) matrix and that the diagonal of the untransposed upper
 *  triangular (U) matrix consists of ones.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to matrix.
 *  RHS  <input>  (RealVector)
 *      RHS is the input data array, the right hand
 *      side. This data is undisturbed and may be reused for other solves.
 *      This vector is only the real portion if the matrix is complex and
 *      spSEPARATED_COMPLEX_VECTORS is set true.
 *  Solution  <output>  (RealVector)
 *      Solution is the real portion of the output data array. This routine
 *      is constructed such that RHS and Solution can be the same array.
 *      This vector is only the real portion if the matrix is complex and
 *      spSEPARATED_COMPLEX_VECTORS is set true.
 *  iRHS  <input>  (RealVector)
 *      iRHS is the imaginary portion of the input data array, the right
 *      hand side. This data is undisturbed and may be reused for other solves.
 *      If either spCOMPLEX or spSEPARATED_COMPLEX_VECTOR is set false, there
 *      is no need to supply this array.
 *  iSolution  <output>  (RealVector)
 *      iSolution is the imaginary portion of the output data array. This
 *      routine is constructed such that iRHS and iSolution can be
 *      the same array.  If spCOMPLEX or spSEPARATED_COMPLEX_VECTOR is set
 *      false, there is no need to supply this array.
 *
 *  >>> Local variables:
 *  Intermediate  (ComplexVector)
 *      Temporary storage for use in forward elimination and backward
 *      substitution.  Commonly referred to as c, when the LU factorization
 *      equations are given as  Ax = b, Lc = b, Ux = c.  Local version of
 *      Matrix->Intermediate, which was created during
 *      the initial factorization in function spcCreateInternalVectors() in the
 *      matrix factorization module.
 *  pElement  (ElementPtr)
 *      Pointer used to address elements in both the lower and upper triangle
 *      matrices.
 *  pExtOrder  (int *)
 *      Pointer used to sequentially access each entry in IntToExtRowMap
 *      and IntToExtColMap arrays.  Used to quickly scramble and unscramble
 *      RHS and Solution to account for row and column interchanges.
 *  pPivot  (ElementPtr)
 *      Pointer that points to current pivot or diagonal element.
 *  Size  (int)
 *      Size of matrix. Made local to reduce indirection.
 *  Temp  (ComplexNumber)
 *      Temporary storage for entries in arrays.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

static void
SolveComplexTransposedMatrix(Matrix, RHS, Solution IMAG_VECTORS )

MatrixPtr  Matrix;
RealVector  RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  ComplexVector  Intermediate;
register  int  I, *pExtOrder, Size;
register  ComplexVector  ExtVector;
ElementPtr  pPivot;
ComplexNumber  Temp;

/* Begin `SolveComplexTransposedMatrix'. */

    Size = Matrix->Size;
    Intermediate = (ComplexVector)Matrix->Intermediate;

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
    --RHS;      --iRHS;
    --Solution; --iSolution;
#else
    RHS -= 2;   Solution -= 2;
#endif
#endif

/* Initialize Intermediate vector. */
    pExtOrder = &Matrix->IntToExtColMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Size; I > 0; I--)
    {   Intermediate[I].Real = RHS[*(pExtOrder)];
        Intermediate[I].Imag = iRHS[*(pExtOrder--)];
    }
#else
    ExtVector = (ComplexVector)RHS;
    for (I = Size; I > 0; I--)
        Intermediate[I] = ExtVector[*(pExtOrder--)];
#endif

/* Forward elimination. */
    for (I = 1; I <= Size; I++)
    {   Temp = Intermediate[I];

/* This step of the elimination is skipped if Temp equals zero. */
        if ((Temp.Real != 0.0) OR (Temp.Imag != 0.0))
        {   pElement = Matrix->Diag[I]->NextInRow;
            while (pElement != NULL)
            {
/* Cmplx expr: Intermediate[Element->Col] -= Temp * *Element. */
                CMPLX_MULT_SUBT_ASSIGN( Intermediate[pElement->Col],
                                        Temp, *pElement);
                pElement = pElement->NextInRow;
            }
        }
    }

/* Backward Substitution. */
    for (I = Size; I > 0; I--)
    {   pPivot = Matrix->Diag[I];
        Temp = Intermediate[I];
        pElement = pPivot->NextInCol;

        while (pElement != NULL)
        {
/* Cmplx expr: Temp -= Intermediate[Element->Row] * *Element. */
            CMPLX_MULT_SUBT_ASSIGN(Temp,Intermediate[pElement->Row],*pElement);

            pElement = pElement->NextInCol;
        }
/* Cmplx expr: Intermediate = Temp * (1.0 / *pPivot). */
        CMPLX_MULT(Intermediate[I], Temp, *pPivot);
    }

/* Unscramble Intermediate vector while placing data in to Solution vector. */
    pExtOrder = &Matrix->IntToExtRowMap[Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Size; I > 0; I--)
    {   Solution[*(pExtOrder)] = Intermediate[I].Real;
        iSolution[*(pExtOrder--)] = Intermediate[I].Imag;
    }
#else
    ExtVector = (ComplexVector)Solution;
    for (I = Size; I > 0; I--)
        ExtVector[*(pExtOrder--)] = Intermediate[I];
#endif

    return;
}
#endif /* TRANSPOSE AND spCOMPLEX */

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
 *  MATRIX UTILITY MODULE
 *
 *  Author:                     Advising professor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains various optional utility routines.
 *
 *  >>> User accessible functions contained in this file:
 *  spMNA_Preorder
 *  spScale
 *  spMultiply
 *  spMultTransposed
 *  spDeterminant
 *  spStripFills
 *  spDeleteRowAndCol
 *  spPseudoCondition
 *  spCondition
 *  spNorm
 *  spLargestElement
 *  spRoundoff
 *  spErrorMessage
 *
 *  >>> Other functions contained in this file:
 *  CountTwins
 *  SwapCols
 *  ScaleComplexMatrix
 *  ComplexMatrixMultiply
 *  ComplexCondition
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



/*
 *  IMPORTS
 *
 *  >>> Import descriptions:
 *  spConfig.h
 *      Macros that customize the sparse matrix routines.
 *  spMatrix.h
 *      Macros and declarations to be imported by the user.
 *  spDefs.h
 *      Matrix type and macro definitions for the sparse matrix routines.
 */

/*
#define CHECK_IND
*/
#define spINSIDE_SPARSE
#include "spconfig.h"
#include "spmatrix.h"
#include "spdefs.h"

extern ElementPtr *returned_elements;
extern int *num_returned_elements;
void spExpandFormat(MatrixPtr);

/*
 *  Function declarations
 */

static int CountTwins( MatrixPtr, int, ElementPtr*, ElementPtr* );
static void SwapCols( MatrixPtr, ElementPtr, ElementPtr );
static void ScaleComplexMatrix( MatrixPtr, RealVector, RealVector );
#if spSEPARATED_COMPLEX_VECTORS
static void ComplexMatrixMultiply( MatrixPtr, RealVector, RealVector,
                                              RealVector, RealVector );
static void ComplexTransposedMatrixMultiply( MatrixPtr, RealVector, RealVector,
                                                        RealVector, RealVector);
#else
static void ComplexMatrixMultiply( MatrixPtr, RealVector, RealVector );
static void ComplexTransposedMatrixMultiply(MatrixPtr, RealVector, RealVector);
#endif
static RealNumber ComplexCondition( MatrixPtr, RealNumber, int* );
int f_ind (MatrixPtr, int, int);
void spCheckInd (MatrixPtr, char *);
MatrixPtr Matrix;

#ifdef DEBUG
#undef DEBUG
#endif
#define DEBUG

void spCheckInd ( MatrixPtr Matrix , char *msg)
{
    int I, J, K, Size = Matrix->Size;
    int err, any_err;
    ElementPtr pElement;

#ifdef DEBUG
#define Col_debug -1
#define Row_debug -1
    any_err = 0;
    for (I=1 ; I<=Size ; I++) {
      err = 0;
      K = 0;
      if (Matrix->FirstInCol[I]) {
        for (J=0 ; J<Matrix->Indsize ; J++) {
          while (f_ind(Matrix, I, K) < J) K++;
          if (pElement = Matrix->Col_fast[I][J]) {
            if (pElement->Col != I) {
              printf ("Col_fast pointing to bad Column number (%d), correct = %d\n",
                       pElement->Col, I);
              err = 1;
            }
            if (pElement->Row >= K) {
              printf ("Col_fast pointing to too high at Row: (%d,%d), %d >= %d\n",
                       I,J,pElement->Row,K);
              err = 1;
            }
            else {
              if (pElement->NextInCol && pElement->NextInCol->Row < K) {
                printf ("Next element in Col has Row too low: (%d,%d) %d < %d\n",
                         I,J,pElement->NextInCol->Row,K);
                err = 1;
              }
              else {
/*   Looks OK */
              }
            }
          }
          else {
            if (Matrix->FirstInCol[I]->Row < K) {
              printf ("Null Col_fast pointer at: (%d,%d)\n",I,J);
              err = 1;
            }
          }
        }
      }
      if (err | I == Col_debug) {
        if (err)
          any_err = 1;
        K = 0;
        for (J=0 ; J<Matrix->Indsize ; J++) {
          while (f_ind(Matrix, I, K) < J) K++;
          if (Matrix->Col_fast[I][J])
            printf ("%4d :: %4d %4d\n",J,Matrix->Col_fast[I][J]->Row,K);
          else
            printf ("%4d :: NULL %4d\n",J,K);
        }
        for (pElement = Matrix->FirstInCol[I] ; pElement ; pElement = pElement->NextInCol) {
          printf ("Element in Row: %4d\n",pElement->Row);
        }
        spColInd(Matrix,I);
      }
    }
    if (Matrix->RowsLinked) {
      for (I=1 ; I<=Size ; I++) {
        err = 0;
        K = 0;
        if (Matrix->FirstInRow[I]) {
          for (J=0 ; J<Matrix->Indsize ; J++) {
            while (f_ind(Matrix, I, K) < J) K++;
            if (pElement = Matrix->Row_fast[I][J]) {
              if (pElement->Row != I) {
                printf ("Row_fast pointing to bad Row number (%d), correct = %d\n",
                         pElement->Row, I);
                err = 1;
              }
              if (pElement->Col >= K) {
                printf ("Row_fast pointing to too high at Col: (%d,%d), %d >= %d\n",
                         I,J,pElement->Col,K);
                err = 1;
              }
              else {
                if (pElement->NextInRow && pElement->NextInRow->Col < K) {
                  printf ("Next element in Row has Col too low: (%d,%d) %d < %d\n",
                           I,J,pElement->NextInRow->Col,K);
                  err = 1;
                }
                else {
/*     Looks OK */
                }
              }
            }
            else {
              if (Matrix->FirstInRow[I]->Col < K) {
                printf ("Null Row_fast pointer at: (%d,%d)\n",I,J);
                err = 1;
              }
            }
          }
        }
        if (err | I == Row_debug) {
          if (err)
            any_err = 1;
          K = 0;
          for (J=0 ; J<Matrix->Indsize ; J++) {
            while (f_ind(Matrix, I, K) < J) K++;
            if (Matrix->Row_fast[I][J])
              printf ("%4d :: %4d %4d\n",J,Matrix->Row_fast[I][J]->Col,K);
            else
              printf ("%4d :: NULL %4d\n",J,K);
          }
          for (pElement = Matrix->FirstInRow[I] ; pElement ; pElement = pElement->NextInRow) {
            printf ("Element in Col: %4d\n",pElement->Col);
          }
          spRowInd(Matrix,I);
        }
      }
    }
    if (any_err) {
      printf ("Error(s) found at tag: %s\n",msg);
    }
    else {
      printf ("NO Error found at tag: %s\n",msg);
    }
#endif

    return;
}





#if MODIFIED_NODAL
/*
 *  PREORDER MODIFIED NODE ADMITTANCE MATRIX TO REMOVE ZEROS FROM DIAGONAL
 *
 *  This routine massages modified node admittance matrices to remove
 *  zeros from the diagonal.  It takes advantage of the fact that the
 *  row and column associated with a zero diagonal usually have
 *  structural ones placed symmetricly.  This routine should be used
 *  only on modified node admittance matrices and should be executed
 *  after the matrix has been built but before the factorization
 *  begins.  It should be executed for the initial factorization only
 *  and should be executed before the rows have been linked.  Thus it
 *  should be run before using spScale(), spMultiply(),
 *  spDeleteRowAndCol(), or spNorm().
 *
 *  This routine exploits the fact that the structural ones are placed
 *  in the matrix in symmetric twins.  For example, the stamps for
 *  grounded and a floating voltage sources are
 *  grounded:              floating:
 *  [  x   x   1 ]         [  x   x   1 ]
 *  [  x   x     ]         [  x   x  -1 ]
 *  [  1         ]         [  1  -1     ]
 *  Notice for the grounded source, there is one set of twins, and for
 *  the floating, there are two sets.  We remove the zero from the diagonal
 *  by swapping the rows associated with a set of twins.  For example:
 *  grounded:              floating 1:            floating 2:
 *  [  1         ]         [  1  -1     ]         [  x   x   1 ]
 *  [  x   x     ]         [  x   x  -1 ]         [  1  -1     ]
 *  [  x   x   1 ]         [  x   x   1 ]         [  x   x  -1 ]
 *
 *  It is important to deal with any zero diagonals that only have one
 *  set of twins before dealing with those that have more than one because
 *  swapping row destroys the symmetry of any twins in the rows being
 *  swapped, which may limit future moves.  Consider
 *  [  x   x   1     ]
 *  [  x   x  -1   1 ]
 *  [  1  -1         ]
 *  [      1         ]
 *  There is one set of twins for diagonal 4 and two for diagonal 3.
 *  Dealing with diagonal 4 first requires swapping rows 2 and 4.
 *  [  x   x   1     ]
 *  [      1         ]
 *  [  1  -1         ]
 *  [  x   x  -1   1 ]
 *  We can now deal with diagonal 3 by swapping rows 1 and 3.
 *  [  1  -1         ]
 *  [      1         ]
 *  [  x   x   1     ]
 *  [  x   x  -1   1 ]
 *  And we are done, there are no zeros left on the diagonal.  However, if
 *  we originally dealt with diagonal 3 first, we could swap rows 2 and 3
 *  [  x   x   1     ]
 *  [  1  -1         ]
 *  [  x   x  -1   1 ]
 *  [      1         ]
 *  Diagonal 4 no longer has a symmetric twin and we cannot continue.
 *
 *  So we always take care of lone twins first.  When none remain, we
 *  choose arbitrarily a set of twins for a diagonal with more than one set
 *  and swap the rows corresponding to that twin.  We then deal with any
 *  lone twins that were created and repeat the procedure until no
 *  zero diagonals with symmetric twins remain.
 *
 *  In this particular implementation, columns are swapped rather than rows.
 *  The algorithm used in this function was developed by Ken Kundert and
 *  Tom Quarles.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix to be preordered.
 *
 *  >>> Local variables;
 *  J  (int)
 *      Column with zero diagonal being currently considered.
 *  pTwin1  (ElementPtr)
 *      Pointer to the twin found in the column belonging to the zero diagonal.
 *  pTwin2  (ElementPtr)
 *      Pointer to the twin found in the row belonging to the zero diagonal.
 *      belonging to the zero diagonal.
 *  AnotherPassNeeded  (BOOLEAN)
 *      Flag indicating that at least one zero diagonal with symmetric twins
 *      remain.
 *  StartAt  (int)
 *      Column number of first zero diagonal with symmetric twins.
 *  Swapped  (BOOLEAN)
 *      Flag indicating that columns were swapped on this pass.
 *  Twins  (int)
 *      Number of symmetric twins corresponding to current zero diagonal.
 */

void
spMNA_Preorder( eMatrix )

char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  int  J, Size;
ElementPtr  pTwin1, pTwin2;
int  Twins, StartAt = 1;
BOOLEAN  Swapped, AnotherPassNeeded;

/* Begin `spMNA_Preorder'. */
    spExpandFormat(Matrix);
    ASSERT( IS_VALID(Matrix) AND NOT Matrix->Factored );

    if (Matrix->RowsLinked) return;
    Size = Matrix->Size;
    Matrix->Reordered = YES;

    do
    {   AnotherPassNeeded = Swapped = NO;

/* Search for zero diagonals with lone twins. */
        for (J = StartAt; J <= Size; J++)
        {   if (Matrix->Diag[J] == NULL)
            {   Twins = CountTwins( Matrix, J, &pTwin1, &pTwin2 );
                if (Twins == 1)
                {   /* Lone twins found, swap rows. */
                    SwapCols( Matrix, pTwin1, pTwin2 );
                    Swapped = YES;
                }
                else if ((Twins > 1) AND NOT AnotherPassNeeded)
                {   AnotherPassNeeded = YES;
                    StartAt = J;
                }
            }
        }

/* All lone twins are gone, look for zero diagonals with multiple twins. */
        if (AnotherPassNeeded)
        {   for (J = StartAt; NOT Swapped AND (J <= Size); J++)
            {   if (Matrix->Diag[J] == NULL)
                {   Twins = CountTwins( Matrix, J, &pTwin1, &pTwin2 );
                    SwapCols( Matrix, pTwin1, pTwin2 );
                    Swapped = YES;
                }
            }
        }
    } while (AnotherPassNeeded);
    return;
}




/*
 *  COUNT TWINS
 *
 *  This function counts the number of symmetric twins associated with
 *  a zero diagonal and returns one set of twins if any exist.  The
 *  count is terminated early at two.
 */

static int
CountTwins( Matrix, Col, ppTwin1, ppTwin2 )

MatrixPtr Matrix;
int Col;
ElementPtr *ppTwin1, *ppTwin2;
{
int Row, Twins = 0;
ElementPtr pTwin1, pTwin2;

/* Begin `CountTwins'. */

    pTwin1 = Matrix->FirstInCol[Col];
    while (pTwin1 != NULL)
    {   if (ABS(pTwin1->Real) == 1.0)
        {   Row = pTwin1->Row;
            pTwin2 = Matrix->FirstInCol[Row];
            while ((pTwin2 != NULL) AND (pTwin2->Row != Col))
                pTwin2 = pTwin2->NextInCol;
            if ((pTwin2 != NULL) AND (ABS(pTwin2->Real) == 1.0))
            {   /* Found symmetric twins. */
                if (++Twins >= 2) return Twins;
#ifdef notdef
                (*ppTwin1 = pTwin1)/*->Col = Col XXX */;
                (*ppTwin2 = pTwin2)/*->Col = Row XXX */;
#else
                (*ppTwin1 = pTwin1)->Col = Col;
                (*ppTwin2 = pTwin2)->Col = Row;
#endif
            }
        }
        pTwin1 = pTwin1->NextInCol;
    }
    return Twins;
}




/*
 *  SWAP COLUMNS
 *
 *  This function swaps two columns and is applicable before the rows are
 *  linked.
 */

static void
SwapCols( Matrix, pTwin1, pTwin2 )

MatrixPtr Matrix;
ElementPtr pTwin1, pTwin2;
{
int Col1 = pTwin1->Col, Col2 = pTwin2->Col;

/* Begin `SwapCols'. */

#ifdef notdef
ElementPtr e; /*XXX*/
    /* XXX Update column numbers */
    for (e = Matrix->FirstInCol[Col1]; e != NULL; e = e->NextInCol)
	e->Col = Col2;
    for (e = Matrix->FirstInCol[Col2]; e != NULL; e = e->NextInCol)
	e->Col = Col1;
#endif

    SWAP (ElementPtr, Matrix->FirstInCol[Col1], Matrix->FirstInCol[Col2]);
    SWAP (int, Matrix->IntToExtColMap[Col1], Matrix->IntToExtColMap[Col2]);
#if TRANSLATE
    Matrix->ExtToIntColMap[Matrix->IntToExtColMap[Col2]]=Col2;
    Matrix->ExtToIntColMap[Matrix->IntToExtColMap[Col1]]=Col1;
#endif

    Matrix->Diag[Col1] = pTwin2;
    Matrix->Diag[Col2] = pTwin1;
    Matrix->NumberOfInterchangesIsOdd = NOT Matrix->NumberOfInterchangesIsOdd;
    return;
}
#endif /* MODIFIED_NODAL */









#if SCALING
/*
 *  SCALE MATRIX
 *
 *  This function scales the matrix to enhance the possibility of
 *  finding a good pivoting order.  Note that scaling enhances accuracy
 *  of the solution only if it affects the pivoting order, so it makes
 *  no sense to scale the matrix before spFactor().  If scaling is
 *  desired it should be done before spOrderAndFactor().  There
 *  are several things to take into account when choosing the scale
 *  factors.  First, the scale factors are directly multiplied against
 *  the elements in the matrix.  To prevent roundoff, each scale factor
 *  should be equal to an integer power of the number base of the
 *  machine.  Since most machines operate in base two, scale factors
 *  should be a power of two.  Second, the matrix should be scaled such
 *  that the matrix of element uncertainties is equilibrated.  Third,
 *  this function multiplies the scale factors by the elements, so if
 *  one row tends to have uncertainties 1000 times smaller than the
 *  other rows, then its scale factor should be 1024, not 1/1024.
 *  Fourth, to save time, this function does not scale rows or columns
 *  if their scale factors are equal to one.  Thus, the scale factors
 *  should be normalized to the most common scale factor.  Rows and
 *  columns should be normalized separately.  For example, if the size
 *  of the matrix is 100 and 10 rows tend to have uncertainties near
 *  1e-6 and the remaining 90 have uncertainties near 1e-12, then the
 *  scale factor for the 10 should be 1/1,048,576 and the scale factors
 *  for the remaining 90 should be 1.  Fifth, since this routine
 *  directly operates on the matrix, it is necessary to apply the scale
 *  factors to the RHS and Solution vectors.  It may be easier to
 *  simply use spOrderAndFactor() on a scaled matrix to choose the
 *  pivoting order, and then throw away the matrix.  Subsequent
 *  factorizations, performed with spFactor(), will not need to have
 *  the RHS and Solution vectors descaled.  Lastly, this function
 *  should not be executed before the function spMNA_Preorder.
 *
 *  >>> Arguments:
 *  eMatrix  <input> (char *)
 *      Pointer to the matrix to be scaled.
 *  SolutionScaleFactors  <input>  (RealVector)
 *      The array of Solution scale factors.  These factors scale the columns.
 *      All scale factors are real valued.
 *  RHS_ScaleFactors  <input>  (RealVector)
 *      The array of RHS scale factors.  These factors scale the rows.
 *      All scale factors are real valued.
 *
 *  >>> Local variables:
 *  lSize  (int)
 *      Local version of the size of the matrix.
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.
 *  pExtOrder  (int *)
 *      Pointer into either IntToExtRowMap or IntToExtColMap vector. Used to
 *      compensate for any row or column swaps that have been performed.
 *  ScaleFactor  (RealNumber)
 *      The scale factor being used on the current row or column.
 */

void
spScale( eMatrix, RHS_ScaleFactors, SolutionScaleFactors )

char *eMatrix;
register  RealVector  RHS_ScaleFactors, SolutionScaleFactors;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr  pElement;
register int  I, lSize, *pExtOrder;
RealNumber  ScaleFactor;
void ScaleComplexMatrix();

/* Begin `spScale'. */
    spExpandFormat(Matrix);
    ASSERT( IS_VALID(Matrix) AND NOT Matrix->Factored );
    if (NOT Matrix->RowsLinked) spcLinkRows( Matrix );

#if spCOMPLEX
    if (Matrix->Complex)
    {   ScaleComplexMatrix( Matrix, RHS_ScaleFactors, SolutionScaleFactors );
        return;
    }
#endif

#if REAL
    lSize = Matrix->Size;

/* Correct pointers to arrays for ARRAY_OFFSET */
#if NOT ARRAY_OFFSET
    --RHS_ScaleFactors;
    --SolutionScaleFactors;
#endif

/* Scale Rows */
    pExtOrder = &Matrix->IntToExtRowMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = RHS_ScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInRow[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement = pElement->NextInRow;
            }
        }
    }

/* Scale Columns */
    pExtOrder = &Matrix->IntToExtColMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = SolutionScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInCol[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement = pElement->NextInCol;
            }
        }
    }
    return;

#endif /* REAL */
}
#endif /* SCALING */


#ifdef CHILE
void spBackup(MatrixPtr  Matrix)
{
    int I, Size;
    ElementPtr pElement;

    Size = Matrix->Size;
    for (I = 1; I <= Size; I++) {
       pElement = Matrix->FirstInRow[I];
       while (pElement != NULL) {
          pElement->RealBackup = pElement->Real;
          pElement = pElement->NextInRow;
       }
    }
    return;
}
#endif






#if spCOMPLEX AND SCALING
/*
 *  SCALE COMPLEX MATRIX
 *
 *  This function scales the matrix to enhance the possibility of
 *  finding a good pivoting order.  Note that scaling enhances accuracy
 *  of the solution only if it affects the pivoting order, so it makes
 *  no sense to scale the matrix before spFactor().  If scaling is
 *  desired it should be done before spOrderAndFactor().  There
 *  are several things to take into account when choosing the scale
 *  factors.  First, the scale factors are directly multiplied against
 *  the elements in the matrix.  To prevent roundoff, each scale factor
 *  should be equal to an integer power of the number base of the
 *  machine.  Since most machines operate in base two, scale factors
 *  should be a power of two.  Second, the matrix should be scaled such
 *  that the matrix of element uncertainties is equilibrated.  Third,
 *  this function multiplies the scale factors by the elements, so if
 *  one row tends to have uncertainties 1000 times smaller than the
 *  other rows, then its scale factor should be 1024, not 1/1024.
 *  Fourth, to save time, this function does not scale rows or columns
 *  if their scale factors are equal to one.  Thus, the scale factors
 *  should be normalized to the most common scale factor.  Rows and
 *  columns should be normalized separately.  For example, if the size
 *  of the matrix is 100 and 10 rows tend to have uncertainties near
 *  1e-6 and the remaining 90 have uncertainties near 1e-12, then the
 *  scale factor for the 10 should be 1/1,048,576 and the scale factors
 *  for the remaining 90 should be 1. Fifth, since this routine
 *  directly operates on the matrix, it is necessary to apply the scale
 *  factors to the RHS and Solution vectors.  It may be easier to
 *  simply use spOrderAndFactor() on a scaled matrix to choose the
 *  pivoting order, and then throw away the matrix.  Subsequent
 *  factorizations, performed with spFactor(), will not need to have
 *  the RHS and Solution vectors descaled.  Lastly, this function
 *  should not be executed before the function spMNA_Preorder.
 *
 *  >>> Arguments:
 *  Matrix  <input> (char *)
 *      Pointer to the matrix to be scaled.
 *  SolutionScaleFactors  <input>  (RealVector)
 *      The array of Solution scale factors.  These factors scale the columns.
 *      All scale factors are real valued.
 *  RHS_ScaleFactors  <input>  (RealVector)
 *      The array of RHS scale factors.  These factors scale the rows.
 *      All scale factors are real valued.
 *
 *  >>> Local variables:
 *  lSize  (int)
 *      Local version of the size of the matrix.
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.
 *  pExtOrder  (int *)
 *      Pointer into either IntToExtRowMap or IntToExtColMap vector. Used to
 *      compensate for any row or column swaps that have been performed.
 *  ScaleFactor  (RealNumber)
 *      The scale factor being used on the current row or column.
 */

static void
ScaleComplexMatrix( Matrix, RHS_ScaleFactors, SolutionScaleFactors )

MatrixPtr  Matrix;
register  RealVector  RHS_ScaleFactors, SolutionScaleFactors;
{
register ElementPtr  pElement;
register int  I, lSize, *pExtOrder;
RealNumber  ScaleFactor;

/* Begin `ScaleComplexMatrix'. */
    lSize = Matrix->Size;

/* Correct pointers to arrays for ARRAY_OFFSET */
#if NOT ARRAY_OFFSET
    --RHS_ScaleFactors;
    --SolutionScaleFactors;
#endif

/* Scale Rows */
    pExtOrder = &Matrix->IntToExtRowMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = RHS_ScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInRow[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement->Imag *= ScaleFactor;
                pElement = pElement->NextInRow;
            }
        }
    }

/* Scale Columns */
    pExtOrder = &Matrix->IntToExtColMap[1];
    for (I = 1; I <= lSize; I++)
    {   if ((ScaleFactor = SolutionScaleFactors[*(pExtOrder++)]) != 1.0)
        {   pElement = Matrix->FirstInCol[I];

            while (pElement != NULL)
            {   pElement->Real *= ScaleFactor;
                pElement->Imag *= ScaleFactor;
                pElement = pElement->NextInCol;
            }
        }
    }
    return;
}
#endif /* SCALING AND spCOMPLEX */








#if MULTIPLICATION
/*
 *  MATRIX MULTIPLICATION
 *
 *  Multiplies matrix by solution vector to find source vector.
 *  Assumes matrix has not been factored.  This routine can be used
 *  as a test to see if solutions are correct.  It should not be used
 *  before spMNA_Preorder().
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 *  RHS  <output>  (RealVector)
 *      RHS is the right hand side. This is what is being solved for.
 *  Solution  <input>  (RealVector)
 *      Solution is the vector being multiplied by the matrix.
 *  iRHS  <output>  (RealVector)
 *      iRHS is the imaginary portion of the right hand side. This is
 *      what is being solved for.  This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *  iSolution  <input>  (RealVector)
 *      iSolution is the imaginary portion of the vector being multiplied
 *      by the matrix. This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

void
spMultiply( eMatrix, RHS, Solution IMAG_VECTORS )

char *eMatrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  RealVector  Vector;
register  RealNumber  Sum;
register  int  I, *pExtOrder;
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
extern void ComplexMatrixMultiply();

/* Begin `spMultiply'. */
    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE( Matrix ) AND NOT Matrix->Factored );
    if (NOT Matrix->RowsLinked)
	spcLinkRows(Matrix);
    if (NOT Matrix->InternalVectorsAllocated)
	spcCreateInternalVectors( Matrix );

#if spCOMPLEX
    if (Matrix->Complex)
    {   ComplexMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS );
        return;
    }
#endif

#if REAL
#if NOT ARRAY_OFFSET
/* Correct array pointers for ARRAY_OFFSET. */
    --RHS;
    --Solution;
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = Solution[*(pExtOrder--)];

    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInRow[I];
        Sum = 0.0;

        while (pElement != NULL)
        {   Sum += pElement->Real * Vector[pElement->Col];
            pElement = pElement->NextInRow;
        }
        RHS[*pExtOrder--] = Sum;
    }
    return;
#endif /* REAL */
}
#endif /* MULTIPLICATION */







#if spCOMPLEX AND MULTIPLICATION
/*
 *  COMPLEX MATRIX MULTIPLICATION
 *
 *  Multiplies matrix by solution vector to find source vector.
 *  Assumes matrix has not been factored.  This routine can be  used
 *  as a test to see if solutions are correct.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to the matrix.
 *  RHS  <output>  (RealVector)
 *      RHS is the right hand side. This is what is being solved for.
 *      This is only the real portion of the right-hand side if the matrix
 *      is complex and spSEPARATED_COMPLEX_VECTORS is set true.
 *  Solution  <input>  (RealVector)
 *      Solution is the vector being multiplied by the matrix. This is only
 *      the real portion if the matrix is complex and
 *      spSEPARATED_COMPLEX_VECTORS is set true.
 *  iRHS  <output>  (RealVector)
 *      iRHS is the imaginary portion of the right hand side. This is
 *      what is being solved for.  This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *  iSolution  <input>  (RealVector)
 *      iSolution is the imaginary portion of the vector being multiplied
 *      by the matrix. This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

static void
ComplexMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS )

MatrixPtr  Matrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  ComplexVector  Vector;
ComplexNumber  Sum;
register  int  I, *pExtOrder;

/* Begin `ComplexMatrixMultiply'. */

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
    --RHS;              --iRHS;
    --Solution;         --iSolution;
#else
    RHS -= 2;           Solution -= 2;
#endif
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = (ComplexVector)Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Matrix->Size; I > 0; I--)
    {   Vector[I].Real = Solution[*pExtOrder];
        Vector[I].Imag = iSolution[*(pExtOrder--)];
    }
#else
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = ((ComplexVector)Solution)[*(pExtOrder--)];
#endif

    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInRow[I];
        Sum.Real = Sum.Imag = 0.0;

        while (pElement != NULL)
        {   /* Cmplx expression : Sum += Element * Vector[Col] */
            CMPLX_MULT_ADD_ASSIGN( Sum, *pElement, Vector[pElement->Col] );
            pElement = pElement->NextInRow;
        }

#if spSEPARATED_COMPLEX_VECTORS
        RHS[*pExtOrder] = Sum.Real;
        iRHS[*pExtOrder--] = Sum.Imag;
#else
        ((ComplexVector)RHS)[*pExtOrder--] = Sum;
#endif
    }
    return;
}
#endif /* spCOMPLEX AND MULTIPLICATION */








#if MULTIPLICATION AND TRANSPOSE
/*
 *  TRANSPOSED MATRIX MULTIPLICATION
 *
 *  Multiplies transposed matrix by solution vector to find source vector.
 *  Assumes matrix has not been factored.  This routine can be used
 *  as a test to see if solutions are correct.  It should not be used
 *  before spMNA_Preorder().
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 *  RHS  <output>  (RealVector)
 *      RHS is the right hand side. This is what is being solved for.
 *  Solution  <input>  (RealVector)
 *      Solution is the vector being multiplied by the matrix.
 *  iRHS  <output>  (RealVector)
 *      iRHS is the imaginary portion of the right hand side. This is
 *      what is being solved for.  This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *  iSolution  <input>  (RealVector)
 *      iSolution is the imaginary portion of the vector being multiplied
 *      by the matrix. This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

void
spMultTransposed( eMatrix, RHS, Solution IMAG_VECTORS )

char *eMatrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  RealVector  Vector;
register  RealNumber  Sum;
register  int  I, *pExtOrder;
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
extern void ComplexTransposedMatrixMultiply();

/* Begin `spMultTransposed'. */
    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE( Matrix ) AND NOT Matrix->Factored );
    if (NOT Matrix->InternalVectorsAllocated)
	spcCreateInternalVectors( Matrix );

#if spCOMPLEX
    if (Matrix->Complex)
    {   ComplexTransposedMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS );
        return;
    }
#endif

#if REAL
#if NOT ARRAY_OFFSET
/* Correct array pointers for ARRAY_OFFSET. */
    --RHS;
    --Solution;
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = Solution[*(pExtOrder--)];

    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInCol[I];
        Sum = 0.0;

        while (pElement != NULL)
        {   Sum += pElement->Real * Vector[pElement->Row];
            pElement = pElement->NextInCol;
        }
        RHS[*pExtOrder--] = Sum;
    }
    return;
#endif /* REAL */
}
#endif /* MULTIPLICATION AND TRANSPOSE */







#if spCOMPLEX AND MULTIPLICATION AND TRANSPOSE
/*
 *  COMPLEX TRANSPOSED MATRIX MULTIPLICATION
 *
 *  Multiplies transposed matrix by solution vector to find source vector.
 *  Assumes matrix has not been factored.  This routine can be  used
 *  as a test to see if solutions are correct.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to the matrix.
 *  RHS  <output>  (RealVector)
 *      RHS is the right hand side. This is what is being solved for.
 *      This is only the real portion of the right-hand side if the matrix
 *      is complex and spSEPARATED_COMPLEX_VECTORS is set true.
 *  Solution  <input>  (RealVector)
 *      Solution is the vector being multiplied by the matrix. This is only
 *      the real portion if the matrix is complex and
 *      spSEPARATED_COMPLEX_VECTORS is set true.
 *  iRHS  <output>  (RealVector)
 *      iRHS is the imaginary portion of the right hand side. This is
 *      what is being solved for.  This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *  iSolution  <input>  (RealVector)
 *      iSolution is the imaginary portion of the vector being multiplied
 *      by the matrix. This is only necessary if the matrix is
 *      complex and spSEPARATED_COMPLEX_VECTORS is true.
 *
 *  >>> Obscure Macros
 *  IMAG_VECTORS
 *      Replaces itself with `, iRHS, iSolution' if the options spCOMPLEX and
 *      spSEPARATED_COMPLEX_VECTORS are set, otherwise it disappears
 *      without a trace.
 */

static void
ComplexTransposedMatrixMultiply( Matrix, RHS, Solution IMAG_VECTORS )

MatrixPtr  Matrix;
RealVector RHS, Solution IMAG_VECTORS;
{
register  ElementPtr  pElement;
register  ComplexVector  Vector;
ComplexNumber  Sum;
register  int  I, *pExtOrder;

/* Begin `ComplexMatrixMultiply'. */

/* Correct array pointers for ARRAY_OFFSET. */
#if NOT ARRAY_OFFSET
#if spSEPARATED_COMPLEX_VECTORS
    --RHS;              --iRHS;
    --Solution;         --iSolution;
#else
    RHS -= 2;           Solution -= 2;
#endif
#endif

/* Initialize Intermediate vector with reordered Solution vector. */
    Vector = (ComplexVector)Matrix->Intermediate;
    pExtOrder = &Matrix->IntToExtRowMap[Matrix->Size];

#if spSEPARATED_COMPLEX_VECTORS
    for (I = Matrix->Size; I > 0; I--)
    {   Vector[I].Real = Solution[*pExtOrder];
        Vector[I].Imag = iSolution[*(pExtOrder--)];
    }
#else
    for (I = Matrix->Size; I > 0; I--)
        Vector[I] = ((ComplexVector)Solution)[*(pExtOrder--)];
#endif

    pExtOrder = &Matrix->IntToExtColMap[Matrix->Size];
    for (I = Matrix->Size; I > 0; I--)
    {   pElement = Matrix->FirstInCol[I];
        Sum.Real = Sum.Imag = 0.0;

        while (pElement != NULL)
        {   /* Cmplx expression : Sum += Element * Vector[Row] */
            CMPLX_MULT_ADD_ASSIGN( Sum, *pElement, Vector[pElement->Row] );
            pElement = pElement->NextInCol;
        }

#if spSEPARATED_COMPLEX_VECTORS
        RHS[*pExtOrder] = Sum.Real;
        iRHS[*pExtOrder--] = Sum.Imag;
#else
        ((ComplexVector)RHS)[*pExtOrder--] = Sum;
#endif
    }
    return;
}
#endif /* spCOMPLEX AND MULTIPLICATION AND TRANSPOSE */








#if DETERMINANT
/*
 *  CALCULATE DETERMINANT
 *
 *  This routine in capable of calculating the determinant of the
 *  matrix once the LU factorization has been performed.  Hence, only
 *  use this routine after spFactor() and before spClear().
 *  The determinant equals the product of all the diagonal elements of
 *  the lower triangular matrix L, except that this product may need
 *  negating.  Whether the product or the negative product equals the
 *  determinant is determined by the number of row and column
 *  interchanges performed.  Note that the determinants of matrices can
 *  be very large or very small.  On large matrices, the determinant
 *  can be far larger or smaller than can be represented by a floating
 *  point number.  For this reason the determinant is scaled to a
 *  reasonable value and the logarithm of the scale factor is returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      A pointer to the matrix for which the determinant is desired.
 *  pExponent  <output>  (int *)
 *      The logarithm base 10 of the scale factor for the determinant.  To find
 *      the actual determinant, Exponent should be added to the exponent of
 *      Determinant.
 *  pDeterminant  <output>  (RealNumber *)
 *      The real portion of the determinant.   This number is scaled to be
 *      greater than or equal to 1.0 and less than 10.0.
 *  piDeterminant  <output>  (RealNumber *)
 *      The imaginary portion of the determinant.  When the matrix is real
 *      this pointer need not be supplied, nothing will be returned.   This
 *      number is scaled to be greater than or equal to 1.0 and less than 10.0.
 *
 *  >>> Local variables:
 *  Norm  (RealNumber)
 *      L-infinity norm of a complex number.
 *  Size  (int)
 *      Local storage for Matrix->Size.  Placed in a register for speed.
 *  Temp  (RealNumber)
 *      Temporary storage for real portion of determinant.
 */

#if spCOMPLEX
void
spDeterminant( eMatrix, pExponent, pDeterminant, piDeterminant )
RealNumber *piDeterminant;
#else
void
spDeterminant( eMatrix, pExponent, pDeterminant )
#endif

char *eMatrix;
register  RealNumber *pDeterminant;
int  *pExponent;
{
register MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register int I, Size;
RealNumber Norm, nr, ni;
ComplexNumber Pivot, cDeterminant;

#define  NORM(a)     (nr = ABS((a).Real), ni = ABS((a).Imag), MAX (nr,ni))

/* Begin `spDeterminant'. */
    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE( Matrix ) AND IS_FACTORED(Matrix) );
    *pExponent = 0;

    if (Matrix->Error == spSINGULAR)
    {   *pDeterminant = 0.0;
#if spCOMPLEX
        if (Matrix->Complex) *piDeterminant = 0.0;
#endif
        return;
    }

    Size = Matrix->Size;
    I = 0;

#if spCOMPLEX
    if (Matrix->Complex)        /* Complex Case. */
    {   cDeterminant.Real = 1.0;
        cDeterminant.Imag = 0.0;

        while (++I <= Size)
        {   CMPLX_RECIPROCAL( Pivot, *Matrix->Diag[I] );
            CMPLX_MULT_ASSIGN( cDeterminant, Pivot );

/* Scale Determinant. */
            Norm = NORM( cDeterminant );
            if (Norm != 0.0)
            {   while (Norm >= 1.0e12)
                {   cDeterminant.Real *= 1.0e-12;
                    cDeterminant.Imag *= 1.0e-12;
                    *pExponent += 12;
                    Norm = NORM( cDeterminant );
                }
                while (Norm < 1.0e-12)
                {   cDeterminant.Real *= 1.0e12;
                    cDeterminant.Imag *= 1.0e12;
                    *pExponent -= 12;
                    Norm = NORM( cDeterminant );
                }
            }
        }

/* Scale Determinant again, this time to be between 1.0 <= x < 10.0. */
        Norm = NORM( cDeterminant );
        if (Norm != 0.0)
        {   while (Norm >= 10.0)
            {   cDeterminant.Real *= 0.1;
                cDeterminant.Imag *= 0.1;
                (*pExponent)++;
                Norm = NORM( cDeterminant );
            }
            while (Norm < 1.0)
            {   cDeterminant.Real *= 10.0;
                cDeterminant.Imag *= 10.0;
                (*pExponent)--;
                Norm = NORM( cDeterminant );
            }
        }
        if (Matrix->NumberOfInterchangesIsOdd)
            CMPLX_NEGATE( cDeterminant );
        
        *pDeterminant = cDeterminant.Real;
        *piDeterminant = cDeterminant.Imag;
    }
#endif /* spCOMPLEX */
#if REAL AND spCOMPLEX
    else
#endif
#if REAL
    {   /* Real Case. */
        *pDeterminant = 1.0;

        while (++I <= Size)
        {   *pDeterminant /= Matrix->Diag[I]->Real;

/* Scale Determinant. */
            if (*pDeterminant != 0.0)
            {   while (ABS(*pDeterminant) >= 1.0e12)
                {   *pDeterminant *= 1.0e-12;
                    *pExponent += 12;
                }
                while (ABS(*pDeterminant) < 1.0e-12)
                {   *pDeterminant *= 1.0e12;
                    *pExponent -= 12;
                }
            }
        }

/* Scale Determinant again, this time to be between 1.0 <= x < 10.0. */
        if (*pDeterminant != 0.0)
        {   while (ABS(*pDeterminant) >= 10.0)
            {   *pDeterminant *= 0.1;
                (*pExponent)++;
            }
            while (ABS(*pDeterminant) < 1.0)
            {   *pDeterminant *= 10.0;
                (*pExponent)--;
            }
        }
        if (Matrix->NumberOfInterchangesIsOdd)
            *pDeterminant = -*pDeterminant;
    }
#endif /* REAL */
}
#endif /* DETERMINANT */









/*
 *  STRIP FILL-INS FROM MATRIX
 *
 *  Strips the matrix of all fill-ins.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to the matrix to be stripped.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      Pointer that is used to step through the matrix.
 *  ppElement  (ElementPtr *)
 *      Pointer to the location of an ElementPtr.  This location will be
 *      updated if a fill-in is stripped from the matrix.
 *  pFillin  (ElementPtr)
 *      Pointer used to step through the lists of fill-ins while marking them.
 *  pLastFillin  (ElementPtr)
 *      A pointer to the last fill-in in the list.  Used to terminate a loop.
 *  pListNode  (struct  FillinListNodeStruct *)
 *      A pointer to a node in the FillinList linked-list.
 */

void print_col(MatrixPtr,int), print_row(MatrixPtr,int);
void
spStripFills( char *eMatrix, int Step )
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
struct FillinListNodeStruct  *pListNode;
int fast_ind, fast_ind_done;

/* Begin `spStripFills'. */
    ASSERT( IS_SPARSE( Matrix ) );
    if (Matrix->Fillins == 0) return;
    Matrix->NeedsOrdering = YES;

#ifdef CHECK_IND
    spCheckInd(Matrix, "spStripFills top");
#endif
/* Unlink fill-ins by searching for elements marked with Fillin >= Step */
    {   register  ElementPtr pElement, *ppElement, last_el;
        register  int  I, Size = Matrix->Size;

/* Unlink fill-ins in all columns. */
        for (I = (Step >= 1) ? Step : 1; I <= Size; I++) {
            fast_ind_done = -1;
            last_el = NULL;
            ppElement = &(Matrix->FirstInCol[I]);
            while ((pElement = *ppElement) != NULL) {   
              if (pElement->Col != I) {
                printf ("Internal Error: Trouble on traversing Column (%d,%d) in spStripFills\n",
                  pElement->Col, I);
              }
              if (pElement->Fillin >= Step && pElement->Fillin != 0)
              {   *ppElement = pElement->NextInCol;  /* Unlink fill-in. */
                  if (Matrix->Diag[pElement->Col] == pElement)
                      Matrix->Diag[pElement->Col] = NULL;
                  Matrix->Fillins--;
              }
              else
              {
                  fast_ind = f_ind(Matrix, I, pElement->Row);
                  while (fast_ind_done < fast_ind) {
                    Matrix->Col_fast[I][++fast_ind_done] = last_el;
                  }
                  last_el = pElement;
                  ppElement = &pElement->NextInCol;  /* Skip element. */
              }
            }
            while (fast_ind_done < Matrix->Indsize)
              Matrix->Col_fast[I][f_ind(Matrix, I,++fast_ind_done)] = last_el;
        }

/* Unlink fill-ins in all rows. */
        for (I = 1; I <= Size; I++) {
            fast_ind_done = -1;
            last_el = NULL;
            ppElement = &(Matrix->FirstInRow[I]);
            while ((pElement = *ppElement) != NULL) {
              if (pElement->Row != I) {
                printf ("Internal Error: Trouble on traversing Row (%d,%d) in spStripFills\n",
                  pElement->Row,I);
              }
              if (pElement->Col >= Step && pElement->Fillin >= Step && pElement->Fillin != 0) {
                  *ppElement = pElement->NextInRow;  /* Unlink fill-in. */
                  pElement->NextInCol = returned_elements[pElement->Row];
                  returned_elements[pElement->Row] = pElement;
                  num_returned_elements[pElement->Row]++;
              }
              else {
                  fast_ind = f_ind(Matrix, I, pElement->Col);
                  while (fast_ind_done < fast_ind) {
                    Matrix->Row_fast[I][++fast_ind_done] = last_el;
                  }
                  last_el = pElement;
                  ppElement = &pElement->NextInRow;  /* Skip element. */
              }
            }
            while (fast_ind_done < Matrix->Indsize)
              Matrix->Row_fast[I][++fast_ind_done] = last_el;
        }
    }
    spSetIndex(Matrix);
#ifdef CHECK_IND
    spCheckInd(Matrix, "spStripFills bottom");
#endif
    return;
}

int f_ind(MatrixPtr Matrix, int i, int j)
{
    int rval, diff, del, offset;

    del = (int) sqrt(fabs((double) IND_DENSITY*(i-j)));
    offset = (int) sqrt((double) IND_DENSITY*i);
    if (i > j)
      rval = offset-del;
    else
      rval = offset+del;

    return rval;
}



#if TRANSLATE AND DELETE
/*
 *  DELETE A ROW AND COLUMN FROM THE MATRIX
 *
 *  Deletes a row and a column from a matrix.
 *
 *  Sparse will abort if an attempt is made to delete a row or column that
 *  doesn't exist.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix in which the row and column are to be deleted.
 *  Row  <input>  (int)
 *      Row to be deleted.
 *  Col  <input>  (int)
 *      Column to be deleted.
 *
 *  >>> Local variables:
 *  ExtCol  (int)
 *      The external column that is being deleted.
 *  ExtRow  (int)
 *      The external row that is being deleted.
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.  Used when scanning rows and
 *      columns in order to eliminate elements from the last row or column.
 *  ppElement  (ElementPtr *)
 *      Pointer to the location of an ElementPtr.  This location will be
 *      filled with a NULL pointer if it is the new last element in its row
 *      or column.
 *  pElement  (ElementPtr)
 *      Pointer to an element in the last row or column of the matrix.
 *  Size  (int)
 *      The local version Matrix->Size, the size of the matrix.
 */

void
spDeleteRowAndCol( eMatrix, Row, Col )

char *eMatrix;
int  Row, Col;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement, *ppElement, pLastElement;
int  Size, ExtRow, ExtCol;
ElementPtr  spcFindElementInCol();

/* Begin `spDeleteRowAndCol'. */

    ASSERT( IS_SPARSE(Matrix) AND Row > 0 AND Col > 0 );
    ASSERT( Row <= Matrix->ExtSize AND Col <= Matrix->ExtSize );

    Size = Matrix->Size;
    ExtRow = Row;
    ExtCol = Col;
    if (NOT Matrix->RowsLinked) spcLinkRows( Matrix );

    Row = Matrix->ExtToIntRowMap[Row];
    Col = Matrix->ExtToIntColMap[Col];
    ASSERT( Row > 0 AND Col > 0 );

/* Move Row so that it is the last row in the matrix. */
    if (Row != Size) spcRowExchange( Matrix, Row, Size );

/* Move Col so that it is the last column in the matrix. */
    if (Col != Size) spcColExchange( Matrix, Col, Size );

/* Correct Diag pointers. */
    if (Row == Col)
        SWAP( ElementPtr, Matrix->Diag[Row], Matrix->Diag[Size] )
    else
    {   Matrix->Diag[Row] = spcFindElementInCol( Matrix, Matrix->FirstInCol+Row,
                                                 Row, Row, NO );
        Matrix->Diag[Col] = spcFindElementInCol( Matrix, Matrix->FirstInCol+Col,
                                                 Col, Col, NO );
    }

/*
 * Delete last row and column of the matrix.
 */
/* Break the column links to every element in the last row. */
    pLastElement = Matrix->FirstInRow[ Size ];
    while (pLastElement != NULL)
    {   ppElement = &(Matrix->FirstInCol[ pLastElement->Col ]);
        while ((pElement = *ppElement) != NULL)
        {   if (pElement == pLastElement)
                *ppElement = NULL;  /* Unlink last element in column. */
            else
                ppElement = &pElement->NextInCol;  /* Skip element. */
        }
        pLastElement = pLastElement->NextInRow;
    }

/* Break the row links to every element in the last column. */
    pLastElement = Matrix->FirstInCol[ Size ];
    while (pLastElement != NULL)
    {   ppElement = &(Matrix->FirstInRow[ pLastElement->Row ]);
        while ((pElement = *ppElement) != NULL)
        {   if (pElement == pLastElement)
                *ppElement = NULL;  /* Unlink last element in row. */
            else
                ppElement = &pElement->NextInRow;  /* Skip element. */
        }
        pLastElement = pLastElement->NextInCol;
    }

/* Clean up some details. */
    Matrix->Size = Size - 1;
    Matrix->Diag[Size] = NULL;
    Matrix->FirstInRow[Size] = NULL;
    Matrix->FirstInCol[Size] = NULL;
    Matrix->CurrentSize--;
    Matrix->ExtToIntRowMap[ExtRow] = -1;
    Matrix->ExtToIntColMap[ExtCol] = -1;
    Matrix->NeedsOrdering = YES;

    return;
}
#endif








#if PSEUDOCONDITION
/*
 *  CALCULATE PSEUDOCONDITION
 *
 *  Computes the magnitude of the ratio of the largest to the smallest
 *  pivots.  This quantity is an indicator of ill-conditioning in the
 *  matrix.  If this ratio is large, and if the matrix is scaled such
 *  that uncertainties in the RHS and the matrix entries are
 *  equilibrated, then the matrix is ill-conditioned.  However, a small
 *  ratio does not necessarily imply that the matrix is
 *  well-conditioned.  This routine must only be used after a matrix has
 *  been factored by spOrderAndFactor() or spFactor() and before it is
 *  cleared by spClear() or spInitialize().  The pseudocondition is
 *  faster to compute than the condition number calculated by
 *  spCondition(), but is not as informative.
 *
 *  >>> Returns:
 *  The magnitude of the ratio of the largest to smallest pivot used during
 *  previous factorization.  If the matrix was singular, zero is returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 */

RealNumber
spPseudoCondition( eMatrix )

char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register int I;
register ArrayOfElementPtrs Diag;
RealNumber MaxPivot, MinPivot, Mag;

/* Begin `spPseudoCondition'. */

    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE(Matrix) AND IS_FACTORED(Matrix) );
    if (Matrix->Error == spSINGULAR OR Matrix->Error == spZERO_DIAG)
        return 0.0;

    Diag = Matrix->Diag;
    MaxPivot = MinPivot = ELEMENT_MAG( Diag[1] );
    for (I = 2; I <= Matrix->Size; I++)
    {   Mag = ELEMENT_MAG( Diag[I] );
        if (Mag > MaxPivot)
            MaxPivot = Mag;
        else if (Mag < MinPivot)
            MinPivot = Mag;
    }
    ASSERT( MaxPivot > 0.0);
    return MaxPivot / MinPivot;
}
#endif








#if CONDITION
/*
 *  ESTIMATE CONDITION NUMBER
 *
 *  Computes an estimate of the condition number using a variation on
 *  the LINPACK condition number estimation algorithm.  This quantity is
 *  an indicator of ill-conditioning in the matrix.  To avoid problems
 *  with overflow, the reciprocal of the condition number is returned.
 *  If this number is small, and if the matrix is scaled such that
 *  uncertainties in the RHS and the matrix entries are equilibrated,
 *  then the matrix is ill-conditioned.  If the this number is near
 *  one, the matrix is well conditioned.  This routine must only be
 *  used after a matrix has been factored by spOrderAndFactor() or
 *  spFactor() and before it is cleared by spClear() or spInitialize().
 *
 *  Unlike the LINPACK condition number estimator, this routines
 *  returns the L infinity condition number.  This is an artifact of
 *  Sparse placing ones on the diagonal of the upper triangular matrix
 *  rather than the lower.  This difference should be of no importance.
 *
 *  References:
 *  A.K. Cline, C.B. Moler, G.W. Stewart, J.H. Wilkinson.  An estimate
 *  for the condition number of a matrix.  SIAM Journal on Numerical
 *  Analysis.  Vol. 16, No. 2, pages 368-375, April 1979.
 *  
 *  J.J. Dongarra, C.B. Moler, J.R. Bunch, G.W. Stewart.  LINPACK
 *  User's Guide.  SIAM, 1979.
 *  
 *  Roger G. Grimes, John G. Lewis.  Condition number estimation for
 *  sparse matrices.  SIAM Journal on Scientific and Statistical
 *  Computing.  Vol. 2, No. 4, pages 384-388, December 1981.
 *  
 *  Dianne Prost O'Leary.  Estimating matrix condition numbers.  SIAM
 *  Journal on Scientific and Statistical Computing.  Vol. 1, No. 2,
 *  pages 205-209, June 1980.
 *
 *  >>> Returns:
 *  The reciprocal of the condition number.  If the matrix was singular,
 *  zero is returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 *  NormOfMatrix  <input>  (RealNumber)
 *      The L-infinity norm of the unfactored matrix as computed by
 *      spNorm().
 *  pError  <output>  (int *)
 *      Used to return error code.
 *
 *  >>> Possible errors:
 *  spSINGULAR
 *  spNO_MEMORY
 */

RealNumber
spCondition( eMatrix, NormOfMatrix, pError )

char *eMatrix;
RealNumber NormOfMatrix;
int *pError;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
register RealVector T, Tm;
register int I, K, Row;
ElementPtr pPivot;
int Size;
RealNumber E, Em, Wp, Wm, ASp, ASm, ASw, ASy, ASv, ASz, MaxY, ScaleFactor;
RealNumber Linpack, OLeary, InvNormOfInverse, ComplexCondition();
#define SLACK   1e4

/* Begin `spCondition'. */

    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE(Matrix) AND IS_FACTORED(Matrix) );
    *pError = Matrix->Error;
    if (Matrix->Error >= spFATAL) return 0.0;
    if (NormOfMatrix == 0.0)
    {   *pError = spSINGULAR;
        return 0.0;
    }

#if spCOMPLEX
    if (Matrix->Complex)
        return ComplexCondition( Matrix, NormOfMatrix, pError );
#endif

#if REAL
    Size = Matrix->Size;
    T = Matrix->Intermediate;
#if spCOMPLEX
    Tm = Matrix->Intermediate + Size;
#else
    Tm = ALLOC( RealNumber, Size+1 );
    if (Tm == NULL)
    {   *pError = spNO_MEMORY;
        return 0.0;
    }
#endif
    for (I = Size; I > 0; I--) T[I] = 0.0;

/*
 * Part 1.  Ay = e.
 * Solve Ay = LUy = e where e consists of +1 and -1 terms with the sign
 * chosen to maximize the size of w in Lw = e.  Since the terms in w can
 * get very large, scaling is used to avoid overflow.
 */

/* Forward elimination. Solves Lw = e while choosing e. */
    E = 1.0;
    for (I = 1; I <= Size; I++)
    {   pPivot = Matrix->Diag[I];
        if (T[I] < 0.0) Em = -E; else Em = E;
        Wm = (Em + T[I]) * pPivot->Real;
        if (ABS(Wm) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(Wm) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            E *= ScaleFactor;
            Em *= ScaleFactor;
            Wm = (Em + T[I]) * pPivot->Real;
        }
        Wp = (T[I] - Em) * pPivot->Real;
        ASp = ABS(T[I] - Em);
        ASm = ABS(Em + T[I]);

/* Update T for both values of W, minus value is placed in Tm. */
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   Row = pElement->Row;
            Tm[Row] = T[Row] - (Wm * pElement->Real);
            T[Row] -= (Wp * pElement->Real);
            ASp += ABS(T[Row]);
            ASm += ABS(Tm[Row]);
            pElement = pElement->NextInCol;
        }

/* If minus value causes more growth, overwrite T with its values. */
        if (ASm > ASp)
        {   T[I] = Wm;
            pElement = pPivot->NextInCol;
            while (pElement != NULL)
            {   T[pElement->Row] = Tm[pElement->Row];
                pElement = pElement->NextInCol;
            }
        }
        else T[I] = Wp;
    }

/* Compute 1-norm of T, which now contains w, and scale ||T|| to 1/SLACK. */
    for (ASw = 0.0, I = Size; I > 0; I--) ASw += ABS(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASw);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) T[I] *= ScaleFactor;
        E *= ScaleFactor;
    }

/* Backward Substitution. Solves Uy = w.*/
    for (I = Size; I >= 1; I--)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   T[I] -= pElement->Real * T[pElement->Col];
            pElement = pElement->NextInRow;
        }
        if (ABS(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(T[I]) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            E *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains y, and scale ||T|| to 1/SLACK. */
    for (ASy = 0.0, I = Size; I > 0; I--) ASy += ABS(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASy);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) T[I] *= ScaleFactor;
        ASy = 1.0 / SLACK;
        E *= ScaleFactor;
    }

/* Compute infinity-norm of T for O'Leary's estimate. */
    for (MaxY = 0.0, I = Size; I > 0; I--)
        if (MaxY < ABS(T[I])) MaxY = ABS(T[I]);

/*
 * Part 2.  A* z = y where the * represents the transpose.
 * Recall that A = LU implies A* = U* L*.
 */

/* Forward elimination, U* v = y. */
    for (I = 1; I <= Size; I++)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   T[pElement->Col] -= T[I] * pElement->Real;
            pElement = pElement->NextInRow;
        }
        if (ABS(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(T[I]) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains v, and scale ||T|| to 1/SLACK. */
    for (ASv = 0.0, I = Size; I > 0; I--) ASv += ABS(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASv);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) T[I] *= ScaleFactor;
        ASy *= ScaleFactor;
    }

/* Backward Substitution, L* z = v. */
    for (I = Size; I >= 1; I--)
    {   pPivot = Matrix->Diag[I];
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   T[I] -= pElement->Real * T[pElement->Row];
            pElement = pElement->NextInCol;
        }
        T[I] *= pPivot->Real;
        if (ABS(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), ABS(T[I]) );
            for (K = Size; K > 0; K--) T[K] *= ScaleFactor;
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains z. */
    for (ASz = 0.0, I = Size; I > 0; I--) ASz += ABS(T[I]);

#if NOT spCOMPLEX
    FREE( Tm );
#endif

    Linpack = ASy / ASz;
    OLeary = E / MaxY;
    InvNormOfInverse = MIN( Linpack, OLeary );
    return InvNormOfInverse / NormOfMatrix;
#endif /* REAL */
}





#if spCOMPLEX
/*
 *  ESTIMATE CONDITION NUMBER
 *
 *  Complex version of spCondition().
 *
 *  >>> Returns:
 *  The reciprocal of the condition number.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *  NormOfMatrix  <input>  (RealNumber)
 *      The L-infinity norm of the unfactored matrix as computed by
 *      spNorm().
 *  pError  <output>  (int *)
 *      Used to return error code.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

static RealNumber
ComplexCondition( Matrix, NormOfMatrix, pError )

MatrixPtr Matrix;
RealNumber NormOfMatrix;
int *pError;
{
register ElementPtr pElement;
register ComplexVector T, Tm;
register int I, K, Row;
ElementPtr pPivot;
int Size;
RealNumber E, Em, ASp, ASm, ASw, ASy, ASv, ASz, MaxY, ScaleFactor;
RealNumber Linpack, OLeary, InvNormOfInverse;
ComplexNumber Wp, Wm;

/* Begin `ComplexCondition'. */

    Size = Matrix->Size;
    T = (ComplexVector)Matrix->Intermediate;
    Tm = ALLOC( ComplexNumber, Size+1 );
    if (Tm == NULL)
    {   *pError = spNO_MEMORY;
        return 0.0;
    }
    for (I = Size; I > 0; I--) T[I].Real = T[I].Imag = 0.0;

/*
 * Part 1.  Ay = e.
 * Solve Ay = LUy = e where e consists of +1 and -1 terms with the sign
 * chosen to maximize the size of w in Lw = e.  Since the terms in w can
 * get very large, scaling is used to avoid overflow.
 */

/* Forward elimination. Solves Lw = e while choosing e. */
    E = 1.0;
    for (I = 1; I <= Size; I++)
    {   pPivot = Matrix->Diag[I];
        if (T[I].Real < 0.0) Em = -E; else Em = E;
        Wm = T[I];
        Wm.Real += Em;
        ASm = CMPLX_1_NORM( Wm );
        CMPLX_MULT_ASSIGN( Wm, *pPivot );
        if (CMPLX_1_NORM(Wm) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(Wm) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            E *= ScaleFactor;
            Em *= ScaleFactor;
            ASm *= ScaleFactor;
            SCLR_MULT_ASSIGN( Wm, ScaleFactor );
        }
        Wp = T[I];
        Wp.Real -= Em;
        ASp = CMPLX_1_NORM( Wp );
        CMPLX_MULT_ASSIGN( Wp, *pPivot );

/* Update T for both values of W, minus value is placed in Tm. */
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   Row = pElement->Row;
            /* Cmplx expr: Tm[Row] = T[Row] - (Wp * *pElement). */
            CMPLX_MULT_SUBT( Tm[Row], Wm, *pElement, T[Row] );
            /* Cmplx expr: T[Row] -= Wp * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[Row], Wm, *pElement );
            ASp += CMPLX_1_NORM(T[Row]);
            ASm += CMPLX_1_NORM(Tm[Row]);
            pElement = pElement->NextInCol;
        }

/* If minus value causes more growth, overwrite T with its values. */
        if (ASm > ASp)
        {   T[I] = Wm;
            pElement = pPivot->NextInCol;
            while (pElement != NULL)
            {   T[pElement->Row] = Tm[pElement->Row];
                pElement = pElement->NextInCol;
            }
        }
        else T[I] = Wp;
    }

/* Compute 1-norm of T, which now contains w, and scale ||T|| to 1/SLACK. */
    for (ASw = 0.0, I = Size; I > 0; I--) ASw += CMPLX_1_NORM(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASw);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) SCLR_MULT_ASSIGN( T[I], ScaleFactor );
        E *= ScaleFactor;
    }

/* Backward Substitution. Solves Uy = w.*/
    for (I = Size; I >= 1; I--)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   /* Cmplx expr: T[I] -= T[pElement->Col] * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[I], T[pElement->Col], *pElement );
            pElement = pElement->NextInRow;
        }
        if (CMPLX_1_NORM(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(T[I]) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            E *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains y, and scale ||T|| to 1/SLACK. */
    for (ASy = 0.0, I = Size; I > 0; I--) ASy += CMPLX_1_NORM(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASy);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) SCLR_MULT_ASSIGN( T[I], ScaleFactor );
        ASy = 1.0 / SLACK;
        E *= ScaleFactor;
    }

/* Compute infinity-norm of T for O'Leary's estimate. */
    for (MaxY = 0.0, I = Size; I > 0; I--)
        if (MaxY < CMPLX_1_NORM(T[I])) MaxY = CMPLX_1_NORM(T[I]);

/*
 * Part 2.  A* z = y where the * represents the transpose.
 * Recall that A = LU implies A* = U* L*.
 */

/* Forward elimination, U* v = y. */
    for (I = 1; I <= Size; I++)
    {   pElement = Matrix->Diag[I]->NextInRow;
        while (pElement != NULL)
        {   /* Cmplx expr: T[pElement->Col] -= T[I] * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[pElement->Col], T[I], *pElement );
            pElement = pElement->NextInRow;
        }
        if (CMPLX_1_NORM(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(T[I]) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains v, and scale ||T|| to 1/SLACK. */
    for (ASv = 0.0, I = Size; I > 0; I--) ASv += CMPLX_1_NORM(T[I]);
    ScaleFactor = 1.0 / (SLACK * ASv);
    if (ScaleFactor < 0.5)
    {   for (I = Size; I > 0; I--) SCLR_MULT_ASSIGN( T[I], ScaleFactor );
        ASy *= ScaleFactor;
    }

/* Backward Substitution, L* z = v. */
    for (I = Size; I >= 1; I--)
    {   pPivot = Matrix->Diag[I];
        pElement = pPivot->NextInCol;
        while (pElement != NULL)
        {   /* Cmplx expr: T[I] -= T[pElement->Row] * *pElement. */
            CMPLX_MULT_SUBT_ASSIGN( T[I], T[pElement->Row], *pElement );
            pElement = pElement->NextInCol;
        }
        CMPLX_MULT_ASSIGN( T[I], *pPivot );
        if (CMPLX_1_NORM(T[I]) > SLACK)
        {   ScaleFactor = 1.0 / MAX( SQR( SLACK ), CMPLX_1_NORM(T[I]) );
            for (K = Size; K > 0; K--) SCLR_MULT_ASSIGN( T[K], ScaleFactor );
            ASy *= ScaleFactor;
        }
    }

/* Compute 1-norm of T, which now contains z. */
    for (ASz = 0.0, I = Size; I > 0; I--) ASz += CMPLX_1_NORM(T[I]);

    FREE( Tm );

    Linpack = ASy / ASz;
    OLeary = E / MaxY;
    InvNormOfInverse = MIN( Linpack, OLeary );
    return InvNormOfInverse / NormOfMatrix;
}
#endif /* spCOMPLEX */





/*
 *  L-INFINITY MATRIX NORM 
 *
 *  Computes the L-infinity norm of an unfactored matrix.  It is a fatal
 *  error to pass this routine a factored matrix.
 *
 *  One difficulty is that the rows may not be linked.
 *
 *  >>> Returns:
 *  The largest absolute row sum of matrix.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 */

RealNumber
spNorm( eMatrix )

char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
register int I;
RealNumber Max = 0.0, AbsRowSum;

/* Begin `spNorm'. */
    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE(Matrix) AND NOT IS_FACTORED(Matrix) );
    if (NOT Matrix->RowsLinked) spcLinkRows( Matrix );

/* Compute row sums. */
#if REAL
    if (NOT Matrix->Complex)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInRow[I];
            AbsRowSum = 0.0;
            while (pElement != NULL)
            {   AbsRowSum += ABS( pElement->Real );
                pElement = pElement->NextInRow;
            }
            if (Max < AbsRowSum) Max = AbsRowSum;
        }
    }
#endif
#if spCOMPLEX
    if (Matrix->Complex)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInRow[I];
            AbsRowSum = 0.0;
            while (pElement != NULL)
            {   AbsRowSum += CMPLX_1_NORM( *pElement );
                pElement = pElement->NextInRow;
            }
            if (Max < AbsRowSum) Max = AbsRowSum;
        }
    }
#endif
    return Max;
}
#endif /* CONDITION */






#if STABILITY
/*
 *  STABILITY OF FACTORIZATION
 *
 *  The following routines are used to gauge the stability of a
 *  factorization.  If the factorization is determined to be too unstable,
 *  then the matrix should be reordered.  The routines compute quantities
 *  that are needed in the computation of a bound on the error attributed
 *  to any one element in the matrix during the factorization.  In other
 *  words, there is a matrix E = [e_ij] of error terms such that A+E = LU.
 *  This routine finds a bound on |e_ij|.  Erisman & Reid [1] showed that
 *  |e_ij| < 3.01 u rho m_ij, where u is the machine rounding unit,
 *  rho = max a_ij where the max is taken over every row i, column j, and
 *  step k, and m_ij is the number of multiplications required in the
 *  computation of l_ij if i > j or u_ij otherwise.  Barlow [2] showed that
 *  rho < max_i || l_i ||_p max_j || u_j ||_q where 1/p + 1/q = 1.
 *
 *  The first routine finds the magnitude on the largest element in the
 *  matrix.  If the matrix has not yet been factored, the largest
 *  element is found by direct search.  If the matrix is factored, a
 *  bound on the largest element in any of the reduced submatrices is
 *  computed using Barlow with p = oo and q = 1.  The ratio of these
 *  two numbers is the growth, which can be used to determine if the
 *  pivoting order is adequate.  A large growth implies that
 *  considerable error has been made in the factorization and that it
 *  is probably a good idea to reorder the matrix.  If a large growth
 *  in encountered after using spFactor(), reconstruct the matrix and
 *  refactor using spOrderAndFactor().  If a large growth is
 *  encountered after using spOrderAndFactor(), refactor using
 *  spOrderAndFactor() with the pivot threshold increased, say to 0.1.
 *
 *  Using only the size of the matrix as an upper bound on m_ij and
 *  Barlow's bound, the user can estimate the size of the matrix error
 *  terms e_ij using the bound of Erisman and Reid.  The second routine
 *  computes a tighter bound (with more work) based on work by Gear
 *  [3], |e_ij| < 1.01 u rho (t c^3 + (1 + t)c^2) where t is the
 *  threshold and c is the maximum number of off-diagonal elements in
 *  any row of L.  The expensive part of computing this bound is
 *  determining the maximum number of off-diagonals in L, which changes
 *  only when the order of the matrix changes.  This number is computed
 *  and saved, and only recomputed if the matrix is reordered.
 *
 *  [1] A. M. Erisman, J. K. Reid.  Monitoring the stability of the
 *      triangular factorization of a sparse matrix.  Numerische
 *      Mathematik.  Vol. 22, No. 3, 1974, pp 183-186.
 *
 *  [2] J. L. Barlow.  A note on monitoring the stability of triangular
 *      decomposition of sparse matrices.  "SIAM Journal of Scientific
 *      and Statistical Computing."  Vol. 7, No. 1, January 1986, pp 166-168.
 *
 *  [3] I. S. Duff, A. M. Erisman, J. K. Reid.  "Direct Methods for Sparse
 *      Matrices."  Oxford 1986. pp 99.
 */

/*
 *  LARGEST ELEMENT IN MATRIX
 *
 *  >>> Returns:
 *  If matrix is not factored, returns the magnitude of the largest element in
 *  the matrix.  If the matrix is factored, a bound on the magnitude of the
 *  largest element in any of the reduced submatrices is returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 */

RealNumber
spLargestElement( eMatrix )

char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register int I;
RealNumber Mag, AbsColSum, Max = 0.0, MaxRow = 0.0, MaxCol = 0.0;
RealNumber Pivot;
ComplexNumber cPivot;
register ElementPtr pElement, pDiag;

/* Begin `spLargestElement'. */
    ASSERT( IS_SPARSE(Matrix) );
    spExpandFormat(Matrix);
#if REAL
    if (Matrix->Factored AND NOT Matrix->Complex)
    {   if (Matrix->Error == spSINGULAR) return 0.0;

/* Find the bound on the size of the largest element over all factorization. */
        for (I = 1; I <= Matrix->Size; I++)
        {   pDiag = Matrix->Diag[I];

/* Lower triangular matrix. */
            Pivot = 1.0 / pDiag->Real;
            Mag = ABS( Pivot );
            if (Mag > MaxRow) MaxRow = Mag;
            pElement = Matrix->FirstInRow[I];
            while (pElement != pDiag)
            {   Mag = ABS( pElement->Real );
                if (Mag > MaxRow) MaxRow = Mag;
                pElement = pElement->NextInRow;
            }

/* Upper triangular matrix. */
            pElement = Matrix->FirstInCol[I];
            AbsColSum = 1.0;  /* Diagonal of U is unity. */
            while (pElement != pDiag)
            {   AbsColSum += ABS( pElement->Real );
                pElement = pElement->NextInCol;
            }
            if (AbsColSum > MaxCol) MaxCol = AbsColSum;
        }
    }
    else if (NOT Matrix->Complex)
    {   for (I = 1; I <= Matrix->Size; I++)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   Mag = ABS( pElement->Real );
                if (Mag > Max) Max = Mag;
                pElement = pElement->NextInCol;
            }
        }
        return Max;
    }
#endif
#if spCOMPLEX
    if (Matrix->Factored AND Matrix->Complex)
    {   if (Matrix->Error == spSINGULAR) return 0.0;

/* Find the bound on the size of the largest element over all factorization. */
        for (I = 1; I <= Matrix->Size; I++)
        {   pDiag = Matrix->Diag[I];

/* Lower triangular matrix. */
            CMPLX_RECIPROCAL( cPivot, *pDiag );
            Mag = CMPLX_1_NORM( cPivot );
            if (Mag > MaxRow) MaxRow = Mag;
            pElement = Matrix->FirstInRow[I];
            while (pElement != pDiag)
            {   Mag = CMPLX_1_NORM( *pElement );
                if (Mag > MaxRow) MaxRow = Mag;
                pElement = pElement->NextInRow;
            }

/* Upper triangular matrix. */
            pElement = Matrix->FirstInCol[I];
            AbsColSum = 1.0;  /* Diagonal of U is unity. */
            while (pElement != pDiag)
            {   AbsColSum += CMPLX_1_NORM( *pElement );
                pElement = pElement->NextInCol;
            }
            if (AbsColSum > MaxCol) MaxCol = AbsColSum;
        }
    }
    else if (Matrix->Complex)
    {   for (I = 1; I <= Matrix->Size; I++)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   Mag = CMPLX_1_NORM( *pElement );
                if (Mag > Max) Max = Mag;
                pElement = pElement->NextInCol;
            }
        }
        return Max;
    }
#endif
    return MaxRow * MaxCol;
}




/*
 *  MATRIX ROUNDOFF ERROR
 *
 *  >>> Returns:
 *  Returns a bound on the magnitude of the largest element in E = A - LU.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to the matrix.
 *  Rho  <input>  (RealNumber)
 *      The bound on the magnitude of the largest element in any of the
 *      reduced submatrices.  This is the number computed by the function
 *      spLargestElement() when given a factored matrix.  If this number is
 *      negative, the bound will be computed automatically.
 */

RealNumber
spRoundoff( eMatrix, Rho )

char *eMatrix;
RealNumber Rho;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
register int Count, I, MaxCount = 0;
RealNumber Reid, Gear;

/* Begin `spRoundoff'. */
    spExpandFormat(Matrix);
    ASSERT( IS_SPARSE(Matrix) AND IS_FACTORED(Matrix) );

/* Compute Barlow's bound if it is not given. */
    if (Rho < 0.0) Rho = spLargestElement( eMatrix );

/* Find the maximum number of off-diagonals in L if not previously computed. */
    if (Matrix->MaxRowCountInLowerTri < 0)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInRow[I];
            Count = 0;
            while (pElement->Col < I)
            {   Count++;
                pElement = pElement->NextInRow;
            }
            if (Count > MaxCount) MaxCount = Count;
        }
        Matrix->MaxRowCountInLowerTri = MaxCount;
    }
    else MaxCount = Matrix->MaxRowCountInLowerTri;

/* Compute error bound. */
    Gear = 1.01*((MaxCount + 1) * Matrix->RelThreshold + 1.0) * SQR(MaxCount);
    Reid = 3.01 * Matrix->Size;

    if (Gear < Reid)
        return (MACHINE_RESOLUTION * Rho * Gear);
    else
        return (MACHINE_RESOLUTION * Rho * Reid);
}
#endif








#if DOCUMENTATION
/*
 *  SPARSE ERROR MESSAGE
 *
 *  This routine prints a short message to a stream describing the error
 *  error state of sparse.  No message is produced if there is no error.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *	Matrix for which the error message is to be printed.
 *  Stream  <input>  (FILE *)
 *	Stream to which the error message is to be printed.
 *  Originator  <input>  (char *)
 *	Name of originator of error message.  If NULL, `sparse' is used.
 *	If zero-length string, no originator is printed.
 */

void
spErrorMessage( eMatrix, Stream, Originator )

char *eMatrix, *Originator;
FILE *Stream;
{
int Row, Col, Error;

/* Begin `spErrorMessage'. */
    if (eMatrix == NULL)
	Error = spNO_MEMORY;
    else
    {   ASSERT(((MatrixPtr)eMatrix)->ID == SPARSE_ID);
	Error = ((MatrixPtr)eMatrix)->Error;
    }

    if (Error == spOKAY) return;
    if (Originator == NULL) Originator = "sparse";
    if (Originator[0] != '\0') fprintf( Stream, "%s: ", Originator);
    if (Error >= spFATAL)
	fprintf( Stream, "fatal error, ");
    else
	fprintf( Stream, "warning, ");
/*
 * Print particular error message.
 * Do not use switch statement because error codes may not be unique.
 */
    if (Error == spPANIC)
	fprintf( Stream, "Sparse called improperly.\n");
    else if (Error == spNO_MEMORY)
	fprintf( Stream, "insufficient memory available.\n");
    else if (Error == spSINGULAR)
    {   spWhereSingular( eMatrix, &Row, &Col );
	fprintf( Stream, "singular matrix detected at row %d and column %d.\n",
		 Row, Col);
    }
    else if (Error == spZERO_DIAG)
    {   spWhereSingular( eMatrix, &Row, &Col );
	fprintf( Stream, "zero diagonal detected at row %d and column %d.\n",
		 Row, Col);
    }
    else if (Error == spSMALL_PIVOT)
    {   fprintf( Stream,
	    "unable to find a pivot that is larger than absolute threshold.\n");
    }
    else ABORT();
    return;
}
#endif /* DOCUMENTATION */

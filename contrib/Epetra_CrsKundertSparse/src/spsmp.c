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
 *  Spice3 COMPATIBILITY MODULE
 *
 *  Author:                     Advising professor:
 *     Kenneth S. Kundert           Alberto Sangiovanni-Vincentelli
 *     UC Berkeley
 *
 *  This module contains routines that make Sparse1.3 a direct
 *  replacement for the SMP sparse matrix package in Spice3c1 or Spice3d1.
 *  Sparse1.3 is in general a faster and more robust package than SMP.
 *  These advantages become significant on large circuits.
 *
 *  >>> User accessible functions contained in this file:
 *  SMPaddElt
 *  SMPmakeElt
 *  SMPcClear
 *  SMPclear
 *  SMPcLUfac
 *  SMPluFac
 *  SMPcReorder
 *  SMPreorder
 *  SMPcaSolve
 *  SMPcSolve
 *  SMPsolve
 *  SMPmatSize
 *  SMPnewMatrix
 *  SMPdestroy
 *  SMPpreOrder
 *  SMPprint
 *  SMPgetError
 *  SMPcProdDiag
 *  LoadGmin
 *  SMPfindElt
 */

/*
 *  To replace SMP with Sparse, rename the file spSpice3.h to
 *  spMatrix.h and place Sparse in a subdirectory of SPICE called
 *  `sparse'.  Then on UNIX compile Sparse by executing `make spice'.
 *  If not on UNIX, after compiling Sparse and creating the sparse.a
 *  archive, compile this file (spSMP.c) and spSMP.o to the archive,
 *  then copy sparse.a into the SPICE main directory and rename it
 *  SMP.a.  Finally link SPICE.
 *
 *  To be compatible with SPICE, the following Sparse compiler options
 *  (in spConfig.h) should be set as shown below:
 *
 *      REAL                            YES
 *      EXPANDABLE                      YES
 *      TRANSLATE                       NO
 *      INITIALIZE                      NO or YES, YES for use with test prog.
 *      DIAGONAL_PIVOTING               YES
 *      ARRAY_OFFSET                    YES
 *      MODIFIED_MARKOWITZ              NO
 *      DELETE                          NO
 *      MODIFIED_NODAL                  YES
 *      QUAD_ELEMENT                    NO
 *      TRANSPOSE                       YES
 *      SCALING                         NO
 *      DOCUMENTATION                   YES
 *      MULTIPLICATION                  NO
 *      DETERMINANT                     YES
 *      STABILITY                       NO
 *      CONDITION                       NO
 *      PSEUDOCONDITION                 NO
 *      FORTRAN                         NO
 *      DEBUG                           YES
 *      spCOMPLEX                       1
 *      spSEPARATED_COMPLEX_VECTORS     1
 *      spCOMPATIBILITY                 0
 *
 *      spREAL  double
 */

/*
 *  Revision and copyright information.
 *
 *  Copyright (c) 1985,86,87,88,89,90
 *  by Kenneth S. Kundert and the University of California.
 *
 *  Permission to use, copy, modify, and distribute this software and its
 *  documentation for any purpose and without fee is hereby granted, provided
 *  that the above copyright notice appear in all copies and supporting
 *  documentation and that the authors and the University of California
 *  are properly credited.  The authors and the University of California
 *  make no representations as to the suitability of this software for
 *  any purpose.  It is provided `as is', without express or implied warranty.
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
 *  spMatrix.h
 *     Sparse macros and declarations.
 *  SMPdefs.h
 *     Spice3's matrix macro definitions.
 */

#include "spice.h"
#include <stdio.h>
#include "spmatrix.h"
#include "smpdefs.h"
#include "spdefs.h"

/* #define NO   0 */
/* #define YES  1 */

#ifdef __STDC__
static void LoadGmin( char * /*eMatrix*/, double /*Gmin*/ );
#else
static void LoadGmin( );
#endif

/*
 * SMPaddElt()
 */
int
SMPaddElt( Matrix, Row, Col, Value )
SMPmatrix *Matrix;
int Row, Col;
double Value;
{
    *spGetElement( (char *)Matrix, Row, Col ) = Value;
    return spError( (char *)Matrix );
}

/*
 * SMPmakeElt 
 */
double* SMPmakeElt( Matrix, Row, Col )
SMPmatrix *Matrix;
int Row, Col;
{
  MatrixPtr csr_matrix = (MatrixPtr)Matrix;
  int ptr;

  if( csr_matrix->Cscflag == 0 ){
      return spGetElement( (char *)Matrix, Row, Col );
  }
  else{
    if( (Col>0) && (Row>0) ){
      for (ptr = csr_matrix->Column_pointers[Col-1]; 
        ptr < csr_matrix->Column_pointers[Col]; ptr++){
        if( csr_matrix->Row_indices[ptr] == Row-1 ) 
          return &csr_matrix->Values[ptr];
      }
    } 
    else{

      /* Ignore the electrical ground */

      return &csr_matrix->TrashCan.Real;
    }
  }
  printf("SMPmakeElt:Error matrix entry not found\n");
  return NULL;
}
/*
double *
SMPmakeElt( Matrix, Row, Col )
SMPmatrix *Matrix;
int Row, Col;
{
    return spGetElement( (char *)Matrix, Row, Col );
}
*/
/*
 * SMPcClear()
 */
void
SMPcClear( Matrix )
SMPmatrix *Matrix;
{
    spClear( (char *)Matrix );
}

/*
 * SMPclear()
 */
void
SMPclear( Matrix )
SMPmatrix *Matrix;
{
    spClear( (char *)Matrix );
}

/*
 * SMPcLUfac()
 */
/*ARGSUSED*/
int
SMPcLUfac( Matrix, PivTol )
SMPmatrix *Matrix;
double PivTol;
{
    spSetComplex( (char *)Matrix );
    return spFactor( (char *)Matrix );
}

/*
 * SMPluFac()
 */
/*ARGSUSED*/
int
SMPluFac( Matrix, PivTol, Gmin )
SMPmatrix *Matrix;
double PivTol, Gmin;
{
    spSetReal( (char *)Matrix );
/*
    LoadGmin( (char *)Matrix, Gmin );
*/
    return spFactor( (char *)Matrix );
}

/*
 * SMPcReorder()
 */
int
SMPcReorder( Matrix, PivTol, PivRel, NumSwaps )
SMPmatrix *Matrix;
double PivTol, PivRel;
int *NumSwaps;
{
    *NumSwaps = 1;
    spSetComplex( (char *)Matrix );
    return spOrderAndFactor( (char *)Matrix, (spREAL*)NULL,
                             (spREAL)PivRel, (spREAL)PivTol, YES );
}

/*
 * SMPreorder()
 */
int
SMPreorder( Matrix, PivTol, PivRel, Gmin )
SMPmatrix *Matrix;
double PivTol, PivRel, Gmin;
{
    spSetReal( (char *)Matrix );
/*
    LoadGmin( (char *)Matrix, Gmin );
*/
    return spOrderAndFactor( (char *)Matrix, (spREAL*)NULL,
                             (spREAL)PivRel, (spREAL)PivTol, YES );
}

/*
 * SMPcaSolve()
 */
void
SMPcaSolve( Matrix, RHS, iRHS, Spare, iSpare)
SMPmatrix *Matrix;
double RHS[], iRHS[], Spare[], iSpare[];
{
    spSolveTransposed( (char *)Matrix, RHS, RHS, iRHS, iRHS );
}

/*
 * SMPcSolve()
 */
int
SMPcSolve( Matrix, RHS, iRHS, Spare, iSpare)
SMPmatrix *Matrix;
double RHS[], iRHS[], Spare[], iSpare[];
{
    return (spSolve( (char *)Matrix, RHS, RHS, iRHS, iRHS ));
}

/*
 * SMPsolve()
 */
int
SMPsolve( Matrix, RHS, Spare )
SMPmatrix *Matrix;
double RHS[], Spare[];
{
    return (spSolve( (char *)Matrix, RHS, RHS, (spREAL*)NULL, (spREAL*)NULL ));
}

/*
 * SMPmatSize()
 */
int
SMPmatSize( Matrix )
SMPmatrix *Matrix;
{
    return spGetSize( (char *)Matrix, 1 );
}

/*
 * SMPnewMatrix()
 */
int
SMPnewMatrix( pMatrix )
SMPmatrix **pMatrix;
{
int Error;
    *pMatrix = (SMPmatrix *)spCreate( 0, 1, &Error );
    return Error;
}

/*
 * SMPdestroy()
 */
void
SMPdestroy( Matrix )
SMPmatrix *Matrix;
{
    spDestroy( (char *)Matrix );
}

/*
 * SMPpreOrder()
 */
int
SMPpreOrder( Matrix )
SMPmatrix *Matrix;
{
    spMNA_Preorder( (char *)Matrix );
    return spError( (char *)Matrix );
}


/*
 * SMPprint()
 */
/*ARGSUSED*/
void
SMPprint( Matrix, filename)
SMPmatrix *Matrix;
char *filename;
{
    char  tmp[16];
    sprintf(tmp,"%s.1",filename);
    spFileMatrix( (char *) Matrix, tmp      , tmp     ,0,1,1);
 
    sprintf(tmp,"%s.2",filename);
    spFileMatrix( (char *) Matrix, tmp      , tmp     ,1,1,1);
 
    sprintf(tmp,"%s.stat",filename);
    spFileStats ( (char *) Matrix, tmp      ,tmp      );
}

/*
 * SMPgetError()
 */
void
SMPgetError( Matrix, Col, Row)
SMPmatrix *Matrix;
int *Row, *Col;
{
    spWhereSingular( (char *)Matrix, Row, Col );
}

/*
 * SMPcProdDiag()
 *    note: obsolete for Spice3d2 and later
 */
int
SMPcProdDiag( Matrix, pMantissa, pExponent)
SMPmatrix *Matrix;
SPcomplex *pMantissa;
int *pExponent;
{
    spDeterminant( (char *)Matrix, pExponent, &(pMantissa->real),
                                              &(pMantissa->imag) );
    return spError( (char *)Matrix );
}

/*
 * SMPcDProd()
 */
int
SMPcDProd( Matrix, pMantissa, pExponent)
SMPmatrix *Matrix;
SPcomplex *pMantissa;
int *pExponent;
{
    double	re, im, x, y, z;
    int		p;

    spDeterminant( (char *)Matrix, &p, &re, &im);

#ifndef M_LN2
#define M_LN2   0.69314718055994530942
#endif
#ifndef M_LN10
#define M_LN10  2.30258509299404568402
#endif

#ifdef debug_print
    printf("Determinant 10: (%20g,%20g)^%d\n", re, im, p);
#endif

    /* Convert base 10 numbers to base 2 numbers, for comparison */
    y = p * M_LN10 / M_LN2;
    x = (int) y;
    y -= x;

    /* ASSERT
     *	x = integral part of exponent, y = fraction part of exponent
     */

    /* Fold in the fractional part */
#ifdef debug_print
    printf(" ** base10 -> base2 int =  %g, frac = %20g\n", x, y);
#endif
    z = pow(2.0, y);
    re *= z;
    im *= z;
#ifdef debug_print
    printf(" ** multiplier = %20g\n", z);
#endif

    /* Re-normalize (re or im may be > 2.0 or both < 1.0 */
    if (re != 0.0) {
	y = logb(re);
	if (im != 0.0)
	    z = logb(im);
	else
	    z = 0;
    } else if (im != 0.0) {
	z = logb(im);
	y = 0;
    } else {
	/* Singular */
	/*printf("10 -> singular\n");*/
	y = 0;
	z = 0;
    }

#ifdef debug_print
    printf(" ** renormalize changes = %g,%g\n", y, z);
#endif
    if (y < z)
	y = z;

    *pExponent = x + y;
    x = scalb(re, (int) -y);
    z = scalb(im, (int) -y);
#ifdef debug_print
    printf(" ** values are: re %g, im %g, y %g, re' %g, im' %g\n",
	    re, im, y, x, z);
#endif
    pMantissa->real = scalb(re, (int) -y);
    pMantissa->imag = scalb(im, (int) -y);

#ifdef debug_print
    printf("Determinant 10->2: (%20g,%20g)^%d\n", pMantissa->real,
	pMantissa->imag, *pExponent);
#endif
    return spError( (char *)Matrix );
}



/*
 *  The following routines need internal knowledge of the Sparse data
 *  structures.
 */

/*
 *  LOAD GMIN
 *
 *  This routine adds Gmin to each diagonal element.  Because Gmin is
 *  added to the current diagonal, which may bear little relation to
 *  what the outside world thinks is a diagonal, and because the
 *  elements that are diagonals may change after calling spOrderAndFactor,
 *  use of this routine is not recommended.  It is included here simply
 *  for compatibility with Spice3.
 */

static void
LoadGmin( eMatrix, Gmin )
char *eMatrix;
register double Gmin;
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;
register int I;
register ArrayOfElementPtrs Diag;
register ElementPtr diag;

/* Begin `LoadGmin'. */
    ASSERT( IS_SPARSE( Matrix ) );

    if (Gmin != 0.0) {
	Diag = Matrix->Diag;
	for (I = Matrix->Size; I > 0; I--) {
	    if (diag = Diag[I])
		diag->Real += Gmin;
	}
    }
    return;
}




/*
 *  FIND ELEMENT
 *
 *  This routine finds an element in the matrix by row and column number.
 *  If the element exists, a pointer to it is returned.  If not, then NULL
 *  is returned unless the CreateIfMissing flag is true, in which case a
 *  pointer to the new element is returned.
 */

SMPelement *
SMPfindElt( eMatrix, Row, Col, CreateIfMissing )

SMPmatrix *eMatrix;
int Row, Col;
int CreateIfMissing;
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;
ElementPtr Element;

/* Begin `SMPfindElt'. */
    ASSERT( IS_SPARSE( Matrix ) );
    if (Matrix->ExtToIntRowMap)
	Row = Matrix->ExtToIntRowMap[Row];
    if (Matrix->ExtToIntColMap)
	Col = Matrix->ExtToIntColMap[Col];
    Element = Matrix->FirstInCol[Col];
    Element = spcFindElementInCol(Matrix, &Element, Row, Col, CreateIfMissing);
    return (SMPelement *)Element;
}

/* XXX The following should probably be implemented in spUtils */

/*
 * SMPcZeroCol()
 */
int
SMPcZeroCol( Matrix, Col )
MatrixPtr Matrix;
int Col;
{
    ElementPtr	Element;

    if (Matrix->ExtToIntColMap)
      Col = Matrix->ExtToIntColMap[Col];

    for (Element = Matrix->FirstInCol[Col];
	Element != NULL;
	Element = Element->NextInCol)
    {
	Element->Real = 0.0;
	Element->Imag = 0.0;
    }

    return spError( (char *)Matrix );
}

/*
 * SMPcAddCol()
 */
int
SMPcAddCol( Matrix, Accum_Col, Addend_Col )
MatrixPtr Matrix;
int Accum_Col, Addend_Col;
{
    ElementPtr	Accum, Addend, *Prev;

    if (Matrix->ExtToIntColMap) {
      Accum_Col = Matrix->ExtToIntColMap[Accum_Col];
      Addend_Col = Matrix->ExtToIntColMap[Addend_Col];
    }

    Addend = Matrix->FirstInCol[Addend_Col];
    Prev = &Matrix->FirstInCol[Accum_Col];
    Accum = *Prev;;

    while (Addend != NULL) {
	while (Accum && Accum->Row < Addend->Row) {
	    Prev = &Accum->NextInCol;
	    Accum = *Prev;
	}
	if (!Accum || Accum->Row > Addend->Row) {
	    Accum = spcCreateElement(Matrix, Addend->Row, Accum_Col, Prev, 0);
	}
	Accum->Real += Addend->Real;
	Accum->Imag += Addend->Imag;
	Addend = Addend->NextInCol;
    }

    return spError( (char *)Matrix );
}

/*
 * SMPzeroRow()
 */
int
SMPzeroRow( Matrix, Row )
MatrixPtr Matrix;
int Row;
{
    ElementPtr	Element;

    if (Matrix->ExtToIntRowMap)
      Row = Matrix->ExtToIntColMap[Row];

    if (Matrix->RowsLinked == NO)
	spcLinkRows(Matrix);

#if spCOMPLEX
    if (Matrix->PreviousMatrixWasComplex OR Matrix->Complex) {
	for (Element = Matrix->FirstInRow[Row];
	    Element != NULL;
	    Element = Element->NextInRow)
	{
	    Element->Real = 0.0;
	    Element->Imag = 0.0;
	}
    } else
#endif
    {
	for (Element = Matrix->FirstInRow[Row];
	    Element != NULL;
	    Element = Element->NextInRow)
	{
            Element->Real = 0.0;
	}
    }

    return spError( (char *)Matrix );
}
/*
 * SMPaddDiagConstraint ()
 */
int
SMPaddDiagConstraint ( Matrix, Row, Value, Rvalue )
MatrixPtr Matrix;
int Row;
double Value, *Rvalue;
{
    ElementPtr	Element;
    double diag_val;

    if (Matrix->ExtToIntRowMap)
      Row = Matrix->ExtToIntColMap[Row];

    if (Matrix->RowsLinked == NO)
      spcLinkRows(Matrix);
    Element = Matrix->Diag[Row];
    if (Element->Real == 0) {
       diag_val = 1;
    }
    else {
       diag_val = Value*Element->Real;
    }

    Element->Real += diag_val;
    *Rvalue = diag_val;

    return 0;
}

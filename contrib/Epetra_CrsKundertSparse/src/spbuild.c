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

#include "spice.h"
#ifdef SHARED_MEM
#include "shared_mem.h"
#endif

/*
#define CHECK_IND
*/

/*
 *  MATRIX BUILD MODULE
 *
 *  Author:                     Advising professor:
 *     Kenneth S. Kundert           Alberto Sangiovanni-Vincentelli
 *     UC Berkeley
 *
 *  This file contains the routines associated with clearing, loading and
 *  preprocessing the matrix for the sparse matrix routines.
 *
 *  >>> User accessible functions contained in this file:
 *  spClear
 *  spGetElement
 *  spGetAdmittance
 *  spGetQuad
 *  spGetOnes
 *  spInstallInitInfo
 *  spGetInitInfo
 *  spInitialize
 *
 *  >>> Other functions contained in this file:
 *  spcFindElementInCol
 *  Translate
 *  spcCreateElement
 *  spcLinkRows
 *  EnlargeMatrix
 *  ExpandTranslationArrays
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





/*
 *  Function declarations
 */

static void Translate( MatrixPtr, int*, int* );
static void EnlargeMatrix( MatrixPtr, int );
static void ExpandTranslationArrays( MatrixPtr, int );
extern int f_ind(MatrixPtr, int, int);
void spCheckInd (MatrixPtr, char *);
void spSetIndex (MatrixPtr);
void spColInd (MatrixPtr, int);
void spRowInd (MatrixPtr, int);
void add_fast_col_index (MatrixPtr, int, int, ElementPtr);
void add_fast_row_index (MatrixPtr, int, int, ElementPtr);



/*
 *  CLEAR MATRIX
 *
 *  Sets every element of the matrix to zero and clears the error flag.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *     Pointer to matrix that is to be cleared.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *     A pointer to the element being cleared.
 */

void
spClear( eMatrix )

char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement;
register  int  I;

/* Begin `spClear'. */
    ASSERT( IS_SPARSE( Matrix ) );
    Matrix->Format = FORMAT_SPARSE;

/* Clear matrix. */
#if spCOMPLEX
    if (Matrix->PreviousMatrixWasComplex OR Matrix->Complex)
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   pElement->Real = 0.0;
                pElement->Imag = 0.0;
                pElement = pElement->NextInCol;
            }
        }
    }
    else
#endif
    {   for (I = Matrix->Size; I > 0; I--)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            {   pElement->Real = 0.0;
                pElement = pElement->NextInCol;
            }
        }
    }

/* Empty the trash. */
    Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
    Matrix->TrashCan.Imag = 0.0;
#endif

    Matrix->Error = spOKAY;
    Matrix->Factored = NO;
    Matrix->SingularCol = 0;
    Matrix->SingularRow = 0;
    Matrix->PreviousMatrixWasComplex = Matrix->Complex;
    return;
}

#ifdef SHARED_MEM
int
spClear_p( eMatrix, rhs, size, ind_d)

char *eMatrix;
double *rhs;
int size, ind_d;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
register  ElementPtr  pElement;
register  int  i, I, n_inc;

/* Begin `spClear'. */
    ASSERT( IS_SPARSE( Matrix ) );

    if (ind_d == 0) return -1;
    for (i=0 ; i<Bins ; i++) {
      num_ind[i] = 0;
    }
    
/* Clear matrix. */
#if spCOMPLEX
    if (Matrix->PreviousMatrixWasComplex OR Matrix->Complex)
    {   for (I = size; I >= 1; I--)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL)
            { n_inc = (((long) &pElement->Real) & Binmask) >> BIN_MASK_SHIFT;
              if (num_ind[n_inc] >= ind_d) return -1;
              ind[n_inc][num_ind[n_inc]++] = &pElement->Real;
              n_inc = (((long) &pElement->Imag) & Binmask) >> BIN_MASK_SHIFT;
              if (num_ind[n_inc] >= ind_d) return -1;
              ind[n_inc][num_ind[n_inc]++] = &pElement->Imag;
              pElement = pElement->NextInCol;
            }
        }
    }
    else
#endif
    {
    for (I = size; I >= 1; I--)
        {   pElement = Matrix->FirstInCol[I];
            while (pElement != NULL) { 
              n_inc = pElement->pe;
              if (n_inc < 0 || n_inc >= num_pes_in_smp)
                n_inc = (((long) &pElement->Real) & Binmask) >> BIN_MASK_SHIFT; 
              if (num_ind[n_inc] >= ind_d) return -1;
              ind[n_inc][num_ind[n_inc]++] = &pElement->Real;
              pElement = pElement->NextInCol;
            }
        }
    }
    return 0;
}

void clean_p (eMatrix)
char *eMatrix;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
/* Empty the trash. */
    Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
    Matrix->TrashCan.Imag = 0.0;
#endif

    Matrix->Error = spOKAY;
    Matrix->Factored = NO;
    Matrix->SingularCol = 0;
    Matrix->SingularRow = 0;
    Matrix->PreviousMatrixWasComplex = Matrix->Complex;
    return;
}
#endif /* SHARED_MEM */












/*
 *  SINGLE ELEMENT ADDITION TO MATRIX BY INDEX
 *
 *  Finds element [Row,Col] and returns a pointer to it.  If element is
 *  not found then it is created and spliced into matrix.  This routine
 *  is only to be used after spCreate() and before spMNA_Preorder(),
 *  spFactor() or spOrderAndFactor().  Returns a pointer to the
 *  Real portion of a MatrixElement.  This pointer is later used by
 *  spADD_xxx_ELEMENT to directly access element.
 *
 *  >>> Returns:
 *  Returns a pointer to the element.  This pointer is then used to directly
 *  access the element during successive builds.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *     Pointer to the matrix that the element is to be added to.
 *  Row  <input>  (int)
 *     Row index for element.  Must be in the range of [0..Size] unless
 *     the options EXPANDABLE or TRANSLATE are used. Elements placed in
 *     row zero are discarded.  In no case may Row be less than zero.
 *  Col  <input>  (int)
 *     Column index for element.  Must be in the range of [0..Size] unless
 *     the options EXPANDABLE or TRANSLATE are used. Elements placed in
 *     column zero are discarded.  In no case may Col be less than zero.
 *
 *  >>> Local variables:
 *  pElement  (RealNumber *)
 *     Pointer to the element.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 *  Error is not cleared in this routine.
 */

RealNumber *
spGetElement( eMatrix, Row, Col )

char *eMatrix;
int  Row, Col;
{
MatrixPtr  Matrix = (MatrixPtr)eMatrix;
RealNumber  *pElement;
ElementPtr spcFindElementInCol();
void  Translate();

/* Begin `spGetElement'. */
    ASSERT( IS_SPARSE( Matrix ) AND Row >= 0 AND Col >= 0 );
    if (Matrix->Format != FORMAT_SPARSE)
      spExpandFormat (Matrix);

    if ((Row == 0) OR (Col == 0))
        return &Matrix->TrashCan.Real;

#if NOT TRANSLATE
    ASSERT(Matrix->NeedsOrdering);
#endif

#if TRANSLATE
    Translate( Matrix, &Row, &Col );
    if (Matrix->Error == spNO_MEMORY) return NULL;
#endif

#if NOT TRANSLATE
#if NOT EXPANDABLE
    ASSERT(Row <= Matrix->Size AND Col <= Matrix->Size);
#endif

#if EXPANDABLE
/* Re-size Matrix if necessary. */
    if ((Row > Matrix->Size) OR (Col > Matrix->Size))
        EnlargeMatrix( Matrix, MAX(Row, Col) );
    if (Matrix->Error == spNO_MEMORY) return NULL;
#endif
#endif

/*
 * The condition part of the following if statement tests to see if the
 * element resides along the diagonal, if it does then it tests to see
 * if the element has been created yet (Diag pointer not NULL).  The
 * pointer to the element is then assigned to Element after it is cast
 * into a pointer to a RealNumber.  This casting makes the pointer into
 * a pointer to Real.  This statement depends on the fact that Real
 * is the first record in the MatrixElement structure.
 */

    if ((Row != Col) OR ((pElement = (RealNumber *)Matrix->Diag[Row]) == NULL))
    {
/*
 * Element does not exist or does not reside along diagonal.  Search
 * column for element.  As in the if statement above, the pointer to the
 * element which is returned by spcFindElementInCol is cast into a
 * pointer to Real, a RealNumber.
 */
        pElement = (RealNumber*)spcFindElementInCol( Matrix,
                                                     &(Matrix->FirstInCol[Col]),
                                                     Row, Col, YES );
    }
    return pElement;
}











/*
 *  FIND ELEMENT BY SEARCHING COLUMN
 *
 *  Searches column starting at element specified at PtrAddr and finds element
 *  in Row. If Element does not exists, it is created. The pointer to the
 *  element is returned.
 *
 *  >>> Returned:
 *  A pointer to the desired element:
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to Matrix.
 *  LastAddr  <input-output>  (ElementPtr *)
 *      Address of pointer that initially points to the element in Col at which
 *      the search is started.  The pointer in this location may be changed if
 *      a fill-in is required in and adjacent element. For this reason it is
 *      important that LastAddr be the address of a FirstInCol or a NextInCol
 *      rather than a temporary variable.
 *  Row  <input>  (int)
 *      Row being searched for.
 *  Col  (int)
 *      Column being searched.
 *  CreateIfMissing  <input>  (BOOLEAN)
 *      Indicates what to do if element is not found, create one or return a
 *      NULL pointer.
 *
 *  Local variables:
 *  pElement  (ElementPtr)
 *      Pointer used to search through matrix.
 */

ElementPtr
spcFindElementInCol( Matrix, LastAddr, Row, Col, CreateIfMissing )

MatrixPtr Matrix;
register ElementPtr *LastAddr;
register int  Row;
int  Col;
BOOLEAN  CreateIfMissing;
{
register  ElementPtr  pElement, pLastElement;
ElementPtr  spcCreateElement();
int i, nfill;

    pElement = Matrix->Col_fast[Col][f_ind(Matrix, Col, Row)];
    if (pElement == NULL || pElement->Row >= Row || pElement->Col != Col) {
      pElement = Matrix->FirstInCol[Col];
      pLastElement = NULL;
    }
    else {
      pLastElement = pElement;
      pElement = pElement->NextInCol;
    }
    while (pElement != NULL && pElement->Row < Row)
    {
      pLastElement = pElement;
      pElement = pElement->NextInCol;
    }

/* Begin `spcFindElementInCol'. */
    if (pElement == NULL || pElement->Row > Row || pElement->Col != Col)
      pElement = *LastAddr; 

/* Search for element. */
    while (pElement != NULL)
    {   if (pElement->Row < Row)
        {
/* Have not reached element yet. */
            LastAddr = &(pElement->NextInCol);
            pElement = pElement->NextInCol;
        }
        else if (pElement->Row == Row)
        {
/* Reached element. */
            return pElement;
        }
        else break;  /* while loop */
    }

/* Element does not exist and must be created. */
    if (CreateIfMissing)
       return spcCreateElement( Matrix, Row, Col, LastAddr, NO );
    else return NULL;
}








#if TRANSLATE

/*
 *  TRANSLATE EXTERNAL INDICES TO INTERNAL
 *
 *  Convert external row and column numbers to internal row and column numbers.
 *  Also updates Ext/Int maps.
 *
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  Row  <input/output>  (int *)
 *     Upon entry Row is either a external row number or an external node
 *     number.  Upon return, the internal equivalent is supplied.
 *  Col  <input/output>  (int *)
 *     Upon entry Column is either a external column number or an external node
 *     number.  Upon return, the internal equivalent is supplied.
 *
 *  >>> Local variables:
 *  ExtCol  (int)
 *     Temporary variable used to hold the external column or node number
 *     during the external to internal column number translation.
 *  ExtRow  (int)
 *     Temporary variable used to hold the external row or node number during
 *     the external to internal row number translation.
 *  IntCol  (int)
 *     Temporary variable used to hold the internal column or node number
 *     during the external to internal column number translation.
 *  IntRow  (int)
 *     Temporary variable used to hold the internal row or node number during
 *     the external to internal row number translation.
 */

static void
Translate( Matrix, Row, Col )

MatrixPtr Matrix;
int  *Row, *Col;
{
register int IntRow, IntCol, ExtRow, ExtCol;

/* Begin `Translate'. */
    ExtRow = *Row;
    ExtCol = *Col;

/* Expand translation arrays if necessary. */
    if ((ExtRow > Matrix->AllocatedExtSize) OR
        (ExtCol > Matrix->AllocatedExtSize))
    {
        ExpandTranslationArrays( Matrix, MAX(ExtRow, ExtCol) );
        if (Matrix->Error == spNO_MEMORY) return;
    }

/* Set ExtSize if necessary. */
    if ((ExtRow > Matrix->ExtSize) OR (ExtCol > Matrix->ExtSize))
        Matrix->ExtSize = MAX(ExtRow, ExtCol);

/* Translate external row or node number to internal row or node number. */
    if ((IntRow = Matrix->ExtToIntRowMap[ExtRow]) == -1)
    {   Matrix->ExtToIntRowMap[ExtRow] = ++Matrix->CurrentSize;
        Matrix->ExtToIntColMap[ExtRow] = Matrix->CurrentSize;
        IntRow = Matrix->CurrentSize;

#if NOT EXPANDABLE
        ASSERT(IntRow <= Matrix->Size);
#endif

#if EXPANDABLE
/* Re-size Matrix if necessary. */
        if (IntRow > Matrix->Size)
            EnlargeMatrix( Matrix, IntRow );
        if (Matrix->Error == spNO_MEMORY) return;
#endif

        Matrix->IntToExtRowMap[IntRow] = ExtRow;
        Matrix->IntToExtColMap[IntRow] = ExtRow;
    }

/* Translate external column or node number to internal column or node number.*/
    if ((IntCol = Matrix->ExtToIntColMap[ExtCol]) == -1)
    {   Matrix->ExtToIntRowMap[ExtCol] = ++Matrix->CurrentSize;
        Matrix->ExtToIntColMap[ExtCol] = Matrix->CurrentSize;
        IntCol = Matrix->CurrentSize;

#if NOT EXPANDABLE
        ASSERT(IntCol <= Matrix->Size);
#endif

#if EXPANDABLE
/* Re-size Matrix if necessary. */
        if (IntCol > Matrix->Size)
            EnlargeMatrix( Matrix, IntCol );
        if (Matrix->Error == spNO_MEMORY) return;
#endif

        Matrix->IntToExtRowMap[IntCol] = ExtCol;
        Matrix->IntToExtColMap[IntCol] = ExtCol;
    }

    *Row = IntRow;
    *Col = IntCol;
    return;
}
#endif






#if QUAD_ELEMENT
/*
 *  ADDITION OF ADMITTANCE TO MATRIX BY INDEX
 *
 *  Performs same function as spGetElement except rather than one
 *  element, all four Matrix elements for a floating component are
 *  added.  This routine also works if component is grounded.  Positive
 *  elements are placed at [Node1,Node2] and [Node2,Node1].  This
 *  routine is only to be used after spCreate() and before
 *  spMNA_Preorder(), spFactor() or spOrderAndFactor().
 *
 *  >>> Returns:
 *  Error code.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *     Pointer to the matrix that component is to be entered in.
 *  Node1  <input>  (int)
 *     Row and column indices for elements. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Node zero is the
 *     ground node.  In no case may Node1 be less than zero.
 *  Node2  <input>  (int)
 *     Row and column indices for elements. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Node zero is the
 *     ground node.  In no case may Node2 be less than zero.
 *  Template  <output>  (struct spTemplate *)
 *     Collection of pointers to four elements that are later used to directly
 *     address elements.  User must supply the template, this routine will
 *     fill it.
 *
 *  Possible errors:
 *  spNO_MEMORY
 *  Error is not cleared in this routine.
 */

int
spGetAdmittance( Matrix, Node1, Node2, Template )

char  *Matrix;
int  Node1, Node2;
struct  spTemplate  *Template;
{

/* Begin `spGetAdmittance'. */
    Template->Element1 = spGetElement(Matrix, Node1, Node1 );
    Template->Element2 = spGetElement(Matrix, Node2, Node2 );
    Template->Element3Negated = spGetElement( Matrix, Node2, Node1 );
    Template->Element4Negated = spGetElement( Matrix, Node1, Node2 );
    if
    (   (Template->Element1 == NULL)
        OR (Template->Element2 == NULL)
        OR (Template->Element3Negated == NULL)
        OR (Template->Element4Negated == NULL)
    )   return spNO_MEMORY;

    if (Node1 == 0)
        SWAP( RealNumber*, Template->Element1, Template->Element2 );

    return spOKAY;
}
#endif /* QUAD_ELEMENT */









#if QUAD_ELEMENT
/*
 *  ADDITION OF FOUR ELEMENTS TO MATRIX BY INDEX
 *
 *  Similar to spGetAdmittance, except that spGetAdmittance only
 *  handles 2-terminal components, whereas spGetQuad handles simple
 *  4-terminals as well.  These 4-terminals are simply generalized
 *  2-terminals with the option of having the sense terminals different
 *  from the source and sink terminals.  spGetQuad adds four
 *  elements to the matrix.  Positive elements occur at Row1,Col1
 *  Row2,Col2 while negative elements occur at Row1,Col2 and Row2,Col1.
 *  The routine works fine if any of the rows and columns are zero.
 *  This routine is only to be used after spCreate() and before
 *  spMNA_Preorder(), spFactor() or spOrderAndFactor()
 *  unless TRANSLATE is set true.
 *
 *  >>> Returns:
 *  Error code.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *     Pointer to the matrix that component is to be entered in.
 *  Row1  <input>  (int)
 *     First row index for elements. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Row1 be less than zero.
 *  Row2  <input>  (int)
 *     Second row index for elements. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Row2 be less than zero.
 *  Col1  <input>  (int)
 *     First column index for elements. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground column.  In no case may Col1 be less than zero.
 *  Col2  <input>  (int)
 *     Second column index for elements. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground column.  In no case may Col2 be less than zero.
 *  Template  <output>  (struct spTemplate *)
 *     Collection of pointers to four elements that are later used to directly
 *     address elements.  User must supply the template, this routine will
 *     fill it.
 *  Real  <input>  (RealNumber)
 *     Real data to be added to elements.
 *  Imag  <input>  (RealNumber)
 *     Imag data to be added to elements.  If matrix is real, this argument
 *     may be deleted.
 *
 *  Possible errors:
 *  spNO_MEMORY
 *  Error is not cleared in this routine.
 */

int
spGetQuad( Matrix, Row1, Row2, Col1, Col2, Template )

char  *Matrix;
int  Row1, Row2, Col1, Col2;
struct  spTemplate  *Template;
{
/* Begin `spGetQuad'. */
    Template->Element1 = spGetElement( Matrix, Row1, Col1);
    Template->Element2 = spGetElement( Matrix, Row2, Col2 );
    Template->Element3Negated = spGetElement( Matrix, Row2, Col1 );
    Template->Element4Negated = spGetElement( Matrix, Row1, Col2 );
    if
    (   (Template->Element1 == NULL)
        OR (Template->Element2 == NULL)
        OR (Template->Element3Negated == NULL)
        OR (Template->Element4Negated == NULL)
    )   return spNO_MEMORY;

    if (Template->Element1 == &((MatrixPtr)Matrix)->TrashCan.Real)
        SWAP( RealNumber *, Template->Element1, Template->Element2 );

    return spOKAY;
}
#endif /* QUAD_ELEMENT */









#if QUAD_ELEMENT
/*
 *  ADDITION OF FOUR STRUCTURAL ONES TO MATRIX BY INDEX
 *
 *  Performs similar function to spGetQuad() except this routine is
 *  meant for components that do not have an admittance representation.
 *
 *  The following stamp is used:
 *         Pos  Neg  Eqn
 *  Pos  [  .    .    1  ]
 *  Neg  [  .    .   -1  ]
 *  Eqn  [  1   -1    .  ]
 *
 *  >>> Returns:
 *  Error code.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *     Pointer to the matrix that component is to be entered in.
 *  Pos  <input>  (int)
 *     See stamp above. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Pos be less than zero.
 *  Neg  <input>  (int)
 *     See stamp above. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Neg be less than zero.
 *  Eqn  <input>  (int)
 *     See stamp above. Must be in the range of [0..Size]
 *     unless the options EXPANDABLE or TRANSLATE are used. Zero is the
 *     ground row.  In no case may Eqn be less than zero.
 *  Template  <output>  (struct spTemplate *)
 *     Collection of pointers to four elements that are later used to directly
 *     address elements.  User must supply the template, this routine will
 *     fill it.
 *
 *  Possible errors:
 *  spNO_MEMORY
 *  Error is not cleared in this routine.
 */

int
spGetOnes(Matrix, Pos, Neg, Eqn, Template)

char  *Matrix;
int  Pos, Neg, Eqn;
struct  spTemplate  *Template;
{
/* Begin `spGetOnes'. */
    Template->Element4Negated = spGetElement( Matrix, Neg, Eqn );
    Template->Element3Negated = spGetElement( Matrix, Eqn, Neg );
    Template->Element2 = spGetElement( Matrix, Pos, Eqn );
    Template->Element1 = spGetElement( Matrix, Eqn, Pos );
    if
    (   (Template->Element1 == NULL)
        OR (Template->Element2 == NULL)
        OR (Template->Element3Negated == NULL)
        OR (Template->Element4Negated == NULL)
    )   return spNO_MEMORY;

    spADD_REAL_QUAD( *Template, 1.0 );
    return spOKAY;
}
#endif /* QUAD_ELEMENT */







/*
 *
 *  CREATE AND SPLICE ELEMENT INTO MATRIX
 *
 *  This routine is used to create new matrix elements and splice them into the
 *  matrix.
 *
 *  >>> Returned:
 *  A pointer to the element that was created is returned.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *  Row  <input>  (int)
 *      Row index for element.
 *  Col  <input>  (int)
 *      Column index for element.
 *  LastAddr  <input-output>  (ElementPtr *)
 *      This contains the address of the pointer to the element just above the
 *      one being created. It is used to speed the search and it is updated with
 *      address of the created element.
 *  Fillin  <input>  (BOOLEAN)
 *      Flag that indicates if created element is to be a fill-in.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix. It is used to refer to the newly
 *      created element and to restring the pointers of the element's row and
 *      column.
 *  pLastElement  (ElementPtr)
 *      Pointer to the element in the matrix that was just previously pointed
 *      to by pElement. It is used to restring the pointers of the element's
 *      row and column.
 *  pCreatedElement  (ElementPtr)
 *      Pointer to the desired element, the one that was just created.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

ElementPtr
spcCreateElement( Matrix, Row, Col, LastAddr, Fillin )

MatrixPtr Matrix;
int  Row;
register int  Col;
register ElementPtr  *LastAddr;
BOOLEAN Fillin;
{
register  ElementPtr  pElement, pLastElement;
ElementPtr  pCreatedElement, spcGetElement(),  spcGetFillin();
int i, nfill;
#ifdef CHECK_IND
char msg[40];
#endif

/* Begin `spcCreateElement'. */

#ifdef CHECK_IND
    sprintf (msg,"Top createdElement (%d,%d)",Row,Col);
    spCheckInd(Matrix, msg);
#endif
#ifdef SHARED_MEM
    new_element = 1;
#endif
    if (Matrix->RowsLinked)
    {
/* Row pointers cannot be ignored. */
        if (Fillin)
        {   pElement = spcGetFillin( Matrix, Row, Col );
            Matrix->Fillins++;
            pElement->Fillin = 1;
        }
        else
        {   pElement = spcGetElement( Matrix, Row, Col );
            Matrix->NeedsOrdering = YES;
            Matrix->NeedsScale = YES;
            Matrix->Elements++;
        }
        if (pElement == NULL) return NULL;

/* If element is on diagonal, store pointer in Diag. */
        if (Row == Col) Matrix->Diag[Row] = pElement;

/* Initialize Element. */
        pCreatedElement = pElement;
        pElement->Row = Row;
        pElement->Col = Col;
        pElement->Real = 0.0;
#if spCOMPLEX
        pElement->Imag = 0.0;
#endif
#if INITIALIZE
        pElement->pInitInfo = NULL;
#endif

/* Splice element into column. */
        pElement->NextInCol = *LastAddr;
        *LastAddr = pElement;

 /* Search row for proper element position. */
        pElement = Matrix->Row_fast[Row][f_ind(Matrix, Row, Col)];
        if (pElement == NULL || pElement->Col >= Col || pElement->Row != Row) {
          pElement = Matrix->FirstInRow[Row];
          pLastElement = NULL;
        }
        else {
          pLastElement = pElement;
          pElement = pElement->NextInRow;
        }
        while (pElement != NULL && pElement->Col < Col)
        {
          pLastElement = pElement;
          pElement = pElement->NextInRow;
        }

/* Splice element into row. */
        pElement = pCreatedElement;
        if (pLastElement == NULL)
        {
/* Element is first in row. */
            pElement->NextInRow = Matrix->FirstInRow[Row];
            Matrix->FirstInRow[Row] = pElement;
        }
        else
/* Element is not first in row. */
        {
            pElement->NextInRow = pLastElement->NextInRow;
            pLastElement->NextInRow = pElement;
        }
        add_fast_row_index (Matrix, Row, Col, pCreatedElement);
    }
    else
    {
/*
 * Matrix has not been factored yet.  Thus get element rather than fill-in.
 * Also, row pointers can be ignored.
 */

/* Allocate memory for Element. */
        pElement = spcGetElement( Matrix, Row, Col );
        if (pElement == NULL) return NULL;

/* If element is on diagonal, store pointer in Diag. */
        if (Row == Col) Matrix->Diag[Row] = pElement;

/* Initialize Element. */
        pCreatedElement = pElement;
        pElement->Row = Row;
#if DEBUG
        pElement->Col = Col;
#endif
        pElement->Real = 0.0;
#if spCOMPLEX
        pElement->Imag = 0.0;
#endif
#if INITIALIZE
        pElement->pInitInfo = NULL;
#endif

/* Splice element into column. */
        pElement->NextInCol = *LastAddr;
        *LastAddr = pElement;
        pCreatedElement->Fillin = 0;
        Matrix->Elements++;
    }

    add_fast_col_index (Matrix, Row, Col, pCreatedElement);

#ifdef CHECK_IND
    sprintf (msg,"Bottom createdElement (%d,%d)",Row,Col);
    spCheckInd(Matrix, msg);
#endif
    return pCreatedElement;
}








/*
 *
 *  LINK ROWS
 *
 *  This routine is used to generate the row links.  The spGetElement()
 *  routines do not create row links, which are needed by the spFactor()
 *  routines.
 *
 *  >>>  Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to the matrix.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      Pointer to an element in the matrix.
 *  FirstInRowEntry  (ElementPtr *)
 *      A pointer into the FirstInRow array.  Points to the FirstInRow entry
 *      currently being operated upon.
 *  FirstInRowArray  (ArrayOfElementPtrs)
 *      A pointer to the FirstInRow array.  Same as Matrix->FirstInRow but
 *      resides in a register and requires less indirection so is faster to
 *      use.
 *  Col  (int)
 *      Column currently being operated upon.
 */

void
spcLinkRows( Matrix )

MatrixPtr Matrix;
{
register  ElementPtr  pElement, *FirstInRowEntry;
register  ArrayOfElementPtrs  FirstInRowArray;
register  int  Col;

/* Begin `spcLinkRows'. */
    FirstInRowArray = Matrix->FirstInRow;
    for (Col = Matrix->Size; Col >= 1; Col--)
    {
/* Generate row links for the elements in the Col'th column. */
        pElement = Matrix->FirstInCol[Col];

        while (pElement != NULL)
        {   pElement->Col = Col;
            FirstInRowEntry = &FirstInRowArray[pElement->Row];
            pElement->NextInRow = *FirstInRowEntry;
            *FirstInRowEntry = pElement;
            pElement = pElement->NextInCol;
        }
    }
    Matrix->RowsLinked = YES;
    spSetIndex(Matrix);
    return;
}








/*
 *  ENLARGE MATRIX
 *
 *  Increases the size of the matrix.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  NewSize  <input>  (int)
 *     The new size of the matrix.
 *
 *  >>> Local variables:
 *  OldAllocatedSize  (int)
 *     The allocated size of the matrix before it is expanded.
 */

static void
EnlargeMatrix( Matrix, NewSize )

MatrixPtr Matrix;
register int  NewSize;
{
int I, OldAllocatedSize = Matrix->AllocatedSize;
ElementPtr  pElement, pLastElement;
int i, j, old_indsize, Row, Col;

/* Begin `EnlargeMatrix'. */
    Matrix->Size = NewSize;

    if (NewSize <= OldAllocatedSize)
        return;

/* Expand the matrix frame. */
    NewSize = MAX( NewSize, EXPANSION_FACTOR * OldAllocatedSize );
    Matrix->AllocatedSize = NewSize;
    old_indsize = Matrix->Indsize;
    Matrix->Indsize = (int) sqrt((double) (IND_DENSITY*2*Matrix->AllocatedSize))+3;

    if (( REALLOC(Matrix->IntToExtColMap, int, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC(Matrix->IntToExtRowMap, int, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC(Matrix->Diag, ElementPtr, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC(Matrix->FirstInCol, ElementPtr, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC(Matrix->FirstInRow, ElementPtr, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC( Matrix->Col_fast, ElementPtr *, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC( Matrix->Row_fast, ElementPtr *, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }

#define CONSERVE_IND

    for (I=1 ; I<=OldAllocatedSize ; I++) {
#ifdef CONSERVE_IND
      REALLOC(Matrix->Col_fast[I], ElementPtr, Matrix->Indsize+1);
      REALLOC(Matrix->Row_fast[I], ElementPtr, Matrix->Indsize+1);
      for (j=old_indsize ; j<= Matrix->Indsize ; j++) {
        Matrix->Col_fast[I][j] = Matrix->Col_fast[I][old_indsize-1];
        Matrix->Row_fast[I][j] = Matrix->Row_fast[I][old_indsize-1];
      }
#else
      FREE(Matrix->Col_fast[I]);
      FREE(Matrix->Row_fast[I]);
#endif
    }
    for (I=OldAllocatedSize+1 ; I<=NewSize ; I++) {
      Matrix->FirstInCol[I] = NULL;
      Matrix->FirstInRow[I] = NULL;
#ifdef CONSERVE_IND
      CALLOC( Matrix->Col_fast[I], ElementPtr, Matrix->Indsize+1);
      CALLOC( Matrix->Row_fast[I], ElementPtr, Matrix->Indsize+1);
#endif
      Matrix->Diag[I] = NULL;
    }
#ifndef CONSERVE_IND
    for (I = 1 ; I<=NewSize ; I++) {
      CALLOC( Matrix->Col_fast[I], ElementPtr, Matrix->Indsize+1);
      if (Matrix->Col_fast[I] == NULL)
      {   Matrix->Error = spNO_MEMORY;
          return;
      }
      CALLOC( Matrix->Row_fast[I], ElementPtr, Matrix->Indsize+1);
      if (Matrix->Row_fast[I] == NULL)
      {   Matrix->Error = spNO_MEMORY;
          return;
      }
    }
    spSetIndex(Matrix);
#endif

/*
 * Destroy the Markowitz and Intermediate vectors, they will be recreated
 * in spOrderAndFactor().
 */
    FREE( Matrix->MarkowitzRow );
    FREE( Matrix->MarkowitzCol );
    FREE( Matrix->MarkowitzProd );
    FREE( Matrix->Nc );
    FREE( Matrix->Nm );
    FREE( Matrix->No );
    FREE( Matrix->DoRealDirect );
    FREE( Matrix->DoCmplxDirect );
    FREE( Matrix->Intermediate );
    FREE( Matrix->Intermediate2 );
    FREE( Matrix->Intermediate3 );
    FREE( Matrix->Intermediate4 );
    Matrix->InternalVectorsAllocated = NO;

/* Initialize the new portion of the vectors. */
    for (I = OldAllocatedSize+1; I <= NewSize; I++)
    {   Matrix->IntToExtColMap[I] = I;
        Matrix->IntToExtRowMap[I] = I;
        Matrix->Diag[I] = NULL;
        Matrix->FirstInRow[I] = NULL;
        Matrix->FirstInCol[I] = NULL;
    }

    return;
}

void spSetIndex (MatrixPtr Matrix)
{
    int i, Row, Col;
    ElementPtr  pElement, pLastElement;

    for (Col=1 ; Col<=Matrix->Size ; Col++) {
      spColInd(Matrix, Col);
    }
    if (Matrix->RowsLinked) {
      for (Row=1 ; Row<=Matrix->Size ; Row++) {
        spRowInd(Matrix, Row);
      }
    }
    return;
}
void spColInd (MatrixPtr Matrix, int Col)
{
    int i;
    ElementPtr  pElement, pLastElement;

    i = 0;
    pElement=Matrix->FirstInCol[Col];
    if (pElement) {
      while (f_ind(Matrix, Col, pElement->Row) >= i)
        Matrix->Col_fast[Col][i++] = NULL;
      pLastElement = pElement;
      for (pElement=pElement->NextInCol ; pElement ; pElement=pElement->NextInCol){
        while (f_ind(Matrix, Col, pElement->Row) >= i) {
          Matrix->Col_fast[Col][i++] = pLastElement;
        }
        pLastElement = pElement;
      }
      while (i<Matrix->Indsize)
        Matrix->Col_fast[Col][i++] = pLastElement;
    }
    return;
}

void spRowInd (MatrixPtr Matrix, int Row)
{
    int i;
    ElementPtr  pElement, pLastElement;

    i = 0;
    pElement=Matrix->FirstInRow[Row];
    if (pElement) {
      while (f_ind(Matrix, Row, pElement->Col) >= i)
        Matrix->Row_fast[Row][i++] = NULL;
      pLastElement = pElement;
      for (pElement=pElement->NextInRow ; pElement ; pElement=pElement->NextInRow){
        while (f_ind(Matrix, Row, pElement->Col) >= i) {
          Matrix->Row_fast[Row][i++] = pLastElement;
        }
        pLastElement = pElement;
      }
      while (i<Matrix->Indsize)
        Matrix->Row_fast[Row][i++] = pLastElement;
    }
    return;
}


#if TRANSLATE

/*
 *  EXPAND TRANSLATION ARRAYS
 *
 *  Increases the size arrays that are used to translate external to internal
 *  row and column numbers.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  NewSize  <input>  (int)
 *     The new size of the translation arrays.
 *
 *  >>> Local variables:
 *  OldAllocatedSize  (int)
 *     The allocated size of the translation arrays before being expanded.
 */

static void
ExpandTranslationArrays( Matrix, NewSize )

MatrixPtr Matrix;
register int  NewSize;
{
register int I, OldAllocatedSize = Matrix->AllocatedExtSize;

/* Begin `ExpandTranslationArrays'. */
    Matrix->ExtSize = NewSize;

    if (NewSize <= OldAllocatedSize)
        return;

/* Expand the translation arrays ExtToIntRowMap and ExtToIntColMap. */
    NewSize = MAX( NewSize, EXPANSION_FACTOR * OldAllocatedSize );
    Matrix->AllocatedExtSize = NewSize;

    if (( REALLOC(Matrix->ExtToIntRowMap, int, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }
    if (( REALLOC(Matrix->ExtToIntColMap, int, NewSize+1)) == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }

/* Initialize the new portion of the vectors. */
    for (I = OldAllocatedSize+1; I <= NewSize; I++)
    {   Matrix->ExtToIntRowMap[I] = -1;
        Matrix->ExtToIntColMap[I] = -1;
    }

    return;
}
#endif









#if INITIALIZE
/*
 *   INITIALIZE MATRIX
 *
 *   With the INITIALIZE compiler option (see spConfig.h) set true,
 *   Sparse allows the user to keep initialization information with each
 *   structurally nonzero matrix element.  Each element has a pointer
 *   that is set and used by the user.  The user can set this pointer
 *   using spInstallInitInfo and may be read using spGetInitInfo.  Both
 *   may be used only after the element exists.  The function
 *   spInitialize() is a user customizable way to initialize the matrix.
 *   Passed to this routine is a function pointer.  spInitialize() sweeps
 *   through every element in the matrix and checks the pInitInfo
 *   pointer (the user supplied pointer).  If the pInitInfo is NULL,
 *   which is true unless the user changes it (almost always true for
 *   fill-ins), then the element is zeroed.  Otherwise, the function
 *   pointer is called and passed the pInitInfo pointer as well as the
 *   element pointer and the external row and column numbers.  If the
 *   user sets the value of each element, then spInitialize() replaces
 *   spClear().
 *
 *   The user function is expected to return a nonzero integer if there
 *   is a fatal error and zero otherwise.  Upon encountering a nonzero
 *   return code, spInitialize() terminates and returns the error code.
 *
 *   >>> Arguments:
 *   Matrix  <input>  (char *)
 *       Pointer to matrix.
 *
 *   >>> Possible Errors:
 *   Returns nonzero if error, zero otherwise.
 */

void
spInstallInitInfo( pElement, pInitInfo )

RealNumber *pElement;
char *pInitInfo;
{
/* Begin `spInstallInitInfo'. */
    ASSERT(pElement != NULL);

    ((ElementPtr)pElement)->pInitInfo = pInitInfo;
}


char *
spGetInitInfo( pElement )

RealNumber *pElement;
{
/* Begin `spGetInitInfo'. */
    ASSERT(pElement != NULL);

    return (char *)((ElementPtr)pElement)->pInitInfo;
}


int
spInitialize( eMatrix, pInit )

char *eMatrix;
int (*pInit)();
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;
register ElementPtr pElement;
int J, Error, Col;

/* Begin `spInitialize'. */
    ASSERT( IS_SPARSE( Matrix ) );

#if spCOMPLEX
/* Clear imaginary part of matrix if matrix is real but was complex. */
    if (Matrix->PreviousMatrixWasComplex AND NOT Matrix->Complex)
    {   for (J = Matrix->Size; J > 0; J--)
        {   pElement = Matrix->FirstInCol[J];
            while (pElement != NULL)
            {   pElement->Imag = 0.0;
                pElement = pElement->NextInCol;
            }
        }
    }
#endif /* spCOMPLEX */

/* Initialize the matrix. */
    for (J = Matrix->Size; J > 0; J--)
    {   pElement = Matrix->FirstInCol[J];
        Col = Matrix->IntToExtColMap[J];
        while (pElement != NULL)
        {   if (pElement->pInitInfo == NULL)
            {   pElement->Real = 0.0;
#               if spCOMPLEX
                    pElement->Imag = 0.0;
#               endif
            }
            else
            {   Error = (*pInit)((RealNumber *)pElement, pElement->pInitInfo,
                                 Matrix->IntToExtRowMap[pElement->Row], Col);
                if (Error)
                {   Matrix->Error = spFATAL;
                    return Error;
                }

            }
            pElement = pElement->NextInCol;
        }
    }

/* Empty the trash. */
    Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
    Matrix->TrashCan.Imag = 0.0;
#endif

    Matrix->Error = spOKAY;
    Matrix->Factored = NO;
    Matrix->SingularCol = 0;
    Matrix->SingularRow = 0;
    Matrix->PreviousMatrixWasComplex = Matrix->Complex;
    return 0;
}
#endif /* INITIALIZE */

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
 *  MATRIX ALLOCATION MODULE
 *
 *  Author:                     Advising professor:
 *      Kenneth S. Kundert          Alberto Sangiovanni-Vincentelli
 *      UC Berkeley
 *
 *  This file contains the allocation and deallocation routines for the
 *  sparse matrix routines.
 *
 *  >>> User accessible functions contained in this file:
 *  spCreate
 *  spDestroy
 *  spError
 *  spWhereSingular
 *  spGetSize
 *  spSetReal
 *  spSetComplex
 *  spFillinCount
 *  spElementCount
 *
 *  >>> Other functions contained in this file:
 *  spcGetElement
 *  InitializeElementBlocks
 *  spcGetFillin
 *  RecordAllocation
 *  AllocateBlockOfAllocationList
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
 *      Macros that customize the sparse matrix routines.
 *  spMatrix.h
 *      Macros and declarations to be imported by the user.
 *  spDefs.h
 *      Matrix type and macro definitions for the sparse matrix routines.
 */

#define spINSIDE_SPARSE
#include "spconfig.h"
#include "spmatrix.h"
#include "spdefs.h"

ElementPtr *returned_elements;
int num_return_cols, *num_returned_elements;

/*
 *  Function declarations
 */

#ifdef __STDC__
static void InitializeElementBlocks( MatrixPtr, int, int );
static void RecordAllocation( MatrixPtr, char* );
static void AllocateBlockOfAllocationList( MatrixPtr );
extern int f_ind(MatrixPtr, int, int);
#else /* __STDC__ */
static void InitializeElementBlocks();
static void RecordAllocation();
static void AllocateBlockOfAllocationList();
#endif /* __STDC__ */




/*
 *  MATRIX ALLOCATION
 *
 *  Allocates and initializes the data structures associated with a matrix.
 *
 *  >>> Returned:
 *  A pointer to the matrix is returned cast into the form of a pointer to
 *  a character.  This pointer is then passed and used by the other matrix
 *  routines to refer to a particular matrix.  If an error occurs, the NULL
 *  pointer is returned.
 *
 *  >>> Arguments:
 *  Size  <input>  (int)
 *      Size of matrix or estimate of size of matrix if matrix is EXPANDABLE.
 *  Complex  <input>  (int)
 *      Type of matrix.  If Complex is 0 then the matrix is real, otherwise
 *      the matrix will be complex.  Note that if the routines are not set up
 *      to handle the type of matrix requested, then a spPANIC error will occur.
 *      Further note that if a matrix will be both real and complex, it must
 *      be specified here as being complex.
 *  pError  <output>  (int *)
 *      Returns error flag, needed because function spError() will not work
 *      correctly if spCreate() returns NULL.
 *
 *  >>> Local variables:
 *  AllocatedSize  (int)
 *      The size of the matrix being allocated.
 *  Matrix  (MatrixPtr)
 *      A pointer to the matrix frame being created.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 *  spPANIC
 *  Error is cleared in this routine.
 */

char *
spCreate( Size, Complex, pError )

int  Size, *pError;
BOOLEAN  Complex;
{
register  unsigned  SizePlusOne;
register  MatrixPtr  Matrix;
register  int  I;
int  AllocatedSize;

/* Begin `spCreate'. */
/* Clear error flag. */
    *pError = spOKAY;

/* Test for valid size. */
    if ((Size < 0) OR (Size == 0 AND NOT EXPANDABLE))
    {   *pError = spPANIC;
        return NULL;
    }

/* Test for valid type. */
#if NOT spCOMPLEX
    if (Complex)
    {   *pError = spPANIC;
        return NULL;
    }
#endif
#if NOT REAL
    if (NOT Complex)
    {   *pError = spPANIC;
        return NULL;
    }
#endif

/* Create Matrix. */
    AllocatedSize = MAX( Size, MINIMUM_ALLOCATED_SIZE );
    SizePlusOne = (unsigned)(AllocatedSize + 1);

    if ((Matrix = ALLOC(struct MatrixFrame, 1)) == NULL)
    {   *pError = spNO_MEMORY;
        return NULL;
    }

/* Initialize matrix */
    Matrix->ID = SPARSE_ID;
    Matrix->Complex = Complex;
    Matrix->PreviousMatrixWasComplex = Complex;
    Matrix->Factored = NO;
    Matrix->Elements = 0;
    Matrix->Error = *pError;
    Matrix->Fillins = 0;
    Matrix->Reordered = NO;
    Matrix->NeedsOrdering = YES;
    Matrix->NeedsScale = YES;
    Matrix->NumberOfInterchangesIsOdd = NO;
    Matrix->Partitioned = NO;
    Matrix->RowsLinked = NO;
    Matrix->InternalVectorsAllocated = NO;
    Matrix->SingularCol = 0;
    Matrix->SingularRow = 0;
    Matrix->Size = Size;
    Matrix->AllocatedSize = AllocatedSize;
    Matrix->ExtSize = Size;
    Matrix->AllocatedExtSize = AllocatedSize;
    Matrix->CurrentSize = 0;
    Matrix->ExtToIntColMap = NULL;
    Matrix->ExtToIntRowMap = NULL;
    Matrix->IntToExtColMap = NULL;
    Matrix->IntToExtRowMap = NULL;
    Matrix->MarkowitzRow = NULL;
    Matrix->MarkowitzCol = NULL;
    Matrix->MarkowitzProd = NULL;
    Matrix->Nc = NULL;
    Matrix->Nm = NULL;
    Matrix->No = NULL;
    Matrix->DoCmplxDirect = NULL;
    Matrix->DoRealDirect = NULL;
    Matrix->Intermediate = NULL;
    Matrix->Intermediate2 = NULL;
    Matrix->Intermediate3 = NULL;
    Matrix->Intermediate4 = NULL;
    Matrix->RelThreshold = DEFAULT_THRESHOLD;
    Matrix->AbsThreshold = 0.0;
    Matrix->Indsize = (int) sqrt((double) (IND_DENSITY*2*AllocatedSize))+3;
    Matrix->Pivots_d = 0;
    Matrix->Pivots = (double *) NULL;
    Matrix->NumReals = NULL;
    Matrix->RealDim = NULL;
    Matrix->DiagPos = NULL;
    Matrix->RealIndex = NULL;
    Matrix->RealValues = NULL;
    Matrix->Format = FORMAT_SPARSE;
    Matrix->DensePointers = 0;

    Matrix->TopOfAllocationList = NULL;
    Matrix->RecordsRemaining = 0;
#ifdef SHARED_MEM
    Matrix->Avgpiv_ratios = (double *) sm_calloc(num_pes_in_smp, sizeof(double));
    Matrix->Minpivs = (double *) sm_calloc(num_pes_in_smp, sizeof(double));
    Matrix->Minpiv_ratios = (double *) sm_calloc(num_pes_in_smp, sizeof(double));
    Matrix->RUpdate = 1;
    Matrix->NextAvailElement = (int *) sm_malloc(num_pes_in_smp*sizeof(int));
    Matrix->ElementsRemaining = (ElementPtr *) sm_malloc(num_pes_in_smp*sizeof(ElementPtr));
#else
    Matrix->ElementsRemaining = 0;
#endif
    Matrix->FillinsRemaining = 0;

    RecordAllocation( Matrix, (char *)Matrix );
    if (Matrix->Error == spNO_MEMORY) goto MemoryError;

#ifdef SHARED_MEM
    Matrix->strips = (struct strip_out **) sm_malloc((num_pes_in_smp*MAX_STRIPS)*sizeof(struct strip_out *));
    Matrix->strips[0] = (struct strip_out *) NULL;
    Matrix->rowStrips = (int *) sm_malloc((num_pes_in_smp*MAX_STRIPS+1)*sizeof(int));
    Matrix->Strips = 0;
    Matrix->MyStuff = (struct context_m *) sm_calloc(num_pes_in_smp, sizeof(struct context_m));
#else
    Matrix->MyStuff = (struct context_m *) tmalloc(sizeof(struct context_m));
#endif
    Matrix->Dense = 0;

/* Take out the trash. */
    Matrix->TrashCan.Real = 0.0;
#if spCOMPLEX
    Matrix->TrashCan.Imag = 0.0;
#endif
    Matrix->TrashCan.Row = 0;
    Matrix->TrashCan.Col = 0;
    Matrix->TrashCan.NextInRow = NULL;
    Matrix->TrashCan.NextInCol = NULL;
#ifdef CHILE
    Matrix->TrashCan.RealDense = &(Matrix->TrashCan.Real);
#endif
#if INITIALIZE
    Matrix->TrashCan.pInitInfo = NULL;
#endif

/* Allocate space in memory for Diag pointer vector. */
    CALLOC( Matrix->Diag, ElementPtr, SizePlusOne);
    if (Matrix->Diag == NULL)
        goto MemoryError;

/* Allocate space in memory for FirstInCol pointer vector. */
    CALLOC( Matrix->FirstInCol, ElementPtr, SizePlusOne);
    if (Matrix->FirstInCol == NULL)
        goto MemoryError;

    CALLOC( Matrix->Col_fast, ElementPtr *, SizePlusOne);
    if (Matrix->Col_fast == NULL)
        goto MemoryError;

    CALLOC( Matrix->Row_fast, ElementPtr *, SizePlusOne);
    if (Matrix->Row_fast == NULL)
        goto MemoryError;

    for (I=1 ; I<=AllocatedSize ; I++) {
      CALLOC( Matrix->Col_fast[I], ElementPtr, Matrix->Indsize+1);
      if (Matrix->Col_fast[I] == NULL)
        goto MemoryError;
      CALLOC( Matrix->Row_fast[I], ElementPtr, Matrix->Indsize+1);
      if (Matrix->Row_fast[I] == NULL)
        goto MemoryError;
    }

/* Allocate space in memory for FirstInRow pointer vector. */
    CALLOC( Matrix->FirstInRow, ElementPtr, SizePlusOne);
    if (Matrix->FirstInRow == NULL)
        goto MemoryError;

/* Allocate space in memory for IntToExtColMap vector. */
    if (( Matrix->IntToExtColMap = ALLOC(int, SizePlusOne)) == NULL)
        goto MemoryError;

/* Allocate space in memory for IntToExtRowMap vector. */
    if (( Matrix->IntToExtRowMap = ALLOC(int, SizePlusOne)) == NULL)
        goto MemoryError;

/* Initialize MapIntToExt vectors. */
    for (I = 1; I <= AllocatedSize; I++)
    {   Matrix->IntToExtRowMap[I] = I;
        Matrix->IntToExtColMap[I] = I;
    }

#if TRANSLATE
/* Allocate space in memory for ExtToIntColMap vector. */
    if (( Matrix->ExtToIntColMap = ALLOC(int, SizePlusOne)) == NULL)
        goto MemoryError;

/* Allocate space in memory for ExtToIntRowMap vector. */
    if (( Matrix->ExtToIntRowMap = ALLOC(int, SizePlusOne)) == NULL)
        goto MemoryError;

/* Initialize MapExtToInt vectors. */
    for (I = 1; I <= AllocatedSize; I++)
    {   Matrix->ExtToIntColMap[I] = -1;
        Matrix->ExtToIntRowMap[I] = -1;
    }
    Matrix->ExtToIntColMap[0] = 0;
    Matrix->ExtToIntRowMap[0] = 0;
#endif

/* Allocate space for fill-ins and initial set of elements. */
    InitializeElementBlocks( Matrix, SPACE_FOR_ELEMENTS*AllocatedSize,
                                     SPACE_FOR_FILL_INS*AllocatedSize );
    if (Matrix->Error == spNO_MEMORY)
        goto MemoryError;

    return (char *)Matrix;

MemoryError:

/* Deallocate matrix and return no pointer to matrix if there is not enough
   memory. */
    *pError = spNO_MEMORY;
    spDestroy( (char *)Matrix);
    return NULL;
}









/*
 *  ELEMENT ALLOCATION
 *
 *  This routine allocates space for matrix elements. It requests large blocks
 *  of storage from the system and doles out individual elements as required.
 *  This technique, as opposed to allocating elements individually, tends to
 *  speed the allocation process.
 *
 *  >>> Returned:
 *  A pointer to an element.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      A pointer to the first element in the group of elements being allocated.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

ElementPtr
spcGetElement(MatrixPtr Matrix, int Row, int Col)
{
ElementPtr  pElement, rElement, *ppElement;
int i, SelectAlloc, Select, pe_memory;

/* Begin `spcGetElement'. */

    Select = Row;
    if (num_return_cols < Select) {
      SelectAlloc = Select+1000;
      num_returned_elements = REALLOC(num_returned_elements, int, SelectAlloc+1);
      returned_elements = REALLOC(returned_elements, ElementPtr, SelectAlloc+1);
      for (i=num_return_cols+1 ; i<=SelectAlloc ; i++) {
        num_returned_elements[i] = 0;
        returned_elements[i] = NULL;
      }
      num_return_cols = SelectAlloc;
    }
    if (num_returned_elements[Select] > 0)
      return (spcGetFillin(Matrix, Row, Col));

#ifdef SHARED_MEM
/* Allocate block of MatrixElements if necessary. */
    pe_memory = 0;
/* If it were possible to allocate memory on a particular PE, then we would use this, and
   pe sensitve calls like sm_malloc_pe()
    pe_memory = get_pe_number(Row, Col);
*/
    if (Matrix->ElementsRemaining[pe_memory] == 0) {
      pElement = (ElementPtr) ALLOC(char, padsize*(ELEMENTS_PER_ALLOCATION+ELEMENTS_PER_CACHE));
      RecordAllocation( Matrix, (char *) pElement );
      if (Matrix->Error == spNO_MEMORY) return NULL;
      Matrix->NextAvailElement[pe_memory] = (ElementPtr) (((((long) pElement) >>
           padshift+LOG_ELEMENTS_PER_CACHE) + 1) << padshift+LOG_ELEMENTS_PER_CACHE);
      Matrix->ElementsRemaining[pe_memory] = ELEMENTS_PER_ALLOCATION;
    }

/* Update Element counter and return pointer to Element. */
    rElement = Matrix->NextAvailElement[pe_memory];
    returned_elements[Select] = NULL;
    ppElement = &returned_elements[Select];
    for (i=0 ; i<ELEMENTS_PER_CACHE ; i++) {
      if (--Matrix->ElementsRemaining[pe_memory] == 0) break;
      pElement = Matrix->NextAvailElement[pe_memory];
      if (i>0) {
        pElement->NextInCol = *ppElement;
        *ppElement = pElement;
        ppElement = &pElement->NextInCol;
        num_returned_elements[Select]++;
      }

      Matrix->NextAvailElement[pe_memory] = (ElementPtr) (((long) Matrix->NextAvailElement[pe_memory]) + padsize);
#else
/* Allocate block of MatrixElements if necessary. */
    if (Matrix->ElementsRemaining == 0) {
      pElement = (ElementPtr) ALLOC(char, padsize*(ELEMENTS_PER_ALLOCATION+ELEMENTS_PER_CACHE));
      RecordAllocation( Matrix, (char *) pElement );
      if (Matrix->Error == spNO_MEMORY) return NULL;
      Matrix->NextAvailElement = (ElementPtr) (((((long) pElement) >>
           padshift+LOG_ELEMENTS_PER_CACHE) + 1) << padshift+LOG_ELEMENTS_PER_CACHE);
      Matrix->ElementsRemaining = ELEMENTS_PER_ALLOCATION;
    }

/* Update Element counter and return pointer to Element. */
    rElement = Matrix->NextAvailElement;
    returned_elements[Select] = NULL;
    ppElement = &returned_elements[Select];
    for (i=0 ; i<ELEMENTS_PER_CACHE ; i++) {
      if (--Matrix->ElementsRemaining == 0) break;
      pElement = Matrix->NextAvailElement;
      if (i>0) {
        pElement->NextInCol = *ppElement;
        *ppElement = pElement;
        ppElement = &pElement->NextInCol;
        num_returned_elements[Select]++;
      }

      Matrix->NextAvailElement = (ElementPtr) (((long) Matrix->NextAvailElement) + padsize);
#endif

/* This alternate method puts as many elements as possible from the chunk into the available
   column.  This saves a small amount of memory, but may make the code slower, especially
   for parallel runs
      Matrix->NextAvailElement++;
      if ((long) Matrix->NextAvailElement > ((long) rElement) + ELEMENTS_PER_CACHE*padsize) {
        Matrix->NextAvailElement = (ElementPtr) (((long) rElement) + ELEMENTS_PER_CACHE*padsize);
        break;
      }
*/

    }
    memset (rElement, 0, sizeof(struct MatrixElement));
#ifdef SHARED_MEM
    rElement->pe = -1;
#endif
    rElement->Row = Row;
    rElement->Col = Col;
#ifdef CHILE
    rElement->RealDense = &(rElement->Real);
#endif
    return rElement;
}








/*
 *  ELEMENT ALLOCATION INITIALIZATION
 *
 *  This routine allocates space for matrix fill-ins and an initial set of
 *  elements.  Besides being faster than allocating space for elements one
 *  at a time, it tends to keep the fill-ins physically close to the other
 *  matrix elements in the computer memory.  This keeps virtual memory paging
 *  to a minimum.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  InitialNumberOfElements  <input> (int)
 *      This number is used as the size of the block of memory, in
 *      MatrixElements, reserved for elements. If more than this number of
 *      elements are generated, then more space is allocated later.
 *  NumberOfFillinsExpected  <input> (int)
 *      This number is used as the size of the block of memory, in
 *      MatrixElements, reserved for fill-ins. If more than this number of
 *      fill-ins are generated, then more space is allocated, but they may
 *      not be physically close in computer's memory.
 *
 *  >>> Local variables:
 *  pElement  (ElementPtr)
 *      A pointer to the first element in the group of elements being allocated.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

static void
InitializeElementBlocks( Matrix, InitialNumberOfElements,
                         NumberOfFillinsExpected )

MatrixPtr Matrix;
int  InitialNumberOfElements, NumberOfFillinsExpected;
{
ElementPtr  pElement;
int i, ColAlloc;

/* Begin `InitializeElementBlocks'. */

    i = 0;
    padsize = sizeof(struct MatrixElement);
    while (padsize) {
      padsize >>= 1;
      i++;
    }
    padshift = i;
    padsize = 1;
    while (i>0) {
      i--;
      padsize <<= 1;
    }

    ColAlloc = 250;
    num_returned_elements = ALLOC(int, ColAlloc+1);
    returned_elements = ALLOC(ElementPtr, ColAlloc+1);
    for (i=num_return_cols+1 ; i<=ColAlloc ; i++) {
      num_returned_elements[i] = 0;
      returned_elements[i] = NULL;
    }

#ifdef SHARED_MEM
    for (i=0 ; i<num_pes_in_smp ; i++) {
      Matrix->ElementsRemaining[i] = 0;
      Matrix->NextAvailElement[i] = NULL;
    }
#else
    Matrix->ElementsRemaining = 0;
    Matrix->NextAvailElement = NULL;
#endif

    return;
}










/*
 *  FILL-IN ALLOCATION
 *
 *  This routine allocates space for matrix fill-ins. It requests large blocks
 *  of storage from the system and doles out individual elements as required.
 *  This technique, as opposed to allocating elements individually, tends to
 *  speed the allocation process.
 *
 *  >>> Returned:
 *  A pointer to the fill-in.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (MatrixPtr)
 *      Pointer to matrix.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

ElementPtr
spcGetFillin(MatrixPtr Matrix, int Row, int Col )
{
    struct FillinListNodeStruct *pListNode;
    ElementPtr  pFillins, rval;
    int Select;

/* Begin `spcGetFillin'. */

    Select = Row;
    if (returned_elements[Select] != NULL) {
      rval = returned_elements[Select];
      returned_elements[Select] = rval->NextInCol;
      num_returned_elements[Select]--;
      if (num_returned_elements[Select] == 0) {
        if (returned_elements[Select] != NULL) {
          printf ("Pointer not Null when count reached zero in spcGetFillin\n");
        }
      }
      memset (rval, 0, sizeof(struct MatrixElement));
#ifdef SHARED_MEM
      rval->pe = -1;
#endif
      rval->Row = Row;
      rval->Col = Col;
      return rval;
    }
    else {
      if (num_returned_elements[Select] != 0) {
          printf ("No returned elements found with num_returned_elements = %d\n",
          num_returned_elements[Select]);
      }
    }

    return spcGetElement( Matrix, Row, Col );
}









/*
 *  RECORD A MEMORY ALLOCATION
 *
 *  This routine is used to record all memory allocations so that the memory
 *  can be freed later.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *  AllocatedPtr  <input>  (char *)
 *      The pointer returned by malloc or calloc.  These pointers are saved in
 *      a list so that they can be easily freed.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

static void
RecordAllocation( Matrix, AllocatedPtr )

MatrixPtr Matrix;
char  *AllocatedPtr;
{
/* Begin `RecordAllocation'. */
/*
 * If Allocated pointer is NULL, assume that malloc returned a NULL pointer,
 * which indicates a spNO_MEMORY error.
 */
    if (AllocatedPtr == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }

/* Allocate block of MatrixElements if necessary. */
    if (Matrix->RecordsRemaining == 0)
    {   AllocateBlockOfAllocationList( Matrix );
        if (Matrix->Error == spNO_MEMORY)
        {   FREE(AllocatedPtr);
            return;
        }
    }

/* Add Allocated pointer to Allocation List. */
    (++Matrix->TopOfAllocationList)->AllocatedPtr = AllocatedPtr;
    Matrix->RecordsRemaining--;
    return;

}








/*
 *  ADD A BLOCK OF SLOTS TO ALLOCATION LIST
 *
 *  This routine increases the size of the allocation list.
 *
 *  >>> Arguments:
 *  Matrix  <input>    (MatrixPtr)
 *      Pointer to the matrix.
 *
 *  >>> Local variables:
 *  ListPtr  (AllocationListPtr)
 *      Pointer to the list that contains the pointers to segments of memory
 *      that were allocated by the operating system for the current matrix.
 *
 *  >>> Possible errors:
 *  spNO_MEMORY
 */

static void
AllocateBlockOfAllocationList( Matrix )

MatrixPtr Matrix;
{
register  int  I;
register  AllocationListPtr  ListPtr;

/* Begin `AllocateBlockOfAllocationList'. */
/* Allocate block of records for allocation list. */
    ListPtr = ALLOC(struct AllocationRecord, (ELEMENTS_PER_ALLOCATION+1));
    if (ListPtr == NULL)
    {   Matrix->Error = spNO_MEMORY;
        return;
    }

/* String entries of allocation list into singly linked list.  List is linked
   such that any record points to the one before it. */

    ListPtr->NextRecord = Matrix->TopOfAllocationList;
    Matrix->TopOfAllocationList = ListPtr;
    ListPtr += ELEMENTS_PER_ALLOCATION;
    for (I = ELEMENTS_PER_ALLOCATION; I > 0; I--)
    {    ListPtr->NextRecord = ListPtr - 1;
         ListPtr--;
    }

/* Record allocation of space for allocation list on allocation list. */
    Matrix->TopOfAllocationList->AllocatedPtr = (char *)ListPtr;
    Matrix->RecordsRemaining = ELEMENTS_PER_ALLOCATION;

    return;
}








/*
 *  MATRIX DEALLOCATION
 *
 *  Deallocates pointers and elements of Matrix.
 *
 *  >>> Arguments:
 *  Matrix  <input>  (char *)
 *      Pointer to the matrix frame which is to be removed from memory.
 *
 *  >>> Local variables:
 *  ListPtr  (AllocationListPtr)
 *      Pointer into the linked list of pointers to allocated data structures.
 *      Points to pointer to structure to be freed.
 *  NextListPtr  (AllocationListPtr)
 *      Pointer into the linked list of pointers to allocated data structures.
 *      Points to the next pointer to structure to be freed.  This is needed
 *      because the data structure to be freed could include the current node
 *      in the allocation list.
 */

void
spDestroy( eMatrix )

register char *eMatrix;
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;
register  AllocationListPtr  ListPtr, NextListPtr;
int I;


/* Begin `spDestroy'. */
    ASSERT( IS_SPARSE( Matrix ) );

/* Deallocate the vectors that are located in the matrix frame. */
    FREE( Matrix->IntToExtColMap );
    FREE( Matrix->IntToExtRowMap );
    FREE( Matrix->ExtToIntColMap );
    FREE( Matrix->ExtToIntRowMap );
    FREE( Matrix->Diag );
    FREE( Matrix->FirstInRow );
    FREE( Matrix->FirstInCol );
    FREE( Matrix->MarkowitzRow );
    FREE( Matrix->MarkowitzCol );
    FREE( Matrix->MarkowitzProd );
    FREE( Matrix->Nc );
    FREE( Matrix->Nm );
    FREE( Matrix->No );
    FREE( Matrix->DoCmplxDirect );
    FREE( Matrix->DoRealDirect );
    FREE( Matrix->Intermediate );
    FREE( Matrix->Intermediate2 );
    FREE( Matrix->Intermediate3 );
    for (I=1 ; I<=Matrix->Size ; I++) {
      FREE ( Matrix->Col_fast[I]);
      FREE ( Matrix->Row_fast[I]);
    }
    FREE (Matrix->Col_fast );
    FREE (Matrix->Row_fast );

/* Sequentially step through the list of allocated pointers freeing pointers
 * along the way. */
    ListPtr = Matrix->TopOfAllocationList;
    while (ListPtr != NULL)
    {   NextListPtr = ListPtr->NextRecord;
	if ((char *) ListPtr == ListPtr->AllocatedPtr)
	{
	    FREE( ListPtr );
	}
	else
	{
	    FREE( ListPtr->AllocatedPtr );
	}
        ListPtr = NextListPtr;
    }
    return;
}







/*
 *  RETURN MATRIX ERROR STATUS
 *
 *  This function is used to determine the error status of the given matrix.
 *
 *  >>> Returned:
 *      The error status of the given matrix.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      The matrix for which the error status is desired.
 */

int
spError( eMatrix )

char  *eMatrix;
{
/* Begin `spError'. */

    if (eMatrix != NULL)
    {   ASSERT(((MatrixPtr)eMatrix)->ID == SPARSE_ID);
        return ((MatrixPtr)eMatrix)->Error;
    }
    else return spNO_MEMORY;   /* This error may actually be spPANIC,
                                * no way to tell. */
}









/*
 *  WHERE IS MATRIX SINGULAR
 *
 *  This function returns the row and column number where the matrix was
 *  detected as singular or where a zero was detected on the diagonal.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      The matrix for which the error status is desired.
 *  pRow  <output>  (int *)
 *      The row number.
 *  pCol  <output>  (int *)
 *      The column number.
 */

void
spWhereSingular( eMatrix, pRow, pCol )

char *eMatrix;
int *pRow, *pCol;
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;

/* Begin `spWhereSingular'. */
    ASSERT( IS_SPARSE( Matrix ) );

    if (Matrix->Error == spSINGULAR OR Matrix->Error == spZERO_DIAG)
    {   *pRow = Matrix->SingularRow;
        *pCol = Matrix->SingularCol;
    }
    else *pRow = *pCol = 0;
    return;
}






/*
 *  MATRIX SIZE
 *
 *  Returns the size of the matrix.  Either the internal or external size of
 *  the matrix is returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to matrix.
 *  External  <input>  (BOOLEAN)
 *      If External is set true, the external size , i.e., the value of the
 *      largest external row or column number encountered is returned.
 *      Otherwise the true size of the matrix is returned.  These two sizes
 *      may differ if the TRANSLATE option is set true.
 */

int
spGetSize( eMatrix, External )

char  *eMatrix;
BOOLEAN  External;
{
MatrixPtr Matrix = (MatrixPtr)eMatrix;

/* Begin `spGetSize'. */
    ASSERT( IS_SPARSE( Matrix ) );

#if TRANSLATE
    if (External)
        return Matrix->ExtSize;
    else
        return Matrix->Size;
#else
    return Matrix->Size;
#endif
}








/*
 *  SET MATRIX COMPLEX OR REAL
 *
 *  Forces matrix to be either real or complex.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to matrix.
 */

void
spSetReal( eMatrix )

char *eMatrix;
{
/* Begin `spSetReal'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) AND REAL);
    ((MatrixPtr)eMatrix)->Complex = NO;
    return;
}


void
spSetComplex( eMatrix )

char  *eMatrix;
{
/* Begin `spSetComplex'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) AND spCOMPLEX);
    ((MatrixPtr)eMatrix)->Complex = YES;
    return;
}









/*
 *  ELEMENT OR FILL-IN COUNT
 *
 *  Two functions used to return simple statistics.  Either the number
 *  of total elements, or the number of fill-ins can be returned.
 *
 *  >>> Arguments:
 *  eMatrix  <input>  (char *)
 *      Pointer to matrix.
 */

int
spFillinCount( eMatrix )

char *eMatrix;
{
/* Begin `spFillinCount'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) );
    return ((MatrixPtr)eMatrix)->Fillins;
}


int
spElementCount( eMatrix )

char  *eMatrix;
{
/* Begin `spElementCount'. */

    ASSERT( IS_SPARSE( (MatrixPtr)eMatrix ) );
    return ((MatrixPtr)eMatrix)->Elements;
}

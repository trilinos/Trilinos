#ifndef BOOLEAN
#define BOOLEAN int
#endif
#ifndef RealNumber
#define RealNumber double
#endif
#include "sppars.h"

/*
 *  MATRIX ELEMENT DATA STRUCTURE
 *
 *  Every nonzero element in the matrix is stored in a dynamically allocated
 *  MatrixElement structure.  These structures are linked together in an
 *  orthogonal linked list.  Two different MatrixElement structures exist.
 *  One is used when only real matrices are expected, it is missing an entry
 *  for imaginary data.  The other is used if complex matrices are expected.
 *  It contains an entry for imaginary data.
 *
 *  >>> Structure fields:
 *  Real  (RealNumber)
 *      The real portion of the value of the element.  Real must be the first
 *      field in this structure.
 *  Imag  (RealNumber)
 *      The imaginary portion of the value of the element. If the matrix
 *      routines are not compiled to handle complex matrices, then this
 *      field does not exist.  If it exists, it must follow immediately after
 *      Real.
 *  Row  (int)
 *      The row number of the element.
 *  Col  (int)
 *      The column number of the element.
 *  NextInRow  (struct MatrixElement *)
 *      NextInRow contains a pointer to the next element in the row to the
 *      right of this element.  If this element is the last nonzero in the
 *      row then NextInRow contains NULL.
 *  NextInCol  (struct MatrixElement *)
 *      NextInCol contains a pointer to the next element in the column below
 *      this element.  If this element is the last nonzero in the column then
 *      NextInCol contains NULL.
 *  pInitInfo  (char *)
 *      Pointer to user data used for initialization of the matrix element.
 *      Initialized to NULL.
 *
 *  >>> Type definitions:
 *  ElementPtr
 *      A pointer to a MatrixElement.
 *  ArrayOfElementPtrs
 *      An array of ElementPtrs.  Used for FirstInRow, FirstInCol and
 *      Diag pointer arrays.
 */

/* Begin `MatrixElement'. */

/* Note: map_SM (in main.c) requires that 'Real' be the first entry
   in this structure. */
struct  MatrixElement
{   RealNumber   Real;
#if spCOMPLEX
    RealNumber   Imag;
#endif
#ifdef CHILE
    RealNumber   RealBackup;
    RealNumber   RealCopy;
    RealNumber  *RealDense;
#endif
#ifdef SHARED_MEM
    int          pe;
#endif
    int          Fillin;
    int          Row;
    int          Col;
    struct MatrixElement  *NextInRow;
    struct MatrixElement  *NextInCol;
#if INITIALIZE
    char        *pInitInfo;
#endif
};

typedef  struct MatrixElement  *ElementPtr;
typedef  ElementPtr  *ArrayOfElementPtrs;








/*
 *  ALLOCATION DATA STRUCTURE
 *
 *  The sparse matrix routines keep track of all memory that is allocated by
wapTask **
 *  saving the pointers to all the chunks of memory that are allocated to a
 *  particular matrix in an allocation list.  That list is organized as a
 *  linked list so that it can grow without a priori bounds.
 *
 *  >>> Structure fields:
 *  AllocatedPtr  (char *)
 *      Pointer to chunk of memory that has been allocated for the matrix.
 *  NextRecord  (struct  AllocationRecord *)
 *      Pointer to the next allocation record.
 */

/* Begin `AllocationRecord'. */
struct AllocationRecord
{   char  *AllocatedPtr;
    struct  AllocationRecord  *NextRecord;
};

typedef  struct  AllocationRecord  *AllocationListPtr;









/*
 *  FILL-IN LIST DATA STRUCTURE
 *
 *  The sparse matrix routines keep track of all fill-ins separately from
 *  user specified elements so they may be removed by spStripFills().  Fill-ins
 *  are allocated in bunched in what is called a fill-in lists.  The data
 *  structure defined below is used to organize these fill-in lists into a
 *  linked-list.
 *
 *  >>> Structure fields:
 *  pFillinList  (ElementPtr)
 *      Pointer to a fill-in list, or a bunch of fill-ins arranged contiguously
 *      in memory.
 *  NumberOfFillinsInList  (int)
 *      Seems pretty self explanatory to me.
 *  Next  (struct  FillinListNodeStruct *)
 *      Pointer to the next fill-in list structures.
 */

/* Begin `FillinListNodeStruct'. */
struct FillinListNodeStruct
{   ElementPtr  pFillinList;
    int         NumberOfFillinsInList;
    struct      FillinListNodeStruct  *Next;
};










/*
 *  MATRIX FRAME DATA STRUCTURE
 *
 *  This structure contains all the pointers that support the orthogonal
 *  linked list that contains the matrix elements.  Also included in this
 *  structure are other numbers and pointers that are used globally by the
 *  sparse matrix routines and are associated with one particular matrix.
 *
 *  >>> Type definitions:
 *  MatrixPtr
 *      A pointer to MatrixFrame.  Essentially, a pointer to the matrix.
 *
 *  >>> Structure fields:
 *  AbsThreshold  (RealNumber)
 *      The absolute magnitude an element must have to be considered as a
 *      pivot candidate, except as a last resort.
 *  AllocatedExtSize  (int)
 *      The allocated size of the arrays used to translate external row and
 *      column numbers to their internal values.
 *  AllocatedSize  (int)
 *      The currently allocated size of the matrix; the size the matrix can
 *      grow to when EXPANDABLE is set true and AllocatedSize is the largest
 *      the matrix can get without requiring that the matrix frame be
 *      reallocated.
 *  Complex  (BOOLEAN)
 *      The flag which indicates whether the matrix is complex (true) or
 *      real.
 *  CurrentSize  (int)
 *      This number is used during the building of the matrix when the
 *      TRANSLATE option is set true.  It indicates the number of internal
 *      rows and columns that have elements in them.
 *  Diag  (ArrayOfElementPtrs)
 *      Array of pointers that points to the diagonal elements.
 *  DoCmplxDirect  (BOOLEAN *)
 *      Array of flags, one for each column in matrix.  If a flag is true
 *      then corresponding column in a complex matrix should be eliminated
 *      in spFactor() using direct addressing (rather than indirect
 *      addressing).
 *  DoRealDirect  (BOOLEAN *)
 *      Array of flags, one for each column in matrix.  If a flag is true
 *      then corresponding column in a real matrix should be eliminated
 *      in spFactor() using direct addressing (rather than indirect
 *      addressing).
 *  Elements  (int)
 *      The number of original elements (total elements minus fill ins)
 *      present in matrix.
 *  Error  (int)
 *      The error status of the sparse matrix package.
 *  ExtSize  (int)
 *      The value of the largest external row or column number encountered.
 *  ExtToIntColMap  (int [])
 *      An array that is used to convert external columns number to internal
 *      external column numbers.  Present only if TRANSLATE option is set true.
 *  ExtToIntRowMap  (int [])
 *      An array that is used to convert external row numbers to internal
 *      external row numbers.  Present only if TRANSLATE option is set true.
 *  Factored  (BOOLEAN)
 *      Indicates if matrix has been factored.  This flag is set true in
 *      spFactor() and spOrderAndFactor() and set false in spCreate()
 *      and spClear().
 *  Fillins  (int)
 *      The number of fill-ins created during the factorization the matrix.
 *  FirstInCol  (ArrayOfElementPtrs)
 *      Array of pointers that point to the first nonzero element of the
 *      column corresponding to the index.
 *  FirstInRow  (ArrayOfElementPtrs)
 *      Array of pointers that point to the first nonzero element of the row
 *      corresponding to the index.
 *  ID  (unsigned long int)
 *      A constant that provides the sparse data structure with a signature.
 *      When DEBUG is true, all externally available sparse routines check
 *      this signature to assure they are operating on a valid matrix.
 *  Intermediate  (RealVector)
 *      Temporary storage used in the spSolve routines. Intermediate is an
 *      array used during forward and backward substitution.  It is
 *      commonly called y when the forward and backward substitution process is
 *      denoted  Ax = b => Ly = b and Ux = y.
 *  InternalVectorsAllocated  (BOOLEAN)
 *      A flag that indicates whether the Markowitz vectors and the
 *      Intermediate vector have been created.
 *      These vectors are created in spcCreateInternalVectors().
 *  IntToExtColMap  (int [])
 *      An array that is used to convert internal column numbers to external
 *      external column numbers.
 *  IntToExtRowMap  (int [])
 *      An array that is used to convert internal row numbers to external
 *      external row numbers.
 *  MarkowitzCol  (int [])
 *      An array that contains the count of the non-zero elements excluding
 *      the pivots for each column. Used to generate and update MarkowitzProd.
 *  MarkowitzProd  (long [])
 *      The array of the products of the Markowitz row and column counts. The
 *      element with the smallest product is the best pivot to use to maintain
 *      sparsity.
 *  MarkowitzRow  (int [])
 *      An array that contains the count of the non-zero elements excluding
 *      the pivots for each row. Used to generate and update MarkowitzProd.
 *  MaxRowCountInLowerTri  (int)
 *      The maximum number of off-diagonal element in the rows of L, the
 *      lower triangular matrix.  This quantity is used when computing an
 *      estimate of the roundoff error in the matrix.
 *  Max_TS (double)
 *      The maximum time step since last pivot determination.
 *  NeedsOrdering  (BOOLEAN)
 *      This is a flag that signifies that the matrix needs to be ordered
 *      or reordered.  NeedsOrdering is set true in spCreate() and
 *      spGetElement() or spGetAdmittance() if new elements are added to the
 *      matrix after it has been previously factored.  It is set false in
 *      spOrderAndFactor().
 *  NumberOfInterchangesIsOdd  (BOOLEAN)
 *      Flag that indicates the sum of row and column interchange counts
 *      is an odd number.  Used when determining the sign of the determinant.
 *  Partitioned  (BOOLEAN)
 *      This flag indicates that the columns of the matrix have been 
 *      partitioned into two groups.  Those that will be addressed directly
 *      and those that will be addressed indirectly in spFactor().
 *  PivotsOriginalCol  (int)
 *      Column pivot was chosen from.
 *  PivotsOriginalRow  (int)
 *      Row pivot was chosen from.
 *  PivotSelectionMethod  (char)
 *      Character that indicates which pivot search method was successful.
 *  PreviousMatrixWasComplex  (BOOLEAN)
 *      This flag in needed to determine how to clear the matrix.  When
 *      dealing with real matrices, it is important that the imaginary terms
 *      in the matrix elements be zero.  Thus, if the previous matrix was
 *      complex, then the current matrix will be cleared as if it were complex
 *      even if it is real.
 *  RelThreshold  (RealNumber)
 *      The magnitude an element must have relative to others in its row
 *      to be considered as a pivot candidate, except as a last resort.
 *  Reordered  (BOOLEAN)
 *      This flag signifies that the matrix has been reordered.  It
 *      is cleared in spCreate(), set in spMNA_Preorder() and
 *      spOrderAndFactor() and is used in spPrint().
 *  RowsLinked  (BOOLEAN)
 *      A flag that indicates whether the row pointers exist.  The AddByIndex
 *      routines do not generate the row pointers, which are needed by some
 *      of the other routines, such as spOrderAndFactor() and spScale().
 *      The row pointers are generated in the function spcLinkRows().
 *  SingularCol  (int)
 *      Normally zero, but if matrix is found to be singular, SingularCol is
 *      assigned the external column number of pivot that was zero.
 *  SingularRow  (int)
 *      Normally zero, but if matrix is found to be singular, SingularRow is
 *      assigned the external row number of pivot that was zero.
 *  Singletons  (int)
 *      The number of singletons available for pivoting.  Note that if row I
 *      and column I both contain singletons, only one of them is counted.
 *  Size  (int)
 *      Number of rows and columns in the matrix.  Does not change as matrix
 *      is factored.
 *  TrashCan  (MatrixElement)
 *      This is a dummy MatrixElement that is used to by the user to stuff
 *      data related to the zero row or column.  In other words, when the user
 *      adds an element in row zero or column zero, then the matrix returns
 *      a pointer to TrashCan.  In this way the user can have a uniform way
 *      data into the matrix independent of whether a component is connected
 *      to ground.
 *
 *  >>> The remaining fields are related to memory allocation.
 *  TopOfAllocationList  (AllocationListPtr)
 *      Pointer which points to the top entry in a list. The list contains
 *      all the pointers to the segments of memory that have been allocated
 *      to this matrix. This is used when the memory is to be freed on
 *      deallocation of the matrix.
 *  RecordsRemaining  (int)
 *      Number of slots left in the list of allocations.
 *  NextAvailElement  (ElementPtr)
 *      Pointer to the next available element which has been allocated but as
 *      yet is unused. Matrix elements are allocated in groups of
 *      ELEMENTS_PER_ALLOCATION in order to speed element allocation and
 *      freeing.
 *  ElementsRemaining  (int)
 *      Number of unused elements left in last block of elements allocated.
 *  NextAvailFillin  (ElementPtr)
 *      Pointer to the next available fill-in which has been allocated but
 *      as yet is unused.  Fill-ins are allocated in a group in order to keep
 *      them physically close in memory to the rest of the matrix.
 *  FillinsRemaining  (int)
 *      Number of unused fill-ins left in the last block of fill-ins
 *      allocated.
 *  FirstFillinListNode  (FillinListNodeStruct *)
 *      A pointer to the head of the linked-list that keeps track of the
 *      lists of fill-ins.
 *  LastFillinListNode  (FillinListNodeStruct *)
 *      A pointer to the tail of the linked-list that keeps track of the
 *      lists of fill-ins.
 *  Column_pointers   -- added for compressed sparse column representation
 *  Row_indices          See spice3f5/src/lib/sparse/lltocsc.c
 *  Values
 */
#ifdef SHARED_MEM
struct pivcol {
    double pivot;
    int col;
    int col0;
};

struct strip_out {
    ElementPtr *FirstInCol;
    int *update;
    int jcop;
    struct pivcol *pc;
    struct pivcol **pc_out;
    struct pivcol **pc_gen;
    int len;
    int done;
    struct strip_out *my_next;
    struct strip_out *prev;
};
#endif

struct context_m {
    int                          Dsize;
    int                          BufDim;
    int                          BufUsed;
    int                         *ColStart_s;
    int                         *ColDiag;
    int                         *MyI;
    double                      *MyD;
    double                      *Dest;
    unsigned char               *ind_list;
    int                         *ind_list_i;
    int                          ind_list_d;
    int                          ind_list_i_d;
    int                         *ColStart[MAX_STRIPS];
};

/* Begin `MatrixFrame'. */
struct  MatrixFrame
{   RealNumber                   AbsThreshold;
    int                          AllocatedSize;
    int                          AllocatedExtSize;
    double                      *Avgpiv_ratios;
    void                        *Ckt;
    BOOLEAN                      Complex;
    int                          CurrentSize;
    ArrayOfElementPtrs           Diag;
    int                         *DiagPos;
    BOOLEAN                     *DoCmplxDirect;
    BOOLEAN                     *DoRealDirect;
    int                          Elements;
    int                          Error;
    int                          ExtSize;
    int                         *ExtToIntColMap;
    int                         *ExtToIntRowMap;
    BOOLEAN                      Factored;
    int                          Fillins;
    ArrayOfElementPtrs           FirstInCol;
    ArrayOfElementPtrs           FirstInRow;
    int                          Format;
    int                          DensePointers;
    ArrayOfElementPtrs          *Col_fast;
    ArrayOfElementPtrs          *Row_fast;
    int                          Hi_lim;
    int                          Curr_lim;
    unsigned long                ID;
    int                          Indsize;
    double                      *Intermediate;
    double                      *Intermediate2;
    double                      *Intermediate3;
    double                      *Intermediate4;
    BOOLEAN                      InternalVectorsAllocated;
    int                         *IntToExtColMap;
    int                         *IntToExtRowMap;
    int                          LastIterationCount;
    int                         *MarkowitzRow;
    int                         *MarkowitzCol;
    long                        *MarkowitzProd;
    int                          MaxRowCountInLowerTri;
    double                       Max_TS;
    double                      *Minpivs;
    double                      *Minpiv_ratios;
    int                         *Nc;
    int                         *Nm;
    int                         *No;
    BOOLEAN                      NeedsOrdering;
    BOOLEAN                      NeedsScale;
    int                         *NumReals;
    BOOLEAN                      NumberOfInterchangesIsOdd;
    int                          OverflowDanger;
    BOOLEAN                      Partitioned;
    int                          PivotRefining;
    double                      *Pivots;
    int                          Pivots_d;
    int                          PivotsOriginalCol;
    int                          PivotsOriginalRow;
    char                         PivotSelectionMethod;
    BOOLEAN                      PreviousMatrixWasComplex;
    int                         *RealDim;
    int                        **RealIndex;
    double                     **RealValues;
    RealNumber                   RelThreshold;
    BOOLEAN                      Reordered;
    BOOLEAN                      RowsLinked;
    int                          SingularCol;
    int                          SingularRow;
    int                          Singletons;
    int                          SmallTimeStep;
    int                          Size;
    struct MatrixElement         TrashCan;
    int                          Updated;
    struct context_m            *MyStuff;
    int                          Dense;
    int                          NewFlags;
    int                          RePart;
    int                          has_scale_factors;
    int                          scale_factors_d;
    double                      *row_scale_factors;
    double                      *col_scale_factors;

    AllocationListPtr            TopOfAllocationList;
    int                          RecordsRemaining;
    ElementPtr                   NextAvailFillin;
    int                          FillinsRemaining;
    struct FillinListNodeStruct *FirstFillinListNode;
    struct FillinListNodeStruct *LastFillinListNode;
    int                         *Column_pointers;
    int                         *Row_indices;
    double                      *Values;
    int                          Cscflag;
#ifdef SHARED_MEM
    int                          RUpdate;
    struct strip_out           **strips;
    int                         *rowStrips;
    int                          Strips;
    ElementPtr                   *NextAvailElement;
    int                          *ElementsRemaining;
#else
    ElementPtr                   NextAvailElement;
    int                          ElementsRemaining;
#endif
};
typedef  struct MatrixFrame  *MatrixPtr;


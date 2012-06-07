/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#ifndef ML_QR_FIX_H
#define ML_QR_FIX_H
/*
  This header allows ML to handle cases where an aggregate is not large enough
  to support the given number of nullspace vectors.  Take the following *scalar*
  example:

  aggregate   :  DOFs
       0      :  0 1 2 3
       1      :  4 5 6
       2      :  7 8

  with nullspace of size 4.  The tentative prolongator should look something like

          0     x x x x
          1     x x x x
          2     x x x x
          3     x x x x
          4             x x x x
          5             x x x x
          6             x x x x
          7                    x x x x
          8                    x x x x

  Aggregate 0 is large enough that locally 4 linearly independent vectors can
  be represented.  However, aggregate 1 can represent at most 3 of the nullspace vectors and still
  maintain linear independence (LI).  Aggregate 2 can only represent 2 LI vectors.  This lack of
  LI will come out in the QR factorization.

  We have two choices in how to handle this.  First, we could drop columns, but then bookkeeping gets
  complicated.  Second, we can somehow mark columns that are invalid. ML implements the second option
  as follows.

  ML maintains a "dead DOF" data structure for each multigrid level.  The data structure keeps
  the level number, the total number of aggregtes, the total number of dead DOFs, and an array of chars or
  integers that maybe short, unsigned, or long.  The array has one entry per aggregate.  This entry
  is treated as an array of bits, with the first bit (bit 0) being the rightmost.  If nullspace
  vector i can be represented in the current aggregate, then bit i is set to 0.  If vector i *cannot* be
  represented, then bit i is set to 1.

  So for the above example, the array of (char) entries looks like the following

  00000000
  00001000
  00001100

  Why use bits?  First, it's very space efficient.  Second, it's very easy to populate each entry e as follows:

    for myAgg=0:2
      e=0;
      for i= size(myAgg) : size(nullSpace)-1
        e |= (1<<i);
      end
    end

  which will give the three entries above.

  Finally, it'also easy to figure out which DOFs are "dead":

  for i = 0: size(nullspace)-1
     if (dead & (1 << k))
       printf("DOF %d is dead\n",i);
     end
  end

  In fact, these two operations are used in ml_agg_uncoupled.c, ml_agg_MIS.c, and ml_struct.c.

  Observations:
  1) Only the rightmost 4 bits matter here.  There's room in fact to encode information for 4 more
     (for a total of 8) nullspace vectors.

  2) On my Linux workstation, a long int is 64 bits.  This is the largest nullspace that can be encoded
     using this approach, as only integers can be bit shifted.  So if you have nullspace with dimension > 64,
     you're out of luck with the current scheme.


  This encoding scheme is used in two spots.  The first is to zero out (local to an aggregate) NS vectors that
  cannot be represented because the aggregate is too small.  The second is to fix up the coarse grid matrix A
  so that it does not have a zero on the diagonal.

*/

/* If we need more than 16 kernel components, define ML_QR_FIX_TYPE
 * as unsigned int, otherwise use unsigned short int to conserve memory */
#define ML_QR_FIX_TYPE long int

typedef struct ML_qr_fix {

  int                 level;
  int                 numDeadNodDof;
 /* -mb: can later replace the following two with a hash structure */ 
  int                 nDeadNodDof;     /*=number of aggregates*/
  ML_QR_FIX_TYPE     *xDeadNodDof;

} ML_qr_fix;

#ifdef __cplusplus
extern "C" {
  int ML_qr_fix_Create(const int nCoarseNod);

  int ML_qr_fix_Destroy(void);

  int ML_qr_fix_Print(ML_qr_fix* ptr);

  int ML_qr_fix_NumDeadNodDof(void);

  ML_QR_FIX_TYPE ML_qr_fix_getDeadNod(const int inx);

  void ML_qr_fix_setNumDeadNod(int num);

  void ML_qr_fix_setDeadNod( const int inx, ML_QR_FIX_TYPE val);

  int  ML_fixCoarseMtx(
          ML_Operator *Cmat,          /*-up- coarse operator in MSR format   */
          const int    CoarseMtxType  /*-in- coarse-lev mtx storage type     */
  );
 
  int  ML_qr_fix_Bitsize(void);
}
#else

int ML_qr_fix_Create(const int nCoarseNod);

int ML_qr_fix_Destroy(void);

int ML_qr_fix_Print(ML_qr_fix* ptr);

int ML_qr_fix_NumDeadNodDof(void);

ML_QR_FIX_TYPE ML_qr_fix_getDeadNod(const int inx);

void ML_qr_fix_setNumDeadNod(int num);

void ML_qr_fix_setDeadNod( const int inx, ML_QR_FIX_TYPE val);

int  ML_fixCoarseMtx(
        ML_Operator *Cmat,          /*-up- coarse operator in MSR format     */
        const int    CoarseMtxType  /*-in- coarse-lev mtx storage type       */
     );

int  ML_qr_fix_Bitsize(void);

#endif
#endif

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_QR_FIX_H
#define ML_QR_FIX_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
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

  ML maintains a "dead DOF" data structure during each multigrid level's construction.  The data structure keeps
  the total number of aggregates, the total number of dead DOFs, and a std::vector with one entry per aggregate.
  This entry is itself a std::vector of bools.  Due to a quirk in the standard library, std::vector<bool> uses
  bits for storage.  So for example, given std::vector<bool> v(8), v is 8 bits, rather than 8 bytes.

  If nullspace vector i can be represented in the current aggregate, then entry i is set to false.  If vector i *cannot* be
  represented, then bit i is set to true.

  This encoding scheme is used in two spots.  The first is to zero out (local to an aggregate) NS vectors that
  cannot be represented because the aggregate is too small (ml_agg_uncoupled.c and ml_agg_MIS.c).  The second is to fix up
  the coarse grid matrix A so that it does not have a zero on the diagonal (ml_struct.c).

*/

#ifdef __cplusplus
extern "C" {
#endif

   int ML_qr_fix_Create(const int nCoarseNod, const int nullspaceDim);

   int ML_qr_fix_Destroy(void);

   int ML_qr_fix_NumDeadNodDof(void);

   void ML_qr_fix_setNumDeadNod(int num);

   int  ML_fixCoarseMtx(
          ML_Operator *Cmat,          /*-up- coarse operator in MSR format   */
          const int    CoarseMtxType  /*-in- coarse-lev mtx storage type     */
  );

  /* Mark all appropriate coarse DOFs as dead. */
   void ML_qr_fix_markDOFsAsDead(const int aggNum, const int aggSize, const int nullspaceDim);

  /* Returns status of coarse DOF. */
   int ML_qr_fix_isDOFDead(const int aggNum, const int coarseDOF);

  /* Returns 1 if aggregate (node) has dead coarse DOFs associated with it.  Returns 0 otherwise. */
   int ML_qr_fix_nodeHasDeadDOFs(const int aggNum);

#ifdef __cplusplus
}
#endif

#endif /*ifndef ML_QR_FIX_H*/

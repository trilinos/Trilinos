#ifndef EPETRA_CRSKUNDERTSPARSE_H
#define EPETRA_CRSKUNDERTSPARSE_H

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_MultiVector;

//! Epetra_CrsKundertSparse:  Solves an Epetra_LinearProblem using the solver "Sparse" from ChiliSpice.
/*! Epetra_CrsKundertSparse solves a fully-defined Epetra_LinearProblem (see the Epetra user documentation)
    using an optimized version of Ken Kundert's sparse solver called Sparse.  

    A fully-defined Epetra_LinearProblem is an Epetra_LinearProblem instance with the Epetra_RowMatrix and the RHS
    and LHS Epetra_MultiVectors defined and compatible with each other. If one or more of these attributes is not
    defined, indeterminant behavior may occur.

    Kundert's sparse solver was developed as part of the circuit simulation code called Spice.  Sandia National Labs
    has developed an optimized version of this code called ChiliSpice.  As part of the development of ChiliSpice,
    Kundert's solver was also optimized, with most of the work being done by Dave Shirley.  Questions about the C
    code in this directory are probably best answered by Dr. Shirley (dnshirl@sandia.gov).

    Epetra_CrsKundertSparse is a wrapper around Sparse that takes an Epetra_LinearProblem and solves it.
    This class is designed to support repeated calls to the Solve() method where each time Solve() is called the 
    matrix has \e identical structure to the initial matrix that was passed in to the Epetra_CrsKundertSparse
    constructor.  The first call to Solve() tends to be very expensive relative to subsequent calls.  This is because
    solver data structures are being allocated and because the initial (Markowitz) reordering is being computed. 
    Occasionally, a subsequent call may be expensive because Sparse compares the diagonal values of the 
    incoming matrix with values from the matrix for which the reordering was computed.  If the incoming matrix values
    have changed "too much" then matrix reordering is computed again, and again the overhead cost is paid.

    Epetra_CrsKundertSparse  provides access to Sparse because Sparse has been proven a valuable solver for small
    matrices (order 100-5000).  Initial results show that SuperLU may be more effective than Sparse for larger 
    problems.

    Epetra_CrsKundertSparse is a serial-only solver.  A shared memory parallel version of Sparse exists, but we
    do not support it at this time.  It is possible to use Epetra_CrsKundertSparse for distributed memory linear
    problems by first collapsing the problem to one processor (easy to do with an Epetra_Import object).  We will
    add this capability if needed in the future.
*/
class Epetra_CrsKundertSparse {

 public:

  //@{ \name Constructors/destructors.
  //! Construct an Epetra_CrsKundertSparse using a fully-defined Epetra_LinearProblem.
  /*! This constructor takes a fully-defined Epetra_LinearProblem and, optionally, values
      for the relative and absolute thresholds and diagonal pivots.  The optional parameters are defined
      below using an excerpt from the internal documentation for Sparse.

      \param InOut
             Problem - An Epetra_LinearProblem with a well-defined Epetra_RowMatrix, RHS and LHS.
      \param In
             RelThreshold - This number determines what the pivot relative threshold will
	     be.  It should be between zero and one.  If it is one then the
	     pivoting method becomes complete pivoting, which is very slow
	     and tends to fill up the matrix.  If it is set close to zero
	     the pivoting method becomes strict Markowitz with no
	     threshold.  The pivot threshold is used to eliminate pivot
	     candidates that would cause excessive element growth if they
	     were used.  Element growth is the cause of roundoff error.
	     Element growth occurs even in well-conditioned matrices.
	     Setting the RelThreshold large will reduce element growth and
	     roundoff error, but setting it too large will cause execution
	     time to be excessive and will result in a large number of
	     fill-ins.  If this occurs, accuracy can actually be degraded
	     because of the large number of operations required on the
	     matrix due to the large number of fill-ins.  A good value seems
	     to be 0.001.  The default is chosen by giving a value larger
	     than one or less than or equal to zero.  This value should be
	     increased and the matrix resolved if growth is found to be
	     excessive.  Changing the pivot threshold does not improve
	     performance on matrices where growth is low, as is often the
	     case with ill-conditioned matrices.  Once a valid threshold is
	     given, it becomes the new default.  The default value of
	     RelThreshold was choosen for use with nearly diagonally
	     dominant matrices such as node- and modified-node admittance
	     matrices.  For these matrices it is usually best to use
	     diagonal pivoting.  For matrices without a strong diagonal, it
	     is usually best to use a larger threshold, such as 0.01 or
	     0.1.
    \param In 
           AbsThreshold - The absolute magnitude an element must have to be considered
	   as a pivot candidate, except as a last resort.  This number
	   should be set significantly smaller than the smallest diagonal
	   element that is is expected to be placed in the matrix.  If
	   there is no reasonable prediction for the lower bound on these
	   elements, then AbsThreshold should be set to zero.
	   AbsThreshold is used to reduce the possibility of choosing as a
	   pivot an element that has suffered heavy cancellation and as a
	   result mainly consists of roundoff error.  Once a valid
	   threshold is given, it becomes the new default.
    \param In
           DiagPivoting - A flag indicating that pivot selection should be confined to the
	   diagonal if possible.  If DiagPivoting is nonzero and if
	   DIAGONAL_PIVOTING is enabled pivots will be chosen only from
	   the diagonal unless there are no diagonal elements that satisfy
	   the threshold criteria.  Otherwise, the entire reduced
	   submatrix is searched when looking for a pivot.  The diagonal
	   pivoting in Sparse is efficient and well refined, while the
	   off-diagonal pivoting is not.  For symmetric and near symmetric
	   matrices, it is best to use diagonal pivoting because it
	   results in the best performance when reordering the matrix and
	   when factoring the matrix without ordering.  If there is a
	   considerable amount of nonsymmetry in the matrix, then
	   off-diagonal pivoting may result in a better equation ordering
	   simply because there are more pivot candidates to choose from.
	   A better ordering results in faster subsequent factorizations.
	   However, the initial pivot selection process takes considerably
	   longer for off-diagonal pivoting.
  */
  Epetra_CrsKundertSparse( Epetra_LinearProblem * Problem,
			   const double RelThreshold = 0.001,
			   const double AbsThreshold = 1.0E-13,
			   const int DiagPivoting = 1);

  ~Epetra_CrsKundertSparse();
  //@}

  //@{ \name Solve methods.
  //! Solves the current Epetra_LinearProblem using the Sparse solver from ChiliSpice.
  /*! Given a well-defined Epetra_LinearProblem, this method will return the solution in the
      LHS multivector (currently only works with a multivector containing one vector).
      Typically Solve() is called many times, each time with a new Epetra_LinearProblem where the
      pattern of the matrix is identical but the values have changed.  The first call to Solve() is
      typically very expensive relative to subsequent calls, because a one-time reordering is performed.
  */
  int Solve();
  //@}

 private:

  void deleteArrays();
  double RelThreshold_;
  double AbsThreshold_;
  int DiagPivoting_;

  Epetra_CrsMatrix * A_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  int NumMyRows_;
  int NumMyCols_;
  int NumGlobalRows_;
  int NumGlobalCols_;


  char *Matrix_;
  double **addr_list_;

  bool FirstSolve_;
};

#endif // EPETRA_CRSKUNDERTSPARSE_H

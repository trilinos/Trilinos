#ifndef _PETRA_RDP_LINEARPROBLEM_H_
#define _PETRA_RDP_LINEARPROBLEM_H_

//! Petra_RDP_LinearProblem:  The Petra Linear Problem Class.
/*! The Petra_RDP_LinearProblem class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
  Currently it accepts a Petra matrix, initial guess and RHS and
  returns the solution.
  the elapsed time for each calling processor.
*/

#include "Petra_Petra.h"
#include "Petra_Comm.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_RowMatrix.h"
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#ifdef PETRA_MPI
#include "mpi.h"
#endif
enum ProblemDifficultyLevel {easy, moderate, hard, unsure};
#endif

class Petra_RDP_LinearProblem {
    
  public:
  //!  Petra_RDP_LinearProblem Default Constructor.
  /*! Creates an empty Petra_RDP_LinearProblem instance. The operator A, left-hand-side X
      and right-hand-side B must be set use the SetOperator(), SetLHS() and SetRHS()
      methods respectively.
  */
  Petra_RDP_LinearProblem(void);

  //!  Petra_RDP_LinearProblem Constructor.
  /*! Creates a Petra_RDP_LinearProblem instance. 
  */
  Petra_RDP_LinearProblem(Petra_RDP_RowMatrix * A, Petra_RDP_MultiVector * X,
			 Petra_RDP_MultiVector * B);

  //! Petra_RDP_LinearProblem Copy Constructor.
  /*! Makes copy of an existing Petra_RDP_LinearProblem instance.
  */
  Petra_RDP_LinearProblem(const Petra_RDP_LinearProblem& Problem);

  // Solver strategy assertions

  void AssertSymmetric(){OperatorSymmetric_ = true;};
#ifdef DOXYGEN_SHOULD_SKIP_THIS
  enum ProblemDifficultyLevel {easy, moderate, hard, unsure};
#endif
  //! Set problem difficulty level.
  /*! Sets Aztec options and parameters based on a definition of easy moderate or hard problem.
      Relieves the user from explicitly setting a large number of individual parameter values.
      This function can be used in conjunction with the SetOptions() and SetParams() functions.
  */
  void SetPDL(ProblemDifficultyLevel PDL) {PDL_ = PDL;};

  //! Set Operator A of linear problem AX = B.
  /*! Sets a pointer to a Petra_RDP_RowMatrix.  No copy of the operator is made.
  */
  void SetOperator(Petra_RDP_RowMatrix * A) {A_ = A;};

  //! Set left-hand-side X of linear problem AX = B.
  /*! Sets a pointer to a Petra_RDP_MultiVector.  No copy of the object is made.
  */
  void SetLHS(Petra_RDP_MultiVector * X) {X_ = X;};

  //! Set right-hand-side B of linear problem AX = B.
  /*! Sets a pointer to a Petra_RDP_MultiVector.  No copy of the object is made.
  */
  void SetRHS(Petra_RDP_MultiVector * B) {B_ = B;};


  //! Check input parameters for size consistency.
  /*! Returns 0 if all input parameters are valid */
  int CheckInput() const;

  //! Perform left scaling of a linear problem.
  /*! Applies the scaling vector D to the left side of the matrix A() and
    to the right hand side B().
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the ith row of A() and B().
    \return Integer error code, set to 0 if successful.
  */
  int LeftScale(const Petra_RDP_Vector & D);

  //! Perform right scaling of a linear problem.
  /*! Applies the scaling vector D to the right side of the matrix A().
      Apply the inverse of D to the initial guess.
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the ith row of A().  1/D[i] will be applied to the
               ith row of B().
    \return Integer error code, set to 0 if successful.
  */
  int RightScale(const Petra_RDP_Vector & D);

  //! Petra_RDP_LinearProblem Destructor.
  /*! Completely deletes a Petra_RDP_LinearProblem object.  
  */
  virtual ~Petra_RDP_LinearProblem(void);

  //@{ \name Accessor methods
  //! Get a pointer to the operator A.
  Petra_RDP_RowMatrix * GetOperator() const {return(A_);};
  //! Get a pointer to the left-hand-side X.
  Petra_RDP_MultiVector * GetLHS() const {return(X_);};
  //! Get a pointer to the right-hand-side B.
  Petra_RDP_MultiVector * GetRHS() const {return(B_);};
  //! Get problem difficulty level.
  ProblemDifficultyLevel GetPDL() const {return(PDL_);};
  //! Get operator symmetry bool.
  bool IsOperatorSymmetric() const {return(OperatorSymmetric_);};
  //@}

 private:

  Petra_RDP_RowMatrix * A_;
  Petra_RDP_MultiVector * X_;
  Petra_RDP_MultiVector * B_;

  bool OperatorSymmetric_;
  ProblemDifficultyLevel PDL_;
  bool LeftScaled_;
  bool RightScaled_;
  Petra_RDP_Vector * LeftScaleVector_;
  Petra_RDP_Vector * RightScaleVector_;
    
};

#endif /* _PETRA_RDP_LINEARPROBLEM_H_ */


/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_LINEARPROBLEM_H_
#define _EPETRA_LINEARPROBLEM_H_

#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#ifndef DOXYGEN_SHOULD_SKIP_THIS
enum ProblemDifficultyLevel {easy, moderate, hard, unsure};
#endif

//! Epetra_LinearProblem:  The Epetra Linear Problem Class.
/*! The Epetra_LinearProblem class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
  Currently it accepts a Epetra matrix, initial guess and RHS and
  returns the solution.
  the elapsed time for each calling processor.
*/


class Epetra_LinearProblem {
    
  public:
  //!  Epetra_LinearProblem Default Constructor.
  /*! Creates an empty Epetra_LinearProblem instance. The operator A, left-hand-side X
      and right-hand-side B must be set use the SetOperator(), SetLHS() and SetRHS()
      methods respectively.
  */
  Epetra_LinearProblem(void);

  //!  Epetra_LinearProblem Constructor.
  /*! Creates a Epetra_LinearProblem instance. 
  */
  Epetra_LinearProblem(Epetra_RowMatrix * A, Epetra_MultiVector * X,
			 Epetra_MultiVector * B);

  //! Epetra_LinearProblem Copy Constructor.
  /*! Makes copy of an existing Epetra_LinearProblem instance.
  */
  Epetra_LinearProblem(const Epetra_LinearProblem& Problem);

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
  /*! Sets a pointer to a Epetra_RowMatrix.  No copy of the operator is made.
  */
  void SetOperator(Epetra_RowMatrix * A) {A_ = A;};

  //! Set left-hand-side X of linear problem AX = B.
  /*! Sets a pointer to a Epetra_MultiVector.  No copy of the object is made.
  */
  void SetLHS(Epetra_MultiVector * X) {X_ = X;};

  //! Set right-hand-side B of linear problem AX = B.
  /*! Sets a pointer to a Epetra_MultiVector.  No copy of the object is made.
  */
  void SetRHS(Epetra_MultiVector * B) {B_ = B;};


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
  int LeftScale(const Epetra_Vector & D);

  //! Perform right scaling of a linear problem.
  /*! Applies the scaling vector D to the right side of the matrix A().
      Apply the inverse of D to the initial guess.
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the ith row of A().  1/D[i] will be applied to the
               ith row of B().
    \return Integer error code, set to 0 if successful.
  */
  int RightScale(const Epetra_Vector & D);

  //! Epetra_LinearProblem Destructor.
  /*! Completely deletes a Epetra_LinearProblem object.  
  */
  virtual ~Epetra_LinearProblem(void);

  //@{ \name Accessor methods
  //! Get a pointer to the operator A.
  Epetra_RowMatrix * GetOperator() const {return(A_);};
  //! Get a pointer to the left-hand-side X.
  Epetra_MultiVector * GetLHS() const {return(X_);};
  //! Get a pointer to the right-hand-side B.
  Epetra_MultiVector * GetRHS() const {return(B_);};
  //! Get problem difficulty level.
  ProblemDifficultyLevel GetPDL() const {return(PDL_);};
  //! Get operator symmetry bool.
  bool IsOperatorSymmetric() const {return(OperatorSymmetric_);};
  //@}

 private:

  Epetra_RowMatrix * A_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  bool OperatorSymmetric_;
  ProblemDifficultyLevel PDL_;
  bool LeftScaled_;
  bool RightScaled_;
  Epetra_Vector * LeftScaleVector_;
  Epetra_Vector * RightScaleVector_;
    
};

#endif /* _EPETRA_LINEARPROBLEM_H_ */

#ifndef _PETRA_BLAS_EASY_DLU_PROBLEM_H_
#define _PETRA_BLAS_EASY_DLU_PROBLEM_H_

//! Petra_BLAS_Easy_DLU_Problem: A class for solving "straight-forward" linear problems.

/*! The Petra_BLAS_Easy_DLU_Problem class supports the construction and solution of straight-forward linear problems of the form:
    \f[ A X = B, \f]
    where \f$A\f$ and \f$B\f$ are known matrices and \f$X\f$ is unknown.
    The Petra_BLAS_Easy_DLU_Problem class is intended to provide basic support for solving linear 
    problems for general dense rectangular (or square) matrices.  It is written on top of BLAS and LAPACK and thus has excellent
    performance and numerical capabilities.  Using this class, one can perform simple factorizations and solves.  To solve more 
    difficult problems, one should use the Petra_BLAS_Hard_DLU_Problem or Petra_BLAS_DSV_Problem classes.

<b>Petra_BLAS_Easy_DLU_Problem vs. Petra_LAPACK</b>

The Petra_LAPACK class provides access to most of the same functionality as Petra_BLAS_Easy_DLU_Problem.
The primary difference is that Petra_LAPACK is a "thin" layer on top of LAPACK and Petra_BLAS_Easy_DLU_Problem
attempts to provide object oriented access to solving dense linear systems.
<ul>
<li> When you should use Petra_LAPACK:  If you are simply looking for a convenient wrapper around the Fortran LAPACK
     routines and you have a well-conditioned problem, you should probably use Petra_LAPACK directly.
<li> When you should use Petra_BLAS_Easy_DLU_Problem: If you want to (or potentially want to) 
     work with a more object-oriented interface, you should probably use Petra_BLAS_Easy_DLU_Problem.
     
</ul>

<b>Constructing Petra_BLAS_Easy_DLU_Problem Objects</b>

There is a single (default) constructor for this class, taking no arguments.  In order to use methods in this class, one
must use the SetOperator() method to set the matrix. The argument A must be a Petra_BLAS_DGE_Matrix.  Some methods, e.g., Factor(), 
require only A to be set.  Others require that X and B also be set.

<b>Setting vectors used for linear solves</b>

Setting the X and B vectors (which are Petra_BLAS_DGE_Matrix objects) used for solving linear systems is 
done separately from the constructor.  
This allows a single matrix factor to be used for multiple solves.  Similar to the constructor, 
the vectors X and B can be copied or viewed using the Petra_DataAccess argument.

<b>Extracting Data from Petra_BLAS_Easy_DLU_Problem Objects</b>

Once a Petra_BLAS_Easy_DLU_Problem is constructed, it is possible to view the data via access functions.

\warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.


<b>Vector and Utility Functions</b>

Once a Petra_BLAS_Easy_DLU_Problem is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Factorizations.
  <li> Solves.
  <li> Condition estimates.
  <li> Norms.
</ul>

The final useful function is Flops().  Each Petra_BLAS_Easy_DLU_Problem object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Petra_Time class, one can get accurate parallel performance
numbers.

<b>Strategies for Solving Linear Systems</b>
In many cases, linear systems can be accurately solved by simply computing the LU factorization
of the matrix and then performing a forward back solve with a given set of right hand side vectors.  The 
Petra_BLAS_Easy_DLU_Problem class is intended for these problems. However,
in some instances, the factorization may be very poorly conditioned and this simple approach may not work.  In
these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
the factorization. Petra_BLAS_Hard_DLU_Problem will use equilibration with the factorization if, once the object
is constructed and \e before it is factored, you call the function FactorWithEquilibration(true) to force 
equilibration to be used.  If you are uncertain if equilibration should be used, you may call the function
ShouldEquilibrate() which will return true if equilibration could possibly help.  ShouldEquilibrate() uses
guidelines specified in the LAPACK User Guide, namely if SCOND < 0.1 and AMAX < Underflow or AMAX > Overflow, to 
determine if equilibration \e might be useful. 
 
Petra_BLAS_Hard_DLU_Problem will use iterative refinement after a forward/back solve if you call
SolveToRefinedSolution(true).  It will also compute forward and backward error estimates if you call
EstimateSolutionErrors(true).  Access to the forward (back) error estimates is available via FERR() (BERR()).

Examples using Petra_BLAS_Easy_DLU_Problem and Petra_BLAS_Hard_DLU_Problem can be found in the Petra test directories.

*/
#include "Petra_Petra.h" 
#include "Petra_Flops.h"
#include "Petra_BLAS.h"
#include "Petra_LAPACK.h"
#include "Petra_BLAS_DGE_Matrix.h"


//=========================================================================
class Petra_BLAS_Easy_DLU_Problem : public Petra_Flops, public Petra_BLAS, public Petra_LAPACK{

  public:
  
  //! Default constructor; defines an empty object.
  /*!
    Petra_BLAS_Easy_DLU_Problem objects defined by the default constructor should be completed with calls to
    SetOperator(), SetLHS() and SetRHS().  
   */
  Petra_BLAS_Easy_DLU_Problem(void);
    
  //! Petra_BLAS_Easy_DLU_Problem copy constructor.
  
  Petra_BLAS_Easy_DLU_Problem(const Petra_BLAS_Easy_DLU_Problem& Source);

  //! Petra_BLAS_Easy_DLU_Problem destructor.  
  virtual ~Petra_BLAS_Easy_DLU_Problem ();

  //! Set the matrix. 
  /* A should contain the matrix that will be factored.  The original values with be over-written (in-place factorization).
     \param A In 
            Matrix as a Petra_BLAS_DGE_Matrix object.
    \return Integer error code, set to 0 if successful.
   */
  int SetOperator(Petra_BLAS_DGE_Matrix & A);

  //! Set the right-hand-side(s). 
  /* B should contain the RHS vectors.  It is possible to have the
     solution returned in B by pasing in the same object for both B and X.
     \param B In 
            Right hand side as a Petra_BLAS_DGE_Matrix object.
    \return Integer error code, set to 0 if successful.
   */
  int SetRHS(const Petra_BLAS_DGE_Matrix & B);

  //! Set the target for the solution vectors 
  /* X should be space allocated for the solution.  It is possible to have the
     solution returned in B by pasing in the same object for both B and X.
     \param X In 
            Left hand side as a Petra_BLAS_DGE_Matrix object.

    \return Integer error code, set to 0 if successful.
   */
  int SetLHS(Petra_BLAS_DGE_Matrix & X);

  //! If Flag is true, causes all subsequent function calls to work with the transpose of \e this matrix, otherwise not.
  void SolveWithTranspose(bool Flag) {Transpose_ = Flag; return;};

  //! Computes the in-place LU factorization of the matrix using the LAPACK routine \e DGETRF. SetOperator() must have been called.
  /*!
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Factor(void);

  //! Computes the solution X to AX = B. SetOperator(), SetLHS() and SetRHS() must have been called first.
  /*!
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Solve(void);

  //! Inverts the \e this matrix.
  /*! Note: This function works a little differently that DPOTRI in that it fills the entire
      matrix with the inverse, independent of the UPLO specification.

    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Invert(void);

  //! Returns the reciprocal of the 1-norm condition number of the \e this matrix.
  /*! 
    \param Value Out
           On return contains the reciprocal of the 1-norm condition number of the \e this matrix.
    
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int ReciprocalConditionEstimate(double & Value);

  //! Returns true if transpose of \e this matrix has and will be used.
  bool Transpose() {return(Transpose_);};

  //! Returns true if matrix is factored (factor available via Operator()).
  bool Factored() {return(Factored_);};

  //! Returns true if matrix inverse has been computed (inverse available via Operator()).
  bool Inverted() {return(Inverted_);};

  //! Returns true if reciprocal of matrix condition number has been computed.
  bool ReciprocalConditionEstimated() {return(ReciprocalConditionEstimated_);};

  //! Returns true if the current set of vectors has been solved.
  bool Solved() {return(Solved_);};

  //! Returns pointer to the \e this matrix.
  Petra_BLAS_DGE_Matrix & Operator()  const {return(*Operator_);};

  //! Returns pointer to current RHS.
  Petra_BLAS_DGE_Matrix & RHS()  const {return(*RHS_);};

  //! Returns pointer to current solution.
  Petra_BLAS_DGE_Matrix & LHS()  const {return(*LHS_);};

  //! Returns pointer to pivot vector (if factorization has been computed), zero otherwise.
  int * IPIV()  const {return(IPIV_);};
  
  void Print (ostream& os) const;

 protected:

  void AllocateWORK(void) {if (WORK_==0) {LWORK_ = 4*Operator().N(); WORK_ = new double[LWORK_];} return;};
  void AllocateIWORK(void) {if (IWORK_==0) IWORK_ = new int[Operator().N()]; return;};
  void DeleteArrays(void);
  bool Transpose_;
  bool Factored_;
  bool Solved_;
  bool Inverted_;
  bool ReciprocalConditionEstimated_;

  int INFO_;
  int LWORK_;

  int * IPIV_;
  int * IWORK_;
  double * WORK_;

  double RCOND_;


  Petra_BLAS_DGE_Matrix * Operator_;

  Petra_BLAS_DGE_Matrix * RHS_;
  Petra_BLAS_DGE_Matrix * LHS_;


};

#endif /* _PETRA_BLAS_EASY_DLU_PROBLEM_H_ */

/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef IFPACK_SERIALTRIDISOLVER_H
#define IFPACK_SERIALTRIDISOLVER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif
class Ifpack_SerialTriDiMatrix;
//class Epetra_SerialDenseVector;
class Epetra_SerialDenseMatrix;

#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_BLAS.h"
#include "Teuchos_LAPACK.hpp"


//! Ifpack_SerialTriDiSolver: A class for solving TriDi linear problems.

/*! The Ifpack_SerialTriDiSolver class enables the definition, in terms of Ifpack_SerialTriDiMatrix
    and Ifpack_SerialTriDiVector objects, of a TriDi linear problem, followed by the solution of that problem via the
    most sophisticated techniques available in LAPACK.

The Ifpack_SerialTriDiSolver class is intended to provide full-featured support for solving linear
problems for general TriDi rectangular (or square) matrices.  It is written on top of BLAS and LAPACK and thus has excellent
performance and numerical capabilities.  Using this class, one can either perform simple factorizations and solves or
apply all the tricks available in LAPACK to get the best possible solution for very ill-conditioned problems.

<b>Ifpack_SerialTriDiSolver vs. Epetra_LAPACK</b>

The Epetra_LAPACK class provides access to most of the same functionality as Ifpack_SerialTriDiSolver.
The primary difference is that Epetra_LAPACK is a "thin" layer on top of LAPACK and Ifpack_SerialTriDiSolver
attempts to provide easy access to the more sophisticated aspects of solving TriDi linear and eigensystems.
<ul>
<li> When you should use Epetra_LAPACK:  If you are simply looking for a convenient wrapper around the Fortran LAPACK
     routines and you have a well-conditioned problem, you should probably use Epetra_LAPACK directly.
<li> When you should use Ifpack_SerialTriDiSolver: If you want to (or potentially want to) solve ill-conditioned
     problems or want to work with a more object-oriented interface, you should probably use Ifpack_SerialTriDiSolver.

</ul>

<b>Constructing Ifpack_SerialTriDiSolver Objects</b>

There is a single Ifpack_SerialTriDiSolver constructor.   However, the matrix, right hand side and solution
vectors must be set prior to executing most methods in this class.

<b>Setting vectors used for linear solves</b>

The matrix A, the left hand side X and the right hand side B (when solving AX = B, for X), can be set by appropriate set
methods.  Each of these three objects must be an Ifpack_SerialTriDiMatrix or and Epetra_SerialDenseVector object.  The
set methods are as follows:
<ul>
<li> SetMatrix() - Sets the matrix.
<li> SetVectors() - Sets the left and right hand side vector(s).
</ul>

<b>Vector and Utility Functions</b>

Once a Ifpack_SerialTriDiSolver is constructed, several mathematical functions can be applied to
the object.  Specifically:
<ul>
  <li> Factorizations.
  <li> Solves.
  <li> Condition estimates.
  <li> Equilibration.
  <li> Norms.
</ul>

<b>Counting floating point operations </b>
The Ifpack_SerialTriDiSolver class has Epetra_CompObject as a base class.  Thus, floating point operations
are counted and accumulated in the Epetra_Flop object (if any) that was set using the SetFlopCounter()
method in the Epetra_CompObject base class.

<b>Strategies for Solving Linear Systems</b>
In many cases, linear systems can be accurately solved by simply computing the LU factorization
of the matrix and then performing a forward back solve with a given set of right hand side vectors.  However,
in some instances, the factorization may be very poorly conditioned and this simple approach may not work.  In
these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
the factorization.

LAPACK does not provide an equilibration functionality for tridiagonal data structures.

Ifpack_SerialTriDiSolver will use iterative refinement after a forward/back solve if you call
SolveToRefinedSolution(true).  It will also compute forward and backward error estimates if you call
EstimateSolutionErrors(true).  Access to the forward (back) error estimates is available via FERR() (BERR()).

Examples using Ifpack_SerialTriDiSolver can be found in the Epetra test directories.
  
*/

//=========================================================================
class Ifpack_SerialTriDiSolver :
  public Epetra_CompObject, public Epetra_BLAS,
  public Epetra_Object    {
  public:
  
    //! @name Constructor/Destructor Methods
  //@{
  //! Default constructor; matrix should be set using SetMatrix(), LHS and RHS set with SetVectors().
  Ifpack_SerialTriDiSolver();


  //! Ifpack_SerialTriDiSolver destructor.
  virtual ~Ifpack_SerialTriDiSolver();
  //@}

  //! @name Set Methods
  //@{

  //! Sets the pointers for coefficient matrix
  int SetMatrix(Ifpack_SerialTriDiMatrix & A);

  //! Sets the pointers for left and right hand side vector(s).
  /*! Row dimension of X must match column dimension of matrix A, row dimension of B
      must match row dimension of A.  X and B must have the same dimensions.
  */
  int SetVectors(Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & B);
  //@}


  //! @name Strategy modifying Methods
  //@{

  //! Causes equilibration to be called just before the matrix factorization as part of the call to Factor.
  /*! This function must be called before the factorization is performed.
   */
  //  void FactorWithEquilibration(bool Flag) {Equilibrate_ = Flag; return;};

  //! If Flag is true, causes all subsequent function calls to work with the transpose of \e this matrix, otherwise not.
  void SolveWithTranspose(bool Flag) {Transpose_ = Flag; if (Flag) TRANS_ = 'T'; else TRANS_ = 'N'; return;};

  //! Causes all solves to compute solution to best ability using iterative refinement.
  void SolveToRefinedSolution(bool Flag) {RefineSolution_ = Flag; return;};

  //! Causes all solves to estimate the forward and backward solution error.
  /*! Error estimates will be in the arrays FERR and BERR, resp, after the solve step is complete.
      These arrays are accessible via the FERR() and BERR() access functions.
  */
  void EstimateSolutionErrors(bool Flag) ;
  //@}

  //! @name Factor/Solve/Invert Methods
  //@{

  //! Computes the in-place LU factorization of the matrix using the LAPACK routine \e DGETRF.
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual int Factor(void);

  //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual int Solve(void);

  //! Inverts the \e this matrix.
  /*!
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int Invert(void);

  //! Apply Iterative Refinement.
  /*!
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int ApplyRefinement(void);

  //! Unscales the solution vectors if equilibration was used to solve the system.
  /*!
    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  // int UnequilibrateLHS(void);

  //! Returns the reciprocal of the 1-norm condition number of the \e this matrix.
  /*!
    \param Value Out
           On return contains the reciprocal of the 1-norm condition number of the \e this matrix.

    \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
  */
  virtual int ReciprocalConditionEstimate(double & Value);
  //@}

  //! @name Query methods
  //@{

  //! Returns true if transpose of \e this matrix has and will be used.
  bool Transpose() {return(Transpose_);};

  //! Returns true if matrix is factored (factor available via AF() and LDAF()).
  bool Factored() {return(Factored_);};

  //! Returns true if forward and backward error estimated have been computed (available via FERR() and BERR()).
  bool SolutionErrorsEstimated() {return(SolutionErrorsEstimated_);};

  //! Returns true if matrix inverse has been computed (inverse available via AF() and LDAF()).
  bool Inverted() {return(Inverted_);};

  //! Returns true if the condition number of the \e this matrix has been computed (value available via ReciprocalConditionEstimate()).
  bool ReciprocalConditionEstimated() {return(ReciprocalConditionEstimated_);};

  //! Returns true if the current set of vectors has been solved.
  bool Solved() {return(Solved_);};

  //! Returns true if the current set of vectors has been refined.
  bool SolutionRefined() {return(SolutionRefined_);};
  //@}

  //! @name Data Accessor methods
  //@{

  //! Returns pointer to current matrix.
   Ifpack_SerialTriDiMatrix * Matrix()  const {return(Matrix_);};

  //! Returns pointer to factored matrix (assuming factorization has been performed).
   Ifpack_SerialTriDiMatrix * FactoredMatrix()  const {return(Factor_);};

  //! Returns pointer to current LHS.
  Epetra_SerialDenseMatrix * LHS()  const {return(LHS_);};

  //! Returns pointer to current RHS.
  Epetra_SerialDenseMatrix * RHS()  const {return(RHS_);};

  //! Returns column dimension of system.
  int N()  const {return(N_);};

  //! Returns pointer to the \e this matrix.
  double * A()  const {return(A_);};

  //! Returns the leading dimension of the \e this matrix.
  int LDA()  const {return(LDA_);};

  //! Returns pointer to current RHS.
  double * B()  const {return(B_);};

  //! Returns the leading dimension of the RHS.
  int LDB()  const {return(LDB_);};

  //! Returns the number of current right hand sides and solution vectors.
  int NRHS()  const {return(NRHS_);};

  //! Returns pointer to current solution.
  double * X()  const {return(X_);};

  //! Returns the leading dimension of the solution.
  int LDX()  const {return(LDX_);};

  //! Returns pointer to the factored matrix (may be the same as A() if factorization done in place).
  double * AF()  const {return(AF_);};

  //! Returns the leading dimension of the factored matrix.
  int LDAF()  const {return(LDAF_);};

  //! Returns pointer to pivot vector (if factorization has been computed), zero otherwise.
  int * IPIV()  const {return(IPIV_);};

  //! Returns the 1-Norm of the \e this matrix (returns -1 if not yet computed).
  double ANORM()  const {return(ANORM_);};

  //! Returns the reciprocal of the condition number of the \e this matrix (returns -1 if not yet computed).
  double RCOND()  const {return(RCOND_);};

  //! Ratio of smallest to largest row scale factors for the \e this matrix (returns -1 if not yet computed).
  /*! If ROWCND() is >= 0.1 and AMAX() is not close to overflow or underflow, then equilibration is not needed.
   */
  double ROWCND()  const {return(ROWCND_);};

  //! Ratio of smallest to largest column scale factors for the \e this matrix (returns -1 if not yet computed).
  /*! If COLCND() is >= 0.1 then equilibration is not needed.
   */
  double COLCND()  const {return(COLCND_);};

  //! Returns the absolute value of the largest entry of the \e this matrix (returns -1 if not yet computed).
  double AMAX()  const {return(AMAX_);};

  //! Returns a pointer to the forward error estimates computed by LAPACK.
  double * FERR()  const {return(FERR_);};

  //! Returns a pointer to the backward error estimates computed by LAPACK.
  double * BERR()  const {return(BERR_);};

  //@}

  //! @name I/O methods
  //@{
  //! Print service methods; defines behavior of ostream << operator.
  virtual void Print(std::ostream& os) const;
  //@}
 protected:

  void AllocateWORK() {if (WORK_==0) {LWORK_ = 4*N_; WORK_ = new double[LWORK_];} return;};
  void AllocateIWORK() {if (IWORK_==0) IWORK_ = new int[N_]; return;};
  void InitPointers();
  void DeleteArrays();
  void ResetMatrix();
  void ResetVectors();

  bool Transpose_;
  bool Factored_;
  bool EstimateSolutionErrors_;
  bool SolutionErrorsEstimated_;
  bool Solved_;
  bool Inverted_;
  bool ReciprocalConditionEstimated_;
  bool RefineSolution_;
  bool SolutionRefined_;

  char TRANS_;

  int N_;
  int Min_MN_;
  int NRHS_;
  int LDA_;
  int LDAF_;
  int LDB_;
  int LDX_;
  int INFO_;
  int LWORK_;

  int * IPIV_;
  int * IWORK_;

  double ANORM_;
  double RCOND_;
  double ROWCND_;
  double COLCND_;
  double AMAX_;

  Ifpack_SerialTriDiMatrix * Matrix_;
  Epetra_SerialDenseMatrix * LHS_;
  Epetra_SerialDenseMatrix * RHS_;
  Ifpack_SerialTriDiMatrix * Factor_;

  double * A_;
  double * FERR_;
  double * BERR_;
  double * AF_;
  double * WORK_;

  double * B_;
  double * X_;


 private:
  Teuchos::LAPACK<int,double> lapack;
  // Ifpack_SerialTriDiSolver copy constructor (put here because we don't want user access)

  Ifpack_SerialTriDiSolver(const Ifpack_SerialTriDiSolver& Source);
  Ifpack_SerialTriDiSolver & operator=(const Ifpack_SerialTriDiSolver& Source);
};

#endif /* EPETRA_SERIALTRIDISOLVER_H */

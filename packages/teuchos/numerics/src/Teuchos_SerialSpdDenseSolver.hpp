
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_SERIALSPDDENSESOLVER_H
#define TEUCHOS_SERIALSPDDENSESOLVER_H

/*! \file Teuchos_SerialSpdDenseSolver.hpp
  \brief Templated class for constructing and using Hermitian positive definite dense matrices.
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_SerialSymDenseMatrix.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

/*! \class Teuchos::SerialSpdDenseSolver
    \brief A class for constructing and using Hermitian positive definite dense matrices.

    The Teuchos::SerialSpdDenseSolver class enables the construction and use of a templated, 
    Hermitian positive definite, double-precision dense matrices.  It is built on the BLAS and 
    LAPACK via the Teuchos::BLAS and Teuchos::LAPACK classes. 

    The Teuchos::SerialSpdDenseSolver class is intended to provide full-featured support for 
    solving linear system problems for Hermitian positive definite matrices.  It is written on 
    top of BLAS and LAPACK and thus has excellent performance and numerical capabilities.  
    Using this class, one can either perform simple factorizations and solves or apply all the 
    tricks available in LAPACK to get the best possible solution for very ill-conditioned problems.

    <b>Teuchos::SerialSpdDenseSolver vs. Teuchos::LAPACK</b>

    The Teuchos::LAPACK class provides access to most of the same functionality as Teuchos::SerialSpdDenseSolver.
    The primary difference is that Teuchos::LAPACK is a "thin" layer on top of LAPACK and 
    Teuchos::SerialSpdDenseSolver attempts to provide easy access to the more sophisticated 
    aspects of solving dense linear and eigensystems.
    <ul>
    <li> When you should use Teuchos::LAPACK:  If you are simply looking for a convenient wrapper 
    around the Fortran LAPACK routines and you have a well-conditioned problem, you should probably 
    use Teuchos::LAPACK directly.
    <li> When you should use Teuchos::SerialSpdDenseSolver: If you want to (or potentially want to) 
    solve ill-conditioned problems or want to work with a more object-oriented interface, you should 
    probably use Teuchos::SerialSpdDenseSolver.
    </ul>

    <b>Extracting Data from Teuchos::SerialSpdDenseSolver Objects</b>
    
    Once a Teuchos::SerialSpdDenseSolver is constructed, it is possible to view the data via access functions.

    \warning Use of these access functions cam be \e extremely dangerous from a data hiding perspective.

    <b>Vector and Utility Functions</b>
    
    Once a Teuchos::SerialSpdDenseSolver is constructed, several mathematical functions can be applied to
    the object.  Specifically:
    <ul>
    <li> Factorizations.
    <li> Solves.
    <li> Condition estimates.
    <li> Equilibration.
    <li> Norms.
    </ul>
    
    <b>Strategies for Solving Linear Systems</b>
    In many cases, linear systems can be accurately solved by simply computing the Cholesky factorization
    of the matrix and then performing a forward back solve with a given set of right hand side vectors.  However,
    in some instances, the factorization may be very poorly conditioned and the simple approach may not work.  In
    these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
    the factorization. 
    
    Teuchos::SerialSpdDenseSolver will use equilibration with the factorization if, once the object
    is constructed and \e before it is factored, you call the function FactorWithEquilibration(true) to force 
    equilibration to be used.  If you are uncertain if equilibration should be used, you may call the function
    ShouldEquilibrate() which will return true if equilibration could possibly help.  ShouldEquilibrate() uses
    guidelines specified in the LAPACK User Guide, namely if SCOND < 0.1 and AMAX < Underflow or AMAX > Overflow, to 
    determine if equilibration \e might be useful. 
    
    Teuchos::SerialSpdDenseSolver will use iterative refinement after a forward/back solve if you call
    SolveToRefinedSolution(true).  It will also compute forward and backward error estimates if you call
    EstimateSolutionErrors(true).  Access to the forward (back) error estimates is available via FERR() (BERR()).
    
    Examples using Teuchos::SerialSpdDenseSolver can be found in the Teuchos test directories.  
*/

namespace Teuchos {

  template<typename OrdinalType, typename ScalarType>
  class SerialSpdDenseSolver : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>,
			    public LAPACK<OrdinalType, ScalarType>
  {
  public:
    
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    
    //! @name Constructor/Destructor Methods
    //@{ 
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    SerialSpdDenseSolver();
    
    
    //! SerialSpdDenseSolver destructor.  
    virtual ~SerialSpdDenseSolver();
    //@}
    
    //! @name Set Methods
    //@{ 
    
    //! Sets the pointers for coefficient matrix.
    int setMatrix(const RCP<SerialSymDenseMatrix<OrdinalType,ScalarType> >& A_in);
    
    //! Sets the pointers for left and right hand side vector(s).
    /*! Row dimension of X must match column dimension of matrix A, row dimension of B 
      must match row dimension of A.  X and B must have the same dimensions.
    */
    int setVectors(const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& X, 
		   const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& B);
    //@}
    
    //! @name Strategy Modifying Methods
    //@{ 
    
    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! This function must be called before the factorization is performed. 
     */
    void factorWithEquilibration(bool flag) {equilibrate_ = flag; return;}
    
    //! Causes all solves to compute solution to best ability using iterative refinement.
    void solveToRefinedSolution(bool flag) {refineSolution_ = flag; return;}
    
    //! Causes all solves to estimate the forward and backward solution error. 
    /*! Error estimates will be in the arrays FERR and BERR, resp, after the solve step is complete.
      These arrays are accessible via the FERR() and BERR() access functions.
    */
    void estimateSolutionErrors(bool flag);
    //@}
    
    //! @name Factor/Solve/Invert Methods
    //@{ 
    
    //! Computes the in-place Cholesky factorization of the matrix using the LAPACK routine \e DPOTRF.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor();
    
    //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve();
    
    //! Inverts the \e this matrix.
    /*! Note: This function works a little differently that DPOTRI in that it fills the entire
      matrix with the inverse, independent of the UPLO specification.
      
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int invert();
    
    //! Computes the scaling vector S(i) = 1/sqrt(A(i,i) of the \e this matrix.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int computeEquilibrateScaling();
    
    //! Equilibrates the \e this matrix.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int equilibrateMatrix();
    
    //! Equilibrates the current RHS.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int equilibrateRHS();
    
    
    //! Apply Iterative Refinement.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int applyRefinement();
    
    //! Unscales the solution vectors if equilibration was used to solve the system.
    /*! 
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int unequilibrateLHS();
    
    //! Returns the reciprocal of the 1-norm condition number of the \e this matrix.
    /*! 
      \param Value Out
      On return contains the reciprocal of the 1-norm condition number of the \e this matrix.
      
      \return Integer error code, set to 0 if successful. Otherwise returns the LAPACK error code INFO.
    */
    int reciprocalConditionEstimate(MagnitudeType & Value);
    //@}
    
    //! @name Query methods
    //@{ 
    
    //! Returns true if transpose of \e this matrix has and will be used.
    bool transpose() {return(transpose_);}
    
    //! Returns true if matrix is factored (factor available via AF() and LDAF()).
    bool factored() {return(factored_);}
    
    //! Returns true if factor is equilibrated (factor available via AF() and LDAF()).
    bool equilibratedA() {return(equilibratedA_);}
    
    //! Returns true if RHS is equilibrated (RHS available via B() and LDB()).
    bool equilibratedB() {return(equilibratedB_);}
    
    //! Returns true if the LAPACK general rules for equilibration suggest you should equilibrate the system.
    bool shouldEquilibrate() {computeEquilibrateScaling(); return(shouldEquilibrate_);}
    
    //! Returns true if forward and backward error estimated have been computed (available via FERR() and BERR()).
    bool solutionErrorsEstimated() {return(solutionErrorsEstimated_);}
    
    //! Returns true if matrix inverse has been computed (inverse available via AF() and LDAF()).
    bool inverted() {return(inverted_);}
    
    //! Returns true if the condition number of the \e this matrix has been computed (value available via ReciprocalConditionEstimate()).
    bool reciprocalConditionEstimated() {return(reciprocalConditionEstimated_);}
    
    //! Returns true if the current set of vectors has been solved.
    bool solved() {return(solved_);}
    
    //! Returns true if the current set of vectors has been refined.
    bool solutionRefined() {return(solutionRefined_);}
    //@}
    
    //! @name Data Accessor methods
    //@{ 
    
    //! Returns pointer to current matrix.
    RCP<SerialSymDenseMatrix<OrdinalType, ScalarType> > getMatrix()  const {return(Matrix_);}
    
    //! Returns pointer to factored matrix (assuming factorization has been performed).
    RCP<SerialSymDenseMatrix<OrdinalType, ScalarType> > getFactoredMatrix()  const {return(Factor_);}
    
    //! Returns pointer to current LHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getLHS() const {return(LHS_);}
    
    //! Returns pointer to current RHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getRHS() const {return(RHS_);}
    
    //! Returns row dimension of system.
    OrdinalType numRows()  const {return(numRowCols_);}
    
    //! Returns column dimension of system.
    OrdinalType numCols()  const {return(numRowCols_);}
    
    //! Returns the 1-Norm of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType ANORM()  const {return(ANORM_);}
    
    //! Returns the reciprocal of the condition number of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType RCOND()  const {return(RCOND_);}
    
    //! Ratio of smallest to largest equilibration scale factors for the \e this matrix (returns -1 if not yet computed).
    /*! If SCOND() is >= 0.1 and AMAX() is not close to overflow or underflow, then equilibration is not needed.
     */
    MagnitudeType SCOND() {return(SCOND_);};
    
    //! Returns the absolute value of the largest entry of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType AMAX()  const {return(AMAX_);}
    
    //! Returns a pointer to the forward error estimates computed by LAPACK.
    std::vector<ScalarType> FERR()  const {return(FERR_);}
    
    //! Returns a pointer to the backward error estimates computed by LAPACK.
    std::vector<ScalarType> BERR()  const {return(BERR_);}
    
    //! Returns a pointer to the row scaling vector used for equilibration.
    std::vector<ScalarType> R()  const {return(R_);}
    
    //@}

    //! @name I/O methods
    //@{
    //! Print service methods; defines behavior of ostream << operator.
    void Print(std::ostream& os) const;
    //@}
    
  protected:
    
    void allocateWORK() { LWORK_ = 4*numRowCols_; WORK_.resize( LWORK_ ); return;}
    void allocateIWORK() { IWORK_.resize( numRowCols_ ); return;}
    void resetMatrix();
    void resetVectors();
    
    bool equilibrate_;
    bool shouldEquilibrate_;
    bool equilibratedA_;
    bool equilibratedB_;
    bool transpose_;
    bool factored_;
    bool estimateSolutionErrors_;
    bool solutionErrorsEstimated_;
    bool solved_;
    bool inverted_;
    bool reciprocalConditionEstimated_;
    bool refineSolution_;
    bool solutionRefined_;
    
    OrdinalType numRowCols_;
    OrdinalType LDA_;
    OrdinalType LDAF_;
    OrdinalType INFO_;
    OrdinalType LWORK_;
    
    std::vector<int> IWORK_;
    
    MagnitudeType ANORM_;
    MagnitudeType RCOND_;
    MagnitudeType SCOND_;
    MagnitudeType AMAX_;
    
    RCP<SerialSymDenseMatrix<OrdinalType, ScalarType> > Matrix_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > LHS_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > RHS_;
    RCP<SerialSymDenseMatrix<OrdinalType, ScalarType> > Factor_;
    
    ScalarType * A_;
    ScalarType * AF_;
    std::vector<ScalarType> FERR_;
    std::vector<ScalarType> BERR_;
    std::vector<ScalarType> WORK_;
    std::vector<ScalarType> R_;
    
  private:
    // SerialSpdDenseSolver copy constructor (put here because we don't want user access)
    
    SerialSpdDenseSolver(const SerialSpdDenseSolver<OrdinalType, ScalarType>& Source);
    SerialSpdDenseSolver & operator=(const SerialSpdDenseSolver<OrdinalType, ScalarType>& Source);
    
  };  

//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialSpdDenseSolver<OrdinalType,ScalarType>::SerialSpdDenseSolver()
  : CompObject(),
    Object("Teuchos::SerialSpdDenseSolver"),
    equilibrate_(false),
    shouldEquilibrate_(false),
    equilibratedA_(false),
    equilibratedB_(false),
    transpose_(false),
    factored_(false),
    estimateSolutionErrors_(false),
    solutionErrorsEstimated_(false),
    solved_(false),
    inverted_(false),
    reciprocalConditionEstimated_(false),
    refineSolution_(false),
    solutionRefined_(false),
    numRowCols_(0),
    LDA_(0),
    LDAF_(0),
    INFO_(0),
    LWORK_(0),
    ANORM_(ScalarTraits<ScalarType>::zero()),
    RCOND_(ScalarTraits<ScalarType>::zero()),
    SCOND_(ScalarTraits<ScalarType>::zero()),
    AMAX_(ScalarTraits<ScalarType>::zero()),
    A_(0),
    AF_(0)
{
  resetMatrix();
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialSpdDenseSolver<OrdinalType,ScalarType>::~SerialSpdDenseSolver()
{}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialSpdDenseSolver<OrdinalType,ScalarType>::resetVectors()
{
  LHS_ = Teuchos::null;
  RHS_ = Teuchos::null;
  reciprocalConditionEstimated_ = false;
  solutionRefined_ = false;
  solved_ = false;
  solutionErrorsEstimated_ = false;
  equilibratedB_ = false;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialSpdDenseSolver<OrdinalType,ScalarType>::resetMatrix()
{
  resetVectors();
  equilibratedA_ = false;
  factored_ = false;
  inverted_ = false;
  numRowCols_ = 0;
  LDA_ = 0;
  LDAF_ = 0;
  ANORM_ = -ScalarTraits<MagnitudeType>::one();
  RCOND_ = -ScalarTraits<MagnitudeType>::one();
  SCOND_ = -ScalarTraits<MagnitudeType>::one();
  AMAX_ = -ScalarTraits<MagnitudeType>::one();
  A_ = 0;
  AF_ = 0;
  INFO_ = 0;
  LWORK_ = 0;
  R_.resize(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::setMatrix(const RCP<SerialSymDenseMatrix<OrdinalType,ScalarType> >& A) 
{
  resetMatrix();
  Matrix_ = A;
  Factor_ = A;
  numRowCols_ = A->numRows();
  LDA_ = A->stride();
  LDAF_ = LDA_;
  A_ = A->values();
  AF_ = A->values();
  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::setVectors(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& X, 
							     const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& B)
{
  // Check that these new vectors are consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(B->numRows()!=X->numRows() || B->numCols() != X->numCols(), std::invalid_argument,
		     "SerialSpdDenseSolver<T>::setVectors: X and B are not the same size!");
  TEUCHOS_TEST_FOR_EXCEPTION(B->values()==0, std::invalid_argument,
		     "SerialSpdDenseSolver<T>::setVectors: B is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(X->values()==0, std::invalid_argument,
		     "SerialSpdDenseSolver<T>::setVectors: X is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(B->stride()<1, std::invalid_argument,
		     "SerialSpdDenseSolver<T>::setVectors: B has an invalid stride!");
  TEUCHOS_TEST_FOR_EXCEPTION(X->stride()<1, std::invalid_argument,
		     "SerialSpdDenseSolver<T>::setVectors: X has an invalid stride!");

  resetVectors(); 
  LHS_ = X;
  RHS_ = B;
  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialSpdDenseSolver<OrdinalType,ScalarType>::estimateSolutionErrors(bool flag) 
{
  estimateSolutionErrors_ = flag;

  // If the errors are estimated, this implies that the solution must be refined
  refineSolution_ = refineSolution_ || flag;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::factor() {

  if (factored()) return(0); // Already factored

  TEUCHOS_TEST_FOR_EXCEPTION(inverted(), std::logic_error,
		     "SerialSpdDenseSolver<T>::factor: Cannot factor an inverted matrix!");

  ANORM_ = Matrix_->normOne(); // Compute 1-Norm of A


  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (refineSolution_ ) {
      Factor_ = rcp( new SerialSymDenseMatrix<OrdinalType,ScalarType>(*Matrix_) );
      AF_ = Factor_->values();
      LDAF_ = Factor_->stride();
    }
  
  int ierr = 0;
  if (equilibrate_) ierr = equilibrateMatrix();

  if (ierr!=0) return(ierr);
  
  INFO_ = 0;
  this->POTRF(Matrix_->UPLO(), numRowCols_, AF_, LDAF_, &INFO_);
  factored_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::solve() {

  // We will call one of four routines depending on what services the user wants and 
  // whether or not the matrix has been inverted or factored already.
  //
  // If the matrix has been inverted, use DGEMM to compute solution.
  // Otherwise, if the user want the matrix to be equilibrated or wants a refined solution, we will
  // call the X interface.
  // Otherwise, if the matrix is already factored we will call the TRS interface.
  // Otherwise, if the matrix is unfactored we will call the SV interface.

  int ierr = 0;
  if (equilibrate_) {
    ierr = equilibrateRHS();
    equilibratedB_ = true;
  }
  if (ierr != 0) return(ierr);  // Can't equilibrate B, so return.

  TEUCHOS_TEST_FOR_EXCEPTION( (equilibratedA_ && !equilibratedB_) || (!equilibratedA_ && equilibratedB_) , 
                     std::logic_error, "SerialSpdDenseSolver<T>::solve: Matrix and vectors must be similarly scaled!");
  TEUCHOS_TEST_FOR_EXCEPTION( RHS_==Teuchos::null, std::invalid_argument, 
                     "SerialSpdDenseSolver<T>::solve: No right-hand side vector (RHS) has been set for the linear system!");
  TEUCHOS_TEST_FOR_EXCEPTION( LHS_==Teuchos::null, std::invalid_argument, 
                     "SerialSpdDenseSolver<T>::solve: No solution vector (LHS) has been set for the linear system!");

  if (shouldEquilibrate() && !equilibratedA_)
    std::cout << "WARNING!  SerialSpdDenseSolver<T>::solve: System should be equilibrated!" << std::endl;

  if (inverted()) {

    TEUCHOS_TEST_FOR_EXCEPTION( RHS_->values() == LHS_->values(), std::invalid_argument, 
                        "SerialSpdDenseSolver<T>::solve: X and B must be different vectors if matrix is inverted.");

    INFO_ = 0;
    this->GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, numRowCols_, RHS_->numCols(), 
	       numRowCols_, 1.0, AF_, LDAF_, RHS_->values(), RHS_->stride(), 0.0, 
	       LHS_->values(), LHS_->stride());
    if (INFO_!=0) return(INFO_);
    solved_ = true;
  }
  else {

    if (!factored()) factor(); // Matrix must be factored
    
    if (RHS_->values()!=LHS_->values()) {
       (*LHS_) = (*RHS_); // Copy B to X if needed
    }
    INFO_ = 0;
    this->POTRS(Matrix_->UPLO(), numRowCols_, RHS_->numCols(), AF_, LDAF_, LHS_->values(), LHS_->stride(), &INFO_);
    if (INFO_!=0) return(INFO_);
    solved_ = true;

  }
  int ierr1=0;
  if (refineSolution_ && !inverted()) ierr1 = applyRefinement();
  if (ierr1!=0) 
    return(ierr1);
  
  if (equilibrate_) ierr1 = unequilibrateLHS();
  return(ierr1);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::applyRefinement()
{
  TEUCHOS_TEST_FOR_EXCEPTION(!solved(), std::logic_error,
		     "SerialSpdDenseSolver<T>::applyRefinement: Must have an existing solution!");
  TEUCHOS_TEST_FOR_EXCEPTION(A_==AF_, std::logic_error,
		     "SerialSpdDenseSolver<T>::applyRefinement: Cannot apply refinement if no original copy of A!");

  OrdinalType NRHS = RHS_->numCols();
  FERR_.resize( NRHS );
  BERR_.resize( NRHS );
  allocateWORK();
  allocateIWORK();
  
  INFO_ = 0;
  this->PORFS(Matrix_->UPLO(), numRowCols_, NRHS, A_, LDA_, AF_, LDAF_,
	      RHS_->values(), RHS_->stride(), LHS_->values(), LHS_->stride(), 
              &FERR_[0], &BERR_[0], &WORK_[0], &IWORK_[0], &INFO_);
  
  solutionErrorsEstimated_ = true;
  reciprocalConditionEstimated_ = true;
  solutionRefined_ = true;
  
  return(INFO_);

}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::computeEquilibrateScaling() 
{
  if (R_.size()!=0) return(0); // Already computed
 
  R_.resize( numRowCols_ );
  
  INFO_ = 0;
  this->POEQU (numRowCols_, AF_, LDAF_, &R_[0], &SCOND_, &AMAX_, &INFO_);
  if (SCOND_<0.1*ScalarTraits<MagnitudeType>::one() || 
      AMAX_ < ScalarTraits<ScalarType>::rmin() || 
      AMAX_ > ScalarTraits<ScalarType>::rmax()) 
    shouldEquilibrate_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::equilibrateMatrix()
{
  OrdinalType i, j;
  int ierr = 0;

  if (equilibratedA_) return(0); // Already done.
  if (R_.size()==0) ierr = computeEquilibrateScaling(); // Compute R if needed.
  if (ierr!=0) return(ierr);     // If return value is not zero, then can't equilibrate.
  if (Matrix_->upper()) {
    if (A_==AF_) {
      ScalarType * ptr;
      for (j=0; j<numRowCols_; j++) {
	ptr = A_ + j*LDA_;
	ScalarType s1 = R_[j];
	for (i=0; i<=j; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	}
      }
    }
    else {
      ScalarType * ptr;
      ScalarType * ptr1;
      for (j=0; j<numRowCols_; j++) {
	ptr = A_ + j*LDA_;
	ptr1 = AF_ + j*LDAF_;
	ScalarType s1 = R_[j];
	for (i=0; i<=j; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	  *ptr1 = *ptr1*s1*R_[i];
	  ptr1++;
	}
      }
    }
  }
  else {
    if (A_==AF_) {
      ScalarType * ptr;
      for (j=0; j<numRowCols_; j++) {
	ptr = A_ + j + j*LDA_;
	ScalarType s1 = R_[j];
	for (i=j; i<numRowCols_; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	}
      }
    }
    else {
      ScalarType * ptr;
      ScalarType * ptr1;
      for (j=0; j<numRowCols_; j++) {
	ptr = A_ + j + j*LDA_;
	ptr1 = AF_ + j + j*LDAF_;
	ScalarType s1 = R_[j];
	for (i=j; i<numRowCols_; i++) {
	  *ptr = *ptr*s1*R_[i];
	  ptr++;
	  *ptr1 = *ptr1*s1*R_[i];
	  ptr1++;
	}
      }
    }
  }
  
  equilibratedA_ = true;
  
  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::equilibrateRHS()
{
  OrdinalType i, j;
  int ierr = 0;

  if (equilibratedB_) return(0); // Already done.
  if (R_.size()==0) ierr = computeEquilibrateScaling(); // Compute R if needed.
  if (ierr!=0) return(ierr);     // Can't count on R being computed.

  OrdinalType LDB = RHS_->stride(), NRHS = RHS_->numCols();
  ScalarType * B = RHS_->values();
  ScalarType * ptr;
  for (j=0; j<NRHS; j++) {
    ptr = B + j*LDB;
    for (i=0; i<numRowCols_; i++) {
      *ptr = *ptr*R_[i];
      ptr++;
    }
  }
  
  equilibratedB_ = true;

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::unequilibrateLHS()
{
  OrdinalType i, j;

  if (!equilibratedB_) return(0); // Nothing to do

  OrdinalType LDX = LHS_->stride(), NLHS = LHS_->numCols();
  ScalarType * X = LHS_->values();
  ScalarType * ptr;
  for (j=0; j<NLHS; j++) {
    ptr = X + j*LDX;
    for (i=0; i<numRowCols_; i++) {
      *ptr = *ptr*R_[i];
      ptr++;
    }
  }

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::invert()
{
  if (!factored()) factor(); // Need matrix factored.

  INFO_ = 0;
  this->POTRI( Matrix_->UPLO(), numRowCols_, AF_, LDAF_, &INFO_);
  
  // Copy lower/upper triangle to upper/lower triangle to make full inverse.
  if (Matrix_->upper())
    Matrix_->setLower();
  else
    Matrix_->setUpper();

  inverted_ = true;
  factored_ = false;
  
  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialSpdDenseSolver<OrdinalType,ScalarType>::reciprocalConditionEstimate(MagnitudeType & Value)
{
  if (reciprocalConditionEstimated()) {
    Value = RCOND_;
    return(0); // Already computed, just return it.
  }

  if ( ANORM_<ScalarTraits<MagnitudeType>::zero() ) ANORM_ = Matrix_->normOne();

  int ierr = 0;
  if (!factored()) ierr = factor(); // Need matrix factored.
  if (ierr!=0) return(ierr);

  allocateWORK();
  allocateIWORK();

  // We will assume a one-norm condition number
  INFO_ = 0;
  this->POCON( Matrix_->UPLO(), numRowCols_, AF_, LDAF_, ANORM_, &RCOND_, &WORK_[0], &IWORK_[0], &INFO_);
  reciprocalConditionEstimated_ = true;
  Value = RCOND_;

  return(INFO_);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialSpdDenseSolver<OrdinalType,ScalarType>::Print(std::ostream& os) const {

  if (Matrix_!=Teuchos::null) os << "Solver Matrix"          << std::endl << *Matrix_ << std::endl;
  if (Factor_!=Teuchos::null) os << "Solver Factored Matrix" << std::endl << *Factor_ << std::endl;
  if (LHS_   !=Teuchos::null) os << "Solver LHS"             << std::endl << *LHS_    << std::endl;
  if (RHS_   !=Teuchos::null) os << "Solver RHS"             << std::endl << *RHS_    << std::endl;

}

} // namespace Teuchos

#endif /* TEUCHOS_SERIALSPDDENSESOLVER_H */

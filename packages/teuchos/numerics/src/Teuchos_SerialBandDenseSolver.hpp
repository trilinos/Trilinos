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

#ifndef _TEUCHOS_SERIALBANDDENSESOLVER_HPP_
#define _TEUCHOS_SERIALBANDDENSESOLVER_HPP_
/// \file Teuchos_SerialBandDenseSolver.hpp
///
/// Declaration and definition of SerialBandDenseSolver,
/// a templated class for solving banded dense linear systems.

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialBandDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Teuchos {

/*! \class SerialBandDenseSolver
  \brief A class for representing and solving banded dense linear systems.
  \tparam OrdinalType The index type used by the linear algebra implementation.
    This should always be \c int.
  \tparam ScalarType The type of entries in the matrix.

  \section Teuchos_SerialBandDenseSolver_Intro Introduction

  This class defines a banded dense matrix, which may have any number
  of rows or columns (not necessarily equal).  It's called "serial"
  because the matrix lives in a single memory space.  Thus, it's the
  kind of matrix that one might give to the BLAS or LAPACK, not a
  distributed matrix like one would give to ScaLAPACK.

  This class also has methods for computing the (banded) LU
  factorization of the matrix, and solving linear systems with the
  matrix.  We use instances of SerialDenseVector to represent the
  right-hand side b or the solution x in the linear system \f$Ax=b\f$.
  The instance of this class used store the banded matrix must contain
  KL extra superdiagonals to store the L and U factors (see details
  below).

  Users have the option to do equilibration before factoring the matrix.
  This may improve accuracy when solving ill-conditioned problems.

  \section Teuchos_SerialBandDenseSolver_LAPACK SerialBandDenseSolver and LAPACK

  Teuchos' LAPACK class wraps LAPACK's LU factorizations, including
  the banded factorizations.  It thus gives access to most of the same
  functionality as SerialBandDenseSolver.  The main difference is that
  SerialBandDenseSolver offers a higher level of abstraction.  It
  hides the details of which LAPACK routines to call.  Furthermore, if
  you have built Teuchos with support for the third-party library
  Eigen, SerialBandDenseSolver lets you solve linear systems for
  <tt>ScalarType</tt> other than the four supported by the LAPACK
  library.

  \section Teuchos_SerialBandDenseSolver_Construct Constructing SerialBandDenseSolver objects

  There is a single Teuchos::SerialBandDenseSolver constructor.
  However, the matrix, right hand side and solution vectors must be
  set prior to executing most methods in this class.

  \subsection Teuchos_SerialBandDenseSolver_Construct_Setting Setting vectors used for linear solves

  The matrix A, the left hand side X and the right hand side B (when
  solving AX = B, for X), can be set by appropriate set methods.  Each
  of these three objects must be a SerialDenseMatrix or a
  SerialDenseVector object.  The set methods are as follows:
  - setMatrix(): Sets the matrix
  - setVectors() - Sets the left and right hand side vector(s)

  \subsection Teuchos_SerialBandDenseSolver_Construct_Format Format of the matrix A

  The SerialBandDenseMatrix must contain KL extra superdiagonals to store the L and U factors, where KL
  is the upper bandwidth. Consider using the non-member conversion routines generalToBanded and bandedToGeneral if the
  full SerialDenseMatrix is already in storage. However, it is more efficient simply to construct the
  SerialBandDenseMatrix with the desired parameters and use the provided matrix access operators so
  that the full rectangular matrix need not be stored. The conversion routine generalToBanded has a flag to store
  the input Teuchos::SerialDenseMatrix in banded format with KL extra superdiagonals so this class can use it. Again,
  it is more efficient to simply construct a Teuchos::SerialBandDenseMatrix object with KL extra superdiagonals than are
  needed for the matrix data and fill the matrix using the matrix access operators.

  See the documentation of Teuchos::SerialBandDenseMatrix for further details on the storage format.

  \section Teuchos_SerialBandDenseSolver_Util Vector and Utility Functions

  Once a Teuchos::SerialBandDenseSolver is constructed, several mathematical functions can be applied to
  the object.  Specifically:
  <ul>
  <li> Conversion between storage formats
  <li> Factorizations
  <li> Solves
  <li> Condition estimates
  <li> Equilibration
  <li> Norms
  </ul>

  \section Teuchos_SerialBandDenseSolver_Strategies Strategies for Solving Linear Systems

  In many cases, linear systems can be accurately solved by simply computing the LU factorization
  of the matrix and then performing a forward back solve with a given set of right hand side vectors.  However,
  in some instances, the factorization may be very poorly conditioned and this simple approach may not work.  In
  these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
  the factorization.

  SerialBandDenseSolver will use equilibration with the factorization if, once the object
  is constructed and \e before it is factored, you call the function factorWithEquilibration(true) to force
  equilibration to be used.  If you are uncertain if equilibration should be used, you may call the function
  shouldEquilibrate() which will return true if equilibration could possibly help.  shouldEquilibrate() uses
  guidelines specified in the LAPACK User Guide, namely if SCOND < 0.1 and AMAX < Underflow or AMAX > Overflow, to
  determine if equilibration \e might be useful.

  SerialBandDenseSolver will use iterative refinement after a forward/back solve if you call
  solveToRefinedSolution(true).  It will also compute forward and backward error estimates if you call
  estimateSolutionErrors(true).  Access to the forward (back) error estimates is available via FERR() (BERR()).

  Examples using SerialBandDenseSolver can be found in the Teuchos test directories.
*/

  template<typename OrdinalType, typename ScalarType>
  class SerialBandDenseSolver : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>,
				public LAPACK<OrdinalType, ScalarType>
  {

  public:

    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

    //! @name Constructor/Destructor Methods
    //@{

    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    SerialBandDenseSolver();

    //! SerialBandDenseSolver destructor.
    virtual ~SerialBandDenseSolver();

    //@}
    //! @name Set Methods
    //@{

    //! Sets the pointer for coefficient matrix
    int setMatrix(const RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> >& AB);

    //! Sets the pointers for left and right hand side vector(s).
    /*! Row dimension of X must match column dimension of matrix A, row dimension of B
      must match row dimension of A.  X and B must have the same dimensions.
    */
    int setVectors(const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& X,
		   const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& B);

    //@}
    //! @name Strategy Modifying Methods
    //@{

    /// Set whether or not to equilibrate just before the matrix factorization.
    /// This function must be called before the factorization is performed.
    void factorWithEquilibration(bool flag) {equilibrate_ = flag; return;}

    /// If \c flag is true, causes all subsequent function calls to work with the transpose of \e this matrix, otherwise not.
    ///
    /// \note This interface will not work correctly for complex-valued linear systems, use solveWithTransposeFlag().
    void solveWithTranspose(bool flag) {transpose_ = flag; if (flag) TRANS_ = Teuchos::TRANS; else TRANS_ = Teuchos::NO_TRANS; return;}

    /// All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS, \c Teuchos::TRANS, and \c Teuchos::CONJ_TRANS).
    /// \note This interface will allow correct behavior for complex-valued linear systems, solveWithTranspose() will not.
    void solveWithTransposeFlag(Teuchos::ETransp trans) {TRANS_ = trans; if (trans != Teuchos::NO_TRANS) {  transpose_ = true; } }

    //! Set whether or not to use iterative refinement to improve solutions to linear systems.
    void solveToRefinedSolution(bool flag) {refineSolution_ = flag; return;}

    //! Causes all solves to estimate the forward and backward solution error.
    /*! Error estimates will be in the arrays FERR and BERR, resp, after the solve step is complete.
      These arrays are accessible via the FERR() and BERR() access functions.
    */
    void estimateSolutionErrors(bool flag);
    //@}

    //! @name Factor/Solve Methods
    //@{

    //! Computes the in-place LU factorization of the matrix using the LAPACK routine \e _GBTRF.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor();

    //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve();

    //! Computes the scaling vector S(i) = 1/sqrt(A(i,i)) of the \e this matrix.
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
    RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> > getMatrix()  const {return(Matrix_);}

    //! Returns pointer to factored matrix (assuming factorization has been performed).
    RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> > getFactoredMatrix()  const {return(Factor_);}

    //! Returns pointer to current LHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getLHS() const {return(LHS_);}

    //! Returns pointer to current RHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getRHS() const {return(RHS_);}

    //! Returns row dimension of system.
    OrdinalType numRows()  const {return(M_);}

    //! Returns column dimension of system.
    OrdinalType numCols()  const {return(N_);}

    //! Returns lower bandwidth of system.
    OrdinalType numLower()  const {return(KL_);}

    //! Returns upper bandwidth of system.
    OrdinalType numUpper()  const {return(KU_);}

    //! Returns pointer to pivot vector (if factorization has been computed), zero otherwise.
    std::vector<OrdinalType> IPIV()  const {return(IPIV_);}

    //! Returns the 1-Norm of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType ANORM()  const {return(ANORM_);}

    //! Returns the reciprocal of the condition number of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType RCOND()  const {return(RCOND_);}

    //! Ratio of smallest to largest row scale factors for the \e this matrix (returns -1 if not yet computed).
    /*! If ROWCND() is >= 0.1 and AMAX() is not close to overflow or underflow, then equilibration is not needed.
     */
    MagnitudeType ROWCND()  const {return(ROWCND_);}

    //! Ratio of smallest to largest column scale factors for the \e this matrix (returns -1 if not yet computed).
    /*! If COLCND() is >= 0.1 then equilibration is not needed.
     */
    MagnitudeType COLCND()  const {return(COLCND_);}

    //! Returns the absolute value of the largest entry of the \e this matrix (returns -1 if not yet computed).
    MagnitudeType AMAX()  const {return(AMAX_);}

    //! Returns a pointer to the forward error estimates computed by LAPACK.
    std::vector<MagnitudeType> FERR()  const {return(FERR_);}

    //! Returns a pointer to the backward error estimates computed by LAPACK.
    std::vector<MagnitudeType> BERR()  const {return(BERR_);}

    //! Returns a pointer to the row scaling vector used for equilibration.
    std::vector<MagnitudeType> R()  const {return(R_);}

    //! Returns a pointer to the column scale vector used for equilibration.
    std::vector<MagnitudeType> C()  const {return(C_);}
    //@}

    //! @name I/O methods
    //@{
    //! Print service methods; defines behavior of ostream << operator.
    void Print(std::ostream& os) const;
    //@}
  protected:

    void allocateWORK() { LWORK_ = 3*N_; WORK_.resize( LWORK_ ); return;}
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
    bool reciprocalConditionEstimated_;
    bool refineSolution_;
    bool solutionRefined_;

    Teuchos::ETransp TRANS_;

    OrdinalType M_;
    OrdinalType N_;
    OrdinalType KL_;
    OrdinalType KU_;
    OrdinalType Min_MN_;
    OrdinalType LDA_;
    OrdinalType LDAF_;
    OrdinalType INFO_;
    OrdinalType LWORK_;

    std::vector<OrdinalType> IPIV_;
    std::vector<int> IWORK_;

    MagnitudeType ANORM_;
    MagnitudeType RCOND_;
    MagnitudeType ROWCND_;
    MagnitudeType COLCND_;
    MagnitudeType AMAX_;

    RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> > Matrix_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > LHS_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > RHS_;
    RCP<SerialBandDenseMatrix<OrdinalType, ScalarType> > Factor_;

    ScalarType * A_;
    ScalarType * AF_;
    std::vector<MagnitudeType> FERR_;
    std::vector<MagnitudeType> BERR_;
    std::vector<ScalarType> WORK_;
    std::vector<MagnitudeType> R_;
    std::vector<MagnitudeType> C_;

  private:

    // SerialBandDenseSolver copy constructor (put here because we don't want user access)
    SerialBandDenseSolver(const SerialBandDenseSolver<OrdinalType, ScalarType>& Source);
    SerialBandDenseSolver & operator=(const SerialBandDenseSolver<OrdinalType, ScalarType>& Source);

  };

  // Helper traits to distinguish work arrays for real and complex-valued datatypes.
  using namespace details;

//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialBandDenseSolver<OrdinalType,ScalarType>::SerialBandDenseSolver()
  : CompObject(),
    Object("Teuchos::SerialBandDenseSolver"),
    equilibrate_(false),
    shouldEquilibrate_(false),
    equilibratedA_(false),
    equilibratedB_(false),
    transpose_(false),
    factored_(false),
    estimateSolutionErrors_(false),
    solutionErrorsEstimated_(false),
    solved_(false),
    reciprocalConditionEstimated_(false),
    refineSolution_(false),
    solutionRefined_(false),
    TRANS_(Teuchos::NO_TRANS),
    M_(0),
    N_(0),
    KL_(0),
    KU_(0),
    Min_MN_(0),
    LDA_(0),
    LDAF_(0),
    INFO_(0),
    LWORK_(0),
    RCOND_(ScalarTraits<MagnitudeType>::zero()),
    ROWCND_(ScalarTraits<MagnitudeType>::zero()),
    COLCND_(ScalarTraits<MagnitudeType>::zero()),
    AMAX_(ScalarTraits<MagnitudeType>::zero()),
    A_(0),
    AF_(0)
{
  resetMatrix();
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialBandDenseSolver<OrdinalType,ScalarType>::~SerialBandDenseSolver()
{}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialBandDenseSolver<OrdinalType,ScalarType>::resetVectors()
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
void SerialBandDenseSolver<OrdinalType,ScalarType>::resetMatrix()
{
  resetVectors();
  equilibratedA_ = false;
  factored_ = false;
  M_ = 0;
  N_ = 0;
  KL_ = 0;
  KU_ = 0;
  Min_MN_ = 0;
  LDA_ = 0;
  LDAF_ = 0;
  RCOND_ = -ScalarTraits<MagnitudeType>::one();
  ROWCND_ = -ScalarTraits<MagnitudeType>::one();
  COLCND_ = -ScalarTraits<MagnitudeType>::one();
  AMAX_ = -ScalarTraits<MagnitudeType>::one();
  A_ = 0;
  AF_ = 0;
  INFO_ = 0;
  LWORK_ = 0;
  R_.resize(0);
  C_.resize(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::setMatrix(const RCP<SerialBandDenseMatrix<OrdinalType,ScalarType> >& AB)
{

  OrdinalType m = AB->numRows();
  OrdinalType n = AB->numCols();
  OrdinalType kl = AB->lowerBandwidth();
  OrdinalType ku = AB->upperBandwidth();

  // Check that the new matrix is consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(AB->values()==0, std::invalid_argument,
		     "SerialBandDenseSolver<T>::setMatrix: A is an empty SerialBandDenseMatrix<T>!");

  resetMatrix();
  Matrix_ = AB;
  Factor_ = AB;
  M_ = m;
  N_ = n;
  KL_ = kl;
  KU_ = ku-kl;
  Min_MN_ = TEUCHOS_MIN(M_,N_);
  LDA_ = AB->stride();
  LDAF_ = LDA_;
  A_ = AB->values();
  AF_ = AB->values();

  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::setVectors(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& X,
							   const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& B)
{
  // Check that these new vectors are consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(B->numRows()!=X->numRows() || B->numCols() != X->numCols(), std::invalid_argument,
		     "SerialBandDenseSolver<T>::setVectors: X and B are not the same size!");
  TEUCHOS_TEST_FOR_EXCEPTION(B->values()==0, std::invalid_argument,
		     "SerialBandDenseSolver<T>::setVectors: B is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(X->values()==0, std::invalid_argument,
		     "SerialBandDenseSolver<T>::setVectors: X is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(B->stride()<1, std::invalid_argument,
		     "SerialBandDenseSolver<T>::setVectors: B has an invalid stride!");
  TEUCHOS_TEST_FOR_EXCEPTION(X->stride()<1, std::invalid_argument,
		     "SerialBandDenseSolver<T>::setVectors: X has an invalid stride!");

  resetVectors();
  LHS_ = X;
  RHS_ = B;
  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialBandDenseSolver<OrdinalType,ScalarType>::estimateSolutionErrors(bool flag)
{
  estimateSolutionErrors_ = flag;

  // If the errors are estimated, this implies that the solution must be refined
  refineSolution_ = refineSolution_ || flag;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::factor() {

  if (factored()) return(0); // Already factored

  ANORM_ = Matrix_->normOne(); // Compute 1-Norm of A

  // If we want to refine the solution, then the factor must
  // be stored separatedly from the original matrix

  if (A_ == AF_)
    if (refineSolution_ ) {
      Factor_ = rcp( new SerialBandDenseMatrix<OrdinalType,ScalarType>(*Matrix_) );
      AF_ = Factor_->values();
      LDAF_ = Factor_->stride();
    }

  int ierr = 0;
  if (equilibrate_) ierr = equilibrateMatrix();

  if (ierr!=0) return(ierr);

  if (IPIV_.size()==0) IPIV_.resize( N_ ); // Allocated Pivot vector if not already done.

  INFO_ = 0;

  this->GBTRF(M_, N_, KL_, KU_, AF_, LDAF_, &IPIV_[0], &INFO_);
  factored_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::solve() {

  // If the user want the matrix to be equilibrated or wants a refined solution, we will
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
		     std::logic_error, "SerialBandDenseSolver<T>::solve: Matrix and vectors must be similarly scaled!");
  TEUCHOS_TEST_FOR_EXCEPTION( RHS_==Teuchos::null, std::invalid_argument,
		     "SerialBandDenseSolver<T>::solve: No right-hand side vector (RHS) has been set for the linear system!");
  TEUCHOS_TEST_FOR_EXCEPTION( LHS_==Teuchos::null, std::invalid_argument,
		     "SerialBandDenseSolver<T>::solve: No solution vector (LHS) has been set for the linear system!");

  if (shouldEquilibrate() && !equilibratedA_)
    std::cout << "WARNING!  SerialBandDenseSolver<T>::solve: System should be equilibrated!" << std::endl;

  if (!factored()) factor(); // Matrix must be factored

  if (RHS_->values()!=LHS_->values()) {
    (*LHS_) = (*RHS_); // Copy B to X if needed
  }
  INFO_ = 0;

  this->GBTRS(ETranspChar[TRANS_], N_, KL_, KU_, RHS_->numCols(), AF_, LDAF_, &IPIV_[0], LHS_->values(), LHS_->stride(), &INFO_);

  if (INFO_!=0) return(INFO_);
  solved_ = true;

  int ierr1=0;
  if (refineSolution_) ierr1 = applyRefinement();
  if (ierr1!=0)
    return(ierr1);

  if (equilibrate_) ierr1 = unequilibrateLHS();

  return(ierr1);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::applyRefinement()
{
  TEUCHOS_TEST_FOR_EXCEPTION(!solved(), std::logic_error,
		     "SerialBandDenseSolver<T>::applyRefinement: Must have an existing solution!");
  TEUCHOS_TEST_FOR_EXCEPTION(A_==AF_, std::logic_error,
		     "SerialBandDenseSolver<T>::applyRefinement: Cannot apply refinement if no original copy of A!");

  OrdinalType NRHS = RHS_->numCols();
  FERR_.resize( NRHS );
  BERR_.resize( NRHS );
  allocateWORK();

  INFO_ = 0;
  std::vector<typename details::lapack_traits<ScalarType>::iwork_type> GBRFS_WORK( N_ );
  this->GBRFS(ETranspChar[TRANS_], N_, KL_, KU_, NRHS, A_+KL_, LDA_, AF_, LDAF_, &IPIV_[0],
	      RHS_->values(), RHS_->stride(), LHS_->values(), LHS_->stride(),
	      &FERR_[0], &BERR_[0], &WORK_[0], &GBRFS_WORK[0], &INFO_);

  solutionErrorsEstimated_ = true;
  reciprocalConditionEstimated_ = true;
  solutionRefined_ = true;

  return(INFO_);

}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::computeEquilibrateScaling()
{
  if (R_.size()!=0) return(0); // Already computed

  R_.resize( M_ );
  C_.resize( N_ );

  INFO_ = 0;
  this->GBEQU (M_, N_, KL_, KU_, AF_+KL_, LDAF_, &R_[0], &C_[0], &ROWCND_, &COLCND_, &AMAX_, &INFO_);

  if (COLCND_<0.1*ScalarTraits<MagnitudeType>::one() ||
      ROWCND_<0.1*ScalarTraits<MagnitudeType>::one() ||
      AMAX_ < ScalarTraits<ScalarType>::rmin() || AMAX_ > ScalarTraits<ScalarType>::rmax() )
    shouldEquilibrate_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::equilibrateMatrix()
{
  OrdinalType i, j;
  int ierr = 0;

  if (equilibratedA_) return(0); // Already done.
  if (R_.size()==0) ierr = computeEquilibrateScaling(); // Compute R and C if needed.
  if (ierr!=0) return(ierr);     // If return value is not zero, then can't equilibrate.

  if (A_==AF_) {

    ScalarType * ptr;
    for (j=0; j<N_; j++) {
      ptr = A_ + KL_ + j*LDA_ + std::max(KU_-j, 0);
      ScalarType s1 = C_[j];
      for (i=std::max(0,j-KU_); i<=std::min(M_-1,j+KL_); i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
      }
    }
  } else {

    ScalarType * ptr;
    ScalarType * ptrL;
    ScalarType * ptrU;
    for (j=0; j<N_; j++) {
      ptr = A_ + KL_ + j*LDA_ + std::max(KU_-j, 0);
      ScalarType s1 = C_[j];
      for (i=std::max(0,j-KU_); i<=std::min(M_-1,j+KL_); i++) {
	*ptr = *ptr*s1*R_[i];
	ptr++;
      }
    }
    for (j=0; j<N_; j++) {
      ptrU = AF_ + j*LDAF_ + std::max(KL_+KU_-j, 0);
      ptrL = AF_ + KL_ + KU_ + 1 + j*LDAF_;
      ScalarType s1 = C_[j];
      for (i=std::max(0,j-(KL_+KU_)); i<=std::min(M_-1,j); i++) {
	*ptrU = *ptrU*s1*R_[i];
	ptrU++;
      }
      for (i=std::max(0,j); i<=std::min(M_-1,j+KL_); i++) {
	*ptrL = *ptrL*s1*R_[i];
	ptrL++;
      }
    }
  }

  equilibratedA_ = true;

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::equilibrateRHS()
{
  OrdinalType i, j;
  int ierr = 0;

  if (equilibratedB_) return(0); // Already done.
  if (R_.size()==0) ierr = computeEquilibrateScaling(); // Compute R and C if needed.
  if (ierr!=0) return(ierr);     // Can't count on R and C being computed.

  MagnitudeType * R_tmp = &R_[0];
  if (transpose_) R_tmp = &C_[0];

  OrdinalType LDB = RHS_->stride(), NRHS = RHS_->numCols();
  ScalarType * B = RHS_->values();
  ScalarType * ptr;
  for (j=0; j<NRHS; j++) {
    ptr = B + j*LDB;
    for (i=0; i<M_; i++) {
      *ptr = *ptr*R_tmp[i];
      ptr++;
    }
  }

  equilibratedB_ = true;

  return(0);
}


//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::unequilibrateLHS()
{
  OrdinalType i, j;

  if (!equilibratedB_) return(0); // Nothing to do

  MagnitudeType * C_tmp = &C_[0];
  if (transpose_) C_tmp = &R_[0];

  OrdinalType LDX = LHS_->stride(), NLHS = LHS_->numCols();
  ScalarType * X = LHS_->values();
  ScalarType * ptr;
  for (j=0; j<NLHS; j++) {
    ptr = X + j*LDX;
    for (i=0; i<N_; i++) {
      *ptr = *ptr*C_tmp[i];
      ptr++;
    }
  }

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialBandDenseSolver<OrdinalType,ScalarType>::reciprocalConditionEstimate(MagnitudeType & Value)
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

  // We will assume a one-norm condition number
  INFO_ = 0;
  std::vector<typename details::lapack_traits<ScalarType>::iwork_type> GBCON_WORK( N_ );
  this->GBCON( '1', N_, KL_, KU_, AF_, LDAF_, &IPIV_[0], ANORM_, &RCOND_, &WORK_[0], &GBCON_WORK[0], &INFO_);
  reciprocalConditionEstimated_ = true;
  Value = RCOND_;

  return(INFO_);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialBandDenseSolver<OrdinalType,ScalarType>::Print(std::ostream& os) const {

  if (Matrix_!=Teuchos::null) os << "Solver Matrix"          << std::endl << *Matrix_ << std::endl;
  if (Factor_!=Teuchos::null) os << "Solver Factored Matrix" << std::endl << *Factor_ << std::endl;
  if (LHS_   !=Teuchos::null) os << "Solver LHS"             << std::endl << *LHS_    << std::endl;
  if (RHS_   !=Teuchos::null) os << "Solver RHS"             << std::endl << *RHS_    << std::endl;

}

//=============================================================================


//=============================================================================


} // namespace Teuchos

#endif /* _TEUCHOS_SERIALBANDDENSESOLVER_HPP_ */

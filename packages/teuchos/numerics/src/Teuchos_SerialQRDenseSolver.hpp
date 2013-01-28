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

#ifndef _TEUCHOS_SERIALQRDENSESOLVER_HPP_
#define _TEUCHOS_SERIALQRDENSESOLVER_HPP_
/*! \file Teuchos_SerialQRDenseSolver.hpp
  \brief Templated class for solving dense linear problems.
*/

#include "Teuchos_CompObject.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_TEUCHOS_EIGEN
#include "Eigen/Dense"
#endif

/*! \class Teuchos::SerialQRDenseSolver
  \brief A class for solving dense linear problems.

  The Teuchos::SerialQRDenseSolver class enables the definition, in terms of Teuchos::SerialDenseMatrix
  and Teuchos::SerialDenseVector objects, of a dense linear problem, followed by the solution of that
  problem via the most sophisticated techniques available in LAPACK.

  The Teuchos::SerialQRDenseSolver class is intended to provide full-featured support for solving linear
  problems for general dense rectangular (or square) matrices.  It is written on top of BLAS and LAPACK
  and thus has excellent performance and numerical capabilities.  Using this class, one can either perform
  simple factorizations and solves or apply all the tricks available in LAPACK to get the best possible
  solution for very ill-conditioned problems.

  <b>Teuchos::SerialQRDenseSolver vs. Teuchos::LAPACK</b>

  The Teuchos::LAPACK class provides access to most of the same functionality as Teuchos::SerialQRDenseSolver.
  The primary difference is that Teuchos::LAPACK is a "thin" layer on top of LAPACK and Teuchos::SerialQRDenseSolver
  attempts to provide easy access to the more sophisticated aspects of solving dense linear and eigensystems.
  <ul>
  <li> When you should use Teuchos::LAPACK:  If you are simply looking for a convenient wrapper around
  the Fortran LAPACK routines and you have a well-conditioned problem, you should probably use Teuchos::LAPACK directly.
  <li> When you should use Teuchos::SerialQRDenseSolver: If you want to (or potentially want to) solve
  ill-conditioned problems or want to work with a more object-oriented interface, you should probably use
  Teuchos::SerialQRDenseSolver.
  </ul>

  <b>Constructing Teuchos::SerialQRDenseSolver Objects</b>

  There is a single Teuchos::SerialQRDenseSolver constructor.   However, the matrix, right hand side and solution
  vectors must be set prior to executing most methods in this class.

  <b>Setting vectors used for linear solves</b>

  The matrix A, the left hand side X and the right hand side B (when solving AX = B, for X), can be set by appropriate set
  methods.  Each of these three objects must be an Teuchos::SerialDenseMatrix or and Teuchos::SerialDenseVector object.  The
  set methods are as follows:
  <ul>
  <li> setMatrix()  - Sets the matrix.
  <li> setVectors() - Sets the left and right hand side vector(s).
  </ul>

  <b>Vector and Utility Functions</b>

  Once a Teuchos::SerialQRDenseSolver is constructed, several mathematical functions can be applied to
  the object.  Specifically:
  <ul>
  <li> Factorizations.
  <li> Solves.
  <li> Equilibration.
  <li> Norms.
  </ul>

  <b>Strategies for Solving Linear Systems</b>
  In many cases, linear least squares systems can be accurately solved by simply computing the QR factorization
  of the matrix and then performing a forward back solve with a given set of right hand side vectors.  However,
  in some instances, the factorization may be very poorly conditioned and this simple approach may not work.  In
  these situations, equilibration and iterative refinement may improve the accuracy, or prevent a breakdown in
  the factorization.

  Teuchos::SerialQRDenseSolver will use equilibration with the factorization if, once the object
  is constructed and \e before it is factored, you call the function factorWithEquilibration(true) to force
  equilibration to be used.  If you are uncertain if equilibration should be used, you may call the function
  shouldEquilibrate() which will return true if equilibration could possibly help.  shouldEquilibrate() uses
  guidelines specified in the LAPACK User Guide to determine if equilibration \e might be useful.

  Examples using Teuchos::SerialQRDenseSolver can be found in the Teuchos test directories.
*/

namespace Teuchos {

  template<typename OrdinalType, typename ScalarType>
  class SerialQRDenseSolver : public CompObject, public Object, public BLAS<OrdinalType, ScalarType>,
			      public LAPACK<OrdinalType, ScalarType>
  {

  public:

    typedef typename ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
#ifdef HAVE_TEUCHOS_EIGEN
    typedef typename Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> EigenMatrix;
    typedef typename Eigen::Matrix<ScalarType,Eigen::Dynamic,1> EigenVector;
    typedef typename Eigen::InnerStride<Eigen::Dynamic> EigenInnerStride;
    typedef typename Eigen::OuterStride<Eigen::Dynamic> EigenOuterStride;
    typedef typename Eigen::Map<EigenVector,0,EigenInnerStride > EigenVectorMap;
    typedef typename Eigen::Map<const EigenVector,0,EigenInnerStride > EigenConstVectorMap;
    typedef typename Eigen::Map<EigenMatrix,0,EigenOuterStride > EigenMatrixMap;
    typedef typename Eigen::Map<const EigenMatrix,0,EigenOuterStride > EigenConstMatrixMap;
    typedef typename Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,OrdinalType> EigenPermutationMatrix;
    typedef typename Eigen::Array<OrdinalType,Eigen::Dynamic,1> EigenOrdinalArray;
    typedef typename Eigen::Map<EigenOrdinalArray> EigenOrdinalArrayMap;
    typedef typename Eigen::Array<ScalarType,Eigen::Dynamic,1> EigenScalarArray;
    typedef typename Eigen::Map<EigenScalarArray> EigenScalarArrayMap;
#endif

    //! @name Constructor/Destructor Methods
    //@{
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    SerialQRDenseSolver();

    //! SerialQRDenseSolver destructor.
    virtual ~SerialQRDenseSolver();
    //@}

    //! @name Set Methods
    //@{

    //! Sets the pointers for coefficient matrix
    /*! Row dimension of A must be greater than or equal to the column dimension of A.
    */
    int setMatrix(const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& A);

    //! Sets the pointers for left and right hand side vector(s).
    /*! Row dimension of X must match column dimension of matrix A, row dimension of B
      must match row dimension of A.
    */
    int setVectors(const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& X,
		   const RCP<SerialDenseMatrix<OrdinalType, ScalarType> >& B);
    //@}

    //! @name Strategy Modifying Methods
    //@{

    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! \note This method must be called before the factorization is performed, otherwise it will have no effect.
   */
    void factorWithEquilibration(bool flag) {equilibrate_ = flag; return;}

    //! If \c flag is true, causes all subsequent function calls to work with the adjoint of \e this matrix, otherwise not.
    void solveWithTranspose(bool flag) {transpose_ = flag; if (flag) TRANS_ = Teuchos::CONJ_TRANS; else TRANS_ = Teuchos::NO_TRANS; return;}

    //! All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS or Teuchos::CONJ_TRANS).
    void solveWithTransposeFlag(Teuchos::ETransp trans) {TRANS_ = trans; if (trans != Teuchos::NO_TRANS) {  transpose_ = true; } }

    //@}

    //! @name Factor/Solve/Invert Methods
    //@{

    //! Computes the in-place QR factorization of the matrix using the LAPACK routine \e _GETRF or the Eigen class \e HouseholderQR
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor();

    //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve();

    //! Determines if \e this matrix should be scaled.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int computeEquilibrateScaling();

    //! Equilibrates the \e this matrix.
    /*!
      \note This method will be called automatically in solve() method if factorWithEquilibration( true ) is called.
      \return Integer error code, set to 0 if successful.
    */
    int equilibrateMatrix();

    //! Equilibrates the current RHS.
    /*!
      \note This method will be called automatically in solve() method if factorWithEquilibration( true ) is called.
      \return Integer error code, set to 0 if successful.
    */
    int equilibrateRHS();

    //! Unscales the solution vectors if equilibration was used to solve the system.
    /*!
      \note This method will be called automatically in solve() method if factorWithEquilibration( true ) is called.
      \return Integer error code, set to 0 if successful.
    */
    int unequilibrateLHS();

    //! Explicitly forms the unitary matrix Q.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int formQ();

    //! Explicitly forms the upper triangular matrix R.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int formR();

    //! Left multiply the input matrix by the unitary matrix Q or its adjoint.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int multiplyQ (ETransp transq, SerialDenseMatrix<OrdinalType, ScalarType> &C);

    //! Solve input matrix on the left with the upper triangular matrix R or its adjoint.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solveR (ETransp transr, SerialDenseMatrix<OrdinalType, ScalarType> &C);
    //@}

    //! @name Query methods
    //@{

    //! Returns true if adjoint of \e this matrix has and will be used.
    bool transpose() {return(transpose_);}

    //! Returns true if matrix is factored (factor available via getFactoredMatrix()).
    bool factored() {return(factored_);}

    //! Returns true if factor is equilibrated (factor available via getFactoredMatrix()).
    bool equilibratedA() {return(equilibratedA_);}

    //! Returns true if RHS is equilibrated (RHS available via getRHS()).
    bool equilibratedB() {return(equilibratedB_);}

    //! Returns true if the LAPACK general rules for equilibration suggest you should equilibrate the system.
    bool shouldEquilibrate() {computeEquilibrateScaling(); return(shouldEquilibrate_);}

    //! Returns true if the current set of vectors has been solved.
    bool solved() {return(solved_);}

    //! Returns true if Q has been formed explicitly.
    bool formedQ() {return(formedQ_);}

    //! Returns true if R has been formed explicitly.
    bool formedR() {return(formedR_);}

    //@}

    //! @name Data Accessor methods
    //@{

    //! Returns pointer to current matrix.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getMatrix()  const {return(Matrix_);}

    //! Returns pointer to factored matrix (assuming factorization has been performed).
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getFactoredMatrix()  const {return(Factor_);}

    //! Returns pointer to Q (assuming factorization has been performed).
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getQ()  const {return(FactorQ_);}

    //! Returns pointer to R (assuming factorization has been performed).
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getR()  const {return(FactorR_);}

    //! Returns pointer to current LHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getLHS() const {return(LHS_);}

    //! Returns pointer to current RHS.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getRHS() const {return(RHS_);}

    //! Returns row dimension of system.
    OrdinalType numRows()  const {return(M_);}

    //! Returns column dimension of system.
    OrdinalType numCols()  const {return(N_);}

    //! Returns pointer to pivot vector (if factorization has been computed), zero otherwise.
    std::vector<ScalarType> tau()  const {return(TAU_);}

    //! Returns the absolute value of the largest element of \e this matrix (returns -1 if not yet computed).
    MagnitudeType ANORM()  const {return(ANORM_);}

    //@}

    //! @name I/O methods
    //@{
    //! Print service methods; defines behavior of ostream << operator.
    void Print(std::ostream& os) const;
    //@}
  protected:

    void allocateWORK() { LWORK_ = 4*N_; WORK_.resize( LWORK_ ); return;}
    void resetMatrix();
    void resetVectors();


    bool equilibrate_;
    bool shouldEquilibrate_;
    bool equilibratedA_;
    bool equilibratedB_;
    OrdinalType equilibrationModeA_;
    OrdinalType equilibrationModeB_;
    bool transpose_;
    bool factored_;
    bool solved_;
    bool matrixIsZero_;
    bool formedQ_;
    bool formedR_;

    Teuchos::ETransp TRANS_;

    OrdinalType M_;
    OrdinalType N_;
    OrdinalType Min_MN_;
    OrdinalType LDA_;
    OrdinalType LDAF_;
    OrdinalType LDQ_;
    OrdinalType LDR_;
    OrdinalType INFO_;
    OrdinalType LWORK_;

    std::vector<ScalarType> TAU_;

    MagnitudeType ANORM_;
    MagnitudeType BNORM_;

    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > Matrix_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > LHS_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > TMP_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > RHS_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > Factor_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > FactorQ_;
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > FactorR_;

    ScalarType * A_;
    ScalarType * AF_;
    ScalarType * Q_;
    ScalarType * R_;
    std::vector<ScalarType> WORK_;
#ifdef HAVE_TEUCHOS_EIGEN
    Eigen::HouseholderQR<EigenMatrix> qr_;
#endif

  private:
    // SerialQRDenseSolver copy constructor (put here because we don't want user access)

    SerialQRDenseSolver(const SerialQRDenseSolver<OrdinalType, ScalarType>& Source);
    SerialQRDenseSolver & operator=(const SerialQRDenseSolver<OrdinalType, ScalarType>& Source);

  };

  namespace details {

    // Helper traits to distinguish work arrays for real and complex-valued datatypes.
    template<typename T>
    struct lapack_traits {
      typedef int iwork_type;
    };

    // Complex-valued specialization
    template<typename T>
    struct lapack_traits<std::complex<T> > {
      typedef typename ScalarTraits<T>::magnitudeType iwork_type;
    };

  } // end namespace details

//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialQRDenseSolver<OrdinalType,ScalarType>::SerialQRDenseSolver()
  : CompObject(),
    Object("Teuchos::SerialQRDenseSolver"),
    equilibrate_(false),
    shouldEquilibrate_(false),
    equilibratedA_(false),
    equilibratedB_(false),
    equilibrationModeA_(0),
    equilibrationModeB_(0),
    transpose_(false),
    factored_(false),
    solved_(false),
    matrixIsZero_(false),
    formedQ_(false),
    formedR_(false),
    TRANS_(Teuchos::NO_TRANS),
    M_(0),
    N_(0),
    Min_MN_(0),
    LDA_(0),
    LDAF_(0),
    LDQ_(0),
    LDR_(0),
    INFO_(0),
    LWORK_(0),
    ANORM_(ScalarTraits<MagnitudeType>::zero()),
    BNORM_(ScalarTraits<MagnitudeType>::zero()),
    A_(0),
    AF_(0),
    Q_(0),
    R_(0)
{
  resetMatrix();
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
SerialQRDenseSolver<OrdinalType,ScalarType>::~SerialQRDenseSolver()
{}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialQRDenseSolver<OrdinalType,ScalarType>::resetVectors()
{
  LHS_ = Teuchos::null;
  TMP_ = Teuchos::null;
  RHS_ = Teuchos::null;
  solved_ = false;
  equilibratedB_ = false;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialQRDenseSolver<OrdinalType,ScalarType>::resetMatrix()
{
  resetVectors();
  equilibratedA_ = false;
  equilibrationModeA_ = 0;
  equilibrationModeB_ = 0;
  factored_ = false;
  matrixIsZero_ = false;
  formedQ_ = false;
  formedR_ = false;
  M_ = 0;
  N_ = 0;
  Min_MN_ = 0;
  LDA_ = 0;
  LDAF_ = 0;
  LDQ_ = 0;
  LDR_ = 0;
  ANORM_ = -ScalarTraits<MagnitudeType>::one();
  BNORM_ = -ScalarTraits<MagnitudeType>::one();
  A_ = 0;
  AF_ = 0;
  Q_ = 0;
  R_ = 0;
  INFO_ = 0;
  LWORK_ = 0;
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::setMatrix(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& A)
{
  TEUCHOS_TEST_FOR_EXCEPTION(A->numRows() < A->numCols(), std::invalid_argument,
		     "SerialQRDenseSolver<T>::setMatrix: the matrix A must have A.numRows() >= A.numCols()!");

  resetMatrix();
  Matrix_ = A;
  Factor_ = A;
  FactorQ_ = A;
  FactorR_ = A;
  M_ = A->numRows();
  N_ = A->numCols();
  Min_MN_ = TEUCHOS_MIN(M_,N_);
  LDA_ = A->stride();
  LDAF_ = LDA_;
  LDQ_ = LDA_;
  LDR_ = N_;
  A_ = A->values();
  AF_ = A->values();
  Q_ = A->values();
  R_ = A->values();

  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::setVectors(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& X,
							   const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& B)
{
  // Check that these new vectors are consistent
  TEUCHOS_TEST_FOR_EXCEPTION(B->numCols() != X->numCols(), std::invalid_argument,
		     "SerialQRDenseSolver<T>::setVectors: X and B have different numbers of columns!");
  TEUCHOS_TEST_FOR_EXCEPTION(B->values()==0, std::invalid_argument,
		     "SerialQRDenseSolver<T>::setVectors: B is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(X->values()==0, std::invalid_argument,
		     "SerialQRDenseSolver<T>::setVectors: X is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(B->stride()<1, std::invalid_argument,
		     "SerialQRDenseSolver<T>::setVectors: B has an invalid stride!");
  TEUCHOS_TEST_FOR_EXCEPTION(X->stride()<1, std::invalid_argument,
		     "SerialQRDenseSolver<T>::setVectors: X has an invalid stride!");

  resetVectors();
  LHS_ = X;
  RHS_ = B;

  return(0);
}
//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::factor() {

  if (factored()) return(0);

  // Equilibrate matrix if necessary
  int ierr = 0;
  if (equilibrate_) ierr = equilibrateMatrix();
  if (ierr!=0) return(ierr);

  allocateWORK();
  if ((int)TAU_.size()!=Min_MN_) TAU_.resize( Min_MN_ );

  INFO_ = 0;

  // Factor
#ifdef HAVE_TEUCHOS_EIGEN
  EigenMatrixMap aMap( AF_, M_, N_, EigenOuterStride(LDAF_) );
  EigenScalarArray tau;
  EigenScalarArrayMap tauMap(&TAU_[0],TEUCHOS_MIN(M_,N_));
  qr_.compute(aMap);
  tau = qr_.hCoeffs();
  for (OrdinalType i=0; i<tauMap.innerSize(); i++) {
    tauMap(i) = tau(i);
  }
  EigenMatrix qrMat = qr_.matrixQR();
  for (OrdinalType j=0; j<aMap.outerSize(); j++) {
    for (OrdinalType i=0; i<aMap.innerSize(); i++) {
      aMap(i,j) = qrMat(i,j);
    }
  }
#else
  this->GEQRF(M_, N_, AF_, LDAF_, &TAU_[0], &WORK_[0], LWORK_, &INFO_);
#endif

  factored_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::solve() {

  // Check if the matrix is zero
  if (matrixIsZero_) {
    LHS_->putScalar(ScalarTraits<ScalarType>::zero());
    return(0);
  }

  // Equilibrate RHS if necessary
  int ierr = 0;
  if (equilibrate_) {
    ierr = equilibrateRHS();
    equilibratedB_ = true;
  }
  if (ierr != 0) return(ierr);

  TEUCHOS_TEST_FOR_EXCEPTION( (equilibratedA_ && !equilibratedB_) || (!equilibratedA_ && equilibratedB_) ,
  		     std::logic_error, "SerialQRDenseSolver<T>::solve: Matrix and vectors must be similarly scaled!");
  TEUCHOS_TEST_FOR_EXCEPTION( RHS_==Teuchos::null, std::invalid_argument,
		     "SerialQRDenseSolver<T>::solve: No right-hand side vector (RHS) has been set for the linear system!");
  TEUCHOS_TEST_FOR_EXCEPTION( LHS_==Teuchos::null, std::invalid_argument,
		     "SerialQRDenseSolver<T>::solve: No solution vector (LHS) has been set for the linear system!");
  if ( TRANS_ == Teuchos::NO_TRANS ) {
    TEUCHOS_TEST_FOR_EXCEPTION( LHS_->numRows() != N_, std::invalid_argument,
   		     "SerialQRDenseSolver<T>::solve: No transpose: must have LHS_->numRows() = N_");
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( LHS_->numRows() != M_, std::invalid_argument,
   		     "SerialQRDenseSolver<T>::solve: Transpose: must have LHS_->numRows() = M_");
  }

  if (shouldEquilibrate() && !equilibratedA_)
    std::cout << "WARNING!  SerialQRDenseSolver<T>::solve: System should be equilibrated!" << std::endl;

  // Matrix must be factored
  if (!factored()) factor();

  TMP_ = rcp( new SerialDenseMatrix<OrdinalType,ScalarType>(M_, RHS_->numCols()) );
  for (OrdinalType j=0; j<RHS_->numCols(); j++) {
    for (OrdinalType i=0; i<RHS_->numRows(); i++) {
      (*TMP_)(i,j) = (*RHS_)(i,j);
    }
  }

  INFO_ = 0;

  // Solve
  if ( TRANS_ == Teuchos::NO_TRANS ) {

    // Solve A lhs = rhs
    this->multiplyQ( Teuchos::CONJ_TRANS, *TMP_ );
    this->solveR( Teuchos::NO_TRANS, *TMP_ );

  } else {

    // Solve A**H lhs = rhs
    this->solveR( Teuchos::CONJ_TRANS, *TMP_ );
    for (OrdinalType j = 0; j < RHS_->numCols(); j++ ) {
      for (OrdinalType i = N_; i < M_; i++ ) {
	(*TMP_)(i, j) = ScalarTraits<ScalarType>::zero();
      }
    }
    this->multiplyQ( Teuchos::NO_TRANS, *TMP_ );

  }
  for (OrdinalType j = 0; j < LHS_->numCols(); j++ ) {
    for (OrdinalType i = 0; i < LHS_->numRows(); i++ ) {
      (*LHS_)(i, j) = (*TMP_)(i,j);
    }
  }

  solved_ = true;

  // Unequilibrate LHS if necessary
  int ierr1=0;
  if (equilibrate_) ierr1 = unequilibrateLHS();
  if (ierr != 0) return(ierr);

  return(INFO_);

}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::computeEquilibrateScaling()
{
  MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
  MagnitudeType prec = ScalarTraits<ScalarType>::prec();
  ScalarType sZero = ScalarTraits<ScalarType>::zero();
  ScalarType sOne  = ScalarTraits<ScalarType>::one();
  MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);

  MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin/prec);
  MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);

  // Compute maximum absolute value of matrix entries
  OrdinalType i, j;
  MagnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  for (j = 0; j < N_; j++) {
    for (i = 0; i < M_; i++) {
      anorm = TEUCHOS_MAX( anorm, ScalarTraits<ScalarType>::magnitude((*Matrix_)(i,j)) );
    }
  }

  ANORM_ = anorm;

  if (ANORM_ > mZero && ANORM_ < smlnum) {
    // Scale matrix norm up to smlnum
    shouldEquilibrate_ = true;
  } else if (ANORM_ > bignum) {
    // Scale matrix norm down to bignum
    shouldEquilibrate_ = true;
  } else if (ANORM_ == mZero) {
    // Matrix all zero. Return zero solution
    matrixIsZero_ = true;
  }

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::equilibrateMatrix()
{
  if (equilibratedA_) return(0);

  MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
  MagnitudeType prec = ScalarTraits<ScalarType>::prec();
  ScalarType sZero = ScalarTraits<ScalarType>::zero();
  ScalarType sOne  = ScalarTraits<ScalarType>::one();
  MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);

  MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin/prec);
  MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);
  OrdinalType BW = 0;

  // Compute maximum absolute value of matrix entries
  OrdinalType i, j;
  MagnitudeType anorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  for (j = 0; j < N_; j++) {
    for (i = 0; i < M_; i++) {
      anorm = TEUCHOS_MAX( anorm, ScalarTraits<ScalarType>::magnitude((*Matrix_)(i,j)) );
    }
  }

  ANORM_ = anorm;
  int ierr1 = 0;
  if (ANORM_ > mZero && ANORM_ < smlnum) {
    // Scale matrix norm up to smlnum
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, ANORM_, smlnum, M_, N_, A_, LDA_, &INFO_);
    equilibrationModeA_ = 1;
  } else if (ANORM_ > bignum) {
    // Scale matrix norm down to bignum
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, ANORM_, bignum, M_, N_, A_, LDA_, &INFO_);
    equilibrationModeA_ = 2;
  } else if (ANORM_ == mZero) {
    // Matrix all zero. Return zero solution
    matrixIsZero_ = true;
  }

  equilibratedA_ = true;

  return(ierr1);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::equilibrateRHS()
{
  if (equilibratedB_) return(0);

  MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
  MagnitudeType prec = ScalarTraits<ScalarType>::prec();
  ScalarType sZero = ScalarTraits<ScalarType>::zero();
  ScalarType sOne  = ScalarTraits<ScalarType>::one();
  MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);

  MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin/prec);
  MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);
  OrdinalType BW = 0;

  // Compute maximum absolute value of rhs entries
  OrdinalType i, j;
  MagnitudeType bnorm = ScalarTraits<ScalarType>::magnitude(ScalarTraits<ScalarType>::zero());
  for (j = 0; j <RHS_->numCols(); j++) {
    for (i = 0; i < RHS_->numRows(); i++) {
      bnorm = TEUCHOS_MAX( bnorm, ScalarTraits<ScalarType>::magnitude((*RHS_)(i,j)) );
    }
  }

  BNORM_ = bnorm;

  int ierr1 = 0;
  if (BNORM_ > mZero && BNORM_ < smlnum) {
    // Scale RHS norm up to smlnum
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, BNORM_, smlnum, RHS_->numRows(), RHS_->numCols(),
     		RHS_->values(), RHS_->stride(), &INFO_);
    equilibrationModeB_ = 1;
  } else if (BNORM_ > bignum) {
    // Scale RHS norm down to bignum
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, BNORM_, bignum, RHS_->numRows(), RHS_->numCols(),
     		RHS_->values(), RHS_->stride(), &INFO_);
    equilibrationModeB_ = 2;
  }

  equilibratedB_ = true;

  return(ierr1);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::unequilibrateLHS()
{
  if (!equilibratedB_) return(0);

  MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
  MagnitudeType prec = ScalarTraits<ScalarType>::prec();
  ScalarType sZero = ScalarTraits<ScalarType>::zero();
  ScalarType sOne  = ScalarTraits<ScalarType>::one();
  MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);

  MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin/prec);
  MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);
  OrdinalType BW = 0;

  int ierr1 = 0;
  if (equilibrationModeA_ == 1) {
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, ANORM_, smlnum, LHS_->numRows(), LHS_->numCols(),
		LHS_->values(), LHS_->stride(), &INFO_);
  } else if (equilibrationModeA_ == 2) {
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, ANORM_, bignum, LHS_->numRows(), LHS_->numCols(),
		LHS_->values(), LHS_->stride(), &INFO_);
  }
  if (equilibrationModeB_ == 1) {
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, smlnum, BNORM_, LHS_->numRows(), LHS_->numCols(),
		LHS_->values(), LHS_->stride(), &INFO_);
  } else if (equilibrationModeB_ == 2) {
    this->LASCL(Teuchos::ETypeChar[Teuchos::FULL], BW, BW, bignum, BNORM_, LHS_->numRows(), LHS_->numCols(),
		LHS_->values(), LHS_->stride(), &INFO_);
  }

  return(ierr1);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::formQ() {

  // Matrix must be factored first
  if (!factored()) factor();

  // Store Q separately from factored matrix
  if (AF_ == Q_) {
    FactorQ_ = rcp( new SerialDenseMatrix<OrdinalType,ScalarType>(*Factor_) );
    Q_ = FactorQ_->values();
    LDQ_ = FactorQ_->stride();
  }

  INFO_ = 0;

  // Form Q
#ifdef HAVE_TEUCHOS_EIGEN
  EigenMatrixMap qMap( Q_, M_, N_, EigenOuterStride(LDQ_) );
  EigenMatrix qMat = qr_.householderQ();
  for (OrdinalType j=0; j<qMap.outerSize(); j++) {
    for (OrdinalType i=0; i<qMap.innerSize(); i++) {
      qMap(i,j) = qMat(i,j);
    }
  }
#else
  this->UNGQR(M_, N_, N_, Q_, LDQ_, &TAU_[0], &WORK_[0], LWORK_, &INFO_);
#endif

  if (INFO_!=0) return(INFO_);

  formedQ_ = true;

  return(INFO_);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int SerialQRDenseSolver<OrdinalType,ScalarType>::formR() {

  // Matrix must be factored first
  if (!factored()) factor();

  // Store R separately from factored matrix
  if (AF_ == R_) {
    FactorR_ = rcp( new SerialDenseMatrix<OrdinalType,ScalarType>(N_, N_) );
    R_ = FactorR_->values();
    LDR_ = FactorR_->stride();
  }

  // Form R
  for (OrdinalType j=0; j<N_; j++) {
    for (OrdinalType i=0; i<=j; i++) {
      (*FactorR_)(i,j) = (*Factor_)(i,j);
    }
  }

  formedR_ = true;

  return(0);
}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int  SerialQRDenseSolver<OrdinalType, ScalarType>::multiplyQ(ETransp transq, SerialDenseMatrix<OrdinalType, ScalarType> &C)
{

  // Check that the matrices are consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(C.numRows()!=M_, std::invalid_argument,
		     "SerialQRDenseSolver<T>::multiplyQ: C.numRows() != M_!");
  TEUCHOS_TEST_FOR_EXCEPTION(C.values()==0, std::invalid_argument,
		     "SerialQRDenseSolver<T>::multiplyQ: C is an empty SerialDenseMatrix<T>!");

  // Matrix must be factored
  if (!factored()) factor();

  INFO_ = 0;

  // Multiply
#ifdef HAVE_TEUCHOS_EIGEN
  EigenMatrixMap cMap( C.values(), C.numRows(), C.numCols(), EigenOuterStride(C.stride()) );
  if ( transq == Teuchos::NO_TRANS ) {
    // C = Q * C
    cMap = qr_.householderQ() * cMap;
  } else {
    // C = Q**H * C
    cMap = qr_.householderQ().adjoint() * cMap;
    for (OrdinalType j = 0; j < C.numCols(); j++) {
      for (OrdinalType i = N_; i < C.numRows(); i++ ) {
	cMap(i, j) = ScalarTraits<ScalarType>::zero();
      }
    }
  }
#else
  Teuchos::ETransp NO_TRANS = Teuchos::NO_TRANS;
  Teuchos::ETransp TRANS = (Teuchos::ScalarTraits<ScalarType>::isComplex) ? Teuchos::CONJ_TRANS : Teuchos::TRANS;
  Teuchos::ESide SIDE = Teuchos::LEFT_SIDE;

  if ( transq == Teuchos::NO_TRANS ) {

    // C = Q * C
    this->UNMQR(ESideChar[SIDE], ETranspChar[NO_TRANS], M_, C.numCols(), N_, AF_, LDAF_,
		&TAU_[0], C.values(), C.stride(), &WORK_[0], LWORK_, &INFO_);

  } else {

    // C = Q**H * C
    this->UNMQR(ESideChar[SIDE], ETranspChar[TRANS], M_, C.numCols(), N_, AF_, LDAF_,
		&TAU_[0], C.values(), C.stride(), &WORK_[0], LWORK_, &INFO_);

    for (OrdinalType j = 0; j < C.numCols(); j++) {
      for (OrdinalType i = N_; i < C.numRows(); i++ ) {
	C(i, j) = ScalarTraits<ScalarType>::zero();
      }
    }
  }
#endif

  return(INFO_);

}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
int  SerialQRDenseSolver<OrdinalType, ScalarType>::solveR(ETransp transr, SerialDenseMatrix<OrdinalType, ScalarType> &C)
{

  // Check that the matrices are consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(C.numRows()<N_ || C.numRows()>M_, std::invalid_argument,
		     "SerialQRDenseSolver<T>::solveR: must have N_ < C.numRows() < M_!");
  TEUCHOS_TEST_FOR_EXCEPTION(C.values()==0, std::invalid_argument,
		     "SerialQRDenseSolver<T>::solveR: C is an empty SerialDenseMatrix<T>!");

  // Matrix must be factored
  if (!factored()) factor();

  INFO_ = 0;

  // Solve
#ifdef HAVE_TEUCHOS_EIGEN
  EigenMatrixMap cMap( C.values(), N_, C.numCols(), EigenOuterStride(C.stride()) );
  // Check for singularity first like TRTRS
  for (OrdinalType j=0; j<N_; j++) {
    if ((qr_.matrixQR())(j,j) == ScalarTraits<ScalarType>::zero()) {
      INFO_ = j+1;
      return(INFO_);
    }
  }
  if ( transr == Teuchos::NO_TRANS ) {
    // C = R \ C
    qr_.matrixQR().topRows(N_).template triangularView<Eigen::Upper>().solveInPlace(cMap);
  } else {
    // C = R**H \ C
    qr_.matrixQR().topRows(N_).template triangularView<Eigen::Upper>().adjoint().solveInPlace(cMap);
  }
#else
  Teuchos::ETransp NO_TRANS = Teuchos::NO_TRANS;
  Teuchos::ETransp TRANS = (Teuchos::ScalarTraits<ScalarType>::isComplex) ? Teuchos::CONJ_TRANS : Teuchos::TRANS;
  Teuchos::EUplo UPLO = Teuchos::UPPER_TRI;
  Teuchos::EDiag DIAG = Teuchos::NON_UNIT_DIAG;

  ScalarType* RMAT = (formedR_) ? R_ : AF_;
  OrdinalType LDRMAT = (formedR_) ? LDR_ : LDAF_;

  if ( transr == Teuchos::NO_TRANS ) {

    // C = R \ C
    this->TRTRS(EUploChar[UPLO], ETranspChar[NO_TRANS], EDiagChar[DIAG], N_, C.numCols(),
		RMAT, LDRMAT, C.values(), C.stride(), &INFO_);

  } else {

    // C = R**H \ C
    this->TRTRS(EUploChar[UPLO], ETranspChar[TRANS], EDiagChar[DIAG], N_, C.numCols(),
		RMAT, LDRMAT, C.values(), C.stride(), &INFO_);

  }
#endif

  return(INFO_);

}

//=============================================================================

template<typename OrdinalType, typename ScalarType>
void SerialQRDenseSolver<OrdinalType,ScalarType>::Print(std::ostream& os) const {

  if (Matrix_!=Teuchos::null) os << "Solver Matrix"          << std::endl << *Matrix_ << std::endl;
  if (Factor_!=Teuchos::null) os << "Solver Factored Matrix" << std::endl << *Factor_ << std::endl;
  if (Q_!=Teuchos::null) os << "Solver Factor Q" << std::endl << *Q_ << std::endl;
  if (LHS_   !=Teuchos::null) os << "Solver LHS"             << std::endl << *LHS_    << std::endl;
  if (RHS_   !=Teuchos::null) os << "Solver RHS"             << std::endl << *RHS_    << std::endl;

}

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALQRDENSESOLVER_HPP_ */

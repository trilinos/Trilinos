// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
#ifndef BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP
#define BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP

/*! \file BelosTeuchosDenseAdapter.hpp
  \brief Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix
  with ordinal type int and arbitrary scalar type.
*/

#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"

#include "BelosDenseMatTraits.hpp" 
#include "BelosDenseSolver.hpp"

namespace Belos {

  //! Full specialization of Belos::DenseSolver for Teuchos::SerialDenseMatrix<int,ST>.
  template<class ScalarType>
  class TeuchosDenseSolver : public DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType> >
  {
  public:

    //! @name Constructor/Destructor Methods
    //@{
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    TeuchosDenseSolver() {}

    //! TeuchosDenseSolver destructor.
    virtual ~TeuchosDenseSolver() {}
    //@}

    //! @name Set Methods
    //@{
    //! Sets the pointers for coefficient matrix
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::setMatrix;

    //! Sets the pointers for left and right hand side vector(s).
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::setVectors;
    //@}

    //! @name Strategy Modifying Methods
    //@{

    //! Set if dense matrix is symmetric positive definite
    //! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::setSPD;

    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    */
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::factorWithEquilibration;

    //! All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS, \c Teuchos::TRANS, and \c Teuchos::CONJ_TRANS).
    /*! \note This interface will set correct transpose flag for matrix, including complex-valued linear systems.
    */
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::solveWithTransposeFlag;
    //@}

    //! @name Factor/Solve Methods
    //@{

    //! Computes the in-place LU factorization of the matrix.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor() 
    {
       //std::cout << "Calling DenseMatrix<Teuchos> factor!" << std::endl;
       // Create solver object if one doesn't exist.
       if ( (dSolver_==Teuchos::null) && (spdSolver_==Teuchos::null) )
       {
         if (spd_)
           spdSolver_ = Teuchos::rcp( new Teuchos::SerialSpdDenseSolver<int,ScalarType>() );
         else
           dSolver_ = Teuchos::rcp( new Teuchos::SerialDenseSolver<int,ScalarType>() );
      }

      int info = 0;

      // Set the matrix and factor
      if (spd_)
      {
        if (newMatrix_)
        {
          spdMatrix_ = Teuchos::rcp( new Teuchos::SerialSymDenseMatrix<int,ScalarType>(Teuchos::View, true, A_->values(),
                                                                                       A_->numRows(), A_->numCols()) );
        } 
        spdSolver_->setMatrix( spdMatrix_ );
        spdSolver_->factorWithEquilibration( equilibrate_ );
        info = spdSolver_->factor();
      }
      else
      {
        dSolver_->setMatrix( A_ );
        dSolver_->solveWithTransposeFlag( TRANS_ );
        dSolver_->factorWithEquilibration( equilibrate_ );
        info = dSolver_->factor();
      }     

      if (!info) 
        newMatrix_ = false;
      return info;
    }

    //! Computes the solution X to AX = B for the \e this matrix and the B provided.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve() 
    {
      //std::cout << "Calling DenseMatrix<Teuchos> solve!" << std::endl;
      // Check if this is a new matrix that has not been factored.
      if (newMatrix_)
        factor();

      int info = 0;

      if (spd_)
      {
        spdSolver_->setVectors( X_, B_ );
        info = spdSolver_->solve();
      }
      else
      {
        dSolver_->setVectors( X_, B_ );
        info = dSolver_->solve();
      }

      return info;
    }

    //@}

  private:

    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::newMatrix_;
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::equilibrate_;
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::TRANS_;
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::spd_;

    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::A_;
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::X_;
    using DenseSolver<ScalarType,Teuchos::SerialDenseMatrix<int,ScalarType>>::B_;

    // Pointer to dense solver for general linear systems
    Teuchos::RCP<Teuchos::SerialDenseSolver<int,ScalarType>> dSolver_;

    // Pointer to dense solver for SPD linear systems
    Teuchos::RCP<Teuchos::SerialSymDenseMatrix<int,ScalarType>> spdMatrix_;
    Teuchos::RCP<Teuchos::SerialSpdDenseSolver<int,ScalarType>> spdSolver_;

  };


  //! Full specialization of Belos::DenseMatTraits for Teuchos::SerialDenseMatrix<int,ST>.
  template<class ScalarType>
  class DenseMatTraits<ScalarType, Teuchos::SerialDenseMatrix<int,ScalarType>>{
  public:

    using ST = Teuchos::ScalarTraits<ScalarType>;
    using MagnitudeType = typename ST::magnitudeType;
    using DM = Teuchos::SerialDenseMatrix<int,ScalarType>;
    using MDM = Teuchos::SerialDenseMatrix<int, MagnitudeType>;
    using MDMT = DenseMatTraits<MagnitudeType, MDM>;
    
    //@{ \name Creation methods

    /*! \brief Creates a new empty \c DM with no dimension.

    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> Create() {
      return Teuchos::rcp(new DM());
    }     

    /*! \brief Creates a new empty \c DM containing \c numvecs columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> Create( const int numrows, const int numcols, bool initZero = true) {
      return Teuchos::rcp(new DM(numrows,numcols,initZero));
    }     

    /*! \brief Create a new copy \c DM, possibly transposed.
   
       \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM> CreateCopy(const DM & dm, bool transpose=false) {
      if (transpose)    
        return Teuchos::rcp(new DM(dm, Teuchos::CONJ_TRANS));
      else 
        return Teuchos::rcp(new DM(Teuchos::Copy, dm));
    }

    //! \brief Returns a raw pointer to the (non-const) data on the host.
    static ScalarType* GetRawHostPtr( DM & dm )
    { return dm.values(); }     

    //! \brief Returns a raw pointer to const data on the host.
    static ScalarType const * GetConstRawHostPtr(const DM & dm )
    { return const_cast<ScalarType const *>(dm.values()); }     

    //! \brief Marks host data modified to avoid device sync errors. 
    /// \note Belos developers must call this function after EVERY
    ///   call to LAPACK!!!
    //static void RawPtrDataModified(DM & dm ) { }

    //! \brief Returns an RCP to a DM which has a subview of the given DM.
    //        Row and column indexing is zero-based.
    //        Source  - Reference to another dense matrix from which values are to be copied.
    //        numRows - The number of rows in this matrix.
    //        numCols - The number of columns in this matrix.
    //        startRow  - The row of Source from which the submatrix copy should start.
    //        startCol  - The column of Source from which the submatrix copy should start.
    //  
    //        Should ints be const? Should they be ints or some other ordinal type?
    static Teuchos::RCP<DM> Subview(DM & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { return Teuchos::rcp(new DM(Teuchos::View, source, numRows, numCols, startRow, startCol)); }

    static Teuchos::RCP<const DM> SubviewConst(
                              const DM & source, int numRows, int numCols, int startRow=0, int startCol=0)
    { return Teuchos::rcp(new const DM(Teuchos::View, source, numRows, numCols, startRow, startCol)); }

    //! \brief Returns a deep copy of the requested subview.
    static Teuchos::RCP<DM> SubviewCopy( const DM& source, int numRows, int numCols, int startRow=0, int startCol=0)
    { return Teuchos::rcp(new DM(Teuchos::Copy, source, numRows, numCols, startRow, startCol)); }
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const DM& dm )
    { return dm.numRows(); }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const DM& dm )
    { return dm.numCols(); }     

    //! \brief Obtain the stride between the columns of \c dm.
    static int GetStride( const DM& dm )
    { return dm.stride(); }    

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm to have \c numrows rows and \c numcols columns.
     *        All values will be initialized to zero if the final argument is true. 
     *        If the final argument is fale, the previous entries in 
     *        the matrix will be maintained. For new entries that did not exist in the previous matrix, values will
     *        contain noise from memory. 
    */
    static void Reshape( DM& dm, const int numrows, const int numcols, bool initZero = false) {
      if(initZero){
        int err =  dm.shape(numrows,numcols); 
        if(err != 0){throw std::runtime_error ("Error in DenseMatrixTraits::shape. Teuchos::SerialDenseMatrix.shape failed.");}
      }
      else{
        int err =  dm.reshape(numrows,numcols); 
        if(err != 0){throw std::runtime_error ("Error in DenseMatrixTraits::reshape. Teuchos::SerialDenseMatrix.reshape failed.");}
      }
    }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static ScalarType & Value( DM& dm, const int i, const int j )
    { 
      return dm(i,j);
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const ScalarType & ValueConst( const DM& dm, const int i, const int j )
    { 
      return dm(i,j);
    }

    static void SyncDeviceToHost(DM &){ }

    static void SyncHostToDevice(DM &){ }
    //@}

    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( DM& thisDM, const DM& sourceDM){
      thisDM += sourceDM;
    }

    //!  \brief Add source to diagonal entries of dest
    static void AddDiag( DM& dest, ScalarType source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumCols(dest));
      for (int i = 0; i < GetNumRows(dest); i++) {
          Value(dest, i, i) += source;
      }
    }

    //!  \brief Fill all entries with \c value. Value is zero if not specified.
    static void PutScalar( DM & dm, ScalarType value = Teuchos::ScalarTraits<ScalarType>::zero()){
      dm.putScalar(value);
    }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( DM& dm, ScalarType value){
      dm.scale(value);
    }

    //!  \brief Multiply two dense matrices. C = beta*C + alpha*op(A)*op(B)
    static void Multiply( bool transposeA, bool transposeB, ScalarType alpha, const Teuchos::SerialDenseMatrix<int,ScalarType>& A, const Teuchos::SerialDenseMatrix<int,ScalarType> &B, ScalarType beta, Teuchos::SerialDenseMatrix<int,ScalarType>& C)
    {
      Teuchos::BLAS<int, ScalarType> blas;

      auto transATeuchos = transposeA ? Teuchos::TRANS : Teuchos::NO_TRANS;
      auto transBTeuchos = transposeB ? Teuchos::TRANS : Teuchos::NO_TRANS;
      blas.GEMM( transATeuchos, transBTeuchos,
                transposeA ? GetNumCols(A) : GetNumRows(A),
                transposeB ? GetNumRows(B) : GetNumCols(B),
                transposeA ? GetNumRows(A) : GetNumCols(A),
                alpha,
                GetConstRawHostPtr(A), GetStride(A),
                GetConstRawHostPtr(B), GetStride(B),
                beta, GetRawHostPtr(C), GetStride(C));
    }

    //!  \brief Fill the DM with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy). 
    //TODO What to do here? Kinda needs random synced version??
    static void Randomize( DM& dm){
      dm.random();
    }

    //!  \brief Copies entries of source to dest (deep copy). 
    static void Assign( DM& dest, const DM& source){
      dest.assign(source);
    }

    //!  \brief Assign source to diagonal entries of dest
    static void AssignDiag( DM& dest, ScalarType source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumCols(dest));
      for (int i = 0; i < GetNumRows(dest); i++) {
          Value(dest, i, i) = source;
      }
    }

    //!  \brief Assign source to diagonal entries of dest
    static void AssignDiag( DM& dest, const DM& source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumRows(source));
      TEUCHOS_ASSERT(GetNumCols(dest) == GetNumRows(source));
      TEUCHOS_ASSERT(GetNumCols(source) == 1);
      for (int i = 0; i < GetNumRows(source); i++) {
          Value(dest, i, i) = ValueConst(source, i, 0);
      }
    }

    //!  \brief Assign upper triangular entries of source to dest
    static void AssignUpperTri( DM& dest, const DM& source){
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumRows(source));
      TEUCHOS_ASSERT(GetNumCols(dest) == GetNumCols(source));
      for (int i = 0; i < GetNumRows(source); i++) {
        for (int j = i; j < GetNumCols(source); j++)
          Value(dest, i, j) = ValueConst(source, i, j);
      }
    }

    //!  \brief Returns the Frobenius norm of the dense matrix.
    static typename Teuchos::ScalarTraits<ScalarType>::magnitudeType NormFrobenius(const DM& dm) {
      return dm.normFrobenius(); 
    }

    //!  \brief Returns the one-norm of the dense matrix.
    static typename Teuchos::ScalarTraits<ScalarType>::magnitudeType NormOne(const DM& dm) {
      return dm.normOne();
    }
    //@}

    //@{ \name Solver methods 
    
    //!  \brief Returns a dense solver object for the dense matrix.
    static Teuchos::RCP<DenseSolver<ScalarType, DM>>
      createDenseSolver() { 
   
      Teuchos::RCP<DenseSolver<ScalarType, DM>> newSolver
        = Teuchos::rcp( new TeuchosDenseSolver<ScalarType>() );
      return newSolver;
    }
    //@}

  static void trsm(const char side[], const char uplo[], const char trans[], const char diag[],
                   const ScalarType& alpha, const DM& A, DM& B) {
    Teuchos::BLAS<int, ScalarType> blas;

    auto sideTeuchos = (*side == 'L') ? Teuchos::LEFT_SIDE : Teuchos::RIGHT_SIDE;
    auto uploTeuchos = (*uplo == 'U') ? Teuchos::UPPER_TRI : Teuchos::LOWER_TRI;
    auto transTeuchos = (*trans == 'N') ? Teuchos::NO_TRANS : Teuchos::TRANS;
    auto diagTeuchos = (*diag == 'N') ? Teuchos::NON_UNIT_DIAG : Teuchos::UNIT_DIAG;

    blas.TRSM( sideTeuchos, uploTeuchos, transTeuchos,
              diagTeuchos, GetNumRows(A), GetNumCols(B), alpha,
              GetConstRawHostPtr(A), GetStride(A),
              GetRawHostPtr(B), GetStride(B) );
  }

  static void nrm2(std::vector<MagnitudeType>&R, const DM &X) {
    Teuchos::BLAS<int, ScalarType> blas;
    for (int k = 0; k<GetNumCols(X); ++k)
      R[k] = blas.NRM2(GetNumRows(X), &ValueConst(X, 0, k), 1);
  }

  static void geqrf(DM &A, DM &tau) {
    Teuchos::LAPACK<int,ScalarType> lapack;

    // Step #1: Perform workspace size query for QR
    // factorization of HP.  On input, lwork must be -1.
    // _GEQRF will put the workspace size in work_[0].
    int lwork = -1;
    int info = 0;
    std::vector<ScalarType> work;
    work.resize(1);
    lapack.GEQRF (GetNumRows(A), GetNumCols(A), GetRawHostPtr(A),
                  GetStride(A), GetRawHostPtr(tau), &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, " LAPACK's _GEQRF failed to compute a workspace size.");

    // Step #2: Compute QR factorization of HP
    //
    // NOTE (mfh 17 Apr 2014) LAPACK promises that the value of
    // work_[0] after the workspace query will fit in int.  This
    // justifies the cast.  We call real() first because
    // static_cast from std::complex to int doesn't work.
    lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));

    work.resize (lwork); // Allocate workspace for the QR factorization
    lapack.GEQRF (GetNumRows(A), GetNumCols(A), GetRawHostPtr(A),
                  GetStride(A), GetRawHostPtr(tau), &work[0], lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, "LAPACK's _GEQRF failed to compute a QR factorization.");
  }

  static void ungqr(const int k, DM& A, const DM& tau) {
    Teuchos::LAPACK<int,ScalarType> lapack;

    int lwork = -1;
    int info = 0;
    std::vector<ScalarType> work;
    work.resize(1);
    lapack.UNGQR (GetNumRows(A), GetNumCols(A), k,
                  GetRawHostPtr(A), GetStride(A), GetConstRawHostPtr(tau), &work[0],
                  lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, "LAPACK's _UNGQR failed to construct the Q factor.");

    lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));

    work.resize (lwork);
    lapack.UNGQR (GetNumRows(A), GetNumCols(A), k,
                  GetRawHostPtr(A), GetStride(A), GetConstRawHostPtr(tau), &work[0],
                  lwork, &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, "LAPACK's _UNGQR failed to construct the Q factor.");
  }

  static void updateLSQR(DM &H, DM &z, Teuchos::RCP<MDM> cs, Teuchos::RCP<DM> sn, Teuchos::RCP<DM> beta, int dim, int blockSize) {
    Teuchos::BLAS<int, ScalarType> blas;
    const ScalarType zero = ST::zero();

    if (blockSize == 1) {
      for (int i=0; i<dim; i++) {
        //
        // Apply previous Givens rotations to new column of Hessenberg matrix
        //
        blas.ROT( 1, &Value(H,i,dim), 1, &Value(H,i+1, dim), 1, &MDMT::Value(*cs, i, 0), &Value(*sn, i, 0) );
      }

      // Calculate new Givens rotation
      blas.ROTG( &Value(H,dim,dim), &Value(H,dim+1,dim), &MDMT::Value(*cs, dim, 0), &Value(*sn, dim, 0) );
      Value(H,dim+1,dim) = zero;

      // Update RHS w/ new transformation
      blas.ROT( 1, &Value(z,dim,0), 1, &Value(z,dim+1,0), 1, &MDMT::Value(*cs, dim, 0), &Value(*sn, dim, 0) );
    } else {
      int maxidx;
      ScalarType sigma, mu, vscale, maxelem;
      //
      // QR factorization of Least-Squares system with Householder reflectors
      //
      for (int j=0; j<blockSize; j++) {
        //
        // Apply previous Householder reflectors to new block of Hessenberg matrix
        //
        for (int i=0; i<dim+j; i++) {
          sigma = blas.DOT( blockSize, &Value(H,i+1,i), 1, &Value(H,i+1,dim+j), 1);
          sigma += ValueConst(H,i,dim+j);
          sigma *= ST::conjugate(ValueConst(*beta, i, 0));
          blas.AXPY(blockSize, ScalarType(-sigma), &Value(H,i+1,i), 1, &Value(H,i+1,dim+j), 1);
          Value(H,i,dim+j) -= sigma;
        }
        //
        // Compute new Householder reflector
        //
        maxidx = blas.IAMAX( blockSize+1, &Value(H,dim+j,dim+j), 1 );
        maxelem = ST::magnitude(Value(H,dim+j+maxidx-1,dim+j));
        for (int i=0; i<blockSize+1; i++)
          Value(H,dim+j+i,dim+j) /= maxelem;
        sigma = blas.DOT( blockSize, &Value(H,dim+j+1,dim+j), 1,
                          &Value(H,dim+j+1,dim+j), 1 );
        MagnitudeType sign_Rjj = -ST::real(Value(H,dim+j,dim+j)) /
                 ST::magnitude(ST::real((Value(H,dim+j,dim+j))));
        if (sigma == zero) {
          Value(*beta, dim + j, 0) = zero;
        } else {
          mu = ST::squareroot(ST::conjugate(Value(H,dim+j,dim+j))*Value(H,dim+j,dim+j)+sigma);
          vscale = ValueConst(H,dim+j,dim+j) - Teuchos::as<ScalarType>(sign_Rjj)*mu;
          Value(*beta, dim + j, 0) = -Teuchos::as<ScalarType>(sign_Rjj) * vscale / mu;
          Value(H,dim+j,dim+j) = Teuchos::as<ScalarType>(sign_Rjj)*maxelem*mu;
          for (int i=0; i<blockSize; i++)
            Value(H,dim+j+1+i,dim+j) /= vscale;
        }
        //
        // Apply new Householder reflector to rhs
        //
        for (int i=0; i<blockSize; i++) {
          sigma = blas.DOT( blockSize, &Value(H,dim+j+1,dim+j),
                            1, &Value(z,dim+j+1,i), 1);
          sigma += ValueConst(z,dim+j,i);
          sigma *= ST::conjugate(ValueConst(*beta, dim+j, 0));
          blas.AXPY(blockSize, ScalarType(-sigma), &Value(H,dim+j+1,dim+j),
                    1, &Value(z,dim+j+1,i), 1);
          Value(z,dim+j,i) -= sigma;
        }
      }
    }
  }
};

} // namespace Belos

#endif // end file BELOS_TEUCHOS_DENSE_MAT_TRAITS_HPP

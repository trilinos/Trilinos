// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_SERIALQRDENSESOLVER_MP_VECTOR_HPP_
#define _TEUCHOS_SERIALQRDENSESOLVER_MP_VECTOR_HPP_
/*! \file Teuchos_SerialQRDenseSolver.hpp
  \brief Templated class for solving dense linear problems.
*/

#include <new>
#include "Teuchos_SerialQRDenseSolver.hpp"

#include "Sacado_MP_Vector.hpp"

/*! \class Teuchos::SerialQRDenseSolver
  \brief Specialization for Sacado::MP::Vector< Storage<...> >
*/

namespace Teuchos {

  template<typename OrdinalType, typename Storage>
  class SerialQRDenseSolver<OrdinalType, Sacado::MP::Vector<Storage> >
    : public CompObject,
      public Object
  {

  public:

    typedef Sacado::MP::Vector<Storage> ScalarType;
    typedef typename ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

    //! @name Constructor/Destructor Methods
    //@{
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    SerialQRDenseSolver();

    //! SerialQRDenseSolver destructor.
    virtual ~SerialQRDenseSolver() {}
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
    void factorWithEquilibration(bool flag) {
      base_QR_.factorWithEquilibration(flag);
    }

    //! If \c flag is true, causes all subsequent function calls to work with the adjoint of \e this matrix, otherwise not.
    void solveWithTranspose(bool flag) {
      base_QR_.solveWithTranspose(flag);
    }

    //! All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS or Teuchos::CONJ_TRANS).
    void solveWithTransposeFlag(Teuchos::ETransp trans) {
      base_QR_.solveWithTranspose(trans);
    }

    //@}

    //! @name Factor/Solve/Invert Methods
    //@{

    //! Computes the in-place QR factorization of the matrix using the LAPACK routine \e _GETRF or the Eigen class \e HouseholderQR
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor() { return base_QR_.factor(); }

    //! Computes the solution X to AX = B for the \e this matrix and the B provided to SetVectors()..
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve() { return base_QR_.solve(); }

    //! Determines if \e this matrix should be scaled.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int computeEquilibrateScaling() {
      return base_QR_.computeEquilibrateScaling();
    }

    //! Equilibrates the \e this matrix.
    /*!
      \note This method will be called automatically in solve() method if factorWithEquilibration( true ) is called.
      \return Integer error code, set to 0 if successful.
    */
    int equilibrateMatrix() { return base_QR_.equilibrateMatrix(); }

    //! Equilibrates the current RHS.
    /*!
      \note This method will be called automatically in solve() method if factorWithEquilibration( true ) is called.
      \return Integer error code, set to 0 if successful.
    */
    int equilibrateRHS() { return base_QR_.equilibrateRHS(); }

    //! Unscales the solution vectors if equilibration was used to solve the system.
    /*!
      \note This method will be called automatically in solve() method if factorWithEquilibration( true ) is called.
      \return Integer error code, set to 0 if successful.
    */
    int unequilibrateLHS() { return base_QR_.unequilibrateLHS(); }

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
    bool transpose() { return base_QR_.transpose(); }

    //! Returns true if matrix is factored (factor available via getFactoredMatrix()).
    bool factored() { return base_QR_.factored(); }

    //! Returns true if factor is equilibrated (factor available via getFactoredMatrix()).
    bool equilibratedA() { return base_QR_.equilibratedA(); }

    //! Returns true if RHS is equilibrated (RHS available via getRHS()).
    bool equilibratedB() { return base_QR_.equilibratedB(); }

    //! Returns true if the LAPACK general rules for equilibration suggest you should equilibrate the system.
    bool shouldEquilibrate() { return base_QR_.shouldEquilibrate(); }

    //! Returns true if the current set of vectors has been solved.
    bool solved() { return base_QR_.solved(); }

    //! Returns true if Q has been formed explicitly.
    bool formedQ() { return base_QR_.formedQ(); }

    //! Returns true if R has been formed explicitly.
    bool formedR() { return base_QR_.formedR();}

    //@}

    //! @name Data Accessor methods
    //@{

    //! Returns pointer to current matrix.
    RCP<SerialDenseMatrix<OrdinalType, ScalarType> > getMatrix() const {return(Matrix_);}

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
    std::vector<ScalarType> tau()  const { return base_QR_.tau(); }

    //! Returns the absolute value of the largest element of \e this matrix (returns -1 if not yet computed).
    MagnitudeType ANORM()  const { return base_QR_.ANORM(); }

    //@}

    //! @name I/O methods
    //@{
    //! Print service methods; defines behavior of ostream << operator.
    void Print(std::ostream& os) const;
    //@}
  protected:

    typedef typename ScalarType::value_type BaseScalarType;
    typedef SerialQRDenseSolver<OrdinalType,BaseScalarType> BaseQRType;
    typedef SerialDenseMatrix<OrdinalType, BaseScalarType> BaseMatrixType;
    typedef SerialDenseMatrix<OrdinalType, ScalarType> MatrixType;

    void resetMatrix();
    void resetVectors();
    RCP<BaseMatrixType> createBaseMatrix(const RCP<MatrixType>& mat) const;
    RCP<MatrixType> createMatrix(const RCP<BaseMatrixType>& base_mat) const;
    BaseQRType base_QR_;

    OrdinalType M_;
    OrdinalType N_;
    OrdinalType SacadoSize_;

    RCP<MatrixType> Matrix_;
    RCP<MatrixType> LHS_;
    RCP<MatrixType> RHS_;
    RCP<MatrixType> Factor_;
    RCP<MatrixType> FactorQ_;
    RCP<MatrixType> FactorR_;

    RCP<BaseMatrixType> Base_Matrix_;
    RCP<BaseMatrixType> Base_LHS_;
    RCP<BaseMatrixType> Base_RHS_;
    RCP<BaseMatrixType> Base_Factor_;
    RCP<BaseMatrixType> Base_FactorQ_;
    RCP<BaseMatrixType> Base_FactorR_;

  private:

    // SerialQRDenseSolver copy constructor (put here because we don't want user access)
    SerialQRDenseSolver(const SerialQRDenseSolver& Source);
    SerialQRDenseSolver & operator=(const SerialQRDenseSolver& Source);

  };

  namespace details {

    // Helper traits class for converting between arrays of
    // Sacado::MP::Vector and its corresponding value type
    template <typename Storage> struct MPVectorArrayHelper;

    template <typename Ordinal, typename Value, typename Device>
    struct MPVectorArrayHelper< Stokhos::DynamicStorage<Ordinal,Value,Device> >
    {
      typedef Stokhos::DynamicStorage<Ordinal,Value,Device> Storage;
      typedef Sacado::MP::Vector<Storage> Vector;

      static Value* getValueArray(Vector* v, Ordinal len) {
        if (len == 0)
          return 0;
        return v[0].coeff(); // Assume data is contiguous
      }

      static Vector* getVectorArray(Value *v, Ordinal len, Ordinal vector_size){
        if (len == 0)
          return 0;
        Vector* v_vec = static_cast<Vector*>(operator new(len*sizeof(Vector)));
        for (Ordinal i=0; i<len; ++i)
          new (v_vec+i) Vector(vector_size, v+i*vector_size, false);
        return v_vec;
      }

    };

    template <typename Ordinal, typename Value, typename Device, int Num>
    struct MPVectorArrayHelper< Stokhos::StaticFixedStorage<Ordinal,Value,Num,Device> > {
      typedef Stokhos::StaticFixedStorage<Ordinal,Value,Num,Device> Storage;
      typedef Sacado::MP::Vector<Storage> Vector;

      static Value* getValueArray(Vector* v, Ordinal len)
      {
        return reinterpret_cast<Value*>(v);
      }

      static Vector* getVectorArray(Value *v, Ordinal len, Ordinal vector_size)
      {
        return reinterpret_cast<Vector*>(v);
      }
    };

  }

  // Helper traits to distinguish work arrays for real and complex-valued datatypes.
  using namespace details;

//=============================================================================

template<typename OrdinalType, typename Storage>
SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
SerialQRDenseSolver()
  : CompObject(),
    Object("Teuchos::SerialQRDenseSolver"),
    M_(0),
    N_(0)
{
  resetMatrix();
}

//=============================================================================

template<typename OrdinalType, typename Storage>
void SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
resetVectors()
{
  LHS_ = Teuchos::null;
  RHS_ = Teuchos::null;
}
//=============================================================================

template<typename OrdinalType, typename Storage>
void SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
resetMatrix()
{
  resetVectors();
  M_ = 0;
  N_ = 0;
}
//=============================================================================

template<typename OrdinalType, typename Storage>
RCP<SerialDenseMatrix<OrdinalType, typename Sacado::MP::Vector<Storage>::value_type> >
SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
createBaseMatrix(
  const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& mat) const
{
  BaseScalarType* base_ptr =
    MPVectorArrayHelper<Storage>::getValueArray(
      mat->values(), mat->stride()*mat->numCols());
  RCP<BaseMatrixType> base_mat =
    rcp(new BaseMatrixType(Teuchos::View,
                           base_ptr,
                           mat->stride()*SacadoSize_,
                           mat->numRows()*SacadoSize_,
                           mat->numCols()));
  return base_mat;
}
//=============================================================================

template<typename OrdinalType, typename Storage>
RCP<SerialDenseMatrix<OrdinalType, Sacado::MP::Vector<Storage> > >
SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
createMatrix(
  const RCP<SerialDenseMatrix<OrdinalType,BaseScalarType> >& base_mat) const
{
  ScalarType* ptr =
    MPVectorArrayHelper<Storage>::getVectorArray(
      base_mat->values(), base_mat->stride()*base_mat->numCols(), SacadoSize_);
  RCP<MatrixType> mat =
    rcp(new MatrixType(Teuchos::View,
                       ptr,
                       base_mat->stride()/SacadoSize_,
                       base_mat->numRows()/SacadoSize_,
                       base_mat->numCols()));
  return mat;
}
//=============================================================================

template<typename OrdinalType, typename Storage>
int SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
setMatrix(const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& A)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A->numRows() < A->numCols(), std::invalid_argument,
    "SerialQRDenseSolver<T>::setMatrix: the matrix A must have A.numRows() >= A.numCols()!");

  resetMatrix();
  Matrix_ = A;
  Factor_ = A;
  FactorQ_ = A;
  FactorR_ = A;
  M_ = A->numRows();
  N_ = A->numCols();
  if (Storage::is_static)
    SacadoSize_ = Storage::static_size;
  else if (M_ > 0 && N_ > 0)
    SacadoSize_ = (*A)(0,0).size();
  else
    SacadoSize_ = 0;

  return base_QR_.setMatrix( createBaseMatrix(A) );
}
//=============================================================================

template<typename OrdinalType, typename Storage>
int SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
setVectors(
  const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& X,
  const RCP<SerialDenseMatrix<OrdinalType,ScalarType> >& B)
{
  // Check that these new vectors are consistent
  TEUCHOS_TEST_FOR_EXCEPTION(
    B->numCols() != X->numCols(), std::invalid_argument,
    "SerialQRDenseSolver<T>::setVectors: X and B have different numbers of columns!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    B->values()==0, std::invalid_argument,
    "SerialQRDenseSolver<T>::setVectors: B is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X->values()==0, std::invalid_argument,
    "SerialQRDenseSolver<T>::setVectors: X is an empty SerialDenseMatrix<T>!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    B->stride()<1, std::invalid_argument,
    "SerialQRDenseSolver<T>::setVectors: B has an invalid stride!");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X->stride()<1, std::invalid_argument,
    "SerialQRDenseSolver<T>::setVectors: X has an invalid stride!");

  resetVectors();
  LHS_ = X;
  RHS_ = B;

  return base_QR_.setVectors( createBaseMatrix(X), createBaseMatrix(B) );
}

//=============================================================================

template<typename OrdinalType, typename Storage>
int SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
formQ() {
  int ret = base_QR_.formQ();
  FactorQ_ = createMatrix( base_QR_.getQ() );
  return(ret);
}

//=============================================================================

template<typename OrdinalType, typename Storage>
int SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
formR() {
  int ret = base_QR_.formR();
  FactorR_ = createMatrix( base_QR_.getR() );
  Factor_ = createMatrix( base_QR_.getFactoredMatrix() );
  return(ret);
}

//=============================================================================

template<typename OrdinalType, typename Storage>
int  SerialQRDenseSolver<OrdinalType, Sacado::MP::Vector<Storage> >::
multiplyQ(ETransp transq, SerialDenseMatrix<OrdinalType, ScalarType> &C)
{
  return base_QR_.multiplyQ(transq, createBaseMatrix(C));
}

//=============================================================================

template<typename OrdinalType, typename Storage>
int  SerialQRDenseSolver<OrdinalType, Sacado::MP::Vector<Storage> >::
solveR(ETransp transr, SerialDenseMatrix<OrdinalType, ScalarType> &C)
{
  return base_QR_.solveR(transr, createBaseMatrix(C));
}

//=============================================================================

template<typename OrdinalType, typename Storage>
void SerialQRDenseSolver<OrdinalType,Sacado::MP::Vector<Storage> >::
Print(std::ostream& os) const {

  if (Matrix_!=Teuchos::null) os << "Solver Matrix"          << std::endl << *Matrix_ << std::endl;
  if (Factor_!=Teuchos::null) os << "Solver Factored Matrix" << std::endl << *Factor_ << std::endl;
  if (LHS_   !=Teuchos::null) os << "Solver LHS"             << std::endl << *LHS_    << std::endl;
  if (RHS_   !=Teuchos::null) os << "Solver RHS"             << std::endl << *RHS_    << std::endl;

}

} // namespace Teuchos

#endif /* _TEUCHOS_SERIALQRDENSESOLVER_MP_VECTOR_HPP_ */

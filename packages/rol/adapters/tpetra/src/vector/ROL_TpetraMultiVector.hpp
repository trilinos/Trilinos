// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TPETRAMULTIVECTOR_HPP
#define ROL_TPETRAMULTIVECTOR_HPP

/** \class ROL::TpetraMultiVector
    \brief Implements the ROL::Vector interface for a Tpetra_MultiVector.
    \author Created by Greg von Winckel
*/



#include "ROL_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "MatrixMarket_Tpetra.hpp"

namespace ROL {


namespace TPMultiVector {


  // Locally define a Kokkos wrapper functor for UnaryFunction
  template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type,
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
  struct unaryFunc {
    typedef typename Tpetra::MultiVector<Real,LO,GO,Node>::dual_view_type::t_dev ViewType;
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    const Elementwise::UnaryFunction<Real>* const f_;

    unaryFunc(ViewType& X, const Elementwise::UnaryFunction<Real>* const f)
      : X_(X), f_(f) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      const int M = X_.extent(1);
      for(int j=0;j<M;++j) {
        X_(i,j) = f_->apply(X_(i,j));
      }
    }
  };

  // Locally define a Kokkos wrapper functor for BinaryFunction
  template <class Real,
            class LO=Tpetra::Map<>::local_ordinal_type,
            class GO=Tpetra::Map<>::global_ordinal_type,
            class Node=Tpetra::Map<>::node_type >
   struct binaryFunc {
    typedef typename Tpetra::MultiVector<Real,LO,GO,Node>::dual_view_type::t_dev ViewType;
    typedef typename Tpetra::MultiVector<const Real,LO,GO,Node>::dual_view_type::t_dev ConstViewType;
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    ConstViewType Y_;
    const Elementwise::BinaryFunction<Real>* const f_;

    binaryFunc(ViewType& X, const ConstViewType& Y, const Elementwise::BinaryFunction<Real>* const f)
      : X_(X), Y_(Y), f_(f) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      const int M = X_.extent(1);
      for(int j=0;j<M;++j) {
        X_(i,j) = f_->apply(X_(i,j),Y_(i,j));
      }
    }
  };

  // Locally define a Kokkos wrapper functor for ReductionOp
  template <class Real,
            class LO=Tpetra::Map<>::local_ordinal_type,
            class GO=Tpetra::Map<>::global_ordinal_type,
            class Node=Tpetra::Map<>::node_type >
  struct reduceFunc {
    typedef typename Tpetra::MultiVector<const Real,LO,GO,Node>::dual_view_type::t_dev ConstViewType;
    typedef typename ConstViewType::execution_space execution_space;
    ConstViewType X_;
    const Elementwise::ReductionOp<Real>* const r_;

    reduceFunc(const ConstViewType& X, const Elementwise::ReductionOp<Real>* const r)
      : X_(X), r_(r) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, Real &rval) const {
      const int M = X_.extent(1);

      for(int j=0;j<M;++j) {
        r_->reduce(X_(i,j),rval);
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init(Real &rval) const {
      rval = r_->initialValue();
    }

    KOKKOS_INLINE_FUNCTION
    void join(Real &globalVal, const Real &localVal) const {
      r_->reduce(localVal,globalVal);
    }

  };


} // end namespace TPMultiVector


// Template on the Real/Scalar type, Local Ordinal, Global Ordinal, and Node
template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type,
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class TpetraMultiVector : public Vector<Real> {
private:
  const Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > tpetra_vec_;
  const Ptr<const Tpetra::Map<LO,GO,Node> > map_;
  const Ptr<const Teuchos::Comm<int> > comm_;

protected:
  const Ptr<const Tpetra::Map<LO,GO,Node> > getMap(void) const {
    return map_;
  }

public:
  virtual ~TpetraMultiVector() {}

  TpetraMultiVector(const Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec)
    : tpetra_vec_(tpetra_vec), map_(tpetra_vec_->getMap()), comm_(map_->getComm()) {}

  /** \brief Assign \f$y \leftarrow x \f$ where \f$y = \mbox{*this}\f$.
  */
  void set(const Vector<Real> &x) {
    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );
    const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
    tpetra_vec_->assign(*ex.getVector());
  }

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus(const Vector<Real> &x) {
    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );
    Real one(1);
    const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
    tpetra_vec_->update(one,*ex.getVector(),one);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );
    Real one(1);
    const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
    tpetra_vec_->update(alpha,*ex.getVector(),one);
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) {
    tpetra_vec_->scale( alpha );
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  virtual Real dot( const Vector<Real> &x ) const {
    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );
    const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
    size_t n = tpetra_vec_->getNumVectors();
    // Perform Euclidean dot between *this and x for each vector
    Teuchos::Array<Real> val(n,0);
    tpetra_vec_->dot( *ex.getVector(), val.view(0,n) );
    // Combine dots for each vector to get a scalar
    Real xy(0);
    for (size_t i = 0; i < n; ++i) {
      xy += val[i];
    }
    return xy;
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    return std::sqrt(dot(*this));
  }

  /** \brief Clone to make a new (uninitialized) vector.
  */
  virtual Ptr<Vector<Real> > clone() const {
    size_t n = tpetra_vec_->getNumVectors();
    return makePtr<TpetraMultiVector>(
           makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(map_,n));
  }

  /**  \brief Set to zero vector.
  */
  void zero() {
    Real zero(0);
    tpetra_vec_->putScalar(zero);
  }

  void setScalar(const Real C) {
    tpetra_vec_->putScalar(static_cast<double>(C));
  }

  void randomize(const Real l=0.0, const Real u = 1.0) {
    tpetra_vec_->randomize(static_cast<double>(l),static_cast<double>(u));
  }

  Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > getVector() const {
    return tpetra_vec_;
  }

  Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > getVector() {
    return tpetra_vec_;
  }

  Ptr<Vector<Real> > basis( const int i ) const {
    ROL_TEST_FOR_EXCEPTION( i >= dimension() || i<0,
                                std::invalid_argument,
                                "Error: Basis index must be between 0 and vector dimension." );
    Ptr<Vector<Real> > e = clone();
    const size_t n = tpetra_vec_->getNumVectors();
    if ( (map_ != nullPtr) && map_->isNodeGlobalElement(static_cast<GO>(i))) {
      for (size_t j = 0; j < n; ++j) {
        (staticPtrCast<TpetraMultiVector<Real,LO,GO,Node> >(e)->getVector())->replaceGlobalValue(i, j, ScalarTraits<Real>::one());
      }
    }
    return e;
  }

  int dimension() const {
    int nVecs = static_cast<int>(tpetra_vec_->getNumVectors());
    int dim   = static_cast<int>(tpetra_vec_->getGlobalLength());
    return nVecs*dim;
  }

/*****************************************************************************/
/************** ELEMENTWISE DEFINITIONS AND IMPLEMENTATIONS ******************/
/*****************************************************************************/
private:
  typedef typename Tpetra::MultiVector<Real,LO,GO,Node>::dual_view_type::t_dev ViewType;
  typedef typename Tpetra::MultiVector<const Real,LO,GO,Node>::dual_view_type::t_dev ConstViewType;
public:
  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    ViewType v_lcl =  tpetra_vec_->getLocalViewDevice(Tpetra::Access::ReadWrite);

    int lclDim = tpetra_vec_->getLocalLength();
    TPMultiVector::unaryFunc<Real,LO,GO,Node> func(v_lcl, &f);

    Kokkos::parallel_for(lclDim,func);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

   const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
   Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > xp = ex.getVector();

    ViewType v_lcl = tpetra_vec_->getLocalViewDevice(Tpetra::Access::ReadWrite);
    ConstViewType x_lcl = xp->getLocalViewDevice(Tpetra::Access::ReadOnly);

    int lclDim = tpetra_vec_->getLocalLength();

    TPMultiVector::binaryFunc<Real,LO,GO,Node> func(v_lcl,x_lcl,&f);

    Kokkos::parallel_for(lclDim,func);

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    ConstViewType v_lcl = tpetra_vec_->getLocalViewDevice(Tpetra::Access::ReadOnly);

    int lclDim = tpetra_vec_->getLocalLength();
    TPMultiVector::reduceFunc<Real,LO,GO,Node> func(v_lcl, &r);

    Real lclValue;

    // Reduce for this MPI process
    Kokkos::parallel_reduce(lclDim,func,lclValue);

    Real gblValue;

    auto reductionType = static_cast<Teuchos::EReductionType>(r.reductionType());

    // Reduce over MPI processes
    Teuchos::reduceAll<int,Real>(*comm_,reductionType,lclValue,Teuchos::outArg(gblValue));

    return gblValue;
  }

  void print( std::ostream &outStream ) const {
    //tpetra_vec_->print(outStream);
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<>> vecWriter;
    vecWriter.writeDense(outStream, *tpetra_vec_);
    vecWriter.writeMap(outStream, *tpetra_vec_->getMap());
  }

}; // class TpetraMultiVector


} // namespace ROL

#endif

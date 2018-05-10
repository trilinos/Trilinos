// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_TPETRAMULTIVECTOR_HPP
#define ROL_TPETRAMULTIVECTOR_HPP

/** \class ROL::TpetraMultiVector
    \brief Implements the ROL::Vector interface for a Tpetra_MultiVector.
    \author Created by Greg von Winckel
*/



#include "ROL_Vector.hpp"
#include "Tpetra_MultiVector_def.hpp"

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
      const int M = X_.dimension_1();
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
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    ViewType Y_;
    const Elementwise::BinaryFunction<Real>* const f_;

    binaryFunc(ViewType& X, const ViewType& Y, const Elementwise::BinaryFunction<Real>* const f)
      : X_(X), Y_(Y), f_(f) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      const int M = X_.dimension_1();
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
    typedef typename Tpetra::MultiVector<Real,LO,GO,Node>::dual_view_type::t_dev ViewType;
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    const Elementwise::ReductionOp<Real>* const r_;

    reduceFunc(const ViewType& X, const Elementwise::ReductionOp<Real>* const r)
      : X_(X), r_(r) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i, Real &rval) const {
      const int M = X_.dimension_1();
  
      for(int j=0;j<M;++j) {
        r_->reduce(X_(i,j),rval);
      }
    } 

    KOKKOS_INLINE_FUNCTION
    void init(Real &rval) const {
      rval = r_->initialValue();
    } 

    KOKKOS_INLINE_FUNCTION
    void join(volatile Real &globalVal, const volatile Real &localVal) const {
      r_->reduce(localVal,globalVal);
    } 
       
  };



}


// Template on the Real/Scalar type, Local Ordinal, Global Ordinal, and Node
template <class Real, 
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type, 
          class Node=Tpetra::Map<>::node_type >
class TpetraMultiVector : public Vector<Real> {
private:
  const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > tpetra_vec_;
  const ROL::Ptr<const Tpetra::Map<LO,GO,Node> > map_;
  const ROL::Ptr<const Teuchos::Comm<int> > comm_;

protected:
  const ROL::Ptr<const Tpetra::Map<LO,GO,Node> > getMap(void) const {
    return map_;
  }

public:
  virtual ~TpetraMultiVector() {}

  TpetraMultiVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec)
    : tpetra_vec_(tpetra_vec), map_(tpetra_vec_->getMap()), comm_(map_->getComm()) {}

  /** \brief Assign \f$y \leftarrow x \f$ where \f$y = \mbox{*this}\f$.
  */
  void set(const Vector<Real> &x) {
    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );
    const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
    tpetra_vec_->assign(*ex.getVector());
  }

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus(const Vector<Real> &x) {
    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );
    Real one(1);
    const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
    tpetra_vec_->update(one,*ex.getVector(),one);
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
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
    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
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
  virtual ROL::Ptr<Vector<Real> > clone() const {
    size_t n = tpetra_vec_->getNumVectors();
    return ROL::makePtr<TpetraMultiVector>(
           ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(map_,n));
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

  ROL::Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > getVector() const {
    return tpetra_vec_;
  }

  ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > getVector() {
    return tpetra_vec_;
  }

  ROL::Ptr<Vector<Real> > basis( const int i ) const {
    TEUCHOS_TEST_FOR_EXCEPTION( i >= dimension() || i<0,
                                std::invalid_argument,
                                "Error: Basis index must be between 0 and vector dimension." );
    const size_t n = tpetra_vec_->getNumVectors();
    ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > e
      = ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(map_,n);
    if ( (map_ != ROL::nullPtr) && map_->isNodeGlobalElement(static_cast<GO>(i))) {
      for (size_t j = 0; j < n; ++j) {
        e->replaceGlobalValue (i, j, Teuchos::ScalarTraits<Real>::one());
      }
    }
    return ROL::makePtr<TpetraMultiVector>(e);
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

public:
  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    ViewType v_lcl =  tpetra_vec_->getDualView().d_view;
        
    int lclDim = tpetra_vec_->getLocalLength();
    TPMultiVector::unaryFunc<Real,LO,GO,Node> func(v_lcl, &f);

    Kokkos::parallel_for(lclDim,func);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {

    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

   const TpetraMultiVector &ex = dynamic_cast<const TpetraMultiVector&>(x);
   ROL::Ptr<const Tpetra::MultiVector<Real,LO,GO,Node> > xp = ex.getVector();

    ViewType v_lcl = tpetra_vec_->getDualView().d_view;
    ViewType x_lcl = xp->getDualView().d_view;
       
    int lclDim = tpetra_vec_->getLocalLength();

    TPMultiVector::binaryFunc<Real,LO,GO,Node> func(v_lcl,x_lcl,&f);

    Kokkos::parallel_for(lclDim,func);

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    ViewType v_lcl = tpetra_vec_->getDualView().d_view;

    int lclDim = tpetra_vec_->getLocalLength();
    TPMultiVector::reduceFunc<Real,LO,GO,Node> func(v_lcl, &r);

    Real lclValue;

    // Reduce for this MPI process
    Kokkos::parallel_reduce(lclDim,func,lclValue);

    Real gblValue;

    // Reduce over MPI processes
    Teuchos::reduceAll<int,Real>(*comm_,r.reductionType(),lclValue,Teuchos::outArg(gblValue));

    return gblValue; 
  }

  void print( std::ostream &outStream ) const {
    tpetra_vec_->print(outStream);
  }

}; // class TpetraMultiVector



} // namespace ROL

#endif

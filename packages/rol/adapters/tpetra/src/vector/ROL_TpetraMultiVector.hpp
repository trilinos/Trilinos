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

// Template on the Real/Scalar type, Local Ordinal, Global Ordinal, and Node
template <class Real, 
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type, 
          class Node=Tpetra::Map<>::node_type >
class TpetraMultiVector : public Vector<Real> {

  typedef Vector<Real>                         V;
  typedef Teuchos::RCP<V>                      VP;
  typedef Tpetra::MultiVector<Real,LO,GO,Node> MV;
  typedef Teuchos::RCP<MV>                     MVP;
  typedef Teuchos::RCP<const MV>               CMVP;
  typedef TpetraMultiVector                    TMV;
  typedef Teuchos::RCP<TMV>                    TMVP; 
  typedef Elementwise::UnaryFunction<Real>     UF;
  typedef Elementwise::BinaryFunction<Real>    BF; 
  typedef Elementwise::ReductionOp<Real>       RO;

  typedef typename MV::dual_view_type::t_dev    ViewType;
           

  // Locally define a Kokkos wrapper functor for UnaryFunction  
  struct unaryFunc {
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    const UF* const f_;
  
    unaryFunc(ViewType& X, const UF* const f) : X_(X), f_(f) {}  

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      const int M = X_.dimension_1();
      for(int j=0;j<M;++j) {
        X_(i,j) = f_->apply(X_(i,j));
      } 
    }
  };


  // Locally define a Kokkos wrapper functor for BinaryFunction  
  struct binaryFunc {
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    ViewType Y_;
    const BF* const f_;

    binaryFunc(ViewType& X, const ViewType& Y, const BF* const f) : X_(X), Y_(Y), f_(f) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int i) const {
      const int M = X_.dimension_1();
      for(int j=0;j<M;++j) {
        X_(i,j) = f_->apply(X_(i,j),Y_(i,j));    
      }
    } 
  };

  // Locally define a Kokkos wrapper functor for ReductionOp
  struct reduceFunc {
    typedef typename ViewType::execution_space execution_space;
    ViewType X_;
    const RO* const r_;

    reduceFunc(const ViewType& X, const RO* const r) : X_(X), r_(r) {}

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

  private:
    MVP tpetra_vec_;

    Teuchos::RCP<const Tpetra::Comm<int> > comm_;  

  public:
    virtual ~TpetraMultiVector() {}

    TpetraMultiVector(const MVP &tpetra_vec) : tpetra_vec_(tpetra_vec),
      comm_(tpetra_vec_->getMap()->getComm()) {}
       
    /** \brief Assign \f$y \leftarrow x \f$ where \f$y = \mbox{*this}\f$.
    */
    void set(const Vector<Real> &x) {

    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

      const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
      tpetra_vec_->assign(*ex.getVector()); 
    }
       
    /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
    */
    void plus(const V &x) {

    TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

      const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
      tpetra_vec_->update(1.0,*ex.getVector(),1.0); 
    }


    void axpy( const Real alpha, const Vector<Real> &x ) {

      TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );

      const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
      tpetra_vec_->update(alpha,*ex.getVector(),1.0); 

    }

    /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
    */
    void scale( const Real alpha ) {
      tpetra_vec_->scale( alpha ); 
    }

    /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
    */
    virtual Real dot( const V &x ) const {

      TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );

      Real v; // Need this to make a 1-element ArrayView
      Teuchos::ArrayView<Real> val(&v,1);
      const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
      tpetra_vec_->dot( *ex.getVector(), val );
      return val[0];
    }

    /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
    */
    Real norm() const {
      Real v; // Need this to make a 1-element ArrayView
      Teuchos::ArrayView<Real> val(&v,1);
      tpetra_vec_->norm2(val);
      return val[0];
    } 

    /** \brief Clone to make a new (uninitialized) vector.
    */
    virtual VP clone() const {
      using Teuchos::rcp; 
      size_t n = tpetra_vec_->getNumVectors();  
      return rcp( new TMV(rcp( new MV(tpetra_vec_->getMap(),n,false)) ));
    }

    /**  \brief Set to zero vector.
    */
    void zero() {
      tpetra_vec_->putScalar(0.0);
    }

    CMVP getVector() const {
      return this->tpetra_vec_;
    }

    MVP getVector() {
      return this->tpetra_vec_;
    }

    VP basis( const int i ) const {

      TEUCHOS_TEST_FOR_EXCEPTION( i >= dimension() || i<0,
                                  std::invalid_argument,
                                  "Error: Basis index must be between 0 and vector dimension." );

      using Teuchos::RCP; 
      using Teuchos::rcp; 
      typedef typename MV::map_type map_type;

      const size_t n = tpetra_vec_->getNumVectors();  
          
      RCP<const map_type> map = tpetra_vec_->getMap ();
      MVP e = rcp (new MV (map,n));  

      if (! map.is_null () && map->isNodeGlobalElement (static_cast<GO> (i))) {
        for (size_t j = 0; j < n; ++j) {
          e->replaceGlobalValue (i, j, Teuchos::ScalarTraits<Real>::one ());
        }
      }
      return rcp(new TMV(e) );  
    }

    int dimension() const { return tpetra_vec_->getGlobalLength(); }


    void applyUnary( const UF &f ) {
      ViewType v_lcl =  tpetra_vec_->getDualView().d_view;
          
      int lclDim = tpetra_vec_->getLocalLength();
      unaryFunc func(v_lcl, &f);

      Kokkos::parallel_for(lclDim,func);
    }

    void applyBinary( const BF &f, const V &x ) {

      TEUCHOS_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
 
     const TMV &ex = Teuchos::dyn_cast<const TMV>(x);
      CMVP xp = ex.getVector();

      ViewType v_lcl = tpetra_vec_->getDualView().d_view;
      ViewType x_lcl = xp->getDualView().d_view;
         
      int lclDim = tpetra_vec_->getLocalLength();

      binaryFunc func(v_lcl,x_lcl,&f);

      Kokkos::parallel_for(lclDim,func);

    }

    Real reduce( const RO &r ) const {
      ViewType v_lcl = tpetra_vec_->getDualView().d_view;

      int lclDim = tpetra_vec_->getLocalLength();
      reduceFunc func(v_lcl, &r);

      Real lclValue;
 
      // Reduce for this MPI process
      Kokkos::parallel_reduce(lclDim,func,lclValue);

      Real gblValue;

      // Reduce over MPI processes
      using Teuchos::outArg;
      using Teuchos::reduceAll;
      reduceAll<int,Real>(*comm_,r.reductionType(),lclValue,outArg(gblValue));

      return gblValue; 
    }

}; // class TpetraMultiVector



} // namespace ROL

#endif

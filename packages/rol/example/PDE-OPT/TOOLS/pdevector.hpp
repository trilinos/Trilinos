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

#ifndef ROL_PDEOPT_PDEVECTOR_HPP
#define ROL_PDEOPT_PDEVECTOR_HPP

#include "ROL_TpetraMultiVector.hpp"

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalSimVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualSimVector;

template <class Real, class LO, class GO, class Node>
class PDE_PrimalSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_DualSimVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_PrimalSimVector() {}

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap1(pde);
    }

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDERieszMap1(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PDE_PrimalSimVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = Teuchos::rcp(new PDE_DualSimVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                    assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDERieszMap1(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }

}; // class PrimalScaledTpetraMultiVector

template <class Real, class LO, class GO, class Node>
class PDE_DualSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_PrimalSimVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_DualSimVector() {}

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap1(pde);
    }

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDEInverseRieszMap1(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new PDE_DualSimVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = Teuchos::rcp(new PDE_PrimalSimVector<Real,LO,GO,Node>(
                      Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                      assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDEInverseRieszMap1(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }

}; // class DualScaledTpetraMultiVector

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalOptVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualOptVector;

template <class Real, class LO, class GO, class Node>
class PDE_PrimalOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_DualOptVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_PrimalOptVector() {}

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap2(pde);
    }

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDERieszMap2(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PDE_PrimalOptVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = Teuchos::rcp(new PDE_DualOptVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                    assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDERieszMap2(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }

}; // class PrimalScaledTpetraMultiVector

template <class Real, class LO, class GO, class Node>
class PDE_DualOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const Teuchos::RCP<Assembler<Real> > assembler_;

    mutable Teuchos::RCP<PDE_PrimalOptVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

  public:
    virtual ~PDE_DualOptVector() {}

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {
      assembler_->assemblePDERieszMap2(pde);
    }

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), assembler_(assembler),
        isDualInitialized_(false) {}

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      assembler_->applyPDEInverseRieszMap2(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new PDE_DualOptVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             assembler_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = Teuchos::rcp(new PDE_PrimalOptVector<Real,LO,GO,Node>(
                      Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                      assembler_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      assembler_->applyPDEInverseRieszMap2(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }

}; // class DualScaledTpetraMultiVector

template<class Real> 
class PDE_OptVector : public ROL::Vector<Real> {
private:
  Teuchos::RCP<ROL::TpetraMultiVector<Real> > vec1_;
  Teuchos::RCP<ROL::StdVector<Real> > vec2_;
  mutable Teuchos::RCP<ROL::Vector<Real> > dual_vec1_;
  mutable Teuchos::RCP<ROL::Vector<Real> > dual_vec2_;
  mutable Teuchos::RCP<PDE_OptVector<Real> > dual_vec_;

public:
  PDE_OptVector(const Teuchos::RCP<ROL::TpetraMultiVector<Real> > &vec1,
                const Teuchos::RCP<ROL::StdVector<Real> > &vec2 ) 
    : vec1_(vec1), vec2_(vec2) {
    dual_vec1_ = (vec1_->dual()).clone();
    dual_vec2_ = (vec2_->dual()).clone();
  }

  PDE_OptVector(const Teuchos::RCP<ROL::TpetraMultiVector<Real> > &vec)
    : vec1_(vec), vec2_(Teuchos::null), dual_vec2_(Teuchos::null) {
    dual_vec1_ = (vec1_->dual()).clone();
  }

  PDE_OptVector(const Teuchos::RCP<ROL::StdVector<Real> > &vec)
    : vec1_(Teuchos::null), vec2_(vec), dual_vec1_(Teuchos::null) {
    dual_vec2_ = (vec2_->dual()).clone();
  }

  void plus( const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->plus(*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->plus(*(xs.getParameter()));
    }
  }

  void scale( const Real alpha ) {
    if ( vec1_ != Teuchos::null ) {
      vec1_->scale(alpha);
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->scale(alpha);
    }
  }

  void axpy( const Real alpha, const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->axpy(alpha,*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->axpy(alpha,*(xs.getParameter()));
    }
  }

  Real dot( const ROL::Vector<Real> &x ) const {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    Real val(0);
    if ( vec1_ != Teuchos::null ) {
      val += vec1_->dot(*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      val += vec2_->dot(*(xs.getParameter()));
    }
    return val;
  }

  Real norm() const {
    Real val(0);
    if ( vec1_ != Teuchos::null ) {
      Real norm1 = vec1_->norm();
      val += norm1*norm1;
    }
    if ( vec2_ != Teuchos::null ) {
      Real norm2 = vec2_->norm();
      val += norm2*norm2;
    }
    return std::sqrt(val);
  } 

  Teuchos::RCP<ROL::Vector<Real> > clone() const {
    if ( vec1_ != Teuchos::null && vec2_ != Teuchos::null ) {
      return Teuchos::rcp( new PDE_OptVector(vec1_->clone(),vec2_->clone()) );
    }
    if ( vec1_ != Teuchos::null && vec2_ == Teuchos::null ) {
      return Teuchos::rcp( new PDE_OptVector(vec1_->clone(),Teuchos::null) );
    }
    if ( vec1_ == Teuchos::null && vec2_ != Teuchos::null ) {
      return Teuchos::rcp( new PDE_OptVector(Teuchos::null,vec2_->clone()) );
    }
    return Teuchos::null;
  }

  const ROL::Vector<Real> & dual(void) const {
    if ( vec1_ != Teuchos::null ) {
      dual_vec1_->set(vec1_->dual());
    }
    if ( vec2_ != Teuchos::null ) {
      dual_vec2_->set(vec2_->dual());
    }
    dual_vec_ = Teuchos::rcp( new PDE_OptVector<Real>(dual_vec1_,dual_vec2_) ); 
    return *dual_vec_;
  }

  Teuchos::RCP<ROL::Vector<Real> > basis( const int i )  const {
    if ( vec1_ != Teuchos::null && vec2_ != Teuchos::null ) {
      int n1 = vec1_->dimension();
      if ( i < n1 ) {
        Teuchos::RCP<ROL::Vector<Real> > e1 = vec1_->basis(i);
        Teuchos::RCP<ROL::Vector<Real> > e2 = vec2_->clone(); e2->zero();
        Teuchos::RCP<ROL::Vector<Real> > e  = Teuchos::rcp(new PDE_OptVector(e1,e2));
        return e;
      }
      else {
        Teuchos::RCP<ROL::Vector<Real> > e1 = vec1_->clone(); e1->zero();
        Teuchos::RCP<ROL::Vector<Real> > e2 = vec2_->basis(i-n1);
        Teuchos::RCP<ROL::Vector<Real> > e  = Teuchos::rcp(new PDE_OptVector<Real>(e1,e2));
        return e;
      }
    }
    if ( vec1_ != Teuchos::null && vec2_ == Teuchos::null ) {
      int n1 = vec1_->dimension();
      if ( i < n1 ) {
        Teuchos::RCP<ROL::Vector<Real> > e1 = vec1_->basis(i);
        Teuchos::RCP<ROL::Vector<Real> > e  = Teuchos::rcp(new PDE_OptVector(e1,Teuchos::null));
        return e;
      }
    }
    if ( vec1_ == Teuchos::null && vec2_ != Teuchos::null ) {
      int n2 = vec2_->dimension();
      if ( i < n2 ) {
        Teuchos::RCP<ROL::Vector<Real> > e2 = vec2_->basis(i);
        Teuchos::RCP<ROL::Vector<Real> > e  = Teuchos::rcp(new PDE_OptVector(Teuchos::null,e2));
        return e;
      }
    }
  }

  void applyUnary( const ROL::Elementwise::UnaryFunction<Real> &f ) {
    if ( vec1_ != Teuchos::null ) {
      vec1_->applyUnary(f);
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->applyUnary(f);
    }
  }

  void applyBinary( const ROL::Elementwise::BinaryFunction<Real> &f, const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->applyBinary(f,*xs.getField());
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->applyBinary(f,*xs.getParameter());
    }
  }

  Real reduce( const ROL::Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    if ( vec1_ != Teuchos::null ) {
      r.reduce(vec1_->reduce(r),result);
    }
    if ( vec2_ != Teuchos::null ) {
      r.reduce(vec2_->reduce(r),result);
    }
    return result;
  }

  int dimension() const {
    int dim(0);
    if ( vec1_ != Teuchos::null ) {
      dim += vec1_->dimension();
    }
    if ( vec2_ != Teuchos::null ) {
      dim += vec2_->dimension();
    }
    return dim;
  }

  Teuchos::RCP<const ROL::Vector<Real> > getField() const { 
    return vec1_;
  }

  Teuchos::RCP<const ROL::Vector<Real> > getParameter() const { 
    return vec2_; 
  }

  Teuchos::RCP<ROL::Vector<Real> > getField() { 
    return vec1_;
  }

  Teuchos::RCP<ROL::Vector<Real> > getParameter() { 
    return vec2_; 
  }

  void set_1(const ROL::Vector<Real>& vec) { 
    vec1_->set(vec);
  }
  
  void set_2(const ROL::Vector<Real>& vec) { 
    vec2_->set(vec); 
  }
};

#endif

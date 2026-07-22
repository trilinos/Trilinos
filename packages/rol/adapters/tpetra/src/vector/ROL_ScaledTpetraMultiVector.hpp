// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALEDTPETRAMULTIVECTOR_HPP
#define ROL_SCALEDTPETRAMULTIVECTOR_HPP

/** \class ROL::ScaledTpetraMultiVector
    \brief Implements the ROL::Vector interface for a Tpetra_MultiVector
           with scaled dot function.
    \author Created by Drew Kouri
*/

#include "ROL_TpetraMultiVector.hpp"

namespace ROL {

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PrimalScaledTpetraMultiVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class DualScaledTpetraMultiVector;

template <class Real, class LO, class GO, class Node>
class PrimalScaledTpetraMultiVector : public TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const ROL::Ptr<const Tpetra::Vector<Real,LO,GO,Node> > scale_vec_;
    mutable ROL::Ptr<DualScaledTpetraMultiVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(Tpetra::MultiVector<Real,LO,GO,Node> &out,
                const Tpetra::MultiVector<Real,LO,GO,Node> &in) const {
      Real zero(0), one(1);
      out.elementWiseMultiply(one,*scale_vec_,in,zero);
    }

  public:
    virtual ~PrimalScaledTpetraMultiVector() {}

    PrimalScaledTpetraMultiVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                                  const ROL::Ptr<const Tpetra::Vector<Real,LO,GO,Node> > &scale_vec)
      : TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        scale_vec_(scale_vec), isDualInitialized_(false) {}

    Real dot( const Vector<Real> &x ) const {
      ROL_TEST_FOR_EXCEPTION( (TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Tpetra::MultiVector<Real,LO,GO,Node> &ex
        = *(dynamic_cast<const TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector());
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Tpetra::MultiVector<Real,LO,GO,Node> wex(TpetraMultiVector<Real>::getMap(), n);
      applyScaling(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    ROL::Ptr<Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return ROL::makePtr<PrimalScaledTpetraMultiVector<Real,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(TpetraMultiVector<Real>::getMap(),n),
             scale_vec_);
    }

    const Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = ROL::makePtr<DualScaledTpetraMultiVector<Real,LO,GO,Node>>(
                    ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(TpetraMultiVector<Real>::getMap(),n),
                    scale_vec_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyScaling(*(dual_vec_->getVector()),*(TpetraMultiVector<Real,LO,GO,Node>::getVector()));
      return *dual_vec_;
    }

    Real apply(const Vector<Real> &x) const {
      const DualScaledTpetraMultiVector<Real,LO,GO,Node> & ex = dynamic_cast<const DualScaledTpetraMultiVector<Real,LO,GO,Node>&>(x);
      return TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }

}; // class PrimalScaledTpetraMultiVector

template <class Real, class LO, class GO, class Node>
class DualScaledTpetraMultiVector : public TpetraMultiVector<Real,LO,GO,Node> {
  private:
    const ROL::Ptr<const Tpetra::Vector<Real,LO,GO,Node> > scale_vec_;
    mutable ROL::Ptr<PrimalScaledTpetraMultiVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(Tpetra::MultiVector<Real,LO,GO,Node> &out,
                const Tpetra::MultiVector<Real,LO,GO,Node> &in) const {
      Real zero(0), one(1);
      Tpetra::Vector<Real,LO,GO,Node> W(scale_vec_->getMap());
      W.reciprocal(*scale_vec_);
      out.elementWiseMultiply(one,W,in,zero);
    }

  public:
    virtual ~DualScaledTpetraMultiVector() {}

    DualScaledTpetraMultiVector(const ROL::Ptr<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                                const ROL::Ptr<const Tpetra::Vector<Real,LO,GO,Node> > &scale_vec)
      : TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        scale_vec_(scale_vec), isDualInitialized_(false) {}

    Real dot( const Vector<Real> &x ) const {
      ROL_TEST_FOR_EXCEPTION( (TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Tpetra::MultiVector<Real,LO,GO,Node> &ex
        = *(dynamic_cast<const TpetraMultiVector<Real,LO,GO,Node>&>(x).getVector());
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Tpetra::MultiVector<Real,LO,GO,Node> wex(TpetraMultiVector<Real>::getMap(), n);
      applyScaling(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    ROL::Ptr<Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return ROL::makePtr<DualScaledTpetraMultiVector<Real,LO,GO,Node>>(
             ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(TpetraMultiVector<Real>::getMap(),n),
             scale_vec_);
    }

    const Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = ROL::makePtr<PrimalScaledTpetraMultiVector<Real,LO,GO,Node>>(
                      ROL::makePtr<Tpetra::MultiVector<Real,LO,GO,Node>>(TpetraMultiVector<Real>::getMap(),n),
                      scale_vec_);
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyScaling(*(primal_vec_->getVector()),*(TpetraMultiVector<Real,LO,GO,Node>::getVector()));
      return *primal_vec_;
    }

    Real apply(const Vector<Real> &x) const {
      const PrimalScaledTpetraMultiVector<Real,LO,GO,Node> & ex = dynamic_cast<const PrimalScaledTpetraMultiVector<Real,LO,GO,Node>&>(x);
      return TpetraMultiVector<Real,LO,GO,Node>::dot(ex);
    }

}; // class DualScaledTpetraMultiVector

} // namespace ROL

#endif

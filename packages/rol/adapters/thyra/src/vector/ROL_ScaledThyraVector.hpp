// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALEDTHYRAVECTOR_HPP
#define ROL_SCALEDTHYRAVECTOR_HPP

/** \class ROL::ScaledThyraVector
    \brief Implements the ROL::Vector interface for a Thyra::Vector
           with scaled dot function.
    \author Created by Eric Cyr
*/

#include "ROL_ThyraVector.hpp"

namespace ROL {

template <class Real>
class PrimalScaledThyraVector;

template <class Real>
class DualScaledThyraVector;

template <class Real>
class PrimalScaledThyraVector : public ThyraVector<Real> {
  private:
    const Teuchos::RCP<const Thyra::VectorBase<Real> >  scaling_vec_;

    mutable Teuchos::RCP<DualScaledThyraVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(const Teuchos::RCP<Thyra::VectorBase<Real> > &out,
                      const Teuchos::RCP<const Thyra::VectorBase<Real> > & in) const 
    {
      Real one(1);
      ::Thyra::assign(out.ptr(),0.0);
      ::Thyra::ele_wise_prod(one,*in,*scaling_vec_,out.ptr());
    }

  public:
    virtual ~PrimalScaledThyraVector() {}

    PrimalScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > & thyra_vec,
                            const Teuchos::RCP<const Thyra::VectorBase<Real> > & scale_vec)
      : ThyraVector<Real>(thyra_vec),
        scaling_vec_(scale_vec), isDualInitialized_(false) {}

    Real dot( const Vector<Real> &x ) const {
      Teuchos::RCP<const Thyra::VectorBase<Real> > ex 
          = dynamic_cast<const ThyraVector<Real>&>(x).getVector();

      // compute a scaled version of the "this" vector
      Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
      Teuchos::RCP<Thyra::VectorBase<Real> > scaled_vec = vec->clone_v();
      applyScaling(scaled_vec,vec);

      // compute this vector but scaled
      return ::Thyra::dot<Real>(*scaled_vec, *ex);
    }

    Teuchos::RCP<Vector<Real> > clone() const {
      Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
      return Teuchos::rcp( new PrimalScaledThyraVector<Real>(vec->clone_v(),scaling_vec_) );
    }

    const Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
        dual_vec_ = Teuchos::rcp( new DualScaledThyraVector<Real>(vec->clone_v(),scaling_vec_) );
        isDualInitialized_ = true;
      }
      // Scale this with scaling_vec_ and place in dual vector
      applyScaling(dual_vec_->getVector(),ThyraVector<Real>::getVector());
      return *dual_vec_;
    }

    Real apply( const Vector<Real> &x ) const {
      const DualScaledThyraVector<Real> &ex = dynamic_cast<const DualScaledThyraVector<Real>&>(x);
      return ThyraVector<Real>::dot(ex);
    }

}; // class PrimalScaledThyraVector

template <class Real>
class DualScaledThyraVector : public ThyraVector<Real> {
  private:
    const Teuchos::RCP<const Thyra::VectorBase<Real> > scaling_vec_;
    mutable Teuchos::RCP<PrimalScaledThyraVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(const Teuchos::RCP<Thyra::VectorBase<Real> > &out,
                      const Teuchos::RCP<const Thyra::VectorBase<Real> > & in) const 
    {
      Real one(1);
      ::Thyra::assign(out.ptr(),0.0);
      ::Thyra::ele_wise_divide(one,*in,*scaling_vec_,out.ptr());
    }

  public:
    virtual ~DualScaledThyraVector() {}

    DualScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > & thyra_vec,
                          const Teuchos::RCP<const Thyra::VectorBase<Real> > & scale_vec)
      : ThyraVector<Real>(thyra_vec),
        scaling_vec_(scale_vec), isDualInitialized_(false) {}

    Real dot( const Vector<Real> &x ) const {
      Teuchos::RCP<const Thyra::VectorBase<Real> > ex 
          = dynamic_cast<const ThyraVector<Real>&>(x).getVector();

      // compute a scaled version of the "this" vector
      Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
      Teuchos::RCP<Thyra::VectorBase<Real> > scaled_vec = vec->clone_v();
      applyScaling(scaled_vec,vec);

      // compute this vector but scaled
      return ::Thyra::dot<Real>(*scaled_vec, *ex);
    }

    Teuchos::RCP<Vector<Real> > clone() const {
      Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
      return Teuchos::rcp( new DualScaledThyraVector<Real>(vec->clone_v(),scaling_vec_) );
    }

    const Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
        primal_vec_ = Teuchos::rcp( new PrimalScaledThyraVector<Real>(vec->clone_v(),scaling_vec_) );
        isDualInitialized_ = true;
      }
      // Scale this with scaling_vec_ and place in dual vector
      applyScaling(primal_vec_->getVector(),ThyraVector<Real>::getVector());
      return *primal_vec_;
    }

    Real apply( const Vector<Real> &x ) const {
      const PrimalScaledThyraVector<Real> &ex = dynamic_cast<const PrimalScaledThyraVector<Real>&>(x);
      return ThyraVector<Real>::dot(ex);
    }

}; // class DualScaledThyraVector

} // namespace ROL

#endif

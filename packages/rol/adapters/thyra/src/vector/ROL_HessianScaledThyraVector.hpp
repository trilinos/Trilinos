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

#ifndef ROL_HESSIANSCALEDTHYRAVECTOR_HPP
#define ROL_HESSIANSCALEDTHYRAVECTOR_HPP

/** \class ROL::HessianScaledThyraVector
    \brief Implements the ROL::Vector interface for a Thyra::Vector
           with Hessian dot function.
*/

#include "ROL_ThyraVector.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_SolveInverseFactory.hpp"

namespace ROL {

template <class Real>
class PrimalHessianScaledThyraVector;

template <class Real>
class DualHessianScaledThyraVector;

template <class Real>
class PrimalHessianScaledThyraVector : public ThyraVector<Real> {
   private:
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > hessian_;
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invHessian_;

    const bool removeMeanOfTheRHS_;

    mutable Teuchos::RCP<DualHessianScaledThyraVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(const Teuchos::RCP<Thyra::VectorBase<Real> > &out,
                      const Teuchos::RCP<const Thyra::VectorBase<Real> > &in) const {
        if(Teuchos::is_null(hessian_))
            Thyra::copy(*in, out.ptr());
        else
            Thyra::apply(*hessian_, Thyra::NOTRANS, *in, out.ptr(), (Real)1, (Real)0);
    }

   public:
    virtual ~PrimalHessianScaledThyraVector() {}

    PrimalHessianScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > &thyra_vec,
                                   const Teuchos::RCP<const Thyra::LinearOpBase<Real> > hessian,
                                   const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invHessian,
                                   const bool removeMeanOfTheRHS)
        : ThyraVector<Real>(thyra_vec),
          hessian_(hessian),
          invHessian_(invHessian),
          removeMeanOfTheRHS_(removeMeanOfTheRHS),
          isDualInitialized_(false) {}

    PrimalHessianScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > &thyra_vec,
                                   const Teko::BlockedLinearOp hessian,
                                   const Teko::BlockedLinearOp invHessian,
                                   const bool removeMeanOfTheRHS)
        : ThyraVector<Real>(thyra_vec),
          hessian_(Teko::toLinearOp(hessian)),
          invHessian_(Teko::toLinearOp(invHessian)),
          removeMeanOfTheRHS_(removeMeanOfTheRHS),
          isDualInitialized_(false) {}

    Real dot(const Vector<Real> &x) const {
        Teuchos::RCP<const Thyra::VectorBase<Real> > ex = dynamic_cast<const ThyraVector<Real> &>(x).getVector();

        // compute a scaled version of the "this" vector
        Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
        Teuchos::RCP<Thyra::VectorBase<Real> > scaled_vec = vec->clone_v();
        applyScaling(scaled_vec, vec);

        // compute this vector but scaled
        return ::Thyra::dot<Real>(*scaled_vec, *ex);
    }

    Teuchos::RCP<Vector<Real> > clone() const {
        Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
        return Teuchos::rcp(new PrimalHessianScaledThyraVector<Real>(vec->clone_v(), hessian_, invHessian_, removeMeanOfTheRHS_));
    }

    const Vector<Real> &dual() const {
        if (!isDualInitialized_) {
            // Create new memory for dual vector
            Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
            dual_vec_ = Teuchos::rcp(new DualHessianScaledThyraVector<Real>(vec->clone_v(), hessian_, invHessian_, removeMeanOfTheRHS_));
            isDualInitialized_ = true;
        }
        // Scale this with scaling_vec_ and place in dual vector
        applyScaling(dual_vec_->getVector(), ThyraVector<Real>::getVector());
        return *dual_vec_;
    }

};  // class PrimalHessianScaledThyraVector

template <class Real>
class DualHessianScaledThyraVector : public ThyraVector<Real> {
   private:
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > hessian_;
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invHessian_;

    const bool removeMeanOfTheRHS_;

    mutable Teuchos::RCP<PrimalHessianScaledThyraVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(const Teuchos::RCP<Thyra::VectorBase<Real> > &out,
                      const Teuchos::RCP<const Thyra::VectorBase<Real> > &in) const {
        if(Teuchos::is_null(invHessian_))
            Thyra::copy(*in, out.ptr());
        else {
            Teuchos::RCP<Thyra::VectorBase<Real> > in2 = in->clone_v();
            if(removeMeanOfTheRHS_) {
                const Real meanValue = -Thyra::sum(*in2)/(in2->space()->dim());
                Thyra::add_scalar(meanValue, in2.ptr());
            }
            Thyra::apply(*invHessian_, Thyra::NOTRANS, *in2, out.ptr(), (Real)1, (Real)0);
        }
    }

   public:
    virtual ~DualHessianScaledThyraVector() {}

    DualHessianScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > &thyra_vec,
                                 const Teuchos::RCP<const Thyra::LinearOpBase<Real> > hessian,
                                 const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invHessian,
                                 const bool removeMeanOfTheRHS)
        : ThyraVector<Real>(thyra_vec),
          hessian_(hessian),
          invHessian_(invHessian),
          removeMeanOfTheRHS_(removeMeanOfTheRHS),
          isDualInitialized_(false) {}

    DualHessianScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > &thyra_vec,
                                 const Teko::BlockedLinearOp hessian,
                                 const Teko::BlockedLinearOp invHessian,
                                 const bool removeMeanOfTheRHS)
        : ThyraVector<Real>(thyra_vec),
          hessian_(Teko::toLinearOp(hessian)),
          invHessian_(Teko::toLinearOp(invHessian)),
          removeMeanOfTheRHS_(removeMeanOfTheRHS),
          isDualInitialized_(false) {}

    Real dot(const Vector<Real> &x) const {
        Teuchos::RCP<const Thyra::VectorBase<Real> > ex = dynamic_cast<const ThyraVector<Real> &>(x).getVector();

        // compute a scaled version of the "this" vector
        Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
        Teuchos::RCP<Thyra::VectorBase<Real> > scaled_vec = vec->clone_v();
        applyScaling(scaled_vec, vec);

        // compute this vector but scaled
        return ::Thyra::dot<Real>(*scaled_vec, *ex);
    }

    Teuchos::RCP<Vector<Real> > clone() const {
        Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
        return Teuchos::rcp(new DualHessianScaledThyraVector<Real>(vec->clone_v(), hessian_, invHessian_, removeMeanOfTheRHS_));
    }

    const Vector<Real> &dual() const {
        if (!isDualInitialized_) {
            // Create new memory for dual vector
            Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
            primal_vec_ = Teuchos::rcp(new PrimalHessianScaledThyraVector<Real>(vec->clone_v(), hessian_, invHessian_, removeMeanOfTheRHS_));
            isDualInitialized_ = true;
        }
        // Scale this with scaling_vec_ and place in dual vector
        applyScaling(primal_vec_->getVector(), ThyraVector<Real>::getVector());
        return *primal_vec_;
    }

};  // class DualHessianScaledThyraVector

}  // namespace ROL

#endif

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

#ifndef ROL_LINEAROPSCALEDTHYRAVECTOR_HPP
#define ROL_LINEAROPSCALEDTHYRAVECTOR_HPP

/** \class ROL::LinearOpScaledThyraVector
    \brief Implements the ROL::Vector interface for a Thyra::Vector
           with Operator based dot function.
*/

#include "ROL_ThyraVector.hpp"
#include "Thyra_LinearOpBase_decl.hpp"

namespace ROL {

template <class Real>
class PrimalLinearOpScaledThyraVector;

template <class Real>
class DualLinearOpScaledThyraVector;

template <class Real>
class PrimalLinearOpScaledThyraVector : public ThyraVector<Real> {
   private:
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > operator_;
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invOperator_;

    mutable Teuchos::RCP<DualLinearOpScaledThyraVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void applyScaling(const Teuchos::RCP<Thyra::VectorBase<Real> > &out,
                      const Teuchos::RCP<const Thyra::VectorBase<Real> > &in) const {
        if(Teuchos::is_null(operator_))
            Thyra::copy(*in, out.ptr());
        else
            Thyra::apply(*operator_, Thyra::NOTRANS, *in, out.ptr(), (Real)1, (Real)0);
    }

   public:
    virtual ~PrimalLinearOpScaledThyraVector() {}

    PrimalLinearOpScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > &thyra_vec,
                                   const Teuchos::RCP<const Thyra::LinearOpBase<Real> > op,
                                   const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invOp)
        : ThyraVector<Real>(thyra_vec),
          operator_(op),
          invOperator_(invOp),
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
        return Teuchos::rcp(new PrimalLinearOpScaledThyraVector<Real>(vec->clone_v(), operator_, invOperator_));
    }

    const Vector<Real> &dual() const {
        if (!isDualInitialized_) {
            // Create new memory for dual vector
            Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
            dual_vec_ = Teuchos::rcp(new DualLinearOpScaledThyraVector<Real>(vec->clone_v(), operator_, invOperator_));
            isDualInitialized_ = true;
        }
        // Scale this with scaling_vec_ and place in dual vector
        applyScaling(dual_vec_->getVector(), ThyraVector<Real>::getVector());
        return *dual_vec_;
    }

};  // class PrimalLinearOpScaledThyraVector

template <class Real>
class DualLinearOpScaledThyraVector : public ThyraVector<Real> {
   private:
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > operator_;
    const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invOperator_;

    mutable Teuchos::RCP<PrimalLinearOpScaledThyraVector<Real> > primal_vec_;
    mutable bool isPrimalInitialized_;

    void applyScaling(const Teuchos::RCP<Thyra::VectorBase<Real> > &out,
                      const Teuchos::RCP<const Thyra::VectorBase<Real> > &in) const {
        if(Teuchos::is_null(invOperator_))
            Thyra::copy(*in, out.ptr());
        else {
            Thyra::apply(*invOperator_, Thyra::NOTRANS, *in, out.ptr(), (Real)1, (Real)0);
        }
    }

   public:
    virtual ~DualLinearOpScaledThyraVector() {}

    DualLinearOpScaledThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > &thyra_vec,
                                 const Teuchos::RCP<const Thyra::LinearOpBase<Real> > op,
                                 const Teuchos::RCP<const Thyra::LinearOpBase<Real> > invOp)
        : ThyraVector<Real>(thyra_vec),
          operator_(op),
          invOperator_(invOp),
          isPrimalInitialized_(false) {}
          
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
        return Teuchos::rcp(new DualLinearOpScaledThyraVector<Real>(vec->clone_v(), operator_, invOperator_));
    }

    const Vector<Real> &dual() const {
        if (!isPrimalInitialized_) {
            // Create new memory for dual vector
            Teuchos::RCP<const Thyra::VectorBase<Real> > vec = ThyraVector<Real>::getVector();
            primal_vec_ = Teuchos::rcp(new PrimalLinearOpScaledThyraVector<Real>(vec->clone_v(), operator_, invOperator_));
            isPrimalInitialized_ = true;
        }
        // Scale this with scaling_vec_ and place in dual vector
        applyScaling(primal_vec_->getVector(), ThyraVector<Real>::getVector());
        return *primal_vec_;
    }

};  // class DualLinearOpScaledThyraVector

}  // namespace ROL

#endif

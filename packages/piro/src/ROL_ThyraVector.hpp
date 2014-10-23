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

#ifndef ROL_THYRAVECTOR_H
#define ROL_THYRAVECTOR_H

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_Tuple.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"


#include "ROL_Vector.hpp"

/** \class ROL::ThyraVector
    \brief Implements the ROL::Vector interface for an Thyra Vector.
*/

namespace ROL {

template <class Real>
class ThyraVector : public Vector<Real> {
private:

  Teuchos::RCP< Thyra::VectorBase<Real> >  thyra_vec_;

public:
  ~ThyraVector() {}

  ThyraVector(const Teuchos::RCP<Thyra::VectorBase<Real> > & thyra_vec) : thyra_vec_(thyra_vec) {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus( const Vector<Real> &x ) {
    axpy(1.0,x);
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) { 
    ::Thyra::scale(alpha, outArg(*thyra_vec_));
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  Real dot( const Vector<Real> &x ) const {
    Real val[1];
    const ThyraVector &ex = Teuchos::dyn_cast<const ThyraVector>(x);
    return ::Thyra::dot<Real>(*thyra_vec_, *ex.thyra_vec_);
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    return ::Thyra::norm_2<Real>(*thyra_vec_);
  } 

  /** \brief Clone to make a new (uninitialized) vector.
  */
  Teuchos::RCP<Vector<Real> > clone() const{
    Teuchos::RCP<Thyra::VectorBase<Real> > tv = thyra_vec_->clone_v();
    return Teuchos::rcp(new ThyraVector(tv)); 
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void axpy( const Real alpha, const Vector<Real> &x ) {
    const ThyraVector &ex = Teuchos::dyn_cast<const ThyraVector>(x);

  ::Thyra::linear_combination<Real>(
      Teuchos::tuple<Real>(alpha)(),
      Teuchos::tuple<Teuchos::Ptr<const ::Thyra::VectorBase<Real> > >(ex.getVector().ptr())(),
      1.0,
      outArg(*thyra_vec_)
      );
  }

  /**  \brief Set to zero vector.
  */
  void zero() {
    ::Thyra::put_scalar(0.0, outArg(*thyra_vec_));
  }

  void putScalar(Real alpha) {
      ::Thyra::put_scalar(alpha, outArg(*thyra_vec_));
    }

  Teuchos::RCP<const Thyra::VectorBase<Real> > getVector() const {
    return thyra_vec_;
  }

  Teuchos::RCP<Thyra::VectorBase<Real> > getNonConstVector()  {
    return thyra_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<Thyra::VectorBase<Real> > basisThyraVec = thyra_vec_->clone_v(); 
    ::Thyra::put_scalar(0.0, outArg(*basisThyraVec));
    ::Thyra::set_ele(i,1.0,outArg(*basisThyraVec));

    return Teuchos::rcp(new ThyraVector(basisThyraVec));
  }

  int dimension() const {return thyra_vec_->space()->dim();}


}; // class ThyraVector

} // namespace ROL

#endif


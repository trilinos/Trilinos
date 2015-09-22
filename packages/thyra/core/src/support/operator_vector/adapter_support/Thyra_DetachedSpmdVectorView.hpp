// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DETACHED_SPMD_VECTOR_VIEW_HPP
#define THYRA_DETACHED_SPMD_VECTOR_VIEW_HPP


#include "Thyra_SpmdVectorBase.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


/** \brief Create an explicit detached non-mutable (const) view of all of the
 * local elements on this process of an <tt>VectorBase</tt> object.
 *
 * The default constructor, copy constructor and assignment operators
 * are not allowed.
 *
 * \ingroup Thyra_Op_Vec_spmd_adapters_grp
 */
template<class Scalar>
class ConstDetachedSpmdVectorView {
public:
  /** \brief . */
  ConstDetachedSpmdVectorView(const Teuchos::RCP<const VectorBase<Scalar> > &v)
    {
      using Teuchos::rcp_dynamic_cast;
      if (!is_null(v)) {
        const RCP<const SpmdVectorBase<Scalar> > spmd_v =
          rcp_dynamic_cast<const SpmdVectorBase<Scalar> >(v, true);
        v_ = spmd_v;
        sv_ = spmd_v->getLocalSubVector();
      }
      else {
        v_ = Teuchos::null;
        sv_ = RTOpPack::ConstSubVectorView<Scalar>();
      }
    }
  /** \brief . */
  ~ConstDetachedSpmdVectorView()
    {}
  /** \brief . */
  const RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const
    { if (!is_null(v_)) return v_->spmdSpace(); return Teuchos::null; }
  /** \brief . */
  const RTOpPack::ConstSubVectorView<Scalar>& sv() const { return sv_; }
  /** \brief . */
  Teuchos_Ordinal globalOffset() const { return sv_.globalOffset(); }
  /** \brief . */
  Teuchos_Ordinal subDim() const { return sv_.subDim(); }
  /** \brief . */
  const ArrayRCP<const Scalar> values() const { return sv_.values(); }
  /** \brief . */
  ptrdiff_t stride() const { return sv_.stride(); }
  /** \brief . */
  const Scalar& operator[](Teuchos_Ordinal i) const { return sv_[i]; }
  /** \brief . */
  const Scalar& operator()(Teuchos_Ordinal i) const { return sv_(i); }
private:
  Teuchos::RCP<const SpmdVectorBase<Scalar> > v_;
  RTOpPack::ConstSubVectorView<Scalar>  sv_;
  // Not defined and not to be called
  ConstDetachedSpmdVectorView();
  ConstDetachedSpmdVectorView(const ConstDetachedSpmdVectorView<Scalar>&);
  ConstDetachedSpmdVectorView<Scalar>& operator==(
    const ConstDetachedSpmdVectorView<Scalar>&);
};


/** \brief Create an explicit detached mutable (non-const) view of all of the
 * local elements on this process of an <tt>VectorBase</tt> object.
 *
 * The default constructor, copy constructor and assignment operators
 * are not allowed.
 *
 * \ingroup Thyra_Op_Vec_spmd_adapters_grp
 */
template<class Scalar>
class DetachedSpmdVectorView {
public:
  /** \brief . */
  DetachedSpmdVectorView(const Teuchos::RCP<VectorBase<Scalar> > &v)
    {
      using Teuchos::rcp_dynamic_cast;
      if (!is_null(v)) {
        const RCP<SpmdVectorBase<Scalar> > spmd_v =
          rcp_dynamic_cast<SpmdVectorBase<Scalar> >(v, true);
        v_ = spmd_v;
        sv_ = spmd_v->getNonconstLocalSubVector();
      }
      else {
        v_ = Teuchos::null;
        sv_ = RTOpPack::SubVectorView<Scalar>();
      }
    }
  /** \brief . */
  ~DetachedSpmdVectorView()
    {}
  /** \brief . */
  const RCP<const SpmdVectorSpaceBase<Scalar> > spmdSpace() const
    { if (!is_null(v_)) return v_->spmdSpace(); return Teuchos::null; }
  /** \brief . */
  const RTOpPack::SubVectorView<Scalar>& sv() const { return sv_; }
  /** \brief . */
  Teuchos_Ordinal globalOffset() const { return sv_.globalOffset(); }
  /** \brief . */
  Teuchos_Ordinal subDim() const { return sv_.subDim(); }
  /** \brief . */
  const ArrayRCP<const Scalar> values() const { return sv_.values(); }
  /** \brief . */
  ptrdiff_t stride() const { return sv_.stride(); }
  /** \brief . */
  Scalar& operator[](Teuchos_Ordinal i) const { return sv_[i]; }
  /** \brief . */
  Scalar& operator()(Teuchos_Ordinal i) const { return sv_(i); }
private:
  Teuchos::RCP<SpmdVectorBase<Scalar> > v_;
  RTOpPack::SubVectorView<Scalar>  sv_;
  // Not defined and not to be called
  DetachedSpmdVectorView();
  DetachedSpmdVectorView(const DetachedSpmdVectorView<Scalar>&);
  DetachedSpmdVectorView<Scalar>& operator==(
    const DetachedSpmdVectorView<Scalar>&);
};


} // namespace Thyra


#endif // THYRA_DETACHED_SPMD_VECTOR_VIEW_HPP

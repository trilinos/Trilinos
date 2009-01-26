// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
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
        const RCP<const SpmdVectorSpaceBase<Scalar> > spmd_vs =
          spmd_v->spmdSpace();
        const Scalar *localValues = 0;
        Index l_stride = 0;
        spmd_v->getLocalData(&localValues, &l_stride);
        TEUCHOS_ASSERT_EQUALITY(l_stride, 1);
        const Index localSubDim = spmd_vs->localSubDim();
        sv_.initialize(
          spmd_vs->localOffset(), // globalOffset?
          localSubDim,
          Teuchos::arcp(localValues, 0, localSubDim, false),
          1 // stride
          );
      }
      else {
        v_ = Teuchos::null;
        sv_ = RTOpPack::ConstSubVectorView<Scalar>();
      }
    }
  /** \brief . */
  ~ConstDetachedSpmdVectorView()
    {
      if (!is_null(v_)) {
        v_->freeLocalData(sv_.values().get());
      }
    }
  /** \brief . */
  const RTOpPack::ConstSubVectorView<Scalar>& sv() const { return sv_; }
  /** \brief . */
  Teuchos_Index globalOffset() const { return sv_.globalOffset(); }
  /** \brief . */
  Teuchos_Index subDim() const { return sv_.subDim(); }
  /** \brief . */
  const ArrayRCP<const Scalar> values() const { return sv_.values(); }
  /** \brief . */
  ptrdiff_t stride() const { return sv_.stride(); }
  /** \brief . */
  const Scalar& operator[](Teuchos_Index i) const { return sv_[i]; }
  /** \brief . */
  const Scalar& operator()(Teuchos_Index i) const { return sv_(i); }
private:
  Teuchos::RCP<const SpmdVectorBase<Scalar> > v_;
  RTOpPack::ConstSubVectorView<Scalar>  sv_;
  // Not defined and not to be called
  ConstDetachedSpmdVectorView();
  ConstDetachedSpmdVectorView(const ConstDetachedSpmdVectorView<Scalar>&);
  ConstDetachedSpmdVectorView<Scalar>& operator==(
    const ConstDetachedSpmdVectorView<Scalar>&);
};


} // namespace Thyra


#endif // THYRA_DETACHED_SPMD_VECTOR_VIEW_HPP

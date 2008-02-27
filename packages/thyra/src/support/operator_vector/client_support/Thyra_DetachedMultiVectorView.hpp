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

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_AssertOp.hpp"

#ifndef THYRA_EXPLICIT_MULTI_VECTOR_VIEW_HPP
#define THYRA_EXPLICIT_MULTI_VECTOR_VIEW_HPP

namespace Thyra {

/** \brief Create an explicit non-mutable (const) view of a <tt>MultiVectorBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class ConstDetachedMultiVectorView {
public:
  /** \brief . */
  ConstDetachedMultiVectorView(
    const MultiVectorBase<Scalar>& mv, const Range1D &rowRng = Range1D(), const Range1D &colRng = Range1D()
    )
    : mv_(mv) { mv_.acquireDetachedView(rowRng,colRng,&smv_); }
  /** \brief . */
  ~ConstDetachedMultiVectorView() { mv_.releaseDetachedView(&smv_); }
  /** \brief . */
  const RTOpPack::ConstSubMultiVectorView<Scalar>& smv() const { return smv_; }
  /** \brief . */
  Teuchos_Index globalOffset() const { return smv_.globalOffset(); }
  /** \brief . */
  Teuchos_Index subDim() const { return smv_.subDim(); }
  /** \brief . */
  Teuchos_Index colOffset() const { return smv_.colOffset(); }
  /** \brief . */
  Teuchos_Index numSubCols() const { return smv_.numSubCols(); }
  /** \brief . */
  const Scalar* values() const { return smv_.values().get(); }
  /** \brief . */
  Teuchos_Index leadingDim() const { return smv_.leadingDim(); }
  /// Zero-based indexing: Preconditions: <tt>values()!=NULL && (0<=i<subDim()) && (0<=j<numSubCols())</tt>
  const Scalar& operator()(Teuchos_Index i,Teuchos_Index j) const { return smv_(i,j); }
private:
  const MultiVectorBase<Scalar> &mv_;
  RTOpPack::ConstSubMultiVectorView<Scalar> smv_;
  // Not defined and not to be called
  ConstDetachedMultiVectorView();
  ConstDetachedMultiVectorView(const ConstDetachedMultiVectorView<Scalar>&);
  ConstDetachedMultiVectorView<Scalar>& operator==(const ConstDetachedMultiVectorView<Scalar>&);
};
 
/** \brief Create an explicit mutable (non-const) view of a <tt>MultiVectorBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DetachedMultiVectorView {
public:
  /** \brief . */
  DetachedMultiVectorView(
    MultiVectorBase<Scalar>& mv, const Range1D &rowRng = Range1D(), const Range1D &colRng = Range1D()
    )
    : mv_(mv) { mv_.acquireDetachedView(rowRng,colRng,&smv_); }
  /** \brief . */
  ~DetachedMultiVectorView() { mv_.commitDetachedView(&smv_); }
  /** \brief . */
  const RTOpPack::SubMultiVectorView<Scalar>& smv() const { return smv_; }
  /** \brief . */
  Teuchos_Index globalOffset() const { return smv_.globalOffset(); }
  /** \brief . */
  Teuchos_Index subDim() const { return smv_.subDim(); }
  /** \brief . */
  Teuchos_Index colOffset() const { return smv_.colOffset(); }
  /** \brief . */
  Teuchos_Index numSubCols() const { return smv_.numSubCols(); }
  /** \brief . */
  Scalar* values() const { return smv_.values().get(); }
  /** \brief . */
  Teuchos_Index leadingDim() const { return smv_.leadingDim(); }
  /// Zero-based indexing: Preconditions: <tt>values()!=NULL && (0<=i<subDim()) && (0<=j<numSubCols())</tt>
  Scalar& operator()(Teuchos_Index i,Teuchos_Index j) { return smv_(i,j); }
private:
  MultiVectorBase<Scalar> &mv_;
  RTOpPack::SubMultiVectorView<Scalar> smv_;
  // Not defined and not to be called
  DetachedMultiVectorView();
  DetachedMultiVectorView(const DetachedMultiVectorView<Scalar>&);
  DetachedMultiVectorView<Scalar>& operator==(const DetachedMultiVectorView<Scalar>&);
};

/** \brief Do an explicit multi-vector adjoint.
 *
 * \relates DetachedMultiVectorView
 */
template<class Scalar>
void doExplicitMultiVectorAdjoint(
  const MultiVectorBase<Scalar>& mvIn, MultiVectorBase<Scalar>* mvTransOut
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==mvTransOut);
  THYRA_ASSERT_VEC_SPACES("doExplicitMultiVectorAdjoint(...)",
    *mvIn.domain(), *mvTransOut->range()
    );
  THYRA_ASSERT_VEC_SPACES("doExplicitMultiVectorAdjoint(...)",
    *mvIn.range(), *mvTransOut->domain()
    );
#endif
  ConstDetachedMultiVectorView<Scalar> dMvIn(mvIn);
  DetachedMultiVectorView<Scalar> dMvTransOut(*mvTransOut);
  const int m = dMvIn.subDim();
  const int n = dMvIn.numSubCols();
  for ( int j = 0; j < n; ++j ) {
    for ( int i = 0; i < m; ++i ) {
      dMvTransOut(j,i) = ST::conjugate(dMvIn(i,j));
    }
  }
}


} // namespace Thyra

#endif // THYRA_EXPLICIT_MULTI_VECTOR_VIEW_HPP

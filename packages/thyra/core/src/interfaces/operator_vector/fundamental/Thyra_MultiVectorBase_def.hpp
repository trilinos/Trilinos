// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_BASE_HPP
#define THYRA_MULTI_VECTOR_BASE_HPP

#include "Thyra_MultiVectorBase_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps_decl.hpp"

namespace Thyra {


// Provide access to the columns as VectorBase objects


template<class Scalar>
RCP<const VectorBase<Scalar> >
MultiVectorBase<Scalar>::colImpl(Ordinal j) const
{
  return const_cast<MultiVectorBase*>(this)->nonconstColImpl(j);
}


// Overridden methods from LinearOpBase


template<class Scalar>
RCP<const LinearOpBase<Scalar> >
MultiVectorBase<Scalar>::clone() const
{
  return this->clone_mv();
}

// Overridden methods from RowStatLinearOpBase

template<class Scalar>
bool MultiVectorBase<Scalar>::
rowStatIsSupportedImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat) const
{
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      return true;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  TEUCHOS_UNREACHABLE_RETURN(false);
}

template<class Scalar>
void MultiVectorBase<Scalar>::
getRowStatImpl(const RowStatLinearOpBaseUtils::ERowStat rowStat,
               const Ptr<VectorBase<Scalar> > &rowStatVec) const
{
  switch (rowStat) {
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM:
      absRowSum(rowStatVec);
      ::Thyra::reciprocal<Scalar>(*rowStatVec,rowStatVec.ptr());
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM:
      // compute absolute row sum
      absRowSum(rowStatVec);
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_INV_COL_SUM:
      absColSum(rowStatVec);
      ::Thyra::reciprocal<Scalar>(*rowStatVec,rowStatVec.ptr());
      break;
    case RowStatLinearOpBaseUtils::ROW_STAT_COL_SUM:
      // compute absolute row sum
      absColSum(rowStatVec);
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}

// Overridden methods from ScaledLinearOpBase

template<class Scalar>
bool MultiVectorBase<Scalar>::
supportsScaleLeftImpl() const
{
  return true;
}
 	
template<class Scalar>
bool MultiVectorBase<Scalar>::
supportsScaleRightImpl() const
{
  return true;
}
 	
template<class Scalar>
void MultiVectorBase<Scalar>::
scaleLeftImpl(const VectorBase< Scalar > &row_scaling)
{
  // loop over each column applying the row scaling
  for(Ordinal i=0;i<this->domain()->dim();i++)
    ::Thyra::ele_wise_scale<Scalar>(row_scaling,this->col(i).ptr());
}
 	
template<class Scalar>
void MultiVectorBase<Scalar>::
scaleRightImpl(const VectorBase< Scalar > &col_scaling)
{
  // this is probably incorrect if the domain is distrbuted
  // but if it is on every processor its probably fine...

  RTOpPack::SubVectorView<Scalar> view;
  col_scaling.acquireDetachedView(Thyra::Range1D(),&view);

  Teuchos::ArrayRCP<const Scalar> col_scaling_vec = view.values();

  // check to make sure things match up
  TEUCHOS_ASSERT(this->domain()->dim()==col_scaling_vec.size());

  for(Ordinal i=0;i<this->domain()->dim();i++)
    ::Thyra::scale<Scalar>(col_scaling_vec[i],this->col(i).ptr());
}

// helper methods

template<class Scalar>
void MultiVectorBase<Scalar>::
absRowSum(const Teuchos::Ptr<Thyra::VectorBase<Scalar> > & output) const
{
  using Teuchos::RCP;
  using Teuchos::ptrFromRef;
  using Teuchos::tuple;

  // compute absolute value of multi-vector
  RCP<MultiVectorBase<Scalar> > abs_mv = createMembers(this->range(),this->domain());
  for (Ordinal i = 0; i < abs_mv->domain()->dim(); ++i)
    abs_mv->col(i)->abs(*this->col(i));

  // compute sum over all rows
  RCP<VectorBase<Scalar> > ones = Thyra::createMember(this->domain());
  ::Thyra::put_scalar<Scalar>(Teuchos::ScalarTraits<Scalar>::one(),ones.ptr());
  ::Thyra::apply<Scalar>(*abs_mv,Thyra::NOTRANS,*ones,output);
}

template<class Scalar>
void MultiVectorBase<Scalar>::
absColSum(const Teuchos::Ptr<Thyra::VectorBase<Scalar> > & output) const
{ 
  using Teuchos::tuple; 
  using Teuchos::ptrInArg; 
  using Teuchos::null;
  using Teuchos::Array;
  using Teuchos::ArrayView;

  RTOpPack::SubVectorView<Scalar> view;
  output->acquireDetachedView(Thyra::Range1D(),&view);
  Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(view.values().size());
  this->norms_1(norms());
  for (Ordinal i = 0; i < norms.size(); ++i)
    view[i] = Teuchos::as<Scalar>(norms[i]);
  output->commitDetachedView(&view);
}


} // end namespace Thyra


#endif // THYRA_MULTI_VECTOR_BASE_HPP

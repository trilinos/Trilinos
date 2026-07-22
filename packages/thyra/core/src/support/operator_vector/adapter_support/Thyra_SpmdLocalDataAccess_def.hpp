// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SPMD_LOCAL_DATA_ACCESS_DEF_HPP
#define THYRA_SPMD_LOCAL_DATA_ACCESS_DEF_HPP


#include "Thyra_SpmdLocalDataAccess_decl.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"


template<class Scalar>
RTOpPack::SubVectorView<Scalar>
Thyra::getNonconstLocalSubVectorView(
  const RCP<VectorBase<Scalar> > &vec)
{
  const RCP<SpmdVectorBase<Scalar> > spmd_v = 
    Teuchos::rcp_dynamic_cast<SpmdVectorBase<Scalar> >(vec);
  if (nonnull(spmd_v)) {
    return spmd_v->getNonconstLocalSubVector();
  }
  const RCP<ProductVectorBase<Scalar> > p_v = 
    Teuchos::rcp_dynamic_cast<ProductVectorBase<Scalar> >(vec, true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(p_v->productSpace()->numBlocks() == 1), 
    std::logic_error,
    "Error, the function getNonconstLocalSubVectorView() can only return"
    " a contiguous view of local SPMD data from a product vector with a single"
    " block (which also must be able to give up a product view.");
  return getNonconstLocalSubVectorView<Scalar>(p_v->getNonconstVectorBlock(0));
}


template<class Scalar>
RTOpPack::ConstSubVectorView<Scalar>
Thyra::getLocalSubVectorView(
  const RCP<const VectorBase<Scalar> > &vec)
{
  const RCP<const SpmdVectorBase<Scalar> > spmd_v = 
    Teuchos::rcp_dynamic_cast<const SpmdVectorBase<Scalar> >(vec);
  if (nonnull(spmd_v)) {
    return spmd_v->getLocalSubVector();
  }
  const RCP<const ProductVectorBase<Scalar> > p_v = 
    Teuchos::rcp_dynamic_cast<const ProductVectorBase<Scalar> >(vec, true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(p_v->productSpace()->numBlocks() == 1), 
    std::logic_error,
    "Error, the function getLocalSubVectorView() can only return"
    " a contiguous view of local SPMD data from a product vector with a single"
    " block (which also must be able to give up a product view.");
  return getLocalSubVectorView<Scalar>(p_v->getVectorBlock(0));
}


template<class Scalar>
RTOpPack::SubMultiVectorView<Scalar>
Thyra::getNonconstLocalSubMultiVectorView(
  const RCP<MultiVectorBase<Scalar> > &multivec)
{
  const RCP<SpmdMultiVectorBase<Scalar> > spmd_mv = 
    Teuchos::rcp_dynamic_cast<SpmdMultiVectorBase<Scalar> >(multivec);
  if (nonnull(spmd_mv)) {
    return spmd_mv->getNonconstLocalSubMultiVector();
  }
  const RCP<ProductMultiVectorBase<Scalar> > p_mv = 
    Teuchos::rcp_dynamic_cast<ProductMultiVectorBase<Scalar> >(multivec, true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(p_mv->productSpace()->numBlocks() == 1), 
    std::logic_error,
    "Error, the function getNonconstLocalSubVectorView() can only return"
    " a contiguous view of local SPMD data from a product vector with a single"
    " block (which also must be able to give up a product view.");
  return getNonconstLocalSubMultiVectorView<Scalar>(p_mv->getNonconstMultiVectorBlock(0));
}


template<class Scalar>
RTOpPack::ConstSubMultiVectorView<Scalar>
Thyra::getLocalSubMultiVectorView(
  const RCP<const MultiVectorBase<Scalar> > &multivec)
{
  const RCP<const SpmdMultiVectorBase<Scalar> > spmd_mv = 
    Teuchos::rcp_dynamic_cast<const SpmdMultiVectorBase<Scalar> >(multivec);
  if (nonnull(spmd_mv)) {
    return spmd_mv->getLocalSubMultiVector();
  }
  const RCP<const ProductMultiVectorBase<Scalar> > p_mv = 
    Teuchos::rcp_dynamic_cast<const ProductMultiVectorBase<Scalar> >(multivec, true);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(p_mv->productSpace()->numBlocks() == 1), 
    std::logic_error,
    "Error, the function getLocalSubVectorView() can only return"
    " a contiguous view of local SPMD data from a product vector with a single"
    " block (which also must be able to give up a product view.");
  return getLocalSubMultiVectorView<Scalar>(p_mv->getMultiVectorBlock(0));
}


//
// Explicit instantiation macro
//


#define THYRA_SPMD_LOCAL_DATA_ACCESS_INSTANT(SCALAR) \
   \
  template RTOpPack::SubVectorView<SCALAR> \
  getNonconstLocalSubVectorView(const RCP<VectorBase<SCALAR> > &vec); \
  \
  template RTOpPack::ConstSubVectorView<SCALAR> \
  getLocalSubVectorView(const RCP<const VectorBase<SCALAR> > &vec); \
  \
  template RTOpPack::SubMultiVectorView<SCALAR> \
  getNonconstLocalSubMultiVectorView(const RCP<MultiVectorBase<SCALAR> > &vec); \
  \
  template RTOpPack::ConstSubMultiVectorView<SCALAR> \
  getLocalSubMultiVectorView(const RCP<const MultiVectorBase<SCALAR> > &vec); \


#endif // THYRA_SPMD_LOCAL_DATA_ACCESS_DEF_HPP

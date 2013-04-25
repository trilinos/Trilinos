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
    "Error, the function getNonconstLocalSubVectorView() can only return"
    " a contiguous view of local SPMD data from a product vector with a single"
    " block (which also must be able to give up a product view.");
  return getLocalSubVectorView<Scalar>(p_v->getVectorBlock(0));
}


template<class Scalar>
RTOpPack::SubMultiVectorView<Scalar>
Thyra::getNonconstLocalSubMultiVectorView(
  const RCP<MultiVectorBase<Scalar> > &vec)
{
  return Teuchos::rcp_dynamic_cast<SpmdMultiVectorBase<Scalar> >(vec, true)
    ->getNonconstLocalSubMultiVector();
}


template<class Scalar>
RTOpPack::ConstSubMultiVectorView<Scalar>
Thyra::getLocalSubMultiVectorView(
  const RCP<const MultiVectorBase<Scalar> > &vec)
{
  return Teuchos::rcp_dynamic_cast<const SpmdMultiVectorBase<Scalar> >(vec, true)
    ->getLocalSubMultiVector();
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

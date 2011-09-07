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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_COLUMNWISE_MULTI_VECTOR_DEF_HPP
#define THYRA_DEFAULT_COLUMNWISE_MULTI_VECTOR_DEF_HPP

#include "Thyra_DefaultColumnwiseMultiVector_decl.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_as.hpp"

namespace Thyra {


// Constructors/Initializers


template<class Scalar>
DefaultColumnwiseMultiVector<Scalar>::DefaultColumnwiseMultiVector()
{}


template<class Scalar>
DefaultColumnwiseMultiVector<Scalar>::DefaultColumnwiseMultiVector(
  const RCP<VectorBase<Scalar> > &col_vec
  )
{
  this->initialize(col_vec);
}


template<class Scalar>
DefaultColumnwiseMultiVector<Scalar>::DefaultColumnwiseMultiVector(
  const RCP<const VectorSpaceBase<Scalar> > &range_in,
  const RCP<const VectorSpaceBase<Scalar> > &domain_in,
  const ArrayView<const RCP<VectorBase<Scalar> > > &col_vecs_in
  )
{
  this->initialize(range_in, domain_in, col_vecs_in);
}


template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::initialize(
  const RCP<VectorBase<Scalar> > &col_vec
  )
{
#ifdef TEUCHOS_DEBUG
  const std::string err_msg =
    "DefaultColumnwiseMultiVector<Scalar>::initialize(...): Error!";
  TEST_FOR_EXCEPT_MSG( is_null(col_vec), err_msg ); 
  TEST_FOR_EXCEPT_MSG( is_null(col_vec->space()), err_msg ); 
#endif
  range_  = col_vec->space();
  domain_ = range_->smallVecSpcFcty()->createVecSpc(1);
  col_vecs_.resize(1);
  col_vecs_[0] = col_vec;
}

  
template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::initialize(
  const RCP<const VectorSpaceBase<Scalar> > &range_in,
  const RCP<const VectorSpaceBase<Scalar> > &domain_in,
  const ArrayView<const RCP<VectorBase<Scalar> > > &col_vecs
  )
{
#ifdef TEUCHOS_DEBUG
  const std::string err_msg =
    "DefaultColumnwiseMultiVector<Scalar>::initialize(...): Error!";
  TEST_FOR_EXCEPT_MSG( is_null(range_in), err_msg ); 
  TEST_FOR_EXCEPT_MSG( is_null(domain_in), err_msg ); 
  TEST_FOR_EXCEPT_MSG( range_in->dim()  == 0, err_msg ); 
  TEST_FOR_EXCEPT_MSG( domain_in->dim() == 0, err_msg );
  // ToDo: Check the compatibility of the vectors in col_vecs!
#endif
  const int domainDim = domain_in->dim();
  range_ = range_in;
  domain_ = domain_in;
  col_vecs_.clear();
  col_vecs_.reserve(domainDim);
  if (col_vecs.size()) {
    for( Ordinal j = 0; j < domainDim; ++j )
      col_vecs_.push_back(col_vecs[j]);
  }
  else {
    for( Ordinal j = 0; j < domainDim; ++j )
      col_vecs_.push_back(createMember(range_));
  }
}


template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::uninitialize()
{
  col_vecs_.resize(0);
  range_ = Teuchos::null;
  domain_ = Teuchos::null;
}


// Overridden from OpBase


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::range() const
{
  return range_;
}


template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::domain() const
{
  return domain_;
}


// Overridden from LinearOpBase


template<class Scalar>
bool DefaultColumnwiseMultiVector<Scalar>::opSupportedImpl(EOpTransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ( ST::isComplex ? ( M_trans==NOTRANS || M_trans==CONJTRANS ) : true );
}


template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::applyImpl(
  const EOpTransp M_trans,
  const MultiVectorBase<Scalar> &X,
  const Ptr<MultiVectorBase<Scalar> > &Y,
  const Scalar alpha,
  const Scalar beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "MultiVectorBase<Scalar>::apply()", *this, M_trans, X, &*Y);
#endif
  const Ordinal nc = this->domain()->dim();
  const Ordinal m = X.domain()->dim();
  for (Ordinal col_j = 0; col_j < m; ++col_j) {
    const RCP<const VectorBase<Scalar> > x_j = X.col(col_j);
    const RCP<VectorBase<Scalar> > y_j = Y->col(col_j);
    // y_j *= beta
    Vt_S(y_j.ptr(), beta);
    // y_j += alpha*op(M)*x_j
    if(M_trans == NOTRANS) {
      //
      // y_j += alpha*M*x_j = alpha*M.col(0)*x_j(0) + ... + alpha*M.col(nc-1)*x_j(nc-1)
      //
      // Extract an explicit view of x_j
      RTOpPack::ConstSubVectorView<Scalar> x_sub_vec;               
      x_j->acquireDetachedView(Range1D(), &x_sub_vec);
      // Loop through and add the multiple of each column
      for (Ordinal j = 0; j < nc; ++j )
        Vp_StV( y_j.ptr(), Scalar(alpha*x_sub_vec(j)), *this->col(j) );
      // Release the view of x
      x_j->releaseDetachedView(&x_sub_vec);
    }
    else {
      //
      //                        [ alpha*dot(M.col(0),x_j)    ]
      // y_j += alpha*M^T*x_j = [ alpha*dot(M.col(1),x_j)    ]
      //                        [ ...                        ]
      //                        [ alpha*dot(M.col(nc-1),x_j) ]
      //
      // Extract an explicit view of y_j
      RTOpPack::SubVectorView<Scalar> y_sub_vec;               
      y_j->acquireDetachedView(Range1D(), &y_sub_vec);
      // Loop through and add to each element in y_j
      for (Ordinal j = 0; j < nc; ++j )
        y_sub_vec(j) += alpha*dot(*this->col(j), *x_j);
      // Commit explicit view of y
      y_j->commitDetachedView(&y_sub_vec);
    }
  }
}


// Overridden from MultiVectorBase


template<class Scalar>
RCP<VectorBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::nonconstColImpl(Ordinal j)
{
  using Teuchos::as;
  const int num_cols = col_vecs_.size();
  TEST_FOR_EXCEPTION(
    !(  0 <= j  && j < num_cols ), std::logic_error
    ,"Error, j = " << j << " does not fall in the range [0,"<<(num_cols-1)<< "]!"
    );
  return col_vecs_[j];
}


template<class Scalar>
RCP<MultiVectorBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::nonconstContigSubViewImpl(
  const Range1D& col_rng_in
  )
{
  const Ordinal numCols = domain_->dim();
  const Range1D col_rng = Teuchos::full_range(col_rng_in,0,numCols-1);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !( col_rng.ubound() < numCols ), std::logic_error
    ,"DefaultColumnwiseMultiVector<Scalar>::subView(col_rng):"
    "Error, the input range col_rng = ["<<col_rng.lbound()
    <<","<<col_rng.ubound()<<"] "
    "is not in the range [0,"<<(numCols-1)<<"]!"
    );
#endif
  return Teuchos::rcp(
    new DefaultColumnwiseMultiVector<Scalar>(
      range_,
      domain_->smallVecSpcFcty()->createVecSpc(col_rng.size()),
      col_vecs_(col_rng.lbound(),col_rng.size())
      )
    );
}

  
} // end namespace Thyra


#endif // THYRA_DEFAULT_COLUMNWISE_MULTI_VECTOR_DEF_HPP

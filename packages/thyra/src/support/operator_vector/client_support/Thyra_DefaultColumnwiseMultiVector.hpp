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

#ifndef THYRA_MULTI_VECTOR_COLS_HPP
#define THYRA_MULTI_VECTOR_COLS_HPP

#include "Thyra_DefaultColumnwiseMultiVectorDecl.hpp"
#include "Thyra_MultiVectorDefaultBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_SingleRhsLinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Teuchos_TestForException.hpp"

namespace Thyra {

// Constructors/Initializers

template<class Scalar>
DefaultColumnwiseMultiVector<Scalar>::DefaultColumnwiseMultiVector()
{}

template<class Scalar>
DefaultColumnwiseMultiVector<Scalar>::DefaultColumnwiseMultiVector(
  const Teuchos::RefCountPtr<VectorBase<Scalar> > &col_vec
  )
{
  this->initialize(col_vec);
}

template<class Scalar>
DefaultColumnwiseMultiVector<Scalar>::DefaultColumnwiseMultiVector(
  const  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >   &range
  ,const  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &domain
  ,const Teuchos::RefCountPtr<VectorBase<Scalar> >              col_vecs[]
  )
{
  this->initialize(range,domain,col_vecs);
}

template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<VectorBase<Scalar> > &col_vec
  )
{
#ifdef TEUCHOS_DEBUG
  const char err_msg[] = "DefaultColumnwiseMultiVector<Scalar>::initialize(...): Error!";
  TEST_FOR_EXCEPTION( col_vec.get() == NULL,           std::invalid_argument, err_msg ); 
  TEST_FOR_EXCEPTION( col_vec->space().get() == NULL,  std::invalid_argument, err_msg ); 
#endif
  range_  = col_vec->space();
  domain_ = range_->smallVecSpcFcty()->createVecSpc(1);
  num_cols_ = 1;
  col_vecs_.resize(1);
  col_vecs_[0] = col_vec;
}
  
template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::initialize(
  const  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >   &range
  ,const  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &domain
  ,const Teuchos::RefCountPtr<VectorBase<Scalar> >              col_vecs[]
  )
{
#ifdef TEUCHOS_DEBUG
  const char err_msg[] = "DefaultColumnwiseMultiVector<Scalar>::initialize(...): Error!";
  TEST_FOR_EXCEPTION( range.get()   == NULL, std::invalid_argument, err_msg ); 
  TEST_FOR_EXCEPTION( domain.get()  == NULL, std::invalid_argument, err_msg ); 
  TEST_FOR_EXCEPTION( range->dim()  == 0,    std::invalid_argument, err_msg ); 
  TEST_FOR_EXCEPTION( domain->dim() == 0,    std::invalid_argument, err_msg );
  // ToDo: Check the compatibility of the vectors in col_vecs!
#endif
  range_ = range;
  domain_ = domain;
  num_cols_ = domain->dim();
  col_vecs_.resize(num_cols_);
  if(col_vecs) {
    for( Index j = 0; j < num_cols_; ++j )
      col_vecs_[j] = col_vecs[j];
  }
  else {
    for( Index j = 0; j < num_cols_; ++j )
      col_vecs_[j] = createMember(range_);
  }
}

template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::set_uninitialized()
{
  col_vecs_.resize(0);
  range_ = Teuchos::null;
  domain_ = Teuchos::null;
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::range() const
{
  return range_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::domain() const
{
  return domain_;
}

// Overridden from SingleRhsLinearOpBase

template<class Scalar>
bool DefaultColumnwiseMultiVector<Scalar>::opSupported(ETransp M_trans) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return ( ST::isComplex ? ( M_trans==NOTRANS || M_trans==CONJTRANS ) : true );
}

template<class Scalar>
void DefaultColumnwiseMultiVector<Scalar>::apply(
  const ETransp                M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  ,const Scalar                beta
  ) const
{
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("MultiVectorBase<Scalar>::apply()",*this,M_trans,x,y);
#endif
  const Index nc = this->domain()->dim();
  // y *= beta
  Vt_S(y,beta);
  // y += alpha*op(M)*x
  if(M_trans == NOTRANS) {
    //
    // y += alpha*M*x = alpha*M.col(0)*x(0) + ... + alpha*M.col(nc-1)*x(nc-1)
    //
    // Extract an explicit view of x
    RTOpPack::ConstSubVectorView<Scalar> x_sub_vec;               
    x.acquireDetachedView(Range1D(),&x_sub_vec);
    // Loop through and add the multiple of each column
    for(Index j = 0; j < nc; ++j )
      Vp_StV( y, Scalar(alpha*x_sub_vec(j)), *this->col(j) );
    // Release the view of x
    x.releaseDetachedView(&x_sub_vec);
  }
  else {
    //
    //                   [ alpha*dot(M.col(0),x)    ]
    // y += alpha*M'*x = [ alpha*dot(M.col(1),x)    ]
    //                   [ ...                      ]
    //                   [ alpha*dot(M.col(nc-1),x) ]
    //
    // Extract an explicit view of y
    RTOpPack::SubVectorView<Scalar> y_sub_vec;               
    y->acquireDetachedView(Range1D(),&y_sub_vec);
    // Loop through and add to each element in y
    for(Index j = 0; j < nc; ++j )
      y_sub_vec(j) += alpha*dot(*this->col(j),x);
    // Commit explicit view of y
    y->commitDetachedView(&y_sub_vec);
  }
}

// Overridden from MultiVectorBase

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::col(Index j)
{
  TEST_FOR_EXCEPTION(
    !(  0 <= j  && j < num_cols_ ), std::logic_error
    ,"Error, j = " << j << " does not fall in the range [0,"<<(num_cols_-1)<< "]!"
    );
  return col_vecs_[j];
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultColumnwiseMultiVector<Scalar>::subView( const Range1D& col_rng_in )
{
  const Index numCols = domain_->dim();
  const Range1D col_rng = Teuchos::full_range(col_rng_in,0,numCols-1);
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !( col_rng.ubound() < numCols ), std::logic_error
    ,"DefaultColumnwiseMultiVector<Scalar>::subView(col_rng): Error, the input range col_rng = ["<<col_rng.lbound()<<","<<col_rng.ubound()<<"] "
    "is not in the range [0,"<<(numCols-1)<<"]!"
    );
#endif
  return Teuchos::rcp(
    new DefaultColumnwiseMultiVector<Scalar>(
      range_
      ,domain_->smallVecSpcFcty()->createVecSpc(col_rng.size())
      ,&col_vecs_[col_rng.lbound()]
      )
    );
}
  
} // end namespace Thyra

#endif // THYRA_MULTI_VECTOR_COLS_HPP

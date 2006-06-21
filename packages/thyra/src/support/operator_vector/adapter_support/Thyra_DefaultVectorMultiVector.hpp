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

#ifndef THYRA_VECTOR_MULTI_VECTOR_HPP
#define THYRA_VECTOR_MULTI_VECTOR_HPP

// Define to make some verbose output
//#define THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT

#include "Thyra_DefaultVectorMultiVectorDecl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Teuchos_Workspace.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DefaultVectorMultiVector<Scalar>::DefaultVectorMultiVector()
{}

template<class Scalar>
DefaultVectorMultiVector<Scalar>::DefaultVectorMultiVector(
  const Teuchos::RefCountPtr<MultiVectorBase<Scalar> > &mv
  )
{
  initialize(mv);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::initialize(
  const Teuchos::RefCountPtr<MultiVectorBase<Scalar> > &mv
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION( mv.get() == NULL, std::invalid_argument, "Error!" );
  TEST_FOR_EXCEPTION( mv->domain().get() == NULL, std::invalid_argument, "Error!" );
  TEST_FOR_EXCEPTION( mv->domain()->dim() != 1, std::invalid_argument, "Error!" );
#endif
  mv_ = mv;
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::uninitialize(
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > *mv
  )
{
  if(mv) *mv = mv_;
  mv_ = Teuchos::null;
}

// Overridden from LinearOpBase (forwarded to this->mv())

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultVectorMultiVector<Scalar>::range() const
{
  return (mv_.get() ? mv_->range() : Teuchos::null );
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultVectorMultiVector<Scalar>::domain() const
{
  return (mv_.get() ? mv_->domain() : Teuchos::null );
}

template<class Scalar>
bool DefaultVectorMultiVector<Scalar>::opSupported(ETransp M_trans) const
{
  return mv_->opSupported(M_trans);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  mv_->apply(M_trans,X,Y,alpha,beta);
}

// Overridden from MultiVectorBase (forwarded to this->mv())

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::col(Index j)
{
  return mv_->col(j);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::clone_mv() const
{
  return mv_->clone_mv();
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::subView( const Range1D& col_rng ) const
{
  return mv_->subView(col_rng);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::subView( const Range1D& col_rng )
{
  return mv_->subView(col_rng);
}

template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::subView( const int numCols, const int cols[] ) const
{
  return mv_->subView(numCols,cols);
}

template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
DefaultVectorMultiVector<Scalar>::subView( const int numCols, const int cols[] )
{
  return mv_->subView(numCols,cols);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>   &primary_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget*        reduct_objs[]
  ,const Index                    primary_first_ele
  ,const Index                    primary_sub_dim
  ,const Index                    primary_global_offset
  ,const Index                    secondary_first_ele
  ,const Index                    secondary_sub_dim
  ) const
{
  mv_->applyOp(
    primary_op,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs,reduct_objs
    ,primary_first_ele,primary_sub_dim,primary_global_offset,secondary_first_ele,secondary_sub_dim
    );
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>   &primary_op
  ,const RTOpPack::RTOpT<Scalar>  &secondary_op
  ,const int                      num_multi_vecs
  ,const MultiVectorBase<Scalar>* multi_vecs[]
  ,const int                      num_targ_multi_vecs
  ,MultiVectorBase<Scalar>*       targ_multi_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const Index                    primary_first_ele
  ,const Index                    primary_sub_dim
  ,const Index                    primary_global_offset
  ,const Index                    secondary_first_ele
  ,const Index                    secondary_sub_dim
  ) const
{
  mv_->applyOp(
    primary_op,secondary_op,num_multi_vecs,multi_vecs,num_targ_multi_vecs,targ_multi_vecs,reduct_obj
    ,primary_first_ele,primary_sub_dim,primary_global_offset,secondary_first_ele,secondary_sub_dim
    );
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::acquireDetachedView(
  const Range1D                       &rowRng
  ,const Range1D                      &colRng
  ,RTOpPack::ConstSubMultiVectorView<Scalar>  *sub_mv
  ) const
{
  mv_->acquireDetachedView(rowRng,colRng,sub_mv);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::releaseDetachedView( RTOpPack::ConstSubMultiVectorView<Scalar>* sub_mv ) const
{
  mv_->releaseDetachedView(sub_mv);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::acquireDetachedView(
  const Range1D                                &rowRng
  ,const Range1D                               &colRng
  ,RTOpPack::SubMultiVectorView<Scalar>    *sub_mv
  )
{
  mv_->acquireDetachedView(rowRng,colRng,sub_mv);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::commitDetachedView( RTOpPack::SubMultiVectorView<Scalar>* sub_mv )
{
  mv_->commitDetachedView(sub_mv);
}

// Overridden from VectorBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultVectorMultiVector<Scalar>::space() const
{
#ifdef THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVectorMultiVector<Scalar>::space() const called!\n";
#endif
  return ( mv_.get() ? mv_->range() : Teuchos::null );
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::applyOp(
  const RTOpPack::RTOpT<Scalar>   &op
  ,const int                      num_vecs
  ,const VectorBase<Scalar>*      vecs[]
  ,const int                      num_targ_vecs
  ,VectorBase<Scalar>*            targ_vecs[]
  ,RTOpPack::ReductTarget         *reduct_obj
  ,const Index                    first_ele
  ,const Index                    sub_dim
  ,const Index                    global_offset
  ) const
{
#ifdef THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVectorMultiVector<Scalar>::applyOp() const called!\n";
#endif
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  // Get MultiVectorBase arguments
  Workspace<const MultiVectorBase<Scalar>*> multi_vecs(wss,num_vecs,false);
  for( int k = 0; k < num_vecs; ++k ) {
    const DefaultVectorMultiVector<Scalar> *v_mv = dynamic_cast<const DefaultVectorMultiVector<Scalar>*>(vecs[k]);
    multi_vecs[k] = ( v_mv ? &*v_mv->mv() : vecs[k] );
  }
  Workspace<MultiVectorBase<Scalar>*> targ_multi_vecs(wss,num_targ_vecs,false);
  for( int k = 0; k < num_targ_vecs; ++k ) {
    DefaultVectorMultiVector<Scalar> *v_mv = dynamic_cast<DefaultVectorMultiVector<Scalar>*>(targ_vecs[k]);
    targ_multi_vecs[k] = ( v_mv ? &*v_mv->mv() : targ_vecs[k] );
  }
  RTOpPack::ReductTarget* reduct_objs[] = { reduct_obj };
  Thyra::applyOp(
    op
    ,num_vecs,num_vecs?&multi_vecs[0]:NULL
    ,num_targ_vecs,num_targ_vecs?&targ_multi_vecs[0]:NULL
    ,reduct_objs
    ,first_ele,sub_dim,global_offset
    ,1,0
    );
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::acquireDetachedView( const Range1D& rng, RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const
{
#ifdef THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVectorMultiVector<Scalar>::acquireDetachedView() const called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_vec==NULL);
#endif
  RTOpPack::ConstSubMultiVectorView<Scalar> sub_mv;
  mv_->acquireDetachedView(rng,Range1D(),&sub_mv);
  sub_vec->initialize(sub_mv.globalOffset(),sub_mv.subDim(),sub_mv.values(),1);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::releaseDetachedView( RTOpPack::ConstSubVectorView<Scalar>* sub_vec ) const
{
#ifdef THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVectorMultiVector<Scalar>::releaseDetachedView() const called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_vec==NULL);
#endif
  RTOpPack::ConstSubMultiVectorView<Scalar> sub_mv(sub_vec->globalOffset(),sub_vec->subDim(),0,1,sub_vec->values(),sub_vec->subDim());
  mv_->releaseDetachedView(&sub_mv);
  sub_vec->set_uninitialized();
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::acquireDetachedView( const Range1D& rng, RTOpPack::SubVectorView<Scalar>* sub_vec )
{
#ifdef THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVectorMultiVector<Scalar>::acquireDetachedView() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_vec==NULL);
#endif
  RTOpPack::SubMultiVectorView<Scalar> sub_mv;
  mv_->acquireDetachedView(rng,Range1D(),&sub_mv);
  sub_vec->initialize(sub_mv.globalOffset(),sub_mv.subDim(),sub_mv.values(),1);
}

template<class Scalar>
void DefaultVectorMultiVector<Scalar>::commitDetachedView( RTOpPack::SubVectorView<Scalar>* sub_vec )
{
#ifdef THYRA_VECTOR_MULTI_VECTOR_VERBOSE_TO_ERROR_OUT
  std::cerr << "\nVectorMultiVector<Scalar>::commitDetachedView() called!\n";
#endif
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(sub_vec==NULL);
#endif
  RTOpPack::SubMultiVectorView<Scalar> sub_mv(sub_vec->globalOffset(),sub_vec->subDim(),0,1,sub_vec->values(),sub_vec->subDim());
  mv_->commitDetachedView(&sub_mv);
  sub_vec->set_uninitialized();
}

} // end namespace Thyra

#endif // THYRA_VECTOR_MULTI_VECTOR_HPP

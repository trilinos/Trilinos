// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_LINEAR_COMBINATION_DEF_HPP
#define RTOPPACK_TOP_LINEAR_COMBINATION_DEF_HPP


#include "Teuchos_Workspace.hpp"


namespace RTOpPack {


template<class Scalar>
TOpLinearCombination<Scalar>::TOpLinearCombination(
  const ArrayView<const Scalar> &alpha_in,
  const Scalar &beta_in
  )
  :beta_(beta_in)
{
  if (alpha_in.size())
    this->alpha(alpha_in);
  this->setOpNameBase("TOpLinearCombination");
}



template<class Scalar>
void TOpLinearCombination<Scalar>::alpha(
  const ArrayView<const Scalar> &alpha_in )
{
  TEUCHOS_TEST_FOR_EXCEPT( alpha_in.size() == 0 );
  alpha_ = alpha_in;
}


template<class Scalar>
const ArrayView<const Scalar>
TOpLinearCombination<Scalar>::alpha() const
{ return alpha_; }


template<class Scalar>
void TOpLinearCombination<Scalar>::beta( const Scalar& beta_in ) { beta_ = beta_in; }


template<class Scalar>
Scalar TOpLinearCombination<Scalar>::beta() const { return beta_; }


template<class Scalar>
int TOpLinearCombination<Scalar>::num_vecs() const { return alpha_.size(); }


// Overridden from RTOpT


template<class Scalar>
void TOpLinearCombination<Scalar>::apply_op_impl(
  const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
  const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
  const Ptr<ReductTarget> &reduct_obj_inout
  ) const
{

  using Teuchos::as;
  using Teuchos::Workspace;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;
  typedef typename Teuchos::ArrayRCP<const Scalar>::iterator const_iter_t;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

#ifdef TEUCHOS_DEBUG
  validate_apply_op<Scalar>(*this, as<int>(alpha_.size()), 1, false,
    sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
  (void)reduct_obj_inout;
#endif

  const int l_num_vecs = alpha_.size();

  // Get iterators to local data
  const RTOpPack::index_type subDim = targ_sub_vecs[0].subDim();
  iter_t z0_val = targ_sub_vecs[0].values().begin();
  const ptrdiff_t z0_s = targ_sub_vecs[0].stride();
  Workspace<const_iter_t> v_val(wss,l_num_vecs);
  Workspace<ptrdiff_t> v_s(wss,l_num_vecs,false);
  for( int k = 0; k < l_num_vecs; ++k ) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT( sub_vecs[k].subDim() != subDim );
    TEUCHOS_TEST_FOR_EXCEPT( sub_vecs[k].globalOffset() != targ_sub_vecs[0].globalOffset() );
#endif					
    v_val[k] = sub_vecs[k].values().begin();
    v_s[k] = sub_vecs[k].stride();
  }

  //
  // Perform the operation and specialize the cases for l_num_vecs = 1 and 2
  // in order to get good performance.
  //
  if( l_num_vecs == 1 ) {
    //
    // z0 = alpha*v0 + beta*z0
    //
    const Scalar l_alpha = alpha_[0], l_beta = beta_;
    const_iter_t v0_val = v_val[0];
    const ptrdiff_t v0_s = v_s[0]; 
    if( l_beta==ST::zero() ) {
      // z0 = alpha*v0
      if( z0_s==1 && v0_s==1 ) {
        for( int j = 0; j < subDim; ++j )
          (*z0_val++) = l_alpha * (*v0_val++);
      }
      else {
        for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s )
          (*z0_val) = l_alpha * (*v0_val);
      }
    }
    else if( l_beta==ST::one() ) {
      //
      // z0 = alpha*v0 + z0
      //
      if( z0_s==1 && v0_s==1 ) {
        for( int j = 0; j < subDim; ++j )
          (*z0_val++) += l_alpha * (*v0_val++);
      }
      else {
        for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s )
          (*z0_val) += l_alpha * (*v0_val);
      }
    }
    else {
      // z0 = alpha*v0 + beta*z0
      if( z0_s==1 && v0_s==1 ) {
        for( int j = 0; j < subDim; ++j, ++z0_val )
          (*z0_val) = l_alpha * (*v0_val++) + l_beta*(*z0_val);
      }
      else {
        for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s )
          (*z0_val) = l_alpha * (*v0_val) + l_beta*(*z0_val);
      }
    }
  }
  else if( l_num_vecs == 2 ) {
    //
    // z0 = alpha0*v0 + alpha1*v1 + beta*z0
    //
    const Scalar alpha0 = alpha_[0], alpha1=alpha_[1], l_beta = beta_;
    const_iter_t v0_val = v_val[0];
    const ptrdiff_t v0_s = v_s[0]; 
    const_iter_t v1_val = v_val[1];
    const ptrdiff_t v1_s = v_s[1]; 
    if( l_beta==ST::zero() ) {
      if( alpha0 == ST::one() ) {
        if( alpha1 == ST::one() ) {
          // z0 = v0 + v1
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j )
              (*z0_val++) = (*v0_val++) + (*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = (*v0_val) + (*v1_val);
          }
        }
        else {
          // z0 = v0 + alpha1*v1
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j )
              (*z0_val++) = (*v0_val++) + alpha1*(*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = (*v0_val) + alpha1*(*v1_val);
          }
        }
      }
      else {
        if( alpha1 == ST::one() ) {
          // z0 = alpha0*v0 + v1
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j )
              (*z0_val++) = alpha0*(*v0_val++) + (*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = alpha0*(*v0_val) + (*v1_val);
          }
        }
        else {
          // z0 = alpha0*v0 + alpha1*v1
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j )
              (*z0_val++) = alpha0*(*v0_val++) + alpha1*(*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = alpha0*(*v0_val) + alpha1*(*v1_val);
          }
        }
      }
    }
    else if( l_beta==ST::one() ) {
      if( alpha0 == ST::one() ) {
        if( alpha1 == ST::one() ) {
          // z0 = v0 + v1 + z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) += (*v0_val++) + (*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) += (*v0_val) + (*v1_val);
          }
        }
        else {
          // z0 = v0 + alpha1*v1 + z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) += (*v0_val++) + alpha1*(*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) += (*v0_val) + alpha1*(*v1_val);
          }
        }
      }
      else {
        if( alpha1 == ST::one() ) {
          // z0 = alpha0*v0 + v1 + z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) += alpha0*(*v0_val++) + (*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) += alpha0*(*v0_val) + (*v1_val);
          }
        }
        else {
          // z0 = alpha0*v0 + alpha1*v1 + z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) += alpha0*(*v0_val++) + alpha1*(*v1_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) += alpha0*(*v0_val) + alpha1*(*v1_val);
          }
        }
      }
    }
    else {
      if( alpha0 == ST::one() ) {
        if( alpha1 == ST::one() ) {
          // z0 = v0 + v1 + beta*z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) = (*v0_val++) + (*v1_val++) + l_beta*(*z0_val);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = (*v0_val) + (*v1_val) + l_beta*(*z0_val);
          }
        }
        else {
          // z0 = v0 + alpha1*v1 + beta*z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) = (*v0_val++) + alpha1*(*v1_val++) + l_beta*(*z0_val);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = (*v0_val) + alpha1*(*v1_val) + l_beta*(*z0_val);
          }
        }
      }
      else {
        if( alpha1 == ST::one() ) {
          // z0 = alpha0*v0 + v1 + beta*z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) = alpha0*(*v0_val++) + (*v1_val++) + l_beta*(*z0_val);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = alpha0*(*v0_val) + (*v1_val) + l_beta*(*z0_val);
          }
        }
        else {
          // z0 = alpha0*v0 + alpha1*v1 + beta*z0
          if( z0_s==1 && v0_s==1 && v1_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) = alpha0*(*v0_val++) + alpha1*(*v1_val++) + l_beta*(*z0_val);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
              (*z0_val) = alpha0*(*v0_val) + alpha1*(*v1_val) + l_beta*(*z0_val);
          }
        }
      }
    }
  }
  else {
    //
    // Totally general implementation (but least efficient)
    //
    // z0 *= beta
    if( beta_ == ST::zero() ) {
      for( int j = 0; j < subDim; ++j, z0_val += z0_s )
        (*z0_val) = ST::zero();
    }
    else if( beta_ != ST::one() ) {
      for( int j = 0; j < subDim; ++j, z0_val += z0_s )
        (*z0_val) *= beta_;
    }
    // z0 += sum( alpha[k]*v[k], k=0...l_num_vecs-1)
    z0_val = targ_sub_vecs[0].values().begin();
    for( int j = 0; j < subDim; ++j, z0_val += z0_s ) {
      for( int k = 0; k < l_num_vecs; ++k ) {
        const Scalar
          &alpha_k = alpha_[k],
          &v_k_val = *v_val[k];
        (*z0_val) += alpha_k * v_k_val;
        v_val[k] += v_s[k];
      }
    }
  }
}


} // namespace RTOpPack


#endif // RTOPPACK_TOP_LINEAR_COMBINATION_DEF_HPP

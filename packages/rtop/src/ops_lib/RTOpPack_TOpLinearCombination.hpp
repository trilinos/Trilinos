// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
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

#ifndef RTOPPACK_TOP_LINEAR_COMBINATION_HPP
#define RTOPPACK_TOP_LINEAR_COMBINATION_HPP

#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_Workspace.hpp"

namespace RTOpPack {

/** \brief Linear combination transformation operator: <tt>z0[i] = beta*z0[i]
 * + sum( alpha[k]*v[k][i], k=0...num_vecs-1 ), i=0...n-1</tt>.
 *
 * This transformation operator only accepts <tt>num_targ_vec==1</tt>
 * but accepts any <tt>num_vecs > 0</tt>.
 *
 * Warning! this class can only be used in SPMD mode and not
 * client/server or master/slave.  You know what needs to happen for
 * this to work!
 */
template<class Scalar>
class TOpLinearCombination : public RTOpT<Scalar> {
public:
  /** \brief . */
  TOpLinearCombination(
    const int       num_vecs  = 0
    ,const Scalar   alpha[]   = NULL
    ,const Scalar   &beta     = Teuchos::ScalarTraits<Scalar>::zero()
    )
    :RTOpT<Scalar>("TOpLinearCombination")
     ,beta_(beta)
    { if(num_vecs) this->alpha(num_vecs,alpha); }
  /** \brief . */
  void beta( const Scalar& beta ) { beta_ = beta; }
  /** \brief . */
  Scalar beta() const { return beta_; }
  /** \brief . */
  void alpha( 
    const int       num_vecs
    ,const Scalar   alpha[]    ///< Array length <tt>num_vecs</tt>
    )
    {
      TEST_FOR_EXCEPT( num_vecs<=0 || alpha==NULL );
      alpha_.resize(0);
      alpha_.insert(alpha_.begin(),alpha,alpha+num_vecs);
    }
  /** \brief . */
  int num_vecs() const { return alpha_.size(); }
  /** \brief . */
  const Scalar* alpha() const { return &alpha_[0]; }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
      using Teuchos::Workspace;
      Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
      // Validate input
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPT( static_cast<int>(alpha_.size()) != num_vecs );
      TEST_FOR_EXCEPT( sub_vecs == NULL );
      TEST_FOR_EXCEPT( num_targ_vecs != 1 );
      TEST_FOR_EXCEPT( targ_sub_vecs == NULL );
#endif
      // Get pointers to local data
      const RTOpPack::index_type    subDim  = targ_sub_vecs[0].subDim();
      Scalar                        *z0_val = targ_sub_vecs[0].values();
      const ptrdiff_t               z0_s    = targ_sub_vecs[0].stride();
      Workspace<const Scalar*> v_val(wss,num_vecs,false);
      Workspace<ptrdiff_t>     v_s(wss,num_vecs,false);
      for( int k = 0; k < num_vecs; ++k ) {
#ifdef TEUCHOS_DEBUG
        TEST_FOR_EXCEPT( sub_vecs[k].subDim() != subDim );
        TEST_FOR_EXCEPT( sub_vecs[k].globalOffset() != targ_sub_vecs[0].globalOffset() );
#endif					
        v_val[k] = sub_vecs[k].values();
        v_s[k]   = sub_vecs[k].stride();
      }
      //
      // Perform the operation and specialize the cases for num_vecs = 1 and 2
      // in order to get good performance.
      //
      if( num_vecs == 1 ) {
        //
        // z0 = alpha*v0 + beta*z0
        //
        const Scalar alpha = alpha_[0], beta = beta_;
        const Scalar       *v0_val = v_val[0];
        const ptrdiff_t    v0_s    = v_s[0]; 
        if( beta==ST::zero() ) {
          // z0 = alpha*v0
          if( z0_s==1 && v0_s==1 ) {
            for( int j = 0; j < subDim; ++j )
              (*z0_val++) = alpha * (*v0_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s )
              (*z0_val) = alpha * (*v0_val);
          }
        }
        else if( beta==ST::one() ) {
          //
          // z0 = alpha*v0 + z0
          //
          if( z0_s==1 && v0_s==1 ) {
            for( int j = 0; j < subDim; ++j )
              (*z0_val++) += alpha * (*v0_val++);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s )
              (*z0_val) += alpha * (*v0_val);
          }
        }
        else {
          // z0 = alpha*v0 + beta*z0
          if( z0_s==1 && v0_s==1 ) {
            for( int j = 0; j < subDim; ++j, ++z0_val )
              (*z0_val) = alpha * (*v0_val++) + beta*(*z0_val);
          }
          else {
            for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s )
              (*z0_val) = alpha * (*v0_val) + beta*(*z0_val);
          }
        }
      }
      else if( num_vecs == 2 ) {
        //
        // z0 = alpha0*v0 + alpha1*v1 + beta*z0
        //
        const Scalar alpha0 = alpha_[0], alpha1=alpha_[1], beta = beta_;
        const Scalar     *v0_val = v_val[0];
        const ptrdiff_t  v0_s    = v_s[0]; 
        const Scalar     *v1_val = v_val[1];
        const ptrdiff_t  v1_s    = v_s[1]; 
        if( beta==ST::zero() ) {
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
        else if( beta==ST::one() ) {
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
                  (*z0_val) = (*v0_val++) + (*v1_val++) + beta*(*z0_val);
              }
              else {
                for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
                  (*z0_val) = (*v0_val) + (*v1_val) + beta*(*z0_val);
              }
            }
            else {
              // z0 = v0 + alpha1*v1 + beta*z0
              if( z0_s==1 && v0_s==1 && v1_s==1 ) {
                for( int j = 0; j < subDim; ++j, ++z0_val )
                  (*z0_val) = (*v0_val++) + alpha1*(*v1_val++) + beta*(*z0_val);
              }
              else {
                for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
                  (*z0_val) = (*v0_val) + alpha1*(*v1_val) + beta*(*z0_val);
              }
            }
          }
          else {
            if( alpha1 == ST::one() ) {
              // z0 = alpha0*v0 + v1 + beta*z0
              if( z0_s==1 && v0_s==1 && v1_s==1 ) {
                for( int j = 0; j < subDim; ++j, ++z0_val )
                  (*z0_val) = alpha0*(*v0_val++) + (*v1_val++) + beta*(*z0_val);
              }
              else {
                for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
                  (*z0_val) = alpha0*(*v0_val) + (*v1_val) + beta*(*z0_val);
              }
            }
            else {
              // z0 = alpha0*v0 + alpha1*v1 + beta*z0
              if( z0_s==1 && v0_s==1 && v1_s==1 ) {
                for( int j = 0; j < subDim; ++j, ++z0_val )
                  (*z0_val) = alpha0*(*v0_val++) + alpha1*(*v1_val++) + beta*(*z0_val);
              }
              else {
                for( int j = 0; j < subDim; ++j, z0_val+=z0_s, v0_val+=v0_s, v1_val+=v1_s )
                  (*z0_val) = alpha0*(*v0_val) + alpha1*(*v1_val) + beta*(*z0_val);
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
        // z0 += sum( alpha[k]*v[k], k=0...num_vecs-1)
        z0_val = targ_sub_vecs[0].values();
        for( int j = 0; j < subDim; ++j, z0_val += z0_s ) {
          for( int k = 0; k < num_vecs; ++k ) {
            const Scalar
              &alpha_k = alpha_[k],
              &v_k_val = *v_val[k];
            (*z0_val) += alpha_k * v_k_val;
            v_val[k] += v_s[k];
          }
        }
      }
    }
  //@}
private:
  Scalar                  beta_;
  std::vector<Scalar>  alpha_;
}; // class TOpLinearCombination

} // namespace RTOpPack

#endif // RTOPPACK_TOP_LINEAR_COMBINATION_HPP

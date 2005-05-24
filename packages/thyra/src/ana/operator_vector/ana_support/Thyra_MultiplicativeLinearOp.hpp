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

#ifndef THYRA_MULTIPLICATIVE_LINEAR_OP_HPP
#define THYRA_MULTIPLICATIVE_LINEAR_OP_HPP

#include "Thyra_MultiplicativeLinearOpDecl.hpp"
#include "Thyra_AssertOp.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
MultiplicativeLinearOp<Scalar>::MultiplicativeLinearOp()
{}

template<class Scalar>
MultiplicativeLinearOp<Scalar>::MultiplicativeLinearOp(
  const int                                                   numOps
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    Ops[]
  ,const Scalar                                               &gamma
  )
{
  initialize(numOps,Ops,gamma);
}

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::initialize(
  const int                                                   numOps
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    Ops[]
  ,const Scalar                                               &gamma
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( numOps <= 0 || Ops == NULL );
  for( int k = 0; k < numOps; ++k ) {
    TEST_FOR_EXCEPT( Ops[k].get() == NULL );
    if( k < numOps-1 ) {
      THYRA_ASSERT_VEC_SPACES(
        "MultiplicativeLinearOp<Scalar>::initialize(...)"
        ,*Ops[k]->domain(), *Ops[k+1]->range()
        );
    }
  }
#endif
  Ops_.resize(numOps);
  std::copy( Ops, Ops+numOps, Ops_.begin() );
  gamma_ = gamma;
}

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::uninitialize(
  const int                                             numOps
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >    Ops[]
  ,Scalar                                               *gamma
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( Ops!=NULL && numOps!=this->numOps() );
#endif
  if(Ops) std::copy( Ops_.begin(), Ops_.end(), Ops );
  if(gamma) *gamma = gamma_;
  Ops_.resize(0);
}


// Overridden from Teuchos::Describable
                                                
template<class Scalar>
std::string MultiplicativeLinearOp<Scalar>::describe() const
{
  assertInitialized();
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  std::ostringstream oss;
  oss << "MultiplicativeLinearOp<" << ST::name() << ">{numOps = "<<numOps()<<"}";
  return oss.str();
}

template<class Scalar>
std::ostream& MultiplicativeLinearOp<Scalar>::describe(
  std::ostream                         &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ,const std::string                   li
  ,const std::string                   is
  ) const
{
  assertInitialized();
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  const int numOps = Ops_.size();
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << this->describe() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      out << li << is << "type = \'MultiplicativeLinearOp<" << ST::name() << ">\', "
          << "rangeDim = " << this->range()->dim() << ", domainDim = " << this->domain()->dim() << std::endl
          << li << is << is << "numOps="<< numOps << std::endl
          << li << is << is << "Constituent LinearOpBase objects for M = Op[0]*...*Op[numOps-1]:\n";
      const std::string new_li = li+is+is+is;
      for( int k = 0; k < numOps; ++k ) {
        out << new_li << "Op["<<k<<"] =\n" << Teuchos::describe(*Ops_[k],verbLevel,new_li,is);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true);
  }
  return out;
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
MultiplicativeLinearOp<Scalar>::range() const
{
  assertInitialized();
  return Ops_[0]->range();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
MultiplicativeLinearOp<Scalar>::domain() const
{
  assertInitialized();
  return Ops_[numOps()-1]->domain();
}

template<class Scalar>
bool MultiplicativeLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
  bool opSupported = true;
  for( int k = 0; k < static_cast<int>(Ops_.size()); ++k )
    if(!Ops_[k]->opSupported(M_trans)) opSupported = false;
  return opSupported;
  // ToDo: Cache these?
}

// Overridden from LinearOpBase

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::apply(
  const ETransp                M_trans
  ,const VectorBase<Scalar>    &x
  ,VectorBase<Scalar>          *y
  ,const Scalar                alpha
  ,const Scalar                beta
  ) const
{
  this->apply(
    M_trans
    ,static_cast<const MultiVectorBase<Scalar>&>(x)
    ,static_cast<MultiVectorBase<Scalar>*>(y)
    ,alpha
    ,beta
    );
}

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  const int numOps = Ops_.size();
  const Index m = X.domain()->dim();
  if( real_trans(M_trans)==NOTRANS ) {
    //
    // Y = alpha * M * X + beta*Y
    // =>
    // Y = (alpha*gamma) * op(Op[0]) * op(Op[1]) * ... * op(Op[numOps-1]) * X + beta*Y
    //
    RefCountPtr<MultiVectorBase<Scalar> > T_kp1, T_k; // Temporary propagated between loops 
    for( int k = numOps-1; k >= 0; --k ) {
      RefCountPtr<MultiVectorBase<Scalar> >         Y_k;
      RefCountPtr<const MultiVectorBase<Scalar> >   X_k;
      if(k==0)        Y_k = rcp(Y,false);  else Y_k = T_k = createMembers(Ops_[k]->range(),m);
      if(k==numOps-1) X_k = rcp(&X,false); else X_k = T_kp1;
      if( k > 0 )
        Ops_[k]->apply(M_trans,*X_k,&*Y_k);
      else
        Ops_[k]->apply(M_trans,*X_k,&*Y_k,(alpha*gamma_),beta);
      T_kp1 = T_k;
    }
  }
  else {
    //
    // Y = alpha * M' * X + beta*Y
    // =>
    // Y = (alpha*gamma) * ( Op[numOps-1]' * ( .... ( Op[1]' * ( Op[0]' * X ) ) ... ) ) + beta * Y
    //
    RefCountPtr<MultiVectorBase<Scalar> > T_km1, T_k; // Temporary propagated between loops 
    for( int k = 0; k <= numOps-1; ++k ) {
      RefCountPtr<MultiVectorBase<Scalar> >         Y_k;
      RefCountPtr<const MultiVectorBase<Scalar> >   X_k;
      if(k==numOps-1)   Y_k = rcp(Y,false);  else Y_k = T_k = createMembers(Ops_[k]->domain(),m);
      if(k==0)          X_k = rcp(&X,false); else X_k = T_km1;
      if( k < numOps-1 )
        Ops_[k]->apply(M_trans,*X_k,&*Y_k);
      else
        Ops_[k]->apply(M_trans,*X_k,&*Y_k,(alpha*gamma_),beta);
      T_km1 = T_k;
    }
  }
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
MultiplicativeLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

}	// end namespace Thyra

#endif	// THYRA_MULTIPLICATIVE_LINEAR_OP_HPP

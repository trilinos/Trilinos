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

#ifndef THYRA_SCALED_ADJOINT_LINEAR_OP_HPP
#define THYRA_SCALED_ADJOINT_LINEAR_OP_HPP

#include "Thyra_ScaledAdjointLinearOpDecl.hpp"

namespace Thyra {

//Constructors/initializers/accessors

template<class Scalar>
void ScaledAdjointLinearOp<Scalar>::initialize(
  const Scalar                                               &scalar
  ,const ETransp                                             &transp
  ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> >   &Op
  )
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(
    Op.get()==NULL, std::invalid_argument
    ,"ScaledAdjointLinearOp<"<<ST::name()<<">::initialize(scalar,transp,Op): Error!, Op.get()==NULL is not allowed!"
    );
#endif // _DEBUG
  Teuchos::RefCountPtr<const ScaledAdjointLinearOp<Scalar> >
    saOp = Teuchos::rcp_dynamic_cast<const ScaledAdjointLinearOp<Scalar> >(Op);
  if(saOp.get()) {
    origOp_            = saOp->origOp_;
    overallScalar_     = saOp->overallScalar_*scalar;
    overallTransp_     = trans_trans(saOp->overallTransp_,transp) ;
    my_index_          = saOp->my_index_ + 1;
    allScalarETransp_  = saOp->allScalarETransp_;
  }
  else {
    origOp_            = Op;
    overallScalar_     = scalar;
    overallTransp_     = transp;
    my_index_          = 0;
    allScalarETransp_  = Teuchos::rcp(new allScalarETransp_t());
  }
  allScalarETransp_->push_back(ScalarETransp<Scalar>(scalar,transp));
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
ScaledAdjointLinearOp<Scalar>::getOp() const
{
  assertInitialized();
  if( my_index_ > 0 ) {
    const ScalarETransp<Scalar> &scalar_transp = allScalarETransp_->at(my_index_);
    Teuchos::RefCountPtr<ScaledAdjointLinearOp<Scalar> >
      Op = Teuchos::rcp(new ScaledAdjointLinearOp<Scalar>());
    Op->origOp_            = origOp_;
    Op->overallScalar_     = overallScalar_/scalar_transp.scalar;
    Op->overallTransp_     = trans_trans(overallTransp_,scalar_transp.transp);
    Op->my_index_          = my_index_ - 1;
    Op->allScalarETransp_  = allScalarETransp_;
    return Op;
  }
  return origOp_;
}

template<class Scalar>
void ScaledAdjointLinearOp<Scalar>::uninitialize(
  Scalar                                              *scalar
  ,ETransp                                            *transp
  ,Teuchos::RefCountPtr<const LinearOpBase<Scalar> >  *Op
  )
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;

  assertInitialized();

  const ScalarETransp<Scalar> &scalar_transp = allScalarETransp_->at(my_index_);
  if(scalar) *scalar = scalar_transp.scalar;
  if(transp) *transp = scalar_transp.transp;
  if(Op)     *Op     = getOp();

  origOp_           = Teuchos::null;
  overallScalar_    = ST::zero();
  overallTransp_    = NOTRANS;
  allScalarETransp_ = Teuchos::null;

}

// Overridden from Teuchos::Describable
                                                
template<class Scalar>
std::string ScaledAdjointLinearOp<Scalar>::description() const
{
  assertInitialized();
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  std::ostringstream oss;
  oss << "ScaledAdjointLinearOp<" << ST::name() << ">{overallScalar="
      << overallScalar() << ",overallTransp="<<toString(overallTransp())<<",origOp="
      << origOp_->description() << "}";
  return oss.str();
}

template<class Scalar>
std::ostream& ScaledAdjointLinearOp<Scalar>::describe(
  std::ostream                         &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ,const std::string                   li
  ,const std::string                   is
  ) const
{
  assertInitialized();
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      out << li << is << "type = \'ScaledAdjointLinearOp<" << ST::name() << ">\', "
          << "rangeDim = " << this->range()->dim() << ", domainDim = " << this->domain()->dim() << std::endl
          << li << is << is << "overallScalar="<< overallScalar() << std::endl
          << li << is << is << "overallTransp="<<toString(overallTransp()) << std::endl
          << li << is << is << "Constituent transformations:\n";
      for( int i = 0; i <= my_index_; ++i ) {
        const ScalarETransp<Scalar> &scalar_transp = (*allScalarETransp_)[my_index_-i];
        out << li << is << is;
        for( int j = 0; j <= i; ++j ) out << is;
        if(scalar_transp.scalar != ST::one() && scalar_transp.transp != NOTRANS)
          out << "scalar="<<scalar_transp.scalar<<",transp="<<toString(scalar_transp.transp)<<std::endl;
        else if(scalar_transp.scalar != ST::one())
          out << "scalar="<<scalar_transp.scalar<<std::endl;
        else if( scalar_transp.transp != NOTRANS )
          out << "transp="<<toString(scalar_transp.transp)<<std::endl;
        else
          out << "no-transformation\n";
      }
      std::ostringstream new_indent;
      new_indent << li << is << is << is;
      for( int i = 0; i <= my_index_; ++i ) new_indent << is;
      out << new_indent.str() << "origOp =\n" << Teuchos::describe(*origOp_,verbLevel,new_indent.str(),is);
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
ScaledAdjointLinearOp<Scalar>::range() const
{
  assertInitialized();
  return ( real_trans(this->overallTransp())==NOTRANS ? this->getOrigOp()->range() : this->getOrigOp()->domain() );
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
ScaledAdjointLinearOp<Scalar>::domain() const
{
  assertInitialized();
  return ( real_trans(this->overallTransp())==NOTRANS ? this->getOrigOp()->domain() : this->getOrigOp()->range() );
}

template<class Scalar>
bool ScaledAdjointLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
  assertInitialized();
  return Thyra::opSupported(*this->getOrigOp(),trans_trans(this->overallTransp(),M_trans));
}

// Overridden from LinearOpBase

template<class Scalar>
void ScaledAdjointLinearOp<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  assertInitialized();
  Thyra::apply(*this->getOrigOp(),trans_trans(M_trans,this->overallTransp()),X,Y,Scalar(this->overallScalar()*alpha),beta);
}

} // namespace Thyra

// Helper functions for creating scaled/adjoined linear operators.

template<class Scalar>
void Thyra::unwrap(
  const LinearOpBase<Scalar>      &Op
  ,Scalar                         *scalar
  ,ETransp                        *transp
  ,const LinearOpBase<Scalar>*    *origOp
  )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT( scalar==NULL );
  TEST_FOR_EXCEPT( transp==NULL );
  TEST_FOR_EXCEPT( origOp==NULL );
#endif
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  const ScaledAdjointLinearOp<Scalar>
    *saOp = dynamic_cast<const ScaledAdjointLinearOp<Scalar>*>(&Op);
  if(saOp) {
    *scalar = saOp->overallScalar();
    *transp = saOp->overallTransp();
    *origOp = &*saOp->getOrigOp();
  }
  else {
    *scalar = ST::one();
    *transp = NOTRANS;
    *origOp = &Op;
  }
}

#endif	// THYRA_SCALED_ADJOINT_LINEAR_OP_HPP

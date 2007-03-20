// @HEADER
// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
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

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_getConst.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

namespace Thyra {

// Constructors / initializers / accessors

EpetraLinearOp::EpetraLinearOp()
{}

EpetraLinearOp::EpetraLinearOp(
  const Teuchos::RefCountPtr<Epetra_Operator>                        &op
  ,ETransp                                                           opTrans
  ,EApplyEpetraOpAs                                                  applyAs
  ,EAdjointEpetraOp                                                  adjointSupport
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdRange
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdDomain
  )
{
  initialize(op,opTrans,applyAs,adjointSupport,spmdRange,spmdDomain);
}

void EpetraLinearOp::initialize(
  const Teuchos::RefCountPtr<Epetra_Operator>                        &op
  ,ETransp                                                           opTrans
  ,EApplyEpetraOpAs                                                  applyAs
  ,EAdjointEpetraOp                                                  adjointSupport
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdRange
  ,const Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >   &spmdDomain
  )
{

  // Validate input, allocate spaces, validate ...
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION( op.get()==NULL, std::invalid_argument, "EpetraLinearOp::initialize(...): Error!" );
#endif
  Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> > range, domain;
  if(spmdRange.get())
    range = spmdRange;
  else
    range = ( applyAs==EPETRA_OP_APPLY_APPLY ? allocateRange(op,opTrans)  : allocateDomain(op,opTrans) );
  if(spmdDomain.get())
    domain = spmdDomain;
  else
    domain  = ( applyAs==EPETRA_OP_APPLY_APPLY ? allocateDomain(op,opTrans) : allocateRange(op,opTrans)  );
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<EpetraLinearOp::Scalar> >
    sp_range = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<EpetraLinearOp::Scalar> >(range),
    sp_domain = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<EpetraLinearOp::Scalar> >(domain);

  // Set data (no exceptions should be thrown now)
  op_      = op;
  opTrans_ = opTrans;
  applyAs_ = applyAs;
  adjointSupport_ = adjointSupport;
  range_ = range;
  domain_ = domain;
  sp_range_ = sp_range;
  sp_domain_ = sp_domain;

}

void EpetraLinearOp::uninitialize(
  Teuchos::RefCountPtr<Epetra_Operator>                       *op
  ,ETransp                                                    *opTrans
  ,EApplyEpetraOpAs                                           *applyAs
  ,EAdjointEpetraOp                                           *adjointSupport
  ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >  *spmdRange
  ,Teuchos::RefCountPtr< const SpmdVectorSpaceBase<Scalar> >  *spmdDomain
  )
{

  if(op)      *op      = op_;
  if(opTrans) *opTrans = opTrans_;
  if(applyAs) *applyAs = applyAs_;
  if(adjointSupport) *adjointSupport = adjointSupport_;
  if(spmdRange) *spmdRange = range_;
  if(spmdDomain) *spmdDomain = domain_;

  op_      = Teuchos::null;
  opTrans_ = NOTRANS;
  applyAs_ = EPETRA_OP_APPLY_APPLY;
  adjointSupport_ = EPETRA_OP_ADJOINT_SUPPORTED;
  range_   = Teuchos::null;
  domain_  = Teuchos::null;
  sp_range_   = Teuchos::null;
  sp_domain_  = Teuchos::null;

}

Teuchos::RefCountPtr< const SpmdVectorSpaceBase<EpetraLinearOp::Scalar> >
EpetraLinearOp::spmdRange() const
{
  return range_;
}

Teuchos::RefCountPtr< const SpmdVectorSpaceBase<EpetraLinearOp::Scalar> >
EpetraLinearOp::spmdDomain() const
{
  return domain_;
}

Teuchos::RefCountPtr<Epetra_Operator>
EpetraLinearOp::epetra_op() 
{
  return op_;
}

Teuchos::RefCountPtr<const Epetra_Operator>
EpetraLinearOp::epetra_op() const 
{
  return op_;
}

// Overridden from EpetraLinearOpBase

void EpetraLinearOp::getEpetraOpView(
  Teuchos::RefCountPtr<Epetra_Operator>   *epetraOp
  ,ETransp                                *epetraOpTransp
  ,EApplyEpetraOpAs                       *epetraOpApplyAs
  ,EAdjointEpetraOp                       *epetraOpAdjointSupport
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(epetraOp==NULL);
  TEST_FOR_EXCEPT(epetraOpTransp==NULL);
  TEST_FOR_EXCEPT(epetraOpApplyAs==NULL);
  TEST_FOR_EXCEPT(epetraOpAdjointSupport==NULL);
#endif
  *epetraOp = op_;
  *epetraOpTransp = opTrans_;
  *epetraOpApplyAs = applyAs_;
  *epetraOpAdjointSupport = adjointSupport_;
}

void EpetraLinearOp::getEpetraOpView(
  Teuchos::RefCountPtr<const Epetra_Operator>   *epetraOp
  ,ETransp                                      *epetraOpTransp
  ,EApplyEpetraOpAs                             *epetraOpApplyAs
  ,EAdjointEpetraOp                             *epetraOpAdjointSupport
  ) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(epetraOp==NULL);
  TEST_FOR_EXCEPT(epetraOpTransp==NULL);
  TEST_FOR_EXCEPT(epetraOpApplyAs==NULL);
  TEST_FOR_EXCEPT(epetraOpAdjointSupport==NULL);
#endif
  *epetraOp = op_;
  *epetraOpTransp = opTrans_;
  *epetraOpApplyAs = applyAs_;
  *epetraOpAdjointSupport = adjointSupport_;
}

// Overridden from SingleScalarLinearOpBase

bool EpetraLinearOp::opSupported(ETransp M_trans) const
{
  return ( M_trans == NOTRANS ? true : adjointSupport_==EPETRA_OP_ADJOINT_SUPPORTED );
}

// Overridden from EuclideanLinearOpBase

Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<EpetraLinearOp::Scalar> >
EpetraLinearOp::rangeScalarProdVecSpc() const
{
  return sp_range_;
}

Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<EpetraLinearOp::Scalar> >
EpetraLinearOp::domainScalarProdVecSpc() const
{
  return sp_domain_;
}

void EpetraLinearOp::euclideanApply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X_in
  ,MultiVectorBase<Scalar>          *Y_inout
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  const ETransp real_M_trans = real_trans(M_trans);
#ifdef TEUCHOS_DEBUG
  // ToDo: Assert vector spaces!
  TEST_FOR_EXCEPTION(
    real_M_trans==TRANS && adjointSupport_==EPETRA_OP_ADJOINT_UNSUPPORTED
    ,Exceptions::OpNotSupported
    ,"EpetraLinearOp::apply(...): *this was informed that adjoints are not supported when initialized." 
    );
#endif
  //
  // Get Epetra_MultiVector objects for the arguments
  //
  Teuchos::RefCountPtr<const Epetra_MultiVector>
    X = get_Epetra_MultiVector(
      real_M_trans==NOTRANS ? getDomainMap() : getRangeMap()
      ,Teuchos::rcp(&X_in,false)
      );
  Teuchos::RefCountPtr<Epetra_MultiVector>
    Y;
  if( beta == 0 ) {
    Y = get_Epetra_MultiVector(
      real_M_trans==NOTRANS ? getRangeMap() : getDomainMap()
      ,Teuchos::rcp(Y_inout,false)
      );
  }
  //
  // Set the operator mode
  //
  /* We need to save the transpose state here, and then reset it after 
   * application. The reason for this is that if we later apply the 
   * operator outside Thyra (in Aztec, for instance), it will remember
   * the transpose flag set here. */
  bool oldState = op_->UseTranspose();
  op_->SetUseTranspose( real_trans(trans_trans(opTrans_,M_trans)) == NOTRANS ? false : true );
  //
  // Perform the operation
  //
  if( beta == 0.0 ) {
    // Y = M * X
    if( applyAs_ == EPETRA_OP_APPLY_APPLY )
      op_->Apply( *X, *Y );
    else if( applyAs_ == EPETRA_OP_APPLY_APPLY_INVERSE )
      op_->ApplyInverse( *X, *Y );
    else
      TEST_FOR_EXCEPT(true);
    // Y = alpha * Y
    if( alpha != 1.0 ) Y->Scale(alpha);
  }
  else {
    // Y_inout = beta * Y_inout
    if(beta != 0.0) scale( beta, Y_inout );
    else assign( Y_inout, 0.0 );
    // T = M * X
    Epetra_MultiVector T(
      real_M_trans == NOTRANS ? op_->OperatorRangeMap() : op_->OperatorDomainMap()
      ,X_in.domain()->dim()
      ,false
      );
    if( applyAs_ == EPETRA_OP_APPLY_APPLY )
      op_->Apply( *X, T );
    else if( applyAs_ == EPETRA_OP_APPLY_APPLY_INVERSE )
      op_->ApplyInverse( *X, T );
    else
      TEST_FOR_EXCEPT(true);
    // Y_inout += alpha * T
    update(
      alpha
      ,*create_MultiVector(
        Teuchos::rcp(&Teuchos::getConst(T),false)
        ,Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar> >(Y_inout->range(),true)
        ,Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<Scalar> >(Y_inout->domain(),true)
        )
      ,Y_inout
      );
  }
  // Reset the transpose state
  op_->SetUseTranspose(oldState);
}

// Overridden from LinearOpBase

Teuchos::RefCountPtr<const LinearOpBase<EpetraLinearOp::Scalar> >
EpetraLinearOp::clone() const
{
  assert(0); // ToDo: Implement when needed
  return Teuchos::null;
}

// Overridden from Teuchos::Describable

std::string EpetraLinearOp::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{";
  if(op_.get()) {
    oss << "op=\'"<<typeName(*op_)<<"\'";
    oss << ",dimRange="<<this->range()->dim();
    oss << ",dimDomain="<<this->domain()->dim();
  }
  else {
    oss << "op=NULL";
  }
  oss << "}";
  return oss.str();
}

void EpetraLinearOp::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::describe;
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim() << ",domainDim=" << this->domain()->dim() << "}\n";
      OSTab tab(out);
      if(op_.get()) {
        *out << "op="<<typeName(*op_)<<"\n";
        *out << "opTrans="<<toString(opTrans_)<<"\n";
        *out << "applyAs="<<toString(applyAs_)<<"\n";
        *out << "adjointSupport="<<toString(adjointSupport_)<<"\n";
      }
      else {
        *out << "op=NULL"<<"\n";
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

// protected

// Allocators for domain and range spaces

Teuchos::RefCountPtr< const SpmdVectorSpaceBase<EpetraLinearOp::Scalar> > 
EpetraLinearOp::allocateDomain(
  const Teuchos::RefCountPtr<Epetra_Operator>  &op 
  ,ETransp                                     op_trans 
  )  const
{
  return create_VectorSpace(Teuchos::rcp(&op->OperatorDomainMap(),false));
  // ToDo: What about the transpose argument???, test this!!!
}

Teuchos::RefCountPtr< const SpmdVectorSpaceBase<EpetraLinearOp::Scalar> > 
EpetraLinearOp::allocateRange(
  const Teuchos::RefCountPtr<Epetra_Operator>  &op 
  ,ETransp                                     op_trans 
  )  const
{
  return create_VectorSpace(Teuchos::rcp(&op->OperatorRangeMap(),false));
  // ToDo: What about the transpose argument???, test this!!!
}

// private

const Epetra_Map& EpetraLinearOp::getRangeMap() const
{
  return ( applyAs_ == EPETRA_OP_APPLY_APPLY ? op_->OperatorRangeMap() : op_->OperatorDomainMap() );
  // ToDo: What about the transpose argument???, test this!!!
}

const Epetra_Map& EpetraLinearOp::getDomainMap() const
{
  return ( applyAs_ == EPETRA_OP_APPLY_APPLY ? op_->OperatorDomainMap() : op_->OperatorRangeMap() );
  // ToDo: What about the transpose argument???, test this!!!
}

}	// end namespace Thyra

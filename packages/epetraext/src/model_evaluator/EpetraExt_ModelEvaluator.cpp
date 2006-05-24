// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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

#include "EpetraExt_ModelEvaluator.h"

namespace {

Teuchos::RefCountPtr<Epetra_Operator>
getLinearOp(
  const std::string                                                    &modelEvalDescription
  ,const EpetraExt::ModelEvaluator::Derivative                         &deriv
  ,const std::string                                                   &derivName
  )
{
  TEST_FOR_EXCEPTION(
    deriv.getDerivativeMultiVector().getMultiVector().get() != NULL, std::logic_error
    ,"For model \'" << modelEvalDescription << "\' the derivative \'"
    << derivName << "\' is of type Epetra_MultiVector and not of type Epetra_Operator!"
    );
  return deriv.getLinearOp();
}

Teuchos::RefCountPtr<Epetra_MultiVector>
getMultiVector(
  const std::string                                                    &modelEvalDescription
  ,const EpetraExt::ModelEvaluator::Derivative                         &deriv
  ,const std::string                                                   &derivName
  ,const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation  mvOrientation
  )
{
  TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get() != NULL, std::logic_error
    ,"For model \'" << modelEvalDescription << "\' the derivative \'"
    << derivName << "\' is of type Epetra_Operator and not of type Epetra_MultiVector!"
    );
  Teuchos::RefCountPtr<Epetra_MultiVector>
    mv = deriv.getDerivativeMultiVector().getMultiVector();
  if(mv.get()) {
    TEST_FOR_EXCEPTION(
      deriv.getDerivativeMultiVector().getOrientation()!=mvOrientation, std::logic_error
      ,"For model \'" << modelEvalDescription << "\' the derivative \'"
      << derivName << "\' if not the orientation \'" << toString(mvOrientation)
      << "\'"
      );
  }
  return mv;
}

} // namespace

namespace EpetraExt {

//
// ModelEvaluator::InArgs
//

ModelEvaluator::InArgs::InArgs()
{
  std::fill_n(&supports_[0],NUM_E_IN_ARGS_MEMBERS,false);
  t_     = 0.0;
  alpha_ = 0.0;
  beta_  = 0.0;
}

bool ModelEvaluator::InArgs::supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  return supports_[arg];
}

void ModelEvaluator::InArgs::_setSupports( EInArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\':Error, arg="<<arg<<" is invalid!"
    );
  supports_[arg] = supports;
}

void ModelEvaluator::InArgs::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(arg): model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

void ModelEvaluator::InArgs::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_l(l): model = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter l = " << l << " is not in the range [0,"<<Np()-1<<"]!"
    );
}

//
// ModelEvaluator::OutArgs
//

ModelEvaluator::OutArgs::OutArgs()
{
  std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false);
}

bool ModelEvaluator::OutArgs::supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  return supports_[arg];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  return supports_DfDp_[l];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  return supports_DgDx_[j];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDp arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_[ j*Np() + l ];
}

void ModelEvaluator::OutArgs::_setModelEvalDescription( const std::string &modelEvalDescription )
{
  modelEvalDescription_ = modelEvalDescription;
}

void ModelEvaluator::OutArgs::_set_Np_Ng(int Np, int Ng)
{
  if(Np) {
    supports_DfDp_.resize(Np);
    DfDp_.resize(Np);                 std::fill_n(DfDp_.begin(),Np,Derivative());
    DfDp_properties_.resize(Np);      std::fill_n(DfDp_properties_.begin(),Np,DerivativeProperties());
  }
  if(Ng) {
    g_.resize(Ng);
    supports_DgDx_.resize(Ng);
    DgDx_.resize(Ng);                 std::fill_n(DgDx_.begin(),Ng,Derivative());
    DgDx_properties_.resize(Ng);      std::fill_n(DgDx_properties_.begin(),Ng,DerivativeProperties());
  }
  if(Np && Ng) {
    const int NpNg = Np*Ng;
    supports_DgDp_.resize(NpNg);
    DgDp_.resize(NpNg);                 std::fill_n(DgDp_.begin(),NpNg,Derivative());
    DgDp_properties_.resize(NpNg);      std::fill_n(DgDp_properties_.begin(),NpNg,DerivativeProperties());
  }
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  supports_[arg] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& supports )
{
  assert_l(l);
  supports_DfDp_[l] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_[j] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_[ j*Np() + l ] = supports;
}

void ModelEvaluator::OutArgs::_set_W_properties( const DerivativeProperties &W_properties )
{
  W_properties_ = W_properties;
}

void ModelEvaluator::OutArgs::_set_DfDp_properties( int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_properties_[l] = properties;
}

void ModelEvaluator::OutArgs::_set_DgDx_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_properties_[j] = properties;
}

void ModelEvaluator::OutArgs::_set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_properties_[ j*Np() + l ] = properties;
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DfDp_[l].none(), std::logic_error
    ,"Thyra::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp(l) with index l = " << l << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_[j].none(), std::logic_error
    ,"Thyra::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDp arg, int j, int l) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_[ j*Np() + l ].none(), std::logic_error
    ,"Thyra::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    Np()==0, std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_l(l): model = \'"<<modelEvalDescription_<<"\':  Error, "
    "no auxiliary parameters subvectors p(l) are supported!!"
    );
  TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluator::OutArgs::assert_l(l): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter subvector p(l) index l = " << l << " is not in the range [0,"<<Np()-1<<"]!"
    );
}

void ModelEvaluator::OutArgs::assert_j(int j) const
{
  TEST_FOR_EXCEPTION(
    Ng()==0, std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_j(j): model = \'"<<modelEvalDescription_<<"\':  Error, "
    "no auxiliary functions g(j) are supported!!"
    );
  TEST_FOR_EXCEPTION(
    !( 0 <= j && j < Ng() ), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_j(j): model = \'"<<modelEvalDescription_<<"\':  Error, "
    "The auxiliary function g(j) index j = " << j << " is not in the range [0,"<<Ng()-1<<"]!"
    );
}

//
// ModelEvaluator
//

// Destructor

ModelEvaluator::~ModelEvaluator()
{}

// Vector maps

Teuchos::RefCountPtr<const Epetra_Map>
ModelEvaluator::get_p_map(int l) const
{ return Teuchos::null; }

Teuchos::RefCountPtr<const Epetra_Map>
ModelEvaluator::get_g_map(int j) const
{ return Teuchos::null; }

// Initial guesses for variables/parameters

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_x_init() const
{ return Teuchos::null; }

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_x_dot_init() const
{ return Teuchos::null; }

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_p_init(int l) const
{ return Teuchos::null; }

double ModelEvaluator::get_t_init() const
{ return 0.0; }

// Bounds for variables/parameters

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_x_lower_bounds() const
{ return Teuchos::null; }

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_x_upper_bounds() const
{ return Teuchos::null; }

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_p_lower_bounds(int l) const
{ return Teuchos::null; }

Teuchos::RefCountPtr<const Epetra_Vector>
ModelEvaluator::get_p_upper_bounds(int l) const
{ return Teuchos::null; }

double ModelEvaluator::get_t_lower_bound() const
{ return 0.0; }

double ModelEvaluator::get_t_upper_bound() const
{ return 0.0; }

// Factory functions for creating derivative objects

Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_W() const
{ return Teuchos::null; }

Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_DfDp_op(int l) const
{ return Teuchos::null; }

Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_DgDx_op(int j) const
{ return Teuchos::null; }

Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_DgDp_op( int j, int l ) const
{ return Teuchos::null; }

} // namespace EpetraExt

//
// Helper functions
//

std::string EpetraExt::toString( ModelEvaluator::EDerivativeMultiVectorOrientation orientation )
{
  switch(orientation) {
    case ModelEvaluator::DERIV_MV_BY_COL:
      return "DERIV_MV_BY_COL";
    case ModelEvaluator::DERIV_TRANS_MV_BY_ROW:
      return "DERIV_TRANS_MV_BY_ROW";
    default:
      TEST_FOR_EXCEPT(true);
  }
  return ""; // Should never be called
}

Teuchos::RefCountPtr<Epetra_Operator>
EpetraExt::get_DfDp_op(
  const int                                                            l
  ,const ModelEvaluator::OutArgs                                       &outArgs
  )
{
  std::ostringstream derivName; derivName << "DfDp("<<l<<")";
  return getLinearOp(
    outArgs.modelEvalDescription()
    ,outArgs.get_DfDp(l)
    ,derivName.str()
    );
}

Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraExt::get_DfDp_mv(
  const int                                                            l
  ,const ModelEvaluator::OutArgs                                       &outArgs
  )
{
  std::ostringstream derivName; derivName << "DfDp("<<l<<")";
  return getMultiVector(
    outArgs.modelEvalDescription()
    ,outArgs.get_DfDp(l)
    ,derivName.str()
    ,ModelEvaluator::DERIV_MV_BY_COL
    );
}

Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraExt::get_DgDx_mv(
  const int                                                            j
  ,const ModelEvaluator::OutArgs                                       &outArgs
  ,const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation  mvOrientation
  )
{
  std::ostringstream derivName; derivName << "DgDx("<<j<<")";
  return getMultiVector(
    outArgs.modelEvalDescription()
    ,outArgs.get_DgDx(j)
    ,derivName.str()
    ,mvOrientation
    );
}

Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraExt::get_DgDp_mv(
  const int                                                            j
  ,const int                                                           l
  ,const ModelEvaluator::OutArgs                                       &outArgs
  ,const EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation  mvOrientation
  )
{
  std::ostringstream derivName; derivName << "DgDp("<<j<<","<<l<<")";
  return getMultiVector(
    outArgs.modelEvalDescription()
    ,outArgs.get_DgDp(j,l)
    ,derivName.str()
    ,mvOrientation
    );
}

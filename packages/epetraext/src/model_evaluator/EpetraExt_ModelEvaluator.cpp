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
#include "Teuchos_implicit_cast.hpp"


namespace EpetraExt {


//
// ModelEvaluator::InArgs
//


ModelEvaluator::InArgs::InArgs()
  :modelEvalDescription_("WARNING!  THIS INARGS OBJECT IS UNINITALIZED!")
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
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}

bool ModelEvaluator::InArgs::supports(EInArgs_p_sg arg, int l) const
{
  assert_l(l);
  return supports_p_sg_[l];
}

bool ModelEvaluator::InArgs::supports(EInArgs_p_mp arg, int l) const
{
  assert_l(l);
  return supports_p_mp_[l];
}

void ModelEvaluator::InArgs::_setSupports( EInArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\':Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supports;
}

void ModelEvaluator::InArgs::_setSupports(EInArgs_p_sg arg, int l, bool supports)
{
  assert_l(l);
  supports_p_sg_[l] = supports;
}

void ModelEvaluator::InArgs::_setSupports(EInArgs_p_mp arg, int l, bool supports)
{
  assert_l(l);
  supports_p_mp_[l] = supports;
}


void ModelEvaluator::InArgs::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(arg): model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

void ModelEvaluator::InArgs::assert_supports(EInArgs_p_sg arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    !supports_p_sg_[l], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(IN_ARG_p_sg,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument p_sg(l) with index l = " << l << " is not supported!"
    );
}

void ModelEvaluator::InArgs::assert_supports(EInArgs_p_mp arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    !supports_p_mp_[l], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(IN_ARG_p_mp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument p_mp(l) with index l = " << l << " is not supported!"
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
  :modelEvalDescription_("WARNING!  THIS OUTARGS OBJECT IS UNINITALIZED!"),
  isFailed_( false )
{
  std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false);
}


bool ModelEvaluator::OutArgs::supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
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
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_dot arg, int j) const
{
  assert_j(j);
  return supports_DgDx_dot_[j];
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

bool ModelEvaluator::OutArgs::supports(EOutArgs_g_sg arg, int j) const
{
  assert_j(j);
  return supports_g_sg_[j];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDfDp_sg arg, int l) const
{
  assert_l(l);
  return supports_DfDp_sg_[l];
}


const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_dot_sg arg, int j) const
{
  assert_j(j);
  return supports_DgDx_dot_sg_[j];
}


const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_sg arg, int j) const
{
  assert_j(j);
  return supports_DgDx_sg_[j];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDp_sg arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_sg_[ j*Np() + l ];
}

bool ModelEvaluator::OutArgs::supports(EOutArgs_g_mp arg, int j) const
{
  assert_j(j);
  return supports_g_mp_[j];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDfDp_mp arg, int l) const
{
  assert_l(l);
  return supports_DfDp_mp_[l];
}


const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_dot_mp arg, int j) const
{
  assert_j(j);
  return supports_DgDx_dot_mp_[j];
}


const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_mp arg, int j) const
{
  assert_j(j);
  return supports_DgDx_mp_[j];
}

const ModelEvaluator::DerivativeSupport&
ModelEvaluator::OutArgs::supports(EOutArgsDgDp_mp arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_mp_[ j*Np() + l ];
}


bool ModelEvaluator::OutArgs::funcOrDerivesAreSet(EOutArgsMembers arg) const
{
  using Teuchos::implicit_cast;
  bool areSet = false;
  switch(arg) {
    case OUT_ARG_f: {
      if (!is_null(f_)) areSet = true;
      if (!is_null(W_)) areSet = true;
      for ( int l = 0; l < implicit_cast<int>(DfDp_.size()); ++l )
        if(!DfDp_[l].isEmpty()) areSet = true;
      break;
    }
    default:
      TEST_FOR_EXCEPTION(true,std::logic_error,
        "ModelEvaluator::OutArgs::funcOrDerivesAreSet(arg): Error, we can not handle"
        " the argument " << toString(arg) << "yet!");
  }
  return areSet;
}

void ModelEvaluator::OutArgs::setFailed() const
{
  isFailed_ = true;
  // TODO: Set objects to NaN?
}

bool ModelEvaluator::OutArgs::isFailed() const
{
  return isFailed_;
}


void ModelEvaluator::OutArgs::_setModelEvalDescription( const std::string &modelEvalDescription )
{
  modelEvalDescription_ = modelEvalDescription;
}


void ModelEvaluator::OutArgs::_set_Np_Ng(int Np, int Ng)
{
  if(Np) {
    supports_DfDp_.resize(Np);
    DfDp_.resize(Np);
    std::fill_n(DfDp_.begin(),Np,Derivative());
    DfDp_properties_.resize(Np);
    std::fill_n(DfDp_properties_.begin(),Np,DerivativeProperties());

    supports_DfDp_sg_.resize(Np);
    DfDp_sg_.resize(Np);
    std::fill_n(DfDp_sg_.begin(),Np,SGDerivative());
    DfDp_sg_properties_.resize(Np);
    std::fill_n(DfDp_sg_properties_.begin(),Np,DerivativeProperties());

    supports_DfDp_mp_.resize(Np);
    DfDp_mp_.resize(Np);
    std::fill_n(DfDp_mp_.begin(),Np,MPDerivative());
    DfDp_mp_properties_.resize(Np);
    std::fill_n(DfDp_mp_properties_.begin(),Np,DerivativeProperties());
  }
  if(Ng) {
    g_.resize(Ng);
    supports_DgDx_dot_.resize(Ng);
    DgDx_dot_.resize(Ng);
    std::fill_n(DgDx_dot_.begin(),Ng,Derivative());
    DgDx_dot_properties_.resize(Ng);
    std::fill_n(DgDx_dot_properties_.begin(),Ng,DerivativeProperties());
    supports_DgDx_.resize(Ng);
    DgDx_.resize(Ng);
    std::fill_n(DgDx_.begin(),Ng,Derivative());
    DgDx_properties_.resize(Ng);
    std::fill_n(DgDx_properties_.begin(),Ng,DerivativeProperties());

    g_sg_.resize(Ng);
    supports_g_sg_.resize(Ng);
    supports_DgDx_dot_sg_.resize(Ng);
    DgDx_dot_sg_.resize(Ng);
    std::fill_n(DgDx_dot_sg_.begin(),Ng,SGDerivative());
    DgDx_dot_sg_properties_.resize(Ng);
    std::fill_n(DgDx_dot_sg_properties_.begin(),Ng,DerivativeProperties());
    supports_DgDx_sg_.resize(Ng);
    DgDx_sg_.resize(Ng);
    std::fill_n(DgDx_sg_.begin(),Ng,SGDerivative());
    DgDx_sg_properties_.resize(Ng);
    std::fill_n(DgDx_sg_properties_.begin(),Ng,DerivativeProperties());

    g_mp_.resize(Ng);
    supports_g_mp_.resize(Ng);
    supports_DgDx_dot_mp_.resize(Ng);
    DgDx_dot_mp_.resize(Ng);
    std::fill_n(DgDx_dot_mp_.begin(),Ng,MPDerivative());
    DgDx_dot_mp_properties_.resize(Ng);
    std::fill_n(DgDx_dot_mp_properties_.begin(),Ng,DerivativeProperties());
    supports_DgDx_mp_.resize(Ng);
    DgDx_mp_.resize(Ng);
    std::fill_n(DgDx_mp_.begin(),Ng,MPDerivative());
    DgDx_mp_properties_.resize(Ng);
    std::fill_n(DgDx_mp_properties_.begin(),Ng,DerivativeProperties());
  }
  if(Np && Ng) {
    const int NpNg = Np*Ng;
    supports_DgDp_.resize(NpNg);
    DgDp_.resize(NpNg);
    std::fill_n(DgDp_.begin(),NpNg,Derivative());
    DgDp_properties_.resize(NpNg);
    std::fill_n(DgDp_properties_.begin(),NpNg,DerivativeProperties());

    supports_DgDp_sg_.resize(NpNg);
    DgDp_sg_.resize(NpNg);
    std::fill_n(DgDp_sg_.begin(),NpNg,SGDerivative());
    DgDp_sg_properties_.resize(NpNg);
    std::fill_n(DgDp_sg_properties_.begin(),NpNg,DerivativeProperties());

    supports_DgDp_mp_.resize(NpNg);
    DgDp_mp_.resize(NpNg);
    std::fill_n(DgDp_mp_.begin(),NpNg,MPDerivative());
    DgDp_mp_properties_.resize(NpNg);
    std::fill_n(DgDp_mp_properties_.begin(),NpNg,DerivativeProperties());
  }
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& supports )
{
  assert_l(l);
  supports_DfDp_[l] = supports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dot arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_dot_[j] = supports;
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

void ModelEvaluator::OutArgs::_setSupports( EOutArgs_g_sg arg, int j, bool supports )
{
  assert_j(j);
  supports_g_sg_[j] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp_sg arg, int l, const DerivativeSupport& supports )
{
  assert_l(l);
  supports_DfDp_sg_[l] = supports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dot_sg arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_dot_sg_[j] = supports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_sg arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_sg_[j] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDp_sg arg, int j, int l, const DerivativeSupport& supports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_sg_[ j*Np() + l ] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgs_g_mp arg, int j, bool supports )
{
  assert_j(j);
  supports_g_mp_[j] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp_mp arg, int l, const DerivativeSupport& supports )
{
  assert_l(l);
  supports_DfDp_mp_[l] = supports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_dot_mp_[j] = supports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_mp arg, int j, const DerivativeSupport& supports )
{
  assert_j(j);
  supports_DgDx_mp_[j] = supports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& supports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_mp_[ j*Np() + l ] = supports;
}


void ModelEvaluator::OutArgs::_set_W_properties( const DerivativeProperties &W_properties )
{
  W_properties_ = W_properties;
}

void ModelEvaluator::OutArgs::_set_WPrec_properties( const DerivativeProperties &WPrec_properties )
{
  WPrec_properties_ = WPrec_properties;
}

void ModelEvaluator::OutArgs::_set_DfDp_properties( int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_properties_[l] = properties;
}


void ModelEvaluator::OutArgs::_set_DgDx_dot_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  DgDx_dot_properties_[j] = properties;
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

void ModelEvaluator::OutArgs::_set_DfDp_sg_properties( int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DfDp_sg,l);
  DfDp_sg_properties_[l] = properties;
}


void ModelEvaluator::OutArgs::_set_DgDx_dot_sg_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_dot_sg,j);
  DgDx_dot_sg_properties_[j] = properties;
}


void ModelEvaluator::OutArgs::_set_DgDx_sg_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_sg,j);
  DgDx_sg_properties_[j] = properties;
}

void ModelEvaluator::OutArgs::_set_DgDp_sg_properties( int j, int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDp_sg,j,l);
  DgDp_sg_properties_[ j*Np() + l ] = properties;
}


void ModelEvaluator::OutArgs::_set_DfDp_mp_properties( int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  DfDp_mp_properties_[l] = properties;
}


void ModelEvaluator::OutArgs::_set_DgDx_dot_mp_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  DgDx_dot_mp_properties_[j] = properties;
}


void ModelEvaluator::OutArgs::_set_DgDx_mp_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  DgDx_mp_properties_[j] = properties;
}

void ModelEvaluator::OutArgs::_set_DgDp_mp_properties( int j, int l, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  DgDp_mp_properties_[ j*Np() + l ] = properties;
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DfDp_[l].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dot arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_dot_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dot,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dot(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_[j].none(), std::logic_error
    ,"TEpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDp arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DgDp_[ j*Np() + l ].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgs_g_sg arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    !supports_g_sg_[j], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_g_sg,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument g_sg(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp_sg arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DfDp_sg_[l].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp_sg,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp_sg(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dot_sg arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_dot_sg_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dot_sg,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dot_sg(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_sg arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_sg_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_sg,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_sg(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDp_sg arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DgDp_sg_[ j*Np() + l ].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp_sg,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp_sg(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgs_g_mp arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    !supports_g_mp_[j], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_g_mp,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument g_mp(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp_mp arg, int l) const
{
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DfDp_mp_[l].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp_mp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp_mp(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dot_mp arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_dot_mp_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dot_mp,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dot_mp(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_mp arg, int j) const
{
  assert_j(j);
  TEST_FOR_EXCEPTION(
    supports_DgDx_mp_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_mp,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_mp(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDp_mp arg, int j, int l) const
{
  assert_j(j);
  assert_l(l);
  TEST_FOR_EXCEPTION(
    supports_DgDp_mp_[ j*Np() + l ].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp_mp,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp_mp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
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
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_l(l): "
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

Teuchos::RefCountPtr<const Teuchos::Array<std::string> >
ModelEvaluator::get_p_names(int l) const
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


double ModelEvaluator::getInfBound() const
{
  return 1e+50;
}


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

Teuchos::RefCountPtr<EpetraExt::ModelEvaluator::Preconditioner>
ModelEvaluator::create_WPrec() const
{ return Teuchos::null; }

Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_DfDp_op(int l) const
{ return Teuchos::null; }

Teuchos::RefCountPtr<Epetra_Operator>
ModelEvaluator::create_DgDx_dot_op(int j) const
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


std::string EpetraExt::toString(
   ModelEvaluator::EDerivativeMultiVectorOrientation orientation
   )
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


std::string EpetraExt::toString( ModelEvaluator::EInArgsMembers inArg )
{
  switch(inArg) {
    case ModelEvaluator::IN_ARG_x_dot:
      return "IN_ARG_x_dot";
    case ModelEvaluator::IN_ARG_x:
      return "IN_ARG_x";
    case ModelEvaluator::IN_ARG_x_dot_poly:
      return "IN_ARG_x_dot_poly";
    case ModelEvaluator::IN_ARG_x_poly:
      return "IN_ARG_x_poly";
    case ModelEvaluator::IN_ARG_x_dot_sg:
      return "IN_ARG_x_dot_sg";
    case ModelEvaluator::IN_ARG_x_sg:
      return "IN_ARG_x_sg";
    case ModelEvaluator::IN_ARG_x_dot_mp:
      return "IN_ARG_x_dot_mp";
    case ModelEvaluator::IN_ARG_x_mp:
      return "IN_ARG_x_mp";
    case ModelEvaluator::IN_ARG_t:
      return "IN_ARG_t";
    case ModelEvaluator::IN_ARG_alpha:
      return "IN_ARG_alpha";
    case ModelEvaluator::IN_ARG_beta:
      return "IN_ARG_beta";
    default:
      TEST_FOR_EXCEPT("Invalid inArg!");
  }
  return ""; // Will never be executed!
}


std::string EpetraExt::toString( ModelEvaluator::EOutArgsMembers outArg )
{
  switch(outArg) {
    case  ModelEvaluator::OUT_ARG_f:
      return "OUT_ARG_f";
    case ModelEvaluator::OUT_ARG_W:
      return "OUT_ARG_W";
    case ModelEvaluator::OUT_ARG_WPrec:
      return "OUT_ARG_WPrec";
    case ModelEvaluator::OUT_ARG_f_poly:
      return "OUT_ARG_f_poly";
    case ModelEvaluator::OUT_ARG_f_sg:
      return "OUT_ARG_f_sg";
    case ModelEvaluator::OUT_ARG_W_sg:
      return "OUT_ARG_W_sg";
    case ModelEvaluator::OUT_ARG_f_mp:
      return "OUT_ARG_f_mp";
    case ModelEvaluator::OUT_ARG_W_mp:
      return "OUT_ARG_W_mp";
    default:
      TEST_FOR_EXCEPT("Invalid outArg!");
  }
  return ""; // Will never be executed!
}


Teuchos::RefCountPtr<Epetra_Operator>
EpetraExt::getLinearOp(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName
  )
{
  TEST_FOR_EXCEPTION(
    deriv.getMultiVector().get() != NULL, std::logic_error
    ,"For model \'" << modelEvalDescription << "\' the derivative \'"
    << derivName << "\' is of type Epetra_MultiVector and not of type Epetra_Operator!"
    );
  return deriv.getLinearOp();
}


Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraExt::getMultiVector(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  )
{
  TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get() != NULL, std::logic_error
    ,"For model \'" << modelEvalDescription << "\' the derivative \'"
    << derivName << "\' is of type Epetra_Operator and not of type Epetra_MultiVector!"
    );
  Teuchos::RefCountPtr<Epetra_MultiVector>
    mv = deriv.getMultiVector();
  if(mv.get()) {
    TEST_FOR_EXCEPTION(
      deriv.getMultiVectorOrientation()!=mvOrientation, std::logic_error
      ,"For model \'" << modelEvalDescription << "\' the derivative \'"
      << derivName << "\' if not the orientation \'" << toString(mvOrientation)
      << "\'"
      );
  }
  return mv;
}


Teuchos::RefCountPtr<Epetra_Operator>
EpetraExt::get_DfDp_op(
  const int l,
  const ModelEvaluator::OutArgs &outArgs
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
  const int l,
  const ModelEvaluator::OutArgs &outArgs
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
EpetraExt::get_DgDx_dot_mv(
  const int j,
  const ModelEvaluator::OutArgs &outArgs,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  )
{
  std::ostringstream derivName; derivName << "DgDx_dot("<<j<<")";
  return getMultiVector(
    outArgs.modelEvalDescription(),
    outArgs.get_DgDx_dot(j),
    derivName.str(),
    mvOrientation
    );
}


Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraExt::get_DgDx_mv(
  const int j,
  const ModelEvaluator::OutArgs &outArgs,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  )
{
  std::ostringstream derivName; derivName << "DgDx("<<j<<")";
  return getMultiVector(
    outArgs.modelEvalDescription(),
    outArgs.get_DgDx(j),
    derivName.str(),
    mvOrientation
    );
}


Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraExt::get_DgDp_mv(
  const int j,
  const int l,
  const ModelEvaluator::OutArgs &outArgs,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  )
{
  std::ostringstream derivName; derivName << "DgDp("<<j<<","<<l<<")";
  return getMultiVector(
    outArgs.modelEvalDescription(),
    outArgs.get_DgDp(j,l),
    derivName.str(),
    mvOrientation
    );
}

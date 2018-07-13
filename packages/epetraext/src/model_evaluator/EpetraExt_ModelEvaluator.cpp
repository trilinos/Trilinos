//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//@HEADER

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
  omega_ = 0.0;
  beta_  = 0.0;
}


bool ModelEvaluator::InArgs::supports(EInArgsMembers arg) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
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

void ModelEvaluator::InArgs::_setSupports( EInArgsMembers arg, bool supportsIt )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\':Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supportsIt;
}

void ModelEvaluator::InArgs::_setSupports(EInArgs_p_sg arg, int l, bool supportsIt)
{
  assert_l(l);
  supports_p_sg_[l] = supportsIt;
}

void ModelEvaluator::InArgs::_setSupports(EInArgs_p_mp arg, int l, bool supportsIt)
{
  assert_l(l);
  supports_p_mp_[l] = supportsIt;
}


void ModelEvaluator::InArgs::assert_supports(EInArgsMembers arg) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(arg): model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

void ModelEvaluator::InArgs::assert_supports(EInArgs_p_sg arg, int l) const
{
  assert_l(l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_p_sg_[l], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(IN_ARG_p_sg,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument p_sg(l) with index l = " << l << " is not supported!"
    );
}

void ModelEvaluator::InArgs::assert_supports(EInArgs_p_mp arg, int l) const
{
  assert_l(l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_p_mp_[l], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(IN_ARG_p_mp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument p_mp(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::InArgs::assert_l(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
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
  TEUCHOS_TEST_FOR_EXCEPTION(
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
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_dotdot arg, int j) const
{
  assert_j(j);
  return supports_DgDx_dotdot_[j];
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
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_dotdot_sg arg, int j) const
{
  assert_j(j);
  return supports_DgDx_dotdot_sg_[j];
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
ModelEvaluator::OutArgs::supports(EOutArgsDgDx_dotdot_mp arg, int j) const
{
  assert_j(j);
  return supports_DgDx_dotdot_mp_[j];
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
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
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


void ModelEvaluator::OutArgs::_setModelEvalDescription( const std::string &theModelEvalDescription )
{
  modelEvalDescription_ = theModelEvalDescription;
}


void ModelEvaluator::OutArgs::_set_Np_Ng(int Np_in, int Ng_in)
{
  if(Np_in) {
    supports_DfDp_.resize(Np_in);
    DfDp_.resize(Np_in);
    std::fill_n(DfDp_.begin(),Np_in,Derivative());
    DfDp_properties_.resize(Np_in);
    std::fill_n(DfDp_properties_.begin(),Np_in,DerivativeProperties());

    supports_DfDp_sg_.resize(Np_in);
    DfDp_sg_.resize(Np_in);
    std::fill_n(DfDp_sg_.begin(),Np_in,SGDerivative());
    DfDp_sg_properties_.resize(Np_in);
    std::fill_n(DfDp_sg_properties_.begin(),Np_in,DerivativeProperties());

    supports_DfDp_mp_.resize(Np_in);
    DfDp_mp_.resize(Np_in);
    std::fill_n(DfDp_mp_.begin(),Np_in,MPDerivative());
    DfDp_mp_properties_.resize(Np_in);
    std::fill_n(DfDp_mp_properties_.begin(),Np_in,DerivativeProperties());
  }
  if(Ng_in) {
    g_.resize(Ng_in);
    supports_DgDx_dot_.resize(Ng_in);
    DgDx_dot_.resize(Ng_in);
    std::fill_n(DgDx_dot_.begin(),Ng_in,Derivative());
    DgDx_dot_properties_.resize(Ng_in);
    std::fill_n(DgDx_dot_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_dotdot_.resize(Ng_in);
    DgDx_dotdot_.resize(Ng_in);
    std::fill_n(DgDx_dotdot_.begin(),Ng_in,Derivative());
    DgDx_dotdot_properties_.resize(Ng_in);
    std::fill_n(DgDx_dotdot_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_.resize(Ng_in);
    DgDx_.resize(Ng_in);
    std::fill_n(DgDx_.begin(),Ng_in,Derivative());
    DgDx_properties_.resize(Ng_in);
    std::fill_n(DgDx_properties_.begin(),Ng_in,DerivativeProperties());

    g_sg_.resize(Ng_in);
    supports_g_sg_.resize(Ng_in);
    supports_DgDx_dot_sg_.resize(Ng_in);
    DgDx_dot_sg_.resize(Ng_in);
    std::fill_n(DgDx_dot_sg_.begin(),Ng_in,SGDerivative());
    DgDx_dot_sg_properties_.resize(Ng_in);
    std::fill_n(DgDx_dot_sg_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_dotdot_sg_.resize(Ng_in);
    DgDx_dotdot_sg_.resize(Ng_in);
    std::fill_n(DgDx_dotdot_sg_.begin(),Ng_in,SGDerivative());
    DgDx_dotdot_sg_properties_.resize(Ng_in);
    std::fill_n(DgDx_dotdot_sg_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_sg_.resize(Ng_in);
    DgDx_sg_.resize(Ng_in);
    std::fill_n(DgDx_sg_.begin(),Ng_in,SGDerivative());
    DgDx_sg_properties_.resize(Ng_in);
    std::fill_n(DgDx_sg_properties_.begin(),Ng_in,DerivativeProperties());

    g_mp_.resize(Ng_in);
    supports_g_mp_.resize(Ng_in);
    supports_DgDx_dot_mp_.resize(Ng_in);
    DgDx_dot_mp_.resize(Ng_in);
    std::fill_n(DgDx_dot_mp_.begin(),Ng_in,MPDerivative());
    DgDx_dot_mp_properties_.resize(Ng_in);
    std::fill_n(DgDx_dot_mp_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_dotdot_mp_.resize(Ng_in);
    DgDx_dotdot_mp_.resize(Ng_in);
    std::fill_n(DgDx_dotdot_mp_.begin(),Ng_in,MPDerivative());
    DgDx_dotdot_mp_properties_.resize(Ng_in);
    std::fill_n(DgDx_dotdot_mp_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_mp_.resize(Ng_in);
    DgDx_mp_.resize(Ng_in);
    std::fill_n(DgDx_mp_.begin(),Ng_in,MPDerivative());
    DgDx_mp_properties_.resize(Ng_in);
    std::fill_n(DgDx_mp_properties_.begin(),Ng_in,DerivativeProperties());
  }
  if(Np_in && Ng_in) {
    const int NpNg = Np_in*Ng_in;
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

void ModelEvaluator::OutArgs::_setSupports( EOutArgsMembers arg, bool supportsIt )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_<<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supportsIt;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& theSupports )
{
  assert_l(l);
  supports_DfDp_[l] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dot arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_dot_[j] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dotdot arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_dotdot_[j] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_[j] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& theSupports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_[ j*Np() + l ] = theSupports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgs_g_sg arg, int j, bool supportsIt )
{
  assert_j(j);
  supports_g_sg_[j] = supportsIt;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp_sg arg, int l, const DerivativeSupport& theSupports )
{
  assert_l(l);
  supports_DfDp_sg_[l] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dot_sg arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_dot_sg_[j] = theSupports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dotdot_sg arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_dotdot_sg_[j] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_sg arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_sg_[j] = theSupports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDp_sg arg, int j, int l, const DerivativeSupport& theSupports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_sg_[ j*Np() + l ] = theSupports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgs_g_mp arg, int j, bool supportsIt )
{
  assert_j(j);
  supports_g_mp_[j] = supportsIt;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDfDp_mp arg, int l, const DerivativeSupport& theSupports )
{
  assert_l(l);
  supports_DfDp_mp_[l] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_dot_mp_[j] = theSupports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_dotdot_mp arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_dotdot_mp_[j] = theSupports;
}


void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDx_mp arg, int j, const DerivativeSupport& theSupports )
{
  assert_j(j);
  supports_DgDx_mp_[j] = theSupports;
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& theSupports )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_mp_[ j*Np() + l ] = theSupports;
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

void ModelEvaluator::OutArgs::_set_DgDx_dotdot_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_dotdot,j);
  DgDx_dotdot_properties_[j] = properties;
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

void ModelEvaluator::OutArgs::_set_DgDx_dotdot_sg_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_dotdot_sg,j);
  DgDx_dotdot_sg_properties_[j] = properties;
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

void ModelEvaluator::OutArgs::_set_DgDx_dotdot_mp_properties( int j, const DerivativeProperties &properties )
{
  assert_supports(OUT_ARG_DgDx_dotdot_mp,j);
  DgDx_dotdot_mp_properties_[j] = properties;
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp arg, int l) const
{
  assert_l(l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DfDp_[l].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dot arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDx_dot_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dot,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dot(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dotdot arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDx_dotdot_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dotdot,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dotdot(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDp_[ j*Np() + l ].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgs_g_sg arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_g_sg_[j], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_g_sg,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument g_sg(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp_sg arg, int l) const
{
  assert_l(l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DfDp_sg_[l].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp_sg,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp_sg(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dot_sg arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDx_dot_sg_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dot_sg,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dot_sg(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dotdot_sg arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDx_dotdot_sg_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dotdot_sg,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dotdot_sg(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_sg arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDp_sg_[ j*Np() + l ].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp_sg,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp_sg(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgs_g_mp arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_g_mp_[j], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_g_mp,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument g_mp(j) with index j = " << j << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsDfDp_mp arg, int l) const
{
  assert_l(l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DfDp_mp_[l].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DfDp_mp,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DfDp_mp(l) with index l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dot_mp arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDx_dot_mp_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dot_mp,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dot_mp(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_dotdot_mp arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDx_dotdot_mp_[j].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDx_dotdot_mp,j): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDx_dotdot_mp(j) with index j = " << j << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_supports(EOutArgsDgDx_mp arg, int j) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    supports_DgDp_mp_[ j*Np() + l ].none(), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(OUT_ARG_DgDp_mp,j,l): "
    "model = \'"<<modelEvalDescription_<<"\': Error,"
    "The argument DgDp_mp(j,l) with indexes j = " << j << " and l = " << l << " is not supported!"
    );
}


void ModelEvaluator::OutArgs::assert_l(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Np()==0, std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_l(l): model = \'"<<modelEvalDescription_<<"\':  Error, "
    "no auxiliary parameters subvectors p(l) are supported!!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_l(l): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter subvector p(l) index l = " << l << " is not in the range [0,"<<Np()-1<<"]!"
    );
}


void ModelEvaluator::OutArgs::assert_j(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    Ng()==0, std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_j(j): model = \'"<<modelEvalDescription_<<"\':  Error, "
    "no auxiliary functions g(j) are supported!!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
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


Teuchos::RCP<const Epetra_Map>
ModelEvaluator::get_p_map(int l) const
{ return Teuchos::null; }

Teuchos::RCP<const Teuchos::Array<std::string> >
ModelEvaluator::get_p_names(int l) const
{ return Teuchos::null; }

Teuchos::RCP<const Epetra_Map>
ModelEvaluator::get_g_map(int j) const
{ return Teuchos::null; }

Teuchos::ArrayView<const std::string>
ModelEvaluator::get_g_names(int j) const
{ return Teuchos::null; }


// Initial guesses for variables/parameters


Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_x_init() const
{ return Teuchos::null; }

Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_x_dot_init() const
{ return Teuchos::null; }

Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_x_dotdot_init() const
{ return Teuchos::null; }

Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_p_init(int l) const
{ return Teuchos::null; }

double ModelEvaluator::get_t_init() const
{ return 0.0; }


// Bounds for variables/parameters


double ModelEvaluator::getInfBound() const
{
  return 1e+50;
}


Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_x_lower_bounds() const
{ return Teuchos::null; }


Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_x_upper_bounds() const
{ return Teuchos::null; }


Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_p_lower_bounds(int l) const
{ return Teuchos::null; }


Teuchos::RCP<const Epetra_Vector>
ModelEvaluator::get_p_upper_bounds(int l) const
{ return Teuchos::null; }


double ModelEvaluator::get_t_lower_bound() const
{ return 0.0; }


double ModelEvaluator::get_t_upper_bound() const
{ return 0.0; }


// Factory functions for creating derivative objects


Teuchos::RCP<Epetra_Operator>
ModelEvaluator::create_W() const
{ return Teuchos::null; }

Teuchos::RCP<EpetraExt::ModelEvaluator::Preconditioner>
ModelEvaluator::create_WPrec() const
{ return Teuchos::null; }

Teuchos::RCP<Epetra_Operator>
ModelEvaluator::create_DfDp_op(int l) const
{ return Teuchos::null; }

Teuchos::RCP<Epetra_Operator>
ModelEvaluator::create_DgDx_dot_op(int j) const
{ return Teuchos::null; }

Teuchos::RCP<Epetra_Operator>
ModelEvaluator::create_DgDx_dotdot_op(int j) const
{ return Teuchos::null; }

Teuchos::RCP<Epetra_Operator>
ModelEvaluator::create_DgDx_op(int j) const
{ return Teuchos::null; }

Teuchos::RCP<Epetra_Operator>
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
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  TEUCHOS_UNREACHABLE_RETURN(""); // Should never be called
}


std::string EpetraExt::toString( ModelEvaluator::EInArgsMembers inArg )
{
  switch(inArg) {
    case ModelEvaluator::IN_ARG_x_dot:
      return "IN_ARG_x_dot";
    case ModelEvaluator::IN_ARG_x_dotdot:
      return "IN_ARG_x_dotdot";
    case ModelEvaluator::IN_ARG_x:
      return "IN_ARG_x";
    case ModelEvaluator::IN_ARG_x_dot_poly:
      return "IN_ARG_x_dot_poly";
    case ModelEvaluator::IN_ARG_x_dotdot_poly:
      return "IN_ARG_x_dotdot_poly";
    case ModelEvaluator::IN_ARG_x_poly:
      return "IN_ARG_x_poly";
    case ModelEvaluator::IN_ARG_x_dot_sg:
      return "IN_ARG_x_dot_sg";
    case ModelEvaluator::IN_ARG_x_dotdot_sg:
      return "IN_ARG_x_dotdot_sg";
    case ModelEvaluator::IN_ARG_x_sg:
      return "IN_ARG_x_sg";
    case ModelEvaluator::IN_ARG_x_dot_mp:
      return "IN_ARG_x_dot_mp";
    case ModelEvaluator::IN_ARG_x_dotdot_mp:
      return "IN_ARG_x_dotdot_mp";
    case ModelEvaluator::IN_ARG_x_mp:
      return "IN_ARG_x_mp";
    case ModelEvaluator::IN_ARG_t:
      return "IN_ARG_t";
    case ModelEvaluator::IN_ARG_alpha:
      return "IN_ARG_alpha";
    case ModelEvaluator::IN_ARG_omega:
      return "IN_ARG_omega";
    case ModelEvaluator::IN_ARG_beta:
      return "IN_ARG_beta";
    case ModelEvaluator::IN_ARG_step_size:
      return "IN_ARG_step_size";
    case ModelEvaluator::IN_ARG_stage_number:
      return "IN_ARG_stage_number";
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid inArg!");
  }
  TEUCHOS_UNREACHABLE_RETURN(""); // Will never be executed!
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
      TEUCHOS_TEST_FOR_EXCEPT("Invalid outArg!");
  }
  TEUCHOS_UNREACHABLE_RETURN(""); // Will never be executed!
}


Teuchos::RCP<Epetra_Operator>
EpetraExt::getLinearOp(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    deriv.getMultiVector().get() != NULL, std::logic_error
    ,"For model \'" << modelEvalDescription << "\' the derivative \'"
    << derivName << "\' is of type Epetra_MultiVector and not of type Epetra_Operator!"
    );
  return deriv.getLinearOp();
}


Teuchos::RCP<Epetra_MultiVector>
EpetraExt::getMultiVector(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    deriv.getLinearOp().get() != NULL, std::logic_error
    ,"For model \'" << modelEvalDescription << "\' the derivative \'"
    << derivName << "\' is of type Epetra_Operator and not of type Epetra_MultiVector!"
    );
  Teuchos::RCP<Epetra_MultiVector>
    mv = deriv.getMultiVector();
  if(mv.get()) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      deriv.getMultiVectorOrientation()!=mvOrientation, std::logic_error
      ,"For model \'" << modelEvalDescription << "\' the derivative \'"
      << derivName << "\' if not the orientation \'" << toString(mvOrientation)
      << "\'"
      );
  }
  return mv;
}


Teuchos::RCP<Epetra_Operator>
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


Teuchos::RCP<Epetra_MultiVector>
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


Teuchos::RCP<Epetra_MultiVector>
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


Teuchos::RCP<Epetra_MultiVector>
EpetraExt::get_DgDx_dotdot_mv(
  const int j,
  const ModelEvaluator::OutArgs &outArgs,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  )
{
  std::ostringstream derivName; derivName << "DgDx_dotdot("<<j<<")";
  return getMultiVector(
    outArgs.modelEvalDescription(),
    outArgs.get_DgDx_dotdot(j),
    derivName.str(),
    mvOrientation
    );
}


Teuchos::RCP<Epetra_MultiVector>
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


Teuchos::RCP<Epetra_MultiVector>
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

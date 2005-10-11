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
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  return supports_[arg];
}

void ModelEvaluator::InArgs::_setSupports( EInArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"*this = \'"<<modelEvalDescription_<<"\':Error, arg="<<arg<<" is invalid!"
    );
  supports_[arg] = supports;
}

void ModelEvaluator::InArgs::assert_supports(EInArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_supports(arg): *this = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

void ModelEvaluator::InArgs::assert_l(int l) const
{
  TEST_FOR_EXCEPTION(
    !( 1 <= l && l <= Np() ), std::logic_error
    ,"EpetraExt::ModelEvaluator::InArgs::assert_l(l): *this = \'"<<modelEvalDescription_<<"\': Error, "
    "The parameter l = " << l << " is not in the range [1,"<<Np()<<"]!"
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
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  return supports_[arg];
}

void ModelEvaluator::OutArgs::_setSupports( EOutArgsMembers arg, bool supports )
{
  TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"*this = \'"<<modelEvalDescription_<<"\': Error, arg="<<arg<<" is invalid!"
    );
  supports_[arg] = supports;
}

void ModelEvaluator::OutArgs::assert_supports(EOutArgsMembers arg) const
{
  TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_supports(arg): "
    "*this = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << arg << " is not supported!"
    );
}

void ModelEvaluator::OutArgs::assert_j(int j) const
{
  TEST_FOR_EXCEPTION(
    !( 1 <= j && j <= Ng() ), std::logic_error
    ,"EpetraExt::ModelEvaluator::OutArgs::assert_j(j): *this = \'"<<modelEvalDescription_<<"\':  Error, "
    "The auxiliary function g(j) index j = " << j << " is not in the range [1,"<<Ng()<<"]!"
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

} // namespace EpetraExt

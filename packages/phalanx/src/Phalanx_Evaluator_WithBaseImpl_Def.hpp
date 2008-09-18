// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_WITHBASEIMPL_DEF_H
#define PHX_EVALUATOR_WITHBASEIMPL_DEF_H

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::
EvaluatorWithBaseImpl(const std::string& evaluator_name) :
  name_(evaluator_name)
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::EvaluatorWithBaseImpl() :
  name_("???")
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::~EvaluatorWithBaseImpl()
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::FieldTag& ft)
{ 
  PHX::FTPredRef pred(ft);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(evaluated_.begin(), evaluated_.end(), pred);
  
  if ( test == evaluated_.end() )
    evaluated_.push_back(ft.clone());
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::Field<DataT>& f)
{ 
  this->addEvaluatedField(f.fieldTag());
}

//**********************************************************************
template<typename Traits>
template<typename DataT, PHX::ArrayOrder Order,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::MDField<DataT,Order,Tag0,Tag1,Tag2,Tag3,
		  Tag4,Tag5,Tag6,Tag7>& f)
{ 
  this->addEvaluatedField(f.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::FieldTag& ft)
{
  PHX::FTPredRef pred(ft);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(required_.begin(), required_.end(), pred);
  
  if ( test == required_.end() )
    required_.push_back(ft.clone());
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::Field<DataT>& v)
{
  this->addDependentField(v.fieldTag());
}

//**********************************************************************
template<typename Traits>
template<typename DataT, PHX::ArrayOrder Order,
	 typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	 typename Tag4, typename Tag5, typename Tag6, typename Tag7>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::MDField<DataT,Order,Tag0,Tag1,Tag2,Tag3,
		  Tag4,Tag5,Tag6,Tag7>& f)
{
  this->addDependentField(f.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
setName(const std::string& name)
{ name_ = name; }

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::evaluatedFields() const
{ return evaluated_; }

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::dependentFields() const
{ return required_; }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ }

//**********************************************************************
template<typename Traits>
const std::string& PHX::EvaluatorWithBaseImpl<Traits>::
getName() const
{return name_;}

//**********************************************************************

#endif

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

#ifndef PHX_SCALAR_CONTAINER_DEF_HPP
#define PHX_SCALAR_CONTAINER_DEF_HPP

#include "Teuchos_TestForException.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_TypeStrings.hpp"
#include <sstream>

// *************************************************************************
template <typename EvalT, typename Traits>
PHX::EvaluationContainer<EvalT, Traits>::EvaluationContainer() :
  post_registration_setup_called_(false)
{
  this->vp_manager_.setEvaluationTypeName( PHX::typeAsString<EvalT>() );
  this->data_container_template_manager_.buildObjects();
}

// *************************************************************************
template <typename EvalT, typename Traits> 
PHX::EvaluationContainer<EvalT, Traits>::~EvaluationContainer()
{

}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
requireField(const PHX::FieldTag& f)
{
  this->vp_manager_.requireField(f);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  this->vp_manager_.registerEvaluator(p);
}

// *************************************************************************
template <typename EvalT, typename Traits> 
void PHX::EvaluationContainer<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  // Figure out all evaluator dependencies
  if ( !(this->vp_manager_.sortingCalled()) )
    this->vp_manager_.sortAndOrderEvaluators();

  // Determine total amount of memory for all variables
  allocator_.reset();
  
  const std::vector< Teuchos::RCP<PHX::FieldTag> >& var_list = 
    this->vp_manager_.getFieldTags();

  std::vector< Teuchos::RCP<PHX::FieldTag> >::const_iterator  var = 
    var_list.begin();
  for (; var != var_list.end(); ++var) {

    typename DCTM::iterator it = data_container_template_manager_.begin();
    for (; it != data_container_template_manager_.end(); ++it) {
      
      if ((*var)->dataTypeInfo() == it->dataTypeInfo()) {

	std::size_t size_of_data_type = it->getSizeOfDataType();

	std::size_t num_elements = (*var)->dataLayout().size();

	allocator_.addRequiredChunk(size_of_data_type, num_elements);
      }
    }
  }

  allocator_.setup();

  // Allocate field data arrays
  //std::vector<PHX::FieldTag>::const_iterator  var = var_list.begin();
  for (var = var_list.begin(); var != var_list.end(); ++var) {

    typename DCTM::iterator it = data_container_template_manager_.begin();
    for (; it != data_container_template_manager_.end(); ++it) {
      
      if ((*var)->dataTypeInfo() == it->dataTypeInfo())
	it->allocateField(*var, allocator_);

    }
  }

  // Allow field evaluators to grab pointers to relevant field data
  this->vp_manager_.postRegistrationSetup(d,fm);

  post_registration_setup_called_ = true;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the evaluateFields() method for that type!");
#endif

  this->vp_manager_.evaluateFields(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the preEvaluate() method for that type!");
#endif

  this->vp_manager_.preEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION( !(this->setupCalled()) , std::logic_error,
		      "You must call post registration setup for each evaluation type before calling the postEvaluate() method for that type!");
#endif

  this->vp_manager_.postEvaluate(d);
}

// *************************************************************************
template <typename EvalT, typename Traits> template <typename DataT>
Teuchos::ArrayRCP<DataT> 
PHX::EvaluationContainer<EvalT, Traits>::getFieldData(const PHX::FieldTag& f)
{
  Teuchos::ArrayRCP<DataT> r = 
    data_container_template_manager_.template getAsObject<DataT>()->
    getFieldData(f);
  return r;
}


// *************************************************************************
template <typename EvalT, typename Traits>
bool PHX::EvaluationContainer<EvalT, Traits>::setupCalled() const
{
  return post_registration_setup_called_;
}

// *************************************************************************
template <typename EvalT, typename Traits>
void PHX::EvaluationContainer<EvalT, Traits>::print(std::ostream& os) const
{
  std::string type = PHX::typeAsString<EvalT>();

  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Starting PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << this->vp_manager_ << std::endl;
  typename DCTM::const_iterator it = data_container_template_manager_.begin();
  for (; it != data_container_template_manager_.end(); ++it)
    os << *it << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Finished PHX::EvaluationContainer Output" << std::endl;
  os << "Evaluation Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << std::endl;
}

// *************************************************************************

#endif 

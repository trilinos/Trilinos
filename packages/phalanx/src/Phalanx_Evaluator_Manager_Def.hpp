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

#ifndef PHX_FIELD_EVALUATOR_MANAGER_DEF_HPP
#define PHX_FIELD_EVALUATOR_MANAGER_DEF_HPP

#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "Teuchos_TestForException.hpp"
#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"

//=======================================================================
template<typename Traits>
PHX::EvaluatorManager<Traits>::
EvaluatorManager(const std::string& evaluation_type_name) :
  evaluation_type_name_(evaluation_type_name),
  sorting_called_(false)
{ }

//=======================================================================
template<typename Traits>
PHX::EvaluatorManager<Traits>::~EvaluatorManager()
{ }

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
requireField(const PHX::FieldTag& t)
{
  FTPredRef pred(t);
  std::vector< Teuchos::RCP<PHX::FieldTag> >::iterator i = 
    find_if(fields_.begin(), fields_.end(), pred);
  
  if (i == fields_.end())
    fields_.push_back(t.clone());
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  varProviders.push_back(p);
  providerVariables.push_back(p->evaluatedFields());
  providerRequirements.push_back(p->dependentFields());
  providerNames.push_back(p->getName());

#ifdef CHARONDEBUG
  char const* const methodID =
    "charon:EvaluatorManager::registerEvaluator";
  if (DO_DEBUG_OUTPUT(methodID,10)) {
    std::cout << "Registered provider: " << p->getName()
              << std::endl
              << "Element block: " << blockID << std::endl
              << "Index in vector: " << varProviders.size() - 1
              << std::endl << "Evaluates:" << std::endl;
    for (uint i=0; i < p->evaluatedFields().size(); ++i) {
      std::cout << "   " << p->evaluatedFields()[i].name() << std::endl;
    }
  }
#endif

  /*!
    \todo RPP: need to add a check to make sure multiple providers
    can't supply the same variable.
  */
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
sortAndOrderEvaluators()
{
  if (sorting_called_) {
    std::string msg = "Setup was already called.  ";
    msg += "Don't call setup more than once!";
    TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  // Construct the order in which providers need to be called
  this->createProviderEvaluationOrder();

  /*
   * After we have figured out which providers to call, we need to
   * ensure that all variables a provider will provide are added to
   * the varaible list.  For example, if someone writes a provider
   * that evaluates both DENSITY and VISCOSITY, but the physics only
   * requires DENSITY, an exception would be thrown when the provider
   * tries to get an index for the VISCOSITY variable from the
   * VariableArray.  We need to ensure that all provided variables
   * have storage allocated in the array - i.e. register the VISCOSITY
   * as a variable if it was not.
   */
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++) {
    std::size_t k = providerEvalOrderIndex[i];
    for (std::size_t j = 0; j <providerVariables[k].size(); j++)
      this->requireField(*(providerVariables[k][j]));
  }
  
  sorting_called_ = true;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  // Call each providers' post registration setup
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->postRegistrationSetup(vm);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::createProviderEvaluationOrder()
{
  // Before sorting provider order, we need to add any intermediate
  // variables to the fields_ that are not specified by the operators.
  bool done = false;
  while (!done) {
    bool addedVariables = false;
    
    for (std::size_t i = 0; i < fields_.size(); i++) {
      const PHX::FieldTag& v = *(fields_[i]);
      
      // Loop over providers and add any requirements as variables.
      for (std::size_t prov = 0; prov < providerVariables.size(); prov++) {
	for (std::size_t var = 0; var < providerVariables[prov].size(); var++) {
	  if (*(providerVariables[prov][var]) == v) {
	    // Loop over requirements to see if they are in the variable list.
	    for (std::size_t r = 0; r < providerRequirements[prov].size(); r++) {
	      bool isVariable = false;
	      for (std::size_t j = 0; j < fields_.size(); j++) {
		if (*(fields_[j]) == *(providerRequirements[prov][r]))
		  isVariable = true;
	      }
	      if (!isVariable) {
		fields_.push_back(providerRequirements[prov][r]);
		addedVariables = true;
	      }
	    }
	  }
	}
      }
    }
    if (!addedVariables)
      done = true;
  }
  
  std::vector<Teuchos::RCP<PHX::FieldTag> > tmpList = fields_;
  std::vector<Teuchos::RCP<PHX::FieldTag> > tmpProvided;
  
  // Loop over variable list until it is empty or we fail to remove var
  while (tmpList.size() > 0) {
    
    bool removedVariable = false;
    
    // Loop over all varibles still in the list until we find a
    // Provider that can remove a varible
    bool foundProvider = false;
    int providerIndex = -1;
    for (std::size_t var = 0; var < tmpList.size(); var++) {
      
      foundProvider = false;
      providerIndex = -1;
      
      // Loop over variable providers to find one that supplies this variable
      for (std::size_t prov = 0; prov < varProviders.size(); prov++) {
	
	// Loop over provided variable names in provider[prov]
	for (std::size_t i = 0; i < providerVariables[prov].size(); i++) {
	  
	  if (*(tmpList[var]) == *(providerVariables[prov][i])) {
	    foundProvider = true;
	    providerIndex = prov;
	    break;
	  }
	  
	}
	
	if (foundProvider)
	  break;
      }
      

      // Make sure requirements are satisfied for this provider
      bool requirementsSatisfied = true;
      
      if (foundProvider) {
	if (providerRequirements[providerIndex].size() > 0) {
	  
	  for (std::size_t req = 0;
	       req < providerRequirements[providerIndex].size();
	       req++) {
	    bool requiredVariableFound = false;
	    for (std::size_t j = 0; j < tmpProvided.size(); j++) {
	      if (*(providerRequirements[providerIndex][req]) == 
		  *(tmpProvided[j]))
		requiredVariableFound = true;
	    }
	    if (!requiredVariableFound) {
	      requirementsSatisfied = false;
	      break;
	    }
	    
	  }
	}
      }
      
      if (foundProvider && requirementsSatisfied) {
	
	// Remove the variable and exit loop
	std::vector<Teuchos::RCP<PHX::FieldTag> >::iterator p = 
	  tmpList.begin();
	tmpList.erase(p+var);
	// Add all vars to provided list and remove all variables
	// that this provider adds
	for (std::size_t i = 0; i < providerVariables[providerIndex].size(); i++) {
	  tmpProvided.push_back(providerVariables[providerIndex][i]);
	  for (std::size_t j = 0; j < tmpList.size(); j++) {
	    if (*(providerVariables[providerIndex][i]) == *(tmpList[j])) {
	      std::vector<Teuchos::RCP<PHX::FieldTag> >::iterator a = 
		tmpList.begin();
	      tmpList.erase(a+j);
	      break;
	    }
	  }
	}
	providerEvalOrderIndex.push_back(providerIndex);
	removedVariable = true;
	break;
      }

    }  // for (std::size_t var = 0; var < tmpList.size(); var++) {

    if (!removedVariable) {
      std::string msg = "EvaluatorManager";
      msg += evaluation_type_name_;
      msg += " \nCould not meet dependencies!\n";
      msg += "The following variables either have no provider or have a\n";
      msg += "provider but could not satisfy provider requirements:\n\n";
      std::ostringstream ost;
      for (std::size_t i = 0; i < tmpList.size(); i++)
	ost << *(tmpList[i]) << std::endl;
      msg += ost.str();
      msg += "\nPrinting EvaluatorManager:\n";
      std::ostringstream ost2;
      ost2 << *this << std::endl;
      msg += ost2.str();
      TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    
  } // While tmpList.size() != 0
  
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
evaluateFields(typename Traits::EvalData d)
{
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->evaluateFields(d);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->preEvaluate(d);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  for (std::size_t i = 0; i < providerEvalOrderIndex.size(); i++)
    (varProviders[providerEvalOrderIndex[i]])->postEvaluate(d);
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::
setEvaluationTypeName(const std::string& evaluation_type_name)
{
  evaluation_type_name_ = evaluation_type_name;
}

//=======================================================================
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >& 
PHX::EvaluatorManager<Traits>::getFieldTags()
{
  return fields_;
}

//=======================================================================
template<typename Traits>
bool PHX::EvaluatorManager<Traits>::sortingCalled() const
{
  return sorting_called_;
}

//=======================================================================
template<typename Traits>
void PHX::EvaluatorManager<Traits>::print(std::ostream& os) const
{
  os << "******************************************************" << std::endl;
  os << "PHX::EvaluatorManager" << std::endl;
  os << "Evaluation Type = " << evaluation_type_name_ << std::endl;
  os << "******************************************************" << std::endl;

  os << "\n** Starting Required Field List" << std::endl;
  for (std::size_t i = 0; i < fields_.size(); i++) {
    os << *(this->fields_[i]) << std::endl;
  }
  os << "** Finished Required Field List" << std::endl;

  os << "\n** Starting Registered Field Evaluators" << std::endl;
  for (std::size_t i = 0; i < varProviders.size(); i++) {
    os << "Evaluator[" << i << "]: " << providerNames[i] << std::endl;
    os << "  *Evaluates:" << std::endl;
    for (std::size_t j = 0; j < providerVariables[i].size(); j++)
      os << "    " << *((this->providerVariables[i])[j]) << std::endl;
    os << "  *Dependencies:";
    if (providerRequirements[i].size() == 0) {
      os << " None!" << std::endl;
    }
    else {
      os << std::endl;
      for (std::size_t j = 0; j < providerRequirements[i].size(); j++)
	os << "    " << *((this->providerRequirements[i])[j]) << std::endl;
    }
  }
  os << "** Finished Registered Field Evaluators" << std::endl;


  os << "\n** Starting Evaluator Order" << std::endl;
  for (std::size_t k = 0; k < providerEvalOrderIndex.size(); k++) {
    os << k << "    " << providerEvalOrderIndex[k] << std::endl;
  }
  os << "\nDetails:\n";
  for (std::size_t k = 0; k < providerEvalOrderIndex.size(); k++) {
    int i = providerEvalOrderIndex[k];
    os << "Evaluator[" << i << "]: " << providerNames[i] << std::endl;
    os << "  *Evaluates:" << std::endl;
    for (std::size_t j = 0; j < providerVariables[i].size(); j++)
      os << "    " << *((this->providerVariables[i])[j]) << std::endl;
    os << "  *Dependencies:";
    if (providerRequirements[i].size() == 0) {
      os << " None!" << std::endl;
    }
    else {
      os << std::endl;
      for (std::size_t j = 0; j < providerRequirements[i].size(); j++)
	os << "    " << *((this->providerRequirements[i])[j]) << std::endl;
    }
  }
  os << "** Finished Provider Evaluation Order" << std::endl;

  os << "******************************************************" << std::endl;
  os << "Finished PHX::EvaluatorManager" << std::endl;
  os << "Evaluation Type = " << evaluation_type_name_ << std::endl;
  os << "******************************************************" << std::endl;

}

//=======================================================================
template<typename Traits>
std::ostream&
PHX::operator<<(std::ostream& os, const PHX::EvaluatorManager<Traits>& m)
{
  m.print(os);
  return os;
}

//=======================================================================

#endif

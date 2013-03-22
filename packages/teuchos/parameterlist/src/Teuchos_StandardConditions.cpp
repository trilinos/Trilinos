// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// @HEADER



#include "Teuchos_StandardConditions.hpp"

namespace Teuchos{

ParameterCondition::ParameterCondition(RCP<const ParameterEntry> parameter):
  parameterEntry_(parameter)
{
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(parameter),
    InvalidConditionException,
    "Parameter conditions can't be given a null parameter" <<
    std::endl << std::endl);
}

Dependency::ConstParameterEntryList 
ParameterCondition::getAllParameters() const
{
  Dependency::ConstParameterEntryList toReturn;
  toReturn.insert(getParameter());
  return toReturn;
}

BoolLogicCondition::BoolLogicCondition(ConstConditionList& conditions):
  conditions_(conditions)
{
  TEUCHOS_TEST_FOR_EXCEPTION(conditions_.size() ==0,
    InvalidConditionException,
    "You must provide at least one condition "
    "when you're constructing a BoolLogicCondition. "
    << std::endl << std::endl <<
    "Error: Empty condition list given to a BoolLogicCondition "
    "constructor.");
}


void BoolLogicCondition::addCondition(RCP<const Condition> toAdd){
  conditions_.append(toAdd);
}

bool BoolLogicCondition::isConditionTrue() const{
  ConstConditionList::const_iterator it = conditions_.begin();
  bool toReturn = (*it)->isConditionTrue();
  ++it;
  for(;it != conditions_.end(); ++it){
    toReturn = applyOperator(toReturn,(*it)->isConditionTrue());
  }
  return toReturn;
}

bool BoolLogicCondition::containsAtLeasteOneParameter() const{
  for(
    ConstConditionList::const_iterator it=conditions_.begin();
    it!=conditions_.end();
    ++it)
  {
    if((*it)->containsAtLeasteOneParameter()){
      return true;
    }
  }
  return false;
}

Dependency::ConstParameterEntryList 
BoolLogicCondition::getAllParameters() const{
  Dependency::ConstParameterEntryList toReturn;
  Dependency::ConstParameterEntryList currentList;
  for(
    ConstConditionList::const_iterator it = conditions_.begin();
    it != conditions_.end();
    ++it)
  {
    currentList = (*it)->getAllParameters();
    toReturn.insert(currentList.begin(), currentList.end());
  }
  return toReturn;
}

OrCondition::OrCondition(ConstConditionList& conditions):
  BoolLogicCondition(conditions){}

bool OrCondition::applyOperator(bool op1, bool op2) const{
  return op1 || op2;
}

RCP<OrCondition> DummyObjectGetter<OrCondition>::getDummyObject(){
  Condition::ConstConditionList dummyList;
  dummyList.append(DummyObjectGetter<BoolCondition>::getDummyObject());
  return rcp(new OrCondition(dummyList));
}

AndCondition::AndCondition(ConstConditionList& conditions):
  BoolLogicCondition(conditions){}

bool AndCondition::applyOperator(bool op1, bool op2) const{
  return op1 && op2;
}

RCP<AndCondition> DummyObjectGetter<AndCondition>::getDummyObject(){
  Condition::ConstConditionList dummyList;
  dummyList.append(DummyObjectGetter<BoolCondition>::getDummyObject());
  return rcp(new AndCondition(dummyList));
}

EqualsCondition::EqualsCondition(ConstConditionList& conditions):
  BoolLogicCondition(conditions){}

bool EqualsCondition::applyOperator(bool op1, bool op2) const{
  return op1 == op2;
}

RCP<EqualsCondition> DummyObjectGetter<EqualsCondition>::getDummyObject(){
  Condition::ConstConditionList dummyList;
  dummyList.append(DummyObjectGetter<BoolCondition>::getDummyObject());
  return rcp(new EqualsCondition(dummyList));
}

NotCondition::NotCondition(RCP<const Condition> childCondition):
  childCondition_(childCondition)
{
  TEUCHOS_TEST_FOR_EXCEPTION(childCondition_.is_null(),
    InvalidConditionException,
    "OOOOOOOOPppppps! Looks like you tried "
    "to give me "
    "a null pointer when you were making a not conditon. "
    "That's a no no. Go back and "
    "checkout your not conditions and make sure you didn't "
    "give any of them a null pointer "
    "as an argument to the constructor." << std::endl << std::endl <<
    "Error: Null pointer given to NotCondition constructor.");
}

bool NotCondition::isConditionTrue() const{
  return (!childCondition_->isConditionTrue());
}

bool NotCondition::containsAtLeasteOneParameter() const{
  return childCondition_->containsAtLeasteOneParameter();
}

Dependency::ConstParameterEntryList NotCondition::getAllParameters() const{
  return childCondition_->getAllParameters();
}

RCP<NotCondition> DummyObjectGetter<NotCondition>::getDummyObject(){
  return rcp(new NotCondition(
    DummyObjectGetter<BoolCondition>::getDummyObject()));
}

StringCondition::StringCondition(
  RCP<const ParameterEntry> parameter,
  std::string value):
  ParameterCondition(parameter), 
  values_(ValueList(1,value))
{
  checkParameterType();
}

StringCondition::StringCondition(
  RCP<const ParameterEntry> parameter,
  ValueList values):
  ParameterCondition(parameter), 
  values_(values)
{
  checkParameterType();
}

void StringCondition::checkParameterType(){
  TEUCHOS_TEST_FOR_EXCEPTION(!getParameter()->isType<std::string>(),
    InvalidConditionException,
    "The parameter of a String Condition "
    "must be of type string." << std::endl << 
    "Expected type: " << TypeNameTraits<std::string>::name() << std::endl <<
    "Actual type: " << getParameter()->getAny().typeName() << 
    std::endl << std::endl);
}
  

bool StringCondition::evaluateParameter() const{
  return  find(
    values_.begin(), values_.end(), 
    getValue<std::string>(*getParameter())) != values_.end();
}

RCP<StringCondition> DummyObjectGetter<StringCondition>::getDummyObject(){
  std::string empty = "";
  return rcp(new StringCondition(rcp(new ParameterEntry(empty)), empty));
}

BoolCondition::BoolCondition(RCP<const ParameterEntry> parameter):
  ParameterCondition(parameter)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!getParameter()->isType<bool>(),
    InvalidConditionException,
    "The parameter of a Bool Condition "
    "must be of type " << TypeNameTraits<bool>::name() << std::endl <<
    "Expected type: Bool" << std::endl <<
    "Actual type: " << getParameter()->getAny().typeName() <<
    std::endl << std::endl);
}

bool BoolCondition::evaluateParameter() const{
  return getValue<bool>(*getParameter());
}

RCP<BoolCondition> DummyObjectGetter<BoolCondition>::getDummyObject(){
  return rcp(new BoolCondition(rcp(new ParameterEntry(true))));
}


} //namespace Teuchos


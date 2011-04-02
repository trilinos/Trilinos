// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef Teuchos_STANDARD_FUNCTION_OBJECTS_H
#define Teuchos_STANDARD_FUNCTION_OBJECTS_H

/*! \file Teuchos_FunctionObject.hpp
    \brief An object representation of a function
*/

#include "Teuchos_FunctionObject.hpp"


namespace Teuchos{


template<class ArgType, class ReturnType>
class SingleArguementFunctionObject : public FunctionObject<ReturnType>{

public:
  SingleArguementFunctionObject(){}

  SingleArguementFunctionObject(ArgType arguementValue):
    arguementValue_(arguementValue){}

  inline ArgType getArguementValue() const{
    return arguementValue_;
  }

  inline void setArguementValue(ArgType arguementValue){
    arguementValue_ = arguementValue;
  }

private:

  ArgType arguementValue_;

};

template<class OperandType>
class SingleOperatorFunction : 
  public SingleArguementFunctionObject<OperandType, OperandType>
{
public:
  SingleOperatorFunction(OperandType modifyingOperand):
    SingleArguementFunctionObject<OperandType, OperandType>(),
    _modifyingOperand(modifyingOperand){}

  SingleOperatorFunction(
    OperandType modifyingOperand,
    OperandType arguementValue):
    SingleArguementFunctionObject<OperandType, OperandType>(arguementValue),
    _modifyingOperand(modifyingOperand){}

  inline OperandType getModifiyingOperand() const{
    return _modifyingOperand;
  }

  inline OperandType setModifyingOperand(OperandType newOperand){
    _modifyingOperand = newOperand;
  }

private:
  OperandType _modifyingOperand;
};

template<class OperandType>
class SubtractionFunction :
  public SingleOperatorFunction<OperandType>
{
public:
  SubtractionFunction(OperandType amountToSubtract):
    SingleOperatorFunction<OperandType>(amountToSubtract){}

  SubtractionFunction(
    OperandType amountToSubtract,
    OperandType arguementValue):
    SingleOperatorFunction<OperandType>(amountToSubtract, arguementValue){}

  inline OperandType runFunction() const{
    return 
      SingleArguementFunctionObject<OperandType, OperandType>::getArguementValue()
      -
      SingleOperatorFunction<OperandType>::getModifiyingOperand();
  }
};

template<class OperandType>
class AdditionFunction :
  public SingleOperatorFunction<OperandType>
{
public:
  AdditionFunction(OperandType amountToAdd):
    SingleOperatorFunction<OperandType>(amountToAdd){}

  AdditionFunction(
    OperandType amountToAdd,
    OperandType arguementValue):
    SingleOperatorFunction<OperandType>(amountToAdd, arguementValue){}

  inline OperandType runFunction() const{
    return 
      SingleArguementFunctionObject<OperandType, OperandType>::getArguementValue()
      +
      SingleOperatorFunction<OperandType>::getModifiyingOperand();
  }
};

} // namespace Teuchos


#endif

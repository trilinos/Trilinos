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

/**
 * \brief A simple function object that
 * applies a given operand to a spcified arguement using a specific operator.
 */
template<class OperandType>
class SimpleFunctionObject : public FunctionObject
{
public:

  //! @name Constructors/Destructors
  //@{

  /**
   * \brief Constructs a SimpleFunctionObject.
   *
   */
  SimpleFunctionObject():
    FunctionObject(){}
  
  /**
   * \brief Constructs a SimpleFunctionObject.
   *
   * @param modifyingOperand The operand that will be 
   * modifying the arguement given in the runFuction function.
   */
  SimpleFunctionObject(OperandType modifyingOperand):
    FunctionObject(),
    _modifyingOperand(modifyingOperand){}

  //@}

  /**
   * Runs the desired function on the arguement and returns the
   * result.
   *
   * @param arguement Arguement on which the function should be run.
   * @return The result of running the function on the give arguement.
   */
  virtual OperandType runFunction(OperandType arguement) const=0;

  //! @name Getters/Setters
  //@{

  /**
   * \brief Returns the modifying operand.
   */
  inline OperandType getModifiyingOperand() const{
    return _modifyingOperand;
  }

  /**
   * \brief Sets the modifyingOperand.
   *
   * @param newOperand The new modifyingOperand
   * to use.
   */
  inline OperandType setModifyingOperand(OperandType newOperand){
    _modifyingOperand = newOperand;
  }

  //@}

private:
  //! @name Private Members
  //@{

  /**
   * \brief The modifying operand.
   */
  OperandType _modifyingOperand;
  
  //@}
};

/**
 * \brief A simple function object that subtracts
 * a specififed value from the given arguement
 * in the runFunction function.
 *
 * Please see the SubtractionFunctionXMLConverter for details on the XML
 * representation of this fucntion object.
 */
template<class OperandType>
class SubtractionFunction :
  public SimpleFunctionObject<OperandType>
{
public:
  //! @name Constructors/Destructors
  //@{
  
  /**
   * \brief Constructs a SubtractionFunction.
   */
  SubtractionFunction():SimpleFunctionObject<OperandType>(){}

  /**
   * \brief Constructs a SubtractionFunction.
   *
   * @param amountToSubtract The amound to be subtracted
   * from a given arguement in the runFunction function.
   */
  SubtractionFunction(OperandType amountToSubtract):
    SimpleFunctionObject<OperandType>(amountToSubtract){}
  
  //@}

  //! @name Overridden from SimpleFunctionObject
  //@{

  /** \brief. */
  OperandType runFunction(OperandType arguement) const{
    return 
      arguement
      -
      SimpleFunctionObject<OperandType>::getModifiyingOperand();
  }

  //@}
  
  //! @name Overridden from FunctionObject
  //@{

  /** \brief. */
  std::string getTypeAttributeValue() const{
    return 
      "SubtractionFunction(" 
      + TypeNameTraits<OperandType>::name()
      +")";
  }

  //@}
};

/**
 * \brief A simple function object that adds
 * a specififed value from the given arguement
 * in the runFunction function.
 *
 * Please see the AdditionFunctionXMLConverter for details on the XML
 * representation of this fucntion object.
 */
template<class OperandType>
class AdditionFunction :
  public SimpleFunctionObject<OperandType>
{
public:
  //! @name Constructors/Destructors
  //@{

  /**
   * \brief Constructs a AdditionFunction.
   */
  AdditionFunction():SimpleFunctionObject<OperandType>(){}

  /**
   * \brief Constructs a AdditionFunction.
   *
   * @param amountToAdd The amound to be added
   * to a given arguement in the runFunction function.
   */
  AdditionFunction(OperandType amountToAdd):
    SimpleFunctionObject<OperandType>(amountToAdd){}

  //@}

  //! @name Overridden from SimpleFunctionObject
  //@{

  /** \brief. */
  OperandType runFunction(OperandType arguement) const{
    return 
      arguement
      +
      SimpleFunctionObject<OperandType>::getModifiyingOperand();
  }

  //@}
  
  //! @name Overridden from FunctionObject
  //@{

  /** \brief. */
  std::string getTypeAttributeValue() const{
    return 
      "AdditionFunction("
      + TypeNameTraits<OperandType>::name()
      +")";
  }
  
  //@}
};

/**
 * \brief A simple function object that multiplys
 * a specififed value from the given arguement
 * in the runFunction function.
 *
 * Please see the MultiplicationFunctionXMLConverter for details on the XML
 * representation of this fucntion object.
 */
template<class OperandType>
class MultiplicationFunction :
  public SimpleFunctionObject<OperandType>
{
public:
  //! @name Constructors/Destructors
  //@{

  /**
   * \brief Constructs a MultiplicationFunction.
   *
   */
  MultiplicationFunction():SimpleFunctionObject<OperandType>(){}

  /**
   * \brief Constructs a MultiplicationFunction.
   *
   * @param amountToMultiplyBy The amound to be by
   * which the given arguement in the runFunction
   * function should be multiplied.
   */
  MultiplicationFunction(OperandType amountToMultiplyBy):
    SimpleFunctionObject<OperandType>(amountToMultiplyBy){}

  //@}
  
  //! @name Overridden from SimpleFunctionObject
  //@{

  /** \brief. */
  OperandType runFunction(OperandType arguement) const{
    return 
      arguement
      *
      SimpleFunctionObject<OperandType>::getModifiyingOperand();
  }

  //@}
  
  //! @name Overridden from FunctionObject
  //@{

  /** \brief. */
  std::string getTypeAttributeValue() const{
    return "MultiplicationFunction(" + 
      TypeNameTraits<OperandType>::name()
      +")";
  }

  //@}
};

/**
 * \brief A simple function object that divides
 * a specififed value from the given arguement
 * in the runFunction function.
 *
 * Please see the DivisionFunctionXMLConverter for details on the XML
 * representation of this fucntion object.
 */
template<class OperandType>
class DivisionFunction :
  public SimpleFunctionObject<OperandType>
{
public:

  //! @name Constructors/Destructors
  //@{
 
  /**
   * \brief Constructs a DivisionFunction.
   *
   */
  DivisionFunction():SimpleFunctionObject<OperandType>(){}

  /**
   * \brief Constructs a DivisionFunction.
   *
   * @param amoundToDivideBy The amound to be by
   * which the given arguement in the runFunction
   * function should be divided.
   */
  DivisionFunction(OperandType amountToDivideBy):
    SimpleFunctionObject<OperandType>(amountToDivideBy){}

  //@}
  
  //! @name Overridden from SimpleFunctionObject
  //@{

  /** \brief. */
  OperandType runFunction(OperandType arguement) const{
    return 
      arguement
      /
      SimpleFunctionObject<OperandType>::getModifiyingOperand();
  }

  //@}
  
  //! @name Overridden from FunctionObject
  //@{

  /** \brief. */
  std::string getTypeAttributeValue() const{
    return 
      "DivisionFunction(" 
      + TypeNameTraits<OperandType>::name()
      +")";
  }

  //@}
  
};

} // namespace Teuchos


#endif

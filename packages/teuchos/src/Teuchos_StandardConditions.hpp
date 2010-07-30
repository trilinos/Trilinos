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



#ifndef TEUCHOS_STANDARDCONDITION_HPP_
#define TEUCHOS_STANDARDCONDITION_HPP_
#include "Teuchos_Condition.hpp"
#include "Teuchos_InvalidConditionException.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Teuchos{

/**
 * An abstract parent class for all Binary Logic Conditions.
 * Binary Logic Conditions return the result of performing some
 * Logical operation on a set of conditions. Note that although the
 * name implies the evaluation of two conditions, Binary Logic Conditions
 * can actually evaluate and arbiturary number of conditions.
 */
class BinaryLogicalCondition : public Condition{
public:
	/**
	 * Constructs a BinaryLogicCondition
	 *
	 * @param conditions The conditions to be evaluated.
	 */
	BinaryLogicalCondition(ConditionList& conditions);

	/**
	 * Deconstructor for a BinaryLogicCondition
	 */
	virtual ~BinaryLogicalCondition(){}
	
	/**
	 * Adds a Condition to the list of conditions that will
	 * be evaluated by this Binary Logica Condition.
	 *
	 * @param toAdd The condition to be added to the list of
	 * conditions this Binary Logic Condition will evaluate.
	 */
	virtual void addCondition(Teuchos::RCP<Condition> toAdd);

	virtual bool isConditionTrue()=0;

	bool containsAtLeasteOneParameter();

	Dependency::ParameterParentMap getAllParameters();

protected:
	/*
	 * A list of conditions on which to perform some logic operation.
	 */
	ConditionList conditions;
};

/**
 * A Binary Logic Condition that returns the result
 * or perfroming a logical OR on the conditions.
 */
class OrCondition : public BinaryLogicalCondition{
public:
	/**
	 * Constructs an Or Condition
	 *
	 * @param conditions The conditions to be evaluated.
	 */
	OrCondition(ConditionList& conditions);

	/**
	 * Deconstructs an Or Condition.
	 */
	virtual ~OrCondition(){}

	bool isConditionTrue();
};

/**
 * A Binary Logic Condition that returns the result
 * or perfroming a logical AND on the conditions.
 */
class AndCondition : public BinaryLogicalCondition{
public:
	/**
	 * Constructs an And Condition
	 *
	 * @param conditions The conditions to be evaluated.
	 */
	AndCondition(ConditionList& conditions);

	/**
	 * Deconstructs an And Condition.
	 */
	virtual ~AndCondition(){}

	bool isConditionTrue();
};

/**
 * A Binary Logic Condition that returns the result
 * or perfroming a logical EQUALS on the conditions.
 */
class EqualsCondition : public BinaryLogicalCondition{
public:
	/**
	 * Constructs an Equals Condition
	 *
	 * @param conditions The conditions to be evaluated.
	 */
	EqualsCondition(ConditionList& conditions);

	/**
	 * Deconstructs an Equals Condition.
	 */
	virtual ~EqualsCondition(){}

	bool isConditionTrue();
};

/**
 * A Not condition returns the result of
 * performing a logical NOT on a given
 * condition.
 */
class NotCondition : public Condition{
public:
	/**
	 * Constructs a Not Condition
	 *
	 * @param condition The condition to be evaluated.
	 */
	NotCondition(Teuchos::RCP<Condition> condition);

	/**
	 * Deconstructs a Not Condition.
	 */
	virtual ~NotCondition(){}

	bool isConditionTrue();

	bool containsAtLeasteOneParameter();

	Dependency::ParameterParentMap getAllParameters();

private:
	/**
	 * The condition on which to perfrom the logical NOT.
	 */
	Teuchos::RCP<Condition> condition;
};

/**
 * An Abstract Base class for all ParameterConditions.
 * A Parmaeter Condition examines the value of a given
 * parameter and returns a bool based on the condition of
 * that value.
 */
class ParameterCondition : public Condition{
public:

   /**
 	* Constructs a Parameter Condition.
	*
	* @param parameterName The name of the parameter to be evaluated.
	* @param parentList The parent Parameter List of the parameter to be evaluated.
	* @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
	* results in a true or when the evaluation results in a false. When set to true, if the parameter
	* evaluates to true then the condition will evaluate to true. If set to false if the parameter
	* evaluates to false, then the condition will evaluatae to true.
 	*/
	ParameterCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue);

	virtual ~ParameterCondition(){}

	virtual bool isConditionTrue()=0;

	bool containsAtLeasteOneParameter();

	Dependency::ParameterParentMap getAllParameters();

protected:
	/**
	 * Name of parameter to be evaluated.
	 */
	std::string parameterName;

	/**
	 * Parent List of the parameter to be evaluated.
	 */
	Teuchos::RCP<Teuchos::ParameterList> parentList;

	/**
	 * Wether or not the condition should evaluate to true if the parameter evaluated to true.
	 */
	bool whenParamEqualsValue;

	/**
	 * A pointer to the actual parameter to be evaluated.
	 */
	Teuchos::ParameterEntry* parameter;
};

/**
 * A String Condition is a Parameter Condition that evaluates
 * whether or not a string parameter has taken on a particular
 * value or set of values.
 */
class StringCondition : public ParameterCondition{
public:
	/**
	 * Convience typedef representing an array of strings.
	 */
	typedef Teuchos::Array<std::string> ValueList; 

   /**
 	* Constructs a String Condition.
	*
	* @param parameterName The name of the parameter to be evaluated.
	* @param parentList The parent Parameter List of the parameter to be evaluated.
	* #param value The value to compare the parameter's value against.
	* @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
	* results in a true or when the evaluation results in a false. When set to true, if the parameter
	* evaluates to true then the condition will evaluate to true. If set to false if the parameter
	* evaluates to false, then the condition will evaluatae to true.
 	*/
	StringCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, std::string value, bool whenParamEqualsValue=true);

   /**
 	* Constructs a String Condition.
	*
	* @param parameterName The name of the parameter to be evaluated.
	* @param parentList The parent Parameter List of the parameter to be evaluated.
	* #param values The values to compare the parameter's value against.
	* @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
	* results in a true or when the evaluation results in a false. When set to true, if the parameter
	* evaluates to true then the condition will evaluate to true. If seet to false if the parameter
	* evaluates to false, then the condition will evaluatae to true.
 	*/
	StringCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, ValueList values, bool whenParamEqualsValue=true);

	virtual ~StringCondition(){}

	bool isConditionTrue();

private:
	/**
	 * A list of values against which to evaluate the parameter's value.
	 */
	ValueList values;
};

/**
 * A Number Condition is a Parameter Condition that evaluates
 * whether or not a number parameter is greater 0. If the parameter is
 * greater than 0 this is interperted as the condition being "true".
 * Otherwise the oncidiont is interperted as false.
 */
template<class T>
class NumberCondition : public ParameterCondition{
public:
	/**
 	* Constructs a Number Condition.
	*
	* @param parameterName The name of the parameter to be evaluated.
	* @param parentList The parent Parameter List of the parameter to be evaluated.
	* @param func A function to run the value of the parameter through. If the function returns a value
	* greater than 0, this will be interperted as the condition being "true". If the 
	* function returns a value of 0 or less, this will be interperted as the condition being false.
	* @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
	* results in a true or when the evaluation results in a false. When set to true, if the parameter
	* evaluates to true then the condition will evaluate to true. If seet to false if the parameter
	* evaluates to false, then the condition will evaluatae to true.
 	*/
	NumberCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, T (*func)(T)=0, bool whenParamEqualsValue=true):
		ParameterCondition(parameterName, parentList, whenParamEqualsValue), func(func)
	{
		if(!parameter->isType<int>()
		&& !parameter->isType<short>()
		&& !parameter->isType<double>()
		&& !parameter->isType<float>()){
			std::string expectedTypeName = typeid(T).name();
			throw InvalidConditionException("The parameter of a "
			"Number Condition must be of a supported number type!\n"
			"Expected type: " + expectedTypeName + "\n" 
			"Actual type: " + parameter->getAny().typeName());
		}
	}
		
	/**
 	* Constructs a Number Condition.
	*
	* @param parameterName The name of the parameter to be evaluated.
	* @param parentList The parent Parameter List of the parameter to be evaluated.
	* @param func A function to run the value of the parameter through. If the function returns a value
	* greater than 0, this will be interperted as the parameter's current state being "true". If the 
	* function returns a value of 0 or less, this will be interperted as the parameter's current state
	* being "false".
	* @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
	* results in a true or when the evaluation results in a fals. When set to true, if the parameter
	* evaluates to true then the condition will evaluate to true. If seet to false if the parameter
	* evaluates to false, then the condition will evaluatae to true.
 	*/
	NumberCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue=true):
		ParameterCondition(parameterName, parentList, whenParamEqualsValue), func(0){}

	virtual ~NumberCondition(){}

	bool isConditionTrue(){
		bool toReturn = (runFunction(Teuchos::getValue<T>(*parameter)) > 0);
		return whenParamEqualsValue ? toReturn : !toReturn;
	}

private:
	T (*func)(T); 	

	/**
	 * Runs the function associated with this condition and
	 * returns the result.
	 *
	 * @param argument The value upon which to run the function.
	 */
	T runFunction(T argument) const{
		if(func !=0)
			return (*func)(argument);
		else
			return argument;
	}	
};

/**
 * A Bool Condition is a Parameter Condition that evaluates
 * whether or not a Boolean parameter is ture.
 * */
class BoolCondition : public ParameterCondition{
public:
   /**
 	* Constructs a Bool Condition.
	*
	* @param parameterName The name of the parameter to be evaluated.
	* @param parentList The parent Parameter List of the parameter to be evaluated.
	* @param whenParamEqualsValue Indicates whether the condition should be true when the evaluation
	* results in a true or when the evaluation results in a false. When set to true, if the parameter
	* evaluates to true then the condition will evaluate to true. If set to false if the parameter
	* evaluates to false, then the condition will evaluatae to true.
 	*/
	BoolCondition(std::string parameterName, Teuchos::RCP<Teuchos::ParameterList> parentList, bool whenParamEqualsValue=true);

	virtual ~BoolCondition(){}

	bool isConditionTrue();
};

}
#endif //TEUCHOS_STANDARDCONDITION_HPP_

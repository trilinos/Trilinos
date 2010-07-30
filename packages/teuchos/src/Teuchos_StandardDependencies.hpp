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



#ifndef TEUCHOS_STANDARDDEPENDCIES_HPP_
#define TEUCHOS_STANDARDDEPENDCIES_HPP_
#include "Teuchos_Dependency.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Condition.hpp"

namespace Teuchos{

/**
 * An abstract parent class for all visual dependencies.
 *
 * IMPORTANT NOTE:
 * If a parameter becomes hidden, it's validity will not be checked by the GUI. This means that it
 * is indeed possible for the GUI to return a non-valid ParameterList. Make sure that you program code takes
 * this into account.
 */
class VisualDependency : public Dependency{
public:
	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 */
	VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf=true);

	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 */
	VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, bool showIf=true);

	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependees A map containing all the dependee Parameters associated with their parent ParameterLists.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 */
	VisualDependency(ParameterParentMap dependees, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf=true);

	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependees A map containing all the dependee Parameters associated with their parent ParameterLists.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 */
	VisualDependency(ParameterParentMap dependees, ParameterParentMap dependents, bool showIf=true);

	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 */
	VisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents);

	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependees A map containing all the dependee Parameters associated with their parent ParameterLists.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 */
	VisualDependency(ParameterParentMap dependees, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList);

	/**
	 * Constructs a VisualDependency.
	 *
	 * @param dependees A map containing all the dependee Parameters associated with their parent ParameterLists.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 */
	VisualDependency(ParameterParentMap dependees, ParameterParentMap dependents);

	/**
	 * Desctructor
	 *
	 * Simply declaring the descrutor as virtual.
	 */
	virtual ~VisualDependency(){}

	virtual void evaluate() = 0;

	/**
	 * Determines whether or not the dependent is currently visible.
	 */
	bool isDependentVisible();

protected:

	/**
	 * Whether or not the dependent is currently visible.
	 */
	bool dependentVisible;

	/**
	 * Whether or not to show the dependent if the dependee is set to the value.
	 */
	bool showIf;

private:
	virtual void validateDep() = 0;
};


/**
 * An abstract base class for all validator dependencies.
 */
class ValidatorDependency : public Dependency{
public:
	/**
	 * Constructs a ValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 */
	ValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList);

	/**
	 * Constructs a ValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 */
	ValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents);

	/**
	 * Desctructor
	 *
	 * Simply declaring the descrutor as virtual.
	 */
	virtual ~ValidatorDependency(){}

	virtual void evaluate() = 0;

private:
	virtual void validateDep() = 0;

};

/**
 * A string visual depdencies says the following about the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a particular value, the dependent may or may not be displayed to the user in a GUI.
 * 
 * The dependee of a StringVisualDependency must be of type string and can't be an array. The dependent may be any type of
 * parameter or parameter list.
 */
class StringVisualDependency : public VisualDependency{
public:
	/**
	 * Convience typedef representing an array of strings.
	 */
	typedef Teuchos::Array<std::string> ValueList; 

	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param value The value of the depndee that affects the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to the same value as specified by the value parameter.
	 * If false, the dependent will be shown only when the dependee is set to a value other than the one specified by the value parameter.
	 */
	StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, std::string value, bool showIf=true);

	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param value The value of the depndee that affects the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to the same value as specified by the value parameter.
	 * If false, the dependent will be shown only when the dependee is set to a value other than the one specified by the value parameter.
	 */
	StringVisualDependency(std::string dependeeName, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> parentList, 
	std::string value, bool showIf=true);

	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param values The values of the depndee that affect the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to one of the values specified by the values parameter.
	 * If false, the dependent will be shown only when the dependee is set to a value other than the ones specified by the values parameter.
	 */
	StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, const ValueList& values, bool showIf=true);

	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param values The values of the depndee that affect the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to the one of the values as specified by the values parameter.
	 * If false, the dependent will be shown only when the dependee is set to a value other than the ones specified by the values parameter.
	 */
	StringVisualDependency(std::string dependeeName, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> parentList, 
	const ValueList& values, bool showIf=true);

	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param value The value of the depndee that affects the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to the same value as specified by the value parameter.
	 * If false, the dependent will be shown only when the dependee is set to a value other than the one specified by the value parameter.
	 */
	StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, std::string value, bool showIf=true);

	/**
	 * Constructs a StringVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param values The values of the depndee that affect the visiblity of the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is set to one of the values specified by the values parameter.
	 * If false, the dependent will be shown only when the dependee is set to a value other than the ones specified by the values parameter.
	 */
	StringVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, const ValueList& values, bool showIf=true);

	void evaluate();

private:
	/**
	 * The value used to deteremine the visiblity of the dependent.
	 */
	const ValueList values;

	void validateDep();
};

/**
 * A bool visual dependency says the following about the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee is true or false, the dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a BoolVisualDependency must be of type bool and can't be an array. The dependent may be any type of parameter
 * or parameter list.
 */
class BoolVisualDependency : public VisualDependency{
public:
	/**
	 * Constructs a BoolVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 * If false, the dependent will be shown only when the dependee is false.
	 */
	BoolVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf=true);

	/**
	 * Constructs a BoolVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 * If false, the dependent will be shown only when the dependee is false.
	 */
	BoolVisualDependency(std::string dependeeName, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> parentList, 
	bool showIf=true);

	/**
	 * Constructs a BoolVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param showIf When true, the depndent will be be shown if the dependee is true.
	 * If false, the dependent will be shown only when the dependee is false.
	 */
	BoolVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, bool showIf=true);

	void evaluate();

	void validateDep();
};

/**
 * A condition visual dependency says the following about the relationship between elements in a Parameter List:
 * Depending on whether or not the dependee(s) statisfy a particual condition, the dependent may or may not be displayed to the user in a GUI.
 *
 * Conditoin Visual Dependencies are unique in that via the Condition class, they allow for multiple dependees.
 * The dependee(s) of a ConditionVisualDependency must be expressed as a Condition and are subject to the consquential constraints. The dependent may be any type of parameter
 * or parameter list.
 */
class ConditionVisualDependency : public VisualDependency{
public:
	/**
	 * Constructs a ConditionVisualDependency.
	 *
	 *
	 * @param condition The condition that must be satisfied in order to display the dependent
	 * parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param showIf When true, the depndent will be be shown if the condition is true.
	 * If false, the dependent will be shown only when the dependee is false.
	 */
	ConditionVisualDependency(Teuchos::RCP<Condition> condition,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, bool showIf=true);

	/**
	 * Constructs a BoolVisualDependency.
	 *
	 * @param condition The condition that must be satisfied in order to display the dependent
	 * parameter.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param showIf When true, the depndent will be be shown if the condition is true.
	 * If false, the dependent will be shown only when the dependee is false.
	 */
	ConditionVisualDependency(Teuchos::RCP<Condition> condition, ParameterParentMap dependents, bool showIf=true);

	void evaluate();

	void validateDep();

private:
	/**
	 * The Condition to determine whether or not the dependent is displayed.
	 */
	Teuchos::RCP<Condition> condition;
};

/**
 * A number visual dependency says the following about the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a certain value, the dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a NumberVisualDependency must be a number type and can't be an array. The dependent may be any type of parameter
 * or parameter list.
 */
template <class S>
class NumberVisualDependency : public VisualDependency{
public:
	/**
	 * Constructs a NumberVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param func A function that takes the dependees value, does some calculations on it, and then
	 * returns a value. If this value is greater than 0, the dependent is show. If the value returned is
	 * less than or equal to zero, the dependent is not shown. If no fuction is specified, the direct
	 * value of the dependee will be used to determine the dependents visibility in a similar fashion (postive
	 * numbers causing the dependent to be displayed and 0 or negative numbers cuasing the dependent to be
	 * hidden).
	 */
	NumberVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, S (*func)(S) =0)
		:VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList), func(func)
	{
		validateDep();
	}

	/**
	 * Constructs a NumberVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param func A function that takes the dependees value, does some calculations on it, and then
	 * returns a value. If this value is greater than 0, the dependent is show. If the value returned is
	 * less than or equal to zero, the dependent is not shown. If no fuction is specified, the direct
	 * value of the dependee will be used to determine the dependents visibility in a similar fashion (postive
	 * numbers causing the dependent to be displayed and 0 or negative numbers cuasing the dependent to be
	 * hidden).
	 */
	NumberVisualDependency(std::string dependeeName, std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> parentList, S 
	(*func)(S) =0)
		:VisualDependency(dependeeName, parentList, dependentName, parentList), func(func)
	{
		validateDep();
	}

	/**
	 * Constructs a NumberVisualDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param func A function that takes the dependees value, does some calculations on it, and then
	 * returns a value. If this value is greater than 0, the dependent is show. If the value returned is
	 * less than or equal to zero, the dependent is not shown. If no fuction is specified, the direct
	 * value of the dependee will be used to determine the dependents visibility in a similar fashion (postive
	 * numbers causing the dependent to be displayed and 0 or negative numbers cuasing the dependent to be
	 * hidden).
	 */
	NumberVisualDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, S (*func)(S) =0)
		:VisualDependency(dependeeName, dependeeParentList, dependents), func(func)
	{
		validateDep();
	}

	void evaluate(){
		S functionVal = runFunction(getFirstDependeeValue<S>());
		functionVal > 0 ? dependentVisible = true : dependentVisible = false;
	}

private:
	/**
	 * the function used to determine the
	 * visibility of the dependent.
	 */
	S (*func)(S);

	/**
	 * Run the function on the argument and returns the value of the fucntion. If no function is specified,
	 * the argument is simple returned.
	 *
	 * @param argument the value to use as an argument for the function.
	 * @return the result of running the function with the value. If no function is specified,
	 * the argument is simple returned.
	 */
	S runFunction(S argument) const{
		if(func !=0)
			return (*func)(argument);
		else
			return argument;
	}	

	void validateDep(){
		/*
		 * This error should never get thrown, unless someone
		 * is doing something wonky in a sublcass.
		 */
		if(dependees.size() != 1){
			std::stringstream out;
			out << dependees.size();
			throw InvalidDependencyException("Uh oh. Looks like you tried to make a " 
			"Number Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
			"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
			"be doing some subclassing you sly-dog ;)\n\n" 
			"Error: A Number Visual Dependency must have exactly 1 dependee. " 
			"You have tried to assign it " + out.str()+ " dependees.\n" 
			"Dependees: " + getDependeeNamesString() + "\n" 
			"Dependents: " + getDependentNamesString());
		}
		const Teuchos::ParameterEntry* dependee = getFirstDependee();
		if(!dependee->isType<int>()
		&& !dependee->isType<short>()
		&& !dependee->isType<double>()
		&& !dependee->isType<float>()){
			throw InvalidDependencyException("The dependee of a "
			"Number Visual Dependency must be of a supported number type!\n"
			"Problem dependee: " + getFirstDependeeName() + "\n"
			"Actual type: " + dependee->getAny().typeName() + "\n"
			"Dependees: " + getDependeeNamesString() + "\n"
			"Dependents: " + getDependentNamesString());
		}
	}
};

/**
 * A NumberValidatorAspectDependency says the following about the relationship between two parameters in a dependent parameter list:
 * depending of the value of the dependee a particular aspect of the dependents validator will have a certain value.
 *
 * A NumberValidatorAspectDependency must have the following characteristics:
 * <ul>
 * 	<li>The dependee must be of a number type and can't be an array.</li>
 * 	<li>The dependent must have an enhancednumberparameter validator or an arrayNumberValidator.</li>
 * 	<li>The validator specified in the constructor must be the same validator being used on the dependent parameter.</li>
 * 	<li>The template type, dependee type, and dependent type must all be the same.</li>
 * <ul>
 *
 * This dependency can have an interesting effect on your program. The class will modifies
 * dependents validator. If that same validator is used more than once, every other 
 * parameter that used that validator will also see the change. Make sure that you are
 * aware of this when making your validators.
 */
template <class S>
class NumberValidatorAspectDependency : public Dependency{
public:
	/**
	 * An enum specifying the aspect of the
	 * validator that should be modified
	 */
	enum ValidatorAspect{
		Min,
		Max,
		Step
	};

	/**
	 * Constructs a NumberValidatorDependency
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, 
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberValidatorAspectDep),
		aspect(aspect),
		validator(validator),
		func(func)

	{
		validateDep();
	}

	/**
	 * Constructs a NumberValidatorDependency
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, std::string dependentName, 
	Teuchos::RCP<Teuchos::ParameterList> parentList, 
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:NumberValidatorAspectDependency(dependeeName, parentList, dependentName, parentList, validator, aspect, func) {}

	/**
	 * Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, 
	Teuchos::RCP<Teuchos::ArrayNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberValidatorAspectDep),
		aspect(aspect),
		validator(validator->getPrototype()),
		func(func)
	{
		validateDep();
	}

	/**
	 * Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, std::string dependentName, 
	Teuchos::RCP<Teuchos::ParameterList> parentList,
	Teuchos::RCP<Teuchos::ArrayNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, parentList, dependentName, parentList, Dependency::NumberValidatorAspectDep),
		aspect(aspect),
		validator(validator),
		func(func)

	{
		validateDep();
	}

	/**
	 * Constructs a NumberValidatorDependency
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, 
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, dependeeParentList, dependents, Dependency::NumberValidatorAspectDep),
		aspect(aspect),
		validator(validator),
		func(func)

	{
		validateDep();
	}

	/**
	 * Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param validator The validator whose aspect will change. 
	 * @param aspect The aspect of the validator that should change.
	 * @param func A function specifying how the value of the validators
	 * aspect should be calculated from the dependees value.
	 */
	NumberValidatorAspectDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents,
	Teuchos::RCP<Teuchos::ArrayNumberValidator<S> > validator,
	ValidatorAspect aspect, S (*func)(S) =0)
		:Dependency(dependeeName, dependeeParentList, dependents, Dependency::NumberValidatorAspectDep),
		aspect(aspect),
		validator(validator->getPrototype()),
		func(func)
	{
		validateDep();
	}

	void evaluate(){
		S newAspectValue = runFunction(getFirstDependeeValue<S>());
		switch(aspect){
			case NumberValidatorAspectDependency<S>::Min:
				validator->setMin(newAspectValue);
				break;
			case NumberValidatorAspectDependency<S>::Max:
				validator->setMax(newAspectValue);
				break;
			case NumberValidatorAspectDependency<S>::Step:
				validator->setStep(newAspectValue);
				break;
		}
	}

private:
	/**
	 * The aspect of the validator to be modified.
	 */
	ValidatorAspect aspect;

	/**
	 * The validator to be modified.
	 */
	Teuchos::RCP<Teuchos::EnhancedNumberValidator<S> > validator;

	/**
	 * The function used to calculate the new value of the
	 * aspect of the validator.
	 */
	S (*func)(S);
	
	/**
	 * Runs the dependency's function on the given argument and returns
	 * the value that function returns.
	 *
	 * @param The value to run the function on.
	 * @return The value the function returned.
	 */
	S runFunction(S argument) const{
		if(func !=0)
			return (*func)(argument);
		else
			return argument;
	}	

	void validateDep(){
		/*
		 * This error should never get thrown, unless someone
		 * is doing something wonky in a sublcass.
		 */
		if(dependees.size() != 1){
			std::stringstream out;
			out << dependees.size();
			throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
			"Number Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
			"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
			"be doing some subclassing you sly-dog ;)\n\n" 
			"Error: A Number Visual Dependency must have exactly 1 dependee. " 
			"You have tried to assign it "+ out.str() + " dependees.\n" 
			"Dependees: " + getDependeeNamesString() + "\n" 
			"Dependents: " + getDependentNamesString());
		}

		const Teuchos::ParameterEntry* dependee = getFirstDependee();
		std::string dependeeName = getFirstDependeeName();
		if(!dependee->isType<int>()
		&& !dependee->isType<short>()
		&& !dependee->isType<double>()
		&& !dependee->isType<float>()){
			throw InvalidDependencyException("The dependee of a "
			"Number Validator Aspect Dependency must be of a supported number type!\n"
			"Problem dependee: " + dependeeName + "\n"
			"Actual type: " + dependee->getAny().typeName() + "\n"
			"Dependents: " + getDependentNamesString());
		}

		if(typeid(S) != dependee->getAny().type()){
			throw InvalidDependencyException("The dependee type and EnhancedNumberValidator "
			"template type must all be the same for a Number Validator Aspect Dependency.\n"
			"Problem Dependee: " + dependeeName + "\n"
			"Dependee Type: " + dependee->getAny().typeName() + "\n"
			"Validator Template Type: " + typeid(S).name());
		}

		typename ParameterParentMap::const_iterator it;
		Teuchos::ParameterEntry *currentDependent;
		for(it = dependents.begin(); it != dependents.end(); ++it){ 
			currentDependent = it->second->getEntryPtr(it->first);
			if(currentDependent->validator() == Teuchos::null){
				throw InvalidDependencyException("The dependent of an "
				"Number Validator Aspect Dependency must have an EnhancedNumberValidator "
				"or an ArrayNumberValidator\n"
				"Problem dependent: " + it->first + "\n" 
				"Dependees: " + getDependeeNamesString() + "\n"
				"Dependents: " + getDependentNamesString());
			}
			if(validator != currentDependent->validator()){
				throw InvalidDependencyException("The dependent's validator and the validator specified "
				"in the constructor must be the same for a Number Validator Aspect Dependency!\n"
				"Problem dependent: " + it->first + "\n" 
				"Dependees: " + getDependeeNamesString() + "\n"
				"Dependents: " + getDependentNamesString());
			}
			if(typeid(S) != currentDependent->getAny().type()){
				throw InvalidDependencyException("The dependent type and EnhancedNumberValidator "
				"template type must all be the same for a Number Validator Aspect Dependency.\n"
				"Dependent: " + it->first + "\n" 
				"Dependent Type: " + currentDependent->getAny().typeName() + "\n"
				"Validator Template Type: " + typeid(S).name());
			}
		}
	}
};

/**
 * An ArrayLengthDependency says the following about the relationship between two parameters:
 * The length of the dependent's array depends on the value of the dependee.
 *
 * An ArrayLengthDependency must have the following characteristics:
 * <ul>
 * 	<li>The dependee must be either of type int or type short.</li>
 * 	<li>The dependent must be an array</li>
 * </ul>
 */
class NumberArrayLengthDependency : public Dependency{
public:
	/**
	 * Constructs an ArrayLengthDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param func A function specifying how the arrays length 
	 * should be calculated from the dependees value.
	 */
	NumberArrayLengthDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, int (*func)(int) = 0);

	/**
	 * Constructs an ArrayLengthDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param func A function specifying how the arrays length 
	 * should be calculated from the dependees value.
	 */
	NumberArrayLengthDependency(std::string dependeeName, std::string dependentName, 
	Teuchos::RCP<Teuchos::ParameterList> parentList, int (*func)(int) = 0);

	/**
	 * Constructs an ArrayLengthDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param func A function specifying how the arrays length 
	 * should be calculated from the dependees value.
	 */
	NumberArrayLengthDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, int (*func)(int) = 0);
	
	void evaluate();

private:
	/**
	 * The function used to calculate the new value of the
	 * arrays length.
	 */
	int (*func)(int);
	
	/**
	 * Runs the dependency's function on the given argument and returns
	 * the value that function returns.
	 *
	 * @param The value to run the function on.
	 * @return The value the function returned.
	 */
	int runFunction(int argument) const;

	/**
	 * Modifies the length of an array.
	 *
	 * @param newLength The new length the array should be.
	 * @param dependentValue The index of the dependent array that is going to be changed.
	 */
	template <class S>
	void modifyArrayLength(int newLength, Teuchos::ParameterEntry* dependentToModify);

	void validateDep();
};

/**
 * A StringValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A StringValidatorDependency must have the following characterisitics:
 * <ul>
 * 	<li>The dependee must be of type string</li>
 * 	<li>The dependent must already have a validator assigned to it.</li>
 * 	<li>All of the validators that may be applied to the dependent must match the type of the
 * 	validator that is currently on the dependent.</li>
 * </ul>
 */
class StringValidatorDependency : public ValidatorDependency{
public:
	/**
	 * Conveniece typedef
	 */
	typedef std::map<std::string, Teuchos::RCP<const Teuchos::ParameterEntryValidator> > ValueToValidatorMap;

	/**
	 * Constructs a StringValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param valuesAndValidators A map associating string values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that is not in the valuesAndValidators map,
	 * this is the validator that will be assigned to the dependent.
	 */
	StringValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,  ValueToValidatorMap valuesAndValidators, 
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator);

	/**
	 * Constructs a StringValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param valuesAndValidators A map associating string values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that is not in the valuesAndValidators map,
	 * this is the validator that will be assigned to the dependent.
	 */
	StringValidatorDependency(std::string dependeeName, std::string dependentName, 
	Teuchos::RCP<Teuchos::ParameterList> parentList,  ValueToValidatorMap valuesAndValidators, 
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator);

	/**
	 * Constructs a StringValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param valuesAndValidators A map associating string values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that is not in the valuesAndValidators map,
	 * this is the validator that will be assigned to the dependent.
	 */
	StringValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents,  ValueToValidatorMap valuesAndValidators, 
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator);


	void evaluate();

private:
	/**
	 * The default validator to be used if a request is made for a value that does not
	 * appear in the valuesAndValidators map.
	 */
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator;

	/**
	 * A map associating particular dependee values with validators that
	 * could be placed on the dependent.
	 */
	ValueToValidatorMap valuesAndValidators;
	
	void validateDep();
};

/**
 * A BoolValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A StringValidatorDependency must have the following characterisitics:
 * <ul>
 * 	<li>The dependee must be of type bool</li>
 * 	<li>The dependent must already have a validator assigned to it.</li>
 * 	<li>The "true" and "false" validators that may be applied to the dependent must match the type of the
 * 	validator that is currently on the dependent.</li>
 * </ul>
 */
class BoolValidatorDependency : public ValidatorDependency{
public:
	/**
	 * Constructs a BoolValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
	 * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
	 */
	BoolValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator);

	/**
	 * Constructs a BoolValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
	 * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
	 */
	BoolValidatorDependency(std::string dependeeName, 
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> parentList,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator);

	/**
	 * Constructs a BoolValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
	 * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
	 */
	BoolValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator,
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> falseValidator);

	void evaluate();

private:
	/**
	 * The validators to be used when the dependee is either
	 * true or false.
	 */
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> trueValidator, falseValidator;

	void validateDep();
};


/**
 * A RangeValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A RangeValidatorDependency achieves this by associating ranges of numerical values with validators.
 * If the dependees value falls within the one of the ranges, the validator associated with the range is
 * used on the dependent. All ranges are inclusive.
 *
 * A RangeValidatorDependency must have the following characterisitics:
 * <ul>
 * 	<li>The dependee must be a number type</li>
 * 	<li>The dependent must already have a validator assigned to it.</li>
 * 	<li>All of the validators that may be applied to the dependent must match the type of the
 * 	validator that is currently on the dependent.</li>
 * </ul>
 */
template<class S>
class RangeValidatorDependency : public ValidatorDependency{
public:
	/**
	 * Convenience typedef
	 */
	typedef std::pair<S,S> Range;
	typedef std::map<Range, Teuchos::RCP<const Teuchos::ParameterEntryValidator> > RangeToValidatorMap;

	/**
	 * Constructs a RangeValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependentName The name of the dependent parameter.
	 * @param dependentParentList The ParameterList containing the dependent.
	 * @param rangesAndValidators A map associating ranges of values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that does not fall within any of the ranges in
	 * the rangesAndValidators map, this is the validator that will be assigned to the dependent.
	 */
	RangeValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	std::string dependentName, Teuchos::RCP<Teuchos::ParameterList> dependentParentList, RangeToValidatorMap rangesAndValidators,
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
		:ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList),
		defaultValidator(defaultValidator),
		rangesAndValidators(rangesAndValidators)
	{
		validateDep();
	}

	/**
	 * Constructs a RangeValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependentName The name of the dependent parameter.
	 * @param parentList The ParameterList containing the dependent and the dependee.
	 * @param rangesAndValidators A map associating ranges of values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that does not fall within any of the ranges in
	 * the rangesAndValidators map, this is the validator that will be assigned to the dependent.
	 */
	RangeValidatorDependency(std::string dependeeName, std::string dependentName, 
	Teuchos::RCP<Teuchos::ParameterList> parentList, RangeToValidatorMap rangesAndValidators,
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
		:ValidatorDependency(dependeeName, parentList, dependentName, parentList),
		defaultValidator(defaultValidator),
		rangesAndValidators(rangesAndValidators)
	{
		validateDep();
	}

	/**
	 * Constructs a RangeValidatorDependency.
	 *
	 * @param dependeeName The name of the dependee parameter.
	 * @param dependeeParentList The ParameterList containing the dependee.
	 * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
	 * @param rangesAndValidators A map associating ranges of values with ParameterEntryValidators. This will be used
	 * to deteremine what type of validator should be applied to the dependent based on the dependees value.
	 * @param defaultValidator If a value is entered in the dependee that does not fall within any of the ranges in
	 * the rangesAndValidators map, this is the validator that will be assigned to the dependent.
	 */
	RangeValidatorDependency(std::string dependeeName, Teuchos::RCP<Teuchos::ParameterList> dependeeParentList,
	ParameterParentMap dependents, RangeToValidatorMap rangesAndValidators,
	Teuchos::RCP<Teuchos::ParameterEntryValidator> defaultValidator)
		:ValidatorDependency(dependeeName, dependeeParentList, dependents),
		defaultValidator(defaultValidator),
		rangesAndValidators(rangesAndValidators)
	{
		validateDep();
	}

	void evaluate(){
		typename RangeToValidatorMap::const_iterator it;
		S dependeeValue = getFirstDependeeValue<S>();
		for(it = rangesAndValidators.begin(); it != rangesAndValidators.end(); ++it){
			S min = it->first.first;
			S max = it->first.second;
			if(dependeeValue >= min && dependeeValue <=max){
				typename ParameterParentMap::iterator it2;
				Teuchos::ParameterEntry *currentDependent;
				for(it2 = dependents.begin(); it2 != dependents.end(); ++it2){ 
					currentDependent = it2->second->getEntryPtr(it2->first);
					currentDependent->setValidator(it->second);
				}
				return;
			}
		}

	}
private:
	/**
	 * The default validator
	 */
	Teuchos::RCP<const Teuchos::ParameterEntryValidator> defaultValidator;

	/**
	 * A map associating ranges with validators.
	 */
	RangeToValidatorMap rangesAndValidators;
	
	void validateDep(){
		/*
		 * This error should never get thrown, unless someone
		 * is doing something wonky in a sublcass.
		 */
		if(dependees.size() != 1){
			std::stringstream out;
			out << dependees.size();
			throw InvalidDependencyException("Uh oh. Looks like you tried to make a "
			"Number Visual Dependency doesn't have exactly one dependee. This is kind of a problem. "
			"You should probably take a look into it. I'm actually amazed you even threw this error. You must "
			"be doing some subclassing you sly-dog ;)\n\n" 
			"Error: A Number Visual Dependency must have exactly 1 dependee. "
			"You have tried to assign it "+ out.str() + " dependees.\n"
			"Dependees: " + getDependeeNamesString() + "\n" 
			"Dependents: " + getDependentNamesString());
		}

		const Teuchos::ParameterEntry* dependee = getFirstDependee();
		std::string dependeeName = getFirstDependeeName();
		if(!dependee->isType<int>()
		&& !dependee->isType<short>()
		&& !dependee->isType<double>()
		&& !dependee->isType<float>()){
			throw InvalidDependencyException("The dependee of a "
			"Range Validator Dependency must be of a supported number type!\n"
			"Problem dependee: " + dependeeName + "\n"
			"Actual type: " + dependee->getAny().typeName() + "\n"
			"Dependents: " + getDependentNamesString());
		}

		typename RangeToValidatorMap::const_iterator it;
		for(it = rangesAndValidators.begin(); it != rangesAndValidators.end(); ++it){
			typename ParameterParentMap::const_iterator it2;
			Teuchos::ParameterEntry *currentDependent;
			for(it2 = dependents.begin(); it2 != dependents.end(); ++it2){ 
				currentDependent = it2->second->getEntryPtr(it2->first);
				if(typeid(*(currentDependent->validator().get())) != typeid(*(it->second.get()))){
					throw InvalidDependencyException("The validator of a dependent of a "
					"Range Validator Dependency must be the same type as all of the validators "
					"in the rangesAndValidators map.\n"
					"Note this means that the dependent must have an initial validator.\n"
					"Problem dependent: " + it2->first + "\n"
					"Validator Type: " + typeid(*(currentDependent->validator())).name() + "\n"
					"One of the validators in the rangesAndValidators map is of type: " + typeid(*(it->second)).name());
				}
			}
		}
	}
};


}
#endif //TEUCHOS_STANDARDDEPENDCIES_HPP_

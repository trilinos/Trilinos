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

/*! \file Teuchos_StandardDependencies.hpp
    \brief A collection of standard dependencies.
*/

#include "Teuchos_Dependency.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Condition.hpp"

namespace Teuchos{

/**
 * \brief An abstract parent class for all visual dependencies.
 *
 * IMPORTANT NOTE:
 * If a parameter becomes hidden, it's validity will not necessarily be checked. This means that it
 * is indeed possible for a non-valid ParameterList to occur. Make sure that you program code takes
 * this into account.
 */
class VisualDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependees A map containing all the dependee Parameters associated with their parent ParameterLists.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   */
  VisualDependency(ParameterParentMap dependees, std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependees A map containing all the dependee Parameters associated with their parent ParameterLists.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   */
  VisualDependency(ParameterParentMap dependees, ParameterParentMap dependents, bool showIf=true);

  /**
   * \brief Desctructor
   *
   * Simply declaring the descrutor as virtual.
   */
  virtual ~VisualDependency(){}
  
  //@}

  //! @name Attribute/Query Methods 
  //@{

  /**
   * \brief Get the state of the dependee in order to evaluate the
   * dependency.
   *
   * @return The state of the dependee.
   */
  virtual bool getDependeeState() const = 0;
  
  /**
   * \brief Determines whether or not the dependent is currently visible.
   */
  bool isDependentVisible() const;

  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  virtual void evaluate(){
    if((getDependeeState() && showIf_) || (!getDependeeState() && !showIf_)){
      dependentVisible_ = true;
    }
    else{
      dependentVisible_ = false;
    }
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief Whether or not the dependent is currently visible.
   */
  bool dependentVisible_;

  /**
   * \brief Whether or not to show the dependent if the dependee is set to the value.
   */
  bool showIf_;
  
  //@}

};


/**
 * \brief An abstract base class for all validator dependencies.
 */
class ValidatorDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a ValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   */
  ValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList);

  /**
   * \brief Constructs a ValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   */
  ValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents);

  /**
   * \brief Desctructor
   *
   * Simply declaring the descrutor as virtual.
   */
  virtual ~ValidatorDependency(){}
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  virtual void evaluate() = 0;
  
  //@}

};

/**
 * \brief A string visual depdencies says the following about the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a particular value, the dependent may or may not be displayed to the user in a UI.
 * 
 * The dependee of a StringVisualDependency must be of type string and can't be an array. The dependent may be any type of
 * parameter or parameter list.
 */
class StringVisualDependency : public VisualDependency{

public:

  /** \name Public types */
  //@{

  /**
   * Convience typedef representing an array of strings.
   */
  typedef Array<std::string> ValueList; 
  
  //@}

  /** \name Constructors/Destructor */
  //@{

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
  StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, std::string value, bool showIf=true);

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
  StringVisualDependency(std::string dependeeName, std::string dependentName, RCP<ParameterList> parentList, 
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
  StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, const ValueList& values, bool showIf=true);

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
  StringVisualDependency(std::string dependeeName, std::string dependentName, RCP<ParameterList> parentList, 
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
  StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
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
  StringVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, const ValueList& values, bool showIf=true);
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  bool getDependeeState() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * The value used to deteremine the visiblity of the dependent.
   */
  const ValueList values_;
  
  //@}
  
  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const;
  
  //@}
  
};

/**
 * \brief A bool visual dependency says the following about the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee is true or false, the dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a BoolVisualDependency must be of type bool and can't be an array. The dependent may be any type of parameter
 * or parameter list.
 */
class BoolVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

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
  BoolVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf=true);

  /**
   * Constructs a BoolVisualDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependentName The name of the dependent parameter.
   * @param parentList The ParameterList containing the dependent and the dependee.
   * @param showIf When true, the depndent will be be shown if the dependee is true.
   * If false, the dependent will be shown only when the dependee is false.
   */
  BoolVisualDependency(std::string dependeeName, std::string dependentName, RCP<ParameterList> parentList, 
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
  BoolVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, bool showIf=true);
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  bool getDependeeState() const;
  
  //@}

private:

  /** \name Overridden from Dependency */
  //@{
  
  /** \brief . */
  void validateDep() const;
  
  //@}

};

/**
 * \brief A condition visual dependency says the following about the relationship between elements in a Parameter List:
 * Depending on whether or not the dependee(s) statisfy a particual condition, the dependent may or may not be displayed to the user in a GUI.
 *
 * Conditoin Visual Dependencies are unique in that via the Condition class, they allow for multiple dependees.
 * The dependee(s) of a ConditionVisualDependency must be expressed as a Condition and are subject to the consquential constraints. The dependent may be any type of parameter
 * or parameter list.
 */
class ConditionVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

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
  ConditionVisualDependency(RCP<Condition> condition,
  std::string dependentName, RCP<ParameterList> dependentParentList, bool showIf=true);

  /**
   * Constructs a ConditionVisualDependency.
   *
   * @param condition The condition that must be satisfied in order to display the dependent
   * parameter.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param showIf When true, the depndent will be be shown if the condition is true.
   * If false, the dependent will be shown only when the dependee is false.
   */
  ConditionVisualDependency(RCP<Condition> condition, ParameterParentMap dependents, bool showIf=true);
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  bool getDependeeState() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The Condition to determine whether or not the dependent is displayed.
   */
  RCP<Condition> condition_;
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const;
  
  //@}

};

/**
 * \brief A number visual dependency says the following about the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a certain value, the dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a NumberVisualDependency must be a number type and can't be an array. The dependent may be any type of parameter
 * or parameter list.
 */
template <class T>
class NumberVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a NumberVisualDependency.
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
  NumberVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, T (*func)(T) =0)
    :VisualDependency(dependeeName, dependeeParentList, dependentName, dependentParentList), func_(func)
  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberVisualDependency.
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
  NumberVisualDependency(std::string dependeeName, std::string dependentName, RCP<ParameterList> parentList, 
  T (*func)(T) =0)
    :VisualDependency(dependeeName, parentList, dependentName, parentList), func_(func)
  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberVisualDependency.
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
  NumberVisualDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, T (*func)(T) =0)
    :VisualDependency(dependeeName, dependeeParentList, dependents), func_(func)
  {
    validateDep();
  }
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  bool getDependeeState() const{
    return runFunction(getFirstDependeeValue<T>()) > 0 ? true : false;
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief the function used to determine the
   * visibility of the dependent.
   */
  T (*func_)(T);

  /**
   * \brief Run the function on the argument and returns the value of the fucntion. If no function is specified,
   * the argument is simple returned.
   *
   * @param argument the value to use as an argument for the function.
   * @return the result of running the function with the value. If no function is specified,
   * the argument is simple returned.
   */
  T runFunction(T argument) const{
    if(func_ !=0)
      return (*func_)(argument);
    else
      return argument;
  }  
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const{
    /*
     * This error should never get thrown, unless someone
     * is doing something wonky in a sublcass.
     */
    TEST_FOR_EXCEPTION(getDependees().size() != 1,
      InvalidDependencyException,
      "Uh oh. Looks like you tried to make a " 
      "Number Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
      "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
      "be doing some subclassing you sly-dog ;)\n\n" 
      "Error: A Number Visual Dependency must have exactly 1 dependee. " 
      "You have tried to assign it " << getDependees().size() << " dependees.\n" 
      "Dependees: " << getDependeeNamesString() << "\n" 
      "Dependents: " << getDependentNamesString());

    const ParameterEntry* dependee = getFirstDependee();
    TEST_FOR_EXCEPTION(
      !dependee->isType<int>()
      && !dependee->isType<short>()
      && !dependee->isType<double>()
      && !dependee->isType<float>(),
      InvalidDependencyException,
      "The dependee of a "
      "Number Visual Dependency must be of a supported number type!\n"
      "Problem dependee: " << getFirstDependeeName() << "\n"
      "Actual type: " << dependee->getAny().typeName() << "\n"
      "Dependees: " << getDependeeNamesString() << "\n"
      "Dependents: " << getDependentNamesString());
  }
  
  //@}
  
};

/**
 * \brief A NumberValidatorAspectDependency says the following about the relationship between two parameters in a dependent parameter list:
 * depending of the value of the dependee a particular aspect of the dependents validator will have a certain value.
 *
 * A NumberValidatorAspectDependency must have the following characteristics:
 *
 *   \li The dependee must be of a number type and can't be an array.
 *
 *   \li The dependent must have an enhancednumberparameter validator or an arrayNumberValidator.
 *
 *   \li The validator specified in the constructor must be the same validator being used on the dependent parameter.
 *
 *   \li The template type, dependee type, and dependent type must all be the same.
 *
 * This dependency can have an interesting effect on your program. The class will modifies
 * dependents validator. If that same validator is used more than once, every other 
 * parameter that used that validator will also see the change. Make sure that you are
 * aware of this when making your validators.
 */
template <class T>
class NumberValidatorAspectDependency : public Dependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief An enum specifying the aspect of the
   * validator that should be modified
   */
  enum ValidatorAspect{
    Min,
    Max,
    Step
  };
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a NumberValidatorDependency
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
  NumberValidatorAspectDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, 
  RCP<EnhancedNumberValidator<T> > validator,
  ValidatorAspect aspect, T (*func)(T) =0)
    :Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberValidatorAspectDep),
    aspect_(aspect),
    validator_(validator),
    func_(func)

  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberValidatorDependency
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
  RCP<ParameterList> parentList, 
  RCP<EnhancedNumberValidator<T> > validator,
  ValidatorAspect aspect, T (*func)(T) =0)
    :Dependency(dependeeName, parentList, dependentName, parentList, Dependency::NumberValidatorAspectDep),
    aspect_(aspect),
    validator_(validator),
    func_(func)
  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
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
  NumberValidatorAspectDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, 
  RCP<ArrayNumberValidator<T> > validator,
  ValidatorAspect aspect, T (*func)(T) =0)
    :Dependency(dependeeName, dependeeParentList, dependentName, dependentParentList, Dependency::NumberValidatorAspectDep),
    aspect_(aspect),
    validator_(validator->getPrototype()),
    func_(func)
  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
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
  RCP<ParameterList> parentList,
  RCP<ArrayNumberValidator<T> > validator,
  ValidatorAspect aspect, T (*func)(T) =0)
    :Dependency(dependeeName, parentList, dependentName, parentList, Dependency::NumberValidatorAspectDep),
    aspect_(aspect),
    validator_(validator),
    func_(func)

  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberValidatorDependency
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param validator The validator whose aspect will change. 
   * @param aspect The aspect of the validator that should change.
   * @param func A function specifying how the value of the validators
   * aspect should be calculated from the dependees value.
   */
  NumberValidatorAspectDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, 
  RCP<EnhancedNumberValidator<T> > validator,
  ValidatorAspect aspect, T (*func)(T) =0)
    :Dependency(dependeeName, dependeeParentList, dependents, Dependency::NumberValidatorAspectDep),
    aspect_(aspect),
    validator_(validator),
    func_(func)

  {
    validateDep();
  }

  /**
   * \brief Constructs a NumberValidatorDependency. Conveniece Constructor for ArrayNumberValidators
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param validator The validator whose aspect will change. 
   * @param aspect The aspect of the validator that should change.
   * @param func A function specifying how the value of the validators
   * aspect should be calculated from the dependees value.
   */
  NumberValidatorAspectDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents,
  RCP<ArrayNumberValidator<T> > validator,
  ValidatorAspect aspect, T (*func)(T) =0)
    :Dependency(dependeeName, dependeeParentList, dependents, Dependency::NumberValidatorAspectDep),
    aspect_(aspect),
    validator_(validator->getPrototype()),
    func_(func)
  {
    validateDep();
  }

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate(){
    T newAspectValue = runFunction(getFirstDependeeValue<T>());
    switch(aspect_){
      case NumberValidatorAspectDependency<T>::Min:
        validator_->setMin(newAspectValue);
        break;
      case NumberValidatorAspectDependency<T>::Max:
        validator_->setMax(newAspectValue);
        break;
      case NumberValidatorAspectDependency<T>::Step:
        validator_->setStep(newAspectValue);
        break;
    }
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The aspect of the validator to be modified.
   */
  ValidatorAspect aspect_;

  /**
   * \brief The validator to be modified.
   */
  RCP<EnhancedNumberValidator<T> > validator_;

  /**
   * \brief The function used to calculate the new value of the
   * aspect of the validator.
   */
  T (*func_)(T);
  
  /**
   * \brief Runs the dependency's function on the given argument and returns
   * the value that function returns.
   *
   * @param The value to run the function on.
   * @return The value the function returned.
   */
  T runFunction(T argument) const{
    if(func_ !=0)
      return (*func_)(argument);
    else
      return argument;
  }  
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const{
    /*
     * This error should never get thrown, unless someone
     * is doing something wonky in a sublcass.
     */
    TEST_FOR_EXCEPTION(getDependees().size() != 1,
      InvalidDependencyException,
      "Uh oh. Looks like you tried to make a "
      "Number Visual Dependency doesn't have exactly one dependee. This is kind of a problem. " 
      "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
      "be doing some subclassing you sly-dog ;)\n\n" 
      "Error: A Number Visual Dependency must have exactly 1 dependee. " 
      "You have tried to assign it "<< getDependees().size() << " dependees.\n" 
      "Dependees: " << getDependeeNamesString() << "\n" 
      "Dependents: " << getDependentNamesString());

    const ParameterEntry* dependee = getFirstDependee();
    std::string dependeeName = getFirstDependeeName();
    TEST_FOR_EXCEPTION(
      !dependee->isType<int>()
      && !dependee->isType<short>()
      && !dependee->isType<double>()
      && !dependee->isType<float>(),
      InvalidDependencyException,
      "The dependee of a "
      "Number Validator Aspect Dependency must be of a supported number type!\n"
      "Problem dependee: " << dependeeName << "\n"
      "Actual type: " << dependee->getAny().typeName() << "\n"
      "Dependents: " << getDependentNamesString());

    TEST_FOR_EXCEPTION(typeid(T) != dependee->getAny().type(),
      InvalidDependencyException,
    "The dependee type and EnhancedNumberValidator "
      "template type must all be the same for a Number Validator Aspect Dependency.\n"
      "Problem Dependee: " << dependeeName << "\n"
      "Dependee Type: " << dependee->getAny().typeName() << "\n"
      "Validator Template Type: " << typeid(T).name());

    typename ParameterParentMap::const_iterator it;
    ParameterEntry *currentDependent;
    for(it = getDependents().begin(); it != getDependents().end(); ++it){ 
      currentDependent = it->second->getEntryPtr(it->first);
      TEST_FOR_EXCEPTION(currentDependent->validator() == null,
        InvalidDependencyException,
    "The dependent of an "
        "Number Validator Aspect Dependency must have an EnhancedNumberValidator "
        "or an ArrayNumberValidator\n"
        "Problem dependent: " << it->first << "\n" 
        "Dependees: " << getDependeeNamesString() << "\n"
        "Dependents: " << getDependentNamesString());

      TEST_FOR_EXCEPTION(validator_ != currentDependent->validator(),
        InvalidDependencyException,
    "The dependent's validator and the validator specified "
        "in the constructor must be the same for a Number Validator Aspect Dependency!\n"
        "Problem dependent: " << it->first << "\n" 
        "Dependees: " << getDependeeNamesString() << "\n"
        "Dependents: " << getDependentNamesString());

      TEST_FOR_EXCEPTION(typeid(T) != currentDependent->getAny().type(),
        InvalidDependencyException,
    "The dependent type and EnhancedNumberValidator "
        "template type must all be the same for a Number Validator Aspect Dependency.\n"
        "Dependent: " << it->first << "\n" 
        "Dependent Type: " << currentDependent->getAny().typeName() << "\n"
        "Validator Template Type: " << typeid(T).name());
    }
  }
  
  //@}
  
};

/**
 * \brief A NumberArrayLengthDependency says the following about the relationship between two parameters:
 * The length of the dependent's array depends on the value of the dependee.
 *
 * A NumberArrayLengthDependency must have the following characteristics:
 *
 *   \li The dependee must be either of type int or type short.
 *
 *   \li The dependent must be an array.
 *
 */
class NumberArrayLengthDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a NumberArrayLengthDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, int (*func)(int) = 0);

  /**
   * \brief Constructs an NumberArrayLengthDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependentName The name of the dependent parameter.
   * @param parentList The ParameterList containing the dependent and the dependee.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(std::string dependeeName, std::string dependentName, 
  RCP<ParameterList> parentList, int (*func)(int) = 0);

  /**
   * \brief Constructs an NumberArrayLengthDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, int (*func)(int) = 0);
  
  //@}
  
  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The function used to calculate the new value of the
   * arrays length.
   */
  int (*func_)(int);
  
  /**
   * \brief Runs the dependency's function on the given argument and returns
   * the value that function returns.
   *
   * @param The value to run the function on.
   * @return The value the function returned.
   */
  int runFunction(int argument) const;

  /**
   * \brief Modifies the length of an array.
   *
   * @param newLength The new length the array should be.
   * @param dependentValue The index of the dependent array that is going to be changed.
   */
  template <class T>
  void modifyArrayLength(int newLength, ParameterEntry* dependentToModify);
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}
  
};

/**
 * \brief A StringValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A StringValidatorDependency must have the following characterisitics:
 * 
 *   \li The dependee must be of type string
 *
 *   \li The dependent must already have a validator assigned to it.
 *
 *   \li All of the validators that may be applied to the dependent must match the type of the
 *   validator that is currently on the dependent.
 *
 */
class StringValidatorDependency : public ValidatorDependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Conveniece typedef
   */
  typedef std::map<std::string, RCP<const ParameterEntryValidator> > ValueToValidatorMap;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a StringValidatorDependency.
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
  StringValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList,  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator);

  /**
   * \brief Constructs a StringValidatorDependency.
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
  RCP<ParameterList> parentList,  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator);

  /**
   * \brief Constructs a StringValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param valuesAndValidators A map associating string values with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied to the dependent based on the dependees value.
   * @param defaultValidator If a value is entered in the dependee that is not in the valuesAndValidators map,
   * this is the validator that will be assigned to the dependent.
   */
  StringValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents,  ValueToValidatorMap valuesAndValidators, 
  RCP<ParameterEntryValidator> defaultValidator);
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The default validator to be used if a request is made for a value that does not
   * appear in the valuesAndValidators map.
   */
  RCP<ParameterEntryValidator> defaultValidator_;

  /**
   * \brief A map associating particular dependee values with validators that
   * could be placed on the dependent.
   */
  ValueToValidatorMap valuesAndValidators_;
  
  //@}
  
  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}

};

/**
 * \brief A BoolValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A StringValidatorDependency must have the following characterisitics:
 *
 *   \li The dependee must be of type bool
 *
 *   \li The dependent must already have a validator assigned to it.
 *
 *   \li The "true" and "false" validators that may be applied to the dependent must match the type of the
 *   validator that is currently on the dependent.
 *
 */
class BoolValidatorDependency : public ValidatorDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependentName The name of the dependent parameter.
   * @param dependentParentList The ParameterList containing the dependent.
   * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
   * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
   */
  BoolValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator);

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependentName The name of the dependent parameter.
   * @param parentList The ParameterList containing the dependent and the dependee.
   * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
   * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
   */
  BoolValidatorDependency(std::string dependeeName, 
  std::string dependentName, RCP<ParameterList> parentList,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator);

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param trueValidator The validator to be used on the dependent if the dependee is set to true.
   * @param falseValidator The validator to be used on the dependent if the dependee is set to false.
   */
  BoolValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents,
  RCP<const ParameterEntryValidator> trueValidator,
  RCP<const ParameterEntryValidator> falseValidator);
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  void evaluate();
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The validators to be used when the dependee is either
   * true or false.
   */
  RCP<const ParameterEntryValidator> trueValidator_, falseValidator_;
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}
  //
};


/**
 * \brief A RangeValidatorDependency says the following about the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a particular validator from
 * a given set of validators.
 *
 * A RangeValidatorDependency achieves this by associating ranges of numerical values with validators.
 * If the dependees value falls within the one of the ranges, the validator associated with the range is
 * used on the dependent. All ranges are inclusive.
 *
 * A RangeValidatorDependency must have the following characterisitics:
 *
 *   \li The dependee must be a number type
 *
 *   \li The dependent must already have a validator assigned to it.
 *
 *   \li All of the validators that may be applied to the dependent must match the type of the
 *   validator that is currently on the dependent.
 *
 */
template<class S>
class RangeValidatorDependency : public ValidatorDependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Convenience typedef
   */
  typedef std::pair<S,S> Range;

  /**
   * \brief Convenience typedef
   */
  typedef std::map<Range, RCP<const ParameterEntryValidator> > RangeToValidatorMap;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a RangeValidatorDependency.
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
  RangeValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  std::string dependentName, RCP<ParameterList> dependentParentList, RangeToValidatorMap rangesAndValidators,
  RCP<ParameterEntryValidator> defaultValidator)
    :ValidatorDependency(dependeeName, dependeeParentList, dependentName, dependentParentList),
    defaultValidator_(defaultValidator),
    rangesAndValidators_(rangesAndValidators)
  {
    validateDep();
  }

  /**
   * \brief Constructs a RangeValidatorDependency.
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
  RCP<ParameterList> parentList, RangeToValidatorMap rangesAndValidators,
  RCP<ParameterEntryValidator> defaultValidator)
    :ValidatorDependency(dependeeName, parentList, dependentName, parentList),
    defaultValidator_(defaultValidator),
    rangesAndValidators_(rangesAndValidators)
  {
    validateDep();
  }

  /**
   * \brief Constructs a RangeValidatorDependency.
   *
   * @param dependeeName The name of the dependee parameter.
   * @param dependeeParentList The ParameterList containing the dependee.
   * @param dependents A map containing dependent Parameters associated with their paraent ParameterLists.
   * @param rangesAndValidators A map associating ranges of values with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied to the dependent based on the dependees value.
   * @param defaultValidator If a value is entered in the dependee that does not fall within any of the ranges in
   * the rangesAndValidators map, this is the validator that will be assigned to the dependent.
   */
  RangeValidatorDependency(std::string dependeeName, RCP<ParameterList> dependeeParentList,
  ParameterParentMap dependents, RangeToValidatorMap rangesAndValidators,
  RCP<ParameterEntryValidator> defaultValidator)
    :ValidatorDependency(dependeeName, dependeeParentList, dependents),
    defaultValidator_(defaultValidator),
    rangesAndValidators_(rangesAndValidators)
  {
    validateDep();
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate(){
    typename RangeToValidatorMap::const_iterator it;
    S dependeeValue = getFirstDependeeValue<S>();
    for(it = rangesAndValidators_.begin(); it != rangesAndValidators_.end(); ++it){
      S min = it->first.first;
      S max = it->first.second;
      if(dependeeValue >= min && dependeeValue <=max){
        typename ParameterParentMap::const_iterator it2;
        ParameterEntry *currentDependent;
        for(it2 = getDependents().begin(); it2 != getDependents().end(); ++it2){ 
          currentDependent = it2->second->getEntryPtr(it2->first);
          currentDependent->setValidator(it->second);
        }
        return;
      }
    }
  }
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The default validator
   */
  RCP<const ParameterEntryValidator> defaultValidator_;

  /**
   * \brief A map associating ranges with validators.
   */
  RangeToValidatorMap rangesAndValidators_;
  
  //@}

  /** \name Overridden from Dependency */
  //@{
  
  /** \brief . */
  void validateDep() const{
    /*
     * This error should never get thrown, unless someone
     * is doing something wonky in a sublcass.
     */
    TEST_FOR_EXCEPTION(getDependees().size() != 1,
      InvalidDependencyException,
    "Uh oh. Looks like you tried to make a "
      "Number Visual Dependency doesn't have exactly one dependee. This is kind of a problem. "
      "You should probably take a look into it. I'm actually amazed you even threw this error. You must "
      "be doing some subclassing you sly-dog ;)\n\n" 
      "Error: A Number Visual Dependency must have exactly 1 dependee. "
      "You have tried to assign it "<< getDependees().size() << " dependees.\n"
      "Dependees: " << getDependeeNamesString() << "\n" 
      "Dependents: " << getDependentNamesString());

    const ParameterEntry* dependee = getFirstDependee();
    std::string dependeeName = getFirstDependeeName();
    TEST_FOR_EXCEPTION(
    !dependee->isType<int>() 
    && !dependee->isType<short>() 
    && !dependee->isType<double>() 
    && !dependee->isType<float>(),
      InvalidDependencyException,
    "The dependee of a "
      "Range Validator Dependency must be of a supported number type!\n"
      "Problem dependee: " << dependeeName << "\n"
      "Actual type: " << dependee->getAny().typeName() << "\n"
      "Dependents: " << getDependentNamesString());

    typename RangeToValidatorMap::const_iterator it;
    for(it = rangesAndValidators_.begin(); it != rangesAndValidators_.end(); ++it){
      typename ParameterParentMap::const_iterator it2;
      ParameterEntry *currentDependent;
      for(it2 = getDependents().begin(); it2 != getDependents().end(); ++it2){ 
        currentDependent = it2->second->getEntryPtr(it2->first);
        TEST_FOR_EXCEPTION(typeid(*(currentDependent->validator().get())) != typeid(*(it->second.get())),
          InvalidDependencyException,
      "The validator of a dependent of a "
          "Range Validator Dependency must be the same type as all of the validators "
          "in the rangesAndValidators map.\n"
          "Note this means that the dependent must have an initial validator.\n"
          "Problem dependent: " << it2->first << "\n"
          "Validator Type: " << typeid(*(currentDependent->validator())).name() << "\n"
          "One of the validators in the rangesAndValidators map is of type: " << typeid(*(it->second)).name());
      }
    }
  }
  
  //@}

};


}
#endif //TEUCHOS_STANDARDDEPENDCIES_HPP_

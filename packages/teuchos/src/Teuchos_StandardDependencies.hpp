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
#include "Teuchos_StandardConditions.hpp"
#include "Teuchos_StandardFunctionObjects.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"


namespace Teuchos{


/**
 * \brief An abstract parent class for all visual dependencies.
 *
 * IMPORTANT NOTE:
 * If a parameter becomes hidden, it's validity will not necessarily 
 * be checked. This means that it
 * is indeed possible for a non-valid ParameterList to occur. Make 
 * sure that you program code takes
 * this into account.
 */
class VisualDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if the dependee is 
   * true.
   */
  VisualDependency(
    RCP<const ParameterEntry> dependee, 
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependent parameters.
   * @param showIf When true, the depndent will be be shown if the dependee is 
   * true.
   */
  VisualDependency(
    RCP<const ParameterEntry> dependee,
    ParameterEntryList dependents,
    bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependees The dependees.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if the dependee is 
   * true.
   */
  VisualDependency(
    ConstParameterEntryList dependees, 
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * \brief Constructs a VisualDependency.
   *
   * @param dependees The dependees.
   * @param dependents The dependets.
   * @param showIf When true, the depndent will be be shown if the dependee is 
   * true.
   */
  VisualDependency(
    ConstParameterEntryList dependees,
    ParameterEntryList dependents,
    bool showIf=true);

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

  /** \brief Get's the value of the showIf variable. */
  bool getShowIf() const;

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
   * \brief Whether or not the dependent is currently visible.
   */
  bool dependentVisible_;

  /**
   * \brief Whether or not to show the dependent if the dependee is set to the 
   * value.
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
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   */
  ValidatorDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent);

  /**
   * \brief Constructs a ValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   */
  ValidatorDependency(
    RCP<const ParameterEntry> dependee, 
    ParameterEntryList dependents);

  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  virtual void evaluate() = 0;
  
  //@}

};

/**
 * \brief A string visual depdencies says the following about the 
 * relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a particular value, 
 * the dependent may or may not be displayed to the user in a UI.
 * 
 * The dependee of a StringVisualDependency must be of type string and 
 * can't be an array. The dependent may be any type of
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
   * @param dependee The dependee paramter.
   * @parame dependent The dependent parameter.
   * @param value The value of the depndee that affects the visiblity 
   * of the dependent.
   * @param showIf When true, the depndent will be be shown 
   * if the dependee 
   * is set to the same value as specified by the value parameter.
   * If false, the dependent will be shown only when the dependee is 
   * set to a value other than the one specified by the value parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    std::string value,
    bool showIf=true);

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param values The values of the depndee that affect the 
   * visiblity of the dependent.
   * @param showIf When true, the depndent will be be shown if 
   * the dependee is set to one of the values specified by the 
   * values parameter.
   * If false, the dependent will be shown only when the dependee is set 
   * to a value other than the ones specified by the values parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    const ValueList& values,
    bool showIf=true);

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents
   * @param value The value of the depndee that affects the visiblity 
   * of the dependent.
   * @param showIf When true, the depndent will be be shown if 
   * the dependee is set to one of the values specified by the values 
   * parameter. If false, the dependent will be shown only when the 
   * dependee is set to a value other than the ones specified by the 
   * values parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents, 
    const std::string& value,
    bool showIf=true);

  /**
   * Constructs a StringVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents
   * @param values The values of the depndee that affect the 
   * visiblity of the dependent.
   * @param showIf When true, the depndent will be be 
   * shown if the dependee is 
   * set to one of the values specified by the values parameter.
   * If false, the dependent will be shown only when the dependee 
   * is set to a value other than the ones specified by 
   * the values parameter.
   */
  StringVisualDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents, 
    const ValueList& values,
    bool showIf=true);
  
  //@}

  /** \name Attribute/Query Functions */
  //@{

  /** \brief get the List of values the dependee will be checked against. */
  const ValueList& getValues() const;

  //@}

  /** \name Overridden from VisualDependency */
  //@{

  /** \brief . */
  bool getDependeeState() const;
  
  //@}
  
  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * The value used to deteremine the visiblity of the dependent.
   */
  const ValueList values_;
  
  //@}
  
};


/** \brief Speicialized class for retrieving a dummy object of type
 * StringVisualDependency.
 *
 * \relates StringVisualDependency
 */
template<>
class DummyObjectGetter<StringVisualDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * StringVisualDependency.
  */
  static RCP<StringVisualDependency> getDummyObject();
  
  //@}
  
};


/**
 * \brief A bool visual dependency says the following about the 
 * relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee is true or false, the 
 * dependent may or may not be displayed to the user in a GUI.
 *
 * The dependee of a BoolVisualDependency must be of type bool and can't 
 * be an array. The dependent may be any type of parameter
 * or parameter list.
 */
class BoolVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * Constructs a BoolVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if the dependee is 
   * true.
   * If false, the dependent will be shown only when the dependee is false.
   */
  BoolVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * Constructs a BoolVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependent parameters.
   * @param showIf When true, the depndent will be be shown if the dependee is 
   * true.
   * If false, the dependent will be shown only when the dependee is false.
   */
  BoolVisualDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents, 
    bool showIf=true);
  
  //@}

  /** \name Overridden from VisualDependency */
  //@{

  /** \brief . */
  bool getDependeeState() const;
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{
  
  /** \brief . */
  void validateDep() const;
  
  //@}

};


/** \brief Speicialized class for retrieving a dummy object of type
 * BoolVisualDependency.
 *
 * \relates BoolVisualDependency
 */
template<>
class DummyObjectGetter<BoolVisualDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * BoolVisualDependency.
  */
  static RCP<BoolVisualDependency> getDummyObject();
  
  //@}
  
};

/**
 * \brief A condition visual dependency says the following about the 
 * relationship between elements in a Parameter List:
 * Depending on whether or not the dependee(s) statisfy 
 * a particual condition, the dependent may or may not be displayed to 
 * the user in a UI.
 *
 * Condition Visual Dependencies are unique in that via 
 * the Condition class, they allow for multiple dependees.
 * The dependee(s) of a ConditionVisualDependency must be expressed as a 
 * Condition and are subject to the consquential constraints. The 
 * dependent may be any type of parameter or parameter list.
 */
class ConditionVisualDependency : public VisualDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * Constructs a ConditionVisualDependency.
   *
   *
   * @param condition The condition that must be satisfied in 
   * order to display the dependent parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown if 
   * the condition is true. If false, the dependent will be shown 
   * only when the condition is false.
   */
  ConditionVisualDependency(
    RCP<const Condition> condition,
    RCP<ParameterEntry> dependent,
    bool showIf=true);

  /**
   * Constructs a ConditionVisualDependency.
   *
   * @param condition The condition that must be satisfied in 
   * order to display the dependent parameter.
   * @param dependents The dependent parameters.
   * @param showIf When true, the depndent will be be shown if 
   * the condition is true. If false, the dependent will be shown 
   * only when the condition is false.
   */
  ConditionVisualDependency(
    RCP<const Condition> condition, 
    Dependency::ParameterEntryList dependents,
    bool showIf=true);
  
  //@}
  
  /** \name Getter Functions */
  //@{
  
  /** \brief Gets the condition being used in this dependency */
  RCP<const Condition> getCondition() const;

  /** \name Overridden from VisualDependency */
  //@{

  /** \brief . */
  bool getDependeeState() const;
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const {}
  
  //@}


private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The Condition to determine whether or not the dependent is displayed.
   */
  RCP<const Condition> condition_;
  
  //@}

};


/** \brief Speicialized class for retrieving a dummy object of type
 * ConditionVisualDependency.
 *
 * \relates ConditionVisualDependency
 */
template<>
class DummyObjectGetter<ConditionVisualDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type ConditionVisualDependency.
  */
  static RCP<ConditionVisualDependency> getDummyObject();

  //@}
  
};


/**
 * \brief A number visual dependency says the following about 
 * the relationship between two elements in a Parameter List:
 * Depending on whether or not the dependee has a certain value, 
 * the dependent may or may not be displayed to the user in a UI.
 *
 * The dependee of a NumberVisualDependency must 
 * be a number type and can't be an array. The dependent may be 
 * any type of parameter or parameter list.
 *
 * If no function is provided, then the value of the Dependee
 * is simply compared to 0. If it is greater than 0, the
 * dependency evaluates to true. Otherwise it evaluates to false.
 *
 * If a function is provided, then the value of the Dependee
 * is first ran through that function. The result of that
 * function is then compared to 0 using the same
 * criteria as above.
 */
template <class T>
class NumberVisualDependency : public VisualDependency{

public:

  /** \name Constructors*/
  //@{

  /**
   * \brief Constructs a NumberVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param showIf When true, the depndent will be be shown 
   * if the dependee 
   * is set to the same value as specified by the value parameter.
   * If false, the dependent will be shown only when the dependee is 
   * set to a value other than the one specified by the value parameter.
   * @param func A function that takes the dependees value, does some 
   * calculations on it, and then returns a value. If this value is 
   * greater than 0, the dependent is show. If the value returned is
   * less than or equal to zero, the dependent is not shown. If showIf is set
   * to valse then these visibility results will be reversed. If no 
   * fuction is specified, the direct value of the dependee will be used 
   * to determine the dependents visibility in a similar fashion (postive
   * numbers causing the dependent to be displayed and 0 or 
   * negative numbers causing the dependent to be hidden).
   */
  NumberVisualDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    bool showIf=true,
    RCP<SimpleFunctionObject<T> > func=null);

  /**
   * \brief Constructs a NumberVisualDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param showIf When true, the depndent will be be shown 
   * if the dependee 
   * is set to the same value as specified by the value parameter.
   * If false, the dependent will be shown only when the dependee is 
   * set to a value other than the one specified by the value parameter.
   * @param func A function that takes the dependees value, does some 
   * calculations on it, and then returns a value. If this value is 
   * greater than 0, the dependent is show. If the value returned is
   * less than or equal to zero, the dependent is not shown. If showIf is set
   * to false, then these visibility results will be reversed. If no 
   * fuction is specified, the direct value of the dependee will be used 
   * to determine the dependents visibility in a similar fashion (postive
   * numbers causing the dependent to be displayed and 0 or 
   * negative numbers causing the dependent to be hidden).
   */
  NumberVisualDependency(
    RCP<const ParameterEntry> dependee,
    ParameterEntryList dependents,
    bool showIf=true,
    RCP<SimpleFunctionObject<T> > func=null);

  //@}

  /** \name Overridden from VisualDependency */
  //@{

  bool getDependeeState() const;
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}
  
  /** \name Getter Functions */
  //@{
  
  /** \brief Const version of function getter. */
  RCP<const SimpleFunctionObject<T> > getFunctionObject() const;

  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void validateDep() const;
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief the function used to determine the
   * visibility of the dependent.
   */
    RCP<SimpleFunctionObject<T> > func_;
  
  //@}
  //
};

template<class T>
NumberVisualDependency<T>::NumberVisualDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  bool showIf,
  RCP<SimpleFunctionObject<T> > func)
  :VisualDependency(dependee, dependent, showIf),
  func_(func)
{
  validateDep();
}

template<class T>
NumberVisualDependency<T>::NumberVisualDependency(
  RCP<const ParameterEntry> dependee,
  ParameterEntryList dependents,
  bool showIf,
  RCP<SimpleFunctionObject<T> > func)
  :VisualDependency(dependee, dependents, showIf),
  func_(func)
{
  validateDep();
}

template<class T>
bool NumberVisualDependency<T>::getDependeeState() const{
  T value = getFirstDependeeValue<T>();
  if(!func_.is_null()){
    value = func_->runFunction(value);
  }
  return value > ScalarTraits<T>::zero() ? true : false;
}

template<class T>
std::string NumberVisualDependency<T>::getTypeAttributeValue() const{
  return "NumberVisualDependency(" + TypeNameTraits<T>::name() + ")";
}

template<class T>
RCP<const SimpleFunctionObject<T> > 
  NumberVisualDependency<T>::getFunctionObject() const
{
  return func_.getConst();
}

template<class T>
void NumberVisualDependency<T>::validateDep() const{
  RCP<const ParameterEntry> dependee = getFirstDependee();
  TEST_FOR_EXCEPTION(
    !dependee->isType<T>(),
    InvalidDependencyException,
    "The dependee of a " <<
    "Number Visual Dependency must the same type as the dependency's " <<
    "template type!" << std::endl <<
    "Type Encountered: " << dependee->getAny().typeName() << std::endl <<
    "Template Type: " << TypeNameTraits<T>::name() << std::endl <<
    std::endl);
}
  

/** \brief Speicialized class for retrieving a dummy object of type
 * NumberVisualDependency.
 *
 * \relates NumberVisualDependency
 */
template<class T>
class DummyObjectGetter<NumberVisualDependency<T> >{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * NumberVisualDependency.
  */
  static RCP<NumberVisualDependency<T> > getDummyObject();
  
  //@}
  
};

template<class T> 
RCP<NumberVisualDependency<T> >
  DummyObjectGetter<NumberVisualDependency<T> >::getDummyObject()
{
  return rcp(new NumberVisualDependency<T>(
    rcp(new ParameterEntry(ScalarTraits<T>::zero())),
    DummyObjectGetter<ParameterEntry>::getDummyObject()));
}

/**
 * \brief A NumberArrayLengthDependency says the following about the 
 * relationship between two parameters:
 * The length of the dependent's array depends on the value 
 * of the dependee.
 *
 * A NumberArrayLengthDependency must have the following characteristics:
 *
 *   \li The dependee type must be an ordinal.
 *
 *   \li The dependent must be an array.
 *   
 *   \li When supplying template parameters for this class, the dependent type
 *   should be the type which the dependent array is templated on. For 
 *   example: if the dependent is of type Array<int> then the
 *   NumberArrayLengthDependency's dependent template parameter should be set
 *   to int.
 *
 */
template<class DependeeType, class DependentType>
  class NumberArrayLengthDependency : public Dependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a NumberArrayLengthDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RCP<SimpleFunctionObject<DependeeType> > func=null);


  /**
   * \brief Constructs a NumberArrayLengthDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param func A function specifying how the arrays length 
   * should be calculated from the dependees value.
   */
  NumberArrayLengthDependency(
    RCP<const ParameterEntry> dependee,
    ParameterEntryList dependent,
    RCP<SimpleFunctionObject<DependeeType> > func=null);

  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

  /** \name Getters */
  //@{

  /** \brief Const version of function getter. */
  RCP<const SimpleFunctionObject<DependeeType> > 
    getFunctionObject() const;

  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}
  //
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief The function used to calculate the new value of the
   * arrays length.
   */
    RCP<SimpleFunctionObject<DependeeType> > func_;
  
  /**
   * \brief Modifies the length of an array.
   *
   * @param newLength The new length the array should be.
   * @param dependentValue The index of the dependent array that is going to be changed.
   */
  void modifyArrayLength(
    DependeeType newLength, RCP<ParameterEntry> dependentToModify);
  
  //@}
  
};

template<class DependeeType, class DependentType>
NumberArrayLengthDependency<DependeeType, DependentType>::NumberArrayLengthDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  RCP<SimpleFunctionObject<DependeeType> > func):
  Dependency(dependee, dependent),
  func_(func)
{
  validateDep();
}

template<class DependeeType, class DependentType>
NumberArrayLengthDependency<DependeeType, DependentType>::NumberArrayLengthDependency(
  RCP<const ParameterEntry> dependee,
  ParameterEntryList dependents,
  RCP<SimpleFunctionObject<DependeeType> > func):
  Dependency(dependee, dependents),
  func_(func)
{
  validateDep();
}


template<class DependeeType, class DependentType>
std::string 
NumberArrayLengthDependency<DependeeType, DependentType>::getTypeAttributeValue()
const
{
  return "NumberArrayLengthDependency(" +
    TypeNameTraits<DependeeType>::name() + ", " +
    TypeNameTraits<DependentType>::name() +")";
}

template<class DependeeType, class DependentType>
RCP<const SimpleFunctionObject<DependeeType> >
NumberArrayLengthDependency<DependeeType, DependentType>::getFunctionObject()
const
{
  return func_.getConst();
}

template <class DependeeType, class DependentType>
void 
NumberArrayLengthDependency<DependeeType, DependentType>::modifyArrayLength(
  DependeeType newLength, RCP<ParameterEntry> dependentToModify)
{
  const Array<DependentType> originalArray = 
    any_cast<Array<DependentType> >(dependentToModify->getAny()); 
  Array<DependentType> newArray(newLength);
  DependeeType i;
  for(
    i=OrdinalTraits<DependeeType>::zero(); 
    i<originalArray.size() && i<newLength; 
    ++i)
  {
    newArray[i] = originalArray[i];
  }

  dependentToModify->setValue(newArray,
    false, dependentToModify->docString(), dependentToModify->validator());
}

template<class DependeeType, class DependentType>
void 
NumberArrayLengthDependency<DependeeType, DependentType>::evaluate(){
  DependeeType newLength = getFirstDependeeValue<DependeeType>();
  if(!func_.is_null()){
    newLength = func_->runFunction(newLength);
  }

  TEST_FOR_EXCEPTION(newLength < OrdinalTraits<DependeeType>::zero(),
    Exceptions::InvalidParameterValue,
    "Ruh Roh Shaggy! Looks like a dependency tried to set the length "
    "of the Array(s) to a negative number. Silly. You can't have "
    "an Array with a negative length!" << std::endl << std::endl <<
    "Error:" << std::endl <<
    "An attempt was made to set the length of an Array to a negative "
    "number by a NumberArrayLengthDependency" << std::endl << std::endl);
  for(
    ParameterEntryList::iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    modifyArrayLength(newLength, *it);
  }
}

template<class DependeeType, class DependentType>
void 
NumberArrayLengthDependency<DependeeType, DependentType>::validateDep() 
  const
{
  TEST_FOR_EXCEPTION(
    typeid(DependeeType) != getFirstDependee()->getAny().type(),
    InvalidDependencyException,
    "Ay no! The dependee parameter types don't match." << std::endl <<
    "Dependee Template Type: " << TypeNameTraits<DependeeType>::name() <<
    std::endl <<
    "Dependee Parameter Type: " << getFirstDependee()->getAny().typeName()
    << std::endl << std::endl);

  for(
    ConstParameterEntryList::const_iterator it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    TEST_FOR_EXCEPTION(
      typeid(Teuchos::Array<DependentType>) != (*it)->getAny().type(),
        InvalidDependencyException,
        "Ay no! The dependent parameter types don't match." << std::endl <<
        "Dependent Template Type: " << 
        TypeNameTraits<DependentType>::name() << std::endl <<
        "Dependent Parameter Type: " << 
        (*it)->getAny().typeName() << std::endl << std::endl);
  }
}

/** \brief Speicialized class for retrieving a dummy object of type
 * NumberArrayLengthDependency.
 *
 * \relates NumberArrayLengthDependency
 */
template<class DependeeType, class DependentType>
class DummyObjectGetter<NumberArrayLengthDependency<DependeeType, DependentType> >{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * NumberArrayLengthDependency.
  */
  static RCP<NumberArrayLengthDependency<DependeeType, DependentType> >
    getDummyObject();
  
  //@}
  
};

template<class DependeeType, class DependentType>
RCP<NumberArrayLengthDependency<DependeeType, DependentType> >
  DummyObjectGetter<NumberArrayLengthDependency<DependeeType, DependentType> >::getDummyObject()
{
  return rcp(
    new NumberArrayLengthDependency<DependeeType, DependentType>(
    rcp(new ParameterEntry(ScalarTraits<DependeeType>::zero())),
    rcp(new ParameterEntry(Array<DependentType>(1)))));
}

/**
 * \brief A StringValidatorDependency says the following about 
 * the relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should 
 * use a particular validator from
 * a given set of validators.
 *
 *
 * A StringValidatorDependency must have the following characterisitics:
 * 
 *   \li The dependee must be of type string
 *
 *   \li The validators in the ValueToValidatorMap must all be the same
 *   type.
 *
 *   \li If a default validator is specified, it must be of the same
 *   type as the validators found in the ValueToValidatorMap.
 *
 * If the dependee takes on a value not in the valuesAndValidators
 * map, then the default validator is assigned to the dependent.
 */
class StringValidatorDependency : public ValidatorDependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Conveniece typedef
   */
  typedef std::map<std::string, RCP<const ParameterEntryValidator> > 
    ValueToValidatorMap;

  /**
   * \brief Conveniece typedef
   */
  typedef std::pair<std::string, RCP<const ParameterEntryValidator> > 
    ValueToValidatorPair;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a StringValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param valuesAndValidators A map associating string values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should 
   * be applied to the dependent based on the dependees value.
   * @param defaultValidator If a value is entered in the 
   * dependee that is not in the valuesAndValidators map,
   * this is the validator that will be assigned to the dependent.
   */
  StringValidatorDependency(
    RCP<const ParameterEntry> dependee, 
    RCP<ParameterEntry> dependent,
    ValueToValidatorMap valuesAndValidators, 
    RCP<ParameterEntryValidator> defaultValidator=null);

  /**
   * \brief Constructs a StringValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param valuesAndValidators A map associating string values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied to 
   * the dependent based on the dependees value.
   * @param defaultValidator If a value is entered in the dependee 
   * that is not in the valuesAndValidators map,
   * this is the validator that will be assigned to the dependent.
   */
  StringValidatorDependency(
    RCP<const ParameterEntry> dependee, 
    Dependency::ParameterEntryList dependents,
    ValueToValidatorMap valuesAndValidators, 
    RCP<ParameterEntryValidator> defaultValidator = null);

  //@}

  /** \name Getters */
  //@{

  /** \brief retrieve a const reference to the ValueToValidator map being 
   * used by this StringValidatorDependency */
  const ValueToValidatorMap& getValuesAndValidators() const;

  /** \brief . */
  RCP<const ParameterEntryValidator> getDefaultValidator() const;

  //@}
  
  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief A map associating particular dependee values with validators 
   * that could be placed on the dependent.
   */
  ValueToValidatorMap valuesAndValidators_;

  /**
   * \brief The default validator to be used if a request is made 
   * for a value that does not
   * appear in the valuesAndValidators map.
   */
  RCP<ParameterEntryValidator> defaultValidator_;
  
  //@}
  
};


/** \brief Speicialized class for retrieving a dummy object of type
 * StringValidatorDependency.
 *
 * \relates StringValidatorDependency
 */
template<>
class DummyObjectGetter<StringValidatorDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type StringValidatorDependency.
  */
  static RCP<StringValidatorDependency > getDummyObject();
  
  //@}
  
};

/**
 * \brief A BoolValidatorDependency says the following about the 
 * relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should use a 
 * particular validator from a given set of validators.
 *
 * A BoolValidatorDependency must have the following characterisitics:
 *
 *   \li The dependee must be of type bool
 *
 *   \li The false and true validators must be the same type.
 *
 */
class BoolValidatorDependency : public ValidatorDependency{

public:

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param trueValidator The validator to be used on the dependent 
   * if the dependee is set to true.
   * @param falseValidator The validator to be used on the 
   * dependent if the dependee is set to false.
   */
  BoolValidatorDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RCP<const ParameterEntryValidator> trueValidator,
    RCP<const ParameterEntryValidator> falseValidator=null);

  /**
   * \brief Constructs a BoolValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param trueValidator The validator to be used on the dependent 
   * if the dependee is set to true.
   * @param falseValidator The validator to be used on the dependent 
   * if the dependee is set to false.
   */
  BoolValidatorDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RCP<const ParameterEntryValidator> trueValidator,
    RCP<const ParameterEntryValidator> falseValidator=null);

  //@}

  /** \name Overridden from Dependency */
  //@{

  void evaluate();
  
  //@}

  /** \name Getters */
  //@{
    
  /** \brief Gets the true validator */
  RCP<const ParameterEntryValidator> getTrueValidator() const;

  /** \brief Gets the false validator */
  RCP<const ParameterEntryValidator> getFalseValidator() const;
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{

  void validateDep() const;
  
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

};

/** \brief Speicialized class for retrieving a dummy object of type
 * BoolValidatorDependency.
 *
 * \relates BoolValidatorDependency
 */
template<>
class DummyObjectGetter<BoolValidatorDependency>{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type BoolValidatorDependency.  */
  static RCP<BoolValidatorDependency > getDummyObject();
  
  //@}
  
};

/**
 * \brief A RangeValidatorDependency says the following about the
 * relationship between two parameters:
 * Dependening on the value of the dependee, the dependent should 
 * use a particular validator from a given set of validators.
 *
 * A RangeValidatorDependency achieves this by associating ranges of 
 * values with validators.
 * If the dependees value falls within the one of the ranges, 
 * the validator associated with the range is
 * used on the dependent. If the value doesn't fall within
 * any of the ranges, the dependent's validator is set to null.
 * All ranges are inclusive.
 *
 * A RangeValidatorDependency must have the following characterisitics:
 *
 *   \li The dependee type must be the same as the template type.
 *
 *   \li All the validators in the rangesAndValidators_ map must be
 *   the same type.
 */
template<class T>
class RangeValidatorDependency : public ValidatorDependency{

public:

  /** \name Public types */
  //@{

  /**
   * \brief Convenience typedef
   */
  typedef std::pair<T,T> Range;

  /**
   * \brief Convenience typedef
   */
  typedef std::map<Range, RCP<const ParameterEntryValidator> > 
    RangeToValidatorMap;

  /**
   * \brief Convenience typedef
   */
  typedef std::pair<Range, RCP<const ParameterEntryValidator> > 
    RangeValidatorPair;
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /**
   * \brief Constructs a RangeValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependent The dependent parameter.
   * @param rangesAndValidators A map associating ranges of values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied 
   * to the dependent based on the dependees value.
   */
  RangeValidatorDependency(
    RCP<const ParameterEntry> dependee,
    RCP<ParameterEntry> dependent,
    RangeToValidatorMap rangesAndValidators);

  /**
   * \brief Constructs a RangeValidatorDependency.
   *
   * @param dependee The dependee parameter.
   * @param dependents The dependents.
   * @param rangesAndValidators A map associating ranges of values 
   * with ParameterEntryValidators. This will be used
   * to deteremine what type of validator should be applied 
   * to the dependent based on the dependees value.
   */
  RangeValidatorDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RangeToValidatorMap rangesAndValidators);

  //@}

  /** \name Getters */
  //@{

  /** \brief . */
  const RangeToValidatorMap& getRangeToValidatorMap() const{
    return rangesAndValidators_;
  }
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  void evaluate();
  
  //@}

  /** \name Overridden from Dependency */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const;
  
  //@}

protected:

  /** \name Overridden from Dependency */
  //@{
  
  /** \brief . */
  void validateDep() const;
  
  //@}

  
private:

  /** \name Private Members */
  //@{
  
  /**
   * \brief A map associating ranges with validators.
   */
  RangeToValidatorMap rangesAndValidators_;
  
  void setDependentsToValidator(RCP<const ParameterEntryValidator> toSet);

  //@}

};

template<class T>
RangeValidatorDependency<T>::RangeValidatorDependency(
  RCP<const ParameterEntry> dependee,
  RCP<ParameterEntry> dependent,
  RangeToValidatorMap rangesAndValidators)
  :ValidatorDependency(dependee, dependent),
  rangesAndValidators_(rangesAndValidators)
{
  validateDep();
}

template<class T>
RangeValidatorDependency<T>::RangeValidatorDependency(
  RCP<const ParameterEntry> dependee,
  Dependency::ParameterEntryList dependents,
  RangeToValidatorMap rangesAndValidators)
  :ValidatorDependency(dependee, dependents),
  rangesAndValidators_(rangesAndValidators)
{
  validateDep();
}

template<class T>
std::string RangeValidatorDependency<T>::getTypeAttributeValue() const
{
  return "RangeValidatorDependency(" + TypeNameTraits<T>::name() + ")";
}
  

template<class T>
void RangeValidatorDependency<T>::evaluate(){
  typename RangeToValidatorMap::const_iterator it;
  T dependeeValue = getFirstDependeeValue<T>();
  for(
    it = rangesAndValidators_.begin(); 
    it != rangesAndValidators_.end(); 
    ++it)
  {
    T min = it->first.first;
    T max = it->first.second;
    if(dependeeValue >= min && dependeeValue <=max){
       setDependentsToValidator(it->second);
      return;
    }
  }
  setDependentsToValidator(null); 
}

template<class T>
void RangeValidatorDependency<T>::validateDep() const{
  RCP<const ParameterEntry> dependee = getFirstDependee();
  TEST_FOR_EXCEPTION(dependee->getAny().type() != typeid(T),
    InvalidDependencyException,
    "The dependee of a RangeValidatorDependency must be the same type as " <<
    "The RangeValidatorDependency template type!" << std::endl <<
    "Dependee Type: " << dependee->getAny().typeName() << std::endl <<
    "Templated Type: " << TypeNameTraits<T>::name() << std::endl << std::endl);
  
  TEST_FOR_EXCEPTION(
    rangesAndValidators_.size() < 1,
    InvalidDependencyException,
    "The rangesAndValidators map RangeValidatorDependency "
    "must have at least one entry!" << std::endl << std::endl);

  typename RangeToValidatorMap::const_iterator it = 
    rangesAndValidators_.begin();
  RCP<const ParameterEntryValidator> firstValidator = it->second;
  ++it;
  for(; it!=rangesAndValidators_.end(); ++it){
    TEST_FOR_EXCEPTION( typeid(*firstValidator) != typeid(*(it->second)),
      InvalidDependencyException,
      "Ay no! All of the validators in a RangeValidatorDependency "
      "must have the same type.");
    TEST_FOR_EXCEPTION(
      it->first.first > it->first.second,
      InvalidDependencyException,
      "The Range " << it->first.first << " to " << it->first.second <<
      " is invalid. The min can't be greater than the max, you silly goose!"
      );
  }

    
}

template<class T>
void RangeValidatorDependency<T>::setDependentsToValidator(
  RCP<const ParameterEntryValidator> toSet)
{
  typename ParameterEntryList::const_iterator it;
  for(
    it = getDependents().begin(); 
    it != getDependents().end(); 
    ++it)
  {
    (*it)->setValidator(toSet);
  }
}

/** \brief Speicialized class for retrieving a dummy object of type
 * RangeValidatorDependency.
 *
 * \relates RangeValidatorDependency
 */
template<class T>
class DummyObjectGetter<RangeValidatorDependency<T> >{

public:

  /** \name GetterFunctions */
  //@{

  /** \brief Retrieves a dummy object of type
  * RangeValidatorDependency.
  */
  static RCP<RangeValidatorDependency<T> > getDummyObject();
  
  //@}
  
};

template<class T>
RCP<RangeValidatorDependency<T> > 
  DummyObjectGetter<RangeValidatorDependency<T> >::getDummyObject()
{
  typename RangeValidatorDependency<T>::RangeToValidatorMap dummyMap;
  typename RangeValidatorDependency<T>::Range dummyRange(
    ScalarTraits<T>::zero(), ScalarTraits<T>::one());
  RCP<FileNameValidator> dummyValidator = 
    DummyObjectGetter<FileNameValidator>::getDummyObject();
  dummyMap.insert(typename RangeValidatorDependency<T>::RangeValidatorPair(
    dummyRange, dummyValidator));
  return rcp(new RangeValidatorDependency<T>(
    rcp(new ParameterEntry(ScalarTraits<T>::zero())),
    DummyObjectGetter<ParameterEntry>::getDummyObject(),
    dummyMap));
}


} //namespace Teuchos
#endif //TEUCHOS_STANDARDDEPENDCIES_HPP_

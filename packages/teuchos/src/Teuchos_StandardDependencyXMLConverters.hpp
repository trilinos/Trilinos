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

#ifndef TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardDependencyXMLConverters.hpp
 * \brief A collection of standard DependencyXMLConverters.
*/


#include "Teuchos_DependencyXMLConverter.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_XMLDependencyExceptions.hpp"


namespace Teuchos {


/** \brief An xml converter for VisualDepenencies
 */
class VisualDependencyConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific visual dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @return An XMLObject representing the VisualDependency.
   */
  virtual void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj) const = 0;

  /** \brief Converts any special aspects of a
   * specific visual dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependees The dependees of the visual dependency.
   * @param dependents The dependents of the visua dependency.
   * @param showIf The showIf attribute of the visual dependency.
   * @return The converted VisualDependency.
   */
  virtual RCP<VisualDependency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj) const;
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  static const std::string& getShowIfAttributeName(){
    static const std::string showIfAttributeName = "showIf";
    return showIfAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for ValidatorDependencies.
 */
class ValidatorDependencyConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific validator dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @return An XMLObject representing the ValidatorDepenency.
   */
  virtual void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj) const = 0;

  /** \brief Converts any special aspects of a
   * specific validator dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependee The dependees of the validator dependency.
   * @param dependents The dependents of the validator dependency.
   * @return The converted ValidatorDependency.
   */
  virtual RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj) const;
  
  //@}
  
};

/** \brief An xml converter for StringVisualDepenencies
 */
class StringVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents,
    bool showIf) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "stringVisual";
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  static const std::string& getStringValuesTagName(){
    static const std::string stringValuesTagName = "StringValues";
    return stringValuesTagName;
  }

  static const std::string& getStringTagName(){
    static const std::string stringTagName = "string";
    return stringTagName;
  }

  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for BoolVisualDepenencies
 */
class BoolVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents,
    bool showIf) const;

  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "boolVisual";
  }
  
  //@}

};

/** \brief An xml converter for NumberVisualDependencies
 */
template<class T>
class NumberVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents,
    bool showIf) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name() + "NumberVisual";
  }
  
  //@}

};

template<class T>
void NumberVisualDependencyConverter<T>::convertSpecialVisualAttributes(
  RCP<const VisualDependency> dependency,
  XMLObject& xmlObj) const
{
  
}
  
template<class T>
RCP<VisualDependency> 
NumberVisualDependencyConverter<T>::convertSpecialVisualAttributes(
  XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A NumberValidatorDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  return rcp(new NumberVisualDependency<T>(
    *(dependees.begin()), dependents));
}

/** \brief An xml converter for ConditionVisualDependencies
 */
class ConditionVisualDependencyConverter : public VisualDependencyConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents,
    bool showIf) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "conditionVisual";
  }
  
  //@}

};


/** \brief An xml converter for StringValidatorDependencies
 */
class StringValidatorDependencyConverter : public ValidatorDependencyConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    XMLObject& xmlObj,
    const RCP<const ParameterEntry> dependees,
    const Dependency::ParameterEntryList dependents) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "stringValidator";
  }
  
  //@}

private:
  /** \name Private Members */
  //@{
 
  /** \brief . */
  static const std::string& getValuesAndValidatorsTag(){
    static const std::string valuesAndValidatorsTag = "ValuesAndValidators";
    return valuesAndValidatorsTag;
  }

  /** \brief . */
  static const std::string& getPairTag(){
    static const std::string pairTag = "Pair";
    return pairTag;
  }

  /** \brief . */
  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  /** \brief . */
  static const std::string& getValidatorIDAttributeName(){
    static const std::string validatorIDAttributeName = "validatorID";
    return validatorIDAttributeName;
  }

  /** \brief . */
  static const std::string& getDefaultValidatorIDAttributeName(){
    static const std::string defaultValidatorIDAttributeName = "defaultValidatorID";
    return defaultValidatorIDAttributeName;
  }
  
  //@}
  
};

/** \brief An xml converter for BoolValidatorDependencies
 */
class BoolValidatorDependencyConverter : public ValidatorDependencyConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "boolValidator";
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
    
    /** \brief . */
    static const std::string& getFalseValidatorIDAttributeName(){
      static const std::string falseValidatorIDAttributeName = 
        "falseValidatorID";
      return falseValidatorIDAttributeName;
    }

    /** \brief . */
    static const std::string& getTrueValidatorIDAttributeName(){
      static const std::string trueValidatorIDAttributeName = 
        "trueValidatorID";
      return trueValidatorIDAttributeName;
    }
  
  //@}
  
};

/** \brief An xml converter for RangeValidatorDependencies
 */
template<class T>
class RangeValidatorDependencyConverter : public ValidatorDependencyConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    XMLObject& xmlObj,
    const RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents) const;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name() + "RangeVisual";
  }
  
  //@}

private:
  /** \name Private Members */
  //@{
 
  /** \brief . */
  static const std::string& getRangesAndValidatorsTag(){
    static const std::string rangesAndValidatorsTag = "RangesAndValidators";
    return rangesAndValidatorsTag;
  }

  /** \brief . */
  static const std::string& getPairTag(){
    static const std::string pairTag = "Pair";
    return pairTag;
  }

  /** \brief . */
  static const std::string& getMinAttributeName(){
    static const std::string minAttributeName = "min";
    return minAttributeName;
  }

  /** \brief . */
  static const std::string& getMaxAttributeName(){
    static const std::string maxAttributeName = "max";
    return maxAttributeName;
  }


  /** \brief . */
  static const std::string& getValidatorIDAttributeName(){
    static const std::string validatorIDAttributeName = "validatorID";
    return validatorIDAttributeName;
  }

  //@}
  
};

template<class T>
void
RangeValidatorDependencyConverter<T>::convertSpecialValidatorAttributes(
  RCP<const ValidatorDependency> dependency,
  XMLObject& xmlObj) const
{
  RCP<const RangeValidatorDependency<T> > castedDependency = 
    rcp_dynamic_cast<const RangeValidatorDependency<T> >(dependency, true);

  XMLObject rangesAndValidatorsTag(getRangesAndValidatorsTag());

  castedDependency->getRangeToValidatorMap();
  for(
    typename RangeValidatorDependency<T>::RangeToValidatorMap::const_iterator 
      it = castedDependency->getRangeToValidatorMap().begin();
    it != castedDependency->getRangeToValidatorMap().end();
    ++it)
  {
    T min = it->first.first;
    T max = it->first.second;
    ParameterEntryValidator::ValidatorID validatorID = 
      ParameterEntryValidator::getValidatorID(it->second);
    XMLObject pairTag(getPairTag());
    pairTag.addAttribute(getMinAttributeName(), min);
    pairTag.addAttribute(getMaxAttributeName(), max);
    pairTag.addAttribute(getValidatorIDAttributeName(), validatorID);
    rangesAndValidatorsTag.addChild(pairTag);
  }
  xmlObj.addChild(rangesAndValidatorsTag);
}

template<class T>
RCP<ValidatorDependency> 
RangeValidatorDependencyConverter<T>::convertSpecialValidatorAttributes(
  XMLObject& xmlObj,
  const RCP<const ParameterEntry> dependee,
  const Dependency::ParameterEntryList dependents) const
{

  int result = xmlObj.findFirstChild(getRangesAndValidatorsTag()); 
  TEST_FOR_EXCEPTION(result == -1,
    MissingRangesAndValidatorsTagException,
    "Error: All RangeValidatorDependencies must have a " << 
    getRangesAndValidatorsTag() << " tag!" << std::endl << std::endl);

  XMLObject rangesAndValidatorsTag = xmlObj.getChild(result);

  typename RangeValidatorDependency<T>::RangeToValidatorMap rangesAndValidators;
  for(int i = 0 ; i < rangesAndValidatorsTag.numChildren(); ++i){
    XMLObject child = rangesAndValidatorsTag.getChild(i);
    T min = child.getRequired<T>(getMinAttributeName());
    T max = child.getRequired<T>(getMinAttributeName());
    RCP<ParameterEntryValidator> validator = 
      ParameterEntryValidator::getValidator(
        child.getRequired<ParameterEntryValidator::ValidatorID>(
          getValidatorIDAttributeName()));
   
    rangesAndValidators.insert(
      typename RangeValidatorDependency<T>::RangeValidatorPair(
        typename RangeValidatorDependency<T>::Range(min, max), validator));
  }

  return rcp(new RangeValidatorDependency<T>(
    dependee, dependents, rangesAndValidators));
}


/** \brief An xml converter for NumberValidatorAspectDependencies.
 */
template<class T>
class NumberValidatorAspectDependencyConverter : 
  public DependencyXMLConverter{

public:

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj) const;
  
  //@}
  
  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name() + "NumberValidatorAspect";
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  static const std::string& getAspectAttributeName(){
    static const std::string aspectAttributeName = "aspect";
    return aspectAttributeName;
  }

  static const std::string& getValidatorIDAttributeName(){
    static const std::string aspectAttributeName = "validatorID";
    return aspectAttributeName;
  }

  //@}
  
};

template<class T>
RCP<Dependency> NumberValidatorAspectDependencyConverter<T>::convertXML(
  const XMLObject& xmlObj, 
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A NumberValidatorAspectDependency can only have 1 dependee!" <<
    std::endl << std::endl);

  typename NumberValidatorAspectDependency<T>::ValidatorAspect aspect = 
    NumberValidatorAspectDependency<T>::getAspectStringEnum(
      xmlObj.getAttribute(getAspectAttributeName()));

  RCP<EnhancedNumberValidator<T> > validator = 
    rcp_dynamic_cast<EnhancedNumberValidator<T> >(
      ParameterEntryValidator::getValidator(
        xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
          getValidatorIDAttributeName())), true);
  return rcp(new NumberValidatorAspectDependency<T>(
    *(dependees.begin()), dependents, aspect, validator)); 
}

template<class T>
void 
NumberValidatorAspectDependencyConverter<T>::convertDependency(
  const RCP<const Dependency> dependency, 
  XMLObject& xmlObj) const
{
  RCP<const NumberValidatorAspectDependency<T> > castedDependency = 
    rcp_dynamic_cast<const NumberValidatorAspectDependency<T> >(
      dependency, true);
  std::string aspectString = 
    NumberValidatorAspectDependency<T>::getAspectString(
      castedDependency->getAspect());
  xmlObj.addAttribute(getAspectAttributeName(), aspectString);
}

/** \brief An xml converter for NumberArrayLengthDependencies.
 */
class NumberArrayLengthDependencyConverter : public DependencyXMLConverter{

public:

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependents) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj) const;
  
  //@}
  
  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  std::string getTypeAttributeValue() const{
    return "NumberArrayLength";
  }
  
  //@}

};



} // namespace Teuchos


#endif // TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP


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

#ifndef TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardDependencyXMLConverters.hpp
 * \brief A collection of standard DependencyXMLConverters.
*/


#include "Teuchos_DependencyXMLConverter.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_XMLDependencyExceptions.hpp"
#include "Teuchos_FunctionObjectXMLConverterDB.hpp"


namespace Teuchos {


/** \brief An xml converter for VisualDepenencies
 */
class TEUCHOS_LIB_DLL_EXPORT VisualDependencyXMLConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific visual dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   * @return An XMLObject representing the VisualDependency.
   */
  virtual void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const = 0;

  /** \brief Converts any special aspects of a
   * specific visual dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependees The dependees of the visual dependency.
   * @param dependents The dependents of the visua dependency.
   * @param showIf The showIf attribute of the visual dependency.
   * @param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   * @return The converted VisualDependency.
   */
  virtual RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  //@}
  
private:

  /** \name Private Members */
  //@{
  
  /** \brief Gets the name of the showif attribute */
  static const std::string& getShowIfAttributeName(){
    static const std::string showIfAttributeName = "showIf";
    return showIfAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for ValidatorDependencies.
 */
class TEUCHOS_LIB_DLL_EXPORT ValidatorDependencyXMLConverter : public DependencyXMLConverter{

public:

  /** \name Special converter methods */
  //@{

  /** \brief Converts any special aspects of a
   * specific validator dependency to xml.
   *
   * @param dependency The dependency being converted.
   * @param xmlObj The XMLObject to which the dependency is
   * being converted.
   * @param validatorIDsMap A map containing validators and their associated
   * IDs.
   * @return An XMLObject representing the ValidatorDepenency.
   */
  virtual void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const = 0;

  /** \brief Converts any special aspects of a
   * specific validator dependency from xml.
   *
   * @param xmlObj The xml being converted.
   * @param dependee The dependees of the validator dependency.
   * @param dependents The dependents of the validator dependency.
   * @param validatorIDsMap A map containing validators and their associated
   * IDs.
   * @return The converted ValidatorDependency.
   */
  virtual RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const = 0;
  
  //@}

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  //@}
  
};

/** \brief An xml converter for StringVisualDepenencies
 *
 * The valid XML representation of a StringVisualDependency is:
 * \code
    <Dependency showIf="showIf value" type="StringVisualDependency">
      <Dependee parameterId="Id of dependee parameter"/>
      <Dependent parameterId="Id of dependent parameter"/>
      ...Any addiditional Dependents...
      <StringValues>
        <String value="First value"/>
        <String value="Second value"/>
        ...Other Values...
      </StringValues>
    </Dependency>
  \endcode
  The "showIf" XML attribute is optional and if not present will be considered
  true.
 */
class TEUCHOS_LIB_DLL_EXPORT StringVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  /** \brief Gets the StringValues Tag */
  static const std::string& getStringValuesTagName(){
    static const std::string stringValuesTagName = "StringValues";
    return stringValuesTagName;
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief Gets the String Tag */
  static const std::string& getStringTagName(){
    static const std::string stringTagName = "String";
    return stringTagName;
  }

  /** \brief Gets the Value attribute name */
  static const std::string& getValueAttributeName(){
    static const std::string valueAttributeName = "value";
    return valueAttributeName;
  }

  //@}
  
};

/** \brief An xml converter for BoolVisualDepenencies
 *
 * The valid XML representation of a BoolVisualDependency is:
 * \code
    <Dependency showIf="showIf value" type="BoolVisualDependency">
      <Dependee parameterId="Id of dependee parameter"/>
      <Dependent parameterId="Id of dependent parameter"/>
      ...Any other dependents...
    </Dependency>
  \endcode
 *  The "showIf" XML attribute is optional and if not present will be considered
 *  true.
 */
class TEUCHOS_LIB_DLL_EXPORT BoolVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  //@}

};

/** \brief An xml converter for NumberVisualDependencies
 *
 * The valid XML representation of a NumberVisualDependency is:
 * \code
   <Dependency showIf="showIf value"
    type="NumberVisualDependency(number_type_of_dependee)"
   >
    <Dependee parameterId="Id of dependee parameter"/>
    <Dependent parameterId="Id of dependent parameter"/>
    ...Any additional Dependents...
    ...Optional function tag...
  </Dependency>
  \endcode
 *  The "showIf" XML attribute is optional and if not present will be considered
 *  true.
 */
template<class T>
class NumberVisualDependencyXMLConverter : public VisualDependencyXMLConverter{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;
  
  //@}

};

template<class T>
void NumberVisualDependencyXMLConverter<T>::convertSpecialVisualAttributes(
  RCP<const VisualDependency> dependency,
  XMLObject& xmlObj,
  const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const
{
  RCP<const NumberVisualDependency<T> > castedDependency = 
    rcp_dynamic_cast<const NumberVisualDependency<T> >(dependency);
  RCP<const SimpleFunctionObject<T> > functionObject = 
    castedDependency->getFunctionObject();

  if(functionObject != null){
    XMLObject functionXML = 
      FunctionObjectXMLConverterDB::convertFunctionObject(functionObject);
    xmlObj.addChild(functionXML);
  }
    
}
  
template<class T>
RCP<VisualDependency> 
NumberVisualDependencyXMLConverter<T>::convertSpecialVisualAttributes(
  const XMLObject& xmlObj,
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  bool showIf,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A NumberVisualDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  int functionIndex = xmlObj.findFirstChild(FunctionObject::getXMLTagName());
  RCP<SimpleFunctionObject<T> > functionObject = null;
  if(functionIndex != -1){
    functionObject = rcp_dynamic_cast<SimpleFunctionObject<T> >(
      FunctionObjectXMLConverterDB::convertXML(xmlObj.getChild(functionIndex)));
  }
  return rcp(new NumberVisualDependency<T>(
    *(dependees.begin()), dependents, showIf, functionObject));
      
}

/** \brief An xml converter for ConditionVisualDependencies
 *
 * The valid XML representation of a ConditionVisualDependency is:
 * \code
    <Dependency showIf="showIf value" type="ConditionVisualDependency">
      <Dependee parameterId="Id of first dependee parameter"/>
      <Dependee parameterId="id of second dependee parameter"/>
      ...Any additional dependees...
      <Dependent parameterId="Id of dependent"/>
      ...Any additional dependents...
      ...Condition Tag and it's children...
    </Dependency>
  \endcode
 *  The "showIf" XML attribute is optional and if not present will be considered
 *  true.
 */
class TEUCHOS_LIB_DLL_EXPORT ConditionVisualDependencyXMLConverter : 
  public VisualDependencyXMLConverter
{

public:

  /** \name Overridden from VisualDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialVisualAttributes(
    RCP<const VisualDependency> dependency,
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  RCP<VisualDependency> convertSpecialVisualAttributes(
    const XMLObject& xmlObj,
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    bool showIf,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;
  
  //@}

};


/** \brief An xml converter for StringValidatorDependencies
 *
 * The valid XML representation of a StringValidatorDependency is:
 * \code
   <Dependency type="StringValidatorDependency"
    defaultValidatorId="value of default validator"
   />
     <Dependee parameterId="Id of dependee parameter"/>
     <Dependent parameterId="Id of dependent parameter"/>
     ...Additional Dependents...
     <ValuesAndValidators>
       <Pair value="Value 1" validatorId="Id of first mapped validator"/>
       <Pair value="Value 2" validatorId="Id of second mapped validator"/>
       ...Other value-to-validator mappings
     </ValuesAndValidators>
   </Dependency>
  \endcode
 * The "defaultValidatorId" XML attribute is optional.
 */
class TEUCHOS_LIB_DLL_EXPORT StringValidatorDependencyXMLConverter : 
  public ValidatorDependencyXMLConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  static const std::string& getValuesAndValidatorsTag(){
    static const std::string valuesAndValidatorsTag = "ValuesAndValidators";
    return valuesAndValidatorsTag;
  }
  
  //@}

private:
  /** \name Private Members */
  //@{
 
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
  static const std::string& getValidatorIdAttributeName(){
    static const std::string validatorIDAttributeName = "validatorId";
    return validatorIDAttributeName;
  }

  /** \brief . */
  static const std::string& getDefaultValidatorIdAttributeName(){
    static const std::string defaultValidatorIdAttributeName = 
      "defaultValidatorId";
    return defaultValidatorIdAttributeName;
  }
  
  //@}
  
};

/** \brief An xml converter for BoolValidatorDependencies
 *
 * The valid XML representation of a BoolValidatorDependency is:
 * \code
    <Dependency type="BoolValidatorDependency"
      trueValidatorId="Id of the true validator"
      falseValidatorId="Id of the false validator"
    />
      <Dependee parameterId="Id of dependee parameter"/>
      <Dependent parameterId="Id of dependent parameter"/>
      ...Any other dependent parameters...
    </Dependency>
 \endcode
 * You don't have to include both a "trueValidatorId" and "falseValidatorId" 
 * XML attribute, but you must include at least one of them.
 */
class TEUCHOS_LIB_DLL_EXPORT BoolValidatorDependencyXMLConverter : public ValidatorDependencyXMLConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const;
  
  //@}

private:

  /** \name Private Members */
  //@{
    
    /** \brief . */
    static const std::string& getFalseValidatorIdAttributeName(){
      static const std::string falseValidatorIdAttributeName = 
        "falseValidatorId";
      return falseValidatorIdAttributeName;
    }

    /** \brief . */
    static const std::string& getTrueValidatorIdAttributeName(){
      static const std::string trueValidatorIdAttributeName = 
        "trueValidatorId";
      return trueValidatorIdAttributeName;
    }
  
  //@}
  
};

/** \brief An xml converter for RangeValidatorDependencies
 *
 * The valid XML representation of a RangeValidatorDependency is:
 * \code
    <Dependency type="RangeValidatorDependency(number_type_of_dependee)"
      defaultValidatoId="id of default validator"
    />
      <Dependee parameterId="Id of dependee parameter"/>
      <Dependent parameterId="Id of dependent parameter"/>
      ...Any other dependent parameters...
      <RangesAndValidators>
        <Pair min="min value" max="max value" 
          validatorId="Id of first mapped validator"/>
        <Pair min="min value" max="max value" 
          validatorId="Id of second mapped validator"/>
        ...Other range-to-validator mappings...
      </RangesAndValidators>
      ...Optional function tag...
    </Dependency>
  \endcode
 * The "defaultValidatorId" XML attribute is optional.
 */
template<class T>
class RangeValidatorDependencyXMLConverter : 
  public ValidatorDependencyXMLConverter{

public:

  /** \name Overridden from ValidatorDependencyConverter */
  //@{
  
  /** \brief . */
  void convertSpecialValidatorAttributes(
    RCP<const ValidatorDependency> dependency,
    XMLObject& xmlObj,
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief . */
  RCP<ValidatorDependency> convertSpecialValidatorAttributes(
    const XMLObject& xmlObj,
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  static const std::string& getRangesAndValidatorsTag(){
    static const std::string rangesAndValidatorsTag = "RangesAndValidators";
    return rangesAndValidatorsTag;
  }
  
  //@}

private:
  /** \name Private Members */
  //@{
 

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
  static const std::string& getValidatorIdAttributeName(){
    static const std::string validatorIdAttributeName = "validatorId";
    return validatorIdAttributeName;
  }

  /** \brief . */
  static const std::string& getDefaultValidatorIdAttributeName(){
    static const std::string defaultValidatorIdAttributeName = 
      "defaultValidatorId";
    return defaultValidatorIdAttributeName;
  }

  //@}
  
};

template<class T>
void
RangeValidatorDependencyXMLConverter<T>::convertSpecialValidatorAttributes(
  RCP<const ValidatorDependency> dependency,
  XMLObject& xmlObj,
  ValidatortoIDMap& validatorIDsMap) const
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
    if(validatorIDsMap.find(it->second) == validatorIDsMap.end()){
      validatorIDsMap.insert(it->second);
    }
    ParameterEntryValidator::ValidatorID validatorID = 
      validatorIDsMap.find(it->second)->second;
    XMLObject pairTag(getPairTag());
    pairTag.addAttribute(getMinAttributeName(), min);
    pairTag.addAttribute(getMaxAttributeName(), max);
    pairTag.addAttribute(getValidatorIdAttributeName(), validatorID);
    rangesAndValidatorsTag.addChild(pairTag);
  }
  xmlObj.addChild(rangesAndValidatorsTag);
  RCP<const ParameterEntryValidator> defaultValidator = 
    castedDependency->getDefaultValidator();
  if(nonnull(defaultValidator)){
    if(validatorIDsMap.find(defaultValidator) == validatorIDsMap.end()){
      validatorIDsMap.insert(defaultValidator);
    }
    xmlObj.addAttribute(
      getDefaultValidatorIdAttributeName(),
      validatorIDsMap.find(defaultValidator)->second);
  }
}

template<class T>
RCP<ValidatorDependency> 
RangeValidatorDependencyXMLConverter<T>::convertSpecialValidatorAttributes(
  const XMLObject& xmlObj,
  RCP<const ParameterEntry> dependee,
  const Dependency::ParameterEntryList dependents,
  const IDtoValidatorMap& validatorIDsMap) const
{

  int result = xmlObj.findFirstChild(getRangesAndValidatorsTag()); 
  TEST_FOR_EXCEPTION(result == -1,
    MissingRangesAndValidatorsTagException,
    "Error: All RangeValidatorDependencies must have a " << 
    getRangesAndValidatorsTag() << " tag!" << std::endl << std::endl);

  XMLObject rangesAndValidatorsTag = xmlObj.getChild(result);

  typename RangeValidatorDependency<T>::RangeToValidatorMap 
    rangesAndValidators;
  for(int i = 0 ; i < rangesAndValidatorsTag.numChildren(); ++i){
    XMLObject child = rangesAndValidatorsTag.getChild(i);
    T min = child.getRequired<T>(getMinAttributeName());
    T max = child.getRequired<T>(getMaxAttributeName());
    ParameterEntryValidator::ValidatorID currentID =
      child.getRequired<ParameterEntryValidator::ValidatorID>(
          getValidatorIdAttributeName());
      
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(currentID) == validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find validator in given ValidatorIDsMap! " << std::endl <<
      std::endl);
    RCP<ParameterEntryValidator> validator = 
      validatorIDsMap.find(currentID)->second;
   
    rangesAndValidators.insert(
      typename RangeValidatorDependency<T>::RangeValidatorPair(
        typename RangeValidatorDependency<T>::Range(min, max), validator));
  }

  RCP<ParameterEntryValidator> defaultValidator = null;
  if(xmlObj.hasAttribute(getDefaultValidatorIdAttributeName())){
    ParameterEntryValidator::ValidatorID defaultValiID = 
      xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
        getDefaultValidatorIdAttributeName());
    TEST_FOR_EXCEPTION(
      validatorIDsMap.find(defaultValiID) == validatorIDsMap.end(),
      MissingValidatorException,
      "Could not find a validator (for the default validator) " <<
      "corresponding to the ID " << defaultValiID << 
      " in the given validatorIDsMap!" << std::endl << std::endl);
    defaultValidator = validatorIDsMap.find(defaultValiID)->second;
  }

  return rcp(new RangeValidatorDependency<T>(
    dependee, dependents, rangesAndValidators, defaultValidator));
}

/**
 * \brief A converter used to convert ArrayModifierDepdencies to and from
 * xml.
 */
template<class DependeeType, class DependentType>
class ArrayModifierDependencyXMLConverter : public DependencyXMLConverter{

public:

  /** \name Overridden from DependencyXMLConverter */
  //@{

  /** \brief . */
  RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  //@}

protected:

  /**
   * \brief Obtains a concrete ArrayModifierDependency given a
   * dependee, dependtns, and a funciton object.
   *
   * Because ArrayModifierDependency is an abstact class with pure virtual
   * methods we need to be able to get a concrete object to actually 
   * return. This is the reponsibility of any classes subclassing this one.
   *
   * @param dependee The dependee to be used in the construction of the
   * concrete dependency.
   * @param dependents The dependts to be used in the construction of the
   * concrete dependency.
   * @param function The function object to be used in the construction of the
   * concrete dependency.
   * @return A concrete dependency object.which subclasses 
   * ArrayModifierDependency.
   */
  virtual RCP<ArrayModifierDependency<DependeeType, DependentType> > 
  getConcreteDependency(
    RCP<const ParameterEntry> dependee,
    const Dependency::ParameterEntryList dependents,
    RCP<const SimpleFunctionObject<DependeeType> > function) const = 0;

};

template<class DependeeType, class DependentType>
RCP<Dependency> 
ArrayModifierDependencyXMLConverter<DependeeType, DependentType>::convertXML(
  const XMLObject& xmlObj, 
  const Dependency::ConstParameterEntryList dependees,
  const Dependency::ParameterEntryList dependents,
  const XMLParameterListReader::EntryIDsMap& entryIDsMap,
  const IDtoValidatorMap& validatorIDsMap) const
{
  TEST_FOR_EXCEPTION(dependees.size() > 1,
    TooManyDependeesException,
    "A ArrayModifierDependency can only have 1 dependee!" <<
    std::endl << std::endl);
  RCP<SimpleFunctionObject<DependeeType> > functionObject = null;
  int functionIndex = xmlObj.findFirstChild(FunctionObject::getXMLTagName());
  if(functionIndex != -1){
    functionObject = rcp_dynamic_cast<SimpleFunctionObject<DependeeType> >(
      FunctionObjectXMLConverterDB::convertXML(xmlObj.getChild(functionIndex)));
  }
  return 
    getConcreteDependency(*(dependees.begin()), dependents, functionObject);
}

template<class DependeeType, class DependentType>
void
ArrayModifierDependencyXMLConverter<DependeeType, DependentType>::convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const ArrayModifierDependency<DependeeType, DependentType> > castedDep =
    rcp_dynamic_cast<const ArrayModifierDependency<DependeeType, DependentType> >(
      dependency);
  RCP<const SimpleFunctionObject<DependeeType> > functionObject = 
    castedDep->getFunctionObject();
  if(functionObject != null){
    XMLObject functionXML = FunctionObjectXMLConverterDB::convertFunctionObject(
      functionObject);
    xmlObj.addChild(functionXML);
  }
}


/** \brief An xml converter for NumberArrayLengthDependencies.
 *
 * The valid XML representation of a NumberArrayLengthDependency is:
 * \code
    <Dependency 
      type="NumberArrayLengthDependency(dependee_number_type, type_of_array_values)"
    >
      <Dependee parameterId="Id of dependee parameter"/>
      <Dependent parameterId="Id of dependent parameter"/>
      ...Any other dependent parameters...
      ...Optional Function tag...
    </Dependency>
 \endcode
 */
template<class DependeeType, class DependentType>
class NumberArrayLengthDependencyXMLConverter : 
  public ArrayModifierDependencyXMLConverter<DependeeType, DependentType>{

protected:

  /** \name Overridden from ArrayModifierDependency */
  //@{
  virtual RCP<ArrayModifierDependency<DependeeType, DependentType> > 
  getConcreteDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RCP<const SimpleFunctionObject<DependeeType> > function) const;
  //@}

};

template<class DependeeType, class DependentType>
RCP<ArrayModifierDependency<DependeeType, DependentType> > 
NumberArrayLengthDependencyXMLConverter<DependeeType, DependentType>::getConcreteDependency(
  RCP<const ParameterEntry> dependee,
  Dependency::ParameterEntryList dependents,
  RCP<const SimpleFunctionObject<DependeeType> > function) const
{
return rcp(
    new NumberArrayLengthDependency<DependeeType, DependentType>(
      dependee, dependents, function));
}


template<class DependeeType, class DependentType>
class TwoDRowDependencyXMLConverter :
  public ArrayModifierDependencyXMLConverter<DependeeType, DependentType>
{

protected:

  /** \name Overridden from ArrayModifierDependency */
  //@{
  virtual RCP<ArrayModifierDependency<DependeeType, DependentType> > 
  getConcreteDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RCP<const SimpleFunctionObject<DependeeType> > function) const;
  //@}

};

template<class DependeeType, class DependentType>
RCP<ArrayModifierDependency<DependeeType, DependentType> > 
TwoDRowDependencyXMLConverter<DependeeType, DependentType>::getConcreteDependency(
  RCP<const ParameterEntry> dependee,
  Dependency::ParameterEntryList dependents,
  RCP<const SimpleFunctionObject<DependeeType> > function) const
{
return rcp(
    new TwoDRowDependency<DependeeType, DependentType>(
      dependee, dependents, function));
}

template<class DependeeType, class DependentType>
class TwoDColDependencyXMLConverter :
  public ArrayModifierDependencyXMLConverter<DependeeType, DependentType>
{

protected:

  /** \name Overridden from ArrayModifierDependency */
  //@{
  virtual RCP<ArrayModifierDependency<DependeeType, DependentType> > 
  getConcreteDependency(
    RCP<const ParameterEntry> dependee,
    Dependency::ParameterEntryList dependents,
    RCP<const SimpleFunctionObject<DependeeType> > function) const;
  //@}

};

template<class DependeeType, class DependentType>
RCP<ArrayModifierDependency<DependeeType, DependentType> > 
TwoDColDependencyXMLConverter<DependeeType, DependentType>::getConcreteDependency(
  RCP<const ParameterEntry> dependee,
  Dependency::ParameterEntryList dependents,
  RCP<const SimpleFunctionObject<DependeeType> > function) const
{
return rcp(
    new TwoDColDependency<DependeeType, DependentType>(
      dependee, dependents, function));
}



} // namespace Teuchos


#endif // TEUCHOS_STANDARDDEPENDENCYXMLCONVERTERS_HPP

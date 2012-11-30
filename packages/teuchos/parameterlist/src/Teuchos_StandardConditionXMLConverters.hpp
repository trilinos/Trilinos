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

#ifndef TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP

/*! \file Teuchos_StandardConditionXMLConverters.hpp
 * \brief A collection of standard ConditionXMLConverters.
*/

#include "Teuchos_FunctionObjectXMLConverterDB.hpp"
#include "Teuchos_ConditionXMLConverter.hpp"
#include "Teuchos_StandardConditions.hpp"


namespace Teuchos {


/** \brief An xml converter for BoolLogicConditions
 */
class BoolLogicConditionConverter : public ConditionXMLConverter{

public:

  /** \name Special Converter Functions */
  //@{

  /** 
   * \brief Gets the specific BoolLogicCondition to be returned 
   * by this conveter when converting from XML.
   *
   * @param conditions The condition list for the BoolLogic converter 
   * that is being converted.
   */
  virtual RCP<BoolLogicCondition> getSpecificBoolLogicCondition(
    Condition::ConstConditionList& conditions) const = 0;
 
  //@}

  /** \name Overridden from ConditionXMLConverter */
  //@{

  /** \brief . */
  virtual RCP<Condition> convertXML(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  void convertCondition(
    const RCP<const Condition> condition, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;
  
  //@}

};

/** \brief An xml converter for OrConditions
 *
 * The valid XML represntation for an OrCondition is:
 * \code
    <Condition type="OrCondition">
      ...Conditions to "or" together...
    </Condition>
 \endcode
 */
class OrConditionConverter : public BoolLogicConditionConverter{

public:

  /** \name Overridden from BoolLogicConditionConverter */
  //@{

  /** \brief . */
  RCP<BoolLogicCondition> getSpecificBoolLogicCondition(
    Condition::ConstConditionList& conditions) const;
 
  //@}

};

/** \brief An xml converter for AndConditions
 *
 * The valid XML represntation for an AndCondition is:
 * \code
    <Condition type="AndCondition">
      ...Conditions to "and" together...
    </Condition>
 \endcode
 */
class AndConditionConverter : public BoolLogicConditionConverter{

public:

  /** \name Overridden from BoolLogicConditionConverter */
  //@{

  /** \brief . */
  RCP<BoolLogicCondition> getSpecificBoolLogicCondition(
    Condition::ConstConditionList& conditions) const;
 
  //@}

};


/** \brief An xml converter for EqualsConditions
 *
 * The valid XML represntation for an EqualsCondition is:
 * \code
    <Condition type="EqualsCondition">
      ...Conditions to "equals" together...
    </Condition>
 \endcode
 */
class EqualsConditionConverter : public BoolLogicConditionConverter{

public:

  /** \name Overridden from BoolLogicConditionConverter */
  //@{

  /** \brief . */
  RCP<BoolLogicCondition> getSpecificBoolLogicCondition(
    Condition::ConstConditionList& conditions) const;
 
  //@}

};

/** \brief An xml converter for NotConditions
 *
 * The valid XML represntation for an NotCondition is:
 * \code
    <Condition type="NotCondition">
      ...Condition to negate...
    </Condition>
 \endcode
 */
class NotConditionConverter : public ConditionXMLConverter{

public:

  /** \name Overridden from ConditionXMLConverter */
  //@{

  /** \brief . */
  virtual RCP<Condition> convertXML(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  void convertCondition(
    const RCP<const Condition> condition, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;

  //@}

};

/** \brief An xml converter for ParameterConditions
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterConditionConverter : public ConditionXMLConverter{

public:

  /** \name Special Converter Functions */
  //@{

  /** \brief Gets the specific ParameterCondition to be returned 
   * by this conveter
   * when converting from XML.
   *
   * @param xmlObj The xml object from which the ParameterCondition 
   * is being derived.
   * @param parameterEntry The ParameterEntry the ParameterCondition evaluates.
   */
  virtual RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry) const = 0;
  
  /** \brief Adds specific xml traits to the xmlObj for a particular
   * ParmaterCondtion
   *
   * @param condition The ParameterCondition to be converted.
   * @param xmlObj The XMLObject to which the specific traits should be 
   * added.
   */
  virtual void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const = 0;
 
  //@}

  /** \name Overridden from ConditionXMLConverter */
  //@{

  /** \brief . */
  virtual RCP<Condition> convertXML(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap) const;

  /** \brief . */
  void convertCondition(
    const RCP<const Condition> condition, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap) const;
  
  //@}
 
private:

  /** \name Private Members */
  //@{
 
  /** \brief gets the ParameterEntryID attribute name. */
  static const std::string& getParameterEntryIdAttributeName(){
    static const std::string parameterEntryIdAttributeName = "parameterId";
    return parameterEntryIdAttributeName;
  }
  
  //@}
  
};

/** \brief An xml converter for StringConditions
 * The valid XML represntation for a StringCondition is:
 * \code
    <Condition type="StringCondition" parameterId="Id of parameter">
      <Values>
        <String value="First string value"/>
        ...Additional String Values...
      </Values>
    </Condition>
 \endcode
 */
class StringConditionConverter : public ParameterConditionConverter{

public:

  /** \name Overridden from ParameterConditionConverter */
  //@{

  /** \brief . */
  RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry) const;
 
  /** \brief . */
  void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const;
 
  //@}

private:

  /** \brief . */
  static const std::string& getValuesTagName(){
    static const std::string valuesTagName = "Values";
    return valuesTagName;
  }

  /** \brief . */
  static const std::string& getStringTagName(){
    static const std::string stringTagName = "String";
    return stringTagName;
  }

  /** \brief . */
  static const std::string& getStringValueAttributeName(){
    static const std::string stringValueAttributeName = "value";
    return stringValueAttributeName;
  }


};

/** \brief An xml converter for BoolConditions
 *
 * The valid XML represntation for a BoolCondition is:
 * \code
    <Condition type="BoolCondition" parameterId="Id of parameter"/>
   \endcode
 */
class BoolConditionConverter : public ParameterConditionConverter{

public:

  /** \name Overridden from ParameterConditionConverter */
  //@{

  /** \brief . */
  RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry) const;
 
  /** \brief . */
  void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const;

  //@}

};

/** \brief An xml converter for NumberConditions
 * The valid XML represntation for a NumberCondition is:
 * \code
   <Condition type="NumberCondition(number_type_of_parameter)"
      parameterId="id of parameter">
     ...Optional Function Tag...
   </Condition>
 \endcode
 */
template<class T>
class NumberConditionConverter : public ParameterConditionConverter{

public:

  /** \name Overridden from ParameterConditionConverter */
  //@{

  /** \brief . */
  RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry) const;

  /** \brief . */
  void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const;
 
  //@}

};

template<class T>
RCP<ParameterCondition> 
NumberConditionConverter<T>::getSpecificParameterCondition(
  const XMLObject& xmlObj,
  RCP<ParameterEntry> parameterEntry) const
{
  int functionTag = xmlObj.findFirstChild(FunctionObject::getXMLTagName());
  if(functionTag == -1){
    return rcp(new NumberCondition<T>(parameterEntry));
  }
  else{
    RCP<FunctionObject> functionObj = 
      FunctionObjectXMLConverterDB::convertXML(xmlObj.getChild(functionTag));
    RCP<SimpleFunctionObject<T> > castedFunction = 
      rcp_dynamic_cast<SimpleFunctionObject<T> >(functionObj);
    return rcp(new NumberCondition<T>(parameterEntry, castedFunction));
  }
}

template<class T>
void NumberConditionConverter<T>::addSpecificXMLTraits(
  RCP<const ParameterCondition> condition, XMLObject& xmlObj) const
{
  RCP<const NumberCondition<T> > castedCondition = 
    rcp_dynamic_cast<const NumberCondition<T> >(condition);
  RCP<const SimpleFunctionObject<T> > functionObject =
    castedCondition->getFunctionObject();
  if(!functionObject.is_null()){
    XMLObject functionXML = 
      FunctionObjectXMLConverterDB::convertFunctionObject(functionObject);
    xmlObj.addChild(functionXML);
  }
}



} // namespace Teuchos


#endif // TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP


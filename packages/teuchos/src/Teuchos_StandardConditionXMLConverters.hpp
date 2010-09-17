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

#ifndef TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP

/*! \file Teuchos_StandardConditionXMLConverters.hpp
 * \brief A collection of standard ConditionXMLConverters.
*/


#include "Teuchos_ConditionXMLConverter.hpp"
#include "Teuchos_StandardConditions.hpp"


namespace Teuchos {


/** \brief An xml converter for BinaryLogicalConditions
 */
class BinaryLogicalConditionConverter : public ConditionXMLConverter{

public:

  /** \name Special Converter Functions */
  //@{

  /** \brief 
   * Gets the specific BinaryLogicalCondition to be returned 
   * by this conveter when converting from XML.
   *
   * @param conditions The condition list for the BinaryLogical converter 
   * that is being converted.
   */
  virtual RCP<BinaryLogicalCondition> getSpecificBinaryLogicalCondition(
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
 */
class OrConditionConverter : public BinaryLogicalConditionConverter{

public:

  /** \name Overridden from BinaryLogicalConditionConverter */
  //@{

  /** \brief . */
  RCP<BinaryLogicalCondition> getSpecificBinaryLogicalCondition(
    Condition::ConstConditionList& conditions) const;
 
  //@}

};

/** \brief An xml converter for AndConditions
 */
class AndConditionConverter : public BinaryLogicalConditionConverter{

public:

  /** \name Overridden from BinaryLogicalConditionConverter */
  //@{

  /** \brief . */
  RCP<BinaryLogicalCondition> getSpecificBinaryLogicalCondition(
    Condition::ConstConditionList& conditions) const;
 
  //@}

};


/** \brief An xml converter for EqualsConditions
 */
class EqualsConditionConverter : public BinaryLogicalConditionConverter{

public:

  /** \name Overridden from BinaryLogicalConditionConverter */
  //@{

  /** \brief . */
  RCP<BinaryLogicalCondition> getSpecificBinaryLogicalCondition(
    Condition::ConstConditionList& conditions) const;
 
  //@}

};

/** \brief An xml converter for NotConditions
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
class ParameterConditionConverter : public ConditionXMLConverter{

public:

  /** \name Special Converter Functions */
  //@{

  /** \brief Gets the specific ParameterCondition to be returned 
   * by this conveter
   * when converting from XML.
   *
   * @param xmlObj The xml object from which the ParameterCondition 
   * is being derived.
   * @param whenParamEqualsValue The whenParamEqualsValue for the 
   * ParameterCondition to be created
   */
  virtual RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry,
    bool whenParamEqualsValue) const = 0;
  
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
 
  static const std::string& getParameterEntryIDAttributeName(){
    static const std::string parameterEntryIDAttributeName = "parameterID";
    return parameterEntryIDAttributeName;
  }

  static const std::string& getWhenParamEqualsValueAttributeName(){
    static const std::string whenParamEqualsValueAttributeName = 
      "whenParamEqualsValue";
    return whenParamEqualsValueAttributeName;
  }

  
  //@}
  
};

/** \brief An xml converter for StringConditions
 */
class StringConditionConverter : public ParameterConditionConverter{

public:

  /** \name Overridden from ParameterConditionConverter */
  //@{

  /** \brief . */
  RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry,
    bool whenParamEqualsValue) const;
 
  /** \brief . */
  void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const;
 
  //@}

private:
  static const std::string& getValuesTagName(){
    static const std::string valuesTagName = "Values";
    return valuesTagName;
  }

  static const std::string& getStringTagName(){
    static const std::string stringTagName = "String";
    return stringTagName;
  }

  static const std::string& getStringValueAttributeName(){
    static const std::string stringValueAttributeName = "value";
    return stringValueAttributeName;
  }


};

/** \brief An xml converter for BoolConditions
 */
class BoolConditionConverter : public ParameterConditionConverter{

public:

  /** \name Overridden from ParameterConditionConverter */
  //@{

  /** \brief . */
  RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry,
    bool whenParamEqualsValue) const;
 
  /** \brief . */
  void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const;
  //@}

  /** \name Overridden from ConditionXMLConverter */
  //@{

};

/** \brief An xml converter for NumberConditions
 */
template<class T>
class NumberConditionConverter : public ParameterConditionConverter{

public:

  /** \name Overridden from ParameterConditionConverter */
  //@{

  /** \brief . */
  RCP<ParameterCondition> getSpecificParameterCondition(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> parameterEntry,
    bool whenParamEqualsValue) const
  {
    return rcp(new NumberCondition<T>(
      parameterEntry, null, whenParamEqualsValue));
  }

  /** \brief . */
  void addSpecificXMLTraits(
    RCP<const ParameterCondition> condition, XMLObject& xmlObj) const
  {}
 
  //@}

};




} // namespace Teuchos


#endif // TEUCHOS_STANDARDCONDITIONXMLCONVERTERS_HPP


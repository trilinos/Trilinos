// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardParameterEntryXMLConverters.hpp
 * \brief A collection of standard ParameterEntryXMLConverters.
*/


#include "Teuchos_ParameterEntryXMLConverter.hpp"


namespace Teuchos {


/** \brief A last resort converter for when no others will do.
 *
 * Writes out a raw string representation to xml and sets ParameterEntryValues
 * as strings when they are read back in.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT AnyParameterEntryConverter : public ParameterEntryXMLConverter{

public:

  /** \name Overridden from ParameterEntryXMLConverter */
  //@{

  /** \brief . */
  const std::string getTypeAttributeValue() const;

  /** \brief . */
  const std::string getValueAttributeValue(
    RCP<const ParameterEntry> entry) const;

  /** \brief . */
  any getAny(const XMLObject& xmlObj) const;
  //@}

};


/** \brief A standard ParameterEntryXMLConverter for most data types.
 *
 * This converter is appropriate for most data types.
 */
template<class T>
class StandardTemplatedParameterConverter :
  public ParameterEntryXMLConverter
{

public:

  /** \name Overridden from ParameterEntryXMLConverter */
  //@{

  /** \brief . */
  virtual const std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name();
  }

  /** \brief . */
  virtual const std::string getValueAttributeValue(
    RCP<const ParameterEntry> entry) const {
    return toString(any_cast<T>(entry->getAny(false)));
  }

  /** \brief . */
  any getAny(const XMLObject& xmlObj) const{
    return any(xmlObj.getRequired<T>(getValueAttributeName()));
  }

  //@}

};



} // namespace Teuchos


#endif // TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

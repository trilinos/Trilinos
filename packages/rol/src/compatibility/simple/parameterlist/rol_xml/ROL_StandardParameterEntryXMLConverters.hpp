// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
#define ROL_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

/*! \file ROL_StandardParameterEntryXMLConverters.hpp
 * \brief A collection of standard ParameterEntryXMLConverters.
*/

#include "ROL_ParameterEntryXMLConverter.hpp"
#include "ROL_toString.hpp"

namespace ROL {


/** \brief A last resort converter for when no others will do.
 *
 * Writes out a raw string representation to xml and sets ParameterEntryValues
 * as strings when they are read back in.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT AnyParameterEntryConverter : public ParameterEntryXMLConverter{

public:

  /** \name Overridden from ParameterEntryXMLConverter */
  //@{

    /* Destructor */
  virtual ~AnyParameterEntryConverter() {}

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

  /* Destructor */
  virtual ~StandardTemplatedParameterConverter() {}

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



} // namespace ROL


#endif // ROL_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

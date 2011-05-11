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
class TEUCHOS_LIB_DLL_EXPORT AnyParameterEntryConverter : public ParameterEntryXMLConverter{

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
class TEUCHOS_LIB_DLL_EXPORT StandardTemplatedParameterConverter : 
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

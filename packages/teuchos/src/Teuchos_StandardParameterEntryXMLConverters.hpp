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
*/

#include "Teuchos_ParameterEntryXMLConverter.hpp"

namespace Teuchos {

/**
 * \brief A last resort converter for when no others will do. Writes out
 * a raw string representation to xml and sets ParameterEntryValues as
 * strings when they are read back in.
 */
class AnyParameterEntryConverter : public ParameterEntryXMLConverter{

public:

  const std::string getTypeAttributeValue() const;

  const std::string getValueAttributeValue(const ParameterEntry &entry) const;

  void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;

  bool isAppropriateConverter(const ParameterEntry& entry) const;
};

/**
 * \brief A standard ParameterEntryXMLConverter for most data types.
 *
 * This converter is appropriate for most data types.
 */
template<class T>
class StandardTemplatedParameterConverter : public ParameterEntryXMLConverter{

public:

  /** \brief */
  virtual const std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name();
  }

  /** \brief */
  virtual const std::string getValueAttributeValue(const ParameterEntry& entry) const{
    return toString(any_cast<T>(entry.getAny(false)));
  }

  /** \brief */
  virtual void setEntryValue(ParameterEntry& entry, const XMLObject& xmlObj, bool isDefault) const{
    entry.setValue<T>(xmlObj.getRequired<T>(getValueAttributeName()), isDefault);
  }

  /** \brief */
  virtual bool isAppropriateConverter(const ParameterEntry& entry) const{
    return entry.isType<T>();  
  }

};

/**
 * \brief A standard ParameterEntryXMLConverter for array data types.
 *
 * This converter is appropriate for most array data types.
 */
template<class T>
class ArrayTemplatedParameterConverter : public ParameterEntryXMLConverter{

public:

  /** \brief */
  virtual const std::string getTypeAttributeValue() const{
    return TypeNameTraits<Array<T> >::name();
  }

  /** \brief */
  virtual const std::string getValueAttributeValue(const ParameterEntry& entry) const{
    return toString(any_cast<Array<T> >(entry.getAny(false)));
  }

  /** \brief */
  virtual void setEntryValue(ParameterEntry& entry, const XMLObject& xmlObj, bool isDefault) const{
    std::string arrayString = xmlObj.getRequired(getValueAttributeName());
    Array<T> convertedArray;
    convertedArray = fromStringToArray<T>(arrayString);
    entry.setValue<Array<T> >(convertedArray, isDefault);
  }

  /** \brief */
  virtual bool isAppropriateConverter(const ParameterEntry& entry) const{
    return entry.isType<Array<T> >();  
  }

};

}

#endif // TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

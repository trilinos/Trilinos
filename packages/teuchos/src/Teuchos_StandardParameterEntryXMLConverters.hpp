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

#include "Teuchos_ParameterEntryXMLConverter.hpp"

#ifndef TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardParameterEntryXMLConverters.hpp
*/
namespace Teuchos {

class AnyParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};


class IntParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ShortParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class LongParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class FloatParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class DoubleParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class StringParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class CharParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class BoolParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayIntParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayShortParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayLongParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayFloatParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayDoubleParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayStringParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

class ArrayCharParameterEntryConverter : public ParameterEntryXMLConverter{
public:
	const std::string getTypeAttributeValue() const;
	const std::string getValueAttributeValue(const ParameterEntry &entry) const;
	void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const;
	bool isAppropriateConverter(const ParameterEntry& entry) const;
};

}
#endif // TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

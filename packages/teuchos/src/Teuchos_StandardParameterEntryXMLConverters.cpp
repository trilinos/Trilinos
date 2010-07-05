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
#include "Teuchos_StandardParameterEntryXMLConverters.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos{
	std::string DefaultParameterEntryConverter::getTypeAttributeValue() const{
		return "any";
	}

	std::string DefaultParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
		return toString(entry.getAny(false));
	}

	void DefaultParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<std::string>(xmlObj.getRequired("value"), isDefault);
	}

	bool DefaultParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return false;
	}


	std::string IntParameterEntryConverter::getTypeAttributeValue() const{
		return "int";
	}

	std::string IntParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<int>(entry.getAny(false)));
	}

	void IntParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<int>(xmlObj.getRequiredInt("value"), isDefault);
	}

	bool IntParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<int>();
	}

	std::string ShortParameterEntryConverter::getTypeAttributeValue() const{
		return "short";
	}

	std::string ShortParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<short>(entry.getAny(false)));
	}

	void ShortParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<short>(xmlObj.getRequiredShort("value"), isDefault);
	}

	bool ShortParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<short>();
	}

	std::string DoubleParameterEntryConverter::getTypeAttributeValue() const{
		return "double";
	}

	std::string DoubleParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<double>(entry.getAny(false)));
	}

	void DoubleParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<double>(xmlObj.getRequiredDouble("value"), isDefault);
	}

	bool DoubleParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<double>();
	}

	std::string FloatParameterEntryConverter::getTypeAttributeValue() const{
		return "float";
	}

	std::string FloatParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<float>(entry.getAny(false)));
	}

	void FloatParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<float>(xmlObj.getRequiredFloat("value"), isDefault);
	}

	bool FloatParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<float>();
	}

	std::string StringParameterEntryConverter::getTypeAttributeValue() const{
		return "string";
	}

	std::string StringParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<std::string>(entry.getAny(false)));
	}

	void StringParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<string>(xmlObj.getRequired("value"), isDefault);
	}

	bool StringParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<std::string>();
	}

	std::string CharParameterEntryConverter::getTypeAttributeValue() const{
		return "char";
	}

	std::string CharParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<char>(entry.getAny(false)));
	}

	void CharParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		char *value = (char*)(xmlObj.getRequired("value").c_str());
		entry.setValue<char>(value[0], isDefault);
	}

	bool CharParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<char>();
	}

	std::string BoolParameterEntryConverter::getTypeAttributeValue() const{
		return "bool";
	}

	std::string BoolParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<bool>(entry.getAny(false)));
	}

	void BoolParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<bool>(xmlObj.getRequiredBool("value"), isDefault);
	}

	bool BoolParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<bool>();
	}

	std::string ArrayIntParameterEntryConverter::getTypeAttributeValue() const{
		return "Array int";
	}

	std::string ArrayIntParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<int> >(entry.getAny(false)).toString();
	}

	void ArrayIntParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<int> >(fromStringToArray<int>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayIntParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<int> >();
	}

	std::string ArrayShortParameterEntryConverter::getTypeAttributeValue() const{
		return "Array short";
	}

	std::string ArrayShortParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<short> >(entry.getAny(false)).toString();
	}

	void ArrayShortParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<short> >(fromStringToArray<short>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayShortParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<short> >();
	}

	std::string ArrayDoubleParameterEntryConverter::getTypeAttributeValue() const{
		return "Array double";
	}

	std::string ArrayDoubleParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<double> >(entry.getAny(false)).toString();
	}

	void ArrayDoubleParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<double> >(fromStringToArray<double>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayDoubleParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<double> >();
	}

	std::string ArrayFloatParameterEntryConverter::getTypeAttributeValue() const{
		return "Array float";
	}

	std::string ArrayFloatParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<float> >(entry.getAny(false)).toString();
	}

	void ArrayFloatParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<float> >(fromStringToArray<float>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayFloatParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<float> >();
	}

	std::string ArrayStringParameterEntryConverter::getTypeAttributeValue() const{
		return "Array string";
	}

	std::string ArrayStringParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<std::string> >(entry.getAny(false)).toString();
	}

	void ArrayStringParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<std::string> >(fromStringToArray<std::string>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayStringParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<std::string> >();
	}

}


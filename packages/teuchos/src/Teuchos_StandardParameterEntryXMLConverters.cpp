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
	const std::string AnyParameterEntryConverter::getTypeAttributeValue() const{
		return "any";
	}

	const std::string AnyParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
		return toString(entry.getAny(false));
	}

	void AnyParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<std::string>(xmlObj.getRequired("value"), isDefault);
	}

	bool AnyParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return true;
	}

	const std::string IntParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<int>::name();
	}

	const std::string IntParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<int>(entry.getAny(false)));
	}

	void IntParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<int>(xmlObj.getRequiredInt("value"), isDefault);
	}

	bool IntParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<int>();
	}

	const std::string ShortParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<short>::name();
	}

	const std::string ShortParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<short>(entry.getAny(false)));
	}

	void ShortParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<short>(xmlObj.getRequired<short>("value"), isDefault);
	}

	bool ShortParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<short>();
	}

	const std::string LongParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<long>::name();
	}

	const std::string LongParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<long>(entry.getAny(false)));
	}

	void LongParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<long>(xmlObj.getRequired<long>("value"), isDefault);
	}

	bool LongParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<long>();
	}

	const std::string DoubleParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<double>::name();
	}

	const std::string DoubleParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<double>(entry.getAny(false)));
	}

	void DoubleParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<double>(xmlObj.getRequiredDouble("value"), isDefault);
	}

	bool DoubleParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<double>();
	}

	const std::string FloatParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<float>::name();
	}

	const std::string FloatParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<float>(entry.getAny(false)));
	}

	void FloatParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<float>(xmlObj.getRequired<float>("value"), isDefault);
	}

	bool FloatParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<float>();
	}

	const std::string StringParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<std::string>::name();
	}

	const std::string StringParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<std::string>(entry.getAny(false)));
	}

	void StringParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<string>(xmlObj.getRequired("value"), isDefault);
	}

	bool StringParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<std::string>();
	}

	const std::string CharParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<char>::name();
	}

	const std::string CharParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<char>(entry.getAny(false)));
	}

	void CharParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<char>(xmlObj.getRequired<char>("value"), isDefault);
	}

	bool CharParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<char>();
	}

	const std::string BoolParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<bool>::name();
	}

	const std::string BoolParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return toString(any_cast<bool>(entry.getAny(false)));
	}

	void BoolParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<bool>(xmlObj.getRequiredBool("value"), isDefault);
	}

	bool BoolParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<bool>();
	}

	const std::string ArrayIntParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<int> >::name();
	}

	const std::string ArrayIntParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<int> >(entry.getAny(false)).toString();
	}

	void ArrayIntParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<int> >(fromStringToArray<int>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayIntParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<int> >();
	}

	const std::string ArrayShortParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<short> >::name();
	}

	const std::string ArrayShortParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<short> >(entry.getAny(false)).toString();
	}

	void ArrayShortParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<short> >(fromStringToArray<short>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayShortParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<short> >();
	}

	const std::string ArrayLongParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<long> >::name();
	}

	const std::string ArrayLongParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<long> >(entry.getAny(false)).toString();
	}

	void ArrayLongParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<long> >(fromStringToArray<long>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayLongParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<long> >();
	}

	const std::string ArrayDoubleParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<double> >::name();
	}

	const std::string ArrayDoubleParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<double> >(entry.getAny(false)).toString();
	}

	void ArrayDoubleParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<double> >(fromStringToArray<double>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayDoubleParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<double> >();
	}

	const std::string ArrayFloatParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<float> >::name();
	}

	const std::string ArrayFloatParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<float> >(entry.getAny(false)).toString();
	}

	void ArrayFloatParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<float> >(fromStringToArray<float>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayFloatParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<float> >();
	}

	const std::string ArrayStringParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<std::string> >::name();
	}

	const std::string ArrayStringParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<std::string> >(entry.getAny(false)).toString();
	}

	void ArrayStringParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<std::string> >(fromStringToArray<std::string>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayStringParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<std::string> >();
	}

	const std::string ArrayCharParameterEntryConverter::getTypeAttributeValue() const{
		return TypeNameTraits<Array<char> >::name();
	}

	const std::string ArrayCharParameterEntryConverter::getValueAttributeValue(const ParameterEntry &entry) const{
      return any_cast<Array<char> >(entry.getAny(false)).toString();
	}

	void ArrayCharParameterEntryConverter::setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const{
		entry.setValue<Array<char> >(fromStringToArray<char>(xmlObj.getRequired("value")), isDefault);
	}

	bool ArrayCharParameterEntryConverter::isAppropriateConverter(const ParameterEntry& entry) const{
		return entry.isType<Array<char> >();
	}

}


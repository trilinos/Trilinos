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

namespace Teuchos{

	std::string& IntParameterEntryConverter::getTypeAttributeValue(ParameterEntry &entry) const{
		return "int";
	}

	std::string& IntParameterEntryConverter::getValueAttributeValue(ParameterEntry &entry) const{
      return toString(any_cast<int>(entry.getAny(false)));
	}

	void IntParameterEntryConverter::setEntryValue(ParameterEntry &entry, XMLObject, &xmlObj, bool isDefault) const{
		entry.setValue<int>(xmlObj.getRequiredInt("value"), isDefault);
	}

	std::string& ShortParameterEntryConverter::getTypeAttributeValue(ParameterEntry &entry) const{
		return "short";
	}

	std::string& ShortParameterEntryConverter::getValueAttributeValue(ParameterEntry &entry) const{
      return toString(any_cast<short>(entry.getAny(false)));
	}

	void ShortParameterEntryConverter::setEntryValue(ParameterEntry &entry, XMLObject, &xmlObj, bool isDefault) const{
		entry.setValue<short>(xmlObj.getRequiredShort("value"), isDefault);
	}

	std::string& DoubleParameterEntryConverter::getTypeAttributeValue(ParameterEntry &entry) const{
		return "double";
	}

	std::string& DoubleParameterEntryConverter::getValueAttributeValue(ParameterEntry &entry) const{
      return toString(any_cast<double>(entry.getAny(false)));
	}

	void DoubleParameterEntryConverter::setEntryValue(ParameterEntry &entry, XMLObject, &xmlObj, bool isDefault) const{
		entry.setValue<double>(xmlObj.getRequiredDouble("value"), isDefault);
	}

	std::string& FloatParameterEntryConverter::getTypeAttributeValue(ParameterEntry &entry) const{
		return "float";
	}

	std::string& FloatParameterEntryConverter::getValueAttributeValue(ParameterEntry &entry) const{
      return toString(any_cast<float>(entry.getAny(false)));
	}

	void FloatParameterEntryConverter::setEntryValue(ParameterEntry &entry, XMLObject, &xmlObj, bool isDefault) const{
		entry.setValue<float>(xmlObj.getRequiredFloat("value"), isDefault);
	}

	std::string& StringParameterEntryConverter::getTypeAttributeValue(ParameterEntry &entry) const{
		return "string";
	}

	std::string& StringParameterEntryConverter::getValueAttributeValue(ParameterEntry &entry) const{
      return toString(any_cast<std::string>(entry.getAny(false)));
	}

	void StringParameterEntryConverter::setEntryValue(ParameterEntry &entry, XMLObject, &xmlObj, bool isDefault) const{
		entry.setValue<string>(xmlObj.getRequired("value"), isDefault);
	}

	std::string& CharParameterEntryConverter::getTypeAttributeValue(ParameterEntry &entry) const{
		return "char";
	}

	std::string& CharParameterEntryConverter::getValueAttributeValue(ParameterEntry &entry) const{
      return toString(any_cast<char>(entry.getAny(false)));
	}

	void CharParameterEntryConverter::setEntryValue(ParameterEntry &entry, XMLObject, &xmlObj, bool isDefault) const{
		entry.setValue<char>(xmlObj.getRequired("value"), isDefault);
	}


}


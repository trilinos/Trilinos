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

#ifndef TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP
#define TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP

/*! \file Teuchos_ParameterEntryXMLCoverter.hpp
*/

#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_XMLObject.hpp"

namespace Teuchos {

class ParameterEntryXMLConverter{
public:
	virtual ParameterEntry fromXMLtoParameterEntry(const XMLObject &xmlObj) const;
	virtual XMLObject fromParameterEntrytoXML(const ParameterEntry &entry, const std::string &name) const;
	virtual std::string getTypeAttributeValue() const=0;
	virtual std::string getValueAttributeValue(const ParameterEntry &entry) const=0;
	virtual void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const=0;
	virtual bool isAppropriateConverter(const ParameterEntry& entry) const=0;
	static const std::string& getTypeAttributeName(){
		return typeAttributeName();
	}
private:
	static std::string& defaultAttributeName(){
		static std::string attributeName;
		if(attributeName == ""){
			attributeName = "isDefault";
		}
		return attributeName;
	}

	static std::string& usedAttributeName(){
		static std::string attributeName;
		if(attributeName == ""){
			attributeName = "isUsed";
		}
		return attributeName;
	}

	static std::string& nameAttributeName(){
		static std::string attributeName;
		if(attributeName == ""){
			attributeName = "name";
		}
		return attributeName;
	}

	static std::string& typeAttributeName(){
		static std::string attributeName;
		if(attributeName == ""){
			attributeName = "type";
		}
		return attributeName;
	}

	static std::string& valueAttributeName(){
		static std::string attributeName;
		if(attributeName == ""){
			attributeName = "value";
		}
		return attributeName;
	}
};

}

#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP

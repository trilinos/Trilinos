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

/* \class Teuchos::ParameterEntryXMLConverter
 * \brief A class used to convert parameter entries to xml and vice versa.
 *
 */
class ParameterEntryXMLConverter{
public:
	/* \brief Converts the given xml into a parameter entry. 
	 *
	 * @param xmlObj The xml to be converted to a parameter entry.
	 * @return A ParameterEntry with the aspects specified by the xml.
	 */
	ParameterEntry fromXMLtoParameterEntry(const XMLObject &xmlObj) const;

	/* \brief Converts the given parameter entry to xml.
	 *
	 * @param entry The parameter entry to convert to xml.
	 * @param name The name associated with the parameter entry.
	 * @return An XMLObject representing the parameter entry.
	 */
	XMLObject fromParameterEntrytoXML(const ParameterEntry &entry, const std::string &name) const;

	/* \brief Gets a string representing the value that should be assigned to the "type" attribute when converting
	 * a parameter entry to xml.
	 *
	 * @return The value to be assigned to the "type" attribute when converting a parameter entry to xml.
	 */
	virtual const std::string getTypeAttributeValue() const=0;

	/* \brief Gets the value to be assigned to the "value" attribute when converting the paramter entry to xml.
	 *
	 * @param entry The entry being converted.
	 * @areturn The value to be assigned to the "value" attribute when converting the parameter entry to xml.
	 */
	virtual const std::string getValueAttributeValue(const ParameterEntry &entry) const=0;

	/* \brief sets the value  */
	virtual void setEntryValue(ParameterEntry &entry, const XMLObject &xmlObj, bool isDefault) const=0;

	/* \brief determines wether or not this converter is appropriate for the given parameter entry.*/
	virtual bool isAppropriateConverter(const ParameterEntry& entry) const=0;

	static const std::string& getTypeAttributeName(){
		static const std::string typeAttributeName_ = "type";
		return typeAttributeName_;
	}

	static const std::string& getValueAttributeName(){
		static const std::string valueAttributeName_ = "value";
		return valueAttributeName_;
	}

private:
	static const std::string& getDefaultAttributeName(){
		static const std::string defaultAttributeName_ = "isDefault";
		return defaultAttributeName_;
	}

	static const std::string& getUsedAttributeName(){
		static const std::string usedAttributeName_ = "isUsed";
		return usedAttributeName_;
	}

	static const std::string& getNameAttributeName(){
		static const std::string nameAttributeName_ = "name";
		return nameAttributeName_;
	}
};

}

#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP

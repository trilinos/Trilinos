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


#ifndef TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
#define TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
/*! \file Teuchos_ParameterEntryXMLCoverter.hpp
*/
#include "Teuchos_ValidatorXMLConverter.hpp"


namespace Teuchos {
class ParameterEntryValidator;

class ValidatorXMLConverterDB{
public:
	static void addConverter(ParameterEntryValidator& validator, RCP<ValidatorXMLConverter> converterToAdd);

	static RCP<const ValidatorXMLConverter> getConverter(const ParameterEntryValidator& validator);

	static RCP<const ValidatorXMLConverter> getConverter(const XMLObject& xmlObject);

private:
	typedef std::map<std::string, RCP<ValidatorXMLConverter> > ConverterMap;
	typedef std::pair<std::string, RCP<ValidatorXMLConverter> > ConverterPair;

	static RCP<ValidatorXMLConverter> getDefaultConverter();

	static ConverterMap& getConverterMap();
};

} // end namespace Teuchos

#endif // TEUCHOS_VALIDATORXMLCONVERTERDB_HPP

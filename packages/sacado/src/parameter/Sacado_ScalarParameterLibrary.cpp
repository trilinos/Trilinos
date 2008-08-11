// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Sacado_ScalarParameterLibrary.hpp"

void
Sacado::ScalarParameterLibrary::
setRealValueForAllTypes(const std::string& name, double value)
{
  FamilyMap::iterator it = library.find(name);
  TEST_FOR_EXCEPTION(
     it == library.end(), 
     std::logic_error,
     std::string("Sacado::ScalararameterLibrary::setRealValueForAllTypes():  ")
     + "Invalid parameter family " + name);
  (*it).second->setRealValueForAllTypes(value);
}

void
Sacado::ScalarParameterLibrary::
fillVector(const Teuchos::Array<std::string>& names,
	   Sacado::ScalarParameterVector& pv)
{
  FamilyMap::iterator it;

  // Fill in parameters
  for (unsigned int i=0; i<names.size(); i++) {
    it = library.find(names[i]);
    TEST_FOR_EXCEPTION(
		   it == library.end(), 
		   std::logic_error,
		   std::string("Sacado::ParameterLibraryBase::fillVector():  ")
		   + "Invalid parameter family " + names[i]);
    pv.addParam((*it).second, 0.0);
    pv[i].baseValue = (*it).second->getValue<double>();
  }
}

// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Parameter_Library.H"
#include "LOCA_Parameter_Entry.H"

LOCA::Parameter::Library::~Library() {
  ParameterMapIterator paramIt;
  ValueTypeMapIterator valueIt;

  // Loop over all parameter entries
  for (paramIt = library.begin(); paramIt != library.end(); paramIt++) {

    // Loop over all value entries
    for (valueIt = paramIt->second->begin(); 
	 valueIt != paramIt->second->end(); 
	 valueIt++) {

      // Delete entry
      delete valueIt->second;

    }

    delete paramIt->second;

  }

}

LOCA::Parameter::Library::Library(const LOCA::Parameter::Library& l) 
  : library(l.library) {}

LOCA::Parameter::Library&
LOCA::Parameter::Library::operator = (const LOCA::Parameter::Library& l) {
  library = l.library;
  return *this;
}

LOCA::Parameter::Library::ParameterMapIterator
LOCA::Parameter::Library::getEntryMapIterator(const string& name) {
  return library.find(name);
}

LOCA::Parameter::Library::ValueTypeMapIterator
LOCA::Parameter::Library::getEntryIterator(
				   const string& valueTypeString, 
				   const ParameterMapIterator& paramIterator) {
  return paramIterator->second->find(valueTypeString);
}

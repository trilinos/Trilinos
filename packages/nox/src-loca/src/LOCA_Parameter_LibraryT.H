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
#include <typeinfo>
#include "LOCA_Parameter_Entry.H"
#include "LOCA_ErrorCheck.H"

template <class ValueType>
void
LOCA::Parameter::Library::setParameterValue(const string& name, 
					    const ValueType& value) {
  // Get entry 
  LOCA::Parameter::Entry<ValueType>& entry = getEntry<ValueType>(name);

  // Call entry's setValue method to set value
  entry.setValue(value);
}

template <class ObjectType, class ValueType>
bool
LOCA::Parameter::Library::addParameterEntry(
				      const string& name, 
				      ObjectType& object, 
				      ValueType ObjectType::* object_val_ptr) {

  typedef LOCA::Parameter::DefaultFunctor<ObjectType, ValueType> FunctorType;

  // Instantiate a default functor
  FunctorType *fctr = new FunctorType(object, object_val_ptr);

  // Instantiate a standard entry
  LOCA::Parameter::StandardEntry<FunctorType, ValueType> *entry = 
    new LOCA::Parameter::StandardEntry<FunctorType, ValueType>(fctr);

  // Attempt to add entry
  bool success = addParameterEntry(name, entry);

  // If unsuccessful, delete entry and functor
  if (!success) {
    delete fctr;
    delete entry;
  }

  return success;
}

template <class FunctorType, class ValueType>
bool
LOCA::Parameter::Library::addParameterEntry(const string& name, 
					    FunctorType* fctr) {
  
  // Instantiate a standard entry
  LOCA::Parameter::StandardEntry<FunctorType, ValueType> *entry = 
    new LOCA::Parameter::StandardEntry<FunctorType, ValueType>(fctr);

  // Attempt to add entry
  bool success = addParameterEntry(name, entry);

  // If unsuccessful, delete entry
  if (!success)
    delete entry;

  return success;
}

template <class ValueType>
bool
LOCA::Parameter::Library::addParameterEntry(
				     const string& name, 
				     LOCA::Parameter::Entry<ValueType>* entry)
{
  // Get string representation of ValueType
  string valueTypeString = getTypeName<ValueType>();

  // Determine if a map already exists for this parameter
  ParameterMapIterator paramIterator = getEntryMapIterator(name);

  // If it doesn not, add it
  if (paramIterator == library.end()) {
    ValueTypeMap* paramMap = new ValueTypeMap;
    entry->setIsInLibrary();
    paramMap->insert(pair<string, AbstractEntry*>(valueTypeString, entry));
    library.insert(pair<string, ValueTypeMap*>(name, paramMap));
  }
  
  else {
    
    // Determine if entry already exists for parameter and type
    ValueTypeMapIterator valueIterator = getEntryIterator(valueTypeString,
							  paramIterator);

    // If it does not, add it
    if (valueIterator == paramIterator->second->end()) {
      entry->setIsInLibrary();
      paramIterator->second->insert(pair<string, AbstractEntry*>(name, entry));
    }

    else {
      return false;
    }
    
  }

  return true;
}

template <class ValueType>
string
LOCA::Parameter::Library::getTypeName() {
  return typeid(ValueType).name();
}   

template <class ValueType>
LOCA::Parameter::Entry<ValueType>&
LOCA::Parameter::Library::getEntry(const string& name) {

  const char *func = "LOCA::Parameter::Library::getEntry()";

  // Find map for this parameter name
  ParameterMapIterator paramIterator = getEntryMapIterator(name);
  if (paramIterator == library.end()) {
    LOCA::ErrorCheck::throwError(
	   func,
	   string("Parameter library has no set of entries for parameter ") + 
	   name);
  }

  // Convert typename ValueType to string
  string valueTypeString = getTypeName<ValueType>();

  // Find entry corresponding to this ValueType
  ValueTypeMapIterator valueIterator = getEntryIterator(valueTypeString, 
							paramIterator);
  if (valueIterator == paramIterator->second->end()) {
    LOCA::ErrorCheck::throwError(
	   func,
	   string("Parameter library has no entry for parameter ") + 
	   name + string(" of type ") + valueTypeString);
  }

  // Cast entry to LOCA::Parameter::Entry<ValueType>
  LOCA::Parameter::Entry<ValueType>* entry = 
    dynamic_cast<LOCA::Parameter::Entry<ValueType>*>(valueIterator->second);
  if (entry == NULL) {
    LOCA::ErrorCheck::throwError(
	   func,
	   string("Parameter entry for parameter ") + name + 
	   string(" and type ") + valueTypeString + 
	   string(" is not of the right type!"));
  }

  return *entry;
}

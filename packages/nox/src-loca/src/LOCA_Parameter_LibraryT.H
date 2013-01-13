// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
#include <typeinfo>
#include "LOCA_Parameter_Entry.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

#include "Teuchos_Assert.hpp" // for errors
#include <utility>

template <class ValueType>
void
LOCA::Parameter::Library::setValue(const std::string& name, 
					    const ValueType& value) {
  // Get entry 
  LOCA::Parameter::Entry<ValueType>& entry = getEntry<ValueType>(name);

  // Call entry's setValue method to set value
  entry.setValue(value);
}

template <class ValueType>
ValueType
LOCA::Parameter::Library::getValue(const std::string& name) const {

  // Get entry 
  const LOCA::Parameter::Entry<ValueType>& entry = getEntry<ValueType>(name);

  // Call entry's getValue method to get value
  return entry.getValue();
}

template <class ObjectType, class ValueType>
bool
LOCA::Parameter::Library::addParameterEntry(
				      const std::string& name, 
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
LOCA::Parameter::Library::addParameterEntry(const std::string& name, 
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
				     const std::string& name, 
				     LOCA::Parameter::Entry<ValueType>* entry)
{
  // Get std::string representation of ValueType
  std::string valueTypeString = getTypeName<ValueType>();

  // Determine if a map already exists for this parameter
  ParameterMapIterator paramIterator = getEntryMapIterator(name);

  // If it doesn not, add it
  if (paramIterator == library.end()) {
    ValueTypeMap* paramMap = new ValueTypeMap;
    entry->setIsInLibrary();
    (*paramMap)[valueTypeString] = entry;
    library[name] = paramMap;
  }
  
  else {
    
    // Determine if entry already exists for parameter and type
    ValueTypeMapIterator valueIterator = getEntryIterator(valueTypeString,
							  paramIterator);

    // If it does not, add it
    if (valueIterator == (*paramIterator).second->end()) {
      entry->setIsInLibrary();
      (*paramIterator).second->insert(std::pair<std::string, AbstractEntry*>(name, entry));
    }

    else {
      return false;
    }
    
  }

  return true;
}

template <class ValueType>
std::string
LOCA::Parameter::Library::getTypeName() const {
  return typeid(ValueType).name();
}   

template <class ValueType>
LOCA::Parameter::Entry<ValueType>&
LOCA::Parameter::Library::getEntry(const std::string& name) {

  const char *func = "LOCA::Parameter::Library::getEntry(): ";

  // Find map for this parameter name
  ParameterMapIterator paramIterator = getEntryMapIterator(name);
  TEUCHOS_TEST_FOR_EXCEPTION(paramIterator == library.end(), 
		     std::invalid_argument,
		     func << "Parameter library has no set of entries for " << 
		     "parameter " << name);

  // Convert typename ValueType to std::string
  std::string valueTypeString = getTypeName<ValueType>();

  // Find entry corresponding to this ValueType
  ValueTypeMapIterator valueIterator = getEntryIterator(valueTypeString, 
							paramIterator);
  TEUCHOS_TEST_FOR_EXCEPTION(valueIterator == (*paramIterator).second->end(), 
		     std::invalid_argument,
		     func << "Parameter library has no entry for parameter " <<
		     name << " of type " << valueTypeString);

  // Cast entry to LOCA::Parameter::Entry<ValueType>
  LOCA::Parameter::Entry<ValueType>* entry = 
    dynamic_cast<LOCA::Parameter::Entry<ValueType>*>((*valueIterator).second);
  TEUCHOS_TEST_FOR_EXCEPTION(entry == NULL, 
		     std::invalid_argument,
		     func << "Parameter entry for parameter " << name <<  
		     " and type " << valueTypeString <<
		     " is not of the right type!");

  return *entry;
}

template <class ValueType>
const LOCA::Parameter::Entry<ValueType>&
LOCA::Parameter::Library::getEntry(const std::string& name) const {

  const char *func = "LOCA::Parameter::Library::getEntry(): ";

  // Find map for this parameter name
  ParameterMapConstIterator paramIterator = getEntryMapIterator(name);
  TEUCHOS_TEST_FOR_EXCEPTION(paramIterator == library.end(), 
		     std::invalid_argument,
		     func << "Parameter library has no set of entries for " << 
		     "parameter " << name);

  // Convert typename ValueType to std::string
  std::string valueTypeString = getTypeName<ValueType>();

  // Find entry corresponding to this ValueType
  ValueTypeMapConstIterator valueIterator = getEntryIterator(valueTypeString, 
							     paramIterator);
  TEUCHOS_TEST_FOR_EXCEPTION(valueIterator == (*paramIterator).second->end(), 
		     std::invalid_argument,
		     func << "Parameter library has no entry for parameter " <<
		     name << " of type " << valueTypeString);

  // Cast entry to LOCA::Parameter::Entry<ValueType>
  const LOCA::Parameter::Entry<ValueType>* entry = 
    dynamic_cast<const LOCA::Parameter::Entry<ValueType>*>((*valueIterator).second);
  TEUCHOS_TEST_FOR_EXCEPTION(entry == NULL, 
		     std::invalid_argument,
		     func << "Parameter entry for parameter " << name <<  
		     " and type " << valueTypeString <<
		     " is not of the right type!");

  return *entry;
}

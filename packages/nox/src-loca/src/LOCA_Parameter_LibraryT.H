// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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

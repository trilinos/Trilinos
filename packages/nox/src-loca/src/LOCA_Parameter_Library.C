// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Parameter_Library.H"
#include "LOCA_Parameter_Entry.H"

LOCA::Parameter::Library::~Library() {
  ParameterMapIterator paramIt;
  ValueTypeMapIterator valueIt;

  // Loop over all parameter entries
  for (paramIt = library.begin(); paramIt != library.end(); paramIt++) {

    // Loop over all value entries
    for (valueIt = (*paramIt).second->begin();
     valueIt != (*paramIt).second->end();
     valueIt++) {

      // Delete entry
      delete (*valueIt).second;

    }

    delete (*paramIt).second;

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
LOCA::Parameter::Library::getEntryMapIterator(const std::string& name) {
  return library.find(name);
}

LOCA::Parameter::Library::ParameterMapConstIterator
LOCA::Parameter::Library::getEntryMapIterator(const std::string& name) const {
  return library.find(name);
}

LOCA::Parameter::Library::ValueTypeMapIterator
LOCA::Parameter::Library::getEntryIterator(
                   const std::string& valueTypeString,
                   const ParameterMapIterator& paramIterator) {
  return (*paramIterator).second->find(valueTypeString);
}

LOCA::Parameter::Library::ValueTypeMapConstIterator
LOCA::Parameter::Library::getEntryIterator(
               const std::string& valueTypeString,
               const ParameterMapConstIterator& paramIterator) const {
  return (*paramIterator).second->find(valueTypeString);
}

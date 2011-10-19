// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_StandardParameterEntryValidators.hpp"

std::string Teuchos::getVerbosityLevelParameterValueName(
  const EVerbosityLevel verbLevel
  )
{
  switch(verbLevel) {
    case VERB_DEFAULT:
      return "default";
    case VERB_NONE:
      return "none";
    case VERB_LOW:
      return "low";
    case VERB_MEDIUM:
      return "medium";
    case VERB_HIGH:
      return "high";
    case VERB_EXTREME:
      return "extreme";
    default:
      TEST_FOR_EXCEPT("Should never get here!");
  }
  return ""; // Never get here!
}


Teuchos::RCP<
  Teuchos::StringToIntegralParameterEntryValidator<Teuchos::EVerbosityLevel>
  >
Teuchos::verbosityLevelParameterEntryValidator(
  std::string const& defaultParameterName
  )
{
  return rcp(
    new StringToIntegralParameterEntryValidator<EVerbosityLevel>(
      tuple<std::string>(
        getVerbosityLevelParameterValueName(VERB_DEFAULT),
        getVerbosityLevelParameterValueName(VERB_NONE),
        getVerbosityLevelParameterValueName(VERB_LOW),
        getVerbosityLevelParameterValueName(VERB_MEDIUM),
        getVerbosityLevelParameterValueName(VERB_HIGH),
        getVerbosityLevelParameterValueName(VERB_EXTREME)
        ),
      tuple<std::string>(
        "Use level set in code",
        "Produce no output",
        "Produce minimal output",
        "Produce a little more output",
        "Produce a higher level of output",
        "Produce the highest level of output"
        ),
      tuple<EVerbosityLevel>(
        VERB_DEFAULT,
        VERB_NONE,
        VERB_LOW,
        VERB_MEDIUM,
        VERB_HIGH,
        VERB_EXTREME
        ),
      defaultParameterName
      )
    );
}


namespace Teuchos {


//
// AnyNumberParameterEntryValidator
//


// Constructors


AnyNumberParameterEntryValidator::AnyNumberParameterEntryValidator()
  : preferredType_(PREFER_DOUBLE), acceptedTypes_(AcceptedTypes())
{
  finishInitialization();
}


AnyNumberParameterEntryValidator::AnyNumberParameterEntryValidator(
  EPreferredType const preferredType, AcceptedTypes const& acceptedTypes
  )
  : preferredType_(preferredType), acceptedTypes_(acceptedTypes)
{
  finishInitialization();
}

//  Local non-virtual validated lookup functions


int AnyNumberParameterEntryValidator::getInt(
  const ParameterEntry &entry, const std::string &paramName,
  const std::string &sublistName, const bool activeQuery
  ) const
{
  const any &anyValue = entry.getAny(activeQuery);
  if( acceptedTypes_.allowInt() && anyValue.type() == typeid(int) )
    return any_cast<int>(anyValue);
  if( acceptedTypes_.allowDouble() && anyValue.type() == typeid(double) )
    return static_cast<int>(any_cast<double>(anyValue));
  if( acceptedTypes_.allowString() && anyValue.type() == typeid(std::string) )
    return std::atoi(any_cast<std::string>(anyValue).c_str());
  throwTypeError(entry,paramName,sublistName);
  return 0; // Will never get here!
}


double AnyNumberParameterEntryValidator::getDouble(
  const ParameterEntry &entry, const std::string &paramName,
  const std::string &sublistName, const bool activeQuery
  ) const
{
  const any &anyValue = entry.getAny(activeQuery);
  if( acceptedTypes_.allowInt() && anyValue.type() == typeid(int) )
    return static_cast<double>(any_cast<int>(anyValue));
  if( acceptedTypes_.allowDouble() && anyValue.type() == typeid(double) )
    return any_cast<double>(anyValue);
  if( acceptedTypes_.allowString() && anyValue.type() == typeid(std::string) )
    return std::atof(any_cast<std::string>(anyValue).c_str());
  throwTypeError(entry,paramName,sublistName);
  return 0.0; // Will never get here!
}


std::string AnyNumberParameterEntryValidator::getString(
  const ParameterEntry &entry, const std::string &paramName,
  const std::string &sublistName, const bool activeQuery
  ) const
{
  const any &anyValue = entry.getAny(activeQuery);
  if( acceptedTypes_.allowInt() && anyValue.type() == typeid(int) )
    return Utils::toString(any_cast<int>(anyValue));
  if( acceptedTypes_.allowDouble() && anyValue.type() == typeid(double) )
    return Utils::toString(any_cast<double>(anyValue));
  if( acceptedTypes_.allowString() && anyValue.type() == typeid(std::string) )
    return any_cast<std::string>(anyValue);
  throwTypeError(entry,paramName,sublistName);
  return ""; // Will never get here!
}


int AnyNumberParameterEntryValidator::getInt(
  ParameterList &paramList, const std::string &paramName,
  const int defaultValue
  ) const
{
  const ParameterEntry *entry = paramList.getEntryPtr(paramName);
  if(entry) return getInt(*entry,paramName,paramList.name(),true);
  return paramList.get(paramName,defaultValue);
}


double AnyNumberParameterEntryValidator::getDouble(
  ParameterList &paramList, const std::string &paramName,
  const double defaultValue
  ) const
{
  const ParameterEntry *entry = paramList.getEntryPtr(paramName);
  if(entry) return getDouble(*entry,paramName,paramList.name(),true);
  return paramList.get(paramName,defaultValue);
}


std::string AnyNumberParameterEntryValidator::getString(
  ParameterList &paramList, const std::string &paramName,
  const std::string &defaultValue
  ) const
{
  const ParameterEntry *entry = paramList.getEntryPtr(paramName);
  if(entry) return getString(*entry,paramName,paramList.name(),true);
  return paramList.get(paramName,defaultValue);
}


bool AnyNumberParameterEntryValidator::isDoubleAllowed() const
{
  return acceptedTypes_.allowDouble();
}


bool AnyNumberParameterEntryValidator::isIntAllowed() const
{
  return acceptedTypes_.allowInt();
}

  
bool AnyNumberParameterEntryValidator::isStringAllowed() const
{
  return acceptedTypes_.allowString();
}


AnyNumberParameterEntryValidator::EPreferredType
AnyNumberParameterEntryValidator::getPreferredType() const
{
  return preferredType_;
}


// Overridden from ParameterEntryValidator


const std::string AnyNumberParameterEntryValidator::getXMLTypeName() const
{
  return "anynumberValidator";
}


void AnyNumberParameterEntryValidator::printDoc(
  std::string  const & docString,
  std::ostream & out
  ) const
{
  StrUtils::printLines(out,"# ",docString);
  out << "#  Accepted types: " << acceptedTypesString_ << ".\n";
}


ParameterEntryValidator::ValidStringsList
AnyNumberParameterEntryValidator::validStringValues() const
{
  return null;
}


void AnyNumberParameterEntryValidator::validate(
  ParameterEntry const& entry,
  std::string const& paramName,
  std::string const& sublistName
  ) const
{
  // Validate (any of the get functions will do!)
  getInt(entry,paramName,sublistName,false);
}


void AnyNumberParameterEntryValidator::validateAndModify(
  std::string const& paramName,
  std::string const& sublistName,
  ParameterEntry * entry
  ) const
{
  TEST_FOR_EXCEPT(0==entry);
  switch(preferredType_) {
    case PREFER_INT:
      entry->setValue(
        getInt(*entry,paramName,sublistName,false),
        false // isDefault
        );
      break;
    case PREFER_DOUBLE:
      entry->setValue(
        getDouble(*entry,paramName,sublistName,false),
        false // isDefault
        );
      break;
    case PREFER_STRING:
      entry->setValue(
        getString(*entry,paramName,sublistName,false),
        false // isDefault
        );
      break;
    default:
      TEST_FOR_EXCEPT("Error, Invalid EPreferredType value!");
  }
}


// private


void AnyNumberParameterEntryValidator::finishInitialization()
{

  std::ostringstream oss;
  bool addedType = false;
  if(acceptedTypes_.allowInt()) {
    oss << "\"int\"";
    addedType = true;
  }
  if(acceptedTypes_.allowDouble()) {
    if(addedType) oss << ", ";
    oss << "\"double\"";
    addedType = true;
  }
  if(acceptedTypes_.allowString()) {
    if(addedType) oss << ", ";
    oss << "\"string\"";
    addedType = true;
  }
  acceptedTypesString_ = oss.str();
}


void AnyNumberParameterEntryValidator::throwTypeError(
  ParameterEntry const& entry,
  std::string const& paramName,
  std::string const& sublistName
  ) const
{
  const std::string &entryName = entry.getAny(false).typeName();
  TEST_FOR_EXCEPTION_PURE_MSG(
    true, Exceptions::InvalidParameterType
    ,"Error, the parameter {paramName=\""<<paramName<<"\""
    ",type=\""<<entryName<<"\"}"
    << "\nin the sublist \"" << sublistName << "\""
    << "\nhas the wrong type."
    << "\n\nThe accepted types are: " << acceptedTypesString_ << "!";
    );
}

RCP<AnyNumberParameterEntryValidator> 
  DummyObjectGetter<AnyNumberParameterEntryValidator>::getDummyObject()
{
  return anyNumberParameterEntryValidator(
    AnyNumberParameterEntryValidator::PREFER_INT, 
    AnyNumberParameterEntryValidator::AcceptedTypes());
}

FileNameValidator::FileNameValidator(bool mustAlreadyExist)
  : ParameterEntryValidator(), mustAlreadyExist_(mustAlreadyExist)
{}

bool FileNameValidator::fileMustExist() const
{
  return mustAlreadyExist_;
}


bool FileNameValidator::setFileMustExist(bool shouldFileExist)
{
  this->mustAlreadyExist_ = shouldFileExist;
  return mustAlreadyExist_;
}


ParameterEntryValidator::ValidStringsList
  FileNameValidator::validStringValues() const
{
  return null;
}


void FileNameValidator::validate(ParameterEntry const &entry, std::string const &paramName,
  std::string const &sublistName) const
{
  const std::string &entryName = entry.getAny(false).typeName();
  any anyValue = entry.getAny(true);
  TEST_FOR_EXCEPTION(!(anyValue.type() == typeid(std::string) ),
    Exceptions::InvalidParameterType,
    "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
    " parameter in the \"" << sublistName << 
    "\" sublist didn't quite work out." << std::endl << std::endl <<
    "No need to fret though. I'm sure it's just a small mistake. "
    "Maybe the information below "
    "can help you figure out what went wrong." << std::endl << std::endl <<
    "Error: The value that you entered was the wrong type." << std::endl <<
    "Parameter: " << paramName << std::endl << 
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << typeid(std::string).name() << 
    std::endl << std::endl);
  if(mustAlreadyExist_){
    std::string fileName = getValue<std::string>(entry);
    TEST_FOR_EXCEPTION(!std::ifstream(fileName.c_str()),
      Exceptions::InvalidParameterValue,
      "Aww shoot! Sorry bud, but it looks like the \"" 
      << paramName << "\"" <<
      " parameter in the \"" << sublistName << 
      "\" sublist didn't quite work out." << std::endl << 
      "No need to fret though. I'm sure it's just a small mistake. " <<
      "Maybe the information below "<<
      "can help you figure out what went wrong." << 
      std::endl << std::endl <<
      "Error: The file must already exists. The value you entered does " <<
      "not corresspond to an existing file name." << std::endl <<
      "Parameter: " << paramName << std::endl <<
      "File name specified: " << fileName << std::endl << std::endl);
  }
}


const std::string FileNameValidator::getXMLTypeName() const
{
  return "FilenameValidator";
}


void FileNameValidator::printDoc(
  std::string const &docString, std::ostream &out) const
{
  StrUtils::printLines(out,"# ",docString);
  out << "#  Validator Used: " << std::endl;
  out << "#  FileName Validator" << std::endl;
}

RCP<FileNameValidator> DummyObjectGetter<FileNameValidator>::getDummyObject(){
  return rcp(new FileNameValidator(true));
}

StringValidator::StringValidator()
  : ParameterEntryValidator(), validStrings_(NULL)
{}

StringValidator::StringValidator(const Array<std::string>& validStrings):
  ParameterEntryValidator(),
  validStrings_(rcp(new Array<std::string>(validStrings)))
{}

ParameterEntryValidator::ValidStringsList
StringValidator::setValidStrings(const Array<std::string>& validStrings)
{
  validStrings_ = rcp(new Array<std::string>(validStrings));
  return validStrings_;
}


ParameterEntryValidator::ValidStringsList
StringValidator::validStringValues() const
{
  return validStrings_;
}


void StringValidator::validate(
  ParameterEntry const &entry, std::string const &paramName,
  std::string const &sublistName) const
{
  any anyValue = entry.getAny(true);
  const std::string &entryName = entry.getAny(false).typeName();
  TEST_FOR_EXCEPTION(!(anyValue.type() == typeid(std::string)) ,
    Exceptions::InvalidParameterType,
    "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
    " parameter in the \"" << sublistName << 
    "\" sublist didn't quite work out." << std::endl <<
    "No need to fret though. I'm sure it's just a small mistake. " <<
    "Maybe the information below "<<
    "can help you figure out what went wrong." << std::endl << std::endl <<
    "Error: The value that you entered was the wrong type." <<
    "Parameter: " << paramName << std::endl <<
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << Teuchos::TypeNameTraits<std::string>::name() <<
    std::endl);
  if(!validStrings_.is_null()){
    Array<std::string>::const_iterator
      it = std::find(validStrings_->begin(),
      validStrings_->end(), getValue<std::string>(entry));
    TEST_FOR_EXCEPTION(it == validStrings_->end(),
      Exceptions::InvalidParameterValue,
      "Aww shoot! Sorry bud, but it looks like the \""
      << paramName << "\"" <<
      " parameter in the \"" << sublistName << 
      "\" sublist didn't quite work out." << std::endl <<
      "No need to fret though. I'm sure it's just a small mistake. "
      "Maybe the information below "
      "can help you figure out what went wrong." <<
      std::endl << std::endl <<
      "Error: The value that was entered doesn't fall with in "
      "the range set by the validator." <<
      "Parameter: " << paramName << std::endl <<
      "Acceptable Values: " << validStrings_ << std::endl <<
      "Value entered: " << getValue<std::string>(entry) << std::endl <<
      std::endl);
  }
}


const std::string StringValidator::getXMLTypeName() const
{
  return "StringValidator";
}


void StringValidator::printDoc(std::string const &docString,
  std::ostream &out) const
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#  Validator Used: " << std::endl;
  out << "#  String Validator" << std::endl;
}

RCP<StringValidator> DummyObjectGetter<StringValidator>::getDummyObject(){
  return rcp(new StringValidator(tuple<std::string>("")));
}


} // namespace Teuchos


// Nonmmeber helper functions

Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
Teuchos::anyNumberParameterEntryValidator()
{
  return rcp(new AnyNumberParameterEntryValidator());
}


Teuchos::RCP<Teuchos::AnyNumberParameterEntryValidator>
Teuchos::anyNumberParameterEntryValidator(
  AnyNumberParameterEntryValidator::EPreferredType const preferredType,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  )
{
  return rcp(
    new AnyNumberParameterEntryValidator(
      preferredType, acceptedTypes
      )
    );
}

void Teuchos::setIntParameter(
  std::string const& paramName,
  int const value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  )
{
  TEST_FOR_EXCEPT(0==paramList);
  const RCP<const ParameterEntryValidator> paramEntryValidator =  
    anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::PREFER_INT, acceptedTypes
      );
  paramList->set(paramName, value, docString, paramEntryValidator);
}


void Teuchos::setDoubleParameter(
  std::string const& paramName,
  double const& value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  )
{
  TEST_FOR_EXCEPT(0==paramList);
  const RCP<const ParameterEntryValidator> paramEntryValidator =  
    anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::PREFER_DOUBLE, acceptedTypes
      );
  paramList->set(paramName, value, docString, paramEntryValidator);
}


void Teuchos::setNumericStringParameter(
  std::string const& paramName,
  std::string const& value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  )
{
  TEST_FOR_EXCEPT(0==paramList);
  const RCP<const ParameterEntryValidator> paramEntryValidator =  
    anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::PREFER_STRING, acceptedTypes
      );
  paramList->set(paramName, value, docString, paramEntryValidator);
}


int Teuchos::getIntParameter(
  ParameterList const& paramList,
  std::string const& paramName
  )
{
  const ParameterEntry &entry = paramList.getEntry(paramName);
  RCP<const AnyNumberParameterEntryValidator>
    anyNumValidator = rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      entry.validator()
      );
  if ( !is_null(anyNumValidator) )
    return anyNumValidator->getInt(entry,paramName,paramList.name());
  if ( typeid(int) == entry.getAny().type() )
    return any_cast<int>(entry.getAny());
  // Try the do the conversion which might fail!
  const AnyNumberParameterEntryValidator myAnyNumValidator;
  return myAnyNumValidator.getInt(entry,paramName,paramList.name());
}


double Teuchos::getDoubleParameter(
  ParameterList const& paramList,
  std::string const& paramName
  )
{
  const ParameterEntry &entry = paramList.getEntry(paramName);
  RCP<const AnyNumberParameterEntryValidator>
    anyNumValidator = rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      entry.validator()
      );
  if ( !is_null(anyNumValidator) )
    return anyNumValidator->getDouble(entry,paramName,paramList.name());
  if ( typeid(double) == entry.getAny().type() )
    return any_cast<double>(entry.getAny());
  // Try the do the conversion which might fail!
  const AnyNumberParameterEntryValidator myAnyNumValidator;
  return myAnyNumValidator.getDouble(entry,paramName,paramList.name());
}


std::string Teuchos::getNumericStringParameter(
  ParameterList const& paramList,
  std::string const& paramName
  )
{
  const ParameterEntry &entry = paramList.getEntry(paramName);
  RCP<const AnyNumberParameterEntryValidator>
    anyNumValidator = rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      entry.validator()
      );
  if ( !is_null(anyNumValidator) )
    return anyNumValidator->getString(entry,paramName,paramList.name());
  if ( typeid(std::string) == entry.getAny().type() )
    return any_cast<std::string>(entry.getAny());
  // Try the do the conversion which might fail!
  const AnyNumberParameterEntryValidator myAnyNumValidator;
  return myAnyNumValidator.getString(entry,paramName,paramList.name());
}



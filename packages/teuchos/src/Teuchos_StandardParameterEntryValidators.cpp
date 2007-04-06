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

#include "Teuchos_StandardParameterEntryValidators.hpp"

Teuchos::RefCountPtr<
  Teuchos::StringToIntegralParameterEntryValidator<Teuchos::EVerbosityLevel>
  >
Teuchos::verbosityLevelParameterEntryValidator(
  std::string const& defaultParameterName
  )
{
  return rcp(
    new StringToIntegralParameterEntryValidator<Teuchos::EVerbosityLevel>(
      tuple<std::string>(
        "default",
        "none",
        "low",
        "medium",
        "high",
        "extreme"
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
        Teuchos::VERB_DEFAULT,
        Teuchos::VERB_NONE,
        Teuchos::VERB_LOW,
        Teuchos::VERB_MEDIUM,
        Teuchos::VERB_HIGH,
        Teuchos::VERB_EXTREME
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


AnyNumberParameterEntryValidator::AnyNumberParameterEntryValidator(
  EPreferredType const preferredType,
  AcceptedTypes const& acceptedTypes
  )
  :preferredType_(preferredType), acceptedTypes_(acceptedTypes)
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


//  Local non-virtual validated lookup functions


int AnyNumberParameterEntryValidator::getInt(
  const ParameterEntry &entry, const std::string &paramName
  ,const std::string &sublistName, const bool activeQuery
  ) const
{
  const any &anyValue = entry.getAny(activeQuery);
  if( acceptedTypes_.allowInt() && anyValue.type() == typeid(int) )
    return any_cast<int>(anyValue);
  if( acceptedTypes_.allowDouble() && anyValue.type() == typeid(double) )
    return static_cast<int>(any_cast<double>(anyValue));
  if( acceptedTypes_.allowString() && anyValue.type() == typeid(std::string) )
    return ::atoi(any_cast<string>(anyValue).c_str());
  throwTypeError(entry,paramName,sublistName);
  return 0; // Will never get here!
}


double AnyNumberParameterEntryValidator::getDouble(
  const ParameterEntry &entry, const std::string &paramName
  ,const std::string &sublistName, const bool activeQuery
  ) const
{
  const any &anyValue = entry.getAny(activeQuery);
  if( acceptedTypes_.allowInt() && anyValue.type() == typeid(int) )
    return static_cast<double>(any_cast<int>(anyValue));
  if( acceptedTypes_.allowDouble() && anyValue.type() == typeid(double) )
    return any_cast<double>(anyValue);
  if( acceptedTypes_.allowString() && anyValue.type() == typeid(std::string) )
    return ::atof(any_cast<string>(anyValue).c_str());
  throwTypeError(entry,paramName,sublistName);
  return 0.0; // Will never get here!
}


std::string AnyNumberParameterEntryValidator::getString(
  const ParameterEntry &entry, const std::string &paramName
  ,const std::string &sublistName, const bool activeQuery
  ) const
{
  const any &anyValue = entry.getAny(activeQuery);
  if( acceptedTypes_.allowInt() && anyValue.type() == typeid(int) )
    return Utils::toString(any_cast<int>(anyValue));
  if( acceptedTypes_.allowDouble() && anyValue.type() == typeid(double) )
    return Utils::toString(any_cast<double>(anyValue));
  if( acceptedTypes_.allowString() && anyValue.type() == typeid(std::string) )
    return any_cast<string>(anyValue);
  throwTypeError(entry,paramName,sublistName);
  return ""; // Will never get here!
}


int AnyNumberParameterEntryValidator::getInt(
  ParameterList &paramList, const std::string &paramName
  ,const int defaultValue
  ) const
{
  const ParameterEntry *entry = paramList.getEntryPtr(paramName);
  if(entry) return getInt(*entry,paramName,paramList.name(),true);
  return paramList.get(paramName,defaultValue);
}


double AnyNumberParameterEntryValidator::getDouble(
  ParameterList &paramList, const std::string &paramName
  ,const double defaultValue
  ) const
{
  const ParameterEntry *entry = paramList.getEntryPtr(paramName);
  if(entry) return getDouble(*entry,paramName,paramList.name(),true);
  return paramList.get(paramName,defaultValue);
}


std::string AnyNumberParameterEntryValidator::getString(
  ParameterList &paramList, const std::string &paramName
  ,const std::string &defaultValue
  ) const
{
  const ParameterEntry *entry = paramList.getEntryPtr(paramName);
  if(entry) return getString(*entry,paramName,paramList.name(),true);
  return paramList.get(paramName,defaultValue);
}

  
// Overridden from ParameterEntryValidator


void AnyNumberParameterEntryValidator::printDoc(
  std::string         const& docString
  ,std::ostream            & out
  ) const
{
  StrUtils::printLines(out,"# ",docString);
  out << "#  Accepted types: " << acceptedTypesString_ << ".\n";
}


Teuchos::RefCountPtr<const Array<std::string> >
AnyNumberParameterEntryValidator::validStringValues() const
{
  return null;
}


void AnyNumberParameterEntryValidator::validate(
  ParameterEntry  const& entry
  ,std::string    const& paramName
  ,std::string    const& sublistName
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


void AnyNumberParameterEntryValidator::throwTypeError(
  ParameterEntry  const& entry
  ,std::string    const& paramName
  ,std::string    const& sublistName
  ) const
{
  TEST_FOR_EXCEPTION_PURE_MSG(
    true, Exceptions::InvalidParameterType
    ,"Error, the parameter {paramName=\""<<paramName<<"\""
    ",type=\""<<entry.getAny(false).typeName()<<"\"}"
    << "\nin the sublist \"" << sublistName << "\""
    << "\nhas the wrong type."
    << "\n\nThe accepted types are: " << acceptedTypesString_ << "!";
    );
}


} // namespace Teuchos


// Nonmmeber helper functions

Teuchos::RefCountPtr<Teuchos::AnyNumberParameterEntryValidator>
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
  paramList->set(
    paramName,value,docString,
    anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::PREFER_INT, acceptedTypes
      )
    );
}


void Teuchos::setDoubleParameter(
  std::string const& paramName,
  double const& value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  )
{
  TEST_FOR_EXCEPT(0==paramList);
  paramList->set(
    paramName,value,docString,
    anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::PREFER_DOUBLE, acceptedTypes
      )
    );
}


void Teuchos::setNumericStringParameter(
  std::string const& paramName,
  std::string const& value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  )
{
  TEST_FOR_EXCEPT(0==paramList);
  paramList->set(
    paramName,value,docString,
    anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::PREFER_STRING, acceptedTypes
      )
    );
}


int Teuchos::getIntParameter(
  ParameterList const& paramList,
  std::string const& paramName
  )
{
  const ParameterEntry &entry = paramList.getEntry(paramName);
  RefCountPtr<const AnyNumberParameterEntryValidator>
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
  RefCountPtr<const AnyNumberParameterEntryValidator>
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
  RefCountPtr<const AnyNumberParameterEntryValidator>
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


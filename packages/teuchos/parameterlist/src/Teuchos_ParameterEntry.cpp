// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterEntry.hpp" // class definition
#include "Teuchos_ParameterList.hpp"	 // for sublists
#include "Teuchos_TwoDArray.hpp"


namespace Teuchos {


ParameterEntry::ParameterEntry() :
  isUsed_(false),
  isDefault_(false)
{}


ParameterEntry::ParameterEntry(const ParameterEntry& source)
{
  operator=(source);
}


ParameterEntry& ParameterEntry::operator=(const ParameterEntry& source)
{
  if (&source == this)
    return *this;

  val_ = source.val_;
  isUsed_ = source.isUsed_;
  isDefault_ = source.isDefault_;
  docString_ = source.docString_;
  validator_ = source.validator_;

  return *this;
}

ParameterEntry::ParameterEntry(ParameterEntry&& other)
  : val_(std::move(other.val_)),
    isUsed_(other.isUsed_),
    isDefault_(other.isDefault_),
    docString_(std::move(other.docString_)),
    validator_(std::move(other.validator_))
{}

ParameterEntry& ParameterEntry::operator=(ParameterEntry&& other)
{
  if(this != &other)
  {
    this->val_ = std::move(other.val_);
    this->isUsed_ = other.isUsed_;
    this->isDefault_ = other.isDefault_;
    this->docString_ = std::move(other.docString_);
    this->validator_ = std::move(other.validator_);
  }
  return *this;
}

void ParameterEntry::setAnyValue(
  const any &value_in, bool isDefault_in
  )
{
  val_ = value_in;
  isDefault_ = isDefault_in;
  validator_ = null;
  isUsed_ = false;
  docString_ = "";
}


void ParameterEntry::setValidator(
  RCP<const ParameterEntryValidator> const& validator_in
  )
{
  validator_ = validator_in;
}


void ParameterEntry::setDocString(const std::string &docString_in)
{
  docString_ = docString_in;
}


ParameterList& ParameterEntry::setList(
  bool isDefault_in, const std::string &docString_in
  )
{
  val_ = ParameterList();
  isDefault_ = isDefault_in;
  isUsed_ = true;
  docString_ = docString_in;
  return any_cast<ParameterList>( val_ );
}


bool ParameterEntry::isList() const
{
  return ( !val_.has_value() ? false : val_.type() == typeid(ParameterList) );
}

std::ostream& ParameterEntry::leftshift(std::ostream& os, bool printFlags) const
{
  if( !this->isList() ) os << val_;

  if(printFlags) {
    if (isDefault_)
      os << "   [default]";
    else if (!isUsed_)
      os << "   [unused]";
  }

  return os;
}

bool ParameterEntry::isTwoDArray() const{
  std::string formatString = getTwoDArrayTypeNameTraitsFormat();
  size_t starPos = formatString.find("*");
  std::string prefix = formatString.substr(0,starPos);
  std::string postfix = formatString.substr(starPos+1);
  std::string valueTypeName = val_.typeName();
  size_t prePos = valueTypeName.find(prefix);
  size_t postPos = valueTypeName.find(postfix);
  return (prePos != std::string::npos) && (prePos==0)
    && (postPos != std::string::npos) && (prePos < postPos);
}

bool ParameterEntry::isArray() const{
  std::string formatString = getArrayTypeNameTraitsFormat();
  size_t starPos = formatString.find("*");
  std::string prefix = formatString.substr(0,starPos);
  std::string postfix = formatString.substr(starPos+1);
  std::string valueTypeName = val_.typeName();
  size_t prePos = valueTypeName.find(prefix);
  size_t postPos = valueTypeName.find(postfix);
  return (prePos != std::string::npos) && (prePos==0)
    && (postPos != std::string::npos) && (prePos < postPos);
}


// private


void ParameterEntry::reset()
{
  //delete val_;
  isUsed_ = false;
  isDefault_ = false;
}


} // namespace Teuchos



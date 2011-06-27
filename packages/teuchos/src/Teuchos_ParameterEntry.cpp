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
  return ( val_.empty() ? false : val_.type() == typeid(ParameterList) );
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
  return (prePos != std::string::npos) 
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
  return (prePos != std::string::npos) 
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



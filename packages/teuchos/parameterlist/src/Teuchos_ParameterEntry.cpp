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



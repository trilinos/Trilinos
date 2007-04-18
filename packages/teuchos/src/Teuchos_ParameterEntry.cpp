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
  const any &value, bool isDefault
  )
{
  val_ = value;
  isDefault_ = isDefault;
  validator_ = null;
  isUsed_ = false;
  docString_ = "";
}


void ParameterEntry::setValidator(
  RefCountPtr<const ParameterEntryValidator> const& validator
  )
{
  validator_ = validator;
}


void ParameterEntry::setDocString(const std::string &docString)
{
  docString_ = docString;
}


ParameterList& ParameterEntry::setList(bool isDefault, const std::string &docString)
{
  val_ = ParameterList();
  isDefault_ = isDefault;
  isUsed_ = true;
  docString_ = docString;
  return any_cast<ParameterList>( val_ );
}


bool ParameterEntry::isList() const
{
  return ( val_.empty() ? false : val_.type() == typeid(ParameterList) );
}


ostream& ParameterEntry::leftshift(ostream& os, bool printFlags) const
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


// private


void ParameterEntry::reset()
{
  //delete val_;
  isUsed_ = false;
  isDefault_ = false;
}


} // namespace Teuchos

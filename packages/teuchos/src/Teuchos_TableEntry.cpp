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

#include "Teuchos_TableEntry.hpp"

using namespace Teuchos;


/* --------- base class methods ---------------------------------------------- */

string TableEntry::toChoppedString(int maxWidth) const
{
  return toString().substr(0, maxWidth);
}



/* --------- DoubleEntry methods -------------------------------------------- */

DoubleEntry::DoubleEntry(const double& value, int precision)
  : TableEntry(), data_(value), precision_(precision)
{}

string DoubleEntry::toString() const 
{
  TeuchosOStringStream toss;
  toss << std::setprecision(precision_) << data_;
  return toss.str();
}



/* --------- IntEntry methods -------------------------------------------- */

IntEntry::IntEntry(int value)
  : TableEntry(), data_(value)
{}

string IntEntry::toString() const 
{
  TeuchosOStringStream toss;
  toss << data_;
  return toss.str();
}



/* --------- StringEntry methods -------------------------------------------- */

StringEntry::StringEntry(string value)
  : TableEntry(), data_(value)
{}

string StringEntry::toString() const 
{
  return data_;
}





/* --------- CompoundEntryWithParentheses methods ------------------------- */

CompoundEntryWithParentheses
::CompoundEntryWithParentheses(const RefCountPtr<TableEntry>& first,
                                const RefCountPtr<TableEntry>& second,
                                bool spaceBeforeParens)
  : TableEntry(), 
    first_(first), 
    second_(second), 
    spaceBeforeParens_(spaceBeforeParens)
{}

string CompoundEntryWithParentheses::toString() const 
{
  TeuchosOStringStream toss;
  
  toss << first_->toString();
  if (spaceBeforeParens_) toss << " ";
  toss << "(" << second_->toString() << ")";
  
  return toss.str();
}







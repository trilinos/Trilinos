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

#include "Teuchos_ParameterList.hpp"	// class definition

/* NOTE: ASCI Red (TFLOP) does not support the i-> function for iterators 
 * in the STL.  Therefore when compiling for the TFLOP we must redefine the 
 * iterator from i-> to (*i). This slows things down on other platforms 
 * so we switch between the two when necessary.
 */
using namespace Teuchos;

ParameterList::ParameterList() {}

ParameterList::ParameterList(const ParameterList& source) 
{
  params_ = source.params_;
}

ParameterList& ParameterList::operator=(const ParameterList& source) 
{
  if (&source == this)
    return *this;

  params_ = source.params_;
  return *this;
}

ParameterList::~ParameterList() 
{
}

void ParameterList::unused(ostream& os) const
{
  for (ConstIterator i = params_.begin(); i != params_.end(); ++i) {
    if (!(entry(i).isUsed())) {
      os << "WARNING: Parameter \"" << name(i) << "\" " << entry(i)
	   << " is unused" << endl;
    }
  }
}

bool ParameterList::isSublist(const string& name) const
{
  ConstIterator i = params_.find(name);

  if (i != params_.end())
    return (entry(i).isList());

  return false;
}

bool ParameterList::isParameter(const string& name) const
{
  return (params_.find(name) != params_.end());
}

ParameterList& ParameterList::sublist(const string& name)
{
  // Find name in list, if it exists.
  Iterator i = params_.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params_.end()) {
     TEST_FOR_EXCEPTION( !entry(i).isList(), std::runtime_error,
	" Parameter " << name << " is not a list!" );
     return getValue<ParameterList>(entry(i));
  }

  // If it does not exist, create a new empty list and return a reference
  return params_[name].setList(true);
}

const ParameterList& ParameterList::sublist(const string& name) const
{
  // Find name in list, if it exists.
  ConstIterator i = params_.find(name);

  // If it does not exist, throw an error
  TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
	" Parameter " << name << " is not a valid list!" );

  // If it does exist and is a list, return the list value.
  TEST_FOR_EXCEPTION( !entry(i).isList(), std::runtime_error,
	" Parameter " << name << " is not a list!" );
  return getValue<ParameterList>(entry(i));
}
  
ostream& ParameterList::print(ostream& os, int indent) const
{
  if (params_.begin() == params_.end()) 
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    os << "[empty list]" << endl;
  }
  else 
    for (ConstIterator i = params_.begin(); i != params_.end(); ++i) 
    {
      for (int j = 0; j < indent; j ++)
	os << ' ';
      if (entry(i).isList()) 
      {
	os << name(i) << " -> " << endl;
	getValue<ParameterList>(entry(i)).print(os, indent + 2);
      }
      else
	os << name(i) << " = " << entry(i) << endl;
    }
  return os;
}


#if defined(TFLOP)

const string& ParameterList::name(ConstIterator i) const
{
  return ((*i).first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return ((*i).second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return ((*i).second);
}

#else

const string& ParameterList::name(ConstIterator i) const
{
  return (i->first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return (i->second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return (i->second);
}

#endif


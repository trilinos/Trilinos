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

#include "Teuchos_VerboseObject1.hpp"
#include "Teuchos_OSTab1.hpp"

namespace Teuchos {

// Private static data members

RefCountPtr<std::ostream> VerboseObject::defaultOStream_ = Teuchos::rcp(&std::cout,false);
EVerbosityLevel VerboseObject::defaultVerbLevel_ = VERB_DEFAULT;
int VerboseObject::defaultTabTag_ = 0;

// Public static member functions

void VerboseObject::setDefaultOStream( const RefCountPtr<std::ostream> &defaultOStream )
{
  defaultOStream_ = defaultOStream;
}

RefCountPtr<std::ostream>
VerboseObject::getDefaultOStream()
{
  return defaultOStream_;
}

void VerboseObject::setDefaultVerbLevel( const EVerbosityLevel defaultVerbLevel)
{
  defaultVerbLevel_ = defaultVerbLevel;
}

EVerbosityLevel VerboseObject::getDefaultVerbLevel()
{
  return defaultVerbLevel_;
}

void VerboseObject::setDefaultTabTag( const int defaultTabTag)
{
  defaultTabTag_ = defaultTabTag;
}

int VerboseObject::getDefaultTabTag()
{
  return defaultTabTag_;
}

// Constructors/Initializers

VerboseObject::VerboseObject(
  const EVerbosityLevel              verbLevel
  ,const RefCountPtr<std::ostream>   &oStream
  ,const int                         tabTag
  )
{
  this->initializeVerboseObject(verbLevel,oStream,tabTag);
}

void VerboseObject::initializeVerboseObject(
  const EVerbosityLevel              verbLevel
  ,const RefCountPtr<std::ostream>   &oStream
  ,const int                         tabTag
  )
{
  thisVerbLevel_ = verbLevel;
  thisOStream_ = oStream;
  thisTabTag_ = tabTag;
}

VerboseObject& VerboseObject::setOStream(const RefCountPtr<std::ostream> &oStream)
{
  thisOStream_ = oStream;
  return *this;
}

VerboseObject& VerboseObject::setVerbLevel(const EVerbosityLevel verbLevel)
{
  thisVerbLevel_ = verbLevel;
  return *this;
}

VerboseObject& VerboseObject::setTabTag(const int tabTag)
{
  thisTabTag_ = tabTag;
  return *this;
}

// Query functions

RefCountPtr<std::ostream>
VerboseObject::getOStream() const
{
  if(!thisOStream_.get())
    return defaultOStream_;
  return thisOStream_;
}

EVerbosityLevel VerboseObject::getVerbLevel() const
{
  if(thisVerbLevel_ == VERB_DEFAULT)
    return defaultVerbLevel_;
  return thisVerbLevel_;
}

int VerboseObject::getTabTag() const
{
  if(thisTabTag_ < 0)
    return defaultTabTag_;
  return thisTabTag_;
}

// Utility functions

OSTab VerboseObject::getOSTab() const
{
  return OSTab(this->getTabTag());
}

} // namespace Teuchos

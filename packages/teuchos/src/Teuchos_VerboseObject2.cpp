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

#include "Teuchos_VerboseObject2.hpp"

namespace Teuchos {

// Private static data members

RefCountPtr<FancyOStream>
VerboseObject::defaultOStream_ = rcp(new FancyOStream(rcp(&std::cout,false)));

EVerbosityLevel
VerboseObject::defaultVerbLevel_ = VERB_DEFAULT;

// Public static member functions

void VerboseObject::setDefaultOStream( const RefCountPtr<FancyOStream> &defaultOStream )
{
  defaultOStream_ = defaultOStream;
}

RefCountPtr<FancyOStream>
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

// Constructors/Initializers

VerboseObject::VerboseObject(
  const EVerbosityLevel              verbLevel
  ,const RefCountPtr<FancyOStream>   &oStream
  )
{
  this->initializeVerboseObject(verbLevel,oStream);
}

void VerboseObject::initializeVerboseObject(
  const EVerbosityLevel              verbLevel
  ,const RefCountPtr<FancyOStream>   &oStream
  )
{
  thisVerbLevel_ = verbLevel;
  thisOStream_ = oStream;
}

VerboseObject& VerboseObject::setOStream(const RefCountPtr<FancyOStream> &oStream)
{
  thisOStream_ = oStream;
  return *this;
}

VerboseObject& VerboseObject::setVerbLevel(const EVerbosityLevel verbLevel)
{
  thisVerbLevel_ = verbLevel;
  return *this;
}

VerboseObject& VerboseObject::setLinePrefix(const std::string &linePrefix)
{
  thisLinePrefix_ = linePrefix;
  return *this;
}

// Query functions

RefCountPtr<FancyOStream>
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

// Utility functions

OSTab VerboseObject::getOSTab(const int tabs,const std::string &linePrefix) const
{
  return OSTab(this->getOStream(),tabs,linePrefix.length() ? linePrefix : thisLinePrefix_);
}

} // namespace Teuchos

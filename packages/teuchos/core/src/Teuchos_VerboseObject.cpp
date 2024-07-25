// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"


namespace Teuchos {


// Private static data members


RCP<FancyOStream>& VerboseObjectBase::privateDefaultOStream()
{
  static RCP<FancyOStream> defaultOStream;
  if (defaultOStream.get()==NULL) {
    defaultOStream = fancyOStream(rcpFromRef(std::cout));
    defaultOStream->setOutputToRootOnly(0);
//    if(GlobalMPISession::getNProc()>1)
//      defaultOStream->setShowProcRank(true);
  }
  return defaultOStream;
}


// Public static member functions


void VerboseObjectBase::setDefaultOStream(
  const RCP<FancyOStream> &defaultOStream
  )
{
  privateDefaultOStream() = defaultOStream;
}


RCP<FancyOStream>
VerboseObjectBase::getDefaultOStream()
{
  return privateDefaultOStream();
}

// Destructor

VerboseObjectBase::~VerboseObjectBase()
{
}

// Constructors/Initializers


VerboseObjectBase::VerboseObjectBase(
  const RCP<FancyOStream>   &oStream
  )
  : thisOverridingOStream_(null)
{
  this->initializeVerboseObjectBase(oStream);
}


void VerboseObjectBase::initializeVerboseObjectBase(
  const RCP<FancyOStream>   &oStream
  )
{
  thisOStream_ = oStream;
}


const VerboseObjectBase&
VerboseObjectBase::setOStream(const RCP<FancyOStream> &oStream) const
{
  thisOStream_ = oStream;
  informUpdatedVerbosityState();
  return *this;
}


const VerboseObjectBase&
VerboseObjectBase::setOverridingOStream(
  const RCP<FancyOStream> &oStream
  ) const
{
  thisOverridingOStream_ = oStream;
  informUpdatedVerbosityState();
  return *this;
}


VerboseObjectBase&
VerboseObjectBase::setLinePrefix(const std::string &linePrefix)
{
  thisLinePrefix_ = linePrefix;
  informUpdatedVerbosityState();
  return *this;
}


// Query functions


RCP<FancyOStream>
VerboseObjectBase::getOStream() const
{
  if(!is_null(thisOverridingOStream_))
    return thisOverridingOStream_;
  if(is_null(thisOStream_))
    return getDefaultOStream();
  return thisOStream_;
}


RCP<FancyOStream>
VerboseObjectBase::getOverridingOStream() const
{
  return thisOverridingOStream_;
}


std::string VerboseObjectBase::getLinePrefix() const
{
  return thisLinePrefix_;
}


// Utility functions


OSTab VerboseObjectBase::getOSTab(
  const int tabs,const std::string &linePrefix
  ) const
{
  return OSTab(
    this->getOStream(), tabs, linePrefix.length()
    ? linePrefix : this->getLinePrefix()
    );
}


// protected


void VerboseObjectBase::informUpdatedVerbosityState() const
{}


} // namespace Teuchos

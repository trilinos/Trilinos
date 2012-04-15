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
  informUpdatedVerbosityState();
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

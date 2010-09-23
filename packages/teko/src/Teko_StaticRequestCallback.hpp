/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#ifndef __Teko_StaticRequestCallback_hpp__
#define __Teko_StaticRequestCallback_hpp__

#include "Teko_RequestCallback.hpp"

namespace Teko {

/** A simple request interface that takes
  * a static bit of data to return to Teko.
  * This is meant primarily as a testing and 
  * early stage development tool. 
  *
  * The constructor takes an object of the same
  * type the class was templated on. It also 
  * takes a string to be used to match the 
  * request.
  */
template <typename DataT>
class StaticRequestCallback : public RequestCallback<DataT> {
public:
   StaticRequestCallback(const std::string & name,DataT data) : name_(name), data_(data) {}

   DataT request(const RequestMesg & rm);
   void preRequest(const RequestMesg & rm);
   bool handlesRequest(const RequestMesg & rm);

private:
   std::string name_;
   DataT data_;
};

template <typename DataT>
DataT StaticRequestCallback<DataT>::request(const RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm));
 
   return data_;   
}

template <typename DataT>
void StaticRequestCallback<DataT>::preRequest(const RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm));
}

template <typename DataT>
bool StaticRequestCallback<DataT>::handlesRequest(const RequestMesg & rm)
{
   return name_==rm.getName();
}

}

#endif

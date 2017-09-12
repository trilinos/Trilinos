// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_ResponseBase_hpp__
#define __Panzer_ResponseBase_hpp__

#include <string>

#include "Panzer_GlobalEvaluationData.hpp"

namespace panzer {

// *** Derive, do not modify this interface ***

/** This object serves as the base class which
  * the user derives from to build there response.
  * Note that by itself the interface for this is 
  * not useful for storing data.
  */
class ResponseBase : public GlobalEvaluationData_Default {
public:

   /** Only available constructor for this object.
     * This gurantees that the reponse object is
     * instantiatied with a name, and gives us basic
     * capability for name lookups.
     */
   ResponseBase(const std::string & responseName) :
    responseName_(responseName) {}

   virtual ~ResponseBase() {}

   /** Get the unmodified name for this response.
     */
   std::string getName() const { return responseName_; }

   /** Get the name of this response useful
     * for looking up data containers.
     */
   std::string getLookupName() const { return buildLookupName(responseName_); }

   /** Static member to build consistent look up names
     * based on a response name.
     */
   static std::string buildLookupName(const std::string & responseName)
   { return "RESPONSE_"+responseName; }

   //! Inherited from GlobalEvaluationData, 
   virtual void ghostToGlobal(int) 
   { scatterResponse(); }

   virtual void initializeData()
   { initializeResponse(); }

   //! Prepare the response for access by the user (do global communication)
   virtual void scatterResponse() = 0;

   virtual void initializeResponse() = 0;

private:

   std::string responseName_; 
 
   // hide these methods
   ResponseBase();
   ResponseBase(const ResponseBase &);
};

}

#endif

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

#ifndef __Teko_RequestMesg_hpp__
#define __Teko_RequestMesg_hpp__

#include <string>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Teko {

class RequestMesg {
public:
   /** Construct a request messages specifing the
     * details of the request.
     *
     * \param[in] name Name of request to be satisfied.
     * \param[in] tag Optional tag describing other information
     *                about this tag.
     */
   explicit RequestMesg(const std::string & name,unsigned int tag=0)
      : name_(name), tag_(tag) {}

   /** Construct a parameter list message. This sets the 
     * name to "Parameter List" and the tag to 0. But now the
     * parameter list can completely describe the request.
     *
     * \param[in] pl Parameter list describing the request
     */
   explicit RequestMesg(const Teuchos::RCP<const Teuchos::ParameterList> & pl)
      : paramList_(pl) 
   {
      fromParameterList(*paramList_);
   }

   //! Simple base class destructor
   virtual ~RequestMesg() {}

   //! Get the name for this request
   std::string getName() const
   { return name_; }

   //! Get the tag for this request
   unsigned int getTag() const
   { return tag_; }

   //! Get parameter list for this request
   const Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const
   { return paramList_.getConst(); }
   
protected:

   void fromParameterList(const Teuchos::ParameterList & pl)
   {
      name_ = "Parameter List";
      tag_ = 0;
      if(pl.isParameter("Name"))
         name_ = pl.get<std::string>("Name");
      if(pl.isParameter("Tag"))
         tag_ = pl.get<unsigned int>("Tag");
   }

   std::string name_;
   unsigned int tag_;
   Teuchos::RCP<const Teuchos::ParameterList> paramList_;
};

// simple stream interface for RequestMesg
inline std::ostream & operator<<(std::ostream & os,const Teko::RequestMesg & rm)
{
   os << "RequestMesg <"
      << "name = \"" << rm.getName() << "\", "
      << "tag = " << rm.getTag() << ">";
 
   return os;
}

} // end namespace Teko

#endif

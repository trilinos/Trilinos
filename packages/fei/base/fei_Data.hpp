/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_Data_hpp_
#define _fei_Data_hpp_

#include <string>

/**
  This is a very simple class for passing stuff around
  in a void pointer. It has the ability to store and query
  a type name, so at least there can be user-enforced
  type safety.

  When setTypeName is called, a char* is created and a copy
  of the input argument is taken. This char* is later destroyed
  by the Data destructor. The void* dataPtr_ member is not
  destroyed, it is just a copy of a pointer.
*/

class Data {
 public:
  /** Default constructor. */
   Data() : typeName_(), dataPtr_(NULL) {}

   /** Default destructor. */
   virtual ~Data() {}

   /** Set a string representing the type of the object stored in
       'getDataPtr()'. */
   void setTypeName(const char* name) { typeName_ = name;}

   /** Query the string representing the type of the object stored in
       'getDataPtr()'. */
   const char* getTypeName() const {return(typeName_.c_str());}

   /** Set the contents of the data pointer. */
   void setDataPtr(void* ptr) {dataPtr_ = ptr;}

  /** Retrieve the contents of the data pointer. */
   void* getDataPtr() const {return(dataPtr_);}

 private:
   std::string typeName_;
   void* dataPtr_;
};

#endif


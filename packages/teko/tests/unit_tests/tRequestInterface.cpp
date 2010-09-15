/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"

#define SS_ECHO(ops) { std::stringstream ss; ss << ops; ECHO(ss.str()) };

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

// Build PreconditionerFactory testing rig
///////////////////////////////////////////////////////////

class TestFactory : public Teko::PreconditionerFactory {
public:
   Teko::LinearOp buildPreconditionerOperator(Teko::LinearOp & lo,
                                              Teko::PreconditionerState & state) const;
   
   mutable double timestep_;
   mutable Teko::LinearOp pcdOp_;
   mutable std::string string_;
};

Teko::LinearOp TestFactory::buildPreconditionerOperator(Teko::LinearOp & lo,
                                                        Teko::PreconditionerState & state) const
{
   Teko::RequestMesg pcdMesg("PCD Op"), tsMesg("timestep"), strMesg("name");

   pcdOp_ = callbackHandler_->request<Teko::LinearOp>(pcdMesg);
   timestep_ = callbackHandler_->request<double>(tsMesg);
   string_ = callbackHandler_->request<std::string>(strMesg);

   return Teuchos::null;
}

// test rigs callbacks
///////////////////////////////////////////////////////////

class PCDCallback : public Teko::RequestCallback<Teko::LinearOp> {
public:
   Teko::LinearOp request(const Teko::RequestMesg & rd)
   {
      TEUCHOS_ASSERT(handlesRequest(rd));
      return Teuchos::null;  
   }

   bool handlesRequest(const Teko::RequestMesg & rd)
   {
      if(rd.getName()=="PCD Op")
         return true;
      return false;
   }

   void preRequest(const Teko::RequestMesg & rm) { }
};

class TSCallback : public Teko::RequestCallback<double> {
public:
   double request(const Teko::RequestMesg & rd)
   {
      TEUCHOS_ASSERT(handlesRequest(rd));
      return 0.1;  
   }

   bool handlesRequest(const Teko::RequestMesg & rd)
   {
      if(rd.getName()=="timestep")
         return true;
      return false;
   }

   void preRequest(const Teko::RequestMesg & rm) { }
};

class StringCallback : public Teko::RequestCallback<std::string> {
public:
   std::string request(const Teko::RequestMesg & rd)
   {
      TEUCHOS_ASSERT(handlesRequest(rd));
      return "the string";
   }

   bool handlesRequest(const Teko::RequestMesg & rd)
   {
      if(rd.getName()=="name")
         return true;
      return false;
   }

   void preRequest(const Teko::RequestMesg & rm) { }
};

///////////////////////////////////////////////////////////

TEUCHOS_UNIT_TEST(tRequestInterface, test_request_interface)
{
   RCP<TestFactory> precFact = rcp(new TestFactory);
   RCP<Teko::PreconditionerState> state = precFact->buildPreconditionerState();

   RCP<Teko::RequestHandler> rh = rcp(new Teko::RequestHandler);

   rh->addRequestCallback(rcp(new PCDCallback));
   rh->addRequestCallback(rcp(new TSCallback));
   rh->addRequestCallback(rcp(new StringCallback));

   {
      Teko::RequestMesg pcdMesg("PCD Op"), tsMesg("timestep"), strMesg("name");
      rh->preRequest<Teko::LinearOp>(pcdMesg);
      rh->preRequest<double>(tsMesg);
      rh->preRequest<std::string>(strMesg);
   }
   precFact->setRequestHandler(rh);

   Teko::LinearOp A;
   Teko::LinearOp result = precFact->buildPreconditionerOperator(A,*state);

   TEST_ASSERT(precFact->timestep_==0.1);
   TEST_ASSERT(precFact->pcdOp_==Teuchos::null);
   TEST_ASSERT(precFact->string_=="the string");

   try {
      Teko::RequestMesg intMesg("Int Op");
      rh->preRequest<int>(intMesg);
      out << "Found <int> with name \"Int Op\" in preRequest" << std::endl;
      TEST_ASSERT(false);
   } catch(std::exception & e) 
   { out << "expected exception = " << e.what() << std::endl; }

   try {
      Teko::RequestMesg intMesg("Int Op");
      int size = rh->request<int>(intMesg);
      out << "Found <int> with name \"Int Op\" value=" << size << std::endl;
      TEST_ASSERT(false);
   } catch(std::exception & e) 
   { out << "expected exception = " << e.what() << std::endl; }

   try {
      Teko::RequestMesg intMesg("PCD Op");
      int lo = rh->request<int>(intMesg);
      out << "Found <int> with name \"PCD Op\" value=" << lo << std::endl;
      TEST_ASSERT(false);
   } catch(std::exception & e) 
   { out << "expected exception = " << e.what() << std::endl; }
}

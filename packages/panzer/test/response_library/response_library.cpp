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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_StrUtils.hpp"

#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseContainer.hpp"
#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_ResponseUtilities.hpp"
#include "Panzer_RLDynamicDispatch.hpp"
#include "Panzer_WorksetContainer.hpp"

#include "TestEvaluators.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <vector>
#include <map>

namespace panzer {

TEUCHOS_UNIT_TEST(response_library, test)
{
  typedef Traits::Residual EvalT;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > comm 
           = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif
 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  int worksetSize = 10;
  Teuchos::ParameterList mainParam;
  mainParam.set<int>("Workset Size", worksetSize);

  // build basic field manager
  PHX::FieldManager<Traits> fm;
  {
     RCP<PHX::Evaluator<Traits> > testEval 
        = rcp(new TestEvaluator<EvalT,Traits>(mainParam));
     fm.registerEvaluator<EvalT>(testEval);
  }

  ResponseId dResp  = buildResponse("Dog","Functional");
  ResponseId hResp  = buildResponse("Horse","Functional");
  RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>(Teuchos::rcp(new panzer::WorksetContainer),Teuchos::null,Teuchos::null));
  rLibrary->defineDefaultAggregators();

  rLibrary->reserveVolumeResponse<EvalT>(dResp,"block_0");
  rLibrary->reserveVolumeResponse<EvalT>(hResp,"block_1");

  // test uninitialized access
  TEST_THROW(rLibrary->getVolumeResponse(dResp,"block_0"),std::logic_error);

  std::vector<std::string> eBlocks;
  rLibrary->getRequiredElementBlocks(eBlocks);
  std::sort(eBlocks.begin(),eBlocks.end());
  TEST_EQUALITY(eBlocks.size(),2);
  TEST_EQUALITY(eBlocks[0],"block_0");
  TEST_EQUALITY(eBlocks[1],"block_1");
}

TEUCHOS_UNIT_TEST(response_library, pl_reader)
{
   Teuchos::ParameterList pl0;
   pl0.sublist("Cat");
      pl0.sublist("Cat").set("Type","None");
      pl0.sublist("Cat").set("Field Name","Cat");
      pl0.sublist("Cat").set("Element Blocks","block_0");
      pl0.sublist("Cat").set("Evaluation Types","Residual,Jacobian");
   pl0.sublist("HorseFunc");
      pl0.sublist("HorseFunc").set("Type","Functional");
      pl0.sublist("HorseFunc").set("Field Name","Horse");
      pl0.sublist("HorseFunc").set("Element Blocks","block_0");
      pl0.sublist("HorseFunc").set("Evaluation Types","Residual");
   pl0.sublist("HorseMax");
      pl0.sublist("HorseMax").set("Type","Maximum");
      pl0.sublist("HorseMax").set("Field Name","Horse");
      pl0.sublist("HorseMax").set("Element Blocks","block_0,block_4");
      pl0.sublist("HorseMax").set("Evaluation Types","Residual,Jacobian");

   Teuchos::ParameterList pl1;
   pl1.sublist("Cat");
      pl1.sublist("Cat").set("Type","Functional");
      pl1.sublist("Cat").set("Field Name","Cat");
      pl1.sublist("Cat").set("Element Blocks"," ");
      pl1.sublist("Cat").set("Evaluation Types"," ");
   pl1.sublist("HorseFunc");
      pl1.sublist("HorseFunc").set("Type","Functional");
      pl1.sublist("HorseFunc").set("Field Name","Horse");
      pl1.sublist("HorseFunc").set("Element Blocks","block_1");
      pl1.sublist("HorseFunc").set("Evaluation Types","Residual");
   pl1.sublist("HorseMax");
      pl1.sublist("HorseMax").set("Type","Maximum");
      pl1.sublist("HorseMax").set("Field Name","Horse");
      pl1.sublist("HorseMax").set("Element Blocks","block_1");
      pl1.sublist("HorseMax").set("Evaluation Types","Residual,Jacobian");

   std::map<std::string,std::pair<ResponseId,std::pair<std::list<std::string>,std::list<std::string> > > > responses;
   TEST_THROW(buildResponseMap(pl1,responses),Teuchos::Exceptions::InvalidParameterValue);
   TEST_NOTHROW(buildResponseMap(pl0,responses));

   TEST_EQUALITY(responses.size(),3);
   TEST_ASSERT(responses.find("Cat")!=responses.end());
   TEST_ASSERT(responses.find("HorseFunc")!=responses.end());
   TEST_ASSERT(responses.find("HorseMax")!=responses.end());

   out << responses["Cat"].first.getString() << std::endl;
   out << responses["HorseMax"].first.getString() << std::endl;
   out << responses["HorseFunc"].first.getString() << std::endl;

   TEST_ASSERT(responses["HorseMax"].first==ResponseId("Horse","Maximum"));
   TEST_EQUALITY(responses["HorseMax"].second.first.size(),2);
   TEST_EQUALITY(responses["HorseMax"].second.second.size(),2);

   TEST_ASSERT(responses["HorseFunc"].first==ResponseId("Horse","Functional"));
   TEST_EQUALITY(responses["HorseFunc"].second.first.size(),1);
   TEST_EQUALITY(responses["HorseFunc"].second.second.size(),1);

   TEST_ASSERT(responses["Cat"].first==ResponseId("Cat","None"));
   TEST_EQUALITY(responses["Cat"].second.first.size(),1);
   TEST_EQUALITY(responses["Cat"].second.second.size(),2);
   

   std::vector<std::string> tokens;
   CommaSeperatedEntryValidator::split(" None :Residual, Jacobian",":",tokens);
    
   TEST_EQUALITY(tokens.size(),2);
   TEST_EQUALITY(tokens[0],"None");
   TEST_EQUALITY(tokens[1],"Residual, Jacobian");

   std::string next = tokens[1];
   CommaSeperatedEntryValidator::split(next,",",tokens);

   TEST_EQUALITY(tokens.size(),2);
   TEST_EQUALITY(tokens[0],"Residual");
   TEST_EQUALITY(tokens[1],"Jacobian");

   CommaSeperatedEntryValidator rev;

   // standard operating
   {
      Teuchos::ParameterEntry entry;

      entry.setValue<std::string>("Residual, Jacobian");
      TEST_NOTHROW(rev.validate(entry,"FuncValue","ResponseList"));

      entry.setValue<std::string>("Residual");
      TEST_NOTHROW(rev.validate(entry,"FuncValue","ResponseList"));
   }

   // failing tests
   {
      Teuchos::ParameterEntry entry;

      entry.setValue<std::string>(" ");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);
   }
}

TEUCHOS_UNIT_TEST(response_library, dyn_dispatch)
{
  typedef Traits::Residual EvalT;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > comm 
           = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif
 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  int worksetSize = 10;
  Teuchos::ParameterList mainParam;
  mainParam.set<int>("Workset Size", worksetSize);

  // build basic field manager
  PHX::FieldManager<Traits> fm;
  {
     RCP<PHX::Evaluator<Traits> > testEval 
        = rcp(new TestEvaluator<EvalT,Traits>(mainParam));
     fm.registerEvaluator<EvalT>(testEval);
  }

  ResponseId dResp  = buildResponse("Dog","Functional");
  ResponseId hResp  = buildResponse("Horse","Functional");
  RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>(Teuchos::rcp(new panzer::WorksetContainer),Teuchos::null,Teuchos::null));
  rLibrary->defineDefaultAggregators();

  rLibrary->reserveVolumeResponse(dResp,"block_0","Residual");
  rLibrary->reserveVolumeResponse(hResp,"block_1","Residual");

  // test uninitialized access
  TEST_THROW(rLibrary->getVolumeResponse(dResp,"block_0"),std::logic_error);

  std::vector<std::string> eBlocks;
  rLibrary->getRequiredElementBlocks(eBlocks);
  std::sort(eBlocks.begin(),eBlocks.end());
  TEST_EQUALITY(eBlocks.size(),2);
  TEST_EQUALITY(eBlocks[0],"block_0");
  TEST_EQUALITY(eBlocks[1],"block_1");
}

#ifdef HAVE_STOKHOS
TEUCHOS_UNIT_TEST(response_library, sg_dyn_dispatch)
{
  typedef Traits::Residual EvalT;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > comm 
           = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif
 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  int worksetSize = 10;
  Teuchos::ParameterList mainParam;
  mainParam.set<int>("Workset Size", worksetSize);

  // build basic field manager
  PHX::FieldManager<Traits> fm;
  {
     RCP<PHX::Evaluator<Traits> > testEval 
        = rcp(new TestEvaluator<EvalT,Traits>(mainParam));
     fm.registerEvaluator<EvalT>(testEval);
  }

  ResponseId dResp  = buildResponse("Dog","Functional");
  ResponseId hResp  = buildResponse("Horse","Functional");
  RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>(Teuchos::rcp(new panzer::WorksetContainer),Teuchos::null,Teuchos::null));
  rLibrary->defineDefaultAggregators();

  rLibrary->reserveVolumeResponse(dResp,"block_0","Residual");
  rLibrary->reserveVolumeResponse(hResp,"block_1","SGResidual");

  // test uninitialized access
  TEST_THROW(rLibrary->getVolumeResponse(dResp,"block_0"),std::logic_error);

  std::vector<std::string> eBlocks;
  rLibrary->getRequiredElementBlocks(eBlocks);
  std::sort(eBlocks.begin(),eBlocks.end());
  TEST_EQUALITY(eBlocks.size(),2);
  TEST_EQUALITY(eBlocks[0],"block_0");
  TEST_EQUALITY(eBlocks[1],"block_1");
}
#endif

}

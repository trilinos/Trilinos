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
        = Teuchos::rcp(new ResponseLibrary<Traits>());
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
   Teuchos::ParameterList pl;
   Teuchos::ParameterList & pl0 = pl.sublist("block_0");
   pl0.set("Cat","None(Cat):Residual,Jacobian");
   pl0.set("HorseFunc","Functional(Horse):Residual");
   pl0.set("HorseMax","Maximum(Horse):Residual,Jacobian");
   Teuchos::ParameterList & pl1 = pl.sublist("block_1");
   pl1.set("Cat","None(Dog):   ");
   pl1.set("Horse","    :Residual,Jacobian");
   pl1.set("Monster","Maximum");

   std::map<std::string,std::pair<ResponseId,std::set<std::string> > > responses;
   TEST_THROW(buildResponseMap(pl.sublist("block_1"),responses),Teuchos::Exceptions::InvalidParameterValue);
   TEST_NOTHROW(buildResponseMap(pl.sublist("block_0"),responses));

   TEST_EQUALITY(responses.size(),3);
   TEST_ASSERT(responses.find("Cat")!=responses.end());
   TEST_ASSERT(responses.find("HorseFunc")!=responses.end());
   TEST_ASSERT(responses.find("HorseMax")!=responses.end());

   out << responses["Cat"].first.getString() << std::endl;
   out << responses["HorseMax"].first.getString() << std::endl;
   out << responses["HorseFunc"].first.getString() << std::endl;

   TEST_ASSERT(responses["HorseMax"].first==ResponseId("Horse","Maximum"));
   TEST_EQUALITY(responses["HorseMax"].second.size(),2);

   TEST_ASSERT(responses["HorseFunc"].first==ResponseId("Horse","Functional"));
   TEST_EQUALITY(responses["HorseFunc"].second.size(),1);

   TEST_ASSERT(responses["Cat"].first==ResponseId("Cat","None"));
   TEST_EQUALITY(responses["Cat"].second.size(),2);
   

   std::vector<std::string> tokens;
   ResponseEntryValidator::split(" None :Residual, Jacobian",":",tokens);
    
   TEST_EQUALITY(tokens.size(),2);
   TEST_EQUALITY(tokens[0],"None");
   TEST_EQUALITY(tokens[1],"Residual, Jacobian");

   std::string next = tokens[1];
   ResponseEntryValidator::split(next,",",tokens);

   TEST_EQUALITY(tokens.size(),2);
   TEST_EQUALITY(tokens[0],"Residual");
   TEST_EQUALITY(tokens[1],"Jacobian");

   ResponseEntryValidator rev;

   // standard operating
   {
      Teuchos::ParameterEntry entry;

      entry.setValue<std::string>("Functional(Horse) : Residual, Jacobian");
      TEST_NOTHROW(rev.validate(entry,"FuncValue","ResponseList"));

      entry.setValue<std::string>("Functional (Monkey) : Residual");
      TEST_NOTHROW(rev.validate(entry,"FuncValue","ResponseList"));
   }

   // failing tests
   {
      Teuchos::ParameterEntry entry;

      entry.setValue<std::string>("Functional : Residual");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);

      entry.setValue<std::string>("Functional () : Residual");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);

      entry.setValue<std::string>("Functional (Monkey) : ");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);

      entry.setValue<std::string>(" : Residual");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);

      entry.setValue<std::string>(" Residual");
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
        = Teuchos::rcp(new ResponseLibrary<Traits>());
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


}

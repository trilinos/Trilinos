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

#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseContainer.hpp"
#include "Panzer_ResponseLibrary.hpp"

#include "TestEvaluators.hpp"


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

}

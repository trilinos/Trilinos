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
#include "Panzer_ResponseContainer_novel.hpp"

#include "TestEvaluators.hpp"


namespace panzer {
namespace novel {

TEUCHOS_UNIT_TEST(response_assembly, test)
{
  typedef Traits::Residual EvalT;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
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

  std::vector<std::string> fields, testFields;
  fields.push_back("Dog");
  fields.push_back("Horse");

  RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>());
  RCP<ResponseContainerBase<Traits> > container
        = Teuchos::rcp(new ResponseContainer<EvalT,Traits>(rLibrary));

  ResponseId dResp  = buildResponse("Dog","Functional");
  ResponseId hResp  = buildResponse("Horse","Functional");
  ResponseId hResp2 = buildResponse("Horse","Max");       // should be able to add
  ResponseId cResp  = buildResponse("Cat","Functional"); // not added
  container->reserve(dResp);
  container->reserve(hResp);
  container->reserve(hResp2);
  container->reserve(hResp); // You can reserve things multiple times, should
                            // not increased number of reserved items

  TEST_EQUALITY(container->getReserved().size(),3);
  TEST_ASSERT(container->contains(dResp));
  TEST_ASSERT(container->contains(hResp));
  TEST_ASSERT(container->contains(hResp2));
  TEST_ASSERT(!container->contains(cResp));

  container->buildResponseDataObjs();
  container->registerResponses(fm,mainParam);
  
  // evaluate on block 0
  {
     Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
     workset->num_cells = worksetSize; workset->block_id = "block_0";

     panzer::Traits::SetupData setupData;
     setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
     setupData.worksets_->push_back(*workset);

     fm.postRegistrationSetup(setupData);
     fm.writeGraphvizFile();
     fm.preEvaluate<EvalT>(0);
     fm.evaluateFields<EvalT>(*workset);
     fm.postEvaluate<EvalT>(0);
  }

  double sum_dog = 0.0;
  double sum_horse = 0.0;
  for(int i=0;i<worksetSize;i++) {
     sum_dog += double(i)+1.0; 
     sum_horse += -double(i)-5.5; 
  }

  {
     Teuchos::RCP<ResponseData<Traits> > data = container->getResponseData("Functional");

     panzer::Response<panzer::Traits> dResp_o(dResp),
                                      hResp_o(hResp);
     data->fillResponse("Dog",dResp_o);
     data->fillResponse("Horse",hResp_o);
     TEST_EQUALITY(dResp_o.getValue(),sum_dog);
     TEST_EQUALITY(hResp_o.getValue(),sum_horse);
  }

  {
     Teuchos::RCP<ResponseData<Traits> > data = container->getResponseData("Functional");
     container->globalReduction(*comm);

     panzer::Response<panzer::Traits> dResp_o(dResp),
                                      hResp_o(hResp);
     data->fillResponse("Dog",dResp_o);
     data->fillResponse("Horse",hResp_o);

     TEST_EQUALITY(dResp_o.getValue(),comm->getSize()*sum_dog);
     TEST_EQUALITY(hResp_o.getValue(),comm->getSize()*sum_horse);
  }
}

}
}

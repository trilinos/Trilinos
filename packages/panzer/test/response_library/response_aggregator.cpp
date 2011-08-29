#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_config.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_AssemblyEngine.hpp"
#include "Panzer_AssemblyEngine_TemplateManager.hpp"
#include "Panzer_AssemblyEngine_TemplateBuilder.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "TestEvaluators.hpp"

#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {

TEUCHOS_UNIT_TEST(response_assembly, test)
{
  typedef panzer::Traits::Residual EvalT;

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
  PHX::FieldManager<panzer::Traits> fm;
  {
     RCP<PHX::Evaluator<panzer::Traits> > testEval 
        = rcp(new TestEvaluator<EvalT,panzer::Traits>(mainParam));
     fm.registerEvaluator<EvalT>(testEval);
  }

  std::vector<std::string> fields, testFields;
  fields.push_back("Dog");
  fields.push_back("Horse");

  Teuchos::ParameterList p;
  RCP<ResponseAggregatorBase<panzer::Traits> > aggregator 
        = Teuchos::rcp(new ResponseFunctional_Aggregator<EvalT,panzer::Traits>(p));

  TEST_ASSERT(aggregator->clone(p)!=Teuchos::null); 

  RCP<ResponseData<panzer::Traits> > data = aggregator->buildResponseData(fields); 
  TEST_ASSERT(data!=Teuchos::null);
  TEST_ASSERT(data->contains("Dog"));
  TEST_ASSERT(data->contains("Horse"));
  TEST_ASSERT(!data->contains("cat"));

  aggregator->registerAndRequireEvaluators(fm,data,mainParam);

  // evaluate on block 0
  {
     Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
     workset->num_cells = worksetSize; workset->block_id = "block_0";

     panzer::Traits::SetupData setupData;
     setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
     setupData.worksets_->push_back(*workset);

     fm.postRegistrationSetup(setupData);
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

     panzer::Response<panzer::Traits> dResp(ResponseId("Dog","Value")), 
                                      hResp(ResponseId("Horse","Value"));
     data->fillResponse("Dog",dResp);
     data->fillResponse("Horse",hResp);
     TEST_EQUALITY(dResp.getValue(),sum_dog);
     TEST_EQUALITY(hResp.getValue(),sum_horse);
  }

  {
     aggregator->globalReduction(*comm,*data);

     panzer::Response<panzer::Traits> dResp(ResponseId("Dog","Value")), 
                                      hResp(ResponseId("Horse","Value"));
     data->fillResponse("Dog",dResp);
     data->fillResponse("Horse",hResp);
     TEST_EQUALITY(dResp.getValue(),comm->getSize()*sum_dog);
     TEST_EQUALITY(hResp.getValue(),comm->getSize()*sum_horse);
  }
}

}

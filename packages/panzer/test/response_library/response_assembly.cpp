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
#include "Panzer_ResponseLibrary.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_BCStrategy_Factory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "TestEvaluators.hpp"

#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

template <typename ScalarT>
Teuchos::RCP<PHX::FieldTag> buildFieldTag(const std::string & name)
{
   Teuchos::RCP<PHX::DataLayout> layout = Teuchos::rcp(new PHX::MDALayout<Cell>(10));
   return Teuchos::rcp(new PHX::Tag<ScalarT>(name,layout));
}

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

  ResponseLibrary rLibrary;

  RCP<PHX::FieldTag> dogTag = buildFieldTag<double>("Dog");
  RCP<PHX::FieldTag> catTag = buildFieldTag<double>("Cat");
  RCP<PHX::FieldTag> hrsTag = buildFieldTag<double>("Horse");
  
  rLibrary.addResponse(*dogTag,"block_0");
  rLibrary.addResponse(*dogTag,"block_1");
  rLibrary.addResponse(*catTag,"block_0");
  rLibrary.addResponse(*catTag,"block_2");
  rLibrary.addResponse(*hrsTag,"block_0");

  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Dog","block_0");
  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Dog","block_1");
  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Horse","block_0");

  out << rLibrary << std::endl;

  PHX::FieldManager<panzer::Traits> fm0, fm1;
  
  int worksetSize = 10;
  // add test evaluator
  {
     Teuchos::ParameterList p;
     p.set("Workset Size",worksetSize);

     RCP<PHX::Evaluator<panzer::Traits> > e 
        = rcp(new panzer::TestEvaluator<EvalT,panzer::Traits>(p));

     fm0.registerEvaluator<EvalT>(e);
     fm1.registerEvaluator<EvalT>(e);

     // typically we would add field tags associated with TestEvaluator
     // but not during this test
  }

  rLibrary.registerReservedResponses<panzer::Traits::Value>("block_0",comm,worksetSize,fm0);
  rLibrary.registerReservedResponses<panzer::Traits::Value>("block_1",comm,worksetSize,fm1);
  
  Teuchos::RCP<panzer::Workset> workset0 = Teuchos::rcp(new panzer::Workset);
  Teuchos::RCP<panzer::Workset> workset1 = Teuchos::rcp(new panzer::Workset);
  workset0->num_cells = worksetSize; workset0->block_id = "block_0";
  workset1->num_cells = worksetSize; workset1->block_id = "block_1";

  // evaluate on block 0
  {
     panzer::Traits::SetupData setupData;
     setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
     setupData.worksets_->push_back(*workset0);

     fm0.postRegistrationSetup(setupData);
     fm0.preEvaluate<EvalT>(0);
     fm0.evaluateFields<EvalT>(*workset0);
     fm0.postEvaluate<EvalT>(0);
  }

  // evaluate on block 1
  {
     panzer::Traits::SetupData setupData;
     setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
     setupData.worksets_->push_back(*workset1);

     fm1.postRegistrationSetup(setupData);
     fm1.preEvaluate<EvalT>(0);
     fm1.evaluateFields<EvalT>(*workset1);
     fm1.postEvaluate<EvalT>(0);
  }


  // check results on block 0
  {
     RCP<const Response<panzer::Traits::Value> > response_dog
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Dog","block_0");
     RCP<const Response<panzer::Traits::Value> > response_horse
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Horse","block_0");

     double sum_dog = 0.0;
     double sum_horse = 0.0;
     for(int i=0;i<worksetSize;i++) {
        sum_dog += double(i)+1.0; 
        sum_horse += -double(i)-5.5; 
     }
     TEST_EQUALITY(response_dog->getValue(),comm->getSize()*sum_dog); 
     TEST_EQUALITY(response_horse->getValue(),comm->getSize()*sum_horse); 
  }

  // check results on block 1
  {
     RCP<const Response<panzer::Traits::Value> > response_dog
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Dog","block_1");
     TEST_THROW(rLibrary.getVolumeResponse<panzer::Traits::Value>("Horse","block_1"),
                std::logic_error);

     double sum_dog = 0.0;
     for(int i=0;i<worksetSize;i++) {
        sum_dog += double(i)+1.0+44.3;  // from extra for this block
     }
     TEST_EQUALITY(response_dog->getValue(),comm->getSize()*sum_dog); 
  }

  // check summation (aggregation) over element blocks results
  {
     RCP<const Response<panzer::Traits::Value> > response_dog0
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Dog","block_0");
     RCP<const Response<panzer::Traits::Value> > response_dog1
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Dog","block_1");
     RCP<const Response<panzer::Traits::Value> > response_dog
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Dog");

     RCP<const Response<panzer::Traits::Value> > response_horse0
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Horse","block_0");
     RCP<const Response<panzer::Traits::Value> > response_horse
        = rLibrary.getVolumeResponse<panzer::Traits::Value>("Horse");

     TEST_EQUALITY(response_dog->getValue(),response_dog0->getValue()+response_dog1->getValue());
     TEST_EQUALITY(response_horse->getValue(),response_horse0->getValue());
  }
}

}

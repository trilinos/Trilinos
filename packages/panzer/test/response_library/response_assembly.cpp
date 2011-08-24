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
  rLibrary.addResponse(*hrsTag,"block_5");
  rLibrary.addResponse(*hrsTag,"block_0");
  rLibrary.addResponse(*hrsTag,"block_4");

  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Dog","block_0");
  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Dog","block_1");
  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Horse","block_0");
  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Horse","block_5");
  rLibrary.reserveVolumeResponse<panzer::Traits::Value>("Horse","block_4");

  out << rLibrary << std::endl;

  PHX::FieldManager<panzer::Traits> fm;
  
  int worksetSize = 10;
  // add test evaluator
  {
     Teuchos::ParameterList p;
     p.set("Workset Size",worksetSize);

     RCP<PHX::Evaluator<panzer::Traits> > e 
        = rcp(new panzer::TestEvaluator<EvalT,panzer::Traits>(p));

     fm.registerEvaluator<EvalT>(e);

     // typically we would add field tags associated with TestEvaluator
     // but not during this test
  }

  rLibrary.registerReservedResponses<panzer::Traits::Value>("block_0",comm,worksetSize,fm);
  
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  workset->num_cells = worksetSize;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);

  fm.postRegistrationSetup(setupData);

  fm.preEvaluate<EvalT>(0);
  fm.evaluateFields<EvalT>(*workset);
  fm.postEvaluate<EvalT>(0);
}

}

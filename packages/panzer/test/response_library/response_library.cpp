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

#include "UnitTest_UniqueGlobalIndexer.hpp"

#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

template <typename ScalarT>
Teuchos::RCP<PHX::FieldTag> buildFieldTag(const std::string & name)
{
   Teuchos::RCP<PHX::DataLayout> layout = Teuchos::rcp(new PHX::MDALayout<Cell>(10));
   return Teuchos::rcp(new PHX::Tag<ScalarT>(name,layout));
}

TEUCHOS_UNIT_TEST(response_library, add_response)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  ResponseLibrary rLibrary;

  RCP<PHX::FieldTag> dogTag = buildFieldTag<double>("Dog");
  RCP<PHX::FieldTag> catTag = buildFieldTag<double>("Cat");
  
  rLibrary.addResponse(*dogTag,"block_0");
  rLibrary.addResponse(*dogTag,"block_1");
  rLibrary.addResponse(*catTag,"block_0");
  rLibrary.addResponse(*catTag,"block_0"); // same element block

  // test sets
  {
     std::vector<std::string> responses;
     rLibrary.getAvailableVolumeResponses(responses);

     TEST_EQUALITY(responses.size(),std::size_t(2));
     TEST_EQUALITY(responses[0],"Cat"); // assuming map sorts!
     TEST_EQUALITY(responses[1],"Dog");
  }

  // test element block sets
  {
     std::vector<std::pair<std::string,std::set<std::string> > > responses;
     rLibrary.getAvailableVolumeResponses(responses);

     TEST_EQUALITY(responses.size(),std::size_t(2));
     TEST_EQUALITY(responses[0].first,"Cat"); // assuming map sorts!
     TEST_EQUALITY(responses[1].first,"Dog");

     TEST_EQUALITY(responses[0].second.size(),1); 
     TEST_EQUALITY(responses[1].second.size(),2); 
  }

  out << rLibrary << std::endl;
  
}

TEUCHOS_UNIT_TEST(response_library, register_response)
{
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

  Teuchos::ParameterList p;
  p.set(dogTag->name(),true);
  p.set(hrsTag->name(),true);

  rLibrary.checkoutResponses(p);

  out << rLibrary << std::endl;

  { 
     std::vector<std::string> names;
     rLibrary.getCheckedOutVolumeResponses(names);

     TEST_EQUALITY(names.size(),2);
     TEST_EQUALITY(names[0],"Dog");
     TEST_EQUALITY(names[1],"Horse");
  }

  { 
     std::vector<std::string> names;

     rLibrary.getCheckedOutVolumeResponses("block_0",names);
     TEST_EQUALITY(names.size(),2);
     TEST_EQUALITY(names[0],"Dog");
     TEST_EQUALITY(names[1],"Horse");

     rLibrary.getCheckedOutVolumeResponses("block_1",names);
     TEST_EQUALITY(names.size(),1);
     TEST_EQUALITY(names[0],"Dog");

     rLibrary.getCheckedOutVolumeResponses("block_2",names);
     TEST_EQUALITY(names.size(),0);

     rLibrary.getCheckedOutVolumeResponses("block_4",names);
     TEST_EQUALITY(names.size(),1);
     TEST_EQUALITY(names[0],"Horse");
  }

  TEST_THROW(rLibrary.checkoutResponses(p),std::logic_error);

  PHX::FieldManager<panzer::Traits> fm;
  rLibrary.registerResponses<panzer::Traits,panzer::Traits::Residual>(comm,10,"block_0",fm);
}

}

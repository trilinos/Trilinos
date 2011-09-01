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

#include "Panzer_ResponseAggregator_Manager.hpp"
#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {

class TestBuilder {
public:
   template <typename EvalT>
   Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > build() const
   { return Teuchos::null; }
};

template < >
Teuchos::RCP<ResponseAggregatorBase<panzer::Traits> > TestBuilder::
build<panzer::Traits::Residual>() const
{ 
   return Teuchos::rcp(new ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits>); 
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
  using Teuchos::rcpFromRef;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  TestBuilder tb;
  ResponseAggregator_Manager<panzer::Traits> aggMngr;
  aggMngr.defineAggregatorTypeFromBuilder("Functional_2",tb);

  TEST_ASSERT(aggMngr.isAggregator("Functional_2"));
  TEST_ASSERT(aggMngr.isAggregator<panzer::Traits::Residual>());
  TEST_ASSERT(aggMngr.isAggregator<panzer::Traits::Residual>("Functional_2"));

  TEST_ASSERT(!aggMngr.isAggregator("None"));
  TEST_ASSERT(!aggMngr.isAggregator<panzer::Traits::Jacobian>());
  TEST_ASSERT(!aggMngr.isAggregator<panzer::Traits::Residual>("None"));
  TEST_ASSERT(!aggMngr.isAggregator<panzer::Traits::Jacobian>("Functional_2"));

  TEST_THROW(aggMngr.getAggregator<panzer::Traits::Residual>("None"),std::logic_error);
  TEST_THROW(aggMngr.getAggregator<panzer::Traits::Jacobian>("Functional_2"),std::logic_error);

  RCP<const ResponseAggregatorBase<panzer::Traits> > aggregator; 
  {
     aggregator = rcpFromRef(aggMngr.getAggregator<panzer::Traits::Residual>("Functional_2"));

     RCP<const ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits> > f_aggregator 
        = rcp_dynamic_cast<const ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits> >(aggregator);
     TEST_ASSERT(f_aggregator!=Teuchos::null);
  }

  ResponseAggregator_Manager<panzer::Traits>::defineDefaultAggregators(aggMngr);

  {
     aggregator = rcpFromRef(aggMngr.getAggregator<panzer::Traits::Residual>("Functional"));

     RCP<const ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits> > f_aggregator 
        = rcp_dynamic_cast<const ResponseFunctional_Aggregator<panzer::Traits::Residual,panzer::Traits> >(aggregator);
     TEST_ASSERT(f_aggregator!=Teuchos::null);
  }
}

}

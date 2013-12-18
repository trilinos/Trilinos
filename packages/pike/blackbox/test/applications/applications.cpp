#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_LinearHeatConductionModelEvaluator.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

namespace pike_test {

  TEUCHOS_UNIT_TEST(app, LinearHeatConduction_BlockJacobi)
  {
    using Teuchos::RCP;

    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator> leftWall = 
      pike_test::linearHeatConductionModelEvaluator(globalComm,"left wall",pike_test::LinearHeatConductionModelEvaluator::T_RIGHT_IS_RESPONSE);

    Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator> middleWall = 
      pike_test::linearHeatConductionModelEvaluator(globalComm,"middle wall",pike_test::LinearHeatConductionModelEvaluator::T_RIGHT_IS_RESPONSE);

    Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator> rightWall = 
      pike_test::linearHeatConductionModelEvaluator(globalComm,"right wall",pike_test::LinearHeatConductionModelEvaluator::Q_IS_RESPONSE);

    

  }

}

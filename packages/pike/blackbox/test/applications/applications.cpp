#include "Teuchos_UnitTestHarness.hpp"
#include "Pike_BlackBox_config.hpp"
#include "Pike_LinearHeatConductionModelEvaluator.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

namespace pike {

  TEUCHOS_UNIT_TEST(app, LinearHeatConduction_BlockJacobi)
  {
    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    typedef pike::LinearHeatConductionModelEvaluator::Mode::Q_IS_RESPONSE Q_IS_RESPONSE;
   
    Teuchos::RCP<pike::LinearHeatConductionModelEvaluator> leftWall = 
      pike::linearHeatConductionModelEvaluator(globalComm,"left wall",T_RIGHT_IS_RESPONSE);

    Teuchos::RCP<pike::LinearHeatConductionModelEvaluator> middleWall = 
      pike::linearHeatConductionModelEvaluator(globalComm,"middle wall",T_RIGHT_IS_RESPONSE);

    Teuchos::RCP<pike::LinearHeatConductionModelEvaluator> rightWall = 
      pike::linearHeatConductionModelEvaluator(globalComm,"right wall",Q_IS_RESPONSE);


  }

}

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

#include "Panzer_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntegrationValues.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_Normals.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"

// for making explicit instantiated tests easier 
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(normals,test2d,TYPE)

namespace panzer {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(normals,test2d,EvalType)
{
  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
     Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
  #endif
 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  // build a dummy workset
  //////////////////////////////////////////////////////////
  typedef Intrepid::FieldContainer<double> FieldArray;
  int numCells = 2, numVerts = 4, dim = 2;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  FieldArray & coords = workset->cell_vertex_coordinates;
  coords.resize(numCells,numVerts,dim);
  coords(0,0,0) = 1.0; coords(0,0,1) = 0.0;
  coords(0,1,0) = 1.0; coords(0,1,1) = 1.0;
  coords(0,2,0) = 0.0; coords(0,2,1) = 1.0;
  coords(0,3,0) = 0.0; coords(0,3,1) = 0.0;

  coords(1,0,0) = 1.0; coords(1,0,1) = 1.0;
  coords(1,1,0) = 2.0; coords(1,1,1) = 2.0;
  coords(1,2,0) = 1.0; coords(1,2,1) = 3.0;
  coords(1,3,0) = 0.0; coords(1,3,1) = 2.0;

  int quadOrder = 5;
  panzer::CellData cellData(2,2,1);
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  out << "num quad points = " << quadRule->num_points << std::endl;
  Teuchos::RCP<panzer::IntegrationValues<double,FieldArray> > quadValues = Teuchos::rcp(new panzer::IntegrationValues<double,FieldArray>);
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  workset->cell_local_ids.push_back(0); workset->cell_local_ids.push_back(1);
  workset->num_cells = numCells;
  workset->block_id = "eblock-0_0";
  workset->ir_degrees = Teuchos::rcp(new std::vector<int>);
  workset->ir_degrees->push_back(quadRule->cubature_degree);
  workset->int_rules.push_back(quadValues);

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
     = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>); 

  // typedef panzer::Traits::Residual EvalType;
  typedef Sacado::ScalarType<typename EvalType::ScalarT> ScalarType;
  typedef Sacado::ScalarValue<typename EvalType::ScalarT> ScalarValue;
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim> > normalsPtr;
  {
     Teuchos::ParameterList p;
     p.set("Name","Norms");
     p.set("IR",quadRule);
     p.set("Side ID",1);
    
     RCP<panzer::Normals<EvalType,panzer::Traits> > normEval  
        = rcp(new panzer::Normals<EvalType,panzer::Traits>(p));
     RCP<PHX::Evaluator<panzer::Traits> > eval = normEval;

     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(normEval->getFieldTag());

     const PHX::FieldTag & ft = normEval->getFieldTag();
     normalsPtr = rcp(new
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim>(ft.name(),quadRule->dl_vector));
  }
  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim> & normals = *normalsPtr;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);

  fm->preEvaluate<EvalType>(0);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(normals);

  TEST_EQUALITY(normals.rank(),3);
  TEST_EQUALITY(normals.size(),numCells*quadRule->num_points*dim);
  normals.print(out,true);
  for(int i=0;i<numCells;i++) {

     // useful for checking if normals are consistent: transformation is
     // affine!
     double nx0 = ScalarValue::eval(normals(i,0,0));
     double ny0 = ScalarValue::eval(normals(i,0,1));

     for(int v=0;v<quadRule->num_points;v++) {
        double nx = ScalarValue::eval(normals(i,v,0)); 
        double ny = ScalarValue::eval(normals(i,v,1)); 
 
        TEST_FLOATING_EQUALITY(nx*nx+ny*ny,1.0,1e-15);

        // check point consistency
        TEST_FLOATING_EQUALITY(nx,nx0,1e-15);
        TEST_FLOATING_EQUALITY(ny,ny0,1e-15);
     }
  }

  // check cell 0
  {
     double nx = ScalarValue::eval(normals(0,0,0));
     double ny = ScalarValue::eval(normals(0,0,1));
   
     TEST_FLOATING_EQUALITY(nx,0.0,1e-15);
     TEST_FLOATING_EQUALITY(ny,1.0,1e-15);
  }

  // check cell 1
  {
     double nx = ScalarValue::eval(normals(1,0,0));
     double ny = ScalarValue::eval(normals(1,0,1));
     double sqrt2 = std::sqrt(2.0);
   
     TEST_FLOATING_EQUALITY(nx,1.0/sqrt2,1e-15);
     TEST_FLOATING_EQUALITY(ny,1.0/sqrt2,1e-15);
  }
}

typedef Traits::Residual ResidualType;
typedef Traits::Jacobian JacobianType;

UNIT_TEST_GROUP(ResidualType)
UNIT_TEST_GROUP(JacobianType)

#ifdef HAVE_STOKHOS
   typedef Traits::SGResidual SGResidualType;
   typedef Traits::SGJacobian SGJacobianType;

   UNIT_TEST_GROUP(SGResidualType)
   UNIT_TEST_GROUP(SGJacobianType)
#endif

}

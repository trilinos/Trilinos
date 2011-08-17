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
#include "Panzer_Integrator_Scalar.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"

#include "UnitValueEvaluator.hpp"

// for making explicit instantiated tests easier 
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(integrator_scalar_side,test2d,TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(integrator_scalar_side,test3d,TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(integrator_scalar,test3d,TYPE)

namespace panzer {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(integrator_scalar_side,test2d,EvalType)
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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell> > integralPtr;
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<panzer::Cell>(quadRule->dl_scalar->dimension(0)));

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     Teuchos::ParameterList p;
     p.set("Integral Name","2_x_Length");
     p.set("Integrand Name","Unit Value");
     p.set("Multiplier",2.0);
     p.set("IR",quadRule);
    
     RCP<panzer::Integrator_Scalar<EvalType,panzer::Traits> > eval 
        = rcp(new panzer::Integrator_Scalar<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(eval->getFieldTag());

     const PHX::FieldTag & ft = eval->getFieldTag();
     integralPtr = rcp(new
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell>(ft.name(),dl_cell));
  }

  {
     Teuchos::ParameterList p;
     p.set("Name","Unit Value");
     p.set("IR",quadRule);
    
     RCP<PHX::Evaluator<panzer::Traits> > eval 
        = rcp(new UnitValueEvaluator<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
  }
  PHX::MDField<typename EvalType::ScalarT,panzer::Cell> & integral = *integralPtr;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);

  fm->preEvaluate<EvalType>(0);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(integral);

  integral.print(out,true);

  TEST_FLOATING_EQUALITY(ScalarValue::eval(integral(0)),2.0,1e-15);
  TEST_FLOATING_EQUALITY(ScalarValue::eval(integral(1)),2.0*std::sqrt(2.0),1e-15);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(integrator_scalar_side,test3d,EvalType)
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
  int numCells = 2, numVerts = 8, dim = 3;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  FieldArray & coords = workset->cell_vertex_coordinates;
  coords.resize(numCells,numVerts,dim);
  coords(0,0,0) = 1.0; coords(0,0,1) = 0.0; coords(0,0,2) = 0.0;
  coords(0,1,0) = 1.0; coords(0,1,1) = 1.0; coords(0,1,2) = 0.0;
  coords(0,2,0) = 0.0; coords(0,2,1) = 1.0; coords(0,2,2) = 0.0;
  coords(0,3,0) = 0.0; coords(0,3,1) = 0.0; coords(0,3,2) = 0.0;
  coords(0,4,0) = 1.0; coords(0,4,1) = 0.0; coords(0,4,2) = 1.0;
  coords(0,5,0) = 1.0; coords(0,5,1) = 1.0; coords(0,5,2) = 1.0;
  coords(0,6,0) = 0.0; coords(0,6,1) = 1.0; coords(0,6,2) = 1.0;
  coords(0,7,0) = 0.0; coords(0,7,1) = 0.0; coords(0,7,2) = 1.0;

  coords(1,0,0) = 0.0; coords(1,0,1) = 0.0; coords(1,0,2) = 0.0;
  coords(1,1,0) =-1.0; coords(1,1,1) = 1.0; coords(1,1,2) = 0.0;
  coords(1,2,0) = 2.0; coords(1,2,1) = 2.0; coords(1,2,2) = 0.0;
  coords(1,3,0) = 1.0; coords(1,3,1) = 1.0; coords(1,3,2) = 0.0;
  coords(1,4,0) = 0.0; coords(1,4,1) = 0.0; coords(1,4,2) = 2.0;
  coords(1,5,0) =-1.0; coords(1,5,1) = 1.0; coords(1,5,2) = 2.0;
  coords(1,6,0) = 2.0; coords(1,6,1) = 2.0; coords(1,6,2) = 2.0;
  coords(1,7,0) = 1.0; coords(1,7,1) = 1.0; coords(1,7,2) = 2.0;

  int quadOrder = 5;
  panzer::CellData cellData(2,dim,0);
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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell> > integralPtr;
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<panzer::Cell>(quadRule->dl_scalar->dimension(0)));

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     Teuchos::ParameterList p;
     p.set("Integral Name","2_x_Length");
     p.set("Integrand Name","Unit Value");
     p.set("Multiplier",2.0);
     p.set("IR",quadRule);
    
     RCP<panzer::Integrator_Scalar<EvalType,panzer::Traits> > eval 
        = rcp(new panzer::Integrator_Scalar<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(eval->getFieldTag());

     const PHX::FieldTag & ft = eval->getFieldTag();
     integralPtr = rcp(new
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell>(ft.name(),dl_cell));
  }

  {
     Teuchos::ParameterList p;
     p.set("Name","Unit Value");
     p.set("IR",quadRule);
    
     RCP<PHX::Evaluator<panzer::Traits> > eval 
        = rcp(new UnitValueEvaluator<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
  }
  PHX::MDField<typename EvalType::ScalarT,panzer::Cell> & integral = *integralPtr;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);

  fm->preEvaluate<EvalType>(0);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(integral);

  integral.print(out,true);

  TEST_FLOATING_EQUALITY(ScalarValue::eval(integral(0)),2.0,1e-15);
  TEST_FLOATING_EQUALITY(ScalarValue::eval(integral(1)),4.0*std::sqrt(2),1e-15);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(integrator_scalar,test3d,EvalType)
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
  int numCells = 2, numVerts = 8, dim = 3;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  FieldArray & coords = workset->cell_vertex_coordinates;
  coords.resize(numCells,numVerts,dim);
  coords(0,0,0) = 1.0; coords(0,0,1) = 0.0; coords(0,0,2) = 0.0;
  coords(0,1,0) = 1.0; coords(0,1,1) = 1.0; coords(0,1,2) = 0.0;
  coords(0,2,0) = 0.0; coords(0,2,1) = 1.0; coords(0,2,2) = 0.0;
  coords(0,3,0) = 0.0; coords(0,3,1) = 0.0; coords(0,3,2) = 0.0;
  coords(0,4,0) = 1.0; coords(0,4,1) = 0.0; coords(0,4,2) = 1.0;
  coords(0,5,0) = 1.0; coords(0,5,1) = 1.0; coords(0,5,2) = 1.0;
  coords(0,6,0) = 0.0; coords(0,6,1) = 1.0; coords(0,6,2) = 1.0;
  coords(0,7,0) = 0.0; coords(0,7,1) = 0.0; coords(0,7,2) = 1.0;

  coords(1,0,0) = 0.0; coords(1,0,1) = 0.0; coords(1,0,2) = 0.0;
  coords(1,1,0) =-1.0; coords(1,1,1) = 1.0; coords(1,1,2) = 0.0;
  coords(1,2,0) = 2.0; coords(1,2,1) = 2.0; coords(1,2,2) = 0.0;
  coords(1,3,0) = 1.0; coords(1,3,1) = 1.0; coords(1,3,2) = 0.0;
  coords(1,4,0) = 0.0; coords(1,4,1) = 0.0; coords(1,4,2) = 2.0;
  coords(1,5,0) =-1.0; coords(1,5,1) = 1.0; coords(1,5,2) = 2.0;
  coords(1,6,0) = 2.0; coords(1,6,1) = 2.0; coords(1,6,2) = 2.0;
  coords(1,7,0) = 1.0; coords(1,7,1) = 1.0; coords(1,7,2) = 2.0;

  int quadOrder = 5;
  panzer::CellData cellData(2,dim);
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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell> > integralPtr;
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<panzer::Cell>(quadRule->dl_scalar->dimension(0)));

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     Teuchos::ParameterList p;
     p.set("Integral Name","2_x_Length");
     p.set("Integrand Name","Unit Value");
     p.set("Multiplier",2.0);
     p.set("IR",quadRule);
    
     RCP<panzer::Integrator_Scalar<EvalType,panzer::Traits> > eval 
        = rcp(new panzer::Integrator_Scalar<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(eval->getFieldTag());

     const PHX::FieldTag & ft = eval->getFieldTag();
     integralPtr = rcp(new
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell>(ft.name(),dl_cell));
  }

  {
     Teuchos::ParameterList p;
     p.set("Name","Unit Value");
     p.set("IR",quadRule);
    
     RCP<PHX::Evaluator<panzer::Traits> > eval 
        = rcp(new UnitValueEvaluator<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
  }
  PHX::MDField<typename EvalType::ScalarT,panzer::Cell> & integral = *integralPtr;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);

  fm->preEvaluate<EvalType>(0);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(integral);

  integral.print(out,true);

  TEST_FLOATING_EQUALITY(ScalarValue::eval(integral(0)),2.0,1e-15);
  TEST_FLOATING_EQUALITY(ScalarValue::eval(integral(1)),8.0,1e-15);
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

typedef Traits::Residual ResidualType;
typedef Traits::Jacobian JacobianType;

UNIT_TEST_GROUP(ResidualType)
UNIT_TEST_GROUP(JacobianType)

#ifdef HAVE_STOKHOS
   typedef Traits::SGResidual SGResidualType;
   typedef Traits::SGJacobian SGJacobianType;

   // these are failing with Stokhos enabled! 8/18/2011
   UNIT_TEST_GROUP(SGResidualType)
   UNIT_TEST_GROUP(SGJacobianType)
#endif

}

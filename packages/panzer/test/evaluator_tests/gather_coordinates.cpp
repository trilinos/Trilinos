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

#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_GatherIntegrationCoordinates.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"

#include "UnitValueEvaluator.hpp"

// for making explicit instantiated tests easier 
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(gather_coordinates,basis,TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(gather_coordinates,integration,TYPE)

namespace panzer {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(gather_coordinates,basis,EvalType)
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

  // build topology, basis, integration rule, and basis layout
  int quadOrder = 5;
  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
  panzer::CellData cellData(2,2,1,topo);
  RCP<PureBasis> basis = rcp(new PureBasis("Q1",2,topo));
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  RCP<BasisIRLayout> basisLayout = rcp(new BasisIRLayout(basis,*quadRule));

  RCP<IntegrationValues<double,FieldArray> > quadValues = rcp(new IntegrationValues<double,FieldArray>);
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  RCP<BasisValues<double,FieldArray> > basisValues = rcp(new BasisValues<double,FieldArray>);
  basisValues->setupArrays(basisLayout);
  basisValues->evaluateValues(quadValues->cub_points,
                              quadValues->jac_inv,
                              quadValues->weighted_measure,
                              quadValues->node_coordinates);

  // setup generic workset stuff
  workset->cell_local_ids.push_back(0); workset->cell_local_ids.push_back(1);
  workset->num_cells = numCells;
  workset->block_id = "eblock-0_0";

  // setup integration rule
  workset->ir_degrees = Teuchos::rcp(new std::vector<int>);
  workset->ir_degrees->push_back(quadRule->cubature_degree);
  workset->int_rules.push_back(quadValues);

  // setup basis functions
  workset->basis_names = Teuchos::rcp(new std::vector<std::string>);
  workset->basis_names->push_back(basis->name());
  workset->bases.push_back(basisValues);

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
     = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>); 

  // typedef panzer::Traits::Residual EvalType;
  typedef Sacado::ScalarType<typename EvalType::ScalarT> ScalarType;
  typedef Sacado::ScalarValue<typename EvalType::ScalarT> ScalarValue;
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim> > fmCoordsPtr;
  Teuchos::RCP<PHX::DataLayout> dl_coords = basis->coordinates;

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     RCP<panzer::GatherBasisCoordinates<EvalType,panzer::Traits> > eval 
        = rcp(new panzer::GatherBasisCoordinates<EvalType,panzer::Traits>(*basis));

     const std::vector<RCP<PHX::FieldTag> > & evalFields = eval->evaluatedFields();
     const std::vector<RCP<PHX::FieldTag> > & depFields = eval->dependentFields();

     TEST_EQUALITY(evalFields.size(),1);
     TEST_EQUALITY(depFields.size(), 0);
     std::string ref_fieldName = GatherBasisCoordinates<EvalType,panzer::Traits>::fieldName(basis->name());
     TEST_EQUALITY(ref_fieldName,evalFields[0]->name());

     fm->registerEvaluator<EvalType>(eval);



     const PHX::FieldTag & ft = *evalFields[0];
     fm->requireField<EvalType>(ft);
     fmCoordsPtr = rcp(new
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim>(ft.name(),dl_coords));
  }

  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim> & fmCoords = *fmCoordsPtr;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(fmCoords);

  fmCoords.print(out,true);

  for(int i=0;i<fmCoords.size();i++)
     TEST_EQUALITY(ScalarValue::eval(fmCoords[i]),coords[i]);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(gather_coordinates,integration,EvalType)
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

  // build topology, basis, integration rule, and basis layout
  int quadOrder = 5;
  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
  panzer::CellData cellData(2,2,1,topo);
  RCP<PureBasis> basis = rcp(new PureBasis("Q1",2,topo));
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  RCP<BasisIRLayout> basisLayout = rcp(new BasisIRLayout(basis,*quadRule));

  RCP<IntegrationValues<double,FieldArray> > quadValues = rcp(new IntegrationValues<double,FieldArray>);
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  RCP<BasisValues<double,FieldArray> > basisValues = rcp(new BasisValues<double,FieldArray>);
  basisValues->setupArrays(basisLayout);
  basisValues->evaluateValues(quadValues->cub_points,
                              quadValues->jac_inv,
                              quadValues->weighted_measure,
                              quadValues->node_coordinates);

  // setup generic workset stuff
  workset->cell_local_ids.push_back(0); workset->cell_local_ids.push_back(1);
  workset->num_cells = numCells;
  workset->block_id = "eblock-0_0";

  // setup integration rule
  workset->ir_degrees = Teuchos::rcp(new std::vector<int>);
  workset->ir_degrees->push_back(quadRule->cubature_degree);
  workset->int_rules.push_back(quadValues);

  // setup basis functions
  workset->basis_names = Teuchos::rcp(new std::vector<std::string>);
  workset->basis_names->push_back(basis->name());
  workset->bases.push_back(basisValues);

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
     = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>); 

  // typedef panzer::Traits::Residual EvalType;
  typedef Sacado::ScalarType<typename EvalType::ScalarT> ScalarType;
  typedef Sacado::ScalarValue<typename EvalType::ScalarT> ScalarValue;
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim> > fmCoordsPtr;
  Teuchos::RCP<PHX::DataLayout> dl_coords = quadRule->dl_vector;

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     RCP<panzer::GatherIntegrationCoordinates<EvalType,panzer::Traits> > eval 
        = rcp(new panzer::GatherIntegrationCoordinates<EvalType,panzer::Traits>(*quadRule));

     const std::vector<RCP<PHX::FieldTag> > & evalFields = eval->evaluatedFields();
     const std::vector<RCP<PHX::FieldTag> > & depFields = eval->dependentFields();

     TEST_EQUALITY(evalFields.size(),1);
     TEST_EQUALITY(depFields.size(), 0);
     std::string ref_fieldName = GatherIntegrationCoordinates<EvalType,panzer::Traits>::fieldName(quadRule->cubature_degree);
     TEST_EQUALITY(ref_fieldName,evalFields[0]->name());

     fm->registerEvaluator<EvalType>(eval);

     const PHX::FieldTag & ft = *evalFields[0];
     fm->requireField<EvalType>(ft);
     fmCoordsPtr = rcp(new
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim>(ft.name(),dl_coords));
  }

  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point,panzer::Dim> & fmCoords = *fmCoordsPtr;

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(fmCoords);

  fmCoords.print(out,true);

  for(int i=0;i<fmCoords.size();i++)
     TEST_EQUALITY(ScalarValue::eval(fmCoords[i]),quadValues->ip_coordinates[i]);
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

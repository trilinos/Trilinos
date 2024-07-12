// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#include "Kokkos_View_Fad.hpp"
#include "PanzerDiscFE_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_GatherIntegrationCoordinates.hpp"

#include "Phalanx_FieldManager.hpp"

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
     Teuchos::RCP<const Teuchos::MpiComm<int> > eComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  #else
      auto eComm = Teuchos::rcp(Teuchos::DefaultComm<int>::getComm());
  #endif

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  // build a dummy workset
  //////////////////////////////////////////////////////////
  // typedef Kokkos::DynRankView<double,PHX::Device> FieldArray;
  int numCells = 2, numVerts = 4, dim = 2;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);

  MDFieldArrayFactory af("",true);
  workset->cell_node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("coords",numCells,numVerts,dim);
  Workset::CellCoordArray coords = workset->cell_node_coordinates;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int ) {
      coords(0,0,0) = 1.0; coords(0,0,1) = 0.0;
      coords(0,1,0) = 1.0; coords(0,1,1) = 1.0;
      coords(0,2,0) = 0.0; coords(0,2,1) = 1.0;
      coords(0,3,0) = 0.0; coords(0,3,1) = 0.0;

      coords(1,0,0) = 1.0; coords(1,0,1) = 1.0;
      coords(1,1,0) = 2.0; coords(1,1,1) = 2.0;
      coords(1,2,0) = 1.0; coords(1,2,1) = 3.0;
      coords(1,3,0) = 0.0; coords(1,3,1) = 2.0;
    });
  Kokkos::fence();

  // build topology, basis, integration rule, and basis layout
  int quadOrder = 5;
  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
  panzer::CellData cellData(2,1,topo);
  RCP<PureBasis> basis = rcp(new PureBasis("Q1",1,2,topo));
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  RCP<BasisIRLayout> basisLayout = rcp(new BasisIRLayout(basis,*quadRule));

  RCP<IntegrationValues2<double> > quadValues = rcp(new IntegrationValues2<double>("",true));
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  RCP<BasisValues2<double> > basisValues = rcp(new BasisValues2<double>("",true));
  basisValues->setupArrays(basisLayout);
  basisValues->evaluateValues(quadValues->cub_points,
                              quadValues->jac,
                              quadValues->jac_det,
                              quadValues->jac_inv,
                              quadValues->weighted_measure,
                              coords);
                              // quadValues->node_coordinates);

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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim> > fmCoordsPtr;
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
         PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim>(ft.name(),dl_coords));
  }

  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim> & fmCoords = *fmCoordsPtr;

  panzer::Traits::SD setupData;
  {
    auto worksets = rcp(new std::vector<panzer::Workset>);
    worksets->push_back(*workset);
    setupData.worksets_ = worksets;
  }

  std::vector<PHX::index_size_type> derivative_dimensions;
  derivative_dimensions.push_back(4);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(derivative_dimensions);
#endif
  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<EvalType>(fmCoords);

  fmCoords.print(out,true);
  int error = false;
  Kokkos::parallel_reduce(fmCoords.extent_int(0), KOKKOS_LAMBDA (int cell, int &err) {
    for(int pt=0;pt<fmCoords.extent_int(1);++pt)
      for(int d=0;d<fmCoords.extent_int(2);++d)
	err |= Sacado::scalarValue(fmCoords(cell,pt,d))!=coords(cell,pt,d);
    },error);
  TEST_EQUALITY(error, false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(gather_coordinates,integration,EvalType)
{

  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<const Teuchos::MpiComm<int> > eComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  #else
      auto eComm = Teuchos::rcp(Teuchos::DefaultComm<int>::getComm());
  #endif

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  // build a dummy workset
  //////////////////////////////////////////////////////////
  // typedef Kokkos::DynRankView<double,PHX::Device> FieldArray;
  int numCells = 2, numVerts = 4, dim = 2;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);

  MDFieldArrayFactory af("",true);
  workset->cell_node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("coords",numCells,numVerts,dim);
  Workset::CellCoordArray coords = workset->cell_node_coordinates;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int ) {
      coords(0,0,0) = 1.0; coords(0,0,1) = 0.0;
      coords(0,1,0) = 1.0; coords(0,1,1) = 1.0;
      coords(0,2,0) = 0.0; coords(0,2,1) = 1.0;
      coords(0,3,0) = 0.0; coords(0,3,1) = 0.0;

      coords(1,0,0) = 1.0; coords(1,0,1) = 1.0;
      coords(1,1,0) = 2.0; coords(1,1,1) = 2.0;
      coords(1,2,0) = 1.0; coords(1,2,1) = 3.0;
      coords(1,3,0) = 0.0; coords(1,3,1) = 2.0;
    });
  // build topology, basis, integration rule, and basis layout
  int quadOrder = 5;
  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
  panzer::CellData cellData(2,1,topo);
  RCP<PureBasis> basis = rcp(new PureBasis("Q1",1,2,topo));
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  RCP<BasisIRLayout> basisLayout = rcp(new BasisIRLayout(basis,*quadRule));

  RCP<IntegrationValues2<double> > quadValues = rcp(new IntegrationValues2<double>("",true));
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  RCP<BasisValues2<double> > basisValues = rcp(new BasisValues2<double>("",true));
  basisValues->setupArrays(basisLayout);
  basisValues->evaluateValues(quadValues->cub_points,
                              quadValues->jac,
                              quadValues->jac_det,
                              quadValues->jac_inv,
                              quadValues->weighted_measure,
                              coords);
                              // quadValues->node_coordinates);

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

  panzer::Traits::SD setupData;
  {
    auto worksets = rcp(new std::vector<panzer::Workset>);
    worksets->push_back(*workset);
    setupData.worksets_ = worksets;
  }

  std::vector<PHX::index_size_type> derivative_dimensions;
  derivative_dimensions.push_back(4);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(derivative_dimensions);
#endif
  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<EvalType>(fmCoords);

  fmCoords.print(out,true);

  int error = false;
  auto q_coords = quadValues->ip_coordinates;
  Kokkos::parallel_reduce(fmCoords.extent_int(0), KOKKOS_LAMBDA (int cell, int &err) {
    for(int pt=0;pt<fmCoords.extent_int(1);++pt)
      for(int d=0;d<fmCoords.extent_int(2);++d)
	err |= Sacado::scalarValue(fmCoords(cell,pt,d)) != q_coords(cell,pt,d);
    },error);
  TEST_EQUALITY(error, false);

}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

typedef Traits::Residual ResidualType;
typedef Traits::Jacobian JacobianType;

UNIT_TEST_GROUP(ResidualType)
UNIT_TEST_GROUP(JacobianType)

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
typedef Traits::Hessian HessianType;
UNIT_TEST_GROUP(HessianType)
#endif

}

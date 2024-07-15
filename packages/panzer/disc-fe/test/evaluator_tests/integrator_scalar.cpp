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
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Panzer_Normals.hpp"
#include "Panzer_Integrator_Scalar.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_MDField_UnmanagedAllocator.hpp"
#include "Phalanx_Evaluator_UnmanagedFieldDummy.hpp"

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
  auto coords_v = coords.get_static_view();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int) {
      coords_v(0,0,0) = 1.0; coords_v(0,0,1) = 0.0;
      coords_v(0,1,0) = 1.0; coords_v(0,1,1) = 1.0;
      coords_v(0,2,0) = 0.0; coords_v(0,2,1) = 1.0;
      coords_v(0,3,0) = 0.0; coords_v(0,3,1) = 0.0;

      coords_v(1,0,0) = 1.0; coords_v(1,0,1) = 1.0;
      coords_v(1,1,0) = 2.0; coords_v(1,1,1) = 2.0;
      coords_v(1,2,0) = 1.0; coords_v(1,2,1) = 3.0;
      coords_v(1,3,0) = 0.0; coords_v(1,3,1) = 2.0;
    });

  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

  int quadOrder = 5;
  panzer::CellData cellData(2,1,topo);
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  out << "num quad points = " << quadRule->num_points << std::endl;
  Teuchos::RCP<panzer::IntegrationValues2<double> > quadValues = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell> > integralPtr;
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<panzer::Cell>(quadRule->dl_scalar->extent(0)));

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
  std::vector<PHX::index_size_type> hess_derivative_dimensions;
  hess_derivative_dimensions.push_back(4);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(hess_derivative_dimensions);
#endif

  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<EvalType>(integral);

  integral.print(out,true);

  integral.print(out,true);
  auto integral_v = integral.get_static_view();
  auto integral_h = Kokkos::create_mirror_view(integral_v);
  Kokkos::deep_copy(integral_h, integral_v);
  TEST_FLOATING_EQUALITY(Sacado::scalarValue(integral_h(0)),2.0,1e-15);
  TEST_FLOATING_EQUALITY(Sacado::scalarValue(integral_h(1)),2.0*std::sqrt(2.0),1e-15);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(integrator_scalar_side,test3d,EvalType)
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
  int numCells = 2, numVerts = 8, dim = 3;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  MDFieldArrayFactory af("",true);
  workset->cell_node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("coords",numCells,numVerts,dim);
  Workset::CellCoordArray coords = workset->cell_node_coordinates;
  auto coords_v = coords.get_static_view();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int) {
      coords_v(0,0,0) = 1.0; coords_v(0,0,1) = 0.0; coords_v(0,0,2) = 0.0;
      coords_v(0,1,0) = 1.0; coords_v(0,1,1) = 1.0; coords_v(0,1,2) = 0.0;
      coords_v(0,2,0) = 0.0; coords_v(0,2,1) = 1.0; coords_v(0,2,2) = 0.0;
      coords_v(0,3,0) = 0.0; coords_v(0,3,1) = 0.0; coords_v(0,3,2) = 0.0;
      coords_v(0,4,0) = 1.0; coords_v(0,4,1) = 0.0; coords_v(0,4,2) = 1.0;
      coords_v(0,5,0) = 1.0; coords_v(0,5,1) = 1.0; coords_v(0,5,2) = 1.0;
      coords_v(0,6,0) = 0.0; coords_v(0,6,1) = 1.0; coords_v(0,6,2) = 1.0;
      coords_v(0,7,0) = 0.0; coords_v(0,7,1) = 0.0; coords_v(0,7,2) = 1.0;

      coords_v(1,0,0) = 0.0; coords_v(1,0,1) = 0.0; coords_v(1,0,2) = 0.0;
      coords_v(1,1,0) =-1.0; coords_v(1,1,1) = 1.0; coords_v(1,1,2) = 0.0;
      coords_v(1,2,0) = 2.0; coords_v(1,2,1) = 2.0; coords_v(1,2,2) = 0.0;
      coords_v(1,3,0) = 1.0; coords_v(1,3,1) = 1.0; coords_v(1,3,2) = 0.0;
      coords_v(1,4,0) = 0.0; coords_v(1,4,1) = 0.0; coords_v(1,4,2) = 2.0;
      coords_v(1,5,0) =-1.0; coords_v(1,5,1) = 1.0; coords_v(1,5,2) = 2.0;
      coords_v(1,6,0) = 2.0; coords_v(1,6,1) = 2.0; coords_v(1,6,2) = 2.0;
      coords_v(1,7,0) = 1.0; coords_v(1,7,1) = 1.0; coords_v(1,7,2) = 2.0;
  });
  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));

  int quadOrder = 5;
  panzer::CellData cellData(2,0,topo);
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  out << "num quad points = " << quadRule->num_points << std::endl;
  Teuchos::RCP<panzer::IntegrationValues2<double> > quadValues = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell> > integralPtr;
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<panzer::Cell>(quadRule->dl_scalar->extent(0)));

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     Teuchos::ParameterList p;
     p.set("Integral Name","2_x_Length");
     p.set("Integrand Name","Unit Value");
     p.set("Multiplier",2.0);
     p.set("IR",quadRule);
     RCP<const std::vector<std::string>> fms = rcp(new std::vector<std::string>{"Dummy Field"});
     p.set("Field Multipliers",fms);

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

  panzer::Traits::SD setupData;
  {
    auto worksets = rcp(new std::vector<panzer::Workset>);
    worksets->push_back(*workset);
    setupData.worksets_ = worksets;
  }

  std::vector<PHX::index_size_type> derivative_dimensions;
  derivative_dimensions.push_back(4);
  {
    PHX::MDField<typename EvalType::ScalarT,Cell,IP> field_mult =
      PHX::allocateUnmanagedMDField<typename EvalType::ScalarT,Cell,IP>("Dummy Field",quadRule->dl_scalar,derivative_dimensions);
    Kokkos::deep_copy(field_mult.get_static_view(),2.0);
    fm->setUnmanagedField<EvalType>(field_mult);
    auto e = rcp(new PHX::UnmanagedFieldDummy<EvalType,panzer::Traits,PHX::MDField<typename EvalType::ScalarT,Cell,IP>>(field_mult));
    fm->registerEvaluator<EvalType>(e);
  }

  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  std::vector<PHX::index_size_type> hess_derivative_dimensions;
  hess_derivative_dimensions.push_back(4);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(hess_derivative_dimensions);
#endif

  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<EvalType>(integral);

  integral.print(out,true);

  integral.print(out,true);
  auto integral_v = integral.get_static_view();
  auto integral_h = Kokkos::create_mirror_view(integral_v);
  Kokkos::deep_copy(integral_h, integral_v);
  TEST_FLOATING_EQUALITY(Sacado::scalarValue(integral_h(0)),4.0,1e-15);
  TEST_FLOATING_EQUALITY(Sacado::scalarValue(integral_h(1)),8.0*std::sqrt(2),1e-15);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(integrator_scalar,test3d,EvalType)
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
  int numCells = 2, numVerts = 8, dim = 3;
  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  MDFieldArrayFactory af("",true);
  workset->cell_node_coordinates = af.buildStaticArray<double,Cell,NODE,Dim>("coords",numCells,numVerts,dim);
  Workset::CellCoordArray coords = workset->cell_node_coordinates;
  auto coords_v = coords.get_static_view();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int) {
      coords_v(0,0,0) = 1.0; coords_v(0,0,1) = 0.0; coords_v(0,0,2) = 0.0;
      coords_v(0,1,0) = 1.0; coords_v(0,1,1) = 1.0; coords_v(0,1,2) = 0.0;
      coords_v(0,2,0) = 0.0; coords_v(0,2,1) = 1.0; coords_v(0,2,2) = 0.0;
      coords_v(0,3,0) = 0.0; coords_v(0,3,1) = 0.0; coords_v(0,3,2) = 0.0;
      coords_v(0,4,0) = 1.0; coords_v(0,4,1) = 0.0; coords_v(0,4,2) = 1.0;
      coords_v(0,5,0) = 1.0; coords_v(0,5,1) = 1.0; coords_v(0,5,2) = 1.0;
      coords_v(0,6,0) = 0.0; coords_v(0,6,1) = 1.0; coords_v(0,6,2) = 1.0;
      coords_v(0,7,0) = 0.0; coords_v(0,7,1) = 0.0; coords_v(0,7,2) = 1.0;

      coords_v(1,0,0) = 0.0; coords_v(1,0,1) = 0.0; coords_v(1,0,2) = 0.0;
      coords_v(1,1,0) =-1.0; coords_v(1,1,1) = 1.0; coords_v(1,1,2) = 0.0;
      coords_v(1,2,0) = 2.0; coords_v(1,2,1) = 2.0; coords_v(1,2,2) = 0.0;
      coords_v(1,3,0) = 1.0; coords_v(1,3,1) = 1.0; coords_v(1,3,2) = 0.0;
      coords_v(1,4,0) = 0.0; coords_v(1,4,1) = 0.0; coords_v(1,4,2) = 2.0;
      coords_v(1,5,0) =-1.0; coords_v(1,5,1) = 1.0; coords_v(1,5,2) = 2.0;
      coords_v(1,6,0) = 2.0; coords_v(1,6,1) = 2.0; coords_v(1,6,2) = 2.0;
      coords_v(1,7,0) = 1.0; coords_v(1,7,1) = 1.0; coords_v(1,7,2) = 2.0;
    });

  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));

  int quadOrder = 5;
  panzer::CellData cellData(2,topo);
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  out << "num quad points = " << quadRule->num_points << std::endl;
  Teuchos::RCP<panzer::IntegrationValues2<double> > quadValues = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
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
  Teuchos::RCP<PHX::MDField<typename EvalType::ScalarT,panzer::Cell> > integralPtr;
  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<panzer::Cell>(quadRule->dl_scalar->extent(0)));

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
  std::vector<PHX::index_size_type> hess_derivative_dimensions;
  hess_derivative_dimensions.push_back(4);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(hess_derivative_dimensions);
#endif

  fm->postRegistrationSetup(setupData);

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<EvalType>(integral);

  integral.print(out,true);
  auto integral_v = integral.get_static_view();
  auto integral_h = Kokkos::create_mirror_view(integral_v);
  Kokkos::deep_copy(integral_h, integral_v);
  TEST_FLOATING_EQUALITY(Sacado::scalarValue(integral_h(0)),2.0,1e-15);
  TEST_FLOATING_EQUALITY(Sacado::scalarValue(integral_h(1)),8.0,1e-15);
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

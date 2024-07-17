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

#include "PanzerDiscFE_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_BasisValues2.hpp"

#include "Panzer_DOF.hpp"
#include "Panzer_DOF_PointField.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "UnitValueEvaluator.hpp"

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(dof_pointfield,value,TYPE)

namespace panzer {

typedef Kokkos::DynRankView<double,PHX::Device> FieldArray;

//**********************************************************************
template<typename EvalT, typename Traits>
class DummyFieldEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DummyFieldEvaluator(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,Cell,panzer::BASIS> fieldValue;
}; // end of class DummyFieldEvaluator

template<typename EvalT, typename Traits>
DummyFieldEvaluator<EvalT, Traits>::
DummyFieldEvaluator(
  const Teuchos::ParameterList& p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  Teuchos::RCP<panzer::PureBasis> basis
     = p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");

  // grab information from quadrature rule
  fieldValue = PHX::MDField<ScalarT,Cell,BASIS>(name, basis->functional);

  this->addEvaluatedField(fieldValue);

  std::string n = "DummyFieldEvaluator: " + name;
  this->setName(n);
}
template<typename EvalT, typename Traits>
void
DummyFieldEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{
  auto lfieldValue=fieldValue;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int ) {
      int i = 0;
      for(int cell=0;cell<lfieldValue.extent_int(0);cell++) {
	for(int pt=0;pt<lfieldValue.extent_int(1);pt++) {
	  lfieldValue(cell,pt) = 1.0+i;
	  ++i;
	}
      }
    });
}
//**********************************************************************
template<typename EvalT, typename Traits>
class RefCoordEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RefCoordEvaluator(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,panzer::Point,panzer::Dim> fieldValue;
  Teuchos::RCP<panzer::IntegrationValues2<double> > quadValues;
public:
  Teuchos::RCP<PHX::DataLayout> coordsLayout;
}; // end of class RefCoordEvaluator

template<typename EvalT, typename Traits>
RefCoordEvaluator<EvalT, Traits>::
RefCoordEvaluator(
  const Teuchos::ParameterList& p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  quadValues = p.get< Teuchos::RCP<panzer::IntegrationValues2<double> > >("Quad Values");

  // grab information from quadrature rule
  coordsLayout = Teuchos::rcp(new PHX::MDALayout<panzer::Point,panzer::Dim>(quadValues->int_rule->num_points,
                                                                         quadValues->int_rule->spatial_dimension));
  fieldValue = PHX::MDField<ScalarT,panzer::Point,panzer::Dim>(name, coordsLayout);

  this->addEvaluatedField(fieldValue);

  std::string n = "RefCoordEvaluator: " + name;
  this->setName(n);
}
template<typename EvalT, typename Traits>
void
RefCoordEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{
  auto lfieldValue = fieldValue;
  auto lcub_points = quadValues->cub_points;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA (int ) {
      for(int cell=0;cell<lfieldValue.extent_int(0);cell++)
	for(int pt=0;pt<lfieldValue.extent_int(1);pt++)
	  lfieldValue(cell,pt) = lcub_points(cell,pt);
    });
}
//**********************************************************************

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(dof_pointfield,value,EvalType)
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

  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

  // build quadrature values
  int quadOrder = 5;
  panzer::CellData cellData(2,1,topo);
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  Teuchos::RCP<panzer::IntegrationValues2<double> > quadValues = Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  // build basis values
  std::string basisName = "Q2";
  Teuchos::RCP<panzer::PureBasis> pureBasis = Teuchos::rcp(new panzer::PureBasis(basisName,2,numCells,topo));
  Teuchos::RCP<panzer::BasisIRLayout> basisLayout = Teuchos::rcp(new panzer::BasisIRLayout(pureBasis,*quadRule));
  Teuchos::RCP<panzer::BasisValues2<double> > basisValues
     = Teuchos::rcp(new panzer::BasisValues2<double>("",true,true));
  basisValues->setupArrays(basisLayout);
  basisValues->evaluateValues(quadValues->cub_points,quadValues->jac,quadValues->jac_det,quadValues->jac_inv,quadValues->weighted_measure,coords);

  // construct workset
  workset->cell_local_ids.push_back(0); workset->cell_local_ids.push_back(1);
  workset->num_cells = numCells;
  workset->block_id = "eblock-0_0";
  workset->ir_degrees = Teuchos::rcp(new std::vector<int>);
  workset->ir_degrees->push_back(quadRule->cubature_degree);
  workset->int_rules.push_back(quadValues);
  workset->basis_names = Teuchos::rcp(new std::vector<std::string>);
  workset->basis_names->push_back(basisLayout->name());
  workset->bases.push_back(basisValues);

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
     = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     // fill the basis with some dummy values
     Teuchos::ParameterList p;
     p.set("Name","TestField");
     p.set("Basis",pureBasis);
     RCP<panzer::DummyFieldEvaluator<EvalType,panzer::Traits> > eval
        = rcp(new panzer::DummyFieldEvaluator<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
  }

  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point> refField
     = PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point>("TestField",quadRule->dl_scalar);
  {
     // fill the basis with some dummy values
     Teuchos::ParameterList p;
     p.set("Name","TestField");
     p.set("Basis",basisLayout);
     p.set("IR",quadRule);
     RCP<panzer::DOF<EvalType,panzer::Traits> > eval
        = rcp(new panzer::DOF<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(*eval->evaluatedFields()[0]);
  }

  Teuchos::RCP<PHX::DataLayout> coordsLayout;
  {
     // build reference coordinate evaluator
     Teuchos::ParameterList p;
     p.set("Name","RefCoord");
     p.set("Quad Values",quadValues);
     RCP<panzer::RefCoordEvaluator<EvalType,panzer::Traits> > eval
        = rcp(new panzer::RefCoordEvaluator<EvalType,panzer::Traits>(p));

     coordsLayout = eval->coordsLayout;
     fm->registerEvaluator<EvalType>(eval);
  }

  // add evaluator under test
  ///////////////////////////////////////////////////

  // test 0
  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point> dofPointField0
     = PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point>("TestFieldpostfix",quadRule->dl_scalar);
  {
     RCP<panzer::DOF_PointField<EvalType,panzer::Traits> > eval
        = rcp(new panzer::DOF_PointField<EvalType,panzer::Traits>("postfix","TestField",*pureBasis,"RefCoord",coordsLayout,quadRule->dl_scalar));

     Teuchos::RCP<PHX::FieldTag> ft = eval->evaluatedFields()[0];
     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(*ft);

     out << "Requiring \"" << *ft << std::endl;
     TEST_EQUALITY(ft->name(),"TestFieldpostfix");
  }

  // test 1
  PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point> dofPointField1
     = PHX::MDField<typename EvalType::ScalarT,panzer::Cell,panzer::Point>("TestFieldRefCoord",quadRule->dl_scalar);
  {
     RCP<panzer::DOF_PointField<EvalType,panzer::Traits> > eval
        = rcp(new panzer::DOF_PointField<EvalType,panzer::Traits>("TestField",*pureBasis,"RefCoord",coordsLayout,quadRule->dl_scalar,true));

     Teuchos::RCP<PHX::FieldTag> ft = eval->evaluatedFields()[0];
     fm->registerEvaluator<EvalType>(eval);
     fm->requireField<EvalType>(*ft);

     out << "Requiring \"" << *ft << std::endl;
     TEST_EQUALITY(ft->name(),"TestFieldRefCoord");
  }

  panzer::Traits::SD setupData;
  {
    auto worksets = rcp(new std::vector<panzer::Workset>);
    worksets->push_back(*workset);
    setupData.worksets_ = worksets;
  }

  std::vector<PHX::index_size_type> derivative_dimensions;
  derivative_dimensions.push_back(8);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(derivative_dimensions);
#endif
  fm->postRegistrationSetup(setupData);
  fm->writeGraphvizFile();

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<EvalType>(refField);
  fm->getFieldData<EvalType>(dofPointField0);
  fm->getFieldData<EvalType>(dofPointField1);

  // check names to make sure they are still correct
  TEST_EQUALITY(refField.fieldTag().name(),"TestField");
  TEST_EQUALITY(dofPointField0.fieldTag().name(),"TestFieldpostfix");
  TEST_EQUALITY(dofPointField1.fieldTag().name(),"TestFieldRefCoord");

  // check sizes
  TEST_EQUALITY(refField.size(),dofPointField0.size());
  TEST_EQUALITY(refField.size(),dofPointField1.size());

  // check the results
  auto refField_h = Kokkos::create_mirror_view(refField.get_view());
  Kokkos::deep_copy(refField_h, refField.get_view());
  auto dofPointField0_h = Kokkos::create_mirror_view(dofPointField0.get_view());
  Kokkos::deep_copy(dofPointField0_h, dofPointField0.get_view());
  auto dofPointField1_h = Kokkos::create_mirror_view(dofPointField1.get_view());
  Kokkos::deep_copy(dofPointField1_h, dofPointField1.get_view());
  for(int cell=0;cell<refField.extent_int(0);cell++) {
    for(int pt=0;pt<refField.extent_int(1);pt++) {
      TEST_FLOATING_EQUALITY(Sacado::scalarValue(refField_h(cell,pt)),Sacado::scalarValue(dofPointField0_h(cell,pt)),1e-15);
      TEST_FLOATING_EQUALITY(Sacado::scalarValue(refField_h(cell,pt)),Sacado::scalarValue(dofPointField1_h(cell,pt)),1e-15);
      // TEST_EQUALITY(refField_h(cell,pt),dofPointField0_h(cell,pt));
      // TEST_EQUALITY(refField_h(cell,pt),dofPointField1_h(cell,pt));
    }
  }
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

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
#include "Teuchos_ScalarTraits.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Panzer_DOF.hpp"
#include "Panzer_DOF_BasisToBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "UnitValueEvaluator.hpp"

// for making explicit instantiated tests easier
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(dof_pointfield,value,TYPE)
  //TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(dof_pointfield,gradient,TYPE)

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
  auto fieldValue_h = Kokkos::create_mirror_view(fieldValue.get_view());
  fieldValue_h(0,0) = 1.0;
  fieldValue_h(0,1) = 2.0;
  fieldValue_h(0,2) = 2.0;
  fieldValue_h(0,3) = 1.0;

  fieldValue_h(1,0) = 2.0;
  fieldValue_h(1,1) = 3.0;
  fieldValue_h(1,2) = 3.0;
  fieldValue_h(1,3) = 2.0;
  Kokkos::deep_copy(fieldValue.get_view(), fieldValue_h);
}

//**********************************************************************

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(dof_pointfield,value,EvalType)
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // build a dummy workset
  //////////////////////////////////////////////////////////
  int numCells = 2;

  // build basis values
  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

  Teuchos::RCP<panzer::PureBasis> sourceBasis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,numCells,topo));
  Teuchos::RCP<panzer::PureBasis> targetBasis = Teuchos::rcp(new panzer::PureBasis("HGrad",2,numCells,topo));

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
     = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);

  // add in some evaluators
  ///////////////////////////////////////////////////
  {
     // fill the basis with some dummy values
     Teuchos::ParameterList p;
     p.set("Name","Pressure");
     p.set("Basis",sourceBasis);
     RCP<panzer::DummyFieldEvaluator<EvalType,panzer::Traits> > e =
       rcp(new panzer::DummyFieldEvaluator<EvalType,panzer::Traits>(p));

     fm->registerEvaluator<EvalType>(e);
  }

  // add evaluator under test
  ///////////////////////////////////////////////////

  {
    RCP<panzer::DOF_BasisToBasis<EvalType,panzer::Traits> > e =
      rcp(new panzer::DOF_BasisToBasis<EvalType,panzer::Traits>("Pressure",*sourceBasis,*targetBasis));

    fm->registerEvaluator<EvalType>(e);

    Teuchos::RCP<PHX::FieldTag> ft = e->evaluatedFields()[0];
    fm->requireField<EvalType>(*ft);
  }

  Teuchos::RCP<panzer::Workset> workset = Teuchos::rcp(new panzer::Workset);
  workset->num_cells = numCells;

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

  //fm->writeGraphvizFile();

  panzer::Traits::PED preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  typedef typename EvalType::ScalarT ScalarT;

  typename PHX::MDField<ScalarT,Cell,BASIS> s("Pressure",sourceBasis->functional);
  typename PHX::MDField<ScalarT,Cell,BASIS> t("Pressure",targetBasis->functional);

  fm->getFieldData<EvalType>(s);
  fm->getFieldData<EvalType>(t);

  typename Teuchos::ScalarTraits<typename Sacado::ScalarType<ScalarT>::type>::magnitudeType tol =
    100.0 * Teuchos::ScalarTraits<typename Sacado::ScalarType<ScalarT>::type>::eps();

  auto s_h = Kokkos::create_mirror_view(s.get_view());
  Kokkos::deep_copy(s_h, s.get_view());
  auto t_h = Kokkos::create_mirror_view(t.get_view());
  Kokkos::deep_copy(t_h, t.get_view());
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,0)),ScalarT(s_h(0,0)),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,1)),ScalarT(s_h(0,1)),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,2)),ScalarT(s_h(0,2)),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,3)),ScalarT(s_h(0,3)),tol);

  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,0)),ScalarT(s_h(1,0)),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,1)),ScalarT(s_h(1,1)),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,2)),ScalarT(s_h(1,2)),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,3)),ScalarT(s_h(1,3)),tol);

  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,4)),ScalarT(1.5),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,5)),ScalarT(2.0),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,6)),ScalarT(1.5),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,7)),ScalarT(1.0),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(0,8)),ScalarT(1.5),tol);

  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,4)),ScalarT(2.5),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,5)),ScalarT(3.0),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,6)),ScalarT(2.5),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,7)),ScalarT(2.0),tol);
  TEST_FLOATING_EQUALITY(ScalarT(t_h(1,8)),ScalarT(2.5),tol);
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

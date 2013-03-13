// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

#include "Panzer_config.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntegrationValues.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_BasisValues.hpp"

#include "Panzer_DOF.hpp"
#include "Panzer_DOF_BasisToBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"

#include "UnitValueEvaluator.hpp"

// for making explicit instantiated tests easier 
#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(dof_pointfield,value,TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(dof_pointfield,gradient,TYPE)

namespace panzer {

typedef Intrepid::FieldContainer<double> FieldArray;

//**********************************************************************
PHX_EVALUATOR_CLASS(DummyFieldEvaluator)
  PHX::MDField<ScalarT,Cell,panzer::BASIS> fieldValue;
PHX_EVALUATOR_CLASS_END
PHX_EVALUATOR_CTOR(DummyFieldEvaluator,p)
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
PHX_POST_REGISTRATION_SETUP(DummyFieldEvaluator,sd,fm)
{
  this->utils.setFieldData(fieldValue,fm);

  
}

PHX_EVALUATE_FIELDS(DummyFieldEvaluator,workset)
{ 
  fieldValue(0,0) = 1.0;
  fieldValue(0,1) = 2.0;
  fieldValue(0,2) = 2.0;
  fieldValue(0,3) = 1.0;
  
  fieldValue(1,0) = 2.0;
  fieldValue(1,1) = 3.0;
  fieldValue(1,2) = 3.0;
  fieldValue(1,3) = 2.0;
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

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);

  fm->postRegistrationSetup(setupData);

  //fm->writeGraphvizFile();

  panzer::GlobalEvaluationDataContainer preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  typedef typename EvalType::ScalarT ScalarT;

  typename PHX::MDField<ScalarT,Cell,BASIS> s("Pressure",sourceBasis->functional);
  typename PHX::MDField<ScalarT,Cell,BASIS> t("Pressure",targetBasis->functional);

  fm->getFieldData<ScalarT,EvalType>(s);
  fm->getFieldData<ScalarT,EvalType>(t);

  typename Teuchos::ScalarTraits<ScalarT>::magnitudeType tol =
    100.0 * Teuchos::ScalarTraits<ScalarT>::eps();

  TEST_FLOATING_EQUALITY(t(0,0),s(0,0),tol);
  TEST_FLOATING_EQUALITY(t(0,1),s(0,1),tol);
  TEST_FLOATING_EQUALITY(t(0,2),s(0,2),tol);
  TEST_FLOATING_EQUALITY(t(0,3),s(0,3),tol);

  TEST_FLOATING_EQUALITY(t(1,0),s(1,0),tol);
  TEST_FLOATING_EQUALITY(t(1,1),s(1,1),tol);
  TEST_FLOATING_EQUALITY(t(1,2),s(1,2),tol);
  TEST_FLOATING_EQUALITY(t(1,3),s(1,3),tol);

  TEST_FLOATING_EQUALITY(t(0,4),ScalarT(1.5),tol);
  TEST_FLOATING_EQUALITY(t(0,5),ScalarT(2.0),tol);
  TEST_FLOATING_EQUALITY(t(0,6),ScalarT(1.5),tol);
  TEST_FLOATING_EQUALITY(t(0,7),ScalarT(1.0),tol);
  TEST_FLOATING_EQUALITY(t(0,8),ScalarT(1.5),tol);

  TEST_FLOATING_EQUALITY(t(1,4),ScalarT(2.5),tol);
  TEST_FLOATING_EQUALITY(t(1,5),ScalarT(3.0),tol);
  TEST_FLOATING_EQUALITY(t(1,6),ScalarT(2.5),tol);
  TEST_FLOATING_EQUALITY(t(1,7),ScalarT(2.0),tol);
  TEST_FLOATING_EQUALITY(t(1,8),ScalarT(2.5),tol);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(dof_pointfield,gradient,EvalType)
{
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

   UNIT_TEST_GROUP(SGResidualType)
   UNIT_TEST_GROUP(SGJacobianType)
#endif

}

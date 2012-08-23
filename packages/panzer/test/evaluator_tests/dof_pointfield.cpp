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
#include "Panzer_DOF_PointField.hpp"
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
{ this->utils.setFieldData(fieldValue,fm); }
PHX_EVALUATE_FIELDS(DummyFieldEvaluator,workset)
{ 
   for(int i=0;i<fieldValue.size();i++)
      fieldValue[i] = 1.0+i;
}
//**********************************************************************
PHX_EVALUATOR_CLASS(RefCoordEvaluator)
  PHX::MDField<ScalarT,panzer::Point,panzer::Dim> fieldValue;
  Teuchos::RCP<panzer::IntegrationValues<double,FieldArray> > quadValues;
public:
  Teuchos::RCP<PHX::DataLayout> coordsLayout;
PHX_EVALUATOR_CLASS_END
PHX_EVALUATOR_CTOR(RefCoordEvaluator,p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  quadValues = p.get< Teuchos::RCP<panzer::IntegrationValues<double,FieldArray> > >("Quad Values");

  // grab information from quadrature rule
  coordsLayout = Teuchos::rcp(new PHX::MDALayout<panzer::Point,panzer::Dim>(quadValues->int_rule->num_points,
                                                                         quadValues->int_rule->spatial_dimension));
  fieldValue = PHX::MDField<ScalarT,panzer::Point,panzer::Dim>(name, coordsLayout);

  this->addEvaluatedField(fieldValue);
  
  std::string n = "RefCoordEvaluator: " + name;
  this->setName(n);
}
PHX_POST_REGISTRATION_SETUP(RefCoordEvaluator,sd,fm)
{ this->utils.setFieldData(fieldValue,fm); }
PHX_EVALUATE_FIELDS(RefCoordEvaluator,workset)
{ 
   for(int i=0;i<fieldValue.size();i++)
      fieldValue[i] = quadValues->cub_points[i];
}
//**********************************************************************

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(dof_pointfield,value,EvalType)
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

  Teuchos::RCP<shards::CellTopology> topo
    = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

  // build quadrature values
  int quadOrder = 5;
  panzer::CellData cellData(2,2,1,topo);
  Teuchos::RCP<panzer::IntegrationRule> quadRule = Teuchos::rcp(new panzer::IntegrationRule(quadOrder,cellData));
  Teuchos::RCP<panzer::IntegrationValues<double,FieldArray> > quadValues = Teuchos::rcp(new panzer::IntegrationValues<double,FieldArray>);
  quadValues->setupArrays(quadRule);
  quadValues->evaluateValues(coords);

  // build basis values
  std::string basisName = "Q2";
  Teuchos::RCP<panzer::PureBasis> pureBasis = Teuchos::rcp(new panzer::PureBasis(basisName,numCells,topo));
  Teuchos::RCP<panzer::BasisIRLayout> basisLayout = Teuchos::rcp(new panzer::BasisIRLayout(pureBasis,*quadRule));
  Teuchos::RCP<panzer::BasisValues<double,Intrepid::FieldContainer<double> > > basisValues 
     = Teuchos::rcp(new panzer::BasisValues<double,Intrepid::FieldContainer<double> >());
  panzer::IntrepidFieldContainerFactory af;
  basisValues->setupArrays(basisLayout,af);
  basisValues->evaluateValues(quadValues->cub_points,quadValues->jac,quadValues->jac_det,quadValues->jac_inv,quadValues->weighted_measure,coords);

  // construct workset
  workset->cell_local_ids.push_back(0); workset->cell_local_ids.push_back(1);
  workset->num_cells = numCells;
  workset->block_id = "eblock-0_0";
  workset->ir_degrees = Teuchos::rcp(new std::vector<int>);
  workset->ir_degrees->push_back(quadRule->cubature_degree);
  workset->int_rules.push_back(quadValues);
  workset->basis_names = Teuchos::rcp(new std::vector<std::string>);
  workset->basis_names->push_back(basisName);
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

  panzer::Traits::SetupData setupData;
  setupData.worksets_ = rcp(new std::vector<panzer::Workset>);
  setupData.worksets_->push_back(*workset);
  fm->postRegistrationSetup(setupData);
  fm->writeGraphvizFile();

  panzer::GlobalEvaluationDataContainer preEvalData;

  fm->preEvaluate<EvalType>(preEvalData);
  fm->evaluateFields<EvalType>(*workset);
  fm->postEvaluate<EvalType>(0);

  fm->getFieldData<typename EvalType::ScalarT,EvalType>(refField);
  fm->getFieldData<typename EvalType::ScalarT,EvalType>(dofPointField0);
  fm->getFieldData<typename EvalType::ScalarT,EvalType>(dofPointField1);

  // check names to make sure they are still correct
  TEST_EQUALITY(refField.fieldTag().name(),"TestField");
  TEST_EQUALITY(dofPointField0.fieldTag().name(),"TestFieldpostfix");
  TEST_EQUALITY(dofPointField1.fieldTag().name(),"TestFieldRefCoord");

  // check sizes
  TEST_EQUALITY(refField.size(),dofPointField0.size());
  TEST_EQUALITY(refField.size(),dofPointField1.size());

  // check the results
  for(int i=0;i<refField.size();i++) {
     TEST_EQUALITY(refField[i],dofPointField0[i]);
     TEST_EQUALITY(refField[i],dofPointField1[i]);
  }
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

   // these are failing with Stokhos enabled! 8/18/2011
   UNIT_TEST_GROUP(SGResidualType)
   UNIT_TEST_GROUP(SGJacobianType)
#endif

}

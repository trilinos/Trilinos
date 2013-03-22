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

#include "Panzer_FieldManagerBuilder.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_PointValues_Evaluator.hpp"
#include "Panzer_BasisValues_Evaluator.hpp"
#include "Panzer_DOF.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

#include "RandomFieldEvaluator.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Epetra_MpiComm.h"

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb, int integration_order);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);
  Teuchos::RCP<panzer::IntegrationRule> buildIR(std::size_t worksetSize,int cubature_degree);

  TEUCHOS_UNIT_TEST(point_values_evaluator, eval)
  {
    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    panzer::CellData cell_data(workset_size, mesh->getCellTopology("eblock-0_0"));

    const int integration_order = 1;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    const int default_int_order = integration_order;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,*physicsBlock);
    TEST_EQUALITY(work_sets->size(),1);

    int num_points = 3;
    RCP<const panzer::PointRule> point_rule = rcp(new panzer::PointRule("RandomPoints",num_points, cell_data));

    Teuchos::RCP<Intrepid::FieldContainer<double> > userArray 
       = Teuchos::rcp(new Intrepid::FieldContainer<double>(num_points,2));
    Intrepid::FieldContainer<double> & point_coordinates = *userArray;
    point_coordinates(0,0) =  0.0; point_coordinates(0,1) = 0.0; // mid point
    point_coordinates(1,0) =  0.5; point_coordinates(1,1) = 0.5; // mid point of upper left quadrant
    point_coordinates(2,0) = -0.5; point_coordinates(2,1) = 0.0; // mid point of line from center to left side
    

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<const std::vector<Teuchos::RCP<PHX::FieldTag > > > evalResFields;
    Teuchos::RCP<const std::vector<Teuchos::RCP<PHX::FieldTag > > > evalJacFields;
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Residual,panzer::Traits>(point_rule,*userArray));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),6);
       evalResFields = Teuchos::rcpFromRef(evaluator->evaluatedFields());

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,*userArray));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),6);
       evalJacFields = Teuchos::rcpFromRef(evaluator->evaluatedFields());

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SetupData sd;
    fm.postRegistrationSetup(sd);
    fm.print(out);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.ghostedLinContainer = Teuchos::null;
    workset.linContainer = Teuchos::null;
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);
    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    PHX::MDField<panzer::Traits::Residual::ScalarT> 
       point_coords(point_rule->getName()+"_"+"point_coords",point_rule->dl_vector);
    fm.getFieldData<panzer::Traits::Residual::ScalarT,panzer::Traits::Residual>(point_coords);

    typedef panzer::ArrayTraits<double,Intrepid::FieldContainer<double> >::size_type size_type;

    for(size_type c=0;c<basis_q1->numCells();c++) {
       double dx = 0.5;
       double dy = 0.5;
       for(size_type p=0;p<num_points;p++) {
          double x = dx*(point_coordinates(p,0)+1.0)/2.0 + workset.cell_vertex_coordinates(c,0,0); 
          double y = dy*(point_coordinates(p,1)+1.0)/2.0 + workset.cell_vertex_coordinates(c,0,1);
          TEST_FLOATING_EQUALITY(point_coords(c,p,0),x,1e-10);
          TEST_FLOATING_EQUALITY(point_coords(c,p,1),y,1e-10);
       }
    }
  }

  TEUCHOS_UNIT_TEST(basis_values_evaluator, eval)
  {
    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    panzer::CellData cell_data(basis_q1->numCells(), basis_q1->getCellTopology());

    int integration_order = 4;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    const int default_int_order = integration_order;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,*physicsBlock);
    panzer::Workset & workset = (*work_sets)[0];
    TEST_EQUALITY(work_sets->size(),1);

    Teuchos::RCP<panzer::IntegrationRule> point_rule = buildIR(workset_size,integration_order);
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > int_values;
    int_values.setupArrays(point_rule);
    int_values.evaluateValues(workset.cell_vertex_coordinates);

    Teuchos::RCP<Intrepid::FieldContainer<double> > userArray = Teuchos::rcpFromRef(int_values.cub_points);

    Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis_q1,*point_rule));

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    Teuchos::RCP<const std::vector<Teuchos::RCP<PHX::FieldTag > > > evalJacFields;
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,*userArray));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::BasisValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,basis_q1));

       TEST_EQUALITY(evaluator->evaluatedFields().size(),4);
       evalJacFields = Teuchos::rcpFromRef(evaluator->evaluatedFields());

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]);
    }

    panzer::Traits::SetupData sd;
    fm.postRegistrationSetup(sd);
    fm.print(out);

    // run tests
    /////////////////////////////////////////////////////////////

    workset.ghostedLinContainer = Teuchos::null;
    workset.linContainer = Teuchos::null;
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    PHX::MDField<panzer::Traits::Jacobian::ScalarT> 
       basis(basis_q1->name()+"_"+point_rule->getName()+"_"+"basis",layout->basis);
    fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(basis);
    out << basis << std::endl;

    std::size_t basisIndex = panzer::getBasisIndex(layout->name(), workset);
    Teuchos::RCP<panzer::BasisValues<double,Intrepid::FieldContainer<double> > > bases = workset.bases[basisIndex];
    TEST_ASSERT(bases!=Teuchos::null);
    TEST_EQUALITY(bases->basis.size(),basis.size());
    for(int i=0;i<basis.size();i++) {
       panzer::Traits::Jacobian::ScalarT truth = bases->basis[i];
       TEST_FLOATING_EQUALITY(truth,basis[i],1e-10);
    }
  }

  TEUCHOS_UNIT_TEST(dof_point_values_evaluator, eval)
  {
    const std::size_t workset_size = 4;
    const std::string fieldName_q1 = "U";
    const std::string fieldName_qedge1 = "V";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(2,2);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    panzer::CellData cell_data(basis_q1->numCells(), basis_q1->getCellTopology());

    int integration_order = 4;
    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb,integration_order);

    const int default_int_order = integration_order;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,*physicsBlock);
    panzer::Workset & workset = (*work_sets)[0];
    TEST_EQUALITY(work_sets->size(),1);

    Teuchos::RCP<panzer::IntegrationRule> point_rule = buildIR(workset_size,integration_order);
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > int_values;
    int_values.setupArrays(point_rule);
    int_values.evaluateValues(workset.cell_vertex_coordinates);

    Teuchos::RCP<Intrepid::FieldContainer<double> > userArray = Teuchos::rcpFromRef(int_values.cub_points);

    Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis_q1,*point_rule));

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new RandomFieldEvaluator<panzer::Traits::Jacobian,panzer::Traits>("TEMPERATURE",basis_q1->functional));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::PointValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,*userArray));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::BasisValues_Evaluator<panzer::Traits::Jacobian,panzer::Traits>(point_rule,basis_q1));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
    }
    {
       Teuchos::ParameterList p;
       p.set("Name","TEMPERATURE");
       p.set("Basis",layout);
       p.set("IR",point_rule);
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::DOF<panzer::Traits::Jacobian,panzer::Traits>(p));
       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]); // require DOF
    }
    {
       Teuchos::ParameterList p;
       p.set("Name","TEMPERATURE");
       p.set("Basis",layout->getBasis().getConst());
       p.set("Point Rule",Teuchos::rcp_static_cast<panzer::PointRule>(point_rule).getConst()); 
       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::DOF_PointValues<panzer::Traits::Jacobian,panzer::Traits>(p));

       fm.registerEvaluator<panzer::Traits::Jacobian>(evaluator);
       fm.requireField<panzer::Traits::Jacobian>(*evaluator->evaluatedFields()[0]); // require DOF_PointValues
    }

    panzer::Traits::SetupData sd;
    sd.worksets_ = work_sets;
    fm.postRegistrationSetup(sd);
    fm.print(out);

    // run tests
    /////////////////////////////////////////////////////////////

    workset.ghostedLinContainer = Teuchos::null;
    workset.linContainer = Teuchos::null;
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Jacobian>(workset);

    PHX::MDField<panzer::Traits::Jacobian::ScalarT> ref_field("TEMPERATURE",point_rule->dl_scalar);
    PHX::MDField<panzer::Traits::Jacobian::ScalarT> fut_field("TEMPERATURE_"+point_rule->getName(),point_rule->dl_scalar); // "field under test"
    fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(ref_field);
    fm.getFieldData<panzer::Traits::Jacobian::ScalarT,panzer::Traits::Jacobian>(fut_field);

    TEST_EQUALITY(ref_field.size(),fut_field.size());
    for(int i=0;i<ref_field.size();i++) {
       TEST_FLOATING_EQUALITY(fut_field[i],ref_field[i],1e-10);
    }
  }

  Teuchos::RCP<panzer::PureBasis> buildBasis(std::size_t worksetSize,const std::string & basisName)
  { 
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,topo);
     return Teuchos::rcp(new panzer::PureBasis(basisName,1,cellData)); 
  }

  Teuchos::RCP<panzer::IntegrationRule> buildIR(std::size_t workset_size,int cubature_degree)
  {
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     const panzer::CellData cell_data(workset_size, topo);

     return Teuchos::rcp(new panzer::IntegrationRule(cubature_degree, cell_data));
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
    typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;
    typedef panzer_stk::STK_Interface::VectorFieldType CoordinateField;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",1);
    pl->set("Y Blocks",1);
    pl->set("X Elements",elemX);
    pl->set("Y Elements",elemY);
    
    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD); 

    return mesh;
  }

  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb, int integration_order)
  {
    // Physics block
    ipb->setName("test physics");
    {
      Teuchos::ParameterList& p = ipb->sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",integration_order);
    }
    {
      Teuchos::ParameterList& p = ipb->sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","solid");
      p.set("Basis Type","HCurl");
      p.set("Basis Order",1);
      p.set("Integration Order",integration_order);
    }
    
  }

}

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
#include "Panzer_DOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_WorksetContainer.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Epetra_MpiComm.h"

#include "PointEvaluator.hpp"

#include "user_app_EquationSetFactory.hpp"

#include <cstdio> // for get char
#include <vector>
#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

  Teuchos::RCP<panzer::PureBasis> buildBasis(const std::size_t& worksetSize,const std::string & basisName);
  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb);
  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY);
  Teuchos::RCP<panzer::IntegrationRule> buildIR(const std::size_t& worksetSize,const int& cubature_degree);

  class BilinearPointEvaluator : public PointEvaluation<panzer::Traits::Residual::ScalarT> {
  public:
    virtual ~BilinearPointEvaluator() {}
    virtual void evaluateContainer(const PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::IP,panzer::Dim> & points,
                                   PHX::MDField<panzer::Traits::Residual::ScalarT> & field) const
    {
       int num_cells = field.extent(0);
       int num_qp = points.extent(1);

       for(int i=0;i<num_cells;i++) {
          for(int j=0;j<num_qp;j++) {
             double x = points(i,j,0); // just x and y
             double y = points(i,j,1);

             field(i,j,0) = (x+y)*(x+y);
             field(i,j,1) = sin(x+y);
          }
       }
    }
  };

  TEUCHOS_UNIT_TEST(basis_time_vector, residual)
  {

    const std::size_t workset_size = 1;
    const std::string fieldName_q1 = "TEMPERATURE";
    const std::string fieldName_qedge1 = "ION_TEMPERATURE";

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(1,1);

    // build input physics block
    Teuchos::RCP<panzer::PureBasis> basis_q1 = buildBasis(workset_size,"Q1");
    Teuchos::RCP<panzer::PureBasis> basis_qedge1 = buildBasis(workset_size,"QEdge1");

    Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
    testInitialization(ipb);

    const int default_int_order = 1;
    std::string eBlockID = "eblock-0_0";    
    Teuchos::RCP<user_app::MyFactory> eqset_factory = Teuchos::rcp(new user_app::MyFactory);
    panzer::CellData cellData(workset_size,mesh->getCellTopology("eblock-0_0"));
    Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
    Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = 
      Teuchos::rcp(new PhysicsBlock(ipb,eBlockID,default_int_order,cellData,eqset_factory,gd,false));

    Teuchos::RCP<panzer::IntegrationRule> ir = buildIR(workset_size,4);
    Teuchos::RCP<panzer::BasisIRLayout> layout_qedge1 = Teuchos::rcp(new panzer::BasisIRLayout(basis_qedge1,*ir));

    // build connection manager and field manager
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));
    RCP<panzer::DOFManager<int,int> > dofManager = Teuchos::rcp(new panzer::DOFManager<int,int>(conn_manager,MPI_COMM_WORLD));
    dofManager->addField(fieldName_q1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_q1->getIntrepid2Basis())));
    dofManager->addField(fieldName_qedge1,Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_qedge1->getIntrepid2Basis())));
    dofManager->setOrientationsRequired(true);
    dofManager->buildGlobalUnknowns();

    // build worksets
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks;
    physicsBlocks.push_back(physicsBlock); // pushing back singular

    panzer::WorksetContainer wkstContainer; // (Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)),physicsBlocks,workset_size);

    wkstContainer.setFactory(Teuchos::rcp(new panzer_stk::WorksetFactory(mesh)));
    for(size_t i=0;i<physicsBlocks.size();i++) 
      wkstContainer.setNeeds(physicsBlocks[i]->elementBlockID(),physicsBlocks[i]->getWorksetNeeds());
    wkstContainer.setGlobalIndexer(dofManager);
    wkstContainer.setWorksetSize(workset_size);

    const panzer::WorksetDescriptor wd = panzer::blockDescriptor(eBlockID);
    Teuchos::RCP<std::vector<panzer::Workset> > work_sets = wkstContainer.getWorksets(wd);
    TEST_EQUALITY(work_sets->size(),1);

    // setup field manager, add evaluator under test
    /////////////////////////////////////////////////////////////
 
    PHX::FieldManager<panzer::Traits> fm;

    {
       Teuchos::ParameterList pl;
       pl.set("Name","Integrand");
       pl.set("IR",ir);
       pl.set("Is Vector",true);
       pl.set<Teuchos::RCP<const PointEvaluation<panzer::Traits::Residual::ScalarT> > >("Point Evaluator",
                                                                                        Teuchos::rcp(new BilinearPointEvaluator));

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new PointEvaluator<panzer::Traits::Residual,panzer::Traits>(pl));

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
    }

    {
       Teuchos::ParameterList pl;
       pl.set("Residual Name","Residual");
       pl.set("Value Name","Integrand");
//       pl.set("Test Field Name",fieldName_qedge1);
       pl.set("Basis",layout_qedge1);
       pl.set("IR",ir);
       pl.set<double>("Multiplier", 1.0);
       Teuchos::RCP<const std::vector<std::string> > vec
          = Teuchos::rcp(new std::vector<std::string>);
       pl.set("Field Multipliers", vec);

       Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
          = Teuchos::rcp(new panzer::Integrator_BasisTimesVector<panzer::Traits::Residual,panzer::Traits>(pl));

       fm.registerEvaluator<panzer::Traits::Residual>(evaluator);
       fm.requireField<panzer::Traits::Residual>(*evaluator->evaluatedFields()[0]);
    }

    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(8);
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

    panzer::Traits::SD sd;
    sd.worksets_ = work_sets;
    fm.postRegistrationSetup(sd);

    // run tests
    /////////////////////////////////////////////////////////////

    panzer::Workset & workset = (*work_sets)[0];
    workset.alpha = 0.0;
    workset.beta = 0.0;
    workset.time = 0.0;
    workset.evaluate_transient_terms = false;

    fm.evaluateFields<panzer::Traits::Residual>(workset);

    PHX::MDField<panzer::Traits::Residual::ScalarT,panzer::Cell,panzer::BASIS> 
       fieldData_qedge1("Residual",basis_qedge1->functional);

    fm.getFieldData<panzer::Traits::Residual>(fieldData_qedge1);

    TEST_EQUALITY(fieldData_qedge1.extent(0),1);
    TEST_EQUALITY(fieldData_qedge1.extent(1),4);

    // Transformation is [x,y] = F[x_ref,y_ref] = 0.5*[1,1]+0.5*[1,0;0,1]*[x_ref,y_ref]
    // therefore transformation matrix is DF^{-T} = 2*[1,0;0,1]
    // so curl vector u_ref:Ref_coord=>Ref_Vec transforms with 
    //
    //           u(x,y)=DF^{-T}*u_ref(F^{-1}(x,y))

    TEST_FLOATING_EQUALITY(fieldData_qedge1(0,0),5.0/12.0,1e-5);        // 0 edge basis is [(1-y_ref)/4, 0] 
    TEST_FLOATING_EQUALITY(fieldData_qedge1(0,2),3.0/4.0,1e-5);         // 2 edge basis is [(1+y_ref)/4, 0] 

    // these two have sign changes because of the mesh topology!
    TEST_FLOATING_EQUALITY(fieldData_qedge1(0,1),0.428925006266,1e-5);  // 1 edge basis is [(1+x_ref)/4, 0] 
    TEST_FLOATING_EQUALITY(fieldData_qedge1(0,3),0.344719536524,1e-5);  // 3 edge basis is [(1-x_ref)/4, 0] 
  }

  Teuchos::RCP<panzer::IntegrationRule> buildIR(const std::size_t& workset_size, const int& cubature_degree)
  {
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     const panzer::CellData cell_data(workset_size,topo);

     return Teuchos::rcp(new panzer::IntegrationRule(cubature_degree, cell_data));
  }

  Teuchos::RCP<panzer::PureBasis> buildBasis(const std::size_t& worksetSize,
					     const std::string& basisName)
  { 
     Teuchos::RCP<shards::CellTopology> topo = 
        Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));

     panzer::CellData cellData(worksetSize,topo);
     // hard coded to first order
     return Teuchos::rcp(new panzer::PureBasis(basisName,1,cellData)); 
  }

  Teuchos::RCP<panzer_stk::STK_Interface> buildMesh(int elemX,int elemY)
  {
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

  void testInitialization(const Teuchos::RCP<Teuchos::ParameterList>& ipb)
  {
    // Physics block
    ipb->setName("test physics");
    {
      Teuchos::ParameterList& p = ipb->sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",2);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = ipb->sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","solid");
      p.set("Basis Type","HCurl");
      p.set("Basis Order",1);
      p.set("Integration Order",4);
    }
    
  }

}

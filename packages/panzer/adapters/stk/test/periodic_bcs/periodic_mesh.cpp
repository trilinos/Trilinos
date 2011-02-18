#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"

typedef Intrepid::FieldContainer<double> FieldContainer;

namespace panzer {

  template <typename IntrepidType>
  RCP<const panzer::FieldPattern> buildFieldPattern()
  {
     // build a geometric pattern from a single basis
     RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
     RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
     return pattern;
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, add_get_vector)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::CoordMatcher y_matcher(1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher));

    std::vector<RCP<const panzer_stk::PeriodicBC_MatcherBase> > & mod_vec = mesh->getPeriodicBCVector();
    TEST_EQUALITY(mod_vec.size(),2);

    const std::vector<RCP<const panzer_stk::PeriodicBC_MatcherBase> > & const_vec = mesh.getConst()->getPeriodicBCVector();
    TEST_EQUALITY(const_vec.size(),2);
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    int myRank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::CoordMatcher y_matcher(1);
    // mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager<int,int> > connMngr 
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();
    connMngr->buildConnectivity(*fp);

    const int * conn0 = connMngr->getConnectivity(0);
    const int * conn1 = connMngr->getConnectivity(1);

    if(myRank==0) {
       TEST_EQUALITY(conn0[0],4); // this is periodic
       TEST_EQUALITY(conn0[1],1);
       TEST_EQUALITY(conn0[2],6);
       TEST_EQUALITY(conn0[3],9); // this is periodic
   
       TEST_EQUALITY(conn1[0],2);
       TEST_EQUALITY(conn1[1],3);
       TEST_EQUALITY(conn1[2],8);
       TEST_EQUALITY(conn1[3],7);
    }
    else {
       TEST_EQUALITY(conn0[0],1);
       TEST_EQUALITY(conn0[1],2);
       TEST_EQUALITY(conn0[2],7);
       TEST_EQUALITY(conn0[3],6);
   
       TEST_EQUALITY(conn1[0],3);
       TEST_EQUALITY(conn1[1],4);
       TEST_EQUALITY(conn1[2],9);
       TEST_EQUALITY(conn1[3],8);

    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager2)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    int myRank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::CoordMatcher y_matcher(1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager<int,int> > connMngr 
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();
    connMngr->buildConnectivity(*fp);

    const int * conn0 = connMngr->getConnectivity(0);
    const int * conn1 = connMngr->getConnectivity(1);

    if(myRank==0) {
       TEST_EQUALITY(conn0[0],0);
       TEST_EQUALITY(conn0[1],1);
       TEST_EQUALITY(conn0[2],1);
       TEST_EQUALITY(conn0[3],0);
   
       TEST_EQUALITY(conn1[0],2);
       TEST_EQUALITY(conn1[1],3);
       TEST_EQUALITY(conn1[2],3);
       TEST_EQUALITY(conn1[3],2);
    }
    else {
       TEST_EQUALITY(conn0[0],1);
       TEST_EQUALITY(conn0[1],2);
       TEST_EQUALITY(conn0[2],2);
       TEST_EQUALITY(conn0[3],1);
   
       TEST_EQUALITY(conn1[0],3);
       TEST_EQUALITY(conn1[1],4);
       TEST_EQUALITY(conn1[2],4);
       TEST_EQUALITY(conn1[3],3);

    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager3)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    int myRank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::CoordMatcher y_matcher(1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager<int,int> > connMngr 
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();
    connMngr->buildConnectivity(*fp);

    const int * conn0 = connMngr->getConnectivity(0);
    const int * conn1 = connMngr->getConnectivity(1);

    if(myRank==0) {
       TEST_EQUALITY(conn0[0],4); // from left/right periodicity
       TEST_EQUALITY(conn0[1],1);
       TEST_EQUALITY(conn0[2],1);
       TEST_EQUALITY(conn0[3],4); // from left/right periodicity
   
       TEST_EQUALITY(conn1[0],2);
       TEST_EQUALITY(conn1[1],3);
       TEST_EQUALITY(conn1[2],3);
       TEST_EQUALITY(conn1[3],2);
    }
    else {
       TEST_EQUALITY(conn0[0],1);
       TEST_EQUALITY(conn0[1],2);
       TEST_EQUALITY(conn0[2],2);
       TEST_EQUALITY(conn0[3],1);
   
       TEST_EQUALITY(conn1[0],3);
       TEST_EQUALITY(conn1[1],4);
       TEST_EQUALITY(conn1[2],4);
       TEST_EQUALITY(conn1[3],3);

    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager4)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    int myRank = Comm.MyPID(); 

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::CoordMatcher x_matcher(0);
    panzer_stk::CoordMatcher y_matcher(1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("right","left",y_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager<int,int> > connMngr 
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();
    connMngr->buildConnectivity(*fp);

    const int * conn0 = connMngr->getConnectivity(0);
    const int * conn1 = connMngr->getConnectivity(1);

    if(myRank==0) {
       TEST_EQUALITY(conn0[0],0); // from left/right periodicity
       TEST_EQUALITY(conn0[1],1);
       TEST_EQUALITY(conn0[2],1);
       TEST_EQUALITY(conn0[3],0); // from left/right periodicity
   
       TEST_EQUALITY(conn1[0],2);
       TEST_EQUALITY(conn1[1],3);
       TEST_EQUALITY(conn1[2],3);
       TEST_EQUALITY(conn1[3],2);
    }
    else {
       TEST_EQUALITY(conn0[0],1);
       TEST_EQUALITY(conn0[1],2);
       TEST_EQUALITY(conn0[2],2);
       TEST_EQUALITY(conn0[3],1);
   
       TEST_EQUALITY(conn1[0],3);
       TEST_EQUALITY(conn1[1],0);
       TEST_EQUALITY(conn1[2],0);
       TEST_EQUALITY(conn1[3],3);

    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_hex)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    int myRank = Comm.MyPID(); 

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",1);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    mesh->writeToExodus("what.exo");

    panzer_stk::PlaneMatcher top_matcher(0,2);
    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",top_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",side_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager<int,int> > connMngr 
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid::Basis_HGRAD_HEX_C1_FEM<double,FieldContainer> >();
    connMngr->buildConnectivity(*fp);

    const int * conn0 = connMngr->getConnectivity(0);
    const int * conn1 = connMngr->getConnectivity(1);
    const int * conn2 = connMngr->getConnectivity(2);
    const int * conn3 = connMngr->getConnectivity(3);

    if(myRank==0) {
       TEST_EQUALITY(conn0[0],2);
       TEST_EQUALITY(conn0[1],1);
       TEST_EQUALITY(conn0[2],4);
       TEST_EQUALITY(conn0[3],5);
       TEST_EQUALITY(conn0[4],11);
       TEST_EQUALITY(conn0[5],10);
       TEST_EQUALITY(conn0[6],13);
       TEST_EQUALITY(conn0[7],14);

       TEST_EQUALITY(conn1[0],11);
       TEST_EQUALITY(conn1[1],10);
       TEST_EQUALITY(conn1[2],13);
       TEST_EQUALITY(conn1[3],14);
       TEST_EQUALITY(conn1[4],2);
       TEST_EQUALITY(conn1[5],1);
       TEST_EQUALITY(conn1[6],4);
       TEST_EQUALITY(conn1[7],5);

       TEST_EQUALITY(conn2[0],5);
       TEST_EQUALITY(conn2[1],4);
       TEST_EQUALITY(conn2[2],1);
       TEST_EQUALITY(conn2[3],2);
       TEST_EQUALITY(conn2[4],14);
       TEST_EQUALITY(conn2[5],13);
       TEST_EQUALITY(conn2[6],10);
       TEST_EQUALITY(conn2[7],11);

       TEST_EQUALITY(conn3[0],14);
       TEST_EQUALITY(conn3[1],13);
       TEST_EQUALITY(conn3[2],10);
       TEST_EQUALITY(conn3[3],11);
       TEST_EQUALITY(conn3[4],5);
       TEST_EQUALITY(conn3[5],4);
       TEST_EQUALITY(conn3[6],1);
       TEST_EQUALITY(conn3[7],2);
    }
    else {
       TEST_EQUALITY(conn0[0],1);
       TEST_EQUALITY(conn0[1],2);
       TEST_EQUALITY(conn0[2],5);
       TEST_EQUALITY(conn0[3],4);
       TEST_EQUALITY(conn0[4],10);
       TEST_EQUALITY(conn0[5],11);
       TEST_EQUALITY(conn0[6],14);
       TEST_EQUALITY(conn0[7],13);
       // passed

       TEST_EQUALITY(conn1[0],10);
       TEST_EQUALITY(conn1[1],11);
       TEST_EQUALITY(conn1[2],14);
       TEST_EQUALITY(conn1[3],13);
       TEST_EQUALITY(conn1[4],1);
       TEST_EQUALITY(conn1[5],2);
       TEST_EQUALITY(conn1[6],5);
       TEST_EQUALITY(conn1[7],4);

       TEST_EQUALITY(conn2[0],4);
       TEST_EQUALITY(conn2[1],5);
       TEST_EQUALITY(conn2[2],2);
       TEST_EQUALITY(conn2[3],1);
       TEST_EQUALITY(conn2[4],13);
       TEST_EQUALITY(conn2[5],14);
       TEST_EQUALITY(conn2[6],11);
       TEST_EQUALITY(conn2[7],10);

       TEST_EQUALITY(conn3[0],13);
       TEST_EQUALITY(conn3[1],14);
       TEST_EQUALITY(conn3[2],11);
       TEST_EQUALITY(conn3[3],10);
       TEST_EQUALITY(conn3[4],4);
       TEST_EQUALITY(conn3[5],5);
       TEST_EQUALITY(conn3[6],2);
       TEST_EQUALITY(conn3[7],1);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, getSideIdsAndCoords)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // mesh->writeToExodus("output.exo");

    panzer_stk::CoordMatcher x_matcher(0);
    std::pair<RCP<std::vector<std::size_t> >,
         RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top");

    TEST_EQUALITY(idsAndCoords.first->size(),5);
    TEST_EQUALITY(idsAndCoords.second->size(),5);
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, PeriodicBC_Matcher)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    {
       panzer_stk::CoordMatcher matcher(1);
       Teuchos::RCP<const panzer_stk::PeriodicBC_MatcherBase> pMatch 
             = panzer_stk::buildPeriodicBC_Matcher("left","right",matcher);

       RCP<std::vector<std::pair<std::size_t,std::size_t> > > globallyMatchedIds = pMatch->getMatchedPair(*mesh);

       // for testing purposes!
       RCP<std::vector<std::size_t> > locallyRequiredIds = panzer_stk::periodic_helpers::getLocalSideIds(*mesh,"left");

       TEST_EQUALITY(globallyMatchedIds->size(),locallyRequiredIds->size()); 

       // match top & bottom sides
       for(std::size_t i=0;i<globallyMatchedIds->size();i++) {
          std::pair<std::size_t,std::size_t> pair = (*globallyMatchedIds)[i];
          TEST_EQUALITY(pair.first,pair.second-4);
       }
    }
  }

}

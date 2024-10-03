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

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"
#include "Panzer_STK_PeriodicBC_Parser.hpp"
#include "Panzer_STK_PeriodicBC_MatchConditions.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"

#include <string>

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer {

  template <typename Intrepid2Type>
  RCP<const panzer::FieldPattern> buildFieldPattern()
  {
     // build a geometric pattern from a single basis
     RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
     RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
     return pattern;
  }

  TEUCHOS_UNIT_TEST(periodic_parser, test_obj)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using namespace panzer_stk;

    panzer_stk::PeriodicBC_Parser parser;

    // test basic bc construction functionality
    {
       std::string matcher, bndry1, bndry2;

       parser.buildMatcher_Tokenize("x-coord left;right",matcher,bndry1,bndry2);
       TEST_EQUALITY(matcher,"x-coord"); TEST_EQUALITY(bndry1,"left"); TEST_EQUALITY(bndry2,"right");

       parser.buildMatcher_Tokenize("x-edge left;right",matcher,bndry1,bndry2);
       TEST_EQUALITY(matcher,"x-edge"); TEST_EQUALITY(bndry1,"left"); TEST_EQUALITY(bndry2,"right");

       parser.buildMatcher_Tokenize("y-coord beg ; what ",matcher,bndry1,bndry2);
       TEST_EQUALITY(matcher,"y-coord"); TEST_EQUALITY(bndry1,"beg"); TEST_EQUALITY(bndry2,"what");

       parser.buildMatcher_Tokenize("y-edge beg ; what ",matcher,bndry1,bndry2);
       TEST_EQUALITY(matcher,"y-edge"); TEST_EQUALITY(bndry1,"beg"); TEST_EQUALITY(bndry2,"what");

       RCP<const panzer_stk::PeriodicBC_MatcherBase> matcher_obj;

       matcher_obj = parser.buildMatcher("x-coord left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("xy-coord left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<PlaneMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("(xy)z-quarter-coord left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<QuarterPlaneMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("x-edge left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<QuarterPlaneMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("xy-edge left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<QuarterPlaneMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("xy-face left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<QuarterPlaneMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("x-coord 1e-8: left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("xy-edge 1e-8: left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("x-coord 1e-8, relative: left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(matcher_obj));

       matcher_obj = parser.buildMatcher("xy-edge 1e-8, relative: left;right");
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(matcher_obj));

       TEST_THROW(parser.buildMatcher("dog-coord left;right"),std::logic_error);

       TEST_THROW(parser.buildMatcher("dog-edge left;right"),std::logic_error);

       TEST_THROW(parser.buildMatcher("dog-face left;right"),std::logic_error);

       TEST_THROW(parser.buildMatcher("x-face left;right"),std::logic_error);

       TEST_THROW(parser.buildMatcher("xface left;right"),std::logic_error);

       TEST_THROW(parser.buildMatcher("x-coord 1e-8, misspelled-relatiive: left;right"),std::logic_error);
    }

    // test parameter list based construction
    {
       RCP<Teuchos::ParameterList> pl = Teuchos::rcp(new Teuchos::ParameterList("top_list"));

       pl->set("Count",8);
       pl->set("Periodic Condition 1","y-coord left;right");
       pl->set("Periodic Condition 2","x-coord top;bottom");
       pl->set("Periodic Condition 3","yz-coord fake_a;fake_b");
       pl->set("Periodic Condition 4","(zy)x-quarter-coord fake_a;fake_b");
       pl->set("Periodic Condition 5","y-edge left;right");
       pl->set("Periodic Condition 6","x-edge top;bottom");
       pl->set("Periodic Condition 7","yz-edge fake_a;fake_b");
       pl->set("Periodic Condition 8","yz-face fake_a;fake_b");
       parser.setParameterList(pl);

       TEST_EQUALITY(parser.getMatchers().size(),8);
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(parser.getMatchers()[0]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(parser.getMatchers()[1]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<PlaneMatcher> >(parser.getMatchers()[2]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<QuarterPlaneMatcher> >(parser.getMatchers()[3]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(parser.getMatchers()[4]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<CoordMatcher> >(parser.getMatchers()[5]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<PlaneMatcher> >(parser.getMatchers()[6]));
       TEST_NOTHROW(rcp_dynamic_cast<const PeriodicBC_Matcher<PlaneMatcher> >(parser.getMatchers()[7]));
    }
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
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher,"edge"));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher,"edge"));

    std::vector<RCP<const panzer_stk::PeriodicBC_MatcherBase> > & mod_vec = mesh->getPeriodicBCVector();
    TEST_EQUALITY(mod_vec.size(),4);

    const std::vector<RCP<const panzer_stk::PeriodicBC_MatcherBase> > & const_vec = mesh.getConst()->getPeriodicBCVector();
    TEST_EQUALITY(const_vec.size(),4);
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
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn1 = connMngr->getConnectivity(1);

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

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_edge)
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
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher,"edge"));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn3 = connMngr->getConnectivity(3);

    // test that the left-most edge has the same number as the right-most edge
    if(myRank==0) {
       TEST_EQUALITY(conn0[3],conn3[1]);
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
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn1 = connMngr->getConnectivity(1);

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

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager2_edge)
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
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher,"edge"));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn1 = connMngr->getConnectivity(1);

    // test that top edge has same number as bottom edge for each element
    if(myRank==0) {
       TEST_EQUALITY(conn0[0],conn0[2]);
       TEST_EQUALITY(conn1[0],conn1[2]);
    } else {
       TEST_EQUALITY(conn0[0],conn0[2]);
       TEST_EQUALITY(conn1[0],conn1[2]);
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

    mesh->writeToExodus("twod.exo");

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn1 = connMngr->getConnectivity(1);

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

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager3_edge)
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
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher,"edge"));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher,"edge"));

    mesh->writeToExodus("twod.exo");

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn1 = connMngr->getConnectivity(1);
    const auto * conn3 = connMngr->getConnectivity(3);

    if(myRank==0) {
       TEST_EQUALITY(conn0[3],conn3[1]); // from left/right periodicity
       TEST_EQUALITY(conn0[0],conn0[2]); // from top/bottom periodicity
       TEST_EQUALITY(conn1[0],conn1[2]); // from top/bottom periodicity
    } else {
       TEST_EQUALITY(conn0[0],conn0[2]); // from top/bottom periodicity
       TEST_EQUALITY(conn1[0],conn1[2]); // from top/bottom periodicity
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
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    const auto * conn0 = connMngr->getConnectivity(0);
    const auto * conn1 = connMngr->getConnectivity(1);

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
    // mesh->writeToExodus("what.exo");

    panzer_stk::PlaneMatcher top_matcher(0,2);
    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",top_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",side_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    if(myRank==0) {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(1));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(5));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(3));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(7));

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
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(2));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(6));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(4));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(8));

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

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_hex_edge)
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
    // mesh->writeToExodus("what.exo");

    panzer_stk::PlaneMatcher top_matcher(0,2);
    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",top_matcher,"edge"));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher,"edge"));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HCURL_HEX_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    if(myRank==0) {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(1));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(5));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(3));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(7));

       // top/bottom matching
       TEST_EQUALITY(conn0[0],conn2[2]);
       TEST_EQUALITY(conn0[9],conn2[10]);
       TEST_EQUALITY(conn0[4],conn2[6]);
       TEST_EQUALITY(conn0[8],conn2[11]);
       TEST_EQUALITY(conn1[0],conn3[2]);
       TEST_EQUALITY(conn1[9],conn3[10]);
       TEST_EQUALITY(conn1[4],conn3[6]);
       TEST_EQUALITY(conn1[8],conn3[11]);

       // back/front matching
       TEST_EQUALITY(conn0[0],conn1[4]);
       TEST_EQUALITY(conn0[1],conn1[5]);
       TEST_EQUALITY(conn0[2],conn1[6]);
       TEST_EQUALITY(conn0[3],conn1[7]);
       TEST_EQUALITY(conn2[0],conn3[4]);
       TEST_EQUALITY(conn2[1],conn3[5]);
       TEST_EQUALITY(conn2[2],conn3[6]);
       TEST_EQUALITY(conn2[3],conn3[7]);
    }
    else {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(2));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(6));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(4));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(8));

       // top/bottom matching
       TEST_EQUALITY(conn0[0],conn2[2]);
       TEST_EQUALITY(conn0[9],conn2[10]);
       TEST_EQUALITY(conn0[4],conn2[6]);
       TEST_EQUALITY(conn0[8],conn2[11]);
       TEST_EQUALITY(conn1[0],conn3[2]);
       TEST_EQUALITY(conn1[9],conn3[10]);
       TEST_EQUALITY(conn1[4],conn3[6]);
       TEST_EQUALITY(conn1[8],conn3[11]);

       // back/front matching
       TEST_EQUALITY(conn0[0],conn1[4]);
       TEST_EQUALITY(conn0[1],conn1[5]);
       TEST_EQUALITY(conn0[2],conn1[6]);
       TEST_EQUALITY(conn0[3],conn1[7]);
       TEST_EQUALITY(conn2[0],conn3[4]);
       TEST_EQUALITY(conn2[1],conn3[5]);
       TEST_EQUALITY(conn2[2],conn3[6]);
       TEST_EQUALITY(conn2[3],conn3[7]);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_hex_face)
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
    // mesh->writeToExodus("what.exo");

    panzer_stk::PlaneMatcher top_matcher(0,2);
    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",top_matcher,"face"));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher,"face"));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    if(myRank==0) {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(1));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(5));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(3));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(7));

       // top/bottom matching
       TEST_EQUALITY(conn0[0],conn2[2]);
       TEST_EQUALITY(conn0[2],conn2[0]);
       TEST_EQUALITY(conn1[0],conn3[2]);
       TEST_EQUALITY(conn1[2],conn3[0]);

       // back/front matching
       TEST_EQUALITY(conn0[4],conn1[5]);
       TEST_EQUALITY(conn0[5],conn1[4]);
       TEST_EQUALITY(conn2[4],conn3[5]);
       TEST_EQUALITY(conn2[5],conn3[4]);
    }
    else {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(2));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(6));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(4));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(8));

       // top/bottom matching
       TEST_EQUALITY(conn0[0],conn2[2]);
       TEST_EQUALITY(conn0[2],conn2[0]);
       TEST_EQUALITY(conn1[0],conn3[2]);
       TEST_EQUALITY(conn1[2],conn3[0]);

       // back/front matching
       TEST_EQUALITY(conn0[4],conn1[5]);
       TEST_EQUALITY(conn0[5],conn1[4]);
       TEST_EQUALITY(conn2[4],conn3[5]);
       TEST_EQUALITY(conn2[5],conn3[4]);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_hex_quarter)
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

    // mesh->writeToExodus("what.exo");

    panzer_stk::QuarterPlaneMatcher quarter_matcher(0,2,1); // (xz)y
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("back","left",quarter_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);

    if(myRank==0) {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(1));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(5));
       const auto * conn2 = connMngr->getConnectivity(mesh->elementLocalId(3));
       const auto * conn3 = connMngr->getConnectivity(mesh->elementLocalId(7));

       TEST_EQUALITY(conn0[0],0);
       TEST_EQUALITY(conn0[1],9);
       TEST_EQUALITY(conn0[2],12);
       TEST_EQUALITY(conn0[3],3);
       TEST_EQUALITY(conn0[4],9);
       TEST_EQUALITY(conn0[5],10);
       TEST_EQUALITY(conn0[6],13);
       TEST_EQUALITY(conn0[7],12);

       TEST_EQUALITY(conn1[0],9);
       TEST_EQUALITY(conn1[1],10);
       TEST_EQUALITY(conn1[2],13);
       TEST_EQUALITY(conn1[3],12);
       TEST_EQUALITY(conn1[4],18);
       TEST_EQUALITY(conn1[5],19);
       TEST_EQUALITY(conn1[6],22);
       TEST_EQUALITY(conn1[7],21);

       TEST_EQUALITY(conn2[0],3);
       TEST_EQUALITY(conn2[1],12);
       TEST_EQUALITY(conn2[2],15);
       TEST_EQUALITY(conn2[3],6);
       TEST_EQUALITY(conn2[4],12);
       TEST_EQUALITY(conn2[5],13);
       TEST_EQUALITY(conn2[6],16);
       TEST_EQUALITY(conn2[7],15);

       TEST_EQUALITY(conn3[0],12);
       TEST_EQUALITY(conn3[1],13);
       TEST_EQUALITY(conn3[2],16);
       TEST_EQUALITY(conn3[3],15);
       TEST_EQUALITY(conn3[4],21);
       TEST_EQUALITY(conn3[5],22);
       TEST_EQUALITY(conn3[6],25);
       TEST_EQUALITY(conn3[7],24);
    }
    else {
       const auto * conn0 = connMngr->getConnectivity(mesh->elementLocalId(2));
       const auto * conn1 = connMngr->getConnectivity(mesh->elementLocalId(4));

       TEST_EQUALITY(conn0[0],9);
       TEST_EQUALITY(conn0[1],18);
       TEST_EQUALITY(conn0[2],21);
       TEST_EQUALITY(conn0[3],12);
       TEST_EQUALITY(conn0[4],10);
       TEST_EQUALITY(conn0[5],11);
       TEST_EQUALITY(conn0[6],14);
       TEST_EQUALITY(conn0[7],13);
       // passed

       TEST_EQUALITY(conn1[0],12);
       TEST_EQUALITY(conn1[1],21);
       TEST_EQUALITY(conn1[2],24);
       TEST_EQUALITY(conn1[3],15);
       TEST_EQUALITY(conn1[4],13);
       TEST_EQUALITY(conn1[5],14);
       TEST_EQUALITY(conn1[6],17);
       TEST_EQUALITY(conn1[7],16);
    }
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_3d_yz_xy_periodic)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    //int myRank = Comm.MyPID();

    // panzer::pauseToAttach();

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",1);
       pl->set("Y Blocks",2);
       pl->set("Z Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",2);
       pl->set("Y Procs",1);
       pl->set("Z Procs",1);
       pl->set("X0",0.0);
       pl->set("Y0",0.0);
       pl->set("Z0",0.0);
       pl->set("Xf",4.0);
       pl->set("Yf",2.0);
       pl->set("Zf",0.25);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",side_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_3d_yz_xy_periodic_edge)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    //int myRank = Comm.MyPID();

    // panzer::pauseToAttach();

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",1);
       pl->set("Y Blocks",2);
       pl->set("Z Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",2);
       pl->set("Y Procs",1);
       pl->set("Z Procs",1);
       pl->set("X0",0.0);
       pl->set("Y0",0.0);
       pl->set("Z0",0.0);
       pl->set("Xf",4.0);
       pl->set("Yf",2.0);
       pl->set("Zf",0.25);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",side_matcher,"edge"));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher,"edge"));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HCURL_HEX_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, conn_manager_3d_yz_xy_periodic_face)
  {
    using Teuchos::RCP;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);
    //int myRank = Comm.MyPID();

    // panzer::pauseToAttach();

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",1);
       pl->set("Y Blocks",2);
       pl->set("Z Blocks",1);
       pl->set("X Elements",4);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",2);
       pl->set("Y Procs",1);
       pl->set("Z Procs",1);
       pl->set("X0",0.0);
       pl->set("Y0",0.0);
       pl->set("Z0",0.0);
       pl->set("Xf",4.0);
       pl->set("Yf",2.0);
       pl->set("Zf",0.25);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    panzer_stk::PlaneMatcher side_matcher(1,2);
    panzer_stk::PlaneMatcher front_matcher(0,1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",side_matcher,"face"));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("front","back",front_matcher,"face"));

    // connection manager
    /////////////////////////////////////////////
    Teuchos::RCP<panzer::ConnManager> connMngr
          = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    RCP<const panzer::FieldPattern> fp
         = buildFieldPattern<Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::exec_space,double,double> >();
    connMngr->buildConnectivity(*fp);
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

  TEUCHOS_UNIT_TEST(periodic_mesh, getSideIdsAndCoords_edge)
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
         RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top","edge");

    TEST_EQUALITY(idsAndCoords.first->size(),4);
    TEST_EQUALITY(idsAndCoords.second->size(),4);
  }

  TEUCHOS_UNIT_TEST(periodic_mesh, getSideIdsAndCoords_face)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    TEUCHOS_ASSERT(Comm.NumProc()==2);

    panzer_stk::CubeHexMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("Z Blocks",1);
       pl->set("X Elements",2);
       pl->set("Y Elements",1);
       pl->set("Z Elements",1);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // mesh->writeToExodus("output.exo");

    panzer_stk::CoordMatcher x_matcher(0);
    std::pair<RCP<std::vector<std::size_t> >,
         RCP<std::vector<Tuple<double,3> > > > idsAndCoords = panzer_stk::periodic_helpers::getSideIdsAndCoords(*mesh,"top","face");

    TEST_EQUALITY(idsAndCoords.first->size(),4);
    TEST_EQUALITY(idsAndCoords.second->size(),4);
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

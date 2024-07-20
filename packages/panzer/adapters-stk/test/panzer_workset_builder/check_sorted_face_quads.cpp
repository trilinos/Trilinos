// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//

#include <iomanip>

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_IntegrationDescriptor.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_SubcellConnectivity.hpp"


#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(check_sorted_face_quads, basic)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0";
    std::string sideset = "left";

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set<int>   ("X Blocks",     1);
    pl->set<int>   ("X Elements",   2);
    pl->set<double>("X0",         0.0);
    pl->set<double>("Xf",         1.0);
    pl->set<int>   ("X Procs",      1);
    pl->set<int>   ("Y Blocks",     1);
    pl->set<int>   ("Y Elements",   2);
    pl->set<double>("Y0",           0);
    pl->set<double>("Yf",         1.0);
    pl->set<int>   ("Y Procs",      1);
    pl->sublist("Periodic BCs");
    pl->sublist("Periodic BCs").set<std::string>("Periodic Condition 1","y-all 1e-8: left;right");
    pl->sublist("Periodic BCs").set<std::string>("Periodic Condition 2","x-all 1e-8: top;bottom");
    pl->sublist("Periodic BCs").set<int>("Count",2);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

    // build DOF Manager (with a single HDiv basis)
    /////////////////////////////////////////////////////////////

    // build the connection manager
    const RCP<panzer::ConnManager>
      conn_manager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::DOFManager> dof_manager
        = rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // build an intrepid basis and a related field pattern for seeding the DOFManager
    {
       RCP<IntrepidBasis> intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",1,
                                                                                      *mesh->getCellTopology(element_block));
      RCP<panzer::Intrepid2FieldPattern> field_pattern = rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      dof_manager->addField(element_block, "T", field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    // build WorksetContainer
    //////////////////////////////////////////////////////////////

    panzer::IntegrationDescriptor sid(3, panzer::IntegrationDescriptor::SURFACE);
    panzer::BasisDescriptor bd(2, "HGrad");
    std::map<std::string, panzer::WorksetNeeds> wkstRequirements;
    wkstRequirements[element_block].addIntegrator(sid);
    wkstRequirements[element_block].addBasis(bd);

    RCP<panzer_stk::WorksetFactory> wkstFactory
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer(wkstFactory,wkstRequirements));

    wkstContainer->setGlobalIndexer(dof_manager);

    panzer::WorksetDescriptor workset_descriptor(element_block, panzer::WorksetSizeType::ALL_ELEMENTS, true,false);

    auto worksets = wkstContainer->getWorksets(workset_descriptor);

    TEST_ASSERT(worksets->size()==1);

    auto & workset = (*worksets)[0];
    auto rot_matrices = workset.getIntegrationValues(sid).surface_rotation_matrices;
    auto normals = workset.getIntegrationValues(sid).surface_normals;
    auto ip_coordinates = workset.getIntegrationValues(sid).ip_coordinates;

    const panzer::SubcellConnectivity & face_connectivity = workset.getFaceConnectivity();


    const int num_points = normals.extent(1);
    const int num_faces = face_connectivity.numSubcells();
    const int num_faces_per_cell = face_connectivity.numSubcellsOnCellHost(0);
    const int num_points_per_face = num_points / num_faces_per_cell;

    TEST_EQUALITY(num_faces,8); // 3*2-2 vertical edges (2 removed by periodicity)
                                // 3*2-2 horizontal edges (2 removed by periodicty)
    TEST_EQUALITY(num_faces_per_cell,4);
    TEST_EQUALITY(num_points,2*4);

    auto normals_view = PHX::as_view(normals);
    auto normals_h = Kokkos::create_mirror_view(normals_view);
    Kokkos::deep_copy(normals_h, normals_view);

    auto ip_coordinates_view = PHX::as_view(ip_coordinates);
    auto ip_coordinates_h = Kokkos::create_mirror_view(ip_coordinates_view);
    Kokkos::deep_copy(ip_coordinates_h, ip_coordinates_view);

    auto rot_matrices_view = PHX::as_view(rot_matrices);
    auto rot_matrices_h = Kokkos::create_mirror_view(rot_matrices_view);
    Kokkos::deep_copy(rot_matrices_h, rot_matrices_view);

    for(int face=0; face<num_faces; ++face) {
      int cell_l = face_connectivity.cellForSubcellHost(face,0);
      int cell_r = face_connectivity.cellForSubcellHost(face,1);

      int local_face_l = face_connectivity.localSubcellForSubcellHost(face,0);
      int local_face_r = face_connectivity.localSubcellForSubcellHost(face,1);

      out << "Face " << face << ": "
          << "Left "  << cell_l << " ( " << local_face_l << " ), "
          << "Right " << cell_r << " ( " << local_face_r << " )" << std::endl;
    }

    out << "Num Dim = " << normals.extent(2) << std::endl;

    // test a periodic horizontal pair
    for(int point=0;point<num_points_per_face;point++) {
      const int cell_l = 0;
      const int cell_r = 1;
      const int local_face_l = 3; // Left side of left cell
      const int local_face_r = 1; // Right side of right cell
      const int point_l = local_face_l * num_points_per_face + point;
      const int point_r = local_face_r * num_points_per_face + point;

      out << std::endl;
      out << std::scientific;
      out << std::setprecision(16) << "LEFT Normal = [ " << normals_h(cell_l,point_l,0) << ", "
                                                         << normals_h(cell_l,point_l,1) << " ]" << std::endl;
      out << std::setprecision(16) << "RGHT Normal = [ " << normals_h(cell_r,point_r,0) << ", "
                                                         << normals_h(cell_r,point_r,1) << " ]" << std::endl;

      out << std::endl;
      out << std::setprecision(16) << "LEFT Coords = [ " << ip_coordinates_h(cell_l,point_l,0) << ", "
                                                         << ip_coordinates_h(cell_l,point_l,1) << " ]" << std::endl;
      out << std::setprecision(16) << "RGHT Coords = [ " << ip_coordinates_h(cell_r,point_r,0) << ", "
                                                         << ip_coordinates_h(cell_r,point_r,1) << " ]" << std::endl;

      TEST_FLOATING_EQUALITY(ip_coordinates_h(cell_l,point_l,0), 0.0, 1.0e-14);
      TEST_FLOATING_EQUALITY(ip_coordinates_h(cell_r,point_r,0), 1.0, 1.0e-14);
      TEST_FLOATING_EQUALITY(ip_coordinates_h(cell_l,point_l,1), ip_coordinates_h(cell_r,point_r,1), 1.0e-14);

      out << "LEFT rotation" << std::endl;
      out << std::setprecision(16)
          << "  " << rot_matrices_h(cell_l,point_l,0,0) << ", " << rot_matrices_h(cell_l,point_l,0,1) << "  " << rot_matrices_h(cell_l,point_l,0,2) << std::endl
          << "  " << rot_matrices_h(cell_l,point_l,1,0) << ", " << rot_matrices_h(cell_l,point_l,1,1) << "  " << rot_matrices_h(cell_l,point_l,1,2) << std::endl
          << "  " << rot_matrices_h(cell_l,point_l,2,0) << ", " << rot_matrices_h(cell_l,point_l,2,1) << "  " << rot_matrices_h(cell_l,point_l,2,2) << std::endl;
    }

    // test an interior horizontal pair
    for(int point=0;point<num_points_per_face;point++) {
      const int cell_l = 0;
      const int cell_r = 1;
      const int local_face_l = 1; // Right side of left cell
      const int local_face_r = 3; // Left side of right cell
      const int point_l = local_face_l * num_points_per_face + point;
      const int point_r = local_face_r * num_points_per_face + point;

      out << std::endl;
      out << std::scientific;
      out << std::setprecision(16) << "LEFT Normal = [ " << normals_h(cell_l,point_l,0) << ", "
                                                         << normals_h(cell_l,point_l,1) << " ]" << std::endl;
      out << std::setprecision(16) << "RGHT Normal = [ " << normals_h(cell_r,point_r,0) << ", "
                                                         << normals_h(cell_r,point_r,1) << " ]" << std::endl;

      out << std::endl;
      out << std::setprecision(16) << "LEFT Coords = [ " << ip_coordinates_h(cell_l,point_l,0) << ", "
                                                         << ip_coordinates_h(cell_l,point_l,1) << " ]" << std::endl;
      out << std::setprecision(16) << "RGHT Coords = [ " << ip_coordinates_h(cell_r,point_r,0) << ", "
                                                         << ip_coordinates_h(cell_r,point_r,1) << " ]" << std::endl;

      TEST_FLOATING_EQUALITY(ip_coordinates_h(cell_l,point_l,0), ip_coordinates_h(cell_r,point_r,0), 1.0e-14);
      TEST_FLOATING_EQUALITY(ip_coordinates_h(cell_l,point_l,1), ip_coordinates_h(cell_r,point_r,1), 1.0e-14);

      out << "LEFT rotation" << std::endl;
      out << std::setprecision(16)
          << "  " << rot_matrices_h(cell_l,point_l,0,0) << ", " << rot_matrices_h(cell_l,point_l,0,1) << "  " << rot_matrices_h(cell_l,point_l,0,2) << std::endl
          << "  " << rot_matrices_h(cell_l,point_l,1,0) << ", " << rot_matrices_h(cell_l,point_l,1,1) << "  " << rot_matrices_h(cell_l,point_l,1,2) << std::endl
          << "  " << rot_matrices_h(cell_l,point_l,2,0) << ", " << rot_matrices_h(cell_l,point_l,2,1) << "  " << rot_matrices_h(cell_l,point_l,2,2) << std::endl;
    }
  }
} // end namespace panzer

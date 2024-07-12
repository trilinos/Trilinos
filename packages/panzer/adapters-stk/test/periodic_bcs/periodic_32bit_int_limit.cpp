// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Kokkos_Core.hpp"

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"

#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

TEUCHOS_UNIT_TEST(periodic_bcs, 32_bit_int_limit)
{
  // This test requires 64 bit integers for the GIDs of Tpetra
  // objects. If the GIDs are not 64 bit then disable the test.
  if ( not( (std::numeric_limits<panzer::GlobalOrdinal>::max() >= std::numeric_limits<long long int>::max()) || 
            (std::numeric_limits<panzer::GlobalOrdinal>::max() >= std::numeric_limits<unsigned long long int>::max()) ) )
    {
      std::cout << "\nWARNING: Panzer_AdaptersSTK::periodic_32bit_int_limit test is disabled since Tpetra"
                << " was not configured with a 64 bit GID type.\n\n";
    return;
  }

  // To run in parallel we need to set an ioss property to split
  // exodus file across mpi ranks.
  setenv("IOSS_PROPERTIES", "DECOMPOSITION_METHOD=rib", 1);

  Teuchos::RCP<Teuchos::MpiComm<int>> Comm = Teuchos::rcp( new Teuchos::MpiComm<int>(MPI_COMM_WORLD) );

  using topo_RCP = Teuchos::RCP<const shards::CellTopology>;
  using basis_RCP = Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double>>;

  {
    // ==========================================================
    // Create a mesh from the file defined by the user
    // ==========================================================
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory = Teuchos::rcp(new panzer_stk::STK_ExodusReaderFactory());

    std::string input_file_name = "periodic_32bit_int_limit.exo";
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("File Name",input_file_name);

    // Optional: include periodic BCs
    if (true) {
      RCP<Teuchos::ParameterList> pbc = rcp(new Teuchos::ParameterList);
      pbc->set("Count",2);
      pbc->set("Periodic Condition 1","xz-all 1e-12,3D: top;bottom");
      pbc->set("Periodic Condition 2","yz-all 1e-12,3D: left;right");
#ifdef PANZER_HAVE_STKSEARCH
      pbc->set("Use BBox Search",true);
#endif
      pl->sublist("Periodic BCs").setParameters( *pbc );
    }

    mesh_factory->setParameterList(pl);
    Teuchos::RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(*(Comm->getRawMpiComm()));

    mesh_factory->completeMeshConstruction(*mesh,*(Comm->getRawMpiComm()));

    bool is_print_process = false;
    if (Comm->getRank() == 0)
      is_print_process = true;

    if (is_print_process)
      mesh->printMetaData(std::cout);

    vector<string> blocknames;
    mesh->getElementBlockNames(blocknames);

    // ==========================================================
    // Define the connectivity manager
    // ==========================================================

    Teuchos::RCP<panzer::ConnManager> conn = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // ==========================================================
    // Test out a few DOF managers
    // ==========================================================

    {
      if (is_print_process) {
        out << "================================================" << std::endl;
        out << " *** Testing one HGRAD variable *** " << std::endl;
      }
      Teuchos::RCP<panzer::DOFManager> DOF = Teuchos::rcp(new panzer::DOFManager());
      DOF->setConnManager(conn,*(Comm->getRawMpiComm()));
      DOF->setOrientationsRequired(true);

      for (size_t b=0; b<blocknames.size(); b++) {
        topo_RCP cellTopo = mesh->getCellTopology(blocknames[b]);
        basis_RCP basis = Teuchos::rcp(new Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::Device::execution_space,double,double>() );
        Teuchos::RCP<const panzer::Intrepid2FieldPattern> Pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
        DOF->addField(blocknames[b], "T", Pattern, panzer::FieldType::CG);
      }

      DOF->buildGlobalUnknowns();
      if (is_print_process) {
        DOF->printFieldInformation(out);
        out << "================================================" << std::endl << std::endl;
      }
    }

    {
      if (is_print_process) {
        out << "================================================" << std::endl;
        out << " *** Testing one HCURL variable *** " << std::endl;
      }
      Teuchos::RCP<panzer::DOFManager> DOF = Teuchos::rcp(new panzer::DOFManager());
      DOF->setConnManager(conn,*(Comm->getRawMpiComm()));
      DOF->setOrientationsRequired(true);

      for (size_t b=0; b<blocknames.size(); b++) {
        topo_RCP cellTopo = mesh->getCellTopology(blocknames[b]);
        basis_RCP basis = Teuchos::rcp(new Intrepid2::Basis_HCURL_HEX_I1_FEM<PHX::Device::execution_space,double,double>() );
        Teuchos::RCP<const panzer::Intrepid2FieldPattern> Pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
        DOF->addField(blocknames[b], "E", Pattern, panzer::FieldType::CG);
      }

      DOF->buildGlobalUnknowns();
      if (is_print_process) {
        DOF->printFieldInformation(out);
        out << "================================================" << std::endl << std::endl;
      }
    }

    {
      if (is_print_process) {
        out << "================================================" << std::endl;
        out << " *** Testing one HDIV variable *** " << std::endl;
      }
      Teuchos::RCP<panzer::DOFManager> DOF = Teuchos::rcp(new panzer::DOFManager());
      DOF->setConnManager(conn,*(Comm->getRawMpiComm()));
      DOF->setOrientationsRequired(true);

      for (size_t b=0; b<blocknames.size(); b++) {
        topo_RCP cellTopo = mesh->getCellTopology(blocknames[b]);
        basis_RCP basis = Teuchos::rcp(new Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::Device::execution_space,double,double>() );
        Teuchos::RCP<const panzer::Intrepid2FieldPattern> Pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
        DOF->addField(blocknames[b], "B", Pattern, panzer::FieldType::CG);
      }

      DOF->buildGlobalUnknowns();
      if (is_print_process) {
        DOF->printFieldInformation(out);
        out << "================================================" << std::endl << std::endl;
      }
    }

    {
      if (is_print_process) {
        out << "================================================" << std::endl;
        out << " *** Testing one HVOL variable *** " << std::endl;
      }
      Teuchos::RCP<panzer::DOFManager> DOF = Teuchos::rcp(new panzer::DOFManager());
      DOF->setConnManager(conn,*(Comm->getRawMpiComm()));
      DOF->setOrientationsRequired(true);

      for (size_t b=0; b<blocknames.size(); b++) {
        topo_RCP cellTopo = mesh->getCellTopology(blocknames[b]);
        basis_RCP basis = Teuchos::rcp(new Intrepid2::Basis_HVOL_C0_FEM<PHX::Device::execution_space,double,double>(*cellTopo) );
        Teuchos::RCP<const panzer::Intrepid2FieldPattern> Pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis));
        DOF->addField(blocknames[b], "p", Pattern, panzer::FieldType::CG);
      }

      DOF->buildGlobalUnknowns();
      if (is_print_process) {
        DOF->printFieldInformation(out);
        out << "================================================" << std::endl << std::endl;
      }
    }
  }
}

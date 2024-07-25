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

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"

#include "Panzer_WorksetContainer.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_DOFManager.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(hdiv_basis, check_dirac)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0_0";
    int basis_order = 4;
    int workset_size = 10;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Elements",2);
    pl->set("Y Elements",2);
    pl->set("Z Elements",2);

    panzer_stk::CubeHexMeshFactory factory;
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
    RCP<IntrepidBasis> hdiv_intrepid_basis, hcurl_intrepid_basis;
    {
       hdiv_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HDiv",basis_order,
                                                                                      *mesh->getCellTopology(element_block));
       hcurl_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HCurl",basis_order,
                                                                                      *mesh->getCellTopology(element_block));
      RCP<panzer::Intrepid2FieldPattern> hdiv_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hdiv_intrepid_basis));
      RCP<panzer::Intrepid2FieldPattern> hcurl_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hcurl_intrepid_basis));

      dof_manager->addField(element_block, "B", hdiv_field_pattern);
      dof_manager->addField(element_block, "E", hcurl_field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    // build WorksetContainer
    //////////////////////////////////////////////////////////////
    BasisDescriptor hdiv_basis_desc(basis_order,"HDiv");
    BasisDescriptor hcurl_basis_desc(basis_order,"HCurl");

    out << "workset container setup [start]" << std::endl;

    RCP<panzer_stk::WorksetFactory> wkstFactory
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    {
      WorksetNeeds needs;
      needs.addBasis(hdiv_basis_desc);
      needs.addBasis(hcurl_basis_desc);
      needs.addPoint(hdiv_basis_desc.getPointDescriptor());
      wkstContainer->setNeeds(element_block,needs);
    }
    wkstContainer->setGlobalIndexer(dof_manager);
    wkstContainer->setWorksetSize(workset_size);

    out << "workset container setup [complete]" << std::endl;

    // Get worksets
    ///////////////////////////////////////////////////////////////////////
    out << "getting worksets [start]" << std::endl;

    //  this must use this descriptor!
    panzer::WorksetDescriptor workset_descriptor(element_block, panzer::WorksetSizeType::ALL_ELEMENTS, true,true);
    std::vector<Workset> & worksets = *wkstContainer->getWorksets(workset_descriptor);

    out << "getting worksets [complete]" << std::endl;

    // get BasisValues2
    ///////////////////////////////////////////////////////////////////////

    const PointValues2<double> & point_values = worksets[0].getPointValues(hdiv_basis_desc.getPointDescriptor());
    const BasisValues2<double> & hdiv_basis_values = worksets[0].getBasisValues(hdiv_basis_desc,hdiv_basis_desc.getPointDescriptor());
    const BasisValues2<double> & hcurl_basis_values = worksets[0].getBasisValues(hcurl_basis_desc,hdiv_basis_desc.getPointDescriptor());

    auto hdiv_basis_vector = hdiv_basis_values.getVectorBasisValues(false);
    auto hcurl_curl_basis = hcurl_basis_values.getCurlVectorBasis(false);

    // check some sizing stuff
    ///////////////////////////////////////////////////////////////////////

    TEST_EQUALITY(Teuchos::as<int>(hdiv_basis_vector.extent(1)),hdiv_intrepid_basis->getCardinality());
    TEST_EQUALITY(hdiv_basis_vector.extent(1),hdiv_basis_vector.extent(2));

    TEST_EQUALITY(Teuchos::as<int>(hcurl_curl_basis.extent(1)),hcurl_intrepid_basis->getCardinality());
    TEST_EQUALITY(hcurl_curl_basis.extent(2),hdiv_basis_vector.extent(2));

    TEST_EQUALITY(Teuchos::as<int>(point_values.coords_ref.extent(0)),hdiv_intrepid_basis->getCardinality());
    TEST_EQUALITY(point_values.coords_ref.extent(1),3u);

    TEST_EQUALITY(Teuchos::as<int>(point_values.point_coords.extent(1)),hdiv_intrepid_basis->getCardinality());
    TEST_EQUALITY(point_values.point_coords.extent(2),3u);

    // print out basis values
    ///////////////////////////////////////////////////////////////////////

    auto hdiv_basis_vector_view = PHX::as_view(hdiv_basis_vector);
    auto hdiv_basis_vector_h = Kokkos::create_mirror_view(hdiv_basis_vector_view);
    Kokkos::deep_copy(hdiv_basis_vector_h, hdiv_basis_vector_view);

    // print out phi_i(x_j).phi_j(x_j)
    for(size_t c=0;c<hdiv_basis_vector.extent(0);c++) {
      out << "cell " << c << std::endl;
      for(size_t i=0;i<hdiv_basis_vector.extent(1);i++) {
        // compute diagonal magnitude
        double diagonal = 0.0;
        for(size_t d=0;d<hdiv_basis_vector.extent(3);d++)
          diagonal += hdiv_basis_vector_h(c,i,i,d)* hdiv_basis_vector_h(c,i,i,d);

        out << "   ";
        for(size_t j=0;j<hdiv_basis_vector_h.extent(2);j++) {
          double entry = 0.0;

          // loop over dimension
          for(size_t d=0;d<hdiv_basis_vector.extent(3);d++)
            entry += hdiv_basis_vector_h(c,i,j,d)* hdiv_basis_vector_h(c,j,j,d);

          out << std::fixed << std::setprecision(2) << std::setw(8) << entry / diagonal;

          if(i==j) {
            TEST_ASSERT(std::fabs(entry/diagonal)-1.0 <= 1.0e-15);
          }
          else {
            TEST_ASSERT(std::fabs(entry/diagonal) <= 1.0e-15);
          }
        } // end j
        out << std::endl;
      } // end i
    } // end c
  }

  TEUCHOS_UNIT_TEST(hdiv_basis, check_weighted_orientations)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Intrepid2::Basis<PHX::Device::execution_space,double,double> IntrepidBasis;

    std::string element_block = "eblock-0_0_0";
    int basis_order = 1;
    int workset_size = 2;

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Elements",2);
    pl->set("Y Elements",1);
    pl->set("Z Elements",1);

    panzer_stk::CubeHexMeshFactory factory;
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
    RCP<IntrepidBasis> hdiv_intrepid_basis, hcurl_intrepid_basis;
    {
       hdiv_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HDiv",basis_order,
                                                                                      *mesh->getCellTopology(element_block));
       hcurl_intrepid_basis
           = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>("HCurl",basis_order,
                                                                                      *mesh->getCellTopology(element_block));
      RCP<panzer::Intrepid2FieldPattern> hdiv_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hdiv_intrepid_basis));
      RCP<panzer::Intrepid2FieldPattern> hcurl_field_pattern = rcp(new panzer::Intrepid2FieldPattern(hcurl_intrepid_basis));

      dof_manager->addField(element_block, "B", hdiv_field_pattern);
      dof_manager->addField(element_block, "E", hcurl_field_pattern);
    }
    dof_manager->buildGlobalUnknowns();

    // build WorksetContainer
    //////////////////////////////////////////////////////////////
    BasisDescriptor hdiv_basis_desc(basis_order,"HDiv");
    BasisDescriptor hcurl_basis_desc(basis_order,"HCurl");
    IntegrationDescriptor quad_desc(2,IntegrationDescriptor::VOLUME);

    out << "workset container setup [start]" << std::endl;

    RCP<panzer_stk::WorksetFactory> wkstFactory
       = rcp(new panzer_stk::WorksetFactory(mesh)); // build STK workset factory
    RCP<panzer::WorksetContainer> wkstContainer     // attach it to a workset container (uses lazy evaluation)
       = rcp(new panzer::WorksetContainer);
    wkstContainer->setFactory(wkstFactory);
    {
      WorksetNeeds needs;
      needs.bases.push_back(Teuchos::rcp(new panzer::PureBasis(hdiv_basis_desc,mesh->getCellTopology(element_block),workset_size)));
      needs.bases.push_back(Teuchos::rcp(new panzer::PureBasis(hcurl_basis_desc,mesh->getCellTopology(element_block),workset_size)));
      needs.rep_field_name.push_back("B");
      needs.rep_field_name.push_back("E");
      needs.int_rules.push_back(Teuchos::rcp(new panzer::IntegrationRule(quad_desc,mesh->getCellTopology(element_block),workset_size)));
      needs.cellData = CellData(workset_size,mesh->getCellTopology(element_block));

      wkstContainer->setNeeds(element_block,needs);
    }
    wkstContainer->setGlobalIndexer(dof_manager);
    wkstContainer->setWorksetSize(workset_size);

    out << "workset container setup [complete]" << std::endl;

    // Get worksets
    ///////////////////////////////////////////////////////////////////////
    out << "getting worksets [start]" << std::endl;

    //  this must use this descriptor!
    // panzer::WorksetDescriptor workset_descriptor(element_block, panzer::WorksetSizeType::ALL_ELEMENTS, true,true);
    panzer::WorksetDescriptor workset_descriptor(element_block);
    std::vector<Workset> & worksets = *wkstContainer->getWorksets(workset_descriptor);

    out << "getting worksets [complete]" << std::endl;

    // get BasisValues2
    ///////////////////////////////////////////////////////////////////////
    out <<  worksets[0].bases.size() << std::endl;

    const BasisValues2<double> & hdiv_basis_values = *worksets[0].bases[0];
    out << (*worksets[0].basis_names)[0]
        << " " << (*worksets[0].basis_names)[1] << std::endl;

    auto hdiv_weighted_basis_vector = hdiv_basis_values.weighted_basis_vector;
    auto hdiv_weighted_basis_vector_view = PHX::as_view(hdiv_weighted_basis_vector);
    auto hdiv_weighted_basis_vector_h = Kokkos::create_mirror_view(hdiv_weighted_basis_vector_view);
    Kokkos::deep_copy(hdiv_weighted_basis_vector_h, hdiv_weighted_basis_vector_view);

    auto hdiv_basis_vector = hdiv_basis_values.basis_vector;
    auto hdiv_basis_vector_view = PHX::as_view(hdiv_basis_vector);
    auto hdiv_basis_vector_h = Kokkos::create_mirror_view(hdiv_basis_vector_view);
    Kokkos::deep_copy(hdiv_basis_vector_h, hdiv_basis_vector_view);

    // check some sizing stuff
    ///////////////////////////////////////////////////////////////////////

    TEST_EQUALITY(worksets[0].num_cells,2);
    TEST_EQUALITY(Teuchos::as<int>(hdiv_weighted_basis_vector.extent(0)),2);
    TEST_EQUALITY(Teuchos::as<int>(hdiv_weighted_basis_vector.extent(1)),6);
    TEST_EQUALITY(Teuchos::as<int>(hdiv_weighted_basis_vector.extent(2)),8);
    TEST_EQUALITY(Teuchos::as<int>(hdiv_weighted_basis_vector.extent(3)),3);

    // print out basis values
    ///////////////////////////////////////////////////////////////////////

    out << "BASIS VECTOR\n" << std::endl;
    for(size_t c=0;c<2;c++) {
      out << "cell " << c << " = ";
      std::vector<panzer::GlobalOrdinal> gids;
      dof_manager->getElementGIDs(c,gids);

      out << "  gids = ";
      for(auto g : gids)
        out << g << " ";
      out << std::endl;

      for(size_t b=0;b<hdiv_basis_vector.extent(1);b++) {
        out << "    basis " << b << " = ";
        for(size_t q=0;q<hdiv_basis_vector.extent(2);q++) {
          out << "[ ";
          for(size_t d=0;d<hdiv_basis_vector.extent(3);d++) {
            out << hdiv_basis_vector_h(c,b,q,d) << ", ";
          } // end d
          out << "], ";
        } // end q
        out << std::endl;
      } // end b
      out << std::endl;
    } // end c

    double sign_diff =   hdiv_basis_vector_h(0,1,0,0)/std::fabs(hdiv_basis_vector_h(0,1,0,0))
                       - hdiv_basis_vector_h(1,3,0,0)/std::fabs(hdiv_basis_vector_h(1,3,0,0));
    TEST_ASSERT(std::fabs(sign_diff) <= 1e-15);

    out << "WEIGHTED BASIS VECTOR\n" << std::endl;
    for(size_t c=0;c<2;c++) {
      out << "cell " << c << " = ";
      std::vector<panzer::GlobalOrdinal> gids;
      dof_manager->getElementGIDs(c,gids);

      out << "  gids = ";
      for(auto g : gids)
        out << g << " ";
      out << std::endl;

      for(size_t b=0;b<hdiv_weighted_basis_vector.extent(1);b++) {
        out << "    basis " << b << " = ";
        for(size_t q=0;q<hdiv_weighted_basis_vector.extent(2);q++) {
          out << "[ ";
          for(size_t d=0;d<hdiv_weighted_basis_vector.extent(3);d++) {
            out << hdiv_weighted_basis_vector_h(c,b,q,d) << ", ";
          } // end d
          out << "], ";
        } // end q
        out << std::endl;
      } // end b
      out << std::endl;
    } // end c
  }
}

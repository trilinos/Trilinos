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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCPStdSharedPtrConversions.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef PANZER_HAVE_PERCEPT
#include "percept/PerceptMesh.hpp"
#endif

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/CreateAdjacentEntities.hpp"
#include "stk_mesh/base/DumpMeshInfo.hpp"

namespace panzer_stk {

// Acceptance test that deomstrates panzer's use of pamgen/stk
TEUCHOS_UNIT_TEST(tPamgenFactory, acceptance)
{
  using namespace Teuchos;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // IMPORTANT: This test must be run on one processor. Multiple
  // processors split the surface counts.
  TEST_EQUALITY(Teuchos::DefaultComm<int>::getComm()->getSize(),1);

  const bool verbose = false;

  RCP<stk::mesh::MetaData> metaData;
  RCP<stk::mesh::BulkData> bulkData;
  const std::string output_exodus_file_name = "pamgen_acceptance_output.exo";
  std::vector<std::size_t> gold_ss_values{6,3,2,6,3,2,6,3,2,6,3,2};

  // Create a stk::mesh from a pamgen using io broker, then delete the broker
  {
    out << "\nCreating pamgen mesh." << std::endl;
    RCP<stk::io::StkMeshIoBroker> broker = rcp(new stk::io::StkMeshIoBroker(MPI_COMM_WORLD));
    broker->use_simple_fields();
    broker->add_mesh_database("pamgen_test.gen", "pamgen", stk::io::READ_MESH);
    broker->create_input_mesh();
    metaData = Teuchos::rcp(broker->meta_data_ptr());
    std::unique_ptr<stk::mesh::BulkData> bulkUPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create(Teuchos::get_shared_ptr(metaData));
    bulkData = Teuchos::rcp(bulkUPtr.release());
    broker->set_bulk_data(Teuchos::get_shared_ptr(bulkData));
    broker->add_all_mesh_fields_as_input_fields();
    broker->populate_bulk_data();

    if (verbose) {
      stk::mesh::impl::dump_all_meta_info(*metaData, out);
      stk::mesh::impl::dump_all_mesh_info(*bulkData, out);
    }
  }

  // Setup extra mesh information for edges and faces. This is the

  // **************************************************
  // * create_adjacent_entities wipes out sideset ids
  // being written to exodus! Sidesets do seem to still
  // be valid for the stk::mesh object though.
  // **************************************************
  if (true) {
    stk::mesh::PartVector emptyPartVector;
    stk::mesh::create_adjacent_entities(*bulkData,emptyPartVector);
  }

  // Verify the sidesets exist and are not empty in the mesh
  {
    out << "\nChecking sidesets before writing:" << std::endl;
    auto surfaces = metaData->get_surfaces_in_surface_to_block_map();
    std::size_t surface_index = 0;
    for (auto s : surfaces) {
      // out << s->name() << std::endl;
      stk::mesh::Selector side_selector = *s;
      // stk::mesh::Selector sel = side_selector & lo_selector;
      const stk::mesh::BucketVector& bv = bulkData->get_buckets(stk::topology::FACE_RANK,side_selector);
      std::size_t num_faces_in_sideset = 0;
      for (std::size_t b=0; b < bv.size(); ++b) {
        num_faces_in_sideset += bv[b]->size();
      }
      TEST_EQUALITY(gold_ss_values[surface_index],num_faces_in_sideset);
      ++surface_index;
    }
  }

  // Output mesh to exodus file
  {
    out << "\nWriting output file." << std::endl;
    RCP<stk::io::StkMeshIoBroker> broker = rcp(new stk::io::StkMeshIoBroker(MPI_COMM_WORLD));
    broker->set_bulk_data(Teuchos::get_shared_ptr(bulkData));
    auto meshIndex_ = broker->create_output_mesh(output_exodus_file_name,stk::io::PURPOSE_UNKNOWN);

    const stk::mesh::FieldVector& fields = metaData->get_fields();
    for (size_t i(0); i < fields.size(); ++i) {
      // Do NOT add MESH type stk fields to exodus io, but do add everything
      // else. This allows us to avoid having to catch a throw for
      // re-registering coordinates, sidesets, etc... Note that some
      // fields like LOAD_BAL don't always have a role assigned, so for
      // roles that point to null, we need to register them as well.
      auto role = stk::io::get_field_role(*fields[i]);
      if (role != nullptr) {
        if (*role != Ioss::Field::MESH)
          broker->add_field(meshIndex_, *fields[i]);
      } else {
        broker->add_field(meshIndex_, *fields[i]);
      }
    }

    double currentStateTime = 10.0;
    broker->begin_output_step(meshIndex_, currentStateTime);
    broker->write_defined_output_fields(meshIndex_);
    broker->end_output_step(meshIndex_);
  }

  // Delete mesh
  bulkData = Teuchos::null;
  metaData = Teuchos::null;

  // Read in mesh from exodus
  {
    out << "\nReading Exodus file." << std::endl;
    RCP<stk::io::StkMeshIoBroker> broker = rcp(new stk::io::StkMeshIoBroker(MPI_COMM_WORLD));
    broker->use_simple_fields();
    broker->add_mesh_database(output_exodus_file_name, "exodus", stk::io::READ_MESH);
    broker->create_input_mesh();
    metaData = Teuchos::rcp(broker->meta_data_ptr());
    std::unique_ptr<stk::mesh::BulkData> bulkUPtr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create(Teuchos::get_shared_ptr(metaData));
    bulkData = Teuchos::rcp(bulkUPtr.release());
    broker->set_bulk_data(Teuchos::get_shared_ptr(bulkData));
    broker->add_all_mesh_fields_as_input_fields();
    broker->populate_bulk_data();

    if (verbose) {
      stk::mesh::impl::dump_all_meta_info(*metaData, out);
      stk::mesh::impl::dump_all_mesh_info(*bulkData, out);
    }
  }

  // Verify the sidesets exist and are not empty in the mesh
  {
    out << "\nVerifying sidesets from exodus file:" << std::endl;
    auto surfaces = metaData->get_surfaces_in_surface_to_block_map();
    TEST_COMPARE(surfaces.size(),>,0);
    std::size_t surface_index = 0;
    for (auto s : surfaces) {
      stk::mesh::Selector side_selector = *s;
      const stk::mesh::BucketVector& bv = bulkData->get_buckets(stk::topology::FACE_RANK,side_selector);
      std::size_t num_faces_in_sideset = 0;
      for (std::size_t b=0; b < bv.size(); ++b) {
        num_faces_in_sideset += bv[b]->size();
      }
      TEST_EQUALITY(gold_ss_values[surface_index],num_faces_in_sideset);
      ++surface_index;
    }
  }

}

TEUCHOS_UNIT_TEST(tPamgenFactory, basic)
{
  using namespace Teuchos;
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string output_file_name = "pamgen_basic_test_outputcheck.gen";

  // IMPORTANT: This test must be run on one processor. Multiple
  // processors split the surface counts.
  TEST_EQUALITY(Teuchos::DefaultComm<int>::getComm()->getSize(),1);

  std::vector<std::size_t> gold_ss_values{6,3,2,6,3,2,6,3,2,6,3,2};

  // *******************************************
  // Create and write out a stk mesh from pamgen
  // *******************************************
  {
    out << "\nCreating pamgen mesh." << std::endl;
    RCP<ParameterList> p = parameterList();
    const std::string input_file_name = "pamgen_test.gen";
    p->set("File Name",input_file_name);
    p->set("File Type","Pamgen");

    STK_ExodusReaderFactory pamgenFactory;
    pamgenFactory.setParameterList(p);

    auto mesh = pamgenFactory.buildMesh(MPI_COMM_WORLD);
    auto metaData = mesh->getMetaData();
    auto bulkData = mesh->getBulkData();

    // Write sideset names
    std::vector<std::string> sideset_names;
    mesh->getSidesetNames(sideset_names);
    TEST_COMPARE(sideset_names.size(),>,0)
    out << "\nsideset names:\n";
    for (auto s : sideset_names)
      out << s << std::endl;
    out << "\n";

    // Get locally owned sidesets
    const stk::mesh::Part& lo_part = metaData->locally_owned_part();
    stk::mesh::Selector lo_selector = lo_part;

    auto surfaces = metaData->get_surfaces_in_surface_to_block_map();
    std::size_t surface_index = 0;
    for (auto s : surfaces) {
      stk::mesh::Selector side_selector = *s;
      stk::mesh::Selector sel = side_selector & lo_selector;
      const stk::mesh::BucketVector& bv = bulkData->get_buckets(stk::topology::FACE_RANK,sel);
      std::size_t num_faces_in_sideset = 0;
      for (std::size_t b=0; b < bv.size(); ++b) {
        num_faces_in_sideset += bv[b]->size();
      }
      TEST_EQUALITY(gold_ss_values[surface_index],num_faces_in_sideset);
      ++surface_index;
    }
    mesh->writeToExodus(output_file_name);

    if (false) {
      stk::mesh::impl::dump_all_meta_info(*metaData, out);
      stk::mesh::impl::dump_all_mesh_info(*bulkData, out);
    }
  }

  // *******************************************
  // Read back in pamgen mesh from exodus and check sideset exist
  // (i.e. repeat the above code using an exodus reader)
  // *******************************************
  {
    RCP<ParameterList> p = parameterList();
    p->set("File Name",output_file_name);
    p->set("File Type","Exodus");

    STK_ExodusReaderFactory factory;
    factory.setParameterList(p);

    auto mesh = factory.buildMesh(MPI_COMM_WORLD);
    auto metaData = mesh->getMetaData();
    auto bulkData = mesh->getBulkData();

    // Write sideset names
    std::vector<std::string> sideset_names;
    mesh->getSidesetNames(sideset_names);
    TEST_COMPARE(sideset_names.size(),>,0)
    out << "\nsideset names:\n";
    for (auto s : sideset_names)
      out << s << std::endl;
    out << "\n";

    // Get locally owned sidesets
    const stk::mesh::Part& lo_part = metaData->locally_owned_part();
    stk::mesh::Selector lo_selector = lo_part;

    auto surfaces = metaData->get_surfaces_in_surface_to_block_map();
    std::size_t surface_index = 0;
    for (auto s : surfaces) {
      stk::mesh::Selector side_selector = *s;
      stk::mesh::Selector sel = side_selector & lo_selector;
      const stk::mesh::BucketVector& bv = bulkData->get_buckets(stk::topology::FACE_RANK,sel);
      std::size_t num_faces_in_sideset = 0;
      for (std::size_t b=0; b < bv.size(); ++b) {
        num_faces_in_sideset += bv[b]->size();
      }
      TEST_EQUALITY(gold_ss_values[surface_index],num_faces_in_sideset);
      ++surface_index;
    }
    mesh->writeToExodus(output_file_name);
  }

}

TEUCHOS_UNIT_TEST(tPamgenFactory, getMeshDimension)
{
  TEST_EQUALITY(panzer_stk::getMeshDimension("pamgen_test.gen",MPI_COMM_WORLD,"Pamgen"),3);
}


#ifdef PANZER_HAVE_PERCEPT
// test the ability to save the percept mesh hierarchy
TEUCHOS_UNIT_TEST(tPamgenFactory, keepPerceptData)
{
  using namespace Teuchos;
  using Teuchos::RCP;
  using Teuchos::rcp;

  {
    out << "\nCreating pamgen mesh with parent elements." << std::endl;
    RCP<ParameterList> p = parameterList();
    const std::string input_file_name = "pamgen_test.gen";
    p->set("File Name",input_file_name);
    p->set("File Type","Pamgen");
    p->set("Levels of Uniform Refinement",1);
    p->set("Keep Percept Data",true);
    p->set("Keep Percept Parent Elements",true);

    RCP<STK_ExodusReaderFactory> pamgenFactory = rcp(new STK_ExodusReaderFactory());
    TEST_NOTHROW(pamgenFactory->setParameterList(p));

    RCP<STK_Interface> mesh = pamgenFactory->buildMesh(MPI_COMM_WORLD);
    RCP<percept::PerceptMesh> refinedMesh = mesh->getRefinedMesh();

    // check if the Percept data exists and has elements
    TEST_ASSERT(nonnull(refinedMesh));
    TEST_ASSERT(refinedMesh->get_number_elements()>0);

    // check if the parent information is stored
    std::vector<stk::mesh::EntityRank> ranks_to_be_deleted;
    ranks_to_be_deleted.push_back(stk::topology::ELEMENT_RANK);
    ranks_to_be_deleted.push_back(refinedMesh->side_rank());
    if (refinedMesh->get_spatial_dim() == 3)
      ranks_to_be_deleted.push_back(refinedMesh->edge_rank());

    //percept::SetOfEntities parents(*refinedMesh->get_bulk_data());
    for (unsigned irank=0; irank < ranks_to_be_deleted.size()-1; irank++) // don't look for children of nodes
    {
      const stk::mesh::BucketVector & buckets = refinedMesh->get_bulk_data()->buckets( ranks_to_be_deleted[irank] );
      int npar=0;
      int nchild=0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        stk::mesh::Bucket & bucket = **k ;

        const unsigned num_elements_in_bucket = bucket.size();
        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
        {
	  stk::mesh::Entity element = bucket[iElement];
	  if (!refinedMesh->isParentElement(element, false))
	  {
	    ++nchild;
	  }
     	  else
	  {
	    ++npar;
	    //parents.insert(element);
	  }
        }
      }
      TEST_ASSERT(npar>0);
      TEST_ASSERT(nchild>0);
    }
  }


  {
    out << "\nCreating pamgen mesh without parent elements." << std::endl;
    RCP<ParameterList> p = parameterList();
    const std::string input_file_name = "pamgen_test.gen";
    p->set("File Name",input_file_name);
    p->set("File Type","Pamgen");
    p->set("Levels of Uniform Refinement",1);
    p->set("Keep Percept Data",true);

    RCP<STK_ExodusReaderFactory> pamgenFactory = rcp(new STK_ExodusReaderFactory());
    TEST_NOTHROW(pamgenFactory->setParameterList(p));

    RCP<STK_Interface> mesh = pamgenFactory->buildMesh(MPI_COMM_WORLD);
    RCP<percept::PerceptMesh> refinedMesh = mesh->getRefinedMesh();

    // check if the Percept data exists and has elements
    TEST_ASSERT(nonnull(refinedMesh));
    TEST_ASSERT(refinedMesh->get_number_elements()>0);

    // check if the parent information is stored
    std::vector<stk::mesh::EntityRank> ranks_to_be_deleted;
    ranks_to_be_deleted.push_back(stk::topology::ELEMENT_RANK);
    ranks_to_be_deleted.push_back(refinedMesh->side_rank());
    if (refinedMesh->get_spatial_dim() == 3)
      ranks_to_be_deleted.push_back(refinedMesh->edge_rank());

    //percept::SetOfEntities parents(*refinedMesh->get_bulk_data());
    for (unsigned irank=0; irank < ranks_to_be_deleted.size()-1; irank++) // don't look for children of nodes
    {
      const stk::mesh::BucketVector & buckets = refinedMesh->get_bulk_data()->buckets( ranks_to_be_deleted[irank] );
      int npar=0;
      int nchild=0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        stk::mesh::Bucket & bucket = **k ;

        const unsigned num_elements_in_bucket = bucket.size();
        for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
        {
	  stk::mesh::Entity element = bucket[iElement];
	  if (!refinedMesh->isParentElement(element, false))
	  {
	    ++nchild;
	  }
     	  else
	  {
	    ++npar;
	    //parents.insert(element);
	  }
        }
      }
      TEST_EQUALITY(npar,0);
      TEST_ASSERT(nchild>0);
    }
  }
}
#endif

} // namespace panzer_stk

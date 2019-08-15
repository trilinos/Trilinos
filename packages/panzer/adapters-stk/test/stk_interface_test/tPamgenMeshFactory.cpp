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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/CreateAdjacentEntities.hpp"

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
    broker->add_mesh_database("pamgen_test.gen", "pamgen", stk::io::READ_MESH);
    broker->create_input_mesh();
    metaData = broker->meta_data_rcp();
    bulkData = Teuchos::rcp(new stk::mesh::BulkData(*metaData,MPI_COMM_WORLD));
    broker->set_bulk_data(bulkData);
    broker->add_all_mesh_fields_as_input_fields();
    broker->populate_bulk_data();

    if (verbose) {
      metaData->dump_all_meta_info(out);
      bulkData->dump_all_mesh_info(out);
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
    broker->set_bulk_data(bulkData);
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
    broker->add_mesh_database(output_exodus_file_name, "exodus", stk::io::READ_MESH);
    broker->create_input_mesh();
    metaData = broker->meta_data_rcp();
    bulkData = Teuchos::rcp(new stk::mesh::BulkData(*metaData,MPI_COMM_WORLD));
    broker->set_bulk_data(bulkData);
    broker->add_all_mesh_fields_as_input_fields();
    broker->populate_bulk_data();

    if (verbose) {
      metaData->dump_all_meta_info(out);
      bulkData->dump_all_mesh_info(out);
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
      metaData->dump_all_meta_info(out);
      bulkData->dump_all_mesh_info(out);
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
  TEST_EQUALITY(panzer_stk::getMeshDimension("pamgen_test.gen",MPI_COMM_WORLD,false),3);
}

} // namespace panzer_stk

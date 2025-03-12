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

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#include <limits> // for epsilon for tolerance in FP comparisons

#ifdef PANZER_HAVE_IOSS

// for checking test correctness
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_ElementBlock.h"          // for ElementBlock
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, etc

const std::string short_field_name    = "short_element_variable_name";
const std::string long_field_name     = "long_element_variable_name_XXXXXXXXXXXXX";
const std::string too_long_field_name = "too_long_element_variable_name_XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
const std::string block_name_1 = "block_1";

namespace panzer_stk {

  void setCellValues(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                     std::string field_name,
                     double value)
  {
    auto cell_vals = mesh->getCellField(field_name,block_name_1);

    auto meta_data = mesh->getMetaData();
    auto bulk_data = mesh->getBulkData();

    auto cell_buckets = bulk_data->get_buckets(stk::topology::ELEM_RANK,meta_data->locally_owned_part());
    for (size_t bucket_i=0 ; bucket_i<cell_buckets.size() ; ++bucket_i) {
      stk::mesh::Bucket &cell_bucket = *cell_buckets[bucket_i];
      for (size_t cell_i=0 ; cell_i<cell_bucket.size() ; ++cell_i) {
        double * cell_data = stk::mesh::field_data(*cell_vals,cell_bucket.bucket_id(),cell_i);
        *cell_data = value + 0.5;
      }
    }
  }

TEUCHOS_UNIT_TEST(tExodusIossProperties, DefaultNameLength)
{
  int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
  int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

  stk::ParallelMachine pm(MPI_COMM_WORLD);
  const std::string restart_file = "exodus_ioss_properties_32chars.exo";
  const auto tolerance = std::numeric_limits<double>::epsilon() * 100.0;

  {
    STK_ExodusReaderFactory f("meshes/basic.gen");
    auto mesh = f.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addCellField(short_field_name,block_name_1);
    mesh->addCellField(long_field_name,block_name_1);
    mesh->addCellField(too_long_field_name,block_name_1);

    mesh->initialize(pm,true,false);

    f.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    mesh->setupExodusFile(restart_file);
    setCellValues(mesh,short_field_name,0.0);
    setCellValues(mesh,long_field_name,1.0);
    setCellValues(mesh,too_long_field_name,2.0);
    mesh->writeToExodus(0.0);
  }
  // *****************************
  // Read the mesh using ioss directly and confirm that the field names
  // are no longer than 32 characters long
  // *****************************
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", restart_file, Ioss::READ_MODEL, MPI_COMM_WORLD);
    Ioss::Region results(resultsDb);

    // Number of time steps
    TEST_EQUALITY(results.get_property("state_count").get_int(),1);

    // First time step value (ioss steps are 1 based)
    int step = 1;
    double db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,0.0,tolerance);
    results.end_state(step);

    // Cell fields exist
    Ioss::ElementBlock *eb = results.get_element_blocks()[0];
    TEST_EQUALITY(eb->field_count(Ioss::Field::TRANSIENT), 5);
    TEST_ASSERT(true  == eb->field_exists(short_field_name));
    TEST_ASSERT(false == eb->field_exists(long_field_name));
    TEST_ASSERT(false == eb->field_exists(too_long_field_name));
    TEST_ASSERT(true  == eb->field_exists(long_field_name.substr(0,32)));
    TEST_ASSERT(true  == eb->field_exists(too_long_field_name.substr(0,32)));
  }
}

TEUCHOS_UNIT_TEST(tExodusIossProperties, LongNameLength)
{
  int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
  int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

  stk::ParallelMachine pm(MPI_COMM_WORLD);
  const std::string restart_file = "exodus_ioss_properties_64chars.exo";
  const auto tolerance = std::numeric_limits<double>::epsilon() * 100.0;

  {
    STK_ExodusReaderFactory f("meshes/basic.gen");
    auto mesh = f.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addCellField(short_field_name,block_name_1);
    mesh->addCellField(long_field_name,block_name_1);
    mesh->addCellField(too_long_field_name,block_name_1);

    mesh->initialize(pm,true,false);

    f.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    std::vector<Ioss::Property> ioss_properties;
    ioss_properties.push_back(Ioss::Property("MAXIMUM_NAME_LENGTH", 64));

    mesh->setupExodusFile(restart_file, ioss_properties);
    setCellValues(mesh,short_field_name,0.0);
    setCellValues(mesh,long_field_name,1.0);
    setCellValues(mesh,too_long_field_name,2.0);
    mesh->writeToExodus(0.0);
  }
  // *****************************
  // Read the mesh using ioss directly and confirm that the field names
  // are no longer than 64 characters long
  // *****************************
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", restart_file, Ioss::READ_MODEL, MPI_COMM_WORLD);
    Ioss::Region results(resultsDb);

    // Number of time steps
    TEST_EQUALITY(results.get_property("state_count").get_int(),1);

    // First time step value (ioss steps are 1 based)
    int step = 1;
    double db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,0.0,tolerance);
    results.end_state(step);

    // Cell fields exist
    Ioss::ElementBlock *eb = results.get_element_blocks()[0];
    TEST_EQUALITY(eb->field_count(Ioss::Field::TRANSIENT), 5);
    TEST_ASSERT(true  == eb->field_exists(short_field_name));
    TEST_ASSERT(true  == eb->field_exists(long_field_name));
    TEST_ASSERT(false == eb->field_exists(too_long_field_name));
    TEST_ASSERT(true  == eb->field_exists(too_long_field_name.substr(0,64)));
  }
}

} // namespace panzer_stk

#endif

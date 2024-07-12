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

const std::string node_field_name = "DUMMY_NODAL_FIELD";
const std::string cell_field_name = "DUMMY_CELL_FIELD";
const std::string block_name_1 = "block_1";
const std::string block_name_2 = "block_2";

namespace panzer_stk {

  void setValues(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                 double value)
  {
    auto node_vals = mesh->getSolutionField(node_field_name,block_name_1);
    auto cell_vals = mesh->getCellField(cell_field_name,block_name_1);

    auto meta_data = mesh->getMetaData();
    auto bulk_data = mesh->getBulkData();

    auto node_buckets = bulk_data->get_buckets(stk::topology::NODE_RANK,meta_data->locally_owned_part());
    for (size_t bucket_i=0 ; bucket_i<node_buckets.size() ; ++bucket_i) {
      stk::mesh::Bucket &node_bucket = *node_buckets[bucket_i];
      for (size_t node_i=0 ; node_i<node_bucket.size() ; ++node_i) {
        double * node_data = stk::mesh::field_data(*node_vals,node_bucket.bucket_id(),node_i);
        *node_data = value;
      }
    }

    auto cell_buckets = bulk_data->get_buckets(stk::topology::ELEM_RANK,meta_data->locally_owned_part());
    for (size_t bucket_i=0 ; bucket_i<cell_buckets.size() ; ++bucket_i) {
      stk::mesh::Bucket &cell_bucket = *cell_buckets[bucket_i];
      for (size_t cell_i=0 ; cell_i<cell_bucket.size() ; ++cell_i) {
        double * cell_data = stk::mesh::field_data(*cell_vals,cell_bucket.bucket_id(),cell_i);
        *cell_data = value + 0.5;
      }
    }
  }

TEUCHOS_UNIT_TEST(tExodusRestartAndAppend, Restart)
{
  int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
  int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  out << "Running numprocs = " << numprocs << " rank = " << rank << std::endl;

  stk::ParallelMachine pm(MPI_COMM_WORLD);
  const std::string restart_file = "exodus_restart_and_append.exo";
  const auto tolerance = std::numeric_limits<double>::epsilon() * 100.0;

  // *****************************
  // Read in mesh and write out some fields for three time steps
  // *****************************
  {
    STK_ExodusReaderFactory f("meshes/basic.gen");
    auto mesh = f.buildUncommitedMesh(MPI_COMM_WORLD);

    mesh->addSolutionField(node_field_name,block_name_1);
    mesh->addSolutionField(node_field_name,block_name_2);
    mesh->addCellField(cell_field_name,block_name_1);
    mesh->addCellField(cell_field_name,block_name_2);

    mesh->initialize(pm,true,false);

    f.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    mesh->setupExodusFile(restart_file);
    setValues(mesh,0.0);
    mesh->writeToExodus(0.0);
    setValues(mesh,1.0);
    mesh->writeToExodus(1.0);
    setValues(mesh,2.0);
    mesh->writeToExodus(2.0);
  }

  // *****************************
  // Restart from mesh and append more time steps using original
  // construction path. This appends to the end of the file.
  // *****************************
  {
    STK_ExodusReaderFactory f(restart_file);
    auto mesh = f.buildUncommitedMesh(MPI_COMM_WORLD);
    mesh->addSolutionField(node_field_name,block_name_1);
    mesh->addSolutionField(node_field_name,block_name_2);
    mesh->addCellField(cell_field_name,block_name_1);
    mesh->addCellField(cell_field_name,block_name_2);
    mesh->initialize(pm,true,false);
    f.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    bool appendToFile = true;
    mesh->setupExodusFile(restart_file,appendToFile);
    setValues(mesh,3.0);
    mesh->writeToExodus(3.0);
    setValues(mesh,4.0);
    mesh->writeToExodus(4.0);
    setValues(mesh,5.0);
    mesh->writeToExodus(5.0);
  }

  // *****************************
  // Read the final mesh using ioss directly and make sure the time
  // steps were appended instead of overwritten
  // *****************************
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", restart_file, Ioss::READ_MODEL, MPI_COMM_WORLD);
    Ioss::Region results(resultsDb);

    // Number of time steps
    TEST_EQUALITY(results.get_property("state_count").get_int(),6);

    // First time step value (ioss steps are 1 based)
    int step = 1;
    double db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,0.0,tolerance);
    results.end_state(step);

    // Last time step value
    step = 6;
    db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,5.0,tolerance);
    results.end_state(step);

    // Nodal field exists
    Ioss::NodeBlock *nb = results.get_node_blocks()[0];
    TEST_EQUALITY(nb->field_count(Ioss::Field::TRANSIENT), 1);
    TEST_ASSERT(nb->field_exists(node_field_name));

    // Cell field exists
    Ioss::ElementBlock *eb = results.get_element_blocks()[0];
    TEST_EQUALITY(eb->field_count(Ioss::Field::TRANSIENT), 3);
    TEST_ASSERT(eb->field_exists(cell_field_name));
  }

  // *****************************
  // Now test the "append after time value". This allows a user to
  // back up and restart earlier in the run. Instead of appending to
  // the end, this overwrites any time steps after the specified time
  // and appends after that.
  // *****************************
  {
    STK_ExodusReaderFactory f(restart_file);
    auto mesh = f.buildUncommitedMesh(MPI_COMM_WORLD);
    mesh->addSolutionField(node_field_name,block_name_1);
    mesh->addSolutionField(node_field_name,block_name_2);
    mesh->addCellField(cell_field_name,block_name_1);
    mesh->addCellField(cell_field_name,block_name_2);
    mesh->initialize(pm,true,false);
    f.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    bool appendToFile = true;
    bool appendAfterRestartTime = true;
    double restartTime = 3.0;
    mesh->setupExodusFile(restart_file,appendToFile,appendAfterRestartTime,restartTime);
    setValues(mesh,3.3);
    mesh->writeToExodus(3.3);
    setValues(mesh,3.4);
    mesh->writeToExodus(3.4);
    setValues(mesh,3.5);
    mesh->writeToExodus(3.5);
    setValues(mesh,3.6);
    mesh->writeToExodus(3.6);
  }

  // *****************************
  // Verify the "append after time" worked
  // *****************************
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", restart_file, Ioss::READ_MODEL, MPI_COMM_WORLD);
    Ioss::Region results(resultsDb);

    // Number of time steps
    TEST_EQUALITY(results.get_property("state_count").get_int(),8);

    // First time step value (ioss steps are 1 based)
    int step = 1;
    double db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,0.0,tolerance);
    results.end_state(step);

    // Middle time step value
    step = 5;
    db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,3.3,tolerance);
    results.end_state(step);

    // Last time step value
    step = 8;
    db_time = results.begin_state(step);
    TEST_FLOATING_EQUALITY(db_time,3.6,tolerance);
    results.end_state(step);

    // Nodal field exists
    Ioss::NodeBlock *nb = results.get_node_blocks()[0];
    TEST_EQUALITY(nb->field_count(Ioss::Field::TRANSIENT), 1);
    TEST_ASSERT(nb->field_exists(node_field_name));

    // Cell field exists
    Ioss::ElementBlock *eb = results.get_element_blocks()[0];
    TEST_EQUALITY(eb->field_count(Ioss::Field::TRANSIENT), 3);
    TEST_ASSERT(eb->field_exists(cell_field_name));
  }

}

} // namespace panzer_stk

#endif

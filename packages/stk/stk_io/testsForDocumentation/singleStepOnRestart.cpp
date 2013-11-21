#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ioss_SubSystem.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <fieldNameTestUtils.hpp>
#include <restartTestUtils.hpp>

namespace {
  TEST(StkMeshIoBrokerHowTo, singleStepOnRestart)
  {
    // ============================================================
    // INITIALIZATION...

    std::string filename = "single_step.restart";
    MPI_Comm comm = MPI_COMM_WORLD;

    stk::io::StkMeshIoBroker stkIo(comm);
      
    stk::mesh::MetaData &meta = generateMetaData(stkIo);
    stk::mesh::FieldBase *field = declareTriStateNodalField(meta, "disp");
    stkIo.populate_bulk_data();

    putDataOnTriStateField(stkIo.bulk_data(), field, 1.0, 2.0, 3.0);

    {
      //-BEGIN
      // ... Setup deleted
      // ============================================================
      // EXAMPLE USAGE...
      // Create a restart file,
      size_t fh = stkIo.create_output_mesh(filename,
					   stk::io::WRITE_RESTART);
      stkIo.add_field(fh, *field);

      //+ Set the cycle count to 1.  This will result in a maximum
      //+ of one step on the output database -- when a new step is
      //+ added, it will overwrite the existing step.
      stkIo.get_output_io_region(fh)->get_database()->set_cycle_count(1); /*@\label{io:cycle}*/
	
      // Write multiple steps to the restart file.
      stkIo.begin_output_step(fh, 0.0);
      stkIo.write_defined_output_fields(fh);
      stkIo.end_output_step(fh);

      stkIo.begin_output_step(fh, 1.0);
      stkIo.write_defined_output_fields(fh);
      stkIo.end_output_step(fh);

      stkIo.begin_output_step(fh, 2.0);
      stkIo.write_defined_output_fields(fh);
      stkIo.end_output_step(fh);

      //+ At this point, there should only be a single state on the
      //+ restart database. The time of this state should be 2.0.
      // ... Verification deleted
      //-END
    }

    {
      // ============================================================
      //  VERIFICATION:
      //+ There should only be a single state on the restart database.
      //+ The time of this state should be 2.0.
      Ioss::DatabaseIO *iossDb =
	Ioss::IOFactory::create("exodus", filename,
				Ioss::READ_RESTART, comm);
      Ioss::Region region(iossDb);

      EXPECT_EQ(region.get_property("state_count").get_int(), 1);
      EXPECT_EQ(region.get_state_time(1), 2.0);
    }
    
    // =============================================================
    // CLEANUP:
    unlink(filename.c_str());
  }
}

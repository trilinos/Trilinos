/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace {
  // Do the actual reading and writing of the mesh database and
  // creation and population of the MetaData and BulkData.
  void mesh_read_write(const std::string &type,
		       const std::string &working_directory,
		       const std::string &filename,
		       stk::io::MeshData &mesh_data,
		       int db_integer_size)
  {
    std::string file = working_directory;
    file += filename;

    mesh_data.open_mesh_database(file, type);
    mesh_data.create_input_mesh();

    // This is done just to define some fields in stk
    // that can be used later for reading restart data.
    mesh_data.define_input_fields();
    mesh_data.populate_bulk_data();


    // Iterate all fields and set them as restart fields...
    const stk::mesh::FieldVector &fields = mesh_data.meta_data().get_fields();
    for (size_t i=0; i < fields.size(); i++) {
      std::string name = fields[i]->name();
      mesh_data.add_restart_field(*fields[i], name);
    }

    // ========================================================================
    // Create output mesh...  ("generated_mesh.out") ("exodus_mesh.out")
    std::string output_filename = working_directory + type + "_mesh.out";

    // This outputs the "mesh" portion of the output results file.
    // This consists of the coordinates, connectivity, nodesets and sidesets. 
    // Basically all non-transient (time-dependent) data.
    mesh_data.create_output_mesh(output_filename);

    // Create restart output ...  ("generated_mesh.restart") ("exodus_mesh.restart")
    std::string restart_filename = working_directory + type + "_mesh.restart";

    mesh_data.create_restart_output(restart_filename);

    // Determine the names of the global fields on the input
    // mesh. These will be used below to define the same fields on the
    // restart and results output databases.
    std::vector<std::string> global_fields;
    mesh_data.get_global_variable_names(global_fields);

    // For each global field name on the input database, determine the type of the field
    // and define that same global field on both the results and restart output databases.
    for (size_t i=0; i < global_fields.size(); i++) {
      const Ioss::Field &input_field = mesh_data.input_io_region()->get_fieldref(global_fields[i]);

      // Define the global fields that will be written on each timestep.
      mesh_data.add_restart_global(input_field.get_name(), input_field.raw_storage()->name(), input_field.get_type());
      mesh_data.add_results_global(input_field.get_name(), input_field.raw_storage()->name(), input_field.get_type());
    }

    // ========================================================================
    // Begin the transient loop...  All timesteps on the input database are transferred
    // to the results and restart output databases...

    // Determine number of timesteps on input database...
    int timestep_count = mesh_data.input_io_region()->get_property("state_count").get_int();
    for (int step=1; step <= timestep_count; step++) {
      double time = mesh_data.input_io_region()->get_state_time(step);

      // Normally, an app would only process the restart input at a single step and
      // then continue with execution at that point.  Here just for testing, we are
      // reading restart data at each step on the input restart file/mesh and then
      // outputting that data to the restart and results output.

      mesh_data.process_restart_input(step);
      mesh_data.begin_restart_output_at_time(time);
      mesh_data.begin_results_output_at_time(time);

      mesh_data.process_restart_output();
      mesh_data.process_output_request();

      // Transfer all global variables from the input mesh to the
      // restart and results databases
      for (size_t i=0; i < global_fields.size(); i++) {
	std::vector<double> field_values;
        mesh_data.get_global(global_fields[i], field_values);
        mesh_data.write_restart_global(global_fields[i], field_values);
        mesh_data.write_results_global(global_fields[i], field_values);
      }

      mesh_data.end_current_restart_output();
      mesh_data.end_current_results_output();
    }
  }

  void driver(stk::ParallelMachine  comm,
	      const std::string &parallel_io,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type,
	      const std::string &decomp_method,
	      bool compose_output,
	      int  compression_level,
	      bool compression_shuffle,
	      int  db_integer_size)
  {
    stk::io::MeshData mesh_data(comm);

    bool use_netcdf4 = false;
    if (!decomp_method.empty()) {
      mesh_data.m_property_manager.add(Ioss::Property("DECOMPOSITION_METHOD", decomp_method));
    }

    if (compose_output) {
      mesh_data.m_property_manager.add(Ioss::Property("COMPOSE_RESULTS", true));
      mesh_data.m_property_manager.add(Ioss::Property("COMPOSE_RESTART", true));
    }
    if (!parallel_io.empty()) {
      mesh_data.m_property_manager.add(Ioss::Property("PARALLEL_IO_MODE", parallel_io));
    }
    if (compression_level > 0) {
      mesh_data.m_property_manager.add(Ioss::Property("COMPRESSION_LEVEL", compression_level));
      use_netcdf4 = true;
    }
    if (compression_shuffle) {
      mesh_data.m_property_manager.add(Ioss::Property("COMPRESSION_SHUFFLE", 1));
      use_netcdf4 = true;
    }
    if (use_netcdf4) {
      mesh_data.m_property_manager.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }
    if (db_integer_size == 8) {
      mesh_data.m_property_manager.add(Ioss::Property("INTEGER_SIZE_DB", db_integer_size));
      mesh_data.m_property_manager.add(Ioss::Property("INTEGER_SIZE_API", db_integer_size));
    }

    mesh_read_write(type, working_directory, filename, mesh_data, db_integer_size);
  }
}

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  std::string working_directory = "";
  std::string decomp_method = "";
  std::string mesh = "";
  std::string type = "exodusii";
  int compression_level = 0;
  bool compression_shuffle = false;
  int db_integer_size = 4;
  bool compose_output = false;
  std::string parallel_io = "";

  //----------------------------------
  // Process the broadcast command line arguments
  bopt::options_description desc("options");

  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
  desc.add_options()
    ("directory,d",   bopt::value<std::string>(&working_directory),
     "working directory with trailing '/'" )
    ("decomposition,D", bopt::value<std::string>(&decomp_method),
     "decomposition method.  One of: linear, rcb, rib, hsfc, block, cyclic, random, kway, geom_kway, metis_sfc" )
    ("mesh",          bopt::value<std::string>(&mesh),
     "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Can also specify a filename. The generated mesh will be output to the file 'generated_mesh.out'" )
    ("compression_level", bopt::value<int>(&compression_level), "compression level [1..9] to use" )
    ("shuffle", bopt::value<bool>(&compression_shuffle), "use shuffle filter prior to compressing data: true|false" )
    ("compose_output", bopt::value<bool>(&compose_output), "create a single output file: true|false" )
    ("parallel_io_method", bopt::value<std::string>(&parallel_io),
     "Method to use for parallel io. One of mpiio, mpiposix, or pnetcdf")
    ("db_integer_size", bopt::value<int>(&db_integer_size), "use 4 or 8-byte integers on output database" );


  stk::get_options_description().add(desc);

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  if (mesh.empty()) {
    std::cerr << "\nERROR: The --mesh option is required\n";
    std::cerr << "\nApplication " << desc << "\n";
    std::exit(EXIT_FAILURE);
  }

  type = "exodusii";
  if (strncasecmp("gen:", mesh.c_str(), 4) == 0) {
    mesh = mesh.substr(4, mesh.size());
    type = "generated";
  }
  if (strncasecmp("dof:", mesh.c_str(), 4) == 0) {
    mesh = mesh.substr(4, mesh.size());
    type = "dof";
  }
  driver(use_case_environment.m_comm, parallel_io,
	 working_directory, mesh, type, decomp_method, compose_output, 
	 compression_level, compression_shuffle, db_integer_size);

  return 0;
}


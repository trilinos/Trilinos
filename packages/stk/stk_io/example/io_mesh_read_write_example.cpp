// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <string>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ParameterList.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace {
  // Do the actual reading and writing of the mesh database and
  // creation and population of the MetaData and BulkData.
  void mesh_read_write(const std::string &type,
		       const std::string &working_directory,
		       const std::string &filename,
		       stk::io::StkMeshIoBroker &mesh_data,
		       int integer_size,
		       stk::io::HeartbeatType hb_type,
		       int interpolation_intervals)
  {
    if (interpolation_intervals == 0)
      interpolation_intervals = 1;
    
    std::string file = working_directory;
    file += filename;

    size_t input_index = mesh_data.add_mesh_database(file, type, stk::io::READ_MESH);
    mesh_data.set_active_mesh(input_index);
    mesh_data.create_input_mesh();

    // This is done just to define some fields in stk
    // that can be used later for reading restart data.
    stk::io::MeshField::TimeMatchOption tmo = stk::io::MeshField::CLOSEST;
    if (interpolation_intervals > 1) {
      tmo = stk::io::MeshField::LINEAR_INTERPOLATION;
    }
    mesh_data.add_all_mesh_fields_as_input_fields(tmo);

    mesh_data.populate_bulk_data();

    // ========================================================================
    // Create output mesh...  ("generated_mesh.out") ("exodus_mesh.out")
    std::string output_filename = working_directory + type + "_mesh.out";

    // This call adds an output database for results data to mesh_data.
    // No data is written at this time other than verifying that the
    // file can be created on the disk.
    size_t results_index = mesh_data.create_output_mesh(output_filename, stk::io::WRITE_RESULTS);

    // Create restart output ...  ("generated_mesh.restart") ("exodus_mesh.restart")
    std::string restart_filename = working_directory + type + "_mesh.restart";

    size_t restart_index = mesh_data.create_output_mesh(restart_filename, stk::io::WRITE_RESTART);

    // Iterate all fields and set them as restart fields...
    const stk::mesh::FieldVector &fields = mesh_data.meta_data().get_fields();
    for (size_t i=0; i < fields.size(); i++) {
      const Ioss::Field::RoleType* role = stk::io::get_field_role(*fields[i]);
      if ( role && *role == Ioss::Field::TRANSIENT ) {
	mesh_data.add_field(restart_index, *fields[i]); // restart output
	mesh_data.add_field(results_index, *fields[i]); // results output
      }
    }

    // Determine the names of the global fields on the input
    // mesh. These will be used below to define the same fields on the
    // restart and results output databases.
    std::vector<std::string> global_fields;
    mesh_data.get_global_variable_names(global_fields);

    // Create heartbeat file of the specified format...
    size_t heart = 0;
    if (hb_type != stk::io::NONE && !global_fields.empty()) {
      std::string heartbeat_filename = working_directory + type + ".hrt";
      heart = mesh_data.add_heartbeat_output(heartbeat_filename, hb_type);
    }
    
    stk::util::ParameterList parameters;
    
    // For each global field name on the input database, determine the type of the field
    // and define that same global field on the results, restart, history, and heartbeat outputs.
    if (!global_fields.empty()) {
      std::cout << "Adding " << global_fields.size() << " global fields:\n";
    }

    Teuchos::RCP<Ioss::Region> io_region = mesh_data.get_input_io_region();
    STKIORequire(!Teuchos::is_null(io_region));
      
    for (size_t i=0; i < global_fields.size(); i++) {
      const Ioss::Field &input_field = io_region->get_fieldref(global_fields[i]);
      std::cout << "\t" << input_field.get_name() << " of type " << input_field.raw_storage()->name() << "\n";

      if (input_field.raw_storage()->component_count() == 1) {
	double val = 0.0;
	parameters.set_param(input_field.get_name(), val);
      }
      else {
	std::vector<double> vals(input_field.raw_storage()->component_count());
	parameters.set_param(input_field.get_name(), vals);
      }

      // Define the global fields that will be written on each timestep.
      mesh_data.add_global(restart_index, input_field.get_name(),
			   input_field.raw_storage()->name(), input_field.get_type());
      mesh_data.add_global(results_index, input_field.get_name(),
			   input_field.raw_storage()->name(), input_field.get_type());
      if (hb_type != stk::io::NONE) {
	stk::util::Parameter &param = parameters.get_param(input_field.get_name());
	mesh_data.add_heartbeat_global(heart, input_field.get_name(), &param.value, param.type);
      }
    }

    // ========================================================================
    // Begin the transient loop...  All timesteps on the input database are transferred
    // to the results and restart output databases...

    // Determine number of timesteps on input database...
    int timestep_count = io_region->get_property("state_count").get_int();

    if (timestep_count == 0 ) {
      mesh_data.write_output_mesh(results_index);
    }
    else {
      for (int step=1; step <= timestep_count; step++) {
	double time = io_region->get_state_time(step);
	if (step == timestep_count)
	  interpolation_intervals = 1;
	
	int step_end = step < timestep_count ? step+1 : step;
	double tend =  io_region->get_state_time(step_end);
	double tbeg = time;
	double delta = (tend - tbeg) / static_cast<double>(interpolation_intervals);
	
	for (int interval = 0; interval < interpolation_intervals; interval++) {
	  // Normally, an app would only process the restart input at a single step and
	  // then continue with execution at that point.  Here just for testing, we are
	  // reading restart data at each step on the input restart file/mesh and then
	  // outputting that data to the restart and results output.
	  time = tbeg + delta * static_cast<double>(interval);

	  mesh_data.read_defined_input_fields(time);
	  mesh_data.begin_output_step(restart_index, time);
	  mesh_data.begin_output_step(results_index, time);

	  mesh_data.write_defined_output_fields(restart_index);
	  mesh_data.write_defined_output_fields(results_index);

	  // Transfer all global variables from the input mesh to the
	  // restart and results databases
	  stk::util::ParameterMapType::const_iterator i = parameters.begin();
	  stk::util::ParameterMapType::const_iterator iend = parameters.end();
	  for (; i != iend; ++i) {
	    const std::string parameterName = (*i).first;
	    stk::util::Parameter &parameter = parameters.get_param(parameterName);
	    mesh_data.get_global(parameterName, parameter.value, parameter.type);
	  }

	  for (i=parameters.begin(); i != iend; ++i) {
	    const std::string parameterName = (*i).first;
	    stk::util::Parameter parameter = (*i).second;
	    mesh_data.write_global(restart_index, parameterName, parameter.value, parameter.type);
	    mesh_data.write_global(results_index, parameterName, parameter.value, parameter.type);
	  }

	  mesh_data.end_output_step(restart_index);
	  mesh_data.end_output_step(results_index);

	}
	if (hb_type != stk::io::NONE && !global_fields.empty()) {
	  mesh_data.process_heartbeat_output(heart, step, time);
	}
      }
    }
  }

  void driver(const std::string &parallel_io,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type,
	      const std::string &decomp_method,
	      bool compose_output,
	      int  compression_level,
	      bool compression_shuffle,
	      int  integer_size,
	      stk::io::HeartbeatType hb_type,
	      int interpolation_intervals)
  {
    stk::io::StkMeshIoBroker mesh_data(MPI_COMM_WORLD);

    bool use_netcdf4 = false;
    if (!decomp_method.empty()) {
      mesh_data.property_add(Ioss::Property("DECOMPOSITION_METHOD", decomp_method));
    }

    if (compose_output) {
      mesh_data.property_add(Ioss::Property("COMPOSE_RESULTS", true));
      mesh_data.property_add(Ioss::Property("COMPOSE_RESTART", true));
    }
    if (!parallel_io.empty()) {
      mesh_data.property_add(Ioss::Property("PARALLEL_IO_MODE", parallel_io));
    }
    if (compression_level > 0) {
      mesh_data.property_add(Ioss::Property("COMPRESSION_LEVEL", compression_level));
      use_netcdf4 = true;
    }
    if (compression_shuffle) {
      mesh_data.property_add(Ioss::Property("COMPRESSION_SHUFFLE", 1));
      use_netcdf4 = true;
    }
    if (use_netcdf4) {
      mesh_data.property_add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }
    if (integer_size == 8) {
      mesh_data.property_add(Ioss::Property("INTEGER_SIZE_DB", integer_size));
      mesh_data.property_add(Ioss::Property("INTEGER_SIZE_API", integer_size));
    }

    mesh_read_write(type, working_directory, filename, mesh_data, integer_size, hb_type,
		    interpolation_intervals);
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
  int interpolation_intervals = 0;
  bool compression_shuffle = false;
  int integer_size = 4;
  bool compose_output = false;
  std::string parallel_io = "";
  std::string heartbeat_format = "none";
  //----------------------------------
  // Process the broadcast command line arguments
  bopt::options_description desc("options");

  desc.add_options()
    ("help,h", "produce help message")
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
    ("heartbeat_format", bopt::value<std::string>(&heartbeat_format),
     "Format of heartbeat output. One of binary, csv, text, ts_text, spyhis, [none]")
    ("interpolate", bopt::value<int>(&interpolation_intervals), "number of intervals to divide each input time step into")
    ("integer_size", bopt::value<int>(&integer_size), "use 4 or 8-byte integers for input and output" );


  stk::parallel_machine_init(&argc, &argv);

  bopt::variables_map vm;
  bopt::store(bopt::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  bopt::notify(vm);

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

  stk::io::HeartbeatType hb_type = stk::io::NONE; // Default is no heartbeat output
  if (heartbeat_format == "none")
    hb_type = stk::io::NONE;
  else if (heartbeat_format == "binary")
    hb_type = stk::io::BINARY;
  else if (heartbeat_format == "csv")
    hb_type = stk::io::CSV;
  else if (heartbeat_format == "ts_csv")
    hb_type = stk::io::TS_CSV;
  else if (heartbeat_format == "text")
    hb_type = stk::io::TEXT;
  else if (heartbeat_format == "ts_text")
    hb_type = stk::io::TS_TEXT;
  else if (heartbeat_format == "spyhis")
    hb_type = stk::io::SPYHIS;

  driver(parallel_io,
	 working_directory, mesh, type, decomp_method, compose_output, 
	 compression_level, compression_shuffle, integer_size, hb_type,
	 interpolation_intervals);

  stk::parallel_machine_finalize();
  return 0;
}


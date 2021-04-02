// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <stk_util/Version.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/environment/LogWithTimeAndMemory.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/environment/memory_util.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ParameterList.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

class IoMeshDriver
{
public:
  IoMeshDriver(MPI_Comm comm)
  : m_comm(comm), m_hwmAvg_baseline(0)
  {
    set_output_streams();
    size_t hwmMax = 0, hwmMin = 0, hwmAvg = 0;
    stk::get_current_memory_usage_across_processors(m_comm, hwmMax, hwmMin, hwmAvg);
    m_hwmAvg_baseline = hwmAvg;
  }

  void set_output_streams()
  {
    if (stk::parallel_machine_rank(m_comm) != 0) {
      stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    }
    Ioss::Utils::set_output_stream(sierra::Env::outputP0());
  }

  void log_mesh_counts(const stk::mesh::BulkData& mesh)
  {
    std::vector<size_t> globalCounts;
    std::vector<size_t> minGlobalCounts;
    std::vector<size_t> maxGlobalCounts;
    std::vector<size_t> auraGlobalCounts;
    std::vector<size_t> sharedNotOwnedCounts;
    stk::mesh::comm_mesh_counts(mesh, globalCounts, minGlobalCounts, maxGlobalCounts);
    stk::mesh::Selector sharedNotOwned = mesh.mesh_meta_data().globally_shared_part() & !mesh.mesh_meta_data().locally_owned_part();
    stk::mesh::count_entities(sharedNotOwned, mesh, sharedNotOwnedCounts);
    constexpr unsigned numRanks = static_cast<unsigned>(stk::topology::ELEM_RANK+1);
    stk::all_reduce(m_comm, stk::ReduceSum<numRanks>(sharedNotOwnedCounts.data()));
    stk::mesh::Selector aura = mesh.mesh_meta_data().aura_part();
    stk::mesh::count_entities(aura, mesh, auraGlobalCounts);
    stk::all_reduce(m_comm, stk::ReduceSum<numRanks>(auraGlobalCounts.data()));

    stk::log_with_time_and_memory(m_comm, " - Elements: "+std::to_string(globalCounts[stk::topology::ELEM_RANK])
                                          +" total Aura: "+std::to_string(auraGlobalCounts[stk::topology::ELEM_RANK]));
    stk::log_with_time_and_memory(m_comm, " - Nodes: "+std::to_string(globalCounts[stk::topology::NODE_RANK])
                                          +", shared-not-owned: "+std::to_string(sharedNotOwnedCounts[stk::topology::NODE_RANK])
                                          +", total Aura: "+std::to_string(auraGlobalCounts[stk::topology::NODE_RANK]));
    size_t hwmMax = 0, hwmMin = 0, hwmAvg = 0;
    stk::get_current_memory_usage_across_processors(m_comm, hwmMax, hwmMin, hwmAvg);
    size_t totalBytes = mesh.parallel_size() * (hwmAvg - m_hwmAvg_baseline);
    size_t bytesPerElem = totalBytes/globalCounts[stk::topology::ELEM_RANK];
    stk::log_with_time_and_memory(m_comm, "Max HWM per proc: "+stk::human_bytes(hwmMax)
                                         + ", bytes-per-element: " + std::to_string(bytesPerElem));
  }

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

    stk::log_with_time_and_memory(m_comm, "Finished populating input mesh, aura is "
          +std::string((mesh_data.bulk_data().is_automatic_aura_on() ? "on" : "off")));
    log_mesh_counts(mesh_data.bulk_data());

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

    auto io_region = mesh_data.get_input_io_region();
      
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

	// Flush the data.  This is not necessary in a normal
	// application, Just being done here to verify that the
	// function exists and does not core dump.  
	mesh_data.flush_output();
      }
    }

    stk::log_with_time_and_memory(m_comm, "Finished writing output mesh.");
  }

  void driver(const std::string &parallel_io,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type,
	      const std::string &decomp_method,
	      bool compose_output,
	      int  compression_level,
	      bool compression_shuffle,
	      bool lower_case_variable_names,
	      int  integer_size,
	      stk::io::HeartbeatType hb_type,
	      int interpolation_intervals)
  {
    std::string readOrCreate = ((type=="generated" || type=="pamgen") ? "Creating" : "Reading");
    stk::log_with_time_and_memory(m_comm, readOrCreate+" input mesh: "+filename);

    stk::io::StkMeshIoBroker mesh_data(MPI_COMM_WORLD);

    mesh_data.property_add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", lower_case_variable_names));

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

private:
  MPI_Comm m_comm;
  size_t m_hwmAvg_baseline;
};

int main(int argc, const char** argv)
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
  bool lc_names = true;
  std::string parallel_io = "";
  std::string heartbeat_format = "none";

  MPI_Comm comm = stk::parallel_machine_init(&argc, const_cast<char***>(&argv));
  int myRank = stk::parallel_machine_rank(comm);

  //----------------------------------
  // Process the command line arguments
  stk::CommandLineParserParallel cmdLine(comm);

  cmdLine.add_required<std::string>({"mesh", "m", "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Can also specify a filename. The generated mesh will be output to the file 'generated_mesh.out'"});
  cmdLine.add_optional<std::string>({"directory", "d", "working directory with trailing '/'"}, "./");
  cmdLine.add_optional<std::string>({"decomposition", "D", "decomposition method.  One of: linear, rcb, rib, hsfc, block, cyclic, random, kway, geom_kway, metis_sfc"}, "rcb");
  cmdLine.add_optional<int>({"compression_level", "c", "compression level [1..9] to use"}, 0);
  cmdLine.add_optional<std::string>({"shuffle", "s", "use shuffle filter prior to compressing data: true|false"}, "false");
  cmdLine.add_optional<std::string>({"lower_case_variable_names", "l", "convert variable names to lowercase and replace spaces in names with underscore (default is true): true|false"}, "true");
  cmdLine.add_optional<std::string>({"compose_output", "C", "create a single output file: true|false"}, "false");
  cmdLine.add_optional<std::string>({"parallel_io_method", "p", "Method to use for parallel io. One of mpiio, mpiposix, or pnetcdf"}, "pnetcdf");
  cmdLine.add_optional<std::string>({"heartbeat_format", "h", "Format of heartbeat output. One of binary, csv, text, ts_text, spyhis, [none]"}, "none");
  cmdLine.add_optional<int>({"interpolate", "i", "number of intervals to divide each input time step into"}, 0);
  cmdLine.add_optional<int>({"integer_size", "I", "use 4 or 8-byte integers for input and output"}, 4);

  stk::CommandLineParser::ParseState parseState = cmdLine.parse(argc, argv);

  if (parseState != stk::CommandLineParser::ParseComplete) {
    int returnCode = 0;
    switch(parseState) {
    case stk::CommandLineParser::ParseError:
      returnCode = 1;
      break;
    case stk::CommandLineParser::ParseHelpOnly:
      std::cout << cmdLine.get_usage() << std::endl;
      break;
    case stk::CommandLineParser::ParseVersionOnly:
      if (myRank == 0) {
        std::cout << "STK Version: " << stk::version_string() << std::endl;
      }
      break;
    default: break;
    }
    stk::parallel_machine_finalize();
    return returnCode;
  }

  IoMeshDriver ioMeshDriver(comm);

  if (cmdLine.is_option_provided("directory")) {
    working_directory = cmdLine.get_option_value<std::string>("directory");
  }
  if (cmdLine.is_option_parsed("decomposition")) {
    decomp_method = cmdLine.get_option_value<std::string>("decomposition");
  }
  if (cmdLine.is_option_provided("shuffle")) {
    if (cmdLine.get_option_value<std::string>("shuffle") == "true") {
      compression_shuffle = true;
    }
  }
  if (cmdLine.is_option_provided("lower_case_variable_names")) {
    if (cmdLine.get_option_value<std::string>("lower_case_variable_names") == "false") {
      lc_names = false;
    }
  }
  if (cmdLine.is_option_provided("compose_output")) {
    if (cmdLine.get_option_value<std::string>("compose_output") == "true" ||
        cmdLine.get_option_value<std::string>("compose_output") == "yes") {
      compose_output = true;
    }
  }
  if (cmdLine.is_option_provided("compression_level")) {
    compression_level = cmdLine.get_option_value<int>("compression_level");
  }
  if (cmdLine.is_option_provided("parallel_io_method")) {
    parallel_io = cmdLine.get_option_value<std::string>("parallel_io_method");
  }
  if (cmdLine.is_option_provided("heartbeat_format")) {
    heartbeat_format = cmdLine.get_option_value<std::string>("heartbeat_format");
  }
  if (cmdLine.is_option_provided("interpolate")) {
    interpolation_intervals = cmdLine.get_option_value<int>("interpolate");
  }
  if (cmdLine.is_option_provided("integer_size")) {
    integer_size = cmdLine.get_option_value<int>("integer_size");
  }

  mesh = cmdLine.get_option_value<std::string>("mesh");

  type = "exodusii";
  if (strncasecmp("gen:", mesh.c_str(), 4) == 0) {
    mesh = mesh.substr(4, mesh.size());
    type = "generated";
  }
  if (strncasecmp("generated:", mesh.c_str(), 10) == 0) {
    mesh = mesh.substr(10, mesh.size());
    type = "generated";
  }
  else if (strncasecmp("dof:", mesh.c_str(), 4) == 0) {
    mesh = mesh.substr(4, mesh.size());
    type = "dof";
  }
  else if (strncasecmp("cgns:", mesh.c_str(), 5) == 0) {
    mesh = mesh.substr(5, mesh.size());
    type = "cgns";
  }
  else if (strncasecmp("pamgen:", mesh.c_str(), 7) == 0) {
    mesh = mesh.substr(7, mesh.size());
    type = "pamgen";
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

  ioMeshDriver.driver(parallel_io,
	 working_directory, mesh, type, decomp_method, compose_output, 
	 compression_level, compression_shuffle, lc_names, integer_size, hb_type,
	 interpolation_intervals);

  stk::parallel_machine_finalize();
  return 0;
}


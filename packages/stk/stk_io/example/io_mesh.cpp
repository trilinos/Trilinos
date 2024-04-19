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

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <string.h>                                             // for size_t
#include <algorithm>                                            // for copy
#include <iostream>                                             // for opera...
#include <map>                                                  // for _Rb_t...
#include <stk_io/IossBridge.hpp>                                // for get_f...
#include <stk_io/StkMeshIoBroker.hpp>                           // for StkMe...
#include <stk_mesh/base/Comm.hpp>                               // for comm_...
#include <stk_mesh/base/GetEntities.hpp>                        // for count...
#include <stk_mesh/base/MetaData.hpp>                           // for MetaData
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_util/Version.hpp>                                 // for versi...
#include <stk_util/command_line/CommandLineParserParallel.hpp>  // for Comma...
#include <stk_util/environment/Env.hpp>                         // for outputP0
#include <stk_util/environment/EnvData.hpp>                     // for EnvData
#include <stk_util/environment/LogWithTimeAndMemory.hpp>        // for log_w...
#include <stk_util/environment/memory_util.hpp>                 // for get_c...
#include <stk_util/parallel/Parallel.hpp>                       // for paral...
#include <stk_util/parallel/ParallelReduce.hpp>                 // for Reduc...
#include <stk_util/util/ParameterList.hpp>                      // for Param...
#include <stk_util/util/human_bytes.hpp>                        // for human...
#include <stk_util/util/MemoryTracking.hpp>
#include <string>                                               // for string
#include <utility>                                              // for pair
#include <vector>                                               // for vector
#include "Ioss_Field.h"                                         // for Field
#include "Ioss_Property.h"                                      // for Property
#include "Ioss_Region.h"                                        // for Region
#include "Ioss_Utils.h"                                         // for Utils
#include "Ioss_VariableType.h"                                  // for Varia...
#include "mpi.h"                                                // for MPI_Comm
#include "stk_io/DatabasePurpose.hpp"                           // for READ_...
#include "stk_io/Heartbeat.hpp"                                 // for NONE
#include "stk_io/MeshField.hpp"                                 // for MeshF...
#include "stk_mesh/base/BulkData.hpp"                           // for BulkData
#include "stk_mesh/base/MeshBuilder.hpp"                        // for MeshBuilder
#include "stk_mesh/base/Part.hpp"                               // for Part
#include "stk_mesh/base/Selector.hpp"                           // for Selector
#include "stk_mesh/base/Types.hpp"                              // for Field...
#include "stk_topology/topology.hpp"                            // for topology
#include "stk_util/command_line/CommandLineParser.hpp"          // for Comma...
#include "stk_util/util/SimpleArrayOps.hpp"                     // for Sum
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

void set_output_streams(MPI_Comm comm)
{
  if (stk::parallel_machine_rank(comm) != 0) {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
  }
  Ioss::Utils::set_output_stream(sierra::Env::outputP0());
}

class IoMeshDriver
{
public:
  IoMeshDriver(MPI_Comm comm)
  : m_comm(comm), m_curAvg_baseline(0), m_baselineBuffer()
  {
    equilibrate_memory_baseline();
  }

  void equilibrate_memory_baseline()
  {
    size_t now = 0, hwm = 0;
    stk::get_memory_usage(now, hwm);
    if (hwm > now) {
      m_baselineBuffer.resize((hwm-now)/sizeof(double));
    }

    size_t curMax = 0, curMin = 0, curAvg = 0;
    stk::get_memory_high_water_mark_across_processors(m_comm, curMax, curMin, curAvg);
    m_curAvg_baseline = curAvg;
  }

  void log_mesh_counts(const stk::mesh::BulkData& mesh)
  {
    size_t curMax = 0, curMin = 0, curAvg = 0;
    stk::get_memory_high_water_mark_across_processors(m_comm, curMax, curMin, curAvg);
    size_t totalBytes = mesh.parallel_size() * (curAvg - m_curAvg_baseline);

    constexpr unsigned numRanks = static_cast<unsigned>(stk::topology::ELEM_RANK+1);
    std::vector<size_t> globalCounts(numRanks, 0);
    std::vector<size_t> minGlobalCounts(numRanks, 0);
    std::vector<size_t> maxGlobalCounts(numRanks, 0);
    std::vector<size_t> auraGlobalCounts(numRanks, 0);
    std::vector<size_t> sharedNotOwnedCounts(numRanks, 0);
    stk::mesh::comm_mesh_counts(mesh, globalCounts, minGlobalCounts, maxGlobalCounts);
    stk::mesh::Selector sharedNotOwned = mesh.mesh_meta_data().globally_shared_part() & !mesh.mesh_meta_data().locally_owned_part();
    stk::mesh::count_entities(sharedNotOwned, mesh, sharedNotOwnedCounts);
    stk::all_reduce(m_comm, stk::ReduceSum<numRanks>(sharedNotOwnedCounts.data()));
    stk::mesh::Selector aura = mesh.mesh_meta_data().aura_part();
    stk::mesh::count_entities(aura, mesh, auraGlobalCounts);
    stk::all_reduce(m_comm, stk::ReduceSum<numRanks>(auraGlobalCounts.data()));

    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(mesh.mesh_meta_data().entity_rank_count());
    {
      std::ostringstream os;
      os<<std::setw(8)<<" "<<std::setw(12)<<"owned"<<std::setw(14)<<"sh-not-owned"<<std::setw(10)<<"aura"<<std::endl;
      for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; ++rank) {
        os<<std::setw(34+8)<<(mesh.mesh_meta_data().entity_rank_names()[rank])
          <<std::setw(12)<<globalCounts[rank]
          <<std::setw(12)<<sharedNotOwnedCounts[rank]
          <<std::setw(12)<<auraGlobalCounts[rank]
          <<std::endl;
      }
      stk::log_with_time_and_memory(m_comm, os.str());
    }

    size_t totalEntities = 0;
    for(size_t count : globalCounts) totalEntities += count;
    for(size_t count : sharedNotOwnedCounts) totalEntities += count;
    for(size_t count : auraGlobalCounts) totalEntities += count;

    size_t bytesPerEntity = totalEntities>0 ? totalBytes/totalEntities : 0;
    std::string bytesPerEntityStr = totalEntities>0 ? std::to_string(bytesPerEntity) : std::string("N/A");
    stk::log_with_time_and_memory(m_comm, "Total HWM Mesh Memory: "+stk::human_bytes(totalBytes));
    stk::log_with_time_and_memory(m_comm, "Total Mesh Entities: "+std::to_string(totalEntities)
                                         + ", bytes-per-entity: " + bytesPerEntityStr);

    for(stk::mesh::EntityRank rank=stk::topology::NODE_RANK; rank<endRank; ++rank) {
      size_t localBucketCapacity = 0, localBucketSize = 0;
      const stk::mesh::BucketVector& buckets = mesh.buckets(rank);
      for(const stk::mesh::Bucket* bptr : buckets) {
        localBucketCapacity += bptr->capacity();
        localBucketSize += bptr->size();
      }
      
      size_t globalNumBuckets = stk::get_global_sum(m_comm, buckets.size());
      size_t globalBucketCapacity = stk::get_global_sum(m_comm, localBucketCapacity);
      size_t globalBucketSize = stk::get_global_sum(m_comm, localBucketSize);
      std::ostringstream os;
      os<<globalNumBuckets<<" "<<rank<<" buckets, total size/capacity: "<<globalBucketSize<<" / "<<globalBucketCapacity;
      if (globalNumBuckets > 0) {
        double totalSz = globalBucketSize;
        double proportion = totalSz / globalBucketCapacity;
        os<<"; "<<(100*proportion)<<"%";
      }
      stk::log_with_time_and_memory(m_comm, os.str());
    }
    
#ifdef STK_MEMORY_TRACKING
    size_t localBytes = stk::get_total_bytes_currently_allocated();
    size_t globalBytes = stk::get_global_sum(m_comm, localBytes);
    size_t localPtrs = stk::get_current_num_ptrs();
    size_t globalPtrs = stk::get_global_sum(m_comm, localPtrs);
    stk::log_with_time_and_memory(m_comm, "Total tracked bytes: "+stk::human_bytes(globalBytes)+", num ptrs: "+std::to_string(globalPtrs));
    size_t localHWMBytes = stk::get_high_water_mark_in_bytes();
    size_t globalHWMBytes = stk::get_global_sum(m_comm, localHWMBytes);
    size_t localHWMPtrs = stk::get_high_water_mark_in_ptrs();
    size_t globalHWMPtrs = stk::get_global_sum(m_comm, localHWMPtrs);
    stk::log_with_time_and_memory(m_comm, "Total HWM tracked bytes: "+std::to_string(globalHWMBytes)+", HWM num ptrs: "+std::to_string(globalHWMPtrs));
#endif
  }

  void mesh_read(const std::string &type,
                 const std::string &working_directory,
                 const std::string &filename,
                 stk::io::StkMeshIoBroker &ioBroker,
                 int integer_size,
                 stk::io::HeartbeatType hb_type,
                 int interpolation_intervals)
  {
#ifdef STK_MEMORY_TRACKING
    stk::reset_high_water_mark_in_bytes();
    stk::reset_high_water_mark_in_ptrs();
#endif

//    constexpr unsigned N = 1000000;
//    std::vector<double> v(N, 0.0);
    if (interpolation_intervals == 0)
      interpolation_intervals = 1;
    
    std::string file = working_directory;
    file += filename;

    size_t input_index = ioBroker.add_mesh_database(file, type, stk::io::READ_MESH);
    ioBroker.set_active_mesh(input_index);
    ioBroker.create_input_mesh();

    // This is done just to define some fields in stk
    // that can be used later for reading restart data.
    stk::io::MeshField::TimeMatchOption tmo = stk::io::MeshField::CLOSEST;
    if (interpolation_intervals > 1) {
      tmo = stk::io::MeshField::LINEAR_INTERPOLATION;
    }
    ioBroker.add_all_mesh_fields_as_input_fields(tmo);

    ioBroker.populate_bulk_data();

    if (m_addEdges) {
      stk::mesh::create_edges(ioBroker.bulk_data());
    }

    if (m_addFaces) {
      stk::mesh::create_faces(ioBroker.bulk_data());
    }
  }

  void mesh_write(const std::string &type,
		       const std::string &working_directory,
		       const std::string &filename,
		       stk::io::StkMeshIoBroker &ioBroker,
		       int integer_size,
		       stk::io::HeartbeatType hb_type,
		       int interpolation_intervals)
  {
    if (interpolation_intervals == 0)
      interpolation_intervals = 1;
    
    std::string file = working_directory;
    file += filename;

    // ========================================================================
    // Create output mesh...  ("generated_mesh.out") ("exodus_mesh.out")
    std::string output_filename = working_directory + type + "_mesh.out";

    // This call adds an output database for results data to ioBroker.
    // No data is written at this time other than verifying that the
    // file can be created on the disk.
    size_t results_index = ioBroker.create_output_mesh(output_filename, stk::io::WRITE_RESULTS);

    // Create restart output ...  ("generated_mesh.restart") ("exodus_mesh.restart")
    std::string restart_filename = working_directory + type + "_mesh.restart";

    size_t restart_index = ioBroker.create_output_mesh(restart_filename, stk::io::WRITE_RESTART);

    // Iterate all fields and set them as restart fields...
    const stk::mesh::FieldVector &fields = ioBroker.meta_data().get_fields();
    for (size_t i=0; i < fields.size(); i++) {
      const Ioss::Field::RoleType* role = stk::io::get_field_role(*fields[i]);
      if ( role && *role == Ioss::Field::TRANSIENT ) {
	ioBroker.add_field(restart_index, *fields[i]); // restart output
	ioBroker.add_field(results_index, *fields[i]); // results output
      }
    }

    // Determine the names of the global fields on the input
    // mesh. These will be used below to define the same fields on the
    // restart and results output databases.
    std::vector<std::string> global_fields;
    ioBroker.get_global_variable_names(global_fields);

    // Create heartbeat file of the specified format...
    size_t heart = 0;
    if (hb_type != stk::io::NONE && !global_fields.empty()) {
      std::string heartbeat_filename = working_directory + type + ".hrt";
      heart = ioBroker.add_heartbeat_output(heartbeat_filename, hb_type);
    }
    
    stk::util::ParameterList parameters;
    
    // For each global field name on the input database, determine the type of the field
    // and define that same global field on the results, restart, history, and heartbeat outputs.
    if (!global_fields.empty()) {
      sierra::Env::outputP0() << "Adding " << global_fields.size() << " global fields:\n";
    }

    auto io_region = ioBroker.get_input_ioss_region();
      
    for (size_t i=0; i < global_fields.size(); i++) {
      const Ioss::Field &input_field = io_region->get_fieldref(global_fields[i]);
      sierra::Env::outputP0() << "\t" << input_field.get_name() << " of type " << input_field.raw_storage()->name() << "\n";

      if (input_field.raw_storage()->component_count() == 1) {
	double val = 0.0;
	parameters.set_param(input_field.get_name(), val);
      }
      else {
	std::vector<double> vals(input_field.raw_storage()->component_count());
	parameters.set_param(input_field.get_name(), vals);
      }

      // Define the global fields that will be written on each timestep.
      ioBroker.add_global(restart_index, input_field.get_name(),
			   input_field.raw_storage()->name(), input_field.get_type());
      ioBroker.add_global(results_index, input_field.get_name(),
			   input_field.raw_storage()->name(), input_field.get_type());
      if (hb_type != stk::io::NONE) {
          stk::util::Parameter &param = parameters.get_param(input_field.get_name());
          ioBroker.add_heartbeat_global(heart, input_field.get_name(), param);
      }
    }

    // ========================================================================
    // Begin the transient loop...  All timesteps on the input database are transferred
    // to the results and restart output databases...

    // Determine number of timesteps on input database...
    int timestep_count = io_region->get_property("state_count").get_int();

    if (timestep_count == 0 ) {
      ioBroker.write_output_mesh(results_index);
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

	  ioBroker.read_defined_input_fields(time);
	  ioBroker.begin_output_step(restart_index, time);
	  ioBroker.begin_output_step(results_index, time);

	  ioBroker.write_defined_output_fields(restart_index);
	  ioBroker.write_defined_output_fields(results_index);

	  // Transfer all global variables from the input mesh to the
	  // restart and results databases
	  stk::util::ParameterMapType::const_iterator i = parameters.begin();
	  stk::util::ParameterMapType::const_iterator iend = parameters.end();
	  for (; i != iend; ++i) {
	    const std::string parameterName = (*i).first;
	    stk::util::Parameter &parameter = parameters.get_param(parameterName);
	    ioBroker.get_global(parameterName, parameter);
	  }

	  for (i=parameters.begin(); i != iend; ++i) {
	    const std::string parameterName = (*i).first;
	    stk::util::Parameter parameter = (*i).second;
	    ioBroker.write_global(restart_index, parameterName, parameter.value, parameter.type);
	    ioBroker.write_global(results_index, parameterName, parameter.value, parameter.type);
	  }

	  ioBroker.end_output_step(restart_index);
	  ioBroker.end_output_step(results_index);

	}
	if (hb_type != stk::io::NONE && !global_fields.empty()) {
	  ioBroker.process_heartbeat_output(heart, step, time);
	}

	// Flush the data.  This is not necessary in a normal
	// application, Just being done here to verify that the
	// function exists and does not core dump.  
	ioBroker.flush_output();
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
    stk::mesh::BulkData::AutomaticAuraOption aura = m_auraOption ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA;
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_aura_option(aura);
    builder.set_upward_connectivity(m_upwardConnectivity);
    builder.set_initial_bucket_capacity(m_initialBucketCapacity);
    builder.set_maximum_bucket_capacity(m_maximumBucketCapacity);
    stk::log_with_time_and_memory(m_comm, "Creating MetaData/BulkData objects");
    std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

    stk::log_with_time_and_memory(m_comm, "Creating StkMeshIoBroker object");

    stk::io::StkMeshIoBroker ioBroker(MPI_COMM_WORLD);
    set_io_properties(ioBroker, lower_case_variable_names, decomp_method, compose_output, parallel_io, compression_level, compression_shuffle, integer_size);

    stk::log_with_time_and_memory(m_comm, "Setting memory baseline");
    equilibrate_memory_baseline();
    stk::log_with_time_and_memory(m_comm, "Finished setting memory baseline");

    stk::log_with_time_and_memory(m_comm, "Reading input mesh: "+filename);

    ioBroker.set_bulk_data(bulk);
    mesh_read(type, working_directory, filename, ioBroker, integer_size, hb_type, interpolation_intervals);

    stk::log_with_time_and_memory(m_comm, "Finished reading input mesh");

    log_mesh_counts(*bulk);

    mesh_write(type, working_directory, filename, ioBroker, integer_size, hb_type, interpolation_intervals);
  }

  void set_io_properties(stk::io::StkMeshIoBroker& ioBroker,
                         bool lower_case_variable_names,
                         const std::string& decomp_method,
                         bool compose_output,
                         const std::string& parallel_io,
                         int compression_level,
                         bool compression_shuffle,
                         int integer_size)
  {
    ioBroker.property_add(Ioss::Property("LOWER_CASE_VARIABLE_NAMES", lower_case_variable_names));

    if (!decomp_method.empty()) {
      ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", decomp_method));
    }
    else {
      sierra::Env::outputP0()<<"decomposition not specified, defaulting to file-per-processor mode for mesh-read."<<std::endl;
    }

    if (compose_output) {
      ioBroker.property_add(Ioss::Property("COMPOSE_RESULTS", true));
      ioBroker.property_add(Ioss::Property("COMPOSE_RESTART", true));
    }

    if (!parallel_io.empty()) {
      ioBroker.property_add(Ioss::Property("PARALLEL_IO_MODE", parallel_io));
    }

    bool use_netcdf4 = false;
    if (compression_level > 0) {
      ioBroker.property_add(Ioss::Property("COMPRESSION_LEVEL", compression_level));
      use_netcdf4 = true;
    }
    if (compression_shuffle) {
      ioBroker.property_add(Ioss::Property("COMPRESSION_SHUFFLE", 1));
      use_netcdf4 = true;
    }
    if (use_netcdf4) {
      ioBroker.property_add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    if (integer_size == 8) {
      ioBroker.property_add(Ioss::Property("INTEGER_SIZE_DB", integer_size));
      ioBroker.property_add(Ioss::Property("INTEGER_SIZE_API", integer_size));
    }
  }

  void set_add_edges(bool trueOrFalse) { m_addEdges = trueOrFalse; }
  void set_add_faces(bool trueOrFalse) { m_addFaces = trueOrFalse; }
  void set_upward_connectivity(bool trueOrFalse) { m_upwardConnectivity = trueOrFalse; }
  void set_aura_option(bool trueOrFalse) { m_auraOption = trueOrFalse; }
  void set_initial_bucket_capacity(int initialCapacity) { m_initialBucketCapacity = initialCapacity; }
  void set_maximum_bucket_capacity(int maximumCapacity) { m_maximumBucketCapacity = maximumCapacity; }

private:
  MPI_Comm m_comm;
  size_t m_curAvg_baseline;
  bool m_addEdges;
  bool m_addFaces;
  bool m_upwardConnectivity;
  bool m_auraOption;
  int m_initialBucketCapacity;
  int m_maximumBucketCapacity;
  std::vector<double> m_baselineBuffer;
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
  bool addEdges = false;
  bool addFaces = false;
  bool upwardConnectivity = true;
  bool auraOption = false;
  std::string parallel_io = "";
  std::string heartbeat_format = "none";
  int initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity();
  int maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity();

  MPI_Comm comm = stk::parallel_machine_init(&argc, const_cast<char***>(&argv));

  set_output_streams(comm);

  //----------------------------------
  // Process the command line arguments
  stk::CommandLineParserParallel cmdLine(comm);

  cmdLine.add_required<std::string>({"mesh", "m", "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Can also specify a filename. The generated mesh will be output to the file 'generated_mesh.out'"});
  cmdLine.add_optional<std::string>({"add_edges", "e", "create all internal edges in the mesh: true|false"}, "false");
  cmdLine.add_optional<std::string>({"add_faces", "f", "create all internal faces in the mesh: true|false"}, "false");
  cmdLine.add_optional<std::string>({"upward_connectivity", "u", "create upward connectivity/adjacency in the mesh (default is true): true|false"}, "true");
  cmdLine.add_optional<std::string>({"ghost_aura", "g", "create aura ghosting around each MPI rank (default is false): true|false"}, "false");
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
  cmdLine.add_optional<int>({"initial_bucket_capacity", "b", "initial bucket capacity"}, stk::mesh::get_default_initial_bucket_capacity());
  cmdLine.add_optional<int>({"maximum_bucket_capacity", "B", "maximum bucket capacity"}, stk::mesh::get_default_maximum_bucket_capacity());

  stk::CommandLineParser::ParseState parseState = cmdLine.parse(argc, argv);

  if (parseState != stk::CommandLineParser::ParseComplete) {
    int returnCode = 0;
    switch(parseState) {
    case stk::CommandLineParser::ParseError:
      returnCode = 1;
      break;
    case stk::CommandLineParser::ParseHelpOnly:
      sierra::Env::outputP0() << cmdLine.get_usage() << std::endl;
      break;
    case stk::CommandLineParser::ParseVersionOnly:
      sierra::Env::outputP0() << "STK Version: " << stk::version_string() << std::endl;
      break;
    default: break;
    }
    stk::parallel_machine_finalize();
    return returnCode;
  }

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
  if (cmdLine.is_option_provided("add_edges")) {
    if (cmdLine.get_option_value<std::string>("add_edges") == "true") {
      addEdges = true;
    }
  }
  if (cmdLine.is_option_provided("add_faces")) {
    if (cmdLine.get_option_value<std::string>("add_faces") == "true") {
      addFaces = true;
    }
  }
  if (cmdLine.is_option_provided("upward_connectivity")) {
    if (cmdLine.get_option_value<std::string>("upward_connectivity") == "false") {
      upwardConnectivity = false;
    }
  }
  if (cmdLine.is_option_provided("ghost_aura")) {
    if (cmdLine.get_option_value<std::string>("ghost_aura") == "true") {
      auraOption = true;
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
  if (cmdLine.is_option_provided("initial_bucket_capacity")) {
    initialBucketCapacity = cmdLine.get_option_value<int>("initial_bucket_capacity");
  }
  if (cmdLine.is_option_provided("maximum_bucket_capacity")) {
    maximumBucketCapacity = cmdLine.get_option_value<int>("maximum_bucket_capacity");
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

  IoMeshDriver ioMeshDriver(comm);

  ioMeshDriver.set_add_edges(addEdges);
  ioMeshDriver.set_add_faces(addFaces);
  ioMeshDriver.set_upward_connectivity(upwardConnectivity);
  ioMeshDriver.set_aura_option(auraOption);
  ioMeshDriver.set_initial_bucket_capacity(initialBucketCapacity);
  ioMeshDriver.set_maximum_bucket_capacity(maximumBucketCapacity);

  ioMeshDriver.driver(parallel_io,
	 working_directory, mesh, type, decomp_method, compose_output, 
	 compression_level, compression_shuffle, lc_names, integer_size, hb_type,
	 interpolation_intervals);

  stk::parallel_machine_finalize();
  return 0;
}


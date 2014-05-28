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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace {
  void driver(stk::ParallelMachine  comm,
              size_t dimension,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type,
	      const std::string &decomp_method,
	      int  compression_level,
	      bool compression_shuffle,
	      int  db_integer_size);
}

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  std::string working_directory = "";
  std::string decomp_method = "";
  std::string mesh = "";
  std::string type = "exodusii";
  size_t spatial_dimension = 3;
  int compression_level = 0;
  bool compression_shuffle = false;
  int db_integer_size = 4;


  //----------------------------------
  // Process the broadcast command line arguments
  bopt::options_description desc("options");

  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
  desc.add_options()
    ("directory,d",   bopt::value<std::string>(&working_directory),
     "working directory with trailing '/'" )
    ("decomposition,D", bopt::value<std::string>(&decomp_method),
     "decomposition method" )
    ("mesh",          bopt::value<std::string>(&mesh),
     "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Can also specify a filename. The generated mesh will be output to the file 'generated_mesh.out'" )
    ("dimension", bopt::value<size_t>(&spatial_dimension), "problem spatial dimension" )
    ("compression_level", bopt::value<int>(&compression_level), "compression level [1..9] to use" )
    ("shuffle", bopt::value<bool>(&compression_shuffle), "use shuffle filter prior to compressing data" )
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
  driver(use_case_environment.m_comm, spatial_dimension,
	 working_directory, mesh, type, decomp_method,
	 compression_level, compression_shuffle, db_integer_size);

  return 0;
}

namespace {
  void driver(stk::ParallelMachine  comm,
              size_t spatial_dimension,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type,
	      const std::string &decomp_method,
	      int  compression_level,
	      bool compression_shuffle,
	      int  db_integer_size)
  {

    // Initialize IO system.  Registers all element types and storage
    // types and the exodusII default database type.
    Ioss::Init::Initializer init_db;

    stk::mesh::fem::FEMMetaData fem_meta_data( spatial_dimension );
    stk::mesh::MetaData &meta_data = fem_meta_data.get_meta_data( fem_meta_data );
    stk::io::MeshData mesh_data;

    bool use_netcdf4 = false;
    if (!decomp_method.empty()) {
      mesh_data.m_property_manager.add(Ioss::Property("DECOMPOSITION_METHOD", decomp_method));
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
    }
      
    std::string file = working_directory;
    file += filename;
    stk::io::create_input_mesh(type, file, comm, fem_meta_data, mesh_data);
    stk::io::define_input_fields(mesh_data, fem_meta_data);

    fem_meta_data.commit();
    stk::mesh::BulkData bulk_data(meta_data , comm);
    stk::io::populate_bulk_data(bulk_data, mesh_data);

    //------------------------------------------------------------------
    // Create output mesh...  ("generated_mesh.out") ("exodus_mesh.out")
    std::string output_filename = working_directory + type + "_mesh.out";
    stk::io::create_output_mesh(output_filename, comm, bulk_data, mesh_data);
    stk::io::define_output_fields(mesh_data, fem_meta_data);

    // Determine number of timesteps on input database...
    int timestep_count = mesh_data.m_input_region->get_property("state_count").get_int();
    for (int step=1; step <= timestep_count; step++) {
      double time = mesh_data.m_input_region->get_state_time(step);
      stk::io::process_input_request(mesh_data, bulk_data, step);
      stk::io::process_output_request(mesh_data, bulk_data, time);
    }
  }
}

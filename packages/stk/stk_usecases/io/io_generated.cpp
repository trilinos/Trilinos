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
	      const std::string &type);
}

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  std::string working_directory = "";
  std::string mesh = "";
  std::string type = "exodusii";
  size_t spatial_dimension = 3;


  //----------------------------------
  // Process the broadcast command line arguments
  bopt::options_description desc("options");

  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
  desc.add_options()
    ("directory,d",   bopt::value<std::string>(&working_directory),
     "working directory with trailing '/'" )
    ("mesh",          bopt::value<std::string>(&mesh),
     "mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Can also specify a filename. The generated mesh will be output to the file 'generated_mesh.out'" )
    ("dimension", bopt::value<size_t>(&spatial_dimension), "problem spatial dimension" );

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
  driver(use_case_environment.m_comm, spatial_dimension,
	 working_directory, mesh, type);

  return 0;
}

namespace {
  void driver(stk::ParallelMachine  comm,
              size_t spatial_dimension,
	      const std::string &working_directory,
	      const std::string &filename,
	      const std::string &type)
  {

    // Initialize IO system.  Registers all element types and storage
    // types and the exodusII default database type.
    Ioss::Init::Initializer init_db;

    stk::mesh::fem::FEMMetaData fem_meta_data( spatial_dimension );
    stk::mesh::MetaData &meta_data = fem_meta_data.get_meta_data( fem_meta_data );
    stk::io::MeshData mesh_data;

    std::string file = working_directory;
    file += filename;
    stk::io::create_input_mesh(type, file, comm, fem_meta_data, mesh_data);

    fem_meta_data.commit();
    stk::mesh::BulkData bulk_data(meta_data , comm);
    stk::io::populate_bulk_data(bulk_data, mesh_data);

    //------------------------------------------------------------------
    // Create output mesh...  ("generated_mesh.out") ("exodus_mesh.out")
    std::string output_filename = working_directory + type + "_mesh.out";
    stk::io::create_output_mesh(output_filename, comm, bulk_data, mesh_data);
  }
}

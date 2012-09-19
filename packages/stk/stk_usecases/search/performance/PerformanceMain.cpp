/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string>
#include <iostream>

#include <search/performance/Performance.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{
  stk::search::Options range;
  stk::search::Options domain;

  std::string working_directory = "";
  bool performance              = false;
  
  /** \todo IMPLEMENT -- what should the syntax and options for this look like? */
  std::string search_type       = "";

  
  //----------------------------------
  // Process the broadcast command line arguments
  
  bopt::options_description desc("options");
    
  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
  desc.add_options()
    ("directory,d",   bopt::value<std::string>(&working_directory),
     "working directory with trailing '/'" )
    ("search_type",   bopt::value<std::string>(&search_type),
     "type of search to run. Valid options are:\n"
     "  point_in_box: \tdomain consists of an axis-aligned bounding box for each 'domain_entity' in the mesh,"
     "range is a PointBoundingBox3D at the centroid of each'range_entity'\n"
     "  overlap:      \tdomain consists of an axis-aligned bounding box for each 'domain_entity' in the mesh,"
     "range is also an axis-aligned bounding box of each 'range_entity'.")
    ("range_mesh",    bopt::value<std::string>(&range.mesh_filename),
     "range mesh file.\n  \tUse name of form 'gen:NxMxL' to generate a hex mesh of size N by M by L intervals.\n\tUse --helpmesh for more detailed help on the mesh options.")
    ("range_entity",  bopt::value<std::string>(&range.entity)->default_value("element"),
     "Entity type to use for range:\n"
     "   node:    \tAll nodes in the model,\n"
     "   nodeset: \tAll nodes in a nodeset,\n"
     "   edge:    \tAll edges in the model (probably none),\n"
     "   face:    \tSkin the model and use all exposed faces,\n"
     "   faceset: \tFaces, but only those faces in a faceset,\n"
     "   element: \tAll elements in the model.")
    ("range_offset",
     bopt::value<double>(&range.offset)->default_value(0.0),
     "offset to be applied to domain axis-aligned bounding boxes. The bounding box extent will be increased by this amount in each direction." )
    ("range_scale",
     bopt::value<double>(&range.scale)->default_value(0.0),
     "scale factor to be applied to domain axis-aligned bounding boxes. The bounding box extent will be increased by max_d*scale+offset where max_d is max of max_i-min_i for i=x,y,z)" )
    ("domain_mesh",   bopt::value<std::string>(&domain.mesh_filename),
     "domain mesh file."
     "\n  \tUse name of form 'gen:NxMxL' to generate a hex mesh of size N by M by L intervals."
     "\n  \tUse --helpmesh for more detailed help on the mesh options." )
    ("domain_entity", bopt::value<std::string>(&domain.entity)->default_value("element"),
     "Entity type to use for domain:\n"
     "   node:    \tAll nodes in the model,\n"
     "   nodeset: \tAll nodes in a nodeset,\n"
     "   edge:    \tAll edges in the model (probably none),\n"
     "   face:    \tSkin the model and use all exposed faces,\n"
     "   faceset: \tFaces, but only those faces in a faceset,\n"
     "   element: \tAll elements in the model.")
    ("domain_offset",
     bopt::value<double>(&domain.offset)->default_value(0.0),
     "offset to be applied to domain axis-aligned bounding boxes. The bounding box extent will be increased by this amount in each direction." )
    ("domain_scale",
     bopt::value<double>(&domain.scale)->default_value(0.0),
     "scale factor to be applied to domain axis-aligned bounding boxes. The bounding box extent will be increased by max_d*scale+offset where max_d is max of max_i-min_i for i=x,y,z)" )
    ("performance",  "Run to measure performance; disable output that may affect runtime.")
    ("helpmesh", "Print detailed description of mesh options and then exit.");
  stk::get_options_description().add(desc);

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  bopt::variables_map &vm = stk::get_variables_map();  

  //----------------------------------

  if (vm.count("helpmesh")) {
    stk::io::show_mesh_help();
    std::exit(EXIT_SUCCESS);
  }
  
  if (vm.count("performance"))
    performance = true;
  
  if (search_type.empty()) {
    std::cerr << "\nOPTION ERROR: A --search_type must be specified.\n\n";
    std::cerr << stk::get_options_description() << "\n";
    std::exit(EXIT_FAILURE);
  }

  if (range.mesh_filename.empty()) {
    std::cerr << "\nOPTION ERROR: The '--range_mesh <filename>' option is required for the use cases!\n\n";
    std::cerr << stk::get_options_description() << "\n";
    std::exit(EXIT_FAILURE);
  }

  if (domain.mesh_filename.empty()) {
    domain.mesh_filename = range.mesh_filename;
    domain.mesh_type = range.mesh_type;
  }

  if (strncasecmp("gen:", range.mesh_filename.c_str(), 4) == 0) {
    // Strip off the 'gen:' prefix and set the type to "generated"
    range.mesh_filename = range.mesh_filename.substr(4, range.mesh_filename.size());
    range.mesh_type = "generated";
  }

  if (strncasecmp("gen:", domain.mesh_filename.c_str(), 4) == 0) {
    // Strip off the 'gen:' prefix and set the type to "generated"
    domain.mesh_filename = domain.mesh_filename.substr(4, domain.mesh_filename.size());
    domain.mesh_type = "generated";
  }

  performance_driver(use_case_environment.m_comm,
				  working_directory,
				  search_type, 
				  range,
				  domain,
				  performance);
  
  return 0;
}

/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string>
#include <iostream>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/ProgramOptions.hpp>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

void use_case_1_driver(stk_classic::ParallelMachine  comm,
		       const std::string &working_directory,
		       const std::string &range_mesh,
		       const std::string &range_mesh_type,
		       const std::string &range_entity,
		       const std::string &domain_mesh,
		       const std::string &domain_mesh_type,
		       const std::string &domain_entity);
void use_case_4_driver(stk_classic::ParallelMachine  comm,
		       const std::string &working_directory,
		       const std::string &range_mesh,
		       const std::string &range_mesh_type,
		       const std::string &range_entity,
		       const std::string &domain_mesh,
		       const std::string &domain_mesh_type,
		       const std::string &domain_entity);

namespace bopt = boost::program_options;

int main(int argc, char** argv)
{

  std::string working_directory = "";
  std::string range_mesh        = "";
  std::string range_filetype    = "exodusii";
  std::string range_entity      = "";
  std::string domain_mesh       = "";
  std::string domain_filetype   = "exodusii";
  std::string domain_entity     = "";

  //----------------------------------
  // Process the broadcast command line arguments
  
  bopt::options_description desc("Search use case options");
    
  // NOTE: Options --directory --output-log --runtest are handled/defined in UseCaseEnvironment
  desc.add_options()
    ("range_mesh",    bopt::value<std::string>(&range_mesh),
     "range mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Use 'gears' to generate the gears mesh." )
    ("domain_mesh",   bopt::value<std::string>(&domain_mesh),
     "domain mesh file. Use name of form 'gen:NxMxL' to internally generate a hex mesh of size N by M by L intervals. See GeneratedMesh documentation for more options. Use 'gears' to generate the gears mesh." )
    ("range_entity",    bopt::value<std::string>(&range_entity)->default_value("element"),
     "node, edge, face, element")
    ("domain_entity",   bopt::value<std::string>(&domain_entity)->default_value("element"),
     "node, edge, face, element")
    ("use_case_1",   "search use case 1" )
    ("use_case_4",   "search use case 4" );

  stk_classic::get_options_description().add(desc);

  use_case::UseCaseEnvironment use_case_environment(&argc, &argv);

  bopt::variables_map &vm = stk_classic::get_variables_map();  

  //----------------------------------

  if (range_mesh.empty()) {
    std::cerr << "OPTION ERROR: The '--range_mesh <filename>' option is required for the use cases!\n";
    std::cerr << stk_classic::get_options_description() << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (domain_mesh.empty()) {
    domain_mesh = range_mesh;
    domain_filetype = range_filetype;
  }

  if (strncasecmp("gen:", range_mesh.c_str(), 4) == 0) {
    // Strip off the 'gen:' prefix and set the type to "generated"
    range_mesh = range_mesh.substr(4, range_mesh.size());
    range_filetype = "generated";
  }

  if (strncasecmp("gears:", range_mesh.c_str(), 5) == 0) {
    // Strip off the 'gears:' prefix and set the type to "gears"
    range_mesh = range_mesh.substr(6, range_mesh.size());
    range_filetype = "gears";
  }

  if (strncasecmp("gen:", domain_mesh.c_str(), 4) == 0) {
    // Strip off the 'gen:' prefix and set the type to "generated"
    domain_mesh = domain_mesh.substr(4, domain_mesh.size());
    domain_filetype = "generated";
  }

  if (strncasecmp("gears:", domain_mesh.c_str(), 6) == 0) {
    // Strip off the 'gears:' prefix and set the type to "gears"
    domain_mesh = domain_mesh.substr(6, domain_mesh.size());
    domain_filetype = "gears";
  }

  if (!vm.count("use_case_1") && !vm.count("use_case_4")) {
    std::cout << "OPTION ERROR: At least one of '--use_case_1' or '--use_case_4' must be specified.\n";
    std::exit(EXIT_FAILURE);
  }

  if (vm.count("use_case_1")) {
    use_case_1_driver(use_case_environment.m_comm,
		      working_directory,
		      range_mesh,  range_filetype,  range_entity,
		      domain_mesh, domain_filetype, domain_entity);
  }
  
  if (vm.count("use_case_4")) {
    use_case_4_driver(use_case_environment.m_comm,
		      working_directory,
		      range_mesh,  range_filetype,  range_entity,
		      domain_mesh, domain_filetype, domain_entity);
  }

  // if we've made it this far, the use case has passed
  use_case::print_status(use_case_environment.m_comm, true);

  return 0;
}

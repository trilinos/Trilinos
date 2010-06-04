/*--------------------------------------------------------------------*/
/*    Copyright 2000, 2008, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_CodeTypes.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#ifndef IOSS_STANDALONE
#define OUTPUT sierra::Env::outputP0()
#include <Slib_Env.h>
#else
#define OUTPUT std::cerr
#include <assert.h>
#define ThrowRequire assert
#endif

#include <Ioss_ConcreteVariableType.h>
#include <exodusII/Ioex_SuperElement.h>

using namespace Ioss;

#ifndef IOSS_STANDALONE
namespace {

void bootstrap()
{
  // Add my command line options to the option descriptions.
  boost::program_options::options_description desc("Use case options");
  desc.add_options()
    ("input-deck,i", boost::program_options::value<std::string>(), "Analysis input file")
    ("restart-time,r", boost::program_options::value<std::string>(), "Restart time")
    ("parser-database,p", boost::program_options::value<std::string>(), "Parser database");
  
  stk::get_options_description().add(desc);
}

stk::Bootstrap x(&bootstrap);

} // namespace <unnamed>
#endif

int main(int argc, char *argv[])
{
  StorageInitializer initialize_storage;
#ifndef IOSS_STANDALONE
  sierra::Env::Startup startup__(&argc, &argv, "Utst_superelement-1.1", __DATE__ " " __TIME__); //, opts);

  std::string input_file(sierra::Env::get_param("input-deck"));
#else
  std::string input_file = std::string(argv[1]);
#endif
  if (input_file == "") {
    OUTPUT << "Error: No input file specified\n";
    return (EXIT_FAILURE);
  }

  input_file = Ioss::Utils::local_filename(input_file, "exodusII");
  Ioex::SuperElement se(input_file, "superelement");

  size_t numDOF = se.get_property("numDOF").get_int();
  size_t numEIG = se.get_property("numEIG").get_int();
  size_t numCon = se.get_property("numConstraints").get_int();

  std::cerr << "DOF count = " << numDOF << "\n";
  std::cerr << "EIG count = " << numEIG << "\n";
  std::cerr << "Constraint dof count = " << numCon << "\n";
  ThrowRequire(numCon == numDOF - numEIG);

  // List the fields on the superelement...
  if (!se.field_exists("Kr")) {
    std::cerr << "ERROR: Stiffness matrix field 'Kr' not found\n";
  } else {
    Ioss::Field kr = se.get_field("Kr");
    ThrowRequire(kr.raw_count() == numDOF * numDOF);
  }
  
  if (!se.field_exists("Mr")) {
    std::cerr << "ERROR: Mass matrix field 'Mr' not found\n";
  } else {
    Ioss::Field mr = se.get_field("Mr");
    ThrowRequire(mr.raw_count() == numDOF * numDOF);
  }

  std::vector<double> stiff_mat(numDOF * numDOF);
  std::vector<double> mass_mat(numDOF * numDOF);
  size_t kr_size = se.get_field_data("Kr", stiff_mat);
  size_t mr_size = se.get_field_data("Mr", mass_mat);
  ThrowRequire(kr_size == numDOF * numDOF);
  ThrowRequire(mr_size == numDOF * numDOF);

  OUTPUT << "\nSIERRA execution successful." << std::endl;
  return EXIT_SUCCESS;
}


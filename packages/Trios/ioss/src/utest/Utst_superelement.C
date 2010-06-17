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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#define OUTPUT std::cerr
#undef NDEBUG
#include <assert.h>
#include <Ioss_ConcreteVariableType.h>
#include <exodusII/Ioex_SuperElement.h>

using namespace Ioss;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  StorageInitializer initialize_storage;

  std::string input_file = std::string(argv[argc-1]);
  if (input_file == "") {
    OUTPUT << "Error: No input file specified\n";
    return (EXIT_FAILURE);
  }

  std::string cwd = "";

  input_file = Ioss::Utils::local_filename(input_file, "text", cwd);

  Ioex::SuperElement se(input_file, "superelement");

  size_t numDOF = se.get_property("numDOF").get_int();
  size_t numEIG = se.get_property("numEIG").get_int();
  size_t numCon = se.get_property("numConstraints").get_int();

  std::cerr << "DOF count = " << numDOF << "\n";
  std::cerr << "EIG count = " << numEIG << "\n";
  std::cerr << "Constraint dof count = " << numCon << "\n";
  assert(numCon == numDOF - numEIG);

  // List the fields on the superelement...
  if (!se.field_exists("Kr")) {
    std::cerr << "ERROR: Stiffness matrix field 'Kr' not found\n";
  } else {
    Ioss::Field kr = se.get_field("Kr");
    assert(kr.raw_count() == numDOF * numDOF);
  }
  
  if (!se.field_exists("Mr")) {
    std::cerr << "ERROR: Mass matrix field 'Mr' not found\n";
  } else {
    Ioss::Field mr = se.get_field("Mr");
    assert(mr.raw_count() == numDOF * numDOF);
  }

  std::vector<double> stiff_mat(numDOF * numDOF);
  std::vector<double> mass_mat(numDOF * numDOF);
  size_t kr_size = se.get_field_data("Kr", stiff_mat);
  size_t mr_size = se.get_field_data("Mr", mass_mat);
  assert(kr_size == numDOF * numDOF);
  assert(mr_size == numDOF * numDOF);

  OUTPUT << "\nSIERRA execution successful." << std::endl;
  return EXIT_SUCCESS;
}


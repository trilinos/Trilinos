// Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <iostream>
#include <string>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

#undef NDEBUG
#include "Ioss_ConcreteVariableType.h"
#include "exodus/Ioex_SuperElement.h"
#include <cassert>
#include <stdlib.h>

#include "Ioss_Field.h"
#include "Ioss_Property.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_Utils.h"

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  Ioss::StorageInitializer initialize_storage;

  std::string input_file = std::string(argv[argc - 1]);
  if (input_file.empty()) {
    std::cerr << "Error: No input file specified\n";
    return (EXIT_FAILURE);
  }

  std::string cwd;

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
  }
  else {
    Ioss::Field kr = se.get_field("Kr");
    assert(kr.raw_count() == numDOF * numDOF);
  }

  if (!se.field_exists("Mr")) {
    std::cerr << "ERROR: Mass matrix field 'Mr' not found\n";
  }
  else {
    Ioss::Field mr = se.get_field("Mr");
    assert(mr.raw_count() == numDOF * numDOF);
  }

  std::vector<double> stiff_mat(numDOF * numDOF);
  std::vector<double> mass_mat(numDOF * numDOF);
  size_t              kr_size = se.get_field_data("Kr", stiff_mat);
  size_t              mr_size = se.get_field_data("Mr", mass_mat);
  assert(kr_size == numDOF * numDOF);
  assert(mr_size == numDOF * numDOF);

  std::cerr << "\nSIERRA execution successful." << '\n';
  return EXIT_SUCCESS;
}

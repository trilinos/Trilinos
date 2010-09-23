// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}


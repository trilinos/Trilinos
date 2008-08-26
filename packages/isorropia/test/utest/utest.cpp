//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// ************************************************************************
//@HEADER

#include <Isorropia_ConfigDefs.hpp>
#include <ispatest_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <utst_serial.hpp>

int main(int argc, char** argv) {

  int numProcs = 1;
  int localProc = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

  bool verbose = ispatest::set_verbose(localProc, argc, argv);

  if (localProc == 0) {
    try {
      run_serial_utests(verbose);
    }
    catch(std::exception& exc) {
      std::cout << "utest main caught exception from run_serial_utests: "
            << exc.what() << std::endl;
      return(-1);
    }
  }

  if (verbose) {
    std::cout << "utest main: tests passed."<<std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}


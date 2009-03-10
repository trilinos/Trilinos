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
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraProber.hpp>

#include <Teuchos_CommandLineProcessor.hpp>

#include <ispatest_utils.hpp>
#include <ispatest_epetra_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Vector.h>
#endif

#ifdef HAVE_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#endif

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#ifdef HAVE_EPETRA
Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int local_num_rows = 360;
  int nnz_per_row = local_num_rows/4+1;
  int global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_num_rows/2;

  //adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  //if numProcs is not an even number, undo the local_num_rows adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);

  //create a matrix
  Epetra_CrsMatrix* input_matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_matrix(*input_matrix, nnz_per_row, verbose);
  if (err != 0) {
    delete input_matrix;
    return(0);
  }

  return(input_matrix);
}



bool probing_test(Epetra_CrsMatrix & in_mat, bool build_list){
  const Epetra_CrsGraph & graph=in_mat.Graph();

  Teuchos::ParameterList main, zoltan;

  if(build_list){
    zoltan.set("DISTANCE","2");
    main.set("ZOLTAN",zoltan);
  }

  Isorropia::Epetra::Prober prober(Teuchos::rcp<const Epetra_CrsGraph>(&graph,false),main);
  Epetra_CrsMatrix out_mat(Copy,graph);
  int rv=prober.probe(in_mat,out_mat);
  if(rv!=0) {printf("ERROR: probing failed\n");return false;}
  
#ifdef HAVE_EPETRAEXT
  EpetraExt::MatrixMatrix::Add(in_mat,false,1,out_mat,-1);
  double nrm=out_mat.NormInf()/in_mat.NormInf();
  if(!in_mat.Comm().MyPID())
    printf("diff norm = %22.16e\n",nrm);
  if(nrm < 1e-12) return true;
  else return false;
#endif
  return true;
}










#endif

int main(int argc, char** argv) {


  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif
  
  try {
    verbose = ispatest::set_verbose(localProc, argc, argv);
  }
  catch(std::exception& exc) {
    std::cout << "err, setting verbosity: " << exc.what() << std::endl;
    std::cout << "End Result: TEST FAILED" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif /* HAVE_MPI */
    return(-1);
  }

  // Only print on proc 0
  if (localProc != 0) {
    verbose = false;
  }

  bool test_passed = true;
#ifdef HAVE_EPETRA
  //TEST HERE  
  Epetra_CrsMatrix* mat=create_epetra_test_matrix_1(numProcs,localProc,true);
  test_passed=probing_test(*mat,true);
  test_passed=test_passed && probing_test(*mat,false);

  // Cleanup
  delete mat;
  
#else
  std::cout << "probing main: currently can only test "
         << "rebalancing with Epetra enabled." << std::endl;
  test_passed = false;
#endif

  if (test_passed && verbose) {
    std::cout << "PROBING main: tests passed."<<std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}



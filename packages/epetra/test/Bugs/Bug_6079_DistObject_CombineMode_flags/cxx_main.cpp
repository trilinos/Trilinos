// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check correct behavior of sadd of Trilinos vectors
// if they have different Epetra maps 

#include <iostream>
#include <vector>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_FEVector.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif

void test ()
{
  int n_proc; 
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
#else
  n_proc = 1;
#endif
  int my_id;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
#else
  my_id = 0;
#endif

  //All processes should own 10 entries
  const int entries_per_process = 10;

  const int begin_index = my_id*entries_per_process;
  const int end_index = (my_id+1)*entries_per_process;

  const int local_begin = std::max(0, begin_index-entries_per_process/2);
  const int local_end = entries_per_process*n_proc;

  //create Epetra maps
  std::vector<unsigned int> ghosted_indices;
  ghosted_indices.reserve(local_end-local_begin);
  for (int i = local_begin; i< local_end; ++i)
    ghosted_indices.push_back(i);
  Epetra_Map map_ghosted
    (-1,
    local_end-local_begin,
    reinterpret_cast<int*>(&ghosted_indices[0]),
    0,
#ifdef HAVE_MPI
    Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Epetra_SerialComm());
#endif
  
  std::vector<unsigned int> distributed_indices;
  distributed_indices.reserve(entries_per_process*n_proc);
  for (int i = begin_index; i< end_index; ++i)
    distributed_indices.push_back(i);
  Epetra_Map map_distributed
  (entries_per_process*n_proc,
   entries_per_process,
   reinterpret_cast<int*>(&distributed_indices[0]),
   0,
#ifdef HAVE_MPI
   Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Epetra_SerialComm());
#endif
 
  Epetra_FEVector v_ghosted(map_ghosted);
  Epetra_FEVector v_distributed(map_distributed);
  
  v_distributed.PutScalar(2.);
  v_ghosted.PutScalar(1.);

  Epetra_Import data_exchange (v_distributed.Map(), v_ghosted.Map());
  int ierr = v_distributed.Import(v_ghosted, data_exchange, Epetra_AddLocalAlso);
 
  std::cout << "Distributed:" << std::endl;
  for (int i=begin_index; i<end_index; ++i)
  {
    int trilinos_i
      = v_distributed.Map().LID(i);
    double value = v_distributed[0][trilinos_i];
    std::cout<<"proc "<<my_id<<" "<< i << ": " << value << std::endl;
	if(value != 3)
		std::cerr << "tests FAILED: value = " << value << std::endl;
  }  
}



int main (int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  test();
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}



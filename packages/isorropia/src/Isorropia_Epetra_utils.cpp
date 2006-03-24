//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_Exception.hpp>
#include <Isorropia_Utils.hpp>
#include <Isorropia_Epetra_utils.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>
#endif

#if defined(HAVE_EPETRA) && defined(HAVE_MPI)
#include <Epetra_MpiComm.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {
namespace Epetra_Utils {

#ifdef HAVE_EPETRA
Epetra_Vector* create_row_weights_nnz(const Epetra_RowMatrix& input_matrix)
{
  const Epetra_BlockMap& input_rowmap = input_matrix.RowMatrixRowMap();
  Epetra_Vector* weights = new Epetra_Vector(input_rowmap);
  double* weights_ptr = 0;
  weights->ExtractView(&weights_ptr);
  int local_num_rows = input_rowmap.NumMyElements();

  for(int i=0; i<local_num_rows; ++i) {
    int nnz;
    int err = input_matrix.NumMyRowEntries(i, nnz);
    if (err != 0) {
      throw Isorropia::Exception("create_row_weights_nnz: error in input_matrix.NumMyRowEntries");
    }

    weights_ptr[i] = 1.0*nnz;
  }

  return( weights );
}

Epetra_Vector* create_row_weights_nnz(const Epetra_CrsGraph& input_graph)
{
  const Epetra_BlockMap& input_rowmap = input_graph.RowMap();
  Epetra_Vector* weights = new Epetra_Vector(input_rowmap);
  double* weights_ptr = 0;
  weights->ExtractView(&weights_ptr);
  int local_num_rows = input_rowmap.NumMyElements();

  for(int i=0; i<local_num_rows; ++i) {
    int nnz = input_graph.NumMyIndices(i);

    weights_ptr[i] = 1.0*nnz;
  }

  return( weights );
}

Epetra_Map create_rowmap_balanced(const Epetra_BlockMap& input_rowmap,
                                  const Epetra_Vector& weights)
{
  //create a dummy map which will be reassigned to the newly-created
  //balanced rowmap later...
  const Epetra_Comm& input_comm = input_rowmap.Comm();
  Epetra_Map rowmap(10, 0, input_comm);

  //next we're going to collect weights onto proc 0.
  int global_num_rows = input_rowmap.NumGlobalElements();
  int myPID = input_comm.MyPID();
  int numProcs = input_comm.NumProc();
  int local_num_rows = myPID == 0 ? global_num_rows : 0;
  Epetra_BlockMap proc0_rowmap(global_num_rows, local_num_rows, 1,0,input_comm);
  Epetra_Vector proc0_weights(proc0_rowmap);

  Epetra_Import importer(proc0_rowmap, input_rowmap);
  proc0_weights.Import(weights, importer, Insert);

#ifdef HAVE_MPI
  int tag = 1212121;
  const Epetra_MpiComm* mpiComm =
    dynamic_cast<const Epetra_MpiComm*>(&input_comm);
  if (mpiComm == 0) {
    throw Isorropia::Exception("dynamic_cast to MpiComm failed.");
  }
  MPI_Comm mpicomm = mpiComm->GetMpiComm();
#endif

  int p;
  double total_weight;
  weights.Norm1(&total_weight);

  std::vector<int> all_proc_old_offsets;
  Isorropia::Epetra_Utils::gather_all_proc_global_offsets(input_rowmap,
                                                     all_proc_old_offsets);
  std::vector<int> all_proc_new_offsets(numProcs+1);

  if (myPID == 0) {
    double avg_weight = total_weight/numProcs;

    double* proc0_weights_ptr;
    proc0_weights.ExtractView(&proc0_weights_ptr);
    int weights_length = proc0_weights.MyLength();

    int offset = 0;
    for(p=0; p<numProcs; ++p) {
      all_proc_new_offsets[p] = offset;

      double tmp_weight = 0.0;
      while(offset < weights_length && tmp_weight < avg_weight) {
        tmp_weight += proc0_weights_ptr[offset++];
      }
    }
    all_proc_new_offsets[numProcs] = offset;

#ifdef HAVE_MPI
    if (numProcs > 1) {
      //now broadcast the new offsets
      input_comm.Broadcast(&all_proc_new_offsets[0], numProcs+1, 0);
    }
#endif
  }
  else { //myPID != 0
#ifdef HAVE_MPI
    input_comm.Broadcast(&all_proc_new_offsets[0], numProcs+1, 0);
#endif
  }

  //Now all processors have two lists, all_proc_old_offsets and
  //all_proc_new_offsets.
  //The map elements that this processor currently holds are
  //given by the range:
  //      [all_proc_old_offsets[myPID]...all_proc_old_offsets[myPID+1]-1]
  //and the map elements that this processor needs to hold are
  //      [all_proc_new_offsets[myPID]...all_proc_new_offsets[myPID+1]-1]
  //
  //So now we need to figure out which elements we need to send/recv
  //to/from neighboring processors.

  std::vector<int> send_info;
  std::vector<int> recv_info;
  Isorropia::Utils::create_comm_plan(myPID, all_proc_old_offsets,
                                     all_proc_new_offsets,
                                     send_info, recv_info);

  int new_num_local = all_proc_new_offsets[myPID+1]-all_proc_new_offsets[myPID];

  const int* old_gids = input_rowmap.MyGlobalElements();

  std::vector<int> new_gids(new_num_local);
#ifdef HAVE_MPI
  unsigned i=0;
  MPI_Request* reqs = recv_info.size() > 0 ?
     new MPI_Request[recv_info.size()/3] : 0;
  MPI_Status* sts = recv_info.size() > 0 ?
     new MPI_Status[recv_info.size()/3] : 0;
  while(i<recv_info.size()) {
    int proc = recv_info[i];
    int recv_offset = recv_info[i+1];
    int num_recv = recv_info[i+2];
    
    MPI_Irecv(&new_gids[recv_offset], num_recv, MPI_INT, proc,
             tag, mpicomm, &reqs[i/3]);
    i += 3;
  }

  i=0;
  while(i<send_info.size()) {
    int proc = send_info[i];
    int send_offset = send_info[i+1];
    int num_send = send_info[i+2];

    MPI_Send((void*)&old_gids[send_offset], num_send, MPI_INT,
             proc, tag, mpicomm);

    i += 3;
  }
#endif

  //now copy any overlapping elements from old_gids into new_gids.
  int old_start = all_proc_old_offsets[myPID];
  int new_start = all_proc_new_offsets[myPID];
  int old_end = all_proc_old_offsets[myPID+1]-1;
  int new_end = all_proc_new_offsets[myPID+1]-1;

  int overlap_start = new_start > old_start ? new_start : old_start;
  int overlap_end = new_end < old_end ? new_end : old_end;

  int copy_src_offset = overlap_start - old_start;
  int copy_dest_offset = overlap_start - new_start;

  int num_copy = overlap_end - overlap_start + 1;

  for(int j=0; j<num_copy; ++j) {
    new_gids[copy_dest_offset++] = old_gids[copy_src_offset++];
  }

#ifdef HAVE_MPI
  //now make sure the recvs are finished...
  if (recv_info.size() > 0) {
    MPI_Waitall(recv_info.size()/3, reqs, sts);
  }
#endif

  Epetra_Map tmp_map(global_num_rows, new_num_local, &new_gids[0],
                     0, input_comm);
  rowmap = tmp_map;

  return(rowmap);
}

void gather_all_proc_global_offsets(const Epetra_BlockMap& blkmap,
                                    std::vector<int>& all_proc_offsets)
{
  const Epetra_Comm& comm = blkmap.Comm();
  int numProcs = comm.NumProc();
  int myPID = comm.MyPID();

  all_proc_offsets.resize(numProcs+1);
  for(int i=0; i<numProcs+1; ++i) {
    all_proc_offsets[0] = 0;
  }

  //first put num-local-elements in position myPID, and gather-all so
  //that each proc has all entries.
  all_proc_offsets[myPID] = blkmap.NumMyElements();
  comm.GatherAll(&all_proc_offsets[myPID], &all_proc_offsets[0], 1);

  //now run the list and turn the local-sizes into global offsets.
  int offset = 0;
  for(int p=0; p<numProcs; ++p) {
    int tmp = all_proc_offsets[p];
    all_proc_offsets[p] = offset;
    offset += tmp;
  }
  all_proc_offsets[numProcs] = offset;
}

void import_matrix(const Epetra_CrsMatrix& input_matrix,
                   Epetra_CrsMatrix& target_matrix)
{
  Epetra_Import importer(target_matrix.RowMap(), input_matrix.RowMap());
  target_matrix.Import(input_matrix, importer, Insert);
}

void import_graph(const Epetra_CrsGraph& input_graph,
                   Epetra_CrsGraph& target_graph)
{
  Epetra_Import importer(target_graph.RowMap(), input_graph.RowMap());
  target_graph.Import(input_graph, importer, Insert);
}

#endif //HAVE_EPETRA

}//namespace Epetra_Utils
}//namespace Isorropia


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

************************************************************************
*/
//@HEADER

#include <Isorropia_Exception.hpp>
#include <Isorropia_Utils.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

#if defined(HAVE_EPETRA) && defined(HAVE_MPI)
#include <Epetra_MpiComm.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace Isorropia {

namespace Epetra {

#ifdef HAVE_EPETRA
Teuchos::RCP<Partitioner>
create_partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_graph, paramlist));
  return(partitioner);
}

Teuchos::RCP<Partitioner>
create_partitioner(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		   Teuchos::RCP<CostDescriber> costs,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_graph, costs, paramlist));
  return(partitioner);
}

Teuchos::RCP<Partitioner>
create_partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_matrix, paramlist));
  return(partitioner);
}

Teuchos::RCP<Partitioner>
create_partitioner(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
		   Teuchos::RCP<CostDescriber> costs,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_matrix, costs, paramlist));
  return(partitioner);
}


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

double compute_imbalance(int nprocs, std::vector<int> &offsets, double *wgts, double target)
{
  double imbalance = 1.0;

  for (int p=0; p < nprocs; p++){

    double pweight = 0.0;

    for (int row=offsets[p]; row < offsets[p+1]; row++){
      pweight += wgts[row];
    }

    double ib = 1.0;

    if (pweight <= target)
      ib += ((target - pweight) / target);
    else
      ib += ((pweight - target) / target);

    if (ib > imbalance) imbalance = ib;
  }

  return imbalance;
}

int
repartition(const Epetra_BlockMap& input_map,
	    const Epetra_Vector& weights,
	    std::vector<int>& myNewElements,
	    int& exportsSize,
	    std::vector<int>& imports)
//             std::map<int,int>& exports,
//             std::map<int,int>& imports)
{
  if (!input_map.PointSameAs(weights.Map())) {
    std::string str1("Epetra::repartition ERROR, input_map not ");
    std::string str2("equivalent size/layout to weights.Map()");
    throw Isorropia::Exception(str1+str2);
  }

  const Epetra_Comm& input_comm = input_map.Comm();

  imports.clear();

  //first we're going to collect weights onto proc 0.
  int myPID = input_comm.MyPID();
  int numProcs = input_comm.NumProc();
  int global_num_rows = input_map.NumGlobalElements();
  int local_num_rows = myPID == 0 ? global_num_rows : 0;
  Epetra_BlockMap proc0_rowmap(global_num_rows, local_num_rows, 1,0,input_comm);
  Epetra_Vector proc0_weights(proc0_rowmap);

  Epetra_Import importer(proc0_rowmap, input_map);
  proc0_weights.Import(weights, importer, Insert);

  double total_weight;
  weights.Norm1(&total_weight);

  std::vector<int> all_proc_old_offsets;
  gather_all_proc_global_offsets(input_map, all_proc_old_offsets);
  std::vector<int> all_proc_new_offsets(numProcs+1);

  if (myPID == 0) {
    double avg_weight = total_weight/numProcs;

    double* proc0_weights_ptr;
    proc0_weights.ExtractView(&proc0_weights_ptr);
    int weights_length = proc0_weights.MyLength();

    double old_imbalance = 
      compute_imbalance(numProcs, all_proc_old_offsets, proc0_weights_ptr, avg_weight);

    int offset = 0;
    for(int p=0; p<numProcs; ++p) {
      all_proc_new_offsets[p] = offset;

      double tmp_weight = 0.0;

      while(offset < weights_length && tmp_weight < avg_weight) {
        tmp_weight += proc0_weights_ptr[offset++];
      }
    }
    all_proc_new_offsets[numProcs] = weights_length;

    double new_imbalance = 
      compute_imbalance(numProcs, all_proc_new_offsets, proc0_weights_ptr, avg_weight);

    // Because this is a quick and dirty partitioning, it is possible that
    // if the balance was good to begin with, that we have just created a
    // slightly worse balance.  In that case, return to the old partitioning.

    if (new_imbalance > old_imbalance){
      for (int proc=0; proc <= numProcs; proc++){
        all_proc_new_offsets[proc] = all_proc_old_offsets[proc];
      }
    }

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

  //Now we need to figure out which elements we need to send/recv
  //to/from neighboring processors.

  std::vector<int> send_info;
  std::vector<int> recv_info;
  Isorropia::Utils::create_comm_plan(myPID, all_proc_old_offsets,
                                     all_proc_new_offsets,
                                     send_info, recv_info);

  //Create the list to hold local elements for the new map.
  myNewElements.assign(all_proc_old_offsets[myPID+1]-all_proc_old_offsets[myPID], myPID);
//   myNewElements.resize(new_num_local);

  const int* old_gids = input_map.MyGlobalElements();

#ifdef HAVE_MPI

  int tag = 1212121;
  const Epetra_MpiComm* mpiComm =
    dynamic_cast<const Epetra_MpiComm*>(&input_comm);
  if (mpiComm == 0) {
    throw Isorropia::Exception("dynamic_cast to MpiComm failed.");
  }
  int num_reqs = recv_info.size()/3;

  MPI_Comm mpicomm = mpiComm->GetMpiComm();
  MPI_Request* reqs = num_reqs > 0 ? new MPI_Request[num_reqs] : 0;

  unsigned i=0;
  unsigned int base_offset = (num_reqs > 0)?(recv_info[1]):0;

  int num_import = 0;
  while(i<recv_info.size()) {
    num_import += recv_info[i+2];
    i+=3;
  }
  imports.resize(num_import);

  for (i = 0 ;  i<recv_info.size(); i+=3) {
    int proc = recv_info[i];
    int recv_offset = recv_info[i+1];
    int num_recv = recv_info[i+2];

    MPI_Irecv(&imports[recv_offset-base_offset], num_recv, MPI_INT, proc,
             tag, mpicomm, &reqs[i/3]);
  }

  exportsSize = 0;
  for (i = 0 ;  i<send_info.size(); i+=3) {
    int proc = send_info[i];
    int send_offset = send_info[i+1];
    int num_send = send_info[i+2];

    MPI_Send((void*)&old_gids[send_offset], num_send, MPI_INT,
             proc, tag, mpicomm);

    exportsSize += num_send;
    for(int j=0; j<num_send; ++j) {
      myNewElements[send_offset+j] = proc;
    }
  }
  //make sure the recvs are finished...
  if (recv_info.size() > 0) {
    MPI_Waitall(recv_info.size()/3, reqs, MPI_STATUSES_IGNORE);
  }
  delete [] reqs;

#else /* HAVE_MPI */

  imports.clear();
  unsigned int i=0;
  while(i<send_info.size()) {
    int proc = send_info[i];
    int send_offset = send_info[i+1];
    int num_send = send_info[i+2];

    for(int j=0; j<num_send; ++j) {
      myNewElements[send_offset+j] = proc;
    }

    i += 3;
  }

  exportsSize = 0;
  for (i = 0 ;  i<send_info.size(); i+=3)
    exportsSize += send_info[i+2];

  int num_import = 0;
  while(i<recv_info.size()) {
    num_import += recv_info[i+2];
    i+=3;
  }
  imports.resize(num_import);

  int num_reqs = recv_info.size()/3;
  unsigned int base_offset = (num_reqs > 0)?(recv_info[1]):0;
  for (i = 0 ;  i<recv_info.size(); i+=3) {
    int proc = recv_info[i];
    int recv_offset = recv_info[i+1];
    int num_recv = recv_info[i+2];

    for (int j = 0 ;  i<send_info.size(); i+=3) {
      if (recv_info[i] != myPID)
	continue;
      memcpy(&imports[recv_offset-base_offset], &old_gids[send_info[j+1]],
	     num_recv*sizeof(int));
    }
  }

#endif /* HAVE_MPI */

  return(0);
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
  int tmpOffset = all_proc_offsets[myPID];
  comm.GatherAll(&tmpOffset, &all_proc_offsets[0], 1);

  //now run the list and turn the local-sizes into global offsets.
  int offset = 0;
  for(int p=0; p<numProcs; ++p) {
    int tmp = all_proc_offsets[p];
    all_proc_offsets[p] = offset;
    offset += tmp;
  }
  all_proc_offsets[numProcs] = offset;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix)
{
  CostDescriber costs;
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  vwgts.release();
  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_RowMatrix> matrixPtr=
    Teuchos::rcp(&(input_matrix), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(matrixPtr, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    rd.redistribute(input_matrix);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  // release of ownership of vwgts, so destruction of stack objects
  // does not do an extra destruction of row_weights and cause
  // a fatal fault
  vwgts.release();  

  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> input_graph =
    Teuchos::rcp(&(input_matrix.Graph()), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_graph, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    rd.redistribute(input_matrix);

  return balanced_matrix;
}
Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    create_balanced_copy(input_graph, costs, paramlist);

  return balanced_graph;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  // release of ownership of vwgts, so destruction of stack objects
  // does not do an extra destruction of row_weights and cause
  // a fatal fault
  vwgts.release();  

  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    create_balanced_copy(input_graph, costs, paramlist);

  return balanced_graph;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    create_balanced_copy(input_graph, costs, paramlist);

  return balanced_graph;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> graphPtr=
    Teuchos::rcp(&(input_graph), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(graphPtr, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    rd.redistribute(input_graph);

  return balanced_graph;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    create_balanced_copy(input_problem, costs, paramlist);

  return linprob;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  // release of ownership of vwgts, so destruction of stack objects
  // does not do an extra destruction of row_weights and cause
  // a fatal fault
  vwgts.release();  
  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    create_balanced_copy(input_problem, costs, paramlist);

  return linprob;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    create_balanced_copy(input_problem, costs, paramlist);

  return linprob;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_RowMatrix> rowmat =
    Teuchos::rcp(input_problem.GetMatrix(), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(rowmat, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    rd.redistribute(*input_problem.GetMatrix());

  Teuchos::RCP<Epetra_MultiVector> balanced_rhs =
    rd.redistribute(*input_problem.GetRHS());

  Teuchos::RCP<Epetra_MultiVector> x=
    Teuchos::rcp(new Epetra_MultiVector(*input_problem.GetLHS()));

  // prevent these from being deallocated on return from create_balanced_copy
  balanced_matrix.release(); 
  balanced_rhs.release();
  x.release();

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    Teuchos::rcp(new Epetra_LinearProblem(balanced_matrix.get(), x.get(), balanced_rhs.get()));

  return( linprob );
}
#endif
}//namespace Epetra
}//namespace Isorropia

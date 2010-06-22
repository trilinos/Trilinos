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
#include <Isorropia_Tpetra.hpp>
#include <Teuchos_Comm.hpp>
#include <Isorropia_TpetraPartitioner.hpp>
#include <Isorropia_TpetraRedistributor.hpp>

#ifdef HAVE_ISORROPIA_TPETRA

#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace Isorropia {

namespace Tpetra {

#ifdef HAVE_ISORROPIA_TPETRA

template <typename Node>
::Tpetra::MultiVector<double,int,int,Node>* create_row_weights_nnz(const ::Tpetra::RowMatrix<double,int,int,Node>& input_matrix)
{
  int stride;
  const ::Tpetra::Map<int,int,Node>& input_rowmap = input_matrix.RowMatrixRowMap();
  ::Tpetra::MultiVector<double,int,int,Node>* weights = new ::Tpetra::MultiVector<double,int,int,Node>(input_rowmap, 1);
  double* weights_ptr = 0;
  weights->ExtractView(&weights_ptr, &stride);
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

template <typename Node>
::Tpetra::MultiVector<double,int,int,Node>* create_row_weights_nnz(const ::Tpetra::CrsGraph<int,int,Node>& input_graph)
{
  int stride;
  const ::Tpetra::Map<int,int,Node>& input_rowmap = input_graph.RowMap();
  ::Tpetra::MultiVector<double,int,int,Node>* weights = new ::Tpetra::MultiVector<double,int,int,Node>(input_rowmap, 1);
  double* weights_ptr = 0;
  weights->ExtractView(&weights_ptr, &stride);
  int local_num_rows = input_rowmap.NumMyElements();

  for(int i=0; i<local_num_rows; ++i) {
    int nnz = input_graph.NumMyIndices(i);

    weights_ptr[i] = 1.0*nnz;
  }

  return( weights );
}

template <typename Node>
::Tpetra::MultiVector<double,int,int,Node>* create_unit_weights(const ::Tpetra::MultiVector<double,int,int,Node>& input_coords)
{
  int stride;
  double* weights_ptr = 0;

  ::Tpetra::MultiVector<double,int,int,Node>* weights = new ::Tpetra::MultiVector<double,int,int,Node>(input_coords.Map(), 1);
  weights->ExtractView(&weights_ptr, &stride);

  for(int i=0; i<input_coords.MyLength(); ++i) {
    weights_ptr[i] = 1.0;
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

template <typename Node>
int repartition(const ::Tpetra::Map<int,int,Node>& input_map,
	       const ::Tpetra::MultiVector<double,int,int,Node>& weights,
	       std::vector<int>& newPartition,    // new partition for each of my objects
	       int& exportsSize,                  // how many of my objects are exports
	       std::vector<int>& imports)         // list of gids I will import
{
  if (!input_map.PointSameAs(weights.Map())) {
    std::string str1("Tpetra::repartition ERROR, input_map not ");
    std::string str2("equivalent size/layout to weights.Map()");
    throw Isorropia::Exception(str1+str2);
  }
  int stride=0;
  double *vals = NULL;

  const Teuchos::Comm<int> & comm = input_map.Comm();
  int base = input_map.IndexBase();
  int myPID = comm.getRank();
  int numProcs = comm.getSize();

  int globalSize = input_map.NumGlobalElements();
  int localSize = input_map.NumMyElements();
  int mySize = ((myPID==0) ? globalSize : 0);
  int myCountSize = ((myPID==0) ? numProcs: 0);

  // Create some maps.  The weight vector is not necessarily contiguous by
  // process rank, so we need to reorganize it.

  // Distributed map of objects that are contiguous by process rank, and those
  //  objects gathered onto process zero

  ::Tpetra::Map<int,int,Node> distMap(globalSize, localSize, 1, base, comm);
  ::Tpetra::Map<int,int,Node> procZeroMap(globalSize, mySize, 1, base, comm);

  // Map of counts for each process, and counts gathered on process 0

  ::Tpetra::Map<int,int,Node> distCountMap(numProcs, 1, 1, base, comm);
  ::Tpetra::Map<int,int,Node> procZeroCountMap(numProcs, myCountSize, 1, base, comm);

  // Some importers

  ::Tpetra::Import<int,int,Node> DistToProcZero(procZeroMap, distMap);
  ::Tpetra::Import<int,int,Node> ProcZeroToDist(distMap, procZeroMap);
  ::Tpetra::Import<int,int,Node> CountsToProcZero(procZeroCountMap, distCountMap);
  ::Tpetra::Import<int,int,Node> CountsToProcs(distCountMap, procZeroCountMap);

  // Move all weights to process 0

  //weights.ExtractView(&vals, &stride);
  Teuchos::ArrayRCP<double> valArray = weights.get1dView();

  ::Tpetra::Vector<double,int,int,Node> distWeights(distMap, valArray());
  ::Tpetra::Vector<double,int,int,Node> weightVal(procZeroMap);
  weightVal.doImport(distWeights, DistToProcZero, ::Tpetra::INSERT);

  // Move the count of weights for each process to process 0

  ::Tpetra::Vector<double,int,int,Node> distCount(distCountMap);
  distCount[0] = localSize;
  ::Tpetra::Vector<double,int,int,Node> weightCount(procZeroCountMap);
  weightCount.doImport(distCount, CountsToProcZero, ::Tpetra::INSERT);

  // Move all global IDs to process 0

  double *myGIDs = new double [localSize];
  for (int i=0; i<localSize; i++){
    myGIDs[i] = input_map.GID(i);
  }

  

  Teuchos::ArrayView<double>::ArrayView   myGIDView(myGIDs, localSize);


  ::Tpetra::Vector<double,int,int,Node> localGIDs(distMap, myGIDView);


  ::Tpetra::Vector<double,int,int,Node> gids(procZeroMap);
  gids.Import(localGIDs, DistToProcZero, ::Tpetra::INSERT);
  delete [] myGIDs;

  // Process zero computes a new balance

  double *numImports = new double [numProcs];
  int totalImports=0;
  double *importGIDs = NULL;
  double *ids = NULL;

  if (myPID == 0){

    double total_weight=0.0; 
    weightVal.ExtractView(&vals);
    for (int i=0; i<globalSize; i++){
      total_weight += vals[i];
    }

    double avg_weight = total_weight/numProcs;

    std::vector<int> all_proc_old_offsets(numProcs+1, globalSize);
    std::vector<int> all_proc_new_offsets(numProcs+1, globalSize);

    for (int i=numProcs-1; i >= 0 ; i--){
      all_proc_old_offsets[i] = all_proc_old_offsets[i+1] - static_cast<int>(weightCount[i]);
    }

    double old_imbalance =
      compute_imbalance(numProcs, all_proc_old_offsets, vals, avg_weight);

    int offset = 0;
    for(int p=0; p<numProcs; ++p) {
      all_proc_new_offsets[p] = offset;

      double tmp_weight = 0.0;

      while(offset < globalSize && tmp_weight < avg_weight) {
        tmp_weight += vals[offset++];
      }
    }

    double new_imbalance =
      compute_imbalance(numProcs, all_proc_new_offsets, vals, avg_weight);

    // Because this is a quick and dirty partitioning, it is possible that
    // if the balance was good to begin with, that we have just created a
    // slightly worse balance.  In that case, return to the old partitioning.

    if (new_imbalance >= old_imbalance){
      for (int proc=0; proc <= numProcs; proc++){
        all_proc_new_offsets[proc] = all_proc_old_offsets[proc];
      }
    }

    // Fill weight vector with the new process ID (partition ID) for each object

    for (int proc=0; proc < numProcs; proc++){
      int len = all_proc_new_offsets[proc+1] - all_proc_new_offsets[proc];
      for (int j=0; j<len; j++){
        *vals++ = static_cast<double>(proc);
      }
    }

    // Count the number of imports for each process, and create a list of
    // gids to be imported to each process

    totalImports = 0;
    importGIDs = new double [globalSize];
    double *curr = importGIDs;
    gids.ExtractView(&ids);

    for (int proc=0; proc < numProcs; proc++){
      numImports[proc] = 0;
      int from = all_proc_new_offsets[proc];
      int to =  all_proc_new_offsets[proc+1];
      int old_from = all_proc_old_offsets[proc];
      int old_to =  all_proc_old_offsets[proc+1];

      for (int i=from; i<to ; i++){
        if ((i< old_from) || (i >= old_to)){
          *curr++ = ids[i];
          numImports[proc]++;
        }
      }
      totalImports += static_cast<int>(numImports[proc]);
    }
  }

  // The new partition IDs are in the weightVal vector on process zero

  distWeights.doImport(weightVal, ProcZeroToDist, ::Tpetra::INSERT);

  newPartition.resize(localSize);
  for (int i=0; i<localSize; i++){
    newPartition[i] = static_cast<int>(distWeights[i]);
  }

  // Get the count of my imports - reuse the weightCount vector

  int *indices = new int [numProcs];
  for (int i=0; i<numProcs; i++){
    indices[i] = i + base;
  }

  weightCount.ReplaceGlobalValues(((myPID==0) ? numProcs : 0), numImports, indices);
  distCount.doImport(weightCount, CountsToProcs, ::Tpetra::INSERT);
  delete [] indices;
  if (numImports) delete [] numImports;

  // Get the list of my imports

  int numMyImports = static_cast<int>(distCount[0]);
  double n1;
  distCount.Norm1(&n1);
  totalImports = static_cast<int>(n1);
  int myImportSize = ((myPID==0) ? totalImports : 0);

  ::Tpetra::Map<int,int,Node> distImportMap(totalImports, numMyImports, 1, base, comm);
  ::Tpetra::Map<int,int,Node> proc0ImportMap(totalImports, myImportSize, 1, base, comm);
  ::Tpetra::Import<int,int,Node> importer(distImportMap, proc0ImportMap);

  ::Tpetra::Vector<double,int,int,Node> proc0imports(proc0ImportMap, importGIDs);
  ::Tpetra::Vector<double,int,int,Node> distImports(distImportMap);

  distImports.doImport(proc0imports, importer, ::Tpetra::INSERT);

  if (importGIDs) delete [] importGIDs;

  imports.resize(numMyImports);
  for (int i=0; i<numMyImports; i++){
    imports[i] = static_cast<int>(distImports[i]);
  }

  // Finally, the number of exports

  exportsSize = 0;
  for (int i=0; i<localSize; i++){
    if (newPartition[i] != myPID){
      exportsSize++;
    }
  }

  return 0;
}

/* New createBalancedCopy functions */

template <typename Node>
::Tpetra::MultiVector<double,int,int,Node> * 
createBalancedCopy(const ::Tpetra::MultiVector<double,int,int,Node> &coords)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(coords, paramlist);
}

template <typename Node>
::Tpetra::MultiVector<double,int,int,Node> * 
createBalancedCopy(const ::Tpetra::MultiVector<double,int,int,Node> &coords,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > coordRcp = Teuchos::rcp(&coords, false);

  Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> > partitioner 
    = Teuchos::rcp(new Isorropia::Tpetra::Partitioner<Node>(coordRcp, paramlist));

  Isorropia::Tpetra::Redistributor<Node> rd(partitioner);

  Teuchos::RCP< ::Tpetra::MultiVector<double,int,int,Node> > newVec = rd.redistribute(coords);

  newVec.release();

  return newVec.get();
}

template <typename Node>
::Tpetra::CrsGraph<int,int,Node> *
createBalancedCopy(const ::Tpetra::CrsGraph<int,int,Node>& input_graph)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(input_graph, paramlist);
}

template <typename Node>
::Tpetra::CrsGraph<int,int,Node> *
createBalancedCopy(const ::Tpetra::CrsGraph<int,int,Node>& input_graph,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > rcp_input_graph =
    Teuchos::rcp(&(input_graph), false);

  Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> > partitioner =
    Teuchos::rcp(new Isorropia::Tpetra::Partitioner<Node> (rcp_input_graph, paramlist));

  Isorropia::Tpetra::Redistributor<Node> rd(partitioner);

  Teuchos::RCP< ::Tpetra::CrsGraph<int,int,Node> > balanced_graph =
    rd.redistribute(input_graph);

  balanced_graph.release();

  return balanced_graph.get();
}

template <typename Node>
::Tpetra::CrsMatrix<double,int,int,Node> *
createBalancedCopy(const ::Tpetra::CrsMatrix<double,int,int,Node>& input_matrix)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(input_matrix, paramlist);
}

template <typename Node>
::Tpetra::CrsMatrix<double,int,int,Node> *
createBalancedCopy(const ::Tpetra::CrsMatrix<double,int,int,Node>& input_matrix,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph =
    Teuchos::rcp(&(input_matrix.Graph()), false);

  Teuchos::RCP<Isorropia::Tpetra::Partitioner<Node> > partitioner =
    Teuchos::rcp(new Isorropia::Tpetra::Partitioner<Node>(input_graph, paramlist));

  Isorropia::Tpetra::Redistributor<Node> rd(partitioner);

  Teuchos::RCP< ::Tpetra::CrsMatrix<double,int,int,Node> > balanced_matrix =
    rd.redistribute(input_matrix);

  balanced_matrix.release();
  return balanced_matrix.get();
}


#endif
}//namespace Tpetra
}//namespace Isorropia

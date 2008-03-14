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

#include <ispatest_lbeval_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_BlockMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#endif

#include <map>
#include <set>
#include <iostream>

/** ispatest is the namespace that contains isorropia's test-utilities
*/
namespace ispatest {

#ifdef HAVE_EPETRA
static double compute_balance(const Epetra_Comm &comm, int myRows, int nwgts, float *wgts);
int compute_graph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl)
{
  const Epetra_Comm &comm  = graph.Comm();
  int myProc = comm.MyPID();
  int myRows = graph.NumMyRows();
  int rc;
  int *vgid = NULL;
  float *vwgt = NULL;

  // Compute the balance

  int numVWgts = costs.getNumVertices();

  if ((numVWgts > 0) && (numVWgts != myRows)){
    std::cout << "length of vertex weights array is not equal to number of my vertices";
    std::cout << std::endl;
    return -1;
  }

  if (numVWgts > 0){
    vgid = new int [numVWgts];
    vwgt = new float [numVWgts];

    costs.getVertexWeights(numVWgts, vgid, vwgt);

    delete [] vgid;
  }

  balance = compute_balance(comm, myRows, numVWgts, vwgt);

  if (vwgt) delete [] vwgt;

  // Compute the measures based on cut edges

  int haveEdgeWeights = costs.haveGraphEdgeWeights();

  int localNumCuts = 0;
  double localCutWgt = 0.0;
  double localCutn = 0.0;
  double localCutl = 0.0;

  int maxEdges = graph.MaxNumIndices();

  if (maxEdges > 0){
    int *nborID  = new int [maxEdges];
    int numEdges;
    std::map<int, float> weightMap;

    const Epetra_BlockMap &rowmap = graph.RowMap();
    const Epetra_BlockMap &colmap = graph.ColMap();

    // Get the processes owning my vertices neighbors

    int numCols = colmap.NumMyElements();
    const int *colGIDs = colmap.MyGlobalElements();
    int *nborProc_GID = new int [numCols];
    int *nborRow_LID = new int [numCols];

    rc = rowmap.RemoteIDList(numCols, colGIDs, nborProc_GID, nborRow_LID);

    if (rc != 0){
      std::cout << "Error obtaining remote process ID list";
      std::cout << std::endl;
      delete [] nborProc_GID;
      delete [] nborRow_LID;
      delete [] nborID;
      return -1;
    }

    std::map<int, int> colProc;
    std::map<int, int>::iterator procIter;

    for (int j=0; j<numCols; j++){

      // map from column GID to process owning row with that GID 
      //   (matrix is square)

      colProc[colGIDs[j]] = nborProc_GID[j];
    }
    delete [] nborProc_GID;
    delete [] nborRow_LID;

    for (int i=0; i < rowmap.NumMyElements(); i++){
      int vtxGID = rowmap.GID(i);

      if (haveEdgeWeights){
        costs.getGraphEdgeWeights(vtxGID, weightMap);
      }

      int numEdges = graph.NumMyIndices(i);

      if (numEdges > 0){

        // get neighboring vertices

        graph.ExtractMyRowCopy(i, maxEdges, numEdges, nborID);

        for (int j=0; j<numEdges; j++){
          nborID[j] = colGIDs[nborID[j]];  // convert to global ID
        }
        
        // get processes that own my neighbors
  
        std::set<int> nbors;
        float heWeight = 0.0;

        for (int j=0; j < numEdges; j++){

          if (nborID[j] == vtxGID) continue;  // skip self edges

          procIter = colProc.find(nborID[j]);
          if (procIter == colProc.end()){
            std::cout << "process owning column is missing";
            std::cout << std::endl;
            delete [] nborID;
            return -1;
          }
          int procNum = procIter->second;

          float wgt = 1.0;
          if (haveEdgeWeights){
            std::map<int, float>::iterator curr = weightMap.find(nborID[j]);
            if (curr == weightMap.end()){
              std::cout << "Graph edge weights do not match matrix";
              std::cout << std::endl;
              delete [] nborID;
              return -1;
            }
            wgt = curr->second;
          }
    
          if (procNum != myProc){
            localNumCuts++;            // number of graph edges that are cut 
            nbors.insert(procNum);     // count number of neighboring processes
            localCutWgt += wgt;        // sum of weights of cut edges
          }
          heWeight += wgt;             // implied hyperedge weight
        }
        int numNbors = nbors.size();

        if (numNbors > 0){
          // sum of the implied hyperedge weights of cut hyperedges
          localCutn += heWeight;   

          // sum of (number of partitions - 1) weighted by the 
          // implied hyperedge weight
          localCutl += (numNbors * heWeight);
        }
      }
    } // next vertex in my partition

    delete [] nborID;
  }

  double lval[4], gval[4];

  lval[0] = (double)localNumCuts;
  lval[1] = localCutWgt;
  lval[2] = localCutn;
  lval[3] = localCutl;

  comm.SumAll(lval, gval, 4);

  numCuts = (int)gval[0];
  cutWgt = gval[1];
  cutn   = gval[2];
  cutl   = gval[3];

  return 0;
}
int compute_hypergraph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &balance, double &cutn, double &cutl)  // output
{
  const Epetra_Comm &comm  = graph.Comm();
#ifdef HAVE_MPI
  const Epetra_MpiComm* mpiComm =
    dynamic_cast<const Epetra_MpiComm*>(&comm);

  MPI_Comm mcomm = mpiComm->Comm();
#endif
  int nProcs = comm.NumProc();
  int myProc = comm.MyPID();
  int myRows = graph.NumMyRows();
  int rc;
  int *vgid = NULL;
  float *vwgt = NULL;

  int numVWgts = costs.getNumVertices();

  if ((numVWgts > 0) && (numVWgts != myRows)){
    std::cout << "length of row (vertex) weights array is not equal to number of rows";
    std::cout << std::endl;
    return -1;
  }

  if (numVWgts > 0){
    vgid = new int [numVWgts];
    vwgt = new float [numVWgts];

    costs.getVertexWeights(numVWgts, vgid, vwgt);

    delete [] vgid;
  }

  balance = compute_balance(comm, myRows, numVWgts, vwgt);

  if (vwgt) delete [] vwgt;

  /* Compute cutl and cutn. 
   */

  int totalHEWeights = 0; 
  int numCols = graph.NumGlobalCols();

  int numHEWeights = costs.getNumHypergraphEdgeWeights();

  comm.SumAll(&numHEWeights, &totalHEWeights, 1);
 
  if ((totalHEWeights > 0) && (totalHEWeights !=  numCols)){
    if (myProc == 0)
      std::cerr << "Must supply either no h.e. weights or else supply one for each column" << std::endl;
      return -1;
  }

  int *heGIDs = NULL;
  float *heWeights = NULL;

  if (numHEWeights){
    heGIDs = new int [numHEWeights];
    heWeights = new float [numHEWeights];

    costs.getHypergraphEdgeWeights(numHEWeights, heGIDs, heWeights);
  }

  // Create a map from column global IDs to edge weight.  We assume each
  // edge weight is supplied by only one process.  We don't do the
  // ZOLTAN_EDGE_WEIGHT_OP operation.  TODO

  std::map<int, double> heWgt;
  std::map<int, double>::iterator heWgtIter;

  if (numHEWeights){
    for (int j=0; j<numHEWeights; j++){
      heWgt[heGIDs[j]] = heWeights[j];
    }
    delete [] heGIDs;
    delete [] heWeights;
  }

  // Create a set containing all the columns in my rows.  We assume all
  // the rows are in the same partition.

  const Epetra_BlockMap &colMap = graph.ColMap();
  int numMyCols = colMap.NumMyElements();

  std::set<int> colGIDS;
  std::set<int>::iterator gidIter;

  for (int j=0; j<numMyCols; j++){
    colGIDS.insert(colMap.GID(j));
  }
  
  /* Divide columns among processes, then each process computes its
   * assigned columns' cutl and cutn.
   */
  int ncols = numCols / nProcs;
  int leftover = numCols - (nProcs * ncols);
  std::vector<int> colCount(nProcs, 0);
  for (int i=0; i<nProcs; i++){
    colCount[i] = ncols;
    if (i < leftover) colCount[i]++;
  }
  int *colTotals = NULL;
  double *colWeights = NULL;
  if (colCount[myProc] > 0){
    colTotals = new int [colCount[myProc]];
    if (totalHEWeights > 0){
      colWeights = new double [colCount[myProc]];
    } 
  }
  int *colLocal= new int [ncols + 1];
  double *localWeights = NULL;
  if (totalHEWeights > 0){
    localWeights = new double [ncols + 1];
  }

  int base = colMap.IndexBase();
  int colStart = base;

  for (int i=0; i<nProcs; i++){

    // All processes send info to the process reponsible
    // for the next group of columns

    int ncols = colCount[i];
    int colEnd = colStart + ncols;
    for (int j=colStart,k=0; j < colEnd; j++,k++){
      gidIter = colGIDS.find(j);
      if (gidIter != colGIDS.end()){
        colLocal[k] = 1;     // column j has rows in my partition
      }
      else{
        colLocal[k] = 0;
      }
      if (totalHEWeights > 0){
        heWgtIter = heWgt.find(j);
        if (heWgtIter != heWgt.end()){
          // I have the edge weight for column j
          localWeights[k] = heWgtIter->second;
        }
        else{
          localWeights[k] = 0.0;
        }
      }
      
    }
#ifdef HAVE_MPI
    rc = MPI_Reduce(colLocal, colTotals, ncols, MPI_INT, MPI_SUM, i, mcomm);
    if (totalHEWeights > 0){
      rc = MPI_Reduce(localWeights, colWeights, ncols, MPI_DOUBLE, MPI_SUM, i, mcomm);
    }
#else
    memcpy(colTotals, colLocal, ncols * sizeof(int));
    if (totalHEWeights > 0){
      memcpy(colWeights, localWeights, ncols * sizeof(double));
    }
#endif
    colStart = colEnd;
  }

  delete [] colLocal;
  if (localWeights) delete [] localWeights;

  double localCutN=0;
  double localCutL=0;
  double ewgt = 1.0;

  for (int j=0; j<colCount[myProc]; j++){
    if (totalHEWeights > 0){
      ewgt = colWeights[j];
    }
    if (colTotals[j] > 1){
      localCutL += (colTotals[j] - 1) * ewgt; // # of cuts in columns/edges
      localCutN += ewgt;                      // # of cut columns/edges
    }
  }
  if (colTotals) delete [] colTotals;
  if (colWeights) delete [] colWeights;

  comm.SumAll(&localCutN, &cutn, 1);
  comm.SumAll(&localCutL, &cutl, 1);

  return 0;
}
static double compute_balance(const Epetra_Comm &comm, int myRows, int nwgts, float *wgts)
{
  int nProcs = comm.NumProc();
  int myProc = comm.MyPID();
  double weightTotal, balance;

  /* Proportion of weight desired in each partition.  For now we
     have no interface to specify unequal partitions.
   */
  double partSize = 1.0 / nProcs;
  std::vector<double> partSizes(nProcs, partSize);

  /* Sum of my row weights.  
   */
  double weightLocal = 0.0;

  if (nwgts > 0){
    for (int i=0; i<myRows; i++){
      weightLocal += wgts[i];
    }
  }
  else{
    weightLocal += myRows;   // default weight of each vertex is 1.0
  }

  comm.SumAll(&weightLocal, &weightTotal, 1);

  /* My degree of imbalance
   */
  double goalWeight = partSizes[myProc] * weightTotal;
  double imbalance = 1.0;
  if (weightLocal >= goalWeight)
    imbalance += (weightLocal - goalWeight) / goalWeight;
  else
    imbalance += (goalWeight - weightLocal) / goalWeight;

  comm.MaxAll(&imbalance, &balance, 1);

  return balance;
}
// Print out the matrix showing the partitioning, for debugging purposes.  
// This only works for small example matrices and 10 or fewer processes.

void show_matrix(const char *txt, const Epetra_CrsGraph &graph, const Epetra_Comm &comm)
{
  int me = comm.MyPID();

  if (comm.NumProc() > 10){
    if (me == 0){
      std::cout << txt << std::endl;
      std::cout << "Printed matrix format only works for 10 or fewer processes" << std::endl;
    }
    return;
  }

  const Epetra_BlockMap &rowmap = graph.RowMap();
  const Epetra_BlockMap &colmap = graph.ColMap();

  int myRows = rowmap.NumMyElements();
  int numRows = graph.NumGlobalRows();
  int numCols = graph.NumGlobalCols();
  int base = rowmap.IndexBase();

  int *myA = new int [numRows * numCols];
  int *A = new int [numRows * numCols];
  memset(myA, 0, sizeof(int) * numRows * numCols);

  int *myIndices;

  int *myRowGIDs = rowmap.MyGlobalElements();

  for (int i=0; i< myRows; i++){
    int myRowLID = rowmap.LID(myRowGIDs[i]);

    int numEntries = graph.NumMyIndices(myRowLID);

    if (numEntries > 0){
      int rc = graph.ExtractMyRowView(myRowLID, numEntries, myIndices);
      if (rc){
        std::cout << txt << std::endl;
        std::cout << "extract graph error" << std::endl;
        return;
      }

      int *row = myA + (numCols * (myRowGIDs[i] - base));

      for (int j=0; j < numEntries; j++){
        int gid = colmap.GID(myIndices[j]);
        row[gid-base] = me+1;
      }
    }
  }

  comm.SumAll(myA, A, numRows * numCols);

  if (me == 0){
    std::cout << txt << std::endl;

    std::cout << "  ";
    for (int j=0; j<numCols; j++){
      std::cout << j%10 ;
    }
    std::cout << std::endl;

    int *row = A;
 
    for (int i=0; i < numRows; i++, row += numCols){
      std::cout << i%10 << " ";
      for (int j=0; j < numCols; j++){
        if (row[j] > 0){
          std::cout << row[j]-1;
        }
        else{
          std::cout << " ";
        }
      }
      std::cout << " " << i%10 << std::endl;
    }
    std::cout << "  ";
    for (int j=0; j<numCols; j++){
      std::cout << j%10 ;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

  delete [] myA;
  delete [] A;
}

#endif // HAVE_EPETRA
}

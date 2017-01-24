//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//************************************************************************
//@HEADER

#include <ispatest_lbeval_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
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

namespace ispatest {

#ifdef HAVE_EPETRA

static void printMatrix(const char *txt, int *myA, int *myX, int *myB,
                        int numRows, int numCols, const Epetra_Comm &comm);
static int make_my_A(const Epetra_RowMatrix &matrix, int *myA, const Epetra_Comm &comm);

static int compute_hypergraph_metrics(const Epetra_BlockMap &rowmap,
            const Epetra_BlockMap &colmap,
            int numGlobalColumns,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, double &cutn, double &cutl);

static int compute_graph_metrics(const Epetra_BlockMap &rowmap,
            const Epetra_BlockMap &colmap,
            std::vector<std::vector<int> > &rows,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl);

/******************************************************************
  Compute weight balance
******************************************************************/
int compute_balance(const Epetra_Vector &wgts, double myGoalWeight,
                      double &min, double &max, double &avg)
{
  if ((myGoalWeight < 0) || (myGoalWeight > 1.0)){
    std::cerr << "compute_balance: Goal weight should be in the range [0, 1]" << std::endl;
    return -1;
  }

  double weightTotal;
  wgts.Norm1(&weightTotal);

  double weightLocal = 0.0;
  for (int i=0; i < wgts.MyLength(); i++){
    weightLocal += wgts[i];
  }

  /* My degree of imbalance.
   * If myGoalWeight is zero, I'm in perfect balance since I got what I wanted.
   */
  double goalWeight = myGoalWeight * weightTotal;
  double imbalance = 1.0;

  if (myGoalWeight > 0.0){
    if (weightLocal >= goalWeight)
      imbalance += (weightLocal - goalWeight) / goalWeight;
    else
      imbalance += (goalWeight - weightLocal) / goalWeight;
  }

  const Epetra_Comm &comm = wgts.Comm();

  comm.MaxAll(&imbalance, &max, 1);
  comm.MinAll(&imbalance, &min, 1);
  comm.SumAll(&imbalance, &avg, 1);

  avg /= comm.NumProc();

  return 0;
}

/******************************************************************
  Compute graph metrics
******************************************************************/

int compute_graph_metrics(const Epetra_RowMatrix &matrix,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl)
{

  const Epetra_Map &rmap = matrix.RowMatrixRowMap();
  const Epetra_Map &cmap = matrix.RowMatrixColMap();

  int maxEdges = matrix.MaxNumEntries();

  std::vector<std::vector<int> > myRows(rmap.NumMyElements());

  if (maxEdges > 0){
    int numEdges = 0;
    int *nborLID  = new int [maxEdges];
    double *tmp = new double [maxEdges];

    for (int i=0; i<rmap.NumMyElements(); i++){

      matrix.ExtractMyRowCopy(i, maxEdges, numEdges, tmp, nborLID);
      std::vector<int> cols(numEdges);
      for (int j=0; j<numEdges; j++){
        cols[j] = nborLID[j];
      }
      myRows[i] = cols;
    }
    delete [] nborLID;
    delete [] tmp;
  }

  return compute_graph_metrics(rmap, cmap, myRows, costs, myGoalWeight,
                               balance, numCuts, cutWgt, cutn, cutl);

}

int compute_graph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl)
{
  const Epetra_BlockMap & rmap = graph.RowMap();
  const Epetra_BlockMap & cmap = graph.ColMap();

  int maxEdges = cmap.NumMyElements();

  std::vector<std::vector<int> > myRows(rmap.NumMyElements());

  if (maxEdges > 0){
    int numEdges = 0;
    int *nborLID  = new int [maxEdges];
    for (int i=0; i<rmap.NumMyElements(); i++){

      graph.ExtractMyRowCopy(i, maxEdges, numEdges, nborLID);
      std::vector<int> cols(numEdges);
      for (int j=0; j<numEdges; j++){
        cols[j] = nborLID[j];
      }
      myRows[i] = cols;
    }
    delete [] nborLID;
  }

  return compute_graph_metrics(rmap, cmap, myRows, costs, myGoalWeight,
                               balance, numCuts, cutWgt, cutn, cutl);

}

static int compute_graph_metrics(const Epetra_BlockMap &rowmap, const Epetra_BlockMap &colmap,
            std::vector<std::vector<int> > &rows,
            Isorropia::Epetra::CostDescriber &costs, double &myGoalWeight,
            double &balance, int &numCuts, double &cutWgt, double &cutn, double &cutl)
{
  const Epetra_Comm &comm  = rowmap.Comm();
  int myProc = comm.MyPID();
  int myCols = colmap.NumMyElements();
  double min, avg;
  std::map<int, float> vertexWeights;
  std::map<int, std::map<int, float > > graphEdgeWeights;
  std::map<int, float> hyperEdgeWeights;

  costs.getCosts(vertexWeights,  // vertex global ID -> weight
           graphEdgeWeights,    //  vertex global ID -> map from neighbor global ID to edge weight
           hyperEdgeWeights);   // hyperedge global ID -> weight

  // Compute the balance

  Epetra_Vector vwgt(rowmap);
  int numVWgts = vertexWeights.size();

  if (numVWgts > 0){
    double *wvals = new double [numVWgts];
    int *gids = new int [numVWgts];

    std::map<int, float>::iterator vnext = vertexWeights.begin();
    int i=0;
    while (vnext != vertexWeights.end()){
      wvals[i] =  vnext->second;
      gids[i] = vnext->first;
      vnext++;
      i++;
    }

    vwgt.ReplaceGlobalValues(i, wvals, gids);

    delete [] wvals;
    delete [] gids;
  }
  else{
    vwgt.PutScalar(1.0);   // default to unit weights
  }

  compute_balance(vwgt, myGoalWeight, min, balance, avg);

  if (balance < 0){
    return 1;
  }

  // Compute the measures based on cut edges

  int *procID = new int [myCols];
  int *GID = new int [myCols];
  int *tmp = new int [myCols];

  for (int i=0; i < myCols; i++){
    GID[i] = colmap.GID(i);
  }

  rowmap.RemoteIDList(myCols, GID, procID, tmp);  // matrix is square

  delete [] tmp;

  int haveEdgeWeights = graphEdgeWeights.size();

  int localNumCuts = 0;
  double localCutWgt = 0.0;
  double localCutn = 0.0;
  double localCutl = 0.0;

  for (int i=0; i < rowmap.NumMyElements(); i++){
    int vtxGID = rowmap.GID(i);

    int numEdges = rows[i].size();

    if (numEdges > 0){

      std::map<int, std::map<int, float> >::iterator wnext;
      if (haveEdgeWeights){
        wnext = graphEdgeWeights.find(vtxGID);
        if (wnext == graphEdgeWeights.end()){
          std::cerr << "Graph edge weights are missing for vertex " << vtxGID;
          std::cerr << std::endl;
          return -1;
        }
      }

      double heWeight = 0.0;
      std::set<int> nbors;

      for (int j=0; j < numEdges; j++){

        int colGID = GID[rows[i][j]];

        int nborProc = procID[rows[i][j]];

        if (colGID == vtxGID) continue;  // skip self edges

        float wgt = 1.0;
        if (haveEdgeWeights){
          std::map<int, float>::iterator curr = (wnext->second).find(colGID);
          if (curr == (wnext->second).end()){
            std::cerr << "Graph edge weights do not match matrix";
            std::cerr << std::endl;
            return -1;
          }
          wgt = curr->second;
        }

        if (nborProc != myProc){
          localNumCuts++;            // number of graph edges that are cut
          nbors.insert(nborProc);     // count number of neighboring processes
          localCutWgt += wgt;        // sum of weights of cut edges
        }
        heWeight += wgt;             // implied hyperedge weight, sum all edges
      }
      int numNbors = nbors.size();

      if (numNbors > 0){
        // implied hyperedge is vertex and neighbors, if cut, add in its he weight
        localCutn += heWeight;

        // sum of (number of partitions - 1) weighted by the
        // implied hyperedge weight
        localCutl += (numNbors * heWeight);
      }
    }
  } // next vertex in my partition

  delete [] GID;
  delete [] procID;

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

/******************************************************************
  Compute hypergraph metrics
******************************************************************/

int compute_hypergraph_metrics(const Epetra_RowMatrix &matrix,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, double &cutn, double &cutl)  // output
{
  const Epetra_BlockMap &rmap =
    static_cast<const Epetra_BlockMap &>(matrix.RowMatrixRowMap());

  const Epetra_BlockMap &cmap =
    static_cast<const Epetra_BlockMap &>(matrix.RowMatrixColMap());

  return compute_hypergraph_metrics(rmap, cmap,
                                     matrix.NumGlobalCols(),
                                     costs,
                                     myGoalWeight,
                                     balance, cutn, cutl);

}
int compute_hypergraph_metrics(const Epetra_CrsGraph &graph,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, double &cutn, double &cutl)  // output
{
  return compute_hypergraph_metrics(graph.RowMap(), graph.ColMap(),
                                     graph.NumGlobalCols(), costs,
                                     myGoalWeight,
                                     balance, cutn, cutl);
}
static int compute_hypergraph_metrics(const Epetra_BlockMap &rowmap,
            const Epetra_BlockMap &colmap,
            int numGlobalColumns,
            Isorropia::Epetra::CostDescriber &costs,
            double &myGoalWeight,
            double &balance, double &cutn, double &cutl)  // output
{
  const Epetra_Comm &comm  = rowmap.Comm();
#ifdef HAVE_MPI
  const Epetra_MpiComm* mpiComm =
    dynamic_cast<const Epetra_MpiComm*>(&comm);

  MPI_Comm mcomm = mpiComm->Comm();
#endif
  int nProcs = comm.NumProc();
  int myProc = comm.MyPID();
  double min, avg;
  std::map<int, float> vertexWeights;
  std::map<int, std::map<int, float > > graphEdgeWeights;
  std::map<int, float> hyperEdgeWeights;

  costs.getCosts(vertexWeights,  // vertex global ID -> weight
           graphEdgeWeights,    //  vertex global ID -> map from neighbor global ID to edge weight
           hyperEdgeWeights);   // hyperedge global ID -> weight

  Epetra_Vector vwgt(rowmap);
  int numVWgts = vertexWeights.size();

  if (numVWgts > 0){
    double *wvals = new double [numVWgts];
    int *gids = new int [numVWgts];

    std::map<int, float>::iterator vnext = vertexWeights.begin();
    int i=0;
    while (vnext != vertexWeights.end()){
      wvals[i] =  vnext->second;
      gids[i] = vnext->first;
      vnext++;
      i++;
    }

    vwgt.ReplaceGlobalValues(i, wvals, gids);

    delete [] wvals;
    delete [] gids;
  }
  else{
    vwgt.PutScalar(1.0);   // default to unit weights
  }

  compute_balance(vwgt, myGoalWeight, min, balance, avg);

  if (balance < 0){
    return 1;
  }

  /* Compute cutl and cutn.
   */

  int totalHEWeights = 0;

  int numHEWeights = hyperEdgeWeights.size();

  comm.SumAll(&numHEWeights, &totalHEWeights, 1);

  if ((totalHEWeights > 0) && (totalHEWeights <  numGlobalColumns)){
    if (myProc == 0)
      std::cerr << "Must supply either no h.e. weights or else supply at least one for each column" << std::endl;
    return -1;
  }

  std::map<int, float>::iterator heWgtIter;

  // Create a set containing all the columns in my rows.  We assume all
  // the rows are in the same partition.

  int numMyCols = colmap.NumMyElements();

  std::set<int> colGIDS;
  std::set<int>::iterator gidIter;

  for (int j=0; j<numMyCols; j++){
    colGIDS.insert(colmap.GID(j));
  }

  /* Divide columns among processes, then each process computes its
   * assigned columns' cutl and cutn.
   *  TODO - numGlobalColumns can be less than nprocs

   * Fix this when a process is assigned no columns. TODO
   */
  int ncols = numGlobalColumns / nProcs;
  int leftover = numGlobalColumns - (nProcs * ncols);
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

  int base = colmap.IndexBase();
  int colStart = base;

  for (int i=0; i<nProcs; i++){

    // All processes send info to the process reponsible
    // for the next group of columns

    int numCols = colCount[i];
    int colEnd = colStart + numCols;
    for (int j=colStart,k=0; j < colEnd; j++,k++){
      gidIter = colGIDS.find(j);
      if (gidIter != colGIDS.end()){
        colLocal[k] = 1;     // column j has rows in my partition
      }
      else{
        colLocal[k] = 0;
      }
      if (totalHEWeights > 0){
        std::map<int, float>::iterator theHeWgtIter = hyperEdgeWeights.find(j);
        if (theHeWgtIter != hyperEdgeWeights.end()){
          // I have the edge weight for column j
          localWeights[k] = theHeWgtIter->second;
        }
        else{
          localWeights[k] = 0.0;
        }
      }

    }
#ifdef HAVE_MPI
    MPI_Reduce(colLocal, colTotals, numCols, MPI_INT, MPI_SUM, i, mcomm);
    if (totalHEWeights > 0){
      MPI_Reduce(localWeights, colWeights, numCols, MPI_DOUBLE, MPI_SUM, i, mcomm);
    }
    // TODO handle possible MPI error
#else
    memcpy(colTotals, colLocal, numCols * sizeof(int));
    if (totalHEWeights > 0){
      memcpy(colWeights, localWeights, numCols * sizeof(double));
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

// Print out the graph or linear system showing the partitioning, for debugging or
// educational purposes.
// The display only makes sense for 10 or fewer processes, and fairly small
// graphs or linear systems.

void show_matrix(const char *txt, const Epetra_CrsGraph &graph, const Epetra_Comm &comm)
{
  int me = comm.MyPID();

  if (comm.NumProc() > 10){
    if (me == 0){
      std::cerr << txt << std::endl;
      std::cerr << "Printed matrix format only works for 10 or fewer processes" << std::endl;
    }
    return;
  }

  const Epetra_BlockMap &rowmap = graph.RowMap();
  const Epetra_BlockMap &colmap = graph.ColMap();

  int myRows = rowmap.NumMyElements();
  int numRows = graph.NumGlobalRows();
  int numCols = graph.NumGlobalCols();
  int base = rowmap.IndexBase();

  if ((numRows > 200) || (numCols > 500)){
    if (me == 0){
      std::cerr << txt << std::endl;
      std::cerr << "show_matrix: problem is too large to display" << std::endl;
    }
    return;
  }

  int *myA = new int [numRows * numCols];
  memset(myA, 0, sizeof(int) * numRows * numCols);

  int *myIndices;

  int *myRowGIDs = rowmap.MyGlobalElements();

  for (int i=0; i< myRows; i++){
    int myRowLID = rowmap.LID(myRowGIDs[i]);

    int numEntries = graph.NumMyIndices(myRowLID);

    if (numEntries > 0){
      int rc = graph.ExtractMyRowView(myRowLID, numEntries, myIndices);
      if (rc){
        std::cerr << txt << std::endl;
        std::cerr << "extract graph error" << std::endl;
        return;
      }

      int *row = myA + (numCols * (myRowGIDs[i] - base));

      for (int j=0; j < numEntries; j++){
        int gid = colmap.GID(myIndices[j]);
        row[gid-base] = me+1;
      }
    }
  }

  printMatrix(txt, myA, NULL, NULL, numRows, numCols, comm);

  delete [] myA;
}

void show_matrix(const char *txt, const Epetra_RowMatrix &matrix, const Epetra_Comm &comm)
{
  int me = comm.MyPID();
  if (comm.NumProc() > 10){
    if (me == 0){
      std::cout << txt << std::endl;
      std::cout << "Printed matrix format only works for 10 or fewer processes" << std::endl;
    }
    return;
  }

  int numRows = matrix.NumGlobalRows();
  int numCols = matrix.NumGlobalCols();

  if ((numRows > 200) || (numCols > 500)){
    if (me == 0){
      std::cerr << txt << std::endl;
      std::cerr << "show_matrix: problem is too large to display" << std::endl;
    }
    return;
  }

  int *myA = new int [numRows * numCols];

  make_my_A(matrix, myA, comm);

  printMatrix(txt, myA, NULL, NULL, numRows, numCols, comm);

  delete [] myA;
}
void show_matrix(const char *txt, const Epetra_LinearProblem &problem, const Epetra_Comm &comm)
{
  int me = comm.MyPID();

  if (comm.NumProc() > 10){
    if (me == 0){
      std::cout << txt << std::endl;
      std::cout << "Printed matrix format only works for 10 or fewer processes" << std::endl;
    }
    return;
  }

  Epetra_RowMatrix *matrix = problem.GetMatrix();
  Epetra_MultiVector *lhs = problem.GetLHS();
  Epetra_MultiVector *rhs = problem.GetRHS();

  int numRows = matrix->NumGlobalRows();
  int numCols = matrix->NumGlobalCols();

  if ((numRows > 200) || (numCols > 500)){
    if (me == 0){
      std::cerr << txt << std::endl;
      std::cerr << "show_matrix: problem is too large to display" << std::endl;
    }
    return;
  }

  int *myA = new int [numRows * numCols];

  make_my_A(*matrix, myA, comm);

  int *myX = new int [numCols];
  int *myB = new int [numRows];

  memset(myX, 0, sizeof(int) * numCols);
  memset(myB, 0, sizeof(int) * numRows);

  const Epetra_BlockMap &lhsMap = lhs->Map();
  const Epetra_BlockMap &rhsMap = rhs->Map();

  int base = lhsMap.IndexBase();

  for (int j=0; j < lhsMap.NumMyElements(); j++){
    int colGID = lhsMap.GID(j);
    myX[colGID - base] = me + 1;
  }

  for (int i=0; i < rhsMap.NumMyElements(); i++){
    int rowGID = rhsMap.GID(i);
    myB[rowGID - base] = me + 1;
  }

  printMatrix(txt, myA, myX, myB, numRows, numCols, comm);

  delete [] myA;
  delete [] myX;
  delete [] myB;
}

static int make_my_A(const Epetra_RowMatrix &matrix, int *myA, const Epetra_Comm &comm)
{
  int me = comm.MyPID();

  const Epetra_Map &rowmap = matrix.RowMatrixRowMap();
  const Epetra_Map &colmap = matrix.RowMatrixColMap();

  int myRows = matrix.NumMyRows();
  int numRows = matrix.NumGlobalRows();
  int numCols = matrix.NumGlobalCols();
  int base = rowmap.IndexBase();
  int maxRow = matrix.MaxNumEntries();

  memset(myA, 0, sizeof(int) * numRows * numCols);

  int *myIndices = new int [maxRow];
  double *tmp = new double [maxRow];

  int rowLen = 0;

  for (int i=0; i< myRows; i++){

    int rc = matrix.ExtractMyRowCopy(i, maxRow, rowLen, tmp, myIndices);

    if (rc){
      if (me == 0){
        std::cout << "Error in make_my_A" << std::endl;
      }
       return 1;
    }

    int *row = myA + (numCols * (rowmap.GID(i) - base));

    for (int j=0; j < rowLen; j++){

      int colGID = colmap.GID(myIndices[j]);

      row[colGID - base] = me + 1;
    }
  }

  if (maxRow){
    delete [] myIndices;
    delete [] tmp;
  }
  return 0;
}

static void printMatrix(const char *txt, int *myA, int *myX, int *myB,
                        int numRows, int numCols, const Epetra_Comm &comm)
{
  int me = comm.MyPID();

  int *A = new int [numRows * numCols];
  int *x = NULL;
  int *b = NULL;

  comm.SumAll(myA, A, numRows * numCols);

  if (myX){
    x = new int [numCols];
    comm.SumAll(myX, x, numCols);
  }
  if (myB){
    b = new int [numRows];
    comm.SumAll(myB, b, numRows);
  }

  if (me == 0){
    std::cout << txt << std::endl;

    std::cout << "  ";
    for (int j=0; j<numCols; j++){
      std::cout << j%10 ;
    }
    if (x)
      std::cout << "    LHS";
    if (b)
      std::cout << "        RHS";

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
      std::cout << " " << i%10 ;

      if (x && (i < numCols)){
        std::cout << "   " << x[i]-1;
      }
      if ((i == 0) && b){
        std::cout << "    =  [";
        for (int j=0; j<numRows; j++){
          std::cout << b[j]-1;
        }
        std::cout << "]";
      }
      std::cout << std::endl;

    }
    std::cout << "  ";
    for (int j=0; j<numCols; j++){
      std::cout << j%10 ;
    }

    int columnsRemaining = numCols - numRows;
    int next_x = numRows;

    if ((columnsRemaining > 0) && x){
      std::cout << "     " << x[next_x++] - 1 << std::endl;
      columnsRemaining--;
      int pad = numCols + 7;
      while(columnsRemaining){
        for( int i=0; i < pad; i++){
          std::cout << " ";
        }
        std::cout << x[next_x++] - 1 << std::endl;
        columnsRemaining--;
      }
    }
    std::cout << std::endl;
  }

  delete [] A;
  if (x) delete [] x;
  if (b) delete [] b;
}

#endif // HAVE_EPETRA
}

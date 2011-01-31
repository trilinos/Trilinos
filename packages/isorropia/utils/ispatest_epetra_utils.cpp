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

#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <vector>

#include <Isorropia_Exception.hpp>

#ifdef HAVE_EPETRA

#include <ispatest_epetra_utils.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Import.h>


namespace ispatest {

/***************************************************************
  Epetra_RowMatrix functions
****************************************************************/

int fill_matrix(Epetra_CrsMatrix& matrix,
                int numNonzerosPerRow,
                bool verbose)
{
  int err = 0;
  const Epetra_Map& rowmap = matrix.RowMap();
  int num_my_rows = rowmap.NumMyElements();
  int num_global_rows = rowmap.NumGlobalElements();

  std::vector<int> indices(numNonzerosPerRow);
  std::vector<double> coefs(numNonzerosPerRow);

  for(int i=0; i<num_my_rows; ++i) {
    int global_row = rowmap.GID(i);
    int first_col = global_row - numNonzerosPerRow/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (num_global_rows - numNonzerosPerRow)) {
      first_col = num_global_rows - numNonzerosPerRow;
    }

    for(int j=0; j<numNonzerosPerRow; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    err = matrix.InsertGlobalValues(global_row, numNonzerosPerRow,
                                    &coefs[0], &indices[0]);
    if (err < 0) {
      err = matrix.ReplaceGlobalValues(global_row, numNonzerosPerRow,
                                      &coefs[0], &indices[0]);
      if (err < 0) {
        return(err);
      }
    }
  }

  err = matrix.FillComplete();
  return(err);
}

bool test_matrix_vector_multiply(Epetra_CrsMatrix &A)
{
  const Epetra_Map &xmap = A.DomainMap();
  const Epetra_Map &ymap = A.RangeMap();

  int myLen = xmap.NumMyElements();
  double *val = NULL;
  if (myLen > 0){
    val = new double [myLen];
    for (int i=0; i < myLen; i+=2){
      val[i] = 1.0;
    }
    for (int i=1; i < myLen; i+=2){
      val[i] = -1.0;
    }
  }
  Epetra_Vector x(Copy, xmap, val);

  if (val){
    delete [] val;
  }

  Epetra_Vector y(ymap, true);

  // See if Ax=y completes without error

  int fail = A.Multiply(false, x, y);

  if (!fail){
    // Try again with A transpose
    fail = A.Multiply(true, y, x);
  }

  return (!fail);
}

bool test_row_matrix_vector_multiply(Epetra_RowMatrix &A)
{
  const Epetra_Map &xmap = A.OperatorDomainMap();
  const Epetra_Map &ymap = A.OperatorRangeMap(); // same as A.RowMatrixRowMap()

  int myLen = xmap.NumMyElements();
  double *val = NULL;
  if (myLen > 0){
    val = new double [myLen];
    for (int i=0; i < myLen; i+=2){
      val[i] = 1.0;
    }
    for (int i=1; i < myLen; i+=2){
      val[i] = -1.0;
    }
  }
  Epetra_Vector x(Copy, xmap, val);

  if (val){
    delete [] val;
  }

  Epetra_Vector y(ymap, true);

  // See if Ax=y completes without error

  int fail = A.Multiply(false, x, y);

  if (!fail){
    // Try again with A transpose
    fail = A.Multiply(true, y, x);
  }

  return (!fail);
}

/***************************************************************
  Epetra_CrsGraph functions
****************************************************************/

int fill_graph(Epetra_CrsGraph& graph,
                int numNonzerosPerRow,
                bool verbose)
{
  int err = 0;
  const Epetra_BlockMap& rowmap = graph.RowMap();
  int num_my_rows = rowmap.NumMyElements();
  int num_global_rows = rowmap.NumGlobalElements();

  std::vector<int> indices(numNonzerosPerRow);
  std::vector<double> coefs(numNonzerosPerRow);

  for(int i=0; i<num_my_rows; ++i) {
    int global_row = rowmap.GID(i);
    int first_col = global_row - numNonzerosPerRow/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (num_global_rows - numNonzerosPerRow)) {
      first_col = num_global_rows - numNonzerosPerRow;
    }

    for(int j=0; j<numNonzerosPerRow; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    err = graph.InsertGlobalIndices(global_row, numNonzerosPerRow,
                                     &indices[0]);
    if (err < 0) {
      return(err);
    }
  }

  err = graph.FillComplete();
  return(err);
}

bool test_matrix_vector_multiply(Epetra_CrsGraph &G)
{
  Epetra_CrsMatrix A(Copy, G);
  A.PutScalar(1.0);
  Epetra_Map *domainMap = map_from_blockmap(G.DomainMap());
  Epetra_Map *rangeMap = map_from_blockmap(G.RangeMap());
  A.FillComplete(*domainMap, *rangeMap);

  delete domainMap;
  delete rangeMap;

  return test_matrix_vector_multiply(A);
}
Epetra_Map *map_from_blockmap(const Epetra_BlockMap &b)
{
  int base = b.IndexBase();
  int size = b.NumGlobalElements();
  int mysize = b.NumMyElements();
  int *elts = b.MyGlobalElements();
  const Epetra_Comm &comm = b.Comm();

  Epetra_Map *map = new Epetra_Map(size, mysize, elts, base, comm);

  return map;
} 
/***************************************************************
  Epetra_LinearProblem functions
****************************************************************/
bool test_matrix_vector_multiply(Epetra_LinearProblem &LP)
{
  // Note this test multiplies the matrix by a made up vector,
  // not a vector that is part of the LP

  Epetra_RowMatrix *A = LP.GetMatrix();
  return test_row_matrix_vector_multiply(*A);
}

/***************************************************************
  Epetra_MultiVector functions
****************************************************************/

double unitWeights(const int id, const int me, const int nids, const int nprocs)
{
  return 1.0;
}
double veeWeights(const int id, const int me, const int nids, const int nprocs)
{
  double w = 1.0;

  int mid = nids / 2;
  int diff = id - mid;
  if (id < mid) diff = mid - id;

  w += diff;
  
  return w;
}
double alternateWeights(const int id, const int me, const int nids, const int nprocs)
{
  if (me %2)
    return 1.0;
  else
    return 2.0;
}
Epetra_MultiVector *makeWeights(const Epetra_BlockMap &map, double (*wFunc)(const int, const int, const int, const int))
{
  Epetra_MultiVector *w = new Epetra_MultiVector(map, 1, false);

  if (map.NumMyElements() > 0){
    double *v;
    int stride;
    w->ExtractView(&v, &stride);

    for (int i=0; i<map.NumMyElements(); i++){
      v[i] = wFunc(map.GID(i), map.Comm().MyPID(), 
                   map.NumGlobalElements(), map.Comm().NumProc());
    }
  }

  return w;
}
int readCoordFile(const std::string &fname, 
         std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
{
  std::ifstream f(fname.c_str());

  if (!f){
    std::cerr << "Error trying to read " << fname << std::endl;
    return 0;
  }

  int dim = 0;
  double vals[3], v;
  char line[128];

  x.clear();
  y.clear();
  z.clear();

  while (f.getline(line, 128)){

    // skip blank lines
    std::string checkstring(line);
    if (checkstring.find_first_not_of(" \t\n") == std::string::npos)
      continue;

    std::istringstream is(line);

    is.setf(std::ios::skipws);

    if (!dim){
      is >> v;
      if (!is.fail()){
        is.clear();
        vals[dim++] = v;
        is >> v;
        if (!is.fail()){
          is.clear();
          vals[dim++] = v;
          is >> v;
          if (!is.fail()){
            is.clear();
            vals[dim++] = v;
          }
        }
      }
      else{
        std::cerr << fname << " does not seem to contain coordinates" << std::endl;
        return 0;
      }
    }
    else{
      for (int i=0; i<dim; i++){
        is >> vals[i] ;
      }
    }
    x.push_back(vals[0]);
    if (dim > 1){
      y.push_back(vals[1]);
      if (dim > 2){
        z.push_back(vals[2]);
      }
    }
  }
  return dim;
}

Epetra_MultiVector *file2multivector(const Epetra_Comm &comm, const std::string &fname)
{
  int sizes[2];
  int dim=0, globalSize, mySize;
  std::vector<double> x, y, z;

  if (comm.MyPID() == 0){

    dim = readCoordFile(fname, x, y, z);

    if (dim < 0){
      sizes[0] = sizes[1] = -1;
    }
    else{
      sizes[0] = dim;
      sizes[1] = x.size();
    }

    comm.Broadcast(sizes, 2, 0);
    if (sizes[0] < 0) return NULL;

    globalSize = mySize = x.size();
  }
  else{
    comm.Broadcast(sizes, 2, 0);
    if (sizes[0] < 0) return NULL;

    dim = sizes[0];
    globalSize = sizes[1];
    mySize = 0;
  }

  Epetra_BlockMap map(globalSize, mySize, 1, 0, comm);
  Epetra_MultiVector coords(map, dim, false);

  if (comm.MyPID() == 0){
    double *v;
    int stride, stride2;

    coords.ExtractView(&v, &stride);
    stride2 = 2 * stride;

    for (int i=0; i<globalSize; i++){
      v[0] = x[i];
      if (dim > 1){
        v[stride] = y[i];
        if (dim > 2){
          v[stride2] = z[i];
        }
      }

      v++;
    }

    x.clear();
    y.clear();
    z.clear();
  }

  Epetra_BlockMap newmap(globalSize, 1, 0, comm);  // uniform linear distribution
  Epetra_MultiVector *newcoords = new Epetra_MultiVector(newmap, dim, false);

  Epetra_Import importer(newmap, map);

  newcoords->Import(coords, importer, Insert);

  return newcoords;
}

int printMultiVector(const Epetra_MultiVector &mv, std::ostream &os, const char *s, int max)
{
  const Epetra_Comm &comm = mv.Comm();
  int me = comm.MyPID();
  int nprocs = comm.NumProc();
  int dim = mv.NumVectors();

  const Epetra_BlockMap &map = mv.Map();
  int localSize=0; 
  int globalSize = map.NumGlobalElements();
  int base = map.IndexBase();

  if (globalSize > max){
    if (me == 0){
      os << std::endl << s << std::endl;
      os << "  " << globalSize << " values" << std::endl;
    }

    return 0;
  }

  std::list<int> *gids = NULL;

  // Process 0 needs to know which GIDs are owned by which processes

  int *gidList = new int [globalSize+base];
  int *pidList = new int [globalSize+base];
  int *lidList = new int [globalSize+base];

  if (me > 0){

    map.RemoteIDList(0, gidList, pidList, lidList);
    delete [] lidList;
    delete [] gidList;
    delete [] pidList;
  }
  else {
    gids = new std::list<int> [nprocs];
  
    for (int i = base; i < base+globalSize; i++){
      gidList[i] = i;
    }
  
    map.RemoteIDList(globalSize, gidList+base, pidList+base, lidList+base);

    delete [] lidList;

    for (int i=base; i<globalSize+base; i++){
      gids[pidList[i]].push_back(gidList[i]);
    }
  
    delete [] gidList;
    delete [] pidList;

  }

  // Get all multivector values onto process 0

  if (me == 0){
    localSize = globalSize;
  }

  Epetra_BlockMap newmap(globalSize, localSize, 1, base, comm); 
  Epetra_MultiVector *newmv = new Epetra_MultiVector(newmap, dim, false);
  Epetra_Import importer(newmap, map);
  newmv->Import(mv, importer, Insert);

  if (me > 0){
    return 0;
  }
  
  // Print out multivector info

  double *v, x=0, y=0, z=0;
  int stride, stride2;
  newmv->ExtractView(&v, &stride);
  stride2 = 2 * stride;

  os << std::endl << s << std::endl;

  for (int i=0; i<nprocs; i++){
    std::list<int> procGids = gids[i];
    std::list<int>::iterator curr = procGids.begin();
    std::list<int>::iterator last = procGids.end();

    int nvals = procGids.size();
    int countVals = 0;

    os << "Process " << i << ", " << nvals << " values" << std::endl << "  ";

    while (curr != last){
      int gid = *curr;
      int lid = newmv->Map().LID(gid);

      if (dim > 0){
        x = v[lid]; 
        if (dim > 1){
          y = v[lid+stride]; 
          if (dim > 2){
            z = v[lid+stride2]; 
          }
        }
      }
      if (countVals && (countVals%10==0)) os << std::endl << "  ";
      os << gid << " ";
      if (dim > 0){
        os << "(" << x ;
        if (dim > 1){
          os << " " << y;
          if (dim > 2){
            os << " " << z;
          }
        }
        os << ") ";
      }
      curr++;
      countVals++;
    }
    os << std::endl;
  }
  delete [] gids;

  return 0;
}

int printRowMatrix(const Epetra_RowMatrix &m, std::ostream &os, const char *s, 
                     bool withGraphCuts, int max)
{
  const Epetra_Comm &Comm = m.Comm();
  int me = Comm.MyPID();
  int ncols = m.NumGlobalCols();
  int nrows = m.NumGlobalRows();
  int localRows = m.NumMyRows();
  int len;
  const Epetra_Map &rowmap = m.RowMatrixRowMap();
  const Epetra_Map &colmap = m.RowMatrixColMap();
  int col_base = colmap.IndexBase();
  int row_base = rowmap.IndexBase();

  if (nrows > max){   // printing doesn't make sense for huge matrices
    if (me == 0){
      os << std::endl << s << std::endl;
      os << "  " << nrows << " rows" << std::endl;
    }

    return 0;
  }
  int *pins = new int [ncols * nrows]; 
  int *owners = new int [nrows];

  for (int i=0; i < ncols * nrows; i++){
    pins[i] = 0;
    if (i < nrows) owners[i] = 0;
  }

  int maxNnz = m.MaxNumEntries();
  double *val = new double [maxNnz];
  int *col = new int [maxNnz];

  for (int i=0; i < localRows; i++){

    int rowGID = rowmap.GID(i) - row_base;
    m.ExtractMyRowCopy(i, maxNnz, len, val, col);

    owners[rowGID] = me;

    int *rowptr = pins + (rowGID * ncols);

    for (int j=0; j < len; j++){

      int colGID = colmap.GID(col[j]) - col_base;

      rowptr[colGID] = 1;
    }
  }

  delete [] val;
  delete [] col;

  int *gpins = new int [ncols*nrows];
  int *gowners = new int [nrows];

  Comm.MaxAll(pins, gpins, ncols*nrows);
  Comm.MaxAll(owners, gowners, nrows);

  delete [] pins;
  delete [] owners;

  Comm.Barrier();

  if (me == 0){

    int *pinptr = gpins;
    int xcount=0;

    os << std::endl;
    if (s) os << s << std::endl;
    for (int i=0; i < nrows; i++){
      os << "Proc " << gowners[i] << ", Row " << i+row_base << ": ";
      for (int j=0; j < ncols; j++){
        int pin = *pinptr++;
        if (pin > 0){
          if (!withGraphCuts || (gowners[j] == gowners[i])){
            os << "1 ";
          }
          else{
            os << "X ";
            xcount++;
          }
        }
        else{
          os << "- ";
        }
      }
      os << std::endl;
    }
    os << std::endl;
    if (xcount > 0) os << "Number of cuts: " << xcount << std::endl;
  }

  delete [] gowners;
  delete [] gpins;

  Comm.Barrier();

  return 0;
}

}//namespace ispatest

#endif //HAVE_EPETRA


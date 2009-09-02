//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2009) Sandia Corporation
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

//--------------------------------------------------------------------
// This is a serial program.  It creates a sparse lower triangular
// matrix L, and computes a partitioning of the rows of L into level sets.
// In solving Lx = y for x, the rows in level set k can not be computed
// until the rows in level set k-1 have completed. 
//
// Two optional arguments:
//     --dim=100  --sparsity=.1
//
//  dim - the number of rows in L, default is 50
//  sparsity - the proportion of non-zeros, default is 1/10
//
//--------------------------------------------------------------------


#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_EpetraLevelScheduler.hpp>
#include <Teuchos_RCP.hpp>

#include <string>
#include <map>
#include <list>
#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace std;

#define DEFAULT_DIM 50
#define DEFAULT_SPARSITY .1

#ifdef HAVE_EPETRA

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialComm.h>

Epetra_SerialComm comm;
Epetra_CrsMatrix *L=NULL;
Epetra_Vector *y=NULL;
Epetra_Vector *x=NULL;
Epetra_Vector *lhs=NULL;
Epetra_Vector *rhs=NULL;

vector<int> levelCounts;
map<int, int> forwardPermute;
map<int, int> reversePermute;

void parseArgs(int argc, char *argv[], map<string, float> &argvals);
void make_L(int nrows, float sparsity);
void make_problem();
void compute_level_scheduling();
void print_func(Epetra_CrsMatrix &M, map<int, int> &rowOrder, vector<int> &counts,  bool diag,
                Epetra_Vector &lhs_vec, Epetra_Vector &rhs_vec);

void usage()
{
  cout << "Usage:" << endl;
  cout << "--dim={nrows} --sparsity={proportion of nonzeros}" << endl;
  cout << "Dimension (number of rows) defaults to " << DEFAULT_DIM << endl;
  cout << "Sparsity (proportion of nonzeros) defaults to " << DEFAULT_SPARSITY << endl << endl;
  cout << "--verbose=0    Don't try to print out problem" << endl;
  exit(0);
}

int main(int argc, char *argv[]) 
{
int dim = DEFAULT_DIM;
float sparsity = DEFAULT_SPARSITY;
bool verbose = true;
bool includeDiagonal = true;
map<string, float> args;
vector<int> dummyVector;
Epetra_BlockMap blkmap(1, 1, 0, comm);
Epetra_Vector dummyEpetraVector(blkmap);
map<int, int> dummyPermutation;

  // parse arguments

  parseArgs(argc, argv, args);

  map<string, float>::iterator arg = args.find("help");
  if (arg != args.end()){
    usage();
  }

  arg = args.find("verbose");
  if (arg != args.end()){
    if ((int)(arg->second) == 0)
      verbose = false;
  }

  arg = args.find("dim");
  if (arg != args.end()){
    dim = (int)(arg->second);
    if (dim < 2) usage();
  }

  arg = args.find("sparsity");
  if (arg != args.end()){
    sparsity = arg->second;
    if ((sparsity < 0) || (sparsity > 1.0)) usage();
  }

  // create a sparse lower triangular matrix (write L)

  make_L(dim, sparsity);

  if (verbose){
    cout << "Randomly created sparse lower triangular square matrix:" << endl;
    print_func(*L, dummyPermutation, dummyVector, includeDiagonal,
                   dummyEpetraVector, dummyEpetraVector);
  }

  // compute a level scheduling of L 
  // (write forwardPermute, reversePermute, and levelCounts)

  compute_level_scheduling();

  if (verbose){
    cout << "The partitioning of rows into level sets" << endl;
    print_func(*L, reversePermute , levelCounts, !includeDiagonal,
                   dummyEpetraVector, dummyEpetraVector);
  }

  // Create a problem to be solved. 

  make_problem();

  if (verbose){
    cout << "Lx = y, we'll use different methods to solve for x" << endl;
    print_func(*L, dummyPermutation, dummyVector, includeDiagonal, *lhs, *rhs);
  }

  // Solve L x = y using the levels of L.
  //
  // Let D be the diagonal matrix where d(i,i) = 1/L(i,i).  It scales L
  // so that L has a unit diagonal.
  //
  //  DLx = Dy
  //
  // Subtract off the diagonal of L:
  //
  //  (DL - I)x = Dy - Ix
  //
  // Permute (DL-I) so it is arranged according to the levels:
  //
  //  P(DL - I)x = P(Dy - Ix)   or     L'x = y'
  //
  // The form of L' is such that each level contains a sparse matrix and we
  // can solve the smaller m(i)x(i) = y'(i).

  int nrows = DEFAULT_DIM;

  double *entries;
  int *idx;
  int numEntries;

  Epetra_Vector D(rhs->Map());
  Epetra_Vector I(rhs->Map());

  double value;
  int j;
  for (int i=0; i < nrows; i++){
    L->ExtractMyRowView(i, numEntries, entries, idx);
    for (j=0; j < numEntries; j++){
      if (idx[j] == i) {
        value = entries[j];   // value on L's diagonal
        break;
      }
    }
    if (j == numEntries){
      cerr << "FIX THIS EXAMPLE" << endl;
      exit(1);
    }
    D[i] = 1.0 / value;
    I[i] = 1.0;
  }

  idx = new int [nrows];
  for (int i=0; i < nrows; i++){
    idx[i] = 1;
  }

  Epetra_CrsMatrix P(Copy, L->RowMap(), idx, true);

  map<int, int>::iterator curr = reversePermute.begin();
  map<int, int>::iterator end = reversePermute.end();

  value = 1.0;

  while (curr != end){
    int row = curr->first;   // this row in the permuted matrix
    int col = curr->second;  // gets this row in the original matrix
    P.InsertGlobalValues(row, 1, &value, &col);
    curr++;
  }
  delete [] idx;
  P.FillComplete();

  delete L;
}

void parseArgs(int argc, char *argv[], map<string, float> &argvals)
{
string argString;
list<string::size_type> pos;

  for (int i=1; i < argc; i++){
    argString = argString + string(argv[i]);
  }
  if (argString.size() == 0){
    return;
  }

  // location of arguments in argument string

  string::size_type loc = argString.find("dim");
  if (loc != string::npos)
    pos.push_back(loc);

  loc = argString.find("spar");
  if (loc != string::npos)
    pos.push_back(loc);

  loc = argString.find("help");
  if (loc != string::npos)
    pos.push_back(loc);

  loc = argString.find("verbose");
  if (loc != string::npos)
    pos.push_back(loc);

  pos.push_back(argString.size());

  pos.sort();

  list<string::size_type>::iterator pcurr = pos.begin();
  list<string::size_type>::iterator pend = pos.end();

  // get numeric value for each argument, if it's supposed to have one

  while (1){

    string::size_type argStart = *pcurr;
    pcurr++;
    if (pcurr == pend) break;
    string::size_type argEnd = *pcurr;

    string arg = argString.substr(argStart, argEnd - argStart);

    if (arg.find("help") != string::npos){
      argvals["help"] = 1;
    }
    else{

      string::size_type numStart = arg.find_first_of(".,0123456789");
      string::size_type numEnd = arg.find_last_of(".,0123456789");
  
      string num = arg.substr(numStart, numEnd - numStart + 1);
  
      istringstream is(num);
      float numValue;
      is >> numValue;
  
      if (arg.find("dim") != string::npos){
        argvals["dim"] = numValue;
      }
      else if (arg.find("spar") != string::npos){
        argvals["sparsity"] = numValue;
      }
      else if (arg.find("verbose") != string::npos){
        argvals["verbose"] = numValue;
      }
      else{
        usage();
      }
    }
  }
}
void compute_level_scheduling()
{
  // compute a level scheduling of L 

  Teuchos::RCP<const Epetra_CrsGraph> graph = Teuchos::rcp(&L->Graph());

  Isorropia::Epetra::LevelScheduler level(graph);

  graph.release();

  int numLevels = level.numLevels();

  levelCounts.assign(numLevels, 0);
  int **rows = new int * [numLevels];

  for (int i=0; i < numLevels; i++){
    levelCounts[i] = level.numElemsWithLevel(i);
    rows[i] = new int [levelCounts[i]];
    level.elemsWithLevel(i, rows[i], levelCounts[i]);
  }

  // Create permutation maps.  The forward permutation map maps each existing
  // row to its new row.  The reverse permutation map maps the reverse.

  int nextRow = 0;
  for (int i=0; i < numLevels; i++){
    for (int j=0; j < levelCounts[i]; j++){
      int rowID = rows[i][j];  // a row in level i
      int newRowID = nextRow++;
      forwardPermute[rowID] = newRowID;
      reversePermute[newRowID] = rowID;
    }
  }
}
//----------------------------------------------------------------------
// Create a problem.
//
// Create a simple lhs, calculate rhs.  We'll solve for a lhs later.
//----------------------------------------------------------------------
void make_problem()
{
  int nrows = L->NumGlobalRows();
  Epetra_BlockMap blkmap(nrows, nrows, 1, 0, comm);

  if (lhs){
    delete lhs;
    lhs = NULL;
  }

  lhs = new Epetra_Vector(blkmap);

  double one = 1.0;
  double minus_one = -1.0;

  for (int i=0; i < nrows; i+=2)
    lhs->ReplaceGlobalValues(1, &one, &i);

  for (int i=1; i < nrows; i+=2)
    lhs->ReplaceGlobalValues(1, &minus_one, &i);

  if (rhs){
    delete rhs;
    rhs = NULL;
  }

  rhs = new Epetra_Vector(blkmap);

  L->Multiply(false, *lhs, *rhs);
}
void make_L(int nrows, float sparsity)
{
  int id = getpid();
  srand(id);
  int cutoff = static_cast<int>(static_cast<float>(RAND_MAX) * sparsity);

  char **nonZeros = new char * [nrows];
  int *numNonZeros = new int [nrows];
  for (int i=0; i < nrows; i++){
    nonZeros[i] = new char [i+1];
  }

  int maxNonZeros = 0;
  for (int i=0; i<nrows; i++){
    numNonZeros[i] = 1;
    for (int j=0; j < i; j++){
      int randVal = rand();
      if (randVal < cutoff){
        nonZeros[i][j] = 1;
        numNonZeros[i]++;
      }
      else{
        nonZeros[i][j] = 0;
      }
    }
    nonZeros[i][i] = 1;

    if (numNonZeros[i] > maxNonZeros) maxNonZeros = numNonZeros[i];
  }

  if (L){
    delete L;
    L = NULL;
  }

  Epetra_Map m(nrows, nrows, 0, comm);
  L = new Epetra_CrsMatrix(Copy, m, numNonZeros, true);

  int *row = new int [maxNonZeros];
  double *values = new double [maxNonZeros];
  for (int j=0; j< maxNonZeros; j++){
    values[j] = 1.0;
  }

  for (int i=0; i < nrows; i++){

    int *r = row;

    for (int j=0; j <= i; j++){
      if (static_cast<int>(nonZeros[i][j]) > 0){
        *r = j;
        r++;
      }
    }
    L->InsertGlobalValues(i, numNonZeros[i], values, row);
  }

  for (int i=0; i < nrows; i++){
    delete [] nonZeros[i];
  }
  delete [] nonZeros; 
  delete [] numNonZeros; 
  delete [] row; 
  delete [] values; 

  L->FillComplete();
}

//--------------------------------------------------------------------------
// print out small examples for debugging
//--------------------------------------------------------------------------

void print_func(Epetra_CrsMatrix &M,
                map<int, int> &rowOrder, vector<int> &counts,  bool diag,
                Epetra_Vector &lhs_vec, Epetra_Vector &rhs_vec)
{
  int nrows = M.NumGlobalRows();
  if (nrows > 50){
    cout << "Matrix is too large to print.  It has " << nrows << " rows" << endl;
  }

  int *rowvals = new int [nrows];
  bool print_lhs = false;
  bool print_rhs = false;
  double *vals= NULL;
  int *pos= NULL;
  int *lineBreak = NULL;
  int numNonZeros;
  int nextBreak = -1;
  int nbreaks = counts.size() - 1;

  if (nbreaks > 0){
    lineBreak = new int [nbreaks];
    lineBreak[0] = counts[0] - 1;
    for (int i=1; i < nbreaks; i++){
      lineBreak[i] = lineBreak[i-1] + counts[i];
    }
    nextBreak = 0;
  }

  if (lhs_vec.GlobalLength() == nrows) print_lhs = true;
  if (rhs_vec.GlobalLength() == nrows) print_rhs = true;

  for (int i = 0 ; i < nrows; i++){
    int row = i;

    if (rowOrder.size()){
      row = rowOrder[i];
    }

    M.ExtractMyRowView(row, numNonZeros, vals, pos);
    memset(rowvals, 0, sizeof(int) * nrows);
    for (int j=0; j < numNonZeros; j++){
      rowvals[pos[j]] = vals[j]; 
    }
    if (!diag && (row <= i)){
      rowvals[row] = 0;
    }

    cout << setw(2) << right << row << " ";

    for (int j=0; j <= i; j++){
      // Assuming M's elements are 1 digit - else make precision/width arguments
      cout << setw(2) << setprecision(0) << fixed << right << rowvals[j];
    }

    if (print_lhs || print_rhs){
      int pad = 2 * (nrows - i); 
      for (int k=0; k < pad; k++) cout << " ";
      if (print_lhs) cout << setprecision(3) << setw(12) << scientific << lhs_vec[row];
      else           cout << setw(12) << " ";
      if (print_rhs) cout << setprecision(3) << setw(12) << scientific << rhs_vec[row];
      else           cout << setw(12) << " ";
    }

    cout << endl;

    if (lineBreak && (i == lineBreak[nextBreak])){
      cout << endl;
      if (nextBreak == nbreaks - 1){
        delete [] lineBreak;
        lineBreak = NULL;
      }
      else{
        nextBreak++;
      }
    }
  }
  cout << endl;

  delete [] rowvals;
}
#else
int main( int argc, char *argv[])
{
  cout << "This example requires epetra" << endl;
}
#endif

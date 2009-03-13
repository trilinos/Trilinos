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
// until the rows in level set k-1 have completed.  And the rows in
// level set k can be computed concurrently.
//
// The motivation for computing the level sets is to partition the 
// solution into calculations that can be done concurrently, perhaps
// using multiple cores that share memory.
//
// Two optional arguments:
//     --dim=100  --sparsity=.1
//
//  dim - the number of rows in L, default is 100
//  sparsity - the proportion of non-zeros, default is 1/10
//
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
Epetra_CrsMatrix *L;

void parseArgs(int argc, char *argv[], map<string, float> &argvals);
void make_L(int nrows, float sparsity);
void print_L(map <int, int> &rowOrder, vector<int> &levelSizes, bool printDiag);

void usage()
{
  cout << "Usage:" << endl;
  cout << "--dim={nrows} --sparsity={proportion of nonzeros}" << endl;
  cout << "Dimension (number of rows) defaults to " << DEFAULT_DIM << endl;
  cout << "Sparsity (proportion of nonzeros) defaults to " << DEFAULT_SPARSITY << endl;
  exit(0);
}


int main(int argc, char *argv[]) 
{
int dim = DEFAULT_DIM;
float sparsity = DEFAULT_SPARSITY;
map<string, float> args;
vector<int> dummy;

  // parse arguments

  parseArgs(argc, argv, args);

  map<string, float>::iterator arg = args.find("help");
  if (arg != args.end()){
    usage();
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

  // create a sparse lower triangular matrix

  make_L(dim, sparsity);

  cout << "Initial sparse matrix:" << endl;
  map<int, int> m;
  print_L(m, dummy, true);

  // compute a level scheduling of L

  Teuchos::RCP<const Epetra_CrsGraph> graph = Teuchos::rcp(&L->Graph());

  Isorropia::Epetra::LevelScheduler level(graph);

  graph.release();

  int numLevels = level.numLevels();

  vector<int> levelCounts(numLevels, 0);
  int **rows = new int * [numLevels];

  for (int i=0; i < numLevels; i++){
    levelCounts[i] = level.numElemsWithLevel(i);
    rows[i] = new int [levelCounts[i]];
    level.elemsWithLevel(i, rows[i], levelCounts[i]);
  }

  // Create permutation maps.  The forward permutation map maps each existing
  // row to its new row.  The reverse permutation map maps the reverse.

  map <int, int> forward;
  map <int, int> reverse;
  vector<int> nextRow(numLevels, 0);
  for (int i=1; i < numLevels; i++){
    nextRow[i] = nextRow[i-1] + levelCounts[i-1];
  }

  for (int i=0; i < numLevels; i++){
    for (int j=0; j < levelCounts[i]; j++){
      int rowID = rows[i][j];  // a row in level i
      int newRowID = nextRow[i]++;
      forward[rowID] = newRowID;
      reverse[newRowID] = rowID;
    }
  }

  // Print the rows of L in groups of level sets

  cout << "L reordered into level sets:" << endl;
  print_L(reverse, levelCounts, false);

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
      else{
        usage();
      }
    }
  }
}
void print_L(map<int, int> &rowOrder, vector<int> &counts,  bool diag=true)
{
  int nrows = L->NumGlobalRows();
  if (nrows > 50){
    cout << "Matrix is too large to print.  It has " << nrows << " rows" << endl;
  }

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

  for (int i = 0 ; i < nrows; i++){
    int row = i;

    if (rowOrder.size()){
      row = rowOrder[i];
    }

    L->ExtractMyRowView(row, numNonZeros, vals, pos);

    int nextpos = 0;

    if (row < 10) cout << " " ;
    cout << row << ": ";

    for (int j=0; j <= i; j++){
      if ((nextpos < numNonZeros) && (pos[nextpos] == j)){
        if (diag || (row != j))
          cout << vals[nextpos] << " ";
        else
          cout << "0 ";
        
        nextpos++;
      }
      else{
        cout << "0 ";
      }
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
}
void make_L(int nrows, float sparsity)
{
  int id = getpid();
  srand(id);
  int cutoff = static_cast<int>(static_cast<float>(RAND_MAX) * sparsity);

  Epetra_Map map(nrows, nrows, 0, comm);

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

  L = new Epetra_CrsMatrix(Copy, map, numNonZeros, true);

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
#else
int main( int argc, char *argv[])
{
  cout << "This example requires epetra" << endl;
}
#endif

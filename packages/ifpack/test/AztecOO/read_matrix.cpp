#include <fstream>

#include "Teuchos_CommHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

//I/O for Harwell-Boeing files
#include <iohb.h>

Epetra_CrsMatrix*
read_matrix_mm(const std::string& mm_file,
               const Epetra_Comm& comm)
{
  int my_proc = comm.MyPID();

  int num_global_rows = 0;
  int nnz_per_row = 0;

  std::ifstream* infile = NULL;
  if (my_proc == 0) {
    infile = new std::ifstream(mm_file.c_str());
    if (infile == NULL || !*infile) {
      throw std::runtime_error("Failed to open file "+mm_file);
    }

    std::ifstream& in = *infile;

    //first skip over the file header, which has
    //lines beginning with '%'.
    std::string line;
    do {
      getline(in, line);
    } while(line[0] == '%');

    //now get the matrix dimensions.

    int numrows, numcols, nnz;
    std::istringstream isstr(line);
    isstr >> numrows >> numcols >> nnz;

    //make sure we successfully read the three ints from that line.
    if (isstr.fail()) {
      throw std::runtime_error("Failed to parse matrix-market header.");
    }

    num_global_rows = numrows;
    nnz_per_row = nnz/numrows;
  }

  comm.Broadcast(&num_global_rows, 1, 0);
  comm.Broadcast(&nnz_per_row, 1, 0);

  const int indexBase = 0;
  Epetra_Map rowmap(num_global_rows, indexBase, comm);

  Epetra_CrsMatrix* A = new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  if (my_proc == 0) {
    Teuchos::Array<int> col(1,0);
    Teuchos::Array<double> coef(1,0);

    int irow=0, icol=0;
    double val=0;

    std::string line;
    std::ifstream& in = *infile;
    while(!in.eof()) {
      getline(in, line);
      std::istringstream isstr(line);
      isstr >> irow >> icol >> val;
    
      if (isstr.fail()) continue;

      int g_row = irow-1;
      col[0] = icol-1;
      coef[0] = val;

      A->InsertGlobalValues(g_row, 1, &coef[0], &col[0] );
    }
  }

  A->FillComplete();

  return A;
}

Epetra_MultiVector*
read_vector_mm(const std::string& mm_file,
               const Epetra_Comm& comm)
{
  int my_proc = comm.MyPID();

  int num_global_rows = 0;

  std::ifstream* infile = NULL;
  if (my_proc == 0) {
    infile = new std::ifstream(mm_file.c_str());
    if (infile == NULL || !*infile) {
      throw std::runtime_error("Failed to open file "+mm_file);
    }

    std::ifstream& in = *infile;

    //first skip over the file header, which has
    //lines beginning with '%'.
    std::string line;
    do {
      getline(in, line);
    } while(line[0] == '%');

    //now get the matrix dimensions.

    int numrows, numcols;
    std::istringstream isstr(line);
    isstr >> numrows >> numcols;

    //make sure we successfully read the ints from that line.
    if (isstr.fail()) {
      throw std::runtime_error("Failed to parse matrix-market header.");
    }

    num_global_rows = numrows;
  }

  comm.Broadcast(&num_global_rows, 1, 0);

  const int indexBase = 0;
  Epetra_Map rowmap(num_global_rows, indexBase, comm);

  Epetra_MultiVector* b = new Epetra_MultiVector(rowmap, 1);

  if (my_proc == 0) {
    int irow=0, icol=0;
    double val=0;

    std::string line;
    std::ifstream& in = *infile;
    while(!in.eof()) {
      getline(in, line);
      std::istringstream isstr(line);
      isstr >> val;
    
      if (isstr.fail()) continue;

      b->ReplaceGlobalValue(irow++, icol, val);
    }
  }

  return b;
}


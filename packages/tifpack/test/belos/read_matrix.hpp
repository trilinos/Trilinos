#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <fstream>

#include "Tpetra_CrsMatrix.hpp"

//I/O for Harwell-Boeing files
#include <iohb.h>

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_matrix_mm(const std::string& mm_file,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
 
  std::ifstream infile(mm_file.c_str());
  if (!infile) {
    throw std::runtime_error("Failed to open file "+mm_file);
  }

  //first skip over the file header, which has
  //lines beginning with '%'.
  std::string line;
  do {
    getline(infile, line);
  } while(line[0] == '%');

  //now get the matrix dimensions.

  int numrows, numcols, nnz;
  std::istringstream isstr(line);
  isstr >> numrows >> numcols >> nnz;

  //make sure we successfully read the three ints from that line.
  if (isstr.fail()) {
    throw std::runtime_error("Failed to parse matrix-market header.");
  }

  const GlobalOrdinal num_global_rows = numrows;
  const LocalOrdinal num_local_rows = numrows;//currently hard-coded for 1 proc!!!
  const LocalOrdinal indexBase = 0;
  Teuchos::RCP<const TMap> rowmap = Teuchos::rcp(new TMap(num_global_rows, num_local_rows, indexBase, comm));

  LocalOrdinal nnz_per_row = nnz/numrows;
  Teuchos::RCP<TCRS> A = Teuchos::rcp(new TCRS(rowmap, nnz_per_row));

  Teuchos::Array<GlobalOrdinal> col(1,0);
  Teuchos::Array<Scalar> coef(1,0);

  int irow=0, icol=0;
  double val=0;

  while(!infile.eof()) {
    getline(infile, line);
    std::istringstream isstr(line);
    isstr >> irow >> icol >> val;
    
    GlobalOrdinal g_row = irow-1;
    col[0] = icol-1;
    coef[0] = val;

    A->insertGlobalValues(g_row, col(), coef() );
  }

  A->fillComplete();

  return A;
}

#endif


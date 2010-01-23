#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <iostream>
#include <fstream>

#include "Teuchos_Time.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"

//I/O for Harwell-Boeing files
#include <iohb.h>

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_matrix_mm(const std::string& mm_file,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  Teuchos::Time timer("read_matrix");
  timer.start();

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
 
  int my_proc = comm->getRank();

  GlobalOrdinal num_global_rows = 0;
  LocalOrdinal nnz_per_row = 0;

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

  Teuchos::broadcast<int,GlobalOrdinal>(*comm, (int)0, (int)1, &num_global_rows);
  Teuchos::broadcast<int,LocalOrdinal>(*comm, (int)0, (int)1, &nnz_per_row);

  const LocalOrdinal indexBase = 0;
  Teuchos::RCP<const TMap> rowmap = Teuchos::rcp(new TMap(num_global_rows, indexBase, comm));

  Teuchos::RCP<TCRS> A = Teuchos::rcp(new TCRS(rowmap, nnz_per_row));

  if (my_proc == 0) {
    Teuchos::Array<GlobalOrdinal> col;
    Teuchos::Array<Scalar> coef;

    GlobalOrdinal g_row=0;
    int last_row=-1;
    int irow=0, icol=0;
    double val=0;

    std::string line;
    std::ifstream& in = *infile;
    while(!in.eof()) {
      getline(in, line);
      std::istringstream isstr(line);
      isstr >> irow >> icol >> val;
    
      g_row = irow-1;
      if (g_row != last_row) {
        if (col.size() > 0) {
          A->insertGlobalValues(last_row, col(), coef() );
          col.clear();
          coef.clear();
        }
        last_row = g_row;
      }
      col.push_back(icol-1);
      coef.push_back(val);
    }

    if (col.size() > 0) {
      A->insertGlobalValues(g_row, col(), coef() );
    }

    delete infile;
  }

  A->fillComplete();

  timer.stop();
  if (my_proc==0) {
    std::cout << "proc 0 time to read and fill matrix: " << timer.totalElapsedTime() << std::endl;
  }

  return A;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_vector_mm(const std::string& mm_file,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  Teuchos::Time timer("read_vector");
  timer.start();

  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TMV;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
 
  int my_proc = comm->getRank();

  GlobalOrdinal num_global_rows = 0;

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

    //now get the dimensions.

    int numrows, numcols;
    std::istringstream isstr(line);
    isstr >> numrows >> numcols;

    //make sure we successfully read the ints from that line.
    if (isstr.fail()) {
      throw std::runtime_error("Failed to parse matrix-market header.");
    }

    num_global_rows = numrows;
  }

  Teuchos::broadcast<int,GlobalOrdinal>(*comm, (int)0, (int)1, &num_global_rows);

  const LocalOrdinal indexBase = 0;
  Teuchos::RCP<const TMap> rowmap = Teuchos::rcp(new TMap(num_global_rows, indexBase, comm));

  Teuchos::RCP<TMV> b = Teuchos::rcp(new TMV(rowmap, 1));

  if (my_proc == 0) {
    Teuchos::Array<GlobalOrdinal> col;
    Teuchos::Array<Scalar> coef;

    LocalOrdinal l_row=0;
    double val=0;

    std::string line;
    std::ifstream& in = *infile;
    while(!in.eof()) {
      getline(in, line);
      std::istringstream isstr(line);
      isstr >> val;
      if (isstr.fail()) continue;
    
      Scalar sval = val;
      b->replaceGlobalValue(l_row, 0, sval);
      ++l_row;
    }

    delete infile;
  }

  timer.stop();
  if (my_proc==0) {
    std::cout << "proc 0 time to read and fill vector: " << timer.totalElapsedTime() << std::endl;
  }

  return b;
}

#endif


#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <iostream>
#include <fstream>

#include <Teuchos_Time.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>
// #include <Tpetra_MatrixIO.hpp>

// template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
// Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
// read_matrix_hb(const std::string& hb_file,
//                const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
//                const Teuchos::RCP<Node> &node,
//                int thisNodeWeight)
// {
//   Teuchos::Time timer("read_matrix");
//   timer.start();
// 
//   typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
// 
//   Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A;
//   
//   Teuchos::RCP<const TMap> rowMap;
//   int matrixSize;
//   if (comm->getRank() == 0) {
//     Teuchos::ArrayRCP<char> type;
//     int N, nz, Nrhs;
//     Tpetra::Utils::readHBInfo(hb_file,matrixSize,N,nz,type,Nrhs);
//     // std::cout << matrixSize << " by " << N << ", " << nz << " nonzeros, type==" << type << ", " << Nrhs << " RHS." << std::endl;
//   }
//   Teuchos::broadcast<int,int>(*comm, 0, 1, &matrixSize);
//   rowMap = Tpetra::createWeightedContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(thisNodeWeight,matrixSize,comm,node);
//   // std::cout << *rowMap << std::endl;
//   Tpetra::Utils::readHBMatrix(hb_file,comm,node,A,rowMap);
// 
//   timer.stop();
// 
//   int my_proc = comm->getRank();
// 
//   if (my_proc==0) {
//     std::cout << "proc 0 time to read and fill matrix: " << timer.totalElapsedTime() << std::endl;
//   }
//   Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//   A->describe( *fos, Teuchos::VERB_MEDIUM );
// 
//   return A;
// }

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
void 
read_matrix_mm(const std::string& mm_file,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Node> &node,
               int thisNodeWeight, int maxnumnnz, 
               Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &A,
               Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &dvec)
{
  Teuchos::Time timer("read_matrix");
  timer.start();

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   TCRS;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      TV;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                TMap;
  typedef Teuchos::ScalarTraits<Scalar> SCT;
 
  int my_proc = comm->getRank();

  GlobalOrdinal num_global_rows = 0;
  LocalOrdinal avg_nnz_per_row = 0;

  std::ifstream in(mm_file.c_str());
  if (!in) {
    throw std::runtime_error("Failed to open file "+mm_file);
  }
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
  if (isstr.fail()) {
    throw std::runtime_error("Failed to parse matrix-market header.");
  }
  if (my_proc == 0) {
    std::cout << "Matrix size: " << numrows << " by " << numcols << " with " << nnz << " nonzeros." << std::endl;
  }

  num_global_rows = numrows;
  avg_nnz_per_row = nnz/numrows;
  // make sure we agree, just in case we're parsing different pieces of the same matrix in different files
  Teuchos::broadcast<int,GlobalOrdinal>(*comm, (int)0, (int)1, &num_global_rows);
  Teuchos::broadcast<int,LocalOrdinal>(*comm, (int)0, (int)1, &avg_nnz_per_row);

  Teuchos::RCP<const TMap> rowmap = Tpetra::createWeightedContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(thisNodeWeight,num_global_rows,comm,node);
  if (maxnumnnz > 0) {
    A = Teuchos::rcp(new TCRS(rowmap, maxnumnnz, Tpetra::StaticProfile));
  }
  else {
    A = Teuchos::rcp(new TCRS(rowmap, avg_nnz_per_row, Tpetra::DynamicProfile));
  }

  Teuchos::Array<GlobalOrdinal> col;
  Teuchos::Array<Scalar> coef;

  // this assumes that entries are listed by row, in order

  GlobalOrdinal g_row=0;
  int last_row=-1;
  int irow=0, icol=0;
  double val=0;
  Teuchos::ArrayRCP<Scalar> diags = Teuchos::arcp<Scalar>(rowmap->getNodeNumElements());
  std::fill(diags.begin(), diags.end(), SCT::zero());

  while (!in.eof()) {
    getline(in, line);
    std::istringstream isstr(line);
    isstr >> irow >> icol >> val;
    g_row = irow-1;
    // is it my row?
    if (rowmap->isNodeGlobalElement(g_row)) {
      if (g_row != last_row && col.size() > 0) {  // is it a new row? if so, dump the cache into the matrix.
        // std::cout << "Inserting into row " << last_row << " with maximum colind " << *std::max_element(col.begin(),col.end()) << std::endl;
        A->insertGlobalValues(last_row, col(), coef() );
        col.clear();
        coef.clear();
      }
      last_row = g_row;
      col.push_back(icol-1);
      coef.push_back(val);
      if (g_row == icol-1) {
        diags[rowmap->getLocalElement(g_row)] = SCT::one() / val;
      }
    }
  }

  if (col.size() > 0) { // dump the cache
    // std::cout << "Inserting into row " << g_row << " with maximum colind " << *std::max_element(col.begin(),col.end()) << std::endl;
    A->insertGlobalValues(last_row, col(), coef() );
  }
  // construct Tpetra::Vector from diag values
  dvec = rcp(new TV(rowmap,diags()));
  diags = Teuchos::null;

  in.close();

  A->fillComplete(Tpetra::DoOptimizeStorage);

  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  A->describe( *fos, Teuchos::VERB_MEDIUM );

  timer.stop();
  if (my_proc==0) {
    std::cout << "proc 0 time to read and fill matrix: " << timer.totalElapsedTime() << std::endl;
  }
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_vector_mm(const std::string& mm_file,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               const Teuchos::RCP<Node> &node,
               const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &rowmap)
{
  Teuchos::Time timer("read_vector");
  timer.start();

  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>  TMV;
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
  TEST_FOR_EXCEPTION( rowmap->getGlobalNumElements() != num_global_rows, std::runtime_error, 
      "read_vector_mm: specified row map was not appropriate for this RHS." );

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


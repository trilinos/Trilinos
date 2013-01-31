#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <algorithm>
#include <iostream>
#include <fstream>

#include "Teuchos_Time.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"


template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_matrix_hb(const std::string& hb_file,
               const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
               Teuchos::RCP<Node> node)
{
  Teuchos::Time timer("read_matrix");
  timer.start();

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A;
  Tpetra::Utils::readHBMatrix(hb_file,comm,node,A);

  timer.stop();

  int my_proc = comm->getRank();

  if (my_proc==0) {
    std::cout << "proc 0 time to read and fill matrix: " << timer.totalElapsedTime() << std::endl;
  }

  return A;
}


/// \fn read_matrix_mm
/// \brief Read a sparse matrix from a Matrix Market file
///
/// \param mm_file [in] Path of a Matrix Market file.  To be opened
///   and read only by Process 0 of the given communicator.
/// \param comm [in] Communicator object, over which to distribute the
///   sparse matrix to return.
/// \param node [in] Kokkos Node instance to be used by the returned
///   sparse matrix.
///
/// \return The sparse matrix, distributed over the given communicator.
///
/// \note This defers to Tpetra::MatrixMarket::Reader for reading the
///   sparse matrix from the Matrix Market file.
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_matrix_mm (const std::string& mm_file,
		const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		const Teuchos::RCP<Node>& node)
{
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;
  using std::cout;
  using std::endl;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;

  RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewCounter ("read_matrix");
  RCP<crs_matrix_type> A;
  {
    Teuchos::TimeMonitor timeMon (*timer);
    A = reader_type::readSparseFile (mm_file, comm, node); 
  }

  if (comm->getRank () == 0) {
    cout << "Proc 0: Time in seconds to read the Matrix Market - format sparse "
	 << "matrix and finish fillComplete(): " 
	 << timer->totalElapsedTime () << endl;
  }
  return rcp_const_cast<const crs_matrix_type> (A);
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


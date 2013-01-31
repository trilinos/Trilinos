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


/// \brief Is the given character c not a white space character?
///
/// Is the given character c not a white space character?  Here "white
/// space" means a space or tab.
static bool
notWhiteSpaceCharacter (const char c)
{
  return c != ' ' && c != '\t';
}

/// Is the given line a line of white space (or an empty line)?
///
/// \note This function works if the input line was retrieved using
/// getline(), because getline() discards the end-of-line terminator.
static bool 
isWhiteSpace (const std::string& line)
{
  if (line.size() == 0)
    return true; // An empty line is "white space" too
  else 
    {
      // If the line does _not_ contain a _non_-whitespace character,
      // it is a white space line.  Sorry about the double negative...
      if (std::find_if (line.begin(), line.end(), notWhiteSpaceCharacter) == line.end())
	return true;
      else
	return false;
    }
}

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
///
/// \param comm [in] Communicator object
///
/// \return The sparse matrix, distributed using the communicator
///
/// \warning This is a fragile stopgap implementation.  It only
///   supports real-valued matrices, and it does not check the Matrix
///   Market header for whether the data represents a sparse matrix
///   (Matrix Market format also supports dense matrices).  It also
///   does not support Matrix Market features like symmetric storage
///   or pattern-only sparse matrices.  For example, even if the
///   Matrix Market file claims that it's using symmetric storage,
///   you'll only get the data that's in the file (which could be
///   either the upper or lower triangle).
///
/// \note The fact that Ifpack2 is under development becomes apparent
///   to new users when attempting to input a Matrix Market format
///   matrix.  A design principle of Tpetra and Ifpack2 is to avoid
///   coupling to Epetra and EpetraExt.  Certain Epetra* components
///   will be reimplemented in Tpetra*, a task that has not yet begun.
///   The input files and this parser are coupled.
///
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


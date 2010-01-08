#ifndef TPETRA_MATRIX_IO
#define TPETRA_MATRIX_IO

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace Tpetra {
  namespace Utils {

    bool parseIfmt(Teuchos::ArrayRCP<char> fmt, int &perline, int &width);
    bool parseRfmt(Teuchos::ArrayRCP<char> fmt, int &perline, int &width, int &prec, char &flag);
    void readHBInfo(const std::string &filename, int &M, int &N, int &nz, Teuchos::ArrayRCP<char> &Type, int &Nrhs);

    void readHBHeader(std::ifstream &in_file, Teuchos::ArrayRCP<char> &Title, Teuchos::ArrayRCP<char> &Key, Teuchos::ArrayRCP<char> &Type, 
        int &Nrow, int &Ncol, int &Nnzero, int &Nrhs,
        Teuchos::ArrayRCP<char> &Ptrfmt, Teuchos::ArrayRCP<char> &Indfmt, Teuchos::ArrayRCP<char> &Valfmt, Teuchos::ArrayRCP<char> &Rhsfmt, 
        int &Ptrcrd, int &Indcrd, int &Valcrd, int &Rhscrd, Teuchos::ArrayRCP<char> &Rhstype);

    void readHBMatDouble(const std::string &filename, int &M, int &N, int &nonzeros, std::string &Type, Teuchos::ArrayRCP<int> &colptr, Teuchos::ArrayRCP<int> &rowind, Teuchos::ArrayRCP<double> &val);

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
    void
    readHBMatrix(const std::string &filename, 
                         const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                         const Teuchos::RCP<Node> &node,
                         Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A);
  } // end of Tpetra::Utils namespace
} // end of Tpetra namespace

#ifndef HIDE_TPETRA_INOUT_IMPLEMENTATIONS

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
Tpetra::Utils::readHBMatrix(const std::string &filename, 
                             const Teuchos::RCP<const Teuchos::Comm<int> > &comm, 
                             const Teuchos::RCP<Node> &node,
                             Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > &A)
{
  const int myRank = comm->getRank();
  int numRows,numCols,numNZ;
  Teuchos::ArrayRCP<Scalar> svals;
  Teuchos::ArrayRCP<GlobalOrdinal> colinds;
  Teuchos::ArrayRCP<int>           rowptrs;
  Teuchos::ArrayRCP<size_t>        nnzPerRow;
  Teuchos::ArrayRCP<char>          type;
  int fail = 0;
  if (myRank == 0) {
    bool isSymmetric;
    Teuchos::ArrayRCP<double> dvals;
    Teuchos::ArrayRCP<int> colptrs, rowinds;
    std::string type;
    Tpetra::Utils::readHBMatDouble(filename,numRows,numCols,numNZ,type,colptrs,rowinds,dvals);
    TEST_FOR_EXCEPT(type.size() != 3);
    if (type[0] != 'R' && type[0] != 'r') {
      // only real matrices right now
      fail = 1;
    }
    if (fail == 0 && numNZ > 0) {
      if (type[1] == 'S' || type[1] == 's') {
        isSymmetric = true;
      }
      else {
        isSymmetric = false;
      }
    }
    if (fail == 0 && numNZ > 0) {
      // find num non-zero per row
      nnzPerRow = Teuchos::arcp<size_t>(numRows);
      std::fill(nnzPerRow.begin(), nnzPerRow.end(), 0);
      for (Teuchos::ArrayRCP<int>::const_iterator ri=rowinds.begin(); ri != rowinds.end(); ++ri) {
        // count each row index towards its row
        ++nnzPerRow[*ri-1];
      }
      if (isSymmetric) {
        // count each column toward the corresponding row as well
        for (int c=0; c < numCols; ++c) {
          // the diagonal was already counted; neglect it, if it exists
          for (int i=colptrs[c]-1; i != colptrs[c+1]-1; ++i) {
            if (rowinds[i] != c+1) {
              ++nnzPerRow[c];
              ++numNZ;
            }
          }
        }
      }
      // allocate/set new matrix data
      svals = Teuchos::arcp<Scalar>(numNZ);
      colinds = Teuchos::arcp<GlobalOrdinal>(numNZ);
      rowptrs = Teuchos::arcp<int>(numRows+1);
      rowptrs[0] = 0;
#ifdef HAVE_TPETRA_DEBUG
      Teuchos::ArrayRCP<size_t> nnzPerRow_debug(nnzPerRow.size());
      std::copy(nnzPerRow.begin(), nnzPerRow.end(), nnzPerRow_debug.begin());
#endif
      for (int j=1; j <= numRows; ++j) {
        rowptrs[j] = rowptrs[j-1] + nnzPerRow[j-1];
        nnzPerRow[j-1] = 0;
      }
      // translate from column-oriented to row-oriented
      for (int col=0; col<numCols; ++col) {
        for (int i=colptrs[col]-1; i != colptrs[col+1]-1; ++i) {
          const int row = rowinds[i]-1;
          // add entry to (row,col), with value dvals[i]
          const size_t entry = rowptrs[row] + nnzPerRow[row];
          svals[entry] = Teuchos::as<Scalar>(dvals[i]);
          colinds[entry] = Teuchos::as<GlobalOrdinal>(col);
          ++nnzPerRow[row];
          if (isSymmetric && row != col) {
            // add entry to (col,row), with value dvals[i]
            const size_t symentry = rowptrs[col] + nnzPerRow[col];
            svals[symentry] = Teuchos::as<Scalar>(dvals[i]);
            colinds[symentry] = Teuchos::as<GlobalOrdinal>(row);
            ++nnzPerRow[col];
          }
        }
      }
#ifdef HAVE_TPETRA_DEBUG
      {
        bool isequal = true;
        Teuchos::ArrayRCP<size_t>::const_iterator it1, it2;
        for (it1 = nnzPerRow.begin(), it2 = nnzPerRow_debug.begin(); it1 != nnzPerRow.end(); ++it1, ++it2) {
          if (*it1 != *it2) {
            isequal = false; 
            break;
          }
        }
        TEST_FOR_EXCEPTION(!isequal || nnzPerRow.size() != nnzPerRow_debug.size(), std::logic_error,
            "Tpetra::Utils::readHBMatrix(): Logic error.");
      }
#endif
    }
    // std::cout << "Matrix " << filename << " of type " << type << ": " << numRows << " by " << numCols << ", " << numNZ << " nonzeros" << std::endl;
  }
  // check for read errors
  broadcast(*comm,0,&fail);
  TEST_FOR_EXCEPTION(fail == 1, std::runtime_error, "Tpetra::Utils::readHBMatrix() can only read Real matrices.");
  // distribute global matrix info
  broadcast(*comm,0,&numRows);
  broadcast(*comm,0,&numCols);
  // create map with uniform partitioning
  Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap;
  rowMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((global_size_t)numRows,(GlobalOrdinal)0,comm,GloballyDistributed,node));
  Teuchos::ArrayRCP<size_t> myNNZ;
  if (rowMap->getNodeNumElements()) {
    myNNZ = Teuchos::arcp<size_t>(rowMap->getNodeNumElements());
  }
  if (myRank == 0) {
    size_t numRowsAlreadyDistributed = rowMap->getNodeNumElements();
    std::copy(nnzPerRow.begin(), nnzPerRow.begin()+numRowsAlreadyDistributed,myNNZ);
    for (int p=1; p < Teuchos::size(*comm); ++p) {
      size_t numRowsForP;
      Teuchos::receive(*comm,p,&numRowsForP);
      Teuchos::send<int,size_t>(*comm,numRowsForP,nnzPerRow(numRowsAlreadyDistributed,numRowsForP).getRawPtr(),p);
      numRowsAlreadyDistributed += numRowsForP;
    }
  }
  else {
    const size_t numMyRows = rowMap->getNodeNumElements();
    Teuchos::send(*comm,numMyRows,0);
    Teuchos::receive<int,size_t>(*comm,0,numMyRows,myNNZ(0,numMyRows).getRawPtr());
  }
  nnzPerRow = Teuchos::null;
  // create column map
  Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domMap;
  if (numRows == numCols) {
    domMap = rowMap;
  }
  else {
    domMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>((global_size_t)numCols,(GlobalOrdinal)0,comm,GloballyDistributed,node));
  }
  A = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowMap,myNNZ,Tpetra::StaticProfile));
  // free this locally, A will keep it allocated as long as it is needed by A (up until allocation of nonzeros)
  myNNZ = Teuchos::null;
  if (myRank == 0 && numNZ > 0) {
    for (int r=0; r < numRows; ++r) {
      const size_t nnz = rowptrs[r+1] - rowptrs[r];
      if (nnz > 0) {
        Teuchos::ArrayView<const GlobalOrdinal> inds = colinds(rowptrs[r],nnz);
        Teuchos::ArrayView<const        Scalar> vals = svals(  rowptrs[r],nnz);
        A->insertGlobalValues(r, inds, vals);
      }
    }
  }
  // don't need these anymore
  colinds = Teuchos::null;
  svals   = Teuchos::null;
  rowptrs = Teuchos::null;
  A->fillComplete(domMap,rowMap,Tpetra::DoOptimizeStorage);
}

#endif

#endif

/**
  \file   Amesos2_TpetraCrsMatrixAdapter_def.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Fri May 28 15:04:38 CDT 2010
  
  \brief  Amesos2::MatrixAdapter specialization for the
          Tpetra::CrsMatrix class.
*/

#ifndef AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DEF_HPP

#include <Teuchos_Ptr.hpp>

#include "Amesos2_MatrixAdapter.hpp"

using Tpetra::CrsMatrix;

namespace Amesos {


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::MatrixAdapter()
      : mat_(Teuchos::null)
    { }


/// Copy constructor
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::MatrixAdapter(
    const MatrixAdapter<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >& adapter)
    : mat_(adapter.mat_)
  { }


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::MatrixAdapter(
    const Teuchos::RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >& m)
    : mat_(m)
  { }


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::getCrs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<Tpetra::global_size_t> rowptr,
    size_t& nnz,
    bool local,
    int root)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();

  // Could we have numrows be LocalOrdinal or GlobalOrdinal depending on $local
  GlobalOrdinal numrows = 0;
  if ( local ){
    numrows = mat_->getNodeNumRows();
    nnz = mat_->getNodeNumEntries();
  } else {
    numrows = mat_->getGlobalNumRows();
    nnz = mat_->getGlobalNumEntries();
  }
  TEST_FOR_EXCEPTION( ((rank == root)
      && (Teuchos::as<size_t>(nzval.size()) < nnz)) ||
    (local && (Teuchos::as<size_t>(nzval.size()) < nnz)),
    std::runtime_error,
    "nzval not large enough to hold nonzero values");
  TEST_FOR_EXCEPTION( ((rank == root)
      && (Teuchos::as<size_t>(colind.size()) < nnz)) ||
    (local && (Teuchos::as<size_t>(colind.size()) < nnz)),
    std::runtime_error,
    "colind not large enough to hold index values");
  TEST_FOR_EXCEPTION( ((rank == root)
      && (Teuchos::as<GlobalOrdinal>(rowptr.size()) < numrows + 1)) ||
    (local && (Teuchos::as<GlobalOrdinal>(rowptr.size()) < numrows + 1)),
    std::runtime_error,
    "rowptr not large enough to hold row index values");
  // TODO: Is there possibly a better way to test locality?
  if ( !mat_->getRowMap()->isDistributed() || local ){ // This matrix is *not*
    // distributed
    global_size_type rowInd = 0;
    for( GlobalOrdinal row = 0; row < numrows; ++row ){
      rowptr[row] = rowInd;
      GlobalOrdinal myRow = mat_->getRowMap()->getGlobalElement(row);
      size_t rowNNZ = mat_->getNumEntriesInLocalRow(myRow);
      size_t nnzRet = 0;
      mat_->getGlobalRowCopy(myRow,
        colind.view(rowInd,rowNNZ),
        nzval.view(rowInd,rowNNZ),
        nnzRet);
      TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
        std::runtime_error,
        "Number of values returned different from \
                          number of values reported");
      rowInd += rowNNZ;
    }
    rowptr[numrows] = nnz;
  } else {                      // This matrix is distributed
    using Teuchos::Array;
    // using Teuchos::Ptr;
    using Teuchos::OrdinalTraits;

    /* Each processor has an array the length of globalNumRows, each index
     * initialized to the max ordinal.  For each global index that the
     * processor owns, they write their rank in that spot, however the
     * processor that has rank == root write a -1.  We do a global array MIN
     * reduction on those arrays.  The root writes -1 because we want to
     * assure that we do local operations as much as we can.  It may be
     * possible that two nodes have a copy of the same row, but we want to
     * prefer the root's copy in all cases.
     * 
     * Now the processor with rank == root goes through the globally reduced
     * array: For indices that have -1, then it does a local row copy, for
     * others it gets the row copy from that processor indicated.  Meanwhile,
     * other processors are sending the rows they've been assigned (the global
     * rows corresponding to the indices in the reduced array that == their
     * rank) to the root.
     */

    Array<int> minOwner(numrows); // Will hold global reduction
    Array<int> myRows(numrows,OrdinalTraits<int>::max());
    for( GlobalOrdinal row = 0; row < numrows; ++row ){
      // if global index is found on this node
      if( mat_->getRowMap()->isNodeGlobalElement(row) ){
        myRows[row] = (rank == root) ? -1 : rank;
      }
    }
    // now perform global reduce of these arrays
    Teuchos::MinValueReductionOp<int,int> op;
    Teuchos::reduceAll(*comm, op, (int)myRows.size(),
      myRows.getRawPtr(), minOwner.getRawPtr());

    if( rank == root ){
      global_size_type rowInd = 0;
      size_t rowNNZ, nnzRet;
      for( GlobalOrdinal row = 0; row < numrows; ++row ){
        int from = minOwner[row];
        if( from == -1 ){       // Copy from myself
          rowptr[row] = rowInd;
          rowNNZ = mat_->getNumEntriesInGlobalRow(row);
          nnzRet = 0;
          mat_->getGlobalRowCopy(row,
            colind.view(rowInd,rowNNZ),
            nzval.view(rowInd,rowNNZ),
            nnzRet);
          TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
            std::runtime_error,
            "Number of values returned different from number of values reported");
          rowInd += rowNNZ;
        } else {
          rowptr[row] = rowInd;
          Teuchos::receive(*comm, from, &rowNNZ);
          Teuchos::receive(*comm, from, (int)rowNNZ,
            colind.view(rowInd, rowNNZ).getRawPtr());
          Teuchos::receive(*comm, from, (int)rowNNZ,
            nzval.view(rowInd, rowNNZ).getRawPtr());
          rowInd += rowNNZ;
        }
        rowptr[numrows] = nnz;
      }
    } else {                  // I am not the root
      size_t localNNZ = mat_->getNodeNumEntries();
      Teuchos::Array<GlobalOrdinal> colind_l(colind);
      Teuchos::Array<Scalar> nzval_l(nzval);
      colind_l.resize(localNNZ);
      nzval_l.resize(localNNZ);
      GlobalOrdinal rowInd = 0;
      size_t rowNNZ, nnzRet;
      for( GlobalOrdinal row = 0; row < numrows; ++row ){
        int from = minOwner[row];
        if( from == rank ){   // I am responsible for this row
          rowNNZ = mat_->getNumEntriesInGlobalRow(row);
          nnzRet = 0;
          mat_->getGlobalRowCopy(row,
            colind_l.view(rowInd, rowNNZ),
            nzval_l.view(rowInd, rowNNZ),
            nnzRet);
          TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
            std::runtime_error,
            "Number of values returned different from number of values reported");
          // Send values, column indices, and row pointers to root
          Teuchos::send(*comm, rowNNZ, root);
          Teuchos::send(*comm, (int)rowNNZ,
            colind_l.view(rowInd,rowNNZ).getRawPtr(), root);
          Teuchos::send(*comm, (int)rowNNZ,
            nzval_l.view(rowInd,rowNNZ).getRawPtr(), root);
          rowInd += rowNNZ;   // This is purely local to this node
        }
      }
    }
  }
}


/* We left in the root parameter just for advanced use.  It may in some cases
 * affect the runtime speed of the getCrs method depending on the distribution
 * of the matrix across node.
 */
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::getCrsAll(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<Tpetra::global_size_t> rowptr,
    size_t& nnz,
    int root)
{
  using Teuchos::ArrayView;
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();
  // Root gather data
  if ( rank == root ){
    getCrs(nzval, colind, rowptr, nnz);
  }
  // Then we broadcast to all
  Teuchos::broadcast(*comm, root, (ArrayView<Scalar>)nzval);
  Teuchos::broadcast(*comm, root, (ArrayView<Scalar>)colind);
  Teuchos::broadcast(*comm, root, (ArrayView<Scalar>)rowptr);
  Teuchos::broadcast(*comm, root, Teuchos::ptrFromRef(nnz));
}


/* We left in the root parameter just for advanced use.  It may in some cases
 * affect the runtime speed of the getCrs method depending on the distribution
 * of the matrix across node.
 */
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::getCcsAll(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<Tpetra::global_size_t> colptr,
    size_t& nnz,
    int root)
{
  using Teuchos::ArrayView;
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();
  // Root gather data
  if ( rank == root ){
    getCcs(nzval, rowind, colptr, nnz);
  }
  // Then we broadcast to all
  Teuchos::broadcast(*comm, root, (ArrayView<Scalar>)nzval);
  Teuchos::broadcast(*comm, root, (ArrayView<Scalar>)rowind);
  Teuchos::broadcast(*comm, root, (ArrayView<Scalar>)colptr);
  Teuchos::broadcast(*comm, root, Teuchos::ptrFromRef(nnz));
}


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::getCcs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<Tpetra::global_size_t> colptr,
    size_t& nnz,
    bool local,
    int root)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();

  global_size_type numCols = 0;
  global_size_type numRows = 0;
  
  if ( local ){
    numCols = mat_->getNodeNumCols();
    numRows = mat_->getNodeNumRows();
  } else {
    numCols = mat_->getGlobalNumCols();
    numRows = mat_->getGlobalNumRows();
  }
  
  Array<scalar_type> nzval_tmp(nzval.size());
  Array<GlobalOrdinal> colind(nzval.size());
  Array<global_size_type> rowptr(numRows + 1);

  this->getCrs(nzval_tmp, colind, rowptr, nnz, local, root);
  /* We now have a compressed-row storage format of this matrix.  We transform
   * this into a compressed-column format using a distribution-counting sort
   * algorithm, which is described by D.E. Knuth in TAOCP Vol 3, 2nd ed pg 78.
   
   * Check the conditional execution here.  Also, we were having trouble with
   * Tpetra::CrsMatrix not reporting the correct global number of columns.
   */
  if( local || rank == root ){
    // Count the number of entries in each column
    Array<global_size_type> count(numCols + 1, 0);
    typename Array<GlobalOrdinal>::iterator ind_it;
    for( ind_it = colind.begin(); ind_it != colind.end(); ++ind_it ){
      ++(count[(*ind_it) + 1]);
    }

    // Accumulate
    typename Array<global_size_type>::iterator cnt_it;
    for( cnt_it = count.begin() + 1; cnt_it != count.end(); ++cnt_it ){
      *cnt_it = *cnt_it + *(cnt_it - 1);
    }
    // This becomes the array of column pointers
    colptr.assign(count);

    /* Move the nonzero values into their final place in nzval, based on the
     * counts found previously.
     *
     * This sequence deviates from Knuth's algorithm a bit, following more
     * closely the description presented in Gustavson, Fred G. "Two Fast
     * Algorithms for Sparse Matrices: Multiplication and Permuted
     * Transposition" ACM Trans. Math. Softw. volume 4, number 3, 1978, pages
     * 250--269, http://doi.acm.org/10.1145/355791.355796.
     *
     * The row indices end up in sorted order
     */

    for( global_size_type row = 0; row < numRows; ++row ){
      GlobalOrdinal u = rowptr[row];
      GlobalOrdinal v = rowptr[row + 1];
      for( GlobalOrdinal j = u; j < v; ++j ){
        GlobalOrdinal i = count[colind[j]];
        nzval[i] = nzval_tmp[j];
        rowind[i] = row;
        ++(count[colind[j]]);
      }
    }
  }
}


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::updateValuesCrs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> colind,
    const Teuchos::ArrayView<global_size_type> rowptr)
{
  // TODO: implement!
}
  

template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::updateValuesCcs(
    const Teuchos::ArrayView<Scalar> nzval,
    const Teuchos::ArrayView<GlobalOrdinal> rowind,
    const Teuchos::ArrayView<global_size_type> colptr)
{
  // TODO: implement!
}



/// Get a short description of this adapter class
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
std::string
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::description() const
{
  std::ostringstream oss;
  oss << "Amesos2::MatrixAdapter wrapping: ";
  oss << mat_->description();
  return oss.str();
}


/// Print a description of this adapter to the Fancy Output Stream.
template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
void
MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::describe(
    Teuchos::FancyOStream& os,
    const Teuchos::EVerbosityLevel verbLevel) const
{}


template< typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          class    Node,
          class    LocalMatOps >
const char* MatrixAdapter<
  CrsMatrix<
    Scalar,
    LocalOrdinal,
    GlobalOrdinal,
    Node,
    LocalMatOps
    >
  >::name
= "Amesos2 adapter for Tpetra::CrsMatrix";


} // end namespace Amesos

#endif // AMESOS2_TPETRA_CRSMATRIX_ADAPTER_DEF_HPP

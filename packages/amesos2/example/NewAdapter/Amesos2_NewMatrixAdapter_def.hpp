// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_NewMatrixAdapter_def.hpp
  \author John Doe <jd@sandia.gov>
  \date   today
  
  \brief  Amesos2::MatrixAdapter specialization for the
          NewMatrix class.
*/

#ifndef AMESOS2_NEWMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_NEWMATRIX_ADAPTER_DEF_HPP


namespace Amesos {


MatrixAdapter<NewMatrix>::MatrixAdapter()
  : mat_(Teuchos::null)
{ }


MatrixAdapter<NewMatrix>::MatrixAdapter(
  const MatrixAdapter<NewMatrix>& adapter)
  : mat_(adapter.mat_)
{ }


MatrixAdapter<NewMatrix>::MatrixAdapter(
  const Teuchos::RCP<NewMatrix>& m)
  : mat_(m)
{ }



bool MatrixAdapter<NewMatrix>::isLocal() const
{
  // TODO: implement!
}


const Teuchos::RCP<const Teuchos::Comm<int> >
MatrixAdapter<NewMatrix>::getComm() const
{
  // TODO: implement!
}



size_t MatrixAdapter<NewMatrix>::getLocalNumRows() const
{
  // TODO: implement!
}



size_t MatrixAdapter<NewMatrix>::getLocalNumCols() const
{
  // TODO: implement!
}



global_size_type MatrixAdapter<NewMatrix>::getGlobalNumRows() const
{
  // TODO: implement!
}



global_size_type MatrixAdapter<NewMatrix>::getGlobalNumCols() const
{
  // TODO: implement!
}



size_t MatrixAdapter<NewMatrix>::getLocalNNZ() const
{
  // TODO: implement!
}



global_size_type MatrixAdapter<NewMatrix>::getGlobalNNZ() const
{
  // TODO: implement!
}



size_t MatrixAdapter<NewMatrix>::getMaxNNZ() const
{
  // TODO: implement!
}



Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
MatrixAdapter<NewMatrix>::getRowMap() const
{
  // TODO: implement!
}




Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
MatrixAdapter<NewMatrix>::getColMap() const
{
  // TODO: implement!
}


void
MatrixAdapter<NewMatrix>::getCrs(
  const Teuchos::ArrayView<Scalar> nzval,
  const Teuchos::ArrayView<GlobalOrdinal> colind,
  const Teuchos::ArrayView<global_size_type> rowptr,
  size_t& nnz,
  bool local,
  int root)
{
  // TODO: implement!
}


/*
 * TODO: The following definition assumes that the matrix being adapted lends
 * itself to easier access in the row direction.  If this is not the case, and
 * mat_ is more easily access in column-major fashion, then it may make sense
 * to "swap" this definition to getCrs(), and define getCcs() (using getCrs)
 * instead.
 *
 * This definitions uses getCrs() to first get a CRS representation, then uses
 * an adaption of one of Knuth's algorithms to change it into a CCS format.
 */
void
MatrixAdapter<NewMatrix>::getCcs(
  const Teuchos::ArrayView<Scalar> nzval,
  const Teuchos::ArrayView<GlobalOrdinal> rowind,
  const Teuchos::ArrayView<global_size_type> colptr,
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
    numCols = getLocalNumCols();
    numRows = getLocalNumRows();
  } else {
    numCols = getGlobalNumCols();
    numRows = getGlobalNumRows();
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

/* TODO: The following two definitions should not need to be changed. */

/* We left in the root parameter just for advanced use.  It may in some cases
 * affect the runtime speed of the getCrs method depending on the distribution
 * of the matrix across node.
 */
void
MatrixAdapter<NewMatrix>::getCrsAll(
  const Teuchos::ArrayView<Scalar> nzval,
  const Teuchos::ArrayView<GlobalOrdinal> colind,
  const Teuchos::ArrayView<global_size_type> rowptr,
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
void
MatrixAdapter<NewMatrix>::getCcsAll(
  const Teuchos::ArrayView<Scalar> nzval,
  const Teuchos::ArrayView<GlobalOrdinal> rowind,
  const Teuchos::ArrayView<global_size_type> colptr,
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


void
MatrixAdapter<NewMatrix>::updateValuesCrs(
  const Teuchos::ArrayView<Scalar> nzval,
  const Teuchos::ArrayView<GlobalOrdinal> colind,
  const Teuchos::ArrayView<global_size_type> rowptr)
{
  // TODO: implement!
}
  

void
MatrixAdapter<NewMatrix>::updateValuesCcs(
  const Teuchos::ArrayView<Scalar> nzval,
  const Teuchos::ArrayView<GlobalOrdinal> rowind,
  const Teuchos::ArrayView<global_size_type> colptr)
{
  // TODO: implement!
}


/* TODO: adapt to your needs */
std::string
MatrixAdapter<NewMatrix>::description() const
{
  std::ostringstream oss;
  oss << "Amesos2::MatrixAdapter wrapping: NewMatrix";
  return oss.str();
}


void
MatrixAdapter<NewMatrix>::describe(
  Teuchos::FancyOStream& os,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  // TODO: implement!
}


const char* MatrixAdapter<NewMatrix>::name
= "Amesos2 adapter for NewMatrix";


} // end namespace Amesos

#endif // AMESOS2_NEWMATRIX_ADAPTER_DEF_HPP

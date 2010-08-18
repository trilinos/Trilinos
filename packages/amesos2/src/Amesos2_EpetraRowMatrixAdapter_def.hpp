// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_EpetraRowMatrixAdapter_def.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Tue Jul 20 00:05:08 CDT 2010
  
  \brief  Amesos2::MatrixAdapter specialization for the
          Epetra_RowMatrix interface.
*/

#ifndef AMESOS2_EPETRA_ROWMATRIX_ADAPTER_DEF_HPP
#define AMESOS2_EPETRA_ROWMATRIX_ADAPTER_DEF_HPP


namespace Amesos {

using Teuchos::RCP;
using Teuchos::rcp;


MatrixAdapter<Epetra_RowMatrix>::MatrixAdapter()
  : mat_(Teuchos::null)
{ }


MatrixAdapter<Epetra_RowMatrix>::MatrixAdapter(
  const MatrixAdapter<Epetra_RowMatrix>& adapter)
  : mat_(adapter.mat_)
{ }


MatrixAdapter<Epetra_RowMatrix>::MatrixAdapter(
  const RCP<Epetra_RowMatrix>& m)
  : mat_(m)
{ }


bool MatrixAdapter<Epetra_RowMatrix>::isLocal() const
{
  return mat_->Map().DistributedGlobal();
}


size_t MatrixAdapter<Epetra_RowMatrix>::getLocalNumRows() const
{
  return Teuchos::as<size_t>(mat_->NumMyRows());
}


size_t MatrixAdapter<Epetra_RowMatrix>::getLocalNumCols() const
{
  return Teuchos::as<size_t>(mat_->NumMyCols());
}


MatrixAdapter<Epetra_RowMatrix>::global_size_type
MatrixAdapter<Epetra_RowMatrix>::getGlobalNumRows() const
{
  return Teuchos::as<global_size_type>(mat_->NumGlobalRows());
}


MatrixAdapter<Epetra_RowMatrix>::global_size_type
MatrixAdapter<Epetra_RowMatrix>::getGlobalNumCols() const
{
  return Teuchos::as<global_size_type>(mat_->NumGlobalCols());
}


size_t MatrixAdapter<Epetra_RowMatrix>::getLocalNNZ() const
{
  return Teuchos::as<size_t>(mat_->NumMyNonzeros());
}


MatrixAdapter<Epetra_RowMatrix>::global_size_type
MatrixAdapter<Epetra_RowMatrix>::getGlobalNNZ() const
{
  return Teuchos::as<global_size_type>(mat_->NumGlobalNonzeros());
}


size_t MatrixAdapter<Epetra_RowMatrix>::getMaxNNZ() const
{
  return Teuchos::as<size_t>(mat_->MaxNumEntries());
}


const RCP<const Teuchos::Comm<int> >
MatrixAdapter<Epetra_RowMatrix>::getComm() const
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(rcp(mat_->Comm().Clone()));
  if( serialEpetraComm.get() ) {
    RCP<const Teuchos::SerialComm<int> >
      serialComm = rcp(new Teuchos::SerialComm<int>());
    set_extra_data( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }

#ifdef HAVE_MPI
  
  RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm>(rcp(mat_->Comm().Clone()));
  if( mpiEpetraComm.get() ) {
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
  set_extra_data( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg(rawMpiComm) );
  RCP<const Teuchos::MpiComm<int> >
    mpiComm = rcp(new Teuchos::MpiComm<int>(rawMpiComm));
  return mpiComm; 
}

#endif // HAVE_MPI
  
  // If you get here then the conversion failed!
  return Teuchos::null;
}


RCP<const Tpetra::Map<MatrixAdapter<Epetra_RowMatrix>::local_ordinal_type,
                      MatrixAdapter<Epetra_RowMatrix>::global_ordinal_type,
                      MatrixAdapter<Epetra_RowMatrix>::node_type> >
MatrixAdapter<Epetra_RowMatrix>::getRowMap() const
{
  Teuchos::Array<int> local_element_list(mat_->RowMatrixRowMap().NumMyElements());
  // Get element list with mat_->MyGlobalElements(), then use that list, the
  // number of global elements, and the index base to create a Tpetra::Map
  mat_->RowMatrixRowMap().MyGlobalElements( local_element_list.getRawPtr() );
  int num_global_elements = mat_->RowMatrixRowMap().NumGlobalElements();
  int index_base = mat_->RowMatrixRowMap().IndexBase();

  return( rcp(new Tpetra::Map<local_ordinal_type,global_ordinal_type>(
        num_global_elements,
        local_element_list,
        index_base,
        this->getComm())) );
}




RCP<const Tpetra::Map<MatrixAdapter<Epetra_RowMatrix>::local_ordinal_type,
                      MatrixAdapter<Epetra_RowMatrix>::global_ordinal_type,
                      MatrixAdapter<Epetra_RowMatrix>::node_type> >
MatrixAdapter<Epetra_RowMatrix>::getColMap() const
{
  Teuchos::Array<int> local_element_list(mat_->RowMatrixColMap().NumMyElements());
  // Get element list with mat_->MyGlobalElements(), then use that list, the
  // number of global elements, and the index base to create a Tpetra::Map
  mat_->RowMatrixColMap().MyGlobalElements( local_element_list.getRawPtr() );
  int num_global_elements = mat_->RowMatrixColMap().NumGlobalElements();
  int index_base = mat_->RowMatrixColMap().IndexBase();

  return( rcp(new Tpetra::Map<local_ordinal_type,global_ordinal_type>(
        num_global_elements,
        local_element_list,
        index_base,
        this->getComm())) );
}


void
MatrixAdapter<Epetra_RowMatrix>::getCrs(
  const Teuchos::ArrayView<scalar_type> nzval,
  const Teuchos::ArrayView<global_ordinal_type> colind,
  const Teuchos::ArrayView<global_size_type> rowptr,
  size_t& nnz,
  bool local)
{
  RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();
  const int root = 0;

  // Could we have numrows be local_ordinal_type or global_ordinal_type
  // depending on $local
  global_ordinal_type numrows = 0;
  if ( local ){
    numrows = mat_->NumMyRows();
    nnz = mat_->NumMyNonzeros();
  } else {
    numrows = mat_->NumGlobalRows();
    nnz = mat_->NumGlobalNonzeros();
  }
  TEST_FOR_EXCEPTION( ((rank == root)
      && (Teuchos::as<size_t>(nzval.size()) < nnz)) ||
    (local && (Teuchos::as<size_t>(nzval.size()) < nnz)),
    std::length_error,
    "nzval not large enough to hold nonzero values");
  TEST_FOR_EXCEPTION( ((rank == root)
      && (Teuchos::as<size_t>(colind.size()) < nnz)) ||
    (local && (Teuchos::as<size_t>(colind.size()) < nnz)),
    std::length_error,
    "colind not large enough to hold index values");
  TEST_FOR_EXCEPTION( ((rank == root)
      && (Teuchos::as<global_ordinal_type>(rowptr.size()) < numrows + 1)) ||
    (local && (Teuchos::as<global_ordinal_type>(rowptr.size()) < numrows + 1)),
    std::length_error,
    "rowptr not large enough to hold row index values");
  // TODO: Is there possibly a better way to test locality?
  if ( !mat_->RowMatrixRowMap().DistributedGlobal() || local ){ // This matrix is *not*
    // distributed
    global_size_type rowInd = 0;
    for( global_ordinal_type row = 0; row < numrows; ++row ){
      rowptr[row] = rowInd;
      int rowNNZ;
      mat_->NumMyRowEntries(row,rowNNZ);
      int nnzRet = 0;
      TEST_FOR_EXCEPTION(
        mat_->ExtractMyRowCopy(
          row,
          rowNNZ,
          nnzRet,
          nzval.view(rowInd,rowNNZ).getRawPtr(),
          colind.view(rowInd,rowNNZ).getRawPtr()) != 0, // Method returns 0 on success
        std::runtime_error,
        "Error extracting row copy from Epetra_RowMatrix");
      rowInd += rowNNZ;
    }
    rowptr[numrows] = nnz;
  } else {                      // This matrix is distributed
    using Teuchos::Array;
    using Teuchos::OrdinalTraits;

    /* Each processor has an array the length of globalNumRows, each index
     * initialized to the max ordinal.  For each global index that the
     * processor owns, they write their rank in that spot, however the
     * processor that has rank == root writes a -1.  We do a global array MIN
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
    for( global_ordinal_type row = 0; row < numrows; ++row ){
      // if global index is found on this node
      if( mat_->RowMatrixRowMap().MyGID(row) ){
        myRows[row] = (rank == root) ? -1 : rank;
      }
    }
    // now perform global reduce of these arrays
    Teuchos::MinValueReductionOp<int,int> op;
    Teuchos::reduceAll(*comm, op, (int)myRows.size(),
      myRows.getRawPtr(), minOwner.getRawPtr());

    if( rank == root ){
      global_size_type rowInd = 0;
      int rowNNZ, nnzRet;
      for( global_ordinal_type row = 0; row < numrows; ++row ){
        int from = minOwner[row];
        if( from == -1 ){       // Copy from myself
          rowptr[row] = rowInd;
          mat_->NumMyRowEntries(row,rowNNZ);
          nnzRet = 0;
          TEST_FOR_EXCEPTION(
            mat_->ExtractMyRowCopy(
              row,
              rowNNZ,
              nnzRet,
              nzval.view(rowInd,rowNNZ).getRawPtr(),
              colind.view(rowInd,rowNNZ).getRawPtr()) != 0, // Method returns 0 on success
            std::runtime_error,
            "Error extracting row copy from Epetra_RowMatrix");
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
      int localNNZ = mat_->NumMyNonzeros();
      Teuchos::Array<global_ordinal_type> colind_l(colind);
      Teuchos::Array<scalar_type> nzval_l(nzval);
      colind_l.resize(localNNZ);
      nzval_l.resize(localNNZ);
      global_ordinal_type rowInd = 0;
      int rowNNZ, nnzRet;
      for( global_ordinal_type row = 0; row < numrows; ++row ){
        int from = minOwner[row];
        if( from == rank ){   // I am responsible for this row
          int my_row = mat_->RowMatrixRowMap().LID(row);
          mat_->NumMyRowEntries(my_row,rowNNZ);
          nnzRet = 0;
          TEST_FOR_EXCEPTION(
            mat_->ExtractMyRowCopy(
              my_row,
              rowNNZ,
              nnzRet,
              nzval_l.view(rowInd, rowNNZ).getRawPtr(),
              colind_l.view(rowInd, rowNNZ).getRawPtr()) != 0,
            std::runtime_error,
            "Error extraction row copy from Epetra_RowMatrix");
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


/*
 * This definitions uses getCrs() to first get a CRS representation, then uses
 * an adaption of one of Knuth's algorithms to change it into a CCS format.
 */
void
MatrixAdapter<Epetra_RowMatrix>::getCcs(
  const Teuchos::ArrayView<scalar_type> nzval,
  const Teuchos::ArrayView<global_ordinal_type> rowind,
  const Teuchos::ArrayView<global_size_type> colptr,
  size_t& nnz,
  bool local)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;

  RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();
  const int root = 0;

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
  Array<global_ordinal_type> colind(nzval.size());
  Array<global_size_type> rowptr(numRows + 1);

  getCrs(nzval_tmp, colind, rowptr, nnz, local);
  /* We now have a compressed-row storage format of this matrix.  We transform
   * this into a compressed-column format using a distribution-counting sort
   * algorithm, which is described by D.E. Knuth in TAOCP Vol 3, 2nd ed pg 78.
   
   * Check the conditional execution here.  Also, we were having trouble with
   * Tpetra::CrsMatrix not reporting the correct global number of columns.
   */
  if( local || rank == root ){
    // Count the number of entries in each column
    Array<global_size_type> count(numCols + 1, 0);
    Array<global_ordinal_type>::iterator ind_it;
    for( ind_it = colind.begin(); ind_it != colind.end(); ++ind_it ){
      ++(count[(*ind_it) + 1]);
    }

    // Accumulate
    Array<global_size_type>::iterator cnt_it;
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
      global_ordinal_type u = rowptr[row];
      global_ordinal_type v = rowptr[row + 1];
      for( global_ordinal_type j = u; j < v; ++j ){
        global_ordinal_type i = count[colind[j]];
        nzval[i] = nzval_tmp[j];
        rowind[i] = row;
        ++(count[colind[j]]);
      }
    }
  }
}


/* We left in the root parameter just for advanced use.  It may in some cases
 * affect the runtime speed of the getCrs method depending on the distribution
 * of the matrix across node.
 */
void
MatrixAdapter<Epetra_RowMatrix>::getCrsAll(
  const Teuchos::ArrayView<scalar_type> nzval,
  const Teuchos::ArrayView<global_ordinal_type> colind,
  const Teuchos::ArrayView<global_size_type> rowptr,
  size_t& nnz)
{
  using Teuchos::ArrayView;
  RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();
  const int root = 0;

  // Root gather data
  if ( rank == root ){
    getCrs(nzval, colind, rowptr, nnz);
  }
  // Then we broadcast to all
  Teuchos::broadcast(*comm, root, (ArrayView<scalar_type>)nzval);
  Teuchos::broadcast(*comm, root, (ArrayView<global_ordinal_type>)colind);
  Teuchos::broadcast(*comm, root, (ArrayView<global_size_type>)rowptr);
  Teuchos::broadcast(*comm, root, Teuchos::ptrFromRef(nnz));
}


/* We left in the root parameter just for advanced use.  It may in some cases
 * affect the runtime speed of the getCrs method depending on the distribution
 * of the matrix across node.
 */
void
MatrixAdapter<Epetra_RowMatrix>::getCcsAll(
  const Teuchos::ArrayView<scalar_type> nzval,
  const Teuchos::ArrayView<global_ordinal_type> rowind,
  const Teuchos::ArrayView<global_size_type> colptr,
  size_t& nnz)
{
  using Teuchos::ArrayView;
  RCP<const Teuchos::Comm<int> > comm = getComm();
  const int rank = comm->getRank();
  const int root = 0;

  // Root gather data
  if ( rank == root ){
    getCcs(nzval, rowind, colptr, nnz);
  }
  // Then we broadcast to all
  Teuchos::broadcast(*comm, root, (ArrayView<scalar_type>)nzval);
  Teuchos::broadcast(*comm, root, (ArrayView<global_ordinal_type>)rowind);
  Teuchos::broadcast(*comm, root, (ArrayView<global_size_type>)colptr);
  Teuchos::broadcast(*comm, root, Teuchos::ptrFromRef(nnz));
}


void
MatrixAdapter<Epetra_RowMatrix>::updateValuesCrs(
  const Teuchos::ArrayView<scalar_type> nzval,
  const Teuchos::ArrayView<global_ordinal_type> colind,
  const Teuchos::ArrayView<global_size_type> rowptr)
{
  // TODO: implement!
}
  

void
MatrixAdapter<Epetra_RowMatrix>::updateValuesCcs(
  const Teuchos::ArrayView<scalar_type> nzval,
  const Teuchos::ArrayView<global_ordinal_type> rowind,
  const Teuchos::ArrayView<global_size_type> colptr)
{
  // TODO: implement!
}


std::string
MatrixAdapter<Epetra_RowMatrix>::description() const
{
  std::ostringstream oss;
  oss << "Amesos2::MatrixAdapter wrapping: ";
  oss << mat_->Label();
  return oss.str();
}


void
MatrixAdapter<Epetra_RowMatrix>::describe(
  Teuchos::FancyOStream& os,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  // TODO: implement!
}


const char* MatrixAdapter<Epetra_RowMatrix>::name
= "Amesos2 adapter for Epetra_RowMatrix";


} // end namespace Amesos

#endif // AMESOS2_EPETRA_ROWMATRIX_ADAPTER_DEF_HPP

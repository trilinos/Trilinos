//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef TPETRA_MATRIXMATRIX_DEF_HPP
#define TPETRA_MATRIXMATRIX_DEF_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_map.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MMHelpers_def.hpp"

#ifdef DOXYGEN_USE_ONLY
  //#include "Tpetra_MMMultiply_decl.hpp"
#endif

/*! \file Tpetra_MMMultiply_def.hpp 

    The implementations for the members of class Tpetra::MatrixMatrixMultiply and related non-member constructors.
 */
namespace Tpetra {

template<class Scalar, class LocalOrdinal>
Scalar sparsedot(Teuchos::ArrayRCP<Scalar> u, Teuchos::ArrayRCP<LocalOrdinal> u_ind, 
     Teuchos::ArrayRCP<Scalar> v, Teuchos::ArrayRCP<LocalOrdinal> v_ind)
{
  Scalar result = Teuchos::ScalarTraits<Scalar>::zero();

  LocalOrdinal v_idx = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal u_idx = Teuchos::OrdinalTraits<LocalOrdinal>::zero();

  while(v_idx < v.size() && u_idx < u.size()) {
    LocalOrdinal ui = u_ind[u_idx];
    LocalOrdinal vi = v_ind[v_idx];

    if (ui < vi) {
      ++u_idx;
    }
    else if (ui > vi) {
      ++v_idx;
    }
    else {
      result += u[u_idx++]*v[v_idx++];
    }
  }

  return(result);
}

//kernel method for computing the local portion of C = A*B

template<class Scalar, 
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatVec, 
  class SpMatSlv>
int mult_A_B(
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Aview, 
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C)
{

  LocalOrdinal C_firstCol = Bview->colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Bview->colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

  //int* bcols = Bview.colMap->MyGlobalElements();
  Teuchos::ArrayView<const GlobalOrdinal> bcols =Bview->colMap->getNodeElementList();
  //int* bcols_import = NULL;
  Teuchos::ArrayView<const GlobalOrdinal> bcols_import = Teuchos::null;
  if (Bview->importColMap != Teuchos::null) {
    C_firstCol_import = Bview->importColMap->getMinLocalIndex();
    C_lastCol_import = Bview->importColMap->getMaxLocalIndex();

    bcols_import = Bview->importColMap->getNodeElementList();
  }

  size_t C_numCols = C_lastCol - C_firstCol + Teuchos::OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + Teuchos::OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;
  Teuchos::ArrayRCP<Scalar> dwork = Teuchos::ArrayRCP<Scalar>(C_numCols);
  //int* iwork = new int[C_numCols];
  Teuchos::ArrayRCP<GlobalOrdinal> iwork = Teuchos::ArrayRCP<GlobalOrdinal>(C_numCols);

  Teuchos::ArrayRCP<Scalar> C_row_i = dwork;
  //int* C_cols = iwork;
  Teuchos::ArrayRCP<GlobalOrdinal> C_cols = iwork;

  //int C_row_i_length, i, j, k;
  size_t C_row_i_length, i, j, k;

  //To form C = A*B we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc.

  bool C_filled = C.isFillComplete();

  //loop over the rows of A.
  for(i=0; i<Aview->numRows; ++i) {

    //only navigate the local portion of Aview... (It's probable that we
    //imported more of A than we need for A*B, because other cases like A^T*B 
    //need the extra rows.)
    if (Aview->remote[i]) {
      continue;
    }

    Teuchos::ArrayRCP<const LocalOrdinal> Aindices_i = Aview->indices[i];
    Teuchos::ArrayRCP<const Scalar> Aval_i  = Aview->values[i];

    GlobalOrdinal global_row = Aview->rowMap->getGlobalElement(i);


    //loop across the i-th row of A and for each corresponding row
    //in B, loop across colums and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i. In other words,
    //as we stride across B(k,:) we're calculating updates for row i of the
    //result matrix C.

    for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<Aview->numEntriesPerRow[i]; ++k) {
      LocalOrdinal Ak = Bview->rowMap->getLocalElement(Aview->colMap->getGlobalElement(Aindices_i[k]));
      Scalar Aval = Aval_i[k];

      Teuchos::ArrayRCP<const LocalOrdinal> Bcol_inds = Bview->indices[Ak];
      Teuchos::ArrayRCP<const Scalar> Bvals_k = Bview->values[Ak];

      C_row_i_length = Teuchos::OrdinalTraits<size_t>::zero();

      if (Bview->remote[Ak]) {
        for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Bview->numEntriesPerRow[Ak]; ++j) {
          C_row_i[C_row_i_length] = Aval*Bvals_k[j];
          C_cols[C_row_i_length++] = bcols_import[Bcol_inds[j]];
        }
      }
      else {
        for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Bview->numEntriesPerRow[Ak]; ++j) {
          C_row_i[C_row_i_length] = Aval*Bvals_k[j];
          C_cols[C_row_i_length++] = bcols[Bcol_inds[j]];
        }
      }
	  /*
	  std::cout << "About to insert row: " << global_row << std::endl;
	  Teuchos::ArrayView<const Scalar> C_row_iView = C_row_i.view(Teuchos::OrdinalTraits<size_t>::zero(), C_row_i_length);
      typename Teuchos::ArrayView<const Scalar>::iterator it = C_row_iView.begin();
      for(; it != C_row_iView.end(); ++it){
        std::cout << *it << ", ";
      }
      std::cout << std::endl;*/

      //
      //Now put the C_row_i values into C.
      //

      C_filled ?
        C.sumIntoGlobalValues(global_row, C_cols.view(Teuchos::OrdinalTraits<size_t>::zero(), C_row_i_length), C_row_i.view(Teuchos::OrdinalTraits<size_t>::zero(), C_row_i_length))
        :
        C.insertGlobalValues(global_row, C_cols.view(Teuchos::OrdinalTraits<size_t>::zero(), C_row_i_length), C_row_i.view(Teuchos::OrdinalTraits<size_t>::zero(), C_row_i_length));

    }
  }

  //delete [] dwork;
  //delete [] iwork;

  return(0);
}

//kernel method for computing the local portion of C = A*B^T
template<
  class Scalar,
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node,
  class SpMatVec, 
  class SpMatSlv>
int mult_A_Btrans(
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Aview, 
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C)
{
  size_t i, j, k;
  int returnValue = 0;

  size_t maxlen = Teuchos::OrdinalTraits<size_t>::zero();
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Aview->numRows; ++i) {
    if (Aview->numEntriesPerRow[i] > maxlen) maxlen = Aview->numEntriesPerRow[i];
  }
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Bview->numRows; ++i) {
    if (Bview->numEntriesPerRow[i] > maxlen) maxlen = Bview->numEntriesPerRow[i];
  }

  //cout << "Aview: " << endl;
  //dumpCrsMatrixStruct(Aview);

  //cout << "Bview: " << endl;
  //dumpCrsMatrixStruct(Bview);

  //int numBcols = Bview.colMap->NumMyElements();
  size_t numBcols = Bview->colMap->getNodeNumElements();
  size_t numBrows = Bview->numRows;

  size_t iworklen = maxlen*2 + numBcols;
  //int* iwork = new int[iworklen];
  Teuchos::ArrayRCP<LocalOrdinal> iwork = Teuchos::ArrayRCP<LocalOrdinal>(iworklen);

  //int* bcols = iwork+maxlen*2;
  Teuchos::ArrayRCP<GlobalOrdinal> bcols = iwork+maxlen*2;
  //int* bgids = Bview.colMap->MyGlobalElements();
  Teuchos::ArrayView<const GlobalOrdinal> bgids = Bview->colMap->getNodeElementList();
  Teuchos::ArrayRCP<Scalar> bvals = Teuchos::ArrayRCP<Scalar>(maxlen*2);
  Teuchos::ArrayRCP<Scalar> avals = bvals+maxlen;

  GlobalOrdinal max_all_b = Bview->colMap->getMaxAllGlobalIndex();
  GlobalOrdinal min_all_b = Bview->colMap->getMinAllGlobalIndex();

  //bcols will hold the GIDs from B's column-map for fast accMinAllGIDess
  //during the computations below
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<numBcols; ++i) {
    LocalOrdinal blid = Bview->colMap->getLocalElement(bgids[i]);
    bcols[blid] = bgids[i];
  }

  //next create arrays indicating the first and last column-index in
  //each row of B, so that we can know when to skip certain rows below.
  //This will provide a large performance gain for banded matrices, and
  //a somewhat smaller gain for *most* other matrices.
  //int* b_firstcol = new int[2*numBrows];
  Teuchos::ArrayRCP<GlobalOrdinal> b_firstcol = Teuchos::ArrayRCP<GlobalOrdinal>(2*numBrows);
  //int* b_lastcol = b_firstcol+numBrows;
  Teuchos::ArrayRCP<GlobalOrdinal> b_lastcol = b_firstcol+numBrows;
  //int temp;
  GlobalOrdinal temp;
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<numBrows; ++i) {
    b_firstcol[i] = max_all_b;
    b_lastcol[i] = min_all_b;

    size_t Blen_i = Bview->numEntriesPerRow[i];
    if (Blen_i < Teuchos::OrdinalTraits<size_t>::one()) continue;
    Teuchos::ArrayRCP<const LocalOrdinal> Bindices_i = Bview->indices[i];

    if (Bview->remote[i]) {
      for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<Blen_i; ++k) {
        temp = Bview->importColMap->getGlobalElement(Bindices_i[k]);
        if (temp < b_firstcol[i]) b_firstcol[i] = temp;
        if (temp > b_lastcol[i]) b_lastcol[i] = temp;
      }
    }
    else {
      for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<Blen_i; ++k) {
        temp = bcols[Bindices_i[k]];
        if (temp < b_firstcol[i]) b_firstcol[i] = temp;
        if (temp > b_lastcol[i]) b_lastcol[i] = temp;
      }
    }
  }

  //Epetra_Util util;

  //int* Aind = iwork;
  Teuchos::ArrayRCP<GlobalOrdinal> Aind = iwork;
  //int* Bind = iwork+maxlen;
  Teuchos::ArrayRCP<GlobalOrdinal> Bind = iwork+maxlen;

  bool C_filled = C.isFillComplete();

  //To form C = A*B^T, we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(j,k) )
  //
  //This is the easiest case of all to code (easier than A*B, A^T*B, A^T*B^T).
  //But it requires the use of a 'sparsedot' function (we're simply forming
  //dot-products with row A_i and row B_j for all i and j).

  //loop over the rows of A.
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Aview->numRows; ++i) {
    if (Aview->remote[i]) {
      continue;
    }

    Teuchos::ArrayRCP<const LocalOrdinal> Aindices_i = Aview->indices[i];
    Teuchos::ArrayRCP<const Scalar> Aval_i  = Aview->values[i];
    size_t A_len_i = Aview->numEntriesPerRow[i];
    if (A_len_i < Teuchos::OrdinalTraits<size_t>::one()) {
      continue;
    }

    for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<A_len_i; ++k) {
      Aind[k] = Aview->colMap->getGlobalElement(Aindices_i[k]);
      avals[k] = Aval_i[k];
    }

    //util.Sort(true, A_len_i, Aind, 1, &avals, 0, NULL);
    sort2(Aind.begin(), Aind.end(), avals.begin());

    //int mina = Aind[0];
    GlobalOrdinal mina = Aind[0];
  //int maxa = Aind[A_len_i-1];
    GlobalOrdinal maxa = Aind[A_len_i-1];

    if (mina > max_all_b || maxa < min_all_b) {
      continue;
    }

    GlobalOrdinal global_row = Aview->rowMap->getGlobalElement(i);

    //loop over the rows of B and form results C_ij = dot(A(i,:),B(j,:))
    for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Bview->numRows; ++j) {
      if (b_firstcol[j] > maxa || b_lastcol[j] < mina) {
        continue;
      }

      Teuchos::ArrayRCP<const LocalOrdinal> Bindices_j = Bview->indices[j];
      size_t B_len_j = Bview->numEntriesPerRow[j];
      if (B_len_j < Teuchos::OrdinalTraits<size_t>::one()) {
        continue;
      }

      //int tmp, Blen = 0;
      GlobalOrdinal tmp;
      size_t Blen = Teuchos::OrdinalTraits<size_t>::zero();



      if (Bview->remote[j]) {
        for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<B_len_j; ++k) {
          tmp = Bview->importColMap->getGlobalElement(Bindices_j[k]);
          if (tmp < mina || tmp > maxa) {
            continue;
          }

          bvals[Blen] = Bview->values[j][k];
          Bind[Blen++] = tmp;
        }
      }
      else {
        for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<B_len_j; ++k) {
          tmp = bcols[Bindices_j[k]];
          if (tmp < mina || tmp > maxa) {
            continue;
          }

          bvals[Blen] = Bview->values[j][k];
          Bind[Blen++] = tmp;
        }
      }

      if (Blen < Teuchos::OrdinalTraits<size_t>::one()) {
        continue;
      }

      //util.Sort(true, Blen, Bind, 1, &bvals, 0, NULL);
      sort2(Bind.begin(), Bind.end(), bvals.begin());

      Teuchos::ArrayRCP<Scalar> C_ij = Teuchos::ArrayRCP<Scalar>(1, sparsedot(avals, Aind, bvals, Bind ));

      if (C_ij[0] == Teuchos::ScalarTraits<Scalar>::zero()) {
        continue;
      }
      Teuchos::ArrayRCP<GlobalOrdinal> global_col = Teuchos::ArrayRCP<GlobalOrdinal>(1, Bview->rowMap->getGlobalElement(j));

      C_filled ?
        C.sumIntoGlobalValues(global_row, global_col(), C_ij())
        :
        C.insertGlobalValues(global_row, global_col(), C_ij());

    }
  }

  //delete [] iwork;
  //delete [] bvals;
  //delete [] b_firstcol;

  return(returnValue);
}

//kernel method for computing the local portion of C = A^T*B
template<
  class Scalar, 
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatVec, 
  class SpMatSlv>
int mult_Atrans_B(
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Aview, 
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>&  C)
{
  LocalOrdinal C_firstCol = Bview->colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Bview->colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

  if (Bview->importColMap != Teuchos::null) {
    C_firstCol_import = Bview->importColMap->getMinLocalIndex();
    C_lastCol_import = Bview->importColMap->getMaxLocalIndex();
  }

  size_t C_numCols = C_lastCol - C_firstCol + Teuchos::OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + Teuchos::OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  Teuchos::ArrayRCP<Scalar> C_row_i = Teuchos::ArrayRCP<Scalar>(C_numCols);
  Teuchos::ArrayRCP<GlobalOrdinal> C_colInds = Teuchos::ArrayRCP<GlobalOrdinal>(C_numCols);

  size_t i, j, k;

  for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<C_numCols; ++j) {
    C_row_i[j] = Teuchos::ScalarTraits<Scalar>::zero();
    C_colInds[j] = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();
  }

  //To form C = A^T*B, compute a series of outer-product updates.
  //
  // for (ith column of A^T) { 
  //   C_i = outer product of A^T(:,i) and B(i,:)
  // Where C_i is the ith matrix update,
  //       A^T(:,i) is the ith column of A^T, and
  //       B(i,:) is the ith row of B.
  // }
  //

  //dumpCrsMatrixStruct(Aview->;
  //dumpCrsMatrixStruct(Bview);
  //int localProc = Bview.colMap->Comm().MyPID();
  int localProc = Bview->colMap->getComm()->getRank();

  Teuchos::ArrayView<const GlobalOrdinal> Arows = Aview->rowMap->getNodeElementList();

  bool C_filled = C.isFillComplete();

  //loop over the rows of A (which are the columns of A^T).
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Aview->numRows; ++i) {

    Teuchos::ArrayRCP<const LocalOrdinal> Aindices_i = Aview->indices[i];
    Teuchos::ArrayRCP<const Scalar> Aval_i  = Aview->values[i];

    //we'll need to get the row of B corresponding to Arows[i],
    //where Arows[i] is the GID of A's ith row.
    LocalOrdinal Bi = Bview->rowMap->getLocalElement(Arows[i]);
    TEST_FOR_EXCEPTION(Bi<Teuchos::OrdinalTraits<LocalOrdinal>::zero(), std::runtime_error,
      "mult_Atrans_B ERROR, proc "<<localProc<<" needs row "
     <<Arows[i]<<" of matrix B, but doesn't have it.");

    Teuchos::ArrayRCP<const LocalOrdinal> Bcol_inds = Bview->indices[Bi];
    Teuchos::ArrayRCP<const Scalar> Bvals_i = Bview->values[Bi];

    //for each column-index Aj in the i-th row of A, we'll update
    //global-row GID(Aj) of the result matrix C. In that row of C,
    //we'll update column-indices given by the column-indices in the
    //ith row of B that we're now holding (Bcol_inds).

    //First create a list of GIDs for the column-indices
    //that we'll be updating.

    size_t Blen = Bview->numEntriesPerRow[Bi];
    if (Bview->remote[Bi]) {
      for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Blen; ++j) {
        C_colInds[j] = Bview->importColMap->getGlobalElement(Bcol_inds[j]);
      }
    }
    else {
      for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Blen; ++j) {
        C_colInds[j] = Bview->colMap->getGlobalElement(Bcol_inds[j]);
      }
    }

    //loop across the i-th row of A (column of A^T)
    for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Aview->numEntriesPerRow[i]; ++j) {

      LocalOrdinal Aj = Aindices_i[j];
      Scalar Aval = Aval_i[j];

      GlobalOrdinal global_row;
      if (Aview->remote[i]) {
        global_row = Aview->importColMap->getGlobalElement(Aj);
      }
      else {
        global_row = Aview->colMap->getGlobalElement(Aj);
      }

      if (!C.getRowMap()->isNodeGlobalElement(global_row)) {
        continue;
      }

      for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<Blen; ++k) {
        C_row_i[k] = Aval*Bvals_i[k];
      }

      //
      //Now add this row-update to C.
      //

      C_filled ?
        C.sumIntoGlobalValues(global_row, C_colInds(), C_row_i() )
        :
        C.insertGlobalValues(global_row, C_colInds(), C_row_i());

    }
  }

  //delete [] C_row_i;
  //delete [] C_colInds;

  return(0);
}

//kernel method for computing the local portion of C = A^T*B^T
template<
  class Scalar,
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatVec, 
  class SpMatSlv>
int mult_Atrans_Btrans(
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Aview, 
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C)
{
  LocalOrdinal C_firstCol = Aview->rowMap->getMinLocalIndex();
  LocalOrdinal C_lastCol = Aview->rowMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = Teuchos::OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

  if (Aview->importColMap != Teuchos::null) {
    C_firstCol_import = Aview->importColMap->getMinLocalIndex();
    C_lastCol_import = Aview->importColMap->getMaxLocalIndex();
  }

  size_t C_numCols = C_lastCol - C_firstCol + Teuchos::OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import + Teuchos::OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  Teuchos::ArrayRCP<Scalar> dwork = Teuchos::ArrayRCP<Scalar>(C_numCols);
  Teuchos::ArrayRCP<GlobalOrdinal> iwork = Teuchos::ArrayRCP<GlobalOrdinal>(C_numCols);

  Teuchos::ArrayRCP<Scalar> C_col_j = dwork;
  Teuchos::ArrayRCP<GlobalOrdinal> C_inds = iwork;

  //cout << "Aview: " << endl;
  //dumpCrsMatrixStruct(Aview);

  //cout << "Bview: " << endl;
  //dumpCrsMatrixStruct(Bview);


  size_t i, j, k;

  for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<C_numCols; ++j) {
    C_col_j[j] = Teuchos::ScalarTraits<Scalar>::zero();
    C_inds[j] = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
  }

  Teuchos::ArrayView<const GlobalOrdinal> A_col_inds = Aview->colMap->getNodeElementList();
  Teuchos::ArrayView<const GlobalOrdinal> A_col_inds_import = Aview->importColMap == Teuchos::null ?
    Aview->importColMap->getNodeElementList() 
	:
	Teuchos::null;

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > Crowmap = C.getRowMap();

  //To form C = A^T*B^T, we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(k,i)*B(j,k) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc. In other words, we avoid
  //column-wise operations like the plague...

  Teuchos::ArrayView<const GlobalOrdinal> Brows = Bview->rowMap->getNodeElementList();

  //loop over the rows of B
  for(j=Teuchos::OrdinalTraits<size_t>::zero(); j<Bview->numRows; ++j) {
    Teuchos::ArrayRCP<const LocalOrdinal> Bindices_j = Bview->indices[j];
    Teuchos::ArrayRCP<const Scalar> Bvals_j = Bview->values[j];

    //GlobalOrdinal global_col = Brows[j];
    Teuchos::ArrayView<const GlobalOrdinal> global_col = Brows.view(j,1);

    //loop across columns in the j-th row of B and for each corresponding
    //row in A, loop across columns and accumulate product
    //A(k,i)*B(j,k) into our partial sum quantities in C_col_j. In other
    //words, as we stride across B(j,:), we use selected rows in A to
    //calculate updates for column j of the result matrix C.

    for(k=Teuchos::OrdinalTraits<size_t>::zero(); k<Bview->numEntriesPerRow[j]; ++k) {

      LocalOrdinal bk = Bindices_j[k];
      Scalar Bval = Bvals_j[k];

      GlobalOrdinal global_k;
      if (Bview->remote[j]) {
        global_k = Bview->importColMap->getGlobalElement(bk);
      }
      else {
        global_k = Bview->colMap->getGlobalElement(bk);
      }

      //get the corresponding row in A
      LocalOrdinal ak = Aview->rowMap->getLocalElement(global_k);
      if (ak<Teuchos::OrdinalTraits<LocalOrdinal>::zero()) {
        continue;
      }

      Teuchos::ArrayRCP<const LocalOrdinal> Aindices_k = Aview->indices[ak];
      Teuchos::ArrayRCP<const Scalar> Avals_k = Aview->values[ak];

      size_t C_len = Teuchos::OrdinalTraits<size_t>::zero();

      if (Aview->remote[ak]) {
        for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Aview->numEntriesPerRow[ak]; ++i) {
          C_col_j[C_len] = Avals_k[i]*Bval;
          C_inds[C_len++] = A_col_inds_import[Aindices_k[i]];
        }
      }
      else {
        for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Aview->numEntriesPerRow[ak]; ++i) {
          C_col_j[C_len] = Avals_k[i]*Bval;
          C_inds[C_len++] = A_col_inds[Aindices_k[i]];
        }
      }

      //Now loop across the C_col_j values and put non-zeros into C.

      for(i=Teuchos::OrdinalTraits<size_t>::zero(); i < C_len ; ++i) {
        if (C_col_j[i] == Teuchos::ScalarTraits<Scalar>::zero()) continue;

        GlobalOrdinal global_row = C_inds[i];
        if (!Crowmap->isNodeGlobalElement(global_row)) {
          continue;
        }

/*  int err = C.SumIntoGlobalValues(global_row, 1, &(C_col_j[i]), &global_col);

  if (err < 0) {
    return(err);
  }
  else {
          if (err > 0) {
      err = C.InsertGlobalValues(global_row, 1, &(C_col_j[i]), &global_col);
      if (err < 0) {
              return(err);
            }
    }
  }*/

        try{
          //C.sumIntoGlobalValues(global_row, global_col, C_col_j[i]);
          C.sumIntoGlobalValues(global_row, global_col, C_col_j.view(i,1));
        }
        catch(std::runtime_error){
          //C.insertGlobalValues(global_row, global_col, C_col_j[i]);
          C.insertGlobalValues(global_row, global_col, C_col_j.view(i,1));
        }
      }
    }
  }

  //delete [] dwork;
  //delete [] iwork;

  return(0);
}

template<
  class Scalar, 
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatVec, 
  class SpMatSlv>
int import_and_extract_views(
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& targetMap,
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> >& Mview)
{
  //The goal of this method is to populate the 'Mview' struct with views of the
  //rows of M, including all rows that correspond to elements in 'targetMap'.
  //
  //If targetMap includes local elements that correspond to remotely-owned rows
  //of M, then those remotely-owned rows will be imported into
  //'Mview.importMatrix', and views of them will be included in 'Mview'.

  Mview->deleteContents();

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > Mrowmap = M->getRowMap();

  int numProcs = Mrowmap->getComm()->getSize();

  Mview->numRows = targetMap->getNodeNumElements();

  Teuchos::ArrayView<const GlobalOrdinal> Mrows = targetMap->getNodeElementList();

  if (Mview->numRows > Teuchos::OrdinalTraits<size_t>::zero()) {
    Mview->numEntriesPerRow = Teuchos::ArrayRCP<size_t>(Mview->numRows);
    Mview->indices = Teuchos::ArrayRCP<Teuchos::ArrayRCP<const LocalOrdinal> >(Mview->numRows, Teuchos::null);
    Mview->values = Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >(Mview->numRows, Teuchos::null);
    Mview->remote = Teuchos::ArrayRCP<bool>(Mview->numRows);
  }

  Mview->numRemote = Teuchos::OrdinalTraits<global_size_t>::zero();

  size_t i;
  for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Mview->numRows; ++i) {
    LocalOrdinal mlid = Mrowmap->getLocalElement(Mrows[i]);
    if (mlid < Teuchos::OrdinalTraits<LocalOrdinal>::zero()) {
      Mview->remote[i] = true;
      ++Mview->numRemote;
    }
    else {
  //    EPETRA_CHK_ERR( M->ExtractMyRowView(mlid, Mview.numEntriesPerRow[i],
   //        Mview.values[i], Mview.indices[i]) );
      M->getLocalRowView(mlid, Mview->indices[i], Mview->values[i]);
	  Mview->numEntriesPerRow[i] = Mview->indices[i].size();
      Mview->remote[i] = false;
    }
  }


  Mview->origRowMap = M->getRowMap();
  Mview->rowMap = targetMap;
  Mview->colMap = M->getColMap();
  Mview->domainMap = M->getDomainMap();
  Mview->importColMap = Teuchos::null;

  if (numProcs < 2) {

    TEST_FOR_EXCEPTION(Mview->numRemote > Teuchos::OrdinalTraits<global_size_t>::zero(), std::runtime_error,
      "MatrixMatrix::Multiply ERROR, numProcs < 2 but attempting to import remote matrix rows." <<std::endl);
      return(-1);
    

    //If only one processor we don't need to import any remote rows, so return.
    return(0);
  }

  //
  //Now we will import the needed remote rows of M, if the global maximum
  //value of numRemote is greater than 0.
  //

  //int globalMaxNumRemote = 0;
  global_size_t globalMaxNumRemote = Teuchos::OrdinalTraits<global_size_t>::zero();
  //Mrowmap.Comm().MaxAll(&Mview.numRemote, &globalMaxNumRemote, 1);
  Teuchos::reduceAll(*(Mrowmap->getComm()) , Teuchos::REDUCE_MAX, 1, &Mview->numRemote, &globalMaxNumRemote);

  if(globalMaxNumRemote > Teuchos::OrdinalTraits<global_size_t>::zero()){
    //Create a map that describes the remote rows of M that we need.

    //int* MremoteRows = Mview.numRemote>0 ? new int[Mview.numRemote] : NULL;
    Teuchos::ArrayRCP<GlobalOrdinal> MremoteRows = Mview->numRemote>Teuchos::OrdinalTraits<global_size_t>::zero() ? Teuchos::ArrayRCP<GlobalOrdinal>(Mview->numRemote) : Teuchos::null;

    //int offset = 0;
    global_size_t offset = Teuchos::OrdinalTraits<global_size_t>::zero();
    for(i=0; i<Mview->numRows; ++i) {
      if (Mview->remote[i]) {
        MremoteRows[offset++] = Mrows[i];
      }
    }

    //Epetra_Map MremoteRowMap(-1, Mview->numRemote, MremoteRows,
     //      Mrowmap.IndexBase(), Mrowmap.Comm());
    Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > MremoteRowMap = Teuchos::rcp(new Map<LocalOrdinal, GlobalOrdinal, Node>(Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), MremoteRows(), Mrowmap->getIndexBase(), Mrowmap->getComm(), Mrowmap->getNode()));

    //Create an importer with target-map MremoteRowMap and 
    //source-map Mrowmap.
    //Epetra_Import importer(MremoteRowMap, Mrowmap);
    Import<LocalOrdinal, GlobalOrdinal, Node> importer(Mrowmap, MremoteRowMap);

    //Now create a new matrix into which we can import the remote rows of M
    //that we need.
    Mview->importMatrix = Teuchos::rcp(new CrsMatrix<Scalar,LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>( MremoteRowMap, 1));

    //EPETRA_CHK_ERR( Mview.importMatrix->Import(M, importer, Insert) );
    Mview->importMatrix->doImport(*M, importer, INSERT);

    //EPETRA_CHK_ERR( Mview.importMatrix->FillComplete(M.DomainMap(), M.RangeMap()) );
    Mview->importMatrix->fillComplete(M->getDomainMap(), M->getRangeMap());

    //Finally, use the freshly imported data to fill in the gaps in our views
    //of rows of M.
    for(i=Teuchos::OrdinalTraits<size_t>::zero(); i<Mview->numRows; ++i) {
      if (Mview->remote[i]) {
        LocalOrdinal importLID = MremoteRowMap->getLocalElement(Mrows[i]);
        Mview->importMatrix->getLocalRowView(importLID,
          Mview->indices[i],
          Mview->values[i]);
		Mview->numEntriesPerRow[i] = Mview->indices[i].size();
      }
    }

    Mview->importColMap = Mview->importMatrix->getColMap();

    //delete [] MremoteRows;
  }

  return(0);
}

template<class Ordinal, class GlobalOrdinal>
int distribute_list(const Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm,
  size_t lenSendList,
  const Teuchos::ArrayRCP<GlobalOrdinal>& sendList,
  size_t& maxSendLen,
  Teuchos::ArrayRCP<GlobalOrdinal>& recvList)
{
  maxSendLen = Teuchos::OrdinalTraits<size_t>::zero() ; 
  //Comm.MaxAll(&lenSendList, &maxSendLen, 1);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::OrdinalTraits<Ordinal>::one(), &lenSendList, &maxSendLen);
  int numProcs = comm->getSize();
  //recvList = new int[numProcs*maxSendLen];
  recvList = Teuchos::ArrayRCP<GlobalOrdinal>(numProcs*maxSendLen);
  Teuchos::ArrayRCP<GlobalOrdinal> send = Teuchos::ArrayRCP<GlobalOrdinal>(maxSendLen);
  for(size_t i=Teuchos::OrdinalTraits<size_t>::zero(); i<lenSendList; ++i) {
    send[i] = sendList[i];
  }

  Teuchos::gatherAll(*comm, (Ordinal)maxSendLen, send.getRawPtr(), (Ordinal)maxSendLen, recvList.getRawPtr());
  //delete [] send;

  return(0);
}

template<
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > 
create_map_from_imported_rows(
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > map,
  size_t totalNumSend,
  Teuchos::ArrayRCP<GlobalOrdinal> sendRows,
  int numProcs,
  Teuchos::ArrayRCP<size_t> numSendPerProc)
{
  //Perform sparse all-to-all communication to send the row-GIDs
  //in sendRows to appropriate processors according to offset
  //information in numSendPerProc.
  //Then create and return a map containing the rows that we
  //received on the local processor.

  Teuchos::RCP<Distributor> distributor = rcp(new Distributor(map->getComm()));

  Teuchos::ArrayRCP<int> sendPIDs = totalNumSend>0 ? Teuchos::ArrayRCP<int>(totalNumSend) : Teuchos::null;
  int offset = 0;
  for(int i=0; i<numProcs; ++i) {
    for(size_t j=0; j<numSendPerProc[i]; ++j) {
      sendPIDs[offset++] = i;
    }
  }

  //int numRecv = 0;
  //int err = distributor->CreateFromSends(totalNumSend, sendPIDs,
      //     true, numRecv);
  size_t numRecv = distributor->createFromSends(sendPIDs());
  //assert( err == 0 );

  //char* c_recv_objs = numRecv>0 ? new char[numRecv*sizeof(int)] : NULL;
  Teuchos::ArrayRCP<GlobalOrdinal> recv_rows = 
    numRecv>Teuchos::OrdinalTraits<size_t>::zero() ? Teuchos::ArrayRCP<GlobalOrdinal>(numRecv) : Teuchos::null;
  //int num_c_recv = numRecv*(int)sizeof(int);

  //err = distributor->Do(reinterpret_cast<char*>(sendRows),
      //(int)sizeof(int), num_c_recv, c_recv_objs);
  //ArrayView<GlobalOrdinal> sendRowsView = sendRows();

  distributor->doPostsAndWaits(sendRows.getConst()(), numRecv, recv_rows());
  //assert( err == 0 );

  //int* recvRows = reinterpret_cast<int*>(c_recv_objs);

  //Now create a map with the rows we've received from other processors.
  //Epetra_Map* import_rows = new Epetra_Map(-1, numRecv, recvRows,
   //          map->IndexBase(), map->Comm());
  Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > import_rows = 
    Teuchos::rcp(new Map<LocalOrdinal, GlobalOrdinal, Node>(
    Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), recv_rows(), map->getIndexBase(), map->getComm()));

  //delete [] c_recv_objs;
  //delete [] sendPIDs;

  //delete distributor;

  return( import_rows );
}

template<
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node>
int form_map_union(
  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > map1,
  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > map2,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& mapunion)
{
  //form the union of two maps

  if (map1 == Teuchos::null) {
    //mapunion = new Epetra_Map(*map2);
    mapunion = map2;
    return(0);
  }

  if (map2 == Teuchos::null) {
    //mapunion = new Epetra_Map(*map1);
    mapunion = map1;
    return(0);
  }

  size_t map1_len = map1->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> map1_elements = map1->getNodeElementList();
  size_t map2_len = map2->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> map2_elements = map2->getNodeElementList();

  Teuchos::ArrayRCP<GlobalOrdinal> union_elements = Teuchos::ArrayRCP<GlobalOrdinal>(map1_len+map2_len);

  //int map1_offset = 0, map2_offset = 0, union_offset = 0;
  global_size_t map1_offset = Teuchos::OrdinalTraits<global_size_t>::zero(), map2_offset = Teuchos::OrdinalTraits<global_size_t>::zero(), union_offset = Teuchos::OrdinalTraits<global_size_t>::zero();

  while(map1_offset < map1_len && map2_offset < map2_len) {
    GlobalOrdinal map1_elem = map1_elements[map1_offset];
    GlobalOrdinal map2_elem = map2_elements[map2_offset];

    if (map1_elem < map2_elem) {
      union_elements[union_offset++] = map1_elem;
      ++map1_offset;
    }
    else if (map1_elem > map2_elem) {
      union_elements[union_offset++] = map2_elem;
      ++map2_offset;
    }
    else {
      union_elements[union_offset++] = map1_elem;
      ++map1_offset;
      ++map2_offset;
    }
  }

  global_size_t i;
  for(i=map1_offset; i<map1_len; ++i) {
    union_elements[union_offset++] = map1_elements[i];
  }

  for(i=map2_offset; i<map2_len; ++i) {
    union_elements[union_offset++] = map2_elements[i];
  }

  //mapunion = new Epetra_Map(-1, union_offset, union_elements,
   //       map1->IndexBase(), map1->Comm());
  mapunion = Teuchos::rcp(new Map<LocalOrdinal, GlobalOrdinal, Node>(Teuchos::OrdinalTraits<global_size_t>::invalid(), union_elements(), map1->getIndexBase(), map1->getComm(), map1->getNode()));

  //delete [] union_elements;

  return(0);
}

template<
  class Scalar,
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatVec, 
  class SpMatSlv >
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
find_rows_containing_cols(
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > colmap)
{
  //The goal of this function is to find all rows in the matrix M that contain
  //column-indices which are in 'colmap'. A map containing those rows is
  //returned.

  int numProcs = colmap->getComm()->getSize();
  int localProc = colmap->getComm()->getRank();

  if (numProcs < 2) {
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > result_map = Teuchos::null;

    int err = form_map_union(
	  M->getRowMap(), 
	  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::null), 
	  result_map);
    if (err != 0) {
      return(Teuchos::null);
    }
    return(result_map);
  }

  size_t MnumRows = M->getNodeNumRows();
  size_t numCols = colmap->getNodeNumElements();

  /*int* iwork = new int[numCols+2*numProcs+numProcs*MnumRows];
  int iworkOffset = 0;*/
  //ArrayRCP<GlobalOrdinal> iwork = Teuchos::ArrayRCP<GlobalOrdinal>(numCols+2*numProcs + numProcs*MnumRows);
  //global_size_t iworkOffset = Teuchos::OrdinalTraits<global_size_t>::zero();

 // int* cols = &(iwork[iworkOffset]); iworkOffset += numCols;
  Teuchos::ArrayRCP<GlobalOrdinal> cols = Teuchos::ArrayRCP<GlobalOrdinal>(numCols + Teuchos::OrdinalTraits<size_t>::one());
  cols[0] = numCols;
  ///cols[1] = colmap->getNodeElementList();
  cols.view(1,numCols).assign(colmap->getNodeElementList());

  //cols are not necessarily sorted at this point, so we'll make sure
  //they are sorted.
  //Epetra_Util util;
  //util.Sort(true, numCols, &(cols[1]), 0, NULL, 0, NULL);
  //sort2(cols.begin(), cols.end(), Teuchos::null);
  std::sort(cols.begin(), cols.end());

// int* all_proc_cols = NULL;
  Teuchos::ArrayRCP<GlobalOrdinal> all_proc_cols = Teuchos::null;
  
  size_t max_num_cols;
  distribute_list(colmap->getComm(), numCols+Teuchos::OrdinalTraits<size_t>::one(), cols, max_num_cols, all_proc_cols);

  Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > Mgraph = M->getCrsGraph();
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > Mrowmap = M->getRowMap();
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > Mcolmap = M->getColMap();
  LocalOrdinal MminMyLID = Mrowmap->getMinLocalIndex();

  /*int* procNumCols = &(iwork[iworkOffset]); iworkOffset += numProcs;
  int* procNumRows = &(iwork[iworkOffset]); iworkOffset += numProcs;
  int* procRows_1D = &(iwork[iworkOffset]);
  int** procCols = new int*[numProcs];
  int** procRows = new int*[numProcs];*/
  Teuchos::ArrayRCP<size_t> procNumCols(numProcs);
  Teuchos::ArrayRCP<size_t> procNumRows(numProcs);
  Teuchos::ArrayRCP<GlobalOrdinal> procRows_1D(numProcs*MnumRows);
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > procCols (numProcs);
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> > procRows = Teuchos::ArrayRCP<Teuchos::ArrayRCP<GlobalOrdinal> >(numProcs);
  int i, err;
  //int offset = 0;
  size_t offset = Teuchos::OrdinalTraits<size_t>::zero();
  for(i=0; i<numProcs; ++i) {
    procNumCols[i] = all_proc_cols[offset];
    procCols[i] = all_proc_cols.persistingView(offset+1, max_num_cols);
    offset += max_num_cols;

    procNumRows[i] = Teuchos::OrdinalTraits<size_t>::zero();
    procRows[i] = procRows_1D.persistingView(i*MnumRows, MnumRows);
  }

  Teuchos::ArrayRCP<const LocalOrdinal> Mindices;

  for(size_t row=Teuchos::OrdinalTraits<size_t>::zero(); row<MnumRows; ++row) {
    LocalOrdinal localRow = MminMyLID+row;
    GlobalOrdinal globalRow = Mrowmap->getGlobalElement(localRow);
    //int MnumCols;
    //err = Mgraph.ExtractMyRowView(localRow, MnumCols, Mindices);
    Mindices = Mgraph->getLocalRowView(localRow);
    /*if (err != 0) {
      cerr << "proc "<<localProc<<", error in Mgraph.ExtractMyRowView, row "
           <<localRow<<endl;
      return(NULL);
    }*/

    for(LocalOrdinal j=Teuchos::OrdinalTraits<size_t>::zero(); j<Mindices.size(); ++j) {
      GlobalOrdinal colGID = Mcolmap->getGlobalElement(Mindices[j]);

      for(int p=0; p<numProcs; ++p) {
        if (p==localProc) continue;

        /*int insertPoint;
        int foundOffset = Epetra_Util_binary_search(colGID, procCols[p],
                                                    procNumCols[p], insertPoint);*/
        typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator result = binary_serach(procCols[p].begin(), procCols[p].end(), colGID);      
        if (result != procCols[p].end()) {
          size_t numRowsP = procNumRows[p];
          Teuchos::ArrayRCP<GlobalOrdinal> prows = procRows[p];
          if (numRowsP < Teuchos::OrdinalTraits<size_t>::one() || prows[numRowsP-Teuchos::OrdinalTraits<size_t>::one()] < globalRow) {
            prows[numRowsP] = globalRow;
            procNumRows[p]++;
          }
        }
      }
    }
  }

  //Now make the contents of procRows occupy a contiguous section
  //of procRows_1D.
  offset = procNumRows[0];
  for(i=Teuchos::OrdinalTraits<global_size_t>::one(); i<numProcs; ++i) {
    for(size_t j=Teuchos::OrdinalTraits<size_t>::one(); j<procNumRows[i]; ++j) {
      procRows_1D[offset++] = procRows[i][j];
    }
  }

  size_t totalNumSend = offset;
  //Next we will do a sparse all-to-all communication to send the lists of rows
  //to the appropriate processors, and create a map with the rows we've received
  //from other processors.
  /*Epetra_Map* recvd_rows =
    create_map_from_imported_rows(&Mrowmap, totalNumSend,
                                  procRows_1D, numProcs, procNumRows);*/
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > recvd_rows = create_map_from_imported_rows(Mrowmap, totalNumSend, procRows_1D, numProcs, procNumRows);
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > result_map = Teuchos::null;

  err = form_map_union(M->getRowMap(), recvd_rows, result_map);
  if (err != 0) {
    return(Teuchos::null);
  }

  /*delete [] iwork;
  delete [] procCols;
  delete [] procRows;
  delete [] all_proc_cols;
  delete recvd_rows;*/

  return(result_map);
}


template <
  class Scalar, 
  class LocalOrdinal,
  class GlobalOrdinal,
  class Node,
  class SpMatVec,
  class SpMatSlv >
int MatrixMatrix::Multiply(
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > A,
  bool transposeA,
  Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > B,
  bool transposeB,
  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > C,
  bool call_FillComplete_on_result)
{
  //
  //This method forms the matrix-matrix product C = op(A) * op(B), where
  //op(A) == A   if transposeA is false,
  //op(A) == A^T if transposeA is true,
  //and similarly for op(B).
  //

  //A and B should already be Filled.
  //(Should we go ahead and call FillComplete() on them if necessary?
  // or error out? For now, we choose to error out.)
  TEST_FOR_EXCEPTION(!A->isFillComplete(), std::runtime_error,
    "Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
    "a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
    "\n\n MatrixMatrix::Multiply(): Matrix A is not fill complete.");
  TEST_FOR_EXCEPTION(!B->isFillComplete(), std::runtime_error,
    "Uh oh. Looks like there's a bit of a problem here. No worries though. We'll help you figure it out. You're "
    "a fantastic programer and this just a minor bump in the road! Maybe the information below can help you out a bit."
    "\n\n MatrixMatrix::Multiply(): Matrix B is not fill complete.");

  //We're going to refer to the different combinations of op(A) and op(B)
  //as scenario 1 through 4.

  int scenario = 1;//A*B
  if (transposeB && !transposeA) scenario = 2;//A*B^T
  if (transposeA && !transposeB) scenario = 3;//A^T*B
  if (transposeA && transposeB)  scenario = 4;//A^T*B^T

  //now check size compatibility
  global_size_t Aouter = transposeA ? A->getGlobalNumCols() : A->getGlobalNumRows();
  global_size_t Bouter = transposeB ? B->getGlobalNumRows() : B->getGlobalNumCols();
  global_size_t Ainner = transposeA ? A->getGlobalNumRows() : A->getGlobalNumCols();
  global_size_t Binner = transposeB ? B->getGlobalNumCols() : B->getGlobalNumRows();
  TEST_FOR_EXCEPTION(!A->isFillComplete(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    <<Aouter<<"x"<<Ainner << ", op(B) is "<<Binner<<"x"<<Bouter<<std::endl);

  //The result matrix C must at least have a row-map that reflects the
  //correct row-size. Don't check the number of columns because rectangular
  //matrices which were constructed with only one map can still end up
  //having the correct capacity and dimensions when filled.
  TEST_FOR_EXCEPTION(Aouter > C->getGlobalNumRows(), std::runtime_error,
    "MatrixMatrix::Multiply: ERROR, dimensions of result C must "
    "match dimensions of op(A) * op(B). C has "<<C->getGlobalNumRows()
     << " rows, should have at least "<<Aouter << std::endl);

  //It doesn't matter whether C is already Filled or not. If it is already
  //Filled, it must have space allocated for the positions that will be
  //referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  //we'll error out later when trying to store result values.

  //We're going to need to import remotely-owned sections of A and/or B
  //if more than 1 processor is performing this run, depending on the scenario.
  int numProcs = A->getComm()->getSize();

  //If we are to use the transpose of A and/or B, we'll need to be able to 
  //access, on the local processor, all rows that contain column-indices in
  //the domain-map.
//  const Epetra_Map* domainMap_A = &(A.DomainMap());
//  const Epetra_Map* domainMap_B = &(B.DomainMap());

  //const Epetra_Map* rowmap_A = &(A.RowMap());
  //const Epetra_Map* rowmap_B = &(B.RowMap());

  //Declare some 'work-space' maps which may be created depending on
  //the scenario, and which will be deleted before exiting this function.
  //const Epetra_Map* workmap1 = NULL;
  //const Epetra_Map* workmap2 = NULL;
  //const Epetra_Map* mapunion1 = NULL;

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > workmap1 = Teuchos::null;
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > workmap2 = Teuchos::null;
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > mapunion1 = Teuchos::null;

  //Declare a couple of structs that will be used to hold views of the data
  //of A and B, to be used for fast access during the matrix-multiplication.
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > Aview = Teuchos::rcp(new CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>);
  Teuchos::RCP<CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv> > Bview = Teuchos::rcp(new CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatVec, SpMatSlv>);

  //const Epetra_Map* targetMap_A = rowmap_A;
  //const Epetra_Map* targetMap_B = rowmap_B;

  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap_A = A->getRowMap();
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap_B = B->getRowMap();

  if (numProcs > 1) {
    //If op(A) = A^T, find all rows of A that contain column-indices in the
    //local portion of the domain-map. (We'll import any remote rows
    //that fit this criteria onto the local processor.)
    if (transposeA) {
      workmap1 = find_rows_containing_cols(A, A->getDomainMap());
      targetMap_A = workmap1;
    }
  }

  //Now import any needed remote rows and populate the Aview struct.
  //EPETRA_CHK_ERR( import_and_extract_views(A, *targetMap_A, Aview) );
  import_and_extract_views(A, targetMap_A, Aview);

    

  //We will also need local access to all rows of B that correspond to the
  //column-map of op(A).
  if (numProcs > 1) {
    //const Epetra_Map* colmap_op_A = NULL;
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > colmap_op_A = Teuchos::null;
    if (transposeA) {
      colmap_op_A = targetMap_A;
    }
    else {
      colmap_op_A = A->getColMap(); 
    }

    targetMap_B = colmap_op_A;

    //If op(B) = B^T, find all rows of B that contain column-indices in the
    //local-portion of the domain-map, or in the column-map of op(A).
    //We'll import any remote rows that fit this criteria onto the
    //local processor.
    if (transposeB) {
      form_map_union(colmap_op_A, B->getDomainMap(), mapunion1);
      workmap2 = find_rows_containing_cols(B, mapunion1);
      targetMap_B = workmap2;
    }
  }

  //Now import any needed remote rows and populate the Bview struct.
  import_and_extract_views(B, targetMap_B, Bview);


  //If the result matrix C is not already FillComplete'd, we will do a
  //preprocessing step to create the nonzero structure, then call FillComplete,
  if (!C->isFillComplete()) {
    CrsWrapper_GraphBuilder<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsgraphbuilder(C->getRowMap());

    //pass the graph-builder object to the multiplication kernel to fill in all
    //the nonzero positions that will be used in the result matrix.
    switch(scenario) {
    case 1:    mult_A_B(Aview, Bview, crsgraphbuilder);
      break;
    case 2:    mult_A_Btrans(Aview, Bview, crsgraphbuilder);
      break;
    case 3:    mult_Atrans_B(Aview, Bview, crsgraphbuilder);
      break;
    case 4:    mult_Atrans_Btrans(Aview, Bview, crsgraphbuilder);
      break;
    }

    //now insert all of the nonzero positions into the result matrix.
    insert_matrix_locations(crsgraphbuilder, C);

  /*  if (call_FillComplete_on_result) {
      Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domainmap = transposeB ? B->getRangeMap() : B->getDomainMap();

      Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > rangemap = transposeA ? A->getDomainMap() : A->getRangeMap();

      C->fillComplete(domainmap, rangemap, DoNotOptimizeStorage);
      call_FillComplete_on_result = false;
    }*/
  }

  //Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
  //C->describe(*out, Teuchos::VERB_EXTREME);

  //Pre-zero the result matrix:
  //C->setAllToScalar(Teuchos::ScalarTraits<Scalar>::zero());

  //Now call the appropriate method to perform the actual multiplication.

  CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

  switch(scenario) {
  case 1:    mult_A_B(Aview, Bview, crsmat);
    break;
  case 2:    mult_A_Btrans(Aview, Bview, crsmat);
    break;
  case 3:    mult_Atrans_B(Aview, Bview, crsmat);
    break;
  case 4:    mult_Atrans_Btrans(Aview, Bview, crsmat);
    break;
  }

  if (call_FillComplete_on_result) {
    //We'll call FillComplete on the C matrix before we exit, and give
    //it a domain-map and a range-map.
    //The domain-map will be the domain-map of B, unless
    //op(B)==transpose(B), in which case the range-map of B will be used.
    //The range-map will be the range-map of A, unless
    //op(A)==transpose(A), in which case the domain-map of A will be used.
    if (!C->isFillComplete()) {
      Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > domainmap = transposeB ? B->getRangeMap() : B->getDomainMap();

      Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > rangemap = transposeA ? A->getDomainMap() : A->getRangeMap();
      //C->fillComplete(transposeB ? B->getRangeMap() : B->getDomainMap(), transposeA ? A->getDomainMap() : B->getRangeMap());
      C->fillComplete(domainmap, rangemap);
    }
  }

  //Finally, delete the objects that were potentially created
  //during the course of importing remote sections of A and B.

  //delete mapunion1; mapunion1 = NULL;
  //delete workmap1; workmap1 = NULL;
  //delete workmap2; workmap2 = NULL;

  return(0);
}
/*
int MatrixMatrix::Add(Teuchos::RCP<const CrsMatrixType> A,
                      bool transposeA,
                      double scalarA,
                      Teuchos::RCP<CrsMatrixType> B,
                      double scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  //A should already be Filled. It doesn't matter whether B is
  //already Filled, but if it is, then its graph must already contain
  //all nonzero locations that will be referenced in forming the
  //sum.

  TEST_FOR_EXCEPTION(!A->isFillComplete(), std::runtime_error,
    "MatrixMatrix::Add ERROR, input matrix A.isFillComplete() is false, it is required to be true. (Result matrix B is not required to be isFillComplete()).");

  //explicit tranpose A formed as necessary
  Epetra_CrsMatrix * Aprime = 0;
  Teuchos::RCP<CrsMatrixType> Aprime = Teuchos::null;
  EpetraExt::RowMatrix_Transpose * Atrans = 0;
  if( transposeA )
  {
    Atrans = new EpetraExt::RowMatrix_Transpose();
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  int MaxNumEntries = EPETRA_MAX( A.MaxNumEntries(), B.MaxNumEntries() );
  int A_NumEntries, B_NumEntries;
  int * A_Indices = new int[MaxNumEntries];
  double * A_Values = new double[MaxNumEntries];
  int* B_Indices;
  double* B_Values;

  int NumMyRows = B.NumMyRows();
  int Row, err;
  int ierr = 0;

  if( scalarA )
  {
    //Loop over B's rows and sum into
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = B.GRID(i);
      EPETRA_CHK_ERR( Aprime->ExtractGlobalRowCopy( Row, MaxNumEntries, A_NumEntries, A_Values, A_Indices ) );

      if (scalarB != 1.0) {
        if (!B.Filled()) {
          EPETRA_CHK_ERR( B.ExtractGlobalRowView( Row, B_NumEntries,
                                                  B_Values, B_Indices));
        }
        else {
          EPETRA_CHK_ERR( B.ExtractMyRowView( i, B_NumEntries,
                                              B_Values, B_Indices));
        }

        for(int jj=0; jj<B_NumEntries; ++jj) {
          B_Values[jj] = scalarB*B_Values[jj];
        }
      }

      if( scalarA != 1.0 ) {
        for( int j = 0; j < A_NumEntries; ++j ) A_Values[j] *= scalarA;
      }

      if( B.Filled() ) {//Sum In Values
        err = B.SumIntoGlobalValues( Row, A_NumEntries, A_Values, A_Indices );
        assert( err >= 0 );
        if (err < 0) ierr = err;
      }
      else {
        err = B.InsertGlobalValues( Row, A_NumEntries, A_Values, A_Indices );
        assert( err == 0 || err == 1 || err == 3 );
        if (err < 0) ierr = err;
      }
    }
  }
  else {
    EPETRA_CHK_ERR( B.Scale(scalarB) );
  }

  delete [] A_Indices;
  delete [] A_Values;

  if( Atrans ) delete Atrans;

  return(ierr);
}

int MatrixMatrix::Add(const Epetra_CrsMatrix& A,
                      bool transposeA,
                      double scalarA,
                      const Epetra_CrsMatrix & B,
                      bool transposeB,
                      double scalarB,
                      Epetra_CrsMatrix * & C)
{
  //
  //This method forms the matrix-matrix sum C = scalarA * op(A) + scalarB * op(B), where

  //A and B should already be Filled. C should be an empty pointer.

  if (!A.Filled() || !B.Filled() ) {
     std::cerr << "EpetraExt::MatrixMatrix::Add ERROR, input matrix A.Filled() or B.Filled() is false,"
               << "they are required to be true. (Result matrix C should be an empty pointer)" << std::endl;
     EPETRA_CHK_ERR(-1);
  }

  Epetra_CrsMatrix * Aprime = 0, * Bprime=0;
  EpetraExt::RowMatrix_Transpose * Atrans = 0,* Btrans = 0;

  //explicit tranpose A formed as necessary
  if( transposeA ) {
     Atrans = new EpetraExt::RowMatrix_Transpose();
     Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
     Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  //explicit tranpose B formed as necessary
  if( transposeB ) {
     Btrans = new EpetraExt::RowMatrix_Transpose();
     Bprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Btrans)(const_cast<Epetra_CrsMatrix&>(B)))));
  }
  else
     Bprime = const_cast<Epetra_CrsMatrix*>(&B);

  // allocate or zero the new matrix
  if(C!=0)
     C->PutScalar(0.0);
  else
     C = new Epetra_CrsMatrix(Copy,Aprime->RowMap(),0);

  // build arrays  for easy resuse
  int ierr = 0;
  Epetra_CrsMatrix * Mat[] = { Aprime,Bprime};
  double scalar[] = { scalarA, scalarB};

  // do a loop over each matrix to add: A reordering might be more efficient
  for(int k=0;k<2;k++) {
     int MaxNumEntries = Mat[k]->MaxNumEntries();
     int NumEntries;
     int * Indices = new int[MaxNumEntries];
     double * Values = new double[MaxNumEntries];
   
     int NumMyRows = Mat[k]->NumMyRows();
     int Row, err;
     int ierr = 0;
   
     //Loop over rows and sum into C
     for( int i = 0; i < NumMyRows; ++i ) {
        Row = Mat[k]->GRID(i);
        EPETRA_CHK_ERR( Mat[k]->ExtractGlobalRowCopy( Row, MaxNumEntries, NumEntries, Values, Indices));
   
        if( scalar[k] != 1.0 )
           for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalar[k];
   
        if(C->Filled()) { // Sum in values
           err = C->SumIntoGlobalValues( Row, NumEntries, Values, Indices );
           if (err < 0) ierr = err;
        } else { // just add it to the unfilled CRS Matrix
           err = C->InsertGlobalValues( Row, NumEntries, Values, Indices );
           if (err < 0) ierr = err;
        }
     }

     delete [] Indices;
     delete [] Values;
  }

  if( Atrans ) delete Atrans;
  if( Btrans ) delete Btrans;

  return(ierr);
}

*/

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_MATRIXMATRIX_INSTANT(SCALAR,LO,GO,NODE,SPMATVEC,SPMATSLV) \
  \
  template<> \
  int MatrixMatrix::Multiply( \
  Teuchos::RCP<const CrsMatrix< SCALAR , LO , GO , NODE , SPMATVEC , SPMATSLV > >& A, \
  bool transposeA, \
  Teuchos::RCP<const CrsMatrix< SCALAR , LO , GO , NODE , SPMATVEC , SPMATSLV > >& B, \
  bool transposeB, \
  Teuchos::RCP<CrsMatrix< SCALAR , LO , GO , NODE , SPMATVEC , SPMATSLV > >& C, \
  bool call_FillComplete_on_result)
}

#endif // TPETRA_MATRIXMATRIX_DEF_HPP

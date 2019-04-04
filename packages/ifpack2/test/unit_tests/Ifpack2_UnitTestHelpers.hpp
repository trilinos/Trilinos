/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
//@HEADER
*/

#ifndef IFPACK2_UNITTESTHELPERS_HPP
#define IFPACK2_UNITTESTHELPERS_HPP

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Ifpack2_Details_RowMatrix.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace tif_utest {
using Tpetra::global_size_t;

typedef Tpetra::Map<>::node_type Node;

inline
Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
  return Tpetra::getDefaultComm();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
create_tpetra_map (const LocalOrdinal num_elements_per_proc)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  RCP<const map_type> tmap =
    rcp (new map_type (INVALID, num_elements_per_proc, indexBase, comm));

  return tmap;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
create_tridiag_graph (const LocalOrdinal num_rows_per_proc)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  RCP<const map_type> rowmap =
    rcp (new map_type (INVALID, num_rows_per_proc, indexBase, comm));

  RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph =
    rcp (new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> (rowmap, 0));

  Teuchos::Array<GlobalOrdinal> cols;
  const LocalOrdinal lclNumRows =
    static_cast<LocalOrdinal> (rowmap->getNodeNumElements ());
  for (LocalOrdinal l_row = 0; l_row < lclNumRows; ++l_row) {
    GlobalOrdinal g_row = rowmap->getGlobalElement(l_row);
    if (g_row == rowmap->getMinAllGlobalIndex() ||
        g_row == rowmap->getMaxAllGlobalIndex()) {
      cols.resize(2);
    }
    else {
      cols.resize(3);
    }

    size_t coloffset = 0;
    if (g_row > rowmap->getMinAllGlobalIndex()) {
      cols[coloffset++] = g_row-1;
    }
    cols[coloffset++] = g_row;
    if (g_row < rowmap->getMaxAllGlobalIndex()) {
      cols[coloffset++] = g_row+1;
    }
    crsgraph->insertGlobalIndices(g_row, cols());
  }

  crsgraph->fillComplete();

  return crsgraph;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
create_test_graph (const LocalOrdinal num_rows_per_proc)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_rows_per_proc, indexBase, comm));

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowmap, 0));

  size_t global_size = rowmap->getGlobalNumElements();
  Teuchos::Array<GlobalOrdinal> cols;
  for(LocalOrdinal l_row = 0; (size_t) l_row<rowmap->getNodeNumElements(); l_row++) {
    GlobalOrdinal g_row = rowmap->getGlobalElement(l_row);
    if (g_row == rowmap->getMinAllGlobalIndex()) {
      cols.resize(global_size);
      for(size_t i=0; i<global_size; ++i) {
        GlobalOrdinal gcol = g_row + i;
        cols[i] = gcol;
      }
    }
    else {
      cols.resize(2);
      cols[0] = rowmap->getMinAllGlobalIndex();
      cols[1] = g_row;
    }

    crsgraph->insertGlobalIndices(g_row, cols());
  }

  crsgraph->fillComplete();

  return crsgraph;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > create_banded_graph(LocalOrdinal num_rows_per_proc, size_t rbandwidth)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_rows_per_proc, indexBase, comm));

  const size_t maxNumEntPerRow = 5;
  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowmap, maxNumEntPerRow, Tpetra::StaticProfile));

  // Fill the graph with only diagonal entries
  Teuchos::Array<GlobalOrdinal> gblColInds (maxNumEntPerRow);

  for (LocalOrdinal lclRowInd = rowmap->getMinLocalIndex ();
       lclRowInd <= rowmap->getMaxLocalIndex (); ++lclRowInd) {
    const GlobalOrdinal gblRowInd = rowmap->getGlobalElement (lclRowInd);
    size_t jvalid = 0;
    for (size_t i = 0; i < 3; ++i)
    {
      LocalOrdinal lclCol = lclRowInd - 1 + i;
      if (lclCol >= 0 && lclCol < num_rows_per_proc)
      {
        gblColInds[jvalid] = rowmap->getGlobalElement (lclCol);
        jvalid += 1;
      }
    }
    for (size_t i = 0; i < 2; ++i)
    {
      LocalOrdinal lclCol = lclRowInd - rbandwidth + i*(2*rbandwidth);
      if (lclCol >= 0 && lclCol < num_rows_per_proc)
      {
        gblColInds[jvalid] = rowmap->getGlobalElement (lclCol);
        jvalid += 1;
      }

    }
    crsgraph->insertGlobalIndices (gblRowInd, gblColInds.view(0, jvalid));
  }
  crsgraph->fillComplete ();

  return crsgraph;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > create_diagonal_graph(LocalOrdinal num_rows_per_proc)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_rows_per_proc, indexBase, comm));

  const size_t maxNumEntPerRow = 1;
  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowmap, maxNumEntPerRow, Tpetra::StaticProfile));

  // Fill the graph with only diagonal entries
  Teuchos::Array<GlobalOrdinal> gblColInds (maxNumEntPerRow);
  for (LocalOrdinal lclRowInd = rowmap->getMinLocalIndex ();
       lclRowInd <= rowmap->getMaxLocalIndex (); ++lclRowInd) {
    const GlobalOrdinal gblRowInd = rowmap->getGlobalElement (lclRowInd);
    gblColInds[0] = gblRowInd;
    crsgraph->insertGlobalIndices (gblRowInd, gblColInds ());
  }
  crsgraph->fillComplete ();

  return crsgraph;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > create_dense_local_graph(LocalOrdinal num_rows_per_proc)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  const LocalOrdinal indexBase = 0;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(INVALID, num_rows_per_proc, indexBase, comm));

  const size_t maxNumEntPerRow = num_rows_per_proc;
  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowmap, maxNumEntPerRow, Tpetra::StaticProfile));

  // Fill the graph with only diagonal entries
  Teuchos::Array<GlobalOrdinal> gblColInds (maxNumEntPerRow);
  for (LocalOrdinal lclRowInd = rowmap->getMinLocalIndex ();
       lclRowInd <= rowmap->getMaxLocalIndex (); ++lclRowInd) {
    const GlobalOrdinal gblRowInd = rowmap->getGlobalElement (lclRowInd);
    for (size_t k = 0; k < maxNumEntPerRow; ++k)
    {
      LocalOrdinal lclColInd = rowmap->getMinLocalIndex() + k;
      const GlobalOrdinal gblColInd = rowmap->getGlobalElement(lclColInd);
      gblColInds[k] = gblColInd;
    }
    crsgraph->insertGlobalIndices(gblRowInd,gblColInds());
  }
  crsgraph->fillComplete ();

  return crsgraph;
}


// ///////////////////////////////////////////////////////////////////////

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_test_matrix(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap, Scalar val=Teuchos::ScalarTraits<Scalar>::zero())
{
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 3/*tri-diagonal matrix*/));

  Teuchos::Array<GlobalOrdinal> col(3);
  Teuchos::Array<Scalar> coef(3);

  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two = one + one;

  size_t upper_limit;

  for(LocalOrdinal l_row = 0; (size_t) l_row<rowmap->getNodeNumElements(); l_row++) {
    GlobalOrdinal g_row = rowmap->getGlobalElement(l_row);

    if (g_row == rowmap->getMinAllGlobalIndex()) {
      upper_limit=2;
      col[0] = g_row;
      col[1] = g_row+1;
      coef[0] = two;
      coef[1] = val;
    }
    else if (g_row == rowmap->getMaxAllGlobalIndex()) {
      upper_limit=2;
      col[0] = g_row-1;
      col[1] = g_row;
      coef[0] = val;
      coef[1] = two;
    }
    else {
      upper_limit=3;
      col[0] = g_row-1;
      col[1] = g_row;
      col[2] = g_row+1;
      coef[0] = val;
      coef[1] = two;
      coef[2] = val;
    }

    crsmatrix->insertGlobalValues(g_row, col(0,upper_limit), coef(0,upper_limit) );
  }

  crsmatrix->fillComplete();
  return crsmatrix;
} //create_test_matrix

// ///////////////////////////////////////////////////////////////////////

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_test_matrix2(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap)
{
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 3/*tri-diagonal matrix*/));

  Teuchos::Array<GlobalOrdinal> col(1);
  Teuchos::Array<Scalar> coef(1);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  magnitude mag_one = Teuchos::ScalarTraits<Scalar>::magnitude(one);
  magnitude mag_two = mag_one*2.0;
  magnitude mag_ten = mag_one*10.0;
  Scalar two = one*mag_two;
  const Scalar onetenth = one / mag_ten;

  for(LocalOrdinal l_row = 0; (size_t) l_row<rowmap->getNodeNumElements(); l_row++) {
    GlobalOrdinal g_row = rowmap->getGlobalElement(l_row);
    if (g_row == rowmap->getMinAllGlobalIndex()) {
      col.resize(2);
      coef.resize(2);
      col[0] = g_row;
      col[1] = g_row+1;
      coef[0] = two;
      coef[1] = onetenth;
    }
    else if (g_row == rowmap->getMaxAllGlobalIndex()) {
      col.resize(2);
      coef.resize(2);
      col[0] = g_row-1;
      col[1] = g_row;
      coef[0] = onetenth;
      coef[1] = two;
    }
    else {
      col.resize(3);
      coef.resize(3);
      col[0] = g_row-1;
      col[1] = g_row;
      col[2] = g_row+1;
      coef[0] = onetenth;
      coef[1] = two;
      coef[2] = onetenth;
    }

    crsmatrix->insertGlobalValues(g_row, col(), coef() );
  }

  crsmatrix->fillComplete();

  return crsmatrix;
} //create_test_matrix2

// ///////////////////////////////////////////////////////////////////////

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_test_matrix3(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap)
{
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 3/*tri-diagonal matrix*/));

  /*
     NOTE:  this utility creates a matrix whose column map is equal to its row map.  At processor boundaries,
     the row stencil is truncated so that it only contains local columns.  This is done by using
     getMinGlobalIndex() & getMinGlobalIndex(), rather than getMinAllGlobalIndex() & getMaxAllGlobalIndex().
     "
  */

  Teuchos::Array<GlobalOrdinal> col(1);
  Teuchos::Array<Scalar> coef(1);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  magnitude mag_one = Teuchos::ScalarTraits<Scalar>::magnitude(one);
  magnitude mag_two = mag_one*2.0;
  magnitude mag_ten = mag_one*10.0;
  Scalar two = one*mag_two;
  const Scalar onetenth = one / mag_ten;

  for(LocalOrdinal l_row = 0; (size_t) l_row<rowmap->getNodeNumElements(); l_row++) {
    GlobalOrdinal g_row = rowmap->getGlobalElement(l_row);
    if (g_row == rowmap->getMinGlobalIndex()) {
      col.resize(2);
      coef.resize(2);
      col[0] = g_row;
      col[1] = g_row+1;
      coef[0] = two;
      coef[1] = onetenth;
    }
    else if (g_row == rowmap->getMaxGlobalIndex()) {
      col.resize(2);
      coef.resize(2);
      col[0] = g_row-1;
      col[1] = g_row;
      coef[0] = onetenth;
      coef[1] = two;
    }
    else {
      col.resize(3);
      coef.resize(3);
      col[0] = g_row-1;
      col[1] = g_row;
      col[2] = g_row+1;
      coef[0] = onetenth;
      coef[1] = two;
      coef[2] = onetenth;
    }

    crsmatrix->insertGlobalValues(g_row, col(), coef() );
  }

  crsmatrix->fillComplete();

  return crsmatrix;
} //create_test_matrix3

// ///////////////////////////////////////////////////////////////////////

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_banded_matrix(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap, const GlobalOrdinal bw)
{
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 5));

  Teuchos::Array<GlobalOrdinal> col(1);
  Teuchos::Array<Scalar> coef(1);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  magnitude mag_one = Teuchos::ScalarTraits<Scalar>::magnitude(one);
  magnitude mag_four = mag_one*4.0;
  Scalar four = one*mag_four;

  size_t ne;
  for(LocalOrdinal l_row = 0; (size_t) l_row<rowmap->getNodeNumElements(); l_row++) {
    GlobalOrdinal g_row = rowmap->getGlobalElement(l_row);
    if (g_row == rowmap->getMinGlobalIndex()) {
      ne=2;
      col.resize(2);
      coef.resize(2);
      col[0] = g_row;
      col[1] = g_row+1;
      coef[0] = four;
      coef[1] = -one;
    }
    else if (g_row == rowmap->getMaxGlobalIndex()) {
      ne=2;
      col.resize(2);
      coef.resize(2);
      col[0] = g_row-1;
      col[1] = g_row;
      coef[0] = -one;
      coef[1] = four;
    }
    else {
      ne=3;
      col.resize(3);
      coef.resize(3);
      col[0] = g_row-1;
      col[1] = g_row;
      col[2] = g_row+1;
      coef[0] = -one;
      coef[1] = four;
      coef[2] = -one;
    }

    // upper band
    if (g_row <= rowmap->getMaxGlobalIndex() - bw) {
      ne++;
      col.resize(ne);
      coef.resize(ne);
      col[ne-1] = g_row+bw;
      coef[ne-1] = -one;
    }

    // lower band
    if (g_row >= rowmap->getMinGlobalIndex() + bw) {
      ne++;
      col.resize(ne);
      coef.resize(ne);
      col[ne-1] = g_row-bw;
      coef[ne-1] = -one;
    }

    crsmatrix->insertGlobalValues(g_row, col(), coef() );
  }

  crsmatrix->fillComplete();

  return crsmatrix;
} //create_banded_matrix

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  create_block_diagonal_matrix(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& graph, const int blockSize)
{

  Teuchos::RCP<Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bcrsmatrix =
      Teuchos::rcp(new Tpetra::Experimental::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*graph, blockSize));

  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& meshRowMap = *bcrsmatrix->getRowMap();

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two = one+one;
  const Scalar three = two+one;

  Teuchos::Array<Scalar> basematrix(blockSize*blockSize, zero);
  basematrix[0] = two;
  basematrix[2] = three;
  basematrix[3] = three;
  basematrix[4] = two;
  basematrix[7] = three;
  basematrix[8] = two;

  Teuchos::Array<LocalOrdinal> lclColInds(1);
  for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex (); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd)
  {
    lclColInds[0] = lclRowInd;
    bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &basematrix[0], 1);
  }

  return bcrsmatrix;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,bool LowerTriangular>
Teuchos::RCP<const Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  create_triangular_matrix(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& graph, const int blockSize)
{

  Teuchos::RCP<Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bcrsmatrix =
      Teuchos::rcp(new Tpetra::Experimental::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*graph, blockSize));

  const LocalOrdinal numLocalRows = bcrsmatrix->getNodeNumRows();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& meshRowMap = *bcrsmatrix->getRowMap();

  Teuchos::Array<GlobalOrdinal> col(1);
  Teuchos::Array<Scalar> coef(1);

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two = one+one;
  const Scalar three = one+one+one;

  Teuchos::Array<Scalar> basematrix(numLocalRows*numLocalRows, zero);
  basematrix[0] = two;
  basematrix[3] = three;
  basematrix[4] = two;
  basematrix[6] = one;
  basematrix[7] = three;
  basematrix[8] = two;

  Teuchos::Array<Scalar> assembleMatrix(blockSize*blockSize, zero);

  Teuchos::Array<LocalOrdinal> lclColInds(1);

  for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex (); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd)
  {
    const LocalOrdinal rowOffset = lclRowInd-meshRowMap.getMinLocalIndex();

    for (LocalOrdinal i = 0; i < numLocalRows; ++i)
    {

      for (int j = 0; j < blockSize*blockSize; ++j)
        assembleMatrix[j] = zero;
      size_t indexBaseMatrix = 0;
      if (LowerTriangular)
        indexBaseMatrix = numLocalRows*rowOffset+i;
      else
        indexBaseMatrix = i*numLocalRows+rowOffset;

      for (int j = 0; j < blockSize; ++j)
        assembleMatrix[(blockSize+1)*j] = basematrix[indexBaseMatrix];

      lclColInds[0] = meshRowMap.getMinLocalIndex () + i;
      bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &assembleMatrix[0], 1);

    }
  }

  return bcrsmatrix;
}

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  create_banded_block_matrix(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& graph, const int blockSize, const size_t rbandwidth)
{
  Teuchos::RCP<Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bcrsmatrix
    = Teuchos::rcp(new Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(*graph, blockSize));
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& meshRowMap = *bcrsmatrix->getRowMap();

  const int blockMatSize = blockSize*blockSize;
  const LocalOrdinal maxNumEntries = 2*rbandwidth+1;
  Teuchos::Array<GlobalOrdinal> col(maxNumEntries);
  Teuchos::Array<Scalar> coef(maxNumEntries*blockMatSize);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  magnitude mag_one = Teuchos::ScalarTraits<Scalar>::magnitude(one);
  magnitude mag_four = mag_one*4.0;
  Scalar four = one*mag_four;

  for(LocalOrdinal l_row = 0; (size_t) l_row < meshRowMap.getNodeNumElements(); ++l_row)
  {

    const LocalOrdinal * inds;
    Scalar * vals;
    LocalOrdinal numInd;
    bcrsmatrix->getLocalRowView(l_row, inds, vals, numInd);
    for (LocalOrdinal j = 0; j < numInd; ++j)
    {
      const LocalOrdinal lcl_col = inds[j];
      if (lcl_col == l_row)
        vals[j] = four;
      else
        vals[j] = -one;
    }
  }

  return bcrsmatrix;
}
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  create_full_local_matrix(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& graph, const int blockSize)
{

  Teuchos::RCP<Tpetra::Experimental::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bcrsmatrix =
      Teuchos::rcp(new Tpetra::Experimental::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*graph, blockSize));

  const LocalOrdinal numLocalRows = bcrsmatrix->getNodeNumRows();
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>& meshRowMap = *bcrsmatrix->getRowMap();

  Teuchos::Array<GlobalOrdinal> col(1);
  Teuchos::Array<Scalar> coef(1);

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two = one+one;
  const Scalar three = one+one+one;

  Teuchos::Array<Scalar> basematrix(numLocalRows*numLocalRows, zero);
  basematrix[0] = two*three;
  basematrix[1] = three;
  basematrix[2] = one;

  basematrix[3] = three*two;
  basematrix[4] = three*two;
  basematrix[5] = three;

  basematrix[6] = two;
  basematrix[7] = two*three;
  basematrix[8] = two*three;

  Teuchos::Array<Scalar> assembleMatrix(blockSize*blockSize, zero);

  Teuchos::Array<LocalOrdinal> lclColInds(1);

  for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex (); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd)
  {
    const LocalOrdinal rowOffset = lclRowInd-meshRowMap.getMinLocalIndex();

    for (LocalOrdinal i = 0; i < numLocalRows; ++i)
    {

      for (int j = 0; j < blockSize*blockSize; ++j)
        assembleMatrix[j] = zero;
      const size_t indexBaseMatrix = numLocalRows*rowOffset+i;

      for (int j = 0; j < blockSize; ++j)
        assembleMatrix[blockSize*j+blockSize-1-j] = basematrix[indexBaseMatrix];

      lclColInds[0] = meshRowMap.getMinLocalIndex () + i;
      bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &assembleMatrix[0], 1);

    }
  }

  return bcrsmatrix;
}


// ///////////////////////////////////////////////////////////////////////
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_test_matrix_variable_blocking(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap) {
  // Basically each processor gets this 5x5 block lower-triangular matrix:
  //
  // [ 2 -1  0  0  0 ;...
  // [-1  2  0  0  0 ;...
  // [ 0 -1  3 -1  0 ;...
  // [ 0  0 -1  3 -1 ;...
  // [ 0  0  0 -1  2  ];
  //
  // Beyond this the matrix is diagonal.  We also explicitly stick in some hard zeros for the line partitioning test

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 0));

  Teuchos::Array<GlobalOrdinal> indices(4);
  Teuchos::Array<Scalar> values(4);
  GlobalOrdinal rb = rowmap->getGlobalElement(0);

  /*** Fill Matrix ****/
  // Row 0
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2;
  values[0] =2;  values[1] =-1;   values[2] = 0;
  crsmatrix->insertGlobalValues(rb, indices(0,3), values(0,3));

  // Row 1
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2;
  values[0] =-1; values[1] =2;    values[2] = 0;
  crsmatrix->insertGlobalValues(rb+1, indices(0,3), values(0,3));

  // Row 2
  indices[0]=rb+1; indices[1]=rb+2; indices[2]=rb+3; indices[3]=rb+4;
  values[0] =-1;   values[1] = 3;   values[2] =-1;   values[3] = 0;
  crsmatrix->insertGlobalValues(rb+2, indices(0,4), values(0,4));

  // Row 3
  indices[0]=rb+2; indices[1]=rb+3; indices[2]=rb+4; indices[3]=rb+1;
  values[0] =-1;   values[1] = 3;   values[2] =-1;   values[3] =0.0;
  crsmatrix->insertGlobalValues(rb+3, indices(0,4), values(0,4));

  // Row 4
  indices[0]=rb+3; indices[1]=rb+4;
  values[0] =-1;   values[1] = 2;
  crsmatrix->insertGlobalValues(rb+4, indices(0,2), values(0,2));

  // All other rows
  for(size_t i=5; i<rowmap->getNodeNumElements(); i++) {
    indices[0]=rb+i;
    values[0] = 1;
    crsmatrix->insertGlobalValues(rb+i, indices(0,1), values(0,1));
  }
  crsmatrix->fillComplete();

  return crsmatrix;
}

// ///////////////////////////////////////////////////////////////////////
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_test_matrix_variable_banded(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap) {
  // Basically each processor gets this 5x5 block banded matrix:
  //
  // [  3  -1   0  -1   0 ;...
  // [ -1   3   0   0  -1 ;...
  // [ -1  -1   4  -1   0 ;...
  // [  0  -1  -1   4  -1 ;...
  // [  0   0  -1  -1   3 ];
  //
  // Beyond this the matrix is diagonal.

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 0));

  Teuchos::Array<GlobalOrdinal> indices(5);
  Teuchos::Array<Scalar> values(5);
  GlobalOrdinal rb = rowmap->getGlobalElement(0);

  /*** Fill Matrix ****/
  // Row 0
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3;
  values[0] =3;  values[1] =-1;   values[2] = 0;   values[3] = -1;
  crsmatrix->insertGlobalValues(rb, indices(0,4), values(0,4));

  // Row 1
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3 ;indices[4]=rb+4;
  values[0] =-1; values[1] =3;    values[2] = 0;   values[3] = 0;   values[4] = -1;
  crsmatrix->insertGlobalValues(rb+1, indices(0,5), values(0,5));

  // Row 2
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3; indices[4]=rb+4;
  values[0] =-1;   values[1] =-1;  values[2] =4;   values[3] = -1;  values[4] = 0;
  crsmatrix->insertGlobalValues(rb+2, indices(0,5), values(0,5));

  // Row 3
  indices[0]=rb+1; indices[1]=rb+2; indices[2]=rb+3; indices[3]=rb+4;
  values[0] =-1;   values[1] = -1;   values[2]=4;   values[3]  =-1;
  crsmatrix->insertGlobalValues(rb+3, indices(0,4), values(0,4));

  // Row 4
  indices[0]=rb+2; indices[1]=rb+3; indices[2]=rb+4;
  values[0] =-1;   values[1] = -1;   values[2]= 3;
  crsmatrix->insertGlobalValues(rb+4, indices(0,3), values(0,3));

  // All other rows
  for(size_t i=5; i<rowmap->getNodeNumElements(); i++) {
    indices[0]=rb+i;
    values[0] = 1;
    crsmatrix->insertGlobalValues(rb+i, indices(0,1), values(0,1));
  }
  crsmatrix->fillComplete();

  return crsmatrix;
}

// ///////////////////////////////////////////////////////////////////////
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > create_test_matrix_banded_variable_blocking(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowmap) {
  // Basically each processor gets this 5x5 block lower-triangular matrix:
  //
  // [ 2  -1   0   0   0 ;...
  // [-1   2   0   0   0 ;...
  // [-1  -1   4  -1   0 ;...
  // [ 0  -1  -1   4  -1 ;...
  // [ 0   0  -1  -1   3  ];
  //
  // Beyond this the matrix is diagonal.  We also explicitly stick in some hard zeros for the line partitioning test

  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = Teuchos::rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap, 0));

  Teuchos::Array<GlobalOrdinal> indices(5);
  Teuchos::Array<Scalar> values(5);
  GlobalOrdinal rb = rowmap->getGlobalElement(0);

  /*** Fill Matrix ****/
  // Row 0
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2;
  values[0] =2;  values[1] =-1;   values[2] = 0;
  crsmatrix->insertGlobalValues(rb, indices(0,3), values(0,3));

  // Row 1
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2;
  values[0] =-1; values[1] =2;    values[2] = 0;
  crsmatrix->insertGlobalValues(rb+1, indices(0,3), values(0,3));

  // Row 2
  indices[0]=rb; indices[1]=rb+1; indices[2]=rb+2; indices[3]=rb+3; indices[4]=rb+4;
  values[0] =-1;   values[1] = -1;  values[2] =4;   values[3] = -1; values[4] = 0;
  crsmatrix->insertGlobalValues(rb+2, indices(0,5), values(0,5));

  // Row 3
  indices[0]=rb+2; indices[1]=rb+3; indices[2]=rb+4; indices[3]=rb+1; indices[4]=rb;
  values[0] =-1;   values[1] = 4;   values[2] =-1;   values[3] =-1;  values[4] =0.0;
  crsmatrix->insertGlobalValues(rb+3, indices(0,5), values(0,5));

  // Row 4
  indices[0]=rb+3; indices[1]=rb+4; indices[2]=rb+2;
  values[0] =-1;   values[1] = 3;   values[2] = -1;
  crsmatrix->insertGlobalValues(rb+4, indices(0,3), values(0,3));

  // All other rows
  for(size_t i=5; i<rowmap->getNodeNumElements(); i++) {
    indices[0]=rb+i;
    values[0] = 1;
    crsmatrix->insertGlobalValues(rb+i, indices(0,1), values(0,1));
  }
  crsmatrix->fillComplete();

  return crsmatrix;
}

  template<class Scalar = Tpetra::RowMatrix<>::scalar_type,
           class LocalOrdinal =
             typename Tpetra::RowMatrix<Scalar>::local_ordinal_type,
           class GlobalOrdinal =
             typename Tpetra::RowMatrix<Scalar, LocalOrdinal>::global_ordinal_type,
           class Node =
             typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class NotCrsMatrix :
    public Ifpack2::Details::RowMatrix<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > {
  public:
    NotCrsMatrix (Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A) : A_(A){;}
    virtual ~NotCrsMatrix(){;}
    virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {return A_->getComm();}
    virtual Teuchos::RCP<Node> getNode() const {return A_->getNode();}
    virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRangeMap () const {return A_->getRangeMap();}
    virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getDomainMap () const {return A_->getDomainMap();}
    virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRowMap () const {return A_->getRowMap();}
    virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getColMap () const {return A_->getColMap();}
    virtual Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >
    getGraph () const {return A_->getGraph();}
    virtual global_size_t getGlobalNumRows() const {return A_->getGlobalNumRows();}
    virtual global_size_t getGlobalNumCols() const {return A_->getGlobalNumCols();}
    virtual size_t getNodeNumRows() const {return A_->getNodeNumRows();}
    virtual size_t getNodeNumCols() const {return A_->getNodeNumCols();}
    virtual GlobalOrdinal getIndexBase() const {return A_->getIndexBase();}
    virtual global_size_t getGlobalNumEntries() const {return A_->getGlobalNumEntries();}
    virtual size_t getNodeNumEntries() const {return A_->getNodeNumEntries();}
    virtual size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const {return A_->getNumEntriesInGlobalRow(globalRow);}
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {return A_->getNumEntriesInLocalRow(localRow);}
    virtual size_t getGlobalMaxNumRowEntries() const {return A_->getGlobalMaxNumRowEntries();}
    virtual size_t getNodeMaxNumRowEntries() const {return A_->getNodeMaxNumRowEntries();}
    virtual bool hasColMap() const {return A_->hasColMap();}
    virtual bool isLocallyIndexed() const {return A_->isLocallyIndexed();}
    virtual bool isGloballyIndexed() const {return A_->isGloballyIndexed();}
    virtual bool isFillComplete() const {return A_->isFillComplete();}
    virtual bool supportsRowViews() const {return A_->supportsRowViews();}

    virtual void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                      const Teuchos::ArrayView<Scalar> &Values,
                      size_t &NumEntries) const {A_->getGlobalRowCopy(GlobalRow,Indices,Values,NumEntries);}

    virtual void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const Teuchos::ArrayView<LocalOrdinal> &Indices,
                     const Teuchos::ArrayView<Scalar> &Values,
                     size_t &NumEntries) const {A_->getLocalRowCopy(LocalRow,Indices,Values,NumEntries);}
    virtual void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      Teuchos::ArrayView<const GlobalOrdinal> &indices,
                      Teuchos::ArrayView<const Scalar> &values) const {A_->getGlobalRowView(GlobalRow,indices,values);}

    virtual void
    getLocalRowView (LocalOrdinal LocalRow,
                     Teuchos::ArrayView<const LocalOrdinal> &indices,
                     Teuchos::ArrayView<const Scalar> &values) const {A_->getLocalRowView(LocalRow,indices,values);}

    virtual void getLocalDiagCopy (Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const {A_->getLocalDiagCopy(diag);}
    virtual void leftScale (const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {A_->leftScale(x);}
    virtual void rightScale (const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {A_->rightScale(x);}
    virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {return A_->getFrobeniusNorm();}

    virtual Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    add (const Scalar& alpha,
         const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Scalar& beta,
         const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap=Teuchos::null,
         const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap=Teuchos::null,
         const Teuchos::RCP<Teuchos::ParameterList>& params=Teuchos::null) const
    {return A_->add(alpha,A,beta,domainMap,rangeMap,params);}

    virtual void
    apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
    {A_->apply(X,Y,mode,alpha,beta);}

    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<char>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Tpetra::Distributor& distor) const
    {A_->pack(exportLIDs,exports,numPacketsPerLID,constantNumPackets,distor);}

   private:
   Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  }; // class NotCrsMatrix



}//namespace tif_utest

#endif


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

#ifndef TPETRA_MATRIXMATRIX_DECL_HPP
#define TPETRA_MATRIXMATRIX_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MMHelpers.hpp"

/*! \file Tpetra_MatrixMatrix_decl.hpp 

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Tpetra {

// CGB: you included the header file above. you don't need a forward declaration.
// #ifndef DOXYGEN_SHOULD_SKIP_THIS  
//   // forward declaration
// template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
// class CrsMatrix;
// #endif

  /** Collection of matrix-matrix operations. This class basically
      functions as a namespace, containing only static methods.
      See the program epetraext/test/MatrixMatrix/cxx_main.cpp for
      a usage example.
   */
  
namespace MatrixMatrix {


    /** Given CrsMatrix objects A, B and C, form the product C = A*B.
  In a parallel setting, A and B need not have matching distributions,
  but C needs to have the same row-map as A.
  At this time C=AT*B and C=A*BT are known to not work. However,
  C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on C or not. If it has,
       then C's graph must already contain all nonzero locations that
       will be produced when forming the product A*B. On exit,
       C.FillComplete() will have been called, unless the last argument
             to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

     */
template <class Scalar, 
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class SpMatOps >
void Multiply(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
  bool transposeA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
  bool transposeB,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& C,
  bool call_FillComplete_on_result=true);

    /** Given CrsMatrix objects A and B, form the sum B = a*A + b*B
     * Currently not functional.

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on B or not. If it has,
       then B's graph must already contain all nonzero locations that
       will be produced when forming the sum.
    @param scalarB Input, scalar multiplier for matrix B.

     */
template <class Scalar, 
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class SpMatOps >
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
  bool transposeA,
  Scalar scalarA,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
  Scalar scalarB );

    /** Given CrsMatrix objects A and B, form the sum C = a*A + b*B

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param scalarB Input, scalar multiplier for matrix B.
    @param C Result. On entry to this method, C can be NULL or a pointer
             to an unfilled or filled matrix. If C is NULL then a new
             object is allocated and must be deleted by the user.
             If C is not NULL and FillComplete has already been called then the sparsity pattern is assumed to be fixed and compatible  with the sparsity of A+B. If FillComplete has not been called then the sum is completed and the function
             returns without calling FillComplete on C.

     */
template <class Scalar, 
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class SpMatOps>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
  bool transposeA,
  Scalar scalarA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
  bool transposeB,
  Scalar scalarB,
  RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> > C);
} //End Namespace MatrixMatrxix

namespace MMdetails{

template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node,
         class SpMatOps>
void
getGlobalRowFromLocalIndex(
  LocalOrdinal localRow,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Mview, 
  Array<GlobalOrdinal>& indices,
  ArrayView<Scalar>& values);

template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > find_rows_containing_cols(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& M,
  const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > colmap);


template<class Scalar, 
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node, 
         class SpMatOps>
void mult_A_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Aview, 
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Bview, 
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C);
  
template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node,
         class SpMatOps>
void mult_A_Btrans(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node> & C);

template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node,
         class SpMatOps>
void mult_Atrans_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Aview, 
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node> & C);

template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node,
         class SpMatOps>
void mult_Atrans_Btrans(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Aview, 
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C);

template<class Scalar,
         class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node,
         class SpMatOps>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& M,
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Mview);

template<class Ordinal,
         class GlobalOrdinal>
void distribute_list(
  RCP<const Comm<Ordinal> > comm,
  size_t lenSendList,
  const Array<GlobalOrdinal>& sendList,
  size_t& maxSendLen,
  Array<GlobalOrdinal>& recvList);

template<class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node>
RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > create_map_from_imported_rows(
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
  const size_t& totalNumSend,
  ArrayView<GlobalOrdinal> sendRows,
  const int& numProcs,
  ArrayView<size_t> numSendPerProc);

template<class LocalOrdinal, 
         class GlobalOrdinal, 
         class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > form_map_union(
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > map1,
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > map2);

  /**
   * Method for internal use... sparsedot forms a dot-product between two
   * sparsely-populated 'vectors'.
   * Important assumption: assumes the indices in u_ind and v_ind are sorted.
   */
  template <class Scalar, class GlobalOrdinal>
  Scalar sparsedot(
    const Teuchos::ArrayView<Scalar> &u, 
    const Teuchos::ArrayView<GlobalOrdinal> &u_ind, 
    const Teuchos::ArrayView<Scalar> &v, 
    const Teuchos::ArrayView<GlobalOrdinal> &v_ind);
  
}//end namespace MMdetails




  

} // end of Tpetra namespace

#endif // TPETRA_MATRIXMATRIX_DECL_HPP


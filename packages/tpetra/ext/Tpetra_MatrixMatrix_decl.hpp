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
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MMHelpers.hpp"

/*! \file Tpetra_MMMultiply_decl.hpp 

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
  
template <
  class Scalar, 
  class LocalOrdinal  = int, 
  class GlobalOrdinal = LocalOrdinal, 
  class Node = Kokkos::DefaultNode::DefaultNodeType, 
  class SpMatOps = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >
class MatrixMatrix {
public:

    typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map_t;
    typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node, SpMatOps> 
      CrsGraph_t;
    typedef CrsMatrix<
      Scalar, 
      LocalOrdinal,
      GlobalOrdinal,
      Node,
      SpMatOps> CrsMatrix_t;

    typedef CrsMatrixStruct<
      Scalar, 
      LocalOrdinal,
      GlobalOrdinal,
      Node,
      SpMatOps> CrsMatrixStruct_t;

    /** destructor */
  virtual ~MatrixMatrix(){}

    /** Given CrsMatrix objects A, B and C, form the product C = A*B.
  In a parallel setting, A and B need not have matching distributions,
  but C needs to have the same row-map as A.

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

    @return error-code, 0 if successful. non-zero returns may result if A or
             B are not already Filled, or if errors occur in putting values
             into C, etc.
     */
static int Multiply(
  Teuchos::RCP<const CrsMatrix_t > A,
  bool transposeA,
  Teuchos::RCP<const CrsMatrix_t > B,
  bool transposeB,
  Teuchos::RCP<CrsMatrix_t > C,
  bool call_FillComplete_on_result=true);

    /** Given CrsMatrix objects A and B, form the sum B = a*A + b*B

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on B or not. If it has,
       then B's graph must already contain all nonzero locations that
       will be produced when forming the sum.
    @param scalarB Input, scalar multiplier for matrix B.

    @return error-code, 0 if successful. non-zero returns may result if A is
             not already Filled, or if errors occur in putting values
             into B, etc.
     */
  static int Add(
    Teuchos::RCP<CrsMatrix_t > A,
    bool transposeA,
    Scalar scalarA,
    Teuchos::RCP<CrsMatrix_t > B,
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

    @return error-code, 0 if successful. non-zero returns may result if A or is
             not already Filled, or if errors occur in putting values
             into C, etc.
     */
  static int Add(
    Teuchos::RCP<const CrsMatrix_t > A,
    bool transposeA,
    Scalar scalarA,
    Teuchos::RCP<const CrsMatrix_t > B,
    bool transposeB,
    Scalar scalarB,
    Teuchos::RCP<CrsMatrix_t > C);

  static Teuchos::RCP<const Map_t >
  find_rows_containing_cols(
    Teuchos::RCP<const CrsMatrix_t > M,
    Teuchos::RCP<const Map_t > colmap);

private:

static int mult_A_B(
  Teuchos::RCP<CrsMatrixStruct_t >& Aview, 
  Teuchos::RCP<CrsMatrixStruct_t >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C);
  
static void mult_A_Btrans(
  const Teuchos::RCP<const CrsMatrixStruct_t> & Aview, 
  const Teuchos::RCP<const CrsMatrixStruct_t> & Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C);

static int mult_Atrans_B(
  Teuchos::RCP<CrsMatrixStruct_t >& Aview, 
  Teuchos::RCP<CrsMatrixStruct_t >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>&  C);

static int mult_Atrans_Btrans(
  Teuchos::RCP<CrsMatrixStruct_t >& Aview, 
  Teuchos::RCP<CrsMatrixStruct_t >& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C);

static int import_and_extract_views(
  Teuchos::RCP<const CrsMatrix_t >& M,
  Teuchos::RCP<const Map_t >& targetMap,
  Teuchos::RCP<CrsMatrixStruct_t >& Mview);

template<class Ordinal>
static int distribute_list(
  const Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm,
  size_t lenSendList,
  const Teuchos::Array<GlobalOrdinal>& sendList,
  size_t& maxSendLen,
  Teuchos::Array<GlobalOrdinal>& recvList);

static
Teuchos::RCP<const Map_t > 
create_map_from_imported_rows(
  Teuchos::RCP<const Map_t > map,
  size_t totalNumSend,
  const Teuchos::ArrayView<const GlobalOrdinal> &sendRows,
  int numProcs,
  const Teuchos::ArrayView<const size_t> &numSendPerProc);

static 
Teuchos::RCP<const Map_t > 
form_map_union(
  Teuchos::RCP<const Map_t > map1,
  Teuchos::RCP<const Map_t > map2);
  
};//class MatrixMatrix


// CGB: If it isn't for public consumption, then hide it in a different namespace.
namespace MMdebug {

    RCP<Teuchos::FancyOStream> debug_stream;
    Teuchos::EVerbosityLevel   debug_level;

}

namespace MMdetails {

  /**
   * Method for internal use... sparsedot forms a dot-product between two
   * sparsely-populated 'vectors'.
   * Important assumption: assumes the indices in u_ind and v_ind are sorted.
   */
  template <class Scalar, class LocalOrdinal>
  Scalar sparsedot(const Teuchos::ArrayView<const Scalar> &u, const Teuchos::ArrayView<const LocalOrdinal> &u_ind, 
                   const Teuchos::ArrayView<const Scalar> &v, const Teuchos::ArrayView<const LocalOrdinal> &v_ind);
  
} // end of MMdetails namespace

} // end of Tpetra namespace

#endif // TPETRA_MATRIXMATRIX_DECL_HPP


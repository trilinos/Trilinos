// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_DECL_HPP
#define TPETRA_MATRIXMATRIX_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "TpetraExt_MMHelpers.hpp"

/*! \file Tpetra_MatrixMatrix_decl.hpp

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Tpetra {

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


/// \brief Compute the sparse matrix sum <tt>C = scalarA * Op(A) +
///   scalarB * Op(B)</tt>, where Op(X) is either X or its transpose.
///
/// \pre If A and B are both fill complete, then they must have the
///   same domain and range Maps.
///
/// \param scalarA [in] Scalar multiplier for A in the sum.
/// \param transposeA [in] If true, use the transpose of A.
/// \param A [in] The first input matrix.
///
/// \param scalarB [in] Scalar multiplier for B in the sum.
/// \param transposeB [in] If true, use the transpose of B.
/// \param B [in] The second input matrix.
///
/// \param domainMap [in] Domain Map of C.
/// \param rangeMap [in] Range Map of C.
/// \param params [in/out] Same as the parameters of RowMatrix::add.
///
/// For \c domainMap and \c rangeMap, see the documentation of
/// RowMatrix::add.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class SpMatOps >
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >
add (const Scalar& alpha,
     const bool transposeA,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& A,
     const Scalar& beta,
     const bool transposeB,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& B,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap=Teuchos::null,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap=Teuchos::null,
     const Teuchos::RCP<Teuchos::ParameterList>& params=Teuchos::null);

/// \brief Compute the sparse matrix sum <tt>C = scalarA * Op(A) +
///   scalarB * Op(B)</tt>, where Op(X) is either X or its transpose.
///
/// \pre Both input matrices A and B must be fill complete.  That is,
///   their fillComplete() method must have been called at least once,
///   without an intervening call to resumeFill().
///
/// \param A [in] The first input matrix.
/// \param transposeA [in] If true, use the transpose of A.
/// \param scalarA [in] Scalar multiplier for A in the sum.
///
/// \param B [in] The second input matrix.
/// \param transposeB [in] If true, use the transpose of B.
/// \param scalarB [in] Scalar multiplier for B in the sum.
///
/// \param C [in/out] On entry, C may be either null or a valid
///   matrix.  If C is null on input, this function will allocate a
///   new CrsMatrix to contain the sum.  If C is not null and is fill
///   complete, then this function assumes that the sparsity pattern
///   of the sum is fixed and compatible with the sparsity pattern of
///   A + B.  If C is not null and is not fill complete, then this
///   function returns without calling fillComplete on C.
///
/// \warning The case where C == null on input does not actually work.
///   In order for it to work, we would need to change the interface
///   of this function to pass in C as a Ptr<RCP<CrsMatrix> >.
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
} // namespace MatrixMatrix

namespace MMdetails{

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class SpMatOps>
void mult_A_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  bool onlyCalculateStructure=false);

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class SpMatOps>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& M,
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Mview);

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class SpMatOps>
void setMaxNumEntriesPerRow(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>& Mview);

}//end namespace MMdetails

} // end of Tpetra namespace

#endif // TPETRA_MATRIXMATRIX_DECL_HPP


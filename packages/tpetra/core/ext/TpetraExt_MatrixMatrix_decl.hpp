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

#include <string>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "TpetraExt_MMHelpers.hpp"
#include "KokkosKernels_Handle.hpp"


/*! \file TpetraExt_MatrixMatrix_decl.hpp

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Tpetra {

namespace MatrixMatrix {

/// \brief Sparse matrix-matrix multiply
///
/// Given CrsMatrix instances A and B, compute the product C = A*B,
/// overwriting an existing CrsMatrix instance C with the result.
///
/// \pre All three matrices A, B, and C must have uniquely owned row
///   Maps.
/// \pre On input, C must have the same row Map as A.
/// \pre A and B must be fill complete.
/// \pre If C has a range Map on input, then A and C must have the
///   same range Maps.
/// \pre If C has a domain Map on input, then B and C must have the
///   same domain Maps.
///
/// For the latter two preconditions, recall that a matrix does not
/// have a domain or range Map unless fillComplete has been called on
/// it at least once.
///
/// \param A [in] fill-complete sparse matrix.
/// \param transposeA [in] Whether to use transpose of matrix A.
/// \param B [in] fill-complete sparse matrix.
/// \param transposeB [in] Whether to use transpose of matrix B.
/// \param C [in/out] On entry to this method, if C is fill complete,
///   then C's graph must have the correct structure, that is, its
///   structure must equal the structure of A*B.  On exit, C will be
///   fill complete, unless the last argument to this function is
///   false.
/// \param call_FillComplete_on_result [in] Optional argument;
///   defaults to true.  If false, C will <i>not</i> be fill complete
///   on output.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Multiply(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  bool call_FillComplete_on_result = true,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

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
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  Scalar scalarB );

/// \brief Compute the sparse matrix sum <tt>C = scalarA * Op(A) +
///   scalarB * Op(B)</tt>, where Op(X) is either X or its transpose.
///
/// This version of sparse matrix-matrix add returns a new CrsMatrix
/// instance, rather than using an existing instance for the result.
/// The returned matrix is always fill complete, with the given domain and
/// range Maps.  It is correct (though less efficient) for A and B to
/// have different row Maps; the returned matrix will have the same
/// row Map as the row Map of B.
///
/// \pre A and B must both be fillComplete and have matching domain and
///   range Maps.
///
/// \param scalarA [in] Scalar multiplier for A in the sum.
/// \param transposeA [in] If true, use the transpose of A.
/// \param A [in] The first input matrix.
///
/// \param scalarB [in] Scalar multiplier for B in the sum.
/// \param transposeB [in] If true, use the transpose of B.
/// \param B [in] The second input matrix.
///
/// \param domainMap [in] Domain Map of C (on output).  If null or not
///   provided, this defaults to the domain map of Op(B).
/// \param rangeMap [in] Range Map of C (on output).  If null or not
///   provided, this defaults to the range map of Op(B).
/// \param params [in/out] Same as the parameters of RowMatrix::add.
///
/// See the documentation of RowMatrix::add for a more detailed
/// discussion of the optional parameters.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
add (const Scalar& alpha,
     const bool transposeA,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
     const Scalar& beta,
     const bool transposeB,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap = Teuchos::null,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap  = Teuchos::null,
     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


/// \brief Compute the sparse matrix sum <tt>C = scalarA * Op(A) +
///   scalarB * Op(B)</tt>, where Op(X) is either X or its transpose.
///
/// This version of sparse matrix-matrix add returns a new CrsMatrix
/// instance, rather than using an existing instance for the result.
/// The returned matrix is always locally indexed and
/// is fill complete by default (with the given
/// domain and range Maps, or using those of Op(B) if null).
/// If params->get("Call fillComplete") = false,
/// then the resulting matrix will not be fill complete but it will have
/// an immutable sparsity pattern (only scalar values can be changed).
/// This is because C is allocated with the exact amount of storage to hold the sum.
///
/// It is correct (though less efficient) for A and B to
/// have different row Maps; the returned matrix will have the same
/// row Map as the row Map of B.
///
/// \pre A and B must both be fillComplete and have matching domain and
///   range Maps.
///
/// \param scalarA [in] Scalar multiplier for A in the sum.
/// \param transposeA [in] If true, use the transpose of A.
/// \param A [in] The first input matrix.
/// \param scalarB [in] Scalar multiplier for B in the sum.
/// \param transposeB [in] If true, use the transpose of B.
/// \param B [in] The second input matrix.
/// \param C [out] The result matrix, which we expect to be 'new' (no entries inserted) on input.
/// \param domainMap [in] Domain Map of C (on output).  If null or not
///   provided, this defaults to the domain map of B (or range map, if transposeB).
/// \param rangeMap [in] Range Map of C (on output).  If null or not
///   provided, this defaults to the range map of B (or domain map, if transposeB).
/// \param params [in/out] Same as the parameters of RowMatrix::add.
///
/// See the documentation of RowMatrix::add for a more detailed
/// discussion of the optional parameters.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
add (const Scalar& alpha,
     const bool transposeA,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
     const Scalar& beta,
     const bool transposeB,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & C,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap = Teuchos::null,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap  = Teuchos::null,
     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


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
///   of this function (for example, to pass in C as a (pointer or
///   nonconst reference) to a Teuchos::RCP).  Please use add() (which
///   see) if you want matrix-matrix add to return a new instance of
///   CrsMatrix.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  Scalar scalarB,
  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > C);


  /** Given CrsMatrix objects A, B and C, and Vector Dinv, form the product C = (I-omega * Dinv A)*B
      In a parallel setting, A and B need not have matching distributions,
      but C needs to have the same row-map as A.

    @param omega Input, scalar multiplier for Dinverse A
    @param Dinv Input, Vector representing a diagonal matrix, must match A->getRowMap()
    @param A Input, must already have had 'fillComplete()' called.
    @param B Input, must already have had 'fillComplete()' called.
    @param C Result. On entry to this method, it doesn't matter whether
             fillComplete() has already been called on C or not. If it has,
             then C's graph must already contain all nonzero locations that
             will be produced when forming the product A*B. On exit,
             C.FillComplete() will have been called, unless the last argument
             to this function is specified to be false.
    @param call_fillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.fillComplete. (It is often useful
           to allow this function to call C.fillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)
  */
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  void Jacobi(Scalar omega,
              const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
              const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
              const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
              bool call_FillComplete_on_result = true,
                   const std::string& label = std::string(),
              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

} // namespace MatrixMatrix

namespace MMdetails{

  // Matrix multiplication functions
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_AT_B_newmatrix(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B_newmatrix(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);




template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B_reuse(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);




  // Matrix jacobi functions
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void jacobi_A_B_newmatrix(
  Scalar omega,
  const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void jacobi_A_B_reuse(
  Scalar omega,
  const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  // Other functions
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview,
  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal, Node> > prototypeImporter = Teuchos::null,
  bool userAssertsThereAreNoRemotes = false,
  const std::string& label = std::string(),
  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void setMaxNumEntriesPerRow(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview);


  // MMM Kernel wrappers struct
  // Because C++ doesn't support partial template specialization of functions.
  template<class Scalar,
	   class LocalOrdinal,
	   class GlobalOrdinal,
	   class Node,
           class LocalOrdinalViewType>
  struct KernelWrappers {
    static inline void mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
							 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
							 const LocalOrdinalViewType & Acol2Brow,
							 const LocalOrdinalViewType & Acol2Irow,
							 const LocalOrdinalViewType & Bcol2Ccol,
							 const LocalOrdinalViewType & Icol2Ccol,
							 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
							 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                         const std::string& label = std::string(),
							 const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                     const LocalOrdinalViewType & Acol2Brow,
                                                     const LocalOrdinalViewType & Acol2Irow,
                                                     const LocalOrdinalViewType & Bcol2Ccol,
                                                     const LocalOrdinalViewType & Icol2Ccol,
                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                     const std::string& label = std::string(),
                                                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
  };

  // Jacobi Kernel wrappers struct
  // Because C++ doesn't support partial template specialization of functions.

  template<class Scalar,
	   class LocalOrdinal,
	   class GlobalOrdinal,
	   class Node,
           class LocalOrdinalViewType>
  struct KernelWrappers2 {
    static inline void jacobi_A_B_newmatrix_kernel_wrapper(Scalar omega,
                                                           const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Dinv,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                           const LocalOrdinalViewType & Acol2Brow,
                                                           const LocalOrdinalViewType & Acol2Irow,
                                                           const LocalOrdinalViewType & Bcol2Ccol,
                                                           const LocalOrdinalViewType & Icol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                           const std::string& label = std::string(),
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void jacobi_A_B_reuse_kernel_wrapper(Scalar omega,
                                                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Dinv,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                       const LocalOrdinalViewType & Acol2Brow,
                                                       const LocalOrdinalViewType & Acol2Irow,
                                                           const LocalOrdinalViewType & Bcol2Ccol,
                                                       const LocalOrdinalViewType & Icol2Ccol,
                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                       Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                       const std::string& label = std::string(),
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
  };


  // Triple-Product Kernel wrappers struct
  // Because C++ doesn't support partial template specialization of functions.
  template<class Scalar,
	   class LocalOrdinal,
	   class GlobalOrdinal,
	   class Node,
           class LocalOrdinalViewType>
  struct KernelWrappers3 {
    static inline void mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                           const LocalOrdinalViewType & Acol2Prow,
                                                           const LocalOrdinalViewType & Acol2PIrow,
                                                           const LocalOrdinalViewType & Pcol2Ccol,
                                                           const LocalOrdinalViewType & PIcol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                           const std::string& label = std::string(),
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void mult_R_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                           const LocalOrdinalViewType & Acol2Prow,
                                                           const LocalOrdinalViewType & Acol2PIrow,
                                                           const LocalOrdinalViewType & Pcol2Ccol,
                                                           const LocalOrdinalViewType & PIcol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                           const std::string& label = std::string(),
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    static inline void mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                           const LocalOrdinalViewType & Acol2Prow,
                                                           const LocalOrdinalViewType & Acol2PIrow,
                                                           const LocalOrdinalViewType & Pcol2Ccol,
                                                           const LocalOrdinalViewType & PIcol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                           const std::string& label = std::string(),
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void mult_PT_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                           const LocalOrdinalViewType & Acol2Prow,
                                                           const LocalOrdinalViewType & Acol2PIrow,
                                                           const LocalOrdinalViewType & Pcol2Ccol,
                                                           const LocalOrdinalViewType & PIcol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                           const std::string& label = std::string(),
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
  };

  // This only merges matrices that look like B & Bimport, aka, they have no overlapping rows
  template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalOrdinalViewType>
  inline const typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type 
  merge_matrices(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                 const LocalOrdinalViewType & Acol2Brow,
                 const LocalOrdinalViewType & Acol2Irow,
                 const LocalOrdinalViewType & Bcol2Ccol,
                 const LocalOrdinalViewType & Icol2Ccol,  
                 const size_t mergedNodeNumCols);




template<class CrsMatrixType>
size_t C_estimate_nnz(CrsMatrixType & A, CrsMatrixType &B);

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
struct AddKernels
{
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
  typedef typename Node::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef typename device_type::memory_space memory_space;
  typedef typename crs_matrix_type::impl_scalar_type impl_scalar_type;
  typedef typename crs_matrix_type::local_matrix_type KCRS;
  typedef typename KCRS::values_type::non_const_type values_array;
  typedef typename KCRS::row_map_type::non_const_type row_ptrs_array;
  typedef typename KCRS::row_map_type row_ptrs_array_const;
  typedef typename KCRS::index_type::non_const_type col_inds_array;
  typedef typename map_type::local_map_type local_map_type;
  typedef typename Kokkos::View<GlobalOrdinal*, device_type> global_col_inds_array;
  typedef Kokkos::RangePolicy<execution_space> range_type;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_t, LocalOrdinal, impl_scalar_type,
              execution_space, memory_space, memory_space> KKH;

  /// \brief Given two matrices in CRS format, return their sum
  /// \pre A and B must both have column indices sorted within each row
  /// \param Avals Values array for A
  /// \param Arowptrs Row pointers array for A
  /// \param Acolinds Column indices array for A
  /// \param scalarA Scaling factor for A
  /// \param Bvals Values array for B
  /// \param Browptrs Row pointers array for B
  /// \param Bcolinds Column indices array for B
  /// \param scalarB Scaling factor for B
  /// \param[Out] Cvals Values array for C (allocated inside function)
  /// \param[Out] Crowptrs Row pointers array for C (allocated inside function)
  /// \param[Out] Ccolinds Column indices array for C (allocated inside function)
  static void addSorted(
    const values_array& Avals,
    const row_ptrs_array_const& Arowptrs,
    const col_inds_array& Acolinds,
    const impl_scalar_type scalarA,
    const values_array& Bvals,
    const row_ptrs_array_const& Browptrs,
    const col_inds_array& Bcolinds,
    const impl_scalar_type scalarB,
    values_array& Cvals,
    row_ptrs_array& Crowptrs,
    col_inds_array& Ccolinds);

  /// \brief Given two matrices in CRS format, return their sum
  /// \pre A and B don't need to be sorted, and column indices are given as global indices
  /// \param A A (local) matrix
  /// \param scalarA Scaling factor for A
  /// \param B B (local) matrix
  /// \param scalarB Scaling factor for B
  /// \param minGlobalCol The minimum global index owned by this processor in the column map
  /// \param globalNumCols The global size of the column map
  /// \param[Out] Cvals Values array for C (allocated inside function)
  /// \param[Out] Crowptrs Row pointers array for C (allocated inside function)
  /// \param[Out] Ccolinds Column indices array for C (allocated inside function)
  static void convertToGlobalAndAdd(
    const KCRS& A,
    const impl_scalar_type scalarA,
    const KCRS& B,
    const impl_scalar_type scalarB,
    const local_map_type& AcolMap,
    const local_map_type& BcolMap,
    values_array& Cvals,
    row_ptrs_array& Crowptrs,
    global_col_inds_array& Ccolinds);

  /// \brief Given two matrices in CRS format, return their sum
  /// \pre A and B don't need to be sorted
  /// \param Avals Values array for A
  /// \param Arowptrs Row pointers array for A
  /// \param Acolinds Column indices array for A
  /// \param scalarA Scaling factor for A
  /// \param Bvals Values array for B
  /// \param Browptrs Row pointers array for B
  /// \param Bcolinds Column indices array for B
  /// \param scalarB Scaling factor for B
  /// \param globalNumCols The global size of the column map
  /// \param[Out] Cvals Values array for C (allocated inside function)
  /// \param[Out] Crowptrs Row pointers array for C (allocated inside function)
  /// \param[Out] Ccolinds Column indices array for C (allocated inside function)
  static void addUnsorted(
    const values_array& Avals,
    const row_ptrs_array_const& Arowptrs,
    const col_inds_array& Acolinds,
    const impl_scalar_type scalarA,
    const values_array& Bvals,
    const row_ptrs_array_const& Browptrs,
    const col_inds_array& Bcolinds,
    const impl_scalar_type scalarB,
    GlobalOrdinal numGlobalCols,
    values_array& Cvals,
    row_ptrs_array& Crowptrs,
    col_inds_array& Ccolinds);
};

}//end namespace MMdetails

} // end of Tpetra namespace

#endif // TPETRA_MATRIXMATRIX_DECL_HPP


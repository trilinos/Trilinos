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

#ifndef TPETRA_TRIPLEMATRIXMULTIPLY_DECL_HPP
#define TPETRA_TRIPLEMATRIXMULTIPLY_DECL_HPP

#include <string>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "TpetraExt_MMHelpers.hpp"


/*! \file TpetraExt_TripleMatrixMultiply_decl.hpp

  The declarations for the class Tpetra::TripleMatrixMultiply and related non-member constructors.
*/

namespace Tpetra {

  namespace TripleMatrixMultiply {

    /// \brief Sparse matrix-matrix multiply
    ///
    /// Given CrsMatrix instances R, A and P, compute the product Ac = R*A*P,
    /// overwriting an existing CrsMatrix instance Ac with the result.
    ///
    /// \pre All four matrices R, A, P, and Ac must have uniquely owned row
    ///   Maps.
    /// \pre On input, Ac must have the same row Map as R.
    /// \pre R, A and P must be fill complete.
    /// \pre If Ac has a range Map on input, then Ac and R must have the
    ///   same range Maps.
    /// \pre If Ac has a domain Map on input, then P and Ac must have the
    ///   same domain Maps.
    ///
    /// For the latter two preconditions, recall that a matrix does not
    /// have a domain or range Map unless fillComplete has been called on
    /// it at least once.
    ///
    /// \param R [in] fill-complete sparse matrix.
    /// \param transposeR [in] Whether to use transpose of matrix R.
    /// \param A [in] fill-complete sparse matrix.
    /// \param transposeR [in] Whether to use transpose of matrix A.
    /// \param P [in] fill-complete sparse matrix.
    /// \param transposeB [in] Whether to use transpose of matrix P.
    /// \param Ac [in/out] On entry to this method, if Ac is fill complete,
    ///   then Ac's graph must have the correct structure, that is, its
    ///   structure must equal the structure of R*A*P. (This is currently
    ///   not implemented.)  On exit, C will be fill complete, unless the
    ///   last argument to this function is false.
    /// \param call_FillComplete_on_result [in] Optional argument;
    ///   defaults to true.  If false, C will <i>not</i> be fill complete
    ///   on output.
    template <class Scalar,
              class LocalOrdinal,
              class GlobalOrdinal,
              class Node>
    void MultiplyRAP(
                     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R,
                     bool transposeR,
                     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                     bool transposeA,
                     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& P,
                     bool transposeP,
                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                     bool call_FillComplete_on_result = true,
                     const std::string& label = std::string(),
                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


  } // namespace TripleMatrixMultiply

  namespace MMdetails{

    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_R_A_P_newmatrix(
                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                              const std::string& label = std::string(),
                              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_PT_A_P_newmatrix(
                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                               const std::string& label = std::string(),
                               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_R_A_P_newmatrix_kernel_wrapper(
                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                             CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                             const std::string& label = std::string(),
                                             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_PT_A_P_newmatrix_kernel_wrapper(
                                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                              const std::string& label = std::string(),
                                              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    // This struct contains the kernels that have a separate specialization for Node = OpenMP
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    struct KernelWrappers3MMM_Specialized {
      static inline void mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                             const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                             const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                             const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                             const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                             CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                             Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                             const std::string& label = std::string(),
                                                             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

      static inline void mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                              const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                              const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                              const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                              const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                              Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                              const std::string& label = std::string(),
                                                              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    };

    // Kernel wrappers struct (for non-specialized kernels)
    // Because C++ doesn't support partial template specialization of functions.
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    struct KernelWrappers3MMM {

      static inline void mult_PT_A_P_newmatrix_kernel_wrapper_2pass(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                    const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                                    const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                                    const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                                    const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                                    CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                    const std::string& label = std::string(),
                                                                    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
    };

  }//end namespace MMdetails

} // end of Tpetra namespace

#endif // TPETRA_TRIPLEMATRIXMULTIPLY_DECL_HPP

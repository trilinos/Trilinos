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

#ifndef TPETRA_MATRIXMATRIX_EXTRAKERNELS_DECL_HPP
#define TPETRA_MATRIXMATRIX_EXTRAKERNELS_DECL_HPP
#include "TpetraExt_MatrixMatrix_decl.hpp"


namespace Tpetra {

namespace MatrixMatrix {

  // This guy allows us to easily get an Unmanaged Kokkos View from a ManagedOne
  template <typename View>
  using UnmanagedView = Kokkos::View< typename View::data_type
                                      , typename View::array_layout
                                      , typename View::device_type
                                      , typename Kokkos::MemoryTraits< Kokkos::Unmanaged>
                                      >;

  namespace ExtraKernels {

    template<class CrsMatrixType>
    size_t C_estimate_nnz_per_row(CrsMatrixType & A, CrsMatrixType &B);

    template<class InRowptrArrayType, class InColindArrayType, class InValsArrayType,
             class OutRowptrType, class OutColindType, class OutValsType>
    void copy_out_from_thread_memory(const InRowptrArrayType & Inrowptr, const InColindArrayType &Incolind, const InValsArrayType & Invals,
                                       size_t m, double thread_chunk,
                                       OutRowptrType & Outrowptr, OutColindType &Outcolind, OutValsType & Outvals);

    /***************************** OpenMP Only Kernels *****************************/
#ifdef HAVE_TPETRA_INST_OPENMP
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void mult_A_B_newmatrix_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                   const LocalOrdinalViewType & Acol2Brow,
                                                                   const LocalOrdinalViewType & Acol2Irow,
                                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                                   const std::string& label,
                                                                   const Teuchos::RCP<Teuchos::ParameterList>& params);

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void mult_A_B_reuse_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                   const LocalOrdinalViewType & Acol2Brow,
                                                                   const LocalOrdinalViewType & Acol2Irow,
                                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                                   const std::string& label,
                                                                   const Teuchos::RCP<Teuchos::ParameterList>& params);

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void jacobi_A_B_newmatrix_LowThreadGustavsonKernel(Scalar omega,
                                                                     const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                     const LocalOrdinalViewType & Acol2Brow,
                                                                     const LocalOrdinalViewType & Acol2Irow,
                                                                     const LocalOrdinalViewType & Bcol2Ccol,
                                                                     const LocalOrdinalViewType & Icol2Ccol,
                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                                     const std::string& label,
                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params);

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void jacobi_A_B_reuse_LowThreadGustavsonKernel(Scalar omega,
                                                                     const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                     const LocalOrdinalViewType & Acol2Brow,
                                                                     const LocalOrdinalViewType & Acol2Irow,
                                                                     const LocalOrdinalViewType & Bcol2Ccol,
                                                                     const LocalOrdinalViewType & Icol2Ccol,
                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                                     const std::string& label,
                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params);
#endif

    /***************************** Generic Kernels *****************************/
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalOrdinalViewType>
    static inline void jacobi_A_B_newmatrix_MultiplyScaleAddKernel(Scalar omega,
                                                                   const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node> & Dinv,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                                   const LocalOrdinalViewType & Acol2Brow,
                                                                   const LocalOrdinalViewType & Acol2Irow,
                                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                                   const std::string& label,
                                                                   const Teuchos::RCP<Teuchos::ParameterList>& params);

  }// ExtraKernels
}//MatrixMatrix
}//Tpetra



#endif

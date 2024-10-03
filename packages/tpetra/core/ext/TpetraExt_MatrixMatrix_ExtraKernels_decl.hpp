// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

    // 2019 Apr 10 JJE:
    // copies data from thread local chunks into a unified CSR structure
    // 'const' on the inCol and inVals array is a lie.  The routine will deallocate
    // the thread local storage. Maybe they shouldn't be const. Or mark, non-const
    // and have a helper function for the actual copies that takes these as const
    // . The point of const is that we want the loops to optimize assuming the
    // RHS is unchanging
    template<class InColindArrayType,
             class InValsArrayType,
             class OutRowptrType,
             class OutColindType,
             class OutValsType>
    void copy_out_from_thread_memory(const OutColindType& thread_total_nnz,
                                     const InColindArrayType& Incolind,
                                     const InValsArrayType& Invals,
                                     const size_t m,
                                     const double thread_chunk,
                                     OutRowptrType& Outrowptr,
                                     OutColindType& Outcolind,
                                     OutValsType& Outvals);

    /***************************** Matrix-Matrix OpenMP Only Kernels *****************************/
#ifdef HAVE_TPETRA_INST_OPENMP
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void mult_A_B_newmatrix_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                                   const LocalOrdinalViewType & Acol2Brow,
                                                                   const LocalOrdinalViewType & Acol2Irow,
                                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                                   const std::string& label,
                                                                   const Teuchos::RCP<Teuchos::ParameterList>& params);

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void mult_A_B_reuse_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                                   const LocalOrdinalViewType & Acol2Brow,
                                                                   const LocalOrdinalViewType & Acol2Irow,
                                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                                   const std::string& label,
                                                                   const Teuchos::RCP<Teuchos::ParameterList>& params);

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void jacobi_A_B_newmatrix_LowThreadGustavsonKernel(Scalar omega,
                                                                     const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> & Dinv,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                                     const LocalOrdinalViewType & Acol2Brow,
                                                                     const LocalOrdinalViewType & Acol2Irow,
                                                                     const LocalOrdinalViewType & Bcol2Ccol,
                                                                     const LocalOrdinalViewType & Icol2Ccol,
                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                                     const std::string& label,
                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params);

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void jacobi_A_B_reuse_LowThreadGustavsonKernel(Scalar omega,
                                                                     const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> & Dinv,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                                     const LocalOrdinalViewType & Acol2Brow,
                                                                     const LocalOrdinalViewType & Acol2Irow,
                                                                     const LocalOrdinalViewType & Bcol2Ccol,
                                                                     const LocalOrdinalViewType & Icol2Ccol,
                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                                     const std::string& label,
                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params);
#endif

    /***************************** Matrix-Matrix Generic Kernels *****************************/
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalOrdinalViewType>
    static inline void jacobi_A_B_newmatrix_MultiplyScaleAddKernel(Scalar omega,
                                                                   const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node> & Dinv,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                                   const LocalOrdinalViewType & Acol2rrow,
                                                                   const LocalOrdinalViewType & Acol2Irow,
                                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                                   const std::string& label,
                                                                   const Teuchos::RCP<Teuchos::ParameterList>& params);


    /***************************** Triple Product OpenMP Only Kernels *****************************/
#ifdef HAVE_TPETRA_INST_OPENMP
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
    static inline void mult_R_A_P_newmatrix_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                                     const LocalOrdinalViewType & Acol2Prow,
                                                                     const LocalOrdinalViewType & Acol2PIrow,
                                                                     const LocalOrdinalViewType & Pcol2Accol,
                                                                     const LocalOrdinalViewType & PIcol2Accol,
                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                                     const std::string& label = std::string(),
                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
#endif


  }// ExtraKernels
}//MatrixMatrix
}//Tpetra



#endif

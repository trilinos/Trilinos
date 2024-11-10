// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_Array.hpp>
#include "Xpetra_EpetraCrsMatrix.hpp"

namespace Xpetra {

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraCrsMatrixT<int, Xpetra::EpetraNode>;
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraCrsMatrixT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraCrsMatrixT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraCrsMatrixT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type>;
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type>;
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraCrsMatrixT<int, default_node_type>;
#endif  // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraCrsMatrixT<long long, Xpetra::EpetraNode>;
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraCrsMatrixT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraCrsMatrixT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraCrsMatrixT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type>;
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type>;
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraCrsMatrixT<long long, default_node_type>;
#endif  // HAVE_XPETRA_TPETRA
#endif

}  // namespace Xpetra

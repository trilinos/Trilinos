// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_EpetraIntMultiVector.hpp"
#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"

namespace Xpetra {

// TODO: move that elsewhere
template <class GlobalOrdinal, class Node>
Epetra_IntMultiVector &toEpetra(MultiVector<int, int, GlobalOrdinal, Node> &x) {
  XPETRA_DYNAMIC_CAST(EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_IntMultiVector();
}

template <class GlobalOrdinal, class Node>
const Epetra_IntMultiVector &toEpetra(const MultiVector<int, int, GlobalOrdinal, Node> &x) {
  XPETRA_DYNAMIC_CAST(const EpetraIntMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_IntMultiVector();
}
//

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraIntMultiVectorT<int, Xpetra::EpetraNode>;
template Epetra_IntMultiVector &toEpetra<int, Xpetra::EpetraNode>(MultiVector<int, int, int, Xpetra::EpetraNode> &);
template const Epetra_IntMultiVector &toEpetra<int, Xpetra::EpetraNode>(const MultiVector<int, int, int, Xpetra::EpetraNode> &);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraIntMultiVectorT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template Epetra_IntMultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(MultiVector<int, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_IntMultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const MultiVector<int, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraIntMultiVectorT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template Epetra_IntMultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(MultiVector<int, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_IntMultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const MultiVector<int, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraIntMultiVectorT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template Epetra_IntMultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(MultiVector<int, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_IntMultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const MultiVector<int, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraIntMultiVectorT<int, default_node_type>;
template Epetra_IntMultiVector &toEpetra<int, default_node_type>(MultiVector<int, int, int, default_node_type> &);
template const Epetra_IntMultiVector &toEpetra<int, default_node_type>(const MultiVector<int, int, int, default_node_type> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraIntMultiVectorT<int, default_node_type>;
template Epetra_IntMultiVector &toEpetra<int, default_node_type>(MultiVector<int, int, int, default_node_type> &);
template const Epetra_IntMultiVector &toEpetra<int, default_node_type>(const MultiVector<int, int, int, default_node_type> &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraIntMultiVectorT<int, default_node_type>;
template Epetra_IntMultiVector &toEpetra<int, default_node_type>(MultiVector<int, int, int, default_node_type> &);
template const Epetra_IntMultiVector &toEpetra<int, default_node_type>(const MultiVector<int, int, int, default_node_type> &);
#endif  // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraIntMultiVectorT<long long, Xpetra::EpetraNode>;
template Epetra_IntMultiVector &toEpetra<long long, Xpetra::EpetraNode>(MultiVector<int, int, long long, Xpetra::EpetraNode> &);
template const Epetra_IntMultiVector &toEpetra<long long, Xpetra::EpetraNode>(const MultiVector<int, int, long long, Xpetra::EpetraNode> &);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraIntMultiVectorT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template Epetra_IntMultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(MultiVector<int, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_IntMultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const MultiVector<int, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraIntMultiVectorT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template Epetra_IntMultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(MultiVector<int, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_IntMultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const MultiVector<int, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraIntMultiVectorT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template Epetra_IntMultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(MultiVector<int, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_IntMultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const MultiVector<int, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraIntMultiVectorT<long long, default_node_type>;
template Epetra_IntMultiVector &toEpetra<long long, default_node_type>(MultiVector<int, int, long long, default_node_type> &);
template const Epetra_IntMultiVector &toEpetra<long long, default_node_type>(const MultiVector<int, int, long long, default_node_type> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraIntMultiVectorT<long long, default_node_type>;
template Epetra_IntMultiVector &toEpetra<long long, default_node_type>(MultiVector<int, int, long long, default_node_type> &);
template const Epetra_IntMultiVector &toEpetra<long long, default_node_type>(const MultiVector<int, int, long long, default_node_type> &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraIntMultiVectorT<long long, default_node_type>;
template Epetra_IntMultiVector &toEpetra<long long, default_node_type>(MultiVector<int, int, long long, default_node_type> &);
template const Epetra_IntMultiVector &toEpetra<long long, default_node_type>(const MultiVector<int, int, long long, default_node_type> &);
#endif  // HAVE_XPETRA_TPETRA
#endif

}  // namespace Xpetra

// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_EpetraVector.hpp"

// TODO: replace double -> Scalar etc.

namespace Xpetra {

// TODO: move that elsewhere
template <class GlobalOrdinal, class Node>
Epetra_Vector &toEpetra(Vector<double, int, GlobalOrdinal, Node> &x) {
  XPETRA_DYNAMIC_CAST(EpetraVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_Vector();
}

template <class GlobalOrdinal, class Node>
const Epetra_Vector &toEpetra(const Vector<double, int, GlobalOrdinal, Node> &x) {
  XPETRA_DYNAMIC_CAST(const EpetraVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_Vector();
}
//

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraVectorT<int, Xpetra::EpetraNode>;
template Epetra_Vector &toEpetra<int, Xpetra::EpetraNode>(Vector<double, int, int, Xpetra::EpetraNode> &);
template const Epetra_Vector &toEpetra<int, Xpetra::EpetraNode>(const Vector<double, int, int, Xpetra::EpetraNode> &);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraVectorT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template Epetra_Vector &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(Vector<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_Vector &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Vector<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraVectorT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template Epetra_Vector &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(Vector<double, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_Vector &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Vector<double, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraVectorT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template Epetra_Vector &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(Vector<double, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_Vector &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Vector<double, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraVectorT<int, default_node_type>;
template Epetra_Vector &toEpetra<int, default_node_type>(Vector<double, int, int, default_node_type> &);
template const Epetra_Vector &toEpetra<int, default_node_type>(const Vector<double, int, int, default_node_type> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraVectorT<int, default_node_type>;
template Epetra_Vector &toEpetra<int, default_node_type>(Vector<double, int, int, default_node_type> &);
template const Epetra_Vector &toEpetra<int, default_node_type>(const Vector<double, int, int, default_node_type> &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraVectorT<int, default_node_type>;
template Epetra_Vector &toEpetra<int, default_node_type>(Vector<double, int, int, default_node_type> &);
template const Epetra_Vector &toEpetra<int, default_node_type>(const Vector<double, int, int, default_node_type> &);
#endif  // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraVectorT<long long, Xpetra::EpetraNode>;
template Epetra_Vector &toEpetra<long long, Xpetra::EpetraNode>(Vector<double, int, long long, Xpetra::EpetraNode> &);
template const Epetra_Vector &toEpetra<long long, Xpetra::EpetraNode>(const Vector<double, int, long long, Xpetra::EpetraNode> &);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraVectorT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template Epetra_Vector &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(Vector<double, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_Vector &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Vector<double, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraVectorT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template Epetra_Vector &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(Vector<double, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_Vector &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Vector<double, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraVectorT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template Epetra_Vector &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(Vector<double, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_Vector &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Vector<double, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraVectorT<long long, default_node_type>;
template Epetra_Vector &toEpetra<long long, default_node_type>(Vector<double, int, long long, default_node_type> &);
template const Epetra_Vector &toEpetra<long long, default_node_type>(const Vector<double, int, long long, default_node_type> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraVectorT<long long, default_node_type>;
template Epetra_Vector &toEpetra<long long, default_node_type>(Vector<double, int, long long, default_node_type> &);
template const Epetra_Vector &toEpetra<long long, default_node_type>(const Vector<double, int, long long, default_node_type> &);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraVectorT<long long, default_node_type>;
template Epetra_Vector &toEpetra<long long, default_node_type>(Vector<double, int, long long, default_node_type> &);
template const Epetra_Vector &toEpetra<long long, default_node_type>(const Vector<double, int, long long, default_node_type> &);
#endif  // HAVE_XPETRA_TPETRA
#endif

}  // namespace Xpetra

// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_EPETRA

#include "Xpetra_EpetraMap.hpp"

namespace Xpetra {

template <class GlobalOrdinal, class Node>
const Epetra_Map &toEpetra(const Map<int, GlobalOrdinal, Node> &map) {
  const EpetraMapT<GlobalOrdinal, Node> &epetraMap = dynamic_cast<const EpetraMapT<GlobalOrdinal, Node> &>(*map.getMap());
  return epetraMap.getEpetra_Map();
}

template <class GlobalOrdinal, class Node>
const Epetra_Map &toEpetra(const RCP<const Map<int, GlobalOrdinal, Node> > &map) {
  XPETRA_RCP_DYNAMIC_CAST(const EpetraMapT<GlobalOrdinal XPETRA_COMMA Node>, map->getMap(), epetraMap, "toEpetra");
  return epetraMap->getEpetra_Map();
}

template <class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<int, GlobalOrdinal, Node> > toXpetra(const Epetra_BlockMap &map) {
  RCP<const Epetra_BlockMap> m = rcp(new Epetra_BlockMap(map));
  return rcp(new EpetraMapT<GlobalOrdinal, Node>(m));
}

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template const RCP<const Map<int, int, Xpetra::EpetraNode> > toXpetra<int, Xpetra::EpetraNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, Xpetra::EpetraNode>(const RCP<const Map<int, int, Xpetra::EpetraNode> > &map);
template const Epetra_Map &toEpetra<int, Xpetra::EpetraNode>(const Map<int, int, Xpetra::EpetraNode> &map);
#endif
#ifdef HAVE_TPETRA_INST_SERIAL
// template class EpetraMapT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode >;
template const RCP<const Map<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const RCP<const Map<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > &map);
template const Epetra_Map &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Map<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &map);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
// template class EpetraMapT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template const RCP<const Map<int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const RCP<const Map<int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > &map);
template const Epetra_Map &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Map<int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &map);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
// template class EpetraMapT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >;
template const RCP<const Map<int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const RCP<const Map<int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > &map);
template const Epetra_Map &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Map<int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &map);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
// template class EpetraMapT<int, default_node_type >;
template const RCP<const Map<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, default_node_type>(const RCP<const Map<int, int, default_node_type> > &map);
template const Epetra_Map &toEpetra<int, default_node_type>(const Map<int, int, default_node_type> &map);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
// template class EpetraMapT<int, default_node_type >;
template const RCP<const Map<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, default_node_type>(const RCP<const Map<int, int, default_node_type> > &map);
template const Epetra_Map &toEpetra<int, default_node_type>(const Map<int, int, default_node_type> &map);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template const RCP<const Map<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<int, default_node_type>(const RCP<const Map<int, int, default_node_type> > &map);
template const Epetra_Map &toEpetra<int, default_node_type>(const Map<int, int, default_node_type> &map);
#endif  // HAVE_XPETRA_TPETRA

#endif  // XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template const RCP<const Map<int, long long, Xpetra::EpetraNode> > toXpetra<long long, Xpetra::EpetraNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, Xpetra::EpetraNode>(const RCP<const Map<int, long long, Xpetra::EpetraNode> > &map);
template const Epetra_Map &toEpetra<long long, Xpetra::EpetraNode>(const Map<int, long long, Xpetra::EpetraNode> &map);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
// template class EpetraMapT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode >;
template const RCP<const Map<int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const RCP<const Map<int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > &map);
template const Epetra_Map &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Map<int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &map);
#endif

#ifdef HAVE_TPETRA_INST_PTHREAD
// template class EpetraMapT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template const RCP<const Map<int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const RCP<const Map<int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > &map);
template const Epetra_Map &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Map<int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &map);
#endif

#ifdef HAVE_TPETRA_INST_OPENMP
// template class EpetraMapT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode >;
template const RCP<const Map<int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const RCP<const Map<int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > &map);
template const Epetra_Map &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Map<int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &map);
#endif

#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
// template class EpetraMapT<long long, default_node_type >;
template const RCP<const Map<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, default_node_type>(const RCP<const Map<int, long long, default_node_type> > &map);
template const Epetra_Map &toEpetra<long long, default_node_type>(const Map<int, long long, default_node_type> &map);
#endif

#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
// template class EpetraMapT<long long, default_node_type >;
template const RCP<const Map<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, default_node_type>(const RCP<const Map<int, long long, default_node_type> > &map);
template const Epetra_Map &toEpetra<long long, default_node_type>(const Map<int, long long, default_node_type> &map);
#endif

#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template const RCP<const Map<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_BlockMap &map);
template const Epetra_Map &toEpetra<long long, default_node_type>(const RCP<const Map<int, long long, default_node_type> > &map);
template const Epetra_Map &toEpetra<long long, default_node_type>(const Map<int, long long, default_node_type> &map);
#endif  // HAVE_XPETRA_TPETRA

#endif  // HAVE_XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

}  // namespace Xpetra

#endif  // HAVE_XPETRA_EPETRA

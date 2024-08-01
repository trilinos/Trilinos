// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class GlobalOrdinal, class Node>
RCP<const Export<int, GlobalOrdinal, Node> > toXpetra(const Epetra_Export *exp) {
  if (exp != NULL) {
    RCP<const Epetra_Export> eexp = rcp(new Epetra_Export(*exp));  // NOTE: non consitent: return pointer, take ref
    return rcp(new Xpetra::EpetraExportT<GlobalOrdinal, Node>(eexp));
  }

  return Teuchos::null;
}

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraExportT<int, Xpetra::EpetraNode>;
template RCP<const Export<int, int, Xpetra::EpetraNode> > toXpetra<int, Xpetra::EpetraNode>(const Epetra_Export *);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraExportT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template RCP<const Export<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraExportT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template RCP<const Export<int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraExportT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template RCP<const Export<int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraExportT<int, default_node_type>;
template RCP<const Export<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraExportT<int, default_node_type>;
template RCP<const Export<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_Export *);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraExportT<int, default_node_type>;
template RCP<const Export<int, int, default_node_type> > toXpetra<int, default_node_type>(const Epetra_Export *);
#endif  // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraExportT<long long, Xpetra::EpetraNode>;
template RCP<const Export<int, long long, Xpetra::EpetraNode> > toXpetra<long long, Xpetra::EpetraNode>(const Epetra_Export *);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraExportT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template RCP<const Export<int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraExportT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template RCP<const Export<int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraExportT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template RCP<const Export<int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode default_node_type;
template class EpetraExportT<long long, default_node_type>;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_Export *);
#endif
#ifdef HAVE_TPETRA_INST_HIP
typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode default_node_type;
template class EpetraExportT<long long, default_node_type>;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_Export *);
#endif
#else
// Tpetra is disabled and Kokkos not available: use dummy node type
typedef Xpetra::EpetraNode default_node_type;
template class EpetraExportT<long long, default_node_type>;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long, default_node_type>(const Epetra_Export *);
#endif  // HAVE_XPETRA_TPETRA
#endif

}  // namespace Xpetra

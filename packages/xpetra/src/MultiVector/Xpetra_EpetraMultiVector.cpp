// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_EpetraMultiVector.hpp"

#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_EpetraVector.hpp"

#include "Epetra_SerialComm.h"

namespace Xpetra {

// specialization for GO=int and NO=EpetraNode
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
Teuchos::RCP<const Vector<double, int, int, EpetraNode> > EpetraMultiVectorT<int, EpetraNode>::getVector(size_t j) const {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new Xpetra::EpetraVectorT<int, EpetraNode>(vec_, j));  // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}

//! Return a Vector which is a nonconst view of column j.
Teuchos::RCP<Vector<double, int, int, EpetraNode> > EpetraMultiVectorT<int, EpetraNode>::getVectorNonConst(size_t j) {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new EpetraVectorT<int, EpetraNode>(vec_, j));  // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}
#endif

// specialization for GO=long long and NO=EpetraNode
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
Teuchos::RCP<const Vector<double, int, long long, EpetraNode> > EpetraMultiVectorT<long long, EpetraNode>::getVector(size_t j) const {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new Xpetra::EpetraVectorT<long long, EpetraNode>(vec_, j));  // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}

//! Return a Vector which is a nonconst view of column j.
Teuchos::RCP<Vector<double, int, long long, EpetraNode> > EpetraMultiVectorT<long long, EpetraNode>::getVectorNonConst(size_t j) {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new EpetraVectorT<long long, EpetraNode>(vec_, j));  // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}
#endif

// TODO: move that elsewhere
template <class GlobalOrdinal, class Node>
const Epetra_MultiVector &toEpetra(const MultiVector<double, int, GlobalOrdinal, Node> &x) {
  XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_MultiVector();
}

template <class GlobalOrdinal, class Node>
Epetra_MultiVector &toEpetra(MultiVector<double, int, GlobalOrdinal, Node> &x) {
  XPETRA_DYNAMIC_CAST(EpetraMultiVectorT<GlobalOrdinal XPETRA_COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_MultiVector();
}
//

template <class GlobalOrdinal, class Node>
RCP<MultiVector<double, int, GlobalOrdinal, Node> > toXpetra(RCP<Epetra_MultiVector> vec) {
  if (!vec.is_null())
    return rcp(new EpetraMultiVectorT<GlobalOrdinal, Node>(vec));

  return Teuchos::null;
}

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraMultiVectorT<int, Xpetra::EpetraNode>;
template RCP<MultiVector<double, int, int, Xpetra::EpetraNode> > toXpetra<int, Xpetra::EpetraNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, Xpetra::EpetraNode>(MultiVector<double, int, int, Xpetra::EpetraNode> &);
template const Epetra_MultiVector &toEpetra<int, Xpetra::EpetraNode>(const MultiVector<double, int, int, Xpetra::EpetraNode> &);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraMultiVectorT<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template RCP<MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraMultiVectorT<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template RCP<MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraMultiVectorT<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template RCP<MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
template class EpetraMultiVectorT<int, Tpetra::KokkosCompat::KokkosCudaWrapperNode>;
template RCP<MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosCudaWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosCudaWrapperNode>(MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode> &);
template const Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosCudaWrapperNode>(const MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosCudaWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
template class EpetraMultiVectorT<int, Tpetra::KokkosCompat::KokkosHIPWrapperNode>;
template RCP<MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > toXpetra<int, Tpetra::KokkosCompat::KokkosHIPWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosHIPWrapperNode>(MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode> &);
template const Epetra_MultiVector &toEpetra<int, Tpetra::KokkosCompat::KokkosHIPWrapperNode>(const MultiVector<double, int, int, Tpetra::KokkosCompat::KokkosHIPWrapperNode> &);
#endif
#else  // Tpetra is disabled
typedef Xpetra::EpetraNode default_node_type;
template class EpetraMultiVectorT<int, default_node_type>;
template RCP<MultiVector<double, int, int, default_node_type> > toXpetra<int, default_node_type>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<int, default_node_type>(MultiVector<double, int, int, default_node_type> &);
template const Epetra_MultiVector &toEpetra<int, default_node_type>(const MultiVector<double, int, int, default_node_type> &);
#endif
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#include "TpetraCore_config.h"
#if ((defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_OPENMP)) || \
     (!defined(EPETRA_HAVE_OMP) && !defined(HAVE_TPETRA_INST_SERIAL)))
template class EpetraMultiVectorT<long long, Xpetra::EpetraNode>;
template RCP<MultiVector<double, int, long long, Xpetra::EpetraNode> > toXpetra<long long, Xpetra::EpetraNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, Xpetra::EpetraNode>(MultiVector<double, int, long long, Xpetra::EpetraNode> &);
template const Epetra_MultiVector &toEpetra<long long, Xpetra::EpetraNode>(const MultiVector<double, int, long long, Xpetra::EpetraNode> &);
#endif

#ifdef HAVE_TPETRA_INST_SERIAL
template class EpetraMultiVectorT<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>;
template RCP<MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
template const Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(const MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosSerialWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_PTHREAD
template class EpetraMultiVectorT<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>;
template RCP<MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
template const Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode>(const MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_OPENMP
template class EpetraMultiVectorT<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
template RCP<MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
template const Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>(const MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_CUDA
template class EpetraMultiVectorT<long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode>;
template RCP<MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode>(MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode> &);
template const Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode>(const MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosCudaWrapperNode> &);
#endif
#ifdef HAVE_TPETRA_INST_HIP
template class EpetraMultiVectorT<long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode>;
template RCP<MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > toXpetra<long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode>(MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode> &);
template const Epetra_MultiVector &toEpetra<long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode>(const MultiVector<double, int, long long, Tpetra::KokkosCompat::KokkosHIPWrapperNode> &);
#endif
#else  // Tpetra is disabled
typedef Xpetra::EpetraNode default_node_type;
template class EpetraMultiVectorT<long long, default_node_type>;
template RCP<MultiVector<double, int, long long, default_node_type> > toXpetra<long long, default_node_type>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector &toEpetra<long long, default_node_type>(MultiVector<double, int, long long, default_node_type> &);
template const Epetra_MultiVector &toEpetra<long long, default_node_type>(const MultiVector<double, int, long long, default_node_type> &);
#endif
#endif

}  // namespace Xpetra

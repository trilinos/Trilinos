// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
Teuchos::RCP< const Vector< double, int, int, EpetraNode > > EpetraMultiVectorT<int, EpetraNode>::getVector(size_t j) const {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new Xpetra::EpetraVectorT<int,EpetraNode>(vec_, j)); // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}

//! Return a Vector which is a nonconst view of column j.
Teuchos::RCP< Vector< double, int, int, EpetraNode > > EpetraMultiVectorT<int,EpetraNode>::getVectorNonConst(size_t j) {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new EpetraVectorT<int,EpetraNode>(vec_, j)); // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}
#endif

// specialization for GO=long long and NO=EpetraNode
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
Teuchos::RCP< const Vector< double, int, long long, EpetraNode > > EpetraMultiVectorT<long long, EpetraNode>::getVector(size_t j) const {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new Xpetra::EpetraVectorT<long long,EpetraNode>(vec_, j)); // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}

//! Return a Vector which is a nonconst view of column j.
Teuchos::RCP< Vector< double, int, long long, EpetraNode > > EpetraMultiVectorT<long long,EpetraNode>::getVectorNonConst(size_t j) {
  XPETRA_MONITOR("EpetraMultiVectorT::getVector");
  return rcp(new EpetraVectorT<long long,EpetraNode>(vec_, j)); // See constructor EpetraVectorT(const RCP<EpetraMultiVectorT> &mv, size_t j) for more info
}
#endif

// TODO: move that elsewhere
template<class GlobalOrdinal,class Node>
const Epetra_MultiVector & toEpetra(const MultiVector<double, int, GlobalOrdinal,Node> & x) {
  XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_MultiVector();
}

template<class GlobalOrdinal,class Node>
Epetra_MultiVector & toEpetra(MultiVector<double, int, GlobalOrdinal,Node> & x) {
  XPETRA_DYNAMIC_CAST(      EpetraMultiVectorT<GlobalOrdinal COMMA Node>, x, tX, "toEpetra");
  return *tX.getEpetra_MultiVector();
}
//

template<class GlobalOrdinal,class Node>
RCP<MultiVector<double, int, GlobalOrdinal,Node> > toXpetra(RCP<Epetra_MultiVector> vec) {
  if (!vec.is_null())
    return rcp(new EpetraMultiVectorT<GlobalOrdinal,Node>(vec));

  return Teuchos::null;
}



#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraMultiVectorT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP<MultiVector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosSerialWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<int,Kokkos::Compat::KokkosSerialWrapperNode >(MultiVector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> &);
template const Epetra_MultiVector & toEpetra<int, Kokkos::Compat::KokkosSerialWrapperNode >(const MultiVector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraMultiVectorT<int, Kokkos::Compat::KokkosThreadsWrapperNode >;
template RCP<MultiVector<double, int, int, Kokkos::Compat::KokkosThreadsWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosThreadsWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<int,Kokkos::Compat::KokkosThreadsWrapperNode >(MultiVector<double, int, int, Kokkos::Compat::KokkosThreadsWrapperNode> &);
template const Epetra_MultiVector & toEpetra<int, Kokkos::Compat::KokkosThreadsWrapperNode >(const MultiVector<double, int, int, Kokkos::Compat::KokkosThreadsWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraMultiVectorT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP<MultiVector<double, int, int, Kokkos::Compat::KokkosOpenMPWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosOpenMPWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<int,Kokkos::Compat::KokkosOpenMPWrapperNode >(MultiVector<double, int, int, Kokkos::Compat::KokkosOpenMPWrapperNode> &);
template const Epetra_MultiVector & toEpetra<int, Kokkos::Compat::KokkosOpenMPWrapperNode >(const MultiVector<double, int, int, Kokkos::Compat::KokkosOpenMPWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_CUDA
template class EpetraMultiVectorT<int, Kokkos::Compat::KokkosCudaWrapperNode >;
template RCP<MultiVector<double, int, int, Kokkos::Compat::KokkosCudaWrapperNode > > toXpetra<int, Kokkos::Compat::KokkosCudaWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<int,Kokkos::Compat::KokkosCudaWrapperNode >(MultiVector<double, int, int,Kokkos::Compat::KokkosCudaWrapperNode> &);
template const Epetra_MultiVector & toEpetra<int, Kokkos::Compat::KokkosCudaWrapperNode >(const MultiVector<double, int, int, Kokkos::Compat::KokkosCudaWrapperNode > &);
#endif
#else
typedef Kokkos::Compat::KokkosSerialWrapperNode default_node_type;
template class EpetraMultiVectorT<int, default_node_type >;
template RCP<MultiVector<double, int, int, default_node_type > > toXpetra<int, default_node_type>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<int,default_node_type >(MultiVector<double, int, int,default_node_type> &);
template const Epetra_MultiVector & toEpetra<int, default_node_type >(const MultiVector<double, int, int, default_node_type > &);
#endif // HAVE_XPETRA_SERIAL
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraMultiVectorT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP<MultiVector<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosSerialWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<long long,Kokkos::Compat::KokkosSerialWrapperNode >(MultiVector<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode> &);
template const Epetra_MultiVector & toEpetra<long long, Kokkos::Compat::KokkosSerialWrapperNode >(const MultiVector<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraMultiVectorT<long long,Kokkos::Compat::KokkosThreadsWrapperNode >;
template RCP<MultiVector<double, int, long long, Kokkos::Compat::KokkosThreadsWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosThreadsWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<long long,Kokkos::Compat::KokkosThreadsWrapperNode >(MultiVector<double, int, long long,Kokkos::Compat::KokkosThreadsWrapperNode> &);
template const Epetra_MultiVector & toEpetra<long long, Kokkos::Compat::KokkosThreadsWrapperNode >(const MultiVector<double, int, long long, Kokkos::Compat::KokkosThreadsWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraMultiVectorT<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP<MultiVector<double, int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosOpenMPWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<long long,Kokkos::Compat::KokkosOpenMPWrapperNode >(MultiVector<double, int, long long,Kokkos::Compat::KokkosOpenMPWrapperNode> &);
template const Epetra_MultiVector & toEpetra<long long, Kokkos::Compat::KokkosOpenMPWrapperNode >(const MultiVector<double, int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode > &);
#endif
#ifdef HAVE_XPETRA_CUDA
template class EpetraMultiVectorT<long long, Kokkos::Compat::KokkosCudaWrapperNode >;
template RCP<MultiVector<double, int, long long, Kokkos::Compat::KokkosCudaWrapperNode > > toXpetra<long long, Kokkos::Compat::KokkosCudaWrapperNode>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<long long,Kokkos::Compat::KokkosCudaWrapperNode >(MultiVector<double, int, long long, Kokkos::Compat::KokkosCudaWrapperNode> &);
template const Epetra_MultiVector & toEpetra<long long, Kokkos::Compat::KokkosCudaWrapperNode >(const MultiVector<double, int, long long, Kokkos::Compat::KokkosCudaWrapperNode > &);
#endif
#else
typedef Kokkos::Compat::KokkosSerialWrapperNode default_node_type;
template class EpetraMultiVectorT<long long, default_node_type >;
template RCP<MultiVector<double, int, long long, default_node_type > > toXpetra<long long, default_node_type>(RCP<Epetra_MultiVector>);
template Epetra_MultiVector & toEpetra<long long,default_node_type >(MultiVector<double, int, long long,default_node_type> &);
template const Epetra_MultiVector & toEpetra<long long, default_node_type >(const MultiVector<double, int, long long, default_node_type > &);
#endif // HAVE_XPETRA_SERIAL
#endif

} // namespace Xpetra

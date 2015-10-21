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
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template<class GlobalOrdinal, class Node>
  RCP<const Export<int, GlobalOrdinal, Node> > toXpetra(const Epetra_Export *exp) {
    if (exp != NULL) {
      RCP<const Epetra_Export> eexp = rcp(new Epetra_Export(*exp)); //NOTE: non consitent: return pointer, take ref
      return rcp(new Xpetra::EpetraExportT<GlobalOrdinal, Node>(eexp));
    }

    return Teuchos::null;
  }

  template<class EpetraGlobalOrdinal, class Node>
  EpetraExportT<EpetraGlobalOrdinal, Node>::EpetraExportT(const Teuchos::RCP<const map_type> & source, const Teuchos::RCP<const map_type> & target)
    : export_(rcp(new Epetra_Export(toEpetra<EpetraGlobalOrdinal,Node>(source), toEpetra<EpetraGlobalOrdinal,Node>(target)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  //
  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal, Node>::getExportPIDs() const { XPETRA_MONITOR("EpetraExportT::getExportImageIDs"); return ArrayView<const int> (export_->ExportPIDs(),export_->NumExportIDs()); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal, Node>::getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getPermuteFromLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getExportImageIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal, Node>::getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getPermuteToLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getPermuteToLIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  size_t EpetraExportT<EpetraGlobalOrdinal, Node>::getNumRemoteIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumRemoteIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getNumRemoteIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal, Node>::getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getRemoteLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getRemoteLIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  size_t EpetraExportT<EpetraGlobalOrdinal, Node>::getNumExportIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumExportIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getNumExportIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal, Node>::getExportLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getExportLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getExportLIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraExportT<EpetraGlobalOrdinal, Node>::print(std::ostream &os) const { XPETRA_MONITOR("EpetraExportT::");  export_->Print(os); }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraExportT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP<const Export<int,int,Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<int,Kokkos::Compat::KokkosSerialWrapperNode >(const Epetra_Export *);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraExportT<int, Kokkos::Compat::KokkosThreadsWrapperNode>;
template RCP<const Export<int,int,Kokkos::Compat::KokkosThreadsWrapperNode> > toXpetra<int,Kokkos::Compat::KokkosThreadsWrapperNode >(const Epetra_Export *);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraExportT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP<const Export<int, int, Kokkos::Compat::KokkosOpenMPWrapperNode> > toXpetra<int,Kokkos::Compat::KokkosOpenMPWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::View<int>::HostMirror::execution_space default_host_execution_space;
typedef Kokkos::Compat::KokkosDeviceWrapperNode<host_execution_space, Kokkos::HostSpace> default_node_type;
template class EpetraExportT<int, default_node_type >;
template RCP<const Export<int, int, default_node_type> > toXpetra<int,default_node_type>(const Epetra_Export *);
#endif
#else
  // TODO What, if Tpetra is disabled? Use fake Kokkos thing?
#endif // HAVE_XPETRA_TPETRA
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraExportT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP<const Export<int,long long,Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<long long,Kokkos::Compat::KokkosSerialWrapperNode >(const Epetra_Export *);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraExportT<long long, Kokkos::Compat::KokkosThreadsWrapperNode>;
template RCP<const Export<int,long long,Kokkos::Compat::KokkosThreadsWrapperNode> > toXpetra<long long,Kokkos::Compat::KokkosThreadsWrapperNode >(const Epetra_Export *);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraExportT<long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP<const Export<int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode> > toXpetra<long long,Kokkos::Compat::KokkosOpenMPWrapperNode>(const Epetra_Export *);
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::View<long long>::HostMirror::execution_space default_host_execution_space;
typedef Kokkos::Compat::KokkosDeviceWrapperNode<host_execution_space, Kokkos::HostSpace> default_node_type;
template class EpetraExportT<long long, default_node_type >;
template RCP<const Export<int, long long, default_node_type> > toXpetra<long long,default_node_type>(const Epetra_Export *);
#endif
#else
  // TODO: What, if Tpetra is disabled? Use fake Kokkos thing?
#endif // HAVE_XPETRA_TPETRA
#endif

} // Xpetra namespace


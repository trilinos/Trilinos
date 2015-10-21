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
#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template<class EpetraGlobalOrdinal, class Node>
  EpetraImportT<EpetraGlobalOrdinal, Node>::EpetraImportT(const Teuchos::RCP<const map_type > & source, const Teuchos::RCP<const map_type > & target)
    : import_(rcp(new Epetra_Import(toEpetra<EpetraGlobalOrdinal,Node>(target), toEpetra<EpetraGlobalOrdinal,Node>(source)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  // //! copy constructor.
  // EpetraImportT<EpetraGlobalOrdinal>::EpetraImportT(const Import<int,GlobalOrdinal> & import) { // TODO: refactoring
  //   XPETRA_DYNAMIC_CAST(const EpetraImportT, import, tImport, "Xpetra::EpetraImportT copy constructors only accept Xpetra::EpetraImportT as input arguments.");
  //   import_ = rcp(new Epetra_Import(*tImport.getEpetra_Import()));
  // }

  // TODO: move that elsewhere
  //   template<class GlobalOrdinal>
  //   const Epetra_Import & toEpetra(const Import<int, GlobalOrdinal> &import) {
  //     // TODO: throw exception
  //     const EpetraImportT & tpetraImport = dynamic_cast<const EpetraImportT<GlobalOrdinal> &>(import);
  //     return *tpetraImport.getEpetra_Import();
  //   }

  template<class GlobalOrdinal, class Node>
  RCP< const Import<int, GlobalOrdinal, Node > > toXpetra(const Epetra_Import *import) {
    if (import != NULL) {
      RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*import)); //NOTE: non consitent: return pointer, take ref
      return rcp(new Xpetra::EpetraImportT<GlobalOrdinal, Node>(imp));
    }

    return Teuchos::null;
  }
  //

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal, Node>::getExportPIDs() const { XPETRA_MONITOR("EpetraImportT::getExportImageIDs"); return ArrayView<const int> (import_->ExportPIDs(),import_->NumExportIDs()); }

  template<class EpetraGlobalOrdinal, class Node>
  size_t EpetraImportT<EpetraGlobalOrdinal, Node>::getNumRemoteIDs() const { XPETRA_MONITOR("EpetraImportT::getNumRemoteIDs"); return import_->NumRemoteIDs(); }
  template<class EpetraGlobalOrdinal, class Node>
  size_t EpetraImportT<EpetraGlobalOrdinal, Node>::getNumExportIDs() const { XPETRA_MONITOR("EpetraImportT::getNumExportIDs"); return import_->NumExportIDs(); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal, Node>::getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getPermuteFromLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getExportImageIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal, Node>::getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getPermuteToLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getPermuteToLIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal, Node>::getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getRemoteLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getRemoteLIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal, Node>::getRemotePIDs() const {
    XPETRA_MONITOR("EpetraImportT::getRemotePIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getRemotePIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal, Node>::getExportLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getExportLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getExportLIDs not implemented"); }

  template<class EpetraGlobalOrdinal, class Node>
  void EpetraImportT<EpetraGlobalOrdinal, Node>::print(std::ostream &os) const { XPETRA_MONITOR("EpetraImportT::print"); import_->Print(os); }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraImportT<int, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP<const Import<int,int,Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<int,Kokkos::Compat::KokkosSerialWrapperNode >(const Epetra_Import *);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraImportT<int, Kokkos::Compat::KokkosThreadsWrapperNode>;
template RCP<const Import<int,int,Kokkos::Compat::KokkosThreadsWrapperNode> > toXpetra<int,Kokkos::Compat::KokkosThreadsWrapperNode >(const Epetra_Import *);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraImportT<int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP<const Import<int, int, Kokkos::Compat::KokkosOpenMPWrapperNode> > toXpetra<int,Kokkos::Compat::KokkosOpenMPWrapperNode>(const Epetra_Import *);
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::View<int>::HostMirror::execution_space default_host_execution_space;
typedef Kokkos::Compat::KokkosDeviceWrapperNode<host_execution_space, Kokkos::HostSpace> default_node_type;
template class EpetraImportT<int, default_node_type >;
template RCP<const Import<int, int, default_node_type> > toXpetra<int,default_node_type>(const Epetra_Import *);
#endif
#else
  // TODO What, if Tpetra is disabled? Use fake Kokkos thing?
#endif // HAVE_XPETRA_TPETRA
#endif // XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_SERIAL
template class EpetraImportT<long long, Kokkos::Compat::KokkosSerialWrapperNode >;
template RCP<const Import<int,long long,Kokkos::Compat::KokkosSerialWrapperNode > > toXpetra<long long,Kokkos::Compat::KokkosSerialWrapperNode >(const Epetra_Import *);
#endif
#ifdef HAVE_XPETRA_PTHREAD
template class EpetraImportT<long long, Kokkos::Compat::KokkosThreadsWrapperNode>;
template RCP<const Import<int,long long,Kokkos::Compat::KokkosThreadsWrapperNode> > toXpetra<long long,Kokkos::Compat::KokkosThreadsWrapperNode >(const Epetra_Import *);
#endif
#ifdef HAVE_XPETRA_OPENMP
template class EpetraImportT<long, Kokkos::Compat::KokkosOpenMPWrapperNode >;
template RCP<const Import<int, long long, Kokkos::Compat::KokkosOpenMPWrapperNode> > toXpetra<long long,Kokkos::Compat::KokkosOpenMPWrapperNode>(const Epetra_Import *);
#endif
#ifdef HAVE_XPETRA_CUDA
typedef Kokkos::View<long long>::HostMirror::execution_space default_host_execution_space;
typedef Kokkos::Compat::KokkosDeviceWrapperNode<host_execution_space, Kokkos::HostSpace> default_node_type;
template class EpetraImportT<long long, default_node_type >;
template RCP<const Import<int, long long, default_node_type> > toXpetra<long long,default_node_type>(const Epetra_Import *);
#endif
#else
  // What, if Tpetra is disabled? Use fake Kokkos thing?
#endif // HAVE_XPETRA_TPETRA

// TODO
//template class EpetraImportT<long long, typename Xpetra::Map<int, long long>::node_type>;
//template RCP< const Import<int, long long> > toXpetra<long long>(const Epetra_Import *);
#endif

} // Xpetra namespace


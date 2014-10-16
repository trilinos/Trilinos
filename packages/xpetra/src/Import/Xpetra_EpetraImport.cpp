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

  template<class EpetraGlobalOrdinal>
  EpetraImportT<EpetraGlobalOrdinal>::EpetraImportT(const Teuchos::RCP<const map_type > & source, const Teuchos::RCP<const map_type > & target)
    : import_(rcp(new Epetra_Import(toEpetra(target), toEpetra(source)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

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

  template<class GlobalOrdinal>
  RCP< const Import<int, GlobalOrdinal > > toXpetra(const Epetra_Import *import) {
    if (import != NULL) {
      RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*import)); //NOTE: non consitent: return pointer, take ref
      return rcp(new Xpetra::EpetraImportT<GlobalOrdinal>(imp));
    }

    return Teuchos::null;
  }
  //

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal>::getExportPIDs() const { XPETRA_MONITOR("EpetraImportT::getExportImageIDs"); return ArrayView<const int> (import_->ExportPIDs(),import_->NumExportIDs()); }

  template<class EpetraGlobalOrdinal>
  size_t EpetraImportT<EpetraGlobalOrdinal>::getNumRemoteIDs() const { XPETRA_MONITOR("EpetraImportT::getNumRemoteIDs"); return import_->NumRemoteIDs(); }
  template<class EpetraGlobalOrdinal>
  size_t EpetraImportT<EpetraGlobalOrdinal>::getNumExportIDs() const { XPETRA_MONITOR("EpetraImportT::getNumExportIDs"); return import_->NumExportIDs(); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal>::getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getPermuteFromLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getExportImageIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal>::getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getPermuteToLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getPermuteToLIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal>::getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getRemoteLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getRemoteLIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal>::getRemotePIDs() const {
    XPETRA_MONITOR("EpetraImportT::getRemotePIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getRemotePIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraImportT<EpetraGlobalOrdinal>::getExportLIDs() const {
    XPETRA_MONITOR("EpetraImportT::getExportLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImportT<EpetraGlobalOrdinal>::getExportLIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  void EpetraImportT<EpetraGlobalOrdinal>::print(std::ostream &os) const { XPETRA_MONITOR("EpetraImportT::print"); import_->Print(os); }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template class EpetraImportT<int>;
template RCP< const Import<int, int> > toXpetra<int>(const Epetra_Import *);
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template class EpetraImportT<long long>;
template RCP< const Import<int, long long> > toXpetra<long long>(const Epetra_Import *);
#endif

} // Xpetra namespace


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

  template<class GlobalOrdinal>
  RCP<const Export<int, GlobalOrdinal> > toXpetra(const Epetra_Export *exp) {
    if (exp != NULL) {
      RCP<const Epetra_Export> eexp = rcp(new Epetra_Export(*exp)); //NOTE: non consitent: return pointer, take ref
      return rcp(new Xpetra::EpetraExportT<GlobalOrdinal>(eexp));
    }

    return Teuchos::null;
  }

  template<class EpetraGlobalOrdinal>
  EpetraExportT<EpetraGlobalOrdinal>::EpetraExportT(const Teuchos::RCP<const Map<int,GlobalOrdinal> > & source, const Teuchos::RCP<const Map<int,GlobalOrdinal> > & target)
    : export_(rcp(new Epetra_Export(toEpetra(source), toEpetra(target)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  //
  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal>::getExportPIDs() const { XPETRA_MONITOR("EpetraExportT::getExportImageIDs"); return ArrayView<const int> (export_->ExportPIDs(),export_->NumExportIDs()); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal>::getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getPermuteFromLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getExportImageIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal>::getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getPermuteToLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getPermuteToLIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  size_t EpetraExportT<EpetraGlobalOrdinal>::getNumRemoteIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumRemoteIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getNumRemoteIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal>::getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getRemoteLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getRemoteLIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  size_t EpetraExportT<EpetraGlobalOrdinal>::getNumExportIDs() const {
    XPETRA_MONITOR("EpetraExportT::getNumExportIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getNumExportIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  ArrayView< const int > EpetraExportT<EpetraGlobalOrdinal>::getExportLIDs() const {
    XPETRA_MONITOR("EpetraExportT::getExportLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExportT<EpetraGlobalOrdinal>::getExportLIDs not implemented"); }

  template<class EpetraGlobalOrdinal>
  void EpetraExportT<EpetraGlobalOrdinal>::print(std::ostream &os) const { XPETRA_MONITOR("EpetraExportT::");  export_->Print(os); }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template class EpetraExportT<int>;
template RCP<const Export<int, int> > toXpetra<int>(const Epetra_Export *);
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template class EpetraExportT<long long>;
template RCP<const Export<int, long long> > toXpetra<long long>(const Epetra_Export *);
#endif

} // Xpetra namespace


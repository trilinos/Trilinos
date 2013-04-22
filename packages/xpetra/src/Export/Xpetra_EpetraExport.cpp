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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  EpetraExport::EpetraExport(const Teuchos::RCP<const Map<int,int> > & source, const Teuchos::RCP<const Map<int,int> > & target)
    : export_(rcp(new Epetra_Export(toEpetra(source), toEpetra(target)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  //
  RCP< const Export<int, int > > toXpetra(const Epetra_Export *import) {
    RCP<const Epetra_Export> imp = rcp(new Epetra_Export(*import)); //NOTE: non consitent: return pointer, take ref
    return rcp ( new Xpetra::EpetraExport(imp) );
  }
  //

  ArrayView< const int > EpetraExport::getExportPIDs() const { XPETRA_MONITOR("EpetraExport::getExportImageIDs"); return ArrayView<const int> (export_->ExportPIDs(),export_->NumExportIDs()); }

  ArrayView< const int > EpetraExport::getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraExport::getPermuteFromLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getExportImageIDs not implemented"); }

  ArrayView< const int > EpetraExport::getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraExport::getPermuteToLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getPermuteToLIDs not implemented"); }

  size_t EpetraExport::getNumRemoteIDs() const {
    XPETRA_MONITOR("EpetraExport::getNumRemoteIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getNumRemoteIDs not implemented"); }

  ArrayView< const int > EpetraExport::getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraExport::getRemoteLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getRemoteLIDs not implemented"); }

  size_t EpetraExport::getNumExportIDs() const {
    XPETRA_MONITOR("EpetraExport::getNumExportIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getNumExportIDs not implemented"); }

  ArrayView< const int > EpetraExport::getExportLIDs() const {
    XPETRA_MONITOR("EpetraExport::getExportLIDs");
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getExportLIDs not implemented"); }

  void EpetraExport::print(std::ostream &os) const { XPETRA_MONITOR("EpetraExport::");  export_->Print(os); }

} // Xpetra namespace


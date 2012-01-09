#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  EpetraExport::EpetraExport(const Teuchos::RCP<const Map<int,int> > & source, const Teuchos::RCP<const Map<int,int> > & target)
    : export_(rcp(new Epetra_Export(toEpetra(source), toEpetra(target)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  RCP< const Export<int, int > > toXpetra(const Epetra_Export *import) {
    RCP<const Epetra_Export> imp = rcp(new Epetra_Export(*import)); //NOTE: non consitent: return pointer, take ref
    return rcp ( new Xpetra::EpetraExport(imp) );
  }
  //

  ArrayView< const int > EpetraExport::getExportImageIDs() const { return ArrayView<const int> (export_->ExportPIDs(),export_->NumExportIDs()); }

  ArrayView< const int > EpetraExport::getPermuteFromLIDs() const {
         TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getExportImageIDs not implemented"); }

  ArrayView< const int > EpetraExport::getPermuteToLIDs() const {
         TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getPermuteToLIDs not implemented"); }

  size_t EpetraExport::getNumRemoteIDs() const {
         TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getNumRemoteIDs not implemented"); }

  ArrayView< const int > EpetraExport::getRemoteLIDs() const {
         TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getRemoteLIDs not implemented"); }

  size_t EpetraExport::getNumExportIDs() const {
         TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getNumExportIDs not implemented"); }

  ArrayView< const int > EpetraExport::getExportLIDs() const {
         TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraExport::getExportLIDs not implemented"); }

  void EpetraExport::print(std::ostream &os) const {export_->Print(os);}

} // Xpetra namespace


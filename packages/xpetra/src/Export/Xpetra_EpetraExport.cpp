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

  ArrayView< const int > EpetraExport::getExportImageIDs() const { XPETRA_MONITOR("EpetraExport::getExportImageIDs"); return ArrayView<const int> (export_->ExportPIDs(),export_->NumExportIDs()); }

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


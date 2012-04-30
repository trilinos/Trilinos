#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  EpetraImport::EpetraImport(const Teuchos::RCP<const Map<int,int> > & source, const Teuchos::RCP<const Map<int,int> > & target)
    : import_(rcp(new Epetra_Import(toEpetra(target), toEpetra(source)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)
  
  // //! copy constructor. 
  // EpetraImport::EpetraImport(const Import<int,int> & import) { // TODO: refactoring
  //   XPETRA_DYNAMIC_CAST(const EpetraImport, import, tImport, "Xpetra::EpetraImport copy constructors only accept Xpetra::EpetraImport as input arguments.");
  //   import_ = rcp(new Epetra_Import(*tImport.getEpetra_Import()));
  // }
  
  // TODO: move that elsewhere
  //   const Epetra_Import & toEpetra(const Import<int, int> &import) {
  //     // TODO: throw exception
  //     const EpetraImport & tpetraImport = dynamic_cast<const EpetraImport &>(import);
  //     return *tpetraImport.getEpetra_Import();
  //   }

  RCP< const Import<int, int > > toXpetra(const Epetra_Import *import) {
    RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*import)); //NOTE: non consitent: return pointer, take ref
    return rcp ( new Xpetra::EpetraImport(imp) );
  }
  //

  ArrayView< const int > EpetraImport::getExportImageIDs() const { XPETRA_MONITOR("EpetraImport::getExportImageIDs"); return ArrayView<const int> (import_->ExportPIDs(),import_->NumExportIDs()); }

  ArrayView< const int > EpetraImport::getPermuteFromLIDs() const {
    XPETRA_MONITOR("EpetraImport::getPermuteFromLIDs"); 
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImport::getExportImageIDs not implemented"); }

  ArrayView< const int > EpetraImport::getPermuteToLIDs() const {
    XPETRA_MONITOR("EpetraImport::getPermuteToLIDs"); 
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImport::getPermuteToLIDs not implemented"); }

  size_t EpetraImport::getNumRemoteIDs() const {
    XPETRA_MONITOR("EpetraImport::getNumRemoteIDs"); 
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImport::getNumRemoteIDs not implemented"); }

  ArrayView< const int > EpetraImport::getRemoteLIDs() const {
    XPETRA_MONITOR("EpetraImport::getRemoteLIDs"); 
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImport::getRemoteLIDs not implemented"); }

  size_t EpetraImport::getNumExportIDs() const {
    XPETRA_MONITOR("EpetraImport::getNumExportIDs"); 
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImport::getNumExportIDs not implemented"); }

  ArrayView< const int > EpetraImport::getExportLIDs() const {
    XPETRA_MONITOR("EpetraImport::getExportLIDs"); 
    TEUCHOS_TEST_FOR_EXCEPTION(1, Xpetra::Exceptions::NotImplemented, "TODO EpetraImport::getExportLIDs not implemented"); }

  void EpetraImport::print(std::ostream &os) const { XPETRA_MONITOR("EpetraImport::print"); import_->Print(os); }


} // Xpetra namespace


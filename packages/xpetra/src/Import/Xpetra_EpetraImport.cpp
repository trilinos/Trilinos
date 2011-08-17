#include "Cthulhu_EpetraImport.hpp"

namespace Cthulhu {

  EpetraImport::EpetraImport(const Teuchos::RCP<const Map<int,int> > & source, const Teuchos::RCP<const Map<int,int> > & target)
    : import_(rcp(new Epetra_Import(toEpetra(target), toEpetra(source)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)
  
  // //! copy constructor. 
  // EpetraImport::EpetraImport(const Import<int,int> & import) { // TODO: refactoring
  //   CTHULHU_DYNAMIC_CAST(const EpetraImport, import, tImport, "Cthulhu::EpetraImport copy constructors only accept Cthulhu::EpetraImport as input arguments.");
  //   import_ = rcp(new Epetra_Import(*tImport.getEpetra_Import()));
  // }
  
  // TODO: move that elsewhere
  //   const Epetra_Import & toEpetra(const Import<int, int> &import) {
  //     // TODO: throw exception
  //     const EpetraImport & tpetraImport = dynamic_cast<const EpetraImport &>(import);
  //     return *tpetraImport.getEpetra_Import();
  //   }

  RCP< const Import<int, int > > toCthulhu(const Epetra_Import *import) {
    RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*import)); //NOTE: non consitent: return pointer, take ref
    return rcp ( new Cthulhu::EpetraImport(imp) );
  }
  //

} // Cthulhu namespace


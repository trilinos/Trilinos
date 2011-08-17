#include "Xpetra_EpetraImport.hpp"

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

} // Xpetra namespace


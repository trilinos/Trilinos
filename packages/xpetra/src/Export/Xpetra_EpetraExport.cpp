#include "Xpetra_EpetraExport.hpp"

namespace Xpetra {

  EpetraExport::EpetraExport(const Teuchos::RCP<const Map<int,int> > & source, const Teuchos::RCP<const Map<int,int> > & target)
    : export_(rcp(new Epetra_Export(toEpetra(target), toEpetra(source)))) { } // Warning: Epetra(Target, Source) vs. Tpetra(Source, Target)

  RCP< const Export<int, int > > toXpetra(const Epetra_Export *import) {
    RCP<const Epetra_Export> imp = rcp(new Epetra_Export(*import)); //NOTE: non consitent: return pointer, take ref
    return rcp ( new Xpetra::EpetraExport(imp) );
  }
  //

} // Xpetra namespace


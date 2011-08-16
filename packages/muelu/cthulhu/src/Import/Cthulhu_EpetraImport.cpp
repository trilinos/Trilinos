#include <Cthulhu_EpetraImport.hpp>

namespace Cthulhu {

  // TODO: move that elsewhere
//   template <class LocalOrdinal, class GlobalOrdinal, class Node>
//   const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & toTpetra(const Import<LocalOrdinal,GlobalOrdinal,Node> &import) {
//     // TODO: throw exception
//     const TpetraImport<LocalOrdinal,GlobalOrdinal,Node> & tpetraImport = dynamic_cast<const TpetraImport<LocalOrdinal,GlobalOrdinal,Node> &>(import);
//     return *tpetraImport.getTpetra_Import();
//   }

// RCP< const Import<int, int > > toCthulhu(const Epetra_Import &import) {
//   RCP<const Epetra_Import> imp = rcp(new Epetra_Import(import)); //NOTE: non consitent: return pointer, take ref
//   return rcp ( new Cthulhu::EpetraImport(imp) );
// }

RCP< const Import<int, int > > toCthulhu(const Epetra_Import *import) {
  RCP<const Epetra_Import> imp = rcp(new Epetra_Import(*import)); //NOTE: non consitent: return pointer, take ref
  return rcp ( new Cthulhu::EpetraImport(imp) );
}
//

} // Cthulhu namespace

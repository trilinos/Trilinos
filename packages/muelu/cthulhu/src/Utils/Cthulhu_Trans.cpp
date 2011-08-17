#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#include "Cthulhu_Trans.hpp"

namespace Cthulhu {
  
  //! Convert a Teuchos::ETransp to a boolean (for Epetra).
  bool toEpetra(Teuchos::ETransp trans) { 
    if (trans == Teuchos::NO_TRANS)
      return false;
    else if (trans == Teuchos::TRANS)
      return true;
    else { 
      TEST_FOR_EXCEPTION((trans != Teuchos::NO_TRANS) && (trans == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cannot convert Teuchos::ETransp to a boolean.");
    }
    
    return false; // to skip compilation warning msg.
  }

} // Cthulhu namespace

#endif

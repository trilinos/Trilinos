#ifndef CTHULHU_EPETRATRANS_HPP
#define CTHULHU_EPETRATRANS_HPP

#include <Teuchos_BLAS_types.hpp>

//   enum ETransp { 	
//     NO_TRANS,	/*!< Not transposed */ 
//     TRANS, 		/*!< Transposed */
//     CONJ_TRANS 	/*!< Conjugate transposed */
//   };

#include "Cthulhu_Exceptions.hpp"
#include "Cthulhu_Debug.hpp"

//! Convert a Teuchos::ETransp to a boolean (for Epetra).
bool Teuchos2Epetra_Trans(Teuchos::ETransp trans) { CTHULHU_DEBUG_ME;
  if (trans == Teuchos::NO_TRANS)
    return false;
  else if (trans == Teuchos::TRANS)
    return true;
  else { 
    TEST_FOR_EXCEPTION((trans != Teuchos::NO_TRANS) && (trans == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cannot convert Teuchos::ETransp to a boolean.");
  }

  return false; // to skip compilation warning msg.
}

#endif // CTHULHU_EPETRATRANS_HPP

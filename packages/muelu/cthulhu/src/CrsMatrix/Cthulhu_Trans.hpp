#ifndef CTHULHU_TRANS_HPP 
#define CTHULHU_TRANS_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#include <Teuchos_BLAS_types.hpp>

#include "Cthulhu_Exceptions.hpp"

//! Convert a Teuchos::ETransp to a boolean (for Epetra).
bool Teuchos2Epetra_Trans(Teuchos::ETransp trans);

#endif // HAVE_CTHULHU_EPETRA

#endif // CTHULHU_EPETRATRANS_HPP

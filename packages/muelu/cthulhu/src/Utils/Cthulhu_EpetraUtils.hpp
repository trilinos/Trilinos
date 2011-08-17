#ifndef CTHULHU_COMM_HPP
#define CTHULHU_COMM_HPP

//! Conversion between Epetra and Teuchos objects

#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

// header files for comm objects conversion
#include <Teuchos_Comm.hpp>
#include <Epetra_Comm.h>

// header file for Teuchos::ETransp
#include <Teuchos_BLAS_types.hpp>

namespace Cthulhu {

  using Teuchos::RCP;
  
  //! Convert a Teuchos_Comm to an Epetra_Comm.
  const RCP<const Epetra_Comm> toEpetra(const RCP<const Teuchos::Comm<int> > & comm);
  
  //! Convert an Epetra_Comm.to a Teuchos_Comm.
  const RCP<const Teuchos::Comm<int> > toCthulhu(const Epetra_Comm & comm);

  //! Convert a Teuchos::ETransp to an Epetra boolean.
  bool toEpetra(Teuchos::ETransp);
  
}
#endif // HAVE_CTHULHU_EPETRA

#endif // CTHULHU_EPETRACOMM_HPP

// Note: no support for Epetra_MpiSmpComm
// TODO: remove return RCP for toEpetra?

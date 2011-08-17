#ifndef CTHULHU_COMM_HPP
#define CTHULHU_COMM_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#ifdef HAVE_MPI
#include <mpi.h>
#include <Teuchos_OpaqueWrapper.hpp>
#endif

#include <Teuchos_Comm.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#endif
#include <Teuchos_DefaultSerialComm.hpp>

#include <Epetra_Comm.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>

#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {

  using Teuchos::RCP;
  
  //! Convert a Teuchos_Comm to an Epetra_Comm.
  const RCP<const Epetra_Comm> Teuchos2Epetra_Comm(const RCP<const Teuchos::Comm<int> > & comm);
  
  //! Convert an Epetra_Comm.to a Teuchos_Comm
  const RCP<const Teuchos::Comm<int> > Epetra2Teuchos_Comm(RCP<const Epetra_Comm> & comm);

  const RCP<const Teuchos::Comm<int> > toCthulhu(const Epetra_Comm & comm);
  
}
#endif // HAVE_CTHULHU_EPETRA

#endif // CTHULHU_EPETRACOMM_HPP

// JG Note: At present, there is no support for Epetra_MpiSmpComm but if
// needed, Cthulhu object XXX with such internal Comm can be created
// using constructors Cthulhu::EpetraXXX(myEpetra_Map) which
// directly wrap an Epetra_XXX.

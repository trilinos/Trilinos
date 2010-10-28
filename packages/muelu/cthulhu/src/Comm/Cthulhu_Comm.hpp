#ifndef CTHULHU_EPETRACOMM_HPP
#define CTHULHU_EPETRACOMM_HPP

#ifdef HAVE_MPI
#include <mpi.h>
#include <Teuchos_OpaqueWrapper.hpp>
#endif

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include "Cthulhu_Exceptions.hpp"
#include "Cthulhu_Debug.hpp"

using Teuchos::RCP;

//! Convert a Teuchos_Comm to an Epetra_Comm.
const RCP<const Epetra_Comm> Teuchos2Epetra_Comm(const RCP<const Teuchos::Comm<int> > & comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_MPI
  const RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
  if (mpiComm != Teuchos::null) {
    return Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()));
  }  else
#endif
    if ((Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<int> >(comm) != Teuchos::null))
      return Teuchos::rcp(new Epetra_SerialComm());
    else
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot convert a Teuchos::Comm to an Epetra_Comm: The exact type of the Teuchos::Comm object is unknown"); 
}

//! Convert an Epetra_Comm.to a Teuchos_Comm
const RCP<const Teuchos::Comm<int> > Epetra2Teuchos_Comm(RCP<const Epetra_Comm> & comm) { CTHULHU_DEBUG_ME;
#ifdef HAVE_MPI
  const RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(comm);
  if (mpiComm != Teuchos::null) 
    return Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm->Comm()))); 
  else 
#endif
    if (Teuchos::rcp_dynamic_cast<const Epetra_SerialComm>(comm) != Teuchos::null)
      return Teuchos::rcp(new Teuchos::SerialComm<int>());
    else
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot convert an Epetra_Comm to a Teuchos::Comm: The exact type of the Epetra_Comm object is unknown");
}

#endif // CTHULHU_EPETRACOMM_HPP

// JG Note: At present, there is no support for Epetra_MpiSmpComm but if
// needed, Cthulhu object XXX with such internal Comm can be created
// using constructors Cthulhu::EpetraXXX(myEpetra_Map) which
// directly wrap an Epetra_XXX.

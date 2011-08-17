#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#include "Cthulhu_EpetraUtils.hpp"

// header files for comm objects conversion
#ifdef HAVE_MPI
#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#endif
#include <Teuchos_DefaultSerialComm.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>

#include "Cthulhu_Exceptions.hpp"


namespace Cthulhu {

  using Teuchos::RCP;

  const RCP<const Epetra_Comm> toEpetra(const RCP<const Teuchos::Comm<int> > & comm) { 
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

  const RCP<const Teuchos::Comm<int> > toCthulhu(const Epetra_Comm & comm) {
#ifdef HAVE_MPI
    try {
      const Epetra_MpiComm& mpiComm = dynamic_cast<const Epetra_MpiComm&>(comm);
      return Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm.Comm()))); 
    } catch (std::bad_cast & b) {}
#endif
    try {
      const Epetra_SerialComm& serialComm = dynamic_cast<const Epetra_SerialComm&>(comm);
      serialComm.NumProc(); // avoid compilation warning
      return Teuchos::rcp(new Teuchos::SerialComm<int>());
    } catch (std::bad_cast & b) {
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot convert an Epetra_Comm to a Teuchos::Comm: The exact type of the Epetra_Comm object is unknown");
    }
  }

  bool toEpetra(Teuchos::ETransp trans) { 
    if (trans == Teuchos::NO_TRANS)
      return false;
    else if (trans == Teuchos::TRANS)
      return true;
    else { 
      TEST_FOR_EXCEPTION((trans != Teuchos::NO_TRANS) && (trans == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cannot convert Teuchos::ETransp to a boolean.");
    }
    
    return false; // to skip a compilation warning msg.
  }

}

#endif

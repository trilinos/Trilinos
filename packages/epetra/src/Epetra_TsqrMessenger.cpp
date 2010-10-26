#include <Epetra_TsqrMessenger.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_KOKKOS_TSQR

namespace TSQR {
  namespace Epetra { 

#ifdef EPETRA_MPI
    std::pair< MPI_Comm, bool >
    extractRawMpiComm (const Teuchos::RCP< const Epetra_Comm >& pComm)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;

      MPI_Comm rawMpiComm = MPI_COMM_NULL;
      bool haveMpiComm = false;
      
      RCP< const Epetra_MpiComm > pMpiComm = 
	rcp_dynamic_cast< const Epetra_MpiComm > (pComm, false);
      if (pMpiComm.get() == NULL)
	{
	  // Teuchos::rcp_dynamic_cast doesn't seem to work with
	  // Epetra_MpiSmpComm.
#if 0 
	  // See if the input Epetra_Comm is really an
	  // Epetra_MpiSmpComm.  If so, pull out its raw MPI_Comm
	  // object and use that.
	  RCP< const Epetra_MpiSmpComm > pMpiSmpComm = 
	    rcp_dynamic_cast< const Epetra_MpiSmpComm > (pComm, false);
	  if (pMpiSmpComm.get() != NULL)
	    {
	      rawMpiComm = pMpiSmpComm->Comm();
	      haveMpiComm = true;
	    }
#else
	  haveMpiComm = false;
#endif // 0
	}
      else
	{
	  rawMpiComm = pMpiComm->Comm();
	  haveMpiComm = true;
	}
      return std::make_pair (rawMpiComm, haveMpiComm);
    }
#endif // EPETRA_MPI

  } // namespace Epetra
} // namespace TSQR

#endif // HAVE_KOKKOS_TSQR

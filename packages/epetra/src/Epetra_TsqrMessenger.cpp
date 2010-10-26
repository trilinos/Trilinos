#include <Epetra_TsqrMessenger.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Epetra { 

    // We only build this routine if Epetra was built with MPI
    // support.  Otherwise, we don't even have an MPI_Comm definition.
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
	haveMpiComm = false;
      else
	{
	  rawMpiComm = pMpiComm->Comm();
	  haveMpiComm = true;
	}
      return std::make_pair (rawMpiComm, haveMpiComm);
    }
#else
#  error BADNESS IN THE EXTREME!
#endif // EPETRA_MPI

  } // namespace Epetra
} // namespace TSQR


#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultSpmdVectorSpace_decl.hpp"
#include "Thyra_DefaultSpmdVector_decl.hpp"
#include "Thyra_MultiVectorBase_decl.hpp"
#include "Thyra_ScalarProdVectorSpaceBase_decl.hpp"
#include "Thyra_DefaultSpmdMultiVector_decl.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include <iostream> // std::cerr, std::endl
#include <sstream>
#include <string>

// If Thyra is compiled with TEUCHOS_DEBUG defined then the following
// wil hang in a collective MPI communication when run on four
// processors.
TEUCHOS_UNIT_TEST( ThyraEpetraMultiVector, HangingInParallelDebug )
{
   using std::cerr;
   using std::endl;

   int myRank = 0;
   int numProcs = 1;
#ifdef HAVE_MPI
   (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
   (void) MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
#endif // HAVE_MPI

   std::string prefix;
   {
     std::ostringstream os;
     os << "(Process " << myRank << ") ";
     prefix = os.str ();
   }

   {
     std::ostringstream os;
     os << prefix << "Creating Epetra_Comm" << endl;
     cerr << os.str ();
   }
#ifdef HAVE_MPI
   const Epetra_MpiComm epetra_comm (MPI_COMM_WORLD);
#else
   const Epetra_SerialComm epetra_comm ();
#endif
   {
     std::ostringstream os;
     os << prefix << "Creating Teuchos::Comm" << endl;
     cerr << os.str ();
   }
   Teuchos::RCP<const Teuchos::Comm<Teuchos_Ordinal> > comm =
     Teuchos::DefaultComm<Teuchos_Ordinal>::getComm ();

   // Some processors have to have data and some not.
   const int localDim  = myRank % 2;
   const int globalDim = numProcs / 2;
   Teuchos::RCP<const Epetra_Map> epetra_map;
   {
     std::ostringstream os;
     os << prefix << "Creating Epetra_Map: localDim=" << localDim << ", globalDim=" << globalDim << endl;
     cerr << os.str ();
   }
   epetra_map = Teuchos::rcp(new Epetra_Map(globalDim,localDim,0,epetra_comm));
   {
     std::ostringstream os;
     os << prefix << "Creating Thyra::DefaultSpmdVectorSpace" << endl;
     cerr << os.str ();
   }
   Teuchos::RCP<Thyra::DefaultSpmdVectorSpace<double>> SPMD = Thyra::DefaultSpmdVectorSpace<double>::create();
   SPMD->initialize(comm, localDim, globalDim);
   {
     std::ostringstream os;
     os << prefix << "Creating Thyra::MultiVectorBase" << endl;
     cerr << os.str ();
   }
   Teuchos::RCP<const Thyra::MultiVectorBase<double>> spmd =
      Teuchos::rcp(
         new Thyra::DefaultSpmdMultiVector<double>(
            SPMD,
            Teuchos::rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<double>>(
               SPMD->smallVecSpcFcty()->createVecSpc(1),true)
            )
        );
   {
     std::ostringstream os;
     os << prefix << "Creating Thyra::get_Epetra_MultiVector" << endl;
     cerr << os.str ();
   }
   Thyra::get_Epetra_MultiVector(*epetra_map,*spmd);
   {
     std::ostringstream os;
     os << prefix << "Done with test on this process" << endl;
     cerr << os.str ();
   }
}


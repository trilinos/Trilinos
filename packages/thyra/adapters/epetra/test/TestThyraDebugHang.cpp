// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
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
   using Teuchos::outArg;
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::REDUCE_MIN;
   using Teuchos::reduceAll;
   using std::cerr;
   using std::endl;
   int lclSuccess = 1; // to be revised below
   int gblSuccess = 1; // to be revised below
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
   RCP<const Teuchos::Comm<Teuchos_Ordinal> > comm =
     Teuchos::DefaultComm<Teuchos_Ordinal>::getComm ();
   // Make sure that everything is OK on all processes.
   TEST_ASSERT( ! comm.is_null () );
   lclSuccess = success ? 1 : 0;
   reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
   TEST_EQUALITY( gblSuccess, 1 );
   if (gblSuccess != 1) {
     out << "FAILED; some process(es) have a null Teuchos::Comm" << endl;
     return;
   }

   // Some processors have to have data and some not.
   const int localDim  = myRank % 2;
   const int globalDim = numProcs / 2;
   RCP<const Epetra_Map> epetra_map;
   {
     std::ostringstream os;
     os << prefix << "Creating Epetra_Map: localDim=" << localDim << ", globalDim=" << globalDim << endl;
     cerr << os.str ();
   }
   epetra_map = rcp (new Epetra_Map (globalDim, localDim, 0, epetra_comm));
   // Make sure that everything is OK on all processes.
   TEST_ASSERT( ! epetra_map.is_null () );
   lclSuccess = success ? 1 : 0;
   reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
   TEST_EQUALITY( gblSuccess, 1 );
   if (gblSuccess != 1) {
     out << "FAILED; some process(es) have a null Epetra_Map" << endl;
     return;
   }

   {
     std::ostringstream os;
     os << prefix << "Creating Thyra::DefaultSpmdVectorSpace" << endl;
     cerr << os.str ();
   }
   RCP<Thyra::DefaultSpmdVectorSpace<double> > SPMD =
     Thyra::DefaultSpmdVectorSpace<double>::create();

   // Make sure that everything is OK on all processes.
   TEST_ASSERT( ! epetra_map.is_null () );
   lclSuccess = success ? 1 : 0;
   reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
   TEST_EQUALITY( gblSuccess, 1 );
   if (gblSuccess != 1) {
     out << "FAILED; some process(es) have a null SPMD" << endl;
     return;
   }

   SPMD->initialize(comm, localDim, globalDim);

   {
     std::ostringstream os;
     os << prefix << "Creating Thyra::MultiVectorBase" << endl;
     cerr << os.str ();
   }
   RCP<const Thyra::MultiVectorBase<double> > spmd =
      rcp (new Thyra::DefaultSpmdMultiVector<double> (
        SPMD,
        rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<double> > (
          SPMD->smallVecSpcFcty()->createVecSpc(1),true)
        )
      );
   // Make sure that everything is OK on all processes.
   TEST_ASSERT( ! spmd.is_null () );
   lclSuccess = success ? 1 : 0;
   reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
   TEST_EQUALITY( gblSuccess, 1 );
   if (gblSuccess != 1) {
     out << "FAILED; some process(es) have a null Thyra::MultiVectorBase"
         << endl;
     return;
   }

   {
     std::ostringstream os;
     os << prefix << "Calling Thyra::get_Epetra_MultiVector "
       "(const overload; see #1941)" << endl;
     cerr << os.str ();
   }
   // Make sure that we invoke the const overload.
   RCP<const Epetra_MultiVector> mv_c =
     Thyra::get_Epetra_MultiVector (*epetra_map,
       const_cast<const Thyra::MultiVectorBase<double>& > (*spmd));
   // Make sure that everything is OK on all processes.
   TEST_ASSERT( ! mv_c.is_null () );
   lclSuccess = success ? 1 : 0;
   reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
   TEST_EQUALITY( gblSuccess, 1 );
   if (gblSuccess != 1) {
     out << "FAILED; some process(es) have a null const Epetra_MultiVector"
         << endl;
     return;
   }

   {
     std::ostringstream os;
     os << prefix << "Calling Thyra::get_Epetra_MultiVector "
       "(nonconst overload; see #2061)" << endl;
     cerr << os.str ();
   }
   // Make sure that we invoke the nonconst overload.
   RCP<Epetra_MultiVector> mv_nc =
     Thyra::get_Epetra_MultiVector (*epetra_map,
       const_cast<Thyra::MultiVectorBase<double>& > (*spmd));
   // Make sure that everything is OK on all processes.
   TEST_ASSERT( ! mv_nc.is_null () );
   lclSuccess = success ? 1 : 0;
   reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
   TEST_EQUALITY( gblSuccess, 1 );
   if (gblSuccess != 1) {
     out << "FAILED; some process(es) have a null nonconst Epetra_MultiVector"
         << endl;
     return;
   }

   {
     std::ostringstream os;
     os << prefix << "Done with test on this process" << endl;
     cerr << os.str ();
   }
}


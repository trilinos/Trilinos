// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Phalanx_KokkosUtilities.hpp"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "Thyra_EpetraThyraWrappers.hpp"

#include <vector>

namespace panzer {
  
  TEUCHOS_UNIT_TEST(epetra_test, maptest)
  {
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    int veca[] = {0, 1, 2, 3, 4, 5};
    int vecb[] = {6, 7, 8, 3, 4, 5, 0, 1, 2, 9, 10, 11, 15, 16, 17, 12, 13, 14, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    int vecc[] = {8, 7, 6, 9, 11, 10, 12, 13, 14};
    int vecd[] = {27};
    int empty = 0;

    int numPs = comm.NumProc();
    int myPid = comm.MyPID();

    int num0, num1 = 0;
    int * ind0, * ind1;
    switch(myPid) {
    case 0:

// 0 passes, 1 fails
#if 1 
      num0 = 0;
      ind0 = &empty;
#else
      num0 = 1;
      ind0 = vecd;
#endif

      num1 = 6;
      ind1 = veca;
      break;
    case 1:
      num0 = 6;
      ind0 = vecb;

      num1 = 9;
      ind1 = vecc;
      break;
    }

    Teuchos::RCP<Epetra_Map> map0 = Teuchos::rcp(new Epetra_Map(-1,num0,ind0,0,comm));
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs0 = Thyra::create_VectorSpace(map0);

    Teuchos::RCP<Epetra_Map> map1 = Teuchos::rcp(new Epetra_Map(-1,num1,ind1,0,comm));
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs1 = Thyra::create_VectorSpace(map1);
  }

}

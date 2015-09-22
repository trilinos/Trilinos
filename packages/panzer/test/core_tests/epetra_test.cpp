// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
    PHX::KokkosDeviceSession session;
    
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

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
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_config.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_MpiComm.h"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

TEUCHOS_UNIT_TEST(tEpetra_GlbEvalData, basic)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Thyra::VectorBase<double> Thyra_Vector;
  typedef Thyra::SpmdVectorBase<double> Thyra_SpmdVec;
  typedef Thyra::SpmdVectorSpaceBase<double> Thyra_SpmdVecSpace;

  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // This is required
  TEST_ASSERT(comm.NumProc()==2);

  std::vector<int> ghosted(5);
  std::vector<int> unique(3);

  if(comm.MyPID()==0) {
    unique[0] = 0;
    unique[1] = 1;
    unique[2] = 2;

    ghosted[0] = 0;
    ghosted[1] = 1;
    ghosted[2] = 2;
    ghosted[3] = 3;
    ghosted[4] = 4;
  }
  else {
    unique[0] = 3;
    unique[1] = 4;
    unique[2] = 5;

    ghosted[0] = 1;
    ghosted[1] = 2;
    ghosted[2] = 3;
    ghosted[3] = 4;
    ghosted[4] = 5;
  }

  RCP<const Epetra_Map> uniqueMap = rcp(new Epetra_Map(-1,unique.size(),&unique[0],0,comm));
  RCP<const Epetra_Map> ghostedMap = rcp(new Epetra_Map(-1,ghosted.size(),&ghosted[0],0,comm));
  RCP<const Epetra_Import> importer = rcp(new Epetra_Import(*ghostedMap,*uniqueMap));

  EpetraVector_ReadOnly_GlobalEvaluationData ged;
 
  TEST_ASSERT(!ged.isInitialized());

  ged.initialize(importer,ghostedMap,uniqueMap);

  TEST_ASSERT(ged.isInitialized());

  // test the ghosted vector sizing (we don't care what the entries are!)
  { 
    RCP<Epetra_Vector> ghostedVecE = ged.getGhostedVector_Epetra();
    RCP<Thyra_Vector>  ghostedVecT = ged.getGhostedVector();

    TEST_ASSERT(ghostedVecE!=Teuchos::null); 
    TEST_ASSERT(ghostedVecT!=Teuchos::null); 

    RCP<const Thyra::SpmdVectorSpaceBase<double> > ghostedSpace 
        = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(ghostedVecT->space());
    
    TEST_EQUALITY(ghostedMap->NumMyElements(),ghostedVecE->MyLength());
    TEST_EQUALITY(ghostedMap->NumGlobalElements(),ghostedVecE->GlobalLength());

    TEST_EQUALITY(ghostedSpace->isLocallyReplicated(),false);
    TEST_EQUALITY(ghostedSpace->localSubDim(),ghostedVecE->MyLength());
  }

  // test setting a unique vector
  {
    RCP<Epetra_Vector> uniqueVec = rcp(new Epetra_Vector(*uniqueMap));

    if(comm.MyPID()==0) {
      (*uniqueVec)[0] = 3.14;
      (*uniqueVec)[1] = 1.82;
      (*uniqueVec)[2] = -.91;
    }
    else {
      (*uniqueVec)[0] = 2.72;
      (*uniqueVec)[1] = 6.23;
      (*uniqueVec)[2] = -.17;
    }

    // set the unique vector, assure that const can be used
    ged.setUniqueVector_Epetra(uniqueVec.getConst());
  }

  // test the unique vector sizing and thyra entries
  { 
    const Epetra_Vector & uniqueVecE = *ged.getUniqueVector_Epetra();
    RCP<const Thyra_Vector>  uniqueVecT = ged.getUniqueVector();

    TEST_ASSERT(uniqueVecT!=Teuchos::null); 

    RCP<const Thyra::SpmdVectorSpaceBase<double> > uniqueSpace 
        = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(uniqueVecT->space());
    
    TEST_EQUALITY(uniqueMap->NumMyElements(),uniqueVecE.MyLength());
    TEST_EQUALITY(uniqueMap->NumGlobalElements(),uniqueVecE.GlobalLength());

    TEST_EQUALITY(uniqueSpace->isLocallyReplicated(),false);
    TEST_EQUALITY(uniqueSpace->localSubDim(),uniqueVecE.MyLength());

    RCP<const Thyra_SpmdVec> spmdVec = rcp_dynamic_cast<const Thyra_SpmdVec>(uniqueVecT);

    Teuchos::ArrayRCP<const double> thyraVec;
    spmdVec->getLocalData(Teuchos::ptrFromRef(thyraVec));

    TEST_EQUALITY(thyraVec.size(),uniqueVecE.MyLength());

    if(comm.MyPID()==0) {
      TEST_EQUALITY(uniqueVecE[0],3.14);
      TEST_EQUALITY(uniqueVecE[1],1.82);
      TEST_EQUALITY(uniqueVecE[2],-.91);
    }
    else {
      TEST_EQUALITY(uniqueVecE[0],2.72);
      TEST_EQUALITY(uniqueVecE[1],6.23);
      TEST_EQUALITY(uniqueVecE[2],-.17);
    }

    TEST_EQUALITY(uniqueVecE[0],thyraVec[0]);
    TEST_EQUALITY(uniqueVecE[1],thyraVec[1]);
    TEST_EQUALITY(uniqueVecE[2],thyraVec[2]);
  }

  // actually do something...
  ged.initializeData();
  ged.globalToGhost(0);

  {
    const Epetra_Vector & ghostedVecE = *ged.getGhostedVector_Epetra();
    RCP<Thyra_Vector>  ghostedVecT = ged.getGhostedVector();

    RCP<const Thyra_SpmdVec> spmdVec = rcp_dynamic_cast<const Thyra_SpmdVec>(ghostedVecT);

    Teuchos::ArrayRCP<const double> thyraVec;
    spmdVec->getLocalData(Teuchos::ptrFromRef(thyraVec));

    TEST_EQUALITY(thyraVec.size(),ghostedVecE.MyLength());

    if(comm.MyPID()==0) {
      TEST_EQUALITY(ghostedVecE[0],3.14);
      TEST_EQUALITY(ghostedVecE[1],1.82);
      TEST_EQUALITY(ghostedVecE[2],-.91);
      TEST_EQUALITY(ghostedVecE[3],2.72);
      TEST_EQUALITY(ghostedVecE[4],6.23);
    }
    else {
      TEST_EQUALITY(ghostedVecE[0],1.82);
      TEST_EQUALITY(ghostedVecE[1],-.91);
      TEST_EQUALITY(ghostedVecE[2],2.72);
      TEST_EQUALITY(ghostedVecE[3],6.23);
      TEST_EQUALITY(ghostedVecE[4],-.17);
    }

    TEST_EQUALITY(ghostedVecE[0],thyraVec[0]);
    TEST_EQUALITY(ghostedVecE[1],thyraVec[1]);
    TEST_EQUALITY(ghostedVecE[2],thyraVec[2]);
    TEST_EQUALITY(ghostedVecE[3],thyraVec[3]);
    TEST_EQUALITY(ghostedVecE[4],thyraVec[4]);
  }
}

}

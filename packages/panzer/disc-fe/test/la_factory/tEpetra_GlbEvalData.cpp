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

#include "PanzerDiscFE_config.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

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
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_ProductVectorBase.hpp"

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
    RCP<const Thyra_Vector>  uniqueVecT = ged.getUniqueVector();

    TEST_ASSERT(uniqueVecT!=Teuchos::null); 

    RCP<const Thyra::SpmdVectorSpaceBase<double> > uniqueSpace 
        = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(uniqueVecT->space());
    
    TEST_EQUALITY(uniqueSpace->isLocallyReplicated(),false);

    RCP<const Thyra_SpmdVec> spmdVec = rcp_dynamic_cast<const Thyra_SpmdVec>(uniqueVecT);

    Teuchos::ArrayRCP<const double> thyraVec;
    spmdVec->getLocalData(Teuchos::ptrFromRef(thyraVec));

    TEST_EQUALITY(thyraVec.size(),3);

    if(comm.MyPID()==0) {
      TEST_EQUALITY(thyraVec[0],3.14);
      TEST_EQUALITY(thyraVec[1],1.82);
      TEST_EQUALITY(thyraVec[2],-.91);
    }
    else {
      TEST_EQUALITY(thyraVec[0],2.72);
      TEST_EQUALITY(thyraVec[1],6.23);
      TEST_EQUALITY(thyraVec[2],-.17);
    }
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

TEUCHOS_UNIT_TEST(tEpetra_GlbEvalData, blocked)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Thyra::SpmdVectorBase<double> Thyra_SpmdVec;

  Teuchos::RCP<Teuchos::MpiComm<Thyra::Ordinal> > tComm = Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(MPI_COMM_WORLD));
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

  RCP<EpetraVector_ReadOnly_GlobalEvaluationData> ged_a, ged_b;
  ged_a = rcp(new EpetraVector_ReadOnly_GlobalEvaluationData(importer,ghostedMap,uniqueMap));
  ged_b = rcp(new EpetraVector_ReadOnly_GlobalEvaluationData(importer,ghostedMap,uniqueMap));
  std::vector<RCP<ReadOnlyVector_GlobalEvaluationData> > gedBlocks;
  gedBlocks.push_back(ged_a);
  gedBlocks.push_back(ged_b);

  RCP<const Thyra::VectorSpaceBase<double> > uniqueSpace_ab = Thyra::defaultSpmdVectorSpace<double>(tComm,3,-1);
  RCP<const Thyra::VectorSpaceBase<double> > ghostedSpace_ab = ged_a->getGhostedVector()->space();

  RCP<Thyra::DefaultProductVectorSpace<double> > uniqueSpace = Thyra::productVectorSpace<double>(uniqueSpace_ab,2);
  RCP<Thyra::DefaultProductVectorSpace<double> > ghostedSpace = Thyra::productVectorSpace<double>(ghostedSpace_ab,2);

  BlockedVector_ReadOnly_GlobalEvaluationData ged;
  ged.initialize(ghostedSpace,uniqueSpace,gedBlocks);

  RCP<Thyra::VectorBase<double> > uniqueVec = Thyra::createMember(*uniqueSpace);

  {
    RCP<Thyra::ProductVectorBase<double> > vec = Thyra::castOrCreateNonconstProductVectorBase(uniqueVec);

    TEST_ASSERT(vec->productSpace()->numBlocks()==2);

    if(comm.MyPID()==0) {
      Teuchos::ArrayRCP<double> thyraVec;

      rcp_dynamic_cast<Thyra_SpmdVec>(vec->getNonconstVectorBlock(0))->getNonconstLocalData(Teuchos::ptrFromRef(thyraVec));
      thyraVec[0] = 3.14;
      thyraVec[1] = 1.82;
      thyraVec[2] = -.91;

      rcp_dynamic_cast<Thyra_SpmdVec>(vec->getNonconstVectorBlock(1))->getNonconstLocalData(Teuchos::ptrFromRef(thyraVec));
      thyraVec[0] = 3.14+9.0;
      thyraVec[1] = 1.82+9.0;
      thyraVec[2] = -.91+9.0;
    }
    else {
      Teuchos::ArrayRCP<double> thyraVec;

      rcp_dynamic_cast<Thyra_SpmdVec>(vec->getNonconstVectorBlock(0))->getNonconstLocalData(Teuchos::ptrFromRef(thyraVec));
      thyraVec[0] = 2.72;
      thyraVec[1] = 6.23;
      thyraVec[2] = -.17;

      rcp_dynamic_cast<Thyra_SpmdVec>(vec->getNonconstVectorBlock(1))->getNonconstLocalData(Teuchos::ptrFromRef(thyraVec));
      thyraVec[0] = 2.72+7.0;
      thyraVec[1] = 6.23+7.0;
      thyraVec[2] = -.17+7.0;
    }
  }

  ged.setUniqueVector(uniqueVec);

  ged.initializeData();
  ged.globalToGhost(0);

  {
    RCP<Thyra::ProductVectorBase<double> > ghostedVec = Thyra::castOrCreateNonconstProductVectorBase(ged.getGhostedVector());

    if(comm.MyPID()==0) {
      Teuchos::ArrayRCP<const double> thyraVec;
      rcp_dynamic_cast<const Thyra_SpmdVec>(ghostedVec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(thyraVec));

      TEST_EQUALITY(thyraVec[0],3.14);
      TEST_EQUALITY(thyraVec[1],1.82);
      TEST_EQUALITY(thyraVec[2],-.91);
      TEST_EQUALITY(thyraVec[3],2.72);
      TEST_EQUALITY(thyraVec[4],6.23);

      rcp_dynamic_cast<const Thyra_SpmdVec>(ghostedVec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(thyraVec));

      TEST_EQUALITY(thyraVec[0],3.14+9.0);
      TEST_EQUALITY(thyraVec[1],1.82+9.0);
      TEST_EQUALITY(thyraVec[2],-.91+9.0);
      TEST_EQUALITY(thyraVec[3],2.72+7.0);
      TEST_EQUALITY(thyraVec[4],6.23+7.0);
    }
    else {
      Teuchos::ArrayRCP<const double> thyraVec;
      rcp_dynamic_cast<const Thyra_SpmdVec>(ghostedVec->getVectorBlock(0))->getLocalData(Teuchos::ptrFromRef(thyraVec));

      TEST_EQUALITY(thyraVec[0],1.82);
      TEST_EQUALITY(thyraVec[1],-.91);
      TEST_EQUALITY(thyraVec[2],2.72);
      TEST_EQUALITY(thyraVec[3],6.23);
      TEST_EQUALITY(thyraVec[4],-.17);

      rcp_dynamic_cast<const Thyra_SpmdVec>(ghostedVec->getVectorBlock(1))->getLocalData(Teuchos::ptrFromRef(thyraVec));

      TEST_EQUALITY(thyraVec[0],1.82+9.0);
      TEST_EQUALITY(thyraVec[1],-.91+9.0);
      TEST_EQUALITY(thyraVec[2],2.72+7.0);
      TEST_EQUALITY(thyraVec[3],6.23+7.0);
      TEST_EQUALITY(thyraVec[4],-.17+7.0);
    }
  }
}

TEUCHOS_UNIT_TEST(tEpetra_GlbEvalData, filtered_dofs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // This is required
  TEST_ASSERT(comm.NumProc()==2);

  std::vector<int> ghosted(4);
  std::vector<int> unique(2);

  // This is a line with 6 notes (numbered 0-5). The boundaries
  // are removed at 0 and 5.
  if(comm.MyPID()==0) {
    unique[0] = 1;
    unique[1] = 2;

    ghosted[0] = 0;
    ghosted[1] = 1;
    ghosted[2] = 2;
    ghosted[3] = 3;
  }
  else {
    unique[0] = 3;
    unique[1] = 4;

    ghosted[0] = 2;
    ghosted[1] = 3;
    ghosted[2] = 4;
    ghosted[3] = 5;
  }

  RCP<const Epetra_Map> uniqueMap = rcp(new Epetra_Map(-1,unique.size(),&unique[0],0,comm));
  RCP<const Epetra_Map> ghostedMap = rcp(new Epetra_Map(-1,ghosted.size(),&ghosted[0],0,comm));
  RCP<const Epetra_Import> importer = rcp(new Epetra_Import(*ghostedMap,*uniqueMap));

  EpetraVector_ReadOnly_GlobalEvaluationData ged;
 
  std::vector<int> constIndex(1);
 
  // setup filtered values
  constIndex[0] = 0;
  ged.useConstantValues(constIndex,2.0);

  constIndex[0] = 5;
  ged.useConstantValues(constIndex,3.0);

  ged.initialize(importer,ghostedMap,uniqueMap);

  TEST_THROW(ged.useConstantValues(constIndex,4.0),std::logic_error);

  {
    RCP<Epetra_Vector> uniqueVec = rcp(new Epetra_Vector(*uniqueMap));

    if(comm.MyPID()==0) {
      (*uniqueVec)[0] = 3.14;
      (*uniqueVec)[1] = 1.82;
    }
    else {
      (*uniqueVec)[0] = 2.72;
      (*uniqueVec)[1] = 6.23;
    }

    // set the unique vector, assure that const can be used
    ged.setUniqueVector_Epetra(uniqueVec.getConst());
  }

  // actually do something...
  ged.initializeData();
  ged.globalToGhost(0);

  // check values making sure that the constants are there
  {
    const Epetra_Vector & ghostedVecE = *ged.getGhostedVector_Epetra();

    if(comm.MyPID()==0) {
      TEST_EQUALITY(ghostedVecE[0],2.0);   // <= Replaced constant value
      TEST_EQUALITY(ghostedVecE[1],3.14);
      TEST_EQUALITY(ghostedVecE[2],1.82);
      TEST_EQUALITY(ghostedVecE[3],2.72);
    }
    else {
      TEST_EQUALITY(ghostedVecE[0],1.82);
      TEST_EQUALITY(ghostedVecE[1],2.72);
      TEST_EQUALITY(ghostedVecE[2],6.23);
      TEST_EQUALITY(ghostedVecE[3],3.0);   // <= Replaced constant value
    }
  }
}

} // end namespace panzer

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

#include "Phalanx_KokkosUtilities.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

TEUCHOS_UNIT_TEST(tTpetra_GlbEvalData, basic)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Thyra::VectorBase<double> Thyra_Vector;
  typedef Thyra::SpmdVectorBase<double> Thyra_SpmdVec;
  typedef Thyra::SpmdVectorSpaceBase<double> Thyra_SpmdVecSpace;

  typedef Tpetra::Vector<double,int,panzer::Ordinal64> Tpetra_Vector;
  typedef Tpetra::Map<int,panzer::Ordinal64> Tpetra_Map;
  typedef Tpetra::Import<int,panzer::Ordinal64> Tpetra_Import;


  Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

  // This is required
  TEST_ASSERT(comm->getSize()==2);

  std::vector<panzer::Ordinal64> ghosted(5);
  std::vector<panzer::Ordinal64> unique(3);

  if(comm->getRank()==0) {
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

  RCP<const Tpetra_Map> uniqueMap = rcp(new Tpetra_Map(-1,unique,0,comm));
  RCP<const Tpetra_Map> ghostedMap = rcp(new Tpetra_Map(-1,ghosted,0,comm));
  RCP<const Tpetra_Import> importer = rcp(new Tpetra_Import(uniqueMap,ghostedMap));

  TpetraVector_ReadOnly_GlobalEvaluationData<double,int,panzer::Ordinal64> ged;
 
  TEST_ASSERT(!ged.isInitialized());

  ged.initialize(importer,ghostedMap,uniqueMap);

  TEST_ASSERT(ged.isInitialized());

  // test the ghosted vector sizing (we don't care what the entries are!)
  { 
    RCP<Tpetra_Vector> ghostedVecTp = ged.getGhostedVector_Tpetra();
    RCP<Thyra_Vector>  ghostedVecT = ged.getGhostedVector();

    TEST_ASSERT(ghostedVecTp!=Teuchos::null); 
    TEST_ASSERT(ghostedVecT!=Teuchos::null); 

    RCP<const Thyra::SpmdVectorSpaceBase<double> > ghostedSpace 
        = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(ghostedVecT->space());
    
    TEST_EQUALITY(ghostedMap->getNodeNumElements(),ghostedVecTp->getLocalLength());
    TEST_EQUALITY(ghostedMap->getGlobalNumElements(),ghostedVecTp->getGlobalLength());

    TEST_EQUALITY(ghostedSpace->isLocallyReplicated(),false);
    TEST_EQUALITY(Teuchos::as<size_t>(ghostedSpace->localSubDim()),ghostedVecTp->getLocalLength());
  }

  // test setting a unique vector
  {
    RCP<Tpetra_Vector> uniqueVec_tp = rcp(new Tpetra_Vector(uniqueMap));
    auto uv_2d = uniqueVec_tp->getLocalView<Kokkos::HostSpace> ();
    auto uniqueVec = Kokkos::subview (uv_2d, Kokkos::ALL (), 0);

    if(comm->getRank()==0) {
      uniqueVec(0) = 3.14;
      uniqueVec(1) = 1.82;
      uniqueVec(2) = -.91;
    }
    else {
      uniqueVec(0) = 2.72;
      uniqueVec(1) = 6.23;
      uniqueVec(2) = -.17;
    }

    // set the unique vector, assure that const can be used
    ged.setUniqueVector_Tpetra(uniqueVec_tp.getConst());
  }

  // test the unique vector sizing and thyra entries
  { 
    const Tpetra_Vector & uniqueVecTp = *ged.getUniqueVector_Tpetra();
    auto uv_2d = uniqueVecTp.getLocalView<Kokkos::HostSpace> ();
    auto uniqueVecTpKv = Kokkos::subview (uv_2d, Kokkos::ALL (), 0);

    RCP<const Thyra_Vector>  uniqueVecT = ged.getUniqueVector();

    TEST_ASSERT(uniqueVecT!=Teuchos::null); 

    RCP<const Thyra::SpmdVectorSpaceBase<double> > uniqueSpace 
        = rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(uniqueVecT->space());
    
    TEST_EQUALITY(uniqueMap->getNodeNumElements(),uniqueVecTp.getLocalLength());
    TEST_EQUALITY(uniqueMap->getGlobalNumElements(),uniqueVecTp.getGlobalLength());

    TEST_EQUALITY(uniqueSpace->isLocallyReplicated(),false);
    TEST_EQUALITY(Teuchos::as<size_t>(uniqueSpace->localSubDim()),uniqueVecTp.getLocalLength());

    RCP<const Thyra_SpmdVec> spmdVec = rcp_dynamic_cast<const Thyra_SpmdVec>(uniqueVecT);

    Teuchos::ArrayRCP<const double> thyraVec;
    spmdVec->getLocalData(Teuchos::ptrFromRef(thyraVec));

    TEST_EQUALITY(Teuchos::as<size_t>(thyraVec.size()),uniqueVecTp.getLocalLength());

    if(comm->getRank()==0) {
      TEST_EQUALITY(uniqueVecTpKv(0),3.14);
      TEST_EQUALITY(uniqueVecTpKv(1),1.82);
      TEST_EQUALITY(uniqueVecTpKv(2),-.91);
    }
    else {
      TEST_EQUALITY(uniqueVecTpKv(0),2.72);
      TEST_EQUALITY(uniqueVecTpKv(1),6.23);
      TEST_EQUALITY(uniqueVecTpKv(2),-.17);
    }

    TEST_EQUALITY(uniqueVecTpKv(0),thyraVec[0]);
    TEST_EQUALITY(uniqueVecTpKv(1),thyraVec[1]);
    TEST_EQUALITY(uniqueVecTpKv(2),thyraVec[2]);
  }

  // actually do something...
  ged.initializeData();
  ged.globalToGhost(0);

  {
    const Tpetra_Vector & ghostedVecTp = *ged.getGhostedVector_Tpetra();
    auto uv_2d = ghostedVecTp.getLocalView<Kokkos::HostSpace> ();
    auto ghostedVecKv = Kokkos::subview (uv_2d, Kokkos::ALL (), 0);

    RCP<Thyra_Vector>  ghostedVecT = ged.getGhostedVector();
    RCP<const Thyra_SpmdVec> spmdVec = rcp_dynamic_cast<const Thyra_SpmdVec>(ghostedVecT);

    Teuchos::ArrayRCP<const double> thyraVec;
    spmdVec->getLocalData(Teuchos::ptrFromRef(thyraVec));

    TEST_EQUALITY(Teuchos::as<size_t>(thyraVec.size()),ghostedVecTp.getLocalLength());

    if(comm->getRank()==0) {
      TEST_EQUALITY(ghostedVecKv(0),3.14);
      TEST_EQUALITY(ghostedVecKv(1),1.82);
      TEST_EQUALITY(ghostedVecKv(2),-.91);
      TEST_EQUALITY(ghostedVecKv(3),2.72);
      TEST_EQUALITY(ghostedVecKv(4),6.23);
    }
    else {
      TEST_EQUALITY(ghostedVecKv(0),1.82);
      TEST_EQUALITY(ghostedVecKv(1),-.91);
      TEST_EQUALITY(ghostedVecKv(2),2.72);
      TEST_EQUALITY(ghostedVecKv(3),6.23);
      TEST_EQUALITY(ghostedVecKv(4),-.17);
    }

    TEST_EQUALITY(ghostedVecKv(0),thyraVec[0]);
    TEST_EQUALITY(ghostedVecKv(1),thyraVec[1]);
    TEST_EQUALITY(ghostedVecKv(2),thyraVec[2]);
    TEST_EQUALITY(ghostedVecKv(3),thyraVec[3]);
    TEST_EQUALITY(ghostedVecKv(4),thyraVec[4]);
  }
}

TEUCHOS_UNIT_TEST(tTpetra_GlbEvalData, blocked)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Thyra::VectorBase<double> Thyra_Vector;
  typedef Thyra::SpmdVectorBase<double> Thyra_SpmdVec;
  typedef Thyra::SpmdVectorSpaceBase<double> Thyra_SpmdVecSpace;

  typedef Tpetra::Vector<double,int,panzer::Ordinal64> Tpetra_Vector;
  typedef Tpetra::Map<int,panzer::Ordinal64> Tpetra_Map;
  typedef Tpetra::Import<int,panzer::Ordinal64> Tpetra_Import;


  Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

  // This is required
  TEST_ASSERT(comm->getSize()==2);

  std::vector<Ordinal64> ghosted(5);
  std::vector<Ordinal64> unique(3);

  if(comm->getRank()==0) {
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

  RCP<const Tpetra_Map> uniqueMap = rcp(new Tpetra_Map(-1,unique,0,comm));
  RCP<const Tpetra_Map> ghostedMap = rcp(new Tpetra_Map(-1,ghosted,0,comm));
  RCP<const Tpetra_Import> importer = rcp(new Tpetra_Import(uniqueMap,ghostedMap));

  RCP<TpetraVector_ReadOnly_GlobalEvaluationData<double,int,Ordinal64> > ged_a, ged_b;
  ged_a = rcp(new TpetraVector_ReadOnly_GlobalEvaluationData<double,int,Ordinal64>(importer,ghostedMap,uniqueMap));
  ged_b = rcp(new TpetraVector_ReadOnly_GlobalEvaluationData<double,int,Ordinal64>(importer,ghostedMap,uniqueMap));
  std::vector<RCP<ReadOnlyVector_GlobalEvaluationData> > gedBlocks;
  gedBlocks.push_back(ged_a);
  gedBlocks.push_back(ged_b);

  RCP<const Thyra::VectorSpaceBase<double> > uniqueSpace_ab = Thyra::tpetraVectorSpace<double,int,Ordinal64>(uniqueMap);
  RCP<const Thyra::VectorSpaceBase<double> > ghostedSpace_ab = ged_a->getGhostedVector()->space();

  RCP<Thyra::DefaultProductVectorSpace<double> > uniqueSpace = Thyra::productVectorSpace<double>(uniqueSpace_ab,2);
  RCP<Thyra::DefaultProductVectorSpace<double> > ghostedSpace = Thyra::productVectorSpace<double>(ghostedSpace_ab,2);

  BlockedVector_ReadOnly_GlobalEvaluationData ged;
  ged.initialize(ghostedSpace,uniqueSpace,gedBlocks);

  RCP<Thyra::VectorBase<double> > uniqueVec = Thyra::createMember(*uniqueSpace);

  {
    RCP<Thyra::ProductVectorBase<double> > vec = Thyra::castOrCreateNonconstProductVectorBase(uniqueVec);

    TEST_ASSERT(vec->productSpace()->numBlocks()==2);

    if(comm->getRank()==0) {
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

    if(comm->getRank()==0) {
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

TEUCHOS_UNIT_TEST(tTpetra_GlbEvalData, filtered_dofs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Thyra::VectorBase<double> Thyra_Vector;
  typedef Thyra::SpmdVectorBase<double> Thyra_SpmdVec;
  typedef Thyra::SpmdVectorSpaceBase<double> Thyra_SpmdVecSpace;

  typedef Tpetra::Vector<double,int,panzer::Ordinal64> Tpetra_Vector;
  typedef Tpetra::Map<int,panzer::Ordinal64> Tpetra_Map;
  typedef Tpetra::Import<int,panzer::Ordinal64> Tpetra_Import;


  Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

  // This is required
  TEST_ASSERT(comm->getSize()==2);

  std::vector<panzer::Ordinal64> ghosted(4);
  std::vector<panzer::Ordinal64> unique(2);

  // This is a line with 6 notes (numbered 0-5). The boundaries
  // are removed at 0 and 5.
  if(comm->getRank()==0) {
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

  RCP<const Tpetra_Map> uniqueMap = rcp(new Tpetra_Map(-1,unique,0,comm));
  RCP<const Tpetra_Map> ghostedMap = rcp(new Tpetra_Map(-1,ghosted,0,comm));
  RCP<const Tpetra_Import> importer = rcp(new Tpetra_Import(uniqueMap,ghostedMap));

  TpetraVector_ReadOnly_GlobalEvaluationData<double,int,panzer::Ordinal64> ged;
 
  std::vector<panzer::Ordinal64> constIndex(1);
 
  // setup filtered values
  constIndex[0] = 0;
  ged.useConstantValues(constIndex,2.0);

  constIndex[0] = 5;
  ged.useConstantValues(constIndex,3.0);

  ged.initialize(importer,ghostedMap,uniqueMap);

  TEST_THROW(ged.useConstantValues(constIndex,4.0),std::logic_error);

  {
    RCP<Tpetra_Vector> uniqueVecTp = rcp(new Tpetra_Vector(uniqueMap));
    auto uv_2d = uniqueVecTp->getLocalView<Kokkos::HostSpace> ();
    auto uniqueVec = Kokkos::subview (uv_2d, Kokkos::ALL (), 0);

    if(comm->getRank()==0) {
      uniqueVec(0) = 3.14;
      uniqueVec(1) = 1.82;
    }
    else {
      uniqueVec(0) = 2.72;
      uniqueVec(1) = 6.23;
    }

    // set the unique vector, assure that const can be used
    ged.setUniqueVector_Tpetra(uniqueVecTp.getConst());
  }

  // actually do something...
  ged.initializeData();
  ged.globalToGhost(0);

  // check values making sure that the constants are there
  {
    const Tpetra_Vector & ghostedVecTp = *ged.getGhostedVector_Tpetra();
    auto uv_2d = ghostedVecTp.getLocalView<Kokkos::HostSpace> ();
    auto ghostedVecKv = Kokkos::subview (uv_2d, Kokkos::ALL (), 0);

    if(comm->getRank()==0) {
      TEST_EQUALITY(ghostedVecKv(0),2.0);   // <= Replaced constant value
      TEST_EQUALITY(ghostedVecKv(1),3.14);
      TEST_EQUALITY(ghostedVecKv(2),1.82);
      TEST_EQUALITY(ghostedVecKv(3),2.72);
    }
    else {
      TEST_EQUALITY(ghostedVecKv(0),1.82);
      TEST_EQUALITY(ghostedVecKv(1),2.72);
      TEST_EQUALITY(ghostedVecKv(2),6.23);
      TEST_EQUALITY(ghostedVecKv(3),3.0);   // <= Replaced constant value
    }
  }
}

} // end namespace panzer

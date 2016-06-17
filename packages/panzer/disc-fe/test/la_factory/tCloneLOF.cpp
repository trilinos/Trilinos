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
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_LinearObjFactory_Utilities.hpp"

#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"

#include "UnitTest_ConnManager.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

namespace panzer {

template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<double,FieldContainer> > basis = rcp(new Intrepid2Type);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

TEUCHOS_UNIT_TEST(tCloneLOF, epetra)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager<int>(myRank,numProc));

   RCP<const FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   RCP<panzer::DOFManager<int,int> > indexer = rcp(new panzer::DOFManager<int,int>());
   indexer->setConnManager(connManager,MPI_COMM_WORLD);
   indexer->addField("U",patternC1);
   indexer->addField("V",patternC1);
   indexer->buildGlobalUnknowns();

   RCP<panzer::DOFManager<int,int> > control_indexer = rcp(new panzer::DOFManager<int,int>());
   control_indexer->setConnManager(connManager,MPI_COMM_WORLD);
   control_indexer->addField("Z",patternC1);
   control_indexer->buildGlobalUnknowns();
 
   // setup factory
   RCP<EpetraLinearObjFactory<Traits,int> > ep_lof
         = Teuchos::rcp(new EpetraLinearObjFactory<Traits,int>(tComm.getConst(),indexer));
   
   // this is the member we are testing!
   RCP<const LinearObjFactory<Traits> > control_lof = cloneWithNewDomain(*ep_lof,control_indexer);
   RCP<const EpetraLinearObjFactory<Traits,int> > ep_control_lof 
       = rcp_dynamic_cast<const EpetraLinearObjFactory<Traits,int> >(control_lof);

   std::vector<int> control_owned;
   control_indexer->getOwnedIndices(control_owned);

   TEST_ASSERT(ep_control_lof->getMap(0)->SameAs(*ep_lof->getMap(0)));
   TEST_EQUALITY(ep_control_lof->getColMap(0)->NumMyElements(),Teuchos::as<int>(control_owned.size()));
}

TEUCHOS_UNIT_TEST(tCloneLOF, blocked_epetra)
{
   typedef Thyra::ProductVectorBase<double> PVector;
   typedef Thyra::BlockedLinearOpBase<double> BLinearOp;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));


   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager<int>(myRank,numProc));

   RCP<const FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();
   RCP<const FieldPattern> patternC2
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<double,FieldContainer> >();

   RCP<panzer::BlockedDOFManager<int,int> > indexer = rcp(new panzer::BlockedDOFManager<int,int>());
   {
     std::vector<std::vector<std::string> > fieldOrder(2);

     indexer->setConnManager(connManager,MPI_COMM_WORLD);
     indexer->addField("U",patternC1);
     indexer->addField("V",patternC1);

     fieldOrder[0].push_back("U");
     fieldOrder[1].push_back("V");
     indexer->setFieldOrder(fieldOrder);
     indexer->buildGlobalUnknowns();
   }

   RCP<panzer::BlockedDOFManager<int,int> > control_indexer = rcp(new panzer::BlockedDOFManager<int,int>());
   {
     std::vector<std::vector<std::string> > fieldOrder(1);
     control_indexer->setConnManager(connManager,MPI_COMM_WORLD);
     control_indexer->addField("Z",patternC1);
     fieldOrder[0].push_back("Z");
     indexer->setFieldOrder(fieldOrder);
     control_indexer->buildGlobalUnknowns();

     patternC2->print(out);

     {
     std::vector<int> gids;
     control_indexer->getFieldDOFManagers()[0]->getOwnedIndices(gids);
     out << "GIDs 0 = " << gids.size() << std::endl;
     }
   }

   // setup factory
   RCP<BlockedEpetraLinearObjFactory<Traits,int> > ep_lof
         = Teuchos::rcp(new BlockedEpetraLinearObjFactory<Traits,int>(tComm,indexer));
   
   // this is the member we are testing!
   RCP<const LinearObjFactory<Traits> > control_lof = cloneWithNewDomain(*ep_lof,control_indexer);
   RCP<const BlockedEpetraLinearObjFactory<Traits,int> > ep_control_lof 
       = rcp_dynamic_cast<const BlockedEpetraLinearObjFactory<Traits,int> >(control_lof);

   RCP<BLinearOp> mat = rcp_dynamic_cast<BLinearOp>(ep_control_lof->getThyraMatrix()); 
   RCP<BLinearOp> gmat = rcp_dynamic_cast<BLinearOp>(ep_control_lof->getGhostedThyraMatrix()); 
   RCP<PVector>   x   = rcp_dynamic_cast<PVector>(ep_control_lof->getThyraDomainVector()); 
   RCP<PVector>   gx   = rcp_dynamic_cast<PVector>(ep_control_lof->getGhostedThyraDomainVector()); 
   RCP<PVector>   f   = rcp_dynamic_cast<PVector>(ep_control_lof->getThyraRangeVector()); 
   RCP<PVector>   gf   = rcp_dynamic_cast<PVector>(ep_control_lof->getGhostedThyraRangeVector()); 

   TEST_EQUALITY(x->productSpace()->numBlocks(),1);
   TEST_EQUALITY(x->productSpace()->dim(),18);
   TEST_EQUALITY(gx->productSpace()->numBlocks(),1);
   TEST_EQUALITY(gx->productSpace()->dim(),10+15);

   TEST_EQUALITY(f->productSpace()->numBlocks(),2);
   TEST_EQUALITY(f->productSpace()->dim(),36);
   TEST_EQUALITY(gf->productSpace()->numBlocks(),2);
   TEST_EQUALITY(gf->productSpace()->dim(),50);


   TEST_EQUALITY(mat->productRange()->numBlocks(),2);
   TEST_EQUALITY(mat->productRange()->dim(),36);
   TEST_EQUALITY(mat->productDomain()->numBlocks(),1);
   TEST_EQUALITY(mat->productDomain()->dim(),18);

   TEST_EQUALITY(gmat->productRange()->numBlocks(),2);
   TEST_EQUALITY(gmat->productRange()->dim(),50);
   TEST_EQUALITY(gmat->productDomain()->numBlocks(),1);
   TEST_EQUALITY(gmat->productDomain()->dim(),10+15);
}

TEUCHOS_UNIT_TEST(tCloneLOF, blocked_epetra_nonblocked_domain)
{
   typedef Thyra::ProductVectorBase<double> PVector;
   typedef Thyra::BlockedLinearOpBase<double> BLinearOp;
   typedef Thyra::VectorBase<double> Vector;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));


   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager<int>(myRank,numProc));

   RCP<const FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();
   RCP<const FieldPattern> patternC2
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<double,FieldContainer> >();

   RCP<panzer::BlockedDOFManager<int,int> > indexer = rcp(new panzer::BlockedDOFManager<int,int>());
   {
     std::vector<std::vector<std::string> > fieldOrder(2);

     indexer->setConnManager(connManager,MPI_COMM_WORLD);
     indexer->addField("U",patternC1);
     indexer->addField("V",patternC1);

     fieldOrder[0].push_back("U");
     fieldOrder[1].push_back("V");
     indexer->setFieldOrder(fieldOrder);
     indexer->buildGlobalUnknowns();
   }

   RCP<panzer::DOFManager<int,int> > control_indexer = rcp(new panzer::DOFManager<int,int>());
   {
     control_indexer->setConnManager(connManager,MPI_COMM_WORLD);
     control_indexer->addField("Z",patternC1);
     control_indexer->buildGlobalUnknowns();

     patternC2->print(out);
   }

   // setup factory
   out << "build lof" << std::endl;
   RCP<BlockedEpetraLinearObjFactory<Traits,int> > ep_lof
         = Teuchos::rcp(new BlockedEpetraLinearObjFactory<Traits,int>(tComm,indexer));
   
   // this is the member we are testing!
   out << "cloning lof" << std::endl;
   RCP<const LinearObjFactory<Traits> > control_lof = cloneWithNewDomain(*ep_lof,control_indexer);

   out << "casting lof" << std::endl;
   RCP<const BlockedEpetraLinearObjFactory<Traits,int> > ep_control_lof 
       = rcp_dynamic_cast<const BlockedEpetraLinearObjFactory<Traits,int> >(control_lof,true);

   out << "using casted lof" << std::endl;
   RCP<BLinearOp> mat  = rcp_dynamic_cast<BLinearOp>(ep_control_lof->getThyraMatrix(),true); 
   RCP<BLinearOp> gmat = rcp_dynamic_cast<BLinearOp>(ep_control_lof->getGhostedThyraMatrix(),true); 
   RCP<Vector>    x    = ep_control_lof->getThyraDomainVector();
   RCP<Vector>    gx   = ep_control_lof->getGhostedThyraDomainVector();
   RCP<PVector>   f    = rcp_dynamic_cast<PVector>(ep_control_lof->getThyraRangeVector(),true); 
   RCP<PVector>   gf   = rcp_dynamic_cast<PVector>(ep_control_lof->getGhostedThyraRangeVector(),true); 

   TEST_EQUALITY(x->space()->dim(),18);
   TEST_EQUALITY(gx->space()->dim(),10+15);

   TEST_EQUALITY(f->productSpace()->numBlocks(),2);
   TEST_EQUALITY(f->productSpace()->dim(),36);
   TEST_EQUALITY(gf->productSpace()->numBlocks(),2);
   TEST_EQUALITY(gf->productSpace()->dim(),50);


   TEST_EQUALITY(mat->productRange()->numBlocks(),2);
   TEST_EQUALITY(mat->productRange()->dim(),36);
   TEST_EQUALITY(mat->productDomain()->numBlocks(),1);
   TEST_EQUALITY(mat->productDomain()->dim(),18);

   TEST_EQUALITY(gmat->productRange()->numBlocks(),2);
   TEST_EQUALITY(gmat->productRange()->dim(),50);
   TEST_EQUALITY(gmat->productDomain()->numBlocks(),1);
   TEST_EQUALITY(gmat->productDomain()->dim(),10+15);
}

}

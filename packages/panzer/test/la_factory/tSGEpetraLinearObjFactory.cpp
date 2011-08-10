#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_SGEpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

// for testing gather/scatter construction
#include "Panzer_Basis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_TensorProductQuadrature.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#include "Thyra_EpetraThyraWrappers.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > buildExpansion(int numDim,int order)
{
   Teuchos::Array<Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(numDim);
   for(int i=0;i<numDim;i++)
      bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(order));
   Teuchos::RCP<Stokhos::ProductBasis<int,double> > basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

   // build Cijk and "expansion"
   int kExpOrder = basis->size();
   // if(!fullExpansion)
   //    kExpOrder = numDim+1;
   Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor(kExpOrder);
   Teuchos::RCP<Stokhos::Quadrature<int,double> > quadrature = Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

   return Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis,Cijk,quadrature));
}

TEUCHOS_UNIT_TEST(tSGEpetraLinearObjFactory, basic)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   // panzer::pauseToAttach();

   RCP<Stokhos::OrthogPolyExpansion<int,double> > sgExpansion = buildExpansion(3,5);
   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new unit_test::UniqueGlobalIndexer(myRank,numProc));
 
   // setup factory
   RCP<panzer::EpetraLinearObjFactory<panzer::Traits,short> > epetraFactory
         = rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));
   RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,short> > la_factory
         = rcp(new panzer::SGEpetraLinearObjFactory<panzer::Traits,short>(epetraFactory,sgExpansion));
   
   RCP<panzer::LinearObjContainer> ghostedContainer = la_factory->buildGhostedLinearObjContainer(); 
   RCP<panzer::LinearObjContainer> container = la_factory->buildLinearObjContainer(); 

   TEST_NOTHROW(rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(ghostedContainer,true));
   TEST_NOTHROW(rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(container,true));

   RCP<panzer::SGEpetraLinearObjContainer> sgGhostedContainer 
         = rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(ghostedContainer,true);
   RCP<panzer::SGEpetraLinearObjContainer> sgContainer 
         = rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(container,true);

   // make sure we have correct number of sub containers
   TEST_EQUALITY(sgExpansion->size(),sgGhostedContainer->end()-sgGhostedContainer->begin());
   TEST_EQUALITY(sgExpansion->size(),sgContainer->end()-sgContainer->begin());

   // make sure all "sub-containers" are epetra containers
   panzer::SGEpetraLinearObjContainer::const_iterator itr;
   for(itr=sgGhostedContainer->begin();itr!=sgGhostedContainer->end();++itr)
      TEST_NOTHROW(rcp_dynamic_cast<EpetraLinearObjContainer>(*itr));
   for(itr=sgContainer->begin();itr!=sgContainer->end();++itr)
      TEST_NOTHROW(rcp_dynamic_cast<EpetraLinearObjContainer>(*itr));

   la_factory->ghostToGlobalContainer(*ghostedContainer,*container);
   la_factory->globalToGhostContainer(*container,*ghostedContainer);
}

TEUCHOS_UNIT_TEST(tSGEpetraLinearObjFactory, initializeContainer)
{
   typedef EpetraLinearObjContainer ELOC;

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   // panzer::pauseToAttach();

   RCP<Stokhos::OrthogPolyExpansion<int,double> > sgExpansion = buildExpansion(3,5);
   RCP<panzer::UniqueGlobalIndexer<short,int> > indexer 
         = rcp(new unit_test::UniqueGlobalIndexer(myRank,numProc));

   std::vector<int> ownedIndices, ownedAndSharedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndSharedIndices(ownedAndSharedIndices);
 
   // setup factory
   RCP<panzer::EpetraLinearObjFactory<panzer::Traits,short> > epetraFactory
         = rcp(new panzer::EpetraLinearObjFactory<panzer::Traits,short>(eComm.getConst(),indexer));
   RCP<panzer::SGEpetraLinearObjFactory<panzer::Traits,short> > la_factory
         = rcp(new panzer::SGEpetraLinearObjFactory<panzer::Traits,short>(epetraFactory,sgExpansion));
   
   RCP<panzer::LinearObjContainer> ghostedContainer = la_factory->buildGhostedLinearObjContainer(); 
   RCP<panzer::LinearObjContainer> container = la_factory->buildLinearObjContainer(); 

   RCP<panzer::SGEpetraLinearObjContainer> sgGhostedContainer 
         = rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(ghostedContainer,true);
   RCP<panzer::SGEpetraLinearObjContainer> sgContainer 
         = rcp_dynamic_cast<panzer::SGEpetraLinearObjContainer>(container,true);

   la_factory->initializeContainer(ELOC::X | ELOC::DxDt,*container);
   la_factory->initializeGhostedContainer(ELOC::X | ELOC::DxDt,*ghostedContainer);

   // make sure all "sub-containers" are epetra containers
   panzer::SGEpetraLinearObjContainer::const_iterator itr;
   for(itr=sgContainer->begin();itr!=sgContainer->end();++itr) {
      TEST_ASSERT((*itr)->x!=Teuchos::null);
      TEST_ASSERT((*itr)->dxdt!=Teuchos::null);
      TEST_EQUALITY((*itr)->f,Teuchos::null);
      TEST_EQUALITY((*itr)->A,Teuchos::null);
      TEST_EQUALITY((*itr)->x->MyLength(),(int) ownedIndices.size());
      TEST_EQUALITY((*itr)->dxdt->MyLength(),(int) ownedIndices.size());
   }
   for(itr=sgGhostedContainer->begin();itr!=sgGhostedContainer->end();++itr) {
      TEST_ASSERT((*itr)->x!=Teuchos::null);
      TEST_ASSERT((*itr)->dxdt!=Teuchos::null);
      TEST_EQUALITY((*itr)->f,Teuchos::null);
      TEST_EQUALITY((*itr)->A,Teuchos::null);
      TEST_EQUALITY((*itr)->x->MyLength(),(int) ownedAndSharedIndices.size());
      TEST_EQUALITY((*itr)->dxdt->MyLength(),(int) ownedAndSharedIndices.size());
   }

   la_factory->initializeContainer(ELOC::Mat | ELOC::F,*container);
   la_factory->initializeGhostedContainer(ELOC::Mat | ELOC::F,*ghostedContainer);

   // make sure all "sub-containers" are epetra containers
   for(itr=sgContainer->begin();itr!=sgContainer->end();++itr) {
      TEST_ASSERT((*itr)->f!=Teuchos::null);
      TEST_ASSERT((*itr)->A!=Teuchos::null);
      TEST_EQUALITY((*itr)->x,Teuchos::null);
      TEST_EQUALITY((*itr)->dxdt,Teuchos::null);
      TEST_EQUALITY((*itr)->f->MyLength(),(int) ownedIndices.size());
      TEST_EQUALITY((*itr)->A->NumMyRows(),(int) ownedIndices.size());
   }
   for(itr=sgGhostedContainer->begin();itr!=sgGhostedContainer->end();++itr) {
      TEST_ASSERT((*itr)->f!=Teuchos::null);
      TEST_ASSERT((*itr)->A!=Teuchos::null);
      TEST_EQUALITY((*itr)->x,Teuchos::null);
      TEST_EQUALITY((*itr)->dxdt,Teuchos::null);
      TEST_EQUALITY((*itr)->f->MyLength(),(int) ownedAndSharedIndices.size());
      TEST_EQUALITY((*itr)->A->NumMyRows(),(int) ownedAndSharedIndices.size());
   }

}

}

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
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

// for testing gather/scatter construction
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_ConnManager.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include "Tpetra_Operator.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_FECrsMatrix.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef double ScalarT;
typedef int LocalOrdinalT;
typedef panzer::GlobalOrdinal GlobalOrdinalT;
typedef panzer::TpetraNodeType NodeT;

typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
typedef Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> OperatorType;
typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

typedef Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> ThyraLinearOp;

typedef panzer::BlockedTpetraLinearObjContainer<double,int,panzer::GlobalOrdinal> BLOC;
typedef panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> BLOFact;


namespace panzer {

Teuchos::RCP<CrsMatrixType> getSubBlock_tp(int i,int j,Thyra::LinearOpBase<double> & lo)
{
   Thyra::BlockedLinearOpBase<double> & blo = Teuchos::dyn_cast<Thyra::BlockedLinearOpBase<double> >(lo);
   Teuchos::RCP<OperatorType> t_blo = 
       rcp_dynamic_cast<Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(blo.getNonconstBlock(i,j),true)->getTpetraOperator();

   return rcp_dynamic_cast<CrsMatrixType>(t_blo);
}

Teuchos::RCP<const CrsMatrixType> getSubBlock_tp(int i,int j,const Thyra::LinearOpBase<double> & lo)
{
   const Thyra::BlockedLinearOpBase<double> & blo = Teuchos::dyn_cast<const Thyra::BlockedLinearOpBase<double> >(lo);
   Teuchos::RCP<const OperatorType> t_blo = 
       rcp_dynamic_cast<const Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(blo.getBlock(i,j),true)->getConstTpetraOperator();

   return rcp_dynamic_cast<const CrsMatrixType>(t_blo);
}

void resumeFill(CrsMatrixType & A)
{
 A.resumeFill();
}

void fillComplete(CrsMatrixType & A)
{
  A.fillComplete(A.getDomainMap(),A.getRangeMap());
}

void putScalar(ScalarT s,CrsMatrixType & A)
{
  A.setAllToScalar(s);
}

template <typename Intrepid2Type>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // build a geometric pattern from a single basis
  RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
  RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
  return pattern;
}

Teuchos::RCP<const panzer::BlockedDOFManager> buildBlockedIndexer64(int myRank,int numProc,int numBlocks)
{
  std::string names[] = {"U","V","W","X"};

  Teuchos::RCP<const FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
  Teuchos::RCP<ConnManager> connManager = rcp(new unit_test::ConnManager(myRank,numProc));
  Teuchos::RCP<panzer::BlockedDOFManager> indexer = rcp(new panzer::BlockedDOFManager());

  indexer->setConnManager(connManager,MPI_COMM_WORLD);

  std::vector<std::vector<std::string> > fieldOrder(numBlocks);
  for(int i=0;i<numBlocks;i++) {
    indexer->addField(names[i],patternC1);

    fieldOrder[i].push_back(names[i]);
  }
  indexer->setFieldOrder(fieldOrder);
  indexer->buildGlobalUnknowns();

  return indexer;
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, intializeContainer_tpetra)
{

   panzer::BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> container;

   TEST_ASSERT(container.checkCompatibility());
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, tpetra_factory_tests)
{

   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // pauseToAttach();

   typedef LinearObjContainer LOC;

   int numBlocks = 3;
   int myRank = comm->getRank();
   int numProc = comm->getSize();

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(myRank,numProc,numBlocks);
   BLOFact factory(comm,blkIndexer);

   RCP<LinearObjContainer> container = factory.buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = factory.buildGhostedLinearObjContainer();
   TEST_ASSERT(container!=Teuchos::null);
   TEST_ASSERT(ghosted!=Teuchos::null);


   RCP<BLOC> bContainer = rcp_dynamic_cast<BLOC>(container);
   RCP<BLOC> b_ghosted = rcp_dynamic_cast<BLOC>(ghosted);
   TEST_ASSERT(bContainer!=Teuchos::null);
   TEST_ASSERT(b_ghosted!=Teuchos::null);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      factory.initializeContainer(LOC::X,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::DxDt,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::F,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      factory.initializeContainer(LOC::F | LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      factory.initializeContainer(LOC::X | LOC::DxDt,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      // everything
      factory.initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // Epetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      factory.initializeContainer(LOC::X,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::DxDt,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::F,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::Mat,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      factory.initializeContainer(LOC::F | LOC::Mat,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      factory.initializeContainer(LOC::X | LOC::DxDt,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      // everything
      factory.initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   }
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, ghostToGlobal)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = comm->getRank();
   int numProc = comm->getSize();

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(myRank,numProc,numBlocks);
   Teuchos::RCP<BLOFact> la_factory
         = Teuchos::rcp(new BLOFact(comm,blkIndexer));

   Teuchos::RCP<LinearObjContainer> global  = la_factory->buildLinearObjContainer();
   Teuchos::RCP<LinearObjContainer> ghosted = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeContainer(LinearObjContainer::Mat,*global);
   la_factory->initializeGhostedContainer(LinearObjContainer::Mat,*ghosted);

   Teuchos::rcp_dynamic_cast<BLOC>(ghosted)->initializeMatrix(1.0);

   la_factory->ghostToGlobalContainer(*ghosted,*global,LinearObjContainer::Mat);

   Teuchos::RCP<Thyra::LinearOpBase<double> > th_A = Teuchos::rcp_dynamic_cast<BLOC>(global)->get_A();
   Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blk_A = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(th_A);
  
   Teuchos::RCP<OperatorType> cA_00 = Teuchos::rcp_dynamic_cast<ThyraLinearOp>(blk_A->getNonconstBlock(0,0),true)->getTpetraOperator();
   Teuchos::RCP<OperatorType> cA_01 = Teuchos::rcp_dynamic_cast<ThyraLinearOp>(blk_A->getNonconstBlock(0,1),true)->getTpetraOperator();
   Teuchos::RCP<OperatorType> cA_10 = Teuchos::rcp_dynamic_cast<ThyraLinearOp>(blk_A->getNonconstBlock(1,0),true)->getTpetraOperator();
   Teuchos::RCP<OperatorType> cA_11 = Teuchos::rcp_dynamic_cast<ThyraLinearOp>(blk_A->getNonconstBlock(1,1),true)->getTpetraOperator();

   Teuchos::RCP<CrsMatrixType> A_00 = Teuchos::rcp_dynamic_cast<CrsMatrixType>(cA_00);
   Teuchos::RCP<CrsMatrixType> A_01 = Teuchos::rcp_dynamic_cast<CrsMatrixType>(cA_01);

   Teuchos::RCP<CrsMatrixType> A_10 = Teuchos::rcp_dynamic_cast<CrsMatrixType>(cA_10);
   Teuchos::RCP<CrsMatrixType> A_11 = Teuchos::rcp_dynamic_cast<CrsMatrixType>(cA_11);

   out << "A_00 = \n";
   A_00->print(out);

   out << "A_01 = \n";
   A_01->print(out);

   out << "A_10 = \n";
   A_10->print(out);

   out << "A_11 = \n";
   A_11->print(out);
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, graph_constr)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = comm->getRank();
   int numProc = comm->getSize();

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(myRank,numProc,numBlocks);
   Teuchos::RCP<BLOFact> la_factory = Teuchos::rcp(new BLOFact(comm,blkIndexer));

   Teuchos::RCP<CrsMatrixType> A_00 = la_factory->getGhostedTpetraMatrix(0,0); putScalar(1.0,*A_00);
   Teuchos::RCP<CrsMatrixType> A_01 = la_factory->getGhostedTpetraMatrix(0,1); putScalar(1.0,*A_01);

   Teuchos::RCP<CrsMatrixType> A_10 = la_factory->getGhostedTpetraMatrix(1,0); putScalar(1.0,*A_10);
   Teuchos::RCP<CrsMatrixType> A_11 = la_factory->getGhostedTpetraMatrix(1,1); putScalar(1.0,*A_11);

   out << "A_00 = \n";
   A_00->print(out);

   out << "A_01 = \n";
   A_01->print(out);

   out << "A_10 = \n";
   A_10->print(out);

   out << "A_11 = \n";
   A_11->print(out);
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, adjustDirichlet)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 3;
   int myRank = comm->getRank();
   int numProc = comm->getSize();
 
   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(myRank,numProc,numBlocks);
   Teuchos::RCP<BLOFact> la_factory
         = Teuchos::rcp(new BLOFact(comm,blkIndexer));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<BLOC> b_0   = rcp_dynamic_cast<BLOC>(ghosted_0);
   RCP<BLOC> b_1   = rcp_dynamic_cast<BLOC>(ghosted_1);
   RCP<BLOC> b_sys = rcp_dynamic_cast<BLOC>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(b_0->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_1->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_A()));

   Thyra::assign(b_0->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_1->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_sys->get_f().ptr(),-3.0); // put some garbage in the systems

   // b_sys->get_A()->PutScalar(-3.0);
   for(int i=0;i<numBlocks;i++) {
      for(int j=0;j<numBlocks;j++) {
         RCP<CrsMatrixType> M = getSubBlock_tp(i,j,*b_sys->get_A());
         M->setAllToScalar(-3.0);
      }
   }

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      for(int i=0;i<numBlocks;i++) {
         RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(i);
         RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(i);

         Teuchos::ArrayRCP<double> data_0,data_1;
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         // case 0
         data_0[0] = 1.0; // GID = 0
         data_1[0] = 1.0;
   
         // case 1
         data_0[2] = 1.0; // GID = 2
         data_1[2] = 2.0;
   
         // case 2
         data_1[5] = 2.0; // GID = 5
      }
   }
   else if(myRank==1) {
      for(int i=0;i<numBlocks;i++) {
         RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(i);
         RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(i);

         Teuchos::ArrayRCP<double> data_0,data_1;
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         // case 0
         data_0[3] = 1.0; // GID = 9
         data_1[3] = 1.0;
   
         // case 1
         data_0[0] = 1.0; // GID =2
         data_1[0] = 2.0;
   
         // case 2
         data_1[6] = 2.0; // GID = 4
      }
   }
   else 
      TEUCHOS_ASSERT(false);

   out << "LOCAL " << std::endl;
   b_0->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;
   out << "GLOBAL " << std::endl;
   b_1->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   std::size_t sz = 20;
   std::size_t numEntries = 0;
   typename CrsMatrixType::nonconst_local_inds_host_view_type indices("indices", sz);
   typename CrsMatrixType::nonconst_values_host_view_type values("values", sz);

   if(myRank==0) {   
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      for(int i=0;i<numBlocks;i++) {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);     // case 0
         TEST_EQUALITY(data[2],-3.0/2.0); // case 1
         TEST_EQUALITY(data[5],0.0);      // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const CrsMatrixType> subA = getSubBlock_tp(i,j,*A);

            subA->getLocalRowCopy(0,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);

            subA->getLocalRowCopy(2,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);

            subA->getLocalRowCopy(5,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else if(myRank==1) {
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      for(int i=0;i<numBlocks;i++) {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[3],-3.0);     // case 0
         TEST_EQUALITY(data[0],-3.0/2.0); // case 1
         TEST_EQUALITY(data[6],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const CrsMatrixType> subA = getSubBlock_tp(i,j,*A);

            subA->getLocalRowCopy(3,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
   
            subA->getLocalRowCopy(0,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);
   
            subA->getLocalRowCopy(6,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, node_cell)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = comm->getRank();
   int numProc = comm->getSize();
 
   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(myRank,numProc,numBlocks);
   Teuchos::RCP<BLOFact> la_factory = Teuchos::rcp(new BLOFact(comm,blkIndexer));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::F,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<BLOC> b_0   = rcp_dynamic_cast<BLOC>(ghosted_0);
   RCP<BLOC> b_1   = rcp_dynamic_cast<BLOC>(ghosted_1);
   RCP<BLOC> b_sys = rcp_dynamic_cast<BLOC>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(b_0->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_1->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_A()));

   Thyra::assign(b_0->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_1->get_f().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_sys->get_f().ptr(),-3.0); // put some garbage in the systems

   for(int i=0;i<numBlocks;i++)
      for(int j=0;j<numBlocks;j++)
         putScalar(-3.0,*getSubBlock_tp(i,j,*b_sys->get_A()));

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(0);
      RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(0);

      Teuchos::ArrayRCP<double> data_0,data_1;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

      // case 0
      data_0[0] = 1.0; // GID = 0
      data_1[0] = 1.0;

      // case 1
      data_0[2] = 1.0; // GID = 2
      data_1[2] = 2.0;

      // case 2
      data_1[5] = 2.0; // GID = 5

      {
         x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(1);
         x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(1);

         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         data_1[0] = 2.0;
      }
   }
   else if(myRank==1) {
      RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_f())->getNonconstVectorBlock(0);
      RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_f())->getNonconstVectorBlock(0);

      Teuchos::ArrayRCP<double> data_0,data_1;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

      // case 0
      data_0[3] = 1.0; // GID = 9
      data_1[3] = 1.0;

      // case 1
      data_0[0] = 1.0; // GID =2
      data_1[0] = 2.0;

      // case 2
      data_1[6] = 2.0; // GID = 4
   }
   else 
      TEUCHOS_ASSERT(false);

   out << "LOCAL " << std::endl;
   b_0->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;
   out << "GLOBAL " << std::endl;
   b_1->get_f()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   std::size_t sz = 20;
   std::size_t numEntries = 0;
   typename CrsMatrixType::nonconst_local_inds_host_view_type indices("indices", sz);
   typename CrsMatrixType::nonconst_values_host_view_type values("values", sz);

   if(myRank==0) {   
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      int i = 0;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);     // case 0
         TEST_EQUALITY(data[2],-3.0/2.0); // case 1
         TEST_EQUALITY(data[5],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const CrsMatrixType> subA = getSubBlock_tp(i,j,*A);

            subA->getLocalRowCopy(0,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);

            subA->getLocalRowCopy(2,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);

            subA->getLocalRowCopy(5,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }

      i = 1;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],0.0);

         for(int j=0;j<numBlocks;j++) {
            RCP<const CrsMatrixType> subA = getSubBlock_tp(i,j,*A);

            subA->getLocalRowCopy(0,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else if(myRank==1) {
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      int i = 0;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[3],-3.0);     // case 0
         TEST_EQUALITY(data[0],-3.0/2.0); // case 1
         TEST_EQUALITY(data[6],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const CrsMatrixType> subA = getSubBlock_tp(i,j,*A);

            subA->getLocalRowCopy(3,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
   
            subA->getLocalRowCopy(0,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);
   
            subA->getLocalRowCopy(6,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }

      i = 1;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);

         for(int j=0;j<numBlocks;j++) {
            RCP<const CrsMatrixType> subA = getSubBlock_tp(i,j,*A);

            subA->getLocalRowCopy(0,indices,values,numEntries);
            for(std::size_t k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
         }
      }
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, exclusion)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // pauseToAttach();

   int numBlocks = 3;
   int myRank = comm->getRank();
   int numProc = comm->getSize();

   RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(myRank,numProc,numBlocks);
   BLOFact factory(comm,blkIndexer);
 
   // exclude some pairs
   std::vector<std::pair<int,int> > exPairs;
   exPairs.push_back(std::make_pair(0,2));
   exPairs.push_back(std::make_pair(2,1));
   factory.addExcludedPairs(exPairs);
   factory.addExcludedPair(1,2);

   RCP<LinearObjContainer> container = factory.buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = factory.buildGhostedLinearObjContainer();
   TEST_ASSERT(container!=Teuchos::null);
   TEST_ASSERT(ghosted!=Teuchos::null);

   RCP<BLOC> bContainer = rcp_dynamic_cast<BLOC>(container);
   RCP<BLOC> b_ghosted = rcp_dynamic_cast<BLOC>(ghosted);
   TEST_ASSERT(bContainer!=Teuchos::null);
   TEST_ASSERT(b_ghosted!=Teuchos::null);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
      factory.initializeContainer(LinearObjContainer::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);

      RCP<Thyra::BlockedLinearOpBase<double> > blo 
         = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(bContainer->get_A());

      TEST_ASSERT(!blo->getNonconstBlock(0,0).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(0,1).is_null());
      TEST_ASSERT(blo->getNonconstBlock(0,2).is_null());

      TEST_ASSERT(!blo->getNonconstBlock(1,0).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(1,1).is_null());
      TEST_ASSERT(blo->getNonconstBlock(1,2).is_null());

      TEST_ASSERT(!blo->getNonconstBlock(2,0).is_null());
      TEST_ASSERT(blo->getNonconstBlock(2,1).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(2,2).is_null());
   }
}


// FE assembly tests
/////////////////////////////////////////////////////////////////////

typedef Tpetra::FECrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> FECrsMatrixType;

// Run one AssemblyEngine-style cycle: sum 1.0 into every (row,col) entry the element
// traversal touches -- the same sparsity a scatter evaluator produces -- for every block of
// the ghosted matrix, then migrate ghost -> owned. Returns the owned row sums, block by
// block, which is a value-level summary that must agree between the classic and FE paths.
std::vector<std::vector<double> >
assembleOnesAndCollectOwnedRowSums(const Teuchos::RCP<BLOFact> & la_factory,
                                   const Teuchos::RCP<const panzer::BlockedDOFManager> & blkIndexer,
                                   int numBlocks,
                                   const Teuchos::RCP<LinearObjContainer> & global,
                                   const Teuchos::RCP<LinearObjContainer> & ghosted)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // AssemblyEngine's ordering: beginFill(ghosted,global), scatter, ghostToGlobal,
  // beginFill(global), endFill(global), endFill(ghosted).
  //
  // The ghosted operator must be read AFTER beginFill: under FE assembly the ghosted
  // container holds none of its own and borrows the global container's there.
  la_factory->beginFill(*ghosted,*global);

  RCP<Thyra::LinearOpBase<double> > gh_A = rcp_dynamic_cast<BLOC>(ghosted)->get_A();

  std::vector<std::string> elementBlockIds;
  blkIndexer->getFieldDOFManagers()[0]->getElementBlockIds(elementBlockIds);

  for(int i=0;i<numBlocks;i++) {
    for(int j=0;j<numBlocks;j++) {
      RCP<CrsMatrixType> mat = getSubBlock_tp(i,j,*gh_A);
      if(mat==Teuchos::null)
        continue;

      RCP<const panzer::GlobalIndexer> rowProvider = blkIndexer->getFieldDOFManagers()[i];
      RCP<const panzer::GlobalIndexer> colProvider = blkIndexer->getFieldDOFManagers()[j];

      for(std::size_t b=0;b<elementBlockIds.size();b++) {
        const std::vector<LocalOrdinalT> & elements =
            blkIndexer->getFieldDOFManagers()[0]->getElementBlock(elementBlockIds[b]);

        std::vector<GlobalOrdinalT> row_gids, col_gids;
        for(std::size_t e=0;e<elements.size();e++) {
          rowProvider->getElementGIDs(elements[e],row_gids);
          colProvider->getElementGIDs(elements[e],col_gids);

          Teuchos::Array<GlobalOrdinalT> cols(col_gids.begin(),col_gids.end());
          Teuchos::Array<double> vals(col_gids.size(),1.0);
          for(std::size_t r=0;r<row_gids.size();r++)
            mat->sumIntoGlobalValues(row_gids[r],cols,vals);
        }
      }
    }
  }

  la_factory->ghostToGlobalContainer(*ghosted,*global,LinearObjContainer::Mat);
  la_factory->beginFill(*global);
  la_factory->endFill(*global);
  la_factory->endFill(*ghosted,*global);

  RCP<Thyra::LinearOpBase<double> > g_A = rcp_dynamic_cast<BLOC>(global)->get_A();

  std::vector<std::vector<double> > result(numBlocks);
  for(int i=0;i<numBlocks;i++) {
    RCP<const MapType> ownedMap = la_factory->getMap(i);
    std::vector<double> rowSums(ownedMap->getLocalNumElements(),0.0);

    for(int j=0;j<numBlocks;j++) {
      RCP<CrsMatrixType> mat = getSubBlock_tp(i,j,*g_A);
      if(mat==Teuchos::null)
        continue;

      for(size_t lr=0;lr<ownedMap->getLocalNumElements();lr++) {
        GlobalOrdinalT gid = ownedMap->getGlobalElement(lr);
        size_t numEnt = mat->getNumEntriesInGlobalRow(gid);
        typename CrsMatrixType::nonconst_global_inds_host_view_type inds("inds",numEnt);
        typename CrsMatrixType::nonconst_values_host_view_type vals("vals",numEnt);
        size_t n = 0;
        mat->getGlobalRowCopy(gid,inds,vals,n);
        for(size_t k=0;k<n;k++)
          rowSums[lr] += vals(k);
      }
    }
    result[i] = rowSums;
  }

  return result;
}

// The flag is opt-in and it changes what the matrix getters hand back.
//
// Both paths allocate a fresh, disjoint matrix per call. What differs is the ghosted side:
// classic gives out a separate ghosted matrix per block, while FE has no ghosted matrix at
// all -- the ghosted container borrows the owned one at beginFill(ghosted,owned), which is
// what lets endAssembly() migrate ghost rows onto owned rows in place.
TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, fe_opt_in)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   int numBlocks = 2;
   RCP<const panzer::BlockedDOFManager> blkIndexer
         = buildBlockedIndexer64(comm->getRank(),comm->getSize(),numBlocks);

   // default: unchanged classic behavior
   {
      RCP<BLOFact> classic = rcp(new BLOFact(comm,blkIndexer));
      TEST_ASSERT(!classic->useFEAssembly());

      for(int i=0;i<numBlocks;i++) {
         for(int j=0;j<numBlocks;j++) {
            RCP<CrsMatrixType> owned   = classic->getTpetraMatrix(i,j);
            RCP<CrsMatrixType> ghosted = classic->getGhostedTpetraMatrix(i,j);

            TEST_ASSERT(owned.get()!=ghosted.get());
            TEST_ASSERT(rcp_dynamic_cast<FECrsMatrixType>(owned)==Teuchos::null);

            // each call allocates anew
            TEST_ASSERT(classic->getTpetraMatrix(i,j).get()!=owned.get());
         }
      }
   }

   // opt in
   {
      RCP<BLOFact> fe = rcp(new BLOFact(comm,blkIndexer,true));
      TEST_ASSERT(fe->useFEAssembly());

      for(int i=0;i<numBlocks;i++) {
         for(int j=0;j<numBlocks;j++) {
            RCP<CrsMatrixType> owned = fe->getTpetraMatrix(i,j);
            TEST_ASSERT(rcp_dynamic_cast<FECrsMatrixType>(owned)!=Teuchos::null);

            // a fresh object per call, exactly like the classic getter
            TEST_ASSERT(fe->getTpetraMatrix(i,j).get()!=owned.get());

            // ...but they share the cached graph, so only values are reallocated
            TEST_ASSERT(fe->getFEGraph(i,j).get()==fe->getFEGraph(i,j).get());

            // There is no ghosted matrix to hand out under FE assembly.
            TEST_THROW(fe->getGhostedTpetraMatrix(i,j),std::logic_error);

            // The graph rests in its OWNED view once endAssembly() has run, so this is the
            // owned row map, not the ghosted one.
            TEST_ASSERT(fe->getFEGraph(i,j)->getRowMap()->isSameAs(*fe->getMap(i)));
         }
      }

      // The owned operator must declare stable spaces. The matrix's own getRangeMap()/
      // getDomainMap() cannot be used for this: they follow whichever view is active, which
      // flips across the assembly cycle (and a freshly built FECrsMatrix starts in
      // owned+shared, not owned). So the factory states the spaces explicitly instead of
      // letting Thyra deduce them from the matrix.
      RCP<Thyra::BlockedLinearOpBase<double> > ownedOp
         = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(fe->getThyraMatrix(),true);

      for(int i=0;i<numBlocks;i++) {
         for(int j=0;j<numBlocks;j++) {
            TEST_ASSERT(ownedOp->getBlock(i,j)->range()->isCompatible(*fe->getThyraRangeSpace(i)));
            TEST_ASSERT(ownedOp->getBlock(i,j)->domain()->isCompatible(*fe->getThyraDomainSpace(j)));
         }
      }

      // ...and there is no ghosted operator at all.
      TEST_THROW(fe->getGhostedThyraMatrix(),std::logic_error);
   }
}

// The two paths must produce the same owned matrix. The pre-existing blocked ghostToGlobal
// and graph_constr tests only print their matrices, so this is also the first value-level
// check of the blocked ghost->global migration.
TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, fe_matches_classic_assembly)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   int numBlocks = 2;
   RCP<const panzer::BlockedDOFManager> blkIndexer
         = buildBlockedIndexer64(comm->getRank(),comm->getSize(),numBlocks);

   std::vector<std::vector<double> > classicSums, feSums;

   {
      RCP<BLOFact> la_factory = rcp(new BLOFact(comm,blkIndexer));
      RCP<LinearObjContainer> global  = la_factory->buildLinearObjContainer();
      RCP<LinearObjContainer> ghosted = la_factory->buildGhostedLinearObjContainer();
      la_factory->initializeContainer(LinearObjContainer::Mat,*global);
      la_factory->initializeGhostedContainer(LinearObjContainer::Mat,*ghosted);
      rcp_dynamic_cast<BLOC>(ghosted)->initializeMatrix(0.0);

      classicSums = assembleOnesAndCollectOwnedRowSums(la_factory,blkIndexer,numBlocks,global,ghosted);
   }

   {
      RCP<BLOFact> la_factory = rcp(new BLOFact(comm,blkIndexer,true));
      RCP<LinearObjContainer> global  = la_factory->buildLinearObjContainer();
      RCP<LinearObjContainer> ghosted = la_factory->buildGhostedLinearObjContainer();
      la_factory->initializeContainer(LinearObjContainer::Mat,*global);
      la_factory->initializeGhostedContainer(LinearObjContainer::Mat,*ghosted);
      rcp_dynamic_cast<BLOC>(ghosted)->initializeMatrix(0.0);

      feSums = assembleOnesAndCollectOwnedRowSums(la_factory,blkIndexer,numBlocks,global,ghosted);
   }

   TEST_EQUALITY(classicSums.size(),feSums.size());
   for(std::size_t blk=0;blk<classicSums.size();blk++) {
      TEST_EQUALITY(classicSums[blk].size(),feSums[blk].size());
      for(std::size_t r=0;r<classicSums[blk].size();r++) {
         out << "  block " << blk << " owned row " << r
             << ": classic=" << classicSums[blk][r] << " fe=" << feSums[blk][r] << std::endl;
         TEST_FLOATING_EQUALITY(feSums[blk][r],classicSums[blk][r],1e-14);
      }
   }
}

// Blocked counterpart of tTpetraLinearObjFactory::fe_repeated_assembly_is_not_polluted.
//
// Assembling twice with the matrix zeroed in between must give the same answer both times.
// For an FECrsMatrix that is not automatic: its owned view aliases only the leading chunk of
// the owned+shared values, so the pre-assembly initializeMatrix(0.0) reaches owned rows only
// and would otherwise leave ghost rows holding the previous cycle's contributions. Each
// block needs its ghost rows cleared, which is what the factory's beginFill() does.
TEUCHOS_UNIT_TEST(tBlockedTpetraLinearObjFactory, fe_repeated_assembly_is_not_polluted)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Teuchos::MpiComm<int> > comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
   #else
      NOPE_PANZER_DOESNT_SUPPORT_SERIAL
   #endif

   int numBlocks = 2;
   RCP<const panzer::BlockedDOFManager> blkIndexer
         = buildBlockedIndexer64(comm->getRank(),comm->getSize(),numBlocks);

   RCP<BLOFact> la_factory = rcp(new BLOFact(comm,blkIndexer,true));

   RCP<LinearObjContainer> global  = la_factory->buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = la_factory->buildGhostedLinearObjContainer();
   la_factory->initializeContainer(LinearObjContainer::Mat,*global);
   la_factory->initializeGhostedContainer(LinearObjContainer::Mat,*ghosted);

   rcp_dynamic_cast<BLOC>(ghosted)->initializeMatrix(0.0);
   std::vector<std::vector<double> > firstPass
      = assembleOnesAndCollectOwnedRowSums(la_factory,blkIndexer,numBlocks,global,ghosted);

   // exactly what panzer::ModelEvaluator does before the next Jacobian evaluation
   rcp_dynamic_cast<BLOC>(ghosted)->initializeMatrix(0.0);
   std::vector<std::vector<double> > secondPass
      = assembleOnesAndCollectOwnedRowSums(la_factory,blkIndexer,numBlocks,global,ghosted);

   TEST_EQUALITY(firstPass.size(),secondPass.size());
   for(std::size_t blk=0;blk<firstPass.size();blk++) {
      TEST_EQUALITY(firstPass[blk].size(),secondPass[blk].size());
      for(std::size_t r=0;r<firstPass[blk].size();r++) {
         out << "  block " << blk << " owned row " << r
             << ": first=" << firstPass[blk][r] << " second=" << secondPass[blk][r] << std::endl;
         TEST_FLOATING_EQUALITY(secondPass[blk][r],firstPass[blk][r],1e-14);
      }
   }
}

}

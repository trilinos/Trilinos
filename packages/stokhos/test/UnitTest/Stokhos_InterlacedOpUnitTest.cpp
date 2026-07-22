// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCP.hpp>

#include "Stokhos_InterlacedTestSupport.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "Stokhos_SGModelEvaluator_Interlaced.hpp"
#include "Stokhos_InterlacedOperator.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "EpetraExt_BlockUtility.h"
#include "EpetraExt_RowMatrixOut.h"

TEUCHOS_UNIT_TEST(interlaced_op, test)
{
#ifdef HAVE_MPI
   Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
   Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

   //int rank = comm->MyPID();
   int numProc = comm->NumProc();

   int num_KL = 1;
   int porder = 5;
   bool full_expansion = false;

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
   Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data;
   Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion; 
   {
      if(full_expansion)
         Cijk = basis->computeTripleProductTensor();
      else
        Cijk = basis->computeLinearTripleProductTensor();
   
      Teuchos::ParameterList parallelParams;
      parallelParams.set("Number of Spatial Processors", numProc);
      sg_parallel_data = Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, comm,
                                                                parallelParams));

      expansion = Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis,
		 							             Cijk));
   }
   Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = sg_parallel_data->getMultiComm();

   // determinstic PDE graph
   Teuchos::RCP<Epetra_Map> determRowMap = Teuchos::rcp(new Epetra_Map(-1,10,0,*comm));
   Teuchos::RCP<Epetra_CrsGraph> determGraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*determRowMap,1));
   for(int row=0;row<determRowMap->NumMyElements();row++) {
      int gid = determRowMap->GID(row);
      determGraph->InsertGlobalIndices(gid,1,&gid);
   }
   for(int row=1;row<determRowMap->NumMyElements()-1;row++) {
      int gid = determRowMap->GID(row);
      int indices[2] = {gid-1,gid+1};
      determGraph->InsertGlobalIndices(gid,2,indices);
   }
   determGraph->FillComplete();
   
   Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
   params->set("Scale Operator by Inverse Basis Norms", false);
   params->set("Include Mean", true);
   params->set("Only Use Linear Terms", false);

   Teuchos::RCP<Stokhos::EpetraSparse3Tensor> epetraCijk = 
         Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis,Cijk,sg_comm));
   Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly> W_sg_blocks = 
     Teuchos::rcp(new Stokhos::EpetraOperatorOrthogPoly(basis, epetraCijk->getStochasticRowMap(), determRowMap, determRowMap, sg_comm));
   for(int i=0; i<W_sg_blocks->size(); i++) {
      Teuchos::RCP<Epetra_CrsMatrix> crsMat = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*determGraph));
      crsMat->PutScalar(1.0 + i);
      W_sg_blocks->setCoeffPtr(i,crsMat); // allocate a bunch of matrices   
   }

   Teuchos::RCP<const Epetra_Map> sg_map = 
     Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		    *determRowMap, *(epetraCijk->getStochasticRowMap()), 
		    *(epetraCijk->getMultiComm())));

   // build an interlaced operator (object under test) and a benchmark
   // fully assembled operator
   ///////////////////////////////////////////////////////////////////////

   Stokhos::InterlacedOperator op(sg_comm,basis,epetraCijk,determGraph,params);
   op.PutScalar(0.0);
   op.setupOperator(W_sg_blocks);  

   Stokhos::FullyAssembledOperator full_op(sg_comm,basis,epetraCijk,determGraph,sg_map,sg_map,params);
   full_op.PutScalar(0.0);
   full_op.setupOperator(W_sg_blocks);  

   // here we test interlaced operator against the fully assembled operator
   ///////////////////////////////////////////////////////////////////////
   bool result = true;
   for(int i=0;i<100;i++) {
      // build vector for fully assembled operator (blockwise)
      Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> x_vec_blocks = 
            Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis,epetraCijk->getStochasticRowMap(),determRowMap,epetraCijk->getMultiComm()));
      Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> f_vec_blocks = 
            Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis,epetraCijk->getStochasticRowMap(),determRowMap,epetraCijk->getMultiComm()));
      Teuchos::RCP<Epetra_Vector> x_vec_blocked = x_vec_blocks->getBlockVector(); 
      Teuchos::RCP<Epetra_Vector> f_vec_blocked = f_vec_blocks->getBlockVector(); 
      x_vec_blocked->Random();       // build an initial vector
      f_vec_blocked->PutScalar(0.0);

      // build interlaced vectors
      Teuchos::RCP<Epetra_Vector> x_vec_inter = Teuchos::rcp(new Epetra_Vector(op.OperatorDomainMap()));
      Teuchos::RCP<Epetra_Vector> f_vec_inter = Teuchos::rcp(new Epetra_Vector(op.OperatorRangeMap()));
      Teuchos::RCP<Epetra_Vector> f_vec_blk_inter = Teuchos::rcp(new Epetra_Vector(op.OperatorRangeMap()));
      Stokhos::SGModelEvaluator_Interlaced::copyToInterlacedVector(*x_vec_blocks,*x_vec_inter); // copy random x to 
      f_vec_inter->PutScalar(0.0);

      full_op.Apply(*x_vec_blocked,*f_vec_blocked); 
      op.Apply(*x_vec_inter,*f_vec_inter); 

      // copy blocked action to interlaced for comparison
      Stokhos::SGModelEvaluator_Interlaced::copyToInterlacedVector(*f_vec_blocks,*f_vec_blk_inter); 

      // compute norm
      double error = 0.0;
      double true_norm = 0.0;
      f_vec_blk_inter->NormInf(&true_norm);
      f_vec_blk_inter->Update(-1.0,*f_vec_inter,1.0);
      f_vec_blk_inter->NormInf(&error);

      out << "rel error = " << error/true_norm << " ( " << true_norm << " ), ";
      result &= (error/true_norm < 1e-14);
   }
   out << std::endl;

   TEST_ASSERT(result);
}

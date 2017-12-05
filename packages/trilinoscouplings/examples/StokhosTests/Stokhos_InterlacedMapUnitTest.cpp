/*
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
*/

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCP.hpp>

#include "Stokhos_InterlacedTestSupport.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "Stokhos_SGModelEvaluator_Interlaced.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "EpetraExt_BlockVector.h"

TEUCHOS_UNIT_TEST(map_test, uniform_buildInterlacedMap)
{
#ifdef HAVE_MPI
   Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
   Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

   int rank = comm->MyPID();

   int stochaUnks = 5; 
   Epetra_LocalMap stocha_map(stochaUnks,0,*comm);
   Epetra_Map determ_map(-1,(rank+1)*10,0,*comm);

   Teuchos::RCP<Epetra_Map> obj_ut = Stokhos::SGModelEvaluator_Interlaced::buildInterlaceMap(determ_map,stocha_map);

   TEST_EQUALITY(obj_ut->NumMyElements(),determ_map.NumMyElements()*stochaUnks);
   TEST_EQUALITY(obj_ut->NumGlobalElements(),determ_map.NumGlobalElements()*stochaUnks);

   bool result = true;
   for(int s=0;s<stocha_map.NumMyElements();s++) {
      for(int d=0;d<determ_map.NumMyElements();d++) {
         result &= (obj_ut->GID(stochaUnks*d+s)==(stochaUnks*determ_map.GID(d)+s));
      }
   } 
   TEST_ASSERT(result);
}

TEUCHOS_UNIT_TEST(map_test, copyToInterlace)
{
#ifdef HAVE_MPI
   Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
   Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

   //int rank = comm->MyPID();
   int numProc = comm->NumProc();

   int num_KL = 1;
   int porder = 1;
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

   Teuchos::RCP<Stokhos::EpetraSparse3Tensor> epetraCijk = 
         Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis,Cijk,sg_comm));

   Epetra_Map determRowMap(-1,3,0,*comm);
   Teuchos::RCP<Epetra_Map> determRowMap_rcp = Teuchos::rcpFromRef(determRowMap);
   Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> x_vec_blocks = 
         Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis,epetraCijk->getStochasticRowMap(),determRowMap_rcp,epetraCijk->getMultiComm()));
   for(int b=0;b<porder+1;b++) 
      x_vec_blocks->getBlockVector()->GetBlock(b)->Print(std::cout);   
//      for(int i=0;i<3;i++)
//         (*x_vec_blocks->getBlockVector()->GetBlock(b))[i] = 1.0+i + (b+1)*3;
   x_vec_blocks->getBlockVector()->Print(std::cout);   
}

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
#include "Stokhos_BasisInteractionGraph.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Stokhos_ParallelData.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
#endif

TEUCHOS_UNIT_TEST(basis_interaction_graph, test_square)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
   #endif

   int numProc = comm->NumProc();

   int num_KL = 3;
   int porder = 3;

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor();
   Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data;

   {
      Teuchos::ParameterList parallelParams;
      parallelParams.set("Number of Spatial Processors", numProc);
      sg_parallel_data = Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, comm,
                                                                 parallelParams));
   }
   Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = sg_parallel_data->getMultiComm();
   Teuchos::RCP<Stokhos::EpetraSparse3Tensor> epetraCijk = 
         Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis,Cijk,sg_comm));

   Stokhos::BasisInteractionGraph big(*basis);   

   // for grins print out the graph
   out << "Size = " << big.rowCount() << " x " << big.colCount() << std::endl;

   for(std::size_t i=0;i<big.rowCount();i++) {
      for(std::size_t j=0;j<big.colCount();j++)
         out << " " << (big(i,j) ? "*" : " " ) << " ";
      out << std::endl;
   }
   out << std::endl;

   // check the graph for correctness
   Teuchos::RCP<const Epetra_CrsGraph> graph = epetraCijk->getStochasticGraph();
   TEST_EQUALITY(graph->NumMyRows(),graph->NumGlobalRows());
   TEST_EQUALITY(int(big.rowCount()),graph->NumMyRows());
   TEST_EQUALITY(int(big.colCount()),graph->NumMyCols());

   for(int i=0;i<graph->NumGlobalRows();i++) {
      bool rowPassed = true;
      int count = 0;
      std::vector<bool> row(graph->NumGlobalRows(),false);
      std::vector<int> activeRow(graph->NumGlobalRows());
      graph->ExtractGlobalRowCopy(i,graph->NumGlobalRows(),count,&activeRow[0]);

      // get truth graph
      for(int j=0;j<count;j++) 
         row[activeRow[j]] = true;

      // check active row
      for(std::size_t j=0;j<row.size();j++) 
         rowPassed = rowPassed && (row[j]==big(i,j));

      TEST_ASSERT(rowPassed);
   }
}

TEUCHOS_UNIT_TEST(basis_interaction_graph, test_isotropic_rect)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
   #endif

   int num_KL = 2;
   int porder = 3;

   Teuchos::RCP<const Stokhos::ProductBasis<int,double> > masterBasis = buildBasis(num_KL,porder);
   Teuchos::RCP<const Stokhos::ProductBasis<int,double> > rowBasis = buildBasis(num_KL,2);
   Teuchos::RCP<const Stokhos::ProductBasis<int,double> > colBasis = buildBasis(num_KL,3);

   out << "Master Array Basis = \n";
   for(int i=0;i<masterBasis->size();i++) {
     Stokhos::MultiIndex<int> masterArray = masterBasis->term(i);
      for(int i=0;i<num_KL;i++) { 
         out << masterArray[i] << " ";
      }
      out << std::endl;
   }

   out << "Row Array Basis = \n";
   for(int i=0;i<rowBasis->size();i++) {
     Stokhos::MultiIndex<int> rowArray = rowBasis->term(i);
      for(int i=0;i<num_KL;i++) { 
         out << rowArray[i] << " ";
      }
      out << std::endl;
   }

   Stokhos::BasisInteractionGraph masterBig(*masterBasis);   
   Stokhos::BasisInteractionGraph rectBig(*masterBasis,*rowBasis,*colBasis);   

   out << "rowBasis.size = " << rowBasis->size() << std::endl;
   out << "colBasis.size = " << colBasis->size() << std::endl;

   // for grins print out the graph
   out << "Size = " << rectBig.rowCount() << " x " << rectBig.colCount() << std::endl;
   for(std::size_t i=0;i<rectBig.rowCount();i++) {
      for(std::size_t j=0;j<rectBig.colCount();j++)
         out << " " << (rectBig(i,j) ? "*" : " " ) << " ";
      out << std::endl;
   }
   out << std::endl;

   out << "Size = " << masterBig.rowCount() << " x " << masterBig.colCount() << std::endl;
   for(std::size_t i=0;i<masterBig.rowCount();i++) {
      for(std::size_t j=0;j<masterBig.colCount();j++)
         out << " " << (masterBig(i,j) ? "*" : " " ) << " ";
      out << std::endl;
   }
   out << std::endl;

   // loop over rectangle graph making sure it matches with the master graph
   bool graphs_match = true;
   for(std::size_t i=0;i<rectBig.rowCount();i++) {
      std::size_t masterI = masterBasis->index(rowBasis->term(i));
      for(std::size_t j=0;j<rectBig.colCount();j++) {
         std::size_t masterJ = masterBasis->index(colBasis->term(j));
         graphs_match &= (rectBig(i,j)==masterBig(masterI,masterJ));
      }
   }

   TEST_ASSERT(graphs_match);
}

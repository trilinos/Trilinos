// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Zoltan2_componentMetrics.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <iostream>
#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>


/////////////////////////////////////////////////////////////////////////////
// Test program for componentMetrics code.
// Usage:
//     a.out 
// Karen Devine, 2017
/////////////////////////////////////////////////////////////////////////////

typedef Tpetra::Map<zlno_t, zgno_t> zmap_t;
typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t> zmatrix_t;
typedef Tpetra::CrsGraph<zlno_t, zgno_t> zgraph_t;

typedef Zoltan2::XpetraCrsMatrixAdapter<zmatrix_t> matrixAdapter_t;
typedef Zoltan2::XpetraCrsGraphAdapter<zgraph_t> graphAdapter_t;

/////////////////////////////////////////////////////////////////////////////
int test_no_graph(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  std::string name("no_graph ");
  int ierr = 0;

  // Create a default map
  const size_t gNvtx = 25;

  Teuchos::RCP<const zmap_t> map = rcp(new zmap_t(gNvtx, 0, comm));
  size_t nVtx = map->getNodeNumElements();

  // Create a Tpetra::Matrix with no edges
  size_t maxRowLen = 1;
  Teuchos::RCP<const zmatrix_t> matrix = rcp(new zmatrix_t(map, maxRowLen));

  // Create a Zoltan2 XpetraGraphAdapter
  matrixAdapter_t ia(matrix, 0);

  // Call connected components utility
  Zoltan2::perProcessorComponentMetrics<matrixAdapter_t> cc(ia, *comm);

  // Check result:  
  // With no edges, every vertex should be a component
  if (nVtx != cc.getNumComponents()) {
    std::cout << name << "Invalid number of components "
              << cc.getNumComponents() << " should be " << nVtx << std::endl;
    ierr++;
  }

  // With no edges, very component should have size one
  if (cc.getMaxComponentSize() != 1) {
    std::cout << name << "Maximum component size "
              << cc.getMaxComponentSize() << " should be one" << std::endl;
    ierr++;
  }

  if (cc.getMinComponentSize() != 1) {
    std::cout << name << "Minimum component size "
              << cc.getMinComponentSize() << " should be one" << std::endl;
    ierr++;
  }

  if (cc.getAvgComponentSize() != 1.) {
    std::cout << name << "Average component size "
              << cc.getAvgComponentSize() << " should be one" << std::endl;
    ierr++;
  }

  return ierr;
}

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  // Establish session.
  Teuchos::GlobalMPISession mpiSession(&narg, &arg, NULL);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
           Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  int me = comm->getRank();

  int testReturn = 0;
//testReturn += test_empty_procs(comm);      // all data on rank 0
  testReturn += test_no_graph(comm);         // no edges in graph
//testReturn += test_single_component(comm); // one component per rank
//testReturn += test_every_third(comm);      // every 3rd vtx in same component

  int gtestReturn = 0;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 1,
                               &testReturn, &gtestReturn);
  if (me == 0) {
    if (gtestReturn) std::cout << "FAIL" << std::endl;
    else             std::cout << "PASS" << std::endl;
  }

  return gtestReturn;
}

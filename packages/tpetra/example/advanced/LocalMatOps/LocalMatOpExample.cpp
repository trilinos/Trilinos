/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER
*/

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <examples/KokkosExamples_EmptySparseKernelClass.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix.hpp"

/** \file LocalMatOpsExample.cpp
    \brief A file testing the build of the KokkosExamples::EmptySparseKernel and illustrating a custom sparse mat-vec with Tpetra::CrsMatrix.
 */

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  typedef Tpetra::DefaultPlatform::DefaultPlatformType              Platform;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType    Node;
  typedef KokkosExamples::EmptySparseKernel<void,Node>              SparseOps;
  typedef Tpetra::Map<int,int,Node>                                 Map;
  typedef Tpetra::CrsMatrix<float,int,int,Node,SparseOps>           Matrix;
  typedef Tpetra::MultiVector<float,int,int,Node>                   MultiVector;

  // 
  // Get the default communicator and node
  //
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();
  Teuchos::RCP<Node>             node = platform.getNode();

  std::cout << "Note, this class doesn't actually do anything. We are only testing that it compiles." << std::endl;

  // create the matrix with a custom sparse mat-vec kernel
  const size_t numRows = 5;
  const Tpetra::global_size_t numGlobalRows = numRows*comm->getSize();
  Teuchos::RCP<const Map> map = Tpetra::createUniformContigMapWithNode<int,int,Node>(numGlobalRows, comm, node);
  Teuchos::RCP<Matrix> matrix = Teuchos::rcp( new Matrix(map,1,Tpetra::DynamicProfile) );
  matrix->fillComplete(Teuchos::null);

  Teuchos::RCP<MultiVector> X = Tpetra::createMultiVector<float>(map, 2),
                            Y = Tpetra::createMultiVector<float>(map, 2);
  matrix->apply(*X, *Y);

  std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;
}

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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Assert.hpp>

#include <Kokkos_TPINode.hpp>

// I/O for Harwell-Boeing files
#include "Tpetra_MatrixIO.hpp"

#include "Tpetra_Version.hpp"
#include "CRSTiming.hpp"

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  //
  // Get example parameters from command-line processor
  //
  int numThreads = -1;
  std::string filename("bcsstk14.hb");
  int verbose = 1;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("num-threads",&numThreads,"Number of threads.");
  cmdp.setOption("verbose",&verbose,"Verbose (zero for silent).");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  //
  // Say hello, print some communicator info
  //
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
  if (comm->getRank () == 0) {
    std::cout << "\n" << Tpetra::version() << std::endl << std::endl;
    std::cout << argv[0] << filename << std::endl;
    std::cout << "Comm info: " << *comm;
  }

  typedef KokkosClassic::TPINode Node;
  Teuchos::ParameterList params;
  params.set<int>("Num Threads",numThreads);
  params.set<int>("Verbose",verbose);
  Teuchos::RCP<Node> node = Teuchos::rcp(new Node(params));

  //
  // Read Tpetra::CrsMatrix from file
  //
  Teuchos::RCP< Tpetra::CrsMatrix<double,int,int,Node> > A;
  Tpetra::Utils::readHBMatrix(filename,comm,node,A);
  if (comm->getRank() == 0 && verbose) {
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  CRSTiming<Node>(A);

  return 0;
}



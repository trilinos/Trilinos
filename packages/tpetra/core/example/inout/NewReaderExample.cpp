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
// ************************************************************************
// @HEADER
*/

// Small example showing how to use the new MatrixMarket reader to read
// a file into a Tpetra::CrsMatrix.

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"


int main(int narg, char *arg[]) 
{
  // Create Tpetra scope (calls Kokkos::Initialize and Tpetra::Initialize
  Tpetra::ScopeGuard scope(&narg, &arg);

  // Get a default communicator:  MPI_COMM_WORLD
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Command-line options for this example:  filename, distribution, etc.
  // You can specify options any way you like
  Teuchos::CommandLineProcessor cmdp(false,true);

  std::string filename = "";
  cmdp.setOption("file", &filename,
                 "Path and filename of the matrix to be read.");

  std::string distribution ="1D";
  cmdp.setOption("distribution", &distribution,
                 "Parallel distribution to use: "
                 "1D, 2D, LowerTriangularBlock, MMFile");

  bool randomize = false;
  cmdp.setOption("randomize", "norandomize", &randomize,
                 "Randomly permute the matrix rows and columns");

  bool symmetrize = false;  // Not yet tested
  cmdp.setOption("symmetrize", "nosymmetrize", &symmetrize,
                 "Symmetrize the matrix as it is being read");

  std::string diagonal = "";
  cmdp.setOption("diagonal", &diagonal,
                 "How to manipulate the matrix diagonal entries, if necessary");

  bool binary = false;  
  cmdp.setOption("binary", "mtx", &binary,
                 "Reading a binary file instead of a matrix market file");

  int chunkSize = 500;
  cmdp.setOption("chunksize", &chunkSize,
		 "Number of edges to be read and broadcasted at once");

  if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Load the options into a Teuchos::Parameter list
  Teuchos::ParameterList params;
  params.set("distribution", distribution);
  params.set("randomize", randomize);
  params.set("symmetrize", symmetrize);
  params.set("diagonal", diagonal);
  params.set("binary", binary);
  params.set("chunkSize", (size_t)chunkSize);

  // Call readSparseFile to read the file
  using scalar_t = Tpetra::CrsMatrix<>::scalar_type;
  using matrix_t = Tpetra::CrsMatrix<scalar_t>;
  Teuchos::RCP<matrix_t> Amat;
  try {
    using reader_t = Tpetra::MatrixMarket::Reader<matrix_t>;
    Amat = reader_t::readSparseFile(filename, comm, params);
  }
  catch (std::exception &e) {
    if (comm->getRank() == 0) {
        std::cout << ":  matrix reading failed " << filename << std::endl;
        std::cout << e.what() << std::endl;
      }
    throw e;
  }

  // The resulting matrix is ready to use
  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));
  Amat->describe(foo, Teuchos::VERB_EXTREME);

  // The end
  if (comm->getRank() == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}


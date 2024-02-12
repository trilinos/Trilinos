// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_MatrixTraits.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"

#ifndef GALERI_TEST_USE_LONGLONG_GO
#define GO int
#else
#define GO long long
#endif
#define Scalar int
#define LO int
#define Node Tpetra::KokkosClassic::DefaultNode::DefaultNodeType

using namespace Galeri;

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef Tpetra::Map<LO, GO, Node> Tpetra_Map;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> Tpetra_CrsMatrix;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Create comm
  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  std::string mapType = "Cartesian2D";
  auto mapParameters = Teuchos::ParameterList("Tpetra::Map");
  // here we specify the global dimension of the problem
  int nx = 5 * comm->getSize();
  int ny = 5 * comm->getSize();
  mapParameters.set("nx", nx);
  mapParameters.set("ny", ny);

  try
  {
    // Creation of the map
    auto map = RCP{Galeri::Xpetra::CreateMap<Scalar, GO, Tpetra_Map>(mapType, comm, mapParameters)};

    // Creates a diagonal matrix with 1's on the diagonal
    auto matrix = Galeri::Xpetra::BigStar2D<Scalar, LO, GO, Tpetra_Map, Tpetra_CrsMatrix>(map, nx, ny, 20, -8, -8, -8, -8, 2, 2, 2, 2, 1, 1, 1, 1);

    // print out the matrix
    auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    matrix->describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
  }
  catch (Exception &rhs)
  {
    if (comm->getRank() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();

#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return (EXIT_FAILURE);
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}

// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/**
  \file   SimpleSolve.cpp
  \author Eric Bavier <etbavie@sandia.gov>
  \date   Sat Jul 17 10:35:39 2010

  \brief  Simple example of Amesos2 usage.

  This example solves a simple sparse system of linear equations using the
  Amesos2 interface to the Superlu solver.
  */

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard mpiSession(&argc, &argv);

  typedef Tpetra::MultiVector<> MV;
  typedef MV::scalar_type Scalar;
  typedef MV::local_ordinal_type LO;
  typedef MV::global_ordinal_type GO;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> MAT;

  using std::endl;
  using Teuchos::tuple;
  using Tpetra::global_size_t;

  std::ostream &out              = std::cout;
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(rcpFromRef(out));

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  size_t myRank                       = comm->getRank();

  out << "Amesos2 stand-alone test" << endl
      << endl;

  const size_t numVectors = 1;

  int numGlobalElements                         = 1000;
  RCP<const Tpetra::Map<> > map                 = Tpetra::createUniformContigMap<LO, GO>(numGlobalElements, comm);
  const size_t numMyElements                    = map->getLocalNumElements();
  Teuchos::ArrayView<const GO> myGlobalElements = map->getLocalElementList();

  RCP<MAT> A = Tpetra::createCrsMatrix<Scalar>(map, 3);

  // 1D Laplace
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i],
                            tuple<GO>(myGlobalElements[i], myGlobalElements[i] + 1),
                            tuple<Scalar>(2.0, -1.0));
    } else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i],
                            tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i]),
                            tuple<Scalar>(-1.0, 2.0));
    } else {
      A->insertGlobalValues(myGlobalElements[i],
                            tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i], myGlobalElements[i] + 1),
                            tuple<Scalar>(-1.0, 2.0, -1.0));
    }
  }
  A->fillComplete();

  /* Create X */
  RCP<MV> X = rcp(new MV(map, numVectors));
  Teuchos::ScalarTraits<Scalar>::seedrandom(846930886);
  X->randomize();

  /* Print X norm */
  {
    Teuchos::Array<ST::magnitudeType> norms(1);
    X->norm2(norms);
    if (myRank == 0)
      *fos << "||X_true|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << endl;
  }

  /* Create B  */
  RCP<MV> B = rcp(new MV(map, numVectors));
  A->apply(*X, *B, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);

  /* Reset X */
  X->putScalar((Scalar)0.0);

  // Create solver interface to Superlu through Amesos::Factory
  RCP<Amesos2::Solver<MAT, MV> > solver = Amesos2::create<MAT, MV>("Superlu", A, X, B);

  // Solve
  solver->symbolicFactorization().numericFactorization().solve();

  /* Print X norm */
  {
    Teuchos::Array<ST::magnitudeType> norms(1);
    X->norm2(norms);
    if (myRank == 0)
      *fos << "||X_directSolve|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << endl;
  }

  //   /* Print the solution */
  //   RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(rcpFromRef(out));

  //   *fos << "Solution :" << endl;
  //   X->describe(*fos,Teuchos::VERB_EXTREME);
  //   *fos << endl;

  // We are done.
  return 0;
}

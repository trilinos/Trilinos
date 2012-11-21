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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

// This file is intended to house all the tests for MueLu_Utilities.hpp.

namespace MueLuTests {

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  TEUCHOS_UNIT_TEST(Utilities,MatMatMult_EpetraVsTpetra)
  {
    out << "version: " << MueLu::Version() << std::endl;
    out << "This test compares the matrix matrix multiply between Tpetra and Epetra" << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Teuchos::Array<ST::magnitudeType> normResult1(1);

    //Calculate result = (Op*Op)*X for Tpetra
    int nx = 37*comm->getSize();
    int ny=nx;
    RCP<Matrix> Op = TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build2DPoisson(nx,ny,Xpetra::UseEpetra);
    RCP<Matrix> OpOp = Utils::TwoMatrixMultiply(Op,false,Op,false);
    RCP<MultiVector> result = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    RCP<MultiVector> X = MultiVectorFactory::Build(OpOp->getDomainMap(),1);
    Teuchos::Array<ST::magnitudeType> xnorm(1);
    X->setSeed(8675309);
    X->randomize(true);
    X->norm2(xnorm);
    OpOp->apply(*X,*result,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    result->norm2(normResult1);

    //Calculate result = (Op*Op)*X for Epetra
    Op = TestHelpers::TestFactory<SC, LO, GO, NO, LMO>::Build2DPoisson(nx,ny,Xpetra::UseTpetra);
    OpOp = Utils::TwoMatrixMultiply(Op,false,Op,false);
    result = MultiVectorFactory::Build(OpOp->getRangeMap(),1);
    X = MultiVectorFactory::Build(OpOp->getDomainMap(),1);
    X->setSeed(8675309);
    X->randomize(true);
    X->norm2(xnorm);
    OpOp->apply(*X,*result,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    Teuchos::Array<ST::magnitudeType> normResult2(1);
    result->norm2(normResult2);

    TEST_FLOATING_EQUALITY(normResult1[0], normResult2[0], 1e-12);

  } //EpetraVersusTpetra
#endif

}//namespace MueLuTests


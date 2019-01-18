// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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

#include <string>

#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Core.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2.hpp"

using std::cout;
using std::endl;
using std::string;

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::tuple;
using Teuchos::ScalarTraits;
using Teuchos::OrdinalTraits;
using Teuchos::CONJ_TRANS;
using Teuchos::TRANS;
using Teuchos::NO_TRANS;


using Tpetra::global_size_t;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using Tpetra::Map;

int main(int argc, char *argv[]){
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  typedef float SCALAR;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef std::complex<SCALAR> cmplx;
  typedef CrsMatrix<cmplx,LO,GO> MAT;
  typedef ScalarTraits<cmplx> ST;
  typedef MultiVector<cmplx,LO,GO> MV;
  typedef ST::magnitudeType Mag;
  typedef ScalarTraits<Mag> MT;
  const size_t numVecs = 1;

  std::ostream &out = std::cout;
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  RCP<MAT> A =
    Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../test/matrices/amesos2_test_mat3.mtx",comm);

  RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  RCP<const Map<LO,GO> > rngmap = A->getRangeMap();

  // Create the know-solution vector
  std::map<GO,cmplx> xValues;
  xValues[0] = cmplx(-0.58657, 0.10646);
  xValues[1] = cmplx(0.86716, 0.75421);
  xValues[2] = cmplx(0.58970, 0.29876);

  std::map<GO,cmplx>::iterator it;
  RCP<MV> X = rcp(new MV(dmnmap, numVecs));
  X->setObjectLabel("X");

  for( it = xValues.begin(); it != xValues.end(); ++it ){
    if( rngmap->isNodeGlobalElement( (*it).first ) ){
      X->replaceGlobalValue( (*it).first, 0, (*it).second );
    }
  }

  RCP<MV> Xhat = rcp(new MV(dmnmap, numVecs));
  RCP<MV> B    = rcp(new MV(rngmap, numVecs));
  X->setObjectLabel("X");
  B->setObjectLabel("B");
  Xhat->setObjectLabel("Xhat");

  A->apply(*X,*B,Teuchos::CONJ_TRANS); // use conjugate-transpose

  Xhat->randomize();

  // Solve A*Xhat = B for Xhat using the Superlu solver
  RCP<Amesos2::Solver<MAT,MV> > solver
    = Amesos2::create<MAT,MV>("SuperLU", A, Xhat, B );

  Teuchos::ParameterList amesos2_params;
  // Setting the following will cause Amesos2 to complain
  amesos2_params.sublist("SuperLU").set("Trans","CONJ","Solve with conjugate-transpose");

  solver->setParameters( rcpFromRef(amesos2_params) );
  solver->symbolicFactorization().numericFactorization().solve();

  B->describe(*fos, Teuchos::VERB_EXTREME);
  Xhat->describe(*fos, Teuchos::VERB_EXTREME);
  X->describe(*fos, Teuchos::VERB_EXTREME);

  // Check result of solve
  Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
  Xhat->norm2(xhatnorms());
  X->norm2(xnorms());
}

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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <iterator>


namespace {

  //
  // UNIT TESTS
  //

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, basic, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar scalar = 100. * (me+1);
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos

  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  std::cout << me << " DEFAULT MAP" << std::endl;
  defaultMap->describe(foo, Teuchos::VERB_EXTREME);

  // Create vector
  vector_t defaultVec(defaultMap);
  defaultVec.putScalar(scalar);

  // Check result; all vector entries should be 0
  auto data = defaultVec.getLocalViewHost(Tpetra::Access::ReadOnly());

  for (size_t i = 0; i < defaultVec.getLocalLength(); i++) {
    if (data(i,0) != scalar) { 
      ierr++;
      std::cout << "Expected: " << scalar << ", got: "<< data(i, 0) << std::endl;
    }
  }
    

  if (ierr > 0) 
    std::cout << "TEST FAILED:  DEFAULT-TO-DEFAULT TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}



//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, basic             , LO, GO, SCALAR, NODE ) \

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

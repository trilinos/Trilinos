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
// // Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <MueLu_config.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_AmalgamationInfo.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(AmalgamationInfo, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  LO fullblocksize    = 1;              // block dim for fixed size blocks
  GO offset           = 0;              // global offset of dof gids
  LO blockid          = -1;             // block id in strided map
  LO nStridedOffset   = 0;              // DOF offset for strided block id "blockid" (default = 0)
  LO stridedblocksize = fullblocksize;  // size of strided block id "blockid" (default = fullblocksize, only if blockid!=-1 stridedblocksize <= fullblocksize)

  const GO nx        = 199;
  using test_factory = TestHelpers::TestFactory<SC, LO, GO, NO>;
  RCP<Matrix> A      = test_factory::Build1DPoisson(nx);

  RCP<AmalgamationInfo> amalgamationData;
  RCP<Array<LO> > theRowTranslation = rcp(new Array<LO>);
  RCP<Array<LO> > theColTranslation = rcp(new Array<LO>);

  RCP<AmalgamationInfo> amalgamationInfo = rcp(new AmalgamationInfo(theRowTranslation,
                                                                    theColTranslation,
                                                                    A->getRowMap(),
                                                                    A->getColMap(),
                                                                    A->getColMap(),
                                                                    fullblocksize,
                                                                    offset,
                                                                    blockid,
                                                                    nStridedOffset,
                                                                    stridedblocksize));
  TEST_INEQUALITY(amalgamationInfo, Teuchos::null);

  TEST_EQUALITY(amalgamationInfo->description(), "AmalgamationInfo");
  LO fullBlockSizeExpected;
  LO blockIDExpected;
  LO stridingOffsetExpected;
  LO stridedBlockSizeExpected;
  GO indexBaseExpected;
  amalgamationInfo->GetStridingInformation(fullBlockSizeExpected, blockIDExpected, stridingOffsetExpected, stridedBlockSizeExpected, indexBaseExpected);
  TEST_EQUALITY(fullBlockSizeExpected == fullblocksize, true);
  TEST_EQUALITY(blockIDExpected == blockid, true);
  TEST_EQUALITY(stridingOffsetExpected == nStridedOffset, true);
  TEST_EQUALITY(stridedBlockSizeExpected == stridedblocksize, true);
  TEST_EQUALITY(indexBaseExpected == A->getColMap()->getIndexBase(), true);
  amalgamationInfo->print(*Teuchos::VerboseObjectBase::getDefaultOStream(), MueLu::Extreme);

}  // Constructor

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(AmalgamationInfo, Constructor, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests

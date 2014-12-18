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
#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_SerialPlatform.hpp"
#ifdef HAVE_TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#endif

#include <Kokkos_DefaultNode.hpp>

namespace {

  using Teuchos::RCP;
  using Teuchos::Comm;
  using Teuchos::rcp;
  using Tpetra::SerialPlatform;
#ifdef HAVE_TPETRA_MPI
  using Tpetra::MpiPlatform;
#endif
  typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;

  template <class PLAT>
  RCP<PLAT> getPlatform() {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Platform type " << Teuchos::TypeNameTraits<PLAT>::name() << " not defined.");
  }

  RCP<node_type> snode;

  template <>
  RCP<SerialPlatform<node_type> > getPlatform() {
    if (snode == Teuchos::null) {
      Teuchos::ParameterList pl;
      snode = rcp(new node_type(pl));
    }
    return rcp(new SerialPlatform<node_type>(snode));
  }

#ifdef HAVE_TPETRA_MPI
  template <>
  RCP<MpiPlatform<node_type> > getPlatform() {
    if (snode == Teuchos::null) {
      Teuchos::ParameterList pl;
      snode = rcp(new node_type(pl));
    }
    return rcp(new MpiPlatform<node_type>(snode));
  }
#endif

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Platform, basic, PlatformType )
  {
    out << "Testing " << Teuchos::TypeNameTraits<PlatformType>::name() << std::endl;
    typedef typename PlatformType::NodeType            N;
    // create a platform
    RCP<PlatformType> platform = getPlatform<PlatformType>();
    platform->setObjectLabel("not the default label");
    // get the comm for this platform
    RCP<const Comm<int> > comm = platform->getComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    TEST_EQUALITY( myImageID < numImages, true );
    TEST_EQUALITY_CONST( comm != Teuchos::null, true );
    RCP<N> node  = platform->getNode();
    (void)node;
  }


  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

typedef SerialPlatform<node_type> SP;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, SP)
#ifdef HAVE_TPETRA_MPI
typedef MpiPlatform<node_type> MP;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Platform, basic, MP )
#endif

}



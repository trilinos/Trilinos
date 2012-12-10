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

#ifndef TPETRA_TESTINGUTILITIES_HPP_
#define TPETRA_TESTINGUTILITIES_HPP_

#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_ETIHelperMacros.h>
#include <Tpetra_DefaultPlatform.hpp>

namespace Tpetra {

  namespace TestingUtilities {

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Tpetra::DefaultPlatform;
    using Teuchos::Comm;

    bool testMpi = true;

    using Kokkos::SerialNode;
    RCP<SerialNode> snode;
#ifdef HAVE_KOKKOSCLASSIC_TBB
    using Kokkos::TBBNode;
    RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
    using Kokkos::TPINode;
    RCP<TPINode> tpinode;
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
    using Kokkos::OpenMPNode;
    RCP<OpenMPNode> ompnode;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
    using Kokkos::ThrustGPUNode;
    RCP<ThrustGPUNode> thrustnode;
#endif

    RCP<const Comm<int> > getDefaultComm()
    {
      if (testMpi) {
        return DefaultPlatform::getDefaultPlatform().getComm();
      }
      return rcp(new Teuchos::SerialComm<int>());
    }

    template <class Node>
      RCP<Node> getNode() {
        assert(false);
      }

    template <>
      RCP<SerialNode> getNode<SerialNode>() {
        if (snode == null) {
          Teuchos::ParameterList pl;
          snode = rcp(new SerialNode(pl));
        }
        return snode;
      }

#ifdef HAVE_KOKKOSCLASSIC_TBB
    template <>
      RCP<TBBNode> getNode<TBBNode>() {
        if (tbbnode == null) {
          Teuchos::ParameterList pl;
          pl.set<int>("Num Threads",0);
          tbbnode = rcp(new TBBNode(pl));
        }
        return tbbnode;
      }
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
    template <>
      RCP<TPINode> getNode<TPINode>() {
        if (tpinode == null) {
          Teuchos::ParameterList pl;
          pl.set<int>("Num Threads",0);
          tpinode = rcp(new TPINode(pl));
        }
        return tpinode;
      }
#endif

#ifdef HAVE_KOKKOSCLASSIC_OPENMP
    template <>
      RCP<OpenMPNode> getNode<OpenMPNode>() {
        if (ompnode == null) {
          Teuchos::ParameterList pl;
          pl.set<int>("Num Threads",0);
          ompnode = rcp(new OpenMPNode(pl));
        }
        return ompnode;
      }
#endif

#ifdef HAVE_KOKKOSCLASSIC_THRUST
    template <>
      RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
        if (thrustnode == null) {
          Teuchos::ParameterList pl;
          pl.set<int>("Num Threads",0);
          pl.set<int>("Verbose",1);
          thrustnode = rcp(new ThrustGPUNode(pl));
        }
        return thrustnode;
      }
#endif

  } // end namespace Tpetra::TestingUtilities

} // end namespace Tpetra

#endif // TPETRA_TESTINGUTILITIES_HPP_

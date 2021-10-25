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

// Computes norms in both host and device space to exercise calls 
// to KokkosBlas in both cases.

#include "Tpetra_Core.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"

namespace {


//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug9856, HostNorm, Scalar,LO,GO,Node)
{
  // Compute the norm on the host

  auto comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using mv_t = Tpetra::MultiVector<Scalar,LO,GO,Node>;

  // use small vector length to keep norm computation on host
  const size_t nGlobalEntries = 100;  
  const size_t nVecs = 100;

  Teuchos::RCP<const map_t> map = rcp(new map_t(nGlobalEntries, 0, comm));
  mv_t mv(map, nVecs);
  mv.randomize();

  // grab host view to ensure device is not most up-to-date
  {
    auto hostview = mv.getLocalViewHost(Tpetra::Access::ReadWrite);
  }

  // compute the norm
  std::vector<Scalar> norm(nVecs);
  Teuchos::ArrayView<Scalar> normView(norm);
  TEST_NOTHROW( mv.norm2(normView) );
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug9856, DeviceNorm, Scalar,LO,GO,Node)
{
  // This test reproduces the issue in issue #9856
  // Compute the norm on the device

  auto comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using mv_t = Tpetra::MultiVector<Scalar,LO,GO,Node>;

  // use small vector length 
  const size_t nGlobalEntries = 100;  
  const size_t nVecs = 100;

  Teuchos::RCP<const map_t> map = rcp(new map_t(nGlobalEntries, 0, comm));
  mv_t mv(map, nVecs);
  mv.randomize();

  // grab device view to ensure device is most up-to-date
  {
    auto deviceview = mv.getLocalViewDevice(Tpetra::Access::ReadWrite);
  }

  // compute the norm
  std::vector<Scalar> norm(nVecs);
  Teuchos::ArrayView<Scalar> normView(norm);
  TEST_NOTHROW( mv.norm2(normView) );
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug9856, LongDeviceNorm, Scalar,LO,GO,Node)
{
  // Compute the norm on the device because the vector is long

  auto comm = Tpetra::getDefaultComm();

  using map_t = Tpetra::Map<LO,GO,Node>;
  using mv_t = Tpetra::MultiVector<Scalar,LO,GO,Node>;

  // use small vector length 
  const size_t threshold = 
               Tpetra::Details::Behavior::multivectorKernelLocationThreshold();

std::cout << "KDDKDD threshold " << threshold << std::endl;

  const size_t nGlobalEntries = comm->getSize() * 1.1 * threshold;
  const size_t nVecs = 1;

  Teuchos::RCP<const map_t> map = rcp(new map_t(nGlobalEntries, 0, comm));
  mv_t mv(map, nVecs);
  mv.randomize();

  // grab host view to ensure device is not most up-to-date;
  // norm should run on device anyway because its length warrants a copy
  {
    auto hostview = mv.getLocalViewHost(Tpetra::Access::ReadWrite);
  }

  // compute the norm
  std::vector<Scalar> norm(nVecs);
  Teuchos::ArrayView<Scalar> normView(norm);
  TEST_NOTHROW( mv.norm2(normView) );
}

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug9856, HostNorm, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug9856, DeviceNorm, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug9856, LongDeviceNorm, SCALAR, LO, GO, NODE) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)


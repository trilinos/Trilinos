// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, TeuchosArray, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  using vector_t = Tpetra::MultiVector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar scalar = 100. * (me+1);
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos
  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vector
  vector_t defaultVec(defaultMap, 2);
  defaultVec.putScalar(scalar);

  bool shouldThrow = ! Kokkos::SpaceAccessibility<
                               Kokkos::Serial,
                               typename Node::memory_space>::accessible;
  {
    int threw = false;
    auto dataHost = defaultVec.get1dView();
    try {
      auto dataDevice = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
    } 
    catch (...) {
      threw = true;
    }
    ierr += (threw == shouldThrow) ? 0 : 1;
  }

  {
    int threw = false;
    auto dataHost = defaultVec.get2dView();
    try {
      auto dataDevice = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
    } 
    catch (...) {
      threw = true;
    }
    ierr += (threw == shouldThrow) ? 0 : 1;
  }

  {
    int threw = false;
    auto dataHost = defaultVec.getData(0);
    try {
      auto dataDevice = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
    } 
    catch (...) {
      threw = true;
    }
    ierr += (threw == shouldThrow) ? 0 : 1;
  }

  {
    int threw = false;
    auto dataHost = defaultVec.getDataNonConst(1);
    try {
      auto dataDevice = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
    } 
    catch (...) {
      threw = true;
    }
    ierr += (threw == shouldThrow) ? 0 : 1;
  }

  if (ierr)
    std::cout << "TEST FAILED:  TeuchosArray test had " << ierr 
              << " failures on rank " << me << std::endl;
  
  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, HostView, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar scalar = 100. * (me+1);
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos
  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vector
  vector_t defaultVec(defaultMap);
  defaultVec.putScalar(scalar);

  // Check result; all vector entries should be scalar
  auto data = defaultVec.getLocalViewHost(Tpetra::Access::ReadOnly);

  for (size_t i = 0; i < defaultVec.getLocalLength(); i++) {
    if (data(i,0) != scalar) { 
      ierr++;
      std::cout << "Expected: " << scalar << ", got: "<< data(i, 0) << std::endl;
    }
  }

  if (ierr > 0) 
    std::cout << "TEST FAILED:  HOSTVIEW TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DeviceView, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using device_t = typename vector_t::device_type;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar scalar = 100. * (me+1);
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos
  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vector
  vector_t defaultVec(defaultMap);
  defaultVec.putScalar(scalar);

  // Check result; all vector entries should be the same
  auto data = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto data_old = defaultVec.template getLocalView<device_t>(Tpetra::Access::ReadOnly);

  if (data != data_old) {
    ierr++;
  }
  if (ierr > 0) 
    std::cout << "TEST FAILED:  DeviceView TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//when holding a host view, requesting a device view should fail
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, HostDeviceView, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar scalar = 100. * (me+1);
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos
  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vector
  vector_t defaultVec(defaultMap);
  defaultVec.putScalar(scalar);

  // Check result; all vector entries should be the same
  auto data = defaultVec.getLocalViewHost(Tpetra::Access::ReadOnly);
  bool shouldThrow = ! Kokkos::SpaceAccessibility<Kokkos::Serial, typename Node::memory_space>::accessible;
  int threw = false;
  try {
    auto data_old = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
  } catch (...) {
    threw = true;
  }

  int ierr = (threw == shouldThrow) ? 0 : 1;
  if (ierr)
    std::cout << "TEST FAILED:  HostDeviceView TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//when holding a device view, requesting a local view should fail
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DeviceHostView, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar scalar = 100. * (me+1);
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos
  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vector
  vector_t defaultVec(defaultMap);
  defaultVec.putScalar(scalar);

  // Check result; all vector entries should be the same
  auto data = defaultVec.getLocalViewDevice(Tpetra::Access::ReadOnly);
  bool shouldThrow = !Kokkos::SpaceAccessibility<Kokkos::Serial, typename Node::memory_space>::accessible;
  int threw = false;
  try {
    auto data_old = defaultVec.getLocalViewHost(Tpetra::Access::ReadOnly);
  } catch (...) {
    threw = true;
  }

  int ierr = (threw == shouldThrow) ? 0 : 1;
  if (ierr) 
    std::cout << "TEST FAILED:  DeviceHostView TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, TemplatedGetLocalView, LO, GO, Scalar , Node )
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int np = comm->getSize();

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;
  using WDV = typename vector_t::wrapped_dual_view_type;
  using device_t = typename WDV::DeviceType;
  using host_t = typename WDV::HostType;

  const size_t nGlobalEntries = 8 * np;
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos
  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vector
  vector_t x(defaultMap);

  {
    // Check that getLocalView<device_t> produces a view with the correct device type
    auto deviceView = x.template getLocalView<device_t>(Tpetra::Access::ReadWrite);
    bool correctType = std::is_same<typename decltype(deviceView)::device_type, device_t>::value;
    TEST_ASSERT(correctType);
  }
  constexpr bool needsSyncPath = !std::is_same<Kokkos::HostSpace, typename device_t::memory_space>::value;
  if(needsSyncPath)
  {
    TEST_ASSERT(x.need_sync_host());
    TEST_ASSERT(!x.need_sync_device());
  }
  // Assuming device/host device types aren't the same, make sure getting the host view also works
  {
    auto hostView = x.template getLocalView<host_t>(Tpetra::Access::ReadWrite);
    bool correctType = std::is_same<typename decltype(hostView)::device_type, host_t>::value;
    TEST_ASSERT(correctType);
    if(needsSyncPath)
    {
      TEST_ASSERT(!x.need_sync_host());
      TEST_ASSERT(x.need_sync_device());
    }
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, TeuchosArray, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, HostView, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DeviceView, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, HostDeviceView, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DeviceHostView, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, TemplatedGetLocalView, LO, GO, SCALAR, NODE ) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

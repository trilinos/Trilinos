// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_TestHelpers.hpp
    \brief common code used by tests
*/
#ifndef ZOLTAN2_TESTHELPERS_HPP
#define ZOLTAN2_TESTHELPERS_HPP

#include <Teuchos_UnitTestHarness.hpp>
#include <Zoltan2_Util.hpp>
#include <iostream>

#include <Tpetra_Map.hpp>
    typedef Tpetra::Map<>::node_type znode_t;

// The path to the directory of test data

#define STR_VALUE(path) #path
#define PATH_NAME(path) STR_VALUE(path)

#ifdef Z2_DATA_DIR
std::string testDataFilePath(PATH_NAME(Z2_DATA_DIR));
#else
std::string testDataFilePath(".");
#endif

// The path to the Zoltan1 test directory.  We use
// some of their data for testing.

#ifdef Z1_TEST_DIR
std::string zoltanTestDirectory(PATH_NAME(Z1_TEST_DIR));
#else
std::string zoltanTestDirectory(".");
#endif

//////////////////////////////////////////////////////////////////////////
//
// If Tpetra is compiled with explicit instantiation,
// we have to choose data types that are compiled into Tpetra.
//
// Epetra uses (scalar/lno/gno) == (double/int/int) data types.  If we
// are using these data types, we can test Epetra user input.

// TODO:  KDD 8/13/14
// Global definitions of types gno_t, lno_t, zgid_t and
// scalar_t can cause bugs in the code.  If a class fails to define these
// types, but this file is included before the class file, the types
// from Zoltan2_TestHelpers.hpp will be used in the class.  Compilation in
// user programs (without Zoltan2_TestHelpers.hpp) would then fail.  An
// example of this bug was in the GeometricGenerator class, which used
// scalar_t without defining it.
// In this "fix," I changed gno_t, lno_t, zgid_t, scalar_t, and node_t to
// zgno_t, zlno_t, zzgid_t, zscalar_t and znode_t in Zoltan2_TestHelpers.hpp.
// This change is not the best fix; a better fix would remove the global
// definitions, but that would require more work.  (An even better change
// would use the Teuchos test framework to cycle through various options,
// but that would require even more work and should be delayed until we
// revamp the testing.)

#include <TpetraCore_config.h>

typedef int zpart_t; // added this for consistency but needs further discussion

typedef Tpetra::Map<>::local_ordinal_type zlno_t;
typedef Tpetra::Map<>::global_ordinal_type zgno_t;

using Teuchos::compareArrays;

#ifdef HAVE_TPETRA_DOUBLE
typedef double zscalar_t;
#define HAVE_EPETRA_SCALAR_TYPE
#else
typedef float zscalar_t;
#endif

#if defined HAVE_TPETRA_INT_INT
#if defined HAVE_EPETRA_SCALAR_TYPE
#define HAVE_EPETRA_DATA_TYPES
#endif
#endif

#ifndef HAVE_ZOLTAN2_EPETRA
#undef HAVE_EPETRA_SCALAR_TYPE
#undef HAVE_EPETRA_DATA_TYPES
#endif

//////////////////////////////////////////////////////////////////////////

#define MEMORY_CHECK(iPrint, msg)                                              \
  if (iPrint) {                                                                \
    long kb = Zoltan2::getProcessKilobytes();                                  \
    std::cout.width(10);                                                       \
    std::cout.fill('*');                                                       \
    std::cout << kb << " KB, " << msg << std::endl;                            \
    std::cout.width(0);                                                        \
    std::cout.fill(' ');                                                       \
  }

#define Z2_TEST(TEST)                                                          \
  {                                                                            \
    Teuchos::RCP<Teuchos::FancyOStream> fout =                                 \
        Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));                 \
    auto &out = *fout;                                                         \
    bool success = true;                                                       \
    try {                                                                      \
      TEST;                                                                    \
    } catch (...) {                                                            \
      out << "Test failed.";                                                   \
    }                                                                          \
    if (!success) {                                                            \
      throw std::runtime_error(#TEST " FAIL");                                 \
    }                                                                          \
  }

#define Z2_TEST_THROW(code, ExceptType) Z2_TEST(TEST_THROW(code, ExceptType))
#define Z2_TEST_NOTHROW(code) Z2_TEST(TEST_NOTHROW(code))
#define Z2_TEST_EQUALITY(val1, val2) Z2_TEST(TEST_EQUALITY(val1, val2))
#define Z2_TEST_INEQUALITY(val1, val2) Z2_TEST(TEST_INEQUALITY(val1, val2))
#define Z2_TEST_ASSERT(expr) Z2_TEST(TEST_ASSERT(expr))
#define Z2_TEST_EQUALITY_CONST(val1, val2)                                     \
  Z2_TEST(TEST_EQUALITY_CONST(val1, val2))
#define Z2_TEST_INEQUALITY_CONST(val1, val2)                                   \
  Z2_TEST(TEST_INEQUALITY_CONST(val1, val2))
#define Z2_TEST_COMPARE(val1, comp, val2)                                      \
  Z2_TEST(TEST_COMPARE(val1, comp, val2))
#define Z2_TEST_COMPARE_ARRAYS(val1, val2)                                     \
  Z2_TEST(TEST_COMPARE_ARRAYS(val1, val2))
#define Z2_TEST_COMPARE_FLOATING_ARRAYS(val1, val2, tol)                       \
  Z2_TEST(TEST_COMPARE_FLOATING_ARRAYS(val1, val2, tol))
#define Z2_TEST_FLOATING_EQUALITY(val1, val2, tol)                             \
  Z2_TEST(TEST_FLOATING_EQUALITY(val1, val2, tol))

inline void PrintFromRoot(const std::string &message) {
  if (Tpetra::getDefaultComm()->getRank() == 0) {
    printf("%s \n", message.c_str());
  }
}

template <typename DeviceType, typename HostType>
void TestDeviceHostView(const DeviceType &deviceView,
                        const HostType &hostView) {
  // Should we test for more dimensions?
  for (int dim = 0; dim <= 2; ++dim) {
    Z2_TEST_EQUALITY(deviceView.extent(dim), hostView.extent(dim));
  }

  const auto mirrorDevice = Kokkos::create_mirror_view(deviceView);
  Kokkos::deep_copy(mirrorDevice, deviceView);

  // Compare the values element-wise
  Z2_TEST_COMPARE_ARRAYS(hostView, mirrorDevice);
}

#define Z2_TEST_DEVICE_HOST_VIEWS(deviceView, hostView)                        \
                                                                               \
  {                                                                            \
    for (int dim = 0; dim <= 2; ++dim) {                                       \
      Z2_TEST_EQUALITY(deviceView.extent(dim), hostView.extent(dim));          \
    }                                                                          \
                                                                               \
    const auto mirrorDevice = Kokkos::create_mirror_view(deviceView);          \
    Kokkos::deep_copy(mirrorDevice, deviceView);                               \
                                                                               \
    Z2_TEST_COMPARE_ARRAYS(hostView, mirrorDevice);                            \
  }

#include <ErrorHandlingForTests.hpp>
#include <PrintData.hpp>
#include <UserInputForTests.hpp>

#endif

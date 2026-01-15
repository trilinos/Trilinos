// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   RealSpaceToolsTests.cpp
    \brief  Tests against the RealSpaceTools class.
    \author Nathan V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <Kokkos_Core.hpp>

namespace
{
  using namespace Intrepid2;

  template <ordinal_type D>
  KOKKOS_INLINE_FUNCTION
  double diagonalDeterminant()
  {
    typedef Kokkos::DynRankView<double, Kokkos::MemoryUnmanaged> jac_type;
        
    double jac[D*D];
    jac_type jac_view(jac,D,D);

    for (int d1=0; d1<D; d1++)
      for (int d2=0; d2<D; d2++)
        jac_view(d1,d2) = (d1 == d2) ? 2.0 : 0.0;
    
    double det=Intrepid2::RealSpaceTools<>::Serial::det<D,jac_type>(jac_view);
    return det;
  }

  KOKKOS_INLINE_FUNCTION
  double diagonalDeterminant(const Kokkos::DynRankView<double,DefaultTestDeviceType> &jac_view)
  {
    auto dim = jac_view.extent_int(0);
    typedef Kokkos::DynRankView<double,DefaultTestDeviceType> jac_type;
    
    for (int d1=0; d1<dim; d1++)
      for (int d2=0; d2<dim; d2++)
        jac_view(d1,d2) = (d1 == d2) ? 2.0 : 0.0;
    
    double det=Intrepid2::RealSpaceTools<>::Serial::det<jac_type>(jac_view);
    return det;
  }

  template<int D>
  double getDiagonalDeterminant()
  {
    using ExecSpace = DefaultTestDeviceType::execution_space;
    Kokkos::DynRankView<double,DefaultTestDeviceType> detView("detView",1);
    
    Kokkos::RangePolicy<ExecSpace> policy(ExecSpace(), 0, 1);
    Kokkos::parallel_for("determinant invocation", policy,
    KOKKOS_LAMBDA(int) {
      detView(0) = diagonalDeterminant<D>();
    });
    auto detViewHost = getHostCopy(detView);
    return detViewHost(0);
  }

  double getDiagonalDeterminant(int dim)
  {
    using ExecSpace = DefaultTestDeviceType::execution_space;
    
    Kokkos::DynRankView<double,DefaultTestDeviceType> detView("detView",1);
    Kokkos::DynRankView<double,DefaultTestDeviceType> jacView("jacView",dim,dim);
    
    Kokkos::RangePolicy<ExecSpace> policy(ExecSpace(), 0, 1);
    
    Kokkos::parallel_for("determinant invocation", policy,
    KOKKOS_LAMBDA(int) {
      detView(0) = diagonalDeterminant(jacView);
    });
    auto detViewHost = getHostCopy(detView);
    return detViewHost(0);
  }

  template<int D>
  void testDiagonalDeterminant(Teuchos::FancyOStream &out, bool &success)
  {
    double relTol = 1.e-14;
    double expectedDeterminant = 1.0;
    for (int d=0; d<D; d++)
    {
      expectedDeterminant *= 2.0;
    }
    
    double det = getDiagonalDeterminant(D);
    TEST_FLOATING_EQUALITY(expectedDeterminant, det, relTol);
    
    det = getDiagonalDeterminant<D>();
    TEST_FLOATING_EQUALITY(expectedDeterminant, det, relTol);
  }

  TEUCHOS_UNIT_TEST( RealSpaceTools, Det1D )
  {
    const int DIMENSION = 1;
    testDiagonalDeterminant<DIMENSION>(out, success);
  }

  TEUCHOS_UNIT_TEST( RealSpaceTools, Det2D )
  {
    const int DIMENSION = 2;
    testDiagonalDeterminant<DIMENSION>(out, success);
  }

  TEUCHOS_UNIT_TEST( RealSpaceTools, Det3D )
  {
    const int DIMENSION = 3;
    testDiagonalDeterminant<DIMENSION>(out, success);
  }
} // namespace

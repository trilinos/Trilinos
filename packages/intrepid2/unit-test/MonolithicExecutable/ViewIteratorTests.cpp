// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   ViewIteratorTests.cpp
    \brief  Tests to verify ViewIterator.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_TestUtils.hpp"
#include "Intrepid2_ViewIterator.hpp"

namespace
{
  using namespace Intrepid2;

  void testIterationCountMatchesEntryCount(Teuchos::FancyOStream &out, bool &success)
  {
    using namespace Intrepid2;
    using Scalar = double;
    
    // note that this test does not involve any access of View data; therefore, it should work fine regardless of the memory space of the View
    
    using DeviceType = DefaultTestDeviceType;
    using ViewIteratorScalar = ViewIterator<ViewType<Scalar,DeviceType>, Scalar>;
    
    // check that the increment operator works to give us the right number of entries
    // we'll use trivial fields so as to factor out problems in the tensor product logic
    int num_fields = 2;
    int num_points = 64;
    ViewType<Scalar,DeviceType> view("view to iterate over",num_fields,num_points);
    ViewIteratorScalar view_iterator(view);
    int entry_count = 0;
    do
    {
      entry_count++;
    } while (view_iterator.increment() >= 0);
    if (entry_count != num_fields * num_points)
    {
      out << "TEST FAILURE: expected to iterate over " << num_fields * num_points << " entries; ";
      out << "instead iterated over " << entry_count << std::endl;
    }
  }
  
  TEUCHOS_UNIT_TEST( ViewIterator, IterationCountMatchesEntryCount )
  {
    testIterationCountMatchesEntryCount(out, success);
  }
  
} // namespace

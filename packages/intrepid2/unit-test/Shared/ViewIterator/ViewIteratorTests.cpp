// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
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
    
    using ViewIteratorScalar = ViewIterator<ViewType<Scalar>, Scalar>;
    
    // check that the increment operator works to give us the right number of entries
    // we'll use trivial fields so as to factor out problems in the tensor product logic
    int num_fields = 2;
    int num_points = 64;
    ViewType<Scalar> view("view to iterate over",num_fields,num_points);
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

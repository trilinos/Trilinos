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

/** \file   TensorBasisTests.cpp
    \brief  Tests to verify the TensorBasis class and free methods defined in Intrepid2_TensorBasis.hpp.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_TensorBasis.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;

  // TensorBasis::getDkCardinality() requires an operator; this one just takes the operator order
  template<int spaceDim>
  int getDkCardinality(const int n)
  {
    return (spaceDim==1) ? 1 :
           (spaceDim==2) ?  n + 1 :
           (spaceDim==3) ? (n + 1) * (n + 2) / 2 :
           (spaceDim==4) ? (n + 1) * (n + 2) * (n + 3) / 6 :
           (spaceDim==5) ? (n + 1) * (n + 2) * (n + 3) * (n + 4) / 24 :
           (spaceDim==6) ? (n + 1) * (n + 2) * (n + 3) * (n + 4) * (n + 5) / 120 :
                           (n + 1) * (n + 2) * (n + 3) * (n + 4) * (n + 5) * (n + 6) / 720;
  }

  template<int spaceDim>
  void testDkEnumeration(Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::Array<int,spaceDim> entries;
    for (int opOrder=1; opOrder <= 10; opOrder++)
    {
      const int dkCardinality = getDkCardinality<spaceDim>(opOrder);
      entries[0] = opOrder;
      for (int d=1; d<spaceDim; d++)
      {
        entries[d] = 0;
      }
      for (int dkEnum=0; dkEnum<dkCardinality; dkEnum++)
      {
        const int dkEnumActual = getDkEnumeration<spaceDim>(entries);
        TEST_EQUALITY(dkEnum, dkEnumActual);
        
        if (dkEnum < dkCardinality - 1)
        {
          // proceed to next entry
          // find the first dimension we should decrement.
          int dimToDecrement = spaceDim - 2;
          int totalOrder = 0;
          for (; dimToDecrement >= 0; dimToDecrement--)
          {
            if (entries[dimToDecrement] == 0) continue; // can't decrement from 0!
            totalOrder = 0; // for dimToDecrement through spaceDim
            for (int d=dimToDecrement; d<spaceDim; d++)
            {
              totalOrder += entries[d];
            }
            if (totalOrder == 0) continue; // try next-lowest dimension, then
            if (totalOrder == entries[dimToDecrement] + entries[spaceDim-1])
            {
              break; // dimToDecrement found
            }
          }
          entries[dimToDecrement]--;
          entries[dimToDecrement+1] = totalOrder - entries[dimToDecrement];
          for (int d=dimToDecrement+2; d<spaceDim; d++)
          {
            entries[d] = 0;
          }
        }
      }
    }
  }

  template<int spaceDim>
  void testDkEnumerationInverse(Teuchos::FancyOStream &out, bool &success)
  {
    for (int opOrder=1; opOrder <= 10; opOrder++)
    {
      const int dkCardinality = getDkCardinality<spaceDim>(opOrder);
      Kokkos::Array<int,spaceDim> entries;
      for (int dkEnum=0; dkEnum<dkCardinality; dkEnum++)
      {
        // fill entries using getDkEnumerationInverse:
        getDkEnumerationInverse<spaceDim>(entries, dkEnum, opOrder);
        const int dkEnumActual = getDkEnumeration<spaceDim>(entries);
        TEST_EQUALITY(dkEnum, dkEnumActual);
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( TensorBasis, DkEnumeration )
  {
    testDkEnumeration<1>(out, success);
    testDkEnumeration<2>(out, success);
    testDkEnumeration<3>(out, success);
    testDkEnumeration<4>(out, success);
    testDkEnumeration<5>(out, success);
    testDkEnumeration<6>(out, success);
    testDkEnumeration<7>(out, success);
  }

  TEUCHOS_UNIT_TEST( TensorBasis, DkEnumerationInverse )
  {
    testDkEnumerationInverse<1>(out, success);
    testDkEnumerationInverse<2>(out, success);
    testDkEnumerationInverse<3>(out, success);
    testDkEnumerationInverse<4>(out, success);
    testDkEnumerationInverse<5>(out, success);
    testDkEnumerationInverse<6>(out, success);
    testDkEnumerationInverse<7>(out, success);
  }
} // namespace

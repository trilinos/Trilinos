/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <stk_util/util/Foreach.hpp>    // for stk_foreach
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "boost/foreach.hpp"            // for auto_any_base
#include "boost/range/iterator_range_core.hpp"  // for operator==, etc
#include "boost/range/iterator_range_io.hpp"  // for operator<<
#include "boost/range/sub_range.hpp"    // for sub_range, operator==



TEST( UnitTestRange, range )
{
  typedef std::vector<int> IntVector;
  typedef boost::sub_range<IntVector> IntRange;

  IntVector array(10,0);

  IntRange a_range = std::make_pair(array.begin(),array.end());
  IntRange b_range(std::make_pair(array.begin(),array.end()));
  IntRange c_range(array.begin(),array.end() );
  IntRange d_range(array);

  EXPECT_EQ( a_range, b_range );
  EXPECT_EQ( a_range, c_range );
  EXPECT_EQ( a_range, d_range );

  stk_foreach( int & i, a_range)
  {
    EXPECT_EQ( i, 0 );
    i = 1;
  }

  stk_foreach( int i, c_range)
  {
    EXPECT_EQ( i, 1 );
  }

  EXPECT_TRUE( a_range == array);

}



/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/


#include <sstream>
#include <stdexcept>
#include <impl/KokkosArray_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

void assert_shapes_are_equal_throw(
  const std::type_info & x_layout ,
  const std::type_info & x_value_type ,
  const unsigned x_rank , const unsigned x_stride ,
  const unsigned x_N0 , const unsigned x_N1 ,
  const unsigned x_N2 , const unsigned x_N3 ,
  const unsigned x_N4 , const unsigned x_N5 ,
  const unsigned x_N6 , const unsigned x_N7 ,

  const std::type_info & y_layout ,
  const std::type_info & y_value_type ,
  const unsigned y_rank , const unsigned y_stride ,
  const unsigned y_N0 , const unsigned y_N1 ,
  const unsigned y_N2 , const unsigned y_N3 ,
  const unsigned y_N4 , const unsigned y_N5 ,
  const unsigned y_N6 , const unsigned y_N7 )
{
  std::ostringstream msg ;

  msg << "KokkosArray::Impl::assert_shape_are_equal_throw( {"
      << " layout(" << x_layout.name()
      << ") value(" << x_value_type.name()
      << ") rank(" << x_rank
      << ") stride(" << x_stride
      << ") dimension("
      << " " << x_N0
      << " " << x_N1
      << " " << x_N2
      << " " << x_N3
      << " " << x_N4
      << " " << x_N5
      << " " << x_N6
      << " " << x_N7
      << " ) } != { "
      << " layout(" << y_layout.name()
      << ") value(" << y_value_type.name()
      << ") rank(" << y_rank
      << ") stride(" << y_stride
      << ") dimension("
      << " " << y_N0
      << " " << y_N1
      << " " << y_N2
      << " " << y_N3
      << " " << y_N4
      << " " << y_N5
      << " " << y_N6
      << " " << y_N7
      << " ) } )" ;

  throw std::runtime_error( msg.str() );
}

}
}


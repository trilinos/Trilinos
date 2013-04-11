// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_Geometry_i_h)
#define Intrepid_MiniTensor_Geometry_i_h


namespace Intrepid {

//
// Helper functions for determining the type of element
//
namespace {

inline
ELEMENT::Type
find_type_1D(Index const nodes)
{
  switch (nodes) {
    case 2:
      return ELEMENT::SEGMENTAL;
      break;
    default:
      break;
  }
  return ELEMENT::UNKNOWN;
}

inline
ELEMENT::Type
find_type_2D(Index const nodes)
{
  switch (nodes) {
    case 3:
      return ELEMENT::TRIANGULAR;
      break;
    case 4:
      return ELEMENT::QUADRILATERAL;
      break;
    default:
      break;
  }
  return ELEMENT::UNKNOWN;
}

inline
ELEMENT::Type
find_type_3D(Index const nodes)
{
  switch (nodes) {
    case 4:
      return ELEMENT::TETRAHEDRAL;
      break;
    case 8:
      return ELEMENT::HEXAHEDRAL;
      break;
    default:
      break;
  }
  return ELEMENT::UNKNOWN;
}

} // anonymous namespace


//
//
//
inline
ELEMENT::Type
find_type(Index const dimension, Index const number_nodes)
{

  ELEMENT::Type
  type = ELEMENT::UNKNOWN;

  switch (dimension) {

    case 1:
      type = find_type_1D(number_nodes);
      break;

    case 2:
      type = find_type_2D(number_nodes);
      break;

    case 3:
      type = find_type_3D(number_nodes);
      break;

    default:
      break;

  }

  if (type == ELEMENT::UNKNOWN) {
    std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
    std::cerr << std::endl;
    std::cerr << "Unknown element type." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Spatial dimension: ";
    std::cerr << dimension << std::endl;
    std::cerr << "Vertices per element: ";
    std::cerr << number_nodes << std::endl;
    exit(1);
  }

  return type;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Geometry_i_h


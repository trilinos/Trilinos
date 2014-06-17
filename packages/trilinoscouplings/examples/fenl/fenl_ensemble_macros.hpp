/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_FENL_ENSEMBLE_MACROS_HPP
#define KOKKOS_EXAMPLE_FENL_ENSEMBLE_MACROS_HPP

#include <fenl_ensemble.hpp>
#include <HexElement.hpp>
#include <fenl_impl.hpp>

namespace Kokkos {
namespace Example {
namespace FENL {

#define INST_DEVICE_SIZE( Device , VectorSize )                         \
  typedef Stokhos::StaticFixedStorage<int,double,VectorSize,Device>     \
    Storage_ ## Device ## _ ## VectorSize;                              \
  typedef Sacado::MP::Vector< Storage_ ## Device ## _ ## VectorSize >   \
    Scalar_ ## Device ## _ ## VectorSize;                               \
  typedef ElementComputationKLCoefficient< Scalar_ ## Device ## _ ## VectorSize,double,Device> \
    KL_ ## Device ## _ ## VectorSize;                                   \
                                                                        \
  INST_FENL( Scalar_ ## Device ## _ ## VectorSize ,                     \
             Device , BoxElemPart::ElemLinear ,                         \
             KL_ ## Device ## _ ## VectorSize ,                         \
             TrivialManufacturedSolution )                              \
                                                                        \
  INST_KL( Scalar_ ## Device ## _ ## VectorSize , double , Device )

#define INST_DEVICE_DOUBLE( Device )                                    \
  typedef ElementComputationKLCoefficient<double,double,Device>         \
    KL_Scalar_ ## Device;                                               \
  INST_FENL( double , Device , BoxElemPart::ElemLinear ,                \
             KL_Scalar_ ## Device , TrivialManufacturedSolution )       \
  INST_FENL( double , Device , BoxElemPart::ElemQuadratic ,             \
             KL_Scalar_ ## Device , TrivialManufacturedSolution )       \
  INST_KL( double , double , Device )

#if defined(__MIC__)
#define INST_DEVICE_HOST( Device ) \
  INST_DEVICE_SIZE( Device, 16 )   \
  INST_DEVICE_SIZE( Device, 32 )   \
  INST_DEVICE_DOUBLE( Device )
#else
#define INST_DEVICE_HOST( Device ) \
  INST_DEVICE_SIZE( Device,  4 )   \
  INST_DEVICE_SIZE( Device, 16 )   \
  INST_DEVICE_SIZE( Device, 32 )   \
  INST_DEVICE_DOUBLE( Device )
#endif

#define INST_DEVICE_GPU( Device )  \
  INST_DEVICE_SIZE( Device, 16 )   \
  INST_DEVICE_SIZE( Device, 32 )   \
  INST_DEVICE_DOUBLE( Device )

}
}
}

#endif /* #ifndef KOKKOS_EXAMPLE_FENL_ENSEMBLE_MACROS_HPP */

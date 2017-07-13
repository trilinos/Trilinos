/*
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER
*/


#ifndef INTREPID2_FIELDCONTAINER_KOKKOS_HPP
#define INTREPID2_FIELDCONTAINER_KOKKOS_HPP

#include "Kokkos_Core.hpp"
#include "Sacado.hpp"
#include <impl/Kokkos_Timer.hpp>

#include <random>
#include <time.h>
#include <stdlib.h>
#include <Kokkos_Random.hpp>
#include "Intrepid2_KokkosRank.hpp"


namespace Intrepid2{
class none{};
class FieldContainer_Kokkos_Ptr;
template <class Scalar,class MemoryLayout=Kokkos::LayoutRight,class ExecutionSpace=Kokkos::DefaultExecutionSpace>
class FieldContainer_Kokkos;



template<class Scalar>
struct initFieldContKokkos{
Scalar* a;
Scalar initValue;
initFieldContKokkos(Scalar initValue_, Scalar* a_): a(a_),initValue(initValue_)
{}
KOKKOS_INLINE_FUNCTION
void operator()(const index_type i)const{
a[i]=initValue;
}

};


template<class FadType, class Layout, class Device, class Scalar>
struct Return_Type< FieldContainer_Kokkos<FadType, Layout, Device>, Scalar>{
      typedef FadType& return_type;
      typedef FadType  const_return_type;
};

template<class FadType, class Layout, class Device, class Scalar>
struct Return_Type<const FieldContainer_Kokkos<FadType, Layout, Device>, Scalar>{
      typedef FadType& return_type;
      typedef FadType  const_return_type;
};

template<class ScalarT, class  Layout,  class Device>
struct CheckType<FieldContainer_Kokkos<ScalarT, Layout, Device> >{
static const bool value = true;
};

template<class ScalarT, class  Layout,  class Device>
struct CheckType<const FieldContainer_Kokkos<ScalarT, Layout, Device> >{
static const bool value = true;
};


}
#include "Intrepid2_FieldContainer_Kokkos_CUDA_Left.hpp"
#include "Intrepid2_FieldContainer_Kokkos_CUDA_Right.hpp"
#include "Intrepid2_FieldContainer_Kokkos_OpenMP_Left.hpp"
#include "Intrepid2_FieldContainer_Kokkos_OpenMP_Right.hpp"
#include "Intrepid2_FieldContainer_Kokkos_PThreads_Left.hpp"
#include "Intrepid2_FieldContainer_Kokkos_PThreads_Right.hpp"
#include "Intrepid2_FieldContainer_Kokkos_Serial_Left.hpp"
#include "Intrepid2_FieldContainer_Kokkos_Serial_Right.hpp"
#endif

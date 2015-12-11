
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

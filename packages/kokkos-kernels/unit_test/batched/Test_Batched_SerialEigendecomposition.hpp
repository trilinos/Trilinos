/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

//#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Eigendecomposition_Decl.hpp"
#include "KokkosBatched_Eigendecomposition_Serial_Impl.hpp"

#include "KokkosKernels_TestUtils.hpp"

using namespace KokkosBatched;

namespace Test {

  template<typename DeviceType,
           typename ViewRank3Type,
           typename ViewRank2Type>
  struct Functor_TestBatchedSerialEigendecomposition {
    ViewRank3Type _A; 
    ViewRank2Type _Er, _Ei;
    ViewRank3Type _UL, _UR;
    ViewRank2Type _W;
    
    KOKKOS_INLINE_FUNCTION
    Functor_TestBatchedSerialEigendecomposition(const ViewRank3Type A,
                                                const ViewRank2Type Er,
                                                const ViewRank2Type Ei,
                                                const ViewRank3Type UL,
                                                const ViewRank3Type UR,
                                                const ViewRank2Type W)
      : _A(A), _Er(Er), _Ei(Ei), _UL(UL), _UR(UR), _W(W)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int k) const {
      auto A  = Kokkos::subview(_A,  k, Kokkos::ALL(), Kokkos::ALL());
      auto er = Kokkos::subview(_Er, k, Kokkos::ALL());
      auto ei = Kokkos::subview(_Ei, k, Kokkos::ALL());
      auto UL = Kokkos::subview(_UL, k, Kokkos::ALL(), Kokkos::ALL());
      auto UR = Kokkos::subview(_UR, k, Kokkos::ALL(), Kokkos::ALL());
      auto W  = Kokkos::subview(_W,  k, Kokkos::ALL());

      SerialEigendecomposition::invoke(A, er, ei, UL, UR, W);        
    }
    
    inline
    void run() {
      typedef typename ViewRank3Type::value_type value_type;
      std::string name_region("KokkosBatched::Test::SerialEigendecomposition");
      std::string name_value_type = ( std::is_same<value_type,float>::value ? "::Float" : 
                                      std::is_same<value_type,double>::value ? "::Double" :
                                      std::is_same<value_type,Kokkos::complex<float> >::value ? "::ComplexFloat" :
                                      std::is_same<value_type,Kokkos::complex<double> >::value ? "::ComplexDouble" : "::UnknownValueType" );                               
      std::string name = name_region + name_value_type;
      Kokkos::Profiling::pushRegion( name.c_str() );
      Kokkos::RangePolicy<DeviceType> policy(0, _A.extent(0));
      Kokkos::parallel_for(name.c_str(), policy, *this);            
      Kokkos::Profiling::popRegion();
    }
  };
    
  template<typename DeviceType,
           typename ValueType,
           typename LayoutType>
  void impl_test_batched_eigendecomposition(const int N, const int m) {
    typedef ValueType value_type;
    typedef Kokkos::View<value_type***,LayoutType,DeviceType> ViewRank3Type;
    typedef Kokkos::View<value_type**,LayoutType,DeviceType> ViewRank2Type;

    /// input
    ViewRank3Type A("A", N, m, m);
    ViewRank2Type W("W", N, 2*m*m+m*5);    

    /// output
    ViewRank2Type Er("Er", N, m);
    ViewRank2Type Ei("Ei", N, m);
    ViewRank3Type UL("UL", N, m, m);
    ViewRank3Type UR("UR", N, m, m);

    
    Kokkos::Random_XorShift64_Pool<typename DeviceType::execution_space> random(13718);
    Kokkos::fill_random(A, random, value_type(1.0));
    Kokkos::fence();

    /// test body
    Functor_TestBatchedSerialEigendecomposition
      <DeviceType,ViewRank3Type,ViewRank2Type>(A, Er, Ei, UL, UR, W).run();
    Kokkos::fence();
  }
}

template<typename DeviceType, 
         typename ValueType>
int test_batched_eigendecomposition() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) 
  {
    Test::impl_test_batched_eigendecomposition<DeviceType,ValueType,Kokkos::LayoutLeft>(0, 10);
    for (int i=0;i<10;++i) {                                                                   
      Test::impl_test_batched_eigendecomposition<DeviceType,ValueType,Kokkos::LayoutLeft>(0, 1);                     
    }
  }
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) 
  {
    Test::impl_test_batched_eigendecomposition<DeviceType,ValueType,Kokkos::LayoutRight>(0, 10);
    for (int i=0;i<10;++i) {                                                                   
      Test::impl_test_batched_eigendecomposition<DeviceType,ValueType,Kokkos::LayoutRight>(0, 1);
    }
  }
#endif
  
  return 0;
}

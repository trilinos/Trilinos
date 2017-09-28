#ifndef __KOKKOSBATCHED_VECTOR_HPP__
#define __KOKKOSBATCHED_VECTOR_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
  namespace Experimental {
    template<typename T>
    class Vector;

    template<typename T, typename SpT, int l>
    struct is_vector<Vector<VectorTag<SIMD<T,SpT>,l> > > {
      static const bool value = true;
    };

    template<typename T, typename SpT, int l>
    struct is_vector<Vector<VectorTag<AVX<T,SpT>,l> > > {
      static const bool value = true;
    };

  }
}

#include "KokkosBatched_Vector_SIMD.hpp"
#include "KokkosBatched_Vector_AVX256D.hpp"
#include "KokkosBatched_Vector_AVX512D.hpp"

namespace Kokkos {
  namespace Details {

    using namespace KokkosBatched::Experimental;

    template<typename T, typename SpT, int l>
    class ArithTraits<Vector<VectorTag<SIMD<T,SpT>,l> > > { 
    public:
      typedef typename ArithTraits<T>::val_type val_type;
      typedef typename ArithTraits<T>::mag_type mag_type;
      
      static const bool is_specialized = ArithTraits<T>::is_specialized;
      static const bool is_signed = ArithTraits<T>::is_signed;
      static const bool is_integer = ArithTraits<T>::is_integer;
      static const bool is_exact = ArithTraits<T>::is_exact;
      static const bool is_complex = ArithTraits<T>::is_complex;
    };
    template<typename T, typename SpT, int l>
    class ArithTraits<Vector<VectorTag<AVX<T,SpT>,l> > > { 
    public:
      typedef typename ArithTraits<T>::val_type val_type;
      typedef typename ArithTraits<T>::mag_type mag_type;
      
      static const bool is_specialized = ArithTraits<T>::is_specialized;
      static const bool is_signed = ArithTraits<T>::is_signed;
      static const bool is_integer = ArithTraits<T>::is_integer;
      static const bool is_exact = ArithTraits<T>::is_exact;
      static const bool is_complex = ArithTraits<T>::is_complex;
    };

  }
}

#endif

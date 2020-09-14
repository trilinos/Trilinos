#ifndef __TACHO_COPY_ON_DEVICE_HPP__
#define __TACHO_COPY_ON_DEVICE_HPP__


/// \file  Tacho_COPY_OnDevice.hpp
/// \brief COPY
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<>
  struct Copy<Algo::OnDevice> {
    template<typename MemberType,
             typename ViewTypeA,
             typename ViewTypeB>
    inline
    static int
    invoke(MemberType &member,
           const ViewTypeA &A,
           const ViewTypeB &B) {      
      const auto exec_instance = member;
      Kokkos::deep_copy(exec_instance, A, B);

      return 0;
    }
  };

}
#endif

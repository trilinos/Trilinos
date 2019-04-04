#ifndef _NgpUnitTestUtils_hpp_
#define _NgpUnitTestUtils_hpp_

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

namespace ngp_unit_test_utils {

template<typename DualViewType>
DualViewType create_dualview(const std::string& name, unsigned size)
{
  DualViewType result(name, size);

  Kokkos::deep_copy(result.h_view, 0);
  result.template modify<typename DualViewType::host_mirror_space>();
  result.template sync<typename DualViewType::execution_space>();

  return result;
}

} // ngp_unit_test_utils

#endif


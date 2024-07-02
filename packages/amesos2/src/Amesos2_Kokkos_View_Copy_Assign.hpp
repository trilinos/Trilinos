// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_Kokkos_View_Copy_Assign.hpp
  \author
  \date   Fri Sept 13 6:00:00 2019

  \brief  Copy or assign views based on memory spaces
*/

#ifndef AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP
#define AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP

namespace Amesos2 {

// allocate dst size if necessary - 2 methods handle 1d and 2d
template<class dst_t, class src_t> // version for 1d view
typename std::enable_if<static_cast<int>(dst_t::rank) == 1>::type
update_dst_size(dst_t & dst, const src_t & src) {
  if(dst.extent(0) != src.extent(0)) { // templated just for 1D
    dst = dst_t(Kokkos::ViewAllocateWithoutInitializing("dst"),
      src.extent(0));
  }
}

template<class dst_t, class src_t> // version for 2d view
typename std::enable_if<static_cast<int>(dst_t::rank) == 2>::type
update_dst_size(dst_t & dst, const src_t & src) {  // templated just for 2d
  if(dst.extent(0) != src.extent(0) || dst.extent(1) != src.extent(1)) {
    dst = dst_t(Kokkos::ViewAllocateWithoutInitializing("dst"),
      src.extent(0), src.extent(1));
  }
}

// now handle type mismatch for same memory space - here types are same
// bInitialize:
//   If bInitialize is false, then the data needs to be allocated but not initialized.
//   If we are about to solve into x we don't care about setting the original values.
//   In this case, we are assigning the view directly so bInitialize does not matter.
// bAssigned:
//   bAssigned tells the caller if the data was simply assigned, so it is set true in this case.
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_same_mem_check_types(
  [[maybe_unused]] bool bInitialize,
  dst_t & dst,
  const src_t & src,
  bool & bAssigned
) {
  dst = src; // just assign the ptr - no need to copy
  bAssigned = true;
}


// deep-copy version (no checking)
// bInitialize:
//   If bInitialize is false, then the data needs to be allocated but not initialized.
//   If we are about to solve into x we don't care about setting the original values.
//   In this case, we are allocating so we first make the memory via update_dst_size.
//   Then we only copy from the source if bInitialize is true.
// bAssigned:
//   bAssigned tells the caller if the data was simply assigned, so it is set false in this case.
template<class dst_t, class src_t> // actual implementation
void deep_copy_only(bool bInitialize, dst_t & dst, const src_t & src, bool & bAssigned) {
  update_dst_size(dst, src); // allocates if necessary
  if(bInitialize) { // bInitialize false would be for solver getting x, where the actual values are not needed
    Kokkos::deep_copy(dst, src); // full copy
  }
  bAssigned = false;
}

template<class dst_t, class src_t> // actual implementation
void deep_copy_only(dst_t & dst, const src_t & src) {
  bool bAssigned;
  deep_copy_only(true, dst, src, bAssigned);
}

// now handle type mismatch for same memory space - now types are different
// bInitialize:
//   If bInitialize is false, then the data needs to be allocated but not initialized.
//   If we are about to solve into x we don't care about setting the original values.
//   In this case, we are allocating so we first make the memory via update_dst_size.
//   Then we only copy from the source if bInitialize is true.
// bAssigned:
//   bAssigned tells the caller if the data was simply assigned, so it is set false in this case.
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<!std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_same_mem_check_types(bool bInitialize, dst_t & dst, const src_t & src, bool & bAssigned) {
  update_dst_size(dst, src); // allocates if necessary
  if(bInitialize) { // bInitialize false would be for solver getting x, where the actual values are not needed
    Kokkos::deep_copy(dst, src); // full copy
  }
  bAssigned = false;
}

// implement_copy_or_assign has 2 versions for matched memory and
// mismatched memory. Right now we just check the memory space.
// a layout mismatch is going to compile fail so probably reflects an error
// in the initial setup.
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
deep_copy_or_assign_view(bool bInitialize, dst_t & dst, const src_t & src, bool & bAssigned) {
  implement_copy_or_assign_same_mem_check_types(bInitialize, dst, src, bAssigned);
}

// for convenience this version does not take bInitialize input and bAssigned ouput
// then it's assumed you want bInitialize true and don't need to know bAssigned
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
deep_copy_or_assign_view(dst_t & dst, const src_t & src) {
  bool bAssigned; // output not needed
  implement_copy_or_assign_same_mem_check_types(true, dst, src, bAssigned);
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_diff_mem_check_types(bool bInitialize, dst_t & dst, const src_t & src, bool & bAssigned) {
  update_dst_size(dst, src); // allocates if necessary
  if(bInitialize) { // bInitialize false would be for solver getting x, where the actual values are not needed
    Kokkos::deep_copy(dst, src); // full copy
  }
  bAssigned = false;
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<static_cast<int>(dst_t::rank) == 1>::type
implement_copy_or_assign_diff_mem_diff_types_check_dim(dst_t & dst, const src_t & src) {
  Kokkos::View<typename dst_t::value_type*, typename src_t::execution_space>
    intermediate(Kokkos::ViewAllocateWithoutInitializing("intermediate"), src.extent(0));
  Kokkos::deep_copy(intermediate, src); // to dst type
  Kokkos::deep_copy(dst, intermediate); // to dst mem
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<static_cast<int>(dst_t::rank) == 2>::type
implement_copy_or_assign_diff_mem_diff_types_check_dim(dst_t & dst, const src_t & src) {
  Kokkos::View<typename dst_t::value_type**, Kokkos::LayoutLeft, typename src_t::execution_space>
    intermediate(Kokkos::ViewAllocateWithoutInitializing("intermediate"), src.extent(0), src.extent(1));
  Kokkos::deep_copy(intermediate, src); // to dst type
  Kokkos::deep_copy(dst, intermediate); // to dst mem
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<!std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_diff_mem_check_types(bool bInitialize, dst_t & dst, const src_t & src, bool & bAssigned) {
  update_dst_size(dst, src); // allocates if necessary
  bAssigned = false;
  if(bInitialize) { // bInitialize false would be for solver getting x, where the actual values are not needed
    // since mem space and types are different, we specify the order of operations
    // Kokkos::deep_copy won't do both since it would be a hidden deep_copy
    implement_copy_or_assign_diff_mem_diff_types_check_dim(dst, src);
  }
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<!std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
deep_copy_or_assign_view(bool bInitialize, dst_t & dst, const src_t & src, bool & bAssigned) {
  implement_copy_or_assign_diff_mem_check_types(bInitialize, dst, src, bAssigned); // full copy
}

// for convenience this version does not take bInitialize input and bAssigned ouput
// then it's assumed you want bInitialize true and don't need to know bAssigned
template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<!std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
deep_copy_or_assign_view(dst_t & dst, const src_t & src) {
  bool bAssigned; // output not needed
  implement_copy_or_assign_diff_mem_check_types(true, dst, src, bAssigned); // full copy
}

} // end namespace Amesos2

#endif  // AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP

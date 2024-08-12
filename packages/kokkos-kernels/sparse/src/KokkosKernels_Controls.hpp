//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef _KOKKOSKERNEL_CONTROLS_HPP
#define _KOKKOSKERNEL_CONTROLS_HPP
/// \file  KokkosKernels_Controls.hpp
/// \brief Mechanism to control internal behavior of kernels
/// \author Luc Berger-Vergiat (lberge@sandia.gov)

#include <string>
#include <unordered_map>
#include <string>
#include "KokkosKernels_config.h"
#include "KokkosKernels_tpl_handles_decl.hpp"

// TPLS headers
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include "cublas_v2.h"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#include <rocsparse/rocsparse.h>
#endif

namespace KokkosKernels {
namespace Experimental {

// Declaration of Controls class
class Controls {
 public:
  using key_type    = std::string;
  using mapped_type = std::string;
  using value_type  = std::pair<const key_type, mapped_type>;

  // Constructor
  Controls() = default;
  Controls(std::initializer_list<value_type> init) : kernel_parameters(init) {}

  // set a new parameter
  void setParameter(const std::string& name, const std::string& value) { kernel_parameters[name] = value; }

  // check if a parameter is already set
  bool isParameter(const std::string& name) const { return kernel_parameters.end() != kernel_parameters.find(name); }

  /// \brief get the value associated with \c name, or \c default if not present
  ///
  /// \param name the name of the parameter to retrieve
  /// \param orUnset (default \c "" ) the value to return if \c name is not set
  key_type getParameter(const std::string& name, const std::string& orUnset = "") const {
    auto search = kernel_parameters.find(name);
    if (kernel_parameters.end() == search) {
      return orUnset;
    } else {
      return search->second;
    }
  }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
  mutable cublasHandle_t cublasHandle = 0;

  cublasHandle_t getCublasHandle() const {
    if (cublasHandle == 0) {
      KokkosBlas::Impl::CudaBlasSingleton& s = KokkosBlas::Impl::CudaBlasSingleton::singleton();
      cublasHandle                           = s.handle;
    }
    return cublasHandle;
  }

  void setCublasHandle(const cublasHandle_t userCublasHandle) { cublasHandle = userCublasHandle; }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  mutable cusparseHandle_t cusparseHandle = 0;

  cusparseHandle_t getCusparseHandle() const {
    if (cusparseHandle == 0) {
      KokkosKernels::Impl::CusparseSingleton& s = KokkosKernels::Impl::CusparseSingleton::singleton();
      cusparseHandle                            = s.cusparseHandle;
    }
    return cusparseHandle;
  }

  void setCusparseHandle(const cusparseHandle_t userCusparseHandle) { cusparseHandle = userCusparseHandle; }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  mutable rocsparse_handle rocsparseHandle = 0;

  rocsparse_handle getRocsparseHandle() const {
    if (rocsparseHandle == 0) {
      KokkosKernels::Impl::RocsparseSingleton& s = KokkosKernels::Impl::RocsparseSingleton::singleton();
      rocsparseHandle                            = s.rocsparseHandle;
    }
    return rocsparseHandle;
  }

  void setRocsparseHandle(const rocsparse_handle userRocsparseHandle) { rocsparseHandle = userRocsparseHandle; }
#endif

 private:
  // storage for kernel parameters
  std::unordered_map<key_type, mapped_type> kernel_parameters;
};

}  // namespace Experimental
}  // namespace KokkosKernels

#endif  // _KOKKOSKERNEL_CONTROLS_HPP

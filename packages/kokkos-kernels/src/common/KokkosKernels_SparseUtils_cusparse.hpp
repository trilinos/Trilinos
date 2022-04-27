/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP
#define _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"

namespace KokkosSparse {
namespace Impl {

inline void cusparse_internal_error_throw(cusparseStatus_t cusparseStatus,
                                          const char* name, const char* file,
                                          const int line) {
  std::ostringstream out;
#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
  out << name << " error( " << cusparseGetErrorName(cusparseStatus)
      << "): " << cusparseGetErrorString(cusparseStatus);
#else
  out << name << " error( ";
  switch (cusparseStatus) {
    case CUSPARSE_STATUS_NOT_INITIALIZED:
      out << "CUSPARSE_STATUS_NOT_INITIALIZED): cusparse handle was not "
             "created correctly.";
      break;
    case CUSPARSE_STATUS_ALLOC_FAILED:
      out << "CUSPARSE_STATUS_ALLOC_FAILED): you might tried to allocate too "
             "much memory";
      break;
    case CUSPARSE_STATUS_INVALID_VALUE:
      out << "CUSPARSE_STATUS_INVALID_VALUE)";
      break;
    case CUSPARSE_STATUS_ARCH_MISMATCH:
      out << "CUSPARSE_STATUS_ARCH_MISMATCH)";
      break;
    case CUSPARSE_STATUS_MAPPING_ERROR:
      out << "CUSPARSE_STATUS_MAPPING_ERROR)";
      break;
    case CUSPARSE_STATUS_EXECUTION_FAILED:
      out << "CUSPARSE_STATUS_EXECUTION_FAILED)";
      break;
    case CUSPARSE_STATUS_INTERNAL_ERROR:
      out << "CUSPARSE_STATUS_INTERNAL_ERROR)";
      break;
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      out << "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED)";
      break;
    case CUSPARSE_STATUS_ZERO_PIVOT:
      out << "CUSPARSE_STATUS_ZERO_PIVOT)";
      break;
    default: out << "unrecognized error code): this is bad!"; break;
  }
#endif  // CUSPARSE_VERSION
  if (file) {
    out << " " << file << ":" << line;
  }
  throw std::runtime_error(out.str());
}

inline void cusparse_internal_safe_call(cusparseStatus_t cusparseStatus,
                                        const char* name,
                                        const char* file = nullptr,
                                        const int line   = 0) {
  if (CUSPARSE_STATUS_SUCCESS != cusparseStatus) {
    cusparse_internal_error_throw(cusparseStatus, name, file, line);
  }
}

// The macro below defines is the public interface for the safe cusparse calls.
// The functions themselves are protected by impl namespace.
#define KOKKOS_CUSPARSE_SAFE_CALL(call)                                  \
  KokkosSparse::Impl::cusparse_internal_safe_call(call, #call, __FILE__, \
                                                  __LINE__)

}  // namespace Impl

}  // namespace KokkosSparse

#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#endif  // _KOKKOSKERNELS_SPARSEUTILS_CUSPARSE_HPP

/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/
#ifndef TPETRA_DETAILS_DEEP_COPY_COUNTER_HPP
#define TPETRA_DETAILS_DEEP_COPY_COUNTER_HPP

/// \file Tpetra_Details_DeepCopyCounter.hpp
/// \brief Declaration of Tpetra::Details::DeepCopyCounter, a class that
///  uses Kokkos' profiling library to count deep copies between memory spaces

#include <cstddef>

namespace Tpetra {
namespace Details {

/// \brief Counter for Kokkos::deep_copy's between memory spaces.
class DeepCopyCounter {
public:
  /// \brief Start the deep_copy counter
  static void start();

  /// \brief Reset the deep_copy counter
  static void reset();

  /// \brief Stop the deep_copy counter
  static size_t stop();

  /// \brief Query the deep_copy counter
  static size_t get_count();


  static size_t count;
  static bool count_active;

};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_DEEP_COPY_COUNTER_HPP

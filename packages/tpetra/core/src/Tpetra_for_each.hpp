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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_FOR_EACH_HPP
#define TPETRA_FOR_EACH_HPP

#include "TpetraCore_config.h"
#include <type_traits>

/// \file Tpetra_for_each.hpp
/// \brief Declaration and definition of Tpetra::for_each;
///   declaration of helper classes for users to specialize.

namespace Tpetra {

  /// \brief Apply a function entrywise to each local entry of a
  ///   Tpetra global data structure, analogously to std::for_each.
  ///
  /// Generically, the function f takes the current entry by nonconst
  /// reference, in read-write fashion.  Overloads for specific Tpetra
  /// data structures may permit functions that take other arguments.
  ///
  /// \param kernelLabel [in] Kernel label for Kokkos Profiling.
  /// \param execSpace [in] Kokkos execution space on which to run.
  /// \param X [in/out] Tpetra global data structure to modify.
  /// \param f [in] Function to apply to each entry of X.
  ///
  /// There is no overload without a kernel label.  You WILL label
  /// your kernels.  If you bother to use Kokkos, that means you care
  /// about performance, so we need to be able to measure it.
  ///
  /// To implement this for new GlobalDataStructure types, specialize
  /// Details::ForEach (see below).
  template<class ExecutionSpace,
           class GlobalDataStructure,
           class UserFunctionType>
  void
  for_each (const char kernelLabel[],
            ExecutionSpace execSpace,
            GlobalDataStructure& X,
            UserFunctionType f);

  /// \brief Overload of for_each (see above) that runs on X's
  ///   default Kokkos execution space.
  ///
  /// \param kernelLabel [in] Kernel label for Kokkos Profiling.
  /// \param X [in/out] Tpetra global data structure to modify.
  /// \param f [in] Function to apply entrywise to X (could have
  ///   different signatures; see above).
  ///
  /// There is no overload without a kernel label.  You WILL label
  /// your kernels.  If you bother to use Kokkos, that means you care
  /// about performance, so we need to be able to measure it.
  ///
  /// To implement this for new GlobalDataStructure types, specialize
  /// Details::ForEach (see below).
  template<class GlobalDataStructure,
           class UserFunctionType>
  void
  for_each (const char kernelLabel[],
            GlobalDataStructure& X,
            UserFunctionType f);

  namespace Details {

    /// \brief Specialize this class to implement Tpetra::for_each for
    ///   specific GlobalDataStructure types.
    template<class ExecutionSpace,
             class GlobalDataStructure,
             class UserFunctionType>
    struct ForEach {
      static void
      for_each (const char kernelLabel[],
                ExecutionSpace execSpace,
                GlobalDataStructure& X,
                UserFunctionType f);
    };

  } // namespace Details

  template<class ExecutionSpace,
           class GlobalDataStructure,
           class UserFunctionType>
  void
  for_each (const char kernelLabel[],
            ExecutionSpace execSpace,
            GlobalDataStructure& X,
            UserFunctionType f)
  {
    using impl_type = Details::ForEach<
      ExecutionSpace,
      typename std::remove_const<GlobalDataStructure>::type,
      UserFunctionType>;

    impl_type::for_each (kernelLabel, execSpace, X, f);
  }

  template<class GlobalDataStructure,
           class UserFunctionType>
  void
  for_each (const char kernelLabel[],
            GlobalDataStructure& X,
            UserFunctionType f)
  {
    using execution_space =
      typename GlobalDataStructure::device_type::execution_space;
    using impl_type = Details::ForEach<
      execution_space,
      typename std::remove_const<GlobalDataStructure>::type,
      UserFunctionType>;

    execution_space execSpace;
    impl_type::for_each (kernelLabel, execSpace, X, f);
  }

} // namespace Tpetra

#endif // TPETRA_FOR_EACH_HPP


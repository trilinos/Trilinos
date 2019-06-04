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

#ifndef TPETRA_TRANSFORM_HPP
#define TPETRA_TRANSFORM_HPP

#include "TpetraCore_config.h"

/// \file Tpetra_transform.hpp
/// \brief Declaration and definition of Tpetra::transform;
///   declaration of helper classes for users to specialize.

namespace Tpetra {

  /// \brief For each local entry input_i of input, assign f(input_i)
  ///   to the corresponding local entry output_i of output,
  ///   analogously to std::transform.
  ///
  /// \param kernelLabel [in] Kernel label for Kokkos Profiling.
  /// \param execSpace [in] Kokkos execution space on which to run.
  /// \param input [in] Tpetra global data structure to read.
  /// \param output [out] Tpetra global data structure to write.
  /// \param f [in] Entrywise function that takes an input object's
  ///   entry by const reference, and returns a value assignable to an
  ///   entry of the output object.
  ///
  /// There is no overload without a kernel label.  You WILL label
  /// your kernels.  If you bother to use Kokkos, that means you care
  /// about performance, so we need to be able to measure it.
  ///
  /// To implement this for new GlobalDataStructure types, specialize
  /// Details::Transform (see below).
  template<class ExecutionSpace,
           class GlobalDataStructure,
           class UnaryFunctionType>
  void
  transform (const char kernelLabel[],
             ExecutionSpace execSpace,
             GlobalDataStructure& input,
             GlobalDataStructure& output,
             UnaryFunctionType f);

  /// \brief Overload of transform (see above) that runs on the output
  ///   object's default Kokkos execution space.
  ///
  /// \param kernelLabel [in] Kernel label for Kokkos Profiling.
  /// \param input [in] Tpetra global data structure to read.
  /// \param output [out] Tpetra global data structure to write.
  /// \param f [in] Entrywise function that takes an input object's
  ///   entry by const reference, and returns a value assignable to an
  ///   entry of the output object.
  ///
  /// There is no overload without a kernel label.  You WILL label
  /// your kernels.  If you bother to use Kokkos, that means you care
  /// about performance, so we need to be able to measure it.
  ///
  /// To implement this for new GlobalDataStructure types, specialize
  /// Details::ForEach (see below).
  template<class GlobalDataStructure,
           class UnaryFunctionType>
  void
  transform (const char kernelLabel[],
             GlobalDataStructure& input,
             GlobalDataStructure& output,
             UnaryFunctionType f);

  namespace Details {

    /// \brief Specialize this class to implement Tpetra::transform
    ///   for specific GlobalDataStructure types.
    template<class ExecutionSpace,
             class GlobalDataStructure>
    struct Transform {
      //! Unary transform: output_i = f(input_i).
      template<class UnaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 GlobalDataStructure& input,
                 GlobalDataStructure& output,
                 UnaryFunctionType f);

      //! Binary transform: output_i = f(input1_i, input2_i).
      template<class BinaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 GlobalDataStructure& input1,
                 GlobalDataStructure& input2,
                 GlobalDataStructure& output,
                 BinaryFunctionType f);
    };

  } // namespace Details

  /// \brief Unary transform on the given execution space instance
  ///   execSpace: output_i = f(input_i).
  template<class ExecutionSpace,
           class GlobalDataStructure,
           class UnaryFunctionType>
  void
  transform (const char kernelLabel[],
             ExecutionSpace execSpace,
             GlobalDataStructure& input,
             GlobalDataStructure& output,
             UnaryFunctionType f)
  {
    using impl_type =
      Details::Transform<ExecutionSpace, GlobalDataStructure>;
    using UFT = UnaryFunctionType;

    impl_type::template transform<UFT> (kernelLabel, execSpace,
                                        input, output, f);
  }

  /// \brief Unary transform on output's default execution space:
  ///   output_i = f(input_i).
  template<class GlobalDataStructure,
           class UnaryFunctionType>
  void
  transform (const char kernelLabel[],
             GlobalDataStructure& input,
             GlobalDataStructure& output,
             UnaryFunctionType f)
  {
    using execution_space =
      typename GlobalDataStructure::device_type::execution_space;
    using impl_type =
      Details::Transform<execution_space, GlobalDataStructure>;
    using UFT = UnaryFunctionType;

    execution_space execSpace;
    impl_type::template transform<UFT> (kernelLabel, execSpace,
                                        input, output, f);
  }

  /// \brief Binary transform on the given execution space execSpace:
  ///   output_i = f(input1_i, input2_i).
  template<class ExecutionSpace,
           class GlobalDataStructure,
           class BinaryFunctionType>
  void
  transform (const char kernelLabel[],
             ExecutionSpace execSpace,
             GlobalDataStructure& input1,
             GlobalDataStructure& input2,
             GlobalDataStructure& output,
             BinaryFunctionType f)
  {
    using impl_type =
      Details::Transform<ExecutionSpace, GlobalDataStructure>;
    using BFT = BinaryFunctionType;

    impl_type::template transform<BFT> (kernelLabel, execSpace,
                                        input1, input2, output, f);
  }

} // namespace Tpetra

#endif // TPETRA_TRANSFORM_HPP

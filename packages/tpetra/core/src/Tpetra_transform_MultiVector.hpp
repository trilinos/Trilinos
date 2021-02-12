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

#ifndef TPETRA_TRANSFORM_MULTIVECTOR_HPP
#define TPETRA_TRANSFORM_MULTIVECTOR_HPP

#include "Tpetra_transform.hpp"
#include "Tpetra_withLocalAccess_MultiVector.hpp"
#include "Tpetra_for_each_MultiVector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_TestForException.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>
#include <sstream>

/// \file Tpetra_transform_MultiVector.hpp
/// \brief Include this file to make Tpetra::transform work with
///   Tpetra::MultiVector and Tpetra::Vector.

namespace Tpetra {
  namespace Details {

    // Kokkos::parallel_for functor that implements unary
    // Tpetra::transform for MultiVector objects.
    //
    // The input and output Views may be the same object (locally).
    template<class InputViewType,
             class OutputViewType,
             class UnaryFunctionType,
             class LocalIndexType>
    class MultiVectorUnaryTransformLoopBody {
    private:
      static_assert (static_cast<int> (InputViewType::Rank) == 2,
                     "InputViewType must be a rank-2 Kokkos::View.");
      static_assert (static_cast<int> (OutputViewType::Rank) == 2,
                     "OutputViewType must be a rank-2 Kokkos::View.");

    public:
      MultiVectorUnaryTransformLoopBody (const InputViewType& in,
                                         const OutputViewType& out,
                                         UnaryFunctionType f) :
        in_ (in), out_ (out), f_ (f)
      {}

      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i, const LocalIndexType j) const {
        out_(i,j) = f_ (in_(i,j));
      }

    private:
      InputViewType in_;
      OutputViewType out_;
      UnaryFunctionType f_;
    };

    // Kokkos::parallel_for functor that implements binary
    // Tpetra::transform for MultiVector objects.
    //
    // The input and output Views may be the same object (locally).
    template<class InputViewType1,
             class InputViewType2,
             class OutputViewType,
             class BinaryFunctionType,
             class LocalIndexType>
    class MultiVectorBinaryTransformLoopBody {
    private:
      static_assert (static_cast<int> (InputViewType1::Rank) == 2,
                     "InputViewType1 must be a rank-2 Kokkos::View.");
      static_assert (static_cast<int> (InputViewType2::Rank) == 2,
                     "InputViewType2 must be a rank-2 Kokkos::View.");
      static_assert (static_cast<int> (OutputViewType::Rank) == 2,
                     "OutputViewType must be a rank-2 Kokkos::View.");

    public:
      MultiVectorBinaryTransformLoopBody (const InputViewType1& in1,
                                          const InputViewType2& in2,
                                          const OutputViewType& out,
                                          BinaryFunctionType f) :
        in1_ (in1), in2_ (in2), out_ (out), f_ (f)
      {}

      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i, const LocalIndexType j) const {
        out_(i,j) = f_ (in1_(i,j), in2_(i,j));
      }

    private:
      InputViewType1 in1_;
      InputViewType2 in2_;
      OutputViewType out_;
      BinaryFunctionType f_;
    };

    /// \brief Implementation of Tpetra::transform for
    ///   Tpetra::MultiVector.
    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT>
    struct Transform<ExecutionSpace,
                     ::Tpetra::MultiVector<SC, LO, GO, NT> >
    {
    private:
      using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
      using IST = typename MV::impl_scalar_type;

      template<class UnaryFunctionType>
      static void
      transform_unary
        (const char kernelLabel[],
         ExecutionSpace execSpace,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& input,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
         UnaryFunctionType f)
      {
        auto input_lcl = input.template getLocalView<ExecutionSpace>(Access::ReadOnly);
        auto output_lcl = output.template getLocalView<ExecutionSpace>(Access::ReadWrite);

        using functor_type = MultiVectorUnaryTransformLoopBody<
          decltype(input_lcl), decltype(output_lcl), UnaryFunctionType, LO>;

        functor_type g (input_lcl, output_lcl, f);

        const LO lclNumRows = input_lcl.extent (0);
        const LO numVecs = input_lcl.extent (1);
        using range_type = Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<2>>;
        range_type range (execSpace, {0, 0}, {lclNumRows, numVecs});

        Kokkos::parallel_for (kernelLabel, range, g);
      }

      template<class BinaryFunctionType>
      static void
      transform_binary
        (const char kernelLabel[],
         ExecutionSpace execSpace,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& input1,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& input2,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
         BinaryFunctionType f)
      {
        // Generic lambdas need C++14, so we need typedefs here.

        auto input1_lcl = input1.template getLocalView<ExecutionSpace>(Access::ReadOnly);
        auto input2_lcl = input2.template getLocalView<ExecutionSpace>(Access::ReadOnly);
        auto output_lcl = output.template getLocalView<ExecutionSpace>(Access::ReadWrite);

        using functor_type = MultiVectorBinaryTransformLoopBody<
          decltype(input1_lcl), decltype(input2_lcl), decltype(output_lcl), BinaryFunctionType, LO>;

        functor_type g (input1_lcl, input2_lcl, output_lcl, f);

        const LO lclNumRows = input1_lcl.extent (0);
        const LO numVecs = input1_lcl.extent (1);
        using range_type = Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<2>>;
        range_type range (execSpace, {0, 0}, {lclNumRows, numVecs});

        Kokkos::parallel_for (kernelLabel, range, g);
      }

    public:
      template<class UnaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 ::Tpetra::MultiVector<SC, LO, GO, NT>& input,
                 ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
                 UnaryFunctionType f)
      {
        using Teuchos::TypeNameTraits;
        using std::endl;

        const int myRank = output.getMap ()->getComm ()->getRank ();
        const bool verbose = ::Tpetra::Details::Behavior::verbose ();
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank << ": Tpetra::transform:" << endl
             << " kernelLabel: " << kernelLabel << endl
             << " ExecutionSpace: "
             << TypeNameTraits<ExecutionSpace>::name () << endl;
          std::cerr << os.str ();
        }

        const size_t numVecs = output.getNumVectors ();
        TEUCHOS_TEST_FOR_EXCEPTION
          (input.getNumVectors () != numVecs, std::invalid_argument,
           "input.getNumVectors() = " << input.getNumVectors () <<
           " != output.getNumVectors() = " << numVecs << ".");

        const bool constStride = output.isConstantStride () &&
          input.isConstantStride ();

        if (!constStride) {
          for (size_t j = 0; j < numVecs; ++j) {
            auto output_j = output.getVectorNonConst (j);
            auto input_j = input.getVectorNonConst (j);

            transform_unary
              (kernelLabel,
               execSpace,
               *input_j,
               *output_j,
               f);
          }
        }
        else {
          transform_unary
            (kernelLabel,
             execSpace,
             input,
             output,
             f);
        }
      }

      // Implementation of binary transform on MultiVectors.
      template<class BinaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 ::Tpetra::MultiVector<SC, LO, GO, NT>& input1,
                 ::Tpetra::MultiVector<SC, LO, GO, NT>& input2,
                 ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
                 BinaryFunctionType f)
      {
        using Teuchos::TypeNameTraits;
        using std::endl;
        const char prefix[] = "Tpetra::transform (binary): ";

        const int myRank = output.getMap ()->getComm ()->getRank ();
        const bool verbose = ::Tpetra::Details::Behavior::verbose ();
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank << ": " << prefix << endl
             << " Tpetra::MultiVector<" << TypeNameTraits<SC>::name ()
             << ", " << TypeNameTraits<LO>::name () << ", "
             << TypeNameTraits<GO>::name () << ", "
             << TypeNameTraits<NT>::name () << ">" << endl
             << " kernelLabel: " << kernelLabel << endl
             << " ExecutionSpace: "
             << TypeNameTraits<ExecutionSpace>::name () << endl;
          std::cerr << os.str ();
        }

        const size_t numVecs = output.getNumVectors ();
        TEUCHOS_TEST_FOR_EXCEPTION
          (input1.getNumVectors () != numVecs, std::invalid_argument,
           prefix << "input1.getNumVectors() = " << input1.getNumVectors ()
           << " != output.getNumVectors() = " << numVecs << ".");
        TEUCHOS_TEST_FOR_EXCEPTION
          (input2.getNumVectors () != numVecs, std::invalid_argument,
           prefix << "input2.getNumVectors() = " << input2.getNumVectors ()
           << " != output.getNumVectors() = " << numVecs << ".");

        const bool constStride = output.isConstantStride () &&
          input1.isConstantStride () && input2.isConstantStride ();

        const LO lclNumRows = static_cast<LO> (output.getLocalLength ());
        using range_type = Kokkos::RangePolicy<ExecutionSpace, LO>;
        range_type range (execSpace, 0, lclNumRows);

        if (! constStride) { // operate on Vectors
          for (size_t j = 0; j < numVecs; ++j) {
            auto output_j = output.getVectorNonConst (j);
            auto input1_j = input1.getVectorNonConst (j);
            auto input2_j = input2.getVectorNonConst (j);

            transform_binary
              (kernelLabel,
               execSpace,
               *input1_j,
               *input2_j,
               *output_j,
               f);
          }
        }
        else { // operate on MultiVectors
          transform_binary
            (kernelLabel,
             execSpace,
             input1,
             input2,
             output,
             f);
        }
      }
    };

    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT>
    struct Transform<ExecutionSpace,
                     ::Tpetra::Vector<SC, LO, GO, NT> >
    {
      using MV = Tpetra::MultiVector<SC, LO, GO, NT>;

      template<class UnaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 ::Tpetra::Vector<SC, LO, GO, NT>& input,
                 ::Tpetra::Vector<SC, LO, GO, NT>& output,
                 UnaryFunctionType f)
      {
        Transform<ExecutionSpace, MV>::transform(
            kernelLabel, execSpace,
            (MV&) input, (MV&) output, f);
      }
    
      template<class BinaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 ::Tpetra::Vector<SC, LO, GO, NT>& input1,
                 ::Tpetra::Vector<SC, LO, GO, NT>& input2,
                 ::Tpetra::Vector<SC, LO, GO, NT>& output,
                 BinaryFunctionType f)
      {
        Transform<ExecutionSpace, MV>::transform(
            kernelLabel, execSpace,
            (MV&) input1, (MV&) input2, (MV&) output, f);
      }
    };
  } // namespace Details

} // namespace Tpetra

#endif // TPETRA_TRANSFORM_MULTIVECTOR_HPP

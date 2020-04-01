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
      operator () (const LocalIndexType i) const {
        using LO = LocalIndexType;
        const LO numCols = static_cast<LO> (in_.extent (1));
        for (LO j = 0; j < numCols; ++j) {
          out_(i,j) = f_ (in_(i,j));
        }
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
      operator () (const LocalIndexType i) const {
        using LO = LocalIndexType;
        const LO numCols = static_cast<LO> (in1_.extent (1));
        for (LO j = 0; j < numCols; ++j) {
          out_(i,j) = f_ (in1_(i,j), in2_(i,j));
        }
      }

    private:
      InputViewType1 in1_;
      InputViewType2 in2_;
      OutputViewType out_;
      BinaryFunctionType f_;
    };

    // Kokkos::parallel_for functor that implements unary
    // Tpetra::transform for Vector objects.
    //
    // The input and output Views may be the same object (locally).
    template<class InputViewType,
             class OutputViewType,
             class UnaryFunctionType,
             class LocalIndexType>
    class VectorUnaryTransformLoopBody {
    private:
      static_assert (static_cast<int> (InputViewType::Rank) == 1,
                     "InputViewType must be a rank-1 Kokkos::View.");
      static_assert (static_cast<int> (OutputViewType::Rank) == 1,
                     "OutputViewType must be a rank-1 Kokkos::View.");

    public:
      VectorUnaryTransformLoopBody (const InputViewType& in,
                                    const OutputViewType& out,
                                    UnaryFunctionType f) :
        in_ (in), out_ (out), f_ (f)
      {}

      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i) const {
        out_(i) = f_ (in_(i));
      }

    private:
      InputViewType in_;
      OutputViewType out_;
      UnaryFunctionType f_;
    };

    // Kokkos::parallel_for functor that implements binary
    // Tpetra::transform for Vector objects.
    //
    // The input and output Views may be the same object (locally).
    template<class InputViewType1,
             class InputViewType2,
             class OutputViewType,
             class BinaryFunctionType,
             class LocalIndexType>
    class VectorBinaryTransformLoopBody {
    private:
      static_assert (static_cast<int> (InputViewType1::Rank) == 1,
                     "InputViewType1 must be a rank-1 Kokkos::View.");
      static_assert (static_cast<int> (InputViewType1::Rank) == 1,
                     "InputViewType1 must be a rank-1 Kokkos::View.");
      static_assert (static_cast<int> (OutputViewType::Rank) == 1,
                     "OutputViewType must be a rank-1 Kokkos::View.");

    public:
      VectorBinaryTransformLoopBody (const InputViewType1& in1,
                                     const InputViewType2& in2,
                                     const OutputViewType& out,
                                     BinaryFunctionType f) :
        in1_ (in1), in2_ (in2), out_ (out), f_ (f)
      {}

      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i) const {
        out_(i) = f_ (in1_(i), in2_(i));
      }

    private:
      InputViewType1 in1_;
      InputViewType2 in2_;
      OutputViewType out_;
      BinaryFunctionType f_;
    };

    // CUDA 9.2 doesn't like it when you call lambdas in private or
    // protected methods (of Transform, in this case).  Thus, I've
    // broken out Transform::transform_sameObject into a separate
    // functor and nonmember function.
    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT,
             class UnaryFunctionType>
    class UnaryTransformSameMultiVector {
    private:
      using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
      using IST = typename MV::impl_scalar_type;

    public:
      UnaryTransformSameMultiVector (UnaryFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION void operator() (IST& X_ij) const {
        // User function has the form IST(const IST&) suitable for
        // transform, so we have to convert it to a for_each function
        // of the form void(IST&).
        X_ij = f_(X_ij);
      }

    private:
      UnaryFunctionType f_;
    };

    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT,
             class UnaryFunctionType>
    void
    unaryTransformSameMultiVector (const char kernelLabel[],
                                   ExecutionSpace execSpace,
                                   ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
                                   UnaryFunctionType f)
    {
      using functor_type = UnaryTransformSameMultiVector<ExecutionSpace,
        SC, LO, GO, NT, UnaryFunctionType>;
      ::Tpetra::for_each (kernelLabel, execSpace, output, functor_type (f));
    }

    /// \brief Implementation of Tpetra::transform for
    ///   Tpetra::MultiVector.
    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT>
    struct Transform<ExecutionSpace,
                     ::Tpetra::MultiVector<SC, LO, GO, NT> >
    {
    private:
      // Given a Kokkos execution space on which the user wants to run
      // the transform, and a memory space in which the MultiVector's
      // data live, determine the memory space that transform should
      // use in its withLocalAccess call.
      template<class MemorySpace>
      using transform_memory_space =
        typename std::conditional<
          Kokkos::SpaceAccessibility<
            ExecutionSpace,
            typename MemorySpace::memory_space>::accessible,
          typename MemorySpace::memory_space,
          typename ExecutionSpace::memory_space>::type;

      using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
      using preferred_memory_space =
        typename MV::device_type::memory_space;
      using memory_space =
        transform_memory_space<preferred_memory_space>;
      using IST = typename MV::impl_scalar_type;

      // This is not the same as "aliases" -- we actually want to know
      // if input and output are the same object (locally), so that we
      // can sync correctly.  The result of transform is undefined if
      // input and output partially alias one another.
      static bool
      sameObject (const ::Tpetra::MultiVector<SC, LO, GO, NT>& input,
                  const ::Tpetra::MultiVector<SC, LO, GO, NT>& output)
      {
        return &input == &output ||
          input.getLocalViewHost ().data () ==
          output.getLocalViewHost ().data () ||
          input.getLocalViewDevice ().data () ==
          output.getLocalViewDevice ().data ();
      }

      template<class UnaryFunctionType>
      static void
      transform_vec_notSameObject
        (const char kernelLabel[],
         ExecutionSpace execSpace,
         ::Tpetra::Vector<SC, LO, GO, NT>& input,
         ::Tpetra::Vector<SC, LO, GO, NT>& output,
         UnaryFunctionType f)
      {
        memory_space memSpace;
        // Generic lambdas need C++14, so we need a typedef here.
        using input_view_type =
          with_local_access_function_argument_type<
            decltype (readOnly (input).on (memSpace). at(execSpace))>;
        using output_view_type =
          with_local_access_function_argument_type<
            decltype (writeOnly (output).on (memSpace). at(execSpace))>;

        withLocalAccess
          ([=] (const input_view_type& input_lcl,
                const output_view_type& output_lcl) {
             using functor_type = VectorUnaryTransformLoopBody<
               input_view_type, output_view_type, UnaryFunctionType, LO>;
             functor_type g (input_lcl, output_lcl, f);

             const LO lclNumRows = static_cast<LO> (input_lcl.extent (0));
             using range_type = Kokkos::RangePolicy<ExecutionSpace, LO>;
             range_type range (execSpace, 0, lclNumRows);

             Kokkos::parallel_for (kernelLabel, range, g);
           },
           readOnly (input).on (memSpace).at (execSpace),
           writeOnly (output).on (memSpace).at (execSpace));
      }

      template<class UnaryFunctionType>
      static void
      transform_mv_notSameObject
        (const char kernelLabel[],
         ExecutionSpace execSpace,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& input,
         ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
         UnaryFunctionType f)
      {
        memory_space memSpace;
        // Generic lambdas need C++14, so we need typedefs here.
        using input_view_type =
          with_local_access_function_argument_type<
            decltype (readOnly (input).on (memSpace). at(execSpace))>;
        using output_view_type =
          with_local_access_function_argument_type<
            decltype (writeOnly (output).on (memSpace). at(execSpace))>;

        withLocalAccess
          ([=] (const input_view_type& input_lcl,
                const output_view_type& output_lcl) {
            using functor_type = MultiVectorUnaryTransformLoopBody<
            input_view_type, output_view_type, UnaryFunctionType, LO>;
            functor_type g (input_lcl, output_lcl, f);

            const LO lclNumRows = static_cast<LO> (input_lcl.extent (0));
            using range_type = Kokkos::RangePolicy<ExecutionSpace, LO>;
            range_type range (execSpace, 0, lclNumRows);

            Kokkos::parallel_for (kernelLabel, range, g);
          },
          readOnly (input).on (memSpace).at (execSpace),
          writeOnly (output).on (memSpace).at (execSpace));
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

        memory_space memSpace;
        if (numVecs == size_t (1) || ! constStride) {
          for (size_t j = 0; j < numVecs; ++j) {
            auto output_j = output.getVectorNonConst (j);
            auto input_j = input.getVectorNonConst (j);

            // Check for aliasing here, since it's possible for only
            // some columns of input & output to alias.  Aliasing is a
            // correctness issue (e.g., for sync'ing).
            if (sameObject (*output_j, *input_j)) {
              unaryTransformSameMultiVector (kernelLabel, execSpace,
                                             *output_j, f);
            }
            else {
              transform_vec_notSameObject (kernelLabel, execSpace,
                                           *input_j, *output_j, f);
            }
          }
        }
        else {
          if (sameObject (output, input)) {
            unaryTransformSameMultiVector (kernelLabel, execSpace,
                                           output, f);
          }
          else {
            transform_mv_notSameObject (kernelLabel, execSpace,
                                        input, output, f);
          }
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
        memory_space memSpace;

        const LO lclNumRows = static_cast<LO> (output.getLocalLength ());
        using range_type = Kokkos::RangePolicy<ExecutionSpace, LO>;
        range_type range (execSpace, 0, lclNumRows);

        if (numVecs == size_t (1) || ! constStride) { // operate on Vectors
          for (size_t j = 0; j < numVecs; ++j) {
            auto output_j = output.getVectorNonConst (j);
            auto input1_j = input1.getVectorNonConst (j);
            auto input2_j = input2.getVectorNonConst (j);

            // Check for aliasing here, since it's possible for only
            // some columns of input & output to alias.  Aliasing is a
            // correctness issue (e.g., for sync'ing).
            const bool outin1same = sameObject (*output_j, *input1_j);
            const bool outin2same = sameObject (*output_j, *input2_j);
            // Don't double-view input1.
            const bool in1in2same = sameObject (*input1_j, *input2_j);
            const bool allsame = outin1same && outin2same; // by transitivity

            // Help GCC 4.9.3 deduce the types of *output_j,
            // *input1_j, and *input2_j.  See discussion here:
            // https://github.com/trilinos/Trilinos/pull/5115
            using vec_type = ::Tpetra::Vector<SC, LO, GO, NT>;
            vec_type& input1_j_ref = *input1_j;
            vec_type& input2_j_ref = *input2_j;
            vec_type& output_j_ref = *output_j;

            // Once we get C++14 generic lambdas, we can get rid of
            // these typedefs and use "const auto&" as the argument(s)
            // for the withLocalAccess lambdas below.
            using input1_view_type =
              with_local_access_function_argument_type<
                decltype (readOnly (input1_j_ref).on (memSpace). at(execSpace))>;
            using input2_view_type =
              with_local_access_function_argument_type<
                decltype (readOnly (input2_j_ref).on (memSpace). at(execSpace))>;
            using rw_output_view_type =
              with_local_access_function_argument_type<
                decltype (readWrite (output_j_ref).on (memSpace). at(execSpace))>;
            using wo_output_view_type =
              with_local_access_function_argument_type<
                decltype (writeOnly (output_j_ref).on (memSpace). at(execSpace))>;

            if (allsame) {
              withLocalAccess
                ([=] (const rw_output_view_type& output_lcl) {
                  using functor_type = VectorBinaryTransformLoopBody<
                    typename rw_output_view_type::const_type,
                    typename rw_output_view_type::const_type,
                    rw_output_view_type,
                    BinaryFunctionType, LO>;
                  functor_type functor (output_lcl, output_lcl, output_lcl, f);
                  Kokkos::parallel_for (kernelLabel, range, functor);
                },
                readWrite (output_j_ref).on (memSpace).at (execSpace));
            }
            else if (in1in2same) { // and not same as output
              withLocalAccess
                ([=] (const input1_view_type& input1_lcl,
                      const wo_output_view_type& output_lcl) {
                  using functor_type = VectorBinaryTransformLoopBody<
                    input1_view_type,
                    input1_view_type,
                    wo_output_view_type,
                    BinaryFunctionType, LO>;
                  functor_type functor (input1_lcl, input1_lcl, output_lcl, f);
                  Kokkos::parallel_for (kernelLabel, range, functor);
                },
                readOnly (input1_j_ref).on (memSpace).at (execSpace),
                writeOnly (output_j_ref).on (memSpace).at (execSpace));
            }
            else if (outin1same) { // and input1 not same as input2
              withLocalAccess
                ([=] (const input2_view_type& input2_lcl,
                      const rw_output_view_type& output_lcl) {
                  using functor_type = VectorBinaryTransformLoopBody<
                    typename rw_output_view_type::const_type,
                    input2_view_type,
                    rw_output_view_type,
                    BinaryFunctionType, LO>;
                  functor_type functor (output_lcl, input2_lcl, output_lcl, f);
                  Kokkos::parallel_for (kernelLabel, range, functor);
                },
                readOnly (input2_j_ref).on (memSpace).at (execSpace),
                readWrite (output_j_ref).on (memSpace).at (execSpace));
            }
            else if (outin2same) { // and input1 not same as input2
              withLocalAccess
                ([=] (const input1_view_type& input1_lcl,
                      const rw_output_view_type& output_lcl) {
                  using functor_type = VectorBinaryTransformLoopBody<
                    input1_view_type,
                    typename rw_output_view_type::const_type,
                    rw_output_view_type,
                    BinaryFunctionType, LO>;
                  functor_type functor (input1_lcl, output_lcl, output_lcl, f);
                  Kokkos::parallel_for (kernelLabel, range, functor);
                },
                readOnly (input1_j_ref).on (memSpace).at (execSpace),
                readWrite (output_j_ref).on (memSpace).at (execSpace));
            }
            else { // output, input1, and input2 all differ
              withLocalAccess
                ([=] (const input1_view_type& input1_lcl,
                      const input2_view_type& input2_lcl,
                      const wo_output_view_type& output_lcl) {
                  using functor_type = VectorBinaryTransformLoopBody<
                    input1_view_type,
                    input2_view_type,
                    wo_output_view_type,
                    BinaryFunctionType, LO>;
                  functor_type functor (input1_lcl, input2_lcl, output_lcl, f);
                  Kokkos::parallel_for (kernelLabel, range, functor);
                },
                readOnly (input1_j_ref).on (memSpace).at (execSpace),
                readOnly (input2_j_ref).on (memSpace).at (execSpace),
                writeOnly (output_j_ref).on (memSpace).at (execSpace));
            }
          }
        }
        else { // operate on MultiVectors
          // Check for aliasing here, since it's possible for only
          // some columns of input & output to alias.  Aliasing is a
          // correctness issue (e.g., for sync'ing).
          const bool outin1same = sameObject (output, input1);
          const bool outin2same = sameObject (output, input2);
          // Don't double-view input1.
          const bool in1in2same = sameObject (input1, input2);
          const bool allsame = outin1same && outin2same; // by transitivity

          // Once we get C++14 generic lambdas, we can get rid of
          // these typedefs and use "const auto&" as the argument(s)
          // for the withLocalAccess lambdas below.
          using input1_view_type =
            with_local_access_function_argument_type<
              decltype (readOnly (input1).on (memSpace). at(execSpace))>;
          using input2_view_type =
            with_local_access_function_argument_type<
              decltype (readOnly (input2).on (memSpace). at(execSpace))>;
          using rw_output_view_type =
            with_local_access_function_argument_type<
              decltype (readWrite (output).on (memSpace). at(execSpace))>;
          using wo_output_view_type =
            with_local_access_function_argument_type<
              decltype (writeOnly (output).on (memSpace). at(execSpace))>;

          if (allsame) {
            withLocalAccess
              ([=] (const rw_output_view_type& output_lcl) {
                using functor_type = MultiVectorBinaryTransformLoopBody<
                  typename rw_output_view_type::const_type,
                  typename rw_output_view_type::const_type,
                  rw_output_view_type,
                  BinaryFunctionType, LO>;
                functor_type functor (output_lcl, output_lcl, output_lcl, f);
                Kokkos::parallel_for (kernelLabel, range, functor);
              },
              readWrite (output).on (memSpace).at (execSpace));
          }
          else if (in1in2same) { // and not same as output
            withLocalAccess
              ([=] (const input1_view_type& input1_lcl,
                    const wo_output_view_type& output_lcl) {
                using functor_type = MultiVectorBinaryTransformLoopBody<
                  input1_view_type,
                  input1_view_type,
                  wo_output_view_type,
                  BinaryFunctionType, LO>;
                functor_type functor (input1_lcl, input1_lcl, output_lcl, f);
                Kokkos::parallel_for (kernelLabel, range, functor);
              },
              readOnly (input1).on (memSpace).at (execSpace),
              writeOnly (output).on (memSpace).at (execSpace));
          }
          else if (outin1same) { // and input1 not same as input2
            withLocalAccess
              ([=] (const input2_view_type& input2_lcl,
                    const rw_output_view_type& output_lcl) {
                using functor_type = MultiVectorBinaryTransformLoopBody<
                  typename rw_output_view_type::const_type,
                  input2_view_type,
                  rw_output_view_type,
                  BinaryFunctionType, LO>;
                functor_type functor (output_lcl, input2_lcl, output_lcl, f);
                Kokkos::parallel_for (kernelLabel, range, functor);
              },
              readOnly (input2).on (memSpace).at (execSpace),
              readWrite (output).on (memSpace).at (execSpace));
          }
          else if (outin2same) { // and input1 not same as input2
            withLocalAccess
              ([=] (const input1_view_type& input1_lcl,
                    const rw_output_view_type& output_lcl) {
                using functor_type = MultiVectorBinaryTransformLoopBody<
                  input1_view_type,
                  typename rw_output_view_type::const_type,
                  rw_output_view_type,
                  BinaryFunctionType, LO>;
                functor_type functor (input1_lcl, output_lcl, output_lcl, f);
                Kokkos::parallel_for (kernelLabel, range, functor);
              },
              readOnly (input1).on (memSpace).at (execSpace),
              readWrite (output).on (memSpace).at (execSpace));
          }
          else { // output, input1, and input2 all differ
            withLocalAccess
              ([=] (const input1_view_type& input1_lcl,
                    const input2_view_type& input2_lcl,
                    const wo_output_view_type& output_lcl) {
                using functor_type = MultiVectorBinaryTransformLoopBody<
                  input1_view_type,
                  input2_view_type,
                  wo_output_view_type,
                  BinaryFunctionType, LO>;
                functor_type functor (input1_lcl, input2_lcl, output_lcl, f);
                Kokkos::parallel_for (kernelLabel, range, functor);
              },
              readOnly (input1).on (memSpace).at (execSpace),
              readOnly (input2).on (memSpace).at (execSpace),
              writeOnly (output).on (memSpace).at (execSpace));
          }
        }
      }
    };

    /// \brief Implementation of Tpetra::transform for Tpetra::Vector.
    ///
    /// Even though Tpetra::Vector is a subclass of
    /// Tpetra::MultiVector, this needs to exist, since partial
    /// specializations don't recognize subclasses of their
    /// specialized template arguments.
    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT>
    struct Transform<ExecutionSpace,
                     ::Tpetra::Vector<SC, LO, GO, NT> >
    {
      // Implementation of unary transform on Vectors.
      template<class UnaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 ::Tpetra::Vector<SC, LO, GO, NT>& input,
                 ::Tpetra::Vector<SC, LO, GO, NT>& output,
                 UnaryFunctionType f)
      {
        using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
        using impl_type = Transform<ExecutionSpace, MV>;
        using UFT = UnaryFunctionType;

        impl_type::template transform<UFT> (kernelLabel, execSpace,
                                            input, output, f);
      }

      // Implementation of binary transform on Vectors.
      template<class BinaryFunctionType>
      static void
      transform (const char kernelLabel[],
                 ExecutionSpace execSpace,
                 ::Tpetra::Vector<SC, LO, GO, NT>& input1,
                 ::Tpetra::Vector<SC, LO, GO, NT>& input2,
                 ::Tpetra::Vector<SC, LO, GO, NT>& output,
                 BinaryFunctionType f)
      {
        using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
        using impl_type = Transform<ExecutionSpace, MV>;
        using BFT = BinaryFunctionType;

        impl_type::template transform<BFT> (kernelLabel, execSpace,
                                            input1, input2, output, f);
      }
    };

  } // namespace Details

} // namespace Tpetra

#endif // TPETRA_TRANSFORM_MULTIVECTOR_HPP

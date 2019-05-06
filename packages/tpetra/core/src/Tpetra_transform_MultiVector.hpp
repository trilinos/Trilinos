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

    // We implement unary transform as a for_each, where
    // InnerLoopBodyType (a function f(X_ij,i,j)) just does X_ij =
    // g(Y(i,j)), for the user's unary function g and input View Y.
    // That is, Tpetra::transform will pass
    // MultiVectorUnaryTransformLoopBody(Y,g) as the "user function"
    // to Tpetra::for_each, after using withLocalAccess to get Y.
    //
    // The input View Y and output View X may be the same object
    // (locally).  Tpetra::transform must check for this case before
    // invoking Tpetra::for_each, so that it doesn't incorrectly sync.
    //
    // This suggests a general approach, based on expression
    // templates, for chaining arithmetic expressions involving
    // Tpetra::MultiVector or Tpetra::Vector.  We won't take the
    // general approach for now, but see
    // MultiVectorBinaryTransformLoopBody.
    //
    // OutputScalarRefType is normally impl_scalar_type&, but
    // could be a proxy reference type -- anything that
    // Kokkos::View<T>::operator() returns for nonconst T.

    template<class InputViewType,
             class OutputScalarRefType,
             class InnerLoopBodyType,
             class LocalIndexType>
    struct MultiVectorUnaryTransformLoopBody {
      static_assert (static_cast<int> (InputViewType::Rank) == 2,
                     "InputViewType must be a rank-2 Kokkos::View.");

      MultiVectorUnaryTransformLoopBody (const InputViewType& Y,
                                         InnerLoopBodyType g) :
        Y_ (Y), g_ (g)
      {}
      KOKKOS_INLINE_FUNCTION void
      operator () (OutputScalarRefType X_ij,
                   const LocalIndexType i,
                   const LocalIndexType j) const
      {
        X_ij = g_ (Y_(i,j));
      };
      InputViewType Y_;
      InnerLoopBodyType g_;
    };

    template<class InputViewType,
             class OutputScalarRefType,
             class InnerLoopBodyType,
             class LocalIndexType>
    struct VectorUnaryTransformLoopBody {
      static_assert (static_cast<int> (InputViewType::Rank) == 1,
                     "InputViewType must be a rank-1 Kokkos::View.");

      VectorUnaryTransformLoopBody (const InputViewType& Y,
                                    InnerLoopBodyType g) :
        Y_ (Y), g_ (g)
      {}
      KOKKOS_INLINE_FUNCTION void
      operator () (OutputScalarRefType X_ij,
                   const LocalIndexType i,
                   const LocalIndexType /* j */) const
      {
        X_ij = g_ (Y_(i));
      };
      InputViewType Y_;
      InnerLoopBodyType g_;
    };

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
      // can defer to for_each in that case.  The result of transform
      // is undefined if input and output partially alias one another.
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
            decltype (readOnly (input).on (memSpace))>;
        // Lambda needs to be mutable, else it makes output const.
        withLocalAccess
          ([=] (const input_view_type& input_lcl) mutable {
             using output_scalar_ref_type = IST&;
             using loop_body_type = VectorUnaryTransformLoopBody<
               input_view_type, IST&, UnaryFunctionType, LO>;

             loop_body_type g (input_lcl, f);
             ::Tpetra::for_each (kernelLabel, execSpace, output, g);
           }, readOnly (input).on (memSpace));
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
        // Generic lambdas need C++14, so we need a typedef here.
        using input_view_type =
          with_local_access_function_argument_type<
            decltype (readOnly (input).on (memSpace))>;
        // Lambda needs to be mutable, else it makes output const.
        withLocalAccess
          ([=] (const input_view_type& input_lcl) mutable {
             using output_scalar_ref_type = IST&;
             using loop_body_type = MultiVectorUnaryTransformLoopBody<
               input_view_type, IST&, UnaryFunctionType, LO>;
             loop_body_type g (input_lcl, f);

             ::Tpetra::for_each (kernelLabel, execSpace, output, g);
           }, readOnly (input).on (memSpace));
      }

      template<class UnaryFunctionType>
      static void
      transform_sameObject (const char kernelLabel[],
                            ExecutionSpace execSpace,
                            ::Tpetra::MultiVector<SC, LO, GO, NT>& output,
                            UnaryFunctionType f)
      {
        // User function has the form IST(const IST&) suitable for
        // transform, so we have to convert it to a for_each function.
        ::Tpetra::for_each
            (kernelLabel, execSpace, output,
             KOKKOS_LAMBDA (IST& X_ij) { X_ij = f(X_ij); });
      }

    public:
      // Implementation of unary Transform on
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
              transform_sameObject (kernelLabel, execSpace, *output_j, f);
            }
            else {
              transform_vec_notSameObject (kernelLabel, execSpace,
                                           *input_j, *output_j, f);
            }
          }
        }
        else {
          if (sameObject (output, input)) {
            transform_sameObject (kernelLabel, execSpace, output, f);
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

        if (numVecs == size_t (1) || ! constStride) { // operate on Vectors
          for (size_t j = 0; j < numVecs; ++j) {
            auto output_j = output.getVectorNonConst (j);
            auto input1_j = input1.getVectorNonConst (j);
            auto input2_j = input2.getVectorNonConst (j);

            // Generic lambdas need C++14, so we need typedefs here.
            using input1_view_type =
              with_local_access_function_argument_type<
                decltype (readOnly (*input1_j).on (memSpace))>;
            using input2_view_type =
              with_local_access_function_argument_type<
                decltype (readOnly (*input2_j).on (memSpace))>;

            // Check for aliasing here, since it's possible for only
            // some columns of input & output to alias.  Aliasing is a
            // correctness issue (e.g., for sync'ing).
            const bool outin1same = sameObject (*output_j, *input1_j);
            const bool outin2same = sameObject (*output_j, *input2_j);
            // Don't double-view input1.
            const bool in1in2same = sameObject (*input1_j, *input2_j);
            const bool allsame = outin1same && outin2same; // by transitivity

            // If input1 or input2 is the same as output, then we can
            // rely on Tpetra::for_each accessing output as readWrite.
            // If either input1 or input2 differs from output, we can
            // use withLocalAccess with readOnly to access it.  The
            // ensuing lambda needs to be mutable, else it makes
            // output const.

            if (allsame) {
              ::Tpetra::for_each
                (kernelLabel, execSpace, *output_j,
                 KOKKOS_LAMBDA (IST& out_i) {
                  out_i = f (out_i, out_i);
                });
            }
            else if (outin1same) {
              withLocalAccess
                ([=] (const input2_view_type& input2_lcl) mutable {
                   ::Tpetra::for_each
                       (kernelLabel, execSpace, *output_j,
                        KOKKOS_LAMBDA (IST& out, const LO i) {
                         out = f (out, input2_lcl(i));
                       });
                 }, readOnly (*input2_j).on (memSpace));
            }
            else if (outin2same) {
              withLocalAccess
                ([=] (const input1_view_type& input1_lcl) mutable {
                   ::Tpetra::for_each
                       (kernelLabel, execSpace, *output_j,
                        KOKKOS_LAMBDA (IST& out, const LO i) {
                         out = f (input1_lcl(i), out);
                       });
                 }, readOnly (*input1_j).on (memSpace));
            }
            else if (in1in2same) {
              withLocalAccess
                ([=] (const input1_view_type& input1_lcl) mutable {
                   ::Tpetra::for_each
                    (kernelLabel, execSpace, *output_j,
                     KOKKOS_LAMBDA (IST& out, const LO i) {
                      out = f (input1_lcl(i), input1_lcl(i));
                    });
                 }, readOnly (*input1_j).on (memSpace));
            }
            else {
              withLocalAccess
                ([=] (const input1_view_type& input1_lcl) mutable {
                   withLocalAccess
                     ([=] (const input2_view_type& input2_lcl) mutable {
                        ::Tpetra::for_each
                         (kernelLabel, execSpace, *output_j,
                          KOKKOS_LAMBDA (IST& out, const LO i) {
                           out = f (input1_lcl(i), input2_lcl(i));
                         });
                      }, readOnly (*input2_j).on (memSpace));
                 }, readOnly (*input1_j).on (memSpace));
              
              // withLocalAccess
              //   ([=] (const input1_view_type& input1_lcl,
              //         const input2_view_type& input2_lcl) mutable {
              //      ::Tpetra::for_each
              //       (kernelLabel, execSpace, *output_j,
              //        KOKKOS_LAMBDA (IST& out, const LO i) {
              //         out = f (input1_lcl(i), input2_lcl(i));
              //       });
              //    },
              //    readOnly (*input1_j).on (memSpace),
              //    readOnly (*input2_j).on (memSpace));
            }
          }
        }
        else { // operate on MultiVectors
          // Generic lambdas need C++14, so we need typedefs here.
          using input1_view_type =
            with_local_access_function_argument_type<
              decltype (readOnly (input1).on (memSpace))>;
          using input2_view_type =
            with_local_access_function_argument_type<
              decltype (readOnly (input2).on (memSpace))>;

          // Check for aliasing here, since it's possible for only
          // some columns of input & output to alias.  Aliasing is a
          // correctness issue (e.g., for sync'ing).
          const bool outin1same = sameObject (output, input1);
          const bool outin2same = sameObject (output, input2);
          // Don't double-view input1.
          const bool in1in2same = sameObject (input1, input2);
          const bool allsame = outin1same && outin2same; // by transitivity

          // If input1 or input2 is the same as output, then we can
          // rely on Tpetra::for_each accessing output as readWrite.
          // If either input1 or input2 differs from output, we can
          // use withLocalAccess with readOnly to access it.  The
          // ensuing lambda needs to be mutable, else it makes output
          // const.

          if (allsame) {
            ::Tpetra::for_each
              (kernelLabel, execSpace, output,
               KOKKOS_LAMBDA (IST& X_ij) {
                X_ij = f (X_ij, X_ij);
              });
          }
          else if (outin1same) {
            withLocalAccess
              ([=] (const input2_view_type& input2_lcl) mutable {
                 ::Tpetra::for_each
                  (kernelLabel, execSpace, output,
                   KOKKOS_LAMBDA (IST& out, const LO i, const LO j) {
                    out = f (out, input2_lcl(i, j));
                  });
               }, readOnly (input2).on (memSpace));
          }
          else if (outin2same) {
            withLocalAccess
              ([=] (const input1_view_type& input1_lcl) mutable {
                 ::Tpetra::for_each
                  (kernelLabel, execSpace, output,
                   KOKKOS_LAMBDA (IST& out, const LO i, const LO j) {
                    out = f (input1_lcl(i, j), out);
                  });
               }, readOnly (input1).on (memSpace));
          }
          else if (in1in2same) {
            withLocalAccess
              ([=] (const input1_view_type& input1_lcl) mutable {
                 ::Tpetra::for_each
                  (kernelLabel, execSpace, output,
                   KOKKOS_LAMBDA (IST& out, const LO i, const LO j) {
                    out = f (input1_lcl(i, j), input1_lcl(i, j));
                  });
               }, readOnly (input1).on (memSpace));
          }
          else {
            withLocalAccess
              ([=] (const input1_view_type& input1_lcl) mutable {
                 withLocalAccess
                   ([=] (const input2_view_type& input2_lcl) mutable {
                      ::Tpetra::for_each
                       (kernelLabel, execSpace, output,
                        KOKKOS_LAMBDA (IST& out, const LO i, const LO j) {
                         out = f (input1_lcl(i, j), input2_lcl(i, j));
                       });
                    }, readOnly (input2).on (memSpace));
               }, readOnly (input1).on (memSpace));
            
            // withLocalAccess
            //   ([=] (const input1_view_type& input1_lcl,
            //         const input2_view_type& input2_lcl) mutable {
            //           ::Tpetra::for_each
            //            (kernelLabel, execSpace, output,
            //             KOKKOS_LAMBDA (IST& out, const LO i, const LO j) {
            //              out = f (input1_lcl(i, j), input2_lcl(i, j));
            //            });
            //    },
            //    readOnly (input1).on (memSpace),
            //    readOnly (input2).on (memSpace));
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

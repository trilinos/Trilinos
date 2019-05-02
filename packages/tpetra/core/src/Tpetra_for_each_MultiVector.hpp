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

#ifndef TPETRA_FOR_EACH_MULTIVECTOR_HPP
#define TPETRA_FOR_EACH_MULTIVECTOR_HPP

#include "Tpetra_withLocalAccess_MultiVector.hpp"
#include "Tpetra_for_each.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>
#include <sstream>

/// \file Tpetra_for_each_MultiVector.hpp
/// \brief Include this file to make Tpetra::for_each work with
///   Tpetra::MultiVector and Tpetra::Vector.
///
/// The overload of Tpetra::for_each for Tpetra::MultiVector resp.
/// Tpetra::Vector applies a function entrywise to each entry of the
/// MultiVector resp. Vector.  It works analogously to std::for_each.
///
/// The function f may have any of the following forms:
/// <ul>
/// <li> Takes the current entry in read-write fashion as
///   <tt>impl_scalar_type&</tt>, the local row index as LO, and the
///   local column index as LO; </li>
/// <li> Takes the current entry in read-write fashion as
///   <tt>impl_scalar_type</tt> and the local row index as LO; or
///   </li>
/// <li> Takes the current entry in read-write fashion as
///   <tt>impl_scalar_type</tt>. </li>
/// </ul>
///
/// <tt>impl_scalar_type</tt> is a public typedef in
/// Tpetra::MultiVector and Tpetra::Vector.  When <tt>scalar_type</tt>
/// is <tt>std::complex<T></tt>, then <tt>impl_scalar_type</tt> is
/// <tt>Kokkos::complex<T></tt>.  Otherwise, <tt>scalar_type</tt> and
/// <tt>impl_scalar_type</tt> are the same.

namespace Tpetra {
  namespace Details {

    // for_each uses this for the loop body of a parallel_for or
    // parallel_reduce over the rows of a Tpetra::MultiVector with
    // constant stride and multiple columns.
    template<class ViewType,
             class InnerLoopBodyType,
             class LocalIndexType>
    struct MultiVectorOuterForEachLoopBody {
      static_assert (static_cast<int> (ViewType::Rank) == 2,
                     "ViewType must be a rank-2 Kokkos::View.");
      MultiVectorOuterForEachLoopBody (const ViewType& X_lcl,
                                       InnerLoopBodyType f) :
        X_lcl_ (X_lcl), f_ (f)
      {}
      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i) const
      {
        const LocalIndexType numCols =
          static_cast<LocalIndexType> (X_lcl_.extent (1));
        for (LocalIndexType j = 0; j < numCols; ++j) {
          f_ (X_lcl_(i,j), i, j);
        }
      };
      ViewType X_lcl_;
      InnerLoopBodyType f_;
    };

    // for_each uses this for the loop body of a parallel_for or
    // parallel_reduce over the rows of a Tpetra::Vector (or each
    // column of a Tpetra::MultiVector, if the Tpetra::MultiVector has
    // nonconstant stride or only a single column).
    template<class ViewType,
             class InnerLoopBodyType,
             class LocalIndexType>
    struct VectorOuterForEachLoopBody {
      static_assert (static_cast<int> (ViewType::Rank) == 1,
                     "ViewType must be a rank-1 Kokkos::View.");
      VectorOuterForEachLoopBody (const ViewType& X_lcl,
                                  InnerLoopBodyType f) :
        X_lcl_ (X_lcl), f_ (f)
      {}
      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i) const
      {
        f_ (X_lcl_(i), i, LocalIndexType (0));
      };
      ViewType X_lcl_;
      InnerLoopBodyType f_;
    };

    // Distinguish between functions that take (scalar&, index,
    // index), (scalar&, index), and (scalar&).
    enum class EMultiVectorForEachFuncArgs {
      SCALAR,
      SCALAR_ROWINDEX,
      SCALAR_ROWINDEX_COLINDEX,
      ERROR
    };

    // Whether UserFunctionType is convertible to
    // std::function<void(ImplScalarType&, const LocalIndexType, const
    // LocalIndexType)> or std::function<void(ImplScalarType&,
    // LocalIndexType, LocalIndexType)>.
    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType>
    constexpr bool
    isScalarIndexIndexFunction ()
    {
      using UFT = UserFunctionType;
      using IST = ImplScalarType;
      using LIT = LocalIndexType;
      using func_type_1 =
        std::function<void (IST&, const LIT, const LIT)>;
      using func_type_2 = std::function<void (IST&, LIT, LIT)>;

      return std::is_convertible<UFT, func_type_1>::value ||
        std::is_convertible<UFT, func_type_2>::value;
    }

    // Whether UserFunctionType is convertible to
    // std::function<void(ImplScalarType&, const LocalIndexType)> or
    // std::function<void(ImplScalarType&, LocalIndexType)>.
    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType>
    constexpr bool
    isScalarIndexFunction ()
    {
      using UFT = UserFunctionType;
      using IST = ImplScalarType;
      using LIT = LocalIndexType;
      using func_type_1 = std::function<void (IST&, const LIT)>;
      using func_type_2 = std::function<void (IST&, LIT)>;

      return std::is_convertible<UFT, func_type_1>::value ||
        std::is_convertible<UFT, func_type_2>::value;
    }

    // Whether UserFunctionType is convertible to
    // std::function<void(ImplScalarType&)>.
    template<class UserFunctionType, class ImplScalarType>
    constexpr bool
    isScalarFunction ()
    {
      using func_type_1 = std::function<void (ImplScalarType&)>;

      return std::is_convertible<UserFunctionType, func_type_1>::value;
    }

    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType>
    constexpr EMultiVectorForEachFuncArgs
    getMultiVectorForEachFuncArgs ()
    {
      using UFT = UserFunctionType;
      using IST = ImplScalarType;
      using LIT = LocalIndexType;

      return isScalarIndexIndexFunction<UFT, IST, LIT> () ?
        EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX_COLINDEX :
        (isScalarIndexFunction<UFT, IST, LIT> () ?
         EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX :
         (isScalarFunction<UFT, IST> () ?
          EMultiVectorForEachFuncArgs::SCALAR :
          EMultiVectorForEachFuncArgs::ERROR));
    }

    // Functor that MultiVectorOuterForEachLoopBody or
    // VectorOuterForEachLoopBody uses.  This functor in turn wraps
    // the user's function given to for_each.  We have different cases
    // for whether the user's function takes (scalar&, local row
    // index, local column index), (scalar&, local row index), or
    // (scalar&).
    //
    // The point of MultiVectorInnerForEachLoopBody is so that the
    // MultiVector / Vector specialization of ForEach only needs one
    // implementation, but can work with three different kinds of user
    // functions: (scalar&, local row index, local column index),
    // (scalar&, local row index), and (scalar&).
    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType,
             const EMultiVectorForEachFuncArgs argsType =
               getMultiVectorForEachFuncArgs<UserFunctionType,
                                             ImplScalarType,
                                             LocalIndexType> ()>
    struct MultiVectorInnerForEachLoopBody {
      static_assert (argsType != EMultiVectorForEachFuncArgs::ERROR,
                     "Please report this bug to the Tpetra developers.");
    };

    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType>
    struct MultiVectorInnerForEachLoopBody<
      UserFunctionType, ImplScalarType, LocalIndexType,
      EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX_COLINDEX>
    {
      MultiVectorInnerForEachLoopBody (UserFunctionType f) :
        f_ (f) {}
      KOKKOS_INLINE_FUNCTION void
      operator () (ImplScalarType& x_ij,
                   const LocalIndexType i,
                   const LocalIndexType j) const
      {
        f_ (x_ij, i, j);
      };
      UserFunctionType f_;
    };

    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType>
    struct MultiVectorInnerForEachLoopBody<
      UserFunctionType, ImplScalarType, LocalIndexType,
      EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX>
    {
      MultiVectorInnerForEachLoopBody (UserFunctionType f) :
        f_ (f) {}
      KOKKOS_INLINE_FUNCTION void
      operator () (ImplScalarType& x_ij,
                   const LocalIndexType i,
                   const LocalIndexType /* j */) const
      {
        f_ (x_ij, i);
      };
      UserFunctionType f_;
    };

    template<class UserFunctionType,
             class ImplScalarType,
             class LocalIndexType>
    struct MultiVectorInnerForEachLoopBody<
      UserFunctionType, ImplScalarType, LocalIndexType,
      EMultiVectorForEachFuncArgs::SCALAR>
    {
      MultiVectorInnerForEachLoopBody (UserFunctionType f) :
        f_ (f) {}
      KOKKOS_INLINE_FUNCTION void
      operator () (ImplScalarType& x_ij,
                   const LocalIndexType /* i */,
                   const LocalIndexType /* j */) const
      {
        f_ (x_ij);
      };
      UserFunctionType f_;
    };

    // The MultiVector specialization of ForEach uses the result of
    // makeMultiVectorForEachLoopBody or makeVectorForEachLoopBody as
    // the functor in a parallel_for over the local rows of the
    // Tpetra::(Multi)Vector.

    template<class ViewType,
             class UserFunctionType,
             class LocalIndexType>
    MultiVectorOuterForEachLoopBody<
      ViewType,
      MultiVectorInnerForEachLoopBody<
        UserFunctionType,
        typename ViewType::non_const_value_type,
        LocalIndexType>,
      LocalIndexType>
    makeMultiVectorForEachLoopBody (const ViewType& X_lcl,
                                    UserFunctionType f,
                                    const LocalIndexType /* numCols */)
    {
      using return_type = typename ViewType::non_const_value_type;
      using inner_loop_body_type =
        MultiVectorInnerForEachLoopBody<
          UserFunctionType, return_type, LocalIndexType>;
      using outer_loop_body_type =
        MultiVectorOuterForEachLoopBody<
          ViewType, inner_loop_body_type, LocalIndexType>;
      return outer_loop_body_type (X_lcl, inner_loop_body_type (f));
    }

    template<class ViewType,
             class UserFunctionType,
             class LocalIndexType>
    VectorOuterForEachLoopBody<
      ViewType,
      MultiVectorInnerForEachLoopBody<
        UserFunctionType,
        typename ViewType::non_const_value_type,
        LocalIndexType>,
      LocalIndexType>
    makeVectorForEachLoopBody (const ViewType& X_lcl,
                               UserFunctionType f,
                               const LocalIndexType /* numCols */)
    {
      using return_type = typename ViewType::non_const_value_type;
      using inner_loop_body_type =
        MultiVectorInnerForEachLoopBody<
          UserFunctionType, return_type, LocalIndexType>;
      using outer_loop_body_type =
        VectorOuterForEachLoopBody<
          ViewType, inner_loop_body_type, LocalIndexType>;
      return outer_loop_body_type (X_lcl, inner_loop_body_type (f));
    }

    /// \brief Implementation of Tpetra::for_each for
    ///   Tpetra::MultiVector.
    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT,
             class UserFunctionType>
    struct ForEach<ExecutionSpace,
                   ::Tpetra::MultiVector<SC, LO, GO, NT>,
                   UserFunctionType>
    {
    private:
      // Given a Kokkos execution space on which the user wants to run
      // the for_each, and a memory space in which the MultiVector's
      // data live, determine the memory space that for_each should
      // use in its withLocalAccess call.
      template<class MemorySpace>
      using for_each_memory_space =
        typename std::conditional<
          Kokkos::SpaceAccessibility<
            ExecutionSpace,
            typename MemorySpace::memory_space>::accessible,
          typename MemorySpace::memory_space,
          typename ExecutionSpace::memory_space>::type;

    public:
      static void
      for_each (const char kernelLabel[],
                ExecutionSpace execSpace,
                ::Tpetra::MultiVector<SC, LO, GO, NT>& X,
                UserFunctionType f)
      {
        using Teuchos::TypeNameTraits;
        using std::endl;
        using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
        using preferred_memory_space =
          typename MV::device_type::memory_space;
        using memory_space = for_each_memory_space<preferred_memory_space>;
        using range_type = Kokkos::RangePolicy<ExecutionSpace, LO>;

        const int myRank = X.getMap ()->getComm ()->getRank ();
        const bool verbose = ::Tpetra::Details::Behavior::verbose ();
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank << ": Tpetra::for_each:" << endl
             << " kernelLabel: " << kernelLabel << endl
             << " ExecutionSpace: "
             << TypeNameTraits<ExecutionSpace>::name () << endl
             << " memory_space: "
             << TypeNameTraits<memory_space>::name () << endl;
          std::cerr << os.str ();
        }

        memory_space memSpace;
        if (X.getNumVectors () == size_t (1) || ! X.isConstantStride ()) {
          const size_t numVecs = X.getNumVectors ();
          for (size_t j = 0; j < numVecs; ++j) {
            auto X_j = X.getVectorNonConst (j);
            // Generic lambdas need C++14, so we must use arg_type here.
            using arg_type =
              with_local_access_function_argument_type<
                decltype (readWrite (*X_j).on (memSpace))>;
            withLocalAccess
              ([=] (const arg_type& X_j_lcl) {
                auto loopBody =
                  Details::makeVectorForEachLoopBody (X_j_lcl, f, LO (1));
                Kokkos::parallel_for
                  (kernelLabel,
                   range_type (execSpace, 0, X_j_lcl.extent (0)),
                   loopBody);
              }, readWrite (*X_j).on (memSpace));
          }
        }
        else {
          // Generic lambdas need C++14, so we must use arg_type here.
          using arg_type =
            with_local_access_function_argument_type<
              decltype (readWrite (X).on (memSpace))>;
          withLocalAccess
            ([=] (const arg_type& X_lcl) {
              const LO numCols = static_cast<LO> (X_lcl.extent (1));
              auto loopBody =
                Details::makeMultiVectorForEachLoopBody (X_lcl, f, numCols);
              Kokkos::parallel_for
                (kernelLabel,
                 range_type (execSpace, 0, X_lcl.extent (0)),
                 loopBody);
            }, readWrite (X).on (memSpace));
        }
      }
    };

    //! Implementation of Tpetra::for_each for Tpetra::Vector.
    template<class ExecutionSpace,
             class SC, class LO, class GO, class NT,
             class UserFunctionType>
    struct ForEach<ExecutionSpace,
                   ::Tpetra::Vector<SC, LO, GO, NT>,
                   UserFunctionType>
    {
      static void
      for_each (const char kernelLabel[],
                ExecutionSpace execSpace,
                ::Tpetra::Vector<SC, LO, GO, NT>& X,
                UserFunctionType f)
      {
        using MV = ::Tpetra::MultiVector<SC, LO, GO, NT>;
        using impl_type = ForEach<ExecutionSpace, MV, UserFunctionType>;
        impl_type::for_each (kernelLabel, execSpace, X, f);
      }
    };

  } // namespace Details

} // namespace Tpetra

#endif // TPETRA_FOR_EACH_MULTIVECTOR_HPP


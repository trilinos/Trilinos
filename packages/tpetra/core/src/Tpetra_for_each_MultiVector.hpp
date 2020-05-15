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

#include "Tpetra_for_each.hpp"
#include "Tpetra_withLocalAccess_MultiVector.hpp"
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
/// The input function f takes the current entry in read-write fashion
/// as <tt>impl_scalar_type&</tt>.  <tt>impl_scalar_type</tt> is a
/// public typedef in Tpetra::MultiVector and Tpetra::Vector.  When
/// <tt>scalar_type</tt> is <tt>std::complex<T></tt>, then
/// <tt>impl_scalar_type</tt> is <tt>Kokkos::complex<T></tt>.
/// Otherwise, <tt>scalar_type</tt> and <tt>impl_scalar_type</tt> are
/// the same.

namespace Tpetra {
  namespace Details {

    // for_each uses this for the loop body of a parallel_for or
    // parallel_reduce over the rows of a Tpetra::MultiVector with
    // constant stride and multiple columns.
    template<class ViewType,
             class UserFunctionType,
             class LocalIndexType>
    class MultiVectorForEachLoopBody {
    private:
      static_assert (static_cast<int> (ViewType::Rank) == 2,
                     "ViewType must be a rank-2 Kokkos::View.");
    public:
      MultiVectorForEachLoopBody (const ViewType& X_lcl,
                                  UserFunctionType f) :
        X_lcl_ (X_lcl), f_ (f)
      {}
      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i) const {
        const LocalIndexType numCols =
          static_cast<LocalIndexType> (X_lcl_.extent (1));
        for (LocalIndexType j = 0; j < numCols; ++j) {
          f_ (X_lcl_(i,j));
        }
      }
    private:
      ViewType X_lcl_;
      UserFunctionType f_;
    };

    // for_each uses this for the loop body of a parallel_for or
    // parallel_reduce over the rows of a Tpetra::Vector (or each
    // column of a Tpetra::MultiVector, if the Tpetra::MultiVector has
    // nonconstant stride or only a single column).
    template<class ViewType,
             class UserFunctionType,
             class LocalIndexType>
    class VectorForEachLoopBody {
    private:
      static_assert (static_cast<int> (ViewType::Rank) == 1,
                     "ViewType must be a rank-1 Kokkos::View.");
    public:
      VectorForEachLoopBody (const ViewType& X_lcl,
                             UserFunctionType f) :
        X_lcl_ (X_lcl), f_ (f)
      {}
      KOKKOS_INLINE_FUNCTION void
      operator () (const LocalIndexType i) const {
        f_ (X_lcl_(i));
      }
    private:
      ViewType X_lcl_;
      UserFunctionType f_;
    };

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

        // Generic lambdas need C++14, so for now, we use
        // with_local_access_function_argument_type to get a named
        // typedef for our withLocalAccess lambda argument(s).

        range_type range (execSpace, 0, X.getLocalLength ());
        memory_space memSpace;
        if (X.getNumVectors () == size_t (1) || ! X.isConstantStride ()) {
          const size_t numVecs = X.getNumVectors ();
          for (size_t j = 0; j < numVecs; ++j) {
            auto X_j = X.getVectorNonConst (j);
            // Help GCC 4.9.3 deduce the type of *X_j.
            // See discussion here:
            // https://github.com/trilinos/Trilinos/pull/5115
            Tpetra::Vector<SC, LO, GO, NT>& X_j_ref = *X_j;
            using read_write_view_type =
              with_local_access_function_argument_type<
                decltype (readWrite (X_j_ref).on (memSpace). at(execSpace))>;
            withLocalAccess
              ([=] (const read_write_view_type& X_j_lcl) {
                using functor_type = VectorForEachLoopBody<
                  read_write_view_type, UserFunctionType, LO>;
                Kokkos::parallel_for (kernelLabel, range,
                                      functor_type (X_j_lcl, f));
              }, readWrite (X_j_ref).on (memSpace). at(execSpace));
          }
        }
        else {
          // Generic lambdas need C++14, so we need a typedef here.
          using read_write_view_type =
            with_local_access_function_argument_type<
              decltype (readWrite (X).on (memSpace). at(execSpace))>;
          withLocalAccess
            ([=] (const read_write_view_type& X_lcl) {
              using functor_type = MultiVectorForEachLoopBody<
                read_write_view_type, UserFunctionType, LO>;
              Kokkos::parallel_for (kernelLabel, range,
                                    functor_type (X_lcl, f));
            }, readWrite (X).on (memSpace). at(execSpace));
        }
      }
    };

    /// \brief Implementation of Tpetra::for_each for Tpetra::Vector.
    ///
    /// Even though Tpetra::Vector is a subclass of
    /// Tpetra::MultiVector, this needs to exist, since partial
    /// specializations don't recognize subclasses of their
    /// specialized template arguments.
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


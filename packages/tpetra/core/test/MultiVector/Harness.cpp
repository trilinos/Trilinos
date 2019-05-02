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

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_withLocalAccess.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <functional>

namespace Tpetra {

  ////////////////////////////////////////////////////////////
  // Specializations for Tpetra::MultiVector
  ////////////////////////////////////////////////////////////

  namespace Details {
    // Specialization of GetMasterLocalObject for Tpetra::MultiVector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const AccessMode am>
    struct GetMasterLocalObject<
      LocalAccess<
        Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am> > {
    public:
      using local_access_type =
        LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;
    private:
      using global_object_type =
        typename local_access_type::global_object_type;
      using memory_space = typename local_access_type::memory_space;
      static constexpr AccessMode access_mode =
        local_access_type::access_mode;
      using non_const_value_type =
        typename global_object_type::impl_scalar_type;
      using value_type = typename std::conditional<
          access_mode == AccessMode::ReadOnly,
          const non_const_value_type,
          non_const_value_type
        >::type;

      // FIXME (mfh 22 Oct 2018, 25 Apr 2019) Need to make sure
      // execution space matches.  If not, we would need to allocate a
      // new View, and then we should actually make the
      // std::unique_ptr's destructor "copy back."  This is why
      // master_local_object_type is a std::unique_ptr<view_type>, not
      // just a view_type.
      //
      // mfh 01 May 2019: For now, we avoid allocation and copy back,
      // by using only the Views available in the MV's DualView.
      using dual_view_type = typename global_object_type::dual_view_type;
      static constexpr bool is_host =
        std::is_same<memory_space, Kokkos::HostSpace>::value;
      using result_device_type = typename std::conditional<
        is_host,
        typename dual_view_type::t_host::device_type,
        typename dual_view_type::t_dev::device_type>::type;
      using view_type = Kokkos::View<
        value_type**,
        typename dual_view_type::t_dev::array_layout,
        result_device_type>;

    public:
      using master_local_object_type = std::unique_ptr<view_type>;

      static master_local_object_type
      get (local_access_type LA)
      {
        if (LA.isValid ()) {
          if (access_mode == Details::AccessMode::WriteOnly) {
            LA.G_.clear_sync_state ();
          }

          // The various templated methods want an execution space
          // rather than a memory space.  Otherwise, DualView of
          // CudaUVMSpace complains that HostSpace is not one of its
          // two memory spaces.  (Both the device and the host Views
          // of a DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)
          using execution_space = typename memory_space::execution_space;

          if (LA.G_.template need_sync<execution_space> ()) {
            LA.G_.template sync<execution_space> ();
          }
          if (access_mode != Details::AccessMode::ReadOnly) {
            LA.G_.template modify<execution_space> ();
          }

          // See note about "copy back" above.
          auto G_lcl_2d = LA.G_.template getLocalView<execution_space> ();
          // This converts the View to const if applicable.
          return std::unique_ptr<view_type> (new view_type (G_lcl_2d));
        }
        else { // invalid; return "null" Kokkos::View
          return std::unique_ptr<view_type> (new view_type ());
        }
      }
    };

    // Specialization of GetMasterLocalObject for Tpetra::Vector.
    template<class SC, class LO, class GO, class NT,
             class MemorySpace,
             const Details::AccessMode am>
    struct GetMasterLocalObject<
      LocalAccess<
        Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am> > {
    public:
      using local_access_type =
        LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am>;
    private:
      using global_object_type =
        typename local_access_type::global_object_type;
      using memory_space = typename local_access_type::memory_space;
      static constexpr AccessMode access_mode =
        local_access_type::access_mode;
      using non_const_value_type =
        typename global_object_type::impl_scalar_type;
      using value_type = typename std::conditional<
          access_mode == AccessMode::ReadOnly,
          const non_const_value_type,
          non_const_value_type
        >::type;

      // FIXME (mfh 22 Oct 2018, 25 Apr 2019) Need to make sure
      // execution space matches.  If not, we would need to allocate a
      // new View, and then we should actually make the
      // std::unique_ptr's destructor "copy back."  This is why
      // master_local_object_type is a std::unique_ptr<view_type>, not
      // just a view_type.
      //
      // mfh 01 May 2019: For now, we avoid allocation and copy back,
      // by using only the Views available in the MV's DualView.
      using dual_view_type = typename global_object_type::dual_view_type;
      static constexpr bool is_host =
        std::is_same<memory_space, Kokkos::HostSpace>::value;
      using result_device_type = typename std::conditional<
        is_host,
        typename dual_view_type::t_host::device_type,
        typename dual_view_type::t_dev::device_type>::type;
      using view_type = Kokkos::View<
        value_type*,
        typename dual_view_type::t_dev::array_layout,
        result_device_type>;

    public:
      using master_local_object_type = std::unique_ptr<view_type>;

      static master_local_object_type
      get (local_access_type LA)
      {
        if (LA.isValid ()) {
          if (access_mode == Details::AccessMode::WriteOnly) {
            LA.G_.clear_sync_state ();
          }

          // The various templated methods want an execution space
          // rather than a memory space.  Otherwise, DualView of
          // CudaUVMSpace complains that HostSpace is not one of its
          // two memory spaces.  (Both the device and the host Views
          // of a DualView of CudaUVMSpace have memory_space =
          // CudaUVMSpace.)
          using execution_space = typename memory_space::execution_space;

          if (LA.G_.template need_sync<execution_space> ()) {
            LA.G_.template sync<execution_space> ();
          }
          if (access_mode != Details::AccessMode::ReadOnly) {
            LA.G_.template modify<execution_space> ();
          }

          // See note about "copy back" above.
          auto G_lcl_2d = LA.G_.template getLocalView<execution_space> ();
          auto G_lcl_1d = Kokkos::subview (G_lcl_2d, Kokkos::ALL (), 0);
          // This converts the View to const if applicable.
          return std::unique_ptr<view_type> (new view_type (G_lcl_1d));
        }
        else { // invalid; return "null" Kokkos::View
          return std::unique_ptr<view_type> (new view_type ());
        }
      }
    };

    // Specialization of GetNonowningLocalObject for Kokkos::View.
    template<class DataType,
             class LayoutType,
             class MemorySpace>
    struct GetNonowningLocalObject<
      std::unique_ptr<
        Kokkos::View<DataType, LayoutType, MemorySpace>>>
    {
    private:
      using input_view_type =
        Kokkos::View<DataType, LayoutType, MemorySpace>;
      using output_view_type =
        Kokkos::View<DataType,
                     LayoutType,
                     MemorySpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    public:
      using master_local_object_type = std::unique_ptr<input_view_type>;
      using nonowning_local_object_type = output_view_type;

      static nonowning_local_object_type
      get (const master_local_object_type& M)
      {
        input_view_type* viewPtr = M.get ();
        return viewPtr == nullptr ?
          nonowning_local_object_type () :
          nonowning_local_object_type (*viewPtr);
      }
    };

    ////////////////////////////////////////////////////////////
    // Implementation details of for_each
    ////////////////////////////////////////////////////////////

    // Implementation detail of for_each (see below).  Given a Kokkos
    // execution space on which the user wants to run the for_each,
    // and a memory space in which the MultiVector's data live,
    // determine the memory space that for_each should use in its
    // withLocalAccess call.
    template<class ExecutionSpace, class MemorySpace>
    using for_each_memory_space =
      typename std::conditional<
        Kokkos::SpaceAccessibility<
          ExecutionSpace,
          typename MemorySpace::memory_space>::accessible,
        typename MemorySpace::memory_space,
        typename ExecutionSpace::memory_space>::type;

    // for_each uses this for the loop body of a parallel_for or
    // parallel_reduce over the rows of a Tpetra::MultiVector with
    // constant stride and multiple columns.
    template<class ViewType,
             class InnerLoopBodyType,
             class IndexType>
    struct MultiVectorOuterLoopBody {
      static_assert (static_cast<int> (ViewType::Rank) == 2,
                     "ViewType must be a rank-2 Kokkos::View.");
      MultiVectorOuterLoopBody (const ViewType& X_lcl, InnerLoopBodyType f) :
        X_lcl_ (X_lcl), f_ (f)
      {}
      KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
        const IndexType numCols = static_cast<IndexType> (X_lcl_.extent (1));
        for (IndexType j = 0; j < numCols; ++j) {
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
             class IndexType>
    struct VectorOuterLoopBody {
      static_assert (static_cast<int> (ViewType::Rank) == 1,
                     "ViewType must be a rank-1 Kokkos::View.");
      VectorOuterLoopBody (const ViewType& X_lcl, InnerLoopBodyType f) :
        X_lcl_ (X_lcl), f_ (f)
      {}
      KOKKOS_INLINE_FUNCTION void operator () (const IndexType i) const {
        f_ (X_lcl_(i), i, IndexType (0));
      };
      ViewType X_lcl_;
      InnerLoopBodyType f_;
    };

    // Distinguish between functions that take (scalar, index, index),
    // (scalar, index), and (scalar).
    enum class EMultiVectorForEachFuncArgs {
      SCALAR,
      SCALAR_ROWINDEX,
      SCALAR_ROWINDEX_COLINDEX,
      ERROR
    };

    template<class FunctionType, class ReturnType, class IndexType>
    constexpr bool
    isScalarIndexIndexFunction ()
    {
      using func_type_1 =
        std::function<void (ReturnType&, const IndexType, const IndexType)>;
      using func_type_2 =
        std::function<void (ReturnType&, IndexType, IndexType)>;

      return std::is_convertible<FunctionType, func_type_1>::value ||
        std::is_convertible<FunctionType, func_type_2>::value;
    }

    template<class FunctionType, class ReturnType, class IndexType>
    constexpr bool
    isScalarIndexFunction ()
    {
      using func_type_1 = std::function<void (ReturnType&, const IndexType)>;
      using func_type_2 = std::function<void (ReturnType&, IndexType)>;

      return std::is_convertible<FunctionType, func_type_1>::value ||
        std::is_convertible<FunctionType, func_type_2>::value;
    }

    template<class FunctionType, class ReturnType>
    constexpr bool
    isScalarFunction ()
    {
      using func_type_1 = std::function<void (ReturnType&)>;

      return std::is_convertible<FunctionType, func_type_1>::value;
    }

    template<class FunctionType, class ReturnType, class IndexType>
    constexpr EMultiVectorForEachFuncArgs
    getMultiVectorForEachFuncArgs ()
    {
      return isScalarIndexIndexFunction<FunctionType, ReturnType, IndexType> () ?
        EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX_COLINDEX :
        (isScalarIndexFunction<FunctionType, ReturnType, IndexType> () ?
         EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX :
         (isScalarFunction<FunctionType, ReturnType> () ?
          EMultiVectorForEachFuncArgs::SCALAR :
          EMultiVectorForEachFuncArgs::ERROR));
    }

    // Functor that MultiVectorOuterLoopBody or VectorOuterLoopBody
    // uses.  This functor in turn wraps the user's function given to
    // for_each.  We have different cases for whether the user's
    // function takes (scalar&, row index, column index), (scalar&,
    // row index), or (scalar&).
    //
    // The point of InnerLoopBody is so that Tpetra::for_each only
    // needs one implementation, but can work with three different
    // kinds of user functions: (scalar&, local row index, local
    // column index), (scalar&, local row index), and (scalar&).
    template<class UserFunctionType,
             class ReturnType,
             class IndexType,
             const EMultiVectorForEachFuncArgs argsType =
               getMultiVectorForEachFuncArgs<UserFunctionType,
                                             ReturnType,
                                             IndexType> ()>
    struct InnerLoopBody {
      static_assert (argsType != EMultiVectorForEachFuncArgs::ERROR,
                     "Please report this bug to the Tpetra developers.");
    };

    template<class UserFunctionType,
             class ReturnType,
             class IndexType>
    struct InnerLoopBody<
      UserFunctionType, ReturnType, IndexType,
      EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX_COLINDEX>
    {
      InnerLoopBody (UserFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION void
      operator () (ReturnType& x_ij,
                   const IndexType i,
                   const IndexType j) const
      {
        f_ (x_ij, i, j);
      };

      UserFunctionType f_;
    };

    template<class UserFunctionType,
             class ReturnType,
             class IndexType>
    struct InnerLoopBody<
      UserFunctionType, ReturnType, IndexType,
      EMultiVectorForEachFuncArgs::SCALAR_ROWINDEX>
    {
      InnerLoopBody (UserFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION void
      operator () (ReturnType& x_ij,
                   const IndexType i,
                   const IndexType /* j */) const
      {
        f_ (x_ij, i);
      };

      UserFunctionType f_;
    };

    template<class UserFunctionType,
             class ReturnType,
             class IndexType>
    struct InnerLoopBody<
      UserFunctionType, ReturnType, IndexType,
      EMultiVectorForEachFuncArgs::SCALAR>
    {
      InnerLoopBody (UserFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION void
      operator () (ReturnType& x_ij,
                   const IndexType /* i */,
                   const IndexType /* j */) const
      {
        f_ (x_ij);
      };

      UserFunctionType f_;
    };

    // The implementation of for_each uses the result of
    // makeMultiVectorLoopBody or makeVectorLoopBody as the functor in
    // a parallel_for over the local rows of the
    // Tpetra::(Multi)Vector.

    template<class ViewType,
             class UserFunctionType,
             class IndexType>
    MultiVectorOuterLoopBody<
      ViewType,
      InnerLoopBody<
        UserFunctionType,
        typename ViewType::non_const_value_type,
        IndexType>,
      IndexType>
    makeMultiVectorLoopBody (const ViewType& X_lcl,
                             UserFunctionType f,
                             const IndexType /* numCols */)
    {
      using return_type = typename ViewType::non_const_value_type;
      using inner_loop_body_type =
        InnerLoopBody<UserFunctionType, return_type, IndexType>;
      using outer_loop_body_type =
        MultiVectorOuterLoopBody<ViewType, inner_loop_body_type, IndexType>;
      return outer_loop_body_type (X_lcl, inner_loop_body_type (f));
    }

    template<class ViewType,
             class UserFunctionType,
             class IndexType>
    VectorOuterLoopBody<
      ViewType,
      InnerLoopBody<
        UserFunctionType,
        typename ViewType::non_const_value_type,
        IndexType>,
      IndexType>
    makeVectorLoopBody (const ViewType& X_lcl,
                        UserFunctionType f,
                        const IndexType /* numCols */)
    {
      using return_type = typename ViewType::non_const_value_type;
      using inner_loop_body_type =
        InnerLoopBody<UserFunctionType, return_type, IndexType>;
      using outer_loop_body_type =
        VectorOuterLoopBody<ViewType, inner_loop_body_type, IndexType>;
      return outer_loop_body_type (X_lcl, inner_loop_body_type (f));
    }

    // Implementation of for_each (see below).
    template<class ExecutionSpace,
             class TpetraMultiVectorType,
             class UserFunctionType>
    struct ForEach {
      static void
      for_each (const char kernelLabel[],
                ExecutionSpace execSpace,
                TpetraMultiVectorType& X,
                UserFunctionType f)
      {
        using Teuchos::TypeNameTraits;
        using std::endl;
        using MV = TpetraMultiVectorType;
        using preferred_memory_space =
          typename MV::device_type::memory_space;
        using memory_space = Details::for_each_memory_space<
          ExecutionSpace, preferred_memory_space>;
        using LO = typename MV::local_ordinal_type;
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
                  Details::makeVectorLoopBody (X_j_lcl, f, LO (1));
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
                Details::makeMultiVectorLoopBody (X_lcl, f, numCols);
              Kokkos::parallel_for
                (kernelLabel,
                 range_type (execSpace, 0, X_lcl.extent (0)),
                 loopBody);
            }, readWrite (X).on (memSpace));
        }
      }
    };
  } // namespace Details

  /// \brief Apply a function entrywise to each entry of a
  ///   Tpetra::MultiVector or Tpetra::Vector, analogously to
  ///   std::for_each.
  ///
  /// The function f may have any of the following forms:
  /// <ul>
  /// <li> Takes the current entry in read-write fashion as
  ///   <tt>impl_scalar_type&</tt>, the local row index as LO, and the
  ///   local column index as LO; </li>
  ///
  /// <li> Takes the current entry in read-write fashion as
  ///   <tt>impl_scalar_type</tt> and the local row index as LO; or
  ///   </li>
  ///
  /// <li> Takes the current entry in read-write fashion as
  ///   <tt>impl_scalar_type</tt>. </li>
  /// </ul>
  ///
  /// <tt>impl_scalar_type</tt> is a public typedef in
  /// Tpetra::MultiVector and Tpetra::Vector.  When
  /// <tt>scalar_type</tt> is <tt>std::complex<T></tt>, then
  /// <tt>impl_scalar_type</tt> is <tt>Kokkos::complex<T></tt>.
  /// Otherwise, <tt>scalar_type</tt> and <tt>impl_scalar_type</tt>
  /// are the same.
  ///
  /// \param kernellabel [in] Kernel label for Kokkos Profiling.
  /// \param execSpace [in] Kokkos execution space on which to run.
  /// \param X [in/out] MultiVector or Vector to modify.
  /// \param f [in] Function to apply to each entry of X.
  template<class SC, class LO, class GO, class NT,
           class UserFunctionType,
           class ExecutionSpace>
  void
  for_each (const char kernelLabel[],
            ExecutionSpace execSpace,
            Tpetra::MultiVector<SC, LO, GO, NT>& X,
            UserFunctionType f)
  {
    using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
    using impl_type =
      Details::ForEach<ExecutionSpace, MV, UserFunctionType>;
    impl_type::for_each (kernelLabel, execSpace, X, f);
  }

  /// \brief Overload of for_each (see above) that runs on X's
  ///   default Kokkos execution space.
  ///
  /// \param kernellabel [in] Kernel label for Kokkos Profiling.
  /// \param X [in/out] MultiVector to modify.
  /// \param f [in] Function to apply entrywise to X (could have
  ///   different signatures; see above).
  template<class SC, class LO, class GO, class NT,
           class UserFunctionType>
  void
  for_each (const char kernelLabel[],
            Tpetra::MultiVector<SC, LO, GO, NT>& X,
            UserFunctionType f)
  {
    using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
    using execution_space = typename MV::device_type::execution_space;
    using impl_type =
      Details::ForEach<execution_space, MV, UserFunctionType>;
    execution_space execSpace;
    impl_type::for_each (kernelLabel, execSpace, X, f);
  }

} // namespace Tpetra

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<>;
  using multivec_type = Tpetra::MultiVector<>;
  using vec_type = Tpetra::Vector<>;
  using GO = map_type::global_ordinal_type;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( VectorHarness, GetLocalObject )
  {
    using Tpetra::Details::getMasterLocalObject;
    using Tpetra::Details::getNonowningLocalObject;
    using Tpetra::readOnly;
    using Tpetra::readWrite;
    using Tpetra::writeOnly;
    const bool debug = ::Tpetra::Details::Behavior::debug ();

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Tpetra::Details::{getMasterLocalObject, "
      "getNonowningLocalObject} for MultiVector and Vector" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type mvec (map, numVecs);
    TEST_EQUALITY( mvec.getNumVectors (), numVecs );

    myOut << "Create a Vector, and make sure that "
      "it has exactly one vector (column)" << endl;
    vec_type vec (map);
    TEST_EQUALITY_CONST(vec.getNumVectors (), size_t (1));

    // Test read-only nonowning MultiVector access.
    {
      auto X_lcl_ro_owning = getMasterLocalObject (readOnly (mvec));
      auto X_lcl_ro = getNonowningLocalObject (X_lcl_ro_owning);
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
      //static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }

    // Test whether read-only access works with a const MultiVector&.
    {
      using MV = multivec_type;
      auto X_lcl_ro_owning = getMasterLocalObject (readOnly (mvec));
      auto X_lcl_ro = getNonowningLocalObject (X_lcl_ro_owning);
      // Make sure X_lcl_ro can be assigned to the type we expect it to
      // be.  It doesn't have to be that type, it just has to be
      // assignable to that type.
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
      // Errors look like this:
      //
      // error: ‘__T0’ has not been declared
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }

    // Test write-only nonowning MultiVector access.
    {
      auto X_lcl_wo_owning = getMasterLocalObject (writeOnly (mvec));
      auto X_lcl_wo = getNonowningLocalObject (X_lcl_wo_owning);
      Kokkos::View<double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_wo2 = X_lcl_wo;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_wo)::Rank == 2, "Rank is not 2");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_wo.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_wo.extent (1)) == numVecs );
    }

    // Test read-write nonowning MultiVector access.
    {
      auto X_lcl_rw_owning = getMasterLocalObject (readWrite (mvec));
      auto X_lcl_rw = getNonowningLocalObject (X_lcl_rw_owning);
      Kokkos::View<double**,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_rw2 = X_lcl_rw;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_rw)::Rank == 2, "Rank is not 2");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_rw.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_rw.extent (1)) == numVecs );
    }

    // Test read-write nonowning Vector access.
    {
      auto X_lcl_1d_ro_owning = getMasterLocalObject (readOnly (vec));
      auto X_lcl_1d_ro = getNonowningLocalObject (X_lcl_1d_ro_owning);
      Kokkos::View<const double*,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_1d_ro2 = X_lcl_1d_ro;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_1d_ro)::Rank == 1, "Rank is not 1");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_1d_ro.extent (0)) == numLocal );
    }

    // Test write-only nonowning Vector access.
    {
      auto X_lcl_1d_wo_owning = getMasterLocalObject (writeOnly (vec));
      auto X_lcl_1d_wo = getNonowningLocalObject (X_lcl_1d_wo_owning);
      Kokkos::View<double*,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_1d_wo2 = X_lcl_1d_wo;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_1d_wo)::Rank == 1, "Rank is not 1");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_1d_wo.extent (0)) == numLocal );
    }

    // Test read-write nonowning Vector access.
    {
      auto X_lcl_1d_wr_owning = getMasterLocalObject (readWrite (vec));
      auto X_lcl_1d_wr = getNonowningLocalObject (X_lcl_1d_wr_owning);
      Kokkos::View<double*,
                   Kokkos::LayoutLeft,
                   multivec_type::device_type,
                   Kokkos::MemoryUnmanaged> X_lcl_1d_wr2 = X_lcl_1d_wr;
      // mfh 30 Apr 2019: Commented out due to possible NVCC bug.
#ifndef KOKKOS_ENABLE_CUDA
      static_assert (decltype (X_lcl_1d_wr)::Rank == 1, "Rank is not 1");
#endif // KOKKOS_ENABLE_CUDA
      TEST_ASSERT( size_t (X_lcl_1d_wr.extent (0)) == numLocal );
    }

    //
    // Examples of using the result of getNonowningLocalObject in
    // Kokkos::parallel_for kernels.
    //

    {
      using execution_space = vec_type::device_type::execution_space;
      using memory_space = vec_type::device_type::memory_space;
      using LO = vec_type::local_ordinal_type;
      using range_type = Kokkos::RangePolicy<execution_space, LO>;

      auto X_lcl_1d_wo_owning =
        getMasterLocalObject (writeOnly (vec).on (memory_space ()));
      auto X_lcl_1d_wo = getNonowningLocalObject (X_lcl_1d_wo_owning);
      static_assert
        (std::is_same<
           decltype (X_lcl_1d_wo)::device_type::execution_space,
           vec_type::dual_view_type::t_dev::execution_space>::value,
         "Wrong execution space");
      Kokkos::parallel_for (
        "Device kernel for write-only Tpetra::Vector",
        range_type (0, LO (numLocal)),
        KOKKOS_LAMBDA (const LO lclRow) {
          X_lcl_1d_wo(lclRow) = 42.0;
        });
    }

    {
      using host_execution_space =
        vec_type::dual_view_type::t_host::execution_space;
      using LO = vec_type::local_ordinal_type;
      using range_type = Kokkos::RangePolicy<host_execution_space, LO>;

      auto X_lcl_1d_wo_owning =
        getMasterLocalObject (writeOnly (vec).on (Kokkos::HostSpace ()));
      auto X_lcl_1d_wo = getNonowningLocalObject (X_lcl_1d_wo_owning);
      static_assert
        (std::is_same<
           decltype (X_lcl_1d_wo)::device_type::execution_space,
           vec_type::dual_view_type::t_host::execution_space>::value,
         "Wrong execution space");
      // test with some not-device function
      Kokkos::parallel_for (
        "Host kernel for write-only Tpetra::Vector",
        range_type (0, LO (numLocal)),
        [=] (const LO lclRow) {
          std::pair<double, double> p {3.0, 4.0};
          X_lcl_1d_wo(lclRow) = p.first * p.second;
        });

      // Just plain modify some entries, in some sequential order.
      // Just in case LO is unsigned, start at +1.
      for (LO lclRowPlusOne = LO (numLocal);
           lclRowPlusOne > LO (0); --lclRowPlusOne) {
        const LO lclRow = lclRowPlusOne - LO (1);
        // operator[] for 1-D Kokkos::View does the same thing as
        // operator().
        X_lcl_1d_wo[lclRow] = double (lclRow) + 42.0;
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                         outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  TEUCHOS_UNIT_TEST( VectorHarness, ForEach )
  {
    using Tpetra::for_each;
    using Kokkos::ALL;
    using std::endl;
    using device_execution_space =
      typename multivec_type::device_type::execution_space;
    using LO = typename multivec_type::local_ordinal_type;
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    int lclSuccess = 0;
    int gblSuccess = 0;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Tpetra::for_each" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type X (map, numVecs);
    TEST_EQUALITY( X.getNumVectors (), numVecs );

    out << "Test for_each(MV, void(double&)): Set entries to 418" << endl;
    for_each ("X_ij=418", X, KOKKOS_LAMBDA (double& X_ij) {
        X_ij = 418.0;
      });
    {
      X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 418.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    out << "Test for_each(DefaultHostExecutionSpace, MV, "
      "void(double&)): Set entries to 777" << endl;
    for_each ("X_ij=777", Kokkos::DefaultHostExecutionSpace (), X,
              KOKKOS_LAMBDA (double& X_ij) {
                X_ij = 777.0;
              });
    {
      //X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 777.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << std::endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    out << "X.sync_device(); X.putScalar(666.0);" << endl;
    X.sync_device ();
    X.putScalar (666.0);
    // out << "Test for_each(device_execution_space (), "
    //   "MultiVector, void(double&))" << endl;
    // for_each ("X_ij=666", device_execution_space (), X,
    //           KOKKOS_LAMBDA (double& X_ij) {
    //     X_ij = 666.0;
    //   });
    {
      X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 666.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << std::endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    out << "Test for_each(DefaultHostExecutionSpace, MV, "
      "void(double&)): Set entries to 44" << endl;
    for_each ("X_ij=44", Kokkos::DefaultHostExecutionSpace (), X,
              KOKKOS_LAMBDA (double& X_ij) {
                X_ij = 44.0;
              });
    {
      //X.sync_host (); // Doesn't help with CUDA_LAUNCH_BLOCKING unset
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 44.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << std::endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    out << "Test for_each(MV, void(double&)): Set entries to 31" << endl;
    //Kokkos::fence (); // Doesn't help with CUDA_LAUNCH_BLOCKING unset
    for_each ("X_ij=31", X, KOKKOS_LAMBDA (double& X_ij) {
        X_ij = 31.0;
      });
    {
      X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 31.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    out << "Test for_each(MV, void(double&,LO)): Set entries to 93" << endl;
    for_each ("X_ij=93", X, KOKKOS_LAMBDA (double& X_ij, LO /* i */) {
        X_ij = 93.0;
      });
    {
      X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 93.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << std::endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    out << "Test for_each(MV, void(double&,LO,LO))" << endl;
    for_each ("X_ij=777", X, KOKKOS_LAMBDA (double& X_ij,
                                const LO /* i */,
                                const LO /* j */) {
                X_ij = 777.0;
              });
    {
      X.sync_host ();
      auto X_lcl = X.getLocalViewHost ();
      for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
        out << "Column " << j << std::endl;
        bool ok = true;
        for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
          const double expectedVal = 777.0;
          if (X_lcl(i,j) != expectedVal) {
            out << "X_lcl(" << i << "," << j << ") = " << X_lcl(i,j)
                << " != " << expectedVal << std::endl;
            ok = false;
          }
        }
        TEST_ASSERT( ok );
      }
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    myOut << "Create a Vector, and make sure that "
      "it has exactly one vector (column)" << endl;
    vec_type vec (map);
    TEST_EQUALITY_CONST(vec.getNumVectors (), size_t (1));

    // Exercise overload of for_each that runs on X's default
    // execution space, and whose function takes (scalar&, LO, LO)
    // arguments.  Exercise it for a Vector.
    for_each ("X_ij+=(i+1)+(j+1)", vec,
        KOKKOS_LAMBDA (double& X_ij, const LO i, const LO j) {
          X_ij += double (i+1.0) + double (j+1.0);
        });

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    // Exercise overload of for_each that runs on X's default
    // execution space, and whose function takes (scalar&, LO)
    // arguments.  Exercise it for a Vector.
    for_each ("X_ij+=i+1", vec, KOKKOS_LAMBDA (double& X_ij, const LO i) {
        X_ij += double (i+1.0);
      });

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    // Exercise overload of for_each that runs on X's default
    // execution space, and whose function takes (scalar&).
    // Exercise it for a Vector.
    for_each ("X_ij=42", vec, KOKKOS_LAMBDA (double& X_ij) {
        X_ij = 42.0;
      });

    {
      vec.sync_host ();
      auto vec_lcl = subview (vec.getLocalViewHost (), ALL (), 0);
      bool ok = true;
      for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
        if (vec_lcl(i) != 42.0) {
          ok = false;
        }
      }
      TEST_ASSERT( ok );
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    for_each ("X_ij+=1", vec, KOKKOS_LAMBDA (double& X_ij) {
        X_ij += 1.0;
      });

    {
      vec.sync_host ();
      auto vec_lcl = subview (vec.getLocalViewHost (), ALL (), 0);
      bool ok = true;
      for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
        if (vec_lcl(i) != 43.0) {
          out << "vec_lcl(" << i << ") = " << vec_lcl(i) << " != 43.0"
              << std::endl;
          ok = false;
        }
      }
      TEST_ASSERT( ok );
    }

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  TEUCHOS_UNIT_TEST( VectorHarness, WithLocalAccess )
  {
    using LO = typename multivec_type::local_ordinal_type;
    const bool debug = ::Tpetra::Details::Behavior::debug ();

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Tpetra::withLocalAccess" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type X (map, numVecs);
    TEST_EQUALITY( X.getNumVectors (), numVecs );

    using const_lcl_mv_type =
      Kokkos::View<const double**,
                   Kokkos::LayoutLeft,
                   typename multivec_type::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using nonconst_lcl_mv_type =
      Kokkos::View<double**,
                   Kokkos::LayoutLeft,
                   typename multivec_type::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    // using lcl_vec_type =
    //   Kokkos::View<double*,
    //                Kokkos::LayoutLeft,
    //                Kokkos::MemoryUnmanaged,
    //                typename multivec_type::device_type>;
    using Tpetra::readOnly;
    using Tpetra::withLocalAccess;

    withLocalAccess
      ([&] (const const_lcl_mv_type& X_lcl) {
         TEST_EQUALITY( static_cast<size_t> (X_lcl.extent (0)),
                        static_cast<size_t> (X.getLocalLength ()) );
       },
       readOnly (X));
    withLocalAccess
      ([&] (const nonconst_lcl_mv_type& X_lcl) {
         TEST_EQUALITY( static_cast<size_t> (X_lcl.extent (0)),
                        static_cast<size_t> (X.getLocalLength ()) );
       },
       writeOnly (X));
    withLocalAccess
      ([&] (const nonconst_lcl_mv_type& X_lcl) {
         TEST_EQUALITY( static_cast<size_t> (X_lcl.extent (0)),
                        static_cast<size_t> (X.getLocalLength ()) );
       },
       readWrite (X));
  }
} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}

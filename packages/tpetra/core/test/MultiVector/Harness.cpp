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
#include <Tpetra_Details_Behavior.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <functional>

namespace Tpetra {

  ////////////////////////////////////////////////////////////
  // Generic classes needed to implement withLocalAccess
  ////////////////////////////////////////////////////////////

  namespace Details {
    /// \brief Access intent.
    enum class AccessMode {
      ReadOnly,
      WriteOnly,
      ReadWrite
    };

    /// \brief Given a global object, get its default memory space
    ///   (both the type and the default instance thereof).
    template<class GlobalObjectType>
    struct DefaultMemorySpace {
      using type = typename GlobalObjectType::device_type::memory_space;

      // Given a global object, get its (default) memory space instance.
      static type space (const GlobalObjectType& /* G */) {
        // This stub just assumes that 'type' is default constructible.
        // In Kokkos, default-constructing a memory space instance just
        // gives the default memory space.
        return type ();
      }
    };

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template<class GlobalObjectType,
             class MemorySpace,
             const AccessMode am>
    class LocalAccess; // forward declaration
#endif // DOXYGEN_SHOULD_SKIP_THIS

    /// \brief Mapping from LocalAccess to the "master" local object type.
    ///
    /// The latter gets the local data from a global object, and holds
    /// on to it until after the user's function (input to
    /// withLocalAccess) returns.  "Holds on to it" is key: the master
    /// object owns (governs the lifetime of) any data accessible from
    /// any resulting nonowning local object (see
    /// GetNonowningLocalObject below).
    ///
    /// Specializations require the following two public features:
    /// <ul>
    /// <li> <tt>master_local_object_type</tt> typedef </li>
    /// <li> <tt>master_local_object_type get(LocalAccessType)</tt>
    ///      static method </li>
    /// </ul>
    template<class LocalAccessType>
    struct GetMasterLocalObject {};

    /// \brief Given a LocalAccess instance (which has a reference to
    ///   a global object), get an instance of its master local object.
    ///
    /// This may be a heavyweight operation.  In fact, implementations
    /// for particular global object types may make this an MPI
    /// collective.  For example, a reasonable generalization of
    /// LocalAccess with
    /// <tt>GlobalObjectType=Tpetra::MultiVector</tt>, is LocalAccess
    /// with <tt>GlobalObjectType=pair<Tpetra::MultiVector,
    /// Tpetra::Import></tt>.  At that point, it may make sense to
    /// have a method on LocalAccess that annotates it with an Import
    /// (or Export).
    ///
    /// Developers should not need to overload this function.  The
    /// right way is to specialize GetMasterLocalObject::get (see
    /// above).
    template<class LocalAccessType>
    typename GetMasterLocalObject<LocalAccessType>::master_local_object_type
    getMasterLocalObject (LocalAccessType LA) {
      return GetMasterLocalObject<LocalAccessType>::get (LA);
    }

    /// \brief Mapping from "master" local object type to the
    ///   nonowning "local view" type that users see (as arguments to
    ///   the function that they give to withLocalAccess).
    ///
    /// The master local object may encode the memory space and access
    /// mode, but the mapping to local view type may also need
    /// run-time information.
    ///
    /// Specializations require the following two public features:
    /// <ul>
    /// <li> <tt>nonowning_local_object_type</tt> typedef </li>
    /// <li> <tt>nonowning_local_object_type get(const MasterLocalObjectType&)</tt>
    ///      static method </li>
    /// </ul>
    template<class MasterLocalObjectType>
    struct GetNonowningLocalObject {};

    /// \brief Given a master local object, get an instance of a
    ///   nonowning local object.
    ///
    /// Users only ever see the nonowning local object, and subviews
    /// (slices) thereof.  This is supposed to be a lightweight
    /// operation.
    ///
    /// Developers should not need to overload this function.  The
    /// right way is to specialize GetNonowningLocalObject::get (see
    /// above).
    template<class MasterLocalObjectType>
    typename GetNonowningLocalObject<MasterLocalObjectType>::
    nonowning_local_object_type
    getNonowningLocalObject (const MasterLocalObjectType& master) {
      return GetNonowningLocalObject<MasterLocalObjectType>::get (master);
    }

    /// \brief Compile-time mapping from LocalAccess specialization
    ///   type to the corresponding withLocalAccess function argument.
    ///
    /// Developers: Please don't specialize this.  It's fine just as
    /// it is.
    ///
    /// Use the LocalAccess type as the template parameter to determine
    /// the type of the nonowning local view to the global object's data.
    /// This only works if GetMasterLocalObject has been specialized for
    /// these template parameters, and if GetNonowningLocalObject has
    /// been specialized for the resulting "master" local object type.
    template<class LocalAccessType>
    class LocalAccessFunctionArgument {
    private:
      using gmlo = GetMasterLocalObject<LocalAccessType>;
      using master_local_object_type =
        typename gmlo::master_local_object_type;
      using gnlo = GetNonowningLocalObject<master_local_object_type>;
    public:
      using type = typename gnlo::nonowning_local_object_type;
    };
  } // namespace Details

  /// \brief Type of the local object, that is an argument to the
  ///   function the user gives to withLocalAccess.
  ///
  /// \tparam LocalAccessType Specialization of LocalAccess.
  ///
  /// withLocalAccess takes zero or more global objects, each wrapped
  /// in a LocalAccess struct.  Wrapping happens indirectly, through
  /// e.g., the readOnly, writeOnly, or readWrite functions.
  /// withLocalAccess also takes a function, functor, or lambda with
  /// the same number of arguments (the "local objects") as the number
  /// of global objects.  If users don't have C++14 (generic lambdas),
  /// then they will need to know the types of the local objects.
  /// This alias maps directly from the LocalAccess struct type, to
  /// the corresponding local object type.
  template<class LocalAccessType>
  using with_local_access_function_argument_type =
    typename Details::LocalAccessFunctionArgument<LocalAccessType>::type;

  //////////////////////////////////////////////////////////////////////
  // Users call readOnly, writeOnly, and readWrite, in order to declare
  // how they intend to access a global object's local data.
  //////////////////////////////////////////////////////////////////////

  /// \brief Declare that you want to access the given global object's
  ///   local data in read-only mode, in the object's default memory
  ///   space.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    typename Details::DefaultMemorySpace<GlobalObjectType>::type,
    Details::AccessMode::ReadOnly>
  readOnly (GlobalObjectType&);

  /// \brief Declare that you want to access the given global object's
  ///   local data in read-only mode (overload for const
  ///   GlobalObjectType), in the object's default memory space.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    typename Details::DefaultMemorySpace<GlobalObjectType>::type,
    Details::AccessMode::ReadOnly>
  readOnly (const GlobalObjectType&);

  /// \brief Declare that you want to access the given global object's
  ///   local data in write-only mode, in the object's default memory
  ///   space.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    typename Details::DefaultMemorySpace<GlobalObjectType>::type,
    Details::AccessMode::WriteOnly>
  writeOnly (GlobalObjectType&);

  /// \brief Declare that you want to access the given global object's
  ///   local data in read-and-write mode, in the object's default
  ///   memory space.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    typename Details::DefaultMemorySpace<GlobalObjectType>::type,
    Details::AccessMode::ReadWrite>
  readWrite (GlobalObjectType&);

  ////////////////////////////////////////////////////////////
  // LocalAccess struct
  ////////////////////////////////////////////////////////////

  namespace Details {
    /// \brief Declaration of access intent for a global object.
    ///
    /// Users aren't supposed to make instances of this class.  They
    /// should use readOnly, writeOnly, or readWrite instead, then
    /// call instance methods like on() and valid() on the resulting
    /// LocalAccess instance.
    template<class GlobalObjectType,
             class MemorySpace,
             const AccessMode am>
    class LocalAccess {
    public:
      using global_object_type = GlobalObjectType;
      using memory_space = typename MemorySpace::memory_space;
      static constexpr AccessMode access_mode = am;

    private:
      using canonical_this_type = LocalAccess<global_object_type,
                                              memory_space,
                                              access_mode>;
    public:
      /// \brief Constructor.
      ///
      /// Users must NOT call the LocalAccess constructor directly.
      /// They should instead start by calling readOnly, writeOnly, or
      /// readWrite above.  They may then use instance methods like
      /// on() or valid() (see below).
      ///
      /// G is a reference, because we only access it in a delimited
      /// scope.  G is nonconst, because even read-only local access
      /// may modify G.  For example, G may give access to its local
      /// data via lazy allocation of a data structure that differs
      /// from its normal internal storage format.
      ///
      /// Memory spaces should behave like Kokkos memory spaces.
      /// Default construction should work and should get the default
      /// instance of the space.  Otherwise, it may make sense to get
      /// the default memory space from G.
      LocalAccess (global_object_type& G,
                   memory_space space = memory_space (),
                   const bool isValid = true) :
        G_ (G),
        space_ (space),
        valid_ (isValid)
      {}

      /// \brief Type that users see, that's an argument to the
      ///   function that they give to withLocalAccess.
      using function_argument_type =
        with_local_access_function_argument_type<canonical_this_type>;

    public:
      /// \brief Declare at run time whether you actually want to
      ///   access the object.
      ///
      /// \param isValid [in] If false, then the caller promises that
      ///   they won't actually access the object.
      ///
      /// If isValid is false, implementations should not spend any
      /// effort getting the master local object.  This may save time
      /// on allocating temporary space, copying from device to host,
      /// etc.  This implies that implementations must be able to
      /// construct "null" / empty master local objects.
      LocalAccess<GlobalObjectType, MemorySpace, am>
      valid (const bool isValid) const {
        return {this->G_, this->space_, isValid};
      }

      /// \brief Declare intent to access this object's local data in
      ///   a specific (Kokkos) memory space (instance).
      template<class NewMemorySpace>
      LocalAccess<GlobalObjectType, NewMemorySpace, am>
      on (NewMemorySpace space) const {
        return {this->G_, space, this->valid_};
      }

      //! Is access supposed to be valid?  (See valid() above.)
      bool isValid () const { return this->valid_; }

      /// \brief Memory space instance in which the user will access
      ///   local data.
      memory_space getSpace () const { return space_; }

    public:
      /// \brief Reference to the global object whose data the user
      ///   will access.
      ///
      /// Keep by reference, because this struct is only valid in a
      /// delimited scope.
      global_object_type& G_;

      /// \brief Memory space instance in which the user will access
      ///   local data.
      ///
      /// We assume that Kokkos memory spaces have shallow-copy
      /// semantics.
      memory_space space_;
      //! Will I actually need to access this object?
      bool valid_;

    private:
      // Nonmember "constructors"; see above for declarations.  This are
      // friends, because they are the only ways that users are supposed
      // to construct LocalAccess instances.
      template<class GOT> friend
      LocalAccess<GOT, typename Details::DefaultMemorySpace<GOT>::type,
                  AccessMode::ReadOnly> readOnly (GOT&);
      template<class GOT> friend
      LocalAccess<GOT, typename Details::DefaultMemorySpace<GOT>::type,
                  AccessMode::ReadOnly> readOnly (const GOT&);
      template<class GOT> friend
      LocalAccess<GOT, typename Details::DefaultMemorySpace<GOT>::type,
                  AccessMode::WriteOnly> writeOnly (GOT&);
      template<class GOT> friend
      LocalAccess<GOT, typename Details::DefaultMemorySpace<GOT>::type,
                  AccessMode::ReadWrite> readWrite (GOT&);
    };
  } // namespace Details

  ////////////////////////////////////////////////////////////
  // Implementations of readOnly, writeOnly, and readWrite
  ////////////////////////////////////////////////////////////

  template<class GOT>
  Details::LocalAccess<
    GOT,
    typename Details::DefaultMemorySpace<GOT>::type,
    Details::AccessMode::ReadOnly>
  readOnly (GOT& G)
  {
    return {G, Details::DefaultMemorySpace<GOT>::space (G), true};
  }

  template<class GOT>
  Details::LocalAccess<
    GOT,
    typename Details::DefaultMemorySpace<GOT>::type,
    Details::AccessMode::ReadOnly>
  readOnly (const GOT& G)
  {
    GOT& G_nc = const_cast<GOT&> (G);
    return {G_nc, Details::DefaultMemorySpace<GOT>::space (G_nc), true};
  }

  template<class GOT>
  Details::LocalAccess<
    GOT,
    typename Details::DefaultMemorySpace<GOT>::type,
    Details::AccessMode::WriteOnly>
  writeOnly (GOT& G)
  {
    return {G, Details::DefaultMemorySpace<GOT>::space (G), true};
  }

  template<class GOT>
  Details::LocalAccess<
    GOT,
    typename Details::DefaultMemorySpace<GOT>::type,
    Details::AccessMode::ReadWrite>
  readWrite (GOT& G)
  {
    return {G, Details::DefaultMemorySpace<GOT>::space (G), true};
  }

  ////////////////////////////////////////////////////////////
  // Implementation of withLocalAccess
  ////////////////////////////////////////////////////////////

  namespace Details {

    /////////////////////////////////////////////////////////////////
    // Use std::tuple as a compile-time list (Greenspun's 10th Rule)
    /////////////////////////////////////////////////////////////////

    // cons<T, std::tuple<Args...>>::type is std::tuple<T, Args...>.
    // This is the usual Lisp CONS function, but for std::tuple.  I
    // got the idea from
    // https://stackoverflow.com/questions/18701798/building-and-accessing-a-list-of-types-at-compile-time
    // but without "struct Void {};", and with head of new list
    // ordered first instead of last (hence "cons").
    template<class ...> struct cons;

    // (CONS SomeType NIL)
    template<class T, template <class ...> class List>
    struct cons<T, List<>> {
      using type = List<T>;
    };

    // (CONS SomeType ListOfTypes)
    template <class T, template <class ...> class List, class ...Types>
    struct cons<T, List<Types...>>
    {
      typedef List<T, Types...> type;
    };

    ////////////////////////////////////////////////////////////
    // Map from std::tuple<Ts...> to std::function<void (Ts...)>.
    ////////////////////////////////////////////////////////////

    // I got inspiration from
    // https://stackoverflow.com/questions/15418841/how-do-i-strip-a-tuple-back-into-a-variadic-template-list-of-types
    //
    // The only significant change change was from "using type =
    // T<Ts...>;", to "using type = std::function<void (Ts...)>;".
    template<class T>
    struct tuple_to_function_type { };

    template<typename... Ts>
    struct tuple_to_function_type<std::tuple<Ts...> >
    {
      using type = std::function<void (Ts...)>;
    };

    // Map from a list of zero or more LocalAccess types, to the
    // corresponding list of arguments for the user function to give to
    // withLocalAccess.
    template<class ... Args>
    struct ArgsToFunction {};

    template<>
    struct ArgsToFunction<> {
      using arg_list_type = std::tuple<>;

      // Implementers should only use this.
      using type = std::function<void ()>;
    };

    template<class FirstLocalAccessType, class ... Rest>
    struct ArgsToFunction<FirstLocalAccessType, Rest...> {
      using head_arg_type =
        with_local_access_function_argument_type<FirstLocalAccessType>;
      using tail_arg_list_type =
        typename ArgsToFunction<Rest...>::arg_list_type;
      using arg_list_type =
        typename cons<head_arg_type, tail_arg_list_type>::type;

      // Implementers should only use this.
      using type = typename tuple_to_function_type<arg_list_type>::type;
    };

    /// \brief Implementation of withLocalAccess.
    ///
    /// \tparam LocalAccessTypes Zero or more possibly different
    ///   specializations of LocalAccess.
    template<class ... LocalAccessTypes>
    struct WithLocalAccess {
      using current_user_function_type =
        typename Details::ArgsToFunction<LocalAccessTypes...>::type;

      static void
      withLocalAccess (LocalAccessTypes...,
                       typename Details::ArgsToFunction<LocalAccessTypes...>::type);
    };

    /// \brief Specialization of withLocalAccess that implements the
    ///   "base class" of the user providing no GlobalObject
    ///   arguments, and a function that takes no arguments.
    template<>
    struct WithLocalAccess<> {
      using current_user_function_type =
        typename Details::ArgsToFunction<>::type;

      static void
      withLocalAccess (current_user_function_type userFunction)
      {
        userFunction ();
      }
    };

    /// \brief Specialization of withLocalAccess that implements the
    ///   "recursion case."
    ///
    /// \tparam FirstLocalAccessType Specialization of LocalAccess.
    /// \tparam Rest Zero or more possibly different specializations
    ///   of LocalAccess.
    template<class FirstLocalAccessType, class ... Rest>
    struct WithLocalAccess<FirstLocalAccessType, Rest...> {
      using current_user_function_type =
        typename Details::ArgsToFunction<FirstLocalAccessType, Rest...>::type;

      static void
      withLocalAccess (FirstLocalAccessType first,
                       Rest... rest,
                       current_user_function_type userFunction)
      {
        // The "master" local object is the scope guard for local
        // data.  Its constructor may allocate temporary storage, copy
        // data to the desired memory space, etc.  Its destructor will
        // put everything back.  "Put everything back" could be a
        // no-op, or it could copy data back so where they need to go
        // and/or free temporary storage.
        //
        // Users define this function and the type it returns by
        // specializing GetMasterLocalObject for LocalAccess
        // specializations.
        auto first_lcl_master = getMasterLocalObject (first);

        // The "nonowning" local object is a nonowning view of the
        // "master" local object.  This is the only local object that
        // users see, and they see it as input to their function.
        // Subsequent slices / subviews view this nonowning local
        // object.  All such nonowning views must have lifetime
        // contained within the lifetime of the master local object.
        //
        // Users define this function and the type it returns by
        // specializing GetNonowningLocalObject on the "master local
        // object" type (the type of first_lcl_master).
        //
        // Constraining the nonowning views' lifetime to this scope
        // means that master local object types may use low-cost
        // ownership models, like that of std::unique_ptr.  There
        // should be no need for reference counting (in the manner of
        // std::shared_ptr) or Herb Sutter's deferred_heap.
        auto first_lcl_view = getNonowningLocalObject (first_lcl_master);

        // Curry the user's function by fixing the first argument.

        // The commented-out implementation requires C++14, because it
        // uses a generic lambda (the special case where parameters
        // are "auto").  We do have the types of the arguments,
        // though, from ArgsToFunction, so we don't need this feature.

        // WithLocalAccess<Rest...>::withLocalAccess
        //   (rest...,
        //    [=] (auto ... args) {
        //      userFunction (first_lcl_view, args...);
        //    });

        WithLocalAccess<Rest...>::withLocalAccess
          (rest...,
           [=] (ArgsToFunction<Rest>... args) {
             userFunction (first_lcl_view, args...);
           });
      }
    };

    // Implementation detail of transform (see below).  Given a Kokkos
    // execution space on which the user wants to run the transform,
    // and a memory space in which the MultiVector's data live,
    // determine the memory space that transform should use in its
    // withLocalAccess call.
    template<class ExecutionSpace, class MemorySpace>
    using transform_memory_space =
      typename std::conditional<
        Kokkos::SpaceAccessibility<
          ExecutionSpace,
          typename MemorySpace::memory_space>::accessible,
        typename MemorySpace::memory_space,
        typename ExecutionSpace::memory_space>::type;

  } // namespace Details

  ////////////////////////////////////////////////////////////
  // withLocalAccess function declaration and definition
  ////////////////////////////////////////////////////////////

  // User's function goes first, followed by the list of zero or more
  // global objects to view, annotated by LocalAccess annotations (e.g.,
  // readOnly, writeOnly, readWrite).
  //
  // I would have preferred the LocalAccess arguments to go first.
  // However, C++ needs the user function has to go first, else C++
  // can't deduce the template parameters.  See the en.cppreference.com
  // article on "parameter pack."
  template<class ... LocalAccessTypes>
  void
  withLocalAccess
    (typename Details::ArgsToFunction<LocalAccessTypes...>::type userFunction,
     LocalAccessTypes... localAccesses)
  {
    using impl_type = Details::WithLocalAccess<LocalAccessTypes...>;
    impl_type::withLocalAccess (localAccesses..., userFunction);
  }

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
    // Implementation details of transform
    ////////////////////////////////////////////////////////////

    // transform uses this for the loop body of a parallel_for or
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
          X_lcl_(i,j) = f_ (X_lcl_(i,j), i, j);
        }
      };
      ViewType X_lcl_;
      InnerLoopBodyType f_;
    };

    // transform uses this for the loop body of a parallel_for or
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
        X_lcl_(i) = f_ (X_lcl_(i), i, IndexType (0));
      };
      ViewType X_lcl_;
      InnerLoopBodyType f_;
    };

    // Distinguish between functions that take (scalar, index, index),
    // (scalar, index), and (scalar).
    enum class EMultiVectorTransformFuncArgs {
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
        std::function<ReturnType (const ReturnType&, const IndexType, const IndexType)>;
      using func_type_2 =
        std::function<ReturnType (ReturnType, IndexType, IndexType)>;

      return std::is_convertible<FunctionType, func_type_1>::value ||
        std::is_convertible<FunctionType, func_type_2>::value;
    }

    template<class FunctionType, class ReturnType, class IndexType>
    constexpr bool
    isScalarIndexFunction ()
    {
      using func_type_1 =
        std::function<ReturnType (const ReturnType&, const IndexType)>;
      using func_type_2 =
        std::function<ReturnType (ReturnType, IndexType)>;

      return std::is_convertible<FunctionType, func_type_1>::value ||
        std::is_convertible<FunctionType, func_type_2>::value;
    }

    template<class FunctionType, class ReturnType>
    constexpr bool
    isScalarFunction ()
    {
      using func_type_1 =
        std::function<ReturnType (const ReturnType&)>;
      using func_type_2 =
        std::function<ReturnType (ReturnType)>;

      return std::is_convertible<FunctionType, func_type_1>::value ||
        std::is_convertible<FunctionType, func_type_2>::value;
    }

    template<class FunctionType, class ReturnType, class IndexType>
    constexpr EMultiVectorTransformFuncArgs
    getMultiVectorTransformFuncArgs ()
    {
      return isScalarIndexIndexFunction<FunctionType, ReturnType, IndexType> () ?
        EMultiVectorTransformFuncArgs::SCALAR_ROWINDEX_COLINDEX :
        (isScalarIndexFunction<FunctionType, ReturnType, IndexType> () ?
         EMultiVectorTransformFuncArgs::SCALAR_ROWINDEX :
         (isScalarFunction<FunctionType, ReturnType> () ?
          EMultiVectorTransformFuncArgs::SCALAR :
          EMultiVectorTransformFuncArgs::ERROR));
    }

    // Functor that MultiVectorOuterLoopBody or VectorOuterLoopBody
    // uses.  This functor in turn wraps the user's function given to
    // transform.  We have different cases for whether the user's
    // function takes (scalar, row index, column index), (scalar, row
    // index), or (scalar).
    template<class UserFunctionType,
             class ReturnType,
             class IndexType,
             const EMultiVectorTransformFuncArgs argsType =
               getMultiVectorTransformFuncArgs<UserFunctionType,
                                               ReturnType,
                                               IndexType> ()>
    struct InnerLoopBody {
      static_assert (argsType != EMultiVectorTransformFuncArgs::ERROR,
                     "Please report this bug to the Tpetra developers.");
    };

    template<class UserFunctionType,
             class ReturnType,
             class IndexType>
    struct InnerLoopBody<
      UserFunctionType, ReturnType, IndexType,
      EMultiVectorTransformFuncArgs::SCALAR_ROWINDEX_COLINDEX>
    {
      InnerLoopBody (UserFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION ReturnType
      operator () (const ReturnType& x_ij,
                   const IndexType i,
                   const IndexType j) const
      {
        return f_ (x_ij, i, j);
      };

      UserFunctionType f_;
    };

    template<class UserFunctionType,
             class ReturnType,
             class IndexType>
    struct InnerLoopBody<
      UserFunctionType, ReturnType, IndexType,
      EMultiVectorTransformFuncArgs::SCALAR_ROWINDEX>
    {
      InnerLoopBody (UserFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION ReturnType
      operator () (const ReturnType& x_ij,
                   const IndexType i,
                   const IndexType /* j */) const
      {
        return f_ (x_ij, i);
      };

      UserFunctionType f_;
    };

    template<class UserFunctionType,
             class ReturnType,
             class IndexType>
    struct InnerLoopBody<
      UserFunctionType, ReturnType, IndexType,
      EMultiVectorTransformFuncArgs::SCALAR>
    {
      InnerLoopBody (UserFunctionType f) : f_ (f) {}

      KOKKOS_INLINE_FUNCTION ReturnType
      operator () (const ReturnType& x_ij,
                   const IndexType /* i */,
                   const IndexType /* j */) const
      {
        return f_ (x_ij);
      };

      UserFunctionType f_;
    };

    // The implementation of transform uses the result of
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

    // Implementation of transform (see below).
    template<class ExecutionSpace,
             class TpetraMultiVectorType,
             class UserFunctionType>
    struct Transform {
      static void
      transform (const char debugLabel[],
                 ExecutionSpace execSpace,
                 TpetraMultiVectorType& X,
                 UserFunctionType f)
      {
        using Teuchos::TypeNameTraits;
        using std::endl;
        using MV = TpetraMultiVectorType;
        using preferred_memory_space =
          typename MV::device_type::memory_space;
        using memory_space = Details::transform_memory_space<
          ExecutionSpace, preferred_memory_space>;
        using LO = typename MV::local_ordinal_type;
        using range_type = Kokkos::RangePolicy<ExecutionSpace, LO>;

        const int myRank = X.getMap ()->getComm ()->getRank ();
        const bool verbose = ::Tpetra::Details::Behavior::verbose ();
        if (verbose) {
          std::ostringstream os;
          os << "Proc " << myRank << ": Tpetra::transform:" << endl
             << " debugLabel: " << debugLabel << endl
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
                  ("Tpetra::transform(Vector)",
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

              if (verbose) {
                std::ostringstream os;
                os << "Proc " << myRank << ": Contiguous MV case: X_lcl: "
                   << X_lcl.extent (0) << " x "
                   << X_lcl.extent (1) << endl;
                std::cerr << os.str ();
              }

              Kokkos::parallel_for
                ("Tpetra::transform(MultiVector)",
                 range_type (execSpace, 0, X_lcl.extent (0)),
                 loopBody);

            }, readWrite (X).on (memSpace));
        }
      }
    };
  } // namespace Details

  /// \brief Apply a function entrywise to each entry of a
  ///   Tpetra::MultiVector.
  ///
  /// X := f(X) entrywise, where X is a Tpetra::MultiVector and f
  /// has one of the following forms:
  ///
  /// <ul>
  /// <li> Takes the current entry as <tt>impl_scalar_type</tt>, the
  ///   local row index as LO, and the local column index as LO, and
  ///   returns <tt>impl_scalar_type</tt>; </li>
  ///
  /// <li> Takes the current entry as <tt>impl_scalar_type</tt> and
  ///   the local row index as LO, and returns
  ///   <tt>impl_scalar_type</tt>; </li>
  ///
  /// <li> Takes the current entry as <tt>impl_scalar_type</tt>, and
  ///   returns <tt>impl_scalar_type</tt>; </li>
  /// </ul>
  ///
  /// \param execSpace [in] Kokkos execution space on which to run.
  /// \param X [in/out] MultiVector to modify.
  /// \param f [in] Function to apply to each entry of X.
  template<class SC, class LO, class GO, class NT,
           class UserFunctionType,
           class ExecutionSpace>
  void
  transform (ExecutionSpace execSpace,
             Tpetra::MultiVector<SC, LO, GO, NT>& X,
             UserFunctionType f)
  {
    using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
    using impl_type =
      Details::Transform<ExecutionSpace, MV, UserFunctionType>;
    impl_type::transform ("transform(execSpace,MV,f)", execSpace, X, f);
  }

  /// \brief Overload of transform (see above) that runs on X's
  ///   default Kokkos execution space.
  ///
  /// \param X [in/out] MultiVector to modify.
  /// \param f [in] Function to apply entrywise to X (could have
  ///   different signatures; see above).
  template<class SC, class LO, class GO, class NT,
           class UserFunctionType>
  void
  transform (Tpetra::MultiVector<SC, LO, GO, NT>& X,
             UserFunctionType f)
  {
    using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
    using execution_space = typename MV::device_type::execution_space;
    using impl_type =
      Details::Transform<execution_space, MV, UserFunctionType>;
    execution_space execSpace;
    impl_type::transform ("transform(MV,f)", execSpace, X, f);
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

  TEUCHOS_UNIT_TEST( VectorHarness, Transform )
  {
    using Tpetra::transform;
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

    myOut << "Test Harness::transform" << endl;
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

    out << "Test transform(MV, double(double)): Set entries to 418" << endl;
    transform (X, KOKKOS_LAMBDA (const double /* X_ij */) {
        return double (418.0);
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

    out << "Test transform(DefaultHostExecutionSpace, MV, "
      "double(double)): Set entries to 777" << endl;
    transform (Kokkos::DefaultHostExecutionSpace (),
               X, KOKKOS_LAMBDA (const double /* X_ij */) {
        return double (777.0);
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
    // out << "Test transform(device_execution_space (), "
    //   "MultiVector, double(double))" << endl;
    // transform (device_execution_space (), X,
    //            KOKKOS_LAMBDA (const double /* X_ij */) {
    //     return double (666.0);
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

    out << "Test transform(DefaultHostExecutionSpace, MV, "
      "double(double)): Set entries to 44" << endl;
    transform (Kokkos::DefaultHostExecutionSpace (),
               X, KOKKOS_LAMBDA (const double /* X_ij */) {
        return double (44.0);
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

    out << "Test transform(MV, double(double)): Set entries to 31" << endl;
    //Kokkos::fence (); // Doesn't help with CUDA_LAUNCH_BLOCKING unset
    transform (X, KOKKOS_LAMBDA (const double /* X_ij */) {
        return double (31.0);
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

    out << "Test transform(MV, double(double,LO)): Set entries to 93" << endl;
    transform (X, KOKKOS_LAMBDA (double /* X_ij */, LO /* i */) {
        return double (93.0);
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

    out << "Test transform(MultiVector, double(double,LO,LO))" << endl;
    transform (X, KOKKOS_LAMBDA (const double /* X_ij */,
                                 const LO /* i */,
                                 const LO /* j */) {
        return double (777.0);
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

    // Exercise overload of transform that runs on X's default
    // execution space, and whose function takes (SC, LO, LO)
    // arguments.  Exercise it for a Vector.
    transform (vec, KOKKOS_LAMBDA (const double X_ij, const LO i, const LO j) {
        return X_ij + double (i+1.0) + double (j+1.0);
      });

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    // Exercise overload of transform that runs on X's default
    // execution space, and whose function takes (SC, LO)
    // arguments.  Exercise it for a Vector.
    transform (vec, KOKKOS_LAMBDA (const double X_ij, const LO i) {
        return X_ij + double (i+1.0);
      });

    lclSuccess = success ? 1 : 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (gblSuccess != 1) {
      out << "Returning early" << endl;
      return;
    }

    // Exercise overload of transform that runs on X's default
    // execution space, and whose function takes (SC) arguments.
    // Exercise it for a Vector.
    transform (vec, KOKKOS_LAMBDA (const double /* X_ij */) {
        return 42.0;
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

    transform (vec, KOKKOS_LAMBDA (const double X_ij) {
        return X_ij + 1.0;
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

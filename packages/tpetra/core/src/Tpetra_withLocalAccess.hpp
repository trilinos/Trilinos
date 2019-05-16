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

#ifndef TPETRA_WITHLOCALACCESS_HPP
#define TPETRA_WITHLOCALACCESS_HPP

#include "TpetraCore_config.h"
#include <functional>
#include <type_traits>

/// \file Tpetra_withLocalAccess.hpp
/// \brief Declaration and definition of Tpetra::withLocalAccess;
///   declaration of helper classes for users to specialize.

namespace Tpetra {

  ////////////////////////////////////////////////////////////
  // Generic classes needed to implement withLocalAccess
  ////////////////////////////////////////////////////////////

  namespace Details {
    //! Access intent.
    enum class AccessMode {
      ReadOnly,
      WriteOnly,
      ReadWrite
    };

    /// \brief Given a global object, get its default memory space
    ///   (both the type and the default instance thereof).
    ///
    /// This generic version should suffice for most Tpetra classes,
    /// and for any class with a public <tt>device_type</tt> typedef
    /// that is a Kokkos::Device specialization.
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
    /// \tparam LocalAccessType Specialization of LocalAccess.
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
    /// \tparam LocalAccessType Specialization of LocalAccess.
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
    /// \tparam LocalAccessType Specialization of LocalAccess, the
    ///   same as that used by the corresponding GetMasterLocalObject
    ///   specialization.
    ///
    /// The master local object may encode the memory space and access
    /// mode, but the mapping to local view type may also need
    /// run-time information.
    ///
    /// Specializations require the following public features:
    /// <ul>
    /// <li> <tt>nonowning_local_object_type</tt> typedef </li>
    /// <li> <tt>nonowning_local_object_type get(LocalAccessType,
    ///            const MasterLocalObjectType&)</tt>
    ///      static method </li>
    /// </ul>
    ///
    /// Implementations of get() need not use the first
    /// LocalAccessType argument.  It exists mainly to help the
    /// compiler deduce template parameters of the nonmember
    /// getNonowningLocalObject function below.
    template<class LocalAccessType>
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
    template<class LocalAccessType>
    typename GetNonowningLocalObject<LocalAccessType>::
    nonowning_local_object_type
    getNonowningLocalObject (LocalAccessType LA,
      const typename GetMasterLocalObject<LocalAccessType>::
        master_local_object_type& master)
    {
      return GetNonowningLocalObject<LocalAccessType>::get (LA, master);
    }
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
    typename Details::GetNonowningLocalObject<LocalAccessType>::
      nonowning_local_object_type;

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
    // The only significant change was from "using type = T<Ts...>;",
    // to "using type = std::function<void (Ts...)>;".
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
      withLocalAccess (current_user_function_type userFunction,
                       FirstLocalAccessType first,
                       Rest... rest)
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
        // specializing GetNonowningLocalObject (see above).
        //
        // Constraining the nonowning views' lifetime to this scope
        // means that master local object types may use low-cost
        // ownership models, like that of std::unique_ptr.  There
        // should be no need for reference counting (in the manner of
        // std::shared_ptr) or Herb Sutter's deferred_heap.
        auto first_lcl_view = getNonowningLocalObject (first, first_lcl_master);

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
          ([=] (typename Rest::function_argument_type... args) {
             userFunction (first_lcl_view, args...);
           },
           rest...);
      }
    };
  } // namespace Details

  ////////////////////////////////////////////////////////////
  // withLocalAccess function declaration and definition
  ////////////////////////////////////////////////////////////

  /// \brief Get access to a Tpetra global object's local data.
  ///
  /// \param userFunction [in] Lambda, closure, or function that takes
  ///   the local data structure by const reference.
  ///
  /// \param localAccesses [in] Zero or more specializations of
  ///   Details::LocalAccess.  You do not create these directly; they
  ///   result from applying the readOnly, writeOnly, or readWrite
  ///   functions (see above) to Tpetra global objects.
  ///
  /// A "Tpetra global object" is a Tpetra::DistObject subclass
  /// representing a data structure distributed over one or more MPI
  /// processes.  Examples include Tpetra::MultiVector,
  /// Tpetra::Vector, Tpetra::CrsGraph, and Tpetra::CrsMatrix.
  ///
  /// The type of the local data structure depends on the type of the
  /// Tpetra global object.  For example, Tpetra::MultiVector's
  /// local data structure is a rank-2 unmanaged (nonowning)
  /// Kokkos::View, and Tpetra::Vector's local data structure is a
  /// rank-2 unmanaged Kokkos::View.
  ///
  /// You may use the with_local_access_function_argument_type alias
  /// (see above) to figure out the type of the local data structure,
  /// as a function of the type of the corresponding "local access"
  /// argument.  Alternately, if your function is a lambda and you
  /// have C++14, you may use <tt>const auto&</tt> ("generic
  /// lambdas").
  ///
  /// \note I would have preferred the LocalAccess arguments to go
  ///   first.  However, C++ needs the user function has to go first,
  ///   else C++ can't deduce the template parameters.  See the
  ///   en.cppreference.com article on "parameter pack."
  template<class ... LocalAccessTypes>
  void
  withLocalAccess
    (typename Details::ArgsToFunction<LocalAccessTypes...>::type userFunction,
     LocalAccessTypes... localAccesses)
  {
    using impl_type = Details::WithLocalAccess<LocalAccessTypes...>;
    impl_type::withLocalAccess (userFunction, localAccesses...);
  }

} // namespace Tpetra

#endif // TPETRA_WITHLOCALACCESS_HPP


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
    /// \brief Enum for declaring access intent.
    ///
    /// This is not for users; it's an implementation detail of
    /// functions readOnly, writeOnly, and readWrite (see below).
    enum class EAccess {
      ReadOnly,
      WriteOnly,
      ReadWrite
    };

    /// \brief Tag class for declaring access intent.
    ///
    /// This is a class, not an enum, so that LocalAccess (see below)
    /// can have a parameter pack.  You can't mix class types and
    /// constexpr values in a parameter pack.  Compare to
    /// Kokkos::MemoryTraits or Kokkos::Schedule.
    template<const EAccess accessMode>
    class Access {
      /// \brief Type alias to help identifying the tag class.
      ///
      /// We follow Kokkos in identifying a tag class by checking
      /// whether a type alias inside it has the same type as the tag
      /// class itself.
      using access_mode = Access<accessMode>;
    };

    template<class T> struct is_access_mode : public std::false_type {};
    template<EAccess am>
    struct is_access_mode<Access<am>> : public std::true_type {};

    using read_only = Access<EAccess::ReadOnly>;
    using write_only = Access<EAccess::WriteOnly>;
    using read_write = Access<EAccess::ReadWrite>;

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

    /// \brief Given a global object, get its default execution space
    ///   (both the type and the default instance thereof).
    ///
    /// This generic version should suffice for most Tpetra classes,
    /// and for any class with a public <tt>device_type</tt> typedef
    /// that is a Kokkos::Device specialization.
    template<class GlobalObjectType>
    struct DefaultExecutionSpace {
      using type = typename GlobalObjectType::device_type::execution_space;

      // Given a global object, get its (default) execution space instance.
      static type space (const GlobalObjectType& /* G */) {
        // This stub just assumes that 'type' is default constructible.
        // In Kokkos, default-constructing a execution space instance just
        // gives the default execution space.
        return type ();
      }
    };

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template<class GlobalObjectType, class ... Args>
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
    /// <li> <tt>master_local_object_type</tt> type alias </li>
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
    typename GetMasterLocalObject<LocalAccessType>::
      master_local_object_type
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
    getNonowningLocalObject(LocalAccessType LA,
      const typename GetMasterLocalObject<LocalAccessType>::
        master_local_object_type& master)
    {
      using impl_type = GetNonowningLocalObject<LocalAccessType>;
      return impl_type::get(LA, master);
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
  ///   local data in read-only mode.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::read_only>
  readOnly (GlobalObjectType&);

  /// \brief Declare that you want to access the given global object's
  ///   local data in read-only mode (overload for const
  ///   GlobalObjectType).
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::read_only>
  readOnly (const GlobalObjectType&);

  /// \brief Declare that you want to access the given global object's
  ///   local data in write-only mode.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::write_only>
  writeOnly (GlobalObjectType&);

  /// \brief Declare that you want to access the given global object's
  ///   local data in read-and-write mode.
  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::read_write>
  readWrite (GlobalObjectType&);

  ////////////////////////////////////////////////////////////
  // LocalAccess struct
  ////////////////////////////////////////////////////////////

  namespace Details {
    //! Tag indicating an unspecified type in LocalAccessTraits.
    struct unspecified_type {};

    /// \brief Deduce types from parameter pack of LocalAccess.
    ///
    /// Deduce the following types from Args:
    /// <ul>
    /// <li>execution_space (Kokkos execution space)</li>
    /// <li>memory_space (Kokkos memory space)</li>
    /// <li>access_mode (Access type)</li>
    /// </ul>
    template<class ... Args>
    struct LocalAccessTraits {};

    // An empty list of template arguments constrains nothing.
    template<>
    struct LocalAccessTraits<> {
      using execution_space = unspecified_type;
      using memory_space = unspecified_type;
      using access_mode = unspecified_type;
    };

    // Skip over any instances of unspecified_type.  This makes it
    // easy to define the return types of LocalAccess::on and
    // LocalAccess::at (see below).
    template<class ... Rest>
    struct LocalAccessTraits<unspecified_type, Rest...> {
    private:
      using rest_type = LocalAccessTraits<Rest...>;
    public:
      using execution_space = typename rest_type::execution_space;
      using memory_space = typename rest_type::memory_space;
      using access_mode = typename rest_type::access_mode;
    };

    template<class First, class ... Rest>
    struct LocalAccessTraits<First, Rest...> {
    private:
      using rest_type = LocalAccessTraits<Rest...>;
    public:
      using execution_space = typename std::conditional<
        Kokkos::Impl::is_execution_space<First>::value,
        First,
        typename rest_type::execution_space>::type;
      using memory_space = typename std::conditional<
        Kokkos::Impl::is_memory_space<First>::value,
        First,
        typename rest_type::memory_space>::type;
      using access_mode = typename std::conditional<
        is_access_mode<First>::value,
        First,
        typename rest_type::access_mode>::type;
    };

    template<class GlobalObjectType,
             class Traits,
             bool is_execution_space_specified = ! std::is_same<
               typename Traits::execution_space,
               unspecified_type>::value,
             bool is_memory_space_specified = ! std::is_same<
               typename Traits::memory_space,
               unspecified_type>::value>
    struct SpaceTypes {};

    template<class GlobalObjectType, class Traits>
    struct SpaceTypes<GlobalObjectType, Traits, true, true> {
      using execution_space = typename Traits::execution_space;
      using memory_space = typename Traits::memory_space;
    };

    // If only memory_space was specified, then get execution_space
    // from memory_space's default execution space.
    template<class GlobalObjectType, class Traits>
    struct SpaceTypes<GlobalObjectType, Traits, false, true> {
      using execution_space =
        typename Traits::memory_space::execution_space;
      using memory_space = typename Traits::memory_space;
    };

    // If only execution_space was specified, then get memory_space
    // from execution_space's default memory space.
    template<class GlobalObjectType, class Traits>
    struct SpaceTypes<GlobalObjectType, Traits, true, false> {
      using execution_space = typename Traits::execution_space;
      using memory_space =
        typename Traits::execution_space::memory_space;
    };

    // If neither execution_space nor memory_space were specified,
    // then get them both from their defaults, that depend on
    // GlobalObjectType.
    template<class GlobalObjectType, class Traits>
    struct SpaceTypes<GlobalObjectType, Traits, false, false> {
      using execution_space =
        typename DefaultExecutionSpace<GlobalObjectType>::type;
      using memory_space =
        typename DefaultMemorySpace<GlobalObjectType>::type;
    };

    /// \brief Declaration of access intent for a global object.
    ///
    /// \tparam GlobalObjectType Type of the global object whose local
    ///   data you want to access.
    ///
    /// \tparam Args Zero or more types, each of which may be a Kokkos
    ///   execution space, a Kokkos memory space, or Access (see
    ///   above).
    ///
    /// Users must not make instances of this class directly.  The
    /// only way they may create an instance of LocalAccess is by
    /// calling readOnly(), writeOnly(), readWrite(), or LocalAccess
    /// instance methods like at(), on(), and valid().
    ///
    /// For X a Tpetra::MultiVector with
    /// <tt>memory_space=Kokkos::CudaUVMSpace</tt> and
    /// <tt>execution_space=Kokkos::Cuda</tt>, the following two code
    /// examples must have the same effect:
    ///
    /// \code
    /// readWrite(X).on(CudaUVMSpace()).at(DefaultHostExecutionSpace());
    /// \endcode
    ///
    /// and
    ///
    /// \code
    /// readWrite(X).at(DefaultHostExecutionSpace()).on(CudaUVMSpace());
    /// \endcode
    ///
    /// That effect should be: "I intend to view X's local data in
    /// read-and-write fashion at host execution, on a UVM
    /// allocation."  Given Tpetra's current design (as of 26 Jan
    /// 2020), if X needs sync to host, then that implies a fence, to
    /// ensure that any device writes to X's local data are done.
    ///
    /// This means that LocalAccess needs to be able to know whether
    /// the user has explicitly specified an execution space in which
    /// to access local data.  Here is how <tt>on</tt> behaves:
    ///
    /// 1. If the input LocalAccess has no execution space explicitly
    ///    specified, then do not constrain it.  withLocalAccess will
    ///    assign a default execution space.  (This makes the "at
    ///    after on" example above work correctly.)
    ///
    /// 2. Else, use the input LocalAccess' execution_space.  The
    ///    caller is responsible for knowing whether NewMemorySpace is
    ///    accessible from execution_space.
    ///
    /// Here is how <tt>at</tt> behaves:
    ///
    /// 1. If the input LocalAccess has no memory space explicitly
    ///    specified, then do not constrain it.  withLocalAccess will
    ///    assign a default memory space.  (This makes the "on after
    ///    at" example above work correctly.)
    ///
    /// 2. Else, use the input LocalAccess' memory_space.  The caller
    ///    is responsible for knowing whether NewExecutionSpace can
    ///    access memory_space.
    ///
    /// The general behavior is that at() and on() both only constrain
    /// the one thing that the user specified.  Only the "final"
    /// LocalAccess object determines the behavior of withLocalAccess.
    ///
    /// For these reasons, LocalAccess needs to be able to distinguish
    /// between "the user explicitly assigned an {execution, memory}
    /// space," and "use the default {execution, memory} space."  This
    /// is why LocalAccess has a template parameter pack.
    template<class GlobalObjectType, class ... Args>
    class LocalAccess {
    private:
      using this_type = LocalAccess<GlobalObjectType, Args...>;
      using traits = LocalAccessTraits<Args...>;
      static constexpr bool is_execution_space_specified =
        ! std::is_same<typename traits::execution_space,
                       unspecified_type>::value;
      static constexpr bool is_memory_space_specified =
        ! std::is_same<typename traits::memory_space,
                       unspecified_type>::value;
      static constexpr bool is_access_mode_specified =
        ! std::is_same<typename traits::access_mode,
                       unspecified_type>::value;
      static_assert(is_access_mode_specified, "LocalAccess requires "
        "that you specify the access mode.  You may do this by "
        "always starting construction of a LocalAccess instance "
        "with readOnly, writeOnly, or readWrite.");

    public:
      using global_object_type = GlobalObjectType;

      using execution_space =
        typename SpaceTypes<global_object_type, traits>::execution_space;
      static_assert(! is_execution_space_specified ||
        Kokkos::Impl::is_execution_space<execution_space>::value,
        "Specified execution space is not a Kokkos execution space.");
      static_assert(is_execution_space_specified ||
        Kokkos::Impl::is_execution_space<execution_space>::value,
        "Default execution space is not a Kokkos execution space.");

      using memory_space =
        typename SpaceTypes<global_object_type, traits>::memory_space;
      static_assert(! is_memory_space_specified ||
        Kokkos::Impl::is_memory_space<memory_space>::value,
        "Specified memory space is not a Kokkos memory space.");
      static_assert(is_memory_space_specified ||
        Kokkos::Impl::is_memory_space<memory_space>::value,
        "Default memory space is not a Kokkos memory space.");

      using access_mode = typename traits::access_mode;

      /// \brief Constructor that specifies the global object, the
      ///   execution memory space instance on which to view its local
      ///   data, the memory space instance on which to view its local
      ///   data, and (optionally) whether the global object is valid.
      ///
      /// Users must NOT call the LocalAccess constructor directly.
      /// They should instead start by calling readOnly, writeOnly, or
      /// readWrite above.  They may then use instance methods like
      /// at(), on(), or valid() (see below).
      ///
      /// G is a reference, because we only access it in a delimited
      /// scope.  G is nonconst, because even read-only local access
      /// may modify G.  For example, G may give access to its local
      /// data via lazy allocation of a data structure that differs
      /// from its normal internal storage format.
      LocalAccess (global_object_type& G,
                   const execution_space& execSpace,
                   const memory_space& memSpace,
                   const bool viewIsValid = true) :
        G_ (G),
        execSpace_ (execSpace),
        memSpace_ (memSpace),
        valid_ (viewIsValid)
      {}

      /// \brief Constructor that specifies the global object and
      ///   (optionally) whether it is valid.
      LocalAccess (global_object_type& G,
                   const bool viewIsValid = true) :
        LocalAccess (G,
          DefaultExecutionSpace<global_object_type>::space (G),
          DefaultMemorySpace<global_object_type>::space (G),
          viewIsValid)
      {}

      /// \brief Constructor that specifies the global object, the
      ///   execution space instance at which to view its local data,
      ///   and (optionally) whether the global object is valid.
      LocalAccess (global_object_type& G,
                   const execution_space& execSpace,
                   const bool viewIsValid = true) :
        LocalAccess (G,
          execSpace,
          DefaultMemorySpace<global_object_type>::space (G),
          viewIsValid)
      {}

      /// \brief Constructor that specifies the global object, the
      ///   memory space instance on which to view its local data, and
      ///   (optionally) whether the global object is valid.
      LocalAccess (global_object_type& G,
                   const memory_space& memSpace,
                   const bool viewIsValid = true) :
        LocalAccess (G,
          DefaultExecutionSpace<global_object_type>::space (G),
          memSpace,
          viewIsValid)
      {}

      /// \brief Is access supposed to be valid?
      ///
      /// If not, then you may not access the global object's local
      /// data inside a withLocalAccess scope governed by this
      /// LocalAccess instance.
      ///
      /// See <tt>valid(const bool)</tt> below.
      bool isValid () const { return valid_; }

      /// \brief Execution space instance, at which the user will
      ///   access local data.
      execution_space getExecutionSpace () const { return execSpace_; }

      /// \brief Memory space instance, on which the user will access
      ///   local data.
      memory_space getMemorySpace () const { return memSpace_; }

      /// \brief Declare at run time whether you actually want to
      ///   access the object.
      ///
      /// \param thisIsValid [in] If false, then the caller promises that
      ///   they won't actually access the object.
      ///
      /// If thisIsValid is false, implementations should not spend any
      /// effort getting the master local object.  This may save time
      /// on allocating temporary space, copying from device to host,
      /// etc.  This implies that implementations must be able to
      /// construct "null" / empty master local objects.
      this_type
      valid (const bool thisIsValid) const {
        return {G_, getExecutionSpace(), getMemorySpace(),
                thisIsValid};
      }

      /// \brief Declare intent to access this object's local data on
      ///   a specific (Kokkos) memory space (instance).
      template<class NewMemorySpace>
      typename std::conditional<
        is_execution_space_specified,
        LocalAccess<global_object_type, execution_space,
                    NewMemorySpace, access_mode>,
        LocalAccess<global_object_type, NewMemorySpace, access_mode>
        >::type
      on(const NewMemorySpace& memSpace) const {
        using Kokkos::Impl::is_memory_space;
        static_assert(is_memory_space<NewMemorySpace>::value,
          "NewMemorySpace must be a Kokkos memory space.");

        // We can't use std::conditional here, because we need the
        // "select" method.
        using alt_execution_space =
          typename LocalAccess<global_object_type, NewMemorySpace,
                               access_mode>::execution_space;
        auto execSpace = Kokkos::Impl::if_c<
          is_execution_space_specified,
          execution_space,
          alt_execution_space>::select(
            getExecutionSpace(),
            alt_execution_space());
        return {G_, execSpace, memSpace, isValid()};
      }

      /// \brief Declare intent to access this object's local data at
      ///   a specific (Kokkos) execution space (instance).
      template<class NewExecutionSpace>
      typename std::conditional<
        is_memory_space_specified,
        LocalAccess<global_object_type, NewExecutionSpace,
                    memory_space, access_mode>,
        LocalAccess<global_object_type, NewExecutionSpace, access_mode>
        >::type
      at(const NewExecutionSpace& execSpace) const {
        using Kokkos::Impl::is_execution_space;
        static_assert(is_execution_space<NewExecutionSpace>::value,
          "NewExecutionSpace must be a Kokkos execution space.");

        // We can't use std::conditional here, because we need the
        // "select" method.
        using alt_memory_space =
          typename LocalAccess<global_object_type, NewExecutionSpace,
                               access_mode>::memory_space;
        auto memSpace = Kokkos::Impl::if_c<
          is_memory_space_specified,
          memory_space,
          alt_memory_space>::select(
            getMemorySpace(),
            alt_memory_space());
        return {G_, execSpace, memSpace, isValid()};
      }

    public:
      /// \brief Reference to the global object whose data the user
      ///   will access.
      ///
      /// Keep by reference, because this struct is only valid in a
      /// delimited scope.
      global_object_type& G_;

      /// \brief Execution space instance, with which the user will
      ///   access local data.
      ///
      /// We assume that Kokkos execution spaces have shallow-copy
      /// semantics.
      execution_space execSpace_;

      /// \brief Memory space instance, on which the user will access
      ///   local data.
      ///
      /// We assume that Kokkos memory spaces have shallow-copy
      /// semantics.
      memory_space memSpace_;

      //! Will I actually need to access this object?
      bool valid_;
    };
  } // namespace Details

  ////////////////////////////////////////////////////////////
  // Implementations of readOnly, writeOnly, and readWrite
  ////////////////////////////////////////////////////////////

  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::read_only>
  readOnly (GlobalObjectType& G)
  {
    return {G};
  }

  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::read_only>
  readOnly (const GlobalObjectType& G)
  {
    GlobalObjectType& G_nc = const_cast<GlobalObjectType&> (G);
    return {G_nc};
  }

  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::write_only>
  writeOnly (GlobalObjectType& G)
  {
    return {G};
  }

  template<class GlobalObjectType>
  Details::LocalAccess<
    GlobalObjectType,
    Details::read_write>
  readWrite (GlobalObjectType& G)
  {
    return {G};
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
        typename ArgsToFunction<LocalAccessTypes...>::type;

      static void
      withLocalAccess (LocalAccessTypes...,
                       typename ArgsToFunction<LocalAccessTypes...>::type);
    };

    /// \brief Specialization of WithLocalAccess that implements the
    ///   "base class" of the user providing no GlobalObject
    ///   arguments, and a function that takes no arguments.
    template<>
    struct WithLocalAccess<> {
      using current_user_function_type =
        typename ArgsToFunction<>::type;

      static void
      withLocalAccess (current_user_function_type userFunction)
      {
        userFunction ();
      }
    };

    /// \brief Specialization of WithLocalAccess that implements the
    ///   "recursion case."
    ///
    /// \tparam FirstLocalAccessType Specialization of LocalAccess.
    /// \tparam Rest Zero or more possibly different specializations
    ///   of LocalAccess.
    template<class FirstLocalAccessType, class ... Rest>
    struct WithLocalAccess<FirstLocalAccessType, Rest...> {
      using current_user_function_type =
        typename ArgsToFunction<FirstLocalAccessType, Rest...>::type;

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
          ([=] (with_local_access_function_argument_type<Rest>... args) {
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

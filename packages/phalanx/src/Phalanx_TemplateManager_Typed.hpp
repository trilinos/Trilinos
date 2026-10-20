// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_TEMPLATE_MANAGER_TYPED_HPP
#define PHALANX_TEMPLATE_MANAGER_TYPED_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <typeinfo>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_TypeList.hpp"

namespace PHX {

  namespace detail {
    // A builder with a member function build<T>() returns an RCP to the base class.
    template <typename B, typename T>
    concept LegacyTemplateBuilder = requires (const B& b) { b.template build<T>(); };

    // A builder that is a generic lambda: []<typename T>() { return rcp(new Foo<T>); }
    template <typename B, typename T>
    concept TemplateBuilder = requires (const B& b) { b.template operator()<T>(); };

    template <typename L, template <typename> class ObjectT>
    struct object_tuple;

    template <typename... Ts, template <typename> class ObjectT>
    struct object_tuple<PHX::mpl::type_list<Ts...>, ObjectT> {
      using type = std::tuple<Teuchos::RCP<ObjectT<Ts>>...>;
    };
  }

  /** \brief Holds one ObjectT<T> for each type T in a type list.

      This is PHX::TemplateManager when Phalanx is configured with
      Phalanx_HIDE_DEPRECATED_CODE=ON. Include Phalanx_TemplateManager.hpp
      rather than this file. Objects are stored with their
      full static type, so get<T>() needs no cast, and asking for a type
      that is not in the list is a compile-time error.

      \tparam TypeSeq A PHX::mpl::type_list of types (usually Traits::EvalTypes).
      \tparam ObjectT The class template to hold, e.g. panzer::EquationSet.
                      Use an alias template for templates with more than
                      one parameter.
      \tparam BaseT   Optional common non-template base class of every
                      ObjectT<T>. When given, objects can also be accessed
                      through it with getAsBase<T>() and forEachBase().

      Objects are built for every enabled type with one of:
      \code
      // Same constructor arguments for every type:
      c.build<MyEquationSet>(params, order, cell_data, global_data, transient);

      // Any construction logic, as a generic lambda:
      c.buildObjects([&]<typename EvalT>() { return Teuchos::rcp(new Foo<EvalT>(...)); });

      // A builder with a build<T>() member:
      c.buildObjects(MyLegacyBuilder(...));
      \endcode

      <b>Migrating from the deprecated PHX::TemplateManager</b>

      | Deprecated                             | Replacement                        |
      |----------------------------------------|------------------------------------|
      | TemplateManager<Seq, Base, Foo<_>>     | TemplateManager<Seq, Foo, Base>    |
      | getAsObject<T>()                       | get<T>()                           |
      | disableType<T>()                       | disable<T>()                       |
      | deleteType<T>()                        | reset<T>()                         |
      | for (it = begin(); it != end(); ++it)  | forEach() or forEachBase()         |

      Differences from the deprecated class:
      - A builder's result is checked to be an ObjectT<T> when it is
        built, not when getAsObject() is first called.
      - ObjectT<T> must derive from BaseT. If the objects only share BaseT
        through a cross cast, use an alias template that ignores its
        argument as ObjectT, so that each slot holds a BaseT.
      - forEach() and forEachBase() skip types with no object (disabled,
        reset, or a builder that returned null). The iterators visited
        them.
  */
  template <typename TypeSeq, template <typename> class ObjectT, typename BaseT = void>
  class TemplateManager {
  public:
    using type_list = PHX::mpl::to_type_list_t<TypeSeq>;
    using base_type = BaseT;
    template <typename T> using object_type = ObjectT<T>;

    static constexpr std::size_t size = PHX::mpl::size_v<type_list>;

  private:
    typename detail::object_tuple<type_list,ObjectT>::type objects_;
    std::array<bool,size> disabled_{};

    template <typename T>
    static constexpr std::size_t index = PHX::mpl::index_v<T,type_list>;

    template <typename T>
    static void checkBase()
    {
      if constexpr (!std::is_void_v<BaseT>)
        static_assert(std::derived_from<ObjectT<T>,BaseT>,
                      "Every ObjectT<T> must derive from BaseT.");
    }

    // Kept out of the lambda in buildObjects(), where the builder's type
    // is not dependent and the discarded branch would still be checked.
    template <typename T, typename BuilderT>
    static Teuchos::RCP<ObjectT<T>> invokeBuilder(const BuilderT& builder)
    {
      if constexpr (detail::LegacyTemplateBuilder<BuilderT,T>) {
        auto obj = builder.template build<T>();
        using ReturnT = typename decltype(obj)::element_type;
        if constexpr (std::derived_from<ReturnT,ObjectT<T>>) {
          // Already an ObjectT<T>, no runtime check needed.
          return obj;
        }
        else {
          // Usually returns a base class pointer, so the type is checked at runtime.
          static_assert(std::derived_from<ObjectT<T>,ReturnT>,
                        "A builder's build<T>() must return an RCP to ObjectT<T>, a class derived from it, or one of its base classes.");
          auto typed = Teuchos::rcp_dynamic_cast<ObjectT<T>>(obj);
          TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::nonnull(obj) && Teuchos::is_null(typed), std::runtime_error,
            "PHX::TemplateManager::buildObjects(): the builder returned an object for type \""
            << typeid(T).name() << "\" that is not an ObjectT<T>.");
          return typed;
        }
      }
      else {
        static_assert(detail::TemplateBuilder<BuilderT,T>,
                      "A builder must provide build<T>() or be callable as builder.template operator()<T>().");
        return builder.template operator()<T>();
      }
    }

  public:

    TemplateManager() = default;

    /** \brief Builds an object for each enabled type.

        The builder is either a generic lambda (or other function
        object) called as `builder.template operator()<T>()`, or a
        class with a `build<T>()` member.
        Either may return null to leave a type empty.
    */
    template <typename BuilderT>
    void buildObjects(const BuilderT& builder)
    {
      PHX::mpl::for_each_type<type_list>([&]<typename T>() {
        checkBase<T>();
        if (disabled_[index<T>])
          return;
        std::get<index<T>>(objects_) = invokeBuilder<T>(builder);
      });
    }

    /** \brief Builds a Concrete<T>(args...) for each enabled type T.

        Concrete<T> must derive from ObjectT<T>. The arguments are
        passed by const reference to every constructor.
    */
    template <template <typename> class Concrete = ObjectT, typename... Args>
    void build(const Args&... args)
    {
      buildObjects([&]<typename T>() -> Teuchos::RCP<ObjectT<T>> {
        return Teuchos::rcp(new Concrete<T>(args...));
      });
    }

    //! Stores obj for type T, replacing any existing object.
    template <typename T>
    void set(const Teuchos::RCP<ObjectT<T>>& obj)
    { std::get<index<T>>(objects_) = obj; }

    //! Returns the object for type T (null if not built).
    template <typename T>
    Teuchos::RCP<ObjectT<T>> get()
    { return std::get<index<T>>(objects_); }

    //! Returns the object for type T (null if not built).
    template <typename T>
    Teuchos::RCP<const ObjectT<T>> get() const
    { return std::get<index<T>>(objects_); }

    //! Returns the object for type T as a BaseT (null if not built).
    template <typename T> requires (!std::is_void_v<BaseT>)
    Teuchos::RCP<BaseT> getAsBase()
    { return std::get<index<T>>(objects_); }

    //! Returns the object for type T as a BaseT (null if not built).
    template <typename T> requires (!std::is_void_v<BaseT>)
    Teuchos::RCP<const BaseT> getAsBase() const
    { return std::get<index<T>>(objects_); }

    //! Returns true if an object is stored for type T.
    template <typename T>
    bool has() const
    { return Teuchos::nonnull(std::get<index<T>>(objects_)); }

    //! Prevents later calls to build() and buildObjects() from building type T. Does not remove an existing object.
    template <typename T>
    void disable()
    { disabled_[index<T>] = true; }

    //! Returns false if disable<T>() was called.
    template <typename T>
    bool isEnabled() const
    { return !disabled_[index<T>]; }

    //! Removes the object for type T.
    template <typename T>
    void reset()
    { std::get<index<T>>(objects_) = Teuchos::null; }

    /** \brief Calls `f.template operator()<T>(obj)` with `ObjectT<T>& obj`
        for each type T that has an object, in type list order.
    */
    template <typename F>
    void forEach(F&& f)
    {
      PHX::mpl::for_each_type<type_list>([&]<typename T>() {
        if (auto& p = std::get<index<T>>(objects_); Teuchos::nonnull(p))
          f.template operator()<T>(*p);
      });
    }

    //! Const version of forEach(); obj is a `const ObjectT<T>&`.
    template <typename F>
    void forEach(F&& f) const
    {
      PHX::mpl::for_each_type<type_list>([&]<typename T>() {
        if (const auto& p = std::get<index<T>>(objects_); Teuchos::nonnull(p))
          f.template operator()<T>(static_cast<const ObjectT<T>&>(*p));
      });
    }

    //! Calls `f(BaseT&)` for each type that has an object, in type list order.
    template <typename F> requires (!std::is_void_v<BaseT>)
    void forEachBase(F&& f)
    {
      forEach([&]<typename T>(ObjectT<T>& obj) { f(static_cast<BaseT&>(obj)); });
    }

    //! Calls `f(const BaseT&)` for each type that has an object, in type list order.
    template <typename F> requires (!std::is_void_v<BaseT>)
    void forEachBase(F&& f) const
    {
      forEach([&]<typename T>(const ObjectT<T>& obj) { f(static_cast<const BaseT&>(obj)); });
    }
  };

}

#endif

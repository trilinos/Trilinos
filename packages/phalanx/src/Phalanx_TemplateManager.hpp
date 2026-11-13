// $Id$
// $Source$
// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_TEMPLATEMANAGER_HPP
#define PHX_TEMPLATEMANAGER_HPP

#include "Phalanx_config.hpp"

// The deprecated form below takes its object type as a Sacado placeholder
// expression and stores base class pointers. Configuring with
// Phalanx_HIDE_DEPRECATED_CODE=ON replaces it with the typed form.
// Code using get(), set(), has(), build(), buildObjects(), disable(),
// isEnabled(), reset(), getAsBase(), forEach() and forEachBase() compiles
// with either one.
#ifdef PHALANX_HIDE_DEPRECATED_CODE

#include "Phalanx_TemplateManager_Typed.hpp"

#else

#include <typeinfo>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "Sacado_mpl_apply.hpp"

#include "Phalanx_TypeList.hpp"

namespace PHX {

  // Forward declaration
  template <typename TypeSeq, typename BaseT, typename ObjectT>
  class TemplateIterator;
  template <typename TypeSeq, typename BaseT, typename ObjectT>
  class ConstTemplateIterator;

  //! Container class to manager template instantiations of a template class
  /*!
   * This class evaluatedFields a generic container class for managing multiple
   * instantiations of another class ObjectT.  It assumes each class
   * ObjectT<ScalarT> is derived from a non-template base class BaseT.  It
   * stores a vector of reference counted pointers to objects of type
   * BaseT corresponding to each instantiation of ObjectT.  The instantiations
   * ObjectT for each ScalarT are provided by a builder class, passed through
   * the buildObjects() method (see DefaultBuilderOp for an example builder).
   * An iterator is provided for accessing each template instantiation, and
   * non-templated virtual methods of BaseT can be called by dereferencing
   * the iterator.  Finally, template methods are provided to access the
   * stored objects as either objects of type BaseT (getAsBase()) or objects
   * of type ObjectT<ScalarT> (getAsObject()).  A map using RTTI is used to
   * map a typename to an index into the vector corresponding to the object
   * of that type.
   *
   * Template managers for specific types should derive from this class
   * instantiated on those types.  In most cases, a builder class will also
   * need to be created to instantiate the objects for each scalar type.
   */
  template <typename TypeSeq, typename BaseT, typename ObjectT>
  class TemplateManager {

    //! Implementation of < for type_info objects
    struct type_info_less {
      bool operator() (const std::type_info* a, const std::type_info* b) {
	return a->before(*b);
      }
    };

    template <typename BuilderOpT, typename T>
    static auto invokeBuilder(const BuilderOpT& builder, int)
      -> decltype(builder.template build<T>())
    { return builder.template build<T>(); }

    template <typename BuilderOpT, typename T>
    static auto invokeBuilder(const BuilderOpT& builder, long)
      -> decltype(builder.template operator()<T>())
    { return builder.template operator()<T>(); }

    template <typename BuilderOpT>
    struct BuildObject {
      std::vector< Teuchos::RCP<BaseT> >& objects;
      std::vector< bool >& disabled;
      const BuilderOpT& builder;
      BuildObject(std::vector< Teuchos::RCP<BaseT> >& objects_,
                  std::vector<bool>& disabled_,
		  const BuilderOpT& builder_) :
	objects(objects_),disabled(disabled_),builder(builder_) {}
      template <typename T> void operator()() const {
	constexpr std::size_t idx = PHX::mpl::index_v<T,TypeSeq>;
        if (!disabled[idx])
          objects[idx] = invokeBuilder<BuilderOpT,T>(builder,0);
      }
    };

    //! Declare TemplateIterator as a friend class
    friend class TemplateIterator<TypeSeq,BaseT,ObjectT>;

  public:

    //! Typedef for iterator
    typedef TemplateIterator<TypeSeq,BaseT,ObjectT> iterator;

    //! Typedef for const_iterator
    typedef ConstTemplateIterator<TypeSeq,BaseT,ObjectT> const_iterator;

    //! The default builder class for building objects for each ScalarT
    struct DefaultBuilderOp {

      //! Returns a new rcp for an object of type ObjectT<ScalarT>
      template<class ScalarT>
      Teuchos::RCP<BaseT> build() const {
	typedef typename Sacado::mpl::apply<ObjectT,ScalarT>::type type;
	return Teuchos::rcp( dynamic_cast<BaseT*>( new type ) );
      }

    };

    //! Default constructor
    TemplateManager();

    //! Destructor
    ~TemplateManager();

    /** \brief Build objects for each ScalarT

        The builder either has a member `build<T>()` or is a generic lambda
        called as `builder.template operator()<T>()`.
    */
    template <typename BuilderOpT>
    void buildObjects(const BuilderOpT& builder);

    //! Build objects for each ScalarT using default builder
    void buildObjects();

    //! Get RCP to object corrensponding to ScalarT as BaseT
    template<typename ScalarT>
    Teuchos::RCP<BaseT> getAsBase();

    //! Get RCP to object corrensponding to ScalarT as BaseT
    template<typename ScalarT>
    Teuchos::RCP<const BaseT> getAsBase() const;

    //! Get RCP to object corrensponding to ScalarT as ObjectT<ScalarT>
    template<typename ScalarT>
    PHALANX_DEPRECATED_MSG("Use get<T>(). See packages/phalanx/scripts/migrate_template_manager.py.")
    Teuchos::RCP< typename Sacado::mpl::apply<ObjectT,ScalarT>::type > getAsObject();

    //! Get RCP to object corrensponding to ScalarT as ObjectT<ScalarT>
    template<typename ScalarT>
    PHALANX_DEPRECATED_MSG("Use get<T>(). See packages/phalanx/scripts/migrate_template_manager.py.")
    Teuchos::RCP< const typename Sacado::mpl::apply<ObjectT,ScalarT>::type > getAsObject() const;

    //! Return an iterator that points to the first type object
    PHALANX_DEPRECATED_MSG("Use forEach() or forEachBase(). See packages/phalanx/scripts/migrate_template_manager.py.")
    typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator begin();

    //! Return an iterator that points to the first type object
    PHALANX_DEPRECATED_MSG("Use forEach() or forEachBase(). See packages/phalanx/scripts/migrate_template_manager.py.")
    typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::const_iterator
    begin() const;

    //! Return an iterator that points one past the last type object
    PHALANX_DEPRECATED_MSG("Use forEach() or forEachBase(). See packages/phalanx/scripts/migrate_template_manager.py.")
    typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::iterator end();

    //! Return an iterator that points one past the last type object
    PHALANX_DEPRECATED_MSG("Use forEach() or forEachBase(). See packages/phalanx/scripts/migrate_template_manager.py.")
    typename PHX::TemplateManager<TypeSeq,BaseT,ObjectT>::const_iterator
    end() const;

    //! Delete the underlying type. Used to clean out unused types.
    template <typename BuilderOpT>
    PHALANX_DEPRECATED_MSG("Use reset<T>(). See packages/phalanx/scripts/migrate_template_manager.py.")
    void deleteType();

    //! Disable the type so that it is not allocated.
    template <typename BuilderOpT>
    PHALANX_DEPRECATED_MSG("Use disable<T>(). See packages/phalanx/scripts/migrate_template_manager.py.")
    void disableType();

    //! Builds a Concrete<T>(args...) for each enabled type T.
    template <template <typename> class Concrete, typename... Args>
    void build(const Args&... args)
    {
      buildObjects([&]<typename T>() -> Teuchos::RCP<BaseT> {
        return Teuchos::rcp(new Concrete<T>(args...));
      });
    }

    //! Stores obj for type T, replacing any existing object.
    template <typename T>
    void set(const Teuchos::RCP<typename Sacado::mpl::apply<ObjectT,T>::type>& obj)
    { objects[PHX::mpl::index_v<T,TypeSeq>] = obj; }

    //! Returns the object for type T (null if not built).
    template <typename T>
    Teuchos::RCP<typename Sacado::mpl::apply<ObjectT,T>::type> get()
    {
      return Teuchos::rcp_dynamic_cast<typename Sacado::mpl::apply<ObjectT,T>::type>(
        objects[PHX::mpl::index_v<T,TypeSeq>], true);
    }

    //! Returns the object for type T (null if not built).
    template <typename T>
    Teuchos::RCP<const typename Sacado::mpl::apply<ObjectT,T>::type> get() const
    {
      return Teuchos::rcp_dynamic_cast<const typename Sacado::mpl::apply<ObjectT,T>::type>(
        objects[PHX::mpl::index_v<T,TypeSeq>], true);
    }

    //! Returns true if an object is stored for type T.
    template <typename T>
    bool has() const
    { return Teuchos::nonnull(objects[PHX::mpl::index_v<T,TypeSeq>]); }

    //! Prevents later calls to build() and buildObjects() from building type T. Does not remove an existing object.
    template <typename T>
    void disable()
    { disabled[PHX::mpl::index_v<T,TypeSeq>] = true; }

    //! Returns false if disable<T>() was called.
    template <typename T>
    bool isEnabled() const
    { return !disabled[PHX::mpl::index_v<T,TypeSeq>]; }

    //! Removes the object for type T.
    template <typename T>
    void reset()
    { objects[PHX::mpl::index_v<T,TypeSeq>] = Teuchos::null; }

    /** \brief Calls `f.template operator()<T>(obj)` with the object for
        each type T that has one, in type list order.
    */
    template <typename F>
    void forEach(F&& f)
    {
      PHX::mpl::for_each_type<TypeSeq>([&]<typename T>() {
        if (auto p = this->template get<T>(); Teuchos::nonnull(p))
          f.template operator()<T>(*p);
      });
    }

    //! Const version of forEach(); the object is passed as const.
    template <typename F>
    void forEach(F&& f) const
    {
      PHX::mpl::for_each_type<TypeSeq>([&]<typename T>() {
        if (auto p = this->template get<T>(); Teuchos::nonnull(p))
          f.template operator()<T>(*p);
      });
    }

    //! Calls `f(BaseT&)` for each type that has an object, in type list order.
    template <typename F>
    void forEachBase(F&& f)
    {
      for (auto& p : objects)
        if (Teuchos::nonnull(p))
          f(*p);
    }

    //! Calls `f(const BaseT&)` for each type that has an object, in type list order.
    template <typename F>
    void forEachBase(F&& f) const
    {
      for (const auto& p : objects)
        if (Teuchos::nonnull(p))
          f(static_cast<const BaseT&>(*p));
    }

  private:

    //! Stores array of rcp's to objects of each type
    std::vector< Teuchos::RCP<BaseT> > objects;

    //! Set to true if the type should not be allocated
    std::vector<bool> disabled;
  };

}

// Include template definitions
#include "Phalanx_TemplateManager_Def.hpp"

// Include template iterator definition
#include "Phalanx_TemplateIterator.hpp"

#endif // PHALANX_HIDE_DEPRECATED_CODE

#endif

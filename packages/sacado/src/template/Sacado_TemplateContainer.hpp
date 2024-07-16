// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TEMPLATE_CONTAINER_HPP
#define SACADO_TEMPLATE_CONTAINER_HPP

#include <type_traits>

// While this code does not directly use C++11 features, it uses mpl::vector,
// which does
#include "Sacado_ConfigDefs.h"

#include "Sacado_mpl_size.hpp"
#include "Sacado_mpl_find.hpp"
#include "Sacado_mpl_for_each.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_deref.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_enable_if.hpp"

namespace Sacado {

  namespace Impl {

    // Forward declaration
    template <typename TypeSeq,
              typename ObjectT,
              typename Iter1 = typename mpl::begin<TypeSeq>::type,
              typename Iter2 =typename  mpl::end<TypeSeq>::type>
    struct TupleSeq;

    // Forward declaration
    template <typename T,
              typename TypeSeq,
              typename ObjectT,
              typename Iter1 = typename mpl::begin<TypeSeq>::type,
              typename Iter2 = typename mpl::end<TypeSeq>::type,
              typename Enabled = void>
    struct GetTupleSeq;

  } // namespace Impl


  //! Container class to manager template instantiations of a template class
  /*!
   * This class provides a generic container class for managing multiple
   * instantiations of a class ObjectT<T> with the values of T specified
   * through a type sequence.  Objects of type ObjectT<T> must have value
   * semantics, be default constructable and copyable/assignable (for objects
   * that don't satisfy these requirements, one would typically wrap them in
   * a (smart) pointer.  Methods are provided to access the object of type
   * ObjectT<T> for a given type T, as well as initialize them.  One would
   * typically operate on the contained objects using mpl::for_each(), e.g.,
   *
   *    template <class T> class MyClass;
   *    typedef mpl::vector< double, Fad::DFad<double> > MyTypes;
   *    typedef TemplateContainer< MyTypes, MyClass<_> > MyObjects;
   *    MyObjects myObjects;
   *    mpl::for_each<MyObjects>([](auto x) {
   *       typedef decltype(x) T;
   *       myObjects.get<T>() = T;
   *    });
   *
   * (Note:  This requires generalized lambda's introduced in C++14.  For
   * C++11, provide a functor that implements the lambda body).
   */
  template <typename TypeSeq, typename ObjectT>
  class TemplateContainer {

    //! Our container for storing each object
    typedef Impl::TupleSeq<TypeSeq, ObjectT> tuple_type;

    //! Helper class for building objects in container
    template <typename BuilderOpT>
    struct BuildObject {
      tuple_type& objects;
      const BuilderOpT& builder;
      BuildObject(tuple_type& objects_,
                  const BuilderOpT& builder_) :
        objects(objects_), builder(builder_) {}
      template <typename T> void operator()(T x) const {
        Impl::GetTupleSeq<T,TypeSeq,ObjectT>::apply(objects) = builder(x);
      }
    };

  public:

    //! Type sequence containing template types
    typedef TypeSeq types;

    //! The default builder class for building objects for each ScalarT
    struct DefaultBuilderOp {

      //! Returns a new object of type ObjectT<ScalarT>
      template<class T>
      typename Sacado::mpl::apply<ObjectT,T>::type operator() (T) const {
        return typename Sacado::mpl::apply<ObjectT,T>::type();
      }

    };

    //! Default constructor
    TemplateContainer() {}

    //! Destructor
    ~TemplateContainer() {}

    //! Get object corrensponding ObjectT<T>
    template<typename T>
    typename Sacado::mpl::apply<ObjectT,T>::type& get() {
      return Impl::GetTupleSeq<T,TypeSeq,ObjectT>::apply(objects);
    }

    //! Get object corrensponding ObjectT<T>
    template<typename T>
    const typename Sacado::mpl::apply<ObjectT,T>::type& get() const  {
      return Impl::GetTupleSeq<T,TypeSeq,ObjectT>::apply(objects);
    }

    //! Build objects for each type T
    template <typename BuilderOpT = DefaultBuilderOp>
    void build(const BuilderOpT& builder) {
      Sacado::mpl::for_each_no_kokkos<TypeSeq>(
        BuildObject<BuilderOpT>(objects,builder));
    }

  private:

    //! Stores type of objects of each type
    tuple_type objects;

  };

  // Wrap call to mpl::for_each so you don't have to specify the container
  // or type sequence
  template <typename TypeSeq, typename ObjectT, typename FunctorT>
  void container_for_each(TemplateContainer<TypeSeq,ObjectT>& container,
                          const FunctorT& op) {
    typedef TemplateContainer<TypeSeq,ObjectT> Container;
    Sacado::mpl::for_each<Container> f(op);
  }

  // Wrap call to mpl::for_each so you don't have to specify the container
  // or type sequence
  template <typename TypeSeq, typename ObjectT, typename FunctorT>
  void container_for_each_no_kokkos(TemplateContainer<TypeSeq,ObjectT>& container,
                                    const FunctorT& op) {
    typedef TemplateContainer<TypeSeq,ObjectT> Container;
    Sacado::mpl::for_each_no_kokkos<Container> f(op);
  }

  namespace mpl {

    // Give TemplateContainer begin<> and end<> iterators for for_each

    template <typename TypeSeq, typename ObjectT>
    struct begin< TemplateContainer<TypeSeq,ObjectT> > {
      typedef typename begin<TypeSeq>::type type;
    };

    template <typename TypeSeq, typename ObjectT>
    struct end< TemplateContainer<TypeSeq,ObjectT> > {
      typedef typename end<TypeSeq>::type type;
    };

  }

  namespace Impl {

    // Container class to store our tuple sequence
    template <typename TypeSeq,
              typename ObjectT,
              typename Iter1,
              typename Iter2>
    struct TupleSeq :
      TupleSeq<TypeSeq, ObjectT, typename mpl::next<Iter1>::type, Iter2>
    {
      typedef typename mpl::apply<ObjectT,
                                  typename mpl::deref<Iter1>::type>::type type;
      type tail;
    };

    template <typename TypeSeq,
              typename ObjectT,
              typename Iter1>
    struct TupleSeq<TypeSeq, ObjectT, Iter1, Iter1> {};

    // Helper class to get a value out of the tuple sequence from a given type
    template <typename T,
              typename TypeSeq,
              typename ObjectT,
              typename Iter1,
              typename Iter2,
              typename Enabled>
    struct GetTupleSeq {};

    template <typename T,
              typename TypeSeq,
              typename ObjectT,
              typename Iter1,
              typename Iter2>
    struct GetTupleSeq< T,
                        TypeSeq,
                        ObjectT,
                        Iter1,
                        Iter2,
                        typename mpl::enable_if_c<
                          std::is_same< T, typename mpl::deref<Iter1>::type
                                        >::value
                          >::type > {
      static typename TupleSeq<TypeSeq,ObjectT,Iter1,Iter2>::type&
      apply(TupleSeq<TypeSeq,ObjectT,Iter1,Iter2>& t) {
        return t.tail;
      }

      static const typename TupleSeq<TypeSeq,ObjectT,Iter1,Iter2>::type&
      apply(const TupleSeq<TypeSeq,ObjectT,Iter1,Iter2>& t) {
        return t.tail;
      }
    };

    template <typename T,
              typename TypeSeq,
              typename ObjectT,
              typename Iter1,
              typename Iter2>
    struct GetTupleSeq< T,
                        TypeSeq,
                        ObjectT,
                        Iter1,
                        Iter2,
                        typename mpl::enable_if_c<
                          !std::is_same< T, typename mpl::deref<Iter1>::type
                                         >::value
                          >::type > :
      GetTupleSeq< T, TypeSeq, ObjectT, typename mpl::next<Iter1>::type, Iter2>
    {};

    template <typename T,
              typename TypeSeq,
              typename ObjectT,
              typename Iter1>
    struct GetTupleSeq< T, TypeSeq, ObjectT, Iter1, Iter1, void> {};


  } // namespace Impl

}

#endif

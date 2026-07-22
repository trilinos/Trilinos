// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TEMPLATEITERATOR_HPP
#define SACADO_TEMPLATEITERATOR_HPP

#include <vector>
#include <iterator>

namespace Sacado {

  /*! 
   * Iterator for traversing through template instantiations stored by
   * the TemplateManager class.
   */
  /*!
   * This class implements a standard forward iterator for the 
   * TemplateManager.
   */
  template <typename BaseT>
  class TemplateIterator : public std::iterator<std::input_iterator_tag,
                                                BaseT> {
  public:

    //! Constructor
    TemplateIterator(
	    typename std::vector< Teuchos::RCP<BaseT> >::iterator p) :
      object_iterator(p) {}
    
    // No default constructor
    // Use default copy constructor and copy assignment

    //! Equal operator
    bool operator==(const TemplateIterator& t) const {
      return object_iterator == t.objectIterator; 
    }

    //! Not equal operator
    bool operator!=(const TemplateIterator& t) const {
      return object_iterator != t.object_iterator; 
    }

    //! Dereference operator
    typename Sacado::TemplateIterator<BaseT>::reference 
    operator*() const {
      return *(*object_iterator);
    }

    //! -> operator
    typename Sacado::TemplateIterator<BaseT>::pointer 
    operator->() const {
      return &(*(*object_iterator));
    }

    //! Prefix ++
    TemplateIterator& operator++() {
      ++object_iterator;
      return *this;
    }

    //! Postfix ++
    TemplateIterator operator++(int) {
      TemplateIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    //! Returns a reference counted pointer object
    Teuchos::RCP<BaseT> rcp() const {
      return *object_iterator;
    }

  private:

    //! Underlying object iterator
    typename std::vector< Teuchos::RCP<BaseT> >::iterator object_iterator;
    
  };

  /*! 
   * Iterator for traversing through template instantiations stored by
   * the TemplateManager class.
   */
  /*!
   * This class implements a standard forward iterator for the 
   * TemplateManager.
   */
  template <typename BaseT>
  class ConstTemplateIterator : public std::iterator<std::input_iterator_tag,
                                                     BaseT> {
  public:

    //! Constructor
    ConstTemplateIterator(
       typename std::vector< Teuchos::RCP<BaseT> >::const_iterator p) :
      object_iterator(p) {}
    
    // No default constructor
    // Use default copy constructor and copy assignment

    //! Equal operator
    bool operator==(const ConstTemplateIterator& t) const {
      return object_iterator == t.objectIterator; 
    }

    //! Not equal operator
    bool operator!=(const ConstTemplateIterator& t) const {
      return object_iterator != t.object_iterator; 
    }

    //! Dereference operator
    typename Sacado::ConstTemplateIterator<BaseT>::reference 
    operator*() const {
      return *(*object_iterator);
    }

    //! -> operator
    typename Sacado::ConstTemplateIterator<BaseT>::pointer 
    operator->() const {
      return &(*(*object_iterator));
    }

    //! Prefix ++
    ConstTemplateIterator& operator++() {
      ++object_iterator;
      return *this;
    }

    //! Postfix ++
    ConstTemplateIterator operator++(int) {
      ConstTemplateIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    //! Returns a reference counted pointer object
    Teuchos::RCP<BaseT> rcp() const {
      return *object_iterator;
    }

  private:

    //! Underlying object iterator
    typename std::vector< Teuchos::RCP<BaseT> >::const_iterator object_iterator;
    
  };

}

#endif

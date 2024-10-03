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

#ifndef PHX_TEMPLATEITERATOR_HPP
#define PHX_TEMPLATEITERATOR_HPP

#include <iterator>
#include <vector>
#include "Phalanx_config.hpp"
#include "Phalanx_TemplateManager.hpp"

namespace PHX {

  /*!
   * Iterator for traversing through template instantiations stored by
   * the TemplateManager class.
   */
  /*!
   * This class implements a standard forward iterator for the
   * TemplateManager.
   */
  template <typename TypeSeq, typename BaseT, typename ObjectT>
  class TemplateIterator {
  public:

    // Iterator requirements
    using iterator_category = std::input_iterator_tag;
    using value_type = BaseT;
    using differnce_type = std::ptrdiff_t;
    using pointer = BaseT*;
    using reference = BaseT&;

    //! Constructor
    TemplateIterator(
	    PHX::TemplateManager<TypeSeq,BaseT,ObjectT>& m,
	    typename std::vector< Teuchos::RCP<BaseT> >::iterator p) :
      manager(&m), object_iterator(p) {}

    // No default constructor
    // Use default copy constructor and copy assignment

    //! Equal operator
    bool operator==(const TemplateIterator& t) const {
      return object_iterator == t.object_iterator;
    }

    //! Not equal operator
    bool operator!=(const TemplateIterator& t) const {
      return object_iterator != t.object_iterator;
    }

    //! Dereference operator
    typename PHX::TemplateIterator<TypeSeq, BaseT, ObjectT>::reference
    operator*() const {
      return *(*object_iterator);
    }

    //! -> operator
    typename PHX::TemplateIterator<TypeSeq, BaseT, ObjectT>::pointer
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

    //! Underlying template manager
    PHX::TemplateManager<TypeSeq,BaseT,ObjectT>* manager;

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
  template <typename TypeSeq, typename BaseT, typename ObjectT>
  class ConstTemplateIterator {
  public:

    // Iterator requirements
    using iterator_category = std::input_iterator_tag;
    using value_type = BaseT;
    using differnce_type = std::ptrdiff_t;
    using pointer = BaseT*;
    using reference = BaseT&;

    //! Constructor
    ConstTemplateIterator(
       const PHX::TemplateManager<TypeSeq,BaseT,ObjectT>& m,
       typename std::vector< Teuchos::RCP<BaseT> >::const_iterator p) :
      manager(&m), object_iterator(p) {}

    // No default constructor
    // Use default copy constructor and copy assignment

    //! Equal operator
    bool operator==(const ConstTemplateIterator& t) const {
      return object_iterator == t.object_iterator;
    }

    //! Not equal operator
    bool operator!=(const ConstTemplateIterator& t) const {
      return object_iterator != t.object_iterator;
    }

    //! Dereference operator
    typename PHX::ConstTemplateIterator<TypeSeq, BaseT, ObjectT>::reference
    operator*() const {
      return *(*object_iterator);
    }

    //! -> operator
    typename PHX::ConstTemplateIterator<TypeSeq, BaseT, ObjectT>::pointer
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

    //! Underlying template manager
    const PHX::TemplateManager<TypeSeq,BaseT,ObjectT>* manager;

    //! Underlying object iterator
    typename std::vector< Teuchos::RCP<BaseT> >::const_iterator object_iterator;

  };

}

#endif

// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
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
    const typename Sacado::ConstTemplateIterator<BaseT>::reference 
    operator*() const {
      return *(*object_iterator);
    }

    //! -> operator
    const typename Sacado::ConstTemplateIterator<BaseT>::pointer 
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

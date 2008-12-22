// $Id$ 
// $Source$ 
// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_TEMPLATEITERATOR_HPP
#define PHX_TEMPLATEITERATOR_HPP

#include <iterator>
#include <vector>
#include "Phalanx_ConfigDefs.hpp"
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
  class TemplateIterator : public std::iterator<std::input_iterator_tag,
						BaseT> {
  public:

    //! Constructor
    TemplateIterator(
	    PHX::TemplateManager<TypeSeq,BaseT,ObjectT>& m,
	    typename std::vector< Teuchos::RCP<BaseT> >::iterator p) :
      manager(&m), object_iterator(p) {}
    
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
  class ConstTemplateIterator : public std::iterator<std::input_iterator_tag,
						     BaseT> {
  public:

    //! Constructor
    ConstTemplateIterator(
       const PHX::TemplateManager<TypeSeq,BaseT,ObjectT>& m,
       typename std::vector< Teuchos::RCP<BaseT> >::const_iterator p) :
      manager(&m), object_iterator(p) {}
    
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
    const typename PHX::ConstTemplateIterator<TypeSeq, BaseT, ObjectT>::reference 
    operator*() const {
      return *(*object_iterator);
    }

    //! -> operator
    const typename PHX::ConstTemplateIterator<TypeSeq, BaseT, ObjectT>::pointer 
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

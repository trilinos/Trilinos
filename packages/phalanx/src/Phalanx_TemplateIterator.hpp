// $Id$ 
// $Source$ 
// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

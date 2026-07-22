// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_VectorTemplateIterator_hpp__
#define __Panzer_VectorTemplateIterator_hpp__

#include <iterator>
#include <vector>

#include "Teuchos_RCP.hpp"

namespace panzer {

/*! 
 * Iterator for traversing through template instantiations stored by
 * the TemplateManager class.
 */
/*!
 * This class implements a standard forward iterator for the 
 * TemplateManager.
 */
template <typename TypeSeq, typename BaseT, typename ObjectT>
class VectorTemplateIterator {
public:

  // Iterator requirements
  using iterator_category = std::input_iterator_tag;
  using value_type = std::vector<Teuchos::RCP<BaseT>>;
  using differnce_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;
  
   //! Constructor
   VectorTemplateIterator(panzer::VectorTemplateManager<TypeSeq,BaseT,ObjectT> & m,
                    typename std::vector<std::vector< Teuchos::RCP<BaseT> > >::iterator p) 
      : manager_(&m), object_iterator_(p) {}

   //! Equal operator
   bool operator==(const VectorTemplateIterator& t) const 
   { return object_iterator_ == t.objectIterator_; }

   //! Not equal operator
   bool operator!=(const VectorTemplateIterator& t) const 
   { return object_iterator_ != t.object_iterator_; }

   //! Dereference operator
   typename panzer::VectorTemplateIterator<TypeSeq, BaseT, ObjectT>::reference
   operator*() const {
     return *(object_iterator_);
   }

   //! -> operator
   typename panzer::VectorTemplateIterator<TypeSeq, BaseT, ObjectT>::pointer
   operator->() const {
     return &(*(*object_iterator_));
   }

   //! Prefix ++
   VectorTemplateIterator& operator++() {
     ++object_iterator_;
     return *this;
   }

   //! Postfix ++
   VectorTemplateIterator operator++(int) {
     VectorTemplateIterator tmp = *this;
     ++(*this);
     return tmp;
   }

private:
   panzer::VectorTemplateManager<TypeSeq,BaseT,ObjectT> * manager_;
   typename std::vector<std::vector< Teuchos::RCP<BaseT> > >::iterator object_iterator_;
};

template <typename TypeSeq, typename BaseT, typename ObjectT>
class ConstVectorTemplateIterator {
public:

  // Iterator requirements
  using iterator_category = std::input_iterator_tag;
  using value_type = const std::vector<Teuchos::RCP<BaseT>>;
  using differnce_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  //! Constructor
   ConstVectorTemplateIterator(const panzer::VectorTemplateManager<TypeSeq,BaseT,ObjectT> & m,
                    typename std::vector<std::vector< Teuchos::RCP<BaseT> > >::const_iterator p) 
      : manager_(&m), object_iterator_(p) {}

   //! Equal operator
   bool operator==(const ConstVectorTemplateIterator& t) const 
   { return object_iterator_ == t.objectIterator_; }

   //! Not equal operator
   bool operator!=(const ConstVectorTemplateIterator& t) const 
   { return object_iterator_ != t.object_iterator_; }

   //! Dereference operator
   typename panzer::ConstVectorTemplateIterator<TypeSeq, BaseT, ObjectT>::reference
   operator*() const {
     return *(object_iterator_);
   }

   //! -> operator
   typename panzer::ConstVectorTemplateIterator<TypeSeq, BaseT, ObjectT>::pointer
   operator->() const {
     return &(*(*object_iterator_));
   }

   //! Prefix ++
   ConstVectorTemplateIterator& operator++() {
     ++object_iterator_;
     return *this;
   }

   //! Postfix ++
   ConstVectorTemplateIterator operator++(int) {
     ConstVectorTemplateIterator tmp = *this;
     ++(*this);
     return tmp;
   }

private:
   const panzer::VectorTemplateManager<TypeSeq,BaseT,ObjectT> * manager_;
   typename std::vector<std::vector< Teuchos::RCP<BaseT> > >::const_iterator object_iterator_;
};

}

#endif

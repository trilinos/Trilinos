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
class VectorTemplateIterator : public std::iterator<std::input_iterator_tag,
                                              std::vector<Teuchos::RCP<BaseT> > > {
public:
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
class ConstVectorTemplateIterator : public std::iterator<std::input_iterator_tag,
                                                         const std::vector<Teuchos::RCP<BaseT> > > {
public:
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

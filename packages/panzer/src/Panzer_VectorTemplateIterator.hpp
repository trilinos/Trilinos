// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

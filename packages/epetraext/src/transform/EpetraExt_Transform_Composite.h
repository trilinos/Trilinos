//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

//-----------------------------------------------------------------------
// EpetraExt_Transform_Composite.h
//-----------------------------------------------------------------------

#ifndef EPETRAEXT_TRANSFORM_COMPOSITE_H
#define EPETRAEXT_TRANSFORM_COMPOSITE_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <EpetraExt_Transform.h>

#include <list>

namespace EpetraExt {

//! Composition Class for Epetra Transform SameType Operators.
/*! This class allows <tt>SameType</tt> Transforms to be composed as
  a single Transform.
 */
template<typename T>
class Transform_Composite : public SameTypeTransform<T>
{

 public:

  typedef SameTypeTransform<T> * TransformTypePtr;

  //! EpetraExt::Transform_Composite Constructor
  Transform_Composite() {}

  //! EpetraExt::Transform_Composite Destructor
  virtual ~Transform_Composite();

  //! Transform Addition
  /*! Add <tt>SameType</tt> Transform to composition.
      Order of Addition == Order of Application
    */
  void addTransform( TransformTypePtr new_trans );

  //! Analysis phase generates plan and check feasibility
  /*! Analysis of transform operation on original object and construction
      of new object.
     
      \return Returns a pointer to the newly created object of type
      NewTypePtr.  The Transform object maintains ownership of this
      new object and deletes as a part of it's destruction.
    */
  virtual
  typename Transform<T,T>::NewTypeRef
  operator()
  ( typename Transform<T,T>::OriginalTypeRef orig );

  //! Forward Data Transfer
  /*! Forward transfer of data from <tt>orig</tt> object input in the
      <tt>operator()</tt> method call to the new object created in this
      same call.  Returns <tt>true</tt> is
      operation is successful.
    */
  virtual bool fwd();

  //!
  /*! Reverse transfer of data from new object created in the
      <tt>operator()</tt> method call to the <tt>orig</tt> object input
      to this same method. Returns <tt>true</tt> if operation is successful.
    */
  virtual bool rvs();

 protected:

  typedef typename std::list<TransformTypePtr>::iterator         TransListIter;
  typedef typename std::list<TransformTypePtr>::reverse_iterator TransListRvsIter;

  std::list<TransformTypePtr> transList_;

}; // end class Tranform_Composite

template<typename T>
Transform_Composite<T>::
~Transform_Composite()
{
  TransListIter iter = transList_.begin();
  TransListIter end = transList_.end();
  for( ; iter != end; ++iter ) delete *iter;
}

template<typename T>
void
Transform_Composite<T>::
addTransform( TransformTypePtr new_trans )
{
  transList_.push_back( new_trans );
}

template<typename T>
typename Transform<T,T>::NewTypeRef
Transform_Composite<T>::
operator()
( typename Transform<T,T>::OriginalTypeRef orig )
{
  this->origObj_ = &orig;
  this->newObj_ = &orig;

  TransListIter iter = transList_.begin();
  TransListIter end = transList_.end();
  for( ; iter != end; ++iter )
    this->newObj_ = &((**iter)( *(this->newObj_) ));

  return *(this->newObj_);
}

template<typename T>
bool
Transform_Composite<T>::
fwd()
{
  bool success = true;

  TransListIter iter = transList_.begin();
  TransListIter end = transList_.end();
  for( ; iter != end; ++iter )
    if( !(**iter).fwd() ) return false;

  return success;
}

template<typename T>
bool
Transform_Composite<T>::
rvs()
{
  bool success = true;

  TransListRvsIter iter = transList_.rbegin();
  TransListRvsIter end = transList_.rend();
  for( ; iter != end; ++iter )
    if( !(**iter).rvs() ) return false;

  return success;
}

} //namespace EpetraExt
  
#endif //EPETRAEXT_TRANSFORM_COMPOSITE_H

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
// EpetraExt_Transform.h
//-----------------------------------------------------------------------

#ifndef EPETRAEXT_TRANSFORM_H
#define EPETRAEXT_TRANSFORM_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <EpetraExt_ConfigDefs.h>

#include <Teuchos_RCP.hpp>

namespace EpetraExt {

///
/** Base Class for all Epetra Transform Operators.
 *
 * This is the abstract definition for all Epetra Transform Operators.
 * Depending on the type of Transform, several specializations are
 * available: Structural, SameType, InPlace, View.
 */

template<typename T, typename U>
class Transform
{
 public:

  /** @name Typedefs for templated classes */
  //@{

  typedef T  OriginalType;
  typedef T* OriginalTypePtr;
  typedef Teuchos::RCP<T> OriginalTypeRCP;
  typedef T& OriginalTypeRef;

  typedef U  NewType;
  typedef U* NewTypePtr;
  typedef Teuchos::RCP<U> NewTypeRCP;
  typedef U& NewTypeRef;

  //@}

  ///
  virtual ~Transform() {}

  /** @name Pure Virtual Methods which must be implemented by subclasses */
  //@{

  ///
  /** Analysis of transform operation on original object and construction
    * of new object.
    *
    * Preconditions:<ul>
    * </ul>
    *
    * Invariants:<ul>
    * </ul>
    *
    * Postconditions:<ul>
    * </ul>
    *
    * @return Returns a pointer to the newly created object of type
    * NewTypeRef.  The Transform object maintains ownership of this
    * new object and deletes as a part of it's destruction.
    */
  virtual NewTypeRef operator()( OriginalTypeRef orig ) = 0;

  ///
  /** Forward transfer of data from <tt>orig</tt> object input in the
    * <tt>operator()</tt> method call to the new object created in this
    * same call.  Returns <tt>true</tt> is
    * operation is successful.
    *
    * Preconditions:<ul>
    * </ul>
    *
    * Invariants:<ul>
    * </ul>
    *
    * Postconditions:<ul>
    * </ul>
    *
    */
  virtual bool fwd() = 0;

  ///
  /** Reverse transfer of data from new object created in the
    * <tt>operator()</tt> method call to the <tt>orig</tt> object input
    * to this same method. Returns <tt>true</tt> if operation is successful.
    *
    * Preconditions:<ul>
    * </ul>
    *
    * Invariants:<ul>
    * </ul>
    *
    * Postconditions:<ul>
    * </ul>
    *
    */
  virtual bool rvs() = 0;

  //@}

  /** @name Virtual functions with default implements allowing for optional
    * implementation by the Transform developer */
  //@{

  ///
  /** Initial analysis phase of transform.  Returns <tt>true</tt> if
    * the transform is possible allowing methods <tt>construct()</tt>,
    * <tt>fwd()</tt> and <tt>rvs()</tt> to be successfully utilized.
    *
    * Preconditions:<ul>
    * </ul>
    *
    * Invariants:<ul>
    * </ul>
    *
    * Postconditions:<ul>
    * </ul>
    *
    * The default implementation calls method <tt>operator()</tt>
    * and stores the resulting object in an internal attribute <tt>newObj_</tt>.
    */
  virtual bool analyze( OriginalTypeRef orig );

  ///
  /** Construction of new object as a result of the transform.
    *
    * Preconditions:<ul>
    * </ul>
    *
    * Invariants:<ul>
    * </ul>
    *
    * Postconditions:<ul>
    * </ul>
    *
    * The default implementation returns internal attribute <tt>newObj_</tt>.
    */
  virtual NewTypeRef construct();

  ///
  /** Check for whether transformed object has been constructed
    *
    * Preconditions:<ul>
    * </ul>
    *
    * Invariants:<ul>
    * </ul>
    *
    * Postconditions:<ul>
    * </ul>
    *
    * The default implementation returns <tt>true</tt> if
    * <tt>newObj_</tt> != 0.
    */
  virtual bool isConstructed();

  //@}

 protected:

  ///
  /** Default constructor, protected to allow only derived classes to use.
    *
    * Initializes attributes <tt>origObj_</tt> and <tt>newObj_</tt> to 0.
    *
    */
  Transform()
  : origObj_(0),
    newObj_(0)
  {}

  OriginalTypePtr origObj_;

  NewTypePtr      newObj_;

 private:
  Transform(const Transform<T,U>& src)
    :origObj_(src.origObj_), newObj_(src.newObj_) {}

  Transform<T,U>& operator=(const Transform<T,U>& src)
    {
      //not currently supported
      abort();
      return(*this);
    }

}; // end class Transform

template<typename T,typename U>
bool
Transform<T,U>::
analyze( OriginalTypeRef orig )
{
  origObj_ = &orig;
  newObj_ = &((*this)( *origObj_ ));
  return true;
}

template<typename T,typename U>
typename Transform<T,U>::NewTypeRef
Transform<T,U>::
construct()
{
  return *newObj_;
}

template<typename T,typename U>
bool
Transform<T,U>::
isConstructed()
{
  return ( newObj_ != 0 );
}

template<typename T, typename U>
class StructuralTransform : public Transform<T,U>
{
 public:
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~StructuralTransform() {}
};

template<typename T>
class SameTypeTransform : public Transform<T,T>
{
 public:
  typedef T  TransformType;
  typedef T* TransformTypePtr;
  typedef T& TransformTypeRef;

  virtual ~SameTypeTransform() {}
};

template<typename T>
class StructuralSameTypeTransform : public SameTypeTransform<T>
{
 public:
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~StructuralSameTypeTransform() {}
};

template<typename T>
class InPlaceTransform : public SameTypeTransform<T>
{
 public:
  typename Transform<T,T>::NewTypeRef
  operator()
  ( typename Transform<T,T>::OriginalTypeRef orig )
  { this->origObj_ = &orig;
    this->newObj_ = &orig;
    return orig;
  }

  virtual ~InPlaceTransform() {}
};

template<typename T>
class ViewTransform : public SameTypeTransform<T>
{
 public:
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~ViewTransform() {}
};

} //namespace EpetraExt
  
#endif //EPETRAEXT_TRANSFORM_H

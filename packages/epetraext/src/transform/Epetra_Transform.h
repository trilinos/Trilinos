//-----------------------------------------------------------------------
// Epetra_Transform.h
//-----------------------------------------------------------------------

#ifndef EPETRA_TRANSFORM_H
#define EPETRA_TRANSFORM_H

#ifdef HAVE_CONFIG_H

#undef PACKAGE

#include <EpetraExt_config.h>

#endif

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
  typedef T& OriginalTypeRef;

  typedef U  NewType;
  typedef U* NewTypePtr;
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

}; // end class Tranform

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
struct StructuralTransform : public Transform<T,U>
{
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~StructuralTransform() {}
};

template<class T>
struct SameTypeTransform : public Transform<T,T>
{
  typedef T  TransformType;
  typedef T* TransformTypePtr;
  typedef T& TransformTypeRef;

  virtual ~SameTypeTransform() {}
};

template<typename T>
struct StructuralSameTypeTransform : public SameTypeTransform<T>
{
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~StructuralSameTypeTransform() {}
};

template<class T>
struct InPlaceTransform : public SameTypeTransform<T>
{
  typename Transform<T,T>::NewTypeRef
  operator()
  ( typename Transform<T,T>::OriginalTypeRef orig )
  { origObj_ = &orig;
    newObj_ = &orig;
    return orig;
  }

  virtual ~InPlaceTransform() {}
};

template<class T>
struct ViewTransform : public SameTypeTransform<T>
{
  bool fwd() { return true; }
  bool rvs() { return true; }

  virtual ~ViewTransform() {}
};

} //namespace EpetraExt
  
#endif //EPETRA_TRANSFORM_H

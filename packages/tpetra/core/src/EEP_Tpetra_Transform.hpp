#ifndef TPETRA_TRANSFORM_H
#define TPETRA_TRANSFORM_H

#include <Teuchos_RCP.hpp>

namespace Tpetra {

///
/** Base Class for all Tpetra Transform Operators.
 *
 * This is the abstract definition for all Tpetra Transform Operators.
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
    :origObj_(src.origObj_), newObj_(src.newObj_)
  {
    std::cout << "Passing through tpetra/code/src/Tpetra_Transform.h Transform::copy_constructor()..." << std::endl;
  }

  Transform<T,U>& operator=(const Transform<T,U>& src)
    {
      std::cout << "Passing through tpetra/code/src/Tpetra_Transform.h Transform::assingment_opeeator()..." << std::endl;
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

} //namespace Tpetra

#endif

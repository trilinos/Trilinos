//-----------------------------------------------------------------------
// EpetraExt_Transform_Composite.h
//-----------------------------------------------------------------------

#ifndef EPETRAEXT_TRANSFORM_COMPOSITE_H
#define EPETRAEXT_TRANSFORM_COMPOSITE_H

#include <EpetraExt_Transform.h>

#include <list>

namespace EpetraExt {

///
/** Composition Class for Epetra Transform SameType Operators.
 *
 * This class allows <tt>SameType</tt> Transforms to be composed as
 * a single Transform.
 */

template<typename T>
class Transform_Composite : public SameTypeTransform<T>
{

 public:

  typedef SameTypeTransform<T> * TransformTypePtr;

  ///
  /** Default Constructor
    */
  Transform_Composite() {}

  ///
  /** Destructor
    */
  virtual ~Transform_Composite();

  ///
  /** Add <tt>SameType</tt> Transform to composition.
    * Order of Addition == Order of Application
    */
  void addTransform( TransformTypePtr new_trans );

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
    * NewTypePtr.  The Transform object maintains ownership of this
    * new object and deletes as a part of it's destruction.
    */
  virtual
  typename Transform<T,T>::NewTypeRef
  operator()
  ( typename Transform<T,T>::OriginalTypeRef orig );

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
  virtual bool fwd();

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
  virtual bool rvs();

 protected:

  typedef typename list<TransformTypePtr>::iterator         TransListIter;
  typedef typename list<TransformTypePtr>::reverse_iterator TransListRvsIter;

  list<TransformTypePtr> transList_;

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
  origObj_ = &orig;
  newObj_ = &orig;

  TransListIter iter = transList_.begin();
  TransListIter end = transList_.end();
  for( ; iter != end; ++iter )
    newObj_ = &((**iter)( *newObj_ ));

  return *newObj_;
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

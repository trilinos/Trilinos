//-----------------------------------------------------------------------
// Epetra_Transform_Composite.h
//-----------------------------------------------------------------------

#ifndef EPETRA_TRANSFORM_COMPOSITE_H
#define EPETRA_TRANSFORM_COMPOSITE_H

#include <Epetra_Transform.h>

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
  virtual NewTypeRef operator()( OriginalTypeRef orig );

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

  typedef std::list<TransformTypePtr>::iterator         TransListIter;
  typedef std::list<TransformTypePtr>::reverse_iterator TransListRvsIter;

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
Transform_Composite<T>::NewTypeRef
Transform_Composite<T>::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  NewTypeRef intermediateObj = *origObj_;

  TransListIter iter = transList_.begin();
  TransListIter end = transList_.end();
  for( ; iter != end; ++iter )
    intermediateObj = (**iter)( intermediateObj );

  newObj_ = &intermediateObj;

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
  
#endif //EPETRA_TRANSFORM_COMPOSITE_H

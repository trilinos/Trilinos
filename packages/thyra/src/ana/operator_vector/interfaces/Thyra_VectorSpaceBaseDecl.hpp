// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_SPACE_BASE_DECL_HPP
#define THYRA_VECTOR_SPACE_BASE_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_Describable.hpp"

namespace Thyra {


/** \defgroup Thyra_Op_Vec_createMember_grp User callable creational functions for VectorBase and MultiVectorBase.

These functions must be used by clients to create all <tt>VectorBase</tt> and
<tt>MultiVectorBase</tt> objects in order to ensure that the created objects
can live past the life time of the <tt>VectorSpaceBase</tt> object that
created them.

\ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp

*/

/** \brief Helper function that clones a <tt>VectorSpaceBase</tt> object
 * if the <tt>RefCountPtr</tt> does not have ownership.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
makeHaveOwnership( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs );

/** \brief Create a vector member from the vector space.
 *
 * Calls <tt>VectorSpaceBase::createMember()</tt> on <tt>vs</tt> but
 * the returned <tt>VectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr< VectorBase<Scalar> >
createMember( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs );

/** \brief Calls <tt>createMember(Teuchos::rcp(&vs,false))</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr< VectorBase<Scalar> >
createMember( const VectorSpaceBase<Scalar> &vs );

/** \brief Create a set of vector members (a <tt>MultiVectorBase</tt>) from the vector space.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr< MultiVectorBase<Scalar> >
createMembers( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, int numMembers );

/** \brief  Calls <tt>createMembers(Teuchos::rcp(&vs,false),numMembers)</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr< MultiVectorBase<Scalar> >
createMembers( const VectorSpaceBase<Scalar> &vs, int numMembers );

/** \brief Create a vector member that is a non-<tt>const</tt> view of raw data.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
createMemberView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::MutableSubVectorT<Scalar> &raw_v );

/** \brief Calls <tt>createMemberView(Teuchos::rcp(&vs,false),raw_v)</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
createMemberView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::MutableSubVectorT<Scalar> &raw_v );

/** \brief Create a vector member that is a <tt>const</tt> view of raw data.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
createMemberView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::SubVectorT<Scalar> &raw_v );

/** \brief Calls <tt>createMembersView(Teuchos::rcp(&vs,false),raw_v)</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const VectorBase<Scalar> >
createMemberView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::SubVectorT<Scalar> &raw_v );

/** \brief Create a multi-vector member that is a non-<tt>const</tt> view of raw data.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
createMembersView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv );

/** \brief Calls <tt>createMembersView(Teuchos::rcp(&vs,false),raw_mv)</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
createMembersView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv );

/** \brief Create a multi-vector member that is a <tt>const</tt> view of raw data.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
createMembersView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::SubMultiVectorT<Scalar> &raw_mv );

/** \brief Calls <tt>createMembersView(Teuchos::rcp(&vs,false),raw_mv)</tt>.
 *
 * \ingroup Thyra_Op_Vec_createMember_grp
 */
template<class Scalar>
Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
createMembersView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::SubMultiVectorT<Scalar> &raw_mv );


/** \brief Abstract interface for objects that represent a space for vectors.
 *
 * This interface acts primarily as an "Abstract Factory" interface for
 * creating <tt>VectorBase</tt> objects using the <tt>createMember()</tt>
 * method.  A <tt>%VectorSpaceBase</tt> can also create <tt>
 * MultiVectorBase</tt> objects which represent a compact collection of
 * vectors.  A secondary role for <tt>%VectorSpaceBase</tt> objects is to test
 * for compatibility of vector spaces (and the vectors and linear operators
 * using those spaces) objects using the <tt>isCompatible()</tt> method.
 *
 * Clients can not directly create <tt>%VectorBase</tt> and
 * <tt>%MultiVectorBase</tt> objects using the member functions
 * <tt>createMember()</tt>, <tt>createMembers()</tt>,
 * <tt>createMemberView()</tt>, and <tt>createMembersView()</tt> but instead
 * must use the \ref Thyra_Op_Vec_createMember_grp.
 *
 * A <tt>%VectorSpaceBase</tt> object can exist independent from any
 * individual <tt>VectorBase</tt> (or <tt>MultiVector</tt>) object; Or, a
 * <tt>%VectorSpaceBase</tt> object can have a lifetime that is dependent on a
 * single <tt>VectorBase</tt> ( or <tt>MultiVectorBase</tt>) object.  The same
 * interface serves both roles.
 *
 * Note that a <tt>%VectorSpaceBase</tt> object must also be able to create
 * <tt>MultiVectorBase</tt> objects with any number of column vectors using
 * the <tt>createMembers()</tt> function, and <tt>LinearOpBase::domain()</tt>
 * of these returned <tt>%MultiVectorBase</tt> objects returns a vector space
 * of that dimension.  An interesting side effect of this design is that the
 * creation of a multi-vector provides a way for clients to create vector
 * spaces of any arbitrary (although small usually) dimension.  In order to
 * give the client the same ability without having to create a full
 * multi-vector object first, the method <tt>smallVecSpcFcty()</tt> is
 * included.  The method <tt>smallVecSpcFcty()</tt> returns a
 * <tt>VectorSpaceFactoryBase</tt> object that can create (serial) vector
 * spaces of any small dimension that are compatible with the domain spaces of
 * so created <tt>%MultiVectorBase</tt> objects.
 *
 * A vector space is also where the scalar product for the space is defined
 * which is computed by the <tt>scalarProd()</tt> method.  A scalar product
 * allows the vector space to introduce scaling into many different types of
 * numerical algorithms.
 *
 * If the underlying object is not initialized, then <tt>dim()==0</tt> will be
 * true and none of the other methods should be called or exceptions will be
 * thrown.
 *
 * <b>Notes to subclass developers</b>
 *
 * A subclass is only required to override three methods:
 * <tt>dim()</tt>, <tt>isCompatible()</tt> and
 * <tt>createMember()</tt>.  Note that implementing the
 * <tt>createMember()</tt> method also entails defining a concrete
 * <tt>VectorBase</tt> subclass.
 *
 * If a subclass can support specialized multi-vectors, then the
 * <tt>createMembers()</tt> should be overridden as well.  Note that
 * implementing the <tt>createMembers()</tt> also entails defining a
 * concrete <tt>MultiVectorBase</tt> subclass.  For some types of concrete
 * <tt>MultiVectorBase</tt> subclass implementations (e.g. serial
 * multi-vectors), the same default <tt>VectorSpaceFactoryBase</tt> typed
 * object returned from the default <tt>smallVecSpcFcty()</tt> method
 * can be used.  However, more specialized <tt>MultiVectorBase</tt>
 * subclasses (e.g. distributed memory parallel) will require an
 * override of the <tt>smallVecSpcFcty()</tt> method to return a
 * specialized type of <tt>VectorSpaceFactoryBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_fundamental_interfaces_code_grp
 */
template<class Scalar>
class VectorSpaceBase : virtual public Teuchos::Describable {
public:

  /** @name Public pure virtual functions that must be overridden */
  //@{

  /** \brief Return the dimension of the vector space.
   *
   * If the underlying object is not initialized, then <tt>dim()==0</tt>
   * will be true and none of the other methods should be called.
   */
  virtual Index dim() const = 0;

  /** \brief Compare the compatibility of two vector spaces.
   *
   * If this function returns <tt>true</tt>, then vectors created
   * from either of the vector spaces will be compatible and can be
   * combined in vector operations.
   *
   * Invariants:<ul>
   * <li> [<tt>this->isCompatible(vecSpc) == true</tt>] <tt>vecSpc.isCompatible(*this) == true</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>this->dim() != vecSpc.dim()</tt>] <tt>return == false</tt>
   * </ul>
   */
  virtual bool isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const = 0;

  /** \brief Return a <tt>VectorSpaceFactoryBase</tt> object for the creation
   * of (serial) vector spaces with a small dimension.
   */
  virtual Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const = 0;

  /** \brief Return the scalar product of two vectors in the vector space.
   *
   * Preconditions:<ul>
   * <li><tt>x.space()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li><tt>y.space()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li>The scalar product is returned.
   * </ul>
   */
  virtual Scalar scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const = 0;

  /** \brief Return the scalar product of each column in two multi-vectors in the vector space.
   *
   * @param  X            [in] Multi-vector.
   * @param  Y            [in] Multi-vector.
   * @param  scalar_prod  [out] Array (length <tt>X.domain()->dim()</tt>) containing the
   *                      scalar products <tt>scalar_prod[j-1] = this->scalarProd(*X.col(j),*Y.col(j))</tt>,
   *                      for <tt>j = 1 ... X.domain()->dim()</tt>.
   *
   * Preconditions:<ul>
   * <li><tt>X.range()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li><tt>Y.range()->isCompatible(*this)</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * <li><tt>X.domain()->isCompatible(*Y.domain())</tt> (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>scalar_prod[j-1] = this->scalarProd(*X.col(j),*Y.col(j))</tt>, for <tt>j = 1 ... X.domain()->dim()</tt>
   * </ul>
   */
  virtual void scalarProds( const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y, Scalar scalar_prods[] ) const = 0;

  //@}

  /** @name Public virtual functions with default implementations */
  //@{

  /** \brief Returns if all of the vector elements are cheaply accessible
   * on this processor.
   *
   * The default implementation returns <tt>false</tt>.
   */
  virtual bool isInCore() const;

  /** \brief Clone this object (if supported).
   *
   * It is allowed for <tt>return.get()==NULL</tt> which means that this
   * capability is optional.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>.
   */
  virtual Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > clone() const;

  //@}

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED

  /** \name Public friends */
  //@{

  friend Teuchos::RefCountPtr< VectorBase<Scalar> >
  createMember<>( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs );

  friend Teuchos::RefCountPtr< MultiVectorBase<Scalar> >
  createMembers<>( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, int numMembers );

  friend Teuchos::RefCountPtr<VectorBase<Scalar> >
  createMemberView<>( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::MutableSubVectorT<Scalar> &raw_v );

  friend Teuchos::RefCountPtr<const VectorBase<Scalar> >
  createMemberView<>( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::SubVectorT<Scalar> &raw_v );

  friend Teuchos::RefCountPtr<MultiVectorBase<Scalar> >
  createMembersView<>( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv );

  friend Teuchos::RefCountPtr<const MultiVectorBase<Scalar> >
  createMembersView<>( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::SubMultiVectorT<Scalar> &raw_mv );

  //@}

#endif

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED
protected:
#endif

  /** @name Protected pure virtual functions that must be overridden */
  //@{

  /** \brief Create a vector member from the vector space.
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * <li> <tt>return->dim() == this->dim()</tt>
   * <li> <tt>return->space()->isCompatible(*this) == true</tt>
   * </ul>
   *
   * Warning, the returned <tt>VectorBase</tt> object is only guaranteed to be
   * valid while <tt>*this</tt> <tt>VectorSpaceBase</tt> object is still in
   * existence.
   *
   * @return Returns a smart reference counted pointer to a
   * dynamically allocated vector object.  After construction the
   * values in the vector <tt>*return</tt> are unspecified
   * (uninitialized).  This allows for faster execution times.  Note
   * that <tt>return->space().get() == this</tt> need not be true.
   */
  virtual Teuchos::RefCountPtr< VectorBase<Scalar> > createMember() const = 0;

  /** \brief Create a set of vector members (a <tt>MultiVectorBase</tt>) from the vector space.
   *
   * Preconditions:<ul>
   * <li> <tt>num_vecs >= 1</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return->range()->isCompatible(*this) == true</tt>
   * <li> <tt>return->domain()->dim() == numMembers</tt>
   * </ul>
   *
   * Warning, the returned <tt>MultiVectorBase</tt> object is only guaranteed
   * to be valid while <tt>*this</tt> <tt>VectorSpaceBase</tt> object is still
   * in existence.
   *
   * @return This method returns a smart reference counted pointer
   * to a dynamically allocated multi-vector object.  After
   * construction, the values in <tt>*return</tt> are unspecified
   * (uninitialized).  This allows for faster execution times.  Note
   * that <tt>return->range().get()==this</tt> does not have to be true
   * but will be in may cases.
   */
  virtual Teuchos::RefCountPtr< MultiVectorBase<Scalar> > createMembers(int numMembers) const = 0;

  /** \brief Create a vector member that is a non-<tt>const</tt> view of raw data.
   *
   * @param  raw_v  [in] On input contains pointer (i.e. <tt>raw_v.values()</tt>)
   *                to array that the returned <tt>VectorBase</tt> will be a view of.
   *                The data pointed to by <tt>raw_v.values()</tt> must remain
   *                valid until the returned <tt>VectorBase</tt> object is destroyed.
   *
   * Preconditions:<ul>
   * <li><tt>raw_v</tt> has been initialized to memory (i.e.
   *     <tt>raw_v.subDim()!=0 && raw_v.values()!=NULL</tt>).
   * <li><tt>raw_v</tt> is *consistent* with the local storage
   *     of a vector's data.  This precondition is purposefully vague since
   *     this function can be used an variety of specialized use-cases.
   * </ul>
   *
   * Postconditions:<ul>
   * <li>See <tt>createMember()</tt>
   * </ul>
   *
   * It is stated here that the client can not expect that the values
   * pointed to by <tt>raw_v.values()</tt> to be changed until
   * the smart pointer returned goes out of scope.  This is to allow a
   * default implementation that temporarily copies data into and out
   * of a <tt>VectorBase</tt> object using explicit vector access.
   */
  virtual Teuchos::RefCountPtr<VectorBase<Scalar> > createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const = 0;

  /** \brief Create a vector member that is a <tt>const</tt> view of raw data.
   *
   * @param  raw_v  [in] On input contains pointer (i.e. <tt>raw_v.values()</tt>)
   *                to array that the returned <tt>VectorBase</tt> will be a view of.
   *                The data pointed to by <tt>raw_v.values()</tt> must remain
   *                valid until the returned <tt>VectorBase</tt> object is destroyed.
   *
   * This function works exactly the same as the version that takes
   * a <tt>RTOpPack::MutableSubVectorT</tt> object except that this
   * version takes a <tt>RTOpPack::SubVectorT</tt> and returns a
   * smart pointer to a <tt>const</tt> <tt>VectorBase</tt> object.
   *
   * Preconditions:<ul>
   * <li>See the <tt>RTOpPack::MutableSubVectorT</tt> version of this function.
   * </ul>
   *
   * Postconditions:<ul>
   * <li>See <tt>createMember()</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<const VectorBase<Scalar> > createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const = 0;

  /** \brief Create a multi-vector member that is a non-<tt>const</tt> view of raw data.
   *
   * @param  raw_mv  [in] On input contains pointer (i.e. <tt>raw_mv.values()</tt>)
   *                 to array that the returned <tt>MultiVectorBase</tt> will be a view of.
   *
   * Preconditions:<ul>
   * <li><tt>raw_mv</tt> has been initialized to memory (i.e.
   *     <tt>raw_mv.subDim()!=0 && raw_mv.values()!=NULL</tt>).
   * <li><tt>raw_mv</tt> is *consistent* with the local storage
   *     of a vector's data.  This precondition is purposefully vague since
   *     this function can be used an variety of specialized use-cases.
   * </ul>
   *
   * Postconditions:<ul>
   * <li>See <tt>createMembers()</tt> where <tt>numMembers==raw_mv.numSubCols()</tt>
   * </ul>
   *
   * It is stated here that the client can not expect that the values
   * pointed to by <tt>raw_mv.values()</tt> to be changed until
   * the smart pointer returned goes out of scope.  This is to allow a
   * default implementation that temporarily copies data into and out
   * of a <tt>MultiVectorBase</tt> object using explicit vector access.
   */
  virtual Teuchos::RefCountPtr<MultiVectorBase<Scalar> > createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const = 0;

  /** \brief Create a multi-vector member that is a <tt>const</tt> view of raw data.
   *
   * @param  raw_mv  [in] On input contains pointer (i.e. <tt>raw_mv.values()</tt>)
   *                 to array that the returned <tt>MultiVectorBase</tt> will be a view of.
   *                 The data pointed to by <tt>raw_mv.values()</tt> must remain
   *                 valid until the returned <tt>MultiVectorBase</tt> object is destroyed.
   *
   * This function works exactly the same as the version that takes
   * a <tt>RTOpPack::MutableSubMultiVectorT</tt> object except that this
   * version takes a <tt>RTOpPack::SubMultiVectorT</tt> and returns a
   * smart pointer to a <tt>const</tt> <tt>MultiVectorBase</tt> object.
   *
   * Preconditions:<ul>
   * <li>See the <tt>RTOpPack::MutableSubMultiVectorT</tt> version of this function.
   * </ul>
   *
   * Postconditions:<ul>
   * <li>See <tt>createMember()</tt>
   * </ul>
   */
  virtual Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const = 0;

  //@}

}; // end class VectorSpaceBase

} // end namespace Thyra

// Inline members

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
Thyra::makeHaveOwnership( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs )
{
  if( vs.has_ownership() ) return vs;
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > _vs = vs->clone();
  TEST_FOR_EXCEPTION(
    _vs.get() == NULL, std::logic_error
    ,"Thyra::makeHaveOwnership(vs): Error, the concrete VectorSpaceBase object identified as \'"
    << vs->description() << "\' does not support the clone() function!"
    );
  return _vs;
}

template<class Scalar>
inline Teuchos::RefCountPtr< Thyra::VectorBase<Scalar> >
Thyra::createMember( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs )
{
  Teuchos::RefCountPtr<VectorBase<Scalar> > v = vs->createMember();
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &v );
  return v;
}
  
template<class Scalar>
inline Teuchos::RefCountPtr< Thyra::VectorBase<Scalar> >
Thyra::createMember( const VectorSpaceBase<Scalar> &vs )
{
  return createMember(Teuchos::rcp(&vs,false));
}

template<class Scalar>
inline Teuchos::RefCountPtr< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, int numMembers )
{
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > mv = vs->createMembers(numMembers);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &mv );
  return mv;
}

template<class Scalar>
inline Teuchos::RefCountPtr< Thyra::MultiVectorBase<Scalar> >
Thyra::createMembers( const VectorSpaceBase<Scalar> &vs, int numMembers )
{
  return createMembers(Teuchos::rcp(&vs,false),numMembers);
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
Thyra::createMemberView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::MutableSubVectorT<Scalar> &raw_v )
{
  Teuchos::RefCountPtr<VectorBase<Scalar> > v = vs->createMemberView(raw_v);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &v );
  return v;
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
Thyra::createMemberView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::MutableSubVectorT<Scalar> &raw_v )
{
  return createMemberView(Teuchos::rcp(&vs,false),raw_v);
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
Thyra::createMemberView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::SubVectorT<Scalar> &raw_v )
{
  Teuchos::RefCountPtr<const VectorBase<Scalar> > v = vs->createMemberView(raw_v);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &v );
  return v;
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::VectorBase<Scalar> >
Thyra::createMemberView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::SubVectorT<Scalar> &raw_v )
{
  return createMemberView(Teuchos::rcp(&vs,false),raw_v);
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(  const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv )
{
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > mv = vs->createMembersView(raw_mv);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &mv );
  return mv;
}

template<class Scalar>
inline Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView(  const VectorSpaceBase<Scalar> &vs, const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv )
{
  return createMembersView(Teuchos::rcp(&vs,false),raw_mv);
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView( const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > &vs, const RTOpPack::SubMultiVectorT<Scalar> &raw_mv )
{
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > mv = vs->createMembersView(raw_mv);
  Teuchos::set_extra_data( makeHaveOwnership(vs), "VectorSpaceBase", &mv );
  return mv;
}

template<class Scalar>
inline Teuchos::RefCountPtr<const Thyra::MultiVectorBase<Scalar> >
Thyra::createMembersView( const VectorSpaceBase<Scalar> &vs, const RTOpPack::SubMultiVectorT<Scalar> &raw_mv )
{
  return createMembersView(Teuchos::rcp(&vs,false),raw_mv);
}

#endif  // THYRA_VECTOR_SPACE_BASE_DECL_HPP

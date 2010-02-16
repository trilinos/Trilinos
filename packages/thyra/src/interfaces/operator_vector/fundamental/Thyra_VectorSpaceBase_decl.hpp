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


/** \brief Helper function that clones a <tt>VectorSpaceBase</tt> object if
 * the <tt>RCP</tt> does not have ownership.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<const VectorSpaceBase<Scalar> >
makeHaveOwnership( const RCP<const VectorSpaceBase<Scalar> > &vs );


/** \brief Create a vector member from the vector space.
 *
 * Calls <tt>VectorSpaceBase::createMember()</tt> on <tt>vs</tt> but
 * the returned <tt>VectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP< VectorBase<Scalar> >
createMember(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const std::string &label=""
  );


/** \brief Calls <tt>createMember(Teuchos::rcp(&vs,false))</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP< VectorBase<Scalar> >
createMember( const VectorSpaceBase<Scalar> &vs, const std::string &label="" );


/** \brief Create a set of vector members (a <tt>MultiVectorBase</tt>) from the vector space.
 *
 * Calls <tt>VectorSpaceBase::createMembers()</tt> on <tt>vs</tt> but the
 * returned <tt>MultiVectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP< MultiVectorBase<Scalar> >
createMembers(
  const RCP<const VectorSpaceBase<Scalar> > &vs, int numMembers,
  const std::string &label=""
  );


/** \brief Create a set of vector members (a <tt>MultiVectorBase</tt>) from the vector space.
 *
 * Calls <tt>VectorSpaceBase::createMembers()</tt> on <tt>vs</tt> but the
 * returned <tt>MultiVectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP< MultiVectorBase<Scalar> >
createMembers(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RCP<const VectorSpaceBase<Scalar> > &domain,
  const std::string &label=""
  );


/** \brief Calls <tt>createMembers(Teuchos::rcp(&vs,false),numMembers)</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP< MultiVectorBase<Scalar> >
createMembers(
  const VectorSpaceBase<Scalar> &vs, int numMembers,
  const std::string &label=""
  );


/** \brief Create a vector member that is a non-<tt>const</tt> view of raw data.
 *
 * Calls <tt>VectorSpaceBase::createMemberView()</tt> on <tt>vs</tt> but the
 * returned <tt>VectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<VectorBase<Scalar> >
createMemberView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::SubVectorView<Scalar> &raw_v,
  const std::string &label=""
  );


/** \brief Calls <tt>createMemberView(Teuchos::rcp(&vs,false),raw_v)</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<VectorBase<Scalar> >
createMemberView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::SubVectorView<Scalar> &raw_v,
  const std::string &label=""
  );


/** \brief Create a vector member that is a <tt>const</tt> view of raw data.
 *
 * Calls <tt>VectorSpaceBase::createMemberView()</tt> on <tt>vs</tt> but the
 * returned <tt>VectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<const VectorBase<Scalar> >
createMemberView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
  const std::string &label=""
  );


/** \brief Calls <tt>createMemberView(Teuchos::rcp(&vs,false),raw_v)</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<const VectorBase<Scalar> >
createMemberView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
  const std::string &label=""
  );


/** \brief Create a multi-vector member that is a non-<tt>const</tt> view of raw data.
 *
 * Calls <tt>VectorSpaceBase::createMembersView()</tt> on <tt>vs</tt> but the
 * returned <tt>MultiVectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<MultiVectorBase<Scalar> >
createMembersView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
  const std::string &label=""
  );


/** \brief Calls <tt>createMembersView(Teuchos::rcp(&vs,false),raw_mv)</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<MultiVectorBase<Scalar> >
createMembersView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
  const std::string &label=""
  );


/** \brief Create a multi-vector member that is a <tt>const</tt> view of raw data.
 *
 * Calls <tt>VectorSpaceBase::createMembersView()</tt> on <tt>vs</tt> but the
 * returned <tt>MultiVectorBase</tt> object can live past <tt>vs</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
createMembersView(
  const RCP<const VectorSpaceBase<Scalar> > &vs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
  const std::string &label=""
  );


/** \brief Calls <tt>createMembersView(Teuchos::rcp(&vs,false),raw_mv)</tt>.
 *
 * \relates VectorSpaceBase
 */
template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
createMembersView(
  const VectorSpaceBase<Scalar> &vs,
  const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
  const std::string &label=""
  );


/** \brief Abstract interface for objects that represent a space for vectors.
 *
 * This interface acts primarily as an "Abstract Factory" interface for
 * creating <tt>VectorBase</tt> objects using the non-member function
 * <tt>Thyra::createMember()</tt>.  A <tt>%VectorSpaceBase</tt> can also
 * create <tt> MultiVectorBase</tt> objects which represent a compact
 * collection of vectors using the non-member function
 * <tt>Thyra::createMembers()</tt>.  A secondary role for
 * <tt>%VectorSpaceBase</tt> objects is to test for compatibility of vector
 * space objects using the <tt>isCompatible()</tt> method and to apply the
 * space's scalar (inner) product.
 *
 * Clients can not directly create <tt>%VectorBase</tt> and
 * <tt>%MultiVectorBase</tt> objects using the member functions
 * <tt>createMember()</tt>, <tt>createMembers()</tt>,
 * <tt>createMemberView()</tt>, and <tt>createMembersView()</tt> but instead
 * must use the non-member \ref Thyra_Op_Vec_createMember_grp.
 *
 * Note that a <tt>%VectorSpaceBase</tt> object must also be able to create
 * <tt>MultiVectorBase</tt> objects with any number of column vectors using
 * the <tt>Thyra::createMembers()</tt> function.  Once created, the
 * <tt>LinearOpBase::domain()</tt> function supported by a created
 * <tt>%MultiVectorBase</tt> object returns a vector space of dimension equal
 * to the number of columns in the multi-vector.  An interesting implication
 * of this design is that the creation of a multi-vector provides a way for
 * clients to create vector spaces of any arbitrary (although small usually)
 * dimension.  In order to give the client the same ability to create smaller
 * vector spaces without having to create a dummy multi-vector object first,
 * the method <tt>smallVecSpcFcty()</tt> is included.  The method
 * <tt>smallVecSpcFcty()</tt> returns a <tt>VectorSpaceFactoryBase</tt> object
 * that can create (typically serial) vector spaces of any small dimension
 * that are compatible with the domain spaces of <tt>%MultiVectorBase</tt>
 * objects created by the vector space object.
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
 * <b>Notes for subclass developers</b>
 *
 * This is a fairly bare-bones interface class without much in the way of
 * default function implementations.  The subclass
 * <tt>VectorSpaceDefaultBase</tt> provides a default multi-vector
 * implementation and should the the first choice for subclass
 * implementations.
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
  virtual Ordinal dim() const = 0;

  /** \brief Compare the compatibility of two vector spaces.
   *
   * If this function returns <tt>true</tt>, then vectors created
   * from either of the vector spaces will be compatible and can be
   * combined in vector operations.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->dim() > 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li> [<tt>this->dim() != vecSpc.dim()</tt>] <tt>returnVal == false</tt>
   *
   * </ul>
   *
   * <b>Invariants:</b><ul>
   *
   * <li> [<tt>this->isCompatible(vecSpc) == true</tt>]
   * <tt>vecSpc.isCompatible(*this) == true</tt>
   *
   * </ul>
   */
  virtual bool isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const = 0;

  /** \brief Return a <tt>VectorSpaceFactoryBase</tt> object for the creation
   * of (usually serial) vector spaces with a small dimension.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->dim() > 0</tt>
   * </ul>
   */
  virtual RCP< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const = 0;

  /** \brief Return the scalar product of two vectors in the vector space.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>x.space()->isCompatible(*this)</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li><tt>y.space()->isCompatible(*this)</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * </ul>
   */
  virtual Scalar scalarProd(
    const VectorBase<Scalar>& x, const VectorBase<Scalar>& y
    ) const = 0;

  /** \brief Return the scalar product of each column in two multi-vectors in
   * the vector space.
   *
   * @param X [in] Multi-vector.
   *
   * @param Y [in] Multi-vector.
   *
   * @param scalarProds_out [out] Array (length <tt>X.domain()->dim()</tt>)
   * containing the scalar products <tt>scalarProds_out[j] =
   * this->scalarProds_out(*X.col(j),*Y.col(j))</tt>, for <tt>j = 0
   * ... X.domain()->dim()-1</tt>.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>X.range()->isCompatible(*this)</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li><tt>Y.range()->isCompatible(*this)</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * <li><tt>X.domain()->isCompatible(*Y.domain())</tt> (throw
   * <tt>Exceptions::IncompatibleVectorSpaces</tt>)
   *
   * </ul>
   */
  void scalarProds(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out
    ) const
    {
      scalarProdsImpl(X, Y, scalarProds_out);
    }

  //@}

  /** @name Public virtual functions with default implementations */
  //@{

  /** \brief Return if this vector space has a Euclidean (identity) basis in
   * which case the scalar product is the same as the dot product.
   *
   * The default implementation returns <tt>false</tt> (even though on average
   * the Euclidean scalar product is used).
   */
  virtual bool isEuclidean() const;

  /** \brief Returns <tt>true</tt> if <tt>this->acquireDetachedView(rng,...)</tt> returns
   * a direct view of the range of data requested.
   *
   * \param rng [in] The range of elements for the view (see
   * <tt>acquireDetachedView()</tt>).  The default value is <tt>Range1D()</tt>
   * (i.e. All of the elements in the vector).
   *
   * \param viewType [in] The type of view allowed.
   *
   * \param strideType [in] The type of stride the view is allowed to be.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>this->dim() > 0</tt>
   *
   * <li><tt>full_range(rng,0,this->dim()-1).ubound < this->dim()</tt>
   *
   * </ul>
   *
   * There are three different questions about the behavior of the <tt>acquireDetachedView()</tt>
   * that this query function can answer:
   *
   * <ul>
   *
   * <li>The elements in <tt>rng</tt> are fairly cheaply accessble in local
   * (i.e. in-core) memory if <tt>this->hasInCoreView(rng)==true</tt>.  Note
   * that this also allows for detached temporary copies of data.
   *
   * <li>A direct view of the elements in <tt>rng</tt> is available in local
   * (i.e. in-core) memory if
   * <tt>this->hasInCoreView(rng,VIEW_TYPE_DIRECT)==true</tt>.  No copy of
   * data is allowed here.
   *
   * <li>A direct view of the elements in <tt>rng</tt> with unit stride is
   * available in local (i.e. in-core) memory if
   * <tt>this->hasInCoreView(rng,VIEW_TYPE_DIRECT,STRIDE_TYPE_UNIT)==true</tt>
   * No copy of data is allowed here.
   *
   * </ul>
   *
   * The default implementation returns <tt>false</tt> (i.e. by default we do
   * not assume that any direct and/or contiguous views of any range of data
   * are provided).
   */
  virtual bool hasInCoreView(
    const Range1D &rng = Range1D(),
    const EViewType viewType = VIEW_TYPE_DETACHED,
    const EStrideType strideType = STRIDE_TYPE_NONUNIT
    ) const;

  /** \brief Clone this object (if supported).
   *
   * It is allowed for <tt>returnVal.get()==NULL</tt> which means that this
   * capability is optional.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->dim() > 0</tt>
   * </ul>
   *
   * The default implementation returns <tt>returnVal.get()==NULL</tt>.
   */
  virtual RCP< const VectorSpaceBase<Scalar> > clone() const;

  //@}

#ifndef DOXYGEN_COMPILE

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED

  /** \name Public friend functions */
  //@{

  friend RCP< VectorBase<Scalar> >
  createMember<>( 
    const RCP<const VectorSpaceBase<Scalar> > &vs,
    const std::string &label
    );

  friend RCP< MultiVectorBase<Scalar> >
  createMembers<>(
    const RCP<const VectorSpaceBase<Scalar> > &vs,
    int numMembers, const std::string &label
    );

  friend RCP<VectorBase<Scalar> >
  createMemberView<>(
    const RCP<const VectorSpaceBase<Scalar> > &vs,
    const RTOpPack::SubVectorView<Scalar> &raw_v,
    const std::string &label
    );

  friend RCP<const VectorBase<Scalar> >
  createMemberView<>(
    const RCP<const VectorSpaceBase<Scalar> > &vs,
    const RTOpPack::ConstSubVectorView<Scalar> &raw_v,
    const std::string &label
    );

  friend RCP<MultiVectorBase<Scalar> >
  createMembersView<>(
    const RCP<const VectorSpaceBase<Scalar> > &vs,
    const RTOpPack::SubMultiVectorView<Scalar> &raw_mv,
    const std::string &label
    );

  friend RCP<const MultiVectorBase<Scalar> >
  createMembersView<>(
    const RCP<const VectorSpaceBase<Scalar> > &vs,
    const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv,
    const std::string &label
    );

  //@}

#endif // DOXYGEN_COMPILE

#endif // TEMPLATE_FRIENDS_NOT_SUPPORTED

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED
protected:
#endif

  /** @name Protected pure virtual functions that must be overridden */
  //@{

  /** \brief Create a vector member from the vector space.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>this->dim() > 0</tt>
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li> <tt>returnVal.get() != NULL</tt>
   *
   * <li> <tt>returnVal->space()->isCompatible(*this) == true</tt>
   *
   * </ul>
   *
   * <b>Note:</b> This function is not to be called directly since it is
   * protected!  See the \ref Thyra_Op_Vec_createMember_grp.
   *
   * \returns A smart reference counted pointer to a dynamically allocated
   * vector object.  After construction the values in the vector
   * <tt>*returnVal</tt> are unspecified (uninitialized).  This allows for
   * faster execution times.  Note that <tt>returnVal->space().get() ==
   * this</tt> need not be true.
   */
  virtual RCP< VectorBase<Scalar> > createMember() const = 0;

  /** \brief Create a set of vector members (a <tt>MultiVectorBase</tt>) from
   * the vector space.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>this->dim() > 0</tt>
   *
   * <li> <tt>num_vecs >= 1</tt>
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li> <tt>returnVal->range()->isCompatible(*this) == true</tt>
   *
   * <li> <tt>returnVal->domain()->dim() == numMembers</tt>
   *
   * </ul>
   *
   * \returns A smart reference-counted pointer to a dynamically allocated
   * multi-vector object.  After construction, the values in
   * <tt>*returnVal</tt> are unspecified (uninitialized).  This allows for
   * faster execution times.  Note that
   * <tt>returnVal->range().get()==this</tt> does not have to be true but will
   * be in may cases.
   */
  virtual RCP< MultiVectorBase<Scalar> >
  createMembers(int numMembers) const = 0;

  /** \brief Create a vector member that is a non-<tt>const</tt> view of raw vector data.
   *
   * @param raw_v [in] On input contains pointer
   * (i.e. <tt>raw_v.values()</tt>) to array that the returned
   * <tt>VectorBase</tt> will be a view of.  The data pointed to by
   * <tt>raw_v.values()</tt> must remain valid until the returned
   * <tt>VectorBase</tt> object is destroyed.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>raw_v</tt> has been initialized to memory (i.e.
   * <tt>raw_v.subDim()!=0 && raw_v.values()!=NULL</tt>).
   *
   * <li><tt>raw_v</tt> is <em>consistent</em> with the local storage of this
   * vector spaces vector data.  This precondition is purposefully vague since
   * this function can be used an variety of specialized use-cases.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>See <tt>this->createMember()</tt>
   *
   * </ul>
   *
   * It is stated here that the client can not expect that the values pointed
   * to by <tt>raw_v.values()</tt> to be changed until the smart pointer
   * <tt>returnVal</tt> goes out of scope.  This is to allow an implementation
   * that temporarily copies data into and out of a <tt>VectorBase</tt> object
   * using explicit vector access.
   */
  virtual RCP<VectorBase<Scalar> >
  createMemberView( const RTOpPack::SubVectorView<Scalar> &raw_v ) const = 0;

  /** \brief Create a vector member that is a <tt>const</tt> view of raw vector data.
   *
   * @param raw_v [in] On input contains pointer
   * (i.e. <tt>raw_v.values()</tt>) to array that the returned
   * <tt>VectorBase</tt> will be a view of.  The data pointed to by
   * <tt>raw_v.values()</tt> must remain valid until the returned
   * <tt>VectorBase</tt> object is destroyed.
   *
   * This function works exactly the same as the previous version that takes a
   * <tt>RTOpPack::SubVectorView</tt> object except that this version
   * takes a <tt>RTOpPack::ConstSubVectorView</tt> and returns a smart pointer to a
   * <tt>const</tt> <tt>VectorBase</tt> object.
   *
   * <b>Preconditions:</b><ul>
   * <li>See the previous <tt>RTOpPack::SubVectorView</tt> version of this function.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>See <tt>this->createMember()</tt>
   * </ul>
   */
  virtual RCP<const VectorBase<Scalar> >
  createMemberView( const RTOpPack::ConstSubVectorView<Scalar> &raw_v ) const = 0;

  /** \brief Create a multi-vector member that is a non-<tt>const</tt> view of
   * raw multi-vector data.
   *
   * @param raw_mv [in] On input contains pointer
   * (i.e. <tt>raw_mv.values()</tt>) to array that the returned
   * <tt>MultiVectorBase</tt> will be a view of.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li><tt>raw_mv</tt> has been initialized to memory (i.e.
   * <tt>raw_mv.subDim()!=0 && raw_mv.values()!=NULL</tt>).
   *
   * <li><tt>raw_mv</tt> is <em>consistent</em> with the local storage of this
   * spaces vector data.  This precondition is purposefully vague since this
   * function can be used an variety of specialized use-cases.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>See <tt>this->createMembers()</tt> where
   * <tt>numMembers==raw_mv.numSubCols()</tt>
   *
   * </ul>
   *
   * It is stated here that the client can not expect that the values pointed
   * to by <tt>raw_mv.values()</tt> to be changed until the smart pointer
   * <tt>returnVal</tt> goes out of scope.  This is to allow for an
   * implementation that temporarily copies data into and out of a
   * <tt>MultiVectorBase</tt> object using explicit vector access.
   */
  virtual RCP<MultiVectorBase<Scalar> >
  createMembersView( const RTOpPack::SubMultiVectorView<Scalar> &raw_mv ) const = 0;

  /** \brief Create a multi-vector member that is a <tt>const</tt> view of raw
   * multi-vector data.
   *
   * @param raw_mv [in] On input contains pointer
   * (i.e. <tt>raw_mv.values()</tt>) to array that the returned
   * <tt>MultiVectorBase</tt> will be a view of.  The data pointed to by
   * <tt>raw_mv.values()</tt> must remain valid until the returned
   * <tt>MultiVectorBase</tt> object is destroyed.
   *
   * This function works exactly the same as the previous version that takes a
   * <tt>RTOpPack::SubMultiVectorView</tt> object except that this version
   * takes a <tt>RTOpPack::ConstSubMultiVectorView</tt> object and returns a smart
   * pointer to a <tt>const</tt> <tt>MultiVectorBase</tt> object.
   *
   * <b>Preconditions:</b><ul>
   *
   * <li>See the previous <tt>RTOpPack::SubMultiVectorView</tt> version of
   * this function.
   *
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   *
   * <li>See <tt>this->createMember()</tt>
   *
   * </ul>
   */
  virtual RCP<const MultiVectorBase<Scalar> >
  createMembersView(
    const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv ) const = 0;

  //@}

protected:

  /** \name Protected virtual funtions. */

  /** \brief . */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds
    ) const = 0;

public:

  /** \name Deprecated . */
  //@{

  /** \brief Deprecated . */
  void scalarProds(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    Scalar scalarProds[]
    ) const;

  //@}

private:
  
  // Not defined and not to be called
  VectorSpaceBase<Scalar>&
  operator=(const VectorSpaceBase<Scalar>&);

};


} // end namespace Thyra


#endif  // THYRA_VECTOR_SPACE_BASE_DECL_HPP

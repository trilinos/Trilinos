// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// @HEADER

#ifndef THYRA_EPETRA_THYRA_WRAPPERS_HPP
#define THYRA_EPETRA_THYRA_WRAPPERS_HPP


#include "Thyra_EpetraTypes.hpp"


namespace Teuchos { template<class Ordinal> class Comm; }


namespace Thyra {


/** \brief Given an <tt>Epetra_Comm</tt> object, return an equivalent
 * <tt>Teuchos::Comm</tt> object.
 *
 * If a successful conversion could not be performed then
 * <tt>return.get()==NULL</tt>.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Teuchos::Comm<Ordinal> >
create_Comm( const RCP<const Epetra_Comm> &epetraComm );


/** \brief Create an <tt>VectorSpaceBase</tt> object given an
 * <tt>Epetra_Map</tt> object.
 *
 * \param epetra_map [in] The Epetra map defining the partitioning of elements
 * to processors.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>epetra_map.get() != NULL</tt>
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li> <tt>return.get() != NULL</tt>
 * <li> The <tt>RCP</tt> object <tt>epetra_map</tt> is copied into
 *      the <tt>return</tt> object and therefore a memory of epetra_map is
 *      kept.
 * </ul>
 *
 * This uses an <tt>Epetra_Map</tt> object to initialize a compatible
 * <tt>VectorSpaceBase</tt> object.
 *
 * The fact that this function only accepts an <tt>Epetra_Map</tt> object
 * means that only maps that have elements of size one can be used to define a
 * vector space.  General <tt>Epetra_BlockMap</tt>s can not be used.  This is
 * not a serious limitation since <tt>Epetra_Operator</tt>'s domain and range
 * maps are of type <tt>Epetra_Map</tt>.
 *
 * This function works properly even if Epetra is not compiled with support
 * for SPMD (i.e. <tt>HAVE_MPI</tt> is not defined when compiling and
 * linking).  If SPMD support is not compiled into Epetra, then a serial
 * implementation of the communication is used instead.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const VectorSpaceBase<double> >
create_VectorSpace(
  const RCP<const Epetra_Map> &epetra_map
  );


/** \brief Create a <tt>VectorSpaceBase</tt> object that creates locally
 * replicated vector objects.
 *
 * \param parentSpace [in] The vector space that will be used to create the
 * smaller locally-replicated vector space.
 *
 * \param dim [in] The dimension of the locally replicated vector space.
 *
 * Note: This vector space will be compatible with the domain space of a
 * multivector. which has the range space <tt>parentSpace</tt>.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const VectorSpaceBase<double> >
create_LocallyReplicatedVectorSpace(
  const RCP<const VectorSpaceBase<double> > &parentSpace,
  const int dim
  );


/** \brief Create a non-<tt>const</tt> <tt>VectorBase</tt> object from a
 * non-<tt>const</tt> <tt>Epetra_Vector</tt> object.
 *
 * \param  epetra_v  [in] Smart pointer to the <tt>Epetra_Vector</tt> object to wrap.
 * \param  space     [in] The vector space that is compatible with <tt>epetra_v->Map()</tt>.
 *
 * <b>Precondiitions:</b><ul>
 * <li>[<tt>epetra_v.get()!=NULL</tt>] <tt>space.get()!=NULL</tt>
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li>[<tt>epetra_v.get()==NULL</tt>] <tt>return.get()==NULL<tt>
 * <li>[<tt>epetra_v.get()!=NULL</tt>] <tt>return.get()!=NULL<tt>
 * </ul>
 *
 * \return The returned <tt>RCP</tt> object contains a copy of the
 * input <tt>RCP<Epetra_Vector></tt> wrapped <tt>Epetra_Vector</tt>
 * object.  It is also stated that <tt>*epetra_v</tt> will only be guaranteed
 * to be modifed after the last <tt>RCP</tt> to the returned
 * <tt>VectorBase</tt> is destroyed.  In addition, <tt>*return</tt> is only
 * valid as long as one <tt>RefCoutPtr</tt> wrapper object still exits.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<VectorBase<double> >
create_Vector(
  const RCP<Epetra_Vector> &epetra_v,
  const RCP<const VectorSpaceBase<double> > &space
  );


/** \brief Create an <tt>const</tt> <tt>VectorBase</tt> wrapper object for
 * a <tt>const</tt> <tt>Epetra_Vector</tt> object.
 *
 * \param  epetra_v  [in] Smart pointer to the <tt>Epetra_Vector</tt> object to wrap.
 * \param  space     [in] The vector space that is compatible with <tt>epetra_v->Map()</tt>.
 *
 * <b>Precondiitions:</b><ul>
 * <li>[<tt>epetra_v.get()!=NULL</tt>] <tt>space.get()!=NULL</tt>
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li>[<tt>epetra_v.get()==NULL</tt>] <tt>return.get()==NULL<tt>
 * <li>[<tt>epetra_v.get()!=NULL</tt>] <tt>return.get()!=NULL<tt>
 * </ul>
 *
 * \return The returned <tt>RCP</tt> object contains a copy of the
 * input <tt>RCP<Epetra_Vector></tt> wrapped
 * <tt>Epetra_Vector</tt> object.  In addition, <tt>*return</tt> is only
 * valid as long as one <tt>RefCoutPtr</tt> wrapper object still exits.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const VectorBase<double> >
create_Vector(
  const RCP<const Epetra_Vector> &epetra_v,
  const RCP<const VectorSpaceBase<double> > &space
  );


/** \brief Create a non-<tt>const</tt> <tt>MultiVectorBase</tt> object from a
 * non-<tt>const</tt> <tt>Epetra_MultiVector</tt> object.
 *
 * \param  epetra_mv  [in] Smart pointer to the <tt>Epetra_MultiVector</tt> object to wrap.
 * \param  range      [in] The vector space that is compatible with <tt>epetra_mv->Map()</tt>.
 * \param  domain     [in] The vector space that is compatible with <tt>epetra_mv.NumVectors</tt>.
 *                    If <tt>domain.get()==NULL</tt>, then a space will be created internally.
 *
 * <b>Precondiitions:</b><ul>
 * <li><tt>epetra_mv.get()!=NULL</tt>
 * <li><tt>range.get()!=NULL</tt>
 * </ul>
 *
 * \return The returned <tt>RCP</tt> object contains a copy of the
 * input <tt>RCP<Epetra_MultiVector></tt> wrapped
 * <tt>Epetra_MultiVector</tt> object.  It is also stated that
 * <tt>*epetra_mv</tt> will only be guaranteed to be modifed after the last
 * <tt>RCP</tt> to the returned <tt>MultiVectorBase</tt> is destroyed.
 * In addition, <tt>*return</tt> is only valid as long as one
 * <tt>RefCoutPtr</tt> wrapper object still exits.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<MultiVectorBase<double> >
create_MultiVector(
  const RCP<Epetra_MultiVector> &epetra_mv,
  const RCP<const VectorSpaceBase<double> > &range,
  const RCP<const VectorSpaceBase<double> > &domain = Teuchos::null
  );


/** \brief Create an <tt>const</tt> <tt>MultiVectorBase</tt> wrapper object
 * for a <tt>const</tt> <tt>Epetra_MultiVector</tt> object.
 *
 * \param epetra_mv [in] Smart pointer to the <tt>Epetra_MultiVector</tt>
 * object to wrap.
 *
 * \param range [in] The vector space that is compatible with
 * <tt>epetra_mv->Map()</tt>.
 *
 * \param domain [in] The vector space that is compatible with
 * <tt>epetra_mv.NumVectors</tt>.  If <tt>domain.get()==NULL</tt>, then a
 * space will be created internally.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>epetra_mv.get()!=NULL</tt>
 * <li><tt>range.get()!=NULL</tt>
 * </ul>
 *
 * \return The returned <tt>RCP</tt> object contains a copy of the
 * input <tt>RCP<Epetra_MultiVector></tt> wrapped
 * <tt>Epetra_MultiVector</tt> object.  In addition, <tt>*return</tt> is only
 * valid as long as one <tt>RefCoutPtr</tt> wrapper object still exits.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const MultiVectorBase<double> >
create_MultiVector(
  const RCP<const Epetra_MultiVector> &epetra_mv,
  const RCP<const VectorSpaceBase<double> > &range,
  const RCP<const VectorSpaceBase<double> > &domain = Teuchos::null
  );


/** \brief Get (or create) and <tt>Epetra_Comm</tt> given a
 * <tt>Teuchos::Comm</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Epetra_Comm>
get_Epetra_Comm(const Teuchos::Comm<Ordinal>& comm);


/** \brief Get (or create) an <tt>Epetra_Map</tt> object given an
 * <tt>VectorSpaceBase</tt> object an optionally an extra <tt>Epetra_Comm</tt>
 * object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Epetra_Map>
get_Epetra_Map(const VectorSpaceBase<double>& vs, const RCP<const Epetra_Comm>& comm);


/** \brief Get a non-<tt>const</tt> <tt>Epetra_Vector</tt> view from a
 * non-<tt>const</tt> <tt>VectorBase</tt> object if possible.
 *
 * Preconditions:<ul>
 * <li> <tt>v.get()!=NULL</tt>
 * <li> <tt>map</tt> must be compatible with <tt>*v.space()</tt>
 * </ul>
 *
 * If a <tt>RCP<Epetra_Vector></tt> object is already
 * attached to the node of the smart pointer for <tt>mv</tt> then this is
 * returned directly.  If not, then a view of the data in <tt>*v</tt> is
 * created and returned.  In the latter case the smart pointer <tt>v</tt> is
 * copied and attached to the returned RCP object.  Therefore, a temporary
 * <tt>VectorBase</tt> object can be created in the call to this function and
 * the view in <tt>return</tt> will persist until all of the RCP objects to
 * the returned <tt>Epetra_Vector</tt> object go away.
 *
 * Note: the <tt>v</tt> object is not guaranteed to be modified until the last
 * smart pointer to the returned <tt>Epetra_Vector</tt> object is destroyed.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<Epetra_Vector>
get_Epetra_Vector(
  const Epetra_Map &map,
  const RCP<VectorBase<double> > &v
  );


/** \brief Get a <tt>const</tt> <tt>Epetra_Vector</tt> view from a
 * <tt>const</tt> <tt>VectorBase</tt> object if possible.
 *
 * Preconditions:<ul>
 * <li> <tt>v.get()!=NULL</tt>
 * <li> <tt>map</tt> must be compatible with <tt>*v->space()</tt>
 * </ul>
 *
 * If a <tt>RCP<Epetra_Vector></tt> object is already
 * attached to the node of the smart pointer for <tt>mv</tt> then this is
 * returned directly.  If not, then a view of the data in <tt>*v</tt> is
 * created and returned.  In the latter case the smart pointer <tt>v</tt> is
 * copied and attached to the returned RCP object.  Therefore, a temporary
 * <tt>VectorBase</tt> object can be created in the call to this function and
 * the view in <tt>return</tt> will persist until all of the RCP objects to
 * the returned <tt>Epetra_Vector</tt> object go away.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Epetra_Vector>
get_Epetra_Vector(
  const Epetra_Map &map,
  const RCP<const VectorBase<double> > &v
  );


/** \brief Get a non-<tt>const</tt> <tt>Epetra_MultiVector</tt> view from a
 * non-<tt>const</tt> <tt>MultiVectorBase</tt> object if possible.
 *
 * <b>Preconditions:</b><ul>
 * <li> <tt>mv.get()!=NULL</tt>
 * <li> <tt>map</tt> must be compatible with <tt>*mv->range()</tt>
 * </ul>
 *
 * If a <tt>RCP<Epetra_MultiVector></tt> object is already
 * attached to the node of the smart pointer for <tt>mv</tt> then this is
 * returned directly.  If not, then a view of the data in <tt>*mv</tt> is
 * created and returned.  In the latter case the smart pointer <tt>mv</tt> is
 * copied and attached to the returned RCP object.  Therefore, a temporary
 * <tt>MultiVectorBase</tt> object can be created in the call to this function
 * and the view in <tt>return</tt> will persist until all of the RCP objects
 * to the returned <tt>Epetra_MultiVector</tt> object go away.
 *
 * Note: the <tt>mv</tt> object is not guaranteed to be modified until
 * the last smart pointer to the returned <tt>Epetra_MultiVector</tt>
 * object is destroyed.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<Epetra_MultiVector>
get_Epetra_MultiVector(
  const Epetra_Map &map,
  const RCP<MultiVectorBase<double> > &mv
  );


/** \brief Get a <tt>const</tt> <tt>Epetra_MultiVector</tt> view from a
 * <tt>const</tt> <tt>MultiVectorBase</tt> object if possible.
 *
 * <b>Preconditions:</b><ul>
 * <li> <tt>mv.get()!=NULL</tt>
 * <li> <tt>map</tt> must be compatible with <tt>*mv.range()</tt>
 * </ul>
 *
 * If a <tt>RCP<const Epetra_MultiVector></tt> object is
 * already attached to the node of the smart pointer for <tt>mv</tt> then this
 * is returned directly.  If not, then a view of the data in <tt>*mv</tt> is
 * created and returned.  In the latter case the smart pointer <tt>mv</tt> is
 * copied and attached to the returned RCP object.  Therefore, a temporary
 * <tt>MultiVectorBase</tt> object can be created in the call to this function
 * and the view in <tt>return</tt> will persist until all of the RCP objects
 * to the returned <tt>Epetra_MultiVector</tt> object go away.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Epetra_MultiVector>
get_Epetra_MultiVector(
  const Epetra_Map &map, 
  const RCP<const MultiVectorBase<double> > &mv
  );


/** \brief Get a non-<tt>const</tt> <tt>Epetra_MultiVector</tt> view from a
 * non-<tt>const</tt> <tt>MultiVectorBase</tt> object if possible where the
 * client must maintain the memory of the input multivector.
 *
 * <b>Preconditions:</b><ul>
 * <li> <tt>map</tt> must be compatible with <tt>*mv.range()</tt>
 * </ul>
 *
 * This function trys to dynamic cast some some known interfaces classes where
 * data can be directly accessed and no RCP magic needs to be used.  This
 * results in improved performance in time-critical use cases (like when
 * called for <tt>EpetraLinearOp</tt> in the inner loop of a Krylov solver).
 *
 * If this function can not dynamic cast to the direct data access interfaces
 * it punts and calls the more general (but more expensive)
 * <tt>get_Epetra_MultiVector()</tt> function.
 *
 * Note: the <tt>mv</tt> object is not guaranteed to be modified until
 * the last smart pointer to the returned <tt>Epetra_MultiVector</tt>
 * object is destroyed.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
Teuchos::RCP<Epetra_MultiVector>
get_Epetra_MultiVector(
  const Epetra_Map &map,
  MultiVectorBase<double> &mv
  );


/** \brief Get a <tt>const</tt> <tt>Epetra_MultiVector</tt> view from a
 * <tt>const</tt> <tt>MultiVectorBase</tt> object if possible where the client
 * must maintain the memory of the input multivector.
 *
 * <b>Preconditions:</b><ul>
 * <li> <tt>map</tt> must be compatible with <tt>*mv.range()</tt>
 * </ul>
 *
 * This function trys to dynamic cast some some known interfaces classes where
 * data can be directly accessed and no RCP magic needs to be used.  This
 * results in improved performance in time-critical use cases (like when
 * called for <tt>EpetraLinearOp</tt> in the inner loop of a Krylov solver).
 *
 * If this function can not dynamic cast to the direct data access interfaces
 * it punts and calls the more general (but more expensive)
 * <tt>get_Epetra_MultiVector()</tt> function.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
Teuchos::RCP<const Epetra_MultiVector>
get_Epetra_MultiVector(
  const Epetra_Map &map,
  const MultiVectorBase<double> &mv
  );


} // namespace Thyra


#endif // THYRA_EPETRA_THYRA_WRAPPERS_HPP

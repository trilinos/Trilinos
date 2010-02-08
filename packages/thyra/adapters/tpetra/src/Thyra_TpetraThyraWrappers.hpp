
#ifndef THYRA_TPETRA_THYRA_WRAPPERS_HPP
#define THYRA_TPETRA_THYRA_WRAPPERS_HPP

#ifndef THYRA_TPETRA_CANNOT_COME_OUT_TO_PLAY_RIGHT_NOW

#include "Thyra_TpetraTypes.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#  include "Tpetra_MpiComm.hpp"
#  include "Teuchos_DefaultMpiComm.hpp"
#endif

namespace Thyra {

/** \defgroup Thyra_Tpetra_Thyra_Wrappers_grp  Collection of functions for wrapping and unwrapping Tpetra objects

\ingroup Tpetra_Thyra_Op_Vec_adapters_grp

This set of functions provides some general utility code for wrapping %Tpetra
objects in %Thyra object and for getting %Tpetra views of %Thyra objects.

*/

/** \brief Given a <tt>Tpetra::Comm</tt> object, return an equivalent
 * <tt>Teuchos::Comm>/tt> object.
 *
 * If a successful conversion could not be performed then
 * <tt>return.get()==NULL</tt>.
 */
template<typename Ordinal, typename Packet>
RCP<const Teuchos::Comm<Ordinal> >
create_Comm(
  const RCP<const Tpetra::Comm<Ordinal,Packet> > &tpetraComm
  );

/** \brief Concrete an <tt>SpmdVectorSpaceBase</tt> object given an
 * <tt>Tpetra::VectorSpace</tt> object.
 *
 * \param  tpetra_vs  [in] The Tpetra vector space defining the partitioning of elements
 *                     to processors.
 *
 * Preconditions:<ul>
 * <li><tt>tpetra_vs.get() != NULL</tt>
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li> <tt>return.get() != NULL</tt>
 * <li> The <tt>RCP</tt> object <tt>tpetra_vs</tt> is copied into
 *      the <tt>return</tt> object and therefore a memory of tpetra_vs is
 *      kept.
 * </ul>
 *
 * This uses an <tt>Tpetra::VectorSpace</tt> object to initialize a compatible
 * <tt>DefaultSpmdVectorSpace</tt> object.
 *
 * \ingroup Thyra_Tpetra_Thyra_Wrappers_grp
 */
template<class Ordinal, class Scalar>
RCP<const SpmdVectorSpaceDefaultBase<Scalar> >
create_VectorSpace(
  const RCP<const Tpetra::VectorSpace<Ordinal,Scalar> > &tpetra_vs
  );

/** \brief Create a non-<tt>const</tt> <tt>SpmdVectorBase</tt> object from
 * a <tt>const> <tt>Tpetra::Vector</tt> object.
 *
 * @param  tpetra_v  [in] Smart pointer to the <tt>Tpetra::Vector</tt> object to wrap.
 * @param  space     [in] The vector space that is compatible with <tt>tpetra_v->vectorSpace()</tt>.
 *                   If <tt>space.get()==NULL</tt> then this space will be created from
 *                   <tt>tpetra_v->vectorSpace()</tt>.
 *
 * <b>Postconditions:</b><ul>
 * <li>[<tt>tpetra_v.get()==NULL</tt>] <tt>return.get()==NULL<tt>
 * <li>[<tt>tpetra_v.get()!=NULL</tt>] <tt>return.get()!=NULL<tt>
 * </ul>
 *
 * \return The returned <tt>RCP</tt> object contains a copy of the
 * input <tt>RCP<Tpetra::Vector<Ordinal,Scalar> ></tt> wrapped
 * <tt>Tpetra::Vector</tt> object.  It is also stated that
 * <tt>*tpetra_v</tt> will only be guaranteed to be modifed after the last
 * <tt>RCP</tt> to the returned <tt>SpmdVectorBase</tt> object is
 * destroyed.  In addition, <tt>*return</tt> is only valid as long as
 * one <tt>RefCoutPtr</tt> wrapper object still exits.
 *
 * \ingroup Thyra_Tpetra_Thyra_Wrappers_grp
 */
template<class Ordinal, class Scalar>
RCP<SpmdVectorBase<Scalar> >
create_Vector(
  const RCP<Tpetra::Vector<Ordinal,Scalar> >      &tpetra_v
  ,const RCP<const SpmdVectorSpaceBase<Scalar> >   &space       = Teuchos::null
  );

/** \brief Create an <tt>const</tt> <tt>SpmdVectorBase</tt> wrapper object
 * for a <tt>const</tt> <tt>Tpetra::Vector</tt> object.
 *
 * @param  tpetra_v  [in] Smart pointer to the <tt>Tpetra::Vector</tt> object to wrap.
 * @param  space     [in] The vector space that is compatible with <tt>tpetra_v->vectorSpace()</tt>.
 *                   If <tt>space.get()==NULL</tt> then this space will be created from
 *                   <tt>tpetra_v->vectorSpace()</tt>.
 *
 * <b>Postconditions:</b><ul>
 * <li>[<tt>tpetra_v.get()==NULL</tt>] <tt>return.get()==NULL<tt>
 * <li>[<tt>tpetra_v.get()!=NULL</tt>] <tt>return.get()!=NULL<tt>
 * </ul>
 *
 * \return The returned <tt>RCP</tt> object contains a copy of the
 * input <tt>RCP<Tpetra::Vector<Ordinal,Scalar> ></tt> wrapped
 * <tt>Tpetra::Vector</tt> object.  In addition, <tt>*return</tt> is only
 * valid as long as one <tt>RefCoutPtr</tt> wrapper object still exits.
 *
 * \ingroup Thyra_Tpetra_Thyra_Wrappers_grp
 */
template<class Ordinal, class Scalar>
RCP<const SpmdVectorBase<Scalar> >
create_Vector(
  const RCP<const Tpetra::Vector<Ordinal,Scalar> >  &tpetra_v
  ,const RCP<const SpmdVectorSpaceBase<Scalar> >     &space       = Teuchos::null
  );

/** \brief Get a non-<tt>const</tt> <tt>Tpetra::Vector</tt> view from a
 * non-<tt>const</tt> <tt>VectorBase</tt> object if possible.
 *
 * <b>Preconditions:</b><ul>
 * <li> <tt>v.get()!=NULL</tt>
 * <li> <tt>tpetra_vs</tt> must be compatible with <tt>*v.space()</tt>
 * </ul>
 *
 * If a <tt>RCP<Tpetra::Vector<Ordinal,Scalar> ></tt> object is already
 * attached to the node of the smart pointer for <tt>v</tt> then this is
 * returned directly.  If not, then a view of the data in <tt>*v</tt> is
 * created and returned.  In the latter case the smart pointer <tt>v</tt> is
 * copied and attached to the returned RCP object.  Therefore, a temporary
 * <tt>VectorBase</tt> object can be created in the call to this function and
 * the view in <tt>*return</tt> will persist until all of the RCP objects to
 * the returned <tt>Tpetra::Vector</tt> object go away.
 *
 * Note: The <tt>v</tt> object is not guaranteed to be modified until the last
 * smart pointer to the returned <tt>Tpetra::Vector</tt> object is destroyed.
 *
 * \ingroup Thyra_Tpetra_Thyra_Wrappers_grp
 */
template<class Ordinal, class Scalar>
RCP<Tpetra::Vector<Ordinal,Scalar> >
get_Tpetra_Vector(
  const Tpetra::VectorSpace<Ordinal,Scalar>           &tpetra_vs
  ,const RCP<VectorBase<Scalar> >    &v
  );

/** \brief Get a <tt>const</tt> <tt>Tpetra::Vector</tt> view from a
 * <tt>const</tt> <tt>VectorBase</tt> object if possible.
 *
 * <b>Preconditions:</b><ul>
 * <li> <tt>v.get()!=NULL</tt>
 * <li> <tt>tpetra_vs</tt> must be compatible with <tt>*v.space()</tt>
 * </ul>
 *
 * If a <tt>RCP<Tpetra::Vector<Ordinal,Scalar> ></tt> object is already
 * attached to the node of the smart pointer for <tt>v</tt> then this is
 * returned directly.  If not, then a view of the data in <tt>*v</tt> is
 * created and returned.  In the latter case the smart pointer <tt>v</tt> is
 * copied and attached to the returned RCP object.  Therefore, a temporary
 * <tt>VectorBase</tt> object can be created in the call to this function and
 * the view in <tt>*return</tt> will persist until all of the RCP objects to
 * the returned <tt>Tpetra::Vector</tt> object go away.
 *
 * \ingroup Thyra_Tpetra_Thyra_Wrappers_grp
 */
template<class Ordinal, class Scalar>
RCP<const Tpetra::Vector<Ordinal,Scalar> >
get_Tpetra_Vector(
  const Tpetra::VectorSpace<Ordinal,Scalar>                &tpetra_vs
  ,const RCP<const VectorBase<Scalar> >   &v
  );
/** \brief Get smart pointer to non-<tt>const</tt>
 * <tt>Tpetra::Operator</tt> object from reference to a
 * non-<tt>const</tt> <tt>TpetraLinearOp</tt> accessed through its
 * <tt>LinearOpBase</tt> interface.
 *
 * @param op [in] Reference to operator to extract <tt>Tpetra::Operator</tt> out of.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>dynamic_cast<TpetraLinearOp<Ordinal,Scalar>*>(&op) != NULL</tt>
 * </ul>
 *
 * This function is designed to provide an easy way for non-C++ experts
 * to get at the <tt>Tpetra::Operator</tt> object that was stuck into
 * an <tt>TpetraLinearOp</tt> object.
 *
 * If the dynamic cast fails then a <tt>std::bad_cast</tt> exception
 * is thrown containing a very detailed error message as to why the
 * cast failed.
 *
 * This function is simple enough and developers can see what needs to
 * be done to accomplish this type of access by looking at the source
 * code by clicking on:
 *
 * \ingroup Tpetra_Thyra_Op_Vec_get_Tpetra_Operator_grp
 */
template<class Ordinal, class Scalar>
RCP<Tpetra::Operator<Ordinal,Scalar> >
get_Tpetra_Operator( LinearOpBase<Scalar> &op );

/** \brief Get smart pointer to <tt>const</tt>
 * <tt>Tpetra::Operator</tt> object from reference to a <tt>const</tt>
 * <tt>TpetraLinearOp</tt> accessed through its <tt>LinearOpBase</tt>
 * interface.
 *
 * @param op [in] Reference to operator to extract <tt>Tpetra::Operator</tt> out of.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>dynamic_cast<const TpetraLinearOp<Ordinal,Scalar>*>(&op) != NULL</tt>
 * </ul>
 *
 * This function is designed to provide an easy way for non-C++ experts
 * to get at the <tt>Tpetra::Operator</tt> object that was stuck into
 * an <tt>TpetraLinearOp</tt> object.
 *
 * If the dynamic cast fails then a <tt>std::bad_cast</tt> exception
 * is thrown containing a very detailed error message as to why the
 * cast failed.
 *
 * This function is simple enough and developers can see what needs to
 * be done to accomplish this type of access by looking at the source
 * code by clicking on:
 *
 * \ingroup Tpetra_Thyra_Op_Vec_get_Tpetra_Operator_grp
 */
template<class Ordinal, class Scalar>
RCP<const Tpetra::Operator<Ordinal,Scalar> >
get_Tpetra_Operator( const LinearOpBase<Scalar> &op );

} // namespace Thyra

// //////////////////////////////
// Implementations

#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_VectorSpace.hpp"

namespace Thyra {
namespace TpetraUtils {

// Utility class to copy data back from a Tpetra vector to a Thyra vector.
template<class Ordinal, class Scalar>
class CopyFromTpetraToThyraVector {
public:
  CopyFromTpetraToThyraVector(
    Tpetra::Vector<Ordinal,Scalar>                          const& tpetra_v
    ,RCP<DetachedVectorView<Scalar> > const& detachedView
    )
    :tpetra_v_(tpetra_v),detachedView_(detachedView)
    {}
  ~CopyFromTpetraToThyraVector()
    {
      for( Ordinal i = 0; i < detachedView_->subDim(); ++i )
        (*detachedView_)[i] = tpetra_v_[i];
    }
private:
  Tpetra::Vector<Ordinal,Scalar> const& tpetra_v_; // Can not be a RCP due to circular references!
  RCP<DetachedVectorView<Scalar> > detachedView_;
  // Not defined and not to be called
  CopyFromTpetraToThyraVector();
  CopyFromTpetraToThyraVector(const CopyFromTpetraToThyraVector&);
  CopyFromTpetraToThyraVector& operator=(const CopyFromTpetraToThyraVector&);
};

} // namespace TpetraUtils
} // namespace Thyra

template<typename Ordinal, typename Packet>
RCP<const Teuchos::Comm<Thyra::Ordinal> >
Thyra::create_Comm(
  const RCP<const Tpetra::Comm<Ordinal,Packet> > &tpetraComm
  )
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  RCP<const Tpetra::SerialComm<Ordinal,Packet> >
    serialTpetraComm = rcp_dynamic_cast<const Tpetra::SerialComm<Ordinal,Packet> >(tpetraComm);
  if( serialTpetraComm.get() ) {
    RCP<const Teuchos::SerialComm<Ordinal> >
      serialComm = rcp(new Teuchos::SerialComm<Ordinal>());
    set_extra_data( serialTpetraComm, "serialTpetraComm", Teuchos::outArg(serialComm) );
    return serialComm;
  }

#ifdef HAVE_MPI
  
  RCP<const Tpetra::MpiComm<Ordinal,Packet> >
    mpiTpetraComm = rcp_dynamic_cast<const Tpetra::MpiComm<Ordinal,Packet> >(tpetraComm);
  if( mpiTpetraComm.get() ) {
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiTpetraComm->getMpiComm());
    set_extra_data( mpiTpetraComm, "mpiTpetraComm", Teuchos::outArg(rawMpiComm) );
    RCP<const Teuchos::MpiComm<Ordinal> >
      mpiComm = rcp(new Teuchos::MpiComm<Ordinal>(rawMpiComm));
    return mpiComm;
  }

#endif // HAVE_MPI
  
  // If you get here then the failed!
  return Teuchos::null;

}

template<class Ordinal, class Scalar>
RCP<const Thyra::SpmdVectorSpaceDefaultBase<Scalar> >
Thyra::create_VectorSpace(
  const RCP<const Tpetra::VectorSpace<Ordinal,Scalar> > &tpetra_vs
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION( !tpetra_vs.get(), std::invalid_argument, "create_VectorSpace::initialize(...): Error!" );
#endif // TEUCHOS_DEBUG
  RCP<const Teuchos::Comm<Ordinal> >
    comm = create_Comm(Teuchos::rcp(&tpetra_vs->comm(),false)).assert_not_null();
  Teuchos::set_extra_data( tpetra_vs, "tpetra_vs", Teuchos::outArg(comm) );
  const Ordinal localSubDim = tpetra_vs->getNumMyEntries();
  RCP<DefaultSpmdVectorSpace<Scalar> >
    vs = Teuchos::rcp(
      new DefaultSpmdVectorSpace<Scalar>(
        comm
        ,localSubDim
        ,tpetra_vs->getNumGlobalEntries()
        )
      );
#ifndef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        vs->dim() != tpetra_vs->getNumGlobalEntries(), std::logic_error
        ,"create_VectorSpace(...): Error, vs->dim() = "<<vs->dim()<<" != "
        "tpetra_vs->getNumGlobalEntries() = "<<tpetra_vs->getNumGlobalEntries()<<"!"
        );
#endif		
  Teuchos::set_extra_data( tpetra_vs, "tpetra_vs", Teuchos::outArg(vs) );
  return vs;
}

template<class Ordinal, class Scalar>
RCP<Thyra::SpmdVectorBase<Scalar> >
Thyra::create_Vector(
  const RCP<Tpetra::Vector<Ordinal,Scalar> >             &tpetra_v
  ,const RCP<const Thyra::SpmdVectorSpaceBase<Scalar> >   &space_in
  )
{
  if(!tpetra_v.get()) return Teuchos::null;
  // Create the space if it is missing
  RCP<const Thyra::SpmdVectorSpaceBase<Scalar> >
    space = space_in;
  if(!space.get())
    space =
      create_VectorSpace<Ordinal,Scalar>(
        Teuchos::rcp(new Tpetra::VectorSpace<Ordinal,Scalar>(tpetra_v->vectorSpace()))
        );
  // New local view of raw data
  Scalar *localValues = &(*tpetra_v)[0]; // This points to contiguous memory!
  // Build the Vector with a view of the data
  RCP<SpmdVectorBase<Scalar> >
    v = Teuchos::rcp(
      new DefaultSpmdVector<Scalar>(
        space,
        Teuchos::arcp(localValues,0,tpetra_v->vectorSpace().getNumMyEntries(),false),
        1
        )
      );
  Teuchos::set_extra_data( tpetra_v, "Tpetra_Vector", Teuchos::outArg(v) );
  return v;
}

template<class Ordinal, class Scalar>
RCP<const Thyra::SpmdVectorBase<Scalar> >
Thyra::create_Vector(
  const RCP<const Tpetra::Vector<Ordinal,Scalar> >       &tpetra_v
  ,const RCP<const Thyra::SpmdVectorSpaceBase<Scalar> >   &space_in
  )
{
  if(!tpetra_v.get()) return Teuchos::null;
  // Create the space if it is missing
  RCP<const Thyra::SpmdVectorSpaceBase<Scalar> >
    space = space_in;
  if(!space.get())
    space =
      create_VectorSpace<Ordinal,Scalar>(
        Teuchos::rcp(new Tpetra::VectorSpace<Ordinal,Scalar>(tpetra_v->vectorSpace()))
        );
  // New local view of raw data
  const Scalar *localValues = &(*tpetra_v)[0]; // This points to contiguous memory!
  // Build the Vector with a view of the data
  RCP<const SpmdVectorBase<Scalar> >
    v = Teuchos::rcp(
      new DefaultSpmdVector<Scalar>(
        space,
        Teuchos::arcp(const_cast<Scalar*>(localValues),0,tpetra_v->vectorSpace().getNumMyEntries(),false),
        1
        )
      );
  Teuchos::set_extra_data( tpetra_v, "Tpetra_Vector", Teuchos::outArg(v) );
  return v;
}

template<class Ordinal, class Scalar>
RCP<Tpetra::Vector<Ordinal,Scalar> >
Thyra::get_Tpetra_Vector(
  const Tpetra::VectorSpace<Ordinal,Scalar>                  &tpetra_vs
  ,const RCP<Thyra::VectorBase<Scalar> >    &v
  )
{
  
  using Teuchos::get_optional_extra_data;
  
  //
  // Warning! There is some advanced stuff going on here using RCP and extra
  // data.  This code is not fragile but you have to know what you are doing
  // if you are going to modify it!
  //
#ifdef TEUCHOS_DEBUG
  RCP<const VectorSpaceBase<Scalar> >
    thyra_vs = create_VectorSpace<Ordinal,Scalar>(Teuchos::rcp(&tpetra_vs,false));
  THYRA_ASSERT_VEC_SPACES( "Thyra::get_Tpetra_Vector(tpetra_vs,v)", *thyra_vs, *v->space() );
#endif
  //
  // First, try to grab the Tpetra::Vector straight out of the
  // RCP since this is the fastest way.
  //
  const Ptr<const RCP<Tpetra::Vector<Ordinal,Scalar> > > tpetra_v_ptr =
    get_optional_extra_data<RCP<Tpetra::Vector<Ordinal,Scalar> > >(v,"Tpetra::Vector");
  if(!is_null(tpetra_v_ptr)) {
    return *tpetra_v_ptr;
  }
  //
  // The assumption that we (rightly) make here is that if the vector spaces
  // are compatible, that either the vectors are both in-core or the vector
  // spaces are both derived from SpmdVectorSpaceBase and have compatible
  // distributions.
  // 
  const VectorSpaceBase<Scalar>  &vs = *v->range();
  const SpmdVectorSpaceBase<Scalar> *mpi_vs = dynamic_cast<const SpmdVectorSpaceBase<Scalar>*>(&vs);
  const Ordinal localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Ordinal localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
  //
  // Here we will extract a view of the local elements in the underlying
  // VectorBase object.  In most cases, no data will be allocated or copied
  // and only some small objects will be created and a few virtual functions
  // will be called so the overhead should be low and fixed.  Note that 'v' is
  // "remembered" by 'detachedView' so that 'v' will not go away until the //
  // detached view is copied back.
  //
  RCP<DetachedVectorView<Scalar> >
    detachedView = Teuchos::rcp(
      new DetachedVectorView<Scalar>(
        v,Range1D(localOffset,localOffset+localSubDim-1)
        )
      );
  // Create a temporary Tpetra::Vector object and copy the local data into it.
  RCP<Tpetra::Vector<Ordinal,Scalar> >
    tpetra_v = Teuchos::rcp(new Tpetra::Vector<Ordinal,Scalar>(tpetra_vs));
  for( Ordinal i = 0; i < detachedView->subDim(); ++i )
    (*tpetra_v)[i] = (*detachedView)[i];
  // Create a utility object that will copy back the values in the
  // Tpetra::Vector into the Thyra vector when this temp Tpetra::Vector is
  // destroyed
  Teuchos::set_extra_data(
    Teuchos::rcp(new TpetraUtils::CopyFromTpetraToThyraVector<Ordinal,Scalar>(*tpetra_v,detachedView)),
    "CopyFromTpetraToThyraVector",
    Teuchos::outArg(tpetra_v),
    Teuchos::PRE_DESTROY
    );
  // Note that when the last RCP to tpetra_v is destroyed, just before the
  // Tpetra::Vector object itself is destroyed, the
  // CopyFromTpetraToThyraVector object will be destroyed and its destructor
  // will copy back the data to the detached view.  When the detached view is
  // destroyed (after the CopyFromTpetraToThyraVector is destoryed and the
  // data is copied back) the view will be commited back to the underlying
  // Thyra vector v.
  return tpetra_v;
  // Note that when Tpetra::Vector supports a view mode, we will initalize it
  // with the view in detachedView (with the stride forced to be unit) and we
  // will have to set detachedView as extra data to the RCP for tpetra_v!
}

template<class Ordinal, class Scalar>
RCP<const Tpetra::Vector<Ordinal,Scalar> >
Thyra::get_Tpetra_Vector(
  const Tpetra::VectorSpace<Ordinal,Scalar>                       &tpetra_vs
  ,const RCP<const Thyra::VectorBase<Scalar> >   &v
  )
{
  using Teuchos::get_optional_extra_data;
#ifdef TEUCHOS_DEBUG
  RCP<const VectorSpaceBase<Scalar> >
    thyra_vs = create_VectorSpace<Ordinal,Scalar>(Teuchos::rcp(&tpetra_vs,false));
  THYRA_ASSERT_VEC_SPACES( "Thyra::get_Tpetra_Vector(tpetra_vs,v)", *thyra_vs, *v->space() );
#endif
  //
  // First, try to grab the Tpetra::Vector straight out of the
  // RCP since this is the fastest way.
  //
  const Ptr<const RCP<const Tpetra::Vector<Ordinal,Scalar> > > tpetra_v_ptr =
    get_optional_extra_data<RCP<const Tpetra::Vector<Ordinal,Scalar> > >(v,"Tpetra::Vector");
  if(!is_null(tpetra_v_ptr))
    return *tpetra_v_ptr;
  const Ptr<const RCP<Tpetra::Vector<Ordinal,Scalar> > > tpetra_nonconst_v_ptr =
    get_optional_extra_data<RCP<Tpetra::Vector<Ordinal,Scalar> > >(v,"Tpetra::Vector");
  if(!is_null(tpetra_nonconst_v_ptr))
    return *tpetra_nonconst_v_ptr;
  //
  // The assumption that we (rightly) make here is that if the vector spaces
  // are compatible, that either the vectors are both in-core or the vector
  // spaces are both derived from SpmdVectorSpaceBase and have compatible
  // distributions.
  // 
  const VectorSpaceBase<Scalar>  &vs = *v->range();
  const SpmdVectorSpaceBase<Scalar> *mpi_vs = dynamic_cast<const SpmdVectorSpaceBase<Scalar>*>(&vs);
  const Ordinal localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Ordinal localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
  //
  // Here we will extract a view of the local elements in the underlying
  // VectorBase object.  In most cases, no data will be allocated or copied
  // and only some small objects will be created and a few virtual functions
  // will be called so the overhead should be low and fixed.  Note that 'v' is
  // "remembered" by 'detachedView' so that 'v' will not go away until the
  // detached view is finished being used and is destroyed.
  //
  RCP<ConstDetachedVectorView<Scalar> >
    detachedView = Teuchos::rcp(
      new ConstDetachedVectorView<Scalar>(
        v,Range1D(localOffset,localOffset+localSubDim-1)
        )
      );
  // Create a temporary Tpetra::Vector object and copy the local data into it.
  RCP<Tpetra::Vector<Ordinal,Scalar> >
    tpetra_v = Teuchos::rcp(new Tpetra::Vector<Ordinal,Scalar>(tpetra_vs));
  for( Ordinal i = 0; i < detachedView->subDim(); ++i )
    (*tpetra_v)[i] = (*detachedView)[i];
  return tpetra_v;
  // Note that when Tpetra::Vector supports a view mode, we will initalize it
  // with the view in detachedView (with the stride forced to be unit) and we
  // will have to set detachedView as extra data to the RCP for tpetra_v!
  // This is why I have detachedView managed as a RCP object even through it
  // goes away in this function.
}

template<class Ordinal, class Scalar>
RCP<Tpetra::Operator<Ordinal,Scalar> >
Thyra::get_Tpetra_Operator( Thyra::LinearOpBase<Scalar> &op )
{
  return Teuchos::dyn_cast<TpetraLinearOp<Ordinal,Scalar> >(op).getNonconstTpetraOp();
}

template<class Ordinal, class Scalar>
RCP<const Tpetra::Operator<Ordinal,Scalar> >
Thyra::get_Tpetra_Operator( const Thyra::LinearOpBase<Scalar> &op )
{
  return Teuchos::dyn_cast<const TpetraLinearOp<Ordinal,Scalar> >(op).getTpetraOp();
}

#endif

#endif // THYRA_TPETRA_THYRA_WRAPPERS_HPP

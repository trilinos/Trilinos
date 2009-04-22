// @HEADER
// ***********************************************************************
// 
//               Thyra: Trilinos Solver Framework Core
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

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"


Teuchos::RCP<const Teuchos::Comm<Thyra::Index> >
Thyra::create_Comm( const RCP<const Epetra_Comm> &epetraComm )
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(epetraComm);
  if( serialEpetraComm.get() ) {
    RCP<const Teuchos::SerialComm<Index> >
      serialComm = rcp(new Teuchos::SerialComm<Index>());
    set_extra_data( serialEpetraComm, "serialEpetraComm", Teuchos::inOutArg(serialComm) );
    return serialComm;
  }

#ifdef HAVE_MPI
  
  RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm>(epetraComm);
  if( mpiEpetraComm.get() ) {
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
    set_extra_data( mpiEpetraComm, "mpiEpetraComm", Teuchos::inOutArg(rawMpiComm) );
    RCP<const Teuchos::MpiComm<Index> >
      mpiComm = rcp(new Teuchos::MpiComm<Index>(rawMpiComm));
    return mpiComm;
  }

#endif // HAVE_MPI
  
  // If you get here then the failed!
  return Teuchos::null;

}


Teuchos::RCP<const Thyra::VectorSpaceBase<double> >
Thyra::create_VectorSpace(
  const RCP<const Epetra_Map> &epetra_map
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    !epetra_map.get(), std::invalid_argument,
    "create_VectorSpace::initialize(...): Error!" );
#endif // TEUCHOS_DEBUG
  RCP<const Teuchos::Comm<Index> >
    comm = create_Comm(Teuchos::rcp(&epetra_map->Comm(),false)).assert_not_null();
  Teuchos::set_extra_data( epetra_map, "epetra_map", Teuchos::inOutArg(comm) );
  const Index localSubDim = epetra_map->NumMyElements();
  RCP<DefaultSpmdVectorSpace<double> > vs =
    defaultSpmdVectorSpace<double>(
      comm, localSubDim, epetra_map->NumGlobalElements());
#ifndef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    vs->dim() != epetra_map->NumGlobalElements(), std::logic_error
    ,"create_VectorSpace(...): Error, vs->dim() = "<<vs->dim()<<" != "
    "epetra_map->NumGlobalElements() = "<<epetra_map->NumGlobalElements()<<"!"
    );
#endif		
  Teuchos::set_extra_data( epetra_map, "epetra_map", Teuchos::inOutArg(vs) );
  return vs;
}


Teuchos::RCP<const Thyra::VectorSpaceBase<double> >
Thyra::create_LocallyReplicatedVectorSpace(
  const RCP<const VectorSpaceBase<double> > &parentSpace,
  const int dim
  )
{
  using Teuchos::rcp_dynamic_cast;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(parentSpace.get()==NULL);
  Teuchos::dyn_cast<const SpmdVectorSpaceBase<double> >(*parentSpace);
  TEST_FOR_EXCEPT(dim <= 0);
#endif
  return parentSpace->smallVecSpcFcty()->createVecSpc(dim);
}


Teuchos::RCP<Thyra::VectorBase<double> >
Thyra::create_Vector(
  const RCP<Epetra_Vector> &epetra_v,
  const RCP<const VectorSpaceBase<double> > &space_in
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(space_in.get()==NULL);
#endif
  RCP<const SpmdVectorSpaceBase<double> >
    space = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
      space_in,true);
  if(!epetra_v.get())
    return Teuchos::null;
  // New local view of raw data
  double *localValues;
  epetra_v->ExtractView( &localValues );
  // Build the Vector
  RCP<SpmdVectorBase<double> >
    v = Teuchos::rcp(
      new DefaultSpmdVector<double>(
        space,
        Teuchos::arcp(localValues,0,epetra_v->Map().NumMyElements(),false),
        1
        )
      );
  Teuchos::set_extra_data<RCP<Epetra_Vector> >( epetra_v, "Epetra_Vector",
    Teuchos::inOutArg(v) );
  return v;
}


Teuchos::RCP<const Thyra::VectorBase<double> >
Thyra::create_Vector(
  const RCP<const Epetra_Vector> &epetra_v,
  const RCP<const VectorSpaceBase<double> > &space_in
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(space_in.get()==NULL);
#endif
  RCP<const SpmdVectorSpaceBase<double> >
    space = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
      space_in,true);
  if(!epetra_v.get())
    return Teuchos::null;
  // New local view of raw data
  double *localValues;
  epetra_v->ExtractView( &localValues );
  // Build the Vector
  RCP<const SpmdVectorBase<double> >
    v = Teuchos::rcp(
      new DefaultSpmdVector<double>(
        space,
        Teuchos::arcp(localValues,0,epetra_v->Map().NumMyElements(),false),
        1
        )
      );
  Teuchos::set_extra_data<RCP<const Epetra_Vector> >( epetra_v, "Epetra_Vector",
    Teuchos::inOutArg(v) );
  return v;
}


Teuchos::RCP<Thyra::MultiVectorBase<double> >
Thyra::create_MultiVector(
  const RCP<Epetra_MultiVector> &epetra_mv,
  const RCP<const VectorSpaceBase<double> > &range_in,
  const RCP<const VectorSpaceBase<double> > &domain_in
  )
{
  using Teuchos::rcp_dynamic_cast;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(range_in.get()==NULL);
#endif
  RCP<const SpmdVectorSpaceBase<double> >
    range = rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
      range_in,true);
  RCP<const ScalarProdVectorSpaceBase<double> >
    domain = rcp_dynamic_cast<const ScalarProdVectorSpaceBase<double> >(
      domain_in,true);
  if (!epetra_mv.get() )
    return Teuchos::null;
  if ( is_null(domain) ) {
    domain = rcp_dynamic_cast<const ScalarProdVectorSpaceBase<double> >(
      create_LocallyReplicatedVectorSpace(range,epetra_mv->NumVectors())
      );
  }
  // New local view of raw data
  double *localValues; int leadingDim;
  if( epetra_mv->ConstantStride() ) {
    epetra_mv->ExtractView( &localValues, &leadingDim );
  }
  else {
    TEST_FOR_EXCEPT(true); // ToDo: Implement views of non-contiguous mult-vectors!
  }
  // Build the MultiVector
  RCP<SpmdMultiVectorBase<double> >
    mv = Teuchos::rcp(
      new DefaultSpmdMultiVector<double>(
        range,
        domain,
        Teuchos::arcp(localValues,0,leadingDim*epetra_mv->NumVectors(),false),
        leadingDim
        )
      );
  Teuchos::set_extra_data<RCP<Epetra_MultiVector> >(
    epetra_mv, "Epetra_MultiVector", Teuchos::inOutArg(mv) );
  return mv;
}


Teuchos::RCP<const Thyra::MultiVectorBase<double> >
Thyra::create_MultiVector(
  const RCP<const Epetra_MultiVector> &epetra_mv,
  const RCP<const VectorSpaceBase<double> > &range_in,
  const RCP<const VectorSpaceBase<double> > &domain_in
  )
{
  using Teuchos::rcp_dynamic_cast;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(range_in.get()==NULL);
#endif
  RCP<const SpmdVectorSpaceBase<double> >
    range = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
      range_in,true);
  RCP<const ScalarProdVectorSpaceBase<double> >
    domain = Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<double> >(
      domain_in,true);
  if (!epetra_mv.get())
    return Teuchos::null;
  if ( is_null(domain) ) {
    domain = rcp_dynamic_cast<const ScalarProdVectorSpaceBase<double> >(
      create_LocallyReplicatedVectorSpace(range,epetra_mv->NumVectors())
      );
  }
  // New local view of raw data
  double *localValues; int leadingDim;
  if( epetra_mv->ConstantStride() ) {
    epetra_mv->ExtractView( &localValues, &leadingDim );
  }
  else {
    TEST_FOR_EXCEPT(true); // ToDo: Implement views of non-contiguous mult-vectors!
  }
  // Build the MultiVector
  RCP<const SpmdMultiVectorBase<double> >
    mv = Teuchos::rcp(
      new DefaultSpmdMultiVector<double>(
        range,
        domain,
        Teuchos::arcp(localValues,0,leadingDim*epetra_mv->NumVectors(),false),
        leadingDim
        )
      );
  Teuchos::set_extra_data<RCP<const Epetra_MultiVector> >(
    epetra_mv, "Epetra_MultiVector", Teuchos::inOutArg(mv) );
  return mv;
}


Teuchos::RCP<Epetra_Vector>
Thyra::get_Epetra_Vector(
  const Epetra_Map &map,
  const RCP<VectorBase<double> > &v
  )
{
  using Teuchos::get_optional_extra_data;
#ifdef TEUCHOS_DEBUG
  RCP<const VectorSpaceBase<double> >
    epetra_vs = create_VectorSpace(Teuchos::rcp(&map,false));
  THYRA_ASSERT_VEC_SPACES( 
    "Thyra::get_Epetra_Vector(map,v)", *epetra_vs, *v->space() );
#endif
  //
  // First, try to grab the Epetra_Vector straight out of the
  // RCP since this is the fastest way.
  //
  const Ptr<const RCP<Epetra_Vector> >
    epetra_v_ptr = get_optional_extra_data<RCP<Epetra_Vector> >(
      v,"Epetra_Vector");
  if(!is_null(epetra_v_ptr)) {
    return *epetra_v_ptr;
  }
  //
  // The assumption that we (rightly) make here is that if the vector spaces
  // are compatible, that either the multi-vectors are both in-core or the
  // vector spaces are both derived from SpmdVectorSpaceBase and have
  // compatible maps.
  // 
  const VectorSpaceBase<double>  &vs = *v->range();
  const SpmdVectorSpaceBase<double> *mpi_vs
    = dynamic_cast<const SpmdVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
  //
  // Here we will extract a view of the local elements in the underlying
  // VectorBase object.  In most cases, no data will be allocated or copied
  // and only some small objects will be created and a few virtual functions
  // will be called so the overhead should be low and fixed.
  //
  // Create a *mutable* view of the local elements, this view will be set on
  // the RCP that is returned.  As a result, this view will be relased
  // when the returned Epetra_Vector is released.
  //
  // Note that the input vector 'v' will be remembered through this detached
  // view!
  //
  RCP<DetachedVectorView<double> >
    emvv = Teuchos::rcp(
      new DetachedVectorView<double>(
        v
        ,Range1D(localOffset,localOffset+localSubDim-1)
        ,true // forceContiguous
        )
      );
  // Create a temporary Epetra_Vector object and give it
  // the above local view
  RCP<Epetra_Vector>
    epetra_v = Teuchos::rcp(
      new Epetra_Vector(
        ::View                                 // CV
        ,map                                   // Map
        ,const_cast<double*>(emvv->values())   // V
        )
      );
  // Give the explict view object to the above Epetra_Vector smart pointer
  // object.  In this way, when the client is finished with the Epetra_Vector
  // view the destructor from the object in emvv will automatically commit the
  // changes to the elements in the input v VectorBase object (reguardless of
  // its implementation).  This is truly an elegant result!
  Teuchos::set_extra_data( emvv, "emvv", Teuchos::inOutArg(epetra_v),
    Teuchos::PRE_DESTROY );
  // We are done!
  return epetra_v;
}


Teuchos::RCP<const Epetra_Vector>
Thyra::get_Epetra_Vector(
  const Epetra_Map &map, 
  const RCP<const VectorBase<double> > &v
  )
{
  using Teuchos::get_optional_extra_data;
#ifdef TEUCHOS_DEBUG
  RCP<const VectorSpaceBase<double> >
    epetra_vs = create_VectorSpace(Teuchos::rcp(&map,false));
  THYRA_ASSERT_VEC_SPACES(
    "Thyra::get_Epetra_Vector(map,v)", *epetra_vs, *v->space() );
#endif
  //
  // First, try to grab the Epetra_Vector straight out of the
  // RCP since this is the fastest way.
  //
  const Ptr<const RCP<const Epetra_Vector> >
    epetra_v_ptr = get_optional_extra_data<RCP<const Epetra_Vector> >(
      v,"Epetra_Vector");
  if(!is_null(epetra_v_ptr))
    return *epetra_v_ptr;
  const Ptr<const RCP<Epetra_Vector> >
    epetra_nonconst_v_ptr = get_optional_extra_data<RCP<Epetra_Vector> >(
      v,"Epetra_Vector");
  if(!is_null(epetra_nonconst_v_ptr))
    return *epetra_nonconst_v_ptr;
  //
  // Same as above function except as stated below
  //
  const VectorSpaceBase<double> &vs = *v->range();
  const SpmdVectorSpaceBase<double> *mpi_vs
    = dynamic_cast<const SpmdVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
  // Get an explicit *non-mutable* view of all of the elements in the multi
  // vector.  Note that 'v' will be remembered by this view!
  RCP<ConstDetachedVectorView<double> >
    evv = Teuchos::rcp(
      new ConstDetachedVectorView<double>(
        v
        ,Range1D(localOffset,localOffset+localSubDim-1)
        ,true // forceContiguous
        )
      );
  // Create a temporary Epetra_Vector object and give it
  // the above view
  RCP<Epetra_Vector>
    epetra_v = Teuchos::rcp(
      new Epetra_Vector(
        ::View                                  // CV
        ,map                                    // Map
        ,const_cast<double*>(evv->values())     // V
        )
      );
  // This next statement will cause the destructor to free the view if
  // needed (see above function).  Since this view is non-mutable,
  // only a releaseDetachedView(...) and not a commit will be called.
  // This is the whole reason there is a seperate implementation for
  // the const and non-const cases.
  Teuchos::set_extra_data( evv, "evv", Teuchos::inOutArg(epetra_v),
    Teuchos::PRE_DESTROY );
  // We are done!
  return epetra_v;
}


Teuchos::RCP<Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
  const Epetra_Map &map,
  const RCP<MultiVectorBase<double> > &mv
  )
{
  using Teuchos::get_optional_extra_data;
#ifdef TEUCHOS_DEBUG
  RCP<const VectorSpaceBase<double> >
    epetra_vs = create_VectorSpace(Teuchos::rcp(&map,false));
  THYRA_ASSERT_VEC_SPACES(
    "Thyra::get_Epetra_MultiVector(map,mv)", *epetra_vs, *mv->range() );
#endif
  //
  // First, try to grab the Epetra_MultiVector straight out of the
  // RCP since this is the fastest way.
  //
  const Ptr<const RCP<Epetra_MultiVector> >
    epetra_mv_ptr = get_optional_extra_data<RCP<Epetra_MultiVector> >(
      mv,"Epetra_MultiVector");
  if(!is_null(epetra_mv_ptr)) {
    return *epetra_mv_ptr;
  }
  //
  // The assumption that we (rightly) make here is that if the vector spaces
  // are compatible, that either the multi-vectors are both in-core or the
  // vector spaces are both derived from SpmdVectorSpaceBase and have
  // compatible maps.
  // 
  const VectorSpaceBase<double> &vs = *mv->range();
  const SpmdVectorSpaceBase<double> *mpi_vs
    = dynamic_cast<const SpmdVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
  //
  // Here we will extract a view of the local elements in the underlying
  // MultiVectorBase object.  In most cases, no data will be allocated or
  // copied and only some small objects will be created and a few virtual
  // functions will be called so the overhead should be low and fixed.
  //
  // Create a *mutable* view of the local elements, this view will be set on
  // the RCP that is returned.  As a result, this view will be relased
  // when the returned Epetra_MultiVector is released.
  //
  RCP<DetachedMultiVectorView<double> >
    emmvv = Teuchos::rcp(
      new DetachedMultiVectorView<double>(
        *mv
        ,Range1D(localOffset,localOffset+localSubDim-1)
        )
      );
  // Create a temporary Epetra_MultiVector object and give it
  // the above local view
  RCP<Epetra_MultiVector>
    epetra_mv = Teuchos::rcp(
      new Epetra_MultiVector(
        ::View                                  // CV
        ,map                                    // Map
        ,const_cast<double*>(emmvv->values())   // A
        ,emmvv->leadingDim()                    // MyLDA
        ,emmvv->numSubCols()                    // NumVectors
        )
      );
  // Give the explict view object to the above Epetra_MultiVector
  // smart pointer object.  In this way, when the client is finished
  // with the Epetra_MultiVector view the destructor from the object
  // in emmvv will automatically commit the changes to the elements in
  // the input mv MultiVectorBase object (reguardless of its
  // implementation).  This is truly an elegant result!
  Teuchos::set_extra_data( emmvv, "emmvv", Teuchos::inOutArg(epetra_mv),
    Teuchos::PRE_DESTROY );
  // Also set the mv itself as extra data just to be safe
  Teuchos::set_extra_data( mv, "mv", Teuchos::inOutArg(epetra_mv) );
  // We are done!
  return epetra_mv;
}


Teuchos::RCP<const Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
  const Epetra_Map &map,
  const RCP<const MultiVectorBase<double> > &mv
  )
{
  using Teuchos::get_optional_extra_data;
#ifdef TEUCHOS_DEBUG
  RCP<const VectorSpaceBase<double> >
    epetra_vs = create_VectorSpace(Teuchos::rcp(&map,false));
  THYRA_ASSERT_VEC_SPACES(
    "Thyra::get_Epetra_MultiVector(map,mv)",
    *epetra_vs, *mv->range() );
#endif
  //
  // First, try to grab the Epetra_MultiVector straight out of the
  // RCP since this is the fastest way.
  //
  const Ptr<const RCP<const Epetra_MultiVector> >
    epetra_mv_ptr
    = get_optional_extra_data<RCP<const Epetra_MultiVector> >(
      mv,"Epetra_MultiVector" );
  if(!is_null(epetra_mv_ptr)) {
    return *epetra_mv_ptr;
  }
  //
  // Same as above function except as stated below
  //
  const VectorSpaceBase<double> &vs = *mv->range();
  const SpmdVectorSpaceBase<double> *mpi_vs
    = dynamic_cast<const SpmdVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
  // Get an explicit *non-mutable* view of all of the elements in
  // the multi vector.
  RCP<ConstDetachedMultiVectorView<double> >
    emvv = Teuchos::rcp(
      new ConstDetachedMultiVectorView<double>(
        *mv
        ,Range1D(localOffset,localOffset+localSubDim-1)
        )
      );
  // Create a temporary Epetra_MultiVector object and give it
  // the above view
  RCP<Epetra_MultiVector>
    epetra_mv = Teuchos::rcp(
      new Epetra_MultiVector(
        ::View                                  // CV
        ,map                                    // Map
        ,const_cast<double*>(emvv->values())    // A
        ,emvv->leadingDim()                     // MyLDA
        ,emvv->numSubCols()                     // NumVectors
        )
      );
  // This next statement will cause the destructor to free the view if
  // needed (see above function).  Since this view is non-mutable,
  // only a releaseDetachedView(...) and not a commit will be called.
  // This is the whole reason there is a seperate implementation for
  // the const and non-const cases.
  Teuchos::set_extra_data( emvv, "emvv", Teuchos::inOutArg(epetra_mv),
    Teuchos::PRE_DESTROY );
  // Also set the mv itself as extra data just to be safe
  Teuchos::set_extra_data( mv, "mv", Teuchos::inOutArg(epetra_mv) );
  // We are done!
  return epetra_mv;
}


Teuchos::RCP<Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
  const Epetra_Map &map,
  MultiVectorBase<double> &mv
  )
{
  using Teuchos::rcp;
  double *mvData = 0;
  Index mvLeadingDim = -1;
  SpmdMultiVectorBase<double> *mvSpmdMv = 0;
  SpmdVectorBase<double> *mvSpmdV = 0;
  if ((mvSpmdMv = dynamic_cast<SpmdMultiVectorBase<double>*>(&mv))) {
    mvSpmdMv->getLocalData(&mvData,&mvLeadingDim);
  }
  else if ((mvSpmdV = dynamic_cast<SpmdVectorBase<double>*>(&mv))) {
    Index mvStride = -1;
    mvSpmdV->getLocalData(&mvData,&mvStride);
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_EQUALITY(mvStride,1);
#endif
    mvLeadingDim = mvSpmdV->spmdSpace()->localSubDim();
  }
  if (mvData) {
    rcp(
      new Epetra_MultiVector(
        ::View, map, mvData, mvLeadingDim, mv.domain()->dim()
        )
      );
  }
  return ::Thyra::get_Epetra_MultiVector(map,rcp(&mv,false));
}


Teuchos::RCP<const Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
  const Epetra_Map &map,
  const MultiVectorBase<double> &mv
  )
{
  using Teuchos::rcp;
  const double *mvData = 0;
  Index mvLeadingDim = -1;
  const SpmdMultiVectorBase<double> *mvSpmdMv = 0;
  const SpmdVectorBase<double> *mvSpmdV = 0;
  if ((mvSpmdMv = dynamic_cast<const SpmdMultiVectorBase<double>*>(&mv))) {
    mvSpmdMv->getLocalData(&mvData,&mvLeadingDim);
  }
  else if ((mvSpmdV = dynamic_cast<const SpmdVectorBase<double>*>(&mv))) {
    Index mvStride = -1;
    mvSpmdV->getLocalData(&mvData,&mvStride);
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_EQUALITY(mvStride,1);
#endif
    mvLeadingDim = mvSpmdV->spmdSpace()->localSubDim();
  }
  if (mvData) {
    rcp(
      new Epetra_MultiVector(
        ::View,map, const_cast<double*>(mvData), mvLeadingDim, mv.domain()->dim()
        )
      );
  }
  return ::Thyra::get_Epetra_MultiVector(map,rcp(&mv,false));
}

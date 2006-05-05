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
#include "Thyra_DefaultMPIVectorSpace.hpp"
#include "Thyra_DefaultMPIMultiVector.hpp"
#include "Thyra_DefaultMPIVector.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Teuchos_dyn_cast.hpp"

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#ifdef RTOp_USE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

Teuchos::RefCountPtr<const Thyra::MPIVectorSpaceDefaultBase<double> >
Thyra::create_MPIVectorSpaceBase(
	const Teuchos::RefCountPtr<const Epetra_Map> &epetra_map
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !epetra_map.get(), std::invalid_argument, "create_MPIVectorSpaceBase::initialize(...): Error!" );
#endif // _DEBUG
	MPI_Comm mpiComm;
#ifdef RTOp_USE_MPI
	const Epetra_MpiComm
		*epetra_mpi_comm = dynamic_cast<const Epetra_MpiComm*>(&epetra_map->Comm());
	if(epetra_mpi_comm) {
    //std::cout << "EpetraVectorSpace::initialize(...): Using an Epetra_MpiComm!\n";
		mpiComm = epetra_mpi_comm->Comm();
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			mpiComm == MPI_COMM_NULL, std::logic_error
			,"EpetraVectorSpace::initialize(...), Error, if using Epetra_MpiComm then "
			"the associated MPI_Comm object can not be MPI_COMM_NULL!"
			);
#endif // _DEBUG
     //std::cout << "EpetraVectorSpace::initialize(...): mpiComm = " << mpiComm << std::endl;
	}
	else {
		//std::cout << "EpetraVectorSpace::initialize(...): Not using an Epetra_MpiComm!\n";
		mpiComm = MPI_COMM_NULL;
	}
#else // RTOp_USE_MPI
  //std::cout << "EpetraVectorSpace::initialize(...): Not using an Epetra_MpiComm!\n";
	mpiComm = MPI_COMM_NULL;
#endif // RTOp_USE_MPI
	const Index localSubDim = epetra_map->NumMyElements();
	Teuchos::RefCountPtr<DefaultMPIVectorSpace<double> >
		vs = Teuchos::rcp(
			new DefaultMPIVectorSpace<double>(
				mpiComm
				,localSubDim
				,epetra_map->NumGlobalElements()
				)
			);
#ifndef _DEBUG
			TEST_FOR_EXCEPTION(
				vs->dim() != epetra_map->NumGlobalElements(), std::logic_error
				,"create_MPIVectorSpaceBase(...): Error, vs->dim() = "<<vs->dim()<<" != "
        "epetra_map->NumGlobalElements() = "<<epetra_map->NumGlobalElements()<<"!"
        );
#endif		
	Teuchos::set_extra_data( epetra_map, "epetra_map", &vs );
	return vs;
}

Teuchos::RefCountPtr<Thyra::MPIVectorBase<double> >
Thyra::create_MPIVectorBase(
	const Teuchos::RefCountPtr<Epetra_Vector>                                &epetra_v
	,const Teuchos::RefCountPtr<const MPIVectorSpaceBase<double> >           &space
	)
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(space.get()==NULL);
#endif
  if(!epetra_v.get()) return Teuchos::null;
  // New local view of raw data
  double *localValues;
  epetra_v->ExtractView( &localValues );
  // Build the Vector
	Teuchos::RefCountPtr<MPIVectorBase<double> >
    v = Teuchos::rcp(new DefaultMPIVector<double>(space,Teuchos::rcp(localValues,false),1));
  Teuchos::set_extra_data( epetra_v, "Epetra_Vector", &v );
  return v;
}

Teuchos::RefCountPtr<const Thyra::MPIVectorBase<double> >
Thyra::create_MPIVectorBase(
	const Teuchos::RefCountPtr<const Epetra_Vector>                          &epetra_v
	,const Teuchos::RefCountPtr<const MPIVectorSpaceBase<double> >           &space
	)
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(space.get()==NULL);
#endif
  if(!epetra_v.get()) return Teuchos::null;
  // New local view of raw data
  double *localValues;
  epetra_v->ExtractView( &localValues );
  // Build the Vector
	Teuchos::RefCountPtr<const MPIVectorBase<double> >
    v = Teuchos::rcp(new DefaultMPIVector<double>(space,Teuchos::rcp(localValues,false),1));
  Teuchos::set_extra_data( epetra_v, "Epetra_Vector", &v );
  return v;
}

Teuchos::RefCountPtr<Thyra::MPIMultiVectorBase<double> >
Thyra::create_MPIMultiVectorBase(
	const Teuchos::RefCountPtr<Epetra_MultiVector>                           &epetra_mv
	,const Teuchos::RefCountPtr<const MPIVectorSpaceBase<double> >           &range
	,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<double> >    &domain
	)
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
#endif
  if(!epetra_mv.get()) return Teuchos::null;
  // New local view of raw data
  double *localValues; int leadingDim;
	if( epetra_mv->ConstantStride() ) {
		epetra_mv->ExtractView( &localValues, &leadingDim );
	}
	else {
		TEST_FOR_EXCEPT(true); // ToDo: Implement!
	}
  // Build the MultiVector
	Teuchos::RefCountPtr<MPIMultiVectorBase<double> >
    mv = Teuchos::rcp(new DefaultMPIMultiVector<double>(range,domain,Teuchos::rcp(localValues,false),leadingDim));
  Teuchos::set_extra_data( epetra_mv, "Epetra_MultiVector", &mv );
  return mv;
}

Teuchos::RefCountPtr<const Thyra::MPIMultiVectorBase<double> >
Thyra::create_MPIMultiVectorBase(
	const Teuchos::RefCountPtr<const Epetra_MultiVector>                     &epetra_mv
	,const Teuchos::RefCountPtr<const MPIVectorSpaceBase<double> >           &range
	,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<double> >    &domain
	)
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(range.get()==NULL);
  TEST_FOR_EXCEPT(domain.get()==NULL);
#endif
  if(!epetra_mv.get()) return Teuchos::null;
  // New local view of raw data
  double *localValues; int leadingDim;
	if( epetra_mv->ConstantStride() ) {
		epetra_mv->ExtractView( &localValues, &leadingDim );
	}
	else {
		TEST_FOR_EXCEPT(true); // ToDo: Implement!
	}
  // Build the MultiVector
	Teuchos::RefCountPtr<const MPIMultiVectorBase<double> >
    mv = Teuchos::rcp(new DefaultMPIMultiVector<double>(range,domain,Teuchos::rcp(localValues,false),leadingDim));
  Teuchos::set_extra_data( epetra_mv, "Epetra_MultiVector", &mv );
  return mv;
}

Teuchos::RefCountPtr<Epetra_Vector>
Thyra::get_Epetra_Vector(
	const Epetra_Map                                    &map
	,const Teuchos::RefCountPtr<VectorBase<double> >    &v
	)
{
#ifdef _DEBUG
  Teuchos::RefCountPtr<const VectorSpaceBase<double> >
    epetra_vs = create_MPIVectorSpaceBase(Teuchos::rcp(&map,false));
	THYRA_ASSERT_VEC_SPACES( "Thyra::get_Epetra_Vector(map,v)", *epetra_vs, *v->space() );
#endif
  //
  // First, try to grab the Epetra_Vector straight out of the
  // RCP since this is the fastest way.
  //
  const Teuchos::RefCountPtr<Epetra_Vector>
    *epetra_v_ptr = Teuchos::get_optional_extra_data<Teuchos::RefCountPtr<Epetra_Vector> >(v,"Epetra_Vector");
  if(epetra_v_ptr) {
    return *epetra_v_ptr;
  }
	//
	// The assumption that we (rightly) make here is that if the vector spaces
	// are compatible, that either the multi-vectors are both in-core or the
	// vector spaces are both derived from MPIVectorSpaceBase and have
	// compatible maps.
	// 
  const VectorSpaceBase<double>  &vs = *v->range();
  const MPIVectorSpaceBase<double> *mpi_vs = dynamic_cast<const MPIVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
	//
	// Here we will extract a view of the local elements in the underlying
	// VectorBase object.  In most cases, no data will be allocated or copied
	// and only some small objects will be created and a few virtual functions
	// will be called so the overhead should be low and fixed.
	//
	// Create a *mutable* view of the local elements, this view will be set on
	// the RefCountPtr that is returned.  As a result, this view will be relased
	// when the returned Epetra_Vector is released.
	//
	Teuchos::RefCountPtr<DetachedVectorView<double> >
		emvv = Teuchos::rcp(
			new DetachedVectorView<double>(
				v
				,Range1D(localOffset,localOffset+localSubDim-1)
        ,true // forceContiguous
				)
			);
	// Create a temporary Epetra_Vector object and give it
	// the above local view
	Teuchos::RefCountPtr<Epetra_Vector>
		epetra_v = Teuchos::rcp(
			new Epetra_Vector(
				::View                                 // CV
				,map                                   // Map
				,const_cast<double*>(emvv->values())   // V
				)
			);
	// Give the explict view object to the above Epetra_Vector
	// smart pointer object.  In this way, when the client is finished
	// with the Epetra_Vector view the destructor from the object
	// in emvv will automatically commit the changes to the elements in
	// the input v MultiVectorBase object (reguardless of its
	// implementation).  This is truly an elegant result!
	Teuchos::set_extra_data( emvv, "emvv", &epetra_v, Teuchos::PRE_DESTROY );
  // We are done!
	return epetra_v;
}

Teuchos::RefCountPtr<const Epetra_Vector>
Thyra::get_Epetra_Vector(
	const Epetra_Map                                         &map 
	,const Teuchos::RefCountPtr<const VectorBase<double> >   &v
	)
{
#ifdef _DEBUG
  Teuchos::RefCountPtr<const VectorSpaceBase<double> >
    epetra_vs = create_MPIVectorSpaceBase(Teuchos::rcp(&map,false));
	THYRA_ASSERT_VEC_SPACES( "Thyra::get_Epetra_Vector(map,v)", *epetra_vs, *v->space() );
#endif
  //
  // First, try to grab the Epetra_Vector straight out of the
  // RCP since this is the fastest way.
  //
  const Teuchos::RefCountPtr<const Epetra_Vector>
    *epetra_v_ptr = Teuchos::get_optional_extra_data<Teuchos::RefCountPtr<const Epetra_Vector> >(v,"Epetra_Vector");
  if(epetra_v_ptr)
    return *epetra_v_ptr;
  const Teuchos::RefCountPtr<Epetra_Vector>
    *epetra_nonconst_v_ptr = Teuchos::get_optional_extra_data<Teuchos::RefCountPtr<Epetra_Vector> >(v,"Epetra_Vector");
  if(epetra_nonconst_v_ptr)
    return *epetra_nonconst_v_ptr;
	//
	// Same as above function except as stated below
	//
  const VectorSpaceBase<double> &vs = *v->range();
  const MPIVectorSpaceBase<double> *mpi_vs = dynamic_cast<const MPIVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
	// Get an explicit *non-mutable* view of all of the elements in
	// the multi vector.
	Teuchos::RefCountPtr<ConstDetachedVectorView<double> >
		evv = Teuchos::rcp(
			new ConstDetachedVectorView<double>(
				*v
				,Range1D(localOffset,localOffset+localSubDim-1)
        ,true // forceContiguous
				)
			);
	// Create a temporary Epetra_Vector object and give it
	// the above view
	Teuchos::RefCountPtr<Epetra_Vector>
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
	Teuchos::set_extra_data( evv, "evv", &epetra_v, Teuchos::PRE_DESTROY );
  // Also set the v itself as extra data just to be safe
	Teuchos::set_extra_data( v, "v", &epetra_v );
  // We are done!
	return epetra_v;
}

Teuchos::RefCountPtr<Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
	const Epetra_Map                                         &map
	,const Teuchos::RefCountPtr<MultiVectorBase<double> >    &mv
	)
{
#ifdef _DEBUG
  Teuchos::RefCountPtr<const VectorSpaceBase<double> >
    epetra_vs = create_MPIVectorSpaceBase(Teuchos::rcp(&map,false));
	THYRA_ASSERT_VEC_SPACES( "Thyra::get_Epetra_MultiVector(map,mv)", *epetra_vs, *mv->range() );
#endif
  //
  // First, try to grab the Epetra_MultiVector straight out of the
  // RCP since this is the fastest way.
  //
  const Teuchos::RefCountPtr<Epetra_MultiVector>
    *epetra_mv_ptr = Teuchos::get_optional_extra_data<Teuchos::RefCountPtr<Epetra_MultiVector> >(mv,"Epetra_MultiVector");
  if(epetra_mv_ptr) {
    return *epetra_mv_ptr;
  }
	//
	// The assumption that we (rightly) make here is that if the vector spaces
	// are compatible, that either the multi-vectors are both in-core or the
	// vector spaces are both derived from MPIVectorSpaceBase and have
	// compatible maps.
	// 
  const VectorSpaceBase<double> &vs = *mv->range();
  const MPIVectorSpaceBase<double> *mpi_vs = dynamic_cast<const MPIVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
	//
	// Here we will extract a view of the local elements in the underlying
	// MultiVectorBase object.  In most cases, no data will be allocated or
	// copied and only some small objects will be created and a few virtual
	// functions will be called so the overhead should be low and fixed.
	//
	// Create a *mutable* view of the local elements, this view will be set on
	// the RefCountPtr that is returned.  As a result, this view will be relased
	// when the returned Epetra_MultiVector is released.
	//
	Teuchos::RefCountPtr<DetachedMultiVectorView<double> >
		emmvv = Teuchos::rcp(
			new DetachedMultiVectorView<double>(
				*mv
				,Range1D(localOffset,localOffset+localSubDim-1)
				)
			);
	// Create a temporary Epetra_MultiVector object and give it
	// the above local view
	Teuchos::RefCountPtr<Epetra_MultiVector>
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
	Teuchos::set_extra_data( emmvv, "emmvv", &epetra_mv, Teuchos::PRE_DESTROY );
  // Also set the mv itself as extra data just to be safe
	Teuchos::set_extra_data( mv, "mv", &epetra_mv );
  // We are done!
	return epetra_mv;
}

Teuchos::RefCountPtr<const Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
	const Epetra_Map                                              &map
	,const Teuchos::RefCountPtr<const MultiVectorBase<double> >   &mv
	)
{
#ifdef _DEBUG
  Teuchos::RefCountPtr<const VectorSpaceBase<double> >
    epetra_vs = create_MPIVectorSpaceBase(Teuchos::rcp(&map,false));
	THYRA_ASSERT_VEC_SPACES( "Thyra::get_Epetra_MultiVector(map,mv)", *epetra_vs, *mv->range() );
#endif
  //
  // First, try to grab the Epetra_MultiVector straight out of the
  // RCP since this is the fastest way.
  //
  const Teuchos::RefCountPtr<const Epetra_MultiVector>
    *epetra_mv_ptr = Teuchos::get_optional_extra_data<Teuchos::RefCountPtr<const Epetra_MultiVector> >(mv,"Epetra_MultiVector");
  if(epetra_mv_ptr) {
    return *epetra_mv_ptr;
  }
	//
	// Same as above function except as stated below
	//
  const VectorSpaceBase<double> &vs = *mv->range();
  const MPIVectorSpaceBase<double> *mpi_vs = dynamic_cast<const MPIVectorSpaceBase<double>*>(&vs);
  const Index localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Index localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
	// Get an explicit *non-mutable* view of all of the elements in
	// the multi vector.
	Teuchos::RefCountPtr<ConstDetachedMultiVectorView<double> >
		emvv = Teuchos::rcp(
			new ConstDetachedMultiVectorView<double>(
				*mv
				,Range1D(localOffset,localOffset+localSubDim-1)
				)
			);
	// Create a temporary Epetra_MultiVector object and give it
	// the above view
	Teuchos::RefCountPtr<Epetra_MultiVector>
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
	Teuchos::set_extra_data( emvv, "emvv", &epetra_mv, Teuchos::PRE_DESTROY );
  // Also set the mv itself as extra data just to be safe
	Teuchos::set_extra_data( mv, "mv", &epetra_mv );
  // We are done!
	return epetra_mv;
}

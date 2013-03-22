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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
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

//
// Helpers
//


namespace {


Teuchos::RCP<const Thyra::VectorSpaceBase<double> >
unwrapSingleProductVectorSpace(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > &vs_in
  )
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::ProductVectorSpaceBase;
  const RCP<const ProductVectorSpaceBase<double> > pvs =
    rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(vs_in);
  if (nonnull(pvs)) {
    TEUCHOS_ASSERT_EQUALITY( pvs->numBlocks(), 1 );
    return pvs->getBlock(0);
  }
  return vs_in;
}


} // namespace


//
// Implementations of user function
//


Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >
Thyra::create_Comm( const RCP<const Epetra_Comm> &epetraComm )
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

  RCP<const Epetra_SerialComm>
    serialEpetraComm = rcp_dynamic_cast<const Epetra_SerialComm>(epetraComm);
  if( serialEpetraComm.get() ) {
    RCP<const Teuchos::SerialComm<Ordinal> >
      serialComm = rcp(new Teuchos::SerialComm<Ordinal>());
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
    RCP<const Teuchos::MpiComm<Ordinal> >
      mpiComm = rcp(new Teuchos::MpiComm<Ordinal>(rawMpiComm));
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
  using Teuchos::as; using Teuchos::inoutArg;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    !epetra_map.get(), std::invalid_argument,
    "create_VectorSpace::initialize(...): Error!" );
#endif // TEUCHOS_DEBUG
  RCP<const Teuchos::Comm<Ordinal> > comm =
    create_Comm(Teuchos::rcpFromRef(epetra_map->Comm())).assert_not_null();
  Teuchos::set_extra_data(epetra_map, "epetra_map", inoutArg(comm));
  const Ordinal localSubDim = epetra_map->NumMyElements();
  RCP<DefaultSpmdVectorSpace<double> > vs =
    defaultSpmdVectorSpace<double>(
      comm, localSubDim, epetra_map->NumGlobalElements());
  TEUCHOS_ASSERT_EQUALITY(vs->dim(), as<Ordinal>(epetra_map->NumGlobalElements()));
  // NOTE: It is impossible to trigger the above exception unless
  // NumGlobalElemenets() overflows 'int'.  However, this is a nice sanity
  // check to stop the code early in case we seen an overflow in practice.
  // Because this assert will only likely trigger in a non-debug build, we
  // will unguard the assert since it is very cheap to perform.
  Teuchos::set_extra_data( epetra_map, "epetra_map", inoutArg(vs) );
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
  TEUCHOS_TEST_FOR_EXCEPT(parentSpace.get()==NULL);
  Teuchos::dyn_cast<const SpmdVectorSpaceBase<double> >(*parentSpace);
  TEUCHOS_TEST_FOR_EXCEPT(dim <= 0);
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
  TEUCHOS_TEST_FOR_EXCEPT(space_in.get()==NULL);
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
  TEUCHOS_TEST_FOR_EXCEPT(space_in.get()==NULL);
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
  TEUCHOS_TEST_FOR_EXCEPT(range_in.get()==NULL);
#endif
  const RCP<const SpmdVectorSpaceBase<double> > range =
    Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
      unwrapSingleProductVectorSpace(range_in),
      true
      );
  RCP<const ScalarProdVectorSpaceBase<double> > domain =
    Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<double> >(
      unwrapSingleProductVectorSpace(domain_in),
      true
      );
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
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement views of non-contiguous mult-vectors!
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
  TEUCHOS_TEST_FOR_EXCEPT(range_in.get()==NULL);
#endif
  const RCP<const SpmdVectorSpaceBase<double> > range =
    Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
      unwrapSingleProductVectorSpace(range_in),
      true
      );
  RCP<const ScalarProdVectorSpaceBase<double> > domain =
    Teuchos::rcp_dynamic_cast<const ScalarProdVectorSpaceBase<double> >(
      unwrapSingleProductVectorSpace(domain_in),
      true
      );
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
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement views of non-contiguous mult-vectors!
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


Teuchos::RCP<const Epetra_Comm>
Thyra::get_Epetra_Comm(const Teuchos::Comm<Ordinal>& comm_in)
{

  using Teuchos::rcp;
  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;
  using Teuchos::SerialComm;
#ifdef HAVE_MPI
  using Teuchos::MpiComm;
#endif

  const Ptr<const Teuchos::Comm<Ordinal> > comm = Teuchos::ptrFromRef(comm_in);

  const Ptr<const SerialComm<Ordinal> > serialComm =
    ptr_dynamic_cast<const SerialComm<Ordinal> >(comm);

  RCP<const Epetra_Comm> epetraComm;

#ifdef HAVE_MPI

  const Ptr<const MpiComm<Ordinal> > mpiComm =
    ptr_dynamic_cast<const MpiComm<Ordinal> >(comm);

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(mpiComm) && is_null(serialComm),
    std::runtime_error, 
    "SPMD std::vector space has a communicator that is "
    "neither a serial comm nor an MPI comm");

  if (nonnull(mpiComm)) {
    epetraComm = rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()()));
  }
  else {
    epetraComm = rcp(new Epetra_SerialComm());
  }

#else

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(serialComm), std::runtime_error, 
    "SPMD std::vector space has a communicator that is "
    "neither a serial comm nor an MPI comm");

  epetraComm = rcp(new Epetra_SerialComm());
  
#endif
  
  TEUCHOS_TEST_FOR_EXCEPTION(is_null(epetraComm), std::runtime_error,
    "null communicator created");
  
  return epetraComm;

}


Teuchos::RCP<const Epetra_Map>
Thyra::get_Epetra_Map(const VectorSpaceBase<double>& vs_in,
  const RCP<const Epetra_Comm>& comm)
{

  using Teuchos::rcpFromRef;
  using Teuchos::rcpFromPtr;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ptrFromRef;
  using Teuchos::ptr_dynamic_cast;

  const Ptr<const VectorSpaceBase<double> > vs_ptr = ptrFromRef(vs_in);
  
  const Ptr<const SpmdVectorSpaceBase<double> > spmd_vs =
    ptr_dynamic_cast<const SpmdVectorSpaceBase<double> >(vs_ptr);

  const Ptr<const ProductVectorSpaceBase<double> > &prod_vs = 
    ptr_dynamic_cast<const ProductVectorSpaceBase<double> >(vs_ptr);

  TEUCHOS_TEST_FOR_EXCEPTION( is_null(spmd_vs) && is_null(prod_vs), std::logic_error,
    "Error, the concrete VectorSpaceBase object of type "
    +Teuchos::demangleName(typeid(vs_in).name())+" does not support the"
    " SpmdVectorSpaceBase or the ProductVectorSpaceBase interfaces!" );

  const int numBlocks = (nonnull(prod_vs) ? prod_vs->numBlocks() : 1);

  // Get an array of SpmdVectorBase objects for the blocks
  
  Array<RCP<const SpmdVectorSpaceBase<double> > > spmd_vs_blocks;
  if (nonnull(prod_vs)) {
    for (int block_i = 0; block_i < numBlocks; ++block_i) {
      const RCP<const SpmdVectorSpaceBase<double> > spmd_vs_i =
        rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(
          prod_vs->getBlock(block_i), true);
      spmd_vs_blocks.push_back(spmd_vs_i);
    }
  }
  else {
    spmd_vs_blocks.push_back(rcpFromPtr(spmd_vs));
  }
  
  // Find the number of local elements, summed over all blocks

  int myLocalElements = 0;
  for (int block_i = 0; block_i < numBlocks; ++block_i) {
    myLocalElements += spmd_vs_blocks[block_i]->localSubDim();
  }
  
  // Find the GIDs owned by this processor, taken from all blocks
  
  int count=0;
  int blockOffset = 0;
  Array<int> myGIDs(myLocalElements);
  for (int block_i = 0; block_i < numBlocks; ++block_i) {
    const RCP<const SpmdVectorSpaceBase<double> > spmd_vs_i = spmd_vs_blocks[block_i];
    const int lowGIDInBlock = spmd_vs_i->localOffset();
    const int numLocalElementsInBlock = spmd_vs_i->localSubDim();
    for (int i=0; i < numLocalElementsInBlock; ++i, ++count) {
      myGIDs[count] = blockOffset + lowGIDInBlock + i;
    }
    blockOffset += spmd_vs_i->dim();
  }
  
  const int globalDim = vs_in.dim();

  return Teuchos::rcp(
    new Epetra_Map(globalDim, myLocalElements, &(myGIDs[0]), 0, *comm));

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
  const Ordinal localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Ordinal localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
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
  const Ordinal localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Ordinal localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
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
  const Ordinal localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Ordinal localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
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
  const Ordinal localOffset = ( mpi_vs ? mpi_vs->localOffset() : 0 );
  const Ordinal localSubDim = ( mpi_vs ? mpi_vs->localSubDim() : vs.dim() );
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
  using Teuchos::rcpWithEmbeddedObj;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  ArrayRCP<double> mvData;
  Ordinal mvLeadingDim = -1;
  SpmdMultiVectorBase<double> *mvSpmdMv = 0;
  SpmdVectorBase<double> *mvSpmdV = 0;
  if ((mvSpmdMv = dynamic_cast<SpmdMultiVectorBase<double>*>(&mv))) {
    mvSpmdMv->getNonconstLocalData(outArg(mvData), outArg(mvLeadingDim));
  }
  else if ((mvSpmdV = dynamic_cast<SpmdVectorBase<double>*>(&mv))) {
    mvSpmdV->getNonconstLocalData(outArg(mvData));
    mvLeadingDim = mvSpmdV->spmdSpace()->localSubDim();
  }
  if (nonnull(mvData)) {
    return rcpWithEmbeddedObj(
      new Epetra_MultiVector(
        ::View, map, mvData.getRawPtr(), mvLeadingDim, mv.domain()->dim()
        ),
      mvData
      );
  }
  return ::Thyra::get_Epetra_MultiVector(map, rcpFromRef(mv));
}


Teuchos::RCP<const Epetra_MultiVector>
Thyra::get_Epetra_MultiVector(
  const Epetra_Map &map,
  const MultiVectorBase<double> &mv
  )
{
  using Teuchos::rcpWithEmbeddedObj;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  ArrayRCP<const double> mvData;
  Ordinal mvLeadingDim = -1;
  const SpmdMultiVectorBase<double> *mvSpmdMv = 0;
  const SpmdVectorBase<double> *mvSpmdV = 0;
  if ((mvSpmdMv = dynamic_cast<const SpmdMultiVectorBase<double>*>(&mv))) {
    mvSpmdMv->getLocalData(outArg(mvData), outArg(mvLeadingDim));
  }
  else if ((mvSpmdV = dynamic_cast<const SpmdVectorBase<double>*>(&mv))) {
    mvSpmdV->getLocalData(outArg(mvData));
    mvLeadingDim = mvSpmdV->spmdSpace()->localSubDim();
  }
  if (nonnull(mvData)) {
    return rcpWithEmbeddedObj(
      new Epetra_MultiVector(
        ::View,map, const_cast<double*>(mvData.getRawPtr()), mvLeadingDim, mv.domain()->dim()
        ),
      mvData
      );
  }
  return ::Thyra::get_Epetra_MultiVector(map, rcpFromRef(mv));
}

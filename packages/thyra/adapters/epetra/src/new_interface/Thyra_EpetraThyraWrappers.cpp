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

#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

#include "Thyra_EpetraVectorSpace.hpp"
#include "Thyra_EpetraMultiVector.hpp"
#include "Thyra_EpetraVector.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"

#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

//
// Helpers
//

namespace {


Teuchos::RCP<const Thyra::EpetraVectorSpace>
getOrCreateEpetraVectorSpace(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space,
  const Epetra_BlockMap& epetraBlockMap
  )
{
  Teuchos::RCP<const Thyra::EpetraVectorSpace> epetraSpace;
  if (Teuchos::nonnull(space)) {
    epetraSpace = Teuchos::rcp_dynamic_cast<const Thyra::EpetraVectorSpace>(space, true);
  } else {
    // Epetra maps use an in-house reference counting, and hide all data in 'Epetra_BlockMapData',
    // so making a copy is inexpensive
    Teuchos::RCP<const Epetra_Map> epetraMapRCP( new Epetra_Map(static_cast<const Epetra_Map&>(epetraBlockMap)) );
    epetraSpace = Thyra::epetraVectorSpace(epetraMapRCP);
  }
  return epetraSpace;
}


Teuchos::RCP<const Thyra::ScalarProdVectorSpaceBase<double>>
getOrCreateLocallyReplicatedEpetraVectorSpace(
  const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> space,
  const Epetra_Comm& epetraComm,
  const int numCols
  )
{
  Teuchos::RCP<const Thyra::EpetraVectorSpace> epetraSpace;
  if (Teuchos::nonnull(space)) {
    epetraSpace = Teuchos::rcp_dynamic_cast<const Thyra::EpetraVectorSpace>(space, true);
  }
  else {
    epetraSpace = Thyra::epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(numCols,0,epetraComm) ));
  }
  return epetraSpace;
}

} // anonymous namespace

namespace Thyra
{

Teuchos::RCP<const Teuchos::Comm<Ordinal>>
convertEpetraToThyraComm( const Epetra_Comm& epetraComm )
{
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<const Epetra_Comm> epetraCommPtr = Teuchos::rcpFromRef(epetraComm);

#ifdef HAVE_MPI
  const Teuchos::RCP<const Epetra_MpiComm> epetraMpiComm = 
    rcp_dynamic_cast<const Epetra_MpiComm>(epetraCommPtr);
  if (nonnull(epetraMpiComm)) {
    return Teuchos::createMpiComm<Ordinal>(Teuchos::opaqueWrapper(epetraMpiComm->Comm()));
  }
#endif // HAVE_MPI

  // Assert conversion to Epetra_SerialComm as a last resort (or throw)
  rcp_dynamic_cast<const Epetra_SerialComm>(epetraCommPtr, true);

  return Teuchos::createSerialComm<Ordinal>();

  // NOTE: Above will throw if the type is not Epetra_SerialComm.  In this
  // case, the type could not be converted. 
}


Teuchos::RCP<const VectorSpaceBase<double>>
createVectorSpace(
  const RCP<const Epetra_Map>& epetraMap
  )
{
  return epetraVectorSpace(epetraMap);
}


Teuchos::RCP<VectorBase<double>>
createVector(
  const RCP<Epetra_Vector>& epetraVector_in,
  const RCP<const VectorSpaceBase<double>> space_in
  )
{
  return epetraVector(
    getOrCreateEpetraVectorSpace(space_in, epetraVector_in->Map()),
    epetraVector_in
    );
}


Teuchos::RCP<const VectorBase<double>>
createConstVector(
  const RCP<const Epetra_Vector>& epetraVector_in,
  const RCP<const VectorSpaceBase<double>> space
  )
{
  return constEpetraVector(
    getOrCreateEpetraVectorSpace(space, epetraVector_in->Map()),
    epetraVector_in
    );
}


Teuchos::RCP<MultiVectorBase<double>>
createMultiVector(
  const RCP<Epetra_MultiVector>& epetraMultiVector_in,
  const RCP<const VectorSpaceBase<double>> rangeSpace,
  const RCP<const VectorSpaceBase<double>> domainSpace
  )
{
  return epetraMultiVector(
    getOrCreateEpetraVectorSpace(rangeSpace, epetraMultiVector_in->Map()),
    getOrCreateLocallyReplicatedEpetraVectorSpace(
      domainSpace, epetraMultiVector_in->Comm(),
      epetraMultiVector_in->NumVectors()
      ),
    epetraMultiVector_in
    );
}


Teuchos::RCP<const MultiVectorBase<double>>
createConstMultiVector(
  const RCP<const Epetra_MultiVector>& epetraMultiVector_in,
  const RCP<const VectorSpaceBase<double>> rangeSpace,
  const RCP<const VectorSpaceBase<double>> domainSpace
  )
{
  return constEpetraMultiVector(
    getOrCreateEpetraVectorSpace(rangeSpace, epetraMultiVector_in->Map()),
    getOrCreateLocallyReplicatedEpetraVectorSpace(
      domainSpace, epetraMultiVector_in->Comm(),
      epetraMultiVector_in->NumVectors()
      ),
    epetraMultiVector_in
    );
}


Teuchos::RCP<LinearOpBase<double>>
createLinearOp(
  const RCP<Epetra_Operator>& epetraOperator_in,
  const RCP<const VectorSpaceBase<double>> rangeSpace,
  const RCP<const VectorSpaceBase<double>> domainSpace
  )
{
  return epetraLinearOp(
    getOrCreateEpetraVectorSpace(rangeSpace, epetraOperator_in->OperatorRangeMap()),
    getOrCreateEpetraVectorSpace(domainSpace, epetraOperator_in->OperatorDomainMap()),
    epetraOperator_in
    );
}


Teuchos::RCP<const LinearOpBase<double>>
createConstLinearOp(
  const RCP<const Epetra_Operator> &epetraOperator_in,
  const RCP<const VectorSpaceBase<double>> rangeSpace,
  const RCP<const VectorSpaceBase<double>> domainSpace
  )
{
  return constEpetraLinearOp(
    getOrCreateEpetraVectorSpace(rangeSpace, epetraOperator_in->OperatorRangeMap()),
    getOrCreateEpetraVectorSpace(domainSpace, epetraOperator_in->OperatorDomainMap()),
    epetraOperator_in
    );
}



RCP<Epetra_Vector>
EpetraOperatorVectorExtraction::
getEpetraVector(const RCP<VectorBase<double>> &v)
{
  return Teuchos::rcp_dynamic_cast<EpetraVector>(v, true)->getEpetraVector();
}


RCP<const Epetra_Vector>
EpetraOperatorVectorExtraction::
getConstEpetraVector(const RCP<const VectorBase<double>> &v)
{
  return Teuchos::rcp_dynamic_cast<const EpetraVector>(v, true)->getConstEpetraVector();
}


RCP<Epetra_MultiVector>
EpetraOperatorVectorExtraction::
getEpetraMultiVector(const RCP<MultiVectorBase<double>> &mv)
{

#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(mv));
#endif

  using Teuchos::rcp_dynamic_cast;
  
  const RCP<EpetraMultiVector> emv =
    rcp_dynamic_cast<EpetraMultiVector>(mv);
  if (nonnull(emv)) {
    return emv->getEpetraMultiVector();
  }
  
  const RCP<EpetraVector> ev =
    rcp_dynamic_cast<EpetraVector>(mv);
  if (nonnull(ev)) {
    return ev->getEpetraVector();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error, the input mv = " << mv->description() << " does not support the"
    " EpetraMultiVector or the EpetraVector interfaces!");

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);

}


RCP<const Epetra_MultiVector>
EpetraOperatorVectorExtraction::
getConstEpetraMultiVector(const RCP<const MultiVectorBase<double>> &mv)
{

#ifdef THYRA_DEBUG
  TEUCHOS_ASSERT(nonnull(mv));
#endif

  using Teuchos::rcp_dynamic_cast;
  
  const RCP<const EpetraMultiVector> emv =
    rcp_dynamic_cast<const EpetraMultiVector>(mv);
  if (nonnull(emv)) {
    return emv->getConstEpetraMultiVector();
  }
  
  const RCP<const EpetraVector> ev =
    rcp_dynamic_cast<const EpetraVector>(mv);
  if (nonnull(ev)) {
    return ev->getConstEpetraVector();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error, the input mv = " << mv->description() << " does not support the"
    " EpetraMultiVector or the EpetraVector interfaces!");

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);

}


RCP<Epetra_Operator>
EpetraOperatorVectorExtraction::
getEpetraOperator(const RCP<LinearOpBase<double>> &op)
{
  return Teuchos::rcp_dynamic_cast<EpetraLinearOp>(op, true)->getEpetraOperator();
}


RCP<const Epetra_Operator>
EpetraOperatorVectorExtraction::
getConstEpetraOperator(const RCP<const LinearOpBase<double>> &op)
{
  return Teuchos::rcp_dynamic_cast<const EpetraLinearOp>(op, true)->getConstEpetraOperator();
}

///////////////////////////////////////////////////
//             getOrCreate methods               //
///////////////////////////////////////////////////

RCP<const Epetra_Map>
EpetraOperatorVectorExtraction::
getOrCreateEpetraMap (const RCP<const VectorSpaceBase<double>> &vs)
{
  // Best case scenario: input vs is an EpetraVectorSpace
  auto tmp = Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(vs,false);
  if (!tmp.is_null()) {
    return tmp->getEpetraMap();
  }

  // Bad luck. We need to create an Epetra map given the vs information.
  // There are two acceptable cases: vs is an SpmdVectorSpaceBase, or
  // vs is a ProductVectorSpaceBase, each block being an SpmdVectorSpaceBase.
  auto spmd_vs = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double>>(vs,false);
  auto prod_vs = Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double>>(vs,false);

  // Check if we are out of luck
  TEUCHOS_TEST_FOR_EXCEPTION(spmd_vs.is_null() && prod_vs.is_null(), std::runtime_error,
                             "Error! Input vector space is neither SpmdVectorSpaceBase nor ProductVectorSpaceBase.\n");

  // NOTE: we are using 32bits indices!
  using GO = int;

  int numBlocks = 0;
  Teuchos::Array<Teuchos::RCP<const SpmdVectorSpaceBase<double>>> blocks;
  int numMyElements = 0;
  if (!spmd_vs.is_null()) {
    numBlocks = 1;
    blocks.resize(1,spmd_vs);
    numMyElements = spmd_vs->localSubDim();
  } else {
    numBlocks = prod_vs->numBlocks();
    TEUCHOS_TEST_FOR_EXCEPTION(numBlocks<1, std::logic_error, "Error! There are no blocks in the Product vector space.\n");
    blocks.resize(numBlocks);
    for (int i=0; i<numBlocks; ++i) {
      blocks[i] = Teuchos::rcp_dynamic_cast<const SpmdVectorSpaceBase<double>>(prod_vs->getBlock(i),true);
      numMyElements += blocks[i]->localSubDim();
    }
  }

  // Ok, now we know how many elements this rank owns. Now it's time to fill the GIDs.
  // NOTE: we assume that, globally, the global map will contain all the gids of block 0,
  //       then those of block 1, ..., then those of block N. This means that the GIDs
  //       on this rank may not be contiguous. For instance, if there are 2 ranks, and
  //       2 blocks in the Spmd VS, each of which has, globally, 4 elements, linearly
  //       distributed, then the owned GIDs will be:
  //          p0: [0, 1, 4, 5]
  //          p1: [2, 3, 6, 7]
  //       where 0-3 are the GIDs of the 1st block, while 4-7 are the GIDs of the 2nd block.
  Teuchos::Array<GO> indices(numMyElements);
  GO offset = 0;
  for (int iblock=0, k=0; iblock<numBlocks; ++iblock) {
    const GO lowest_gid_in_block = blocks[iblock]->localOffset();
    const int localDim = blocks[iblock]->localSubDim();
    for (int i=0; i<localDim; ++i, ++k) {
      indices[k] = offset + i + lowest_gid_in_block;
    }
    offset += blocks[iblock]->dim();
  }

  // Create an Epetra_Comm
  RCP<const Epetra_Comm> comm = createEpetraComm(blocks[0]->getComm());

  // Now we can create the map
  // NOTE: the map's constructor clones the comm, so don't worry about comm going out of scope
  // NOTE: we don't take the global dimension from the input vs, in case the map is overlapped.
  //       Yes, doing vector operations like norm/dot on overlapped partitions is wrong, but
  //       we don't know what the vector was meant to be for, so maybe it is correct.
  //       If you really are concerned about this, you should put an error in the Norm/Dot methods
  //       of vector and multivector.
  return Teuchos::rcp( new Epetra_Map(-1,numMyElements,indices().getConst().getRawPtr(),0,*comm) );
}

RCP<Epetra_Vector>
EpetraOperatorVectorExtraction::
getOrCreateEpetraVector(const RCP<VectorBase<double> > &v)
{
  // Best case scenario: input vs is an EpetraVector
  auto tmp = Teuchos::rcp_dynamic_cast<EpetraVector>(v);
  if (!tmp.is_null()) {
    return tmp->getEpetraVector();
  }

  // Bad luck. First, create an Epetra_Map
  auto map = getOrCreateEpetraMap(v->space());

  // Then create a EpetraVectorSpace
  RCP<const EpetraVectorSpace> evs = epetraVectorSpace(map);

  // Then create a EpetraVector
  RCP<EpetraVector> ev = epetraVector(evs, Teuchos::rcp( new Epetra_Vector(*map) ));

  // Then use assign to copy values (will rely on RTOp)
  ev->assign(*v);

  return ev->getEpetraVector();
}

RCP<const Epetra_Vector>
EpetraOperatorVectorExtraction::
getOrCreateConstEpetraVector(const RCP<const VectorBase<double> > &v)
{
  // Best case scenario: input vs is an EpetraVector
  auto tmp = Teuchos::rcp_dynamic_cast<const EpetraVector>(v);
  if (!tmp.is_null()) {
    return tmp->getConstEpetraVector();
  }

  // Bad luck. First, create an Epetra_Map
  auto map = getOrCreateEpetraMap(v->space());

  // Then create a EpetraVectorSpace
  RCP<const EpetraVectorSpace> evs = epetraVectorSpace(map);

  // Then create a EpetraVector
  RCP<EpetraVector> ev = epetraVector(evs, Teuchos::rcp( new Epetra_Vector(*map) ));

  // Then use assign to copy values (will rely on RTOp)
  ev->assign(*v);

  return ev->getConstEpetraVector();
}

RCP<Epetra_MultiVector>
EpetraOperatorVectorExtraction::
getOrCreateEpetraMultiVector(const RCP<MultiVectorBase<double> > &mv)
{
  // Best case scenario: input vs is an EpetraMultiVector
  auto tmp = Teuchos::rcp_dynamic_cast<EpetraMultiVector>(mv);
  if (!tmp.is_null()) {
    return tmp->getEpetraMultiVector();
  }

  // Bad luck. First, create an Epetra_Map
  auto map = getOrCreateEpetraMap(mv->range());

  // Then create the range and domain EpetraVectorSpace
  RCP<const EpetraVectorSpace> range  = epetraVectorSpace(map);
  RCP<const EpetraVectorSpace> domain = epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(static_cast<int>(mv->domain()->dim()),0,map->Comm()) ));

  // Then create a EpetraMultiVector
  RCP<EpetraMultiVector> emv = epetraMultiVector(range, domain, Teuchos::rcp( new Epetra_MultiVector(*map,mv->domain()->dim()) ));

  // Then use assign to copy values (will rely on RTOp)
  emv->assign(*mv);

  return emv->getEpetraMultiVector();
}

RCP<const Epetra_MultiVector>
EpetraOperatorVectorExtraction::
getOrCreateConstEpetraMultiVector(const RCP<const MultiVectorBase<double> > &mv)
{
  // Best case scenario: input vs is an EpetraMultiVector
  auto tmp = Teuchos::rcp_dynamic_cast<const EpetraMultiVector>(mv);
  if (!tmp.is_null()) {
    return tmp->getConstEpetraMultiVector();
  }

  // Bad luck. First, create an Epetra_Map
  auto map = getOrCreateEpetraMap(mv->range());

  // Then create the range and domain EpetraVectorSpace
  RCP<const EpetraVectorSpace> range  = epetraVectorSpace(map);
  RCP<const EpetraVectorSpace> domain = epetraVectorSpace(Teuchos::rcp( new Epetra_LocalMap(static_cast<int>(mv->domain()->dim()),0,map->Comm()) ));

  // Then create a EpetraMultiVector
  RCP<EpetraMultiVector> emv = epetraMultiVector(range, domain, Teuchos::rcp( new Epetra_MultiVector(*map,mv->domain()->dim()) ));

  // Then use assign to copy values (will rely on RTOp)
  emv->assign(*mv);

  return emv->getConstEpetraMultiVector();
}

RCP<Epetra_Operator>
EpetraOperatorVectorExtraction::
getOrCreateEpetraOperator(const RCP<LinearOpBase<double> > &op)
{
  // Best case scenario: input vs is an EpetraLinearOp
  auto tmp = Teuchos::rcp_dynamic_cast<EpetraLinearOp>(op);
  if (!tmp.is_null()) {
    return tmp->getEpetraOperator();
  }

  // Bad luck. Wrap inside EpetraOperatorWrapper
  return Teuchos::rcp( new EpetraOperatorWrapper(op) );
}

RCP<const Epetra_Operator>
EpetraOperatorVectorExtraction::
getOrCreateConstEpetraOperator(const RCP<const LinearOpBase<double> > &op)
{
  // Best case scenario: input vs is an EpetraLinearOp
  auto tmp = Teuchos::rcp_dynamic_cast<const EpetraLinearOp>(op);
  if (!tmp.is_null()) {
    return tmp->getConstEpetraOperator();
  }

  // Bad luck. Wrap inside EpetraOperatorWrapper
  return Teuchos::rcp( new EpetraOperatorWrapper(op) );
}

RCP<const Epetra_Comm>
EpetraOperatorVectorExtraction::
createEpetraComm(const RCP<const Teuchos::Comm<Ordinal>> comm)
{
  RCP<const Epetra_Comm> epetraComm;
#ifdef HAVE_MPI
  const RCP<const Teuchos::SerialComm<Ordinal>> serialComm =
    Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<Ordinal>>(comm);

  const RCP<const Teuchos::MpiComm<Ordinal>> mpiComm =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<Ordinal>>(comm);

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(mpiComm) && is_null(serialComm),
    std::runtime_error,
    "SPMD std::vector space has a communicator that is "
    "neither a serial comm nor an MPI comm");

  if (nonnull(mpiComm)) {
    epetraComm = Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()()));
  }
  else {
    epetraComm = Teuchos::rcp(new Epetra_SerialComm());
  }

#else
  const RCP<const Teuchos::SerialComm<Ordinal>> serialComm =
    Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<Ordinal>>(comm);

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(serialComm), std::runtime_error,
    "SPMD std::vector space has a communicator that is "
    "not a serial comm and MPI is not enabled (so I can't check)");

  epetraComm = Teuchos::rcp(new Epetra_SerialComm());

#endif

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(epetraComm), std::runtime_error,
    "null communicator created");

  return epetraComm;
}

} // namespace Thyra

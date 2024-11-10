// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_EpetraThyraConverter.hpp"

// Teuchos includes
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"

// Thyra includes
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include <iostream>
#include <vector>

using Teuchos::null;
using Teuchos::Ptr;
using Teuchos::ptr_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

namespace Teko {
namespace Epetra {

// const Teuchos::RCP<const Thyra::MultiVectorBase<double> >
// blockEpetraToThyra(int numVectors,const double * epetraData,int leadingDim,const
// Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs,int & localDim)

void blockEpetraToThyra(int numVectors, const double* epetraData, int leadingDim,
                        const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& mv, int& localDim) {
  localDim = 0;

  // check the base case
  const Ptr<Thyra::ProductMultiVectorBase<double> > prodMV =
      ptr_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(mv);
  if (prodMV == Teuchos::null) {
    // VS object must be a SpmdMultiVector object
    const Ptr<Thyra::SpmdMultiVectorBase<double> > spmdX =
        ptr_dynamic_cast<Thyra::SpmdMultiVectorBase<double> >(mv, true);
    const RCP<const Thyra::SpmdVectorSpaceBase<double> > spmdVS = spmdX->spmdSpace();

    int localSubDim = spmdVS->localSubDim();

    Thyra::Ordinal thyraLeadingDim = 0;
    // double * thyraData=0;
    // spmdX->getLocalData(&thyraData,&thyraLeadingDim);

    Teuchos::ArrayRCP<double> thyraData_arcp;
    Teuchos::ArrayView<double> thyraData;
    spmdX->getNonconstLocalData(Teuchos::outArg(thyraData_arcp), Teuchos::outArg(thyraLeadingDim));
    thyraData = thyraData_arcp();  // build array view

    for (int i = 0; i < localSubDim; i++) {
      // copy each vector
      for (int v = 0; v < numVectors; v++)
        thyraData[i + thyraLeadingDim * v] = epetraData[i + leadingDim * v];
    }

    // spmdX->commitLocalData(&thyraData[0]);

    // set the local dimension
    localDim = localSubDim;

    return;
  }

  // this keeps track of current location in the epetraData vector
  const double* localData = epetraData;

  // loop over all the blocks in the vector space
  for (int blkIndex = 0; blkIndex < prodMV->productSpace()->numBlocks(); blkIndex++) {
    int subDim = 0;
    const RCP<Thyra::MultiVectorBase<double> > blockVec =
        prodMV->getNonconstMultiVectorBlock(blkIndex);

    // perorm the recusive copy
    blockEpetraToThyra(numVectors, localData, leadingDim, blockVec.ptr(), subDim);

    // shift to the next block
    localData += subDim;

    // account for the size of this subblock
    localDim += subDim;
  }
}

// Convert a Epetra_MultiVector with assumed block structure dictated by the
// vector space into a Thyra::MultiVectorBase object.
// const Teuchos::RCP<const Thyra::MultiVectorBase<double> > blockEpetraToThyra(const
// Epetra_MultiVector & e,const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs)
void blockEpetraToThyra(const Epetra_MultiVector& epetraX,
                        const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& thyraX) {
  TEUCHOS_ASSERT(thyraX->range()->dim() == epetraX.GlobalLength());

  // extract local information from the Epetra_MultiVector
  int leadingDim = 0, numVectors = 0, localDim = 0;
  double* epetraData = 0;
  epetraX.ExtractView(&epetraData, &leadingDim);

  numVectors = epetraX.NumVectors();

  blockEpetraToThyra(numVectors, epetraData, leadingDim, thyraX.ptr(), localDim);

  TEUCHOS_ASSERT(localDim == epetraX.MyLength());
}

void blockThyraToEpetra(int numVectors, double* epetraData, int leadingDim,
                        const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& tX,
                        int& localDim) {
  localDim = 0;

  // check the base case
  const RCP<const Thyra::ProductMultiVectorBase<double> > prodX =
      rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(tX);
  if (prodX == Teuchos::null) {
    // the base case

    // VS object must be a SpmdMultiVector object
    RCP<const Thyra::SpmdMultiVectorBase<double> > spmdX =
        rcp_dynamic_cast<const Thyra::SpmdMultiVectorBase<double> >(tX, true);
    RCP<const Thyra::SpmdVectorSpaceBase<double> > spmdVS = spmdX->spmdSpace();

    int localSubDim = spmdVS->localSubDim();

    Thyra::Ordinal thyraLeadingDim = 0;
    // const double * thyraData=0;
    // spmdX->getLocalData(&thyraData,&thyraLeadingDim);

    Teuchos::ArrayView<const double> thyraData;
    Teuchos::ArrayRCP<const double> thyraData_arcp;
    spmdX->getLocalData(Teuchos::outArg(thyraData_arcp), Teuchos::outArg(thyraLeadingDim));
    thyraData = thyraData_arcp();  // grab the array view

    for (int i = 0; i < localSubDim; i++) {
      // copy each vector
      for (int v = 0; v < numVectors; v++)
        epetraData[i + leadingDim * v] = thyraData[i + thyraLeadingDim * v];
    }

    // set the local dimension
    localDim = localSubDim;

    return;
  }

  const RCP<const Thyra::ProductVectorSpaceBase<double> > prodVS = prodX->productSpace();

  // this keeps track of current location in the epetraData vector
  double* localData = epetraData;

  // loop over all the blocks in the vector space
  for (int blkIndex = 0; blkIndex < prodVS->numBlocks(); blkIndex++) {
    int subDim = 0;

    // construct the block vector
    blockThyraToEpetra(numVectors, localData, leadingDim, prodX->getMultiVectorBlock(blkIndex),
                       subDim);

    // shift to the next block
    localData += subDim;

    // account for the size of this subblock
    localDim += subDim;
  }

  return;
}

// Convert a Thyra::MultiVectorBase object to a Epetra_MultiVector object with
// the map defined by the Epetra_Map.
// const Teuchos::RCP<const Epetra_MultiVector>
// blockThyraToEpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,const RCP<const
// Epetra_Map> & map)
void blockThyraToEpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& thyraX,
                        Epetra_MultiVector& epetraX) {
  // build an Epetra_MultiVector object
  int numVectors = thyraX->domain()->dim();

  // make sure the number of vectors are the same
  TEUCHOS_ASSERT(numVectors == epetraX.NumVectors());
  TEUCHOS_ASSERT(thyraX->range()->dim() == epetraX.GlobalLength());

  // extract local information from the Epetra_MultiVector
  int leadingDim = 0, localDim = 0;
  double* epetraData = 0;
  epetraX.ExtractView(&epetraData, &leadingDim);

  // perform recursive copy
  blockThyraToEpetra(numVectors, epetraData, leadingDim, thyraX, localDim);

  // sanity check
  TEUCHOS_ASSERT(localDim == epetraX.Map().NumMyElements());
}

void thyraVSToEpetraMap(std::vector<int>& myIndicies, int blockOffset,
                        const Thyra::VectorSpaceBase<double>& vs, int& localDim) {
  // zero out set local dimension
  localDim = 0;

  const RCP<const Thyra::ProductVectorSpaceBase<double> > prodVS =
      rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(rcpFromRef(vs));

  // is more recursion needed?
  if (prodVS == Teuchos::null) {
    // base case

    // try to cast to an SPMD capable vector space
    const RCP<const Thyra::SpmdVectorSpaceBase<double> > spmdVS =
        rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(rcpFromRef(vs));
    TEUCHOS_TEST_FOR_EXCEPTION(spmdVS == Teuchos::null, std::runtime_error,
                               "thyraVSToEpetraMap requires all subblocks to be SPMD");

    // get local data storage information
    int localOffset = spmdVS->localOffset();
    int localSubDim = spmdVS->localSubDim();

    // add indicies to matrix
    for (int i = 0; i < localSubDim; i++) myIndicies.push_back(blockOffset + localOffset + i);

    localDim += localSubDim;

    return;
  }

  // loop over all the blocks in the vector space
  for (int blkIndex = 0; blkIndex < prodVS->numBlocks(); blkIndex++) {
    int subDim = 0;

    // construct the block vector
    thyraVSToEpetraMap(myIndicies, blockOffset, *prodVS->getBlock(blkIndex), subDim);

    blockOffset += prodVS->getBlock(blkIndex)->dim();

    // account for the size of this subblock
    localDim += subDim;
  }
}

// From a Thyra vector space create a compatable Epetra_Map
const RCP<Epetra_Map> thyraVSToEpetraMap(const Thyra::VectorSpaceBase<double>& vs,
                                         const RCP<const Epetra_Comm>& comm) {
  int localDim = 0;
  std::vector<int> myGIDs;

  // call recursive routine that constructs the mapping
  thyraVSToEpetraMap(myGIDs, 0, vs, localDim);

  TEUCHOS_ASSERT(myGIDs.size() == (unsigned int)localDim);

  // create the map
  return rcp(new Epetra_Map(vs.dim(), myGIDs.size(), &(myGIDs[0]), 0, *comm));
}

}  // end namespace Epetra
}  // end namespace Teko

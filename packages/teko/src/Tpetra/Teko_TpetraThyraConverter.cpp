// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_TpetraThyraConverter.hpp"
#include "Tpetra_Core.hpp"

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
namespace TpetraHelpers {

// const Teuchos::RCP<const Thyra::MultiVectorBase<double> >
// blockTpetraToThyra(int numVectors,const double * tpetraData,int leadingDim,const
// Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs,int & localDim)

void blockTpetraToThyra(int numVectors, Teuchos::ArrayRCP<const ST> tpetraData, int leadingDim,
                        const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& mv, int& localDim) {
  localDim = 0;

  // check the base case
  const Ptr<Thyra::ProductMultiVectorBase<ST> > prodMV =
      ptr_dynamic_cast<Thyra::ProductMultiVectorBase<ST> >(mv);
  if (prodMV == Teuchos::null) {
    // VS object must be a SpmdMultiVector object
    const Ptr<Thyra::SpmdMultiVectorBase<ST> > spmdX =
        ptr_dynamic_cast<Thyra::SpmdMultiVectorBase<ST> >(mv, true);
    const RCP<const Thyra::SpmdVectorSpaceBase<ST> > spmdVS = spmdX->spmdSpace();

    int localSubDim = spmdVS->localSubDim();

    Thyra::Ordinal thyraLeadingDim = 0;

    Teuchos::ArrayRCP<ST> thyraData_arcp;
    Teuchos::ArrayView<ST> thyraData;
    spmdX->getNonconstLocalData(Teuchos::outArg(thyraData_arcp), Teuchos::outArg(thyraLeadingDim));
    thyraData = thyraData_arcp();  // build array view

    for (int i = 0; i < localSubDim; i++) {
      // copy each vector
      for (int v = 0; v < numVectors; v++)
        thyraData[i + thyraLeadingDim * v] = tpetraData[i + leadingDim * v];
    }

    // set the local dimension
    localDim = localSubDim;

    return;
  }

  // this keeps track of current location in the tpetraData vector
  Teuchos::ArrayRCP<const ST> localData = tpetraData;

  // loop over all the blocks in the vector space
  for (int blkIndex = 0; blkIndex < prodMV->productSpace()->numBlocks(); blkIndex++) {
    int subDim                                      = 0;
    const RCP<Thyra::MultiVectorBase<ST> > blockVec = prodMV->getNonconstMultiVectorBlock(blkIndex);

    // perform the recusive copy
    blockTpetraToThyra(numVectors, localData, leadingDim, blockVec.ptr(), subDim);

    // shift to the next block
    localData += subDim;

    // account for the size of this subblock
    localDim += subDim;
  }
}

// Convert a Tpetra_MultiVector with assumed block structure dictated by the
// vector space into a Thyra::MultiVectorBase object.
// const Teuchos::RCP<const Thyra::MultiVectorBase<double> > blockTpetraToThyra(const
// Tpetra_MultiVector & e,const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs)
void blockTpetraToThyra(const Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX,
                        const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyraX) {
  TEUCHOS_ASSERT((Tpetra::global_size_t)thyraX->range()->dim() == tpetraX.getGlobalLength());

  // extract local information from the Tpetra_MultiVector
  LO leadingDim = 0, localDim = 0;
  leadingDim                             = tpetraX.getStride();
  Teuchos::ArrayRCP<const ST> tpetraData = tpetraX.get1dView();

  int numVectors = tpetraX.getNumVectors();

  blockTpetraToThyra(numVectors, tpetraData, leadingDim, thyraX.ptr(), localDim);

  TEUCHOS_ASSERT((size_t)localDim == tpetraX.getLocalLength());
}

void blockThyraToTpetra(LO numVectors, Teuchos::ArrayRCP<ST> tpetraData, LO leadingDim,
                        const Teuchos::RCP<const Thyra::MultiVectorBase<ST> >& tX, LO& localDim) {
  localDim = 0;

  // check the base case
  const RCP<const Thyra::ProductMultiVectorBase<ST> > prodX =
      rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<ST> >(tX);
  if (prodX == Teuchos::null) {
    // the base case

    // VS object must be a SpmdMultiVector object
    RCP<const Thyra::SpmdMultiVectorBase<ST> > spmdX =
        rcp_dynamic_cast<const Thyra::SpmdMultiVectorBase<ST> >(tX, true);
    RCP<const Thyra::SpmdVectorSpaceBase<ST> > spmdVS = spmdX->spmdSpace();

    Thyra::Ordinal thyraLeadingDim = 0;
    Teuchos::ArrayView<const ST> thyraData;
    Teuchos::ArrayRCP<const ST> thyraData_arcp;
    spmdX->getLocalData(Teuchos::outArg(thyraData_arcp), Teuchos::outArg(thyraLeadingDim));
    thyraData = thyraData_arcp();  // grab the array view

    LO localSubDim = spmdVS->localSubDim();
    for (LO i = 0; i < localSubDim; i++) {
      // copy each vector
      for (LO v = 0; v < numVectors; v++) {
        tpetraData[i + leadingDim * v] = thyraData[i + thyraLeadingDim * v];
      }
    }

    // set the local dimension
    localDim = localSubDim;

    return;
  }

  const RCP<const Thyra::ProductVectorSpaceBase<ST> > prodVS = prodX->productSpace();

  // this keeps track of current location in the tpetraData vector
  Teuchos::ArrayRCP<ST> localData = tpetraData;

  // loop over all the blocks in the vector space
  for (int blkIndex = 0; blkIndex < prodVS->numBlocks(); blkIndex++) {
    int subDim = 0;

    // construct the block vector
    blockThyraToTpetra(numVectors, localData, leadingDim, prodX->getMultiVectorBlock(blkIndex),
                       subDim);

    // shift to the next block
    localData += subDim;

    // account for the size of this subblock
    localDim += subDim;
  }

  return;
}

// Convert a Thyra::MultiVectorBase object to a Tpetra_MultiVector object with
// the map defined by the Tpetra_Map.
// const Teuchos::RCP<const Tpetra_MultiVector>
// blockThyraToTpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & tX,const RCP<const
// Tpetra_Map> & map)
void blockThyraToTpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<ST> >& thyraX,
                        Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX) {
  // build an Tpetra_MultiVector object
  LO numVectors = thyraX->domain()->dim();

  // make sure the number of vectors are the same
  TEUCHOS_ASSERT((size_t)numVectors == tpetraX.getNumVectors());
  TEUCHOS_ASSERT((Tpetra::global_size_t)thyraX->range()->dim() == tpetraX.getGlobalLength());

  // extract local information from the Tpetra_MultiVector
  LO leadingDim = 0, localDim = 0;
  leadingDim                       = tpetraX.getStride();
  Teuchos::ArrayRCP<ST> tpetraData = tpetraX.get1dViewNonConst();

  // perform recursive copy
  blockThyraToTpetra(numVectors, tpetraData, leadingDim, thyraX, localDim);

  // sanity check
  TEUCHOS_ASSERT((size_t)localDim == tpetraX.getLocalLength());
}

void thyraVSToTpetraMap(std::vector<GO>& myIndicies, int blockOffset,
                        const Thyra::VectorSpaceBase<ST>& vs, int& localDim) {
  // zero out set local dimension
  localDim = 0;

  const RCP<const Thyra::ProductVectorSpaceBase<ST> > prodVS =
      rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<ST> >(rcpFromRef(vs));

  // is more recursion needed?
  if (prodVS == Teuchos::null) {
    // base case

    // try to cast to an SPMD capable vector space
    const RCP<const Thyra::SpmdVectorSpaceBase<ST> > spmdVS =
        rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<ST> >(rcpFromRef(vs));
    TEUCHOS_TEST_FOR_EXCEPTION(spmdVS == Teuchos::null, std::runtime_error,
                               "thyraVSToTpetraMap requires all subblocks to be SPMD");

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
    thyraVSToTpetraMap(myIndicies, blockOffset, *prodVS->getBlock(blkIndex), subDim);

    blockOffset += prodVS->getBlock(blkIndex)->dim();

    // account for the size of this subblock
    localDim += subDim;
  }
}

// From a Thyra vector space create a compatable Tpetra_Map
const RCP<Tpetra::Map<LO, GO, NT> > thyraVSToTpetraMap(
    const Thyra::VectorSpaceBase<ST>& vs,
    const RCP<const Teuchos::Comm<Thyra::Ordinal> >& /* comm */) {
  int localDim = 0;
  std::vector<GO> myGIDs;

  // call recursive routine that constructs the mapping
  thyraVSToTpetraMap(myGIDs, 0, vs, localDim);

  TEUCHOS_ASSERT(myGIDs.size() == (size_t)localDim);

  // FIXME (mfh 12 Jul 2018) This ignores the input comm, so it can't
  // be right.

  // create the map
  return rcp(new Tpetra::Map<LO, GO, NT>(vs.dim(), Teuchos::ArrayView<const GO>(myGIDs), 0,
                                         Tpetra::getDefaultComm()));
}

}  // namespace TpetraHelpers
}  // end namespace Teko

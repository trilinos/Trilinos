/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/HexBeamCR.hpp>

HexBeamCR::HexBeamCR(int W, int D, int DofPerNode,
		 int decomp, int numProcs, int localProc)
  : HexBeam(W, D, DofPerNode, decomp, numProcs, localProc)
{
  totalNumElems_ = W*W*D;
  totalNumNodes_ = (W+1)*(W+1)*(D+1) + (W+1)*(W+1)*(2*numProcs_-1);

  numGlobalDOF_ = totalNumNodes_*dofPerNode_;

  numLocalSlices_ = D/numProcs;
  int remainder = D%numProcs;
  firstLocalSlice_ = localProc_*numLocalSlices_;

  switch(decomp) {
  case HexBeamCR::OneD:
    if (D < numProcs) {
      fei::console_out() << "HexBeamCR: too many processors." << FEI_ENDL;
      inErrorState_ = true;
      break;
    }
    if (localProc < remainder) {
      ++numLocalSlices_;
      ++firstLocalSlice_;
    }

    localNumNodes_ = numNodesPerSlice_*(numLocalSlices_+2);

    localNumElems_ = numElemsPerSlice_*numLocalSlices_;
    numLocalDOF_ = localNumNodes_*dofPerNode_;

    if (localProc > 0) {
      firstLocalElem_ = localProc*numLocalSlices_*numElemsPerSlice_;

      firstLocalNode_ = localProc*localNumNodes_;

      if (remainder <= localProc && remainder > 0) {
	firstLocalElem_ += remainder*numElemsPerSlice_;
	firstLocalNode_ += remainder*numNodesPerSlice_;
      }
    }

    break;

  case HexBeamCR::TwoD:
  case HexBeamCR::ThreeD:
  default:
    fei::console_out() << "HexBeamCR: invalid decomp option: " << decomp
	 <<" aborting." << FEI_ENDL;
    std::abort();
  }

  localCRslice_ = firstLocalSlice_ + numLocalSlices_/2;
  numLocalCRs_ = numNodesPerSlice_;
  if (localProc_ < numProcs_-1) {
    numLocalCRs_ += numNodesPerSlice_;
  }

  numNodesPerCR_ = 2;
}

HexBeamCR::~HexBeamCR()
{
}

int HexBeamCR::getElemConnectivity(int elemID, int* nodeIDs)
{
  if (elemID < firstLocalElem_ || elemID > firstLocalElem_+localNumElems_) {
    return(-1);
  }

  int whichGlobalSlice = elemID/numElemsPerSlice_;
  int elemX = elemID%W_;
  int elemY = (elemID%(W_*W_))/W_;
  //FEI_COUT << "whichGlobalSlice: " << whichGlobalSlice << FEI_ENDL;
  int firstElemNode = (whichGlobalSlice + localProc_*2)*numNodesPerSlice_
                     + elemY*(W_+1) + elemX;

  if (whichGlobalSlice >= localCRslice_) {
    firstElemNode += numNodesPerSlice_;
  }

  nodeIDs[0] = firstElemNode;
  nodeIDs[1] = firstElemNode+1;
  nodeIDs[2] = firstElemNode+W_+1;
  nodeIDs[3] = nodeIDs[2]+1;

  nodeIDs[4] = nodeIDs[0]+numNodesPerSlice_;
  nodeIDs[5] = nodeIDs[1]+numNodesPerSlice_;
  nodeIDs[6] = nodeIDs[2]+numNodesPerSlice_;
  nodeIDs[7] = nodeIDs[3]+numNodesPerSlice_;

  return(0);
}

int HexBeamCR::getCRNodes(int** nodeIDs)
{
  int offset = 0;
  int firstCRnode = (localCRslice_+localProc_*2)*numNodesPerSlice_;
  int i;
  for(i=0; i<numNodesPerSlice_; ++i) {
    nodeIDs[offset][0] = firstCRnode+i;
    nodeIDs[offset++][1] = firstCRnode+i+numNodesPerSlice_;
  }

  if (localProc_ >= numProcs_-1) return(0);

  int nextCRnode = firstLocalNode_ + localNumNodes_ - numNodesPerSlice_;
  for(i=0; i<numNodesPerSlice_; ++i) {
    nodeIDs[offset][0] = nextCRnode+i;
    nodeIDs[offset++][1] = nextCRnode+i+numNodesPerSlice_;    
  }

  return(0);
}

int HexBeamCR::getElemStiffnessMatrix(int elemID, double* elemMat)
{
  if (elemID < firstLocalElem_ || elemID > firstLocalElem_+localNumElems_) {
    return(-1);
  }

  int i, len = nodesPerElem_*dofPerNode_*nodesPerElem_*dofPerNode_;

  for(i=0; i<len; ++i) {
    elemMat[i] = 0.0;
  }

  //Should set up some semi-realistic stiffness-matrix coefficients here...
  //For now just use arbitrary numbers and set it up so the matrix won't be
  //too ill-conditioned. (This is intended for an assembly test more than
  //a solver test.)

  //Now set the diagonal to 4.0
  len = nodesPerElem_*dofPerNode_;
  for(i=0; i<len; ++i) {
    int offset = i*len+i;
    elemMat[offset] = 4.0;
  }

  //Now set some off-diagonals
  for(i=0; i<len; ++i) {
    int offset = i*len+i;
    if (i>1) {
      elemMat[offset-2] = -0.5;
    }

    if (i<len-2) {
      elemMat[offset+2] = -0.5;
    }

    if (i>3) {
      elemMat[offset-4] = -0.1;
    }
    if (i<len-4) {
      elemMat[offset+4] = -0.1;
    }
  }

  return(0);
}

int HexBeamCR::getElemLoadVector(int elemID, double* elemVec)
{
  if (elemID < firstLocalElem_ || elemID > firstLocalElem_+localNumElems_) {
    return(-1);
  }

  int i, len = nodesPerElem_*dofPerNode_;
  for(i=0; i<len; ++i) {
    elemVec[i] = 1.0;
  }

  return(0);
}

int HexBeamCR::getNumBCNodes()
{
  int numBCNodes = (numLocalSlices_+1)*(W_+1);
  return( numBCNodes );
}

int HexBeamCR::getBCNodes(int numNodes, int* nodeIDs)
{
  if (numNodes != getNumBCNodes()) {
    return(-1);
  }

  int firstBCNode = firstLocalNode_ + W_;

  for(int i=0; i<numNodes; ++i) {
    nodeIDs[i] = firstBCNode + W_+1;
  }

  return(0);
}

int HexBeamCR::getBCGammaValues(int numBCDofs, double* gamma)
{
  if (numBCDofs != getNumBCNodes()*dofPerNode_) {
    return(-1);
  }

  for(int i=0; i<numBCDofs; ++i) {
    gamma[i] = 2.0;
  }

  return(0);
}

int HexBeamCR::getNumSharedNodes()
{
  if (numProcs_ < 2) return(0);

  int numSharedNodes = numNodesPerSlice_;
  if (localProc_ > 0 && localProc_ < numProcs_-1) {
    numSharedNodes += numNodesPerSlice_;
  }

  return(numSharedNodes);
}

int HexBeamCR::getSharedNodes(int numSharedNodes,
			    int*& sharedNodes,
			    int*& numSharingProcsPerNode,
			    int**& sharingProcs)
{
  if (numProcs_ < 2) return(0);

  if (numSharedNodes != getNumSharedNodes()) {
    return(-1);
  }

  sharedNodes = new int[numSharedNodes];
  numSharingProcsPerNode = new int[numSharedNodes];
  sharingProcs = new int*[numSharedNodes];
  int* sharingProcVals = new int[numSharedNodes];
  if (sharedNodes == NULL || numSharingProcsPerNode == NULL ||
      sharingProcs == NULL || sharingProcVals == NULL) {
    return(-1);
  }

  int i;
  for(i=0; i<numSharedNodes; ++i) {
    numSharingProcsPerNode[i] = 1;
    sharingProcs[i] = &(sharingProcVals[i]);
  }

  int firstSharedNode = firstLocalNode_+numNodesPerSlice_*(numLocalSlices_+2);
  int offset = 0;
  //FEI_COUT << localProc_ << ": firstSharedNode: " << firstSharedNode << FEI_ENDL;
  if (localProc_ < numProcs_ - 1) {
    for(i=0; i<numNodesPerSlice_; ++i) {
      sharedNodes[offset] = firstSharedNode+i;
      sharingProcs[offset++][0] = localProc_+1;
    }
  }

  firstSharedNode = firstLocalNode_;
  //FEI_COUT << localProc_ << ":+1 firstSharedNode: " << firstSharedNode << FEI_ENDL;
  if (localProc_ > 0) {
    for(i=0; i<numNodesPerSlice_; ++i) {
      sharedNodes[offset] = firstSharedNode+i;
      sharingProcs[offset++][0] = localProc_-1;
    }
  }

  return(0);
}

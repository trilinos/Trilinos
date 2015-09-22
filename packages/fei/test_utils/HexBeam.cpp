/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <test_utils/HexBeam.hpp>

#undef fei_file
#define fei_file "HexBeam.cpp"

HexBeam::HexBeam(int W, int D, int DofPerNode,
		 int decomp, int numProcs, int localProc)
  : W_(W),
    D_(D),
    decomp_(decomp),
    numProcs_(numProcs),
    localProc_(localProc),
    firstLocalElem_(0),
    firstLocalNode_(0),
    inErrorState_(false),
    nodesPerElem_(8),
    dofPerNode_(DofPerNode)
{
  totalNumElems_ = W*W*D;
  totalNumNodes_ = (W+1)*(W+1)*(D+1);
  numElemsPerSlice_ = W*W;
  numNodesPerSlice_ = (W+1)*(W+1);

  numGlobalDOF_ = totalNumNodes_*dofPerNode_;

  numLocalSlices_ = D/numProcs;
  int remainder = D%numProcs;

  switch(decomp) {
  case HexBeam::OneD:
    if (D < numProcs) {
      throw std::runtime_error("HexBeam: size D must be greater or equal num-procs.");
    }
    if (localProc < remainder) {
      ++numLocalSlices_;
    }

    localNumNodes_ = numNodesPerSlice_*(numLocalSlices_+1);
    localNumElems_ = numElemsPerSlice_*numLocalSlices_;
    numLocalDOF_ = localNumNodes_*dofPerNode_;

    if (localProc > 0) {
      firstLocalElem_ = localProc*numLocalSlices_*numElemsPerSlice_;
      firstLocalNode_ = localProc*numLocalSlices_*numNodesPerSlice_;
      if (remainder <= localProc && remainder > 0) {
	firstLocalElem_ += remainder*numElemsPerSlice_;
	firstLocalNode_ += remainder*numNodesPerSlice_;
      }
    }

    break;

  case HexBeam::TwoD:
  case HexBeam::ThreeD:
  default:
    fei::console_out() << "HexBeam: invalid decomp option: " << decomp
	 <<" aborting." << FEI_ENDL;
    std::abort();
  }
}

HexBeam::~HexBeam()
{
}

int HexBeam::getElemConnectivity(int elemID, int* nodeIDs)
{
  if (elemID < firstLocalElem_ || elemID > firstLocalElem_+localNumElems_) {
    return(-1);
  }

  int whichGlobalSlice = elemID/numElemsPerSlice_;
  int elemX = elemID%W_;
  int elemY = (elemID%(W_*W_))/W_;

  int firstElemNode = whichGlobalSlice*numNodesPerSlice_
                     + elemY*(W_+1) + elemX;

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

int HexBeam::getElemStiffnessMatrix(int elemID, double* elemMat)
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

int HexBeam::getElemLoadVector(int elemID, double* elemVec)
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

int HexBeam::getNumBCNodes()
{
  int numBCNodes = (numLocalSlices_+1)*(W_+1);
  return( numBCNodes );
}

int HexBeam::getBCNodes(int numNodes, int* nodeIDs)
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

int HexBeam::getBCValues(int numBCNodes, int* offsetsIntoField, double* vals)
{
  if (numBCNodes != getNumBCNodes()) {
    return(-1);
  }

  for(int i=0; i<numBCNodes; ++i) {
    offsetsIntoField[i] = 0;
    vals[i] = 2.0;
  }

  return(0);
}

int HexBeam::getNumSharedNodes()
{
  if (numProcs_ < 2) return(0);

  int numSharedNodes = numNodesPerSlice_;
  if (localProc_ > 0 && localProc_ < numProcs_-1) {
    numSharedNodes += numNodesPerSlice_;
  }

  return(numSharedNodes);
}

int HexBeam::getSharedNodes(int numSharedNodes,
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

  int firstSharedNode = firstLocalNode_+numNodesPerSlice_*numLocalSlices_;
  int offset = 0;

  if (localProc_ < numProcs_ - 1) {
    for(i=0; i<numNodesPerSlice_; ++i) {
      sharedNodes[offset] = firstSharedNode+i;
      sharingProcs[offset++][0] = localProc_+1;
    }
  }

  firstSharedNode = firstLocalNode_;
  if (localProc_ > 0) {
    for(i=0; i<numNodesPerSlice_; ++i) {
      sharedNodes[offset] = firstSharedNode+i;
      sharingProcs[offset++][0] = localProc_-1;
    }
  }

  return(0);
}

namespace HexBeam_Functions {

int print_cube_data(HexBeam& hexcube, int numProcs, int localProc)
{
  FEI_COUT << localProc << ": num elems: " << hexcube.numLocalElems() << FEI_ENDL;
  int i;
  int* nodeIDs = new int[hexcube.numNodesPerElem()];
  int firstLocalElem = hexcube.firstLocalElem();

  for(i=0; i<hexcube.numLocalElems(); ++i) {
    hexcube.getElemConnectivity(firstLocalElem+i, nodeIDs);
    FEI_COUT << localProc << ": elem " << firstLocalElem+i << ", nodes: ";
    for(int j=0; j<hexcube.numNodesPerElem(); ++j) {
      FEI_COUT << nodeIDs[j] << " ";
    }
    FEI_COUT << FEI_ENDL;
  }

  delete [] nodeIDs;

  return(0);
}


int init_elem_connectivities(FEI* fei, HexBeam& hexcube)
{
  int numLocalElems = hexcube.numLocalElems();
  int firstLocalElem = hexcube.firstLocalElem();
  int nodesPerElem = hexcube.numNodesPerElem();
  int fieldID = 0;

  int** fieldIDsTable = new int*[nodesPerElem];
  int* numFieldsPerNode = new int[nodesPerElem];

  for(int j=0; j<nodesPerElem; ++j) {
    numFieldsPerNode[j] = 1;
    fieldIDsTable[j] = new int[numFieldsPerNode[j]];
    for(int k=0; k<numFieldsPerNode[j]; ++k) {
      fieldIDsTable[j][k] = fieldID;
    }
  }

  int blockID = 0;
  CHK_ERR( fei->initElemBlock(blockID,
			      numLocalElems,
			      nodesPerElem,
			      numFieldsPerNode,
			      fieldIDsTable,
			      0, // no element-centered degrees-of-freedom
			      NULL, //null list of elem-dof fieldIDs
			      FEI_NODE_MAJOR) );


  int* nodeIDs = new int[nodesPerElem];
  if (nodeIDs == NULL) return(-1);

  for(int i=0; i<numLocalElems; ++i) {
    CHK_ERR( hexcube.getElemConnectivity(firstLocalElem+i, nodeIDs) );

    CHK_ERR( fei->initElem(blockID, firstLocalElem+i, nodeIDs) );
  }

  delete [] nodeIDs;
  delete [] numFieldsPerNode;
  for(int jj=0; jj<nodesPerElem; ++jj) {
    delete [] fieldIDsTable[jj];
  }
  delete [] fieldIDsTable;

  return(0);
}

int init_shared_nodes(FEI* fei, HexBeam& hexcube)
{
  int numSharedNodes = hexcube.getNumSharedNodes();
  if (numSharedNodes == 0) {
    return(0);
  }

  int* sharedNodes = NULL;
  int* numSharingProcsPerNode = NULL;
  int** sharingProcs = NULL;
  if (numSharedNodes > 0) {
    CHK_ERR( hexcube.getSharedNodes(numSharedNodes,
				    sharedNodes, numSharingProcsPerNode,
				    sharingProcs) );
  }

  CHK_ERR( fei->initSharedNodes(numSharedNodes, sharedNodes,
			      numSharingProcsPerNode, sharingProcs) );

  delete [] sharedNodes;
  delete [] numSharingProcsPerNode;
  delete [] sharingProcs[0];
  delete [] sharingProcs;

  return(0);
}

int init_constraints(FEI* fei, HexBeam& hexcube, int& firstLocalCRID)
{
  int numCRs = hexcube.getNumCRs();
  if (numCRs < 1) {
    return(0);
  }

  int numNodesPerCR = hexcube.getNumNodesPerCR();
  int* crnodes_1d = new int[numCRs*numNodesPerCR];
  int** crNodes = new int*[numCRs];
  int i, offset = 0;
  for(i=0; i<numCRs; ++i) {
    crNodes[i] = &(crnodes_1d[offset]);
    offset += numNodesPerCR;
  }

  CHK_ERR( hexcube.getCRNodes(crNodes) );

  int crID;
  int* fieldIDs = new int[numNodesPerCR];
  for(i=0; i<numNodesPerCR; ++i) fieldIDs[i] = 0;

  for(i=0; i<numCRs; ++i) {
    CHK_ERR( fei->initCRMult(numNodesPerCR, crNodes[i], fieldIDs, crID) );
//     FEI_COUT << "crID: " << crID << ", nodes: ";
//     for(int j=0; j<numNodesPerCR; ++j) {
//       FEI_COUT << crNodes[i][j] << " ";
//     }
//     FEI_COUT << FEI_ENDL;

    if (i == 0) {
      firstLocalCRID = crID;
    }
  }

  delete [] crnodes_1d;
  delete [] crNodes;
  delete [] fieldIDs;

  return(0);
}

int load_constraints(FEI* fei, HexBeam& hexcube, int firstLocalCRID)
{
  int numCRs = hexcube.getNumCRs();
  if (numCRs < 1) {
    return(0);
  }

  int numNodesPerCR = hexcube.getNumNodesPerCR();
  int* crnodes_1d = new int[numCRs*numNodesPerCR];
  int** crNodes = new int*[numCRs];
  int i, offset = 0;
  for(i=0; i<numCRs; ++i) {
    crNodes[i] = &(crnodes_1d[offset]);
    offset += numNodesPerCR;
  }

  CHK_ERR( hexcube.getCRNodes(crNodes) );

  int* fieldIDs = new int[numNodesPerCR];
  for(i=0; i<numNodesPerCR; ++i) fieldIDs[i] = 0;

  int fieldSize = hexcube.numDofPerNode();
  double* weights = new double[fieldSize*numNodesPerCR];

  for(i=0; i<fieldSize*numNodesPerCR; ++i) weights[i] = 0.0;
  weights[0] = -1.0;
  weights[fieldSize] = 1.0;
  double rhsValue = 0.0;

  for(i=0; i<numCRs; ++i) {
    CHK_ERR( fei->loadCRMult(firstLocalCRID+i,
			     numNodesPerCR, crNodes[i], fieldIDs,
			     weights, rhsValue) );
  }

  delete [] crnodes_1d;
  delete [] crNodes;
  delete [] fieldIDs;
  delete [] weights;

  return(0);
}

int load_elem_data(FEI* fei, HexBeam& hexcube)
{
  int blockID = 0;
  int numLocalElems = hexcube.localNumElems_;
  int firstLocalElem = hexcube.firstLocalElem_;
  int nodesPerElem = hexcube.numNodesPerElem();
  int fieldSize = hexcube.numDofPerNode();

  int len = nodesPerElem*fieldSize;
  double* elemMat = new double[len*len];
  double** elemMat2D = new double*[len];
  double* elemVec = new double[len];

  int* nodeIDs = new int[nodesPerElem];

  if (elemMat == NULL || elemMat2D == NULL || elemVec == NULL ||
      nodeIDs == NULL) {
    return(-1);
  }

  for(int j=0; j<len; ++j) {
    elemMat2D[j] = &(elemMat[j*len]);
  }

  CHK_ERR( hexcube.getElemStiffnessMatrix(firstLocalElem, elemMat) );
  CHK_ERR( hexcube.getElemLoadVector(firstLocalElem, elemVec) );

  for(int i=0; i<numLocalElems; ++i) {
    CHK_ERR( hexcube.getElemConnectivity(firstLocalElem+i, nodeIDs) );

    CHK_ERR( fei->sumInElemMatrix(blockID, firstLocalElem+i,
				  nodeIDs, elemMat2D, FEI_DENSE_ROW) );

    CHK_ERR( fei->sumInElemRHS(blockID, firstLocalElem+i, nodeIDs, elemVec) );
  }

  delete [] elemMat;
  delete [] elemMat2D;
  delete [] elemVec;
  delete [] nodeIDs;

  return(0);
}

int load_BC_data(FEI* fei, HexBeam& hexcube)
{
  int numBCNodes = hexcube.getNumBCNodes();
  if (numBCNodes == 0) {
    return(0);
  }

  int* nodeIDs = new int[numBCNodes];

  int fieldID = 0;

  int* offsetsIntoField = new int[numBCNodes];
  double* prescribed_vals = new double[numBCNodes];

  CHK_ERR( hexcube.getBCNodes(numBCNodes, nodeIDs) );

  CHK_ERR( hexcube.getBCValues(numBCNodes, offsetsIntoField, prescribed_vals) );

  CHK_ERR( fei->loadNodeBCs(numBCNodes, nodeIDs,
			    fieldID, offsetsIntoField, prescribed_vals) );

  delete [] nodeIDs;
  delete [] offsetsIntoField;
  delete [] prescribed_vals;

  return(0);
}

int init_elem_connectivities(fei::MatrixGraph* matrixGraph, HexBeam& hexcube)
{
  int numLocalElems = hexcube.localNumElems_;
  int firstLocalElem = hexcube.firstLocalElem_;
  int nodesPerElem = hexcube.numNodesPerElem();
  int fieldID = 0;
//  int fieldSize = hexcube.numDofPerNode();
  int nodeIDType = 0;


  int patternID = 0;
//  if (fieldSize > 1) {
    patternID = matrixGraph->definePattern(nodesPerElem,
			       nodeIDType, fieldID);
//  }
//  else {
//    //if fieldSize == 1, let's not even bother associating a field with
//    //our mesh-nodes. fei:: objects assume that identifiers without an
//    //associated field always have exactly one degree-of-freedom.
//    //
//    patternID = matrixGraph->definePattern(nodesPerElem, nodeIDType);
//  }

  int blockID = 0;
  CHK_ERR( matrixGraph->initConnectivityBlock(blockID, numLocalElems, patternID));

  int* nodeIDs = new int[nodesPerElem];
  if (nodeIDs == NULL) return(-1);

  for(int i=0; i<numLocalElems; ++i) {
    CHK_ERR( hexcube.getElemConnectivity(firstLocalElem+i, nodeIDs) );

    CHK_ERR( matrixGraph->initConnectivity(blockID, firstLocalElem+i, nodeIDs) );
  }

  delete [] nodeIDs;

  return(0);
}

int init_shared_nodes(fei::MatrixGraph* matrixGraph, HexBeam& hexcube)
{
  int numSharedNodes = hexcube.getNumSharedNodes();
  if (numSharedNodes == 0) {
    return(0);
  }

  int* sharedNodes = NULL;
  int* numSharingProcsPerNode = NULL;
  int** sharingProcs = NULL;
  if (numSharedNodes > 0) {
    CHK_ERR( hexcube.getSharedNodes(numSharedNodes,
				    sharedNodes, numSharingProcsPerNode,
				    sharingProcs) );
  }

  fei::SharedPtr<fei::VectorSpace> nodeSpace = matrixGraph->getRowSpace();

  int nodeIDType = 0;

  CHK_ERR( nodeSpace->initSharedIDs(numSharedNodes, nodeIDType,
				    sharedNodes, numSharingProcsPerNode,
				    sharingProcs) );

  delete [] sharedNodes;
  delete [] numSharingProcsPerNode;
  delete [] sharingProcs[0];
  delete [] sharingProcs;

  return(0);
}

int init_constraints(fei::MatrixGraph* matrixGraph, HexBeam& hexcube,
		     int localProc, int& firstLocalCRID)
{
  int numCRs = hexcube.getNumCRs();
  if (numCRs < 1) {
    return(0);
  }

  int numNodesPerCR = hexcube.getNumNodesPerCR();
  int* crnodes_1d = new int[numCRs*numNodesPerCR];
  int** crNodes = new int*[numCRs];
  int i, offset = 0;
  for(i=0; i<numCRs; ++i) {
    crNodes[i] = &(crnodes_1d[offset]);
    offset += numNodesPerCR;
  }

  CHK_ERR( hexcube.getCRNodes(crNodes) );

  int crID = localProc*100000;
  firstLocalCRID = crID;

  int nodeIDType = 0;
  int crIDType = 1;

  int* fieldIDs = new int[numNodesPerCR];
  int* idTypes = new int[numNodesPerCR];
  for(i=0; i<numNodesPerCR; ++i) {
    fieldIDs[i] = 0;
    idTypes[i] = nodeIDType;
  }

  for(i=0; i<numCRs; ++i) {
    CHK_ERR( matrixGraph->initLagrangeConstraint(crID+i, crIDType,
						 numNodesPerCR,
						 idTypes, crNodes[i],
						 fieldIDs) );
//     FEI_COUT << "crID: " << crID << ", nodes: ";
//     for(int j=0; j<numNodesPerCR; ++j) {
//       FEI_COUT << crNodes[i][j] << " ";
//     }
//     FEI_COUT << FEI_ENDL;

  }

  delete [] crnodes_1d;
  delete [] crNodes;
  delete [] fieldIDs;

  return(0);
}

int init_slave_constraints(fei::MatrixGraph* matrixGraph, HexBeam& hexcube)
{
  int numCRs = hexcube.getNumCRs();
  if (numCRs < 1) {
    return(0);
  }

  int numNodesPerCR = hexcube.getNumNodesPerCR();
  int* crnodes_1d = new int[numCRs*numNodesPerCR];
  int** crNodes = new int*[numCRs];
  int i, offset = 0;
  for(i=0; i<numCRs; ++i) {
    crNodes[i] = &(crnodes_1d[offset]);
    offset += numNodesPerCR;
  }

  CHK_ERR( hexcube.getCRNodes(crNodes) );

  int nodeIDType = 0;

  int* fieldIDs = new int[numNodesPerCR];
  int* idTypes = new int[numNodesPerCR];
  for(i=0; i<numNodesPerCR; ++i) {
    fieldIDs[i] = 0;
    idTypes[i] = nodeIDType;
  }

  int fieldSize = hexcube.numDofPerNode();
  double* weights = new double[fieldSize*numNodesPerCR];

  for(i=0; i<fieldSize*numNodesPerCR; ++i) weights[i] = 0.0;
  weights[0] = -1.0;
  weights[fieldSize] = 1.0;
  double rhsValue = 0.0;
  int offsetOfSlave = 0;
  int offsetIntoSlaveField = 0;

  for(i=0; i<numCRs; ++i) {
    CHK_ERR( matrixGraph->initSlaveConstraint(numNodesPerCR,
					      idTypes,
					      crNodes[i],
					      fieldIDs,
					      offsetOfSlave,
					      offsetIntoSlaveField,
					      weights,
					      rhsValue) );
  }

  delete [] crnodes_1d;
  delete [] crNodes;
  delete [] fieldIDs;
  delete [] weights;

  return(0);
}

int load_elem_data(fei::MatrixGraph* matrixGraph,
		   fei::Matrix* mat,
		   fei::Vector* rhs,
		   HexBeam& hexcube)
{
  int blockID = 0;
  int numLocalElems = hexcube.localNumElems_;
  int firstLocalElem = hexcube.firstLocalElem_;
  int nodesPerElem = hexcube.numNodesPerElem();
  int fieldSize = hexcube.numDofPerNode();

  int len = nodesPerElem*fieldSize;
  double* elemMat = new double[len*len];
  double** elemMat2D = new double*[len];
  double* elemVec = new double[len];

  if (elemMat == NULL || elemMat2D == NULL || elemVec == NULL) {
    return(-1);
  }

  for(int j=0; j<len; ++j) {
    elemMat2D[j] = &(elemMat[j*len]);
  }

  CHK_ERR( hexcube.getElemStiffnessMatrix(firstLocalElem, elemMat) );
  CHK_ERR( hexcube.getElemLoadVector(firstLocalElem, elemVec) );

  bool block_matrix = mat->usingBlockEntryStorage();
  std::vector<int> indices(len);

  if (block_matrix) {
    mat->getMatrixGraph()->setIndicesMode(fei::MatrixGraph::POINT_ENTRY_GRAPH);
  }

  for(int i=0; i<numLocalElems; ++i) {
    CHK_ERR( mat->getMatrixGraph()->getConnectivityIndices(blockID,
                                                           firstLocalElem+i,
                                                           len, &indices[0],
                                                           len) );
    CHK_ERR( mat->sumIn(len, &indices[0], len, &indices[0],
                        elemMat2D, FEI_DENSE_COL) );
    CHK_ERR( rhs->sumIn(len, &indices[0], elemVec, 0) );
  }

  if (block_matrix) {
    mat->getMatrixGraph()->setIndicesMode(fei::MatrixGraph::BLOCK_ENTRY_GRAPH);
  }
  delete [] elemMat;
  delete [] elemMat2D;
  delete [] elemVec;

  return(0);
}

int load_constraints(fei::LinearSystem* linSys, HexBeam& hexcube,
		     int firstLocalCRID)
{
  int numCRs = hexcube.getNumCRs();
  if (numCRs < 1) {
    return(0);
  }

  int numNodesPerCR = hexcube.getNumNodesPerCR();

  int fieldSize = hexcube.numDofPerNode();
  double* weights = new double[fieldSize*numNodesPerCR];

  int i;
  for(i=0; i<fieldSize*numNodesPerCR; ++i) weights[i] = 0.0;
  weights[0] = -1.0;
  weights[fieldSize] = 1.0;
  double rhsValue = 0.0;

  for(i=0; i<numCRs; ++i) {
    CHK_ERR( linSys->loadLagrangeConstraint(firstLocalCRID+i,
					    weights, rhsValue) );
  }

  delete [] weights;

  return(0);
}

int load_BC_data(fei::LinearSystem* linSys, HexBeam& hexcube)
{
  int numBCNodes = hexcube.getNumBCNodes();
  if (numBCNodes == 0) {
    return(0);
  }

  int* nodeIDs = new int[numBCNodes];

  int fieldID = 0;
  int nodeIDType = 0;

  int* offsetsIntoField = new int[numBCNodes];
  double* prescribed_vals = new double[numBCNodes];

  CHK_ERR( hexcube.getBCNodes(numBCNodes, nodeIDs) );

  CHK_ERR( hexcube.getBCValues(numBCNodes, offsetsIntoField, prescribed_vals) );

  CHK_ERR( linSys->loadEssentialBCs(numBCNodes, nodeIDs,
				    nodeIDType, fieldID,
				    offsetsIntoField, prescribed_vals) );

  delete [] offsetsIntoField;
  delete [] prescribed_vals;
  delete [] nodeIDs;

  return(0);
}

}//namespace HexBeam_Functions

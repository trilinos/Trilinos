/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include "fei_sstream.hpp"
#include "fei_fstream.hpp"

#include <limits>
#include <cmath>
#include <assert.h>

#include "fei_defs.h"

#include "fei_TemplateUtils.hpp"
#include <fei_CommUtils.hpp>
#include "snl_fei_Constraint.hpp"
typedef snl_fei::Constraint<GlobalID> ConstraintType;

#include "fei_NodeDescriptor.hpp"
#include "fei_NodeCommMgr.hpp"
#include "fei_NodeDatabase.hpp"

#include "fei_SlaveVariable.hpp"

#include "fei_BlockDescriptor.hpp"

#include "snl_fei_PointBlockMap.hpp"
#include "fei_ProcEqns.hpp"
#include "fei_EqnBuffer.hpp"
#include <fei_FillableMat.hpp>
#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>
#include "fei_EqnCommMgr.hpp"

#include "fei_Lookup.hpp"
#include "fei_ConnectivityTable.hpp"
#include "snl_fei_Utils.hpp"
#include "SNL_FEI_Structure.hpp"

#undef fei_file
#define fei_file "SNL_FEI_Structure.cpp"
#include "fei_ErrMacros.hpp"

//-----Constructor-------------------------------------------------------------
SNL_FEI_Structure::SNL_FEI_Structure(MPI_Comm comm)
 : comm_(comm),
   localProc_(0),
   masterProc_(0),
   numProcs_(1),
   fieldIDs_(),
   fieldSizes_(),
   fieldDatabase_(new std::map<int,int>),
   fieldDofMap_(),
   workarray_(),
   blockIDs_(),
   blocks_(),
   connTables_(),
   nodeDatabase_(NULL),
   activeNodesInitialized_(false),
   globalNodeOffsets_(),
   globalEqnOffsets_(),
   globalBlkEqnOffsets_(),
   slaveVars_(NULL),
   slaveEqns_(NULL),
   slvEqnNumbers_(NULL),
   numSlvs_(0),
   lowestSlv_(0),
   highestSlv_(0),
   slaveMatrix_(NULL),
   globalNumNodesVanished_(),
   localVanishedNodeNumbers_(),
   nodeCommMgr_(NULL),
   eqnCommMgr_(NULL),
   slvCommMgr_(NULL),
   numGlobalEqns_(0),
   numLocalEqns_(0),
   localStartRow_(0),
   localEndRow_(0),
   numLocalNodalEqns_(0),
   numLocalElemDOF_(0),
   numLocalMultCRs_(0),
   reducedStartRow_(0),
   reducedEndRow_(0),
   numLocalReducedRows_(0),
   Kid_(NULL),
   Kdi_(NULL),
   Kdd_(NULL),
   tmpMat1_(),
   tmpMat2_(),
   reducedEqnCounter_(0),
   reducedRHSCounter_(0),
   rSlave_(),
   cSlave_(),
   work_nodePtrs_(),
   structureFinalized_(false),
   generateGraph_(true),
   sysMatIndices_(NULL),
   blockMatrix_(false),
   numGlobalEqnBlks_(0),
   numLocalEqnBlks_(0),
   numLocalReducedEqnBlks_(0),
   localBlkOffset_(0),
   localReducedBlkOffset_(0),
   globalMaxBlkSize_(0),
   firstLocalNodeNumber_(-1),
   numGlobalNodes_(0),
   sysBlkMatIndices_(NULL),
   matIndicesDestroyed_(false),
   workSpace_(),
   blkEqnMapper_(new snl_fei::PointBlockMap()),
   multCRs_(),
   penCRs_(),
   checkSharedNodes_(false),
   name_(),
   outputLevel_(0),
   debugOutput_(false),
   dbgPath_(),
   dbgOStreamPtr_(NULL),
   setDbgOutCalled_(true)
{
  numProcs_ = 1, localProc_ = 0, masterProc_ = 0;

  localProc_ = fei::localProc(comm_);
  numProcs_ = fei::numProcs(comm_);
  masterProc_ = 0;

  slaveVars_ = new std::vector<SlaveVariable*>;
  slaveEqns_ = new EqnBuffer;

  nodeCommMgr_ = new NodeCommMgr(comm_, *this);

  eqnCommMgr_ = new EqnCommMgr(comm_);
  eqnCommMgr_->setNumRHSs(1);

  nodeDatabase_ = new NodeDatabase(fieldDatabase_, nodeCommMgr_);

  Kid_ = new fei::FillableMat;
  Kdi_ = new fei::FillableMat;
  Kdd_ = new fei::FillableMat;
}

//-----------------------------------------------------------------------------
int SNL_FEI_Structure::parameters(int numParams, const char*const* paramStrings)
{
  const char* param = NULL;

  param = snl_fei::getParamValue("outputLevel",numParams,paramStrings);
  if (param != NULL){
    std::string str(param);
    FEI_ISTRINGSTREAM isstr(str);
    isstr >> outputLevel_;
  }

  param = snl_fei::getParam("debugOutput",numParams,paramStrings);
  if (param != NULL){
    debugOutput_ = true;
  }

  param = snl_fei::getParam("debugOutputOff",numParams,paramStrings);
  if (param != NULL){
    debugOutput_ = false;
  }

  param = snl_fei::getParam("FEI_CHECK_SHARED_IDS", numParams,paramStrings);
  if (param != NULL){
    checkSharedNodes_ = true;
  }

  param = snl_fei::getParamValue("sharedNodeOwnership",
					numParams,paramStrings);
  if (param != NULL){
    if (!strcmp(param, "LowNumberedProc")) {
      nodeCommMgr_->setSharedOwnershipRule(NodeCommMgr::STRICTLY_LOW_PROC);
    }
    if (!strcmp(param, "ProcWithLocalElem")) {
      nodeCommMgr_->setSharedOwnershipRule(NodeCommMgr::PROC_WITH_LOCAL_ELEM);
    }
  }

  return(0);
}

//-----Destructor--------------------------------------------------------------
SNL_FEI_Structure::~SNL_FEI_Structure()
{
  for(size_t j=0; j<slaveVars_->size(); j++) {
    delete (*slaveVars_)[j];
  }
  delete slaveVars_;

  delete slaveEqns_;
  delete slaveMatrix_;

  destroyBlockRoster();
  destroyConnectivityTables();

  delete nodeCommMgr_;
  delete eqnCommMgr_;
  delete blkEqnMapper_;

  destroyMatIndices();

  deleteMultCRs();

  int numPCRs = penCRs_.size();
  if (numPCRs > 0) {
    fei::destroyValues(penCRs_);
    penCRs_.clear();
  }

  delete nodeDatabase_;
  delete fieldDatabase_;

  delete Kid_;
  delete Kdi_;
  delete Kdd_;
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::setDbgOut(FEI_OSTREAM& ostr,
				 const char* path, const char* feiName)
{
  dbgOStreamPtr_ = &ostr;
  setDbgOutCalled_ = true;
  if (path != NULL) {
    dbgPath_ = path;
  }
  else {
    dbgPath_ = ".";
  }
  if (feiName != NULL) {
    name_ = feiName;
  }

  debugOutput_ = true;
  return(0);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::destroyBlockRoster()
{
  for(size_t i=0; i<blockIDs_.size(); i++) delete blocks_[i];
  blocks_.resize(0);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::destroyConnectivityTables()
{
  for(size_t i=0; i<blockIDs_.size(); i++) {
    delete connTables_[i];
  }
  connTables_.resize(0);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::destroyMatIndices()
{
  if (matIndicesDestroyed_ == true) return;

  delete [] sysMatIndices_;
  sysMatIndices_ = NULL;

  delete [] sysBlkMatIndices_;
  sysBlkMatIndices_ = NULL;

  matIndicesDestroyed_ = true;
}

//------------------------------------------------------------------------------
const int* SNL_FEI_Structure::getNumFieldsPerNode(GlobalID blockID)
{
  BlockDescriptor* block = NULL;
  int err = getBlockDescriptor(blockID, block);
  if (err) return(NULL);

  return(block->fieldsPerNodePtr());
}

//------------------------------------------------------------------------------
const int* const* SNL_FEI_Structure::getFieldIDsTable(GlobalID blockID)
{
  BlockDescriptor* block = NULL;
  int err = getBlockDescriptor(blockID, block);
  if (err) return(NULL);

  return(block->fieldIDsTablePtr());
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getElemBlockInfo(GlobalID blockID,
                         int& interleaveStrategy, int& lumpingStrategy,
                         int& numElemDOF, int& numElements,
                         int& numNodesPerElem, int& numEqnsPerElem)
{
  BlockDescriptor* block = NULL;
  int err = getBlockDescriptor(blockID, block);
  if (err) {
      interleaveStrategy = -1;
      lumpingStrategy = -1;
      numElemDOF = -1;
      numElements = -1;
      numNodesPerElem = -1;
      numEqnsPerElem = -1;
   }

  interleaveStrategy = block->getInterleaveStrategy();
  lumpingStrategy = block->getLumpingStrategy();
  numElemDOF = block->getNumElemDOFPerElement();
  numElements = block->getNumElements();
  numNodesPerElem = block->getNumNodesPerElement();
  numEqnsPerElem = block->getNumEqnsPerElement();
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getEqnNumber(int nodeNumber, int fieldID)
{
  //this function is only used by clients of the Lookup interface, who are
  //expecting equations in the 'reduced' equation space.

  int eqnNumber;

  const NodeDescriptor* node = NULL;
  CHK_ERR( nodeDatabase_->getNodeWithNumber(nodeNumber, node) );

  bool hasField = node->getFieldEqnNumber(fieldID, eqnNumber);
  if (!hasField) {
    fei::console_out() << "SNL_FEI_Structure::getEqnNumber: ERROR, node with nodeNumber "
	 << nodeNumber << " doesn't have fieldID " << fieldID << FEI_ENDL;
    ERReturn(-1);
  }

  int reducedEqn = -1;
  bool isSlave = translateToReducedEqn(eqnNumber, reducedEqn);
  if (isSlave) return(-1);

  return(reducedEqn);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getOwnerProcForEqn(int eqn)
{
  if (eqn < 0) return(-1);

  int len = globalEqnOffsets_.size(); // len is numProcs+1...

  for(int i=0; i<len-1; i++) {
    if (eqn >= globalEqnOffsets_[i] && eqn < globalEqnOffsets_[i+1]) return(i);
  }

  return(-1);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initFields(int numFields, 
				 const int *fieldSizes, 
				 const int *fieldIDs,
         const int *fieldTypes)
{
  // store the incoming solution fields
  //
  if (debugOutput_) {
    FEI_OSTREAM& ostr = dbgOut();
    ostr << "FEI: initFields" << FEI_ENDL
	 << "#num-fields" << FEI_ENDL << numFields << FEI_ENDL;
    int nf;
    ostr << "#field-sizes" << FEI_ENDL;
    for(nf=0; nf<numFields; nf++) {
      ostr <<fieldSizes[nf] << " ";
    }
    ostr << FEI_ENDL << "#field-ids" << FEI_ENDL;
    for(nf=0; nf<numFields; nf++) {
      ostr << fieldIDs[nf] << " ";
    }
    if (fieldTypes != NULL) {
      ostr << FEI_ENDL << "#field-types" << FEI_ENDL;
      for(nf=0; nf<numFields; nf++) {
        ostr << fieldTypes[nf] << " ";
      }
    }
    ostr<<FEI_ENDL;
  }

  for (int i=0; i<numFields; i++) {
    fieldDatabase_->insert(std::pair<int,int>(fieldIDs[i], fieldSizes[i]));
    int offs = fei::sortedListInsert(fieldIDs[i], fieldIDs_); 
    if (offs >= 0) fieldSizes_.insert(fieldSizes_.begin()+offs, fieldSizes[i]);

    if (fieldIDs[i] >= 0) {
      if (fieldTypes != NULL) {
        fieldDofMap_.add_field(fieldIDs[i], fieldSizes[i], fieldTypes[i]);
      }
      else {
        fieldDofMap_.add_field(fieldIDs[i], fieldSizes[i]);
      }
    }
  }

  return(FEI_SUCCESS);
}
 
//------------------------------------------------------------------------------
int SNL_FEI_Structure::initElemBlock(GlobalID elemBlockID,
				    int numElements,
				    int numNodesPerElement,
				    const int* numFieldsPerNode,
				    const int* const* nodalFieldIDs,
				    int numElemDofFieldsPerElement,
				    const int* elemDofFieldIDs,
				    int interleaveStrategy)
{
  int nn, nf;
  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "FEI: initElemBlock" << FEI_ENDL << "#elemBlockID" << FEI_ENDL
       << (int)elemBlockID << FEI_ENDL;
    os << "#numElements"<< FEI_ENDL << numElements << FEI_ENDL;
    os << "#numNodesPerElement"<< FEI_ENDL <<numNodesPerElement << FEI_ENDL;
    os << "#numFieldsPerNode -- one entry per node " << FEI_ENDL;
    for(nn=0; nn<numNodesPerElement; nn++) os << numFieldsPerNode[nn]<<" ";
    os << FEI_ENDL << "#nodalFieldIDs -- one row per node" << FEI_ENDL;
    for(nn=0; nn<numNodesPerElement; ++nn) {
      for(nf=0; nf<numFieldsPerNode[nn]; ++nf) os << nodalFieldIDs[nn][nf] << " ";
      os << FEI_ENDL;
    }
    os << "#numElemDofFieldsPerElement" << FEI_ENDL
       << numElemDofFieldsPerElement<<FEI_ENDL;
    os << "#elemDofFieldIDs -- 'numElemDofFieldsPerElement' entries" << FEI_ENDL;
    for(nn=0; nn<numElemDofFieldsPerElement; ++nn) {
      os << elemDofFieldIDs[nn] << " ";
    }
    if (numElemDofFieldsPerElement > 0) os << FEI_ENDL;
    os << "#interleaveStrategy" << FEI_ENDL << interleaveStrategy << FEI_ENDL;
  }
   int j;

   CHK_ERR( addBlock(elemBlockID) );
   BlockDescriptor* block = NULL;
   CHK_ERR( getBlockDescriptor(elemBlockID, block) );

   CHK_ERR( block->setNumNodesPerElement(numNodesPerElement) );
   block->setNumElements(numElements);
   block->setElemDofFieldIDs(numElemDofFieldsPerElement, elemDofFieldIDs);
   block->setInterleaveStrategy(interleaveStrategy);

   int *fieldsPerNodePtr = block->fieldsPerNodePtr();

// construct the list of nodal solution cardinalities for this block

   int numNodalEqns = 0;
   int countDOF;
   std::vector<int> distinctFields(0);

   for(j=0; j<numNodesPerElement; j++) {

      countDOF = 0;
      for(int k = 0; k < numFieldsPerNode[j]; k++) {
        fei::sortedListInsert(nodalFieldIDs[j][k], distinctFields);

         int fieldSize = getFieldSize(nodalFieldIDs[j][k]);
         if (fieldSize < 0) {
	   fei::console_out() << "SNL_FEI_Structure::initElemBlock ERROR: fieldID " << 
	     nodalFieldIDs[j][k] << " has negative size. " << FEI_ENDL;
	   ERReturn(-1);
	 }
         countDOF += fieldSize;
      }

      fieldsPerNodePtr[j] = numFieldsPerNode[j];
      numNodalEqns += countDOF;
   }

   block->setNumDistinctFields(distinctFields.size());

   int numElemDOFPerElement = 0;
   for(j=0; j<numElemDofFieldsPerElement; j++) {
      int fieldSize = getFieldSize(elemDofFieldIDs[j]);
      if (fieldSize < 0) {
	fei::console_out() << "SNL_FEI_Structure::initElemBlock ERROR: elemDoffieldID " << 
	  elemDofFieldIDs[j] << " has negative size. " << FEI_ENDL;
	ERReturn(-1);
      }
      numElemDOFPerElement += fieldSize;
   }

   block->setNumElemDOFPerElement(numElemDOFPerElement);
   block->setNumEqnsPerElement(numNodalEqns + numElemDOFPerElement);
   block->setNumBlkEqnsPerElement(numNodesPerElement + numElemDOFPerElement);

// cache a copy of the element fields array for later use...

   CHK_ERR( block->allocateFieldIDsTable() );
   int **fieldIDsTablePtr = block->fieldIDsTablePtr();

   for (j = 0; j < numNodesPerElement; j++) {
      for(int k = 0; k < numFieldsPerNode[j]; k++) {
         fieldIDsTablePtr[j][k] = nodalFieldIDs[j][k];
      }
   }

//  create data structures for storage of element ID and topology info

   if (numElements > 0) {
      CHK_ERR( allocateBlockConnectivity(elemBlockID) );
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initElem(GlobalID elemBlockID,
			       GlobalID elemID,
			       const GlobalID* elemConn)
{
  if (debugOutput_ && outputLevel_ > 2) {
    FEI_OSTREAM& os = dbgOut();
    os << "FEI: initElem"<<FEI_ENDL;
    os << "#blkID"<<FEI_ENDL<<(int)elemBlockID<<FEI_ENDL<<"#elmID"<<FEI_ENDL
       <<(int)elemID<< FEI_ENDL;
  }

  //first get the block-descriptor for this elemBlockID...

  BlockDescriptor* block = NULL;
  CHK_ERR( getBlockDescriptor(elemBlockID, block) );

  ConnectivityTable& connTable = getBlockConnectivity(elemBlockID);

  std::map<GlobalID,int>& elemIDList = connTable.elemIDs;
  GlobalID* conn = &((*connTable.elem_conn_ids)[0]);

  int elemIndex = elemIDList.size();
  std::map<GlobalID,int>::iterator
    iter = elemIDList.find(elemID);
  
  bool redundantElem = false;

  if (iter != elemIDList.end()) {
    elemIndex = iter->second;
    redundantElem = true;
  }
  else {
    elemIDList.insert(std::make_pair(elemID,elemIndex));
  }

  int numNodes = block->getNumNodesPerElement();

  if (debugOutput_ && outputLevel_ > 2) {
    FEI_OSTREAM& os = dbgOut();
    os << "#n-nodes"<<FEI_ENDL<<numNodes<<FEI_ENDL<<"#nodeIDs"<<FEI_ENDL;
    for(int nn=0; nn<numNodes; nn++) os << (int)elemConn[nn] << " ";
    os << FEI_ENDL;
  }

  if (redundantElem) {
    //redundantElem == true means this element has been initialized before.
    //So we'll simply make sure the connectivity is the same as it was, and
    //if it isn't return -1.

    int offset = elemIndex*numNodes;
    for(int j=0; j<numNodes; j++) {
      if ( conn[offset+j] != elemConn[j]) {
	fei::console_out() << "SNL_FEI_Structure::initElem ERROR, elemID " << (int)elemID
	     << " registered more than once, with differing connectivity."
	     << FEI_ENDL;
	return(-1);
      }
    }
  }
  else {
    int offset = elemIndex*numNodes;
    for(int j = 0; j < numNodes; j++) {
      conn[offset+j] = elemConn[j];
    }

    CHK_ERR( nodeDatabase_->initNodeIDs(&(conn[offset]), numNodes) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initSlaveVariable(GlobalID slaveNodeID, 
					 int slaveFieldID,
					 int offsetIntoSlaveField,
					 int numMasterNodes,
					 const GlobalID* masterNodeIDs,
					 const int* masterFieldIDs,
					 const double* weights,
					 double rhsValue)
{
  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "FEI: initSlaveVariable" << FEI_ENDL;
    os << "#slvNodeID" << FEI_ENDL << (int)slaveNodeID << FEI_ENDL
       << "#slvFieldID"<< FEI_ENDL << slaveFieldID << FEI_ENDL
       << "#offsetIntoSlvField" << FEI_ENDL << offsetIntoSlaveField << FEI_ENDL
       << "#numMasterNodes" << FEI_ENDL << numMasterNodes << FEI_ENDL
       << "#masterNodeIDs" << FEI_ENDL;
    int nn;
    for(nn=0; nn<numMasterNodes; ++nn) {
      os <<(int)masterNodeIDs[nn]<<" ";
    }
    os << FEI_ENDL << "#masterFieldIDs" << FEI_ENDL;
    for(nn=0; nn<numMasterNodes; ++nn) {
      os <<masterFieldIDs[nn] << " ";
    }
    os << FEI_ENDL << "#field-sizes" << FEI_ENDL;
    for(nn=0; nn<numMasterNodes; ++nn) {
      int size = getFieldSize(masterFieldIDs[nn]);
      os << size << " ";
    }
    os << FEI_ENDL << "#weights" << FEI_ENDL;
    int offset = 0;
    for(nn=0; nn<numMasterNodes; ++nn) {
      int size = getFieldSize(masterFieldIDs[nn]);
      for(int j=0; j<size; ++j) {
	os << weights[offset++] << " ";
      }
    }
    os << FEI_ENDL << "#rhsValue" << FEI_ENDL << rhsValue << FEI_ENDL;
  }

  //Create and initialize a slave-variable object with the incoming data,
  //and add it to our list.

  SlaveVariable* svar = new SlaveVariable;
  svar->setNodeID(slaveNodeID);
  svar->setFieldID(slaveFieldID);
  svar->setFieldOffset(offsetIntoSlaveField);

  int woffset = 0;

  for(int i=0; i<numMasterNodes; i++) {
    svar->addMasterNodeID(masterNodeIDs[i]);
    svar->addMasterField(masterFieldIDs[i]);
    int fieldSize = getFieldSize(masterFieldIDs[i]);
    if (fieldSize < 0) ERReturn(-1);

    for(int j=0; j<fieldSize; j++) {
      svar->addWeight(weights[woffset++]);
    }
  }

  addSlaveVariable(svar);

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::deleteMultCRs()
{
  int numMCRs = multCRs_.size();
  if (numMCRs > 0) {
    fei::destroyValues(multCRs_);
    multCRs_.clear();
  }
  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initSharedNodes(int numSharedNodes,
				      const GlobalID *sharedNodeIDs,  
				      const int* numProcsPerNode, 
				      const int *const *sharingProcIDs)
{
  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "FEI: initSharedNodes"<<FEI_ENDL;
    os << "#n-nodes"<<FEI_ENDL<<numSharedNodes<<FEI_ENDL;
    os << "#num-procs-per-node"<<FEI_ENDL;
    int nn;
    for(nn=0; nn<numSharedNodes; ++nn) {
      os << numProcsPerNode[nn] << " ";
    }
    os << FEI_ENDL << "#following lines: nodeID   sharing-procs" << FEI_ENDL;
    for(nn=0; nn<numSharedNodes; nn++) {
      os << (int)sharedNodeIDs[nn]<<"   ";
      for(int np=0; np<numProcsPerNode[nn]; np++) os<<sharingProcIDs[nn][np]<<" ";
      os << FEI_ENDL;
    }
    os << "#end shared nodes"<<FEI_ENDL;
  }

  //  In this function we simply accumulate the incoming data into the
  //  NodeCommMgr object.
  //
  CHK_ERR( nodeCommMgr_->addSharedNodes(sharedNodeIDs, numSharedNodes,
					sharingProcIDs, numProcsPerNode) );

  if (activeNodesInitialized_) {
    //if active nodes have already been initialized, then we won't be
    //re-running the element-connectivities during the next call to
    //initComplete(), which means we won't have an opportunity to call the
    //NodeCommMgr::informLocal() method for nodes that reside locally. So we
    //need to check now for any locally residing shared nodes and make that
    //call. This is expensive, but it's a case that only arises when constraints
    //change between solves and nothing else changes. In general, constraints
    //are the only structural features allowed to change without requiring a
    //complete destruction and re-calculation of the FEI structure.
    //
    for(int i=0; i<numSharedNodes; ++i) {
      for(size_t nc=0; nc<connTables_.size(); ++nc) {
	if (connTables_[nc]->elem_conn_ids == NULL) continue;
	int len = connTables_[nc]->elem_conn_ids->size();
	if (len < 1) continue;
	GlobalID* conn_ids = &((*connTables_[nc]->elem_conn_ids)[0]);
	NodeDescriptor** nodes = &((*connTables_[nc]->elem_conn_ptrs)[0]);

	for(int j=0; j<len; ++j) {
	  if (conn_ids[j] == sharedNodeIDs[i]) {
	    CHK_ERR( nodeCommMgr_->informLocal(*(nodes[j])));
	    break;
	  }
	}
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initCRMult(int numCRNodes,
				 const GlobalID* CRNodes,
				 const int *CRFields,
				 int& CRID)
{
  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "FEI: initCRMult" << FEI_ENDL << "# numCRNodes: "<<FEI_ENDL<<numCRNodes<<FEI_ENDL;
    os << "#CRNodes:"<<FEI_ENDL;
    int nn;
    for(nn=0; nn<numCRNodes; nn++) {
      os << (int)CRNodes[nn]<<" ";
    }
    os << FEI_ENDL<<"#fields:"<<FEI_ENDL;
    for(nn=0; nn<numCRNodes; nn++) {
      os << CRFields[nn]<<" ";
    }
    os<<FEI_ENDL;
  }
  //  tasks: add local constraint equations, determine sparsity
  //         patterns for the new constraint relation
  //

  CRID = localProc_*100000 + multCRs_.size();
  ConstraintType* multCRPtr = NULL;
  addCR(CRID, multCRPtr, multCRs_);

  ConstraintType& multCR = *multCRPtr;

  multCR.setConstraintID(CRID);

  std::vector<GlobalID>& CRNodeArray = multCR.getMasters();
  std::vector<int>& CRFieldArray = multCR.getMasterFieldIDs();

  for(int j = 0; j < numCRNodes; j++) {
    CRNodeArray.push_back(CRNodes[j]);
    CRFieldArray.push_back(CRFields[j]);
  }

  if (debugOutput_) dbgOut() << "#(output) CRID:"<<FEI_ENDL << CRID << FEI_ENDL;

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initCRPen(int numCRNodes,
				const GlobalID* CRNodes,
				const int *CRFields,
				int& CRID)
{
  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "initCRPen, numCRNodes: "<<numCRNodes<<FEI_ENDL;
    for(int nn=0; nn<numCRNodes; nn++) {
      os << "   crNodeID: "<<(int)CRNodes[nn]<<", field: "<<CRFields[nn]<<FEI_ENDL;
    }
  }

  CRID = localProc_*100000 + penCRs_.size();

  ConstraintType* penCRPtr = NULL;
  addCR(CRID, penCRPtr, penCRs_);

  ConstraintType& penCR = *penCRPtr;

  penCR.setConstraintID(CRID);
  penCR.setIsPenalty(true);

  std::vector<GlobalID>& CRNodesArray = penCR.getMasters();

  std::vector<int>& CRFieldArray = penCR.getMasterFieldIDs();

  for(int i = 0; i < numCRNodes; i++) {
    CRNodesArray.push_back(CRNodes[i]);
    CRFieldArray.push_back(CRFields[i]);
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
bool SNL_FEI_Structure::isInLocalElement(int nodeNumber)
{
  //I'm going to make an assumption: if we're running in serial, then there
  //are no remote anythings, and the node in question must be in a local
  //element.
  if (numProcs_ < 2) return(true);

  const NodeDescriptor* node = NULL;
  int err = nodeDatabase_->getNodeWithNumber(nodeNumber, node);
  if (err != 0) return(false);

  GlobalID nodeID = node->getGlobalNodeID();

  //now let's find out if this node is a shared node.
  int shIndx = nodeCommMgr_->getSharedNodeIndex(nodeID);
  if (shIndx < 0) {
    //if shIndx < 0, then this isn't a shared node. For now, we will assume
    //that it must be a local node.
    return(true);
  }

  //If we reach this line, then the node is a shared node. Let's now ask if
  //it appears locally...
  std::vector<GlobalID>& localNodes = nodeCommMgr_->getLocalNodeIDs();
  int index = fei::binarySearch(nodeID, &localNodes[0], localNodes.size());
  if (index >= 0) return(true);

  //If we reach this line, then the node is shared, but doesn't appear in the
  //"localNodeIDs" held by NodeCommMgr. This means it is not in a local element,
  //so we'll return false.
  return(false);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initComplete(bool generateGraph)
{
  //This is the most significant function in SNL_FEI_Structure. This is where
  //the underlying matrix structure is calculated from all the data that has
  //been provided since construction. Calculating the equation space and
  //forming the matrix structure is a multi-step process, which proceeds like
  //this:
  //
  //1. finalizeActiveNodes()
  //  - makes sure that all shared nodes have been identified to the
  //    NodeDatabase, then has the NodeDatabase allocate its internal list of
  //    NodeDescriptor objects.
  //  - runs the element-connectivity tables and sets the nodal field and block
  //    info on the NodeDescriptors in the NodeDatabase. IMPORTANT: At this
  //    point, the nodeIDs in the connectivitiy tables are replaced with
  //    indices into the nodeDatabase to speed future lookups.
  //  - As the connectivity tables are being run, the NodeCommMgr is given the
  //    nodeID of each node that is connected to a local element, and the node's
  //    owner is set to the initial value of localProc_.
  //1a. finalizeNodeCommMgr()
  //  - NodeCommMgr::initComplete is called, at which time the ownership of all
  //    shared nodes is determined.
  //
  //2. Next, all local nodal degrees of freedom can be counted. This, together
  //  with the number of lagrange multipliers and element-centered DOFs are used
  //  by the function calcGlobalEqnInfo(), which determines global equation
  //  counts and offsets.
  //
  //3. setNodalEqnInfo() runs all locally active nodes and for each one that's
  //  locally owned, sets global equation-numbers on it, as well as its
  //  numNodalDOF and block-entry equation-number.
  //
  //4. setElemDOFEqnInfo() and setMultCREqnInfo() associate global equation-
  //  numbers with each of the lagrange multipliers and element-dofs.
  //
  //5. NodeDatabase::synchronize() is called, which associates nodeNumbers with
  //  each node, and also in turn calls NodeCommMgr::exchangeEqnInfo() to obtain
  //  equation-numbers and field/block info for remotely-owned shared nodes.
  //
  //6. setNumNodesAndEqnsPerBlock() then runs the active nodes again and counts
  //  active nodes and equations per element-block. This couldn't be done before
  //  because NodeCommMgr::exchangeEqnInfo() may have found out about new blocks
  //  that are associated with shared nodes on other processors.
  //
  //7. calculateSlaveEqns() associates equation-numbers with slave-nodes and
  //  master-nodes that have been identified locally, then gathers slave
  //  equation info from all processors onto all other processors, so that all
  //  processors know about all slave equations.
  //
  //8. initializeEqnCommMgr() runs the shared-nodes from NodeCommMgr and lets
  //  EqnCommMgr know which equations will be sent to other processors, and
  //  which will be received from other processors.
  //
  //9. translate localStartRow_ and localEndRow_ into the reduced equation
  //  space, where they are called reducedStartRow_ and reducedEndRow_, after
  //  which we know how many local equations there are not including any
  //  slave equations. At this point we can allocate the sysMatIndices_ array
  //  of arrays, which will hold the point-entry matrix structure.
  //
  //10. formMatrixStructure()
  //  - initElemBlockStructure() for each element in each element-block,
  //    obtains scatter-indices and makes the corresponding structurally
  //    symmetric contributions to the matrix structure. This process includes
  //    calculations associated with reducing out slave equations, and with
  //    putting off-processor contributions (element contributions with shared
  //    nodes) into the EqnCommMgr object.
  //  - initMultCRStructure() and initPenCRStructure() follow the same general
  //    routine as for the initElemBlockStructure() function.
  //  - EqnCommMgr::exchangeIndices() is called, which sends all shared
  //    contributions to the processors that own the corresponding equations,
  //    after which the receiving processors insert those contributions into
  //    the local matrix structure.
  //
  //11. initializeBlkEqnMapper() run all nodes, elem-dof, multiplier constraints
  //  and pass global point-equation-numbers with corresponding block-equation
  //  numbers to the snl_fei::PointBlockMap object which maintains a mapping
  //  between the point-entry and block-entry equation spaces.
  //
  //12. EqnCommMgr::exchangePtToBlkInfo() ensures that the point-to-block
  //  mapping data is globally consistent across all processors.
  //
  //13. The sysBlkMatIndices_ array is allocated, which is the array of
  //  arrays that will hold the block-entry matrix structure. This block-entry
  //  matrix structure is not filled, however, until later if/when the
  //  method getMatrixStructure is called specifically requesting block-entry
  //  structure data.
  //
  //14. Finally, if the debugOutputLevel_ is greater than 0, a file is written
  //  which contains a mapping from equations to nodes. This file is often very
  //  useful for post-mortem debugging operations, and is also necessary if the
  //  matrix produced from a serial run is to be compared with a matrix 
  //  produced from a similar run on multiple processors. The equation-to-node
  //  mapping is the only way to reconcile the equation-orderings between the
  //  two matrices.
  //

  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "FEI: initComplete" << FEI_ENDL;
  }
  // marshall all of the nodes we received in element connectivities into an
  // active node list, and set the fields, blocks and owning processors 
  // associated with those nodes.

  //IMPORTANT!!! Note that at this point, the nodeIDs in the
  //connectivitiy tables are being replaced with indices from the nodeDatabase.

  CHK_ERR( finalizeActiveNodes() );

  CHK_ERR( finalizeNodeCommMgr() );

  numLocalElemDOF_ = calcTotalNumElemDOF();
  numLocalMultCRs_ = calcNumMultCREqns();

  // ok, now the active equation calculations -- we need to count how many
  // local equations there are, then obtain the global equation info, then
  // set equation numbers on all of our active nodes.

  int numLocalNodes     = nodeDatabase_->countLocalNodeDescriptors(localProc_);
  numLocalNodalEqns_ = nodeDatabase_->countLocalNodalEqns(localProc_);

  numLocalEqns_ = numLocalNodalEqns_ + numLocalElemDOF_ + numLocalMultCRs_;
  numLocalEqnBlks_ = numLocalNodes   + numLocalElemDOF_ + numLocalMultCRs_;

  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    os << "#   numMultCRs: " << numLocalMultCRs_ << ", numElemDOF: "
       << numLocalElemDOF_ << ", numLocalNodalEqns: " <<numLocalNodalEqns_
       << FEI_ENDL;
    os << "#   numLocalEqns_: " << numLocalEqns_ << FEI_ENDL;
  }

  calcGlobalEqnInfo(numLocalNodes, numLocalEqns_, numLocalEqnBlks_);

  CHK_ERR( setNodalEqnInfo() );

  //now we need to set equation numbers for the element-DOF's, and for the
  //lagrange-multipler constraints.
  //(exactly where to put these element DOF is an interesting issue... here, 
  //just toss them at the end of the nodal active eqns, which may or may not
  //be such a great choice.)

  setElemDOFEqnInfo();

  CHK_ERR( setMultCREqnInfo() );

  //now have the NodeDatabase run through the list of active nodes and
  //do whatever final initializations it needs to do. This function also
  //calls nodeCommMgr::exchangeEqnInfo.
  CHK_ERR( nodeDatabase_->synchronize(firstLocalNodeNumber_, localStartRow_,
				      localProc_, comm_) );

  //now run the active nodes again and set the number of
  //active nodes and equations per block. We couldn't do this
  //before, because the nodeCommMgr::initComplete() and
  //nodeCommMgr::exchangeEqnInfo functions may both have found out about new
  //blocks that are associated with shared nodes.

  setNumNodesAndEqnsPerBlock();

  //at this point we'll run through any SlaveVariable records that were set,
  //and establish an EqnBuffer of slave equation numbers, and gather the slave
  //equations initialized on any other processors.

  slvCommMgr_ = new EqnCommMgr(comm_);
  slvCommMgr_->setNumRHSs(1);

  CHK_ERR( calculateSlaveEqns(comm_) );

  //From this point on, all computations involving equation numbers will assume
  //that those numbers are in the reduced equation space. So we'll start by
  //translating the contents of eqnCommMgr_ to the reduced space.
  //CHK_ERR( translateToReducedEqns(*eqnCommMgr_) );

  initializeEqnCommMgr();

  translateToReducedEqn(localEndRow_, reducedEndRow_);
  int startRow = localStartRow_;
  while(isSlaveEqn(startRow)) startRow++;
  translateToReducedEqn(startRow, reducedStartRow_);
  numLocalReducedRows_ = reducedEndRow_ - reducedStartRow_ + 1;

  generateGraph_ = generateGraph;

  if (generateGraph_) {
    sysMatIndices_ = new fei::ctg_set<int>[numLocalReducedRows_];
  }

  //Now we'll run through all the data structures that have been initialized
  //and form the matrix structure (lists of column-indices).

  CHK_ERR( formMatrixStructure() );

  CHK_ERR( initializeBlkEqnMapper() );

  if (globalMaxBlkSize_ > 1) {
    CHK_ERR( eqnCommMgr_->exchangePtToBlkInfo(*blkEqnMapper_) );
    if (numSlvs_ > 0) {
      try {
      CHK_ERR( slvCommMgr_->exchangeIndices() );
      }
      catch (std::runtime_error& exc) {
        fei::console_out() << exc.what() << FEI_ENDL;
        ERReturn(-1);
      }

      CHK_ERR( slvCommMgr_->exchangePtToBlkInfo(*blkEqnMapper_) );
    }
  }

  localReducedBlkOffset_ = ptEqnToBlkEqn(reducedStartRow_);
  int lastLocalReducedBlkEqn = ptEqnToBlkEqn(reducedEndRow_);
  numLocalReducedEqnBlks_ = lastLocalReducedBlkEqn - localReducedBlkOffset_ + 1;

  if (globalMaxBlkSize_ > 1 && generateGraph_) {
    sysBlkMatIndices_ = numLocalReducedEqnBlks_>0 ?
      new fei::ctg_set<int>[numLocalReducedEqnBlks_] : NULL;
  }

  delete slvCommMgr_;

  blockMatrix_ = true;

  structureFinalized_ = true;
  matIndicesDestroyed_ = false;

  //while we're here, let's write an equation -> nodeNumber map into a file
  //for possible debugging purposes... it will come in handy if we need to
  //compare matrices/vectors from different parallel runs.
  if (debugOutput_) writeEqn2NodeMap();

  if (debugOutput_) dbgOut() << "#leaving initComplete" << FEI_ENDL;
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::formMatrixStructure()
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) os << "#   formMatrixStructure" << FEI_ENDL;

  //do our local equation profile calculations (i.e., determine
  //how long each row is). use the element scatter arrays to determine 
  //the sparsity pattern

  CHK_ERR( initElemBlockStructure() );

  // next, handle the matrix structure imposed by the constraint relations...
  //

  CHK_ERR( initMultCRStructure() );

  // we also need to accomodate penalty constraints, so let's process them
  // now...

  CHK_ERR( initPenCRStructure() );

  //we've now processed all of the local data that produces matrix structure.
  //so now let's have the equation comm manager exchange indices with the
  //other processors so we can account for remote contributions to our
  //matrix structure.
  try {
  CHK_ERR( eqnCommMgr_->exchangeIndices(&os) );
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  //so now the remote contributions should be available, let's get them out
  //of the eqn comm mgr and put them into our local matrix structure.

  int numRecvEqns = eqnCommMgr_->getNumLocalEqns();
  std::vector<int>& recvEqnNumbers = eqnCommMgr_->localEqnNumbers();
  std::vector<fei::CSVec*>& recvEqns = eqnCommMgr_->localEqns();
  int i;
  if (debugOutput_) {
    os << "#     after eqnCommMgr_->exchangeIndices, numRecvEqns: "
       << numRecvEqns << FEI_ENDL;
  }

  for(i=0; i<numRecvEqns; i++) {
    int eqn = recvEqnNumbers[i];
    if ((reducedStartRow_ > eqn) || (reducedEndRow_ < eqn)) {
      fei::console_out() << "SNL_FEI_Structure::initComplete: ERROR, recvEqn " << eqn
	   << " out of range. (reducedStartRow_: " << reducedStartRow_
	   << ", reducedEndRow_: " << reducedEndRow_ << ", localProc_: "
	   << localProc_ << ")" << FEI_ENDL;
      return(-1);
    }

    for(size_t j=0; j<recvEqns[i]->size(); j++) {
      CHK_ERR( createMatrixPosition(eqn, recvEqns[i]->indices()[j],
				    "frmMatStr") );
    }
  }

  if (debugOutput_) {
    os << "#  leaving formMatrixStructure" << FEI_ENDL;
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initElemBlockStructure()
{
  int numBlocks = getNumElemBlocks();

  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) {
    os << "#     initElemBlockStructure" << FEI_ENDL;
    os << "#        numElemBlocks: " << numBlocks << FEI_ENDL;
  }

  for (int bIndex = 0; bIndex < numBlocks; bIndex++) {
    BlockDescriptor* block = NULL;
    CHK_ERR( getBlockDescriptor_index(bIndex, block) );

    int numEqns = block->getNumEqnsPerElement();
    int interleave = block->getInterleaveStrategy();
    std::vector<int> scatterIndices(numEqns);

    GlobalID elemBlockID = block->getGlobalBlockID();
    ConnectivityTable& connTable = getBlockConnectivity(elemBlockID);

    std::map<GlobalID,int>& elemIDList = connTable.elemIDs;
    block->setNumElements(elemIDList.size());
    int numBlockElems = block->getNumElements();

    //loop over all the elements, determining the elemental (both from nodes 
    //and from element DOF) contributions to the sparse matrix structure
    if (debugOutput_) {
      os << "#        block " << bIndex << ", numElems: " << numBlockElems<<FEI_ENDL;
    }
    for(int elemIndex = 0; elemIndex < numBlockElems; elemIndex++) {

      getScatterIndices_index(bIndex, elemIndex, interleave,
                              &scatterIndices[0]);

      //now, store the structure that will arise from this contribution,
      //after first performing any manipulations associated with slave
      //reduction.

      CHK_ERR( createSymmEqnStructure(scatterIndices) );
    }
  }

  if (reducedEqnCounter_ > 0) CHK_ERR( assembleReducedStructure() );

  if (debugOutput_) os << "#     leaving initElemBlockStructure" << FEI_ENDL;
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getMatrixRowLengths(std::vector<int>& rowLengths)
{
  if (!structureFinalized_) return(-1);

  rowLengths.resize(numLocalReducedRows_);

  for(int i=0; i<numLocalReducedRows_; i++) {
    rowLengths[i] = sysMatIndices_[i].size();
  }
  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getMatrixStructure(int** indices,
					 std::vector<int>& rowLengths)
{
  if (!structureFinalized_) return(-1);

  rowLengths.resize(numLocalReducedRows_);

  for(int i=0; i<numLocalReducedRows_; i++) {
    rowLengths[i] = sysMatIndices_[i].size();
    fei::copySetToArray(sysMatIndices_[i], rowLengths[i], indices[i]);
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getMatrixStructure(int** ptColIndices,
					 std::vector<int>& ptRowLengths,
					 int** blkColIndices,
					  int* blkIndices_1D,
					 std::vector<int>& blkRowLengths,
					 std::vector<int>& numPtRowsPerBlkRow)
{
  int err = getMatrixStructure(ptColIndices, ptRowLengths);
  if (err != 0) return(-1);

  if (globalMaxBlkSize_ == 1) {
    //No block-equations have more than one point-equation, so we'll just assign
    //the block-structure arrays to be the same as the point-structure arrays.
    int numRows = ptRowLengths.size();
    blkRowLengths.resize(numRows);
    numPtRowsPerBlkRow.assign(numRows, 1);

    blkRowLengths.resize(ptRowLengths.size());
    for(size_t ii=0; ii<ptRowLengths.size(); ++ii) {
      blkRowLengths[ii] = ptRowLengths[ii];
    }

    for(int i=0; i<numRows; i++) {
      blkColIndices[i] = ptColIndices[i];
    }
  }
  else {
    //There are block-equations with more than 1 point-equation, so let's
    //calculate the actual block-structure.

    std::map<int,int>* ptEqns = blkEqnMapper_->getPtEqns();
    int numPtEqns = ptEqns->size();

    std::map<int,int>::const_iterator
      pteq = ptEqns->begin(),
      pteq_end = ptEqns->end();

    int lastBlkRow = -1;

    for(int jj=0; jj<numPtEqns; ++jj, ++pteq) {
      int ptEqn = (*pteq).first;
      int localPtEqn = ptEqn - reducedStartRow_;
      if (localPtEqn < 0 || localPtEqn >= numLocalReducedRows_) continue;

      int rowLength = sysMatIndices_[localPtEqn].size();

      int blkRow = blkEqnMapper_->eqnToBlkEqn(ptEqn);
      if (blkRow < 0) {
	ERReturn(-1);
      }

      if (blkRow == lastBlkRow) continue;

      int localBlkRow = blkRow - localReducedBlkOffset_;
      if (localBlkRow < 0 || localBlkRow >= numLocalReducedEqnBlks_) {
	ERReturn(-1);
      }

      fei::ctg_set<int>& sysBlkIndices = sysBlkMatIndices_[localBlkRow];

      int* theseColIndices = ptColIndices[localPtEqn];

      for(int colj=0; colj<rowLength; colj++) {
	int blkCol = blkEqnMapper_->eqnToBlkEqn(theseColIndices[colj]);
	if (blkCol < 0) {
	  fei::console_out() << localProc_
	       <<"SNL_FEI_Structure::getMatrixStructure ERROR pt row "
	       << ptEqn << ", pt col "
	       << ptColIndices[localPtEqn][colj]
	       << " doesn't have a corresponding block" << FEI_ENDL;
	  blkCol = blkEqnMapper_->eqnToBlkEqn(theseColIndices[colj]);
	  std::abort();
	}

	sysBlkIndices.insert2(blkCol);
      }

      lastBlkRow = blkRow;
    }

    //now we're ready to set the block-sizes...

    numPtRowsPerBlkRow.resize(numLocalReducedEqnBlks_);
    blkRowLengths.resize(numLocalReducedEqnBlks_);
    int offset = 0;
    for(int i=0; i<numLocalReducedEqnBlks_; i++) {
      blkRowLengths[i] = sysBlkMatIndices_[i].size();
      fei::copySetToArray(sysBlkMatIndices_[i],blkRowLengths[i],
				       &(blkIndices_1D[offset]));
      offset += blkRowLengths[i];

      int blkEqn = localReducedBlkOffset_ + i;
      numPtRowsPerBlkRow[i] = blkEqnMapper_->getBlkEqnSize(blkEqn);
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
bool SNL_FEI_Structure::nodalEqnsAllSlaves(const NodeDescriptor* node,
					   std::vector<int>& slaveEqns)
{
  int numFields = node->getNumFields();
  const int* fieldIDs = node->getFieldIDList();
  const int* fieldEqns= node->getFieldEqnNumbers();

  for(int i=0; i<numFields; ++i) {
    int fieldSize = getFieldSize(fieldIDs[i]);
    for(int eqn=0; eqn<fieldSize; ++eqn) {
      int thisEqn = fieldEqns[i] + eqn;
      if (fei::binarySearch(thisEqn, slaveEqns) < 0) {
	//found a nodal eqn that's not a slave eqn, so return false.
	return(false);
      }
    }
  }

  return(true);
}

//------------------------------------------------------------------------------
NodeDescriptor& SNL_FEI_Structure::findNodeDescriptor(GlobalID nodeID)
{
//
//This function returns a NodeDescriptor reference if nodeID is an active node.
//
  NodeDescriptor* node = NULL;
  int err = nodeDatabase_->getNodeWithID(nodeID, node);

  if (err != 0) {
    fei::console_out() << "ERROR, findNodeDescriptor unable to find node " << (int)nodeID
	 << FEI_ENDL;
    std::abort();
  }

  return( *node );
}

//------------------------------------------------------------------------------
NodeDescriptor* SNL_FEI_Structure::findNode(GlobalID nodeID)
{
  NodeDescriptor* node = NULL;
  nodeDatabase_->getNodeWithID(nodeID, node);
  return node;
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initMultCRStructure() 
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) os << "#    initMultCRStructure" << FEI_ENDL;
  //
  //Private SNL_FEI_Structure function, to be called from initComplete.
  //
  //Records the system matrix structure arising from Lagrange Multiplier
  //Constraint Relations.
  //
  // since at initialization all we are being passed is the
  // general form of the constraint equations, we can't check to see if any
  // of the weight terms are zeros at this stage of the game.  Hence, we
  // have to reserve space for all the nodal weight vectors, even though
  // they might turn out to be zeros during the load step....

  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = multCRs_.begin(),
    cr_end  = multCRs_.end();

  while(cr_iter != cr_end) {
    ConstraintType& multCR = *((*cr_iter).second);

    int lenList = multCR.getMasters().size();

    std::vector<GlobalID>& CRNode_vec = multCR.getMasters();
    GlobalID *CRNodePtr = &CRNode_vec[0];
    std::vector<int>& CRField_vec = multCR.getMasterFieldIDs();
    int* CRFieldPtr = &CRField_vec[0];

    int crEqn = multCR.getEqnNumber();
    int reducedCrEqn = 0;
    translateToReducedEqn(crEqn, reducedCrEqn);

    createMatrixPosition(reducedCrEqn, reducedCrEqn, "initMCRStr");

    for(int j=0; j<lenList; j++) {
      GlobalID nodeID = CRNodePtr[j];
      int fieldID = CRFieldPtr[j];

      NodeDescriptor *nodePtr = findNode(nodeID);
      if(nodePtr==NULL) continue;
      NodeDescriptor& node = *nodePtr;

      //first, store the column indices associated with this node, in
      //the constraint's equation. i.e., position (crEqn, node)

      storeNodalColumnIndices(crEqn, node, fieldID);

      //now we need to store the transpose of the above contribution,
      //i.e., position(s) (node, crEqn)

      if (node.getOwnerProc() == localProc_) {
	//if we own this node, we will simply store its equation
	//numbers in local rows (equations) of the matrix.

	storeNodalRowIndices(node, fieldID, crEqn);
      }
      else {
	//if we don't own it, then we need to register with the
	//eqn comm mgr that we'll be sending contributions to
	//column crEqn of the remote equations associated with this
	//node.

	storeNodalSendIndex(node, fieldID, crEqn);
      }
    }
    ++cr_iter;
  }
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initPenCRStructure()
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) os << "#     initPenCRStructure" << FEI_ENDL;
//
//This function oversees the putting in of any matrix structure that results
//from Penalty constraint relations.
//
// note that penalty constraints don't generate new equations
// (unlike Lagrange constraints), but they do add terms to the system
// stiffness matrix that couple all the nodes that contribute to the
// penalty constraint.  In addition, each submatrix is defined by the pair
// of nodes that creates its contribution, hence a submatrix can be defined
// in terms of two weight vectors (of length p and q) instead of the
// generally larger product matrix (of size pq)

// the additional terms take the form of little submatrices that look a lot 
// like an element stiffness and load matrix, where the nodes in the
// constraint list take on the role of the nodes associated with an
// element, and the individual matrix terms arise from the outer products
// of the constraint nodal weight vectors 

// rather than get elegant and treat this task as such an elemental energy
// term, we'll use some brute force to construct these penalty contributions
// on the fly, primarily to simplify -reading- this thing, so that the 
// derivations in the annotated implementation document are more readily
// followed...
  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = penCRs_.begin(),
    cr_end = penCRs_.end();

  while(cr_iter != cr_end) {
    ConstraintType& penCR = *((*cr_iter).second);

    int lenList = penCR.getMasters().size();
    std::vector<GlobalID>& CRNode_vec = penCR.getMasters();
    GlobalID* CRNodesPtr = &CRNode_vec[0];

    std::vector<int>& CRField_vec = penCR.getMasterFieldIDs();
    int* CRFieldPtr = &CRField_vec[0];

    // each constraint equation generates a set of nodal energy terms, so
    // we have to process a matrix of nodes for each constraint

    for(int i = 0; i < lenList; i++) {
      GlobalID iNodeID = CRNodesPtr[i];
      int iField = CRFieldPtr[i];

      NodeDescriptor* iNodePtr = findNode(iNodeID);
      if(!iNodePtr) continue;
      NodeDescriptor& iNode = *iNodePtr;

      for(int j = 0; j < lenList; j++) {
	GlobalID jNodeID = CRNodesPtr[j];
	int jField = CRFieldPtr[j];

	NodeDescriptor* jNodePtr = findNode(jNodeID);
        if(!jNodePtr) continue;
        NodeDescriptor &jNode = *jNodePtr;

	if (iNode.getOwnerProc() == localProc_) {
	  //if iNode is local, we'll put equations into the local 
	  //matrix structure.

	  storeLocalNodeIndices(iNode, iField, jNode, jField);
	}
	else {
	  //if iNode is remotely owned, we'll be registering equations
	  //to send to its owning processor.

	  storeNodalSendIndices(iNode, iField, jNode, jField);
	}
      }
    }   //   end i loop
    ++cr_iter;
   }   //   end while loop
   return(FEI_SUCCESS);
}
 
//------------------------------------------------------------------------------
void SNL_FEI_Structure::storeNodalSendIndex(NodeDescriptor& node, int fieldID,
					    int col)
{
  //
  //This is a private SNL_FEI_Structure function. We can safely assume that it
  //will only be called with a node that is not locally owned.
  //
  //This function tells the eqn comm mgr that we'll be sending contributions
  //to column 'col' for the equations associated with 'fieldID', on 'node', on
  //node's owning processor.
  //
  int proc = node.getOwnerProc();

  int eqnNumber = -1;
  if (!node.getFieldEqnNumber(fieldID, eqnNumber)) voidERReturn;

  int numEqns = getFieldSize(fieldID);
  if (numEqns < 1) {
    fei::console_out() << "FEI error, attempt to store indices for field ("<<fieldID
	 <<") with size "<<numEqns<<FEI_ENDL;
    voidERReturn;
  }

  for(int i=0; i<numEqns; i++) {
    eqnCommMgr_->addRemoteIndices(eqnNumber+i, proc, &col, 1);
  }
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::storeNodalSendIndices(NodeDescriptor& iNode, int iField,
                                  NodeDescriptor& jNode, int jField){
//
//This function will register with the eqn comm mgr the equations associated
//with iNode, field 'iField' having column indices that are the equations
//associated with jNode, field 'jField', to be sent to the owner of iNode.
//
   int proc = iNode.getOwnerProc();

   int iEqn = -1, jEqn = -1;
   if (!iNode.getFieldEqnNumber(iField, iEqn)) return; //voidERReturn;
   if (!jNode.getFieldEqnNumber(jField, jEqn)) return; //voidERReturn;

   int iNumParams = getFieldSize(iField);
   int jNumParams = getFieldSize(jField);
   if (iNumParams < 1 || jNumParams < 1) {
     fei::console_out() << "FEI ERROR, attempt to store indices for field with non-positive size"
	  << " field "<<iField<<", size "<<iNumParams<<", field "<<jField<<", size "
	  << jNumParams<<FEI_ENDL;
     voidERReturn;
   }

   for(int i=0; i<iNumParams; i++) {
      int eqn = iEqn + i;

      for(int j=0; j<jNumParams; j++) {
         int col = jEqn + j;
         eqnCommMgr_->addRemoteIndices(eqn, proc, &col, 1);
      }
   }
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::storeLocalNodeIndices(NodeDescriptor& iNode, int iField,
					     NodeDescriptor& jNode, int jField)
{
//
//This function will add to the local matrix structure the equations associated
//with iNode at iField, having column indices that are the equations associated
//with jNode at jField.
//
   int iEqn = -1, jEqn = -1;
   if (!iNode.getFieldEqnNumber(iField, iEqn)) return; //voidERReturn;
   if (!jNode.getFieldEqnNumber(jField, jEqn)) return; //voidERReturn;

   int iNumParams = getFieldSize(iField);
   int jNumParams = getFieldSize(jField);
   if (iNumParams < 1 || jNumParams < 1) {
     fei::console_out() << "FEI ERROR, attempt to store indices for field with non-positive size"
	  << " field "<<iField<<", size "<<iNumParams<<", field "<<jField<<", size "
	  << jNumParams<<FEI_ENDL;
     voidERReturn;
   }

   for(int i=0; i<iNumParams; i++) {
      int row = iEqn + i;
      int reducedRow = -1;
      bool isSlave = translateToReducedEqn(row, reducedRow);
      if (isSlave) continue;
      
      for(int j=0; j<jNumParams; j++) {
         int col = jEqn + j;
	 int reducedCol = -1;
	 isSlave = translateToReducedEqn(col, reducedCol);
	 if (isSlave) continue;

         createMatrixPosition(reducedRow, reducedCol,
			      "storeLocNdInd");
      }
   }
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::storeNodalColumnIndices(int eqn, NodeDescriptor& node,
					       int fieldID)
{
  //
  //This function stores the equation numbers associated with 'node' at
  //'fieldID' as column indices in row 'eqn' of the system matrix structure.
  //
  if ((localStartRow_ > eqn) || (eqn > localEndRow_)) {
    return;
  }

  int colNumber = -1;
  if (!node.getFieldEqnNumber(fieldID, colNumber)) voidERReturn;

  int numParams = getFieldSize(fieldID);
  if (numParams < 1) {
    fei::console_out() << "FEI error, attempt to store indices for field ("<<fieldID
	 <<") with size "<<numParams<<FEI_ENDL;
    voidERReturn;
  }

  int reducedEqn = -1;
  bool isSlave = translateToReducedEqn(eqn, reducedEqn);
  if (isSlave) return;

  for(int j=0; j<numParams; j++) {
    int reducedCol = -1;
    isSlave = translateToReducedEqn(colNumber+j, reducedCol);
    if (isSlave) continue;

    createMatrixPosition(reducedEqn, reducedCol,
			 "storeNdColInd");
  }
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::storeNodalRowIndices(NodeDescriptor& node,
					    int fieldID, int eqn)
{
  //
  //This function stores column 'eqn' in equation numbers associated with
  //'node' at 'fieldID' in the system matrix structure. 
  //
  int eqnNumber = -1;
  if (!node.getFieldEqnNumber(fieldID, eqnNumber)) voidERReturn;

  int numParams = getFieldSize(fieldID);
  if (numParams < 1) {
    fei::console_out() << "FEI error, attempt to store indices for field ("<<fieldID
	 <<") with size "<<numParams<<FEI_ENDL;
    voidERReturn;
  }

  int reducedEqn = -1;
  bool isSlave = translateToReducedEqn(eqn, reducedEqn);
  if (isSlave) return;

  for(int j=0; j<numParams; j++) {
    int reducedRow = -1;
    isSlave = translateToReducedEqn(eqnNumber+j, reducedRow);
    if (isSlave) continue;

    createMatrixPosition(reducedRow, reducedEqn, "storeNdRowInd");
  }
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::createSymmEqnStructure(std::vector<int>& scatterIndices)
{
  //scatterIndices are both the rows and the columns for a structurally
  //symmetric 2-dimensional contribution to the matrix.
  //

  //if there aren't any slaves in the whole problem, store the scatter-indices
  //using a fast function (no translations-to-reduced-space) and return.
  if (numSlaveEquations() == 0) {
    CHK_ERR( storeElementScatterIndices_noSlaves(scatterIndices) );
    return(FEI_SUCCESS);
  }

  try {

  int len = scatterIndices.size();
  bool anySlaves = false;
  rSlave_.resize(len);
  for(int is=0; is<len; is++) { 
    rSlave_[is] = isSlaveEqn(scatterIndices[is]) ? 1 : 0;
    if (rSlave_[is] == 1) anySlaves = true;
  }

  //if there aren't any slaves in this contribution, then just store it
  //and return
  if (!anySlaves) {
    CHK_ERR( storeElementScatterIndices(scatterIndices) );
    return(FEI_SUCCESS);
  }

  int* scatterPtr = &scatterIndices[0];

  workSpace_.resize(len);
  for(int j=0; j<len; j++) {
    translateToReducedEqn(scatterPtr[j], workSpace_[j]);
  }
  
  for(int i=0; i<len; i++) {
    int row = scatterPtr[i];
    if (rSlave_[i] == 1) {
      reducedEqnCounter_++;
      //
      //'row' is a slave equation, so add this row to Kdi_. But as we do that,
      //watch for columns that are slave equations and add them to Kid_.
      //
      for(int jj=0; jj<len; jj++) {
	int col = scatterPtr[jj];

	if (rSlave_[jj] == 1) {
	  //'col' is a slave equation, so add this column to Kid_.
	  for(int ii=0; ii<len; ii++) {
	    int rowi = scatterPtr[ii];

	    //only add the non-slave rows for this column.
	    if (rSlave_[ii] == 1) continue;

	    Kid_->createPosition(rowi, col);
	  }
	  //now skip to the next column 
	  continue;
	}

	Kdi_->createPosition(row, col);
      }

      //finally, add all slave columns in this slave row to Kdd_.
      for(int kk=0; kk<len; kk++) {
	int colk = scatterPtr[kk];

	if (rSlave_[kk] != 1) continue;

	Kdd_->createPosition(row, colk);
      }
    }
    else {
      //this is not a slave row, so we will loop across it, creating matrix
      //positions for all non-slave columns in it.
      int reducedRow = -1;
      translateToReducedEqn(row, reducedRow);

      bool rowIsLocal = true;
      int owner = localProc_;
      if (reducedStartRow_ > reducedRow || reducedEndRow_ < reducedRow) {
	rowIsLocal = false;
	owner = getOwnerProcForEqn(row);
      }

      for(int j=0; j<len; j++) {
	if (rSlave_[j] == 1) continue;

	int reducedCol = workSpace_[j];

	if (rowIsLocal) {
	  CHK_ERR( createMatrixPosition(reducedRow, reducedCol,
					"crtSymmEqnStr") );
	}
	else {
	  eqnCommMgr_->addRemoteIndices(reducedRow, owner, &reducedCol, 1);
	}
      }
    }
  }

  if (reducedEqnCounter_ > 300) CHK_ERR( assembleReducedStructure() );

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::createBlkSymmEqnStructure(std::vector<int>& scatterIndices)
{
  //scatterIndices are both the rows and the columns for a structurally
  //symmetric 2-dimensional contribution to the matrix.
  //

  //if there aren't any slaves in the whole problem, store the scatter-indices
  //using a fast function (no translations-to-reduced-space) and return.
  if (numSlaveEquations() == 0) {
    return( storeElementScatterBlkIndices_noSlaves(scatterIndices) );
  }

  try {

  int len = scatterIndices.size();
  bool anySlaves = false;
  rSlave_.resize(len);
  for(int is=0; is<len; is++) { 
    rSlave_[is] = isSlaveEqn(scatterIndices[is]) ? 1 : 0;
    if (rSlave_[is] == 1) anySlaves = true;
  }

  //if there aren't any slaves in this contribution, then just store it
  //and return
  if (!anySlaves) {
    CHK_ERR( storeElementScatterIndices(scatterIndices) );
    return(FEI_SUCCESS);
  }

  int* scatterPtr = &scatterIndices[0];

  workSpace_.resize(len);
  for(int j=0; j<len; j++) {
    translateToReducedEqn(scatterPtr[j], workSpace_[j]);
  }
  
  for(int i=0; i<len; i++) {
    int row = scatterPtr[i];
    if (rSlave_[i] == 1) {
      reducedEqnCounter_++;
      //
      //'row' is a slave equation, so add this row to Kdi_. But as we do that,
      //watch for columns that are slave equations and add them to Kid_.
      //
      for(int jj=0; jj<len; jj++) {
	int col = scatterPtr[jj];

	if (rSlave_[jj] == 1) {
	  //'col' is a slave equation, so add this column to Kid_.
	  for(int ii=0; ii<len; ii++) {
	    int rowi = scatterPtr[ii];

	    //only add the non-slave rows for this column.
	    if (rSlave_[ii] == 1) continue;

	    Kid_->createPosition(rowi, col);
	  }
	  //now skip to the next column 
	  continue;
	}

	Kdi_->createPosition(row, col);
      }

      //finally, add all slave columns in this slave row to Kdd_.
      for(int kk=0; kk<len; kk++) {
	int colk = scatterPtr[kk];

	if (rSlave_[kk] != 1) continue;

	Kdd_->createPosition(row, colk);
      }
    }
    else {
      //this is not a slave row, so we will loop across it, creating matrix
      //positions for all non-slave columns in it.
      int reducedRow = -1;
      translateToReducedEqn(row, reducedRow);

      bool rowIsLocal = true;
      int owner = localProc_;
      if (reducedStartRow_ > reducedRow || reducedEndRow_ < reducedRow) {
	rowIsLocal = false;
	owner = getOwnerProcForEqn(row);
      }

      for(int j=0; j<len; j++) {
	if (rSlave_[j] == 1) continue;

	int reducedCol = workSpace_[j];

	if (rowIsLocal) {
	  CHK_ERR( createMatrixPosition(reducedRow, reducedCol,
					"crtSymmEqnStr") );
	}
	else {
	  eqnCommMgr_->addRemoteIndices(reducedRow, owner, &reducedCol, 1);
	}
      }
    }
  }

  if (reducedEqnCounter_ > 300) CHK_ERR( assembleReducedStructure() );

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::storeElementScatterIndices(std::vector<int>& indices)
{
//This function takes a list of equation numbers, as input, and for each
//one, if it is a local equation number, goes to that row of the sysMatIndices
//structure and stores the whole list of equation numbers as column indices
//in that row of the matrix structure. If the equation number is not local,
//then it is given to the EqnCommMgr.
//
   int numIndices = indices.size();
   int* indPtr = &indices[0];

   workSpace_.resize(numIndices);
   int* workPtr = &workSpace_[0];
   int reducedEqn = -1;
   for(size_t j=0; j<indices.size(); j++) {
     bool isSlave = translateToReducedEqn(indPtr[j], reducedEqn);
     if (isSlave) ERReturn(-1);
     workPtr[j] = reducedEqn;
   }

   for(int i=0; i<numIndices; i++) {
      int row = indPtr[i];
      int rrow = workPtr[i];

      if ((localStartRow_ > row) || (row > localEndRow_)) {
	
	int owner = getOwnerProcForEqn(row);
	eqnCommMgr_->addRemoteIndices(rrow, owner, workPtr, numIndices);
	continue;
      }

      CHK_ERR( createMatrixPositions(rrow, numIndices, workPtr,
				     "storeElemScttrInd") );
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::storeElementScatterIndices_noSlaves(std::vector<int>& indices)
{
  //This function takes a list of equation numbers as input and for each one,
  //if it is a local equation number, goes to that row of the sysMatIndices
  //structure and stores the whole list of equation numbers as column indices
  //in that row of the matrix structure. If the equation number is not local,
  //then it is given to the EqnCommMgr.
  //
  int i, numIndices = indices.size();
  int* indPtr = &indices[0];

  for(i=0; i<numIndices; i++) {
    int row = indPtr[i];
    if (row < localStartRow_ || row > localEndRow_) {
      int owner = getOwnerProcForEqn(row);
      eqnCommMgr_->addRemoteIndices(row, owner, indPtr, numIndices);
      continue;
    }

    if (generateGraph_) {
      fei::ctg_set<int>& thisRow =
	sysMatIndices_[row - reducedStartRow_];

      for(int j=0; j<numIndices; ++j) {
	if (debugOutput_ && outputLevel_ > 2) {
	  dbgOut() << "#       storeElemScttrInds_ns crtMatPoss("
		   << row << "," << indPtr[j] << ")"<<FEI_ENDL;
	}

	thisRow.insert2(indPtr[j]);
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::storeElementScatterBlkIndices_noSlaves(std::vector<int>& indices)
{
  //This function takes a list of equation numbers as input and for each one,
  //if it is a local equation number, goes to that row of the sysMatIndices
  //structure and stores the whole list of equation numbers as column indices
  //in that row of the matrix structure. If the equation number is not local,
  //then it is given to the EqnCommMgr.
  //
  int i, numIndices = indices.size();
  int* indPtr = &indices[0];

  for(i=0; i<numIndices; i++) {
    int row = indPtr[i];
    if (row < localStartRow_ || row > localEndRow_) {
      int owner = getOwnerProcForEqn(row);
      eqnCommMgr_->addRemoteIndices(row, owner, indPtr, numIndices);
      continue;
    }

    if (generateGraph_) {
      fei::ctg_set<int>& thisRow =
	sysMatIndices_[row - reducedStartRow_];

      for(int j=0; j<numIndices; ++j) {
	if (debugOutput_ && outputLevel_ > 2) {
	  dbgOut() << "#       storeElemScttrInds_ns crtMatPoss("
		   << row << "," << indPtr[j] << ")"<<FEI_ENDL;
	}

	thisRow.insert2(indPtr[j]);
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::assembleReducedStructure()
{
  //This function performs the appropriate manipulations (matrix-matrix-products,
  //matrix-vector-products) on the Kid,Kdi,Kdd matrices and then assembles the
  //results into the global matrix structure. This function also adds any
  //resulting contributions to off-processor portions of the matrix structure to
  //the EqnCommMgr.
  //
  fei::FillableMat* D = getSlaveDependencies();

  csrD = *D;
  csrKid = *Kid_;
  csrKdi = *Kdi_;
  csrKdd = *Kdd_;

  //form tmpMat1_ = Kid_*D
  fei::multiply_CSRMat_CSRMat(csrKid, csrD, tmpMat1_);

  //form tmpMat2_ = D^T*Kdi_
  fei::multiply_trans_CSRMat_CSRMat(csrD, csrKdi, tmpMat2_);

  CHK_ERR( translateMatToReducedEqns(tmpMat1_) );
  CHK_ERR( translateMatToReducedEqns(tmpMat2_) );
  if (numProcs_ > 1) {
    CHK_ERR( eqnCommMgr_->addRemoteEqns(tmpMat1_, true) );
    CHK_ERR( eqnCommMgr_->addRemoteEqns(tmpMat2_, true) );
  }

  CHK_ERR( createMatrixPositions(tmpMat1_) );
  CHK_ERR( createMatrixPositions(tmpMat2_) );

  //form tmpMat1_ = D^T*Kdd_
  fei::multiply_trans_CSRMat_CSRMat(csrD, csrKdd, tmpMat1_);

  //form tmpMat2_ = tmpMat1_*D = D^T*Kdd_*D
  fei::multiply_CSRMat_CSRMat(tmpMat1_, csrD, tmpMat2_);


  CHK_ERR( translateMatToReducedEqns(tmpMat2_) );
  if (numProcs_ > 1) {
    CHK_ERR( eqnCommMgr_->addRemoteEqns(tmpMat2_, true) );
  }
  CHK_ERR( createMatrixPositions(tmpMat2_) );

  Kdi_->clear();
  Kid_->clear();
  Kdd_->clear();
  reducedEqnCounter_ = 0;

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::translateToReducedEqns(EqnCommMgr& eqnCommMgr)
{
  if (numSlvs_ < 1) return(0);

  CHK_ERR( translateToReducedEqns(*(eqnCommMgr.getRecvProcEqns())) );
  CHK_ERR( translateToReducedEqns(*(eqnCommMgr.getRecvEqns())) );
  CHK_ERR( translateToReducedEqns(*(eqnCommMgr.getSendProcEqns())) );
  CHK_ERR( translateToReducedEqns(*(eqnCommMgr.getSendEqns())) );

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::translateToReducedEqns(EqnBuffer& eqnBuf)
{
  int numEqns = eqnBuf.getNumEqns();
  int* eqnNumbers = &(eqnBuf.eqnNumbers()[0]);
  std::vector<fei::CSVec*>& eqnArray = eqnBuf.eqns();
  for(int i=0; i<numEqns; ++i) {
    int reducedEqn;
    translateToReducedEqn(eqnNumbers[i], reducedEqn);
    eqnNumbers[i] = reducedEqn;

    int* indicesPtr = &(eqnArray[i]->indices()[0]);
    int numIndices = eqnArray[i]->size();
    for(int j=0; j<numIndices; ++j) {
      translateToReducedEqn(indicesPtr[j], reducedEqn);
      indicesPtr[j] = reducedEqn;
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::translateToReducedEqns(ProcEqns& procEqns)
{
  std::vector<std::vector<int>*>& eqnTable = procEqns.procEqnNumbersPtr();
  int dim1 = eqnTable.size();
  for(int i=0; i<dim1; ++i) {
    int dim2 = eqnTable[i]->size();
    int* indicesPtr = &(*(eqnTable[i]))[0];
    for(int j=0; j<dim2; ++j) {
      int reducedEqn;
      translateToReducedEqn(indicesPtr[j], reducedEqn);
      indicesPtr[j] = reducedEqn;
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::translateMatToReducedEqns(fei::CSRMat& mat)
{
  //Given a matrix in unreduced numbering, convert all of its contents to the
  //"reduced" equation space. If any of the row-numbers or column-indices in
  //this matrix object are slave equations, they will not be referenced. In
  //this case, a positive warning-code will be returned.

  int reducedEqn = -1;
  int foundSlave = 0;

  fei::SparseRowGraph& srg = mat.getGraph();

  std::vector<int>& rowNumbers = srg.rowNumbers;
  for(size_t i=0; i<rowNumbers.size(); ++i) {
    bool isSlave = translateToReducedEqn(rowNumbers[i], reducedEqn);
    if (isSlave) foundSlave = 1;
    else rowNumbers[i] = reducedEqn;
  }

  std::vector<int>& colIndices = srg.packedColumnIndices;
  for(size_t i=0; i<colIndices.size(); ++i) {
    bool isSlave = translateToReducedEqn(colIndices[i], reducedEqn);
    if (isSlave) foundSlave = 1;
    else colIndices[i] = reducedEqn;
  }

  return(foundSlave);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::createMatrixPositions(fei::CSRMat& mat)
{
  if (!generateGraph_) {
    return(0);
  }

  //This function must be called with a CSRMat object that already has its
  //contents numbered in "reduced" equation-numbers.
  //
  //This function has one simple task -- for each row,col pair stored in 'mat',
  //call 'createMatrixPosition' to make an entry in the global matrix structure
  //if it doesn't already exist.
  //
  const std::vector<int>& rowNumbers = mat.getGraph().rowNumbers;
  const std::vector<int>& rowOffsets = mat.getGraph().rowOffsets;
  const std::vector<int>& pckColInds = mat.getGraph().packedColumnIndices;

  for(size_t i=0; i<rowNumbers.size(); ++i) {
    int row = rowNumbers[i];
    int offset = rowOffsets[i];
    int rowlen = rowOffsets[i+1]-offset;
    const int* indices = &pckColInds[offset];

    for(int j=0; j<rowlen; j++) {
      CHK_ERR( createMatrixPosition(row, indices[j], "crtMatPos(CSRMat)") );
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::createMatrixPosition(int row, int col,
					    const char* callingFunction)
{
  if (!generateGraph_) {
    return(0);
  }

  //
  //This function inserts the column index 'col' in the system matrix structure
  //in 'row', if it isn't already there.
  //
  //This is a private SNL_FEI_Structure function. These row/col numbers must be
  //global, 0-based equation numbers in the "reduced" equation space..
  //

  //if the equation isn't local, just return. (The assumption being that
  //this call is also being made on the processor that owns the equation.)
  if (row < reducedStartRow_ || row > reducedEndRow_) {
    if (debugOutput_) {
      dbgOut() << "         " << callingFunction << " passed " <<row<<","<<col
	       << ", not local." << FEI_ENDL;
    }

    return(0);
  }

  if (debugOutput_ && outputLevel_ > 2) {
    dbgOut() << "#         " << callingFunction << " crtMatPos("
	     << row << "," << col << ")"<<FEI_ENDL;
  }

  fei::ctg_set<int>& thisRow = sysMatIndices_[row - reducedStartRow_];

  thisRow.insert2(col);

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::createMatrixPositions(int row, int numCols, int* cols,
					    const char* callingFunction)
{
  if (!generateGraph_) {
    return(0);
  }

  //
  //This function inserts the column indices 'cols' in the system matrix structure
  //in 'row', if they aren't already there.
  //
  //This is a private SNL_FEI_Structure function. These row/col numbers must be
  //global, 0-based equation numbers in the "reduced" equation space..
  //

  //if the equation isn't local, just return. (The assumption being that
  //this call is also being made on the processor that owns the equation.)
  if (row < reducedStartRow_ || row > reducedEndRow_) {
    if (debugOutput_) {
      dbgOut() << "         " << callingFunction << " passed " <<row
	       << ", not local." << FEI_ENDL;
    }

    return(0);
  }

  fei::ctg_set<int>& thisRow = sysMatIndices_[row - reducedStartRow_];

  for(int i=0; i<numCols; i++) {
    if (debugOutput_ && outputLevel_ > 2) {
      dbgOut() << "#         " << callingFunction << " crtMatPoss("
	       << row << "," << cols[i] << ")"<<FEI_ENDL;
    }

    thisRow.insert2(cols[i]);
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::initializeBlkEqnMapper()
{
  //Run the nodes, elem-dof's and multiplier constraints, establishing the
  //point-equation to block-equation mappings.  Note that this data is 
  //initialized as non-reduced (i.e., not accounting for any slave equations)
  //equation numbers, then translated to the reduced space at the end of this
  //function.

  if (blkEqnMapper_ != NULL) delete blkEqnMapper_;
  blkEqnMapper_ = new snl_fei::PointBlockMap();
  snl_fei::PointBlockMap& blkEqnMapper = getBlkEqnMapper();

  //First, if there's only one (non-negative) fieldID and that fieldID's size
  //is 1, then we don't need to execute the rest of this function.
  int numFields = fieldDatabase_->size();
  const int* fieldIDs = getFieldIDsPtr();
  const int* fieldSizes = fieldIDs+numFields;
  int numNonNegativeFields = 0;
  int maxFieldSize = 0;
  for(int f=0; f<numFields; f++) {
    if (fieldIDs[f] >= 0) {
      numNonNegativeFields++;
      if (fieldSizes[f] > maxFieldSize) maxFieldSize = fieldSizes[f];
    }
  }

  if (numNonNegativeFields == 1 && maxFieldSize == 1) {
    globalMaxBlkSize_ = 1;
    blkEqnMapper.setPtEqualBlk();
    return(0);
  }

  int localVanishedNodeAdjustment = 0;
  for(int i=0; i<localProc_; ++i) {
    localVanishedNodeAdjustment += globalNumNodesVanished_[i];
  }

  int eqnNumber, blkEqnNumber;

  std::map<GlobalID,int>& nodeIDmap = nodeDatabase_->getNodeIDs();
  std::map<GlobalID,int>::iterator
    iter = nodeIDmap.begin(), iter_end = nodeIDmap.end();

  for(; iter!=iter_end; ++iter) {
    NodeDescriptor* node = NULL;
    nodeDatabase_->getNodeAtIndex(iter->second, node);

    // If the node doesn't exist, skip.
    if (node==NULL) continue;

    // If the node doesn't have any fields, skip
    int numFields  = node->getNumFields();
    if(numFields == 0) continue;

    int firstPtEqn = node->getFieldEqnNumbers()[0];
    int numNodalDOF = node->getNumNodalDOF();
    blkEqnNumber = node->getBlkEqnNumber() - localVanishedNodeAdjustment;

    int blkSize = numNodalDOF;

    for(int j=0; j<numNodalDOF; ++j) {
      bool isSlave = translateToReducedEqn(firstPtEqn+j, eqnNumber);
      if (isSlave) {
	--blkSize;
      }
      else {
	CHK_ERR( blkEqnMapper.setEqn(eqnNumber, blkEqnNumber) );
      }
    }

    if (blkSize > 0) {
      CHK_ERR( blkEqnMapper.setBlkEqnSize(blkEqnNumber, blkSize) );
    }
    else {
      ++localVanishedNodeAdjustment;
    }
  }

  //
  //Now the element-dofs...
  //
  int numBlocks = getNumElemBlocks();
  int numLocalNodes = globalNodeOffsets_[localProc_+1]-globalNodeOffsets_[localProc_];
  eqnNumber = localStartRow_ + numLocalNodalEqns_;
  blkEqnNumber = localBlkOffset_ + numLocalNodes;

  for(int i = 0; i < numBlocks; i++) {
    BlockDescriptor* block = NULL;
    CHK_ERR( getBlockDescriptor_index(i, block) );

    int numElemDOF = block->getNumElemDOFPerElement();

    if (numElemDOF > 0) {
      int numElems = block->getNumElements();

      for(int j=0; j<numElems; j++) {

	for(int ede=0; ede<numElemDOF; ede++) {
	  blkEqnMapper.setEqn(eqnNumber+ede, blkEqnNumber+ede);
	  blkEqnMapper.setBlkEqnSize(blkEqnNumber+ede, 1);
	}

	eqnNumber += numElemDOF;
	blkEqnNumber += numElemDOF;
      }
    }
  }

  //Now we'll set mappings for all local multiplier constraint-relations,
  //and also gather the equation numbers for other processors' multipliers
  //to record in our snl_fei::PointBlockMap object.
  //
  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = multCRs_.begin(),
    cr_end  = multCRs_.end();

  std::vector<int> localMultEqns;
  localMultEqns.reserve(multCRs_.size()*2);
  while(cr_iter != cr_end) {
    ConstraintType* multCR = (*cr_iter).second;
    eqnNumber = multCR->getEqnNumber();
    blkEqnNumber = multCR->getBlkEqnNumber();

    CHK_ERR( blkEqnMapper_->setEqn(eqnNumber, blkEqnNumber) );
    CHK_ERR( blkEqnMapper_->setBlkEqnSize(blkEqnNumber, 1) );

    localMultEqns.push_back(eqnNumber);
    localMultEqns.push_back(blkEqnNumber);
    ++cr_iter;
  }

  //Now gather and store the equation-numbers arising from other
  //processors' multiplier constraints. (We obviously only need to do this if
  //there is more than one processor, and we can skip it if the problem
  //only contains 1 scalar equation for each block-equation.)

  int localMaxBlkEqnSize = blkEqnMapper_->getMaxBlkEqnSize();
  globalMaxBlkSize_ = 0;
  CHK_ERR( fei::GlobalMax(comm_, localMaxBlkEqnSize, globalMaxBlkSize_) );

  blkEqnMapper_->setMaxBlkEqnSize(globalMaxBlkSize_);

  if (globalMaxBlkSize_ == 1) {
    //If the globalMax block-size is 1 then tell the blkEqnMapper, and it will
    //reclaim the memory it allocated in arrays, etc.
    blkEqnMapper_->setPtEqualBlk();
    return(0);
  }

  if (numProcs_ > 1) {
    std::vector<int> recvLengths, globalMultEqns;
    CHK_ERR(fei::Allgatherv(comm_, localMultEqns, recvLengths, globalMultEqns));

    int offset = 0;
    for(int p=0; p<numProcs_; p++) {
      if (p == localProc_) { offset += recvLengths[p]; continue; }

      for(int j=0; j<recvLengths[p]/2; j++) {
	blkEqnMapper_->setEqn(globalMultEqns[offset], globalMultEqns[offset+1]);
	blkEqnMapper_->setBlkEqnSize(globalMultEqns[offset+1], 1);
	offset += 2;
      }
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::setMultCREqnInfo()
{
  //We'll set equation-numbers on all local multiplier constraint-relations,
  //and also gather the equation numbers for other processors' multipliers
  //to record in our snl_fei::PointBlockMap object.
  //
  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = multCRs_.begin(),
    cr_end  = multCRs_.end();

  int eqnNumber = localStartRow_ + numLocalNodalEqns_ + numLocalElemDOF_;
  int blkEqnNum = localBlkOffset_ + numLocalEqnBlks_ - numLocalMultCRs_;

  while(cr_iter != cr_end) {
    ConstraintType* multCR = (*cr_iter).second;

    multCR->setEqnNumber(eqnNumber);

    multCR->setBlkEqnNumber(blkEqnNum);

    ++eqnNumber;
    ++blkEqnNum;
    ++cr_iter;
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::writeEqn2NodeMap()
{
  FEI_OSTRINGSTREAM osstr;
  osstr << dbgPath_ << "/FEI" << name_ << "_nID_nNum_blkEq_ptEq."
	<< numProcs_ << "." << localProc_;

  FEI_OFSTREAM e2nFile(osstr.str().c_str());

  if (!e2nFile) {
    fei::console_out() << "SNL_FEI_Structure::writeEqn2NodeMap: ERROR, couldn't open file "
         << osstr.str() << FEI_ENDL;
    ERReturn(-1);
  }

  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = multCRs_.begin(),
    cr_end  = multCRs_.end();

  while(cr_iter != cr_end) {
    ConstraintType* multCR = (*cr_iter).second;
    int eqnNumber = multCR->getEqnNumber();
    int blkEqnNumber = multCR->getBlkEqnNumber();

    e2nFile << "-999 -999 " << blkEqnNumber << " " << eqnNumber << FEI_ENDL;
    ++cr_iter;
  }

  std::map<GlobalID,int>& nodeIDmap = nodeDatabase_->getNodeIDs();
  std::map<GlobalID,int>::iterator
    iter = nodeIDmap.begin(), iter_end = nodeIDmap.end();

  for(; iter!=iter_end; ++iter) {
    NodeDescriptor* node = NULL;
    nodeDatabase_->getNodeAtIndex(iter->second, node);

    if (node==NULL || node->getOwnerProc() != localProc_) {
      continue;
    }

    const int* fieldIDList = node->getFieldIDList();
    const int* eqnNumbers = node->getFieldEqnNumbers();
    int nodeNum = node->getNodeNumber();
    int blkEqnNumber = node->getBlkEqnNumber();
    int numFields = node->getNumFields();

    for(int j=0; j<numFields; j++) {
      int numFieldParams = getFieldSize(fieldIDList[j]);
      assert( numFieldParams > 0 );

      for(int k=0; k<numFieldParams; k++) {
        e2nFile << (int)node->getGlobalNodeID() << " " << nodeNum
                << " " << blkEqnNumber << " " << eqnNumbers[j]+k<<FEI_ENDL;
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::calcTotalNumElemDOF()
{
   int numBlocks = getNumElemBlocks();
   int totalNumElemDOF = 0;

   for(int i = 0; i < numBlocks; i++) {
     BlockDescriptor* block = NULL;
     CHK_ERR( getBlockDescriptor_index(i, block) );
      int numElemDOF = block->getNumElemDOFPerElement();
      if (numElemDOF > 0) {
         int numElems = block->getNumElements();

         totalNumElemDOF += numElems*numElemDOF;
      }
   }

   return(totalNumElemDOF);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::calcNumMultCREqns()
{
  int totalNumMultCRs = 0;
  int numMCRecords = getNumMultConstRecords();

  totalNumMultCRs = numMCRecords;

  return( totalNumMultCRs );
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::calcGlobalEqnInfo(int numLocallyOwnedNodes, 
					 int numLocalEqns,
					 int numLocalEqnBlks) 
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) os << "#  entering calcGlobalEqnInfo" << FEI_ENDL;

//
//This function does the inter-process communication necessary to obtain,
//on each processor, the global number of equations, and the local starting
//and ending equation numbers.
//While we're here, we do the same thing for node-numbers.
//

  //first, get each processor's local number of equations on the master proc.

  std::vector<int> numProcNodes(numProcs_);
  std::vector<int> numProcEqns(numProcs_);
  std::vector<int> numProcEqnBlks(numProcs_);

  if (numProcs_ > 1) {
#ifndef FEI_SER
    int glist[3];
    std::vector<int> recvList(3*numProcs_);

    glist[0] = numLocallyOwnedNodes;
    glist[1] = numLocalEqns;
    glist[2] = numLocalEqnBlks;
    if (MPI_Gather(&glist[0], 3, MPI_INT, &recvList[0], 3,
		   MPI_INT, masterProc_, comm_) != MPI_SUCCESS) {
      fei::console_out() << "SNL_FEI_Structure::calcGlobalEqnInfo: ERROR in MPI_Gather" << FEI_ENDL;
      MPI_Abort(comm_, -1);
    }

    for(int p=0; p<numProcs_; p++) {
      numProcNodes[p]   = recvList[p*3];
      numProcEqns[p]    = recvList[p*3+1];
      numProcEqnBlks[p] = recvList[p*3+2];
    }
#endif
  }
  else {
    numProcNodes[0] = numLocallyOwnedNodes;
    numProcEqns[0] = numLocalEqns;
    numProcEqnBlks[0] = numLocalEqnBlks;
  }

  //compute offsets for all processors (starting index for local equations)

  globalNodeOffsets_.resize(numProcs_+1);
  globalEqnOffsets_.resize(numProcs_+1);
  globalBlkEqnOffsets_.resize(numProcs_+1);

  globalNodeOffsets_[0] = 0; // 0-based node-numbers
  globalEqnOffsets_[0] = 0;    // we're going to start rows & cols at 0 (global)
  globalBlkEqnOffsets_[0] = 0;

  if (localProc_ == masterProc_) {
    for(int i=1;i<numProcs_;i++) {
      globalNodeOffsets_[i] = globalNodeOffsets_[i-1] +numProcNodes[i-1];
      globalEqnOffsets_[i] = globalEqnOffsets_[i-1] + numProcEqns[i-1];
      globalBlkEqnOffsets_[i] = globalBlkEqnOffsets_[i-1] +
                                  numProcEqnBlks[i-1];
    }

    globalNodeOffsets_[numProcs_] = globalNodeOffsets_[numProcs_-1] +
                                   numProcNodes[numProcs_-1];
    globalEqnOffsets_[numProcs_] = globalEqnOffsets_[numProcs_-1] +
                                  numProcEqns[numProcs_-1];
    globalBlkEqnOffsets_[numProcs_] = globalBlkEqnOffsets_[numProcs_-1] +
                                     numProcEqnBlks[numProcs_-1];
  }

  //now, broadcast vector of eqnOffsets to all processors

  if (numProcs_ > 1) {
#ifndef FEI_SER
    std::vector<int> blist(3*(numProcs_+1));

    if (localProc_ == masterProc_) {
      int offset = 0;
      for(int i=0; i<numProcs_+1; i++) {
	blist[offset++] = globalNodeOffsets_[i];
	blist[offset++] = globalEqnOffsets_[i];
	blist[offset++] = globalBlkEqnOffsets_[i];
      }
    }

    if (MPI_Bcast(&blist[0], 3*(numProcs_+1), MPI_INT,
		  masterProc_, comm_) != MPI_SUCCESS) {
      fei::console_out() << "SNL_FEI_Structure::calcGlobalEqnInfo: ERROR in MPI_Bcast" << FEI_ENDL;
      MPI_Abort(comm_, -1);
    }

    if (localProc_ != masterProc_) {
      int offset = 0;
      for(int i=0; i<numProcs_+1; i++) {
	globalNodeOffsets_[i] = blist[offset++];
	globalEqnOffsets_[i] = blist[offset++];
	globalBlkEqnOffsets_[i] = blist[offset++];
      }
    }
#endif
  }

  localStartRow_ = globalEqnOffsets_[localProc_];
  localEndRow_ = globalEqnOffsets_[localProc_+1]-1;
  numGlobalEqns_ = globalEqnOffsets_[numProcs_];

  firstLocalNodeNumber_ = globalNodeOffsets_[localProc_];
  numGlobalNodes_ = globalNodeOffsets_[numProcs_];

  localBlkOffset_ = globalBlkEqnOffsets_[localProc_];
  numGlobalEqnBlks_ = globalBlkEqnOffsets_[numProcs_];

  if (debugOutput_) {
    os << "#     localStartRow_: " << localStartRow_ << FEI_ENDL;
    os << "#     localEndRow_: "   << localEndRow_ << FEI_ENDL;
    os << "#     numGlobalEqns_: " << numGlobalEqns_ << FEI_ENDL;
    os << "#     numGlobalEqnBlks_: " << numGlobalEqnBlks_ << FEI_ENDL;
    os << "#     localBlkOffset_: " << localBlkOffset_ << FEI_ENDL;
    os << "#     firstLocalNodeNumber_: " << firstLocalNodeNumber_ << FEI_ENDL;
    size_t n;
    os << "#     numGlobalNodes_: " << numGlobalNodes_ << FEI_ENDL;
    os << "#     " << globalNodeOffsets_.size() << " globalNodeOffsets_: ";
    for(size_t nn=0; nn<globalNodeOffsets_.size(); nn++) 
      os << globalNodeOffsets_[nn] << " ";
    os << FEI_ENDL;
    os << "#     " << globalEqnOffsets_.size() << " globalEqnOffsets_: ";
    for(n=0; n<globalEqnOffsets_.size(); n++) 
      os << globalEqnOffsets_[n] << " ";
    os << FEI_ENDL;
    os << "#     " << globalBlkEqnOffsets_.size() << " globalBlkEqnOffsets_: ";
    for(n=0; n<globalBlkEqnOffsets_.size(); n++) 
      os << globalBlkEqnOffsets_[n] << " ";
    os << FEI_ENDL;
    os << "#  leaving calcGlobalEqnInfo" << FEI_ENDL;
  }
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::setNodalEqnInfo()
{
//
//Private SNL_FEI_Structure function, to be called in initComplete, after
//'calcGlobalEqnInfo()'.
//
//This function sets global equation numbers on the nodes.
//
//The return value is an error-code.
//
  int eqn = localStartRow_;

  //localBlkOffset_ is 0-based, and so is blkEqnNumber.
  int blkEqnNumber = localBlkOffset_;

  std::map<GlobalID,int>& nodeIDmap = nodeDatabase_->getNodeIDs();
  std::map<GlobalID,int>::iterator
    iter = nodeIDmap.begin(), iter_end = nodeIDmap.end();

  for(; iter!=iter_end; ++iter) {
    NodeDescriptor* node = NULL;
    nodeDatabase_->getNodeAtIndex(iter->second, node);

    if (node==NULL) continue;

    int numFields = node->getNumFields();
    const int* fieldIDs = node->getFieldIDList();

    if (node->getOwnerProc() == localProc_) {
      int numNodalDOF = 0;

      for(int j=0; j<numFields; j++) {
	node->setFieldEqnNumber(fieldIDs[j], eqn);
	int fieldSize = getFieldSize(fieldIDs[j]);
	eqn += fieldSize;
	numNodalDOF += fieldSize;
      }

      node->setNumNodalDOF(numNodalDOF);
      node->setBlkEqnNumber(blkEqnNumber++);
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::setElemDOFEqnInfo()
{
  //
  //Private SNL_FEI_Structure function, to be called from initComplete after
  //setBlkEqnInfo() has been called.
  //
  //Sets global equation numbers for all element-DOF.
  //
  int numBlocks = getNumElemBlocks();

  int eqnNumber = localStartRow_ + numLocalNodalEqns_;

  for(int i = 0; i < numBlocks; i++) {
    BlockDescriptor* block = NULL;
    int err = getBlockDescriptor_index(i, block);
    if (err) voidERReturn;

    int numElemDOF = block->getNumElemDOFPerElement();

    if (numElemDOF > 0) {
      std::vector<int>& elemDOFEqns = block->elemDOFEqnNumbers();
      std::map<GlobalID,int>& elemIDs = connTables_[i]->elemIDs;

      std::map<GlobalID,int>::const_iterator
        iter = elemIDs.begin(),
        iter_end = elemIDs.end();

      for(; iter != iter_end; ++iter) {
	elemDOFEqns[iter->second] = eqnNumber;

	eqnNumber += numElemDOF;
      }
    }
  }
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::addBlock(GlobalID blockID) {
//
//Append a blockID to our (unsorted) list of blockIDs if it isn't
//already present. Also, add a block-descriptor if blockID wasn't
//already present, and a ConnectivityTable to go with it.
//
  int insertPoint = -1;
  int found = fei::binarySearch(blockID, blockIDs_, insertPoint);

   if (found<0) {
      blockIDs_.insert(blockIDs_.begin()+insertPoint, blockID);

      BlockDescriptor* block = new BlockDescriptor;
      block->setGlobalBlockID(blockID);

      blocks_.insert(blocks_.begin()+insertPoint, block);

      ConnectivityTable* newConnTable = new ConnectivityTable;
      connTables_.insert(connTables_.begin()+insertPoint, newConnTable);
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getBlockDescriptor(GlobalID blockID,
                                         BlockDescriptor*& block)
{
  int index = fei::binarySearch(blockID, blockIDs_);

   if (index < 0) {
      fei::console_out() << "SNL_FEI_Structure::getBlockDescriptor: ERROR, blockID "
           << (int)blockID << " not found." << FEI_ENDL;
      return(-1);
   }

   block = blocks_[index];

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getIndexOfBlock(GlobalID blockID) const
{
  int index = fei::binarySearch(blockID, blockIDs_);
  return(index);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getBlockDescriptor_index(int index,
					       BlockDescriptor*& block)
{
  if (index < 0 || index >= (int)blockIDs_.size()) ERReturn(FEI_FATAL_ERROR);

  block = blocks_[index];

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::allocateBlockConnectivity(GlobalID blockID) {

   int index = fei::binarySearch(blockID, blockIDs_);

   if (index < 0) {
      fei::console_out() << "SNL_FEI_Structure::allocateConnectivityTable: ERROR, blockID "
           << (int)blockID << " not found. Aborting." << FEI_ENDL;
      MPI_Abort(comm_, -1);
   }

   connTables_[index]->numNodesPerElem = 
                         blocks_[index]->getNumNodesPerElement();

   int numRows = blocks_[index]->getNumElements();
   int numCols = connTables_[index]->numNodesPerElem;

   if ((numRows <= 0) || (numCols <= 0)) {
      fei::console_out() << "SNL_FEI_Structure::allocateConnectivityTable: ERROR, either "
           << "numElems or numNodesPerElem not yet set for blockID "
           << (int)blockID << ". Aborting." << FEI_ENDL;
      MPI_Abort(comm_, -1);
   }

   connTables_[index]->elemNumbers.resize(numRows);

   connTables_[index]->elem_conn_ids = new std::vector<GlobalID>(numRows*numCols);

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
ConnectivityTable& SNL_FEI_Structure::getBlockConnectivity(GlobalID blockID) {

   int index = fei::binarySearch(blockID, blockIDs_);

   if (index < 0) {
      fei::console_out() << "SNL_FEI_Structure::getBlockConnectivity: ERROR, blockID "
           << (int)blockID << " not found. Aborting." << FEI_ENDL;
      MPI_Abort(comm_, -1);
   }  

   return(*(connTables_[index]));
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::finalizeActiveNodes()
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) {
    os << "#  finalizeActiveNodes" << FEI_ENDL;
  }
  //First, make sure that the shared-node IDs are in the nodeDatabase, since
  //they might not have been added yet...

  std::vector<GlobalID>& shNodeIDs = nodeCommMgr_->getSharedNodeIDs();
  for(unsigned i=0; i<shNodeIDs.size(); i++) {
    CHK_ERR( nodeDatabase_->initNodeID(shNodeIDs[i]) );
  }

  if (activeNodesInitialized_) return(0);
  //
  //Now run through the connectivity tables and set the nodal field and block
  //info. Also, replace the nodeID in the element-connectivity lists with an 
  //index from the NodeDatabase object. This will save us a lot of time when
  //looking up nodes later.
  //
  //One other thing we're going to do is give all elements an element-number.
  //This element-number will start at zero locally (meaning that it's not
  //globally unique).
  //
  int elemNumber = 0;
  int numBlocks = blockIDs_.size();
  for(int bIndex=0; bIndex<numBlocks; bIndex++) {

    GlobalID blockID = blocks_[bIndex]->getGlobalBlockID();
    ConnectivityTable& conn = *(connTables_[bIndex]);

    GlobalID* elemConn = NULL;
    NodeDescriptor** elemNodeDescPtrs = NULL;

    int numElems = conn.elemIDs.size();
    if (numElems > 0) {
      elemConn = &((*conn.elem_conn_ids)[0]);
      if (!activeNodesInitialized_) {
        int elemConnLen = conn.elem_conn_ids->size();
        conn.elem_conn_ptrs = new std::vector<NodeDescriptor*>(elemConnLen);
      }
      elemNodeDescPtrs = &((*conn.elem_conn_ptrs)[0]);
    }
    int nodesPerElem = conn.numNodesPerElem;

    int* fieldsPerNodePtr = blocks_[bIndex]->fieldsPerNodePtr();
    int** fieldIDsTablePtr = blocks_[bIndex]->fieldIDsTablePtr();

    //
    //Now we want to go through the connectivity table, and for each node,
    //add its fields and this block to the correct NodeDescriptor from the
    //nodeDatabase_.
    //(Note that the addField and addBlock functions only add the input if
    //it isn't already present in that NodeDescriptor.)
    //
    //We also need to set its owner proc to localProc_ as a default, and 
    //inform the nodeCommMgr that the node appears in local connectivities.
    //

    conn.elemNumbers.resize(numElems);

    int numDistinctFields = blocks_[bIndex]->getNumDistinctFields();
    int fieldID = fieldIDsTablePtr[0][0];

    int elem;
    for(elem=0; elem<numElems; elem++) {
      conn.elemNumbers[elem] = elemNumber++;
      GlobalID* elemNodes = &(elemConn[elem*nodesPerElem]);
      NodeDescriptor** elemNodePtrs = &(elemNodeDescPtrs[elem*nodesPerElem]);

      for(int n=0; n<nodesPerElem; n++) {
        NodeDescriptor* node = NULL;
        int index = nodeDatabase_->getNodeWithID( elemNodes[n], node );
        if (index < 0) {
          fei::console_out() << "ERROR in SNL_FEI_Structure::initializeActiveNodes, "
            << FEI_ENDL << "failed to find node "
            << (int)(elemNodes[n]) << FEI_ENDL;
        }

        if (numDistinctFields == 1) {
          node->addField(fieldID);
        }
        else {
          for(int i=0; i<fieldsPerNodePtr[n]; i++) {
            node->addField(fieldIDsTablePtr[n][i]);
          }
        }

        int blk_idx = getIndexOfBlock(blockID);
        if (blk_idx >= 0) {
          node->addBlockIndex(blk_idx);
        }

        //now store this NodeDescriptor pointer, for fast future lookups
        elemNodePtrs[n] = node;

        node->setOwnerProc(localProc_);
        if (numProcs_ > 1) nodeCommMgr_->informLocal(*node);
      }
    }
  }

  //
  //we also need to run through the penalty constraints and inform the
  //nodeCommMgr that the constrained nodes are local.
  //

  if (numProcs_ > 1) {
    std::map<GlobalID,ConstraintType*>::const_iterator
      cr_iter = getPenConstRecords().begin(),
      cr_end = getPenConstRecords().end();

    while(cr_iter != cr_end) {
      ConstraintType& cr = *((*cr_iter).second);
      std::vector<GlobalID>& nodeID_vec = cr.getMasters();
      GlobalID* nodeIDs = &nodeID_vec[0];
      int numNodes = cr.getMasters().size();

      NodeDescriptor* node = NULL;
      for(int k=0; k<numNodes; ++k) {
	 int err = nodeDatabase_->getNodeWithID(nodeIDs[k], node) ;
         if ( err != -1 ) {
           nodeCommMgr_->informLocal(*node);
         }
	//CHK_ERR( nodeDatabase_->getNodeWithID(nodeIDs[k], node) );
	//CHK_ERR( nodeCommMgr_->informLocal(*node) );
      }
      ++cr_iter;
    }
  }

  activeNodesInitialized_ = true;

  if (debugOutput_) os << "#  leaving finalizeActiveNodes" << FEI_ENDL;
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::finalizeNodeCommMgr()
{
  //call initComplete() on the nodeCommMgr, so that it can
  //finalize shared-node ownerships, etc.

  //request the safetyCheck if the user requested it via the
  //parameters function.
  bool safetyCheck = checkSharedNodes_;

  if (safetyCheck && localProc_==0 && numProcs_>1 && outputLevel_>0) {
    FEI_COUT << "FEI Info: A consistency-check of shared-node data will be "
	 << "performed, which involves all-to-all communication. This check is "
	 << "done only if explicitly requested by parameter "
         << "'FEI_CHECK_SHARED_IDS'."<<FEI_ENDL;
  }

  CHK_ERR( nodeCommMgr_->initComplete(*nodeDatabase_, safetyCheck) );

  if (debugOutput_) {
    FEI_OSTREAM& os = dbgOut();
    int numSharedNodes = nodeCommMgr_->getNumSharedNodes();
    os << "#     numSharedNodes: " << numSharedNodes << FEI_ENDL;
    for(int ns=0; ns<numSharedNodes; ns++) {
      NodeDescriptor& node = nodeCommMgr_->getSharedNodeAtIndex(ns);
      GlobalID nodeID = node.getGlobalNodeID();
      int proc = node.getOwnerProc();
      os << "#     shNodeID " << (int)nodeID << ", owned by proc "<<proc<<FEI_ENDL;
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::setNumNodesAndEqnsPerBlock()
{
  //This function will count how 
  //many active nodes there are for each block.

   int numBlocks = blockIDs_.size();
   std::vector<int> nodesPerBlock(numBlocks);
   std::vector<int> eqnsPerBlock(numBlocks);

   int j;
   for(j=0; j<numBlocks; j++) {
      nodesPerBlock[j] = 0;
      eqnsPerBlock[j] = 0;
   }

   int numNodes = nodeDatabase_->getNumNodeDescriptors();

   for(j=0; j<numNodes; j++) {
     NodeDescriptor* node = NULL;
     nodeDatabase_->getNodeAtIndex(j, node);
     if (node==NULL) continue;

     const int* fieldIDList = node->getFieldIDList();

     int numFields = node->getNumFields();

     const std::vector<unsigned>& blkIndices = node->getBlockIndexList();
     for(std::vector<unsigned>::const_iterator b=blkIndices.begin(), bend=blkIndices.end();
         b!=bend; ++b) {
       nodesPerBlock[*b]++;

       for(int fld=0; fld<numFields; fld++) {
         if (blocks_[*b]->containsField(fieldIDList[fld])) {
           int fSize = getFieldSize(fieldIDList[fld]);
           assert(fSize >= 0);
           eqnsPerBlock[*b] += fSize;
         }
       }
     }
   }

   for(j=0; j<numBlocks; j++) {
     blocks_[j]->setNumActiveNodes(nodesPerBlock[j]);
   }

   //now we need to add the elem-dof to the eqn-count for each block,
   //and then set the total number of eqns on each block.

   for(j=0; j<numBlocks; j++) {
     eqnsPerBlock[j] += blocks_[j]->getNumElemDOFPerElement() *
                        blocks_[j]->getNumElements();

     blocks_[j]->setTotalNumEqns(eqnsPerBlock[j]);
   }
   return(0);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::initializeEqnCommMgr()
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) os << "#  initializeEqnCommMgr" << FEI_ENDL;

  //This function will take information from the (internal member) nodeCommMgr
  //and use it to tell the eqnCommMgr which equations we can expect other
  //processors to contribute to, and also which equations we'll be sending to
  //other processors.
  //
  //This function can not be called until after initComplete() has been called
  //on the nodeCommMgr.
  //
  int numSharedNodes = nodeCommMgr_->getNumSharedNodes();

  for(int i=0; i<numSharedNodes; i++) {
    NodeDescriptor& node = nodeCommMgr_->getSharedNodeAtIndex(i);

    int numFields = node.getNumFields();
    const int* nodeFieldIDList = node.getFieldIDList();
    const int* nodeEqnNums = node.getFieldEqnNumbers();

    int owner = node.getOwnerProc();
    if (owner == localProc_) {
      std::vector<int>& proc = nodeCommMgr_->getSharedNodeProcs(i);
      int numProcs = proc.size();

      for(int p=0; p<numProcs; p++) {

	if (proc[p] == localProc_) continue;

	for(int j=0; j<numFields; j++) {
	  int numEqns = getFieldSize(nodeFieldIDList[j]);
	  assert(numEqns >= 0);

	  int eqn;
	  for(eqn=0; eqn<numEqns; eqn++) {
	    int reducedEqn = -1;
	    bool isSlave = translateToReducedEqn(nodeEqnNums[j]+eqn,
						 reducedEqn);
	    if (isSlave) continue;

	    eqnCommMgr_->addLocalEqn(reducedEqn, proc[p]);
	  }
	}
      }
    }
    else {
      for(int j=0; j<numFields; j++) {
	int numEqns = getFieldSize(nodeFieldIDList[j]);
	assert(numEqns >= 0);
	for(int eqn=0; eqn<numEqns; eqn++) {
	  int reducedEqn = -1;
	  bool isSlave = translateToReducedEqn(nodeEqnNums[j]+eqn, reducedEqn);
	  if (isSlave) continue;

	  eqnCommMgr_->addRemoteIndices(reducedEqn, owner, NULL, 0);
	}
      }
    }
  }

  if (debugOutput_) os << "#  leaving initializeEqnCommMgr" << FEI_ENDL;
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getEqnInfo(int& numGlobalEqns, int& numLocalEqns,
                   int& localStartRow, int& localEndRow){

   numGlobalEqns = numGlobalEqns_;
   numLocalEqns = numLocalEqns_;
   localStartRow = localStartRow_;
   localEndRow = localEndRow_;
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getEqnNumbers(GlobalID ID,
                                     int idType, int fieldID,
                                     int& numEqns, int* eqnNumbers)
{
  //for now, only allow node ID. allowance for element ID is coming soon!!!
  if (idType != FEI_NODE) ERReturn(-1);

  NodeDescriptor* node = NULL;
  CHK_ERR( nodeDatabase_->getNodeWithID(ID, node) );

  numEqns = getFieldSize(fieldID);
  if (numEqns < 0) ERReturn(-1);

  int nodeFieldEqnNumber = -1;
  if (!node->getFieldEqnNumber(fieldID, nodeFieldEqnNumber)) {
    ERReturn(-1);
  }

  for(int i=0; i<numEqns; i++) eqnNumbers[i] = nodeFieldEqnNumber + i;

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getEqnNumbers(int numIDs,
				     const GlobalID* IDs,
				     int idType, int fieldID,
				     int& numEqns, int* eqnNumbers)
{
    //for now, only allow node ID. allowance for element ID is coming soon!!!
    if (idType != FEI_NODE) ERReturn(-1);

    int fieldSize = getFieldSize(fieldID);
    if (fieldSize < 0) {
      ERReturn(-1);
    }

    int offset = 0;
    for(int i=0; i<numIDs; ++i) {
      NodeDescriptor* node = NULL;

      if ( nodeDatabase_->getNodeWithID(IDs[i], node) != 0 ) {
        // fei::console_out() << "SNL_FEI_Structure::getEqnNumbers: ERROR getting node " << IDs[i] << FEI_ENDL;
        for(int ii=0; ii<fieldSize; ii++) {
          eqnNumbers[offset++] = -1;
        }
        continue;
        // ERReturn(-1);
      }

      int nodeFieldEqnNumber = -1;

      if ( !node->getFieldEqnNumber(fieldID, nodeFieldEqnNumber) ) {
        //if we didn't find a field eqn-number, set eqnNumbers to -1. These
        //positions will then be skipped by code that passes the data on down to
        //underlying solver objects.
        for(int ii=0; ii<fieldSize; ii++) {
          eqnNumbers[offset++] = -1;
        }
      }
      else {
        //else we did find valid eqn-numbers...
        for(int ii=0; ii<fieldSize; ii++) {
          eqnNumbers[offset++] = nodeFieldEqnNumber + ii;
        }
      }
    }

    numEqns = fieldSize*numIDs;

    return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::translateToReducedNodeNumber(int nodeNumber, int proc)
{
  if (proc != localProc_) {
    return(nodeNumber - globalNumNodesVanished_[proc]);
  }

  int insertPoint = -1;
  int index =
    fei::binarySearch(nodeNumber, localVanishedNodeNumbers_, insertPoint);

  int localAdjustment = index < 0 ? insertPoint : index + 1;

  return(nodeNumber - localAdjustment - globalNumNodesVanished_[proc]);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getEqnBlkInfo(int& numGlobalEqnBlks,
                                     int& numLocalEqnBlks,
                                     int& localBlkOffset) {

   numGlobalEqnBlks = numGlobalEqnBlks_;
   numLocalEqnBlks = numLocalEqnBlks_;
   localBlkOffset = localBlkOffset_;
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::calculateSlaveEqns(MPI_Comm comm)
{
  FEI_OSTREAM& os = dbgOut();
  if (debugOutput_) os << "#  calculateSlaveEqns" << FEI_ENDL;

  if (eqnCommMgr_ != NULL) delete eqnCommMgr_;
  eqnCommMgr_ = new EqnCommMgr(comm_);
  eqnCommMgr_->setNumRHSs(1);

  std::vector<int> eqns;
  std::vector<int> mEqns;
  std::vector<double> mCoefs;

  for(size_t i=0; i<slaveVars_->size(); i++) {
    int numEqns;
    SlaveVariable* svar = (*slaveVars_)[i];

    eqns.resize( getFieldSize(svar->getFieldID()));
    CHK_ERR( getEqnNumbers(svar->getNodeID(), FEI_NODE, svar->getFieldID(),
			   numEqns, &eqns[0]));

    int slaveEqn = eqns[svar->getFieldOffset()];

    const std::vector<GlobalID>* mNodes = svar->getMasterNodeIDs();
    const std::vector<int>* mFields = svar->getMasterFields();
    const std::vector<double>* mWeights = svar->getWeights();
    const std::vector<double>& mWeightsRef = *mWeights;
    int mwOffset = 0;

    for(size_t j=0; j<mNodes->size(); j++) {
      int mfSize = getFieldSize((*mFields)[j]);

      eqns.resize(mfSize);
      GlobalID mNodeID = (*mNodes)[j];
      int mFieldID = (*mFields)[j];
      CHK_ERR( getEqnNumbers( mNodeID, FEI_NODE, mFieldID,
			      mfSize, &eqns[0]));

      double fei_eps = 1.e-49;

      for(int k=0; k<mfSize; k++) {
	if (std::abs(mWeightsRef[mwOffset++]) > fei_eps) {
	  mEqns.push_back(eqns[k]);
	  mCoefs.push_back(mWeightsRef[mwOffset-1]);
	}
      }
    }

    CHK_ERR( slaveEqns_->addEqn(slaveEqn, &mCoefs[0], &mEqns[0],
		       mEqns.size(), false) );
    mEqns.resize(0);
    mCoefs.resize(0);
  }

#ifndef FEI_SER
  int numLocalSlaves = slaveVars_->size();
  int globalMaxSlaves = 0;
  CHK_ERR( fei::GlobalMax(comm_, numLocalSlaves, globalMaxSlaves) );

  if (globalMaxSlaves > 0) {
    CHK_ERR( gatherSlaveEqns(comm, eqnCommMgr_, slaveEqns_) );
  }
#endif

  globalNumNodesVanished_.resize(numProcs_+1, 0);

  slvEqnNumbers_ = &(slaveEqns_->eqnNumbers());
  numSlvs_ = slvEqnNumbers_->size();
  if (numSlvs_ > 0) {
    //first, let's remove any 'couplings' among the slave equations. A coupling
    //is where a slave depends on a master which is also a slave that depends on
    //other masters.

    if (debugOutput_) {
      os << "#  slave-equations:" << FEI_ENDL;
      os << *slaveEqns_;
      os << "#  leaving calculateSlaveEqns" << FEI_ENDL;
    }

    int levelsOfCoupling;
    CHK_ERR( removeCouplings(*slaveEqns_, levelsOfCoupling) );

    if (debugOutput_) {
      os << "#     SNL_FEI_Structure: " << levelsOfCoupling
	 << " level(s) of coupling removed: " << FEI_ENDL;
    }

    lowestSlv_ = (*slvEqnNumbers_)[0];
    highestSlv_ = (*slvEqnNumbers_)[numSlvs_-1];

    //For each slave equation, we need to find out if we (the local proc) either
    //own or share the node from which that equation arises. If we own or share
    //a slave node, then we will need access to the solution for each of the
    //associated master equations, and all other processors that share the slave
    //will also need access to all of the master solutions.
    //So,
    // 1. for each slave node that we share but don't own, we need to add the
    //   master equations to the eqnCommMgr_ object using addRemoteIndices, if
    //   they're not local.
    // 2. for each slave node that we own, we need to add the master equations
    // to the eqnCommMgr_ object using addLocalEqn for each processor that
    // shares the slave node.

    std::vector<int>& slvEqns = *slvEqnNumbers_;
    std::vector<fei::CSVec*>& mstrEqns = slaveEqns_->eqns();

    //keep track of the number of locally owned nodes that vanish due to the
    //fact that all equations at that node are slave equations.
    int numLocalNodesVanished = 0;

    GlobalID lastNodeID = -1;

    for(int i=0; i<numSlvs_; i++) {
      const NodeDescriptor* node = NULL;
      int reducedSlaveEqn;
      translateToReducedEqn(slvEqns[i], reducedSlaveEqn);
      int numMasters = mstrEqns[i]->size();

      int err = nodeDatabase_->getNodeWithEqn(slvEqns[i], node);
      if (err != 0) {
        if (debugOutput_) {
          os << "#  no local node for slave eqn " << slvEqns[i] << FEI_ENDL;
        }

        continue;
      }

      if (node->getGlobalNodeID() != lastNodeID &&
          node->getOwnerProc() == localProc_) {
        if (nodalEqnsAllSlaves(node, slvEqns)) {
          numLocalNodesVanished++;
          if (fei::sortedListInsert(node->getNodeNumber(), localVanishedNodeNumbers_)
              == -2) {
            ERReturn(-1);
          }
        }
        lastNodeID = node->getGlobalNodeID();
      }

      GlobalID slvNodeID = node->getGlobalNodeID();
      int shIndex = nodeCommMgr_->getSharedNodeIndex(slvNodeID);
      if (shIndex < 0) {
        continue;
      }

      std::vector<int>& sharingProcs = nodeCommMgr_->getSharedNodeProcs(shIndex);

      for(int j=0; j<numMasters; j++) {
        int masterEquation = mstrEqns[i]->indices()[j];
        int owner = getOwnerProcForEqn( masterEquation );

        int reducedMasterEqn;
        translateToReducedEqn(masterEquation, reducedMasterEqn);

        if (owner == localProc_) {
          int numSharing = sharingProcs.size();
          for(int sp=0; sp<numSharing; sp++) {
            if (sharingProcs[sp] == localProc_) continue;

            if (debugOutput_) {
              os << "#     slave node " << (int)slvNodeID << ",eqn "<<slvEqns[i]
                 << ", master eqn " << masterEquation << " eqnCommMgr_->addLocal "
                 << reducedMasterEqn << ", proc " << sharingProcs[sp] << FEI_ENDL;
            }
            eqnCommMgr_->addLocalEqn(reducedMasterEqn, sharingProcs[sp]);
	    slvCommMgr_->addRemoteIndices(reducedSlaveEqn, sharingProcs[sp],
					  &reducedMasterEqn, 1);
          }
        }
        else {
          if (debugOutput_) {
            os << "#     slave node " << (int)slvNodeID << ",eqn "<<slvEqns[i]
               << ", master eqn " << masterEquation << " eqnCommMgr_->addRemote "
               << reducedMasterEqn << ", proc " << owner << FEI_ENDL;
          }
          eqnCommMgr_->addRemoteIndices(reducedMasterEqn, owner,
                                        &reducedMasterEqn, 1);
	  slvCommMgr_->addLocalEqn(reducedSlaveEqn, owner);
        }
      }
    }

    std::vector<int> tmp(1), tmp2(numProcs_);
    tmp[0] = numLocalNodesVanished;
    CHK_ERR( fei::Allgatherv(comm_, tmp, tmp2, globalNumNodesVanished_) );

    if ((int)globalNumNodesVanished_.size() != numProcs_) {
      ERReturn(-1);
    }

    globalNumNodesVanished_.resize(numProcs_+1);
    int tmpNumVanished = 0;
    for(int proc=0; proc<numProcs_; ++proc) {
      int temporary = tmpNumVanished;
      tmpNumVanished += globalNumNodesVanished_[proc];
      globalNumNodesVanished_[proc] = temporary;
    }
    globalNumNodesVanished_[numProcs_] = tmpNumVanished;
  }

  if (slaveMatrix_ != NULL) delete slaveMatrix_;
  slaveMatrix_ = new fei::FillableMat(*slaveEqns_);

  if (debugOutput_) {
    os << "#  slave-equations:" << FEI_ENDL;
    os << *slaveEqns_;
    os << "#  leaving calculateSlaveEqns" << FEI_ENDL;
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::removeCouplings(EqnBuffer& eqnbuf, int& levelsOfCoupling)
{
  std::vector<int>& eqnNumbers = eqnbuf.eqnNumbers();
  std::vector<fei::CSVec*>& eqns = eqnbuf.eqns();

  std::vector<double> tempCoefs;
  std::vector<int> tempIndices;

  levelsOfCoupling = 0;
  bool finished = false;
  while(!finished) {
    bool foundCoupling = false;
    for(size_t i=0; i<eqnNumbers.size(); i++) {
      int rowIndex = eqnbuf.isInIndices(eqnNumbers[i]);

      if (rowIndex == (int)i) {
	fei::console_out() <<" SNL_FEI_Structure::removeCouplings ERROR,"
	     << " illegal master-slave constraint coupling. Eqn "
	     << eqnNumbers[i] << " is both master and slave. " << FEI_ENDL;
	ERReturn(-1);
      }

      while(rowIndex >= 0) {
	foundCoupling = true;
	double coef = 0.0;
	CHK_ERR( eqnbuf.getCoefAndRemoveIndex( eqnNumbers[rowIndex],
					       eqnNumbers[i], coef) );

	std::vector<int>& indicesRef = eqns[i]->indices();
	std::vector<double>& coefsRef = eqns[i]->coefs();

	int len = indicesRef.size();
	tempCoefs.resize(len);
	tempIndices.resize(len);

	double* tempCoefsPtr = &tempCoefs[0];
	int* tempIndicesPtr = &tempIndices[0];
	double* coefsPtr = &coefsRef[0];
	int* indicesPtr = &indicesRef[0];

	for(int j=0; j<len; ++j) {
	  tempIndicesPtr[j] = indicesPtr[j];
	  tempCoefsPtr[j] = coef*coefsPtr[j];
	}

	CHK_ERR( eqnbuf.addEqn(eqnNumbers[rowIndex], tempCoefsPtr,
			       tempIndicesPtr, len, true) );

	rowIndex = eqnbuf.isInIndices(eqnNumbers[i]);
      }
    }
    if (foundCoupling) ++levelsOfCoupling;
    else finished = true;

    if (levelsOfCoupling>1 && !finished) {
      fei::console_out() <<" SNL_FEI_Structure::removeCouplings ERROR,"
	   << " too many (>1) levels of master-slave constraint coupling. "
	   << "Hint: coupling is considered infinite if two slaves depend on "
	   << "each other. This may or may not be the case here." << FEI_ENDL;
      ERReturn(-1);
    }
  }

  return(0);
}

#ifndef FEI_SER
//------------------------------------------------------------------------------
int SNL_FEI_Structure::gatherSlaveEqns(MPI_Comm comm,
				       EqnCommMgr* eqnCommMgr,
				       EqnBuffer* slaveEqns)
{
  int numProcs = 0;
  if (MPI_Comm_size(comm, &numProcs) != MPI_SUCCESS) return(-1);
  if (numProcs == 1) return(0);
  int localProc;
  if (MPI_Comm_rank(comm, &localProc) != MPI_SUCCESS) return(-1);

  //We're going to send all of our slave equations to all other processors, and
  //receive the slave equations from all other processors.
  //So we'll first fill a ProcEqns object with all of our eqn/proc pairs,
  //then use the regular exchange functions from EqnCommMgr.
  ProcEqns localProcEqns, remoteProcEqns;
  std::vector<int>& slvEqnNums = slaveEqns->eqnNumbers();
  fei::CSVec** slvEqnsPtr = NULL;
  if (slaveEqns->eqns().size() > 0) slvEqnsPtr = &(slaveEqns->eqns()[0]);

  for(size_t i=0; i<slvEqnNums.size(); i++) {
    for(int p=0; p<numProcs; p++) {
      if (p == localProc) continue;

      localProcEqns.addEqn(slvEqnNums[i], slvEqnsPtr[i]->size(), p);
    }
  }

  CHK_ERR( eqnCommMgr->mirrorProcEqns(localProcEqns, remoteProcEqns) );

  CHK_ERR( eqnCommMgr->mirrorProcEqnLengths(localProcEqns,
					    remoteProcEqns));

  EqnBuffer remoteEqns;
  CHK_ERR( EqnCommMgr::exchangeEqnBuffers(comm, &localProcEqns, slaveEqns,
					  &remoteProcEqns, &remoteEqns,
					  false) );

  //Now the remoteEqns buffer should hold the slave equations from all other
  //processors, so let's add them to the ones we already have.
  CHK_ERR( slaveEqns->addEqns(remoteEqns, false) );

  return(0);
}
#endif

//------------------------------------------------------------------------------
bool SNL_FEI_Structure::isSlaveEqn(int eqn)
{
  if (numSlvs_ == 0) return(false);

  std::vector<int>& slvEqns = slaveEqns_->eqnNumbers();
  int insertPoint = -1;
  int foundOffset = fei::binarySearch(eqn, slvEqns, insertPoint);

  if (foundOffset >= 0) return(true);
  else return(false);
}

//------------------------------------------------------------------------------
bool SNL_FEI_Structure::translateToReducedEqn(int eqn, int& reducedEqn)
{
  if (numSlvs_ == 0) { reducedEqn = eqn; return(false); }

  if (eqn < lowestSlv_) {reducedEqn = eqn; return(false); }
  if (eqn > highestSlv_) {reducedEqn = eqn - numSlvs_; return(false); }

  int index = 0;
  int foundOffset = fei::binarySearch(eqn, *slvEqnNumbers_, index);

  bool isSlave = false;

  if (foundOffset < 0) {
    reducedEqn = eqn - index;
  }
  else {
    isSlave = true; reducedEqn = eqn - (foundOffset+1);
  }

  return(isSlave);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::translateFromReducedEqn(int reducedEqn)
{
  int numSlvs = slaveEqns_->getNumEqns();

  if (numSlvs == 0) return(reducedEqn);

  const int* slvEqns = &(slaveEqns_->eqnNumbers()[0]);

  if (reducedEqn < slvEqns[0]) return(reducedEqn);

  int eqn = reducedEqn;

  for(int i=0; i<numSlvs; i++) {
    if (eqn < slvEqns[i]) return(eqn);
    eqn++;
  }

  return(eqn);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getMasterEqnNumbers(int slaveEqn,
					  std::vector<int>*& masterEqns)
{
  if (slaveEqns_->getNumEqns() == 0) {
    masterEqns = NULL;
    return(0);
  }

  std::vector<int>& slvEqns = slaveEqns_->eqnNumbers();
  int index = 0;
  int foundOffset = fei::binarySearch(slaveEqn, slvEqns, index);

  if (foundOffset >= 0) {
    masterEqns = &(slaveEqns_->eqns()[foundOffset]->indices());
  }
  else {
    masterEqns = NULL;
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getMasterEqnCoefs(int slaveEqn,
					std::vector<double>*& masterCoefs)
{
  if (slaveEqns_->getNumEqns() == 0) {
    masterCoefs = NULL;
    return(0);
  }

  std::vector<int>& slvEqns = slaveEqns_->eqnNumbers();
  int index = 0;
  int foundOffset = fei::binarySearch(slaveEqn, slvEqns, index);

  if (foundOffset >= 0) {
    masterCoefs = &(slaveEqns_->eqns()[foundOffset]->coefs());
  }
  else {
    masterCoefs = NULL;
  }

  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getMasterEqnRHS(int slaveEqn,
				      double& rhsValue)
{
  if (slaveEqns_->getNumEqns() == 0) {
    return(0);
  }

  std::vector<int>& slvEqns = slaveEqns_->eqnNumbers();
  int index = 0;
  int foundOffset = fei::binarySearch(slaveEqn, slvEqns, index);

  if (foundOffset >= 0) {
    std::vector<double>* rhsCoefsPtr = (*(slaveEqns_->rhsCoefsPtr()))[foundOffset];
    rhsValue = (*rhsCoefsPtr)[0];
  }

  return(0);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getScatterIndices_ID(GlobalID blockID, GlobalID elemID, 
                                            int interleaveStrategy,
                                            int* scatterIndices)
{
   int index = fei::binarySearch(blockID, blockIDs_);

   if (index < 0) {
      fei::console_out() << "SNL_FEI_Structure::getScatterIndices_ID: ERROR, blockID "
           << (int)blockID << " not found. Aborting." << FEI_ENDL;
      std::abort();
   }

   std::map<GlobalID,int>& elemIDs = connTables_[index]->elemIDs;

   std::map<GlobalID,int>::const_iterator
     iter = elemIDs.find(elemID);

   if (iter == elemIDs.end()) {
      fei::console_out() << "SNL_FEI_Structure::getScatterIndices_ID: ERROR, blockID: " 
           << (int)blockID << ", elemID "
           << (int)elemID << " not found. Aborting." << FEI_ENDL;
      std::abort();
   }

   int elemIndex = iter->second;

   getScatterIndices_index(index, elemIndex,
                           interleaveStrategy, scatterIndices);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getScatterIndices_ID(GlobalID blockID, GlobalID elemID, 
					     int interleaveStrategy,
					     int* scatterIndices,
					     int* blkScatterIndices,
					     int* blkSizes)
{
   int index = fei::binarySearch(blockID, blockIDs_);

   if (index < 0) {
      fei::console_out() << "SNL_FEI_Structure::getScatterIndices_ID: ERROR, blockID "
           << (int)blockID << " not found. Aborting." << FEI_ENDL;
      std::abort();
   }

   std::map<GlobalID,int>& elemIDs = connTables_[index]->elemIDs;

   std::map<GlobalID,int>::const_iterator
     iter = elemIDs.find(elemID);

   if (iter == elemIDs.end()) {
      fei::console_out() << "SNL_FEI_Structure::getScatterIndices_ID: ERROR, blockID: " 
           << (int)blockID << ", elemID "
           << (int)elemID << " not found. Aborting." << FEI_ENDL;
      std::abort();
   }

   int elemIndex = iter->second;

   getScatterIndices_index(index, elemIndex,
                           interleaveStrategy, scatterIndices,
			   blkScatterIndices, blkSizes);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getBlkScatterIndices_index(int blockIndex,
						  int elemIndex,
						  int* scatterIndices)
{
  BlockDescriptor& block = *(blocks_[blockIndex]);
  int numNodes = block.getNumNodesPerElement();
  work_nodePtrs_.resize(numNodes);
  NodeDescriptor** nodes = &work_nodePtrs_[0];
  int err = getElemNodeDescriptors(blockIndex, elemIndex, nodes);
  if (err) {
    fei::console_out() << "SNL_FEI_Structure::getBlkScatterIndices_index: ERROR getting"
	 << " node descriptors." << FEI_ENDL;
    ERReturn(-1);
  }

  int offset = 0;
  return( getNodeBlkIndices(nodes, numNodes, scatterIndices, offset) );
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getScatterIndices_index(int blockIndex, int elemIndex,
						int interleaveStrategy,
						int* scatterIndices)
{
//On input, scatterIndices, is assumed to be allocated by the calling code,
// and be of length the number of equations per element.
//
   BlockDescriptor& block = *(blocks_[blockIndex]);
   int numNodes = block.getNumNodesPerElement();
   int* fieldsPerNode = block.fieldsPerNodePtr();
   int** fieldIDs = block.fieldIDsTablePtr();

   work_nodePtrs_.resize(numNodes);
   NodeDescriptor** nodes = &work_nodePtrs_[0];

   int err = getElemNodeDescriptors(blockIndex, elemIndex, nodes);
   if (err) {
      fei::console_out() << "SNL_FEI_Structure::getScatterIndices_index: ERROR getting"
           << " node descriptors." << FEI_ENDL;
      std::abort();
   }

   int offset = 0;
   if (fieldDatabase_->size() == 1) {
     err = getNodeIndices_simple(nodes, numNodes, fieldIDs[0][0],
				    scatterIndices, offset);
     if (err) fei::console_out() << "ERROR in getNodeIndices_simple." << FEI_ENDL;
   }
   else {
     switch (interleaveStrategy) {
     case 0:
       err = getNodeMajorIndices(nodes, numNodes, fieldIDs, fieldsPerNode,
				 scatterIndices, offset);
       if (err) fei::console_out() << "ERROR in getNodeMajorIndices." << FEI_ENDL;
       break;

     case 1:
       err = getFieldMajorIndices(nodes, numNodes, fieldIDs, fieldsPerNode,
				  scatterIndices, offset);
       if (err) fei::console_out() << "ERROR in getFieldMajorIndices." << FEI_ENDL;
       break;

     default:
       fei::console_out() << "ERROR, unrecognized interleaveStrategy." << FEI_ENDL;
       break;
     }
   }

   //now the element-DOF.
   int numElemDOF = blocks_[blockIndex]->getNumElemDOFPerElement();
   std::vector<int>& elemDOFEqns = blocks_[blockIndex]->elemDOFEqnNumbers();

   for(int i=0; i<numElemDOF; i++) {
      scatterIndices[offset++] = elemDOFEqns[elemIndex] + i;
   }
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::getScatterIndices_index(int blockIndex, int elemIndex,
						int interleaveStrategy,
						int* scatterIndices,
						int* blkScatterIndices,
						int* blkSizes)
{
//On input, scatterIndices, is assumed to be allocated by the calling code,
// and be of length the number of equations per element.
//
   BlockDescriptor& block = *(blocks_[blockIndex]);
   int numNodes = block.getNumNodesPerElement();
   int* fieldsPerNode = block.fieldsPerNodePtr();
   int** fieldIDs = block.fieldIDsTablePtr();

   work_nodePtrs_.resize(numNodes);
   NodeDescriptor** nodes = &work_nodePtrs_[0];

   int err = getElemNodeDescriptors(blockIndex, elemIndex, nodes);
   if (err) {
      fei::console_out() << "SNL_FEI_Structure::getScatterIndices_index: ERROR getting"
           << " node descriptors." << FEI_ENDL;
      std::abort();
   }

   int offset = 0, blkOffset = 0;
   if (fieldDatabase_->size() == 1) {
     err = getNodeIndices_simple(nodes, numNodes, fieldIDs[0][0],
				 scatterIndices, offset,
				 blkScatterIndices, blkSizes, blkOffset);
     if (err) fei::console_out() << "ERROR in getNodeIndices_simple." << FEI_ENDL;
   }
   else {
     switch (interleaveStrategy) {
     case 0:
       err = getNodeMajorIndices(nodes, numNodes, fieldIDs, fieldsPerNode,
				 scatterIndices, offset,
				 blkScatterIndices, blkSizes, blkOffset);
       if (err) fei::console_out() << "ERROR in getNodeMajorIndices." << FEI_ENDL;
       break;

     case 1:
       err = getFieldMajorIndices(nodes, numNodes, fieldIDs, fieldsPerNode,
				  scatterIndices, offset);
       if (err) fei::console_out() << "ERROR in getFieldMajorIndices." << FEI_ENDL;
       break;

     default:
       fei::console_out() << "ERROR, unrecognized interleaveStrategy." << FEI_ENDL;
       break;
     }
   }

   //now the element-DOF.
   int numElemDOF = blocks_[blockIndex]->getNumElemDOFPerElement();
   std::vector<int>& elemDOFEqns = blocks_[blockIndex]->elemDOFEqnNumbers();

   for(int i=0; i<numElemDOF; i++) {
      scatterIndices[offset++] = elemDOFEqns[elemIndex] + i;
      if (interleaveStrategy == 0) {
	blkSizes[blkOffset] = 1;
	blkScatterIndices[blkOffset++] = elemDOFEqns[elemIndex] + i;
      }
   }
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getElemNodeDescriptors(int blockIndex, int elemIndex,
                                             NodeDescriptor** nodes)
{
  //Private function, called by 'getScatterIndices_index'. We can safely
  //assume that blockIndex and elemIndex are valid in-range indices.
  //
  //This function's task is to obtain the NodeDescriptor objects, from the
  //nodeDatabase, that are connected to the specified element.
  //
  int numNodes = connTables_[blockIndex]->numNodesPerElem;

  if (activeNodesInitialized_) {
    NodeDescriptor** elemNodePtrs =
      &((*connTables_[blockIndex]->elem_conn_ptrs)[elemIndex*numNodes]);
    for(int i=0; i<numNodes; ++i) nodes[i] = elemNodePtrs[i];
  }
  else {
    const GlobalID* elemConn = 
      &((*connTables_[blockIndex]->elem_conn_ids)[elemIndex*numNodes]);
    for(int i=0; i<numNodes; ++i) {
      CHK_ERR( nodeDatabase_->getNodeWithID(elemConn[i], nodes[i]));
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getNodeIndices_simple(NodeDescriptor** nodes,
					     int numNodes,
					     int fieldID,
					     int* scatterIndices,
					     int& offset)
{
  int fieldSize = getFieldSize(fieldID);

  for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
    NodeDescriptor& node = *(nodes[nodeIndex]);
    const int* eqnNumbers = node.getFieldEqnNumbers();
    int eqn = eqnNumbers[0];
    scatterIndices[offset++] = eqn;
    if (fieldSize > 1) {
      for(int i=1; i<fieldSize; i++) {
	scatterIndices[offset++] = eqn+i;
      }
    }
  }
  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getNodeIndices_simple(NodeDescriptor** nodes,
					     int numNodes,
					     int fieldID,
					     int* scatterIndices,
					     int& offset,
					     int* blkScatterIndices,
					     int* blkSizes,
					     int& blkOffset)
{
  int fieldSize = getFieldSize(fieldID);

  for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
    NodeDescriptor& node = *(nodes[nodeIndex]);
    const int* eqnNumbers = node.getFieldEqnNumbers();
    int eqn = eqnNumbers[0];
    scatterIndices[offset++] = eqn;
    if (fieldSize > 1) {
      for(int i=1; i<fieldSize; i++) {
	scatterIndices[offset++] = eqn+i;
      }
    }
    blkSizes[blkOffset] = node.getNumNodalDOF();
    blkScatterIndices[blkOffset++] = node.getBlkEqnNumber();
  }
  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getNodeMajorIndices(NodeDescriptor** nodes, int numNodes,
                                          int** fieldIDs, int* fieldsPerNode,
                                          int* scatterIndices, int& offset)
{
  for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {

      NodeDescriptor& node = *(nodes[nodeIndex]);

      const int* nodeFieldIDList = node.getFieldIDList();
      const int* nodeEqnNums = node.getFieldEqnNumbers();
      int numFields = node.getNumFields();

      int* fieldID_ind = fieldIDs[nodeIndex];

      for(int j=0; j<fieldsPerNode[nodeIndex]; j++) {
         int numEqns = getFieldSize(fieldID_ind[j]);
         assert(numEqns >= 0);

         int insert = -1;
         int nind = fei::binarySearch(fieldID_ind[j], nodeFieldIDList,
                                             numFields, insert);

         if (nind >= 0) {
	   int eqn = nodeEqnNums[nind];

	   if (eqn < 0) {
	     FEI_COUT << "SNL_FEI_Structure::getNodeMajorIndices: ERROR, node "
		  << (int)node.getGlobalNodeID()
		  << ", field " << nodeFieldIDList[nind]
		  << " has equation number " << eqn << FEI_ENDL;
	     std::abort();
	   }

	   for(int jj=0; jj<numEqns; jj++) {
	     scatterIndices[offset++] = eqn+jj;
	   }
         }
         else {
	   if (outputLevel_ > 2) {
	     fei::console_out() << "WARNING, field ID " << fieldIDs[nodeIndex][j]
		  << " not found for node "
		  << (int)(node.getGlobalNodeID()) << FEI_ENDL;
	   }
         }
      }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getNodeBlkIndices(NodeDescriptor** nodes,
					 int numNodes,
					 int* scatterIndices,
					 int& offset)
{
  for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
    NodeDescriptor* node = nodes[nodeIndex];
    scatterIndices[offset++] = node->getBlkEqnNumber();
  }
  return(0);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getNodeMajorIndices(NodeDescriptor** nodes, int numNodes,
					   int** fieldIDs, int* fieldsPerNode,
					   int* scatterIndices, int& offset,
					   int* blkScatterIndices,
					   int* blkSizes,
					   int& blkOffset)
{
  for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {

      NodeDescriptor& node = *(nodes[nodeIndex]);

      const int* nodeFieldIDList = node.getFieldIDList();
      const int* nodeEqnNums = node.getFieldEqnNumbers();
      int numFields = node.getNumFields();

      blkSizes[blkOffset] = node.getNumNodalDOF();
      blkScatterIndices[blkOffset++] = node.getBlkEqnNumber();

      int* fieldID_ind = fieldIDs[nodeIndex];

      for(int j=0; j<fieldsPerNode[nodeIndex]; j++) {
         int numEqns = getFieldSize(fieldID_ind[j]);
         assert(numEqns >= 0);

         int insert = -1;
         int nind = fei::binarySearch(fieldID_ind[j], nodeFieldIDList,
                                             numFields, insert);

         if (nind >= 0) {
	   int eqn = nodeEqnNums[nind];
	   if (eqn < 0) {
	     FEI_COUT << "SNL_FEI_Structure::getNodeMajorIndices: ERROR, node "
		  << (int)node.getGlobalNodeID()
		  << ", field " << nodeFieldIDList[nind]
		  << " has equation number " << eqn << FEI_ENDL;
	     std::abort();
	   }

	   for(int jj=0; jj<numEqns; jj++) {
	     scatterIndices[offset++] = eqn+jj;
	   }
         }
         else {
	   if (outputLevel_ > 2) {
	     fei::console_out() << "WARNING, field ID " << fieldIDs[nodeIndex][j]
		  << " not found for node "
		  << (int)(node.getGlobalNodeID()) << FEI_ENDL;
	   }
         }
      }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getNodeMajorIndices(NodeDescriptor** nodes, int numNodes,
                                          std::vector<int>* fieldIDs,
                                          std::vector<int>& fieldsPerNode,
                                          std::vector<int>& scatterIndices)
{
   int offset = 0;
   scatterIndices.resize(0);

   for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {

      NodeDescriptor& node = *(nodes[nodeIndex]);

      const int* nodeFieldIDList = node.getFieldIDList();
      const int* nodeEqnNums = node.getFieldEqnNumbers();
      int numFields = node.getNumFields();

      int* fieldID_ind = &(fieldIDs[nodeIndex][0]);

      for(int j=0; j<fieldsPerNode[nodeIndex]; j++) {
         int numEqns = getFieldSize(fieldID_ind[j]);
         assert(numEqns >= 0);

         int insert = -1;
         int nind = fei::binarySearch(fieldID_ind[j], nodeFieldIDList,
                                             numFields, insert);

         if (nind >= 0) {
	   int eqn = nodeEqnNums[nind];

	   if (eqn < 0) {
	     FEI_COUT << "SNL_FEI_Structure::getNodeMajorIndices: ERROR, node "
		  << (int)node.getGlobalNodeID()
		  << ", field " << nodeFieldIDList[nind]
		  << " has equation number " << eqn << FEI_ENDL;
	     MPI_Abort(comm_, -1);
	   }

	   scatterIndices.resize(offset+numEqns);
	   int* inds = &scatterIndices[0];

	   for(int jj=0; jj<numEqns; jj++) {
	     inds[offset+jj] = eqn+jj;
	   }
	   offset += numEqns;
         }
         else {
	   if (outputLevel_ > 2) {
	     fei::console_out() << "WARNING, field ID " << fieldID_ind[j]
		  << " not found for node "
		  << (int)node.getGlobalNodeID() << FEI_ENDL;
           }
         }
      }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getFieldMajorIndices(NodeDescriptor** nodes, int numNodes,
                                          int** fieldIDs, int* fieldsPerNode,
                                          int* scatterIndices, int& offset)
{
  //In this function we want to group equations by field, but
  //in what order should the fields be?
  //Let's just run through the fieldIDs table, and add the fields to a
  //flat list, in the order we encounter them, but making sure no fieldID
  //gets added more than once.

  std::vector<int> fields;
  for(int ii=0; ii<numNodes; ii++) {
    for(int j=0; j<fieldsPerNode[ii]; j++) {
      if (std::find(fields.begin(), fields.end(), fieldIDs[ii][j]) == fields.end()) {
        fields.push_back(fieldIDs[ii][j]);
      }
    }
  }

  int* fieldsPtr = fields.size()>0 ? &fields[0] : NULL;

  //ok, we've got our flat list of fields, so let's proceed to get the
  //scatter indices.

  for(size_t i=0; i<fields.size(); i++) {
    int field = fieldsPtr[i];

    for(int nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex) {
      int fidx = fei::searchList(field, fieldIDs[nodeIndex],
				    fieldsPerNode[nodeIndex]);
      if (fidx < 0) {
	continue;
      }

      NodeDescriptor* node = nodes[nodeIndex];

      const int* nodeFieldIDList = node->getFieldIDList();
      const int* nodeEqnNums = node->getFieldEqnNumbers();
      int numFields = node->getNumFields();

      int numEqns = getFieldSize(field);
      assert(numEqns >= 0);

      int insert = -1;
      int nind = fei::binarySearch(field, nodeFieldIDList,
					  numFields, insert);

      if (nind > -1) {
	for(int jj=0; jj<numEqns; ++jj) {
	  scatterIndices[offset++] = nodeEqnNums[nind]+jj;
	}
      }
      else {
	ERReturn(-1);
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int SNL_FEI_Structure::getFieldMajorIndices(NodeDescriptor** nodes, int numNodes,
                                          std::vector<int>* fieldIDs,
                                          std::vector<int>& fieldsPerNode,
                                          std::vector<int>& scatterIndices)
{
   //In this function we want to group equations by field, but
   //in what order should the fields be?
   //Let's just run through the fieldIDs table, and add the fields to a
   //flat list, in the order we encounter them, but making sure no fieldID
   //gets added more than once.

   std::vector<int> fields;
   for(int ii=0; ii<numNodes; ii++) {
      for(int j=0; j<fieldsPerNode[ii]; j++) {
         std::vector<int>& fldIDs = fieldIDs[ii];
         if (std::find(fields.begin(), fields.end(), fldIDs[j]) == fields.end()) {
           fields.push_back(fldIDs[j]);
         }
      }
   }

   //ok, we've got our flat list of fields, so let's proceed to get the
   //scatter indices.

   int offset = 0;
   scatterIndices.resize(0);

   for(size_t i=0; i<fields.size(); i++) {
      for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {

         const int* nodeFieldIDList = nodes[nodeIndex]->getFieldIDList();
         const int* nodeEqnNums = nodes[nodeIndex]->getFieldEqnNumbers();
         int numFields = nodes[nodeIndex]->getNumFields();

         int numEqns = getFieldSize(fields[i]);
         assert(numEqns >= 0);

         int insert = -1;
         int nind = fei::binarySearch(fields[i], nodeFieldIDList,
                                             numFields, insert);

         if (nind >= 0) {
            for(int jj=0; jj<numEqns; jj++) {
               scatterIndices.push_back(nodeEqnNums[nind]+jj);
	       offset++;
            }
         }
         else {
	   if (outputLevel_ > 2) {
	     fei::console_out() << "WARNING, field ID " << fields[i]
		  << " not found for node "
		  << (int)nodes[nodeIndex]->getGlobalNodeID() << FEI_ENDL;
	   }
         }
      }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
void SNL_FEI_Structure::addCR(int CRID,
			     snl_fei::Constraint<GlobalID>*& cr,
	       std::map<GlobalID,snl_fei::Constraint<GlobalID>* >& crDB)
{
  std::map<GlobalID,snl_fei::Constraint<GlobalID>* >::iterator
    cr_iter = crDB.find(CRID);

  if (cr_iter == crDB.end()) {
    cr = new ConstraintType;
    crDB.insert(std::pair<GlobalID,snl_fei::Constraint<GlobalID>*>(CRID, cr));
  }
}

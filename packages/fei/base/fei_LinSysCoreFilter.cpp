/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_sstream.hpp"

#include <limits>
#include <cmath>
#include <assert.h>

#include "fei_utils.hpp"
#include <fei_impl_utils.hpp>

#include "fei_defs.h"

#include <fei_ostream_ops.hpp>
#include "fei_CommUtils.hpp"
#include "fei_TemplateUtils.hpp"
#include "snl_fei_Constraint.hpp"
typedef snl_fei::Constraint<GlobalID> ConstraintType;

#include "fei_DirichletBCManager.hpp"
#include "fei_FillableMat.hpp"
#include "fei_CSVec.hpp"
#include "FEI_Implementation.hpp"
#include "fei_EqnCommMgr.hpp"
#include "fei_ConnectivityTable.hpp"
#include "fei_NodeDescriptor.hpp"
#include "fei_NodeDatabase.hpp"
#include "fei_BlockDescriptor.hpp"
#include "SNL_FEI_Structure.hpp"
#include "snl_fei_Utils.hpp"
#include "fei_Data.hpp"
#include "fei_LinearSystemCore.hpp"
#include "fei_LinSysCore_flexible.hpp"

#include "fei_LinSysCoreFilter.hpp"

#undef fei_file
#define fei_file "fei_LinSysCoreFilter.cpp"
#include "fei_ErrMacros.hpp"

#define ASSEMBLE_PUT 0
#define ASSEMBLE_SUM 1


//------------------------------------------------------------------------------
LinSysCoreFilter::LinSysCoreFilter(FEI_Implementation* owner,
                                   MPI_Comm comm,
                                   SNL_FEI_Structure* probStruct,
                                   LinearSystemCore* lsc,
                                   int masterRank)
 : Filter(probStruct),
   timesInitializeCalled_(0),
   lsc_(lsc),
   useLookup_(true),
   newMatrixData_(false),
   newVectorData_(false),
   newConstraintData_(false),
   newBCData_(false),
   connectivitiesInitialized_(false),
   firstRemEqnExchange_(true),
   needToCallMatrixLoadComplete_(false),
   resolveConflictRequested_(false),
   localStartRow_(0),             //
   localEndRow_(0),               //Initialize all private variables here,
   numLocalEqns_(0),              //in the order that they are declared.
   numGlobalEqns_(0),             //(Don't bother initializing those that will
   numLocallyOwnedNodes_(0),      //be initialized in the body of the 
   numGlobalNodes_(0),            //constructor below.)
   firstLocalNodeNumber_(-1),     //
   blockMatrix_(false),
   tooLateToChooseBlock_(false),
   numLocalEqnBlks_(0),
   localReducedBlkOffset_(0),
   numLocalReducedEqnBlks_(0),
   iterations_(0),
   currentRHS_(0),
   rhsIDs_(),
   outputLevel_(0),
   comm_(comm),
   masterRank_(masterRank),
   problemStructure_(probStruct),
   matrixAllocated_(false),
   rowIndices_(),
   rowColOffsets_(0),
   colIndices_(0),
   nodeIDType_(0),
   eqnCommMgr_(NULL),
   eqnCommMgr_put_(NULL),
   maxElemRows_(0),
   scatterIndices_(),
   blkScatterIndices_(),
   iworkSpace_(),
   iworkSpace2_(),
   dworkSpace_(),
   dworkSpace2_(),
   eStiff_(NULL),
   eStiff1D_(NULL),
   eLoad_(NULL)
{
    localRank_ = fei::localProc(comm_);
    numProcs_ = fei::numProcs(comm_);

    internalFei_ = 0;

    numRHSs_ = 1;
    rhsIDs_.resize(numRHSs_);
    rhsIDs_[0] = 0;

    bcManager_ = new fei::DirichletBCManager(probStruct);
    eqnCommMgr_ = problemStructure_->getEqnCommMgr().deepCopy();
    int err = createEqnCommMgr_put();

   //We need to get the parameters from the owning FEI_Implementation...
   int numParams = -1;
   char** paramStrings = NULL;
   err = owner->getParameters(numParams, paramStrings);

   //Now let's pass them into our own parameter-handling mechanism.
   err = parameters(numParams, paramStrings);
   if (err != 0) {
     fei::console_out() << "LinSysCoreFilter::LinSysCoreFilter ERROR, parameters failed." << FEI_ENDL;
     MPI_Abort(comm_, -1);
   }

   Kid_ = new fei::FillableMat;
   Kdi_ = new fei::FillableMat;
   Kdd_ = new fei::FillableMat;
   reducedEqnCounter_ = 0;
   reducedRHSCounter_ = 0;

   return;
}

//------------------------------------------------------------------------------
LinSysCoreFilter::~LinSysCoreFilter() {
//
//  Destructor function. Free allocated memory, etc.
//
  numRHSs_ = 0;

  delete bcManager_;
  delete eqnCommMgr_;
  delete eqnCommMgr_put_;

  delete [] eStiff_;
  delete [] eStiff1D_;
  delete [] eLoad_;

  delete Kid_;
  delete Kdi_;
  delete Kdd_;
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::initialize()
{
// Determine final sparsity pattern for setting the structure of the
// underlying sparse matrix.
//
  debugOutput("#LinSysCoreFilter::initialize");

  numLocalEqns_ = problemStructure_->getNumLocalEqns();
  numLocalEqnBlks_ = problemStructure_->getNumLocalEqnBlks();
  numLocalReducedEqnBlks_ = problemStructure_->getNumLocalReducedEqnBlks();

  // now, obtain the global equation info, such as how many equations there
  // are globally, and what the local starting and ending row-numbers are.

  // let's also get the number of global nodes, and a first-local-node-number.
  // node-number is a globally 0-based number we are assigning to nodes.
  // node-numbers are contiguous on a processor -- i.e., a processor owns a
  // contiguous block of node-numbers. This provides an easier-to-work-with
  // node numbering than the application-supplied nodeIDs which may not be
  // assumed to be contiguous or 0-based, or anything else.

  std::vector<int>& eqnOffsets = problemStructure_->getGlobalEqnOffsets();
  localStartRow_ = eqnOffsets[localRank_];
  localEndRow_ = eqnOffsets[localRank_+1]-1;
  numGlobalEqns_ = eqnOffsets[numProcs_];

  std::vector<int>& nodeOffsets = problemStructure_->getGlobalNodeOffsets();
  firstLocalNodeNumber_ =  nodeOffsets[localRank_];
  numGlobalNodes_ = nodeOffsets[numProcs_];

  //--------------------------------------------------------------------------
  //  ----- end active equation calculations -----

  if (eqnCommMgr_ != NULL) delete eqnCommMgr_;
  eqnCommMgr_ = NULL;
  if (eqnCommMgr_put_ != NULL) delete eqnCommMgr_put_;
  eqnCommMgr_put_ = NULL;

  eqnCommMgr_ = problemStructure_->getEqnCommMgr().deepCopy();
  if (eqnCommMgr_ == NULL) ERReturn(-1);

  int err = createEqnCommMgr_put();
  if (err != 0) ERReturn(err);

  //(we need to set the number of RHSs in the eqn comm manager)
  eqnCommMgr_->setNumRHSs(numRHSs_);

  //let's let the underlying linear system know what the global offsets are.
  //While we're dealing with global offsets, we'll also calculate the starting
  //and ending 'reduced' rows, etc.
  CHK_ERR( initLinSysCore() );

  setLinSysCoreCREqns();

  if (timesInitializeCalled_ == 0) {
    //
    // let's prepare some arrays for handing the matrix structure to
    // the linear system.

    std::vector<int> rowLengths;
    CHK_ERR( problemStructure_->getMatrixRowLengths(rowLengths) );

    int numReducedEqns = problemStructure_->getNumLocalReducedEqns();
    int maxBlkSize = problemStructure_->getGlobalMaxBlkSize();
    std::vector<int> blkSizes(numLocalReducedEqnBlks_, 1);

    int numNonzeros = 0;
    for(size_t ii=0; ii<rowLengths.size(); ++ii) {
      numNonzeros += rowLengths[ii];
    }

    std::vector<int> indices_1D(numNonzeros);
    std::vector<int*> indices(numReducedEqns);

    int offset = 0;
    for(size_t ii=0; ii<rowLengths.size(); ++ii) {
      indices[ii] = &(indices_1D[offset]);
      offset += rowLengths[ii];
    }

    if (maxBlkSize == 0) ERReturn(-1);

    if (maxBlkSize == 1) {
      CHK_ERR( problemStructure_->getMatrixStructure(&indices[0], rowLengths) );

      debugOutput("#LinSysCoreFilter calling point lsc_->setMatrixStructure");
      CHK_ERR( lsc_->setMatrixStructure(&indices[0], &rowLengths[0],
                                        &indices[0], &rowLengths[0], &blkSizes[0]) );
    }
    else {
      std::vector<int> blkRowLengths;
      int* blkIndices_1D = numNonzeros > 0 ? new int[numNonzeros] : NULL;
      int** blkIndices = numLocalReducedEqnBlks_ ?
        new int*[numLocalReducedEqnBlks_] : NULL;
      if (blkIndices == NULL && numLocalReducedEqnBlks_ != 0) ERReturn(-1);

      CHK_ERR( problemStructure_->getMatrixStructure(&indices[0], rowLengths,
                                                     blkIndices, blkIndices_1D,
                                                     blkRowLengths, blkSizes) );

      offset = 0;
      for(int ii=0; ii<numLocalReducedEqnBlks_; ++ii) {
        blkIndices[ii] = &(blkIndices_1D[offset]);
        offset += blkRowLengths[ii];
      }

      debugOutput("#LinSysCoreFilter calling block lsc_->setMatrixStructure");
      CHK_ERR( lsc_->setMatrixStructure(&indices[0], &rowLengths[0],
                                        blkIndices, &blkRowLengths[0], &blkSizes[0]) );

      if (numLocalReducedEqnBlks_ != 0) delete [] blkIndices;
      if (numNonzeros > 0) delete [] blkIndices_1D;
    }
  }

  matrixAllocated_ = true;

  debugOutput("#leaving LinSysCoreFilter::initialize");
  ++timesInitializeCalled_;
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::createEqnCommMgr_put()
{
  if (eqnCommMgr_put_ != NULL) return(0);

  eqnCommMgr_put_  = eqnCommMgr_->deepCopy();
  if (eqnCommMgr_put_ == NULL) ERReturn(-1);

  eqnCommMgr_put_->resetCoefs();
  eqnCommMgr_put_->accumulate_ = false;
  return(0);
}

//==============================================================================
int LinSysCoreFilter::initLinSysCore()
{
  debugOutput("#LinSysCoreFilter calling lsc_->setLookup");
  int err = lsc_->setLookup(*problemStructure_);

  if (err != 0) {
    useLookup_ = false;
  }

  std::vector<int>& globalNodeOffsets = problemStructure_->getGlobalNodeOffsets();
  std::vector<int>& globalEqnOffsets = problemStructure_->getGlobalEqnOffsets();
  std::vector<int>& globalBlkEqnOffsets =
    problemStructure_->getGlobalBlkEqnOffsets();

  int startRow = localStartRow_;
  while(problemStructure_->isSlaveEqn(startRow)) startRow++;
  int endRow = localEndRow_;
  while(problemStructure_->isSlaveEqn(endRow)) endRow--;

  problemStructure_->translateToReducedEqn(startRow, reducedStartRow_);
  problemStructure_->translateToReducedEqn(endRow, reducedEndRow_);
  numReducedRows_ = reducedEndRow_ - reducedStartRow_ + 1;

  std::vector<int> reducedEqnOffsets(globalEqnOffsets.size());
  std::vector<int> reducedBlkEqnOffsets(globalBlkEqnOffsets.size());
  std::vector<int> reducedNodeOffsets(globalNodeOffsets.size());

  int numSlaveEqns = problemStructure_->numSlaveEquations();

  int reducedNodeNum =
    problemStructure_->getAssociatedNodeNumber(reducedStartRow_);

  std::vector<int> tmpSend(2), tmpRecv, tmpRecvLengths;
  tmpSend[0] = reducedStartRow_;
  tmpSend[1] = reducedNodeNum;
  CHK_ERR( fei::Allgatherv(comm_, tmpSend, tmpRecvLengths, tmpRecv) );

  for(size_t ii=0; ii<globalEqnOffsets.size()-1; ++ii) {
    reducedEqnOffsets[ii] = tmpRecv[ii*2];
    reducedNodeOffsets[ii] = tmpRecv[ii*2+1];

    //Major assumption: we're assuming that if there are slave equations, then
    //there are no lagrange-multiplier constraints. Because if there are no
    //lagrange-multiplier constraints, then blkEqn == nodeNum.
    if (numSlaveEqns > 0) {
      reducedBlkEqnOffsets[ii] = reducedNodeOffsets[ii];
    }
    else {
      reducedBlkEqnOffsets[ii] = globalBlkEqnOffsets[ii];
    }
  }

  if (localRank_ == numProcs_-1) {
    reducedNodeNum = problemStructure_->
      translateToReducedNodeNumber(globalNodeOffsets[numProcs_]-1, localRank_);
    int blkEqn = globalBlkEqnOffsets[numProcs_]-1;
    if (numSlaveEqns > 0) {
      blkEqn = reducedNodeNum;
    }

    tmpSend.resize(3);
    tmpSend[0] = reducedEndRow_;
    tmpSend[1] = reducedNodeNum;
    tmpSend[2] = blkEqn;
  }
  else {
    tmpSend.resize(3);
  }

  CHK_ERR( fei::Bcast(comm_, tmpSend, numProcs_-1) );
  reducedEqnOffsets[numProcs_] = tmpSend[0]+1;
  reducedNodeOffsets[numProcs_] = tmpSend[1]+1;
  reducedBlkEqnOffsets[numProcs_] = tmpSend[2]+1;

  debugOutput("#LinSysCoreFilter calling lsc_->setGlobalOffsets");
  CHK_ERR( lsc_->setGlobalOffsets(numProcs_+1,
                                  &reducedNodeOffsets[0],
                                  &reducedEqnOffsets[0],
                                  &reducedBlkEqnOffsets[0]) );

  if (connectivitiesInitialized_) return(0);

  int numBlocks = problemStructure_->getNumElemBlocks();
  NodeDatabase& nodeDB     = problemStructure_->getNodeDatabase();
  NodeCommMgr& nodeCommMgr = problemStructure_->getNodeCommMgr();

  int numNodes = nodeDB.getNumNodeDescriptors();
  int numRemoteNodes = nodeCommMgr.getSharedNodeIDs().size() -
                         nodeCommMgr.getLocalNodeIDs().size();
  numNodes -= numRemoteNodes;

  std::vector<int> numElemsPerBlock(numBlocks);
  std::vector<int> numNodesPerElem(numBlocks);
  std::vector<int> numDofPerNode(0,1);

  for(int blk=0; blk<numBlocks; blk++) {
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor_index(blk, block) );

    numElemsPerBlock[blk] = block->getNumElements();
    numNodesPerElem[blk]  = block->getNumNodesPerElement();

    int* fieldsPerNode = block->fieldsPerNodePtr();
    int** fieldIDsTable = block->fieldIDsTablePtr();

    for(int nn=0; nn<numNodesPerElem[blk]; nn++) {
      if (fieldsPerNode[nn] <= 0) ERReturn(-1);
      numDofPerNode.push_back(0);
      int indx = numDofPerNode.size()-1;
      
      for(int nf=0; nf<fieldsPerNode[nn]; nf++) {
        numDofPerNode[indx] += problemStructure_->getFieldSize(fieldIDsTable[nn][nf]);
      }
    }
  }

  for(int i=0; i<numBlocks; i++) {
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor_index(i, block) );

    if (block->getNumElements() == 0) continue;

    ConnectivityTable& ctbl =
      problemStructure_->getBlockConnectivity(block->getGlobalBlockID());

    std::vector<int> cNodeList(block->getNumNodesPerElement());

    int* fieldsPerNode = block->fieldsPerNodePtr();
    int** fieldIDsTable = block->fieldIDsTablePtr();

    numDofPerNode.resize(0);
    for(int nn=0; nn<numNodesPerElem[i]; nn++) {
      if (fieldsPerNode[nn] <= 0) ERReturn(-1);
      numDofPerNode.push_back(0);
      int indx = numDofPerNode.size()-1;

      for(int nf=0; nf<fieldsPerNode[nn]; nf++) {
        numDofPerNode[indx] += problemStructure_->getFieldSize(fieldIDsTable[nn][nf]);
      }
    }

    int nodesPerElement = block->getNumNodesPerElement();
    NodeDescriptor** elemNodePtrs = &((*ctbl.elem_conn_ptrs)[0]);
    int offset = 0;
    for(int j=0; j<block->getNumElements(); j++) {

      for(int k=0; k<nodesPerElement; k++) {
        NodeDescriptor* node = elemNodePtrs[offset++];
        cNodeList[k] = node->getNodeNumber();
      }

      int elemID = ctbl.elemNumbers[j];
      int* nodeNumbers = &cNodeList[0];
      debugOutput("#LinSysCoreFilter calling lsc->setConnectivities");
      CHK_ERR( lsc_->setConnectivities(i, 1, nodesPerElement,
                                       &elemID, &nodeNumbers) );
    }
  }

  connectivitiesInitialized_ = true;

  return(FEI_SUCCESS);
}

//==============================================================================
void LinSysCoreFilter::setLinSysCoreCREqns()
{
   int err, i=0;

   std::vector<int> iwork;

   std::map<GlobalID,ConstraintType*>::const_iterator
     cr_iter = problemStructure_->getMultConstRecords().begin(),
     cr_end = problemStructure_->getMultConstRecords().end();

   while(cr_iter != cr_end) {
      ConstraintType& constraint = *((*cr_iter).second);
      int numNodesPerCR = constraint.getMasters().size();
      int meqn = constraint.getEqnNumber();

      std::vector<GlobalID>& nodeID_vec = constraint.getMasters();
      GlobalID* nodeIDPtr = &nodeID_vec[0];

      if ((int)iwork.size() < 2*numNodesPerCR) {
        iwork.resize(2*numNodesPerCR);
      }

      int* nodeList = &(iwork[0]);

      int* eqnList = nodeList+numNodesPerCR;

      std::vector<int>& fieldIDs_vec = constraint.getMasterFieldIDs();
      int* fieldIDs = &fieldIDs_vec[0];

      for(int k=0; k<numNodesPerCR; k++) {
        const NodeDescriptor *node = Filter::findNode(nodeIDPtr[k]);
        if(node == NULL)
        {
          nodeList[k] = -1; // Indicates that the node wasn't found
        }
        else
        {
          nodeList[k] = node->getNodeNumber();
        }

        int eqn = -1; // Indicates that the equation wasn't found.
        if ( node ) {
          node->getFieldEqnNumber(fieldIDs[k], eqn);
        }
        eqnList[k] = eqn;
      }

      int crMultID = constraint.getConstraintID() + i++;
      if (Filter::logStream() != NULL) {
        FEI_OSTREAM& os = *logStream();
        os << "#LinSysCoreFilter calling lsc_->setMultCREqns"<<FEI_ENDL;
        os << "#  multiplier eqn: " << meqn << ", columns: ";
        for(int j=0; j<numNodesPerCR; ++j) os << eqnList[j] << " ";
        os << FEI_ENDL;
      }

      err = lsc_->setMultCREqns(crMultID, 1, numNodesPerCR,
                                &nodeList, &eqnList,
                                fieldIDs, &meqn);
      if (err) voidERReturn;
      ++cr_iter;
   }

   LinSysCore_flexible* lscf = dynamic_cast<LinSysCore_flexible*>(lsc_);
   if (lscf != NULL) {
     debugOutput("LinSysCoreFilter calling lscf->setMultCRComplete");
     err = lscf->setMultCRComplete();
     if (err != 0) {
       fei::console_out() << "LinSysCoreFilter::setLinSysCoreCREqns ERROR returned from "
            << "lscf->setMultCRComplete()" << FEI_ENDL;
     }
   }

   //
   //Now the penalty CRs...
   //
   cr_iter = problemStructure_->getPenConstRecords().begin();
   cr_end = problemStructure_->getPenConstRecords().end();

   while(cr_iter != cr_end) {
      ConstraintType& crset = *((*cr_iter).second);
      int numNodesPerCR = crset.getMasters().size();

      std::vector<GlobalID>& nodeIDs_vec = crset.getMasters();
      GlobalID* nodeIDsPtr = &nodeIDs_vec[0];

      if ((int)iwork.size() < 2*numNodesPerCR) {
        iwork.resize(2*numNodesPerCR);
      }

      int* nodeList = &(iwork[0]);

      int* eqnList = nodeList+numNodesPerCR;

      std::vector<int>& fieldIDs_vec = crset.getMasterFieldIDs();
      int* fieldIDs = &fieldIDs_vec[0];

      for(int k=0; k<numNodesPerCR; k++) {
        const NodeDescriptor& node = Filter::findNodeDescriptor(nodeIDsPtr[k]);
        nodeList[k] = node.getNodeNumber();

        int eqn = -1;
        if (!node.getFieldEqnNumber(fieldIDs[k], eqn)) voidERReturn;

        eqnList[k] = eqn;
      }

      int crPenID = crset.getConstraintID() + i;
      err = lsc_->setPenCREqns(crPenID, 1, numNodesPerCR,
                               &nodeList, &eqnList,
                               fieldIDs);
      if (err) voidERReturn;
      ++cr_iter;
   }
}

//==============================================================================
int LinSysCoreFilter::storeNodalColumnCoefs(int eqn, const NodeDescriptor& node,
                                            int fieldID, int fieldSize,
                                            double* coefs)
{
  //
  //This function stores the coeficients for 'node' at 'fieldID' at the correct
  //column indices in row 'eqn' of the system matrix.
  //
  if ((localStartRow_ > eqn) || (eqn > localEndRow_)) return(0);

  int eqnNumber = -1;
  if (!node.getFieldEqnNumber(fieldID, eqnNumber)) ERReturn(FEI_ID_NOT_FOUND);

  int numParams = fieldSize;
  if (numParams < 1) {
    return(0);
  }

  if ((int)iworkSpace2_.size() < numParams) {
    iworkSpace2_.resize(numParams);
  }

  int* cols = &iworkSpace2_[0];

  for(int j=0; j<numParams; j++) {
    cols[j] = eqnNumber + j;
  }

  CHK_ERR( giveToMatrix(1, &eqn, numParams, cols, &coefs, ASSEMBLE_SUM) );

  return(FEI_SUCCESS);
}

//==============================================================================
int LinSysCoreFilter::storeNodalRowCoefs(const NodeDescriptor& node,
                                         int fieldID, int fieldSize,
                                         double* coefs, int eqn)
{
  //
  //This function stores coeficients in the equations for 'node', 'fieldID' at
  //column index 'eqn' of the system matrix.
  //
  int eqnNumber = -1;
  if (!node.getFieldEqnNumber(fieldID, eqnNumber)) ERReturn(FEI_ID_NOT_FOUND);

  if ((localStartRow_ > eqnNumber) || (eqnNumber > localEndRow_)) return(0);

  int numParams = fieldSize;

  if (numParams < 1) {
    return(0);
  }

  if ((int)iworkSpace2_.size() < numParams) {
    iworkSpace2_.resize(numParams);
  }

  int* ptRows = &iworkSpace2_[0];

  if ((int)dworkSpace2_.size() < numParams) {
    dworkSpace2_.resize(numParams);
  }

  const double* * values = &dworkSpace2_[0];

  for(int j=0; j<numParams; j++) {
    ptRows[j] = eqnNumber + j;
    values[j] = &(coefs[j]);
  }

  CHK_ERR( giveToMatrix(numParams, ptRows, 1, &eqn, values, ASSEMBLE_SUM) );

  return(FEI_SUCCESS);
}

//==============================================================================
void LinSysCoreFilter::storeNodalSendEqn(const NodeDescriptor& node,
                                         int fieldID, int col,
                                         double* coefs)
{
  //
  //This is a private LinSysCoreFilter function. We can safely assume that
  //it will only be called with a node that is not locally owned.
  //
  //This function tells the eqn comm mgr that we'll be sending contributions
  //to column 'col' for the equations associated with 'fieldID', on 'node', on
  //node's owning processor.
  //
  int proc = node.getOwnerProc();

  int eqnNumber = -1;
  if (!node.getFieldEqnNumber(fieldID, eqnNumber)) voidERReturn;

  int numEqns = problemStructure_->getFieldSize(fieldID);

  for(int i=0; i<numEqns; i++) {
    eqnCommMgr_->addRemoteEqn(eqnNumber+i, proc, &coefs[i], &col, 1);
  }
}

//==============================================================================
void LinSysCoreFilter::storePenNodeSendData(const NodeDescriptor& iNode,
                                            int iField, int iFieldSize,
                                            double* iCoefs,
                                            const NodeDescriptor& jNode,
                                            int jField, int jFieldSize,
                                            double* jCoefs,
                                            double penValue, double CRValue)
{
//
//This function will register with the eqn comm mgr the equations associated
//with iNode, field 'iField' having column indices that are the equations
//associated with jNode, field 'jField', to be sent to the owner of iNode.
//
   int proc = iNode.getOwnerProc();

   int iEqn = -1, jEqn = -1;
   if (!iNode.getFieldEqnNumber(iField, iEqn)) voidERReturn;
   if (!jNode.getFieldEqnNumber(jField, jEqn)) voidERReturn;

   int iNumParams = iFieldSize;
   int jNumParams = jFieldSize;
   if (iNumParams < 1 || jNumParams < 1) {
     fei::console_out() << "FEI ERROR, attempt to store indices for field with non-positive size"
          << " field "<<iField<<", size "<<iNumParams<<", field "<<jField<<", size "
          << jNumParams<<FEI_ENDL;
     voidERReturn;
   }

   if ((int)dworkSpace_.size() < jNumParams) {
     dworkSpace_.resize(jNumParams);
   }

   double* coefs = &dworkSpace_[0];

   if ((int)iworkSpace2_.size() < jNumParams) {
     iworkSpace2_.resize(jNumParams);
   }

   int* cols = &iworkSpace2_[0];

   for(int i=0; i<iNumParams; i++) {
      for(int j=0; j<jNumParams; j++) {
         cols[j] = jEqn + j;
         coefs[j] = penValue*iCoefs[i]*jCoefs[j];
      }

      int row = iEqn + i;
      eqnCommMgr_->addRemoteEqn(row, proc, coefs, cols, jNumParams);

      double rhsValue = penValue*iCoefs[i]*CRValue;
      eqnCommMgr_->addRemoteRHS(row, proc, currentRHS_, rhsValue);
   }
}

//==============================================================================
int LinSysCoreFilter::storePenNodeData(const NodeDescriptor& iNode,
                                       int iField, int iFieldSize,
                                       double* iCoefs,
                                       const NodeDescriptor& jNode,
                                       int jField, int jFieldSize,
                                       double* jCoefs,
                                       double penValue, double CRValue){
//
//This function will add to the local matrix the penalty constraint equations
//associated with iNode at iField, having column indices that are the
//equations associated with jNode at jField.
//Also, add the penalty contribution to the RHS vector.
//
   int iEqn = -1, jEqn = -1;
   if (!iNode.getFieldEqnNumber(iField, iEqn)) ERReturn(FEI_ID_NOT_FOUND);
   if (!jNode.getFieldEqnNumber(jField, jEqn)) ERReturn(FEI_ID_NOT_FOUND);

   int iNumParams = iFieldSize;
   int jNumParams = jFieldSize;
   if (iNumParams < 1 || jNumParams < 1) {
     fei::console_out() << "FEI ERROR, attempt to store indices for field with non-positive size"
          << " field "<<iField<<", size "<<iNumParams<<", field "<<jField<<", size "
          << jNumParams<<FEI_ENDL;
     ERReturn(-1);
   }

   if ((int)dworkSpace2_.size() < iNumParams) {
     dworkSpace2_.resize(iNumParams);
   }

   const double* * coefs = &dworkSpace2_[0];

   if ((int)dworkSpace_.size() < iNumParams * jNumParams) {
     dworkSpace_.resize(iNumParams * jNumParams);
   }

   double* coefList = &dworkSpace_[0];

   if ((int)iworkSpace2_.size() < jNumParams+iNumParams) {
     iworkSpace2_.resize(jNumParams+iNumParams);
   }

   int* cols = &iworkSpace2_[0];
   int* rows = &iworkSpace2_[jNumParams];


   for(int i=0; i<iNumParams; i++) {
      double* coefPtr = coefList + i*jNumParams;
      coefs[i] = coefPtr;

      rows[i] = iEqn + i;

      for(int j=0; j<jNumParams; j++) {
         if (i==0) cols[j] = jEqn + j;
         coefPtr[j] = penValue*iCoefs[i]*jCoefs[j];
      }

      double rhsValue = penValue*iCoefs[i]*CRValue;
      CHK_ERR( giveToRHS(1, &rhsValue, &rows[i], ASSEMBLE_SUM) );
   }

   CHK_ERR( giveToMatrix(iNumParams, rows,
                         jNumParams, cols, coefs, ASSEMBLE_SUM) );

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resetSystem(double s)
{
  //
  //  This puts the value s throughout both the matrix and the vector.
  //
  if (Filter::logStream() != NULL) {
    (*logStream()) << "FEI: resetSystem" << FEI_ENDL << s << FEI_ENDL;
  }

  CHK_ERR( resetTheMatrix(s) );
  CHK_ERR( resetTheRHSVector(s) );

  debugOutput("#LinSysCoreFilter  leaving resetSystem");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::deleteMultCRs()
{
  debugOutput("#LinSysCoreFilter::deleteMultCRs");

  LinSysCore_flexible* lscf = dynamic_cast<LinSysCore_flexible*>(lsc_);
  if (lscf == NULL) {
//    fei::console_out() << "FEI::LinSysCoreFilter: ERROR deleteMultCRs requested, but "
//         << "underlying solver doesn't support this operation." << FEI_ENDL;
    return(-1);
  }

  int err = lscf->resetConstraints(0.0);

  debugOutput("#LinSysCoreFilter leaving deleteMultCRs");

  return(err);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resetTheMatrix(double s)
{
  CHK_ERR( lsc_->resetMatrix(s) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resetTheRHSVector(double s)
{
  CHK_ERR( lsc_->resetRHSVector(s) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resetMatrix(double s) {
//
//  This puts the value s throughout both the matrix and the vector.
//

    debugOutput("FEI: resetMatrix");

    CHK_ERR( resetTheMatrix(s) )

   //clear away any boundary condition data.
   bcManager_->clearAllBCs();

    eqnCommMgr_->resetCoefs();

    debugOutput("#LinSysCoreFilter leaving resetMatrix");

    return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resetRHSVector(double s) {
//
//  This puts the value s throughout the rhs vector.
//

    debugOutput("FEI: resetRHSVector");

    CHK_ERR( resetTheRHSVector(s) )

    //clear away any boundary condition data.
    bcManager_->clearAllBCs();

    eqnCommMgr_->resetCoefs();

    debugOutput("# LinSysCoreFilter leaving resetRHSVector");

    return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resetInitialGuess(double s) {
//
//  This puts the value s throughout the initial guess (solution) vector.
//
  if (Filter::logStream() != NULL) {
    FEI_OSTREAM& os = *logStream();
    os << "FEI: resetInitialGuess" << FEI_ENDL;
    os << "#value to which initial guess is to be set" << FEI_ENDL;
    os << s << FEI_ENDL;
  }

  int* eqns = new int[numReducedRows_];
  double* values = new double[numReducedRows_];
  if (eqns == NULL || values == NULL) return(-1);

  for(int i=0; i<numReducedRows_; ++i) {
    eqns[i] = reducedStartRow_ + i;
    values[i] = s;
  }

  CHK_ERR( lsc_->putInitialGuess(eqns, values, numReducedRows_) );

  delete [] eqns;
  delete [] values;

  debugOutput("# LinSysCoreFilter leaving resetInitialGuess");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::loadNodeBCs(int numNodes,
                                  const GlobalID *nodeIDs,
                                  int fieldID,
                                  const int* offsetsIntoField,
                                  const double* prescribedValues)
{
  //
  //  load boundary condition information for a given set of nodes
  //
  int size = problemStructure_->getFieldSize(fieldID);
  if (size < 1) {
    fei::console_out() << "FEI Warning: loadNodeBCs called for fieldID "<<fieldID
         <<", which was defined with size "<<size<<" (should be positive)."<<FEI_ENDL;
    return(0);
  }

  if (Filter::logStream() != NULL) {
    (*logStream())<<"FEI: loadNodeBCs"<<FEI_ENDL
                     <<"#num-nodes"<<FEI_ENDL<<numNodes<<FEI_ENDL
                     <<"#fieldID"<<FEI_ENDL<<fieldID<<FEI_ENDL
                     <<"#field-size"<<FEI_ENDL<<size<<FEI_ENDL;
    (*logStream())<<"#following lines: nodeID offsetIntoField value "<<FEI_ENDL;

    for(int j=0; j<numNodes; j++) {
      GlobalID nodeID = nodeIDs[j];
      (*logStream())<<static_cast<int>(nodeID)<<"  ";
      (*logStream())<<offsetsIntoField[j]<<"  "<<prescribedValues[j];
      (*logStream())<<FEI_ENDL;
    }
  }

  if (numNodes > 0) newBCData_ = true;

  bcManager_->addBCRecords(numNodes, nodeIDType_, fieldID, nodeIDs,
                           offsetsIntoField, prescribedValues);

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::loadElemBCs(int numElems,
                                  const GlobalID *elemIDs,
                                  int fieldID,
                                  const double *const *alpha,
                                  const double *const *beta,
                                  const double *const *gamma)
{
   return(-1);
}

//------------------------------------------------------------------------------
void LinSysCoreFilter::allocElemStuff()
{
   int nb = problemStructure_->getNumElemBlocks();

   for(int i=0; i<nb; i++) {
     BlockDescriptor* block = NULL;
     int err = problemStructure_->getBlockDescriptor_index(i, block);
     if (err) voidERReturn;

      int numEqns = block->getNumEqnsPerElement();
      if (maxElemRows_ < numEqns) maxElemRows_ = numEqns;
   }

   scatterIndices_.resize(maxElemRows_);

   if (eStiff_ != NULL) delete [] eStiff_;
   if (eStiff1D_ != NULL) delete [] eStiff1D_;
   if (eLoad_ != NULL) delete [] eLoad_;

   eStiff_ = new double*[maxElemRows_];
   eStiff1D_ = new double[maxElemRows_*maxElemRows_];

   if (eStiff_ == NULL || eStiff1D_ == NULL) voidERReturn

   for(int r=0; r<maxElemRows_; r++) {
      eStiff_[r] = eStiff1D_ + r*maxElemRows_;
   }

   eLoad_ = new double[maxElemRows_];

   if (eLoad_ == NULL) voidERReturn
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::sumInElem(GlobalID elemBlockID,
                                GlobalID elemID,
                                const GlobalID* elemConn,
                                const double* const* elemStiffness,
                                const double* elemLoad,
                                int elemFormat)
{
  if (Filter::logStream() != NULL && outputLevel_ > 2) {
    (*logStream()) << "FEI: sumInElem" << FEI_ENDL <<"#blkID" << FEI_ENDL
                      << static_cast<int>(elemBlockID) << FEI_ENDL
                      << "#elID" << FEI_ENDL << static_cast<int>(elemID) << FEI_ENDL;
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
    int numNodes = block->getNumNodesPerElement();
    (*logStream()) << "#n-nodes" << FEI_ENDL << numNodes << FEI_ENDL;
    (*logStream()) << "#nodes" << FEI_ENDL;
    for(int i=0; i<numNodes; ++i) {
      GlobalID nodeID = elemConn[i];
      (*logStream())<<static_cast<int>(nodeID)<<" ";
    }
    (*logStream())<<FEI_ENDL;
  }

  return(generalElemInput(elemBlockID, elemID, elemConn, elemStiffness,
                          elemLoad, elemFormat));
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::sumInElemMatrix(GlobalID elemBlockID,
                                      GlobalID elemID,
                                      const GlobalID* elemConn,
                                      const double* const* elemStiffness,
                                      int elemFormat)
{
  if (Filter::logStream() != NULL && outputLevel_ > 2) {
    (*logStream()) << "FEI: sumInElemMatrix"<<FEI_ENDL
                      << "#blkID" << FEI_ENDL << static_cast<int>(elemBlockID) << FEI_ENDL
                      << "#elID" << FEI_ENDL << static_cast<int>(elemID) << FEI_ENDL;
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
    int numNodes = block->getNumNodesPerElement();
    (*logStream()) << "#n-nodes" << FEI_ENDL << numNodes << FEI_ENDL;
    (*logStream()) << "#nodes" << FEI_ENDL;
    for(int i=0; i<numNodes; ++i) {
      GlobalID nodeID = elemConn[i];
      (*logStream())<<static_cast<int>(nodeID)<<" ";
    }
    (*logStream())<<FEI_ENDL;
  }

  return(generalElemInput(elemBlockID, elemID, elemConn, elemStiffness,
                          NULL, elemFormat));
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::sumInElemRHS(GlobalID elemBlockID,
                           GlobalID elemID,
                           const GlobalID* elemConn,
                           const double* elemLoad)
{
  if (Filter::logStream() != NULL && outputLevel_ > 2) {
    (*logStream()) << "FEI: sumInElemRHS"<<FEI_ENDL<<"#blID" << FEI_ENDL
                      <<(int)elemBlockID << FEI_ENDL
                      << "#elID" << FEI_ENDL << static_cast<int>(elemID) << FEI_ENDL;
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
    int numNodes = block->getNumNodesPerElement();
    (*logStream()) << "#n-nodes" << FEI_ENDL << numNodes << FEI_ENDL;
    (*logStream()) << "#nodes" << FEI_ENDL;
    for(int i=0; i<numNodes; ++i) {
      GlobalID nodeID = elemConn[i];
      (*logStream())<<static_cast<int>(nodeID)<<" ";
    }
    (*logStream())<<FEI_ENDL;
  }

  return(generalElemInput(elemBlockID, elemID, elemConn, NULL,
                          elemLoad, -1));
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::generalElemInput(GlobalID elemBlockID,
                                       GlobalID elemID,
                                       const GlobalID* elemConn,
                                       const double* const* elemStiffness,
                                       const double* elemLoad,
                                       int elemFormat)
{
  (void)elemConn;
  return(generalElemInput(elemBlockID, elemID, elemStiffness, elemLoad,
                          elemFormat) );
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::generalElemInput(GlobalID elemBlockID,
                                       GlobalID elemID,
                                       const double* const* elemStiffness,
                                       const double* elemLoad,
                                       int elemFormat)
{
  //first get the block-descriptor for this elemBlockID...

  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );

  //now allocate our local stiffness/load copy if we haven't already.

  if (maxElemRows_ <= 0) allocElemStuff();

  int numElemRows = block->getNumEqnsPerElement();
  int numBlkElemRows = block->getNumBlkEqnsPerElement();
  int interleave = block->getInterleaveStrategy();

  //an std::vector.resize operation is free if the size is either shrinking or
  //staying the same.

  scatterIndices_.resize(numElemRows);
  blkScatterIndices_.resize(numBlkElemRows*2);

  const double* const* stiff = NULL;
  if (elemStiffness != NULL) stiff = elemStiffness;

  const double* load = NULL;
  if (elemLoad != NULL) load = elemLoad;

  //we'll make a local dense copy of the element stiffness array
  //if the stiffness array was passed in using elemFormat != FEI_DENSE_ROW
  //AND if the stiffness array is non-null.
  //Note that this is not optimal if elemFormat is one of the
  //diagonal or block-diagonal formats.
  if (elemFormat != FEI_DENSE_ROW && stiff != NULL) {
    if (elemFormat == FEI_BLOCK_DIAGONAL_ROW ||
        elemFormat == FEI_BLOCK_DIAGONAL_COL) {
      fei::console_out() << "LinSysCoreFilter::generalElemInput ERROR, elemFormat="
               << elemFormat << " not supported."<<FEI_ENDL;
      ERReturn(-1);
    }

    Filter::copyStiffness(stiff, numElemRows, elemFormat, eStiff_);
    stiff = eStiff_;
  }

  if (stiff != NULL) newMatrixData_ = true;
  if (load != NULL) newVectorData_ = true;

  if (Filter::logStream() != NULL && outputLevel_ > 2) {
    FEI_OSTREAM& os = *logStream();

    if (stiff != NULL || load != NULL) {
      os << "#numRows"<< FEI_ENDL << numElemRows << FEI_ENDL;
    }

    if (stiff != NULL) {
      os << "#elem-stiff (after copy into dense-row format)" << FEI_ENDL;
      for(int i=0; i<numElemRows; i++) {
        const double* stiff_i = stiff[i];
        for(int j=0; j<numElemRows; j++) {
          os << stiff_i[j] << " ";
        }
        os << FEI_ENDL;
      }
    }

    if (load != NULL) {
      os << "#elem-load" << FEI_ENDL;
      for(int i=0; i<numElemRows; i++) {
        os << load[i] << " ";
      }
      os<<FEI_ENDL;
    }

    if (stiff != NULL) {
      os << "#elemformat" << FEI_ENDL << elemFormat << FEI_ENDL;
    }
  }

  //now let's obtain the scatter indices for assembling the equations
  //into their appropriate places in the global matrix and rhs vectors

  int* indPtr = &scatterIndices_[0];
  int* blkIndPtr = &blkScatterIndices_[0];
  int* blkSizesPtr = blkIndPtr + numBlkElemRows;

  bool useBlkEqns = false;
  if (interleave == 0) {
    //interleave==0 is node-major, so we'll get the block-indices too.
    problemStructure_->getScatterIndices_ID(elemBlockID, elemID,
                                            interleave, indPtr,
                                            blkIndPtr, blkSizesPtr);
    int sumBlkSizes = 0;
    for(int ii=0; ii<numBlkElemRows; ++ii) {
      sumBlkSizes += blkSizesPtr[ii];
    }
    if (sumBlkSizes == numElemRows) {
      useBlkEqns = true;
    }
  }
  else {
    //interleave!=0 is field-major, and we'll only bother with point-indices.
    problemStructure_->getScatterIndices_ID(elemBlockID, elemID,
                                            interleave, indPtr);
  }

  if (stiff != NULL) {
    if (problemStructure_->numSlaveEquations() == 0) {
      //I'm not checking the return-value (error-code) on this call, because
      //I wasn't even calling it until recently, and I'm not sure if all
      //LinearSystemCore implementations even have it implemented.
      lsc_->setStiffnessMatrices(elemBlockID, 1, &elemID,
                                 &stiff, scatterIndices_.size(),
                                 &indPtr);
    }

    int len = scatterIndices_.size();
    if (interleave == 0) {
      CHK_ERR( assembleEqns( len, len, indPtr, indPtr, stiff, true,
                             numBlkElemRows, blkIndPtr,
                             blkSizesPtr, useBlkEqns, ASSEMBLE_SUM ) );
    }
    else {
      CHK_ERR( assembleEqns( len, len, indPtr, indPtr, stiff, true,
                             numBlkElemRows, blkIndPtr,
                             blkSizesPtr, false, ASSEMBLE_SUM ) );
    }
  }

  if (load != NULL) {
    if (problemStructure_->numSlaveEquations() == 0) {
      //I'm not checking the return-value (error-code) on this call, because
      //I wasn't even calling it until recently, and I'm not sure if all
      //LinearSystemCore implementations even have it implemented.
      lsc_->setLoadVectors(elemBlockID, 1, &elemID,
                           &load, scatterIndices_.size(),
                           &indPtr);
    }

    CHK_ERR( assembleRHS(scatterIndices_.size(), indPtr, load, ASSEMBLE_SUM) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putIntoRHS(int IDType,
                          int fieldID,
                          int numIDs,
                          const GlobalID* IDs,
                          const double* rhsEntries)
{
  int fieldSize = problemStructure_->getFieldSize(fieldID);

  rowIndices_.resize(fieldSize*numIDs);
  int checkNumEqns;

  CHK_ERR( problemStructure_->getEqnNumbers(numIDs, IDs, IDType, fieldID,
                                            checkNumEqns,
                                            &rowIndices_[0]));
  if (checkNumEqns != numIDs*fieldSize) {
    ERReturn(-1);
  }

  CHK_ERR( exchangeRemoteEquations() );

  CHK_ERR(assembleRHS(rowIndices_.size(), &rowIndices_[0], rhsEntries, ASSEMBLE_PUT));

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::sumIntoRHS(int IDType,
                          int fieldID,
                          int numIDs,
                          const GlobalID* IDs,
                          const double* rhsEntries)
{
  int fieldSize = problemStructure_->getFieldSize(fieldID);

  rowIndices_.resize(fieldSize*numIDs);
  int checkNumEqns;

  CHK_ERR( problemStructure_->getEqnNumbers(numIDs, IDs, IDType, fieldID,
                                            checkNumEqns,
                                            &rowIndices_[0]));
  if (checkNumEqns != numIDs*fieldSize) {
    ERReturn(-1);
  }

  CHK_ERR(assembleRHS(rowIndices_.size(), &rowIndices_[0], rhsEntries, ASSEMBLE_SUM));

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::implementAllBCs() {
//
// This function will handle the modifications to the stiffness and load
// necessary to enforce nodal boundary conditions.
//
   debugOutput("# implementAllBCs");

   std::vector<int> essEqns;
   std::vector<double> essGamma;

   EqnBuffer bcEqns;

   CHK_ERR( bcManager_->finalizeBCEqns(bcEqns) );

   if (resolveConflictRequested_) {
     CHK_ERR( resolveConflictingCRs(bcEqns) );
   }

   CHK_ERR( eqnCommMgr_->gatherSharedBCs(bcEqns) );

   //now separate the boundary-condition equations into arrays
   fei::FillableMat bcEqns_mat(bcEqns);
   fei::impl_utils::separate_BC_eqns(bcEqns_mat, essEqns, essGamma);

   std::vector<double> essAlpha(essEqns.size(), 1);

   exchangeRemoteBCs(essEqns, essAlpha, essGamma);

   if (essEqns.size() > 0) {
      CHK_ERR( enforceEssentialBCs(&essEqns[0],
                                   &essAlpha[0],
                                   &essGamma[0], essEqns.size()) );
   }

   debugOutput("#LinSysCoreFilter leaving implementAllBCs");
   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::enforceEssentialBCs(const int* eqns,
                                          const double* alpha,
                                          const double* gamma,
                                          int numEqns)
{
  int* cc_eqns = const_cast<int*>(eqns);
  double* cc_alpha = const_cast<double*>(alpha);
  double* cc_gamma = const_cast<double*>(gamma);

  if (problemStructure_->numSlaveEquations() == 0) {
    CHK_ERR( lsc_->enforceEssentialBC(cc_eqns,
                                      cc_alpha, cc_gamma,
                                      numEqns) );
  }
  else {
    std::vector<int> reducedEqns(numEqns);
    for(int i=0; i<numEqns; i++) {
      problemStructure_->translateToReducedEqn(eqns[i], reducedEqns[i]);
    }

    CHK_ERR( lsc_->enforceEssentialBC(&reducedEqns[0], cc_alpha, cc_gamma, numEqns) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::enforceRemoteEssBCs(int numEqns, const int* eqns, 
                           const int* const* colIndices, const int* colIndLens,
                           const double* const* BCcoefs)
{
  //These eqn-numbers were reduced to the slave-eqn space by
  //LinSysCore::exchangeRemoteBCs, which is the only function that calls THIS
  //function, so we can simply pass them straight on in to LinearSystemCore.
  //

  int* cc_eqns = const_cast<int*>(eqns);
  int** cc_colIndices = const_cast<int**>(colIndices);
  int* cc_colIndLens = const_cast<int*>(colIndLens);
  double** cc_BCcoefs = const_cast<double**>(BCcoefs);

  CHK_ERR( lsc_->enforceRemoteEssBCs(numEqns, cc_eqns, cc_colIndices,
                                     cc_colIndLens, cc_BCcoefs) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::resolveConflictingCRs(EqnBuffer& bcEqns)
{
  int numMultCRs = problemStructure_->getNumMultConstRecords();
  if (numMultCRs < 1) {
    return(0);
  }

  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = problemStructure_->getMultConstRecords().begin(),
    cr_end  = problemStructure_->getMultConstRecords().end();

  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  std::vector<int>& bcEqnNumbers = bcEqns.eqnNumbers();

  double coefs[3];
  int indices[3];
  indices[0] = 0;
  indices[1] = 1;
  indices[2] = 2;

  double fei_eps = 1.e-49;

  while(cr_iter != cr_end) {
    ConstraintType& multCR = *((*cr_iter).second);

    int lenList = multCR.getMasters().size();

    std::vector<GlobalID>& CRNode_vec = multCR.getMasters();
    GlobalID *CRNodePtr = &CRNode_vec[0];
    std::vector<int>& CRField_vec = multCR.getMasterFieldIDs();
    int* CRFieldPtr = &CRField_vec[0];
    std::vector<double>& weights_vec = multCR.getMasterWeights();
    double* weights = &weights_vec[0];

    int offset = 0;
    for(int j=0; j<lenList; ++j) {
      int fieldSize = problemStructure_->getFieldSize(CRFieldPtr[j]);
      for(int k=0; k<fieldSize; ++k) {
        if (std::abs(weights[offset++] + 1.0) < fei_eps) {
          const NodeDescriptor* node = NULL;
          CHK_ERR( nodeDB.getNodeWithID(CRNodePtr[j], node) );
          int eqn = 0;
          node->getFieldEqnNumber(CRFieldPtr[j], eqn);
          eqn += k;

          if (fei::binarySearch(eqn, bcEqnNumbers) >= 0) {
            coefs[0] = 1.0;
            coefs[1] = 0.0;
            coefs[2] = 1.0;
            int crEqn = multCR.getEqnNumber();
            CHK_ERR( bcEqns.addEqn(crEqn, coefs, indices, 3, false) );
          }
        }
      }
    }
    ++cr_iter;
  }

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::exchangeRemoteEquations()
{
  //
  // This function is where processors send local contributions to remote
  // equations to the owners of those equations, and receive remote
  // contributions to local equations.
  //
  debugOutput("#LinSysCoreFilter::exchangeRemoteEquations");

  if (reducedEqnCounter_ > 0) CHK_ERR( assembleReducedEqns() );

  if (reducedRHSCounter_ > 0) CHK_ERR( assembleReducedRHS() );

  int len = 4;
  std::vector<int> flags(len), globalFlags(len);
  flags[0] = newMatrixData_ ? 1 : 0;
  flags[1] = newVectorData_ ? 1 : 0;
  flags[2] = newConstraintData_ ? 1 : 0;
  flags[3] = newBCData_ ? 1 : 0;

  CHK_ERR( fei::GlobalMax(comm_, flags, globalFlags) );

  newMatrixData_     = globalFlags[0] > 0 ? true : false;
  newVectorData_     = globalFlags[1] > 0 ? true : false;
  newConstraintData_ = globalFlags[2] > 0 ? true : false;
  newBCData_         = globalFlags[3] > 0 ? true : false;

  if (newMatrixData_ || newVectorData_ || newConstraintData_) { 

    CHK_ERR( eqnCommMgr_->exchangeEqns(logStream()) );

    needToCallMatrixLoadComplete_ = true;

    //so now the remote contributions should be available, let's get them out
    //of the eqn comm mgr and put them into our local matrix structure.

    debugOutput("#   putting remote contributions into linear system...");

    CHK_ERR( unpackRemoteContributions(*eqnCommMgr_, ASSEMBLE_SUM) );

    eqnCommMgr_->resetCoefs();

    newMatrixData_ = false;
    newVectorData_ = false;
    newConstraintData_ = false;
  }

  firstRemEqnExchange_ = false;

  debugOutput("#LinSysCoreFilter leaving exchangeRemoteEquations");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::unpackRemoteContributions(EqnCommMgr& eqnCommMgr,
                                                int assemblyMode)
{
  int numRecvEqns = eqnCommMgr.getNumLocalEqns();
  std::vector<int>& recvEqnNumbers = eqnCommMgr.localEqnNumbers();
  std::vector<fei::CSVec*>& recvEqns = eqnCommMgr.localEqns();
  std::vector<std::vector<double>*>& recvRHSs = *(eqnCommMgr.localRHSsPtr());

  bool newCoefs = eqnCommMgr.newCoefData();
  bool newRHSs = eqnCommMgr.newRHSData();

  int i;
  double** coefs = new double*[numRecvEqns];

  for(i=0; i<numRecvEqns; i++) {
    coefs[i] = &(recvEqns[i]->coefs()[0]);
  }

  for(i=0; i<numRecvEqns; i++) {

    int eqn = recvEqnNumbers[i];
    if ((reducedStartRow_ > eqn) || (reducedEndRow_ < eqn)) {
      fei::console_out() << "LinSysCoreFilter::unpackRemoteContributions: ERROR, recvEqn "
           << eqn << " out of range. (localStartRow_: " << reducedStartRow_
           << ", localEndRow_: " << reducedEndRow_ << ", localRank_: "
           << localRank_ << ")" << FEI_ENDL;
      MPI_Abort(comm_, -1);
    }

    for(size_t ii=0; ii<recvEqns[i]->size(); ii++) {
      if (coefs[i][ii] > 1.e+200) {
        fei::console_out() << localRank_ << ": LinSysCoreFilter::unpackRemoteContributions: "
             << "WARNING, coefs["<<i<<"]["<<ii<<"]: " << coefs[i][ii]
             << FEI_ENDL;
        MPI_Abort(comm_, -1);
      }
    }

    if (recvEqns[i]->size() > 0 && newCoefs) {
      //contribute this equation to the matrix,
      CHK_ERR( giveToLocalReducedMatrix(1, &(recvEqnNumbers[i]),
                                        recvEqns[i]->size(),
                                        &(recvEqns[i]->indices()[0]),
                                        &(coefs[i]), assemblyMode ) );
    }

    //and now the RHS contributions.
    if (newRHSs) {
      for(int j=0; j<numRHSs_; j++) {
        CHK_ERR( giveToLocalReducedRHS(1, &( (*(recvRHSs[i]))[j] ),
                                       &eqn, assemblyMode) );
      }
    }
  }

  delete [] coefs;

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::exchangeRemoteBCs(std::vector<int>& essEqns,
                                        std::vector<double>& essAlpha,
                                        std::vector<double>& essGamma)
{
  //we need to make sure that the right thing happens for essential
  //boundary conditions that get applied to nodes on elements that touch
  //a processor boundary. (Note that for this case, the BC node itself doesn't
  //touch the processor boundary.) For an essential boundary condition, the row
  //and column of the corresponding equation must be diagonalized. If there is
  //a processor boundary on any side of the element that contains the node,
  //then there are column contributions to the matrix on the other processor.
  //That other processor must be notified and told to make the adjustments
  //necessary to enforce the boundary condition.

  std::vector<int>* eqns = &essEqns;

  if (problemStructure_->numSlaveEquations() > 0) {
    int numEqns = essEqns.size();
    eqns = new std::vector<int>(numEqns);
    int* eqnsPtr = &(*eqns)[0];

    for(int ii=0; ii<numEqns; ++ii) {
      problemStructure_->translateToReducedEqn(essEqns[ii], eqnsPtr[ii]);
    }
  }

  FEI_OSTREAM* dbgOut = NULL;
  if (Filter::logStream() != NULL) {
    dbgOut = logStream();
  }

  eqnCommMgr_->exchangeRemEssBCs(&(*eqns)[0], eqns->size(),
                                 &essAlpha[0], &essGamma[0],
                                 comm_, dbgOut);

  int numRemoteEssBCEqns = eqnCommMgr_->getNumRemEssBCEqns();
  if (numRemoteEssBCEqns > 0) {
    std::vector<int>& remEssBCEqnNumbers = eqnCommMgr_->remEssBCEqnNumbersPtr();
    fei::CSVec** remEssBCEqns = &(eqnCommMgr_->remEssBCEqns()[0]);
    std::vector<int> remEssBCEqnLengths(remEssBCEqnNumbers.size());

    int** indices = new int*[numRemoteEssBCEqns];
    double** coefs = new double*[numRemoteEssBCEqns];

    for(int i=0; i<numRemoteEssBCEqns; i++) {
      coefs[i] = &(remEssBCEqns[i]->coefs()[0]);
      indices[i] = &(remEssBCEqns[i]->indices()[0]);
      remEssBCEqnLengths[i] = remEssBCEqns[i]->size();
    }

    CHK_ERR( enforceRemoteEssBCs(numRemoteEssBCEqns,
                                 &remEssBCEqnNumbers[0], indices,
                                 &remEssBCEqnLengths[0], coefs));

    delete [] indices;
    delete [] coefs;
  }

  if (problemStructure_->numSlaveEquations() > 0) {
    delete eqns;
  }

  debugOutput("#LinSysCoreFilter exchangeRemoteBCs");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::loadCRMult(int CRID,
                                 int numCRNodes,
                                 const GlobalID* CRNodes,
                                 const int* CRFields,
                                 const double* CRWeights,
                                 double CRValue)
{
//
// Load Lagrange multiplier constraint relation data
//
//   Question: do we really need to pass CRNodes again?  Here, I'm going
//            to ignore it for now (i.e., not store it, but just check it), 
//            as it got passed during the initialization phase, so all we'll 
//            do here is check for errors...
//
  if (Filter::logStream() != NULL) {
    FEI_OSTREAM& os = *logStream();
    os<<"FEI: loadCRMult"<<FEI_ENDL;
    os<<"#num-nodes"<<FEI_ENDL<<numCRNodes<<FEI_ENDL;
    os<<"#CRNodes:"<<FEI_ENDL;
    int i;
    for(i=0; i<numCRNodes; ++i) {
      GlobalID nodeID = CRNodes[i];
      os << static_cast<int>(nodeID) << " ";
    }
    os << FEI_ENDL << "#fields:"<<FEI_ENDL;
    for(i=0; i<numCRNodes; ++i) os << CRFields[i] << " ";
    os << FEI_ENDL << "#field-sizes:"<<FEI_ENDL;
    for(i=0; i<numCRNodes; ++i) {
      int size = problemStructure_->getFieldSize(CRFields[i]);
      os << size << " ";
    }
    os << FEI_ENDL<<"#weights:"<<FEI_ENDL;
    int offset = 0;
    for(i=0; i<numCRNodes; ++i) {
      int size = problemStructure_->getFieldSize(CRFields[i]);
      for(int j=0; j<size; ++j) {
        os << CRWeights[offset++] << " ";
      }
    }
    os << FEI_ENDL<<"#CRValue:"<<FEI_ENDL<<CRValue<<FEI_ENDL;
  }

  ConstraintType* multCR = NULL;
  CHK_ERR( problemStructure_->getMultConstRecord(CRID, multCR) );

  int i;

  int lenList = multCR->getMasters().size();
  if (lenList < 1) {
    fei::console_out() << "ERROR in FEI, constraint with ID="<<CRID<<" appears to have"
         <<" a constrained-node list of length "<<lenList<<", should be > 0."<<FEI_ENDL;
    ERReturn(-1);
  }

  //  recall the data stored earlier and ensure that the passed data (here, the
  //  node list) agrees with the initialization data

  std::vector<GlobalID>& CRNode_vec = multCR->getMasters();
  GlobalID *CRNodePtr = &CRNode_vec[0];

  for(i=0; i<lenList; i++) {
    if (CRNodePtr[i] != CRNodes[i]) {
      fei::console_out() << "ERROR in FEI, constraint with ID="<<CRID<<" had different node-list"
           << " in initCRMult than it has in loadCRMult."<<FEI_ENDL;
      ERReturn(-1);
    }
  }

  std::vector<int>& CRField_vec = multCR->getMasterFieldIDs();
  int *CRFieldPtr = &CRField_vec[0];

  for (i = 0; i < lenList; i++) {
    if (CRFieldPtr[i] != CRFields[i]) {
      fei::console_out() <<"ERROR in FEI, constraint with CRID="<<CRID<<" had different field-list"
           <<" in initCRMult than it has in loadCRMult."<<FEI_ENDL;
      ERReturn(-1);
    }
  }

  newConstraintData_ = true;

  if ((int)iworkSpace_.size() < lenList) {
    iworkSpace_.resize(lenList);
  }

  int* fieldSizes = &iworkSpace_[0];

  for (i = 0; i < lenList; i++) {
    int numSolnParams = problemStructure_->getFieldSize(CRFields[i]);
    assert(numSolnParams >= 0);
    fieldSizes[i] = numSolnParams;
  }

  std::vector<double>& CRWeightArray = multCR->getMasterWeights();

  int offset = 0;

  try {

  for(i = 0; i < lenList; i++) {
    for(int j = 0; j < fieldSizes[i]; j++) {
      CRWeightArray.push_back(CRWeights[offset++]);
    }
  }

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  multCR->setRHSValue(CRValue);
  double* CRWeightsPtr = &CRWeightArray[0];

//  next, perform assembly of the various terms into the system arrays
//  (this is a good candidate for a separate function...)

  int irow = multCR->getEqnNumber();
  double zero = 0.0;
  double* zeroPtr = &zero;
  CHK_ERR( giveToMatrix(1, &irow, 1, &irow, &zeroPtr, ASSEMBLE_PUT) );

  CHK_ERR( giveToRHS(1, &(CRValue), &irow, ASSEMBLE_PUT));

  offset = 0;
  for(int j = 0; j < lenList; j++) {
    int myFieldID = CRFields[j];

    const NodeDescriptor& node = Filter::findNodeDescriptor(CRNodePtr[j]);

    //first, store the column coeficients for equation irow, the
    //constraint's equation.
    storeNodalColumnCoefs(irow, node, myFieldID, fieldSizes[j],
                          &(CRWeightsPtr[offset]));


    //next, store store the transpose of the above. i.e., column irow,
    //in equations associated with 'node' at 'myFieldID'.

    if (node.getOwnerProc() == localRank_) {

      storeNodalRowCoefs(node, myFieldID, fieldSizes[j],
                         &(CRWeightsPtr[offset]), irow);
    }
    else {

      storeNodalSendEqn(node, myFieldID, irow, &(CRWeightsPtr[offset]));
    }

    offset += fieldSizes[j];
  }
    
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::loadCRPen(int CRID, 
                                int numCRNodes,
                                const GlobalID* CRNodes,
                                const int* CRFields,
                                const double* CRWeights,
                                double CRValue,
                                double penValue)
{
  //
  // Load penalty constraint relation data
  //

  debugOutput("FEI: loadCRPen");

  ConstraintType* penCR = NULL;
  CHK_ERR( problemStructure_->getPenConstRecord(CRID, penCR) );

  int i;
  int lenList = penCR->getMasters().size();
  if (lenList < 1) {
    fei::console_out() << "ERROR in FEI, constraint with ID="<<CRID<<" appears to have"
         <<" a constrained-node list of length "<<lenList<<", should be > 0."<<FEI_ENDL;
    ERReturn(-1);
  }

  // recall the data stored earlier and ensure that the passed data (here,
  // the node list) agrees with the initialization data

  std::vector<GlobalID>& CRNode_vec = penCR->getMasters();
  GlobalID* CRNodePtr = &CRNode_vec[0];
                                  
  for(int j = 0; j < lenList; j++) {
    if (CRNodePtr[j] != CRNodes[j]) {
      fei::console_out() << "ERROR in FEI, constraint with ID="<<CRID<<" had different node-list"
           << " in initCRPen than it has in loadCRPen."<<FEI_ENDL;
      ERReturn(-1);
    }
  }

  newConstraintData_ = true;

  //  store the weights and rhs-value in the constraint records.

  if ((int)iworkSpace_.size() < lenList) {
    iworkSpace_.resize(lenList);
  }

  int* fieldSizes = &iworkSpace_[0];

  for (i = 0; i < lenList; i++) {
    int numSolnParams = problemStructure_->getFieldSize(CRFields[i]);
    assert(numSolnParams >= 0);
    fieldSizes[i] = numSolnParams;
  }

  std::vector<double>& CRWeightArray = penCR->getMasterWeights();

  try {

  int offset = 0;
  for (i = 0; i < lenList; i++) {
    for (int j = 0; j < fieldSizes[i]; j++) {
      CRWeightArray.push_back(CRWeights[offset++]);
    }
  }

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  penCR->setRHSValue(CRValue);

  double* CRWeightPtr = &CRWeightArray[0];

  int ioffset = 0, joffset = 0;
  for(i = 0; i < lenList; i++) {
    GlobalID iNodeID = CRNodePtr[i];
    int iField = CRFields[i];

    const NodeDescriptor& iNode = Filter::findNodeDescriptor(iNodeID);
    double* iweights = &(CRWeightPtr[ioffset]);
    ioffset += fieldSizes[i];

    joffset = 0;
    for (int j = 0; j < lenList; j++) {
      GlobalID jNodeID = CRNodePtr[j];
      int jField = CRFields[j];

      const NodeDescriptor& jNode = Filter::findNodeDescriptor(jNodeID);
      double* jweights = &(CRWeightPtr[joffset]);
      joffset += fieldSizes[j];

      double rhsValue = CRValue;
      if (j < lenList-1) {
        rhsValue = 0.0;
      }

      if (iNode.getOwnerProc() == localRank_) {

        storePenNodeData(iNode, iField, fieldSizes[i], iweights,
                         jNode, jField, fieldSizes[j], jweights,
                         penValue, rhsValue);
      }
      else {
        storePenNodeSendData(iNode, iField, fieldSizes[i], iweights,
                             jNode, jField, fieldSizes[j], jweights,
                             penValue, rhsValue);
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::parameters(int numParams, const char *const* paramStrings) {
//
// this function takes parameters for setting internal things like solver
// and preconditioner choice, etc.
//
   if (numParams == 0 || paramStrings == NULL) {
      debugOutput("#LinSysCoreFilter::parameters--- no parameters.");
   }
   else {
      const char* param1 = snl_fei::getParamValue("AZ_matrix_type",
                                                         numParams,
                                                         paramStrings);

      if (param1 != NULL) {
        if (!strcmp(param1, "AZ_VBR_MATRIX") ||
            !strcmp(param1, "blockMatrix")) {
          blockMatrix_ = true;
        }        
      }
      else {
        param1 = snl_fei::getParamValue("matrixType",
                                               numParams, paramStrings);
        if (param1 != NULL) {
          if (!strcmp(param1, "AZ_VBR_MATRIX") ||
              !strcmp(param1, "blockMatrix")) {
            blockMatrix_ = true;
          }        
        }
      }

      param1 = snl_fei::getParamValue("outputLevel",
                                             numParams,paramStrings);
      if ( param1 != NULL){
        std::string str(param1);
        FEI_ISTRINGSTREAM isstr(str);
        isstr >> outputLevel_;
      }

      param1 = snl_fei::getParam("resolveConflict",numParams,paramStrings);
      if ( param1 != NULL){
         resolveConflictRequested_ = true;
      }

      param1 = snl_fei::getParamValue("internalFei", numParams,paramStrings);
      if ( param1 != NULL ){
        std::string str(param1);
        FEI_ISTRINGSTREAM isstr(str);
        isstr >> internalFei_;
      }

      if (Filter::logStream() != NULL) {

        (*logStream())<<"#LinSysCoreFilter::parameters"<<FEI_ENDL
                         <<"# --- numParams: "<< numParams<<FEI_ENDL;
         for(int i=0; i<numParams; i++){
           (*logStream())<<"#------ paramStrings["<<i<<"]: "
                            <<paramStrings[i];
           if (paramStrings[i][strlen(paramStrings[i])-1] != '\n') {
             (*logStream())<<FEI_ENDL;
           }
         }
      }
   }

   CHK_ERR( Filter::parameters(numParams, paramStrings) );

   debugOutput("#LinSysCoreFilter leaving parameters function");

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::loadComplete()
{
  int len = 4;
  std::vector<int> flags(len), globalFlags(len);
  flags[0] = newMatrixData_ ? 1 : 0;
  flags[1] = newVectorData_ ? 1 : 0;
  flags[2] = newConstraintData_ ? 1 : 0;
  flags[3] = newBCData_ ? 1 : 0;

  CHK_ERR( fei::GlobalMax(comm_, flags, globalFlags) );

  newMatrixData_     = globalFlags[0] > 0 ? true : false;
  newVectorData_     = globalFlags[1] > 0 ? true : false;
  newConstraintData_ = globalFlags[2] > 0 ? true : false;
  newBCData_         = globalFlags[3] > 0 ? true : false;

  bool called_exchange = false;
  if (newMatrixData_ || newVectorData_ || newConstraintData_) {
    CHK_ERR( exchangeRemoteEquations() );
    called_exchange = true;
  }

  bool called_implbcs = false;
  if (newBCData_) {
    CHK_ERR( implementAllBCs() );
    called_implbcs = true;
  }

  if (called_exchange || called_implbcs ||needToCallMatrixLoadComplete_) {
    debugOutput("#LinSysCoreFilter calling LinSysCore matrixLoadComplete");

    CHK_ERR( lsc_->matrixLoadComplete() );
    needToCallMatrixLoadComplete_ = false;
  }

  newMatrixData_ = false;
  newVectorData_ = false;
  newConstraintData_ = false;
  newBCData_ = false;

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::residualNorm(int whichNorm, int numFields,
                           int* fieldIDs, double* norms, double& residTime)
{
  //
  //This function can do 3 kinds of norms: infinity-norm (denoted
  //by whichNorm==0), 1-norm and 2-norm.
  //
  debugOutput("FEI: residualNorm");

  if (whichNorm < 0 || whichNorm > 2) return(-1);

  CHK_ERR( loadComplete() );

  std::vector<double> residValues(numReducedRows_, 0.0);

  double start = fei::utils::cpu_time();

  CHK_ERR( formResidual(&(residValues[0]), numReducedRows_) );

  residTime = fei::utils::cpu_time() - start;

  CHK_ERR( Filter::calculateResidualNorms(whichNorm, numFields,
                                          fieldIDs, norms, residValues) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::formResidual(double* residValues, int numLocalEqns)
{
  CHK_ERR( lsc_->formResidual(residValues, numLocalEqns))

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::solve(int& status, double& sTime) {

   debugOutput("FEI: solve");

   CHK_ERR( loadComplete() );

   debugOutput("#LinSysCoreFilter in solve, calling launchSolver...");
 
   double start = fei::utils::cpu_time();

   CHK_ERR( lsc_->launchSolver(status, iterations_) );

   sTime = fei::utils::cpu_time() - start;

   debugOutput("#LinSysCoreFilter... back from solver");
 
   //now unpack the locally-owned shared entries of the solution vector into
   //the eqn-comm-mgr data structures.
   CHK_ERR( unpackSolution() );

   debugOutput("#LinSysCoreFilter leaving solve");

   if (status != 0) return(1);
   else return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::setNumRHSVectors(int numRHSs, int* rhsIDs){

   if (numRHSs < 0) {
      fei::console_out() << "LinSysCoreFilter::setNumRHSVectors: ERROR, numRHSs < 0." << FEI_ENDL;
      ERReturn(-1);
   }

   numRHSs_ = numRHSs;

   rhsIDs_.resize(numRHSs_);
   for(int i=0; i<numRHSs_; i++) rhsIDs_[i] = rhsIDs[i];

  //(we need to set the number of RHSs in the eqn comm manager)
  eqnCommMgr_->setNumRHSs(numRHSs_);

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::setCurrentRHS(int rhsID)
{
   std::vector<int>::iterator iter =
      std::find( rhsIDs_.begin(), rhsIDs_.end(), rhsID);

   if (iter == rhsIDs_.end()) ERReturn(-1)
 
   int index = iter - rhsIDs_.begin();
   currentRHS_ = index;

   lsc_->setRHSID(rhsID);

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::giveToMatrix_symm_noSlaves(int numPtRows,
                                                 const int* ptRowNumbers,
                                                 const double* const* coefs,
                                                 int mode)
{
  for(int i=0; i<numPtRows; i++) {
    int row = ptRowNumbers[i];
    const double* valptr = coefs[i];
    if (row < localStartRow_ || row > localEndRow_) {
      eqnCommMgr_->addRemoteEqn(row, valptr, ptRowNumbers, numPtRows);
      continue;
    }

    if (mode == ASSEMBLE_SUM) {
      if (Filter::logStream() != NULL && 0) {
        FEI_OSTREAM& os = *logStream();
        os << "#  calling sumIntoSystemMatrix, row: " << ptRowNumbers[i]
           << ", columns: ";
        for(int j=0; j<numPtRows; ++j) os << ptRowNumbers[j] << " ";
        os << FEI_ENDL;
      }

      CHK_ERR( lsc_->sumIntoSystemMatrix(1, &(ptRowNumbers[i]),
                                         numPtRows, ptRowNumbers,
                                         &valptr) );
    }
    else {
      CHK_ERR( lsc_->putIntoSystemMatrix(1, &(ptRowNumbers[i]),
                                         numPtRows, ptRowNumbers,
                                         &valptr) );
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::giveToBlkMatrix_symm_noSlaves(int numPtRows,
                                                    const int* ptRowNumbers,
                                                    int numBlkRows,
                                                    const int* blkRowNumbers,
                                                    const int* blkRowSizes,
                                                    const double* const* coefs,
                                                    int mode)
{
  int i;
  if ((int)dworkSpace2_.size() < numPtRows) {
    dworkSpace2_.resize(numPtRows);
  }
  const double* * valptr = &dworkSpace2_[0];
  for(i=0; i<numPtRows; i++) {
    int row = ptRowNumbers[i];
    valptr[i] = coefs[i];
    if (row < localStartRow_ || row > localEndRow_) {
      eqnCommMgr_->addRemoteEqn(row, valptr[i], ptRowNumbers, numPtRows);
      continue;
    }

    if (mode == ASSEMBLE_PUT) {
       CHK_ERR( lsc_->putIntoSystemMatrix(1, &(ptRowNumbers[i]),
                                         numPtRows, ptRowNumbers,
                                         &(valptr[i])) );
   }
  }

  int offset = 0;
  for(i=0; i<numBlkRows; i++) {
    int row = ptRowNumbers[offset];
    if (row < localStartRow_ || row > localEndRow_) {
      offset += blkRowSizes[i];
      continue;
    }

    if (mode == ASSEMBLE_SUM) {
      if (Filter::logStream() != NULL && 0) {
        FEI_OSTREAM& os = *logStream();
        os << "#  calling sumIntoSystemMatrix, row: " << ptRowNumbers[i]
           << ", columns: ";
        for(int j=0; j<numPtRows; ++j) os << ptRowNumbers[j] << " ";
        os << FEI_ENDL;
      }
      
      CHK_ERR(lsc_->sumIntoSystemMatrix(blkRowSizes[i],&(ptRowNumbers[offset]),
                                        numPtRows, ptRowNumbers,
                                        1, &(blkRowNumbers[i]),
                                        numBlkRows, blkRowNumbers,
                                        &(valptr[offset])) );
    }
    
    offset += blkRowSizes[i];
  }

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::giveToMatrix(int numPtRows, const int* ptRows,
                           int numPtCols, const int* ptCols,
                           const double* const* values, int mode)
{
  try {

  if (problemStructure_->numSlaveEquations() == 0) {
    for(int i=0; i<numPtRows; i++) {
      if (ptRows[i] < localStartRow_ || ptRows[i] > localEndRow_) {
        eqnCommMgr_->addRemoteEqn(ptRows[i], values[i], ptCols, numPtCols);
        continue;
      }

      if (mode == ASSEMBLE_SUM) {
        if (Filter::logStream() != NULL && 0) {
          FEI_OSTREAM& os = *logStream();
          os << "#  calling sumIntoSystemMatrix, row: " << ptRows[i]
             << ", columns: ";
          for(int j=0; j<numPtCols; ++j) os << ptCols[j] << " ";
          os << FEI_ENDL;
        }

        CHK_ERR( lsc_->sumIntoSystemMatrix(1, &(ptRows[i]),
                                           numPtCols, ptCols,
                                           &(values[i])) );
      }
      else {
        CHK_ERR( lsc_->putIntoSystemMatrix(1, &(ptRows[i]),
                                           numPtCols, ptCols,
                                           &(values[i])) );
      }
    }
  }
  else {
    iworkSpace_.resize(numPtCols);
    iworkSpace2_.resize(numPtCols);
    int* iworkPtr = &iworkSpace_[0];
    int* iworkPtr2= &iworkSpace2_[0];
    int offset = 0;
    for(int ii=0; ii<numPtCols; ii++) {
      int reducedEqn = -1;
      bool isSlave = problemStructure_->translateToReducedEqn(ptCols[ii],
                                                              reducedEqn);
      if (isSlave) {
        reducedEqn = -1;
        iworkPtr[ii] = reducedEqn;
      }
      else {
        iworkPtr[ii] = reducedEqn;
        iworkPtr2[offset++] = reducedEqn;
      }
    }
    iworkSpace2_.resize(offset);

    for(int i=0; i<numPtRows; i++) {
      int row = ptRows[i];

      int reducedRow;
      bool isSlave = problemStructure_->translateToReducedEqn(row, reducedRow);
      if (isSlave) continue;

      if (reducedStartRow_ > reducedRow || reducedRow > reducedEndRow_) {

        dworkSpace_.resize(0);
        for(int j=0; j<numPtCols; j++) {
          if (iworkSpace_[j]>=0) {
            if (Filter::logStream() != NULL) {
              (*logStream())<<"#  giveToMatrix remote("<<reducedRow<<","
                            <<iworkSpace_[j]<<","<<values[i][j]<<")"<<FEI_ENDL;
            }

            dworkSpace_.push_back(values[i][j]);
          }
        }

        if (mode == ASSEMBLE_SUM) {
          if (Filter::logStream() != NULL) {
            (*logStream())<<"sum"<<FEI_ENDL;
          }

          eqnCommMgr_->addRemoteEqn(reducedRow,
                                    &dworkSpace_[0],
                                    &iworkSpace2_[0],
                                    iworkSpace2_.size());
        }
        else {
          if (Filter::logStream() != NULL) {
            (*logStream())<<"put"<<FEI_ENDL;
          }

          eqnCommMgr_put_->addRemoteEqn(reducedRow,
                                        &dworkSpace_[0],
                                        &iworkSpace2_[0],
                                        iworkSpace2_.size());
        }

        continue;
      }

      for(int j=0; j<numPtCols; j++) {

        int reducedCol = iworkPtr[j];
        if (reducedCol<0) continue;

        double* tmpCoef = const_cast<double*>(&(values[i][j]));

        if (Filter::logStream() != NULL) {
          (*logStream())<< "#  giveToMatrix local("<<reducedRow
                        <<","<<reducedCol<<","<<*tmpCoef<<")"<<FEI_ENDL;
        }

        if (mode == ASSEMBLE_SUM) {
          if (Filter::logStream() != NULL && 0) {
            FEI_OSTREAM& os = *logStream();
            os << "#  calling sumIntoSystemMatrix, row: " << reducedRow
               << ", columns: " << reducedCol << FEI_ENDL;
          }

          CHK_ERR( lsc_->sumIntoSystemMatrix(1, &reducedRow, 1, &reducedCol,
                                             &tmpCoef ) );
        }
        else {
          CHK_ERR( lsc_->putIntoSystemMatrix(1, &reducedRow, 1, &reducedCol,
                                             &tmpCoef ) );
        }
      }
    }
  }

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::giveToLocalReducedMatrix(int numPtRows, const int* ptRows,
                                       int numPtCols, const int* ptCols,
                                       const double* const* values, int mode)
{
  bool specialCase = (!firstRemEqnExchange_ && newConstraintData_
                      && !newMatrixData_) ? true : false;

  double fei_min = std::numeric_limits<double>::min();

  for(int i=0; i<numPtRows; i++) {

    if (mode == ASSEMBLE_SUM) {
      const double* values_i = values[i];

      for(int j=0; j<numPtCols; ++j) {
        if (specialCase && std::abs(values_i[j]) < fei_min) continue;

        const double* valPtr = &(values_i[j]);
        CHK_ERR( lsc_->sumIntoSystemMatrix(1, &(ptRows[i]), 1, &(ptCols[j]),
                                           &valPtr) );
      }
    }
    else {
      CHK_ERR( lsc_->putIntoSystemMatrix(1, &(ptRows[i]), numPtCols, ptCols,
                                         &(values[i])) );
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::sumIntoMatrix(fei::CSRMat& mat)
{
  const std::vector<int>& rowNumbers = mat.getGraph().rowNumbers;
  const std::vector<int>& rowOffsets = mat.getGraph().rowOffsets;
  const std::vector<int>& pckColInds = mat.getGraph().packedColumnIndices;
  const std::vector<double>& pckCoefs = mat.getPackedCoefs();

  for(size_t i=0; i<rowNumbers.size(); ++i) {
    int row = rowNumbers[i];
    int offset = rowOffsets[i];
    int rowlen = rowOffsets[i+1]-offset;
    const int* indices = &pckColInds[offset];
    const double* coefs = &pckCoefs[offset];

    if (giveToMatrix(1, &row, rowlen, indices, &coefs, ASSEMBLE_SUM) != 0) {
      throw std::runtime_error("fei::impl_utils::add_to_matrix ERROR in matrix.sumIn.");
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getFromMatrix(int numPtRows, const int* ptRows,
                            const int* rowColOffsets, const int* ptCols,
                            int numColsPerRow, double** values)
{
  //This function may be attempting to retrieve matrix rows that are not
  //locally owned. If those rows correspond to finite-element nodes that we
  //share, AND if the owning processor is also making this function call, then
  //we can communicate with that processor and obtain those matrix rows.
  //

  ProcEqns remoteProcEqns;

  //Let's populate this ProcEqns object with the remote equations and the procs
  //that we need to receive the remote equations from.
  for(int re=0; re<numPtRows; re++) {
    int eqn = ptRows[re];
    int owner = problemStructure_->getOwnerProcForEqn(eqn);
    if (owner == localRank_) continue;

    remoteProcEqns.addEqn(eqn, owner);
  }

  //so now we know which of the requested equations are remotely owned, and we
  //know which processors own them.
  //Next we're going to need to know which locally-owned equations are needed
  //by other processors.
  ProcEqns localProcEqns;
  CHK_ERR( eqnCommMgr_->mirrorProcEqns(remoteProcEqns, localProcEqns) )

  //ok, now we know which local equations we'll need to send, so let's extract
  //those from the matrix
  EqnBuffer localEqns;
  CHK_ERR( getEqnsFromMatrix(localProcEqns, localEqns) )

  //now we can set the lengths in localProcEqns.
  std::vector<int>& eqnNumbers = localEqns.eqnNumbers();
  fei::CSVec** localEqnsPtr = (localEqns.eqns().size() ? &(localEqns.eqns()[0]) : 0);
  std::vector<int> eqnLengths(eqnNumbers.size());
  for(size_t i=0; i<eqnNumbers.size(); ++i) {
    eqnLengths[i] = localEqnsPtr[i]->size();
  }

  localProcEqns.setProcEqnLengths(&eqnNumbers[0], &eqnLengths[0],
                                  eqnNumbers.size());

  //now mirror those lengths in the remoteProcEqns objects to get ready for the
  //all-to-all exchange of equation data.
  CHK_ERR( eqnCommMgr_->mirrorProcEqnLengths(localProcEqns, remoteProcEqns) );

  EqnBuffer remoteEqns;
  //we're now ready to do the exchange.
  CHK_ERR( EqnCommMgr::exchangeEqnBuffers(comm_, &localProcEqns, &localEqns,
                                          &remoteProcEqns, &remoteEqns, false));

  std::vector<int>& remEqnNumbers = remoteEqns.eqnNumbers();
  fei::CSVec** remEqnsPtr = (remoteEqns.eqns().size() ? &(remoteEqns.eqns()[0]) : 0);
  std::vector<fei::CSVec*>& remEqns   = remoteEqns.eqns();

  //now we're ready to fill the values array with the remote coefficients.
  for(int i=0; i<numPtRows; i++) {
    int row = ptRows[i];

    int eqnIndex = fei::binarySearch(row, remEqnNumbers);

    //if eqnIndex < 0, this is a local equation, so skip to the next loop iter.
    if (eqnIndex < 0) continue;

    //the equation is remote, so stride across it copying out the coefs.
    //if ptCols is NULL, then we're going to copy all coefficients (the whole
    //row) into 'values'.
    if (ptCols == NULL) {
      for(size_t j=0; j<remEqnsPtr[eqnIndex]->size(); j++) {
        values[i][j] = remEqns[eqnIndex]->coefs()[j];
      }
      continue;
    }

    for(int j=0; j<numColsPerRow; j++) {
      int offset = rowColOffsets[i] + j;
      int colIndex = fei::binarySearch(ptCols[offset], remEqns[eqnIndex]->indices());
      if (colIndex < 0) ERReturn(-1);

      values[i][j] = remEqns[eqnIndex]->coefs()[colIndex];
    }
  }

  //and now, get the local stuff out of the matrix.
  for(int i=0; i<numPtRows; i++) {
    int row = ptRows[i];
    if (row < localStartRow_ || localEndRow_ < row) continue;

    int rowLen = 0, checkRowLen;
    CHK_ERR( lsc_->getMatrixRowLength(row, rowLen) )
      if (rowLen <= 0) ERReturn(-1);

    //for each local row, establish some temp arrays and get the row from
    //the matrix.

    std::vector<double> coefs(rowLen);
    std::vector<int> indices(rowLen);

    CHK_ERR( lsc_->getMatrixRow(row, &coefs[0], &indices[0], rowLen, checkRowLen) );
    if (rowLen != checkRowLen) ERReturn(-1);

    //now stride across the list of requested column-indices, and find the
    //corresponding location in the matrix row. Copy that location into the
    //values array.

    //again, if ptCols is NULL, then we're going to copy all coefficients 
    //(the whole row) into 'values'.
    if (ptCols == NULL) {
      for(int j=0; j<rowLen; j++) {
        values[i][j] = coefs[j];
      }
      continue;
    }

    for(int j=0; j<numColsPerRow; j++) {
      std::vector<int>::iterator iter =
          std::find(indices.begin(), indices.end(), ptCols[rowColOffsets[i]+j]);
      if (iter == indices.end()) {
        ERReturn(-1);
      }

      int index = iter - indices.begin();
      values[i][j] = coefs[index];
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getEqnsFromMatrix(ProcEqns& procEqns, EqnBuffer& eqnData)
{
  //Given a ProcEqns object containing lists of equation-numbers, get the data
  //for those equations from the local portion of the matrix and store that data
  //in the eqnData object.

  std::vector<std::vector<int>*>& eqnNumbers = procEqns.procEqnNumbersPtr();

  for(unsigned p=0; p<eqnNumbers.size(); p++) {
    for(unsigned i=0; i<eqnNumbers[p]->size(); i++) {
      int eqn = (*(eqnNumbers[p]))[i];

      if (localStartRow_ > eqn || localEndRow_ < eqn) continue;

      //if this equation is already in eqnData, then don't put it in again...
      std::vector<int>& eqnDataEqns = eqnData.eqnNumbers();
      if (fei::binarySearch(eqn, eqnDataEqns) >= 0) continue;

      int len = 0;
      CHK_ERR( lsc_->getMatrixRowLength(eqn, len) );

      if (len <= 0) continue;
      std::vector<double> coefs(len);
      std::vector<int> indices(len);
      int outputLen = 0;

      CHK_ERR( lsc_->getMatrixRow(eqn, &coefs[0], &indices[0],
                                  len, outputLen) );
      if (outputLen != len) ERReturn(-1);

      CHK_ERR( eqnData.addEqn(eqn, &coefs[0], &indices[0], len, false) );
    }
  }
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getEqnsFromRHS(ProcEqns& procEqns, EqnBuffer& eqnData)
{
  //Given a ProcEqns object containing lists of equation-numbers, get the data
  //for those equations from the local portion of the RHS vector and store that
  // data in the eqnData object. We're only storing rhs coefs in an EqnBuffer
  //that was designed for also storing equations with column-indices. So we'll
  //put a bogus column-index in with each equation just to make sure the
  //EqnBuffer does the right stuff internally...

  int numSendProcs = procEqns.getNumProcs();
  std::vector<int>& eqnsPerProc = procEqns.eqnsPerProcPtr();
  std::vector<std::vector<int>*>& eqnNumbers = procEqns.procEqnNumbersPtr();

  eqnData.setNumRHSs(1);

  for(int p=0; p<numSendProcs; p++) {
    for(int i=0; i<eqnsPerProc[p]; i++) {
      int reducedEqn;
      problemStructure_->translateToReducedEqn((*(eqnNumbers[p]))[i], reducedEqn);

      if (reducedStartRow_ > reducedEqn || reducedEndRow_ < reducedEqn) continue;

      //if this equation is already in eqnData, then don't put it in again...
      std::vector<int>& eqnDataEqns = eqnData.eqnNumbers();
      if (fei::binarySearch(reducedEqn, eqnDataEqns) >= 0) continue;

      double rhsValue;

      CHK_ERR( lsc_->getFromRHSVector(1, &rhsValue, &reducedEqn) );

      int bogusIndex = 19;
      CHK_ERR( eqnData.addIndices(reducedEqn, &bogusIndex, 1) );
      CHK_ERR( eqnData.addRHS(reducedEqn, 0, rhsValue) );
    }
  }
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::giveToRHS(int num, const double* values,
                        const int* indices, int mode)
{
  if (problemStructure_->numSlaveEquations() == 0) {
    for(int i=0; i<num; i++) {
      if (indices[i] < localStartRow_ || indices[i] > localEndRow_) {
        if (mode == ASSEMBLE_SUM) {
          eqnCommMgr_->addRemoteRHS(indices[i], currentRHS_, values[i]);
        }
        else {
          eqnCommMgr_put_->addRemoteRHS(indices[i], currentRHS_, values[i]);
        }

        continue;
      }

      if (mode == ASSEMBLE_SUM) {
        CHK_ERR( lsc_->sumIntoRHSVector(1, &(values[i]), &(indices[i])) );
      }
      else {
        CHK_ERR( lsc_->putIntoRHSVector(1, &(values[i]), &(indices[i])) );
      }
    }
  }
  else {
    for(int i=0; i<num; i++) {
      int reducedEqn;
      bool isSlave = problemStructure_->
        translateToReducedEqn(indices[i], reducedEqn);
      if (isSlave) continue;

      if (reducedEqn < reducedStartRow_ || reducedEqn > reducedEndRow_) {
        if (mode == ASSEMBLE_SUM) {
          eqnCommMgr_->addRemoteRHS(reducedEqn, currentRHS_, values[i]);
        }
        else {
          eqnCommMgr_put_->addRemoteRHS(reducedEqn, currentRHS_, values[i]);
        }

        continue;
      }

      if (mode == ASSEMBLE_SUM) {
        CHK_ERR( lsc_->sumIntoRHSVector(1, &(values[i]), &reducedEqn) );
      }
      else {
        CHK_ERR( lsc_->putIntoRHSVector(1, &(values[i]), &reducedEqn) );
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::giveToLocalReducedRHS(int num, const double* values,
                                    const int* indices, int mode)
{
  for(int i=0; i<num; i++) {

    if (mode == ASSEMBLE_SUM) {
      CHK_ERR( lsc_->sumIntoRHSVector(1, &(values[i]), &(indices[i])) );
    }
    else {
      CHK_ERR( lsc_->putIntoRHSVector(1, &(values[i]), &(indices[i])) );
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::sumIntoRHS(fei::CSVec& vec)
{
  std::vector<int>& indices = vec.indices();
  std::vector<double>& coefs = vec.coefs();

  CHK_ERR( giveToRHS(indices.size(), &coefs[0], &indices[0], ASSEMBLE_SUM) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getFromRHS(int num, double* values, const int* indices)
{
  //We need to do similar things here as we do in getFromMatrix, with respect to
  //communications to obtain values for equations that are remotely owned.

  ProcEqns remoteProcEqns;

  //Let's populate this ProcEqns object with the remote equations and the procs
  //that we need to receive the remote equations from.
  for(int re=0; re<num; re++) {
    int eqn = indices[re];
    int owner = problemStructure_->getOwnerProcForEqn(eqn);
    if (owner==localRank_) continue;

    remoteProcEqns.addEqn(eqn, owner);
  }

  //so now we know which of the requested equations are remotely owned, and we
  //know which processors own them.
  //Next we're going to need to know which locally-owned equations are needed
  //by other processors.
  ProcEqns localProcEqns;
  CHK_ERR( eqnCommMgr_->mirrorProcEqns(remoteProcEqns, localProcEqns) );

  //ok, now we know which equations we'll need to send, so let's extract
  //them from the rhs vector.
  EqnBuffer localEqns;
  CHK_ERR( getEqnsFromRHS(localProcEqns, localEqns) );

  //now we can set the lengths in localProcEqns.
  std::vector<int>& eqnNumbers = localEqns.eqnNumbers();
  fei::CSVec** localEqnsPtr = &(localEqns.eqns()[0]);
  std::vector<int> eqnLengths(eqnNumbers.size());
  for(size_t i=0; i<eqnNumbers.size(); ++i) {
    eqnLengths[i] = localEqnsPtr[i]->size();
  }

  localProcEqns.setProcEqnLengths(&eqnNumbers[0], &eqnLengths[0],
                                  eqnNumbers.size());

  //now mirror those lengths in the remoteProcEqns objects to get ready for the
  //all-to-all exchange of equation data.
  CHK_ERR( eqnCommMgr_->mirrorProcEqnLengths(localProcEqns, remoteProcEqns) );

  EqnBuffer remoteEqns;
  //we're now ready to do the exchange.
  CHK_ERR( EqnCommMgr::exchangeEqnBuffers(comm_, &localProcEqns, &localEqns,
                                           &remoteProcEqns, &remoteEqns, false))

  //now we're ready to get the rhs data we've received from other processors.
  std::vector<int>& remEqnNumbers = remoteEqns.eqnNumbers();
  std::vector<std::vector<double>*>& remRhsCoefs = *(remoteEqns.rhsCoefsPtr());

  for(int i=0; i<num; i++) {
    int row = indices[i];

    int eqnIndex = fei::binarySearch(row, remEqnNumbers);
    if (eqnIndex < 0) continue;

    values[i] = (*(remRhsCoefs[eqnIndex]))[0];
  }

  //and now get the local stuff.
  for(int i=0; i<num; i++) {
    if (indices[i] < localStartRow_ || indices[i] > localEndRow_) continue;

    CHK_ERR( lsc_->getFromRHSVector(num, values, indices) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getEqnSolnEntry(int eqnNumber, double& solnValue)
{
  //This function's task is to retrieve the solution-value for a global
  //equation-number. eqnNumber may or may not be a slave-equation, and may or
  //may not be locally owned. If it is not locally owned, it should at least
  //be shared.
  //return 0 if the solution is successfully retrieved, otherwise return 1.
  //

  //
  //First and probably most common case: there are no slave equations.
  //
  if (problemStructure_->numSlaveEquations() == 0) {
    if (localStartRow_ > eqnNumber || eqnNumber > localEndRow_) {
      //Dig into the eqn-comm-mgr for the shared-remote solution value.
      CHK_ERR( getSharedRemoteSolnEntry(eqnNumber, solnValue) );
    }
    else {
      //It's local, simply get the solution from the assembled linear system.
      CHK_ERR( getReducedSolnEntry( eqnNumber, solnValue ) );
    }
    return(0);
  }

  //
  //If we reach this point, there are slave equations to account for.
  //So first translate this equation into 'assembled-linear-system'
  //equation-numbers.
  //
  int reducedEqn;
  bool isSlave = problemStructure_->translateToReducedEqn(eqnNumber,reducedEqn);

  if (isSlave) {
    //This is a slave-equation, so construct its solution-value as the linear-
    //combination of the master-equations it is defined in terms of.

    std::vector<int>* masterEqns = NULL;
    std::vector<double>* masterCoefs = NULL;
    CHK_ERR( problemStructure_->getMasterEqnNumbers(eqnNumber, masterEqns) );
    CHK_ERR( problemStructure_->getMasterEqnCoefs(eqnNumber, masterCoefs) );

    int len = masterEqns->size();
    solnValue = 0.0;
    CHK_ERR( problemStructure_->getMasterEqnRHS(eqnNumber, solnValue) );

    double coef = 0.0;
    for(int i=0; i<len; i++) {
      int mEqn = (*masterEqns)[i];
      int mReducedeqn;
      problemStructure_->translateToReducedEqn(mEqn, mReducedeqn);

      if (reducedStartRow_ > mReducedeqn || mReducedeqn > reducedEndRow_) {
        CHK_ERR( getSharedRemoteSolnEntry(mReducedeqn, coef) );
      }
      else {
        CHK_ERR( getReducedSolnEntry(mReducedeqn, coef) );
      }
      solnValue += coef * (*masterCoefs)[i];
    }
  }
  else {
    //This is not a slave-equation, so retrieve the solution from either the
    //assembled linear system or the shared-remote data structures.

    if (reducedStartRow_ > reducedEqn || reducedEqn > reducedEndRow_) {
      CHK_ERR( getSharedRemoteSolnEntry(reducedEqn, solnValue) );
    }
    else {
      CHK_ERR( getReducedSolnEntry(reducedEqn, solnValue) );
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getSharedRemoteSolnEntry(int eqnNumber, double& solnValue)
{
  std::vector<int>& remoteEqnNumbers = eqnCommMgr_->sendEqnNumbersPtr();
  double* remoteSoln = eqnCommMgr_->sendEqnSolnPtr();

  int index = fei::binarySearch(eqnNumber, remoteEqnNumbers);
  if (index < 0) {
    fei::console_out() << "LinSysCoreFilter::getSharedRemoteSolnEntry: ERROR, eqn "
         << eqnNumber << " not found." << FEI_ENDL;
    ERReturn(-1);
  }
  solnValue = remoteSoln[index];
  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getReducedSolnEntry(int eqnNumber, double& solnValue)
{
  //We may safely assume that this function is called with 'eqnNumber' that is
  //local in the underlying assembled linear system. i.e., it isn't a slave-
  //equation, it isn't remotely owned, etc.
  //
  CHK_ERR( lsc_->getSolnEntry(eqnNumber, solnValue) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::unpackSolution()
{
  //
  //This function should be called after the solver has returned,
  //and we know that there is a solution in the underlying vector.
  //This function ensures that any locally-owned shared solution values are
  //available on the sharing processors.
  //
  if (Filter::logStream() != NULL) {
    (*logStream())<< "#  entering unpackSolution, outputLevel: "
                     <<outputLevel_<<FEI_ENDL;
  }

  //what we need to do is as follows.
  //The eqn comm mgr has a list of what it calls 'recv eqns'. These are
  //equations that we own, for which we received contributions from other
  //processors. The solution values corresponding to these equations need
  //to be made available to those remote contributing processors.

  int numRecvEqns = eqnCommMgr_->getNumLocalEqns();
  std::vector<int>& recvEqnNumbers = eqnCommMgr_->localEqnNumbers();

  for(int i=0; i<numRecvEqns; i++) {
    int eqn = recvEqnNumbers[i];

    if ((reducedStartRow_ > eqn) || (reducedEndRow_ < eqn)) {
      fei::console_out() << "LinSysCoreFilter::unpackSolution: ERROR, 'recv' eqn (" << eqn
           << ") out of local range." << FEI_ENDL;
      MPI_Abort(comm_, -1);
    }

    double solnValue = 0.0;

    CHK_ERR( getReducedSolnEntry(eqn, solnValue) );

    eqnCommMgr_->addSolnValues(&eqn, &solnValue, 1);
  }

  eqnCommMgr_->exchangeSoln();

  debugOutput("#LinSysCoreFilter leaving unpackSolution");
  return(FEI_SUCCESS);
}
             
//------------------------------------------------------------------------------
void LinSysCoreFilter::setEqnCommMgr(EqnCommMgr* eqnCommMgr)
{
  delete eqnCommMgr_;
  eqnCommMgr_ = eqnCommMgr;
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getBlockNodeSolution(GlobalID elemBlockID,  
                                   int numNodes, 
                                   const GlobalID *nodeIDs, 
                                   int *offsets,
                                   double *results) {
        
   debugOutput("FEI: getBlockNodeSolution");

   int numActiveNodes = problemStructure_->getNumActiveNodes();
   NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

   if (numActiveNodes <= 0) return(0);

   int numSolnParams = 0;

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );

   //Traverse the node list, checking if nodes are associated with this block.
   //If so, put its 'answers' in the results list.

   int offset = 0;
   for(int i=0; i<numActiveNodes; i++) {
     const NodeDescriptor* node_i = NULL;
     nodeDB.getNodeAtIndex(i, node_i);

      if (offset == numNodes) break;

      GlobalID nodeID = nodeIDs[offset];

      //first let's set the offset at which this node's solution coefs start.
      offsets[offset++] = numSolnParams;

      const NodeDescriptor* node = NULL;
      int err = 0;
      //Obtain the NodeDescriptor of nodeID in the activeNodes list...
      //Don't call the getActiveNodeDesc_ID function unless we have to.

      if (node_i!=NULL && nodeID == node_i->getGlobalNodeID()) {
        node = node_i;
      }
      else {
         err = nodeDB.getNodeWithID(nodeID, node);
      }

      //ok. If err is not 0, meaning nodeID is NOT in the
      //activeNodes list, then skip to the next loop iteration.

      if (err != 0) {
        continue;
      }

      int numFields = node->getNumFields();
      const int* fieldIDs = node->getFieldIDList();

      for(int j=0; j<numFields; j++) {
        if (block->containsField(fieldIDs[j])) {
          int size = problemStructure_->getFieldSize(fieldIDs[j]);
          assert(size >= 0);

          int thisEqn = -1;
          node->getFieldEqnNumber(fieldIDs[j], thisEqn);

          double answer;
          for(int k=0; k<size; k++) {
            CHK_ERR( getEqnSolnEntry(thisEqn+k, answer) )
              results[numSolnParams++] = answer;
          }
        }
      }//for(j<numFields)loop
   }

   offsets[numNodes] = numSolnParams;

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getNodalSolution(int numNodes, 
                                       const GlobalID *nodeIDs, 
                                       int *offsets,
                                       double *results)
{
  debugOutput("FEI: getNodalSolution");

  int numActiveNodes = problemStructure_->getNumActiveNodes();
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  if (numActiveNodes <= 0) return(0);

  int numSolnParams = 0;

  //Traverse the node list, checking if nodes are local.
  //If so, put 'answers' in the results list.

  int offset = 0;
  for(int i=0; i<numActiveNodes; i++) {
    const NodeDescriptor* node_i = NULL;
    nodeDB.getNodeAtIndex(i, node_i);

    if (offset == numNodes) break;

    GlobalID nodeID = nodeIDs[offset];

    //first let's set the offset at which this node's solution coefs start.
    offsets[offset++] = numSolnParams;

    const NodeDescriptor* node = NULL;
    int err = 0;
    //Obtain the NodeDescriptor of nodeID in the activeNodes list...
    //Don't call the getNodeWithID function unless we have to.

    if (node_i!=NULL && nodeID == node_i->getGlobalNodeID()) {
      node = node_i;
    }
    else {
      err = nodeDB.getNodeWithID(nodeID, node);
    }

    //ok. If err is not 0, meaning nodeID is NOT in the
    //activeNodes list, then skip to the next loop iteration.

    if (err != 0) {
      continue;
    }

    int numFields = node->getNumFields();
    const int* fieldIDs = node->getFieldIDList();

    for(int j=0; j<numFields; j++) {
      int size = problemStructure_->getFieldSize(fieldIDs[j]);
      assert(size >= 0);

      int thisEqn = -1;
      node->getFieldEqnNumber(fieldIDs[j], thisEqn);

      double answer;
      for(int k=0; k<size; k++) {
        CHK_ERR( getEqnSolnEntry(thisEqn+k, answer) )
          results[numSolnParams++] = answer;
      }
    }//for(j<numFields)loop
  }

  offsets[numNodes] = numSolnParams;

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getBlockFieldNodeSolution(GlobalID elemBlockID,
                                        int fieldID,
                                        int numNodes, 
                                        const GlobalID *nodeIDs, 
                                        double *results)
{
  //Note: if the user-supplied nodeIDs list containts nodes which are not in
  //the specified element-block, then the corresponding positions in the
  //results array are simply not referenced. This is dangerous behavior that
  //hasn't gotten me into trouble yet.

  debugOutput("FEI: getBlockFieldNodeSolution");

  int numActiveNodes = problemStructure_->getNumActiveNodes();
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  if (numActiveNodes <= 0) return(0);

  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  if (fieldSize <= 0) ERReturn(-1);

  if (!block->containsField(fieldID)) {
    fei::console_out() << "LinSysCoreFilter::getBlockFieldNodeSolution WARNING: fieldID " << fieldID
         << " not contained in element-block " << (int)elemBlockID << FEI_ENDL;
    return(1);
  }

   //Traverse the node list, checking if nodes are associated with this block.
   //If so, put the answers in the results list.

   for(int i=0; i<numNodes; i++) {
     const NodeDescriptor* node_i = NULL;
     nodeDB.getNodeAtIndex(i, node_i);

     GlobalID nodeID = nodeIDs[i];

     const NodeDescriptor* node = NULL;
     int err = 0;
     //Obtain the NodeDescriptor of nodeID in the activeNodes list...
     //Don't call the getNodeWithID function unless we have to. (getNodeWithID
     //does a binary-search over all local nodeIDs, while getNodeAtIndex is
     //a direct lookup.) Often the user supplies a nodeIDs list that is in the
     //"natural" order, so we don't need to call getNodeWithID at all.

     if (node_i!=NULL && nodeID == node_i->getGlobalNodeID()) {
       node = node_i;
     }
     else {
       err = nodeDB.getNodeWithID(nodeID, node);
     }

     //ok. If err is not 0, meaning nodeID is NOT in the
     //activeNodes list, then skip to the next loop iteration.

     if (err != 0) {
       continue;
     }

     int eqnNumber = -1;
     bool hasField = node->getFieldEqnNumber(fieldID, eqnNumber);
     if (!hasField) continue;

     int offset = fieldSize*i;
     for(int j=0; j<fieldSize; j++) {
       double answer = 0.0;
       CHK_ERR( getEqnSolnEntry(eqnNumber+j, answer) );
       results[offset+j] = answer;
     }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getNodalFieldSolution(int fieldID,
                                    int numNodes, 
                                    const GlobalID *nodeIDs, 
                                    double *results)
{
  debugOutput("FEI: getNodalFieldSolution");

  int numActiveNodes = problemStructure_->getNumActiveNodes();
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  if (numActiveNodes <= 0) return(0);

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  if (fieldSize <= 0) ERReturn(-1);

  //Traverse the node list, checking if nodes have the specified field.
  //If so, put the answers in the results list.

  for(int i=0; i<numNodes; i++) {
    const NodeDescriptor* node_i = NULL;
    nodeDB.getNodeAtIndex(i, node_i);

    GlobalID nodeID = nodeIDs[i];

    const NodeDescriptor* node = NULL;
    int err = 0;
    //Obtain the NodeDescriptor of nodeID in the activeNodes list...
    //Don't call the getNodeWithID function unless we have to.

    if (node_i!=NULL && nodeID == node_i->getGlobalNodeID()) {
      node = node_i;
    }
    else {
      err = nodeDB.getNodeWithID(nodeID, node);
    }

    //ok. If err is not 0, meaning nodeID is NOT in the
    //activeNodes list, then skip to the next loop iteration.

    if (err != 0) {
      continue;
    }

    int eqnNumber = -1;
    bool hasField = node->getFieldEqnNumber(fieldID, eqnNumber);

    //If this node doesn't have the specified field, then skip to the
    //next loop iteration.
    if (!hasField) continue;

    int offset = fieldSize*i;
    for(int j=0; j<fieldSize; j++) {
      double answer = 0.0;
      CHK_ERR( getEqnSolnEntry(eqnNumber+j, answer) );
      results[offset+j] = answer;
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putBlockNodeSolution(GlobalID elemBlockID,
                                   int numNodes, 
                                   const GlobalID *nodeIDs, 
                                   const int *offsets,
                                   const double *estimates) {
        
   debugOutput("FEI: putBlockNodeSolution");

   int numActiveNodes = problemStructure_->getNumActiveNodes();

   if (numActiveNodes <= 0) return(0);

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );

   NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

   //traverse the node list, checking for nodes associated with this block
   //when an associated node is found, put its 'answers' into the linear system.

   unsigned blk_idx = problemStructure_->getIndexOfBlock(elemBlockID);

   for(int i=0; i<numNodes; i++) {
     const NodeDescriptor* node = NULL;
     int err = nodeDB.getNodeWithID(nodeIDs[i], node);

      if (err != 0) continue;
   
      if (!node->hasBlockIndex(blk_idx)) continue;

      if (node->getOwnerProc() != localRank_) continue;

      int numFields = node->getNumFields();
      const int* fieldIDs = node->getFieldIDList();
      const int* fieldEqnNumbers = node->getFieldEqnNumbers();

      if (fieldEqnNumbers[0] < localStartRow_ ||
          fieldEqnNumbers[0] > localEndRow_) continue;

      int offs = offsets[i];

      for(int j=0; j<numFields; j++) {
         int size = problemStructure_->getFieldSize(fieldIDs[j]);

         if (block->containsField(fieldIDs[j])) {
            for(int k=0; k<size; k++) {
               int reducedEqn;
               problemStructure_->
                 translateToReducedEqn(fieldEqnNumbers[j]+k, reducedEqn);

               CHK_ERR( lsc_->putInitialGuess(&reducedEqn,
                                              &estimates[offs+k], 1) );
            }
         }
         offs += size;
      }//for(j<numFields)loop
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                        int fieldID, 
                                        int numNodes, 
                                        const GlobalID *nodeIDs, 
                                        const double *estimates)
{
   int fieldSize = problemStructure_->getFieldSize(fieldID);

   if (Filter::logStream() != NULL) {
     FEI_OSTREAM& os = *logStream();
     os << "FEI: putBlockFieldNodeSolution" << FEI_ENDL;
     os << "#blkID" << FEI_ENDL << (int)elemBlockID << FEI_ENDL
        << "#fieldID"<<FEI_ENDL << fieldID << FEI_ENDL
        << "#fieldSize"<<FEI_ENDL << fieldSize << FEI_ENDL
        << "#numNodes"<<FEI_ENDL << numNodes << FEI_ENDL
        << "#nodeIDs" << FEI_ENDL;
     int i;
     for(i=0; i<numNodes; ++i) os << (int)nodeIDs[i] << FEI_ENDL;
     os << "#estimates" << FEI_ENDL;
     for(i=0; i<numNodes*fieldSize; ++i) os << estimates[i] << FEI_ENDL;
   }

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
   if (!block->containsField(fieldID)) return(1);

   NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

   //if we have a negative fieldID, we'll need a list of length numNodes,
   //in which to put nodeNumbers for passing to the solver... 

   std::vector<int> numbers(numNodes);

   //if we have a fieldID >= 0, then our numbers list will hold equation numbers
   //and we'll need fieldSize*numNodes of them.

   std::vector<double> data;

   if (fieldID >= 0) {
      assert(fieldSize >= 0);
      numbers.resize(numNodes*fieldSize);
      data.resize(numNodes*fieldSize);
   }

   int count = 0;

   for(int i=0; i<numNodes; i++) {
     const NodeDescriptor* node = NULL;
     CHK_ERR( nodeDB.getNodeWithID(nodeIDs[i], node) );

      if (fieldID < 0) numbers[count++] = node->getNodeNumber();
      else {
         int eqn = -1;
         if (node->getFieldEqnNumber(fieldID, eqn)) {
           if (eqn >= localStartRow_ && eqn <= localEndRow_) {
             for(int j=0; j<fieldSize; j++) { 
               data[count] = estimates[i*fieldSize + j];
               problemStructure_->translateToReducedEqn(eqn+j, numbers[count++]);
             }
           }
         }
      }
   }

   if (fieldID < 0) {
     CHK_ERR( lsc_->putNodalFieldData(fieldID, fieldSize, 
                                      &numbers[0], numNodes, estimates));
   }
   else {
     CHK_ERR(lsc_->putInitialGuess(&numbers[0], &data[0], count));
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getBlockElemSolution(GlobalID elemBlockID,
                                   int numElems, 
                                   const GlobalID *elemIDs,
                                   int& numElemDOFPerElement,
                                   double *results)
{
//
//  return the elemental solution parameters associated with a
//  particular block of elements
//
   debugOutput("FEI: getBlockElemSolution");

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) )

   numElemDOFPerElement = block->getNumElemDOFPerElement();
   if (numElemDOFPerElement <= 0) return(0);

   ConnectivityTable& ctable =
     problemStructure_->getBlockConnectivity(elemBlockID);
   std::map<GlobalID,int>& elemIDList = ctable.elemIDs;

   std::vector<int>& elemDOFEqnNumbers = block->elemDOFEqnNumbers();
   double answer;

   for(int i=0; i<numElems; i++) {
      std::map<GlobalID,int>::const_iterator
        iter = elemIDList.find(elemIDs[i]);
      if (iter == elemIDList.end()) continue;
      int index = iter->second;

      int offset = i*numElemDOFPerElement;

      for(int j=0; j<numElemDOFPerElement; j++) {
         int eqn = elemDOFEqnNumbers[index] + j;

         CHK_ERR( getEqnSolnEntry(eqn, answer) )

         results[offset+j] = answer;
      }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putBlockElemSolution(GlobalID elemBlockID,
                                   int numElems,
                                   const GlobalID *elemIDs,
                                   int dofPerElem,
                                   const double *estimates)
{
   debugOutput("FEI: putBlockElemSolution");

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) )

   int DOFPerElement = block->getNumElemDOFPerElement();
   assert(DOFPerElement == dofPerElem);
   if (DOFPerElement <= 0) return(0);

   ConnectivityTable& ctable =
     problemStructure_->getBlockConnectivity(elemBlockID);
   std::map<GlobalID,int>& elemIDList = ctable.elemIDs;

   std::vector<int>& elemDOFEqnNumbers = block->elemDOFEqnNumbers();


   for(int i=0; i<numElems; i++) {
     std::map<GlobalID,int>::const_iterator
       iter = elemIDList.find(elemIDs[i]);
     if (iter == elemIDList.end()) continue;

     int index = iter->second;

      for(int j=0; j<DOFPerElement; j++) {
         int reducedEqn;
         problemStructure_->
           translateToReducedEqn(elemDOFEqnNumbers[index] + j, reducedEqn);
         double soln = estimates[i*DOFPerElement + j];

         CHK_ERR( lsc_->putInitialGuess(&reducedEqn, &soln, 1) );
      }
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::getCRMultipliers(int numCRs,
                                       const int* CRIDs,
                                       double* multipliers)
{
  int multCRsLen = problemStructure_->getNumMultConstRecords();
  if (numCRs > multCRsLen) {
    return(-1);
  }

  std::map<GlobalID, ConstraintType*>::const_iterator
    cr_iter = problemStructure_->getMultConstRecords().begin(),
    cr_end  = problemStructure_->getMultConstRecords().end();

  int i = 0;
  while(cr_iter != cr_end && i < numCRs) {
    GlobalID CRID = (*cr_iter).first;
    ConstraintType* multCR = (*cr_iter).second;
    if (CRID != CRIDs[i]) {
      CHK_ERR( problemStructure_->getMultConstRecord(CRIDs[i], multCR) );
    }

    int eqn = multCR->getEqnNumber();

    CHK_ERR( getEqnSolnEntry(eqn, multipliers[i]) );
    ++cr_iter;
    ++i;
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putCRMultipliers(int numMultCRs,
                               const int* CRIDs,
                               const double *multEstimates)
{
  debugOutput("FEI: putCRMultipliers");

  for(int j = 0; j < numMultCRs; j++) {
    ConstraintType* multCR = NULL;
    CHK_ERR( problemStructure_->getMultConstRecord(CRIDs[j], multCR) );

    int eqnNumber = multCR->getEqnNumber();
    if (eqnNumber < localStartRow_ || eqnNumber > localEndRow_) continue;
    CHK_ERR( lsc_->putInitialGuess(&eqnNumber, &(multEstimates[j]), 1));
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putNodalFieldData(int fieldID,
                                        int numNodes,
                                        const GlobalID* nodeIDs,
                                        const double* nodeData)
{
  debugOutput("FEI: putNodalFieldData");

  if (fieldID > -1) {
    return(putNodalFieldSolution(fieldID, numNodes, nodeIDs, nodeData));
  }

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  std::vector<int> nodeNumbers(numNodes);

  for(int i=0; i<numNodes; i++) {
    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithID(nodeIDs[i], node) );

    int nodeNumber = node->getNodeNumber();
    if (nodeNumber < 0) {
      fei::console_out() << "LinSysCoreFilter::putNodalFieldData ERROR, node with ID " 
           << (int)nodeIDs[i] << " doesn't have an associated nodeNumber "
           << "assigned. putNodalFieldData shouldn't be called until after the "
           << "initComplete method has been called." << FEI_ENDL;
      ERReturn(-1);
    }

    nodeNumbers[i] = nodeNumber;
  }

  CHK_ERR( lsc_->putNodalFieldData(fieldID, fieldSize,
                                   &nodeNumbers[0], numNodes, nodeData) );

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::putNodalFieldSolution(int fieldID,
                                int numNodes,
                                const GlobalID* nodeIDs,
                                const double* nodeData)
{
  debugOutput("FEI: putNodalFieldSolution");

  if (fieldID < 0) {
    return(putNodalFieldData(fieldID, numNodes, nodeIDs, nodeData));
  }

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  std::vector<int> eqnNumbers(fieldSize);

  for(int i=0; i<numNodes; i++) {
    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithID(nodeIDs[i], node) );

    int eqn = -1;
    bool hasField = node->getFieldEqnNumber(fieldID, eqn);
    if (!hasField) continue;

    int reducedEqn = -1;
    bool isSlave = problemStructure_->translateToReducedEqn(eqn, reducedEqn);
    if (isSlave) continue;

    if (reducedStartRow_ > reducedEqn || reducedEndRow_ < reducedEqn) continue;

    int localLen = fieldSize;
    for(int j=0; j<fieldSize; j++) {
      int thisEqn = reducedEqn+j;
      if (reducedStartRow_ > thisEqn || reducedEndRow_ <thisEqn) {
        localLen = j;
      }

      eqnNumbers[j] = reducedEqn+j;
    }

    int offset = i*fieldSize;
    CHK_ERR( lsc_->putInitialGuess(&eqnNumbers[0],
                                   &nodeData[offset], localLen) );
  }

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::assembleEqns(int numPtRows, 
                                   int numPtCols,
                                   const int* rowNumbers,
                                   const int* colIndices,
                                   const double* const* coefs,
                                   bool structurallySymmetric,
                                   int numBlkEqns, int* blkEqns,
                                   int* blkSizes, bool useBlkEqns,
                                   int mode)
{
  if (numPtRows == 0) return(FEI_SUCCESS);

  bool anySlaves = false;
  int numSlaveEqns = problemStructure_->numSlaveEquations();
  if (numSlaveEqns > 0) {
    rSlave_.resize(numPtRows);
    cSlave_.resize(0);
    const int* indPtr = colIndices;
    for(int r=0; r<numPtRows; r++) {
      rSlave_[r] = problemStructure_->isSlaveEqn(rowNumbers[r]) ? 1 : 0;
      if (rSlave_[r] == 1) anySlaves = true;

      for(int j=0; j<numPtCols; j++) {
        int isSlave = problemStructure_->isSlaveEqn(indPtr[j]) ? 1 : 0;
        cSlave_.push_back(isSlave);
        if (isSlave == 1) anySlaves = true;
      }

      if (!structurallySymmetric) indPtr += numPtCols;
    }
  }

  if (numSlaveEqns == 0 || !anySlaves) {
    if (numSlaveEqns == 0 && structurallySymmetric) {
      if (useBlkEqns) {
        CHK_ERR( giveToBlkMatrix_symm_noSlaves(numPtRows, rowNumbers,
                                               numBlkEqns, blkEqns, blkSizes,
                                               coefs, mode) );
      }
      else {
        CHK_ERR( giveToMatrix_symm_noSlaves(numPtRows, rowNumbers, coefs, mode) );
      }
    }
    else {
      if ((int)dworkSpace2_.size() < numPtRows) {
        dworkSpace2_.resize(numPtRows);
      }
      const double* * coefPtr = &dworkSpace2_[0];
      for(int i=0; i<numPtRows; i++) {
        coefPtr[i] = coefs[i];
      }

      if (structurallySymmetric) {
        CHK_ERR( giveToMatrix(numPtRows, rowNumbers, numPtRows, rowNumbers,
                              coefPtr, mode) );
      }
      else {
        const int* indPtr = colIndices;
        for(int i=0; i<numPtRows; i++) {
          int row = rowNumbers[i];

          const double* coefPtr1 = coefs[i];

          CHK_ERR(giveToMatrix(1, &row, numPtCols, indPtr, &coefPtr1, mode));
          indPtr += numPtCols;
        }
      }
    }
  }
  else {
    int offset = 0;
    const int* indicesPtr = colIndices;
    for(int i=0; i<numPtRows; i++) {
      int row = rowNumbers[i];

      const double* coefPtr = coefs[i];
      int* colSlave = &cSlave_[offset];
      offset += numPtCols;

      if (rSlave_[i] == 1) {
        //Since this is a slave equation, the non-slave columns of this row go
        //into 'Kdi_', and the slave columns go into 'Kdd_'.
        for(int jj=0; jj<numPtCols; jj++) {
          int col = indicesPtr[jj];
          if (colSlave[jj]) {
            Kdd_->sumInCoef(row, col, coefPtr[jj]);
          }
          else {
            Kdi_->sumInCoef(row, col, coefPtr[jj]);
          }
        }

        //We also need to put the non-slave rows of column 'row' into 'K_id'.
        const int* ii_indicesPtr = colIndices;
        for(int ii=0; ii<numPtRows; ii++) {
          int rowi = rowNumbers[ii];
          if (rSlave_[ii] == 1) continue;

          int index = fei::binarySearch(row, ii_indicesPtr, numPtCols);
          if (index < 0) continue;

          const double* coefs_ii = coefs[ii];

          Kid_->sumInCoef(rowi, row, coefs_ii[index]);

          if (!structurallySymmetric) ii_indicesPtr += numPtCols;
        }

        reducedEqnCounter_++;

        continue;
      }
      else {//row is not a slave eqn...

        //put all non-slave columns from this row into the assembled matrix.

        CHK_ERR( giveToMatrix(1, &row, numPtCols, indicesPtr, &coefPtr, mode) );
      }

      if (!structurallySymmetric) indicesPtr += numPtCols;
    }

    if (reducedEqnCounter_ > 300) CHK_ERR( assembleReducedEqns() );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::assembleReducedEqns()
{
  fei::FillableMat* D = problemStructure_->getSlaveDependencies();

  csrD = *D;
  csrKid = *Kid_;
  csrKdi = *Kdi_;
  csrKdd = *Kdd_;

  //form tmpMat1_ = Kid_*D
  fei::multiply_CSRMat_CSRMat(csrKid, csrD, tmpMat1_);

  //form tmpMat2_ = D^T*Kdi_
  fei::multiply_trans_CSRMat_CSRMat(csrD, csrKdi, tmpMat2_);

  if (Filter::logStream() != NULL) {
    FEI_OSTREAM& os = *Filter::logStream();
    os << "#  tmpMat1_"<<FEI_ENDL << tmpMat1_ << FEI_ENDL;
    os << "#  tmpMat2_"<<FEI_ENDL << tmpMat2_ << FEI_ENDL;
  }

  //accumulate the above two results into the global system matrix.
  CHK_ERR( sumIntoMatrix(tmpMat1_) );
  CHK_ERR( sumIntoMatrix(tmpMat2_) );

  //form tmpMat1_ = D^T*Kdd_
  fei::multiply_trans_CSRMat_CSRMat(csrD, csrKdd, tmpMat1_);

  //form tmpMat2_ = tmpMat1_*D = D^T*Kdd_*D
  fei::multiply_CSRMat_CSRMat(tmpMat1_, csrD, tmpMat2_);

  if (Filter::logStream() != NULL) {
    FEI_OSTREAM& os = *Filter::logStream();
    os << "#  tmpMat2_"<<FEI_ENDL << tmpMat2_ << FEI_ENDL;
  }

  //finally, accumulate tmpMat2_ = D^T*Kdd_*D into the global system matrix.
  CHK_ERR( sumIntoMatrix(tmpMat2_) );

  Kdi_->clear();
  Kid_->clear();
  Kdd_->clear();
  reducedEqnCounter_ = 0;

  return(0);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::assembleRHS(int numValues,
                                  const int* indices,
                                  const double* coefs,
                                  int mode) {
//
//This function hands the data off to the routine that finally
//sticks it into the RHS vector.
//

  if (problemStructure_->numSlaveEquations() == 0) {
    CHK_ERR( giveToRHS(numValues, coefs, indices, mode) );
    return(FEI_SUCCESS);
  }

  for(int i = 0; i < numValues; i++) {
    int eqn = indices[i];

    if (problemStructure_->isSlaveEqn(eqn)) {
      fei::add_entry( fd_, eqn, coefs[i]);
      reducedRHSCounter_++;
      continue;
    }

    CHK_ERR( giveToRHS(1, &(coefs[i]), &eqn, mode ) );
  }

  if (reducedRHSCounter_ > 300) CHK_ERR( assembleReducedRHS() );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int LinSysCoreFilter::assembleReducedRHS()
{
  fei::FillableMat* D = problemStructure_->getSlaveDependencies();

  csrD = *D;

  //now form tmpVec1_ = D^T*fd_.
  fei::multiply_trans_CSRMat_CSVec(csrD, fd_, tmpVec1_);

  CHK_ERR( sumIntoRHS(tmpVec1_) );

  fd_.clear();
  reducedRHSCounter_ = 0;

  return(0);
}

//==============================================================================
void LinSysCoreFilter::debugOutput(const char* mesg) {
   if (Filter::logStream() != NULL) {
     (*logStream())<<mesg<<FEI_ENDL;
   }
}

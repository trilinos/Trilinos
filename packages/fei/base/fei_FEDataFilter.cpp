/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <limits>
#include <cmath>

#include <fei_defs.h>

#include <fei_CommUtils.hpp>
#include <fei_TemplateUtils.hpp>
#include <snl_fei_Constraint.hpp>
typedef snl_fei::Constraint<GlobalID> ConstraintType;

#include <fei_LibraryWrapper.hpp>
#include <SNL_FEI_Structure.hpp>
#include <fei_FiniteElementData.hpp>
#include <fei_Lookup.hpp>
#include <FEI_Implementation.hpp>
#include <fei_EqnCommMgr.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_NodeDatabase.hpp>
#include <fei_NodeCommMgr.hpp>
#include <fei_ProcEqns.hpp>
#include <fei_BlockDescriptor.hpp>
#include <fei_ConnectivityTable.hpp>
#include <snl_fei_Utils.hpp>

#include <fei_FEDataFilter.hpp>

#undef fei_file
#define fei_file "FEDataFilter.cpp"
#include <fei_ErrMacros.hpp>

#define ASSEMBLE_PUT 0
#define ASSEMBLE_SUM 1

//------------------------------------------------------------------------------
void convert_eqns_to_nodenumbers_and_dof_ids(fei::FieldDofMap<int>& fdmap,
                                             const NodeDatabase& nodeDB,
                                             int numEqns,
                                             const int* eqns,
                                             std::vector<int>& nodeNumbers,
                                             std::vector<int>& dof_ids)
{
  nodeNumbers.resize(numEqns);
  dof_ids.resize(numEqns);

  for(int i=0; i<numEqns; ++i) {
    const NodeDescriptor* nodePtr = NULL;
    int err = nodeDB.getNodeWithEqn(eqns[i], nodePtr);
    if (err < 0) {
      nodeNumbers[i] = -1;
      dof_ids[i] = -1;
      continue;
    }

    nodeNumbers[i] = nodePtr->getNodeNumber();

    int fieldID, offset;
    nodePtr->getFieldID(eqns[i], fieldID, offset);
    dof_ids[i] = fdmap.get_dof_id(fieldID, offset);
  }
}

//------------------------------------------------------------------------------
void convert_field_and_nodes_to_eqns(const NodeDatabase& nodeDB,
                                     int fieldID, int fieldSize,
                                     int numNodes, const GlobalID* nodeIDs,
                                     std::vector<int>& eqns)
{
  eqns.assign(numNodes*fieldSize, -1);

  size_t offset = 0;
  for(int i=0; i<numNodes; ++i) {
    const NodeDescriptor* node = NULL;
    int err = nodeDB.getNodeWithID(nodeIDs[i], node);
    if (err < 0) {
      offset += fieldSize;
      continue;
    }

    int eqn = 0;
    node->getFieldEqnNumber(fieldID, eqn);
    for(int j=0; j<fieldSize; ++j) {
      eqns[offset++] = eqn+j;
    }
  }
}

//------------------------------------------------------------------------------
FEDataFilter::FEDataFilter(FEI_Implementation* owner,
                           MPI_Comm comm,
                           SNL_FEI_Structure* probStruct,
                           LibraryWrapper* wrapper,
                           int masterRank)
 : Filter(probStruct),
   wrapper_(wrapper),
   feData_(NULL),
   useLookup_(true),
   internalFei_(0),
   newData_(false),
   localStartRow_(0),
   localEndRow_(0),
   numGlobalEqns_(0),
   reducedStartRow_(0),
   reducedEndRow_(0),
   numReducedRows_(0),
   iterations_(0),
   numRHSs_(0),
   currentRHS_(0),
   rhsIDs_(),
   outputLevel_(0),
   comm_(comm),
   masterRank_(masterRank),
   problemStructure_(probStruct),
   penCRIDs_(),
   rowIndices_(),
   rowColOffsets_(0),
   colIndices_(0),
   eqnCommMgr_(NULL),
   eqnCommMgr_put_(NULL),
   maxElemRows_(0),
   eStiff_(NULL),
   eStiff1D_(NULL),
   eLoad_(NULL),
   numRegularElems_(0),
   constraintBlocks_(0, 16),
   constraintNodeOffsets_(),
   packedFieldSizes_()
{
  localRank_ = fei::localProc(comm_);
  numProcs_ = fei::numProcs(comm_);

  internalFei_ = 0;

  numRHSs_ = 1;
  rhsIDs_.resize(numRHSs_);
  rhsIDs_[0] = 0;

  eqnCommMgr_ = problemStructure_->getEqnCommMgr().deepCopy();
  createEqnCommMgr_put();

  if (wrapper_->haveFiniteElementData()) {
    feData_ = wrapper_->getFiniteElementData();
  }
  else {
    fei::console_out() << "FEDataFilter::FEDataFilter ERROR, must be constructed with a "
         << "FiniteElementData interface. Aborting." << FEI_ENDL;
#ifndef FEI_SER
    MPI_Abort(comm_, -1);
#else
    abort();
#endif
  }

  //We need to get the parameters from the owning FEI_Implementation, if we've
  //been given a non-NULL FEI_Implementation...
  if (owner != NULL) {
    int numParams = -1;
    char** paramStrings = NULL;
    int err = owner->getParameters(numParams, paramStrings);

    //Now let's pass them into our own parameter-handling mechanism.
    err = parameters(numParams, paramStrings);
    if (err != 0) {
      fei::console_out() << "FEDataFilter::FEDataFilter ERROR, parameters failed." << FEI_ENDL;
      MPI_Abort(comm_, -1);
    }
  }

  return;
}

//------------------------------------------------------------------------------
FEDataFilter::FEDataFilter(const FEDataFilter& src)
 : Filter(NULL),
   wrapper_(NULL),
   feData_(NULL),
   useLookup_(true),
   internalFei_(0),
   newData_(false),
   localStartRow_(0),
   localEndRow_(0),
   numGlobalEqns_(0),
   reducedStartRow_(0),
   reducedEndRow_(0),
   numReducedRows_(0),
   iterations_(0),
   numRHSs_(0),
   currentRHS_(0),
   rhsIDs_(),
   outputLevel_(0),
   comm_(0),
   masterRank_(0),
   problemStructure_(NULL),
   penCRIDs_(),
   rowIndices_(),
   rowColOffsets_(0),
   colIndices_(0),
   eqnCommMgr_(NULL),
   eqnCommMgr_put_(NULL),
   maxElemRows_(0),
   eStiff_(NULL),
   eStiff1D_(NULL),
   eLoad_(NULL),
   numRegularElems_(0),
   constraintBlocks_(0, 16),
   constraintNodeOffsets_(),
   packedFieldSizes_()
{
}

//------------------------------------------------------------------------------
FEDataFilter::~FEDataFilter() {
//
//  Destructor function. Free allocated memory, etc.
//
  numRHSs_ = 0;

  delete eqnCommMgr_;
  delete eqnCommMgr_put_;

  delete [] eStiff_;
  delete [] eStiff1D_;
  delete [] eLoad_;
}

//------------------------------------------------------------------------------
int FEDataFilter::initialize()
{
// Determine final sparsity pattern for setting the structure of the
// underlying sparse matrix.
//
  debugOutput("#  initialize");

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

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::createEqnCommMgr_put()
{
  if (eqnCommMgr_put_ != NULL) return(0);

  eqnCommMgr_put_  = eqnCommMgr_->deepCopy();
  if (eqnCommMgr_put_ == NULL) ERReturn(-1);

  eqnCommMgr_put_->resetCoefs();
  eqnCommMgr_put_->accumulate_ = false;
  return(0);
}

//==============================================================================
int FEDataFilter::initLinSysCore()
{
  try {

  int err = wrapper_->getFiniteElementData()->setLookup(*problemStructure_);

  if (err != 0) {
    useLookup_ = false;
  }

  reducedStartRow_ = localStartRow_;
  reducedEndRow_ = localEndRow_;

  int numElemBlocks = problemStructure_->getNumElemBlocks();
  NodeDatabase& nodeDB     = problemStructure_->getNodeDatabase();
  NodeCommMgr& nodeCommMgr = problemStructure_->getNodeCommMgr();

  int numNodes = nodeDB.getNumNodeDescriptors();
  int numRemoteNodes = nodeCommMgr.getSharedNodeIDs().size() -
                         nodeCommMgr.getLocalNodeIDs().size();
  numNodes -= numRemoteNodes;

  int numSharedNodes = nodeCommMgr.getNumSharedNodes();

  std::vector<int> numElemsPerBlock(numElemBlocks);
  std::vector<int> numNodesPerElem(numElemBlocks);
  std::vector<int> elemMatrixSizePerBlock(numElemBlocks);

  for(int blk=0; blk<numElemBlocks; blk++) {
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor_index(blk, block) );

    numElemsPerBlock[blk] = block->getNumElements();
    numNodesPerElem[blk]  = block->getNumNodesPerElement();

    int* fieldsPerNode = block->fieldsPerNodePtr();
    int** fieldIDsTable = block->fieldIDsTablePtr();

    elemMatrixSizePerBlock[blk] = 0;

    for(int nn=0; nn<numNodesPerElem[blk]; nn++) {
      if (fieldsPerNode[nn] <= 0) ERReturn(-1);
      
      for(int nf=0; nf<fieldsPerNode[nn]; nf++) {
        elemMatrixSizePerBlock[blk] +=
          problemStructure_->getFieldSize(fieldIDsTable[nn][nf]);
      }
    }
  }

  //Now we need to run the penalty constraint records and figure out how many
  //extra "element-blocks" to describe. (A penalty constraint will be treated 
  //exactly like an element.) So first, we need to figure out how many different
  //sizes of constraint connectivities there are, because the constraints with
  //the same numbers of constrained nodes will be grouped together in blocks.

  if (problemStructure_==NULL) {
    FEI_COUT << "problemStructrue_ NULL"<<FEI_ENDL;
    ERReturn(-1);
  }

  std::map<GlobalID,ConstraintType*>::const_iterator
    cr_iter = problemStructure_->getPenConstRecords().begin(),
    cr_end  = problemStructure_->getPenConstRecords().end();

  //constraintBlocks will be a sorted list with each "block-id" being the
  //num-nodes-per-constraint for constraints in that block.

  //numConstraintsPerBlock is the same length as constraintBlocks
  std::vector<int> numConstraintsPerBlock;
  std::vector<int> numDofPerConstraint;
  penCRIDs_.resize(problemStructure_->getNumPenConstRecords());

  int counter = 0;
  while(cr_iter != cr_end) {
    penCRIDs_[counter++] = (*cr_iter).first;
    ConstraintType& cr = *((*cr_iter).second);
    int nNodes = cr.getMasters().size();

    int insertPoint = -1;
    int offset = fei::binarySearch(nNodes, constraintBlocks_, insertPoint);

    int nodeOffset = 0;
    if (offset < 0) {
      constraintBlocks_.insert(constraintBlocks_.begin()+insertPoint, nNodes);
      numConstraintsPerBlock.insert(numConstraintsPerBlock.begin()+insertPoint, 1);
      numDofPerConstraint.insert(numDofPerConstraint.begin()+insertPoint, 0);

      if (insertPoint > 0) {
        nodeOffset = constraintNodeOffsets_[insertPoint-1] +
           constraintBlocks_[insertPoint-1];
      }
      constraintNodeOffsets_.insert(constraintNodeOffsets_.begin()+insertPoint, nodeOffset);
      offset = insertPoint;
    }
    else {
      numConstraintsPerBlock[offset]++;
      ++cr_iter;
      continue;
    }

    std::vector<int>& fieldIDsvec = cr.getMasterFieldIDs();
    int* fieldIDs = &fieldIDsvec[0];
    for(int k=0; k<nNodes; ++k) {
      int fieldSize = problemStructure_->getFieldSize(fieldIDs[k]);
      packedFieldSizes_.insert(packedFieldSizes_.begin()+nodeOffset+k, fieldSize);
      numDofPerConstraint[offset] += fieldSize;
    }
    ++cr_iter;
  }

  //now combine the elem-block info with the penalty-constraint info.
  int numBlocksTotal = numElemBlocks + constraintBlocks_.size();
  for(size_t i=0; i<constraintBlocks_.size(); ++i) {
    numElemsPerBlock.push_back(numConstraintsPerBlock[i]);
    numNodesPerElem.push_back(constraintBlocks_[i]);
    elemMatrixSizePerBlock.push_back(numDofPerConstraint[i]);
  }

  int numMultCRs = problemStructure_->getNumMultConstRecords();

  CHK_ERR( feData_->describeStructure(numBlocksTotal,
                                      &numElemsPerBlock[0],
                                      &numNodesPerElem[0],
                                      &elemMatrixSizePerBlock[0],
                                      numNodes,
                                      numSharedNodes,
                                      numMultCRs) );

  numRegularElems_ = 0;
  std::vector<int> numDofPerNode;
  std::vector<int> dof_ids;
  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  for(int i=0; i<numElemBlocks; i++) {
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor_index(i, block) );

    if (block->getNumElements() == 0) continue;

    ConnectivityTable& ctbl =
      problemStructure_->getBlockConnectivity(block->getGlobalBlockID());

    std::vector<int> cNodeList(block->getNumNodesPerElement());

    int* fieldsPerNode = block->fieldsPerNodePtr();
    int** fieldIDsTable = block->fieldIDsTablePtr();

    numDofPerNode.resize(0);
    int total_num_dof = 0;
    for(int nn=0; nn<numNodesPerElem[i]; nn++) {
      if (fieldsPerNode[nn] <= 0) ERReturn(-1);
      numDofPerNode.push_back(0);
      int indx = numDofPerNode.size()-1;

      for(int nf=0; nf<fieldsPerNode[nn]; nf++) {
        numDofPerNode[indx] += problemStructure_->getFieldSize(fieldIDsTable[nn][nf]);
      }
      total_num_dof += numDofPerNode[indx];
    }

    dof_ids.resize(total_num_dof);
    int doffset = 0;
    for(int nn=0; nn<numNodesPerElem[i]; ++nn) {
      for(int nf=0; nf<fieldsPerNode[nn]; ++nf) {
        int fieldSize = problemStructure_->getFieldSize(fieldIDsTable[nn][nf]);
        for(int dof_offset=0; dof_offset<fieldSize; ++dof_offset) {
          dof_ids[doffset++] = fdmap.get_dof_id(fieldIDsTable[nn][nf], dof_offset);
        }
      }
    }

    int nodesPerElement = block->getNumNodesPerElement();
    NodeDescriptor** elemConn = &((*ctbl.elem_conn_ptrs)[0]);
    int offset = 0;
    int numElems = block->getNumElements();
    numRegularElems_ += numElems;
    for(int j=0; j<numElems; j++) {

      for(int k=0; k<nodesPerElement; k++) {
        NodeDescriptor* node = elemConn[offset++];
        cNodeList[k] = node->getNodeNumber();
      }

      CHK_ERR( feData_->setConnectivity(i, ctbl.elemNumbers[j],
                                        block->getNumNodesPerElement(),
                                        &cNodeList[0],
                                        &numDofPerNode[0],
                                        &dof_ids[0]) );
    }
  }

  std::vector<int> nodeNumbers;
  cr_iter = problemStructure_->getPenConstRecords().begin();
  int i = 0;
  while(cr_iter != cr_end) {
    ConstraintType& cr = *((*cr_iter).second);
    std::vector<GlobalID>& nodeIDsvec = cr.getMasters();
    GlobalID* nodeIDs = &nodeIDsvec[0];
    int nNodes = cr.getMasters().size();
    int index = fei::binarySearch(nNodes, constraintBlocks_);
    if (index < 0) {
      ERReturn(-1);
    }

    int total_num_dof = 0;
    std::vector<int>& masterFieldIDs = cr.getMasterFieldIDs();
    for(size_t k=0; k<masterFieldIDs.size(); ++k) {
      total_num_dof += problemStructure_->getFieldSize(masterFieldIDs[k]);
    }

    dof_ids.resize(total_num_dof);
    int doffset = 0;
    for(size_t k=0; k<masterFieldIDs.size(); ++k) {
      int field_size = problemStructure_->getFieldSize(masterFieldIDs[k]);
      for(int dof_offset=0; dof_offset<field_size; ++dof_offset) {
        dof_ids[doffset++] = fdmap.get_dof_id(masterFieldIDs[k], dof_offset);
      }
    }

    int blockNum = numElemBlocks + index;

    nodeNumbers.resize(nNodes);

    for(int k=0; k<nNodes; ++k) {
      const NodeDescriptor* node = Filter::findNode(nodeIDs[k]);
      if(node == NULL)
      {
        nodeNumbers[k] = -1;
      }
      else
      {
        nodeNumbers[k] = node->getNodeNumber();
      }
    }

    int offset = constraintNodeOffsets_[index];
    CHK_ERR( feData_->setConnectivity(blockNum, numRegularElems_+i++, nNodes, &nodeNumbers[0], &packedFieldSizes_[offset], &dof_ids[0]) );
    ++cr_iter;
  }

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::resetSystem(double s)
{
  //
  //  This puts the value s throughout both the matrix and the vector.
  //
  if (Filter::logStream() != NULL) {
    (*logStream()) << "FEI: resetSystem" << FEI_ENDL << s << FEI_ENDL;
  }

  CHK_ERR( feData_->reset() );
 
  debugOutput("#FEDataFilter leaving resetSystem");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::deleteMultCRs()
{

  debugOutput("#FEDataFilter::deleteMultCRs");

  int err = feData_->deleteConstraints();

  debugOutput("#FEDataFilter leaving deleteMultCRs");

  return(err);
}

//------------------------------------------------------------------------------
int FEDataFilter::resetTheMatrix(double s)
{
  //FiniteElementData implementations can't currently reset the matrix without
  //resetting the rhs vector too. 
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::resetTheRHSVector(double s)
{
  //FiniteElementData implementations can't currently reset the rhs vector
  //without resetting the matrix too.
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::resetMatrix(double s)
{
  //
  //  This puts the value s throughout both the matrix and the vector.
  //

  debugOutput("FEI: resetMatrix");

  CHK_ERR( resetTheMatrix(s) );

  eqnCommMgr_->resetCoefs();

  debugOutput("#FEDataFilter leaving resetMatrix");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::resetRHSVector(double s)
{
  //
  //  This puts the value s throughout the rhs vector.
  //

  debugOutput("FEI: resetRHSVector");

  CHK_ERR( resetTheRHSVector(s) );

  eqnCommMgr_->resetCoefs();

  debugOutput("#FEDataFilter leaving resetRHSVector");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::resetInitialGuess(double s)
{
  //
  //  This puts the value s throughout the initial guess (solution) vector.
  //
  if (Filter::logStream() != NULL) {
    FEI_OSTREAM& os = *logStream();
    os << "FEI: resetInitialGuess" << FEI_ENDL;
    os << "#value to which initial guess is to be set" << FEI_ENDL;
    os << s << FEI_ENDL;
  }

  //Actually, the FiniteElementData doesn't currently allow us to alter
  //values in any initial guess or solution vector.

  debugOutput("#FEDataFilter leaving resetInitialGuess");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadNodeBCs(int numNodes,
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
      int nodeID = nodeIDs[j];
      (*logStream())<<nodeID<<"  ";
      (*logStream())<< offsetsIntoField[j]<<" ";
      (*logStream())<< prescribedValues[j]<<FEI_ENDL;
    }
   }

   NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

   std::vector<int> essEqns(numNodes);
   std::vector<double> alpha(numNodes);
   std::vector<double> gamma(numNodes);

   for(int i=0; i<numNodes; ++i) {
     NodeDescriptor* node = NULL;
     nodeDB.getNodeWithID(nodeIDs[i], node);
     if (node == NULL) {
       fei::console_out() << "fei_FEDataFilter::loadNodeBCs ERROR, node " << nodeIDs[i]
           << " not found." << FEI_ENDL;
       ERReturn(-1);
     }

     int eqn = -1;
     if (!node->getFieldEqnNumber(fieldID, eqn)) {
       ERReturn(-1);
     }

     essEqns[i] = eqn + offsetsIntoField[i];
     gamma[i] = prescribedValues[i];
     alpha[i] = 1.0;
   }

   if (essEqns.size() > 0) {
      CHK_ERR( enforceEssentialBCs(&essEqns[0], &alpha[0],
                                   &gamma[0], essEqns.size()) );
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadElemBCs(int numElems,
                              const GlobalID *elemIDs,
                              int fieldID,
                              const double *const *alpha,
                              const double *const *beta,
                              const double *const *gamma)
{
   return(-1);
}

//------------------------------------------------------------------------------
void FEDataFilter::allocElemStuff()
{
   int nb = problemStructure_->getNumElemBlocks();

   for(int i=0; i<nb; i++) {
     BlockDescriptor* block = NULL;
     int err = problemStructure_->getBlockDescriptor_index(i, block);
     if (err) voidERReturn;

      int numEqns = block->getNumEqnsPerElement();
      if (maxElemRows_ < numEqns) maxElemRows_ = numEqns;
   }

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
int FEDataFilter::sumInElem(GlobalID elemBlockID,
                        GlobalID elemID,
                        const GlobalID* elemConn,
                        const double* const* elemStiffness,
                        const double* elemLoad,
                        int elemFormat)
{
  if (Filter::logStream() != NULL) {
    (*logStream()) << "FEI: sumInElem" << FEI_ENDL <<"# elemBlockID " << FEI_ENDL
                      << static_cast<int>(elemBlockID) << FEI_ENDL
                      << "# elemID " << FEI_ENDL << static_cast<int>(elemID) << FEI_ENDL;
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
    int numNodes = block->getNumNodesPerElement();
    (*logStream()) << "#num-nodes" << FEI_ENDL << numNodes << FEI_ENDL;
    (*logStream()) << "#connected nodes" << FEI_ENDL;
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
int FEDataFilter::sumInElemMatrix(GlobalID elemBlockID,
                              GlobalID elemID,
                              const GlobalID* elemConn,
                              const double* const* elemStiffness,
                              int elemFormat)
{
  if (Filter::logStream() != NULL) {
    (*logStream()) << "FEI: sumInElemMatrix"<<FEI_ENDL
                      << "#elemBlockID" << FEI_ENDL << static_cast<int>(elemBlockID)
                      << "# elemID" << FEI_ENDL << static_cast<int>(elemID) << FEI_ENDL;
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
    int numNodes = block->getNumNodesPerElement();
    (*logStream()) << "#num-nodes" << FEI_ENDL << numNodes << FEI_ENDL;
    (*logStream()) << "#connected nodes" << FEI_ENDL;
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
int FEDataFilter::sumInElemRHS(GlobalID elemBlockID,
                           GlobalID elemID,
                           const GlobalID* elemConn,
                           const double* elemLoad)
{
  if (Filter::logStream() != NULL) {
    (*logStream()) << "FEI: sumInElemRHS"<<FEI_ENDL<<"# elemBlockID " << FEI_ENDL
                      <<static_cast<int>(elemBlockID)
                      << "# elemID " << FEI_ENDL << static_cast<int>(elemID) << FEI_ENDL;
    BlockDescriptor* block = NULL;
    CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
    int numNodes = block->getNumNodesPerElement();
    (*logStream()) << "#num-nodes" << FEI_ENDL << numNodes << FEI_ENDL;
    (*logStream()) << "#connected nodes" << FEI_ENDL;
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
int FEDataFilter::generalElemInput(GlobalID elemBlockID,
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
int FEDataFilter::generalElemInput(GlobalID elemBlockID,
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

  //an std::vector.resize operation is free if the size is either shrinking or
  //staying the same.

  const double* const* stiff = NULL;
  if (elemStiffness != NULL) stiff = elemStiffness;

  const double* load = NULL;
  if (elemLoad != NULL) load = elemLoad;

  //we'll make a local dense copy of the element stiffness array
  //if the stiffness array was passed in using one of the "weird"
  //element formats, AND if the stiffness array is non-null.
  if (elemFormat != FEI_DENSE_ROW && stiff != NULL) {
    Filter::copyStiffness(stiff, numElemRows, elemFormat, eStiff_);
    stiff = eStiff_;
  }

  if (stiff != NULL || load != NULL) newData_ = true;

  if (Filter::logStream() != NULL) {
    if (stiff != NULL) {
      (*logStream())
        << "#numElemRows"<< FEI_ENDL << numElemRows << FEI_ENDL
        << "#elem-stiff (after being copied into dense-row format)"
        << FEI_ENDL;
      for(int i=0; i<numElemRows; i++) {
        const double* stiff_i = stiff[i];
        for(int j=0; j<numElemRows; j++) {
          (*logStream()) << stiff_i[j] << " ";
        }
        (*logStream()) << FEI_ENDL;
      }
    }

    if (load != NULL) {
      (*logStream()) << "#elem-load" << FEI_ENDL;
      for(int i=0; i<numElemRows; i++) {
        (*logStream()) << load[i] << " ";
      }
      (*logStream())<<FEI_ENDL;
    }
  }

  //Now we'll proceed to gather the stuff we need to pass the stiffness
  //data through to the FiniteElementData interface...

  int blockNumber = problemStructure_->getIndexOfBlock(elemBlockID);

  ConnectivityTable& connTable = problemStructure_->
    getBlockConnectivity(elemBlockID);

  std::map<GlobalID,int>::iterator
    iter = connTable.elemIDs.find(elemID);
  if (iter == connTable.elemIDs.end()) {
    ERReturn(-1);
  }

  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  int elemIndex = iter->second;

  int elemNumber = connTable.elemNumbers[elemIndex];

  int numNodes = block->getNumNodesPerElement();
  int* fieldsPerNode = block->fieldsPerNodePtr();
  int** fieldIDsTable = block->fieldIDsTablePtr();

  int numDistinctFields = block->getNumDistinctFields();
  int dof_id = 0;
  int fieldSize = 0;
  int total_num_dofs = 0;
  if (numDistinctFields == 1) {
    fieldSize = problemStructure_->getFieldSize(fieldIDsTable[0][0]);
    for(int i=0; i<numNodes; ++i) {
      total_num_dofs += fieldSize*fieldsPerNode[i];
    }
    dof_id = fdmap.get_dof_id(fieldIDsTable[0][0], 0);
  }
  else {
    for(int i=0; i<numNodes; ++i) {
      for(int nf=0; nf<fieldsPerNode[i]; ++nf) {
        total_num_dofs += problemStructure_->getFieldSize(fieldIDsTable[i][nf]);
      }
    }
  }

  static std::vector<int> iwork;
  iwork.resize(2*numNodes+total_num_dofs);

  int* dofsPerNode = &iwork[0];
  int* nodeNumbers = dofsPerNode+numNodes;
  int* dof_ids = nodeNumbers+numNodes;

  for(int i=0; i<numNodes; ++i) {
    dofsPerNode[i] = 0;
  }


  NodeDescriptor** elemNodes =
    &((*connTable.elem_conn_ptrs)[elemIndex*numNodes]);

  int doffset = 0;
  for(int nn=0; nn<numNodes; nn++) {
    NodeDescriptor* node = elemNodes[nn];
    nodeNumbers[nn] = node->getNodeNumber();

    if (numDistinctFields == 1) {
      for(int nf=0; nf<fieldsPerNode[nn]; nf++) {
        dofsPerNode[nn] += fieldSize;
        for(int dof_offset=0; dof_offset<fieldSize; ++dof_offset) {
          dof_ids[doffset++] = dof_id;
        }
      }
    }
    else {
      for(int nf=0; nf<fieldsPerNode[nn]; nf++) {
        int fieldSize = problemStructure_->getFieldSize(fieldIDsTable[nn][nf]);
        int dof_id = fdmap.get_dof_id(fieldIDsTable[nn][nf], 0);
        dofsPerNode[nn] += fieldSize;
        for(int dof_offset=0; dof_offset<fieldSize; ++dof_offset) {
          dof_ids[doffset++] = dof_id + dof_offset;
        }
      }
    }
  }

  if (stiff != NULL) {
    CHK_ERR( feData_->setElemMatrix(blockNumber, elemNumber, numNodes,
                                    nodeNumbers, dofsPerNode, dof_ids, stiff) );
  }

  if (load != NULL) {
    CHK_ERR( feData_->setElemVector(blockNumber, elemNumber, numNodes,
                                    nodeNumbers, dofsPerNode, dof_ids, load) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::putIntoRHS(int IDType,
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
int FEDataFilter::sumIntoRHS(int IDType,
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
int FEDataFilter::sumIntoMatrixDiagonal(int  IDType,
                             int  fieldID,
                             int  numIDs,
                             const GlobalID*  IDs,
                             const double*  coefficients)
{
  int fieldSize = problemStructure_->getFieldSize(fieldID);
  const NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  std::vector<int> eqns;
  convert_field_and_nodes_to_eqns(nodeDB, fieldID, fieldSize, numIDs, IDs, eqns);

  std::vector<int> nodeNumbers, dof_ids;
  convert_eqns_to_nodenumbers_and_dof_ids(problemStructure_->getFieldDofMap(),
                                        nodeDB, eqns.size(), &eqns[0],
                                      nodeNumbers, dof_ids);

  std::vector<int> ones(nodeNumbers.size(), 1);

  CHK_ERR( feData_->sumIntoMatrix(nodeNumbers.size(), &nodeNumbers[0], &dof_ids[0],
                                  &ones[0], &nodeNumbers[0], &dof_ids[0], coefficients));
  return( 0 );
}

//------------------------------------------------------------------------------
int FEDataFilter::enforceEssentialBCs(const int* eqns, 
                                      const double* alpha,
                                      const double* gamma, 
                                      int numEqns)
{
  std::vector<double> values;
  std::vector<int> nodeNumbers;
  std::vector<int> dof_ids;
  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  for(int i=0; i<numEqns; i++) {
    int reducedEqn = -1;
    bool isSlave = problemStructure_->
      translateToReducedEqn(eqns[i], reducedEqn);
    if (isSlave) continue;

    int nodeNumber = problemStructure_->getAssociatedNodeNumber(eqns[i]);

    nodeNumbers.push_back(nodeNumber);

    const NodeDescriptor* node = NULL;
    CHK_ERR( problemStructure_->getNodeDatabase().
             getNodeWithNumber(nodeNumber, node));

    int fieldID, offset;
    node->getFieldID(eqns[i], fieldID, offset);
    dof_ids.push_back( fdmap.get_dof_id(fieldID, offset) );

    values.push_back(gamma[i]/alpha[i]);
  }

  CHK_ERR( feData_->setDirichletBCs(nodeNumbers.size(),
                                    &nodeNumbers[0], &dof_ids[0], &values[0]) );

  newData_ = true;

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadFEDataMultCR(int CRID,
                               int numCRNodes,
                               const GlobalID* CRNodes, 
                               const int* CRFields,
                               const double* CRWeights,
                               double CRValue)
{
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

  if (numCRNodes <= 0) return(0);

  std::vector<int> nodeNumbers;
  std::vector<int> dof_ids;
  std::vector<double> weights;

  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  double fei_min = std::numeric_limits<double>::min();

  int offset = 0;
  for(int i=0; i<numCRNodes; i++) {
    NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithID(CRNodes[i], node) );

    int fieldEqn = -1;
    bool hasField = node->getFieldEqnNumber(CRFields[i], fieldEqn);
    if (!hasField) ERReturn(-1);

    int fieldSize = problemStructure_->getFieldSize(CRFields[i]);
    int dof_id = fdmap.get_dof_id(CRFields[i], 0);

    for(int f=0; f<fieldSize; f++) {
      double weight = CRWeights[offset++];
      if (std::abs(weight) > fei_min) {
        nodeNumbers.push_back(node->getNodeNumber());
        dof_ids.push_back(dof_id+f);
        weights.push_back(weight);
      }
    }
  }

  CHK_ERR( feData_->setMultiplierCR(CRID, nodeNumbers.size(),
                                    &nodeNumbers[0], &dof_ids[0],
                                    &weights[0], CRValue) );
  newData_ = true;

  return(0);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadFEDataPenCR(int CRID,
                              int numCRNodes,
                              const GlobalID* CRNodes, 
                              const int* CRFields,
                              const double* CRWeights,
                              double CRValue, 
                              double penValue)
{
  if (numCRNodes <= 0) return(0);

  std::vector<int> nodeNumbers;
  std::vector<int> dofsPerNode;
  std::vector<int> dof_ids;
  std::vector<double> weights;

  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  int offset = 0;
  for(int i=0; i<numCRNodes; i++) {
    NodeDescriptor* node = NULL; 
    nodeDB.getNodeWithID(CRNodes[i], node);
    if(node == NULL) continue;

    int fieldEqn = -1;
    bool hasField = node->getFieldEqnNumber(CRFields[i], fieldEqn);
    // If a node doesn't have a field, skip it.
    if (!hasField) continue;

    int fieldSize = problemStructure_->getFieldSize(CRFields[i]);

    nodeNumbers.push_back(node->getNodeNumber());
    dofsPerNode.push_back(fieldSize);

    for(int f=0; f<fieldSize; f++) {
      dof_ids.push_back(fdmap.get_dof_id(CRFields[i], f));
      double weight = CRWeights[offset++];
      weights.push_back(weight);
    }
  }

  std::vector<double*> matrixCoefs(weights.size());
  std::vector<double> rhsCoefs(weights.size());
  offset = 0;
  for(size_t i=0; i<weights.size(); ++i) {
    double* coefPtr = new double[weights.size()];
    for(size_t j=0; j<weights.size(); ++j) {
      coefPtr[j] = weights[i]*weights[j]*penValue;
    }
    matrixCoefs[i] = coefPtr;
    rhsCoefs[i] = weights[i]*penValue*CRValue;
  }

  int crIndex = fei::binarySearch(CRID, penCRIDs_);

  int index = fei::binarySearch(numCRNodes, constraintBlocks_);

  int blockNum = problemStructure_->getNumElemBlocks() + index;
  int elemNum = numRegularElems_ + crIndex;

  CHK_ERR( feData_->setElemMatrix(blockNum, elemNum,
                                  nodeNumbers.size(),
                                  &nodeNumbers[0],
                                  &dofsPerNode[0],
                                  &dof_ids[0],
                                  &matrixCoefs[0]) );

  CHK_ERR( feData_->setElemVector(blockNum, elemNum, nodeNumbers.size(),
                                  &nodeNumbers[0], &dofsPerNode[0], &dof_ids[0], &rhsCoefs[0]) );

  newData_ = true;

  for(size_t i=0; i<weights.size(); ++i) {
    delete [] matrixCoefs[i];
  }

  return(0);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadCRMult(int CRID, 
                         int numCRNodes,
                         const GlobalID* CRNodes, 
                         const int* CRFields,
                         const double* CRWeights,
                         double CRValue)
{
//
// Load Lagrange multiplier constraint relation data
//

  //Give the constraint data to the underlying solver using this special function...
  CHK_ERR( loadFEDataMultCR(CRID, numCRNodes, CRNodes, CRFields, CRWeights,
                            CRValue) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadCRPen(int CRID, 
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

   //Give the constraint data to the underlying solver using this special function...
   CHK_ERR( loadFEDataPenCR(CRID, numCRNodes, CRNodes, CRFields, CRWeights,
                            CRValue, penValue) );

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::parameters(int numParams, const char *const* paramStrings)
{
//
// this function takes parameters for setting internal things like solver
// and preconditioner choice, etc.
//
   if (numParams == 0 || paramStrings == NULL) {
      debugOutput("#FEDataFilter::parameters --- no parameters.");
   }
   else {

      snl_fei::getIntParamValue("outputLevel",numParams, paramStrings,
                                outputLevel_);

      snl_fei::getIntParamValue("internalFei",numParams,paramStrings,
                                internalFei_);

      if (Filter::logStream() != NULL) {
        (*logStream())<<"#FEDataFilter::parameters"<<FEI_ENDL
                         <<"# --- numParams: "<< numParams<<FEI_ENDL;
         for(int i=0; i<numParams; i++){
           (*logStream())<<"#------ paramStrings["<<i<<"]: "
                            <<paramStrings[i]<<FEI_ENDL;
         }
      }
   }

   CHK_ERR( Filter::parameters(numParams, paramStrings) );

   debugOutput("#FEDataFilter leaving parameters function");

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::loadComplete()
{
  debugOutput("#FEDataFilter calling FEData matrixLoadComplete");

  CHK_ERR( feData_->loadComplete() );

  newData_ = false;

  return(0);
}

//------------------------------------------------------------------------------
int FEDataFilter::residualNorm(int whichNorm, int numFields,
                           int* fieldIDs, double* norms, double& residTime)
{
//
//This function can do 3 kinds of norms: infinity-norm (denoted
//by whichNorm==0), 1-norm and 2-norm.
//
   debugOutput("FEI: residualNorm");

   CHK_ERR( loadComplete() );

   //for now, FiniteElementData doesn't do residual calculations.

   int fdbNumFields = problemStructure_->getNumFields();
   const int* fdbFieldIDs = problemStructure_->getFieldIDsPtr();

   int i;

   //Since we don't calculate actual residual norms, we'll fill the user's
   //array with norm data that is obviously not real norm data.
   int offset = 0;
   i = 0;
   while(offset < numFields && i < fdbNumFields) {
     if (fdbFieldIDs[i] >= 0) {
       fieldIDs[offset++] = fdbFieldIDs[i];
     }
     ++i;
   }
   for(i=0; i<numFields; ++i) {
      norms[i] = -99.9;
   }

   //fill out the end of the array with garbage in case the user-provided
   //array is longer than the list of fields we have in fieldDB.
   for(i=offset; i<numFields; ++i) {
      fieldIDs[i] = -99;
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::formResidual(double* residValues, int numLocalEqns)
{
  //FiniteElementData implementations can't currently do residuals.
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::solve(int& status, double& sTime) {

   debugOutput("FEI: solve");

   CHK_ERR( loadComplete() );

   debugOutput("#FEDataFilter in solve, calling launchSolver...");
 
   double start = MPI_Wtime();

   CHK_ERR( feData_->launchSolver(status, iterations_) );

   sTime = MPI_Wtime() - start;

   debugOutput("#FEDataFilter... back from solver");
 
   //now unpack the locally-owned shared entries of the solution vector into
   //the eqn-comm-mgr data structures.
   CHK_ERR( unpackSolution() );

   debugOutput("#FEDataFilter leaving solve");

   if (status != 0) return(1);
   else return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::setNumRHSVectors(int numRHSs, int* rhsIDs){

   if (numRHSs < 0) {
      fei::console_out() << "FEDataFilter::setNumRHSVectors: ERROR, numRHSs < 0." << FEI_ENDL;
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
int FEDataFilter::setCurrentRHS(int rhsID)
{
   std::vector<int>::iterator iter =
       std::find( rhsIDs_.begin(), rhsIDs_.end(), rhsID);

   if (iter == rhsIDs_.end()) ERReturn(-1)
   
   currentRHS_ = iter - rhsIDs_.begin();

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::giveToMatrix(int numPtRows, const int* ptRows,
                           int numPtCols, const int* ptCols,
                           const double* const* values, int mode)
{
  //This isn't going to be fast... I need to optimize the whole structure
  //of code that's associated with passing data to FiniteElementData.

  std::vector<int> rowNodeNumbers, row_dof_ids, colNodeNumbers, col_dof_ids;
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  int i;

  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  //First, we have to get nodeNumbers and dof_ids for each of the
  //row-numbers and col-numbers.

  for(i=0; i<numPtRows; i++) {
    int nodeNumber = problemStructure_->getAssociatedNodeNumber(ptRows[i]);
    if (nodeNumber < 0) ERReturn(-1);
    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithNumber(nodeNumber, node) );
    int fieldID, offset;
    node->getFieldID(ptRows[i], fieldID, offset);

    rowNodeNumbers.push_back(nodeNumber);
    row_dof_ids.push_back(fdmap.get_dof_id(fieldID, offset));
  }

  for(i=0; i<numPtCols; i++) {
    int nodeNumber = problemStructure_->getAssociatedNodeNumber(ptCols[i]);
    if (nodeNumber < 0) ERReturn(-1);
    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithNumber(nodeNumber, node) );
    int fieldID, offset;
    node->getFieldID(ptCols[i], fieldID, offset);

    colNodeNumbers.push_back(nodeNumber);
    col_dof_ids.push_back(fdmap.get_dof_id(fieldID, offset));
  }

  //now we have to flatten the colNodeNumbers and col_dof_ids out into
  //an array of length numPtRows*numPtCols, where the nodeNumbers and
  //dof_ids are repeated 'numPtRows' times.

  int len = numPtRows*numPtCols;
  std::vector<int> allColNodeNumbers(len), all_col_dof_ids(len);
  std::vector<double> allValues(len);

  int offset = 0;
  for(i=0; i<numPtRows; i++) {
    for(int j=0; j<numPtCols; j++) {
      allColNodeNumbers[offset] = colNodeNumbers[j];
      all_col_dof_ids[offset] = col_dof_ids[j];
      allValues[offset++] = values[i][j];
    }
  }

  //while we're at it, let's make an array with numPtCols replicated in it
  //'numPtRows' times.
  std::vector<int> numColsPerRow(numPtRows, numPtCols);

  //now we're ready to hand this stuff off to the FiniteElementData
  //instantiation.

  CHK_ERR( feData_->sumIntoMatrix(numPtRows,
                                  &rowNodeNumbers[0],
                                  &row_dof_ids[0],
                                  &numColsPerRow[0],
                                  &allColNodeNumbers[0],
                                  &all_col_dof_ids[0],
                                  &allValues[0]) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::giveToLocalReducedMatrix(int numPtRows, const int* ptRows,
                                       int numPtCols, const int* ptCols,
                                       const double* const* values, int mode)
{
  //This isn't going to be fast... I need to optimize the whole structure
  //of code that's associated with passing data to FiniteElementData.

  std::vector<int> rowNodeNumbers, row_dof_ids, colNodeNumbers, col_dof_ids;
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  int i;

  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  //First, we have to get nodeNumbers and dof_ids for each of the
  //row-numbers and col-numbers.

  for(i=0; i<numPtRows; i++) {
    int nodeNumber = problemStructure_->getAssociatedNodeNumber(ptRows[i]);
    if (nodeNumber < 0) ERReturn(-1);
    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithNumber(nodeNumber, node) );
    int fieldID, offset;
    node->getFieldID(ptRows[i], fieldID, offset);

    rowNodeNumbers.push_back(nodeNumber);
    row_dof_ids.push_back(fdmap.get_dof_id(fieldID, offset));
  }

  for(i=0; i<numPtCols; i++) {
    int nodeNumber = problemStructure_->getAssociatedNodeNumber(ptCols[i]);
    if (nodeNumber < 0) ERReturn(-1);
    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithNumber(nodeNumber, node) );
    int fieldID, offset;
    node->getFieldID(ptCols[i], fieldID, offset);

    colNodeNumbers.push_back(nodeNumber);
    col_dof_ids.push_back(fdmap.get_dof_id(fieldID, offset));
  }

  //now we have to flatten the colNodeNumbers and col_dof_ids out into
  //an array of length numPtRows*numPtCols, where the nodeNumbers and
  //dof_ids are repeated 'numPtRows' times.

  int len = numPtRows*numPtCols;
  std::vector<int> allColNodeNumbers(len), all_col_dof_ids(len);
  std::vector<double> allValues(len);

  int offset = 0;
  for(i=0; i<numPtRows; i++) {
    for(int j=0; j<numPtCols; j++) {
      allColNodeNumbers[offset] = colNodeNumbers[j];
      all_col_dof_ids[offset] = col_dof_ids[j];
      allValues[offset++] = values[i][j];
    }
  }

  //while we're at it, let's make an array with numPtCols replicated in it
  //'numPtRows' times.
  std::vector<int> numColsPerRow(numPtRows, numPtCols);

  //now we're ready to hand this stuff off to the FiniteElementData
  //instantiation.

  CHK_ERR( feData_->sumIntoMatrix(numPtRows,
                                  &rowNodeNumbers[0],
                                  &row_dof_ids[0],
                                  &numColsPerRow[0],
                                  &allColNodeNumbers[0],
                                  &all_col_dof_ids[0],
                                  &allValues[0]) );

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::getFromMatrix(int numPtRows, const int* ptRows,
                            const int* rowColOffsets, const int* ptCols,
                            int numColsPerRow, double** values)
{
  return(-1);

}

//------------------------------------------------------------------------------
int FEDataFilter::getEqnsFromMatrix(ProcEqns& procEqns, EqnBuffer& eqnData)
{
  ERReturn(-1);
}

//------------------------------------------------------------------------------
int FEDataFilter::getEqnsFromRHS(ProcEqns& procEqns, EqnBuffer& eqnData)
{
  ERReturn(-1);
}

//------------------------------------------------------------------------------
int FEDataFilter::giveToRHS(int num, const double* values,
                        const int* indices, int mode)
{
  std::vector<int> workspace(num*2);
  int* rowNodeNumbers = &workspace[0];
  int* row_dof_ids  = rowNodeNumbers+num;
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  for(int i=0; i<num; ++i) {
    const NodeDescriptor* nodeptr = 0;
    int err = nodeDB.getNodeWithEqn(indices[i], nodeptr);
    if (err < 0) { 
        rowNodeNumbers[i] = -1;
        row_dof_ids[i] = -1;
        continue;
    }

    rowNodeNumbers[i] = nodeptr->getNodeNumber();

    int fieldID, offset;
    nodeptr->getFieldID(indices[i], fieldID, offset);

    row_dof_ids[i] = fdmap.get_dof_id(fieldID, offset);
  }

  if (mode == ASSEMBLE_SUM) {
    CHK_ERR( feData_->sumIntoRHSVector(num,
                                       rowNodeNumbers,
                                       row_dof_ids,
                                       values) );
  }
  else {
    CHK_ERR( feData_->putIntoRHSVector(num,
                                       rowNodeNumbers,
                                       row_dof_ids,
                                       values) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::giveToLocalReducedRHS(int num, const double* values,
                                    const int* indices, int mode)
{
  std::vector<int> rowNodeNumbers, row_dof_ids;
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  fei::FieldDofMap<int>& fdmap = problemStructure_->getFieldDofMap();

  for(int i=0; i<num; i++) {
    int nodeNumber = problemStructure_->getAssociatedNodeNumber(indices[i]);
    if (nodeNumber < 0) ERReturn(-1);

    const NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithNumber(nodeNumber, node) );

    int fieldID, offset;
    node->getFieldID(indices[i], fieldID, offset);

    rowNodeNumbers.push_back(nodeNumber);
    row_dof_ids.push_back(fdmap.get_dof_id(fieldID, offset));
  }

  if (mode == ASSEMBLE_SUM) {
    CHK_ERR( feData_->sumIntoRHSVector(rowNodeNumbers.size(),
                                       &rowNodeNumbers[0],
                                       &row_dof_ids[0], values) );
  }
  else {
    CHK_ERR( feData_->putIntoRHSVector(rowNodeNumbers.size(),
                                       &rowNodeNumbers[0],
                                       &row_dof_ids[0], values) );
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::getFromRHS(int num, double* values, const int* indices)
{
  return(-1);
}

//------------------------------------------------------------------------------
int FEDataFilter::getEqnSolnEntry(int eqnNumber, double& solnValue)
{
  //This function's task is to retrieve the solution-value for a global
  //equation-number. eqnNumber may or may not be a slave-equation, and may or
  //may not be locally owned. If it is not locally owned, it should at least
  //be shared.
  //return 0 if the solution is successfully retrieved, otherwise return 1.
  //

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

//------------------------------------------------------------------------------
int FEDataFilter::getSharedRemoteSolnEntry(int eqnNumber, double& solnValue)
{
  std::vector<int>& remoteEqnNumbers = eqnCommMgr_->sendEqnNumbersPtr();
  double* remoteSoln = eqnCommMgr_->sendEqnSolnPtr();

  int index = fei::binarySearch(eqnNumber, remoteEqnNumbers);
  if (index < 0) {
    fei::console_out() << "FEDataFilter::getSharedRemoteSolnEntry: ERROR, eqn "
         << eqnNumber << " not found." << FEI_ENDL;
    ERReturn(-1);
  }
  solnValue = remoteSoln[index];
  return(0);
}

//------------------------------------------------------------------------------
int FEDataFilter::getReducedSolnEntry(int eqnNumber, double& solnValue)
{
  //We may safely assume that this function is called with 'eqnNumber' that is
  //local in the underlying assembled linear system. i.e., it isn't a slave-
  //equation, it isn't remotely owned, etc.
  //

  int nodeNumber = problemStructure_->getAssociatedNodeNumber(eqnNumber);

  //if nodeNumber < 0, it probably means we're trying to look up the
  //node for a lagrange-multiplier (which doesn't exist). In that
  //case, we're just going to ignore the request and return for now...
  if (nodeNumber < 0) {solnValue = -999.99; return(FEI_SUCCESS);}

  const NodeDescriptor* node = NULL;
  problemStructure_->getNodeDatabase().getNodeWithNumber(nodeNumber, node);
  if(node == NULL) {
    // KHP: If a node doesn't exist, we still need to
    // return a solution value....Zero seems like a logical
    // choice however, FEI_SUCCESS seems wrong however I don't
    // want to trip any asserts or other error conditions.
    solnValue = 0.0;
    return FEI_SUCCESS;
  }

  int eqn = problemStructure_->translateFromReducedEqn(eqnNumber);
  int fieldID, offset;
  node->getFieldID(eqn, fieldID, offset);
  int dof_id = problemStructure_->getFieldDofMap().get_dof_id(fieldID, offset);

  bool fetiHasNode = true;
  GlobalID nodeID = node->getGlobalNodeID();
  NodeCommMgr& nodeCommMgr = problemStructure_->getNodeCommMgr();
  std::vector<GlobalID>& shNodeIDs = nodeCommMgr.getSharedNodeIDs();
  int shIndex = fei::binarySearch(nodeID, shNodeIDs);
  if (shIndex >= 0) {
    if (!(problemStructure_->isInLocalElement(nodeNumber)) ) fetiHasNode = false;
  }

  if (fetiHasNode) {
    int err = feData_->getSolnEntry(nodeNumber, dof_id, solnValue);
    if (err != 0) {
      fei::console_out() << "FEDataFilter::getReducedSolnEntry: nodeNumber " << nodeNumber
           << " (nodeID " << node->getGlobalNodeID() << "), dof_id "<<dof_id
           << " couldn't be obtained from FETI on proc " << localRank_ << FEI_ENDL;
      ERReturn(-1);
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::unpackSolution()
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
       fei::console_out() << "FEDataFilter::unpackSolution: ERROR, 'recv' eqn (" << eqn
             << ") out of local range." << FEI_ENDL;
       MPI_Abort(comm_, -1);
     }

     double solnValue = 0.0;

     CHK_ERR( getReducedSolnEntry(eqn, solnValue) );

     eqnCommMgr_->addSolnValues(&eqn, &solnValue, 1);
   }

   eqnCommMgr_->exchangeSoln();

  debugOutput("#FEDataFilter leaving unpackSolution");
  return(FEI_SUCCESS);
}
             
//------------------------------------------------------------------------------
void FEDataFilter::  setEqnCommMgr(EqnCommMgr* eqnCommMgr)
{
  delete eqnCommMgr_;
  eqnCommMgr_ = eqnCommMgr;
}

//------------------------------------------------------------------------------
int FEDataFilter::getBlockNodeSolution(GlobalID elemBlockID,  
                                       int numNodes, 
                                       const GlobalID *nodeIDs, 
                                       int *offsets,
                                       double *results)
{        
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
     NodeDescriptor* node_i = NULL;
     nodeDB.getNodeAtIndex(i, node_i);

      if (offset == numNodes) break;

      GlobalID nodeID = nodeIDs[offset];

      //first let's set the offset at which this node's solution coefs start.
      offsets[offset++] = numSolnParams;

      NodeDescriptor* node = NULL;
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
          if (size < 1) {
            continue;
          }

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
int FEDataFilter::getNodalSolution(int numNodes, 
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
    NodeDescriptor* node_i = NULL;
    nodeDB.getNodeAtIndex(i, node_i);

    if (offset == numNodes) break;

    GlobalID nodeID = nodeIDs[offset];

    //first let's set the offset at which this node's solution coefs start.
    offsets[offset++] = numSolnParams;

    NodeDescriptor* node = NULL;
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
      if (size < 1) {
        continue;
      }

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
int FEDataFilter::getBlockFieldNodeSolution(GlobalID elemBlockID,
                                        int fieldID,
                                        int numNodes, 
                                        const GlobalID *nodeIDs, 
                                        double *results)
{
  debugOutput("FEI: getBlockFieldNodeSolution");

  int numActiveNodes = problemStructure_->getNumActiveNodes();
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  if (numActiveNodes <= 0) return(0);

  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  if (fieldSize <= 0) ERReturn(-1);

  if (!block->containsField(fieldID)) {
    fei::console_out() << "FEDataFilter::getBlockFieldNodeSolution WARNING: fieldID " << fieldID
         << " not contained in element-block " << static_cast<int>(elemBlockID) << FEI_ENDL;
    return(1);
  }

   //Traverse the node list, checking if nodes are associated with this block.
   //If so, put the answers in the results list.

   for(int i=0; i<numNodes; i++) {
     NodeDescriptor* node_i = NULL;
     nodeDB.getNodeAtIndex(i, node_i);

     GlobalID nodeID = nodeIDs[i];

     NodeDescriptor* node = NULL;
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
int FEDataFilter::getNodalFieldSolution(int fieldID,
                                        int numNodes, 
                                        const GlobalID *nodeIDs, 
                                        double *results)
{
  debugOutput("FEI: getNodalFieldSolution");

  int numActiveNodes = problemStructure_->getNumActiveNodes();
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  if (numActiveNodes <= 0) return(0);

  if (problemStructure_->numSlaveEquations() != 0) {
    fei::console_out() << "FEDataFilter::getEqnSolnEntry ERROR FETI-support is not currently"
         << " compatible with the FEI's constraint reduction." << FEI_ENDL;
    ERReturn(-1);
  }

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  if (fieldSize <= 0) {
    ERReturn(-1);
  }

  NodeCommMgr& nodeCommMgr = problemStructure_->getNodeCommMgr();

  //Traverse the node list, checking if nodes have the specified field.
  //If so, put the answers in the results list.

  for(int i=0; i<numNodes; i++) {
    NodeDescriptor* node_i = NULL;
    nodeDB.getNodeAtIndex(i, node_i);

    GlobalID nodeID = nodeIDs[i];

    NodeDescriptor* node = NULL;
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

    int nodeNumber = node->getNodeNumber();

    int eqnNumber = -1;
    bool hasField = node->getFieldEqnNumber(fieldID, eqnNumber);

    //If this node doesn't have the specified field, then skip to the
    //next loop iteration.
    if (!hasField) continue;

    std::vector<GlobalID>& shNodeIDs = nodeCommMgr.getSharedNodeIDs();
    int shIndex = fei::binarySearch(nodeID, &shNodeIDs[0], shNodeIDs.size());
    if (shIndex > -1) {
      if (!(problemStructure_->isInLocalElement(nodeNumber))) continue;
    }

    int dof_id = problemStructure_->getFieldDofMap().get_dof_id(fieldID, 0);

    int offset = fieldSize*i;

    for(int j=0; j<fieldSize; j++) {
      if (localStartRow_ > eqnNumber || eqnNumber > localEndRow_) {
        CHK_ERR( getSharedRemoteSolnEntry(eqnNumber+j, results[offset+j]) );
        continue;
      }

      err = feData_->getSolnEntry(nodeNumber, dof_id+j, results[offset+j]);
      if (err != 0) {
        fei::console_out() << "FEDataFilter::getReducedSolnEntry: nodeNumber " << nodeNumber
             << " (nodeID " << nodeID << "), dof_id "<<dof_id
             << " couldn't be obtained from FETI on proc " << localRank_ << FEI_ENDL;
        ERReturn(-1);
      }
    }
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::putBlockNodeSolution(GlobalID elemBlockID,
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

   int blk_idx = problemStructure_->getIndexOfBlock(elemBlockID);

   for(int i=0; i<numNodes; i++) {
     NodeDescriptor* node = NULL;
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
            }
         }
         offs += size;
      }//for(j<numFields)loop
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                        int fieldID, 
                                        int numNodes, 
                                        const GlobalID *nodeIDs, 
                                        const double *estimates)
{
   debugOutput("FEI: putBlockFieldNodeSolution");

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) );
   if (!block->containsField(fieldID)) return(1);

   int fieldSize = problemStructure_->getFieldSize(fieldID);
   NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

   //if we have a negative fieldID, we'll need a list of length numNodes,
   //in which to put nodeNumbers for passing to the solver... 

   std::vector<int> numbers(numNodes);

   //if we have a fieldID >= 0, then our numbers list will hold equation numbers
   //and we'll need fieldSize*numNodes of them.

   std::vector<double> data;

   if (fieldID >= 0) {
     if (fieldSize < 1) {
       fei::console_out() << "FEI Warning, putBlockFieldNodeSolution called for field "
            << fieldID<<", which has size "<<fieldSize<<FEI_ENDL;
       return(0);
     }
     try {
     numbers.resize(numNodes*fieldSize);
     data.resize(numNodes*fieldSize);
     }
     catch(std::runtime_error& exc) {
       fei::console_out() << exc.what()<<FEI_ENDL;
       ERReturn(-1);
     }
   }

   int count = 0;

   for(int i=0; i<numNodes; i++) {
     NodeDescriptor* node = NULL;
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
     CHK_ERR( feData_->putNodalFieldData(fieldID, fieldSize, 
                                         numNodes, &numbers[0],
                                         estimates));
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::getBlockElemSolution(GlobalID elemBlockID,
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

   std::map<GlobalID,int>& elemIDList = problemStructure_->
                          getBlockConnectivity(elemBlockID).elemIDs;

   int len = block->getNumElements();

   //if the user is only asking for a subset of element-solutions, shrink len.
   if (len > numElems) len = numElems;

   numElemDOFPerElement = block->getNumElemDOFPerElement();
   std::vector<int>& elemDOFEqnNumbers = block->elemDOFEqnNumbers();
   double answer;


   if (numElemDOFPerElement <= 0) return(0);

   std::map<GlobalID,int>::const_iterator
     elemid_end = elemIDList.end(),
     elemid_itr = elemIDList.begin();

   for(int i=0; i<len; i++) {
      int index = i;

      //if the user-supplied elemIDs are out of order, we need the index of
      //the location of this element.
      if (elemid_itr->first != elemIDs[i]) {
         index = elemid_itr->second;
      }

      if (index < 0) continue;

      int offset = i*numElemDOFPerElement;

      for(int j=0; j<numElemDOFPerElement; j++) {
         int eqn = elemDOFEqnNumbers[index] + j;

         CHK_ERR( getEqnSolnEntry(eqn, answer) )

         results[offset+j] = answer;
      }

      ++elemid_itr;
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::putBlockElemSolution(GlobalID elemBlockID,
                                   int numElems,
                                   const GlobalID *elemIDs,
                                   int dofPerElem,
                                   const double *estimates)
{
   debugOutput("FEI: putBlockElemSolution");

   BlockDescriptor* block = NULL;
   CHK_ERR( problemStructure_->getBlockDescriptor(elemBlockID, block) )

   std::map<GlobalID,int>& elemIDList = problemStructure_->
                          getBlockConnectivity(elemBlockID).elemIDs;

   int len = block->getNumElements();
   if (len > numElems) len = numElems;

   int DOFPerElement = block->getNumElemDOFPerElement();
   if (DOFPerElement != dofPerElem) {
     fei::console_out() << "FEI ERROR, putBlockElemSolution called with bad 'dofPerElem' ("
          <<dofPerElem<<"), block "<<elemBlockID<<" should have dofPerElem=="
          <<DOFPerElement<<FEI_ENDL;
     ERReturn(-1);
   }

   std::vector<int>& elemDOFEqnNumbers = block->elemDOFEqnNumbers();

   if (DOFPerElement <= 0) return(0);

   std::map<GlobalID,int>::const_iterator
     elemid_end = elemIDList.end(),
     elemid_itr = elemIDList.begin();

   for(int i=0; i<len; i++) {
      int index = i;
      if (elemid_itr->first != elemIDs[i]) {
         index = elemid_itr->second;
      }

      if (index < 0) continue;

      for(int j=0; j<DOFPerElement; j++) {
         int reducedEqn;
         problemStructure_->
           translateToReducedEqn(elemDOFEqnNumbers[i] + j, reducedEqn);
//         double soln = estimates[i*DOFPerElement + j];

//       if (useLinSysCore_) {
//            CHK_ERR( lsc_->putInitialGuess(&reducedEqn, &soln, 1) );
//          }
      }

      ++elemid_itr;
   }

   return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::getCRMultipliers(int numCRs,
                                   const int* CRIDs,
                                   double* multipliers)
{
  for(int i=0; i<numCRs; i++) {
    //temporarily, FETI's getMultiplierSoln method isn't implemented.
    //CHK_ERR( feData_->getMultiplierSoln(CRIDs[i], multipliers[i]) );
    multipliers[i] = -999.99;
  }

  return(-1);
}

//------------------------------------------------------------------------------
int FEDataFilter::putCRMultipliers(int numMultCRs,
                               const int* CRIDs,
                               const double *multEstimates)
{
  debugOutput("FEI: putCRMultipliers");

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::putNodalFieldData(int fieldID,
                                int numNodes,
                                const GlobalID* nodeIDs,
                                const double* nodeData)
{
  debugOutput("FEI: putNodalFieldData");

  int fieldSize = problemStructure_->getFieldSize(fieldID);
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  std::vector<int> nodeNumbers(numNodes);

  for(int i=0; i<numNodes; i++) {
    NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithID(nodeIDs[i], node) );

    int nodeNumber = node->getNodeNumber();
    if (nodeNumber < 0) {
      GlobalID nodeID = nodeIDs[i];
      fei::console_out() << "FEDataFilter::putNodalFieldData ERROR, node with ID " 
           << static_cast<int>(nodeID) << " doesn't have an associated nodeNumber "
           << "assigned. putNodalFieldData shouldn't be called until after the "
           << "initComplete method has been called." << FEI_ENDL;
      ERReturn(-1);
    }

    nodeNumbers[i] = nodeNumber;
  }

  CHK_ERR( feData_->putNodalFieldData(fieldID, fieldSize,
                                      numNodes, &nodeNumbers[0],
                                      nodeData) );

  return(0);
}

//------------------------------------------------------------------------------
int FEDataFilter::putNodalFieldSolution(int fieldID,
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
    NodeDescriptor* node = NULL;
    CHK_ERR( nodeDB.getNodeWithID(nodeIDs[i], node) );

    int eqn = -1;
    bool hasField = node->getFieldEqnNumber(fieldID, eqn);
    if (!hasField) continue;

  }

  return(0);
}

//------------------------------------------------------------------------------
int FEDataFilter::assembleEqns(int numRows, 
                               int numCols,
                               const int* rowNumbers,
                               const int* colIndices,
                               const double* const* coefs,
                               bool structurallySymmetric,
                               int mode)
{
  if (numRows == 0) return(FEI_SUCCESS);

  const int* indPtr = colIndices;
  for(int i=0; i<numRows; i++) {
    int row = rowNumbers[i];

    const double* coefPtr = coefs[i];

    CHK_ERR(giveToMatrix(1, &row, numCols, indPtr, &coefPtr, mode));

    if (!structurallySymmetric) indPtr += numCols;
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEDataFilter::assembleRHS(int numValues,
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
    if (eqn < 0) continue;

    CHK_ERR( giveToRHS(1, &(coefs[i]), &eqn, mode ) );
  }

  return(FEI_SUCCESS);
}

//==============================================================================
void FEDataFilter::debugOutput(const char* mesg)
{
  if (Filter::logStream() != NULL) {
    (*logStream()) << mesg << FEI_ENDL;
   }
}

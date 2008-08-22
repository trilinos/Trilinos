/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_sstream.hpp>

#include <limits>
#include <cmath>

#include <fei_BCManager.hpp>
#include <fei_defs.h>
#include <fei_BCRecord.hpp>
#include <fei_MatrixGraph.hpp>

#include <feiPoolAllocator.hpp>

#include <fei_mpi.h>
#include <fei_NodeDatabase.hpp>
#include <fei_NodeDescriptor.hpp>
#include <fei_EqnBuffer.hpp>

#include <fei_TemplateUtils.hpp>

#include <fei_VectorSpace.hpp>
#include <fei_Matrix.hpp>

#undef fei_file
#define fei_file "fei_BCManager.cpp"
#include <fei_ErrMacros.hpp>

//==============================================================================
BCManager::BCManager()
 : bcList_(),
   bcAlloc_(NULL),
   coefAlloc_(NULL)
{
  bcList_.reserve(1000);
  bcAlloc_ = new feiPoolAllocator<BCRecord>(1000);
  coefAlloc_ = new feiPoolAllocator<double>(3000);
}

//==============================================================================
BCManager::~BCManager() {
  delete bcAlloc_;
  delete coefAlloc_;
}

//==============================================================================
void BCManager::addBCRecords(int numNodes, const GlobalID* nodeIDs,
			    int fieldID, int fieldSize,
                            const double*const * alpha,
			    const double*const * beta,
                            const double*const * gamma)
{
  size_t oldLength = bcList_.size();
  bcList_.resize(oldLength+numNodes);
  const BCRecord** bcListPtr = &bcList_[0];

  int len = fieldSize*3;
  for(int i=0; i<numNodes; i++) {
    BCRecord* bc = bcAlloc_->alloc();
    double* coefs = coefAlloc_->alloc(len);
    bc->init(nodeIDs[i], fieldID, fieldSize, coefs);
    //fill the new BCRecord with data.
    bc->setAlpha(alpha[i]);
    bc->setBeta(beta[i]);
    bc->setGamma(gamma[i]);
    bcListPtr[oldLength+i] = bc;
  }
}

//==============================================================================
void BCManager::addBCRecords(int idType, int numNodes, const GlobalID* nodeIDs,
			    int fieldID, int fieldSize,
                            const double*const * gamma,
			    const double*const * alpha)
{
  if (numNodes < 1) {
    return;
  }

  double* beta = new double[fieldSize];

  for(int ii=0; ii<fieldSize; ++ii) {
    beta[ii] = 0.0;
  }

  size_t oldLength = bcList_.size();
  bcList_.resize(oldLength+numNodes);
  const BCRecord** bcListPtr = &bcList_[0];

  int len = fieldSize*3;
  for(int i=0; i<numNodes; i++) {
    BCRecord* bc = bcAlloc_->alloc();
    double* coefs = coefAlloc_->alloc(len);
    bc->init(nodeIDs[i], fieldID, fieldSize, coefs);
    bc->setIDType(idType);
    //fill the new BCRecord with data.
    bc->setAlpha(alpha[i]);
    bc->setBeta(beta);
    bc->setGamma(gamma[i]);
    bcListPtr[oldLength+i] = bc;
  }

  delete [] beta;
}

//==============================================================================
void BCManager::addBCRecords(int idType, int numNodes, const GlobalID* nodeIDs,
                              int fieldID, int fieldSize,
                              const double*const * prescribedValues)
{
  if (numNodes < 1) {
    return;
  }

  double* alpha = new double[fieldSize*2];
  double* beta = alpha+fieldSize;

  for(int ii=0; ii<fieldSize; ++ii) {
    alpha[ii] = 1.0;
    beta[ii] = 0.0;
  }

  size_t oldLength = bcList_.size();
  bcList_.resize(oldLength+numNodes);
  const BCRecord** bcListPtr = &bcList_[0];

  int len = fieldSize*3;
  for(int i=0; i<numNodes; i++) {
    BCRecord* bc = bcAlloc_->alloc();
    double* coefs = coefAlloc_->alloc(len);
    bc->init(nodeIDs[i], fieldID, fieldSize, coefs);
    bc->setIDType(idType);
    //fill the new BCRecord with data.
    bc->setAlpha(alpha);
    bc->setBeta(beta);
    bc->setGamma(prescribedValues[i]);
    bcListPtr[oldLength+i] = bc;
  }

  delete [] alpha;
}

//==============================================================================
void BCManager::addBCRecords(int numNodes,
                             const GlobalID* nodeIDs,
                             int fieldID, int fieldSize,
                             const int* offsetsIntoField,
                             const double* prescribedValues)
{
  if (numNodes < 1) {
    return;
  }

  double* alpha = new double[fieldSize*3];
  double* beta = alpha+fieldSize;
  double* gamma = beta+fieldSize;

  for(int ii=0; ii<fieldSize; ++ii) {
    beta[ii] = 0.0;
  }

  size_t oldLength = bcList_.size();
  bcList_.resize(oldLength+numNodes);
  const BCRecord** bcListPtr = &bcList_[0];

  int len = fieldSize*3;
  for(int i=0; i<numNodes; i++) {
    for(int j=0; j<offsetsIntoField[i]; ++j) {
      alpha[j] = 0.0;
      gamma[j] = 0.0;
    }
    gamma[offsetsIntoField[i]] = prescribedValues[i];
    alpha[offsetsIntoField[i]] = 1.0;
    for(int j=offsetsIntoField[i]+1; j<fieldSize; ++j) {
      alpha[j] = 0.0;
      gamma[j] = 0.0;
    }

    BCRecord* bc = bcAlloc_->alloc();
    double* coefs = coefAlloc_->alloc(len);
    bc->init(nodeIDs[i], fieldID, fieldSize, coefs);
    //bc->setIDType(idType);
    //fill the new BCRecord with data.
    bc->setAlpha(alpha);
    bc->setBeta(beta);
    bc->setGamma(gamma);
    bcListPtr[oldLength+i] = bc;
  }

  delete [] alpha;
}

//==============================================================================
int BCManager::finalizeBCEqns(fei::Matrix& matrix,
                              bool throw_if_bc_slave_conflict)
{
  try {
    consolidateBCs();
  }
  catch(fei::Exception& exc) {
    FEI_CERR << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  size_t numBCs = getNumBCs();
  std::vector<const BCRecord*>& BCs = getBCRecords();

  fei::VectorSpace& vecSpace = *(matrix.getMatrixGraph()->getRowSpace());
  fei::SharedPtr<fei::Reducer> reducer = matrix.getMatrixGraph()->getReducer();
  bool haveSlaves = reducer.get()!=NULL;

  int cols[3];
  double coefs[3];
  const double* coefsPtr = coefs;
  //each "equation" will have 3 coefficients: alpha, beta, gamma. And they will
  //carry arbitrary column-indices.
  if (haveSlaves) {
    int ii=0, num=0;
    while(num<3) {
      if (!reducer->isSlaveEqn(ii++)) {
        cols[num++] = ii-1;
      }
    }
  }
  else {
    for(int ii=0; ii<3; ++ii) cols[ii] = ii;
  }

  for(unsigned i=0; i<numBCs; i++) {
    const BCRecord& bc = *(BCs[i]);

    int fieldID = bc.getFieldID();
    int eqn = -1;
    try {
      CHK_ERR( vecSpace.getGlobalIndex(bc.getIDType(), bc.getNodeID(), fieldID, eqn) );
    }
    catch(fei::Exception& exc) {
      FEI_OSTRINGSTREAM osstr;
      osstr << "BCManager::finalizeBCEqns caught exception: " <<exc.what()
	    << " BC node " << bc.getNodeID()<<", fieldID: " << fieldID;
      FEI_CERR << osstr.str()<<FEI_ENDL;
      ERReturn(-1);
    }

    if (haveSlaves) {
      if (reducer->isSlaveEqn(eqn)) {
        if (throw_if_bc_slave_conflict) {
          FEI_OSTRINGSTREAM osstr;
          osstr << "fei BCManager::finalizeBCeqns ERROR, eqn="<<eqn
            << " is both a BC eqn and slave-constraint eqn.";
          throw fei::Exception(osstr.str());
        }
        continue;
      }
    }

    int fieldSize = bc.getFieldSize();

    const double* alpha = bc.pointerToAlpha();
    const double* beta  = bc.pointerToBeta();
    const double* gamma = bc.pointerToGamma();

    for(int k = 0; k < fieldSize; k++) {
      int thisEqn = eqn + k;

      coefs[0] = alpha[k];
      coefs[1] = beta[k];
      coefs[2] = gamma[k];

      CHK_ERR( matrix.copyIn(1, &thisEqn, 3, cols, &coefsPtr) );
    }
  }

  clearAllBCs();

  return(0);
}

//==============================================================================
int BCManager::finalizeBCEqns(NodeDatabase& nodeDB,
			      EqnBuffer& bcEqns)
{
  try {
    consolidateBCs();
  }
  catch(fei::Exception& exc) {
    FEI_CERR << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  size_t numBCs = getNumBCs();
  std::vector<const BCRecord*>& BCs = getBCRecords();

  int cols[3];
  double coefs[3];
  //each "equation" will have 3 coefficients: alpha, beta, gamma. And they will
  //carry the arbitrary column-indices 0,1,2.
  for(int ii=0; ii<3; ii++) cols[ii] = ii;

  for(unsigned i=0; i<numBCs; i++) {
    const BCRecord& bc = *(BCs[i]);

    NodeDescriptor* node = NULL;
    nodeDB.getNodeWithID(bc.getNodeID(), node);
    int fieldID = bc.getFieldID();
    int fieldSize = bc.getFieldSize();

    int eqn = -1;
    if (!node->getFieldEqnNumber(fieldID, eqn)) {
      ERReturn(-1);
    }

    const double* alpha = bc.pointerToAlpha();
    const double* beta  = bc.pointerToBeta();
    const double* gamma = bc.pointerToGamma();

    for(int k = 0; k < fieldSize; k++) {
      int thisEqn = eqn + k;

      coefs[0] = alpha[k];
      coefs[1] = beta[k];
      coefs[2] = gamma[k];

      CHK_ERR( bcEqns.addEqn(thisEqn, coefs, cols, 3, false) );
    }
  }

  clearAllBCs();

  return(0);
}

//==============================================================================
int BCManager::consolidateBCs()
{
  //
  //This function runs all of the BCRecords that have been created, checking for
  //duplicates, which are BCRecords that have the same node and field. If a
  //duplicate is found, consolidation is done as follows:
  //
  // If the duplicate records are imposing essential (alpha != 0, beta == 0) BCs
  // on the node/field, then only the last one is kept (the BCRecords were
  // simply appended to the list as BCs came in from the application, so they
  // are in "order").
  //
  // If the duplicate records are imposing natural or mixed BCs on the
  // node/field, then the alpha/beta/gamma contributions are simply added
  // together in the BCRecord that is kept.
  //
  // If one record imposes an essential BC on a node/field and another record
  // imposes a natural or mixed BC on a node/field, then the essential BC is
  // kept and the natural/mixed is discarded.
  //
  size_t numBCs = getNumBCs();
  if (numBCs == 0) return(0);

  std::vector<const BCRecord*> consolidatedList;
  consolidatedList.reserve(numBCs);

  for(unsigned i=0; i<numBCs; i++) {
    BCRecord* lastBC = const_cast<BCRecord*>(bcList_[i]);

    const BCRecord** listdata = &consolidatedList[0];
    int insertPoint = -1;
    int index = snl_fei::binarySearchPtr<BCRecord>(lastBC, listdata,
				      (int)consolidatedList.size(), insertPoint);

    if (index >= 0) {
      //lastBC specifies the same nodeID/fieldID pair as a BCRecord that's
      //already in the consolidated list.
      BCRecord* bc = const_cast<BCRecord*>(consolidatedList[index]);

      double* bc_alpha = const_cast<double*>(bc->pointerToAlpha());
      double* bc_beta  = const_cast<double*>(bc->pointerToBeta());
      double* bc_gamma = const_cast<double*>(bc->pointerToGamma());

      const double* lastBC_alpha = lastBC->pointerToAlpha();
      const double* lastBC_beta  = lastBC->pointerToBeta();
      const double* lastBC_gamma = lastBC->pointerToGamma();

      int size = bc->getFieldSize();

      double fei_eps = 1.e-49;
      for(int k=0; k<size; k++) {
	if (std::abs(lastBC_alpha[k]) > fei_eps &&
            std::abs(lastBC_beta[k]) < fei_eps) {
	  //it's essential, so 'put' it in
	  bc_alpha[k] = lastBC_alpha[k];
	  bc_beta[k]  = lastBC_beta[k];
	  bc_gamma[k] = lastBC_gamma[k];
	}
	else {
	  //it's natural or mixed, so sum it in, but not if
	  //an essential BC was already in place.
	  if (!((std::abs(bc_alpha[k]) > fei_eps) &&
                (std::abs(bc_beta[k]) < fei_eps))){
	    bc_alpha[k] += lastBC_alpha[k];
	    bc_beta[k] += lastBC_beta[k];
	    bc_gamma[k] += lastBC_gamma[k];
	  }
	}
      }
    }
    else {
      //lastBC wasn't found, so insert it in the consolidated list.

      consolidatedList.insert(consolidatedList.begin()+insertPoint, lastBC);
    }
  }

  //now we can simply replace bcList_ with consolidatedList.
  bcList_ = consolidatedList;

  return(0);
}

//==============================================================================
size_t BCManager::getNumBCs()
{
  return( bcList_.size() );
}

//==============================================================================
void BCManager::clearAllBCs()
{
  bcList_.clear();
  delete bcAlloc_;
  delete coefAlloc_;
  bcAlloc_ = new feiPoolAllocator<BCRecord>(1000);
  coefAlloc_ = new feiPoolAllocator<double>(3000);
}


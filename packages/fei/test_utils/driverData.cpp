/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_fstream.hpp>

#include <FEI.hpp>
#include <fei_defs.h>

#include <snl_fei_Utils.hpp>
#include <test_utils/DataReader.hpp>

#include <test_utils/driverData.hpp>

#include <cstring>

#ifdef CHK_ERR
#undef CHK_ERR
#endif
#define CHK_ERR(a) {int chkerr; if ((chkerr = a) != 0) { \
                    fei::console_out() << "file " << __FILE__ << ", line " << __LINE__ \
                         << ", err " << chkerr << FEI_ENDL; \
                    return(chkerr); } }

driverData::driverData()
  :
  methodNames(),
  temp_(NULL),
  tempLen_(0),
  initFields_fieldSizes_(NULL),
  initFields_fieldIDs_(NULL),
  initElemBlock_ints_(NULL),
  initElemBlock_fieldsPerNode_(NULL),
  initElemBlock_fieldIDs_(NULL),
  initElemBlock_elemDofFieldIDs_(NULL),
  initElems_(),
  initElemCounter_(0),
  sumInElems_(),
  sumInElemCounter_(0),
  sumInElemMatrix_(),
  sumInElemMatrixCounter_(0),
  sumInElemRHS_(),
  sumInElemRHSCounter_(0),
  loadNodeBCs_(),
  loadNodeBCsCounter_(0),
  initCRMult_(),
  initCRMultCounter_(0),
  loadCRMult_(),
  loadCRMultCounter_(0),
  initSharedNodes_(),
  initSharedNodesCounter_(0),
  parameters_(),
  parametersCounter_(0),
  setIDLists_(),
  setIDListsCounter_(0),
  setCurrentMatrix_(),
  setCurrentMatrixCounter_(0),
  setCurrentRHS_(),
  setCurrentRHSCounter_(0),
  putBlockFieldNodeSolution_(),
  putBlockFieldNodeSolutionCounter_(0)
{
  tempLen_ = 512;
  temp_ = new char[tempLen_];
}

driverData::~driverData()
{
  tempLen_ = 0;
  delete [] temp_;

  for(size_t i=0; i<methodNames.size(); i++) {
    delete [] methodNames[i];
  }

  delete [] initFields_fieldSizes_;
  delete [] initFields_fieldIDs_;

  if (initElemBlock_numInts_ > 0 && initElemBlock_ints_ != NULL) {
    int len = initElemBlock_ints_[2];
    for(int j=0; j<len; ++j) delete [] initElemBlock_fieldIDs_[j];
    delete [] initElemBlock_fieldIDs_;
    delete [] initElemBlock_fieldsPerNode_;
    delete [] initElemBlock_ints_;
    initElemBlock_ints_ = NULL;
  }

  for(size_t i=0; i<initElems_.size(); ++i) delete initElems_[i];
  for(size_t i=0; i<sumInElems_.size(); ++i) delete sumInElems_[i];
  for(size_t i=0; i<sumInElemMatrix_.size(); ++i) delete sumInElemMatrix_[i];
  for(size_t i=0; i<sumInElemRHS_.size(); ++i) delete sumInElemRHS_[i];
  for(size_t i=0; i<loadNodeBCs_.size(); ++i) delete loadNodeBCs_[i];
  for(size_t i=0; i<initCRMult_.size(); ++i) delete initCRMult_[i];
  for(size_t i=0; i<loadCRMult_.size(); ++i) delete loadCRMult_[i];
  for(size_t i=0; i<initSharedNodes_.size(); ++i) delete initSharedNodes_[i];
  for(size_t i=0; i<parameters_.size(); ++i) delete parameters_[i];
  for(size_t i=0; i<setIDLists_.size(); ++i) delete setIDLists_[i];
}

int driverData::readData(const char* fileName)
{
  FEI_IFSTREAM* instr = NULL;
  instr = new FEI_IFSTREAM(fileName);

  if (instr->bad()) {
    fei::console_out() << "driverData::readData: ERROR opening " << fileName << FEI_ENDL;
    return(-1);
  }
  FEI_COUT << "driverData reading from " << fileName << FEI_ENDL;
  char* keyword = NULL;

  int err = getKeyword(instr, keyword);
  while (!instr->eof() && !err) {
    FEI_COUT << "driverData read keyword: " << keyword << FEI_ENDL;
    err = readData(instr, keyword);
    if (err != 0) {
      fei::console_out() << "driverData: ERROR reading data for keyword: " << keyword << FEI_ENDL;
      break;
    }
    delete [] keyword;
    err = getKeyword(instr, keyword);
  }

  delete instr;

  return(0);
}

int driverData::call_fei_method(const char* method, FEI* fei)
{
  if (!std::strcmp("setSolveType", method)) {
    return( fei->setSolveType(solveType_) );
  }

  if (!std::strcmp("setIDLists", method)) {
    if (setIDListsCounter_ >= (int)setIDLists_.size()) {
      fei::console_out() << "driverData ERROR, can't call setIDLists again" << FEI_ENDL;
      return(-1);
    }

    setIDLists* sidl = setIDLists_[setIDListsCounter_++];

    return( fei->setIDLists(sidl->numMatrices, sidl->matrixIDs,
			       sidl->numRHSs, sidl->rhsIDs) );
  }

  if (!std::strcmp("setCurrentMatrix", method)) {
    if (setCurrentMatrixCounter_ >= (int)setCurrentMatrix_.size()) {
      fei::console_out() << "driverData ERROR, can't call setCurrentMatrix again" << FEI_ENDL;
      return(-1);
    }

    int matID = setCurrentMatrix_[setCurrentMatrixCounter_++];

    return( fei->setCurrentMatrix(matID) );
  }

  if (!std::strcmp("setCurrentRHS", method)) {
    if (setCurrentRHSCounter_ >= (int)setCurrentRHS_.size()) {
      fei::console_out() << "driverData ERROR, can't call setCurrentRHS again" << FEI_ENDL;
      return(-1);
    }

    int rhsID = setCurrentRHS_[setCurrentRHSCounter_++];

    return( fei->setCurrentMatrix(rhsID) );
  }

  if (!std::strcmp("initFields", method)) {
    return( fei->initFields(initFields_numFields_,
			       initFields_fieldSizes_,
			       initFields_fieldIDs_) );
  }

  if (!std::strcmp("initElemBlock", method)) {
    return( fei->initElemBlock((GlobalID)initElemBlock_ints_[0],
				  initElemBlock_ints_[1],
				  initElemBlock_ints_[2],
				  initElemBlock_fieldsPerNode_,
				  initElemBlock_fieldIDs_,
				  initElemBlock_ints_[3],
				  initElemBlock_elemDofFieldIDs_,
				  initElemBlock_ints_[4]) );

  }

  if (!std::strcmp("parameters", method)) {
    if (parametersCounter_ >= (int)parameters_.size()) {
      fei::console_out() << "driverData ERROR, can't call parameters again" << FEI_ENDL;
      return(-1);
    }

    parameters* param = parameters_[parametersCounter_++];

    return( fei->parameters(param->paramList.size(),
			    &param->paramList[0]) );
  }

  if (!std::strcmp("initCRMult", method)) {
    if (initCRMultCounter_ >= (int)initCRMult_.size()) {
      fei::console_out() << "driverData ERROR, can't call initCRMult again" << FEI_ENDL;
      return(-1);
    }

    initCR* icr = initCRMult_[initCRMultCounter_++];

    return( fei->initCRMult(icr->numNodes, icr->nodeIDs,
			       icr->fieldIDs, icr->CRID) );
  }

  if (!std::strcmp("initSharedNodes", method)) {
    if (initSharedNodesCounter_ >= (int)initSharedNodes_.size()) {
      fei::console_out() << "driverData ERROR, can't call initSharedNodes again" << FEI_ENDL;
      return(-1);
    }

    sharedNodes* sn = initSharedNodes_[initSharedNodesCounter_++];

    return( fei->initSharedNodes(sn->numNodes, sn->nodeIDs,
			       sn->numProcsPerNode, sn->sharedProcIDs) );
  }

  if (!std::strcmp("loadCRMult", method)) {
    if (loadCRMultCounter_ >= (int)loadCRMult_.size()) {
      fei::console_out() << "driverData ERROR, can't call loadCRMult again" << FEI_ENDL;
      return(-1);
    }

    loadCR* lcr = loadCRMult_[loadCRMultCounter_++];

    return( fei->loadCRMult(lcr->CRID, lcr->numNodes, lcr->nodeIDs,
			       lcr->fieldIDs, lcr->weights, lcr->CRValue) );
  }

  if (!std::strcmp("deleteMultCRs", method)) {
    return( fei->deleteMultCRs() );
  }

  if (!std::strcmp("initElem", method)) {
    if (initElemCounter_ >= (int)initElems_.size()) {
      fei::console_out() << "driverData ERROR, can't call initElem again" << FEI_ENDL;
      return(-1);
    }

    initElem* ie = initElems_[initElemCounter_++];

    return( fei->initElem(ie->elemBlockID, ie->elemID, ie->nodeIDs) );
  }

  if (!std::strcmp("initComplete", method)) {
    return( fei->initComplete() );
  }

  if (!std::strcmp("resetSystem", method)) {
    return( fei->resetSystem(resetSystem_) );
  }

  if (!std::strcmp("resetMatrix", method)) {
    return( fei->resetMatrix(resetMatrix_) );
  }

  if (!std::strcmp("resetRHSVector", method)) {
    return( fei->resetRHSVector(resetRHSVector_) );
  }

  if (!std::strcmp("resetInitialGuess", method)) {
    return( fei->resetInitialGuess(resetInitialGuess_) );
  }

  if (!std::strcmp("sumInElem", method)) {
    if (sumInElemCounter_ >= (int)sumInElems_.size()) {
      fei::console_out() << "driverData ERROR, can't call sumInElem again" << FEI_ENDL;
      return(-1);
    }

    sumInElem* sie = sumInElems_[sumInElemCounter_++];

    return( fei->sumInElem(sie->elemBlockID, sie->elemID, sie->nodeIDs,
			      sie->stiffness, sie->load, sie->elemFormat) );
  }

  if (!std::strcmp("sumInElemMatrix", method)) {
    if (sumInElemMatrixCounter_ >= (int)sumInElemMatrix_.size()) {
      fei::console_out() << "driverData ERROR, can't call sumInElemMatrix again" << FEI_ENDL;
      return(-1);
    }

    sumInElem* sie = sumInElemMatrix_[sumInElemMatrixCounter_++];

    return( fei->sumInElemMatrix(sie->elemBlockID, sie->elemID, sie->nodeIDs,
			      sie->stiffness, sie->elemFormat) );
  }

  if (!std::strcmp("sumInElemRHS", method)) {
    if (sumInElemRHSCounter_ >= (int)sumInElemRHS_.size()) {
      fei::console_out() << "driverData ERROR, can't call sumInElemRHS again" << FEI_ENDL;
      return(-1);
    }

    sumInElem* sie = sumInElemRHS_[sumInElemRHSCounter_++];

    return( fei->sumInElemRHS(sie->elemBlockID, sie->elemID, sie->nodeIDs,
			         sie->load) );
  }

  if (!std::strcmp("putBlockFieldNodeSolution", method)) {
    if (putBlockFieldNodeSolutionCounter_ >=
	(int)putBlockFieldNodeSolution_.size()) {
      fei::console_out() << "driverData ERROR, can't call putBlockFieldNodeSolution again"
	   << FEI_ENDL;
      return(-1);
    }

    putBlockFieldNodeSolution* pbfns =
      putBlockFieldNodeSolution_[putBlockFieldNodeSolutionCounter_++];

    return( fei->putBlockFieldNodeSolution(pbfns->elemBlockID,
					      pbfns->fieldID,
					      pbfns->numNodes,
					      pbfns->nodeIDs,
					      pbfns->estimates) );
  }

  if (!std::strcmp("loadNodeBCs", method)) {
    if (loadNodeBCsCounter_ >= (int)loadNodeBCs_.size()) {
      fei::console_out() << "driverData ERROR, can't call loadNodeBCs again" << FEI_ENDL;
      return(-1);
    }

    fei::console_out() << "driverData: ERROR, loadNodeBCs needs to be re-examined..." << FEI_ENDL;
    return( -1 );
  }

  if (!std::strcmp("loadComplete", method)) {
    return( fei->loadComplete() );
  }

  if (!std::strcmp("solve", method)) {
    int status;
    return( fei->solve(status) );
  }

  if (!std::strcmp("getBlockNodeIDList", method) ||
      !std::strcmp("residualNorm",       method) ||
      !std::strcmp("getBlockFieldNodeSolution", method)) {
    return(0);
  }

  fei::console_out() << "driverData: ERROR unrecognized method name '" << method << "'"<<FEI_ENDL;
  return(1);
}

int driverData::readData(FEI_ISTREAM* instr, char* keyword)
{
  if (!std::strcmp("setSolveType", keyword)) {
    CHK_ERR( appendName(keyword) );
    return(readData(instr, solveType_));
  }

  if (!std::strcmp("setIDLists", keyword)) {
    int numMatrices = 0;
    CHK_ERR( readData(instr, numMatrices) );
    setIDLists* sidl = new setIDLists;
    sidl->numMatrices = numMatrices;
    sidl->matrixIDs = new int[numMatrices];
    int i;
    for(i=0; i<numMatrices; ++i) {
      CHK_ERR( readData(instr, sidl->matrixIDs[i]) );
    }
    int numRHSs = 0;
    CHK_ERR( readData(instr, numRHSs) );
    sidl->numRHSs = numRHSs;
    sidl->rhsIDs = new int[numRHSs];
    for(i=0; i<numRHSs; ++i) {
      CHK_ERR( readData(instr, sidl->rhsIDs[i]) );
    }

    setIDLists_.push_back(sidl);
    return( appendName(keyword) );
  }

  if (!std::strcmp("setCurrentMatrix", keyword)) {
    int matID = 0;
    CHK_ERR( readData(instr, matID) );
    setCurrentMatrix_.push_back(matID);
    return( appendName(keyword) );
  }

  if (!std::strcmp("setCurrentRHS", keyword)) {
    int rhsID = 0;
    CHK_ERR( readData(instr, rhsID) );
    setCurrentRHS_.push_back(rhsID);
    return( appendName(keyword) );
  }

  if (!std::strcmp("initFields", keyword)) {
    int i;
    CHK_ERR( readData(instr, initFields_numFields_) );
    initFields_fieldSizes_ = new int[initFields_numFields_];
    initFields_fieldIDs_ = new int[initFields_numFields_];

    for(i=0; i<initFields_numFields_; ++i) {
      CHK_ERR( readData(instr, initFields_fieldSizes_[i]) );
    } 
    for(i=0; i<initFields_numFields_; ++i) {
      CHK_ERR( readData(instr, initFields_fieldIDs_[i]) );
    }

    return( appendName(keyword) );
  }

  if (!std::strcmp("parameters", keyword)) {
    int numParams = 0;
    CHK_ERR( readData(instr, numParams) );
    parameters* param = new parameters;
    param->paramList.resize(numParams);
    CHK_ERR( skipWhite(instr) );
    for(int i=0; i<numParams; ++i) {
      char* line = new char[512];
      instr->getline(line, 512);
      param->paramList[i] = line;
    }
    parameters_.push_back(param);
    return( appendName(keyword) );
  }

  if (!std::strcmp("initElemBlock", keyword)) {
    initElemBlock_numInts_ = 5;
    int i, intOffset = 0;
    initElemBlock_ints_ = new int[initElemBlock_numInts_];
    //elemBlockID
    CHK_ERR( readData(instr, initElemBlock_ints_[intOffset++]) );
    //numElements
    CHK_ERR( readData(instr, initElemBlock_ints_[intOffset++]) );
    //numNodesPerElement
    CHK_ERR( readData(instr, initElemBlock_ints_[intOffset++]) );
    //now loop and read numFieldsPerNode
    int len = initElemBlock_ints_[intOffset-1];
    initElemBlock_fieldsPerNode_ = new int[len];
    initElemBlock_fieldIDs_ = new int*[len];
    for(i=0; i<len; ++i) {
      CHK_ERR( readData(instr, initElemBlock_fieldsPerNode_[i]) );
    }
    //now double-loop and read nodalFieldIDs
    for(i=0; i<len; ++i) {
      int len2 = initElemBlock_fieldsPerNode_[i];
      initElemBlock_fieldIDs_[i] = new int[len2];
      for(int ii=0; ii<len2; ++ii) {
	CHK_ERR( readData(instr, initElemBlock_fieldIDs_[i][ii]) );
      }
    }
    //numElemDOFPerElement
    CHK_ERR( readData(instr, initElemBlock_ints_[intOffset++]) );
    //now loop and read elemDOFFieldIDs
    len = initElemBlock_ints_[intOffset-1];
    if (len > 0) {
      initElemBlock_elemDofFieldIDs_ = new int[len];
      for(i=0; i<len; ++i) {
	CHK_ERR( readData(instr, initElemBlock_elemDofFieldIDs_[i]) );
      }
    }
    //interleaveStrategy
    CHK_ERR( readData(instr, initElemBlock_ints_[intOffset++]) );
    return( appendName(keyword) );
  }

  if (!std::strcmp("initElem", keyword) ) {
    initElem* ie = new initElem;
    int tmp;
    CHK_ERR( readData(instr, tmp) );
    ie->elemBlockID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    ie->elemID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    ie->nodeIDs = new GlobalID[tmp];
    ie->numNodes = tmp;
    for(int i=0; i<ie->numNodes; ++i) {
      CHK_ERR( readData(instr, tmp) );
      ie->nodeIDs[i] = (GlobalID)tmp;
    }
    initElems_.push_back(ie);
    return( appendName(keyword) );
  }

  if (!std::strcmp("initCRMult", keyword) ) {
    initCR* icr = new initCR;
    CHK_ERR( readData(instr, icr->numNodes) );
    if (icr->numNodes > 0) {
      icr->nodeIDs = new GlobalID[icr->numNodes];
      icr->fieldIDs = new int[icr->numNodes];
      int i, tmp;
      //read the nodeIDs
      for(i=0; i<icr->numNodes; ++i) {
	CHK_ERR( readData(instr, tmp) ); icr->nodeIDs[i] = (GlobalID)tmp;
      }
      //read the fieldIDs
      for(i=0; i<icr->numNodes; ++i) {
	CHK_ERR( readData(instr, icr->fieldIDs[i]) );
      }
    }
    //read the CRID
    CHK_ERR( readData(instr, icr->CRID) );

    initCRMult_.push_back(icr);
    return( appendName(keyword) );
  }

  if (!std::strcmp("loadCRMult", keyword) ) {
    loadCR* lcr = new loadCR;
    CHK_ERR( readData(instr, lcr->numNodes) );
    if (lcr->numNodes > 0) {
      lcr->nodeIDs = new GlobalID[lcr->numNodes];
      lcr->fieldIDs = new int[lcr->numNodes];
      lcr->fieldSizes = new int[lcr->numNodes];
      int i, tmp;
      //read the nodeIDs
      for(i=0; i<lcr->numNodes; ++i) {
	CHK_ERR( readData(instr, tmp) ); lcr->nodeIDs[i] = (GlobalID)tmp;
      }
      //read the fieldIDs
      for(i=0; i<lcr->numNodes; ++i) {
	CHK_ERR( readData(instr, lcr->fieldIDs[i]) );
      }
      //read the field-sizes
      tmp = 0;
      for(i=0; i<lcr->numNodes; ++i) {
	CHK_ERR( readData(instr, lcr->fieldSizes[i]) );
	tmp += lcr->fieldSizes[i];
      }
      //read the weights
      lcr->weights = new double[tmp];
      int offset = 0;
      for(i=0; i<lcr->numNodes; ++i) {
	int size = lcr->fieldSizes[i];
	for(int j=0; j<size; ++j) {
	  CHK_ERR( readData(instr, lcr->weights[offset++]) );
	}
      }
      //read the CRValue
      CHK_ERR( readData(instr, lcr->CRValue) );
    }
    //read the CRID
    CHK_ERR( readData(instr, lcr->CRID) );

    loadCRMult_.push_back(lcr);
    return( appendName(keyword) );
  }

  if (!std::strcmp("deleteMultCRs", keyword) ) {
    return( appendName(keyword) );
  }

  if (!std::strcmp("initSharedNodes", keyword) ) {
    sharedNodes* sn = new sharedNodes;
    CHK_ERR( readData(instr, sn->numNodes) );
    if (sn->numNodes > 0) {
      sn->nodeIDs = new GlobalID[sn->numNodes];
      sn->numProcsPerNode = new int[sn->numNodes];
      sn->sharedProcIDs = new int*[sn->numNodes];
      int i, tmp;
      //read the numProcsPerNode list
      for(i=0; i<sn->numNodes; ++i) {
	CHK_ERR( readData(instr, sn->numProcsPerNode[i]) );
	sn->sharedProcIDs[i] = new int[sn->numProcsPerNode[i]];
      }
      //read the nodeIDs and sharing-proc-ids
      for(i=0; i<sn->numNodes; ++i) {
	CHK_ERR( readData(instr, tmp) ); sn->nodeIDs[i] = (GlobalID)tmp;
	for(int j=0; j<sn->numProcsPerNode[i]; ++j) {
	  CHK_ERR( readData(instr, sn->sharedProcIDs[i][j]) );
	}
      }
    }

    initSharedNodes_.push_back(sn);
    return( appendName(keyword) );
  }

  if (!std::strcmp("initComplete", keyword) ) {
    return( appendName(keyword) );
  }

  if (!std::strcmp("sumInElem", keyword) ) {
    sumInElem* sie = new sumInElem;
    int tmp;
    double dtmp;
    CHK_ERR( readData(instr, tmp) );
    sie->elemBlockID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    sie->elemID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    sie->nodeIDs = new GlobalID[tmp];
    sie->numNodes = tmp;
    int i;
    for(i=0; i<sie->numNodes; ++i) {
      CHK_ERR( readData(instr, tmp) );
      sie->nodeIDs[i] = (GlobalID)tmp;
    }

    CHK_ERR( readData(instr, tmp) );
    sie->numRows = tmp;
    sie->stiff1D = new double[tmp*tmp];
    sie->load = new double[tmp];
    sie->stiffness = new double*[tmp];
    int offset = 0;
    for(i=0; i<sie->numRows; ++i) {
      for(int j=0; j<sie->numRows; ++j) {
	CHK_ERR( readData(instr, dtmp) );
	sie->stiff1D[offset++] = dtmp;
      }
      sie->stiffness[i] = &(sie->stiff1D[i*sie->numRows]);
    }

    for(int j=0; j<sie->numRows; ++j) {
      CHK_ERR( readData(instr, dtmp) );
      sie->load[j] = dtmp;
    }

    CHK_ERR( readData(instr, tmp) );
    sie->elemFormat = tmp;

    sumInElems_.push_back(sie);
    return( appendName(keyword) );
  }

  if (!std::strcmp("sumInElemMatrix", keyword) ) {
    sumInElem* sie = new sumInElem;
    int tmp;
    double dtmp;
    CHK_ERR( readData(instr, tmp) );
    sie->elemBlockID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    sie->elemID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    sie->nodeIDs = new GlobalID[tmp];
    sie->numNodes = tmp;
    int i;
    for(i=0; i<sie->numNodes; ++i) {
      CHK_ERR( readData(instr, tmp) );
      sie->nodeIDs[i] = (GlobalID)tmp;
    }

    CHK_ERR( readData(instr, tmp) );
    sie->numRows = tmp;
    sie->stiff1D = new double[tmp*tmp];
    sie->load = new double[tmp];
    sie->stiffness = new double*[tmp];
    int offset = 0;
    for(i=0; i<sie->numRows; ++i) {
      for(int j=0; j<sie->numRows; ++j) {
	CHK_ERR( readData(instr, dtmp) );
	sie->stiff1D[offset++] = dtmp;
      }
      sie->stiffness[i] = &(sie->stiff1D[i*sie->numRows]);
    }

    CHK_ERR( readData(instr, tmp) );
    sie->elemFormat = tmp;

    sumInElemMatrix_.push_back(sie);
    return( appendName(keyword) );
  }

  if (!std::strcmp("sumInElemRHS", keyword) ) {
    sumInElem* sie = new sumInElem;
    int tmp;
    double dtmp;
    CHK_ERR( readData(instr, tmp) );
    sie->elemBlockID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    sie->elemID = (GlobalID)tmp;
    CHK_ERR( readData(instr, tmp) );
    sie->nodeIDs = new GlobalID[tmp];
    sie->numNodes = tmp;
    for(int i=0; i<sie->numNodes; ++i) {
      CHK_ERR( readData(instr, tmp) );
      sie->nodeIDs[i] = (GlobalID)tmp;
    }

    CHK_ERR( readData(instr, tmp) );
    sie->numRows = tmp;

    sie->load = new double[sie->numRows];
    for(int j=0; j<sie->numRows; ++j) {
      CHK_ERR( readData(instr, dtmp) );
      sie->load[j] = dtmp;
    }

    sumInElemRHS_.push_back(sie);
    return( appendName(keyword) );
  }

  if (!std::strcmp("resetSystem", keyword) ) {
    CHK_ERR( readData(instr, resetSystem_) );
    return( appendName(keyword) );
  }

  if (!std::strcmp("resetMatrix", keyword) ) {
    CHK_ERR( readData(instr, resetMatrix_) );
    return( appendName(keyword) );
  }

  if (!std::strcmp("resetRHSVector", keyword) ) {
    CHK_ERR( readData(instr, resetRHSVector_) );
    return( appendName(keyword) );
  }

  if (!std::strcmp("resetInitialGuess", keyword) ) {
    CHK_ERR( readData(instr, resetInitialGuess_) );
    return( appendName(keyword) );
  }

  if (!std::strcmp("putBlockFieldNodeSolution", keyword) ) {
    putBlockFieldNodeSolution* pbfns = new putBlockFieldNodeSolution;
    CHK_ERR( readData(instr, pbfns->elemBlockID) );
    CHK_ERR( readData(instr, pbfns->fieldID) );
    CHK_ERR( readData(instr, pbfns->fieldSize) );
    CHK_ERR( readData(instr, pbfns->numNodes) );
    if (pbfns->numNodes > 0) {
      pbfns->nodeIDs = new GlobalID[pbfns->numNodes];

      int i=0;
      for(i=0; i<pbfns->numNodes; ++i) {
	CHK_ERR( readData(instr, pbfns->nodeIDs[i]) );
      }
      int len = pbfns->numNodes * pbfns->fieldSize;
      pbfns->estimates = new double[len];

      for(i=0; i<pbfns->numNodes * pbfns->fieldSize; ++i) {
	CHK_ERR( readData(instr, pbfns->estimates[i]) );
      }
    }

    putBlockFieldNodeSolution_.push_back(pbfns);
    return( appendName(keyword) );
  }

  if (!std::strcmp("loadNodeBCs", keyword) ) {
    nodeBC* nbc = new nodeBC;
    CHK_ERR( readData(instr, nbc->numNodes) );
    CHK_ERR( readData(instr, nbc->fieldID) );
    CHK_ERR( readData(instr, nbc->fieldSize) );

    if (nbc->numNodes > 0) {
      nbc->nodeIDs = new GlobalID[nbc->numNodes];
      nbc->alpha = new double*[nbc->numNodes];
      nbc->beta = new double*[nbc->numNodes];
      nbc->gamma = new double*[nbc->numNodes];

      int i, j, tmp;
      for(i=0; i<nbc->numNodes; ++i) {
	nbc->alpha[i] = new double[nbc->fieldSize];
	nbc->beta[i] = new double[nbc->fieldSize];
	nbc->gamma[i] = new double[nbc->fieldSize];

	CHK_ERR( readData(instr, tmp) );
	nbc->nodeIDs[i] = (GlobalID)tmp;

	for(j=0; j<nbc->fieldSize; ++j) {
	  CHK_ERR( readData(instr, nbc->alpha[i][j]));
	}
	for(j=0; j<nbc->fieldSize; ++j) {
	  CHK_ERR( readData(instr, nbc->beta[i][j]));
	}
	for(j=0; j<nbc->fieldSize; ++j) {
	  CHK_ERR( readData(instr, nbc->gamma[i][j]));
	}
      }
    }

    loadNodeBCs_.push_back(nbc);
    return( appendName(keyword) );
  }

  if (!std::strcmp("loadComplete",       keyword) ||
      !std::strcmp("solve",              keyword) ||
      !std::strcmp("destructor",         keyword) ||
      !std::strcmp("getBlockNodeIDList", keyword) ||
      !std::strcmp("getBlockFieldNodeSolution", keyword) ||
      !std::strcmp("residualNorm",       keyword)) {
    return( appendName(keyword) );
  }

  return(-1);
}

int driverData::appendName(const char* name)
{
  if (name == NULL) return(-1);
  char* str = new char[strlen(name)+1];
  strcpy(str, name);
  methodNames.push_back(str);
  return(0);
}

int driverData::getKeyword(FEI_ISTREAM* instr, char*& keyword)
{
  int err = skipWhite(instr);
  if (err) return(err);

  for(int i=0; i<tempLen_; i++) temp_[i] = '\0';

  do {
    instr->getline(temp_, tempLen_);
  } while ((strlen(temp_) == 0) && (!instr->eof()));

  if (instr->eof() || strlen(temp_) == 0) return(-1);

  keyword = new char[strlen(temp_)+1];
  const char* temp2 = snl_fei::getParamValue("FEI:", 1, &temp_);

  if (temp2 != NULL) {
    strcpy(keyword, temp2);
    return(0);
  }

  return(-1);
}

//==============================================================================
int driverData::is_reg_char(char c) {
   int i = (int)c;
   if (i<1 || i>126) return(0);

   return(1);
}

//==============================================================================
int driverData::skipWhite(FEI_ISTREAM* instr) {
   char c = '\0';
   instr->get(c);

   if (!is_reg_char(c)) {
      return(-1);
   }

   while(c == '#' || c == '\n' || c == ' ') {
      if (c=='#') {
         char* buf = new char[128];
         for(int i=0; i<128; i++) buf[i] = '\0';
         instr->getline(buf, 128);
         delete [] buf;
      }

      instr->get(c);

      if (instr->eof()) return(1);
      if ((int)c == EOF) return(1);

      if (!is_reg_char(c)) {
         return(-1);
      }
   }

   instr->putback(c);
   return(0);
}

//==============================================================================
int driverData::readData(FEI_ISTREAM* instr, int& n) {
  int err = skipWhite(instr);
  if (err) return(err);
  (*instr) >> n;
  return(0);
}

//==============================================================================
int driverData::readData(FEI_ISTREAM* instr, double& val) {
  int err = skipWhite(instr);
  if (err) return(err);
  (*instr) >> val;
  return(0);
}


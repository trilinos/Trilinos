/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_sstream.hpp>
#include <fei_fstream.hpp>

#include <test_utils/fei_test_utils.hpp>

#include <test_utils/FEI_tester.hpp>

#include <fei_LinearSystemCore.hpp>
#include <fei_LibraryWrapper.hpp>
#include <snl_fei_Utils.hpp>

#include <fei_FEI_Impl.hpp>

#include <test_utils/LibraryFactory.hpp>

#ifdef HAVE_FEI_FETI
#include <FETI_DP_FiniteElementData.h>
#endif

#include <test_utils/DataReader.hpp>
#include <test_utils/SolnCheck.hpp>

#undef fei_file
#define fei_file "FEI_tester.cpp"

#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
FEI_tester::FEI_tester(fei::SharedPtr<DataReader> data_reader,
		       MPI_Comm comm, int localProc, int numProcs, bool useNewFEI)
  : comm_(comm),
    fei_(),
    wrapper_(),
    data_(data_reader),
    localProc_(localProc),
    numProcs_(numProcs),
    numMatrices(1),
    matrixIDs(NULL),
    numRHSs(1),
    rhsIDs(NULL),
    useNewFEI_(useNewFEI)
{
}

//----------------------------------------------------------------------------
FEI_tester::~FEI_tester()
{
  delete [] matrixIDs;
  delete [] rhsIDs;
}

//----------------------------------------------------------------------------
int FEI_tester::testInitialization()
{
  if (data_.get() == NULL) {
    ERReturn(-1);
  }

  CHK_ERR( createFEIinstance(data_->solverLibraryName_.c_str()) );

  const char* feiVersionString;
  CHK_ERR( fei_->version(feiVersionString) );

  FEI_COUT << "FEI version: " << feiVersionString << FEI_ENDL;

  fei_->parameters(data_->numParams_, data_->paramStrings_);

  //Now we check the solveType. A regular Ax=b solve corresponds to 
  //solveType==FEI_SINGLE_SYSTEM. The aggregate stuff (solveType==
  //FEI_AGGREGATE_SUM) is for power users who
  //want to assemble more than one matrix system and solve a linear
  //combination of them, an aggregate system.

  if (data_->solveType_ == FEI_AGGREGATE_SUM) {
    CHK_ERR( setIDlists());
  }

  CHK_ERR( initializationPhase() );

  int numBlkActNodes;
  for(int i=0; i<data_->numElemBlocks_; ++i) {
    ElemBlock& eblk = data_->elemBlocks_[i];
    int elemBlockID = eblk.blockID_;
    CHK_ERR( fei_->getNumBlockActNodes(elemBlockID, numBlkActNodes) );
  }

  //************** Initialization Phase is now complete *****************

  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::testLoading()
{
  CHK_ERR(fei_->resetRHSVector());
  CHK_ERR(fei_->resetMatrix());    //try out all of these reset functions,
  CHK_ERR(fei_->resetSystem());    //just for the sake of coverage...

  // ************ Load Phase (delegated to a function) ***************

  //obviously only one of these load-phases should be
  //performed.

  if (data_->solveType_ == FEI_SINGLE_SYSTEM) {
    CHK_ERR( normalLoadPhase());
  }
  if (data_->solveType_ == FEI_AGGREGATE_SUM) {
    CHK_ERR( aggregateLoadPhase());
  }

  CHK_ERR( fei_->loadComplete() );

  CHK_ERR( exerciseResidualNorm() );

  //**** let's try out the 'put' functions. ******************
  CHK_ERR( exercisePutFunctions() );

  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::testSolve()
{
  //If solverLibraryName is TEST_LSC, then we're not running a real solver so
  //just return.
  std::string sname(data_->solverLibraryName_);
  if (sname == "TEST_LSC") {
    return( 0 );
  }

  int status;
  int err = fei_->solve(status);

  //FEI_COUT << "fei->solve, err = " << err << ", status = " << status << FEI_ENDL;

  if (err != 0 || status != 0) {
    FEI_COUT << "!!!! solve returned: err: "<<err<<", status: "<<status<<FEI_ENDL;
    return(err);
  }

  if (localProc_ == 0) {
    int itersTaken;
    CHK_ERR( fei_->iterations(itersTaken));
    //FEI_COUT << "iterations: " << itersTaken << FEI_ENDL;
  }

  CHK_ERR( exerciseResidualNorm() );

  return(0);
}

//----------------------------------------------------------------------------
void FEI_tester::dumpMatrixFiles()
{
}

//----------------------------------------------------------------------------
void FEI_tester::setParameter(const char*)
{
}

//----------------------------------------------------------------------------
int FEI_tester::testCheckResult()
{
  //If solverLibraryName is TEST_LSC, then we're not running a real solver so
  //just check the matrix.
  std::string sname(data_->solverLibraryName_);
  if (sname == "TEST_LSC") {
    return( lsc_matrix_check() );
  }

  CHK_ERR( save_block_node_soln(*data_, *fei_, data_->solnFileName_.c_str(),
				numProcs_, localProc_, 1));

  CHK_ERR( save_block_elem_soln(*data_, *fei_, data_->solnFileName_.c_str(),
				numProcs_, localProc_, 1));

  CHK_ERR( save_multiplier_soln(*data_, *fei_, data_->solnFileName_.c_str(),
				numProcs_, localProc_, 1));

  int err = SolnCheck::checkSolution(localProc_, numProcs_, data_->solnFileName_.c_str(),
				     data_->checkFileName_.c_str(), "node", 1);

  err += SolnCheck::checkSolution(localProc_, numProcs_, data_->solnFileName_.c_str(),
				  data_->checkFileName_.c_str(), "elem", 1);

  err += SolnCheck::checkSolution(localProc_, numProcs_, data_->solnFileName_.c_str(),
				  data_->checkFileName_.c_str(), "mult", 1);
  int globalErr = err;
#ifndef FEI_SER
  if (MPI_SUCCESS != MPI_Allreduce(&err, &globalErr, 1, MPI_INT, MPI_SUM,
				   comm_)) return(-1);
#endif
  if (globalErr != 0) {
    ERReturn(-1);
  }

  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::createFEIinstance(const char* solverName)
{
  try {
    wrapper_ = fei::create_LibraryWrapper(comm_, solverName);
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what()<<FEI_ENDL;
    return(1);
  }

  if (wrapper_.get() == NULL) ERReturn(-1);

  if (useNewFEI_) {
    fei_.reset(new fei::FEI_Impl(wrapper_, comm_, 0));
  }
  else {
    fei_.reset(new FEI_Implementation(wrapper_, comm_, 0));
  }

  if (fei_.get() == NULL) ERReturn(-1);

  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::setIDlists()
{
  snl_fei::getIntParamValue("numMatrices",
			    data_->numParams_,
			    data_->paramStrings_,
			    numMatrices);

  matrixIDs = new int[numMatrices];
  numRHSs = 1;
  rhsIDs = new int[1];
  rhsIDs[0] = 0;

  for(int i=0; i<numMatrices; i++) {
    matrixIDs[i] = i;
  }

  CHK_ERR(fei_->setIDLists(numMatrices, matrixIDs, numRHSs, rhsIDs));
  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::initializationPhase()
{
  if (data_->solveType_ != FEI_SINGLE_SYSTEM &&
      data_->solveType_ != FEI_AGGREGATE_SUM) {
    FEI_COUT << "FEI_tester: bad solveType: " << data_->solveType_ << FEI_ENDL;
    return(-1);
  }

  CHK_ERR( fei_->setSolveType(data_->solveType_));

  CHK_ERR(fei_->initFields(data_->numFields_, data_->fieldSizes_, data_->fieldIDs_));

  int i;
  for(i=0; i<data_->numElemBlocks_; i++) {
    ElemBlock& block = data_->elemBlocks_[i];

    CHK_ERR(fei_->initElemBlock(block.blockID_,
			       block.numElements_,
			       block.numNodesPerElement_,
			       block.numFieldsPerNode_,
			       block.nodalFieldIDs_,
			       block.numElemDOF_,
			       block.elemDOFFieldIDs_,
			       block.interleaveStrategy_) );

    for(int el=0; el<block.numElements_; el++) {
      CHK_ERR(fei_->initElem(block.blockID_, 
			    block.elemIDs_[el],
			    block.elemConn_[el]));
    }
  }

  for(i=0; i<data_->numSharedNodeSets_; i++) {
    CommNodeSet& shNodeSet = data_->sharedNodeSets_[i];

    CHK_ERR(fei_->initSharedNodes(shNodeSet.numNodes_, shNodeSet.nodeIDs_,
				 shNodeSet.procsPerNode_, shNodeSet.procs_));
  }

  //********* Initialize any slave variables *****************************

  for(i=0; i<data_->numSlaveVars_; i++) {
    CRSet& crSet = data_->slaveVars_[i];

    CHK_ERR( fei_->initSlaveVariable(crSet.slaveNodeID_,
				    crSet.slaveFieldID_,
				    crSet.slaveOffset_,
				    crSet.numNodes_,
				    crSet.nodeIDs_[0],
				    crSet.fieldIDs_,
				    crSet.weights_,
				    crSet.values_[0]) );
  }

  //*********** Initialize Constraint Relation Equations *****************

  for(i=0; i<data_->numCRMultSets_; i++) {
    CRSet& crSet = data_->crMultSets_[i];

    for(int j=0; j<1; j++) {
      CHK_ERR(fei_->initCRMult(crSet.numNodes_,
			      crSet.nodeIDs_[j],
			      crSet.fieldIDs_,
			      crSet.crID_));
    }
  }

  for(i=0; i<data_->numCRPenSets_; i++) {
    CRSet& crSet = data_->crPenSets_[i];

    for(int j=0; j<1; j++) {
      CHK_ERR(fei_->initCRPen(crSet.numNodes_,
			     crSet.nodeIDs_[j],
			     crSet.fieldIDs_,
			     crSet.crID_));
    }
  }

  CHK_ERR(fei_->initComplete());
  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::normalLoadPhase()
{
  int i;

  //*************** Boundary Condition Node Sets *************************

  for(i=0; i<data_->numBCNodeSets_; i++) {
    BCNodeSet& bcSet = data_->bcNodeSets_[i];

    CHK_ERR(fei_->loadNodeBCs(bcSet.numNodes_,
			     bcSet.nodeIDs_,
			     bcSet.fieldID_,
			     bcSet.offsetsIntoField_,
			     bcSet.prescribed_values_));
  }

   for(i=0; i<data_->numElemBlocks_; i++) {
      ElemBlock& block = data_->elemBlocks_[i];

      for(int el=0; el<block.numElements_; el++) {

         CHK_ERR(fei_->sumInElemMatrix(block.blockID_,
                                block.elemIDs_[el],
                                block.elemConn_[el],
                                block.elemStiff_[el],
                                block.elemFormat_));
      }
   }

   for(i=0; i<data_->numElemBlocks_; i++) {
      ElemBlock& block = data_->elemBlocks_[i];

      for(int el=0; el<block.numElements_; el++) {

         CHK_ERR(fei_->sumInElemRHS(block.blockID_,
                                block.elemIDs_[el],
                                block.elemConn_[el],
                                block.elemLoad_[el]));
      }
   }

   //******** Load Constraint Relation Equations ***********************

   for(i=0; i<data_->numCRMultSets_; i++) {
      CRSet& crSet = data_->crMultSets_[i];

      for(int j=0; j<1; j++) {
         CHK_ERR(fei_->loadCRMult(crSet.crID_,
                                 crSet.numNodes_,
                                 crSet.nodeIDs_[j],
                                 crSet.fieldIDs_,
                                 crSet.weights_,
                                 crSet.values_[j]))
      }
   }

   for(i=0; i<data_->numCRPenSets_; i++) {
      CRSet& crSet = data_->crPenSets_[i];

      for(int j=0; j<1; j++) {
         CHK_ERR(fei_->loadCRPen(crSet.crID_,
                                crSet.numNodes_,
                                crSet.nodeIDs_[j],
                                crSet.fieldIDs_,
                                crSet.weights_,
                                crSet.values_[j],
                                crSet.penValues_[j]))
      }
   }
  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::aggregateLoadPhase()
{
   int i;

   for(i=0; i<numMatrices; i++) {
      CHK_ERR(fei_->setCurrentMatrix(matrixIDs[i]))

      for(int j=0; j<data_->numElemBlocks_; j++) {
         ElemBlock& block = data_->elemBlocks_[j];

         for(int el=0; el<block.numElements_; el++) {

            CHK_ERR(fei_->sumInElemMatrix(block.blockID_,
                                   block.elemIDs_[el],
                                   block.elemConn_[el],
                                   block.elemStiff_[el],
                                   block.elemFormat_))
         }
      }
   }

   for(i=0; i<numRHSs; i++) {
      CHK_ERR(fei_->setCurrentRHS(rhsIDs[i]))

      for(int j=0; j<data_->numElemBlocks_; j++) {
         ElemBlock& block = data_->elemBlocks_[j];

         for(int el=0; el<block.numElements_; el++) {
            CHK_ERR(fei_->sumInElemRHS(block.blockID_,
                                   block.elemIDs_[el],
                                   block.elemConn_[el],
                                   block.elemLoad_[el]))
         }
      }
   }

   //*************** Boundary Condition Node Sets *************************

   for(i=0; i<data_->numBCNodeSets_; i++) {
      BCNodeSet& bcSet = data_->bcNodeSets_[i];

      CHK_ERR(fei_->loadNodeBCs(bcSet.numNodes_,
                     bcSet.nodeIDs_,
                     bcSet.fieldID_,
                     bcSet.offsetsIntoField_,
                     bcSet.prescribed_values_))
   }

   double* matScalars = new double[numMatrices];
   for(i=0; i<numMatrices; i++) {
      matScalars[i] = 1.0;
   }

   int rhsScaleID = rhsIDs[0];
   double rhsScalar = 1.0;

   CHK_ERR(fei_->setMatScalars(numMatrices, matrixIDs, matScalars))
   CHK_ERR(fei_->setRHSScalars(1, &rhsScaleID, &rhsScalar))

   delete [] matScalars;
  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::exerciseResidualNorm()
{
  std::string sname(data_->solverLibraryName_);
  if (sname == "TEST_LSC") {
    return( 0 );
  }
  //FEI_COUT << "numFields: " << data_->numFields_ << FEI_ENDL;
  double* norms = new double[data_->numFields_];
  int *fields = new int[data_->numFields_];
  for(int i=0; i<data_->numFields_; ++i) {
    fields[i] = data_->fieldIDs_[i];
  }

  CHK_ERR( fei_->residualNorm(1, data_->numFields_, fields, norms) );
  //int n;
  //for(n=0; n<data_->numFields_; n++) {
  //  FEI_COUT << " field["<<n<<"]: " << fields[n]
  //    << ", 1-norm: " << norms[n] << FEI_ENDL;
  //}

  delete [] fields;
  delete [] norms;
  return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::exercisePutFunctions()
{
   //let's exercise the putNodalFieldData function.

   int numNodes;
   CHK_ERR( fei_->getNumLocalNodes(numNodes) );
   std::vector<int> nodeIDs(numNodes);
   int* nodeIDsPtr = &nodeIDs[0];
   int checkNumNodes;
   CHK_ERR( fei_->getLocalNodeIDList(checkNumNodes, nodeIDsPtr, numNodes) );

   for(int i=0; i<data_->numFields_; ++i) {
     int fieldID = data_->fieldIDs_[i];
     int fieldSize = data_->fieldSizes_[i];
     std::vector<double> data(numNodes*fieldSize, 0.0001);

     CHK_ERR( fei_->putNodalFieldData(fieldID, numNodes, nodeIDsPtr,
				      &data[0]) );
   }

   return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::save_block_node_soln(DataReader& data, FEI& fei,
				     const char* solnFileName, int numProcs,
				     int localProc, int solveCounter)
{
  (void)solveCounter;
   int i;

   int maxNumEqnsPerNode = 0;
   for(i=0; i<data.numFields_; ++i) {
     maxNumEqnsPerNode += data.fieldSizes_[i];
   }

   std::vector<double> soln(maxNumEqnsPerNode);

   int numNodes;
   CHK_ERR( fei.getNumLocalNodes(numNodes) );

   std::vector<GlobalID> nodes(numNodes);
   int* nodesPtr = &nodes[0];

   int checkNumNodes;
   CHK_ERR( fei.getLocalNodeIDList( checkNumNodes, nodesPtr, numNodes) );

   if (checkNumNodes != numNodes) {
     ERReturn(-1);
   }

   FEI_OSTRINGSTREAM fileName;
   fileName << solnFileName<<".node."<<solveCounter<<"."<<numProcs<<"."<<localProc;
   FEI_OFSTREAM outfile(fileName.str().c_str());

   if (!outfile || outfile.bad()) {
      FEI_COUT << "ERROR opening solution output file " << fileName.str() << FEI_ENDL;
      return(1);
   }

   outfile.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);

   std::vector<int> offsets(2);

   for(i=0; i<numNodes; ++i) {
     CHK_ERR( fei.getNodalSolution(1, &(nodesPtr[i]),
				   &offsets[0], &soln[0]) );

     int numDOF = offsets[1];

     outfile << nodesPtr[i] << " " << numDOF << FEI_ENDL;
     for(int j=0; j<numDOF; j++) {
       outfile << soln[j] << " ";
     }
     outfile << FEI_ENDL;
   }

   return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::save_block_elem_soln(DataReader& data, FEI& fei,
				     const char* solnFileName,
				     int numProcs, int localProc,
				     int solveCounter)
{
   int returnValue = 0;
   FEI_OSTRINGSTREAM fileName;
   fileName << solnFileName<<".elem."<<solveCounter<<"."<<numProcs<<"."<<localProc;
   FEI_OFSTREAM outfile(fileName.str().c_str());

   if (!outfile || outfile.bad()) {
      FEI_COUT << "ERROR opening elem-solution output file " << fileName.str() << FEI_ENDL;
      return(1);
   }

   for(int i=0; i<data.numElemBlocks_; i++) {
      if (returnValue != 0) break;

      ElemBlock& eb = data.elemBlocks_[i];

      GlobalID blockID = eb.blockID_;
      int numElems;
      CHK_ERR( fei.getNumBlockElements(blockID, numElems))

      int dofPerElem;
      CHK_ERR( fei.getNumBlockElemDOF(blockID, dofPerElem))
      int totalNumElemDOF = numElems*dofPerElem;

      if (totalNumElemDOF < 1) {
	continue;
      }

      GlobalID* elemIDs = new GlobalID[numElems];
      if (elemIDs==NULL) return(-1);

      int err = fei.getBlockElemIDList(blockID, numElems, elemIDs);
      if (err) returnValue = 1;

      int* offsets = new int[numElems+1];
      if (offsets == NULL) return(-1);

      if (totalNumElemDOF > 0) {
         double* solnValues = new double[totalNumElemDOF];
         if (solnValues == NULL) return(-1);

         err = fei.getBlockElemSolution(blockID, numElems, elemIDs,
                                         dofPerElem, solnValues);
         if (err) returnValue = 1;

         if (!err) {
            for(int j=0; j<numElems; j++) {

               outfile << (int)elemIDs[j] << " " << dofPerElem << FEI_ENDL << "  ";
               for(int k=0; k<dofPerElem; k++) {
                  outfile << solnValues[j*dofPerElem + k] << " ";
               }
               outfile << FEI_ENDL;
            }
         }

         delete [] solnValues;
      }

      delete [] elemIDs;
      delete [] offsets;
   }

   return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::save_multiplier_soln(DataReader& data, FEI& fei,
			 const char* solnFileName,
                          int numProcs, int localProc, int solveCounter)
{
   int numCRs = 0;

   CHK_ERR( fei.getNumCRMultipliers(numCRs) );

   int* globalNumCRs = new int[numProcs];
#ifndef FEI_SER
   if (MPI_Allgather(&numCRs, 1, MPI_INT, globalNumCRs, 1, MPI_INT,
		 comm_) != MPI_SUCCESS) {
     ERReturn(-1);
   }
#endif

   int localCRStart = 0;
   for(int p=0; p<localProc; p++) localCRStart += globalNumCRs[p];

   delete [] globalNumCRs;

   FEI_OSTRINGSTREAM fileName;
   fileName << solnFileName<<".mult."<<solveCounter<<"."<<numProcs<<"."<<localProc;
   FEI_OFSTREAM outfile(fileName.str().c_str());

   if (!outfile || outfile.bad()) {
      FEI_COUT << "ERROR opening mult-solution output file " << fileName.str() << FEI_ENDL;
      return(-1);
   }

   int* CRIDs = numCRs > 0 ? new int[numCRs] : NULL;
   double* results = numCRs > 0 ? new double[numCRs] : NULL;

   if (numCRs > 0 && (CRIDs==NULL || results==NULL)) {
     ERReturn(-1);
   }

   if (numCRs < 1) {
     return(0);
   }

   CHK_ERR( fei.getCRMultIDList(numCRs, CRIDs) );

   std::string sname(data_->solverLibraryName_);
   if (sname == "FETI") {
     for(int ii=0; ii<numCRs; ++ii) results[ii] = -999.99;
   }
   else {
     CHK_ERR( fei.getCRMultipliers(numCRs, CRIDs, results));
   }

   for(int i=0; i<numCRs; i++) {
      outfile << localCRStart++ << " " << 1 << FEI_ENDL;

      outfile << "   " << results[i] << FEI_ENDL;
   }

   delete [] CRIDs;
   delete [] results;

   return(0);
}

//----------------------------------------------------------------------------
int FEI_tester::lsc_matrix_check()
{
   if (localProc_ == 0) {
     char* current_dir = NULL;
     CHK_ERR( fei_test_utils::dirname(data_->solnFileName_.c_str(), current_dir));

     FEI_OSTRINGSTREAM solnMtxName;
     solnMtxName<< current_dir<<"/A_TLSC.mtx";
     fei::FillableMat solnMtx, checkMtx;
     CHK_ERR( SolnCheck::readMatrix(solnMtxName.str().c_str(), numProcs_, solnMtx) );
     CHK_ERR( SolnCheck::readMatrix(data_->checkFileName_.c_str(), numProcs_, checkMtx) );
     int err = SolnCheck::compareMatrices(solnMtx, checkMtx);
     delete [] current_dir;
     if (err == 0) {
       FEI_COUT << "Utst_fei_lsc: TEST PASSED" << FEI_ENDL;
     }
     else {
       FEI_COUT << "Utst_fei_lsc: TEST FAILED" << FEI_ENDL;
       return(-1);
     }
   }

   return(0);
}

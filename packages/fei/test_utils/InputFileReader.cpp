/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <test_utils/InputFileReader.hpp>

#include <fei_sstream.hpp>
#include <fei_fstream.hpp>
#include <fei_Pattern.hpp>
#include <fei_utils.hpp>
#include "fei_ConnectivityBlock.hpp"

#include <test_utils/DataReader.hpp>
#include <test_utils/fei_test_utils.hpp>
#include <test_utils/SolnCheck.hpp>

//--------------------------------------------------------------------
InputFileReader::InputFileReader(MPI_Comm comm,
                                 fei::SharedPtr<fei::Factory> factory)
 : comm_(comm),
   factory_(factory),
   feiVectorSpace_row_(),
   feiMatrixGraph_(),
   initCompleteCalled_(false),
   feiMatrix_(),
   feiVector_soln_(),
   feiVector_rhs_(),
   feiLinearSystem_(),
   solnFileName_(),
   checkFileName_(),
   solveCounter_(0),
   lagrangeConstraints_(),
   crID_(0)
{
  int localProc = 0;
#ifndef FEI_SER
  MPI_Comm_rank(comm_, &localProc);
#endif

  crID_ = localProc*1000;
}

//--------------------------------------------------------------------
InputFileReader::~InputFileReader() {
}

//--------------------------------------------------------------------
fei::SharedPtr<fei::VectorSpace>& InputFileReader::getRowSpace()
{
  if (feiVectorSpace_row_.get() == NULL) {
    feiVectorSpace_row_ = factory_->createVectorSpace(comm_, NULL);
    if (feiParameterSet.get() != NULL) {
      feiVectorSpace_row_->setParameters(*feiParameterSet);
    }
  }

  return(feiVectorSpace_row_);
}

//--------------------------------------------------------------------
fei::SharedPtr<fei::MatrixGraph>& InputFileReader::getMatrixGraph()
{
  if (feiMatrixGraph_.get() == NULL) {
    feiMatrixGraph_ = factory_->createMatrixGraph(getRowSpace(),
						  feiVectorSpace_col,
						  NULL);
    if (feiParameterSet.get() != NULL) {
      feiMatrixGraph_->setParameters(*feiParameterSet);
    }
  }

  return(feiMatrixGraph_);
}

//--------------------------------------------------------------------
fei::SharedPtr<fei::Matrix>& InputFileReader::getMatrix()
{
  if (feiMatrix_.get() != NULL) {
    return(feiMatrix_);
  }

  if (feiMatrixGraph_.get() == NULL) {
    throw std::runtime_error("InputFileReader::getMatrix: don't yet have MatrixGraph");
  }
  else {
    if (!initCompleteCalled_) {
      feiMatrixGraph_->initComplete();
      initCompleteCalled_ = true;
    }

    feiMatrix_ = factory_->createMatrix(feiMatrixGraph_);
  }

  return(feiMatrix_);
}

//--------------------------------------------------------------------
fei::SharedPtr<fei::LinearSystem>& InputFileReader::getLinearSystem()
{
  if (feiLinearSystem_.get() != NULL) {
    return(feiLinearSystem_);
  }

  fei::SharedPtr<fei::Vector> x = getVector(true);
  fei::SharedPtr<fei::Vector> b = getVector(false);
  fei::SharedPtr<fei::Matrix> A = getMatrix();

  fei::SharedPtr<fei::MatrixGraph>& matgraph = getMatrixGraph();
  fei::SharedPtr<fei::LinearSystem> linsys =
    factory_->createLinearSystem(matgraph);

  if (feiParameterSet.get() != NULL) {
    linsys->parameters(*feiParameterSet);
  }

  linsys->setMatrix(A);
  linsys->setSolutionVector(x);
  linsys->setRHS(b);

  feiLinearSystem_ = linsys;

  return(feiLinearSystem_);
}

//--------------------------------------------------------------------
fei::SharedPtr<fei::Vector>& InputFileReader::getVector(bool soln)
{
  if (soln) {
    if (feiVector_soln_.get() != NULL) return(feiVector_soln_);
  }
  else {
    if (feiVector_rhs_.get() != NULL) return(feiVector_rhs_);
  }

  if (feiMatrixGraph_.get() == NULL) {
    if (soln) {
      feiVector_soln_ = factory_->createVector(getRowSpace(), soln, 1);
      feiVector_soln_->putScalar(0.0);
    }
    else {
      feiVector_rhs_ = factory_->createVector(getRowSpace(), soln, 1);
      feiVector_rhs_->putScalar(0.0);
    }
  }
  else {
    if (!initCompleteCalled_) {
      feiMatrixGraph_->initComplete();
      initCompleteCalled_ = true;
    }

    if (soln) {
      feiVector_soln_ = factory_->createVector(feiMatrixGraph_, soln, 1);
      feiVector_soln_->putScalar(0.0);
    }
    else {
      feiVector_rhs_ = factory_->createVector(feiMatrixGraph_, soln, 1);
      feiVector_rhs_->putScalar(0.0);
    }
  }

  if (soln) return(feiVector_soln_);
  else return(feiVector_rhs_);
}

//--------------------------------------------------------------------
int InputFileReader::readInputFile(const char* fileName)
{
  int numProcs = 1;
  int localProc = 0;
#ifndef FEI_SER
  MPI_Comm_size(comm_, &numProcs);
  MPI_Comm_rank(comm_, &localProc);
#endif

  FEI_OSTRINGSTREAM osstr;
  osstr << fileName << "." << numProcs << "." << localProc;
  std::string fullName = osstr.str();

  FEI_IFSTREAM* instr = new FEI_IFSTREAM(fullName.c_str());

  if (instr->bad()) {
    FEI_CERR << "InputFileReader::readInputFile: ERROR opening " << fileName
	     << FEI_ENDL;
    return(-1);
  }

  char* keyword = NULL;

  int err = DataReader::getKeyword(instr, keyword);

  while (!instr->eof() && !err) {
    readData(instr, keyword);
    delete [] keyword;
    keyword = NULL;

    err = DataReader::getKeyword(instr, keyword);
  }

  delete instr;
  delete [] keyword;

  return(0);
}

//--------------------------------------------------------------------
int InputFileReader::readParameters(const char* fileName)
{
  std::vector<std::string> file_strings;
  fei_test_utils::read_input_file(fileName, comm_, file_strings);

  feiParameterSet.reset(new fei::ParameterSet);

  fei::utils::parse_strings(file_strings, " ", *feiParameterSet);
  return(0);
}

//--------------------------------------------------------------------
void InputFileReader::readData(FEI_ISTREAM* instr, const char* keyword)
{
  std::string keystr(keyword);
  bool keyword_recognized = false;

  if (keystr == "FIELD")
    { readField(instr); keyword_recognized = true;}

  if (keystr == "IDTYPE")
    { readIDType(instr); keyword_recognized = true;}

  if (keystr == "PATTERN_SIMPLE")
    { readSimplePattern(instr); keyword_recognized = true;}

  if (keystr == "ELEM_BLOCK_SYMM")
    { readSymmElemBlock(instr); keyword_recognized = true;}

  if (keystr == "CONN_BLOCK_NONSYMM")
    { readNonsymmConnBlock(instr); keyword_recognized = true;}

  if (keystr == "ELEM_CONN")
    { readSymmConnectivity(instr); keyword_recognized = true;}

  if (keystr == "SLAVE_CONSTRAINT")
    { readSlaveConstraint(instr); keyword_recognized = true;}

  if (keystr == "LAGRANGE_CONSTRAINT")
    { readLagrangeConstraint(instr); keyword_recognized = true;}

  if (keystr == "ELEM_STIFF")
    { readElemStiffness(instr); keyword_recognized = true;}

  if (keystr == "ELEM_LOAD")
    { readElemLoad(instr); keyword_recognized = true;}

  if (keystr == "NONSYMM_CONN")
    { readNonsymmConnectivity(instr); keyword_recognized = true;}

  if (keystr == "NONSYMM_COEFS")
    { readNonsymmCoefficients(instr); keyword_recognized = true;}

  if (keystr == "SHARED_IDS")
    { readSharedIDs(instr); keyword_recognized = true;}

  if (!keyword_recognized) {
    std::string msg(": keyword not recognized");
    throw std::runtime_error(keyword+msg);
  }
}

//--------------------------------------------------------------------
void InputFileReader::readField(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::VectorSpace> rowspace = getRowSpace();

  int fieldID = 0;
  int fieldSize = -1;
  DataReader::readData(instr, fieldID);
  DataReader::readData(instr, fieldSize);

  rowspace->defineFields(1, &fieldID, &fieldSize);
}

//--------------------------------------------------------------------
void InputFileReader::readIDType(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::VectorSpace> rowspace = getRowSpace();

  int idType = 0;
  DataReader::readData(instr, idType);

  rowspace->defineIDTypes(1, &idType);
}

//--------------------------------------------------------------------
void InputFileReader::readSimplePattern(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();

  int patternID = 0, numIDs = 0, idType = 0, fieldID = 0;
  DataReader::readData(instr, patternID);
  DataReader::readData(instr, numIDs);
  DataReader::readData(instr, idType);
  DataReader::readData(instr, fieldID);

  mgraph->definePattern(patternID, numIDs, idType, fieldID);
}

//--------------------------------------------------------------------
void InputFileReader::readSymmElemBlock(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();

  int blockID = 0, numElems = 0, patternID = 0;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, numElems);
  DataReader::readData(instr, patternID);

  int errcode = mgraph->initConnectivityBlock(blockID, numElems, patternID);
  if (errcode != 0) {
    throw std::runtime_error("err in mgraph->initConnectivityBlock");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readSharedIDs(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::VectorSpace> vspace = getRowSpace();

  int numSharedIDs = 0, idType = 0;
  DataReader::readData(instr, idType);
  DataReader::readData(instr, numSharedIDs);

  std::vector<int> ids(numSharedIDs);
  std::vector<int> numSharingProcs(numSharedIDs);
  std::vector<int> sharingProcs;

  int temp = 0;
  for(int i=0; i<numSharedIDs; ++i) {
    DataReader::readData(instr, temp);
    ids[i] = temp;
    DataReader::readData(instr, temp);
    numSharingProcs[i] = temp;
    for(int j=0; j<numSharingProcs[i]; ++j) {
      DataReader::readData(instr, temp);
      sharingProcs.push_back(temp);
    }
  }

  int errcode = vspace->initSharedIDs(numSharedIDs, idType,
					&ids[0], &numSharingProcs[0],
					&sharingProcs[0]);
  if (errcode != 0) {
    throw std::runtime_error("err in vspace->initSharedIDs");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readNonsymmConnBlock(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();

  int blockID = 0, numContributions = 0, rowPatternID = 0, colPatternID;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, numContributions);
  DataReader::readData(instr, rowPatternID);
  DataReader::readData(instr, colPatternID);

  int errcode = mgraph->initConnectivityBlock(blockID, numContributions,
					      rowPatternID, colPatternID);
  if (errcode != 0) {
    throw std::runtime_error("err in mgraph->initConnectivityBlock");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readNonsymmConnectivity(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();

  int blockID = 0, connID = 0;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, connID);

  fei::ConnectivityBlock* cblock = mgraph->getConnectivityBlock(blockID);
  if (cblock == NULL) {
    throw std::runtime_error("error getting connectivity block");
  }

  int numRowIDs = cblock->getRowPattern()->getNumIDs();
  int numColIDs = cblock->getColPattern()->getNumIDs();

  std::vector<int> rowIDs(numRowIDs);
  for(int i=0; i<numRowIDs; ++i) {
    int id = 0;
    DataReader::readData(instr, id);
    rowIDs[i] = id;
  }

  std::vector<int> colIDs(numColIDs);
  for(int i=0; i<numColIDs; ++i) {
    int id = 0;
    DataReader::readData(instr, id);
    colIDs[i] = id;
  }

  int errcode = mgraph->initConnectivity(blockID, connID,
					 &rowIDs[0], &colIDs[0]);
  if (errcode != 0) {
    throw std::runtime_error("err in mgraph->initConnectivityBlock");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readNonsymmCoefficients(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  fei::SharedPtr<fei::Matrix> matrix = getMatrix();

  int blockID = 0, connID = 0;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, connID);

  fei::ConnectivityBlock* cblock = mgraph->getConnectivityBlock(blockID);
  if (cblock == NULL) {
    throw std::runtime_error("error getting connectivity block");
  }

  int numRowIDs = cblock->getRowPattern()->getNumIDs();
  int numColIDs = cblock->getColPattern()->getNumIDs();

  std::vector<double> coefs;
  for(int i=0; i<numRowIDs; ++i) {
    for(int j=0; j<numColIDs; ++j) {
      double coef = 0.0;
      DataReader::readData(instr, coef);
      coefs.push_back(coef);
    }
  }

  std::vector<double*> coefs_2d(numRowIDs);
  for(int i=0; i<numRowIDs; ++i) {
    coefs_2d[i] = &coefs[i*numColIDs];
  }

  int errcode = matrix->sumIn(blockID, connID,
			      &coefs_2d[0]);
  if (errcode != 0) {
    throw std::runtime_error("err in matrix->sumIn(non-symm)");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readSymmConnectivity(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();

  int blockID = 0, elemID = 0, numNodes = 0;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, elemID);
  DataReader::readData(instr, numNodes);

  std::vector<int> nodeIDs(numNodes);

  for(int i=0; i<numNodes; ++i) {
    DataReader::readData(instr, nodeIDs[i]);
  }

  int errcode = mgraph->initConnectivity(blockID, elemID, &nodeIDs[0]);
  if (errcode != 0) {
    throw std::runtime_error("err in mgraph->initConnectivity");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readSlaveConstraint(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  fei::SharedPtr<fei::VectorSpace> vspace = mgraph->getRowSpace();

  int slaveID = 0, slaveIDType = 0, slaveField = 0, slaveFieldOffset = 0;
  int numMasters = 0;

  DataReader::readData(instr, slaveID);
  DataReader::readData(instr, slaveIDType);
  DataReader::readData(instr, slaveField);
  DataReader::readData(instr, slaveFieldOffset);
  DataReader::readData(instr, numMasters);

  std::vector<int> IDs(numMasters+1);
  std::vector<int> IDTypes(numMasters+1);
  std::vector<int> Fields(numMasters+1);
  std::vector<double> weights(vspace->getFieldSize(slaveField), 0.0);

  IDs[0] = slaveID;
  IDTypes[0] = slaveIDType;
  Fields[0] = slaveField;
  weights[slaveFieldOffset] = -1.0;

  for(int i=1; i<numMasters+1; ++i) {
    DataReader::readData(instr, IDs[i]);
    DataReader::readData(instr, IDTypes[i]);
    DataReader::readData(instr, Fields[i]);
    int fieldSize = 0;
    DataReader::readData(instr, fieldSize);
    for(int j=0; j<fieldSize; ++j) {
      double weight = 0.0;
      DataReader::readData(instr, weight);
      weights.push_back(weight);
    }
  }

  double rhsValue = 0.0;
  DataReader::readData(instr, rhsValue);

  int errcode = mgraph->initSlaveConstraint(numMasters+1, &IDTypes[0],
					    &IDs[0], &Fields[0],
					    0, slaveFieldOffset,
					    &weights[0], rhsValue);
  if (errcode != 0) {
    throw std::runtime_error("err in mgraph->initSlaveConstraint");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readLagrangeConstraint(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  fei::SharedPtr<fei::VectorSpace> vspace = mgraph->getRowSpace();

  int numIDs = 0, crIDType = 0;

  DataReader::readData(instr, crIDType);
  DataReader::readData(instr, numIDs);

  std::vector<int> idtypes(numIDs);
  std::vector<int> ids(numIDs);
  std::vector<int> fields(numIDs);
  std::vector<double> weights;

  for(int i=0; i<numIDs; ++i) {
    DataReader::readData(instr, ids[i]);
    DataReader::readData(instr, idtypes[i]);
    DataReader::readData(instr, fields[i]);
    int fieldSize = vspace->getFieldSize(fields[i]);

    for(int j=0; j<fieldSize; ++j) {
      double weight = 0.0;
      DataReader::readData(instr, weight);
      weights.push_back(weight);
    }
  }

  double rhsValue = 0.0;
  DataReader::readData(instr, rhsValue);

  std::map<int,ConstraintType*>::iterator iter = lagrangeConstraints_.find(crID_);
  if (iter != lagrangeConstraints_.end()) {
    throw std::runtime_error("error, duplicate constraint");
  }

  ConstraintType* newconstraint = 
    new ConstraintType(crID_, crIDType,
		       false,//isSlave
		       false,//isPenalty
		       numIDs, &idtypes[0], &ids[0], &fields[0],
		       0, 0,//offsetOfSlave & offsetIntoSlaveField
		       &weights[0], rhsValue, vspace.get());

  lagrangeConstraints_[crID_] = newconstraint;

  int errcode = mgraph->initLagrangeConstraint(crID_, crIDType, numIDs,
					       &idtypes[0], &ids[0], &fields[0]);
  ++crID_;

  if (errcode != 0) {
    throw std::runtime_error("err in mgraph->initLagrangeConstraint");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readElemStiffness(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  fei::SharedPtr<fei::Matrix> matrix = getMatrix();

  int blockID = 0, elemID = 0, numIndices = 0;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, elemID);

  numIndices = mgraph->getConnectivityNumIndices(blockID);

  std::vector<double> coefs;
  std::vector<double*> coefs_2d;

  for(int i=0; i<numIndices; ++i) {
    for(int j=0; j<numIndices; ++j) {
      double coef;
      DataReader::readData(instr, coef);
      coefs.push_back(coef);
    }
  }
  for(int i=0; i<numIndices; ++i) {
    coefs_2d.push_back(&coefs[i*numIndices]);
  }

  int errcode = matrix->sumIn(blockID, elemID, &coefs_2d[0]);
  if (errcode != 0) {
    throw std::runtime_error("err in matrix->sumIn");
  }
}

//--------------------------------------------------------------------
void InputFileReader::readElemLoad(FEI_ISTREAM* instr)
{
  fei::SharedPtr<fei::MatrixGraph> mgraph = getMatrixGraph();
  fei::SharedPtr<fei::Vector> rhs = getVector(false);

  int blockID = 0, elemID = 0, numIndices = 0;
  DataReader::readData(instr, blockID);
  DataReader::readData(instr, elemID);

  numIndices = mgraph->getConnectivityNumIndices(blockID);

  std::vector<double> coefs;

  for(int i=0; i<numIndices; ++i) {
    double coef;
    DataReader::readData(instr, coef);
    coefs.push_back(coef);
  }

  std::vector<int> indices(numIndices);
  mgraph->getConnectivityIndices(blockID, elemID, numIndices, &indices[0],
                                 numIndices);
  int errcode = rhs->sumIn(numIndices, &indices[0], &coefs[0], 0);
  if (errcode != 0) {
    throw std::runtime_error("err in rhs->sumIn");
  }
}

//--------------------------------------------------------------------
int InputFileReader::testInitialization()
{
  return(0);
}

//--------------------------------------------------------------------
int InputFileReader::testLoading()
{
  return(0);
}

//--------------------------------------------------------------------
void InputFileReader::dumpMatrixFiles()
{
}

//--------------------------------------------------------------------
void InputFileReader::setParameter(const char*)
{
}

//--------------------------------------------------------------------
void InputFileReader::
loadLagrangeConstraints(fei::SharedPtr<fei::LinearSystem> linsys)
{
  std::map<int,ConstraintType*>::iterator
    iter = lagrangeConstraints_.begin(),
    iter_end = lagrangeConstraints_.end();

  for(; iter != iter_end; ++iter) {
    ConstraintType* cr = (*iter).second;
    int crID = cr->getConstraintID();
    std::vector<double>& weights = cr->getMasterWeights();
    double rhsValue = cr->getRHSValue();

    int errcode = linsys->loadLagrangeConstraint(crID, &weights[0], rhsValue);
    if (errcode != 0) {
      throw std::runtime_error("error in linsys->loadLagrangeConstraint");
    }
  }
}

//--------------------------------------------------------------------
int InputFileReader::testSolve()
{
  fei::SharedPtr<fei::Solver> solver = factory_->createSolver();

  fei::SharedPtr<fei::LinearSystem> linsys = getLinearSystem();

  loadLagrangeConstraints(linsys);

  linsys->loadComplete();

  int status, itersTaken = 0;
  int errcode = solver->solve(linsys.get(), NULL,
			      *feiParameterSet, itersTaken, status);
  errcode += linsys->getSolutionVector()->scatterToOverlap();
  solveCounter_++;
  return(errcode);
}

//--------------------------------------------------------------------
int InputFileReader::testCheckResult()
{
  FEI_COUT << "testCheckResult"<<FEI_ENDL;
  int numProcs = 1, localProc = 0;
#ifndef FEI_SER
  MPI_Comm_size(comm_, &numProcs);
  MPI_Comm_rank(comm_, &localProc);
#endif

  FEI_OSTRINGSTREAM osstr;
  osstr << solnFileName_ << ".idtype0." << solveCounter_<<"."
	<<numProcs << "." << localProc;
  std::string fullSolnFileName = osstr.str();
  const char* solnFile_c_str = fullSolnFileName.c_str();
  save_soln(0, solnFile_c_str);

  int err = SolnCheck::checkSolution(localProc, numProcs,
				     solnFileName_.c_str(),
				     checkFileName_.c_str(),
				     "idtype0", solveCounter_);

  int globalErr = err;
#ifndef FEI_SER
  if (MPI_SUCCESS != MPI_Allreduce(&err, &globalErr, 1, MPI_INT, MPI_SUM,
				   comm_)) return(-1);
#endif
  if (globalErr != 0) return(-1);

  return(0);
}

//--------------------------------------------------------------------
void InputFileReader::save_soln(int idType, const char* fileName)
{
  FEI_COUT << "save_soln opening '"<<fileName<<"'"<<FEI_ENDL;
  FEI_OFSTREAM outfile(fileName);
  if (!outfile || outfile.bad()) {
    throw std::runtime_error("InputFileReader::save_soln ERROR opening outfile.");
  }

  fei::SharedPtr<fei::VectorSpace> vspace = getRowSpace();
  int numIDs = vspace->getNumOwnedAndSharedIDs(idType);
  std::vector<int> ids(numIDs);
  int errcode = vspace->getOwnedAndSharedIDs(idType, numIDs, &ids[0], numIDs);
  if (errcode != 0) {
    throw std::runtime_error("error in getLocalIDs");
  }

  fei::SharedPtr<fei::LinearSystem> linsys = getLinearSystem();
  fei::SharedPtr<fei::Vector> x = linsys->getSolutionVector();

  std::vector<double> coefs;
  std::vector<int> fields;
  int numFields;

  for(int i=0; i<numIDs; ++i) {
    int numDOF = vspace->getNumDegreesOfFreedom(idType, ids[i]);
    outfile << ids[i] << " " << numDOF << FEI_ENDL;

    numFields = vspace->getNumFields(idType, ids[i]);
    fields.resize(numFields);
    vspace->getFields(idType, ids[i], fields);

    for(int j=0; j<numFields; ++j) {
      int fieldSize = vspace->getFieldSize(fields[j]);
      coefs.resize(fieldSize);
      errcode = x->copyOutFieldData(fields[j], idType, 1, &ids[i], &coefs[0]);
      if (errcode != 0) {
	throw std::runtime_error("error in copyOutFieldData");
      }

      for(int k=0; k<fieldSize; ++k) {
	outfile << coefs[k] << " ";
      }
    }
    outfile << FEI_ENDL;
  }
}

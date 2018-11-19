
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef BDDC_SETUPTEST_H
#define BDDC_SETUPTEST_H
#include <gtest/gtest.h>

#include <mpi.h>
#include "ProblemMakerBDDC.hpp"
#include "OperatorStandard.hpp"
#include "ShyLU_DDBDDC_config.h"
#include "shylu_PreconditionerBDDC.hpp"
#include "shylu_KrylovSolver.hpp"
#include "shylu_UtilBDDC.hpp"
#include "shylu_enumsBDDC.hpp"

#if defined(_OPENMP)
#include <omp.h>
#else
#define omp_get_nested() 0
#endif

using Teuchos::RCP;

namespace bddc {

template <class LO, class GO, class SX, class SM> class setupTest
{
 public:

setupTest(RCP<Teuchos::ParameterList> & parametersPM,
	  RCP<Teuchos::ParameterList> & parametersBDDC,
	  RCP<Teuchos::ParameterList> & parametersMueLu,
	  RCP<Teuchos::ParameterList> & parametersNodalAMG,
	  const std::string & meshDataFile,
	  MPI_Comm Comm)
{
  int numProc, myPID;
  MPI_Comm_size(Comm, &numProc);
  MPI_Comm_rank(Comm, &myPID);
  if (myPID == 0) {
    std::cout << "sizeof(LO) = " << sizeof(LO) << std::endl;
    std::cout << "sizeof(GO) = " << sizeof(GO) << std::endl;
  }
  int spatialDim, problemTypeInt;
  std::vector<SM> xFile, yFile, zFile;
  std::vector<GO> nodeGIDs;
  std::vector< std::vector<LO> > elemConn, subElems;
  std::vector<LO> constrainedNodes;
  // read in problem specification from separate file if requested
  readProblem(meshDataFile, myPID, spatialDim, problemTypeInt, 
	      xFile, yFile, zFile, elemConn, subElems, nodeGIDs, 
	      constrainedNodes);
  if (meshDataFile == "") {
    int returnValue = 
      setParameters1(numProc, parametersPM, parametersBDDC, parametersMueLu,
		     parametersNodalAMG);
    if (returnValue != 0) return;
    setParameters2(myPID, parametersPM, parametersBDDC, parametersMueLu,
		   parametersNodalAMG);
    m_Problem = rcp( new bddc::ProblemMaker<LO,GO,SX,SM>(parametersPM, Comm) );
  }
  else {
    setParameters1(numProc, parametersPM, parametersBDDC, parametersMueLu,
		   parametersNodalAMG);
    setParameters2(myPID, parametersPM, parametersBDDC, parametersMueLu,
		   parametersNodalAMG);
    m_Problem = rcp( new bddc::ProblemMaker<LO,GO,SX,SM>
		     (parametersPM, Comm, xFile, yFile, zFile, elemConn, 
		      nodeGIDs, subElems, constrainedNodes) );
  }
  m_Problem->getSubDomainNodeData(m_subNodes, m_subNodeBegin, m_subLocalDofs);

  m_Problem->getSubdomainMatrices(m_subNodes, m_subNodeBegin, m_subRowBegin, 
				  m_subColumns, m_subValues);
  const bool runTPtest(false);
  if (runTPtest) {
  // following is for higher-order elements
    m_Problem->adjustNodalCoords(m_subNodes);
#if defined(TENSORPRODUCTBDDC)
    const int numNodeSub = m_subNodes[0].size();
    std::vector<double> xTP(numNodeSub), yTP(numNodeSub), zTP(numNodeSub);
    for (int i=0; i<numNodeSub; i++) {
      const int node = m_subNodes[0][i];
      xTP[i] = xCoord[node];
      yTP[i] = yCoord[node];
      zTP[i] = zCoord[node];
    }
    bddc::TensorProduct* TP = m_Problem->getTP();
    std::vector<double> KTP, MTP;
    const double mu(1), rho(1);
    spatialDim = parametersPM->get("Spatial Dimension", 0);
    int integrationOption = 0;
    TP->calculateMatrices(xTP.data(), yTP.data(), zTP.data(), mu, rho, 
			  spatialDim, integrationOption, KTP, MTP);
    bddc::UtilBDDC<double,double>::printDenseMatrix(numNodeSub, numNodeSub,
						    KTP.data(), "KGLeg.dat");
    bddc::UtilBDDC<double,double>::printDenseMatrix(numNodeSub, numNodeSub,
						    MTP.data(), "MGLeg.dat");
    integrationOption = 1;
    TP->calculateMatrices(xTP.data(), yTP.data(), zTP.data(), mu, rho, 
			  spatialDim, integrationOption, KTP, MTP);
    bddc::UtilBDDC<double,double>::printDenseMatrix(numNodeSub, numNodeSub,
						    KTP.data(), "KGLob.dat");
    bddc::UtilBDDC<double,double>::printDenseMatrix(numNodeSub, numNodeSub,
						    MTP.data(), "MGLob.dat");
    const LO numRowsSub = subNodeBegin[0][numNodeSub];
    bddc::UtilBDDC<double,double>::printSparseMatrix
      (numRowsSub, m_subRowBegin[0].data(), m_subColumns[0].data(),
       m_subValues[0].data(), "Ksub.dat");
#endif
  }
  // use operator if specified
  const std::string operatorName = 
    parametersBDDC->get("Operator Name", "none");
  const int spatialDim2 = parametersPM->get("Spatial Dimension", 0);
  const int quadTypeTP = parametersPM->get("Quad Type TP", 1);
  const SM* xCoord = getXcoord();
  const SM* yCoord = getYcoord();
  const SM* zCoord = getZcoord();
  initializeOperator(m_Problem, operatorName, m_Operator, m_subNodes, 
		     m_subRowBegin, m_subColumns, m_subValues, spatialDim2, 
		     quadTypeTP, xCoord, yCoord, zCoord);
  if (m_Operator != nullptr) {
    parametersBDDC->set("Avoid Initial Static Condensation", true);
  }
}
    
~setupTest()
{
  delete m_Operator;
}

double getTime()
{
  return MPI_Wtime();
}

OperatorBase<SX>* getOperator()
{
  return m_Operator;
}

void openFiles
  (RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > & Preconditioner,
   const bool resetFile)
{
  Preconditioner->printBanner("dd_solver.dat", resetFile);
  RCP<const Teuchos::Comm<int> > Comm = Preconditioner->getComm();
  LO myPID = Comm->getRank();
  if (myPID == 0) {
    std::ofstream outputFileKrylov;
    if (resetFile) {
      outputFileKrylov.open("krylov_solver.dat", std::ios::out);
    }
    else {
      outputFileKrylov.open("krylov_solver.dat", 
			    std::ios::out | std::ios::app);
    }
  }
}

void printTimings
(RCP< bddc::PreconditionerBDDC<SX,SM,LO,GO> > & Preconditioner, 
 RCP< bddc::KrylovSolver<SX,SM,LO,GO> > & Solver,
 const LO krylovMethod)
{
  RCP<const Teuchos::Comm<int> > Comm = Preconditioner->getComm();
  LO myPID = Comm->getRank(); 
  if (myPID != 0) return;
  std::ofstream timingsFile;
  timingsFile.open("KrylovSolver_timers.dat", std::ios::out);
  if (krylovMethod == 0) {
    timingsFile << "OthogGCR:projections         = "
		<< Solver->getProjectionTime() << std::endl;
    timingsFile << "OthogGCR:orthogonalizations  = "
		<< Solver->getOrthogonalizationTime() << std::endl;
  }
  timingsFile.close();
}

void getRhs(const LO numMyRows,
	    std::vector<SX> & rhs)
{
  rhs.resize(numMyRows);
  srand(7);
  for (LO i=0; i<numMyRows; i++) {
    rhs[i] = 0.7*rand()/RAND_MAX;
  }
}

void updateRhs(SX* rhs, 
	       LO numRows)
{
  for (LO i=0; i<numRows; i++) {
    rhs[i] += 0.07*rand()/RAND_MAX;
  }
}

void setSubPointers(std::vector<LO*> & subRowBeginPtr, 
		    std::vector<LO*> & subColumnsPtr, 
		    std::vector<SX*> & subValuesPtr)
{
  const int numSub = m_subNodes.size();
  subRowBeginPtr.resize(numSub);
  subColumnsPtr.resize(numSub);
  subValuesPtr.resize(numSub);
  for (int i=0; i<numSub; i++) {
    subRowBeginPtr[i] = m_subRowBegin[i].data();
    subColumnsPtr[i] = m_subColumns[i].data();
    subValuesPtr[i] = m_subValues[i].data();
  }
}

const GO* getNodeGlobalIDs()
{
  return m_Problem->getNodeGlobalIDs();
}

const LO* getNodeBegin() const
{
  return m_Problem->getNodeBegin();
}

const LO* getLocalDofs() const
{
  return m_Problem->getLocalDofs();
}

const LO getNumNode() const
{
  return m_Problem->getNumNode();
}

const SM* getXcoord() const
{
  return m_Problem->getXcoord();
}

 const SM* getYcoord() const
 {
   return m_Problem->getYcoord();
 }

 const SM* getZcoord() const
 {
   return m_Problem->getZcoord();
 }

const std::vector< std::vector<LO> > & getSubNodes() const
{
  return m_subNodes;
}

void getProblemData(LO & numNode, 
		    LO* & nodeBegin, 
		    LO* & localDofs, 
		    const GO* & nodeGlobalIDs,
		    const SM* & xCoord, 
		    const SM* & yCoord, 
		    const SM* & zCoord, 
		    std::vector<LO*> & subRowBeginPtr,
		    std::vector<LO*> & subColumnsPtr, 
		    std::vector<SX*> & subValuesPtr,
		    OperatorBase<SX>* & Operator)
{
  numNode = getNumNode();
  nodeBegin = const_cast<LO*>(getNodeBegin());
  localDofs = const_cast<LO*>(getLocalDofs());
  nodeGlobalIDs = getNodeGlobalIDs();
  xCoord = getXcoord();
  yCoord = getYcoord();
  zCoord = getZcoord();
  setSubPointers(subRowBeginPtr, subColumnsPtr, subValuesPtr);
  Operator = getOperator();
}

 private:
 std::vector< std::vector<LO> > m_subNodes, m_subNodeBegin, m_subRowBegin, 
    m_subLocalDofs, m_subColumns;
 std::vector< std::vector<SX> > m_subValues;
 RCP< ProblemMaker<LO,GO,SX,SM> > m_Problem;
 OperatorBase<SX>* m_Operator{nullptr};

void initializeOperator(RCP< bddc::ProblemMaker<LO,GO,SX,SM> > & Problem,
			const std::string & operatorName, 
			bddc::OperatorBase<SX>* & Operator,
			std::vector< std::vector<LO> > & subNodes,
			const std::vector< std::vector<int> > & subRowBegin,
			const std::vector< std::vector<int> > & subColumns,
			const std::vector< std::vector<SX> > & subValues,
			const int spatialDim,
			const int quadTypeTP,
			const SM* xCoord, 
			const SM* yCoord,
			const SM* zCoord)
{
  Operator = nullptr;
  LO numNode = Problem->getNumNode();
  LO* nodeBegin = const_cast<LO*>(Problem->getNodeBegin());
  if (operatorName == "Standard") {
    Operator = new bddc::OperatorStandard<SX>
      (numNode, nodeBegin, subNodes, subRowBegin, subColumns, subValues);
   
  }
  else if (operatorName == "TP") {
#ifdef TENSORPRODUCTBDDC
    LO numDofNode = 1; // change for elasticity
    LO* localDofs = const_cast<LO*>(Problem->getLocalDofs());
    // just get matrix for first element
    const int numNodeSub = subNodes[0].size();
    std::vector<double> xTP(numNodeSub), yTP(numNodeSub), zTP(numNodeSub);
    for (int i=0; i<numNodeSub; i++) {
      const int node = subNodes[0][i];
      xTP[i] = xCoord[node];
      yTP[i] = yCoord[node];
      zTP[i] = zCoord[node];
    }
    bddc::TensorProduct* TP = Problem->getTP();
    std::vector<double> KTP, MTP;
    const double mu(1), rho(1);
    int integrationOption = quadTypeTP;
    TP->calculateMatrices(xTP.data(), yTP.data(), zTP.data(), mu, rho, 
			  spatialDim, integrationOption, KTP, MTP);
    Operator = new bddc::OperatorTP<SX>
      (numNode, numDofNode, nodeBegin, localDofs, subNodes, KTP);
#else
    std::string msg("Error: Tensor product BDDC not available");
    throw msg;
#endif
  }
}

void readProblem(const std::string & meshDataFile, 
		 const int myPID,
		 LO & spatialDim,
		 LO & problemType,
		 std::vector<SM> & xFile, 
		 std::vector<SM> & yFile,
		 std::vector<SM> & zFile, 
		 std::vector< std::vector<LO> > & elemConn, 
		 std::vector< std::vector<LO> > & subElems, 
		 std::vector<GO> & nodeGIDs,
		 std::vector<LO> & constrainedNodes)
{
  if (meshDataFile == "") return;
  std::ifstream fin;
  fin.open(meshDataFile.c_str());
  if ( ! fin.is_open()) {
    std::cout << "Error: could not open the file " << meshDataFile << std::endl;
    return;
  }
  LO numNode, numElem, numNodePerElem;
  char buff[101];
  fin >> spatialDim; fin.getline(buff,100);
  fin >> problemType; fin.getline(buff,100);
  fin >> numNode; fin.getline(buff,100);
  fin >> numElem; fin.getline(buff,100);
  fin >> numNodePerElem; fin.getline(buff,100);
  std::vector<LO>  myElems, mySubs, elemProc(numElem), elemLocalSub(numElem);
  LO eProc, localSub, maxLocalSub(0), numMyElem(0);
  fin.getline(buff,100);
  for (LO i=0; i<numElem; i++) {
    fin >> eProc >> localSub; fin.getline(buff,100);
    elemProc[i] = eProc;
    elemLocalSub[i] = localSub;
    if (elemProc[i] == myPID) {
      if (localSub > maxLocalSub) maxLocalSub = localSub;
      numMyElem++;
    }
  }
  elemConn.resize(numMyElem);
  subElems.resize(maxLocalSub+1);
  std::vector<LO> count(maxLocalSub+1, 0);
  for (LO i=0; i<numElem; i++) {
    if (elemProc[i] == myPID) {
      count[elemLocalSub[i]]++;
    }
  }
  for (LO i=0; i<maxLocalSub+1; i++) {
    subElems[i].resize(count[i]);
    count[i] = 0;
  }
  fin.getline(buff,100);
  numMyElem = 0;
  for (LO i=0; i<numElem; i++) {
    if (elemProc[i] == myPID) {
      const LO localSub = elemLocalSub[i];
      subElems[localSub][count[localSub]] = numMyElem;
      std::vector<LO> & eConn = elemConn[numMyElem];
      eConn.resize(numNodePerElem);
      for (LO j=0; j<numNodePerElem; j++) {
	fin >> eConn[j];
      }
      count[localSub]++;
      numMyElem++;
    }
    fin.getline(buff,100);
  }
  nodeGIDs.resize(numNodePerElem*numMyElem);
  LO numTerms(0);
  for (LO i=0; i<numMyElem; i++) {
    for (LO j=0; j<numNodePerElem; j++) {
      nodeGIDs[numTerms++] = elemConn[i][j];
    }
  }
  bddc::DofManager<LO,GO>::determineUniqueIndices(nodeGIDs);
  std::map<GO, LO> nodeMap;
  const LO numMyNode = nodeGIDs.size();
  for (LO i=0; i<numMyNode; i++) {
    nodeMap[nodeGIDs[i]] = i;
  }
  for (LO i=0; i<numMyElem; i++) {
    for (LO j=0; j<numNodePerElem; j++) {
      const GO node = elemConn[i][j];
      auto iter = nodeMap.find(node);
      BDDC_TEST_FOR_EXCEPTION(iter == nodeMap.end(), std::runtime_error, 
			      "node not found");
      elemConn[i][j] = iter->second;
    }
  }
  SM x, y, z;
  GO globalID;
  LO nodeFlag;
  xFile.resize(numMyNode);
  yFile.resize(numMyNode);
  zFile.resize(numMyNode);
  fin.getline(buff,100);
  for (LO i=0; i<numNode; i++) {
    fin >> globalID >> x >> y >> z >> nodeFlag; fin.getline(buff,100);
    auto iter = nodeMap.find(globalID);
    if (iter != nodeMap.end()) {
      const LO myNode = iter->second;
      xFile[myNode] = x;
      yFile[myNode] = y;
      zFile[myNode] = z;
      if (nodeFlag == 1) {
	constrainedNodes.push_back(myNode);
      }
    }
  }
}

int setParameters1(const int numProc,
		   RCP<Teuchos::ParameterList> & parametersPM,
		   RCP<Teuchos::ParameterList> & parametersBDDC,
		   RCP<Teuchos::ParameterList> & parametersMueLu,
		   RCP<Teuchos::ParameterList> & parametersNodalAMG)
{
  // set parameters
  LO numSubregionsDir1 = parametersPM->get("Number Subregions Direction 1", 1);
  LO numSubregionsDir2 = parametersPM->get("Number Subregions Direction 2", 1);
  LO numSubregionsDir3 = parametersPM->get("Number Subregions Direction 3", 1);
  LO numSubdomainsPerSubregionDir1 = parametersPM->get
    ("Number Subdomains Per Subregion Direction 1", 1);
  LO numSubdomainsPerSubregionDir2 = parametersPM->get
    ("Number Subdomains Per Subregion Direction 2", 1);
  LO numSubdomainsPerSubregionDir3 = parametersPM->get
    ("Number Subdomains Per Subregion Direction 3", 1);
  LO Hh1 = parametersPM->get("Number of Elements Per Subdomain Direction 1", 4);
  LO Hh2 = parametersPM->get("Number of Elements Per Subdomain Direction 2", 4);
  LO Hh3 = parametersPM->get("Number of Elements Per Subdomain Direction 3", 4);
  parametersPM->set("Num Elems Per Sub Dir 1", Hh1);
  parametersPM->set("Num Elems Per Sub Dir 2", Hh2);
  parametersPM->set("Num Elems Per Sub Dir 3", Hh3);
  // problem type and spatial dimension
  const std::string problemTypeString = 
    parametersPM->get("Problem Type String", "Poisson-3D");
  if (parametersMueLu != Teuchos::null) {
    const std::string problemTypeMueLu = 
      parametersMueLu->get("problem: type", "abc");
    if (usingMueLu(parametersBDDC)) {
      BDDC_TEST_FOR_EXCEPTION(problemTypeMueLu != problemTypeString, 
	  std::runtime_error, "inconsistent problem type in MueLu .xml file");
    }
    parametersMueLu->set("problem: type", problemTypeString);
  }
  if (parametersNodalAMG != Teuchos::null) {
    parametersNodalAMG->set("Problem Type", problemTypeString);
  }
  enum bddc::ProblemType problemType = bddc::SCALARPDE;
  if ((problemTypeString == "Elasticity-2D") ||
      (problemTypeString == "Elasticity-3D")) problemType = bddc::ELASTICITY;
  int spatialDim = 3;
  if ((problemTypeString == "Poisson-2D") || 
      (problemTypeString == "Elasticity-2D")) spatialDim = 2;  
  parametersPM->set("Problem Type", problemType);
  parametersPM->set("Spatial Dimension", spatialDim);
  enum bddc::AnalysisType analysisType = bddc::STANDARD;
  LO numElemPerSubDir1 = numSubdomainsPerSubregionDir1*Hh1;
  LO numElemPerSubDir2 = numSubdomainsPerSubregionDir2*Hh2;
  LO numElemPerSubDir3 = numSubdomainsPerSubregionDir3*Hh3;
  if (spatialDim == 2) {
    numSubregionsDir3 = 1;
    numElemPerSubDir3 = 1;
  }
  parametersPM->set("Number of Subdomains Direction 1", numSubregionsDir1);
  parametersPM->set("Number of Subdomains Direction 2", numSubregionsDir2);
  parametersPM->set("Number of Subdomains Direction 3", numSubregionsDir3);
  parametersPM->set("Number of Elements Per Subregion Direction 1",
		    numElemPerSubDir1);
  parametersPM->set("Number of Elements Per Subregion Direction 2",
		    numElemPerSubDir2);
  parametersPM->set("Number of Elements Per Subregion Direction 3",
		    numElemPerSubDir3);
  int numDofPerNode = 1;
  if (problemType == bddc::ELASTICITY) {
    numDofPerNode = spatialDim;
  }
  int numExpectedProc = numSubregionsDir1*numSubregionsDir2*numSubregionsDir3;
  if (numProc < numExpectedProc) {
    std::cout << "Error: number of processors must be at least " 
	      << numExpectedProc << std::endl;
    return 1;
  }
  //  
  parametersPM->set("numDofPerNode", numDofPerNode);
  parametersPM->set("Problem Type", problemType);
  parametersBDDC->set("Problem Type BDDC", problemType);
  parametersBDDC->set("nullSpaceCorrection", false);
  if (parametersNodalAMG != Teuchos::null) {
    if (parametersNodalAMG->get("nullSpaceCorrection", false)) {
      parametersBDDC->set("nullSpaceCorrection", true);
    }
  }
  else {
    if (parametersMueLu != Teuchos::null) {
      if (usingMueLu(parametersBDDC)) {
	parametersBDDC->set("nullSpaceCorrection", true);
      }
    }
  }
  parametersPM->set("Analysis Type", analysisType);
  bool addAsymmetry = parametersPM->get("Add Asymmetry", false);
  if (addAsymmetry == false) {
    parametersBDDC->set("Matrix Type", "Symmetric");
    parametersBDDC->set("Pardiso Matrix Type", 
			"Real And Symmetric Positive Definite");
  }
  else {
    parametersBDDC->set("Matrix Type", "NonSymmetric");
    parametersBDDC->set("Pardiso Matrix Type", "Real And Structurally Symmetric");
  }
  return 0;
}

void setParameters2(const int myPID,
		    RCP<Teuchos::ParameterList> & parametersPM,
		    RCP<Teuchos::ParameterList> & parametersBDDC,
		    RCP<Teuchos::ParameterList> & parametersMueLu,
		    RCP<Teuchos::ParameterList> & parametersNodalAMG)
{
  checkInterfacePreconditionerSolverCompatibility(myPID, *parametersBDDC);

  parametersBDDC->set("Problem Type BDDC", 
		      parametersPM->get("Problem Type", bddc::SCALARPDE));
  parametersBDDC->set("Spatial Dimension", 
		      parametersPM->get("Spatial Dimension", 3));

  std::string weightCornerString = parametersBDDC->get("Corner Weight Type", "stiffness");
  std::string weightEdgeString = parametersBDDC->get("Edge Weight Type", "stiffness");
  std::string weightFaceString = parametersBDDC->get("Face Weight Type", "stiffness");
  enum bddc::WeightType weightTypeCorner = bddc::STIFFNESS;
  enum bddc::WeightType weightTypeEdge = bddc::STIFFNESS;
  enum bddc::WeightType weightTypeFace = bddc::STIFFNESS;
  // no deluxe scaling for now
  if (weightCornerString == "cardinality") {
    weightTypeCorner = bddc::CARDINALITY;
  }
  if (weightEdgeString == "cardinality") {
    weightTypeEdge = bddc::CARDINALITY;
  }
  if (weightFaceString == "cardinality") {
    weightTypeFace = bddc::CARDINALITY;
  }
  parametersBDDC->set("Weight Type Corner", weightTypeCorner);
  parametersBDDC->set("Weight Type Edge", weightTypeEdge);
  parametersBDDC->set("Weight Type Face", weightTypeFace);
  std::string coarseningOption = parametersBDDC->get("Coarsening Option", "Graph");
  if (coarseningOption == "Graph") {
    parametersBDDC->set("Construct Subdomain Adjacency Graph", true);
  }
  else {
    parametersBDDC->set("Construct Subdomain Adjacency Graph", false);
  }
  parametersBDDC->set("Krylov Method", 0); // GCR
  const std::string & krylovMethod = parametersBDDC->get("Krylov Solver", "GCR");
  if (krylovMethod == "PCG") {
    parametersBDDC->set("Krylov Method", 1);
  }

  parametersBDDC->set("MueLu Parameter List", parametersMueLu);
  parametersBDDC->set("NodalAMG Parameter List", parametersNodalAMG);
}

bool usingMueLu(RCP<Teuchos::ParameterList> & paramsBDDC) 
{
  bool useMueLu = false;
  if (paramsBDDC->get("Dirichlet Solver", "a") == "MueLu") useMueLu = true;
  if (paramsBDDC->get("Neumann Solver", "a") == "MueLu") useMueLu = true;
  if (paramsBDDC->get("Coarse Solver", "a") == "MueLu") useMueLu = true;
  return useMueLu;
}

void checkInterfacePreconditionerSolverCompatibility
  (const LO myPID,
   Teuchos::ParameterList & paramListBDDC)
{
  bool interfacePreconditioner = 
    paramListBDDC.get("Interface Preconditioner", true);
  if (interfacePreconditioner) {
    const std::string solverDirichlet =
      paramListBDDC.get("Dirichlet Solver", "SuperLU");
    if ((solverDirichlet == "MueLu") || (solverDirichlet == "NodalAMG")) {
      if (myPID == 0) {
	std::cout << "Warning: cannot use interface preconditioner with\n";
	std::cout << solverDirichlet << " Dirichlet solver, switching to non-interface preconditioner\n";
      }
      paramListBDDC.set("Interface Preconditioner", false);
    }
  }
}

};

} // end namespace

#endif // BDDC_SETUPTEST_H

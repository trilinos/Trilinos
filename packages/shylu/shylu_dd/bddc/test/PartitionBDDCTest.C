
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

#include <gtest/gtest.h>
#include <mpi.h>
#include "ProblemMakerBDDC.h"
#include "shylu_PartitionOfUnityBDDC.h"
#include "shylu_UtilBDDC.h"
#include "shylu_enumsBDDC.h"

using Teuchos::RCP;

namespace {

typedef int LO; // Local Ordinal
typedef int GO; // Global Ordinal
typedef double SX; // floating point data type
typedef double SM; // real (magnitude) for SX
    
void outputSubdomainMatrices(std::vector< std::vector<LO> > & subRowBegin,
			     std::vector< std::vector<LO> > & subColumns,
			     std::vector< std::vector<SX> > & subValues,
			     LO myPID,
			     const char* fnameBase);

TEST(PartitionBDDC, Test1)
{
  int numProc, myPID;
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Comm_size(Comm, &numProc);
  MPI_Comm_rank(Comm, &myPID);
  if (numProc != 4) return;
  double lengthDir1 = 1;
  double lengthDir2 = 1;
  double lengthDir3 = 1;
  LO numSubDir1 = 2;
  LO numSubDir2 = 2;
  LO numSubDir3 = 1;
  LO numSubDir1PerProc = 2;
  LO numSubDir2PerProc = 2;
  LO numSubDir3PerProc = 1;
  LO numElemPerSubDir1 = 3*numSubDir1PerProc;
  LO numElemPerSubDir2 = 3*numSubDir2PerProc;
  LO numElemPerSubDir3 = 1*numSubDir3PerProc;
  
  enum bddc::ProblemType problemType = bddc::SCALARPDE; 
  // ELASTICITY or SCALARPDE
  int spatialDim(2), loadDirection(0);
  enum bddc::AnalysisType analysisType = bddc::STANDARD;
  RCP<Teuchos::ParameterList> Parameters;
  Parameters = Teuchos::rcp( new Teuchos::ParameterList() );
  Parameters->set("Spatial Dimension", spatialDim);
  Parameters->set("Problem Type", problemType);
  Parameters->set("Analysis Type", analysisType);
  Parameters->set("Length Direction 1", lengthDir1);
  Parameters->set("Length Direction 2", lengthDir2);
  Parameters->set("Length Direction 3", lengthDir3);
  Parameters->set("Number of Subdomains Direction 1", numSubDir1);
  Parameters->set("Number of Subdomains Direction 2", numSubDir2);
  Parameters->set("Number of Subdomains Direction 3", numSubDir3);
  Parameters->set("Number of Elements Per Subdomain Direction 1",
		  numElemPerSubDir1);
  Parameters->set("Number of Elements Per Subdomain Direction 2",
		  numElemPerSubDir2);
  Parameters->set("Number of Elements Per Subdomain Direction 3",
		  numElemPerSubDir3);
  
  Parameters->set("Apply Left Side Essential BCs", false);
  Parameters->set("Apply Right Side Essential BCs", false);

  Parameters->set("Load Direction", loadDirection);
  Parameters->set("Artificial Foundation Stiffness", 0);
  Parameters->set("omega", 0.0);
  Parameters->set("Generate Constraint Equations", false);
  Parameters->set("Interface Preconditioner", true);
  Parameters->set("Print Interior Matrices", false);

  assert (numProc == numSubDir1*numSubDir2*numSubDir3);
  if (spatialDim == 2) assert (numProc == 4);
  //
  // generate problem
  //
  RCP< bddc::ProblemMaker<LO,GO,SX,SM> > Problem = 
    rcp( new bddc::ProblemMaker<LO,GO,SX,SM>(Parameters, Comm) );
  std::vector< std::vector<LO> > subNodes, subNodeBegin, subRowBegin, 
    subLocalDofs, subColumns, subElems;
  std::vector< std::vector<SX> > subValues;
  Problem->getSubDomainElements(numSubDir1PerProc, numSubDir2PerProc,
				numSubDir3PerProc, subElems);
  Problem->getSubDomainNodeData(subElems, subNodes, subNodeBegin,
				subLocalDofs);
  Problem->getSubdomainMatrices(subElems, subNodes, subNodeBegin, subRowBegin, 
				subColumns, subValues);
  LO outputProc = -1;
  if (myPID == outputProc) {
    outputSubdomainMatrices(subRowBegin, subColumns, subValues, myPID, "A");
  }
  LO numNode = Problem->getNumNode();
  LO numSub = subRowBegin.size();
  std::vector< LO* > subRowBeginPtr(numSub), subColumnsPtr(numSub);
  std::vector< SX* > subValuesPtr(numSub);
  for (LO i=0; i<numSub; i++) {
    subRowBeginPtr[i] = subRowBegin[i].data();
    subColumnsPtr[i] = subColumns[i].data();
    subValuesPtr[i] = subValues[i].data();
  }
  const GO* nodeGlobalIDs = Problem->getNodeGlobalIDs();
  RCP<const Teuchos::Comm<int> > tComm = rcp( new Teuchos::MpiComm<int>(Comm) );
  RCP< bddc::PartitionOfUnity<SX,SM,LO,GO> > Partition =
    rcp( new bddc::PartitionOfUnity<SX,SM,LO,GO>
	 (numNode, nodeGlobalIDs, subNodeBegin, subNodes, 
	  subRowBeginPtr.data(), subColumnsPtr.data(), subValuesPtr.data(), 
	  spatialDim, Parameters, tComm) );
  const std::vector< std::vector<LO> > & subdomainEquivClasses = 
    Partition->getSubdomainEquivClasses();
  const std::vector< std::vector<LO> > & equivClasses = 
    Partition->getEquivClasses();
  if (spatialDim == 2) {
    if (myPID == 0) {
      EXPECT_EQ(subdomainEquivClasses[0].size(), 5);
      EXPECT_EQ(subdomainEquivClasses[1].size(), 7);
      EXPECT_EQ(subdomainEquivClasses[2].size(), 7);
      EXPECT_EQ(subdomainEquivClasses[3].size(), 8);
      EXPECT_EQ(equivClasses.size(), 16);
      for (int i=0; i<8; i++) {
	EXPECT_EQ(equivClasses[i].size(), 2);
	EXPECT_EQ(equivClasses[i+8].size(), 1);
      }
    }
  }
}

void outputSubdomainMatrices(std::vector< std::vector<LO> > & subRowBegin,
			     std::vector< std::vector<LO> > & subColumns,
			     std::vector< std::vector<SX> > & subValues,
			     LO myPID,
			     const char* fnameBase)
{
  LO numSub = subRowBegin.size();
  SX value(0);
  bool isComplex = bddc::UtilBDDC<SX,SM>::isComplex(value);
  for (LO sub=0; sub<numSub; sub++) {
    std::vector<LO> & rowBegin = subRowBegin[sub];
    std::vector<LO> & columns = subColumns[sub];
    std::vector<SX> & values = subValues[sub];
    char fname[101];
    sprintf(fname, "%s_proc%d_sub%d.dat", fnameBase, myPID, sub);
    std::ofstream fout;
    fout.open(fname);
    LO numRows = rowBegin.size()-1;
    for (LO i=0; i<numRows; i++) {
      for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	fout << i+1 << "  " << columns[j]+1 << " ";
	fout << std::setw(22) << std::setprecision(15);
	value = values[j];
	fout << bddc::UtilBDDC<SX,SM>::real(value);
	if (isComplex == true) {
	  fout << " " << bddc::UtilBDDC<SX,SM>::imag(value);
	}
	fout << std::endl;
      }
    }
    fout.close();
  }
}

} // end namespace

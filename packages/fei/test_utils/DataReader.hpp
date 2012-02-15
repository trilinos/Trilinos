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

#ifndef _DataReader_h_
#define _DataReader_h_


#include <test_utils/BCNodeSet.hpp>
#include <test_utils/CRSet.hpp>
#include <test_utils/CommNodeSet.hpp>
#include <test_utils/ElemBlock.hpp>
#include <test_utils/AccessPattern.hpp>
#include <test_utils/CoefAccess.hpp>
#include <fei_iostream.hpp>
#include <string>

class DataReader {
 public:
  DataReader();
  ~DataReader();

  int readData(const char* fileName);

  int solveType_;

  std::string solverLibraryName_;
  std::string solnFileName_;
  std::string checkFileName_;

  int numFields_;
  int* fieldIDs_;
  int* fieldSizes_;

  int numParams_;
  char** paramStrings_;

  int numElemBlocks_;
  ElemBlock* elemBlocks_; //list of length numElemBlocks_

  int numCoefAccessPatterns_;
  AccessPattern* accessPatterns_;

  int numCoefAccesses_;
  CoefAccess* coefAccesses_;

  int numCRMultSets_;
  CRSet* crMultSets_;

  int numSlaveVars_;
  CRSet* slaveVars_;

   int numCRPenSets_;
   CRSet* crPenSets_;

   int numBCNodeSets_;
   BCNodeSet* bcNodeSets_;

   int numSharedNodeSets_;
   CommNodeSet* sharedNodeSets_;

   int getFieldSize(int fieldID);

   static int getKeyword(FEI_ISTREAM* instr, char*& keyword);
   void readData(FEI_ISTREAM* instr, char* keyword);
   static void readData(FEI_ISTREAM* instr, int& n);
   static void readData(FEI_ISTREAM* instr, double& val);

   static int is_reg_char(char c);
   static int skipWhite(FEI_ISTREAM* instr);

 private:
   void deleteMemory();

   bool numFieldsRead_;
   bool numElemBlocksRead_;
   int currentElemBlockIndex_;
   int currentElemIndex_;

   int currentShIndex_;
   int currentExtIndex_;
   int currentBCIndex_;
};

#endif


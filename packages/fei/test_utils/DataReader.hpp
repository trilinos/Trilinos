#ifndef _DataReader_h_
#define _DataReader_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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


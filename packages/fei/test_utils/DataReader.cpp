/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <cstring>

#include <fei_fstream.hpp>
#include <fei_iostream.hpp>
#include <fei_defs.h>
#include <test_utils/BCNodeSet.hpp>
#include <test_utils/CRSet.hpp>
#include <test_utils/CommNodeSet.hpp>
#include <test_utils/ElemBlock.hpp>
#include <test_utils/AccessPattern.hpp>
#include <test_utils/DataReader.hpp>

//==============================================================================
DataReader::DataReader()
 : 
   solveType_(0),
   solverLibraryName_(),
   numFields_(0),
   fieldIDs_(NULL),
   fieldSizes_(NULL),
   numParams_(0),
   paramStrings_(NULL),
   numElemBlocks_(0),
   elemBlocks_(NULL),
   numCoefAccessPatterns_(0),
   accessPatterns_(NULL),
   numCoefAccesses_(0),
   coefAccesses_(NULL),
   numCRMultSets_(0),
   crMultSets_(NULL),
   numSlaveVars_(0),
   slaveVars_(NULL),
   numCRPenSets_(0),
   crPenSets_(NULL),
   numBCNodeSets_(0),
   bcNodeSets_(NULL),
   numSharedNodeSets_(0),
   sharedNodeSets_(NULL),
   numFieldsRead_(false),
   numElemBlocksRead_(false),
   currentElemBlockIndex_(0),
   currentElemIndex_(0),
   currentShIndex_(0),
   currentExtIndex_(0),
   currentBCIndex_(0)
{
}

//==============================================================================
DataReader::~DataReader() {
   deleteMemory();

   numElemBlocksRead_ = false;
   numFieldsRead_ = false;
}

//==============================================================================
void DataReader::deleteMemory() {
   for(int i=0; i<numParams_; i++) {
      delete [] paramStrings_[i];
   }
   delete [] paramStrings_;
   numParams_ = 0;

   delete [] accessPatterns_;
   numCoefAccessPatterns_ = 0;

   delete [] coefAccesses_;
   numCoefAccesses_ = 0;

   delete [] elemBlocks_;
   numElemBlocks_ = 0;

   delete [] fieldIDs_;
   delete [] fieldSizes_;
   numFields_ = 0;

   delete [] sharedNodeSets_;
   numSharedNodeSets_ = 0;

   delete [] crMultSets_;
   numCRMultSets_ = 0;

   delete [] slaveVars_;
   numSlaveVars_ = 0;

   delete [] crPenSets_;
   numCRPenSets_ = 0;

   delete [] bcNodeSets_;
   numBCNodeSets_ = 0;
}

//==============================================================================
int DataReader::readData(const char* fileName) {

   FEI_IFSTREAM* instr = new FEI_IFSTREAM(fileName);

   if (instr->bad()) {
      fei::console_out() << "DataReader::readData: ERROR opening " << fileName << FEI_ENDL;
      return(1);
   }
// FEI_COUT << "DataReader reading from " << fileName << FEI_ENDL;
   char* keyword = NULL;

   int err = getKeyword(instr, keyword);

   while (!instr->eof() && !err) {
      readData(instr, keyword);
      delete [] keyword;
      keyword = NULL;

      err = getKeyword(instr, keyword);
   }

   delete instr;
   delete [] keyword;

   return(0);
}

//==============================================================================
int DataReader::getKeyword(FEI_ISTREAM* instr, char*& keyword) {
   int err = skipWhite(instr);
   if (err) return(err);

   char *temp = new char[256];
   for(int i=0; i<256; i++) temp[i] = '\0';

   do {
      (*instr) >> temp;
   } while ((std::strlen(temp) == 0) && (!instr->eof()));

   keyword = new char[std::strlen(temp)+1];
   std::strcpy(keyword, temp);

   delete [] temp;

   return(0);
}

//==============================================================================
int DataReader::is_reg_char(char c) {
   int i = (int)c;
   if (i<1 || i>126) return(0);

   return(1);
}

//==============================================================================
int DataReader::skipWhite(FEI_ISTREAM* instr) {
   char c = '\0';
   instr->get(c);

   if (!is_reg_char(c)) {
      return(1);
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
         return(1);
      }
   }

   instr->putback(c);
   return(0);
}

//==============================================================================
void DataReader::readData(FEI_ISTREAM* instr, char* keyword) {

   if (!std::strcmp("solveType", keyword)) {
      readData(instr, solveType_);
      return;
   }

   if (!std::strcmp("parameters", keyword)) {
      int tmp = 0;
      readData(instr, tmp);

      char** newParams = new char*[numParams_ + tmp];
      for(int pp=0; pp<numParams_; pp++) newParams[pp] = paramStrings_[pp];

      for(int i=numParams_; i<numParams_+tmp; i++) {
         char* buf = new char[256];
         for(int j=0; j<256; j++) buf[j] = '\0';

         skipWhite(instr);
         instr->getline(buf, 128);

         newParams[i] = new char[std::strlen(buf)+2];
         std::strcpy(newParams[i], buf);

         delete [] buf;
      }

      delete [] paramStrings_;
      paramStrings_ = newParams;
      numParams_ += tmp;

      return;
   }

   if (!std::strcmp("numFields", keyword)) {
      readData(instr, numFields_);

      fieldSizes_ = new int[numFields_];
      fieldIDs_ = new int[numFields_];

      numFieldsRead_ = true;
      return;
   }

   if (!std::strcmp("fieldIDs", keyword)) {
      for(int i=0; i<numFields_; i++) readData(instr, fieldIDs_[i]);
      return;
   }

   if (!std::strcmp("fieldSizes", keyword)) {
      for(int i=0; i<numFields_; i++) readData(instr, fieldSizes_[i]);
      return;
   }

   if (!std::strcmp("numElemBlocks", keyword)) {
      if (numElemBlocks_ > 0) {
         FEI_COUT << "DataReader: Caution, re-setting numElemBlocks." << FEI_ENDL;
         delete [] elemBlocks_;
      }

      readData(instr, numElemBlocks_);

      elemBlocks_ = new ElemBlock[numElemBlocks_];

      numElemBlocksRead_ = true;
      currentElemBlockIndex_ = -1;

      return;
   }

   if (!std::strcmp("blockID", keyword)) {
      currentElemBlockIndex_++;
      currentElemIndex_ = 0;

      if (!numElemBlocksRead_) {
         FEI_COUT << "DataReader: ERROR, numElemBlocks not read before blockID."
            << FEI_ENDL;
         return;
      }

      ElemBlock& eb1 = elemBlocks_[currentElemBlockIndex_];

      int tmp;
      readData(instr, tmp);
      eb1.blockID_ = (GlobalID)tmp;

      return;
   }

   if (!std::strcmp("interleaveStrategy", keyword)) {
      ElemBlock& eb2 = elemBlocks_[currentElemBlockIndex_];
      int interleave;

      readData(instr, interleave);

      eb2.interleaveStrategy_ = interleave;

      return;
   }

   if (!std::strcmp("numElements", keyword)) {
      ElemBlock& eb3 = elemBlocks_[currentElemBlockIndex_];
      int numElems;

      readData(instr, numElems);

      eb3.numElements_ = numElems;
      eb3.elemIDs_ = new GlobalID[numElems];
      eb3.elemConn_ = new GlobalID*[numElems];
      eb3.elemStiff_ = new double**[numElems];
      eb3.elemLoad_ = new double*[numElems];

      return;
   }

   if (!std::strcmp("numNodesPerElement", keyword)) {
      ElemBlock& eb4 = elemBlocks_[currentElemBlockIndex_];
      int numNodes;

      readData(instr, numNodes);

      eb4.numNodesPerElement_ = numNodes;
      eb4.numFieldsPerNode_ = new int[numNodes];
      eb4.nodalFieldIDs_ = new int*[numNodes];

      for(int i=0; i<eb4.numElements_; i++) {
         eb4.elemConn_[i] = new GlobalID[numNodes];
      }
      
      return;
   }

   if (!std::strcmp("numElemDOF", keyword)) {
      ElemBlock& eb5 = elemBlocks_[currentElemBlockIndex_];
      int edof;
      readData(instr, edof);
      eb5.numElemDOF_ = edof;

      if (edof > 0) {
         eb5.elemDOFFieldIDs_ = new int[edof];
      }

      return;
   }

   if (!std::strcmp("elemDOFFieldIDs", keyword)) {
      ElemBlock& eb6 = elemBlocks_[currentElemBlockIndex_];
      int edof = eb6.numElemDOF_;

      for(int i=0; i<edof; i++) {
         int eDofFieldID;
         readData(instr, eDofFieldID);
         eb6.elemDOFFieldIDs_[i] = eDofFieldID;
      }
      return;
   }

   if (!std::strcmp("elemFormat", keyword)) {
      ElemBlock& eb7 = elemBlocks_[currentElemBlockIndex_];
      int ef;
      readData(instr, ef);

      eb7.elemFormat_ = ef;
      return;
   }

   if (!std::strcmp("numFieldsPerNode", keyword)) {
      ElemBlock& eb8 = elemBlocks_[currentElemBlockIndex_];

      int i;
      for(i=0; i<eb8.numNodesPerElement_; i++) {
         int nf;
         readData(instr, nf);
         eb8.numFieldsPerNode_[i] = nf;
         eb8.nodalFieldIDs_[i] = new int[nf];
      }

      return;
   }

   if (!std::strcmp("nodalFieldIDs", keyword)) {
      ElemBlock& eb9 = elemBlocks_[currentElemBlockIndex_];

      int i, numStiffRows = 0;
      for(i=0; i<eb9.numNodesPerElement_; i++) {
         for(int j=0; j<eb9.numFieldsPerNode_[i]; j++) {
            int nfid;
            readData(instr, nfid);
            eb9.nodalFieldIDs_[i][j] = nfid;

            numStiffRows += getFieldSize(nfid);
         }
      }

      numStiffRows += eb9.numElemDOF_;

      eb9.numStiffRows_ = numStiffRows;

      for(i=0; i<eb9.numElements_; i++) {
         eb9.elemStiff_[i] = new double*[numStiffRows];
         eb9.elemLoad_[i] = new double[numStiffRows];
         for(int j=0; j<numStiffRows; j++) {
            eb9.elemStiff_[i][j] = new double[numStiffRows];
         }
      }

      return;
   }

   if (!std::strcmp("elemID", keyword)) {
      ElemBlock& eb10 = elemBlocks_[currentElemBlockIndex_];

      int i, tmp;
      readData(instr, tmp);
      eb10.elemIDs_[currentElemIndex_] = (GlobalID)tmp;

      int* conn = eb10.elemConn_[currentElemIndex_];

      for(i=0; i<eb10.numNodesPerElement_; i++) {
         readData(instr, conn[i]);
      }

      double** stiff = eb10.elemStiff_[currentElemIndex_];

      for(i=0; i<eb10.numStiffRows_; i++) {
         for(int j=0; j<eb10.numStiffRows_; j++) {
            readData(instr, stiff[i][j]);
         }
      }

      double* load = eb10.elemLoad_[currentElemIndex_];

      for(i=0; i<eb10.numStiffRows_; i++) {
         readData(instr, load[i]);
      }

      currentElemIndex_++;

      return;
   }

   if (!std::strcmp("numSharedNodeSets", keyword)) {
      readData(instr, numSharedNodeSets_);

      if (numSharedNodeSets_ == 0) return;

      sharedNodeSets_ = new CommNodeSet[numSharedNodeSets_];

      currentShIndex_ = 0;

      for(int i=0; i<numSharedNodeSets_; i++) {
         CommNodeSet& shSet = sharedNodeSets_[currentShIndex_++];

         int j, nn;
         readData(instr, nn);
         shSet.numNodes_ = nn;

         shSet.nodeIDs_ = new GlobalID[nn];

         for(j=0; j<nn; j++) {
            int tmp;
            readData(instr, tmp);
            shSet.nodeIDs_[j] = (GlobalID)tmp;
         }

         shSet.procsPerNode_ = new int[nn];

         for(j=0; j<nn; j++) {
            int tmp;
            readData(instr, tmp);
            shSet.procsPerNode_[j] = tmp;
         }

         shSet.procs_ = new int*[nn];

         for(j=0; j<nn; j++) {
            shSet.procs_[j] = new int[shSet.procsPerNode_[j]];
            for(int k=0; k<shSet.procsPerNode_[j]; k++) {
               readData(instr, shSet.procs_[j][k]);
            }
         }
      }

      return;
   }

   if (!std::strcmp("numBCNodeSets", keyword)) {
      readData(instr, numBCNodeSets_);

      if (numBCNodeSets_ == 0) return;

      bcNodeSets_ = new BCNodeSet[numBCNodeSets_];

      currentBCIndex_ = 0;

      for(int i=0; i<numBCNodeSets_; i++) {
         BCNodeSet& bcSet = bcNodeSets_[currentBCIndex_++];

         int j, nn;
         readData(instr, nn);
         bcSet.numNodes_ = nn;

         readData(instr, bcSet.fieldID_);

         int field_offset;
         readData(instr, field_offset);
         
         bcSet.nodeIDs_ = new GlobalID[nn];
         bcSet.offsetsIntoField_ = new int[nn];
         bcSet.prescribed_values_ = new double[nn];

         for(j=0; j<nn; ++j) bcSet.offsetsIntoField_[j] = field_offset;

         for(j=0; j<nn; j++) {
            readData(instr, bcSet.nodeIDs_[j]);
            readData(instr, bcSet.prescribed_values_[j]);
         }
      }

      return;
   }

   if (!std::strcmp("numCRMultSets", keyword)) {
      readData(instr, numCRMultSets_);

      if (numCRMultSets_ == 0) return;

      crMultSets_ = new CRSet[numCRMultSets_];

      for(int i=0; i<numCRMultSets_; i++) {
         CRSet& cr1 = crMultSets_[i];

	 int dummy;
         readData(instr, dummy);//used to be numCRs_
         readData(instr, cr1.numNodes_);

         cr1.nodeIDs_ = new GlobalID*[1];
         cr1.fieldIDs_ = new int[cr1.numNodes_];

         int j;
         for(j=0; j<1; j++) {
            cr1.nodeIDs_[j] = new GlobalID[cr1.numNodes_];

            for(int k=0; k<cr1.numNodes_; k++) {
               readData(instr, cr1.nodeIDs_[j][k]);
            }
         }

         for(j=0; j<cr1.numNodes_; j++) {
            readData(instr, cr1.fieldIDs_[j]);
         }

         int len = 0;
         for(j=0; j<cr1.numNodes_; j++) {
            len += getFieldSize(cr1.fieldIDs_[j]);
         }
         cr1.weights_ = new double[len];

         int offset = 0;
         for(j=0; j<cr1.numNodes_; j++) {
            int size = getFieldSize(cr1.fieldIDs_[j]);

            for(int k=0; k<size; k++) {
               readData(instr, cr1.weights_[offset++]);
            }
         }

         cr1.values_ = new double[1];
         for(j=0; j<1; j++) {
            readData(instr, cr1.values_[j]);
         }
      }

      return;
   }

   if (!std::strcmp("numCoefAccessPatterns", keyword)) {
     readData(instr, numCoefAccessPatterns_);

     if (numCoefAccessPatterns_ == 0) return;

     accessPatterns_ = new AccessPattern[numCoefAccessPatterns_];
     for(int i=0; i<numCoefAccessPatterns_; i++) {
       AccessPattern& patt = accessPatterns_[i];

       readData(instr, patt.ID_);
       readData(instr, patt.numRowIDs_);
       if (patt.numRowIDs_ <= 0) {
	 fei::console_out() << "DataReader ERROR, numRowIDs_ <=0"<<FEI_ENDL;
	 return;
       }

       patt.numFieldsPerRow_ = new int[patt.numRowIDs_];
       int j;
       for(j=0; j<patt.numRowIDs_; j++) {
	 readData(instr, patt.numFieldsPerRow_[j]);
       }

       patt.rowFieldIDs_ = new int*[patt.numRowIDs_];
       for(j=0; j<patt.numRowIDs_; j++) {
	 patt.rowFieldIDs_[j] = new int[patt.numFieldsPerRow_[j]];
       }

       for(int r=0; r<patt.numRowIDs_; r++) {
	 for(int c=0; c<patt.numFieldsPerRow_[r]; c++) {
	   readData(instr, patt.rowFieldIDs_[r][c]);
	 }
       }

       readData(instr, patt.numColIDsPerRow_);
       if (patt.numColIDsPerRow_ <= 0) {
	 fei::console_out() << "DataReader ERROR, numColIDsPerRow_ <=0"<<FEI_ENDL;
	 return;
       }
       
       patt.numFieldsPerCol_ = new int[patt.numColIDsPerRow_];
       for(j=0; j<patt.numColIDsPerRow_; j++) {
	 readData(instr, patt.numFieldsPerCol_[j]);
       }

       patt.colFieldIDs_ = new int*[patt.numColIDsPerRow_];
       for(j=0; j<patt.numColIDsPerRow_; j++) {
	 patt.colFieldIDs_[j] = new int[patt.numFieldsPerCol_[j]];
       }

       for(int rr=0; rr<patt.numColIDsPerRow_; rr++) {
	 for(int c=0; c<patt.numFieldsPerCol_[rr]; c++) {
	   readData(instr, patt.colFieldIDs_[rr][c]);
	 }
       }

       readData(instr, patt.interleaveStrategy_);
     }

     return;
   }

   if (!std::strcmp("coefAccess", keyword)) {
     int i, patternID = -1;
     readData(instr, patternID);

     //find the access-pattern corresponding to this coef-access.
     int index = -1;
     for(i=0; i<numCoefAccessPatterns_; i++) {
       if (accessPatterns_[i].ID_ == patternID) index = i;
     }

     if (index < 0) {
       fei::console_out() << "DataReader ERROR, patternID " << patternID << " not found."<<FEI_ENDL;
       return;
     }

     AccessPattern& patt = accessPatterns_[index];

     //now lengthen the list of coef-accesses.
     CoefAccess* newAccesses = new CoefAccess[numCoefAccesses_+1];
     for(i=0; i<numCoefAccesses_; i++) newAccesses[i] = coefAccesses_[i];

     delete [] coefAccesses_;
     coefAccesses_ = newAccesses;

     CoefAccess& coefAcc = coefAccesses_[numCoefAccesses_];

     coefAcc.patternID_ = patternID;

     coefAcc.numRowIDs_ = patt.numRowIDs_;
     coefAcc.numColIDsPerRow_ = patt.numColIDsPerRow_;

     if (coefAcc.numRowIDs_ <= 0 || coefAcc.numColIDsPerRow_ <= 0) {
       fei::console_out() << "DataReader ERROR, coef-access has 0 rows or cols." << FEI_ENDL;
       return;
     }

     coefAcc.rowIDs_ = new GlobalID[coefAcc.numRowIDs_];
     for(i=0; i<coefAcc.numRowIDs_; i++) {
       readData(instr, coefAcc.rowIDs_[i]);
     }

     int len = coefAcc.numRowIDs_ * coefAcc.numColIDsPerRow_;
     coefAcc.colIDs_ = new GlobalID[len];
     for(i=0; i<len; i++) {
       readData(instr, coefAcc.colIDs_[i]);
     }

     //calculate numRowCoefs.
     coefAcc.numRowCoefs_ = 0;
     for(i=0; i<coefAcc.numRowIDs_; i++) {
       for(int j=0; j<patt.numFieldsPerRow_[i]; j++) {
	 coefAcc.numRowCoefs_ += getFieldSize(patt.rowFieldIDs_[i][j]);
       }
     }

     //calculate numColCoefs.
     coefAcc.numColCoefs_ = 0;
     for(i=0; i<coefAcc.numColIDsPerRow_; i++) {
       for(int j=0; j<patt.numFieldsPerCol_[i]; j++) {
	 coefAcc.numColCoefs_ += getFieldSize(patt.colFieldIDs_[i][j]);
       }
     }

     len = coefAcc.numRowCoefs_*coefAcc.numColCoefs_;
     coefAcc.coefs_ = new double[len];
     int offset = 0;
     //now read the coefficients-table.
     for(i=0; i<len; i++) {
       readData(instr, coefAcc.coefs_[offset++]);
     }

     numCoefAccesses_++;
     return;
   }

   if (!std::strcmp("numSlaveVariables", keyword)) {
      readData(instr, numSlaveVars_);

      if (numSlaveVars_ == 0) return;

      slaveVars_ = new CRSet[numSlaveVars_];

      for(int i=0; i<numSlaveVars_; i++) {
         CRSet& cr2 = slaveVars_[i];

         readData(instr, cr2.slaveNodeID_);
         readData(instr, cr2.slaveFieldID_);
	 readData(instr, cr2.slaveOffset_);

	 readData(instr, cr2.numNodes_);
         cr2.nodeIDs_ = new GlobalID*[1];
	 cr2.nodeIDs_[0] = new GlobalID[cr2.numNodes_];
         cr2.fieldIDs_ = new int[cr2.numNodes_];

         int j;
	 for(int k=0; k<cr2.numNodes_; k++) {
	   readData(instr, cr2.nodeIDs_[0][k]);
	 }

         for(j=0; j<cr2.numNodes_; j++) {
            readData(instr, cr2.fieldIDs_[j]);
         }

         int len = 0;
         for(j=0; j<cr2.numNodes_; j++) {
            len += getFieldSize(cr2.fieldIDs_[j]);
         }
         cr2.weights_ = new double[len];

         int offset = 0;
         for(j=0; j<cr2.numNodes_; j++) {
            int size = getFieldSize(cr2.fieldIDs_[j]);

            for(int k=0; k<size; k++) {
               readData(instr, cr2.weights_[offset++]);
            }
         }

	 cr2.values_ = new double[1];
	 readData(instr, cr2.values_[0]);
      }

      return;
   }

   if (!std::strcmp("numCRPenSets", keyword)) {
      readData(instr, numCRPenSets_);

      if (numCRPenSets_ == 0) return;

      crPenSets_ = new CRSet[numCRPenSets_];

      for(int i=0; i<numCRPenSets_; i++) {
         CRSet& cr3 = crPenSets_[i];

	 int dummy;
         readData(instr, dummy);//used to be numCRs_
         readData(instr, cr3.numNodes_);

         cr3.nodeIDs_ = new GlobalID*[1];
         cr3.fieldIDs_ = new int[cr3.numNodes_];

         int j;
         for(j=0; j<1; j++) {
            cr3.nodeIDs_[j] = new GlobalID[cr3.numNodes_];

            for(int k=0; k<cr3.numNodes_; k++) {
               readData(instr, cr3.nodeIDs_[j][k]);
            }
         }

         for(j=0; j<cr3.numNodes_; j++) {
            readData(instr, cr3.fieldIDs_[j]);
         }

         int len3 = 0;
         for(j=0; j<cr3.numNodes_; j++) {
            len3 += getFieldSize(cr3.fieldIDs_[j]);
         }
         cr3.weights_ = new double[len3];

         int offset3 = 0;
         for(j=0; j<cr3.numNodes_; j++) {
            int size3 = getFieldSize(cr3.fieldIDs_[j]);

            for(int k=0; k<size3; k++) {
	      double dummy3;
	      readData(instr, dummy3);
	      cr3.weights_[offset3++] = dummy3;
            }
         }

         cr3.values_ = new double[1];
         for(j=0; j<1; j++) {
            readData(instr, cr3.values_[j]);
         }

         cr3.penValues_ = new double[1];
         for(j=0; j<1; j++) {
            readData(instr, cr3.penValues_[j]);
         }
      }

      return;
   }
}

//==============================================================================
int DataReader::getFieldSize(int fieldID) {
   for(int i=0; i<numFields_; i++) {
      if (fieldID == fieldIDs_[i]) return(fieldSizes_[i]);
   }

   fei::console_out() << "DataReader: ERROR, trying to find size of non-existent field."
        << FEI_ENDL;
   return(0);
}

//==============================================================================
void DataReader::readData(FEI_ISTREAM* instr, int& n) {
   int err = skipWhite(instr);
   if (err) return;
   (*instr) >> n;
}

//==============================================================================
void DataReader::readData(FEI_ISTREAM* instr, double& val) {
   int err = skipWhite(instr);
   if (err) return;
   (*instr) >> val;
}


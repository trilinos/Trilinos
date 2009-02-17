/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _driverData_h_
#define _driverData_h_

#include <fei_macros.hpp>
#include <fei_defs.h>

class initElem {
 public:
  initElem(){}
  ~initElem()
    {
      delete [] nodeIDs;
    }
  GlobalID elemBlockID;
  GlobalID elemID;
  int numNodes;
  GlobalID* nodeIDs;
};

class sumInElem {
 public:
  sumInElem() : stiff1D(NULL), stiffness(NULL), load(NULL) {}
  ~sumInElem()
    {
      delete [] nodeIDs;
      delete [] stiff1D;
      delete [] stiffness;
      delete [] load;
    }
  GlobalID elemBlockID;
  GlobalID elemID;
  int numNodes;
  GlobalID* nodeIDs;
  int numRows;
  double* stiff1D;
  double** stiffness;
  double* load;
  int elemFormat;
};

class nodeBC {
 public:
  nodeBC():nodeIDs(NULL), alpha(NULL), beta(NULL), gamma(NULL) {}
  ~nodeBC()
    {
      for(int i=0; i<numNodes; ++i) {
	delete [] alpha[i];
	delete [] beta[i];
	delete [] gamma[i];
      }
      delete [] alpha;
      delete [] beta;
      delete [] gamma;
      delete [] nodeIDs;
    }
  int numNodes;
  GlobalID* nodeIDs;
  int fieldID;
  int fieldSize;
  double** alpha;
  double** beta;
  double** gamma;
};

class initCR {
 public:
  initCR(){}
  ~initCR()
    {
      delete [] nodeIDs;
      delete [] fieldIDs;
    }

  int numNodes;
  GlobalID* nodeIDs;
  int* fieldIDs;
  int CRID;
};

class loadCR {
 public:
  loadCR() {}
  ~loadCR()
    {
      delete [] fieldIDs;
      delete [] fieldSizes;
      delete [] weights;
      delete [] nodeIDs;
    }

  int numNodes;
  GlobalID* nodeIDs;
  int* fieldIDs;
  int* fieldSizes;
  double* weights;
  double CRValue;
  double penValue;
  int CRID;  
};

class sharedNodes {
 public:
  sharedNodes() : nodeIDs(NULL), numProcsPerNode(NULL), sharedProcIDs(NULL){}
  ~sharedNodes()
    {
      for(int i=0; i<numNodes; ++i) delete [] sharedProcIDs[i];
      delete [] sharedProcIDs;
      delete [] numProcsPerNode;
      delete [] nodeIDs;
    }

  int numNodes;
  GlobalID* nodeIDs;
  int* numProcsPerNode;
  int** sharedProcIDs;
};

class parameters {
 public:
  parameters() : paramList() {}
  ~parameters()
    {
      for(size_t i=0; i<paramList.size(); ++i) delete [] paramList[i];
    }

  std::vector<char*> paramList;
};

class setIDLists {
 public:
  setIDLists() : matrixIDs(NULL), rhsIDs(NULL) {}
  ~setIDLists(){delete [] matrixIDs; delete [] rhsIDs;}

  int* matrixIDs;
  int numMatrices;
  int* rhsIDs;
  int numRHSs;
};

class putBlockFieldNodeSolution {
 public:
  putBlockFieldNodeSolution() : nodeIDs(NULL), estimates(NULL) {}
  ~putBlockFieldNodeSolution(){delete [] nodeIDs; delete [] estimates;}

  int elemBlockID;
  int fieldID;
  int fieldSize;
  int numNodes;
  GlobalID* nodeIDs;
  double* estimates;
};

class driverData {
 public:
  driverData();
  ~driverData();

  int readData(const char* fileName);

  /** call a named method on the FEI object.
      @return error-code 0 if successful. negative if an error occurs, positive
      if the named method is un-recognized.
  */
  int call_fei_method(const char* method, FEI* fei);

  std::vector<const char*>& get_methodNames() { return( methodNames ); }

 private:
  int readData(FEI_ISTREAM* instr, char* keyword);
  int getKeyword(FEI_ISTREAM* instr, char*& keyword);

  int is_reg_char(char c);
  int skipWhite(FEI_ISTREAM* instr);
  int readData(FEI_ISTREAM* instr, int& n);
  int readData(FEI_ISTREAM* instr, double& val);
  int appendName(const char* name);

  std::vector<const char*> methodNames;
  char* temp_;
  int tempLen_;

  int solveType_;

  int initFields_numFields_;
  int* initFields_fieldSizes_;
  int* initFields_fieldIDs_;

  int initElemBlock_numInts_;
  int* initElemBlock_ints_;
  int* initElemBlock_fieldsPerNode_;
  int** initElemBlock_fieldIDs_;
  int* initElemBlock_elemDofFieldIDs_;

  std::vector<initElem*> initElems_;
  int initElemCounter_;

  std::vector<sumInElem*> sumInElems_;
  int sumInElemCounter_;

  std::vector<sumInElem*> sumInElemMatrix_;
  int sumInElemMatrixCounter_;

  std::vector<sumInElem*> sumInElemRHS_;
  int sumInElemRHSCounter_;

  double resetSystem_;
  double resetMatrix_;
  double resetRHSVector_;
  double resetInitialGuess_;

  std::vector<nodeBC*> loadNodeBCs_;
  int loadNodeBCsCounter_;

  std::vector<initCR*> initCRMult_;
  int initCRMultCounter_;

  std::vector<loadCR*> loadCRMult_;
  int loadCRMultCounter_;

  std::vector<sharedNodes*> initSharedNodes_;
  int initSharedNodesCounter_;

  std::vector<parameters*> parameters_;
  int parametersCounter_;

  std::vector<setIDLists*> setIDLists_;
  int setIDListsCounter_;

  std::vector<int> setCurrentMatrix_;
  int setCurrentMatrixCounter_;

  std::vector<int> setCurrentRHS_;
  int setCurrentRHSCounter_;

  std::vector<putBlockFieldNodeSolution*> putBlockFieldNodeSolution_;
  int putBlockFieldNodeSolutionCounter_;
};

#endif // _driverData_h_

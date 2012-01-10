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

#ifndef _PoissonData_h_
#define _PoissonData_h_

#include <fei_base.hpp>
//
//This is a class for use in exercising FEI implementations.
//
//This class sets up test data for the Poisson equation on a 2D square,
//and provides query functions for obtaining that data.
//
//The calling program (the 'user' of PoissonData) is left
//with the task of calling the FEI functions.
//
//Also note:
//
// 1. This class only provides 1 element-block per processor currently.
// 2. The function calculateBCs() must be called before the boundary
//    condition data is requested.
//
// Alan Williams 12-20-2000
//

class PoissonData {
  public:
    //constructor -- see PoissonData.cpp for descriptions of these
    //parameters.
    PoissonData(int L,
               int numProcs, int localProc, int outputLevel);

    //destructor.
    ~PoissonData();

    int getElemFormat() {return(elemFormat_); };

    //hardwired for only 1 field...
    int getNumFields() { return(1);};
    int* getFieldSizes() { return(&fieldSize_);};
    int* getFieldIDs() { return(&fieldIDs_[0][0]);};

    GlobalID getElemBlockID() { return(elemBlockID_); };

    int getNumLocalElements() { return(numLocalElements_); };
    GlobalID* getLocalElementIDs() { return(elemIDs_); };
    int getNumNodesPerElement() { return(elem_->numElemNodes()); };

    int* getNumFieldsPerNodeList() { return( numFields_ ); };
    int** getNodalFieldIDsTable() { return( fieldIDs_ ); };

    GlobalID* getElementConnectivity(GlobalID elemID);

    double** getElemStiffness(GlobalID elemID);
    double* getElemLoad(GlobalID elemID);

    void addBCNode(GlobalID nodeID, double x, double y);

    void calculateBCs();

    int getNumBCNodes() { return( BCNodeIDs_.size() ); }
    GlobalID* getBCNodeIDs() { return( &BCNodeIDs_[0] ); }
    int getBCFieldID() { return( fieldIDs_[0][0] ); }
    double* getBCValues() { return( &BCValues_[0] ); }
    

    void getLeftSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs);
    void getRightSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs);
    void getTopSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs);
    void getBottomSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs);
  private:
    void check1();
    void calculateDistribution();

    void messageAbort(const char* message);

    void calculateConnectivity(GlobalID* conn, int size, GlobalID elemID);
    void initializeFieldStuff();
    void deleteFieldArrays();

    void printSharedNodes(const char* str,
			  int numShared,
			  GlobalID* nodeIDs,
                          int** shareProcs,
			  int* numShareProcs);

    Poisson_Elem* elem_; //we're only going to have 1 element instance!!
    int numLocalElements_;
    int startElement_;

    int numProcs_;
    int localProc_;
    int outputLevel_;

    int L_;
    int procX_, procY_;
    int maxProcX_, maxProcY_;

    int numElemBlocks_;
    int solveType_;

    int nodesPerElement_;
    int fieldsPerNode_;
    GlobalID elemBlockID_;
    int elemSetID_;
    int elemFormat_;

    //*************** field description variables *********
    int fieldSize_;
    int* numFields_;
    int** fieldIDs_;
    bool fieldArraysAllocated_;

    //************* element IDs and connectivities ********
    GlobalID* elemIDs_;
    bool elemIDsAllocated_;

    //************* boundary condition stuff **************
    std::vector<GlobalID> BCNodeIDs_;
    std::vector<double> BCValues_;
};

int init_elem_connectivities(FEI* fei, PoissonData& poissonData);

int init_elem_connectivities(fei::MatrixGraph* matrixGraph,
			     PoissonData& poissonData);

int set_shared_nodes(FEI* fei, PoissonData& poissonData);

int set_shared_nodes(fei::VectorSpace* nodeSpace, PoissonData& poissonData);

int load_elem_data(FEI* fei, PoissonData& poissonData);

int load_elem_data_putrhs(FEI* fei, PoissonData& poissonData);

int load_elem_data(fei::MatrixGraph* matrixGraph,
		   fei::Matrix* mat, fei::Vector* rhs,
		   PoissonData& poissonData);

int load_BC_data(FEI* fei, PoissonData& poissonData);

int load_BC_data(fei::LinearSystem* linSys, PoissonData& poissonData);

#endif


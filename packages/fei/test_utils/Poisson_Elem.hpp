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

#ifndef __Poisson_Elem_H
#define __Poisson_Elem_H

class Poisson_Elem {

  public:
    Poisson_Elem();
    ~Poisson_Elem();


    GlobalID getElemID() const {return globalElemID_;};
    void setElemID(GlobalID gNID) {globalElemID_ = gNID;
                                   ID_IsSet_ = true;};

    int numElemRows() const {return numElemRows_;};
    void numElemRows(int gNERows) {numElemRows_ = gNERows;};

    int numElemNodes() const {return numElemNodes_;};
    void numElemNodes(int gNodes) {numElemNodes_ = gNodes;};

    double getElemLength() const {return elemLength_;};
    void setElemLength(double len) {elemLength_ = len;
                                    elemLengthIsSet_ = true;};

    double getTotalLength() const {return totalLength_;};
    void setTotalLength(double len) {totalLength_ = len;
                                    totalLengthIsSet_ = true;};

    int allocateInternals(int DOF);
    int allocateLoad(int DOF);
    int allocateStiffness(int DOF);

    GlobalID* getElemConnPtr(int& size);

    void calculateLoad();
    double* getElemLoad(int& size);

    void calculateStiffness();
    double** getElemStiff(int& size);

    double* getNodalX(int& size) {size = numElemNodes_; return(nodalX_);};
    double* getNodalY(int& size) {size = numElemNodes_; return(nodalY_);};

    void calculateCoords();

    void messageAbort(const char* str);

    void deleteMemory();

//$ temporary output for debugging...

    void dumpToScreen();


  private:
    GlobalID globalElemID_; // global ID number for this element
    bool ID_IsSet_;         // whether ID has been set or not.
    int numElemNodes_;      // number of nodes associated with this element
    int numElemRows_;       // number of rows in the element matrices
    
    GlobalID *nodeList_;      // list of nodes associated with this element
    double* nodalX_;          // list of nodal x-coordinates
    double* nodalY_;          // list of nodal y-coordinates

    double **elemStiff_;      // stiffness matrix for this element
    double *elemLoad_;        // load vector for this element
    bool internalsAllocated_;

    double elemLength_;       // length of a side of this 2D element
    double totalLength_;      // total length of a side of the square region
                              // that this element is in.
    bool elemLengthIsSet_;    // indicates whether length has been set.
    bool totalLengthIsSet_;    // indicates whether length has been set.

    bool loadAllocated_;
    bool stiffAllocated_;
};
 
#endif


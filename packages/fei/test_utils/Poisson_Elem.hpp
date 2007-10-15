#ifndef __Poisson_Elem_H
#define __Poisson_Elem_H

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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


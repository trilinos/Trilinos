/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>
#include <fei_defs.h>
#include <test_utils/Poisson_Elem.hpp>

#include <cmath>
#include <cstdlib>
#include <assert.h>

//==============================================================================
Poisson_Elem::Poisson_Elem() {

    globalElemID_ = (GlobalID)(-1);
    ID_IsSet_ = false;
    numElemNodes_ = 4;
    numElemRows_ = 4;

    nodeList_ = NULL;
    nodalX_ = NULL;
    nodalY_ = NULL;

    elemStiff_ = NULL;
    elemLoad_ = NULL;
    internalsAllocated_ = false;

    elemLength_ = 0.0;
    elemLengthIsSet_ = false;
    totalLength_ = 0.0;
    totalLengthIsSet_ = false;

    loadAllocated_ = false;
    stiffAllocated_ = false;
}

//==============================================================================
Poisson_Elem::~Poisson_Elem() {
   deleteMemory();
}

//==============================================================================
void Poisson_Elem::deleteMemory()
{
    if (!internalsAllocated_) return;

    delete [] nodeList_;
    delete [] nodalX_;
    delete [] nodalY_;

    internalsAllocated_ = false;

    if (loadAllocated_) {
       delete [] elemLoad_;
       loadAllocated_ = false;
    }

    if (stiffAllocated_) {
       delete [] elemStiff_[0];
       delete [] elemStiff_;
       stiffAllocated_ = false;
   }
}

//==============================================================================
int Poisson_Elem::allocateInternals(int DOF) {

    assert(DOF==1);

    nodeList_ = new GlobalID[numElemNodes_];
    if (!nodeList_) return(1);

    nodalX_ = new double[numElemNodes_];
    if (!nodalX_) return(1);

    nodalY_ = new double[numElemNodes_];
    if (!nodalY_) return(1);

    for(int i = 0; i < numElemNodes_; i++) {
        nodeList_[i] = (GlobalID)0;
        nodalX_[i] = 0.0;
        nodalY_[i] = 0.0;
    }

    internalsAllocated_ = true;
    return(0);
}

//==============================================================================
int Poisson_Elem::allocateLoad(int DOF)
{
    assert(DOF==1);

    elemLoad_ = new double[numElemNodes_];
    if (!elemLoad_) return(1);

    loadAllocated_ = true;
    return(0);
}

//==============================================================================
int Poisson_Elem::allocateStiffness(int DOF)
{
    assert(DOF==1);

    elemStiff_ = new double*[numElemNodes_];
    if (!elemStiff_) return(1);

    elemStiff_[0] = NULL;
    elemStiff_[0] = new double[numElemNodes_*numElemNodes_];
    if (!elemStiff_[0]) return(1);

    int i;
    for(i=0; i<numElemNodes_*numElemNodes_; i++) elemStiff_[0][i] = 0.0;

    for(i=1; i<numElemNodes_; i++) {
        elemStiff_[i] = elemStiff_[i-1] + numElemNodes_;
    }

    stiffAllocated_ = true;
    return(0);
}

//==============================================================================
GlobalID* Poisson_Elem::getElemConnPtr(int& size) {

    if (internalsAllocated_) {
        size = numElemNodes_;
        return(nodeList_);
    }
    else {
        size = 0;
        return(NULL);
    }
}

//==============================================================================
double* Poisson_Elem::getElemLoad(int& size) {
        
    if (internalsAllocated_) {
        size = numElemNodes_;
        return(elemLoad_);
    }   
    else {
        size = 0;
        return(NULL);
    }
}

//==============================================================================
double** Poisson_Elem::getElemStiff(int& size) {
        
    if (internalsAllocated_) {
        size = numElemNodes_*numElemNodes_;
        return(elemStiff_);
    }   
    else {
        size = 0;
        return(NULL);
    }
}

//==============================================================================
void Poisson_Elem::calculateCoords() {
//
//This function calculates nodal x- and y-coordinates for this element.
//NOTE: element IDs are assumed to be 1-based.
//
    if (!internalsAllocated_)
        messageAbort("calculateCoords: internals not allocated.");
    if (!elemLengthIsSet_)
        messageAbort("calculateCoords: elemLength not set.");
    if (!totalLengthIsSet_)
        messageAbort("calculateCoords: totalLength not set.");
    if (!ID_IsSet_)
        messageAbort("calculateCoords: elemID not set.");
    if (std::abs(elemLength_) < 1.e-49)
        messageAbort("calculateCoords: elemLength == 0.");

    int lowLeft = 0;
    int lowRight = 1;
    int upperRight = 2;
    int upperLeft = 3;

    int elemsPerSide = (int)std::ceil(totalLength_/elemLength_);

    int elemX = (int)globalElemID_%elemsPerSide;
    if (elemX==0) elemX = elemsPerSide;

    int elemY = ((int)globalElemID_ - elemX)/elemsPerSide + 1;

    //elemX and elemY are 1-based coordinates of this element in
    //the global square of elements. The origin is position (1,1),
    //which is at the bottom left of the square.

    nodalX_[lowLeft] = (elemX-1)*elemLength_;
    nodalX_[upperLeft] = nodalX_[lowLeft];
    nodalX_[lowRight] = elemX*elemLength_;
    nodalX_[upperRight] = nodalX_[lowRight];

    nodalY_[lowLeft] = (elemY-1)*elemLength_;
    nodalY_[lowRight] = nodalY_[lowLeft];
    nodalY_[upperLeft] = elemY*elemLength_;
    nodalY_[upperRight] = nodalY_[upperLeft];
}

//==============================================================================
void Poisson_Elem::messageAbort(const char* str) {
    fei::console_out() << "Poisson_Elem: ERROR: " << str << " Aborting." << FEI_ENDL;
    std::abort();
}

//==============================================================================
void Poisson_Elem::calculateLoad() {

    assert (numElemRows_ == 4);
    
    elemLoad_[0] = -2.0*elemLength_*elemLength_;
    elemLoad_[1] = -2.0*elemLength_*elemLength_;
    elemLoad_[2] = -2.0*elemLength_*elemLength_;
    elemLoad_[3] = -2.0*elemLength_*elemLength_;
}

//==============================================================================
void Poisson_Elem::calculateStiffness() {

    assert (numElemRows_ == 4);
    
    elemStiff_[0][0] = 2.0;
    elemStiff_[0][1] = -1.0;
    elemStiff_[0][2] = 0.0;
    elemStiff_[0][3] = -1.0;

    elemStiff_[1][0] = -1.0;
    elemStiff_[1][1] = 2.0;
    elemStiff_[1][2] = -1.0;
    elemStiff_[1][3] = 0.0;

    elemStiff_[2][0] = 0.0;
    elemStiff_[2][1] = -1.0;
    elemStiff_[2][2] = 2.0;
    elemStiff_[2][3] = -1.0;

    elemStiff_[3][0] = -1.0;
    elemStiff_[3][1] = 0.0;
    elemStiff_[3][2] = -1.0;
    elemStiff_[3][3] = 2.0;
}



//-----------------------------------------------------------------------

void Poisson_Elem::dumpToScreen() {

    int i;
    
    FEI_COUT << " globalElemID_ = " << (int)globalElemID_ << FEI_ENDL;
    FEI_COUT << " numElemRows_  = " << numElemRows_ << FEI_ENDL;
    FEI_COUT << " elemLength_  = " << elemLength_ << FEI_ENDL;
    FEI_COUT << " the " << numElemNodes_ << " nodes: ";
    for (i = 0; i < numElemNodes_; ++i) {
    	FEI_COUT << "   " << (int)nodeList_[i];
    }
    FEI_COUT << FEI_ENDL;
    FEI_COUT << " the " << numElemRows_ << " load vector terms: ";
    for (i = 0; i < numElemRows_; ++i) {
    	FEI_COUT << "   " << elemLoad_[i];
    }
    FEI_COUT << FEI_ENDL << FEI_ENDL;
    return;
}


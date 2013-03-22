/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
//This is a class that will exercise FEI implementations.
//

#include <fei_macros.hpp>
#include <fei_defs.h>
#include <fei_CSVec.hpp>

#include <test_utils/Poisson_Elem.hpp>

#include <test_utils/PoissonData.hpp>

#include <cstdlib>
#include <cmath>

static int int_sqrt(int x) {
//a not-safe function for taking the sqrt of a square int.
    return((int)std::ceil(std::sqrt((double)x)));
}

//==============================================================================
PoissonData::PoissonData(int L,
                       int numProcs, int localProc, int outputLevel)
{
  //
  //PoissonData constructor.
  //
  //Arguments:
  //
  // L:                global square size (number-of-elements along side)
  // numProcs:         number of processors participating in this FEI test.
  // localProc:        local processor number.
  // outputLevel:      affects the amount of screen output.
  //

  L_ = L;
    
  startElement_ = 0;
  numLocalElements_ = 0;

  numProcs_ = numProcs;
  localProc_ = localProc;
  outputLevel_ = outputLevel;

  check1();

  elem_ = new Poisson_Elem();
  int err = elem_->allocateInternals(1);
  err += elem_->allocateLoad(1);
  err += elem_->allocateStiffness(1);
  if (err) messageAbort("Allocation error in element.");

  fieldArraysAllocated_ = false;
  elemIDsAllocated_ = false;

  numFields_ = NULL;
  fieldIDs_ = NULL;

  elemIDs_ = NULL;

  calculateDistribution();

  numElemBlocks_ = 1;
  elemBlockID_ = (GlobalID)0;
  elemSetID_ = 0;
  elemFormat_ = 0;

  nodesPerElement_ = 4;
  fieldsPerNode_ = 1;

  initializeFieldStuff();
}

//==============================================================================
PoissonData::~PoissonData() {
//
//Destructor -- delete any allocated memory.
//

    deleteFieldArrays();

    if (elemIDsAllocated_) delete [] elemIDs_;

    elem_->deleteMemory();
    delete elem_;
}

//==============================================================================
void PoissonData::check1() {
//
//Private function to be called from the constructor, simply makes sure that
//the constructor's input arguments were reasonable.
//
//If they aren't, a message is printed on standard err, and abort() is called.
//
    if (L_ <= 0)                 messageAbort("bar length L <= 0.");
    if (numProcs_ <= 0)          messageAbort("numProcs <= 0.");
    if (L_%int_sqrt(numProcs_)) 
        messageAbort("L must be an integer multiple of sqrt(numProcs).");
    if (localProc_ < 0)          messageAbort("localProc < 0.");
    if (localProc_ >= numProcs_) messageAbort("localProc >= numProcs.");
    if (outputLevel_ < 0)        messageAbort("outputLevel < 0.");
}

//==============================================================================
void PoissonData::calculateDistribution() {
//
//Calculate which elements this processor owns. The element domain is a
//square, and we can assume that sqrt(numProcs_) divides evenly into
//L_. We're working with a (logically) 2D processor arrangement.
//Furthermore, the logical processor layout is such that processor 0 is at
//the bottom left corner of a 2D grid, and a side of the grid is of length
//sqrt(numProcs_). The element domain is numbered such that element 1 is at
//the bottom left corner of the square, and element numbers increase from
//left to right. i.e., element 1 is in position (1,1), element L is in
//position (1,L), element L+1 is in position (2,1).
//
//Use 1-based numbering for the elements and the x- and y- coordinates in
//the element grid, but use 0-based numbering for processor IDs and the
//coordinates in the processor grid.
//
    numLocalElements_ = (L_*L_)/numProcs_;

    elemIDs_ = new GlobalID[numLocalElements_];
    if (!elemIDs_) messageAbort("ERROR allocating elemIDs_.");
    elemIDsAllocated_ = true;

    //0-based x-coordinate of this processor in the 2D processor grid.
    procX_ = localProc_%int_sqrt(numProcs_);

    //0-based maximum processor x-coordinate.
    maxProcX_ = int_sqrt(numProcs_) - 1;

    //0-based y-coordinate of this processor in the 2D processor grid.
    procY_ = localProc_/int_sqrt(numProcs_);

    //0-based maximum processor y-coordinate.
    maxProcY_ = int_sqrt(numProcs_) - 1;

    int sqrtElems = int_sqrt(numLocalElements_);
    int sqrtProcs = int_sqrt(numProcs_);

    //1-based first-element-on-this-processor
    startElement_ = 1 + procY_*sqrtProcs*numLocalElements_ + procX_*sqrtElems;

    if (outputLevel_>1) {
        FEI_COUT << localProc_ << ", calcDist.: numLocalElements: " 
             << numLocalElements_ << ", startElement: " << startElement_ 
             << FEI_ENDL;
        FEI_COUT << localProc_ << ", procX: " << procX_ << ", procY_: " << procY_
             << ", maxProcX: " << maxProcX_ << ", maxProcY: " << maxProcY_
             << FEI_ENDL;
    }

    int offset = 0;
    for(int i=0; i<sqrtElems; i++) {
        for(int j=0; j<sqrtElems; j++) {
            elemIDs_[offset] = (GlobalID)(startElement_ + i*L_ + j);
            offset++;
        }
    }
}

//==============================================================================
void PoissonData::messageAbort(const char* message) {
    fei::console_out() << FEI_ENDL << "PoissonData: " << message 
         << FEI_ENDL << "  Aborting." << FEI_ENDL;
    std::abort();
}

//==============================================================================
GlobalID* PoissonData::getElementConnectivity(GlobalID elemID)
{
   //set the elemID on the internal Poisson_Elem instance.
   elem_->setElemID(elemID);
   elem_->setElemLength(1.0/L_);
   elem_->setTotalLength(1.0);

   //now get a pointer to the element's connectivity array and
   //calculate that connectivity (in place).
   int size = 0;
   GlobalID* elemConn = elem_->getElemConnPtr(size);
   if (size == 0) messageAbort("loadElements: bad conn ptr.");

   calculateConnectivity(elemConn, size, elemID);

   return(elemConn);
}

//==============================================================================
double** PoissonData::getElemStiffness(GlobalID elemID)
{
   elem_->setElemID(elemID);
   elem_->setElemLength(1.0/L_);
   elem_->setTotalLength(1.0);

   //now get a pointer to this element's connectivity array and
   //calculate that connectivity (in place).
   int size = 0;
   GlobalID* elemConn = elem_->getElemConnPtr(size);
   if (size == 0) messageAbort("loadElemStiffnesses: bad conn ptr.");

   calculateConnectivity(elemConn, size, elemID);

   elem_->calculateCoords();

   if (outputLevel_>1) {
      double* x = elem_->getNodalX(size);
      double* y = elem_->getNodalY(size);
      FEI_COUT << localProc_ << ", elemID " << elemID << ", nodes: ";
      for(int j=0; j<size; j++) {
         FEI_COUT << elemConn[j] << " ";
         FEI_COUT << "("<<x[j]<<","<<y[j]<<") ";
      }
      FEI_COUT << FEI_ENDL;
   }

   elem_->calculateStiffness();

   return( elem_->getElemStiff(size) );
}

//==============================================================================
double* PoissonData::getElemLoad(GlobalID elemID)
{
   elem_->setElemID(elemID);
   elem_->setElemLength(1.0/L_);
   elem_->setTotalLength(1.0);

   //now get a pointer to this element's connectivity array and
   //calculate that connectivity (in place).
   int size = 0;
   GlobalID* elemConn = elem_->getElemConnPtr(size);
   if (size == 0) messageAbort("loadElemLoads: bad conn ptr.");

   calculateConnectivity(elemConn, size, elemID);

   elem_->calculateCoords();

   if (outputLevel_>1) {
      double* x = elem_->getNodalX(size);
      double* y = elem_->getNodalY(size);
      FEI_COUT << localProc_ << ", elemID " << elemID << ", nodes: ";
      for(int j=0; j<size; j++) {
         FEI_COUT << elemConn[j] << " ";
         FEI_COUT << "("<<x[j]<<","<<y[j]<<") ";
      }
      FEI_COUT << FEI_ENDL;
   }

   elem_->calculateLoad();

   return( elem_->getElemLoad(size));

}

//==============================================================================
void PoissonData::calculateConnectivity(GlobalID* conn, int size,
                                        GlobalID elemID) {
//
//Calculate a single element's connectivity array -- the list of nodes
//that it 'contains'.
//
//Note that we're assuming the element is a 2D square.
//
    //elemX will be the global 'x-coordinate' of this element in the square. The
    //'origin' is the lower-left corner of the bar, which is element 1,
    //and it is in position 1,1.
    int elemX = (int)elemID%L_;
    if (elemX == 0) elemX = L_;

    //elemY will be the global (1-based) 'y-coordinate'.
    int elemY = ((int)elemID - elemX)/L_ + 1;

    //These are the four nodes for this element.
    GlobalID lowerLeft = elemID + (GlobalID)(elemY-1);
    GlobalID lowerRight = lowerLeft + (GlobalID)1;
    GlobalID upperRight = lowerRight + (GlobalID)(L_+1);
    GlobalID upperLeft = upperRight - (GlobalID)1;

    (void)size;

    //now fill the connectivity array. We'll always fill the connectivity
    //array with the lower left node first, and then proceed counter-clockwise.
    conn[0] = lowerLeft;
    conn[1] = lowerRight;
    conn[2] = upperRight;
    conn[3] = upperLeft;
}

//==============================================================================
void PoissonData::initializeFieldStuff() {
//
//Set up the field-descriptor variables that will be passed
//to the FEI's initFields function, beginInitElemBlock function, etc.
//
//Note we've hardwired 1 dof per field.
//
    fieldSize_ = 1;
    numFields_ = new int[nodesPerElement_];
    fieldIDs_ = new int*[nodesPerElement_];
    for(int i=0; i<nodesPerElement_; i++) {
        numFields_[i] = fieldsPerNode_;
        fieldIDs_[i] = new int[fieldsPerNode_];
        for(int j=0; j<fieldsPerNode_; j++) {
            fieldIDs_[i][j] = fieldsPerNode_;
        }
    }
    fieldArraysAllocated_ = true;
}

//==============================================================================
void PoissonData::deleteFieldArrays() {

    if (fieldArraysAllocated_) {

        for(int i=0; i<nodesPerElement_; i++) {
            delete [] fieldIDs_[i];
        }

        delete [] fieldIDs_;
        delete [] numFields_;
    }
    fieldArraysAllocated_ = false;
}

//==============================================================================
void PoissonData::getLeftSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs) {
//
//This function decides whether any of the nodes along the left edge,
//including the top node but not the bottom node, are shared. It also
//decides which processors the nodes are shared with.
//

    if (numProcs_ == 1) {
        numShared = 0;
        return;
    }

    if (procX_ == 0) {
        //if this proc is on the left edge of the square...

        if (procY_ < maxProcY_) {
            //if this proc is not the top left proc...

            numShared = 1;

            int topLeftElemIndex = numLocalElements_ -
                               int_sqrt(numLocalElements_);

            elem_->setElemID(elemIDs_[topLeftElemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            int size = 0;
            GlobalID* elemConn = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort("loadElements: bad conn ptr.");

            calculateConnectivity(elemConn, size, elemIDs_[topLeftElemIndex]);

            sharedNodeIDs[0] = elemConn[3]; //elem's top left node is node 3
            numProcsPerSharedNode[0] = 2;
            sharingProcs[0][0] = localProc_;
            sharingProcs[0][1] = localProc_ + int_sqrt(numProcs_);

            return;
        }
        else {
            //else this proc is the top left proc...
            numShared = 0;
        }
    }
    else {
        //else this proc is not on the left edge of the square...

        numShared = int_sqrt(numLocalElements_);
        int lowerLeftElemIndex = 0;

        int sqrtElems = int_sqrt(numLocalElements_);

        int shOffset = 0;
        for(int i=0; i<sqrtElems; i++){
            //stride up the left edge of the local elements...
            int size=0;

            int elemIndex = lowerLeftElemIndex+i*sqrtElems;

            elem_->setElemID(elemIDs_[elemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");
      
            calculateConnectivity(nodes, size, elemIDs_[elemIndex]);

            //now put in the top left node
            sharedNodeIDs[shOffset] = nodes[3];
            sharingProcs[shOffset][0] = localProc_-1;
            sharingProcs[shOffset][1] = localProc_;
            numProcsPerSharedNode[shOffset++] = 2;
        }

        if (procY_ < maxProcY_) {
            //if this proc isn't on the top edge, the upper left node (the
            //last one we put into the shared node list) is shared by 4 procs.
            shOffset--;
            numProcsPerSharedNode[shOffset] = 4;
            sharingProcs[shOffset][2] = localProc_ + int_sqrt(numProcs_);
            sharingProcs[shOffset][3] = sharingProcs[shOffset][2] - 1;
        }
    }
}

//==============================================================================
void PoissonData::getRightSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs) {
//
//This function decides whether any of the nodes along the right edge,
//including the bottom node but not the top node, are shared. It also
//decides which processors the nodes are shared with.
//

    if (numProcs_ == 1) {
        numShared = 0;
        return;
    }

    if (procX_ == maxProcX_) {
        //if this proc is on the right edge of the square...

        if (procY_ > 0) {
            //if this proc is not the bottom right proc...

            numShared = 1;

            int lowerRightElemIndex = int_sqrt(numLocalElements_) - 1;

            elem_->setElemID(elemIDs_[lowerRightElemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            int size;
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");
      
            calculateConnectivity(nodes, size, elemIDs_[lowerRightElemIndex]);

            sharedNodeIDs[0] = nodes[1]; //elem's bottom right node is node 1
            numProcsPerSharedNode[0] = 2;
            sharingProcs[0][0] = localProc_;
            sharingProcs[0][1] = localProc_ - int_sqrt(numProcs_);

            return;
        }
        else {
            //else this proc is the bottom right proc...
            numShared = 0;
        }
    }
    else {
        //else this proc is not on the right edge of the square...

        numShared = int_sqrt(numLocalElements_);
        int upperRightElemIndex = numLocalElements_ - 1;

        int sqrtElems = int_sqrt(numLocalElements_);

        int shOffset = 0;
        for(int i=0; i<sqrtElems; i++){
            //stride down the right edge of the local elements...
            int size=0;
            int elemIndex = upperRightElemIndex-i*sqrtElems;
            elem_->setElemID(elemIDs_[elemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");

            calculateConnectivity(nodes, size, elemIDs_[elemIndex]);

            //now put in the lower right node
            sharedNodeIDs[shOffset] = nodes[1];
            sharingProcs[shOffset][0] = localProc_+1;
            sharingProcs[shOffset][1] = localProc_;
            numProcsPerSharedNode[shOffset++] = 2;
        }

        if (procY_ > 0) {
            //if this proc isn't on the bottom edge, the lower right node (the
            //last one we put into the shared node list) is shared by 4 procs.
            shOffset--;
            numProcsPerSharedNode[shOffset] = 4;
            sharingProcs[shOffset][2] = localProc_ - int_sqrt(numProcs_);
            sharingProcs[shOffset][3] = sharingProcs[shOffset][2] + 1;
        }
    }
}

//==============================================================================
void PoissonData::getTopSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs) {
//
//This function decides whether any of the nodes along the top edge,
//including the right node but not the left node, are shared. It also
//decides which processors the nodes are shared with.
//

    if (numProcs_ == 1) {
        numShared = 0;
        return;
    }

    if (procY_ == maxProcY_) {
        //if this proc is on the top edge of the square...

        if (procX_ < maxProcX_) {
            //if this proc is not the top right proc...

            numShared = 1;

            int elemIndex = numLocalElements_ - 1;

            elem_->setElemID(elemIDs_[elemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            int size;
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");

            calculateConnectivity(nodes, size, elemIDs_[elemIndex]);

            sharedNodeIDs[0] = nodes[2]; //elem's top right node is node 2
            numProcsPerSharedNode[0] = 2;
            sharingProcs[0][0] = localProc_;
            sharingProcs[0][1] = localProc_ + 1;

            return;
        }
        else {
            //else this proc is the top right proc...
            numShared = 0;
        }
    }
    else {
        //else this proc is not on the top edge of the square...

        numShared = int_sqrt(numLocalElements_);
        int topLeftElemIndex = numLocalElements_ - int_sqrt(numLocalElements_);

        int sqrtElems = int_sqrt(numLocalElements_);

        int shOffset = 0;
        for(int i=0; i<sqrtElems; i++){
            //stride across the top edge of the local elements...
            int size=0;
            int elemIndex = topLeftElemIndex+i;

            elem_->setElemID(elemIDs_[elemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");

            calculateConnectivity(nodes, size, elemIDs_[elemIndex]);

            //now put in the upper right node
            sharedNodeIDs[shOffset] = nodes[2];
            sharingProcs[shOffset][0] = localProc_+int_sqrt(numProcs_);
            sharingProcs[shOffset][1] = localProc_;
            numProcsPerSharedNode[shOffset++] = 2;
        }
        if (procX_ < maxProcX_) {
            //if this proc isn't on the right edge, the top right node (the
            //last one we put into the shared node list) is shared by 4 procs.
            shOffset--;
            numProcsPerSharedNode[shOffset] = 4;
            sharingProcs[shOffset][2] = localProc_ + 1;
            sharingProcs[shOffset][3] = sharingProcs[shOffset][0] + 1;
        }
    }
}

//==============================================================================
void PoissonData::getBottomSharedNodes(int& numShared, GlobalID* sharedNodeIDs,
                                     int* numProcsPerSharedNode,
                                     int** sharingProcs) {
//
//This function decides whether any of the nodes along the bottom edge,
//including the left node but not the right node, are shared. It also
//decides which processors the nodes are shared with.
//

    if (numProcs_ == 1) {
        numShared = 0;
        return;
    }

    if (procY_ == 0) {
        //if this proc is on the bottom edge of the square...

        if (procX_ > 0) {
            //if this proc is not the bottom left proc...

            numShared = 1;

            int elemIndex = 0;

            elem_->setElemID(elemIDs_[elemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            int size;
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");

            calculateConnectivity(nodes, size, elemIDs_[elemIndex]);

            sharedNodeIDs[0] = nodes[0]; //elem's bottom left node is node 0
            numProcsPerSharedNode[0] = 2;
            sharingProcs[0][0] = localProc_;
            sharingProcs[0][1] = localProc_ - 1;

            return;
        }
        else {
            //else this proc is the top right proc...
            numShared = 0;
        }
    }
    else {
        //else this proc is not on the bottom edge of the square...

        numShared = int_sqrt(numLocalElements_);
        int lowerRightElemIndex = int_sqrt(numLocalElements_) - 1;

        int sqrtElems = int_sqrt(numLocalElements_);

        int shOffset = 0;
        for(int i=0; i<sqrtElems; i++){
            //stride across the bottom edge of the local elements, from 
            //right to left...
            int size=0;
            int elemIndex = lowerRightElemIndex-i;

            elem_->setElemID(elemIDs_[elemIndex]);

            //now get a pointer to this element's connectivity array and
            //calculate that connectivity (in place).
            GlobalID* nodes = elem_->getElemConnPtr(size);
            if (size == 0) messageAbort(": bad conn ptr.");

            calculateConnectivity(nodes, size, elemIDs_[elemIndex]);

            //now put in the lower left node
            sharedNodeIDs[shOffset] = nodes[0];
            sharingProcs[shOffset][0] = localProc_ - int_sqrt(numProcs_);
            sharingProcs[shOffset][1] = localProc_;
            numProcsPerSharedNode[shOffset++] = 2;
        }
        if (procX_ > 0) {
            //if this proc isn't on the left edge, the lower left node (the
            //last one we put into the shared node list) is shared by 4 procs.
            shOffset--;
            numProcsPerSharedNode[shOffset] = 4;
            sharingProcs[shOffset][2] = localProc_ - 1;
            sharingProcs[shOffset][3] = sharingProcs[shOffset][0] - 1;
        }
    }
}

//==============================================================================
void PoissonData::printSharedNodes(const char* str,
				   int numShared, GlobalID* nodeIDs,
                                   int** shareProcs, int* numShareProcs)
{
  for(int i=0; i<numShared; i++) {
    FEI_COUT << localProc_ << ", " << str << " node: " << (int) nodeIDs[i];
    FEI_COUT << ", procs: ";
    for(int j=0; j<numShareProcs[i]; j++) {
      FEI_COUT << shareProcs[i][j] << " ";
    }
    FEI_COUT << FEI_ENDL;
  }
}

//==============================================================================
void PoissonData::calculateBCs() {
//
//This function figures out which nodes lie on the boundary. The ones that
//do are added to the BC set, along with appropriate alpha/beta/gamma values.
//
    for(int i=0; i<numLocalElements_; i++) {
       elem_->setElemID(elemIDs_[i]);
       elem_->setElemLength(1.0/L_);
       elem_->setTotalLength(1.0);

       //now get a pointer to this element's connectivity array and
       //calculate that connectivity (in place).
       int size = 0;
       GlobalID* nodeIDs = elem_->getElemConnPtr(size);
       if (size == 0) messageAbort("loadElements: bad conn ptr.");

       calculateConnectivity(nodeIDs, size, elemIDs_[i]);

       elem_->calculateCoords();

       double* xcoord = elem_->getNodalX(size);
       double* ycoord = elem_->getNodalY(size);

       //now loop over the nodes and see if any are on a boundary.
       for(int j=0; j<size; j++) {
          if ((std::abs(xcoord[j]) < 1.e-49) || (std::abs(xcoord[j] - 1.0) < 1.e-49) ||
             (std::abs(ycoord[j]) < 1.e-49) || (std::abs(ycoord[j] - 1.0) < 1.e-49)) {

             addBCNode(nodeIDs[j], xcoord[j], ycoord[j]);
          }
       }
    }
}

//==============================================================================
void PoissonData::addBCNode(GlobalID nodeID, double x, double y){

  std::vector<GlobalID>::iterator
    iter = std::lower_bound(BCNodeIDs_.begin(), BCNodeIDs_.end(), nodeID);

  if (iter == BCNodeIDs_.end() || *iter != nodeID) {
    unsigned offset = iter - BCNodeIDs_.begin();
    BCNodeIDs_.insert(iter, nodeID);

    double bcValue = std::pow(x, 2.0) + std::pow(y, 2.0);

    BCValues_.insert(BCValues_.begin()+offset, bcValue);
  }
}

//==============================================================================
int init_elem_connectivities(FEI* fei, PoissonData& poissonData)
{
  //first load the information that defines this element block, and
  //the topology of each element in this element block.

  GlobalID elemBlockID = poissonData.getElemBlockID();
  int numLocalElements = poissonData.getNumLocalElements();
  int numNodesPerElement = poissonData.getNumNodesPerElement();
  int* numFieldsPerNode = poissonData.getNumFieldsPerNodeList();
  int** fieldIDsTable = poissonData.getNodalFieldIDsTable();

  CHK_ERR( fei->initElemBlock(elemBlockID,
			      numLocalElements,
			      numNodesPerElement,
			      numFieldsPerNode,
			      fieldIDsTable,
			      0, // no element-centered degrees-of-freedom
			      NULL, //null list of elem-dof fieldIDs
			      FEI_NODE_MAJOR) );

  //now let's loop over all of the local elements, giving their 
  //nodal connectivity lists to the FEI.

  GlobalID* elemIDs = poissonData.getLocalElementIDs();

  for(int elem=0; elem<numLocalElements; elem++) {
    GlobalID* elemConnectivity =
      poissonData.getElementConnectivity(elemIDs[elem]);

    CHK_ERR( fei->initElem(elemBlockID, elemIDs[elem], elemConnectivity) );
  }

  return(0);
}

//==============================================================================
int set_shared_nodes(FEI* fei, PoissonData& poissonData)
{
   int numLocalElements = poissonData.getNumLocalElements();
   int maxNumSharedNodes = (int)std::sqrt((double)numLocalElements);
   GlobalID* sharedNodeIDs = new GlobalID[maxNumSharedNodes];
   int* numProcsPerSharedNode = new int[maxNumSharedNodes];
   int** sharingProcs = new int*[maxNumSharedNodes];
   for(int i=0; i<maxNumSharedNodes; i++) sharingProcs[i] = new int[4];

   int numShared;

   //first, get the shared-node data for the left edge of the local block

   poissonData.getLeftSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( fei->initSharedNodes(numShared, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   //now, get the shared-node data for the right edge of the local block

   poissonData.getRightSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( fei->initSharedNodes(numShared, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   //now, get the shared-node data for the bottom edge of the local block

   poissonData.getBottomSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( fei->initSharedNodes(numShared, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   //finally, get the shared-node data for the top edge of the local block

   poissonData.getTopSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( fei->initSharedNodes(numShared, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   for(int j=0; j<maxNumSharedNodes; j++) delete [] sharingProcs[j];
   delete [] sharingProcs;
   delete [] numProcsPerSharedNode;
   delete [] sharedNodeIDs;

   return(0);
}

//==============================================================================
int load_elem_data(FEI* fei, PoissonData& poissonData)
{
  GlobalID elemBlockID = poissonData.getElemBlockID();
  int numLocalElements = poissonData.getNumLocalElements();
  GlobalID* elemIDs = poissonData.getLocalElementIDs();

  for(int elem=0; elem<numLocalElements; elem++) {
    GlobalID* elemConnectivity =
      poissonData.getElementConnectivity(elemIDs[elem]);
    double** elemStiffness = poissonData.getElemStiffness(elemIDs[elem]);

    CHK_ERR( fei->sumInElemMatrix(elemBlockID, elemIDs[elem],
				  elemConnectivity, elemStiffness,
				  poissonData.getElemFormat()));

    double* elemLoad = poissonData.getElemLoad(elemIDs[elem]);

    CHK_ERR( fei->sumInElemRHS(elemBlockID, elemIDs[elem],
			       elemConnectivity, elemLoad));
  }

  return(0);
}

//==============================================================================
int load_elem_data_putrhs(FEI* fei, PoissonData& poissonData)
{
  GlobalID elemBlockID = poissonData.getElemBlockID();
  int numLocalElements = poissonData.getNumLocalElements();
  GlobalID* elemIDs = poissonData.getLocalElementIDs();

  int numIDs = poissonData.getNumNodesPerElement();

  int* fieldID = poissonData.getFieldIDs();

  fei::CSVec rhs;

  for(int elem=0; elem<numLocalElements; elem++) {
    GlobalID* elemConnectivity =
      poissonData.getElementConnectivity(elemIDs[elem]);
    double** elemStiffness = poissonData.getElemStiffness(elemIDs[elem]);

    CHK_ERR( fei->sumInElemMatrix(elemBlockID, elemIDs[elem],
                                  elemConnectivity, elemStiffness,
                                  poissonData.getElemFormat()));

    double* elemLoad = poissonData.getElemLoad(elemIDs[elem]);

    for(int i=0; i<numIDs; ++i) {
      fei::add_entry(rhs, elemConnectivity[i], elemLoad[i]);
    }
  }

  fei->loadComplete();

  fei->putIntoRHS(0, *fieldID, rhs.size(),
                  &(rhs.indices()[0]), &(rhs.coefs()[0]));

  return(0);
}

//==============================================================================
int load_BC_data(FEI* fei, PoissonData& poissonData)
{
  //first, have the data object generate the BC data
  poissonData.calculateBCs();

  int numBCNodes = poissonData.getNumBCNodes();
  GlobalID* nodeIDs = poissonData.getBCNodeIDs();
  int fieldID = poissonData.getBCFieldID();
  double* values = poissonData.getBCValues();

  std::vector<int> offsets(numBCNodes, 0);

  CHK_ERR( fei->loadNodeBCs(numBCNodes, nodeIDs, fieldID,
			    &offsets[0], values) );

  return(0);
}

//==============================================================================
int init_elem_connectivities(fei::MatrixGraph* matrixGraph,
			     PoissonData& poissonData)
{
  //first load the information that defines this element block, and
  //the topology of each element in this element block.

  GlobalID elemBlockID = poissonData.getElemBlockID();
  int numLocalElements = poissonData.getNumLocalElements();
  int numNodesPerElement = poissonData.getNumNodesPerElement();
  int** fieldIDsTable = poissonData.getNodalFieldIDsTable();

  int nodeIDType = 0;

  int patternID =
    matrixGraph->definePattern(numNodesPerElement,
			     nodeIDType, fieldIDsTable[0][0]);

  CHK_ERR( matrixGraph->initConnectivityBlock(elemBlockID,
					    numLocalElements, patternID) );

  //now let's loop over all of the local elements, giving their 
  //nodal connectivity lists to the matrixGraph object.

  GlobalID* elemIDs = poissonData.getLocalElementIDs();

  for(int elem=0; elem<numLocalElements; elem++) {
    GlobalID* elemConnectivity =
      poissonData.getElementConnectivity(elemIDs[elem]);

    CHK_ERR( matrixGraph->initConnectivity(elemBlockID, elemIDs[elem],
					 elemConnectivity) );
  }

  return(0);
}

//==============================================================================
int set_shared_nodes(fei::VectorSpace* nodeSpace, PoissonData& poissonData)
{
   int numLocalElements = poissonData.getNumLocalElements();
   int maxNumSharedNodes = (int)std::sqrt((double)numLocalElements);
   GlobalID* sharedNodeIDs = new GlobalID[maxNumSharedNodes];
   int* numProcsPerSharedNode = new int[maxNumSharedNodes];
   int** sharingProcs = new int*[maxNumSharedNodes];
   for(int i=0; i<maxNumSharedNodes; i++) sharingProcs[i] = new int[4];

   int numShared;

   //first, get the shared-node data for the left edge of the local block

   poissonData.getLeftSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);
   int nodeIDType = 0;

   CHK_ERR( nodeSpace->initSharedIDs(numShared, nodeIDType, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   //now, get the shared-node data for the right edge of the local block

   poissonData.getRightSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( nodeSpace->initSharedIDs(numShared, nodeIDType, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   //now, get the shared-node data for the bottom edge of the local block

   poissonData.getBottomSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( nodeSpace->initSharedIDs(numShared, nodeIDType, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   //finally, get the shared-node data for the top edge of the local block

   poissonData.getTopSharedNodes(numShared, sharedNodeIDs,
                                  numProcsPerSharedNode, sharingProcs);

   CHK_ERR( nodeSpace->initSharedIDs(numShared, nodeIDType, sharedNodeIDs,
                                numProcsPerSharedNode, sharingProcs));

   for(int j=0; j<maxNumSharedNodes; j++) delete [] sharingProcs[j];
   delete [] sharingProcs;
   delete [] numProcsPerSharedNode;
   delete [] sharedNodeIDs;

   return(0);
}

//==============================================================================
int load_elem_data(fei::MatrixGraph* matrixGraph,
		   fei::Matrix* mat, fei::Vector* rhs,
		   PoissonData& poissonData)
{
  GlobalID elemBlockID = poissonData.getElemBlockID();
  int numLocalElements = poissonData.getNumLocalElements();
  GlobalID* elemIDs = poissonData.getLocalElementIDs();

  int numIndices = matrixGraph->getConnectivityNumIndices(elemBlockID);

  std::vector<int> indicesArray(numIndices);
  int* indicesPtr = &indicesArray[0];

  for(int elem=0; elem<numLocalElements; elem++) {
    double** elemStiffness = poissonData.getElemStiffness(elemIDs[elem]);

    int checkNumIndices = 0;
    CHK_ERR( matrixGraph->getConnectivityIndices(elemBlockID, elemIDs[elem],
					       numIndices, indicesPtr,
					       checkNumIndices) );
    if (checkNumIndices != numIndices) return(-1);

    CHK_ERR( mat->sumIn(elemBlockID, elemIDs[elem],
			elemStiffness));

    double* elemLoad = poissonData.getElemLoad(elemIDs[elem]);

    CHK_ERR( rhs->sumIn(numIndices, indicesPtr, elemLoad));
  }

  return(0);
}

//==============================================================================
int load_BC_data(fei::LinearSystem* linSys, PoissonData& poissonData)
{
  //first, have the data object generate the BC data
  poissonData.calculateBCs();

  int numBCNodes = poissonData.getNumBCNodes();
  GlobalID* nodeIDs = poissonData.getBCNodeIDs();
  int fieldID = poissonData.getBCFieldID();
  double* values = poissonData.getBCValues();

  CHK_ERR( linSys->loadEssentialBCs(numBCNodes, nodeIDs, 0, fieldID, 0, values) );

  return(0);
}

/*Paul
01-August-2002 Initial writeup.
05-August-2002 switched to images.
06-August-2002 Completed switch.
*/

namespace Tpetra {

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType elementSize) 
  : Object("Tpetra::BlockElementSpace") 
  ,constantSize_(true) 
  ,elementSize_(elementSize) 
  ,minMySize_(elementSize) 
  ,maxMySize_(elementSize) 
  ,minGlobalSize_(elementSize) 
  ,maxGlobalSize_(elementSize) 
  ,elementSizeList_(0)
  ,pointToElementList_(0)
  ,firstPointList_(0)
  ,ElementSpace_(&ElementSpace)
{
  if(elementSize <= 0)
    throw reportError("elementSize = " + toString(elementSize) + ".  Should be > 0.", -1);

  numGlobalPoints_ = elementSize_ * elementSpace().getNumGlobalElements();
  numMyPoints_ = numGlobalPoints_;
  OrdinalType nME = elementSpace().getNumMyElements();
  if(nME > 0) {
		elementSizeList_ = new OrdinalType[nME];
		for(OrdinalType i = 0; i < nME; i++)
			elementSizeList_[i] = elementSize_;
	}
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType* elementSizeList) 
  : Object("Tpetra::BlockElementSpace") 
  ,constantSize_(false) 
  ,elementSize_(0) 
  ,elementSizeList_(0)
  ,pointToElementList_(0)
  ,firstPointList_(0)
  ,ElementSpace_(&ElementSpace)
{
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 0; i < nME; i++)
    if(elementSizeList[i] <= 0)
      throw reportError("An element in elementSizeList = " + toString(elementSizeList[i]) + ".  Should be > 0.", -1);
  if(nME > 0) {
    elementSizeList_ = new OrdinalType[nME];
    minMySize_ = elementSizeList[0];
    maxMySize_ = elementSizeList[0];
    numMyPoints = 0;
    for(OrdinalType i = 0; i < nME; i++) {
      elementSizeList_[i] = elementSizeList[i];
      minMySize_ = TPETRA_MIN(minMySize_, elementSizeList[i]);
      maxMySize_ = TPETRA_MAX(maxMySize_, elementSizeList[i]);
      numMyPoints_ += elementSizeList[i];
    }
  }
  else {
    minMySize_ = 1;
    maxMySize_ = 1;
    numMyPoints_ = 0;
  }
  if(elementSpace().isGlobal() == false) {
    minGlobalSize_ = minMySize_;
    maxGlobalSize_ = maxMySize_;
  }
  else {
    elementSpace().comm().sumAll(&numMyPoints_, &numGlobalPoints_, 1);
    elementSpace().comm().minAll(&minMySize_, &minGlobalSize_, 1);
    elementSpace().comm().maxAll(&maxMySize_, &maxGlobalSize_, 1);
  }
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(BlockElementSpace<OrdinalType>& BlockElementSpace) 
  : Object(BlockElementSpace.label())
  ,constantSize_(BlockElementSpace.constantSize_)
  ,elementSize_(BlockElementSpace.elementSize_)
  ,numMyPoints_(BlockElementSpace.numMyPoints_)
  ,numGlobalPoints_(BlockElementSpace.numGlobalPoints_)
  ,minMySize_(BlockElementSpace.minMySize_)
  ,maxMySize_(BlockElementSpace.maxMySize_)
  ,minGlobalSize_(BlockElementSpace.minGlobalSize_)
  ,maxGlobalSize_(BlockElementSpace.maxGlobalSize_)
  ,elementSizeList_(0)
  ,pointToElementList_(0)
  ,firstPointList_(0)
  ,ElementSpace_(BlockElementSpace.ElementSpace_)
{
  // If eSL, pTEL, or fPL were non-null in BlockElementSpace, create them here
  if(BlockElementSpace.elementSizeList_ != 0)
    OrdinalType* tmp = getElementSizeList();
  if(BlockElementSpace.pointToElementList_ != 0)
    OrdinalType* tmp = getPointToElementList();
  if(BlockElementSpace.firstPointList_ != 0)
    OrdinalType* tmp = getFirstPointInElementList();
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::~BlockElementSpace() {
  if(elementSizeList_ != 0) {
    delete [] elementSizeList_;
    elementSizeList_ = 0;
  }
  if(pointToElementList_ != 0) {
    delete [] pointToElementList_;
    pointToElementList_ = 0;
  }
  if(firstPointList_ != 0) {
    delete [] firstPointList_;
    firstPointList_ = 0;
  }
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getRemoteIDList(OrdinalType numIDs, const OrdinalType* GIDList, OrdinalType* imageIDList, OrdinalType* LIDList, OrdinalType* elementSizeList) const {
  elementSpace().getRemoteIDList(numIDs, GIDList, imageIDList, LIDList);

  if(isConstantElementSize())
    for(OrdinalType i = 0; i < numIDs; i++)
      elementSizeList[i] = getElementSize();
  else
    for(OrdinalType i = 0; i < numIDs; i++)
      elementSizeList[i] = getElementSize(LIDList[i]);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getLocalElementID(OrdinalType pointID, OrdinalType& elementID, OrdinalType& elementOffset) const {
  pointID -= elementSpace().getIndexBase(); // convert from indexBase-based to zero-based counting.
  if(pointID < 0 || pointID > getNumMyPoints())
    throw reportError("PointID " + toString(pointID) + " was not found on this processor.", 1);
  if(isConstantSize()) {
    elementID = pointID / getElementSize();
    elementOffset = pointID % getElementSize();
  }
  else {
    OrdinalType* tmp = getPointToElementList();
    elementID = tmp[pointID];
    tmp = getFirstPointInElementList();
    elementOffset = pointID - tmp[elementID];
    tmp = 0;
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType BlockElementSpace<OrdinalType>::getElementSize() const {
  if(!isConstantElementSize()) 
    throw reportError("This BlockElementSpace does not have a constant element size.", 3);
  return(elementSize_);
}

//=======================================================================
template<typename OrdinalType>
OrdinalType BlockElementSpace<OrdinalType>::getElementSize(OrdinalType LID) const {
  if(elementSpace().isMyLID(LID) == false)
    throw reportError("Local ID " + toString(LID) + " was not found on this processor.", 2);
  if(isConstantSize())
    return(getElementSize());
  else {
    LID -= elementSpace().getIndexBase(); // convert to zero-based counting.
    return(elementSizeList_[LID]);
  }
}

//=======================================================================
template<typename OrdinalType>
bool BlockElementSpace<OrdinalType>::isSameAs(const BlockElementSpace<OrdinalType>& BlockElementSpace) const {
  if(this == &BlockElementSpace)
    return(true);
  // check to make sure ElementSpaces match.
  if(! elementSpace().isSameAs(BlockElementSpace.elementSpace()))
    return(false);
  if(isConstantElementSize() != BlockElementSpace.isConstantElementSize())
    return(false);

  if(isConstantElementSize()) {
    if(getElementSize() != BlockElementSpace.getElementSize())
      return(false);
    else
      return(true);
  }
  else {
    int mySameBES = 1;
    OrdinalType nME = elementSpace().getNumMyElements();
    for(OrdinalType i = 0; i < nME; i++)
      if(elementSizeList_[i] != BlockElementSpace.elementSizeList_[i])
				mySameBES = 0;
    
    int globalSameBES = 0;
    elementSpace().comm().minAll(&mySameBES, &globalSameBES, 1);
    return(globalSameBES == 1);
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType* BlockElementSpace<OrdinalType>::getFirstPointInElementList() const {
  OrdinalType nME = elementSpace().getNumMyElements();
	if((firstPointList_ == 0) && (nME > 0)) {
    firstPointList_ = new OrdinalType[nME];
    getFirstPointInElementList(firstPointList_);
  }
  return(firstPointList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getFirstPointInElementList(OrdinalType* firstPointInElementList) const {
  OrdinalType iB = elementSpace().getIndexBase();
	firstPointInElementList[0] = iB;
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 1; i < nME; i++)
    firstPointInElementList[i] = firstPointInElementList[i - 1] + elementSizeList_[i - 1];
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getElementSizeList(OrdinalType* elementSizeList) const {
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 0; i < nME; i++)
    elementSizeList[i] = elementSizeList_[i];
}

//=======================================================================
template<typename OrdinalType>
OrdinalType* BlockElementSpace<OrdinalType>::getPointToElementList() const {
	OrdinalType nME = elementSpace().getNumMyElements();
  if((pointToElementList_ == 0) && (nME > 0)) {
    pointToElementList_ = new OrdinalType[nME];
    getPointToElementList(pointToElementList_);
  }
  return(pointToElementList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getPointToElementList(OrdinalType* pointToElementList) const {
  OrdinalType currPos = 0;
  OrdinalType nME = elementSpace().getNumMyElements();
  OrdinalType currLID = elementSpace().getIndexBase();
  OrdinalType currSize;
  for(OrdinalType i = 0; i < nME; i++) {
    currSize = elementSizeList_[i];
    for(OrdinalType j = 0; j < currSize; j++) {
      pointToElementList[currPos] = currLID;
      currPos++;
    }
    currLID++;
  }
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::print(ostream& os) const {
  OrdinalType* elementSizeList1 = getElementSizeList();
  OrdinalType* firstPointList1 = getFirstPointInElementList();
  OrdinalType* pointToElementList1 = getPointToElementList();

  int myImageID = elementSpace().comm().getMyImageID();
  int numImages = elementSpace().comm().getNumImages();

  for(int imageCtr = 0; imageCtr < numImages; imageCtr++) {
    if(myImageID == imageCtr) {
      if(myImageID == 0) {
				os << "\nNumber of Global Points = "; os << getNumGlobalPoints(); os << endl;
				if(isConstantElementSize()) {
					os << "Constant Element Size   = "; os << getElementSize(); os << endl;
				}
				else {
					os << "Global Min Element Size = "; os << getMinElementSize(); os << endl;
					os << "Global Max Element Size = "; os << getMaxElementSize(); os << endl;
				}
      }
      os << endl;
			
      os <<     "Number of Local Points  = "; os << getNumMyPoints(); os << endl;
      if(!isConstantElementSize()) {
				os <<     "Min Local Element Size  = "; os << getMinMyElementSize(); os << endl;
				os <<     "Max Local Element Size  = "; os << getMaxMyElementSize(); os << endl;
      }
			
      os << "\nElementSizeList"; os << endl;
      OrdinalType nME = elementSpace().getNumMyElements();
      for(OrdinalType i = 0; i < nME; i++)
				os << elementSizeList1[i] << " ";
      os << endl;
      os << "\nFirstPointInElementList"; os << endl;
      for(OrdinalType i = 0; i < nME; i++)
				os << firstPointList1[i] << " ";
      os << endl;
      os << "\nPointToElementList"; os << endl;
      for(OrdinalType i = 0; i < numMyPoints_; i++)
				os << pointToElementList1[i] << " ";
      os << endl;
      os << flush;
    }
  }
	//cout << "BES(print): calling gMGE() from tmp pointer" << endl;
	//OrdinalType* tmp = elementSpace().getMyGlobalElements();
	//cout << "BES(print): printing ES" << endl;
  elementSpace().print(os);
  // global ops to let I/O complete done by elementSpace::print.
}

} // Tpetra namespace
//=======================================================================

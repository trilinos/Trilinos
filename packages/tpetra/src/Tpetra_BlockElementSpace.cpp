/*Paul
01-August-2002 Initial writeup.
05-August-2002 switched to images.
06-August-2002 Completed switch.
21-Sept-2002 Platform/Comm switch done.
16-Oct-2002 Updated to use BESData 
*/

namespace Tpetra {

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType elementSize) 
  : Object("Tpetra::BlockElementSpace") 
	, BlockElementSpaceData_(new BlockElementSpaceData<OrdinalType>(ElementSpace, true))
{
	// get ES data
  OrdinalType numMyElements = ElementSpace.getNumMyElements();
	OrdinalType numGlobalElements = ElementSpace.getNumGlobalElements();

	// initial throws
  if(elementSize <= 0)
    throw reportError("elementSize = " + toString(elementSize) + ".  Should be > 0.", -1);

	// initialize elementSizeList
  if(numMyElements > 0) {
		BlockElementSpaceData_->elementSizeList_ = new OrdinalType[numMyElements];
		for(OrdinalType i = 0; i < numMyElements; i++)
			BlockElementSpaceData_->elementSizeList_[i] = elementSize;
	}

	// set sizes
  BlockElementSpaceData_->elementSize_ = elementSize;
  BlockElementSpaceData_->minMySize_ = elementSize;
  BlockElementSpaceData_->maxMySize_ = elementSize;
  BlockElementSpaceData_->minGlobalSize_ = elementSize;
  BlockElementSpaceData_->maxGlobalSize_ = elementSize;

	// compute numGlobalPoints & numMyPoints
  BlockElementSpaceData_->numGlobalPoints_ = elementSize * numGlobalElements;
  BlockElementSpaceData_->numMyPoints_ = elementSize * numMyElements;
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType* elementSizeList) 
  : Object("Tpetra::BlockElementSpace")
	, BlockElementSpaceData_(new BlockElementSpaceData<OrdinalType>(ElementSpace, false))
{
	// get ES data
  OrdinalType numMyElements = ElementSpace.getNumMyElements();

	// initial throws
  for(OrdinalType i = 0; i < numMyElements; i++)
    if(elementSizeList[i] <= 0)
      throw reportError("An element in elementSizeList = " + toString(elementSizeList[i]) + ".  Should be > 0.", -1);

	// initialize elementSizeList and compute minMySize, MaxMySize, & numMyPoints
	if(numMyElements > 0) {
    BlockElementSpaceData_->elementSizeList_ = new OrdinalType[numMyElements];
    OrdinalType minMySize = elementSizeList[0];
    OrdinalType maxMySize = elementSizeList[0];
    OrdinalType numMyPoints = 0;
    for(OrdinalType i = 0; i < numMyElements; i++) {
      BlockElementSpaceData_->elementSizeList_[i] = elementSizeList[i];
      minMySize = TPETRA_MIN(minMySize, elementSizeList[i]);
      maxMySize = TPETRA_MAX(maxMySize, elementSizeList[i]);
      numMyPoints += elementSizeList[i];
		}
		BlockElementSpaceData_->minMySize_ = minMySize;
		BlockElementSpaceData_->maxMySize_ = maxMySize;
		BlockElementSpaceData_->numMyPoints_ = numMyPoints;
	}
  else {
    BlockElementSpaceData_->minMySize_ = 1;
    BlockElementSpaceData_->maxMySize_ = 1;
    BlockElementSpaceData_->numMyPoints_ = 0;
  }

	// compute minGlobalSize & maxGlobalSize
  if(elementSpace().isGlobal() == false) {
    BlockElementSpaceData_->minGlobalSize_ = BlockElementSpaceData_->minMySize_;
    BlockElementSpaceData_->maxGlobalSize_ = BlockElementSpaceData_->maxMySize_;
  }
  else {
    elementSpace().comm().sumAll(&BlockElementSpaceData_->numMyPoints_, &BlockElementSpaceData_->numGlobalPoints_, 1);
    elementSpace().comm().minAll(&BlockElementSpaceData_->minMySize_, &BlockElementSpaceData_->minGlobalSize_, 1);
    elementSpace().comm().maxAll(&BlockElementSpaceData_->maxMySize_, &BlockElementSpaceData_->maxGlobalSize_, 1);
  }
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(BlockElementSpace<OrdinalType>& BlockElementSpace) 
  : Object(BlockElementSpace.label())
	, BlockElementSpaceData_(BlockElementSpace.BlockElementSpaceData_) {}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getRemoteIDList(const OrdinalType numIDs, const OrdinalType* GIDList, OrdinalType* imageIDList, 
																										 OrdinalType* LIDList, OrdinalType* elementSizeList) const {
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
  return(BlockElementSpaceData_->elementSize_);
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
    return(BlockElementSpaceData_->elementSizeList_[LID]);
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
      if(BlockElementSpaceData_->elementSizeList_[i] != BlockElementSpace.BlockElementSpaceData_->elementSizeList_[i])
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
	if((BlockElementSpaceData_->firstPointList_ == 0) && (nME > 0)) {
    BlockElementSpaceData_->firstPointList_ = new OrdinalType[nME];
    getFirstPointInElementList(BlockElementSpaceData_->firstPointList_);
  }
  return(BlockElementSpaceData_->firstPointList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getFirstPointInElementList(OrdinalType* firstPointInElementList) const {
  OrdinalType iB = elementSpace().getIndexBase();
	firstPointInElementList[0] = iB;
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 1; i < nME; i++)
    firstPointInElementList[i] = firstPointInElementList[i - 1] + BlockElementSpaceData_->elementSizeList_[i - 1];
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getElementSizeList(OrdinalType* elementSizeList) const {
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 0; i < nME; i++)
    elementSizeList[i] = BlockElementSpaceData_->elementSizeList_[i];
}

//=======================================================================
template<typename OrdinalType>
OrdinalType* BlockElementSpace<OrdinalType>::getPointToElementList() const {
	OrdinalType numPoints = getNumMyPoints();
  if((BlockElementSpaceData_->pointToElementList_ == 0) && (numPoints > 0)) {
    BlockElementSpaceData_->pointToElementList_ = new OrdinalType[numPoints];
    getPointToElementList(BlockElementSpaceData_->pointToElementList_);
  }
  return(BlockElementSpaceData_->pointToElementList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getPointToElementList(OrdinalType* pointToElementList) const {
  OrdinalType currPos = 0;
  OrdinalType nME = elementSpace().getNumMyElements();
  OrdinalType currLID = elementSpace().getIndexBase();
  OrdinalType currSize;
  for(OrdinalType i = 0; i < nME; i++) {
    currSize = BlockElementSpaceData_->elementSizeList_[i];
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
 
  int myImageID = elementSpace().platform().getMyImageID();
  int numImages = elementSpace().platform().getNumImages();

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
      for(OrdinalType i = 0; i < BlockElementSpaceData_->numMyPoints_; i++)
				os << pointToElementList1[i] << " ";
			os << endl;
      os << flush;
    }
	}

  elementSpace().print(os);

  // global ops to let I/O complete done by elementSpace::print.
}

} // Tpetra namespace
//=======================================================================

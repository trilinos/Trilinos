/*Paul
01-August-2002 Initial writeup.
05-August-2002 switched to images.
06-August-2002 Completed switch.
21-Sept-2002 Platform/Comm switch done.
16-Oct-2002 Updated to use BESData
12-Nov-2002 Updated to use createOrdinalComm() instead of createComm() (nothing changed)
24-Nov-2002 Updated for imageID methods moved back to Comm. Changed to use massive BESData constructor calls.
27-Jan-2003 Updated for .hpp and for new const syntax.
06-Feb-2003 Tweaked a lil, retested.
*/

namespace Tpetra {

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType elementSize) 
  : Object("Tpetra::BlockElementSpace") 
	, BlockElementSpaceData_()
{
	// get ES data
  OrdinalType numMyElements = ElementSpace.getNumMyElements();
	OrdinalType numGlobalElements = ElementSpace.getNumGlobalElements();

	// initial throws
  if(elementSize <= 0)
    throw reportError("elementSize = " + toString(elementSize) + ".  Should be > 0.", -1);

	// initialize elementSizeList
	OrdinalType* elementSizeList = 0;
  if(numMyElements > 0) {
		elementSizeList = new OrdinalType[numMyElements];
		for(OrdinalType i = 0; i < numMyElements; i++)
			elementSizeList[i] = elementSize;
	}

	// set sizes
  OrdinalType minMySize = elementSize;
  OrdinalType maxMySize = elementSize;
  OrdinalType minGlobalSize = elementSize;
  OrdinalType maxGlobalSize = elementSize;

	// compute numGlobalPoints & numMyPoints
  OrdinalType numGlobalPoints = elementSize * numGlobalElements;
  OrdinalType numMyPoints = elementSize * numMyElements;

	// call BESData constructor
	BlockElementSpaceData_.reset(new BlockElementSpaceData<OrdinalType>(ElementSpace, true, elementSize, numMyPoints, 
																																			numGlobalPoints, minMySize, maxMySize, 
																																			minGlobalSize, maxGlobalSize, elementSizeList));
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(ElementSpace<OrdinalType>& ElementSpace, OrdinalType* elementSizeList) 
  : Object("Tpetra::BlockElementSpace")
	, BlockElementSpaceData_()
{
	// get ES data
  OrdinalType numMyElements = ElementSpace.getNumMyElements();

	// initial throws
  for(OrdinalType i = 0; i < numMyElements; i++)
    if(elementSizeList[i] <= 0)
      throw reportError("An element in elementSizeList = " + toString(elementSizeList[i]) + ".  Should be > 0.", -1);

	// initialize elementSizeList and compute minMySize, MaxMySize, & numMyPoints
	//   we copy elementSizeList into our own array because elementSizeList (the user's array)
	//   is not guaranteed to always hold the same values it does now.
	OrdinalType* myElementSizeList = 0;
	OrdinalType minMySize = 1;
	OrdinalType maxMySize = 1;
	OrdinalType numMyPoints = 0;
	if(numMyElements > 0) {
    myElementSizeList = new OrdinalType[numMyElements];
    minMySize = elementSizeList[0];
    maxMySize = elementSizeList[0];
    numMyPoints = 0;
    for(OrdinalType i = 0; i < numMyElements; i++) {
      myElementSizeList[i] = elementSizeList[i];
      minMySize = TPETRA_MIN(minMySize, elementSizeList[i]);
      maxMySize = TPETRA_MAX(maxMySize, elementSizeList[i]);
      numMyPoints += elementSizeList[i];
		}
	}

	// compute minGlobalSize & maxGlobalSize
	OrdinalType numGlobalPoints = numMyPoints;
	OrdinalType minGlobalSize = minMySize;
	OrdinalType maxGlobalSize = maxMySize;
  if(ElementSpace.isGlobal() == true) {
    ElementSpace.comm().sumAll(&numMyPoints, &numGlobalPoints, 1);
    ElementSpace.comm().minAll(&minMySize, &minGlobalSize, 1);
    ElementSpace.comm().maxAll(&maxMySize, &maxGlobalSize, 1);
  }

	// call BESData constructor
	BlockElementSpaceData_.reset(new BlockElementSpaceData<OrdinalType>(ElementSpace, false, 0, numMyPoints, 
																																			numGlobalPoints, minMySize, maxMySize, 
																																			minGlobalSize, maxGlobalSize, myElementSizeList));
}

//=======================================================================
template<typename OrdinalType>
BlockElementSpace<OrdinalType>::BlockElementSpace(BlockElementSpace<OrdinalType> const& BlockElementSpace) 
  : Object(BlockElementSpace.label())
	, BlockElementSpaceData_(BlockElementSpace.BlockElementSpaceData_) {}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getRemoteIDList(OrdinalType numIDs, OrdinalType* GIDList, OrdinalType* imageIDList, 
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
bool BlockElementSpace<OrdinalType>::isSameAs(BlockElementSpace<OrdinalType> const& BlockElementSpace) const {
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
// LID -> size of that element
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getElementSizeList(OrdinalType* elementSizeList) const {
	if(elementSizeList == 0) 
    throw reportError("This pointer does not have a child allocated.", 4);
  OrdinalType nME = elementSpace().getNumMyElements();
  for(OrdinalType i = 0; i < nME; i++)
    elementSizeList[i] = BlockElementSpaceData_->elementSizeList_[i];
}

//=======================================================================
// LID -> lowest PointID contained in that element
template<typename OrdinalType>
OrdinalType const* BlockElementSpace<OrdinalType>::getFirstPointInElementList() const {
	OrdinalType nME = elementSpace().getNumMyElements();
	if((BlockElementSpaceData_->firstPointList_ == 0) && (nME > 0)) {
		OrdinalType* tmpPtr = new OrdinalType[nME];
		getFirstPointInElementList(tmpPtr);
		BlockElementSpaceData_->firstPointList_ = tmpPtr;
	}
	return(BlockElementSpaceData_->firstPointList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getFirstPointInElementList(OrdinalType* firstPointInElementList) const {
	if(firstPointInElementList == 0) 
    throw reportError("This pointer does not have a child allocated.", 4);
  OrdinalType iB = elementSpace().getIndexBase();
	firstPointInElementList[0] = iB;
	OrdinalType nME = elementSpace().getNumMyElements();
	for(OrdinalType i = 1; i < nME; i++)
		firstPointInElementList[i] = firstPointInElementList[i-1] + BlockElementSpaceData_->elementSizeList_[i-1];
}

//=======================================================================
// pointID -> LID containing that point
template<typename OrdinalType>
OrdinalType const* BlockElementSpace<OrdinalType>::getPointToElementList() const {
	OrdinalType numPoints = getNumMyPoints();
	if((BlockElementSpaceData_->pointToElementList_ == 0) && (numPoints > 0)) {
		OrdinalType* tmpPtr = new OrdinalType[numPoints];
		//for(OrdinalType i = 0; i < numPoints; i++)
		//	tmpPtr[i] = i+1;
		getPointToElementList(tmpPtr);
		BlockElementSpaceData_->pointToElementList_ = tmpPtr;
	}
	return(BlockElementSpaceData_->pointToElementList_);
}

//=======================================================================
template<typename OrdinalType>
void BlockElementSpace<OrdinalType>::getPointToElementList(OrdinalType* pointToElementList) const {
	if(pointToElementList == 0) 
    throw reportError("This pointer does not have a child allocated.", 4);
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
	OrdinalType const* elementSizeList1 = getElementSizeList();
  OrdinalType const* firstPointList1 = getFirstPointInElementList(); 
  OrdinalType const* pointToElementList1 = getPointToElementList(); 
 
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

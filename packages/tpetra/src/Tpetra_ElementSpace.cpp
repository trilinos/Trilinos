/*Paul
10-July-2002 SimpleES.cpp - Simple ElementSpace class for testing purposes
17-July-2002 Now the real ElementSpace, being built up in extreme coding style
24-July-2002 all including directory. untemplated. gMGE & print still not const.
27-July-2002 gMGE & print const. Templated for OrdinalType.
05-Aug-2002 switched from PID to imageID
06-Aug-2002 Completed switch to images.
21-Sept-2002 Comm/Platform split
07-Oct-2002 ElementSpaceData move started
13-Oct-2002 Updated constructors to use new ESData constructor.
22-Oct-2002 Modified slightly - ESData constructor now takes Comm* argument
12-Nov-2002 Updated to use createOrdinalComm() instead of createComm()
24-Nov-2002 Updated for imageID methods moved back to Comm
27-Jan-2003 Updated for .hpp and for new const syntax.
*/

#include "Tpetra_Platform.hpp"

namespace Tpetra {

// constructor #1, tpetra contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType indexBase, 
																				Platform<OrdinalType, OrdinalType> const& Platform)
  : Object("Tpetra::ElementSpace")
	, ElementSpaceData_()
{
	// initial throws
	if (numGlobalElements < 0)
    throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= 0.", -1);

	// platform & comm setup
	Comm<OrdinalType, OrdinalType>* comm = Platform.createOrdinalComm();
  OrdinalType numImages = comm->getNumImages();
  OrdinalType myImageID = comm->getMyImageID();
  
	// compute numMyElements
  OrdinalType numMyElements = numGlobalElements / numImages;
  OrdinalType remainder = numGlobalElements % numImages;
  OrdinalType start_index = myImageID * (numMyElements + 1);
  if (myImageID < remainder)
    numMyElements++;
  else
    start_index -= (myImageID - remainder);

	// setup lgmap & glmap
	map<OrdinalType, OrdinalType> lgMap;
	map<OrdinalType, OrdinalType> glMap;
  
	// setup min/maxs
  OrdinalType minAllGID = indexBase;
  OrdinalType maxAllGID = minAllGID + numGlobalElements - 1;
  OrdinalType minMyGID = start_index + indexBase;
  OrdinalType maxMyGID = minMyGID + numMyElements - 1;
	
	// call ESData constructor
	ElementSpaceData_.reset(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
																														minMyGID, maxMyGID, lgMap, glMap, true, Platform, comm));
  
	// initialize directory
  directorySetup();
}

// constructor #2, user contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType indexBase, 
																				Platform<OrdinalType, OrdinalType> const& Platform)
	: Object("Tpetra::ElementSpace")
	, ElementSpaceData_()
{
	// initial throws
  if(numGlobalElements < -1) 
    throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= -1.", -1);
  if(numMyElements < 0) 
    throw reportError("numMyElements = " + toString(numMyElements) + ".  Should be >= 0.", -2);

	// platform & comm setup
	Comm<OrdinalType, OrdinalType>* comm = Platform.createOrdinalComm();
  OrdinalType numImages = comm->getNumImages();
  OrdinalType myImageID = comm->getMyImageID();

	// check for invalid numGlobalElements
	//   Sum up all local element counts to get global count, and then
	//   check to see if user's value for numGlobalElements is either -1 
	//   (in which case we use our computed value) or matches ours.
  OrdinalType global_sum;
  comm->sumAll(&numMyElements, &global_sum, 1);
	if(numGlobalElements == -1)
		numGlobalElements = global_sum;
	else if(numGlobalElements != global_sum) 
		throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
											".  Should equal " + toString(global_sum) + ", or be set to -1 to compute automatically", -3);

	// setup lgmap & glmap
	map<OrdinalType, OrdinalType> lgMap;
	map<OrdinalType, OrdinalType> glMap;
	
	// setup min/maxs
  OrdinalType minAllGID = indexBase;
  OrdinalType maxAllGID = minAllGID + numGlobalElements - 1;
	OrdinalType start_index;
	comm->scanSum(&numMyElements, &start_index, 1);
	start_index -= numMyElements;
	OrdinalType minMyGID = start_index + indexBase;
	OrdinalType maxMyGID = minMyGID + numMyElements - 1;

	// call ESData constructor
	ElementSpaceData_.reset(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
																														minMyGID, maxMyGID, lgMap, glMap, true, Platform, comm));
  
	// initialize directory
  directorySetup();
}

// constructor #3, user non-contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType* elementList, 
																				OrdinalType indexBase, Platform<OrdinalType, OrdinalType> const& Platform)
  : Object("Tpetra::ElementSpace")
	, ElementSpaceData_()
{
	// initial throws
  if(numGlobalElements < -1) 
    throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= -1.", -1);
  if(numMyElements < 0) 
    throw reportError("numMyElements = " + toString(numMyElements) + ".  Should be >= 0.", -2);

	// platform & comm setup
	Comm<OrdinalType, OrdinalType>* comm = Platform.createOrdinalComm();
  OrdinalType numImages = comm->getNumImages();
  OrdinalType myImageID = comm->getMyImageID();

	// check for invalid numGlobalElements
  //   Sum up all local element counts to get global count, and then
	//   check to see if user's value for numGlobalElements is either -1 
	//   (in which case we use our computed value) or matches ours.
  OrdinalType global_sum;
  comm->sumAll(&numMyElements, &global_sum, 1);
	if(numGlobalElements == -1)
		numGlobalElements = global_sum;
	else if(numGlobalElements != global_sum)
		throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
											".  Should equal " + toString(global_sum) + ", or be set to -1 to compute automatically", -3);
  
	// setup lgmap and glmap, and min/maxMyGIDs
	map<OrdinalType, OrdinalType> lgMap;
	map<OrdinalType, OrdinalType> glMap;
	OrdinalType minMyGID = indexBase;
	OrdinalType maxMyGID = indexBase;
	if(numMyElements > 0) {
		for(OrdinalType i = 0; i < numMyElements; i++) {
			lgMap[i + indexBase] = elementList[i]; // lgmap: LID=key, GID=mapped
			glMap[elementList[i]] = (i + indexBase); // glmap: GID=key, LID=mapped
		}
    minMyGID = elementList[0];
    maxMyGID = elementList[numMyElements - 1];
	}

	// set min/maxAllGIDs
	OrdinalType minAllGID;
	OrdinalType maxAllGID;
	comm->minAll(&minMyGID, &minAllGID, 1);
	comm->maxAll(&maxMyGID, &maxAllGID, 1);
  if (minAllGID < indexBase)
    throw reportError("Minimum global element index = " + toString(minAllGID) + 
											" is less than index base = " + toString(indexBase) +".", -4);
	
	// call ESData constructor
	ElementSpaceData_.reset(new ElementSpaceData<OrdinalType>(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID, 
																														minMyGID, maxMyGID, lgMap, glMap, false, Platform, comm));

	// initialize directory
  directorySetup();
}

// copy constructor
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace (ElementSpace<OrdinalType> const& ElementSpace) 
  : Object(ElementSpace.label())
	, ElementSpaceData_(ElementSpace.ElementSpaceData_)
{}

//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::~ElementSpace() {}

//=======================================================================
template<typename OrdinalType>
OrdinalType ElementSpace<OrdinalType>::getLID (OrdinalType GID) const {
  if(!isMyGID(GID)) 
    throw reportError("Global ID " + toString(GID) + " was not found on this processor.", 1);
  else if(isContiguous()) 
    return(GID - getMinMyGID() + getIndexBase()); //compute with offset
  else {
    return((ElementSpaceData_->glMap_.find(GID))->second);
	}
}

//=======================================================================
template<typename OrdinalType>
OrdinalType ElementSpace<OrdinalType>::getGID (OrdinalType LID) const {
  if(!isMyLID(LID))
    throw reportError("Local ID " + toString(LID) + " was not found on this processor.", 2);
  else if(isContiguous()) 
    return(LID + getMinMyGID() + getIndexBase()); //compute with offset
  else {
    return((ElementSpaceData_->lgMap_.find(LID))->second);
	}
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isMyGID (OrdinalType GID) const {
  if(GID < getMinMyGID() || GID > getMaxMyGID())
    return(false);
  else if(isContiguous())
    return(true);
  else {
    return (ElementSpaceData_->glMap_.find(GID) != ElementSpaceData_->glMap_.end());
	}
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isMyLID (OrdinalType LID) const {
  if(LID < getMinLID() || LID > getMaxLID())
    return(false);
  else if(isContiguous())
    return(true);
  else {
    return (ElementSpaceData_->lgMap_.find(LID) != ElementSpaceData_->lgMap_.end());
	}
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::getMyGlobalElements(OrdinalType* elementList) const {
  if(elementList == 0)
    throw reportError("Pointer does not have child allocated.", 3);
  else if(isContiguous()) {
		OrdinalType nME = getNumMyElements();
		OrdinalType minMyGID = getMinMyGID();
    for(OrdinalType i = 0; i < nME; i++)
      elementList[i] = minMyGID + i;
	}
  else { // not contiguous
    map<OrdinalType, OrdinalType>::iterator lgi = ElementSpaceData_->lgMap_.begin();
    map<OrdinalType, OrdinalType>::iterator lgmax = ElementSpaceData_->lgMap_.end();
    for(OrdinalType i = 0; lgi != lgmax; i++) {
      elementList[i] = lgi->second;
      lgi++;
    }
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType* ElementSpace<OrdinalType>::getMyGlobalElements() const {
	OrdinalType nME = getNumMyElements();
  if((ElementSpaceData_->myGlobalElements_ == 0) && (nME > 0)) {
		ElementSpaceData_->myGlobalElements_ = new OrdinalType[nME];
    getMyGlobalElements(ElementSpaceData_->myGlobalElements_);
  }
  return(ElementSpaceData_->myGlobalElements_);
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isSameAs (ElementSpace<OrdinalType> const& ElementSpace) const {
  if (this == &ElementSpace) 
    return(true);
  if (getMinAllGID() != ElementSpace.getMinAllGID() || 
      getMaxAllGID() != ElementSpace.getMaxAllGID() ||
      getNumGlobalElements() != ElementSpace.getNumGlobalElements() || 
      getIndexBase() != ElementSpace.getIndexBase() ||
      isGlobal() != ElementSpace.isGlobal() || 
      isContiguous() != ElementSpace.isContiguous()) 
    return(false);

  // If we get this far, we need to check local properties and then check across
  // all processors to see if local properties are all true
	
  int mySameSpace = 1;
  if(getNumMyElements() != ElementSpace.getNumMyElements()) 
    mySameSpace=0;
	
  if(!isContiguous() && mySameSpace == 1)
    if(ElementSpaceData_->lgMap_ != ElementSpace.ElementSpaceData_->lgMap_)
      mySameSpace=0;

  // Now get min of mySameSpace across all processors
  int globalSameSpace = 0;
  comm().minAll(&mySameSpace, &globalSameSpace, 1);
  return(globalSameSpace==1);
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::print(ostream& os) const {
  OrdinalType* myGlobalElements1 = getMyGlobalElements();
  OrdinalType myImageID = comm().getMyImageID();
  OrdinalType numImages = comm().getNumImages();
	OrdinalType minLID = getMinLID();
	OrdinalType nME = getNumMyElements();
  
  for (int imageCtr = 0; imageCtr < numImages; imageCtr++) {
    if (myImageID == imageCtr) {
      if (myImageID == 0) {
				os <<  "\nNumber of Global Elements  = "; os << getNumGlobalElements(); os << endl;
				os <<    "Maximum of all GIDs        = "; os << getMaxAllGID(); os << endl;
				os <<    "Minimum of all GIDs        = "; os << getMinAllGID(); os << endl;
				os <<    "Index Base                 = "; os << getIndexBase(); os << endl;
      }
      os << endl;

      os <<    "Number of Local Elements   = "; os << getNumMyElements(); os << endl;
      os <<    "Maximum of my GIDs         = "; os << getMaxMyGID(); os << endl;
      os <<    "Minimum of my GIDs         = "; os << getMinMyGID(); os << endl;
      os << endl;

      os.width(14);
      os <<  "   ImageID"; os << "    ";
      os.width(14);
      os <<  "       Local Index "; os << " ";
      os.width(14);
      os <<  "      Global Index "; os << " ";
      os << endl;
    
      for (OrdinalType i = 0, lid = minLID; i < nME; i++, lid++) {
				os.width(14);
				os <<  myImageID; os << "    ";
				os.width(14);
				os << lid; os << "    ";
				os.width(14);
				os <<  myGlobalElements1[i]; os << "    ";
				os << endl;
      }
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    comm().barrier();
    comm().barrier();
    comm().barrier();
  }
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::directorySetup() {
  if(getNumGlobalElements() != 0)
    if(ElementSpaceData_->Directory_ == 0)
      ElementSpaceData_->Directory_ = platform().createDirectory(*this); // Make directory
}

} // namespace Tpetra
//=======================================================================

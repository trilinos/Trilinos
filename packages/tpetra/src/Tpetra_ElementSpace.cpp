/*Paul
10-July-2002 SimpleES.cpp - Simple ElementSpace class for testing purposes
17-July-2002 Now the real ElementSpace, being built up in extreme coding style
24-July-2002 all including directory. untemplated. gMGE & print still not const.
27-July-2002 gMGE & print const. Templated for OrdinalType.
05-Aug-2002 switched from PID to imageID
06-Aug-2002 Completed switch to images.
*/

#include "Tpetra_Comm.h"

namespace Tpetra {

// constructor #1, tpetra contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType indexBase, const Comm<OrdinalType, OrdinalType>& Comm)
  : Object("Tpetra::ElementSpace")
   ,numGlobalElements_(numGlobalElements)
   ,indexBase_(indexBase)
   ,contiguous_(true)
   ,lgMap_()
   ,glMap_()
   ,Comm_(&Comm)
   ,myGlobalElements_(0)
   ,Directory_(0)
{
  if (numGlobalElements_ < 0)
    throw reportError("numGlobalElements = " + toString(numGlobalElements_) + ".  Should be >= 0.", -1);
  
  global_ = checkGlobalness(numGlobalElements_, numMyElements_);

  OrdinalType numImages = comm().getNumImages();
  OrdinalType myImageID = comm().getMyImageID();
  
  numMyElements_ = numGlobalElements_ / numImages;
  OrdinalType remainder = numGlobalElements_ % numImages;
  OrdinalType start_index = myImageID * (numMyElements_ + 1);
  
  if (myImageID < remainder)
    numMyElements_++;
  else
    start_index -= (myImageID - remainder);
  
  minAllGID_ = indexBase_;
  maxAllGID_ = minAllGID_ + numGlobalElements_ - 1;
  minMyGID_ = start_index + indexBase_;
  maxMyGID_ = minMyGID_ + numMyElements_ - 1;
  minLID_ = indexBase_;
  maxLID_ = minLID_ + numMyElements_ - 1;
	
  directorySetup();
}

// constructor #2, user contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType indexBase, 
					const Comm<OrdinalType, OrdinalType>& Comm)
  : Object("Tpetra::ElementSpace")
   ,numGlobalElements_(numGlobalElements)
   ,numMyElements_(numMyElements)
   ,indexBase_(indexBase)
   ,contiguous_(true)
   ,lgMap_()
   ,glMap_()
   ,Comm_(&Comm)
   ,myGlobalElements_(0)
   ,Directory_(0)
{
  if(numGlobalElements_ < -1) 
    throw reportError("numGlobalElements = " + toString(numGlobalElements_) + ".  Should be >= -1.", -1);
  if(numMyElements_ < 0) 
    throw reportError("numMyElements = " + toString(numGlobalElements_) + ".  Should be >= 0.", -2);

  global_ = checkGlobalness(numGlobalElements, numMyElements);

  OrdinalType numImages = comm().getNumImages();
  OrdinalType myImageID = comm().getMyImageID();

  // Locally replicated and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!global_ || numImages==1) {
    numGlobalElements_ = numMyElements_;
    // Check to see if user's value for numGlobalElements is either indexBase-1 
    // (in which case we use our computed value) or matches ours.
    if ((numGlobalElements != -1) && (numGlobalElements != numGlobalElements_)) 
      throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
			".  Should equal " + toString(numGlobalElements_) + 
			", or be set to -1 to compute automatically", -3);
    
    minAllGID_ = indexBase_;
    maxAllGID_ = minAllGID_ + numGlobalElements_ - 1;
    minMyGID_ = indexBase_;
    maxMyGID_ = minMyGID_ + numMyElements_ - 1;
    minLID_ = indexBase_;
    maxLID_ = minLID_ + numMyElements_ - 1;
  }
  else if (numImages > 1) {
    // Sum up all local element counts to get global count
    comm().sumAll(&numMyElements_, &numGlobalElements_, 1);

    // Check to see if user's value for numGlobalElements is either indexBase-1 
    // (in which case we use our computed value) or matches ours.
    if ((numGlobalElements != -1) && (numGlobalElements != numGlobalElements_))
      throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
			".  Should equal " + toString(numGlobalElements_) + 
			", or be set to -1 to compute automatically", -3);
    
    minAllGID_ = indexBase_;
    maxAllGID_ = minAllGID_ + numGlobalElements_ - 1;
    minLID_ = indexBase_;
    maxLID_ = minLID_ + numMyElements_ - 1;

    comm().scanSum(&numMyElements_, &maxMyGID_, 1);

    OrdinalType start_index = maxMyGID_ - numMyElements_;
    minMyGID_ = start_index + indexBase_;
    maxMyGID_ = minMyGID_ + numMyElements_ - 1;
  }
  else 
    throw reportError("Internal Error.  Report to Tpetra developer", -99);
  
  directorySetup();
}

// constructor #3, user non-contig
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace(OrdinalType numGlobalElements, OrdinalType numMyElements, OrdinalType* elementList, 
					OrdinalType indexBase, const Comm<OrdinalType, OrdinalType>& Comm)
  : Object("Tpetra::ElementSpace")
    ,numGlobalElements_(numGlobalElements)
    ,myGlobalElements_(0)
    ,numMyElements_(numMyElements)
    ,indexBase_(indexBase)
    ,contiguous_(false)
    ,lgMap_()
    ,glMap_()
    ,Comm_(&Comm)
    ,Directory_(0)
{
  if (numGlobalElements_ < - 1) 
    throw reportError("numGlobalElements = " + toString(numGlobalElements_) + ".  Should be >= -1.", -1);
  if (numMyElements_ < 0) 
    throw reportError("numMyElements = " + toString(numGlobalElements_) + ".  Should be >= 0.", -2);

  OrdinalType numImages = comm().getNumImages();
  
  if (numMyElements > 0) {
    for(OrdinalType i = 0; i < numMyElements_; i++) {
      lgMap_[i + indexBase_] = elementList[i]; // lgmap: LID=key, GID=mapped
      glMap_[elementList[i]] = (i + indexBase_); // glmap: GID=key, LID=mapped
    }
    minMyGID_ = elementList[0];
    maxMyGID_ = elementList[numMyElements_ - 1];
    minLID_ = indexBase_;
    maxLID_ = minLID_ + numMyElements_ - 1;
  }
  else {
    minMyGID_ = indexBase_;
    maxMyGID_ = indexBase_;
    minLID_ = indexBase_;
    maxLID_ = indexBase_;
  }
  
  global_ = checkGlobalness(numGlobalElements, numMyElements);

  // Local Map and uniprocessor case:  Each processor gets a complete copy of all elements
  if (!global_ || numImages == 1) {
    numGlobalElements_ = numMyElements_;
    // Check to see if user's value for numGlobalElements is either indexBase-1 
    // (in which case we use our computed value) or matches ours.
    if ((numGlobalElements != -1) && (numGlobalElements != numGlobalElements_)) 
      throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
			".  Should equal " + toString(numGlobalElements_) + 
			", or be set to -1 to compute automatically", -3);
      
    minAllGID_ = minMyGID_;
    maxAllGID_ = maxMyGID_;
  }
  else if (numImages > 1) {
    // Sum up all local element counts to get global count
    comm().sumAll(&numMyElements_, &numGlobalElements_, 1);
    // Check to see if user's value for NumGlobalElements is either -1
    // (in which case we use our computed value) or matches ours.
    if ((numGlobalElements != -1) && (numGlobalElements != numGlobalElements_)) 
      throw reportError("Invalid numGlobalElements.  numGlobalElements = " + toString(numGlobalElements) + 
			".  Should equal " + toString(numGlobalElements_) + 
			", or be set to -1 to compute automatically", -3);
      
    // Use the Allreduce function to find min/max GID 
    comm().minAll(&minMyGID_, &minAllGID_, 1);
    comm().maxAll(&maxMyGID_, &maxAllGID_, 1);
  }
  else
    throw reportError("Internal Error.  Report to Tpetra developer", -99);

  if (minAllGID_ < indexBase_)
    throw reportError("Minimum global element index = " + toString(minAllGID_) + " is less than index base = " + toString(indexBase_) +".", -4);

  directorySetup();
}

// copy constructor
//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::ElementSpace (const ElementSpace<OrdinalType>& ElementSpace) 
  : Object(ElementSpace.label())
    ,numGlobalElements_(ElementSpace.numGlobalElements_)
    ,numMyElements_(ElementSpace.numMyElements_)
    ,indexBase_(ElementSpace.indexBase_)
    ,minLID_(ElementSpace.minLID_)
    ,maxLID_(ElementSpace.maxLID_)
    ,minMyGID_(ElementSpace.minMyGID_)
    ,maxMyGID_(ElementSpace.maxMyGID_)
    ,minAllGID_(ElementSpace.minAllGID_)
    ,maxAllGID_(ElementSpace.maxMyGID_)
    ,contiguous_(ElementSpace.contiguous_)
    ,global_(ElementSpace.global_)
    ,lgMap_(ElementSpace.lgMap_)
    ,glMap_(ElementSpace.glMap_)
    ,Comm_(ElementSpace.Comm_)
    ,myGlobalElements_(0)
    ,Directory_(0)
{
  // Create mGE array if ElementSpace had one
  if(ElementSpace.myGlobalElements_ != 0)
    OrdinalType* tmp = getMyGlobalElements();
  // Create directory if ElementSpace had one
  if(ElementSpace.Directory_ != 0)
    directorySetup();
}

//=======================================================================
template<typename OrdinalType>
ElementSpace<OrdinalType>::~ElementSpace() {
  if(Directory_ != 0) {
    delete Directory_;
    Directory_ = 0;
  }
  if(myGlobalElements_ != 0) {
    delete [] myGlobalElements_;
    myGlobalElements_ = 0;
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType ElementSpace<OrdinalType>::getLID (OrdinalType GID) const {
  if(!isMyGID(GID)) 
    throw reportError("Global ID " + toString(GID) + " was not found on this processor.", 1);
  else if(isContiguous()) 
    return(GID - getMinMyGID() + getIndexBase()); //compute with offset
  else {
    return((glMap_.find(GID))->second);
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
    return((lgMap_.find(LID))->second);
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
    return (glMap_.find(GID) != glMap_.end());
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
    return (lgMap_.find(LID) != lgMap_.end());
	}
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::getMyGlobalElements(OrdinalType* elementList) const {
  if(elementList == 0)
    throw reportError("Pointer does not have child allocated.", 3);
  else if(isContiguous())
    for(OrdinalType i = 0; i < numMyElements_; i++)
      elementList[i] = minMyGID_ + i;
  else { // not contiguous
    map<OrdinalType, OrdinalType>::iterator lgi = lgMap_.begin();
    for(OrdinalType i = 0; lgi != lgMap_.end(); i++) {
      elementList[i] = lgi->second;
      lgi++;
    }
  }
}

//=======================================================================
template<typename OrdinalType>
OrdinalType* ElementSpace<OrdinalType>::getMyGlobalElements() const {
	OrdinalType numMyElements = getNumMyElements();
  if((myGlobalElements_ == 0) && (numMyElements > 0)) {
		myGlobalElements_ = new OrdinalType[numMyElements];
    getMyGlobalElements(myGlobalElements_);
  }
  return(myGlobalElements_);
}

//=======================================================================
template<typename OrdinalType>
bool ElementSpace<OrdinalType>::isSameAs (const ElementSpace<OrdinalType>& ElementSpace) const {
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
    if(lgMap_ != ElementSpace.lgMap_)
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
    
      for (OrdinalType i = 0, lid = minLID_; i < numMyElements_; i++, lid++) {
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
bool ElementSpace<OrdinalType>::checkGlobalness(OrdinalType numGlobalElements, OrdinalType numMyElements) {
  bool global = false;
  if(comm().getNumImages() > 1) {
    int localRep = 0;
    int allLocalRep;
    if(numGlobalElements == numMyElements)
      localRep=1;
    comm().minAll(&localRep, &allLocalRep, 1);
    if(allLocalRep != 1)
      global = true;
  }
  return(global);
}

//=======================================================================
template<typename OrdinalType>
void ElementSpace<OrdinalType>::directorySetup() {
  if(getNumGlobalElements() != 0)
    if(Directory_ == 0)
      Directory_ = comm().createDirectory(*this); // Make directory
}


} // namespace Tpetra
//=======================================================================

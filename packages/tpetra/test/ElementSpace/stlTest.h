#include <map>
#include <iostream>
#include <iomanip>

class stlTest {
public:
  stlTest(ordinalType nGE, ordinalType nME, ordinalType* mGE, ordinalType iB);
	~stlTest();
  bool isMyLID(ordinalType LID);
  bool isMyGID(ordinalType GID);
  ordinalType getLID(ordinalType GID);
  ordinalType getGID(ordinalType LID);
	void getMyGlobalElements(ordinalType* list) const;
	ordinalType* getMyGlobalElements() const;
	void print(ostream& os) const;

	ordinalType getNumGlobalElements() const {return(numGlobalElements_);};
	ordinalType getNumMyElements() const {return(numMyElements_);};
	ordinalType getIndexBase() const {return(indexBase_);};
	ordinalType getMinLID() const {return(minLID_);};
	ordinalType getMaxLID() const {return(maxLID_);};
	ordinalType getMinMyGID() const {return(minMyGID_);};
	ordinalType getMaxMyGID() const {return(maxMyGID_);};
	ordinalType getMinAllGID() const {return(minAllGID_);};
	ordinalType getMaxAllGID() const {return(maxAllGID_);};

private:
	ordinalType numGlobalElements_;
	ordinalType numMyElements_;
	ordinalType indexBase_;
	ordinalType minLID_;
	ordinalType maxLID_;
 	ordinalType minMyGID_;
	ordinalType maxMyGID_;
	ordinalType minAllGID_;
	ordinalType maxAllGID_;
	mutable map<ordinalType, ordinalType> lgMap_;
	mutable map<ordinalType, ordinalType> glMap_;
  mutable ordinalType* myGlobalElements_;
};

stlTest::stlTest(ordinalType nGE, ordinalType nME, ordinalType* elementList, ordinalType iB)
	: numGlobalElements_(nGE)
		 ,numMyElements_(nME)
		 ,indexBase_(iB)
		 ,lgMap_()
		 ,glMap_()
		 ,myGlobalElements_(0)
{
	if (numMyElements_ > 0) {
    for(ordinalType i = 0; i < numMyElements_; i++) {
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

	maxAllGID_ = maxMyGID_;
	minAllGID_ = minMyGID_;
}

stlTest::~stlTest() {
	if(myGlobalElements_ != 0) {
    delete myGlobalElements_;
		myGlobalElements_ = 0;
	}
}

void stlTest::print(ostream& os) const {
	cout << "stlTest::print" << endl;
  ordinalType* myGlobalElements1 = getMyGlobalElements();
  
	os <<  "\nNumber of Global Elements  = "; os << getNumGlobalElements(); os << endl;
	os <<    "Maximum of all GIDs        = "; os << getMaxAllGID(); os << endl;
	os <<    "Minimum of all GIDs        = "; os << getMinAllGID(); os << endl;
	os <<    "Index Base                 = "; os << getIndexBase(); os << endl;
  os << endl;

  os <<    "Number of Local Elements   = "; os << getNumMyElements(); os << endl;
  os <<    "Maximum of my GIDs         = "; os << getMaxMyGID(); os << endl;
  os <<    "Minimum of my GIDs         = "; os << getMinMyGID(); os << endl;
  os << endl;

  os.width(14);
  os <<  "       Local Index "; os << " ";
  os.width(14);
  os <<  "      Global Index "; os << " ";
  os << endl;
    
  for (ordinalType i = 0, lid = minLID_; i < numMyElements_; i++, lid++) {
		os.width(14);
		os << lid; os << "    ";
		os.width(14);
		os <<  myGlobalElements1[i]; os << "    ";
		os << endl;
  }
  os << flush;
}

void stlTest::getMyGlobalElements(ordinalType* list) const {
	map<ordinalType, ordinalType>::iterator gi = lgMap_.begin();
	for(ordinalType i = 0; gi != lgMap_.end(); i++) {
	  list[i] = gi->second;
		gi++;
	}
}

ordinalType* stlTest::getMyGlobalElements() const {
	if(myGlobalElements_ == 0 && getNumMyElements() > 0) {
	  myGlobalElements_ = new ordinalType[numMyElements_ + 1];
		getMyGlobalElements(myGlobalElements_);
	}
	return(myGlobalElements_);
}

bool stlTest::isMyLID(ordinalType LID) {
	return(lgMap_.find(LID) != lgMap_.end());
}

bool stlTest::isMyGID(ordinalType GID) {
	return(glMap_.find(GID) != glMap_.end());
}

ordinalType stlTest::getLID(ordinalType GID) {
	if(!isMyGID(GID)) {
		cout << "GID " << GID << " not found." << endl;
		return(0);
	}
	else
		return((glMap_.find(GID))->second);
}

ordinalType stlTest::getGID(ordinalType LID) {
	if(!isMyLID(LID)) {
		cout << "LID " << LID << " not found." << endl;
		return(0);
	}
	else
		return((lgMap_.find(LID))->second);
}

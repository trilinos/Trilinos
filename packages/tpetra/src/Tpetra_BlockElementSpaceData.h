/*Paul
09-Oct-2002 BlockElementSpaceData move started.
14-Oct-2002 Continued.
12-Nov-2002 Updated to use createOrdinalComm() instead of createComm() (nothing changed)
24-Nov-2002 Rewritten with massive constructor call.
*/

#ifndef _TPETRA_BLOCKELEMENTSPACEDATA_H_
#define _TPETRA_BLOCKELEMENTSPACEDATA_H_

namespace Tpetra {

template<typename OrdinalType>
class BlockElementSpaceData : public Object {
	friend class BlockElementSpace<OrdinalType>;
 public:
	BlockElementSpaceData(const ElementSpace<OrdinalType>& ElementSpace, 
												const bool constantSize, 
												const OrdinalType elementSize,
												const OrdinalType numMyPoints,
												const OrdinalType numGlobalPoints,
												const OrdinalType minMySize,
												const OrdinalType maxMySize,
												const OrdinalType minGlobalSize,
												const OrdinalType maxGlobalSize,
												const OrdinalType* elementSizeList) 
		: Object("Tpetra::BlockElementSpaceData")
		, ElementSpace_(&ElementSpace) 
		, constantSize_(constantSize)
		, elementSize_(elementSize)
		, numMyPoints_(numMyPoints)
		, numGlobalPoints_(numGlobalPoints)
		, minMySize_(minMySize)
		, maxMySize_(maxMySize)
		, minGlobalSize_(minGlobalSize)
		, maxGlobalSize_(maxGlobalSize)
		, elementSizeList_(elementSizeList)
		, pointToElementList_(0)
		, firstPointList_(0) 
		{};

	~BlockElementSpaceData() {
		//cout << "BESData destructor called." << endl;
		//cout << "BESData destructor starting firstPointList." << endl;
		if(firstPointList_ != 0) {
			delete[] firstPointList_;
			firstPointList_ = 0;
		}
		//cout << "BESData destructor starting pointToElementList." << endl;
		if(pointToElementList_ != 0) {
			delete[] pointToElementList_;
			pointToElementList_ = 0;
		}
		//cout << "BESData destructor starting elementSizeList." << endl;
		if(elementSizeList_ != 0) {
			//cout << "BESData destructor about to delete eSL." << endl;
			delete[] elementSizeList_;
			elementSizeList_ = 0;
		}
		//cout << "BESData destructor finished." << endl;
	};

 protected:
	const ElementSpace<OrdinalType>* const ElementSpace_;
	const bool constantSize_;
	const OrdinalType elementSize_;
	const OrdinalType numMyPoints_;
	const OrdinalType numGlobalPoints_;
	const OrdinalType minMySize_;
	const OrdinalType maxMySize_;
	const OrdinalType minGlobalSize_;
	const OrdinalType maxGlobalSize_;
	const OrdinalType* elementSizeList_;
	const OrdinalType* pointToElementList_;
	const OrdinalType* firstPointList_;

}; // class BlockElementSpaceData

} // namespace Tpetra

#endif // _TPETRA_BLOCKELEMENTSPACEDATA_H_

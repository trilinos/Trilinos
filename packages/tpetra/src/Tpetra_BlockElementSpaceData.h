/*Paul
09-Oct-2002 BlockElementSpaceData move started.
14-Oct-2002 Continued.
*/

#ifndef _TPETRA_BLOCKELEMENTSPACEDATA_H_
#define _TPETRA_BLOCKELEMENTSPACEDATA_H_

namespace Tpetra {

template<typename OrdinalType>
class BlockElementSpaceData : public Object {
	friend class BlockElementSpace<OrdinalType>;
 public:
	BlockElementSpaceData(ElementSpace<OrdinalType>& ElementSpace, bool constantSize) 
		: Object("Tpetra::BlockElementSpaceData")
		, constantSize_(constantSize)
		, elementSizeList_(0)
		, pointToElementList_(0)
		, firstPointList_(0) 
		, ElementSpace_(&ElementSpace) {};

	~BlockElementSpaceData() {
		if(elementSizeList_ != 0) {
			delete[] elementSizeList_;
			elementSizeList_ = 0;
		}
		if(pointToElementList_ != 0) {
			delete[] pointToElementList_;
			pointToElementList_ = 0;
		}
		if(pointToElementList_ != 0) {
			delete[] pointToElementList_;
			pointToElementList_ = 0;
		}
	};

	protected:

	bool constantSize_;
	OrdinalType elementSize_;
	OrdinalType numMyPoints_;
	OrdinalType numGlobalPoints_;
	OrdinalType minMySize_;
	OrdinalType maxMySize_;
	OrdinalType minGlobalSize_;
	OrdinalType maxGlobalSize_;
	OrdinalType* elementSizeList_;
	OrdinalType* pointToElementList_;
	OrdinalType* firstPointList_;
	ElementSpace<OrdinalType>* ElementSpace_;

}; // class BlockElementSpaceData

} // namespace Tpetra

#endif // _TPETRA_BLOCKELEMENTSPACEDATA_H_

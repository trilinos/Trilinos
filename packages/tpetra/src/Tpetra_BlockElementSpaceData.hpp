/*Paul
09-Oct-2002 BlockElementSpaceData move started.
14-Oct-2002 Continued.
12-Nov-2002 Updated to use createOrdinalComm() instead of createComm() (nothing changed)
24-Nov-2002 Rewritten with massive constructor call.
27-Jan-2003 Updated for .hpp and for new const syntax.
*/

#ifndef _TPETRA_BLOCKELEMENTSPACEDATA_HPP_
#define _TPETRA_BLOCKELEMENTSPACEDATA_HPP_

namespace Tpetra {

template<typename OrdinalType>
class BlockElementSpaceData : public Object {
	friend class BlockElementSpace<OrdinalType>;
 public:
	BlockElementSpaceData(ElementSpace<OrdinalType> const& ElementSpace, 
												bool const constantSize, 
												OrdinalType const elementSize,
												OrdinalType const numMyPoints,
												OrdinalType const numGlobalPoints,
												OrdinalType const minMySize,
												OrdinalType const maxMySize,
												OrdinalType const minGlobalSize,
												OrdinalType const maxGlobalSize,
												OrdinalType const* elementSizeList) 
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
		if(firstPointList_ != 0) {
			delete[] firstPointList_;
			firstPointList_ = 0;
		}
		if(pointToElementList_ != 0) {
			delete[] pointToElementList_;
			pointToElementList_ = 0;
		}
		if(elementSizeList_ != 0) {
			delete[] elementSizeList_;
			elementSizeList_ = 0;
		}
	};

 protected:
	ElementSpace<OrdinalType> const* ElementSpace_;
	bool const constantSize_;
	OrdinalType const elementSize_;
	OrdinalType const numMyPoints_;
	OrdinalType const numGlobalPoints_;
	OrdinalType const minMySize_;
	OrdinalType const maxMySize_;
	OrdinalType const minGlobalSize_;
	OrdinalType const maxGlobalSize_;
	OrdinalType const* elementSizeList_;
	OrdinalType const* pointToElementList_;
	OrdinalType const* firstPointList_;

}; // class BlockElementSpaceData

} // namespace Tpetra

#endif // _TPETRA_BLOCKELEMENTSPACEDATA_HPP_

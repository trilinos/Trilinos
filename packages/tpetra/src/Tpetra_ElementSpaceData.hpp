/*Paul
26-Jan-2003 Updated for .hpp and for new const syntax.
*/

#ifndef _TPETRA_ELEMENTSPACEDATA_HPP_
#define _TPETRA_ELEMENTSPACEDATA_HPP_

namespace Tpetra {

template<typename OrdinalType>
class ElementSpaceData : public Object {
	friend class ElementSpace<OrdinalType>;
 public:
	ElementSpaceData(OrdinalType const indexBase, 
									 OrdinalType const numGlobalElements,
									 OrdinalType const numMyElements,
									 OrdinalType const minAllGID,
									 OrdinalType const maxAllGID,
									 OrdinalType const minMyGID,
									 OrdinalType const maxMyGID,
									 map<OrdinalType, OrdinalType> const lgMap,
									 map<OrdinalType, OrdinalType> const glMap,
									 bool const contiguous,
									 Platform<OrdinalType, OrdinalType> const& platform,
									 Comm<OrdinalType, OrdinalType> const* comm)
		: Object("Tpetra::ElementSpaceData")
		, Platform_(&platform) 
		, Comm_(comm) 
		, numGlobalElements_(numGlobalElements)
		, numMyElements_(numMyElements)
		, indexBase_(indexBase)
		, minLID_(indexBase)
		, maxLID_(indexBase + numMyElements)
		, minMyGID_(minMyGID)
		, maxMyGID_(maxMyGID)
		, minAllGID_(minAllGID)
		, maxAllGID_(maxAllGID)
		, contiguous_(contiguous)
		, global_(checkGlobalness())
		, lgMap_(lgMap)
		, glMap_(glMap)
		, myGlobalElements_(0)
		, Directory_(0) 
		{};

	~ElementSpaceData() {
		if(Directory_ != 0) {
			delete Directory_;
			Directory_ = 0;
		}
		if(myGlobalElements_ != 0) {
			delete [] myGlobalElements_;
			myGlobalElements_ = 0;
		}
		if(Comm_ != 0) {
			delete Comm_;
			Comm_ = 0;
		}
	};

 protected:
	Platform<OrdinalType, OrdinalType> const* Platform_;
	Comm<OrdinalType, OrdinalType> const* Comm_;
	OrdinalType const numGlobalElements_;
	OrdinalType const numMyElements_;
	OrdinalType const indexBase_;
	OrdinalType const minLID_;
	OrdinalType const maxLID_;
	OrdinalType const minMyGID_;
	OrdinalType const maxMyGID_;
	OrdinalType const minAllGID_;
	OrdinalType const maxAllGID_;
	bool const contiguous_;
	bool const global_;
  map<OrdinalType, OrdinalType> lgMap_;
  map<OrdinalType, OrdinalType> const glMap_;
	OrdinalType* myGlobalElements_;
	Directory<OrdinalType>* Directory_;

 private:
	bool checkGlobalness() {
		bool global = false;
		if(Comm_->getNumImages() > 1) {
			int localRep = 0;
			int allLocalRep;
			if(numGlobalElements_ == numMyElements_)
				localRep = 1;
			Comm_->minAll(&localRep, &allLocalRep, 1);
			if(allLocalRep != 1)
				global = true;
		}
		return(global);
	}

}; // class ElementSpaceData

} // namespace Tpetra

#endif // _TPETRA_ELEMENTSPACEDATA_HPP_

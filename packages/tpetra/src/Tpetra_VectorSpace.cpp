/*Paul
05-Feb-2003
*/

#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"
#include "Tpetra_Platform.hpp"
#include "Tpetra_Comm.hpp"

namespace Tpetra {

// ElementSpace constructor
template<typename OrdinalType, typename ScalarType>
VectorSpace<OrdinalType, ScalarType>::VectorSpace(ElementSpace<OrdinalType> const& elementSpace, 
																									Platform<OrdinalType, ScalarType> const& platform)
	: Object("Tpetra::VectorSpace")
	, blockspace_(false)
	, indexBase_(elementSpace.getIndexBase())
	, numMyEntries_(elementSpace.getNumMyElements())
	, numGlobalEntries_(elementSpace.getNumGlobalElements())
	, ElementSpace_(&elementSpace)
	, BlockElementSpace_()
	, Platform_(&platform)
	, Comm_(platform.createScalarComm())
{

}

// BlockElementSpace constructor
//template<typename OrdinalType, typename ScalarType>
//VectorSpace<OrdinalType, ScalarType>::VectorSpace(BlockElementSpace<OrdinalType> const& blockElementSpace, 
//																									Platform<OrdinalType, ScalarType> const& platform) {
//}

// copy constructor.
//template<typename OrdinalType, typename ScalarType>
//VectorSpace<OrdinalType, ScalarType>::VectorSpace(VectorSpace<OrdinalType, ScalarType> const& vectorSpace) {
//}

// destructor.
template<typename OrdinalType, typename ScalarType>
VectorSpace<OrdinalType, ScalarType>::~VectorSpace() {
}

// getLocalIndex
template<typename OrdinalType, typename ScalarType>
OrdinalType VectorSpace<OrdinalType, ScalarType>::getLocalIndex(OrdinalType globalIndex) const {
	if(!blockspace_)
		return(elementSpace().getLID(globalIndex));
}

// getGlobalIndex
template<typename OrdinalType, typename ScalarType>
OrdinalType VectorSpace<OrdinalType, ScalarType>::getGlobalIndex(OrdinalType localIndex) const {
	if(!blockspace_)
		return(elementSpace().getGID(localIndex));
}

// compatibleVector
//template<typename OrdinalType, typename ScalarType>
//bool VectorSpace<OrdinalType, ScalarType>::compatibleVector(VectorSpace<ScalarType, OrdinalType> const& vector) const {
//}

// isSameAs
//template<typename OrdinalType, typename ScalarType>
//bool VectorSpace<OrdinalType, ScalarType>::isSameAs(VectorSpace<ScalarType, OrdinalType> const& vectorSpace) const {
//	if(blockspace_)
//		return(blockElementSpace().isSameAs(vectorSpace.blockElementSpace())); // compare BlockElementSpaces
//	else
//		return(elementSpace().isSameAs(vectorSpace.elementSpace())); // compare ElementSpaces
//}

// print
/*template<typename OrdinalType, typename ScalarType>
void VectorSpace<OrdinalType, ScalarType>::print(ostream& os) const { 
  OrdinalType myImageID = comm().getMyImageID(); 
  OrdinalType numImages = comm().getNumImages(); 
  
  for (int imageCtr = 0; imageCtr < numImages; imageCtr++) { 
    if (myImageID == imageCtr) { 
      if (myImageID == 0) { 
				os << endl << "Number of Global Entries  = " << getNumGlobalEntries() << endl; 
				os <<         "Index Base                = " << getIndexBase() << endl; 
      } 

      os << endl <<   "Number of Local Entries   = " << getNumMyEntries() << endl; 
      os << endl; 
		} 
	} 
	if(blockspace_) 
		blockElementSpace.print(); 
	else 
		elementSpace().print(); 
		}*/

} // namespace Tpetra

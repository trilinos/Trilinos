/*Paul
27-July-2002 Templated for OrdinalType.
*/

#ifndef _TPETRA_SERIALDIRECTORY_H_
#define _TPETRA_SERIALDIRECTORY_H_

#include "Tpetra_Directory.h"
#include "Tpetra_Object.h"

namespace Tpetra {

// forward declaration
template<typename OrdinalType> class ElementSpace;

template<typename OrdinalType>
class SerialDirectory : public Object, public virtual Directory<OrdinalType> {
 public:
	
	// constructor
	SerialDirectory(const ElementSpace<OrdinalType>& ElementSpace);
  
	// copy constructor
 	SerialDirectory(const SerialDirectory<OrdinalType>& Directory);
  
	// destructor.
 	~SerialDirectory();
	
	// query
 	void getDirectoryEntries(OrdinalType numEntries, const OrdinalType* globalEntries, OrdinalType* procs, OrdinalType* localEntries) const;

 private:
	const ElementSpace<OrdinalType>* ElementSpace_;

}; // class SerialDirectory

} // namespace Tpetra

#include "Tpetra_SerialDirectory.cpp"
#endif // _TPETRA_SERIALDIRECTORY_H_

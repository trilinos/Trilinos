/*Paul
27-July-2002 Templated for OrdinalType.
06-August-2002 Changed to images. Documentation added.
21-Sept-2002 Comm/Platform split
*/

#ifndef _TPETRA_SERIALDIRECTORY_H_
#define _TPETRA_SERIALDIRECTORY_H_

#include "Tpetra_Directory.h"
#include "Tpetra_Object.h"

namespace Tpetra {

// forward declaration
template<class OrdinalType> class ElementSpace;

//! Tpetra::SerialDirectory: This class is a serial implementation of Directory.  Its interface allows ElementSpace and BlockElementSpace objects to reference non-local elements.

/*! For ElementSpace objects, a Directory object must be created by a call to
    the Comm createDirectory method.  The Directory is needed to allow referencing
    of non-local elements.
		For BlockElementSpace objects, a Directory should be created and used through the 
		ElementSpace accessor.

		This class currently has two constructors, one that takes an ElementSpace object, and a copy constructor.
*/

template<class OrdinalType>
class SerialDirectory : public Object, public virtual Directory<OrdinalType> {
 public:
	
	//@{ \name Constructors/Destructor.
  //! constructor
  SerialDirectory(const ElementSpace<OrdinalType>& ElementSpace);
  
  //! copy constructor
  SerialDirectory(const SerialDirectory<OrdinalType>& Directory);
  
  //! destructor.
  ~SerialDirectory();
	//@}
	
  //@{ \name Query method.
  //! getDirectoryEntries : Returns image and local id info for non-local ElementSpace entries
  /*! Given a list of Global Entry IDs, this function returns the list of
      image IDs and local IDs on the owning memory image that correspond
      to the list of entries.  If LocalEntries is 0, then local IDs are 
      not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
      element sizes for the requested global entries.
    \param In
           numEntries - Number of Global IDs being passed in.
    \param In
           globalEntries - List of Global IDs being passed in.
    \param InOut
           images - User allocated array of length at least NumEntries.  On return contains list of images
	   owning the Global IDs in question.
    \param InOut
           localEntries - User allocated array of length at least NumEntries.  On return contains the local ID of
	   the global on the owning image. If LocalEntries is zero, no local ID information is returned.
  */
  void getDirectoryEntries(OrdinalType numEntries, const OrdinalType* globalEntries, OrdinalType* images, OrdinalType* localEntries) const;
	//@}

 private:
  const ElementSpace<OrdinalType>* ElementSpace_;

}; // class SerialDirectory

} // namespace Tpetra

#include "Tpetra_SerialDirectory.cpp"
#endif // _TPETRA_SERIALDIRECTORY_H_

/*Paul
27-July-2002 Status: Templated for OrdinalType. Almost trivial SerialDirectory only implementation.
06-August-2002 Changed to images.
*/

#ifndef _TPETRA_DIRECTORY_H_
#define _TPETRA_DIRECTORY_H_

namespace Tpetra {

//! Tpetra::Directory: This class is a pure virtual class whose interface allows ElementSpace and BlockElementSpace objects to reference non-local elements.

/*! For ElementSpace objects, a Directory object must be created by a call to
    the Comm createDirectory method.  The Directory is needed to allow referencing
    of non-local elements.
*/

template<class OrdinalType>
class Directory {
  public:

  //@{ \name Constructors/Destructor.
  //! Directory destructor.
  virtual ~Directory() {};
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

  virtual void getDirectoryEntries(OrdinalType numEntries, const OrdinalType* globalEntries, OrdinalType* images, OrdinalType* localEntries) const = 0;
	
}; // class Directory

} // namespace Tpetra

#endif // _TPETRA_DIRECTORY_H_

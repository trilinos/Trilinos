/*Paul
27-July-2002 Status: Templated for OrdinalType. Almost trivial SerialDirectory only implementation.
*/

#ifndef _TPETRA_DIRECTORY_H_
#define _TPETRA_DIRECTORY_H_

namespace Tpetra {

template<typename OrdinalType>
class Directory {
  public:

  virtual ~Directory() {};

  virtual void getDirectoryEntries(OrdinalType numEntries, const OrdinalType* globalEntries, OrdinalType* procs, OrdinalType* localEntries) const = 0;
	
}; // class Directory

} // namespace Tpetra

#endif // _TPETRA_DIRECTORY_H_

#ifndef _PETRA_DIRECTORY_H_
#define _PETRA_DIRECTORY_H_

//! Petra_Directory: This class allows Petra_Map objects to reference non-local elements.

/*! For Petra_BlockMap objects, a Petra_Directory object must be created to allow referencing
    of non-local elements.  The Petra_Directory produces and contains a uniform linear
    Petra_BlockMap and a ProcList_ allowing blocks of non-local elements to be accessed
    by dereferencing throught the Petra_Directory.

  This class currently has one constructor, taking a Petra_BlockMap object.

*/

#include "Petra_Petra.h"

#ifdef DOXYGEN_SHOULD_SKIP_THIS
#include "Petra_BlockMap.h" // Doxygen needs the header files
#include "Petra_Map.h"
#else
class Petra_BlockMap;  // Compiler needs forward reference
class Petra_Map;
#endif


#ifdef PETRA_MPI
#include "GSComm_Plan.h"
#include "GSComm_Comm.h"
#endif

class Petra_Directory{
    
  public:

  //! Petra_Directory constructor
  /*! Fill this in later 
  */ 
  Petra_Directory( Petra_BlockMap* Map );
  
  //! Petra_Directory copy constructor.
  
  Petra_Directory(const Petra_Directory& Directory);
  
  //! Petra_Directory destructor.
  
  ~Petra_Directory(void);

  //! Returns the Petra_Map containing the directory
  const Petra_Map & DirectoryMap() const {return(*DirectoryMap_);};

  //! GetDirectoryEntries : Returns proc and local id info for non-local map entries
  /*! Given a list of Global Entry IDs, this function returns the list of
      processor IDs and local IDs on the owning processor that correspond
      to the list of entries.  If LocalEntries is 0, then local IDs are 
      not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
      element sizes for the requested global entries.
  */
  int GetDirectoryEntries( const int NumEntries,
			   const int * GlobalEntries,
			   int * Procs,
			   int * LocalEntries, int * EntrySizes);

 protected:

 friend class Petra_BlockMap;

 friend ostream & operator<<( ostream & os, const Petra_Directory & pd );

 private: // These need to be accessible to derived map classes.
  //! Generate: Sets up Directory tables.
  int Generate();


  Petra_BlockMap* Map_;
  Petra_Map* DirectoryMap_;

  int * ProcList_;
  int * LocalIndexList_;

  int * AllMinGIDs_;
  

};

#endif /* _PETRA_DIRECTORY_H_ */

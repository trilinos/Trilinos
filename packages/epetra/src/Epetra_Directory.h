
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_DIRECTORY_H_
#define _EPETRA_DIRECTORY_H_

#include "Epetra_Object.h"

class Epetra_BlockMap;  // Compiler needs forward reference
class Epetra_Map;

//! Epetra_Directory: This class allows Epetra_Map objects to reference non-local elements.

/*! For Epetra_BlockMap objects, a Epetra_Directory object must be created to allow referencing
    of non-local elements.  The Epetra_Directory produces and contains a uniform linear
    Epetra_BlockMap and a ProcList_ allowing blocks of non-local elements to be accessed
    by dereferencing throught the Epetra_Directory.

  This class currently has one constructor, taking a Epetra_BlockMap object.

*/

#include "Epetra_Object.h"

class Epetra_BlockMap;  // Compiler needs forward reference
class Epetra_Map;

class Epetra_Directory{
    
  public:

  //! Epetra_Directory constructor
  /*! Fill this in later 
  */ 
  Epetra_Directory( Epetra_BlockMap* Map );
  
  //! Epetra_Directory copy constructor.
  
  Epetra_Directory(const Epetra_Directory& Directory);
  
  //! Epetra_Directory destructor.
  
  ~Epetra_Directory(void);

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

  //! Print method
  virtual void Print(ostream & os) const;

 protected:

 friend class Epetra_BlockMap;

 private: // These need to be accessible to derived map classes.
  //! Generate: Sets up Directory tables.
  int Generate();

  //! Returns the Epetra_Map containing the directory
  const Epetra_Map & DirectoryMap() const {return(*DirectoryMap_);};



  Epetra_BlockMap* Map_;
  Epetra_Map* DirectoryMap_;

  int * ProcList_;
  int * LocalIndexList_;

  int * AllMinGIDs_;
  

};

#endif /* _EPETRA_DIRECTORY_H_ */

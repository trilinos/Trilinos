
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRA_BASICDIRECTORY_H
#define EPETRA_BASICDIRECTORY_H

#include "Epetra_Object.h"
#include "Epetra_Directory.h"
#include "Epetra_Map.h"

//! Epetra_BasicDirectory: This class allows Epetra_Map objects to reference non-local elements.

/*! For Epetra_BlockMap objects, a Epetra_Directory object must be created to allow referencing
    of non-local elements.  The Epetra_BasicDirectory produces and contains a uniform linear
    Epetra_BlockMap and a ProcList_ allowing blocks of non-local elements to be accessed
    by dereferencing throught the Epetra_BasicDirectory.

  This class currently has one constructor, taking a Epetra_BlockMap object.

*/

class Epetra_BasicDirectory: public virtual Epetra_Directory {
    
  public:

  //@{ \name Constructors/Destructor.
  //! Epetra_BasicDirectory constructor
  Epetra_BasicDirectory(const Epetra_BlockMap & Map );
  
  //! Epetra_BasicDirectory copy constructor.
  
  Epetra_BasicDirectory(const Epetra_BasicDirectory& Directory);
  
  //! Epetra_BasicDirectory destructor.
  
  ~Epetra_BasicDirectory(void);
  //@}
  
  //@{ \name Query method.
  //! GetDirectoryEntries : Returns proc and local id info for non-local map entries
  /*! Given a list of Global Entry IDs, this function returns the list of
      processor IDs and local IDs on the owning processor that correspond
      to the list of entries.  If LocalEntries is 0, then local IDs are 
      not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
      element sizes for the requested global entries.
    \param In
           NumEntries - Number of Global IDs being passed in.
    \param In
           GlobalEntries - List of Global IDs being passed in.
    \param InOut
           Procs - User allocated array of length at least NumEntries.  On return contains list of processors
	   owning the Global IDs in question.
    \param InOut
           LocalEntries - User allocated array of length at least NumEntries.  On return contains the local ID of
	   the global on the owning processor. If LocalEntries is zero, no local ID information is returned.
    \param InOut
           EntrySizes - User allocated array of length at least NumEntries.  On return contains the size of the
	   object associated with this global ID. If LocalEntries is zero, no size information is returned.
	   
    \return Integer error code, set to 0 if successful.
  */
  int GetDirectoryEntries( const Epetra_BlockMap& Map,
			   const int NumEntries,
			   const int * GlobalEntries,
			   int * Procs,
			   int * LocalEntries, int * EntrySizes) const;
  //@}

  //@{ \name I/O Methods.
  //! Print method
  virtual void Print(ostream & os) const;
  //@}

 private:
  Epetra_BasicDirectory& operator=(const Epetra_BasicDirectory& src)
    {
      //not currently supported
      abort();
      return( *this );
    }

 // These need to be accessible to derived map classes.
  //! Generate: Sets up Directory tables.
  int Generate(const Epetra_BlockMap& Map);

  //! Returns the Epetra_Map containing the directory
  const Epetra_Map & DirectoryMap() const {return(*DirectoryMap_);};

  Epetra_Map* DirectoryMap_;

  int * ProcList_;
  int * LocalIndexList_;
  int * SizeList_;
  bool SizeIsConst_;

  int * AllMinGIDs_;
  

};

#endif /* EPETRA_BASICDIRECTORY_H */

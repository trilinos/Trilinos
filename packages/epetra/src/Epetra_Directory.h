
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef EPETRA_DIRECTORY_H
#define EPETRA_DIRECTORY_H

class Epetra_BlockMap;  // Compiler needs forward reference
class Epetra_Map;

//! Epetra_Directory: This class is a pure virtual class whose interface allows Epetra_Map and Epetr_BlockMap objects to reference non-local elements.

/*! For Epetra_BlockMap objects, a Epetra_Directory object must be created by a call to
    the Epetra_Comm CreateDirectory method.  The Directory is needed to allow referencing
    of non-local elements.

*/
class Epetra_Directory {
    
  public:

  //@{ \name Constructors/Destructor.
  //! Epetra_Directory destructor.
  virtual ~Epetra_Directory(){};
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
  virtual int GetDirectoryEntries( const Epetra_BlockMap& Map,
				   const int NumEntries,
				   const int * GlobalEntries,
				   int * Procs,
				   int * LocalEntries,
				   int * EntrySizes) const = 0;

  //@}
};

#endif /* EPETRA_DIRECTORY_H */

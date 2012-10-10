/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_DIRECTORY_H
#define EPETRA_DIRECTORY_H

#include "Epetra_ConfigDefs.h"

class Epetra_BlockMap;  // Compiler needs forward reference
class Epetra_Map;

//! Epetra_Directory: This class is a pure virtual class whose interface allows Epetra_Map and Epetr_BlockMap objects to reference non-local elements.

/*! For Epetra_BlockMap objects, a Epetra_Directory object must be created by a call to
    the Epetra_Comm CreateDirectory method.  The Directory is needed to allow referencing
    of non-local elements.

*/
class Epetra_Directory {
    
  public:

    //! @name Constructors/Destructor
  //@{ 
  //! Epetra_Directory destructor.
  virtual ~Epetra_Directory(){}
  //@}
  
  //! @name Query method
  //@{ 
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
     
    \param In
           high_rank_sharing_procs Optional argument, defaults to true. If any GIDs appear on multiple
     processors (referred to as "sharing procs"), this specifies whether the lowest-rank proc or the
     highest-rank proc is chosen as the "owner".
     
    \return Integer error code, set to 0 if successful.
  */
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  virtual int GetDirectoryEntries( const Epetra_BlockMap& Map,
           const int NumEntries,
           const int * GlobalEntries,
           int * Procs,
           int * LocalEntries,
           int * EntrySizes,
           bool high_rank_sharing_procs=false) const = 0;
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  virtual int GetDirectoryEntries( const Epetra_BlockMap& Map,
           const int NumEntries,
           const long long * GlobalEntries,
           int * Procs,
           int * LocalEntries,
           int * EntrySizes,
           bool high_rank_sharing_procs=false) const = 0;
#endif

  //!GIDsAllUniquelyOwned: returns true if all GIDs appear on just one processor.
  /*! If any GIDs are owned by multiple processors, returns false.
   */
  virtual bool GIDsAllUniquelyOwned() const = 0;
  //@}
};

#endif /* EPETRA_DIRECTORY_H */

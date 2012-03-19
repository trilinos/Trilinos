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

#ifndef EPETRA_BASICDIRECTORY_H
#define EPETRA_BASICDIRECTORY_H

#include "Epetra_Object.h"
#include "Epetra_Directory.h"
#include "Epetra_Map.h"

//! Epetra_BasicDirectory: This class allows Epetra_Map objects to reference non-local elements.

/*! For Epetra_BlockMap objects, a Epetra_Directory object must be created to allow referencing
    of non-local elements.  The Epetra_BasicDirectory produces and contains a uniform linear
    Epetra_BlockMap and a ProcList_ allowing blocks of non-local elements to be accessed
    by dereferencing through the Epetra_BasicDirectory.

  This class currently has one constructor, taking a Epetra_BlockMap object.

*/

class Epetra_BasicDirectory: public virtual Epetra_Directory {
    
  public:

    //! @name Constructors/Destructor
  //@{ 
  //! Epetra_BasicDirectory constructor
  Epetra_BasicDirectory(const Epetra_BlockMap & Map );
  
  //! Epetra_BasicDirectory copy constructor.
  
  Epetra_BasicDirectory(const Epetra_BasicDirectory& Directory);
  
  //! Epetra_BasicDirectory destructor.
  
  virtual ~Epetra_BasicDirectory(void);
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
	   owning the Global IDs in question. If any of the GIDs is shared by more than
	   one processor, then the lowest-numbered processor is listed in this array, unless the optional
	   argument 'high_rank_sharing_procs' is given as true.
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
  int GetDirectoryEntries( const Epetra_BlockMap& Map,
			   const int NumEntries,
			   const int * GlobalEntries,
			   int * Procs,
			   int * LocalEntries,
			   int * EntrySizes,
			   bool high_rank_sharing_procs=false) const;

  int GetDirectoryEntries( const Epetra_BlockMap& Map,
			   const int NumEntries,
			   const long long * GlobalEntries,
			   int * Procs,
			   int * LocalEntries,
			   int * EntrySizes,
			   bool high_rank_sharing_procs=false) const;

  //!GIDsAllUniquelyOwned: returns true if all GIDs appear on just one processor.
  /*! If any GIDs are owned by multiple processors, returns false.
   */
  bool GIDsAllUniquelyOwned() const;
  //@}

  //! @name I/O Methods
  //@{ 
  //! Print method
  virtual void Print(ostream & os) const;
  //@}

 private:
  Epetra_BasicDirectory& operator=(const Epetra_BasicDirectory& src);

  void create_ProcListArrays();

  void addProcToList(int proc, int LID);

  //! Generate: Sets up Directory tables.
  template<typename int_type>
  int Generate(const Epetra_BlockMap& Map);

  //! Returns the Epetra_Map containing the directory
  const Epetra_Map & DirectoryMap() const {return(*DirectoryMap_);};

  Epetra_Map* DirectoryMap_;

  //ProcList_ is a list containing the associated processor for each
  //directory entry. If any directory entry has more than one associated
  //processor, then the corresponding ProcList_ entry will be the lowest-
  //numbered of those processors. In that case, refer to ProcListLists_
  //for more info.

  int * ProcList_;

  //ProcListLists_ will usually be unallocated, and set to NULL. But if
  //at least one directory entry is associcated with more than one proc,
  //then ProcListLists_ is a list of lists -- it holds, for each
  //directory-entry, a list of processors.
  //But even then, it will have a NULL list for all directory entries that
  //are associated with only one processor.
  //
  //Each list's length will be stored in ProcListLens_.
  //Example:
  //
  //if (numProcLists_ > 0) {
  //  int entry_LID = DirectoryMap_->LID(GID);
  //
  //  for(int i=0; i<ProcListLens_[entry_LID]; ++i) {
  //    cout << "entry "<<GID<<" associated with proc "
  //          <<ProcListLists_[entry_LID][i]<<endl;
  //  }
  //}
  int** ProcListLists_;
  int* ProcListLens_;
  int numProcLists_;

  //true if any directory entry appears on multiple processors
  bool entryOnMultipleProcs_;

  int * LocalIndexList_;
  int * SizeList_;
  bool SizeIsConst_;

  int * AllMinGIDs_int_;
  long long * AllMinGIDs_LL_;
  
  template<typename int_type>
  const int_type * AllMinGIDs() const;

  template<> const int * AllMinGIDs() const { return AllMinGIDs_int_; }
  template<> const long long * AllMinGIDs() const { return AllMinGIDs_LL_; }

	template<typename int_type>
	int	GetDirectoryEntries( const Epetra_BlockMap& Map,
						const int NumEntries,
						const int_type * GlobalEntries,
						int * Procs,
						int * LocalEntries,
						int * EntrySizes,
						bool high_rank_sharing_procs) const;

};

#endif /* EPETRA_BASICDIRECTORY_H */

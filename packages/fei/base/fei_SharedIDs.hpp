/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_SharedIDs_hpp_
#define _fei_SharedIDs_hpp_

#include <fei_macros.hpp>
#include <snl_fei_RaggedTable.hpp>

#include <vector>
#include <map>

namespace fei {
  /** Simple container for IDs that are shared among multiple processors. */
  class SharedIDs {
  public:
    /** constructor */
    SharedIDs()
      : sharedIDs_(-1,-1), owningProcs_() { owningProcs_.reserve(128); }

    /** Copy Constructor */
    SharedIDs(const SharedIDs& src)
      : sharedIDs_(src.sharedIDs_), owningProcs_(src.owningProcs_) {}

    /** destructor */
    virtual ~SharedIDs() {}

    /** alias for the type of the internal data container */
    typedef snl_fei::RaggedTable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> >
      table_type;

    /** Add shared IDs with specified sharing processors. */
    int addSharedID(int ID, int numSharingProcs, const int* sharingProcs)
      {
	sharedIDs_.addIndices(ID, numSharingProcs, sharingProcs);
	return(0);
      }

    /** Retrieve the internal container holding shared-ID data. */
    table_type& getSharedIDs() { return( sharedIDs_ ); }

    /** Retrieve a vector holding the owning-processors for the stored
      shared IDs. */
    std::vector<int>& getOwningProcs() { return( owningProcs_ ); }

  private:
    table_type sharedIDs_;
    std::vector<int> owningProcs_;
  };

} //namespace fei

#endif // _fei_SharedIDs_hpp_


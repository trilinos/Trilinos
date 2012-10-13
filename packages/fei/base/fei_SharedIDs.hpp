/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_SharedIDs_hpp_
#define _fei_SharedIDs_hpp_

#include <fei_macros.hpp>

#include <vector>
#include <set>
#include <map>

namespace fei {

  /** Simple container for IDs that are shared among multiple processors. */
  template<typename T>
  class SharedIDs {
  public:
    /** constructor */
    SharedIDs()
      : sharedIDs_(), owningProcs_()
    {}

    /** Copy Constructor */
    SharedIDs(const SharedIDs<T>& src)
      : sharedIDs_(src.sharedIDs_), owningProcs_(src.owningProcs_)
    {}

    /** destructor */
    virtual ~SharedIDs() {}

    /** typedef for the type of the internal data container */
    typedef std::map<T,std::set<int> > map_type;

    /** Associate ID with specified sharing processors.
     If ID is already present in the sharedIDs table, then this method simply adds
     the specified sharing processors to that ID's set of sharing processors (if
     they are not already present).
    */
    void addSharedID(const T& ID, size_t numSharingProcs, const int* sharingProcs)
    {
      typename map_type::iterator iter = sharedIDs_.find(ID);
      if (iter == sharedIDs_.end()) {
        iter = sharedIDs_.insert(std::make_pair(ID,std::set<int>())).first;
      }
      for(size_t i=0; i<numSharingProcs; ++i) {
        iter->second.insert(sharingProcs[i]);
      }
    }

    /** Retrieve the internal container holding shared-ID data. */
    map_type& getSharedIDs() { return( sharedIDs_ ); }

    /** Retrieve the internal container holding shared-ID data. */
    const map_type& getSharedIDs() const { return( sharedIDs_ ); }

    /** Retrieve a vector holding the owning-processors for the stored
      shared IDs. */
    std::vector<int>& getOwningProcs() { return( owningProcs_ ); }

    /** Retrieve a vector holding the owning-processors for the stored
      shared IDs. */
    const std::vector<int>& getOwningProcs() const { return( owningProcs_ ); }

  private:
    map_type sharedIDs_;
    std::vector<int> owningProcs_;
  };

} //namespace fei

#endif // _fei_SharedIDs_hpp_


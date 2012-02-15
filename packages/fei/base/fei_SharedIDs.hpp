/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


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


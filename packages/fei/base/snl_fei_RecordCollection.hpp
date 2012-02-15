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


#ifndef _snl_fei_RecordCollection_hpp_
#define _snl_fei_RecordCollection_hpp_

#include <fei_iosfwd.hpp>
#include <fei_IndexType.hpp>
#include <fei_Pool_alloc.hpp>
#include <fei_FieldMask.hpp>
#include <fei_Record.hpp>

#include <map>
#include <vector>

#undef fei_file
#define fei_file "snl_fei_RecordCollection.hpp"

#include <fei_ErrMacros.hpp>

namespace fei {
  template<typename T> class SharedIDs;
}

namespace snl_fei {

  /** container for Record objects */
  class RecordCollection {
  public:
    /** Constructor */
    RecordCollection(int localProc);

    /** Copy constructor */
    RecordCollection(const RecordCollection& src);

    /** Destructor */
    virtual ~RecordCollection();

    /** initialize records for specified IDs */
    void initRecords(int numIDs,
		     const int* IDs,
		     std::vector<fei::FieldMask*>& fieldMasks,
		     int* recordLocalIDs=NULL);

    /** initialize records for specified IDs with specified fieldID */
    void initRecords(int fieldID,
		    int fieldSize,
		    int numIDs,
		    const int* IDs,
		    std::vector<fei::FieldMask*>& fieldMasks,
		    int* recordLocalIDs=NULL);

    /** set owner-proc attribute for specified IDs, to be the
       lowest-numbered sharing processor */
    void setOwners_lowestSharing(fei::SharedIDs<int>& sharedIDs);

    /** Query the number of records in this collection */
    size_t getNumRecords() const
    {
      return( m_records.size() );
    }

    /** Get the native version of global-to-local ids */
    fei::IndexType<int,int>& getNativeGlobalToLocalMap()
    { return m_global_to_local; }

    const fei::IndexType<int,int>& getNativeGlobalToLocalMap() const
    { return m_global_to_local; }

  private:
    /** If someone has a reference to std::map of our data and we are using 
	something else, sync from that map to us 
     */
    void syncFrom() const {
      if (doesSomeoneHaveMyMap)
	m_global_to_local.resyncFromMap(m_global_to_local_map_);
    }

    /** If someone has a reference to std::map of our data and we are using 
	something else, sync to that map from us 
     */
    void syncTo() {
      if (doesSomeoneHaveMyMap)
	m_global_to_local.resyncToMap(m_global_to_local_map_);
    }
  public:

    /** Get the std::map version of global-to-local ids */
    FEI_DEPRECATED std::map<int,int>& getGlobalToLocalMap() const ;

    FEI_DEPRECATED const std::map<int,int>& getGlobalToLocalMap() ;

    /** Get the vector containing the records */
    std::vector<fei::Record<int> >& getRecords()
    {
      return( m_records );
    }

    /** Get the vector containing the records */
    const std::vector<fei::Record<int> >& getRecords() const
    {
      return( m_records );
    }

    /** Get record with the specified ID. Returns NULL if not found. */
    fei::Record<int>* getRecordWithID(int ID);

    /** Get record with the specified ID. Returns NULL if not found. */
    const fei::Record<int>* getRecordWithID(int ID) const;

    fei::Record<int>* getRecordWithLocalID(int lid)
    { return &m_records[lid]; }

    const fei::Record<int>* getRecordWithLocalID(int lid) const
    { return &m_records[lid]; }

    int getLocalID(int global_id) const
    {
      fei::IndexType<int,int>::const_iterator iter = m_global_to_local.find(global_id);
      if (iter == m_global_to_local.end()) {
        return -1;
      }
      return iter->second;
    }

    /** Get global equation index for specified ID */
    int getGlobalIndex(int ID,
        int fieldID,
        int fieldSize,
        int fieldOffset,
        int whichComponentOfField,
        const int* eqnNumbers);

    /** Get global block-equation index for specified ID */
    int getGlobalBlkIndex(int ID, int& globalBlkIndex);

    /** specify an output-stream for debug information */
    void setDebugOutput(FEI_OSTREAM* dbgOut)
    {
      dbgOut_ = dbgOut;
      debugOutput_ = true;
    }

  private:

    std::vector<fei::Record<int> > m_records;
    /// This has to be mutable to maintain backwards compatability when 
    /// people give out refs to private data which breaks const anyway
    mutable     fei::IndexType<int, int> m_global_to_local;
    /// this is for backwards compatability for an outdated send reference
    mutable     std::map<int, int> m_global_to_local_map_;
    /// This will cause a lot of work once this is set
    bool doesSomeoneHaveMyMap;

    int localProc_;

    bool debugOutput_;
    FEI_OSTREAM* dbgOut_;
  };

} //namespace snl_fei

#undef fei_file

#endif // _snl_fei_RecordCollection_hpp_

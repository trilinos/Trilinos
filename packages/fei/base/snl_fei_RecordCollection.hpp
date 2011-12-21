/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

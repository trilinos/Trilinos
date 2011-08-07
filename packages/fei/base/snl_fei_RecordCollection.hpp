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

    /** Get the map of global-to-local ids */
    std::map<int,int>& getGlobalToLocalMap()
    { return m_global_to_local; }

    const std::map<int,int>& getGlobalToLocalMap() const
    { return m_global_to_local; }

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
      std::map<int,int>::const_iterator iter = m_global_to_local.find(global_id);
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
    std::map<int,int> m_global_to_local;

    int localProc_;

    bool debugOutput_;
    FEI_OSTREAM* dbgOut_;
  };

} //namespace snl_fei

#undef fei_file

#endif // _snl_fei_RecordCollection_hpp_

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
#include <feiPoolAllocator.hpp>
#include <fei_FieldMask.hpp>
#include <fei_Record.hpp>

#include <map>

#undef fei_file
#define fei_file "snl_fei_RecordCollection.hpp"

#include <fei_ErrMacros.hpp>

namespace fei {
  class SharedIDs;
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

    /** alias for map container */
    typedef std::map<int,fei::Record*> map_type;

    /** initialize records for specified IDs */
    void initRecords(int numIDs,
		     const int* IDs,
		     std::vector<fei::FieldMask*>& fieldMasks,
		     fei::Record** records=NULL);

    /** initialize records for specified IDs with specified fieldID */
    void initRecords(int fieldID,
		    int fieldSize,
		    int numInstances,
		    int numIDs,
		    const int* IDs,
		    std::vector<fei::FieldMask*>& fieldMasks,
		    bool skipIDsWithThisField=true);

    /** initialize records for specified IDs with specified fieldID */
    void initRecords(int fieldID,
		    int fieldSize,
		    int numInstances,
		    int numIDs,
		    const int* IDs,
		    std::vector<fei::FieldMask*>& fieldMasks,
		    fei::Record** records,
		    bool skipIDsWithThisField=true);

    /** set owner-proc attribute for specified IDs, to be the
       lowest-numbered sharing processor */
    int setOwners_lowestSharing(fei::SharedIDs* sharedIDs);

    /** Query the number of records in this collection */
    size_t getNumRecords()
      {
	return( records_.size() );
      }

    /** Get the map containing the records */
    map_type& getRecords()
      {
	return( records_ );
      }

    /** Get record with the specified ID. Returns NULL if not found. */
    fei::Record* getRecordWithID(int ID);

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

    /** power-users only */
    void setFirstLocallyOwnedGlobalIndex(int idx)
      { firstLocallyOwnedGlobalIndex_ = idx; }

    /** power-users only */
    void setLastLocallyOwnedGlobalIndex(int idx)
      { lastLocallyOwnedGlobalIndex_ = idx; }

    /** power-users only */
    int getFirstLocallyOwnedGlobalIndex()
      { return(firstLocallyOwnedGlobalIndex_); }

    /** power-users only */
    int getLastLocallyOwnedGlobalIndex()
      { return(lastLocallyOwnedGlobalIndex_); }

  private:

    map_type records_;

    int localProc_;

    feiPoolAllocator<fei::Record > recordPool_;

    bool debugOutput_;
    FEI_OSTREAM* dbgOut_;

    int firstLocallyOwnedGlobalIndex_;
    int lastLocallyOwnedGlobalIndex_;
  };

} //namespace snl_fei

#undef fei_file

#endif // _snl_fei_RecordCollection_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_sstream.hpp"
#include "fei_fstream.hpp"

#include "fei_utils.hpp"

#include "fei_TemplateUtils.hpp"
#include "fei_chk_mpi.hpp"
#include <fei_CommUtils.hpp>
#include <fei_set_shared_ids.hpp>
#include "snl_fei_Utils.hpp"
#include "fei_Record.hpp"
#include "snl_fei_RecordCollection.hpp"
#include "fei_ParameterSet.hpp"
#include "snl_fei_RecordMsgHandler.hpp"
#include "fei_SharedIDs.hpp"
#include "fei_Pattern.hpp"
#include "fei_VectorSpace.hpp"
#include "fei_FieldMask.hpp"
#include "snl_fei_PointBlockMap.hpp"
#include "fei_LogManager.hpp"

#undef fei_file
#define fei_file "fei_VectorSpace.cpp"
#include "fei_ErrMacros.hpp"

namespace fei {
  class RecordAttributeCounter : public Record_Operator<int> {
  public:
    RecordAttributeCounter(int proc)
      : numLocalDOF_(0),
        numLocalIDs_(0),
        numLocallyOwnedIDs_(0),
        numRemoteSharedDOF_(0),
        proc_(proc)
    {}

    virtual ~RecordAttributeCounter(){}

    void operator()(fei::Record<int>& record)
    {
      fei::FieldMask* mask = record.getFieldMask();
      int owner = record.getOwnerProc();

      if (owner != proc_) {
        numRemoteSharedDOF_ += mask->getNumIndices();
        return;
      }
      else {
        ++numLocallyOwnedIDs_;
      }

      ++numLocalIDs_;

      int numDOF = mask->getNumIndices();

      numLocalDOF_ += numDOF;  
    }

    int numLocalDOF_;
    int numLocalIDs_;
    int numLocallyOwnedIDs_;
    int numRemoteSharedDOF_;

  private:
    int proc_;
  };//class Record_Operator

  class BlkIndexAccessor : public Record_Operator<int> {
  public:
    BlkIndexAccessor(int localProc,
                     int lenBlkIndices,
                     int* globalBlkIndices,
                     int* blkSizes)
      : numBlkIndices_(0),
        proc_(localProc),
        lenBlkIndices_(lenBlkIndices),
        globalBlkIndices_(globalBlkIndices),
        blkSizes_(blkSizes)
    {
    }

    BlkIndexAccessor(int lenBlkIndices,
                     int* globalBlkIndices,
                     int* blkSizes)
      : numBlkIndices_(0),
        proc_(-1),
        lenBlkIndices_(lenBlkIndices),
        globalBlkIndices_(globalBlkIndices),
        blkSizes_(blkSizes)
    {
    }

    void operator()(fei::Record<int>& record)
    {
      int owner = record.getOwnerProc();
      if (owner != proc_ && proc_ > -1) {
        return;
      }

      fei::FieldMask* mask = record.getFieldMask();
      int blkSize = mask->getNumIndices();

      if (numBlkIndices_ < lenBlkIndices_) {
        globalBlkIndices_[numBlkIndices_] = record.getNumber();
        blkSizes_[numBlkIndices_] = blkSize;
      }

      ++numBlkIndices_;
    }

    int numBlkIndices_;

  private:
    int proc_;
    int lenBlkIndices_;
    int* globalBlkIndices_;
    int* blkSizes_;
  };//class BlkIndexAccessor
}//namespace snl_fei

//----------------------------------------------------------------------------
fei::SharedPtr<fei::VectorSpace>
fei::VectorSpace::Factory::createVectorSpace(MPI_Comm comm,
                                             const char* name)
{
  fei::SharedPtr<fei::VectorSpace> sptr(new fei::VectorSpace(comm, name));
  return(sptr);
}

//----------------------------------------------------------------------------
fei::VectorSpace::VectorSpace(MPI_Comm comm, const char* name)
  : fieldMasks_(),
    comm_(comm),
    idTypes_(),
    fieldDatabase_(),
    fieldDofMap_(),
    recordCollections_(),
    sharedIDTables_(),
    ownerPatterns_(),
    sharerPatterns_(),
    sharedRecordsSynchronized_(true),
    ptBlkMap_(NULL),
    globalOffsets_(),
    globalIDOffsets_(),
    simpleProblem_(false),
    firstLocalOffset_(-1),
    lastLocalOffset_(-1),
    eqnNumbers_(),
    newInitData_(false),
    name_(),
    dbgprefix_("VecSpc: "),
    checkSharedIDs_(false)
{
  check_version();

//The initializations below are redundant with the ones above in the
//initializer list. For some reason, when compiled for janus with optimization,
//these pointers are not being set to NULL by the above initializers.
//This was causing seg faults later when other VectorSpace methods delete these
//pointers if they aren't NULL.  ABW 5/19/2004

  ptBlkMap_ = NULL;

  int numProcs = fei::numProcs(comm_);
  globalOffsets_.assign(numProcs+1, 0);
  globalIDOffsets_.assign(numProcs+1, 0);

  setName(name);
}

//----------------------------------------------------------------------------
fei::VectorSpace::~VectorSpace()
{
  int i, len = fieldMasks_.size();
  for(i=0; i<len; ++i) delete fieldMasks_[i];

  len = recordCollections_.size();
  for(i=0; i<len; ++i) delete recordCollections_[i];

  std::map<int,fei::comm_map*>::iterator
    o_iter = ownerPatterns_.begin(), o_end = ownerPatterns_.end();
  for(; o_iter != o_end; ++o_iter) {
    delete o_iter->second;
  }

  std::map<int,fei::comm_map*>::iterator
    s_iter = sharerPatterns_.begin(), s_end = sharerPatterns_.end();
  for(; s_iter != s_end; ++s_iter) {
    delete s_iter->second;
  }

  delete ptBlkMap_;
}

//----------------------------------------------------------------------------
MPI_Comm
fei::VectorSpace::getCommunicator() const
{
  return comm_;
}

//----------------------------------------------------------------------------
void fei::VectorSpace::setParameters(const fei::ParameterSet& paramset)
{
  const fei::Param* param = paramset.get("name");
  fei::Param::ParamType ptype = param != NULL ?
    param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    setName(param->getStringValue().c_str());
  }

  param = paramset.get("FEI_OUTPUT_LEVEL");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    fei::LogManager& log_manager = fei::LogManager::getLogManager();
    log_manager.setOutputLevel(param->getStringValue().c_str());
    setOutputLevel(fei::utils::string_to_output_level(param->getStringValue()));
  }

  param = paramset.get("FEI_LOG_EQN");
  ptype =  param != NULL ? param->getType() : fei::Param::BAD_TYPE;  
  if (ptype == fei::Param::INT) {
    addLogEqn(param->getIntValue());
  }

  param = paramset.get("FEI_LOG_ID");
  ptype =  param != NULL ? param->getType() : fei::Param::BAD_TYPE;  
  if (ptype == fei::Param::INT) {
    addLogID(param->getIntValue());
  }

  param = paramset.get("FEI_CHECK_SHARED_IDS");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype != fei::Param::BAD_TYPE) {
    if (ptype == fei::Param::BOOL) {
      checkSharedIDs_ = param->getBoolValue();
    }
    else if (ptype == fei::Param::INT) {
      checkSharedIDs_ = param->getIntValue() > 0 ? true : false;
    }
    else {
      checkSharedIDs_ = true;
    }
  }
  else {
    checkSharedIDs_ = false;
  }
}

//----------------------------------------------------------------------------
void fei::VectorSpace::defineFields(int numFields,
                                    const int* fieldIDs,
                                    const int* fieldSizes,
                                    const int* fieldTypes)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"defineFields ";
    for(int j=0; j<numFields; ++j) {
      os << "{"<<fieldIDs[j] << "," << fieldSizes[j] << "} ";
    }
    os << FEI_ENDL;
  }

  for (int i=0; i<numFields; ++i) {
    fieldDatabase_.insert(std::pair<int,unsigned>(fieldIDs[i], fieldSizes[i]));
    if (fieldIDs[i] >= 0) {
      if (fieldTypes != NULL) {
        fieldDofMap_.add_field(fieldIDs[i], fieldSizes[i], fieldTypes[i]);
      }
      else {
        fieldDofMap_.add_field(fieldIDs[i], fieldSizes[i]);
      }
    }
  }
}

//----------------------------------------------------------------------------
void fei::VectorSpace::defineIDTypes(int numIDTypes,
                                     const int* idTypes)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"defineIDTypes {";
    for(int j=0; j<numIDTypes; ++j) {
      os << idTypes[j] << " ";
    }
    os << "}"<<FEI_ENDL;
  }

  int localProc = fei::localProc(comm_);
  for (int i=0; i<numIDTypes; ++i) {
    int offset = fei::sortedListInsert(idTypes[i], idTypes_);

    if (offset >= 0) {
      recordCollections_.insert(recordCollections_.begin()+offset, new snl_fei::RecordCollection(localProc));
    }
  }
}

//----------------------------------------------------------------------------
int fei::VectorSpace::addDOFs(int fieldID,
                                          int numInstancesOfThisFieldPerID,
                                          int idType,
                                          int numIDs,
                                          const int* IDs)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addDOFs, fID=" << fieldID
       <<", idT="<<idType <<", ninst="
       << numInstancesOfThisFieldPerID << " {";
    for(int j=0; j<numIDs; ++j) {
      os << IDs[j] << " ";
      if (j>0 && j%20==0) os << FEI_ENDL << dbgprefix_;
    }
    os << "}"<<FEI_ENDL;
  }

  if (numIDs <= 0) return(0);

  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) ERReturn(-1);

  unsigned fieldSize = getFieldSize(fieldID);
  recordCollections_[idx]->initRecords(fieldID, fieldSize,
                                       numInstancesOfThisFieldPerID,
                                       numIDs, IDs,
                                       fieldMasks_);
  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::addDOFs(int fieldID,
                                          int numInstancesOfThisFieldPerID,
                                          int idType,
                                          int numIDs,
                                          const int* IDs,
                                          fei::Record<int>** records)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addDOFs*, fID=" << fieldID
       <<", idT="<<idType <<", ninst="
       << numInstancesOfThisFieldPerID << " {";
    for(int j=0; j<numIDs; ++j) {
      os << IDs[j] << " ";
      if (j>0 && j%20==0) os << FEI_ENDL << dbgprefix_;
    }
    os << "}"<<FEI_ENDL;
  }

  if (numIDs <= 0) return(0);

  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::VectorSpace::addDOFs: error, idType " << idType
          << " not recognized. (idTypes need to be initialized via the"
          << " method VectorSpace::defineIDTypes)";
    throw std::runtime_error(osstr.str());
  }

  unsigned fieldSize = getFieldSize(fieldID);
  recordCollections_[idx]->initRecords(fieldID, fieldSize,
                                       numInstancesOfThisFieldPerID,
                                       numIDs, IDs,
                                       fieldMasks_,
                                       records);
  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::addDOFs(int idType,
                                          int numIDs,
                                          const int* IDs)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addDOFs idT=" <<idType<<" {";
    for(int j=0; j<numIDs; ++j) {
      os << IDs[j] << " ";
      if (j>0 && j%20==0) os << FEI_ENDL << dbgprefix_;
    }
    os << "}"<<FEI_ENDL;
  }

  if (numIDs <= 0) return(0);

  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) ERReturn(-1);

  recordCollections_[idx]->initRecords(numIDs, IDs, fieldMasks_);

  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::addDOFs(int idType,
                                          int numIDs,
                                          const int* IDs,
                                          fei::Record<int>** records)
{
  if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addDOFs* idT=" <<idType<<" {";
    for(int j=0; j<numIDs; ++j) {
      os << IDs[j] << " ";
      if (j>0 && j%20==0) os << FEI_ENDL << dbgprefix_;
    }
    os << "}"<<FEI_ENDL;
  }

  if (numIDs <= 0) return(0);

  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) ERReturn(-1);

  recordCollections_[idx]->initRecords(numIDs, IDs,
                                       fieldMasks_, records);

  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::initSharedIDs(int numShared,
                                    int idType,
                                    const int* sharedIDs,
                                    const int* numSharingProcsPerID,
                                    const int* sharingProcs)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"initSharedIDs n=" <<numShared<<", idT="<<idType<< FEI_ENDL;
    int offset = 0;
    for(int ns=0; ns<numShared; ++ns) {
      os << dbgprefix_<<"#sharedID="<<sharedIDs[ns] << ", nprocs=" << numSharingProcsPerID[ns] << ", procs: ";
      for(int sp=0; sp<numSharingProcsPerID[ns]; ++sp) {
        os << sharingProcs[offset++] << " ";
      }
      os << FEI_ENDL;
    }
    os << FEI_ENDL;
  }

  if (numShared == 0) return(0);

  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) ERReturn(-1);

  fei::SharedIDs<int>& shIDs = getSharedIDs_private(idType);

  int offset = 0;
  for(int i=0; i<numShared; ++i) {
    shIDs.addSharedID(sharedIDs[i], numSharingProcsPerID[i],
                                &(sharingProcs[offset]));
    offset += numSharingProcsPerID[i];

    fei::Record<int>* rec = recordCollections_[idx]->getRecordWithID(sharedIDs[i]);
    if (rec == NULL) {
      CHK_ERR( addDOFs(idType, 1, &(sharedIDs[i])) );
    }
  }

  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::initSharedIDs(int numShared,
                                    int idType,
                                    const int* sharedIDs,
                                    const int* numSharingProcsPerID,
                                    const int* const* sharingProcs)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"initSharedIDs n=" <<numShared<<", idT="<<idType<< FEI_ENDL;
    for(int ns=0; ns<numShared; ++ns) {
      os << dbgprefix_<<"#sharedID="<<sharedIDs[ns] << ", nprocs=" << numSharingProcsPerID[ns] << ", procs: ";
      for(int sp=0; sp<numSharingProcsPerID[ns]; ++sp) {
        os << sharingProcs[ns][sp] << " ";
      }
      os << FEI_ENDL;
    }
    os << FEI_ENDL;
  }

  if (numShared == 0) return(0);

  fei::SharedIDs<int>& shIDs = getSharedIDs_private(idType);

  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) ERReturn(-1);

  for(int i=0; i<numShared; ++i) {
    shIDs.addSharedID(sharedIDs[i], numSharingProcsPerID[i],
                                 sharingProcs[i]);

    fei::Record<int>* rec = recordCollections_[idx]->getRecordWithID(sharedIDs[i]);
    if (rec == NULL) {
      CHK_ERR( addDOFs(idType, 1, &(sharedIDs[i])) );
    }
  }

  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::addVectorSpace(fei::VectorSpace* inputSpace)
{
  idTypes_ = inputSpace->idTypes_;

  std::map<int,unsigned>::const_iterator
    f_iter = inputSpace->fieldDatabase_.begin(),
    f_end = inputSpace->fieldDatabase_.end();

  for(; f_iter != f_end; ++f_iter) {
    const std::pair<const int,unsigned>& fpair = *f_iter;
    int fieldsize = fpair.second;
    defineFields(1, &(fpair.first), &fieldsize);
  }

  int i, len = inputSpace->fieldMasks_.size();
  fieldMasks_.resize(len);
  for(i=0; i<len; ++i) {
    fieldMasks_[i] = new fei::FieldMask(*(inputSpace->fieldMasks_[i]));
  }

  len = inputSpace->recordCollections_.size();
  recordCollections_.resize(len);
  for(i=0; i<len; ++i) {
    recordCollections_[i] =
      new snl_fei::RecordCollection(*(inputSpace->recordCollections_[i]));
  }

  sharedIDTables_ = inputSpace->sharedIDTables_;

  newInitData_ = true;
  sharedRecordsSynchronized_ = false;

  return(0);
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getSendProcs(std::vector<int>& sendProcs) const
{
  sendProcs.clear();
  std::set<int> sendSet;

  std::map<int,fei::SharedIDs<int> >::const_iterator
    s_it = sharedIDTables_.begin(),
    s_end= sharedIDTables_.end();
  for(; s_it!=s_end; ++s_it) {
    const std::vector<int>& owners = s_it->second.getOwningProcs();
    for(size_t i=0; i<owners.size(); ++i) {
      sendSet.insert(owners[i]);
    }
  }

  for(std::set<int>::iterator it=sendSet.begin(); it!=sendSet.end(); ++it) {
    if (fei::localProc(comm_) != *it) sendProcs.push_back(*it);
  }
}

//----------------------------------------------------------------------------
fei::SharedIDs<int>& fei::VectorSpace::getSharedIDs_private(int idType)
{
  std::map<int,fei::SharedIDs<int> >::iterator
    iter = sharedIDTables_.find(idType);

  if (iter == sharedIDTables_.end()) {
    iter = sharedIDTables_.insert(std::make_pair(idType,fei::SharedIDs<int>())).first;
  }

  return iter->second;
}

//----------------------------------------------------------------------------
void fei::VectorSpace::compute_shared_ids()
{
  //For each ID-type, call the fei::set_shared_ids function, which
  //takes a RecordCollection as input and fills a SharedIDs object
  //with mappings for which processors share each ID that appears
  //on multiple processors.
  for(size_t i=0; i<idTypes_.size(); ++i) {
    snl_fei::RecordCollection* records = recordCollections_[i];
    fei::SharedIDs<int>& shIDs = getSharedIDs_private(idTypes_[i]);

    fei::set_shared_ids(comm_, *records, shIDs);
  }
}

//----------------------------------------------------------------------------
int fei::VectorSpace::initComplete()
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os <<dbgprefix_<< "initComplete" << FEI_ENDL;
  }

  simpleProblem_ = (fieldMasks_.size()==1) && (fieldMasks_[0]->getNumFields()==1);

  //we need to know if any processor has newInitData_.
  int localInitData = newInitData_ ? 1 : 0;
  int globalInitData = 0;
  CHK_ERR( fei::GlobalMax(comm_, localInitData, globalInitData) );
  newInitData_ = globalInitData > 0 ? true : false;

  if (!newInitData_) {
    return(0);
  }

  //if running on multiple procs and initSharedIDs(...) has not been
  //called, then we need to specifically figure out which IDs are
  //shared and which procs share them.

  bool need_to_compute_shared_ids = false;
  if (fei::numProcs(comm_) > 1 && sharedIDTables_.size() < 1) {
    need_to_compute_shared_ids = true;
  }

  int localNeedToCompute = need_to_compute_shared_ids ? 1 : 0;
  int globalNeedToCompute = 0;
  CHK_ERR( fei::GlobalMin(comm_, localNeedToCompute, globalNeedToCompute) );
  need_to_compute_shared_ids = globalNeedToCompute==1 ? true : false;
  if (need_to_compute_shared_ids) {
    compute_shared_ids();
  }

  //setOwners_lowestSharing is a local operation (no communication), which
  //assumes that each processor holds CORRECT (globally symmetric)
  //shared-id/sharing-proc tables. No correctness-checking is performed here.
  setOwners_lowestSharing();

  //synchronizeSharedRecords ensures that each sharing processor has the same
  //view of the shared records, with respect to the layout of fields, which
  //determines how many DOFs and equation-numbers reside at each ID.
  //This involves inter-processor communication.
  if ( synchronizeSharedRecords() != 0) {
    return(-1);
  }

  //calculateGlobalIndices is also a global operation.
  if ( calculateGlobalIndices() != 0) {
    return(-1);
  }

  //finally we need to exchange global indices for shared records. i.e., the
  //processors that own shared records, need to send the global indices for
  //those records to the sharing-but-not-owning processors.
  if (fei::numProcs(comm_) > 1) {
    if ( exchangeGlobalIndices() != 0) {
      return(-1);
    }
  }

  newInitData_ = false;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalIndex(int idType,
                                     int ID,
                                     int fieldID,
                                     int fieldOffset,
                                     int whichComponentOfField,
                                     int& globalIndex)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(-1);

  unsigned fieldSize = 0;
  if (fieldOffset > 0) {
    fieldSize = getFieldSize(fieldID);
  }

  try {
    globalIndex = recordCollections_[idindex]->getGlobalIndex(ID,
                                                         fieldID,
                                                         fieldSize,
                                                         fieldOffset,
                                                         whichComponentOfField,
                                                         &eqnNumbers_[0]);
  }
  catch (std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "VectorSpace::getGlobalIndex caught exception: " << exc.what();
    fei::console_out() << osstr.str()<<FEI_ENDL;
    ERReturn(-1);
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalIndex(int idType,
                                     int ID,
                                     int fieldID,
                                     int& globalIndex)
{
  return( getGlobalIndex(idType, ID, fieldID, 0, 0, globalIndex) );
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalBlkIndex(int idType,
                                        int ID,
                                        int& globalBlkIndex)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(-1);

  CHK_ERR(recordCollections_[idindex]->getGlobalBlkIndex(ID, globalBlkIndex));

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalIndices(int numIDs,
                                       const int* IDs,
                                       int idType,
                                       int fieldID,
                                       int* globalIndices)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(-1);

  unsigned fieldSize = getFieldSize(fieldID);
  int offset = 0;

  for(int i=0; i<numIDs; ++i) {
    try {
      globalIndices[offset] =
        recordCollections_[idindex]->getGlobalIndex(IDs[i],
                                                    fieldID,
                                                    fieldSize,
                                                    0, 0,
                                                    &eqnNumbers_[0]);
      if (fieldSize>1) {
        int eqn = globalIndices[offset];
        for(unsigned j=1; j<fieldSize; ++j) {
          globalIndices[offset+j] = eqn+j;
        }
      }
    }
    catch (...) {
      for(unsigned j=0; j<fieldSize; ++j) {
        globalIndices[offset+j] = -1;
      }
    }

    offset += fieldSize;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalBlkIndices(int numIDs,
                                          const int* IDs,
                                          int idType,
                                          int* globalBlkIndices)
{
  int err;
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(-1);

  int offset = 0;

  for(int i=0; i<numIDs; ++i) {
    err = recordCollections_[idindex]->getGlobalBlkIndex(IDs[i],
                                                     globalBlkIndices[offset]);
    if (err != 0) {
      globalBlkIndices[offset] = -1;
    }
    ++offset;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalIndices(int numIDs,
                                       const int* IDs,
                                       const int* idTypes,
                                       const int* fieldIDs,
                                       int* globalIndices)
{
  int err;
  int offset = 0;
  for(int i=0; i<numIDs; ++i) {
    unsigned fieldSize = getFieldSize(fieldIDs[i]);
    err = getGlobalIndex(idTypes[i], IDs[i], fieldIDs[i], 0, 0,
                         globalIndices[offset]);
    if (err) {
      for(unsigned j=1; j<fieldSize; ++j) {
        globalIndices[offset+j] = -1;
      }
    }
    else {
      if (fieldSize>1) {
        int eqn = globalIndices[offset];
        for(unsigned j=1; j<fieldSize; ++j) {
          globalIndices[offset+j] = eqn+j;
        }
      }
    }
    offset += fieldSize;
  }

  return(0);
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalBlkIndices(const fei::Pattern* pattern,
                                           const fei::Record<int>*const* records,
                                           std::vector<int>& indices)
{
  int numRecords = pattern->getNumIDs();
  indices.resize(numRecords);
  int numIndices;
  getGlobalBlkIndices(numRecords, records, numRecords, &indices[0],
                      numIndices);
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalIndices(const fei::Pattern* pattern,
                                        const fei::Record<int>*const* records,
                                        std::vector<int>& indices)
{
  int numRecords = pattern->getNumIDs();
  int numIndices = pattern->getNumIndices();
  indices.resize(numIndices);
  int* indices_ptr = &indices[0];

  fei::Pattern::PatternType pType = pattern->getPatternType();

  if (pType == fei::Pattern::GENERAL ||
      pType == fei::Pattern::SINGLE_IDTYPE) {
    const int* numFieldsPerID = pattern->getNumFieldsPerID();
    const int* fieldIDs = pattern->getFieldIDs();
    int totalNumFields = pattern->getTotalNumFields();

    std::vector<int> fieldSizes(totalNumFields);

    for(int j=0; j<totalNumFields; ++j) {
      fieldSizes[j] = getFieldSize(fieldIDs[j]);
    }

    getGlobalIndices(numRecords, records, numFieldsPerID,
                     fieldIDs, &(fieldSizes[0]),
                     numIndices, indices_ptr, numIndices);
  }
  else if (pType == fei::Pattern::SIMPLE) {
    const int* fieldIDs = pattern->getFieldIDs();

    int fieldID = fieldIDs[0];
    unsigned fieldSize = getFieldSize(fieldID);

    getGlobalIndices(numRecords, records,
                     fieldID, fieldSize,
                     numIndices, indices_ptr, numIndices);
  }
  else if (pType == fei::Pattern::NO_FIELD) {
    getGlobalBlkIndices(numRecords, records, numIndices, indices_ptr, numIndices);
  }
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalIndices(int numRecords,
                                        const fei::Record<int>*const* records,
                                        int fieldID,
                                        int fieldSize,
                                        int indicesAllocLen,
                                        int* indices,
                                        int& numIndices)
{
  numIndices = 0;
  int* eqnPtr = &eqnNumbers_[0];

  int len = numRecords;
  if (len*fieldSize >= indicesAllocLen) {
    len = indicesAllocLen/fieldSize;
  }

  if (fieldSize == 1 && simpleProblem_) {
    for(int i=0; i<len; ++i) {
      const fei::Record<int>* record = records[i];
      indices[numIndices++] = *(eqnPtr+record->getOffsetIntoEqnNumbers());
    }    
    return;
  }

  if (fieldSize == 1) {
    for(int i=0; i<len; ++i) {
      const fei::Record<int>* record = records[i];

      int eqnOffset = 0, numInstances = 0;
      int* eqnNumbers = eqnPtr+record->getOffsetIntoEqnNumbers();

      const fei::FieldMask* fieldMask = record->getFieldMask();
      fieldMask->getFieldEqnOffset(fieldID, eqnOffset, numInstances);
      indices[numIndices++] = eqnNumbers[eqnOffset];
    }
  }
  else {
    for(int i=0; i<len; ++i) {
      const fei::Record<int>* record = records[i];

      int* eqnNumbers = eqnPtr+record->getOffsetIntoEqnNumbers();

      int eqnOffset = 0, numInstances = 0;
      if (!simpleProblem_) {
        const fei::FieldMask* fieldMask = record->getFieldMask();
        fieldMask->getFieldEqnOffset(fieldID, eqnOffset, numInstances);
      }

      for(int fs=0; fs<fieldSize; ++fs) {
        indices[numIndices++] = eqnNumbers[eqnOffset+fs];
      }
    }
  }
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalIndices(int numRecords,
                                        const fei::Record<int>*const* records,
                                        const int* numFieldsPerID,
                                        const int* fieldIDs,
                                        const int* fieldSizes,
                                        int indicesAllocLen,
                                        int* indices,
                                        int& numIndices)
{
  numIndices = 0;
  int fld_offset = 0;
  int* eqnPtr = &eqnNumbers_[0];

  for(int i=0; i<numRecords; ++i) {
    const fei::Record<int>* record = records[i];

    const fei::FieldMask* fieldMask = record->getFieldMask();
    int* eqnNumbers = eqnPtr + record->getOffsetIntoEqnNumbers();

    for(int nf=0; nf<numFieldsPerID[i]; ++nf) {
      int eqnOffset = 0, numInstances = 0;
      if (!simpleProblem_) {
        fieldMask->getFieldEqnOffset(fieldIDs[fld_offset], eqnOffset, numInstances);
      }

      for(int fs=0; fs<fieldSizes[fld_offset]; ++fs) {
        indices[numIndices++] = eqnNumbers[eqnOffset+fs];
      }

      ++fld_offset;
    }
  }
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalBlkIndices(int numRecords,
                                           const fei::Record<int>*const* records,
                                           int indicesAllocLen,
                                           int* indices,
                                           int& numIndices)
{
  numIndices = 0;
  for(int i=0; i<numRecords; ++i) {
    if (numIndices < indicesAllocLen) {
      indices[numIndices++] = records[i]->getNumber();
    }
    else ++numIndices;
  }
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalIndex(int idType,
                                     int ID,
                                     int& globalIndex)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(-1);

  fei::Record<int>* record = recordCollections_[idindex]->getRecordWithID(ID);
  if (record == NULL) {
    ERReturn(-1);
  }

  const int* eqnNums = &eqnNumbers_[0]
                     + record->getOffsetIntoEqnNumbers();

  if (eqnNums != NULL) { globalIndex = eqnNums[0]; return(0); }
  return(-1);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumDegreesOfFreedom(int idType,
                                             int ID)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(0);

  fei::Record<int>* record = recordCollections_[idindex]->getRecordWithID(ID);
  if (record == NULL) {
    ERReturn(-1);
  }

  return( record->getFieldMask()->getNumIndices() );
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumFields()
{
  return(fieldDatabase_.size());
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getFields(std::vector<int>& fieldIDs)
{
  unsigned numFields = fieldDatabase_.size();

  fieldIDs.resize(numFields);

  fei::copyKeysToArray<std::map<int,unsigned> >(fieldDatabase_, numFields,
                                                    &fieldIDs[0]);
}

//----------------------------------------------------------------------------
size_t fei::VectorSpace::getNumIDTypes()
{
  return(idTypes_.size());
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getIDTypes(std::vector<int>& idTypes) const
{
  size_t numIDTypes = idTypes_.size();
  idTypes.resize(numIDTypes);
  for(size_t i=0; i<numIDTypes; ++i) idTypes[i] = idTypes_[i];
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumFields(int idType,
                                   int ID)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(0);

  fei::Record<int>* record = recordCollections_[idindex]->getRecordWithID(ID);
  if (record == NULL) {
    ERReturn(-1);
  }

  return( record->getFieldMask()->getNumFields() );
}

//----------------------------------------------------------------------------
bool fei::VectorSpace::isLocal(int idType, int ID)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(false);

  fei::Record<int>* record = recordCollections_[idindex]->getRecordWithID(ID);
  if (record == NULL) {
    return(false);
  }

  return(true);
}

//----------------------------------------------------------------------------
bool fei::VectorSpace::isLocallyOwned(int idType, int ID)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) return(false);

  fei::Record<int>* record = recordCollections_[idindex]->getRecordWithID(ID);
  if (record == NULL) {
    return(false);
  }
  if (record->getOwnerProc() == fei::localProc(comm_)) {
    return(true);
  }

  return(false);
}

//----------------------------------------------------------------------------
unsigned fei::VectorSpace::getFieldSize(int fieldID)
{
  std::map<int,unsigned>::const_iterator
    f_iter = fieldDatabase_.find(fieldID);

  if (f_iter == fieldDatabase_.end()) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::VectorSpace";
    if (name_.length() > 0) {
      osstr << "(name: "<<name_<<")";
    }
    osstr << "::getFieldSize: fieldID " << fieldID << " not found.";
    throw std::runtime_error(osstr.str());
  }

  return((*f_iter).second);
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getFields(int idType, int ID, std::vector<int>& fieldIDs)
{
  int idindex = fei::binarySearch(idType, idTypes_);
  if (idindex < 0) {
    fieldIDs.resize(0);
    return;
  }

  fei::Record<int>* record = recordCollections_[idindex]->getRecordWithID(ID);
  if (record == NULL) {
    voidERReturn;
  }

  fei::FieldMask* fieldMask = record->getFieldMask();
  std::vector<int>& maskFieldIDs = fieldMask->getFieldIDs();
  int numFields = maskFieldIDs.size();
  fieldIDs.resize(numFields);
  for(int i=0; i<numFields; ++i) {
    fieldIDs[i] = maskFieldIDs[i];
  }
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalIndexOffsets(std::vector<int>& globalOffsets) const
{
  globalOffsets = globalOffsets_;
}

//----------------------------------------------------------------------------
void fei::VectorSpace::getGlobalBlkIndexOffsets(std::vector<int>& globalBlkOffsets) const
{
  globalBlkOffsets = globalIDOffsets_;
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getOwnerProcPtIndex(int globalIndex)
{
  if (globalIndex < 0) return(-1);

  unsigned len = globalOffsets_.size();
  for(int i=0; i<(int)(len-1); ++i) {
    if (globalIndex < globalOffsets_[i+1]) {
      return(i);
    }
  }

  return(-1);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getOwnerProcBlkIndex(int globalIndex)
{
  if (globalIndex < 0) return(-1);

  unsigned len = globalOffsets_.size();
  for(int i=0; i<(int)(len-1); ++i) {
    if (globalIndex < globalIDOffsets_[i+1]) {
      return(i);
    }
  }

  return(-1);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumOwnedAndSharedIDs(int idType)
{
  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) return(0);

  return( recordCollections_[idx]->getNumRecords() );
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumOwnedIDs(int idType)
{
  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) return(0);

  fei::RecordAttributeCounter attrCounter(fei::localProc(comm_));
  runRecords(attrCounter);

  return( attrCounter.numLocallyOwnedIDs_ );
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumSharedIDs(int idType, int& numShared)
{
  numShared = sharedIDTables_[idType].getSharedIDs().size();

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getOwnedAndSharedIDs(int idType,
                                      int lenList,
                                      int* IDs,
                                      int& numLocalIDs)
{
  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) return(-1);

  snl_fei::RecordCollection* records = recordCollections_[idx];

  snl_fei::RecordCollection::map_type& rmap = records->getRecords();

  numLocalIDs = rmap.size();

  fei::copyKeysToArray(rmap, lenList, IDs);

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getOwnedIDs(int idType,
                                         int lenList,
                                         int* IDs,
                                         int& numLocalIDs)
{
  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) return(-1);

  snl_fei::RecordCollection* records = recordCollections_[idx];

  snl_fei::RecordCollection::map_type& rmap = records->getRecords();

  numLocalIDs = 0;

  snl_fei::RecordCollection::map_type::iterator
    r_iter = rmap.begin(),
    r_end  = rmap.end();

  for(; r_iter != r_end; ++r_iter) {
    fei::Record<int>& thisrecord = *((*r_iter).second);

    if (thisrecord.getOwnerProc() == fei::localProc(comm_)) {
      if (numLocalIDs < lenList) {
        IDs[numLocalIDs] = thisrecord.getID();
      }
      ++numLocalIDs;
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumIndices_SharedAndOwned() const
{
  return(eqnNumbers_.size());
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getIndices_SharedAndOwned(std::vector<int>& globalIndices) const
{
  if (eqnNumbers_.size() == 0) {
    globalIndices.resize(0);
    return(0);
  }

  size_t numIndices = eqnNumbers_.size();
  const int* indicesPtr = &eqnNumbers_[0];

  globalIndices.resize(numIndices);
  for(size_t i=0; i<numIndices; ++i) {
    globalIndices[i] = indicesPtr[i];
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumBlkIndices_SharedAndOwned(int& numBlkIndices) const
{
  numBlkIndices = 0;
  for(size_t i=0; i<recordCollections_.size(); ++i) {
    numBlkIndices += recordCollections_[i]->getNumRecords();
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getBlkIndices_SharedAndOwned(int lenBlkIndices,
                                                   int* globalBlkIndices,
                                                   int* blkSizes,
                                                   int& numBlkIndices)
{
  if (!sharedRecordsSynchronized_) {
    numBlkIndices = 0;
    return(-1);
  }

  fei::BlkIndexAccessor blkIndAccessor(lenBlkIndices,
                                  globalBlkIndices, blkSizes);
  runRecords(blkIndAccessor);

  numBlkIndices = blkIndAccessor.numBlkIndices_;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalNumIndices() const
{
  if (globalOffsets_.size() < 1) return(0);
  return(globalOffsets_[globalOffsets_.size()-1]);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumIndices_Owned() const
{
  if (!sharedRecordsSynchronized_) {
    return(-1);
  }

  int localProc = fei::localProc(comm_);
  int numIndices = globalOffsets_[localProc+1]-globalOffsets_[localProc];

  return(numIndices);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getIndices_Owned(std::vector<int>& globalIndices) const
{
  if (!sharedRecordsSynchronized_) {
    globalIndices.clear();
    return(-1);
  }

  int localProc = fei::localProc(comm_);
  int numIndices = globalOffsets_[localProc+1]-globalOffsets_[localProc];

  globalIndices.resize(numIndices);

  int firstOffset = globalOffsets_[localProc];
  for(int i=0; i<numIndices; ++i) {
    globalIndices[i] = firstOffset+i;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getIndices_Owned(int lenIndices, int* globalIndices, int& numIndices) const
{
  if (!sharedRecordsSynchronized_) {
    numIndices = 0;
    return(-1);
  }

  int localProc = fei::localProc(comm_);
  numIndices = globalOffsets_[localProc+1]-globalOffsets_[localProc];

  int len = lenIndices >= numIndices ? numIndices : lenIndices;

  int firstOffset = globalOffsets_[localProc];
  for(int i=0; i<len; ++i) {
    globalIndices[i] = firstOffset+i;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getNumBlkIndices_Owned() const
{
  if (!sharedRecordsSynchronized_) {
    return(-1);
  }

  int localProc = fei::localProc(comm_);
  int numBlkIndices = globalIDOffsets_[localProc+1]-globalIDOffsets_[localProc];

  return(numBlkIndices);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getGlobalNumBlkIndices() const
{
  int numBlkIndices = 0;
  unsigned len = globalIDOffsets_.size();
  if (len > 0) {
    numBlkIndices = globalIDOffsets_[len-1];
  }

  return(numBlkIndices);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getBlkIndices_Owned(int lenBlkIndices,
                                          int* globalBlkIndices,
                                          int* blkSizes,
                                          int& numBlkIndices)
{
  if (!sharedRecordsSynchronized_) {
    numBlkIndices = 0;
    return(-1);
  }

  int localProc = fei::localProc(comm_);
  fei::BlkIndexAccessor blkIndAccessor(localProc, lenBlkIndices,
                                           globalBlkIndices, blkSizes);
  runRecords(blkIndAccessor);

  numBlkIndices = blkIndAccessor.numBlkIndices_;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getRecordCollection(int idType,
                                          snl_fei::RecordCollection*& records)
{
  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) return(-1);

  records = recordCollections_[idx];
  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::getRecordCollection(int idType,
                                          const snl_fei::RecordCollection*& records) const
{
  int idx = fei::binarySearch(idType, idTypes_);
  if (idx < 0) return(-1);

  records = recordCollections_[idx];
  return(0);
}

//----------------------------------------------------------------------------
void fei::VectorSpace::setOwners_lowestSharing()
{
  //first, add localProc to each of the sharing-proc lists, in case it wasn't
  //included via initSharedIDs().
  int localProc = fei::localProc(comm_);

  std::map<int, fei::SharedIDs<int> >::iterator
    s_iter = sharedIDTables_.begin(), s_end = sharedIDTables_.end();

  for(; s_iter != s_end; ++s_iter) {
    fei::SharedIDs<int>::map_type& shid_table = s_iter->second.getSharedIDs();

    fei::SharedIDs<int>::map_type::iterator
      t_iter = shid_table.begin(),
      t_end = shid_table.end();

    for(; t_iter != t_end; ++t_iter) {
      std::set<int>& shProcs = t_iter->second;
      shProcs.insert(localProc);
    }
  }

  //now set the owningProcs for the SharedIDs records, and the owning procs on
  //the appropriate records in the recordCollections. Set the owner to be the
  //lowest-numbered sharing proc in all cases.
  s_iter = sharedIDTables_.begin();
  for(; s_iter != s_end; ++s_iter) {
    fei::SharedIDs<int>::map_type& shid_table = s_iter->second.getSharedIDs();

    std::vector<int>& owningProcs = s_iter->second.getOwningProcs();

    int len = shid_table.size();
    owningProcs.resize(len);
    int j = 0;

    fei::SharedIDs<int>::map_type::iterator
      t_iter = shid_table.begin(),
      t_end = shid_table.end();

    for(; t_iter != t_end; ++t_iter, ++j) {
      std::set<int>& shProcs = (*t_iter).second;
      int lowest = *(shProcs.begin());
      owningProcs[j] = lowest;
    }
    
    int idx = fei::binarySearch(s_iter->first, idTypes_);
    if (idx < 0) {
      throw std::runtime_error("internal error in fei::VectorSapce::setOwners_lowestSharing");
    }

    if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
      recordCollections_[idx]->setDebugOutput(output_stream_);
    }

    recordCollections_[idx]->setOwners_lowestSharing(s_iter->second);
  }
}

//----------------------------------------------------------------------------
int fei::VectorSpace::calculateGlobalIndices()
{
  int localProc = fei::localProc(comm_);
  int numProcs = fei::numProcs(comm_);
  std::vector<int> localOffsets(numProcs +1, 0);
  std::vector<int> localIDOffsets(numProcs+1, 0);
  globalOffsets_.resize(numProcs + 1, 0);

  globalIDOffsets_.assign(globalIDOffsets_.size(), 0);

  //first we'll calculate the number of local degrees of freedom, and the
  //number of local identifiers.
  fei::RecordAttributeCounter counter(localProc);
  runRecords(counter);

  int numLocalDOF = counter.numLocalDOF_;
  int numLocalIDs = counter.numLocalIDs_;
  int numRemoteSharedDOF = counter.numRemoteSharedDOF_;

  eqnNumbers_.resize(numLocalDOF+numRemoteSharedDOF);

  localOffsets[localProc] = numLocalDOF;
  CHK_MPI( fei::GlobalMax(comm_, localOffsets, globalOffsets_) );

  localIDOffsets[localProc] = numLocalIDs;
  CHK_MPI( fei::GlobalMax(comm_, localIDOffsets, globalIDOffsets_) );

  //Now globalOffsets_ contains numLocalDOF for proc i, in the i-th position.
  //(and similarly for globalIDOffsets_)
  //So all we need to do is turn that data into global-offsets (i.e., the
  //starting global-offset for each processor).
  int localOffset = 0;
  int localIDOffset = 0;
  for(int p=0; p<numProcs; ++p) {
    numLocalDOF = globalOffsets_[p];
    globalOffsets_[p] = localOffset;
    localOffset += numLocalDOF;
    numLocalIDs = globalIDOffsets_[p];
    globalIDOffsets_[p] = localIDOffset;
    localIDOffset += numLocalIDs;
  }
  globalOffsets_[numProcs] = localOffset;
  globalIDOffsets_[numProcs] = localIDOffset;

  firstLocalOffset_ = globalOffsets_[localProc];
  lastLocalOffset_ = globalOffsets_[localProc+1] - 1;

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os <<dbgprefix_<< "  firstLocalOffset_: " << firstLocalOffset_ << ", lastLocalOffset_: "
      << lastLocalOffset_ << FEI_ENDL;
  }

  //Now we're ready to set the equation-numbers on all local records.
  CHK_ERR( setLocalEqnNumbers() );

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::synchronizeSharedRecords()
{
  if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os<<dbgprefix_<<"#synchronizeSharedRecords num-field-masks: "<<fieldMasks_.size()<<FEI_ENDL;
    for(unsigned fm=0; fm<fieldMasks_.size(); ++fm) {
      os << dbgprefix_<<"#     maskID["<<fm<<"]: " << fieldMasks_[fm]->getMaskID() << FEI_ENDL;
    }
  }

  bool safetyCheck = checkSharedIDs_;

  int numProcs = fei::numProcs(comm_);
  int localProc = fei::localProc(comm_);
  int numShTables = sharedIDTables_.size();

  if (numProcs < 2) {
    sharedRecordsSynchronized_ = true;
    return(0);
  }

  if (safetyCheck) {
    if (localProc == 0 && output_level_ >= fei::BRIEF_LOGS) {
      FEI_COUT << "fei::VectorSpace: global consistency-check of shared ID"
           << " data (involves all-to-all communication). This is done "
           << "only if requested by parameter 'FEI_CHECK_SHARED_IDS true'."
           << FEI_ENDL;
    }

    int totalNumShTables = 0;
    CHK_ERR( fei::GlobalSum(comm_, numShTables, totalNumShTables) );
    if (totalNumShTables != numShTables * numProcs) {
      //not all processors have the same number of shared-id tables. This means
      //that one or more processors is 'empty', which may not be an error. But
      //it does mean that the safety check can't be performed because the safety
      //check involves all-to-all communication and if one or more processors
      //don't participate then we'll hang...
      safetyCheck = false;
    }
  }

  //Create a list of fei::comm_maps which will be the communication-patterns
  //for each of the sharedIDTables. A sharedIDTable maps lists of processors
  //to each shared ID. The communication-pattern will be a transpose of
  //that, mapping lists of IDs to owning or sharing processors.
  //

  int local_err = 0;

  for(size_t i=0; i<idTypes_.size(); ++i) {

    fei::SharedIDs<int>& shared = getSharedIDs_private(idTypes_[i]);
    
    //now create/initialize ownerPatterns which map owning processors to lists
    //of ids that are shared locally, and sharerPatterns which map sharing
    //processors to lists of ids that are owned locally.
    fei::comm_map* ownerPattern = new fei::comm_map(1, numProcs);

    std::map<int,fei::comm_map*>::iterator iter = ownerPatterns_.find(idTypes_[i]);
    if (iter == ownerPatterns_.end()) {
      ownerPatterns_.insert(std::make_pair(idTypes_[i], ownerPattern));
    }
    else {
      delete iter->second;
      iter->second = ownerPattern;
    }

    fei::comm_map* sharerPattern = new fei::comm_map(1, numProcs);

    iter = sharerPatterns_.find(idTypes_[i]);
    if (iter == sharerPatterns_.end()) {
      sharerPatterns_.insert(std::make_pair(idTypes_[i], sharerPattern));
    }
    else {
      delete iter->second;
      iter->second = sharerPattern;
    }

    std::vector<int>& owningProcs = shared.getOwningProcs();

    fei::SharedIDs<int>::map_type& shtable = shared.getSharedIDs();

    fei::SharedIDs<int>::map_type::iterator
      s_iter = shtable.begin(),
      s_end = shtable.end();
    
    int j = 0;
    for(; s_iter != s_end; ++s_iter, ++j) {
      int ID = (*s_iter).first;
      std::set<int>& shProcs = (*s_iter).second;

      int owner = owningProcs[j];

      if (owner == localProc) {
        std::set<int>::const_iterator
          p_iter = shProcs.begin(),
          p_end = shProcs.end();
        for(; p_iter != p_end; ++p_iter) {
          if (*p_iter != localProc) {
            sharerPattern->addIndices(*p_iter, 1, &ID);
          }
        }
      }
      else {
        ownerPattern->addIndices(owner, 1, &ID);
      }
    }

    if (safetyCheck) {
      fei::comm_map* checkPattern = NULL;
      CHK_ERR( fei::mirrorCommPattern(comm_, ownerPattern, checkPattern) );
      fei::Barrier(comm_);
      if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
        FEI_OSTREAM& os = *output_stream_;
        fei::comm_map::map_type& owner_map = ownerPattern->getMap();
        int numKeys = owner_map.size();
        fei::comm_map::map_type::const_iterator
          omap_iter = owner_map.begin(),
          omap_end = owner_map.end();

        os << dbgprefix_<<"#  synchronizeSharedRecords" << FEI_ENDL
           << dbgprefix_<<"#  ownerPattern, num-procs-to-send-to: " << numKeys << FEI_ENDL;
        for(int sk=0; omap_iter != omap_end; ++sk, ++omap_iter) {
          os << dbgprefix_<<"#    sendProc["<<sk<<"]: " << omap_iter->first << " IDs: ";
          fei::comm_map::row_type::const_iterator
            val_iter = omap_iter->second->begin(),
            val_end =  omap_iter->second->end();
          for(; val_iter != val_end; ++val_iter) {
            os << *val_iter << " ";
          }
          os << FEI_ENDL;
        }

        fei::comm_map::map_type& check_map = checkPattern->getMap();
        int numCheckKeys = check_map.size();
        fei::comm_map::map_type::const_iterator
          cmap_iter = check_map.begin(),
          cmap_end = check_map.end();

        os <<dbgprefix_<< "#  synchronizeSharedRecords" << FEI_ENDL
           <<dbgprefix_<< "#  checkPattern (send mirror), num-procs: "
           << numCheckKeys << FEI_ENDL;
        for(int sk=0; cmap_iter != cmap_end; ++sk, ++cmap_iter) {
          os <<dbgprefix_<< "#    proc["<<sk<<"]: " << cmap_iter->first << " IDs: ";
          fei::comm_map::row_type::const_iterator
          val_iter = cmap_iter->second->begin(),
            val_end = cmap_iter->second->end();
          for(; val_iter != val_end; ++val_iter) {
            os << *val_iter << " ";
          }
          os << FEI_ENDL;
        }
      }

      int err = 0;
      bool quiet = false;
      if (!checkPattern->equal(*sharerPattern, quiet)) {
        //fei::console_out() << "shared-ID safety-check failed..."
        //     << FEI_ENDL;
        err = -1;
      }
      int globalErr = 0;
      CHK_ERR( fei::GlobalSum(comm_, err, globalErr) );

      delete checkPattern;

      if (globalErr != 0) {
        return(globalErr);
      }
    }

    local_err += exchangeFieldInfo(ownerPattern, sharerPattern,
                                   recordCollections_[i], fieldMasks_);
  }

  int global_err = 0;
  CHK_ERR( fei::GlobalSum(comm_, local_err, global_err) );

  if (global_err != 0) {
    ERReturn(-1);
  }

  sharedRecordsSynchronized_ = true;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::exchangeGlobalIndices()
{
  for(size_t i=0; i<idTypes_.size(); ++i) {
    snl_fei::RecordMsgHandler recmsghndlr(fei::localProc(comm_),
                                 recordCollections_[i], *ptBlkMap_,
                                          fieldMasks_, eqnNumbers_);
    recmsghndlr.setTask(snl_fei::RecordMsgHandler::_EqnNumbers_);

    recmsghndlr.setSendPattern(sharerPatterns_[idTypes_[i]]);
    recmsghndlr.setRecvPattern(ownerPatterns_[idTypes_[i]]);
    CHK_ERR( fei::exchange(comm_, &recmsghndlr) );
  }

  return(0);
}

//----------------------------------------------------------------------------
void fei::VectorSpace::runRecords(fei::Record_Operator<int>& record_op)
{
  for(size_t rec=0; rec<recordCollections_.size(); ++rec) {
    snl_fei::RecordCollection* records = recordCollections_[rec];

    snl_fei::RecordCollection::map_type& rmap = records->getRecords();

    snl_fei::RecordCollection::map_type::iterator
      r_iter = rmap.begin(),
      r_end  = rmap.end();

    for(; r_iter != r_end; ++r_iter) {
      fei::Record<int>& thisrecord = *((*r_iter).second);

      record_op(thisrecord);
    }
  }
}

//----------------------------------------------------------------------------
int fei::VectorSpace::setLocalEqnNumbers()
{
  int proc = fei::localProc(comm_);
  int eqnNumber = firstLocalOffset_;
  int idNumber = globalIDOffsets_[proc];

  int numProcs = fei::numProcs(comm_);
  int localProc = fei::localProc(comm_);

  if (ptBlkMap_ != NULL) {
    delete ptBlkMap_;
  }
  ptBlkMap_ = new snl_fei::PointBlockMap;

  int maxNumIndices = 0;
  for(unsigned i=0; i<fieldMasks_.size(); ++i) {
    if (fieldMasks_[i]->getNumIndices() > maxNumIndices) {
      maxNumIndices = fieldMasks_[i]->getNumIndices();
    }
  }

  if (maxNumIndices == 1) {
    ptBlkMap_->setPtEqualBlk();
  }

  FEI_OSTREAM* id2eqnStream = NULL;
  if (output_level_ >= fei::BRIEF_LOGS) {
    std::string path = fei::LogManager::getLogManager().getOutputPath();
    if (path == "") path = ".";
    FEI_OSTRINGSTREAM osstr;
    osstr << path << "/fei_ID2Eqn";
    if (name_.size() > 0) osstr << "_" << name_;
    osstr << "." <<numProcs<< "." << localProc;

    id2eqnStream = new FEI_OFSTREAM(osstr.str().c_str(), IOS_OUT);
    FEI_OSTREAM& os = *id2eqnStream;
    os << "# Each line contains:\n#   ID   blk-eqn   eqn" << FEI_ENDL;
  }

  int eqnNumberOffset = 0;
  int maxNumDOF = 0;

  for(size_t rec=0; rec<recordCollections_.size(); ++rec) {
    snl_fei::RecordCollection* records = recordCollections_[rec];

    snl_fei::RecordCollection::map_type& rmap = records->getRecords();

    snl_fei::RecordCollection::map_type::iterator
      r_iter = rmap.begin(),
      r_end = rmap.end();

    int* eqnNumPtr = &eqnNumbers_[0];

    for(; r_iter != r_end; ++r_iter) {
      fei::Record<int>& thisrecord = *((*r_iter).second);

      fei::FieldMask* mask = thisrecord.getFieldMask();

      thisrecord.setOffsetIntoEqnNumbers(eqnNumberOffset);

      int owner = thisrecord.getOwnerProc();
      if (owner == proc) {
        thisrecord.setNumber(idNumber++);
      }

      int* eqnNumbers = eqnNumPtr
                      + thisrecord.getOffsetIntoEqnNumbers();

      int numDOF = mask->getNumIndices();
      eqnNumberOffset += numDOF;

      if (output_level_ >= fei::BRIEF_LOGS) {
        for(int ii=0; ii<numDOF; ++ii) {
          if (isLogEqn(eqnNumber+ii) && output_stream_ != NULL) {
            FEI_OSTREAM& os = *output_stream_;
            os <<dbgprefix_<< "setLocalEqnNumbers: ID "<<thisrecord.getID()
               << " <--> eqn " << eqnNumber+ii<<FEI_ENDL;
          }
        }
      }

      if (owner == proc) {
        int offset = 0;
        for(int n=0; n<numDOF; ++n) {
          eqnNumbers[offset++] = eqnNumber++;
        }
      }

      if (numDOF > maxNumDOF) maxNumDOF = numDOF;

      if (owner == proc) {
        int thiseqn = eqnNumber-numDOF;
        int thisrecordnumber = thisrecord.getNumber();
        if (maxNumIndices > 1) {
          CHK_ERR( ptBlkMap_->setEqn(thiseqn, thisrecordnumber, numDOF) );
          if (numDOF > 1) {
            for(int i=1; i<numDOF; ++i) {
              CHK_ERR( ptBlkMap_->setEqn(thiseqn+i, thisrecordnumber, numDOF) );
            }
          }
        }
      }

      if (id2eqnStream != NULL) {
        if (owner == proc) {
          FEI_OSTREAM& os = *id2eqnStream;
          for(int n=0; n<numDOF; ++n) {
            os << thisrecord.getID() << " " << thisrecord.getNumber() << " "
               << eqnNumber-numDOF+n<<FEI_ENDL;
          }
        }
      }
    }

  }

  ptBlkMap_->setMaxBlkEqnSize(maxNumDOF);

  int globalMaxNumDOF;
  CHK_ERR( fei::GlobalMax(comm_, maxNumDOF, globalMaxNumDOF) );

  if (globalMaxNumDOF == 1) {
    ptBlkMap_->setPtEqualBlk();
  }

  delete id2eqnStream;

  return(0);
}

//----------------------------------------------------------------------------
int fei::VectorSpace::exchangeFieldInfo(fei::comm_map* ownerPattern,
                                        fei::comm_map* sharerPattern,
                                  snl_fei::RecordCollection* recordCollection,
                                      std::vector<fei::FieldMask*>& fieldMasks)
{
  //ownerPattern maps owning processors to lists of IDs that we (the local
  //processor) share.
  //
  //sharerPattern maps sharing processors to lists of IDs that we own.
  //
  //In this function we need to perform these tasks:
  //1. processors exchange and combine their sets of field-masks so that all
  //   processors will have the super-set of field-masks that they need.
  //2. sharing processors send maskIDs for their shared IDs to the owning
  //   processors. The owning processors then combine masks as necessary to
  //   make sure each shared ID has the union of field-masks that are held by
  //   each of the sharing processors. all owning processors now send their
  //   maskIDs out to the sharing processors. At this point all processors
  //   should have the same view of the masks that belong on shared IDs.
  //3. exchange info describing inactive DOF offsets for shared records.
  //

  int numProcs = fei::numProcs(comm_);
  if (numProcs < 2) return(0);

  if (ptBlkMap_ == 0) ptBlkMap_ = new snl_fei::PointBlockMap;

  snl_fei::RecordMsgHandler recMsgHandler(fei::localProc(comm_),
                                          recordCollection,
                                          *ptBlkMap_,
                                          fieldMasks, eqnNumbers_);

  //Step 1a.
  recMsgHandler.setTask(snl_fei::RecordMsgHandler::_FieldMasks_);

  recMsgHandler.setSendPattern(ownerPattern);
  recMsgHandler.setRecvPattern(sharerPattern);
  CHK_ERR( fei::exchange(comm_, &recMsgHandler) );

  //Step 2a.
  recMsgHandler.setTask(snl_fei::RecordMsgHandler::_MaskIDs_);

  recMsgHandler.setSendPattern(ownerPattern);
  recMsgHandler.setRecvPattern(sharerPattern);
  CHK_ERR( fei::exchange(comm_, &recMsgHandler) );

  //Step 1b.
  recMsgHandler.setTask(snl_fei::RecordMsgHandler::_FieldMasks_);

  recMsgHandler.setSendPattern(sharerPattern);
  recMsgHandler.setRecvPattern(ownerPattern);
  CHK_ERR( fei::exchange(comm_, &recMsgHandler) );

  //Step 2b.
  recMsgHandler.setTask(snl_fei::RecordMsgHandler::_MaskIDs_);

  recMsgHandler.setSendPattern(sharerPattern);
  recMsgHandler.setRecvPattern(ownerPattern);
  CHK_ERR( fei::exchange(comm_, &recMsgHandler) );

  return(0);
}

//----------------------------------------------------------------------------
fei::FieldDofMap<int>&
fei::VectorSpace::getFieldDofMap()
{
  return fieldDofMap_;
}

//----------------------------------------------------------------------------
void fei::VectorSpace::setName(const char* name)
{
  if (name == NULL) return;

  if (name_ == name) return;

  name_ = name;
  dbgprefix_ = "VecSpc_"+name_+": ";
}


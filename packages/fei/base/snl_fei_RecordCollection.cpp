/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_sstream.hpp"

#include "fei_FieldMask.hpp"
#include "fei_Record.hpp"
#include "snl_fei_RecordCollection.hpp"
#include "fei_SharedIDs.hpp"

#undef fei_file
#define fei_file "snl_fei_RecordCollection.cpp"
#include "fei_ErrMacros.hpp"

//----------------------------------------------------------------------------
snl_fei::RecordCollection::RecordCollection(int localProc)
  : m_records(),
    m_global_to_local(),
    localProc_(localProc),
    debugOutput_(false),
    dbgOut_(NULL)
{
  m_records.reserve(2000);
}

//----------------------------------------------------------------------------
snl_fei::RecordCollection::RecordCollection(const RecordCollection& src)
  : m_records(src.m_records),
    m_global_to_local(src.m_global_to_local),
    localProc_(src.localProc_),
    debugOutput_(src.debugOutput_),
    dbgOut_(src.dbgOut_)
{
}

//----------------------------------------------------------------------------
snl_fei::RecordCollection::~RecordCollection()
{
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::initRecords(int numIDs, const int* IDs,
                                            std::vector<fei::FieldMask*>& fieldMasks,
                                            int* recordLocalIDs)
{
  int maskID = 0;
  fei::FieldMask* mask = NULL;
  for(unsigned m=0; m<fieldMasks.size(); ++m) {
    if (maskID == fieldMasks[m]->getMaskID()) {
      mask = fieldMasks[m]; break;
    }
  }

  if (mask == NULL) {
    mask = new fei::FieldMask();
    maskID = mask->getMaskID();
    fieldMasks.push_back(mask);
  }

  for(int i=0; i<numIDs; ++i) {
    int local_id;
    std::map<int,int>::iterator iter = m_global_to_local.lower_bound(IDs[i]);
    if (iter == m_global_to_local.end() || iter->first != IDs[i]) {
      //record doesn't exist, so we'll add a new one.
      local_id = m_records.size();
      m_global_to_local.insert(iter, std::make_pair(IDs[i], local_id));
      fei::Record<int> record;
      record.setID(IDs[i]);
      record.setFieldMask(mask);
      record.setOwnerProc(localProc_);
      m_records.push_back(record);
    }
    else {
      local_id = iter->second;
    }

    if (recordLocalIDs != NULL) recordLocalIDs[i] = local_id;
  }
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::initRecords(int fieldID, int fieldSize,
                                              int numIDs, const int* IDs,
                                              std::vector<fei::FieldMask*>& fieldMasks,
                                              int* recordLocalIDs)
{
  int maskID = fei::FieldMask::calculateMaskID(1, &fieldID);
  fei::FieldMask* mask = NULL;
  for(unsigned m=0; m<fieldMasks.size(); ++m) {
    if (maskID == fieldMasks[m]->getMaskID()) {
      mask = fieldMasks[m]; break;
    }
  }

  if (mask == NULL) {
    mask = new fei::FieldMask(1, &fieldID, &fieldSize);
    maskID = mask->getMaskID();
    fieldMasks.push_back(mask);
  }

  int lastMaskID = maskID;
  fei::FieldMask* lastMask = mask;

  for(int i=0; i<numIDs; ++i) {
    int local_id;
    std::map<int,int>::iterator iter = m_global_to_local.lower_bound(IDs[i]);
    if (iter == m_global_to_local.end() || iter->first != IDs[i]) {
      //record doesn't exist, so we'll add a new one.
      local_id = m_records.size();
      m_global_to_local.insert(iter, std::make_pair(IDs[i], local_id));
      fei::Record<int> record;
      record.setID(IDs[i]);
      record.setFieldMask(mask);
      record.setOwnerProc(localProc_);
      m_records.push_back(record);

      if (recordLocalIDs != NULL) {
        recordLocalIDs[i] = local_id;
      }
    }
    else {
      local_id = iter->second;
      fei::Record<int>& record = m_records[local_id];

      if (recordLocalIDs != NULL) {
        recordLocalIDs[i] = local_id;
      }

      int thisMaskID = record.getFieldMask()->getMaskID();

      fei::FieldMask* thisMask = record.getFieldMask();
      if (maskID == thisMaskID || thisMask->hasFieldID(fieldID)) {
        continue;
      }

      if (lastMaskID == thisMaskID) {
        record.setFieldMask(lastMask);
       continue;
      }

      int newMaskID = fei::FieldMask::calculateMaskID(*thisMask, fieldID);
      if (lastMaskID == newMaskID) {
        record.setFieldMask(lastMask);
        continue;
      }

      bool newMaskAlreadyExists = false;
      for(unsigned m=0; m<fieldMasks.size(); ++m) {
        if (newMaskID == fieldMasks[m]->getMaskID()) {
          lastMask = fieldMasks[m];
          lastMaskID = lastMask->getMaskID();
          record.setFieldMask(lastMask);
          newMaskAlreadyExists = true;
          break;
        }
      }

      if (!newMaskAlreadyExists) {
        fei::FieldMask* newmask = new fei::FieldMask(*record.getFieldMask());
        newmask->addField(fieldID, fieldSize);
        record.setFieldMask(newmask);
        fieldMasks.push_back(newmask);
        lastMask = newmask;
        lastMaskID = lastMask->getMaskID();
      }
    }
  }
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::
setOwners_lowestSharing(fei::SharedIDs<int>& sharedIDs)
{
  fei::SharedIDs<int>::map_type::iterator
    s_beg = sharedIDs.getSharedIDs().begin(),
    s_end = sharedIDs.getSharedIDs().end(),
    s_it;

  std::vector<int>& owningProcs = sharedIDs.getOwningProcs();

  int i=0;
  for(i=0, s_it = s_beg; s_it != s_end; ++i, ++s_it) {
    int sh_id = s_it->first;
    fei::Record<int>* record = getRecordWithID(sh_id);
    if (record == NULL) continue;

    int proc = owningProcs[i];

    if (debugOutput_) {
      *dbgOut_ << "#   setting ID " << (int)(record->getID())
               << "'s owner to proc " << proc << FEI_ENDL;
    }

    record->setOwnerProc(proc);
  }
}

fei::Record<int>* snl_fei::RecordCollection::getRecordWithID(int ID)
{
  std::map<int,int>::iterator iter = m_global_to_local.find(ID);

  if (iter == m_global_to_local.end()) {
    return( NULL );
  }

  return(&m_records[iter->second]);
}

const fei::Record<int>* snl_fei::RecordCollection::getRecordWithID(int ID) const
{
  std::map<int,int>::const_iterator iter = m_global_to_local.find(ID);

  if (iter == m_global_to_local.end()) {
    return( NULL );
  }

  return(&m_records[iter->second]);
}

int snl_fei::RecordCollection::getGlobalBlkIndex(int ID, int& globalBlkIndex)
{
  fei::Record<int>* record = getRecordWithID(ID);
  if (record == NULL) {
    globalBlkIndex = -1;
    ERReturn(-1);
  }

  globalBlkIndex = record->getNumber();
  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::RecordCollection::getGlobalIndex(int ID,
                                              int fieldID,
                                              int fieldSize,
                                              int fieldOffset,
                                              int whichComponentOfField,
                                              const int* eqnNumbers)
{
  fei::Record<int>* record = getRecordWithID(ID);
  if (record == NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "snl_fei::RecordCollection::getGlobalIndex ERROR, no record with "
       << "ID=" << ID;
    throw std::runtime_error(osstr.str());
  }

  fei::FieldMask* mask = record->getFieldMask();
  int offset = 0;
  try {
    mask->getFieldEqnOffset(fieldID, offset);
  }
  catch (...) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "failed to get eqn-offset for fieldID " << fieldID
          << " on record with ID " << ID << ".";
    throw std::runtime_error(osstr.str());
  }

  const int* eqnNums = eqnNumbers + record->getOffsetIntoEqnNumbers();
  if (eqnNums == NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "snl_fei::RecordCollection::getGlobalIndex ERROR: null pointer,"
         << " possibly because initComplete() hasn't been called yet?";
    throw std::runtime_error(osstr.str());
  }

  int globalIndex = -1;
  if (fieldOffset > 0) {
    globalIndex = eqnNums[offset + fieldOffset*fieldSize + whichComponentOfField];
  }
  else {
    globalIndex = eqnNums[offset + whichComponentOfField];
  }

  return(globalIndex);
}

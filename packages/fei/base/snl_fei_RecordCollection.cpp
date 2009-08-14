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
  : records_(),
    localProc_(localProc),
    recordPool_(),
    debugOutput_(false),
    dbgOut_(NULL)
{
}

//----------------------------------------------------------------------------
snl_fei::RecordCollection::RecordCollection(const RecordCollection& src)
  : records_(),
    localProc_(src.localProc_),
    recordPool_(),
    debugOutput_(src.debugOutput_),
    dbgOut_(src.dbgOut_)
{
  map_type& srcRecords =
    const_cast<map_type&>(src.records_);
  map_type::iterator
    iter = srcRecords.begin(),
    iter_end = srcRecords.end();

  map_type::iterator records_end = records_.end();

  static fei::Record dummyRecord;

  for(; iter != iter_end; ++iter) {
    map_type::value_type srcpair = *iter;
    int srcID = srcpair.first;
    fei::Record* srcRec = srcpair.second;

    fei::Record* record = recordPool_.allocate(1);
    recordPool_.construct(record,dummyRecord);
    record->deepCopy(*srcRec);
    records_.insert(records_end, map_type::value_type(srcID, record));
  }
}

//----------------------------------------------------------------------------
snl_fei::RecordCollection::~RecordCollection()
{
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::initRecords(int numIDs, const int* IDs,
                                            std::vector<fei::FieldMask*>& fieldMasks,
                                            fei::Record** records)
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

  static fei::Record dummyRecord;

  for(int i=0; i<numIDs; ++i) {
    fei::Record* record = NULL;

    map_type::iterator riter = records_.lower_bound(IDs[i]);
    if (riter != records_.end()) {
      if ((*riter).first != IDs[i]) {
        record = recordPool_.allocate(1);
        recordPool_.construct(record,dummyRecord);

        record->setID(IDs[i]);
        record->setFieldMask(mask);
        record->setOwnerProc(localProc_);

        records_.insert(riter, map_type::value_type(IDs[i], record));
      }
      else {
        record = (*riter).second;

        record->setFieldMask(mask);
      }
    }
    else {
      record = recordPool_.allocate(1);
      recordPool_.construct(record,dummyRecord);

      record->setID(IDs[i]);
      record->setFieldMask(mask);
      record->setOwnerProc(localProc_);

      records_.insert(riter, map_type::value_type(IDs[i], record));
    }

    if (records != NULL) {
      records[i] = record;
    }
  }
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::initRecords(int fieldID, int fieldSize,
                                            int numInstances,
                                            int numIDs, const int* IDs,
                                            std::vector<fei::FieldMask*>& fieldMasks,
                                            bool skipIDsWithThisField)
{
  int maskID = fei::FieldMask::calculateMaskID(1, &fieldID,
                                                      &numInstances);
  fei::FieldMask* mask = NULL;
  for(unsigned m=0; m<fieldMasks.size(); ++m) {
    if (maskID == fieldMasks[m]->getMaskID()) {
      mask = fieldMasks[m]; break;
    }
  }

  if (mask == NULL) {
    mask = new fei::FieldMask(1, &fieldID, &fieldSize, &numInstances);
    maskID = mask->getMaskID();
    fieldMasks.push_back(mask);
  }

  fei::FieldMask* lastMask = mask;
  int lastMaskID = maskID;
  static fei::Record dummyRecord;

  for(int i=0; i<numIDs; ++i) {
    fei::Record* record = NULL;

    map_type::iterator riter = records_.lower_bound(IDs[i]);

    if (riter != records_.end()) {
      if ((*riter).first != IDs[i]) {
        record = recordPool_.allocate(1);
        recordPool_.construct(record,dummyRecord);

        record->setID(IDs[i]);
        record->setFieldMask(mask);
        record->setOwnerProc(localProc_);

        records_.insert(riter, map_type::value_type(IDs[i], record));
      }
      else {
        record = (*riter).second;

        int thisMaskID = record->getFieldMask()->getMaskID();

        if (skipIDsWithThisField) {
          if (maskID == thisMaskID) continue;
          if (record->getFieldMask()->hasFieldID(fieldID)) continue;
        }

        if (lastMaskID == thisMaskID) {
          record->setFieldMask(lastMask);
          continue;
        }

        fei::FieldMask* thisMask = record->getFieldMask();
        int newMaskID = fei::FieldMask::calculateMaskID(*thisMask,
                                                            fieldID, numInstances);
        if (lastMaskID == newMaskID) {
          record->setFieldMask(lastMask);
          continue;
        }

        bool newMaskAlreadyExists = false;
        for(unsigned m=0; m<fieldMasks.size(); ++m) {
          if (newMaskID == fieldMasks[m]->getMaskID()) {
            lastMask = fieldMasks[m];
            lastMaskID = lastMask->getMaskID();
            record->setFieldMask(lastMask);
            newMaskAlreadyExists = true;
            break;
          }
        }

        if (!newMaskAlreadyExists) {
          fei::FieldMask* newmask = new fei::FieldMask(*record->getFieldMask());
          newmask->addField(fieldID, fieldSize, numInstances);
          record->setFieldMask(newmask);
          fieldMasks.push_back(newmask);
          lastMask = newmask;
          lastMaskID = lastMask->getMaskID();
        }
      }
    }
    else {
      record = recordPool_.allocate(1);
      recordPool_.construct(record,dummyRecord);
      record->setID(IDs[i]);
      record->setFieldMask(mask);
      record->setOwnerProc(localProc_);

      records_.insert(riter, map_type::value_type(IDs[i], record));
    }
  }
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::initRecords(int fieldID, int fieldSize,
                                              int numInstances,
                                              int numIDs, const int* IDs,
                                              std::vector<fei::FieldMask*>& fieldMasks,
                                              fei::Record** records,
                                              bool skipIDsWithThisField)
{
  int maskID = fei::FieldMask::calculateMaskID(1, &fieldID,
                                                      &numInstances);
  fei::FieldMask* mask = NULL;
  for(unsigned m=0; m<fieldMasks.size(); ++m) {
    if (maskID == fieldMasks[m]->getMaskID()) {
      mask = fieldMasks[m]; break;
    }
  }

  if (mask == NULL) {
    mask = new fei::FieldMask(1, &fieldID, &fieldSize, &numInstances);
    maskID = mask->getMaskID();
    fieldMasks.push_back(mask);
  }

  fei::FieldMask* lastMask = mask;
  int lastMaskID = maskID;
  static fei::Record dummyRecord;

  map_type::iterator rend = records_.end();

  for(int i=0; i<numIDs; ++i) {
    fei::Record* record = NULL;
    int ID = IDs[i];

    map_type::iterator riter = records_.lower_bound(ID);

    if (riter != rend) {
      const map_type::value_type& rval = *riter;

      if (rval.first != ID) {
        record = recordPool_.allocate(1);
        recordPool_.construct(record,dummyRecord);

        record->setID(ID);
        record->setFieldMask(mask);
        record->setOwnerProc(localProc_);

        records_.insert(riter, map_type::value_type(ID, record));

        records[i] = record;
      }
      else {
        record = rval.second;

        records[i] = record;

        int thisMaskID = record->getFieldMask()->getMaskID();

        if (skipIDsWithThisField) {
          if (maskID == thisMaskID) continue;
          if (record->getFieldMask()->hasFieldID(fieldID)) continue;
        }

        if (lastMaskID == thisMaskID) {
          record->setFieldMask(lastMask);
          continue;
        }

        fei::FieldMask* thisMask = record->getFieldMask();
        int newMaskID = fei::FieldMask::calculateMaskID(*thisMask,
                                                            fieldID, numInstances);
        if (lastMaskID == newMaskID) {
          record->setFieldMask(lastMask);
          continue;
        }

        bool newMaskAlreadyExists = false;
        for(unsigned m=0; m<fieldMasks.size(); ++m) {
          if (newMaskID == fieldMasks[m]->getMaskID()) {
            lastMask = fieldMasks[m];
            lastMaskID = lastMask->getMaskID();
            record->setFieldMask(lastMask);
            newMaskAlreadyExists = true;
            break;
          }
        }

        if (!newMaskAlreadyExists) {
          fei::FieldMask* newmask = new fei::FieldMask(*record->getFieldMask());
          newmask->addField(fieldID, fieldSize, numInstances);
          record->setFieldMask(newmask);
          fieldMasks.push_back(newmask);
          lastMask = newmask;
          lastMaskID = lastMask->getMaskID();
        }
      }
    }
    else {
      record = recordPool_.allocate(1);
      recordPool_.construct(record,dummyRecord);
      record->setID(ID);
      record->setFieldMask(mask);
      record->setOwnerProc(localProc_);

      records_.insert(riter, map_type::value_type(ID, record));

      records[i] = record;
    }
  }
}

//----------------------------------------------------------------------------
void snl_fei::RecordCollection::
setOwners_lowestSharing(fei::SharedIDs* sharedIDs)
{
  fei::SharedIDs::table_type::iterator
    s_beg = sharedIDs->getSharedIDs().getMap().begin(),
    s_end = sharedIDs->getSharedIDs().getMap().end(),
    s_it;

  std::vector<int>& owningProcs = sharedIDs->getOwningProcs();

  map_type::iterator rend = records_.end();

  int i;
  for(i=0, s_it = s_beg; s_it != s_end; ++i, ++s_it) {
    fei::Record* record = NULL;
    int sh_id = (*s_it).first;
    map_type::iterator riter = records_.find(sh_id);
    if (riter == rend) continue;

    record = (*riter).second;

    int proc = owningProcs[i];

    if (debugOutput_) {
      *dbgOut_ << "#   setting ID " << (int)(record->getID())
               << "'s owner to proc " << proc << FEI_ENDL;
    }

    record->setOwnerProc(proc);
  }
}

fei::Record* snl_fei::RecordCollection::getRecordWithID(int ID)
{
  map_type::iterator rend = records_.end();
  map_type::iterator riter = records_.find(ID);

  if (riter == rend) {
    return( NULL );
  }

  return((*riter).second);
}

int snl_fei::RecordCollection::getGlobalBlkIndex(int ID, int& globalBlkIndex)
{
  fei::Record* record = getRecordWithID(ID);
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
  fei::Record* record = getRecordWithID(ID);
  if (record == NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "snl_fei::RecordCollection::getGlobalIndex ERROR, no record with "
       << "ID=" << ID;
    throw std::runtime_error(osstr.str());
  }

  fei::FieldMask* mask = record->getFieldMask();
  int numInstances = 0;
  int offset = 0;
  try {
    mask->getFieldEqnOffset(fieldID, offset, numInstances);
  }
  catch (std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "failed to get eqn-offset for fieldID " << fieldID
          << " on record with ID " << ID << ".";
    throw std::runtime_error(osstr.str());
  }

  if (fieldOffset >= numInstances) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "snl_fei::RecordCollection::getGlobalIndex: fieldOffset ("<<fieldOffset
          << ") should be less than numInstances ("<<numInstances<<").";
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

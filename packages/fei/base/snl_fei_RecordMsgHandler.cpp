/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_sstream.hpp"

#include "snl_fei_RecordMsgHandler.hpp"
#include "snl_fei_RecordCollection.hpp"
#include "snl_fei_PointBlockMap.hpp"
#include "fei_FieldMask.hpp"
#include "fei_TemplateUtils.hpp"

#undef fei_file
#define fei_file "snl_fei_RecordMsgHandler.cpp"
#include "fei_ErrMacros.hpp"

snl_fei::RecordMsgHandler::RecordMsgHandler(int localProc,
					    RecordCollection* recordCollection,
					    snl_fei::PointBlockMap& ptBlkMap,
					    std::vector<fei::FieldMask*>& fieldMasks,
					    std::vector<int>& eqnNumbers)
  :sendPattern_(NULL),
   recvPattern_(NULL),
   recordCollection_(recordCollection),
   ptBlkMap_(ptBlkMap),
   fieldMasks_(fieldMasks),
   whichTask_(_FieldMasks_),
   sendProcs_(0, 32),
   recvProcs_(0, 32),
   eqnNumbers_(eqnNumbers),
   localProc_(localProc)
{
}

snl_fei::RecordMsgHandler::~RecordMsgHandler()
{
}

std::vector<int>& snl_fei::RecordMsgHandler::getSendProcs()
{
  fei::copyKeysToVector(sendPattern_->getMap(), sendProcs_);
  return(sendProcs_);
}

std::vector<int>& snl_fei::RecordMsgHandler::getRecvProcs()
{
  fei::copyKeysToVector(recvPattern_->getMap(), recvProcs_);
  return(recvProcs_);
}

int snl_fei::RecordMsgHandler::getSendMessageLength(int destProc,
						    int& messageLength)
{
  if (sendPattern_ == NULL || recvPattern_ == NULL) ERReturn(-1);

  switch(whichTask_) {
  case _FieldMasks_:
    messageLength = localFieldMaskMessageSize(fieldMasks_);    break;
  case _MaskIDs_:
    messageLength = sendPattern_->getRow(destProc)->size();  break;
  case _EqnNumbers_:
    messageLength = eqnNumbersMsgLength(destProc);  break;
  default:
    std::abort();
  }

  return(0);
}

int snl_fei::RecordMsgHandler::getSendMessage(int destProc,
					      std::vector<int>& message)
{
  switch(whichTask_) {
  case _FieldMasks_:
    CHK_ERR( packLocalFieldMasks(fieldMasks_, message) );    break;
  case _MaskIDs_:
    CHK_ERR( packMaskIDs(destProc, message) );  break;
  case _EqnNumbers_:
    CHK_ERR( packEqnNumbersMsg(destProc, message) );  break;
  default:
    std::abort();
  }

  return(0);
}

int snl_fei::RecordMsgHandler::processRecvMessage(int srcProc,
						  std::vector<int>& message)
{
  switch(whichTask_) {
  case _FieldMasks_:
    CHK_ERR( addFieldMasks(message, fieldMasks_) );    break;
  case _MaskIDs_:
    CHK_ERR( mergeMaskIDs(srcProc, message) );  break;
  case _EqnNumbers_:
    CHK_ERR( storeEqnNumbers(srcProc, message) );  break;
  default:
    std::abort();
  }

  return(0);
}

int snl_fei::RecordMsgHandler::
localFieldMaskMessageSize(std::vector<fei::FieldMask*>& fieldMasks)
{
  //msg-size will be 1+sum(numFields*3+2)
  int numLocalMasks = fieldMasks.size();
  int i, msgSize = 1;
  for(i=0; i<numLocalMasks; ++i) {
    msgSize += (fieldMasks[i]->getNumFields())*3+2;
  }

  return( msgSize );
}

int snl_fei::RecordMsgHandler::
packLocalFieldMasks(std::vector<fei::FieldMask*>& fieldMasks,
		    std::vector<int>& localFieldMasks)
{
  int numLocalMasks = fieldMasks.size();
  int msgSize = 1;
  for(int ii=0; ii<numLocalMasks; ++ii) {
    msgSize += fieldMasks[ii]->getNumFields()*3+2;
  }

  localFieldMasks.resize(msgSize);

  int offset = 0;
  localFieldMasks[offset++] = fieldMasks.size();
  for(unsigned i=0; i<fieldMasks.size(); ++i) {
    int numFields = fieldMasks[i]->getNumFields();
    int numIndices = fieldMasks[i]->getNumIndices();
    localFieldMasks[offset++] = numFields;
    localFieldMasks[offset++] = numIndices;

    std::vector<int>& fieldIDs = fieldMasks[i]->getFieldIDs();
    std::vector<int>& fieldSizes = fieldMasks[i]->getFieldSizes();

    for(int j=0; j<numFields; ++j) {
      localFieldMasks[offset+j] = fieldIDs[j];
      localFieldMasks[offset+numFields+j] = fieldSizes[j];
    }
    offset += numFields*3;
  }

  return(0);
}

int snl_fei::RecordMsgHandler::
addFieldMasks(std::vector<int>& msg, std::vector<fei::FieldMask*>& fieldMasks)
{
  int offset = 0;
  int* msgPtr = &msg[0];
  int numMasks = msgPtr[offset++];
  for(int i=0; i<numMasks; ++i) {
    int numFields = msgPtr[offset++];
    int numIndices = msgPtr[offset++];
    int* fieldIDs = NULL;
    int* fieldSizes = NULL;
    if (numFields > 0) {
      fieldIDs = &(msgPtr[offset]);
      fieldSizes = &(msgPtr[offset+numFields]);
    }

    int maskID = fei::FieldMask::calculateMaskID(numFields, fieldIDs);

    bool maskAlreadyExists = false;
    for(unsigned j=0; j<fieldMasks.size(); ++j) {
      int existingMaskID = fieldMasks[j]->getMaskID();
      if (maskID == existingMaskID) {
	maskAlreadyExists = true; break;
      }
    }

    if (!maskAlreadyExists) {
      fei::FieldMask* newmask = new fei::FieldMask(numFields,
							   fieldIDs, fieldSizes);
      if (numFields < 1) {
	newmask->setNumIndices(numIndices);
      }
      fieldMasks.push_back(newmask);
    }

    offset += 3*numFields;
  }

  return(0);
}

int snl_fei::RecordMsgHandler::packMaskIDs(int destProc, std::vector<int>& msg)
{
  fei::comm_map::row_type* ids = sendPattern_->getRow(destProc);
  int len = ids->size();

  msg.resize(len);

  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();

  int offset = 0;
  int* msgPtr = &msg[0];

  for(; id_iter != id_end; ++id_iter) {
    fei::Record<int>* rec = recordCollection_->getRecordWithID(*id_iter);
    if (rec == NULL) {
      FEI_OSTRINGSTREAM osstr;
      osstr << "RecordMsgHandler::packMaskIDs: proc " << localProc_
	   << " failed to find ID " << *id_iter;
      throw std::runtime_error(osstr.str());
    }

    msgPtr[offset++] = rec->getFieldMask()->getMaskID();
  }

  return(0);
}

int snl_fei::RecordMsgHandler::mergeMaskIDs(int srcProc, std::vector<int>& msg)
{
  fei::comm_map::row_type* ids = recvPattern_->getRow(srcProc);

  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();

  int offset = 0;
  int* msgPtr = &msg[0];

  for(; id_iter != id_end; ++id_iter) {
    int ID = *id_iter;
    fei::Record<int>* rec = recordCollection_->getRecordWithID(ID);
    if (rec == NULL) {
      ERReturn(-1);
    }

    int maskID = msgPtr[offset++];

    if (maskID != rec->getFieldMask()->getMaskID()) {
      //if the masks don't match, we need to add the fields from the
      //incoming-field-mask-id to our record's field-mask.

      //first, find the field-mask for 'maskID'
      fei::FieldMask* mask = NULL;
      for(unsigned fm=0; fm<fieldMasks_.size(); ++fm) {
	if (fieldMasks_[fm]->getMaskID() == maskID) mask = fieldMasks_[fm];
      }
      if (mask == NULL) {
	fei::console_out() << "mergeMaskIDs didn't find mask for maskID="
           << maskID<<FEI_ENDL;
	ERReturn(-1);
      }

      int numFields = mask->getNumFields();
      std::vector<int>& fieldIDs = mask->getFieldIDs();
      std::vector<int>& fieldSizes = mask->getFieldSizes();

      for(int nf=0; nf<numFields; ++nf) {
	recordCollection_->initRecords(fieldIDs[nf], fieldSizes[nf],
				       1, &ID, fieldMasks_);
      }
    }
  }

  return(0);
}

int snl_fei::RecordMsgHandler::eqnNumbersMsgLength(int destProc)
{
  fei::comm_map::row_type* ids = sendPattern_->getRow(destProc);
  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();
  int len = ids->size();

  len *= 3;

  for(; id_iter != id_end; ++id_iter) {
    int ID = *id_iter;
    fei::Record<int>* rec = recordCollection_->getRecordWithID(ID);
    if (rec == NULL) {
      ERReturn(-1);
    }

    len += rec->getFieldMask()->getNumIndices();
  }

  return(len);
}

int snl_fei::RecordMsgHandler::packEqnNumbersMsg(int destProc,
						 std::vector<int>& msg)
{
  fei::comm_map::row_type* ids = sendPattern_->getRow(destProc);
  int len = ids->size()*3;
  msg.resize(len);
  const int* eqnNumPtr = &eqnNumbers_[0];

  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();
  int offset = 0;
  for(; id_iter != id_end; ++id_iter) {
    int ID = *id_iter;
    fei::Record<int>* rec = recordCollection_->getRecordWithID(ID);
    if (rec == NULL) {
      ERReturn(-1);
    }

    len = rec->getFieldMask()->getNumIndices();

    msg[offset++] = ID;
    msg[offset++] = rec->getNumber();
    msg[offset++] = len;

    const int* eqnNumbers = eqnNumPtr+rec->getOffsetIntoEqnNumbers();
    for(int i=0; i<len; ++i) {
      msg.push_back(eqnNumbers[i]);
    }
  }

  return(0);
}

int snl_fei::RecordMsgHandler::storeEqnNumbers(int srcProc, std::vector<int>& msg)
{
  int numIDs = recvPattern_->getRow(srcProc)->size();
  int offset = numIDs*3;
  int* msgPtr = &msg[0];
  int* eqnNumPtr = &eqnNumbers_[0];
  for(int i=0; i<numIDs; ++i) {
    int ID = msgPtr[i*3];
    int recNumber = msgPtr[i*3+1];
    int numEqns = msgPtr[i*3+2];
    fei::Record<int>* rec = recordCollection_->getRecordWithID(ID);
    if (rec == NULL) {
      ERReturn(-1);
    }

    rec->setNumber(recNumber);
    int* eqnNumbers = eqnNumPtr+rec->getOffsetIntoEqnNumbers();
    for(int eq=0; eq<numEqns; ++eq) {
      eqnNumbers[eq] = msgPtr[offset++];
      ptBlkMap_.setEqn(eqnNumbers[eq], recNumber, numEqns);
    }
  }

  return(0);
}


/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <snl_fei_SubdMsgHandler.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <fei_SharedIDs.hpp>
#include <fei_TemplateUtils.hpp>

#undef fei_file
#define fei_file "snl_fei::SubdMsgHandler"
#include <fei_ErrMacros.hpp>

snl_fei::SubdMsgHandler::SubdMsgHandler(RecordCollection* recordCollection,
					fei::SharedIDs<int>* sharedIDTable,
					fei::SharedIDs<int>* subdomainIDTable)
  :sendPattern_(NULL),
   recvPattern_(NULL),
   recordCollection_(recordCollection),
   sharedIDTable_(sharedIDTable),
   subdomainIDTable_(subdomainIDTable),
   sendProcs_(0, 32),
   recvProcs_(0, 32)
{
}

snl_fei::SubdMsgHandler::~SubdMsgHandler()
{
}

std::vector<int>& snl_fei::SubdMsgHandler::getSendProcs()
{
  fei::copyKeysToVector(sendPattern_->getMap(), sendProcs_);
  return(sendProcs_);
}

std::vector<int>& snl_fei::SubdMsgHandler::getRecvProcs()
{
  fei::copyKeysToVector(recvPattern_->getMap(), recvProcs_);
  return(recvProcs_);
}

int snl_fei::SubdMsgHandler::getSendMessageLength(int destProc,
						  int& messageLength)
{
  fei::comm_map::row_type* list = sendPattern_->getRow(destProc);
  if (list == NULL) {
    FEI_COUT << "SubdMsdHandler: destProc: " << destProc << ", list is NULL."<<FEI_ENDL;
    return(-1);
  }

  messageLength = list->size();

  return( 0 );
}

int snl_fei::SubdMsgHandler::getSendMessage(int destProc,
					  std::vector<int>& message)
{
  int len = sendPattern_->getRow(destProc)->size();
  message.resize(len);
  int* msgPtr = &message[0];

  fei::comm_map::row_type* ids = sendPattern_->getRow(destProc);
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

    if (rec->isInLocalSubdomain_) {
      msgPtr[offset++] = 1;
    }
    else {
      msgPtr[offset++] = 0;
    }
  }

  return(0);
}

int snl_fei::SubdMsgHandler::processRecvMessage(int srcProc,
						std::vector<int>& message)
{
  fei::comm_map::row_type* ids = recvPattern_->getRow(srcProc);
  fei::comm_map::row_type::const_iterator
    id_iter = ids->begin(),
    id_end = ids->end();

  int* msgPtr = &message[0];

  if (message.size() != ids->size()) {
    ERReturn(-1);
  }

  int offset = 0;
  for(; id_iter != id_end; ++id_iter) {
    int ID = *id_iter;

    bool isInRemoteSubdomain = msgPtr[offset++] > 1 ? true : false;

    if (isInRemoteSubdomain) {
      subdomainIDTable_->addSharedID(ID, 1, &srcProc);
    }
  }

  return(0);
}


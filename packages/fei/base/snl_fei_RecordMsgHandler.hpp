/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_RecordMsgHandler_hpp_
#define _snl_fei_RecordMsgHandler_hpp_

#include <fei_macros.hpp>
#include <fei_CommUtils.hpp>
#include <fei_fwd.hpp>

namespace fei {
  class FieldMask;
}

namespace snl_fei {
  /** Implementation of MessageHandler specialized for Record objects. */
  class RecordMsgHandler : public fei::MessageHandler<int> {
  public:
    /** constructor */
    RecordMsgHandler(int localProc,
		     RecordCollection* recordCollection,
		     snl_fei::PointBlockMap& ptBlkMap,
		     std::vector<fei::FieldMask*>& fieldMasks,
		     std::vector<int>& eqnNumbers);

    /** destructor */
    virtual ~RecordMsgHandler();

    /** enumeration for operation-types */
    enum {_FieldMasks_ = 0, _MaskIDs_ = 1,
          _EqnNumbers_};

    /** Get list of processors to be sent to. */
    std::vector<int>& getSendProcs();

    /** Get list of processors to be received from. */
    std::vector<int>& getRecvProcs();

    /** Get length of message to be sent to a specified processor. */
    int getSendMessageLength(int destProc, int& messageLength);

    /** Get message data to be sent to specified processor. */
    int getSendMessage(int destProc, std::vector<int>& message);

    /** Process a message received from a specified processor. */
    int processRecvMessage(int srcProc, std::vector<int>& message);

    /** clumsy method for specifying the next operation to be performed. */
    void setTask(int task) { whichTask_ = task; }

    /** Set the pattern that specifies processors to be sent to. */
    void setSendPattern(fei::comm_map* pattern)
      { sendPattern_ = pattern; }

    /** Set the pattern that specifies processors to be received from. */
    void setRecvPattern(fei::comm_map* pattern)
      { recvPattern_ = pattern; }

  private:
    int localFieldMaskMessageSize(std::vector<fei::FieldMask*>& fieldMasks);

    int packLocalFieldMasks(std::vector<fei::FieldMask*>& fieldMasks,
                            std::vector<int>& localFieldMasks);

    int addFieldMasks(std::vector<int>& msg, std::vector<fei::FieldMask*>& fieldMasks);

    int packMaskIDs(int destProc, std::vector<int>& msg);

    int mergeMaskIDs(int srcProc, std::vector<int>& msg);

    int eqnNumbersMsgLength(int destProc);

    int packEqnNumbersMsg(int destProc, std::vector<int>& msg);

    int storeEqnNumbers(int srcProc, std::vector<int>& msg);

    fei::comm_map* sendPattern_;
    fei::comm_map* recvPattern_;
    RecordCollection* recordCollection_;
    snl_fei::PointBlockMap& ptBlkMap_;
    std::vector<fei::FieldMask*>& fieldMasks_;

    int whichTask_;

    std::vector<int> sendProcs_;
    std::vector<int> recvProcs_;

    std::vector<int>& eqnNumbers_;

    int localProc_;
  };//class RecordMsgHandler
}//namespace snl_fei

#endif // _snl_fei_RecordMsgHandler_hpp_


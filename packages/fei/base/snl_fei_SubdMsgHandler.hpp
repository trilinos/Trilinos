/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_SubdMsgHandler_hpp_
#define _snl_fei_SubdMsgHandler_hpp_

#include <fei_macros.hpp>
#include <fei_CommUtils.hpp>
#include <fei_fwd.hpp>

#include <vector>

namespace fei {
  template<typename T> class SharedIDs;
}

namespace snl_fei {
  /** implementation of MessageHandler for subdomain data */
  class SubdMsgHandler : public fei::MessageHandler<int> {
  public:
    /** constructor */
    SubdMsgHandler(RecordCollection* recordCollection,
		   fei::SharedIDs<int>* sharedIDTable,
		   fei::SharedIDs<int>* subdomainIDTable);
    /** destructor */
    virtual ~SubdMsgHandler();

    /** get list of processors to be sent to */
    std::vector<int>& getSendProcs();

    /** get list of processors to be recvd from */
    std::vector<int>& getRecvProcs();

    /** get length of message to be sent to specified proc */
    int getSendMessageLength(int destProc, int& messageLength);

    /** get message to be sent to specified proc */
    int getSendMessage(int destProc, std::vector<int>& message);

    /** process message received from specified recv proc */
    int processRecvMessage(int srcProc, std::vector<int>& message);

    /** set pattern describing procs to be sent to */
    void setSendPattern(fei::comm_map* pattern)
      { sendPattern_ = pattern; }

    /** set pattern describing procs to be recvd from */
    void setRecvPattern(fei::comm_map* pattern)
      { recvPattern_ = pattern; }

  private:
    fei::comm_map* sendPattern_;
    fei::comm_map* recvPattern_;
    RecordCollection* recordCollection_;
    fei::SharedIDs<int>* sharedIDTable_;
    fei::SharedIDs<int>* subdomainIDTable_;

    std::vector<int> sendProcs_;
    std::vector<int> recvProcs_;
  };//class SubdMsgHandler
}//namespace snl_fei

#endif // _snl_fei_SubdMsgHandler_hpp_


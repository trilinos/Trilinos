/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_BlkSizeMsgHandler_hpp_
#define _snl_fei_BlkSizeMsgHandler_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Graph.hpp>
#include <fei_CommUtils.hpp>

namespace snl_fei {

/** MessageHandler implementation for block-size data. */
class BlkSizeMsgHandler : public fei::MessageHandler<int> {
 public:
  /** constructor */
  BlkSizeMsgHandler(fei::VectorSpace* vspace,
		    fei::Graph* graph,
		    MPI_Comm comm);
  /** destructor */
  virtual ~BlkSizeMsgHandler();

  /** clumsy method to launch the data exchange. */
  int do_the_exchange();

  /** Get list of procs to send to. */
  std::vector<int>& getSendProcs();
  /** Get list of procs to recv from. */
  std::vector<int>& getRecvProcs();

  /** Get length of message for specified destination proc. */
  int getSendMessageLength(int destProc, int& messageLength);
  /** Get message to send to specified destination proc. */
  int getSendMessage(int destProc, std::vector<int>& message);
  /** process message received from specified source proc. */
  int processRecvMessage(int srcProc, std::vector<int>& message);

 private:
  fei::comm_map* remote_colIndices_;
  fei::comm_map* local_colIndices_;
  fei::VectorSpace* vecSpace_;
  snl_fei::PointBlockMap* ptBlkMap_;
  fei::Graph* graph_;
  MPI_Comm comm_;
  std::vector<int> sendProcs_;
  std::vector<int> recvProcs_;

  bool firstExchange_;
};

} // namespace snl_fei

#endif


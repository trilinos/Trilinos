/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <snl_fei_BlkSizeMsgHandler.hpp>

#include <fei_utils.hpp>

#include <snl_fei_Utils.hpp>
#include <fei_FieldMask.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_Graph.hpp>
#include <snl_fei_Constraint.hpp>
#include <fei_TemplateUtils.hpp>

#include <fei_EqnBuffer.hpp>
#include <fei_EqnCommMgr.hpp>
#include <SNL_FEI_Structure.hpp>

#undef fei_file
#define fei_file "snl_fei_BlkSizeMsgHandler.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
snl_fei::BlkSizeMsgHandler::BlkSizeMsgHandler(fei::VectorSpace* vspace,
                                              fei::Graph* graph,
                                              MPI_Comm comm)
  : remote_colIndices_(NULL),
    local_colIndices_(NULL),
    vecSpace_(vspace),
    ptBlkMap_(NULL),
    graph_(graph),
    comm_(comm),
    sendProcs_(0, 64),
    recvProcs_(0, 64),
    firstExchange_(true)
{
  remote_colIndices_ = new fei::comm_map(0,1);
  local_colIndices_ = new fei::comm_map(0,1);

  ptBlkMap_ = vspace->getPointBlockMap();
}

//----------------------------------------------------------------------------
snl_fei::BlkSizeMsgHandler::~BlkSizeMsgHandler()
{
  delete remote_colIndices_;
  delete local_colIndices_;
}

//----------------------------------------------------------------------------
int snl_fei::BlkSizeMsgHandler::do_the_exchange()
{
  int local_proc = fei::localProc(comm_);
  if (fei::numProcs(comm_) < 2) {
    return(0);
  }

  fei::Graph::table_type* localgraph = graph_->getLocalGraph();
  fei::Graph::table_type::iterator
    g_iter = localgraph->begin(),
    g_end = localgraph->end();

  //First create a table that maps remote processors to column-indices from our
  //graph.
  //These are remotely-owned column-indices for which we will need block-sizes.

  for(; g_iter != g_end; ++g_iter) {
    fei::Graph::table_type::row_type* row = (*g_iter).second;

    fei::Graph::table_type::row_type::const_iterator
      iter = row->begin(),
      iter_end = row->end();

    int owner;

    for(; iter != iter_end; ++iter) {
      int col = *iter;
      owner = vecSpace_->getOwnerProcBlkIndex(col);

      if (owner != local_proc) {
        remote_colIndices_->addIndices(owner, 1, &col);
      }
    }
  }

  //Next, we need to send our lists of remotely-owned column-indices to the
  //owning processors. After that, those processors can respond by sending us
  //the sizes for those column-indices.
  fei::copyKeysToVector(remote_colIndices_->getMap(), sendProcs_);

  CHK_ERR( fei::mirrorProcs(comm_, sendProcs_, recvProcs_) );

  firstExchange_ = true;

  CHK_ERR( fei::exchange(comm_, this) );

  firstExchange_ = false;

  CHK_ERR( fei::exchange(comm_, this) );

  return(0);
}

//----------------------------------------------------------------------------
std::vector<int>& snl_fei::BlkSizeMsgHandler::getSendProcs()
{
  if (firstExchange_) {
    return(sendProcs_);
  }
  else {
    return(recvProcs_);
  }
}

//----------------------------------------------------------------------------
std::vector<int>& snl_fei::BlkSizeMsgHandler::getRecvProcs()
{
  if (firstExchange_) {
    return(recvProcs_);
  }
  else {
    return(sendProcs_);
  }
}

//----------------------------------------------------------------------------
int snl_fei::BlkSizeMsgHandler::getSendMessageLength(int destProc,
                                                     int& messageLength)
{
  if (firstExchange_) {
    fei::comm_map::row_type* cols = remote_colIndices_->getRow(destProc);
    messageLength = cols->size();
    return(0);
  }
  else {
    fei::comm_map::row_type* cols = local_colIndices_->getRow(destProc);
    messageLength = cols->size()*2;
    return(0);
  }
}

//----------------------------------------------------------------------------
int snl_fei::BlkSizeMsgHandler::getSendMessage(int destProc,
                                               std::vector<int>& message)
{
  if (firstExchange_) {
    fei::comm_map::row_type* cols = remote_colIndices_->getRow(destProc);
    message.resize(cols->size());
    fei::copySetToArray(*cols, message.size(), &message[0]);
    return(0);
  }
  else {
    fei::comm_map::row_type* cols = local_colIndices_->getRow(destProc);

    message.resize(cols->size()*2);

    fei::comm_map::row_type::const_iterator
      iter = cols->begin(),
      iter_end = cols->end();

    int offset = 0;
    for(; iter != iter_end; ++iter) {
      CHK_ERR( ptBlkMap_->getBlkEqnInfo(*iter,
                                        message[offset], message[offset+1]) );
      offset += 2;
    }

    return( 0 );
  }
}

//----------------------------------------------------------------------------
int snl_fei::BlkSizeMsgHandler::processRecvMessage(int srcProc,
                                                   std::vector<int>& message)
{
  if (firstExchange_) {
    for(unsigned i=0; i<message.size(); ++i) {
      local_colIndices_->addIndices(srcProc, 1, &(message[i]));
    }
  }
  else {
    fei::comm_map::row_type* cols = remote_colIndices_->getRow(srcProc);
    fei::comm_map::row_type::const_iterator
      iter = cols->begin(),
      iter_end = cols->end();

    int offset = 0;
    for(; iter != iter_end; ++iter) {
      int ptEqn = message[offset];
      int blkSize = message[offset+1];
      for(int i=0; i<blkSize; ++i) {
        CHK_ERR( ptBlkMap_->setEqn(ptEqn+i, *iter, blkSize) );
      }
      offset += 2;
    }
  }

  return(0);
}


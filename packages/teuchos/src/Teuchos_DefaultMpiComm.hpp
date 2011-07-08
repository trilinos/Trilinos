// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_MPI_COMM_HPP
#define TEUCHOS_MPI_COMM_HPP


#include "Teuchos_Comm.hpp"
#include "Teuchos_CommUtilities.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_MpiReductionOpSetter.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Assert.hpp"
#include "mpi.h"


// This must be defined globally for the whole program!
//#define TEUCHOS_MPI_COMM_DUMP


#ifdef TEUCHOS_MPI_COMM_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif


namespace Teuchos {


#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal, typename T>
void dumpBuffer(
  const std::string &funcName, const std::string &buffName
  ,const Ordinal bytes, const T buff[]
  )
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out
    << "\n" << funcName << "::" << buffName << ":\n";
  tab.incrTab();
  for( Ordinal i = 0; i < bytes; ++i ) {
    *out << buffName << "[" << i << "] = '" << buff[i] << "'\n";
  }
  *out << "\n";
}
#endif // TEUCHOS_MPI_COMM_DUMP


/** \brief . */
class MpiCommRequest : public CommRequest {
public:
  /** \brief . */
  MpiCommRequest( MPI_Request rawMpiRequest )
    :rawMpiRequest_(rawMpiRequest)
    {}
  /** \brief . */
  MPI_Request releaseRawMpiRequest()
    {
      MPI_Request tmp_rawMpiRequest = rawMpiRequest_;
      rawMpiRequest_ = MPI_REQUEST_NULL;
      return tmp_rawMpiRequest;
    }
private:
  MPI_Request rawMpiRequest_;
  MpiCommRequest(); // Not defined
};


/** \brief Nonmember constructor .
 *
 * \relates MpiCommRequest
 */
inline
const RCP<MpiCommRequest>
mpiCommRequest( MPI_Request rawMpiRequest  )
{
  return Teuchos::rcp(new MpiCommRequest(rawMpiRequest));
}


/** \brief Concrete communicator subclass based on MPI.
 *
 * <b>Assertions:</b><ul>
 * <li><tt>getRawMpiComm().get()!=NULL && *getRawMpiComm()!=MPI_COMM_NULL</tt>
 * <li><tt>getSize() > 0</tt>
 * <li><tt>0 <= getRank() && getRank() < getSize()</tt>
 * </ul>
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal>
class MpiComm : public Comm<Ordinal> {
public:

  //! @name Constructors 
  //@{

  /** \brief Construct given a wrapped MPI_Comm oqaque object.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL</tt>
   * </ul>
   */
  MpiComm(
    const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
    );

  /** \brief Return the embedded wrapped opaque <tt>MPI_Comm</tt> object. */
  RCP<const OpaqueWrapper<MPI_Comm> > getRawMpiComm() const
  {return rawMpiComm_;}

  //@}

  //! @name Overridden from Comm 
  //@{

  /** \brief . */
  virtual int getRank() const;
  /** \brief . */
  virtual int getSize() const;
  /** \brief . */
  virtual void barrier() const;
  /** \brief . */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const;
  /** \brief . */
  virtual void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const;
  /** \brief . */
  virtual void reduceAllAndScatter(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvCounts[], char myGlobalReducts[]
    ) const;
  /** \brief . */
	virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> ireceive(
    const ArrayView<char> &Buffer,
    const int sourceRank
    ) const;
  /** \brief . */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest> > &requests
    ) const;
  /** \brief . */
  virtual void wait(
    const Ptr<RCP<CommRequest> > &request
    ) const;

  //@}

  //! @name Overridden from Describable 
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  // These should be private but the PGI compiler requires them be public

  static int const minTag_ = 26000; // These came from Teuchos::MpiComm???
  static int const maxTag_ = 26099; // ""

private:

  static int tagCounter_;

  RCP<const OpaqueWrapper<MPI_Comm> > rawMpiComm_;
  int rank_;
  int size_;
  int tag_;

  void assertRank(const int rank, const std::string &rankName) const;

  // Not defined and not to be called!
  MpiComm();

#ifdef TEUCHOS_MPI_COMM_DUMP
public:
  static bool show_dump;
#endif // TEUCHOS_MPI_COMM_DUMP
	
};


/** \brief Helper function that creates a dynamically allocated
 * <tt>MpiComm</tt> object or returns <tt>Teuchos::null</tt> to correctly
 * represent a null communicator.
 *
 * <b>Postconditions:</b></ul>
 * <li>[<tt>rawMpiComm.get()!=NULL && *rawMpiComm!=MPI_COMM_NULL</tt>]
 *     <tt>return.get()!=NULL</tt>
 * <li>[<tt>rawMpiComm.get()==NULL || *rawMpiComm==MPI_COMM_NULL</tt>]
 *     <tt>return.get()==NULL</tt>
 * </ul>
 *
 * \relates MpiComm
 */
template<typename Ordinal>
RCP<MpiComm<Ordinal> >
createMpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  );


// ////////////////////////
// Implementations


// Static members


template<typename Ordinal>
int MpiComm<Ordinal>::tagCounter_ = MpiComm<Ordinal>::minTag_;


// Constructors


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  TEST_FOR_EXCEPT( rawMpiComm.get()==NULL );
  TEST_FOR_EXCEPT( *rawMpiComm == MPI_COMM_NULL );
  rawMpiComm_ = rawMpiComm;
  MPI_Comm_size(*rawMpiComm_,&size_);
  MPI_Comm_rank(*rawMpiComm_,&rank_);
  if(tagCounter_ > maxTag_)
    tagCounter_ = minTag_;
  tag_ = tagCounter_++;
}


// Overridden from Comm

  
template<typename Ordinal>
int MpiComm<Ordinal>::getRank() const
{
  return rank_;
}

  
template<typename Ordinal>
int MpiComm<Ordinal>::getSize() const
{
  return size_;
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::barrier() const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::barrier()"
    );
  MPI_Barrier(*rawMpiComm_);
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::broadcast(
  const int rootRank, const Ordinal bytes, char buffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::broadcast(...)"
    );
  MPI_Bcast(buffer,bytes,MPI_CHAR,rootRank,*rawMpiComm_);
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::gatherAll(
  const Ordinal sendBytes, const char sendBuffer[],
  const Ordinal recvBytes, char recvBuffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::gatherAll(...)"
    );
  TEUCHOS_ASSERT_EQUALITY((sendBytes*size_), recvBytes );
  MPI_Allgather(
    const_cast<char *>(sendBuffer), sendBytes, MPI_CHAR,
    recvBuffer, sendBytes, MPI_CHAR,
    *rawMpiComm_
    );
  // NOTE: 'sendBytes' is being sent above for the MPI arg recvcount (which is
  // very confusing in the MPI documentation) for MPI_Allgether(...).
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::reduceAll(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::reduceAll(...)"
    );
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  MPI_Allreduce(
    const_cast<char*>(sendBuffer),globalReducts,bytes,MPI_CHAR,op.mpi_op()
    ,*rawMpiComm_
    );
}


template<typename Ordinal>
void MpiComm<Ordinal>::reduceAllAndScatter(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal sendBytes, const char sendBuffer[]
  ,const Ordinal recvCounts[], char myGlobalReducts[]
  ) const
{

  (void)sendBytes; // Ignore if not in debug mode

  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::reduceAllAndScatter(...)"
    );

#ifdef TEUCHOS_DEBUG
  Ordinal sumRecvBytes = 0;
  for( Ordinal i = 0; i < size_; ++i ) {
    sumRecvBytes += recvCounts[i];
  }
  TEST_FOR_EXCEPT(!(sumRecvBytes==sendBytes));
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "sendBuffer", sendBytes, sendBuffer );
    dumpBuffer<Ordinal,Ordinal>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "recvCounts", as<Ordinal>(size_), recvCounts );
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "myGlobalReducts", as<char>(recvCounts[rank_]), myGlobalReducts );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  // Create a new recvCount[] if Ordinal!=int
  WorkspaceStore* wss = get_default_workspace_store().get();
  const bool Ordinal_is_int = typeid(int)==typeid(Ordinal);
  Workspace<int> ws_int_recvCounts(wss,Ordinal_is_int?0:size_);
  const int *int_recvCounts = 0;
  if(Ordinal_is_int) {
    int_recvCounts = reinterpret_cast<const int*>(recvCounts);
    // Note: We must do an reinterpet cast since this must
    // compile even if it is not executed.  I could implement
    // code that would not need to do this using template
    // conditionals but I don't want to bother.
  }
  else {
    std::copy(recvCounts, recvCounts+size_, &ws_int_recvCounts[0]);
    int_recvCounts = &ws_int_recvCounts[0];
  }

  // Perform the operation
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp, false)));
  MPI_Reduce_scatter(
    const_cast<char*>(sendBuffer), myGlobalReducts,
    const_cast<int*>(int_recvCounts),
    MPI_CHAR,
    op.mpi_op(),
    *rawMpiComm_
    );

}


template<typename Ordinal>
void MpiComm<Ordinal>::scan(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::scan(...)"
    );
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  MPI_Scan(
    const_cast<char*>(sendBuffer),scanReducts,bytes,MPI_CHAR,op.mpi_op()
    ,*rawMpiComm_
    );
}


template<typename Ordinal>
void MpiComm<Ordinal>::send(
  const Ordinal bytes, const char sendBuffer[], const int destRank
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::send(...)"
    );
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    ! ( 0 <= destRank && destRank < size_ ), std::logic_error
    ,"Error, destRank = " << destRank << " is not < 0 or is not"
    " in the range [0,"<<size_-1<<"]!"
    );
#endif // TEUCHOS_DEBUG
#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::send(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP
  MPI_Send(
    const_cast<char*>(sendBuffer),bytes,MPI_CHAR,destRank,tag_,*rawMpiComm_
    );
  // ToDo: What about error handling???
}


template<typename Ordinal>
void MpiComm<Ordinal>::readySend(
  const ArrayView<const char> &sendBuffer,
  const int destRank
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::readySend(...)"
    );
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    ! ( 0 <= destRank && destRank < size_ ), std::logic_error
    ,"Error, destRank = " << destRank << " is not < 0 or is not"
    " in the range [0,"<<size_-1<<"]!"
    );
#endif // TEUCHOS_DEBUG
#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::readySend(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP
  MPI_Rsend(
    const_cast<char*>(sendBuffer.getRawPtr()),sendBuffer.size(),MPI_CHAR,destRank,tag_,*rawMpiComm_
    );
  // ToDo: What about error handling???
}


template<typename Ordinal>
int MpiComm<Ordinal>::receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::receive(...)"
    );
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    sourceRank >=0 && !(sourceRank < size_), std::logic_error
    ,"Error, sourceRank = " << sourceRank << " is not < 0 or is not"
    " in the range [0,"<<(size_-1)<<"]!"
    );
#endif // TEUCHOS_DEBUG
  MPI_Status status;
  MPI_Recv(
    recvBuffer,bytes,MPI_CHAR
    ,sourceRank >= 0 ? sourceRank : MPI_ANY_SOURCE
    ,tag_,*rawMpiComm_
    ,&status
    );
#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::receive(...)"
      ,"recvBuffer", bytes, recvBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP
  return status.MPI_SOURCE;
  // ToDo: What about error handling???
}


template<typename Ordinal>
RCP<CommRequest> MpiComm<Ordinal>::isend(
  const ArrayView<const char> &sendBuffer,
  const int destRank
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::isend(...)"
    );
#ifdef TEUCHOS_DEBUG
  assertRank(destRank, "destRank");
#endif // TEUCHOS_DEBUG
  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  MPI_Isend(
    const_cast<char*>(sendBuffer.getRawPtr()), sendBuffer.size(), MPI_CHAR, destRank,
    tag_, *rawMpiComm_, &rawMpiRequest );
  return mpiCommRequest(rawMpiRequest);
  // ToDo: What about MPI error handling???
}


template<typename Ordinal>
RCP<CommRequest> MpiComm<Ordinal>::ireceive(
  const ArrayView<char> &recvBuffer,
  const int sourceRank
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::ireceive(...)"
    );
#ifdef TEUCHOS_DEBUG
  assertRank(sourceRank, "sourceRank");
#endif // TEUCHOS_DEBUG
  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  MPI_Irecv(
    const_cast<char*>(recvBuffer.getRawPtr()), recvBuffer.size(), MPI_CHAR, sourceRank,
    tag_, *rawMpiComm_, &rawMpiRequest );
  return mpiCommRequest(rawMpiRequest);
  // ToDo: What about MPI error handling???
}


template<typename Ordinal>
void MpiComm<Ordinal>::waitAll(
  const ArrayView<RCP<CommRequest> > &requests
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::waitAll(...)"
    );
  const int count = requests.size();
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( requests.size() == 0 );
#endif
  
  Array<MPI_Request> rawMpiRequests(count, MPI_REQUEST_NULL);
  for (int i = 0; i < count; ++i) {
    RCP<CommRequest> &request = requests[i];
    if (!is_null(request)) {
      const RCP<MpiCommRequest> mpiCommRequest =
        rcp_dynamic_cast<MpiCommRequest>(request);
      rawMpiRequests[i] = mpiCommRequest->releaseRawMpiRequest();
    }
    // else already null
    request = null;
  }

  Array<MPI_Status> rawMpiStatuses(count);
  MPI_Waitall( count, rawMpiRequests.getRawPtr(), rawMpiStatuses.getRawPtr() );
  // ToDo: We really should check the status?

}


template<typename Ordinal>
void MpiComm<Ordinal>::wait(
  const Ptr<RCP<CommRequest> > &request
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::wait(...)"
    );
  if (is_null(*request)) {
    return; // Nothing to wait on ...
  }
  const RCP<MpiCommRequest> mpiCommRequest =
    rcp_dynamic_cast<MpiCommRequest>(*request);
  MPI_Request rawMpiRequest = mpiCommRequest->releaseRawMpiRequest();
  MPI_Status status;
  MPI_Wait( &rawMpiRequest, &status );
  // ToDo: We really should check the status?
  *request = null;
}


// Overridden from Describable


template<typename Ordinal>
std::string MpiComm<Ordinal>::description() const
{
  std::ostringstream oss;
  oss
    << typeName(*this)
    << "{"
    << "size="<<size_
    << ",rank="<<rank_
    << ",rawMpiComm="<<static_cast<MPI_Comm>(*rawMpiComm_)
    <<"}";
  return oss.str();
}


#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal>
bool MpiComm<Ordinal>::show_dump = false;
#endif


// private


template<typename Ordinal>
void MpiComm<Ordinal>::assertRank(const int rank, const std::string &rankName) const
{
  TEST_FOR_EXCEPTION(
    ! ( 0 <= rank && rank < size_ ), std::logic_error
    ,"Error, "<<rankName<<" = " << rank << " is not < 0 or is not"
    " in the range [0,"<<size_-1<<"]!"
    );
}


} // namespace Teuchos


template<typename Ordinal>
Teuchos::RCP<Teuchos::MpiComm<Ordinal> >
Teuchos::createMpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  if( rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL )
    return rcp(new MpiComm<Ordinal>(rawMpiComm));
  return Teuchos::null;
}


#endif // TEUCHOS_MPI_COMM_HPP

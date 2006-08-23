// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
#include "mpi.h"

//#define TEUCHOS_MPI_COMM_DUMP

#ifdef TEUCHOS_MPI_COMM_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif

namespace Teuchos {

#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal>
void dumpCharBuffer(
  const std::string &funcName, const std::string &buffName
  ,const Ordinal bytes, const char buff[]
  )
{
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out
    << "\n" << funcName << "::" << buffName << ":";
  tab.incrTab();
  for( Ordinal i = 0; i < bytes; ++i ) {
    *out << buffName << "[" << i << "] = " << buff[i] << "\n";
  }
  *out << "\n";
}
#endif // TEUCHOS_MPI_COMM_DUMP

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
    const RefCountPtr<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
    );

  /** \brief Return the embedded wrapped opaque <tt>MPI_Comm</tt> object. */
  RefCountPtr<const OpaqueWrapper<MPI_Comm> > getRawMpiComm() const;

  //@}

  //! @name Overridden from Comm 
  //@{

  /** \brief . */
  int getRank() const;
  /** \brief . */
  int getSize() const;
  /** \brief . */
  void barrier() const;
  /** \brief . */
  void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const;
  /** \brief . */
  void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const;
  /** \brief . */
  void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const;
  /** \brief . */
  void reduceAllAndScatter(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvCounts[], const Ordinal blockSize, char myGlobalReducts[]
    ) const;
  /** \brief . */
	void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const;

  //@}

  //! @name Overridden from Describable 
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  static int const minTag_ = 26000; // These came from Teuchos::MpiComm???
  static int const maxTag_ = 26099; // ""
  static int tagCounter_;

  RefCountPtr<const OpaqueWrapper<MPI_Comm> > rawMpiComm_;
  int                                         rank_;
  int                                         size_;
  int                                         tag_;

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
RefCountPtr<MpiComm<Ordinal> >
createMpiComm(
  const RefCountPtr<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  );

// ////////////////////////
// Implementations

// Static members

template<typename Ordinal>
int MpiComm<Ordinal>::tagCounter_ = MpiComm<Ordinal>::minTag_;

// Constructors

template<typename Ordinal>
MpiComm<Ordinal>::MpiComm(
  const RefCountPtr<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
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
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::gatherAll(...)"
    );
  TEST_FOR_EXCEPT(!(sendBytes*size_==recvBytes));
  MPI_Allgather(
    const_cast<char *>(sendBuffer),sendBytes,MPI_CHAR
    ,recvBuffer,sendBytes,MPI_CHAR
    ,*rawMpiComm_
    );
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
  ,const Ordinal recvCounts[], const Ordinal blockSize, char myGlobalReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::reduceAllAndScatter(...)"
    );
#ifdef TEUCHOS_DEBUG
  Ordinal sumRecvBytes = 0;
  for( Ordinal i = 0; i < size_; ++i )
    sumRecvBytes += recvCounts[i];
  sumRecvBytes *= blockSize;
  TEST_FOR_EXCEPT(!(sumRecvBytes==sendBytes));
#endif // TEUCHOS_DEBUG
  WorkspaceStore* wss = get_default_workspace_store().get();
  // Create a of recvCount[] if Ordinal!=int
  const bool Ordinal_is_int = typeid(int)==typeid(Ordinal);
  Workspace<int> _recvCounts(wss,Ordinal_is_int?0:size_);
  const int *int_recvCounts = 0;
  if(Ordinal_is_int) {
    int_recvCounts = reinterpret_cast<const int*>(recvCounts);
    // Note: We must do an reinterpet cast since this must
    // compile even if it is not executed.  I could implement
    // code that would not need to do this using template
    // conditionals but I don't want to bother.
  }
  else {
    std::copy(recvCounts,recvCounts+size_,&_recvCounts[0]);
    int_recvCounts = &_recvCounts[0];
  }
  MPI_Datatype _chars_type;
  MPI_Type_contiguous(blockSize,MPI_CHAR,&_chars_type);
  RefCountPtr<const OpaqueWrapper<MPI_Datatype> >
    chars_type = opaqueWrapper(_chars_type,MPI_Type_free);
  // Perform the operation
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  MPI_Reduce_scatter(
    const_cast<char*>(sendBuffer), myGlobalReducts
    ,const_cast<int*>(int_recvCounts)
    ,*chars_type
    ,op.mpi_op()
    ,*rawMpiComm_
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
    dumpCharBuffer<Ordinal>(
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
    dumpCharBuffer<Ordinal>(
      "Teuchos::MpiComm<Ordinal>::receive(...)"
      ,"recvBuffer", bytes, recvBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP
  return status.MPI_SOURCE;
  // ToDo: What about error handling???
}

// Overridden from Describable

template<typename Ordinal>
std::string MpiComm<Ordinal>::description() const
{
  std::ostringstream oss;
  oss
    << "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">"
    << "{"
    << "size="<<size_
    << ",rawMpiComm="<<static_cast<MPI_Comm>(*rawMpiComm_)
    <<"}";
  return oss.str();
}

#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal>
bool MpiComm<Ordinal>::show_dump = false;
#endif // TEUCHOS_MPI_COMM_DUMP

} // namespace Teuchos

template<typename Ordinal>
Teuchos::RefCountPtr<Teuchos::MpiComm<Ordinal> >
Teuchos::createMpiComm(
  const RefCountPtr<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  if( rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL )
    return rcp(new MpiComm<Ordinal>(rawMpiComm));
  return Teuchos::null;
}

#endif // TEUCHOS_MPI_COMM_HPP

#include "Teuchos_MPIComm.hpp"


using namespace Teuchos;

namespace Teuchos
{
	const int MPIComm::INT = 1;
	const int MPIComm::FLOAT = 2;
	const int MPIComm::DOUBLE = 3;
	const int MPIComm::CHAR = 4;

	const int MPIComm::SUM = 5;
	const int MPIComm::MIN = 6;
	const int MPIComm::MAX = 7;
	const int MPIComm::PROD = 8;
}


MPIComm::MPIComm()
	:
#ifdef HAVE_MPI
	comm_(MPI_COMM_WORLD),
#endif
	nProc_(0), myRank_(0)
{
	init();
}

#ifdef HAVE_MPI
MPIComm::MPIComm(MPI_Comm comm)
	: comm_(comm), nProc_(0), myRank_(0)
{
	init();
}
#endif

void MPIComm::init()
{
#ifdef HAVE_MPI

	errCheck(MPI_Comm_rank(comm_, &myRank_), "Comm_rank");
					
	errCheck(MPI_Comm_size(comm_, &nProc_), "Comm_size");
	
#else
	nProc_ = 1;
	myRank_ = 0;
#endif
}

#ifdef USE_MPI_GROUPS /* we're ignoring groups for now */

MPIComm::MPIComm(const MPIComm& parent, const MPIGroup& group)
	:
#ifdef HAVE_MPI
	comm_(MPI_COMM_WORLD), 
#endif
	nProc_(0), myRank_(0)
{
#ifdef HAVE_MPI
	if (group.getNProc()==0)
		{
			rank_ = -1;
			nProc_ = 0;
		}
	else if (parent.containsMe())
		{
			MPI_Comm parentComm = parent.comm_;
			MPI_Group newGroup = group.group_;
			
			errCheck(MPI_Comm_create(parentComm, newGroup, &comm_), 
							 "Comm_create");
			
			if (group.containsProc(parent.getRank()))
				{
					errCheck(MPI_Comm_rank(comm_, &rank_), "Comm_rank");
					
					errCheck(MPI_Comm_size(comm_, &nProc_), "Comm_size");
				}
			else
				{
					rank_ = -1;
					nProc_ = -1;
					return;
				}
		}
	else
		{
			rank_ = -1;
			nProc_ = -1;
		}
#endif
}

#endif /* USE_MPI_GROUPS */

MPIComm& MPIComm::world()
{
	static MPIComm w = MPIComm();
	return w;
}


void MPIComm::synchronize() const 
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		errCheck(::MPI_Barrier(comm_), "Barrier");
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allToAll(void* sendBuf, int sendCount, int sendType,
											 void* recvBuf, int recvCount, int recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = getDataType(sendType);
		MPI_Datatype mpiRecvType = getDataType(recvType);

		errCheck(::MPI_Alltoall(sendBuf, sendCount, mpiSendType,
														recvBuf, recvCount, mpiRecvType,
														comm_), "Alltoall");
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allToAllv(void* sendBuf, int* sendCount, 
												int* sendDisplacements, int sendType,
												void* recvBuf, int* recvCount, 
												int* recvDisplacements, int recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = getDataType(sendType);
		MPI_Datatype mpiRecvType = getDataType(recvType);
		
		errCheck(::MPI_Alltoallv(sendBuf, sendCount, sendDisplacements, mpiSendType,
														 recvBuf, recvCount, recvDisplacements, mpiRecvType,
														 comm_), "Alltoallv");
	}
	//mutex_.unlock();
#endif
}

void MPIComm::gather(void* sendBuf, int sendCount, int sendType,
										 void* recvBuf, int recvCount, int recvType,
										 int root) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = getDataType(sendType);
		MPI_Datatype mpiRecvType = getDataType(recvType);
		
		errCheck(::MPI_Gather(sendBuf, sendCount, mpiSendType,
													recvBuf, recvCount, mpiRecvType,
													root, comm_), "Gather");
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allGather(void* sendBuf, int sendCount, int sendType,
												void* recvBuf, int recvCount, 
												int recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = getDataType(sendType);
		MPI_Datatype mpiRecvType = getDataType(recvType);
		
		errCheck(::MPI_Allgather(sendBuf, sendCount, mpiSendType,
														 recvBuf, recvCount, 
														 mpiRecvType, comm_), 
						 "AllGather");
	}
	//mutex_.unlock();
#endif
}


void MPIComm::allGatherv(void* sendBuf, int sendCount, int sendType,
												 void* recvBuf, int* recvCount, 
												 int* recvDisplacements,
												 int recvType) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiSendType = getDataType(sendType);
		MPI_Datatype mpiRecvType = getDataType(recvType);
		
		errCheck(::MPI_Allgatherv(sendBuf, sendCount, mpiSendType,
															recvBuf, recvCount, recvDisplacements,
															mpiRecvType, 
															comm_), 
						 "AllGatherv");
	}
	//mutex_.unlock();
#endif
}


void MPIComm::bcast(void* msg, int length, int type, int src) const
{
#ifdef HAVE_MPI
	//mutex_.lock();
	{
		MPI_Datatype mpiType = getDataType(type);
		errCheck(::MPI_Bcast(msg, length, mpiType, src, 
												 comm_), "Bcast");
	}
	//mutex_.unlock();
#endif
}

void MPIComm::allReduce(void* input, void* result, int inputCount, 
												int type, int op) const
{
#ifdef HAVE_MPI

	//mutex_.lock();
	{
		MPI_Op mpiOp = getOp(op);
		MPI_Datatype mpiType = getDataType(type);
		
		errCheck(::MPI_Allreduce(input, result, inputCount, mpiType,
														 mpiOp, comm_), 
						 "Allreduce");
	}
	//mutex_.unlock();
#endif
}


#ifdef HAVE_MPI

MPI_Datatype MPIComm::getDataType(int type)
{
  TEST_FOR_EXCEPTION( !(type == INT || type==FLOAT 
                        || type==DOUBLE || type==CHAR),
                      range_error,
                      "invalid type " << type << " in MPIComm::getDataType");

  if(type == INT) return MPI_INT;
  if(type == FLOAT) return MPI_FLOAT;
  if(type == DOUBLE) return MPI_DOUBLE;
  
  return MPI_CHAR;
}


void MPIComm::errCheck(int errCode, const string& methodName)
{
  TEST_FOR_EXCEPTION(errCode != 0, runtime_error,
                     "MPI function MPI_" << methodName 
                     << " returned error code=" << errCode);
}

MPI_Op MPIComm::getOp(int op)
{

  TEST_FOR_EXCEPTION( !(op == SUM || op==MAX 
                        || op==MIN || op==PROD),
                      range_error,
                      "invalid operator " 
                      << op << " in MPIComm::getOp");

  if( op == SUM) return MPI_SUM;
  else if( op == MAX) return MPI_MAX;
  else if( op == MIN) return MPI_MIN;
  return MPI_PROD;
}

#endif

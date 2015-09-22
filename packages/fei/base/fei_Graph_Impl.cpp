/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_Graph_Impl.hpp>
#include <fei_EqnComm.hpp>
#include <fei_CommUtils.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_VectorSpace.hpp>

#undef fei_file
#define fei_file "fei_Graph_Impl.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
fei::Graph_Impl::Graph_Impl(MPI_Comm comm, int firstLocalRow, int lastLocalRow)
  : localGraphData_(NULL),
    remoteGraphData_(),
    eqnComm_(),
    firstLocalRow_(firstLocalRow),
    lastLocalRow_(lastLocalRow),
    localProc_(0),
    numProcs_(1),
    comm_(comm)
{
  localProc_ = fei::localProc(comm_);
  numProcs_  = fei::numProcs(comm_);
  //for remoteGraphData_, we don't know what the range of row-numbers will
  //be, so we'll just construct it with -1,-1
  remoteGraphData_.resize(numProcs_);
  for(int p=0; p<numProcs_; ++p) {
    remoteGraphData_[p] = new remote_table_type(-1, -1);
  }
  eqnComm_.reset(new fei::EqnComm(comm_, lastLocalRow-firstLocalRow+1));
  localGraphData_       = new table_type(firstLocalRow_, lastLocalRow_);
}

//----------------------------------------------------------------------------
fei::Graph_Impl::~Graph_Impl()
{
  delete localGraphData_;
  for(unsigned i=0; i<remoteGraphData_.size(); ++i) {
    delete remoteGraphData_[i];
  }
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::addIndices(int row, int len, const int* indices)
{
  if (row < 0) {
    return(-1);
  }

  if (row < firstLocalRow_ || row > lastLocalRow_) {
    int p = eqnComm_->getOwnerProc(row);
    remoteGraphData_[p]->addIndices(row, len, indices);
  }
  else localGraphData_->addIndices(row, len, indices);

  return(0);
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::addSymmetricIndices(int numIndices, int* indices,
					bool diagonal)
{
  if (diagonal) {
    addDiagonals(numIndices, indices);
    return(0);
  }

  bool all_local = true;
  if (numProcs_ > 1) {
    for(int i=0; i<numIndices; ++i) {
      if (indices[i] < 0) {
	return(-1);
      }

      bool local = true;
      if (indices[i] < firstLocalRow_ || indices[i] > lastLocalRow_) {
	  all_local = false;
	  local = false;
      }

      if (!local) {
        int p = eqnComm_->getOwnerProc(indices[i]);
	remoteGraphData_[p]->addIndices(indices[i], numIndices, indices);
      }
    }
  }

  if (all_local) {
    localGraphData_->addIndices(numIndices, indices,
				  numIndices, indices);
  }
  else {
    for(int i=0; i<numIndices; ++i) {
      if (indices[i] >= firstLocalRow_ && indices[i] <= lastLocalRow_) {
	  localGraphData_->addIndices(indices[i], numIndices, indices);
      }
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
void fei::Graph_Impl::addDiagonals(int numIndices, int* indices)
{
  bool all_local = true;
  int i;
  if (numProcs_ > 1) {
    for(i=0; i<numIndices; ++i) {
      int ind = indices[i];
      if (ind < 0) {
        throw std::runtime_error("fei::Graph_Impl::addDiagonals given negative index");
      }

      bool local = true;
      if (ind < firstLocalRow_ || ind > lastLocalRow_) {
	  all_local = false;
	  local = false;
      }

      if (!local) {
        int p=eqnComm_->getOwnerProc(ind);
	remoteGraphData_[p]->addIndices(ind, 1, &ind);
      }
    }
  }

  if (all_local) {
    localGraphData_->addDiagonals(numIndices, indices);
  }
  else {
    for(i=0; i<numIndices; ++i) {
      int ind = indices[i];
      if (ind >= firstLocalRow_ && ind <= lastLocalRow_) {
	  localGraphData_->addIndices(ind, 1, &ind);
      }
    }
  }
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::writeLocalGraph(FEI_OSTREAM& os, bool debug,
				    bool prefixLinesWithPoundSign)
{
  if (localGraphData_ == NULL) {
    if (prefixLinesWithPoundSign) {
      os << "# fei::Graph_Impl::writeLocalGraph numRows: 0"<<FEI_ENDL;
    }
    else {
      os << "0" << FEI_ENDL;
    }
    return(0);
  }

  if (prefixLinesWithPoundSign) {
    os << "# fei::Graph_Impl::writeLocalGraph numRows: ";
  }
  os << localGraphData_->getMap().size() <<FEI_ENDL;

  if (prefixLinesWithPoundSign) {
    writeToStream(*localGraphData_, os, "# ");
  }
  else {
    writeToStream(*localGraphData_, os, NULL);
  }

  return( 0 );
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::writeRemoteGraph(FEI_OSTREAM& os)
{
  for(unsigned p=0; p<remoteGraphData_.size(); ++p) {
    writeToStream(*remoteGraphData_[p], os, "# ");
  }

  return( 0 );
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::gatherFromOverlap()
{
  if (numProcs_ == 1) return(0);

#ifndef FEI_SER
  //this function gathers shared-but-not-owned data onto the
  //owning processors.

  //iterate the remoteGraphData_ array and create a list of processors
  //that we will be sending data to. (processors which own graph rows
  //that we share.)
  std::vector<int> sendProcs;
  for(unsigned i=0; i<remoteGraphData_.size(); ++i) {
    if ((int)i == localProc_) continue;
    if (remoteGraphData_[i] != NULL) {
      if (remoteGraphData_[i]->getMap().size() == 0) {
        continue;
      }

      sendProcs.push_back((int)i);
    }
  }

  //now we can find out which procs we'll be receiving from.
  std::vector<int> recvProcs;
  fei::mirrorProcs(comm_, sendProcs, recvProcs);

  //next we'll declare arrays to receive into.
  std::vector<std::vector<int> > recv_ints(recvProcs.size());

  //...and an array for the sizes of the recv buffers:
  std::vector<int> recv_sizes(recvProcs.size());
  std::vector<MPI_Request> mpiReqs(recvProcs.size());
  std::vector<MPI_Status> mpiStatuses(recvProcs.size());

  int tag1 = 11113;

  unsigned offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    MPI_Irecv(&recv_sizes[i], 1, MPI_INT, recvProcs[i],
              tag1, comm_, &mpiReqs[i]);
  }

  //now we'll pack our to-be-sent data into buffers, and send the
  //sizes to the receiving procs:
  std::vector<std::vector<int> > send_ints(sendProcs.size());

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];

    fei::packRaggedTable(*(remoteGraphData_[proc]), send_ints[i]);

    int isize = send_ints[i].size();

    MPI_Send(&isize, 1, MPI_INT, proc, tag1, comm_);
  }

  if (mpiReqs.size() > 0) {
    MPI_Waitall(mpiReqs.size(), &mpiReqs[0], &mpiStatuses[0]);
  }

  //now resize our recv buffers, and post the recvs.
  for(size_t i=0; i<recvProcs.size(); ++i) {
    int intsize = recv_sizes[i];

    recv_ints[i].resize(intsize);

    MPI_Irecv(&(recv_ints[i][0]), intsize, MPI_INT, recvProcs[i],
              tag1, comm_, &mpiReqs[i]);
  }

  //now send our packed buffers.
  for(size_t i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];

    MPI_Send(&(send_ints[i][0]), send_ints[i].size(), MPI_INT,
             proc, tag1, comm_);
  }

  if (mpiReqs.size() > 0) {
    MPI_Waitall(mpiReqs.size(), &mpiReqs[0], &mpiStatuses[0]);
  }

  for(unsigned i=0; i<recvProcs.size(); ++i) {
    std::vector<int> recvdata = recv_ints[i];
    int numRows = recvdata[0];
    int* rowNumbers = &recvdata[1];
    int* rowLengths = rowNumbers+numRows;
    int* packedCols = rowLengths+numRows;
    offset = 0;
    for(int r=0; r<numRows; ++r) {
      addIndices(rowNumbers[r], rowLengths[r], &packedCols[offset]);
      offset += rowLengths[r];
    }
  }

#endif
  return(0);
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::getLocalRowLength(int row)
{
  table_row_type* colIndices = localGraphData_->getRow(row);
  if (colIndices == NULL) {
    return(-1);
  }

  return(colIndices->size());
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::getNumLocalRows()
{
  int numLocalRows = localGraphData_->getMap().size();

  return(numLocalRows);
}

//----------------------------------------------------------------------------
int fei::Graph_Impl::getNumLocalNonzeros() const
{
  int numNonzeros = fei::countNonzeros(*localGraphData_);

  return(numNonzeros);
}


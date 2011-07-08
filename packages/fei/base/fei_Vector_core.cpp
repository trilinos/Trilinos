/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_fstream.hpp"

#include "fei_Vector_core.hpp"
#include "fei_VectorSpace.hpp"
#include "fei_CSVec.hpp"
#include "snl_fei_RecordCollection.hpp"
#include "fei_TemplateUtils.hpp"
#include "fei_impl_utils.hpp"

#include <fstream>
#include <sstream>
#undef fei_file
#define fei_file "fei_Vector_core.cpp"

#include "fei_ErrMacros.hpp"

fei::Vector_core::Vector_core(fei::SharedPtr<fei::VectorSpace> vecSpace,
                              int numLocalEqns)
  : eqnComm_(),
    vecSpace_(vecSpace),
    comm_(vecSpace->getCommunicator()),
    firstLocalOffset_(0),
    lastLocalOffset_(0),
    numLocal_(0),
    work_indices_(),
    work_indices2_(),
    haveFEVector_(false),
    remotelyOwned_(),
    sendProcs_(),
    recvProcs_(),
    recv_sizes_(),
    sendRecvProcsNeedUpdated_(true),
    overlapAlreadySet_(false),
    dbgprefix_("Vcore: ")
{
  eqnComm_.reset(new fei::EqnComm(comm_,numLocalEqns));
  remotelyOwned_.resize(fei::numProcs(comm_));
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    remotelyOwned_[i] = new CSVec;
  }

  const std::vector<int>& offsets = eqnComm_->getGlobalOffsets();
  firstLocalOffset_ = offsets[fei::localProc(comm_)];
  lastLocalOffset_ = offsets[fei::localProc(comm_)+1] - 1;
  numLocal_ = lastLocalOffset_ - firstLocalOffset_ + 1;

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os<<dbgprefix_<<" ctor firstLocal="<<firstLocalOffset_<<", lastLocal="
     <<lastLocalOffset_<<FEI_ENDL;
  }
}

fei::Vector_core::~Vector_core()
{
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    delete remotelyOwned_[i];
  }
}

void fei::Vector_core::setOverlap(int numRemoteEqns,
                                  const int* remoteEqns)
{
  if (numRemoteEqns == 0 && remoteEqns == NULL) {
    if (overlapAlreadySet_) return;
  }

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"setOverlap"<<FEI_ENDL;
  }

  int local_proc = fei::localProc(comm_);

  if (numRemoteEqns != 0 && remoteEqns != NULL) {
    for(int i=0; i<numRemoteEqns; ++i) {
      int proc = eqnComm_->getOwnerProc(remoteEqns[i]);
      if (proc == local_proc) continue;

      fei::add_entry(*remotelyOwned_[proc], remoteEqns[i], 0.0);
    }
  }
  else {
    std::vector<int> eqns;
    vecSpace_->getIndices_SharedAndOwned(eqns);

    for(size_t i=0; i<eqns.size(); ++i) {
      int proc = eqnComm_->getOwnerProc(eqns[i]);
      if (proc == local_proc) continue;

      fei::add_entry(*remotelyOwned_[proc], eqns[i], 0.0);
    }
  }

  overlapAlreadySet_ = true;
}

int fei::Vector_core::scatterToOverlap()
{
  if (fei::numProcs(comm_) == 1 || haveFEVector()) {
    return(0);
  }

#ifndef FEI_SER
  if (!overlapAlreadySet_) {
    setOverlap();
  }

  //...and now the overlap is whatever is in our remotelyOwned_ vectors.

//  if (sendRecvProcsNeedUpdated_) {
    //first find out which procs we'll be receiving from.
    std::vector<int> recvProcs;
    for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
      if ((int)i == fei::localProc(comm_)) continue;
      if (remotelyOwned_[i]->size() == 0) continue;
  
      recvProcs.push_back((int)i);
    }
  
    //find out the send-procs.
    std::vector<int> sendProcs;
    fei::mirrorProcs(comm_, recvProcs, sendProcs);
    sendRecvProcsNeedUpdated_ = false;
//  }

  //declare arrays to send from, and corresponding sizes
  std::vector<std::vector<int> > send_ints(sendProcs.size());
  std::vector<std::vector<double> > send_doubles(sendProcs.size());
  std::vector<int> send_sizes(sendProcs.size());

  std::vector<MPI_Request> mpiReqs(sendProcs.size()+recvProcs.size());
  std::vector<MPI_Status> mpiStatuses(sendProcs.size()+recvProcs.size());
  int tag1 = 11111;
  int tag2 = 11112;

  //first, the procs we're going to send to, have to let us know
  //how much data we're supposed to send. So we have to receive
  //sizes and then indices from the "send"-procs.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    MPI_Irecv(&send_sizes[i], 1, MPI_INT, sendProcs[i],
              tag1, comm_, &mpiReqs[i]);
  }

  //now we'll send the sizes of our remotely-owned data to the
  //procs that we will be receiving the data from, and also the
  //indices that we want to receive.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];

    int size = remotelyOwned_[proc]->size();
    MPI_Send(&size, 1, MPI_INT, proc, tag1, comm_);
  }
 
  MPI_Waitall(sendProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

  //now resize our send_ints and send_doubles arrays, and post the recvs
  //for indices that we're supposed to pack.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];
    int size = send_sizes[i];
    send_ints[i].resize(size);
    MPI_Irecv(&(send_ints[i][0]), size, MPI_INT, proc, tag1,
              comm_, &mpiReqs[i]);
    send_doubles[i].resize(size);
  }

  //now send the indices that we want to receive data for.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];
    int size = remotelyOwned_[proc]->size();
    int* indices = &(remotelyOwned_[proc]->indices())[0];
    MPI_Send(indices, size, MPI_INT, proc, tag1, comm_);
  }

  MPI_Waitall(sendProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

  //now post our recvs.
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int proc = recvProcs[i];
    int size = remotelyOwned_[proc]->size();
    double* coefs = &(remotelyOwned_[proc]->coefs())[0];
    MPI_Irecv(coefs, size, MPI_DOUBLE, proc, tag2, comm_, &mpiReqs[i]);
  }

  //now pack and send the coefs that the other procs need from us.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];

    int num = send_sizes[i];
    int err = copyOutOfUnderlyingVector(num, &(send_ints[i][0]),
                                        &(send_doubles[i][0]), 0);
    if (err != 0) {
      FEI_COUT << "fei::Vector_core::scatterToOverlap ERROR getting data to send."<<FEI_ENDL;
      return(err);
    }

    MPI_Send(&(send_doubles[i][0]), num, MPI_DOUBLE, proc,
             tag2, comm_);
  }

  MPI_Waitall(recvProcs.size(), &mpiReqs[0], &mpiStatuses[0]);

#endif  //#ifndef FEI_SER

  return(0);
}

int fei::Vector_core::copyOut(int numValues,
				  const int* indices,
				  double* values,
				  int vectorIndex) const
{
  const std::vector<CSVec*>& remote = remotelyOwned();

  for(int i=0; i<numValues; ++i) {
    int ind = indices[i];

    int local = ind - firstLocalOffset_;
    if (local < 0 || local >= numLocal_) {
      if (ind < 0) {
	continue;
      }

      int proc = eqnComm_->getOwnerProc(ind);

      int insertPoint = -1;
      int idx = fei::binarySearch(ind, remote[proc]->indices(), insertPoint);
      if (idx < 0) {
	fei::console_out() << "fei::Vector_core::copyOut: proc " << fei::localProc(comm_)
	     << ", index " << ind << " not in remotelyOwned_ vec object for proc "
            <<proc<<FEI_ENDL;
	ERReturn(-1);
      }
      else {
	values[i] = remote[proc]->coefs()[idx];
      }
    }
    else {
      CHK_ERR( copyOutOfUnderlyingVector(1, &ind, &(values[i]), vectorIndex) );
    }
  }

  return(0);
}

int fei::Vector_core::giveToVector(int numValues,
				       const int* indices,
				       const double* values,
				       bool sumInto,
				       int vectorIndex)
{
  std::vector<CSVec*>& remote = remotelyOwned();

  for(int i=0; i<numValues; ++i) {
    int ind = indices[i];
    double val = values[i];

    if (ind < 0) {
//      throw std::runtime_error("negative index not allowed");
      //preservation of existing behavior: certain Sierra scenarios
      //involve passing negative indices for positions that should be
      //ignored... so we'll continue rather than throwing.
      continue;
    }
    int proc = eqnComm_->getOwnerProc(ind);
    int local = ind - firstLocalOffset_;
    if (local < 0 || local >= numLocal_) {
      if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
        FEI_OSTREAM& os = *output_stream_;
        os << dbgprefix_<<"giveToVector remote["<<proc<<"]("
         <<ind<<","<<val<<")"<<FEI_ENDL;
      }

      if (sumInto) {
	fei::add_entry( *remote[proc], ind, val);
      }
      else {
	fei::put_entry( *remote[proc], ind, val);
      }
    }
    else {
      int err = giveToUnderlyingVector(1, &ind, &val, sumInto, vectorIndex);
      if (err != 0) {
	fei::console_out() << "giveToVector sumIn ERROR, ind: " << ind
	     << ", val: " << val << FEI_ENDL;
	ERReturn(-1);
      }
    }
  }

  return(0);
}

int fei::Vector_core::assembleFieldData(int fieldID,
					    int idType,
					    int numIDs,
					    const int* IDs,
					    const double* data,
					    bool sumInto,
					    int vectorIndex)
{
  if (vecSpace_.get() == NULL) ERReturn(-1);

  int fieldSize = vecSpace_->getFieldSize(fieldID);

  work_indices_.resize(numIDs*fieldSize);
  int* indicesPtr = &work_indices_[0];

  CHK_ERR( vecSpace_->getGlobalIndices(numIDs, IDs, idType, fieldID,
					indicesPtr) );

  CHK_ERR( giveToVector(numIDs*fieldSize, indicesPtr, data, sumInto, vectorIndex) );

  return(0);
}

void pack_send_buffers(const std::vector<int>& sendProcs,
                       const std::vector<fei::CSVec*>& remotelyOwned,
                       std::vector<std::vector<char> >& send_chars)
{
  for(size_t i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];
    fei::impl_utils::pack_indices_coefs(remotelyOwned[proc]->indices(),
                       remotelyOwned[proc]->coefs(), send_chars[i]);

    fei::set_values(*remotelyOwned[proc], 0.0);
  }
}

int fei::Vector_core::gatherFromOverlap(bool accumulate)
{
  if (fei::numProcs(comm_) == 1 || haveFEVector()) {
    return(0);
  }

//int localproc = fei::localProc(comm_);
//std::ostringstream osstr;
//osstr << "proc."<<localproc;
//std::string str=osstr.str();
//static std::ofstream ofs(str.c_str());

#ifndef FEI_SER
  //declare arrays to hold the buffers we'll be receiving.
  std::vector<std::vector<char> > recv_chars;
  std::vector<std::vector<char> > send_chars(sendProcs_.size());

  std::vector<MPI_Request> mpiReqs;
  std::vector<MPI_Request> mpiReqs2;
  std::vector<MPI_Status> mpiStatuses;
  int tag1 = 11111;

  if (sendRecvProcsNeedUpdated_) {
    sendProcs_.clear();
    vecSpace_->getSendProcs(sendProcs_);
//    //first create the list of procs we'll be sending to.
//    for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
//      if ((int)i == fei::localProc(comm_)) continue;
//      if (remotelyOwned_[i]->size() == 0) continue;
//  
//      sendProcs_.push_back(i);
//    }
  
    recvProcs_.clear();
    fei::mirrorProcs(comm_, sendProcs_, recvProcs_);
    sendRecvProcsNeedUpdated_ = false;
  }
  
//  int num_remote_sends_needed = 0;
//  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
//    if ((int)i==fei::localProc(comm_)) continue;
//    if (remotelyOwned_[i]->size() > 0) num_remote_sends_needed++;
//  }
//  if (num_remote_sends_needed != (int)sendProcs_.size()) {
//    std::cout << "proc " << fei::localProc(comm_)<<", num_remote_sends_needed="<<num_remote_sends_needed<<", sendProcs_.size()="<<sendProcs_.size()<<std::endl;
//  }

  send_chars.resize(sendProcs_.size());

  pack_send_buffers(sendProcs_, remotelyOwned_, send_chars);

  recv_sizes_.resize(recvProcs_.size());
  mpiReqs.resize(recvProcs_.size());
  mpiStatuses.resize(recvProcs_.size());
  mpiReqs2.resize(recvProcs_.size());
  mpiStatuses.resize(recvProcs_.size());
  recv_chars.resize(recvProcs_.size());


  //post the recvs for the sizes.
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    int proc = recvProcs_[i];
    MPI_Irecv(&recv_sizes_[i], 1, MPI_INT, proc,
              tag1, comm_, &mpiReqs[i]);
  }

  //send the sizes of data we'll be sending.
  for(unsigned i=0; i<sendProcs_.size(); ++i) {
    int proc = sendProcs_[i];
    int size = send_chars[i].size();
//std::cout<<"proc "<<fei::localProc(comm_)<<" sending "<<size<<" chars ("<<remotelyOwned_[proc]->indices().size()<<" indices/coefs) to proc "<<proc<<" accum="<<accumulate<<std::endl;
    MPI_Send(&size, 1, MPI_INT, proc, tag1, comm_);
  }

  //now post the recvs for the data.
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    int index;
    MPI_Status status;
    MPI_Waitany(mpiReqs.size(), &mpiReqs[0], &index, &status);

    recv_chars[index].resize(recv_sizes_[index]);
    MPI_Irecv(&(recv_chars[index][0]), recv_sizes_[index], MPI_CHAR, recvProcs_[index],
              tag1, comm_, &mpiReqs2[index]);
  }

  //now send the outgoing data.
  for(size_t i=0; i<sendProcs_.size(); ++i) {
    int proc = sendProcs_[i];
    int size = send_chars[i].size();
    MPI_Send(&(send_chars[i][0]), size, MPI_CHAR, proc, tag1, comm_);
  }

  std::vector<int> indices;
  std::vector<double> coefs;
  int numRecvProcs = recvProcs_.size();
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    int index;
    MPI_Status status;
    MPI_Waitany(numRecvProcs, &mpiReqs2[0], &index, &status);
  }

  //now store the data we've received.
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    fei::impl_utils::unpack_indices_coefs(recv_chars[i], indices, coefs);
    int num = indices.size();
    if (num == 0) continue;
//ofs<<"proc "<<fei::localProc(comm_)<<" recvd "<<num<<"inds/coefs from proc "<<recvProcs_[i]<<" accum="<<accumulate<<std::endl;
//for(int ii=0; ii<num;++ii)ofs<<"   "<<indices[ii]<<" "<<coefs[ii]<<std::endl;
    int err = giveToUnderlyingVector(num, &(indices[0]),
                                     &(coefs[0]), accumulate, 0);
    if (err != 0) {
      FEI_COUT << "fei::Vector_core::gatherFromOverlap ERROR storing recvd data" << FEI_ENDL;
      return(err);
    }
  }

#endif  //#ifndef FEI_SER

  return(0);
}

int fei::Vector_core::copyOutFieldData(int fieldID,
					   int idType,
					   int numIDs,
					   const int* IDs,
					   double* data,
					   int vectorIndex)
{
  if (vecSpace_.get() == NULL) ERReturn(-1);

  int fieldSize = vecSpace_->getFieldSize(fieldID);

  if (haveFEVector_) {
    snl_fei::RecordCollection* collection = NULL;
    CHK_ERR( vecSpace_->getRecordCollection(idType, collection) );
    int nodeNumber;
    int dofOffset;
    int numInstances;
    int foffset;
    std::vector<int>& eqnNums = vecSpace_->getEqnNumbers();
    int* vspcEqnPtr = eqnNums.size() > 0 ? &eqnNums[0] : NULL;

    int offset = 0;
    for(int i=0; i<numIDs; ++i) {
      fei::Record<int>* node = collection->getRecordWithID(IDs[i]);
      if (node == NULL) {
        ERReturn(-1);
      }

      nodeNumber = node->getNumber();
      int* eqnNumbers = vspcEqnPtr+node->getOffsetIntoEqnNumbers();
      node->getFieldMask()->getFieldEqnOffset(fieldID, foffset, numInstances);
      dofOffset = eqnNumbers[foffset] - eqnNumbers[0];
      for(int j=0; j<fieldSize; ++j) {
	CHK_ERR( copyOut_FE(nodeNumber, dofOffset+j, data[offset++]));
      }
    }
  }
  else {
    work_indices_.resize(numIDs*fieldSize*2);
    int* indicesPtr = &work_indices_[0];

    CHK_ERR( vecSpace_->getGlobalIndices(numIDs, IDs, idType,
					 fieldID, indicesPtr) );

    CHK_ERR( copyOut(numIDs*fieldSize, indicesPtr, data) );
  }

  return(0);
}

int fei::Vector_core::writeToFile(const char* filename,
                                    bool matrixMarketFormat)
{
  int numProcs = fei::numProcs(comm_);
  int localProc =fei::localProc(comm_);

  double coef;

  static char mmbanner[] = "%%MatrixMarket matrix array real general";

  for(int p=0; p<numProcs; ++p) {
    fei::Barrier(comm_);
    if (p != localProc) continue;

    FEI_OFSTREAM* outFile = NULL;
    if (p==0) {
      outFile = new FEI_OFSTREAM(filename, IOS_OUT);
      FEI_OFSTREAM& ofref = *outFile;
      if (matrixMarketFormat) {
        ofref << mmbanner << FEI_ENDL;
        ofref << eqnComm_->getGlobalOffsets()[numProcs] << " 1" << FEI_ENDL;
      }
      else {
        ofref << eqnComm_->getGlobalOffsets()[numProcs] << FEI_ENDL;
      }
    }
    else outFile = new FEI_OFSTREAM(filename, IOS_APP);
    FEI_OFSTREAM& ofref = *outFile;
    ofref.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
    ofref.precision(13);

    for(int i=firstLocalOffset_; i<=lastLocalOffset_; ++i) {
      CHK_ERR( copyOut(1, &i, &coef) );
      if (matrixMarketFormat) {
        ofref << " " << coef << FEI_ENDL;
      }
      else {
        ofref << i << " " << coef << FEI_ENDL;
      }
    }
    
    delete outFile;
  } 
    
  return(0);
}

int fei::Vector_core::writeToStream(FEI_OSTREAM& ostrm,
					bool matrixMarketFormat)
{
  int numProcs = fei::numProcs(comm_);
  int local_proc =fei::localProc(comm_);

  double coef;

  static char mmbanner[] = "%%MatrixMarket matrix array real general";

  IOS_FMTFLAGS oldf = ostrm.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
  ostrm.precision(13);

  for(int proc=0; proc<numProcs; ++proc) {
    fei::Barrier(comm_);
    if (proc != local_proc) continue;

    if (proc==0) {
      if (matrixMarketFormat) {
	ostrm << mmbanner << FEI_ENDL;
	ostrm << eqnComm_->getGlobalOffsets()[numProcs] << " 1" << FEI_ENDL;
      }
      else {
	ostrm << eqnComm_->getGlobalOffsets()[numProcs] << FEI_ENDL;
      }
    }

    for(int p=0; p<local_proc; ++p) {
      for(size_t ii=0; ii<remotelyOwned_[p]->size(); ++ii) {
	if (matrixMarketFormat) {
	  ostrm << " " << remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
	else {
	  ostrm << " " << remotelyOwned_[p]->indices()[ii] << " "
		<< remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
      }
    }

    for(int i=firstLocalOffset_; i<=lastLocalOffset_; ++i) {
      CHK_ERR( copyOut(1, &i, &coef) );
      if (matrixMarketFormat) {
	ostrm << " " << coef << FEI_ENDL;
      }
      else {
	ostrm << " " << i << " " << coef << FEI_ENDL;
      }
    }

    for(int p=local_proc+1; p<numProcs; ++p) {
      for(size_t ii=0; ii<remotelyOwned_[p]->size(); ++ii) {
	if (matrixMarketFormat) {
	  ostrm << " " << remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
	else {
	  ostrm << " " << remotelyOwned_[p]->indices()[ii] << " "
		<< remotelyOwned_[p]->coefs()[ii] << FEI_ENDL;
	}
      }
    }
  }

  ostrm.setf(oldf, IOS_FLOATFIELD);

  return(0);
}


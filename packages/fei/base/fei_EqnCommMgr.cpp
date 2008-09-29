/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_defs.h>

#include <fei_EqnCommMgr.hpp>

#include <fei_CommUtils.hpp>

#include <fei_ProcEqns.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_SSVec.hpp>
#include <fei_SSMat.hpp>
#include <feiArray.hpp>
#include <fei_TemplateUtils.hpp>

#include <algorithm>

#undef fei_file
#define fei_file "EqnCommMgr.cpp"
#include <fei_ErrMacros.hpp>

//-----Constructor--------------------------------------------------------------
EqnCommMgr::EqnCommMgr(MPI_Comm comm, bool accumulate)
 : accumulate_(accumulate),
   localProc_(0),
   recvProcEqns_(NULL),
   exchangeIndicesCalled_(false),
   recvEqns_(NULL),
   solnValues_(0, 8000),
   sendProcEqns_(NULL),
   sendEqns_(NULL),
   sendEqnSoln_(0, 2000),
   essBCEqns_(NULL),
   comm_(comm)
{
  localProc_ = fei::localProc(comm_);
  recvEqns_ = new EqnBuffer();
  sendEqns_ = new EqnBuffer();
  recvProcEqns_ = new ProcEqns();
  sendProcEqns_ = new ProcEqns();

  essBCEqns_ = new EqnBuffer();
}

//-----CopyConstructor--------------------------------------------------------------
EqnCommMgr::EqnCommMgr(const EqnCommMgr& src)
 : accumulate_(src.accumulate_),
   localProc_(0),
   recvProcEqns_(NULL),
   exchangeIndicesCalled_(false),
   recvEqns_(NULL),
   solnValues_(0, 8000),
   sendProcEqns_(NULL),
   sendEqns_(NULL),
   sendEqnSoln_(0, 2000),
   essBCEqns_(NULL),
   comm_(src.comm_)
{
  *this = src;
}

//-----Destructor---------------------------------------------------------------
EqnCommMgr::~EqnCommMgr()
{
   delete recvEqns_;
   delete sendEqns_;
   delete recvProcEqns_;
   delete sendProcEqns_;

   delete essBCEqns_;
}

//------------------------------------------------------------------------------
EqnCommMgr& EqnCommMgr::operator=(const EqnCommMgr& src)
{
   delete recvProcEqns_;
   recvProcEqns_ = src.recvProcEqns_->deepCopy();

   exchangeIndicesCalled_ = src.exchangeIndicesCalled_;
   accumulate_ = src.accumulate_;

   delete recvEqns_;
   recvEqns_ = src.recvEqns_->deepCopy();

   int len = src.solnValues_.length();

   solnValues_.resize(len);

   for(int i=0; i<len; i++) {
      solnValues_[i] = src.solnValues_[i];
   }

   delete sendProcEqns_;
   sendProcEqns_ = src.sendProcEqns_->deepCopy();

   delete sendEqns_;
   sendEqns_ = src.sendEqns_->deepCopy();

   len = src.sendEqnSoln_.length();
   sendEqnSoln_.resize(len);

   for(int i=0; i<len; i++) {
      sendEqnSoln_[i] = src.sendEqnSoln_[i];
   }

   delete essBCEqns_;
   essBCEqns_ = src.essBCEqns_->deepCopy();

   return(*this);
}

//------------------------------------------------------------------------------
EqnCommMgr* EqnCommMgr::deepCopy()
{
//
//This function is like a copy-constructor. This function produces not just a
//structural copy of the current object, it copies EVERYTHING, including data
//contents.
//
   EqnCommMgr* dest = new EqnCommMgr(comm_);
   *dest = *this;
   return(dest);
}

//------------------------------------------------------------------------------
void EqnCommMgr::addLocalEqn(int eqnNumber, int srcProc) {
//
//This function adds the eqnNumber, srcProc pair to the recvProcsEqns_
//object, which does the right thing as far as only putting it in lists that
//it isn't already in, preserving order, etc. 
//
   if (srcProc == localProc_) {
      FEI_CERR << "EqnCommMgr::addRecvEqn: ERROR, srcProc == localProc_, "
          << "which is a recipe for a deadlock." << FEI_ENDL;
      std::abort();
   }

   recvProcEqns_->addEqn(eqnNumber, srcProc);
}

//------------------------------------------------------------------------------
void EqnCommMgr::addSolnValues(int* eqnNumbers, double* values, int num)
{
  if (!exchangeIndicesCalled_) {
    FEI_CERR << "EqnCommMgr::addSolnValues: ERROR, you may not call this until"
      " after exchangeIndices has been called." << FEI_ENDL;
    std::abort();
  }

  feiArray<int>& recvEqnNumbers = recvEqns_->eqnNumbersPtr();

  double* solnValuesPtr = solnValues_.dataPtr();
  for(int i=0; i<num; i++) {
    int index = snl_fei::binarySearch(eqnNumbers[i], recvEqnNumbers);

    if (index < 0) continue;

    solnValuesPtr[index] = values[i];
  }
}

//------------------------------------------------------------------------------
int EqnCommMgr::exchangeIndices(FEI_OSTREAM* dbgOut) {
//
//This function performs the communication necessary to exchange remote
//contributions to the matrix structure (indices only, not coefficients)
//among all participating processors.
//
//Most of this function is #ifdef'd according to whether FEI_SER is
//defined...
//
#ifndef FEI_SER

  int numSendEqns = sendEqns_->getNumEqns();
  SSVec** sendEqnsPtr = sendEqns_->eqns().dataPtr();
  feiArray<int>& sendEqnNumbers = sendEqns_->eqnNumbersPtr();
  feiArray<int> sendEqnLengths(numSendEqns);
  for(int i=0; i<numSendEqns; ++i) {
    sendEqnLengths[i] = sendEqnsPtr[i]->length();
  }

  sendProcEqns_->setProcEqnLengths(sendEqnNumbers.dataPtr(),
                                   sendEqnLengths.dataPtr(),
                                   numSendEqns);

  CHK_ERR( mirrorProcEqns(*sendProcEqns_, *recvProcEqns_) );
  CHK_ERR( mirrorProcEqnLengths(*sendProcEqns_, *recvProcEqns_) );

  sendEqns_->setNumRHSs(sendEqns_->getNumRHSs());
  recvEqns_->setNumRHSs(sendEqns_->getNumRHSs());

  size_t numRecvProcs = recvProcEqns_->getNumProcs();
  size_t numSendProcs = sendProcEqns_->getNumProcs();
  std::vector<int>& sendProcs       = sendProcEqns_->procsPtr();

  //while we're here, let's allocate the array into which we will (later)
  //recv the soln values for the remote equations we've contributed to.
  sendEqnSoln_.resize(numSendEqns);

  //First we'll figure out the expected receive-lengths and the expected
  //send lengths.

  std::vector<int> recvProcTotalLengths(numRecvProcs);
  std::vector<int> sendProcTotalLengths(numSendProcs);

  std::vector<int>& recvProcs = recvProcEqns_->procsPtr();
  std::vector<int>& eqnsPerRecvProc = recvProcEqns_->eqnsPerProcPtr();
  std::vector<std::vector<int>*>& recvProcEqnLengths =
    recvProcEqns_->procEqnLengthsPtr();
  std::vector<std::vector<int>*>& recvProcEqnNumbers = 
    recvProcEqns_->procEqnNumbersPtr();

  std::vector<int>** recvProcEqnLengthsPtr = &recvProcEqnLengths[0];
  std::vector<int>** recvProcEqnNumbersPtr = &recvProcEqnNumbers[0];

  for(unsigned i=0; i<numRecvProcs; i++) {

    //first we need the total of eqn-lengths for this recv-proc
    int totalLength = 0;
    for(int j=0; j<eqnsPerRecvProc[i]; j++) {
      totalLength += (*(recvProcEqnLengthsPtr[i]))[j];
    }
    recvProcTotalLengths[i] = totalLength;
  }

  std::vector<int>& eqnsPerSendProc = sendProcEqns_->eqnsPerProcPtr();
  std::vector<std::vector<int>*>& sendProcEqnNumbers = 
    sendProcEqns_->procEqnNumbersPtr();
  std::vector<std::vector<int>*>& sendProcLengths =
    sendProcEqns_->procEqnLengthsPtr();

  std::vector<int>** sendProcEqnNumbersPtr = &sendProcEqnNumbers[0];
  std::vector<int>** sendProcLengthsPtr = &sendProcLengths[0];

  for(unsigned i=0; i<numSendProcs; i++) {
    int totalLength = 0;
    for(int j=0; j<eqnsPerSendProc[i]; j++) {
      totalLength += (*(sendProcLengthsPtr[i]))[j];
    }
    sendProcTotalLengths[i] = totalLength;
  }

  //Before we get too carried away here, lets make sure that the messages we
  //expect to receive from each other processor are the same length as the
  //other processor plans to send. (Many times I've had crashes in MPI_Wait
  //below, which are always mysterious to figure out. They usually result
  //from mis-matched send/recv lengths.)

  CHK_ERR( consistencyCheck("exchangeIndices", recvProcs, recvProcTotalLengths,
                            sendProcs, sendProcTotalLengths) );

  int** recvProcEqnIndices = NULL;

  if (numRecvProcs > 0) {
    recvProcEqnIndices = new int*[numRecvProcs];
  }

  MPI_Request* indRequests = NULL;

  if (numRecvProcs > 0) {
    indRequests = new MPI_Request[numRecvProcs];
  }

  int numRecvsStarted = 0;
  int indTag = 199904;

  //now, let's start the recvs for the incoming indices.
  for(unsigned i=0; i<numRecvProcs; i++) {

    int totalLength = recvProcTotalLengths[i];

    recvProcEqnIndices[i] = new int[totalLength];

    //launch the recvs for the indices now.
    if (MPI_Irecv(recvProcEqnIndices[i], totalLength, MPI_INT,
                  recvProcs[i], indTag, comm_, &indRequests[i]) != MPI_SUCCESS) {
       ERReturn(-1);
    }

    numRecvsStarted++;
  }

  //now we need to build the lists of outgoing indices and send those.
  std::vector<int> indices;

  for(unsigned i=0; i<numSendProcs; i++) {
    int totalLength = sendProcTotalLengths[i];
    int j;

    indices.resize(totalLength);
    int* indicesPtr = &indices[0];

    int offset = 0;

    for(j=0; j<eqnsPerSendProc[i]; j++) {
      int eqnLoc = sendEqns_->getEqnIndex((*(sendProcEqnNumbersPtr[i]))[j]);
      feiArray<int>& sendIndices = sendEqns_->eqns()[eqnLoc]->indices();
      int* sendIndicesPtr = sendIndices.dataPtr();

      for(int k=0; k<(*(sendProcLengthsPtr[i]))[j]; k++) {
        indicesPtr[offset++] = sendIndicesPtr[k];
      }
    }

    if (MPI_Send(indicesPtr, totalLength, MPI_INT, sendProcs[i],
                 indTag, comm_) != MPI_SUCCESS) ERReturn(-1)
  }

  //and finally, we're ready to complete the irecvs for the indices and put
  //them away.
  int numCompleted = 0;
  for(unsigned i=0; i<numRecvProcs; i++) {
    MPI_Status status;
    int index = i;
    MPI_Wait(&indRequests[i], &status);
    numCompleted++;

    int offset = 0;
    for(int j=0; j<eqnsPerRecvProc[index]; j++) {
      int eqn = (*(recvProcEqnNumbersPtr[index]))[j];
      int* indxs = &(recvProcEqnIndices[index][offset]);
      int len = (*(recvProcEqnLengthsPtr[index]))[j];

      recvEqns_->addIndices(eqn, indxs, len);

      offset += len;
    }

    delete [] recvProcEqnIndices[index];
  }

  delete [] recvProcEqnIndices;
  delete [] indRequests;

  if (numRecvsStarted != numCompleted) {
    FEI_CERR << "EqnCommMgr::exchangeIndices: recv-send mismatch; "
          << "numRecvsStarted: " << numRecvsStarted << ", numCompleted: "
         << numCompleted << FEI_ENDL;
    std::abort();
  }

  //allocate the solnValue_ list, which is of size numRecvEqns.
  int numRecvEqns = recvEqns_->getNumEqns();
  solnValues_.resize(numRecvEqns);

  exchangeIndicesCalled_ = true;

  if (dbgOut != NULL) {
    FEI_OSTREAM& os = *dbgOut;
    os << "#ereb exchangeIndices, sendEqns_:"<<FEI_ENDL;
//    os << *sendEqns_<<FEI_ENDL;
  }

#else
  (void)dbgOut;
#endif // #ifndef FEI_SER

   return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::consistencyCheck(const char* caller,
				 std::vector<int>& recvProcs,
				 std::vector<int>& recvProcTotalLengths,
				 std::vector<int>& sendProcs,
				 std::vector<int>& sendProcTotalLengths)
{
  int err = 0;
  //First we'll gather each processors send-lengths onto all other processors.
  std::vector<int> globalProcSendLengths, globalSendProcs;
  std::vector<int> gatherSizes;

  CHK_ERR( fei::Allgatherv(comm_, sendProcTotalLengths,
				 gatherSizes, globalProcSendLengths) );

  CHK_ERR( fei::Allgatherv(comm_, sendProcs,
				 gatherSizes, globalSendProcs) );

  //Now check the consistency of the global send-lengths against local
  //receive-lengths.
  int offset = 0;
  for(unsigned i=0; i<gatherSizes.size(); i++) {
    int size = gatherSizes[i];
    if ((int)i==localProc_) {offset += size; continue; }

    //now we're ready to stride through processor i's sendProcs. Only do this
    //if processor i is one of our recvProcs.
    std::vector<int>::iterator rp_iter = 
      std::lower_bound(recvProcs.begin(), recvProcs.end(), (int)i);
    int rpIndex = -1;
    if (rp_iter != recvProcs.end() && (int)i == *rp_iter) {
      rpIndex = (int)(rp_iter - recvProcs.begin());
    }

    if (rpIndex < 0) {
      // proc i is not one of our recvProcs. Let's make sure
      //that we're not one of proc i's sendProcs.
      for(int j=0; j<size; j++) {
	if (globalSendProcs[offset+j] == localProc_) {
	  FEI_CERR << "EqnCommMgr::"<<caller<<" ERROR: proc " << localProc_
	       << " is not expecting to receive from proc " << i << " but proc "
	       << i << " is expecting to send to proc " << localProc_ << FEI_ENDL;
	  err = -1;
	}
      }
      if (err == -1) break;

      //if err != -1, simply jump to the next for(i... iteration.
      offset += size;
      continue;
    }

    for(int j=0; j<size; j++) {
      if (globalSendProcs[offset+j] == localProc_) {
	int sendLength = globalProcSendLengths[offset+j];
	int recvLength = recvProcTotalLengths[rpIndex];
	if (sendLength != recvLength) {
	  FEI_CERR << "EqnCommMgr::"<<caller<<" ERROR: proc " << localProc_
	       << " is expecting to receive " << recvLength << " indices from "
	       << "proc " << i << " but proc " << i << " is expecting to send "
	       << sendLength << " indices to proc " << localProc_ << FEI_ENDL;
	  err = -1;
	}
      }
    }

    offset += size;
  }

  int globalErr = 0;
  CHK_ERR( fei::GlobalSum(comm_, err, globalErr) );

  return(globalErr);
}

//------------------------------------------------------------------------------
int EqnCommMgr::exchangeEqns(FEI_OSTREAM* dbgOut)
{
  //
  //This function performs the communication necessary to exchange remote
  //equations (both indices and coefficients) among all participating processors.
  //

  recvEqns_->resetCoefs();

  if (dbgOut != NULL) {
    FEI_OSTREAM& os = *dbgOut;
    os << "#ereb exchangeEqns begin, sendEqns_:"<<FEI_ENDL;
//    os << *sendEqns_<<FEI_ENDL;
  }

  CHK_ERR( exchangeEqnBuffers(comm_, sendProcEqns_,
			      sendEqns_, recvProcEqns_, recvEqns_, accumulate_));

  if (dbgOut != NULL) {
    FEI_OSTREAM& os = *dbgOut;
    os << "#ereb exchangeEqns end, sendEqns_:"<<FEI_ENDL;
//    os << *sendEqns_<<FEI_ENDL;
  }

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::exchangeEqnBuffers(MPI_Comm comm, ProcEqns* sendProcEqns,
                              EqnBuffer* sendEqns, ProcEqns* recvProcEqns,
                              EqnBuffer* recvEqns, bool accumulate)
{
//
//This function performs the communication necessary to exchange remote
//equations (both indices and coefficients) among all participating processors.
//
//Most of this function is #ifdef'd according to whether FEI_SER is
//defined.
#ifdef FEI_SER
  (void)comm;
  (void)sendProcEqns;
  (void)sendEqns;
  (void)recvProcEqns;
  (void)recvEqns;
  (void)accumulate;
#else
   int indTag = 19991130, coefTag = 19991131;

   size_t numRecvProcs = recvProcEqns->getNumProcs();
   size_t numSendProcs = sendProcEqns->getNumProcs();
   if ((numRecvProcs == 0) && (numSendProcs == 0)) return(0);

   //we assume that sendProcEqns and recvProcEqns are fully populated with
   //equation-numbers and their lengths, and that this data is consistent with
   //the data in sendEqns...

   MPI_Request* indRequests = NULL;
   MPI_Request* coefRequests = NULL;
   int** recvProcEqnIndices = NULL;
   double** recvProcEqnCoefs = NULL;

   if (numRecvProcs > 0) {
      indRequests = new MPI_Request[numRecvProcs];
      coefRequests = new MPI_Request[numRecvProcs];
      recvProcEqnIndices = new int*[numRecvProcs];
      recvProcEqnCoefs = new double*[numRecvProcs];
   }

   int numRHSs = sendEqns->getNumRHSs();

   //now, let's allocate the space for the incoming equations.
   //each row of recvProcEqnIndices will be of length
   //sum-of-recvProcEqnLengths[i], and each row of recvProcEqnCoefs will be
   //of length sum-of-recvProcEqnLengths[i] + numRHSs*eqnsPerRecvProc[i].

   std::vector<int>& recvProcs = recvProcEqns->procsPtr();
   std::vector<int>& eqnsPerRecvProc = recvProcEqns->eqnsPerProcPtr();
   std::vector<std::vector<int>*>& recvProcEqnNumbers =
     recvProcEqns->procEqnNumbersPtr();
   std::vector<std::vector<int>*>& recvProcEqnLengths =
     recvProcEqns->procEqnLengthsPtr();

   int padding = 2;
   for(unsigned i=0; i<numRecvProcs; i++) {
      int totalLength = 0;

      for(int j=0; j<eqnsPerRecvProc[i]; j++) {
         totalLength += (*(recvProcEqnLengths[i]))[j];
      }

      //in case we're only exchanging rhs coefs, (i.e., there are no
      //column-indices) let's make the length totalLength+padding so we can
      //send the new*Data_ indicators.
      recvProcEqnIndices[i] = new int[totalLength+padding];

      int coefLength = totalLength + numRHSs*eqnsPerRecvProc[i];

      recvProcEqnCoefs[i] = new double[coefLength];
      //let's go ahead and launch the recvs for the indices and coefs now.
      MPI_Irecv(recvProcEqnIndices[i], totalLength+padding, MPI_INT,
                recvProcs[i], indTag, comm, &indRequests[i]);

      MPI_Irecv(recvProcEqnCoefs[i], coefLength, MPI_DOUBLE,
                recvProcs[i], coefTag, comm, &coefRequests[i]);
   }

   //ok, now we need to build the lists of outgoing indices and coefs, and
   //send those.
   std::vector<int>& sendProcs = sendProcEqns->procsPtr();
   std::vector<int>& eqnsPerSendProc = sendProcEqns->eqnsPerProcPtr();
   std::vector<std::vector<int>*>& sendProcEqnNumbers =
     sendProcEqns->procEqnNumbersPtr();
   std::vector<std::vector<int>*>& sendProcEqnLengths =
     sendProcEqns->procEqnLengthsPtr();

   feiArray<feiArray<double>*>& sendRHS = *(sendEqns->rhsCoefsPtr());

   for(unsigned i=0; i<numSendProcs; i++) {
      int totalLength = 0;
      int* sendProcEqnNumbers_i = &(*(sendProcEqnNumbers[i]))[0];
      int* sendProcEqnLengths_i = &(*(sendProcEqnLengths[i]))[0];
      int j;
      for(j=0; j<eqnsPerSendProc[i]; j++)
         totalLength += sendProcEqnLengths_i[j];

      //As with the recv code above, let's make the indices length
      //be totalLength+padding...
      std::vector<int> indices(totalLength+padding);
      int* indicesPtr = &indices[0];
      int coefLength = totalLength + numRHSs*eqnsPerSendProc[i];

      std::vector<double> coefs(coefLength);
      double* coefsPtr = &coefs[0];

      int offset = 0;

      //first pack up the coefs and indices
      for(j=0; j<eqnsPerSendProc[i]; j++) {
         int eqnLoc = sendEqns->getEqnIndex(sendProcEqnNumbers_i[j]);
	 int* sendIndices = &(sendEqns->eqns()[eqnLoc]->indices())[0];
	 double* sendCoefs= &(sendEqns->eqns()[eqnLoc]->coefs())[0];

         for(int k=0; k<sendProcEqnLengths_i[j]; k++) {
            indicesPtr[offset] = sendIndices[k];
            coefsPtr[offset++] = sendCoefs[k];
         }
      }

      //now append the new*Data_ indicators to the end of the indices array
      indicesPtr[offset] = sendEqns->newCoefData_;
      indicesPtr[offset+1] = sendEqns->newRHSData_;

      //now append the RHS coefs to the end of the coefs array
      for(j=0; j<eqnsPerSendProc[i]; j++) {
         int eqnLoc = sendEqns->getEqnIndex(sendProcEqnNumbers_i[j]);

         for(int k=0; k<numRHSs; k++) {
            coefsPtr[offset++] = (*(sendRHS[eqnLoc]))[k];
         }
      }

      MPI_Send(&indices[0], (int)indices.size(), MPI_INT, sendProcs[i],
               indTag, comm);
      MPI_Send(&coefs[0], coefLength, MPI_DOUBLE, sendProcs[i],
	       coefTag, comm);
   }

   //and finally, we're ready to complete the irecvs for the indices and coefs,
   //and put them away.
   for(unsigned i=0; i<numRecvProcs; i++) {
      MPI_Status status;
      int index = i;
      MPI_Wait(&indRequests[index], &status);
      MPI_Wait(&coefRequests[index], &status);

      int j, offset = 0;
      for(j=0; j<eqnsPerRecvProc[index]; j++) {
         int eqn = (*(recvProcEqnNumbers[index]))[j];
         int* indices = &(recvProcEqnIndices[index][offset]);
         double* coefs = &(recvProcEqnCoefs[index][offset]);
         int len = (*(recvProcEqnLengths[index]))[j];

         recvEqns->addEqn(eqn, coefs, indices, len, accumulate);

         offset += len;
      }

      recvEqns->newCoefData_ += recvProcEqnIndices[index][offset];
      recvEqns->newRHSData_  += recvProcEqnIndices[index][offset+1];
      delete [] recvProcEqnIndices[index];
   }

   //now unpack the RHS entries

   recvEqns->setNumRHSs(numRHSs);

   for(unsigned i=0; i<numRecvProcs; i++) {
      int j, offset = 0;
      for(j=0; j<eqnsPerRecvProc[i]; j++) {
         offset += (*(recvProcEqnLengths[i]))[j];
      }

      for(j=0; j<eqnsPerRecvProc[i]; j++) {
         int eqn = (*(recvProcEqnNumbers[i]))[j];

         for(int k=0; k<numRHSs; k++) {
            CHK_ERR( recvEqns->addRHS(eqn, k, recvProcEqnCoefs[i][offset++],
				      accumulate));
         }
      }

      delete [] recvProcEqnCoefs[i];
   }

   delete [] recvProcEqnIndices;
   delete [] recvProcEqnCoefs;
   delete [] indRequests;
   delete [] coefRequests;

#endif //#ifndef FEI_SER

   return(0);
}

//------------------------------------------------------------------------------
void EqnCommMgr::exchangeSoln()
{
  //Most of this function is #ifdef'd according to whether FEI_SER 
  //is defined...
#ifndef FEI_SER

   int solnTag = 199906;

   MPI_Request* solnRequests = NULL;
   double** solnBuffer = NULL;

   size_t numSendProcs = sendProcEqns_->getNumProcs();
   std::vector<int>& sendProcs = sendProcEqns_->procsPtr();
   std::vector<int>& eqnsPerSendProc = sendProcEqns_->eqnsPerProcPtr();

   if (numSendProcs > 0) {
      solnRequests = new MPI_Request[numSendProcs];
      solnBuffer = new double*[numSendProcs];
   }

   MPI_Comm comm = comm_;

   //let's launch the recv's for the incoming solution values.
   for(unsigned i=0; i<numSendProcs; i++) {
      solnBuffer[i] = new double[eqnsPerSendProc[i]];

      MPI_Irecv(solnBuffer[i], eqnsPerSendProc[i], MPI_DOUBLE, sendProcs[i],
                solnTag, comm, &solnRequests[i]);
   }

   size_t numRecvProcs = recvProcEqns_->getNumProcs();
   std::vector<int>& recvProcs = recvProcEqns_->procsPtr();
   std::vector<int>& eqnsPerRecvProc = recvProcEqns_->eqnsPerProcPtr();
   std::vector<std::vector<int>*>& recvProcEqnNumbers =
     recvProcEqns_->procEqnNumbersPtr();

   //now let's send the outgoing solutions.
   for(unsigned i=0; i<numRecvProcs; i++) {
      double* solnBuff = new double[eqnsPerRecvProc[i]];

      for(int j=0; j<eqnsPerRecvProc[i]; j++) {
	int eqnNumber = (*(recvProcEqnNumbers[i]))[j];

	int index = recvEqns_->getEqnIndex(eqnNumber);
	solnBuff[j] = solnValues_[index];
      }

      MPI_Send(solnBuff, eqnsPerRecvProc[i], MPI_DOUBLE, recvProcs[i],
               solnTag, comm);

      delete [] solnBuff;
   }

   std::vector<std::vector<int>*>& sendProcEqnNumbers =
     sendProcEqns_->procEqnNumbersPtr();

   //make sure the sendEqnSoln_ array is big enough...
   sendEqnSoln_.resize(sendEqns_->getNumEqns());

   //ok, complete the above recvs and store the soln values.
   for(unsigned i=0; i<numSendProcs; i++) {
      int index;
      MPI_Status status;
      MPI_Waitany((int)numSendProcs, solnRequests, &index, &status);

      for(int j=0; j<eqnsPerSendProc[index]; j++) {
	int eqnNumber = (*(sendProcEqnNumbers[index]))[j];
	int ind = sendEqns_->getEqnIndex(eqnNumber);

	sendEqnSoln_[ind] = solnBuffer[index][j];
      }

      delete [] solnBuffer[index];
   }

   delete [] solnRequests;
   delete [] solnBuffer;
#endif //#ifndef FEI_SER
}

//------------------------------------------------------------------------------
int EqnCommMgr::mirrorProcEqns(ProcEqns& inProcEqns, ProcEqns& outProcEqns)
{
  //Beginning assumption: we (the local processor) have a populated ProcEqns
  //object (inProcEqns) which contains information pairing certain equations
  //with certain remote processors. These can be equations that we will be
  //receiving in an all-to-all exchange, or equations that we will be sending in
  //an all-to-all exchange. In either case, the "mirror" of that information is
  //needed before the exchange can be performed. i.e., if we know which eqns
  //we'll be sending, then the mirror info concerns the eqns we'll be recv'ing,
  //and vice-versa if we already know which eqns we'll be recv'ing.
  //
  //This function is to obtain that mirror info, and return it in outProcEqns.
  //
  //Given a populated ProcEqns object, we want to populate the mirror ProcEqns
  //object.
  //
  //First figure out how many procs belong in outProcEqns.
  //Then receive, from each of those procs, the list of associated equations.
  //Then we'll have the info necessary to populate the 'outProcEqns' object.

  //Most of this function is #ifdef'd according to whether FEI_SER 
  //is defined.
#ifdef FEI_SER
  (void)inProcEqns;
  (void)outProcEqns;
#else
  int numProcs = fei::numProcs(comm_);

  if (numProcs < 2) return(0);

  std::vector<int> buf(numProcs*2, 0);

  std::vector<int>& inProcs = inProcEqns.procsPtr();
  std::vector<int> outProcs;
  fei::mirrorProcs(comm_, inProcs, outProcs);

  std::vector<int>& eqnsPerInProc = inProcEqns.eqnsPerProcPtr();

  size_t numOutProcs = outProcs.size();

  std::vector<int> recvbuf(numOutProcs, 0);

  //now send a length (the contents of buf[i]) to each "in-proc"
  //(length of the list of equation data that is to follow).
  MPI_Request* requests = new MPI_Request[numOutProcs];
  MPI_Status* statuses = new MPI_Status[numOutProcs];
  int firsttag = 20051014;
  int offset = 0;
  int* outProcsPtr = &outProcs[0];
  for(unsigned i=0; i<numOutProcs; ++i) {
    if (MPI_Irecv(&(recvbuf[i]), 1, MPI_INT, outProcsPtr[i], firsttag,
		  comm_, &requests[offset++]) != MPI_SUCCESS) ERReturn(-1);
  }

  int* inProcsPtr = &inProcs[0];
  size_t numInProcs = inProcs.size();
  for(unsigned i=0; i<numInProcs; ++i) {
    if (MPI_Send(&(eqnsPerInProc[i]), 1, MPI_INT, inProcsPtr[i], firsttag,
		 comm_) != MPI_SUCCESS) ERReturn(-1);
  }

  MPI_Waitall((int)numOutProcs, requests, statuses);

  delete [] requests;
  delete [] statuses;

  std::vector<int> lengths(numOutProcs);

  offset = 0;
  for(unsigned i=0; i<numOutProcs; ++i) {
    if (recvbuf[i] > 0) {
      lengths[offset++] = recvbuf[i];
    }
  }

  //now we need to create 'numOutProcs' lists, into which we'll receive the
  //equation-numbers that those procs send to us.
  std::vector<std::vector<int>*>* outEqns = NULL;
  if (numOutProcs > 0) {
    outEqns = new std::vector<std::vector<int>*>(numOutProcs);
  }

  for(unsigned i=0; i<numOutProcs; i++) {
    (*outEqns)[i] = new std::vector<int>(lengths[i]);
  }

  //finally we're ready to exchange lists of equations.

  CHK_ERR( fei::exchangeData(comm_, inProcs,
				  inProcEqns.procEqnNumbersPtr(),
				  outProcs,
				  true,
				  *outEqns) );

  //now we've completed all the communication, so we're ready to put the data
  //we received into the outProcEqns object.
  for(unsigned i=0; i<numOutProcs; i++) {
    std::vector<int>* eqnArray = (*outEqns)[i];
    int* eqns = &(*eqnArray)[0];
    size_t len = eqnArray->size();
    for(unsigned j=0; j<len; j++) {
      outProcEqns.addEqn(eqns[j], outProcs[i]);
    }
    delete eqnArray;
  }

  delete outEqns;

#endif //#ifndef FEI_SER

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::mirrorProcEqnLengths(ProcEqns& inProcEqns,
				     ProcEqns& outProcEqns)
{
  //Support function to set up information needed for exchanging equation
  //info among processors.
  //This function plays a similar role to that of the above 'mirrorProcEqns'
  //function, but exchanges the length information rather than the
  //equation-number information. THIS FUNCTION ASSUMES that both inProcEqns and
  //outProcEqns are already populated with eqn-number information.

  //Most of this function is #ifdef'd according to whether FEI_SER 
  //is defined.
#ifdef FEI_SER
  (void)inProcEqns;
  (void)outProcEqns;
#else
  if (fei::numProcs(comm_) == 1) return(0);

  CHK_ERR( fei::exchangeData( comm_, inProcEqns.procsPtr(),
				    inProcEqns.procEqnLengthsPtr(),
				    outProcEqns.procsPtr(),
				    true,
				    outProcEqns.procEqnLengthsPtr() ) );
#endif //#ifndef FEI_SER

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::addRemoteEqn(int eqnNumber, int destProc,
                            const double* coefs, const int* indices, int num) {
   (void)destProc;
   sendEqns_->newCoefData_ = 1;

   return(sendEqns_->addEqn(eqnNumber, coefs, indices, num, accumulate_));
}

//------------------------------------------------------------------------------
int EqnCommMgr::addRemoteEqn(int eqnNumber, const double* coefs,
			     const int* indices, int num)
{
   sendEqns_->newCoefData_ = 1;

   return(sendEqns_->addEqn(eqnNumber, coefs, indices, num, accumulate_));
}

//------------------------------------------------------------------------------
void EqnCommMgr::setNumRHSs(int numRHSs) {
   sendEqns_->setNumRHSs(numRHSs);
   recvEqns_->setNumRHSs(numRHSs);
}

//------------------------------------------------------------------------------
int EqnCommMgr::addRemoteRHS(int eqnNumber, int destProc, int rhsIndex,
                            double value)
{
   (void)destProc;
   sendEqns_->newRHSData_ = 1;
   return(sendEqns_->addRHS(eqnNumber, rhsIndex, value));
}

//------------------------------------------------------------------------------
int EqnCommMgr::addRemoteRHS(int eqnNumber, int rhsIndex, double value)
{
   sendEqns_->newRHSData_ = 1;
   return(sendEqns_->addRHS(eqnNumber, rhsIndex, value));
}

//------------------------------------------------------------------------------
void EqnCommMgr::addRemoteIndices(int eqnNumber, int destProc,
                                int* indices, int num)
{
  if (destProc < 0) {
    FEI_CERR << "fei: EqnCommMgr::addRemoteIndices ERROR, destProc < 0" << FEI_ENDL;
    std::abort();
  }

  sendEqns_->addIndices(eqnNumber, indices, num);

  sendProcEqns_->addEqn(eqnNumber, destProc);
}

//------------------------------------------------------------------------------
void EqnCommMgr::resetCoefs() {
   recvEqns_->resetCoefs();
   sendEqns_->resetCoefs();
   essBCEqns_->resetCoefs();

   sendEqns_->newCoefData_ = 0;
   sendEqns_->newRHSData_ = 0;
   recvEqns_->newCoefData_ = 0;
   recvEqns_->newRHSData_ = 0;

   int numRecvEqns = recvEqns_->getNumEqns();
   for(int i=0; i<numRecvEqns; i++) {
      solnValues_[i] = 0.0;
   }
}

//----------------------------------------------------------------------------
int EqnCommMgr::gatherSharedBCs(EqnBuffer& bcEqns)
{
  //Gather boundary-condition equations from all sharing processors to the
  //owning processors.
  //
  //On entry to this function, bcEqns contains all boundary-condition equations
  //that were specified by the finite-element application on this processor.
  //On exit, bcEqns will also include any boundary-condition equations that were
  //specified by the finite-element application on processors that share nodes
  //that are owned by this processor.
  //
  //We'll set up two new equation buffers: sendBCs and recvBCs. From the input
  //bcEqns, we'll put each eqn that is in sendEqns_ into sendBCs because these
  //are the equations that are shared and remotely owned. We'll then do a
  //gather where all processors send their shared-bc eqns to the owner of those
  //equations. The result of this gather will be in recvBCs, which we will then
  //merge into bcEqns.

  if (fei::numProcs(comm_) == 1) return(0);

  int i;
  EqnBuffer sendBCs, recvBCs;
  ProcEqns sendBCProcEqns, recvBCProcEqns;

  //now loop through the equations in bcEqns, checking whether they're in
  //our sendEqns_ object.
  feiArray<int>& bcEqnNumbers = bcEqns.eqnNumbersPtr();
  int numBCeqns = bcEqnNumbers.length();

  feiArray<int>& sendEqnNumbers = sendEqns_->eqnNumbersPtr();
  std::vector<int>& sendProcs = sendProcEqns_->procsPtr();
  std::vector<std::vector<int>*>& sendProcEqnNumbers = 
    sendProcEqns_->procEqnNumbersPtr();
  size_t numSendProcs = sendProcs.size();

  for(i=0; i<numBCeqns; i++) {
    int eqn = bcEqnNumbers[i];

    int index = snl_fei::binarySearch(eqn, sendEqnNumbers);
    if (index<0) continue;

    feiArray<double>& coefs = bcEqns.eqns()[i]->coefs();
    feiArray<int>& indices = bcEqns.eqns()[i]->indices();
    CHK_ERR( sendBCs.addEqn(eqn,
			    coefs.dataPtr(),
			    indices.dataPtr(),
			    indices.length(), false) );

    for(unsigned p=0; p<numSendProcs; p++) {
      if (std::binary_search(sendProcEqnNumbers[p]->begin(),
                             sendProcEqnNumbers[p]->end(), eqn)) {

        sendBCProcEqns.addEqn(eqn, 3, sendProcs[p]);
      }
    }
  }

  //now set up the required mirror info, then perform the exchange among procs.
  CHK_ERR( mirrorProcEqns(sendBCProcEqns, recvBCProcEqns) );
  CHK_ERR( mirrorProcEqnLengths(sendBCProcEqns, recvBCProcEqns) );
  CHK_ERR( exchangeEqnBuffers(comm_, &sendBCProcEqns,
			      &sendBCs, &recvBCProcEqns, &recvBCs, false) );

  //finally, merge the recvBCs into the bcEqns buffer.
  CHK_ERR( bcEqns.addEqns(recvBCs, false) );

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::exchangeRemEssBCs(int* essEqns, int numEssEqns,double* essAlpha,
				  double* essGamma, MPI_Comm comm,
				  FEI_OSTREAM* dbgOut)
{
  delete essBCEqns_;
  essBCEqns_ = new EqnBuffer();

  EqnBuffer* sendEssEqns = new EqnBuffer();
  ProcEqns* essSendProcEqns = new ProcEqns();

  feiArray<SSVec*>& _sendEqns = sendEqns_->eqns();
  SSVec** _sendEqnsPtr = _sendEqns.dataPtr();
  feiArray<int>& _sendEqnNumbers = sendEqns_->eqnNumbersPtr();
  int* _sendEqnNumbersPtr = _sendEqnNumbers.dataPtr();
  int _numSendEqns = sendEqns_->getNumEqns();

  //check to see if any of the essEqns are in the _sendEqns indices.
  //the ones that are, will need to be sent to other processors.

  if (dbgOut != NULL) {
    FEI_OSTREAM& os = *dbgOut;
    os << "#ereb: num-remote-rows: " << _numSendEqns
       << ", numEssEqns: " << numEssEqns << FEI_ENDL;
//    for(int ii=0; ii<numEssEqns; ++ii) {
//      os << "#ereb, essEqns["<<ii<<"]: "<<essEqns[ii]<<FEI_ENDL;
//    }
    os << "#ereb sendEqns_:"<<FEI_ENDL;
    os << *sendEqns_<<FEI_ENDL;
  }

  int i;
  bool accumulate = false;

  feiArray<int> offsets(numEssEqns);
  int* offsetsPtr = offsets.dataPtr();

  for(int j=0; j<_numSendEqns; j++) {

    feiArray<int>& indices = _sendEqnsPtr[j]->indices();

    snl_fei::binarySearch(numEssEqns, essEqns, offsetsPtr,
			  indices.dataPtr(), indices.length());

    int sendEqn_j = _sendEqnNumbersPtr[j];

    int proc = getSendProcNumber(sendEqn_j);

    const int* sendEqnsPtr_j = indices.dataPtr();

    if (dbgOut != NULL) {
      FEI_OSTREAM& os = *dbgOut;
      os << "#ereb sendeqns["<<j<<"].length: "
         <<_sendEqnsPtr[j]->length()<<", numEssEqns: " << numEssEqns<<FEI_ENDL;
    }

    for(i=0; i<numEssEqns; i++) {

      int index = offsetsPtr[i];

      if (index < 0) continue;

      essSendProcEqns->addEqn(sendEqn_j, proc);

      double coef = essGamma[i]/essAlpha[i];
      int essEqns_i = essEqns[i];

      int* essEqns_i_ptr = &essEqns_i;

      sendEssEqns->addEqn(sendEqn_j, &coef,
			  essEqns_i_ptr, 1, accumulate);

      for(int k=0; k<_sendEqnsPtr[j]->length(); k++) {

	int row = sendEqnsPtr_j[k];

	essBCEqns_->addEqn(row, &coef, essEqns_i_ptr, 1, accumulate);
      }
    }
  }

  if (dbgOut != NULL) {
    FEI_OSTREAM& os = *dbgOut;
    os << "#ereb sendEssEqns:"<<FEI_ENDL;
    os << *sendEssEqns<<FEI_ENDL;
  }

  ProcEqns* essRecvProcEqns = new ProcEqns();

  CHK_ERR( mirrorProcEqns(*essSendProcEqns, *essRecvProcEqns) );

  feiArray<int>& eqnNumbers = sendEssEqns->eqnNumbersPtr();
  SSVec** sendEssEqnsPtr = sendEssEqns->eqns().dataPtr();
  feiArray<int> eqnLengths(eqnNumbers.length());
  for(i=0; i<eqnNumbers.length(); ++i) {
    eqnLengths[i] = sendEssEqnsPtr[i]->length();
  }

  essSendProcEqns->setProcEqnLengths(eqnNumbers.dataPtr(),
				     eqnLengths.dataPtr(),
				     eqnNumbers.length());

  CHK_ERR( mirrorProcEqnLengths(*essSendProcEqns, *essRecvProcEqns) );


  CHK_ERR( exchangeEqnBuffers(comm, essSendProcEqns, sendEssEqns,
			      essRecvProcEqns, essBCEqns_, accumulate) );

  if (dbgOut != NULL) {
    FEI_OSTREAM& os = *dbgOut;
    os << "essBCEqns:"<<FEI_ENDL;
    os << *essBCEqns_ << FEI_ENDL;
  }

  delete sendEssEqns;
  delete essSendProcEqns;
  delete essRecvProcEqns;

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::exchangePtToBlkInfo(snl_fei::PointBlockMap& blkEqnMapper)
{
  std::set<int> sendIndices;
  feiArray<SSVec*>& sendeqns = sendEqns_->eqns();
  for(int i=0; i<sendeqns.length(); ++i) {
    feiArray<int>& indices = sendeqns[i]->indices();
    int len = indices.length();
    if (len < 1) continue;
    int* indicesPtr = indices.dataPtr();
    for(int j=0; j<len; ++j) {
      sendIndices.insert(indicesPtr[j]);
    }
  }

  std::set<int> recvIndices;
  feiArray<SSVec*>& recveqns = recvEqns_->eqns();
  for(int i=0; i<recveqns.length(); ++i) {
    feiArray<int>& indices = recveqns[i]->indices();
    int len = indices.length();
    if (len < 1) continue;
    int* indicesPtr = indices.dataPtr();
    for(int j=0; j<len; ++j) {
      recvIndices.insert(indicesPtr[j]);
    }
  }

  std::map<int,int>* ptEqns =  blkEqnMapper.getPtEqns();
  size_t numPtEqns = ptEqns->size();

  std::map<int,int>::const_iterator
    pteq = ptEqns->begin(),
    pteq_end = ptEqns->end();

  std::vector<int> ptBlkInfo(numPtEqns*3);
  int* ptBlkInfoPtr = &ptBlkInfo[0];

  int offset = 0;
  for(; pteq!=pteq_end; ++pteq) {
    int ptEqn = (*pteq).first;
    if (sendIndices.find(ptEqn) == sendIndices.end()) continue;

    int blkEqn = blkEqnMapper.eqnToBlkEqn(ptEqn);
    int blkSize = blkEqnMapper.getBlkEqnSize(blkEqn);

    ptBlkInfoPtr[offset++] = ptEqn;
    ptBlkInfoPtr[offset++] = blkEqn;
    ptBlkInfoPtr[offset++] = blkSize;
  }

  ptBlkInfo.resize(offset);

  size_t numRecvProcs = recvProcEqns_->getNumProcs();
  std::vector<int>& recvProcs = recvProcEqns_->procsPtr();
  size_t numSendProcs = sendProcEqns_->getNumProcs();
  std::vector<int>& sendProcs = sendProcEqns_->procsPtr();

  std::vector<std::vector<int>* > recvData(numRecvProcs);
  for(unsigned i=0; i<numRecvProcs; ++i) {
    recvData[i] = new std::vector<int>;
  }

  std::vector<std::vector<int>* > sendData(numSendProcs);
  for(unsigned i=0; i<numSendProcs; ++i) {
    sendData[i] = &ptBlkInfo;
  }

  CHK_ERR( fei::exchangeData(comm_, sendProcs, sendData,
				    recvProcs, false, recvData) );

  for(unsigned i=0; i<numRecvProcs; ++i) {
    size_t len = recvData[i]->size()/3;
    int* dataPtr = &(*(recvData[i]))[0];

    offset = 0;
    for(unsigned eq=0; eq<len; ++eq) {
      int ptEqn = dataPtr[offset++];

      if (recvIndices.find(ptEqn) == recvIndices.end()) {
        offset += 2; continue;
      }

      int blkEqn = dataPtr[offset++];
      int blkSize = dataPtr[offset++];

      blkEqnMapper.setEqn(ptEqn, blkEqn, blkSize);
    }
    delete recvData[i];
  }

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::addRemoteEqns(SSMat& mat, bool onlyIndices)
{
  int numRows = mat.getRowNumbers().length();
  int* rowNumbers = mat.getRowNumbers().dataPtr();
  feiArray<SSVec*>& rows = mat.getRows();

  for(int i=0; i<numRows; i++) {
    int proc = getSendProcNumber(rowNumbers[i]);
    if (proc == localProc_ || proc < 0) continue;

    if (onlyIndices) {
      addRemoteIndices(rowNumbers[i], proc,
		       rows[i]->indices().dataPtr(),rows[i]->indices().length());
    }
    else {
      CHK_ERR( addRemoteEqn(rowNumbers[i], proc, rows[i]->coefs().dataPtr(),
			    rows[i]->indices().dataPtr(),
			    rows[i]->indices().length()) );
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::addRemoteRHS(SSVec& vec, int rhsIndex)
{
  feiArray<int>& indices = vec.indices();
  feiArray<double>& coefs = vec.coefs();

  for(int i=0; i<indices.length(); i++) {
    int proc = getSendProcNumber(indices[i]);

    if (proc == localProc_ || proc < 0) continue;

    CHK_ERR( addRemoteRHS(indices[i], proc, rhsIndex, coefs[i]) );
  }

  return(0);
}

//------------------------------------------------------------------------------
int EqnCommMgr::getSendProcNumber(int eqn)
{
  size_t numSendProcs = sendProcEqns_->getNumProcs();
  std::vector<int>& sendProcs = sendProcEqns_->procsPtr();
  std::vector<std::vector<int>*>& sendProcEqnNumbers = 
    sendProcEqns_->procEqnNumbersPtr();

  for(unsigned i=0; i<numSendProcs; i++) {
    if (std::binary_search(sendProcEqnNumbers[i]->begin(),
                           sendProcEqnNumbers[i]->end(), eqn)) {

      return(sendProcs[i]);
    }
  }

  return(-1);
}


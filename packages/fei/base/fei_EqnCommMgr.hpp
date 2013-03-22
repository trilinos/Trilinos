#ifndef _fei_EqnCommMgr_hpp_
#define _fei_EqnCommMgr_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_fwd.hpp"
#include "fei_mpi.h"

#include "snl_fei_PointBlockMap.hpp"

#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>
#include <fei_CommUtils.hpp>
#include "fei_ProcEqns.hpp"
#include "fei_EqnBuffer.hpp"

/**
  The EqnCommMgr (Equation communication manager) class is responsible
  for keeping track of equations that require communications. There
  are two types of equations in this class:
  
   1. Local equations which remote processors contribute to (e.g., because
   they share some local active nodes).
  
   2. Remote equations that the local processor contributes to (the mirror of
   case 1.).

  
  Usage Notes:
  
   1. You can't call exchangeEqns until after exchangeIndices has been called.
   2. You can't call addSolnValues until after exchangeIndices has been
      called.
  
  In general, usage will proceed like this:
  
  in snl_fei::Structure::initComplete:
         addLocalEqn  (a bunch of times, probably)
         addRemoteIndices  (also a bunch of times)
  
         exchangeIndices
  
         getNumLocalEqns
         localEqnNumbers
         localIndicesPtr
  
  in Filter::sumInElem/sumInElemMatrix/sumInElemRHS
         addRemoteEqn and/or addRemoteRHS
  
  in Filter::exchangeRemoteEquations
         exchangeEqns
         getNumLocalEqns
         recvEqnNumbersPtr
         localIndicesPtr
         localCoefsPtr
         localRHSsPtr
  
  in Filter::unpackSolution
         getNumLocalEqns
         localEqnNumbers
         addSolnValues
         exchangeSoln
  
  in Filter various getSoln functions
         getNumSendEqns
         sendEqnNumbersPtr
         sendEqnSolnPtr

This class also provides general support functions for exchanging equations
among processors. At the start of an all-to-all exchange, all processors
usually know which equations they need to receive, and which 
processors they'll be receiving from, OR they know which equations they need to
send, and which processors they'll be sending to. Usually they don't know both
the sending and receiving information though.
This eqn-processor pairing information (for either the "recv-equations" or the
"send-equations") is held in a ProcEqns object. Thus, an all-to-all exchange
of equation data requires two ProcEqns objects -- one holding the recv info, the
other holding the send info.

An EqnCommMgr function (mirrorProcEqns) is provided which does the
communications necessary to populate the 'send' ProcEqns object given a
populated 'recv' ProcEqns object, or vice-versa. Note that at this point, the
equation-lengths need not be known by either the sending or the receiving
processors.

Once the ProcEqns objects have been correctly mirrored, another EqnCommMgr
function (mirrorProcEqnLengths) is available for mirroring the eqn-length
data from one ProcEqns object to the other.

The next step is: given equations (with all associated data) that we need to
send, along with previously known equation-numbers that we need to receive, an
EqnCommMgr function (exchangeEqnBuffers) is provided to first send the
equation-lengths to the receiving processors, followed by the actual equation
data. At this point the exchange is complete. The equation data is 
supplied/returned in EqnBuffer objects.
*/

class EqnCommMgr {
 public:
  /** Constructor.
      @param localProc The MPI rank of 'this' processor.
  */
   EqnCommMgr(MPI_Comm comm, bool accumulate = true);

   /** copy constructor */
   EqnCommMgr(const EqnCommMgr& src);

   /** assignment operator */
   EqnCommMgr& operator=(const EqnCommMgr& src);

   /** Destructor. */
   virtual ~EqnCommMgr();

   /** Produce a clone of 'this' object, including all of its internal data.*/
   EqnCommMgr* deepCopy();

   /** return the number of processors that share (contribute to) equations that
       are locally owned.
   */
   size_t getNumSharingProcs() {return(recvProcEqns_->getNumProcs());};
   std::vector<int>& sharingProcsPtr() {return(recvProcEqns_->procsPtr());};

   /** return the number of processors that own equations that we share
       (equations that we contribute to).
   */
   size_t getNumOwnerProcs() {return(sendProcEqns_->getNumProcs());};
   std::vector<int>& ownerProcsPtr() {return(sendProcEqns_->procsPtr());};

   /** add a local equation to which a remote processor will be contributing. */
   void addLocalEqn(int eqnNumber, int srcProc);

   void addSolnValues(int* eqnNumbers, double* values, int num);

#ifdef FEI_HAVE_IOSFWD
   int exchangeIndices(std::ostream* dbgOut=NULL);
   int exchangeEqns(std::ostream* dbgOut=NULL);
#else
   int exchangeIndices(ostream* dbgOut=NULL);
   int exchangeEqns(ostream* dbgOut=NULL);
#endif

   void exchangeSoln();

   /** Support function to set up information needed for exchanging equation
       info among processors.
       Beginning assumption: we (the local processor) have a populated ProcEqns
       object (inProcEqns) which contains information pairing certain equations
       with certain remote processors. These can be equations that we will be
       receiving in an all-to-all exchange, or equations that we will be sending
       in an all-to-all exchange. In either case, the "mirror" of that
       information is needed before the exchange can be performed. i.e., if we
       know which eqns we'll be sending, then the mirror info concerns the eqns
       we'll be recv'ing, and vice-versa if we already know which eqns we'll be
       recv'ing.
       
       This function is to obtain that mirror info, and return it in
       outProcEqns. Given a populated ProcEqns object, we want to populate the
       mirror ProcEqns object. Note that this function IGNORES any equation-
       lengths, if they are present. i.e., the eqn-length info in 'inProcEqns'
       is not referenced, and on completion the eqn-length info in 'outProcEqns'
       is not populated.

       @param inProcEqns Contains the input equation-numbers and associated
       processors.
       @param outProcEqns Output.
   */
   int mirrorProcEqns(ProcEqns& inProcEqns, ProcEqns& outProcEqns);

   /** Support function to set up information needed for exchanging equation
       info among processors. This function plays a similar role to that of the
       above 'mirrorProcEqns' function, but exchanges the length information
       rather than the equation-number information. THIS FUNCTION ASSUMES that
       both inProcEqns and outProcEqns are already populated with eqn-number
       information.
       @param inProcEqns Contains equation-numbers, and their associated
       lengths.
       @param outProcEqns Input/Output. On entry, contains the equation-numbers
       (but not their lengths.
       @return error non-zero if error occurs.
   */
   int mirrorProcEqnLengths(ProcEqns& inProcEqns,
			    ProcEqns& outProcEqns);


   static int exchangeEqnBuffers(MPI_Comm comm, ProcEqns* sendProcEqns,
				 EqnBuffer* sendEqns, ProcEqns* recvProcEqns,
				 EqnBuffer* recvEqns, bool accumulate);

   int getNumLocalEqns() {return(recvEqns_->getNumEqns());};

   std::vector<int>& localEqnNumbers() {return(recvEqns_->eqnNumbers());};
   std::vector<fei::CSVec*>& localEqns(){return(recvEqns_->eqns());};
   std::vector<std::vector<double>*>* localRHSsPtr()
     {return(recvEqns_->rhsCoefsPtr());};

   int addRemoteEqn(int eqnNumber, int destProc, const double* coefs,
                   const int* indices, int num);

   int addRemoteEqn(int eqnNumber, const double* coefs,
		    const int* indices, int num);

   int addRemoteEqns(fei::CSRMat& mat, bool onlyIndices);
   int addRemoteRHS(fei::CSVec& vec, int rhsIndex);

   void setNumRHSs(int numRHSs);

   int addRemoteRHS(int eqnNumber, int destProc, int rhsIndex, double value);

   int addRemoteRHS(int eqnNumber, int rhsIndex, double value);

   void addRemoteIndices(int eqnNumber, int destProc, int* indices, int num);

   int getNumRemoteEqns() {return(sendEqns_->getNumEqns());};

   std::vector<int>& sendEqnNumbersPtr() {return(sendEqns_->eqnNumbers());};

   double* sendEqnSolnPtr() {return(sendEqnSoln_.size()>0? &sendEqnSoln_[0] : NULL);};

   void resetCoefs();

   int gatherSharedBCs(EqnBuffer& bcEqns);

   int exchangeRemEssBCs(int* essEqns, int numEssEqns, double* essAlpha,
			 double* essGamma, MPI_Comm comm,
			 std::ostream* dbgOut = NULL);

   int getNumRemEssBCEqns() {return(essBCEqns_->getNumEqns());};
   std::vector<int>& remEssBCEqnNumbersPtr() {return(essBCEqns_->eqnNumbers());};
   std::vector<fei::CSVec*>& remEssBCEqns() {return(essBCEqns_->eqns());};

   int exchangePtToBlkInfo(snl_fei::PointBlockMap& blkEqnMapper);

   bool newCoefData() {if (recvEqns_->newCoefData_>0) return(true);
                       else return(false);}
   bool newRHSData() {if (recvEqns_->newRHSData_>0) return(true);
                      else return(false);}

   bool accumulate_;

   EqnBuffer* getRecvEqns() { return( recvEqns_ ); }
   EqnBuffer* getSendEqns() { return( sendEqns_ ); }
   ProcEqns* getRecvProcEqns() { return( recvProcEqns_ ); }
   ProcEqns* getSendProcEqns() { return( sendProcEqns_ ); }

 private:
   void deleteEssBCs();
   int getSendProcNumber(int eqn);

   int consistencyCheck(const char* caller,
			std::vector<int>& recvProcs,
			std::vector<int>& recvProcTotalLengths,
			std::vector<int>& sendProcs,
			std::vector<int>& sendProcTotalLengths);

   int localProc_;

   ProcEqns* recvProcEqns_;

   bool exchangeIndicesCalled_; //whether or not the exchangeIndices function
                                //has been called yet.

   EqnBuffer* recvEqns_;

   std::vector<double> solnValues_; //solution values we'll need to return to the
                              //processors that contribute to our equations

   ProcEqns* sendProcEqns_;

   EqnBuffer* sendEqns_;

   std::vector<double> sendEqnSoln_; 
                          //the solution values for the send equations. i.e.,
                          //we'll recv these solution values for the equations
                          //that we contributed to (sent) for other processors.

   EqnBuffer* essBCEqns_;

   MPI_Comm comm_;
};

#endif


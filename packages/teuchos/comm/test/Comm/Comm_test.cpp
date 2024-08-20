// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_as.hpp"


//
// Unit test for Teuchos::Comm
//

template<typename Ordinal>
bool checkSumResult(
  const Teuchos::Comm<Ordinal> &comm,
  const Teuchos::RCP<Teuchos::FancyOStream> &out,
  const bool result
  )
{
  *out << "\nChecking that the above test passed in all processes ...";
  int thisResult = ( result ? 1 : 0 );
  int sumResult = -1;
  reduceAll(comm,Teuchos::REDUCE_SUM,Ordinal(1),&thisResult,&sumResult);
  const bool passed = sumResult==size(comm);
  if(passed)
    *out << " passed\n";
  else
    *out << " (sumResult="<<sumResult<<"!=numProcs) failed\n";
  return passed;
}


template<typename Ordinal, typename Packet>
bool testComm(
  const Teuchos::Comm<Ordinal> &comm,
  const Teuchos::RCP<Teuchos::FancyOStream> &out
  )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::as;

  typedef Teuchos::ScalarTraits<Packet> ST;
  typedef Teuchos::OrdinalTraits<Ordinal> OT;

  OSTab tab(out);

  bool success = true, result;

  *out
    << "\n***"
    << "\n*** testComm<"<<OT::name()<<","<<ST::name()<<">(...)"
    << "\n***\n";

  *out << "\nTesting Comm = " << comm.description() << "\n";

  const int procRank = rank(comm);
  const int numProcs = size(comm);

  *out
    << "\nnumProcs = size(comm) = " << numProcs << "\n"
    << "\nprocRank = rank(comm) = " << procRank << "\n";

  const Ordinal count = numProcs*2;

  Teuchos::Array<Packet> sendBuff(count), recvBuff(count), recvBuff2(count);
  for( int i = 0; i < count; ++i )
    sendBuff[i] = Packet(procRank+1)*Packet(i);

  //
  // send/receive
  //

  if(numProcs > 1) {

#ifdef TEUCHOS_MPI_COMM_DUMP
    Teuchos::MpiComm<Ordinal>::show_dump = true;
#endif

    if(procRank==numProcs-1) {
      *out << "\nSending data from p="<<procRank<<" to the root process (see p=0 output!) ...\n";
      send(comm,count,&sendBuff[0],0);
    }

    if(procRank==0) {
      *out << "\nReceiving data specifically from p="<<numProcs-1<<" ...\n";
      std::fill_n(&recvBuff[0],count,Packet(0));
      const int sourceRank = receive(comm,numProcs-1,count,&recvBuff[0]);
      result = sourceRank ==numProcs-1;
      *out
        << "\nChecking that sourceRank="<<sourceRank<<" == numProcs-1="<<(numProcs-1)
        << " : " << (result ? "passed" : "falied" ) << "\n";
      *out << "\nChecking that recvBuffer[] == numProcs * sendBuffer[] ...";
      result = true;
      for( int i = 0; i < count; ++i ) {
        const Packet expected = Packet(numProcs)*sendBuff[i];
        if( recvBuff[i] != expected ) {
          result = false;
          *out
            << "\n  recvBuffer["<<i<<"]="<<recvBuff[i]
            << " == numProcs*sendBuffer["<<i<<"]="<<expected<<" : failed";
        }
      }
      if(result) {
        *out << " passed\n";
      }
      else {
        *out << "\n";
        success = false;
      }
    }

#ifdef TEUCHOS_MPI_COMM_DUMP
    Teuchos::MpiComm<Ordinal>::show_dump = false;
#endif

  }


  //
  // broadcast/reduceAll(sum)
  //

  if(procRank==0) {
    std::copy(&sendBuff[0],&sendBuff[0]+count,&recvBuff[0]);
    *out << "\nSending broadcast of data from sendBuff[] in root process to recvBuff[] in each process ...\n";
  }
  else {
    std::fill_n(&recvBuff[0],count,Packet(0));
    *out << "\nReceiving broadcast of data from root process into recvBuff[] ...\n";
  }

  broadcast(comm,0,count,&recvBuff[0]);

  *out << "\nSumming broadcasted data recvBuff[] over all processes into recvBuff2[] ...\n";

  reduceAll(comm,Teuchos::REDUCE_SUM,count,&recvBuff[0],&recvBuff2[0]);

  *out << "\nChecking that recvBuff2[i] == numProcs * i ...";
  result = true;
  for( int i = 0; i < count; ++i ) {
    const Packet expected = Packet(numProcs)*Packet(i);
    //*out << "\nexpected["<<i<<"]=numProcs*i="<<Packet(numProcs)<<"*"<<Packet(i)<<"="<<expected<<"\n";
    if( recvBuff2[i] != expected ) {
      result = false;
      *out
        << "\n  recvBuffer2["<<i<<"]="<<recvBuff2[i]
        << " == numProcs*"<<i<<"="<<expected<<" : failed";
    }
  }
  if(result) {
    *out << " passed\n";
  }
  else {
    *out << "\n";
    success = false;
  }

  result = checkSumResult(comm,out,result);
  if(!result) success = false;

  //
  // reduceAll(min)
  //

  if( ST::isComparable ) {

    *out << "\nTaking min of sendBuff[] and putting it in recvBuff[] ...\n";

    reduceAll(comm,Teuchos::REDUCE_MIN,count,&sendBuff[0],&recvBuff[0]);

    *out << "\nChecking that recvBuff[i] == i ...";
    result = true;
    for( int i = 0; i < count; ++i ) {
      const Packet expected = Packet(i);
      //*out << "\nexpected["<<i<<"]=numProcs*i="<<Packet(numProcs)<<"*"<<Packet(i)<<"="<<expected<<"\n";
      if( recvBuff[i] != expected ) {
        result = false;
        *out
          << "\n  recvBuffer["<<i<<"]="<<recvBuff[i]
          << " == "<<i<<"="<<expected<<" : failed";
      }
    }
    if(result) {
      *out << " passed\n";
    }
    else {
      *out << "\n";
      success = false;
    }

    result = checkSumResult(comm,out,result);
    if(!result) success = false;

  }

  //
  // reduceAll(max)
  //

  if( ST::isComparable ) {

    *out << "\nTaking max of sendBuff[] and putting it in recvBuff[] ...\n";

    reduceAll(comm,Teuchos::REDUCE_MAX,count,&sendBuff[0],&recvBuff[0]);

    *out << "\nChecking that recvBuff[i] == numProcs*i ...";
    result = true;
    for( int i = 0; i < count; ++i ) {
      const Packet expected = Packet(numProcs)*Packet(i);
      //*out << "\nexpected["<<i<<"]=numProcs*i="<<Packet(numProcs)<<"*"<<Packet(i)<<"="<<expected<<"\n";
      if( recvBuff[i] != expected ) {
        result = false;
        *out
          << "\n  recvBuffer["<<i<<"]="<<recvBuff[i]
          << " == numProcs*"<<i<<"="<<expected<<" : failed";
      }
    }
    if(result) {
      *out << " passed\n";
    }
    else {
      *out << "\n";
      success = false;
    }

    result = checkSumResult(comm,out,result);
    if(!result) success = false;

  }

  //
  // gatherAll
  //

  *out << "\nGathering all data from sendBuff[] in each process to all processes to allRecvBuff ...\n";

  Teuchos::Array<Packet>
    allRecvBuff(count*numProcs);

  gatherAll(comm,count,&sendBuff[0],Ordinal(allRecvBuff.size()),&allRecvBuff[0]);

  *out << "\nChecking that allRecvBuff[count*k+i] == (k+1) * i ...";
  result = true;
  for( int k = 0; k < numProcs; ++k ) {
    for( int i = 0; i < count; ++i ) {
      const Packet expected = Packet(k+1)*Packet(i);
      if( allRecvBuff[count*k+i] != expected ) {
        result = false;
        *out
          << "\n  allRecvBuff["<<count<<"*"<<k<<"+"<<i<<"]="<<allRecvBuff[count*k+i]
          << " == (k+1)*i="<<expected<<" : failed";
      }
    }
  }
  if(result) {
    *out << " passed\n";
  }
  else {
    *out << "\n";
    success = false;
  }

  result = checkSumResult(comm,out,result);
  if(!result) success = false;

  //
  // scan
  //

  *out << "\nPerforming a scan sum of sendBuff[] into recvBuff[] ...\n";

  std::fill_n(&recvBuff[0],count,Packet(0));

  scan(comm,Teuchos::REDUCE_SUM,count,&sendBuff[0],&recvBuff[0]);

  *out << "\nChecking that recvBuff[i] == sum(k+1,k=0...procRank) * i ...";
  result = true;
  int sumProcRank = 0;
  for( int k = 0; k <= procRank; ++k ) sumProcRank += (k+1);
  for( int i = 0; i < count; ++i ) {
    const Packet expected = Packet(sumProcRank)*Packet(i);
    //*out << "\nexpected["<<i<<"]=sum(k+1,k=0...procRank)*i="<<Packet(sumProcRank)<<"*"<<Packet(i)<<"="<<expected<<"\n";
    if( recvBuff[i] != expected ) {
      result = false;
      *out
        << "\n  recvBuffer["<<i<<"]="<<recvBuff[i]
        << " == sum(k+1,k=0...procRank)*"<<i<<"="<<expected<<" : failed";
    }
  }
  if(result) {
    *out << " passed\n";
  }
  else {
    *out << "\n";
    success = false;
  }

  result = checkSumResult(comm,out,result);
  if(!result) success = false;

  //
  // The End!
  //

  if(success)
    *out << "\nCongratulations, all tests for this Comm check out!\n";
  else
    *out << "\nOh no, at least one of the tests for this Comm failed!\n";

  return success;

}

template<typename Ordinal>
bool masterTestComm(
  const Teuchos::RCP<Teuchos::FancyOStream>    &out
  )
{

  bool success = true, result;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;

  typedef Teuchos::OrdinalTraits<Ordinal> OT;

  OSTab tab(out);

  RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();

#ifdef HAVE_MPI

  // Test that the DefaultComm is really a DefaultMpiComm.
  RCP<const Teuchos::MpiComm<Ordinal> >
    mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<Ordinal> >( comm, false );

  if (mpiComm == Teuchos::null) {
    success = false;
    *out << "\n*** FAILED to cast the Teuchos::DefaultComm<"<< OT::name() << "> to a Teuchos::MpiComm<" << OT::name() << ">!\n";
  }
  else {
    *out
      << "\n***"
      << "\n*** Successfully casted the Teuchos::DefaultComm<"<< OT::name() << "> to a Teuchos::MpiComm<" << OT::name() << ">!"
      << "\n***\n";

    // Now get the raw pointer to the MPI_Comm object
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = mpiComm->getRawMpiComm();

    if (static_cast<MPI_Comm>(*rawMpiComm) == 0) {
      success = false;
      *out << "\n*** FAILED to get the raw MPI_Comm pointer from the Teuchos::MpiComm<" << OT::name() << ">!\n";
    }
    else {
      *out
        << "\n***"
        << "\n*** Successfully got the raw MPI_Comm pointer from the Teuchos::MpiComm<" << OT::name() << ">!"
        << "\n***\n";
    }
  }

#endif

  *out
    << "\n***"
    << "\n*** Created a Comm of type " << comm->description() << " for testing"
    << "\n***\n";

  *out << "\nOrdinal type = "<<OT::name()<<" with an extent of "<<sizeof(Ordinal)<<" bytes\n";

  if( comm->getSize() <= 4 ) {
    result = testComm<Ordinal,char>(*comm,out);
    if(!result) success = false;
  }

  result = testComm<Ordinal,int>(*comm,out);
  if(!result) success = false;

  result = testComm<Ordinal,size_t>(*comm,out);
  if(!result) success = false;

  result = testComm<Ordinal,float>(*comm,out);
  if(!result) success = false;

  result = testComm<Ordinal,double>(*comm,out);
  if(!result) success = false;

#ifdef HAVE_TEUCHOS_COMPLEX

#  ifdef HAVE_TEUCHOS_FLOAT

  result = testComm<Ordinal,std::complex<float> >(*comm,out);
  if(!result) success = false;

#  endif // HAVE_TEUCHOS_FLOAT

  result = testComm<Ordinal,std::complex<double> >(*comm,out);
  if(!result) success = false;

#endif // HAVE_TEUCHOS_COMPLEX

  return success;

}

//
// Main driver program
//

int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::CommandLineProcessor;

  bool success = true, result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    bool   showTimers = true;

    clp.setOption( "show-timers", "no-show-timers", &showTimers, "Determine if timers are shown or not" );

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    RCP<FancyOStream>
      out = VerboseObjectBase::getDefaultOStream();

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl << std::endl;

    result = masterTestComm<short int>(out);
    if(!result) success = false;

    result = masterTestComm<int>(out);
    if(!result) success = false;

    result = masterTestComm<long int>(out);
    if(!result) success = false;

    if(showTimers) {
      Teuchos::TimeMonitor::summarize(
        *out<<"\n"
        ,out->getOutputToRootOnly() < 0 // Show local time or not
        );
    }

    if(success)
      *out << "\nEnd Result: TEST PASSED\n";

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  return ( success ? 0 : 1 );

}

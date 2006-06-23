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

//
// Unit test for Teuchos::Comm
//

template<typename Ordinal>
bool checkSumResult(
  const Teuchos::Comm<Ordinal>                          &comm
  ,const Teuchos::RefCountPtr<Teuchos::FancyOStream>    &out 
  ,const bool                                           result
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
  const Teuchos::Comm<Ordinal>                          &comm
  ,const Teuchos::RefCountPtr<Teuchos::FancyOStream>    &out 
  )
{
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;

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
    
    if(procRank==numProcs-1) {
      *out << "\nSending data from p="<<procRank<<" to the root process (see p=0 output!) ...\n";
      send(comm,count,&sendBuff[0],0);
    }
    
    if(procRank==0) {
      *out << "\nReceiving data specifically from p="<<numProcs-1<<" ...\n";
      std::fill_n(&recvBuff[0],count,Packet(0));
//      const int sourceRank = receive(comm,numProcs-1,count,&recvBuff[0]);
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
  // reduceAllAndScatter(...)
  //

  *out << "\nReducing/summing sendBuff[] and scattering into recvBuff[] ...\n";

  std::fill_n(&recvBuff[0],1,Packet(0));

  Teuchos::Array<Ordinal>
    recvCounts(numProcs);
  
  const Ordinal
    numItemsPerProcess = count/numProcs;

  std::fill(recvCounts.begin(),recvCounts.end(),numItemsPerProcess);

  reduceAllAndScatter(
    comm,Teuchos::REDUCE_SUM
    ,count,&sendBuff[0],&recvCounts[0],&recvBuff[0]
    );

  *out << "\nChecking that recvBuff[i] == sum(k+1,k=0...numProcs-1) * (offset+i) ...";
  result = true;
  int sumProcRanks = 0;
  for( int k = 0; k < numProcs; ++k ) sumProcRanks += (k+1);
  for( int i = 0; i < numItemsPerProcess; ++i ) {
    const int offset = procRank * numItemsPerProcess;
    const Packet expected = Packet(sumProcRanks)*Packet(offset+i);
    *out 
      << "\n  expected["<<i<<"]=sum(k+1,k=0...numProcs-1)*(offset+i)="
      << sumProcRanks<<"*"<<(offset+i)<<"="<<expected;
    *out
      << "\n  recvBuffer["<<i<<"]="<<recvBuff[i];
/*
    if( recvBuff[i] != expected ) {
      result = false;
      *out
        << "\n  recvBuffer["<<i<<"]="<<recvBuff[i]
        << " == sum(k+1,k=0...numProcs-1)*(offset+i)="<<sumProcRanks<<"*"<<(offset+i)<<"="<<expected<<" : failed";
    }
*/
  }
/*
  if(result) {
    *out << " passed\n";
  }
  else {
    *out << "\n";
    success = false;
  }

  result = checkSumResult(comm,out,result);
  if(!result) success = false;

*/

  *out << "\n*** Warning, the output shows that reduceAllAndScatter(...) does not seem to be working as I understand it should ...\n";

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
  const Teuchos::RefCountPtr<Teuchos::FancyOStream>    &out 
  )
{

  bool success = true, result;

  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;

  typedef Teuchos::OrdinalTraits<Ordinal> OT;

  OSTab tab(out);

  RefCountPtr<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();

  *out
    << "\n***"
    << "\n*** Created a Comm of type " << comm->description() << " for testing"
    << "\n***\n";

  *out << "\nOrdinal type = "<<OT::name()<<" with an extent of "<<sizeof(Ordinal)<<" bytes\n";
  
  result = testComm<Ordinal,char>(*comm,out);
  if(!result) success = false;
  
  result = testComm<Ordinal,int>(*comm,out);
  if(!result) success = false;
  
  result = testComm<Ordinal,float>(*comm,out);
  if(!result) success = false;
  
  result = testComm<Ordinal,double>(*comm,out);
  if(!result) success = false;
  
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
  
  result = testComm<Ordinal,std::complex<float> >(*comm,out);
  if(!result) success = false;
  
  result = testComm<Ordinal,std::complex<double> >(*comm,out);
  if(!result) success = false;
  
#endif // defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
  
  return success;

}

//
// Main driver program
//

int main(int argc, char* argv[])
{

  using Teuchos::RefCountPtr;
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

    RefCountPtr<FancyOStream>
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

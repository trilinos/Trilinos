/*
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
*/

#include "Teuchos_UnitTestHarness.hpp"


#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"

#ifdef HAVE_TEUCHOS_QD
#include <qd/dd_real.h>
#endif

namespace std { 


template <typename Packet>
ostream & operator<< ( ostream& os, const pair<Packet, Packet>& arg)
{
  os << "(" << arg.first << "," << arg.second << ")";
  return os;
}


} // namespace std


namespace Teuchos {


template<typename Packet>
struct ScalarTraits<std::pair<Packet,Packet> >
{
  typedef ScalarTraits<Packet> PST;
      typedef  std::pair<typename PST::magnitudeType, typename PST::magnitudeType> magnitudeType;
  static const bool isComplex = PST::isComplex;
  static const bool isComparable = PST::isComparable;
  static const bool hasMachineParameters = PST::hasMachineParameters;
  // Not defined: eps(), sfmin(), base(), prec(), t(), rnd(), emin(), rmin(), emax(), rmax()
  static inline magnitudeType magnitude(std::pair<Packet,Packet> a) { return std::pair<Packet,Packet>( PST::magnitude(a.first), PST::magnitude(a.second) ); }
  static inline std::pair<Packet,Packet> zero()  { return std::pair<Packet,Packet>(PST::zero(),PST::zero()); }
  static inline std::pair<Packet,Packet> one()   { return std::pair<Packet,Packet>(PST::one(), PST::one()); }
  static inline std::pair<Packet,Packet> conjugate(std::pair<Packet,Packet> x) { return std::pair<Packet,Packet>(PST::conjugate(x.first), PST::conjugate(x.second) ); }
  static inline std::pair<Packet,Packet> real(std::pair<Packet,Packet> x) { return std::pair<Packet,Packet>(PST::real(x.first), PST::real(x.second) ); }
  static inline std::pair<Packet,Packet> imag(std::pair<Packet,Packet> x) { return std::pair<Packet,Packet>(PST::imag(x.first), PST::imag(x.second) ); }
  static inline bool isnaninf(std::pair<Packet,Packet> x) { return PST::isnaninf(x.first) || PST::isnaninf(x.second); }
  static inline void seedrandom(unsigned int s) { PST::seedrandom(s); }
  static inline std::pair<Packet,Packet> random() { return std::pair<Packet,Packet>( PST::random(), PST::random() ); }
  static inline std::string name() { return "std::pair<" + Teuchos::TypeNameTraits<Packet>::name() + "," + Teuchos::TypeNameTraits<Packet>::name() + ">"; }
  static inline std::pair<Packet,Packet> squareroot(std::pair<Packet,Packet> x) { return std::pair<Packet,Packet>(PST::squareroot(x.first), PST::squareroot(x.second)); }
  static inline std::pair<Packet,Packet> pow(std::pair<Packet,Packet> x, std::pair<Packet,Packet> y) { return std::pair<Packet,Packet>( PST::pow(x.first,y.first), PST::pow(x.second,y.second) ); }
};

template<class Packet, class ConvertToPacket>
class ValueTypeConversionTraits<std::pair<Packet,Packet>, ConvertToPacket> {
public:
  static std::pair<Packet,Packet> convert( const ConvertToPacket t )
    {
      return std::pair<Packet,Packet>(t,t);
    }
  static std::pair<Packet,Packet> safeConvert( const ConvertToPacket t )
    {
      return std::pair<Packet,Packet>(t,t);
    }
};


} // namespace Teuchos


namespace {


using Teuchos::as;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::GlobalMPISession;
using Teuchos::defaultSmallNumber;
using Teuchos::outArg;


bool testMpi = true;


double errorTolSlack = 1e+1;



TEUCHOS_STATIC_SETUP()
{

  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();

  clp.addOutputSetupOptions(true);

  clp.setOption(
    "test-mpi", "test-serial", &testMpi,
    "Test MPI (if available) or force test of serial.  In a serial build,"
    " this option is ignored and a serial comm is always used." );

  clp.setOption(
    "error-tol-slack", &errorTolSlack,
    "Slack off of machine epsilon used to check test results" );

}


template<class Ordinal>
RCP<const Comm<Ordinal> > getDefaultComm()
{
  if (testMpi) {
    return DefaultComm<Ordinal>::getComm();
  }
  return rcp(new Teuchos::SerialComm<Ordinal>);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMpiComm, basic, Ordinal )
{
  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  out << "comm = " << Teuchos::describe(*comm);
  TEST_EQUALITY( size(*comm), GlobalMPISession::getNProc() );
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, reduceAllAndScatter_1, Ordinal, Packet )
{

  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);

#ifdef TEUCHOS_MPI_COMM_DUMP
  Teuchos::MpiComm<Ordinal>::show_dump = true;
#endif

  Array<Packet> sendBuffer(as<Ordinal>(numProcs));
  for (Ordinal k = 0; k < numProcs; ++k) {
    sendBuffer[k] = as<Packet>(1);
  }

  Array<Ordinal> recvCounts(as<Ordinal>(numProcs), as<Ordinal>(1));

  Array<Packet> myGlobalReducts(1);

  Teuchos::reduceAllAndScatter<Ordinal,Packet>(
    *comm, Teuchos::REDUCE_SUM,
    as<Ordinal>(sendBuffer.size()), &sendBuffer[0],
    &recvCounts[0], &myGlobalReducts[0]
    );

  if (std::numeric_limits<Packet>::is_integer) {
    TEST_EQUALITY( myGlobalReducts[0], as<Packet>(numProcs) );
  }
  else {
    const PacketMag local_errorTolSlack = static_cast<PacketMag>(errorTolSlack);
    TEST_FLOATING_EQUALITY( myGlobalReducts[0], as<Packet>(numProcs),
      as<PacketMag>(defaultSmallNumber<PacketMag>() * local_errorTolSlack / numProcs)
      );
  }
  
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, reduceAllAndScatter_2, Ordinal, Packet )
{

  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  Array<Packet> sendBuffer(as<Ordinal>(numProcs));
  for (Ordinal k = 0; k < numProcs; ++k) {
    sendBuffer[k] = as<Packet>(procRank + k);
  }

  Array<Ordinal> recvCounts(as<Ordinal>(numProcs), as<Ordinal>(1));

  Array<Packet> myGlobalReducts(1);

  Teuchos::reduceAllAndScatter<Ordinal,Packet>(
    *comm, Teuchos::REDUCE_SUM,
    as<Ordinal>(sendBuffer.size()), &sendBuffer[0],
    &recvCounts[0], &myGlobalReducts[0]
    );

  const Packet expectedMyGlobalReduct = as<Packet>(
    numProcs * procRank + ((numProcs - 1) * numProcs)/2 
    );

  if (std::numeric_limits<Packet>::is_integer) {
    TEST_EQUALITY( myGlobalReducts[0], expectedMyGlobalReduct );
  }
  else {
    const PacketMag local_errorTolSlack = static_cast<PacketMag>(errorTolSlack);
    TEST_FLOATING_EQUALITY( myGlobalReducts[0], expectedMyGlobalReduct,
      as<PacketMag>(defaultSmallNumber<PacketMag>() * local_errorTolSlack / numProcs)
      );
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, ReadySend1, Ordinal, Packet )
{

  using Teuchos::broadcast;
  using Teuchos::readySend;
  using Teuchos::wait;
  using Teuchos::as;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  using Teuchos::isend;
  using Teuchos::ireceive;
  using Teuchos::wait;
  using Teuchos::SerialComm;
  using Teuchos::is_null;
  using Teuchos::arcp;
  using Teuchos::arcpClone;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ArrayRCP;
  using Teuchos::ptr;
  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  if (
    numProcs == 1
    &&
    !is_null(rcp_dynamic_cast<const SerialComm<Ordinal> >(comm))
    )
  {
    out << "\nThis is Teuchos::SerialComm which does not support readySend!\n";
    return; // Pass!
  }

  PT::seedrandom(as<unsigned int>(procRank));
  Packet origSendData = PT::random();
  Packet origRecvData = PT::random();
  broadcast<Ordinal, Packet>( *comm, 0, outArg(origSendData) );

  Packet sendData = origSendData;
  Packet recvData = origRecvData;

  RCP<Teuchos::CommRequest<Ordinal> > recvRequest;

  // Post non-block receive on proc 0
  if (procRank == 0) {
    // Post non-blocking receive from proc n-1
    recvRequest = ireceive<Ordinal, Packet>(
        *comm,
        rcp(&recvData,false),
        numProcs-1
        );
  }
  barrier(*comm);

  if (procRank == numProcs-1) {
    // ready send from proc n-1 to proc 0
    // send data in sendData
    readySend<Ordinal, Packet>(
        *comm,
        sendData,
        0
        );
  }
  barrier(*comm);

  if (procRank == 0) {
    // wait for request on 0
    wait( *comm, outArg(recvRequest) );
  }
  barrier(*comm);

  // test that all procs have recvRequest == Teuchos::null
  TEST_EQUALITY_CONST( recvRequest, Teuchos::null );

  // proc 0 should have recvData == sendData
  if (procRank == 0) {
    TEST_EQUALITY( recvData, sendData );
  }
  // other procs should have recvData == origRecvData (i.e., unchanged)
  else {
    TEST_EQUALITY( recvData, origRecvData );
  }
  // all procs should have sendData == origSendData
  TEST_EQUALITY( sendData, origSendData );

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, ReadySend, Ordinal, Packet )
{

  using Teuchos::broadcast;
  using Teuchos::readySend;
  using Teuchos::wait;
  using Teuchos::as;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  using Teuchos::isend;
  using Teuchos::ireceive;
  using Teuchos::wait;
  using Teuchos::SerialComm;
  using Teuchos::is_null;
  using Teuchos::arcp;
  using Teuchos::arcpClone;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ArrayRCP;
  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  if (
    numProcs == 1
    &&
    !is_null(rcp_dynamic_cast<const SerialComm<Ordinal> >(comm))
    )
  {
    out << "\nThis is Teuchos::SerialComm which does not support readySend!\n";
    return; // Pass!
  }

  const int dataLen = 3;

  const ArrayRCP<Packet> origSendData = arcp<Packet>(dataLen);
  const ArrayRCP<Packet> origRecvData = arcp<Packet>(dataLen);
  PT::seedrandom(as<unsigned int>(procRank));
  for (int j = 0; j < dataLen; ++j) {
    origSendData[j] = PT::random();
    origRecvData[j] = PT::random();
  }
  broadcast<Ordinal, Packet>( *comm, 0, origSendData() );

  const ArrayRCP<Packet> sendData = arcpClone<Packet>(origSendData());
  const ArrayRCP<Packet> recvData = arcpClone<Packet>(origRecvData());

  RCP<Teuchos::CommRequest<Ordinal> > recvRequest;

  // both proc 0 and proc n-1 will post non-block receives, into recvData
  // then proc 0 will initiate a ready-send to proc n-1; the latter will initiate a wait
  // a barrier
  // then proc n-1 will initiate a ready-send to proc 0 of the data from recvData; the latter will initiate a wait
  // a barrier
  // now both of these procs should have matching data in sendData and recvData

  // Post non-block receive on both procs
  if (procRank == 0) {
    // Post non-blocking receive from proc n-1
    recvRequest = ireceive<Ordinal, Packet>(
        *comm,
        recvData.persistingView(0, dataLen),
        numProcs-1
        );
  }
  else if (procRank == numProcs-1) {
    // Post non-blocking receive from proc 0
    recvRequest = ireceive<Ordinal, Packet>(
        *comm,
        recvData.persistingView(0, dataLen),
        0
        );
  }
  barrier(*comm);

  if (procRank == 0) {
    // ready send from proc 0 to proc n-1
    // send data in sendData
    readySend<Ordinal, Packet>(
        *comm,
        sendData(),
        numProcs-1
        );
  }
  else if (procRank == numProcs-1) {
    // wait for request on proc n-1
    wait( *comm, outArg(recvRequest) );
  }
  barrier(*comm);

  if (procRank == 0) {
    // wait for request on 0
    wait( *comm, outArg(recvRequest) );
  }
  else if (procRank == numProcs-1) {
    // ready send from proc n-1 to proc 0
    // send data in recvData: THIS IS IMPORTANT: SEE ABOVE
    readySend<Ordinal, Packet>(
        *comm,
        recvData(),
        0
        );
  }
  barrier(*comm);

  // test that all procs (even non-participating) have recvRequest == Teuchos::null
  TEST_EQUALITY_CONST( recvRequest, Teuchos::null );

  // participating procs should have recvData == sendData
  if (procRank == 0 || procRank == numProcs-1) {
    TEST_COMPARE_ARRAYS( recvData, sendData );
  }
  // non-participating procs should have recvData == origRecvData (i.e., unchanged)
  else {
    TEST_COMPARE_ARRAYS( recvData, origRecvData );
  }
  // all procs should have sendData == origSendData
  TEST_COMPARE_ARRAYS( sendData, origSendData );

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, NonblockingSendReceive, Ordinal, Packet )
{

  using Teuchos::as;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  using Teuchos::isend;
  using Teuchos::ireceive;
  using Teuchos::wait;
  using Teuchos::SerialComm;
  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  if (
    numProcs == 1
    &&
    !is_null(rcp_dynamic_cast<const SerialComm<Ordinal> >(comm))
    )
  {
    out << "\nThis is Teuchos::SerialComm which does not yet support isend/ireceive!\n";
    return; // Pass!
  }

  // Only use randomize on one proc and then broacast
  Packet orig_input_data = PT::random();
  broadcast( *comm, 0, &orig_input_data );

  const Packet orig_output_data = as<Packet>(-1);

  const Packet input_data = orig_input_data;
  Packet output_data = orig_output_data;

  RCP<Teuchos::CommRequest<Ordinal> > recvRequest;
  RCP<Teuchos::CommRequest<Ordinal> > sendRequest;

  if (procRank == 0) {
    // Create copy of data to make sure that persisting relationship is
    // maintained!
    sendRequest = isend<Ordinal, Packet>(
      *comm, Teuchos::rcp(new Packet(input_data)), numProcs-1);
  }
  if (procRank == numProcs-1) {
    // We will need to read output_data after wait(...) below
    recvRequest = ireceive<Ordinal, Packet>(
      *comm, rcpFromRef(output_data), 0);
  }

  if (procRank == 0) {
    wait( *comm, outArg(sendRequest) );
  }
  if (procRank == numProcs-1) {
    wait( *comm, outArg(recvRequest) );
  }
  
  TEST_EQUALITY_CONST( sendRequest, Teuchos::null );
  TEST_EQUALITY_CONST( recvRequest, Teuchos::null );

  if (procRank == numProcs-1) {
    TEST_EQUALITY( output_data, input_data );
  }
  else {
    TEST_EQUALITY( output_data, orig_output_data );
  }
  TEST_EQUALITY( input_data, orig_input_data );

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );

}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( DefaultMpiComm, NonblockingSendReceiveSet, Ordinal, Packet )
{

  using Teuchos::as;
  using Teuchos::rcpFromRef;
  using Teuchos::outArg;
  using Teuchos::arcp;
  using Teuchos::arcpClone;
  using Teuchos::ArrayRCP;
  using Teuchos::isend;
  using Teuchos::ireceive;
  using Teuchos::wait;
  using Teuchos::broadcast;
  using Teuchos::SerialComm;
  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Packet> PT;
  typedef typename PT::magnitudeType PacketMag;
  typedef Teuchos::ScalarTraits<PacketMag> PMT;

  RCP<const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  const Ordinal numProcs = size(*comm);
  const Ordinal procRank = rank(*comm);

  if (
    numProcs == 1
    &&
    !is_null(rcp_dynamic_cast<const SerialComm<Ordinal> >(comm))
    )
  {
    out << "\nThis is Teuchos::SerialComm which does not yet support isend/ireceive!\n";
    return; // Pass!
  }

  const int numSendRecv = 4;
  const int sendLen = 3;

  const ArrayRCP<Packet> origInputData = arcp<Packet>(numSendRecv*sendLen);
  const ArrayRCP<Packet> origOutputData = arcp<Packet>(numSendRecv*sendLen);
  {
    int offset = 0;
    for (int i = 0; i < numSendRecv; ++i, offset += sendLen) {
      const ArrayRCP<Packet> origInputData_i =
        origInputData.persistingView(offset, sendLen); 
      const ArrayRCP<Packet> origOutputData_i =
        origOutputData.persistingView(offset, sendLen); 
      for (int j = 0; j < sendLen; ++j) {
        origInputData_i[j] = PT::random();
        origOutputData_i[j] = PT::random();
      }
    }
  }
  broadcast<Ordinal, Packet>( *comm, 0, origInputData() );

  const ArrayRCP<Packet> inputData = arcpClone<Packet>(origInputData());
  const ArrayRCP<Packet> outputData = arcpClone<Packet>(origOutputData());

  Array<RCP<Teuchos::CommRequest<Ordinal> > > recvRequests;
  Array<RCP<Teuchos::CommRequest<Ordinal> > > sendRequests;

  // Send from proc 0 to proc numProcs-1
  if (procRank == 0) {
    // Create copy of data to make sure that persisting relationship is
    // maintained!
    int offset = 0;
    for (int i = 0; i < numSendRecv; ++i, offset += sendLen) {
      sendRequests.push_back(
        isend<Ordinal, Packet>(
          *comm,
          arcpClone<Packet>(inputData(offset, sendLen)),
          numProcs-1
          )
        );
    }
  }

  // Receive from proc 0 on proc numProcs-1
  if (procRank == numProcs-1) {
    // We will need to read output_data after wait(...) below
    int offset = 0;
    for (int i = 0; i < numSendRecv; ++i, offset += sendLen) {
      recvRequests.push_back(
        ireceive<Ordinal, Packet>(
          *comm, outputData.persistingView(offset, sendLen), 0
          )
        );
    }
  }

  if (procRank == 0) {
    waitAll( *comm, sendRequests() );
  }
  if (procRank == numProcs-1) {
    waitAll( *comm, recvRequests() );
  }

  if (!sendRequests.empty()) {
    for (int i = 0; i < numSendRecv; ++i) {
      TEST_EQUALITY_CONST( sendRequests[i], Teuchos::null );
    }
  }

  if (!recvRequests.empty()) {
    for (int i = 0; i < numSendRecv; ++i) {
      TEST_EQUALITY_CONST( recvRequests[i], Teuchos::null );
    }
  }
  // ToDo: Write a test macro for this in one shot!

  if (procRank == numProcs-1) {
    TEST_COMPARE_ARRAYS( outputData, inputData );
  }
  else {
    TEST_COMPARE_ARRAYS( outputData, origOutputData );
  }
  TEST_COMPARE_ARRAYS( inputData, origInputData );

  // All procs fail if any proc fails
  int globalSuccess_int = -1;
  reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultMpiComm, duplicate, Ordinal)
{
  RCP< const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  int initialRank = comm->getRank();
  int initialSize = comm->getSize();

  RCP< const Comm<Ordinal> > newComm = comm->duplicate();
  TEST_EQUALITY(newComm->getSize(), initialSize);
  TEST_EQUALITY(newComm->getRank(), initialRank);

  // TODO Make sure the communication space is distinct.
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultMpiComm, split, Ordinal) {
  RCP< const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  int initialRank = comm->getRank();
  int initialSize = comm->getSize();

  // Partition this communicator into two: one with the odd ranks and one with
  // the even ones. Pass a common key for everyone to maintain the same
  // ordering as in the initial communicator.
  RCP< const Comm<Ordinal> > newComm = comm->split(initialRank % 2, 0);

  // Check the size of the new communicator and my rank within it.
  int halfSize = initialSize / 2;
  int newSize = newComm->getSize();
  int newRank = newComm->getRank();
  if (initialSize % 2 == 0) {
    TEST_EQUALITY(newSize, halfSize);
  }
  else {
    TEST_EQUALITY(newSize, initialRank % 2 == 0 ? halfSize + 1 : halfSize);
  }
  TEST_EQUALITY(newRank, initialRank / 2);

  // Negative color values get a null communicator.
  RCP< const Comm<Ordinal> > shouldBeNull = comm->split(-1, 0);
  TEST_ASSERT(shouldBeNull.is_null());
}

namespace {

template<typename ValueType>
class MonotoneSequence
{
  ValueType currentValue_;
public:
  typedef ValueType value_type;

  MonotoneSequence(const value_type& initialValue) : currentValue_(initialValue)
  {}

  value_type operator()()
  {
    return currentValue_++;
  }
};

} // namepsace

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(DefaultMpiComm, createSubcommunicator, Ordinal) {
  RCP< const Comm<Ordinal> > comm = getDefaultComm<Ordinal>();
  int initialRank = comm->getRank();
  int initialSize = comm->getSize();

  // Create a new communicator that reverses all of the ranks.
  std::vector< int > ranks(initialSize);
  std::generate(ranks.begin(), ranks.end(), MonotoneSequence<int>(0));
  std::reverse(ranks.begin(), ranks.end());
  RCP< const Comm<Ordinal> > newComm = comm->createSubcommunicator(ranks);
  TEST_EQUALITY(newComm->getSize(), initialSize);
  int expectedNewRank = initialSize - initialRank - 1;
  TEST_EQUALITY(newComm->getRank(), expectedNewRank);

  // Processes that aren't in the group get a null communicator.
  std::vector<int> rank0Only(1, 0);
  RCP< const Comm<Ordinal> > rank0Comm = comm->createSubcommunicator(rank0Only);
  // Original rank 0 should be valid, all others should be null.
  if (initialRank == 0) {
    TEST_ASSERT(rank0Comm.is_valid_ptr());
  } else {
    TEST_ASSERT(rank0Comm.is_null());
  }
}

//
// Instantiations
//


#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME, ORDINAL)\
     typedef std::complex<float> ComplexFloat; \
     TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(TEST_GROUP, TEST_NAME, ORDINAL, ComplexFloat)
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME, ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(TEST_GROUP, TEST_NAME, ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME, ORDINAL)
#  define UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME, ORDINAL)
#endif


#define UNIT_TEST_GROUP_ORDINAL_PACKET( ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, reduceAllAndScatter_1, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, reduceAllAndScatter_2, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, NonblockingSendReceive, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, NonblockingSendReceiveSet, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, ReadySend1, ORDINAL, PACKET ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, ReadySend, ORDINAL, PACKET )

#ifdef HAVE_TEUCHOS_QD
#  define UNIT_TEST_GROUP_ORDINAL_QD(ORDINAL) \
     UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, dd_real) \
     UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, qd_real)
#else
#  define UNIT_TEST_GROUP_ORDINAL_QD(ORDINAL)
#endif

#define UNIT_TEST_GROUP_ORDINAL_PAIROFPACKETS( ORDINAL, PAIROFPACKETS ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, NonblockingSendReceive, ORDINAL, PAIROFPACKETS ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, NonblockingSendReceiveSet, ORDINAL, PAIROFPACKETS ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, ReadySend1, ORDINAL, PAIROFPACKETS ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( DefaultMpiComm, ReadySend, ORDINAL, PAIROFPACKETS )

#define UNIT_TEST_GROUP_ORDINAL_SUBCOMMUNICATORS( ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DefaultMpiComm, duplicate, ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DefaultMpiComm, split, ORDINAL ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DefaultMpiComm, createSubcommunicator, ORDINAL )


typedef std::pair<short, short> PairOfShorts;
typedef std::pair<int,int> PairOfInts;
typedef std::pair<float,float> PairOfFloats;
typedef std::pair<double,double> PairOfDoubles;


// Uncomment this for really fast development cycles but make sure to comment
// it back again before checking in so that we can test all the types.
// #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#  define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, double) \
    UNIT_TEST_GROUP_ORDINAL_PAIROFPACKETS(ORDINAL, PairOfDoubles) \

  UNIT_TEST_GROUP_ORDINAL(int)

#else // FAST_DEVELOPMENT_UNIT_TEST_BUILD

#  define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DefaultMpiComm, basic, ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, short) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, int) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, float) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, double) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, reduceAllAndScatter_1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, reduceAllAndScatter_2, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, NonblockingSendReceive, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, ReadySend1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, ReadySend, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, reduceAllAndScatter_1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, reduceAllAndScatter_2, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, NonblockingSendReceive, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, ReadySend1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, ReadySend, ORDINAL) \
    UNIT_TEST_GROUP_ORDINAL_SUBCOMMUNICATORS(ORDINAL)

#  define UNIT_TEST_GROUP_ORDINAL_WITH_PAIRS_AND_QD( ORDINAL ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( DefaultMpiComm, basic, ORDINAL ) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, short)			\
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, int) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, float) \
    UNIT_TEST_GROUP_ORDINAL_PACKET(ORDINAL, double) \
    UNIT_TEST_GROUP_ORDINAL_QD(ORDINAL) \
    UNIT_TEST_GROUP_ORDINAL_PAIROFPACKETS(ORDINAL, PairOfShorts) \
    UNIT_TEST_GROUP_ORDINAL_PAIROFPACKETS(ORDINAL, PairOfInts) \
    UNIT_TEST_GROUP_ORDINAL_PAIROFPACKETS(ORDINAL, PairOfFloats) \
    UNIT_TEST_GROUP_ORDINAL_PAIROFPACKETS(ORDINAL, PairOfDoubles) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, reduceAllAndScatter_1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, reduceAllAndScatter_2, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, NonblockingSendReceive, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, ReadySend1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_FLOAT(DefaultMpiComm, ReadySend, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, reduceAllAndScatter_1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, reduceAllAndScatter_2, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, NonblockingSendReceive, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, ReadySend1, ORDINAL) \
    UNIT_TEST_TEMPLATE_2_INSTANT_COMPLEX_DOUBLE(DefaultMpiComm, ReadySend, ORDINAL)

  typedef short int ShortInt;
  UNIT_TEST_GROUP_ORDINAL(ShortInt)
  UNIT_TEST_GROUP_ORDINAL_WITH_PAIRS_AND_QD(int)
  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINAL(LongInt) // can't do QD with LongInt, one of the tests complains
  
#  ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#  endif

#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD


} // namespace

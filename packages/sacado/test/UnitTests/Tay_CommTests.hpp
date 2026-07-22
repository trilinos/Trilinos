// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"

#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ValueTypeSerializer;

template <typename TayType>
bool checkFadArrays(const Teuchos::Array<TayType>& x, 
		    const Teuchos::Array<TayType>& x2, 
		    const std::string& tag,
		    Teuchos::FancyOStream& out) {

  // Check sizes match
  bool success = (x.size() == x2.size());
  out << tag << " Taylor array size test";
  if (success)
    out << " passed";
  else
    out << " failed";
  out << ":  \n\tExpected:  " << x.size() << ", \n\tGot:       " << x2.size() 
      << "." << std::endl;
  
  // Check coefficients match
  for (int i=0; i<x.size(); i++) {
    bool success2 = Sacado::IsEqual<TayType>::eval(x[i], x2[i]);
    out << tag << " Taylor array comparison test " << i;
    if (success2)
      out << " passed";
    else
	out << " failed";
    out << ":  \n\tExpected:  " << x[i] << ", \n\tGot:       " << x2[i] << "." 
	<< std::endl;
    success = success && success2;
  }

  return success;
}

template<typename Ordinal>
bool checkResultOnAllProcs(
  const Teuchos::Comm<Ordinal> &comm,
  Teuchos::FancyOStream &out,
  const bool result
  )
{
  out << "\nChecking that the above test passed in all processes ...";
  int thisResult = ( result ? 1 : 0 );
  int sumResult = -1;
  Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM,Ordinal(1),&thisResult,
		     &sumResult);
  const bool passed = sumResult==Teuchos::size(comm);
  if(passed)
    out << " passed\n";
  else
    out << " (sumResult="<<sumResult<<"!=numProcs="<<Teuchos::size(comm)<<") failed\n";
  return passed;
}

#define TAY_COMM_TESTS(TayType, TAY)					\
TEUCHOS_UNIT_TEST( TAY##_Comm, Broadcast ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), x2(n), x3(n);				\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, rnd.number());					\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = rnd.number();				\
  }									\
  for (int i=0; i<n; i++) {						\
    x2[i] = TayType(p, 0.0);						\
  }									\
  if (comm->getRank() == 0) {						\
    x2 = x;								\
    x3 = x;								\
  }									\
									\
  Teuchos::broadcast(*comm, 0, n, &x2[0]);				\
  bool success1 = checkFadArrays(					\
    x, x2, std::string(#TAY)+" Broadcast", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
                                                                        \
  Teuchos::broadcast(*comm, tts, 0, n, &x3[0]);				\
  bool success2 = checkFadArrays(					\
    x, x3, std::string(#TAY)+" Broadcast TTS", out);			\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, GatherAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int size = comm->getSize();						\
  int rank = comm->getRank();						\
  int N = n*size;							\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), x2(N), x3(N), x4(N);			\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, (rank+1)*(i+1));					\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = (rank+1)*(i+1)*(j+1);			\
  }									\
  for (int i=0; i<N; i++) {						\
    x2[i] = TayType(p, 0.0);						\
  }									\
  for (int j=0; j<size; j++) {						\
    for (int i=0; i<n; i++) {						\
      x3[n*j+i] = TayType(p, (j+1)*(i+1));				\
      for (int k=0; k<=p; k++)						\
	x3[n*j+i].fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);		\
    }									\
  }									\
									\
  Teuchos::gatherAll(*comm, n, &x[0], N, &x2[0]);			\
  bool success1 = checkFadArrays(					\
    x3, x2, std::string(#TAY)+" Gather All", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::gatherAll(*comm, tts, n, &x[0], N, &x4[0]);			\
  bool success2 = checkFadArrays(					\
    x3, x4, std::string(#TAY)+" Gather All TTS", out);			\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, SumAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int num_proc = comm->getSize();					\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), sums(n), sums2(n), sums3(n);		\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, 1.0*(i+1));					\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  for (int i=0; i<n; i++) {						\
    sums[i] = TayType(p, 1.0*(i+1)*num_proc);				\
    for (int j=0; j<=p; j++)						\
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;			\
  }									\
  for (int i=0; i<n; i++) {						\
    sums2[i] = TayType(p, 0.0);						\
  }									\
									\
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);	\
  bool success1 = checkFadArrays(					\
    sums, sums2, std::string(#TAY)+" Sum All", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::reduceAll(*comm, tts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]); \
  bool success2 = checkFadArrays(					\
    sums, sums3, std::string(#TAY)+" Sum All TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, MaxAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
  int num_proc = comm->getSize();					\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), maxs(n), maxs2(n), maxs3(n);		\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, 1.0*(i+1)*(rank+1));				\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs[i] = TayType(p, 1.0*(i+1)*num_proc);				\
    for (int j=0; j<=p; j++)						\
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;			\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs2[i] = TayType(p, 0.0);						\
  }									\
									\
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);	\
  bool success1 = checkFadArrays(					\
    maxs, maxs2, std::string(#TAY)+" Max All", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::reduceAll(*comm, tts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]); \
  bool success2 = checkFadArrays(					\
    maxs, maxs3, std::string(#TAY)+" Max All TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, MinAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), mins(n), mins2(n), mins3(n);		\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, 1.0*(i+1)*(rank+1));				\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    mins[i] = TayType(p, 1.0*(i+1));					\
    for (int j=0; j<=p; j++)						\
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  for (int i=0; i<n; i++) {						\
    mins2[i] = TayType(p, 0.0);						\
  }									\
									\
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);	\
  bool success1 = checkFadArrays(					\
    mins, mins2, std::string(#TAY)+" Min All", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::reduceAll(*comm, tts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]); \
  bool success2 = checkFadArrays(					\
    mins, mins3, std::string(#TAY)+" Min All TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, ScanSum ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), sums(n), sums2(n), sums3(n);		\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, 1.0*(i+1));					\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  for (int i=0; i<n; i++) {						\
    sums[i] = TayType(p, 1.0*(i+1)*(rank+1));				\
    for (int j=0; j<=p; j++)						\
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    sums2[i] = TayType(p, 0.0);						\
  }									\
									\
  Teuchos::scan(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);	\
  bool success1 = checkFadArrays(					\
    sums, sums2, std::string(#TAY)+" Scan Sum", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::scan(*comm, tts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]);	\
  bool success2 = checkFadArrays(					\
    sums, sums3, std::string(#TAY)+" Scan Sum TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
  									\
  success = success1 && success2;					\
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, ScanMax ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), maxs(n), maxs2(n), maxs3(n);		\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, 1.0*(i+1)*(rank+1));				\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs[i] = TayType(p, 1.0*(i+1)*(rank+1));				\
    for (int j=0; j<=p; j++)						\
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs2[i] = TayType(p, 0.0);						\
  }									\
									\
  Teuchos::scan(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);	\
  bool success1 = checkFadArrays(					\
    maxs, maxs2, std::string(#TAY)+" Scan Max", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::scan(*comm, tts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]);	\
  bool success2 = checkFadArrays(					\
    maxs, maxs3, std::string(#TAY)+" Scan Max TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, ScanMin ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
  ValueTypeSerializer<int,TayType> tts(					\
    rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
  Teuchos::Array<TayType> x(n), mins(n), mins2(n), mins3(n);		\
  for (int i=0; i<n; i++) {						\
    x[i] = TayType(p, 1.0*(i+1)*(rank+1));				\
    for (int j=0; j<=p; j++)						\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    mins[i] = TayType(p, 1.0*(i+1));					\
    for (int j=0; j<=p; j++)						\
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  for (int i=0; i<n; i++) {						\
    mins2[i] = TayType(p, 0.0);						\
  }									\
									\
  Teuchos::scan(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);	\
  bool success1 = checkFadArrays(					\
    mins, mins2, std::string(#TAY)+" Scan Min", out);			\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::scan(*comm, tts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]);	\
  bool success2 = checkFadArrays(					\
    mins, mins3, std::string(#TAY)+" Scan Min TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, SendReceive ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int num_proc = comm->getSize();					\
  if (num_proc > 1) {							\
    int rank = comm->getRank();						\
    int n = 7;								\
    int p = 5;								\
    ValueTypeSerializer<int,TayType> tts(				\
      rcp(new ValueTypeSerializer<int,double>), p+1);			\
									\
    Teuchos::Array<TayType> x(n), x2(n), x3(n);				\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(p, 1.0*(i+1));					\
      for (int j=0; j<=p; j++)						\
	x[i].fastAccessCoeff(j) = 2.0*(i+1)*(j+1);			\
    }									\
    for (int i=0; i<n; i++) {						\
      x2[i] = TayType(p, 0.0);						\
    }									\
    if (rank != 1) {							\
      x2 = x;								\
      x3 = x2;								\
    }									\
    									\
    if (rank == 0) Teuchos::send(*comm, n, &x[0], 1);			\
    if (rank == 1) Teuchos::receive(*comm, 0, n, &x2[0]);		\
    bool success1 = checkFadArrays(					\
      x, x2, std::string(#TAY)+" Send/Receive", out);			\
    success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
    if (rank == 0) Teuchos::send(*comm, tts, n, &x[0], 1);		\
    if (rank == 1) Teuchos::receive(*comm, tts, 0, n, &x3[0]);		\
    bool success2 = checkFadArrays(					\
      x, x3, std::string(#TAY)+" Send/Receive TTS", out);		\
    success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
    success = success1 && success2;					\
  }									\
  else									\
    success = true;							\
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedBroadcast ) {			\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), x2(n), x3(n);			\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, rnd.number());					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = rnd.number();				\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
       TayType g(p1, rnd.number());					\
       for (int k=0; k<=p1; k++)					\
	 g.fastAccessCoeff(k) = rnd.number();				\
       x[i].fastAccessCoeff(j) = g;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    x2[i] = TayTayType(p2, TayType(p1, 0.0));				\
    for (int j=0; j<=p2; j++)						\
      x2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
  if (comm->getRank() == 0) {						\
    x2 = x;								\
    x3 = x;								\
  }									\
  									\
  Teuchos::broadcast(*comm, 0, n, &x2[0]);				\
  bool success1 = checkFadArrays(					\
    x, x2, std::string(#TAY)+"<"+#TAY+"> Broadcast", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::broadcast(*comm, ttts, 0, n, &x3[0]);			\
  bool success2 = checkFadArrays(					\
    x, x3, std::string(#TAY)+"<"+#TAY+"> Broadcast TTS", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedGatherAll ) {			\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  int size = comm->getSize();						\
  int rank = comm->getRank();						\
  int N = n*size;							\
  Teuchos::Array<TayTayType> x(n), x2(N), x3(N), x4(N);			\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, (rank+1)*(i+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = (rank+1)*(i+1)*(k+1);			\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<N; i++) {						\
    x2[i] = TayTayType(p2, TayType(p1, 0.0));				\
    for (int j=0; j<=p2; j++)						\
      x2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
  for (int j=0; j<size; j++) {						\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, (j+1)*(i+1));					\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);			\
      x3[n*j+i] = TayTayType(p2, f);					\
      for (int k=0; k<=p2; k++)						\
	x3[n*j+i].fastAccessCoeff(k) = f;				\
    }									\
  }									\
									\
  Teuchos::gatherAll(*comm, n, &x[0], N, &x2[0]);			\
  bool success1 = checkFadArrays(					\
    x3, x2, std::string(#TAY)+"<"+#TAY+">  Gather All", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::gatherAll(*comm, ttts, n, &x[0], N, &x4[0]);			\
  bool success2 = checkFadArrays(					\
    x3, x4, std::string(#TAY)+"<"+#TAY+">  Gather All FTS", out);	\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedSumAll ) {				\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  int num_proc = comm->getSize();					\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), sums(n), sums2(n), sums3(n);		\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1));						\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*num_proc);					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;			\
    sums[i] = TayTayType(p2, f);					\
    for (int j=0; j<=p2; j++)						\
      sums[i].fastAccessCoeff(j) = f;					\
  }									\
  for (int i=0; i<n; i++) {						\
    sums2[i] = TayTayType(p2, TayType(p1, 0.0));			\
    for (int j=0; j<=p2; j++)						\
      sums2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
									\
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);	\
  bool success1 = checkFadArrays(					\
    sums, sums2, std::string(#TAY)+"<"+#TAY+"> Sum All", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::reduceAll(*comm, ttts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]); \
  bool success2 = checkFadArrays(					\
    sums, sums3, std::string(#TAY)+"<"+#TAY+"> Sum All", out);		\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedMaxAll ) {				\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  int rank = comm->getRank();						\
  int num_proc = comm->getSize();					\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), maxs(n), maxs2(n), maxs3(n);		\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*(rank+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*num_proc);					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;			\
    maxs[i] = TayTayType(p2, f);					\
    for (int j=0; j<=p2; j++)						\
      maxs[i].fastAccessCoeff(j) = f;					\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs2[i] = TayTayType(p2, TayType(p1, 0.0));			\
    for (int j=0; j<=p2; j++)						\
      maxs2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
									\
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);	\
  bool success1 = checkFadArrays(					\
    maxs, maxs2, std::string(#TAY)+"<"+#TAY+"> Max All", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::reduceAll(*comm, ttts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]); \
  bool success2 = checkFadArrays(					\
    maxs, maxs3, std::string(#TAY)+"<"+#TAY+"> Max All FTS", out);	\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedMinAll ) {				\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  int rank = comm->getRank();						\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), mins(n), mins2(n), mins3(n);		\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*(rank+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1));						\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    mins[i] = TayTayType(p2, f);					\
    for (int j=0; j<=p2; j++)						\
      mins[i].fastAccessCoeff(j) = f;					\
  }									\
  for (int i=0; i<n; i++) {						\
    mins2[i] = TayTayType(p2, TayType(p1, 0.0));			\
    for (int j=0; j<=p2; j++)						\
      mins2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
									\
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);	\
  bool success1 = checkFadArrays(					\
    mins, mins2, std::string(#TAY)+"<"+#TAY+"> Min All", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::reduceAll(*comm, ttts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]); \
  bool success2 = checkFadArrays(					\
    mins, mins3, std::string(#TAY)+"<"+#TAY+"> Min All FTS", out);	\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedScanSum ) {			\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  int rank = comm->getRank();						\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), sums(n), sums2(n), sums3(n);		\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1));						\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*(rank+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    sums[i] = TayTayType(p2, f);					\
    for (int j=0; j<=p2; j++)						\
      sums[i].fastAccessCoeff(j) = f;					\
  }									\
  for (int i=0; i<n; i++) {						\
    sums2[i] = TayTayType(p2, TayType(p1, 0.0));			\
    for (int j=0; j<=p2; j++)						\
      sums2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
									\
  Teuchos::scan(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);	\
  bool success1 = checkFadArrays(					\
    sums, sums2, std::string(#TAY)+"<"+#TAY+"> Scan Sum", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::scan(*comm, ttts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]);	\
  bool success2 = checkFadArrays(					\
    sums, sums3, std::string(#TAY)+"<"+#TAY+"> Scan Sum FTS", out);	\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedScanMax ) {			\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  int rank = comm->getRank();						\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), maxs(n), maxs2(n), maxs3(n);		\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*(rank+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*(rank+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    maxs[i] = TayTayType(p2, f);					\
    for (int j=0; j<=p2; j++)						\
      maxs[i].fastAccessCoeff(j) = f;					\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs2[i] = TayTayType(p2, TayType(p1, 0.0));			\
    for (int j=0; j<=p2; j++)						\
      maxs2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
									\
  Teuchos::scan(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);	\
  bool success1 = checkFadArrays(					\
    maxs, maxs2, std::string(#TAY)+"<"+#TAY+"> Scan Max", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::scan(*comm, ttts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]);	\
  bool success2 = checkFadArrays(					\
    maxs, maxs3, std::string(#TAY)+"<"+#TAY+"> Scan Max FTS", out);	\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedScanMin ) {			\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p1 = 5;								\
  int p2 = 5;								\
  int rank = comm->getRank();						\
  RCP< ValueTypeSerializer<int,TayType> > tts =				\
    rcp(new ValueTypeSerializer<int,TayType>(				\
	  rcp(new ValueTypeSerializer<int,double>), p1+1));		\
  ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);			\
									\
  Teuchos::Array<TayTayType> x(n), mins(n), mins2(n), mins3(n);		\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1)*(rank+1));					\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = TayTayType(p2, f);						\
    for (int j=0; j<=p2; j++) {						\
      x[i].fastAccessCoeff(j) = f;					\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    TayType f(p1, 1.0*(i+1));						\
    for (int k=0; k<=p1; k++)						\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    mins[i] = TayTayType(p2, f);					\
    for (int j=0; j<=p2; j++)						\
      mins[i].fastAccessCoeff(j) = f;					\
  }									\
  for (int i=0; i<n; i++) {						\
    mins2[i] = TayTayType(p2, TayType(p1, 0.0));			\
    for (int j=0; j<=p2; j++)						\
      mins2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
  }									\
									\
  Teuchos::scan(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);	\
  bool success1 = checkFadArrays(					\
    mins, mins2, std::string(#TAY)+"<"+#TAY+"> Scan Min", out);		\
  success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
  Teuchos::scan(*comm, ttts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]);	\
  bool success2 = checkFadArrays(					\
    mins, mins3, std::string(#TAY)+"<"+#TAY+"> Scan Min FTS", out);	\
  success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
  success = success1 && success2;                                       \
}									\
									\
TEUCHOS_UNIT_TEST( TAY##_Comm, NestedSendReceive ) {			\
  typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int num_proc = comm->getSize();					\
  if (num_proc > 1) {							\
    int rank = comm->getRank();						\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 5;								\
    RCP< ValueTypeSerializer<int,TayType> > tts =			\
      rcp(new ValueTypeSerializer<int,TayType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1+1));		\
    ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);		\
									\
    Teuchos::Array<TayTayType> x(n), x2(n), x3(n);			\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, 1.0*(i+1));						\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = 2.0*(i+1)*(k+1);				\
      x[i] = TayTayType(p2, f);						\
      for (int j=0; j<=p2; j++)						\
	x[i].fastAccessCoeff(j) = f;					\
    }									\
    for (int i=0; i<n; i++) {						\
      x2[i] = TayTayType(p2, TayType(p1, 0.0));				\
      for (int j=0; j<=p2; j++)						\
	x2[i].fastAccessCoeff(j) = TayType(p1, 0.0);			\
    }									\
    if (rank != 1) {							\
      x2 = x;								\
      x3 = x2;								\
    }									\
									\
    if (rank == 0) Teuchos::send(*comm, n, &x[0], 1);			\
    if (rank == 1) Teuchos::receive(*comm, 0, n, &x2[0]);		\
    bool success1 = checkFadArrays(					\
      x, x2, std::string(#TAY)+"<"+#TAY+"> Send/Receive", out);		\
    success1 = checkResultOnAllProcs(*comm, out, success1);		\
									\
    if (rank == 0) Teuchos::send(*comm, ttts, n, &x[0], 1);		\
    if (rank == 1) Teuchos::receive(*comm, ttts, 0, n, &x3[0]);		\
    bool success2 = checkFadArrays(					\
      x, x3, std::string(#TAY)+"<"+#TAY+"> Send/Receive FTS", out);	\
    success2 = checkResultOnAllProcs(*comm, out, success2);		\
									\
    success = success1 && success2;					\
  }									\
  else									\
    success = true;							\
}

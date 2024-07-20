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

template <typename ArrayType>
bool checkFadArrays(const ArrayType& x,
                    const ArrayType& x2,
                    const std::string& tag,
                    Teuchos::FancyOStream& out) {
  typedef typename ArrayType::value_type FadType;

  // Check sizes match
  bool success = (x.size() == x2.size());
  out << tag << " Fad array size test";
  if (success)
    out << " passed";
  else
    out << " failed";
  out << ":  \n\tExpected:  " << x.size() << ", \n\tGot:       " << x2.size()
      << "." << std::endl;

  // Check Fads match
  const int sz = x.size();
  for (int i=0; i<sz; i++) {
    bool success2 = Sacado::IsEqual<FadType>::eval(x[i], x2[i]);
    out << tag << " Fad array comparison test " << i;
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

#define FAD_BASE_COMM_TESTS(FadType, FAD)                               \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_Broadcast ) {                        \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), x2(n), x3(n);                           \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, rnd.number());                                    \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = rnd.number();                              \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    x2[i] = FadType(p, 0.0);                                            \
  }                                                                     \
  if (comm->getRank() == 0) {                                           \
    x2 = x;                                                             \
    x3 = x;                                                             \
  }                                                                     \
                                                                        \
  Teuchos::broadcast(*comm, 0, n, &x2[0]);                              \
  bool success1 = checkFadArrays(                                       \
    x, x2, std::string(#FAD)+" Broadcast", out);                        \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::broadcast(*comm, fts, 0, n, &x3[0]);                         \
  bool success2 = checkFadArrays(                                       \
    x, x3, std::string(#FAD)+" Broadcast FTS", out);                    \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_GatherAll ) {                        \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int size = comm->getSize();                                           \
  int rank = comm->getRank();                                           \
  int N = n*size;                                                       \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
   Teuchos::Array<FadType> x(n), x2(N), x3(N), x4(N);                   \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, (rank+1)*(i+1));                                  \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = (rank+1)*(i+1)*(j+1);                      \
  }                                                                     \
  for (int i=0; i<N; i++) {                                             \
    x2[i] = FadType(p, 0.0);                                            \
  }                                                                     \
  for (int j=0; j<size; j++) {                                          \
    for (int i=0; i<n; i++) {                                           \
      x3[n*j+i] = FadType(p, (j+1)*(i+1));                              \
      for (int k=0; k<p; k++)                                           \
        x3[n*j+i].fastAccessDx(k) = (j+1)*(i+1)*(k+1);                  \
    }                                                                   \
  }                                                                     \
                                                                        \
  Teuchos::gatherAll(*comm, n, &x[0], N, &x2[0]);                       \
  bool success1 = checkFadArrays(                                       \
    x3, x2, std::string(#FAD)+" Gather All", out);                      \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::gatherAll(*comm, fts, n, &x[0], N, &x4[0]);                  \
  bool success2 = checkFadArrays(                                       \
    x3, x4, std::string(#FAD)+" Gather All FTS", out);                  \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_SumAll ) {                           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int num_proc = comm->getSize();                                       \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), sums(n), sums2(n), sums3(n);            \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, 1.0*(i+1));                                       \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = 2.0*(i+1);                                 \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums[i] = FadType(p, 1.0*(i+1)*num_proc);                           \
    for (int j=0; j<p; j++)                                             \
      sums[i].fastAccessDx(j) = 2.0*(i+1)*num_proc;                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums2[i] = FadType(p, 0.0);                                         \
  }                                                                     \
                                                                        \
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);  \
  bool success1 = checkFadArrays(                                       \
    sums, sums2, std::string(#FAD)+" Sum All", out);                    \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]); \
  bool success2 = checkFadArrays(                                       \
    sums, sums3, std::string(#FAD)+" Sum All FTS", out);                \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_MaxAll ) {                           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  int num_proc = comm->getSize();                                       \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), maxs(n), maxs2(n), maxs3(n);            \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, 1.0*(i+1)*(rank+1));                              \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                        \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs[i] = FadType(p, 1.0*(i+1)*num_proc);                           \
    for (int j=0; j<p; j++)                                             \
      maxs[i].fastAccessDx(j) = 2.0*(i+1)*num_proc;                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs2[i] = FadType(p, 0.0);                                         \
  }                                                                     \
                                                                        \
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);  \
  bool success1 = checkFadArrays(                                       \
    maxs, maxs2, std::string(#FAD)+" Max All", out);                    \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]); \
  bool success2 = checkFadArrays(                                       \
    maxs, maxs3, std::string(#FAD)+" Max All FTS", out);                \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_MinAll ) {                           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), mins(n), mins2(n), mins3(n);            \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, 1.0*(i+1)*(rank+1));                              \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                        \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins[i] = FadType(p, 1.0*(i+1));                                    \
    for (int j=0; j<p; j++)                                             \
      mins[i].fastAccessDx(j) = 2.0*(i+1);                              \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins2[i] = FadType(p, 0.0);                                         \
  }                                                                     \
                                                                        \
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);  \
  bool success1 = checkFadArrays(                                       \
    mins, mins2, std::string(#FAD)+" Min All", out);                    \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]); \
  bool success2 = checkFadArrays(                                       \
    mins, mins3, std::string(#FAD)+" Min All FTS", out);                \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_ScanSum ) {                          \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), sums(n), sums2(n), sums3(n);            \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, 1.0*(i+1));                                       \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = 2.0*(i+1);                                 \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums[i] = FadType(p, 1.0*(i+1)*(rank+1));                           \
    for (int j=0; j<p; j++)                                             \
      sums[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums2[i] = FadType(p, 0.0);                                         \
  }                                                                     \
                                                                        \
  Teuchos::scan(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);       \
  bool success1 = checkFadArrays(                                       \
    sums, sums2, std::string(#FAD)+" Scan Sum", out);                   \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::scan(*comm, fts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]);  \
  bool success2 = checkFadArrays(                                       \
    sums, sums3, std::string(#FAD)+" Scan Sum FTS", out);               \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_ScanMax ) {                          \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), maxs(n), maxs2(n), maxs3(n);            \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, 1.0*(i+1)*(rank+1));                              \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                        \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs[i] = FadType(p, 1.0*(i+1)*(rank+1));                           \
    for (int j=0; j<p; j++)                                             \
      maxs[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs2[i] = FadType(p, 0.0);                                         \
  }                                                                     \
                                                                        \
  Teuchos::scan(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);       \
  bool success1 = checkFadArrays(                                       \
    maxs, maxs2, std::string(#FAD)+" Scan Max", out);                   \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::scan(*comm, fts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]);  \
  bool success2 = checkFadArrays(                                       \
    maxs, maxs3, std::string(#FAD)+" Scan Max FTS", out);               \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_ScanMin ) {                          \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  Teuchos::Array<FadType> x(n), mins(n), mins2(n), mins3(n);            \
  for (int i=0; i<n; i++) {                                             \
    x[i] = FadType(p, 1.0*(i+1)*(rank+1));                              \
    for (int j=0; j<p; j++)                                             \
      x[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                        \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins[i] = FadType(p, 1.0*(i+1));                                    \
    for (int j=0; j<p; j++)                                             \
      mins[i].fastAccessDx(j) = 2.0*(i+1);                              \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins2[i] = FadType(p, 0.0);                                         \
  }                                                                     \
                                                                        \
  Teuchos::scan(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);       \
  bool success1 = checkFadArrays(                                       \
    mins, mins2, std::string(#FAD)+" Scan Min", out);                   \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::scan(*comm, fts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]);  \
  bool success2 = checkFadArrays(                                       \
    mins, mins3, std::string(#FAD)+" Scan Min FTS", out);               \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, Fad_SendReceive ) {                      \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int num_proc = comm->getSize();                                       \
  if (num_proc > 1) {                                                   \
    int rank = comm->getRank();                                         \
    int n = 7;                                                          \
    int p = 5;                                                          \
    ValueTypeSerializer<int,FadType> fts(                               \
      rcp(new ValueTypeSerializer<int,double>), p);                     \
                                                                        \
    Teuchos::Array<FadType> x(n), x2(n), x3(n);                         \
    for (int i=0; i<n; i++) {                                           \
      x[i] = FadType(p, 1.0*(i+1));                                     \
      for (int j=0; j<p; j++)                                           \
        x[i].fastAccessDx(j) = 2.0*(i+1)*(j+1);                         \
    }                                                                   \
    for (int i=0; i<n; i++) {                                           \
      x2[i] = FadType(p, 0.0);                                          \
    }                                                                   \
    if (rank != 1) {                                                    \
      x2 = x;                                                           \
      x3 = x;                                                           \
    }                                                                   \
                                                                        \
    if (rank == 0) Teuchos::send(*comm, n, &x[0], 1);                   \
    if (rank == 1) Teuchos::receive(*comm, 0, n, &x2[0]);               \
    bool success1 = checkFadArrays(                                     \
      x, x2, std::string(#FAD)+" Send/Receive", out);                   \
    success1 = checkResultOnAllProcs(*comm, out, success1);             \
                                                                        \
    if (rank == 0) Teuchos::send(*comm, fts, n, &x[0], 1);              \
    if (rank == 1) Teuchos::receive(*comm, fts, 0, n, &x3[0]);          \
    bool success2 = checkFadArrays(                                     \
      x, x3, std::string(#FAD)+" Send/Receive FTS", out);               \
    success2 = checkResultOnAllProcs(*comm, out, success2);             \
                                                                        \
    success = success1 && success2;                                     \
  }                                                                     \
  else                                                                  \
    success = true;                                                     \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_Broadcast ) {                     \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), x2(n), x3(n);                        \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, rnd.number());                                        \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = rnd.number();                                 \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
       FadType g(p1, rnd.number());                                     \
       for (int k=0; k<p1; k++)                                         \
         g.fastAccessDx(k) = rnd.number();                              \
       x[i].fastAccessDx(j) = g;                                        \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    x2[i] = FadFadType(p2, FadType(p1, 0.0));                           \
    for (int j=0; j<p2; j++)                                            \
      x2[i].fastAccessDx(j) = FadType(p1, 0.0);                         \
  }                                                                     \
  if (comm->getRank() == 0) {                                           \
    x2 = x;                                                             \
    x3 = x;                                                             \
  }                                                                     \
                                                                        \
  Teuchos::broadcast(*comm, 0, n, &x2[0]);                              \
  bool success1 = checkFadArrays(                                       \
    x, x2, std::string(#FAD)+"<"+#FAD+"> Broadcast", out);              \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::broadcast(*comm, ffts, 0, n, &x3[0]);                        \
  bool success2 = checkFadArrays(                                       \
    x, x3, std::string(#FAD)+"<"+#FAD+"> Broadcast FTS", out);          \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_GatherAll ) {                     \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int size = comm->getSize();                                           \
  int rank = comm->getRank();                                           \
  int N = n*size;                                                       \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), x2(N), x3(N), x4(N);                 \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, (rank+1)*(i+1));                                      \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = (rank+1)*(i+1)*(k+1);                         \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<N; i++) {                                             \
    x2[i] = FadFadType(p2, FadType(p1, 0.0));                           \
    for (int j=0; j<p2; j++)                                            \
      x2[i].fastAccessDx(j) = FadType(p1, 0.0);                         \
  }                                                                     \
  for (int j=0; j<size; j++) {                                          \
    for (int i=0; i<n; i++) {                                           \
      FadType f(p1, (j+1)*(i+1));                                       \
      for (int k=0; k<p1; k++)                                          \
        f.fastAccessDx(k) = (j+1)*(i+1)*(k+1);                          \
      x3[n*j+i] = FadFadType(p2, f);                                    \
      for (int k=0; k<p2; k++)                                          \
        x3[n*j+i].fastAccessDx(k) = f;                                  \
    }                                                                   \
  }                                                                     \
                                                                        \
  Teuchos::gatherAll(*comm, n, &x[0], N, &x2[0]);                       \
  bool success1 = checkFadArrays(                                       \
    x3, x2, std::string(#FAD)+"<"+#FAD+">  Gather All", out);           \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::gatherAll(*comm, ffts, n, &x[0], N, &x4[0]);                 \
  bool success2 = checkFadArrays(                                       \
    x3, x4, std::string(#FAD)+"<"+#FAD+">  Gather All FTS", out);       \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_SumAll ) {                        \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int num_proc = comm->getSize();                                       \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), sums(n), sums2(n), sums3(n);         \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1));                                           \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1);                                    \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*num_proc);                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*num_proc;                           \
    sums[i] = FadFadType(p2, f);                                        \
    for (int j=0; j<p2; j++)                                            \
      sums[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums2[i] = FadFadType(p2, FadType(p1, 0.0));                        \
    for (int j=0; j<p2; j++)                                            \
      sums2[i].fastAccessDx(j) = FadType(p1, 0.0);                      \
  }                                                                     \
                                                                        \
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);  \
  bool success1 = checkFadArrays(                                       \
    sums, sums2, std::string(#FAD)+"<"+#FAD+"> Sum All", out);          \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::reduceAll(*comm, ffts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]); \
  bool success2 = checkFadArrays(                                       \
    sums, sums3, std::string(#FAD)+"<"+#FAD+"> Sum All", out);          \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_MaxAll ) {                        \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int rank = comm->getRank();                                           \
  int num_proc = comm->getSize();                                       \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), maxs(n), maxs2(n), maxs3(n);         \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*(rank+1));                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*(rank+1);                           \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*num_proc);                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*num_proc;                           \
    maxs[i] = FadFadType(p2, f);                                        \
    for (int j=0; j<p2; j++)                                            \
      maxs[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs2[i] = FadFadType(p2, FadType(p1, 0.0));                        \
    for (int j=0; j<p2; j++)                                            \
      maxs2[i].fastAccessDx(j) = FadType(p1, 0.0);                      \
  }                                                                     \
                                                                        \
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);  \
  bool success1 = checkFadArrays(                                       \
    maxs, maxs2, std::string(#FAD)+"<"+#FAD+"> Max All", out);          \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::reduceAll(*comm, ffts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]); \
  bool success2 = checkFadArrays(                                       \
    maxs, maxs3, std::string(#FAD)+"<"+#FAD+"> Max All FTS", out);      \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_MinAll ) {                        \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int rank = comm->getRank();                                           \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), mins(n), mins2(n), mins3(n);         \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*(rank+1));                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*(rank+1);                           \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1));                                           \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1);                                    \
    mins[i] = FadFadType(p2, f);                                        \
    for (int j=0; j<p2; j++)                                            \
      mins[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins2[i] = FadFadType(p2, FadType(p1, 0.0));                        \
    for (int j=0; j<p2; j++)                                            \
      mins2[i].fastAccessDx(j) = FadType(p1, 0.0);                      \
  }                                                                     \
                                                                        \
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);  \
  bool success1 = checkFadArrays(                                       \
    mins, mins2, std::string(#FAD)+"<"+#FAD+"> Min All", out);          \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::reduceAll(*comm, ffts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]); \
  bool success2 = checkFadArrays(                                       \
    mins, mins3, std::string(#FAD)+"<"+#FAD+"> Min All FTS", out);      \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_ScanSum ) {                       \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int rank = comm->getRank();                                           \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), sums(n), sums2(n), sums3(n);         \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1));                                           \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1);                                    \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*(rank+1));                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*(rank+1);                           \
    sums[i] = FadFadType(p2, f);                                        \
    for (int j=0; j<p2; j++)                                            \
      sums[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums2[i] = FadFadType(p2, FadType(p1, 0.0));                        \
    for (int j=0; j<p2; j++)                                            \
      sums2[i].fastAccessDx(j) = FadType(p1, 0.0);                      \
  }                                                                     \
                                                                        \
  Teuchos::scan(*comm, Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);       \
  bool success1 = checkFadArrays(                                       \
    sums, sums2, std::string(#FAD)+"<"+#FAD+"> Scan Sum", out);         \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::scan(*comm, ffts, Teuchos::REDUCE_SUM, n, &x[0], &sums3[0]); \
  bool success2 = checkFadArrays(                                       \
    sums, sums3, std::string(#FAD)+"<"+#FAD+"> Scan Sum FTS", out);     \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_ScanMax ) {                       \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int rank = comm->getRank();                                           \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), maxs(n), maxs2(n), maxs3(n);         \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*(rank+1));                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*(rank+1);                           \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*(rank+1));                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*(rank+1);                           \
    maxs[i] = FadFadType(p2, f);                                        \
    for (int j=0; j<p2; j++)                                            \
      maxs[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs2[i] = FadFadType(p2, FadType(p1, 0.0));                        \
    for (int j=0; j<p2; j++)                                            \
      maxs2[i].fastAccessDx(j) = FadType(p1, 0.0);                      \
  }                                                                     \
                                                                        \
  Teuchos::scan(*comm, Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);       \
  bool success1 = checkFadArrays(                                       \
    maxs, maxs2, std::string(#FAD)+"<"+#FAD+"> Scan Max", out);         \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::scan(*comm, ffts, Teuchos::REDUCE_MAX, n, &x[0], &maxs3[0]); \
  bool success2 = checkFadArrays(                                       \
    maxs, maxs3, std::string(#FAD)+"<"+#FAD+"> Scan Max FTS", out);     \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_ScanMin ) {                       \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p1 = 5;                                                           \
  int p2 = 5;                                                           \
  int rank = comm->getRank();                                           \
  RCP< ValueTypeSerializer<int,FadType> > fts =                         \
    rcp(new ValueTypeSerializer<int,FadType>(                           \
          rcp(new ValueTypeSerializer<int,double>), p1));               \
  ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                    \
                                                                        \
  Teuchos::Array<FadFadType> x(n), mins(n), mins2(n), mins3(n);         \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1)*(rank+1));                                  \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1)*(rank+1);                           \
    x[i] = FadFadType(p2, f);                                           \
    for (int j=0; j<p2; j++) {                                          \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    FadType f(p1, 1.0*(i+1));                                           \
    for (int k=0; k<p1; k++)                                            \
      f.fastAccessDx(k) = 2.0*(i+1);                                    \
    mins[i] = FadFadType(p2, f);                                        \
    for (int j=0; j<p2; j++)                                            \
      mins[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins2[i] = FadFadType(p2, FadType(p1, 0.0));                        \
    for (int j=0; j<p2; j++)                                            \
      mins2[i].fastAccessDx(j) = FadType(p1, 0.0);                      \
  }                                                                     \
                                                                        \
  Teuchos::scan(*comm, Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);       \
  bool success1 = checkFadArrays(                                       \
    mins, mins2, std::string(#FAD)+"<"+#FAD+"> Scan Min", out);         \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  Teuchos::scan(*comm, ffts, Teuchos::REDUCE_MIN, n, &x[0], &mins3[0]); \
  bool success2 = checkFadArrays(                                       \
    mins, mins3, std::string(#FAD)+"<"+#FAD+"> Scan Min FTS", out);     \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( FAD##_Comm, FadFad_SendReceive ) {                   \
  typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int num_proc = comm->getSize();                                       \
  if (num_proc > 1) {                                                   \
    int rank = comm->getRank();                                         \
    int n = 7;                                                          \
    int p1 = 5;                                                         \
    int p2 = 5;                                                         \
    RCP< ValueTypeSerializer<int,FadType> > fts =                       \
      rcp(new ValueTypeSerializer<int,FadType>(                         \
            rcp(new ValueTypeSerializer<int,double>), p1));             \
    ValueTypeSerializer<int,FadFadType> ffts(fts, p2);                  \
                                                                        \
    Teuchos::Array<FadFadType> x(n), x2(n), x3(n);                      \
    for (int i=0; i<n; i++) {                                           \
      FadType f(p1, 1.0*(i+1));                                         \
      for (int k=0; k<p1; k++)                                          \
        f.fastAccessDx(k) = 2.0*(i+1)*(k+1);                            \
      x[i] = FadFadType(p2, f);                                         \
      for (int j=0; j<p2; j++)                                          \
        x[i].fastAccessDx(j) = f;                                       \
    }                                                                   \
    for (int i=0; i<n; i++) {                                           \
      x2[i] = FadFadType(p2, FadType(p1, 0.0));                         \
      for (int j=0; j<p2; j++)                                          \
        x2[i].fastAccessDx(j) = FadType(p1, 0.0);                       \
    }                                                                   \
    if (rank != 1) {                                                    \
      x2 = x;                                                           \
      x3 = x;                                                           \
    }                                                                   \
                                                                        \
    if (rank == 0) Teuchos::send(*comm, n, &x[0], 1);                   \
    if (rank == 1) Teuchos::receive(*comm, 0, n, &x2[0]);               \
    bool success1 = checkFadArrays(                                     \
      x, x2, std::string(#FAD)+"<"+#FAD+"> Send/Receive", out);         \
    success1 = checkResultOnAllProcs(*comm, out, success1);             \
                                                                        \
    if (rank == 0) Teuchos::send(*comm, ffts, n, &x[0], 1);             \
    if (rank == 1) Teuchos::receive(*comm, ffts, 0, n, &x3[0]);         \
    bool success2 = checkFadArrays(                                     \
      x, x3, std::string(#FAD)+"<"+#FAD+"> Send/Receive FTS", out);     \
    success2 = checkResultOnAllProcs(*comm, out, success2);             \
                                                                        \
    success = success1 && success2;                                     \
  }                                                                     \
  else                                                                  \
    success = true;                                                     \
}

#if defined(HAVE_SACADO_KOKKOS) && defined(HAVE_SACADO_TEUCHOSKOKKOSCOMM)

#include "Kokkos_Core.hpp"

#define FAD_KOKKOS_COMM_TESTS_DEV(FadType, FAD, Device)                 \
TEUCHOS_UNIT_TEST( FAD##_Comm_Kokkos_##Device, Fad_Broadcast ) {        \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  typedef Kokkos::View<FadType*,Device> ViewType;                       \
  typedef ViewType::HostMirror HostViewType;                            \
  ViewType x("x",n,p+1), x2("x2",n,p+1), x3("x3",n,p+1);                \
  HostViewType h_x = Kokkos::create_mirror_view(x);                     \
  HostViewType h_x2 = Kokkos::create_mirror_view(x2);                   \
  HostViewType h_x3 = Kokkos::create_mirror_view(x3);                   \
  for (int i=0; i<n; i++) {                                             \
    h_x[i] = FadType(p, rnd.number());                                  \
    for (int j=0; j<p; j++)                                             \
      h_x[i].fastAccessDx(j) = rnd.number();                            \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_x2[i] = FadType(p, 0.0);                                          \
  }                                                                     \
  Kokkos::deep_copy(x, h_x);                                            \
  Kokkos::deep_copy(x2, h_x2);                                          \
  if (comm->getRank() == 0) {                                           \
    x2 = x;                                                             \
    x3 = x;                                                             \
    h_x2 = h_x;                                                         \
    h_x3 = h_x;                                                         \
  }                                                                     \
                                                                        \
  /* The Teuchos MPI wrappers know nothing of CUDA nor CUDA-aware MPI*/ \
  /* so only do the communication on the host.  This probably makes  */ \
  /* the deep copy unnecessary.                                      */ \
  const bool accessible =                                               \
    Kokkos::Impl::MemorySpaceAccess<                                    \
      Kokkos::HostSpace,                                                \
      typename Device::memory_space >::accessible;                      \
  if (accessible) {                                                     \
    Teuchos::broadcast(*comm, 0, n, x2);                                \
    Kokkos::deep_copy(h_x2, x2);                                        \
  }                                                                     \
  else                                                                  \
    Teuchos::broadcast(*comm, 0, n, h_x2);                              \
  bool success1 = checkFadArrays(                                       \
    h_x, h_x2, std::string(#FAD)+" Broadcast", out);                    \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  if (accessible) {                                                     \
    Teuchos::broadcast(*comm, fts, 0, n, x3);                           \
    Kokkos::deep_copy(h_x3, x3);                                        \
  }                                                                     \
  else                                                                  \
    Teuchos::broadcast(*comm, fts, 0, n, h_x3);                         \
  bool success2 = checkFadArrays(                                       \
    h_x, h_x3, std::string(#FAD)+" Broadcast FTS", out);                \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
                                                                        \
  success = success1 && success2;                                       \
}                                                                       \
TEUCHOS_UNIT_TEST( FAD##_Comm_Kokkos_##Device, Fad_SumAll ) {           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int num_proc = comm->getSize();                                       \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  typedef Kokkos::View<FadType*,Device> ViewType;                       \
  typedef ViewType::HostMirror HostViewType;                            \
  ViewType x("x",n,p+1), sums("sums",n,p+1),                            \
    sums2("sums2",n,p+1), sums3("sums3",n,p+1);                         \
  HostViewType h_x = Kokkos::create_mirror_view(x);                     \
  HostViewType h_sums = Kokkos::create_mirror_view(sums);               \
  HostViewType h_sums2 = Kokkos::create_mirror_view(sums2);             \
  HostViewType h_sums3 = Kokkos::create_mirror_view(sums3);             \
  for (int i=0; i<n; i++) {                                             \
    h_x[i] = FadType(p, 1.0*(i+1));                                     \
    for (int j=0; j<p; j++)                                             \
      h_x[i].fastAccessDx(j) = 2.0*(i+1);                               \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_sums[i] = FadType(p, 1.0*(i+1)*num_proc);                         \
    for (int j=0; j<p; j++)                                             \
      h_sums[i].fastAccessDx(j) = 2.0*(i+1)*num_proc;                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_sums2[i] = FadType(p, 0.0);                                       \
  }                                                                     \
  Kokkos::deep_copy(x, h_x);                                            \
  Kokkos::deep_copy(sums, h_sums);                                      \
  Kokkos::deep_copy(sums2, h_sums2);                                    \
                                                                        \
  /* The Teuchos MPI wrappers know nothing of CUDA nor CUDA-aware MPI*/ \
  /* so only do the communication on the host.  This probably makes  */ \
  /* the deep copy unnecessary.                                      */ \
  const bool accessible =                                               \
    Kokkos::Impl::MemorySpaceAccess<                                    \
      Kokkos::HostSpace,                                                \
      typename Device::memory_space >::accessible;                      \
  if (accessible) {                                                     \
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, n, x, sums2);        \
    Kokkos::deep_copy(h_sums2, sums2);                                  \
  }                                                                     \
  else                                                                  \
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, n, h_x, h_sums2);    \
  bool success1 = checkFadArrays(                                       \
    h_sums, h_sums2, std::string(#FAD)+" Sum All", out);                \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  if (accessible) {                                                     \
    Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_SUM, n, x, sums3);   \
    Kokkos::deep_copy(h_sums3, sums3);                                  \
  }                                                                     \
  else                                                                  \
    Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_SUM, n, h_x, h_sums3); \
  bool success2 = checkFadArrays(                                       \
    h_sums, h_sums3, std::string(#FAD)+" Sum All FTS", out);            \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
  success = success1 && success2;                                       \
                                                                        \
}                                                                       \
TEUCHOS_UNIT_TEST( FAD##_Comm_Kokkos_##Device, Fad_MaxAll ) {           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  int num_proc = comm->getSize();                                       \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  typedef Kokkos::View<FadType*,Device> ViewType;                       \
  typedef ViewType::HostMirror HostViewType;                            \
  ViewType x("x",n,p+1), maxs("maxs",n,p+1),                            \
    maxs2("maxs2",n,p+1), maxs3("maxs3",n,p+1);                         \
  HostViewType h_x = Kokkos::create_mirror_view(x);                     \
  HostViewType h_maxs = Kokkos::create_mirror_view(maxs);               \
  HostViewType h_maxs2 = Kokkos::create_mirror_view(maxs2);             \
  HostViewType h_maxs3 = Kokkos::create_mirror_view(maxs3);             \
  for (int i=0; i<n; i++) {                                             \
    h_x[i] = FadType(p, 1.0*(i+1)*(rank+1));                            \
    for (int j=0; j<p; j++)                                             \
      h_x[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_maxs[i] = FadType(p, 1.0*(i+1)*num_proc);                         \
    for (int j=0; j<p; j++)                                             \
      h_maxs[i].fastAccessDx(j) = 2.0*(i+1)*num_proc;                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_maxs2[i] = FadType(p, 0.0);                                       \
  }                                                                     \
  Kokkos::deep_copy(x, h_x);                                            \
  Kokkos::deep_copy(maxs, h_maxs);                                      \
  Kokkos::deep_copy(maxs2, h_maxs2);                                    \
                                                                        \
  /* The Teuchos MPI wrappers know nothing of CUDA nor CUDA-aware MPI*/ \
  /* so only do the communication on the host.  This probably makes  */ \
  /* the deep copy unnecessary.                                      */ \
  const bool accessible =                                               \
    Kokkos::Impl::MemorySpaceAccess<                                    \
      Kokkos::HostSpace,                                                \
      typename Device::memory_space >::accessible;                      \
  if (accessible) {                                                     \
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, n, x, maxs2);        \
    Kokkos::deep_copy(h_maxs2, maxs2);                                  \
  }                                                                     \
  else                                                                  \
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, n, h_x, h_maxs2);    \
  bool success1 = checkFadArrays(                                       \
    h_maxs, h_maxs2, std::string(#FAD)+" Max All", out);                \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  if (accessible) {                                                     \
    Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_MAX, n, x, maxs3);   \
    Kokkos::deep_copy(h_maxs3, maxs3);                                  \
  }                                                                     \
  else                                                                  \
    Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_MAX, n, h_x, h_maxs3); \
  bool success2 = checkFadArrays(                                       \
    h_maxs, h_maxs3, std::string(#FAD)+" Max All FTS", out);            \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
  success = success1 && success2;                                       \
                                                                        \
}                                                                       \
TEUCHOS_UNIT_TEST( FAD##_Comm_Kokkos_##Device, Fad_MinAll ) {           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  ValueTypeSerializer<int,FadType> fts(                                 \
    rcp(new ValueTypeSerializer<int,double>), p);                       \
                                                                        \
  typedef Kokkos::View<FadType*,Device> ViewType;                       \
  typedef ViewType::HostMirror HostViewType;                            \
  ViewType x("x",n,p+1), mins("mins",n,p+1),                            \
    mins2("mins2",n,p+1), mins3("mins3",n,p+1);                         \
  HostViewType h_x = Kokkos::create_mirror_view(x);                     \
  HostViewType h_mins = Kokkos::create_mirror_view(mins);               \
  HostViewType h_mins2 = Kokkos::create_mirror_view(mins2);             \
  HostViewType h_mins3 = Kokkos::create_mirror_view(mins3);             \
  for (int i=0; i<n; i++) {                                             \
    h_x[i] = FadType(p, 1.0*(i+1)*(rank+1));                            \
    for (int j=0; j<p; j++)                                             \
      h_x[i].fastAccessDx(j) = 2.0*(i+1)*(rank+1);                      \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_mins[i] = FadType(p, 1.0*(i+1));                                  \
    for (int j=0; j<p; j++)                                             \
      h_mins[i].fastAccessDx(j) = 2.0*(i+1);                            \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    h_mins2[i] = FadType(p, 0.0);                                       \
  }                                                                     \
  Kokkos::deep_copy(x, h_x);                                            \
  Kokkos::deep_copy(mins, h_mins);                                      \
  Kokkos::deep_copy(mins2, h_mins2);                                    \
                                                                        \
  /* The Teuchos MPI wrappers know nothing of CUDA nor CUDA-aware MPI*/ \
  /* so only do the communication on the host.  This probably makes  */ \
  /* the deep copy unnecessary.                                      */ \
  const bool accessible =                                               \
    Kokkos::Impl::MemorySpaceAccess<                                    \
      Kokkos::HostSpace,                                                \
      typename Device::memory_space >::accessible;                      \
  if (accessible) {                                                     \
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, n, x, mins2);        \
    Kokkos::deep_copy(h_mins2, mins2);                                  \
  }                                                                     \
  else                                                                  \
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, n, h_x, h_mins2);    \
  bool success1 = checkFadArrays(                                       \
    h_mins, h_mins2, std::string(#FAD)+" Min All", out);                \
  success1 = checkResultOnAllProcs(*comm, out, success1);               \
                                                                        \
  if (accessible) {                                                     \
    Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_MIN, n, x, mins3);   \
    Kokkos::deep_copy(h_mins3, mins3);                                  \
  }                                                                     \
  else                                                                  \
    Teuchos::reduceAll(*comm, fts, Teuchos::REDUCE_MIN, n, h_x, h_mins3); \
  bool success2 = checkFadArrays(                                       \
    h_mins, h_mins3, std::string(#FAD)+" Min All FTS", out);            \
  success2 = checkResultOnAllProcs(*comm, out, success2);               \
  success = success1 && success2;                                       \
                                                                        \
}

#ifdef KOKKOS_ENABLE_OPENMP
#define FAD_KOKKOS_COMM_TESTS_OPENMP(FadType, FAD)                      \
  using Kokkos::OpenMP;                                                 \
  FAD_KOKKOS_COMM_TESTS_DEV(FadType, FAD, OpenMP)
#else
#define FAD_KOKKOS_COMM_TESTS_OPENMP(FadType, FAD)
#endif

#ifdef KOKKOS_ENABLE_THREADS
#define FAD_KOKKOS_COMM_TESTS_THREADS(FadType, FAD)                      \
  using Kokkos::Threads;                                                 \
  FAD_KOKKOS_COMM_TESTS_DEV(FadType, FAD, Threads)
#else
#define FAD_KOKKOS_COMM_TESTS_THREADS(FadType, FAD)
#endif

#ifdef KOKKOS_ENABLE_CUDA
#define FAD_KOKKOS_COMM_TESTS_CUDA(FadType, FAD)                         \
  using Kokkos::Cuda;                                                    \
  FAD_KOKKOS_COMM_TESTS_DEV(FadType, FAD, Cuda)
#else
#define FAD_KOKKOS_COMM_TESTS_CUDA(FadType, FAD)
#endif

#ifdef KOKKOS_ENABLE_HIP
#define FAD_KOKKOS_COMM_TESTS_HIP(FadType, FAD)                          \
  using Kokkos::HIP;					\
  FAD_KOKKOS_COMM_TESTS_DEV(FadType, FAD, HIP)
#else
#define FAD_KOKKOS_COMM_TESTS_HIP(FadType, FAD)
#endif

#ifdef KOKKOS_ENABLE_SERIAL
#define FAD_KOKKOS_COMM_TESTS_SERIAL(FadType, FAD)                      \
  using Kokkos::Serial;                                                 \
  FAD_KOKKOS_COMM_TESTS_DEV(FadType, FAD, Serial)
#else
#define FAD_KOKKOS_COMM_TESTS_SERIAL(FadType, FAD)
#endif

#define FAD_KOKKOS_COMM_TESTS(FadType, FAD)                             \
  FAD_KOKKOS_COMM_TESTS_OPENMP(FadType, FAD)                            \
  FAD_KOKKOS_COMM_TESTS_THREADS(FadType, FAD)                           \
  FAD_KOKKOS_COMM_TESTS_CUDA(FadType, FAD)                              \
  FAD_KOKKOS_COMM_TESTS_SERIAL(FadType, FAD)

#else

#define FAD_KOKKOS_COMM_TESTS(FadType, FAD)

#endif

#define FAD_COMM_TESTS(FadType, FAD)        \
  FAD_BASE_COMM_TESTS(FadType, FAD)

// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"

#include "Sacado.hpp"
#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

#include <Kokkos_Core.hpp>

//
// Currently this doesn't test:
//   * the device
//   * threaded storage (needs the device)
//   * strided storage with non-trivial stride
//

using Teuchos::RCP;
using Teuchos::rcp;

// Common setup for unit tests
template <typename VecType, typename FadType>
struct UnitTestSetup {
  int sz;

  typedef Teuchos::ValueTypeSerializer<int, VecType> VecSerializerT;
  RCP<VecSerializerT> vec_serializer;

  typedef typename Sacado::mpl::apply<FadType,VecType>::type FadVecType;
  typedef Teuchos::ValueTypeSerializer<int, FadVecType> FadVecSerializerT;
  RCP<FadVecSerializerT> fad_vec_serializer;

  UnitTestSetup() {
    sz = 8;

    // Serializers
    vec_serializer =
      rcp(new VecSerializerT(
            rcp(new Teuchos::ValueTypeSerializer<int,double>()), sz));
    fad_vec_serializer = rcp(new FadVecSerializerT(vec_serializer, 5));
  }
};

template <typename VecType>
bool checkVecArrays(const Teuchos::Array<VecType>& x,
                    const Teuchos::Array<VecType>& x2,
                    const std::string& tag,
                    Teuchos::FancyOStream& out) {

  // Check sizes match
  bool success = (x.size() == x2.size());
  out << tag << " Vec array size test";
  if (success)
    out << " passed";
  else
    out << " failed";
  out << ":  \n\tExpected:  " << x.size() << ", \n\tGot:       " << x2.size()
      << "." << std::endl;

  // Check Fads match
  for (int i=0; i<x.size(); i++) {
    bool success2 = Sacado::IsEqual<VecType>::eval(x[i], x2[i]);
    out << tag << " Vec array comparison test " << i;
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

#define VEC_COMM_TESTS(VecType, FadType, Vec, FAD)                      \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_Broadcast ) {                        \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  Teuchos::Array<VecType> x(n), x2(n);                                  \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = rnd.number();                           \
  }                                                                     \
  if (comm->getRank() == 0)                                             \
    x2 = x;                                                             \
  Teuchos::broadcast(*comm, *setup.vec_serializer, 0, n, &x2[0]);       \
  success = checkVecArrays(x, x2, std::string(#Vec)+" Broadcast", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_GatherAll ) {                        \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int size = comm->getSize();                                           \
  int rank = comm->getRank();                                           \
  int N = n*size;                                                       \
  Teuchos::Array<VecType> x(n), x2(N), x3(N);                           \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = (rank+1)*(i+1)*(j+1);                   \
  }                                                                     \
  for (int j=0; j<size; j++) {                                          \
    for (int i=0; i<n; i++) {                                           \
      x3[n*j+i] = VecType(setup.sz, 0.0);                               \
      for (int k=0; k<setup.sz; k++)                                    \
        x3[n*j+i].fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);               \
    }                                                                   \
  }                                                                     \
  Teuchos::gatherAll(*comm, *setup.vec_serializer,                      \
                     n, &x[0], N, &x2[0]);                              \
  success = checkVecArrays(x3, x2, std::string(#Vec)+" Gather All", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_SumAll ) {                           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int num_proc = comm->getSize();                                       \
                                                                        \
  Teuchos::Array<VecType> x(n), sums(n), sums2(n);                      \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = 2.0*(i+1);                              \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums[i] = VecType(setup.sz, 0.0);                                   \
    for (int j=0; j<setup.sz; j++)                                      \
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;                  \
  }                                                                     \
  Teuchos::reduceAll(*comm, *setup.vec_serializer,                      \
                     Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);         \
  success = checkVecArrays(sums, sums2,                                 \
                           std::string(#Vec)+" Sum All", out);          \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_MaxAll ) {                           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int rank = comm->getRank();                                           \
  int num_proc = comm->getSize();                                       \
                                                                        \
  Teuchos::Array<VecType> x(n), maxs(n), maxs2(n);                      \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs[i] = VecType(setup.sz, 0.0);                                   \
    for (int j=0; j<setup.sz; j++)                                      \
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;                  \
  }                                                                     \
  Teuchos::reduceAll(*comm, *setup.vec_serializer,                      \
                     Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);         \
  success = checkVecArrays(maxs, maxs2,                                 \
                           std::string(#Vec)+" Max All", out);          \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_MinAll ) {                           \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<VecType> x(n), mins(n), mins2(n);                      \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins[i] = VecType(setup.sz, 0.0);                                   \
    for (int j=0; j<setup.sz; j++)                                      \
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);                           \
  }                                                                     \
  Teuchos::reduceAll(*comm, *setup.vec_serializer,                      \
                     Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);         \
  success = checkVecArrays(mins, mins2,                                 \
                           std::string(#Vec)+" Min All", out);          \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_ScanSum ) {                          \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<VecType> x(n), sums(n), sums2(n);                      \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = 2.0*(i+1);                              \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    sums[i] = VecType(setup.sz, 0.0);                                   \
    for (int j=0; j<setup.sz; j++)                                      \
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);                  \
  }                                                                     \
  Teuchos::scan(*comm, *setup.vec_serializer,                           \
                Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);              \
  success = checkVecArrays(sums, sums2,                                 \
                           std::string(#Vec)+" Scan Sum", out);         \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_ScanMax ) {                          \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<VecType> x(n), maxs(n), maxs2(n);                      \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    maxs[i] = VecType(setup.sz, 0.0);                                   \
    for (int j=0; j<setup.sz; j++)                                      \
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);                  \
  }                                                                     \
  Teuchos::scan(*comm, *setup.vec_serializer,                           \
                Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);              \
  success = checkVecArrays(maxs, maxs2,                                 \
                           std::string(#Vec)+" Scan Max", out);         \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_ScanMin ) {                          \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<VecType> x(n), mins(n), mins2(n);                      \
  for (int i=0; i<n; i++) {                                             \
    x[i] = VecType(setup.sz, 0.0);                                      \
    for (int j=0; j<setup.sz; j++)                                      \
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);                     \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    mins[i] = VecType(setup.sz, 0.0);                                   \
    for (int j=0; j<setup.sz; j++)                                      \
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);                           \
  }                                                                     \
  Teuchos::scan(*comm, *setup.vec_serializer,                           \
                Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);              \
  success = checkVecArrays(mins, mins2,                                 \
                           std::string(#Vec)+" Scan Min", out);         \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, Vec_SendReceive ) {                      \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int num_proc = comm->getSize();                                       \
  if (num_proc > 1) {                                                   \
    int rank = comm->getRank();                                         \
    int n = 7;                                                          \
    Teuchos::Array<VecType> x(n), x2(n);                                \
    for (int i=0; i<n; i++) {                                           \
      x[i] = VecType(setup.sz, 0.0);                                    \
      for (int j=0; j<setup.sz; j++)                                    \
        x[i].fastAccessCoeff(j) = 2.0*(i+1)*(j+1);                      \
    }                                                                   \
    if (rank != 1)                                                      \
      x2 = x;                                                           \
    if (rank == 0) Teuchos::send(*comm, *setup.vec_serializer,          \
                                 n, &x[0], 1);                          \
    if (rank == 1) Teuchos::receive(*comm, *setup.vec_serializer,       \
                                    0, n, &x2[0]);                      \
    success = checkVecArrays(x, x2,                                     \
                             std::string(#Vec)+" Send/Receive", out);   \
    success = checkResultOnAllProcs(*comm, out, success);               \
  }                                                                     \
  else                                                                  \
    success = true;                                                     \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_Broadcast ) {                     \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  Teuchos::Array<FadVecType> x(n), x2(n);                               \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = rnd.number();                              \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      VecType g(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        g.fastAccessCoeff(k) = rnd.number();                            \
      x[i].fastAccessDx(j) = g;                                         \
    }                                                                   \
  }                                                                     \
  if (comm->getRank() == 0)                                             \
    x2 = x;                                                             \
  Teuchos::broadcast(*comm, *setup.fad_vec_serializer, 0, n, &x2[0]);   \
  success = checkVecArrays(x, x2,                                       \
                           std::string(#FAD)+"<"+#Vec+"> Broadcast", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_GatherAll ) {                     \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int size = comm->getSize();                                           \
  int rank = comm->getRank();                                           \
  int N = n*size;                                                       \
  Teuchos::Array<FadVecType> x(n), x2(N), x3(N);                        \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = (rank+1)*(i+1)*(k+1);                      \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int j=0; j<size; j++) {                                          \
    for (int i=0; i<n; i++) {                                           \
      VecType f(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        f.fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);                       \
      x3[n*j+i] = FadVecType(p, f);                                     \
      for (int k=0; k<p; k++)                                           \
        x3[n*j+i].fastAccessDx(k) = f;                                  \
    }                                                                   \
  }                                                                     \
  Teuchos::gatherAll(*comm, *setup.fad_vec_serializer,                  \
                     n, &x[0], N, &x2[0]);                              \
  success = checkVecArrays(x3, x2,                                      \
                           std::string(#FAD)+"<"+#Vec+">  Gather All", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_SumAll ) {                        \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int num_proc = comm->getSize();                                       \
                                                                        \
  Teuchos::Array<FadVecType> x(n), sums(n), sums2(n);                   \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1);                                 \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      VecType g(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        g.fastAccessCoeff(k) = 2.0*(i+1);                               \
      x[i].fastAccessDx(j) = g;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;                        \
    sums[i] = FadVecType(p, f);                                         \
    for (int j=0; j<p; j++) {                                           \
      VecType g(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        g.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;                      \
      sums[i].fastAccessDx(j) = g;                                      \
    }                                                                   \
  }                                                                     \
  Teuchos::reduceAll(*comm, *setup.fad_vec_serializer,                  \
                     Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);         \
  success = checkVecArrays(sums, sums2,                                 \
                           std::string(#FAD)+"<"+#Vec+"> Sum All", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_MaxAll ) {                        \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 8;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
  int num_proc = comm->getSize();                                       \
                                                                        \
  Teuchos::Array<FadVecType> x(n), maxs(n), maxs2(n);                   \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);                        \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;                        \
    maxs[i] = FadVecType(p, f);                                         \
    for (int j=0; j<p; j++)                                             \
      maxs[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  Teuchos::reduceAll(*comm, *setup.fad_vec_serializer,                  \
                     Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);         \
  success = checkVecArrays(maxs, maxs2,                                 \
                           std::string(#FAD)+"<"+#Vec+"> Max All", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_MinAll ) {                        \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 8;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<FadVecType> x(n), mins(n), mins2(n);                   \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);                        \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1);                                 \
    mins[i] = FadVecType(p, f);                                         \
    for (int j=0; j<p; j++)                                             \
      mins[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  Teuchos::reduceAll(*comm, *setup.fad_vec_serializer,                  \
                     Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);         \
  success = checkVecArrays(mins, mins2,                                 \
                           std::string(#FAD)+"<"+#Vec+"> Min All", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_ScanSum ) {                       \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<FadVecType> x(n), sums(n), sums2(n);                   \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1);                                 \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);                        \
    sums[i] = FadVecType(p, f);                                         \
    for (int j=0; j<p; j++)                                             \
      sums[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  Teuchos::scan(*comm, *setup.fad_vec_serializer,                       \
                Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);              \
  success = checkVecArrays(sums, sums2,                                 \
                           std::string(#FAD)+"<"+#Vec+"> Scan Sum", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_ScanMax ) {                       \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<FadVecType> x(n), maxs(n), maxs2(n);                   \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);                        \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);                        \
    maxs[i] = FadVecType(p, f);                                         \
    for (int j=0; j<p; j++)                                             \
      maxs[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  Teuchos::scan(*comm, *setup.fad_vec_serializer,                       \
                Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);              \
  success = checkVecArrays(maxs, maxs2,                                 \
                           std::string(#FAD)+"<"+#Vec+"> Scan Max", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_ScanMin ) {                       \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int n = 7;                                                            \
  int p = 5;                                                            \
  int rank = comm->getRank();                                           \
                                                                        \
  Teuchos::Array<FadVecType> x(n), mins(n), mins2(n);                   \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);                        \
    x[i] = FadVecType(p, f);                                            \
    for (int j=0; j<p; j++) {                                           \
      x[i].fastAccessDx(j) = f;                                         \
    }                                                                   \
  }                                                                     \
  for (int i=0; i<n; i++) {                                             \
    VecType f(setup.sz, 0.0);                                           \
    for (int k=0; k<setup.sz; k++)                                      \
      f.fastAccessCoeff(k) = 2.0*(i+1);                                 \
    mins[i] = FadVecType(p, f);                                         \
    for (int j=0; j<p; j++)                                             \
      mins[i].fastAccessDx(j) = f;                                      \
  }                                                                     \
  Teuchos::scan(*comm, *setup.fad_vec_serializer,                       \
                Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);              \
  success = checkVecArrays(mins, mins2,                                 \
                           std::string(#FAD)+"<"+#Vec+"> Scan Min", out); \
  success = checkResultOnAllProcs(*comm, out, success);                 \
}                                                                       \
                                                                        \
TEUCHOS_UNIT_TEST( Vec##_Comm, FadVec_SendReceive ) {                   \
  typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;         \
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >                           \
    comm = Teuchos::DefaultComm<Ordinal>::getComm();                    \
                                                                        \
  int num_proc = comm->getSize();                                       \
  if (num_proc > 1) {                                                   \
    int rank = comm->getRank();                                         \
    int n = 7;                                                          \
    int p = 5;                                                          \
    Teuchos::Array<FadVecType> x(n), x2(n);                             \
    for (int i=0; i<n; i++) {                                           \
      VecType f(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        f.fastAccessCoeff(k) = 2.0*(i+1)*(k+1);                         \
      x[i] = FadVecType(p, f);                                          \
      for (int j=0; j<p; j++)                                           \
        x[i].fastAccessDx(j) = f;                                       \
    }                                                                   \
    if (rank != 1)                                                      \
      x2 = x;                                                           \
    if (rank == 0) Teuchos::send(*comm, *setup.fad_vec_serializer,      \
                                 n, &x[0], 1);                          \
    if (rank == 1) Teuchos::receive(*comm, *setup.fad_vec_serializer,   \
                                    0, n, &x2[0]);                      \
    success = checkVecArrays(x, x2,                                     \
                             std::string(#FAD)+"<"+#Vec+"> Send/Receive", out); \
    success = checkResultOnAllProcs(*comm, out, success);               \
  }                                                                     \
  else                                                                  \
    success = true;                                                     \
}

namespace DynamicVecTest {
  Sacado::Random<double> rnd;
  typedef int Ordinal;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef Stokhos::DynamicStorage<int,double,execution_space> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_COMM_TESTS(vec_type, fad_type, DynamicVector, DFad)
}

namespace DynamicStridedVecTest {
  Sacado::Random<double> rnd;
  typedef int Ordinal;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef Stokhos::DynamicStridedStorage<int,double,execution_space> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_COMM_TESTS(vec_type, fad_type, DynamicStridedVector, DFad)
}

namespace StaticVecTest {
  Sacado::Random<double> rnd;
  typedef int Ordinal;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef Stokhos::StaticStorage<int,double,8,execution_space> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_COMM_TESTS(vec_type, fad_type, StaticVector, DFad)
}

namespace StaticFixedVecTest {
  Sacado::Random<double> rnd;
  typedef int Ordinal;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef Stokhos::StaticFixedStorage<int,double,8,execution_space> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_COMM_TESTS(vec_type, fad_type, StaticFixedVector, DFad)
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}

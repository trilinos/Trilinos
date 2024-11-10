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
#include "Stokhos.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

// Common setup for unit tests
template <typename PCEType, typename FadType>
struct UnitTestSetup {
  RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis;
  RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;

  typedef typename PCEType::cijk_type kokkos_cijk_type;
  kokkos_cijk_type kokkos_cijk;

  typedef Teuchos::ValueTypeSerializer<int, PCEType> PCESerializerT;
  RCP<PCESerializerT> pce_serializer;

  typedef typename Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  typedef Teuchos::ValueTypeSerializer<int, FadPCEType> FadPCESerializerT;
  RCP<FadPCESerializerT> fad_pce_serializer;
  int sz;

  UnitTestSetup() {
    const int d = 2;
    const int p = 7;

    // Create product basis
    Teuchos::Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++)
      bases[i] =
        rcp(new Stokhos::LegendreBasis<int,double>(p));
    basis =
      rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Triple product tensor
    typedef typename PCEType::execution_space execution_space;
    Cijk = basis->computeTripleProductTensor();
    kokkos_cijk =
      Stokhos::create_product_tensor<execution_space>(*basis, *Cijk);

    // Serializers
    pce_serializer =
      rcp(new PCESerializerT(
            kokkos_cijk,
            rcp(new Teuchos::ValueTypeSerializer<int,double>())));
    fad_pce_serializer = rcp(new FadPCESerializerT(pce_serializer, 5));

    sz = basis->size();
  }
};

template <typename PCEType>
bool checkPCEArrays(const Teuchos::Array<PCEType>& x,
                    const Teuchos::Array<PCEType>& x2,
                    const std::string& tag,
                    Teuchos::FancyOStream& out) {

  // Check sizes match
  bool success = (x.size() == x2.size());
  out << tag << " PCE array size test";
  if (success)
    out << " passed";
  else
    out << " failed";
  out << ":  \n\tExpected:  " << x.size() << ", \n\tGot:       " << x2.size()
      << "." << std::endl;

  // Check Fads match
  for (int i=0; i<x.size(); i++) {
    bool success2 = Sacado::IsEqual<PCEType>::eval(x[i], x2[i]);
    out << tag << " PCE array comparison test " << i;
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

typedef int Ordinal;
typedef Sacado::Fad::DFad<double> FadType;
typedef Kokkos::DefaultExecutionSpace execution_space;
typedef Stokhos::DynamicStorage<int,double,execution_space> storage_type;
typedef Sacado::UQ::PCE<storage_type> PCEType;

Sacado::Random<double> rnd;

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_Broadcast ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  Teuchos::Array<PCEType> x(n), x2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = rnd.number();
  }
  if (comm->getRank() == 0)
    x2 = x;
  Teuchos::broadcast(*comm, *setup.pce_serializer, 0, n, &x2[0]);
  success = checkPCEArrays(x, x2, std::string("UQ::PCE")+" Broadcast", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_GatherAll ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int size = comm->getSize();
  int rank = comm->getRank();
  int N = n*size;
  Teuchos::Array<PCEType> x(n), x2(N), x3(N);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = (rank+1)*(i+1)*(j+1);
  }
  for (int j=0; j<size; j++) {
    for (int i=0; i<n; i++) {
      x3[n*j+i].reset(setup.kokkos_cijk);
      for (int k=0; k<setup.sz; k++)
        x3[n*j+i].fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);
    }
  }
  Teuchos::gatherAll(*comm, *setup.pce_serializer,
                     n, &x[0], N, &x2[0]);
  success = checkPCEArrays(x3, x2, std::string("UQ::PCE")+" Gather All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_SumAll ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int num_proc = comm->getSize();

  Teuchos::Array<PCEType> x(n), sums(n), sums2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = 2.0*(i+1);
  }
  for (int i=0; i<n; i++) {
    sums[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;
  }
  Teuchos::reduceAll(*comm, *setup.pce_serializer,
                     Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);
  success = checkPCEArrays(sums, sums2,
                           std::string("UQ::PCE")+" Sum All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_MaxAll ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int rank = comm->getRank();
  int num_proc = comm->getSize();

  Teuchos::Array<PCEType> x(n), maxs(n), maxs2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);
  }
  for (int i=0; i<n; i++) {
    maxs[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;
  }
  Teuchos::reduceAll(*comm, *setup.pce_serializer,
                     Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);
  success = checkPCEArrays(maxs, maxs2,
                           std::string("UQ::PCE")+" Max All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_MinAll ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int rank = comm->getRank();

  Teuchos::Array<PCEType> x(n), mins(n), mins2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);
  }
  for (int i=0; i<n; i++) {
    mins[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);
  }
  Teuchos::reduceAll(*comm, *setup.pce_serializer,
                     Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);
  success = checkPCEArrays(mins, mins2,
                           std::string("UQ::PCE")+" Min All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_ScanSum ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int rank = comm->getRank();

  Teuchos::Array<PCEType> x(n), sums(n), sums2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = 2.0*(i+1);
  }
  for (int i=0; i<n; i++) {
    sums[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);
  }
  Teuchos::scan(*comm, *setup.pce_serializer,
                Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);
  success = checkPCEArrays(sums, sums2,
                           std::string("UQ::PCE")+" Scan Sum", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_ScanMax ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int rank = comm->getRank();

  Teuchos::Array<PCEType> x(n), maxs(n), maxs2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);
  }
  for (int i=0; i<n; i++) {
    maxs[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);
  }
  Teuchos::scan(*comm, *setup.pce_serializer,
                Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);
  success = checkPCEArrays(maxs, maxs2,
                           std::string("UQ::PCE")+" Scan Max", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_ScanMin ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int rank = comm->getRank();

  Teuchos::Array<PCEType> x(n), mins(n), mins2(n);
  for (int i=0; i<n; i++) {
    x[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);
  }
  for (int i=0; i<n; i++) {
    mins[i].reset(setup.kokkos_cijk);
    for (int j=0; j<setup.sz; j++)
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);
  }
  Teuchos::scan(*comm, *setup.pce_serializer,
                Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);
  success = checkPCEArrays(mins, mins2,
                           std::string("UQ::PCE")+" Scan Min", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, PCE_SendReceive ) {
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int num_proc = comm->getSize();
  if (num_proc > 1) {
    int rank = comm->getRank();
    int n = 7;
    Teuchos::Array<PCEType> x(n), x2(n);
    for (int i=0; i<n; i++) {
      x[i].reset(setup.kokkos_cijk);
      for (int j=0; j<setup.sz; j++)
        x[i].fastAccessCoeff(j) = 2.0*(i+1)*(j+1);
    }
    if (rank != 1)
      x2 = x;
    if (rank == 0) Teuchos::send(*comm, *setup.pce_serializer,
                                 n, &x[0], 1);
    if (rank == 1) Teuchos::receive(*comm, *setup.pce_serializer,
                                    0, n, &x2[0]);
    success = checkPCEArrays(x, x2,
                             std::string("UQ::PCE")+" Send/Receive", out);
    success = checkResultOnAllProcs(*comm, out, success);
  }
  else
    success = true;
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_Broadcast ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int p = 5;
  Teuchos::Array<FadPCEType> x(n), x2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = rnd.number();
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      PCEType g(setup.kokkos_cijk);
      for (int k=0; k<setup.sz; k++)
        g.fastAccessCoeff(k) = rnd.number();
      x[i].fastAccessDx(j) = g;
    }
  }
  if (comm->getRank() == 0)
    x2 = x;
  Teuchos::broadcast(*comm, *setup.fad_pce_serializer, 0, n, &x2[0]);
  success = checkPCEArrays(x, x2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Broadcast", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_GatherAll ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int p = 5;
  int size = comm->getSize();
  int rank = comm->getRank();
  int N = n*size;
  Teuchos::Array<FadPCEType> x(n), x2(N), x3(N);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = (rank+1)*(i+1)*(k+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      x[i].fastAccessDx(j) = f;
    }
  }
  for (int j=0; j<size; j++) {
    for (int i=0; i<n; i++) {
      PCEType f(setup.kokkos_cijk);
      for (int k=0; k<setup.sz; k++)
        f.fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);
      x3[n*j+i] = FadPCEType(p, f);
      for (int k=0; k<p; k++)
        x3[n*j+i].fastAccessDx(k) = f;
    }
  }
  Teuchos::gatherAll(*comm, *setup.fad_pce_serializer,
                     n, &x[0], N, &x2[0]);
  success = checkPCEArrays(x3, x2,
                           std::string("DFad")+"<"+"UQ::PCE"+">  Gather All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_SumAll ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int p = 5;
  int num_proc = comm->getSize();

  Teuchos::Array<FadPCEType> x(n), sums(n), sums2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      PCEType g(setup.kokkos_cijk);
      for (int k=0; k<setup.sz; k++)
        g.fastAccessCoeff(k) = 2.0*(i+1);
      x[i].fastAccessDx(j) = g;
    }
  }
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;
    sums[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      PCEType g(setup.kokkos_cijk);
      for (int k=0; k<setup.sz; k++)
        g.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;
      sums[i].fastAccessDx(j) = g;
    }
  }
  Teuchos::reduceAll(*comm, *setup.fad_pce_serializer,
                     Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);
  success = checkPCEArrays(sums, sums2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Sum All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_MaxAll ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 8;
  int p = 5;
  int rank = comm->getRank();
  int num_proc = comm->getSize();

  Teuchos::Array<FadPCEType> x(n), maxs(n), maxs2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      x[i].fastAccessDx(j) = f;
    }
  }
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;
    maxs[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++)
      maxs[i].fastAccessDx(j) = f;
  }
  Teuchos::reduceAll(*comm, *setup.fad_pce_serializer,
                     Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);
  success = checkPCEArrays(maxs, maxs2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Max All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_MinAll ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 8;
  int p = 5;
  int rank = comm->getRank();

  Teuchos::Array<FadPCEType> x(n), mins(n), mins2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      x[i].fastAccessDx(j) = f;
    }
  }
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1);
    mins[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++)
      mins[i].fastAccessDx(j) = f;
  }
  Teuchos::reduceAll(*comm, *setup.fad_pce_serializer,
                     Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);
  success = checkPCEArrays(mins, mins2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Min All", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_ScanSum ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int p = 5;
  int rank = comm->getRank();

  Teuchos::Array<FadPCEType> x(n), sums(n), sums2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      x[i].fastAccessDx(j) = f;
    }
  }
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);
    sums[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++)
      sums[i].fastAccessDx(j) = f;
  }
  Teuchos::scan(*comm, *setup.fad_pce_serializer,
                Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);
  success = checkPCEArrays(sums, sums2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Scan Sum", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_ScanMax ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int p = 5;
  int rank = comm->getRank();

  Teuchos::Array<FadPCEType> x(n), maxs(n), maxs2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      x[i].fastAccessDx(j) = f;
    }
  }
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);
    maxs[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++)
      maxs[i].fastAccessDx(j) = f;
  }
  Teuchos::scan(*comm, *setup.fad_pce_serializer,
                Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);
  success = checkPCEArrays(maxs, maxs2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Scan Max", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_ScanMin ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int n = 7;
  int p = 5;
  int rank = comm->getRank();

  Teuchos::Array<FadPCEType> x(n), mins(n), mins2(n);
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);
    x[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++) {
      x[i].fastAccessDx(j) = f;
    }
  }
  for (int i=0; i<n; i++) {
    PCEType f(setup.kokkos_cijk);
    for (int k=0; k<setup.sz; k++)
      f.fastAccessCoeff(k) = 2.0*(i+1);
    mins[i] = FadPCEType(p, f);
    for (int j=0; j<p; j++)
      mins[i].fastAccessDx(j) = f;
  }
  Teuchos::scan(*comm, *setup.fad_pce_serializer,
                Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);
  success = checkPCEArrays(mins, mins2,
                           std::string("DFad")+"<"+"UQ::PCE"+"> Scan Min", out);
  success = checkResultOnAllProcs(*comm, out, success);
}

TEUCHOS_UNIT_TEST( UQ_PCE_Comm, FadPCE_SendReceive ) {
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();
  UnitTestSetup<PCEType, FadType> setup;

  int num_proc = comm->getSize();
  if (num_proc > 1) {
    int rank = comm->getRank();
    int n = 7;
    int p = 5;
    Teuchos::Array<FadPCEType> x(n), x2(n);
    for (int i=0; i<n; i++) {
      PCEType f(setup.kokkos_cijk);
      for (int k=0; k<setup.sz; k++)
        f.fastAccessCoeff(k) = 2.0*(i+1)*(k+1);
      x[i] = FadPCEType(p, f);
      for (int j=0; j<p; j++)
        x[i].fastAccessDx(j) = f;
    }
    if (rank != 1)
      x2 = x;
    if (rank == 0) Teuchos::send(*comm, *setup.fad_pce_serializer,
                                 n, &x[0], 1);
    if (rank == 1) Teuchos::receive(*comm, *setup.fad_pce_serializer,
                                    0, n, &x2[0]);
    success = checkPCEArrays(x, x2,
                             std::string("DFad")+"<"+"UQ::PCE"+"> Send/Receive", out);
    success = checkResultOnAllProcs(*comm, out, success);
  }
  else
    success = true;
}

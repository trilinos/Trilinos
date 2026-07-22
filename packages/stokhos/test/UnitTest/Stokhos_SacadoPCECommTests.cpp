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
  RCP< Stokhos::AlgebraicOrthogPolyExpansion<int,double> > exp;

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
    Cijk = basis->computeTripleProductTensor();
      
    // Expansion
    exp = 
      rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis, Cijk));

    // Serializers
    pce_serializer = 
      rcp(new PCESerializerT(
	    exp,
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

#define PCE_COMM_TESTS(PCEType, FadType, PCE, FAD)			\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_Broadcast ) {			\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  Teuchos::Array<PCEType> x(n), x2(n);					\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = rnd.number();				\
  }									\
  if (comm->getRank() == 0)						\
    x2 = x;								\
  Teuchos::broadcast(*comm, *setup.pce_serializer, 0, n, &x2[0]);	\
  success = checkPCEArrays(x, x2, std::string(#PCE)+" Broadcast", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_GatherAll ) {			\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int size = comm->getSize();						\
  int rank = comm->getRank();						\
  int N = n*size;							\
  Teuchos::Array<PCEType> x(n), x2(N), x3(N);				\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = (rank+1)*(i+1)*(j+1);			\
  }									\
  for (int j=0; j<size; j++) {						\
    for (int i=0; i<n; i++) {						\
      x3[n*j+i] = PCEType(setup.exp);					\
      for (int k=0; k<setup.sz; k++)					\
	x3[n*j+i].fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);		\
    }									\
  }									\
  Teuchos::gatherAll(*comm, *setup.pce_serializer,			\
		     n, &x[0], N, &x2[0]);				\
  success = checkPCEArrays(x3, x2, std::string(#PCE)+" Gather All", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_SumAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int num_proc = comm->getSize();					\
									\
  Teuchos::Array<PCEType> x(n), sums(n), sums2(n);			\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  for (int i=0; i<n; i++) {						\
    sums[i] = PCEType(setup.exp);					\
    for (int j=0; j<setup.sz; j++)					\
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;			\
  }									\
  Teuchos::reduceAll(*comm, *setup.pce_serializer,			\
		     Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);		\
  success = checkPCEArrays(sums, sums2,					\
			   std::string(#PCE)+" Sum All", out);		\
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_MaxAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int rank = comm->getRank();						\
  int num_proc = comm->getSize();					\
									\
  Teuchos::Array<PCEType> x(n), maxs(n), maxs2(n);			\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs[i] = PCEType(setup.exp);					\
    for (int j=0; j<setup.sz; j++)					\
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*num_proc;			\
  }									\
  Teuchos::reduceAll(*comm, *setup.pce_serializer,			\
		     Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);		\
  success = checkPCEArrays(maxs, maxs2,					\
			   std::string(#PCE)+" Max All", out);		\
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_MinAll ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<PCEType> x(n), mins(n), mins2(n);			\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    mins[i] = PCEType(setup.exp);					\
    for (int j=0; j<setup.sz; j++)					\
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  Teuchos::reduceAll(*comm, *setup.pce_serializer,			\
		     Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);		\
  success = checkPCEArrays(mins, mins2,					\
			   std::string(#PCE)+" Min All", out);		\
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_ScanSum ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<PCEType> x(n), sums(n), sums2(n);			\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  for (int i=0; i<n; i++) {						\
    sums[i] = PCEType(setup.exp);					\
    for (int j=0; j<setup.sz; j++)					\
      sums[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  Teuchos::scan(*comm, *setup.pce_serializer,				\
		Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);		\
  success = checkPCEArrays(sums, sums2,					\
			   std::string(#PCE)+" Scan Sum", out);		\
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_ScanMax ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<PCEType> x(n), maxs(n), maxs2(n);			\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    maxs[i] = PCEType(setup.exp);					\
    for (int j=0; j<setup.sz; j++)					\
      maxs[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  Teuchos::scan(*comm, *setup.pce_serializer,				\
		Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);		\
  success = checkPCEArrays(maxs, maxs2,					\
			   std::string(#PCE)+" Scan Max", out);		\
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_ScanMin ) {				\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<PCEType> x(n), mins(n), mins2(n);			\
  for (int i=0; i<n; i++) {						\
    x[i] = PCEType(setup.exp);						\
    for (int j=0; j<setup.sz; j++)					\
      x[i].fastAccessCoeff(j) = 2.0*(i+1)*(rank+1);			\
  }									\
  for (int i=0; i<n; i++) {						\
    mins[i] = PCEType(setup.exp);					\
    for (int j=0; j<setup.sz; j++)					\
      mins[i].fastAccessCoeff(j) = 2.0*(i+1);				\
  }									\
  Teuchos::scan(*comm, *setup.pce_serializer,				\
		Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);		\
  success = checkPCEArrays(mins, mins2,					\
			   std::string(#PCE)+" Scan Min", out);		\
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, PCE_SendReceive ) {			\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int num_proc = comm->getSize();					\
  if (num_proc > 1) {							\
    int rank = comm->getRank();						\
    int n = 7;								\
    Teuchos::Array<PCEType> x(n), x2(n);				\
    for (int i=0; i<n; i++) {						\
      x[i] = PCEType(setup.exp);					\
      for (int j=0; j<setup.sz; j++)					\
	x[i].fastAccessCoeff(j) = 2.0*(i+1)*(j+1);			\
    }									\
    if (rank != 1)							\
      x2 = x;								\
    if (rank == 0) Teuchos::send(*comm, *setup.pce_serializer,		\
				 n, &x[0], 1);				\
    if (rank == 1) Teuchos::receive(*comm, *setup.pce_serializer,	\
				    0, n, &x2[0]);			\
    success = checkPCEArrays(x, x2,					\
			     std::string(#PCE)+" Send/Receive", out);	\
    success = checkResultOnAllProcs(*comm, out, success);		\
  }									\
  else									\
    success = true;							\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_Broadcast ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  Teuchos::Array<FadPCEType> x(n), x2(n);				\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = rnd.number();				\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      PCEType g(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	g.fastAccessCoeff(k) = rnd.number();				\
      x[i].fastAccessDx(j) = g;						\
    }									\
  }									\
  if (comm->getRank() == 0)						\
    x2 = x;								\
  Teuchos::broadcast(*comm, *setup.fad_pce_serializer, 0, n, &x2[0]);	\
  success = checkPCEArrays(x, x2,					\
			   std::string(#FAD)+"<"+#PCE+"> Broadcast", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_GatherAll ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int size = comm->getSize();						\
  int rank = comm->getRank();						\
  int N = n*size;							\
  Teuchos::Array<FadPCEType> x(n), x2(N), x3(N);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = (rank+1)*(i+1)*(k+1);			\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      x[i].fastAccessDx(j) = f;						\
    }									\
  }									\
  for (int j=0; j<size; j++) {						\
    for (int i=0; i<n; i++) {						\
      PCEType f(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	f.fastAccessCoeff(k) = (j+1)*(i+1)*(k+1);			\
      x3[n*j+i] = FadPCEType(p, f);					\
      for (int k=0; k<p; k++)						\
	x3[n*j+i].fastAccessDx(k) = f;					\
    }									\
  }									\
  Teuchos::gatherAll(*comm, *setup.fad_pce_serializer,			\
		     n, &x[0], N, &x2[0]);				\
  success = checkPCEArrays(x3, x2,					\
			   std::string(#FAD)+"<"+#PCE+">  Gather All", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_SumAll ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int num_proc = comm->getSize();					\
									\
  Teuchos::Array<FadPCEType> x(n), sums(n), sums2(n);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      PCEType g(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	g.fastAccessCoeff(k) = 2.0*(i+1);				\
      x[i].fastAccessDx(j) = g;						\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;			\
    sums[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      PCEType g(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	g.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;			\
      sums[i].fastAccessDx(j) = g;					\
    }									\
  }									\
  Teuchos::reduceAll(*comm, *setup.fad_pce_serializer,			\
		     Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);		\
  success = checkPCEArrays(sums, sums2,					\
			   std::string(#FAD)+"<"+#PCE+"> Sum All", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_MaxAll ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 8;								\
  int p = 5;								\
  int rank = comm->getRank();						\
  int num_proc = comm->getSize();					\
									\
  Teuchos::Array<FadPCEType> x(n), maxs(n), maxs2(n);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      x[i].fastAccessDx(j) = f;						\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*num_proc;			\
    maxs[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++)						\
      maxs[i].fastAccessDx(j) = f;					\
  }									\
  Teuchos::reduceAll(*comm, *setup.fad_pce_serializer,			\
		     Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);		\
  success = checkPCEArrays(maxs, maxs2,					\
			   std::string(#FAD)+"<"+#PCE+"> Max All", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_MinAll ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 8;								\
  int p = 5;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<FadPCEType> x(n), mins(n), mins2(n);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      x[i].fastAccessDx(j) = f;						\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    mins[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++)						\
      mins[i].fastAccessDx(j) = f;					\
  }									\
  Teuchos::reduceAll(*comm, *setup.fad_pce_serializer,			\
		     Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);		\
  success = checkPCEArrays(mins, mins2,					\
			   std::string(#FAD)+"<"+#PCE+"> Min All", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_ScanSum ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<FadPCEType> x(n), sums(n), sums2(n);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      x[i].fastAccessDx(j) = f;						\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    sums[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++)						\
      sums[i].fastAccessDx(j) = f;					\
  }									\
  Teuchos::scan(*comm, *setup.fad_pce_serializer,			\
		Teuchos::REDUCE_SUM, n, &x[0], &sums2[0]);		\
  success = checkPCEArrays(sums, sums2,					\
			   std::string(#FAD)+"<"+#PCE+"> Scan Sum", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_ScanMax ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<FadPCEType> x(n), maxs(n), maxs2(n);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      x[i].fastAccessDx(j) = f;						\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    maxs[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++)						\
      maxs[i].fastAccessDx(j) = f;					\
  }									\
  Teuchos::scan(*comm, *setup.fad_pce_serializer,			\
		Teuchos::REDUCE_MAX, n, &x[0], &maxs2[0]);		\
  success = checkPCEArrays(maxs, maxs2,					\
			   std::string(#FAD)+"<"+#PCE+"> Scan Max", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_ScanMin ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int n = 7;								\
  int p = 5;								\
  int rank = comm->getRank();						\
									\
  Teuchos::Array<FadPCEType> x(n), mins(n), mins2(n);			\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1)*(rank+1);			\
    x[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++) {						\
      x[i].fastAccessDx(j) = f;						\
    }									\
  }									\
  for (int i=0; i<n; i++) {						\
    PCEType f(setup.exp);						\
    for (int k=0; k<setup.sz; k++)					\
      f.fastAccessCoeff(k) = 2.0*(i+1);					\
    mins[i] = FadPCEType(p, f);						\
    for (int j=0; j<p; j++)						\
      mins[i].fastAccessDx(j) = f;					\
  }									\
  Teuchos::scan(*comm, *setup.fad_pce_serializer,			\
		Teuchos::REDUCE_MIN, n, &x[0], &mins2[0]);		\
  success = checkPCEArrays(mins, mins2,					\
			   std::string(#FAD)+"<"+#PCE+"> Scan Min", out); \
  success = checkResultOnAllProcs(*comm, out, success);			\
}									\
									\
TEUCHOS_UNIT_TEST( PCE##_Comm, FadPCE_SendReceive ) {			\
  typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;		\
  Teuchos::RCP<const Teuchos::Comm<Ordinal> >				\
    comm = Teuchos::DefaultComm<Ordinal>::getComm();			\
									\
  int num_proc = comm->getSize();					\
  if (num_proc > 1) {							\
    int rank = comm->getRank();						\
    int n = 7;								\
    int p = 5;								\
    Teuchos::Array<FadPCEType> x(n), x2(n);				\
    for (int i=0; i<n; i++) {						\
      PCEType f(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	f.fastAccessCoeff(k) = 2.0*(i+1)*(k+1);				\
      x[i] = FadPCEType(p, f);						\
      for (int j=0; j<p; j++)						\
	x[i].fastAccessDx(j) = f;					\
    }									\
    if (rank != 1)							\
      x2 = x;								\
    if (rank == 0) Teuchos::send(*comm, *setup.fad_pce_serializer,	\
				 n, &x[0], 1);				\
    if (rank == 1) Teuchos::receive(*comm, *setup.fad_pce_serializer,	\
				    0, n, &x2[0]);			\
    success = checkPCEArrays(x, x2,					\
			     std::string(#FAD)+"<"+#PCE+"> Send/Receive", out);	\
    success = checkResultOnAllProcs(*comm, out, success);		\
  }									\
  else									\
    success = true;							\
}

typedef int Ordinal;
typedef Stokhos::StandardStorage<int,double> storage_type;
typedef Sacado::Fad::DFad<double> fad_type;
namespace PCETest {
  Sacado::Random<double> rnd;
  typedef Sacado::PCE::OrthogPoly<double,storage_type> pce_type;
  UnitTestSetup<pce_type, fad_type> setup;
  PCE_COMM_TESTS(pce_type, fad_type, OrthogPoly, DFad)
}

namespace ETPCETest {
  Sacado::Random<double> rnd;
  typedef Sacado::ETPCE::OrthogPoly<double,storage_type> pce_type;
  UnitTestSetup<pce_type, fad_type> setup;
  PCE_COMM_TESTS(pce_type, fad_type, ETOrthogPoly, DFad)
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}

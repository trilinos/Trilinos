// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestILUT.cpp

\brief Single MPI rank Ifpack2 Unit test for the ILUT template.
*/


#include <iostream>
#include <type_traits>
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_ConfigDefs.hpp>

// Xpetra / Galeri
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

#include <Ifpack2_Version.hpp>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_ILUT.hpp>


namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2ILUT, ParILUT, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  typedef Node NT;

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,NT>   crs_matrix_type;

  typedef Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,NT> XCrsType;
  typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,NT>                XMapType;
  typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,NT>     XMVectorType;

  // Generate the matrix using Galeri.  Galeri wraps it in an Xpetra
  // matrix, so after it finishes, ask it for the Tpetra matrix.
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  Teuchos::CommandLineProcessor clp;

//#define IFPACK2_DEBUG_PARILUT // This duplicates a Kokkos-Kernels unit test.
#ifdef IFPACK2_DEBUG_PARILUT
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> tpetra_map_type;

  GlobalOrdinal nx = 2, ny=2;
  const Teuchos::RCP<const tpetra_map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(4);
  RCP<crs_matrix_type> mtx = Teuchos::rcp(new crs_matrix_type(rowmap,4));
  Teuchos::Array<GlobalOrdinal> inds(4);
  Teuchos::Array<Scalar>        vals(4);

  inds[0]=0; inds[1]=1; inds[2]=2; inds[3]=3;
  vals[0]=1.; vals[1]=6.; vals[2]=4.; vals[3]=7.;
  mtx->insertGlobalValues(0, 4, vals.getRawPtr(), inds.getRawPtr());
  vals[0]=2.; vals[1]=-5.; vals[2]=0.; vals[3]=8.;
  mtx->insertGlobalValues(1, 4, vals.getRawPtr(), inds.getRawPtr());
  vals[0]=0.5; vals[1]=-3.; vals[2]=6.; vals[3]=0.;
  mtx->insertGlobalValues(2, 4, vals.getRawPtr(), inds.getRawPtr());
  vals[0]=0.2; vals[1]=-0.5; vals[2]=-9.; vals[3]=0.;
  mtx->insertGlobalValues(3, 4, vals.getRawPtr(), inds.getRawPtr());
  mtx->fillComplete();
  RCP<crs_matrix_type> A = mtx;
#else
  GlobalOrdinal nx = 30, ny=30, nz=3;
  Galeri::Xpetra::Parameters<GlobalOrdinal> GaleriParameters (clp, nx, ny, nz, "Laplace2D");
  Xpetra::Parameters xpetraParameters (clp);
  ParameterList GaleriList = GaleriParameters.GetParameterList ();
  RCP<XMapType> xmap =
    Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node> (xpetraParameters.GetLib (),
                                             "Cartesian2D", comm, GaleriList);
  RCP<Galeri::Xpetra::Problem<XMapType,XCrsType,XMVectorType> > Pr =
    Galeri::Xpetra::BuildProblem<Scalar,LocalOrdinal,GlobalOrdinal,XMapType,XCrsType,XMVectorType> (std::string("Laplace2D"),
                                                                           xmap, GaleriList);

  RCP<XCrsType> XA = Pr->BuildMatrix ();
  RCP<crs_matrix_type> A = XA->getTpetra_CrsMatrixNonConst ();
#endif
  TEST_INEQUALITY(A, Teuchos::null);


  using prec_type = Ifpack2::ILUT<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >;

  {
    prec_type prec (A);
    Teuchos::ParameterList params;
    params.set("fact: type", "par_ilut");
    params.set("fact: ilut level-of-fill",0.5);
    TEST_THROW( prec.setParameters (params), std::runtime_error);
  }

  prec_type prec (A);
  {
    Teuchos::ParameterList params;
    params.set ("fact: type", "par_ilut"); // valid, not default
    TEST_NOTHROW( prec.setParameters (params) );
  }

  out << "prec.initialize()" << std::endl;
  prec.initialize();
  out << "prec.compute()" << std::endl;
  prec.compute();

#ifdef IFPACK2_DEBUG_PRINT_PARILUT_FACTORS
  std::cout << "========\nA matrix\n========" << std::endl;
  int nnz=nx*ny;
  for (int i=0; i<nnz; ++i) {
    using h_inds_type = typename crs_matrix_type::nonconst_local_inds_host_view_type;
    using h_vals_type = typename crs_matrix_type::nonconst_values_host_view_type;
    h_inds_type indices("A columns", nnz);
    h_vals_type values("A values", nnz);
    size_t numEntries;
    A->getLocalRowCopy(i,indices,values,numEntries);
    std::cout << "row " << i << ":";
    for (int j=0; j<nnz; j++)
      //std::cout << "  " << values[j];
      printf("   %10.7lf",values[j]);
    std::cout << std::endl;
  }

  std::cout << std::endl << std::endl << "========\nL factor\n========" << std::endl;
  auto Lfactor = prec.getL();
  for (int i=0; i<nnz; ++i) {
    using h_inds_type = typename crs_matrix_type::nonconst_local_inds_host_view_type;
    using h_vals_type = typename crs_matrix_type::nonconst_values_host_view_type;
    h_inds_type indices("L columns", nnz);
    h_inds_type inds2print("L columns 2 print", nnz);
    h_vals_type values("L values", nnz);
    h_vals_type vals2print("L values 2 print", nnz);
    size_t numEntries;
    Lfactor->getLocalRowCopy(i,indices,values,numEntries);
    for (size_t j=0; j<numEntries; j++) {
      inds2print[indices[j]]=indices[j];
      vals2print[indices[j]]=values[j];
    }
    std::cout << "row " << i << ":";
    for (int j=0; j<nnz; j++)
      printf("   %10.7lf (%d)",vals2print[j],inds2print[j]);
    std::cout << std::endl;
  }

  std::cout << std::endl << std::endl << "========\nU factor\n========" << std::endl;
  auto Ufactor = prec.getU();
  for (int i=0; i<nnz; ++i) {
    using h_inds_type = typename crs_matrix_type::nonconst_local_inds_host_view_type;
    using h_vals_type = typename crs_matrix_type::nonconst_values_host_view_type;
    h_inds_type indices("U columns", nnz);
    h_vals_type values("U values", nnz);
    h_inds_type inds2print("U columns 2 print", nnz);
    h_vals_type vals2print("U values 2 print", nnz);
    size_t numEntries;
    Ufactor->getLocalRowCopy(i,indices,values,numEntries);
    for (size_t j=0; j<numEntries; j++) {
      inds2print[indices[j]]=indices[j];
      vals2print[indices[j]]=values[j];
    }
    std::cout << "row " << i << ":";
    for (int j=0; j<nnz; j++)
      printf("   %10.7lf (%d)",vals2print[j],inds2print[j]);
    std::cout << std::endl;
  }

#endif //ifdef IFPACK2_DEBUG_PRINT_PARILUT_FACTORS

  const int numVectors=2;
  using STS = Teuchos::ScalarTraits<Scalar>;
  using STM = typename STS::magnitudeType;
  Kokkos::View<STM*, Kokkos::HostSpace> bnorms("Initial norms", numVectors);
  Kokkos::View<STM*, Kokkos::HostSpace> xnorms_final("Initial norms", numVectors);
  Kokkos::View<STM*, Kokkos::HostSpace> xnorms_true("Initial norms", numVectors);
  Kokkos::View<STM*, Kokkos::HostSpace> norms("Initial norms", numVectors);

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(A->getRowMap(),numVectors), b(A->getRowMap(),numVectors);
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> r(A->getRowMap(),numVectors);
  x.randomize();
  x.norm2(xnorms_true);
  out << "before solve, |xact|: " << ": " << xnorms_true(0) << ", " << xnorms_true(1) << std::endl;

  A->apply(x,b);
  x.putScalar(0);

  b.norm2(bnorms);
  Kokkos::View<STM*, Kokkos::HostSpace> lastNorms("previous norms", numVectors);
  prec.apply(b, x);
  A->apply(x,r);
  r.update(1.,b,-1.);
  x.norm2(xnorms_final);
  out << "|e|        = [" << std::abs(xnorms_true(0)-xnorms_final(0)) << ", " << std::abs(xnorms_true(1)-xnorms_final(1)) << "]" << std::endl;
  r.norm2(norms);
  out << "|b|        = [" << bnorms(0) << ", " << bnorms(1) << "]" << std::endl;
  out << "|b-Ax|     = [" << norms(0) << ", " << norms(1) << "]" << std::endl;
  out << "|b-Ax|/|b| = [" << norms(0)/bnorms(0) << ", " << norms(1)/bnorms(1) << "]" << std::endl;
  TEST_COMPARE(norms(0), <, 0.25 * bnorms(0));
  TEST_COMPARE(norms(1), <, 0.25 * bnorms(1));
} //ParILUT


#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2ILUT, ParILUT, Scalar, LocalOrdinal,GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

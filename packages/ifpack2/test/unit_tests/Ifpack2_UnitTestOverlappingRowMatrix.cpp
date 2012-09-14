// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Ifpack2_UnitTestRelaxation.cpp

\brief Ifpack2 Unit test for the OverlappingRowMatrix template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#ifdef HAVE_IFPACK2_QD
#include <qd/dd_real.h>
#endif

// Xpetra / Galeri
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>
#include <Galeri_XpetraMatrixTypes.hpp>


#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_OverlappingRowMatrix.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;


//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2OverlappingRowMatrix, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  // Typedefs
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>       CrsType;
  typedef Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XCrsType;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                    MapType;
  typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                    XMapType;
  typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>          VectorType;

  // Useful stuff
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  Teuchos::CommandLineProcessor clp;
  Xpetra::Parameters xpetraParameters(clp);

  Teuchos::ParameterList GaleriList;
  int nx = 100; 
  size_t numProcs = comm->getSize();
  size_t numElementsPerProc = nx*nx;
  GaleriList.set("nx", (GlobalOrdinal) nx);
  GaleriList.set("ny", (GlobalOrdinal) (nx * numProcs));
  GaleriList.set("n", (GlobalOrdinal) (numElementsPerProc*numProcs));

  // Short circuit --- this test should only be run in parallel.
  if (comm->getSize() == 1) return;

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  Teuchos::RCP<XMapType > xmap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal>::Build(xpetraParameters.GetLib(), INVALID, numElementsPerProc, 0, comm);
  Teuchos::RCP<XCrsType> XA = Galeri::Xpetra::CreateCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,XMapType,XCrsType>(string("Laplace2D"),xmap,GaleriList);
  Teuchos::RCP<CrsType> A = XA->getTpetra_CrsMatrixNonConst();

  TEST_INEQUALITY(A,Teuchos::null);
  
  VectorType X(A->getRowMap()), Y(A->getRowMap());
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > x_ptr = X.get2dViewNonConst();

  int OverlapLevel = 5;  

  // ======================================== //
  // Build the overlapping matrix using class //
  // Ifpack2::OverlappingRowMatrix.             //
  // ======================================== //
  Ifpack2::OverlappingRowMatrix<CrsType> B(A,OverlapLevel);
  
  for (LocalOrdinal i = 0 ; (size_t)i < A->getNodeNumRows() ; ++i) 
    x_ptr[0][i] = 1.0* A->getRowMap()->getGlobalElement(i);
  Y.putScalar(0.0);

  VectorType ExtX_B(B.getRowMap()), ExtY_B(B.getRowMap());
  ExtY_B.putScalar(0.0);





#ifdef OLD_AND_BUSTED
  IFPACK_CHK_ERR(B.ImportMultiVector(X,ExtX_B));
  IFPACK_CHK_ERR(B.Multiply(false,ExtX_B,ExtY_B));
  IFPACK_CHK_ERR(B.ExportMultiVector(ExtY_B,Y,Add));

  double Norm_B;
  Y.Norm2(&Norm_B);
  if (Comm.MyPID() == 0)
    cout << "Norm of Y using B = " << Norm_B << endl;
  
  // ================================================== //
  //Build the overlapping matrix as an Epetra_CrsMatrix //
  // ================================================== //

  Time.ResetStartTime();
  Epetra_CrsMatrix& C = 
    *(Ifpack_CreateOverlappingCrsMatrix(&*A,OverlapLevel));
  if (Comm.MyPID() == 0)
    cout << "Time to create C = " << Time.ElapsedTime() << endl;

  // simple checks on global quantities
  int NumGlobalRowsC = C.NumGlobalRows();
  int NumGlobalNonzerosC = C.NumGlobalNonzeros();
  assert (NumGlobalRowsB == NumGlobalRowsC);
  assert (NumGlobalNonzerosB == NumGlobalNonzerosC);

  Epetra_Vector ExtX_C(C.RowMatrixRowMap());
  Epetra_Vector ExtY_C(C.RowMatrixRowMap());
  ExtY_C.PutScalar(0.0);
  Y.PutScalar(0.0);

  IFPACK_CHK_ERR(C.Multiply(false,X,Y));

  double Norm_C;
  Y.Norm2(&Norm_C);
  if (Comm.MyPID() == 0)
    cout << "Norm of Y using C = " << Norm_C << endl;

  if (IFPACK_ABS(Norm_B - Norm_C) > 1e-5)
    IFPACK_CHK_ERR(-1);

  // ======================= //
  // now localize the matrix //
  // ======================= //

  Ifpack_LocalFilter D(Teuchos::rcp(&B, false));

#endif
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2OverlappingRowMatrix, Test0, Scalar, LocalOrdinal,GlobalOrdinal) 


UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)

#ifdef HAVE_IFPACK2_QD
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>


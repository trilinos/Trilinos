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


/*! \file Ifpack2_UnitTestContainer.cpp

\brief Ifpack2 Unit test for the Container classes
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#ifdef HAVE_IFPACK2_QD
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_SparseContainer.hpp>
#include <Ifpack2_ILUT.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Container, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CRS;
  typedef Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,LocalOrdinal,Node>    > ILUTlo;
  typedef Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   > ILUTgo;

  // The simple joy of tridiagonal matrices
  global_size_t num_rows_per_proc = 5;
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap), y(rowmap), z(rowmap);
  Teuchos::ArrayRCP<Scalar> x_ptr = x.get1dViewNonConst();
  Teuchos::ArrayRCP<Scalar> y_ptr = y.get1dViewNonConst();

  // Fill x for all time
  Teuchos::ScalarTraits<double>::seedrandom(24601);
  x.randomize();

  // ====================================== //
  //        Sparse Container + ILUT         //
  // ====================================== //
  // 
  Ifpack2::SparseContainer<CRS,ILUTlo> MyContainer((size_t)num_rows_per_proc);
  Teuchos::ParameterList params;
  params.set("fact: ilut level-of-fill", 1.0);
  params.set("fact: drop tolerance", 0.0);
  MyContainer.setParameters(params);

  // Set IDs to grab the whole matrix
  MyContainer.initialize();
  for(size_t i=0; i<num_rows_per_proc; i++)
    MyContainer.ID(i) = i;
  MyContainer.compute(crsmatrix);

  // Reference ILUT
  ILUTgo prec(crsmatrix);
  prec.setParameters(params);
  prec.initialize();
  prec.compute();

  // Apply the SparseContainer
  Teuchos::ArrayRCP<Scalar> x_local = MyContainer.getX()->get1dViewNonConst();
  Teuchos::ArrayRCP<Scalar> y_local = MyContainer.getY()->get1dViewNonConst();
  for(size_t i=0; i <(size_t)num_rows_per_proc; i++) x_local[i] = x_ptr[i];
  MyContainer.apply();
  for(size_t i=0; i <(size_t)num_rows_per_proc; i++) y_ptr[i] = y_local[i];
  
  // Apply raw ILUT
  prec.apply(x,z);

  // Diff
  TEST_COMPARE_FLOATING_ARRAYS(y.get1dView(), z.get1dView(), 1e4*Teuchos::ScalarTraits<Scalar>::eps());

}



#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Container, Test0, Scalar, LocalOrdinal,GlobalOrdinal) 

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

#ifdef HAVE_IFPACK2_QD
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

}//namespace <anonymous>


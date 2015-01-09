/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************


/// \file Ifpack2_UnitTestAmesos2solver.cpp
/// \brief Ifpack2 unit test for the Amesos2 wrapper.
///
/// This test only builds nontrivially and executes if Amesos2 is enabled.

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>

//#if defined(HAVE_IFPACK2_AMESOS2) && defined(HAVE_AMESOS2_SUPERLU)
#if defined(HAVE_IFPACK2_AMESOS2)

#include <Amesos2_config.h>
#include <Ifpack2_Details_Amesos2Wrapper.hpp>
#include <iostream>

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
#include <qd/dd_real.h>
#endif

#include <Ifpack2_UnitTestHelpers.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Amesos2Wrapper, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success
  using Teuchos::RCP;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> mv_type;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);

  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> (rowmap);

  Ifpack2::Details::Amesos2Wrapper<crs_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;
  params.set("Amesos2 solver name","superlu");
  Teuchos::ParameterList &sublist = params.sublist("Amesos2");
  (sublist.sublist("SuperLU")).set("ILU_Flag",false); //create SuperLU sublist to get rid of unused variable warnings
  TEST_NOTHROW(prec.setParameters(params));

  //trivial tests to insist that the preconditioner's domain/range maps are
  //identically those of the matrix:
  const map_type* mtx_dom_map_ptr = &*crsmatrix->getDomainMap();
  const map_type* mtx_rng_map_ptr = &*crsmatrix->getRangeMap();

  const map_type* prec_dom_map_ptr = &*prec.getDomainMap();
  const map_type* prec_rng_map_ptr = &*prec.getRangeMap();

  TEST_EQUALITY( prec_dom_map_ptr, mtx_dom_map_ptr );
  TEST_EQUALITY( prec_rng_map_ptr, mtx_rng_map_ptr );

# if !defined(HAVE_AMESOS2_SUPERLU)
  TEST_THROW(prec.initialize(),std::invalid_argument);
# else
  TEST_NOTHROW(prec.initialize());
  prec.compute();

  mv_type x (rowmap,2), y (rowmap,2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  prec.apply (x, y);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());
# endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2Amesos2Wrapper, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success
  using Teuchos::RCP;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> mv_type;

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  RCP<const map_type> rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node> (num_rows_per_proc);

  RCP<const crs_matrix_type> crsmatrix = tif_utest::create_test_matrix3<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::Details::Amesos2Wrapper<crs_matrix_type> prec (crsmatrix);

  Teuchos::ParameterList params;

  params.set("Amesos2 solver name","superlu");
  Teuchos::ParameterList &sublist = params.sublist("Amesos2");
  (sublist.sublist("SuperLU")).set("ILU_Flag",false); //create SuperLU sublist to get rid of unused variable warnings
  TEST_NOTHROW(prec.setParameters(params));

# if !defined(HAVE_AMESOS2_SUPERLU)
  TEST_THROW(prec.initialize(),std::invalid_argument);
# else
  TEST_NOTHROW(prec.initialize());
  prec.compute();

  mv_type x(rowmap,2), y(rowmap,2);
  x.putScalar (Teuchos::ScalarTraits<Scalar>::one ());

  crsmatrix->apply (x, y);
  prec.apply (y, x);

  Teuchos::ArrayRCP<const Scalar> xview = x.get1dView();

  //x should be full of 1's now.

  Teuchos::ArrayRCP<Scalar> ones(num_rows_per_proc*2, Teuchos::ScalarTraits<Scalar>::one ());

  TEST_COMPARE_FLOATING_ARRAYS(xview, ones(), 2*Teuchos::ScalarTraits<Scalar>::eps());
# endif
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Amesos2Wrapper, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2Amesos2Wrapper, Test1, Scalar, LocalOrdinal,GlobalOrdinal)

//typedef std::complex<double> ComplexDouble;
//UNIT_TEST_GROUP_SCALAR_ORDINAL(ComplexDouble, int, int)
UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)
#ifndef HAVE_IFPACK2_EXPLICIT_INSTANTIATION
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, short, int)
#endif

#if defined(HAVE_IFPACK2_QD) && !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
UNIT_TEST_GROUP_SCALAR_ORDINAL(dd_real, int, int)
#endif

} // namespace (anonymous)

#endif // HAVE_IFPACK2_AMESOS2

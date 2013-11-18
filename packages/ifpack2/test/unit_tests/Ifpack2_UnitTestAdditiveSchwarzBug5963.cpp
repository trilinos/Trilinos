/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
//      Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

/// \file Ifpack2_UnitTestAdditiveSchwarzOverlap.cpp
// \brief Unit test exercising different overlap with AdditiveSchwarz.

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>
#include <Ifpack2_Details_DenseSolver.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>


template<>
Ifpack2::Details::DenseSolver<Tpetra::RowMatrix<double, int, int>, false>;
template<>
Ifpack2::AdditiveSchwarz<Tpetra::RowMatrix<double, int, int>,
                         Ifpack2::Details::DenseSolver<Tpetra::RowMatrix<double, int, int>, false> >;


namespace {

// This example exercises the fix for Bug 5963.  It illustrates the
// "Zero" combine mode.  This example only works with exactly 2 MPI
// processes.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(AdditiveSchwarz, ZeroCombineMode, ScalarType, LocalOrdinalType, GlobalOrdinalType)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::toString;
  using std::endl;

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_type;

  typedef Tpetra::global_size_t GST;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vec_type;
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  typedef Ifpack2::Details::DenseSolver<row_matrix_type> local_solver_type;
  typedef Ifpack2::AdditiveSchwarz<row_matrix_type, local_solver_type> global_solver_type;

  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<node_type> node =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)

  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl;

  out << "Creating Maps" << endl;

  //
  // Each of the two processes gets two rows.
  //
  const size_t localNumRows = 2;
  const GST globalNumRows = comm->getSize () * localNumRows;
  const global_ordinal_type indexBase = 0;

  RCP<const map_type> rowMap (new map_type (globalNumRows, localNumRows, indexBase, comm, node));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Creating and filling matrix A" << endl;

  RCP<crs_matrix_type> A (new crs_matrix_type (rowMap, 3));

  Array<scalar_type> val (3);
  Array<global_ordinal_type> ind (3);
  val[0] = as<scalar_type> (-1);
  val[1] = as<scalar_type> (2);
  val[2] = as<scalar_type> (-1);

  if (myRank == 0) {
    ind[1] = indexBase;
    ind[2] = indexBase + 1;
    A->insertGlobalValues (indexBase, ind.view (1, 2), val.view (1, 2));

    ind[0] = indexBase;
    ind[1] = indexBase + 1;
    ind[2] = indexBase + 2;
    A->insertGlobalValues (indexBase, ind.view (0, 3), val.view (0, 3));
  }
  else if (myRank == 1) {
    ind[0] = indexBase + 1;
    ind[1] = indexBase + 2;
    ind[2] = indexBase + 3;
    A->insertGlobalValues (indexBase, ind.view (0, 3), val.view (0, 3));

    ind[0] = indexBase + 2;
    ind[1] = indexBase + 3;
    A->insertGlobalValues (indexBase, ind.view (0, 2), val.view (0, 2));
  }

  out << "Calling fillComplete on A" << endl;
  A->fillComplete (domMap, ranMap);

  out << "Creating x_exact and b for the test problem" << endl;
  // FIXME (mfh 17 Nov 2013) Insert actual x_exact and b.
  vec_type x_exact (domMap);
  {
    ArrayRCP<scalar_type> x_exact_view = x_exact.get1dViewNonConst ();
    x_exact_view[0] = 1;
    x_exact_view[1] = 4;
  }
  out << "x_exact: " << toString ((x_exact.get1dView ()) ()) << endl;

  vec_type b (ranMap);
  A->apply (x_exact, b); // b := A*x_exact
  out << "b: " << toString ((b.get1dView ()) ()) << endl;

  out << "Creating solver" << endl;
  global_solver_type solver (A);

  out << "Calling initialize" << endl;
  solver.initialize ();

  out << "Calling compute" << endl;
  solver.compute ();

  out << "Calling solver's describe() method" << endl;
  solver.describe (out, Teuchos::VERB_EXTREME);

  out << "Calling apply" << endl;
  vec_type x_computed (ranMap);
  solver.apply (b, x_computed); // solve A*x_computed=b for x_computed

  out << "x_computed: " << toString ((x_computed.get1dView ()) ()) << endl;

  // Compute the residual.
  vec_type r (ranMap);
  A->apply (x_computed, r); // r := A*x_computed
  r.update (STS::one (), b, -STS::one ()); // r := b - A*x_computed
  const magnitude_type absResNorm = r.norm2 ();
  out << "\\| b - A*x_computed \\|_2 = " << absResNorm << endl;
  const magnitude_type normB = b.norm2 ();
  const magnitude_type relResNorm = absResNorm / normB;
  out << "\\| b - A*x_computed \\|_2 / \\|b\\|_2 = " << relResNorm << endl;

  // // 'out' only prints on Proc 0.
  // RCP<Teuchos::FancyOStream> errPtr = Teuchos::getFancyOStream (rcpFromRef (std::cerr));
  // Teuchos::FancyOStream& err = *errPtr;
  // err << "A->describe() result:" << endl;
  // A->describe (err, Teuchos::VERB_EXTREME);

  // err << "x_exact.describe() result:" << endl;
  // x_exact.describe (err, Teuchos::VERB_EXTREME);

  // err << "b.describe() result:" << endl;
  // b.describe (err, Teuchos::VERB_EXTREME);

  // err << "x_computed.describe() result:" << endl;
  // x_computed.describe (err, Teuchos::VERB_EXTREME);

  // err << "r.describe() result:" << endl;
  // r.describe (err, Teuchos::VERB_EXTREME);

  // err << "solver.describe() result:" << endl;
  // solver.describe (err, Teuchos::VERB_EXTREME);

  // TEUCHOS_TEST_FOR_EXCEPTION(
  //   relResNorm > 10*STS::eps (), std::logic_error,
  //   "DenseSolver failed to solve the problem to within a small tolerance "
  //   << 10*STS::eps () << ".  Relative residual norm: " << relResNorm << ".");
}

// Define the set of unit tests to instantiate in this file.
// AdditiveSchwarz with a DenseSolver subdomain solver is not
// explicitly instantiated by default, so do that here.
#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( AdditiveSchwarz, ZeroCombineMode, Scalar, LocalOrdinal, GlobalOrdinal)


// Instantiate the unit tests for Scalar=double, LO=int, and GO=int.
// It's not necessary to exercise other Scalar types, as that would
// just be a Teuchos::LAPACK test, not an Ifpack2 test.
UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int, int)

}//namespace <anonymous>


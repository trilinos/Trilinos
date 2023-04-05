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


/*! \file Ifpack2_UnitTestMDF.cpp

\brief Ifpack2 Unit test for the MDF template.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_MDF.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>

#include "Tpetra_BlockCrsMatrix.hpp"
#include <Ifpack2_Experimental_RBILUK.hpp>

#include <type_traits>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2MDF, Test0, Scalar, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowmap = tif_utest::create_tpetra_map<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsmatrix = tif_utest::create_test_matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rowmap);

  Ifpack2::MDF<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec(crsmatrix);

  Teuchos::ParameterList params;
  params.set("fact: mdf level-of-fill", 0.0);
  params.set("Verbosity", 0);
  TEST_NOTHROW(prec.setParameters(params));

  prec.initialize();
  prec.compute();

  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(rowmap,2), y(rowmap,2);
  x.putScalar(1);

  prec.apply(x, y);

  Teuchos::ArrayRCP<const Scalar> yview = y.get1dView();

  //y should be full of 0.5's now.

  Teuchos::ArrayRCP<Scalar> halfs(num_rows_per_proc*2, 0.5);

  TEST_COMPARE_FLOATING_ARRAYS(yview, halfs(), Teuchos::ScalarTraits<Scalar>::eps());

  // Test that permuation arrays are updated not overwritten
  Teuchos::ArrayRCP<const LocalOrdinal> permuations = prec.getPermutations();
  Teuchos::ArrayRCP<const LocalOrdinal> reversePermuations = prec.getReversePermutations();
  prec.compute();
  TEST_EQUALITY(permuations,prec.getPermutations());
  TEST_EQUALITY(reversePermuations,prec.getReversePermutations());
}

template<size_t maxEntrPerRow,typename Scalar,typename LO,typename GO,typename Node>
void test_mdf_reference_problem(
  bool& success,
  Teuchos::FancyOStream& out,
  const LO locNumRow,
  const LO * crs_row_map,
  const LO * crs_col_ind,
  const Scalar * crs_values,
  const LO * known_Perm,
  const Scalar * known_sln)
{
  using row_matrix_t = Tpetra::RowMatrix<Scalar,LO,GO,Node>;
  using crs_matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > rowmap(tif_utest::create_tpetra_map<LO,GO,Node>(locNumRow));
  Teuchos::RCP<crs_matrix_t> crsmatrix(new crs_matrix_t(rowmap,maxEntrPerRow));

  //Fill matrix
  {
    GO rb = rowmap()->getGlobalElement(0);
    GO col_indices[maxEntrPerRow];

    for(LO row_lo = 0; row_lo < locNumRow; ++row_lo)
    {
      const LO & entrStart = crs_row_map[row_lo];
      const LO entriesInRow = crs_row_map[row_lo+1] - entrStart;
      for (LO i =0; i < entriesInRow; ++i)
        col_indices[i] = crs_col_ind[entrStart + i] + rb;

      crsmatrix->insertGlobalValues (row_lo + rb,
                          entriesInRow,
                          &crs_values[entrStart],
                          col_indices);
    }
    crsmatrix->fillComplete();
  }

  //Create prec
  Ifpack2::MDF<row_matrix_t> prec(crsmatrix.getConst());
  {
    Teuchos::ParameterList params;
    params.set("fact: mdf level-of-fill", 0.0);
    params.set("Verbosity",0);
    TEST_NOTHROW(prec.setParameters(params));
  }

  //Compute prec
  prec.initialize();
  prec.compute();

  //Check the permutations are correct
  if (known_Perm != nullptr)
  {
    Teuchos::ArrayRCP<LO> refPermuations(locNumRow);
    refPermuations.assign(known_Perm,known_Perm+locNumRow);
    TEST_COMPARE_ARRAYS(prec.getReversePermutations(),refPermuations);

    Teuchos::ArrayRCP<LO> refPermuationsInv(locNumRow);
    for (LO i=0;i<locNumRow;++i)
      refPermuationsInv[known_Perm[i]] = i;
    TEST_COMPARE_ARRAYS(prec.getPermutations(),refPermuationsInv);
  }

  //Apply prec
  Tpetra::MultiVector<Scalar,LO,GO,Node> x(rowmap,2), yMDF(rowmap,2), yILU(rowmap,2);
  x.putScalar(1);
  prec.apply(x, yMDF);

  // Check result if known
  if (known_sln != nullptr)
  {
    Teuchos::ArrayRCP<Scalar> knownSln(locNumRow*2);
    for (LO i=0;i<locNumRow;++i)
      knownSln[i] = knownSln[i+locNumRow] = known_sln[i];

    Teuchos::ArrayRCP<const Scalar> yview = yMDF.get1dView();

    const Scalar tol = 1e-8;
    TEST_COMPARE_FLOATING_ARRAYS(yview, knownSln, tol);
  }

  // Now apply reordering with AdditiveSchwarz
  Ifpack2::AdditiveSchwarz<row_matrix_t> reorderedPrec(crsmatrix.getConst());
  {
    Teuchos::ParameterList params;
    params.set ("schwarz: overlap level", static_cast<int> (0));
    params.set ("schwarz: combine mode", "add");
    params.set ("inner preconditioner name", "RILUK");
    params.set ("schwarz: zero starting solution", true);
    params.set ("schwarz: num iterations", 1);
    {
      Teuchos::ParameterList innerParams;
      innerParams.set ("fact: iluk level-of-fill", static_cast<int> (0));
      innerParams.set ("fact: iluk level-of-overlap", static_cast<int> (0));

      params.set ("inner preconditioner parameters", innerParams);
    }
    params.set ("schwarz: use reordering", true);
    {
      Teuchos::ParameterList zlist;
      zlist.set ("order_method", "user");
      zlist.set ("order_method_type", "local");
      zlist.set ("user ordering", prec.getPermutations());
      zlist.set ("user reverse ordering", prec.getReversePermutations());

      params.set ("schwarz: reordering list", zlist);
    }
    reorderedPrec.setParameters(params);
  }

  //Compute prec
  reorderedPrec.initialize();
  reorderedPrec.compute();

  //Apply prec
  reorderedPrec.apply(x, yILU);

  //Check if results match mdf impl
  {
    Teuchos::ArrayRCP<const Scalar> yMDFview = yMDF.get1dView();
    Teuchos::ArrayRCP<const Scalar> yILUview = yILU.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(yMDFview, yILUview, 100*Teuchos::ScalarTraits<Scalar>::eps());
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2MDF, Test1, Scalar, LocalOrdinal, GlobalOrdinal)
{
  // Reference problem from paragraph 2.5 of
  //   D'Azevedo, E. F., et al. "Ordering Methods for Preconditioned Conjugate Gradient Methods Applied to
  //     Unstructured Grid Problems." SIAM Journal on Matrix Analysis and Applications, vol. 13, no. 3, 
  //     July 1992, pp. 944-961., https://doi.org/10.1137/0613057. 
  // Known solution determined in Matlab.

  const LocalOrdinal locNumRow = 16;
  const LocalOrdinal maxEntrPerRow = 5;

  const LocalOrdinal known_Perm[] = {0, 3, 12, 15, 1, 2, 4, 8, 7, 11, 13, 14, 5, 6, 9, 10};
  const Scalar known_sln[] = {0.667850672, 0.8357013441, 0.7697237194, 0.6521907175, 0.8357013441, 1.114156321, 1.097870956, 0.8390391506, 0.7697237194, 1.097870956, 1.198198407, 0.7985258586, 0.6521907175, 0.8390391506, 0.7985258586, 0.6492629293};

  const LocalOrdinal crs_row_map[locNumRow + 1] = {0, 3, 7, 11, 14, 18, 23, 28, 32, 36, 41, 46, 50, 53, 57, 61, 64};
  const LocalOrdinal crs_col_ind[] = {0,   1,  4,
                                      0,   1,  2,  5,
                                      1,   2,  3,  6,
                                      2,   3,  7,
                                      0,   4,  5,  8,
                                      1,   4,  5,  6,  9,
                                      2,   5,  6,  7, 10,
                                      3,   6,  7, 11,
                                      4,   8,  9, 12,
                                      5,   8,  9, 10, 13,
                                      6,   9, 10, 11, 14,
                                      7,  10, 11, 15,
                                      8,  12, 13,
                                      9,  12, 13, 14,
                                      10, 13, 14, 15,
                                      11, 14, 15};
  const Scalar crs_values[] = {4, -1, -1,
                              -1,  4, -1, -1,
                              -1,  4, -1, -1,
                              -1,  4, -1,
                              -1,  4, -1, -1,
                              -1, -1,  4, -1, -1,
                              -1, -1,  4, -1, -1,
                              -1, -1,  4, -1,
                              -1,  4, -1, -1,
                              -1, -1,  4, -1, -1,
                              -1, -1,  4, -1, -1,
                              -1, -1,  4, -1,
                              -1,  4, -1,
                              -1, -1,  4, -1,
                              -1, -1,  4, -1,
                              -1, -1,  4};

  test_mdf_reference_problem<maxEntrPerRow,Scalar,LocalOrdinal,GlobalOrdinal,Node>(
    success,
    out,
    locNumRow,
    crs_row_map,
    crs_col_ind,
    crs_values,
    known_Perm,
    known_sln);
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2MDF, Test2, Scalar, LocalOrdinal, GlobalOrdinal)
{
  // Finite difference Advection-difussion problem with uniform spacing, D = 0.1, and v = {1,0}.
  //     0 = D lapl(c) - v div(c)
  // Matrix generated in matlab

  const LocalOrdinal locNumRow = 9;
  const LocalOrdinal maxEntrPerRow = 5;

  const LocalOrdinal crs_row_map[locNumRow + 1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45};
  const LocalOrdinal crs_col_ind[] = {
    0, 1, 2, 3, 6, 
    0, 1, 2, 4, 7, 
    0, 1, 2, 5, 8, 
    0, 3, 4, 5, 6, 
    1, 3, 4, 5, 7, 
    2, 3, 4, 5, 8, 
    0, 3, 6, 7, 8, 
    1, 4, 6, 7, 8, 
    2, 5, 6, 7, 8
  };
  const Scalar crs_values[] = {
    -0.6, -0.1,  0.9, -0.1, -0.1, 
     0.9, -0.6, -0.1, -0.1, -0.1, 
    -0.1,  0.9, -0.6, -0.1, -0.1, 
    -0.1, -0.6, -0.1,  0.9, -0.1, 
    -0.1,  0.9, -0.6, -0.1, -0.1, 
    -0.1, -0.1,  0.9, -0.6, -0.1, 
    -0.1, -0.1, -0.6, -0.1,  0.9, 
    -0.1, -0.1,  0.9, -0.6, -0.1, 
    -0.1, -0.1, -0.1,  0.9, -0.6
  };
  const LocalOrdinal known_Perm[] = {0, 3, 6, 1, 4, 5, 7, 8, 2};
  const Scalar known_sln[] = {5.16247188607349, 5.52439770923273, 6.238246086077, 4.89405015653464, 5.66305162810038, 5.1671248742544, 4.75093559248463, 5.88204162861661, 4.46367801745558};

  test_mdf_reference_problem<maxEntrPerRow,Scalar,LocalOrdinal,GlobalOrdinal,Node>(
    success,
    out,
    locNumRow,
    crs_row_map,
    crs_col_ind,
    crs_values,
    known_Perm,
    known_sln);
}


#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2MDF, Test0, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2MDF, Test1, Scalar, LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2MDF, Test2, Scalar, LocalOrdinal,GlobalOrdinal)
  
#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

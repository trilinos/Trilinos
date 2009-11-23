// ***********************************************************************
// 
//      Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
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


/*! \file Tifpack_UnitTestTemplate.cpp

\brief Tifpack Unit testing template.

This file demonstrates how you create a unit test for template code.

*/


#include <Teuchos_ConfigDefs.hpp>
#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tifpack_Version.hpp>
#include <iostream>

#include <Tifpack_UnitTestHelpers.hpp>
#include <Tifpack_IlukGraph.hpp>

namespace {
using Tpetra::global_size_t;
typedef tif_utest::Node Node;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(TifpackIlukGraph, IlukGraphTest0, LocalOrdinal, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Tifpack::Version();
  out << "Tifpack::Version(): " << version << std::endl;

  global_size_t num_rows_per_proc = 5;

  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > crsgraph = tif_utest::create_test_graph<LocalOrdinal,GlobalOrdinal,Node>(num_rows_per_proc);

  TEUCHOS_TEST_EQUALITY( crsgraph->getMap()->getNodeNumElements(), num_rows_per_proc, out, success)

  LocalOrdinal overlap_levels = 2;

  LocalOrdinal fill_levels = 0;

  Tifpack::IlukGraph<LocalOrdinal,GlobalOrdinal> iluk0_graph(crsgraph, fill_levels, overlap_levels);

  iluk0_graph.ConstructFilledGraph();

  int num_procs = crsgraph->getMap()->getComm()->getSize();
  size_t num_global_elements = num_rows_per_proc * num_procs;

  size_t num_global_rows = iluk0_graph.getL_Graph().getGlobalNumRows();

  TEUCHOS_TEST_EQUALITY( num_global_rows, num_global_elements, out, success)

  //The number of nonzeros in an ILU(0) graph should be the same as the
  //number of nonzeros in the input graph:

  size_t nnz0 = iluk0_graph.getL_Graph().getGlobalNumEntries() +
                iluk0_graph.getU_Graph().getGlobalNumEntries() +
                iluk0_graph.getNumGlobalDiagonals();

  size_t nnz_input = crsgraph->getGlobalNumEntries();

  TEUCHOS_TEST_EQUALITY( nnz0, nnz_input, out, success)

  fill_levels = 2;

  Tifpack::IlukGraph<int,int> iluk2_graph(crsgraph, fill_levels, overlap_levels);

  iluk2_graph.ConstructFilledGraph();

  TEUCHOS_TEST_EQUALITY( num_global_rows, num_global_elements, out, success)

  //The number of nonzeros in an ILU(2) graph should be greater than the
  //number of nonzeros in the ILU(0) graph:

  size_t nnz2 = iluk2_graph.getL_Graph().getGlobalNumEntries() +
                iluk2_graph.getU_Graph().getGlobalNumEntries() +
                iluk2_graph.getNumGlobalDiagonals();

  bool nnz2_greater_than_nnz0 = nnz2 > nnz0;
  TEUCHOS_TEST_EQUALITY( nnz2_greater_than_nnz0, true, out, success)
}

#define UNIT_TEST_GROUP_ORDINAL(LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( TifpackIlukGraph, IlukGraphTest0, LocalOrdinal,GlobalOrdinal)

UNIT_TEST_GROUP_ORDINAL(int, int)

}//namespace <anonymous>



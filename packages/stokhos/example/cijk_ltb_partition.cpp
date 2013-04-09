// $Id$
// $Source$
// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Stokhos.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example generates partitions the Cijk tensor for a lexicographic tree basis.\n");
    int d = 3;
    CLP.setOption("dimension", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    double drop = 1.0e-12;
    CLP.setOption("drop", &drop, "Drop tolerance");
    bool symmetric = true;
    CLP.setOption("symmetric", "asymmetric", &symmetric, "Use basis polynomials with symmetric PDF");
    int a_size = 100;
    CLP.setOption("a_size", &a_size, "Size of a (matrix) partition");
    int x_size = 100;
    CLP.setOption("x_size", &x_size, "Size of x (input vector) partition");
    int y_size = 100;
    CLP.setOption("y_size", &y_size, "Size of y (output vector) partition");

    // Parse arguments
    CLP.parse( argc, argv );

    // Basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    const double alpha = 1.0;
    const double beta = symmetric ? 1.0 : 2.0 ;
    for (int i=0; i<d; i++) {
      bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<int,double>(
                                p, alpha, beta, true));
    }
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<int,double,less_type> basis_type;
    RCP<const basis_type> basis = Teuchos::rcp(new basis_type(bases, drop));

    // Build LTB Cijk
    typedef Stokhos::LTBSparse3Tensor<int,double> Cijk_LTB_type;
    typedef Cijk_LTB_type::CijkNode node_type;
    Teuchos::RCP<Cijk_LTB_type> Cijk =
      computeTripleProductTensorLTB(*basis, symmetric);

    int sz = basis->size();
    std::cout << "basis size = " << sz
              << " num nonzero Cijk entries = " << Cijk->num_entries()
              << std::endl;

    // Setup partitions
    if (a_size > sz) a_size = sz;
    if (x_size > sz) x_size = sz;
    if (y_size > sz) y_size = sz;

    Teuchos::Array< Teuchos::RCP<const node_type> > node_stack;
    Teuchos::Array< int > index_stack;
    node_stack.push_back(Cijk->getHeadNode());
    index_stack.push_back(0);
    Teuchos::RCP<const node_type> node;
    int child_index;
    Teuchos::Array< Teuchos::RCP<const node_type> > partition_stack;
    while (node_stack.size() > 0) {
      node = node_stack.back();
      child_index = index_stack.back();

      // Leaf -- If we got here, just push this node into the partitions
      if (node->is_leaf) {
        partition_stack.push_back(node);
        node_stack.pop_back();
        index_stack.pop_back();
      }

      // Once sizes are small enough, push node onto partition stack
      else if (node->i_size <= y_size &&
               node->j_size <= a_size &&
               node->k_size <= x_size) {
        partition_stack.push_back(node);
        node_stack.pop_back();
        index_stack.pop_back();
      }

      // More children to process -- process them first
      else if (child_index < node->children.size()) {
        ++index_stack.back();
        node = node->children[child_index];
        node_stack.push_back(node);
        index_stack.push_back(0);
      }

      // No more children
      else {
        node_stack.pop_back();
        index_stack.pop_back();
      }

    }

    std::cout << "num partitions = " << partition_stack.size() << std::endl;
    for (int part=0; part<partition_stack.size(); ++part) {
      node = partition_stack[part];
      std::cout << "partition " << part << ":" << std::endl
                << "\ti-range: [" << node->i_begin << ","
                << node->i_begin+node->i_size << ")" << std::endl
                << "\tj-range: [" << node->j_begin << ","
                << node->j_begin+node->j_size << ")" << std::endl
                << "\tk-range: [" << node->k_begin << ","
                << node->k_begin+node->k_size << ")" << std::endl
                << "\tnum_non_zeros = " << node->total_num_entries << std::endl;
    }

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}

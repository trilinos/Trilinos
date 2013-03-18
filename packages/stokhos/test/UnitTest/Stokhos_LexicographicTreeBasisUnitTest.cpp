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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"

namespace TotalOrderBasisUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol, sparse_tol;
    OrdinalType p,d;

    UnitTestSetup() {
      rtol = 1e-12;
      atol = 1e-12;
      sparse_tol = 1e-12;
      d = 3;
      p = 5;
    }

  };

  typedef int ordinal_type;
  typedef double value_type;
  UnitTestSetup<ordinal_type,value_type> setup;

  template <typename ordinal_type>
  bool test_lexicographic_tree_coeffs(
    const Teuchos::Array<ordinal_type>& basis_orders,
    const ordinal_type p,
    const bool isotropic,
    Teuchos::FancyOStream& out) {
    bool success = true;

    // Build total order basis of dimension d
    ordinal_type d = basis_orders.size();
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(d);
    for (ordinal_type i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(basis_orders[i], true));
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,less_type> basis_type;
    Teuchos::RCP<const basis_type> basis = Teuchos::rcp(new basis_type(bases));

    // Build tree basis
    typedef Stokhos::LexicographicTreeBasisNode<ordinal_type> node_type;
    Teuchos::RCP<node_type> head =
      Stokhos::build_lexicographic_basis_tree(basis_orders(), p);

    // Check head block size is equal to the whole basis size
    TEUCHOS_TEST_EQUALITY(head->block_size, basis->size(), out, success);

    // Check tree is consistent
    Teuchos::Array< Teuchos::RCP<node_type> > node_stack;
    Teuchos::Array< ordinal_type > index_stack;
    node_stack.push_back(head);
    index_stack.push_back(0);
    Teuchos::RCP<node_type> node;
    ordinal_type child_index;
    ordinal_type block_size_expected;
    while (node_stack.size() > 0) {
      node = node_stack.back();
      child_index = index_stack.back();

      // Check block size is the sum of each child's block size
      // or the number of terms for leaves
      if (node->children.size() > 0) {
        block_size_expected = 0;
        for (ordinal_type i=0; i<node->children.size(); ++i)
          block_size_expected += node->children[i]->block_size;
      }
      else {
        block_size_expected = node->terms.size();
      }
      TEUCHOS_TEST_EQUALITY(node->block_size, block_size_expected,
                            out, success);

      // Check block size based on formula (only for isotropic)
      if (isotropic) {
        ordinal_type sum_prev = 0;
        Stokhos::MultiIndex<ordinal_type> term_prefix = node->terms[0];
        for (ordinal_type i=0; i<term_prefix.dimension()-1; ++i)
          sum_prev += term_prefix[i];
        ordinal_type d_prev = term_prefix.dimension()-1;
        ordinal_type my_p = std::min(p-sum_prev,basis_orders[d_prev]);
        ordinal_type d_left = d - d_prev;
        ordinal_type block_size_expected2 =
          Stokhos::n_choose_k(my_p+d_left,d_left);
        TEUCHOS_TEST_EQUALITY(node->block_size, block_size_expected2,
                              out, success);
      }

      if (child_index < node->terms.size() && node->children.size() > 0)
        out << node->terms[child_index] << " : block_size = "
            << node->children[child_index]->block_size << std::endl;

      // Leaf -- check global indices
      if (node->children.size() == 0) {
        TEUCHOS_TEST_EQUALITY(node->block_size, node->terms.size(),
                              out, success);
        for (ordinal_type i=0; i<node->terms.size(); ++i) {
          Stokhos::MultiIndex<ordinal_type> term = node->terms[i];
          ordinal_type index = node->index_begin + i;
          ordinal_type index_expected = basis->index(term);

          out << term << " : index = " << index
              << " index expected = " << index_expected << std::endl;
          TEUCHOS_TEST_EQUALITY(index, index_expected, out, success);
        }
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
    return success;
  }

  template <typename ordinal_type, typename value_type>
  bool test_lexicographic_tree_sparse_3_tensor(
    const Teuchos::Array<ordinal_type>& basis_orders,
    const ordinal_type p,
    const bool symmetric,
    const value_type sparse_tol,
    const value_type atol,
    const value_type rtol,
    Teuchos::FancyOStream& out) {

    bool success = true;
    ordinal_type d = basis_orders.size();

    // Build tensor product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(d);
    if (symmetric)
      for (ordinal_type i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(basis_orders[i], true));
    else
      for (ordinal_type i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<ordinal_type,value_type>(basis_orders[i], 1.0, 2.0, true));
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,less_type> basis_type;
    Teuchos::RCP<const basis_type> basis =
      Teuchos::rcp(new basis_type(bases, sparse_tol));

    // Build "standard" Cijk
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk =
      basis->computeTripleProductTensor();

    // Build LTB Cijk
    typedef Stokhos::LTBSparse3Tensor<ordinal_type, value_type> Cijk_LTB_type;
    Teuchos::RCP<Cijk_LTB_type> Cijk_LTB =
      computeTripleProductTensorLTB(*basis, symmetric);

    // Check sizes
    TEUCHOS_TEST_EQUALITY(Cijk->num_entries(), Cijk_LTB->num_entries(),
                          out, success);

    // Check non-zero structure
    typedef typename Cijk_LTB_type::CijkNode node_type;
    typedef Stokhos::ProductBasisUtils::Cijk_1D_Iterator<ordinal_type> Cijk_Iterator;

    Teuchos::Array< Teuchos::RCP<const node_type> > node_stack;
    Teuchos::Array< ordinal_type > index_stack;
    node_stack.push_back(Cijk_LTB->getHeadNode());
    index_stack.push_back(0);
    Teuchos::RCP<const node_type> node;
    ordinal_type child_index;
    while (node_stack.size() > 0) {
      node = node_stack.back();
      child_index = index_stack.back();

      // Leaf -- check Cijk indices and values
      if (node->is_leaf) {
        TEUCHOS_TEST_EQUALITY(node->my_num_entries, node->values.size(),
                              out, success);
        Cijk_Iterator cijk_iterator(node->p_i, node->p_j, node->p_k, symmetric);
        ordinal_type idx = 0;
        bool more = true;
        while (more) {
          value_type cijk = node->values[idx];
          ordinal_type I = node->i_begin + cijk_iterator.i;
          ordinal_type J = node->j_begin + cijk_iterator.j;
          ordinal_type K = node->k_begin + cijk_iterator.k;
          value_type cijk2 = Cijk->getValue(I,J,K);

          value_type tol = atol + std::abs(cijk2)*rtol;
          value_type err = std::abs(cijk-cijk2);
          bool s = err < tol;
          if (!s) {
            out << std::endl
                << "Check: rel_err( C(" << I << "," << J << "," << K << ") )"
                << " = " << "rel_err( " << cijk << ", " << cijk2 << " ) = "
                << err << " <= " << tol << " : ";
            if (s) out << "Passed.";
            else
              out << "Failed!";
            out << std::endl;
          }
          success = success && s;
          more = cijk_iterator.increment();
          ++idx;
        }
        TEUCHOS_TEST_EQUALITY(node->my_num_entries, idx, out, success);
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

    return success;
  }

  template <typename ordinal_type, typename value_type>
  bool test_lexicographic_tree_sparse_3_tensor_block(
    const Teuchos::Array<ordinal_type>& basis_orders,
    const ordinal_type p,
    const bool symmetric,
    const value_type sparse_tol,
    const value_type atol,
    const value_type rtol,
    Teuchos::FancyOStream& out) {

    bool success = true;
    ordinal_type d = basis_orders.size();

    // Build tensor product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(d);
    if (symmetric)
      for (ordinal_type i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(basis_orders[i], true));
    else
      for (ordinal_type i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<ordinal_type,value_type>(basis_orders[i], 1.0, 2.0, true));
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,less_type> basis_type;
    Teuchos::RCP<const basis_type> basis =
      Teuchos::rcp(new basis_type(bases, sparse_tol));

    // Build "standard" Cijk
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk =
      basis->computeTripleProductTensor();

    // Build LTB Cijk
    typedef Stokhos::LTBSparse3Tensor<ordinal_type, value_type> Cijk_LTB_type;
    Teuchos::RCP<Cijk_LTB_type> Cijk_LTB =
      computeTripleProductTensorLTBBlockLeaf(*basis, symmetric);

    // Check non-zero structure
    typedef typename Cijk_LTB_type::CijkNode node_type;
    typedef Stokhos::ProductBasisUtils::Cijk_1D_Iterator<ordinal_type> Cijk_Iterator;

    Teuchos::Array< Teuchos::RCP<const node_type> > node_stack;
    Teuchos::Array< ordinal_type > index_stack;
    node_stack.push_back(Cijk_LTB->getHeadNode());
    index_stack.push_back(0);
    Teuchos::RCP<const node_type> node;
    ordinal_type child_index;
    while (node_stack.size() > 0) {
      node = node_stack.back();
      child_index = index_stack.back();

      // Leaf -- check Cijk indices and values
      if (node->is_leaf) {
        TEUCHOS_TEST_EQUALITY(node->my_num_entries, node->values.size(),
                              out, success);
        ordinal_type idx = 0;
        for (ordinal_type i=0; i<node->i_size; ++i) {
          for (ordinal_type j=0; j<node->j_size; ++j) {
            const ordinal_type k0 = symmetric ? (i+j)%2 : 0;
            const ordinal_type kinc = symmetric ? 2 : 1;
            for (ordinal_type k=k0; k<node->k_size; k+=kinc) {
              value_type cijk = node->values[idx++];
              ordinal_type I = node->i_begin + i;
              ordinal_type J = node->j_begin + j;
              ordinal_type K = node->k_begin + k;
              value_type cijk2 = Cijk->getValue(I,J,K);

              value_type tol = atol + std::abs(cijk2)*rtol;
              value_type err = std::abs(cijk-cijk2);
              bool s = err < tol;
              if (!s) {
                out << std::endl
                    << "Check: rel_err( C(" << I << "," << J << "," 
                    << K << ") )"
                    << " = " << "rel_err( " << cijk << ", " << cijk2 << " ) = "
                    << err << " <= " << tol << " : ";
                if (s) out << "Passed.";
                else
                  out << "Failed!";
                out << std::endl;
              }
              success = success && s;
            }
          }
        }
        TEUCHOS_TEST_EQUALITY(node->my_num_entries, idx, out, success);
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

    return success;
  }

  template <typename ordinal_type, typename value_type>
  bool test_lexicographic_tree_sparse_3_tensor_multiply(
    const Teuchos::Array<ordinal_type>& basis_orders,
    const ordinal_type p,
    const bool symmetric,
    const value_type sparse_tol,
    const value_type atol,
    const value_type rtol,
    Teuchos::FancyOStream& out) {

    bool success = true;
    ordinal_type d = basis_orders.size();

    // Build tensor product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(d);
    if (symmetric)
      for (ordinal_type i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(basis_orders[i], true));
    else
      for (ordinal_type i=0; i<d; i++)
        bases[i] = Teuchos::rcp(new Stokhos::JacobiBasis<ordinal_type,value_type>(basis_orders[i], 1.0, 2.0, true));
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,less_type> basis_type;
    Teuchos::RCP<const basis_type> basis =
      Teuchos::rcp(new basis_type(bases, sparse_tol));

    // Build "standard" Cijk
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk =
      basis->computeTripleProductTensor();

    // Build "standard" algebraic expansion
    Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type, value_type> expn(
      basis, Cijk);

    // Build quadrature expansion for sin/cos
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type,value_type> > quad =
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<ordinal_type,value_type>(basis));
    Stokhos::QuadOrthogPolyExpansion<ordinal_type, value_type> quad_expn(
      basis, Cijk, quad);

    // Build flat LTB 3 tensor
    Teuchos::RCP< Stokhos::FlatLTBSparse3Tensor<ordinal_type,value_type> > 
      flat_Cijk =
      Stokhos::computeFlatTripleProductTensorLTB(*basis, symmetric);

    // Expansions
    Stokhos::OrthogPolyApprox<ordinal_type,value_type>
      x(basis), a(basis), b(basis), c1(basis), c2(basis);

    // Initialize a and b to reasonable values
    x.term(0, 0) = 1.0;
    for (ordinal_type i=0; i<d; i++)
      x.term(i, 1) = 0.1;
    quad_expn.sin(a,x);
    quad_expn.cos(b,x);

    // Do multiplications
    expn.times(c1,a,b);
    Stokhos::flatLTB3TensorMultiply<10>(c2, a, b, *flat_Cijk);

    // Test c1 == c2
    success = Stokhos::comparePCEs(c1, "c1", c2, "c2", rtol, atol, out);

    return success;
  }

  TEUCHOS_UNIT_TEST( LexicographicTreeCoefficients, Isotropic ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d, setup.p);
    success = test_lexicographic_tree_coeffs(basis_orders, setup.p, true, out);
  }

  TEUCHOS_UNIT_TEST( LexicographicTreeCoefficients, Anisotropic ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d);
    for (ordinal_type i=0; i<setup.d; ++i)
      basis_orders[i] = i+1;
    success = test_lexicographic_tree_coeffs(basis_orders, setup.d, false, out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3Tensor, Isotropic_Symmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d, setup.p);
    success = test_lexicographic_tree_sparse_3_tensor(
      basis_orders, setup.p, true, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3Tensor, Anisotropic_Symmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d);
    for (ordinal_type i=0; i<setup.d; ++i)
      basis_orders[i] = i+1;
    success = test_lexicographic_tree_sparse_3_tensor(
      basis_orders, setup.p, true, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3Tensor, Isotropic_Asymmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d, setup.p);
    success = test_lexicographic_tree_sparse_3_tensor(
      basis_orders, setup.p, false, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3Tensor, Anisotropic_Asymmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d);
    for (ordinal_type i=0; i<setup.d; ++i)
      basis_orders[i] = i+1;
    success = test_lexicographic_tree_sparse_3_tensor(
      basis_orders, setup.p, false, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3TensorBlock, Isotropic_Symmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d, setup.p);
    success = test_lexicographic_tree_sparse_3_tensor_block(
      basis_orders, setup.p, true, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3TensorBlockMultiply, Isotropic_Symmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d, setup.p);
    success = test_lexicographic_tree_sparse_3_tensor_multiply(
      basis_orders, setup.p, true, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3TensorBlockMultiply, Isotropic_Asymmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d, setup.p);
    success = test_lexicographic_tree_sparse_3_tensor_multiply(
      basis_orders, setup.p, false, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3TensorBlockMultiply, Anisotropic_Symmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d);
    for (ordinal_type i=0; i<setup.d; ++i)
      basis_orders[i] = i+1;
    success = test_lexicographic_tree_sparse_3_tensor_multiply(
      basis_orders, setup.p, true, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }

  TEUCHOS_UNIT_TEST( LTBSparse3TensorBlockMultiply, Anisotropic_Asymmetric ) {
    Teuchos::Array<ordinal_type> basis_orders(setup.d);
    for (ordinal_type i=0; i<setup.d; ++i)
      basis_orders[i] = i+1;
    success = test_lexicographic_tree_sparse_3_tensor_multiply(
      basis_orders, setup.p, false, setup.sparse_tol, setup.atol, setup.rtol,
      out);
  }
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return res;
}
